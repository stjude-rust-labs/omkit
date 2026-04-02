//! Maps VCF record coordinates from a source genome to a target genome using a
//! chain file.
//!
//! The [`CoordinateMapper`] takes each record's contig and position, constructs
//! a base interval spanning the `REF` allele, converts it to interbase
//! coordinates, and queries the [`chainfile::liftover::Machine`] for a mapping
//! in the target assembly.
//!
//! When exactly one full-span mapping is found, the record's contig and
//! position are updated in place and the `negative_strand` flag is set on the
//! record's state. When the `REF` span is only partially covered by chain
//! alignment blocks (e.g., the variant straddles a gap within a chain or
//! extends beyond chain coverage), each covered chunk is extracted as its own
//! record with the corresponding substring of the original `REF`. A single
//! valid chunk produces an [`Outcome::Accepted`] with a truncated `REF`, while
//! two or more valid chunks produce an [`Outcome::Split`]. When multiple
//! full-span mappings are found (i.e., ambiguous mapping from overlapping
//! chains), the behavior is controlled by [`AmbiguousMapping`]—either reject
//! the record or emit one copy per candidate.
//!
//! # Rejections
//!
//! A record is rejected in the following situations:
//!
//! * **Missing position** ([`RejectionReason::MissingPosition`]). The record
//!   has no `POS` field (e.g., telomeric breakends).
//! * **Missing or empty `REF`**
//!   ([`RejectionReason::MissingRefAllele`]). The `REF` allele is `.` or
//!   empty. This transformer requires a `REF` allele to construct the
//!   liftover interval.
//! * **No target** ([`RejectionReason::NoTarget`]). No chain interval covers
//!   the record's position at all (e.g., the contig does not appear in the
//!   chain file, or the position falls in an alignment gap with no coverage on
//!   either side).
//! * **Ambiguous mapping** ([`RejectionReason::AmbiguousMapping`]). The
//!   liftover produces multiple full-span results from overlapping chains and
//!   [`AmbiguousMapping::Reject`] is configured.
//! * **Straddle with no valid chunks**
//!   ([`RejectionReason::IndelStraddlesMultipleIntervals`]). The `REF` span is
//!   only partially covered by alignment blocks and every resulting chunk
//!   extends beyond the `REF` string.

use std::sync::Arc;
use std::sync::atomic::AtomicUsize;
use std::sync::atomic::Ordering;

use noodles::core::Position;
use omics::coordinate::Contig;
use omics::coordinate::Strand;
use omics::coordinate::base::Coordinate as BaseCoordinate;
use omics::coordinate::interval::base::Interval as BaseInterval;
use omics::coordinate::interval::interbase::Interval;
use omics::coordinate::position::base::Position as BasePosition;
use tracing::debug;
use tracing::trace;

use crate::MISSING;
use crate::Outcome;
use crate::Result;
use crate::Transformer;
use crate::record::VcfRecord;
use crate::rejection::RejectionReason;

/// Controls behavior when a variant maps to multiple locations.
#[derive(Clone, Debug, Default)]
pub enum AmbiguousMapping {
    /// Reject variants that map to multiple locations.
    #[default]
    Reject,

    /// Emit all candidate mappings as separate records.
    EmitAll,
}

/// Running statistics for [`CoordinateMapper`].
#[derive(Debug, Default)]
struct Stats {
    /// Number of records mapped to exactly one location.
    accepted: AtomicUsize,

    /// Number of records rejected because no chain interval covers the
    /// position.
    rejected_no_target: AtomicUsize,

    /// Number of records rejected because the `REF` span straddles
    /// chain alignment block boundaries.
    rejected_straddle: AtomicUsize,

    /// Number of records rejected due to ambiguous mapping.
    rejected_ambiguous: AtomicUsize,

    /// Number of records split into multiple mappings.
    split: AtomicUsize,
}

/// Maps VCF record coordinates using a chain file.
pub struct CoordinateMapper {
    /// The liftover machine.
    machine: Arc<chainfile::liftover::Machine>,

    /// How to handle ambiguous mappings.
    ambiguous: AmbiguousMapping,

    /// Running statistics.
    stats: Stats,
}

impl CoordinateMapper {
    /// Creates a new coordinate mapper.
    pub fn new(
        machine: impl Into<Arc<chainfile::liftover::Machine>>,
        ambiguous: AmbiguousMapping,
    ) -> Self {
        Self {
            machine: machine.into(),
            ambiguous,
            stats: Stats::default(),
        }
    }
}

/// Applies a single liftover result to a record, updating its contig and
/// position.
fn apply_mapping(record: &mut VcfRecord, query: &Interval) {
    let base_interval = query.clone().into_equivalent_base();
    let start = base_interval.start();
    let end = base_interval.end();

    let new_contig = String::from(start.contig().as_str());
    let is_negative = start.strand() == Strand::Negative;

    // NOTE: on the negative strand, `start()` returns the highest
    // coordinate. VCF `POS` is always the left-most (lowest) coordinate,
    // so we take the minimum of the two endpoints.
    let start_pos = start.position().get();
    let end_pos = end.position().get();
    let new_pos = start_pos.min(end_pos) as usize;

    *record.current_mut().reference_sequence_name_mut() = new_contig;
    *record.current_mut().variant_start_mut() = Position::new(new_pos);

    record.state_mut().negative_strand = is_negative;
}

impl Transformer for CoordinateMapper {
    fn transform(&self, record: VcfRecord) -> Result<Outcome> {
        let contig = String::from(record.current().reference_sequence_name());
        let position = match record.current().variant_start() {
            Some(pos) => pos.get(),
            None => {
                trace!(contig = %contig, "record has no position");
                return Ok(Outcome::Rejected(record, RejectionReason::MissingPosition));
            }
        };

        trace!(contig = %contig, position = position, "coordinate mapping record");

        let ref_bases = record.current().reference_bases();
        let ref_len = ref_bases.len();

        // NOTE: this transformer rejects variants without a `REF` allele
        // because one is required to construct the liftover interval.
        if ref_bases == MISSING || ref_bases.is_empty() {
            trace!(contig = %contig, position = position, "REF is missing");
            return Ok(Outcome::Rejected(record, RejectionReason::MissingRefAllele));
        }

        // SAFETY: `vcf_pos` came from `noodles::core::Position`, which wraps
        // `NonZero<usize>`, so it is always >= 1.
        let start_pos = BasePosition::try_new(position as u64).unwrap();

        // SAFETY: the early return above guarantees `ref_len` >= 1, so
        // `vcf_pos + ref_len - 1` >= `vcf_pos` >= 1.
        let end_pos = BasePosition::try_new((position + ref_len - 1) as u64).unwrap();

        // NOTE: VCF coordinates are defined relative to the reference genome,
        // which corresponds to the positive strand in the chain file's
        // coordinate system.
        let start =
            BaseCoordinate::new(Contig::new_unchecked(&contig), Strand::Positive, start_pos);
        let end = BaseCoordinate::new(Contig::new_unchecked(&contig), Strand::Positive, end_pos);

        // SAFETY: `start` and `end` are constructed with the same contig
        // and strand, so `try_new()` will always succeed.
        let interval = BaseInterval::try_new(start, end)
            .unwrap()
            .into_equivalent_interbase();

        let results = match self.machine.liftover(interval) {
            Some(results) if !results.is_empty() => results,
            _ => {
                trace!(contig = %contig, position = position, "no target found");
                self.stats
                    .rejected_no_target
                    .fetch_add(1, Ordering::Relaxed);
                return Ok(Outcome::Rejected(record, RejectionReason::NoTarget));
            }
        };

        // Results are grouped by chain. Multiple results means the query
        // interval is covered by multiple chains (ambiguous mapping).
        if results.len() > 1 {
            match self.ambiguous {
                AmbiguousMapping::Reject => {
                    trace!(
                        contig = %contig,
                        position = position,
                        candidates = results.len(),
                        "rejected due to ambiguous mapping"
                    );
                    self.stats
                        .rejected_ambiguous
                        .fetch_add(1, Ordering::Relaxed);
                    return Ok(Outcome::Rejected(
                        record,
                        RejectionReason::AmbiguousMapping {
                            candidates: results.len(),
                        },
                    ));
                }
                AmbiguousMapping::EmitAll => {
                    trace!(
                        contig = %contig,
                        position = position,
                        candidates = results.len(),
                        "splitting into multiple mappings"
                    );
                    self.stats.split.fetch_add(1, Ordering::Relaxed);
                    let records = results
                        .iter()
                        .flat_map(|result| {
                            result.segments().iter().map(|seg| {
                                let mut r = record.clone();
                                apply_mapping(&mut r, seg.query());
                                r
                            })
                        })
                        .collect();
                    return Ok(Outcome::Split(records));
                }
            }
        }

        // Single chain — segments are non-overlapping blocks from that
        // chain. A single segment covering the full REF is the common
        // case; multiple segments mean the REF straddles a gap.
        let segments = results[0].segments();
        let ref_bases = String::from(record.current().reference_bases());
        let interbase_origin = start_pos.get() - 1;

        let mut accepted = Vec::new();
        let mut rejected = Vec::new();

        for seg in segments {
            let chunk_offset =
                (seg.reference().start().position().get() - interbase_origin) as usize;
            let chunk_len = seg.reference().count_entities() as usize;

            if chunk_offset + chunk_len > ref_bases.len() {
                trace!(
                    contig = %contig,
                    position = position,
                    chunk_offset,
                    chunk_len,
                    "chunk extends beyond REF"
                );
                rejected.push((
                    record.clone(),
                    RejectionReason::IndelStraddlesMultipleIntervals,
                ));
                continue;
            }

            let chunk_ref = &ref_bases[chunk_offset..][..chunk_len];
            let mut r = record.clone();
            apply_mapping(&mut r, seg.query());
            *r.current_mut().reference_bases_mut() = String::from(chunk_ref);
            accepted.push(r);
        }

        match accepted.len() {
            0 => {
                self.stats.rejected_straddle.fetch_add(1, Ordering::Relaxed);
                let (r, reason) = rejected.into_iter().next().unwrap();
                Ok(Outcome::Rejected(r, reason))
            }
            1 => {
                self.stats.accepted.fetch_add(1, Ordering::Relaxed);
                Ok(Outcome::Accepted(accepted.into_iter().next().unwrap()))
            }
            _ => {
                self.stats.split.fetch_add(1, Ordering::Relaxed);
                Ok(Outcome::Split(accepted))
            }
        }
    }

    fn finish(&self) {
        let accepted = self.stats.accepted.load(Ordering::Relaxed);
        let no_target = self.stats.rejected_no_target.load(Ordering::Relaxed);
        let straddle = self.stats.rejected_straddle.load(Ordering::Relaxed);
        let ambiguous = self.stats.rejected_ambiguous.load(Ordering::Relaxed);
        let split = self.stats.split.load(Ordering::Relaxed);
        let total = accepted + no_target + straddle + ambiguous + split;

        debug!(
            total,
            accepted,
            rejected_no_target = no_target,
            rejected_straddle = straddle,
            rejected_ambiguous = ambiguous,
            split,
            "CoordinateMapper finished"
        );
    }
}

#[cfg(test)]
mod tests {
    use std::sync::Arc;

    use chainfile::Reader;
    use chainfile::liftover::machine;
    use noodles::core::Position;

    use super::*;
    use crate::Outcome;
    use crate::Transformer;
    use crate::record::VcfRecord;
    use crate::rejection::RejectionReason;

    /// Builds a [`chainfile::liftover::Machine`] from raw chain file bytes.
    fn build_machine(data: &[u8]) -> Arc<chainfile::liftover::Machine> {
        let reader = Reader::new(data);
        Arc::new(machine::Builder.try_build_from(reader).unwrap())
    }

    /// Builds a [`VcfRecord`] with the given contig, position, and `REF`
    /// allele.
    fn test_record(contig: &str, pos: usize, ref_bases: &str) -> VcfRecord {
        let mut inner = noodles::vcf::variant::RecordBuf::default();
        *inner.reference_sequence_name_mut() = String::from(contig);
        *inner.variant_start_mut() = Some(Position::try_from(pos).unwrap());
        *inner.reference_bases_mut() = String::from(ref_bases);
        VcfRecord::new(inner)
    }

    /// A simple positive-to-positive chain: `seq0` (size 10) → `seq1`
    /// (size 10), covering the full range `[0, 10)`.
    const CHAIN_POS_TO_POS: &[u8] = b"chain 0 seq0 10 + 0 10 seq1 10 + 0 10 0\n10";

    /// A positive-to-negative chain: `seq0` (size 10) → `seq1` (size 10),
    /// covering the full range `[0, 10)`.
    const CHAIN_POS_TO_NEG: &[u8] = b"chain 0 seq0 10 + 0 10 seq1 10 - 0 10 0\n10";

    /// Two overlapping chains on `seq0` that both cover `[0, 10)`, mapping
    /// to `seq1` (score 100) and `seq2` (score 50) respectively.
    const CHAIN_AMBIGUOUS: &[u8] =
        b"chain 100 seq0 10 + 0 10 seq1 10 + 0 10 0\n10\n\nchain 50 seq0 10 + 0 10 seq2 10 + 0 10 1\n10";

    /// A chain with a gap: two aligned blocks of 4 bases each, separated
    /// by a gap of 2 in both reference and query. Covers `seq0` [0, 10) →
    /// `seq1` [0, 10) with a hole at interbase [4, 6).
    const CHAIN_WITH_GAP: &[u8] = b"chain 0 seq0 10 + 0 10 seq1 10 + 0 10 0\n4\t2\t2\n4";

    #[test]
    fn single_mapping_updates_contig_and_position() {
        let mapper =
            CoordinateMapper::new(build_machine(CHAIN_POS_TO_POS), AmbiguousMapping::Reject);
        // VCF position 2 (1-based) → interbase [1, 2) → maps to seq1 interbase
        // [1, 2) → base position 2.
        let record = test_record("seq0", 2, "A");
        let outcome = mapper.transform(record).unwrap();

        match outcome {
            Outcome::Accepted(r) => {
                assert_eq!(r.current().reference_sequence_name(), "seq1");
                assert_eq!(r.current().variant_start(), Position::new(2));
                assert!(!r.state().negative_strand);
            }
            other => panic!("expected `Accepted`, got {other:?}"),
        }
    }

    #[test]
    fn original_record_is_preserved() {
        let mapper =
            CoordinateMapper::new(build_machine(CHAIN_POS_TO_POS), AmbiguousMapping::Reject);
        let record = test_record("seq0", 3, "A");
        let outcome = mapper.transform(record).unwrap();

        match outcome {
            Outcome::Accepted(r) => {
                assert_eq!(r.original_contig(), "seq0");
                assert_eq!(r.original_position(), Position::new(3));
                assert_eq!(r.current().reference_sequence_name(), "seq1");
                assert_eq!(r.current().variant_start(), Position::new(3));
            }
            other => panic!("expected `Accepted`, got {other:?}"),
        }
    }

    #[test]
    fn negative_strand_mapping_sets_flag() {
        let mapper =
            CoordinateMapper::new(build_machine(CHAIN_POS_TO_NEG), AmbiguousMapping::Reject);
        // SNV at position 2 on the positive strand of a 10-base
        // chromosome. The chain maps to the negative strand, so the
        // left-most forward-strand position is 9.
        let record = test_record("seq0", 2, "A");
        let outcome = mapper.transform(record).unwrap();

        match outcome {
            Outcome::Accepted(r) => {
                assert!(r.state().negative_strand);
                assert_eq!(r.current().reference_sequence_name(), "seq1");
                assert_eq!(r.current().variant_start(), Position::new(9));
            }
            other => panic!("expected `Accepted`, got {other:?}"),
        }
    }

    #[test]
    fn missing_position_is_rejected() {
        let mapper =
            CoordinateMapper::new(build_machine(CHAIN_POS_TO_POS), AmbiguousMapping::Reject);
        let mut inner = noodles::vcf::variant::RecordBuf::default();
        *inner.reference_sequence_name_mut() = String::from("seq0");
        *inner.variant_start_mut() = None;
        *inner.reference_bases_mut() = String::from("A");
        let record = VcfRecord::new(inner);

        let outcome = mapper.transform(record).unwrap();
        assert!(matches!(
            outcome,
            Outcome::Rejected(_, RejectionReason::MissingPosition)
        ));
    }

    #[test]
    fn missing_ref_allele_is_rejected() {
        let mapper =
            CoordinateMapper::new(build_machine(CHAIN_POS_TO_POS), AmbiguousMapping::Reject);
        let record = test_record("seq0", 2, ".");
        let outcome = mapper.transform(record).unwrap();
        assert!(matches!(
            outcome,
            Outcome::Rejected(_, RejectionReason::MissingRefAllele)
        ));
    }

    #[test]
    fn empty_ref_allele_is_rejected() {
        let mapper =
            CoordinateMapper::new(build_machine(CHAIN_POS_TO_POS), AmbiguousMapping::Reject);
        let record = test_record("seq0", 2, "");
        let outcome = mapper.transform(record).unwrap();
        assert!(matches!(
            outcome,
            Outcome::Rejected(_, RejectionReason::MissingRefAllele)
        ));
    }

    #[test]
    fn no_target_is_rejected() {
        let mapper =
            CoordinateMapper::new(build_machine(CHAIN_POS_TO_POS), AmbiguousMapping::Reject);
        // `seq9` does not exist in the chain file.
        let record = test_record("seq9", 2, "A");
        let outcome = mapper.transform(record).unwrap();
        assert!(matches!(
            outcome,
            Outcome::Rejected(_, RejectionReason::NoTarget)
        ));
    }

    #[test]
    fn ambiguous_mapping_rejects_by_default() {
        let mapper =
            CoordinateMapper::new(build_machine(CHAIN_AMBIGUOUS), AmbiguousMapping::Reject);
        let record = test_record("seq0", 2, "A");
        let outcome = mapper.transform(record).unwrap();

        match outcome {
            Outcome::Rejected(_, RejectionReason::AmbiguousMapping { candidates }) => {
                assert_eq!(candidates, 2);
            }
            other => panic!("expected `Rejected(AmbiguousMapping)`, got {other:?}"),
        }
    }

    #[test]
    fn ambiguous_mapping_emits_all_highest_score_first() {
        let mapper =
            CoordinateMapper::new(build_machine(CHAIN_AMBIGUOUS), AmbiguousMapping::EmitAll);
        let record = test_record("seq0", 2, "A");
        let outcome = mapper.transform(record).unwrap();

        match outcome {
            Outcome::Split(records) => {
                assert_eq!(records.len(), 2);
                // The highest-scoring chain (score 100 → `seq1`) comes first.
                assert_eq!(records[0].current().reference_sequence_name(), "seq1");
                assert_eq!(records[1].current().reference_sequence_name(), "seq2");
            }
            other => panic!("expected `Split`, got {other:?}"),
        }
    }

    #[test]
    fn ambiguous_mapping_with_different_scores_still_rejected() {
        // Even when one chain scores much higher than the other, if both
        // cover the query interval the mapping is still ambiguous.
        let data =
            b"chain 1000 seq0 10 + 0 10 seq1 10 + 0 10 0\n10\n\nchain 1 seq0 10 + 0 10 seq2 10 + 0 10 1\n10";
        let mapper = CoordinateMapper::new(build_machine(data), AmbiguousMapping::Reject);
        let record = test_record("seq0", 5, "A");
        let outcome = mapper.transform(record).unwrap();

        match outcome {
            Outcome::Rejected(_, RejectionReason::AmbiguousMapping { candidates }) => {
                assert_eq!(candidates, 2);
            }
            other => panic!("expected `Rejected(AmbiguousMapping)`, got {other:?}"),
        }
    }

    #[test]
    fn snv_maps_single_base_interval() {
        let mapper =
            CoordinateMapper::new(build_machine(CHAIN_POS_TO_POS), AmbiguousMapping::Reject);
        // SNV: REF=`A` (1 base) at position 5 → interval spans [5, 5].
        let record = test_record("seq0", 5, "A");
        let outcome = mapper.transform(record).unwrap();

        match outcome {
            Outcome::Accepted(r) => {
                assert_eq!(r.current().reference_sequence_name(), "seq1");
                assert_eq!(r.current().variant_start(), Position::new(5));
            }
            other => panic!("expected `Accepted`, got {other:?}"),
        }
    }

    #[test]
    fn insertion_maps_padding_base_only() {
        let mapper =
            CoordinateMapper::new(build_machine(CHAIN_POS_TO_POS), AmbiguousMapping::Reject);
        // Insertion: REF=`A`, ALT=`ACGT`. The interval is determined solely
        // by REF length (1 base), so this behaves like a SNV for mapping.
        let record = test_record("seq0", 3, "A");
        let outcome = mapper.transform(record).unwrap();

        match outcome {
            Outcome::Accepted(r) => {
                assert_eq!(r.current().reference_sequence_name(), "seq1");
                assert_eq!(r.current().variant_start(), Position::new(3));
            }
            other => panic!("expected `Accepted`, got {other:?}"),
        }
    }

    #[test]
    fn deletion_maps_full_ref_span() {
        let mapper =
            CoordinateMapper::new(build_machine(CHAIN_POS_TO_POS), AmbiguousMapping::Reject);
        // Deletion: REF=`ACGT` (4 bases) at position 2 → interval spans
        // [2, 5]. The mapped position should still be the start.
        let record = test_record("seq0", 2, "ACGT");
        let outcome = mapper.transform(record).unwrap();

        match outcome {
            Outcome::Accepted(r) => {
                assert_eq!(r.current().reference_sequence_name(), "seq1");
                assert_eq!(r.current().variant_start(), Position::new(2));
            }
            other => panic!("expected `Accepted`, got {other:?}"),
        }
    }

    #[test]
    fn mnv_maps_full_ref_span() {
        let mapper =
            CoordinateMapper::new(build_machine(CHAIN_POS_TO_POS), AmbiguousMapping::Reject);
        // MNV: REF=`ACG`, ALT=`TGA` (3 bases each) at position 4 →
        // interval spans [4, 6].
        let record = test_record("seq0", 4, "ACG");
        let outcome = mapper.transform(record).unwrap();

        match outcome {
            Outcome::Accepted(r) => {
                assert_eq!(r.current().reference_sequence_name(), "seq1");
                assert_eq!(r.current().variant_start(), Position::new(4));
            }
            other => panic!("expected `Accepted`, got {other:?}"),
        }
    }

    #[test]
    fn negative_strand_deletion_maps_correctly() {
        let mapper =
            CoordinateMapper::new(build_machine(CHAIN_POS_TO_NEG), AmbiguousMapping::Reject);
        // Deletion on a chain that maps to the negative strand.
        // REF=`ACGT` (4 bases) at position 2 → base interval [2, 5],
        // interbase [1, 5). The chain maps to the negative strand of a
        // 10-base chromosome, so the left-most forward-strand position
        // is 6.
        let record = test_record("seq0", 2, "ACGT");
        let outcome = mapper.transform(record).unwrap();

        match outcome {
            Outcome::Accepted(r) => {
                assert!(r.state().negative_strand);
                assert_eq!(r.current().reference_sequence_name(), "seq1");
                assert_eq!(r.current().variant_start(), Position::new(6));
            }
            other => panic!("expected `Accepted`, got {other:?}"),
        }
    }

    #[test]
    fn negative_strand_snv_at_position_5() {
        let mapper =
            CoordinateMapper::new(build_machine(CHAIN_POS_TO_NEG), AmbiguousMapping::Reject);
        // SNV at position 5. Chain maps seq0:+:[0,10) → seq1:-:[0,10)
        // on a 10-base chromosome. Interbase [4, 5) maps to
        // forward-strand position 6.
        let record = test_record("seq0", 5, "A");
        let outcome = mapper.transform(record).unwrap();

        match outcome {
            Outcome::Accepted(r) => {
                assert!(r.state().negative_strand);
                assert_eq!(r.current().reference_sequence_name(), "seq1");
                assert_eq!(r.current().variant_start(), Position::new(6));
            }
            other => panic!("expected `Accepted`, got {other:?}"),
        }
    }

    #[test]
    fn negative_strand_insertion_maps_padding_base() {
        let mapper =
            CoordinateMapper::new(build_machine(CHAIN_POS_TO_NEG), AmbiguousMapping::Reject);
        // Insertion at position 3 (REF=`A`, 1 base). Interbase [2, 3)
        // maps to forward-strand position 8.
        let record = test_record("seq0", 3, "A");
        let outcome = mapper.transform(record).unwrap();

        match outcome {
            Outcome::Accepted(r) => {
                assert!(r.state().negative_strand);
                assert_eq!(r.current().reference_sequence_name(), "seq1");
                assert_eq!(r.current().variant_start(), Position::new(8));
            }
            other => panic!("expected `Accepted`, got {other:?}"),
        }
    }

    #[test]
    fn negative_strand_mnv_maps_full_ref_span() {
        let mapper =
            CoordinateMapper::new(build_machine(CHAIN_POS_TO_NEG), AmbiguousMapping::Reject);
        // MNV at position 4 with REF=`ACG` (3 bases). Base interval
        // [4, 6], interbase [3, 6). Maps to forward-strand position 5.
        let record = test_record("seq0", 4, "ACG");
        let outcome = mapper.transform(record).unwrap();

        match outcome {
            Outcome::Accepted(r) => {
                assert!(r.state().negative_strand);
                assert_eq!(r.current().reference_sequence_name(), "seq1");
                assert_eq!(r.current().variant_start(), Position::new(5));
            }
            other => panic!("expected `Accepted`, got {other:?}"),
        }
    }

    #[test]
    fn negative_strand_snv_at_first_position() {
        let mapper =
            CoordinateMapper::new(build_machine(CHAIN_POS_TO_NEG), AmbiguousMapping::Reject);
        // SNV at position 1 (left edge). Interbase [0, 1) maps to
        // forward-strand position 10.
        let record = test_record("seq0", 1, "A");
        let outcome = mapper.transform(record).unwrap();

        match outcome {
            Outcome::Accepted(r) => {
                assert!(r.state().negative_strand);
                assert_eq!(r.current().variant_start(), Position::new(10));
            }
            other => panic!("expected `Accepted`, got {other:?}"),
        }
    }

    #[test]
    fn negative_strand_snv_at_last_position() {
        let mapper =
            CoordinateMapper::new(build_machine(CHAIN_POS_TO_NEG), AmbiguousMapping::Reject);
        // SNV at position 10 (right edge). Interbase [9, 10) maps to
        // forward-strand position 1.
        let record = test_record("seq0", 10, "A");
        let outcome = mapper.transform(record).unwrap();

        match outcome {
            Outcome::Accepted(r) => {
                assert!(r.state().negative_strand);
                assert_eq!(r.current().variant_start(), Position::new(1));
            }
            other => panic!("expected `Accepted`, got {other:?}"),
        }
    }

    #[test]
    fn negative_strand_preserves_original() {
        let mapper =
            CoordinateMapper::new(build_machine(CHAIN_POS_TO_NEG), AmbiguousMapping::Reject);
        let record = test_record("seq0", 4, "ACGT");
        let outcome = mapper.transform(record).unwrap();

        match outcome {
            Outcome::Accepted(r) => {
                assert_eq!(r.original_contig(), "seq0");
                assert_eq!(r.original_position(), Position::new(4));
                assert_eq!(r.original_ref_bases(), "ACGT");
                assert_eq!(r.current().reference_sequence_name(), "seq1");
                assert_eq!(r.current().variant_start(), Position::new(4));
            }
            other => panic!("expected `Accepted`, got {other:?}"),
        }
    }

    #[test]
    fn ref_at_chain_boundary_maps_correctly() {
        let mapper =
            CoordinateMapper::new(build_machine(CHAIN_POS_TO_POS), AmbiguousMapping::Reject);
        // The chain covers base positions 1..=10. A SNV at position 10
        // sits right at the boundary.
        let record = test_record("seq0", 10, "A");
        let outcome = mapper.transform(record).unwrap();

        match outcome {
            Outcome::Accepted(r) => {
                assert_eq!(r.current().reference_sequence_name(), "seq1");
                assert_eq!(r.current().variant_start(), Position::new(10));
            }
            other => panic!("expected `Accepted`, got {other:?}"),
        }
    }

    #[test]
    fn ref_spanning_beyond_chain_is_truncated() {
        let mapper =
            CoordinateMapper::new(build_machine(CHAIN_POS_TO_POS), AmbiguousMapping::Reject);
        // The chain covers base positions 1..=10. A deletion whose REF
        // spans positions 8..=12 extends beyond the chain. Only the
        // portion within chain coverage (positions 8..=10) is kept.
        let record = test_record("seq0", 8, "ACGTA");
        let outcome = mapper.transform(record).unwrap();

        match outcome {
            Outcome::Accepted(r) => {
                assert_eq!(r.current().reference_sequence_name(), "seq1");
                assert_eq!(r.current().variant_start(), Position::new(8));
                assert_eq!(r.current().reference_bases(), "ACG");
            }
            other => panic!("expected `Accepted`, got {other:?}"),
        }
    }

    #[test]
    fn ref_with_iupac_ambiguity_codes() {
        let mapper =
            CoordinateMapper::new(build_machine(CHAIN_POS_TO_POS), AmbiguousMapping::Reject);
        // REF=`N` (IUPAC ambiguity code). The mapper only uses REF
        // length, not the actual bases, so this should map normally.
        let record = test_record("seq0", 3, "N");
        let outcome = mapper.transform(record).unwrap();

        match outcome {
            Outcome::Accepted(r) => {
                assert_eq!(r.current().reference_sequence_name(), "seq1");
                assert_eq!(r.current().variant_start(), Position::new(3));
            }
            other => panic!("expected `Accepted`, got {other:?}"),
        }
    }

    #[test]
    fn ref_and_alt_both_missing() {
        let mapper =
            CoordinateMapper::new(build_machine(CHAIN_POS_TO_POS), AmbiguousMapping::Reject);
        // Both REF and ALT are `.` (missing). The mapper rejects on
        // REF=`.` before it ever considers ALT.
        let record = test_record("seq0", 2, ".");
        let outcome = mapper.transform(record).unwrap();
        assert!(matches!(
            outcome,
            Outcome::Rejected(_, RejectionReason::MissingRefAllele)
        ));
    }

    #[test]
    fn symbolic_alt_allele_maps_normally() {
        let mapper =
            CoordinateMapper::new(build_machine(CHAIN_POS_TO_POS), AmbiguousMapping::Reject);
        // Structural variant with symbolic ALT (`<DEL>`). The mapper
        // ignores ALT entirely, so this maps based on REF alone.
        let record = test_record("seq0", 4, "ACGT");
        let outcome = mapper.transform(record).unwrap();

        match outcome {
            Outcome::Accepted(r) => {
                assert_eq!(r.current().reference_sequence_name(), "seq1");
                assert_eq!(r.current().variant_start(), Position::new(4));
            }
            other => panic!("expected `Accepted`, got {other:?}"),
        }
    }

    #[test]
    fn lowercase_ref_uses_correct_length() {
        let mapper =
            CoordinateMapper::new(build_machine(CHAIN_POS_TO_POS), AmbiguousMapping::Reject);
        // Lowercase REF=`acgt`. The mapper uses `len()` which works
        // regardless of case.
        let record = test_record("seq0", 2, "acgt");
        let outcome = mapper.transform(record).unwrap();

        match outcome {
            Outcome::Accepted(r) => {
                assert_eq!(r.current().reference_sequence_name(), "seq1");
                assert_eq!(r.current().variant_start(), Position::new(2));
            }
            other => panic!("expected `Accepted`, got {other:?}"),
        }
    }

    // ---- straddle-split tests (CHAIN_WITH_GAP) ----
    //
    // CHAIN_WITH_GAP has two 4-base aligned blocks separated by a 2-base
    // gap in both reference and query:
    //   Block 1: ref [0, 4) → query [0, 4)   (base positions 1..=4)
    //   Gap:     ref [4, 6), query [4, 6)     (base positions 5..=6)
    //   Block 2: ref [6, 10) → query [6, 10)  (base positions 7..=10)

    #[test]
    fn gap_chain_snv_within_first_block() {
        let mapper = CoordinateMapper::new(build_machine(CHAIN_WITH_GAP), AmbiguousMapping::Reject);
        // SNV at position 2 (within block 1) maps normally.
        let record = test_record("seq0", 2, "A");
        let outcome = mapper.transform(record).unwrap();

        match outcome {
            Outcome::Accepted(r) => {
                assert_eq!(r.current().reference_sequence_name(), "seq1");
                assert_eq!(r.current().variant_start(), Position::new(2));
            }
            other => panic!("expected `Accepted`, got {other:?}"),
        }
    }

    #[test]
    fn gap_chain_snv_within_second_block() {
        let mapper = CoordinateMapper::new(build_machine(CHAIN_WITH_GAP), AmbiguousMapping::Reject);
        // SNV at position 8 (within block 2) maps normally.
        let record = test_record("seq0", 8, "A");
        let outcome = mapper.transform(record).unwrap();

        match outcome {
            Outcome::Accepted(r) => {
                assert_eq!(r.current().reference_sequence_name(), "seq1");
                assert_eq!(r.current().variant_start(), Position::new(8));
            }
            other => panic!("expected `Accepted`, got {other:?}"),
        }
    }

    #[test]
    fn gap_chain_deletion_within_one_block() {
        let mapper = CoordinateMapper::new(build_machine(CHAIN_WITH_GAP), AmbiguousMapping::Reject);
        // Deletion REF=`ACGT` at position 1, spanning base positions
        // 1..=4 — fully within block 1.
        let record = test_record("seq0", 1, "ACGT");
        let outcome = mapper.transform(record).unwrap();

        match outcome {
            Outcome::Accepted(r) => {
                assert_eq!(r.current().reference_sequence_name(), "seq1");
                assert_eq!(r.current().variant_start(), Position::new(1));
                assert_eq!(r.current().reference_bases(), "ACGT");
            }
            other => panic!("expected `Accepted`, got {other:?}"),
        }
    }

    #[test]
    fn gap_chain_variant_in_gap_is_rejected() {
        let mapper = CoordinateMapper::new(build_machine(CHAIN_WITH_GAP), AmbiguousMapping::Reject);
        // SNV at position 5, which falls in the gap between the two
        // blocks. No chain interval covers this position.
        let record = test_record("seq0", 5, "A");
        let outcome = mapper.transform(record).unwrap();
        assert!(matches!(
            outcome,
            Outcome::Rejected(_, RejectionReason::NoTarget)
        ));
    }

    #[test]
    fn gap_chain_ref_extending_into_gap_is_truncated() {
        let mapper = CoordinateMapper::new(build_machine(CHAIN_WITH_GAP), AmbiguousMapping::Reject);
        // Deletion at position 3 with REF=`ACGT` (4 bases), spanning
        // base positions 3..=6. This overlaps block 1 at positions 3..=4
        // and extends into the gap at 5..=6. Only the covered portion
        // (positions 3..=4) is kept.
        let record = test_record("seq0", 3, "ACGT");
        let outcome = mapper.transform(record).unwrap();

        match outcome {
            Outcome::Accepted(r) => {
                assert_eq!(r.current().reference_sequence_name(), "seq1");
                assert_eq!(r.current().variant_start(), Position::new(3));
                assert_eq!(r.current().reference_bases(), "AC");
            }
            other => panic!("expected `Accepted`, got {other:?}"),
        }
    }

    #[test]
    fn gap_chain_ref_spanning_both_blocks_splits() {
        let mapper = CoordinateMapper::new(build_machine(CHAIN_WITH_GAP), AmbiguousMapping::Reject);
        // Deletion at position 1 with REF=`ACGTNNACGT` (10 bases),
        // spanning base positions 1..=10 — across both blocks and the
        // gap. The liftover returns two partial results (one per block),
        // and the variant is split into two records.
        let record = test_record("seq0", 1, "ACGTNNACGT");
        let outcome = mapper.transform(record).unwrap();

        match outcome {
            Outcome::Split(records) => {
                assert_eq!(records.len(), 2);

                // Chunk 1: first 4 bases of the REF, mapped to block 1.
                assert_eq!(records[0].current().reference_bases(), "ACGT");
                assert_eq!(records[0].current().variant_start(), Position::new(1));

                // Chunk 2: last 4 bases of the REF (positions 7..=10),
                // mapped to block 2.
                assert_eq!(records[1].current().reference_bases(), "ACGT");
                assert_eq!(records[1].current().variant_start(), Position::new(7));
            }
            other => panic!("expected `Split`, got {other:?}"),
        }
    }

    #[test]
    fn gap_chain_ref_spanning_both_blocks_preserves_original() {
        let mapper = CoordinateMapper::new(build_machine(CHAIN_WITH_GAP), AmbiguousMapping::Reject);
        let record = test_record("seq0", 1, "ACGTNNACGT");
        let outcome = mapper.transform(record).unwrap();

        match outcome {
            Outcome::Split(records) => {
                for r in &records {
                    assert_eq!(r.original_contig(), "seq0");
                    assert_eq!(r.original_position(), Position::new(1));
                    assert_eq!(r.original_ref_bases(), "ACGTNNACGT");
                }
            }
            other => panic!("expected `Split`, got {other:?}"),
        }
    }

    #[test]
    fn gap_chain_split_extracts_correct_ref_substrings() {
        let mapper = CoordinateMapper::new(build_machine(CHAIN_WITH_GAP), AmbiguousMapping::Reject);
        // REF=`1234XX6789` — distinct characters to verify correct
        // substring extraction. Positions 1..=4 map to block 1,
        // positions 7..=10 map to block 2.
        let record = test_record("seq0", 1, "1234XX6789");
        let outcome = mapper.transform(record).unwrap();

        match outcome {
            Outcome::Split(records) => {
                assert_eq!(records.len(), 2);
                assert_eq!(records[0].current().reference_bases(), "1234");
                assert_eq!(records[1].current().reference_bases(), "6789");
            }
            other => panic!("expected `Split`, got {other:?}"),
        }
    }
}
