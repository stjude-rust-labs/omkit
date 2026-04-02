//! Validates VCF `REF` alleles against a FASTA reference.
//!
//! The [`RefValidator`] queries an indexed FASTA file to check whether
//! each record's `REF` allele matches the reference sequence at the
//! record's contig and position. It operates in one of two modes:
//!
//! In [`Mode::Reject`] mode (typically used pre-liftover), any mismatch
//! between the `REF` allele and the FASTA causes the record to be
//! rejected with [`RejectionReason::MismatchedRefAllele`].
//!
//! In [`Mode::DetectSwap`] mode (typically used post-liftover), a
//! mismatch triggers a check of the `ALT` allele against the reference.
//! If the `ALT` matches, a `REF`/`ALT` swap is flagged on the record's
//! state so that a downstream transformer can perform the actual swap.
//! If neither allele matches, the record is rejected.
//!
//! This transformer assumes that records are biallelic when running in
//! [`Mode::DetectSwap`]. A `SplitMultiAllelic` transformer must run
//! upstream to decompose multi-allelic records; if multiple valid `ALT`
//! alleles are present, this transformer will panic.
//!
//! Records with no position are rejected with
//! [`RejectionReason::MissingPosition`]. Records whose `REF` allele is
//! missing (`.`) or empty are passed through without validation.

use std::io::BufRead;
use std::io::Seek;
use std::sync::Mutex;
use std::sync::atomic::AtomicUsize;
use std::sync::atomic::Ordering;

use noodles::core::Position;
use noodles::core::Region;
use noodles::fasta::io::indexed_reader::IndexedReader;
use tracing::debug;
use tracing::trace;

use crate::MISSING;
use crate::Outcome;
use crate::Result;
use crate::SPANNING_DELETION;
use crate::Transformer;
use crate::record::VcfRecord;
use crate::rejection::RejectionReason;

/// Controls behavior when the `REF` allele does not match the reference.
#[derive(Clone, Debug)]
pub enum Mode {
    /// Reject the record outright.
    Reject,

    /// Detect whether the `REF` and an `ALT` allele have been swapped
    /// relative to the target reference. If a swap is detected, the
    /// record is accepted with `swap_detected` set on its state. If no
    /// allele matches, the record is rejected.
    DetectSwap,
}

/// Running statistics for [`RefValidator`].
#[derive(Debug, Default)]
struct Stats {
    /// Records whose `REF` allele matched the reference.
    accepted: AtomicUsize,

    /// Records where a `REF`/`ALT` swap was detected.
    swaps_detected: AtomicUsize,

    /// Records rejected because no allele matched the reference.
    rejected_mismatched: AtomicUsize,

    /// Records rejected because the contig was not in the reference
    /// FASTA.
    rejected_contig_not_found: AtomicUsize,
}

/// Validates VCF `REF` alleles against a FASTA reference sequence.
///
/// In [`Mode::Reject`] mode, records whose `REF` allele does not match the
/// FASTA are rejected. In [`Mode::DetectSwap`] mode (post-liftover), the
/// validator checks whether an `ALT` allele matches the reference instead,
/// flagging a swap on the record's state. If no allele matches, the record is
/// rejected.
pub struct RefValidator<R> {
    /// The indexed FASTA reader, wrapped in a [`Mutex`] because
    /// [`IndexedReader::query`] requires `&mut self`.
    reader: Mutex<IndexedReader<R>>,

    /// How to handle mismatches.
    mode: Mode,

    /// Running statistics.
    stats: Stats,
}

impl<R> RefValidator<R> {
    /// Creates a new `REF` allele validator.
    pub fn new(reader: IndexedReader<R>, mode: Mode) -> Self {
        Self {
            reader: Mutex::new(reader),
            mode,
            stats: Stats::default(),
        }
    }
}

impl<R> Transformer for RefValidator<R>
where
    R: BufRead + Seek + Send,
{
    fn transform(&self, mut record: VcfRecord) -> Result<Outcome> {
        let contig_name = String::from(record.current().reference_sequence_name());
        let vcf_pos: usize = match record.current().variant_start() {
            Some(pos) => pos.get(),
            None => {
                trace!(contig = %contig_name, "record has no position");
                return Ok(Outcome::Rejected(record, RejectionReason::MissingPosition));
            }
        };
        let ref_bases = record.current().reference_bases();

        // NOTE: unlike `CoordinateMapper`, which rejects records with
        // missing REF (since it needs REF to build the liftover interval),
        // this transformer passes them through because there is nothing
        // to validate.
        if ref_bases == MISSING || ref_bases.is_empty() {
            trace!(contig = %contig_name, position = vcf_pos, "REF is missing, skipping validation");
            self.stats.accepted.fetch_add(1, Ordering::Relaxed);
            return Ok(Outcome::Accepted(record));
        }

        let ref_bases = String::from(ref_bases);
        let ref_len = ref_bases.len();

        // SAFETY: `vcf_pos` came from `noodles::core::Position`, which wraps
        // `NonZero<usize>`, so it is always >= 1.
        let start = Position::try_from(vcf_pos).unwrap();

        // SAFETY: the early return above guarantees `ref_len` >= 1, so
        // `vcf_pos + ref_len - 1` >= `vcf_pos` >= 1.
        let end = Position::try_from(vcf_pos + ref_len - 1).unwrap();
        let region = Region::new(contig_name.as_str(), start..=end);

        let fasta_seq = {
            let mut reader = self.reader.lock().expect("FASTA reader lock poisoned");
            match reader.query(&region) {
                Ok(fasta_record) => fasta_record.sequence().as_ref().to_ascii_uppercase(),
                Err(e) => {
                    trace!(
                        contig = %contig_name,
                        position = vcf_pos,
                        error = %e,
                        "FASTA query failed"
                    );
                    self.stats
                        .rejected_contig_not_found
                        .fetch_add(1, Ordering::Relaxed);
                    return Ok(Outcome::Rejected(
                        record,
                        RejectionReason::ContigNotInReference,
                    ));
                }
            }
        };

        let ref_upper = ref_bases.as_bytes().to_ascii_uppercase();

        if fasta_seq == ref_upper {
            trace!(
                contig = %contig_name,
                position = vcf_pos,
                "REF allele matches reference"
            );
            self.stats.accepted.fetch_add(1, Ordering::Relaxed);
            return Ok(Outcome::Accepted(record));
        }

        match &self.mode {
            Mode::Reject => {
                trace!(
                    contig = %contig_name,
                    position = vcf_pos,
                    ref_allele = %ref_bases,
                    "REF allele does not match reference"
                );
                self.stats
                    .rejected_mismatched
                    .fetch_add(1, Ordering::Relaxed);
                Ok(Outcome::Rejected(
                    record,
                    RejectionReason::MismatchedRefAllele,
                ))
            }
            Mode::DetectSwap => {
                // NOTE: this assumes the record is biallelic. A
                // `SplitMultiAllelic` transformer must run upstream to
                // decompose multi-allelic records before this stage.
                let valid_alts = record
                    .current()
                    .alternate_bases()
                    .as_ref()
                    .iter()
                    .filter(|alt| *alt != MISSING && *alt != SPANNING_DELETION)
                    .collect::<Vec<_>>();

                assert!(
                    valid_alts.len() <= 1,
                    "`DetectSwap` requires biallelic records; found {} `ALT` alleles at \
                     `{contig_name}`:`{vcf_pos}`. Run `SplitMultiAllelic` upstream.",
                    valid_alts.len(),
                );

                // SAFETY: the assert above guarantees at most one valid
                // `ALT` allele, so `first()` is the only candidate.
                let alt = match valid_alts.first() {
                    Some(alt) => (*alt).clone(),
                    None => {
                        trace!(
                            contig = %contig_name,
                            position = vcf_pos,
                            ref_allele = %ref_bases,
                            "REF allele does not match reference and no ALT allele to swap"
                        );
                        self.stats
                            .rejected_mismatched
                            .fetch_add(1, Ordering::Relaxed);
                        return Ok(Outcome::Rejected(
                            record,
                            RejectionReason::MismatchedRefAllele,
                        ));
                    }
                };

                if alt.as_bytes().to_ascii_uppercase() == fasta_seq {
                    trace!(
                        contig = %contig_name,
                        position = vcf_pos,
                        ref_allele = %ref_bases,
                        "REF/ALT swap detected"
                    );
                    record.state_mut().swap_detected = true;
                    self.stats.swaps_detected.fetch_add(1, Ordering::Relaxed);
                    Ok(Outcome::Accepted(record))
                } else {
                    trace!(
                        contig = %contig_name,
                        position = vcf_pos,
                        ref_allele = %ref_bases,
                        "REF allele does not match reference and no swap possible"
                    );
                    self.stats
                        .rejected_mismatched
                        .fetch_add(1, Ordering::Relaxed);
                    Ok(Outcome::Rejected(
                        record,
                        RejectionReason::MismatchedRefAllele,
                    ))
                }
            }
        }
    }

    fn finish(&self) {
        let accepted = self.stats.accepted.load(Ordering::Relaxed);
        let swaps = self.stats.swaps_detected.load(Ordering::Relaxed);
        let mismatched = self.stats.rejected_mismatched.load(Ordering::Relaxed);
        let contig_not_found = self.stats.rejected_contig_not_found.load(Ordering::Relaxed);
        let total = accepted + swaps + mismatched + contig_not_found;

        debug!(
            total,
            accepted,
            swaps_detected = swaps,
            rejected_mismatched = mismatched,
            rejected_contig_not_found = contig_not_found,
            "`RefValidator` finished"
        );
    }
}

#[cfg(test)]
mod tests {
    use std::io::BufReader;
    use std::io::Cursor;

    use noodles::core::Position;
    use noodles::fasta;
    use noodles::fasta::io::indexed_reader::IndexedReader;
    use noodles::vcf::variant::record_buf::AlternateBases;

    use super::*;
    use crate::Outcome;
    use crate::Transformer;
    use crate::record::VcfRecord;
    use crate::rejection::RejectionReason;

    /// Builds an [`IndexedReader`] backed by an in-memory FASTA with a single
    /// contig `sq0` whose sequence is `ACGTACGT`.
    fn test_reader() -> IndexedReader<BufReader<Cursor<Vec<u8>>>> {
        let fasta_data = b">sq0\nACGTACGT\n";
        let index = fasta::fai::Index::from(vec![fasta::fai::Record::new("sq0", 8, 5, 8, 9)]);
        IndexedReader::new(BufReader::new(Cursor::new(fasta_data.to_vec())), index)
    }

    /// Builds a [`VcfRecord`] with the given contig, position, `REF`, and `ALT`
    /// alleles.
    fn test_record(contig: &str, pos: usize, ref_bases: &str, alts: Vec<&str>) -> VcfRecord {
        let mut inner = noodles::vcf::variant::RecordBuf::default();
        *inner.reference_sequence_name_mut() = String::from(contig);
        *inner.variant_start_mut() = Some(Position::try_from(pos).unwrap());
        *inner.reference_bases_mut() = String::from(ref_bases);
        *inner.alternate_bases_mut() =
            AlternateBases::from(alts.into_iter().map(String::from).collect::<Vec<_>>());
        VcfRecord::new(inner)
    }

    #[test]
    fn reject_mode_accepts_matching_ref() {
        let validator = RefValidator::new(test_reader(), Mode::Reject);
        let record = test_record("sq0", 1, "A", vec!["T"]);
        let outcome = validator.transform(record).unwrap();
        assert!(matches!(outcome, Outcome::Accepted(_)));
    }

    #[test]
    fn reject_mode_rejects_mismatched_ref() {
        let validator = RefValidator::new(test_reader(), Mode::Reject);
        let record = test_record("sq0", 1, "T", vec!["A"]);
        let outcome = validator.transform(record).unwrap();
        assert!(matches!(
            outcome,
            Outcome::Rejected(_, RejectionReason::MismatchedRefAllele)
        ));
    }

    #[test]
    fn detect_swap_accepts_matching_ref() {
        let validator = RefValidator::new(test_reader(), Mode::DetectSwap);
        let record = test_record("sq0", 1, "A", vec!["T"]);
        let outcome = validator.transform(record).unwrap();
        match outcome {
            Outcome::Accepted(r) => assert!(!r.state().swap_detected),
            other => panic!("expected `Accepted`, got {other:?}"),
        }
    }

    #[test]
    fn detect_swap_flags_swap() {
        let validator = RefValidator::new(test_reader(), Mode::DetectSwap);
        // FASTA has `A` at position 1; REF=T, ALT=A → swap detected.
        let record = test_record("sq0", 1, "T", vec!["A"]);
        let outcome = validator.transform(record).unwrap();
        match outcome {
            Outcome::Accepted(r) => assert!(r.state().swap_detected),
            other => panic!("expected `Accepted` with swap, got {other:?}"),
        }
    }

    #[test]
    fn detect_swap_rejects_when_no_allele_matches() {
        let validator = RefValidator::new(test_reader(), Mode::DetectSwap);
        // FASTA has `A` at position 1; REF=T, ALT=G → no match.
        let record = test_record("sq0", 1, "T", vec!["G"]);
        let outcome = validator.transform(record).unwrap();
        assert!(matches!(
            outcome,
            Outcome::Rejected(_, RejectionReason::MismatchedRefAllele)
        ));
    }

    #[test]
    fn unknown_contig_is_rejected() {
        let validator = RefValidator::new(test_reader(), Mode::Reject);
        let record = test_record("chrX", 1, "A", vec!["T"]);
        let outcome = validator.transform(record).unwrap();
        assert!(matches!(
            outcome,
            Outcome::Rejected(_, RejectionReason::ContigNotInReference)
        ));
    }

    #[test]
    #[should_panic(expected = "`DetectSwap` requires biallelic records")]
    fn detect_swap_panics_on_multi_allelic() {
        let validator = RefValidator::new(test_reader(), Mode::DetectSwap);
        // FASTA has `A` at position 1; REF=T, ALT=G,C → multi-allelic.
        let record = test_record("sq0", 1, "T", vec!["G", "C"]);
        let _ = validator.transform(record);
    }
}
