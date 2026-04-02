//! A VCF record wrapper that carries accumulated transformation state.

use noodles::core::Position;
use noodles::vcf;

/// A VCF record with accumulated state from prior transformation stages.
///
/// As a record passes through the transformer pipeline, each stage may
/// modify the record and annotate it with metadata (e.g., whether a
/// `REF`/`ALT` swap was detected). This struct carries the current
/// (possibly modified) record, a snapshot of the original contig,
/// position, and `REF` allele, and accumulated state.
#[derive(Clone, Debug)]
pub struct VcfRecord {
    /// The current `noodles` VCF record, modified in place by transformers.
    current: vcf::variant::RecordBuf,

    /// The original contig name before any transformation.
    original_contig: String,

    /// The original 1-based VCF position before any transformation.
    original_position: Option<Position>,

    /// The original `REF` allele before any transformation.
    original_ref_bases: String,

    /// State accumulated by transformers.
    state: State,
}

/// State accumulated by transformers.
#[derive(Clone, Debug, Default)]
pub struct State {
    /// Whether the liftover mapped to the negative strand.
    pub negative_strand: bool,

    /// Whether a `REF`/`ALT` swap was detected.
    pub swap_detected: bool,
}

impl VcfRecord {
    /// Creates a new [`VcfRecord`] from a `noodles` record.
    pub fn new(current: vcf::variant::RecordBuf) -> Self {
        let original_contig = String::from(current.reference_sequence_name());
        let original_position = current.variant_start();
        let original_ref_bases = String::from(current.reference_bases());

        Self {
            current,
            original_contig,
            original_position,
            original_ref_bases,
            state: State::default(),
        }
    }

    /// Returns the original contig name before any transformation.
    pub fn original_contig(&self) -> &str {
        &self.original_contig
    }

    /// Returns the original 1-based VCF position before any
    /// transformation.
    pub fn original_position(&self) -> Option<Position> {
        self.original_position
    }

    /// Returns the original `REF` allele before any transformation.
    pub fn original_ref_bases(&self) -> &str {
        &self.original_ref_bases
    }

    /// Returns a reference to the current `noodles` record.
    pub fn current(&self) -> &vcf::variant::RecordBuf {
        &self.current
    }

    /// Returns a mutable reference to the current `noodles` record.
    pub fn current_mut(&mut self) -> &mut vcf::variant::RecordBuf {
        &mut self.current
    }

    /// Consumes the wrapper and returns the current `noodles` record.
    pub fn into_current(self) -> vcf::variant::RecordBuf {
        self.current
    }

    /// Returns a reference to the accumulated state.
    pub fn state(&self) -> &State {
        &self.state
    }

    /// Returns a mutable reference to the accumulated state.
    pub fn state_mut(&mut self) -> &mut State {
        &mut self.state
    }
}
