//! A VCF record wrapper that carries accumulated transformation state.

use noodles::vcf;
use omics::coordinate::position::base::Position;

/// A VCF record with accumulated state from prior transformation stages.
///
/// As a record passes through the transformer pipeline, each stage may
/// annotate it with metadata (e.g., the lifted coordinate, whether a
/// REF/ALT swap was detected). This struct carries both the underlying
/// VCF record and that accumulated state.
#[derive(Clone, Debug)]
pub struct VcfRecord {
    /// The underlying `noodles` VCF record.
    inner: vcf::variant::RecordBuf,

    /// The state accumulated by transformers.
    state: State,
}

/// Accumulated state from prior transformation stages.
#[derive(Clone, Debug, Default)]
pub struct State {
    /// Whether the liftover mapped to the negative strand.
    ///
    /// Set by `CoordinateMapper`, read by `StrandFlipper`.
    pub negative_strand: bool,

    /// The number of chain alignment blocks the variant's reference
    /// span covers.
    ///
    /// Set by `CoordinateMapper`, read by `IndelStraddleDetector`.
    ///
    /// A value > `1` means the indel straddles block boundaries.
    pub liftover_block_count: usize,

    /// Whether a REF/ALT swap was detected.
    ///
    /// Set by `RefValidator`, read by `SwapTransformer` and
    /// `GenotypeInverter`.
    pub swap_detected: bool,

    /// The original contig name before liftover.
    ///
    /// Set by `CoordinateMapper`, read by `MetadataUpdater` and
    /// rejected record annotation.
    pub original_contig: Option<String>,

    /// The original 1-based VCF position before liftover.
    ///
    /// Set by `CoordinateMapper`, read by `MetadataUpdater` and
    /// rejected record annotation.
    pub original_position: Option<Position>,

    /// The original alleles before any transformation.
    ///
    /// Set by `CoordinateMapper`, read by `MetadataUpdater` and
    /// rejected record annotation.
    pub original_alleles: Option<Vec<String>>,
}

impl VcfRecord {
    /// Creates a new [`VcfRecord`] from a `noodles` record.
    pub fn new(inner: vcf::variant::RecordBuf) -> Self {
        Self {
            inner,
            state: State::default(),
        }
    }

    /// Returns a reference to the underlying `noodles` record.
    pub fn inner(&self) -> &vcf::variant::RecordBuf {
        &self.inner
    }

    /// Returns a mutable reference to the underlying `noodles` record.
    pub fn inner_mut(&mut self) -> &mut vcf::variant::RecordBuf {
        &mut self.inner
    }

    /// Consumes the wrapper and returns the underlying `noodles` record.
    pub fn into_inner(self) -> vcf::variant::RecordBuf {
        self.inner
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
