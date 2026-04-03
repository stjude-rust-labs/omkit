//! A VCF record wrapper that carries accumulated transformation state.

pub mod extensions;

use noodles::vcf;

use extensions::Extensions;

/// A VCF record with accumulated state from prior transformation stages.
///
/// As a record passes through the transformer pipeline, each stage may
/// modify the record and annotate it with metadata (e.g., whether a
/// `REF`/`ALT` swap was detected). This struct carries the current
/// (possibly modified) record and type-erased [`Extensions`] for
/// inter-transformer communication.
#[derive(Clone, Debug)]
pub struct VcfRecord {
    /// The current `noodles` VCF record, modified in place by transformers.
    current: vcf::variant::RecordBuf,

    /// The type-erased extensions set by transformers.
    extensions: Extensions,
}

impl VcfRecord {
    /// Creates a new [`VcfRecord`] from a `noodles` record.
    pub fn new(current: vcf::variant::RecordBuf) -> Self {
        Self {
            current,
            extensions: Extensions::default(),
        }
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

    /// Returns a reference to the extensions map.
    pub fn extensions(&self) -> &Extensions {
        &self.extensions
    }

    /// Returns a mutable reference to the extensions map.
    pub fn extensions_mut(&mut self) -> &mut Extensions {
        &mut self.extensions
    }
}
