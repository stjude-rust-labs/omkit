//! A library for VCF transformations including liftover and normalization.
//!
//! This crate provides a composable transformer pipeline for VCF records.
//! Each transformer is an independent entity that can accept, split, or
//! reject records. Commands compose transformers into pipelines to build
//! tools like liftover and normalization.
//!
//! A [`Pipeline`] drives records through an ordered sequence of
//! [`Transformer`] implementations. Each transformer returns an
//! [`Outcome`] indicating whether the record was accepted (possibly
//! modified), split into multiple records, or rejected with a reason.

pub mod record;
pub mod rejection;
pub mod transformer;

use record::VcfRecord;
use rejection::RejectionReason;

/// The VCF missing value indicator (`.`).
pub const MISSING: &str = ".";

/// The VCF spanning deletion allele (`*`).
pub const SPANNING_DELETION: &str = "*";

/// An error that can occur during transformation.
#[derive(Debug, thiserror::Error)]
pub enum Error {
    /// An I/O error.
    #[error(transparent)]
    Io(#[from] std::io::Error),

    /// A required extension was not found in the record's extensions.
    #[error("`{consumer}` requires `{producer}` to run first")]
    MissingExtension {
        /// The name of the transformer that requires the extension.
        consumer: &'static str,

        /// The name of the transformer that produces the extension.
        producer: &'static str,
    },
}

/// A [`Result`](std::result::Result) with the crate-local [`Error`].
pub type Result<T> = std::result::Result<T, Error>;

/// The result of a single transformer processing a record.
#[derive(Debug)]
pub enum Outcome {
    /// The record was accepted, possibly modified.
    Accepted(VcfRecord),

    /// The record was split into multiple records. Subsequent
    /// transformers run independently on each resulting record.
    Split(Vec<VcfRecord>),

    /// The record was rejected with a reason. Processing stops for this
    /// record and it is written to the rejected output.
    Rejected(VcfRecord, RejectionReason),
}

/// A transformer that processes a single VCF record.
///
/// Implementations must be `Send + Sync` to support future multithreaded
/// pipeline execution. Statistics should use atomic types for thread safety.
/// The default [`finish()`](Transformer::finish) implementation is a no-op;
/// override it to log summary statistics after all records have been
/// processed.
pub trait Transformer: Send + Sync {
    /// Returns the name of this transformer for diagnostics.
    ///
    /// The default implementation uses the type name.
    fn name(&self) -> &'static str {
        std::any::type_name::<Self>()
    }

    /// Processes a single record and returns the outcome.
    fn transform(&self, record: VcfRecord) -> Result<Outcome>;

    /// Called after all records have been processed.
    ///
    /// Implementations should log summary statistics here. The default
    /// implementation is a no-op.
    fn finish(&self) {}
}

/// The results of processing a single input record through the pipeline.
#[derive(Debug)]
pub struct PipelineResults {
    /// Records that made it through all transformers successfully.
    pub accepted: Vec<VcfRecord>,

    /// Records that were rejected, paired with their rejection reasons.
    pub rejected: Vec<(VcfRecord, RejectionReason)>,
}

/// An ordered sequence of transformers that processes VCF records.
pub struct Pipeline {
    /// The ordered transformers.
    transformers: Vec<Box<dyn Transformer>>,
}

impl Pipeline {
    /// Creates a new pipeline from the given transformers.
    pub fn new(transformers: Vec<Box<dyn Transformer>>) -> Self {
        Self { transformers }
    }

    /// Processes a single record through all transformers.
    ///
    /// Records that are [`Outcome::Split`] are fanned out and each copy
    /// continues independently through the remaining transformers.
    /// Records that are [`Outcome::Rejected`] stop processing and are
    /// collected in the rejected set.
    pub fn process(&self, record: VcfRecord) -> Result<PipelineResults> {
        let mut current = vec![record];
        let mut rejected = Vec::new();

        for transformer in &self.transformers {
            let mut next = Vec::new();

            for rec in current {
                match transformer.transform(rec)? {
                    Outcome::Accepted(r) => next.push(r),
                    Outcome::Split(rs) => next.extend(rs),
                    Outcome::Rejected(r, reason) => rejected.push((r, reason)),
                }
            }

            current = next;
        }

        Ok(PipelineResults {
            accepted: current,
            rejected,
        })
    }

    /// Processes all records and calls [`Transformer::finish()`] on each
    /// transformer when complete.
    ///
    /// [`Transformer::finish()`] is called even if a record produces an
    /// error, ensuring that statistics are always logged.
    pub fn process_all(
        &self,
        records: impl IntoIterator<Item = VcfRecord>,
    ) -> Result<PipelineResults> {
        let mut all_accepted = Vec::new();
        let mut all_rejected = Vec::new();
        let mut error = None;

        for record in records {
            match self.process(record) {
                Ok(results) => {
                    all_accepted.extend(results.accepted);
                    all_rejected.extend(results.rejected);
                }
                Err(e) => {
                    error = Some(e);
                    break;
                }
            }
        }

        self.finish();

        if let Some(e) = error {
            return Err(e);
        }

        Ok(PipelineResults {
            accepted: all_accepted,
            rejected: all_rejected,
        })
    }

    /// Calls [`Transformer::finish()`] on each transformer in order.
    pub fn finish(&self) {
        for transformer in &self.transformers {
            transformer.finish();
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    struct PassThrough;

    impl Transformer for PassThrough {
        fn transform(&self, record: VcfRecord) -> Result<Outcome> {
            Ok(Outcome::Accepted(record))
        }
    }

    struct RejectAll;

    impl Transformer for RejectAll {
        fn transform(&self, record: VcfRecord) -> Result<Outcome> {
            Ok(Outcome::Rejected(record, RejectionReason::NoTarget))
        }
    }

    struct SplitTwo;

    impl Transformer for SplitTwo {
        fn transform(&self, record: VcfRecord) -> Result<Outcome> {
            Ok(Outcome::Split(vec![record.clone(), record]))
        }
    }

    #[test]
    fn pipeline_passthrough() {
        let pipeline = Pipeline::new(vec![Box::new(PassThrough)]);
        let record = VcfRecord::new(Default::default());
        let results = pipeline.process(record).unwrap();
        assert_eq!(results.accepted.len(), 1);
        assert!(results.rejected.is_empty());
    }

    #[test]
    fn pipeline_reject_stops_processing() {
        let pipeline = Pipeline::new(vec![Box::new(RejectAll), Box::new(PassThrough)]);
        let record = VcfRecord::new(Default::default());
        let results = pipeline.process(record).unwrap();
        assert!(results.accepted.is_empty());
        assert_eq!(results.rejected.len(), 1);
    }

    #[test]
    fn pipeline_split_fans_out() {
        let pipeline = Pipeline::new(vec![Box::new(SplitTwo), Box::new(PassThrough)]);
        let record = VcfRecord::new(Default::default());
        let results = pipeline.process(record).unwrap();
        assert_eq!(results.accepted.len(), 2);
        assert!(results.rejected.is_empty());
    }

    #[test]
    fn pipeline_split_then_reject() {
        let pipeline = Pipeline::new(vec![Box::new(SplitTwo), Box::new(RejectAll)]);
        let record = VcfRecord::new(Default::default());
        let results = pipeline.process(record).unwrap();
        assert!(results.accepted.is_empty());
        assert_eq!(results.rejected.len(), 2);
    }

    #[test]
    fn process_all_aggregates_and_finishes() {
        let pipeline = Pipeline::new(vec![Box::new(PassThrough)]);
        let records = vec![
            VcfRecord::new(Default::default()),
            VcfRecord::new(Default::default()),
            VcfRecord::new(Default::default()),
        ];
        let results = pipeline.process_all(records).unwrap();
        assert_eq!(results.accepted.len(), 3);
        assert!(results.rejected.is_empty());
    }
}
