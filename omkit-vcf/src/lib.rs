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

use record::VcfRecord;
use rejection::RejectionReason;

/// An error that can occur during transformation.
#[derive(Debug, thiserror::Error)]
pub enum Error {
    /// An I/O error.
    #[error(transparent)]
    Io(#[from] std::io::Error),
}

/// A [`Result`](std::result::Result) with the crate-local [`Error`].
pub type Result<T> = std::result::Result<T, Error>;

/// The result of a single transformer processing a record.
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
pub trait Transformer {
    /// Processes a single record and returns the outcome.
    fn transform(&self, record: VcfRecord) -> Result<Outcome>;
}

/// The results of processing a single input record through the pipeline.
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
}
