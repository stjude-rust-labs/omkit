//! Rejection reasons for VCF records that fail transformation.

use std::fmt;

/// The reason a VCF record was rejected during transformation.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum RejectionReason {
    /// No chain interval covers this position.
    NoTarget,

    /// The variant's reference span crosses a chain alignment block
    /// boundary.
    IndelStraddlesMultipleIntervals,

    /// The lifted `REF` allele doesn't match the target FASTA and no
    /// swap is possible.
    MismatchedRefAllele,

    /// Reverse complement of alleles failed.
    CannotReverseComplement,

    /// Multiple chain intervals matched and ambiguous mapping is
    /// configured to reject.
    AmbiguousMapping {
        /// The number of candidate mappings found.
        candidates: usize,
    },
}

impl RejectionReason {
    /// Returns the `FILTER` string for this rejection reason.
    pub fn as_filter_str(&self) -> &'static str {
        match self {
            Self::NoTarget => "NoTarget",
            Self::IndelStraddlesMultipleIntervals => "IndelStraddlesMultipleIntervals",
            Self::MismatchedRefAllele => "MismatchedRefAllele",
            Self::CannotReverseComplement => "CannotReverseComplement",
            Self::AmbiguousMapping { .. } => "AmbiguousMapping",
        }
    }

    /// Returns a human-readable description for VCF header `##FILTER`
    /// lines.
    pub fn description(&self) -> &'static str {
        match self {
            Self::NoTarget => "No chain interval covers this position",
            Self::IndelStraddlesMultipleIntervals => {
                "Variant reference span crosses a chain alignment block boundary"
            }
            Self::MismatchedRefAllele => {
                "Lifted REF allele does not match the target reference and no swap is possible"
            }
            Self::CannotReverseComplement => "Reverse complement of alleles failed",
            Self::AmbiguousMapping { .. } => {
                "Multiple chain intervals matched and ambiguous mapping is configured to reject"
            }
        }
    }
}

impl fmt::Display for RejectionReason {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.description())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn all_variants() -> Vec<RejectionReason> {
        vec![
            RejectionReason::NoTarget,
            RejectionReason::IndelStraddlesMultipleIntervals,
            RejectionReason::MismatchedRefAllele,
            RejectionReason::CannotReverseComplement,
            RejectionReason::AmbiguousMapping { candidates: 3 },
        ]
    }

    #[test]
    fn as_filter_str_matches_variant_name() {
        assert_eq!(RejectionReason::NoTarget.as_filter_str(), "NoTarget");
        assert_eq!(
            RejectionReason::IndelStraddlesMultipleIntervals.as_filter_str(),
            "IndelStraddlesMultipleIntervals"
        );
        assert_eq!(
            RejectionReason::MismatchedRefAllele.as_filter_str(),
            "MismatchedRefAllele"
        );
        assert_eq!(
            RejectionReason::CannotReverseComplement.as_filter_str(),
            "CannotReverseComplement"
        );
        assert_eq!(
            RejectionReason::AmbiguousMapping { candidates: 5 }.as_filter_str(),
            "AmbiguousMapping"
        );
    }

    #[test]
    fn description_is_nonempty_for_all_variants() {
        for variant in all_variants() {
            assert!(!variant.description().is_empty(), "{variant:?}");
        }
    }

    #[test]
    fn display_uses_description() {
        for variant in all_variants() {
            assert_eq!(variant.to_string(), variant.description());
        }
    }
}
