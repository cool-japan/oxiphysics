// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Error type for fallible articulated-body computations.

/// Errors returned by fallible articulated-body algorithms (e.g. the
/// operational-space inertia in [`crate::osc`]).
#[derive(Debug, Clone, PartialEq)]
pub enum ArticulatedError {
    /// The joint-space mass matrix was not positive-definite: the Cholesky
    /// factorisation met a non-positive pivot.
    NotPositiveDefinite {
        /// Zero-based pivot index where factorisation failed.
        index: usize,
    },
    /// The task-space inertia inverse `J·M⁻¹·Jᵀ` is singular (rank-deficient),
    /// e.g. a 6-row task Jacobian on a robot with fewer than 6 DoF or at a
    /// kinematic singularity. Use a damped variant to obtain a regularised
    /// result.
    SingularTaskInertia,
    /// An input had the wrong length.
    DimensionMismatch {
        /// Which input was malformed (e.g. `"jacobian rows"`).
        what: &'static str,
        /// Expected length.
        expected: usize,
        /// Actual length supplied.
        got: usize,
    },
    /// A numeric parameter was outside its valid domain.
    InvalidParameter {
        /// Name of the offending parameter.
        what: &'static str,
        /// The supplied value.
        value: f64,
    },
}

impl std::fmt::Display for ArticulatedError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::NotPositiveDefinite { index } => write!(
                f,
                "mass matrix is not positive-definite at pivot index {index}"
            ),
            Self::SingularTaskInertia => write!(
                f,
                "task-space inertia J·M⁻¹·Jᵀ is singular (rank-deficient task Jacobian)"
            ),
            Self::DimensionMismatch {
                what,
                expected,
                got,
            } => write!(
                f,
                "dimension mismatch for {what}: expected {expected}, got {got}"
            ),
            Self::InvalidParameter { what, value } => {
                write!(f, "invalid parameter {what} = {value}")
            }
        }
    }
}

impl std::error::Error for ArticulatedError {}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn display_messages_name_the_problem() {
        let e = ArticulatedError::NotPositiveDefinite { index: 3 };
        assert!(e.to_string().contains("pivot index 3"));
        let e = ArticulatedError::DimensionMismatch {
            what: "jacobian rows",
            expected: 6,
            got: 5,
        };
        let msg = e.to_string();
        assert!(msg.contains("jacobian rows") && msg.contains('6') && msg.contains('5'));
        assert!(
            ArticulatedError::SingularTaskInertia
                .to_string()
                .contains("singular")
        );
    }
}
