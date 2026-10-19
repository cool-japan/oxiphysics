// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Error type of the quantum-chemistry MD drivers.

/// Errors reported by the Born-Oppenheimer drivers of
/// [`quantum_chemistry_md`](crate::quantum_chemistry_md).
#[derive(Debug, Clone, PartialEq, thiserror::Error)]
pub enum QuantumChemistryMdError {
    /// Two inputs have inconsistent sizes.
    #[error("dimension mismatch for {what}: expected {expected}, found {found}")]
    DimensionMismatch {
        /// Which quantity was inconsistent.
        what: &'static str,
        /// Expected length.
        expected: usize,
        /// Actual length.
        found: usize,
    },
    /// A nuclear charge is not an element supported by the electronic-structure method.
    #[error("nuclear charge {charge} of atom {index} is not a supported element")]
    UnsupportedNuclearCharge {
        /// Index of the offending atom.
        index: usize,
        /// The nuclear charge that was given.
        charge: f64,
    },
    /// The electron count is odd; only closed-shell (restricted) references are supported.
    #[error("a closed-shell system with an even electron count is required, got {0} electrons")]
    OpenShell(i64),
    /// The electron count is negative or exceeds what the basis can hold.
    #[error("invalid electron count {electrons} (capacity {capacity})")]
    InvalidElectronCount {
        /// Number of electrons implied by the nuclear charges and the total charge.
        electrons: i64,
        /// Maximum number of electrons the orbital space can hold.
        capacity: usize,
    },
    /// The self-consistent-field iteration did not converge.
    #[error("the self-consistent-field iteration did not converge")]
    ScfNotConverged,
    /// The iterative eigensolver did not converge.
    #[error("the Kohn-Sham eigensolver did not converge (residual {residual:e})")]
    EigensolverNotConverged {
        /// Largest residual norm at the last iteration.
        residual: f64,
    },
    /// A scalar parameter is out of its valid range.
    #[error("invalid parameter {name}: {value}")]
    InvalidParameter {
        /// Name of the parameter.
        name: &'static str,
        /// The rejected value.
        value: f64,
    },
}
