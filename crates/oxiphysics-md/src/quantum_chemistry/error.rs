// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Error type of the quantum-chemistry toolkit.

use thiserror::Error;

/// Failure of a quantum-chemistry calculation (post-HF correlation, embedding, ...).
#[derive(Debug, Clone, PartialEq, Error)]
#[non_exhaustive]
pub enum QuantumChemistryError {
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
    /// A required input (integrals, dipoles, basis functions, ...) is missing.
    #[error("missing required input: {0}")]
    MissingData(&'static str),
    /// The reference SCF calculation did not converge.
    #[error("the reference SCF calculation did not converge")]
    ScfNotConverged,
    /// The method requires a closed-shell reference (even electron count).
    #[error("a closed-shell reference with an even electron count is required, got {0} electrons")]
    OpenShellReference(usize),
    /// An orbital-energy denominator ε_a + ε_b − ε_i − ε_j (or ε_a − ε_i) is not
    /// positive, so perturbation theory is undefined (non-aufbau reference).
    #[error("non-positive orbital-energy denominator {0}")]
    NonPositiveDenominator(f64),
    /// An external point charge is non-finite or coincides with a nucleus.
    #[error("external point charge {index} is non-finite or coincides with a nucleus")]
    InvalidExternalCharge {
        /// Index of the offending charge.
        index: usize,
    },
    /// Integrals for this angular momentum are not available for the requested property.
    #[error("integrals for angular momentum l={0} are not available for this property")]
    UnsupportedAngularMomentum(u8),
    /// The overlap matrix is not positive definite (linearly dependent basis).
    #[error("the overlap matrix is not positive definite (linearly dependent basis)")]
    SingularOverlap,
    /// Two nuclei occupy the same position, so no molecular grid can be built.
    #[error("nuclei {first} and {second} coincide")]
    CoincidentNuclei {
        /// Index of the first nucleus.
        first: usize,
        /// Index of the second nucleus.
        second: usize,
    },
}
