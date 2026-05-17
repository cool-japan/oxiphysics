// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Conservation bookkeeping for cross-domain coupling.
//!
//! These utilities aggregate `CouplingReport`s across multiple pairs to
//! track global momentum exchange and diagnose conservation violations.

use crate::coupling::CouplingReport;

/// Sum the total momentum exchanged from A→B across all reports.
///
/// Each report's `momentum_a_to_b` field represents `force × dt` transferred
/// from domain A to domain B during one step.  Summing across all reports
/// gives the net momentum deposited into the B sides.
pub fn total_exchanged_momentum(reports: &[CouplingReport]) -> [f64; 3] {
    reports.iter().fold([0.0_f64; 3], |mut acc, r| {
        acc[0] += r.momentum_a_to_b[0];
        acc[1] += r.momentum_a_to_b[1];
        acc[2] += r.momentum_a_to_b[2];
        acc
    })
}

/// Compute the fractional momentum drift relative to a reference magnitude.
///
/// `delta` is the difference between the expected and actual exchanged
/// momentum vector.  `reference` is the expected momentum vector.
///
/// Returns `|delta| / |reference|`, or `0.0` if the reference is zero.
///
/// This is useful for checking that action–reaction is satisfied to within
/// numerical tolerance: `|Σ momentum_a_to_b + Σ momentum_b_to_a| / |reference|`.
pub fn momentum_drift_fraction(delta: [f64; 3], reference: [f64; 3]) -> f64 {
    let ref_mag =
        (reference[0] * reference[0] + reference[1] * reference[1] + reference[2] * reference[2])
            .sqrt();
    if ref_mag < f64::EPSILON {
        return 0.0;
    }
    let delta_mag = (delta[0] * delta[0] + delta[1] * delta[1] + delta[2] * delta[2]).sqrt();
    delta_mag / ref_mag
}

/// Compute the RMS position error across all reports.
///
/// Useful for monitoring whether interface mismatch is growing over many
/// steps.
pub fn mean_rms_position_error(reports: &[CouplingReport]) -> f64 {
    if reports.is_empty() {
        return 0.0;
    }
    let sum: f64 = reports.iter().map(|r| r.rms_position_error).sum();
    sum / reports.len() as f64
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn total_exchanged_momentum_zero_for_empty() {
        let result = total_exchanged_momentum(&[]);
        assert_eq!(result, [0.0, 0.0, 0.0]);
    }

    #[test]
    fn momentum_drift_fraction_zero_reference() {
        let drift = momentum_drift_fraction([1.0, 0.0, 0.0], [0.0, 0.0, 0.0]);
        assert!((drift - 0.0).abs() < f64::EPSILON);
    }
}
