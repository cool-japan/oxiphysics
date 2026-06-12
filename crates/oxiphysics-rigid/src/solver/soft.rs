//! Soft-constraint (frequency / damping-ratio) parameterization for the
//! TGS-Soft contact solver.
//!
//! Ported from `oxiphysics_constraints::tgs_solver::types::SoftParams`.
//!
//! `oxiphysics-constraints` already depends on `oxiphysics-rigid`, so
//! `oxiphysics-rigid` cannot depend back on it without introducing a crate
//! dependency cycle. The ~12-line Catto / Box2D-v3 derivation is therefore
//! reproduced here verbatim; keep it numerically in sync with the constraints
//! crate (the soft-params formula is identical in both places).
//!
//! # Derivation summary
//!
//! Given natural frequency ω = 2π·hz and damping ratio ζ, with sub-step dt = h:
//!
//! ```text
//! a₁ = 2ζ + h·ω
//! a₂ = h·ω·a₁
//! a₃ = 1 / (1 + a₂)
//!
//! bias_rate     = ω / a₁
//! mass_scale    = a₂ · a₃
//! impulse_scale = a₃
//! ```
//!
//! See: Erin Catto, "Soft Constraints", Box2D v3 source; also
//! Baumgarte (1972), Ascher et al. (1995) for the classical connection.

/// Soft-constraint coefficients derived from a natural frequency and damping
/// ratio (Catto / Box2D v3).
///
/// The coefficients enter the accumulated-impulse update of a constraint row as
///
/// ```text
/// impulse = -(mass_scale · (Jv + bias_rate · c) + impulse_scale · λ) / eff_mass
/// ```
///
/// where `c` is the (signed) position error fed into the soft bias.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SoftParams {
    /// Coefficient that multiplies the position error to produce the soft bias
    /// contribution.
    pub bias_rate: f64,
    /// Factor applied to the velocity error (softens the effective mass).
    pub mass_scale: f64,
    /// Factor applied to the accumulated impulse `λ` (damps it each iteration).
    pub impulse_scale: f64,
}

impl SoftParams {
    /// Create soft params from a frequency (Hz) and damping ratio ζ, given the
    /// current sub-step size `h` (seconds).
    ///
    /// If `hz <= 0.0`, returns [`SoftParams::rigid`] immediately.
    ///
    /// # Arguments
    ///
    /// * `hz`   — Constraint natural frequency in Hz. Typical range: 1–120 Hz.
    /// * `zeta` — Damping ratio ζ. ζ < 1 = under-damped, ζ = 1 = critically
    ///   damped, ζ > 1 = over-damped.
    /// * `h`    — Sub-step size in seconds (e.g. `dt / substeps`).
    pub fn from_frequency(hz: f64, zeta: f64, h: f64) -> Self {
        if hz <= 0.0 {
            return Self::rigid();
        }
        let omega = std::f64::consts::TAU * hz;
        let a1 = 2.0 * zeta + h * omega;
        let a2 = h * omega * a1;
        let a3 = 1.0 / (1.0 + a2);
        Self {
            bias_rate: omega / a1,
            mass_scale: a2 * a3,
            impulse_scale: a3,
        }
    }

    /// Rigid (infinitely stiff) parameters — equivalent to a standard
    /// Baumgarte-stabilised TGS row with no softening: `bias_rate = 0`,
    /// `mass_scale = 1`, `impulse_scale = 0`.
    pub fn rigid() -> Self {
        Self {
            bias_rate: 0.0,
            mass_scale: 1.0,
            impulse_scale: 0.0,
        }
    }

    /// Returns `true` if these parameters represent a rigid (non-soft)
    /// constraint.
    pub fn is_rigid(&self) -> bool {
        self.mass_scale >= 1.0 - 1e-12 && self.impulse_scale <= 1e-12
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn hz_zero_is_rigid() {
        let soft = SoftParams::from_frequency(0.0, 1.0, 1.0 / 60.0);
        assert_eq!(soft, SoftParams::rigid());
        assert!(soft.is_rigid());
    }

    #[test]
    fn rigid_constants() {
        let r = SoftParams::rigid();
        assert_eq!(r.bias_rate, 0.0);
        assert_eq!(r.mass_scale, 1.0);
        assert_eq!(r.impulse_scale, 0.0);
    }

    #[test]
    fn from_frequency_matches_catto_formula() {
        let hz = 30.0;
        let zeta = 1.0;
        let h = 1.0 / 240.0;
        let omega = std::f64::consts::TAU * hz;
        let a1 = 2.0 * zeta + h * omega;
        let a2 = h * omega * a1;
        let a3 = 1.0 / (1.0 + a2);
        let soft = SoftParams::from_frequency(hz, zeta, h);
        assert!((soft.bias_rate - omega / a1).abs() < 1e-12);
        assert!((soft.mass_scale - a2 * a3).abs() < 1e-12);
        assert!((soft.impulse_scale - a3).abs() < 1e-12);
        assert!(!soft.is_rigid());
    }

    #[test]
    fn coefficients_in_unit_range() {
        // mass_scale and impulse_scale are convex combinations bounded to [0, 1].
        for &hz in &[1.0, 10.0, 30.0, 60.0, 120.0] {
            for &zeta in &[0.5, 1.0, 2.0] {
                let soft = SoftParams::from_frequency(hz, zeta, 1.0 / 120.0);
                assert!((0.0..=1.0).contains(&soft.mass_scale), "mass_scale oob");
                assert!(
                    (0.0..=1.0).contains(&soft.impulse_scale),
                    "impulse_scale oob"
                );
                assert!(soft.bias_rate > 0.0, "bias_rate must be positive");
            }
        }
    }

    #[test]
    fn stiffer_frequency_increases_bias_rate() {
        let h = 1.0 / 240.0;
        let lo = SoftParams::from_frequency(10.0, 1.0, h);
        let hi = SoftParams::from_frequency(60.0, 1.0, h);
        assert!(hi.bias_rate > lo.bias_rate);
    }
}
