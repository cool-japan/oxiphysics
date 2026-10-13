//! Sampled unit-step response of a [`PolesZeros`] model.
//!
//! The pole/zero/gain description `G(s) = K prod(s - z_i) / prod(s - p_i)` is
//! expanded into numerator and denominator polynomials, realised in
//! controllable canonical form and simulated with the exact zero-order hold
//! (see `crate::control_systems::response`).  Zeros, the gain, repeated
//! poles, integrators and complex pole pairs are therefore all honoured
//! exactly at the sample instants.

use super::PolesZeros;
use crate::control_systems::response::{C64, SisoRealization, polynomial_from_roots};

impl PolesZeros {
    /// Sampled unit-step response of `G(s) = K prod(s - z_i) / prod(s - p_i)`.
    ///
    /// Returns `steps` pairs `(t_k, y_k)` with `t_k = k * dt`.  Despite the
    /// historical name the samples are not approximations: the model is
    /// realised in state space and discretised with the exact zero-order hold,
    /// so `y_k` equals the true step response at `t_k` up to rounding for every
    /// pole/zero pattern (distinct, repeated, complex, at the origin, in the
    /// right half plane).  For a stable model the samples converge to `G(0)`.
    ///
    /// Each entry of `poles` and `zeros` is one root, so a complex pair must be
    /// listed as both `(re, im)` and `(re, -im)`.  If the lists are not closed
    /// under conjugation the transfer function has complex coefficients and
    /// the real part of its response is returned.
    ///
    /// An empty vector is returned when no sampled step response exists: the
    /// model is improper (more zeros than poles, so the response contains
    /// Dirac impulses), `dt` is not a positive finite number, or the gain or a
    /// root is not finite.
    pub fn step_response_approx(&self, dt: f64, steps: usize) -> Vec<(f64, f64)> {
        let roots_to_complex = |roots: &[(f64, f64)]| -> Vec<C64> {
            roots.iter().map(|&(re, im)| C64::new(re, im)).collect()
        };
        let numerator =
            polynomial_from_roots(C64::new(self.gain, 0.0), &roots_to_complex(&self.zeros));
        let denominator = polynomial_from_roots(C64::new(1.0, 0.0), &roots_to_complex(&self.poles));
        let Some(realization) = SisoRealization::new(&numerator, &denominator) else {
            return Vec::new();
        };
        let Some(outputs) = realization.step_samples(dt, steps) else {
            return Vec::new();
        };
        outputs
            .into_iter()
            .enumerate()
            .map(|(k, y)| (k as f64 * dt, y))
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Largest deviation of the sampled response from `reference(t)`.
    fn max_deviation(response: &[(f64, f64)], reference: impl Fn(f64) -> f64) -> f64 {
        response
            .iter()
            .map(|&(t, y)| (y - reference(t)).abs())
            .fold(0.0_f64, f64::max)
    }

    #[test]
    fn single_real_pole_has_residue_scaled_by_the_pole() {
        // K / (s + a) with K = 2, a = 3: y(t) = (K / a)(1 - e^(-a t)).  The old
        // code returned K (1 - e^(-a t)), i.e. it ignored the 1/a residue factor.
        let pz = PolesZeros::new(2.0, vec![], vec![(-3.0, 0.0)]);
        let response = pz.step_response_approx(0.01, 400);
        assert_eq!(response.len(), 400);
        let dev = max_deviation(&response, |t| (2.0 / 3.0) * (1.0 - (-3.0 * t).exp()));
        assert!(dev < 1e-10, "max deviation {dev}");
    }

    #[test]
    fn time_stamps_are_multiples_of_dt() {
        let pz = PolesZeros::new(1.0, vec![], vec![(-1.0, 0.0)]);
        let response = pz.step_response_approx(0.25, 5);
        for (k, &(t, _)) in response.iter().enumerate() {
            assert!((t - 0.25 * k as f64).abs() < 1e-15);
        }
        assert_eq!(response[0].1, 0.0);
    }

    #[test]
    fn complex_pole_pair_matches_analytic_response_and_dc_gain() {
        // 5 / ((s + 1)^2 + 4):  y(t) = 1 - e^-t (cos 2t + sin(2t) / 2),  G(0) = 1.
        let pz = PolesZeros::new(5.0, vec![], vec![(-1.0, 2.0), (-1.0, -2.0)]);
        let response = pz.step_response_approx(0.01, 800);
        let dev = max_deviation(&response, |t| {
            1.0 - (-t).exp() * ((2.0 * t).cos() + 0.5 * (2.0 * t).sin())
        });
        assert!(dev < 1e-10, "max deviation {dev}");
        let final_value = response[response.len() - 1].1;
        assert!(
            (final_value - pz.dc_gain()).abs() < 1e-3,
            "final {final_value} vs dc gain {}",
            pz.dc_gain()
        );
    }

    #[test]
    fn minimum_phase_zero_shapes_the_response() {
        // (s + 3) / ((s + 1)(s + 2)): y(t) = 3/2 - 2 e^-t + 1/2 e^-2t -> G(0) = 3/2.
        let pz = PolesZeros::new(1.0, vec![(-3.0, 0.0)], vec![(-1.0, 0.0), (-2.0, 0.0)]);
        let response = pz.step_response_approx(0.02, 500);
        let dev = max_deviation(&response, |t| {
            1.5 - 2.0 * (-t).exp() + 0.5 * (-2.0 * t).exp()
        });
        assert!(dev < 1e-10, "max deviation {dev}");
        let final_value = response[response.len() - 1].1;
        assert!((final_value - pz.dc_gain()).abs() < 1e-3);
    }

    #[test]
    fn right_half_plane_zero_gives_initial_undershoot() {
        // K = -1, zero at +1: G = (1 - s) / ((s + 1)(s + 2)),
        // y(t) = 1/2 - 2 e^-t + 3/2 e^-2t with a negative initial slope.
        let pz = PolesZeros::new(-1.0, vec![(1.0, 0.0)], vec![(-1.0, 0.0), (-2.0, 0.0)]);
        let response = pz.step_response_approx(0.01, 600);
        let dev = max_deviation(&response, |t| {
            0.5 - 2.0 * (-t).exp() + 1.5 * (-2.0 * t).exp()
        });
        assert!(dev < 1e-10, "max deviation {dev}");
        assert!(
            response[10].1 < 0.0,
            "non-minimum phase undershoot expected"
        );
    }

    #[test]
    fn repeated_poles_use_the_exact_polynomial_form() {
        // 1 / (s + 1)^2: y(t) = 1 - (1 + t) e^-t.
        let pz = PolesZeros::new(1.0, vec![], vec![(-1.0, 0.0), (-1.0, 0.0)]);
        let response = pz.step_response_approx(0.01, 800);
        let dev = max_deviation(&response, |t| 1.0 - (1.0 + t) * (-t).exp());
        assert!(dev < 1e-10, "max deviation {dev}");
    }

    #[test]
    fn integrator_pole_at_the_origin_grows_linearly() {
        // 1 / (s (s + 1)): y(t) = t - 1 + e^-t.
        let pz = PolesZeros::new(1.0, vec![], vec![(0.0, 0.0), (-1.0, 0.0)]);
        let response = pz.step_response_approx(0.01, 500);
        let dev = max_deviation(&response, |t| t - 1.0 + (-t).exp());
        assert!(dev < 1e-10, "max deviation {dev}");
    }

    #[test]
    fn unstable_pole_grows_exponentially() {
        // 1 / (s - 0.5): y(t) = 2 (e^(0.5 t) - 1).
        let pz = PolesZeros::new(1.0, vec![], vec![(0.5, 0.0)]);
        let response = pz.step_response_approx(0.01, 300);
        let dev = max_deviation(&response, |t| 2.0 * ((0.5 * t).exp() - 1.0));
        assert!(dev < 1e-9, "max deviation {dev}");
    }

    #[test]
    fn biproper_model_starts_at_its_feedthrough() {
        // (s + 1) / (s + 2) = 1 - 1 / (s + 2): y(t) = 1/2 + 1/2 e^-2t, y(0) = 1.
        let pz = PolesZeros::new(1.0, vec![(-1.0, 0.0)], vec![(-2.0, 0.0)]);
        let response = pz.step_response_approx(0.01, 300);
        assert!((response[0].1 - 1.0).abs() < 1e-12);
        let dev = max_deviation(&response, |t| 0.5 + 0.5 * (-2.0 * t).exp());
        assert!(dev < 1e-10, "max deviation {dev}");
    }

    #[test]
    fn unpaired_complex_pole_returns_the_real_part() {
        // G(s) = 1 / (s - p), p = -1 + 2j has complex coefficients.  The step
        // response is (e^(pt) - 1) / p and its real part is
        // (1 - e^-t cos 2t + 2 e^-t sin 2t) / 5.
        let pz = PolesZeros::new(1.0, vec![], vec![(-1.0, 2.0)]);
        let response = pz.step_response_approx(0.01, 300);
        let dev = max_deviation(&response, |t| {
            (1.0 - (-t).exp() * (2.0 * t).cos() + 2.0 * (-t).exp() * (2.0 * t).sin()) / 5.0
        });
        assert!(dev < 1e-10, "max deviation {dev}");
    }

    #[test]
    fn no_poles_and_no_zeros_is_a_constant_gain() {
        let pz = PolesZeros::new(3.0, vec![], vec![]);
        let response = pz.step_response_approx(0.1, 5);
        assert_eq!(response.len(), 5);
        for &(_, y) in &response {
            assert_eq!(y, 3.0);
        }
    }

    #[test]
    fn invalid_models_and_grids_give_an_empty_response() {
        // Improper: a zero but no pole.
        let improper = PolesZeros::new(1.0, vec![(-1.0, 0.0)], vec![]);
        assert!(improper.step_response_approx(0.1, 10).is_empty());
        // Non-positive or non-finite dt.
        let pz = PolesZeros::new(1.0, vec![], vec![(-1.0, 0.0)]);
        for dt in [0.0, -0.1, f64::NAN, f64::INFINITY] {
            assert!(pz.step_response_approx(dt, 10).is_empty(), "dt = {dt}");
        }
        // Non-finite gain or pole.
        let bad_gain = PolesZeros::new(f64::NAN, vec![], vec![(-1.0, 0.0)]);
        assert!(bad_gain.step_response_approx(0.1, 10).is_empty());
        let bad_pole = PolesZeros::new(1.0, vec![], vec![(f64::INFINITY, 0.0)]);
        assert!(bad_pole.step_response_approx(0.1, 10).is_empty());
        // Zero steps, and a step count whose buffer cannot be allocated.
        assert!(pz.step_response_approx(0.1, 0).is_empty());
        assert!(pz.step_response_approx(0.1, usize::MAX).is_empty());
    }
}
