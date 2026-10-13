//! Gain and phase margins of a [`BodeAnalysis`], evaluated on the *continuous*
//! (unwrapped) phase of the frequency response.
//!
//! `atan2` returns phases in `(-180, 180]`, so a lag that passes through
//! `-180` re-enters the interval at `+180`.  Testing the wrapped phase for a
//! `-180` crossing therefore misses the very crossover the gain margin is
//! defined at.  The sweep is first unwrapped, crossings of `-180 - 360 k`
//! degrees are located on the continuous phase, and each crossing is refined
//! by bisection on the exact frequency response.

use super::{BodeAnalysis, BodePoint};

/// Folds an angle difference in degrees into `(-180, 180]`.
fn wrap_to_180(angle_deg: f64) -> f64 {
    let wrapped = angle_deg - 360.0 * (angle_deg / 360.0).round();
    if wrapped <= -180.0 {
        wrapped + 360.0
    } else {
        wrapped
    }
}

/// Continuous phase of a sweep: every sample is shifted by a multiple of 360
/// degrees so that consecutive samples differ by at most 180 degrees.
///
/// A sample whose phase is not finite (a sweep frequency that hits a pole on
/// the imaginary axis exactly) stays NaN and is skipped: the next finite
/// sample is unwrapped against the last finite one, so a single singular
/// sample does not poison the rest of the sweep.
fn unwrap_phase_deg(points: &[BodePoint]) -> Vec<f64> {
    let mut unwrapped = Vec::with_capacity(points.len());
    let mut offset = 0.0_f64;
    let mut previous_finite: Option<f64> = None;
    for point in points {
        if point.phase_deg.is_finite() {
            if let Some(previous) = previous_finite {
                offset -= 360.0 * ((point.phase_deg - previous) / 360.0).round();
            }
            previous_finite = Some(point.phase_deg);
        }
        unwrapped.push(point.phase_deg + offset);
    }
    unwrapped
}

/// Phase level `-180 - 360 k` (degrees) crossed between two consecutive
/// unwrapped samples, if any.
fn phase_crossing_level(phase_before: f64, phase_after: f64) -> Option<f64> {
    if !(phase_before.is_finite() && phase_after.is_finite()) {
        return None;
    }
    let distance_before = phase_before + 180.0;
    let distance_after = phase_after + 180.0;
    let lowest = distance_before.min(distance_after);
    let highest = distance_before.max(distance_after);
    let level = 360.0 * (lowest / 360.0).ceil();
    (level <= highest).then_some(level - 180.0)
}

/// Linear interpolation fraction of `target` between `before` and `after`.
fn crossing_fraction(before: f64, after: f64, target: f64) -> f64 {
    let span = after - before;
    if span == 0.0 {
        0.0
    } else {
        ((target - before) / span).clamp(0.0, 1.0)
    }
}

/// Bisection (in `log omega`) for a sign change of `f` on `[lo, hi]`.
///
/// Returns `None` if `f` does not change sign on the interval or yields NaN.
fn bisect_log_omega(lo: f64, hi: f64, f: impl Fn(f64) -> f64) -> Option<f64> {
    let (mut lo, mut hi) = (lo, hi);
    if !(lo > 0.0 && hi > lo) {
        return None;
    }
    let mut f_lo = f(lo);
    let f_hi = f(hi);
    if f_lo == 0.0 {
        return Some(lo);
    }
    if f_hi == 0.0 {
        return Some(hi);
    }
    if f_lo.is_nan() || f_hi.is_nan() || (f_lo < 0.0) == (f_hi < 0.0) {
        return None;
    }
    for _ in 0..200 {
        let mid = (lo * hi).sqrt();
        let f_mid = f(mid);
        if f_mid.is_nan() {
            return None;
        }
        if f_mid == 0.0 {
            return Some(mid);
        }
        if (f_mid < 0.0) == (f_lo < 0.0) {
            lo = mid;
            f_lo = f_mid;
        } else {
            hi = mid;
        }
        if hi / lo - 1.0 < 1e-15 {
            break;
        }
    }
    Some((lo * hi).sqrt())
}

impl BodeAnalysis {
    /// Magnitude in dB at `omega` (`-300` dB for an exactly vanishing gain).
    fn gain_db_at(&self, omega: f64) -> f64 {
        let (gain, _) = self.eval_at(omega);
        if gain > 1e-300 {
            20.0 * gain.log10()
        } else {
            -300.0
        }
    }

    /// Phase at `omega` continued from a reference sample on the unwrapped
    /// branch: `unwrapped_ref + wrap(phase - wrapped_ref)`.
    fn continuous_phase_at(&self, omega: f64, wrapped_ref: f64, unwrapped_ref: f64) -> f64 {
        let (_, phase) = self.eval_at(omega);
        unwrapped_ref + wrap_to_180(phase - wrapped_ref)
    }

    /// Gain margin in dB at the first phase crossover, `+inf` if none exists.
    fn gain_margin_db(&self, points: &[BodePoint], unwrapped: &[f64]) -> f64 {
        for i in 1..points.len() {
            let (before, after) = (&points[i - 1], &points[i]);
            let Some(target) = phase_crossing_level(unwrapped[i - 1], unwrapped[i]) else {
                continue;
            };
            let refined = bisect_log_omega(before.omega, after.omega, |omega| {
                self.continuous_phase_at(omega, before.phase_deg, unwrapped[i - 1]) - target
            });
            let gain_db = match refined {
                Some(omega) => self.gain_db_at(omega),
                None => {
                    let fraction = crossing_fraction(unwrapped[i - 1], unwrapped[i], target);
                    before.gain_db + fraction * (after.gain_db - before.gain_db)
                }
            };
            return -gain_db;
        }
        f64::INFINITY
    }

    /// Phase margin in degrees at the first gain crossover, `+inf` if none exists.
    fn phase_margin_deg(&self, points: &[BodePoint], unwrapped: &[f64]) -> f64 {
        for i in 1..points.len() {
            let (before, after) = (&points[i - 1], &points[i]);
            // A sample with a non-finite phase sits exactly on a pole: its gain
            // is a sentinel, not a measurement, so it cannot bracket a crossover.
            let usable = before.phase_deg.is_finite() && after.phase_deg.is_finite();
            let crosses_zero_db = usable && before.gain_db * after.gain_db <= 0.0;
            if !crosses_zero_db {
                continue;
            }
            let refined =
                bisect_log_omega(before.omega, after.omega, |omega| self.gain_db_at(omega));
            let phase_deg = match refined {
                Some(omega) => self.continuous_phase_at(omega, before.phase_deg, unwrapped[i - 1]),
                None => {
                    let fraction = crossing_fraction(before.gain_db, after.gain_db, 0.0);
                    unwrapped[i - 1] + fraction * (unwrapped[i] - unwrapped[i - 1])
                }
            };
            return wrap_to_180(180.0 + phase_deg);
        }
        f64::INFINITY
    }

    /// Find gain margin (dB) and phase margin (degrees).
    ///
    /// The frequency response is sampled at `n_points` logarithmically spaced
    /// frequencies in `[omega_min, omega_max]`, its phase is unwrapped along the
    /// sweep, and the first crossings are located and refined on the exact
    /// response:
    ///
    /// - Gain margin: `-|G(jw)|` in dB at the first phase crossover, i.e. the
    ///   first frequency where the continuous phase equals `-180 - 360 k`
    ///   degrees for an integer `k`.
    /// - Phase margin: `180 + phase` in degrees at the first gain crossover
    ///   (`|G(jw)| = 1`), folded into `(-180, 180]` so that the result does not
    ///   depend on which `360 k` branch the sweep started on.
    ///
    /// Returns `(gain_margin_db, phase_margin_deg)`.  A margin whose crossover
    /// does not occur inside the sweep is reported as `f64::INFINITY` (for
    /// example the gain margin of `1/(s+1)` or the phase margin of
    /// `1/(s+1)^3`, whose gain stays below 0 dB).  Fewer than two sweep
    /// points cannot reveal any crossing, and a sweep whose frequencies are not
    /// positive and finite (for example `omega_min <= 0`) is meaningless; both
    /// give `(NaN, NaN)`.
    pub fn margins(&self, omega_min: f64, omega_max: f64, n_points: usize) -> (f64, f64) {
        let points = self.compute_log(omega_min, omega_max, n_points);
        if points.len() < 2
            || points
                .iter()
                .any(|p| !(p.omega.is_finite() && p.omega > 0.0))
        {
            return (f64::NAN, f64::NAN);
        }
        let unwrapped = unwrap_phase_deg(&points);
        (
            self.gain_margin_db(&points, &unwrapped),
            self.phase_margin_deg(&points, &unwrapped),
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn cubic_lag(gain: f64) -> BodeAnalysis {
        // gain / (s + 1)^3
        BodeAnalysis::new(vec![gain], vec![1.0, 3.0, 3.0, 1.0])
    }

    #[test]
    fn wrap_to_180_folds_into_half_open_interval() {
        assert_eq!(wrap_to_180(0.0), 0.0);
        assert_eq!(wrap_to_180(190.0), -170.0);
        assert_eq!(wrap_to_180(-190.0), 170.0);
        assert_eq!(wrap_to_180(180.0), 180.0);
        assert_eq!(wrap_to_180(-180.0), 180.0);
        assert!((wrap_to_180(725.0) - 5.0).abs() < 1e-12);
    }

    #[test]
    fn unwrap_removes_the_wrap_discontinuity() {
        let make = |phase_deg: f64| BodePoint {
            omega: 1.0,
            gain_db: 0.0,
            phase_deg,
        };
        // A lag sweeping through -180: the wrapped phase jumps from -178 to +178.
        let pts = [
            make(-170.0),
            make(-178.0),
            make(178.0),
            make(170.0),
            make(160.0),
        ];
        let u = unwrap_phase_deg(&pts);
        let expected = [-170.0, -178.0, -182.0, -190.0, -200.0];
        assert_eq!(u.len(), expected.len());
        for (got, want) in u.iter().zip(expected.iter()) {
            assert!((got - want).abs() < 1e-12, "u = {u:?}");
        }
    }

    #[test]
    fn unwrap_skips_non_finite_samples_without_poisoning_the_sweep() {
        let make = |phase_deg: f64| BodePoint {
            omega: 1.0,
            gain_db: 0.0,
            phase_deg,
        };
        let pts = [
            make(-170.0),
            make(-178.0),
            make(f64::NAN),
            make(178.0),
            make(170.0),
        ];
        let u = unwrap_phase_deg(&pts);
        assert!(u[2].is_nan());
        for (got, want) in [
            (u[0], -170.0),
            (u[1], -178.0),
            (u[3], -182.0),
            (u[4], -190.0),
        ] {
            assert!((got - want).abs() < 1e-12, "u = {u:?}");
        }
    }

    #[test]
    fn a_sweep_sample_on_an_imaginary_axis_pole_does_not_hide_later_crossovers() {
        // 1e5 / ((s^2 + 1)(s + 10)^5) = 1 / ((s^2 + 1)(s/10 + 1)^5): the 41-point
        // sweep over [0.1, 1000] contains w = 1 exactly, where the denominator
        // vanishes and the phase sample is NaN.  Above the resonance the phase
        // is -180 - 5 atan(w / 10) degrees, which reaches -540 degrees at
        // w = 10 tan(72 deg).  The gain margin there is
        // 20 log10(|1 - w^2| (1 + (w / 10)^2)^(5/2)).
        let plant = BodeAnalysis::new(
            vec![1e5],
            vec![
                1.0, 50.0, 1001.0, 10050.0, 51000.0, 110000.0, 50000.0, 100000.0,
            ],
        );
        let sweep = plant.compute_log(0.1, 1000.0, 41);
        assert_eq!(sweep[10].omega, 1.0);
        assert!(sweep[10].phase_deg.is_nan(), "{:?}", sweep[10]);
        let w = 10.0 * 72.0_f64.to_radians().tan();
        let expected = 20.0 * ((w * w - 1.0) * (1.0 + (w / 10.0).powi(2)).powf(2.5)).log10();
        let (gm, _) = plant.margins(0.1, 1000.0, 41);
        assert!(
            (gm - expected).abs() < 1e-6,
            "gm = {gm}, expected {expected}"
        );
    }

    #[test]
    fn crossing_level_detects_every_odd_half_turn() {
        assert_eq!(phase_crossing_level(-170.0, -190.0), Some(-180.0));
        assert_eq!(phase_crossing_level(-190.0, -170.0), Some(-180.0));
        assert_eq!(phase_crossing_level(-500.0, -540.0), Some(-540.0));
        assert_eq!(phase_crossing_level(170.0, 190.0), Some(180.0));
        assert_eq!(phase_crossing_level(-100.0, -150.0), None);
        assert_eq!(phase_crossing_level(-181.0, -250.0), None);
        // A non-finite sample never brackets a crossing, even next to a level.
        assert_eq!(phase_crossing_level(-180.0, f64::NAN), None);
        assert_eq!(phase_crossing_level(f64::NAN, 180.0), None);
    }

    #[test]
    fn cubic_lag_gain_margin_is_found_across_the_phase_wrap() {
        // 1/(s+1)^3: phase crossover at w = sqrt(3) where |G| = 1/8, so the
        // gain margin is 20 log10(8) = 18.0618 dB.  The wrapped-phase test
        // used to return NaN for this system.
        let expected = 20.0 * 8.0_f64.log10();
        let (gm, pm) = cubic_lag(1.0).margins(0.01, 100.0, 2000);
        assert!(
            (gm - expected).abs() < 1e-6,
            "gm = {gm}, expected {expected}"
        );
        assert!(
            pm.is_infinite() && pm > 0.0,
            "gain never reaches 0 dB, pm = {pm}"
        );
    }

    #[test]
    fn coarse_sweeps_are_refined_on_the_exact_response() {
        let expected = 20.0 * 8.0_f64.log10();
        let (gm, _) = cubic_lag(1.0).margins(0.01, 100.0, 25);
        assert!(
            (gm - expected).abs() < 1e-6,
            "gm = {gm}, expected {expected}"
        );
    }

    #[test]
    fn stable_loop_has_positive_margins_matching_closed_forms() {
        // K = 4: GM = 20 log10(8/4) = 6.0206 dB; the gain crossover solves
        // (1 + w^2)^(3/2) = 4 and PM = 180 - 3 atan(w) degrees.
        let (gm, pm) = cubic_lag(4.0).margins(0.01, 100.0, 400);
        let w_gc = (4.0_f64.powf(2.0 / 3.0) - 1.0).sqrt();
        let expected_pm = 180.0 - 3.0 * w_gc.atan().to_degrees();
        assert!((gm - 20.0 * 2.0_f64.log10()).abs() < 1e-6, "gm = {gm}");
        assert!(
            (pm - expected_pm).abs() < 1e-6,
            "pm = {pm}, expected {expected_pm}"
        );
        assert!(pm > 0.0 && gm > 0.0);
    }

    #[test]
    fn unstable_loop_has_negative_margins() {
        // K = 16: GM = 20 log10(8/16) = -6.0206 dB, PM = 180 - 3 atan(w_gc) < 0.
        let (gm, pm) = cubic_lag(16.0).margins(0.01, 100.0, 400);
        let w_gc = (16.0_f64.powf(2.0 / 3.0) - 1.0).sqrt();
        let expected_pm = 180.0 - 3.0 * w_gc.atan().to_degrees();
        assert!((gm + 20.0 * 2.0_f64.log10()).abs() < 1e-6, "gm = {gm}");
        assert!(
            (pm - expected_pm).abs() < 1e-6,
            "pm = {pm}, expected {expected_pm}"
        );
        assert!(gm < 0.0 && pm < 0.0);
    }

    #[test]
    fn type_two_plant_phase_margin_is_independent_of_the_starting_branch() {
        // 1/(s^2 (s+1)): the phase starts just below -180 degrees, so the wrapped
        // phase starts near +180.  The gain crossover solves x^3 + x^2 = 1 with
        // x = w^2 and PM = -atan(w) < 0 (the loop 1/(s^2(s+1)) is unstable).
        let plant = BodeAnalysis::new(vec![1.0], vec![1.0, 1.0, 0.0, 0.0]);
        let mut x = 0.75_f64;
        for _ in 0..60 {
            x -= (x * x * x + x * x - 1.0) / (3.0 * x * x + 2.0 * x);
        }
        let expected_pm = -x.sqrt().atan().to_degrees();
        let (_, pm) = plant.margins(0.01, 100.0, 500);
        assert!(
            (pm - expected_pm).abs() < 1e-6,
            "pm = {pm}, expected {expected_pm}"
        );
        assert!(pm < 0.0);
    }

    #[test]
    fn first_order_loop_has_unbounded_gain_margin() {
        // 10/(s+1): the phase only approaches -90 degrees; PM = 180 - atan(sqrt(99)).
        let plant = BodeAnalysis::new(vec![10.0], vec![1.0, 1.0]);
        let (gm, pm) = plant.margins(0.01, 1000.0, 300);
        assert!(gm.is_infinite() && gm > 0.0, "gm = {gm}");
        let expected_pm = 180.0 - 99.0_f64.sqrt().atan().to_degrees();
        assert!(
            (pm - expected_pm).abs() < 1e-6,
            "pm = {pm}, expected {expected_pm}"
        );
    }

    #[test]
    fn degenerate_sweeps_return_nan() {
        let plant = cubic_lag(1.0);
        for n in [0usize, 1] {
            let (gm, pm) = plant.margins(0.1, 10.0, n);
            assert!(gm.is_nan() && pm.is_nan(), "n = {n}: ({gm}, {pm})");
        }
    }

    #[test]
    fn invalid_sweep_ranges_return_nan() {
        let plant = cubic_lag(1.0);
        for (lo, hi) in [
            (0.0, 10.0),
            (-1.0, 10.0),
            (0.1, f64::INFINITY),
            (f64::NAN, 10.0),
        ] {
            let (gm, pm) = plant.margins(lo, hi, 20);
            assert!(gm.is_nan() && pm.is_nan(), "({lo}, {hi}): ({gm}, {pm})");
        }
    }

    #[test]
    fn flat_minus_180_phase_does_not_produce_nan() {
        // 1/s^2 has phase exactly -180 everywhere; the crossing fraction must not be 0/0.
        let plant = BodeAnalysis::new(vec![1.0], vec![1.0, 0.0, 0.0]);
        let (gm, _) = plant.margins(0.1, 10.0, 20);
        assert!(!gm.is_nan(), "gm = {gm}");
    }
}
