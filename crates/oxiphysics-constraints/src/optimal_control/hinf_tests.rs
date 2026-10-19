//! Tests of the H-infinity norm and the small gain test of
//! [`RobustControl`](super::RobustControl).
//!
//! Reference values come from dense frequency sweeps refined by golden-section
//! search (numpy / scipy) and from the closed-form peak
//! `1 / (2 ζ √(1 - ζ²))` of a second-order mode.

use super::hinf::StateSpace;
use super::types::RobustControl;
use crate::dense_linalg::real_eigenvalues;

/// `wn² / (s² + 2 ζ wn s + wn²)`.
fn second_order(zeta: f64, wn: f64) -> RobustControl {
    RobustControl::new(
        2,
        1,
        1,
        vec![0.0, 1.0, -wn * wn, -2.0 * zeta * wn],
        vec![0.0, wn * wn],
        vec![1.0, 0.0],
        vec![0.0],
    )
}

fn peak_gain(zeta: f64) -> f64 {
    1.0 / (2.0 * zeta * (1.0 - zeta * zeta).sqrt())
}

/// The bracket `[lower, upper]` must contain `exact`, with `lower` within the
/// advertised relative `3e-6` and `upper` never below the true value.
fn assert_brackets(rc: &RobustControl, exact: f64, what: &str) {
    let b = rc.h_inf_bounds();
    assert!(
        b.lower <= exact * (1.0 + 1e-12),
        "{what}: lower bound {} exceeds the true norm {exact}",
        b.lower
    );
    assert!(
        b.lower >= exact * (1.0 - 3e-6),
        "{what}: lower bound {} is not within 3e-6 of {exact}",
        b.lower
    );
    assert!(
        b.upper >= exact * (1.0 - 1e-12),
        "{what}: upper bound {} is below the true norm {exact}",
        b.upper
    );
    assert!(
        b.upper <= exact * (1.0 + 4e-6),
        "{what}: upper bound {} is not within 4e-6 of {exact}",
        b.upper
    );
    assert_eq!(rc.h_inf_norm(), b.lower);
    assert_eq!(rc.h_inf_upper_bound(), b.upper);
}

#[test]
fn lightly_damped_oscillator_gives_the_resonance_peak() {
    // 1 / (s^2 + 0.1 s + 1): the audit's reference system, peak 10.0125.
    let rc = second_order(0.05, 1.0);
    assert_brackets(&rc, 10.012523486435176, "zeta 0.05");
    assert!((rc.h_inf_norm() - 10.0125).abs() < 1e-3);
    // 1 / (s^2 + 0.02 s + 1): peak 50.0025.
    assert_brackets(&second_order(0.01, 1.0), 50.00250018751562, "zeta 0.01");
    for zeta in [0.3, 0.1, 1e-3, 1e-4, 1e-5] {
        assert_brackets(
            &second_order(zeta, 2.5),
            peak_gain(zeta),
            "second-order mode",
        );
    }
}

#[test]
fn peak_frequency_of_an_oscillator_is_reported() {
    let b = second_order(0.05, 1.0).h_inf_bounds();
    let expected = (1.0_f64 - 2.0 * 0.05 * 0.05).sqrt();
    assert!(
        (b.peak_frequency - expected).abs() < 1e-5 * expected,
        "peak at {} rad/s, expected {expected}",
        b.peak_frequency
    );
    // The scaled mode peaks at wn * sqrt(1 - 2 zeta^2).
    let b = second_order(0.05, 40.0).h_inf_bounds();
    assert!((b.peak_frequency - 40.0 * expected).abs() < 1e-5 * 40.0);
}

#[test]
fn well_damped_lag_peaks_at_dc() {
    // 1 / (s + 1): the supremum is the DC gain.
    let lag = RobustControl::new(1, 1, 1, vec![-1.0], vec![1.0], vec![1.0], vec![0.0]);
    assert_brackets(&lag, 1.0, "first-order lag");
    assert_eq!(lag.h_inf_bounds().peak_frequency, 0.0);
    // A pole very close to (but left of) the axis: 1 / (s + 1e-6) peaks at 1e6.
    let slow = RobustControl::new(1, 1, 1, vec![-1e-6], vec![1.0], vec![1.0], vec![0.0]);
    assert_brackets(&slow, 1e6, "slow lag");
    // Stiff cascade 1e6 / (s + 1e6) * 1 / (s + 1): also 1 at DC.
    let cascade = RobustControl::new(
        2,
        1,
        1,
        vec![-1.0, 0.0, 1e6, -1e6],
        vec![1.0, 0.0],
        vec![0.0, 1.0],
        vec![0.0],
    );
    assert_brackets(&cascade, 1.0, "stiff cascade");
}

#[test]
fn sharp_resonance_between_sweep_grid_points_is_not_missed() {
    // A broad mode at 0.4 rad/s plus a very sharp one (zeta 2e-4) at 3.7 rad/s:
    // the sharp peak is 1250 high but only ~1.5e-3 rad/s wide, so a 500-point
    // log sweep over eight decades (3.7 % spacing) steps right over it.
    let rc = RobustControl::new(
        4,
        1,
        1,
        vec![
            0.0, 1.0, 0.0, 0.0, //
            -0.16, -0.24, 0.0, 0.0, //
            0.0, 0.0, 0.0, 1.0, //
            0.0, 0.0, -13.69, -0.00148,
        ],
        vec![0.0, 0.16, 0.0, 13.69],
        vec![0.7, 0.0, 0.5, 0.0],
        vec![0.0],
    );
    assert_brackets(&rc, 1250.0005643703023, "two-mode plant");
    let peak = rc.h_inf_bounds().peak_frequency;
    assert!((peak - 3.7).abs() < 1e-3, "peak at {peak}");
}

#[test]
fn stiff_plant_with_a_resonance_is_resolved() {
    // Slow lag (-0.01) in parallel with a stiff resonance at 1000 rad/s
    // (zeta 1e-3): entries of A span eight orders of magnitude.
    let rc = RobustControl::new(
        3,
        1,
        1,
        vec![
            -0.01, 0.0, 0.0, //
            0.0, 0.0, 1.0, //
            0.0, -1.0e6, -2.0,
        ],
        vec![1.0, 0.0, 1.0e6],
        vec![0.5, 1.0, 0.0],
        vec![0.0],
    );
    assert_brackets(&rc, 500.00075000044245, "stiff resonance");
}

#[test]
fn mimo_plant_with_feedthrough_matches_the_reference_norm() {
    let a = vec![
        -0.5, 2.0, 0.0, 0.0, //
        -2.0, -0.5, 0.0, 0.0, //
        0.0, 0.0, -0.2, 5.0, //
        0.0, 0.0, -5.0, -0.2,
    ];
    let b = vec![1.0, 0.0, 0.0, 1.0, 1.0, 0.5, 0.0, -1.0];
    let c = vec![1.0, 0.0, 0.5, 0.0, 0.0, 1.0, 0.0, 2.0];
    let d = vec![0.2, -0.1, 0.0, 0.3];
    let rc = RobustControl::new(4, 2, 2, a.clone(), b.clone(), c.clone(), d.clone());
    assert_brackets(&rc, 7.540196816557535, "mimo");
    let peak = rc.h_inf_bounds().peak_frequency;
    assert!((peak - 5.0003).abs() < 2e-3, "peak at {peak}");

    // Rotating the outputs by an orthogonal matrix leaves every singular value
    // of G(jw), and therefore the norm, unchanged.
    let (cs, sn) = (0.7_f64.cos(), 0.7_f64.sin());
    let rotate = |m: &[f64]| -> Vec<f64> {
        let mut out = vec![0.0; m.len()];
        let cols = m.len() / 2;
        for j in 0..cols {
            out[j] = cs * m[j] - sn * m[cols + j];
            out[cols + j] = sn * m[j] + cs * m[cols + j];
        }
        out
    };
    let rotated = RobustControl::new(4, 2, 2, a, b, rotate(&c), rotate(&d));
    assert_brackets(&rotated, 7.540196816557535, "rotated mimo");
}

#[test]
fn static_gain_norm_is_the_largest_singular_value_of_d() {
    let rc = RobustControl::new(0, 2, 2, vec![], vec![], vec![], vec![3.0, 0.0, 0.0, -2.0]);
    let b = rc.h_inf_bounds();
    assert!((b.lower - 3.0).abs() < 1e-12 && (b.upper - 3.0).abs() < 1e-12);
    assert!(b.peak_frequency.is_infinite());
    assert!(rc.small_gain_stable(0.3));
    assert!(!rc.small_gain_stable(0.4));
}

#[test]
fn norm_reached_only_at_infinite_frequency_is_the_feedthrough() {
    // G = -3 + 1 / (s + 1): |G(jw)| rises towards 3 as w -> infinity.
    let rc = RobustControl::new(1, 1, 1, vec![-1.0], vec![1.0], vec![1.0], vec![-3.0]);
    let b = rc.h_inf_bounds();
    assert!((b.lower - 3.0).abs() < 1e-12, "lower {}", b.lower);
    assert!(
        b.upper >= 3.0 && b.upper <= 3.0 * (1.0 + 4e-6),
        "upper {}",
        b.upper
    );
    assert!(b.peak_frequency.is_infinite());
    // G = 3 + 1 / (s + 1) peaks at DC with 4.
    let dc = RobustControl::new(1, 1, 1, vec![-1.0], vec![1.0], vec![1.0], vec![3.0]);
    assert_brackets(&dc, 4.0, "dc dominated");
}

#[test]
fn unstable_marginal_and_malformed_plants_have_infinite_norm() {
    let cases: Vec<(&str, RobustControl)> = vec![
        (
            "unstable pole",
            RobustControl::new(1, 1, 1, vec![0.5], vec![1.0], vec![1.0], vec![0.0]),
        ),
        (
            "double integrator",
            RobustControl::new(
                2,
                1,
                1,
                vec![0.0, 1.0, 0.0, 0.0],
                vec![0.0, 1.0],
                vec![1.0, 0.0],
                vec![0.0],
            ),
        ),
        (
            "undamped oscillator",
            RobustControl::new(
                2,
                1,
                1,
                vec![0.0, 1.0, -1.0, 0.0],
                vec![0.0, 1.0],
                vec![1.0, 0.0],
                vec![0.0],
            ),
        ),
        (
            "zero state matrix",
            RobustControl::new(1, 1, 1, vec![0.0], vec![1.0], vec![1.0], vec![0.0]),
        ),
        (
            "wrong B size",
            RobustControl::new(
                2,
                1,
                1,
                vec![-1.0, 0.0, 0.0, -1.0],
                vec![1.0],
                vec![1.0, 0.0],
                vec![0.0],
            ),
        ),
        (
            "non-finite A",
            RobustControl::new(1, 1, 1, vec![f64::NAN], vec![1.0], vec![1.0], vec![0.0]),
        ),
    ];
    for (what, rc) in cases {
        let b = rc.h_inf_bounds();
        assert!(
            b.lower.is_infinite() && b.upper.is_infinite(),
            "{what}: {b:?}"
        );
        assert!(rc.h_inf_norm().is_infinite(), "{what}");
        assert!(rc.h_inf_upper_bound().is_infinite(), "{what}");
        assert!(
            !rc.small_gain_stable(1e-12),
            "{what} must never be robustly stable"
        );
        assert!(
            !rc.small_gain_stable(0.0),
            "{what} must never be robustly stable"
        );
    }
}

#[test]
fn plant_without_input_or_gain_has_zero_norm() {
    let no_input = RobustControl::new(1, 0, 1, vec![-1.0], vec![], vec![1.0], vec![]);
    let b = no_input.h_inf_bounds();
    assert_eq!((b.lower, b.upper), (0.0, 0.0));
    let silent = RobustControl::new(1, 1, 1, vec![-1.0], vec![0.0], vec![1.0], vec![0.0]);
    let b = silent.h_inf_bounds();
    assert_eq!((b.lower, b.upper), (0.0, 0.0));
    assert!(silent.small_gain_stable(1e9));
}

#[test]
fn norm_scales_linearly_with_the_input_matrix() {
    let base = second_order(0.05, 1.0);
    let mut doubled = second_order(0.05, 1.0);
    for v in doubled.b.iter_mut() {
        *v *= 2.0;
    }
    let ratio = doubled.h_inf_norm() / base.h_inf_norm();
    assert!((ratio - 2.0).abs() < 1e-5, "ratio {ratio}");
}

#[test]
fn small_gain_test_uses_the_peak_and_not_the_dc_gain() {
    // 1 / (s^2 + 0.02 s + 1) has DC gain 1 but peak gain 50.0025: the previous
    // DC-gain proxy declared gamma = 0.5 robustly stable although the loop
    // gain at resonance is 25.
    let rc = second_order(0.01, 1.0);
    assert!(!rc.small_gain_stable(0.5));
    assert!(!rc.small_gain_stable(0.1));
    assert!(rc.small_gain_stable(0.01));
    // The threshold sits at 1 / ||G||_inf.
    let norm = peak_gain(0.01);
    assert!(rc.small_gain_stable(0.999 / norm));
    assert!(!rc.small_gain_stable(1.001 / norm));
    // Negative uncertainty bounds are interpreted by magnitude.
    assert!(rc.small_gain_stable(-0.01));
    assert!(!rc.small_gain_stable(f64::NAN));
}

#[test]
fn sigma_max_matches_closed_forms() {
    let lag = StateSpace {
        n: 1,
        m: 1,
        p: 1,
        a: &[-1.0],
        b: &[1.0],
        c: &[1.0],
        d: &[0.0],
    };
    assert!((lag.sigma_max(0.0) - 1.0).abs() < 1e-14);
    assert!((lag.sigma_max(1.0) - std::f64::consts::FRAC_1_SQRT_2).abs() < 1e-14);
    assert!((lag.sigma_max(1e6) - 1e-6).abs() < 1e-12);
    // diag(1/(s+1), 2/(s+2)): sigma_max is the larger of the two magnitudes.
    let diag = StateSpace {
        n: 2,
        m: 2,
        p: 2,
        a: &[-1.0, 0.0, 0.0, -2.0],
        b: &[1.0, 0.0, 0.0, 2.0],
        c: &[1.0, 0.0, 0.0, 1.0],
        d: &[0.0; 4],
    };
    assert!((diag.sigma_max(0.0) - 1.0).abs() < 1e-14);
    assert!((diag.sigma_max(1.0) - 2.0 / 5.0_f64.sqrt()).abs() < 1e-14);
    // A pole on the imaginary axis makes jwI - A singular.
    let oscillator = StateSpace {
        n: 2,
        m: 1,
        p: 1,
        a: &[0.0, 1.0, -1.0, 0.0],
        b: &[0.0, 1.0],
        c: &[1.0, 0.0],
        d: &[0.0],
    };
    assert!(oscillator.sigma_max(1.0).is_infinite());
}

#[test]
fn hamiltonian_has_imaginary_eigenvalues_exactly_where_the_gain_crosses_the_level() {
    let rc = second_order(0.05, 1.0);
    let sys = StateSpace {
        n: rc.n,
        m: rc.m,
        p: rc.p,
        a: &rc.a,
        b: &rc.b,
        c: &rc.c,
        d: &rc.d,
    };
    // |G(jw)| = 5 at w = 0.90649 and w = 1.08087 (levels below the 10.0125 peak).
    let crossings = sys.crossing_frequencies(5.0).expect("eigenvalues converge");
    let expected = [0.9064892528256281, 1.0808687406487592];
    assert_eq!(crossings.len(), 2, "crossings {crossings:?}");
    for (got, want) in crossings.iter().zip(expected) {
        assert!((got - want).abs() < 1e-8, "crossing {got} vs {want}");
        assert!((sys.sigma_max(*got) - 5.0).abs() < 1e-8);
    }
    // Above the peak the Hamiltonian has no eigenvalue on the axis ...
    assert!(
        sys.crossing_frequencies(10.5)
            .expect("converges")
            .is_empty()
    );
    // ... and below the DC gain the only crossing besides the resonance pair
    // does not exist either: levels in (0, 1) are crossed once at high
    // frequency and once more below the peak.
    let low = sys.crossing_frequencies(0.5).expect("converges");
    assert_eq!(low.len(), 1, "crossings {low:?}");
    assert!((sys.sigma_max(low[0]) - 0.5).abs() < 1e-8);
    // The Hamiltonian spectrum is symmetric about the imaginary axis.
    let h = sys.hamiltonian(5.0).expect("gamma exceeds sigma_max(D)");
    let eigenvalues = real_eigenvalues(&h, 4).expect("converges");
    for (re, im) in &eigenvalues {
        assert!(
            eigenvalues
                .iter()
                .any(|(r2, i2)| (r2 + re).abs() < 1e-8 && (i2 - im).abs() < 1e-8),
            "no mirror image of {re} + {im}i"
        );
    }
    // R = gamma^2 I - D^T D must be positive definite.
    let with_feedthrough = StateSpace {
        n: 1,
        m: 1,
        p: 1,
        a: &[-1.0],
        b: &[1.0],
        c: &[1.0],
        d: &[2.0],
    };
    assert!(with_feedthrough.hamiltonian(2.0).is_none());
    assert!(with_feedthrough.hamiltonian(2.1).is_some());
}

#[test]
fn mimo_hamiltonian_crossings_lie_on_the_level_set() {
    let a = [
        -0.5, 2.0, 0.0, 0.0, //
        -2.0, -0.5, 0.0, 0.0, //
        0.0, 0.0, -0.2, 5.0, //
        0.0, 0.0, -5.0, -0.2,
    ];
    let b = [1.0, 0.0, 0.0, 1.0, 1.0, 0.5, 0.0, -1.0];
    let c = [1.0, 0.0, 0.5, 0.0, 0.0, 1.0, 0.0, 2.0];
    let d = [0.2, -0.1, 0.0, 0.3];
    let sys = StateSpace {
        n: 4,
        m: 2,
        p: 2,
        a: &a,
        b: &b,
        c: &c,
        d: &d,
    };
    let norm = 7.540196816557535;
    let level = 0.5 * norm;
    let crossings = sys.crossing_frequencies(level).expect("converges");
    assert!(!crossings.is_empty());
    for omega in &crossings {
        let sigma = sys.sigma_max(*omega);
        assert!(
            (sigma - level).abs() < 1e-7 * level,
            "sigma_max({omega}) = {sigma} should equal {level}"
        );
    }
    assert!(
        sys.crossing_frequencies(1.01 * norm)
            .expect("converges")
            .is_empty()
    );
}

// ── Randomized cross-check against dense frequency sweeps ───────────────

/// Deterministic pseudo-random numbers (64-bit LCG, Box–Muller for normals).
struct Lcg(u64);

impl Lcg {
    fn uniform(&mut self) -> f64 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        ((self.0 >> 11) as f64) / ((1_u64 << 53) as f64)
    }

    fn normal(&mut self) -> f64 {
        let (u1, u2) = (self.uniform().max(1e-300), self.uniform());
        (-2.0 * u1.ln()).sqrt() * (2.0 * std::f64::consts::PI * u2).cos()
    }

    fn below(&mut self, bound: usize) -> usize {
        ((self.uniform() * bound as f64) as usize).min(bound - 1)
    }
}

#[test]
fn random_stable_systems_are_bracketed_by_dense_sweeps() {
    // Random MIMO systems (up to 7 states and 3 x 3 channels, half of them with
    // a feedthrough), shifted to be stable.  Every value of sigma_max(G(jw)) on
    // a fine log grid is a lower bound of the norm, so the certified upper
    // bound can never be below the largest of them; the reported norm must be
    // within the advertised tolerance of it or above it, the bracket must be
    // tight, and the reported peak frequency must really attain the norm.
    let mut rng = Lcg(0x1234_5678_9abc_def0);
    let (log_lo, log_hi) = (1e-3_f64.ln(), 1e3_f64.ln());
    let steps = 3000;
    for trial in 0..60 {
        let n = 1 + rng.below(7);
        let (m, p) = (1 + rng.below(3), 1 + rng.below(3));
        let mut a: Vec<f64> = (0..n * n).map(|_| rng.normal()).collect();
        let abscissa = real_eigenvalues(&a, n)
            .expect("eigenvalues converge")
            .iter()
            .fold(f64::NEG_INFINITY, |worst, e| worst.max(e.0));
        let margin = 0.02 + 0.6 * rng.uniform();
        for i in 0..n {
            a[i * n + i] -= abscissa + margin;
        }
        let b: Vec<f64> = (0..n * m).map(|_| rng.normal()).collect();
        let c: Vec<f64> = (0..p * n).map(|_| rng.normal()).collect();
        let d: Vec<f64> = if rng.uniform() < 0.5 {
            vec![0.0; p * m]
        } else {
            (0..p * m).map(|_| 0.5 * rng.normal()).collect()
        };
        let sys = StateSpace {
            n,
            m,
            p,
            a: &a,
            b: &b,
            c: &c,
            d: &d,
        };
        let mut sweep = sys.sigma_max(0.0);
        for k in 0..steps {
            let w = log_lo + (log_hi - log_lo) * k as f64 / (steps - 1) as f64;
            sweep = sweep.max(sys.sigma_max(w.exp()));
        }
        let bounds =
            RobustControl::new(n, m, p, a.clone(), b.clone(), c.clone(), d.clone()).h_inf_bounds();
        assert!(
            bounds.lower.is_finite() && bounds.upper.is_finite(),
            "trial {trial}: {bounds:?}"
        );
        assert!(
            bounds.upper >= sweep * (1.0 - 1e-12),
            "trial {trial} (n={n}, m={m}, p={p}): upper bound {} below the sweep {sweep}",
            bounds.upper
        );
        assert!(
            bounds.lower >= sweep * (1.0 - 4e-6),
            "trial {trial}: norm {} below the sweep {sweep}",
            bounds.lower
        );
        assert!(
            bounds.upper <= bounds.lower * (1.0 + 4e-6),
            "trial {trial}: bracket {bounds:?} is not tight"
        );
        if bounds.peak_frequency.is_finite() {
            let attained = sys.sigma_max(bounds.peak_frequency);
            assert!(
                (attained - bounds.lower).abs() <= 1e-9 * bounds.lower,
                "trial {trial}: sigma_max at the reported peak is {attained}, norm {}",
                bounds.lower
            );
        }
    }
}
