//! Regression tests for the MMA sub-problem and the volume-constrained MMA step.
//!
//! The convex separable approximation `p/(U − x) + q/(x − L)` has its minimiser at
//! `(L + r U)/(1 + r)`, `r = √(q/p)`.  Earlier code swapped `L` and `U` (returning
//! the wrong side of the interval) and produced `NaN` for negative sensitivities
//! (`p = 0 ⇒ r = ∞`).

use super::*;

/// The MMA approximation of one variable.
fn approximation(p: f64, q: f64, l: f64, u: f64, x: f64) -> f64 {
    p / (u - x) + q / (x - l)
}

/// Brute-force minimiser of the approximation on `[a, b]`.
fn grid_minimiser(p: f64, q: f64, l: f64, u: f64, a: f64, b: f64) -> f64 {
    let samples = 200_000;
    let mut best_x = a;
    let mut best = f64::INFINITY;
    for k in 0..=samples {
        let x = a + (b - a) * k as f64 / samples as f64;
        let value = approximation(p, q, l, u, x);
        if value < best {
            best = value;
            best_x = x;
        }
    }
    best_x
}

fn fresh_asymptotes(n: usize, x0: f64) -> Vec<MmaAsymptote> {
    (0..n).map(|_| MmaAsymptote::new(x0, 0.5)).collect()
}

// ── single-variable sub-problem ─────────────────────────────────────────────

#[test]
fn single_variable_minimiser_matches_brute_force() {
    // (p, q, L, U, alpha, beta)
    let cases = [
        (1.0, 4.0, 0.0, 3.0, 0.5, 2.9),   // interior, r = 2  ⇒ x* = 2
        (2.0, 2.0, 0.0, 1.0, 0.1, 0.9),   // r = 1  ⇒ midpoint of the asymptotes
        (4.0, 1.0, 0.0, 3.0, 0.2, 2.8),   // r = ½  ⇒ x* = 1
        (1.0, 100.0, 0.0, 1.0, 0.1, 0.6), // clamped at beta
        (100.0, 1.0, 0.0, 1.0, 0.3, 0.9), // clamped at alpha
        (0.3, 0.7, -1.0, 2.0, -0.5, 1.5), // asymmetric asymptotes
    ];
    for (p, q, l, u, a, b) in cases {
        let x = mma_single_variable(p, q, l, u, a, b);
        let reference = grid_minimiser(p, q, l, u, a, b);
        assert!(
            (x - reference).abs() < 1e-4,
            "p={p} q={q} L={l} U={u}: got {x}, brute force {reference}"
        );
        assert!(
            approximation(p, q, l, u, x) <= approximation(p, q, l, u, reference) + 1e-9,
            "p={p} q={q}: the returned point must be a minimiser"
        );
        assert!(x >= a && x <= b);
    }
}

#[test]
fn single_variable_uses_l_plus_r_u_over_one_plus_r() {
    // r = sqrt(q/p) = 2, L = 0, U = 3:  x* = (L + r U) / (1 + r) = 2
    // (the swapped formula (L r + U) / (1 + r) would give 1).
    let x = mma_single_variable(1.0, 4.0, 0.0, 3.0, 0.1, 2.9);
    assert!((x - 2.0).abs() < 1e-12, "x = {x}");
    // equal p and q ⇒ the centre of the asymptotes
    let x = mma_single_variable(5.0, 5.0, -1.0, 3.0, -0.5, 2.5);
    assert!((x - 1.0).abs() < 1e-12, "x = {x}");
    // larger q/p pushes the minimiser away from the lower asymptote
    let lo = mma_single_variable(1.0, 1.0, 0.0, 3.0, 0.1, 2.9);
    let hi = mma_single_variable(1.0, 9.0, 0.0, 3.0, 0.1, 2.9);
    assert!(hi > lo, "{hi} should exceed {lo}");
}

#[test]
fn single_variable_one_sided_sensitivities_hit_the_bounds() {
    // dc < 0 (p = 0): q/(x − L) decreases with x ⇒ upper bound, finite (was NaN).
    let x = mma_single_variable(0.0, 2.0, 0.0, 1.0, 0.2, 0.8);
    assert_eq!(x, 0.8);
    // dc > 0 (q = 0): p/(U − x) increases with x ⇒ lower bound.
    let x = mma_single_variable(2.0, 0.0, 0.0, 1.0, 0.2, 0.8);
    assert_eq!(x, 0.2);
    // no information at all: centre of the asymptotes
    let x = mma_single_variable(0.0, 0.0, 0.0, 1.0, 0.2, 0.8);
    assert_eq!(x, 0.5);
}

#[test]
fn single_variable_is_always_finite_and_inside_its_interval() {
    let coefficients = [
        0.0,
        1e-300,
        1e-31,
        1e-29,
        1e-3,
        1.0,
        1e3,
        1e300,
        f64::INFINITY,
        f64::NAN,
        -1.0,
    ];
    for &p in &coefficients {
        for &q in &coefficients {
            let x = mma_single_variable(p, q, 0.0, 1.0, 0.25, 0.75);
            assert!(
                x.is_finite() && (0.25..=0.75).contains(&x),
                "p={p} q={q}: x = {x}"
            );
        }
    }
    // degenerate asymptote range and an empty interval never panic
    assert!(mma_single_variable(1.0, 1.0, 0.5, 0.5, 0.2, 0.8).is_finite());
    assert_eq!(mma_single_variable(1.0, 1.0, 0.0, 1.0, 0.8, 0.2), 0.2);
}

// ── asymptote update ──────────────────────────────────────────────────────────

#[test]
fn asymptotes_evolve_geometrically_with_the_sign_change_rule() {
    let mut asym = MmaAsymptote::new(0.5, 0.5);
    assert_eq!((asym.lower, asym.upper), (0.0, 1.0));
    // first move: no history, γ = 1
    asym.update_asymptotes(0.6, 0.0, 1.0, 1.2, 0.7);
    assert!((asym.lower - 0.1).abs() < 1e-12 && (asym.upper - 1.1).abs() < 1e-12);
    // same direction again: γ = 1.2 expands the distances
    asym.update_asymptotes(0.7, 0.0, 1.0, 1.2, 0.7);
    assert!((asym.lower - 0.1).abs() < 1e-12 && (asym.upper - 1.3).abs() < 1e-12);
    // reversal: γ = 0.7 contracts them
    asym.update_asymptotes(0.65, 0.0, 1.0, 1.2, 0.7);
    assert!((asym.lower - 0.23).abs() < 1e-12 && (asym.upper - 1.07).abs() < 1e-12);
    assert_eq!((asym.x_prev, asym.x_prev2), (0.65, 0.7));
}

#[test]
fn asymptote_distances_stay_within_svanberg_limits() {
    let mut asym = MmaAsymptote::new(0.5, 0.5);
    let mut x = 0.5;
    for k in 0..200 {
        x += if k % 2 == 0 { 0.01 } else { -0.01 };
        asym.update_asymptotes(x, 0.0, 1.0, 1.2, 0.7);
        assert!(x - asym.lower >= 0.01 - 1e-12 && x - asym.lower <= 10.0 + 1e-12);
        assert!(asym.upper - x >= 0.01 - 1e-12 && asym.upper - x <= 10.0 + 1e-12);
    }
    // persistent oscillation has contracted the asymptotes to the minimum distance
    assert!((x - asym.lower - 0.01).abs() < 1e-9);
    let mut drift = MmaAsymptote::new(0.5, 0.5);
    let mut y = 0.5;
    for _ in 0..200 {
        y += 0.001;
        drift.update_asymptotes(y, 0.0, 1.0, 1.2, 0.7);
    }
    assert!(y - drift.lower <= 10.0 + 1e-12 && drift.upper - y <= 10.0 + 1e-12);
}

// ── multi-variable step with the volume constraint ───────────────────────────

#[test]
fn negative_sensitivities_give_a_finite_increasing_design() {
    let n = 10;
    let x = vec![0.5; n];
    let dc: Vec<f64> = (0..n).map(|i| -(1.0 + i as f64)).collect();
    let mut asym = fresh_asymptotes(n, 0.5);
    let x_new = mma_step(&x, &dc, &mut asym, 0.0, 1.0, 0.2, 1.0);
    for (i, &v) in x_new.iter().enumerate() {
        assert!(v.is_finite(), "x[{i}] = {v}");
        assert!(v > 0.5, "x[{i}] = {v} must increase for dc < 0");
        assert!(v <= 0.7 + 1e-12, "x[{i}] = {v} exceeds the move limit");
    }
}

#[test]
fn positive_sensitivities_give_a_decreasing_design() {
    let n = 6;
    let x = vec![0.5; n];
    let dc: Vec<f64> = (0..n).map(|i| 0.5 + i as f64).collect();
    let mut asym = fresh_asymptotes(n, 0.5);
    let x_new = mma_step(&x, &dc, &mut asym, 0.0, 1.0, 0.2, 1.0);
    for (i, &v) in x_new.iter().enumerate() {
        assert!(
            v.is_finite() && v < 0.5,
            "x[{i}] = {v} must decrease for dc > 0"
        );
        assert!(v >= 0.3 - 1e-12);
    }
}

#[test]
fn mixed_signs_move_every_variable_along_its_descent_direction() {
    let x = vec![0.5; 6];
    let dc = [-1.0, 1.0, -2.0, 2.0, 0.0, -0.5];
    let mut asym = fresh_asymptotes(6, 0.5);
    let x_new = mma_step(&x, &dc, &mut asym, 0.0, 1.0, 0.2, 1.0);
    assert!(x_new[0] > 0.5 && x_new[2] > 0.5 && x_new[5] > 0.5);
    assert!(x_new[1] < 0.5 && x_new[3] < 0.5);
    // zero sensitivity: the variable stays where it is
    assert!((x_new[4] - 0.5).abs() < 1e-9, "x[4] = {}", x_new[4]);
}

#[test]
fn volume_constraint_is_met_exactly_when_active() {
    let n = 40;
    let x = vec![0.5; n];
    // Compliance-like sensitivities (all negative): the unconstrained step would
    // add material everywhere, so the constraint is active.
    let dc: Vec<f64> = (0..n).map(|i| -(1.0 + 0.1 * i as f64)).collect();
    let mut asym = fresh_asymptotes(n, 0.5);
    let x_new = mma_step(&x, &dc, &mut asym, 0.0, 1.0, 0.2, 0.5);
    let vf = volume_fraction(&x_new);
    assert!((vf - 0.5).abs() < 1e-6, "volume fraction {vf}");
    assert!(
        vf <= 0.5 + 1e-12,
        "the constraint must not be violated: {vf}"
    );
    // stiffer (more negative dc) elements receive more material
    assert!(x_new.windows(2).all(|w| w[1] >= w[0] - 1e-12));
    assert!(x_new[n - 1] > x_new[0] + 1e-3);
    // all values stay inside the trust region
    assert!(
        x_new
            .iter()
            .all(|&v| (0.3 - 1e-12..=0.7 + 1e-12).contains(&v))
    );
}

#[test]
fn slack_volume_constraint_leaves_the_unconstrained_step_untouched() {
    let n = 40;
    let x = vec![0.5; n];
    let dc: Vec<f64> = (0..n).map(|i| -(1.0 + 0.1 * i as f64)).collect();
    let mut asym = fresh_asymptotes(n, 0.5);
    let x_new = mma_step(&x, &dc, &mut asym, 0.0, 1.0, 0.2, 0.9);
    // every variable moves to its move limit; mean 0.7 < 0.9
    assert!(x_new.iter().all(|&v| (v - 0.7).abs() < 1e-12));
}

#[test]
fn unreachable_volume_returns_the_lowest_design_within_the_move_limits() {
    let n = 12;
    let x = vec![0.5; n];
    let dc: Vec<f64> = (0..n).map(|i| -(1.0 + i as f64)).collect();
    let mut asym = fresh_asymptotes(n, 0.5);
    let x_new = mma_step(&x, &dc, &mut asym, 0.0, 1.0, 0.2, 0.1);
    assert!(
        x_new
            .iter()
            .all(|&v| v.is_finite() && (v - 0.3).abs() < 1e-12)
    );
}

#[test]
fn step_updates_the_asymptotes_around_the_new_design() {
    let n = 4;
    let x = vec![0.5; n];
    let dc = [-1.0, -2.0, 1.0, 2.0];
    let mut asym = fresh_asymptotes(n, 0.5);
    let x_new = mma_step(&x, &dc, &mut asym, 0.0, 1.0, 0.2, 1.0);
    for (a, &xn) in asym.iter().zip(x_new.iter()) {
        assert_eq!(a.x_prev, xn);
        assert_eq!(a.x_prev2, 0.5);
        assert!(a.lower < xn && xn < a.upper);
    }
}

#[test]
fn step_survives_asymptotes_that_do_not_bracket_the_design() {
    let x = vec![0.5, 0.2];
    let dc = [-1.0, 3.0];
    let mut asym = vec![
        MmaAsymptote {
            lower: 0.6,
            upper: 0.7,
            x_prev: 0.5,
            x_prev2: 0.5,
        },
        MmaAsymptote {
            lower: f64::NAN,
            upper: 0.9,
            x_prev: 0.2,
            x_prev2: 0.2,
        },
    ];
    let x_new = mma_step(&x, &dc, &mut asym, 0.0, 1.0, 0.2, 1.0);
    assert!(x_new.iter().all(|v| v.is_finite()));
    assert!(x_new[0] > 0.5 && x_new[1] < 0.2);
}

#[test]
fn step_handles_degenerate_inputs_without_panicking() {
    // empty design
    let mut none: Vec<MmaAsymptote> = Vec::new();
    assert!(mma_step(&[], &[], &mut none, 0.0, 1.0, 0.2, 0.5).is_empty());
    // inverted bounds leave the design untouched
    let mut asym = fresh_asymptotes(2, 0.5);
    let x = mma_step(&[0.4, 0.6], &[-1.0, 1.0], &mut asym, 1.0, 0.0, 0.2, 0.5);
    assert_eq!(x, vec![0.4, 0.6]);
    // a design outside the bounds yields finite values
    let mut asym = fresh_asymptotes(2, 1.5);
    let x = mma_step(&[1.5, -0.5], &[-1.0, 1.0], &mut asym, 0.0, 1.0, 0.2, 1.0);
    assert!(x.iter().all(|v| v.is_finite()));
}

// ── optimisation of a separable quadratic ─────────────────────────────────────

fn quadratic_problem() -> (Vec<f64>, Vec<f64>) {
    let n = 10;
    let c: Vec<f64> = (0..n).map(|i| 1.0 + i as f64).collect();
    let t: Vec<f64> = (0..n).map(|i| 0.1 + 0.08 * i as f64).collect();
    (c, t)
}

fn quadratic_value(c: &[f64], t: &[f64], x: &[f64]) -> f64 {
    c.iter()
        .zip(t.iter())
        .zip(x.iter())
        .map(|((&ci, &ti), &xi)| ci * (xi - ti).powi(2))
        .sum()
}

fn quadratic_gradient(c: &[f64], t: &[f64], x: &[f64]) -> Vec<f64> {
    c.iter()
        .zip(t.iter())
        .zip(x.iter())
        .map(|((&ci, &ti), &xi)| 2.0 * ci * (xi - ti))
        .collect()
}

#[test]
fn separable_quadratic_objective_decreases() {
    let (c, t) = quadratic_problem();
    let n = c.len();
    let mut x = vec![0.5; n];
    let mut asym = fresh_asymptotes(n, 0.5);
    let f0 = quadratic_value(&c, &t, &x);
    let mut f_after_first_step = f0;
    for k in 0..80 {
        let dc = quadratic_gradient(&c, &t, &x);
        x = mma_step(&x, &dc, &mut asym, 0.0, 1.0, 0.2, 1.0);
        assert!(x.iter().all(|v| v.is_finite() && (0.0..=1.0).contains(v)));
        if k == 0 {
            f_after_first_step = quadratic_value(&c, &t, &x);
        }
    }
    assert!(
        f_after_first_step < f0,
        "first step must reduce the objective: {f_after_first_step} vs {f0}"
    );
    let f_final = quadratic_value(&c, &t, &x);
    assert!(
        f_final < 1.0e-2 * f0,
        "objective only fell from {f0} to {f_final}"
    );
}

#[test]
fn volume_constrained_quadratic_converges_to_the_kkt_point() {
    // minimise Σ cᵢ (xᵢ − tᵢ)²  s.t. mean(x) = vf, 0 ≤ x ≤ 1:
    //   xᵢ = clamp(tᵢ − μ / (2 cᵢ), 0, 1) with μ fixed by the volume.
    let (c, t) = quadratic_problem();
    let n = c.len();
    let vf = 0.35;
    let design = |mu: f64| -> Vec<f64> {
        c.iter()
            .zip(t.iter())
            .map(|(&ci, &ti)| (ti - mu / (2.0 * ci)).clamp(0.0, 1.0))
            .collect()
    };
    let (mut lo, mut hi) = (-100.0_f64, 100.0_f64);
    for _ in 0..200 {
        let mu = 0.5 * (lo + hi);
        if volume_fraction(&design(mu)) > vf {
            lo = mu;
        } else {
            hi = mu;
        }
    }
    let kkt = design(0.5 * (lo + hi));
    assert!((volume_fraction(&kkt) - vf).abs() < 1e-12);

    let mut x = vec![0.5; n];
    let mut asym = fresh_asymptotes(n, 0.5);
    for _ in 0..60 {
        let dc = quadratic_gradient(&c, &t, &x);
        x = mma_step(&x, &dc, &mut asym, 0.0, 1.0, 0.2, vf);
        assert!((volume_fraction(&x) - vf).abs() < 1e-6);
    }
    for i in 0..n {
        assert!(
            (x[i] - kkt[i]).abs() < 1e-6,
            "variable {i}: MMA {} vs KKT {}",
            x[i],
            kkt[i]
        );
    }
}
