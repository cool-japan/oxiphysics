//! Regression tests for the variable-mobility Cahn-Hilliard update.
//!
//! The update used to be `C += M(C) * laplacian(mu)`, which is not `div(M grad mu)`: the total
//! concentration drifted for every non-constant mobility. The degenerate mobility was also
//! written for a `[0, 1]` order parameter (`C (1 - C)` on a clamped `C`) although the double
//! well puts the phases at `C = ±1`, which froze the whole `C < 0` phase and the interface
//! centre. These tests pin the conservative flux form, its second-order accuracy and the
//! corrected mobility functions.

use super::*;

/// Next pseudo-random number in `[-1, 1)` of a 64-bit LCG.
fn next_unit(state: &mut u64) -> f64 {
    *state = state
        .wrapping_mul(6_364_136_223_846_793_005)
        .wrapping_add(1_442_695_040_888_963_407);
    2.0 * ((*state >> 11) as f64 / (1u64 << 53) as f64) - 1.0
}

/// Solver filled with deterministic noise of the given amplitude around `mean`.
fn noisy_solver(mobility_type: MobilityType, mean: f64, noise: f64, seed: u64) -> CahnHilliardLbm {
    let mut solver = CahnHilliardLbm::new(16, 16, 1.0, 0.1, 0.1, mobility_type);
    let mut state = seed;
    for c in solver.concentration.iter_mut() {
        *c = mean + noise * next_unit(&mut state);
    }
    solver.compute_chemical_potential();
    solver
}

fn total(solver: &CahnHilliardLbm) -> f64 {
    solver.concentration.iter().sum()
}

const VARIABLE_TYPES: [MobilityType; 3] = [
    MobilityType::Degenerate,
    MobilityType::ConcentrationDependent,
    MobilityType::Polynomial,
];

const ALL_TYPES: [MobilityType; 4] = [
    MobilityType::Constant,
    MobilityType::Degenerate,
    MobilityType::ConcentrationDependent,
    MobilityType::Polynomial,
];

#[test]
fn degenerate_mobility_is_the_symmetric_interface_peak() {
    let m0 = 0.8;
    for k in -10..=10 {
        let c = f64::from(k) / 10.0;
        let m = eval_mobility(MobilityType::Degenerate, c, m0);
        assert!((m - m0 * (1.0 - c * c) / 4.0).abs() < 1e-15, "M({c}) = {m}");
        let mirrored = eval_mobility(MobilityType::Degenerate, -c, m0);
        assert!((m - mirrored).abs() < 1e-15, "M must be even in C");
    }
    // The -1 phase is mobile at its interface side: the old clamp(0, 1) returned zero here.
    assert!(eval_mobility(MobilityType::Degenerate, -0.5, m0) > 0.1);
    // Outside the physical range the mobility is clamped to zero, never negative.
    assert_eq!(eval_mobility(MobilityType::Degenerate, 1.3, m0), 0.0);
    assert_eq!(eval_mobility(MobilityType::Degenerate, -1.3, m0), 0.0);
}

#[test]
fn variable_mobilities_vanish_at_the_pure_phases_and_are_positive_between() {
    for mobility_type in VARIABLE_TYPES {
        for pure in [-1.0, 1.0] {
            assert!(
                eval_mobility(mobility_type, pure, 1.0).abs() < 1e-15,
                "{mobility_type:?} at C = {pure}"
            );
        }
        for c in [-0.9, -0.4, 0.0, 0.3, 0.8] {
            assert!(
                eval_mobility(mobility_type, c, 1.0) > 0.0,
                "{mobility_type:?} at C = {c}"
            );
        }
    }
    assert_eq!(eval_mobility(MobilityType::Constant, 1.0, 0.3), 0.3);
}

#[test]
fn total_concentration_is_conserved_for_every_mobility_type() {
    for (n, mobility_type) in ALL_TYPES.into_iter().enumerate() {
        let mut solver = noisy_solver(mobility_type, 0.0, 0.05, 11 + n as u64);
        let initial = solver.concentration.clone();
        let total0 = total(&solver);
        solver.run(200);
        let total1 = total(&solver);
        assert!(
            (total1 - total0).abs() < 1e-10,
            "{mobility_type:?}: total concentration {total0} -> {total1}"
        );
        assert!(
            solver
                .concentration
                .iter()
                .all(|c| c.is_finite() && c.abs() < 1.5),
            "{mobility_type:?}: solution left the physical range"
        );
        let moved = solver
            .concentration
            .iter()
            .zip(initial.iter())
            .map(|(a, b)| (a - b).abs())
            .fold(0.0_f64, f64::max);
        assert!(moved > 1e-3, "{mobility_type:?}: state did not evolve");
    }
}

#[test]
fn conservation_holds_around_a_nonzero_mean() {
    // Off-critical composition: the old pointwise update drifted most here.
    for mobility_type in VARIABLE_TYPES {
        let mut solver = noisy_solver(mobility_type, 0.3, 0.05, 5);
        let total0 = total(&solver);
        solver.run(150);
        assert!(
            (total(&solver) - total0).abs() < 1e-10,
            "{mobility_type:?}: drift {}",
            total(&solver) - total0
        );
    }
}

#[test]
fn step_applies_the_link_flux_divergence_of_the_face_mobility() {
    let mut solver = noisy_solver(MobilityType::Degenerate, 0.0, 0.4, 3);
    let (nx, ny) = (solver.nx, solver.ny);
    let c0 = solver.concentration.clone();
    let mu0 = solver.chemical_potential.clone();
    let mobility = |c: f64| eval_mobility(MobilityType::Degenerate, c, solver.mobility);

    // Independent reference: sum the fluxes through the eight links of every node.
    let mut expected = vec![0.0_f64; nx * ny];
    let mut pointwise = vec![0.0_f64; nx * ny];
    for i in 0..nx {
        for j in 0..ny {
            let idx = i * ny + j;
            for q in 1..9 {
                let ni = (i + nx).wrapping_add_signed(CX[q] as isize) % nx;
                let nj = (j + ny).wrapping_add_signed(CY[q] as isize) % ny;
                let nidx = ni * ny + nj;
                let face = 0.5 * (mobility(c0[idx]) + mobility(c0[nidx]));
                expected[idx] += 6.0 * W9[q] * face * (mu0[nidx] - mu0[idx]);
                pointwise[idx] += 6.0 * W9[q] * mobility(c0[idx]) * (mu0[nidx] - mu0[idx]);
            }
        }
    }

    solver.step();
    let mut max_gap_to_pointwise = 0.0_f64;
    for idx in 0..nx * ny {
        let applied = solver.concentration[idx] - c0[idx];
        assert!(
            (applied - expected[idx]).abs() < 1e-14,
            "node {idx}: applied {applied}, link-flux divergence {}",
            expected[idx]
        );
        max_gap_to_pointwise = max_gap_to_pointwise.max((applied - pointwise[idx]).abs());
    }
    assert!(
        max_gap_to_pointwise > 1e-5,
        "the test must distinguish div(M grad mu) from M laplacian(mu): {max_gap_to_pointwise}"
    );
}

#[test]
fn constant_mobility_reduces_to_m_times_the_laplacian_of_mu() {
    let mut solver = noisy_solver(MobilityType::Constant, 0.0, 0.3, 17);
    let (nx, ny) = (solver.nx, solver.ny);
    let c0 = solver.concentration.clone();
    let mu0 = solver.chemical_potential.clone();
    solver.step();
    for i in 0..nx {
        for j in 0..ny {
            let idx = i * ny + j;
            let mut lap = 0.0;
            for q in 1..9 {
                let ni = (i + nx).wrapping_add_signed(CX[q] as isize) % nx;
                let nj = (j + ny).wrapping_add_signed(CY[q] as isize) % ny;
                lap += 6.0 * W9[q] * (mu0[ni * ny + nj] - mu0[idx]);
            }
            let applied = solver.concentration[idx] - c0[idx];
            assert!((applied - solver.mobility * lap).abs() < 1e-14);
        }
    }
}

/// Analytic smooth fields on a periodic square of side `2 pi` and `div(M(C) grad mu)` for them.
struct SmoothFields {
    m0: f64,
}

impl SmoothFields {
    fn c(&self, x: f64, y: f64) -> f64 {
        0.5 * x.cos() * y.sin()
    }

    fn mu(&self, x: f64, y: f64) -> f64 {
        (x + 0.3).sin() * (2.0 * y).cos()
    }

    /// `div(M grad mu)` for the degenerate mobility `M = m0 (1 - c^2) / 4`.
    fn exact_divergence(&self, x: f64, y: f64) -> f64 {
        let c = self.c(x, y);
        let (c_x, c_y) = (-0.5 * x.sin() * y.sin(), 0.5 * x.cos() * y.cos());
        let m = 0.25 * self.m0 * (1.0 - c * c);
        let (m_x, m_y) = (-0.5 * self.m0 * c * c_x, -0.5 * self.m0 * c * c_y);
        let (mu_x, mu_y) = (
            (x + 0.3).cos() * (2.0 * y).cos(),
            -2.0 * (x + 0.3).sin() * (2.0 * y).sin(),
        );
        let laplacian = -5.0 * (x + 0.3).sin() * (2.0 * y).cos();
        m * laplacian + m_x * mu_x + m_y * mu_y
    }
}

/// Largest deviation of the lattice flux divergence from the analytic one on an `n x n` grid.
fn flux_divergence_error(n: usize) -> f64 {
    let fields = SmoothFields { m0: 0.8 };
    let mut solver = CahnHilliardLbm::new(n, n, 1.0, 0.1, fields.m0, MobilityType::Degenerate);
    let h = 2.0 * std::f64::consts::PI / n as f64;
    for i in 0..n {
        for j in 0..n {
            let idx = solver.idx(i, j);
            solver.concentration[idx] = fields.c(i as f64 * h, j as f64 * h);
            solver.chemical_potential[idx] = fields.mu(i as f64 * h, j as f64 * h);
        }
    }
    let mobility = solver.mobility_field();
    let mut worst = 0.0_f64;
    for i in 0..n {
        for j in 0..n {
            let lattice = solver.mobility_flux_divergence(&mobility, i, j) / (h * h);
            let exact = fields.exact_divergence(i as f64 * h, j as f64 * h);
            worst = worst.max((lattice - exact).abs());
        }
    }
    worst
}

#[test]
fn flux_divergence_is_second_order_accurate() {
    let e16 = flux_divergence_error(16);
    let e32 = flux_divergence_error(32);
    let e64 = flux_divergence_error(64);
    assert!(e64 < 1e-2, "error on the finest grid: {e64}");
    for (coarse, fine, label) in [(e16, e32, "16 -> 32"), (e32, e64, "32 -> 64")] {
        let order = (coarse / fine).log2();
        assert!(
            (1.8..2.3).contains(&order),
            "observed order {order} for {label} ({coarse} -> {fine})"
        );
    }
}

#[test]
fn small_perturbation_grows_with_the_exact_discrete_linear_rate() {
    // For c = eps cos(k x) about C = 0 the update is linear:
    // c' = c (1 + M(0) lambda (-A - kappa lambda)) with lambda = 2 (cos k - 1) the eigenvalue of
    // the D2Q9 Laplacian stencil. With the old degenerate mobility M(0) was zero and nothing grew.
    let (a_coeff, kappa, m0) = (1.0, 0.1, 0.1);
    let steps = 30;
    let eps = 1e-6;
    for mobility_type in ALL_TYPES {
        for mode in [1usize, 3] {
            let mut solver = CahnHilliardLbm::new(16, 16, a_coeff, kappa, m0, mobility_type);
            let k = 2.0 * std::f64::consts::PI * mode as f64 / 16.0;
            for i in 0..16 {
                for j in 0..16 {
                    let idx = solver.idx(i, j);
                    solver.concentration[idx] = eps * (k * i as f64).cos();
                }
            }
            let amplitude = |s: &CahnHilliardLbm| {
                let mut acc = 0.0;
                for i in 0..16 {
                    for j in 0..16 {
                        acc += s.concentration[s.idx(i, j)] * (k * i as f64).cos();
                    }
                }
                2.0 * acc / 256.0
            };
            let before = amplitude(&solver);
            solver.run(steps);
            let measured = amplitude(&solver) / before;

            let lambda = 2.0 * (k.cos() - 1.0);
            let mobility_at_zero = eval_mobility(mobility_type, 0.0, m0);
            let factor = 1.0 + mobility_at_zero * lambda * (-a_coeff - kappa * lambda);
            let expected = factor.powi(steps as i32);
            assert!(
                (measured - expected).abs() < 1e-6 * expected.abs(),
                "{mobility_type:?} mode {mode}: measured {measured}, expected {expected}"
            );
            if mobility_at_zero > 0.0 {
                assert!(measured > 1.05, "{mobility_type:?} mode {mode} must grow");
            }
        }
    }
}

#[test]
fn negative_phase_side_is_mobile_with_the_degenerate_mobility() {
    let mut solver = CahnHilliardLbm::new(12, 12, 1.0, 0.1, 0.1, MobilityType::Degenerate);
    for c in solver.concentration.iter_mut() {
        *c = -0.5;
    }
    let bump = solver.idx(6, 6);
    solver.concentration[bump] = -0.3;
    let before = solver.concentration.clone();
    solver.step();
    let east = solver.idx(7, 6);
    assert!(
        (solver.concentration[east] - before[east]).abs() > 1e-5,
        "the C < 0 phase must respond to the chemical-potential gradient"
    );
    assert!(
        (solver.concentration[bump] - before[bump]).abs() > 1e-5,
        "the perturbed node must relax"
    );
}

#[test]
fn pure_phases_are_stationary_for_every_mobility_type() {
    for mobility_type in ALL_TYPES {
        for pure in [-1.0, 1.0] {
            let mut solver = CahnHilliardLbm::new(8, 8, 1.0, 0.1, 0.1, mobility_type);
            for c in solver.concentration.iter_mut() {
                *c = pure;
            }
            solver.run(10);
            assert!(
                solver
                    .concentration
                    .iter()
                    .all(|c| (c - pure).abs() < 1e-14),
                "{mobility_type:?}: uniform C = {pure} must be a fixed point"
            );
        }
    }
}
