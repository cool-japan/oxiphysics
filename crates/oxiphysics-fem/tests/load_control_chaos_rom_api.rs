// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Public-API regression tests for the load-control solvers, the polynomial
//! chaos projection, the Monte Carlo generator, the reduced-order Newmark
//! integrator and the Biot coupling matrix.
//!
//! They use only items exported by the crate, from outside of it, so a change
//! of a path or a signature shows up here.

use oxiphysics_fem::analysis::FemError;
use oxiphysics_fem::poromechanics::{BiotBlock, ConsolidationSolver, PorousElement};
use oxiphysics_fem::probabilistic_fem::{
    MonteCarloFem, NormalQuadrature, SpectralStochasticFem, StochasticFem, pce_project,
};
use oxiphysics_fem::reduced_order::{
    PodBasis, RomSystem, balanced_truncation, hankel_singular_values, solve_lyapunov,
};
use oxiphysics_fem::{
    ArcLengthMethod, LoadStepper, LoadedResidualFn, NewtonRaphson, NonlinearSpring, ResidualFn,
};

/// User-defined snap-through spring `f_int(u) = u (u - 3)²` loaded by `f_ref = 1`.
struct SnapThrough;

impl ResidualFn for SnapThrough {
    fn residual(&self, u: &[f64]) -> Vec<f64> {
        vec![u[0] * (u[0] - 3.0).powi(2) - 1.0]
    }

    fn tangent_stiffness(&self, u: &[f64]) -> Vec<Vec<f64>> {
        vec![vec![3.0 * u[0] * u[0] - 12.0 * u[0] + 9.0]]
    }
}

impl LoadedResidualFn for SnapThrough {
    fn internal_force(&self, u: &[f64]) -> Vec<f64> {
        vec![u[0] * (u[0] - 3.0).powi(2)]
    }

    fn reference_load(&self) -> Vec<f64> {
        vec![1.0]
    }
}

#[test]
fn arc_length_follows_a_user_problem_through_its_limit_load() {
    let arc = ArcLengthMethod::new(0.25, 40, 1e-10);
    let path = arc.trace_path(&SnapThrough, &[0.0], 0.0, 45);
    assert_eq!(path.len(), 45);
    for (u, lambda) in &path {
        assert!(SnapThrough.residual_at(u, *lambda)[0].abs() < 1e-9);
    }
    // The load rises to the limit load 4 and then falls again.
    let first_decrease = path
        .windows(2)
        .position(|pair| pair[1].1 < pair[0].1)
        .expect("the load decreases after the limit point");
    assert!((path[first_decrease].1 - 4.0).abs() < 0.1);
}

#[test]
fn load_stepper_solves_the_scaled_problem_and_reports_failures() {
    let spring = NonlinearSpring {
        k1: 1.0,
        k3: 1.0,
        force: 1.0,
    };
    let results = LoadStepper::new(3, 0.0, 1.0)
        .solve(&spring, &[0.0], &NewtonRaphson::new(50, 1e-12))
        .expect("converges");
    // u + u^3 = 0.5 at half load, not the full-load root 0.6823.
    assert!((results[1].1[0] - 0.4239).abs() < 5e-5);
    assert!((results[2].1[0] - 0.6823).abs() < 5e-5);

    let hard = NonlinearSpring {
        k1: 1.0,
        k3: 1.0,
        force: 10.0,
    };
    let failure = LoadStepper::new(2, 0.5, 1.0).solve(&hard, &[0.0], &NewtonRaphson::new(1, 1e-14));
    assert!(failure.is_err());
}

#[test]
fn chaos_projection_recovers_cross_terms_exactly() {
    // f = 1 + 2 xi_1 + 3 xi_1 xi_2:  coefficients 1, 2, 3 and zero elsewhere.
    let alphas = vec![
        vec![0, 0],
        vec![1, 0],
        vec![0, 1],
        vec![2, 0],
        vec![1, 1],
        vec![0, 2],
    ];
    let rule = NormalQuadrature::tensor(2, 3).expect("tensor rule");
    let coefficients = pce_project(|xi| 1.0 + 2.0 * xi[0] + 3.0 * xi[0] * xi[1], &alphas, &rule)
        .expect("projection");
    let expected = [1.0, 2.0, 0.0, 0.0, 3.0, 0.0];
    for (c, e) in coefficients.iter().zip(expected.iter()) {
        assert!((c - e).abs() < 1e-13, "{coefficients:?}");
    }
}

#[test]
fn spectral_stochastic_fem_projects_on_the_full_grid_and_rejects_impossible_dimensions() {
    let mut pce = SpectralStochasticFem::new(2, 2, 8, 1.0, 1.0, 100.0, 1e6);
    pce.compute_pce_coefficients(0.05).expect("projection");
    assert_eq!(pce.u_tip_pce.len(), 6);
    assert!((pce.pce_mean() - 1e-4).abs() < 1e-6);
    assert!(pce.pce_variance() > 0.0);

    // Five random dimensions cannot be represented on four elements.
    let mut impossible = SpectralStochasticFem::new(5, 1, 4, 1.0, 1.0, 100.0, 1e6);
    assert!(impossible.compute_pce_coefficients(0.05).is_err());
}

#[test]
fn monte_carlo_samples_are_reproducible_and_unbiased() {
    let build = || StochasticFem::new(8, 1.0, 1.0, 100.0, 1e6, 0.05, 3);
    let mut a = MonteCarloFem::with_seed(build(), 4_000, 11);
    let mut b = MonteCarloFem::with_seed(build(), 4_000, 11);
    a.run();
    b.run();
    assert_eq!(a.samples, b.samples);
    // Mean tip displacement F L / (E A) = 1e-4 up to the second-order bias.
    assert!((a.mean_tip() / 1e-4 - 1.0).abs() < 0.01);
    assert!(a.std_tip() > 0.0);
}

#[test]
fn reduced_order_newmark_integrates_a_free_vibration() {
    let basis = PodBasis {
        modes: vec![vec![1.0]],
        singular_values: vec![1.0],
        n_modes: 1,
    };
    let omega = 2.0 * std::f64::consts::PI;
    let mut rom = RomSystem::project_system(&basis, &[vec![1.0]], &[vec![omega * omega]]);
    rom.set_initial_state(&[1.0], &[0.0], &[0.0])
        .expect("initial state");
    let dt = 1.0 / 400.0;
    let mut last = 1.0;
    for _ in 0..100 {
        last = rom.step_newmark(&[0.0], dt, 0.25, 0.5).expect("step")[0];
    }
    // A quarter period: q = cos(pi / 2) = 0 up to the period elongation of the scheme.
    assert!(last.abs() < 1e-3, "q after a quarter period: {last}");
    assert!(matches!(
        rom.step_newmark(&[0.0], 0.0, 0.25, 0.5),
        Err(FemError::InvalidInput(_))
    ));
}

#[test]
fn biot_coupling_matrix_feeds_the_consolidation_solver() {
    // One tetrahedron, all nodes drained and the base held: the coupling block
    // has the shape the solver expects and responds to the element geometry.
    let element = PorousElement::new(vec![0, 1, 2, 3], [0.0; 9], 0.3, 1e7);
    let coords = [
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
    ];
    let q = element
        .compute_coupling_matrix(&coords, 1.0)
        .expect("coupling matrix");
    assert_eq!(q.len(), 12 * 4);
    let u_dofs: Vec<usize> = (0..12).collect();
    let p_dofs: Vec<usize> = (0..4).collect();
    let mut solver = ConsolidationSolver::new(12, 4);
    solver
        .add_element_block(BiotBlock::Coupling, &u_dofs, &p_dofs, &q)
        .expect("scatter");
    // Node 0 couples through all three components (gradient (-1, -1, -1)), nodes 1
    // to 3 through one component each (unit gradients): 3·4 + 3·4 non-zero entries.
    assert_eq!(solver.block(BiotBlock::Coupling).nnz(), 24);
    let bigger = coords.map(|p| [2.0 * p[0], 2.0 * p[1], 2.0 * p[2]]);
    let q_bigger = element
        .compute_coupling_matrix(&bigger, 1.0)
        .expect("coupling matrix");
    assert!((q_bigger[3 * 4] - 4.0 * q[3 * 4]).abs() < 1e-14);
}

#[test]
fn balanced_truncation_is_driven_by_the_lyapunov_gramians() {
    // x' = -2 x + u, y = 3 x:  Wc = 1/4, Wo = 9/4, so the Hankel singular value is 3/4.
    let hsv = hankel_singular_values(&[vec![-2.0]], &[vec![1.0]], &[vec![3.0]]).expect("stable");
    assert_eq!(hsv.len(), 1);
    assert!((hsv[0] - 0.75).abs() < 1e-13, "sigma = {}", hsv[0]);

    // A W + W A^T + I = 0 for the oscillator with poles -0.3 +- 5i has W = I / 0.6.
    let oscillator = vec![vec![-0.3, 5.0], vec![-5.0, -0.3]];
    let identity = vec![vec![1.0, 0.0], vec![0.0, 1.0]];
    let w = solve_lyapunov(&oscillator, &identity).expect("stable");
    assert!((w[0][0] - 1.0 / 0.6).abs() < 1e-12 && (w[1][1] - 1.0 / 0.6).abs() < 1e-12);
    assert!(w[0][1].abs() < 1e-12 && w[1][0].abs() < 1e-12);

    // Two real poles (-1 and -2) driven and observed through both states: the one-state
    // model stays stable and its static gain is within the bound 2 sigma_2 of the full one.
    let a = vec![vec![-1.0, 0.0], vec![0.0, -2.0]];
    let b = vec![vec![1.0], vec![1.0]];
    let c = vec![vec![1.0, 1.0]];
    let sigma = hankel_singular_values(&a, &b, &c).expect("stable");
    let (a_r, b_r, c_r) = balanced_truncation(&a, &b, &c, 1).expect("stable");
    assert_eq!((a_r.len(), b_r.len(), c_r[0].len()), (1, 1, 1));
    assert!(a_r[0][0] < 0.0, "the reduced model must stay stable");
    let full_gain = 1.0 + 0.5;
    let reduced_gain = -c_r[0][0] * b_r[0][0] / a_r[0][0];
    assert!((full_gain - reduced_gain).abs() <= 2.0 * sigma[1] + 1e-12);

    // An unstable state matrix is rejected instead of yielding fabricated Gramians.
    assert!(balanced_truncation(&[vec![0.5]], &[vec![1.0]], &[vec![1.0]], 1).is_err());
}
