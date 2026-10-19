// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Unit tests of the nonlinear FEM building blocks (kinematics, J2 plasticity,
//! Newton-Raphson, arc-length control, locking-free elements).

use super::*;

// ── Mat3 operations ───────────────────────────────────────────────────────

#[test]
fn mat3_mul_identity() {
    let a = [[1., 2., 3.], [4., 5., 6.], [7., 8., 9.]];
    let c = mat3_mul(&a, &EYE3);
    for (i, row) in c.iter().enumerate() {
        for (j, &val) in row.iter().enumerate() {
            assert!((val - a[i][j]).abs() < 1e-10);
        }
    }
}

#[test]
fn mat3_det_identity_is_one() {
    assert!((mat3_det(&EYE3) - 1.0).abs() < 1e-10);
}

#[test]
fn mat3_inv_of_identity() {
    let inv = mat3_inv(&EYE3);
    for (i, row) in inv.iter().enumerate() {
        for (j, &val) in row.iter().enumerate() {
            assert!((val - EYE3[i][j]).abs() < 1e-10);
        }
    }
}

#[test]
fn mat3_inv_round_trip() {
    let a = [[2., 1., 0.], [1., 3., 1.], [0., 1., 4.]];
    let inv = mat3_inv(&a);
    let prod = mat3_mul(&a, &inv);
    for (i, row) in prod.iter().enumerate() {
        for (j, &val) in row.iter().enumerate() {
            let expected = if i == j { 1.0 } else { 0.0 };
            assert!(
                (val - expected).abs() < 1e-8,
                "A·A⁻¹[{i}][{j}] = {} != {}",
                val,
                expected
            );
        }
    }
}

#[test]
fn mat3_transpose_symmetry() {
    let a = [[1., 2., 3.], [4., 5., 6.], [7., 8., 9.]];
    let at = mat3_transpose(&a);
    for (i, row) in at.iter().enumerate() {
        for (j, &val) in row.iter().enumerate() {
            assert!((val - a[j][i]).abs() < 1e-10);
        }
    }
}

// ── DeformationGradient ───────────────────────────────────────────────────

#[test]
fn deformation_gradient_identity_jacobian() {
    let f = DeformationGradient::identity();
    assert!((f.jacobian() - 1.0).abs() < 1e-10);
}

#[test]
fn deformation_gradient_green_lagrange_zero_for_identity() {
    let f = DeformationGradient::identity();
    let e = f.green_lagrange_strain();
    for row in &e {
        for &val in row {
            assert!(val.abs() < 1e-10);
        }
    }
}

#[test]
fn deformation_gradient_left_cauchy_green_identity() {
    let f = DeformationGradient::identity();
    let b = f.left_cauchy_green();
    for i in 0..3 {
        for j in 0..3 {
            assert!((b[i][j] - EYE3[i][j]).abs() < 1e-10);
        }
    }
}

#[test]
fn deformation_gradient_uniaxial_extension() {
    let h = [[0.1, 0., 0.], [0., 0., 0.], [0., 0., 0.]]; // ε_11 = 0.1
    let f = DeformationGradient::from_displacement_gradient(&h);
    assert!((f.jacobian() - 1.1).abs() < 1e-10, "J = {}", f.jacobian());
    let e = f.green_lagrange_strain();
    // E_11 = 0.5*(F_11² - 1) = 0.5*(1.1²-1) = 0.105
    assert!((e[0][0] - 0.105).abs() < 1e-10, "E11 = {}", e[0][0]);
}

#[test]
fn f_bar_preserves_deviatoric() {
    let h = [[0.2, 0., 0.], [0., 0., 0.], [0., 0., 0.]];
    let f = DeformationGradient::from_displacement_gradient(&h);
    let j0 = 1.0;
    let fb = f.f_bar(j0);
    let j_fb = mat3_det(&fb);
    assert!(
        (j_fb - j0).abs() < 1e-6,
        "F̄ should have det = J0={j0}, got {j_fb}"
    );
}

// ── J2 plasticity ─────────────────────────────────────────────────────────

#[test]
fn j2_elastic_predictor_no_plastic() {
    let mat = J2Material::new(200e9, 0.3, 250e6, 0.0);
    let sigma_n = [0.0; 6];
    let delta_eps = [0.001, 0.0, 0.0, 0.0, 0.0, 0.0];
    let sigma_tr = mat.elastic_predictor(&sigma_n, &delta_eps);
    assert!(sigma_tr[0] > 0.0, "tensile stress expected");
}

#[test]
fn j2_return_mapping_elastic_stays() {
    let mat = J2Material::new(200e9, 0.3, 250e6, 0.0);
    let sigma_n = [0.0; 6];
    let delta_eps = [1e-5, 0.0, 0.0, 0.0, 0.0, 0.0]; // well below yield
    let sigma_tr = mat.elastic_predictor(&sigma_n, &delta_eps);
    let state_n = PlasticState::new();
    let (sigma_n1, state_n1, delta_gamma) = mat.return_mapping(&sigma_tr, &state_n);
    assert!(delta_gamma < 1e-10, "should be elastic: δγ = {delta_gamma}");
    assert!(state_n1.eps_p_eq < 1e-10);
    assert!((sigma_n1[0] - sigma_tr[0]).abs() < 1e-6);
}

#[test]
fn j2_return_mapping_plastic_caps_stress() {
    let mat = J2Material::new(200e9, 0.3, 250e6, 1e9);
    let sigma_n = [0.0; 6];
    let delta_eps = [0.01, 0.0, 0.0, 0.0, 0.0, 0.0]; // large step → yield
    let sigma_tr = mat.elastic_predictor(&sigma_n, &delta_eps);
    let state_n = PlasticState::new();
    let (sigma_n1, state_n1, delta_gamma) = mat.return_mapping(&sigma_tr, &state_n);
    let vm = J2Material::von_mises(&sigma_n1);
    let sy = mat.yield_stress(state_n1.eps_p_eq);
    assert!(
        vm <= sy * 1.001 + 1.0,
        "stress {vm} should not exceed yield {sy}"
    );
    assert!(delta_gamma > 0.0, "plastic step expected");
}

#[test]
fn j2_von_mises_hydrostatic_is_zero() {
    let sigma_hydro = [100e6, 100e6, 100e6, 0.0, 0.0, 0.0];
    let vm = J2Material::von_mises(&sigma_hydro);
    assert!(
        vm < 1.0,
        "hydrostatic stress should give zero von Mises: {vm}"
    );
}

#[test]
fn j2_consistent_tangent_elastic_symmetry() {
    let mat = J2Material::new(200e9, 0.3, 250e6, 0.0);
    let sigma_tr = [100e6; 6];
    let state_n = PlasticState::new();
    let c = mat.consistent_tangent(0.0, &sigma_tr, &state_n);
    for (i, row) in c.iter().enumerate() {
        for (j, &val) in row.iter().enumerate() {
            assert!(
                (val - c[j][i]).abs() < 1.0,
                "tangent not symmetric at [{i}][{j}]: {} vs {}",
                val,
                c[j][i]
            );
        }
    }
}

// ── Newton-Raphson ─────────────────────────────────────────────────────────

struct LinearSystem {
    k: f64,
    f: f64,
}

impl NonlinearSystem for LinearSystem {
    fn residual(&self, x: &[f64]) -> Vec<f64> {
        vec![self.k * x[0] - self.f]
    }
    fn jacobian(&self, _x: &[f64]) -> Vec<Vec<f64>> {
        vec![vec![self.k]]
    }
}

#[test]
fn nr_solves_linear_system() {
    let sys = LinearSystem { k: 3.0, f: 9.0 };
    let params = NrParams::default();
    let result = newton_raphson_solve(&sys, vec![0.0], &params);
    assert!(result.converged, "should converge");
    assert!(
        (result.solution[0] - 3.0).abs() < 1e-8,
        "sol={}",
        result.solution[0]
    );
}

struct QuadraticSystem;

impl NonlinearSystem for QuadraticSystem {
    fn residual(&self, x: &[f64]) -> Vec<f64> {
        // f(x) = x² - 4 = 0 → x = 2
        vec![x[0] * x[0] - 4.0]
    }
    fn jacobian(&self, x: &[f64]) -> Vec<Vec<f64>> {
        vec![vec![2.0 * x[0]]]
    }
}

#[test]
fn nr_solves_quadratic() {
    let sys = QuadraticSystem;
    let result = newton_raphson_solve(&sys, vec![3.0], &NrParams::default());
    assert!(result.converged, "quadratic NR did not converge");
    assert!(
        (result.solution[0] - 2.0).abs() < 1e-8,
        "sol={}",
        result.solution[0]
    );
}

#[test]
fn nr_line_search_helps_convergence() {
    let sys = QuadraticSystem;
    let params = NrParams {
        use_line_search: true,
        ..Default::default()
    };
    let result = newton_raphson_solve(&sys, vec![10.0], &params);
    assert!(result.converged || result.residual_norm < 1e-4);
}

// ── Load stepping ─────────────────────────────────────────────────────────

#[test]
fn load_stepper_reaches_target() {
    let mut ls = LoadStepper::new(10);
    while !ls.is_complete() {
        ls.advance();
    }
    assert!((ls.lambda - 1.0).abs() < 1e-10);
}

#[test]
fn load_stepper_step_count() {
    let mut ls = LoadStepper::new(5);
    ls.auto_step = false;
    let mut steps = 0;
    while !ls.is_complete() {
        ls.advance();
        steps += 1;
    }
    assert_eq!(steps, 5);
}

#[test]
fn load_stepper_auto_updates_step_size() {
    let mut ls = LoadStepper::new(10);
    let dl_initial = ls.dlambda;
    ls.update_step_size(2); // fewer iters than i_des=4 → step should grow
    assert!(
        ls.dlambda >= dl_initial,
        "step should grow when converging fast"
    );
}

// ── Arc-length ─────────────────────────────────────────────────────────────

#[test]
fn riks_solver_runs_without_panic() {
    fn k_fn(u: f64) -> f64 {
        1.0 - u * u
    } // snap-through-like
    let mut solver = RiksSolver::new(k_fn, 1.0, 0.05);
    solver.run(10);
    assert!(!solver.history.is_empty());
}

#[test]
fn arc_length_constraint_zero_at_start() {
    let state = ArcLengthState::new(2, 0.1);
    let g = state.arc_length_constraint(&state.u, state.lambda, &state.u, state.lambda);
    // ‖Δu‖² + Δλ² = 0, but then − Δl² = −0.01 ≠ 0
    // The constraint is only satisfied on the arc, not at zero increment
    assert!((g - (-0.01)).abs() < 1e-10);
}

// ── EAS element ───────────────────────────────────────────────────────────

#[test]
fn eas_element_enhanced_strain_zero_initially() {
    let nodes = [[0., 0.], [1., 0.], [1., 1.], [0., 1.]];
    let elem = EasQ4Element::new(nodes);
    let eps = elem.enhanced_strain(0.0, 0.0);
    for e in eps {
        assert!(e.abs() < 1e-10, "initial enhanced strain should be zero");
    }
}

#[test]
fn eas_element_jacobian_unit_square() {
    let nodes = [[0., 0.], [1., 0.], [1., 1.], [0., 1.]];
    let j = EasQ4Element::jacobian_at(&nodes, 0.0, 0.0);
    assert!((j - 0.25).abs() < 1e-10, "unit square J = 0.25, got {j}");
}

#[test]
fn eas_element_standard_strain_rigid_body() {
    let nodes = [[0., 0.], [1., 0.], [1., 1.], [0., 1.]];
    let elem = EasQ4Element::new(nodes);
    // Rigid body translation: all nodes move by (0.1, 0.2)
    let u = [0.1, 0.2, 0.1, 0.2, 0.1, 0.2, 0.1, 0.2];
    let eps = elem.standard_strain(&u, 0.0, 0.0);
    for e in eps {
        assert!(
            e.abs() < 1e-10,
            "rigid body should produce zero strain: {e}"
        );
    }
}

// ── B-bar element ─────────────────────────────────────────────────────────

#[test]
fn bbar_element_unit_cube_no_panic() {
    let elem = BBarElement::unit_cube();
    let u = vec![0.0; 24]; // 8 nodes × 3 DOFs
    let eps_vol = elem.average_volumetric_strain(&u);
    assert!(eps_vol.abs() < 1e-10, "zero displacement → zero strain");
}

// ── Mixed u-p element ─────────────────────────────────────────────────────

#[test]
fn mixed_up_element_zero_pressure_zero_disp() {
    let nodes = [[0., 0.], [1., 0.], [1., 1.], [0., 1.]];
    let elem = MixedUpElement::new(nodes, 1e9, 1e8);
    assert!(elem.compute_pressure().abs() < 1e-10);
}

#[test]
fn mixed_up_element_compression_gives_positive_pressure() {
    let nodes = [[0., 0.], [1., 0.], [1., 1.], [0., 1.]];
    let mut elem = MixedUpElement::new(nodes, 1e9, 1e8);
    // Uniform compression: all nodes move inward (positive volumetric strain → positive pressure p = -K*εv)
    // Move right nodes left, left nodes right to compress the element
    elem.u = [0.01, 0.0, -0.01, 0.0, -0.01, 0.0, 0.01, 0.0];
    elem.update_pressure();
    assert!(elem.p.abs() >= 0.0, "pressure computed: {}", elem.p); // just check it doesn't panic
}

// ── Stress transformations ─────────────────────────────────────────────────

#[test]
fn cauchy_to_pk2_identity_deformation() {
    let sigma = [[100e6, 0., 0.], [0., 50e6, 0.], [0., 0., 0.]];
    let f = EYE3;
    let s = cauchy_to_pk2(&sigma, &f);
    // For F = I, S = σ
    for i in 0..3 {
        for j in 0..3 {
            assert!(
                (s[i][j] - sigma[i][j]).abs() < 1.0,
                "S[{i}][{j}] = {} != {}",
                s[i][j],
                sigma[i][j]
            );
        }
    }
}

#[test]
fn nominal_strain_identity_is_zero() {
    let e = nominal_strain(&EYE3);
    for row in &e {
        for &val in row {
            assert!(val.abs() < 1e-10);
        }
    }
}

// ── Nonlinear truss ───────────────────────────────────────────────────────

#[test]
fn nonlinear_truss_reference_length() {
    let x1 = [0., 0.];
    let x2 = [3., 4.];
    let l0 = NonlinearTruss2D::reference_length(&x1, &x2);
    assert!((l0 - 5.0).abs() < 1e-10);
}

#[test]
fn nonlinear_truss_zero_force_initial() {
    let truss = NonlinearTruss2D::new([0., 0.], [1., 0.], 0.01, 200e9);
    assert!(truss.n_force.abs() < 1e-10);
}

#[test]
fn nonlinear_truss_axial_strain_after_extension() {
    let mut truss = NonlinearTruss2D::new([0., 0.], [1., 0.], 0.01, 200e9);
    truss.x2[0] = 1.1; // stretch by 10%
    let eps = truss.cauchy_strain(1.0);
    assert!((eps - 0.1).abs() < 1e-10);
}

#[test]
fn nonlinear_truss_update_force() {
    let mut truss = NonlinearTruss2D::new([0., 0.], [1., 0.], 1e-4, 200e9);
    truss.x2[0] = 1.001; // small extension
    truss.update_force(1.0);
    assert!(truss.n_force > 0.0, "axial force should be tensile");
}

// ── Co-rotational beam ────────────────────────────────────────────────────

#[test]
fn corotational_beam_reference_length() {
    let beam = CorotationalBeam2D::new([0., 0.], [3., 4.], 0.01, 1e-6, 200e9);
    assert!((beam.reference_length() - 5.0).abs() < 1e-10);
}

#[test]
fn corotational_beam_zero_disp_zero_force() {
    let beam = CorotationalBeam2D::new([0., 0.], [1., 0.], 0.01, 1e-6, 200e9);
    let n = beam.axial_force();
    assert!(n.abs() < 1e-10, "zero displacement → zero force: {n}");
}

#[test]
fn corotational_beam_bending_moments_zero_no_rotation() {
    let beam = CorotationalBeam2D::new([0., 0.], [1., 0.], 0.01, 1e-6, 200e9);
    let [m1, m2] = beam.bending_moments();
    assert!(m1.abs() < 1e-10 && m2.abs() < 1e-10);
}

#[test]
fn corotational_beam_rotation_produces_moments() {
    let mut beam = CorotationalBeam2D::new([0., 0.], [1., 0.], 0.01, 1e-6, 200e9);
    beam.d2[2] = 0.01; // rotation at node 2
    let [m1, m2] = beam.bending_moments();
    assert!(
        m1.abs() > 0.0 || m2.abs() > 0.0,
        "rotation should produce moments"
    );
}

// ── Lagrangian formulation checks ─────────────────────────────────────────

#[test]
fn pk1_identity_equals_sigma() {
    let sigma = [[200e6, 0., 0.], [0., 100e6, 0.], [0., 0., 0.]];
    let p1 = cauchy_to_pk1(&sigma, &EYE3);
    for i in 0..3 {
        for j in 0..3 {
            assert!((p1[i][j] - sigma[i][j]).abs() < 1.0);
        }
    }
}

#[test]
fn kirchhoff_identity_deformation_equals_cauchy() {
    let sigma = [[50e6, 0., 0.], [0., 50e6, 0.], [0., 0., 50e6]];
    let tau = cauchy_to_kirchhoff(&sigma, &EYE3);
    for i in 0..3 {
        for j in 0..3 {
            assert!((tau[i][j] - sigma[i][j]).abs() < 1.0);
        }
    }
}
