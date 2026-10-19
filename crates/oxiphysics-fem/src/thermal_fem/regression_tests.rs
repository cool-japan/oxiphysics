//! Regression tests for the tetrahedral thermal FEM.
//!
//! Covers the exact shape-function gradients of skewed tetrahedra (linear patch
//! test), the full-matrix θ-method, the Newton-linearised radiation boundary and
//! the equivalent nodal thermo-elastic load.

use super::*;
use crate::analysis::FemError;
use std::f64::consts::PI;

// ── helpers ──────────────────────────────────────────────────────────────────

/// Tetrahedron with a non-symmetric Jacobian `J = [p1 | p2 | p3]`.
fn skewed_tet() -> [[f64; 3]; 4] {
    [
        [0.0, 0.0, 0.0],
        [1.0, 0.2, 0.1],
        [0.3, 1.0, 0.4],
        [0.2, 0.5, 1.5],
    ]
}

fn unit_tet() -> [[f64; 3]; 4] {
    [
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
    ]
}

/// `T(x) = g · x + c`.
fn linear_field(g: [f64; 3], c: f64, p: [f64; 3]) -> f64 {
    g[0] * p[0] + g[1] * p[1] + g[2] * p[2] + c
}

/// Tetrahedron volume from the scalar triple product, computed independently of
/// the production code.
fn tet_volume(nodes: &[[f64; 3]; 4]) -> f64 {
    let e = |i: usize| {
        [
            nodes[i][0] - nodes[0][0],
            nodes[i][1] - nodes[0][1],
            nodes[i][2] - nodes[0][2],
        ]
    };
    let (a, b, c) = (e(1), e(2), e(3));
    let triple = a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0])
        + a[2] * (b[0] * c[1] - b[1] * c[0]);
    triple.abs() / 6.0
}

/// 1-D rod conduction matrix `K` and lumped capacity matrix `M` for `n` nodes:
/// conductance `g` and capacity `c_elem` per element.
fn rod_matrices(n: usize, g: f64, c_elem: f64) -> (Vec<Vec<f64>>, Vec<Vec<f64>>) {
    let mut k = vec![vec![0.0_f64; n]; n];
    let mut m = vec![vec![0.0_f64; n]; n];
    for e in 0..n - 1 {
        k[e][e] += g;
        k[e][e + 1] -= g;
        k[e + 1][e] -= g;
        k[e + 1][e + 1] += g;
        m[e][e] += 0.5 * c_elem;
        m[e + 1][e + 1] += 0.5 * c_elem;
    }
    (k, m)
}

/// Direct dense solution of the θ-method exactly as documented:
/// `(C + θ dt K) T₁ = (C − (1−θ) dt K) T₀ + dt f`.
fn dense_theta_reference(
    theta: f64,
    dt: f64,
    temps: &[f64],
    k: &[Vec<f64>],
    m: &[Vec<f64>],
    f: &[f64],
) -> Vec<f64> {
    let n = temps.len();
    let c: Vec<f64> = m.iter().map(|row| row.iter().sum::<f64>()).collect();
    let mut a = vec![vec![0.0_f64; n]; n];
    let mut rhs = vec![0.0_f64; n];
    for i in 0..n {
        for j in 0..n {
            a[i][j] = theta * dt * k[i][j] + if i == j { c[i] } else { 0.0 };
        }
        let kt: f64 = (0..n).map(|j| k[i][j] * temps[j]).sum();
        rhs[i] = c[i] * temps[i] - (1.0 - theta) * dt * kt + dt * f[i];
    }
    solve_thermal_steady(&a, &rhs)
}

/// Four-tetrahedron patch: outer tetrahedron `A B C D` with one interior node `P`
/// connected to all four faces (orientations deliberately mixed).
fn patch_mesh(boundary_temp: impl Fn([f64; 3]) -> f64) -> (Vec<ThermalNode>, Vec<ThermalElement>) {
    let a = [0.0, 0.0, 0.0];
    let b = [2.0, 0.3, 0.1];
    let c = [0.4, 1.8, 0.2];
    let d = [0.3, 0.5, 2.2];
    let w = [0.4, 0.3, 0.2, 0.1];
    let mut p = [0.0_f64; 3];
    for i in 0..3 {
        p[i] = w[0] * a[i] + w[1] * b[i] + w[2] * c[i] + w[3] * d[i];
    }
    let nodes = vec![
        ThermalNode::dirichlet(a, boundary_temp(a)),
        ThermalNode::dirichlet(b, boundary_temp(b)),
        ThermalNode::dirichlet(c, boundary_temp(c)),
        ThermalNode::dirichlet(d, boundary_temp(d)),
        ThermalNode::free(p, 0.0),
    ];
    let kc = [1.0, 2.0, 3.0];
    let elements = vec![
        ThermalElement::new([4, 1, 2, 3], kc),
        ThermalElement::new([0, 4, 2, 3], kc),
        ThermalElement::new([0, 1, 4, 3], kc),
        ThermalElement::new([0, 1, 2, 4], kc),
    ];
    (nodes, elements)
}

// ── F1378: exact tetrahedral shape gradients ─────────────────────────────────

#[test]
fn skewed_tet_gradients_form_a_partition_of_unity_and_reproduce_coordinates() {
    let nodes = skewed_tet();
    let (grads, vol) = tet_shape_gradients(&nodes).expect("non-degenerate tetrahedron");
    assert!((vol - tet_volume(&nodes)).abs() < 1e-14);
    for d in 0..3 {
        let sum: f64 = grads.iter().map(|g| g[d]).sum();
        assert!(sum.abs() < 1e-13, "sum of ∂N/∂x_{d} = {sum}");
    }
    // x_i = Σ_a N_a x_{a,i}  ⇒  Σ_a x_{a,i} ∂N_a/∂x_j = δ_ij
    for i in 0..3 {
        for j in 0..3 {
            let s: f64 = (0..4).map(|a| nodes[a][i] * grads[a][j]).sum();
            let expected = if i == j { 1.0 } else { 0.0 };
            assert!(
                (s - expected).abs() < 1e-13,
                "Σ x_{i} ∂N/∂x_{j} = {s}, expected {expected}"
            );
        }
    }
}

#[test]
fn skewed_tet_heat_flux_recovers_linear_field_exactly() {
    let nodes = skewed_tet();
    // T = x: the finding's counter-example (the old gradient gave [0.966, 0.084, 0.047]).
    let temps: Vec<f64> = nodes.iter().map(|p| p[0]).collect();
    let q = heat_flux_element(&temps, &nodes, [1.0, 1.0, 1.0]);
    assert!((q[0] + 1.0).abs() < 1e-12, "q = {q:?}");
    assert!(q[1].abs() < 1e-12, "q = {q:?}");
    assert!(q[2].abs() < 1e-12, "q = {q:?}");

    // General linear field and anisotropic conductivity: q_i = -k_i g_i.
    let g = [2.0, -1.0, 3.0];
    let kc = [2.0, 3.0, 5.0];
    let temps: Vec<f64> = nodes.iter().map(|&p| linear_field(g, 7.0, p)).collect();
    let q = heat_flux_element(&temps, &nodes, kc);
    for i in 0..3 {
        assert!(
            (q[i] + kc[i] * g[i]).abs() < 1e-11,
            "component {i}: {} vs {}",
            q[i],
            -kc[i] * g[i]
        );
    }
}

#[test]
fn skewed_tet_conductivity_matrix_is_exact_for_linear_fields() {
    let nodes = skewed_tet();
    let kc = [2.0, 3.0, 5.0];
    let ke = tet_thermal_conductivity_matrix(&nodes, kc);
    let scale = ke
        .iter()
        .flat_map(|r| r.iter())
        .fold(0.0_f64, |m, v| m.max(v.abs()));

    for (i, row) in ke.iter().enumerate() {
        let row_sum: f64 = row.iter().sum();
        assert!(row_sum.abs() < 1e-12 * scale, "row {i} sum = {row_sum}");
        for (j, &v) in row.iter().enumerate() {
            assert!((v - ke[j][i]).abs() < 1e-12 * scale, "asymmetry ({i},{j})");
        }
    }

    // T = g·x + c: Tᵀ K_e T = V (kx gx² + ky gy² + kz gz²) exactly.
    let g = [1.5, -0.5, 2.0];
    let t: Vec<f64> = nodes.iter().map(|&p| linear_field(g, 4.0, p)).collect();
    let mut energy = 0.0;
    for i in 0..4 {
        for j in 0..4 {
            energy += t[i] * ke[i][j] * t[j];
        }
    }
    let expected =
        tet_volume(&nodes) * (kc[0] * g[0] * g[0] + kc[1] * g[1] * g[1] + kc[2] * g[2] * g[2]);
    assert!(
        (energy - expected).abs() < 1e-10 * expected,
        "energy {energy} vs {expected}"
    );

    // K_e T = V B^T D g: nodal heat fluxes of a linear field balance pairwise.
    let net_flux: f64 = (0..4)
        .map(|i| (0..4).map(|j| ke[i][j] * t[j]).sum::<f64>())
        .sum();
    assert!(net_flux.abs() < 1e-11 * scale * 10.0, "net flux {net_flux}");
}

#[test]
fn degenerate_tetrahedra_give_zero_results() {
    // Coplanar vertices (z = 0) and coincident vertices.
    let flat = [
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [1.0, 1.0, 0.0],
    ];
    let coincident = [[0.5_f64; 3]; 4];
    for nodes in [flat, coincident] {
        assert!(tet_shape_gradients(&nodes).is_none());
        let ke = tet_thermal_conductivity_matrix(&nodes, [1.0, 2.0, 3.0]);
        assert!(ke.iter().flatten().all(|&v| v == 0.0));
        let q = heat_flux_element(&[1.0, 2.0, 3.0, 4.0], &nodes, [1.0, 1.0, 1.0]);
        assert_eq!(q, [0.0; 3]);
        let te = ThermoElasticCoupling::new(1e-5, 2e11, 0.3);
        let load = te.nodal_thermal_load(&nodes, 25.0);
        assert_eq!(load.len(), 12);
        assert!(load.iter().all(|&v| v == 0.0));
    }
    // A non-finite coordinate is rejected as well.
    let mut bad = skewed_tet();
    bad[2][1] = f64::NAN;
    assert!(tet_shape_gradients(&bad).is_none());
}

#[test]
fn tiny_but_well_shaped_tetrahedron_is_not_rejected() {
    // Nano-scale element (edge 1e-9 m): the degeneracy test is scale invariant.
    let s = 1.0e-9;
    let nodes = [
        [0.0, 0.0, 0.0],
        [s, 0.2 * s, 0.1 * s],
        [0.3 * s, s, 0.4 * s],
        [0.2 * s, 0.5 * s, 1.5 * s],
    ];
    let (grads, vol) = tet_shape_gradients(&nodes).expect("valid tetrahedron");
    assert!(vol > 0.0 && vol < 1e-26);
    // ∇N scales as 1/s and still sums to zero relative to that scale.
    for d in 0..3 {
        let sum: f64 = grads.iter().map(|g| g[d]).sum();
        assert!(sum.abs() * s < 1e-12, "Σ ∂N/∂x_{d} = {sum}");
    }
    let temps: Vec<f64> = nodes.iter().map(|p| p[0] / s).collect();
    let q = heat_flux_element(&temps, &nodes, [1.0, 1.0, 1.0]);
    // T = x / s  ⇒  ∇T = [1/s, 0, 0]
    assert!((q[0] * s + 1.0).abs() < 1e-9, "q = {q:?}");
    assert!(q[1].abs() * s < 1e-9 && q[2].abs() * s < 1e-9, "q = {q:?}");
}

#[test]
fn linear_patch_test_on_irregular_assembled_mesh() {
    let g = [3.0, -2.0, 0.5];
    let c0 = 10.0;
    let field = move |p: [f64; 3]| linear_field(g, c0, p);
    let (nodes, elements) = patch_mesh(field);
    let mut k = assemble_thermal_stiffness(&elements, &nodes);
    let mut f = vec![0.0_f64; nodes.len()];
    apply_thermal_dirichlet(&mut k, &mut f, &nodes);
    let t = solve_thermal_steady(&k, &f);
    let expected = field(nodes[4].position);
    assert!(
        (t[4] - expected).abs() < 1e-9,
        "interior node: {} vs {expected}",
        t[4]
    );
    // Every sub-element sees the same constant gradient.
    for elem in &elements {
        let pos: [[f64; 3]; 4] = [
            nodes[elem.nodes[0]].position,
            nodes[elem.nodes[1]].position,
            nodes[elem.nodes[2]].position,
            nodes[elem.nodes[3]].position,
        ];
        let temps = [
            t[elem.nodes[0]],
            t[elem.nodes[1]],
            t[elem.nodes[2]],
            t[elem.nodes[3]],
        ];
        let q = heat_flux_element(&temps, &pos, elem.conductivity);
        for i in 0..3 {
            assert!(
                (q[i] + elem.conductivity[i] * g[i]).abs() < 1e-9,
                "element flux component {i}: {}",
                q[i]
            );
        }
    }
}

// ── F1380: θ-method solves the full system ───────────────────────────────────

#[test]
fn implicit_theta_steps_preserve_a_uniform_field_on_a_conducting_rod() {
    let (k, m) = rod_matrices(6, 3.0, 0.4);
    let t0 = vec![300.0_f64; 6];
    let f = vec![0.0_f64; 6];
    for solver in [
        TransientThermalSolver::implicit_euler(0.1),
        TransientThermalSolver::crank_nicolson(0.1),
        TransientThermalSolver {
            theta: 0.7,
            dt: 0.25,
        },
    ] {
        let t1 = solver.step(&t0, &k, &m, &f).expect("step");
        for (i, &t) in t1.iter().enumerate() {
            assert!(
                (t - 300.0).abs() < 1e-10,
                "theta = {}: node {i} drifted to {t}",
                solver.theta
            );
        }
    }
}

#[test]
fn implicit_theta_steps_preserve_a_uniform_field_on_a_tet_mesh() {
    let (nodes, elements) = patch_mesh(|_| 0.0);
    let k = assemble_thermal_stiffness(&elements, &nodes);
    let params = ThermalFemParams::isotropic(1.0, 1.0, 1.0);
    let m = assemble_lumped_capacity(&elements, &nodes, &params);
    let t0 = vec![300.0_f64; nodes.len()];
    let f = vec![0.0_f64; nodes.len()];
    for theta in [0.5_f64, 1.0] {
        let solver = TransientThermalSolver { theta, dt: 0.05 };
        let history = solver.integrate(&t0, &k, &m, &f, 5).expect("integrate");
        for (s, state) in history.iter().enumerate() {
            for (i, &t) in state.iter().enumerate() {
                assert!(
                    (t - 300.0).abs() < 1e-10,
                    "theta = {theta}, step {s}, node {i}: {t}"
                );
            }
        }
    }
}

#[test]
fn crank_nicolson_matches_analytic_decay_on_a_fifty_node_rod() {
    // ∂T/∂t = α ∂²T/∂x², T(0) = T(L) = 0, T(x, 0) = sin(π x / L)
    //   ⇒ T(x, t) = sin(π x / L) exp(−α π² t / L²).
    let n = 50;
    let length = 1.0_f64;
    let h = length / (n - 1) as f64;
    let (conductivity, rho_cp) = (2.0_f64, 4.0_f64);
    let alpha = conductivity / rho_cp;
    let (mut k, m) = rod_matrices(n, conductivity / h, rho_cp * h);
    for d in [0, n - 1] {
        for j in 0..n {
            k[d][j] = 0.0;
            k[j][d] = 0.0;
        }
        k[d][d] = 1.0;
    }
    let mut t0: Vec<f64> = (0..n).map(|i| (PI * i as f64 * h / length).sin()).collect();
    t0[n - 1] = 0.0;
    let f = vec![0.0_f64; n];

    let dt = 0.002;
    let n_steps = 100;
    let solver = TransientThermalSolver::crank_nicolson(dt);
    let history = solver
        .integrate(&t0, &k, &m, &f, n_steps)
        .expect("Crank-Nicolson integration");
    assert_eq!(history.len(), n_steps + 1);
    let t_end = dt * n_steps as f64;
    let decay = (-alpha * PI * PI * t_end / (length * length)).exp();
    let last = &history[n_steps];
    let mut max_err = 0.0_f64;
    for (i, &t) in last.iter().enumerate() {
        let exact = (PI * i as f64 * h / length).sin() * decay;
        max_err = max_err.max((t - exact).abs());
    }
    assert!(
        max_err < 0.01 * decay,
        "max error {max_err} exceeds 1% of the amplitude {decay}"
    );
    // The decay is real (a diagonal-only solve damps the mode far too fast).
    let mid = n / 2;
    assert!(last[mid] > 0.5 * decay && last[mid] < 1.1 * decay);
}

#[test]
fn implicit_euler_converges_to_the_steady_state_of_a_driven_rod() {
    let n = 21;
    let h = 1.0 / (n - 1) as f64;
    let (mut k, m) = rod_matrices(n, 2.0 / h, 4.0 * h);
    let mut nodes: Vec<ThermalNode> = (0..n)
        .map(|i| ThermalNode::free([i as f64 * h, 0.0, 0.0], 0.0))
        .collect();
    nodes[0] = ThermalNode::dirichlet([0.0; 3], 0.0);
    nodes[n - 1] = ThermalNode::dirichlet([1.0, 0.0, 0.0], 100.0);
    let mut f = vec![0.0_f64; n];
    apply_thermal_dirichlet(&mut k, &mut f, &nodes);
    let mut t = vec![0.0_f64; n];
    t[n - 1] = 100.0;
    let solver = TransientThermalSolver::implicit_euler(0.5);
    for _ in 0..400 {
        t = solver.step(&t, &k, &m, &f).expect("step");
    }
    for (i, &ti) in t.iter().enumerate() {
        let expected = 100.0 * i as f64 / (n - 1) as f64;
        assert!((ti - expected).abs() < 1e-7, "node {i}: {ti} vs {expected}");
    }
}

#[test]
fn theta_step_equals_the_dense_direct_solution_of_the_documented_system() {
    let (nodes, elements) = patch_mesh(|_| 0.0);
    let k = assemble_thermal_stiffness(&elements, &nodes);
    let params = ThermalFemParams::isotropic(1.0, 1.0, 1.0);
    let m = assemble_lumped_capacity(&elements, &nodes, &params);
    let temps = vec![10.0, 55.0, -4.0, 120.0, 37.5];
    let f = vec![0.5, -1.0, 2.0, 0.0, 3.0];
    for (theta, dt) in [(0.5_f64, 0.05_f64), (1.0, 0.2), (0.3, 0.01), (0.8, 5.0)] {
        let solver = TransientThermalSolver { theta, dt };
        let got = solver.step(&temps, &k, &m, &f).expect("step");
        let want = dense_theta_reference(theta, dt, &temps, &k, &m, &f);
        for i in 0..temps.len() {
            assert!(
                (got[i] - want[i]).abs() < 1e-8 * (1.0 + want[i].abs()),
                "theta {theta}, dt {dt}, node {i}: {} vs {}",
                got[i],
                want[i]
            );
        }
    }
}

#[test]
fn theta_step_handles_a_non_symmetric_conduction_matrix() {
    // Advection-like skew part: the implicit system is solved by BiCGSTAB.
    let (mut k, m) = rod_matrices(8, 2.0, 0.5);
    for i in 0..7 {
        k[i][i + 1] += 0.6;
        k[i + 1][i] -= 0.6;
    }
    let temps: Vec<f64> = (0..8).map(|i| 20.0 + 5.0 * (i as f64).sin()).collect();
    let f = vec![0.25_f64; 8];
    for theta in [0.5_f64, 1.0] {
        let dt = 0.1;
        let solver = TransientThermalSolver { theta, dt };
        let got = solver.step(&temps, &k, &m, &f).expect("step");
        let want = dense_theta_reference(theta, dt, &temps, &k, &m, &f);
        for i in 0..8 {
            assert!(
                (got[i] - want[i]).abs() < 1e-8 * (1.0 + want[i].abs()),
                "theta {theta}, node {i}: {} vs {}",
                got[i],
                want[i]
            );
        }
    }
}

#[test]
fn explicit_theta_step_is_forward_euler_and_integrate_matches_repeated_steps() {
    let (k, m) = rod_matrices(5, 1.5, 0.3);
    let t0 = vec![4.0, 8.0, 15.0, 16.0, 23.0];
    let f = vec![0.0, 0.5, 0.0, 0.0, 1.0];
    let dt = 0.01;
    let explicit = TransientThermalSolver::explicit_euler(dt);
    let t1 = explicit.step(&t0, &k, &m, &f).expect("explicit step");
    let c: Vec<f64> = m.iter().map(|row| row.iter().sum::<f64>()).collect();
    for i in 0..5 {
        let kt: f64 = (0..5).map(|j| k[i][j] * t0[j]).sum();
        let expected = t0[i] + dt * (f[i] - kt) / c[i];
        assert!((t1[i] - expected).abs() < 1e-12, "node {i}");
    }

    let cn = TransientThermalSolver::crank_nicolson(0.05);
    let history = cn.integrate(&t0, &k, &m, &f, 3).expect("integrate");
    assert_eq!(history.len(), 4);
    assert_eq!(history[0], t0);
    let mut t = t0.clone();
    for state in history.iter().skip(1) {
        t = cn.step(&t, &k, &m, &f).expect("step");
        assert_eq!(*state, t);
    }
}

#[test]
fn transient_solver_reports_invalid_input_instead_of_panicking() {
    let eye = vec![vec![1.0_f64]];
    let t = [1.0_f64];
    let f = [0.0_f64];
    let invalid = |r: Result<Vec<f64>, FemError>| matches!(r, Err(FemError::InvalidInput(_)));

    assert!(invalid(
        TransientThermalSolver::implicit_euler(0.0).step(&t, &eye, &eye, &f)
    ));
    assert!(invalid(
        TransientThermalSolver::implicit_euler(-1.0).step(&t, &eye, &eye, &f)
    ));
    assert!(invalid(
        TransientThermalSolver::implicit_euler(f64::NAN).step(&t, &eye, &eye, &f)
    ));
    assert!(invalid(
        TransientThermalSolver {
            theta: 1.5,
            dt: 0.1
        }
        .step(&t, &eye, &eye, &f)
    ));
    // size mismatches
    assert!(invalid(TransientThermalSolver::implicit_euler(0.1).step(
        &[1.0, 2.0],
        &eye,
        &eye,
        &[0.0, 0.0]
    )));
    assert!(invalid(TransientThermalSolver::implicit_euler(0.1).step(
        &t,
        &eye,
        &eye,
        &[0.0, 0.0]
    )));
    // non-finite entries
    let nan_k = vec![vec![f64::NAN]];
    assert!(invalid(
        TransientThermalSolver::implicit_euler(0.1).step(&t, &nan_k, &eye, &f)
    ));
    assert!(invalid(TransientThermalSolver::implicit_euler(0.1).step(
        &[f64::INFINITY],
        &eye,
        &eye,
        &f
    )));
    // C + θ dt K with a non-positive diagonal cannot be solved by CG.
    let negative = vec![vec![-10.0_f64]];
    assert_eq!(
        TransientThermalSolver::implicit_euler(1.0).step(&t, &negative, &eye, &f),
        Err(FemError::SingularSystem)
    );
}

// ── F1381: Newton-linearised radiation boundary ──────────────────────────────

#[test]
fn radiation_newton_iteration_converges_monotonically_to_the_surroundings() {
    let rb = RadiationBoundary::new(1.0, 300.0);
    let mut t = 1000.0_f64;
    for iter in 0..80 {
        let mut k = vec![vec![0.0_f64]];
        let mut f = vec![0.0_f64];
        rb.apply_to_node(&mut k, &mut f, 0, t);
        assert!(k[0][0] > 0.0);
        let t_new = f[0] / k[0][0];
        assert!(
            t_new <= t + 1e-12 && t_new >= 300.0 - 1e-9,
            "iteration {iter}: {t} -> {t_new} is not a monotone approach to 300 K"
        );
        t = t_new;
    }
    assert!((t - 300.0).abs() < 1e-8, "converged to {t}");
}

#[test]
fn radiation_first_solve_at_equilibrium_returns_the_surrounding_temperature() {
    for (emissivity, t_surr) in [(1.0_f64, 300.0_f64), (0.7, 450.0), (0.2, 1200.0)] {
        let rb = RadiationBoundary::new(emissivity, t_surr);
        let mut k = vec![vec![0.0_f64]];
        let mut f = vec![0.0_f64];
        rb.apply_to_node(&mut k, &mut f, 0, t_surr);
        let t = f[0] / k[0][0];
        assert!(
            (t - t_surr).abs() < 1e-12 * t_surr,
            "T0 = T_surr = {t_surr}: first solve gave {t}"
        );
    }
}

#[test]
fn radiation_rhs_is_the_newton_linearisation_h_t_star_minus_q() {
    // T* = 1000 K, T_surr = 300 K, ε = 1: h = 4σT*³, f = h T* − q(T*) ≈ +1.7057e5 W/m²,
    // not the previous (negative) h T_surr + q − h T* ≈ −1.0253e5.
    let rb = RadiationBoundary::new(1.0, 300.0);
    let mut k = vec![vec![0.0_f64]];
    let mut f = vec![0.0_f64];
    rb.apply_to_node(&mut k, &mut f, 0, 1000.0);
    let h = 4.0 * STEFAN_BOLTZMANN * 1000.0_f64.powi(3);
    let q = STEFAN_BOLTZMANN * (1000.0_f64.powi(4) - 300.0_f64.powi(4));
    assert!((k[0][0] - h).abs() < 1e-9);
    assert!((f[0] - (h * 1000.0 - q)).abs() < 1e-6);
    assert!(f[0] > 1.70e5 && f[0] < 1.71e5, "f = {}", f[0]);
    // The linear model reproduces the flux at the linearisation point.
    assert!((k[0][0] * 1000.0 - f[0] - q).abs() < 1e-6);
}

#[test]
fn radiation_with_conduction_satisfies_the_energy_balance_at_convergence() {
    // Node fed by conduction from a 1000 K wall (g = 50 W/(m²K)) and radiating to 300 K.
    let g = 50.0_f64;
    let t_wall = 1000.0_f64;
    let rb = RadiationBoundary::new(0.8, 300.0);
    let mut t = t_wall;
    for _ in 0..100 {
        let mut k = vec![vec![g]];
        let mut f = vec![g * t_wall];
        rb.apply_to_node(&mut k, &mut f, 0, t);
        t = f[0] / k[0][0];
    }
    let residual = g * (t - t_wall) + rb.heat_flux(t);
    assert!(residual.abs() < 1e-6, "energy balance residual {residual}");
    assert!(t > 300.0 && t < t_wall, "T = {t}");
}

#[test]
fn radiation_apply_to_node_ignores_out_of_range_nodes() {
    let rb = RadiationBoundary::new(0.9, 300.0);
    let mut k = vec![vec![1.0_f64]];
    let mut f = vec![2.0_f64];
    rb.apply_to_node(&mut k, &mut f, 5, 800.0);
    assert_eq!(k[0][0], 1.0);
    assert_eq!(f[0], 2.0);
}

// ── F1382: equivalent nodal thermo-elastic load ──────────────────────────────

#[test]
fn thermal_load_of_a_skewed_tet_has_zero_resultant() {
    let te = ThermoElasticCoupling::new(1.2e-5, 2.0e11, 0.3);
    let load = te.nodal_thermal_load(&skewed_tet(), 80.0);
    assert_eq!(load.len(), 12);
    let max_abs = load.iter().fold(0.0_f64, |m, v| m.max(v.abs()));
    assert!(max_abs > 0.0);
    for dir in 0..3 {
        let sum: f64 = (0..4).map(|a| load[3 * a + dir]).sum();
        assert!(
            sum.abs() < 1e-12 * max_abs,
            "resultant in direction {dir} = {sum} (max |f| = {max_abs})"
        );
    }
}

#[test]
fn thermal_load_of_the_unit_tet_matches_the_hand_computed_values() {
    let (alpha, young, nu, d_t) = (1.0e-5, 2.0e11, 0.3, 50.0);
    let te = ThermoElasticCoupling::new(alpha, young, nu);
    let load = te.nodal_thermal_load(&unit_tet(), d_t);
    // σ_th = −E α ΔT / (1 − 2ν) (hydrostatic); f_a = −V σ_th ∇N_a, V = 1/6.
    let sigma = -young * alpha * d_t / (1.0 - 2.0 * nu);
    let v = 1.0 / 6.0;
    let outward = -v * sigma; // magnitude of the force on a node pushed outward
    assert!(outward > 0.0, "heating pushes the nodes outward");
    let expected = [
        [-outward, -outward, -outward], // node 0: ∇N0 = (−1, −1, −1)
        [outward, 0.0, 0.0],            // node 1
        [0.0, outward, 0.0],            // node 2
        [0.0, 0.0, outward],            // node 3
    ];
    for a in 0..4 {
        for i in 0..3 {
            let got = load[3 * a + i];
            assert!(
                (got - expected[a][i]).abs() < 1e-9 * outward,
                "node {a} dof {i}: {got} vs {}",
                expected[a][i]
            );
        }
    }
}

#[test]
fn thermal_load_does_the_free_expansion_work() {
    // Virtual work of the nodal loads on the free thermal expansion field
    // u = α ΔT x (strain ε_th = α ΔT I) is  V ε_thᵀ D ε_th = 3 V E (α ΔT)² / (1 − 2ν).
    let (alpha, young, nu, d_t) = (1.2e-5, 7.0e10, 0.33, 120.0);
    let te = ThermoElasticCoupling::new(alpha, young, nu);
    for nodes in [unit_tet(), skewed_tet()] {
        let load = te.nodal_thermal_load(&nodes, d_t);
        let mut work = 0.0;
        for (a, p) in nodes.iter().enumerate() {
            for i in 0..3 {
                work += load[3 * a + i] * alpha * d_t * p[i];
            }
        }
        let expected = 3.0 * tet_volume(&nodes) * young * (alpha * d_t).powi(2) / (1.0 - 2.0 * nu);
        assert!(
            (work - expected).abs() < 1e-9 * expected,
            "work {work} vs {expected}"
        );
    }
}

#[test]
fn thermal_load_is_antisymmetric_in_delta_t_and_zero_without_heating() {
    let te = ThermoElasticCoupling::new(1.0e-5, 2.0e11, 0.3);
    let nodes = skewed_tet();
    let hot = te.nodal_thermal_load(&nodes, 40.0);
    let cold = te.nodal_thermal_load(&nodes, -40.0);
    let none = te.nodal_thermal_load(&nodes, 0.0);
    for i in 0..12 {
        assert!((hot[i] + cold[i]).abs() < 1e-12 * hot[i].abs().max(1.0));
        assert_eq!(none[i], 0.0);
    }
}
