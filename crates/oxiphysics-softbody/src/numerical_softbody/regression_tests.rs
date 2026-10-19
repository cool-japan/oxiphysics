//! Regression tests for the implicit spring-mass solvers and the Projective
//! Dynamics solver.
//!
//! * F4238: `NewtonRaphsonSolver::minimise` and `BackwardEulerIntegrator::step`
//!   used to be diagonal (Jacobi) gradient descent; the PCG solver, the PCG
//!   parameters and `ImplicitParticle::external_force` were never used.  They
//!   now assemble the sparse Hessian `M / dt^2 + K` and solve each Newton
//!   system with PCG.
//! * F4241: `ProjectiveDynamicsSolver::step` re-derived the (constant) global
//!   diagonal for every particle by scanning all constraints in every iteration
//!   (`O(n * C)` per iteration).  It is now assembled once per step in `O(C)`.

use super::functions::{
    ImplicitStepProblem, assemble_spring_hessian, newton_pcg_minimise, v3_add, v3_dot, v3_norm,
    v3_scale, v3_sub, vec_norm,
};
use super::types::{
    BackwardEulerIntegrator, CsrMatrix, ImplicitParticle, ImplicitSpring, LineSearch,
    NewtonRaphsonSolver, PcgParams, PcgSolver, ProjectiveConstraint, ProjectiveConstraintKind,
    ProjectiveDynamicsSolver,
};

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

fn to_dense(m: &CsrMatrix) -> Vec<Vec<f64>> {
    let mut dense = vec![vec![0.0f64; m.ncols]; m.nrows];
    for (r, dense_row) in dense.iter_mut().enumerate() {
        for k in m.row_ptr[r]..m.row_ptr[r + 1] {
            dense_row[m.col_idx[k]] += m.values[k];
        }
    }
    dense
}

/// Vertical chain hanging below a static particle at the origin.
fn hanging_chain(
    n: usize,
    rest: f64,
    stiffness: f64,
    damping: f64,
) -> (Vec<ImplicitParticle>, Vec<ImplicitSpring>) {
    let mut particles = Vec::with_capacity(n);
    particles.push(ImplicitParticle::new_static([0.0; 3]));
    for i in 1..n {
        particles.push(ImplicitParticle::new([0.0, -(i as f64) * rest, 0.0], 1.0));
    }
    let springs = (0..n - 1)
        .map(|i| ImplicitSpring::new(i, i + 1, rest, stiffness, damping))
        .collect();
    (particles, springs)
}

/// Norm of the stationarity residual of the incremental potential
/// `sum m/(2dt^2)|x - y|^2 + sum spring energy - sum f . x`, computed from the
/// public `ImplicitSpring::force_on_i` (independent of the solver's own
/// gradient).  Valid for particles at rest (the damper then contributes zero).
fn stationarity_residual(
    particles: &[ImplicitParticle],
    springs: &[ImplicitSpring],
    inertia_target: &[[f64; 3]],
    dt: f64,
) -> f64 {
    let mut residual = vec![[0.0f64; 3]; particles.len()];
    for (i, p) in particles.iter().enumerate() {
        if p.is_static() {
            continue;
        }
        let inertial = v3_scale(v3_sub(p.position, inertia_target[i]), p.mass / (dt * dt));
        residual[i] = v3_sub(inertial, p.external_force);
    }
    for s in springs {
        let f = s.force_on_i(particles);
        if !particles[s.i].is_static() {
            residual[s.i] = v3_sub(residual[s.i], f);
        }
        if !particles[s.j].is_static() {
            residual[s.j] = v3_add(residual[s.j], f);
        }
    }
    residual.iter().map(|r| v3_dot(*r, *r)).sum::<f64>().sqrt()
}

fn positions_of(particles: &[ImplicitParticle]) -> Vec<[f64; 3]> {
    particles.iter().map(|p| p.position).collect()
}

/// A small deterministic generator for jittered test geometry.
struct Lcg(u64);

impl Lcg {
    fn next_unit(&mut self) -> f64 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        ((self.0 >> 11) as f64) / ((1u64 << 53) as f64)
    }
}

/// A scattered, over-stretched tetrahedron-like cluster: one static particle,
/// three free ones, and springs that are all strictly longer than their rest
/// length (so the exact spring Hessian has no clamped transverse term).
fn stretched_cluster() -> (Vec<ImplicitParticle>, Vec<ImplicitSpring>) {
    let particles = vec![
        ImplicitParticle::new_static([0.0, 0.0, 0.0]),
        ImplicitParticle::new([1.3, 0.1, -0.2], 2.0),
        ImplicitParticle::new([0.2, 1.2, 0.4], 3.0),
        ImplicitParticle::new([0.7, 0.6, 1.1], 1.5),
    ];
    let springs = vec![
        ImplicitSpring::new(0, 1, 0.8, 50.0, 4.0),
        ImplicitSpring::new(1, 2, 0.9, 30.0, 2.5),
        ImplicitSpring::new(2, 3, 0.5, 70.0, 0.0),
        ImplicitSpring::new(3, 1, 0.7, 20.0, 1.0),
        ImplicitSpring::new(0, 3, 0.4, 40.0, 3.0),
    ];
    (particles, springs)
}

// ---------------------------------------------------------------------------
// F4238: Hessian assembly
// ---------------------------------------------------------------------------

#[test]
fn spring_hessian_has_exact_block_structure_for_a_single_spring() {
    // Two free unit-mass particles; spring along x with l = 2, L0 = 1, k = 10.
    let particles = vec![
        ImplicitParticle::new([0.0, 0.0, 0.0], 1.0),
        ImplicitParticle::new([2.0, 0.0, 0.0], 1.0),
    ];
    let springs = vec![ImplicitSpring::new(0, 1, 1.0, 10.0, 0.0)];
    let dt = 0.5;
    let h = to_dense(&assemble_spring_hessian(&particles, &springs, dt));
    assert_eq!(h.len(), 6);
    let inertia = 1.0 / (dt * dt);
    // d = (-1, 0, 0); transverse stiffness m = 1 - L0/l = 0.5.
    let axial = 10.0;
    let transverse = 10.0 * 0.5;
    let expected_block = [axial, transverse, transverse];
    for axis in 0..3 {
        let kk = expected_block[axis];
        assert!((h[axis][axis] - (inertia + kk)).abs() < 1e-12);
        assert!((h[3 + axis][3 + axis] - (inertia + kk)).abs() < 1e-12);
        assert!((h[axis][3 + axis] + kk).abs() < 1e-12);
        assert!((h[3 + axis][axis] + kk).abs() < 1e-12);
    }
    // No axis coupling for an axis-aligned spring.
    assert!(h[0][1].abs() < 1e-14 && h[1][0].abs() < 1e-14);
}

#[test]
fn spring_hessian_is_symmetric_with_identity_rows_for_static_particles() {
    let (particles, springs) = stretched_cluster();
    let h = to_dense(&assemble_spring_hessian(&particles, &springs, 0.05));
    let n = h.len();
    assert_eq!(n, 12);
    for (r, row) in h.iter().enumerate() {
        for (c, &value) in row.iter().enumerate() {
            assert!(
                (value - h[c][r]).abs() <= 1e-12 * (1.0 + value.abs()),
                "H is not symmetric at ({r},{c}): {value} vs {}",
                h[c][r]
            );
        }
    }
    // Particle 0 is static: identity rows and columns.
    for (r, row) in h.iter().enumerate().take(3) {
        for (c, &value) in row.iter().enumerate() {
            let expected = if r == c { 1.0 } else { 0.0 };
            assert_eq!(value, expected, "static row {r}, column {c}");
            assert_eq!(h[c][r], expected, "static column {r}, row {c}");
        }
    }
}

#[test]
fn spring_hessian_ignores_invalid_and_degenerate_springs() {
    let particles = vec![
        ImplicitParticle::new([0.0, 0.0, 0.0], 1.0),
        ImplicitParticle::new([1.0, 0.0, 0.0], 1.0),
        ImplicitParticle::new([1.0, 0.0, 0.0], 1.0),
    ];
    let springs = vec![
        ImplicitSpring::new(0, 7, 1.0, 100.0, 0.0), // out of range
        ImplicitSpring::new(1, 1, 1.0, 100.0, 0.0), // self loop
        ImplicitSpring::new(1, 2, 1.0, 100.0, 0.0), // zero length
    ];
    let dt = 0.1;
    let h = to_dense(&assemble_spring_hessian(&particles, &springs, dt));
    let expected = 1.0 / (dt * dt);
    for (r, row) in h.iter().enumerate() {
        for (c, &value) in row.iter().enumerate() {
            let want = if r == c { expected } else { 0.0 };
            assert!((value - want).abs() < 1e-12, "entry ({r},{c}) = {value}");
        }
    }
}

#[test]
fn hessian_matches_finite_difference_of_gradient_and_gradient_matches_energy() {
    for with_damping in [false, true] {
        let (particles, springs) = stretched_cluster();
        let x = positions_of(&particles);
        let masses: Vec<f64> = particles.iter().map(|p| p.mass).collect();
        let dt = 0.05;
        let inertia_target = vec![
            [0.0, 0.0, 0.0],
            [1.2, 0.2, -0.1],
            [0.3, 1.1, 0.5],
            [0.6, 0.7, 1.0],
        ];
        let external = vec![
            [0.0; 3],
            [0.3, -1.0, 0.2],
            [0.0, -2.0, 0.0],
            [1.0, 0.0, -0.5],
        ];
        let origin = x.clone();
        let reference = with_damping.then(|| x.clone());
        let problem = ImplicitStepProblem::new(
            masses,
            &springs,
            inertia_target,
            external,
            origin,
            dt,
            reference,
        )
        .expect("consistent problem");

        // Evaluate at a point that differs from the damping reference, so the
        // damper contributes a non-trivial value.
        let eval = vec![
            x[0],
            [1.31, 0.14, -0.18],
            [0.23, 1.17, 0.43],
            [0.66, 0.58, 1.07],
        ];
        let grad = problem.gradient(&eval);
        let hess = to_dense(&problem.hessian(&eval));
        let h = 1e-6;

        // Gradient == finite difference of the energy.
        for (i, &free) in [false, true, true, true].iter().enumerate() {
            for axis in 0..3 {
                let dof = 3 * i + axis;
                if !free {
                    assert_eq!(grad[dof], 0.0, "static dof {dof} must have zero gradient");
                    continue;
                }
                let mut plus = eval.clone();
                let mut minus = eval.clone();
                plus[i][axis] += h;
                minus[i][axis] -= h;
                let fd = (problem.energy(&plus) - problem.energy(&minus)) / (2.0 * h);
                assert!(
                    (grad[dof] - fd).abs() <= 1e-6 * (1.0 + fd.abs()),
                    "damping={with_damping} dof {dof}: grad {} vs fd {fd}",
                    grad[dof]
                );
            }
        }

        // Hessian == finite difference of the gradient on all free dofs.
        for i in 1..4 {
            for axis in 0..3 {
                let col = 3 * i + axis;
                let mut plus = eval.clone();
                let mut minus = eval.clone();
                plus[i][axis] += h;
                minus[i][axis] -= h;
                let gp = problem.gradient(&plus);
                let gm = problem.gradient(&minus);
                for row in 3..12 {
                    let fd = (gp[row] - gm[row]) / (2.0 * h);
                    assert!(
                        (hess[row][col] - fd).abs() <= 1e-5 * (1.0 + fd.abs()),
                        "damping={with_damping} H[{row}][{col}] = {} vs fd {fd}",
                        hess[row][col]
                    );
                }
            }
        }
    }
}

#[test]
fn hessian_of_compressed_springs_is_positive_semi_definite() {
    // Every spring is compressed (l < L0): the exact Hessian has a negative
    // transverse stiffness, which must be clamped so PCG stays applicable.
    let particles = vec![
        ImplicitParticle::new_static([0.0, 0.0, 0.0]),
        ImplicitParticle::new([0.5, 0.1, 0.0], 1.0),
        ImplicitParticle::new([0.1, 0.4, 0.2], 1.0),
    ];
    let springs = vec![
        ImplicitSpring::new(0, 1, 2.0, 1.0e4, 0.0),
        ImplicitSpring::new(1, 2, 2.0, 1.0e4, 0.0),
        ImplicitSpring::new(0, 2, 2.0, 1.0e4, 0.0),
    ];
    let h = assemble_spring_hessian(&particles, &springs, 0.01);
    let mut rng = Lcg(7);
    for _ in 0..50 {
        let v: Vec<f64> = (0..9).map(|_| rng.next_unit() - 0.5).collect();
        let hv = h.matvec(&v);
        let quad: f64 = v.iter().zip(hv.iter()).map(|(a, b)| a * b).sum();
        assert!(quad > 0.0, "v^T H v = {quad} must be positive");
    }
}

// ---------------------------------------------------------------------------
// F4238: Newton-Raphson minimiser
// ---------------------------------------------------------------------------

#[test]
fn newton_matches_closed_form_for_one_free_particle_on_a_spring() {
    // Anchor at the origin, particle hanging below.  With s the distance
    // below the anchor the incremental potential is
    //   E(s) = k/2 (s - L0)^2 + m/(2 dt^2) (s - s0)^2 - f (s - s0),
    // minimised at s = (k L0 + m/dt^2 s0 + f) / (k + m/dt^2).
    let (k, rest, mass, dt) = (1.0e6, 1.0, 2.0, 0.01);
    let weight = mass * 9.81;
    let mut particles = vec![
        ImplicitParticle::new_static([0.0; 3]),
        ImplicitParticle::new([0.0, -1.0, 0.0], mass),
    ];
    particles[1].external_force = [0.0, -weight, 0.0];
    let springs = vec![ImplicitSpring::new(0, 1, rest, k, 0.0)];
    let inertia = positions_of(&particles);
    let solver = NewtonRaphsonSolver::new(20, 1e-9);
    let iterations = solver.minimise(&mut particles, &springs, &inertia, dt);
    assert!((1..=3).contains(&iterations), "iterations = {iterations}");
    let m_dt2 = mass / (dt * dt);
    let expected_s = (k * rest + m_dt2 * 1.0 + weight) / (k + m_dt2);
    let s = -particles[1].position[1];
    assert!(
        (s - expected_s).abs() < 1e-10,
        "s = {s}, expected {expected_s}"
    );
    assert_eq!(particles[0].position, [0.0; 3], "static particle moved");
    assert!(particles[1].position[0].abs() < 1e-12);
    assert!(particles[1].position[2].abs() < 1e-12);
}

#[test]
fn newton_converges_stiff_chain_in_few_iterations() {
    // k = 1e6 chain under gravity: the replaced solver ignored the external
    // force altogether, and a diagonal (Jacobi) gradient descent on the same
    // system is still far from converged (gradient norm ~27) after 50
    // iterations; Newton + PCG converges in a handful.
    let (mut particles, springs) = hanging_chain(21, 0.1, 1.0e6, 0.0);
    for p in particles.iter_mut().skip(1) {
        p.external_force = [0.0, -9.81 * p.mass, 0.0];
    }
    let inertia = positions_of(&particles);
    let dt = 0.01;
    let solver = NewtonRaphsonSolver::new(50, 1e-6);
    let iterations = solver.minimise(&mut particles, &springs, &inertia, dt);
    assert!(
        (1..=5).contains(&iterations),
        "expected <= 5 Newton iterations, got {iterations}"
    );
    let residual = stationarity_residual(&particles, &springs, &inertia, dt);
    assert!(residual < 1e-5, "stationarity residual {residual}");
    // The chain stretches under its own weight.
    assert!(particles[20].position[1] < inertia[20][1]);
    assert_eq!(particles[0].position, [0.0; 3]);
}

#[test]
fn newton_reads_the_pcg_parameters() {
    let build = || {
        let (mut particles, springs) = hanging_chain(21, 0.1, 1.0e6, 0.0);
        for p in particles.iter_mut().skip(1) {
            p.external_force = [0.0, -9.81 * p.mass, 0.0];
        }
        (particles, springs)
    };
    let dt = 0.01;

    let (mut full, springs) = build();
    let inertia = positions_of(&full);
    let precise = NewtonRaphsonSolver::new(60, 1e-6);
    let precise_iterations = precise.minimise(&mut full, &springs, &inertia, dt);

    let (mut crude, springs) = build();
    let mut crippled = NewtonRaphsonSolver::new(60, 1e-6);
    crippled.pcg_params = PcgParams {
        max_iter: 1,
        tolerance: 1e-8,
    };
    let crippled_iterations = crippled.minimise(&mut crude, &springs, &inertia, dt);

    assert!(
        crippled_iterations > precise_iterations,
        "one PCG iteration per Newton step must need more Newton iterations \
         ({crippled_iterations}) than a converged PCG solve ({precise_iterations})"
    );
}

#[test]
fn newton_applies_external_force_and_keeps_static_particles_fixed() {
    // No springs: x = y + dt^2 / m * f for every free particle.
    let dt = 0.5;
    let mut particles = vec![
        ImplicitParticle::new([1.0, 2.0, 3.0], 2.0),
        ImplicitParticle::new_static([5.0, 5.0, 5.0]),
        ImplicitParticle::new([0.0, 0.0, 0.0], 4.0),
    ];
    particles[0].external_force = [4.0, 0.0, -8.0];
    particles[1].external_force = [100.0, 100.0, 100.0];
    particles[2].external_force = [0.0, 16.0, 0.0];
    let inertia = positions_of(&particles);
    let solver = NewtonRaphsonSolver::new(10, 1e-12);
    let iterations = solver.minimise(&mut particles, &[], &inertia, dt);
    assert!(iterations >= 1);
    let expected0 = [1.0 + 0.25 / 2.0 * 4.0, 2.0, 3.0 + 0.25 / 2.0 * -8.0];
    for (position, expected) in particles[0].position.iter().zip(expected0.iter()) {
        assert!((position - expected).abs() < 1e-12);
    }
    assert_eq!(particles[1].position, [5.0, 5.0, 5.0]);
    assert!((particles[2].position[1] - 0.25 / 4.0 * 16.0).abs() < 1e-12);
}

#[test]
fn newton_rejects_inconsistent_inputs_without_touching_particles() {
    let (mut particles, springs) = hanging_chain(4, 1.0, 100.0, 0.0);
    let inertia = positions_of(&particles);
    let before = positions_of(&particles);
    let solver = NewtonRaphsonSolver::new(10, 1e-9);

    for bad_dt in [0.0, -0.1, f64::NAN, f64::INFINITY] {
        assert_eq!(
            solver.minimise(&mut particles, &springs, &inertia, bad_dt),
            0
        );
        assert_eq!(positions_of(&particles), before);
    }
    // Inertia target of the wrong length.
    assert_eq!(
        solver.minimise(&mut particles, &springs, &inertia[..2], 0.1),
        0
    );
    assert_eq!(positions_of(&particles), before);
    // Spring referencing a non-existent particle.
    let bad_springs = vec![ImplicitSpring::new(0, 9, 1.0, 10.0, 0.0)];
    assert_eq!(
        solver.minimise(&mut particles, &bad_springs, &inertia, 0.1),
        0
    );
    assert_eq!(positions_of(&particles), before);
}

#[test]
fn newton_handles_a_two_thousand_particle_chain() {
    // 2000 particles / 1999 springs.  The previous implementation scanned every
    // spring for every particle (O(n * S)) on each residual evaluation.
    let n = 2000;
    let (mut particles, springs) = hanging_chain(n, 0.05, 1.0e5, 0.0);
    for p in particles.iter_mut().skip(1) {
        p.external_force = [0.0, -9.81 * p.mass, 0.0];
    }
    let inertia = positions_of(&particles);
    let dt = 0.02;
    let mut solver = NewtonRaphsonSolver::new(8, 1e-4);
    solver.pcg_params = PcgParams {
        max_iter: 5000,
        tolerance: 1e-10,
    };
    let iterations = solver.minimise(&mut particles, &springs, &inertia, dt);
    assert!(iterations >= 1);
    let residual = stationarity_residual(&particles, &springs, &inertia, dt);
    assert!(residual < 1e-2, "stationarity residual {residual}");
    assert!(
        particles
            .iter()
            .all(|p| p.position.iter().all(|c| c.is_finite()))
    );
}

// ---------------------------------------------------------------------------
// F4238: backward Euler
// ---------------------------------------------------------------------------

#[test]
fn backward_euler_free_particle_matches_the_analytic_implicit_step() {
    // Constant force: v1 = v0 + dt (g + f / m), x1 = x0 + dt v1 (exact for
    // backward Euler), and the external force must be included.
    let mut particles = vec![ImplicitParticle::new([1.0, 2.0, 3.0], 2.0)];
    particles[0].velocity = [0.5, -0.25, 1.0];
    particles[0].external_force = [6.0, 0.0, -4.0];
    let gravity = [0.0, -9.81, 0.0];
    let dt = 0.1;
    let integrator = BackwardEulerIntegrator::new(10, 1e-10);
    let iterations = integrator.step(&mut particles, &[], gravity, dt);
    assert!(iterations >= 1);
    let v0 = [0.5, -0.25, 1.0];
    // Acceleration: g + f / m = (0, -9.81, 0) + (6, 0, -4) / 2.
    let a = [3.0, -9.81, -2.0];
    let expected_v = [v0[0] + dt * a[0], v0[1] + dt * a[1], v0[2] + dt * a[2]];
    let expected_x = [
        1.0 + dt * expected_v[0],
        2.0 + dt * expected_v[1],
        3.0 + dt * expected_v[2],
    ];
    for axis in 0..3 {
        assert!(
            (particles[0].velocity[axis] - expected_v[axis]).abs() < 1e-10,
            "v[{axis}] = {} vs {}",
            particles[0].velocity[axis],
            expected_v[axis]
        );
        assert!((particles[0].position[axis] - expected_x[axis]).abs() < 1e-10);
    }
}

#[test]
fn backward_euler_stiff_chain_converges_to_static_equilibrium() {
    // Unconditionally stable even with omega * dt >> 1, and the fixed point of
    // the iteration is the static equilibrium: spring j carries the weight of
    // everything below it, so its stretch is (n - 1 - j) m g / k.
    let n = 11;
    let (k, rest, dt) = (1.0e6, 0.1, 0.02);
    let (mut particles, springs) = hanging_chain(n, rest, k, 50.0);
    let integrator = BackwardEulerIntegrator::new(20, 1e-9);
    let gravity = [0.0, -9.81, 0.0];
    let mut worst_iterations = 0;
    for _ in 0..600 {
        let iterations = integrator.step(&mut particles, &springs, gravity, dt);
        worst_iterations = worst_iterations.max(iterations);
        assert!(
            particles
                .iter()
                .all(|p| p.position.iter().all(|c| c.is_finite())),
            "simulation blew up"
        );
    }
    assert!(
        worst_iterations <= 10,
        "worst Newton count {worst_iterations}"
    );
    let mut expected_y = 0.0;
    for j in 0..n - 1 {
        let tension = (n - 1 - j) as f64 * 9.81;
        expected_y -= rest + tension / k;
        let y = particles[j + 1].position[1];
        assert!(
            (y - expected_y).abs() < 1e-7,
            "particle {} at {y}, expected {expected_y}",
            j + 1
        );
    }
    let kinetic: f64 = particles
        .iter()
        .map(|p| 0.5 * p.mass * v3_dot(p.velocity, p.velocity))
        .sum();
    assert!(kinetic < 1e-8, "chain did not settle, KE = {kinetic}");
}

#[test]
fn backward_euler_damper_dissipates_energy_and_undamped_spring_does_not_gain() {
    let run = |damping: f64| -> (f64, f64) {
        let mut particles = vec![
            ImplicitParticle::new_static([0.0; 3]),
            ImplicitParticle::new([1.5, 0.0, 0.0], 1.0),
        ];
        let springs = vec![ImplicitSpring::new(0, 1, 1.0, 400.0, damping)];
        let integrator = BackwardEulerIntegrator::new(20, 1e-10);
        let dt = 0.01;
        let energy = |particles: &[ImplicitParticle]| {
            0.5 * v3_dot(particles[1].velocity, particles[1].velocity)
                + springs[0].potential_energy(particles)
        };
        let e0 = energy(&particles);
        let mut max_energy = e0;
        for _ in 0..300 {
            integrator.step(&mut particles, &springs, [0.0; 3], dt);
            max_energy = max_energy.max(energy(&particles));
        }
        (energy(&particles), max_energy / e0)
    };
    let (e_undamped, growth_undamped) = run(0.0);
    let (e_damped, _) = run(8.0);
    // Backward Euler never increases the energy of an undamped spring.
    assert!(
        growth_undamped <= 1.0 + 1e-9,
        "energy grew by {growth_undamped}"
    );
    // The damper removes considerably more energy than the scheme itself.
    assert!(
        e_damped < 0.25 * e_undamped,
        "damped {e_damped} vs undamped {e_undamped}"
    );
}

#[test]
fn backward_euler_uses_per_particle_external_force_and_pcg_parameters() {
    let (mut with_force, springs) = hanging_chain(6, 0.5, 5.0e4, 0.0);
    let (mut without_force, _) = hanging_chain(6, 0.5, 5.0e4, 0.0);
    with_force[5].external_force = [30.0, 0.0, 0.0];
    let integrator = BackwardEulerIntegrator::new(20, 1e-9);
    for _ in 0..20 {
        integrator.step(&mut with_force, &springs, [0.0; 3], 0.01);
        integrator.step(&mut without_force, &springs, [0.0; 3], 0.01);
    }
    assert!(
        with_force[5].position[0] > 1e-4,
        "external force had no effect: x = {}",
        with_force[5].position[0]
    );
    assert!(without_force[5].position[0].abs() < 1e-12);
    assert_eq!(with_force[0].position, [0.0; 3]);

    // The `pcg` field is the linear solver: a crippled one needs more Newton
    // iterations for the same step.
    let (mut a, springs) = hanging_chain(21, 0.1, 1.0e6, 0.0);
    let (mut b, _) = hanging_chain(21, 0.1, 1.0e6, 0.0);
    let good = BackwardEulerIntegrator::new(60, 1e-6);
    let mut poor = BackwardEulerIntegrator::new(60, 1e-6);
    poor.pcg = PcgSolver::new(PcgParams {
        max_iter: 1,
        tolerance: 1e-8,
    });
    let gravity = [0.0, -9.81, 0.0];
    let good_iterations = good.step(&mut a, &springs, gravity, 0.01);
    let poor_iterations = poor.step(&mut b, &springs, gravity, 0.01);
    assert!(
        poor_iterations > good_iterations,
        "poor PCG {poor_iterations} vs good PCG {good_iterations}"
    );
}

/// `side x side` cloth with structural and shear springs, pinned at the two
/// top corners.
fn hanging_cloth(
    side: usize,
    spacing: f64,
    mass: f64,
    stiffness: f64,
    damping: f64,
) -> (Vec<ImplicitParticle>, Vec<ImplicitSpring>) {
    let mut particles = Vec::with_capacity(side * side);
    for r in 0..side {
        for c in 0..side {
            let position = [c as f64 * spacing, -(r as f64) * spacing, 0.0];
            if r == 0 && (c == 0 || c + 1 == side) {
                particles.push(ImplicitParticle::new_static(position));
            } else {
                particles.push(ImplicitParticle::new(position, mass));
            }
        }
    }
    let index = |r: usize, c: usize| r * side + c;
    let diagonal = spacing * std::f64::consts::SQRT_2;
    let mut springs = Vec::new();
    for r in 0..side {
        for c in 0..side {
            if c + 1 < side {
                springs.push(ImplicitSpring::new(
                    index(r, c),
                    index(r, c + 1),
                    spacing,
                    stiffness,
                    damping,
                ));
            }
            if r + 1 < side {
                springs.push(ImplicitSpring::new(
                    index(r, c),
                    index(r + 1, c),
                    spacing,
                    stiffness,
                    damping,
                ));
            }
            if r + 1 < side && c + 1 < side {
                springs.push(ImplicitSpring::new(
                    index(r, c),
                    index(r + 1, c + 1),
                    diagonal,
                    stiffness,
                    damping,
                ));
                springs.push(ImplicitSpring::new(
                    index(r, c + 1),
                    index(r + 1, c),
                    diagonal,
                    stiffness,
                    damping,
                ));
            }
        }
    }
    (particles, springs)
}

#[test]
fn backward_euler_hanging_cloth_converges_in_a_few_newton_iterations() {
    // The diagonal (Jacobi) iteration this replaced exhausts the iteration cap
    // on a coupled 10 x 10 cloth; Newton + PCG needs a handful of iterations
    // per step even for very stiff springs.
    for stiffness in [500.0, 5.0e4, 1.0e6] {
        let (mut particles, springs) = hanging_cloth(10, 0.1, 0.01, stiffness, 0.05);
        let integrator = BackwardEulerIntegrator::new(50, 1e-6);
        let mut worst = 0;
        for _ in 0..100 {
            worst =
                worst.max(integrator.step(&mut particles, &springs, [0.0, -9.81, 0.0], 1.0 / 60.0));
        }
        assert!(
            worst <= 8,
            "k = {stiffness}: worst Newton iteration count {worst}"
        );
        assert!(
            particles
                .iter()
                .all(|p| p.position.iter().all(|c| c.is_finite()))
        );
        // The pinned corners hold the cloth, which hangs about 0.9 m below them.
        let lowest = particles
            .iter()
            .map(|p| p.position[1])
            .fold(f64::INFINITY, f64::min);
        assert!(
            (-1.3..-0.85).contains(&lowest),
            "k = {stiffness}: lowest particle at y = {lowest}"
        );
        assert_eq!(particles[0].position, [0.0, 0.0, 0.0]);
    }
}

#[test]
fn backward_euler_rejects_inconsistent_inputs() {
    let (mut particles, springs) = hanging_chain(3, 1.0, 100.0, 1.0);
    particles[1].velocity = [1.0, 0.0, 0.0];
    let before_x = positions_of(&particles);
    let before_v: Vec<[f64; 3]> = particles.iter().map(|p| p.velocity).collect();
    let integrator = BackwardEulerIntegrator::new(10, 1e-9);
    for bad_dt in [0.0, -1.0, f64::NAN] {
        assert_eq!(
            integrator.step(&mut particles, &springs, [0.0, -9.81, 0.0], bad_dt),
            0
        );
        assert_eq!(positions_of(&particles), before_x);
        let after_v: Vec<[f64; 3]> = particles.iter().map(|p| p.velocity).collect();
        assert_eq!(after_v, before_v);
    }
    let bad_springs = vec![ImplicitSpring::new(0, 5, 1.0, 10.0, 0.0)];
    assert_eq!(
        integrator.step(&mut particles, &bad_springs, [0.0, -9.81, 0.0], 0.1),
        0
    );
    assert_eq!(positions_of(&particles), before_x);
}

#[test]
fn newton_pcg_minimise_terminates_immediately_at_a_stationary_point() {
    // Chain already in equilibrium without forces: zero gradient, zero updates.
    let (particles, springs) = hanging_chain(5, 1.0, 1.0e3, 0.0);
    let masses: Vec<f64> = particles.iter().map(|p| p.mass).collect();
    let x0 = positions_of(&particles);
    let problem = ImplicitStepProblem::new(
        masses,
        &springs,
        x0.clone(),
        vec![[0.0; 3]; 5],
        x0.clone(),
        0.1,
        None,
    )
    .expect("consistent problem");
    let mut x = x0.clone();
    let pcg = PcgSolver::new(PcgParams::default());
    let iterations = newton_pcg_minimise(&problem, &mut x, 10, 1e-9, &pcg, &LineSearch::default());
    assert_eq!(iterations, 0);
    assert_eq!(x, x0);
    assert!(vec_norm(&problem.gradient(&x)) < 1e-12);
}

// ---------------------------------------------------------------------------
// F4241: Projective Dynamics
// ---------------------------------------------------------------------------

/// Verbatim copy of the pre-fix `ProjectiveDynamicsSolver::step`, which
/// re-derived the constant diagonal for every particle by scanning every
/// constraint inside every iteration.
fn reference_projective_step(
    solver: &ProjectiveDynamicsSolver,
    positions: &mut [[f64; 3]],
    velocities: &mut [[f64; 3]],
    gravity: [f64; 3],
) {
    let n = solver.n_particles;
    let dt = solver.dt;
    let y = solver.inertia_target(positions, velocities, gravity);
    let mut q = y.clone();
    for _iter in 0..solver.n_iterations {
        let mut rhs: Vec<[f64; 3]> = (0..n)
            .map(|i| v3_scale(y[i], solver.masses[i] / (dt * dt)))
            .collect();
        for c in &solver.constraints {
            let proj = c.local_step(&q);
            match &c.kind {
                ProjectiveConstraintKind::Spring { i, j, weight, .. } => {
                    let w_dt2 = weight / (dt * dt);
                    rhs[*i] = v3_add(rhs[*i], v3_scale(proj[0], w_dt2));
                    rhs[*j] = v3_add(rhs[*j], v3_scale(proj[1], w_dt2));
                }
                ProjectiveConstraintKind::Anchor { i, weight, .. } => {
                    let w_dt2 = weight / (dt * dt);
                    rhs[*i] = v3_add(rhs[*i], v3_scale(proj[0], w_dt2));
                }
            }
        }
        for i in 0..n {
            let mut diag = solver.masses[i] / (dt * dt);
            for c in &solver.constraints {
                match &c.kind {
                    ProjectiveConstraintKind::Spring {
                        i: ci,
                        j: cj,
                        weight,
                        ..
                    } => {
                        if *ci == i || *cj == i {
                            diag += weight / (dt * dt);
                        }
                    }
                    ProjectiveConstraintKind::Anchor { i: ci, weight, .. } => {
                        if *ci == i {
                            diag += weight / (dt * dt);
                        }
                    }
                }
            }
            if diag > 1e-15 {
                q[i] = v3_scale(rhs[i], 1.0 / diag);
            }
        }
    }
    for i in 0..n {
        if solver.masses[i] > 1e-15 {
            velocities[i] = v3_scale(v3_sub(q[i], positions[i]), 1.0 / dt);
            positions[i] = q[i];
        }
    }
}

/// `side x side` jittered grid with structural + shear springs, the top row
/// pinned by anchors, a few static particles and one self-referencing spring.
fn grid_solver(side: usize, dt: f64, iterations: usize) -> ProjectiveDynamicsSolver {
    let n = side * side;
    let mut masses = vec![1.0f64; n];
    for (i, m) in masses.iter_mut().enumerate() {
        *m = 0.5 + 0.1 * (i % 7) as f64;
    }
    masses[side * side - 1] = 0.0;
    let mut solver = ProjectiveDynamicsSolver::new(n, masses, dt, iterations);
    let idx = |r: usize, c: usize| r * side + c;
    for r in 0..side {
        for c in 0..side {
            if c + 1 < side {
                solver.add_constraint(ProjectiveConstraint::spring(
                    idx(r, c),
                    idx(r, c + 1),
                    1.0,
                    800.0,
                ));
            }
            if r + 1 < side {
                solver.add_constraint(ProjectiveConstraint::spring(
                    idx(r, c),
                    idx(r + 1, c),
                    1.0,
                    800.0,
                ));
            }
            if r + 1 < side && c + 1 < side {
                solver.add_constraint(ProjectiveConstraint::spring(
                    idx(r, c),
                    idx(r + 1, c + 1),
                    std::f64::consts::SQRT_2,
                    250.0,
                ));
            }
        }
    }
    for c in 0..side {
        solver.add_constraint(ProjectiveConstraint::anchor(
            idx(0, c),
            [c as f64, 0.0, 0.0],
            5.0e4,
        ));
    }
    solver.add_constraint(ProjectiveConstraint::spring(3, 3, 0.0, 10.0));
    solver
}

fn grid_state(side: usize, seed: u64) -> (Vec<[f64; 3]>, Vec<[f64; 3]>) {
    let mut rng = Lcg(seed);
    let mut positions = Vec::with_capacity(side * side);
    let mut velocities = Vec::with_capacity(side * side);
    for r in 0..side {
        for c in 0..side {
            positions.push([
                c as f64 + 0.2 * (rng.next_unit() - 0.5),
                -(r as f64) + 0.2 * (rng.next_unit() - 0.5),
                0.2 * (rng.next_unit() - 0.5),
            ]);
            velocities.push([
                rng.next_unit() - 0.5,
                rng.next_unit() - 0.5,
                rng.next_unit() - 0.5,
            ]);
        }
    }
    (positions, velocities)
}

#[test]
fn projective_dynamics_is_bit_identical_to_the_quadratic_reference() {
    let side = 10;
    let solver = grid_solver(side, 1.0 / 60.0, 12);
    let gravity = [0.0, -9.81, 0.0];
    let (mut pos_new, mut vel_new) = grid_state(side, 12345);
    let mut pos_ref = pos_new.clone();
    let mut vel_ref = vel_new.clone();
    for step in 0..8 {
        solver.step(&mut pos_new, &mut vel_new, gravity);
        reference_projective_step(&solver, &mut pos_ref, &mut vel_ref, gravity);
        assert_eq!(pos_new, pos_ref, "positions differ after step {step}");
        assert_eq!(vel_new, vel_ref, "velocities differ after step {step}");
    }
    // The comparison is meaningful: the grid actually moved.
    let (start, _) = grid_state(side, 12345);
    assert_ne!(pos_new, start);
    // Static particle (zero mass) never moves.
    assert_eq!(pos_new[side * side - 1], start[side * side - 1]);
}

#[test]
fn projective_dynamics_local_step_returns_fixed_size_projections() {
    let spring = ProjectiveConstraint::spring(0, 1, 2.0, 10.0);
    let q = vec![[0.0, 0.0, 0.0], [6.0, 0.0, 0.0]];
    let p = spring.local_step(&q);
    assert!((v3_norm(v3_sub(p[0], p[1])) - 2.0).abs() < 1e-12);
    // Symmetric about the midpoint of the current end points.
    let mid = v3_scale(v3_add(p[0], p[1]), 0.5);
    assert!(v3_norm(v3_sub(mid, [3.0, 0.0, 0.0])) < 1e-12);
    assert!(p[0][0] < p[1][0], "p_i must stay on the side of x_i");

    let coincident = vec![[1.0, 1.0, 1.0], [1.0, 1.0, 1.0]];
    let p = spring.local_step(&coincident);
    assert_eq!(p, [[1.0, 1.0, 1.0], [1.0, 1.0, 1.0]]);

    let anchor = ProjectiveConstraint::anchor(1, [7.0, 8.0, 9.0], 3.0);
    assert_eq!(anchor.local_step(&q), [[7.0, 8.0, 9.0], [0.0; 3]]);
}

#[test]
fn projective_dynamics_two_anchors_give_the_weighted_mean() {
    // Anchor projections do not depend on the iterate, so the global solve is
    // exact after one iteration: q = (m y + sum w t) / (m + sum w), where the
    // 1 / dt^2 factor of the inertia and of the weights cancels.  The result
    // therefore depends on the complete (constant) diagonal.
    let dt = 0.1;
    let mut solver = ProjectiveDynamicsSolver::new(1, vec![2.0], dt, 3);
    solver.add_constraint(ProjectiveConstraint::anchor(0, [0.0, 0.0, 0.0], 30.0));
    solver.add_constraint(ProjectiveConstraint::anchor(0, [4.0, 0.0, 0.0], 20.0));
    let mut positions = vec![[1.0, 0.0, 0.0]];
    let mut velocities = vec![[0.0; 3]];
    solver.step(&mut positions, &mut velocities, [0.0; 3]);
    let expected = (2.0 * 1.0 + 30.0 * 0.0 + 20.0 * 4.0) / (2.0 + 30.0 + 20.0);
    assert!(
        (positions[0][0] - expected).abs() < 1e-12,
        "x = {}, expected {expected}",
        positions[0][0]
    );
    assert!((velocities[0][0] - (expected - 1.0) / dt).abs() < 1e-10);
    assert_eq!(positions[0][1], 0.0);
    assert_eq!(positions[0][2], 0.0);
}

#[test]
fn projective_dynamics_scales_linearly_with_constraint_count() {
    // n = 40_000 particles, ~119_000 constraints.  Rebuilding the diagonal by
    // scanning every constraint for every particle (the previous O(n * C)
    // implementation) visits ~5 * 10^9 constraints per iteration: measured at
    // ~2 minutes for this step in an optimised build, against ~0.1 s for the
    // single-pass assembly.  The bound is ~100x the expected run time, so it
    // is robust against machine load.
    let side = 200;
    let solver = grid_solver(side, 1.0 / 60.0, 6);
    assert!(solver.constraints.len() > 119_000);
    let (mut positions, mut velocities) = grid_state(side, 99);
    let started = std::time::Instant::now();
    solver.step(&mut positions, &mut velocities, [0.0, -9.81, 0.0]);
    let elapsed = started.elapsed();
    assert!(
        elapsed < std::time::Duration::from_secs(10),
        "one projective dynamics step took {elapsed:?}"
    );
    assert!(positions.iter().all(|p| p.iter().all(|c| c.is_finite())));
    assert!(velocities.iter().any(|v| v3_norm(*v) > 0.0));
}
