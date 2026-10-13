// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

// Unit tests for the DFSPH solver (`super`), kept in their own file so the
// solver module stays below the 2000-line limit.

use super::*;

fn make_default_config() -> DfsphConfig {
    DfsphConfig {
        rest_density: 1000.0,
        smoothing_length: 0.1,
        support_radius: 0.2,
        max_pressure_iterations: 50,
        max_divergence_iterations: 50,
        error_threshold: 1e-3,
        omega: 0.5,
        viscosity: 1e-6,
        gravity: [0.0, -9.81, 0.0],
        surface_tension: 0.0,
    }
}

#[test]
fn test_wendland_c2_at_zero() {
    let h = 0.1;
    let w = wendland_c2_kernel(0.0, h);
    // At r=0: W(0,h) = 21/(16π h³), the canonical Wendland C2 constant
    // 21/(2π H³) written for the support radius H = 2h.
    let expected = 21.0 / (16.0 * std::f64::consts::PI * h * h * h);
    assert!(
        (w - expected).abs() < 1e-10,
        "W(0,h) should equal {}, got {}",
        expected,
        w
    );
}

#[test]
fn test_wendland_c2_at_support() {
    let h = 0.1;
    let w = wendland_c2_kernel(2.0 * h, h);
    assert!(w.abs() < 1e-14, "W(2h,h) should be zero, got {}", w);
}

#[test]
fn test_wendland_c2_positive() {
    let h = 0.1;
    for r in [0.0, 0.05, 0.1, 0.15, 0.19] {
        let w = wendland_c2_kernel(r, h);
        assert!(w >= 0.0, "W(r,h) should be non-negative for r={}", r);
    }
}

#[test]
fn test_wendland_c4_at_zero() {
    let h = 0.1;
    let w = wendland_c4_kernel(0.0, h);
    // Canonical Wendland C4 constant 495/(32π H³) with support radius H = 2h.
    let expected = 495.0 / (256.0 * std::f64::consts::PI * h * h * h);
    assert!(
        (w - expected).abs() < 1e-9 * expected,
        "W_C4(0,h) should equal {}, got {}",
        expected,
        w
    );
}

#[test]
fn test_wendland_c4_at_support() {
    let h = 0.1;
    let w = wendland_c4_kernel(2.0 * h + 1e-12, h);
    assert!(w.abs() < 1e-14, "W_C4 should be zero beyond support");
}

#[test]
fn test_wendland_c2_gradient_antisymmetric() {
    let h = 0.1;
    let r_vec = [0.05, 0.0, 0.0];
    let r_neg = [-0.05, 0.0, 0.0];
    let g1 = wendland_c2_gradient(r_vec, h);
    let g2 = wendland_c2_gradient(r_neg, h);
    assert!(
        (g1[0] + g2[0]).abs() < 1e-10,
        "gradient should be antisymmetric"
    );
}

#[test]
fn test_dfsph_config_default() {
    let cfg = DfsphConfig::default();
    assert!(cfg.rest_density > 0.0);
    assert!(cfg.smoothing_length > 0.0);
    assert!(cfg.max_pressure_iterations > 0);
}

#[test]
fn test_dfsph_particle_new() {
    let p = DfsphParticle::new([1.0, 2.0, 3.0], 0.001, 1000.0);
    assert_eq!(p.position, [1.0, 2.0, 3.0]);
    assert_eq!(p.mass, 0.001);
    assert_eq!(p.density, 1000.0);
    assert!(!p.is_boundary);
}

#[test]
fn test_dfsph_particle_kinetic_energy() {
    let mut p = DfsphParticle::new([0.0; 3], 1.0, 1000.0);
    p.velocity = [3.0, 4.0, 0.0];
    let ke = p.kinetic_energy();
    assert!(
        (ke - 12.5).abs() < 1e-10,
        "KE should be 0.5*1*25=12.5, got {}",
        ke
    );
}

#[test]
fn test_dfsph_solver_new() {
    let cfg = make_default_config();
    let solver = DfsphSolver::new(cfg);
    assert_eq!(solver.num_particles(), 0);
}

#[test]
fn test_dfsph_solver_add_particle() {
    let cfg = make_default_config();
    let mut solver = DfsphSolver::new(cfg);
    solver.add_particle([0.0, 0.0, 0.0], [0.0, 0.0, 0.0], 0.001);
    assert_eq!(solver.num_particles(), 1);
}

#[test]
fn test_dfsph_solver_kinetic_energy() {
    let cfg = make_default_config();
    let mut solver = DfsphSolver::new(cfg);
    solver.add_particle([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], 2.0);
    let ke = solver.kinetic_energy();
    assert!((ke - 1.0).abs() < 1e-10);
}

#[test]
fn test_dfsph_stats_new() {
    let stats = DfsphStats::new();
    assert!(stats.density_error_history.is_empty());
}

#[test]
fn test_dfsph_stats_record_density() {
    let mut stats = DfsphStats::new();
    stats.record_density_step(0.01, 5);
    assert_eq!(stats.density_error_history.len(), 1);
    assert_eq!(stats.pressure_iterations[0], 5);
}

#[test]
fn test_dfsph_stats_mean_error() {
    let mut stats = DfsphStats::new();
    stats.record_density_step(0.02, 3);
    stats.record_density_step(0.04, 4);
    let mean = stats.mean_density_error();
    assert!((mean - 0.03).abs() < 1e-10, "mean error should be 0.03");
}

#[test]
fn test_wcsph_density_compute() {
    let density_obj = WcSphDensity::new(KernelType::WendlandC2, 0.1);
    let pos_i = [0.0; 3];
    let neighbors = vec![([0.05, 0.0, 0.0], 0.001), ([-0.05, 0.0, 0.0], 0.001)];
    let rho = density_obj.compute_density(pos_i, &neighbors);
    assert!(rho > 0.0, "density should be positive");
}

#[test]
fn test_compute_alpha_zero_neighbors() {
    let p = DfsphParticle::new([0.0; 3], 0.001, 1000.0);
    let alpha = compute_alpha(
        p.density,
        p.position,
        p.mass,
        &[],
        0.1,
        &KernelType::WendlandC2,
    );
    assert_eq!(alpha, 0.0, "alpha with no neighbors should be 0");
}

#[test]
fn test_pressure_solve_density_change_rate() {
    let solver_iter = PressureSolveIter::new(1000.0, 50, 1e-3, 0.5, KernelType::WendlandC2, 0.1);
    let pos_i = [0.0; 3];
    let vel_i = [1.0, 0.0, 0.0];
    let rho_i = 1000.0;
    let neighbors = vec![([0.05, 0.0, 0.0], [0.5, 0.0, 0.0], 0.001, 1000.0)];
    let rate = solver_iter.density_change_rate(pos_i, vel_i, rho_i, &neighbors);
    // Particle i closes in on its neighbour (1.0 > 0.5 m/s toward +x): compression.
    assert!(rate > 0.0, "an approaching neighbour must compress: {rate}");
    // A neighbour running away faster than i follows: expansion.
    let receding = vec![([0.05, 0.0, 0.0], [2.0, 0.0, 0.0], 0.001, 1000.0)];
    let expansion = solver_iter.density_change_rate(pos_i, vel_i, rho_i, &receding);
    assert!(
        expansion < 0.0,
        "a receding neighbour must expand: {expansion}"
    );
}

#[test]
fn test_divergence_solve_velocity_divergence() {
    let div_solver = DivergenceSolveIter::new(50, 1e-3, 0.5, KernelType::WendlandC2, 0.1);
    let pos_i = [0.0; 3];
    let vel_i = [0.0; 3];
    let rho_i = 1000.0;
    let neighbors = vec![([0.05, 0.0, 0.0], [0.1, 0.0, 0.0], 0.001, 1000.0)];
    let div = div_solver.velocity_divergence(pos_i, vel_i, rho_i, &neighbors);
    // The neighbour moves away from the particle at rest: it is being rarefied,
    // so the density-change rate over density is negative ...
    assert!(div < 0.0, "a receding neighbour must expand: {div}");
    // ... and a neighbour moving toward it compresses it (positive).
    let approaching = vec![([0.05, 0.0, 0.0], [-0.1, 0.0, 0.0], 0.001, 1000.0)];
    let compression = div_solver.velocity_divergence(pos_i, vel_i, rho_i, &approaching);
    assert!(
        compression > 0.0,
        "an approaching neighbour must compress: {compression}"
    );
    assert!((compression + div).abs() < 1e-12 * compression.abs());
}

/// Build a symmetric `n_side³` grid of fluid particles with measured SPH
/// densities and neighbor lists, plus an inward (compressing) radial
/// velocity field of magnitude `v_in`. Returns the particles, neighbor
/// lists, and the index of the geometric centre particle.
fn build_compressing_cluster(
    n_side: usize,
    spacing: f64,
    h: f64,
    mass: f64,
    rho0: f64,
    v_in: f64,
    kernel: &KernelType,
) -> (Vec<DfsphParticle>, Vec<Vec<usize>>, usize) {
    let mut particles = Vec::new();
    for ix in 0..n_side {
        for iy in 0..n_side {
            for iz in 0..n_side {
                let pos = [
                    ix as f64 * spacing,
                    iy as f64 * spacing,
                    iz as f64 * spacing,
                ];
                particles.push(DfsphParticle::new(pos, mass, rho0));
            }
        }
    }
    let n = particles.len();
    let support = 2.0 * h;
    let density_obj = WcSphDensity::new(kernel.clone(), h);
    let mut neighbor_lists = vec![Vec::new(); n];
    for i in 0..n {
        for j in 0..n {
            if i == j {
                continue;
            }
            let r = norm3(sub3(particles[i].position, particles[j].position));
            if r < support {
                neighbor_lists[i].push(j);
            }
        }
    }
    // Measure the real SPH density (self-contribution + neighbours).
    for i in 0..n {
        let pos_i = particles[i].position;
        let mut rho = mass * density_obj.eval_kernel(0.0);
        for &j in &neighbor_lists[i] {
            let r = norm3(sub3(pos_i, particles[j].position));
            rho += mass * density_obj.eval_kernel(r);
        }
        particles[i].density = rho;
    }
    // Centre of the cluster and an inward radial velocity that actively
    // compresses it (positive Dρ/Dt at the interior).
    let c = (n_side as f64 - 1.0) * spacing / 2.0;
    let center = [c, c, c];
    let mut center_idx = 0;
    let mut best = f64::INFINITY;
    for (i, p) in particles.iter_mut().enumerate() {
        let to_center = sub3(center, p.position);
        let r = norm3(to_center);
        if r > 1e-9 {
            p.velocity = scale3(to_center, v_in / r);
        }
        if r < best {
            best = r;
            center_idx = i;
        }
    }
    (particles, neighbor_lists, center_idx)
}

/// Rate of density change `Dρ_i/Dt = Σ_j m_j (v_i − v_j) · ∇W_ij` of particle
/// `i` for the given `velocities`: the exact time derivative of the density
/// summation `ρ_i = Σ_j m_j W(|x_i − x_j|)`, positive where `i` is being
/// compressed. Evaluated independently of the solvers under test.
fn density_rate_at(
    particles: &[DfsphParticle],
    velocities: &[[f64; 3]],
    neighbor_lists: &[Vec<usize>],
    i: usize,
    kernel: &KernelType,
    h: f64,
) -> f64 {
    let density_obj = WcSphDensity::new(kernel.clone(), h);
    neighbor_lists[i]
        .iter()
        .filter(|&&j| j != i)
        .map(|&j| {
            let g = density_obj.eval_gradient(sub3(particles[i].position, particles[j].position));
            particles[j].mass * dot3(sub3(velocities[i], velocities[j]), g)
        })
        .sum()
}

/// Predicted density `ρ*_i = ρ_i + dt · Dρ_i/Dt` of particle `i`, the quantity
/// the constant-density solve must drive to the rest density.
fn predicted_density_at(
    particles: &[DfsphParticle],
    velocities: &[[f64; 3]],
    neighbor_lists: &[Vec<usize>],
    i: usize,
    kernel: &KernelType,
    h: f64,
    dt: f64,
) -> f64 {
    particles[i].density + dt * density_rate_at(particles, velocities, neighbor_lists, i, kernel, h)
}

/// Mean compressive density error `mean_i max(ρ*_i − ρ₀, 0) / ρ₀`.
fn mean_compressive_density_error(
    particles: &[DfsphParticle],
    velocities: &[[f64; 3]],
    neighbor_lists: &[Vec<usize>],
    kernel: &KernelType,
    h: f64,
    rho0: f64,
    dt: f64,
) -> f64 {
    let total: f64 = (0..particles.len())
        .map(|i| {
            let rho_star =
                predicted_density_at(particles, velocities, neighbor_lists, i, kernel, h, dt);
            (rho_star - rho0).max(0.0) / rho0
        })
        .sum();
    total / particles.len() as f64
}

/// Mean compressive divergence `mean_i max(Dρ_i/Dt / ρ_i, 0)` (1/s).
fn mean_compressive_divergence(
    particles: &[DfsphParticle],
    velocities: &[[f64; 3]],
    neighbor_lists: &[Vec<usize>],
    kernel: &KernelType,
    h: f64,
) -> f64 {
    let total: f64 = (0..particles.len())
        .map(|i| {
            let rate = density_rate_at(particles, velocities, neighbor_lists, i, kernel, h);
            (rate / particles[i].density).max(0.0)
        })
        .sum();
    total / particles.len() as f64
}

/// Two equal particles 0.06 apart (inside the support `2h = 0.1` of `h = 0.05`)
/// approaching each other at 1 m/s, with the given density.
fn approaching_pair(
    mass: f64,
    density: f64,
) -> (Vec<DfsphParticle>, Vec<Vec<usize>>, Vec<[f64; 3]>) {
    let mut a = DfsphParticle::new([0.0, 0.0, 0.0], mass, density);
    let mut b = DfsphParticle::new([0.06, 0.0, 0.0], mass, density);
    a.velocity = [0.5, 0.0, 0.0];
    b.velocity = [-0.5, 0.0, 0.0];
    let velocities = vec![a.velocity, b.velocity];
    (vec![a, b], vec![vec![1], vec![0]], velocities)
}

#[test]
fn test_density_rate_matches_derivative_of_density_summation() {
    // d/dt Σ_j m_j W(|x_i − x_j|) = Σ_j m_j (v_i − v_j)·∇W_ij: checked against
    // a central finite difference of the summed density along the velocities.
    let kernel = KernelType::WendlandC2;
    let (h, spacing, mass, rho0) = (0.1, 0.05, 0.16, 1000.0);
    let (particles, neighbor_lists, _) =
        build_compressing_cluster(4, spacing, h, mass, rho0, 0.05, &kernel);
    let velocities: Vec<[f64; 3]> = particles.iter().map(|p| p.velocity).collect();
    let density_obj = WcSphDensity::new(kernel.clone(), h);
    let summed_density = |shift: f64, i: usize| -> f64 {
        let moved = |k: usize| add3(particles[k].position, scale3(velocities[k], shift));
        let mut rho = mass * density_obj.eval_kernel(0.0);
        for &j in &neighbor_lists[i] {
            rho += mass * density_obj.eval_kernel(norm3(sub3(moved(i), moved(j))));
        }
        rho
    };
    let eps = 1e-6;
    let mut largest = 0.0f64;
    for i in 0..particles.len() {
        let numeric = (summed_density(eps, i) - summed_density(-eps, i)) / (2.0 * eps);
        let analytic = density_rate_at(&particles, &velocities, &neighbor_lists, i, &kernel, h);
        largest = largest.max(analytic.abs());
        assert!(
            (numeric - analytic).abs() <= 1e-4 * analytic.abs().max(1.0),
            "particle {i}: numeric {numeric}, analytic {analytic}"
        );
    }
    assert!(largest > 10.0, "the cluster must actually be compressing");
}

#[test]
fn test_pressure_solve_iteration_scales_pair_error_by_one_minus_two_omega() {
    // For two equal particles the response of ρ* to its own pressure (the
    // diagonal Newton step, α) and to the partner's pressure are equal, so one
    // iteration with relaxation ω multiplies the compressive error by exactly
    // 1 − 2ω. The wrong (divide-by-α) scaling cannot satisfy this.
    let kernel = KernelType::WendlandC2;
    let (h, rho0, dt) = (0.05, 1000.0, 1e-3);
    for omega in [0.1, 0.25, 0.4] {
        let (particles, neighbor_lists, mut velocities) = approaching_pair(0.3, 1100.0);
        let alphas = AlphaCompute::new(h, kernel.clone()).compute_all(&particles, &neighbor_lists);
        let solver = PressureSolveIter::new(rho0, 1, 1e-3, omega, kernel.clone(), h);
        let mut pressures = vec![0.0; 2];
        let before =
            predicted_density_at(&particles, &velocities, &neighbor_lists, 0, &kernel, h, dt);
        let excess = before - rho0;
        assert!(
            excess > 50.0,
            "the pair must be compressed, excess {excess}"
        );

        let mean = solver.iterate(
            &particles,
            &alphas,
            &mut pressures,
            &mut velocities,
            dt,
            &neighbor_lists,
        );
        assert!(
            (mean - excess / rho0).abs() < 1e-12,
            "returned mean {mean}, expected {}",
            excess / rho0
        );
        let kappa = omega * excess * 1100.0 * alphas[0] / (dt * dt);
        for (i, pressure) in pressures.iter().enumerate() {
            assert!(
                (pressure - kappa).abs() < 1e-9 * kappa,
                "stiffness of particle {i}: {pressure}, expected {kappa}"
            );
            let after =
                predicted_density_at(&particles, &velocities, &neighbor_lists, i, &kernel, h, dt);
            let expected = (1.0 - 2.0 * omega) * excess;
            assert!(
                ((after - rho0) - expected).abs() < 1e-9 * excess,
                "omega {omega}, particle {i}: error {} after one iteration, expected {expected}",
                after - rho0
            );
        }
    }
}

#[test]
fn test_divergence_solve_iteration_scales_pair_divergence_by_one_minus_two_omega() {
    // Same two-particle identity for the divergence solve: the compressive
    // divergence is multiplied by exactly 1 − 2ω per iteration.
    let kernel = KernelType::WendlandC2;
    let (h, dt) = (0.05, 1e-3);
    for omega in [0.1, 0.25, 0.4] {
        let (particles, neighbor_lists, mut velocities) = approaching_pair(0.3, 1100.0);
        let alphas = AlphaCompute::new(h, kernel.clone()).compute_all(&particles, &neighbor_lists);
        let solver = DivergenceSolveIter::new(1, 1e-3, omega, kernel.clone(), h);
        let rate = density_rate_at(&particles, &velocities, &neighbor_lists, 0, &kernel, h);
        let divergence = rate / 1100.0;
        assert!(
            divergence > 1.0,
            "the pair must be compressed: {divergence}"
        );

        let mean = solver.iterate(&particles, &alphas, &mut velocities, dt, &neighbor_lists);
        assert!(
            (mean - divergence).abs() < 1e-12 * divergence,
            "returned mean {mean}, expected {divergence}"
        );
        for i in 0..2 {
            let after =
                density_rate_at(&particles, &velocities, &neighbor_lists, i, &kernel, h) / 1100.0;
            let expected = (1.0 - 2.0 * omega) * divergence;
            assert!(
                (after - expected).abs() < 1e-9 * divergence,
                "omega {omega}, particle {i}: divergence {after} after one iteration, expected {expected}"
            );
        }
    }
}

#[test]
fn test_pressure_solve_removes_compression_in_cluster_at_default_omega() {
    // A strongly compressed, inward-moving cluster: with the default relaxation
    // ω = 0.5 the first iteration drives every compressive predicted-density
    // error to zero, pushes the compressed particles outward and stays finite.
    let kernel = KernelType::WendlandC2;
    let (h, spacing, mass, rho0, v_in, dt) = (0.1, 0.05, 0.16, 1000.0, 0.05, 1e-3);
    let (particles, neighbor_lists, center_idx) =
        build_compressing_cluster(5, spacing, h, mass, rho0, v_in, &kernel);
    assert!(
        (particles[center_idx].density - rho0) / rho0 > 0.1,
        "the centre must be over-dense"
    );
    let alphas = AlphaCompute::new(h, kernel.clone()).compute_all(&particles, &neighbor_lists);
    let solver = PressureSolveIter::new(rho0, 100, 1e-3, 0.5, kernel.clone(), h);
    let mut velocities: Vec<[f64; 3]> = particles.iter().map(|p| p.velocity).collect();
    let mut pressures = vec![0.0f64; particles.len()];

    // Particle (2, 2, 1) sits one layer below the centre and moves inward (+z).
    let probe = (2 * 5 + 2) * 5 + 1;
    let inward_before = velocities[probe][2];
    assert!(inward_before > 0.0);

    let initial = mean_compressive_density_error(
        &particles,
        &velocities,
        &neighbor_lists,
        &kernel,
        h,
        rho0,
        dt,
    );
    assert!(
        initial > 0.01,
        "the cluster must start compressed: {initial}"
    );
    let reported = solver.iterate(
        &particles,
        &alphas,
        &mut pressures,
        &mut velocities,
        dt,
        &neighbor_lists,
    );
    assert!(
        (reported - initial).abs() < 1e-9 * initial,
        "iterate reports the pre-projection error: {reported} vs {initial}"
    );
    assert!(
        pressures[center_idx] > 0.0,
        "compressed centre needs pressure"
    );
    assert!(pressures.iter().all(|&p| p >= 0.0 && p.is_finite()));
    assert!(
        velocities[probe][2] < inward_before,
        "the projection must oppose the inward motion: {} -> {}",
        inward_before,
        velocities[probe][2]
    );
    assert!(velocities.iter().all(|v| v.iter().all(|c| c.is_finite())));

    let remaining = mean_compressive_density_error(
        &particles,
        &velocities,
        &neighbor_lists,
        &kernel,
        h,
        rho0,
        dt,
    );
    assert!(
        remaining < 1e-3 * initial,
        "compression must be removed in one pass: {initial:e} -> {remaining:e}"
    );
    // A second pass finds nothing left to correct.
    let second = solver.iterate(
        &particles,
        &alphas,
        &mut pressures,
        &mut velocities,
        dt,
        &neighbor_lists,
    );
    assert!(second < 1e-9, "second pass error {second:e}");
}

#[test]
fn test_pressure_solve_enforces_incompressibility_end_to_end() {
    // Strongly over-dense (≈ 25 % at the centre), compressing cluster. With the
    // Newton-scaled stiffness the iteration is a contraction for a modest
    // under-relaxation (here ω = 0.1; the single-step gain at the centre is
    // ≈ 4.6 ω), so the predicted-density residual at the centre converges to
    // the rest density. Larger ω overshoots into expansion in one step; the
    // one-sided correction then stops, which is covered by
    // `test_pressure_solve_removes_compression_in_cluster_at_default_omega`.
    let kernel = KernelType::WendlandC2;
    let (h, spacing, mass, rho0, v_in, dt) = (0.1, 0.05, 0.16, 1000.0, 0.05, 1e-3);
    let (particles, neighbor_lists, center_idx) =
        build_compressing_cluster(5, spacing, h, mass, rho0, v_in, &kernel);
    let alphas = AlphaCompute::new(h, kernel.clone()).compute_all(&particles, &neighbor_lists);
    let solver = PressureSolveIter::new(rho0, 200, 1e-3, 0.1, kernel.clone(), h);
    let mut velocities: Vec<[f64; 3]> = particles.iter().map(|p| p.velocity).collect();
    let mut pressures = vec![0.0f64; particles.len()];

    let residual = |velocities: &[[f64; 3]]| {
        let rho_star = predicted_density_at(
            &particles,
            velocities,
            &neighbor_lists,
            center_idx,
            &kernel,
            h,
            dt,
        );
        ((rho_star - rho0) / rho0).abs()
    };
    let initial_residual = residual(&velocities);
    assert!(
        initial_residual > 0.2,
        "initial centre residual should be large, got {initial_residual:.3e}"
    );

    let mut last_mean = f64::INFINITY;
    let mut iterations = 0;
    for it in 0..solver.max_iter {
        last_mean = solver.iterate(
            &particles,
            &alphas,
            &mut pressures,
            &mut velocities,
            dt,
            &neighbor_lists,
        );
        iterations = it + 1;
        assert!(last_mean.is_finite(), "iteration {it} went non-finite");
        let current = residual(&velocities);
        assert!(
            current <= initial_residual,
            "iteration {it}: centre residual grew to {current:.3e}"
        );
        if last_mean < 1e-6 {
            break;
        }
    }
    let final_residual = residual(&velocities);
    assert!(
        last_mean < 1e-6,
        "mean compressive error {last_mean:e} after {iterations} iterations"
    );
    assert!(
        final_residual < 1e-3,
        "constant-density solve must drive the centre residual below 1e-3; \
         initial = {initial_residual:.3e}, final = {final_residual:.3e}"
    );
}

#[test]
fn test_divergence_solve_removes_compressive_divergence_in_cluster() {
    // The compressing cluster has a large compressive divergence; the
    // compression-only divergence solve contracts it geometrically at the
    // default relaxation ω = 0.5 and never blows up (the two-sided variant
    // already diverges at ω = 0.3 on this compact cluster).
    let kernel = KernelType::WendlandC2;
    let (h, spacing, mass, rho0, v_in, dt) = (0.1, 0.05, 0.16, 1000.0, 0.05, 1e-3);
    let (particles, neighbor_lists, _) =
        build_compressing_cluster(5, spacing, h, mass, rho0, v_in, &kernel);
    let alphas = AlphaCompute::new(h, kernel.clone()).compute_all(&particles, &neighbor_lists);
    let solver = DivergenceSolveIter::new(100, 1e-3, 0.5, kernel.clone(), h);
    let mut velocities: Vec<[f64; 3]> = particles.iter().map(|p| p.velocity).collect();

    let initial = mean_compressive_divergence(&particles, &velocities, &neighbor_lists, &kernel, h);
    assert!(initial > 0.5, "initial compressive divergence {initial}");
    let mut previous = f64::INFINITY;
    for it in 0..30 {
        let mean = solver.iterate(&particles, &alphas, &mut velocities, dt, &neighbor_lists);
        assert!(mean.is_finite(), "iteration {it} went non-finite");
        assert!(
            mean < previous,
            "iteration {it}: divergence {mean:e} did not decrease from {previous:e}"
        );
        previous = mean;
    }
    assert!(
        previous < 1e-6,
        "compressive divergence {previous:e} after 30 iterations (initial {initial:e})"
    );
    let measured =
        mean_compressive_divergence(&particles, &velocities, &neighbor_lists, &kernel, h);
    assert!(measured < 1e-6, "independent measurement {measured:e}");
    let fastest = velocities.iter().map(|v| norm3(*v)).fold(0.0f64, f64::max);
    assert!(
        fastest.is_finite() && fastest < 1.0,
        "corrected velocities must stay bounded, got {fastest}"
    );
}

#[test]
fn test_velocity_advection_gravity() {
    let adv = VelocityAdvection::new([0.0, -9.81, 0.0], 1e-6, 0.1, KernelType::WendlandC2);
    let mut particles = vec![DfsphParticle::new([0.0; 3], 1.0, 1000.0)];
    let neighbor_data: Vec<Vec<NeighborEntry>> = vec![vec![]];
    let dt = 0.01;
    adv.advect_velocity(&mut particles, &neighbor_data, dt);
    // After dt = 0.01 with g = -9.81: v_y should be -0.0981
    assert!((particles[0].velocity[1] - (-9.81 * dt)).abs() < 1e-10);
}

#[test]
fn test_dfsph_boundary_floor() {
    let mut boundary = DfsphBoundary::new(BoundaryType::NoSlip, 1000.0, 0.1);
    boundary.create_floor(0.0, 1.0, 0.0, 0.1);
    assert!(boundary.num_particles() > 0, "floor should have particles");
}

#[test]
fn test_dfsph_boundary_density_contribution() {
    let mut boundary = DfsphBoundary::new(BoundaryType::NoSlip, 1000.0, 0.1);
    boundary.add_particle([0.0, 0.0, 0.0], 0.001);
    let pos_i = [0.05, 0.0, 0.0];
    let contrib = boundary.density_contribution(pos_i, &KernelType::WendlandC2);
    assert!(
        contrib > 0.0,
        "boundary density contribution should be positive"
    );
}

#[test]
fn test_dfsph_solver_step() {
    let cfg = DfsphConfig {
        rest_density: 1000.0,
        smoothing_length: 0.1,
        support_radius: 0.2,
        max_pressure_iterations: 5,
        max_divergence_iterations: 5,
        error_threshold: 1e-2,
        omega: 0.5,
        viscosity: 1e-6,
        gravity: [0.0, -9.81, 0.0],
        surface_tension: 0.0,
    };
    let mut solver = DfsphSolver::new(cfg);
    // Add a small cluster of particles
    for i in 0..3 {
        for j in 0..3 {
            solver.add_particle([i as f64 * 0.05, j as f64 * 0.05, 0.0], [0.0; 3], 0.001);
        }
    }
    let (dens_err, div_err) = solver.step(0.001);
    assert!(dens_err >= 0.0, "density error should be non-negative");
    assert!(div_err >= 0.0, "divergence error should be non-negative");
}

#[test]
fn test_dfsph_solver_mean_density() {
    let cfg = make_default_config();
    let mut solver = DfsphSolver::new(cfg);
    solver.add_particle([0.0; 3], [0.0; 3], 0.001);
    // Before update, density equals rest_density
    assert!((solver.mean_density() - 1000.0).abs() < 1e-10);
}

// ---------------------------------------------------------------------------
// Kernel normalisation (Wendland C2 / C4)
// ---------------------------------------------------------------------------

/// Midpoint rule for `∫_0^{2h} 4π r² W(r, h) dr`.
fn radial_volume_integral(kernel: fn(f64, f64) -> f64, h: f64) -> f64 {
    let steps = 20_000;
    let dr = 2.0 * h / steps as f64;
    (0..steps)
        .map(|k| {
            let r = (k as f64 + 0.5) * dr;
            4.0 * std::f64::consts::PI * r * r * kernel(r, h) * dr
        })
        .sum()
}

#[test]
fn test_wendland_kernels_integrate_to_one_over_the_support_sphere() {
    // An 8x over-normalised kernel integrates to 8.
    for h in [0.01, 0.1, 1.0] {
        let c2 = radial_volume_integral(wendland_c2_kernel, h);
        assert!((c2 - 1.0).abs() < 1e-6, "∫W_C2 dV = {c2} for h = {h}");
        let c4 = radial_volume_integral(wendland_c4_kernel, h);
        assert!((c4 - 1.0).abs() < 1e-6, "∫W_C4 dV = {c4} for h = {h}");
    }
}

#[test]
fn test_wendland_gradient_magnitude_is_derivative_of_kernel() {
    type RadialFn = fn(f64, f64) -> f64;
    let pairs: [(RadialFn, RadialFn); 2] = [
        (wendland_c2_kernel, wendland_c2_gradient_mag),
        (wendland_c4_kernel, wendland_c4_gradient_mag),
    ];
    for h in [0.05, 1.0] {
        for (kernel, gradient) in pairs {
            for frac in [0.1, 0.5, 1.0, 1.5, 1.9] {
                let r = frac * h;
                let d = 1e-6 * h;
                let numeric = (kernel(r + d, h) - kernel(r - d, h)) / (2.0 * d);
                let analytic = gradient(r, h);
                assert!(analytic < 0.0, "kernel must decrease at r = {frac} h");
                assert!(
                    (numeric - analytic).abs() < 1e-6 * analytic.abs(),
                    "dW/dr at r = {frac} h, h = {h}: numeric {numeric}, analytic {analytic}"
                );
            }
        }
    }
}

#[test]
fn test_wendland_gradient_vector_points_toward_the_neighbour() {
    // ∇_i W(x_i − x_j) = (dW/dr) r̂ with dW/dr < 0, so it points from i to j.
    let h = 0.1;
    let r_vec = [0.03, -0.04, 0.0];
    let r = norm3(r_vec);
    for gradient in [wendland_c2_gradient, wendland_c4_gradient] {
        let g = gradient(r_vec, h);
        assert!(dot3(g, r_vec) < 0.0, "gradient must oppose x_i - x_j");
        let cross = [
            g[1] * r_vec[2] - g[2] * r_vec[1],
            g[2] * r_vec[0] - g[0] * r_vec[2],
            g[0] * r_vec[1] - g[1] * r_vec[0],
        ];
        assert!(
            norm3(cross) < 1e-9 * norm3(g) * r,
            "gradient not parallel to r"
        );
    }
    assert_eq!(wendland_c2_gradient([0.0; 3], h), [0.0; 3]);
    assert_eq!(wendland_c4_gradient([0.0; 3], h), [0.0; 3]);
}

#[test]
fn test_wendland_density_of_a_lattice_reproduces_the_rest_density() {
    // Summation density of a cubic lattice with spacing dx and m = ρ₀ dx³ is ρ₀
    // (to the lattice-sampling error) for a correctly normalised kernel.
    let (spacing, rho0) = (0.05, 1000.0);
    let h = 1.5 * spacing;
    let mass = rho0 * spacing * spacing * spacing;
    let n_side = 11i32;
    let mid = n_side / 2;
    let mut neighbours = Vec::new();
    for ix in 0..n_side {
        for iy in 0..n_side {
            for iz in 0..n_side {
                let offset = [
                    (ix - mid) as f64 * spacing,
                    (iy - mid) as f64 * spacing,
                    (iz - mid) as f64 * spacing,
                ];
                neighbours.push((offset, mass));
            }
        }
    }
    for kernel in [KernelType::WendlandC2, KernelType::WendlandC4] {
        let density = WcSphDensity::new(kernel.clone(), h).compute_density([0.0; 3], &neighbours);
        assert!(
            (density - rho0).abs() < 0.01 * rho0,
            "{kernel:?}: lattice density {density}, expected about {rho0}"
        );
    }
}

// ---------------------------------------------------------------------------
// DfsphSolver: constant-density and divergence-free solves
// ---------------------------------------------------------------------------

/// Solver holding the 5×5×5 compressing cluster of [`build_compressing_cluster`]
/// (spacing 0.05, mass 0.16, inward speed 0.05 m/s) with neighbour lists,
/// densities and alphas prepared.
fn prepared_cluster_solver(config: DfsphConfig) -> DfsphSolver {
    let (n_side, spacing, mass, v_in) = (5usize, 0.05, 0.16, 0.05);
    let mut solver = DfsphSolver::new(config);
    let c = (n_side as f64 - 1.0) * spacing / 2.0;
    for ix in 0..n_side {
        for iy in 0..n_side {
            for iz in 0..n_side {
                let pos = [
                    ix as f64 * spacing,
                    iy as f64 * spacing,
                    iz as f64 * spacing,
                ];
                let to_center = sub3([c, c, c], pos);
                let r = norm3(to_center);
                let vel = if r > 1e-9 {
                    scale3(to_center, v_in / r)
                } else {
                    [0.0; 3]
                };
                solver.add_particle(pos, vel, mass);
            }
        }
    }
    solver.build_spatial_hash();
    solver.build_neighbor_lists();
    solver.update_densities();
    solver.update_alphas();
    solver
}

/// Current particle velocities of `solver`.
fn solver_velocities(solver: &DfsphSolver) -> Vec<[f64; 3]> {
    solver.particles.iter().map(|p| p.velocity).collect()
}

/// Independent mean compressive density error of the solver's current state.
fn solver_density_error(solver: &DfsphSolver, dt: f64) -> f64 {
    mean_compressive_density_error(
        &solver.particles,
        &solver_velocities(solver),
        &solver.neighbor_lists,
        &KernelType::WendlandC2,
        solver.config.smoothing_length,
        solver.config.rest_density,
        dt,
    )
}

/// Independent mean compressive divergence of the solver's current state.
fn solver_divergence_error(solver: &DfsphSolver) -> f64 {
    mean_compressive_divergence(
        &solver.particles,
        &solver_velocities(solver),
        &solver.neighbor_lists,
        &KernelType::WendlandC2,
        solver.config.smoothing_length,
    )
}

/// Two particles at rest-distance 0.06 with kernel-summed density above ρ₀,
/// approaching each other, prepared for a solve with relaxation `omega`.
fn prepared_pair_solver(omega: f64, max_iterations: usize, threshold: f64) -> DfsphSolver {
    let config = DfsphConfig {
        rest_density: 1000.0,
        smoothing_length: 0.05,
        support_radius: 0.1,
        max_pressure_iterations: max_iterations,
        max_divergence_iterations: max_iterations,
        error_threshold: threshold,
        omega,
        viscosity: 0.0,
        gravity: [0.0; 3],
        surface_tension: 0.0,
    };
    let mut solver = DfsphSolver::new(config);
    solver.add_particle([0.0; 3], [0.5, 0.0, 0.0], 0.3);
    solver.add_particle([0.06, 0.0, 0.0], [-0.5, 0.0, 0.0], 0.3);
    solver.build_spatial_hash();
    solver.build_neighbor_lists();
    solver.update_densities();
    solver.update_alphas();
    solver
}

#[test]
fn test_solver_constant_density_pair_error_follows_one_minus_two_omega() {
    let dt = 1e-3;
    let omega = 0.25;
    let mut solver = prepared_pair_solver(omega, 1, 1e-12);
    let rho0 = solver.config.rest_density;
    assert!(
        solver.particles[0].density > 1.05 * rho0,
        "pair must be over-dense, got {}",
        solver.particles[0].density
    );
    let e0 = solver_density_error(&solver, dt);
    assert!(e0 > 0.05, "initial error {e0}");
    let err = solver.solve_constant_density(dt);
    // One projection scales the error by (1 − 2ω); it is recorded although the
    // tolerance 1e-12 was not met.
    assert!(
        (err - (1.0 - 2.0 * omega) * e0).abs() < 1e-9 * e0,
        "error after one iteration {err}, expected {}",
        (1.0 - 2.0 * omega) * e0
    );
    assert!((solver_density_error(&solver, dt) - err).abs() < 1e-12);
    assert_eq!(solver.stats.pressure_iterations, vec![1]);
    assert_eq!(solver.stats.density_error_history, vec![err]);
}

#[test]
fn test_solver_constant_density_pair_converges_geometrically() {
    // Error after k iterations is e0 (1 − 2ω)^k; the solve stops at the first
    // iteration whose pre-projection error is below the threshold.
    let dt = 1e-3;
    let omega = 0.25;
    let threshold = 1e-6;
    let mut solver = prepared_pair_solver(omega, 200, threshold);
    let e0 = solver_density_error(&solver, dt);
    let expected_iterations = (0..200)
        .find(|&k| e0 * (1.0 - 2.0 * omega).powi(k) < threshold)
        .map(|k| k as usize + 1)
        .expect("the geometric sequence must reach the threshold");
    let err = solver.solve_constant_density(dt);
    assert!(err < threshold, "final error {err:e}");
    assert_eq!(solver.stats.pressure_iterations, vec![expected_iterations]);
    assert!(expected_iterations > 10 && expected_iterations < 200);
    assert!((solver.stats.density_error_history[0] - err).abs() < 1e-18);
}

#[test]
fn test_solver_divergence_free_pair_error_follows_one_minus_two_omega() {
    let dt = 1e-3;
    let omega = 0.25;
    let mut solver = prepared_pair_solver(omega, 1, 1e-12);
    let d0 = solver_divergence_error(&solver);
    assert!(d0 > 1.0, "initial compressive divergence {d0}");
    let err = solver.solve_divergence_free(dt);
    assert!(
        (err - (1.0 - 2.0 * omega) * d0).abs() < 1e-9 * d0,
        "divergence after one iteration {err}, expected {}",
        (1.0 - 2.0 * omega) * d0
    );
    assert!((solver_divergence_error(&solver) - err).abs() < 1e-12);
    assert_eq!(solver.stats.divergence_iterations, vec![1]);
    assert_eq!(solver.stats.divergence_error_history, vec![err]);
    // The divergence field holds the signed Dρ/Dt / ρ of the final state, which
    // is still compressive for ω < 0.5.
    let expected_field = (1.0 - 2.0 * omega) * d0;
    for p in &solver.particles {
        assert!(
            (p.divergence - expected_field).abs() < 1e-9 * d0,
            "divergence field {}, expected {expected_field}",
            p.divergence
        );
    }
}

#[test]
fn test_solver_constant_density_cluster_converges_and_keeps_state_intact() {
    let dt = 1e-3;
    let mut solver = prepared_cluster_solver(make_default_config());
    let lists_before = solver.neighbor_lists.clone();
    let positions_before: Vec<[f64; 3]> = solver.particles.iter().map(|p| p.position).collect();
    let densities_before: Vec<f64> = solver.particles.iter().map(|p| p.density).collect();

    let initial = solver_density_error(&solver, dt);
    assert!(initial > 0.01, "cluster must start compressed: {initial}");
    let err = solver.solve_constant_density(dt);
    let threshold = solver.config.error_threshold;
    assert!(err < threshold, "error {err:e} not below {threshold:e}");
    assert!(solver_density_error(&solver, dt) < threshold);

    // Statistics: exactly one entry, below the iteration cap (the frozen-density
    // implementation could never converge and recorded nothing).
    assert_eq!(solver.stats.pressure_iterations.len(), 1);
    let iterations = solver.stats.pressure_iterations[0];
    assert!(
        (1..solver.config.max_pressure_iterations).contains(&iterations),
        "iterations {iterations}"
    );
    assert_eq!(solver.stats.density_error_history, vec![err]);

    // The solve only edits velocities and leaves neighbour lists untouched.
    assert_eq!(solver.neighbor_lists, lists_before);
    for (p, x0) in solver.particles.iter().zip(positions_before.iter()) {
        assert_eq!(p.position, *x0);
    }
    for (p, rho0) in solver.particles.iter().zip(densities_before.iter()) {
        assert_eq!(p.density, *rho0);
    }
    assert!(
        solver
            .particles
            .iter()
            .all(|p| p.velocity.iter().all(|c| c.is_finite()))
    );
    let fastest = solver
        .particles
        .iter()
        .map(|p| norm3(p.velocity))
        .fold(0.0f64, f64::max);
    assert!(fastest < 100.0, "velocities stay bounded, got {fastest}");
}

#[test]
fn test_solver_divergence_free_cluster_converges_and_keeps_state_intact() {
    let dt = 1e-3;
    let mut solver = prepared_cluster_solver(make_default_config());
    let lists_before = solver.neighbor_lists.clone();
    let initial = solver_divergence_error(&solver);
    assert!(initial > 0.5, "initial compressive divergence {initial}");
    let err = solver.solve_divergence_free(dt);
    let threshold = solver.config.error_threshold;
    assert!(err < threshold, "error {err:e} not below {threshold:e}");
    assert!(solver_divergence_error(&solver) < threshold);
    assert_eq!(solver.stats.divergence_iterations.len(), 1);
    let iterations = solver.stats.divergence_iterations[0];
    assert!(
        (2..solver.config.max_divergence_iterations).contains(&iterations),
        "iterations {iterations}"
    );
    assert_eq!(solver.stats.divergence_error_history, vec![err]);
    assert_eq!(solver.neighbor_lists, lists_before);
    let fastest = solver
        .particles
        .iter()
        .map(|p| norm3(p.velocity))
        .fold(0.0f64, f64::max);
    assert!(fastest.is_finite() && fastest < 1.0, "max speed {fastest}");
}

#[test]
fn test_solver_divergence_free_leaves_a_divergence_free_field_untouched() {
    // Particles at rest have Dρ/Dt = 0: nothing to correct, converged at once.
    let dt = 1e-3;
    let mut solver = prepared_cluster_solver(make_default_config());
    for p in solver.particles.iter_mut() {
        p.velocity = [0.0; 3];
    }
    let err = solver.solve_divergence_free(dt);
    assert_eq!(err, 0.0);
    assert_eq!(solver.stats.divergence_iterations, vec![1]);
    assert!(solver.particles.iter().all(|p| p.velocity == [0.0; 3]));
}

#[test]
fn test_solver_records_statistics_when_not_converged() {
    // With a tiny relaxation and a tight tolerance neither solve converges in
    // three iterations, but both must still record the iterations used and the
    // residual they stopped at.
    let dt = 1e-3;
    let config = DfsphConfig {
        omega: 0.05,
        max_pressure_iterations: 3,
        max_divergence_iterations: 3,
        error_threshold: 1e-12,
        ..make_default_config()
    };
    let mut density_solver = prepared_cluster_solver(config.clone());
    let density_before = solver_density_error(&density_solver, dt);
    let density_err = density_solver.solve_constant_density(dt);
    assert!(
        density_err > 1e-12,
        "cannot have converged: {density_err:e}"
    );
    assert!(
        density_err < density_before,
        "{density_before} -> {density_err}"
    );
    assert_eq!(density_solver.stats.pressure_iterations, vec![3]);
    assert_eq!(
        density_solver.stats.density_error_history,
        vec![density_err]
    );
    assert_eq!(density_solver.stats.mean_density_error(), density_err);

    let mut divergence_solver = prepared_cluster_solver(config);
    let divergence_before = solver_divergence_error(&divergence_solver);
    let divergence_err = divergence_solver.solve_divergence_free(dt);
    assert!(divergence_err > 1e-12, "cannot have converged");
    assert!(
        divergence_err < divergence_before,
        "{divergence_before} -> {divergence_err}"
    );
    assert_eq!(divergence_solver.stats.divergence_iterations, vec![3]);
    assert_eq!(
        divergence_solver.stats.divergence_error_history,
        vec![divergence_err]
    );
}

#[test]
fn test_solver_step_on_compressed_cluster_stays_bounded() {
    let config = DfsphConfig {
        gravity: [0.0; 3],
        ..make_default_config()
    };
    let mut solver = prepared_cluster_solver(config);
    let dt = 1e-3;
    for _ in 0..5 {
        let (density_error, divergence_error) = solver.step(dt);
        assert!(density_error.is_finite() && density_error >= 0.0);
        assert!(divergence_error.is_finite() && divergence_error >= 0.0);
    }
    assert_eq!(solver.stats.pressure_iterations.len(), 5);
    assert_eq!(solver.stats.divergence_iterations.len(), 5);
    for p in &solver.particles {
        assert!(p.position.iter().all(|c| c.is_finite()));
        assert!(p.velocity.iter().all(|c| c.is_finite()));
        assert!(p.density.is_finite() && p.density > 0.0);
    }
    assert!(solver.mean_density() < 2000.0, "{}", solver.mean_density());
}

// ---------------------------------------------------------------------------
// Neighbour lists and scaling
// ---------------------------------------------------------------------------

/// Deterministic pseudo-random numbers in `[0, 1)` (64-bit LCG).
fn lcg_unit(state: &mut u64) -> f64 {
    *state = state
        .wrapping_mul(6364136223846793005)
        .wrapping_add(1442695040888963407);
    ((*state >> 11) as f64) / ((1u64 << 53) as f64)
}

#[test]
fn test_build_neighbor_lists_matches_brute_force_on_random_cloud() {
    // Includes negative coordinates (floor-based cell keys) and distances very
    // close to the support radius; each list contains the particle itself.
    let mut state = 0x9E37_79B9_7F4A_7C15_u64;
    let config = make_default_config();
    let support = config.support_radius;
    let mut solver = DfsphSolver::new(config);
    let mut positions = Vec::new();
    for _ in 0..250 {
        let pos = [
            0.8 * lcg_unit(&mut state) - 0.4,
            0.8 * lcg_unit(&mut state) - 0.4,
            0.8 * lcg_unit(&mut state) - 0.4,
        ];
        positions.push(pos);
        solver.add_particle(pos, [0.0; 3], 0.01);
    }
    solver.build_spatial_hash();
    solver.build_neighbor_lists();
    assert_eq!(solver.neighbor_lists.len(), positions.len());
    let mut total = 0usize;
    for (i, list) in solver.neighbor_lists.iter().enumerate() {
        let expected: Vec<usize> = (0..positions.len())
            .filter(|&j| norm_sq3(sub3(positions[i], positions[j])) < support * support)
            .collect();
        let mut got = list.clone();
        got.sort_unstable();
        assert_eq!(got, expected, "neighbour list of particle {i}");
        assert_eq!(got.iter().filter(|&&j| j == i).count(), 1);
        total += got.len();
    }
    assert!(total > 2 * positions.len(), "cloud too sparse: {total}");
}

#[test]
fn test_build_neighbor_lists_tolerates_escaped_particles() {
    // Particles with infinite, NaN or astronomically large coordinates (a
    // diverged simulation) must neither overflow the cell arithmetic nor
    // acquire neighbours, and must leave the healthy lists exact.
    let mut state = 0xA076_1D64_78BD_642F_u64;
    let config = make_default_config();
    let support = config.support_radius;
    let mut solver = DfsphSolver::new(config);
    let mut healthy = Vec::new();
    for _ in 0..60 {
        let pos = [
            0.4 * lcg_unit(&mut state),
            0.4 * lcg_unit(&mut state),
            0.4 * lcg_unit(&mut state),
        ];
        healthy.push(pos);
        solver.add_particle(pos, [0.0; 3], 0.01);
    }
    for pos in [
        [f64::INFINITY, 0.0, 0.0],
        [0.0, f64::NEG_INFINITY, 0.0],
        [f64::NAN, 0.0, 0.0],
        [1.0e300, -1.0e300, 1.0e300],
    ] {
        solver.add_particle(pos, [0.0; 3], 0.01);
    }
    solver.build_spatial_hash();
    solver.build_neighbor_lists();
    assert_eq!(solver.neighbor_lists.len(), healthy.len() + 4);
    // Non-finite particles have no entries at all; the finite but astronomically
    // distant one only lists itself (distance zero).
    let escaped = &solver.neighbor_lists[healthy.len()..];
    for list in &escaped[..3] {
        assert!(
            list.is_empty(),
            "non-finite particle has neighbours {list:?}"
        );
    }
    assert_eq!(escaped[3], vec![healthy.len() + 3]);
    for (i, list) in solver.neighbor_lists.iter().enumerate().take(healthy.len()) {
        let expected: Vec<usize> = (0..healthy.len())
            .filter(|&j| norm_sq3(sub3(healthy[i], healthy[j])) < support * support)
            .collect();
        let mut got = list.clone();
        got.sort_unstable();
        assert_eq!(got, expected, "neighbour list of healthy particle {i}");
    }
}

#[test]
fn test_solvers_scale_to_thousands_of_particles() {
    // 12³ = 1728 particles, ~65 neighbours each. Cloning the whole neighbour
    // table for every particle and iteration (O(N²·k) per iteration) made this
    // take minutes; the Jacobi scratch buffer makes it milliseconds.
    let spacing = 0.05;
    let config = DfsphConfig {
        smoothing_length: 0.065,
        support_radius: 0.13,
        max_pressure_iterations: 20,
        max_divergence_iterations: 20,
        error_threshold: 1e-6,
        gravity: [0.0, -9.81, 0.0],
        ..make_default_config()
    };
    let mut solver = DfsphSolver::new(config);
    let mass = 1000.0 * spacing * spacing * spacing;
    for ix in 0..12 {
        for iy in 0..12 {
            for iz in 0..12 {
                let pos = [
                    ix as f64 * spacing,
                    iy as f64 * spacing,
                    iz as f64 * spacing,
                ];
                // A smooth compressive velocity pattern so both solves have work.
                let vel = [
                    -0.05 * (ix as f64 - 5.5) / 5.5,
                    -0.05 * (iy as f64 - 5.5) / 5.5,
                    -0.05 * (iz as f64 - 5.5) / 5.5,
                ];
                solver.add_particle(pos, vel, mass);
            }
        }
    }
    let (density_error, divergence_error) = solver.step(1e-3);
    assert!(density_error.is_finite() && divergence_error.is_finite());
    assert_eq!(solver.stats.pressure_iterations.len(), 1);
    assert_eq!(solver.stats.divergence_iterations.len(), 1);
    let average_neighbours = solver
        .neighbor_lists
        .iter()
        .map(|list| list.len() as f64)
        .sum::<f64>()
        / solver.neighbor_lists.len() as f64;
    assert!(average_neighbours > 30.0, "{average_neighbours}");
    assert!(
        solver
            .particles
            .iter()
            .all(|p| { p.density.is_finite() && p.velocity.iter().all(|c| c.is_finite()) })
    );
}
