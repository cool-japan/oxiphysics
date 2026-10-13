// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Unit tests for [`super`] (adaptive smoothing lengths and the grad-h correction).

use super::*;

// -----------------------------------------------------------------------
// Hoover criterion tests
// -----------------------------------------------------------------------

#[test]
fn hoover_h_at_reference_density() {
    // At ρ_i = ρ₀, hoover_h should return h₀.
    let h = hoover_h(0.1, 1000.0, 1000.0, 3);
    assert!(
        (h - 0.1).abs() < 1e-14,
        "At reference density, h should equal h₀: {h}"
    );
}

#[test]
fn hoover_h_denser_fluid_gives_smaller_h() {
    // Higher density → smaller h (need fewer radius to find same number of neighbours).
    let h_dense = hoover_h(0.1, 1000.0, 2000.0, 3);
    let h_ref = 0.1_f64;
    assert!(
        h_dense < h_ref,
        "Denser fluid should have smaller h: {h_dense} vs {h_ref}"
    );
}

#[test]
fn hoover_h_rarefied_fluid_gives_larger_h() {
    // Lower density → larger h.
    let h_rare = hoover_h(0.1, 1000.0, 500.0, 3);
    assert!(
        h_rare > 0.1,
        "Rarefied fluid should have larger h: {h_rare}"
    );
}

#[test]
fn hoover_residual_zero_at_solution() {
    let h0 = 0.1;
    let rho0 = 1000.0;
    let rho_i = 800.0;
    let h_sol = hoover_h(h0, rho0, rho_i, 3);
    let res = hoover_residual(h0, h_sol, rho0, rho_i, 3);
    assert!(
        res.abs() < 1e-10,
        "Residual at Hoover solution should be ~0: {res}"
    );
}

// -----------------------------------------------------------------------
// SPH density sum tests
// -----------------------------------------------------------------------

#[test]
fn sph_density_single_neighbour_at_origin() {
    let pos_i = [0.0_f64; 3];
    let pos_j = [[0.0_f64; 3]];
    let mass_j = [1.0_f64];
    let h = 0.1;
    let kernel = CubicSplineKernel;
    let rho = sph_density(pos_i, h, &pos_j, &mass_j);
    let expected = kernel.w(0.0, h);
    assert!(
        (rho - expected).abs() < 1e-14,
        "Self contribution: {rho} vs {expected}"
    );
}

#[test]
fn sph_density_zero_outside_support() {
    let pos_i = [0.0_f64; 3];
    let pos_j = [[5.0_f64, 0.0, 0.0]]; // r = 5 >> 2h
    let mass_j = [1.0_f64];
    let h = 0.1;
    let rho = sph_density(pos_i, h, &pos_j, &mass_j);
    assert!(
        rho.abs() < 1e-14,
        "Particle outside support should contribute 0: {rho}"
    );
}

#[test]
fn sph_density_increases_with_more_neighbours() {
    let pos_i = [0.0_f64; 3];
    let h = 0.1;
    let pos_1 = vec![[0.05_f64, 0.0, 0.0]];
    let pos_2 = vec![[0.05_f64, 0.0, 0.0], [-0.05, 0.0, 0.0]];
    let mass = vec![1.0_f64];
    let mass2 = vec![1.0_f64, 1.0_f64];

    let rho1 = sph_density(pos_i, h, &pos_1, &mass);
    let rho2 = sph_density(pos_i, h, &pos_2, &mass2);
    assert!(rho2 > rho1, "More neighbours → higher density");
}

// -----------------------------------------------------------------------
// Newton-Raphson solver tests
// -----------------------------------------------------------------------

#[test]
fn solve_h_newton_converges_to_known_density() {
    let pos_i = [0.0_f64; 3];
    let h0 = 0.1;
    let positions = vec![[0.05_f64, 0.0, 0.0], [-0.05, 0.0, 0.0], [0.0, 0.05, 0.0]];
    let masses = vec![0.001_f64; 3];

    // Compute the density at h0 to use as target
    let rho_target = sph_density(pos_i, h0, &positions, &masses);

    let h_solved = solve_h_newton(
        pos_i, h0, rho_target, &positions, &masses, 100, 1e-8, 0.01, 1.0,
    );

    assert!(h_solved.is_some(), "Newton solver should converge");
    let h_s = h_solved.unwrap();
    let rho_check = sph_density(pos_i, h_s, &positions, &masses);
    let rel_err = (rho_check - rho_target).abs() / rho_target.max(1e-14);
    assert!(
        rel_err < 1e-6,
        "Solved h gives wrong density: target={rho_target}, got={rho_check}"
    );
}

#[test]
fn solve_h_newton_clamps_to_bounds() {
    let pos_i = [0.0_f64; 3];
    let h_min = 0.05;
    let h_max = 0.5;
    // Use a very low target density (only possible with very large h)
    let positions = vec![[0.05_f64, 0.0, 0.0]];
    let masses = vec![0.0001_f64];
    let h_solved = solve_h_newton(
        pos_i, 0.1, 1e10, // unreachably large target
        &positions, &masses, 50, 1e-6, h_min, h_max,
    );
    // Should either converge within bounds or return None (but if Some, must be in bounds)
    if let Some(h) = h_solved {
        assert!(
            h >= h_min && h <= h_max,
            "Solved h must be in [h_min, h_max]: {h}"
        );
    }
}

// -----------------------------------------------------------------------
// Omega correction tests
// -----------------------------------------------------------------------

#[test]
fn omega_correction_near_one_for_few_neighbours() {
    // With no neighbours, drho/dh = 0 → Ω = 1
    let pos_i = [0.0_f64; 3];
    let omega = omega_correction(pos_i, 0.1, &[], &[], 1000.0, 3.0);
    assert!((omega - 1.0).abs() < 1e-14, "No neighbours → Ω = 1");
}

/// Cubic lattice of `n`³ points with spacing `dx`, centred on the origin
/// (for odd `n` one particle sits exactly at the origin).
fn cubic_lattice(n: usize, dx: f64) -> Vec<[f64; 3]> {
    let centre = (n as f64 - 1.0) * 0.5;
    let mut points = Vec::with_capacity(n * n * n);
    for i in 0..n {
        for j in 0..n {
            for k in 0..n {
                points.push([
                    (i as f64 - centre) * dx,
                    (j as f64 - centre) * dx,
                    (k as f64 - centre) * dx,
                ]);
            }
        }
    }
    points
}

/// Central finite-difference estimate of ∂ρ/∂h from the density sum alone.
fn density_dh_central_difference(
    pos_i: [f64; 3],
    h: f64,
    positions: &[[f64; 3]],
    masses: &[f64],
) -> f64 {
    let delta = 1.0e-6 * h;
    let rho_plus = sph_density(pos_i, h + delta, positions, masses);
    let rho_minus = sph_density(pos_i, h - delta, positions, masses);
    (rho_plus - rho_minus) / (2.0 * delta)
}

#[test]
fn sph_density_dh_matches_central_difference() {
    // Irregular cluster spanning q < 1, 1 <= q < 2 and q >= 2 at every h used.
    let pos_i = [0.0_f64; 3];
    let positions: Vec<[f64; 3]> = vec![
        [0.0, 0.0, 0.0],
        [0.031, 0.0, 0.0],
        [-0.047, 0.012, 0.0],
        [0.0, 0.083, -0.021],
        [0.052, -0.044, 0.019],
        [-0.091, -0.063, 0.027],
        [0.0, 0.0, 0.118],
        [0.139, 0.021, -0.057],
        [-0.171, 0.0, 0.0],
        [0.0, -0.205, 0.0],
    ];
    let masses: Vec<f64> = (0..positions.len())
        .map(|i| 0.001 * (1.0 + 0.1 * i as f64))
        .collect();
    for &h in &[0.04_f64, 0.07, 0.1, 0.13] {
        let analytic = sph_density_dh(pos_i, h, &positions, &masses);
        let numeric = density_dh_central_difference(pos_i, h, &positions, &masses);
        let tol = 1.0e-6 * numeric.abs().max(1.0e-9);
        assert!(
            (analytic - numeric).abs() <= tol,
            "h = {h}: analytic ∂ρ/∂h = {analytic}, finite difference = {numeric}"
        );
    }
}

#[test]
fn sph_density_dh_single_particle_follows_h_minus_three_scaling() {
    // A lone particle contributes W(0, h) = 1/(π h³); ∂/∂h = −3/(π h⁴).
    let h = 0.2_f64;
    let drho = sph_density_dh([0.0; 3], h, &[[0.0; 3]], &[2.5]);
    let expected = -3.0 * 2.5 / (std::f64::consts::PI * h.powi(4));
    assert!(
        (drho - expected).abs() <= 1.0e-12 * expected.abs(),
        "drho/dh = {drho}, expected {expected}"
    );
}

#[test]
fn omega_correction_near_unity_on_regular_lattice() {
    // 7³ lattice at h = 1.2 Δx: ρ is nearly independent of h, hence Ω ≈ 1.
    let dx = 0.1_f64;
    let h = 1.2 * dx;
    let rho0 = 1000.0_f64;
    let positions = cubic_lattice(7, dx);
    let masses = vec![rho0 * dx.powi(3); positions.len()];
    let pos_i = [0.0_f64; 3];
    let rho = sph_density(pos_i, h, &positions, &masses);
    assert!(
        (rho - rho0).abs() / rho0 < 5.0e-3,
        "lattice density {rho} should be close to ρ₀ = {rho0}"
    );

    let omega = omega_correction(pos_i, h, &positions, &masses, rho, 3.0);
    assert!(
        (0.9..=1.1).contains(&omega),
        "Ω at the lattice centre must lie in [0.9, 1.1], got {omega}"
    );

    // Independent check against Ω = 1 + h/(3ρ) ∂ρ/∂h from finite differences.
    let numeric_dh = density_dh_central_difference(pos_i, h, &positions, &masses);
    let omega_fd = 1.0 + h * numeric_dh / (3.0 * rho);
    assert!(
        (omega - omega_fd).abs() < 1.0e-6,
        "Ω = {omega} disagrees with finite-difference value {omega_fd}"
    );
    // Value for the cubic-spline lattice at h = 1.2 Δx (scale invariant).
    assert!(
        (omega - 0.9809).abs() < 1.0e-3,
        "Ω = {omega}, expected ≈ 0.981"
    );
}

#[test]
fn omega_correction_below_one_when_density_falls_with_h() {
    // Two particles close together: ρ(h) is dominated by the self term and
    // decreases with h, so ∂ρ/∂h < 0 and Springel–Hernquist Ω = 1 + … < 1.
    let pos_i = [0.0_f64; 3];
    let positions = vec![[0.0, 0.0, 0.0], [0.02, 0.0, 0.0]];
    let masses = vec![0.01, 0.01];
    let h = 0.1;
    let rho = sph_density(pos_i, h, &positions, &masses);
    assert!(sph_density_dh(pos_i, h, &positions, &masses) < 0.0);
    let omega = omega_correction(pos_i, h, &positions, &masses, rho, 3.0);
    assert!(
        omega < 1.0 && omega > 0.0,
        "Ω must be in (0, 1) when ρ falls with h, got {omega}"
    );
}

#[test]
fn solve_h_newton_converges_for_compressed_lattice_target() {
    // Lattice of unit-mass particles with Δx = 1 has natural density ≈ 1;
    // a target of ρ(h = 0.7) ≈ 1.19 is a compressed state that needs a
    // smaller h.  Newton needs the exact ∂ρ/∂h to converge from h = 1.3.
    let positions = cubic_lattice(7, 1.0);
    let masses = vec![1.0_f64; positions.len()];
    let pos_i = [0.0_f64; 3];
    let h_true = 0.7;
    let rho_target = sph_density(pos_i, h_true, &positions, &masses);
    assert!(
        rho_target > 1.15,
        "target density {rho_target} is compressed"
    );

    let solved = solve_h_newton(
        pos_i, 1.3, rho_target, &positions, &masses, 15, 1.0e-10, 0.1, 5.0,
    );
    let h = solved.expect("Newton iteration must converge with the exact Jacobian");
    assert!(
        (h - h_true).abs() < 1.0e-7,
        "solved h = {h}, expected {h_true}"
    );
    let rho = sph_density(pos_i, h, &positions, &masses);
    assert!((rho - rho_target).abs() / rho_target < 1.0e-9);
}

// -----------------------------------------------------------------------
// Variable-h kernel gradient tests
// -----------------------------------------------------------------------

#[test]
fn variable_h_gradient_zero_at_coincident_positions() {
    let r_ij = [0.0_f64; 3];
    let grad = variable_h_kernel_gradient(r_ij, 0.1, 0.1, 1.0, 1.0);
    let mag = (grad[0] * grad[0] + grad[1] * grad[1] + grad[2] * grad[2]).sqrt();
    assert!(mag < 1e-14, "Zero separation → zero gradient");
}

#[test]
fn variable_h_gradient_symmetric_equal_h() {
    // With h_i = h_j and Ω_i = Ω_j, the symmetric gradient should be the same
    // as the standard kernel gradient.
    let r_ij = [0.05_f64, 0.0, 0.0];
    let h = 0.1;
    let kernel = CubicSplineKernel;
    let r = 0.05;
    let gw = kernel.grad_w(r, h);

    let grad = variable_h_kernel_gradient(r_ij, h, h, 1.0, 1.0);
    let expected = gw / r * r_ij[0];
    assert!(
        (grad[0] - expected).abs() < 1e-12,
        "Equal-h gradient: {}, expected {}",
        grad[0],
        expected
    );
}

// -----------------------------------------------------------------------
// H bounds tests
// -----------------------------------------------------------------------

#[test]
fn apply_h_bounds_clamps_values() {
    let mut hs = vec![0.001, 0.1, 5.0, 0.5];
    apply_h_bounds(&mut hs, 0.01, 1.0);
    assert!((hs[0] - 0.01).abs() < 1e-14, "Too small clamped to h_min");
    assert!((hs[1] - 0.1).abs() < 1e-14, "In range unchanged");
    assert!((hs[2] - 1.0).abs() < 1e-14, "Too large clamped to h_max");
    assert!((hs[3] - 0.5).abs() < 1e-14, "In range unchanged");
}

#[test]
fn h_bounds_from_reference_uses_fractions() {
    let h0 = 0.1;
    let (h_min, h_max) = h_bounds_from_reference(h0);
    assert!((h_min - h0 * H_MIN_FRACTION).abs() < 1e-14);
    assert!((h_max - h0 * H_MAX_FRACTION).abs() < 1e-14);
}

// -----------------------------------------------------------------------
// Neighbour count estimate tests
// -----------------------------------------------------------------------

#[test]
fn neighbour_count_single_at_origin() {
    let pos_i = [0.0_f64; 3];
    let positions = vec![[0.0_f64; 3]]; // the particle itself
    let count = estimate_neighbour_count(pos_i, 0.1, &positions);
    assert_eq!(count, 1, "Should count itself");
}

#[test]
fn neighbour_count_excludes_outside_radius() {
    let pos_i = [0.0_f64; 3];
    let positions = vec![[1.0_f64, 0.0, 0.0]]; // far away
    let count = estimate_neighbour_count(pos_i, 0.1, &positions);
    assert_eq!(count, 0, "Far particle should not be counted");
}

#[test]
fn neighbour_count_uniform_grid() {
    let pos_i = [0.0_f64; 3];
    // 3×3×3 grid of particles at spacing 0.05 → all within 2h = 0.2 (h=0.1)
    let positions: Vec<[f64; 3]> = (0..3)
        .flat_map(|ix| {
            (0..3).flat_map(move |iy| {
                (0..3).map(move |iz| {
                    [
                        (ix as f64 - 1.0) * 0.05,
                        (iy as f64 - 1.0) * 0.05,
                        (iz as f64 - 1.0) * 0.05,
                    ]
                })
            })
        })
        .collect();
    let count = estimate_neighbour_count(pos_i, 0.15, &positions);
    assert!(count > 0, "Should find at least some neighbours: {count}");
}

// -----------------------------------------------------------------------
// Particle split/merge tests
// -----------------------------------------------------------------------

#[test]
fn split_particle_gives_eight_children() {
    let parent = AdaptiveHParticle::new([0.0; 3], [0.0; 3], 1.0, 0.1, 1000.0, 0);
    let children = split_particle(&parent);
    assert_eq!(children.len(), 8, "Should produce 8 children");
}

#[test]
fn split_conserves_mass() {
    let parent_mass = 2.5;
    let parent = AdaptiveHParticle::new([0.0; 3], [0.0; 3], parent_mass, 0.1, 1000.0, 0);
    let children = split_particle(&parent);
    let total_child_mass: f64 = children.iter().map(|c| c.mass).sum();
    assert!(
        (total_child_mass - parent_mass).abs() < 1e-14,
        "Split must conserve mass: {total_child_mass} vs {parent_mass}"
    );
}

#[test]
fn split_centre_of_mass_conserved() {
    let parent_pos = [1.0_f64, 2.0, 3.0];
    let parent = AdaptiveHParticle::new(parent_pos, [0.0; 3], 1.0, 0.1, 1000.0, 0);
    let children = split_particle(&parent);
    let n = children.len() as f64;
    let cm = [
        children.iter().map(|c| c.position[0]).sum::<f64>() / n,
        children.iter().map(|c| c.position[1]).sum::<f64>() / n,
        children.iter().map(|c| c.position[2]).sum::<f64>() / n,
    ];
    for k in 0..3 {
        assert!(
            (cm[k] - parent_pos[k]).abs() < 1e-12,
            "Centre of mass dim {k}: {}, expected {}",
            cm[k],
            parent_pos[k]
        );
    }
}

#[test]
fn merge_particles_conserves_mass() {
    let particles: Vec<AdaptiveHParticle> = (0..8)
        .map(|i| AdaptiveHParticle::new([i as f64 * 0.01; 3], [0.0; 3], 0.125, 0.05, 1000.0, 1))
        .collect();
    let total_before: f64 = particles.iter().map(|p| p.mass).sum();
    let merged = merge_particles(&particles);
    assert!(
        (merged.mass - total_before).abs() < 1e-14,
        "Merge must conserve mass: {} vs {}",
        merged.mass,
        total_before
    );
}

#[test]
fn merge_particles_centre_of_mass() {
    // Two equal-mass particles at x=0 and x=2 → merged at x=1
    let p1 = AdaptiveHParticle::new([0.0, 0.0, 0.0], [0.0; 3], 1.0, 0.05, 1000.0, 1);
    let p2 = AdaptiveHParticle::new([2.0, 0.0, 0.0], [0.0; 3], 1.0, 0.05, 1000.0, 1);
    let merged = merge_particles(&[p1, p2]);
    assert!(
        (merged.position[0] - 1.0).abs() < 1e-14,
        "CoM should be at x=1"
    );
}

#[test]
fn split_then_merge_preserves_mass() {
    let parent = AdaptiveHParticle::new([0.0; 3], [1.0, 0.0, 0.0], 2.0, 0.1, 1000.0, 0);
    let original_mass = parent.mass;
    let children = split_particle(&parent);
    let merged = merge_particles(&children);
    assert!(
        (merged.mass - original_mass).abs() < 1e-14,
        "Round-trip mass: {} vs {}",
        merged.mass,
        original_mass
    );
}

// -----------------------------------------------------------------------
// AdaptiveHManager tests
// -----------------------------------------------------------------------

#[test]
fn manager_update_smoothing_lengths_hoover() {
    let h0 = 0.1;
    let rho0 = 1000.0;
    let mut mgr = AdaptiveHManager::new(rho0, h0, 32.0, 3);

    // Add particle at double the reference density → h should halve (3D)
    mgr.add_particle(AdaptiveHParticle::new(
        [0.0; 3],
        [0.0; 3],
        1.0,
        h0,
        2.0 * rho0,
        0,
    ));
    mgr.update_smoothing_lengths();

    let expected_h = hoover_h(h0, rho0, 2.0 * rho0, 3);
    let actual_h = mgr.particles[0].h;
    assert!(
        (actual_h - expected_h).abs() < 1e-12,
        "Hoover h mismatch: {actual_h} vs {expected_h}"
    );
}

#[test]
fn manager_total_mass_conserved_after_refine() {
    let mut mgr = AdaptiveHManager::new(1000.0, 0.1, 32.0, 3);
    for i in 0..4 {
        mgr.add_particle(AdaptiveHParticle::new(
            [i as f64 * 0.1, 0.0, 0.0],
            [0.0; 3],
            1.0,
            0.5, // Large h → will be split
            1000.0,
            0,
        ));
    }
    let mass_before = mgr.total_mass();
    mgr.refine_coarse_particles(0.3); // threshold < 0.5 → all 4 split
    let mass_after = mgr.total_mass();
    assert!(
        (mass_after - mass_before).abs() < 1e-12,
        "Mass after refine: {mass_after} vs {mass_before}"
    );
}

#[test]
fn manager_omega_all_finite() {
    let mut mgr = AdaptiveHManager::new(1000.0, 0.1, 32.0, 3);
    for i in 0..5 {
        mgr.add_particle(AdaptiveHParticle::new(
            [i as f64 * 0.05, 0.0, 0.0],
            [0.0; 3],
            0.001,
            0.1,
            1000.0,
            0,
        ));
    }
    // Update densities from SPH sum
    let positions: Vec<[f64; 3]> = mgr.particles.iter().map(|p| p.position).collect();
    let masses: Vec<f64> = mgr.particles.iter().map(|p| p.mass).collect();
    for p in &mut mgr.particles {
        p.density = sph_density(p.position, p.h, &positions, &masses).max(1e-14);
    }

    let omegas = mgr.compute_omega();
    for (i, &w) in omegas.iter().enumerate() {
        assert!(w.is_finite(), "Omega[{i}] should be finite: {w}");
    }
}

// -----------------------------------------------------------------------
// Zhu-Fox criterion tests
// -----------------------------------------------------------------------

#[test]
fn zhu_fox_uniform_density_zero() {
    let pos_i = [0.0_f64; 3];
    let rho = 1000.0;
    let h = 0.1;
    // All neighbours at same density → epsilon = 0
    let neighbour_pos = vec![[0.05, 0.0, 0.0], [-0.05, 0.0, 0.0]];
    let neighbour_rho = vec![rho, rho];
    let eps = zhu_fox_criterion(pos_i, rho, h, &neighbour_pos, &neighbour_rho);
    assert!(eps.abs() < 1e-14, "Uniform density → eps=0, got {eps}");
}

#[test]
fn zhu_fox_varying_density_positive() {
    let pos_i = [0.0_f64; 3];
    let rho = 1000.0;
    let h = 0.1;
    let neighbour_pos = vec![[0.05, 0.0, 0.0]];
    let neighbour_rho = vec![1500.0_f64]; // 50% higher
    let eps = zhu_fox_criterion(pos_i, rho, h, &neighbour_pos, &neighbour_rho);
    assert!(
        (eps - 0.5).abs() < 1e-10,
        "Should report 50% variation, got {eps}"
    );
}

#[test]
fn zhu_fox_far_neighbour_ignored() {
    let pos_i = [0.0_f64; 3];
    let rho = 1000.0;
    let h = 0.05; // small h
    // Neighbour is at r = 0.1 > h, so outside the criterion radius
    let neighbour_pos = vec![[0.1, 0.0, 0.0]];
    let neighbour_rho = vec![5000.0_f64];
    let eps = zhu_fox_criterion(pos_i, rho, h, &neighbour_pos, &neighbour_rho);
    assert!(eps.abs() < 1e-14, "Far neighbour should be ignored: {eps}");
}

#[test]
fn should_refine_zhu_fox_at_max_level_returns_false() {
    let pos_i = [0.0_f64; 3];
    let rho = 1000.0;
    let result = should_refine_zhu_fox(pos_i, rho, 0.1, 5, 5, &[], &[], 0.1);
    assert!(!result, "At max_level should never refine");
}

#[test]
fn should_refine_zhu_fox_large_variation_returns_true() {
    let pos_i = [0.0_f64; 3];
    let rho = 1000.0;
    let h = 0.1;
    let neighbour_pos = vec![[0.05, 0.0, 0.0]];
    let neighbour_rho = vec![2000.0_f64]; // 100% variation
    let result = should_refine_zhu_fox(pos_i, rho, h, 0, 5, &neighbour_pos, &neighbour_rho, 0.5);
    assert!(result, "100% variation > threshold 50% → should refine");
}

// -----------------------------------------------------------------------
// sph_drho_dt tests
// -----------------------------------------------------------------------

#[test]
fn drho_dt_zero_for_zero_velocity() {
    let pos_i = [0.0_f64; 3];
    let vel_i = [0.0; 3];
    let rho_i = 1000.0;
    let h = 0.1;
    let neighbour_pos = vec![[0.05, 0.0, 0.0]];
    let neighbour_vel = vec![[0.0_f64; 3]]; // same velocity → zero relative
    let neighbour_mass = vec![0.001_f64];
    let neighbour_rho = vec![1000.0_f64];
    let drho = sph_drho_dt(
        pos_i,
        vel_i,
        rho_i,
        h,
        &neighbour_pos,
        &neighbour_vel,
        &neighbour_mass,
        &neighbour_rho,
    );
    assert!(
        drho.abs() < 1e-14,
        "Zero relative velocity → dρ/dt=0, got {drho}"
    );
}

#[test]
fn drho_dt_finite_for_nonzero_velocity() {
    let pos_i = [0.0_f64; 3];
    let vel_i = [1.0, 0.0, 0.0];
    let rho_i = 1000.0;
    let h = 0.1;
    let neighbour_pos = vec![[0.05, 0.0, 0.0]];
    let neighbour_vel = vec![[0.0_f64; 3]];
    let neighbour_mass = vec![0.001_f64];
    let neighbour_rho = vec![1000.0_f64];
    let drho = sph_drho_dt(
        pos_i,
        vel_i,
        rho_i,
        h,
        &neighbour_pos,
        &neighbour_vel,
        &neighbour_mass,
        &neighbour_rho,
    );
    assert!(drho.is_finite(), "dρ/dt should be finite");
}

// -----------------------------------------------------------------------
// update_h_continuity tests
// -----------------------------------------------------------------------

#[test]
fn update_h_continuity_no_change_for_zero_drho() {
    let h = 0.1;
    let h_new = update_h_continuity(h, 1000.0, 0.0, 0.001, 3.0, 0.01, 1.0);
    assert!((h_new - h).abs() < 1e-14, "Zero dρ/dt → h unchanged");
}

#[test]
fn update_h_continuity_increases_for_expanding_flow() {
    // Expanding flow: dρ/dt < 0 → dh/dt > 0 → h increases
    let h = 0.1;
    let h_new = update_h_continuity(h, 1000.0, -100.0, 0.001, 3.0, 0.01, 1.0);
    assert!(
        h_new > h,
        "Expanding flow should increase h: {h_new} vs {h}"
    );
}

#[test]
fn update_h_continuity_clamps_to_bounds() {
    let h = 0.1;
    // drho/dt very large → h would go very small
    let h_new = update_h_continuity(h, 1000.0, 1e10, 0.001, 3.0, 0.05, 1.0);
    assert!(h_new >= 0.05, "Should be clamped to h_min: {h_new}");
}

// -----------------------------------------------------------------------
// Gather / scatter density tests
// -----------------------------------------------------------------------

#[test]
fn gather_density_equals_standard() {
    let pos_i = [0.0_f64; 3];
    let h = 0.1;
    let positions = vec![[0.05_f64, 0.0, 0.0], [-0.05, 0.0, 0.0]];
    let masses = vec![0.001_f64, 0.001];
    let rho_gather = sph_density_gather(pos_i, h, &positions, &masses);
    let rho_std = sph_density(pos_i, h, &positions, &masses);
    assert!(
        (rho_gather - rho_std).abs() < 1e-12,
        "Gather == standard for uniform h"
    );
}

#[test]
fn scatter_density_positive() {
    let pos_i = [0.0_f64; 3];
    let positions = vec![[0.05_f64, 0.0, 0.0]];
    let masses = vec![0.001_f64];
    let hs = vec![0.1_f64];
    let rho = sph_density_scatter(pos_i, &positions, &masses, &hs);
    assert!(rho > 0.0, "Scatter density should be positive: {rho}");
}

#[test]
fn symmetric_density_between_gather_and_scatter() {
    let pos_i = [0.0_f64; 3];
    let h_i = 0.1;
    let positions = vec![[0.05_f64, 0.0, 0.0]];
    let masses = vec![0.001_f64];
    let hs_j = vec![0.08_f64]; // different h for neighbour
    let rho_g = sph_density_gather(pos_i, h_i, &positions, &masses);
    let rho_s = sph_density_scatter(pos_i, &positions, &masses, &hs_j);
    let rho_sym = sph_density_symmetric(pos_i, h_i, &positions, &masses, &hs_j);
    let expected = 0.5 * (rho_g + rho_s);
    assert!(
        (rho_sym - expected).abs() < 1e-14,
        "Symmetric = avg(gather, scatter)"
    );
}

// -----------------------------------------------------------------------
// Gather-scatter gradient test
// -----------------------------------------------------------------------

#[test]
fn gather_scatter_gradient_zero_at_origin() {
    let r_ij = [0.0_f64; 3];
    let g = gather_scatter_gradient(r_ij, 0.1, 0.1);
    let mag = (g[0] * g[0] + g[1] * g[1] + g[2] * g[2]).sqrt();
    assert!(mag < 1e-14, "Zero separation → zero gradient");
}

#[test]
fn gather_scatter_gradient_symmetric_equal_h() {
    // With h_i = h_j, gather-scatter == standard gradient
    let r_ij = [0.05_f64, 0.0, 0.0];
    let h = 0.1;
    let kernel = CubicSplineKernel;
    let r = 0.05;
    let gw = kernel.grad_w(r, h);
    let g = gather_scatter_gradient(r_ij, h, h);
    let expected_x = gw / r * r_ij[0];
    assert!(
        (g[0] - expected_x).abs() < 1e-12,
        "Equal-h gather-scatter == standard"
    );
}

// -----------------------------------------------------------------------
// 4-particle split tests
// -----------------------------------------------------------------------

#[test]
fn split_4_gives_four_children() {
    let parent = AdaptiveHParticle::new([0.0; 3], [0.0; 3], 1.0, 0.1, 1000.0, 0);
    let children = split_particle_4(&parent);
    assert_eq!(children.len(), 4);
}

#[test]
fn split_4_conserves_mass() {
    let parent = AdaptiveHParticle::new([0.0; 3], [0.0; 3], 3.0, 0.1, 1000.0, 0);
    let children = split_particle_4(&parent);
    let total: f64 = children.iter().map(|c| c.mass).sum();
    assert!((total - 3.0).abs() < 1e-14, "Split-4 should conserve mass");
}

#[test]
fn split_4_increments_level() {
    let parent = AdaptiveHParticle::new([0.0; 3], [0.0; 3], 1.0, 0.1, 1000.0, 2);
    let children = split_particle_4(&parent);
    for c in &children {
        assert_eq!(c.level, 3, "Children should be at level parent.level + 1");
    }
}

// -----------------------------------------------------------------------
// merge_pair test
// -----------------------------------------------------------------------

#[test]
fn merge_pair_conserves_mass() {
    let a = AdaptiveHParticle::new([0.0; 3], [1.0, 0.0, 0.0], 0.5, 0.1, 1000.0, 1);
    let b = AdaptiveHParticle::new([0.1; 3], [0.0, 1.0, 0.0], 0.5, 0.1, 1000.0, 1);
    let merged = merge_pair(&a, &b);
    assert!((merged.mass - 1.0).abs() < 1e-14);
}

#[test]
fn merge_pair_momentum_conserved() {
    let a = AdaptiveHParticle::new([0.0; 3], [2.0, 0.0, 0.0], 1.0, 0.1, 1000.0, 1);
    let b = AdaptiveHParticle::new([0.0; 3], [0.0, 0.0, 0.0], 1.0, 0.1, 1000.0, 1);
    let merged = merge_pair(&a, &b);
    // Momentum: m_a * v_a + m_b * v_b = 1*2 + 1*0 = 2 kg m/s
    // merged mass = 2, merged velocity should be 1.0 m/s
    assert!(
        (merged.velocity[0] - 1.0).abs() < 1e-14,
        "Velocity after merge: {}",
        merged.velocity[0]
    );
}

// -----------------------------------------------------------------------
// max_h_ratio test
// -----------------------------------------------------------------------

#[test]
fn max_h_ratio_uniform_is_one() {
    let positions = vec![[0.0_f64; 3], [0.05, 0.0, 0.0], [-0.05, 0.0, 0.0]];
    let hs = vec![0.1_f64; 3];
    let ratio = max_h_ratio(&positions, &hs, 0.2);
    assert!(
        (ratio - 1.0).abs() < 1e-10,
        "Uniform h → ratio 1.0, got {ratio}"
    );
}

#[test]
fn max_h_ratio_detects_difference() {
    let positions = vec![[0.0_f64; 3], [0.05, 0.0, 0.0]];
    let hs = vec![0.2_f64, 0.1_f64]; // ratio 2
    let ratio = max_h_ratio(&positions, &hs, 0.2);
    assert!(
        (ratio - 2.0).abs() < 1e-10,
        "Expected ratio 2.0, got {ratio}"
    );
}

// -----------------------------------------------------------------------
// kernel_mass_conservation_error test
// -----------------------------------------------------------------------

#[test]
fn kernel_mass_conservation_finite() {
    let pos_i = [0.0_f64; 3];
    let h = 0.1;
    let positions: Vec<[f64; 3]> = (0..8)
        .map(|i| {
            let angle = i as f64 * std::f64::consts::TAU / 8.0;
            [0.07 * angle.cos(), 0.07 * angle.sin(), 0.0]
        })
        .collect();
    let masses = vec![0.001_f64; 8];
    let err = kernel_mass_conservation_error(pos_i, h, &positions, &masses, 1000.0);
    assert!(
        err.is_finite(),
        "Mass conservation error should be finite: {err}"
    );
}

// -----------------------------------------------------------------------
// variable_h_pressure_accel test
// -----------------------------------------------------------------------

#[test]
fn variable_h_pressure_accel_zero_for_zero_pressure() {
    let pos_i = [0.0_f64; 3];
    let accel = variable_h_pressure_accel(
        pos_i,
        VariableHParticleState {
            pressure: 0.0,
            rho: 1000.0,
            h: 0.1,
            omega: 1.0,
        },
        &[[0.05, 0.0, 0.0_f64]],
        &[0.001],
        &[0.0], // neighbour pressure = 0
        &[1000.0],
        &[0.1],
        &[1.0],
    );
    let mag = (accel[0] * accel[0] + accel[1] * accel[1] + accel[2] * accel[2]).sqrt();
    assert!(mag < 1e-14, "Zero pressure → zero force: {mag}");
}

#[test]
fn variable_h_pressure_accel_finite_for_nonzero_pressure() {
    let pos_i = [0.0_f64; 3];
    let accel = variable_h_pressure_accel(
        pos_i,
        VariableHParticleState {
            pressure: 1000.0,
            rho: 1000.0,
            h: 0.1,
            omega: 1.0,
        },
        &[[0.05, 0.0, 0.0_f64]],
        &[0.001],
        &[500.0],
        &[1000.0],
        &[0.1],
        &[1.0],
    );
    for a in accel {
        assert!(
            a.is_finite(),
            "Acceleration component should be finite: {a}"
        );
    }
}
