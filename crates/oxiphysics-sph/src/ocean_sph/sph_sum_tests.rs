//! Regression tests for the Wendland C2 gradient, the neighbour-list SPH sums
//! and the neighbour-list based [`OceanSimulation::step`].

use super::*;

/// Relative-tolerance comparison.
fn rel_close(a: f64, b: f64, tol: f64) -> bool {
    if b.abs() < 1e-30 {
        return a.abs() < tol;
    }
    ((a - b) / b).abs() < tol
}

// ── Wendland C2 gradient normalisation ────────────────────────────────

#[test]
fn test_wendland_c2_grad_matches_finite_difference_of_kernel() {
    for &h in &[1.0_f64, 0.37, 2.5] {
        for &q in &[0.3_f64, 0.9, 1.5, 1.9] {
            let r = q * h;
            let eps = 1.0e-6 * h;
            let numeric = (wendland_c2(r + eps, h) - wendland_c2(r - eps, h)) / (2.0 * eps);
            let analytic = wendland_c2_grad(r, h);
            assert!(
                rel_close(analytic, numeric, 1.0e-6),
                "h={h} q={q}: analytic {analytic} vs finite difference {numeric}"
            );
        }
    }
}

#[test]
fn test_wendland_c2_grad_satisfies_divergence_identity() {
    // For a kernel normalised to one, ∫ r W'(r) 4π r² dr = −3 (integration
    // by parts); a mis-normalised gradient would give −3·(factor).
    let h = 0.8;
    let n = 4000;
    let dr = 2.0 * h / n as f64;
    let integral: f64 = (0..n)
        .map(|i| {
            let r = (i as f64 + 0.5) * dr;
            r * wendland_c2_grad(r, h) * 4.0 * PI * r * r * dr
        })
        .sum();
    assert!((integral + 3.0).abs() < 1.0e-3, "integral = {integral}");
}

#[test]
fn test_wendland_c2_grad_is_negative_inside_and_zero_outside_support() {
    let h = 1.3;
    for i in 1..40 {
        let r = 2.0 * h * i as f64 / 40.0;
        assert!(wendland_c2_grad(r, h) < 0.0, "r={r}");
    }
    assert_eq!(wendland_c2_grad(2.0 * h, h), 0.0);
    assert_eq!(wendland_c2_grad(3.0 * h, h), 0.0);
    assert_eq!(wendland_c2_grad(0.0, h), 0.0);
    assert_eq!(wendland_c2_grad(0.5, 0.0), 0.0);
}

#[test]
fn test_sph_pressure_acc_two_particles_uses_exact_kernel_gradient() {
    // Two particles at r = h along x with p = 1e4, ρ = 1000, m = 2:
    // a_i = −m (2p/ρ²) W'(h) r̂_ij and r̂_ij = −x̂ for i on the left.
    let h = 0.5;
    let (m, p, rho) = (2.0_f64, 1.0e4, 1000.0);
    let neighbours = [([h, 0.0, 0.0], m, p, rho)];
    let acc = sph_pressure_acc([0.0; 3], p, rho, h, &neighbours);
    // dW/dr at q = 1 for the Wendland C2 kernel: σ/h · (2/16 − 2/8 · 3) = σ/h · (−5/8)
    let sigma = 21.0 / (16.0 * PI * h * h * h);
    let dw_dr = sigma / h * (2.0 * 0.5_f64.powi(4) - 2.0 * 0.5_f64.powi(3) * 3.0);
    let expected = m * (2.0 * p / (rho * rho)) * dw_dr; // along +x because r̂_ij = −x̂
    assert!(
        rel_close(acc[0], expected, 1.0e-12),
        "{} vs {expected}",
        acc[0]
    );
    assert!(
        acc[0] < 0.0,
        "particle i is pushed away from j (towards −x)"
    );
    assert_eq!(acc[1], 0.0);
    assert_eq!(acc[2], 0.0);
}

// ── Neighbour-list SPH sums ───────────────────────────────────────────

/// Deterministic pseudo-random number in [0, 1) (64-bit LCG).
fn lcg_unit(state: &mut u64) -> f64 {
    *state = state
        .wrapping_mul(6_364_136_223_846_793_005)
        .wrapping_add(1_442_695_040_888_963_407);
    ((*state >> 11) as f64) / ((1_u64 << 53) as f64)
}

/// `nx × ny × nz` lattice of spacing `dx` with a deterministic jitter of
/// at most `jitter` along every axis.
fn jittered_lattice(nx: usize, ny: usize, nz: usize, dx: f64, jitter: f64) -> Vec<[f64; 3]> {
    let mut state = 0x9E37_79B9_7F4A_7C15_u64;
    let mut out = Vec::with_capacity(nx * ny * nz);
    for i in 0..nx {
        for j in 0..ny {
            for k in 0..nz {
                let mut pos = [i as f64 * dx, j as f64 * dx, k as f64 * dx];
                for c in &mut pos {
                    *c += jitter * (2.0 * lcg_unit(&mut state) - 1.0);
                }
                out.push(pos);
            }
        }
    }
    out
}

fn all_pairs_density(positions: &[[f64; 3]], masses: &[f64], h: f64, i: usize) -> f64 {
    let neigh: Vec<([f64; 3], f64)> = positions
        .iter()
        .copied()
        .zip(masses.iter().copied())
        .collect();
    sph_density_sum(positions[i], h, &neigh)
}

#[test]
fn test_sph_density_sum_nl_matches_all_pairs_sum() {
    let h = 0.65;
    let positions = jittered_lattice(6, 5, 4, 0.5, 0.1);
    let masses: Vec<f64> = (0..positions.len())
        .map(|i| 100.0 + 0.5 * (i % 7) as f64)
        .collect();
    let lists = SpatialHash3D::all_neighbors(&positions, 2.0 * h);
    for (i, list) in lists.iter().enumerate() {
        let fast = sph_density_sum_nl(i, &positions, &masses, h, list);
        let slow = all_pairs_density(&positions, &masses, h, i);
        assert!(rel_close(fast, slow, 1.0e-12), "i={i}: {fast} vs {slow}");
    }
}

#[test]
fn test_sph_pressure_acc_nl_matches_all_pairs_sum() {
    let h = 0.65;
    let positions = jittered_lattice(5, 5, 5, 0.5, 0.1);
    let n = positions.len();
    let masses = vec![RHO_WATER * 0.125; n];
    let lists = SpatialHash3D::all_neighbors(&positions, 2.0 * h);
    let densities: Vec<f64> = (0..n)
        .map(|i| sph_density_sum_nl(i, &positions, &masses, h, &lists[i]))
        .collect();
    let pressures: Vec<f64> = densities
        .iter()
        .map(|&rho| tait_pressure(rho, RHO_WATER, 40.0, 7.0))
        .collect();
    for i in 0..n {
        let fast =
            sph_pressure_acc_nl(i, &positions, &masses, &pressures, &densities, h, &lists[i]);
        let neigh: Vec<([f64; 3], f64, f64, f64)> = (0..n)
            .map(|j| (positions[j], masses[j], pressures[j], densities[j]))
            .collect();
        let slow = sph_pressure_acc(positions[i], pressures[i], densities[i], h, &neigh);
        for c in 0..3 {
            assert!(
                (fast[c] - slow[c]).abs() <= 1.0e-9 * (1.0 + slow[c].abs()),
                "i={i} c={c}: {} vs {}",
                fast[c],
                slow[c]
            );
        }
    }
}

#[test]
fn test_neighbour_list_sums_ignore_self_and_invalid_indices() {
    let positions = [[0.0; 3], [0.4, 0.0, 0.0]];
    let masses = [3.0, 5.0];
    let h = 0.5;
    let clean = sph_density_sum_nl(0, &positions, &masses, h, &[1]);
    // Self index and out-of-range index must not change the result.
    let noisy = sph_density_sum_nl(0, &positions, &masses, h, &[0, 1, 99]);
    assert_eq!(clean, noisy);
    let expected = 3.0 * wendland_c2(0.0, h) + 5.0 * wendland_c2(0.4, h);
    assert!(rel_close(clean, expected, 1.0e-14));
    assert_eq!(sph_density_sum_nl(7, &positions, &masses, h, &[0]), 0.0);
    let acc = sph_pressure_acc_nl(7, &positions, &masses, &[1.0, 1.0], &[1.0, 1.0], h, &[0]);
    assert_eq!(acc, [0.0; 3]);
}

// ── OceanSimulation::step ─────────────────────────────────────────────

fn quiet_config(dt: f64, cs: f64) -> OceanSimConfig {
    OceanSimConfig {
        dt,
        cs,
        enable_breaking: false,
        enable_foam: false,
        enable_stokes: false,
        enable_tides: false,
        ..Default::default()
    }
}

#[test]
fn test_ocean_step_matches_all_pairs_reference_with_variable_h() {
    let config = quiet_config(0.002, 30.0);
    let dx = 0.5;
    let positions = jittered_lattice(8, 8, 6, dx, 0.08);
    let n = positions.len();
    let mass = RHO_WATER * dx * dx * dx;
    let smoothing: Vec<f64> = (0..n)
        .map(|i| if i % 5 == 0 { 0.9 } else { 0.65 })
        .collect();
    let mut sim = OceanSimulation::new(config.clone(), JonswapParams::new(0.1));
    for (pos, &h) in positions.iter().zip(smoothing.iter()) {
        sim.particles.push(OceanParticle::new(*pos, mass, h));
    }

    // Independent all-pairs reference (the original O(N²) formulation).
    let masses = vec![mass; n];
    let ref_rho: Vec<f64> = (0..n)
        .map(|i| all_pairs_density(&positions, &masses, smoothing[i], i))
        .collect();
    let ref_p: Vec<f64> = ref_rho
        .iter()
        .map(|&rho| tait_pressure(rho, config.rho0, config.cs, 7.0))
        .collect();
    let ref_acc: Vec<[f64; 3]> = (0..n)
        .map(|i| {
            let neigh: Vec<([f64; 3], f64, f64, f64)> = (0..n)
                .map(|j| (positions[j], masses[j], ref_p[j], ref_rho[j]))
                .collect();
            let a = sph_pressure_acc(positions[i], ref_p[i], ref_rho[i], smoothing[i], &neigh);
            [a[0], a[1], a[2] - G]
        })
        .collect();

    sim.step();

    for i in 0..n {
        let p = &sim.particles[i];
        assert!(rel_close(p.rho, ref_rho[i], 1.0e-12), "rho[{i}]");
        assert!(
            (p.pressure - ref_p[i]).abs() <= 1.0e-9 * (1.0 + ref_p[i].abs()),
            "p[{i}]: {} vs {}",
            p.pressure,
            ref_p[i]
        );
        for c in 0..3 {
            assert!(
                (p.acc[c] - ref_acc[i][c]).abs() <= 1.0e-9 * (1.0 + ref_acc[i][c].abs()),
                "acc[{i}][{c}]: {} vs {}",
                p.acc[c],
                ref_acc[i][c]
            );
        }
    }
}

#[test]
fn test_ocean_step_single_particle_falls_freely() {
    let config = quiet_config(0.01, 30.0);
    let mut sim = OceanSimulation::new(config, JonswapParams::new(0.1));
    let h = 0.7;
    sim.particles
        .push(OceanParticle::new([1.0, 2.0, 3.0], 12.0, h));
    sim.step();
    let p = &sim.particles[0];
    assert!(rel_close(p.rho, 12.0 * wendland_c2(0.0, h), 1.0e-14));
    assert_eq!(p.acc, [0.0, 0.0, -G]);
    assert!((p.vel[2] + G * 0.01).abs() < 1.0e-14);
}

#[test]
fn test_ocean_step_without_particles_or_smoothing_length() {
    let mut empty = OceanSimulation::new(quiet_config(0.01, 30.0), JonswapParams::new(0.1));
    empty.step();
    assert_eq!(empty.step, 1);
    assert!((empty.time - 0.01).abs() < 1.0e-15);

    // A non-positive smoothing length means no kernel support: the step must
    // neither panic nor produce non-finite state.
    let mut sim = OceanSimulation::new(quiet_config(0.01, 30.0), JonswapParams::new(0.1));
    sim.particles.push(OceanParticle::new([0.0; 3], 1.0, 0.0));
    sim.particles
        .push(OceanParticle::new([0.1, 0.0, 0.0], 1.0, 0.0));
    sim.step();
    for p in &sim.particles {
        assert!(p.rho == 0.0 && p.pos.iter().all(|c| c.is_finite()));
        assert_eq!(p.acc[0], 0.0);
    }
}

#[test]
fn test_ocean_step_stokes_drift_uses_each_particles_depth() {
    let jonswap = JonswapParams::new(0.1);
    let omega = jonswap.omega_peak;
    let mut sim = OceanSimulation::new(
        OceanSimConfig {
            dt: 0.01,
            ..Default::default()
        },
        jonswap,
    );
    for &(x, depth) in &[
        (0.0, 50.0),
        (10.0, 50.0),
        (20.0, 5.0),
        (30.0, 5.0),
        (40.0, 50.0),
    ] {
        let mut p = OceanParticle::new([x, 0.0, 0.0], 1000.0, 1.0);
        p.depth = depth;
        sim.particles.push(p);
    }
    sim.step();
    for p in &sim.particles {
        let k = dispersion_wavenumber(omega, p.depth.max(0.1));
        let expected = stokes_drift_deep(0.5, omega, k, p.pos[2], [1.0, 0.0]);
        for c in 0..3 {
            assert!(
                (p.stokes_drift[c] - expected[c]).abs() <= 1.0e-14,
                "depth {}: {:?} vs {expected:?}",
                p.depth,
                p.stokes_drift
            );
        }
    }
    assert!(
        (sim.particles[0].stokes_drift[0] - sim.particles[2].stokes_drift[0]).abs() > 1.0e-6,
        "the shallow and deep particles must use different wavenumbers"
    );
}
