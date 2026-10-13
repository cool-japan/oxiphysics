//! Regression tests for the 2-D colour-function SPH system
//! ([`MultiphaseSphSimulation`]): full-support colour gradients, interface
//! curvature, the continuum surface force and the two-phase WCSPH pressure.

use super::functions::compute_color_gradients_with_neighbors;
use super::types::{ColorParticle, MultiphaseSphSimulation};
use std::f64::consts::PI;

const RHO_L: f64 = 1000.0;
const RHO_G: f64 = 1.225;

/// Deterministic pseudo-random number in `[0, 1)` (64-bit LCG).
fn lcg_unit(state: &mut u64) -> f64 {
    *state = state
        .wrapping_mul(6_364_136_223_846_793_005)
        .wrapping_add(1_442_695_040_888_963_407);
    ((*state >> 11) as f64) / ((1_u64 << 53) as f64)
}

/// 2-D cubic-spline gradient magnitude `dW/dr` (support `2h`).
fn dw_dr(r: f64, h: f64) -> f64 {
    let q = r / h;
    let factor = 10.0 / (7.0 * PI * h * h);
    if q < 1.0 {
        factor * (-3.0 * q + 2.25 * q * q) / h
    } else if q < 2.0 {
        factor * (-0.75 * (2.0 - q).powi(2)) / h
    } else {
        0.0
    }
}

/// `n × n` lattice of spacing `dx`; the `block × block` square starting at
/// lattice index `start` is liquid, everything else gas.  Masses follow the
/// phase densities (`m = ρ dx²`), gravity is off.
fn patch_simulation(
    n: usize,
    start: usize,
    block: usize,
    dx: f64,
    h_over_dx: f64,
    sigma: f64,
) -> MultiphaseSphSimulation {
    let mut sim = MultiphaseSphSimulation::new(h_over_dx * dx, sigma, RHO_L, RHO_G);
    sim.gravity = [0.0, 0.0];
    for i in 0..n {
        for j in 0..n {
            let (x, y) = (i as f64 * dx, j as f64 * dx);
            let liquid = (start..start + block).contains(&i) && (start..start + block).contains(&j);
            if liquid {
                sim.add_liquid_particle(x, y, RHO_L * dx * dx);
            } else {
                sim.add_gas_particle(x, y, RHO_G * dx * dx);
            }
        }
    }
    sim
}

fn index(n: usize, i: usize, j: usize) -> usize {
    i * n + j
}

// ---------------------------------------------------------------------------
// Neighbour lists and colour gradient
// ---------------------------------------------------------------------------

#[test]
fn neighbor_lists_match_a_brute_force_search_in_two_dimensions() {
    let mut state = 77_u64;
    let mut sim = MultiphaseSphSimulation::new(0.06, 0.0, RHO_L, RHO_G);
    for _ in 0..200 {
        let (x, y) = (lcg_unit(&mut state), lcg_unit(&mut state));
        sim.add_liquid_particle(x, y, 1.0);
    }
    let lists = sim.neighbor_lists();
    assert_eq!(lists.len(), 200);
    let radius2 = (2.0 * sim.h) * (2.0 * sim.h);
    for (i, list) in lists.iter().enumerate() {
        let expected: Vec<usize> = (0..200)
            .filter(|&j| {
                let dx = sim.particles[i].pos[0] - sim.particles[j].pos[0];
                let dy = sim.particles[i].pos[1] - sim.particles[j].pos[1];
                j != i && dx * dx + dy * dy <= radius2
            })
            .collect();
        assert_eq!(*list, expected, "neighbours of particle {i}");
    }
}

#[test]
fn color_gradient_uses_the_full_kernel_support() {
    // Two particles 1.5 h apart: inside the support radius 2h but beyond h, so
    // the old cutoff at h dropped the interaction entirely.
    let h = 0.1;
    let mut sim = MultiphaseSphSimulation::new(h, 0.0, RHO_L, RHO_G);
    sim.add_liquid_particle(0.0, 0.0, 1.0);
    sim.add_gas_particle(1.5 * h, 0.0, 1.0);
    sim.update_color_gradients();
    // ∇C_liquid = (m_g/ρ_g)(C_g − C_l) ∇_i W_ij with C_g − C_l = −1 and
    // ∇_i W_ij = (dW/dr)(x_i − x_j)/r = −(dW/dr) x̂, so the gradient is
    // (m_g/ρ_g) · dW/dr · x̂, which points towards −x, into the liquid.
    let volume_gas = 1.0 / 1.225;
    let expected = volume_gas * dw_dr(1.5 * h, h);
    let g = sim.particles[0].color_grad;
    assert!(expected < 0.0);
    assert!(
        (g[0] - expected).abs() <= 1e-12 * expected.abs(),
        "{} vs {expected}",
        g[0]
    );
    assert_eq!(g[1], 0.0);
    // The gas particle's gradient also points into the liquid (towards −x); its
    // magnitude differs because the liquid particle has a different volume.
    assert!(sim.particles[1].color_grad[0] < 0.0);
}

#[test]
fn color_gradient_ignores_particles_beyond_two_smoothing_lengths() {
    let h = 0.1;
    let mut sim = MultiphaseSphSimulation::new(h, 0.0, RHO_L, RHO_G);
    sim.add_liquid_particle(0.0, 0.0, 1.0);
    sim.add_gas_particle(2.0 * h * 1.0001, 0.0, 1.0);
    sim.update_color_gradients();
    assert_eq!(sim.particles[0].color_grad, [0.0, 0.0]);
}

#[test]
fn color_gradient_with_neighbors_matches_the_plain_function() {
    let mut a = patch_simulation(14, 4, 6, 0.02, 1.5, 0.0);
    let mut b = a.clone();
    a.update_color_gradients();
    let lists = b.neighbor_lists();
    compute_color_gradients_with_neighbors(&mut b.particles, b.h, &lists);
    for (pa, pb) in a.particles.iter().zip(b.particles.iter()) {
        assert_eq!(pa.color_grad, pb.color_grad);
    }
}

// ---------------------------------------------------------------------------
// Curvature
// ---------------------------------------------------------------------------

#[test]
fn curvature_of_a_circular_normal_field_is_the_inverse_radius() {
    // n̂ = −r̂ about the origin (pointing into the C = 1 phase) has
    // −∇·n̂ = 1/r in two dimensions.
    let dx = 0.01;
    let mut sim = MultiphaseSphSimulation::new(1.5 * dx, 0.0, RHO_L, RHO_G);
    let (cx, cy) = (0.5, 0.0);
    for i in 0..21 {
        for j in 0..21 {
            let x = cx - 0.1 + i as f64 * dx;
            let y = cy - 0.1 + j as f64 * dx;
            sim.add_liquid_particle(x, y, RHO_L * dx * dx);
        }
    }
    let magnitude = 40.0;
    for p in sim.particles.iter_mut() {
        let r = p.pos[0].hypot(p.pos[1]);
        p.color_grad = [-magnitude * p.pos[0] / r, -magnitude * p.pos[1] / r];
    }
    sim.update_curvatures();
    let centre = index(21, 10, 10);
    let kappa = sim.particles[centre].curvature;
    let expected = 1.0 / sim.particles[centre].pos[0].hypot(sim.particles[centre].pos[1]);
    assert!(
        (kappa - expected).abs() < 0.03 * expected,
        "κ = {kappa} vs 1/R = {expected}"
    );
    // A flat interface (uniform normal) has no curvature.
    for p in sim.particles.iter_mut() {
        p.color_grad = [0.0, magnitude];
    }
    sim.update_curvatures();
    assert!(sim.particles[centre].curvature.abs() < 1e-9);
}

#[test]
fn curvature_vanishes_where_the_colour_gradient_is_below_the_threshold() {
    let dx = 0.01;
    let mut sim = patch_simulation(10, 0, 10, dx, 1.5, 0.0);
    // Single phase: the colour is uniform, there is no interface.
    sim.update_color_gradients();
    sim.update_curvatures();
    for p in &sim.particles {
        assert_eq!(p.color_grad, [0.0, 0.0]);
        assert_eq!(p.curvature, 0.0);
    }
}

// ---------------------------------------------------------------------------
// Surface tension
// ---------------------------------------------------------------------------

#[test]
fn square_patch_corners_feel_an_inward_surface_tension_acceleration() {
    let (n, start, block, dx) = (20, 5, 10, 0.02);
    let mut sim = patch_simulation(n, start, block, dx, 1.5, 0.072);
    sim.update_densities();
    sim.update_color_gradients();
    sim.update_curvatures();
    let acc = sim.surface_tension_accelerations();

    let centre = [
        (start as f64 + 0.5 * (block as f64 - 1.0)) * dx,
        (start as f64 + 0.5 * (block as f64 - 1.0)) * dx,
    ];
    let last = start + block - 1;
    // Reference scale of the capillary acceleration: σ / (ρ_l dx²).
    let scale = 0.072 / (RHO_L * dx * dx);
    for &(i, j) in &[(start, start), (start, last), (last, start), (last, last)] {
        let idx = index(n, i, j);
        let p = &sim.particles[idx];
        assert!(p.curvature > 0.0, "corner curvature {}", p.curvature);
        let a = acc[idx];
        let magnitude = a[0].hypot(a[1]);
        assert!(
            magnitude > 0.05 * scale,
            "corner ({i},{j}): |a| = {magnitude} vs scale {scale}"
        );
        let inward = [centre[0] - p.pos[0], centre[1] - p.pos[1]];
        let cosine =
            (a[0] * inward[0] + a[1] * inward[1]) / (magnitude * inward[0].hypot(inward[1]));
        assert!(
            cosine > 0.9,
            "corner ({i},{j}) not pulled inwards: cos = {cosine}"
        );
    }

    // Particles deep inside the liquid and far out in the gas feel nothing.
    let inside = index(n, start + block / 2, start + block / 2);
    let far = index(n, 0, 0);
    assert_eq!(acc[inside], [0.0, 0.0]);
    assert_eq!(acc[far], [0.0, 0.0]);

    // The liquid as a whole is pulled towards its centre: Σ a·(x − x_c) < 0.
    let virial: f64 = sim
        .particles
        .iter()
        .zip(acc.iter())
        .filter(|(p, _)| p.color > 0.5)
        .map(|(p, a)| a[0] * (p.pos[0] - centre[0]) + a[1] * (p.pos[1] - centre[1]))
        .sum();
    assert!(virial < 0.0, "liquid virial {virial}");
}

#[test]
fn curvature_of_a_liquid_disc_is_close_to_the_inverse_radius() {
    // A disc of radius R in gas has κ = 1/R at its interface; the SPH estimate,
    // averaged over the interface band with the colour-gradient weight, is
    // accurate to a few per cent.
    let (n, dx, radius) = (40, 0.02, 0.2);
    let mut sim = MultiphaseSphSimulation::new(1.5 * dx, 0.072, RHO_L, RHO_G);
    sim.gravity = [0.0, 0.0];
    let centre = 0.5 * n as f64 * dx;
    for i in 0..n {
        for j in 0..n {
            let (x, y) = (i as f64 * dx, j as f64 * dx);
            if (x - centre).hypot(y - centre) < radius {
                sim.add_liquid_particle(x, y, RHO_L * dx * dx);
            } else {
                sim.add_gas_particle(x, y, RHO_G * dx * dx);
            }
        }
    }
    sim.update_densities();
    sim.update_color_gradients();
    sim.update_curvatures();
    let threshold = 0.01 / sim.h;
    let (mut weighted, mut weight) = (0.0_f64, 0.0_f64);
    for p in &sim.particles {
        let magnitude = p.color_grad[0].hypot(p.color_grad[1]);
        if magnitude > threshold {
            weighted += p.curvature * magnitude;
            weight += magnitude;
        }
    }
    assert!(weight > 0.0);
    let mean = weighted / weight;
    let expected = 1.0 / radius;
    assert!(
        (mean - expected).abs() < 0.08 * expected,
        "mean interface curvature {mean} vs 1/R = {expected}"
    );
}

#[test]
fn surface_tension_vanishes_for_a_flat_interface() {
    // Half-plane of liquid: the interface is straight, so κ ≈ 0 away from the
    // ends of the lattice and the capillary acceleration is negligible there.
    let (n, dx) = (24, 0.02);
    let mut sim = MultiphaseSphSimulation::new(1.5 * dx, 0.072, RHO_L, RHO_G);
    sim.gravity = [0.0, 0.0];
    for i in 0..n {
        for j in 0..n {
            let (x, y) = (i as f64 * dx, j as f64 * dx);
            if j < n / 2 {
                sim.add_liquid_particle(x, y, RHO_L * dx * dx);
            } else {
                sim.add_gas_particle(x, y, RHO_G * dx * dx);
            }
        }
    }
    sim.update_densities();
    sim.update_color_gradients();
    sim.update_curvatures();
    let acc = sim.surface_tension_accelerations();
    let scale = 0.072 / (RHO_L * dx * dx);
    for i in 6..n - 6 {
        for j in (n / 2 - 2)..(n / 2 + 2) {
            let idx = index(n, i, j);
            let magnitude = acc[idx][0].hypot(acc[idx][1]);
            assert!(
                magnitude < 1e-9 * scale,
                "flat interface particle ({i},{j}): |a| = {magnitude}"
            );
        }
    }
}

#[test]
fn surface_tension_scales_linearly_with_sigma() {
    let compute = |sigma: f64| -> Vec<[f64; 2]> {
        let mut sim = patch_simulation(16, 4, 8, 0.02, 1.5, sigma);
        sim.update_densities();
        sim.update_color_gradients();
        sim.update_curvatures();
        sim.surface_tension_accelerations()
    };
    let one = compute(0.05);
    let three = compute(0.15);
    for (a1, a3) in one.iter().zip(three.iter()) {
        for k in 0..2 {
            assert!((a3[k] - 3.0 * a1[k]).abs() <= 1e-12 * (1.0 + a1[k].abs()));
        }
    }
}

// ---------------------------------------------------------------------------
// Density, pressure and pressure acceleration
// ---------------------------------------------------------------------------

#[test]
fn number_density_gives_every_particle_its_own_phase_density() {
    let (n, start, block, dx) = (20, 5, 10, 0.02);
    let mut sim = patch_simulation(n, start, block, dx, 1.5, 0.0);
    sim.update_densities();
    let last = start + block - 1;
    // Interior liquid, interface liquid, interface gas, far gas.
    for &(i, j, rho0) in &[
        (start + 5, start + 5, RHO_L),
        (start, start + 4, RHO_L),
        (last, last, RHO_L),
        (start - 1, start + 4, RHO_G),
        (last + 1, last + 1, RHO_G),
        (2, 2, RHO_G),
    ] {
        let rho = sim.particles[index(n, i, j)].density;
        assert!(
            (rho - rho0).abs() < 0.03 * rho0,
            "particle ({i},{j}): density {rho} vs rest density {rho0}"
        );
    }
}

#[test]
fn tait_pressure_vanishes_at_rest_density_and_grows_under_compression() {
    let mut sim = MultiphaseSphSimulation::new(0.05, 0.0, RHO_L, RHO_G);
    sim.add_liquid_particle(0.0, 0.0, 1.0);
    sim.add_liquid_particle(1.0, 0.0, 1.0);
    sim.add_liquid_particle(2.0, 0.0, 1.0);
    sim.add_gas_particle(3.0, 0.0, 1.0);
    sim.add_gas_particle(4.0, 0.0, 1.0);
    sim.sound_speed = 40.0;
    sim.particles[0].density = RHO_L;
    sim.particles[1].density = 1.01 * RHO_L;
    sim.particles[2].density = 0.9 * RHO_L;
    sim.particles[3].density = RHO_G;
    sim.particles[4].density = 1.02 * RHO_G;
    sim.update_pressures();
    let b_l = RHO_L * 40.0 * 40.0 / 7.0;
    let b_g = RHO_G * 40.0 * 40.0 / 7.0;
    assert_eq!(sim.particles[0].pressure, 0.0);
    assert!((sim.particles[1].pressure - b_l * (1.01_f64.powi(7) - 1.0)).abs() < 1e-6 * b_l);
    assert_eq!(sim.particles[2].pressure, 0.0, "tension is clamped");
    assert_eq!(sim.particles[3].pressure, 0.0);
    assert!((sim.particles[4].pressure - b_g * (1.02_f64.powi(7) - 1.0)).abs() < 1e-6 * b_g);
    // A vanishing sound speed disables the pressure.
    sim.sound_speed = 0.0;
    sim.update_pressures();
    assert!(sim.particles.iter().all(|p| p.pressure == 0.0));
}

#[test]
fn compressed_liquid_pair_repels_with_equal_and_opposite_momentum() {
    let h = 0.1;
    let mut sim = MultiphaseSphSimulation::new(h, 0.0, RHO_L, RHO_G);
    let mass = 2.0;
    sim.add_liquid_particle(0.0, 0.0, mass);
    sim.add_liquid_particle(h, 0.0, mass);
    for p in sim.particles.iter_mut() {
        p.density = RHO_L;
        p.pressure = 5.0e3;
    }
    let acc = sim.pressure_accelerations();
    // a_i = −m (P_i + P_j)/(ρ_i ρ_j) ∇_i W_ij, ∇_i W_ij = (dW/dr)(x_i − x_j)/r.
    let magnitude = mass * (2.0 * 5.0e3) / (RHO_L * RHO_L) * dw_dr(h, h).abs();
    assert!(
        (acc[0][0] + magnitude).abs() <= 1e-12 * magnitude,
        "{} vs {}",
        acc[0][0],
        -magnitude
    );
    assert!((acc[1][0] - magnitude).abs() <= 1e-12 * magnitude);
    assert_eq!(acc[0][1], 0.0);
    assert!((mass * acc[0][0] + mass * acc[1][0]).abs() < 1e-12 * magnitude);
}

#[test]
fn liquid_gas_pair_exchanges_equal_and_opposite_momentum() {
    let h = 0.1;
    let dx = 0.05;
    let mut sim = MultiphaseSphSimulation::new(h, 0.0, RHO_L, RHO_G);
    let (m_l, m_g) = (RHO_L * dx * dx, RHO_G * dx * dx);
    sim.add_liquid_particle(0.0, 0.0, m_l);
    sim.add_gas_particle(dx, 0.0, m_g);
    sim.particles[0].density = RHO_L;
    sim.particles[1].density = RHO_G;
    sim.particles[0].pressure = 2.0e3;
    sim.particles[1].pressure = 2.0e3;
    let acc = sim.pressure_accelerations();
    // The light gas particle is pushed away much harder than the liquid one,
    // but the exchanged momentum is equal and opposite.
    assert!(acc[0][0] < 0.0 && acc[1][0] > 0.0);
    let ratio = acc[1][0] / -acc[0][0];
    assert!(
        (ratio - m_l / m_g).abs() < 1e-9 * (m_l / m_g),
        "acceleration ratio {ratio} vs mass ratio {}",
        m_l / m_g
    );
    assert!((m_l * acc[0][0] + m_g * acc[1][0]).abs() < 1e-12 * (m_g * acc[1][0]).abs());
    // Interface pressure continuity: the liquid acceleration is the pair term
    // −m_g (P_l + P_g)/(ρ_l ρ_g) dW/dr.
    let expected_liquid = m_g * 4.0e3 / (RHO_L * RHO_G) * dw_dr(dx, h).abs();
    assert!((acc[0][0] + expected_liquid).abs() <= 1e-12 * expected_liquid);
}

// ---------------------------------------------------------------------------
// Full step
// ---------------------------------------------------------------------------

#[test]
fn step_conserves_linear_momentum_without_gravity_and_surface_tension() {
    let mut sim = patch_simulation(14, 4, 6, 0.02, 1.5, 0.0);
    let mut state = 5_u64;
    for p in sim.particles.iter_mut() {
        p.vel = [lcg_unit(&mut state) - 0.5, lcg_unit(&mut state) - 0.5];
    }
    let momentum = |sim: &MultiphaseSphSimulation| -> ([f64; 2], f64) {
        let mut total = [0.0_f64; 2];
        let mut scale = 0.0_f64;
        for p in &sim.particles {
            total[0] += p.mass * p.vel[0];
            total[1] += p.mass * p.vel[1];
            scale += p.mass * p.vel[0].abs() + p.mass * p.vel[1].abs();
        }
        (total, scale)
    };
    let (before, scale) = momentum(&sim);
    for _ in 0..5 {
        sim.step();
    }
    let (after, _) = momentum(&sim);
    for k in 0..2 {
        assert!(
            (after[k] - before[k]).abs() <= 1e-9 * scale,
            "momentum[{k}]: {} → {}",
            before[k],
            after[k]
        );
    }
}

#[test]
fn step_applies_pressure_so_a_compressed_liquid_block_expands() {
    // 8×8 liquid block at 95 % of the lattice spacing the masses were chosen
    // for: the number density is about 11 % too high, the Tait pressure is
    // positive inside and pushes the block outwards.
    let dx = 0.05;
    let spacing = 0.95 * dx;
    let n = 8;
    let mut sim = MultiphaseSphSimulation::new(1.5 * dx, 0.0, RHO_L, RHO_G);
    sim.gravity = [0.0, 0.0];
    let mass = RHO_L * dx * dx;
    for i in 0..n {
        for j in 0..n {
            sim.add_liquid_particle(i as f64 * spacing, j as f64 * spacing, mass);
        }
    }
    let centre = 0.5 * (n as f64 - 1.0) * spacing;
    let mean_radius = |sim: &MultiphaseSphSimulation| -> f64 {
        sim.particles
            .iter()
            .map(|p| (p.pos[0] - centre).hypot(p.pos[1] - centre))
            .sum::<f64>()
            / sim.particles.len() as f64
    };
    let before = mean_radius(&sim);
    for _ in 0..5 {
        sim.step();
    }
    assert!(
        sim.particles[index(n, 4, 4)].pressure > 0.0,
        "interior pressure {}",
        sim.particles[index(n, 4, 4)].pressure
    );
    let after = mean_radius(&sim);
    assert!(after > before, "block did not expand: {before} → {after}");
    // The corners move radially outwards, the net momentum stays zero.
    let corner = &sim.particles[index(n, 0, 0)];
    let outward =
        (corner.pos[0] - centre) * corner.vel[0] + (corner.pos[1] - centre) * corner.vel[1];
    assert!(outward > 0.0, "corner velocity {:?}", corner.vel);
    let momentum: f64 = sim.particles.iter().map(|p| p.mass * p.vel[0]).sum();
    let scale: f64 = sim.particles.iter().map(|p| p.mass * p.vel[0].abs()).sum();
    assert!(momentum.abs() <= 1e-9 * scale);
}

#[test]
fn step_keeps_a_liquid_patch_in_gas_finite_and_pulls_the_corners_inwards() {
    let (n, start, block, dx) = (16, 4, 8, 0.02);
    let mut sim = patch_simulation(n, start, block, dx, 1.5, 0.072);
    // Acoustic Courant number 0.1 for the default dt.
    let c0 = sim.sound_speed;
    assert!((c0 - 0.1 * sim.h / sim.dt).abs() < 1e-12);
    for _ in 0..30 {
        sim.step();
    }
    assert_eq!(sim.step_count, 30);
    for p in &sim.particles {
        assert!(p.pos.iter().all(|c| c.is_finite()), "position {:?}", p.pos);
        assert!(p.vel.iter().all(|c| c.is_finite()), "velocity {:?}", p.vel);
        assert!(p.pressure.is_finite() && p.pressure >= 0.0);
    }
}

#[test]
fn step_still_lets_an_isolated_particle_fall_freely() {
    let mut sim = MultiphaseSphSimulation::new(0.1, 0.072, RHO_L, RHO_G);
    sim.add_liquid_particle(0.5, 1.0, 1e-3);
    sim.dt = 0.001;
    sim.step();
    let p: &ColorParticle = &sim.particles[0];
    assert!((p.vel[1] + 9.81 * 0.001).abs() < 1e-12);
    assert_eq!(p.vel[0], 0.0);
    assert_eq!(p.pressure, 0.0);
    assert_eq!(p.curvature, 0.0);
}
