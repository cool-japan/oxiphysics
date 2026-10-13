// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Regression tests for the LBM boundary / inlet drive and the SPH force
//! pass of the fluid bridge.

use super::sph_forces::{NeighborGrid, compute_densities, cubic_spline_gradient, tait_pressure};
use super::*;

// ===========================================================================
// LBM: set_boundary / set_inlet_velocity must act on the lattice
// ===========================================================================

/// Kinematic viscosity 0.1 gives the BGK relaxation time tau = 0.8.
const LBM_VISCOSITY: f64 = 0.1;

fn lbm_sim(nx: usize, ny: usize, nz: usize) -> WasmLbmSimulation {
    WasmLbmSimulation::new(WasmLbmConfig::d3q19(nx, ny, nz, LBM_VISCOSITY))
}

fn lbm_index(sim: &WasmLbmSimulation, ix: usize, iy: usize, iz: usize) -> usize {
    ix + iy * sim.config.nx + iz * sim.config.nx * sim.config.ny
}

/// Mean x-velocity over the y-z plane at `ix`.
fn plane_mean_ux(sim: &WasmLbmSimulation, ix: usize) -> f64 {
    let (ny, nz) = (sim.config.ny, sim.config.nz);
    let mut sum = 0.0;
    for iz in 0..nz {
        for iy in 0..ny {
            sum += sim.velocity[lbm_index(sim, ix, iy, iz) * 3];
        }
    }
    sum / (ny * nz) as f64
}

fn lbm_max_speed(sim: &WasmLbmSimulation) -> f64 {
    sim.velocity
        .chunks_exact(3)
        .map(|u| (u[0] * u[0] + u[1] * u[1] + u[2] * u[2]).sqrt())
        .fold(0.0, f64::max)
}

#[test]
fn lbm_set_boundary_marks_the_lattice_solid() {
    let mut sim = lbm_sim(6, 4, 4);
    let idx = lbm_index(&sim, 2, 1, 3);
    sim.set_boundary(2, 1, 3, true);
    assert!(sim.boundary[idx]);
    let grid = sim.grid.as_ref().expect("grid exists after construction");
    assert!(grid.solid[idx], "the cache flag never reached the lattice");
    assert_eq!(grid.solid_count(), 1);

    sim.set_boundary(2, 1, 3, false);
    assert!(!sim.boundary[idx]);
    let grid = sim.grid.as_ref().expect("grid");
    assert_eq!(grid.solid_count(), 0);

    // Coordinates outside the grid are ignored.
    sim.set_boundary(6, 0, 0, true);
    sim.set_boundary(0, 4, 0, true);
    sim.set_boundary(0, 0, 4, true);
    assert!(sim.boundary.iter().all(|&solid| !solid));
    assert_eq!(sim.grid.as_ref().expect("grid").solid_count(), 0);
}

#[test]
fn lbm_undriven_box_stays_at_rest() {
    let mut sim = lbm_sim(8, 4, 4);
    sim.set_boundary(3, 1, 1, true);
    for _ in 0..30 {
        sim.step();
    }
    assert!(
        lbm_max_speed(&sim) < 1e-12,
        "an undriven periodic box must stay at rest, got {}",
        lbm_max_speed(&sim)
    );
}

#[test]
fn lbm_inlet_velocity_drives_the_flow() {
    // The old implementation wrote the inlet speed into the velocity cache,
    // which the next step overwrote from the (never driven) lattice.
    let mut sim = lbm_sim(16, 6, 4);
    sim.set_inlet_velocity(0.05).expect("valid inlet velocity");
    for _ in 0..60 {
        sim.step();
    }
    let inlet = plane_mean_ux(&sim, 0);
    assert!(
        (inlet - 0.05).abs() < 1e-9,
        "the inlet face must carry the prescribed velocity, got {inlet}"
    );
    let first = plane_mean_ux(&sim, 1);
    assert!(
        first > 0.02,
        "the inlet must push fluid into the domain, mean ux at x=1 is {first}"
    );
    for iz in 0..sim.config.nz {
        for iy in 0..sim.config.ny {
            let k = lbm_index(&sim, 0, iy, iz) * 3;
            assert!(sim.velocity[k + 1].abs() < 1e-9 && sim.velocity[k + 2].abs() < 1e-9);
        }
    }
}

#[test]
fn lbm_inlet_is_visible_before_the_first_step() {
    let mut sim = lbm_sim(8, 4, 4);
    sim.set_inlet_velocity(0.03).expect("valid inlet velocity");
    let inlet = plane_mean_ux(&sim, 0);
    assert!((inlet - 0.03).abs() < 1e-9, "inlet plane ux = {inlet}");
    assert!(plane_mean_ux(&sim, 4).abs() < 1e-12);
}

#[test]
fn lbm_inlet_rejects_invalid_velocity_and_keeps_state() {
    let mut sim = lbm_sim(8, 4, 4);
    for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY, 1.0, -1.0, 2.5] {
        let err = sim
            .set_inlet_velocity(bad)
            .expect_err("unusable inlet velocity must be rejected");
        assert!(
            matches!(&err, Error::InvalidParameter { name, .. } if name == "inlet_velocity"),
            "{bad}: {err:?}"
        );
        assert!(sim.inlet_velocity.is_none(), "{bad} was stored");
    }
    sim.step();
    assert!(lbm_max_speed(&sim) < 1e-12);
    assert!(
        sim.set_inlet_velocity(-0.05).is_ok(),
        "reverse flow is valid"
    );
    assert_eq!(sim.inlet_velocity, Some(-0.05));
}

#[test]
fn lbm_clear_inlet_velocity_releases_the_inlet_face() {
    // A wall across the box makes the driven flow develop a velocity gradient,
    // so a released inlet face visibly departs from the prescribed speed
    // (without it the box would keep gliding at a uniform 0.05).
    let build = || {
        let mut sim = lbm_sim(12, 4, 4);
        for iz in 0..4 {
            for iy in 0..4 {
                sim.set_boundary(6, iy, iz, true);
            }
        }
        sim.set_inlet_velocity(0.05).expect("valid inlet velocity");
        for _ in 0..20 {
            sim.step();
        }
        sim
    };
    let mut driven = build();
    let mut released = build();
    assert!((plane_mean_ux(&released, 0) - 0.05).abs() < 1e-9);
    released.clear_inlet_velocity();
    assert!(released.inlet_velocity.is_none());
    for _ in 0..5 {
        driven.step();
        released.step();
    }
    assert!(
        (plane_mean_ux(&driven, 0) - 0.05).abs() < 1e-9,
        "the driven inlet face must stay pinned"
    );
    assert!(
        (plane_mean_ux(&released, 0) - 0.05).abs() > 1e-5,
        "the inlet face is still pinned after clear_inlet_velocity"
    );
}

#[test]
fn lbm_solid_wall_blocks_the_driven_flow() {
    let (nx, ny, nz) = (12, 8, 4);
    let wall_x = 6;
    let mut open = lbm_sim(nx, ny, nz);
    let mut walled = lbm_sim(nx, ny, nz);
    for iz in 0..nz {
        for iy in 0..ny {
            walled.set_boundary(wall_x, iy, iz, true);
        }
    }
    open.set_inlet_velocity(0.05).expect("valid inlet velocity");
    walled
        .set_inlet_velocity(0.05)
        .expect("valid inlet velocity");
    for _ in 0..120 {
        open.step();
        walled.step();
    }
    // A solid node reports no velocity, however hard the flow pushes on it.
    for iz in 0..nz {
        for iy in 0..ny {
            let k = lbm_index(&walled, wall_x, iy, iz) * 3;
            assert_eq!(&walled.velocity[k..k + 3], &[0.0, 0.0, 0.0]);
        }
    }
    // The open channel carries the inlet velocity through the whole box ...
    let open_behind = plane_mean_ux(&open, wall_x + 2);
    assert!(
        (open_behind - 0.05).abs() < 1e-3,
        "open channel ux behind x=6: {open_behind}"
    );
    // ... the walled one is choked ahead of the wall (about 0.006 here) and
    // has no through-flow behind it (the periodic box is even drained
    // towards the inlet, so ux turns slightly negative, about -0.014).
    let ahead = plane_mean_ux(&walled, wall_x - 1);
    assert!(
        ahead < 0.3 * open_behind,
        "flow is not choked ahead of the wall: {ahead}"
    );
    let behind = plane_mean_ux(&walled, wall_x + 2);
    assert!(
        behind < 0.25 * open_behind,
        "wall leaks: ux behind it {behind} vs open {open_behind}"
    );
    // The inlet itself is unaffected by the wall.
    assert!((plane_mean_ux(&walled, 0) - 0.05).abs() < 1e-9);
}

#[test]
fn lbm_solid_mask_and_inlet_survive_a_serde_round_trip() {
    let mut sim = lbm_sim(8, 4, 4);
    sim.set_boundary(3, 1, 1, true);
    sim.set_inlet_velocity(0.04).expect("valid inlet velocity");
    let json = serde_json::to_string(&sim).expect("serialise");
    let mut restored: WasmLbmSimulation = serde_json::from_str(&json).expect("deserialise");
    assert!(restored.grid.is_none(), "the lattice is not serialised");
    assert_eq!(restored.inlet_velocity, Some(0.04));
    restored.step();
    let idx = lbm_index(&restored, 3, 1, 1);
    let grid = restored.grid.as_ref().expect("rebuilt by step");
    assert!(grid.solid[idx], "solid mask lost across serde");
    assert!((plane_mean_ux(&restored, 0) - 0.04).abs() < 1e-9);
}

#[test]
fn lbm_json_without_inlet_field_still_loads() {
    let sim = lbm_sim(4, 4, 4);
    let mut value: serde_json::Value = serde_json::to_value(&sim).expect("serialise");
    let removed = value
        .as_object_mut()
        .and_then(|o| o.remove("inlet_velocity"));
    assert!(removed.is_some(), "the field is serialised");
    let restored: WasmLbmSimulation = serde_json::from_value(value).expect("legacy json loads");
    assert!(restored.inlet_velocity.is_none());
}

#[test]
fn lbm_public_boundary_mask_is_resynced_every_step() {
    let mut sim = lbm_sim(8, 4, 4);
    let idx = lbm_index(&sim, 5, 2, 2);
    // Setting the public mask directly must be as effective as set_boundary.
    sim.boundary[idx] = true;
    sim.set_inlet_velocity(0.04).expect("valid inlet velocity");
    for _ in 0..10 {
        sim.step();
    }
    assert!(sim.grid.as_ref().expect("grid").solid[idx]);
    assert_eq!(&sim.velocity[idx * 3..idx * 3 + 3], &[0.0, 0.0, 0.0]);
    // ... and clearing it again frees the node.
    sim.boundary[idx] = false;
    sim.step();
    assert!(!sim.grid.as_ref().expect("grid").solid[idx]);
}

#[test]
fn lbm_short_boundary_mask_does_not_panic() {
    let mut sim = lbm_sim(6, 4, 4);
    sim.boundary.truncate(3);
    sim.set_boundary(5, 3, 3, true);
    sim.step();
    assert_eq!(sim.velocity.len(), 6 * 4 * 4 * 3);
}

// ===========================================================================
// SPH: pressure, viscosity, XSPH, surface tension and boundary handling
// ===========================================================================

/// Configuration with every optional force switched off and no gravity, so a
/// test enables exactly the term it exercises.
fn bare_sph_config(h: f64) -> WasmSphConfig {
    WasmSphConfig {
        particle_radius: 0.5 * h,
        smoothing_length: h,
        gravity: [0.0; 3],
        viscosity: 0.0,
        surface_tension: 0.0,
        xsph_enabled: false,
        ..WasmSphConfig::default()
    }
}

fn sph_with(config: WasmSphConfig, particles: Vec<WasmSphParticle>) -> WasmSphSimulation {
    let mut sim = WasmSphSimulation::new(config);
    sim.particles = particles;
    sim
}

/// Cubic lattice of `counts` particles centred on `centre`.
fn lattice(
    counts: [usize; 3],
    spacing: f64,
    centre: [f64; 3],
    mass: f64,
    boundary: bool,
) -> Vec<WasmSphParticle> {
    let mut out = Vec::new();
    let origin = [
        centre[0] - 0.5 * spacing * (counts[0] - 1) as f64,
        centre[1] - 0.5 * spacing * (counts[1] - 1) as f64,
        centre[2] - 0.5 * spacing * (counts[2] - 1) as f64,
    ];
    for iz in 0..counts[2] {
        for iy in 0..counts[1] {
            for ix in 0..counts[0] {
                let position = [
                    origin[0] + spacing * ix as f64,
                    origin[1] + spacing * iy as f64,
                    origin[2] + spacing * iz as f64,
                ];
                let id = out.len() as u64;
                out.push(if boundary {
                    WasmSphParticle::boundary(id, position, mass)
                } else {
                    WasmSphParticle::new(id, position, mass)
                });
            }
        }
    }
    out
}

// ---- kernel derivatives ----------------------------------------------------

#[test]
fn sph_kernel_gradient_matches_finite_differences() {
    let h = 0.1;
    let eps = 1e-6 * h;
    for q in [0.05, 0.2, 0.35, 0.49, 0.51, 0.7, 0.9, 0.99] {
        let r = q * h;
        let grad = cubic_spline_gradient([r, 0.0, 0.0], r, h);
        let fd =
            (cubic_spline_compact(r + eps, h) - cubic_spline_compact(r - eps, h)) / (2.0 * eps);
        assert!(
            (grad[0] - fd).abs() <= 1e-4 * fd.abs().max(1.0),
            "q={q}: analytic {} vs finite difference {fd}",
            grad[0]
        );
        assert!(grad[0] < 0.0, "the kernel decreases with distance");
        assert_eq!((grad[1], grad[2]), (0.0, 0.0));
    }
}

#[test]
fn sph_kernel_gradient_points_along_the_separation() {
    let h = 0.1;
    let d = [0.03, 0.04, 0.0];
    let r = 0.05;
    let grad = cubic_spline_gradient(d, r, h);
    let radial = cubic_spline_gradient([r, 0.0, 0.0], r, h)[0];
    // grad = W'(r) d / r, i.e. the radial derivative along the unit vector.
    assert!((grad[0] - radial * 0.6).abs() < 1e-9 * radial.abs());
    assert!((grad[1] - radial * 0.8).abs() < 1e-9 * radial.abs());
    assert_eq!(grad[2], 0.0);
}

#[test]
fn sph_kernel_gradient_vanishes_at_origin_and_outside_support() {
    let h = 0.1;
    assert_eq!(cubic_spline_gradient([0.0; 3], 0.0, h), [0.0; 3]);
    assert_eq!(cubic_spline_gradient([h, 0.0, 0.0], h, h), [0.0; 3]);
    assert_eq!(cubic_spline_gradient([0.3, 0.0, 0.0], 0.3, h), [0.0; 3]);
    assert_eq!(cubic_spline_gradient([f64::NAN; 3], f64::NAN, h), [0.0; 3]);
}

// ---- neighbour search --------------------------------------------------------

#[test]
fn sph_neighbor_grid_finds_exactly_the_particles_within_support() {
    let h = 0.1;
    let mut state: u64 = 0x9E37_79B9_7F4A_7C15;
    let mut next = move || {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        (state >> 11) as f64 / (1u64 << 53) as f64
    };
    let particles: Vec<WasmSphParticle> = (0..150)
        .map(|i| {
            let pos = [next() * 0.6 - 0.3, next() * 0.6 - 0.3, next() * 0.6 - 0.3];
            WasmSphParticle::new(i as u64, pos, 1.0)
        })
        .collect();
    let grid = NeighborGrid::build(&particles, h);
    let mut candidates_total = 0usize;
    for (i, pi) in particles.iter().enumerate() {
        let mut from_grid = Vec::new();
        grid.for_each_candidate(pi.position, |j| {
            candidates_total += 1;
            let d = [
                pi.position[0] - particles[j].position[0],
                pi.position[1] - particles[j].position[1],
                pi.position[2] - particles[j].position[2],
            ];
            if d[0] * d[0] + d[1] * d[1] + d[2] * d[2] < h * h {
                from_grid.push(j);
            }
        });
        let brute: Vec<usize> = (0..particles.len())
            .filter(|&j| {
                let d = [
                    pi.position[0] - particles[j].position[0],
                    pi.position[1] - particles[j].position[1],
                    pi.position[2] - particles[j].position[2],
                ];
                d[0] * d[0] + d[1] * d[1] + d[2] * d[2] < h * h
            })
            .collect();
        from_grid.sort_unstable();
        assert_eq!(from_grid, brute, "neighbours of particle {i}");
    }
    assert!(
        candidates_total < particles.len() * particles.len() / 2,
        "the cell list should cull most pairs ({candidates_total} candidates)"
    );
}

#[test]
fn sph_cell_list_densities_match_the_brute_force_reference() {
    let h = 0.1;
    let mut state: u64 = 0xD1B5_4A32_D192_ED03;
    let mut next = move || {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        (state >> 11) as f64 / (1u64 << 53) as f64
    };
    let particles: Vec<WasmSphParticle> = (0..120)
        .map(|i| {
            let pos = [next() * 0.4 - 0.2, next() * 0.4 - 0.2, next() * 0.4 - 0.2];
            let mass = 0.05 + 0.4 * next();
            if i % 5 == 0 {
                WasmSphParticle::boundary(i as u64, pos, mass)
            } else {
                WasmSphParticle::new(i as u64, pos, mass)
            }
        })
        .collect();
    let grid = NeighborGrid::build(&particles, h);
    let fast = compute_densities(&particles, h, &grid);
    let reference = compute_sph_densities(&particles, h);
    assert_eq!(fast.len(), reference.len());
    let mut nontrivial = 0;
    for (i, (a, b)) in fast.iter().zip(&reference).enumerate() {
        assert!(
            (a - b).abs() <= 1e-12 * b.abs().max(1.0),
            "particle {i}: cell list {a} vs brute force {b}"
        );
        // Count particles whose density exceeds their own contribution, i.e.
        // that have at least one neighbour inside the support.
        if *b > 1.0001 * particles[i].mass * cubic_spline_compact(0.0, h) {
            nontrivial += 1;
        }
    }
    assert!(
        nontrivial > 20,
        "too few particles have neighbours: {nontrivial}"
    );
}

#[test]
fn sph_densities_ignore_particles_with_non_finite_coordinates() {
    let h = 0.1;
    let particles = vec![
        WasmSphParticle::new(0, [f64::NAN, 0.0, 0.0], 1.0),
        WasmSphParticle::new(1, [0.0, 0.0, 0.0], 1.0),
        WasmSphParticle::new(2, [0.05, 0.0, 0.0], 1.0),
        WasmSphParticle::new(3, [f64::INFINITY, 0.0, 0.0], 1.0),
    ];
    let grid = NeighborGrid::build(&particles, h);
    let rho = compute_densities(&particles, h, &grid);
    let expected = cubic_spline_compact(0.0, h) + cubic_spline_compact(0.05, h);
    assert!((rho[1] - expected).abs() < 1e-9 * expected, "{}", rho[1]);
    assert!((rho[2] - expected).abs() < 1e-9 * expected, "{}", rho[2]);
}

#[test]
fn sph_neighbor_grid_tolerates_extreme_and_nan_positions() {
    let h = 0.1;
    let mut particles = vec![
        WasmSphParticle::new(0, [f64::NAN, 0.0, 0.0], 1.0),
        WasmSphParticle::new(1, [1.0e300, -1.0e300, f64::INFINITY], 1.0),
        WasmSphParticle::new(2, [0.0; 3], 1.0),
    ];
    let grid = NeighborGrid::build(&particles, h);
    let mut visited = 0;
    for p in &particles {
        grid.for_each_candidate(p.position, |_| visited += 1);
    }
    assert!(visited > 0);
    // A full step must not panic either.
    particles[0].position = [0.0, 0.0, 0.0];
    let mut sim = sph_with(bare_sph_config(h), particles);
    sim.step(1e-3);
    assert_eq!(sim.step_count, 1);
}

// ---- equation of state -------------------------------------------------------

#[test]
fn sph_tait_pressure_is_clamped_at_zero() {
    let cfg = WasmSphConfig::default();
    assert_eq!(tait_pressure(cfg.rest_density, &cfg), 0.0);
    assert_eq!(tait_pressure(0.5 * cfg.rest_density, &cfg), 0.0);
    assert_eq!(tait_pressure(0.0, &cfg), 0.0);
    assert_eq!(tait_pressure(f64::NAN, &cfg), 0.0);
    let compressed = tait_pressure(1.1 * cfg.rest_density, &cfg);
    let expected = cfg.stiffness * (1.1_f64.powf(cfg.tait_exponent) - 1.0);
    assert!(
        (compressed - expected).abs() < 1e-9 * expected,
        "{compressed}"
    );
    assert!(tait_pressure(1.2 * cfg.rest_density, &cfg) > compressed);
}

// ---- pressure force ----------------------------------------------------------

#[test]
fn sph_overlapping_particles_repel_and_conserve_momentum() {
    // Two particles at 0.3 h compress each other well above the rest density.
    // The old step computed the pressure and then ignored it, so nothing
    // happened.
    let h = 0.1;
    let mass = 0.4;
    let mut sim = sph_with(
        bare_sph_config(h),
        vec![
            WasmSphParticle::new(0, [0.0, 0.0, 0.0], mass),
            WasmSphParticle::new(1, [0.3 * h, 0.0, 0.0], mass),
        ],
    );
    sim.step(1e-3);
    let (a, b) = (&sim.particles[0], &sim.particles[1]);
    assert!(a.density > sim.config.rest_density, "density {}", a.density);
    assert!(a.pressure > 0.0);
    assert!(
        a.velocity[0] < -0.1,
        "particle 0 must be pushed to -x: {:?}",
        a.velocity
    );
    assert!(
        b.velocity[0] > 0.1,
        "particle 1 must be pushed to +x: {:?}",
        b.velocity
    );
    assert!(
        (a.mass * a.velocity[0] + b.mass * b.velocity[0]).abs() < 1e-9 * b.velocity[0].abs(),
        "pressure forces must conserve momentum"
    );
    for p in [a, b] {
        assert_eq!(p.velocity[1], 0.0);
        assert_eq!(p.velocity[2], 0.0);
    }
    assert!(
        b.position[0] - a.position[0] > 0.3 * h,
        "the pair moved apart"
    );
}

#[test]
fn sph_particles_beyond_the_support_do_not_interact() {
    let h = 0.1;
    let mut sim = sph_with(
        WasmSphConfig {
            smoothing_length: h,
            particle_radius: 0.5 * h,
            ..WasmSphConfig::default()
        },
        vec![
            WasmSphParticle::new(0, [0.0, 0.0, 0.0], 0.4),
            WasmSphParticle::new(1, [2.0 * h, 0.0, 0.0], 0.4),
        ],
    );
    let dt = 1e-3;
    sim.step(dt);
    for p in &sim.particles {
        assert!(
            (p.velocity[1] + 9.81 * dt).abs() < 1e-12,
            "{:?}",
            p.velocity
        );
        assert_eq!(p.velocity[0], 0.0);
        assert_eq!(p.velocity[2], 0.0);
    }
}

// ---- boundary particles ------------------------------------------------------

#[test]
fn sph_boundary_particles_repel_fluid_but_never_move() {
    let h = 0.1;
    let mass = 0.4;
    let mut sim = sph_with(
        bare_sph_config(h),
        vec![
            WasmSphParticle::boundary(0, [0.0, 0.0, 0.0], mass),
            WasmSphParticle::new(1, [0.0, 0.4 * h, 0.0], mass),
        ],
    );
    sim.step(1e-3);
    let (wall, fluid) = (&sim.particles[0], &sim.particles[1]);
    assert!(
        wall.density > 0.0 && wall.pressure > 0.0,
        "boundary density/pressure are computed"
    );
    assert!(
        fluid.velocity[1] > 0.1,
        "fluid must be pushed off the wall: {:?}",
        fluid.velocity
    );
    assert_eq!(wall.position, [0.0, 0.0, 0.0]);
    assert_eq!(wall.velocity, [0.0; 3]);
    assert_eq!(wall.acceleration, [0.0; 3]);
}

#[test]
fn sph_boundary_particle_api_adds_static_particles() {
    let mut sim = WasmSphSimulation::new(WasmSphConfig::default());
    let fluid = sim.add_particle([0.0, 1.0, 0.0], 0.1);
    let wall = sim.add_boundary_particle([0.0, 0.0, 0.0], 0.1);
    assert_eq!((fluid, wall), (0, 1));
    sim.add_boundary_particles(&[1.0, 0.0, 0.0, 2.0, 0.0, 0.0, 9.0], 0.1);
    assert_eq!(
        sim.particles.len(),
        4,
        "the trailing incomplete triple is ignored"
    );
    assert_eq!(sim.fluid_count(), 1);
    assert!(sim.particles[1..].iter().all(|p| p.is_boundary));
    let ids: Vec<u64> = sim.particles.iter().map(|p| p.id).collect();
    assert_eq!(ids, vec![0, 1, 2, 3], "ids stay unique");
    assert_eq!(sim.add_boundary_particle_js(5.0, 0.0, 0.0, 0.1), 4);
    sim.add_boundary_particles_js(&[6.0, 0.0, 0.0], 0.1);
    assert_eq!(sim.particles.len(), 6);
}

/// Outcome of letting a column of water settle on a two-layer boundary floor.
struct SettledColumn {
    /// Lowest fluid particle (the top layer of the floor sits at `y = 0`).
    y_min: f64,
    /// Highest fluid particle.
    y_max: f64,
    /// Fastest fluid particle.
    v_max: f64,
    /// Whether every coordinate of every particle is finite.
    finite: bool,
    /// The floor particles are exactly where they started, at rest.
    floor_untouched: bool,
    /// At least one boundary particle has a computed (non-default) density.
    boundary_density_computed: bool,
}

/// A 3x3x6 column of water (lattice spacing 0.05 m) on a two-layer boundary
/// floor, advanced for 0.15 s with `particle_radius` as the contact distance.
fn settle_column(particle_radius: f64) -> SettledColumn {
    let dx = 0.05;
    let mass = 1000.0 * dx * dx * dx;
    let cfg = WasmSphConfig {
        particle_radius,
        smoothing_length: 2.0 * dx,
        stiffness: 1.0e5,
        ..WasmSphConfig::default()
    };
    let mut sim = WasmSphSimulation::new(cfg);
    for layer in 0..2 {
        for ix in -4..=4 {
            for iz in -4..=4 {
                sim.add_boundary_particle(
                    [ix as f64 * dx, -(layer as f64) * dx, iz as f64 * dx],
                    mass,
                );
            }
        }
    }
    for iy in 1..=6 {
        for ix in -1..=1 {
            for iz in -1..=1 {
                sim.add_particle([ix as f64 * dx, iy as f64 * dx, iz as f64 * dx], mass);
            }
        }
    }
    let floor: Vec<[f64; 3]> = sim
        .particles
        .iter()
        .filter(|p| p.is_boundary)
        .map(|p| p.position)
        .collect();
    for _ in 0..500 {
        sim.step(3.0e-4);
    }
    let fluid: Vec<&WasmSphParticle> = sim.particles.iter().filter(|p| !p.is_boundary).collect();
    let walls: Vec<&WasmSphParticle> = sim.particles.iter().filter(|p| p.is_boundary).collect();
    SettledColumn {
        y_min: fluid
            .iter()
            .map(|p| p.position[1])
            .fold(f64::INFINITY, f64::min),
        y_max: fluid
            .iter()
            .map(|p| p.position[1])
            .fold(f64::NEG_INFINITY, f64::max),
        v_max: fluid.iter().map(|p| p.speed()).fold(0.0, f64::max),
        finite: sim.particles.iter().all(|p| {
            p.position
                .iter()
                .chain(p.velocity.iter())
                .all(|c| c.is_finite())
        }),
        floor_untouched: walls.len() == floor.len()
            && walls
                .iter()
                .zip(&floor)
                .all(|(wall, start)| &wall.position == start && wall.velocity == [0.0; 3]),
        boundary_density_computed: walls.iter().any(|w| (w.density - 1000.0).abs() > 1.0),
    }
}

#[test]
fn sph_column_rests_on_the_boundary_layer() {
    // Without the pressure force the column free-falls: after 0.15 s its
    // lowest particle is at y = -0.06 m, moving at 1.5 m/s.  With it, the
    // particles come to rest on the floor at about y = 0.044 m (the usual gap
    // of dynamic boundary particles) and keep moving slowly.
    let dx = 0.05;
    let column = settle_column(0.5 * dx);
    assert!(
        column.y_min > 0.5 * dx,
        "fluid sank into the boundary layer: y_min = {}",
        column.y_min
    );
    assert!(
        column.y_max - column.y_min > 2.5 * dx,
        "column collapsed: {}..{}",
        column.y_min,
        column.y_max
    );
    assert!(
        column.v_max < 1.0,
        "column is still falling at {} m/s",
        column.v_max
    );
    assert!(column.finite);
    // The boundary layer neither moved nor picked up velocity, and its
    // densities were computed (they used to stay at the 1000 default).
    assert!(column.floor_untouched);
    assert!(column.boundary_density_computed);
}

#[test]
fn sph_column_with_the_default_particle_radius_settles_smoothly() {
    // The default particle_radius (0.05 m) is as large as the resting gap
    // above the floor, so the contact response is active while the column
    // settles.  It must act as a clean contact layer (y >= ~0.05), not make
    // the column jitter or sink.
    let column = settle_column(WasmSphConfig::default().particle_radius);
    assert!(
        column.y_min > 0.045,
        "fluid sank below the contact layer: y_min = {}",
        column.y_min
    );
    assert!(
        column.y_max - column.y_min > 0.125,
        "column collapsed: {}..{}",
        column.y_min,
        column.y_max
    );
    assert!(
        column.v_max < 1.0,
        "column is not settling: {} m/s",
        column.v_max
    );
    assert!(column.finite && column.floor_untouched);
}

// ---- viscosity and XSPH --------------------------------------------------------

/// Two close particles moving in opposite y-directions, no pressure.
fn shearing_pair(config: WasmSphConfig) -> WasmSphSimulation {
    let h = config.smoothing_length;
    let mut a = WasmSphParticle::new(0, [0.0, 0.0, 0.0], 0.1);
    let mut b = WasmSphParticle::new(1, [0.5 * h, 0.0, 0.0], 0.1);
    a.velocity = [0.0, 1.0, 0.0];
    b.velocity = [0.0, -1.0, 0.0];
    sph_with(config, vec![a, b])
}

#[test]
fn sph_viscosity_damps_relative_motion_and_conserves_momentum() {
    let h = 0.1;
    let inviscid = WasmSphConfig {
        stiffness: 0.0,
        ..bare_sph_config(h)
    };
    let viscous = WasmSphConfig {
        viscosity: 5.0,
        ..inviscid.clone()
    };
    let mut free = shearing_pair(inviscid);
    free.step(1e-4);
    assert_eq!(free.particles[0].velocity[1], 1.0);
    assert_eq!(free.particles[1].velocity[1], -1.0);

    let mut sim = shearing_pair(viscous);
    sim.step(1e-4);
    let (va, vb) = (sim.particles[0].velocity[1], sim.particles[1].velocity[1]);
    assert!(va < 1.0 && va > 0.0, "particle 0 slowed to {va}");
    assert!(vb > -1.0 && vb < 0.0, "particle 1 slowed to {vb}");
    assert!((va + vb).abs() < 1e-12, "momentum: {va} + {vb}");
    assert!(va - vb < 2.0, "relative speed must shrink");
}

#[test]
fn sph_xsph_blends_neighbour_velocities_and_conserves_momentum() {
    let h = 0.1;
    let base = WasmSphConfig {
        stiffness: 0.0,
        ..bare_sph_config(h)
    };
    let mut off = shearing_pair(base.clone());
    off.step(1e-6);
    assert_eq!(
        off.particles[0].velocity[1], 1.0,
        "XSPH disabled: untouched"
    );

    let smoothed = WasmSphConfig {
        xsph_enabled: true,
        xsph_epsilon: 0.5,
        ..base
    };
    let mut sim = shearing_pair(smoothed);
    sim.step(1e-6);
    let (va, vb) = (sim.particles[0].velocity[1], sim.particles[1].velocity[1]);
    // Both particles have rho = m (W(0) + W(r)), so the correction is
    // eps * m / rho * W(r) * (v_j - v_i).
    let w0 = cubic_spline_compact(0.0, h);
    let wr = cubic_spline_compact(0.5 * h, h);
    let expected = 1.0 + 0.5 * (0.1 / (0.1 * (w0 + wr))) * wr * (-2.0);
    assert!((va - expected).abs() < 1e-9, "{va} vs {expected}");
    assert!((vb + expected).abs() < 1e-9);
    assert!((va + vb).abs() < 1e-12, "XSPH must conserve momentum");
    assert!(va < 1.0 && va > -1.0);
}

// ---- surface tension ---------------------------------------------------------

/// Droplet used by the surface-tension tests: a 5x5x5 lattice cube.
fn cube_droplet(dx: f64, mass: f64) -> Vec<WasmSphParticle> {
    lattice([5, 5, 5], dx, [0.0; 3], mass, false)
}

/// Extra acceleration of every particle caused by surface tension `gamma`:
/// the difference of two otherwise identical single steps.
fn surface_tension_accelerations(gamma: f64) -> (Vec<WasmSphParticle>, Vec<[f64; 3]>, Vec<f64>) {
    let dx = 0.05;
    let h = 2.0 * dx;
    let mass = 1000.0 * dx * dx * dx;
    let calm = WasmSphConfig {
        stiffness: 0.0,
        ..bare_sph_config(h)
    };
    let tense = WasmSphConfig {
        surface_tension: gamma,
        ..calm.clone()
    };
    let mut reference = sph_with(calm, cube_droplet(dx, mass));
    let mut sim = sph_with(tense, cube_droplet(dx, mass));
    reference.step(1e-9);
    sim.step(1e-9);
    let extra: Vec<[f64; 3]> = sim
        .particles
        .iter()
        .zip(&reference.particles)
        .map(|(p, q)| {
            [
                p.acceleration[0] - q.acceleration[0],
                p.acceleration[1] - q.acceleration[1],
                p.acceleration[2] - q.acceleration[2],
            ]
        })
        .collect();
    let colour = sim.particles.iter().map(|p| p.colour).collect();
    (sim.particles, extra, colour)
}

#[test]
fn sph_surface_tension_pulls_convex_corners_and_edges_inward() {
    let dx = 0.05;
    let (particles, extra, colour) = surface_tension_accelerations(2000.0);
    let norm = |v: [f64; 3]| (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
    let radial = |i: usize| {
        let (p, e) = (&particles[i].position, &extra[i]);
        p[0] * e[0] + p[1] * e[1] + p[2] * e[2]
    };
    let lattice_coords = |i: usize| -> [i32; 3] {
        let p = &particles[i].position;
        [
            (p[0] / dx).round() as i32,
            (p[1] / dx).round() as i32,
            (p[2] / dx).round() as i32,
        ]
    };
    let centre = 62; // (2, 2, 2) of the 5x5x5 lattice
    assert_eq!(particles[centre].position, [0.0; 3]);
    assert_eq!(
        extra[centre], [0.0; 3],
        "an interior particle feels no surface tension"
    );

    let mut corner_strength = 0.0_f64;
    let mut strongest_face_centre = 0.0_f64;
    let mut inward_edges_and_corners = 0;
    let mut net_radial = 0.0;
    for (i, e) in extra.iter().enumerate() {
        let k = lattice_coords(i);
        let on_shell = k.iter().filter(|c| c.abs() == 2).count();
        let on_axis = k.iter().filter(|c| **c == 0).count();
        net_radial += radial(i);
        match (on_shell, on_axis) {
            // corner of the cube
            (3, _) => {
                corner_strength = corner_strength.max(norm(*e));
                assert!(radial(i) < 0.0, "corner {k:?} pushed outward");
                inward_edges_and_corners += 1;
            }
            // edge of the cube
            (2, _) => {
                assert!(radial(i) < 0.0, "edge particle {k:?} pushed outward");
                inward_edges_and_corners += 1;
            }
            // centre of a flat face
            (1, 2) => {
                strongest_face_centre = strongest_face_centre.max(norm(*e));
            }
            _ => {}
        }
    }
    assert_eq!(inward_edges_and_corners, 8 + 36);
    assert!(corner_strength > 100.0, "corner pull {corner_strength}");
    assert!(
        strongest_face_centre < 0.05 * corner_strength,
        "a flat face must feel (almost) no surface tension: {strongest_face_centre} vs {corner_strength}"
    );
    assert!(net_radial < 0.0, "the droplet is pulled inward overall");

    // The force is linear in the surface tension coefficient.
    let (_, doubled, _) = surface_tension_accelerations(4000.0);
    for (single, twice) in extra.iter().zip(&doubled) {
        for k in 0..3 {
            assert!((2.0 * single[k] - twice[k]).abs() <= 1e-9 * twice[k].abs().max(1.0));
        }
    }
    // The colour field is 1 in the bulk and stays positive and bounded.
    assert!(
        (colour[centre] - 1.0).abs() < 1e-9,
        "bulk colour {}",
        colour[centre]
    );
    assert!(colour.iter().all(|c| *c > 0.5 && *c < 1.5));
}

#[test]
fn sph_zero_surface_tension_adds_no_force() {
    let dx = 0.05;
    let h = 2.0 * dx;
    let mass = 1000.0 * dx * dx * dx;
    // Without pressure, gravity or viscosity the only possible acceleration
    // is the surface-tension force.
    let accelerations = |gamma: f64| {
        let cfg = WasmSphConfig {
            surface_tension: gamma,
            stiffness: 0.0,
            ..bare_sph_config(h)
        };
        let mut sim = sph_with(cfg, lattice([3, 3, 3], dx, [0.0; 3], mass, false));
        sim.step(1e-6);
        sim.particles
            .iter()
            .map(|p| p.acceleration)
            .collect::<Vec<_>>()
    };
    assert!(accelerations(0.0).iter().all(|a| *a == [0.0; 3]));
    assert!(
        accelerations(500.0).iter().any(|a| *a != [0.0; 3]),
        "surface tension had no effect at all"
    );
}

// ---- boundary restitution ------------------------------------------------------

/// A fluid particle falling straight onto a boundary particle, with every
/// other force disabled so only the contact response acts.
fn contact_sim(restitution: f64, fluid_velocity: [f64; 3]) -> WasmSphSimulation {
    let cfg = WasmSphConfig {
        particle_radius: 0.02,
        smoothing_length: 0.1,
        stiffness: 0.0,
        boundary_restitution: restitution,
        ..bare_sph_config(0.1)
    };
    let mut fluid = WasmSphParticle::new(1, [0.0, 0.05, 0.0], 0.1);
    fluid.velocity = fluid_velocity;
    sph_with(
        cfg,
        vec![WasmSphParticle::boundary(0, [0.0; 3], 0.1), fluid],
    )
}

#[test]
fn sph_boundary_contact_reflects_with_the_configured_restitution() {
    for restitution in [0.0, 0.3, 1.0] {
        let mut sim = contact_sim(restitution, [0.0, -5.0, 0.0]);
        sim.step(0.008);
        let fluid = &sim.particles[1];
        assert!(
            (fluid.position[1] - 0.02).abs() < 1e-12,
            "e={restitution}: pushed back to the contact distance, got y = {}",
            fluid.position[1]
        );
        assert!(
            (fluid.velocity[1] - restitution * 5.0).abs() < 1e-12,
            "e={restitution}: vy = {}",
            fluid.velocity[1]
        );
        assert_eq!((fluid.velocity[0], fluid.velocity[2]), (0.0, 0.0));
        assert_eq!(sim.particles[0].position, [0.0; 3]);
        assert_eq!(sim.particles[0].velocity, [0.0; 3]);
    }
}

#[test]
fn sph_boundary_contact_oblique_impact_reflects_only_the_normal_component() {
    let restitution = 0.5;
    let velocity = [1.0, -5.0, 0.0];
    let dt = 0.008;
    let mut sim = contact_sim(restitution, velocity);
    sim.step(dt);
    // Position after the free flight: (0.008, 0.01, 0), inside the contact
    // distance of 0.02 around the boundary particle at the origin.
    let flown = [velocity[0] * dt, 0.05 + velocity[1] * dt, 0.0];
    let r = (flown[0] * flown[0] + flown[1] * flown[1]).sqrt();
    let n = [flown[0] / r, flown[1] / r, 0.0];
    let v_n = velocity[0] * n[0] + velocity[1] * n[1];
    assert!(v_n < 0.0);
    let fluid = &sim.particles[1];
    for k in 0..3 {
        let expected_v = velocity[k] - (1.0 + restitution) * v_n * n[k];
        assert!(
            (fluid.velocity[k] - expected_v).abs() < 1e-12,
            "component {k}: {} vs {expected_v}",
            fluid.velocity[k]
        );
        assert!(
            (fluid.position[k] - 0.02 * n[k]).abs() < 1e-12,
            "component {k}: position {} vs {}",
            fluid.position[k],
            0.02 * n[k]
        );
    }
    // The tangential component is untouched, the normal one is scaled by -e.
    let after_n = fluid.velocity[0] * n[0] + fluid.velocity[1] * n[1];
    assert!((after_n + restitution * v_n).abs() < 1e-12);
}

#[test]
fn sph_boundary_contact_ignores_separating_and_distant_particles() {
    // Moving away from the wall: only the position is corrected.
    let mut sim = contact_sim(0.5, [0.0, 1.0, 0.0]);
    sim.particles[1].position = [0.0, 0.01, 0.0];
    sim.step(1e-3);
    assert_eq!(sim.particles[1].velocity[1], 1.0);
    assert!((sim.particles[1].position[1] - 0.02).abs() < 1e-12);
    // Far from every boundary particle nothing happens.
    let mut sim = contact_sim(0.5, [0.0, -1.0, 0.0]);
    sim.step(1e-3);
    assert!((sim.particles[1].position[1] - (0.05 - 1e-3)).abs() < 1e-12);
    assert_eq!(sim.particles[1].velocity[1], -1.0);
}

#[test]
fn sph_boundary_contact_handles_coincident_centres_and_odd_restitution() {
    // The particle flies exactly onto the boundary particle (2e-4 - 2 * 1e-4
    // is exactly zero), so the contact normal comes from the reversed velocity.
    let mut sim = contact_sim(f64::NAN, [0.0, -2.0, 0.0]);
    sim.particles[1].position = [0.0, 2.0e-4, 0.0];
    sim.step(1.0e-4);
    let fluid = &sim.particles[1];
    assert!(fluid.position.iter().all(|c| c.is_finite()));
    assert!(
        (fluid.position[1] - 0.02).abs() < 1e-12,
        "{:?}",
        fluid.position
    );
    // NaN restitution is treated as 0 (a perfectly inelastic contact) ...
    assert!(fluid.velocity[1].abs() < 1e-12, "{:?}", fluid.velocity);
    // ... a particle resting exactly on the boundary particle is pushed up ...
    let mut sim = contact_sim(0.5, [0.0; 3]);
    sim.particles[1].position = [0.0; 3];
    sim.step(1.0e-4);
    assert!((sim.particles[1].position[1] - 0.02).abs() < 1e-12);
    assert_eq!(sim.particles[1].velocity, [0.0; 3]);
    // ... and restitution above 1 is clamped (no energy gain).
    let mut sim = contact_sim(3.0, [0.0, -5.0, 0.0]);
    sim.step(0.008);
    assert!((sim.particles[1].velocity[1] - 5.0).abs() < 1e-12);
}

// ---- robustness ----------------------------------------------------------------

#[test]
fn sph_invalid_kernel_config_falls_back_to_gravity_only() {
    for (h, rho0) in [
        (0.0, 1000.0),
        (-0.1, 1000.0),
        (f64::NAN, 1000.0),
        (0.1, 0.0),
    ] {
        let cfg = WasmSphConfig {
            smoothing_length: h,
            rest_density: rho0,
            ..WasmSphConfig::default()
        };
        let mut sim = sph_with(
            cfg,
            vec![
                WasmSphParticle::new(0, [0.0, 1.0, 0.0], 0.1),
                WasmSphParticle::boundary(1, [0.0, 0.0, 0.0], 0.1),
            ],
        );
        sim.step(0.01);
        let p = &sim.particles[0];
        assert!(
            (p.velocity[1] + 9.81 * 0.01).abs() < 1e-12,
            "h={h} rho0={rho0}: {:?}",
            p.velocity
        );
        assert!(p.position.iter().all(|c| c.is_finite()));
        assert_eq!(sim.particles[1].position, [0.0; 3]);
    }
}

#[test]
fn sph_zero_mass_particles_do_not_poison_the_simulation() {
    let h = 0.1;
    let mut sim = sph_with(
        WasmSphConfig {
            smoothing_length: h,
            particle_radius: 0.5 * h,
            ..WasmSphConfig::default()
        },
        vec![
            // Isolated and massless: its density is exactly zero.
            WasmSphParticle::new(0, [5.0, 0.0, 0.0], 0.0),
            WasmSphParticle::new(1, [0.0, 0.0, 0.0], 0.4),
            WasmSphParticle::new(2, [0.02, 0.0, 0.0], 0.4),
            // Massless but embedded among the others.
            WasmSphParticle::new(3, [0.01, 0.0, 0.0], 0.0),
        ],
    );
    for _ in 0..5 {
        sim.step(1e-3);
    }
    for p in &sim.particles {
        assert!(
            p.position
                .iter()
                .chain(p.velocity.iter())
                .all(|c| c.is_finite()),
            "{p:?}"
        );
    }
    let isolated = &sim.particles[0];
    assert!(
        (isolated.velocity[1] + 9.81 * 5e-3).abs() < 1e-12,
        "{:?}",
        isolated.velocity
    );
}

#[test]
fn sph_step_is_deterministic() {
    let dx = 0.05;
    let mass = 1000.0 * dx * dx * dx;
    let build = || {
        let cfg = WasmSphConfig {
            particle_radius: 0.5 * dx,
            smoothing_length: 2.0 * dx,
            stiffness: 1.0e5,
            xsph_enabled: true,
            ..WasmSphConfig::default()
        };
        let mut particles = lattice([5, 1, 5], dx, [0.0, -dx, 0.0], mass, true);
        particles.extend(lattice([3, 3, 3], dx, [0.0, 2.0 * dx, 0.0], mass, false));
        sph_with(cfg, particles)
    };
    let (mut a, mut b) = (build(), build());
    for _ in 0..25 {
        a.step(2e-4);
        b.step(2e-4);
    }
    for (p, q) in a.particles.iter().zip(&b.particles) {
        for k in 0..3 {
            assert_eq!(p.position[k].to_bits(), q.position[k].to_bits());
            assert_eq!(p.velocity[k].to_bits(), q.velocity[k].to_bits());
        }
    }
}
