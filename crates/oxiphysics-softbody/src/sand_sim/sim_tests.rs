// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! API-level tests of [`SandSimulation`]: parameter validation, sub-stepping,
//! inert particles, state ownership, boundaries and conservation laws.

use super::mpm::{Material, mat_from_col_major};
use super::test_support::{Lcg, assert_close, assert_vec_close, max_abs_diff, slab_column};
use super::{MAX_SUBSTEPS_PER_STEP, SandParticle, SandSimulation};
use crate::Error;

fn block(grid_size: f64, spacing: f64) -> SandSimulation {
    let mut sim = SandSimulation::new(grid_size, [0.0, -9.81, 0.0]);
    sim.cfl = 0.6;
    let filled = sim.fill_box([-0.1, 0.0, -0.1], [0.1, 0.1, 0.1], spacing, 1600.0);
    assert_eq!(
        filled.ok(),
        Some(((0.2f64 / spacing).round() as usize).pow(2) * (0.1f64 / spacing).round() as usize)
    );
    sim
}

/// Bit-exact equality of particle states (NaN-safe).
fn bits(p: &SandParticle) -> Vec<u64> {
    let mut v: Vec<u64> = p.position.iter().map(|x| x.to_bits()).collect();
    v.extend(p.velocity.iter().map(|x| x.to_bits()));
    v.push(p.mass.to_bits());
    v.push(p.volume.to_bits());
    v.extend(p.deformation_gradient.iter().map(|x| x.to_bits()));
    v.push(p.plastic_strain_vol.to_bits());
    v.extend(p.affine_velocity.iter().map(|x| x.to_bits()));
    v
}

fn assert_same_states(a: &[SandParticle], b: &[SandParticle], what: &str) {
    assert_eq!(a.len(), b.len(), "{what}: length");
    for (i, (p, q)) in a.iter().zip(b).enumerate() {
        assert_eq!(bits(p), bits(q), "{what}: particle {i}");
    }
}

// ---------------------------------------------------------------------------
// Construction and parameters
// ---------------------------------------------------------------------------

#[test]
fn new_simulation_has_the_documented_defaults() {
    let sim = SandSimulation::new(0.05, [0.0, -9.81, 0.0]);
    assert!(sim.particles.is_empty());
    assert_eq!(sim.grid_size, 0.05);
    assert_eq!(sim.gravity, [0.0, -9.81, 0.0]);
    assert_eq!(sim.youngs_modulus_pa, 3.537e5);
    assert_eq!(sim.poisson_ratio, 0.3);
    assert_eq!(sim.friction_angle_deg, 30.0);
    assert_eq!(sim.floor_friction, 1.0);
    assert_eq!(sim.cfl, 0.3);
    assert_eq!(sim.plane_strain_z, None);
    assert_eq!(sim.total_mass(), 0.0);
    assert_eq!(sim.kinetic_energy(), 0.0);
    assert_eq!(sim.total_momentum(), [0.0; 3]);
    assert!(sim.stable_time_step().is_ok_and(f64::is_infinite));
}

#[test]
fn invalid_simulation_parameters_are_rejected_and_change_nothing() {
    type Mutation = Box<dyn Fn(&mut SandSimulation)>;
    let cases: Vec<(&str, Mutation)> = vec![
        ("zero grid", Box::new(|s| s.grid_size = 0.0)),
        ("negative grid", Box::new(|s| s.grid_size = -0.05)),
        ("NaN grid", Box::new(|s| s.grid_size = f64::NAN)),
        ("infinite grid", Box::new(|s| s.grid_size = f64::INFINITY)),
        ("subnormal grid", Box::new(|s| s.grid_size = 1.0e-320)),
        ("NaN gravity", Box::new(|s| s.gravity[1] = f64::NAN)),
        (
            "infinite gravity",
            Box::new(|s| s.gravity[0] = f64::INFINITY),
        ),
        ("zero E", Box::new(|s| s.youngs_modulus_pa = 0.0)),
        ("NaN E", Box::new(|s| s.youngs_modulus_pa = f64::NAN)),
        ("nu = 0.5", Box::new(|s| s.poisson_ratio = 0.5)),
        ("negative nu", Box::new(|s| s.poisson_ratio = -0.1)),
        ("phi = 90", Box::new(|s| s.friction_angle_deg = 90.0)),
        ("negative phi", Box::new(|s| s.friction_angle_deg = -5.0)),
        ("NaN phi", Box::new(|s| s.friction_angle_deg = f64::NAN)),
        (
            "negative floor friction",
            Box::new(|s| s.floor_friction = -0.1),
        ),
        (
            "NaN floor friction",
            Box::new(|s| s.floor_friction = f64::NAN),
        ),
        ("zero cfl", Box::new(|s| s.cfl = 0.0)),
        ("cfl above one", Box::new(|s| s.cfl = 1.5)),
        ("NaN cfl", Box::new(|s| s.cfl = f64::NAN)),
        (
            "reversed slab",
            Box::new(|s| s.plane_strain_z = Some([0.1, 0.0])),
        ),
        (
            "empty slab",
            Box::new(|s| s.plane_strain_z = Some([0.1, 0.1])),
        ),
        (
            "NaN slab",
            Box::new(|s| s.plane_strain_z = Some([f64::NAN, 0.1])),
        ),
    ];
    for (what, mutate) in cases {
        let mut sim = SandSimulation::new(0.05, [0.0, -9.81, 0.0]);
        sim.add_particle(SandParticle::new([0.0, 0.5, 0.0], 1.0e-3, 1.0e-6));
        mutate(&mut sim);
        let before = sim.particles.clone();
        assert!(
            matches!(sim.try_step(0.01), Err(Error::InvalidParameter(_))),
            "{what}: try_step must fail"
        );
        assert!(
            matches!(sim.stable_time_step(), Err(Error::InvalidParameter(_))),
            "{what}: stable_time_step must fail"
        );
        sim.step(0.01);
        assert_same_states(&sim.particles, &before, what);
    }
    // A sticky floor (infinite friction) is valid.
    let mut sticky = SandSimulation::new(0.05, [0.0, -9.81, 0.0]);
    sticky.floor_friction = f64::INFINITY;
    assert!(sticky.try_step(0.01).is_ok());
    // So are the boundary values of the material ranges.
    let mut edge = SandSimulation::new(0.05, [0.0; 3]);
    edge.friction_angle_deg = 0.0;
    edge.poisson_ratio = 0.0;
    edge.cfl = 1.0;
    assert!(edge.try_step(0.01).is_ok());
}

#[test]
fn invalid_time_steps_are_rejected_and_zero_is_a_no_op() {
    let mut sim = block(0.05, 0.025);
    let before = sim.particles.clone();
    for dt in [-0.01, f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert!(
            matches!(sim.try_step(dt), Err(Error::InvalidParameter(_))),
            "dt={dt}"
        );
        sim.step(dt);
        assert_same_states(&sim.particles, &before, "invalid dt");
    }
    assert!(sim.try_step(0.0).is_ok());
    assert_same_states(&sim.particles, &before, "zero dt");
}

// ---------------------------------------------------------------------------
// fill_box
// ---------------------------------------------------------------------------

#[test]
fn fill_box_builds_a_regular_lattice_with_consistent_mass_and_volume() {
    let mut sim = SandSimulation::new(0.05, [0.0, -9.81, 0.0]);
    let n = sim.fill_box([0.0, 0.0, 0.0], [0.1, 0.05, 0.075], 0.025, 1600.0);
    assert_eq!(n.ok(), Some(4 * 2 * 3));
    assert_eq!(sim.particles.len(), 24);
    for p in &sim.particles {
        assert_close(p.volume, 0.025f64.powi(3), 1e-12, "cell volume");
        assert_close(p.mass, 1600.0 * p.volume, 1e-12, "mass");
        assert!(p.position[0] > 0.0 && p.position[0] < 0.1);
        assert!(p.position[1] > 0.0 && p.position[1] < 0.05);
        assert!(p.position[2] > 0.0 && p.position[2] < 0.075);
    }
    // The cells tile the box: total volume and mass are those of the box.
    assert_close(
        sim.particles.iter().map(|p| p.volume).sum(),
        0.1 * 0.05 * 0.075,
        1e-12,
        "volume",
    );
    assert_close(sim.total_mass(), 1600.0 * 0.1 * 0.05 * 0.075, 1e-12, "mass");
    // First particle sits in the centre of the first cell.
    assert_vec_close(
        sim.particles[0].position,
        [0.0125, 0.0125, 0.0125],
        1e-12,
        "first",
    );
    // A spacing that does not divide the extent stretches the lattice to fit.
    let mut stretched = SandSimulation::new(0.05, [0.0; 3]);
    assert_eq!(
        stretched
            .fill_box([0.0; 3], [0.1, 0.1, 0.1], 0.03, 1000.0)
            .ok(),
        Some(27)
    );
    assert_close(
        stretched.total_mass(),
        1000.0 * 1.0e-3,
        1e-12,
        "stretched mass",
    );
    // Filling appends.
    let before = sim.particles.len();
    assert!(
        sim.fill_box([1.0, 0.0, 0.0], [1.05, 0.05, 0.05], 0.025, 1600.0)
            .is_ok()
    );
    assert_eq!(sim.particles.len(), before + 8);
}

#[test]
fn fill_box_rejects_bad_input_and_adds_nothing() {
    let mut sim = SandSimulation::new(0.05, [0.0; 3]);
    let bad: Vec<([f64; 3], [f64; 3], f64, f64)> = vec![
        ([0.0; 3], [0.0, 1.0, 1.0], 0.1, 1000.0),
        ([0.0; 3], [-1.0, 1.0, 1.0], 0.1, 1000.0),
        ([f64::NAN, 0.0, 0.0], [1.0; 3], 0.1, 1000.0),
        ([0.0; 3], [f64::INFINITY, 1.0, 1.0], 0.1, 1000.0),
        ([0.0; 3], [1.0; 3], 0.0, 1000.0),
        ([0.0; 3], [1.0; 3], -0.1, 1000.0),
        ([0.0; 3], [1.0; 3], f64::NAN, 1000.0),
        ([0.0; 3], [1.0; 3], 0.1, 0.0),
        ([0.0; 3], [1.0; 3], 0.1, f64::NAN),
        // 10^4 cells per axis = 10^12 particles.
        ([0.0; 3], [1.0e3; 3], 0.1, 1000.0),
    ];
    for (min, max, spacing, density) in bad {
        let result = sim.fill_box(min, max, spacing, density);
        assert!(
            matches!(result, Err(Error::InvalidParameter(_))),
            "{min:?} {max:?} {spacing} {density}"
        );
        assert!(sim.particles.is_empty());
    }
    // A box thinner than half a spacing holds no particles (not an error).
    assert_eq!(
        sim.fill_box([0.0; 3], [1.0, 1.0, 0.004], 0.1, 1000.0).ok(),
        Some(0)
    );
    assert!(sim.particles.is_empty());
}

// ---------------------------------------------------------------------------
// Time stepping
// ---------------------------------------------------------------------------

#[test]
fn stable_time_step_matches_the_material_wave_speed() {
    let sim = block(0.05, 0.025);
    let mat = Material::new(
        sim.youngs_modulus_pa,
        sim.poisson_ratio,
        sim.friction_angle_deg,
    )
    .expect("valid");
    let expected = sim.cfl * sim.grid_size / (mat.p_wave_modulus() / 1600.0).sqrt();
    let got = sim.stable_time_step().expect("valid parameters");
    assert_close(got, expected, 1e-12, "stable dt");
    // A stiffer sand needs a smaller step, a larger CFL a larger one.
    let mut stiff = sim.clone();
    stiff.youngs_modulus_pa *= 100.0;
    assert!(stiff.stable_time_step().is_ok_and(|dt| dt < 0.11 * got));
    let mut generous = sim.clone();
    generous.cfl = 1.0;
    assert_close(
        generous.stable_time_step().expect("valid"),
        got / 0.6,
        1e-12,
        "cfl scaling",
    );
}

#[test]
fn large_steps_are_sub_stepped_and_agree_with_small_steps() {
    let mut once = block(0.05, 0.025);
    let mut chunked = once.clone();
    assert!(
        once.stable_time_step().is_ok_and(|dt| dt < 0.002),
        "the step below must need many sub-steps"
    );
    once.try_step(0.2).expect("sub-stepped step");
    for _ in 0..8 {
        chunked.try_step(0.025).expect("small step");
    }
    let mean = |s: &SandSimulation| {
        s.particles.iter().map(|p| p.position[1]).sum::<f64>() / s.particles.len() as f64
    };
    assert_close(mean(&once), mean(&chunked), 2e-3, "mean height");
    assert!(
        once.particles
            .iter()
            .chain(&chunked.particles)
            .all(|p| p.position.iter().all(|x| x.is_finite()))
    );
    // A one-second step in a single call (dt ≫ dt_stable) must neither explode nor
    // produce non-finite values.
    let mut wild = block(0.05, 0.025);
    wild.try_step(1.0).expect("one second in one call");
    assert!(wild.kinetic_energy().is_finite() && wild.kinetic_energy() < 1.0);
}

#[test]
fn a_step_that_would_need_too_many_sub_steps_fails_after_partial_progress() {
    let mut sim = SandSimulation::new(0.05, [0.0, -9.81, 0.0]);
    sim.cfl = 1.0;
    sim.add_particle(SandParticle::new([0.0, 5.0, 0.0], 1.6e-3, 1.0e-6));
    let before = sim.particles.clone();
    // c ≈ 17 m/s → dt_stable ≈ 3e-3 s: 20 000 sub-steps cover ≈ 60 s.
    let result = sim.try_step(1.0e4);
    assert!(
        matches!(&result, Err(Error::InvalidParameter(msg)) if msg.contains(&MAX_SUBSTEPS_PER_STEP.to_string())),
        "{result:?}"
    );
    assert_ne!(
        sim.particles[0].position, before[0].position,
        "partial progress stays applied"
    );
    assert!(sim.particles[0].position[1] >= 0.0);
}

#[test]
fn step_is_the_infallible_form_of_try_step() {
    let mut a = block(0.05, 0.025);
    let mut b = a.clone();
    a.step(0.03);
    b.try_step(0.03).expect("valid");
    assert_same_states(&a.particles, &b.particles, "step vs try_step");
}

// ---------------------------------------------------------------------------
// Inert particles and state ownership
// ---------------------------------------------------------------------------

#[test]
fn inert_particles_are_untouched_and_do_not_disturb_the_others() {
    let mut clean = block(0.05, 0.025);
    let mut dirty = clean.clone();
    let inert: Vec<SandParticle> = {
        let mut list = Vec::new();
        let mut nan_pos = SandParticle::new([0.0, 0.3, 0.0], 1.0e-3, 1.0e-6);
        nan_pos.position[0] = f64::NAN;
        list.push(nan_pos);
        list.push(SandParticle::new([0.0, 0.3, 0.0], 0.0, 1.0e-6));
        list.push(SandParticle::new([0.0, 0.3, 0.0], 1.0e-3, -1.0e-6));
        list.push(SandParticle::new([1.0e9, 0.3, 0.0], 1.0e-3, 1.0e-6));
        let mut nan_f = SandParticle::new([0.01, 0.2, 0.0], 1.0e-3, 1.0e-6);
        nan_f.deformation_gradient[0] = f64::NAN;
        list.push(nan_f);
        let mut inf_v = SandParticle::new([0.0, 0.1, 0.01], 1.0e-3, 1.0e-6);
        inf_v.velocity[1] = f64::INFINITY;
        list.push(inf_v);
        list
    };
    // Interleave the inert particles with the active ones.
    let active = dirty.particles.clone();
    dirty.particles.clear();
    let mut inert_iter = inert.iter().cloned();
    for (i, p) in active.into_iter().enumerate() {
        if i % 20 == 3
            && let Some(extra) = inert_iter.next()
        {
            dirty.particles.push(extra);
        }
        dirty.particles.push(p);
    }
    dirty.particles.extend(inert_iter);
    let inert_before: Vec<Vec<u64>> = dirty
        .particles
        .iter()
        .filter(|p| !super::mpm::is_active(p, 1.0 / dirty.grid_size))
        .map(bits)
        .collect();
    assert_eq!(inert_before.len(), 6);

    for _ in 0..6 {
        clean.try_step(0.02).expect("clean");
        dirty.try_step(0.02).expect("dirty");
    }
    let inert_after: Vec<Vec<u64>> = dirty
        .particles
        .iter()
        .filter(|p| !super::mpm::is_active(p, 1.0 / dirty.grid_size))
        .map(bits)
        .collect();
    assert_eq!(inert_before, inert_after, "inert particles must not change");
    let active_dirty: Vec<SandParticle> = dirty
        .particles
        .iter()
        .filter(|p| super::mpm::is_active(p, 1.0 / dirty.grid_size))
        .cloned()
        .collect();
    assert_same_states(&clean.particles, &active_dirty, "active particles");
}

#[test]
fn all_state_lives_in_the_particles() {
    // No hidden per-particle data: continuing from a copy of the particles in a
    // brand-new simulation reproduces the original run bit for bit.
    let mut original = block(0.05, 0.025);
    for _ in 0..5 {
        original.try_step(0.02).expect("step");
    }
    assert!(
        original
            .particles
            .iter()
            .any(|p| p.affine_velocity.iter().any(|c| *c != 0.0)),
        "APIC state must have developed"
    );
    let mut copy = SandSimulation::new(original.grid_size, original.gravity);
    copy.cfl = original.cfl;
    copy.particles = original.particles.clone();
    for _ in 0..5 {
        original.try_step(0.02).expect("original");
        copy.try_step(0.02).expect("copy");
    }
    assert_same_states(&original.particles, &copy.particles, "continued run");
    // Removing a particle between steps is allowed.
    copy.particles.pop();
    assert!(copy.try_step(0.02).is_ok());
}

#[test]
fn identical_runs_are_bitwise_reproducible() {
    let mut a = block(0.05, 0.025);
    let mut b = a.clone();
    for _ in 0..8 {
        a.try_step(0.015).expect("a");
        b.try_step(0.015).expect("b");
    }
    assert_same_states(&a.particles, &b.particles, "reproducibility");
}

// ---------------------------------------------------------------------------
// Boundaries
// ---------------------------------------------------------------------------

#[test]
fn particles_below_the_floor_are_lifted_onto_it() {
    let mut sim = SandSimulation::new(0.05, [0.0; 3]);
    let mut p = SandParticle::new([0.0, -0.3, 0.0], 1.0e-3, 1.0e-6);
    p.velocity = [0.0, -4.0, 0.0];
    sim.add_particle(p);
    sim.try_step(0.001).expect("step");
    assert!(sim.particles[0].position[1] >= 0.0);
    assert!(sim.particles[0].velocity[1] >= 0.0);
}

#[test]
fn plane_strain_slab_keeps_the_motion_two_dimensional() {
    let h = 0.025;
    let mut sim = slab_column(30.0, 0.1, 0.12, 0.5 * h);
    assert_eq!(sim.plane_strain_z, Some([0.0, h]));
    // One particle starts outside the slab and is pulled in.
    sim.particles
        .push(SandParticle::new([0.0, 0.05, 0.3], 1.0e-3, 1.0e-6));
    sim.particles
        .push(SandParticle::new([0.0, 0.05, -0.2], 1.0e-3, 1.0e-6));
    for _ in 0..10 {
        sim.try_step(0.03).expect("step");
    }
    for p in &sim.particles {
        assert!(
            p.position[2] >= 0.0 && p.position[2] <= h,
            "z = {}",
            p.position[2]
        );
    }
    // The slab particles stay exactly on their z planes with zero z velocity.
    let (slab, strays) = sim.particles.split_at(sim.particles.len() - 2);
    assert!(slab.iter().all(|p| p.velocity[2] == 0.0), "v_z must vanish");
    assert!(
        slab.iter()
            .all(|p| p.position[2] > 0.0 && p.position[2] < h)
    );
    assert!(strays.iter().all(|p| (0.0..=h).contains(&p.position[2])));
    // The collapse itself happened in the x-y plane.
    let half_width = slab.iter().map(|p| p.position[0].abs()).fold(0.0, f64::max);
    assert!(
        half_width > 0.07,
        "the column should have spread: {half_width}"
    );
}

// ---------------------------------------------------------------------------
// Conservation and dissipation
// ---------------------------------------------------------------------------

#[test]
fn mass_is_conserved_and_a_dropped_block_comes_to_rest_on_the_floor() {
    let mut sim = SandSimulation::new(0.05, [0.0, -9.81, 0.0]);
    sim.cfl = 0.6;
    sim.fill_box([-0.1, 0.1, -0.1], [0.1, 0.25, 0.1], 0.025, 1600.0)
        .expect("fill");
    let mass0 = sim.total_mass();
    let mut peak = 0.0f64;
    for _ in 0..40 {
        sim.try_step(0.025).expect("step");
        peak = peak.max(sim.kinetic_energy());
        assert!(
            sim.particles.iter().all(|p| p.position[1] >= 0.0),
            "tunnelled through the floor"
        );
    }
    assert_eq!(sim.total_mass(), mass0, "mass changed");
    assert!(
        peak > 1.0,
        "the block must have fallen and gathered speed: {peak}"
    );
    assert!(
        sim.kinetic_energy() < 0.01 * peak,
        "no rest: {} vs peak {peak}",
        sim.kinetic_energy()
    );
    // It rests as a heap on the floor: the lowest particles are near y = 0.
    let lowest = sim
        .particles
        .iter()
        .map(|p| p.position[1])
        .fold(f64::INFINITY, f64::min);
    assert!(lowest < 0.02, "lowest particle {lowest}");
    assert!(
        sim.particles
            .iter()
            .all(|p| p.position.iter().all(|x| x.is_finite()))
    );
    assert!(sim.particles.iter().all(|p| p.plastic_strain_vol >= 0.0));
}

#[test]
fn free_flight_conserves_total_momentum_exactly() {
    // No gravity, far from the floor: internal forces only. The grid forces
    // of a particle sum to zero (Σ_i ∇w_ip = 0), so Σ m v is a constant of the
    // motion even while the block deforms.
    let mut sim = SandSimulation::new(0.05, [0.0; 3]);
    let mut rng = Lcg::new(77);
    sim.cfl = 0.6;
    sim.fill_box([-0.1, 1.0, -0.1], [0.1, 1.1, 0.1], 0.025, 1600.0)
        .expect("fill");
    for p in &mut sim.particles {
        p.velocity = [
            rng.range(-1.0, 1.0),
            rng.range(-1.0, 1.0),
            rng.range(-1.0, 1.0),
        ];
    }
    let p0 = sim.total_momentum();
    let mut peak = 0.0f64;
    for _ in 0..20 {
        sim.try_step(0.01).expect("step");
        peak = peak.max(sim.kinetic_energy());
        assert_vec_close(sim.total_momentum(), p0, 1e-9, "momentum");
    }
    // The random velocity field really deformed the block (plastic work).
    assert!(sim.particles.iter().any(|p| {
        let f = mat_from_col_major(&p.deformation_gradient);
        max_abs_diff(&f, &super::mpm::IDENTITY) > 1e-3
    }));
}

#[test]
fn a_three_dimensional_drop_stays_symmetric_in_x_and_z() {
    let mut sim = SandSimulation::new(0.05, [0.0, -9.81, 0.0]);
    sim.cfl = 0.6;
    sim.fill_box([-0.1, 0.0, -0.1], [0.1, 0.12, 0.1], 0.025, 1600.0)
        .expect("fill");
    let n = sim.particles.len() as f64;
    for _ in 0..6 {
        sim.try_step(0.05).expect("step");
    }
    let centre: Vec<f64> = (0..3)
        .map(|d| sim.particles.iter().map(|p| p.position[d]).sum::<f64>() / n)
        .collect();
    assert!(
        centre[0].abs() < 1e-9 && centre[2].abs() < 1e-9,
        "centre of mass drifted: {centre:?}"
    );
    let half_x = sim
        .particles
        .iter()
        .map(|p| p.position[0].abs())
        .fold(0.0, f64::max);
    let half_z = sim
        .particles
        .iter()
        .map(|p| p.position[2].abs())
        .fold(0.0, f64::max);
    assert_close(half_x, half_z, 1e-9, "x/z symmetry of the footprint");
    assert!(half_x > 0.105, "the slumping block must widen: {half_x}");
    let mean_height = sim.particles.iter().map(|p| p.position[1]).sum::<f64>() / n;
    assert!(
        mean_height > 0.04,
        "the block must still stand on the floor: {mean_height}"
    );
}

#[test]
fn totals_helpers_sum_over_the_particles() {
    let mut sim = SandSimulation::new(0.05, [0.0; 3]);
    let mut a = SandParticle::new([0.0, 1.0, 0.0], 2.0, 1.0e-3);
    let mut b = SandParticle::new([0.2, 1.0, 0.0], 3.0, 1.0e-3);
    a.velocity = [1.0, 0.0, 0.0];
    b.velocity = [0.0, -2.0, 0.0];
    sim.add_particle(a);
    sim.add_particle(b);
    assert_eq!(sim.total_mass(), 5.0);
    assert_eq!(sim.total_momentum(), [2.0, -6.0, 0.0]);
    assert_close(
        sim.kinetic_energy(),
        0.5 * 2.0 + 0.5 * 3.0 * 4.0,
        1e-15,
        "kinetic energy",
    );
}
