// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Regression tests for the wind coupling of [`HairSimulation`].
//!
//! The aerodynamic drag is a force (N).  It used to be added directly to the
//! particle positions (m), without `dt^2 / m`, and, because the previous
//! positions were not touched, it also injected a spurious velocity into the
//! Verlet scheme.  Drag now enters as an acceleration `F / m` that
//! [`HairStrand::integrate_verlet`] integrates over `dt^2`.

use super::*;

const AIR_DENSITY: f64 = 1.225;
const DRAG_COEFFICIENT: f64 = 1.0;
const RADIUS: f64 = 5.0e-5;
const PARTICLE_MASS: f64 = 1.0e-6;
const STRAND_LENGTH: f64 = 0.3;
const SEGMENTS: usize = 10;
const SEGMENT_LENGTH: f64 = STRAND_LENGTH / SEGMENTS as f64;
const ROOT: [f64; 3] = [0.0, 1.0, 0.0];
const PI_F64: f64 = std::f64::consts::PI;

fn simulation() -> HairSimulation {
    let mut sim = HairSimulation::new(AIR_DENSITY, DRAG_COEFFICIENT, RADIUS)
        .with_particle_mass(PARTICLE_MASS);
    sim.add_strand(ROOT, STRAND_LENGTH, SEGMENTS);
    sim
}

/// Drag acceleration (m/s²) of a particle that is at rest in a wind of `speed`.
fn drag_acceleration(speed: f64, mass: f64) -> f64 {
    0.5 * AIR_DENSITY * DRAG_COEFFICIENT * (2.0 * RADIUS * SEGMENT_LENGTH) * speed * speed / mass
}

#[test]
fn wind_drag_is_returned_as_an_acceleration_along_the_wind() {
    let sim = simulation();
    let positions_before = sim.strands[0].positions.clone();
    let acc = sim.apply_wind_forces([4.0, 0.0, 0.0]);

    assert_eq!(acc.len(), 1, "one entry per strand");
    assert_eq!(acc[0].len(), SEGMENTS + 1, "one entry per particle");
    assert_eq!(acc[0][0], [0.0; 3], "the pinned root feels no wind");
    let expected = drag_acceleration(4.0, PARTICLE_MASS);
    assert!(
        expected > 1.0 && expected < 1.0e3,
        "m/s^2 scale: {expected}"
    );
    for (i, a) in acc[0].iter().enumerate().skip(1) {
        assert!(
            (a[0] - expected).abs() < 1.0e-9 * expected,
            "particle {i}: a_x = {} vs {expected}",
            a[0]
        );
        assert_eq!(a[1], 0.0);
        assert_eq!(a[2], 0.0);
    }
    // Querying the load does not move anything.
    assert_eq!(sim.strands[0].positions, positions_before);
}

#[test]
fn drag_scales_with_speed_squared_and_inverse_particle_mass() {
    let sim = simulation();
    let slow = sim.apply_wind_forces([2.0, 0.0, 0.0])[0][5][0];
    let fast = sim.apply_wind_forces([4.0, 0.0, 0.0])[0][5][0];
    assert!(
        (fast / slow - 4.0).abs() < 1.0e-12,
        "v^2 law: {slow} -> {fast}"
    );

    let heavy = simulation().with_particle_mass(2.0 * PARTICLE_MASS);
    let heavy_acc = heavy.apply_wind_forces([4.0, 0.0, 0.0])[0][5][0];
    assert!(
        (fast / heavy_acc - 2.0).abs() < 1.0e-12,
        "a = F / m: {fast} vs {heavy_acc}"
    );
}

#[test]
fn drag_uses_the_velocity_relative_to_the_air() {
    let mut sim = simulation();
    let wind = [4.0, 0.0, 0.0];
    // Particles moving with the wind feel no drag ...
    for v in sim.strands[0].velocities.iter_mut() {
        *v = wind;
    }
    for a in sim.apply_wind_forces(wind)[0].iter() {
        assert_eq!(*a, [0.0; 3]);
    }
    // ... a particle against the wind sees the larger relative speed ...
    sim.strands[0].velocities[3] = [-2.0, 0.0, 0.0];
    let against = sim.apply_wind_forces(wind)[0][3][0];
    let expected = drag_acceleration(6.0, PARTICLE_MASS);
    assert!((against - expected).abs() < 1.0e-9 * expected);
    // ... and one that outruns the wind is decelerated.
    sim.strands[0].velocities[3] = [8.0, 0.0, 0.0];
    let outrunning = sim.apply_wind_forces(wind)[0][3][0];
    let expected = -drag_acceleration(4.0, PARTICLE_MASS);
    assert!((outrunning - expected).abs() < 1.0e-9 * expected.abs());
}

#[test]
fn invalid_particle_mass_switches_the_wind_off() {
    for mass in [0.0, -1.0e-6, f64::NAN, f64::INFINITY] {
        let sim = simulation().with_particle_mass(mass);
        for strand in sim.apply_wind_forces([5.0, 1.0, -2.0]) {
            for a in strand {
                assert_eq!(a, [0.0; 3], "mass {mass}");
            }
        }
    }
}

#[test]
fn default_particle_mass_is_a_keratin_segment() {
    let sim = HairSimulation::new(AIR_DENSITY, DRAG_COEFFICIENT, RADIUS);
    let expected = 1300.0 * PI_F64 * RADIUS * RADIUS * 0.01;
    assert!((sim.particle_mass - expected).abs() < 1.0e-12 * expected);
    assert!(sim.particle_mass > 1.0e-8 && sim.particle_mass < 1.0e-6);
}

#[test]
fn verlet_integrates_extra_acceleration_over_dt_squared() {
    let dt = 2.0e-3;
    let gravity = [0.0, -9.81, 0.0];
    let mut extra = vec![[0.0_f64; 3]; SEGMENTS + 1];
    extra[3] = [5.0, -1.0, 2.0];
    extra[7] = [-3.0, 0.5, 0.25];

    let mut strand = HairStrand::new(ROOT, STRAND_LENGTH, SEGMENTS);
    let before = strand.positions.clone();
    strand.integrate_verlet(dt, gravity, &extra);
    assert_eq!(strand.positions[0], before[0], "the root stays pinned");
    for i in 1..=SEGMENTS {
        for d in 0..3 {
            // From rest: x_new = x + (g + a_extra) dt^2.
            let expected = before[i][d] + (gravity[d] + extra[i][d]) * dt * dt;
            assert!(
                (strand.positions[i][d] - expected).abs() < 1.0e-15,
                "particle {i} component {d}"
            );
        }
    }

    // A short (or empty) slice means "no extra acceleration" for the rest.
    let mut short = HairStrand::new(ROOT, STRAND_LENGTH, SEGMENTS);
    short.integrate_verlet(dt, gravity, &extra[..4]);
    let mut gravity_only = HairStrand::new(ROOT, STRAND_LENGTH, SEGMENTS);
    gravity_only.integrate_verlet(dt, gravity, &[]);
    assert_eq!(short.positions[3], strand.positions[3]);
    for i in 4..=SEGMENTS {
        assert_eq!(
            short.positions[i], gravity_only.positions[i],
            "particle {i}"
        );
    }
}

#[test]
fn one_step_wind_displacement_scales_with_dt_squared() {
    let wind = [4.0, 0.0, 0.0];
    let tip_displacement = |dt: f64| {
        let mut sim = simulation();
        sim.step(dt, [0.0; 3], wind);
        sim.strands[0].tip_position()[0]
    };
    let dt = 1.0e-3;
    let single = tip_displacement(dt);
    let doubled = tip_displacement(2.0 * dt);

    // x = a dt^2 for a particle starting at rest ...
    let expected = drag_acceleration(4.0, PARTICLE_MASS) * dt * dt;
    assert!(
        (single - expected).abs() < 0.01 * expected,
        "one-step displacement {single} vs a dt^2 = {expected}"
    );
    // ... so doubling dt quadruples it.
    assert!(
        (doubled / single - 4.0).abs() < 0.04,
        "dt^2 law: {single} -> {doubled}"
    );
}

#[test]
fn wind_pushes_a_hanging_strand_downwind() {
    let mut sim = simulation();
    let gravity = [0.0, -9.81, 0.0];
    let mut max_tip_x = 0.0_f64;
    for _ in 0..200 {
        sim.step(1.0e-3, gravity, [1.5, 0.0, 0.0]);
        let tip = sim.strands[0].tip_position();
        assert!(tip.iter().all(|c| c.is_finite()), "tip must stay finite");
        max_tip_x = max_tip_x.max(tip[0]);
    }
    assert!(
        max_tip_x > 1.0e-3,
        "the strand must swing downwind: {max_tip_x}"
    );
    let strand = &sim.strands[0];
    assert!(
        strand.tip_position()[1] < ROOT[1],
        "the tip stays below the root"
    );
    assert!(
        (strand.arc_length() - STRAND_LENGTH).abs() < 0.02 * STRAND_LENGTH,
        "length constraints keep the strand intact: {}",
        strand.arc_length()
    );
    // No wind: the strand just hangs.
    let mut calm = simulation();
    for _ in 0..200 {
        calm.step(1.0e-3, gravity, [0.0; 3]);
    }
    assert!(calm.strands[0].tip_position()[0].abs() < 1.0e-12);
}
