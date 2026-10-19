// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Behavioural tests of the sand simulation: things the previous ballistic
//! free-fall integrator could not do (stand on the floor, hold a slope, flow
//! when the slope is too steep, spread and settle, slide with friction).
//!
//! All runs are plane-strain slabs one cell thick (a 2-D problem at a fraction
//! of the 3-D cost) on a grid of 25 mm cells with two particles per cell edge.
//! Note that the Drucker-Prager cone is calibrated to the Mohr-Coulomb friction
//! angle in *triaxial compression*; in plane strain it admits steeper heaps than
//! `φ` (a heap of `φ = 30°` sand stands up to about 45°), which is why the
//! stability tests compare against measured outcomes at slopes well away from
//! the boundary rather than against `tan φ`.

use super::test_support::{assert_close, half_width, max_height, mean_height, slab_column, wedge};

const SPACING: f64 = 0.0125;
const HALF_BASE: f64 = 0.1;
const WEDGE_TIME: f64 = 0.6;

/// Height retained by a pre-built heap after `WEDGE_TIME`: final / initial
/// maximum height (1 = unmoved, → 0 = flowed out flat).
fn retained_height(friction_angle_deg: f64, slope_deg: f64) -> f64 {
    let mut sim = wedge(friction_angle_deg, slope_deg, HALF_BASE, SPACING);
    let initial = max_height(&sim);
    assert!(initial > 0.0 && sim.particles.len() > 30);
    let mut elapsed = 0.0;
    while elapsed < WEDGE_TIME - 1e-12 {
        assert!(sim.try_step(0.05).is_ok());
        elapsed += 0.05;
    }
    assert!(
        sim.particles
            .iter()
            .all(|p| p.position.iter().all(|x| x.is_finite()))
    );
    max_height(&sim) / initial
}

// ---------------------------------------------------------------------------
// Standing on the floor
// ---------------------------------------------------------------------------

#[test]
fn a_flat_block_rests_on_the_floor() {
    // Ballistic free fall would drop every particle to y = 0 (mean height 0).
    let mut sim = slab_column(30.0, 0.2, 0.04, SPACING);
    let initial_mean = mean_height(&sim);
    let initial_top = max_height(&sim);
    for _ in 0..10 {
        assert!(sim.try_step(0.05).is_ok());
    }
    assert!(
        mean_height(&sim) > 0.98 * initial_mean,
        "block sank: {}",
        mean_height(&sim) / initial_mean
    );
    assert!(max_height(&sim) > 0.97 * initial_top);
    assert!(sim.kinetic_energy() < 1.0e-5, "KE {}", sim.kinetic_energy());
    assert!(
        half_width(&sim) < 0.11,
        "a wide flat block barely spreads: {}",
        half_width(&sim)
    );
    assert!(sim.particles.iter().all(|p| p.position[1] >= 0.0));
}

// ---------------------------------------------------------------------------
// Slope stability depends on the friction angle
// ---------------------------------------------------------------------------

#[test]
fn heaps_well_below_the_limit_stand() {
    // φ = 30° at 20°, φ = 40° at 30°: 10° below the friction angle.
    for (phi, slope) in [(30.0, 20.0), (40.0, 30.0)] {
        let retained = retained_height(phi, slope);
        assert!(
            retained > 0.98,
            "phi={phi} slope={slope}: retained {retained}"
        );
    }
}

#[test]
fn a_heap_slumps_when_the_friction_angle_is_too_small_for_its_slope() {
    // Same 40° heap: it flows for φ = 15° and stands for φ = 30°.
    let weak = retained_height(15.0, 40.0);
    let strong = retained_height(30.0, 40.0);
    assert!(weak < 0.8, "φ=15° must slump: retained {weak}");
    assert!(strong > 0.97, "φ=30° must stand: retained {strong}");
}

#[test]
fn retained_height_grows_with_the_friction_angle() {
    let retained: Vec<f64> = [15.0, 20.0, 25.0, 30.0]
        .iter()
        .map(|&phi| retained_height(phi, 40.0))
        .collect();
    for pair in retained.windows(2) {
        assert!(pair[1] > pair[0] + 0.03, "not increasing: {retained:?}");
    }
}

#[test]
fn a_steeper_heap_than_any_friction_angle_supports_collapses() {
    // 60° is beyond the plane-strain limit of φ = 30° ...
    let collapsed = retained_height(30.0, 60.0);
    assert!(
        collapsed < 0.85,
        "60° heap of φ=30° sand must flow: retained {collapsed}"
    );
    // ... but φ = 40° sand still carries it.
    let stands = retained_height(40.0, 60.0);
    assert!(stands > 0.95, "φ=40° sand must hold 60°: retained {stands}");
}

// ---------------------------------------------------------------------------
// Collapse of a column
// ---------------------------------------------------------------------------

#[test]
fn a_column_collapses_spreads_and_settles() {
    let mut sim = slab_column(30.0, 0.12, 0.18, SPACING);
    let mass = sim.total_mass();
    let initial_half_width = half_width(&sim);
    let initial_top = max_height(&sim);
    let mut peak = 0.0f64;
    for _ in 0..12 {
        assert!(sim.try_step(0.05).is_ok());
        peak = peak.max(sim.kinetic_energy());
    }
    assert_eq!(sim.total_mass(), mass);
    assert!(
        max_height(&sim) < 0.6 * initial_top,
        "top {} of {initial_top}",
        max_height(&sim)
    );
    assert!(
        half_width(&sim) > 2.5 * initial_half_width,
        "half-width {}",
        half_width(&sim)
    );
    assert!(peak > 0.01, "the collapse must release energy: {peak}");
    assert!(
        sim.kinetic_energy() < 1.0e-3 * peak,
        "not at rest: {} vs {peak}",
        sim.kinetic_energy()
    );
    assert!(sim.particles.iter().all(|p| p.position[1] >= 0.0));
    // Spreading loosened the sand: dilation was recorded and never decreases.
    assert!(sim.particles.iter().any(|p| p.plastic_strain_vol > 0.0));
    assert!(sim.particles.iter().all(|p| p.plastic_strain_vol >= 0.0));
}

#[test]
fn runout_decreases_and_the_heap_grows_taller_with_the_friction_angle() {
    let mut outcomes = Vec::new();
    for phi in [20.0, 30.0, 40.0] {
        let mut sim = slab_column(phi, 0.12, 0.18, SPACING);
        for _ in 0..12 {
            assert!(sim.try_step(0.05).is_ok());
        }
        outcomes.push((half_width(&sim), max_height(&sim)));
    }
    for pair in outcomes.windows(2) {
        assert!(
            pair[1].0 < pair[0].0 - 0.015,
            "runout must shrink: {outcomes:?}"
        );
        assert!(
            pair[1].1 > pair[0].1 + 0.005,
            "heap must grow taller: {outcomes:?}"
        );
    }
}

// ---------------------------------------------------------------------------
// Floor friction
// ---------------------------------------------------------------------------

#[test]
fn a_sliding_block_decelerates_at_mu_times_g() {
    // A block of high-friction sand sliding on a floor with Coulomb friction µ:
    // the floor force is µ M g, so the mean velocity falls as v0 − µ g t.
    let (mu, v0) = (0.4, 1.0);
    let mut sim = slab_column(45.0, 0.2, 0.05, SPACING);
    sim.floor_friction = mu;
    for p in &mut sim.particles {
        p.velocity = [v0, 0.0, 0.0];
    }
    let mass = sim.total_mass();
    let mut t = 0.0;
    for _ in 0..4 {
        assert!(sim.try_step(0.05).is_ok());
        t += 0.05;
        let mean_velocity = sim.total_momentum()[0] / mass;
        assert_close(
            mean_velocity,
            v0 - mu * 9.81 * t,
            0.012,
            &format!("mean velocity at t={t}"),
        );
    }
    // It comes to rest and stays there.
    for _ in 0..4 {
        assert!(sim.try_step(0.05).is_ok());
    }
    assert!((sim.total_momentum()[0] / mass).abs() < 0.01);
}

#[test]
fn a_frictionless_floor_does_not_slow_a_sliding_block_and_a_sticky_one_stops_it() {
    let run = |floor_friction: f64| {
        let mut sim = slab_column(45.0, 0.2, 0.05, SPACING);
        sim.floor_friction = floor_friction;
        for p in &mut sim.particles {
            p.velocity = [1.0, 0.0, 0.0];
        }
        let mass = sim.total_mass();
        for _ in 0..4 {
            assert!(sim.try_step(0.05).is_ok());
        }
        sim.total_momentum()[0] / mass
    };
    assert!(
        run(0.0) > 0.99,
        "a frictionless floor must not brake the block"
    );
    assert!(
        run(f64::INFINITY).abs() < 0.02,
        "a sticky floor must stop it"
    );
}
