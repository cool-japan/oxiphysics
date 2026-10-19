// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0
//! Regression tests for the `fracture_rigid` generators.
//!
//! * The pseudo-random helper used to be mapped onto `[0, 0.5)`, which confined
//!   every seed and every explosion fragment to a single octant.
//! * `ImpactFracture::fracture_body` must conserve mass and linear momentum and
//!   deposit exactly the prescribed relative kinetic energy.
//! * The Sadovsky blast correlations are in MPa and Pa·s respectively.

use oxiphysics_rigid::fracture_rigid::{ExplosionModel, FragmentSystem, ImpactFracture};

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn norm(a: [f64; 3]) -> f64 {
    (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt()
}

/// `a` and `b` agree to the relative tolerance `rel` (absolute below 1).
fn close(a: f64, b: f64, rel: f64) -> bool {
    (a - b).abs() <= rel * a.abs().max(b.abs()).max(1.0)
}

// ── pseudo-random seeds and directions ──────────────────────────────────────

#[test]
fn seeds_cover_both_signs_of_every_axis() {
    let model = ImpactFracture::new(1.0, 200, 0.5, 1.0);
    let centre = [3.0, -2.0, 0.5];
    let seeds = model.generate_seeds(centre);
    assert_eq!(seeds.len(), 200);
    assert_eq!(seeds[0].position, centre, "first seed sits on the impact");
    let offsets: Vec<[f64; 3]> = seeds.iter().map(|s| sub(s.position, centre)).collect();
    for axis in 0..3 {
        assert!(
            offsets.iter().any(|o| o[axis] > 0.1),
            "no seed on the positive side of axis {axis}"
        );
        assert!(
            offsets.iter().any(|o| o[axis] < -0.1),
            "no seed on the negative side of axis {axis}"
        );
        let mean = offsets.iter().map(|o| o[axis]).sum::<f64>() / offsets.len() as f64;
        assert!(
            mean.abs() < 0.15,
            "seeds are lopsided along axis {axis}: mean offset {mean}"
        );
    }
    for o in &offsets {
        assert!(norm(*o) <= 1.0 + 1.0e-9, "seed outside the scatter radius");
    }
}

#[test]
fn seeds_are_deterministic_for_one_impact_point() {
    let model = ImpactFracture::new(1.0, 12, 0.5, 0.4);
    let a = model.generate_seeds([0.25, -1.5, 2.0]);
    let b = model.generate_seeds([0.25, -1.5, 2.0]);
    assert_eq!(a.len(), b.len());
    for (sa, sb) in a.iter().zip(b.iter()) {
        assert_eq!(sa.position, sb.position);
    }
}

#[test]
fn distinct_impact_points_give_distinct_scatter_patterns() {
    let model = ImpactFracture::new(1.0, 8, 0.5, 1.0);
    // Negative coordinates used to saturate to zero in the seed hash, so these
    // two points produced the very same pattern of offsets.
    let p = [-0.5, 0.0, 0.0];
    let q = [-0.7, 0.0, 0.0];
    let a = model.generate_seeds(p);
    let b = model.generate_seeds(q);
    let differing = a
        .iter()
        .zip(b.iter())
        .skip(1)
        .filter(|(sa, sb)| sub(sa.position, p) != sub(sb.position, q))
        .count();
    assert_eq!(differing, 7, "every scattered seed should differ");

    // Mirror images are different points as well.
    let c = model.generate_seeds([0.5, 0.0, 0.0]);
    let mirrored = a
        .iter()
        .zip(c.iter())
        .skip(1)
        .any(|(sa, sc)| sub(sa.position, p) != sub(sc.position, [0.5, 0.0, 0.0]));
    assert!(mirrored);
}

#[test]
fn zero_seeds_requested_gives_no_seeds_and_no_fragments() {
    let model = ImpactFracture::new(1.0, 0, 0.5, 1.0);
    assert!(model.generate_seeds([0.0; 3]).is_empty());
    let system = model.fracture_body(4.0, [0.0; 3], [1.0; 3], [10.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
    assert!(system.fragments.is_empty());
}

fn explosion_system(count: usize) -> FragmentSystem {
    ExplosionModel::new(1.0, 2440.0, 1.0, 5.0, count).create_fragment_system([0.0, 5.0, 0.0])
}

#[test]
fn explosion_fragments_fly_in_every_direction() {
    let system = explosion_system(2000);
    assert_eq!(system.fragments.len(), 2000);
    for axis in 0..3 {
        let positive = system
            .fragments
            .iter()
            .filter(|f| f.velocity[axis] > 0.0)
            .count();
        let negative = system
            .fragments
            .iter()
            .filter(|f| f.velocity[axis] < 0.0)
            .count();
        // Directions are isotropic, so each half space holds about half of them.
        assert!(
            (800..=1200).contains(&positive) && (800..=1200).contains(&negative),
            "axis {axis}: {positive} fragments forward, {negative} backward"
        );
    }
}

#[test]
fn explosion_fragments_carry_almost_no_net_momentum() {
    let system = explosion_system(2000);
    let mut momentum = [0.0_f64; 3];
    let mut scale = 0.0;
    for f in &system.fragments {
        for (m, v) in momentum.iter_mut().zip(f.velocity.iter()) {
            *m += f.mass * v;
        }
        scale += f.mass * norm(f.velocity);
    }
    let ratio = norm(momentum) / scale;
    // Independent random directions leave a residual of order 1 / sqrt(N); a
    // generator that confines the directions to one octant gives about 0.7.
    assert!(ratio < 0.1, "net momentum / sum(m|v|) = {ratio}");
}

#[test]
fn explosion_mass_distribution_matches_the_casing() {
    let model = ExplosionModel::new(1.0, 2440.0, 2.0, 5.0, 2000);
    let distribution = model.fragment_distribution();
    let total = distribution.total_mass();
    assert!(
        (total - 2.0).abs() < 0.1,
        "fragment masses add up to {total} kg for a 2 kg casing"
    );
}

// ── ImpactFracture conservation ─────────────────────────────────────────────

const MASS: f64 = 12.5;
const IMPULSE: [f64; 3] = [30.0, 120.0, -40.0];
const IMPACT: [f64; 3] = [0.2, 1.0, -0.3];

fn conservation_report(model: &ImpactFracture) -> (usize, f64, [f64; 3], f64) {
    let system = model.fracture_body(MASS, [0.0; 3], [1.0; 3], IMPULSE, IMPACT);
    let bulk = [IMPULSE[0] / MASS, IMPULSE[1] / MASS, IMPULSE[2] / MASS];
    let mut mass = 0.0;
    let mut momentum = [0.0_f64; 3];
    let mut relative_energy = 0.0;
    for f in &system.fragments {
        mass += f.mass;
        for (m, v) in momentum.iter_mut().zip(f.velocity.iter()) {
            *m += f.mass * v;
        }
        let rel = sub(f.velocity, bulk);
        relative_energy += 0.5 * f.mass * norm(rel) * norm(rel);
    }
    (system.fragments.len(), mass, momentum, relative_energy)
}

fn assert_conserved(model: &ImpactFracture, expected_count: usize) {
    let (count, mass, momentum, relative_energy) = conservation_report(model);
    assert_eq!(count, expected_count);
    assert!(
        close(mass, MASS, 1.0e-12),
        "fragment masses add up to {mass}"
    );
    for k in 0..3 {
        assert!(
            (momentum[k] - IMPULSE[k]).abs() < 1.0e-9,
            "momentum {momentum:?} differs from the impulse {IMPULSE:?}"
        );
    }
    let expected = model.energy_transfer * norm(IMPULSE) * norm(IMPULSE) / (2.0 * MASS);
    assert!(
        close(relative_energy, expected, 1.0e-9),
        "relative kinetic energy {relative_energy} differs from {expected}"
    );
}

#[test]
fn four_fragments_conserve_mass_momentum_and_energy() {
    let model = ImpactFracture::new(10.0, 4, 0.5, 0.3);
    assert_conserved(&model, 4);
}

#[test]
fn forty_fragments_conserve_mass_momentum_and_energy() {
    let mut model = ImpactFracture::new(10.0, 40, 0.5, 0.3);
    // No lower bound on the fragment size, so all 40 seeds become fragments.
    model.min_fragment_mass_fraction = 0.0;
    assert_conserved(&model, 40);
}

#[test]
fn minimum_fragment_mass_fraction_caps_the_count_without_inflating_the_mass() {
    // The default floor of 5 % allows at most 20 fragments.  It used to raise
    // every fragment of a 40-seed fracture to 5 %, giving twice the body mass.
    let model = ImpactFracture::new(10.0, 40, 0.5, 0.3);
    assert_conserved(&model, 20);
    let system = model.fracture_body(MASS, [0.0; 3], [1.0; 3], IMPULSE, IMPACT);
    for f in &system.fragments {
        assert!(close(f.mass, MASS / 20.0, 1.0e-12), "mass {}", f.mass);
        assert!(f.mass >= 0.05 * MASS - 1.0e-12);
    }

    let mut loose = ImpactFracture::new(10.0, 40, 0.5, 0.3);
    loose.min_fragment_mass_fraction = 0.1;
    assert_conserved(&loose, 10);
}

#[test]
fn energy_share_follows_the_energy_transfer_fraction() {
    for transfer in [0.0, 0.1, 0.5, 1.0, 2.0] {
        let model = ImpactFracture::new(10.0, 6, transfer, 0.3);
        assert_conserved(&model, 6);
    }
}

#[test]
fn a_single_fragment_moves_with_the_parent() {
    let model = ImpactFracture::new(10.0, 1, 0.5, 0.3);
    let system = model.fracture_body(MASS, [0.0; 3], [1.0; 3], IMPULSE, IMPACT);
    assert_eq!(system.fragments.len(), 1);
    let f = &system.fragments[0];
    assert!(close(f.mass, MASS, 1.0e-12));
    for (v, j) in f.velocity.iter().zip(IMPULSE.iter()) {
        assert!((v - j / MASS).abs() < 1.0e-12);
    }
}

#[test]
fn fragments_start_without_spin_and_inside_the_scatter_region() {
    let model = ImpactFracture::new(10.0, 8, 0.5, 0.3);
    let system = model.fracture_body(MASS, [0.0; 3], [1.0; 3], IMPULSE, IMPACT);
    for f in &system.fragments {
        assert_eq!(f.angular_velocity, [0.0; 3]);
        assert!(norm(sub(f.position, IMPACT)) <= 0.3 + 1.0e-9);
    }
}

#[test]
fn gravity_comes_from_the_model() {
    let mut model = ImpactFracture::new(10.0, 4, 0.5, 0.3);
    assert_eq!(model.gravity, [0.0, -9.81, 0.0]);
    model.gravity = [0.0, 0.0, -1.62];
    let system = model.fracture_body(MASS, [0.0; 3], [1.0; 3], IMPULSE, IMPACT);
    assert_eq!(system.gravity, [0.0, 0.0, -1.62]);
}

#[test]
fn non_positive_mass_or_weak_impulse_gives_no_fragments() {
    let model = ImpactFracture::new(10.0, 4, 0.5, 0.3);
    for mass in [0.0, -3.0, f64::NAN, f64::INFINITY] {
        let system = model.fracture_body(mass, [0.0; 3], [1.0; 3], IMPULSE, IMPACT);
        assert!(system.fragments.is_empty(), "mass {mass}");
    }
    let weak = model.fracture_body(MASS, [0.0; 3], [1.0; 3], [1.0, 0.0, 0.0], IMPACT);
    assert!(weak.fragments.is_empty());
}

// ── Sadovsky blast correlations ─────────────────────────────────────────────

fn charge(mass_kg: f64) -> ExplosionModel {
    ExplosionModel::new(mass_kg, 2440.0, 1.0, 5.0, 10)
}

#[test]
fn overpressure_at_unit_scaled_distance_is_about_one_megapascal() {
    // z = 1: 0.084 + 0.27 + 0.7 = 1.054 MPa.
    let p = charge(1.0).peak_overpressure(1.0);
    assert!(close(p, 1.054e6, 1.0e-9), "p = {p}");
}

#[test]
fn one_kilogram_at_five_metres_gives_about_33_kilopascal() {
    // z = 5: 0.084/5 + 0.27/25 + 0.7/125 = 0.0332 MPa.
    let p = charge(1.0).peak_overpressure(5.0);
    assert!(close(p, 33_200.0, 1.0e-9), "p = {p}");
    assert!((20.0e3..50.0e3).contains(&p));
}

#[test]
fn overpressure_depends_only_on_the_scaled_distance() {
    let reference = charge(1.0).peak_overpressure(5.0);
    // 8 kg at 10 m and 27 kg at 15 m both have z = 5.
    assert!(close(
        charge(8.0).peak_overpressure(10.0),
        reference,
        1.0e-9
    ));
    assert!(close(
        charge(27.0).peak_overpressure(15.0),
        reference,
        1.0e-9
    ));
}

#[test]
fn overpressure_decays_monotonically_with_distance() {
    let model = charge(2.0);
    let mut previous = f64::INFINITY;
    for step in 1..=60 {
        let p = model.peak_overpressure(0.5 * f64::from(step));
        assert!(p < previous, "overpressure must fall with range");
        previous = p;
    }
}

#[test]
fn blast_impulse_follows_the_sadovsky_correlation() {
    // 200 * W^(2/3) / r: one kilogram at five metres is 40 Pa s.
    let i = charge(1.0).blast_impulse(5.0);
    assert!(close(i, 40.0, 1.0e-9), "i = {i}");
}

#[test]
fn blast_impulse_scales_with_the_two_thirds_power_of_the_charge() {
    let one = charge(1.0).blast_impulse(4.0);
    let eight = charge(8.0).blast_impulse(4.0);
    let twenty_seven = charge(27.0).blast_impulse(4.0);
    assert!(close(eight / one, 4.0, 1.0e-9), "ratio {}", eight / one);
    assert!(
        close(twenty_seven / one, 9.0, 1.0e-9),
        "ratio {}",
        twenty_seven / one
    );
    // And it falls off as 1 / r at a fixed charge.
    let far = charge(8.0).blast_impulse(8.0);
    assert!(close(eight / far, 2.0, 1.0e-9));
}

#[test]
fn blast_quantities_vanish_without_range_or_charge() {
    assert_eq!(charge(1.0).peak_overpressure(0.0), 0.0);
    assert_eq!(charge(1.0).blast_impulse(-1.0), 0.0);
    assert_eq!(charge(0.0).peak_overpressure(5.0), 0.0);
    assert_eq!(charge(0.0).blast_impulse(5.0), 0.0);
}
