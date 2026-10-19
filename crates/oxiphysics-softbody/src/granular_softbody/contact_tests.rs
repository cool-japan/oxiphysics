// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Regression tests for the Hertz-Mindlin contact kinematics and for the
//! cell-list contact detection of [`DemFemSystem`].
//!
//! * The contact point lies `+r_a n` from grain A and `-r_b n` from grain B for
//!   a contact normal `n` pointing from A to B.  Friction must therefore brake
//!   a spinning grain; with the lever arms flipped it accelerated the spin.
//! * Contact detection has to match the brute-force all-pairs loop exactly
//!   while touching only candidate pairs, and the Mindlin tangential history
//!   may only exist for pairs that are in contact.

use super::*;

const TOL: f64 = 1.0e-9;

fn sand_grain(id: usize, position: Vec3, radius: f64) -> Grain {
    Grain::new(id, position, radius, 2650.0, 70.0e9, 0.3, 50.0, 10.0)
}

/// Deterministic uniform numbers in `[0, 1)` (64-bit LCG).
fn uniform(seed: &mut u64) -> f64 {
    *seed = seed
        .wrapping_mul(6364136223846793005)
        .wrapping_add(1442695040888963407);
    ((*seed >> 11) as f64) / ((1u64 << 53) as f64)
}

fn centred(seed: &mut u64) -> f64 {
    2.0 * uniform(seed) - 1.0
}

#[test]
fn spinning_grain_is_braked_by_friction_against_a_resting_grain() {
    let params = DemParams::default_sand();
    let r = 1.0e-3;
    let omega = 50.0;
    let mut ga = sand_grain(0, Vec3::zero(), r);
    ga.omega = Vec3::new(0.0, 0.0, omega);
    let gb = sand_grain(1, Vec3::new(1.9e-3, 0.0, 0.0), r);
    let overlap = 2.0 * r - 1.9e-3;

    let mut delta_t = Vec3::zero();
    let (fa, ta, fb, tb) = hertz_mindlin_force(&ga, &gb, &mut delta_t, &params);

    // The contact point of A moves in +y (omega x r_a n = z x x), so the
    // friction force on A points to -y, on B to +y.
    assert!(fa.y < 0.0 && fb.y > 0.0, "friction forces: {fa:?} / {fb:?}");
    assert!(
        ta.dot(ga.omega) < 0.0,
        "friction torque {ta:?} must oppose the spin {:?}",
        ga.omega
    );
    // B is dragged along its own contact point: gear-like counter-rotation.
    assert!(tb.z < 0.0, "torque on the dragged grain: {tb:?}");

    // Magnitude: slip speed omega r, stick regime (spring + dashpot).
    let slip = omega * r;
    let g_single = 70.0e9 / (2.0 * (1.0 + 0.3));
    let g_eff = 1.0 / (1.0 / g_single + 1.0 / g_single);
    let r_eff = 0.5 * r;
    let k_t = 8.0 * g_eff * (r_eff * overlap).sqrt();
    let expected_ft = k_t * slip * params.dt + ga.gamma_t * slip;
    assert!(
        (fb.y - expected_ft).abs() < TOL * expected_ft,
        "tangential force {} vs {expected_ft}",
        fb.y
    );
    assert!(
        (ta.z + r * expected_ft).abs() < TOL * r * expected_ft,
        "torque on A {} vs {}",
        ta.z,
        -r * expected_ft
    );
}

#[test]
fn spinning_grain_is_braked_when_the_resting_grain_is_on_the_other_side() {
    let params = DemParams::default_sand();
    let r = 1.0e-3;
    // A resting grain below, a spinning one on top: contact normal -y.
    let ga = sand_grain(0, Vec3::new(0.0, 1.9e-3, 0.0), r);
    let mut gb = sand_grain(1, Vec3::zero(), r);
    gb.omega = Vec3::new(0.0, 0.0, 30.0);
    let mut delta_t = Vec3::zero();
    let (_, _, _, tb) = hertz_mindlin_force(&ga, &gb, &mut delta_t, &params);
    assert!(
        tb.dot(gb.omega) < 0.0,
        "torque {tb:?} on the spinning second grain must oppose its spin {:?}",
        gb.omega
    );
}

#[test]
fn rolling_contact_without_slip_carries_no_friction() {
    let params = DemParams::default_sand();
    let r = 1.0e-3;
    let v = 0.2;
    // Grain A rolls over the resting grain B (normal +x): the contact point of
    // A is at rest when v + omega x (r n) = 0, i.e. omega_z = -v / r.
    let mut ga = sand_grain(0, Vec3::zero(), r);
    ga.velocity = Vec3::new(0.0, v, 0.0);
    ga.omega = Vec3::new(0.0, 0.0, -v / r);
    let gb = sand_grain(1, Vec3::new(1.9e-3, 0.0, 0.0), r);

    let mut delta_t = Vec3::zero();
    let (fa, ta, fb, tb) = hertz_mindlin_force(&ga, &gb, &mut delta_t, &params);
    let normal_force = fb.x;
    assert!(normal_force > 0.0);
    for (what, force) in [("A", fa), ("B", fb)] {
        assert!(
            force.y.abs() < 1.0e-9 * normal_force && force.z.abs() < 1.0e-9 * normal_force,
            "grain {what} feels a friction force although the contact points do not slide: {force:?}"
        );
    }
    assert!(
        ta.norm() < 1.0e-12 * normal_force * r,
        "torque on A: {ta:?}"
    );
    assert!(
        tb.norm() < 1.0e-12 * normal_force * r,
        "torque on B: {tb:?}"
    );
    assert!(
        delta_t.norm() < 1.0e-18,
        "no tangential displacement builds up"
    );
}

#[test]
fn tangential_contact_forces_never_add_energy() {
    let params = DemParams::default_sand();
    let mut seed = 0x5eed_1234_abcd_ef01_u64;
    let mut checked = 0;
    for _ in 0..300 {
        let direction = Vec3::new(centred(&mut seed), centred(&mut seed), centred(&mut seed));
        if direction.norm() < 0.1 {
            continue;
        }
        let n = direction.normalise();
        let ra = 1.0e-3 * (0.5 + uniform(&mut seed));
        let rb = 1.0e-3 * (0.5 + uniform(&mut seed));
        let mut ga = sand_grain(0, Vec3::zero(), ra);
        let mut gb = sand_grain(1, n.scale(0.9 * (ra + rb)), rb);
        ga.velocity = Vec3::new(centred(&mut seed), centred(&mut seed), centred(&mut seed));
        gb.velocity = Vec3::new(centred(&mut seed), centred(&mut seed), centred(&mut seed));
        ga.omega = Vec3::new(
            200.0 * centred(&mut seed),
            200.0 * centred(&mut seed),
            200.0 * centred(&mut seed),
        );
        gb.omega = Vec3::new(
            200.0 * centred(&mut seed),
            200.0 * centred(&mut seed),
            200.0 * centred(&mut seed),
        );

        let mut delta_t = Vec3::zero();
        let (fa, ta, fb, tb) = hertz_mindlin_force(&ga, &gb, &mut delta_t, &params);
        let power_total =
            fa.dot(ga.velocity) + ta.dot(ga.omega) + fb.dot(gb.velocity) + tb.dot(gb.omega);
        // The normal force acts along the line of centres; the tangential force
        // is perpendicular to it, so fb . n is the normal force magnitude.
        let normal_force = fb.dot(n);
        let power_normal = normal_force * n.dot(gb.velocity - ga.velocity);
        let power_tangential = power_total - power_normal;
        let scale = power_total.abs() + power_normal.abs() + 1.0e-30;
        assert!(
            power_tangential <= TOL * scale,
            "friction injects power {power_tangential} (total {power_total}, normal {power_normal})"
        );
        checked += 1;
    }
    assert!(
        checked > 250,
        "too few random contacts were checked: {checked}"
    );
}

/// Random polydisperse packing in the unit cube with random motion.
fn random_system(n: usize, r_min: f64, r_span: f64, seed: u64) -> DemFemSystem {
    let mut seed = seed;
    let mut sys = DemFemSystem::new(DemParams::default_sand(), None);
    for _ in 0..n {
        let radius = r_min + r_span * uniform(&mut seed);
        let position = Vec3::new(uniform(&mut seed), uniform(&mut seed), uniform(&mut seed));
        let mut grain = sand_grain(0, position, radius);
        grain.velocity = Vec3::new(centred(&mut seed), centred(&mut seed), centred(&mut seed));
        grain.omega = Vec3::new(
            50.0 * centred(&mut seed),
            50.0 * centred(&mut seed),
            50.0 * centred(&mut seed),
        );
        sys.add_grain(grain);
    }
    sys
}

/// Number of grain pairs whose spheres overlap (independent of the contact code).
fn overlapping_pairs(grains: &[Grain]) -> usize {
    let mut count = 0;
    for i in 0..grains.len() {
        for j in (i + 1)..grains.len() {
            let dist = (grains[j].position - grains[i].position).norm();
            if dist < grains[i].radius + grains[j].radius {
                count += 1;
            }
        }
    }
    count
}

type Reference = (Vec<Vec3>, Vec<Vec3>, Vec<(usize, usize, Vec3, f64)>);

/// The original all-pairs algorithm with a fresh tangential history.
fn brute_force_reference(sys: &DemFemSystem) -> Reference {
    let n = sys.grains.len();
    let mut forces = vec![Vec3::zero(); n];
    let mut torques = vec![Vec3::zero(); n];
    let mut contacts = Vec::new();
    for i in 0..n {
        for j in (i + 1)..n {
            let mut delta_t = Vec3::zero();
            let (fa, ta, fb, tb) =
                hertz_mindlin_force(&sys.grains[i], &sys.grains[j], &mut delta_t, &sys.params);
            let fn_mag = fa.norm();
            if fn_mag > 0.0 {
                let r = sys.grains[j].position - sys.grains[i].position;
                contacts.push((i, j, r.scale(1.0 / r.norm()), fn_mag));
            }
            forces[i] = forces[i] + fa;
            torques[i] = torques[i] + ta;
            forces[j] = forces[j] + fb;
            torques[j] = torques[j] + tb;
        }
    }
    (forces, torques, contacts)
}

#[test]
fn contact_detection_matches_brute_force_on_a_random_packing() {
    let mut sys = random_system(300, 0.04, 0.03, 0xd15c_0de5_0000_0001);
    let touching = overlapping_pairs(&sys.grains);
    assert!(
        touching > 100,
        "the packing must be dense enough: {touching}"
    );

    let (forces, torques, contacts) = brute_force_reference(&sys);
    sys.compute_contacts();

    assert_eq!(sys.contacts, contacts, "contact list");
    for (i, grain) in sys.grains.iter().enumerate() {
        assert_eq!(grain.force, forces[i], "force on grain {i}");
        assert_eq!(grain.torque, torques[i], "torque on grain {i}");
    }
    assert_eq!(
        sys.tangential_cache.len(),
        touching,
        "the tangential history must hold exactly the touching pairs, not all n(n-1)/2 pairs"
    );
    for &(i, j) in sys.tangential_cache.keys() {
        let dist = (sys.grains[j].position - sys.grains[i].position).norm();
        assert!(i < j && dist < sys.grains[i].radius + sys.grains[j].radius);
    }
}

#[test]
fn contact_detection_is_repeatable_and_keeps_history_of_persistent_contacts() {
    let mut sys = random_system(120, 0.05, 0.03, 0x0123_4567_89ab_cdef);
    sys.compute_contacts();
    let keys: std::collections::BTreeSet<(usize, usize)> =
        sys.tangential_cache.keys().copied().collect();
    assert!(!keys.is_empty(), "the packing must contain contacts");
    let first_history: Vec<Vec3> = keys.iter().map(|k| sys.tangential_cache[k]).collect();
    for grain in &mut sys.grains {
        grain.clear_forces();
    }
    sys.compute_contacts();
    let again: std::collections::BTreeSet<(usize, usize)> =
        sys.tangential_cache.keys().copied().collect();
    assert_eq!(keys, again, "the same pairs are in contact");
    // The spring displacement of a persistent contact keeps accumulating.
    let grown = keys
        .iter()
        .zip(first_history.iter())
        .filter(|(k, first)| sys.tangential_cache[*k].norm() > first.norm())
        .count();
    assert!(
        grown > keys.len() / 2,
        "history grew for {grown} of {} pairs",
        keys.len()
    );
}

#[test]
fn large_sparse_system_keeps_a_cache_of_the_size_of_its_contacts() {
    let mut sys = random_system(2000, 0.008, 0.004, 0xfeed_beef_cafe_0042);
    let touching = overlapping_pairs(&sys.grains);
    assert!(
        touching > 0 && touching < 2000,
        "sparse packing: {touching}"
    );
    sys.compute_contacts();
    assert_eq!(sys.tangential_cache.len(), touching);
}

#[test]
fn tangential_history_lives_only_as_long_as_the_contact() {
    let mut sys = DemFemSystem::new(DemParams::default_sand(), None);
    let mut ga = sand_grain(0, Vec3::zero(), 1.0e-3);
    ga.velocity = Vec3::new(0.0, 0.1, 0.0);
    sys.add_grain(ga);
    let touching_position = Vec3::new(1.9e-3, 0.0, 0.0);
    sys.add_grain(sand_grain(1, touching_position, 1.0e-3));

    sys.compute_contacts();
    assert_eq!(sys.tangential_cache.len(), 1);
    let first = sys.tangential_cache[&(0, 1)];
    assert!(first.norm() > 0.0, "sliding builds up a tangential spring");

    for grain in &mut sys.grains {
        grain.clear_forces();
    }
    sys.compute_contacts();
    let second = sys.tangential_cache[&(0, 1)];
    assert!(
        second.norm() > first.norm(),
        "history accumulates while touching"
    );

    // Separate the grains: the history must be dropped ...
    sys.grains[1].position = Vec3::new(1.0, 0.0, 0.0);
    for grain in &mut sys.grains {
        grain.clear_forces();
    }
    sys.compute_contacts();
    assert!(
        sys.tangential_cache.is_empty(),
        "no history without a contact"
    );
    assert!(sys.contacts.is_empty());

    // ... so that a new contact starts from zero, not from the stale spring.
    sys.grains[1].position = touching_position;
    for grain in &mut sys.grains {
        grain.clear_forces();
    }
    sys.compute_contacts();
    assert_eq!(
        sys.tangential_cache[&(0, 1)],
        first,
        "a fresh contact restarts with the displacement of a single step"
    );
}

#[test]
fn crushed_and_degenerate_grains_make_no_contacts() {
    let mut sys = DemFemSystem::new(DemParams::default_sand(), None);
    sys.compute_contacts();
    assert!(sys.contacts.is_empty() && sys.tangential_cache.is_empty());

    // Overlapping pair with one crushed grain: skipped.
    sys.add_grain(sand_grain(0, Vec3::zero(), 1.0e-3));
    let mut crushed = sand_grain(1, Vec3::new(1.0e-3, 0.0, 0.0), 1.0e-3);
    crushed.crushed = true;
    sys.add_grain(crushed);
    sys.compute_contacts();
    assert!(sys.contacts.is_empty() && sys.tangential_cache.is_empty());
    assert_eq!(sys.grains[0].force, Vec3::zero());

    // Zero-radius grains neither hang the cell list nor create contacts.
    let mut points = DemFemSystem::new(DemParams::default_sand(), None);
    points.add_grain(sand_grain(0, Vec3::zero(), 0.0));
    points.add_grain(sand_grain(1, Vec3::zero(), 0.0));
    points.compute_contacts();
    assert!(points.contacts.is_empty() && points.tangential_cache.is_empty());

    // Non-finite positions are never reported as contacts.
    let mut broken = DemFemSystem::new(DemParams::default_sand(), None);
    broken.add_grain(sand_grain(0, Vec3::zero(), 1.0e-3));
    broken.add_grain(sand_grain(1, Vec3::new(f64::NAN, 0.0, 0.0), 1.0e-3));
    broken.compute_contacts();
    assert!(broken.contacts.is_empty() && broken.tangential_cache.is_empty());
}

#[test]
fn diverged_coordinates_are_skipped_instead_of_overflowing_the_cell_hash() {
    let mut sys = DemFemSystem::new(DemParams::default_sand(), None);
    sys.add_grain(sand_grain(0, Vec3::zero(), 1.0e-3));
    sys.add_grain(sand_grain(1, Vec3::new(1.9e-3, 0.0, 0.0), 1.0e-3));
    // Finite but astronomically distant grains of a blown-up simulation: their
    // cell indices exceed the 32-bit range of the hash, so a neighbour-cell
    // lookup around them would overflow.
    sys.add_grain(sand_grain(2, Vec3::new(1.0e300, 0.0, 0.0), 1.0e-3));
    sys.add_grain(sand_grain(3, Vec3::new(-1.0e30, 1.0e30, 0.0), 1.0e-3));
    sys.add_grain(sand_grain(4, Vec3::new(1.0e300, 0.0, 0.0), 1.0e-3));

    sys.compute_contacts();

    assert_eq!(sys.contacts.len(), 1, "only the healthy pair touches");
    assert_eq!((sys.contacts[0].0, sys.contacts[0].1), (0, 1));
    assert_eq!(sys.tangential_cache.len(), 1);
    for far in 2..5 {
        assert_eq!(sys.grains[far].force, Vec3::zero(), "grain {far}");
    }
    assert!(sys.grains[0].force.norm() > 0.0 && sys.grains[1].force.norm() > 0.0);
}
