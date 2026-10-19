// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Regression tests for the N-body integrator and the Barnes-Hut octree.

use super::{G_CONST, GravBody, NBodyGravity, OctreeNode};
use crate::rng_util::uniform;
use rand::SeedableRng;
use rand::rngs::StdRng;
use std::f64::consts::PI;

fn norm(v: [f64; 3]) -> f64 {
    (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
}

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

fn scale(a: [f64; 3], factor: f64) -> [f64; 3] {
    [a[0] * factor, a[1] * factor, a[2] * factor]
}

/// Bodies scattered uniformly in a cube with random velocities and masses.
fn random_cloud(rng: &mut StdRng, count: usize, half_extent: f64) -> Vec<GravBody> {
    let mut bodies = Vec::with_capacity(count);
    for _ in 0..count {
        let pos: [f64; 3] = std::array::from_fn(|_| uniform(rng, -half_extent, half_extent));
        let vel: [f64; 3] = std::array::from_fn(|_| uniform(rng, -1.0e3, 1.0e3));
        let mass = uniform(rng, 1.0e21, 1.0e23);
        bodies.push(GravBody::new(pos, vel, mass));
    }
    bodies
}

/// Two equal masses on a circular orbit of radius `d` around their centre.
///
/// The circular speed accounts for the Plummer softening `epsilon`, so the
/// pair separation stays exactly `2 d` for the true dynamics.  Returns the
/// bodies and the orbital period.
fn circular_pair(mass: f64, d: f64, epsilon: f64) -> (Vec<GravBody>, f64) {
    let separation = 2.0 * d;
    let pull =
        G_CONST * mass * separation / (separation * separation + epsilon * epsilon).powf(1.5);
    let speed = (d * pull).sqrt();
    let period = 2.0 * PI * d / speed;
    let bodies = vec![
        GravBody::new([-d, 0.0, 0.0], [0.0, -speed, 0.0], mass),
        GravBody::new([d, 0.0, 0.0], [0.0, speed, 0.0], mass),
    ];
    (bodies, period)
}

fn total_energy(sim: &NBodyGravity) -> f64 {
    sim.kinetic_energy() + sim.potential_energy()
}

fn pair_separation(sim: &NBodyGravity) -> f64 {
    norm(sub(sim.bodies[1].pos, sim.bodies[0].pos))
}

fn total_momentum(sim: &NBodyGravity) -> [f64; 3] {
    sim.bodies
        .iter()
        .fold([0.0_f64; 3], |p, body| add(p, scale(body.vel, body.mass)))
}

fn assert_states_identical(a: &NBodyGravity, b: &NBodyGravity, context: &str) {
    assert_eq!(a.bodies.len(), b.bodies.len(), "{context}: body count");
    for (i, (x, y)) in a.bodies.iter().zip(&b.bodies).enumerate() {
        assert_eq!(x.pos, y.pos, "{context}: position of body {i}");
        assert_eq!(x.vel, y.vel, "{context}: velocity of body {i}");
    }
}

// ── forces ───────────────────────────────────────────────────────────────

#[test]
fn barnes_hut_with_zero_opening_angle_reproduces_direct_summation() {
    let mut rng = StdRng::seed_from_u64(11);
    let sim = NBodyGravity::new(random_cloud(&mut rng, 300, 1.0e9));
    let tree = sim.build_tree();
    let direct = sim.accelerations_direct();
    let rms = (direct.iter().map(|a| norm(*a).powi(2)).sum::<f64>() / direct.len() as f64).sqrt();
    for (i, expected) in direct.iter().enumerate() {
        let got = sim.accel_bh(&tree, i, 0.0);
        let err = norm(sub(got, *expected));
        assert!(err <= 1e-12 * rms, "body {i}: error {err} vs rms {rms}");
    }
}

#[test]
fn barnes_hut_half_opening_angle_matches_direct_summation() {
    let mut rng = StdRng::seed_from_u64(2024);
    let mut sim = NBodyGravity::with_barnes_hut(random_cloud(&mut rng, 500, 1.0e9), 0.5);
    sim.epsilon = 1.0e3;
    let tree = sim.build_tree();
    let direct = sim.accelerations_direct();
    let barnes_hut = sim.accelerations();
    assert_eq!(barnes_hut.len(), direct.len());

    let rms = (direct.iter().map(|a| norm(*a).powi(2)).sum::<f64>() / direct.len() as f64).sqrt();
    let mut per_body = Vec::with_capacity(direct.len());
    let mut worst_vs_rms = 0.0_f64;
    let mut squared_error = 0.0_f64;
    let mut squared_reference = 0.0_f64;
    for (i, expected) in direct.iter().enumerate() {
        let err = norm(sub(barnes_hut[i], *expected));
        assert_eq!(barnes_hut[i], sim.accel_bh(&tree, i, 0.5));
        per_body.push(err / norm(*expected));
        worst_vs_rms = worst_vs_rms.max(err / rms);
        squared_error += err * err;
        squared_reference += norm(*expected).powi(2);
    }
    per_body.sort_by(f64::total_cmp);
    let median = per_body[per_body.len() / 2];
    let aggregate = (squared_error / squared_reference).sqrt();
    // The per-body relative error is ill-conditioned for bodies whose net pull
    // nearly cancels (deep inside the cloud), so the worst case is judged
    // against the rms acceleration; the aggregate error covers all bodies at once.
    assert!(
        worst_vs_rms < 1e-2,
        "worst error relative to the rms acceleration: {worst_vs_rms}"
    );
    assert!(median < 3e-3, "median per-body relative error: {median}");
    assert!(
        aggregate < 3e-3,
        "aggregate relative acceleration error: {aggregate}"
    );
}

#[test]
fn cell_multipoles_improve_the_far_field_of_a_compact_cluster() {
    let mut rng = StdRng::seed_from_u64(314);
    let mut bodies = Vec::new();
    for _ in 0..12 {
        let pos: [f64; 3] = std::array::from_fn(|_| uniform(&mut rng, -1.0e7, 1.0e7));
        bodies.push(GravBody::new(
            pos,
            [0.0; 3],
            uniform(&mut rng, 1.0e22, 5.0e22),
        ));
    }
    let probe = bodies.len();
    bodies.push(GravBody::new([1.0e9, 3.0e8, -2.0e8], [0.0; 3], 1.0e20));
    let sim = NBodyGravity::new(bodies);
    let tree = sim.build_tree();

    let exact = sim.accelerations_direct()[probe];
    let got = sim.accel_bh(&tree, probe, 0.5);
    // Monopole-only reference: all cluster mass at its centre of mass.
    let cluster = &sim.bodies[..probe];
    let mass: f64 = cluster.iter().map(|b| b.mass).sum();
    let com = cluster.iter().fold([0.0_f64; 3], |c, body| {
        add(c, scale(body.pos, body.mass / mass))
    });
    let r = sub(sim.bodies[probe].pos, com);
    let rho = (r[0] * r[0] + r[1] * r[1] + r[2] * r[2] + sim.epsilon * sim.epsilon).sqrt();
    let monopole = [
        -G_CONST * mass * r[0] / rho.powi(3),
        -G_CONST * mass * r[1] / rho.powi(3),
        -G_CONST * mass * r[2] / rho.powi(3),
    ];
    let monopole_error = norm(sub(monopole, exact)) / norm(exact);
    let tree_error = norm(sub(got, exact)) / norm(exact);
    assert!(
        monopole_error > 1e-5,
        "cluster must be resolved as one cell: {monopole_error}"
    );
    assert!(
        tree_error < 0.05 * monopole_error,
        "quadrupole should remove the leading error: tree {tree_error}, monopole {monopole_error}"
    );
}

#[test]
fn smaller_opening_angle_is_more_accurate() {
    let mut rng = StdRng::seed_from_u64(5);
    let sim = NBodyGravity::new(random_cloud(&mut rng, 200, 1.0e9));
    let tree = sim.build_tree();
    let direct = sim.accelerations_direct();
    let rms = (direct.iter().map(|a| norm(*a).powi(2)).sum::<f64>() / direct.len() as f64).sqrt();
    let worst = |theta: f64| -> f64 {
        direct
            .iter()
            .enumerate()
            .map(|(i, a)| norm(sub(sim.accel_bh(&tree, i, theta), *a)) / rms)
            .fold(0.0, f64::max)
    };
    let coarse = worst(1.0);
    let medium = worst(0.5);
    let fine = worst(0.2);
    assert!(fine < medium && medium < coarse, "{fine} {medium} {coarse}");
    assert!(fine < 1e-3, "theta 0.2 error {fine}");
}

#[test]
fn two_body_barnes_hut_is_exact_for_any_opening_angle() {
    let bodies = vec![
        GravBody::new([0.0, 0.0, 0.0], [0.0; 3], 6.0e24),
        GravBody::new([3.0e8, 1.0e8, -2.0e8], [0.0; 3], 7.0e22),
    ];
    let exact = NBodyGravity::new(bodies.clone()).accelerations_direct();
    for theta in [0.1, 0.5, 1.0, 50.0] {
        let got = NBodyGravity::with_barnes_hut(bodies.clone(), theta).accelerations();
        for (g, e) in got.iter().zip(&exact) {
            assert!(norm(sub(*g, *e)) <= 1e-14 * norm(*e), "theta {theta}");
        }
    }
}

#[test]
fn direct_forces_obey_newtons_third_law() {
    let mut rng = StdRng::seed_from_u64(77);
    let sim = NBodyGravity::new(random_cloud(&mut rng, 40, 1.0e9));
    let acc = sim.accelerations_direct();
    let mut net = [0.0_f64; 3];
    let mut magnitude = 0.0_f64;
    for (a, body) in acc.iter().zip(&sim.bodies) {
        net = add(net, scale(*a, body.mass));
        magnitude += body.mass * norm(*a);
    }
    assert!(
        norm(net) <= 1e-12 * magnitude,
        "net force {net:?} vs scale {magnitude}"
    );
}

#[test]
fn force_is_the_negative_gradient_of_the_softened_potential() {
    // Softening comparable to the separations makes the old `max(r, eps)`
    // potential disagree with the Plummer force.
    let mut sim = NBodyGravity::new(vec![
        GravBody::new([0.0, 0.0, 0.0], [0.0; 3], 3.0e23),
        GravBody::new([2.0e8, 0.0, 0.0], [0.0; 3], 2.0e23),
        GravBody::new([0.5e8, 1.5e8, 0.4e8], [0.0; 3], 1.0e23),
    ]);
    sim.epsilon = 1.0e8;
    let acc = sim.accelerations_direct();
    let h = 1.0e3;
    for (i, body_acc) in acc.iter().enumerate() {
        for (axis, &accel) in body_acc.iter().enumerate() {
            let mut plus = sim.clone();
            plus.bodies[i].pos[axis] += h;
            let mut minus = sim.clone();
            minus.bodies[i].pos[axis] -= h;
            let gradient = (plus.potential_energy() - minus.potential_energy()) / (2.0 * h);
            let force = sim.bodies[i].mass * accel;
            assert!(
                (force + gradient).abs() <= 1e-6 * force.abs().max(gradient.abs()),
                "body {i} axis {axis}: force {force}, -dU/dx {}",
                -gradient
            );
        }
    }
}

#[test]
fn potential_energy_of_coincident_bodies_is_finite_and_softened() {
    let mut sim = NBodyGravity::new(vec![
        GravBody::new([1.0, 2.0, 3.0], [0.0; 3], 1.0e20),
        GravBody::new([1.0, 2.0, 3.0], [0.0; 3], 2.0e20),
    ]);
    sim.epsilon = 5.0e3;
    let expected = -G_CONST * 1.0e20 * 2.0e20 / 5.0e3;
    assert!((sim.potential_energy() - expected).abs() <= 1e-12 * expected.abs());
    for a in sim.accelerations_direct() {
        assert_eq!(a, [0.0; 3], "coincident bodies exert no softened force");
    }
}

// ── octree ───────────────────────────────────────────────────────────────

fn tree_depth(tree: &[OctreeNode], index: usize) -> usize {
    1 + tree[index]
        .children
        .iter()
        .flatten()
        .map(|&child| tree_depth(tree, child))
        .max()
        .unwrap_or(0)
}

/// Indices of all bodies stored in the subtree rooted at `index`.
fn subtree_bodies(tree: &[OctreeNode], index: usize) -> Vec<usize> {
    let node = &tree[index];
    if node.is_leaf() {
        return node.body_idx.into_iter().collect();
    }
    node.children
        .iter()
        .flatten()
        .flat_map(|&child| subtree_bodies(tree, child))
        .collect()
}

#[test]
fn octree_aggregates_mass_and_centre_of_mass_bottom_up() {
    let mut rng = StdRng::seed_from_u64(99);
    let sim = NBodyGravity::new(random_cloud(&mut rng, 400, 5.0e8));
    let tree = sim.build_tree();
    let total_mass: f64 = sim.bodies.iter().map(|b| b.mass).sum();
    assert!((tree[0].mass - total_mass).abs() <= 1e-12 * total_mass);
    let com = sim.center_of_mass();
    assert!(
        norm(sub(tree[0].com, com)) <= 1e-9 * 5.0e8,
        "{:?} vs {com:?}",
        tree[0].com
    );

    for (index, node) in tree.iter().enumerate() {
        let members = subtree_bodies(&tree, index);
        let mass: f64 = members.iter().map(|&b| sim.bodies[b].mass).sum();
        assert!(
            (node.mass - mass).abs() <= 1e-12 * mass,
            "mass of node {index}"
        );
        let mut expected = [0.0_f64; 6];
        for &b in &members {
            let d = sub(sim.bodies[b].pos, node.com);
            let m = sim.bodies[b].mass;
            let terms = [
                d[0] * d[0],
                d[1] * d[1],
                d[2] * d[2],
                d[0] * d[1],
                d[0] * d[2],
                d[1] * d[2],
            ];
            for (total, term) in expected.iter_mut().zip(terms) {
                *total += m * term;
            }
        }
        let moment_scale = mass * (2.0 * node.half_width).powi(2);
        for (component, (&got, want)) in node.quad.iter().zip(expected).enumerate() {
            assert!(
                (got - want).abs() <= 1e-9 * moment_scale,
                "node {index} quad[{component}]: {got} vs {want}"
            );
        }
        if node.is_leaf() {
            assert_eq!(node.quad, [0.0; 6]);
        }
    }

    let mut seen = vec![0_usize; sim.bodies.len()];
    for (index, node) in tree.iter().enumerate() {
        if node.is_leaf() {
            let body = node
                .body_idx
                .expect("every leaf in a built tree holds a body");
            seen[body] += 1;
            assert_eq!(node.mass, sim.bodies[body].mass);
            assert_eq!(node.com, sim.bodies[body].pos);
            assert!(
                node.contains(sim.bodies[body].pos),
                "leaf {index} must hold its body"
            );
        } else {
            assert!(
                node.body_idx.is_none(),
                "internal node {index} holds no body"
            );
            let child_mass: f64 = node.children.iter().flatten().map(|&c| tree[c].mass).sum();
            assert!((node.mass - child_mass).abs() <= 1e-12 * node.mass);
            for &child in node.children.iter().flatten() {
                assert!(child > index, "children follow their parent in the arena");
                assert_eq!(tree[child].half_width, 0.5 * node.half_width);
            }
        }
    }
    assert!(
        seen.iter().all(|&count| count == 1),
        "each body sits in exactly one leaf"
    );
}

#[test]
fn coincident_and_nearly_coincident_bodies_do_not_blow_up_the_tree() {
    let mut bodies = Vec::new();
    for k in 0..6 {
        bodies.push(GravBody::new(
            [1.0e6, 2.0e6, -3.0e6],
            [0.0; 3],
            1.0e20 * f64::from(k + 1),
        ));
    }
    bodies.push(GravBody::new([-4.0e8, 1.0e8, 0.0], [0.0; 3], 5.0e22));
    bodies.push(GravBody::new([2.0e8, 2.0e8, 2.0e8], [0.0; 3], 5.0e22));
    // Two bodies separated by far less than any sensible cell size.
    bodies.push(GravBody::new([3.0e8, -1.0e8, 5.0e7], [0.0; 3], 1.0e22));
    bodies.push(GravBody::new(
        [3.0e8 + 1.0e-6, -1.0e8, 5.0e7],
        [0.0; 3],
        2.0e22,
    ));
    let total_mass: f64 = bodies.iter().map(|b| b.mass).sum();

    let sim = NBodyGravity::with_barnes_hut(bodies, 0.5);
    let tree = sim.build_tree();
    assert!((tree[0].mass - total_mass).abs() <= 1e-12 * total_mass);
    assert!(tree_depth(&tree, 0) <= 66, "depth {}", tree_depth(&tree, 0));
    let direct = sim.accelerations_direct();
    for (got, expected) in sim.accelerations().iter().zip(&direct) {
        assert!(got.iter().all(|c| c.is_finite()));
        assert!(norm(sub(*got, *expected)) <= 0.05 * norm(*expected));
    }
}

#[test]
fn degenerate_bodies_fall_back_gracefully() {
    let empty = NBodyGravity::with_barnes_hut(Vec::new(), 0.5);
    assert!(empty.build_tree().is_empty());
    assert!(empty.accelerations().is_empty());

    let single =
        NBodyGravity::with_barnes_hut(vec![GravBody::new([1.0; 3], [0.0; 3], 1.0e20)], 0.5);
    assert_eq!(single.build_tree().len(), 1);
    assert_eq!(single.accelerations(), vec![[0.0; 3]]);

    let nan = NBodyGravity::with_barnes_hut(
        vec![
            GravBody::new([f64::NAN, 0.0, 0.0], [0.0; 3], 1.0e20),
            GravBody::new([1.0e6, 0.0, 0.0], [0.0; 3], 1.0e20),
        ],
        0.5,
    );
    assert!(
        nan.build_tree().is_empty(),
        "non-finite positions cannot be bounded"
    );
    assert_eq!(
        nan.accelerations().len(),
        2,
        "falls back to direct summation"
    );

    let tree = single.build_tree();
    assert_eq!(single.accel_bh(&tree, 7, 0.5), [0.0; 3], "invalid index");
    assert_eq!(single.accel_bh(&[], 0, 0.5), [0.0; 3], "empty tree");
}

// ── time integration ─────────────────────────────────────────────────────

#[test]
fn leapfrog_conserves_energy_on_a_circular_two_body_orbit() {
    let mass = 1.0e24;
    let d = 1.0e8;
    for theta in [0.0, 0.5] {
        let (bodies, period) = circular_pair(mass, d, 1.0e3);
        let mut sim = NBodyGravity::with_barnes_hut(bodies, theta);
        let e0 = total_energy(&sim);
        let dt = period / 4000.0;
        let mut worst_energy = 0.0_f64;
        let mut worst_radius = 0.0_f64;
        for _ in 0..1000 {
            sim.step(dt);
            worst_energy = worst_energy.max(((total_energy(&sim) - e0) / e0).abs());
            worst_radius = worst_radius.max((pair_separation(&sim) / (2.0 * d) - 1.0).abs());
        }
        assert!(
            worst_energy < 1e-6,
            "theta {theta}: energy drift {worst_energy}"
        );
        assert!(
            worst_radius < 1e-5,
            "theta {theta}: radius error {worst_radius}"
        );
    }
}

#[test]
fn leapfrog_has_no_secular_energy_drift_over_many_orbits() {
    let (bodies, period) = circular_pair(1.0e24, 1.0e8, 1.0e3);
    let mut sim = NBodyGravity::new(bodies);
    let e0 = total_energy(&sim);
    let steps_per_orbit = 200;
    let dt = period / f64::from(steps_per_orbit);
    let mut early = 0.0_f64;
    let mut late = 0.0_f64;
    for orbit in 0..50 {
        let mut worst = 0.0_f64;
        for _ in 0..steps_per_orbit {
            sim.step(dt);
            worst = worst.max(((total_energy(&sim) - e0) / e0).abs());
        }
        if orbit < 10 {
            early = early.max(worst);
        } else if orbit >= 40 {
            late = late.max(worst);
        }
    }
    assert!(early < 1e-3, "energy error should be small: {early}");
    assert!(
        late <= 1.5 * early,
        "energy error grows secularly: early {early}, late {late}"
    );
}

#[test]
fn leapfrog_is_time_reversible() {
    let mut rng = StdRng::seed_from_u64(31);
    let mut sim = NBodyGravity::new(random_cloud(&mut rng, 8, 1.0e9));
    let start = sim.bodies.clone();
    let dt = 2.0e3;
    for _ in 0..100 {
        sim.step(dt);
    }
    for body in &mut sim.bodies {
        body.vel = [-body.vel[0], -body.vel[1], -body.vel[2]];
    }
    for _ in 0..100 {
        sim.step(dt);
    }
    for (now, then) in sim.bodies.iter().zip(&start) {
        assert!(
            norm(sub(now.pos, then.pos)) <= 1e-3,
            "{:?} vs {:?}",
            now.pos,
            then.pos
        );
        let speed_err = norm([
            now.vel[0] + then.vel[0],
            now.vel[1] + then.vel[1],
            now.vel[2] + then.vel[2],
        ]);
        assert!(speed_err <= 1e-9, "velocity not reversed: {speed_err}");
    }
}

#[test]
fn direct_integration_conserves_momentum_com_motion_and_angular_momentum() {
    let mut rng = StdRng::seed_from_u64(1234);
    let mut sim = NBodyGravity::new(random_cloud(&mut rng, 12, 1.0e9));
    let total_mass: f64 = sim.bodies.iter().map(|b| b.mass).sum();
    let p0 = total_momentum(&sim);
    let l0 = sim.angular_momentum();
    let com0 = sim.center_of_mass();
    let momentum_scale: f64 = sim.bodies.iter().map(|b| b.mass * norm(b.vel)).sum();
    let dt = 1.0e4;
    let steps = 200;
    for _ in 0..steps {
        sim.step(dt);
    }
    let p1 = total_momentum(&sim);
    assert!(
        norm(sub(p1, p0)) <= 1e-10 * momentum_scale,
        "momentum drift"
    );

    let elapsed = dt * f64::from(steps);
    let predicted = [
        com0[0] + p0[0] / total_mass * elapsed,
        com0[1] + p0[1] / total_mass * elapsed,
        com0[2] + p0[2] / total_mass * elapsed,
    ];
    assert!(
        norm(sub(sim.center_of_mass(), predicted)) <= 1e-3,
        "centre of mass must move uniformly"
    );

    let l1 = sim.angular_momentum();
    let l_scale: f64 = sim
        .bodies
        .iter()
        .map(|b| b.mass * norm(b.pos) * norm(b.vel))
        .sum();
    assert!(
        norm(sub(l1, l0)) <= 1e-9 * l_scale,
        "angular momentum drift"
    );
}

#[test]
fn acceleration_cache_never_changes_the_numerics() {
    let mut rng = StdRng::seed_from_u64(8);
    let bodies = random_cloud(&mut rng, 20, 1.0e9);
    let dt = 3.0e3;
    for theta in [0.0, 0.6] {
        let mut cached = NBodyGravity::with_barnes_hut(bodies.clone(), theta);
        let mut fresh = NBodyGravity::with_barnes_hut(bodies.clone(), theta);
        for _ in 0..5 {
            cached.step(dt);
            // A brand-new integrator has no cache and must recompute everything.
            fresh = {
                let mut again = NBodyGravity::with_barnes_hut(fresh.bodies.clone(), theta);
                again.step(dt);
                again
            };
            assert_states_identical(&cached, &fresh, "cached vs fresh integrator");
        }
    }
}

/// A modification of the public state of a system between two steps.
type SystemEdit = fn(&mut NBodyGravity);

#[test]
fn stale_acceleration_cache_is_detected_after_external_edits() {
    let mut rng = StdRng::seed_from_u64(64);
    let bodies = random_cloud(&mut rng, 10, 1.0e9);
    let dt = 2.0e3;

    let edits: [(&str, SystemEdit); 6] = [
        ("move a body", |s| s.bodies[4].pos[1] += 3.0e7),
        ("change a mass", |s| s.bodies[2].mass *= 3.0),
        ("remove a body", |s| {
            s.bodies.pop();
        }),
        ("add a body", |s| {
            s.bodies.push(GravBody::new([2.0e8; 3], [0.0; 3], 4.0e22));
        }),
        ("change softening", |s| s.epsilon = 5.0e7),
        ("change opening angle", |s| s.theta = 0.4),
    ];
    for (label, edit) in &edits {
        let mut sim = NBodyGravity::new(bodies.clone());
        sim.step(dt);
        edit(&mut sim);
        let mut reference = NBodyGravity::new(sim.bodies.clone());
        reference.epsilon = sim.epsilon;
        reference.theta = sim.theta;
        sim.step(dt);
        reference.step(dt);
        assert_states_identical(&sim, &reference, label);
    }

    // Editing only velocities keeps the cached accelerations valid and must
    // be honoured by the integration.
    let mut sim = NBodyGravity::new(bodies.clone());
    sim.step(dt);
    sim.bodies[0].vel = [10.0, -20.0, 30.0];
    let mut reference = NBodyGravity::new(sim.bodies.clone());
    sim.step(dt);
    reference.step(dt);
    assert_states_identical(&sim, &reference, "velocity-only edit");
}

#[test]
fn empty_and_single_body_steps_are_well_defined() {
    let mut empty = NBodyGravity::new(Vec::new());
    empty.step(10.0);
    assert!(empty.bodies.is_empty());
    assert_eq!(empty.center_of_mass(), [0.0; 3]);

    let mut single = NBodyGravity::with_barnes_hut(
        vec![GravBody::new([1.0, 2.0, 3.0], [4.0, 5.0, 6.0], 7.0)],
        0.5,
    );
    single.step(2.0);
    assert_eq!(single.bodies[0].pos, [9.0, 12.0, 15.0]);
    assert_eq!(single.bodies[0].vel, [4.0, 5.0, 6.0]);
}
