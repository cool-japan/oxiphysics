// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0
//! Regression tests for `contact_dynamics::ImpulseResolver`.
//!
//! The resolver used to compute a normal impulse only: the friction
//! coefficient was dropped and the Baumgarte penetration bias that
//! `PositionCorrection` promises "inside the solver" did not exist.

use oxiphysics_rigid::contact_dynamics::{
    ContactManifold, ContactPoint, ImpulseResolver, PositionCorrection, PositionCorrectionMethod,
    SolverBody, contact_tangent_frame,
};

const G: f64 = 9.81;
const DT: f64 = 1.0 / 240.0;

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// Translating block (rotation locked) of mass `mass` resting on ground at the
/// origin, with the given velocity. Body 0 is the block, body 1 the ground.
fn block_on_ground(mass: f64, vel: [f64; 3]) -> Vec<SolverBody> {
    vec![
        SolverBody::dynamic([0.0, 0.5, 0.0], vel, [0.0; 3], mass, [0.0; 3]),
        SolverBody::static_body([0.0, -0.5, 0.0]),
    ]
}

/// One contact at the origin, normal +y (from the ground, body 1, towards the
/// block, body 0).
fn ground_manifold(
    friction_static: f64,
    friction_dynamic: f64,
    restitution: f64,
    penetration: f64,
) -> ContactManifold {
    let mut manifold = ContactManifold::new(0, 1, friction_static, friction_dynamic, restitution);
    manifold.add_point(ContactPoint::new([0.0; 3], [0.0, 1.0, 0.0], penetration));
    manifold
}

/// One simulation frame: gravity on body 0, then a fresh solve (the
/// accumulated impulses are per-solve quantities).
fn frame(
    resolver: &ImpulseResolver,
    bodies: &mut [SolverBody],
    manifolds: &mut [ContactManifold],
    dt: f64,
) {
    bodies[0].vel[1] -= G * dt;
    for manifold in manifolds.iter_mut() {
        for cp in manifold.points.iter_mut() {
            cp.reset_impulses();
        }
    }
    resolver.solve_manifolds(bodies, manifolds, dt);
}

// ── friction ────────────────────────────────────────────────────────────────

#[test]
fn sliding_block_decelerates_at_mu_times_g() {
    let resolver = ImpulseResolver::default();
    let mut bodies = block_on_ground(2.0, [1.0, 0.0, 0.0]);
    let mut manifolds = vec![ground_manifold(0.8, 0.5, 0.0, 0.0)];
    let steps = 24; // 0.1 s
    for _ in 0..steps {
        frame(&resolver, &mut bodies, &mut manifolds, DT);
    }
    let expected = 1.0 - 0.5 * G * (steps as f64) * DT;
    assert!(
        (bodies[0].vel[0] - expected).abs() < 1e-9,
        "vx = {} vs {expected}",
        bodies[0].vel[0]
    );
    assert!(bodies[0].vel[1].abs() < 1e-12, "vy = {}", bodies[0].vel[1]);
    assert!(bodies[0].vel[2].abs() < 1e-12, "vz = {}", bodies[0].vel[2]);
}

#[test]
fn sliding_block_comes_to_rest_exactly_and_then_sticks() {
    let resolver = ImpulseResolver::default();
    let mut bodies = block_on_ground(1.0, [1.0, 0.0, 0.0]);
    let mut manifolds = vec![ground_manifold(0.8, 0.5, 0.0, 0.0)];
    // mu g t = 1  =>  t = 0.2039 s; run for 0.5 s.
    for _ in 0..120 {
        frame(&resolver, &mut bodies, &mut manifolds, DT);
    }
    assert!(
        bodies[0].vel[0].abs() < 1e-12,
        "block must stop, vx = {}",
        bodies[0].vel[0]
    );
    // Friction never reverses the motion.
    assert!(bodies[0].vel[0] >= -1e-12);
}

#[test]
fn friction_impulse_is_the_coulomb_limit_while_sliding() {
    let resolver = ImpulseResolver::default();
    let mut bodies = block_on_ground(2.0, [3.0, 0.0, 4.0]);
    bodies[0].vel[1] = -1.0;
    let mut manifolds = vec![ground_manifold(0.9, 0.5, 0.0, 0.0)];
    resolver.solve_manifolds(&mut bodies, &mut manifolds, DT);

    let cp = &manifolds[0].points[0];
    // Normal impulse stops the approach: m * 1 = 2.
    assert!((cp.accumulated_normal_impulse - 2.0).abs() < 1e-12);
    // Tangential impulse saturates the Coulomb disc: mu * jn = 1.
    let jt = cp
        .accumulated_tangent1_impulse
        .hypot(cp.accumulated_tangent2_impulse);
    assert!((jt - 1.0).abs() < 1e-12, "|jt| = {jt}");
    // The velocity change opposes the slip direction (0.6, 0, 0.8) exactly:
    // dv = jt / m = 0.5.
    let v = bodies[0].vel;
    assert!((v[0] - (3.0 - 0.5 * 0.6)).abs() < 1e-12, "vx = {}", v[0]);
    assert!((v[2] - (4.0 - 0.5 * 0.8)).abs() < 1e-12, "vz = {}", v[2]);
    assert!(v[1].abs() < 1e-12);
}

#[test]
fn accumulated_tangent_impulses_equal_the_momentum_change() {
    let resolver = ImpulseResolver::default();
    let mut bodies = block_on_ground(2.0, [3.0, 0.0, 4.0]);
    bodies[0].vel[1] = -1.0;
    let mut manifolds = vec![ground_manifold(0.9, 0.5, 0.0, 0.0)];
    resolver.solve_manifolds(&mut bodies, &mut manifolds, DT);
    let cp = &manifolds[0].points[0];
    let (t1, t2) = contact_tangent_frame(cp.normal);
    let dp = [
        2.0 * (bodies[0].vel[0] - 3.0),
        2.0 * (bodies[0].vel[1] + 1.0),
        2.0 * (bodies[0].vel[2] - 4.0),
    ];
    // Normal + tangential impulses reconstruct the momentum change.
    for k in 0..3 {
        let expected = cp.accumulated_normal_impulse * cp.normal[k]
            + cp.accumulated_tangent1_impulse * t1[k]
            + cp.accumulated_tangent2_impulse * t2[k];
        assert!(
            (dp[k] - expected).abs() < 1e-12,
            "component {k}: {} vs {expected}",
            dp[k]
        );
    }
    // The tangential impulse has no normal component.
    assert!(dot(t1, cp.normal).abs() < 1e-12 && dot(t2, cp.normal).abs() < 1e-12);
}

#[test]
fn zero_friction_leaves_the_tangential_velocity_alone() {
    let resolver = ImpulseResolver::default();
    let mut bodies = block_on_ground(1.0, [2.0, 0.0, -1.0]);
    bodies[0].vel[1] = -3.0;
    let mut manifolds = vec![ground_manifold(0.0, 0.0, 0.0, 0.0)];
    resolver.solve_manifolds(&mut bodies, &mut manifolds, DT);
    assert!((bodies[0].vel[0] - 2.0).abs() < 1e-12);
    assert!((bodies[0].vel[2] + 1.0).abs() < 1e-12);
    assert!(manifolds[0].points[0].accumulated_tangent1_impulse.abs() < 1e-15);
    assert!(manifolds[0].points[0].accumulated_tangent2_impulse.abs() < 1e-15);
}

#[test]
fn friction_needs_a_normal_load() {
    // The contact separates (vn > 0): no normal impulse, hence no friction.
    let resolver = ImpulseResolver::default();
    let mut bodies = block_on_ground(1.0, [2.0, 0.0, 0.0]);
    bodies[0].vel[1] = 1.0;
    let mut manifolds = vec![ground_manifold(1.0, 1.0, 0.0, 0.0)];
    resolver.solve_manifolds(&mut bodies, &mut manifolds, DT);
    assert!(manifolds[0].points[0].accumulated_normal_impulse.abs() < 1e-15);
    assert!((bodies[0].vel[0] - 2.0).abs() < 1e-12);
    assert!((bodies[0].vel[1] - 1.0).abs() < 1e-12);
}

#[test]
fn static_friction_holds_a_slow_contact() {
    // Slip 0.005 m/s is below the 0.01 m/s threshold: the static coefficient
    // (0.6) bounds the friction impulse, so a load of m * 0.1 stops the slip
    // although the dynamic coefficient (0.001) alone could not.
    let resolver = ImpulseResolver::default();
    let mut bodies = block_on_ground(1.0, [0.005, 0.0, 0.0]);
    bodies[0].vel[1] = -0.1;
    let mut manifolds = vec![ground_manifold(0.6, 0.001, 0.0, 0.0)];
    resolver.solve_manifolds(&mut bodies, &mut manifolds, DT);
    assert!(bodies[0].vel[0].abs() < 1e-12, "vx = {}", bodies[0].vel[0]);
}

#[test]
fn dynamic_friction_applies_to_a_fast_contact() {
    // Same load, but the 0.05 m/s slip is above the threshold: the dynamic
    // coefficient (0.001) bounds the impulse to 0.001 * 0.1.
    let resolver = ImpulseResolver::default();
    let mut bodies = block_on_ground(1.0, [0.05, 0.0, 0.0]);
    bodies[0].vel[1] = -0.1;
    let mut manifolds = vec![ground_manifold(0.6, 0.001, 0.0, 0.0)];
    resolver.solve_manifolds(&mut bodies, &mut manifolds, DT);
    assert!(
        (bodies[0].vel[0] - (0.05 - 0.001 * 0.1)).abs() < 1e-12,
        "vx = {}",
        bodies[0].vel[0]
    );
}

#[test]
fn friction_exchanges_equal_and_opposite_impulses_between_dynamic_bodies() {
    let resolver = ImpulseResolver::default();
    // A (mass 1) slides over B (mass 3); both translate only.
    let mut bodies = vec![
        SolverBody::dynamic([0.0, 0.5, 0.0], [2.0, -1.0, 0.0], [0.0; 3], 1.0, [0.0; 3]),
        SolverBody::dynamic([0.0, -0.5, 0.0], [0.0, 1.0, 0.0], [0.0; 3], 3.0, [0.0; 3]),
    ];
    let mut manifolds = vec![ground_manifold(0.2, 0.2, 0.0, 0.0)];
    resolver.solve_manifolds(&mut bodies, &mut manifolds, DT);

    // Normal: effective mass 0.75, approach speed 2  =>  jn = 1.5.
    let cp = &manifolds[0].points[0];
    assert!((cp.accumulated_normal_impulse - 1.5).abs() < 1e-12);
    // Friction is capped at mu * jn = 0.3 (the full cancellation would need 1.5).
    assert!(
        (bodies[0].vel[0] - 1.7).abs() < 1e-12,
        "{:?}",
        bodies[0].vel
    );
    assert!(
        (bodies[1].vel[0] - 0.1).abs() < 1e-12,
        "{:?}",
        bodies[1].vel
    );
    // Tangential momentum is conserved.
    let momentum = bodies[0].vel[0] + 3.0 * bodies[1].vel[0];
    assert!((momentum - 2.0).abs() < 1e-12, "momentum = {momentum}");
}

#[test]
fn friction_below_the_cone_cancels_the_slip_completely() {
    // Slip 0.2 with load 1.0 and mu 0.5: the needed impulse (0.2) is inside
    // the cone, so the contact point stops sliding exactly.
    let resolver = ImpulseResolver::default();
    let mut bodies = block_on_ground(1.0, [0.2, 0.0, -0.1]);
    bodies[0].vel[1] = -1.0;
    let mut manifolds = vec![ground_manifold(0.5, 0.5, 0.0, 0.0)];
    resolver.solve_manifolds(&mut bodies, &mut manifolds, DT);
    assert!(bodies[0].vel[0].abs() < 1e-12, "vx = {}", bodies[0].vel[0]);
    assert!(bodies[0].vel[2].abs() < 1e-12, "vz = {}", bodies[0].vel[2]);
    let cp = &manifolds[0].points[0];
    let jt = cp
        .accumulated_tangent1_impulse
        .hypot(cp.accumulated_tangent2_impulse);
    assert!((jt - (0.2_f64).hypot(0.1)).abs() < 1e-12, "|jt| = {jt}");
}

#[test]
fn off_centre_friction_spins_the_body() {
    // A block with inertia 1 slides at 1 m/s; the contact lies 0.5 m below the
    // centre of mass. Friction at the contact (mu 0.3, jn 1) drags the bottom
    // back and tips the block forward (negative rotation about z).
    let resolver = ImpulseResolver::default();
    let mut bodies = vec![
        SolverBody::dynamic([0.0, 0.5, 0.0], [1.0, -1.0, 0.0], [0.0; 3], 1.0, [1.0; 3]),
        SolverBody::static_body([0.0, -0.5, 0.0]),
    ];
    let mut manifolds = vec![ground_manifold(0.3, 0.3, 0.0, 0.0)];
    resolver.solve_manifolds(&mut bodies, &mut manifolds, DT);
    let cp = &manifolds[0].points[0];
    assert!((cp.accumulated_normal_impulse - 1.0).abs() < 1e-12);
    // J = 0.3 along -x: dv = -0.3, d(omega_z) = (r x J)_z / I = -0.5 * 0.3.
    assert!((bodies[0].vel[0] - 0.7).abs() < 1e-9, "{:?}", bodies[0].vel);
    assert!(
        (bodies[0].omega[2] + 0.15).abs() < 1e-9,
        "{:?}",
        bodies[0].omega
    );
}

#[test]
fn rolling_without_slip_produces_no_friction_impulse() {
    // Sphere of radius 0.5, mass 2, I = 2/5 m r^2 = 0.2; v = 1 along x rolls
    // with omega_z = -v / r = -2, so the contact point does not slip.
    let resolver = ImpulseResolver::default();
    let mut bodies = vec![
        SolverBody::dynamic(
            [0.0, 0.5, 0.0],
            [1.0, -0.1, 0.0],
            [0.0, 0.0, -2.0],
            2.0,
            [0.2; 3],
        ),
        SolverBody::static_body([0.0, -0.5, 0.0]),
    ];
    let mut manifolds = vec![ground_manifold(0.8, 0.5, 0.0, 0.0)];
    resolver.solve_manifolds(&mut bodies, &mut manifolds, DT);
    let cp = &manifolds[0].points[0];
    assert!(cp.accumulated_normal_impulse > 0.0);
    assert!(cp.accumulated_tangent1_impulse.abs() < 1e-12);
    assert!(cp.accumulated_tangent2_impulse.abs() < 1e-12);
    assert!((bodies[0].vel[0] - 1.0).abs() < 1e-12);
    assert!((bodies[0].omega[2] + 2.0).abs() < 1e-12);
}

#[test]
fn resolve_contact_applies_friction_for_a_single_contact() {
    let resolver = ImpulseResolver::default();
    let mut ba = SolverBody::dynamic([0.0, 0.5, 0.0], [1.0, -1.0, 0.0], [0.0; 3], 1.0, [0.0; 3]);
    let mut bb = SolverBody::static_body([0.0, -0.5, 0.0]);
    let mut cp = ContactPoint::new([0.0; 3], [0.0, 1.0, 0.0], 0.0);
    resolver.resolve_contact(&mut ba, &mut bb, &mut cp, 0.0, 0.4, DT);
    assert!((ba.vel[0] - 0.6).abs() < 1e-12, "vx = {}", ba.vel[0]);
    assert!(ba.vel[1].abs() < 1e-12);
    assert!((cp.accumulated_normal_impulse - 1.0).abs() < 1e-12);
}

#[test]
fn friction_is_independent_of_the_body_order() {
    // Ground first (index 0), block second (index 1): the manifold lists the
    // block as body A with the normal pointing from the ground to the block.
    let resolver = ImpulseResolver::default();
    let mut bodies = vec![
        SolverBody::static_body([0.0, -0.5, 0.0]),
        SolverBody::dynamic([0.0, 0.5, 0.0], [1.0, -1.0, 0.0], [0.0; 3], 1.0, [0.0; 3]),
    ];
    let mut manifold = ContactManifold::new(1, 0, 0.5, 0.5, 0.0);
    manifold.add_point(ContactPoint::new([0.0; 3], [0.0, 1.0, 0.0], 0.0));
    let mut manifolds = vec![manifold];
    resolver.solve_manifolds(&mut bodies, &mut manifolds, DT);
    assert!(
        (bodies[1].vel[0] - 0.5).abs() < 1e-12,
        "{:?}",
        bodies[1].vel
    );
    assert!(bodies[1].vel[1].abs() < 1e-12);
    assert_eq!(bodies[0].vel, [0.0; 3]);
}

// ── Baumgarte bias ──────────────────────────────────────────────────────────

#[test]
fn penetrating_resting_contact_separates_with_the_baumgarte_velocity() {
    let resolver = ImpulseResolver::default();
    let dt = 1.0 / 60.0;
    let mut bodies = block_on_ground(1.0, [0.0; 3]);
    let mut manifolds = vec![ground_manifold(0.5, 0.5, 0.0, 0.05)];
    resolver.solve_manifolds(&mut bodies, &mut manifolds, dt);
    // bias = beta / dt * (depth - slop) = 0.2 * 60 * 0.045.
    let expected = 0.2 * 60.0 * (0.05 - 0.005);
    assert!(
        (bodies[0].vel[1] - expected).abs() < 1e-12,
        "vy = {} vs {expected}",
        bodies[0].vel[1]
    );
    assert!(bodies[0].vel[1] > 0.0);
}

#[test]
fn baumgarte_bias_is_capped_by_the_maximum_correction() {
    let resolver = ImpulseResolver::default();
    let dt = 1.0 / 60.0;
    let mut bodies = block_on_ground(1.0, [0.0; 3]);
    let mut manifolds = vec![ground_manifold(0.5, 0.5, 0.0, 10.0)];
    resolver.solve_manifolds(&mut bodies, &mut manifolds, dt);
    let cap = resolver.position_correction.max_correction / dt;
    assert!(
        (bodies[0].vel[1] - cap).abs() < 1e-9,
        "vy = {}",
        bodies[0].vel[1]
    );
}

#[test]
fn penetration_within_the_slop_is_not_corrected() {
    let resolver = ImpulseResolver::default();
    let mut bodies = block_on_ground(1.0, [0.0; 3]);
    let mut manifolds = vec![ground_manifold(0.5, 0.5, 0.0, 0.004)];
    resolver.solve_manifolds(&mut bodies, &mut manifolds, 1.0 / 60.0);
    assert!(bodies[0].vel[1].abs() < 1e-15);
}

#[test]
fn only_the_baumgarte_method_biases_the_velocity_solver() {
    for method in [
        PositionCorrectionMethod::None,
        PositionCorrectionMethod::PseudoVelocity,
    ] {
        let resolver = ImpulseResolver {
            position_correction: PositionCorrection {
                method,
                ..PositionCorrection::default()
            },
            ..ImpulseResolver::default()
        };
        let mut bodies = block_on_ground(1.0, [0.0; 3]);
        let mut manifolds = vec![ground_manifold(0.5, 0.5, 0.0, 0.05)];
        resolver.solve_manifolds(&mut bodies, &mut manifolds, 1.0 / 60.0);
        assert!(
            bodies[0].vel[1].abs() < 1e-15,
            "{method:?}: vy = {}",
            bodies[0].vel[1]
        );
    }
}

#[test]
fn baumgarte_bias_needs_a_positive_time_step() {
    let pc = PositionCorrection::default();
    assert_eq!(pc.baumgarte_bias(0.1, 0.0), 0.0);
    assert_eq!(pc.baumgarte_bias(0.1, -1.0), 0.0);
    assert_eq!(pc.baumgarte_bias(0.1, f64::NAN), 0.0);
    assert_eq!(pc.baumgarte_bias(0.0, 0.0), 0.0);
    // A zero time step through the solver leaves the resting contact alone.
    let resolver = ImpulseResolver::default();
    let mut bodies = block_on_ground(1.0, [0.0; 3]);
    let mut manifolds = vec![ground_manifold(0.5, 0.5, 0.0, 0.05)];
    resolver.solve_manifolds(&mut bodies, &mut manifolds, 0.0);
    assert!(bodies[0].vel[1].is_finite());
    assert!(bodies[0].vel[1].abs() < 1e-15);
}

// ── restitution ─────────────────────────────────────────────────────────────

#[test]
fn restitution_bounce_uses_the_pre_solve_velocity_over_all_iterations() {
    // Two contact points of one translating block share a single degree of
    // freedom. With 10 iterations the bounce must still be e * v: evaluating
    // the restitution from the already-bounced velocity would erase it again.
    let resolver = ImpulseResolver {
        restitution_velocity_threshold: 0.0,
        ..ImpulseResolver::default()
    };
    let mut bodies = block_on_ground(1.0, [0.0, -2.0, 0.0]);
    let mut manifold = ContactManifold::new(0, 1, 0.5, 0.5, 0.5);
    manifold.add_point(ContactPoint::new([-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], 0.0));
    manifold.add_point(ContactPoint::new([1.0, 0.0, 0.0], [0.0, 1.0, 0.0], 0.0));
    let mut manifolds = vec![manifold];
    resolver.solve_manifolds(&mut bodies, &mut manifolds, DT);
    assert!(
        (bodies[0].vel[1] - 1.0).abs() < 1e-12,
        "vy = {}",
        bodies[0].vel[1]
    );
}

#[test]
fn restitution_and_baumgarte_bias_add_up() {
    let resolver = ImpulseResolver {
        restitution_velocity_threshold: 0.0,
        ..ImpulseResolver::default()
    };
    let dt = 1.0 / 60.0;
    let mut bodies = block_on_ground(1.0, [0.0, -2.0, 0.0]);
    let mut manifolds = vec![ground_manifold(0.5, 0.5, 0.5, 0.05)];
    resolver.solve_manifolds(&mut bodies, &mut manifolds, dt);
    let expected = 0.5 * 2.0 + 0.2 * 60.0 * (0.05 - 0.005);
    assert!(
        (bodies[0].vel[1] - expected).abs() < 1e-12,
        "vy = {} vs {expected}",
        bodies[0].vel[1]
    );
}

#[test]
fn slow_impacts_do_not_bounce() {
    // Approach 0.3 m/s is below the default 0.5 m/s restitution threshold.
    let resolver = ImpulseResolver::default();
    let mut bodies = block_on_ground(1.0, [0.0, -0.3, 0.0]);
    let mut manifolds = vec![ground_manifold(0.5, 0.5, 1.0, 0.0)];
    resolver.solve_manifolds(&mut bodies, &mut manifolds, DT);
    assert!(bodies[0].vel[1].abs() < 1e-12, "vy = {}", bodies[0].vel[1]);
}

#[test]
fn separating_contact_receives_no_impulse() {
    let resolver = ImpulseResolver::default();
    let mut bodies = block_on_ground(1.0, [0.4, 2.0, 0.0]);
    let mut manifolds = vec![ground_manifold(0.9, 0.9, 0.5, 0.0)];
    resolver.solve_manifolds(&mut bodies, &mut manifolds, DT);
    assert_eq!(bodies[0].vel, [0.4, 2.0, 0.0]);
    assert_eq!(manifolds[0].points[0].accumulated_normal_impulse, 0.0);
}

// ── robustness ──────────────────────────────────────────────────────────────

#[test]
fn invalid_manifolds_are_skipped_without_panicking() {
    let resolver = ImpulseResolver::default();
    let mut bodies = block_on_ground(1.0, [1.0, -1.0, 0.0]);
    let before = (bodies[0].vel, bodies[1].vel);

    // Body index out of range.
    let mut missing = ContactManifold::new(0, 7, 0.5, 0.5, 0.0);
    missing.add_point(ContactPoint::new([0.0; 3], [0.0, 1.0, 0.0], 0.1));
    // A body in contact with itself.
    let mut itself = ContactManifold::new(0, 0, 0.5, 0.5, 0.0);
    itself.add_point(ContactPoint::new([0.0; 3], [0.0, 1.0, 0.0], 0.1));
    let mut manifolds = vec![missing, itself];

    resolver.solve_manifolds(&mut bodies, &mut manifolds, DT);
    assert_eq!((bodies[0].vel, bodies[1].vel), before);
    for manifold in &manifolds {
        let cp = &manifold.points[0];
        assert_eq!(cp.accumulated_normal_impulse, 0.0);
        assert_eq!(cp.accumulated_tangent1_impulse, 0.0);
    }

    // The valid manifold next to the invalid ones is still solved.
    manifolds.push(ground_manifold(0.5, 0.5, 0.0, 0.0));
    resolver.solve_manifolds(&mut bodies, &mut manifolds, DT);
    assert!(bodies[0].vel[1].abs() < 1e-12);
    assert!((bodies[0].vel[0] - 0.5).abs() < 1e-12);
}

#[test]
fn two_static_bodies_stay_put() {
    let resolver = ImpulseResolver::default();
    let mut bodies = vec![
        SolverBody::static_body([0.0, 0.5, 0.0]),
        SolverBody::static_body([0.0, -0.5, 0.0]),
    ];
    let mut manifolds = vec![ground_manifold(0.5, 0.5, 0.5, 0.2)];
    resolver.solve_manifolds(&mut bodies, &mut manifolds, DT);
    assert_eq!(bodies[0].vel, [0.0; 3]);
    assert_eq!(bodies[1].vel, [0.0; 3]);
}

#[test]
fn zero_iterations_change_nothing() {
    let resolver = ImpulseResolver {
        iterations: 0,
        ..ImpulseResolver::default()
    };
    let mut bodies = block_on_ground(1.0, [1.0, -1.0, 0.0]);
    let mut manifolds = vec![ground_manifold(0.5, 0.5, 0.0, 0.1)];
    resolver.solve_manifolds(&mut bodies, &mut manifolds, DT);
    assert_eq!(bodies[0].vel, [1.0, -1.0, 0.0]);
}

#[test]
fn friction_never_exceeds_the_cone_for_random_contacts() {
    // A deterministic sweep over slip directions and loads.
    let resolver = ImpulseResolver {
        iterations: 20,
        ..ImpulseResolver::default()
    };
    let mu = 0.35;
    for i in 0..64 {
        let angle = 0.1 * f64::from(i);
        let speed = 0.05 + 0.2 * f64::from(i % 7);
        let approach = 0.1 + 0.3 * f64::from(i % 5);
        let mut bodies = vec![
            SolverBody::dynamic(
                [0.0, 0.5, 0.0],
                [speed * angle.cos(), -approach, speed * angle.sin()],
                [0.3, -0.2, 0.4],
                1.5,
                [0.4, 0.5, 0.6],
            ),
            SolverBody::static_body([0.0, -0.5, 0.0]),
        ];
        let mut manifolds = vec![ground_manifold(mu, mu, 0.0, 0.0)];
        resolver.solve_manifolds(&mut bodies, &mut manifolds, DT);
        let cp = &manifolds[0].points[0];
        let jt = cp
            .accumulated_tangent1_impulse
            .hypot(cp.accumulated_tangent2_impulse);
        assert!(
            jt <= mu * cp.accumulated_normal_impulse + 1e-12,
            "case {i}: |jt| = {jt} exceeds mu * jn = {}",
            mu * cp.accumulated_normal_impulse
        );
        assert!(cp.accumulated_normal_impulse >= 0.0);
    }
}
