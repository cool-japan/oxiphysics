// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0
//! Regression tests for the collider-driven contact pipeline of
//! [`PhysicsWorld`]: contacts come from the attached collider shapes (never
//! from body mass), approaching contacts are resolved while separating ones
//! are left alone, and the collider restitution / friction reach the solver.

use oxiphysics_core::BodyHandle;
use oxiphysics_core::Transform;
use oxiphysics_core::math::Vec3;
use oxiphysics_geometry::{BoxShape, Compound, HeightField, Shape, Sphere, Torus, TriangleMesh};
use oxiphysics_rigid::world::PhysicsWorld;
use oxiphysics_rigid::{BodyState, Collider, RigidBody};
use std::sync::Arc;

/// Add a dynamic body at `pos` with collider `shape` (inertia from the shape).
fn add_dynamic(
    world: &mut PhysicsWorld,
    mass: f64,
    pos: [f64; 3],
    shape: Arc<dyn Shape>,
    restitution: f64,
    friction: f64,
) -> BodyHandle {
    let mut body = RigidBody::new(mass);
    body.transform.position = Vec3::new(pos[0], pos[1], pos[2]);
    body.set_inertia_from_shape(shape.as_ref());
    body.update_world_inertia();
    let h = world.add_rigid_body(body);
    let collider = Collider::new(shape)
        .with_restitution(restitution)
        .with_friction(friction);
    world.add_collider(collider, Some(h));
    h
}

/// Add a static body at `pos` with collider `shape`.
fn add_static(
    world: &mut PhysicsWorld,
    pos: [f64; 3],
    shape: Arc<dyn Shape>,
    restitution: f64,
    friction: f64,
) -> BodyHandle {
    let mut body = RigidBody::new_static();
    body.transform.position = Vec3::new(pos[0], pos[1], pos[2]);
    let h = world.add_rigid_body(body);
    let collider = Collider::new(shape)
        .with_restitution(restitution)
        .with_friction(friction);
    world.add_collider(collider, Some(h));
    h
}

/// Static ground box whose top face is the plane `y = 0`.
fn add_ground(world: &mut PhysicsWorld, restitution: f64, friction: f64) -> BodyHandle {
    add_static(
        world,
        [0.0, -0.5, 0.0],
        Arc::new(BoxShape::new(Vec3::new(20.0, 0.5, 20.0))),
        restitution,
        friction,
    )
}

fn unit_box() -> Arc<dyn Shape> {
    Arc::new(BoxShape::new(Vec3::new(0.5, 0.5, 0.5)))
}

fn body_y(world: &PhysicsWorld, h: BodyHandle) -> f64 {
    world
        .get_body(h)
        .expect("body inserted by the test")
        .transform
        .position
        .y
}

#[test]
fn box_rests_on_static_box_with_four_point_manifold() {
    let mut world = PhysicsWorld::new([0.0, -9.81, 0.0], 1.0 / 120.0);
    world.linear_sleep_threshold = 0.0;
    world.angular_sleep_threshold = 0.0;
    add_ground(&mut world, 0.0, 0.5);
    let b = add_dynamic(&mut world, 1.0, [0.0, 0.8, 0.0], unit_box(), 0.0, 0.5);
    for _ in 0..240 {
        world.step();
    }
    let y = body_y(&world, b);
    assert!((y - 0.5).abs() < 1e-3, "box must rest at y = 0.5, got {y}");
    let contacts: Vec<_> = world
        .last_contacts
        .iter()
        .filter(|c| c.body_a == b || c.body_b == b)
        .collect();
    assert_eq!(contacts.len(), 4, "face-on-face box contact has 4 corners");
    for c in contacts {
        let n_y = if c.body_a == b {
            c.normal[1]
        } else {
            -c.normal[1]
        };
        assert!(
            n_y > 0.99,
            "normal must push the box up, got {:?}",
            c.normal
        );
        assert!(c.depth >= 0.0 && c.depth < 2e-3, "depth {}", c.depth);
    }
    let body = world.get_body(b).expect("box");
    assert!(
        body.velocity.norm() < 1e-2,
        "box at rest: v = {}",
        body.velocity
    );
    assert!(
        body.angular_velocity.norm() < 1e-2,
        "box must not rock: w = {}",
        body.angular_velocity
    );
}

#[test]
fn bodies_without_colliders_never_collide() {
    let mut world = PhysicsWorld::new([0.0, -9.81, 0.0], 1.0 / 60.0);
    // Heavy bodies used to be treated as `mass.cbrt() * 0.5` spheres (5 m for
    // 1000 kg); without colliders they must pass through each other.
    let mut ground = RigidBody::new_static();
    ground.mass = 1000.0;
    world.add_rigid_body(ground);
    let falling = world.add_rigid_body(RigidBody::new(1000.0));
    for _ in 0..30 {
        world.step();
    }
    assert!(world.last_contacts.is_empty(), "no collider, no contact");
    let y = body_y(&world, falling);
    let t = 30.0 / 60.0;
    let expected = -0.5 * 9.81 * t * (t + 1.0 / 60.0);
    assert!(
        (y - expected).abs() < 1e-9,
        "free fall expected y = {expected}, got {y}"
    );
}

#[test]
fn contact_size_comes_from_the_shape_not_the_mass() {
    let mut world = PhysicsWorld::new([0.0, 0.0, 0.0], 1.0 / 60.0);
    let small_sphere: Arc<dyn Shape> = Arc::new(Sphere::new(0.1));
    let a = add_dynamic(&mut world, 1000.0, [0.0, 0.0, 0.0], small_sphere, 0.0, 0.5);
    let b = add_dynamic(
        &mut world,
        1000.0,
        [1.0, 0.0, 0.0],
        Arc::new(Sphere::new(0.1)),
        0.0,
        0.5,
    );
    world.step();
    assert!(
        world.last_contacts.is_empty(),
        "two 0.1 m spheres 1 m apart must not touch"
    );
    // Now give them colliders that do overlap.
    let mut world2 = PhysicsWorld::new([0.0, 0.0, 0.0], 1.0 / 60.0);
    add_dynamic(&mut world2, 0.001, [0.0, 0.0, 0.0], unit_box(), 0.0, 0.5);
    add_dynamic(&mut world2, 0.001, [0.95, 0.0, 0.0], unit_box(), 0.0, 0.5);
    world2.step();
    assert!(
        !world2.last_contacts.is_empty(),
        "overlapping 1 m boxes must touch even when nearly massless"
    );
    let _ = (a, b);
}

#[test]
fn approaching_contact_is_stopped_and_separating_contact_is_untouched() {
    // Separating: overlapping spheres moving apart keep their velocities.
    let mut world = PhysicsWorld::new([0.0, 0.0, 0.0], 1.0 / 60.0);
    let a = add_dynamic(
        &mut world,
        1.0,
        [0.0, 0.0, 0.0],
        Arc::new(Sphere::new(0.5)),
        0.0,
        0.0,
    );
    let b = add_dynamic(
        &mut world,
        1.0,
        [0.9, 0.0, 0.0],
        Arc::new(Sphere::new(0.5)),
        0.0,
        0.0,
    );
    world.solver_config.position_iterations = 0;
    world.get_body_mut(a).expect("a").velocity.x = -1.0;
    world.get_body_mut(b).expect("b").velocity.x = 1.0;
    world.step();
    assert!(!world.last_contacts.is_empty(), "spheres overlap");
    assert!((world.get_body(a).expect("a").velocity.x + 1.0).abs() < 1e-12);
    assert!((world.get_body(b).expect("b").velocity.x - 1.0).abs() < 1e-12);

    // Approaching: inelastic head-on contact removes the closing velocity.
    let mut world = PhysicsWorld::new([0.0, 0.0, 0.0], 1.0 / 60.0);
    let a = add_dynamic(
        &mut world,
        1.0,
        [0.0, 0.0, 0.0],
        Arc::new(Sphere::new(0.5)),
        0.0,
        0.0,
    );
    let b = add_dynamic(
        &mut world,
        1.0,
        [0.99, 0.0, 0.0],
        Arc::new(Sphere::new(0.5)),
        0.0,
        0.0,
    );
    world.get_body_mut(a).expect("a").velocity.x = 2.0;
    world.get_body_mut(b).expect("b").velocity.x = -2.0;
    world.step();
    let va = world.get_body(a).expect("a").velocity.x;
    let vb = world.get_body(b).expect("b").velocity.x;
    assert!(va.abs() < 1e-9 && vb.abs() < 1e-9, "va = {va}, vb = {vb}");
}

#[test]
fn static_first_contact_still_pushes_the_dynamic_body() {
    // The static ground is inserted first so it is pair member A.
    let mut world = PhysicsWorld::new([0.0, -9.81, 0.0], 1.0 / 120.0);
    world.linear_sleep_threshold = 0.0;
    world.angular_sleep_threshold = 0.0;
    let ground = add_ground(&mut world, 0.0, 0.5);
    let ball = add_dynamic(
        &mut world,
        1.0,
        [0.0, 2.0, 0.0],
        Arc::new(Sphere::new(0.5)),
        0.0,
        0.5,
    );
    for _ in 0..360 {
        world.step();
    }
    assert!(
        world.last_contacts.iter().any(|c| c.body_a == ground),
        "the static ground must appear as body A"
    );
    let y = body_y(&world, ball);
    assert!((y - 0.5).abs() < 2e-3, "ball rests on the ground: y = {y}");
}

#[test]
fn collider_restitution_sets_the_bounce_height() {
    // Drop from 1 m above rest onto ground; both colliders e = 0.5 → h' = 0.25.
    let mut world = PhysicsWorld::new([0.0, -9.81, 0.0], 1.0 / 480.0);
    add_ground(&mut world, 0.5, 0.0);
    let ball = add_dynamic(
        &mut world,
        1.0,
        [0.0, 1.5, 0.0],
        Arc::new(Sphere::new(0.5)),
        0.5,
        0.0,
    );
    let mut bounced = false;
    let mut peak = f64::NEG_INFINITY;
    for _ in 0..800 {
        world.step();
        let b = world.get_body(ball).expect("ball");
        if b.velocity.y > 0.0 {
            bounced = true;
        }
        if bounced {
            peak = peak.max(b.transform.position.y);
            if b.velocity.y < 0.0 {
                break;
            }
        }
    }
    let rebound = peak - 0.5;
    assert!(bounced, "ball must bounce");
    assert!(
        (rebound - 0.25).abs() < 0.02,
        "rebound {rebound} must match e² · h = 0.25"
    );
}

#[test]
fn restitution_combines_with_max_rule() {
    // Ground e = 0, ball e = 0.5 → combined e = 0.5.
    let mut world = PhysicsWorld::new([0.0, -9.81, 0.0], 1.0 / 480.0);
    add_ground(&mut world, 0.0, 0.0);
    add_dynamic(
        &mut world,
        1.0,
        [0.0, 0.52, 0.0],
        Arc::new(Sphere::new(0.5)),
        0.5,
        0.0,
    );
    for _ in 0..200 {
        world.step();
        if !world.last_contacts.is_empty() {
            break;
        }
    }
    let c = world
        .last_contacts
        .first()
        .expect("ball reaches the ground");
    assert!((c.restitution - 0.5).abs() < 1e-12);
    assert!(c.friction.abs() < 1e-12, "sqrt(0 · 0) = 0");
}

#[test]
fn coulomb_friction_decelerates_a_sliding_box() {
    // μ_a = μ_b = 0.5 → combined μ = 0.5; deceleration μ g, stop distance v²/(2μg).
    let mut world = PhysicsWorld::new([0.0, -9.81, 0.0], 1.0 / 240.0);
    world.linear_sleep_threshold = 0.0;
    world.angular_sleep_threshold = 0.0;
    add_ground(&mut world, 0.0, 0.5);
    let b = add_dynamic(&mut world, 1.0, [0.0, 0.5, 0.0], unit_box(), 0.0, 0.5);
    for _ in 0..60 {
        world.step();
    }
    let x0 = world.get_body(b).expect("box").transform.position.x;
    world.get_body_mut(b).expect("box").velocity.x = 2.0;
    for _ in 0..480 {
        world.step();
    }
    let body = world.get_body(b).expect("box");
    let travelled = body.transform.position.x - x0;
    let expected = 2.0 * 2.0 / (2.0 * 0.5 * 9.81);
    assert!(
        body.velocity.x.abs() < 1e-3,
        "box must stop: vx = {}",
        body.velocity.x
    );
    assert!(
        (travelled - expected).abs() < 0.1 * expected,
        "slide distance {travelled} vs analytic {expected}"
    );
    assert!(
        (body.transform.position.y - 0.5).abs() < 5e-3,
        "box stays on the ground"
    );
}

#[test]
fn frictionless_box_keeps_sliding() {
    let mut world = PhysicsWorld::new([0.0, -9.81, 0.0], 1.0 / 240.0);
    world.linear_sleep_threshold = 0.0;
    world.angular_sleep_threshold = 0.0;
    add_ground(&mut world, 0.0, 0.0);
    let b = add_dynamic(&mut world, 1.0, [0.0, 0.5, 0.0], unit_box(), 0.0, 0.0);
    world.get_body_mut(b).expect("box").velocity.x = 2.0;
    for _ in 0..240 {
        world.step();
    }
    let vx = world.get_body(b).expect("box").velocity.x;
    assert!((vx - 2.0).abs() < 1e-6, "no friction, vx stays 2: {vx}");
}

#[test]
fn sphere_rests_on_height_field() {
    let mut world = PhysicsWorld::new([0.0, -9.81, 0.0], 1.0 / 120.0);
    world.linear_sleep_threshold = 0.0;
    world.angular_sleep_threshold = 0.0;
    let hf = HeightField::new(vec![0.0; 9 * 9], 9, 9, 1.0, 1.0);
    add_static(&mut world, [-4.0, 0.0, -4.0], Arc::new(hf), 0.0, 0.5);
    let ball = add_dynamic(
        &mut world,
        1.0,
        [0.3, 1.0, 0.2],
        Arc::new(Sphere::new(0.5)),
        0.0,
        0.5,
    );
    for _ in 0..360 {
        world.step();
    }
    let y = body_y(&world, ball);
    assert!((y - 0.5).abs() < 3e-3, "ball on height field: y = {y}");
}

#[test]
fn box_rests_on_triangle_mesh_floor() {
    let mut world = PhysicsWorld::new([0.0, -9.81, 0.0], 1.0 / 120.0);
    world.linear_sleep_threshold = 0.0;
    world.angular_sleep_threshold = 0.0;
    let mesh = TriangleMesh::new(
        vec![
            Vec3::new(-5.0, 0.0, -5.0),
            Vec3::new(5.0, 0.0, -5.0),
            Vec3::new(5.0, 0.0, 5.0),
            Vec3::new(-5.0, 0.0, 5.0),
        ],
        vec![[0, 2, 1], [0, 3, 2]],
    );
    add_static(&mut world, [0.0, 0.0, 0.0], Arc::new(mesh), 0.0, 0.5);
    let b = add_dynamic(&mut world, 1.0, [0.2, 1.0, 0.1], unit_box(), 0.0, 0.5);
    for _ in 0..360 {
        world.step();
    }
    let y = body_y(&world, b);
    assert!((y - 0.5).abs() < 3e-3, "box on mesh floor: y = {y}");
}

#[test]
fn compound_children_all_collide() {
    let mut world = PhysicsWorld::new([0.0, -9.81, 0.0], 1.0 / 120.0);
    world.linear_sleep_threshold = 0.0;
    world.angular_sleep_threshold = 0.0;
    add_ground(&mut world, 0.0, 0.5);
    let dumbbell = Compound::new(vec![
        (
            Transform::from_position(Vec3::new(-1.0, 0.0, 0.0)),
            Arc::new(Sphere::new(0.25)) as Arc<dyn Shape>,
        ),
        (
            Transform::from_position(Vec3::new(1.0, 0.0, 0.0)),
            Arc::new(Sphere::new(0.25)) as Arc<dyn Shape>,
        ),
    ]);
    let h = add_dynamic(
        &mut world,
        1.0,
        [0.0, 1.0, 0.0],
        Arc::new(dumbbell),
        0.0,
        0.5,
    );
    for _ in 0..360 {
        world.step();
    }
    let y = body_y(&world, h);
    assert!(
        (y - 0.25).abs() < 3e-3,
        "dumbbell rests on both balls: y = {y}"
    );
    let touching = world
        .last_contacts
        .iter()
        .filter(|c| c.body_a == h || c.body_b == h)
        .count();
    assert_eq!(touching, 2, "one contact per child sphere");
}

#[test]
fn ball_falls_through_torus_hole_but_rests_on_its_ring() {
    let torus: Arc<dyn Shape> = Arc::new(Torus::new(2.0, 0.5));
    // Through the hole: radius 0.5 ball along the axis clears the 1.5 m hole.
    let mut world = PhysicsWorld::new([0.0, -9.81, 0.0], 1.0 / 120.0);
    add_static(&mut world, [0.0, 0.0, 0.0], Arc::clone(&torus), 0.0, 0.5);
    let ball = add_dynamic(
        &mut world,
        1.0,
        [0.0, 2.0, 0.0],
        Arc::new(Sphere::new(0.5)),
        0.0,
        0.5,
    );
    for _ in 0..120 {
        world.step();
    }
    assert!(
        body_y(&world, ball) < -1.0,
        "ball must fall through the torus hole (convex-hull treatment would stop it)"
    );
    // On the ring: dropped above the tube centre line it comes to rest on top.
    let mut world = PhysicsWorld::new([0.0, -9.81, 0.0], 1.0 / 120.0);
    world.linear_sleep_threshold = 0.0;
    world.angular_sleep_threshold = 0.0;
    add_static(&mut world, [0.0, 0.0, 0.0], torus, 0.0, 1.0);
    let ball = add_dynamic(
        &mut world,
        1.0,
        [2.0, 1.5, 0.0],
        Arc::new(Sphere::new(0.5)),
        0.0,
        1.0,
    );
    for _ in 0..240 {
        world.step();
    }
    let y = body_y(&world, ball);
    // Tube top at y = 0.5 → ball centre at 1.0 (facet tolerance ≈ 3 mm).
    assert!((y - 1.0).abs() < 1e-2, "ball rests on the ring: y = {y}");
}

#[test]
fn sensors_produce_no_contacts() {
    let mut world = PhysicsWorld::new([0.0, -9.81, 0.0], 1.0 / 60.0);
    let ground = world.add_rigid_body(RigidBody::new_static());
    world.add_collider(
        Collider::new(Arc::new(BoxShape::new(Vec3::new(5.0, 0.5, 5.0)))).as_sensor(),
        Some(ground),
    );
    let ball = add_dynamic(
        &mut world,
        1.0,
        [0.0, 0.9, 0.0],
        Arc::new(Sphere::new(0.5)),
        0.0,
        0.5,
    );
    for _ in 0..30 {
        world.step();
    }
    assert!(world.last_contacts.is_empty(), "sensor exerts no force");
    assert!(body_y(&world, ball) < 0.0, "ball falls through the sensor");
}

#[test]
fn sleeping_body_is_woken_by_an_impact() {
    let mut world = PhysicsWorld::new([0.0, -9.81, 0.0], 1.0 / 120.0);
    add_ground(&mut world, 0.0, 0.5);
    let bottom = add_dynamic(&mut world, 1.0, [0.0, 0.5, 0.0], unit_box(), 0.0, 0.5);
    for _ in 0..240 {
        world.step();
    }
    assert_eq!(
        world.get_body(bottom).expect("bottom").state,
        BodyState::Sleeping,
        "resting box falls asleep"
    );
    let top = add_dynamic(&mut world, 1.0, [0.0, 2.5, 0.0], unit_box(), 0.0, 0.5);
    for _ in 0..240 {
        world.step();
    }
    let y_top = body_y(&world, top);
    let y_bottom = body_y(&world, bottom);
    assert!(
        (y_top - 1.5).abs() < 5e-3,
        "top box lands on the bottom box: {y_top}"
    );
    assert!(
        (y_bottom - 0.5).abs() < 5e-3,
        "bottom box holds: {y_bottom}"
    );
}

#[test]
fn small_steps_settle_a_box_on_the_ground() {
    let mut world = PhysicsWorld::new([0.0, -9.81, 0.0], 1.0 / 60.0);
    world.solver_config.substeps = 4;
    world.linear_sleep_threshold = 0.0;
    world.angular_sleep_threshold = 0.0;
    add_ground(&mut world, 0.0, 0.5);
    let b = add_dynamic(&mut world, 1.0, [0.0, 0.7, 0.0], unit_box(), 0.0, 0.5);
    for _ in 0..120 {
        world.step_small_steps();
    }
    let body = world.get_body(b).expect("box");
    let y = body.transform.position.y;
    assert!((y - 0.5).abs() < 1e-2, "box settles on the ground: y = {y}");
    assert!(
        body.velocity.norm() < 0.05,
        "box at rest: {}",
        body.velocity
    );
    assert!(
        body.angular_velocity.norm() < 0.05,
        "no spin: {}",
        body.angular_velocity
    );
}

#[test]
fn speculative_contacts_report_negative_depth_from_collider_distance() {
    let mut world = PhysicsWorld::new([0.0, 0.0, 0.0], 1.0 / 60.0);
    world.solver_config.use_speculative = true;
    world.solver_config.speculative_margin = 0.2;
    add_ground(&mut world, 0.0, 0.5);
    let b = add_dynamic(&mut world, 1.0, [0.0, 0.6, 0.0], unit_box(), 0.0, 0.5);
    world.get_body_mut(b).expect("box").velocity.y = -12.0;
    world.step();
    let y = body_y(&world, b);
    assert!(
        y >= 0.5 - 1e-4,
        "speculative contact prevents tunnelling: {y}"
    );
    let c = world.last_contacts.first().expect("speculative contact");
    assert!(
        (c.depth + 0.1).abs() < 1e-6,
        "gap of 0.1 m reported as depth -0.1, got {}",
        c.depth
    );
}

#[test]
fn three_box_stack_settles_without_drift_or_spin() {
    // Block-solved 4-point manifolds + warm starting keep a symmetric stack
    // symmetric: no rocking, no lateral drift, bounded penetration.
    let mut world = PhysicsWorld::new([0.0, -9.81, 0.0], 1.0 / 100.0);
    add_ground(&mut world, 0.3, 0.5);
    let rest = [0.5_f64, 1.5, 2.5];
    let handles: Vec<_> = rest
        .iter()
        .map(|&y| add_dynamic(&mut world, 1.0, [0.0, y + 0.5, 0.0], unit_box(), 0.3, 0.5))
        .collect();
    for _ in 0..300 {
        world.step();
    }
    for (h, &y_rest) in handles.iter().zip(rest.iter()) {
        let b = world.get_body(*h).expect("stacked box");
        let p = b.transform.position;
        assert!(y_rest - p.y < 2e-3, "penetration too deep: y = {}", p.y);
        assert!((p.x * p.x + p.z * p.z).sqrt() < 1e-4, "lateral drift: {p}");
        assert!(
            0.5 * b.mass * b.velocity.norm_squared() < 1e-3,
            "not at rest"
        );
        assert!(
            b.angular_velocity.norm() < 1e-3,
            "spinning: {}",
            b.angular_velocity
        );
    }
}
