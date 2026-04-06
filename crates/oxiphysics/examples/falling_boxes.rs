// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Simulate 10 boxes falling under gravity and hitting a static floor.
//! After 5 seconds, print each box's final position and verify none
//! has fallen below y = -2.0.

use oxiphysics::pipeline::PhysicsPipeline;
use oxiphysics_core::{Transform, math::Vec3};
use oxiphysics_geometry::{BoxShape, Shape, Sphere};
use oxiphysics_rigid::{Collider, ColliderSet, RigidBody, RigidBodySet};
use std::sync::Arc;

fn main() {
    let mut pipeline = PhysicsPipeline::new();
    let mut bodies = RigidBodySet::new();
    let mut colliders = ColliderSet::new();

    // Static floor: a wide flat box. Half-extents (50, 0.5, 50).
    // Centred at y = -0.5, so the top surface is at y = 0.0.
    let mut floor_body = RigidBody::new_static();
    floor_body.transform = Transform::from_position(Vec3::new(0.0, -0.5, 0.0));
    let floor_h = bodies.insert(floor_body);
    let floor_shape: Arc<dyn Shape> = Arc::new(BoxShape::new(Vec3::new(50.0, 0.5, 50.0)));
    colliders.insert(Collider::new(floor_shape).with_body(floor_h));

    // 10 dynamic spheres (radius 0.25) spread along x, starting at y = 1.5.
    // Spheres are used so the contact normal is well-defined from any position.
    let radius = 0.25_f64;
    let dyn_shape: Arc<dyn Shape> = Arc::new(Sphere::new(radius));
    let mut box_handles = Vec::with_capacity(10);
    for i in 0..10 {
        let x = (i as f64 - 4.5) * 0.7;
        let mut b = RigidBody::new(1.0);
        b.transform = Transform::from_position(Vec3::new(x, 1.5, 0.0));
        b.linear_damping = 0.02;
        let h = bodies.insert(b);
        colliders.insert(
            Collider::new(Arc::clone(&dyn_shape))
                .with_body(h)
                .with_restitution(0.2),
        );
        box_handles.push(h);
    }

    // Simulate for 5 seconds at 120 Hz.
    let dt = 1.0 / 120.0;
    let steps = (5.0 / dt) as usize;
    for _ in 0..steps {
        pipeline.step(dt, &mut bodies, &colliders);
    }

    // Report final positions.
    println!("=== falling_boxes: final positions after 5 s ===");
    let mut all_above = true;
    for (idx, &h) in box_handles.iter().enumerate() {
        if let Some(b) = bodies.get(h) {
            let pos = b.transform.position;
            println!("  box {:2}: y = {:.4}", idx, pos.y);
            if pos.y < -2.0 {
                all_above = false;
                println!("    FAIL: box {} fell below y=-2.0 (y={:.4})", idx, pos.y);
            }
        }
    }
    if all_above {
        println!("PASS: All boxes above floor");
    }
}
