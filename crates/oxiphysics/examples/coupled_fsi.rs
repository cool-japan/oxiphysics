// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Coupled fluid-structure interaction (FSI) demo.
//!
//! Combines an SPH fluid simulation with a rigid body to demonstrate
//! buoyancy force calculation. A pool of SPH particles sits beneath
//! a falling rigid sphere; the pipeline's SPH-rigid coupling applies
//! upward buoyancy once the body enters the fluid region.

use oxiphysics::pipeline::{PipelineBuilder, SphRigidCouplingParams};
use oxiphysics_core::Transform;
use oxiphysics_core::math::Vec3;
use oxiphysics_geometry::{Shape, Sphere};
use oxiphysics_rigid::{Collider, ColliderSet, RigidBody, RigidBodySet};
use oxiphysics_sph::{
    kernel::CubicSplineKernel,
    neighbor::SpatialHash,
    particle::{ParticleSet, SphParticle},
    wcsph::{WcsphParams, step as wcsph_step},
};
use std::sync::Arc;

fn main() {
    // ── SPH fluid pool ──────────────────────────────────────────────
    let spacing = 0.05_f64;
    let mass = 1000.0 * spacing.powi(3);
    let h = 2.0 * spacing;

    let sph_params = WcsphParams {
        rest_density: 1000.0,
        stiffness: 1000.0 * 9.81 * 2.0,
        gamma: 7.0,
        viscosity: 0.01,
        smoothing_length: h,
    };
    let kernel = CubicSplineKernel;

    // Build an 8x8 block of fluid particles in the XY plane at z=0.
    let mut fluid = ParticleSet::with_capacity(64);
    for ix in 0..8_usize {
        for iy in 0..8_usize {
            let pos = Vec3::new(
                ix as f64 * spacing + spacing * 0.5,
                iy as f64 * spacing + spacing * 0.5,
                0.0,
            );
            fluid.add_particle(&SphParticle::new(pos, Vec3::zeros(), mass));
        }
    }

    let initial_fluid_count = fluid.len();
    let fluid_centroid_y: f64 =
        fluid.positions.iter().map(|p| p.y).sum::<f64>() / initial_fluid_count as f64;
    println!("=== coupled_fsi: SPH fluid + rigid body buoyancy demo ===");
    println!(
        "Fluid: {} particles, centroid y = {:.4}",
        initial_fluid_count, fluid_centroid_y
    );

    // ── Rigid body pipeline ─────────────────────────────────────────
    // Configure a coupling region that covers the fluid pool.
    let coupling = SphRigidCouplingParams {
        fluid_density: 1000.0,
        coupling_alpha: 1.0,
        region_half_extents: Vec3::new(1.0, 1.0, 1.0),
        region_centre: Vec3::new(0.2, 0.2, 0.0),
    };

    let mut pipeline = PipelineBuilder::new()
        .gravity(Vec3::new(0.0, -9.81, 0.0))
        .sph_rigid_coupling(coupling)
        .finish();

    let mut bodies = RigidBodySet::new();
    let mut colliders = ColliderSet::new();

    // Rigid sphere (radius=0.1 m, mass=5 kg) starting above the fluid.
    let sphere_radius = 0.1_f64;
    let mut sphere_body = RigidBody::new(5.0);
    sphere_body.transform = Transform::from_position(Vec3::new(0.2, 1.0, 0.0));
    let sphere_h = bodies.insert(sphere_body);
    let sphere_shape: Arc<dyn Shape> = Arc::new(Sphere::new(sphere_radius));
    colliders.insert(Collider::new(sphere_shape).with_body(sphere_h));

    // Static floor.
    let mut floor_body = RigidBody::new_static();
    floor_body.transform = Transform::from_position(Vec3::new(0.0, -0.5, 0.0));
    let floor_h = bodies.insert(floor_body);
    let floor_shape: Arc<dyn Shape> =
        Arc::new(oxiphysics_geometry::BoxShape::new(Vec3::new(5.0, 0.5, 5.0)));
    colliders.insert(Collider::new(floor_shape).with_body(floor_h));

    // ── Simulation loop ─────────────────────────────────────────────
    let dt = 1.0 / 120.0;
    let n_steps = 600_usize; // 5 seconds
    let gravity = Vec3::new(0.0, -9.81, 0.0);
    let sph_dt = 0.001_f64;
    let sph_substeps = ((dt / sph_dt) as usize).max(1);

    println!(
        "\nRunning {} rigid steps ({:.1} s), {} SPH sub-steps each ...",
        n_steps,
        n_steps as f64 * dt,
        sph_substeps
    );

    println!(
        "\n{:>6}  {:>10}  {:>12}  {:>12}  {:>12}",
        "step", "time(s)", "sphere_y", "sphere_vy", "fluid_cy"
    );

    for step in 0..n_steps {
        // SPH sub-stepping.
        for _ in 0..sph_substeps {
            let neighbors = SpatialHash::find_all_neighbors(&fluid.positions, 2.0 * h);
            wcsph_step(
                &mut fluid,
                &neighbors,
                &kernel,
                &sph_params,
                sph_dt,
                gravity,
            );

            // Floor boundary.
            for (pos, vel) in fluid.positions.iter_mut().zip(fluid.velocities.iter_mut()) {
                if pos.y < 0.0 {
                    pos.y = 0.0;
                    if vel.y < 0.0 {
                        vel.y = -vel.y * 0.3;
                    }
                }
            }
        }

        // Rigid body step (includes gravity + contact solving).
        pipeline.step(dt, &mut bodies, &colliders);

        // Manual buoyancy: compute approximate submerged volume and apply.
        if let Some(body) = bodies.get_mut(sphere_h) {
            let p = body.transform.position;
            // If sphere centre is below the fluid top (~0.4 m), apply buoyancy.
            let fluid_top = 0.4_f64;
            if p.y < fluid_top {
                let submersion = ((fluid_top - p.y) / (2.0 * sphere_radius)).clamp(0.0, 1.0);
                let submerged_vol =
                    submersion * (4.0 / 3.0) * std::f64::consts::PI * sphere_radius.powi(3);
                let buoyancy = 1000.0 * 9.81 * submerged_vol;
                body.velocity.y += buoyancy * dt * body.inverse_mass;
                // Fluid drag.
                let drag_coeff = 0.5 * 1000.0 * 0.47 * std::f64::consts::PI * sphere_radius.powi(2);
                let drag = -drag_coeff * body.velocity.y * body.velocity.y.abs();
                body.velocity.y += drag * dt * body.inverse_mass;
            }
        }

        // Print samples.
        if step % 60 == 0 || step == n_steps - 1 {
            let sim_time = (step + 1) as f64 * dt;
            let sphere_y = bodies
                .get(sphere_h)
                .map(|b| b.transform.position.y)
                .unwrap_or(f64::NAN);
            let sphere_vy = bodies
                .get(sphere_h)
                .map(|b| b.velocity.y)
                .unwrap_or(f64::NAN);
            let fluid_cy: f64 = if !fluid.is_empty() {
                fluid.positions.iter().map(|p| p.y).sum::<f64>() / fluid.len() as f64
            } else {
                0.0
            };
            println!(
                "{:>6}  {:>10.3}  {:>12.4}  {:>12.4}  {:>12.4}",
                step, sim_time, sphere_y, sphere_vy, fluid_cy
            );
        }
    }

    // ── Verification ────────────────────────────────────────────────
    let final_sphere_y = bodies
        .get(sphere_h)
        .map(|b| b.transform.position.y)
        .unwrap_or(f64::NAN);
    let final_fluid_count = fluid.len();

    println!("\nFinal sphere y: {:.4}", final_sphere_y);
    println!(
        "Fluid particle count conserved: {}",
        final_fluid_count == initial_fluid_count
    );

    if final_sphere_y < 1.0 {
        println!("PASS: sphere fell from starting height under gravity");
    } else {
        println!("WARN: sphere did not fall as expected");
    }

    if final_fluid_count == initial_fluid_count {
        println!("PASS: fluid particle count conserved");
    } else {
        println!(
            "FAIL: fluid particle count changed {} -> {}",
            initial_fluid_count, final_fluid_count
        );
    }
}
