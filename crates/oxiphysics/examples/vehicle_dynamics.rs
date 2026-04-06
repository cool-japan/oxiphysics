// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Vehicle dynamics demo: Pacejka tire model on a flat road.
//!
//! Creates a 4-wheel raycast vehicle, applies throttle for 5 seconds,
//! then prints speed, tire slip, and lateral force at each step sample.

use oxiphysics_core::math::Vec3;
use oxiphysics_rigid::RigidBody;
use oxiphysics_vehicle::{
    LinearSuspension, PacejkaCoeffs, PacejkaTire, RaycastVehicle, TireModel, flat_ground_query,
};

fn main() {
    // Build a default 4-wheel vehicle.
    let mut vehicle = RaycastVehicle::new_default_4wheel();

    // Configure Pacejka tire with typical street-tire coefficients.
    let tire = PacejkaTire {
        lateral: PacejkaCoeffs {
            b: 10.0,
            c: 1.9,
            d: 1.0,
            e: 0.97,
        },
        longitudinal: PacejkaCoeffs {
            b: 12.0,
            c: 2.3,
            d: 1.0,
            e: 0.97,
        },
        combined_slip_weight: 0.5,
    };
    let suspension = LinearSuspension;

    // Chassis rigid body (1500 kg sedan).
    let mut chassis = RigidBody::new(1500.0);
    chassis.transform.position = Vec3::new(0.0, 0.5, 0.0);

    let gravity = Vec3::new(0.0, -9.81, 0.0);
    let ground = flat_ground_query(0.0);

    let dt = 1.0 / 120.0;
    let total_time = 5.0;
    let steps = (total_time / dt) as usize;
    let sample_interval = steps / 10; // Print 10 samples

    println!("=== vehicle_dynamics: Pacejka tire acceleration demo ===");
    println!(
        "{:>6}  {:>10}  {:>12}  {:>12}  {:>12}",
        "step", "time(s)", "speed(m/s)", "slip_ratio", "lat_force(N)"
    );

    for step in 0..steps {
        // Full throttle, no braking, no steering.
        vehicle.set_inputs(1.0, 0.0, 0.0);
        vehicle.update(dt, &mut chassis, &gravity, &*ground, &suspension, &tire);

        // Update speed from chassis velocity.
        vehicle.state.speed = chassis.velocity.z;

        if step % sample_interval == 0 || step == steps - 1 {
            let sim_time = (step + 1) as f64 * dt;

            // Compute sample tire forces at a representative slip condition.
            // Use wheel 2 (rear-left, driven wheel) state for display.
            let wheel_speed = vehicle.wheel_states[2].spin_velocity * vehicle.wheels[2].radius;
            let chassis_speed = vehicle.state.speed.abs().max(0.01);
            let slip_ratio = (wheel_speed - chassis_speed) / chassis_speed;

            // Query Pacejka for lateral force at a small slip angle (1 degree).
            let slip_angle = 1.0_f64.to_radians();
            let normal_force = 1500.0 * 9.81 / 4.0; // approximate per-wheel
            let friction = 1.0;
            let (lat_force, _long_force) =
                tire.compute_forces(slip_angle, slip_ratio, normal_force, friction);

            println!(
                "{:>6}  {:>10.3}  {:>12.4}  {:>12.6}  {:>12.2}",
                step, sim_time, chassis_speed, slip_ratio, lat_force
            );
        }
    }

    let final_speed = chassis.velocity.norm();
    println!("\nFinal chassis speed: {:.4} m/s", final_speed);
    println!(
        "Final position: ({:.2}, {:.2}, {:.2})",
        chassis.transform.position.x, chassis.transform.position.y, chassis.transform.position.z,
    );

    // Sanity check: vehicle should have moved forward.
    if chassis.transform.position.z.abs() > 0.01 || final_speed > 0.01 {
        println!("PASS: vehicle moved under throttle");
    } else {
        println!("INFO: vehicle barely moved (check drivetrain config)");
    }
}
