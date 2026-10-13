use oxiphysics_core::math::Vec3;
use oxiphysics_rigid::RigidBody;
use oxiphysics_vehicle::{LinearSuspension, LinearTire};

#[path = "zz_baseline/vehicle_baseline.rs"]
#[allow(dead_code, unused_imports)]
mod baseline;

macro_rules! run_impl {
    ($name:ident, $vehicle_ty:ty, $ground:expr) => {
        fn $name(dt: f64, seconds: f64, v0: f64, brake: f64, tag: &str) {
            let mut v = <$vehicle_ty>::new_default_4wheel();
            let mut chassis = RigidBody::new(1500.0);
            chassis.transform.position = Vec3::new(0.0, 0.52, 0.0);
            chassis.linear_damping = 0.0;
            chassis.angular_damping = 0.0;
            chassis.velocity = Vec3::new(0.0, 0.0, v0);
            for ws in v.wheel_states.iter_mut() {
                ws.spin_velocity = v0 / 0.3;
            }
            let gravity = Vec3::new(0.0, -9.81, 0.0);
            let ground = $ground;
            let suspension = LinearSuspension;
            let tire = LinearTire::default();
            v.drivetrain.gearbox.shift(0);
            v.set_inputs(0.0, brake, 0.0);
            let steps = (seconds / dt).round() as usize;
            let mut max_spin_last = 0.0f64;
            let mut min_vz_last = f64::MAX;
            let mut max_vz_last = f64::MIN;
            for i in 0..steps {
                chassis.force_accumulator = Vec3::zeros();
                chassis.torque_accumulator = Vec3::zeros();
                v.update(dt, &mut chassis, &gravity, &ground, &suspension, &tire);
                chassis.integrate_forces(dt, &gravity);
                chassis.integrate_velocity(dt);
                if (i as f64) * dt > seconds - 1.0 {
                    for ws in &v.wheel_states {
                        max_spin_last = max_spin_last.max(ws.spin_velocity.abs());
                    }
                    min_vz_last = min_vz_last.min(chassis.velocity.z);
                    max_vz_last = max_vz_last.max(chassis.velocity.z);
                }
            }
            println!(
                "SCRATCH {tag} dt={dt:.6} v0={v0} brake={brake} max|spin|={max_spin_last:.4}, vz range=[{min_vz_last:.5},{max_vz_last:.5}], z={:.4}, final vz={:.5}",
                chassis.transform.position.z, chassis.velocity.z
            );
        }
    };
}

run_impl!(
    run_new,
    oxiphysics_vehicle::RaycastVehicle,
    oxiphysics_vehicle::flat_ground_query(0.0)
);
run_impl!(
    run_old,
    baseline::RaycastVehicle,
    baseline::flat_ground_query(0.0)
);

#[test]
fn scratch_standstill() {
    for dt in [1.0 / 60.0, 1.0 / 120.0, 1.0 / 240.0, 1.0 / 1000.0] {
        for v0 in [0.0, 0.01, 0.3, 2.0, 10.0] {
            run_old(dt, 10.0, v0, 0.0, "OLD");
            run_new(dt, 10.0, v0, 0.0, "NEW");
        }
    }
}
