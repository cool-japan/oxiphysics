// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Tests for all engine submodules.

#![allow(missing_docs)]
#![cfg(test)]

use super::*;

// ===========================================================================
// Core WasmPhysicsEngine tests
// ===========================================================================

#[test]
fn test_engine_creation() {
    let engine = WasmPhysicsEngine::new(0.0, -9.81, 0.0);
    assert_eq!(engine.get_body_count(), 0);
    assert!((engine.gravity()[1] + 9.81).abs() < 1e-10);
}

#[test]
fn test_add_remove_body() {
    let mut engine = WasmPhysicsEngine::new(0.0, -9.81, 0.0);
    let h = engine.add_dynamic_body(1.0, 0.0, 0.0, 0.0);
    assert_eq!(engine.get_body_count(), 1);
    engine.remove_body(h).expect("should remove");
    assert_eq!(engine.get_body_count(), 0);
}

#[test]
fn test_remove_invalid_handle() {
    let mut engine = WasmPhysicsEngine::new(0.0, -9.81, 0.0);
    assert!(engine.remove_body(99).is_err());
}

#[test]
fn test_add_bodies_and_step() {
    let mut engine = WasmPhysicsEngine::new(0.0, -10.0, 0.0);
    let h0 = engine.add_dynamic_body(1.0, 0.0, 100.0, 0.0);
    let h1 = engine.add_static_body(0.0, 0.0, 0.0);
    assert_eq!(engine.get_body_count(), 2);

    engine.step(1.0);

    // Dynamic body should have moved (below 100)
    let pos = engine.get_position(h0);
    assert!(pos[1] < 100.0, "body should fall, got y={}", pos[1]);

    // Static body should not move
    let pos1 = engine.get_position(h1);
    assert!((pos1[1]).abs() < 1e-10);
}

#[test]
fn test_gravity_integration() {
    let mut engine = WasmPhysicsEngine::new(0.0, -10.0, 0.0);
    let h = engine.add_dynamic_body(1.0, 0.0, 50.0, 0.0);
    // Step with one substep of exactly 0.5s
    let mut cfg = engine.config.clone();
    cfg.fixed_dt = 0.5;
    engine.config = cfg;
    engine.step(0.5);

    let vel = engine.get_velocity(h);
    // v = g * dt = -10 * 0.5 = -5 (approximately, with damping)
    assert!(vel[1] < 0.0, "velocity should be negative");

    let pos = engine.get_position(h);
    // position should be below 50
    assert!(pos[1] < 50.0, "position should decrease");
}

#[test]
fn test_set_velocity() {
    let mut engine = WasmPhysicsEngine::new(0.0, 0.0, 0.0); // no gravity
    let h = engine.add_dynamic_body(1.0, 0.0, 0.0, 0.0);
    engine
        .set_velocity(h, 10.0, 0.0, 0.0)
        .expect("set velocity");
    engine.set_fixed_dt(1.0);
    engine.set_body_linear_damping(h, 0.0).expect("set damping");
    engine.step(1.0);

    let pos = engine.get_position(h);
    assert!(
        (pos[0] - 10.0).abs() < 1e-6,
        "expected x~10, got {}",
        pos[0]
    );
}

#[test]
fn test_apply_force() {
    let mut engine = WasmPhysicsEngine::new(0.0, 0.0, 0.0);
    let h = engine.add_dynamic_body(1.0, 0.0, 0.0, 0.0);
    engine.apply_force(h, 10.0, 0.0, 0.0).expect("apply force");
    engine.set_body_linear_damping(h, 0.0).expect("set damping");
    engine.set_fixed_dt(1.0);
    engine.step(1.0);
    let vel = engine.get_velocity(h);
    // v = F/m * dt = 10 * 1 = 10
    assert!(
        (vel[0] - 10.0).abs() < 1e-6,
        "expected vx~10, got {}",
        vel[0]
    );
}

#[test]
fn test_apply_impulse() {
    let mut engine = WasmPhysicsEngine::new(0.0, 0.0, 0.0);
    let h = engine.add_dynamic_body(2.0, 0.0, 0.0, 0.0);
    engine
        .apply_impulse(h, 10.0, 0.0, 0.0)
        .expect("apply impulse");
    let vel = engine.get_velocity(h);
    // dv = J/m = 10/2 = 5
    assert!(
        (vel[0] - 5.0).abs() < 1e-10,
        "expected vx=5, got {}",
        vel[0]
    );
}

#[test]
fn test_get_all_positions_flat() {
    let mut engine = WasmPhysicsEngine::new(0.0, 0.0, 0.0);
    engine.add_dynamic_body(1.0, 1.0, 2.0, 3.0);
    engine.add_dynamic_body(1.0, 4.0, 5.0, 6.0);

    let all = engine.get_all_positions();
    assert_eq!(all.len(), 6);
    assert!((all[0] - 1.0).abs() < 1e-10);
    assert!((all[1] - 2.0).abs() < 1e-10);
    assert!((all[2] - 3.0).abs() < 1e-10);
}

#[test]
fn test_colliders_sphere() {
    let mut engine = WasmPhysicsEngine::new(0.0, -9.81, 0.0);
    let body = engine.add_dynamic_body(1.0, 0.0, 10.0, 0.0);
    let c0 = engine.add_sphere_collider(body, 0.5);
    let c1 = engine.add_box_collider(body, 1.0, 1.0, 1.0);
    assert_eq!(c0, 0);
    assert_eq!(c1, 1);
}

#[test]
fn test_sphere_plane_contact() {
    let mut engine = WasmPhysicsEngine::new(0.0, -10.0, 0.0);
    // Body starts near the plane (y=1, sphere radius=1.5 -> penetrating plane at y=0)
    let dyn_body = engine.add_dynamic_body(1.0, 0.0, 1.0, 0.0);
    let static_body = engine.add_static_body(0.0, 0.0, 0.0);
    engine.add_sphere_collider(dyn_body, 1.5);
    engine.add_plane_collider(static_body, 0.0, 1.0, 0.0, 0.0);

    engine.step(1.0 / 60.0);
    // Should have at least one contact
    assert!(
        engine.get_contact_count() >= 1,
        "expected contact, got {}",
        engine.get_contact_count()
    );
}

#[test]
fn test_sphere_sphere_contact() {
    let mut engine = WasmPhysicsEngine::new(0.0, 0.0, 0.0);
    let b0 = engine.add_dynamic_body(1.0, 0.0, 0.0, 0.0);
    let b1 = engine.add_dynamic_body(1.0, 1.5, 0.0, 0.0); // distance = 1.5 < r0+r1=2.0
    engine.add_sphere_collider(b0, 1.0);
    engine.add_sphere_collider(b1, 1.0);
    let _ = engine.set_velocity(b0, 2.0, 0.0, 0.0);
    let _ = engine.set_velocity(b1, -2.0, 0.0, 0.0);
    engine.step(1.0 / 60.0);
    assert!(
        engine.get_contact_count() >= 1,
        "expected sphere-sphere contact"
    );
}

#[test]
fn test_reset() {
    let mut engine = WasmPhysicsEngine::new(0.0, -9.81, 0.0);
    engine.add_dynamic_body(1.0, 0.0, 0.0, 0.0);
    engine.step(1.0 / 60.0);
    engine.reset();
    assert_eq!(engine.get_body_count(), 0);
    assert!((engine.time()).abs() < 1e-15);
}

#[test]
fn test_set_gravity() {
    let mut engine = WasmPhysicsEngine::new(0.0, -9.81, 0.0);
    engine.set_gravity(0.0, -1.62, 0.0);
    assert!((engine.gravity()[1] + 1.62).abs() < 1e-10);
}

#[test]
fn test_get_body_state() {
    let mut engine = WasmPhysicsEngine::new(0.0, 0.0, 0.0);
    let h = engine.add_dynamic_body(2.0, 1.0, 2.0, 3.0);
    let state = engine.get_body_state(h).expect("should have state");
    assert_eq!(state.handle, h);
    assert!((state.position[0] - 1.0).abs() < 1e-10);
}

#[test]
fn test_handle_reuse_after_remove() {
    let mut engine = WasmPhysicsEngine::new(0.0, 0.0, 0.0);
    let h0 = engine.add_dynamic_body(1.0, 0.0, 0.0, 0.0);
    engine.remove_body(h0).expect("remove body");
    let h1 = engine.add_dynamic_body(1.0, 5.0, 5.0, 5.0);
    // The free list reuses the slot
    assert_eq!(h1, h0);
    let pos = engine.get_position(h1);
    assert!((pos[0] - 5.0).abs() < 1e-10);
}

#[test]
fn test_raycast_sphere_hit() {
    let mut engine = WasmPhysicsEngine::new(0.0, 0.0, 0.0);
    let b = engine.add_dynamic_body(1.0, 10.0, 0.0, 0.0);
    engine.add_sphere_collider(b, 1.0);
    let result = engine.raycast(0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 100.0);
    assert!(result.hit, "expected ray to hit sphere");
    assert!((result.distance - 9.0).abs() < 1e-6);
}

#[test]
fn test_raycast_no_hit() {
    let mut engine = WasmPhysicsEngine::new(0.0, 0.0, 0.0);
    let b = engine.add_dynamic_body(1.0, 0.0, 100.0, 0.0); // sphere above
    engine.add_sphere_collider(b, 1.0);
    let result = engine.raycast(0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 50.0);
    assert!(!result.hit);
}

#[test]
fn test_get_all_transforms() {
    let mut engine = WasmPhysicsEngine::new(0.0, 0.0, 0.0);
    engine.add_dynamic_body(1.0, 1.0, 2.0, 3.0);
    let t = engine.get_all_transforms();
    assert_eq!(t.len(), 7); // one body -> 7 floats
    assert!((t[0] - 1.0).abs() < 1e-10);
}

#[test]
fn test_invalid_handle_errors() {
    let mut engine = WasmPhysicsEngine::new(0.0, -9.81, 0.0);
    assert!(engine.set_velocity(99, 0.0, 0.0, 0.0).is_err());
    assert!(engine.apply_force(99, 0.0, 0.0, 0.0).is_err());
    assert!(engine.apply_impulse(99, 0.0, 0.0, 0.0).is_err());
}

// ===========================================================================
// Float64View, LBM, SPH, Vehicle, Constraint tests
// ===========================================================================

#[test]
fn test_float64view_basic() {
    let v = Float64View::from_vec(vec![1.0, 2.0, 3.0]);
    assert_eq!(v.len(), 3);
    assert!(!v.is_empty());
    assert_eq!(v.get(1), 2.0);
    assert!(v.get(99).is_nan());
}

#[test]
fn test_float64view_empty() {
    let v = Float64View::empty();
    assert!(v.is_empty());
    assert_eq!(v.len(), 0);
}

#[test]
fn test_float64view_into_vec() {
    let data = vec![4.0, 5.0, 6.0];
    let v = Float64View::from_vec(data.clone());
    assert_eq!(v.into_vec(), data);
}

#[test]
fn test_lbm_creation() {
    let cfg = EngineWasmLbmConfig {
        width: 4,
        height: 4,
        viscosity: 0.1,
    };
    let sim = WasmLbmSim::new(&cfg);
    assert_eq!(sim.width(), 4);
    assert_eq!(sim.height(), 4);
    assert_eq!(sim.step_count(), 0.0);
}

#[test]
fn test_lbm_initial_density() {
    let cfg = EngineWasmLbmConfig {
        width: 4,
        height: 4,
        viscosity: 0.1,
    };
    let sim = WasmLbmSim::new(&cfg);
    let rho = sim.density_at(0, 0);
    assert!((rho - 1.0).abs() < 1e-10, "rho = {}", rho);
}

#[test]
fn test_lbm_step_increments_count() {
    let cfg = EngineWasmLbmConfig {
        width: 4,
        height: 4,
        viscosity: 0.1,
    };
    let mut sim = WasmLbmSim::new(&cfg);
    sim.step();
    sim.step();
    assert_eq!(sim.step_count(), 2.0);
}

#[test]
fn test_lbm_velocity_field_size() {
    let cfg = EngineWasmLbmConfig {
        width: 8,
        height: 6,
        viscosity: 0.05,
    };
    let sim = WasmLbmSim::new(&cfg);
    let vf = sim.velocity_field();
    assert_eq!(vf.len(), 8 * 6 * 2);
}

#[test]
fn test_lbm_density_field_size() {
    let cfg = EngineWasmLbmConfig {
        width: 8,
        height: 6,
        viscosity: 0.05,
    };
    let sim = WasmLbmSim::new(&cfg);
    let df = sim.density_field();
    assert_eq!(df.len(), 8 * 6);
}

#[test]
fn test_sph_add_particles() {
    let mut sim = WasmSphSim::new();
    sim.add_particle(0.0, 0.0, 0.0);
    sim.add_particle(0.1, 0.0, 0.0);
    assert_eq!(sim.particle_count(), 2);
}

#[test]
fn test_sph_step_moves_particle() {
    let mut sim = WasmSphSim::new();
    sim.add_particle(0.0, 1.0, 0.0);
    let y0 = sim.get_position(0)[1];
    sim.step(0.001);
    let y1 = sim.get_position(0)[1];
    assert!(y1 < y0, "particle should fall: {} -> {}", y0, y1);
}

#[test]
fn test_sph_all_positions_length() {
    let mut sim = WasmSphSim::new();
    sim.add_particle(0.0, 0.0, 0.0);
    sim.add_particle(1.0, 0.0, 0.0);
    assert_eq!(sim.all_positions().len(), 6);
}

#[test]
fn test_sph_time_accumulation() {
    let mut sim = WasmSphSim::new();
    sim.add_particle(0.0, 0.0, 0.0);
    sim.step(0.01);
    sim.step(0.01);
    assert!((sim.time() - 0.02).abs() < 1e-12);
}

#[test]
fn test_vehicle_creation() {
    let v = WasmVehicleSim::new();
    let state = v.get_state();
    assert_eq!(state.gear, 1);
    assert!(state.position[0].abs() < 1e-10);
}

#[test]
fn test_vehicle_throttle_accelerates() {
    let mut v = WasmVehicleSim::new();
    v.set_throttle(1.0);
    for _ in 0..100 {
        v.step(0.05);
    }
    let state = v.get_state();
    // After 5 seconds at full throttle the car should be moving
    assert!(
        state.speed > 0.1,
        "expected speed > 0.1, got {}",
        state.speed
    );
}

#[test]
fn test_vehicle_steering() {
    let mut v = WasmVehicleSim::new();
    v.set_throttle(1.0);
    v.set_steering(0.5);
    for _ in 0..60 {
        v.step(1.0 / 60.0);
    }
    let state = v.get_state();
    // Heading should have changed with steering
    assert!(state.heading.abs() > 0.01, "heading = {}", state.heading);
}

#[test]
fn test_vehicle_brake() {
    let mut v = WasmVehicleSim::new();
    v.set_throttle(1.0);
    for _ in 0..100 {
        v.step(0.05);
    }
    let speed_before = v.get_state().speed;
    v.set_throttle(0.0);
    v.set_brake(1.0);
    for _ in 0..60 {
        v.step(0.05);
    }
    let speed_after = v.get_state().speed;
    assert!(
        speed_after < speed_before,
        "expected braking: {} -> {}",
        speed_before,
        speed_after
    );
}

#[test]
fn test_vehicle_time_advances() {
    let mut v = WasmVehicleSim::new();
    v.step(0.1);
    v.step(0.1);
    assert!((v.time() - 0.2).abs() < 1e-12);
}

#[test]
fn test_constraint_solver_creation() {
    let solver = EngineConstraintSolver::new();
    assert_eq!(solver.contact_count(), 0);
    assert_eq!(solver.constraint_count(), 0);
}

#[test]
fn test_constraint_solver_add_contact() {
    let mut solver = EngineConstraintSolver::new();
    solver.add_contact(0, 1, 0.0, 1.0, 0.0, 0.05, 0.3);
    assert_eq!(solver.contact_count(), 1);
    solver.clear_contacts();
    assert_eq!(solver.contact_count(), 0);
}

#[test]
fn test_constraint_solver_add_constraint() {
    let mut solver = EngineConstraintSolver::new();
    solver.add_constraint(0, 1, 0.0, 1.0, 0.0, 0.0);
    assert_eq!(solver.constraint_count(), 1);
    solver.clear_constraints();
    assert_eq!(solver.constraint_count(), 0);
}

#[test]
fn test_constraint_solver_solve() {
    let mut engine = WasmPhysicsEngine::new(0.0, 0.0, 0.0);
    let b0 = engine.add_dynamic_body(1.0, 0.0, 0.0, 0.0);
    let b1 = engine.add_dynamic_body(1.0, 0.1, 0.0, 0.0);

    let mut solver = EngineConstraintSolver::new();
    // Contact: bodies overlapping at depth 0.05
    solver.add_contact(b0, b1, -1.0, 0.0, 0.0, 0.05, 0.0);
    let rows = solver.solve(&mut engine, 1.0 / 60.0);
    assert_eq!(rows, 1, "should have processed 1 constraint row");
}

#[test]
fn test_constraint_solver_iterations() {
    let mut solver = EngineConstraintSolver::new();
    solver.set_iterations(20);
    assert_eq!(solver.iterations, 20);
    solver.set_iterations(0); // clamped to 1
    assert_eq!(solver.iterations, 1);
}

// ===========================================================================
// WasmVec3, WasmTransform, error helpers, bindings tests
// ===========================================================================

#[test]
fn test_wasm_vec3_new() {
    let v = WasmVec3::new(1.0, 2.0, 3.0);
    assert!((v.x() - 1.0).abs() < 1e-12);
    assert!((v.y() - 2.0).abs() < 1e-12);
    assert!((v.z() - 3.0).abs() < 1e-12);
}

#[test]
fn test_wasm_vec3_zero() {
    let v = WasmVec3::zero();
    assert!((v.length()).abs() < 1e-12);
}

#[test]
fn test_wasm_vec3_length() {
    let v = WasmVec3::new(3.0, 4.0, 0.0);
    assert!((v.length() - 5.0).abs() < 1e-10);
}

#[test]
fn test_wasm_vec3_normalized() {
    let v = WasmVec3::new(0.0, 5.0, 0.0);
    let n = v.normalized();
    assert!((n.y() - 1.0).abs() < 1e-10);
    assert!(n.x().abs() < 1e-10);
}

#[test]
fn test_wasm_vec3_dot() {
    let a = WasmVec3::new(1.0, 0.0, 0.0);
    let b = WasmVec3::new(0.0, 1.0, 0.0);
    assert!(a.dot(&b).abs() < 1e-12);
    assert!((a.dot(&a) - 1.0).abs() < 1e-12);
}

#[test]
fn test_wasm_vec3_add_sub() {
    let a = WasmVec3::new(1.0, 2.0, 3.0);
    let b = WasmVec3::new(4.0, 5.0, 6.0);
    let s = a.add(&b);
    assert!((s.x() - 5.0).abs() < 1e-12);
    let d = b.sub(&a);
    assert!((d.x() - 3.0).abs() < 1e-12);
}

#[test]
fn test_wasm_vec3_scale() {
    let v = WasmVec3::new(2.0, 3.0, 4.0);
    let s = v.scale(2.0);
    assert!((s.x() - 4.0).abs() < 1e-12);
    assert!((s.y() - 6.0).abs() < 1e-12);
}

#[test]
fn test_wasm_vec3_distance_to() {
    let a = WasmVec3::new(0.0, 0.0, 0.0);
    let b = WasmVec3::new(3.0, 4.0, 0.0);
    assert!((a.distance_to(&b) - 5.0).abs() < 1e-10);
}

#[test]
fn test_wasm_vec3_to_json() {
    let v = WasmVec3::new(1.0, 2.0, 3.0);
    let j = v.to_json();
    assert!(j.contains("\"x\":1"));
    assert!(j.contains("\"y\":2"));
    assert!(j.contains("\"z\":3"));
}

#[test]
fn test_wasm_vec3_array_roundtrip() {
    let arr = [7.0, 8.0, 9.0];
    let v = WasmVec3::from_array(arr);
    assert_eq!(v.to_array(), arr);
}

#[test]
fn test_wasm_vec3_setters() {
    let mut v = WasmVec3::new(0.0, 0.0, 0.0);
    v.set_x(5.0);
    v.set_y(6.0);
    v.set_z(7.0);
    assert!((v.x() - 5.0).abs() < 1e-12);
    assert!((v.y() - 6.0).abs() < 1e-12);
    assert!((v.z() - 7.0).abs() < 1e-12);
}

#[test]
fn test_wasm_transform_identity() {
    let t = WasmTransform::identity();
    let flat = t.to_flat();
    assert_eq!(flat.len(), 7);
    assert!(flat[0].abs() < 1e-12);
    assert!(flat[1].abs() < 1e-12);
    assert!(flat[2].abs() < 1e-12);
    assert!((flat[6] - 1.0).abs() < 1e-12);
}

#[test]
fn test_wasm_transform_from_position() {
    let t = WasmTransform::from_position(1.0, 2.0, 3.0);
    let p = t.position();
    assert!((p.x() - 1.0).abs() < 1e-12);
    assert!((p.y() - 2.0).abs() < 1e-12);
    assert!((p.z() - 3.0).abs() < 1e-12);
}

#[test]
fn test_wasm_transform_from_flat_roundtrip() {
    let t = WasmTransform::from_position(4.0, 5.0, 6.0);
    let flat = t.to_flat();
    let t2 = WasmTransform::from_flat(&flat).expect("should parse");
    let p = t2.position();
    assert!((p.x() - 4.0).abs() < 1e-12);
    assert!((p.y() - 5.0).abs() < 1e-12);
}

#[test]
fn test_wasm_transform_from_flat_too_short() {
    let short = vec![1.0, 2.0];
    assert!(WasmTransform::from_flat(&short).is_none());
}

#[test]
fn test_wasm_transform_matrix4_identity() {
    let t = WasmTransform::identity();
    let m = t.to_matrix4();
    assert!((m[0] - 1.0).abs() < 1e-10);
    assert!((m[5] - 1.0).abs() < 1e-10);
    assert!((m[10] - 1.0).abs() < 1e-10);
    assert!((m[15] - 1.0).abs() < 1e-10);
}

#[test]
fn test_wasm_transform_to_json() {
    let t = WasmTransform::identity();
    let j = t.to_json();
    assert!(j.contains("\"position\""));
    assert!(j.contains("\"rotation\""));
    assert!(j.contains("\"w\":1"));
}

#[test]
fn test_wasm_transform_position_array() {
    let t = WasmTransform::from_position(10.0, 20.0, 30.0);
    let pa = t.position_array();
    assert!((pa[0] - 10.0).abs() < 1e-12);
    assert!((pa[1] - 20.0).abs() < 1e-12);
    assert!((pa[2] - 30.0).abs() < 1e-12);
}

#[test]
fn test_error_to_js_string() {
    use crate::error::Error;
    let s = error_to_js_string(&format!("{}", Error::InvalidHandle(99)));
    assert!(s.contains("InvalidHandle") || s.contains("error"));
}

#[test]
fn test_result_to_js_ok() {
    let r: crate::error::Result<u32> = Ok(42);
    let js = result_to_js(r);
    assert_eq!(js, std::result::Result::Ok(42));
}

// Gated to wasm32 only: under host (non-wasm32) builds, wasm-bindgen 0.2.120's
// `JsValue::from_str` invokes the `__wbindgen_describe` machinery, which aborts
// the process when no wasm runtime is present. The Ok branch is exercised by
// `test_result_to_js_ok`; for the Err branch we only need the wasm32 path.
#[cfg(target_arch = "wasm32")]
#[test]
fn test_result_to_js_err() {
    let r: crate::error::Result<u32> = Err(Error::InvalidHandle(7));
    let js = result_to_js(r);
    assert!(js.is_err());
}

#[test]
fn test_free_f64_buffer() {
    let buf = vec![1.0f64, 2.0, 3.0];
    free_f64_buffer(buf); // should not panic
}

#[test]
fn test_free_u32_buffer() {
    let buf = vec![1u32, 2, 3];
    free_u32_buffer(buf);
}

#[test]
fn test_free_u8_buffer() {
    let buf = vec![1u8, 2, 3];
    free_u8_buffer(buf);
}

#[test]
fn test_wasm_page_size() {
    assert_eq!(WASM_PAGE_SIZE, 65536);
}

#[test]
fn test_ts_definitions_content() {
    assert!(TS_DEFINITIONS.contains("WasmPhysicsEngine"));
    assert!(TS_DEFINITIONS.contains("bodyCount()"));
    assert!(TS_DEFINITIONS.contains("addRigidBody"));
    assert!(TS_DEFINITIONS.contains("setVelocity"));
    assert!(TS_DEFINITIONS.contains("getAllPositions"));
    assert!(TS_DEFINITIONS.contains("setGravity"));
    assert!(TS_DEFINITIONS.contains("getContacts"));
    assert!(TS_DEFINITIONS.contains("addColliderSphere"));
    assert!(TS_DEFINITIONS.contains("addColliderBox"));
    assert!(TS_DEFINITIONS.contains("webgpuComputePlaceholder"));
}

#[test]
fn test_ts_vec3_def() {
    assert!(TS_VEC3_DEF.contains("Vec3"));
    assert!(TS_VEC3_DEF.contains("number"));
}

#[test]
fn test_ts_transform_def() {
    assert!(TS_TRANSFORM_DEF.contains("Transform"));
    assert!(TS_TRANSFORM_DEF.contains("Float64Array"));
}

#[test]
fn test_ts_contact_def() {
    assert!(TS_CONTACT_DEF.contains("ContactInfo"));
    assert!(TS_CONTACT_DEF.contains("impulse"));
}

#[test]
fn test_wasm_bindings_new() {
    let api = WasmBindings::new(0.0, -9.81, 0.0);
    assert_eq!(api.body_count(), 0);
}

#[test]
fn test_wasm_bindings_add_rigid_body() {
    let mut api = WasmBindings::new(0.0, -9.81, 0.0);
    let id = api.add_rigid_body(0.0, 10.0, 0.0, 1.0);
    assert_eq!(api.body_count(), 1);
    let pos = api.get_position(id);
    assert_eq!(pos.len(), 3);
    assert!((pos[1] - 10.0).abs() < 1e-12);
}

#[test]
fn test_wasm_bindings_add_static_body() {
    let mut api = WasmBindings::new(0.0, -9.81, 0.0);
    let _id = api.add_rigid_body(0.0, 0.0, 0.0, 0.0); // mass=0 -> static
    assert_eq!(api.body_count(), 1);
}

#[test]
fn test_wasm_bindings_step_moves_body() {
    let mut api = WasmBindings::new(0.0, -9.81, 0.0);
    let id = api.add_rigid_body(0.0, 10.0, 0.0, 1.0);
    api.step(1.0 / 60.0);
    let pos = api.get_position(id);
    assert!(pos[1] < 10.0, "body should have fallen");
}

#[test]
fn test_wasm_bindings_set_velocity() {
    let mut api = WasmBindings::new(0.0, 0.0, 0.0);
    let id = api.add_rigid_body(0.0, 0.0, 0.0, 1.0);
    api.set_velocity(id, 10.0, 0.0, 0.0);
    api.step(1.0);
    let pos = api.get_position(id);
    assert!(
        pos[0] > 0.5,
        "body should have moved along X, got {}",
        pos[0]
    );
}

#[test]
fn test_wasm_bindings_get_all_positions() {
    let mut api = WasmBindings::new(0.0, 0.0, 0.0);
    api.add_rigid_body(1.0, 2.0, 3.0, 1.0);
    api.add_rigid_body(4.0, 5.0, 6.0, 1.0);
    let all = api.get_all_positions();
    assert_eq!(all.len(), 6);
    assert!((all[0] - 1.0).abs() < 1e-12);
    assert!((all[3] - 4.0).abs() < 1e-12);
}

#[test]
fn test_wasm_bindings_set_gravity() {
    let mut api = WasmBindings::new(0.0, -9.81, 0.0);
    api.set_gravity(0.0, -1.62, 0.0);
    let g = api.gravity();
    assert!((g[1] + 1.62).abs() < 1e-10);
}

#[test]
fn test_wasm_bindings_add_colliders() {
    let mut api = WasmBindings::new(0.0, -9.81, 0.0);
    let id = api.add_rigid_body(0.0, 5.0, 0.0, 1.0);
    let c0 = api.add_collider_sphere(id, 0.5);
    let c1 = api.add_collider_box(id, 1.0, 1.0, 1.0);
    assert_eq!(c0, 0);
    assert_eq!(c1, 1);
}

#[test]
fn test_wasm_bindings_contacts_flat() {
    let mut api = WasmBindings::new(0.0, 0.0, 0.0);
    let b0 = api.add_rigid_body(0.0, 0.0, 0.0, 1.0);
    let b1 = api.add_rigid_body(1.5, 0.0, 0.0, 1.0);
    api.add_collider_sphere(b0, 1.0);
    api.add_collider_sphere(b1, 1.0);
    api.set_velocity(b0, 2.0, 0.0, 0.0);
    api.set_velocity(b1, -2.0, 0.0, 0.0);
    api.step(1.0 / 60.0);
    // Flat output length should be a multiple of 12
    let flat = api.get_contacts_flat();
    assert_eq!(flat.len() % 12, 0);
    if !flat.is_empty() {
        assert!((flat[0] - b0 as f64).abs() < 1e-12 || (flat[0] - b1 as f64).abs() < 1e-12);
    }
}

#[test]
fn test_wasm_bindings_contacts_json() {
    let mut api = WasmBindings::new(0.0, 0.0, 0.0);
    let b0 = api.add_rigid_body(0.0, 0.0, 0.0, 1.0);
    let b1 = api.add_rigid_body(1.5, 0.0, 0.0, 1.0);
    api.add_collider_sphere(b0, 1.0);
    api.add_collider_sphere(b1, 1.0);
    api.step(1.0 / 60.0);
    let json_list = api.get_contacts_json();
    for j in &json_list {
        assert!(j.contains("bodyA"));
        assert!(j.contains("normal"));
    }
}

#[test]
fn test_wasm_bindings_webgpu_placeholder() {
    let mut api = WasmBindings::new(0.0, -9.81, 0.0);
    api.add_rigid_body(0.0, 1.0, 0.0, 1.0);
    api.step(1.0 / 60.0);
    let result = api.webgpu_compute_placeholder(64);
    assert_eq!(result.len(), 64);
    assert!((result[0] - 1.0).abs() < 1e-10);
    assert!(result[1] > 0.0); // simulation time > 0
}

#[test]
fn test_wasm_bindings_webgpu_placeholder_min_size() {
    let api = WasmBindings::new(0.0, -9.81, 0.0);
    let result = api.webgpu_compute_placeholder(0);
    assert_eq!(result.len(), 2); // clamped to min 2
}

#[test]
fn test_wasm_bindings_get_transform() {
    let mut api = WasmBindings::new(0.0, 0.0, 0.0);
    let id = api.add_rigid_body(1.0, 2.0, 3.0, 1.0);
    let t = api.get_transform(id);
    let p = t.position();
    assert!((p.x() - 1.0).abs() < 1e-12);
    assert!((p.y() - 2.0).abs() < 1e-12);
    assert!((p.z() - 3.0).abs() < 1e-12);
}

#[test]
fn test_wasm_bindings_get_all_transforms_flat() {
    let mut api = WasmBindings::new(0.0, 0.0, 0.0);
    api.add_rigid_body(1.0, 2.0, 3.0, 1.0);
    let t = api.get_all_transforms_flat();
    assert_eq!(t.len(), 7);
    assert!((t[0] - 1.0).abs() < 1e-12);
}

#[test]
fn test_wasm_bindings_state_to_json() {
    let mut api = WasmBindings::new(0.0, -9.81, 0.0);
    api.add_rigid_body(0.0, 5.0, 0.0, 1.0);
    api.step(0.016);
    let json = api.state_to_json();
    assert!(json.contains("\"time\""));
    assert!(json.contains("\"gravity\""));
    assert!(json.contains("\"bodies\""));
    assert!(json.contains("\"pos\""));
}

#[test]
fn test_wasm_bindings_reset() {
    let mut api = WasmBindings::new(0.0, -9.81, 0.0);
    api.add_rigid_body(0.0, 0.0, 0.0, 1.0);
    api.reset();
    assert_eq!(api.body_count(), 0);
}

#[test]
fn test_wasm_bindings_apply_force() {
    let mut api = WasmBindings::new(0.0, 0.0, 0.0);
    let id = api.add_rigid_body(0.0, 0.0, 0.0, 1.0);
    api.apply_force(id, 100.0, 0.0, 0.0);
    api.step(0.1);
    let pos = api.get_position(id);
    assert!(pos[0] > 0.0, "body should move along X from applied force");
}

#[test]
fn test_wasm_bindings_apply_impulse() {
    let mut api = WasmBindings::new(0.0, 0.0, 0.0);
    let id = api.add_rigid_body(0.0, 0.0, 0.0, 2.0);
    api.apply_impulse(id, 4.0, 0.0, 0.0);
    api.step(1.0);
    let pos = api.get_position(id);
    assert!(pos[0] > 0.1, "body should move along X, got {}", pos[0]);
}

#[test]
fn test_contact_list_empty() {
    let list = WasmContactList::from_contacts(&[]);
    assert!(list.is_empty());
    assert_eq!(list.len(), 0);
    assert_eq!(list.to_json(), "[]");
}

#[test]
fn test_contact_list_from_contacts() {
    let mut c = crate::types::ContactResult::new(1, 2);
    c.depth = 0.05;
    c.impulse = 1.5;
    c.normal = [0.0, 1.0, 0.0];
    let list = WasmContactList::from_contacts(&[c]);
    assert_eq!(list.len(), 1);
    let entry = list.get(0).expect("should have entry");
    assert_eq!(entry.body_a, 1);
    assert_eq!(entry.body_b, 2);
    assert!((entry.depth - 0.05).abs() < 1e-12);
}

#[test]
fn test_contact_list_to_json() {
    let mut c = crate::types::ContactResult::new(3, 4);
    c.depth = 0.1;
    c.impulse = 2.0;
    c.normal = [1.0, 0.0, 0.0];
    let list = WasmContactList::from_contacts(&[c]);
    let json = list.to_json();
    assert!(json.starts_with('['));
    assert!(json.ends_with(']'));
    assert!(json.contains("\"bodyA\":3"));
    assert!(json.contains("\"bodyB\":4"));
}

#[test]
fn test_contact_list_to_flat() {
    let mut c = crate::types::ContactResult::new(0, 1);
    c.depth = 0.02;
    c.impulse = 0.5;
    c.is_new = true;
    let list = WasmContactList::from_contacts(&[c]);
    let flat = list.to_flat();
    assert_eq!(flat.len(), 8);
    assert!((flat[0] - 0.0).abs() < 1e-12); // body_a = 0
    assert!((flat[1] - 1.0).abs() < 1e-12); // body_b = 1
    assert!((flat[7] - 1.0).abs() < 1e-12); // is_new = true
}

#[test]
fn test_contact_info_entry_to_json() {
    let entry = ContactInfoEntry {
        body_a: 5,
        body_b: 6,
        normal: [0.0, 1.0, 0.0],
        depth: 0.03,
        impulse: 1.1,
        is_new: false,
    };
    let j = entry.to_json();
    assert!(j.contains("\"bodyA\":5"));
    assert!(j.contains("\"isNew\":false"));
}

#[test]
fn test_wasm_bindings_get_contacts_structured() {
    let mut api = WasmBindings::new(0.0, 0.0, 0.0);
    let b0 = api.add_rigid_body(0.0, 0.0, 0.0, 1.0);
    let b1 = api.add_rigid_body(1.5, 0.0, 0.0, 1.0);
    api.add_collider_sphere(b0, 1.0);
    api.add_collider_sphere(b1, 1.0);
    api.step(1.0 / 60.0);
    let clist = api.get_contacts_structured();
    for i in 0..clist.len() as usize {
        let entry = clist.get(i).expect("should have entry");
        assert!(entry.depth >= 0.0);
    }
}

// ===========================================================================
// WasmEngine tests (performance, debug draw, ray-cast)
// ===========================================================================

#[test]
fn test_get_performance_metrics_empty() {
    let engine = WasmEngine::new(0.0, -9.81, 0.0);
    let m = engine.get_performance_metrics();
    assert_eq!(m.body_count, 0);
    assert_eq!(m.contact_count, 0);
    assert!(m.fps > 0.0, "FPS must be positive");
}

#[test]
fn test_get_performance_metrics_body_count() {
    let mut engine = WasmEngine::new(0.0, -9.81, 0.0);
    engine.add_dynamic_body(1.0, 0.0, 0.0, 0.0);
    engine.add_dynamic_body(2.0, 5.0, 0.0, 0.0);
    let m = engine.get_performance_metrics();
    assert_eq!(m.body_count, 2);
}

#[test]
fn test_get_performance_metrics_sim_time_advances() {
    let mut engine = WasmEngine::new(0.0, -9.81, 0.0);
    engine.add_dynamic_body(1.0, 0.0, 10.0, 0.0);
    engine.step(0.1);
    let m = engine.get_performance_metrics();
    assert!(m.sim_time > 0.0, "sim_time should be > 0 after step");
}

#[test]
fn test_get_performance_metrics_fps_from_step_time() {
    let mut engine = WasmEngine::new(0.0, -9.81, 0.0);
    engine.step(0.016); // ~60 fps
    let m = engine.get_performance_metrics();
    assert!((m.fps - 62.5).abs() < 1.0, "fps ~ 62.5, got {}", m.fps);
}

#[test]
fn test_set_debug_draw_all_enabled() {
    let mut engine = WasmEngine::new(0.0, -9.81, 0.0);
    engine.set_debug_draw(DebugDrawFlags::all_enabled());
    let flags = engine.debug_draw_flags();
    assert!(flags.show_aabbs);
    assert!(flags.show_velocities);
    assert!(flags.show_contacts);
    assert!(flags.any_enabled());
}

#[test]
fn test_set_debug_draw_disabled() {
    let mut engine = WasmEngine::new(0.0, -9.81, 0.0);
    engine.set_debug_draw(DebugDrawFlags::disabled());
    assert!(!engine.debug_draw_flags().any_enabled());
}

#[test]
fn test_toggle_debug_flag() {
    let mut engine = WasmEngine::new(0.0, -9.81, 0.0);
    engine.toggle_debug_flag("aabbs", true);
    assert!(engine.debug_draw_flags().show_aabbs);
    engine.toggle_debug_flag("aabbs", false);
    assert!(!engine.debug_draw_flags().show_aabbs);
}

#[test]
fn test_toggle_unknown_flag_ignored() {
    let mut engine = WasmEngine::new(0.0, -9.81, 0.0);
    engine.toggle_debug_flag("nonexistent_flag", true);
    assert!(!engine.debug_draw_flags().any_enabled());
}

#[test]
fn test_compute_ray_cast_miss() {
    let engine = WasmEngine::new(0.0, -9.81, 0.0);
    let result = engine.compute_ray_cast(0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1000.0);
    assert!(!result.hit, "Ray in empty world should miss");
}

#[test]
fn test_compute_ray_cast_hit_sphere() {
    let mut engine = WasmEngine::new(0.0, 0.0, 0.0);
    let body = engine.add_dynamic_body(1.0, 0.0, 0.0, 10.0);
    engine.add_sphere_collider(body, 1.0);
    let result = engine.compute_ray_cast(0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1000.0);
    assert!(result.hit, "Ray should hit sphere at z=10");
    assert!(
        (result.distance - 9.0).abs() < 0.1,
        "dist ~ 9.0, got {}",
        result.distance
    );
}

#[test]
fn test_compute_ray_cast_flat_length() {
    let engine = WasmEngine::new(0.0, -9.81, 0.0);
    let flat = engine.compute_ray_cast_flat(0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 100.0);
    assert_eq!(flat.len(), 9, "flat ray cast result should have 9 elements");
}

#[test]
fn test_compute_ray_cast_flat_miss_encoding() {
    let engine = WasmEngine::new(0.0, -9.81, 0.0);
    let flat = engine.compute_ray_cast_flat(0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 100.0);
    assert!((flat[0] - 0.0).abs() < 1e-12, "flat[0] = 0 for miss");
}

// ===========================================================================
// Extended simulation control, queries, serialisation tests
// ===========================================================================

#[test]
fn test_step_n_advances_time() {
    let mut e = WasmEngine::new(0.0, -9.81, 0.0);
    e.step_n(5, 0.016);
    assert!(e.elapsed_time() > 0.0, "time should advance after step_n");
}

#[test]
fn test_step_n_moves_body() {
    let mut e = WasmEngine::new(0.0, -9.81, 0.0);
    let h = e.add_dynamic_body(1.0, 0.0, 10.0, 0.0);
    e.step_n(10, 0.016);
    let pos = e.get_body_position(h);
    assert!(pos[1] < 10.0, "body should fall under gravity");
}

#[test]
fn test_add_dynamic_sphere_returns_handle() {
    let mut e = WasmEngine::new(0.0, -9.81, 0.0);
    let h = e.add_dynamic_sphere(1.0, 0.0, 5.0, 0.0, 0.5);
    assert!(e.get_body_snapshot(h).is_some());
}

#[test]
fn test_add_static_plane_returns_valid_handle() {
    let mut e = WasmEngine::new(0.0, -9.81, 0.0);
    let h = e.add_static_plane(0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0);
    assert!(
        e.body_count() >= 1,
        "static plane body must be registered; handle={h}"
    );
}

#[test]
fn test_remove_body_returns_true() {
    let mut e = WasmEngine::new(0.0, -9.81, 0.0);
    let h = e.add_dynamic_body(1.0, 0.0, 0.0, 0.0);
    assert!(
        e.remove_body(h),
        "remove_body should return true for existing handle"
    );
}

#[test]
fn test_remove_body_invalid_returns_false() {
    let mut e = WasmEngine::new(0.0, -9.81, 0.0);
    assert!(
        !e.remove_body(9999),
        "remove_body should return false for unknown handle"
    );
}

#[test]
fn test_set_body_position_changes_position() {
    let mut e = WasmEngine::new(0.0, 0.0, 0.0);
    let h = e.add_dynamic_body(1.0, 0.0, 0.0, 0.0);
    let ok = e.set_body_position(h, 5.0, 5.0, 5.0);
    assert!(ok, "set_body_position should succeed");
    let p = e.get_body_position(h);
    assert!((p[0] - 5.0).abs() < 1e-10, "x should be 5");
    assert!((p[1] - 5.0).abs() < 1e-10, "y should be 5");
    assert!((p[2] - 5.0).abs() < 1e-10, "z should be 5");
}

#[test]
fn test_set_body_velocity_changes_velocity() {
    let mut e = WasmEngine::new(0.0, 0.0, 0.0);
    let h = e.add_dynamic_body(1.0, 0.0, 0.0, 0.0);
    e.set_body_velocity(h, 3.0, 0.0, 0.0);
    let v = e.get_body_velocity(h);
    assert!((v[0] - 3.0).abs() < 1e-10, "vx should be 3");
}

#[test]
fn test_apply_body_force_creates_acceleration() {
    let mut e = WasmEngine::new(0.0, 0.0, 0.0);
    let h = e.add_dynamic_body(1.0, 0.0, 0.0, 0.0);
    e.apply_body_force(h, 100.0, 0.0, 0.0);
    e.step(0.016);
    let v = e.get_body_velocity(h);
    assert!(v[0] > 0.0, "force should produce positive x velocity");
}

#[test]
fn test_apply_body_impulse_creates_velocity() {
    let mut e = WasmEngine::new(0.0, 0.0, 0.0);
    let h = e.add_dynamic_body(1.0, 0.0, 0.0, 0.0);
    e.apply_body_impulse(h, 10.0, 0.0, 0.0);
    let v = e.get_body_velocity(h);
    assert!(v[0] > 0.0, "impulse should create velocity");
}

#[test]
fn test_set_linear_damping_succeeds() {
    let mut e = WasmEngine::new(0.0, 0.0, 0.0);
    let h = e.add_dynamic_body(1.0, 0.0, 0.0, 0.0);
    assert!(e.set_linear_damping(h, 0.1));
}

#[test]
fn test_set_gravity_changes_gravity() {
    let mut e = WasmEngine::new(0.0, -9.81, 0.0);
    e.set_gravity(0.0, -1.0, 0.0);
    let g = e.get_gravity();
    assert!((g[1] + 1.0).abs() < 1e-10, "gravity y should be -1.0");
}

#[test]
fn test_reset_with_new_gravity_empties_world() {
    let mut e = WasmEngine::new(0.0, -9.81, 0.0);
    e.add_dynamic_body(1.0, 0.0, 10.0, 0.0);
    e.reset_with_new_gravity(0.0, -1.62, 0.0);
    assert_eq!(e.body_count(), 0, "reset should clear all bodies");
}

#[test]
fn test_body_snapshot_valid_handle() {
    let mut e = WasmEngine::new(0.0, -9.81, 0.0);
    let h = e.add_dynamic_body(2.0, 1.0, 2.0, 3.0);
    let snap = e.get_body_snapshot(h).expect("snapshot should exist");
    assert!((snap.position[0] - 1.0).abs() < 1e-8, "snapshot x");
    assert!((snap.position[1] - 2.0).abs() < 1e-8, "snapshot y");
    assert!((snap.position[2] - 3.0).abs() < 1e-8, "snapshot z");
}

#[test]
fn test_body_snapshot_invalid_handle_returns_none() {
    let e = WasmEngine::new(0.0, -9.81, 0.0);
    assert!(e.get_body_snapshot(9999).is_none());
}

#[test]
fn test_get_body_speed_at_rest() {
    let mut e = WasmEngine::new(0.0, 0.0, 0.0);
    let h = e.add_dynamic_body(1.0, 0.0, 0.0, 0.0);
    assert!(
        (e.get_body_speed(h)).abs() < 1e-10,
        "speed at rest should be 0"
    );
}

#[test]
fn test_query_aabb_finds_body_inside() {
    let mut e = WasmEngine::new(0.0, 0.0, 0.0);
    let h = e.add_dynamic_body(1.0, 5.0, 5.0, 5.0);
    let q = e.query_aabb([0.0, 0.0, 0.0, 10.0, 10.0, 10.0]);
    assert!(
        q.handles.contains(&h),
        "body at (5,5,5) should be inside AABB"
    );
}

#[test]
fn test_query_aabb_excludes_body_outside() {
    let mut e = WasmEngine::new(0.0, 0.0, 0.0);
    let h = e.add_dynamic_body(1.0, 100.0, 100.0, 100.0);
    let q = e.query_aabb([0.0, 0.0, 0.0, 10.0, 10.0, 10.0]);
    assert!(
        !q.handles.contains(&h),
        "body at (100,100,100) should be outside AABB"
    );
}

#[test]
fn test_query_sphere_overlap_finds_body_inside() {
    let mut e = WasmEngine::new(0.0, 0.0, 0.0);
    let h = e.add_dynamic_body(1.0, 1.0, 0.0, 0.0);
    let q = e.query_sphere_overlap(0.0, 0.0, 0.0, 5.0);
    assert!(q.handles.contains(&h));
}

#[test]
fn test_find_closest_body_empty_returns_none() {
    let e = WasmEngine::new(0.0, -9.81, 0.0);
    assert!(e.find_closest_body(0.0, 0.0, 0.0).is_none());
}

#[test]
fn test_find_closest_body_returns_nearest() {
    let mut e = WasmEngine::new(0.0, 0.0, 0.0);
    let _far = e.add_dynamic_body(1.0, 100.0, 0.0, 0.0);
    let near = e.add_dynamic_body(1.0, 1.0, 0.0, 0.0);
    let closest = e
        .find_closest_body(0.0, 0.0, 0.0)
        .expect("should find a body");
    assert_eq!(closest, near, "nearest body should be at x=1");
}

#[test]
fn test_simulation_summary_empty() {
    let e = WasmEngine::new(0.0, -9.81, 0.0);
    let s = e.simulation_summary();
    assert_eq!(s.active_body_count, 0);
    assert!((s.total_kinetic_energy).abs() < 1e-12);
}

#[test]
fn test_simulation_summary_body_count() {
    let mut e = WasmEngine::new(0.0, -9.81, 0.0);
    e.add_dynamic_body(1.0, 0.0, 10.0, 0.0);
    e.add_dynamic_body(2.0, 1.0, 10.0, 0.0);
    let s = e.simulation_summary();
    assert_eq!(s.active_body_count, 2);
}

#[test]
fn test_is_empty_true_when_empty() {
    let e = WasmEngine::new(0.0, -9.81, 0.0);
    assert!(e.is_empty());
}

#[test]
fn test_is_empty_false_after_add() {
    let mut e = WasmEngine::new(0.0, -9.81, 0.0);
    e.add_dynamic_body(1.0, 0.0, 0.0, 0.0);
    assert!(!e.is_empty());
}

#[test]
fn test_serialise_body_flat_length() {
    let mut e = WasmEngine::new(0.0, -9.81, 0.0);
    let h = e.add_dynamic_body(1.0, 3.0, 4.0, 5.0);
    let flat = e.serialise_body_flat(h);
    assert_eq!(flat.len(), BODY_STATE_FLAT_LEN);
}

#[test]
fn test_serialise_body_flat_position_values() {
    let mut e = WasmEngine::new(0.0, 0.0, 0.0);
    let h = e.add_dynamic_body(1.0, 3.0, 4.0, 5.0);
    let flat = e.serialise_body_flat(h);
    assert!((flat[1] - 3.0).abs() < 1e-8, "px should be 3.0");
    assert!((flat[2] - 4.0).abs() < 1e-8, "py should be 4.0");
    assert!((flat[3] - 5.0).abs() < 1e-8, "pz should be 5.0");
}

#[test]
fn test_serialise_body_flat_invalid_handle_nan() {
    let e = WasmEngine::new(0.0, -9.81, 0.0);
    let flat = e.serialise_body_flat(9999);
    assert!(
        flat[1].is_nan(),
        "invalid handle should produce NAN position"
    );
}

#[test]
fn test_serialise_all_bodies_flat_length() {
    let mut e = WasmEngine::new(0.0, -9.81, 0.0);
    e.add_dynamic_body(1.0, 0.0, 0.0, 0.0);
    e.add_dynamic_body(2.0, 1.0, 0.0, 0.0);
    let flat = e.serialise_all_bodies_flat();
    assert_eq!(flat.len(), 2 * BODY_STATE_FLAT_LEN);
}

#[test]
fn test_deserialise_body_flat_updates_position() {
    let mut e = WasmEngine::new(0.0, 0.0, 0.0);
    let h = e.add_dynamic_body(1.0, 0.0, 0.0, 0.0);
    let mut data = [0.0f64; BODY_STATE_FLAT_LEN];
    data[0] = h as f64;
    data[1] = 7.0;
    data[2] = 8.0;
    data[3] = 9.0; // position
    assert!(e.deserialise_body_flat(&data));
    let p = e.get_body_position(h);
    assert!((p[0] - 7.0).abs() < 1e-8);
}

#[test]
fn test_deserialise_body_flat_too_short_returns_false() {
    let mut e = WasmEngine::new(0.0, 0.0, 0.0);
    assert!(!e.deserialise_body_flat(&[1.0, 2.0]));
}

#[test]
fn test_state_to_json_extended_is_valid_json_prefix() {
    let mut e = WasmEngine::new(0.0, -9.81, 0.0);
    e.add_dynamic_body(1.0, 0.0, 5.0, 0.0);
    let json = e.state_to_json_extended();
    assert!(json.starts_with('{'), "JSON should start with a brace");
    assert!(json.contains("\"time\""), "JSON should contain time field");
    assert!(
        json.contains("\"bodies\""),
        "JSON should contain bodies field"
    );
}

#[test]
fn test_contacts_to_json_empty() {
    let e = WasmEngine::new(0.0, -9.81, 0.0);
    let json = e.contacts_to_json();
    assert_eq!(json, "[]", "empty contact list should be []");
}

#[test]
fn test_get_contacts_flat_extended_empty() {
    let e = WasmEngine::new(0.0, -9.81, 0.0);
    let flat = e.get_contacts_flat_extended();
    assert!(flat.is_empty(), "no contacts expected in empty world");
}

#[test]
fn test_get_all_handles_empty() {
    let e = WasmEngine::new(0.0, -9.81, 0.0);
    assert!(e.get_all_handles().is_empty());
}

#[test]
fn test_get_all_handles_after_add() {
    let mut e = WasmEngine::new(0.0, -9.81, 0.0);
    let h = e.add_dynamic_body(1.0, 0.0, 0.0, 0.0);
    assert!(e.get_all_handles().contains(&h));
}
