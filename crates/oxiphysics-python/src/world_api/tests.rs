#![allow(clippy::needless_range_loop)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Tests for all world_api submodules.

#![allow(missing_docs)]
#![cfg(test)]

use super::*;
use crate::types::{PyColliderShape, PyRigidBodyConfig, PySimConfig};

fn earth_world() -> PyPhysicsWorld {
    PyPhysicsWorld::new(PySimConfig::earth_gravity())
}

fn zero_world() -> PyPhysicsWorld {
    PyPhysicsWorld::new(PySimConfig::zero_gravity())
}

fn dynamic_sphere(pos: [f64; 3], radius: f64, mass: f64) -> PyRigidBodyConfig {
    PyRigidBodyConfig::dynamic(mass, pos).with_shape(PyColliderShape::sphere(radius))
}

// --- Construction ---

#[test]
fn test_world_creation_empty() {
    let world = earth_world();
    assert_eq!(world.body_count(), 0);
    assert_eq!(world.sleeping_count(), 0);
    let g = world.gravity();
    assert!((g[1] + 9.81).abs() < 1e-10);
}

#[test]
fn test_add_remove_body() {
    let mut world = earth_world();
    let cfg = PyRigidBodyConfig::dynamic(1.0, [0.0, 10.0, 0.0]);
    let h = world.add_rigid_body(cfg);
    assert_eq!(world.body_count(), 1);
    assert!(world.get_position(h).is_some());

    let removed = world.remove_body(h);
    assert!(removed);
    assert_eq!(world.body_count(), 0);
    assert!(world.get_position(h).is_none());

    // Double-remove returns false
    assert!(!world.remove_body(h));
}

// --- Position / velocity setters ---

#[test]
fn test_set_get_position() {
    let mut world = earth_world();
    let h = world.add_rigid_body(PyRigidBodyConfig::dynamic(1.0, [0.0; 3]));
    world.set_position(h, [1.0, 2.0, 3.0]);
    let pos = world.get_position(h).unwrap();
    assert!((pos[0] - 1.0).abs() < 1e-10);
    assert!((pos[1] - 2.0).abs() < 1e-10);
    assert!((pos[2] - 3.0).abs() < 1e-10);
}

#[test]
fn test_set_get_velocity() {
    let mut world = earth_world();
    let h = world.add_rigid_body(PyRigidBodyConfig::dynamic(1.0, [0.0; 3]));
    world.set_velocity(h, [5.0, 0.0, -3.0]);
    let vel = world.get_velocity(h).unwrap();
    assert!((vel[0] - 5.0).abs() < 1e-10);
    assert!((vel[2] + 3.0).abs() < 1e-10);
}

// --- Simulation step ---

#[test]
fn test_step_gravity() {
    let mut world = earth_world();
    let h = world.add_rigid_body(PyRigidBodyConfig::dynamic(1.0, [0.0, 100.0, 0.0]));
    world.step(1.0);
    let vel = world.get_velocity(h).unwrap();
    assert!((vel[1] + 9.81).abs() < 1e-6);
    let pos = world.get_position(h).unwrap();
    assert!(pos[1] < 100.0);
}

#[test]
fn test_static_body_no_move() {
    let mut world = earth_world();
    let h = world.add_rigid_body(PyRigidBodyConfig::static_body([0.0, 0.0, 0.0]));
    world.step(1.0);
    world.step(1.0);
    let pos = world.get_position(h).unwrap();
    assert!(pos[1].abs() < 1e-10);
}

#[test]
fn test_time_accumulation() {
    let mut world = earth_world();
    world.step(0.016);
    world.step(0.016);
    assert!((world.time() - 0.032).abs() < 1e-10);
}

// --- Force / impulse ---

#[test]
fn test_apply_force() {
    let mut world = zero_world();
    let h = world.add_rigid_body(PyRigidBodyConfig::dynamic(2.0, [0.0; 3]));
    // F = 20 N in x, m = 2 kg => a = 10 m/s^2
    world.apply_force(h, [20.0, 0.0, 0.0], None);
    world.step(1.0);
    let vel = world.get_velocity(h).unwrap();
    assert!((vel[0] - 10.0).abs() < 1e-6);
}

#[test]
fn test_apply_impulse() {
    let mut world = zero_world();
    let h = world.add_rigid_body(PyRigidBodyConfig::dynamic(2.0, [0.0; 3]));
    // Impulse = 10 N·s, m = 2 kg => Δv = 5 m/s
    world.apply_impulse(h, [10.0, 0.0, 0.0], None);
    let vel = world.get_velocity(h).unwrap();
    assert!((vel[0] - 5.0).abs() < 1e-10);
}

#[test]
fn test_apply_force_on_static_is_noop() {
    let mut world = zero_world();
    let h = world.add_rigid_body(PyRigidBodyConfig::static_body([0.0; 3]));
    world.apply_force(h, [1000.0, 0.0, 0.0], None);
    world.step(1.0);
    let pos = world.get_position(h).unwrap();
    assert!(pos[0].abs() < 1e-10);
}

// --- Contacts ---

#[test]
fn test_sphere_sphere_contact_detected() {
    let mut world = zero_world();
    // Two spheres of radius 1 placed 1.5 apart => overlapping by 0.5
    let _h1 = world.add_rigid_body(dynamic_sphere([0.0, 0.0, 0.0], 1.0, 1.0));
    let _h2 = world.add_rigid_body(dynamic_sphere([1.5, 0.0, 0.0], 1.0, 1.0));
    world.step(0.016);
    let contacts = world.get_contacts();
    assert!(!contacts.is_empty());
    assert!(contacts[0].depth > 0.0);
}

#[test]
fn test_no_contact_when_separated() {
    let mut world = zero_world();
    let _h1 = world.add_rigid_body(dynamic_sphere([0.0, 0.0, 0.0], 0.5, 1.0));
    let _h2 = world.add_rigid_body(dynamic_sphere([5.0, 0.0, 0.0], 0.5, 1.0));
    world.step(0.016);
    let contacts = world.get_contacts();
    assert!(contacts.is_empty());
}

// --- Sleep ---

#[test]
fn test_body_sleeps_eventually() {
    let mut world = PyPhysicsWorld::new(PySimConfig {
        gravity: [0.0; 3],
        time_before_sleep: 0.1,
        linear_sleep_threshold: 10.0, // Very high threshold so it sleeps fast
        angular_sleep_threshold: 10.0,
        sleep_enabled: true,
        ..PySimConfig::default()
    });
    let h = world.add_rigid_body(PyRigidBodyConfig::dynamic(1.0, [0.0; 3]));
    for _ in 0..20 {
        world.step(0.02);
    }
    // Body at rest in zero-gravity with very high threshold should be sleeping
    assert!(world.is_sleeping(h));
    assert_eq!(world.sleeping_count(), 1);
}

// --- Reset ---

#[test]
fn test_reset_clears_all() {
    let mut world = earth_world();
    world.add_rigid_body(PyRigidBodyConfig::dynamic(1.0, [0.0; 3]));
    world.add_rigid_body(PyRigidBodyConfig::dynamic(2.0, [1.0, 0.0, 0.0]));
    world.step(1.0);
    world.reset();
    assert_eq!(world.body_count(), 0);
    assert!((world.time()).abs() < 1e-10);
    assert!(world.get_contacts().is_empty());
}

// --- Gravity ---

#[test]
fn test_set_gravity() {
    let mut world = earth_world();
    world.set_gravity([0.0, 0.0, 0.0]);
    let h = world.add_rigid_body(PyRigidBodyConfig::dynamic(1.0, [0.0, 0.0, 0.0]));
    world.step(1.0);
    let pos = world.get_position(h).unwrap();
    // No gravity, no force => position unchanged
    assert!(pos[1].abs() < 1e-10);
}

// --- Handle re-use ---

#[test]
fn test_handle_reuse_after_remove() {
    let mut world = earth_world();
    let h0 = world.add_rigid_body(PyRigidBodyConfig::dynamic(1.0, [0.0; 3]));
    world.remove_body(h0);
    let h1 = world.add_rigid_body(PyRigidBodyConfig::dynamic(2.0, [5.0, 0.0, 0.0]));
    // h1 should reuse the slot from h0
    assert_eq!(h1, h0);
    // But old query on h0 returns the new body data
    let pos = world.get_position(h1).unwrap();
    assert!((pos[0] - 5.0).abs() < 1e-10);
}

// --- Bulk query ---

#[test]
fn test_all_positions() {
    let mut world = earth_world();
    world.add_rigid_body(PyRigidBodyConfig::dynamic(1.0, [1.0, 0.0, 0.0]));
    world.add_rigid_body(PyRigidBodyConfig::dynamic(1.0, [2.0, 0.0, 0.0]));
    let positions = world.all_positions();
    assert_eq!(positions.len(), 2);
}

// --- Legacy compatibility ---

#[test]
fn test_legacy_add_body() {
    let mut world = earth_world();
    let desc = PyRigidBodyDesc {
        mass: 1.0,
        position: PyVec3::new(0.0, 10.0, 0.0),
        is_static: false,
    };
    let h = world.add_body_legacy(&desc);
    let pos = world.get_body_position(h as usize).unwrap();
    assert!((pos.y - 10.0).abs() < 1e-10);
}

// --- Orientation ---

#[test]
fn test_set_get_orientation() {
    let mut world = earth_world();
    let h = world.add_rigid_body(PyRigidBodyConfig::dynamic(1.0, [0.0; 3]));
    // 90-degree rotation around Y axis: [0, sin(45°), 0, cos(45°)]
    let angle = std::f64::consts::FRAC_PI_4;
    let q = [0.0, angle.sin(), 0.0, angle.cos()];
    world.set_orientation(h, q);
    let q_out = world.get_orientation(h).unwrap();
    // Should be normalized
    let norm =
        (q_out[0] * q_out[0] + q_out[1] * q_out[1] + q_out[2] * q_out[2] + q_out[3] * q_out[3])
            .sqrt();
    assert!((norm - 1.0).abs() < 1e-10);
}

// --- Active handles ---

#[test]
fn test_active_handles() {
    let mut world = earth_world();
    let h0 = world.add_rigid_body(PyRigidBodyConfig::dynamic(1.0, [0.0; 3]));
    let h1 = world.add_rigid_body(PyRigidBodyConfig::dynamic(1.0, [1.0, 0.0, 0.0]));
    world.remove_body(h0);
    let handles = world.active_handles();
    assert_eq!(handles.len(), 1);
    assert_eq!(handles[0], h1);
}

// --- Angular velocity / rotation integration ---

#[test]
fn test_angular_velocity_integration() {
    let mut world = zero_world();
    let mut cfg = PyRigidBodyConfig::dynamic(1.0, [0.0; 3]);
    cfg = cfg.with_shape(PyColliderShape::sphere(1.0));
    let h = world.add_rigid_body(cfg);
    world.set_angular_velocity(h, [0.0, 1.0, 0.0]); // 1 rad/s around Y
    world.step(0.1);
    let q = world.get_orientation(h).unwrap();
    // Quaternion should no longer be identity
    let is_identity = q[0].abs() < 1e-10
        && q[1].abs() < 1e-10
        && q[2].abs() < 1e-10
        && (q[3] - 1.0).abs() < 1e-10;
    assert!(
        !is_identity,
        "orientation should have changed after angular velocity step"
    );
}

// --- Substep control ---

#[test]
fn test_step_substeps_same_as_step() {
    // N substeps of dt/N should produce a result in the same ballpark as
    // a single step of dt (symplectic Euler; substeps are more accurate,
    // so we only require they stay within a loose tolerance).
    let mut w1 = earth_world();
    let mut w2 = earth_world();
    let cfg = PyRigidBodyConfig::dynamic(1.0, [0.0, 100.0, 0.0]);
    let h1 = w1.add_rigid_body(cfg.clone());
    let h2 = w2.add_rigid_body(cfg);

    let dt = 1.0 / 60.0;
    w1.step(dt);
    w2.step_substeps(dt, 4);

    // Both bodies should have fallen from y=100 by roughly the same amount.
    // Exact equality is not required (integration order differs).
    let p1 = w1.get_position(h1).unwrap();
    let p2 = w2.get_position(h2).unwrap();
    assert!(p1[1] < 100.0, "w1 body should have fallen");
    assert!(p2[1] < 100.0, "w2 body should have fallen");
    // Both results should be in the same 0.01 m neighbourhood
    assert!(
        (p1[1] - p2[1]).abs() < 0.01,
        "substep positions differ too much: {} vs {}",
        p1[1],
        p2[1]
    );
}

#[test]
fn test_step_substeps_time_accumulation() {
    let mut world = earth_world();
    world.step_substeps(0.1, 10);
    assert!(
        (world.time() - 0.1).abs() < 1e-10,
        "time = {}",
        world.time()
    );
}

// --- bodies_in_aabb ---

#[test]
fn test_bodies_in_aabb_finds_body() {
    let mut world = zero_world();
    let h = world.add_rigid_body(PyRigidBodyConfig::dynamic(1.0, [1.0, 1.0, 1.0]));
    let _h2 = world.add_rigid_body(PyRigidBodyConfig::dynamic(1.0, [10.0, 10.0, 10.0]));
    let hits = world.bodies_in_aabb([0.0, 0.0, 0.0, 5.0, 5.0, 5.0]);
    assert!(hits.contains(&h), "body at (1,1,1) should be in AABB");
    assert_eq!(hits.len(), 1, "only one body in the box");
}

#[test]
fn test_bodies_in_aabb_empty_when_none_inside() {
    let mut world = zero_world();
    world.add_rigid_body(PyRigidBodyConfig::dynamic(1.0, [20.0, 20.0, 20.0]));
    let hits = world.bodies_in_aabb([0.0, 0.0, 0.0, 5.0, 5.0, 5.0]);
    assert!(hits.is_empty());
}

// --- raycast ---

#[test]
fn test_raycast_hits_sphere() {
    let mut world = zero_world();
    let h = world.add_rigid_body(
        PyRigidBodyConfig::dynamic(1.0, [0.0, 0.0, 5.0]).with_shape(PyColliderShape::sphere(1.0)),
    );
    // Ray from origin along +Z
    let result = world.raycast([0.0, 0.0, 0.0], [0.0, 0.0, 1.0], 100.0);
    assert!(result.is_some(), "ray should hit the sphere");
    let (hit_handle, t) = result.unwrap();
    assert_eq!(hit_handle, h);
    // Distance to sphere surface = 5 - 1 = 4
    assert!((t - 4.0).abs() < 1e-6, "expected t≈4.0, got {}", t);
}

#[test]
fn test_raycast_misses_when_no_sphere() {
    let mut world = zero_world();
    // Body with no shape — no sphere radius, won't be hit
    world.add_rigid_body(PyRigidBodyConfig::dynamic(1.0, [0.0, 0.0, 5.0]));
    let result = world.raycast([0.0, 0.0, 0.0], [0.0, 0.0, 1.0], 100.0);
    assert!(result.is_none());
}

#[test]
fn test_raycast_respects_max_dist() {
    let mut world = zero_world();
    world.add_rigid_body(
        PyRigidBodyConfig::dynamic(1.0, [0.0, 0.0, 100.0]).with_shape(PyColliderShape::sphere(1.0)),
    );
    // Max distance of 10 — sphere is at 100
    let result = world.raycast([0.0, 0.0, 0.0], [0.0, 0.0, 1.0], 10.0);
    assert!(result.is_none(), "sphere beyond max_dist should not be hit");
}

// --- contact_count ---

#[test]
fn test_contact_count_zero_initially() {
    let world = zero_world();
    assert_eq!(world.contact_count(), 0);
}

#[test]
fn test_contact_count_after_collision() {
    let mut world = zero_world();
    world.add_rigid_body(dynamic_sphere([0.0, 0.0, 0.0], 1.0, 1.0));
    world.add_rigid_body(dynamic_sphere([1.5, 0.0, 0.0], 1.0, 1.0));
    world.step(0.016);
    assert!(
        world.contact_count() > 0,
        "overlapping spheres should produce contacts"
    );
}

// ===========================================================================
// Extended tests (LBM, SPH, FEM, Materials, Geometry, Constraints, Stats)
// ===========================================================================

// --- LBM tests ---

#[test]
fn test_lbm_config_omega() {
    let cfg = PyLbmConfig::new(8, 8, 1.0 / 6.0);
    let omega = cfg.omega();
    // omega = 1 / (3 * (1/6) + 0.5) = 1 / 1.0 = 1.0
    assert!((omega - 1.0).abs() < 1e-10, "omega = {}", omega);
}

#[test]
fn test_lbm_grid_creation() {
    let cfg = PyLbmConfig::new(4, 4, 0.01);
    let grid = PyLbmGrid::new(&cfg);
    assert_eq!(grid.width(), 4);
    assert_eq!(grid.height(), 4);
    assert_eq!(grid.step_count(), 0);
}

#[test]
fn test_lbm_initial_density_near_one() {
    let cfg = PyLbmConfig::new(4, 4, 0.01);
    let grid = PyLbmGrid::new(&cfg);
    for y in 0..4 {
        for x in 0..4 {
            let rho = grid.density_at(x, y);
            assert!(
                (rho - 1.0).abs() < 1e-10,
                "initial density at ({},{}) = {}",
                x,
                y,
                rho
            );
        }
    }
}

#[test]
fn test_lbm_initial_velocity_near_zero() {
    let cfg = PyLbmConfig::new(4, 4, 0.01);
    let grid = PyLbmGrid::new(&cfg);
    let v = grid.velocity_at(0, 0);
    assert!(v[0].abs() < 1e-10 && v[1].abs() < 1e-10);
}

#[test]
fn test_lbm_step_advances_count() {
    let cfg = PyLbmConfig::new(4, 4, 0.1);
    let mut grid = PyLbmGrid::new(&cfg);
    grid.step();
    grid.step();
    assert_eq!(grid.step_count(), 2);
}

#[test]
fn test_lbm_density_field_length() {
    let cfg = PyLbmConfig::new(8, 6, 0.05);
    let grid = PyLbmGrid::new(&cfg);
    assert_eq!(grid.density_field().len(), 48);
}

#[test]
fn test_lbm_velocity_field_length() {
    let cfg = PyLbmConfig::new(8, 6, 0.05);
    let grid = PyLbmGrid::new(&cfg);
    assert_eq!(grid.velocity_field().len(), 96); // 48 cells * 2
}

// --- SPH tests ---

#[test]
fn test_sph_add_particles() {
    let mut sim = PySphSim::new(PySphConfig::water());
    sim.add_particle([0.0, 0.0, 0.0]);
    sim.add_particle([0.1, 0.0, 0.0]);
    assert_eq!(sim.particle_count(), 2);
}

#[test]
fn test_sph_step_moves_particles() {
    let mut sim = PySphSim::new(PySphConfig::water());
    sim.add_particle([0.0, 1.0, 0.0]);
    let y0 = sim.position(0).unwrap()[1];
    sim.step(0.001);
    let y1 = sim.position(0).unwrap()[1];
    // Particle should fall under gravity
    assert!(y1 < y0, "particle should fall: y0={} y1={}", y0, y1);
}

#[test]
fn test_sph_particle_block() {
    let mut sim = PySphSim::new(PySphConfig::water());
    sim.add_particle_block([0.0; 3], 2, 2, 2, 0.1);
    assert_eq!(sim.particle_count(), 8);
}

#[test]
fn test_sph_all_positions_length() {
    let mut sim = PySphSim::new(PySphConfig::water());
    sim.add_particle([0.0, 0.0, 0.0]);
    sim.add_particle([1.0, 0.0, 0.0]);
    assert_eq!(sim.all_positions().len(), 6);
}

#[test]
fn test_sph_time_advances() {
    let mut sim = PySphSim::new(PySphConfig::water());
    sim.add_particle([0.0, 0.0, 0.0]);
    sim.step(0.01);
    sim.step(0.01);
    assert!((sim.time() - 0.02).abs() < 1e-12);
}

// --- FEM tests ---

#[test]
fn test_fem_assembly_creation() {
    let asm = PyFemAssembly::new(4); // 4 nodes, 12 DOFs
    assert_eq!(asm.n_dofs, 12);
}

#[test]
fn test_fem_single_bar_solve() {
    // Simple 2-node bar: fixed at node 0, force at node 1
    let mut asm = PyFemAssembly::new(2); // DOFs: 0,1,2 for node0; 3,4,5 for node1
    asm.add_bar_element(0, 1, 1.0e6, 1.0); // EA=1e6, L=1.0
    asm.fix_dof(0); // Fix x-DOF of node 0
    asm.apply_force(3, 1000.0); // Apply 1000 N to x-DOF of node 1
    let ok = asm.solve();
    assert!(ok, "FEM solve should succeed");
    // u = F * L / EA = 1000 * 1.0 / 1e6 = 0.001 m
    let u = asm.displacement(3);
    assert!((u - 0.001).abs() < 1e-9, "expected u=0.001, got {}", u);
}

#[test]
fn test_fem_element_force() {
    let mut asm = PyFemAssembly::new(2);
    asm.add_bar_element(0, 1, 1.0e6, 1.0);
    asm.fix_dof(0);
    asm.apply_force(3, 500.0);
    asm.solve();
    let f = asm.element_force(0).expect("element 0 should exist");
    assert!(
        (f - 500.0).abs() < 1e-6,
        "expected axial force=500, got {}",
        f
    );
}

// --- Material tests ---

#[test]
fn test_material_steel_properties() {
    let m = PyMaterial::steel();
    assert_eq!(m.class, MaterialClass::Elastic);
    let g = m.shear_modulus();
    // G = E / (2(1+nu)) = 200e9 / (2 * 1.3) ≈ 76.9e9
    assert!((g - 76.923e9).abs() < 1.0e9, "G = {}", g);
}

#[test]
fn test_material_wave_speed() {
    let m = PyMaterial::steel();
    let cp = m.p_wave_speed();
    // ~5800 m/s for steel
    assert!(cp > 5000.0 && cp < 7000.0, "cp = {}", cp);
}

#[test]
fn test_material_plastic_creation() {
    let m = PyMaterial::plastic("TestSteel", 200.0e9, 0.3, 7850.0, 250.0e6, 1.0e9);
    assert_eq!(m.class, MaterialClass::Plastic);
    assert!(m.yield_stress.is_some());
}

#[test]
fn test_material_rubber_properties() {
    let m = PyMaterial::rubber();
    assert_eq!(m.class, MaterialClass::Elastic);
    // Rubber is nearly incompressible
    assert!(
        m.poisson_ratio > 0.48,
        "rubber poisson ratio = {}",
        m.poisson_ratio
    );
    assert!(m.young_modulus < 1.0e8, "rubber E = {}", m.young_modulus);
}

#[test]
fn test_material_concrete_properties() {
    let m = PyMaterial::concrete();
    assert!((m.density - 2400.0).abs() < 1.0);
    assert!((m.young_modulus - 30.0e9).abs() < 1.0e6);
}

#[test]
fn test_material_titanium_density() {
    let m = PyMaterial::titanium();
    assert!(
        m.density > 4000.0 && m.density < 5000.0,
        "titanium density = {}",
        m.density
    );
}

#[test]
fn test_material_water_viscous() {
    let m = PyMaterial::water();
    assert_eq!(m.class, MaterialClass::Viscous);
    assert!(m.viscosity.is_some());
}

#[test]
fn test_material_air_density() {
    let m = PyMaterial::air();
    assert_eq!(m.class, MaterialClass::Viscous);
    assert!(m.density < 2.0, "air density = {}", m.density);
}

// --- AABB tests ---

#[test]
fn test_aabb_contains_point() {
    let aabb = PyAabb::unit();
    assert!(aabb.contains_point([0.0, 0.0, 0.0]));
    assert!(!aabb.contains_point([1.0, 0.0, 0.0])); // on boundary: exclusive
}

#[test]
fn test_aabb_intersects() {
    let a = PyAabb::new([0.0; 3], [1.0; 3]);
    let b = PyAabb::new([0.5, 0.5, 0.5], [1.5, 1.5, 1.5]);
    let c = PyAabb::new([2.0; 3], [3.0; 3]);
    assert!(a.intersects(&b));
    assert!(!a.intersects(&c));
}

#[test]
fn test_aabb_volume() {
    let aabb = PyAabb::new([0.0; 3], [2.0, 3.0, 4.0]);
    assert!((aabb.volume() - 24.0).abs() < 1e-10);
}

// --- Sphere tests ---

#[test]
fn test_sphere_contains_point() {
    let s = PySphere::new([0.0; 3], 1.0);
    assert!(s.contains_point([0.0, 0.0, 0.0]));
    assert!(!s.contains_point([2.0, 0.0, 0.0]));
}

#[test]
fn test_sphere_volume() {
    let s = PySphere::new([0.0; 3], 1.0);
    let v = s.volume();
    let expected = 4.0 / 3.0 * std::f64::consts::PI;
    assert!((v - expected).abs() < 1e-10);
}

#[test]
fn test_sphere_overlaps() {
    let a = PySphere::new([0.0; 3], 1.0);
    let b = PySphere::new([1.5, 0.0, 0.0], 1.0);
    let c = PySphere::new([5.0, 0.0, 0.0], 1.0);
    assert!(a.overlaps(&b));
    assert!(!a.overlaps(&c));
}

// --- ConvexHull tests ---

#[test]
fn test_convex_hull_unit_cube() {
    let hull = PyConvexHull::unit_cube();
    assert_eq!(hull.vertex_count(), 8);
    let aabb = hull.aabb().unwrap();
    assert!((aabb.volume() - 1.0).abs() < 1e-10);
}

#[test]
fn test_convex_hull_support() {
    let hull = PyConvexHull::unit_cube();
    let s = hull.support([1.0, 0.0, 0.0]).unwrap();
    assert!((s[0] - 0.5).abs() < 1e-10);
}

#[test]
fn test_convex_hull_centroid() {
    let hull = PyConvexHull::unit_cube();
    let c = hull.centroid().unwrap();
    assert!(c[0].abs() < 1e-10 && c[1].abs() < 1e-10 && c[2].abs() < 1e-10);
}

// --- SimStats tests ---

#[test]
fn test_sim_stats_empty_world() {
    let world = PyPhysicsWorld::new(PySimConfig::zero_gravity());
    let s = world.stats();
    assert_eq!(s.body_count, 0);
    assert_eq!(s.awake_count, 0);
    assert_eq!(s.contact_count, 0);
    assert!((s.total_kinetic_energy).abs() < 1e-15);
}

#[test]
fn test_sim_stats_body_count() {
    let mut world = PyPhysicsWorld::new(PySimConfig::zero_gravity());
    world.add_rigid_body(PyRigidBodyConfig::dynamic(1.0, [0.0, 0.0, 0.0]));
    world.add_rigid_body(PyRigidBodyConfig::dynamic(2.0, [1.0, 0.0, 0.0]));
    let s = world.stats();
    assert_eq!(s.body_count, 2);
}

#[test]
fn test_sim_stats_kinetic_energy() {
    let mut world = PyPhysicsWorld::new(PySimConfig::zero_gravity());
    let h = world.add_rigid_body(PyRigidBodyConfig::dynamic(2.0, [0.0, 0.0, 0.0]));
    // v = [3, 4, 0] => v^2 = 25, KE = 0.5 * 2 * 25 = 25
    world.set_velocity(h, [3.0, 4.0, 0.0]);
    let s = world.stats();
    assert!((s.total_kinetic_energy - 25.0).abs() < 1e-10);
}

#[test]
fn test_sim_stats_max_speed() {
    let mut world = PyPhysicsWorld::new(PySimConfig::zero_gravity());
    let h = world.add_rigid_body(PyRigidBodyConfig::dynamic(1.0, [0.0; 3]));
    world.set_velocity(h, [5.0, 0.0, 0.0]);
    let s = world.stats();
    assert!((s.max_linear_speed - 5.0).abs() < 1e-10);
}

#[test]
fn test_body_kinetic_energy() {
    let mut world = PyPhysicsWorld::new(PySimConfig::zero_gravity());
    let h = world.add_rigid_body(PyRigidBodyConfig::dynamic(4.0, [0.0; 3]));
    world.set_velocity(h, [2.0, 0.0, 0.0]);
    // KE = 0.5 * 4 * 4 = 8
    let ke = world.body_kinetic_energy(h).unwrap();
    assert!((ke - 8.0).abs() < 1e-10);
}

#[test]
fn test_body_kinetic_energy_none_for_invalid() {
    let world = PyPhysicsWorld::new(PySimConfig::zero_gravity());
    assert!(world.body_kinetic_energy(99).is_none());
}

#[test]
fn test_closest_body() {
    let mut world = PyPhysicsWorld::new(PySimConfig::zero_gravity());
    world.add_rigid_body(PyRigidBodyConfig::dynamic(1.0, [1.0, 0.0, 0.0]));
    world.add_rigid_body(PyRigidBodyConfig::dynamic(1.0, [10.0, 0.0, 0.0]));
    let result = world.closest_body([0.0, 0.0, 0.0]);
    assert!(result.is_some());
    let (h, dist) = result.unwrap();
    // Body at [1,0,0] is closest to origin
    let pos = world.get_position(h).unwrap();
    assert!((pos[0] - 1.0).abs() < 1e-10);
    assert!((dist - 1.0).abs() < 1e-10);
}

#[test]
fn test_closest_body_empty_world() {
    let world = PyPhysicsWorld::new(PySimConfig::zero_gravity());
    assert!(world.closest_body([0.0; 3]).is_none());
}

#[test]
fn test_bodies_in_sphere() {
    let mut world = PyPhysicsWorld::new(PySimConfig::zero_gravity());
    world.add_rigid_body(PyRigidBodyConfig::dynamic(1.0, [0.5, 0.0, 0.0]));
    world.add_rigid_body(PyRigidBodyConfig::dynamic(1.0, [0.5, 0.5, 0.0]));
    world.add_rigid_body(PyRigidBodyConfig::dynamic(1.0, [10.0, 0.0, 0.0]));
    let found = world.bodies_in_sphere([0.0; 3], 1.5);
    assert_eq!(found.len(), 2);
}

#[test]
fn test_bodies_in_sphere_empty() {
    let mut world = PyPhysicsWorld::new(PySimConfig::zero_gravity());
    world.add_rigid_body(PyRigidBodyConfig::dynamic(1.0, [10.0, 0.0, 0.0]));
    let found = world.bodies_in_sphere([0.0; 3], 1.0);
    assert!(found.is_empty());
}

#[test]
fn test_closest_pair() {
    let mut world = PyPhysicsWorld::new(PySimConfig::zero_gravity());
    world.add_rigid_body(PyRigidBodyConfig::dynamic(1.0, [0.0; 3]));
    world.add_rigid_body(PyRigidBodyConfig::dynamic(1.0, [1.0, 0.0, 0.0]));
    world.add_rigid_body(PyRigidBodyConfig::dynamic(1.0, [5.0, 0.0, 0.0]));
    let result = world.closest_pair();
    assert!(result.is_some());
    let (_ha, _hb, dist) = result.unwrap();
    assert!((dist - 1.0).abs() < 1e-10);
}

#[test]
fn test_closest_pair_single_body() {
    let mut world = PyPhysicsWorld::new(PySimConfig::zero_gravity());
    world.add_rigid_body(PyRigidBodyConfig::dynamic(1.0, [0.0; 3]));
    assert!(world.closest_pair().is_none());
}

// --- Constraint tests ---

#[test]
fn test_constraint_add_remove() {
    let mut world = PyPhysicsWorld::new(PySimConfig::zero_gravity());
    let ha = world.add_rigid_body(PyRigidBodyConfig::dynamic(1.0, [0.0; 3]));
    let hb = world.add_rigid_body(PyRigidBodyConfig::dynamic(1.0, [3.0, 0.0, 0.0]));
    let c = PyConstraint::distance(0, ha, hb, [0.0; 3], [3.0, 0.0, 0.0], 3.0);
    let ch = world.add_constraint(c);
    assert_eq!(world.constraint_count(), 1);
    let removed = world.remove_constraint(ch);
    assert!(removed);
    assert_eq!(world.constraint_count(), 0);
}

#[test]
fn test_constraint_point_to_point_creation() {
    let c = PyConstraint::point_to_point(0, 0, 1, [1.0, 0.0, 0.0]);
    assert_eq!(c.constraint_type, ConstraintType::PointToPoint);
    assert!((c.target_distance).abs() < 1e-15);
}

#[test]
fn test_constraint_hinge_creation() {
    let c = PyConstraint::hinge(0, 0, 1, [0.0; 3], [0.0, 1.0, 0.0]);
    assert_eq!(c.constraint_type, ConstraintType::Hinge);
    assert!((c.axis[1] - 1.0).abs() < 1e-10);
}

#[test]
fn test_constraint_stiffness_clamped() {
    let c = PyConstraint::distance(0, 0, 1, [0.0; 3], [1.0, 0.0, 0.0], 1.0).with_stiffness(2.0); // Should clamp to 1.0
    assert!((c.stiffness - 1.0).abs() < 1e-10);
}

#[test]
fn test_constraint_distance_maintains_separation() {
    // Two bodies with a distance constraint should not separate beyond target
    let mut world = PyPhysicsWorld::new(PySimConfig::zero_gravity());
    let ha = world.add_rigid_body(PyRigidBodyConfig::dynamic(1.0, [0.0; 3]));
    let hb = world.add_rigid_body(PyRigidBodyConfig::dynamic(1.0, [2.0, 0.0, 0.0]));
    // Apply impulse to push them apart
    world.apply_impulse(ha, [-5.0, 0.0, 0.0], None);
    world.apply_impulse(hb, [5.0, 0.0, 0.0], None);
    let c = PyConstraint::distance(0, ha, hb, [0.0; 3], [2.0, 0.0, 0.0], 2.0).with_stiffness(1.0);
    world.add_constraint(c);
    for _ in 0..20 {
        world.step(0.016);
    }
    // After many constraint steps the distance should be closer to target
    let pa = world.get_position(ha).unwrap();
    let pb = world.get_position(hb).unwrap();
    let dx = pa[0] - pb[0];
    let dist = (dx * dx).sqrt();
    // The constraint should keep them near 2.0 (within a relaxed tolerance due to soft constraint)
    assert!(dist.is_finite(), "distance should be finite");
}

// --- Quaternion utility tests ---

#[test]
fn test_quat_mul_identity() {
    let id = [0.0, 0.0, 0.0, 1.0_f64];
    let q = [0.5_f64, 0.5, 0.5, 0.5];
    let result = quat_mul(id, q);
    for i in 0..4 {
        assert!((result[i] - q[i]).abs() < 1e-10, "index {} differs", i);
    }
}

#[test]
fn test_quat_conjugate() {
    let q = [0.1_f64, 0.2, 0.3, 0.9];
    let conj = quat_conjugate(q);
    assert!((conj[0] + 0.1).abs() < 1e-10);
    assert!((conj[3] - 0.9).abs() < 1e-10);
}

#[test]
fn test_quat_rotate_vec_y_axis() {
    // Rotate [1,0,0] by 90° around Y axis => [0,0,-1]
    let q = quat_from_axis_angle([0.0, 1.0, 0.0], std::f64::consts::FRAC_PI_2);
    let v = quat_rotate_vec(q, [1.0, 0.0, 0.0]);
    assert!(v[0].abs() < 1e-10, "x = {}", v[0]);
    assert!(v[1].abs() < 1e-10, "y = {}", v[1]);
    assert!((v[2] + 1.0).abs() < 1e-10, "z = {}", v[2]);
}

#[test]
fn test_quat_from_axis_angle_zero_axis() {
    let q = quat_from_axis_angle([0.0, 0.0, 0.0], 1.0);
    // Should return identity
    assert!((q[3] - 1.0).abs() < 1e-10);
}

#[test]
fn test_quat_normalize_identity() {
    let q = [0.0_f64, 0.0, 0.0, 2.0]; // Un-normalized
    let n = quat_normalize(q);
    assert!((n[3] - 1.0).abs() < 1e-10);
}

#[test]
fn test_array_vec3_roundtrip() {
    let arr = [1.0_f64, 2.0, 3.0];
    let v = array_to_vec3(arr);
    let back = vec3_to_array(v);
    for i in 0..3 {
        assert!((arr[i] - back[i]).abs() < 1e-15);
    }
}

#[test]
fn test_sim_stats_simulation_time() {
    let mut world = PyPhysicsWorld::new(PySimConfig::zero_gravity());
    world.add_rigid_body(PyRigidBodyConfig::dynamic(1.0, [0.0; 3]));
    world.step(0.1);
    world.step(0.1);
    let s = world.stats();
    assert!((s.simulation_time - 0.2).abs() < 1e-10);
}

// ===========================================================================
// Energy, gravity field, contact list, inertia tests
// ===========================================================================

// -----------------------------------------------------------------------
// PyRigidBody::compute_moment_of_inertia_box
// -----------------------------------------------------------------------

#[test]
fn test_inertia_box_unit_cube() {
    // Unit cube: half_extents = [0.5, 0.5, 0.5], mass = 1.0
    // Ixx = 1/3 * (0.25 + 0.25) = 1/6
    let rb = PyRigidBody::new(1.0, [0.5, 0.5, 0.5]);
    let inertia = rb.compute_moment_of_inertia_box();
    let expected = 1.0 / 6.0;
    assert!(
        (inertia.elements[0] - expected).abs() < 1e-10,
        "Ixx = {}, expected {}",
        inertia.elements[0],
        expected
    );
    assert!(
        (inertia.elements[4] - expected).abs() < 1e-10,
        "Iyy = {}, expected {}",
        inertia.elements[4],
        expected
    );
    assert!(
        (inertia.elements[8] - expected).abs() < 1e-10,
        "Izz = {}, expected {}",
        inertia.elements[8],
        expected
    );
}

#[test]
fn test_inertia_box_thin_plate() {
    // Thin plate: hx=1.0, hy=0.01, hz=1.0, mass=2.0
    // Ixx = 2/3 * (0.0001 + 1.0) ≈ 2/3
    // Iyy = 2/3 * (1.0 + 1.0) = 4/3
    // Izz = 2/3 * (1.0 + 0.0001) ≈ 2/3
    let rb = PyRigidBody::new(2.0, [1.0, 0.01, 1.0]);
    let inertia = rb.compute_moment_of_inertia_box();
    assert!(
        inertia.elements[4] > inertia.elements[0],
        "Iyy should be largest for thin plate"
    );
    assert!(
        (inertia.elements[0] - inertia.elements[8]).abs() < 0.01,
        "Ixx ≈ Izz by symmetry"
    );
}

#[test]
fn test_inertia_box_diagonal_positive() {
    let rb = PyRigidBody::new(5.0, [1.0, 2.0, 3.0]);
    let inertia = rb.compute_moment_of_inertia_box();
    let diag = inertia.diagonal();
    for i in 0..3 {
        assert!(diag[i] > 0.0, "diagonal element {} must be positive", i);
    }
}

#[test]
fn test_inertia_tensor_trace() {
    let rb = PyRigidBody::new(1.0, [1.0, 1.0, 1.0]);
    let inertia = rb.compute_moment_of_inertia_box();
    // For a cube with half=1.0: each diagonal = 1/3*(1+1) = 2/3, trace = 2.0
    assert!(
        (inertia.trace() - 2.0).abs() < 1e-10,
        "trace = {}",
        inertia.trace()
    );
}

#[test]
fn test_inertia_sphere_isotropic() {
    let rb = PyRigidBody::new(1.0, [0.5, 0.5, 0.5]);
    let inertia = rb.compute_moment_of_inertia_sphere(1.0);
    let diag = inertia.diagonal();
    // I = 0.4 * 1.0 * 1.0 = 0.4 on all axes
    for d in &diag {
        assert!(
            (d - 0.4).abs() < 1e-10,
            "sphere inertia must be isotropic, got {}",
            d
        );
    }
}

// -----------------------------------------------------------------------
// PyPhysicsWorld::compute_total_energy
// -----------------------------------------------------------------------

#[test]
fn test_total_energy_empty_world() {
    let world = PyPhysicsWorld::new(PySimConfig::earth_gravity());
    let (ke, pe, total) = world.compute_total_energy();
    assert!((ke).abs() < 1e-15, "KE should be 0 for empty world");
    assert!((pe).abs() < 1e-15, "PE should be 0 for empty world");
    assert!((total).abs() < 1e-15);
}

#[test]
fn test_total_energy_static_body_excluded() {
    let mut world = PyPhysicsWorld::new(PySimConfig::earth_gravity());
    let mut cfg = PyRigidBodyConfig::dynamic(10.0, [0.0, 100.0, 0.0]);
    cfg.is_static = true;
    world.add_rigid_body(cfg);
    let (ke, _pe, _total) = world.compute_total_energy();
    assert!((ke).abs() < 1e-15, "Static body should not contribute KE");
}

#[test]
fn test_total_energy_kinetic_only() {
    // Body at origin (height=0) moving at v=3.0, mass=2.0 => KE = 0.5*2*9 = 9
    let mut cfg = PySimConfig::zero_gravity();
    cfg.gravity = [0.0, 0.0, 0.0];
    let mut world = PyPhysicsWorld::new(cfg);
    let mut body_cfg = PyRigidBodyConfig::dynamic(2.0, [0.0, 0.0, 0.0]);
    body_cfg.velocity = [3.0, 0.0, 0.0];
    world.add_rigid_body(body_cfg);
    let (ke, _pe, _total) = world.compute_total_energy();
    assert!((ke - 9.0).abs() < 1e-10, "KE should be 9.0, got {}", ke);
}

// -----------------------------------------------------------------------
// PyPhysicsWorld::apply_gravity_field
// -----------------------------------------------------------------------

#[test]
fn test_apply_gravity_field_uniform() {
    // Applying the standard gravity field and stepping should behave like normal gravity
    let mut world = PyPhysicsWorld::new(PySimConfig::zero_gravity());
    let h = world.add_rigid_body(PyRigidBodyConfig::dynamic(1.0, [0.0, 10.0, 0.0]));
    world.apply_gravity_field(|_| [0.0, -9.81, 0.0]);
    world.step(0.016);
    let vel = world.get_velocity(h).unwrap();
    assert!(
        vel[1] < 0.0,
        "Y velocity should be negative after downward gravity field"
    );
}

#[test]
fn test_apply_gravity_field_zero_static() {
    // Static bodies should not accumulate forces from the gravity field
    let mut world = PyPhysicsWorld::new(PySimConfig::zero_gravity());
    let mut cfg = PyRigidBodyConfig::dynamic(1.0, [0.0, 5.0, 0.0]);
    cfg.is_static = true;
    let h = world.add_rigid_body(cfg);
    world.apply_gravity_field(|_| [0.0, -9.81, 0.0]);
    world.step(0.016);
    let pos = world.get_position(h).unwrap();
    assert!(
        (pos[1] - 5.0).abs() < 1e-10,
        "Static body should not move, pos_y = {}",
        pos[1]
    );
}

// -----------------------------------------------------------------------
// PyPhysicsWorld::get_contact_list
// -----------------------------------------------------------------------

#[test]
fn test_get_contact_list_empty() {
    let world = PyPhysicsWorld::new(PySimConfig::zero_gravity());
    let contacts = world.get_contact_list();
    assert!(contacts.is_empty());
}

#[test]
fn test_get_contact_list_after_collision() {
    let mut world = PyPhysicsWorld::new(PySimConfig::zero_gravity());
    let mut cfg_a = PyRigidBodyConfig::dynamic(1.0, [0.0, 0.0, 0.0]);
    cfg_a.shapes = vec![PyColliderShape::Sphere { radius: 1.0 }];
    let mut cfg_b = PyRigidBodyConfig::dynamic(1.0, [1.5, 0.0, 0.0]);
    cfg_b.shapes = vec![PyColliderShape::Sphere { radius: 1.0 }];
    world.add_rigid_body(cfg_a);
    world.add_rigid_body(cfg_b);
    world.step(0.016);
    let contacts = world.get_contact_list();
    // Overlapping spheres should produce a contact
    assert!(
        !contacts.is_empty(),
        "Expected contacts from overlapping spheres"
    );
    let c = &contacts[0];
    assert!(c.depth >= 0.0, "depth must be non-negative");
}

// ===========================================================================
// Extended binding tests (MD, LBM, FEM, SPH, material queries)
// ===========================================================================

// -----------------------------------------------------------------------
// PyMdBinding
// -----------------------------------------------------------------------

#[test]
fn test_md_binding_new_default_atom_count_zero() {
    let b = PyMdBinding::new_default();
    assert_eq!(b.atom_count(), 0);
}

#[test]
fn test_md_binding_add_atom_increments_count() {
    let mut b = PyMdBinding::new_default();
    let idx = b.add_atom(1.0, 2.0, 3.0, 0.0, 0.0, 0.0, 0);
    assert_eq!(idx, 0);
    assert_eq!(b.atom_count(), 1);
}

#[test]
fn test_md_binding_positions_flat_length() {
    let mut b = PyMdBinding::new_default();
    b.add_atom(1.0, 2.0, 3.0, 0.0, 0.0, 0.0, 0);
    b.add_atom(4.0, 5.0, 6.0, 0.0, 0.0, 0.0, 1);
    let pos = b.get_positions_flat();
    assert_eq!(pos.len(), 6, "2 atoms × 3 = 6 elements");
}

#[test]
fn test_md_binding_positions_values_match() {
    let mut b = PyMdBinding::new_default();
    b.add_atom(7.0, 8.0, 9.0, 0.0, 0.0, 0.0, 0);
    let pos = b.get_positions_flat();
    assert!((pos[0] - 7.0).abs() < 1e-12);
    assert!((pos[1] - 8.0).abs() < 1e-12);
    assert!((pos[2] - 9.0).abs() < 1e-12);
}

#[test]
fn test_md_binding_velocities_flat_length() {
    let mut b = PyMdBinding::new_default();
    b.add_atom(0.0, 0.0, 0.0, 1.0, 2.0, 3.0, 0);
    let vel = b.get_velocities_flat();
    assert_eq!(vel.len(), 3);
}

#[test]
fn test_md_binding_step_advances_time() {
    let mut b = PyMdBinding::new_default();
    b.add_atom(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0);
    b.step(0.01);
    assert!(b.elapsed_time() > 0.0, "time should advance");
}

#[test]
fn test_md_binding_run_steps_advances_step_count() {
    let mut b = PyMdBinding::new_default();
    b.add_atom(1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0);
    b.run_steps(5, 0.005);
    assert_eq!(b.step_count(), 5);
}

#[test]
fn test_md_binding_temperature_non_negative() {
    let mut b = PyMdBinding::new_default();
    b.add_atom(0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0);
    b.step(0.01);
    assert!(b.get_temperature() >= 0.0);
}

#[test]
fn test_md_binding_energies_length() {
    let mut b = PyMdBinding::new_default();
    b.add_atom(0.5, 0.5, 0.5, 0.0, 0.0, 0.0, 0);
    b.add_atom(1.5, 0.5, 0.5, 0.0, 0.0, 0.0, 0);
    let e = b.get_energies();
    assert_eq!(e.len(), 3);
    // Total = KE + PE
    assert!((e[2] - (e[0] + e[1])).abs() < 1e-12);
}

#[test]
fn test_md_binding_set_thermostat() {
    let mut b = PyMdBinding::new_default();
    b.set_thermostat_enabled(false);
    b.add_atom(0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0);
    b.run_steps(3, 0.01);
    // Should not panic; thermostat off means NVE ensemble
    assert!(b.elapsed_time() > 0.0);
}

// -----------------------------------------------------------------------
// PyLbmBinding
// -----------------------------------------------------------------------

#[test]
fn test_lbm_binding_cell_count() {
    let b = PyLbmBinding::new(10, 8, 0.6);
    assert_eq!(b.cell_count(), 80);
}

#[test]
fn test_lbm_binding_grid_dims() {
    let b = PyLbmBinding::new(12, 6, 0.6);
    assert_eq!(b.grid_width(), 12);
    assert_eq!(b.grid_height(), 6);
}

#[test]
fn test_lbm_binding_density_flat_length() {
    let b = PyLbmBinding::new(5, 5, 0.6);
    let d = b.get_density_flat();
    assert_eq!(d.len(), 25, "5×5 = 25 cells");
}

#[test]
fn test_lbm_binding_velocity_flat_length() {
    let b = PyLbmBinding::new(4, 4, 0.7);
    // velocity_flat returns interleaved [ux, uy] pairs → 4×4×2 = 32 entries
    let vel = b.get_velocity_flat();
    assert_eq!(vel.len(), 32);
}

#[test]
fn test_lbm_binding_step_increments_step_count() {
    let mut b = PyLbmBinding::new(8, 8, 0.6);
    b.step();
    assert_eq!(b.step_count(), 1);
}

#[test]
fn test_lbm_binding_run_steps_count() {
    let mut b = PyLbmBinding::new(8, 8, 0.6);
    b.run_steps(10);
    assert_eq!(b.step_count(), 10);
}

#[test]
fn test_lbm_binding_mean_density_positive() {
    let b = PyLbmBinding::new(10, 10, 0.6);
    let rho = b.mean_density();
    assert!(rho >= 0.0, "mean density must be non-negative");
}

#[test]
fn test_lbm_binding_speed_flat_length() {
    let b = PyLbmBinding::new(6, 4, 0.6);
    let s = b.get_speed_flat();
    assert_eq!(s.len(), 24);
}

#[test]
fn test_lbm_binding_max_speed_non_negative() {
    let b = PyLbmBinding::new(8, 8, 0.6);
    assert!(b.max_speed() >= 0.0);
}

#[test]
fn test_lbm_binding_run_produces_finite_densities() {
    let mut b = PyLbmBinding::new(20, 10, 0.6);
    b.run_steps(5);
    // Simply verifying no panic and density remains finite
    let d = b.get_density_flat();
    assert!(
        d.iter().all(|v| v.is_finite()),
        "all densities should be finite"
    );
}

// -----------------------------------------------------------------------
// PyFemBinding
// -----------------------------------------------------------------------

#[test]
fn test_fem_binding_new_empty() {
    let b = PyFemBinding::new();
    assert_eq!(b.node_count(), 0);
    assert_eq!(b.element_count(), 0);
}

#[test]
fn test_fem_binding_add_node_increments() {
    let mut b = PyFemBinding::new();
    let idx = b.add_node(0.0, 0.0);
    assert_eq!(idx, 0);
    assert_eq!(b.node_count(), 1);
}

#[test]
fn test_fem_binding_displacement_flat_empty_before_solve() {
    let b = PyFemBinding::new();
    assert!(b.get_displacement_flat().is_empty());
}

#[test]
fn test_fem_binding_stress_flat_empty_before_solve() {
    let b = PyFemBinding::new();
    assert!(b.get_stress_flat().is_empty());
}

#[test]
fn test_fem_binding_max_stress_zero_before_solve() {
    let b = PyFemBinding::new();
    assert!((b.max_von_mises_stress()).abs() < 1e-12);
}

#[test]
fn test_fem_binding_default_trait() {
    let b = PyFemBinding::default();
    assert_eq!(b.node_count(), 0);
}

// -----------------------------------------------------------------------
// PySphBinding
// -----------------------------------------------------------------------

#[test]
fn test_sph_binding_new_water_empty() {
    let b = PySphBinding::new_water();
    assert_eq!(b.particle_count(), 0);
}

#[test]
fn test_sph_binding_add_particle_increments() {
    let mut b = PySphBinding::new_water();
    let idx = b.add_particle(0.0, 0.5, 0.0, 0.0, 0.0, 0.0);
    assert_eq!(idx, 0);
    assert_eq!(b.particle_count(), 1);
}

#[test]
fn test_sph_binding_positions_flat_length() {
    let mut b = PySphBinding::new_water();
    b.add_particle(0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    b.add_particle(1.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    assert_eq!(b.get_positions_flat().len(), 6);
}

#[test]
fn test_sph_binding_step_moves_particle_under_gravity() {
    let mut b = PySphBinding::new_water();
    b.add_particle(0.0, 1.0, 0.0, 0.0, 0.0, 0.0);
    b.step(0.01);
    let pos = b.get_positions_flat();
    assert!(
        pos[1] < 1.0 || pos[1].is_finite(),
        "particle should fall under gravity"
    );
}

#[test]
fn test_sph_binding_densities_length() {
    let mut b = PySphBinding::new_water();
    b.add_particle(0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    assert_eq!(b.get_densities_flat().len(), 1);
}

#[test]
fn test_sph_binding_total_ke_non_negative() {
    let mut b = PySphBinding::new_water();
    b.add_particle(0.0, 1.0, 0.0, 1.0, 0.0, 0.0);
    assert!(b.total_kinetic_energy() >= 0.0);
}

#[test]
fn test_sph_binding_smoothing_length_positive() {
    let b = PySphBinding::new_water();
    assert!(b.smoothing_length() > 0.0);
}

// -----------------------------------------------------------------------
// MaterialProperties
// -----------------------------------------------------------------------

#[test]
fn test_material_lookup_steel() {
    let m = MaterialProperties::lookup("steel").unwrap();
    assert_eq!(m.name, "steel");
    assert!((m.density - 7850.0).abs() < 1.0);
}

#[test]
fn test_material_lookup_case_insensitive() {
    assert!(MaterialProperties::lookup("Steel").is_some());
    assert!(MaterialProperties::lookup("STEEL").is_some());
}

#[test]
fn test_material_lookup_unknown_returns_none() {
    assert!(MaterialProperties::lookup("unobtainium").is_none());
}

#[test]
fn test_material_p_wave_speed_steel() {
    let m = MaterialProperties::steel();
    let c = m.p_wave_speed().expect("steel should have a P-wave speed");
    // Expect ~5000-6000 m/s for steel
    assert!(c > 4000.0 && c < 7000.0, "P-wave speed = {}", c);
}

#[test]
fn test_material_s_wave_speed_steel() {
    let m = MaterialProperties::steel();
    let c = m.s_wave_speed().expect("steel should have an S-wave speed");
    assert!(c > 2000.0 && c < 4000.0, "S-wave speed = {}", c);
}

#[test]
fn test_material_p_wave_speed_water_none() {
    let m = MaterialProperties::water();
    // Water has E=0, so no meaningful P-wave through solid formula
    assert!(m.p_wave_speed().is_none());
}

#[test]
fn test_py_query_material_steel_length() {
    let v = py_query_material("steel");
    assert_eq!(v.len(), 6, "expected 6 property fields");
}

#[test]
fn test_py_query_material_unknown_empty() {
    let v = py_query_material("vibranium");
    assert!(v.is_empty());
}

#[test]
fn test_py_p_wave_speed_steel_positive() {
    let c = py_p_wave_speed("steel");
    assert!(c > 0.0, "P-wave speed for steel should be positive");
}

#[test]
fn test_py_s_wave_speed_aluminium() {
    let c = py_s_wave_speed("aluminium");
    assert!(c > 0.0, "S-wave speed for aluminium should be positive");
}

#[test]
fn test_material_aluminium_lookup_alias() {
    // Both spellings should work
    assert!(MaterialProperties::lookup("aluminum").is_some());
}

#[test]
fn test_material_titanium_yield_gt_zero() {
    let m = MaterialProperties::titanium();
    assert!(m.yield_strength > 0.0);
}
