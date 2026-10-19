// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0
//! Regression tests for `deformable_coupling::CoupledSimulation::step`.
//!
//! The step used to apply the penalty coupling force to the rigid body along
//! +z and to the nodes along +x (no action-reaction), never integrated the
//! node displacements and velocities, and debonded by comparing the *total*
//! force on the rigid body - gravity and the external load included - with the
//! bond strength.

use oxiphysics_rigid::deformable_coupling::{CoupledSimulation, DeformationState};

fn no_gravity() -> CoupledSimulation {
    let mut sim = CoupledSimulation::new_default();
    sim.gravity = [0.0; 3];
    sim
}

fn norm(v: [f64; 3]) -> f64 {
    (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
}

fn total_energy(sim: &CoupledSimulation) -> f64 {
    sim.rigid_kinetic_energy() + sim.nodal_kinetic_energy() + sim.coupling_strain_energy()
}

// ── action and reaction ─────────────────────────────────────────────────────

#[test]
fn internal_forces_conserve_the_total_momentum() {
    let mut sim = no_gravity();
    sim.rigid_velocity = [3.0, -2.0, 1.0];
    sim.rigid_position = [0.01, 0.02, -0.01];
    sim.deform.velocities[3] = [0.0, 5.0, 0.0];
    sim.deform.velocities[9] = [-1.0, 0.0, 2.0];
    sim.deform.displacements[5] = [0.0, 0.003, 0.0];
    sim.thermal.set_temperature(343.15);
    let p0 = sim.total_momentum();
    for _ in 0..200 {
        sim.step(0.01, [0.0; 3]);
    }
    assert!(sim.bond.active, "the run must not end in a debond");
    let p1 = sim.total_momentum();
    let scale = norm(p0).max(1.0);
    for k in 0..3 {
        assert!(
            (p1[k] - p0[k]).abs() <= 1e-12 * scale,
            "momentum component {k}: {} -> {}",
            p0[k],
            p1[k]
        );
    }
    // The state really changed: the rigid body and the nodes exchanged momentum.
    assert!((sim.rigid_velocity[0] - 3.0).abs() > 1e-3);
}

#[test]
fn coupling_force_acts_along_the_rigid_displacement_with_equal_reaction() {
    // Rigid body displaced 1 cm along +y from its rest position: every node is
    // pulled along +y with k * 0.01, the rigid body is pulled back by the sum.
    let mut sim = no_gravity();
    sim.rigid_position = [0.0, 0.01, 0.0];
    let dt = 1e-7;
    sim.step(dt, [0.0; 3]);
    let per_node = sim.bond.penalty_stiffness / 16.0 * 0.01;
    for (i, f) in sim.deform.nodal_forces.iter().enumerate() {
        assert!(f[0].abs() < 1e-6 * per_node, "node {i}: fx = {}", f[0]);
        assert!(f[2].abs() < 1e-6 * per_node, "node {i}: fz = {}", f[2]);
        assert!(
            (f[1] - per_node).abs() < 1e-3 * per_node,
            "node {i}: fy = {} vs {per_node}",
            f[1]
        );
    }
    let total = 16.0 * per_node;
    let expected_dv = -total * dt / sim.mass;
    assert!(
        (sim.rigid_velocity[1] - expected_dv).abs() < 1e-3 * expected_dv.abs(),
        "rigid dvy = {} vs {expected_dv}",
        sim.rigid_velocity[1]
    );
    assert!(sim.rigid_velocity[0].abs() < 1e-12 && sim.rigid_velocity[2].abs() < 1e-12);
    // The nodes are accelerated the opposite way.
    for v in &sim.deform.velocities {
        assert!(
            v[1] > 0.0 && v[0].abs() < 1e-12 && v[2].abs() < 1e-12,
            "{v:?}"
        );
    }
}

#[test]
fn relaxed_springs_exert_no_force() {
    let mut sim = no_gravity();
    sim.rigid_position = [1.0, 2.0, 3.0];
    sim.rigid_rest_position = [1.0, 2.0, 3.0];
    sim.step(0.01, [0.0; 3]);
    assert_eq!(sim.rigid_position, [1.0, 2.0, 3.0]);
    assert!(sim.deform.nodal_forces.iter().all(|f| norm(*f) < 1e-12));
    assert!(sim.deform.velocities.iter().all(|v| norm(*v) < 1e-15));
}

#[test]
fn nodes_follow_the_rigid_body_and_share_its_momentum() {
    // The node displacements and velocities used to stay at zero forever.
    let mut sim = no_gravity();
    sim.rigid_velocity = [1.0, 0.0, 0.0];
    for _ in 0..3000 {
        sim.step(1e-3, [0.0; 3]);
    }
    // Perfectly inelastic exchange: v = M v0 / (M + 16 m_node).
    let common = 100.0 / (100.0 + 16.0 * 0.5);
    assert!(
        (sim.rigid_velocity[0] - common).abs() < 1e-6,
        "rigid vx = {}",
        sim.rigid_velocity[0]
    );
    for (i, v) in sim.deform.velocities.iter().enumerate() {
        assert!((v[0] - common).abs() < 1e-6, "node {i}: vx = {}", v[0]);
    }
    for (i, d) in sim.deform.displacements.iter().enumerate() {
        assert!(d[0] > 0.0, "node {i} never moved");
        assert!((d[0] - sim.rigid_position[0]).abs() < 1e-6, "node {i}");
    }
    assert!(sim.deform.nodal_forces.iter().all(|f| norm(*f) < 1e-3));
}

#[test]
fn external_force_changes_the_momentum_by_exactly_force_times_time() {
    let mut sim = no_gravity();
    let steps = 100;
    let dt = 0.01;
    for _ in 0..steps {
        sim.step(dt, [1.0e3, 0.0, 0.0]);
    }
    assert!(sim.bond.active);
    let p = sim.total_momentum();
    let expected = 1.0e3 * steps as f64 * dt;
    assert!((p[0] - expected).abs() < 1e-9 * expected, "px = {}", p[0]);
    assert!(p[1].abs() < 1e-9 && p[2].abs() < 1e-9);
}

#[test]
fn gravity_acts_on_the_nodes_as_well_as_on_the_rigid_body() {
    let mut sim = CoupledSimulation::new_default();
    let dt = 0.01;
    let steps = 100;
    for _ in 0..steps {
        sim.step(dt, [0.0; 3]);
    }
    let total_mass = 100.0 + 16.0 * 0.5;
    let p = sim.total_momentum();
    let expected = -9.81 * total_mass * steps as f64 * dt;
    assert!(
        (p[1] - expected).abs() < 1e-9 * expected.abs(),
        "py = {}",
        p[1]
    );
    // Free fall of the whole system: everything has the same velocity and
    // nothing stretches.
    for v in &sim.deform.velocities {
        assert!((v[1] - sim.rigid_velocity[1]).abs() < 1e-9);
    }
    assert!(sim.coupling_strain_energy() < 1e-12);
}

#[test]
fn fluid_pressure_loads_the_rigid_body_only() {
    let mut sim = no_gravity();
    sim.fsi.fluid_pressure = 1.0e3;
    sim.fsi.surface_area = 0.25;
    sim.fsi.surface_normal = [0.0, 1.0, 0.0];
    let steps = 50;
    for _ in 0..steps {
        sim.step(0.01, [0.0; 3]);
    }
    let p = sim.total_momentum();
    let expected = 1.0e3 * 0.25 * steps as f64 * 0.01;
    assert!((p[1] - expected).abs() < 1e-9 * expected, "py = {}", p[1]);
    assert!(p[0].abs() < 1e-9 && p[2].abs() < 1e-9);
}

// ── thermal expansion ───────────────────────────────────────────────────────

#[test]
fn thermal_expansion_is_an_internal_force_pair_along_the_thermal_axis() {
    let mut sim = no_gravity();
    sim.thermal.set_temperature(393.15);
    sim.thermal_axis = [0.0, 0.0, 2.0];
    let force = sim.thermal.thermal_force();
    assert!(force > 1.0e4);
    for _ in 0..3000 {
        sim.step(1e-3, [0.0; 3]);
    }
    assert!(sim.bond.active);
    // The centre of mass of the system does not move: no net force.
    let p = sim.total_momentum();
    assert!(norm(p) < 1e-9 * force, "momentum {p:?}");
    // Static equilibrium: the nodes sit F / K_bond ahead of the rigid body
    // along +z and no force is left on them.
    let offset = force / sim.bond.penalty_stiffness;
    let rigid_z = sim.rigid_position[2] - sim.rigid_rest_position[2];
    for (i, d) in sim.deform.displacements.iter().enumerate() {
        assert!(d[0].abs() < 1e-9 && d[1].abs() < 1e-9, "node {i}: {d:?}");
        assert!(
            ((d[2] - rigid_z) - offset).abs() < 1e-6 * offset,
            "node {i}: relative offset {} vs {offset}",
            d[2] - rigid_z
        );
    }
    assert!(sim.rigid_position[2] < 0.0, "the rigid body is pushed back");
    assert!(sim.deform.nodal_forces.iter().all(|f| norm(*f) < 1e-3));
    assert!(sim.rigid_position[0].abs() < 1e-12 && sim.rigid_position[1].abs() < 1e-12);
}

#[test]
fn thermal_expansion_at_the_reference_temperature_is_inert() {
    let mut sim = no_gravity();
    for _ in 0..10 {
        sim.step(0.01, [0.0; 3]);
    }
    assert_eq!(sim.rigid_position, [0.0; 3]);
    assert!(sim.deform.displacements.iter().all(|d| norm(*d) == 0.0));
}

// ── debonding ───────────────────────────────────────────────────────────────

#[test]
fn free_fall_and_a_light_load_do_not_debond_a_weak_bond() {
    // The bond carries nothing while the whole system falls freely, so even a
    // 10 N bond holds although the weight is ~1 kN.
    let mut sim = CoupledSimulation::new_default();
    sim.bond.max_force = 10.0;
    for _ in 0..100 {
        sim.step(0.01, [0.0; 3]);
    }
    assert!(sim.bond.active, "free fall must not debond the bond");
    let expected_y = -9.81 * 0.01 * 0.01 * (100.0 * 101.0 / 2.0);
    assert!(
        (sim.rigid_position[1] - expected_y).abs() < 1e-9,
        "y = {} vs {expected_y}",
        sim.rigid_position[1]
    );
    assert!(sim.deform.nodal_forces.iter().all(|f| norm(*f) < 1e-6));

    // A 1 kN pull needs only m_node * a = 4.6 N per node.
    let mut pulled = no_gravity();
    pulled.bond.max_force = 100.0;
    for _ in 0..20 {
        pulled.step(0.01, [1.0e3, 0.0, 0.0]);
    }
    assert!(pulled.bond.active);
}

#[test]
fn bond_debonds_when_a_nodal_coupling_force_exceeds_its_strength() {
    // 100 kN on the rigid body: each node must carry m_node * a = 463 N > 100 N.
    let mut sim = no_gravity();
    sim.bond.max_force = 100.0;
    sim.step(0.01, [1.0e5, 0.0, 0.0]);
    assert!(!sim.bond.active, "the overload must debond the bond");
    // A debonded system is frozen.
    let (position, velocity) = (sim.rigid_position, sim.rigid_velocity);
    sim.step(0.01, [1.0e5, 0.0, 0.0]);
    assert_eq!(sim.rigid_position, position);
    assert_eq!(sim.rigid_velocity, velocity);
}

#[test]
fn unbreakable_bond_never_debonds() {
    let mut sim = no_gravity();
    sim.bond.max_force = f64::INFINITY;
    sim.step(0.01, [1.0e9, 0.0, 0.0]);
    assert!(sim.bond.active);
}

// ── stability ───────────────────────────────────────────────────────────────

#[test]
fn stiff_bond_at_a_large_step_does_not_gain_energy() {
    // Node frequency sqrt(k / m) = 2.5e5 rad/s against a 50 ms step: an
    // explicit integrator diverges immediately.
    let mut sim = no_gravity();
    sim.bond.penalty_stiffness = 1.0e10;
    sim.bond.penalty_damping = 1.0e3;
    sim.bond.max_force = f64::INFINITY;
    sim.node_mass = 0.01;
    sim.rigid_position = [0.0, 0.0, 0.05];
    sim.rigid_velocity = [0.5, 0.0, 0.0];
    let mut energy = total_energy(&sim);
    assert!(energy > 0.0);
    let initial = energy;
    for step in 0..200 {
        sim.step(0.05, [0.0; 3]);
        let next = total_energy(&sim);
        assert!(next.is_finite(), "step {step}: energy {next}");
        assert!(
            next <= energy * (1.0 + 1e-9),
            "step {step}: energy grew from {energy} to {next}"
        );
        energy = next;
    }
    assert!(energy < initial);
}

#[test]
fn default_bond_is_stable_at_the_default_step() {
    // k = 6.25e4 N/m per node against a 0.5 kg node: omega * dt = 3.5 for
    // dt = 0.01, beyond the limit of explicit semi-implicit Euler.
    let mut sim = no_gravity();
    sim.rigid_position = [0.0, 0.0, 0.02];
    let initial = total_energy(&sim);
    for _ in 0..500 {
        sim.step(0.01, [0.0; 3]);
        assert!(total_energy(&sim) <= initial * (1.0 + 1e-9));
    }
    assert!(sim.rigid_position.iter().all(|c| c.is_finite()));
    assert!(sim.rigid_position[2].abs() <= 0.02 + 1e-12);
}

// ── robustness ──────────────────────────────────────────────────────────────

#[test]
fn rigid_body_without_interface_nodes_is_integrated_on_its_own() {
    let mut sim = CoupledSimulation::new_default();
    sim.deform = DeformationState::new_undeformed(0);
    for _ in 0..10 {
        sim.step(0.1, [0.0; 3]);
    }
    assert!((sim.rigid_velocity[1] + 9.81).abs() < 1e-12);
    assert_eq!(
        sim.total_momentum(),
        [0.0, 100.0 * sim.rigid_velocity[1], 0.0]
    );
    assert_eq!(sim.coupling_strain_energy(), 0.0);
}

#[test]
fn inconsistent_deformation_arrays_are_repaired_instead_of_panicking() {
    let mut sim = no_gravity();
    sim.deform.velocities.clear();
    sim.deform.nodal_forces.truncate(3);
    sim.rigid_position = [0.01, 0.0, 0.0];
    sim.step(0.01, [0.0; 3]);
    assert_eq!(sim.deform.velocities.len(), 16);
    assert_eq!(sim.deform.nodal_forces.len(), 16);
    assert!(sim.deform.velocities.iter().all(|v| v[0] > 0.0));
}

#[test]
fn non_positive_and_nan_time_steps_are_ignored() {
    for dt in [0.0, -0.01, f64::NAN] {
        let mut sim = CoupledSimulation::new_default();
        sim.step(dt, [1.0, 2.0, 3.0]);
        assert_eq!(sim.rigid_position, [0.0; 3], "dt {dt}");
        assert_eq!(sim.rigid_velocity, [0.0; 3], "dt {dt}");
        assert!(sim.bond.active);
    }
}

#[test]
fn coupling_strain_energy_is_half_k_times_the_stretch_squared() {
    let mut sim = no_gravity();
    sim.rigid_position = [0.01, 0.0, 0.0];
    // 16 nodes, k_total = 1e6 N/m: 1/2 * 1e6 * 0.01^2 = 50 J.
    assert!((sim.coupling_strain_energy() - 50.0).abs() < 1e-9);
    sim.deform.displacements[0] = [0.01, 0.0, 0.0];
    let k_node = 1.0e6 / 16.0;
    assert!((sim.coupling_strain_energy() - (50.0 - 0.5 * k_node * 1e-4)).abs() < 1e-9);
}
