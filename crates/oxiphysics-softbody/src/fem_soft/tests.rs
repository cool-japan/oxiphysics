// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Tests for FEM soft body modules.

#![allow(clippy::needless_range_loop)]

use super::math_helpers::*;
use super::*;

fn make_material() -> SoftMaterial {
    SoftMaterial {
        youngs_modulus: 1000.0,
        poisson_ratio: 0.3,
        damping: 0.01,
    }
}

fn tet_positions() -> Vec<[f64; 3]> {
    vec![
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
    ]
}

#[test]
fn test_raw_element_rest_volume() {
    let pos = tet_positions();
    let elem = CorotationalElementRaw::new([0, 1, 2, 3], &pos, make_material());
    // Volume of unit tet = 1/6
    assert!((elem.rest_volume - 1.0 / 6.0).abs() < 1e-12);
}

#[test]
fn test_raw_element_zero_force_at_rest() {
    let pos = tet_positions();
    let elem = CorotationalElementRaw::new([0, 1, 2, 3], &pos, make_material());
    let forces = elem.compute_forces(&pos);
    for k in 0..4 {
        let mag = (forces[k][0] * forces[k][0]
            + forces[k][1] * forces[k][1]
            + forces[k][2] * forces[k][2])
            .sqrt();
        assert!(
            mag < 1e-8,
            "force on node {k} should be zero at rest, got {mag}"
        );
    }
}

#[test]
fn test_raw_element_nonzero_force_when_stretched() {
    let pos = tet_positions();
    let elem = CorotationalElementRaw::new([0, 1, 2, 3], &pos, make_material());
    let mut stretched = pos.clone();
    stretched[1] = [2.0, 0.0, 0.0]; // stretch edge 0-1
    let forces = elem.compute_forces(&stretched);
    let mag =
        (forces[1][0] * forces[1][0] + forces[1][1] * forces[1][1] + forces[1][2] * forces[1][2])
            .sqrt();
    assert!(mag > 1e-6, "stretching should produce nonzero force");
}

#[test]
fn test_raw_fem_body_step_moves_nodes() {
    let pos = tet_positions();
    let elem = CorotationalElementRaw::new([0, 1, 2, 3], &pos, make_material());
    let mut body = FemSoftBodyRaw::new(pos.clone(), vec![elem], 1.0);
    body.step(0.01, [0.0, -9.81, 0.0]);
    // Nodes should have moved under gravity
    let moved = body.nodes.iter().zip(pos.iter()).any(|(a, b)| {
        (a[0] - b[0]).abs() > 1e-10 || (a[1] - b[1]).abs() > 1e-10 || (a[2] - b[2]).abs() > 1e-10
    });
    assert!(moved, "nodes should have moved after step");
}

#[test]
fn test_raw_fem_kinetic_energy_zero_initially() {
    let pos = tet_positions();
    let elem = CorotationalElementRaw::new([0, 1, 2, 3], &pos, make_material());
    let body = FemSoftBodyRaw::new(pos, vec![elem], 1.0);
    assert!(body.kinetic_energy().abs() < 1e-15);
}

#[test]
fn test_raw_fem_kinetic_energy_positive_after_step() {
    let pos = tet_positions();
    let elem = CorotationalElementRaw::new([0, 1, 2, 3], &pos, make_material());
    let mut body = FemSoftBodyRaw::new(pos, vec![elem], 1.0);
    body.step(0.01, [0.0, -9.81, 0.0]);
    assert!(body.kinetic_energy() > 0.0);
}

#[test]
fn test_raw_fem_potential_energy() {
    let pos = tet_positions();
    let elem = CorotationalElementRaw::new([0, 1, 2, 3], &pos, make_material());
    let body = FemSoftBodyRaw::new(pos, vec![elem], 1.0);
    let pe = body.potential_energy([0.0, -9.81, 0.0]);
    // Nodes at y=0,0,1,0 => dot with (0,-9.81,0) sum = -1*(-9.81) = 9.81
    assert!(pe > 0.0, "PE should be positive for nodes above reference");
}

#[test]
fn test_compute_rotation_identity() {
    let pos = tet_positions();
    let r = CorotationalElementRaw::compute_rotation(&pos, &[0, 1, 2, 3]);
    // For the identity tet, rotation should be close to identity
    for i in 0..3 {
        for j in 0..3 {
            let expected = if i == j { 1.0 } else { 0.0 };
            assert!(
                (r[i][j] - expected).abs() < 0.1,
                "R[{i}][{j}] = {} expected {expected}",
                r[i][j]
            );
        }
    }
}

#[test]
fn test_soft_material_struct() {
    let m = SoftMaterial {
        youngs_modulus: 500.0,
        poisson_ratio: 0.4,
        damping: 0.05,
    };
    assert!((m.youngs_modulus - 500.0).abs() < 1e-10);
    assert!((m.poisson_ratio - 0.4).abs() < 1e-10);
    assert!((m.damping - 0.05).abs() < 1e-10);
}

// ── Neo-Hookean material tests ──────────────────────────────────────

#[test]
fn test_neohookean_from_young_poisson() {
    let mat = NeoHookeanMaterial::from_young_poisson(1000.0, 0.3);
    // μ = E / (2(1+ν)) = 1000 / 2.6 ≈ 384.6
    let expected_mu = 1000.0 / (2.0 * 1.3);
    assert!((mat.mu - expected_mu).abs() < 1e-6);
}

#[test]
fn test_neohookean_identity_deformation() {
    let mat = NeoHookeanMaterial::from_young_poisson(1000.0, 0.3);
    let identity = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    let p = mat.piola_kirchhoff(identity);
    // At identity deformation, stress should be zero
    for i in 0..3 {
        for j in 0..3 {
            assert!(
                p[i][j].abs() < 1e-10,
                "P[{i}][{j}] should be 0 at identity, got {}",
                p[i][j]
            );
        }
    }
}

#[test]
fn test_neohookean_strain_energy_zero_at_rest() {
    let mat = NeoHookeanMaterial::from_young_poisson(1000.0, 0.3);
    let identity = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    let w = mat.strain_energy(identity);
    assert!(
        w.abs() < 1e-10,
        "Strain energy at rest should be 0, got {w}"
    );
}

#[test]
fn test_neohookean_strain_energy_positive_when_stretched() {
    let mat = NeoHookeanMaterial::from_young_poisson(1000.0, 0.3);
    let stretched = [[2.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    let w = mat.strain_energy(stretched);
    assert!(
        w > 0.0,
        "Strain energy should be positive when stretched, got {w}"
    );
}

#[test]
fn test_neohookean_element_rest_forces_zero() {
    let pos = tet_positions();
    let mat = NeoHookeanMaterial::from_young_poisson(1000.0, 0.3);
    let elem = NeoHookeanElement::new([0, 1, 2, 3], &pos, mat);
    let forces = elem.compute_forces(&pos);
    for k in 0..4 {
        let mag = (forces[k][0] * forces[k][0]
            + forces[k][1] * forces[k][1]
            + forces[k][2] * forces[k][2])
            .sqrt();
        assert!(
            mag < 1e-8,
            "Force on node {k} at rest should be zero, got {mag}"
        );
    }
}

#[test]
fn test_neohookean_element_stretched_forces() {
    let pos = tet_positions();
    let mat = NeoHookeanMaterial::from_young_poisson(1000.0, 0.3);
    let elem = NeoHookeanElement::new([0, 1, 2, 3], &pos, mat);
    let mut stretched = pos.clone();
    stretched[1] = [2.0, 0.0, 0.0];
    let forces = elem.compute_forces(&stretched);
    let mag =
        (forces[1][0] * forces[1][0] + forces[1][1] * forces[1][1] + forces[1][2] * forces[1][2])
            .sqrt();
    assert!(mag > 1e-6, "Stretched element should have nonzero force");
}

#[test]
fn test_neohookean_element_strain_energy() {
    let pos = tet_positions();
    let mat = NeoHookeanMaterial::from_young_poisson(1000.0, 0.3);
    let elem = NeoHookeanElement::new([0, 1, 2, 3], &pos, mat);
    let energy_rest = elem.strain_energy(&pos);
    assert!(
        energy_rest.abs() < 1e-8,
        "Strain energy at rest should be ~0"
    );

    let mut stretched = pos.clone();
    stretched[1] = [2.0, 0.0, 0.0];
    let energy_stretched = elem.strain_energy(&stretched);
    assert!(
        energy_stretched > 0.0,
        "Strain energy when stretched should be > 0"
    );
}

#[test]
fn test_hyperelastic_body_step() {
    let pos = tet_positions();
    let mat = NeoHookeanMaterial::from_young_poisson(1000.0, 0.3);
    let elem = NeoHookeanElement::new([0, 1, 2, 3], &pos, mat);
    let mut body = HyperelasticBody::new(pos.clone(), vec![elem], 1.0, 0.01);
    body.step(0.01, [0.0, -9.81, 0.0]);
    let moved = body.nodes.iter().zip(pos.iter()).any(|(a, b)| {
        (a[0] - b[0]).abs() > 1e-10 || (a[1] - b[1]).abs() > 1e-10 || (a[2] - b[2]).abs() > 1e-10
    });
    assert!(moved, "Nodes should move after step");
}

#[test]
fn test_hyperelastic_body_kinetic_energy() {
    let pos = tet_positions();
    let mat = NeoHookeanMaterial::from_young_poisson(1000.0, 0.3);
    let elem = NeoHookeanElement::new([0, 1, 2, 3], &pos, mat);
    let body = HyperelasticBody::new(pos, vec![elem], 1.0, 0.01);
    assert!(
        body.kinetic_energy().abs() < 1e-15,
        "Initial KE should be 0"
    );
}

#[test]
fn test_hyperelastic_body_ke_after_step() {
    let pos = tet_positions();
    let mat = NeoHookeanMaterial::from_young_poisson(1000.0, 0.3);
    let elem = NeoHookeanElement::new([0, 1, 2, 3], &pos, mat);
    let mut body = HyperelasticBody::new(pos, vec![elem], 1.0, 0.01);
    body.step(0.01, [0.0, -9.81, 0.0]);
    assert!(
        body.kinetic_energy() > 0.0,
        "KE should be positive after gravity step"
    );
}

// ── Stress computation tests ────────────────────────────────────────

#[test]
fn test_von_mises_zero_stress() {
    let sigma = [[0.0; 3]; 3];
    let vm = von_mises_stress(sigma);
    assert!(vm.abs() < 1e-15);
}

#[test]
fn test_von_mises_uniaxial() {
    // Uniaxial stress σ11 = 100, rest zero
    let mut sigma = [[0.0; 3]; 3];
    sigma[0][0] = 100.0;
    let vm = von_mises_stress(sigma);
    // For uniaxial: σ_vm = σ11 = 100
    assert!(
        (vm - 100.0).abs() < 1e-6,
        "Von Mises should be 100, got {vm}"
    );
}

#[test]
fn test_von_mises_hydrostatic() {
    // Hydrostatic: σ11 = σ22 = σ33 = p, von Mises = 0
    let p = 50.0;
    let sigma = [[p, 0.0, 0.0], [0.0, p, 0.0], [0.0, 0.0, p]];
    let vm = von_mises_stress(sigma);
    assert!(
        vm.abs() < 1e-10,
        "Hydrostatic: von Mises should be 0, got {vm}"
    );
}

#[test]
fn test_green_lagrange_identity() {
    let identity = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    let e = green_lagrange_strain(identity);
    for i in 0..3 {
        for j in 0..3 {
            assert!(e[i][j].abs() < 1e-14, "E[{i}][{j}] should be 0 at identity");
        }
    }
}

#[test]
fn test_green_lagrange_stretch() {
    let stretched = [[2.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    let e = green_lagrange_strain(stretched);
    // E11 = 0.5(4 - 1) = 1.5
    assert!((e[0][0] - 1.5).abs() < 1e-10);
    assert!(e[1][1].abs() < 1e-10);
    assert!(e[2][2].abs() < 1e-10);
}

#[test]
fn test_cauchy_stress_at_identity() {
    let identity = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    let p_zero = [[0.0; 3]; 3];
    let sigma = cauchy_stress(p_zero, identity);
    for i in 0..3 {
        for j in 0..3 {
            assert!(sigma[i][j].abs() < 1e-14);
        }
    }
}

// ── Contact handling tests ──────────────────────────────────────────

#[test]
fn test_ground_contact_stops_falling() {
    let mut nodes = vec![[0.0, -0.1, 0.0]]; // below ground
    let mut velocities = vec![[0.0, -1.0, 0.0]];
    apply_ground_contact(&mut nodes, &mut velocities, 0.0, 10000.0, 100.0, 0.01);
    assert!(nodes[0][1] >= 0.0, "Node should be at or above ground");
}

#[test]
fn test_ground_contact_no_effect_above() {
    let mut nodes = vec![[0.0, 1.0, 0.0]]; // above ground
    let mut velocities = vec![[0.0, -1.0, 0.0]];
    let vel_before = velocities[0][1];
    apply_ground_contact(&mut nodes, &mut velocities, 0.0, 10000.0, 100.0, 0.01);
    assert_eq!(
        velocities[0][1], vel_before,
        "No contact → velocity unchanged"
    );
}

// ── Matrix helper tests ─────────────────────────────────────────────

#[test]
fn test_det3x3_identity() {
    let identity = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    assert!((det3x3(identity) - 1.0).abs() < 1e-14);
}

#[test]
fn test_inv3x3_identity() {
    let identity = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    let inv = inv3x3(identity);
    for i in 0..3 {
        for j in 0..3 {
            let expected = if i == j { 1.0 } else { 0.0 };
            assert!((inv[i][j] - expected).abs() < 1e-14);
        }
    }
}

#[test]
fn test_transpose3x3() {
    let m = [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]];
    let t = transpose3x3(m);
    assert_eq!(t[0][1], 4.0);
    assert_eq!(t[1][0], 2.0);
    assert_eq!(t[2][0], 3.0);
}

#[test]
fn test_mul3x3_identity() {
    let identity = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    let m = [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]];
    let result = mul3x3(identity, m);
    for i in 0..3 {
        for j in 0..3 {
            assert!((result[i][j] - m[i][j]).abs() < 1e-14);
        }
    }
}

// ── CorotFem tests ───────────────────────────────────────────────────────

fn make_corot_nodes() -> Vec<CorotFemNode> {
    vec![
        CorotFemNode {
            position: [0.0, 0.0, 0.0],
            velocity: [0.0; 3],
            force: [0.0; 3],
            mass: 1.0,
        },
        CorotFemNode {
            position: [1.0, 0.0, 0.0],
            velocity: [0.0; 3],
            force: [0.0; 3],
            mass: 1.0,
        },
        CorotFemNode {
            position: [0.0, 1.0, 0.0],
            velocity: [0.0; 3],
            force: [0.0; 3],
            mass: 1.0,
        },
        CorotFemNode {
            position: [0.0, 0.0, 1.0],
            velocity: [0.0; 3],
            force: [0.0; 3],
            mass: 1.0,
        },
    ]
}

#[test]
fn test_corot_polar_decomp_of_rotation() {
    // Pure rotation by 90° about Z axis: [[0,-1,0],[1,0,0],[0,0,1]]
    let rot = [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]];
    let r = CorotFemTet::polar_decompose_rotation(rot);
    // Should recover the same rotation
    for i in 0..3 {
        for j in 0..3 {
            assert!(
                (r[i][j] - rot[i][j]).abs() < 1e-6,
                "R[{i}][{j}]={} expected {}",
                r[i][j],
                rot[i][j]
            );
        }
    }
}

#[test]
fn test_corot_deformation_gradient_identity_at_rest() {
    let nodes = make_corot_nodes();
    let tet = CorotFemTet::new([0, 1, 2, 3], &nodes);
    let f = tet.compute_deformation_gradient(&nodes);
    // At rest, F = I
    for i in 0..3 {
        for j in 0..3 {
            let expected = if i == j { 1.0 } else { 0.0 };
            assert!(
                (f[i][j] - expected).abs() < 1e-10,
                "F[{i}][{j}]={} expected {expected}",
                f[i][j]
            );
        }
    }
}

#[test]
fn test_corot_elastic_force_zero_at_rest() {
    let nodes = make_corot_nodes();
    let tet = CorotFemTet::new([0, 1, 2, 3], &nodes);
    let forces = tet.compute_elastic_forces(&nodes, 500.0, 300.0);
    for k in 0..4 {
        let mag = (forces[k][0] * forces[k][0]
            + forces[k][1] * forces[k][1]
            + forces[k][2] * forces[k][2])
            .sqrt();
        assert!(
            mag < 1e-8,
            "force at rest should be ~0 for node {k}, got {mag}"
        );
    }
}

#[test]
fn test_corot_elastic_force_newton_third_law() {
    // Forces should sum to zero
    let nodes = make_corot_nodes();
    let tet = CorotFemTet::new([0, 1, 2, 3], &nodes);
    // Deform node 1 away from rest
    let mut deformed = nodes.clone();
    deformed[1].position = [2.0, 0.0, 0.0];
    let forces = tet.compute_elastic_forces(&deformed, 500.0, 300.0);
    let mut sum = [0.0f64; 3];
    for k in 0..4 {
        for d in 0..3 {
            sum[d] += forces[k][d];
        }
    }
    for d in 0..3 {
        assert!(
            sum[d].abs() < 1e-8,
            "force sum[{d}]={} should be 0 (Newton 3rd law)",
            sum[d]
        );
    }
}

#[test]
fn test_corot_sim_step_moves_nodes() {
    let nodes = make_corot_nodes();
    let mut sim_nodes = nodes.clone();
    // Give node 1 an initial force via stretching
    sim_nodes[1].position = [1.5, 0.0, 0.0];
    let tet = CorotFemTet::new([0, 1, 2, 3], &nodes); // rest shape from original
    let sim = CorotFemSim::new(sim_nodes.clone(), vec![tet], 500.0, 300.0, 0.001);
    let mut sim2 = sim.clone();
    sim2.step();
    let moved = sim2.nodes.iter().zip(sim_nodes.iter()).any(|(a, b)| {
        (a.position[0] - b.position[0]).abs() > 1e-15
            || (a.position[1] - b.position[1]).abs() > 1e-15
            || (a.position[2] - b.position[2]).abs() > 1e-15
    });
    assert!(moved, "nodes should move after sim step");
}

#[test]
fn test_corot_total_strain_energy_zero_at_rest() {
    let nodes = make_corot_nodes();
    let tet = CorotFemTet::new([0, 1, 2, 3], &nodes);
    let sim = CorotFemSim::new(nodes, vec![tet], 500.0, 300.0, 0.001);
    let energy = sim.total_strain_energy();
    assert!(
        energy.abs() < 1e-8,
        "strain energy at rest should be ~0, got {energy}"
    );
}

#[test]
fn test_corot_total_strain_energy_positive_when_deformed() {
    let nodes = make_corot_nodes();
    let tet = CorotFemTet::new([0, 1, 2, 3], &nodes);
    let mut sim_nodes = nodes;
    sim_nodes[1].position = [2.0, 0.0, 0.0];
    let sim = CorotFemSim::new(sim_nodes, vec![tet], 500.0, 300.0, 0.001);
    let energy = sim.total_strain_energy();
    assert!(
        energy > 0.0,
        "strain energy when deformed should be positive, got {energy}"
    );
}

#[test]
fn test_corot_tet_alias_matches_corot_fem_tet() {
    // CorotTet should behave identically to CorotFemTet
    let nodes = make_corot_nodes();
    let tet = CorotTet::new([0, 1, 2, 3], &nodes);
    let f = tet.compute_deformation_gradient(&nodes);
    // At rest, F should be identity
    for i in 0..3 {
        for j in 0..3 {
            let expected = if i == j { 1.0 } else { 0.0 };
            assert!((f[i][j] - expected).abs() < 1e-10);
        }
    }
}

#[test]
fn test_newmark_beta_rest_no_motion() {
    // A body at rest with no forces should not move
    let nodes = make_corot_nodes();
    let tet = CorotFemTet::new([0, 1, 2, 3], &nodes);
    let mut sim = NewmarkBetaSim::new(nodes.clone(), vec![tet], 500.0, 300.0, 0.01);
    let initial_positions: Vec<[f64; 3]> = sim.nodes.iter().map(|n| n.position).collect();
    sim.step();
    for (i, n) in sim.nodes.iter().enumerate() {
        for d in 0..3 {
            assert!(
                (n.position[d] - initial_positions[i][d]).abs() < 1e-10,
                "node {i} dim {d} should not move: {} vs {}",
                n.position[d],
                initial_positions[i][d]
            );
        }
    }
}

#[test]
fn test_newmark_beta_deformed_returns_to_rest() {
    // After deforming, elastic forces should push nodes back
    let nodes = make_corot_nodes();
    let tet = CorotFemTet::new([0, 1, 2, 3], &nodes);
    let mut deformed = nodes.clone();
    deformed[1].position = [1.5, 0.0, 0.0];
    let mut sim = NewmarkBetaSim::new(deformed, vec![tet], 500.0, 300.0, 0.001);
    // Initial elastic energy
    let e0 = sim.total_strain_energy();
    assert!(e0 > 0.0, "deformed body should have positive strain energy");
    // Step forward - the body should respond
    sim.step();
    // Nodes should have moved
    let moved = sim.nodes[1].velocity.iter().any(|&v| v.abs() > 1e-15);
    assert!(moved, "nodes should have nonzero velocity after step");
}

#[test]
fn test_newmark_beta_energy_conservation_undamped() {
    // Undamped Newmark-beta (beta=0.25, gamma=0.5) should be energy-conserving
    let nodes = make_corot_nodes();
    let tet = CorotFemTet::new([0, 1, 2, 3], &nodes);
    let mut deformed = nodes.clone();
    deformed[1].position = [1.2, 0.0, 0.0];
    let mut sim = NewmarkBetaSim::new(deformed, vec![tet], 1000.0, 500.0, 0.001);
    let e0 = sim.total_mechanical_energy();
    for _ in 0..50 {
        sim.step();
    }
    let e1 = sim.total_mechanical_energy();
    // Energy should be within 10% after 50 steps (implicit Newmark is stable)
    assert!(
        (e1 - e0).abs() / (e0.abs() + 1e-30) < 0.10,
        "energy drift too large: e0={e0}, e1={e1}"
    );
}

#[test]
fn test_newmark_beta_boundary_pinned_nodes_no_move() {
    let nodes = make_corot_nodes();
    let tet = CorotFemTet::new([0, 1, 2, 3], &nodes);
    let mut deformed = nodes.clone();
    deformed[1].position = [1.5, 0.0, 0.0];
    let mut sim = NewmarkBetaSim::new(deformed, vec![tet], 500.0, 300.0, 0.01);
    // Pin node 0
    sim.pinned[0] = true;
    let pos0 = sim.nodes[0].position;
    sim.step();
    for d in 0..3 {
        assert!(
            (sim.nodes[0].position[d] - pos0[d]).abs() < 1e-14,
            "pinned node 0 should not move in dim {d}"
        );
    }
}

#[test]
fn test_corot_tet_deformation_gradient_uniform_stretch() {
    // Uniform stretch along x by factor 2: F should be diag(2,1,1)
    let mut nodes = make_corot_nodes();
    nodes[1].position = [2.0, 0.0, 0.0]; // was [1,0,0]
    // tet was constructed at rest [0,0,0],[1,0,0],[0,1,0],[0,0,1]
    let rest_nodes = make_corot_nodes();
    let tet = CorotFemTet::new([0, 1, 2, 3], &rest_nodes);
    let f = tet.compute_deformation_gradient(&nodes);
    // F[0][0] should be 2, rest identity
    assert!(
        (f[0][0] - 2.0).abs() < 1e-10,
        "F[0][0] = {} expected 2",
        f[0][0]
    );
    assert!(
        (f[1][1] - 1.0).abs() < 1e-10,
        "F[1][1] = {} expected 1",
        f[1][1]
    );
    assert!(
        (f[2][2] - 1.0).abs() < 1e-10,
        "F[2][2] = {} expected 1",
        f[2][2]
    );
}

#[test]
fn test_corot_tet_polar_decomp_pure_rotation_z90() {
    // 90° rotation about Z: [[0,-1,0],[1,0,0],[0,0,1]]
    let rot = [[0.0_f64, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]];
    let r = CorotFemTet::polar_decompose_rotation(rot);
    for i in 0..3 {
        for j in 0..3 {
            assert!(
                (r[i][j] - rot[i][j]).abs() < 1e-6,
                "R[{i}][{j}] = {} expected {}",
                r[i][j],
                rot[i][j]
            );
        }
    }
}

#[test]
fn test_newmark_beta_total_strain_energy_zero_at_rest() {
    let nodes = make_corot_nodes();
    let tet = CorotFemTet::new([0, 1, 2, 3], &nodes);
    let sim = NewmarkBetaSim::new(nodes, vec![tet], 500.0, 300.0, 0.001);
    assert!(
        sim.total_strain_energy().abs() < 1e-8,
        "strain energy at rest should be 0"
    );
}

// ── CorotationalTet tests ────────────────────────────────────────────────

fn unit_tet_positions() -> Vec<[f64; 3]> {
    vec![
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
    ]
}

#[test]
#[allow(non_snake_case)]
fn test_corot_tet_new_fields() {
    let pos = unit_tet_positions();
    let tet = CorotationalTet::new([0, 1, 2, 3], &pos);
    assert!(
        (tet.V0 - 1.0 / 6.0).abs() < 1e-12,
        "rest volume should be 1/6, got {}",
        tet.V0
    );
    // Dm columns should be the edge vectors
    let e1 = [1.0, 0.0, 0.0];
    for d in 0..3 {
        assert!((tet.Dm[d][0] - e1[d]).abs() < 1e-12, "Dm[:,0] should be e1");
    }
}

#[test]
fn test_corot_tet_deformation_gradient_identity() {
    let pos = unit_tet_positions();
    let tet = CorotationalTet::new([0, 1, 2, 3], &pos);
    let f = tet.compute_deformation_gradient(&pos);
    for i in 0..3 {
        for j in 0..3 {
            let expected = if i == j { 1.0 } else { 0.0 };
            assert!(
                (f[i][j] - expected).abs() < 1e-10,
                "F[{i}][{j}] should be I at rest, got {}",
                f[i][j]
            );
        }
    }
}

#[test]
fn test_corot_tet_polar_decompose_identity() {
    let identity = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    let (r, s) = CorotationalTet::polar_decompose(&identity);
    // R should be identity, S should be identity
    for i in 0..3 {
        for j in 0..3 {
            let expected = if i == j { 1.0 } else { 0.0 };
            assert!((r[i][j] - expected).abs() < 1e-6, "R[{i}][{j}] should be I");
            assert!((s[i][j] - expected).abs() < 1e-6, "S[{i}][{j}] should be I");
        }
    }
}

#[test]
fn test_corot_tet_elastic_forces_zero_at_rest() {
    let pos = unit_tet_positions();
    let tet = CorotationalTet::new([0, 1, 2, 3], &pos);
    let forces = tet.elastic_forces(&pos, 384.6, 576.9);
    for k in 0..4 {
        let mag = (forces[k][0] * forces[k][0]
            + forces[k][1] * forces[k][1]
            + forces[k][2] * forces[k][2])
            .sqrt();
        assert!(
            mag < 1e-8,
            "elastic force at rest should be ~0 for node {k}, got {mag}"
        );
    }
}

#[test]
fn test_corot_tet_elastic_forces_nonzero_when_stretched() {
    let pos = unit_tet_positions();
    let tet = CorotationalTet::new([0, 1, 2, 3], &pos);
    let mut stretched = pos.clone();
    stretched[1] = [2.0, 0.0, 0.0];
    let forces = tet.elastic_forces(&stretched, 384.6, 576.9);
    let mag1 =
        (forces[1][0] * forces[1][0] + forces[1][1] * forces[1][1] + forces[1][2] * forces[1][2])
            .sqrt();
    assert!(
        mag1 > 1e-6,
        "stretched tet should have nonzero force, got {mag1}"
    );
}

#[test]
fn test_corot_tet_elastic_forces_newton_third_law() {
    let pos = unit_tet_positions();
    let tet = CorotationalTet::new([0, 1, 2, 3], &pos);
    let mut stretched = pos.clone();
    stretched[1] = [1.5, 0.0, 0.0];
    let forces = tet.elastic_forces(&stretched, 384.6, 576.9);
    let mut sum = [0.0f64; 3];
    for k in 0..4 {
        for d in 0..3 {
            sum[d] += forces[k][d];
        }
    }
    for d in 0..3 {
        assert!(
            sum[d].abs() < 1e-8,
            "force sum[{d}]={} should be 0 (Newton 3rd law)",
            sum[d]
        );
    }
}

#[test]
fn test_corot_tet_stiffness_matrix_size() {
    let pos = unit_tet_positions();
    let tet = CorotationalTet::new([0, 1, 2, 3], &pos);
    let k = tet.stiffness_matrix(&pos, 384.6, 576.9);
    assert_eq!(
        k.len(),
        144,
        "12×12 stiffness matrix should have 144 entries"
    );
}

#[test]
fn test_corot_tet_stiffness_matrix_symmetry() {
    let pos = unit_tet_positions();
    let tet = CorotationalTet::new([0, 1, 2, 3], &pos);
    let mut stretched = pos.clone();
    stretched[1] = [1.3, 0.0, 0.0];
    let k = tet.stiffness_matrix(&stretched, 384.6, 576.9);
    // Stiffness matrix should be approximately symmetric
    for r in 0..12 {
        for c in 0..12 {
            let diff = (k[r * 12 + c] - k[c * 12 + r]).abs();
            assert!(
                diff < 1e-3,
                "K[{r},{c}]={} vs K[{c},{r}]={}: not symmetric (diff={diff})",
                k[r * 12 + c],
                k[c * 12 + r]
            );
        }
    }
}

// ── HighLevelFemBody tests ───────────────────────────────────────────────

fn make_high_level_body() -> HighLevelFemBody {
    let pos = unit_tet_positions();
    let masses = vec![1.0; 4];
    let mut body = HighLevelFemBody::new(pos, masses, 1000.0, 0.3);
    body.add_tet(0, 1, 2, 3);
    body
}

#[test]
fn test_high_level_fem_energy_zero_at_rest() {
    let body = make_high_level_body();
    let e = body.total_elastic_energy();
    assert!(
        e.abs() < 1e-8,
        "elastic energy at rest should be ~0, got {e}"
    );
}

#[test]
fn test_high_level_fem_energy_positive_when_deformed() {
    let mut body = make_high_level_body();
    body.positions[1] = [2.0, 0.0, 0.0];
    let e = body.total_elastic_energy();
    assert!(
        e > 0.0,
        "elastic energy when deformed should be > 0, got {e}"
    );
}

#[test]
fn test_high_level_fem_step_moves_nodes() {
    let mut body = make_high_level_body();
    body.positions[1] = [1.5, 0.0, 0.0];
    let init = body.positions.clone();
    body.step(0.001);
    let moved = body
        .positions
        .iter()
        .zip(init.iter())
        .any(|(a, b)| (0..3).any(|d| (a[d] - b[d]).abs() > 1e-15));
    assert!(moved, "nodes should move after step");
}

#[test]
fn test_high_level_fem_add_tet() {
    let pos = unit_tet_positions();
    let masses = vec![1.0; 4];
    let mut body = HighLevelFemBody::new(pos, masses, 1000.0, 0.3);
    assert_eq!(body.elements.len(), 0);
    body.add_tet(0, 1, 2, 3);
    assert_eq!(body.elements.len(), 1);
}

#[cfg(test)]
mod fem_verlet_tests {
    use super::BoundaryConditions;

    use crate::fem_soft::CorotFemElement4;
    use crate::fem_soft::CorotFemTet;
    use crate::fem_soft::FemSoftBodyVerlet;
    use crate::fem_soft::LumpedMassMatrix;
    use crate::fem_soft::RayleighDamping;
    use crate::fem_soft::VelocityVerletIntegrator;
    use crate::fem_soft::assemble_internal_forces;

    fn unit_tet_pos() -> Vec<[f64; 3]> {
        vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ]
    }

    // ── LumpedMassMatrix ─────────────────────────────────────────────────

    #[test]
    fn test_lumped_mass_total() {
        let density = 1000.0;
        let elements = vec![([0usize, 1, 2, 3], 1.0 / 6.0)];
        let lm = LumpedMassMatrix::from_elements(4, &elements, density);
        let expected = density / 6.0;
        assert!(
            (lm.total_mass() - expected).abs() < 1e-10,
            "total lumped mass = {}, expected {expected}",
            lm.total_mass()
        );
    }

    #[test]
    fn test_lumped_mass_equal_distribution() {
        let elements = vec![([0usize, 1, 2, 3], 1.0 / 6.0)];
        let lm = LumpedMassMatrix::from_elements(4, &elements, 1200.0);
        // Each node should get 1/4 of total mass
        let m_per_node = 1200.0 / 6.0 * 0.25;
        for i in 0..4 {
            assert!(
                (lm.masses[i] - m_per_node).abs() < 1e-10,
                "node {i} mass = {}, expected {m_per_node}",
                lm.masses[i]
            );
        }
    }

    #[test]
    fn test_lumped_mass_inv_mass() {
        let elements = vec![([0usize, 1, 2, 3], 0.5)];
        let lm = LumpedMassMatrix::from_elements(4, &elements, 2.0);
        for i in 0..4 {
            let m = lm.masses[i];
            let inv_m = lm.inv_mass(i);
            assert!((inv_m - 1.0 / m).abs() < 1e-12, "inv_mass mismatch at {i}");
        }
    }

    #[test]
    fn test_lumped_mass_apply_inv() {
        let elements = vec![([0usize, 1, 2, 3], 1.0 / 6.0)];
        let lm = LumpedMassMatrix::from_elements(4, &elements, 1.0);
        let forces = vec![[1.0, 0.0, 0.0]; 4];
        let accels = lm.apply_inv(&forces);
        // All accelerations should be inv_mass * 1.0
        for i in 0..4 {
            let expected = lm.inv_mass(i);
            assert!(
                (accels[i][0] - expected).abs() < 1e-10,
                "accel[{i}] mismatch"
            );
        }
    }

    // ── RayleighDamping ──────────────────────────────────────────────────

    #[test]
    fn test_rayleigh_mass_damping_force() {
        let rd = RayleighDamping::new(2.0, 0.0);
        let force = rd.mass_damping_force(3.0, [1.0, 0.0, 0.0]);
        // f = -alpha * m * v = -2.0 * 3.0 * 1.0 = -6.0
        assert!(
            (force[0] - (-6.0)).abs() < 1e-12,
            "mass damping force = {}",
            force[0]
        );
        assert!(force[1].abs() < 1e-15);
        assert!(force[2].abs() < 1e-15);
    }

    #[test]
    fn test_rayleigh_stiffness_damping_force() {
        let rd = RayleighDamping::new(0.0, 0.5);
        let elastic_f = [10.0, 0.0, 0.0];
        let force = rd.stiffness_damping_force(elastic_f);
        // f = -beta_K * f_elastic = -0.5 * 10 = -5.0
        assert!((force[0] - (-5.0)).abs() < 1e-12);
    }

    #[test]
    fn test_rayleigh_zero_damping_no_change() {
        let rd = RayleighDamping::new(0.0, 0.0);
        let mut forces = vec![[1.0_f64, 2.0, 3.0]; 4];
        let velocities = vec![[0.1_f64; 3]; 4];
        let masses = vec![1.0_f64; 4];
        let orig = forces.clone();
        rd.apply(&mut forces, &velocities, &masses);
        for i in 0..4 {
            for d in 0..3 {
                assert!(
                    (forces[i][d] - orig[i][d]).abs() < 1e-15,
                    "zero damping should not change forces"
                );
            }
        }
    }

    #[test]
    fn test_rayleigh_damping_ratio() {
        let rd = RayleighDamping::new(2.0, 0.01);
        let omega = 10.0; // rad/s
        let xi = rd.damping_ratio(omega);
        let expected = 2.0 / (2.0 * omega) + 0.01 * omega / 2.0;
        assert!((xi - expected).abs() < 1e-12);
    }

    #[test]
    fn test_rayleigh_from_two_modes() {
        let xi1 = 0.05;
        let omega1 = 10.0;
        let xi2 = 0.05;
        let omega2 = 100.0;
        let rd = RayleighDamping::from_two_modes(xi1, omega1, xi2, omega2);
        // Check that the damping ratio at each mode is approximately correct
        let xi1_check = rd.damping_ratio(omega1);
        let xi2_check = rd.damping_ratio(omega2);
        assert!(
            (xi1_check - xi1).abs() < 1e-8,
            "xi at omega1: got {xi1_check}, expected {xi1}"
        );
        assert!(
            (xi2_check - xi2).abs() < 1e-8,
            "xi at omega2: got {xi2_check}, expected {xi2}"
        );
    }

    // ── VelocityVerletIntegrator ─────────────────────────────────────────

    #[test]
    fn test_velocity_verlet_no_force_constant_velocity() {
        // With no forces, a particle should move at constant velocity
        let mut integrator = VelocityVerletIntegrator::new(1, 0.01);
        let mut positions = vec![[0.0_f64, 0.0, 0.0]];
        let mut velocities = vec![[1.0_f64, 0.0, 0.0]];
        let masses = vec![1.0];
        let tets: Vec<CorotFemTet> = vec![];
        let pinned = vec![false];
        integrator.step_sim(
            &mut positions,
            &mut velocities,
            &masses,
            &tets,
            500.0,
            300.0,
            [0.0; 3],
            None,
            &pinned,
        );
        // x should have moved by roughly dt * v = 0.01
        assert!(
            (positions[0][0] - 0.01).abs() < 1e-10,
            "position should be 0.01, got {}",
            positions[0][0]
        );
    }

    #[test]
    fn test_velocity_verlet_gravity_free_fall() {
        // Free fall: single node under gravity
        let mut integrator = VelocityVerletIntegrator::new(1, 0.01);
        let mut positions = vec![[0.0_f64, 10.0, 0.0]];
        let mut velocities = vec![[0.0_f64; 3]];
        let masses = vec![1.0];
        let tets: Vec<CorotFemTet> = vec![];
        let pinned = vec![false];

        for _ in 0..100 {
            integrator.step_sim(
                &mut positions,
                &mut velocities,
                &masses,
                &tets,
                500.0,
                300.0,
                [0.0, -9.81, 0.0],
                None,
                &pinned,
            );
        }
        // Should have fallen (y < 10)
        assert!(
            positions[0][1] < 10.0,
            "node should have fallen, y = {}",
            positions[0][1]
        );
        // Velocity should be downward
        assert!(
            velocities[0][1] < 0.0,
            "velocity should be downward, vy = {}",
            velocities[0][1]
        );
    }

    #[test]
    fn test_velocity_verlet_pinned_node_no_movement() {
        let mut integrator = VelocityVerletIntegrator::new(2, 0.01);
        let mut positions = vec![[0.0_f64; 3], [1.0, 0.0, 0.0]];
        let mut velocities = vec![[0.0_f64; 3]; 2];
        let masses = vec![1.0; 2];
        let tets: Vec<CorotFemTet> = vec![];
        let pinned = vec![true, false];

        integrator.step_sim(
            &mut positions,
            &mut velocities,
            &masses,
            &tets,
            500.0,
            300.0,
            [0.0, -9.81, 0.0],
            None,
            &pinned,
        );
        // Pinned node 0 should not move
        for d in 0..3 {
            assert!(
                positions[0][d].abs() < 1e-15,
                "pinned node should not move in dim {d}"
            );
        }
    }

    // ── BoundaryConditions ───────────────────────────────────────────────

    #[test]
    fn test_bc_pin_and_enforce() {
        let mut bc = BoundaryConditions::new();
        bc.pin_at(2, [3.0, 4.0, 5.0]);

        let mut positions = vec![[0.0_f64; 3]; 4];
        let mut velocities = vec![[1.0_f64; 3]; 4];

        bc.enforce(&mut positions, &mut velocities);

        // Node 2 should be reset
        assert_eq!(positions[2], [3.0, 4.0, 5.0]);
        assert_eq!(velocities[2], [0.0; 3]);

        // Other nodes untouched
        for i in [0usize, 1, 3] {
            for d in 0..3 {
                assert!(
                    (velocities[i][d] - 1.0).abs() < 1e-15,
                    "node {i} vel should be 1"
                );
            }
        }
    }

    #[test]
    fn test_bc_unpin() {
        let mut bc = BoundaryConditions::new();
        bc.pin_at(1, [0.0; 3]);
        assert!(bc.is_pinned(1));
        bc.unpin(1);
        assert!(!bc.is_pinned(1));
        assert_eq!(bc.len(), 0);
    }

    #[test]
    fn test_bc_pinned_mask() {
        let mut bc = BoundaryConditions::new();
        bc.pin_at(0, [0.0; 3]);
        bc.pin_at(3, [0.0; 3]);
        let mask = bc.pinned_mask(5);
        assert!(mask[0]);
        assert!(!mask[1]);
        assert!(!mask[2]);
        assert!(mask[3]);
        assert!(!mask[4]);
    }

    #[test]
    fn test_bc_is_empty() {
        let mut bc = BoundaryConditions::new();
        assert!(bc.is_empty());
        bc.pin_at(0, [0.0; 3]);
        assert!(!bc.is_empty());
    }

    // ── CorotFemElement4 ─────────────────────────────────────────────────

    #[test]
    fn test_corot_fem4_rest_volume() {
        let pos = unit_tet_pos();
        let elem = CorotFemElement4::new([0, 1, 2, 3], &pos, 1000.0, 0.3);
        assert!(
            (elem.rest_volume - 1.0 / 6.0).abs() < 1e-12,
            "rest volume should be 1/6, got {}",
            elem.rest_volume
        );
    }

    #[test]
    fn test_corot_fem4_zero_force_at_rest() {
        let pos = unit_tet_pos();
        let elem = CorotFemElement4::new([0, 1, 2, 3], &pos, 1000.0, 0.3);
        let forces = elem.elastic_forces(&pos);
        for k in 0..4 {
            let mag = (forces[k][0] * forces[k][0]
                + forces[k][1] * forces[k][1]
                + forces[k][2] * forces[k][2])
                .sqrt();
            assert!(
                mag < 1e-8,
                "elastic force at rest should be ~0 for node {k}, got {mag}"
            );
        }
    }

    #[test]
    fn test_corot_fem4_nonzero_force_when_stretched() {
        let pos = unit_tet_pos();
        let elem = CorotFemElement4::new([0, 1, 2, 3], &pos, 1000.0, 0.3);
        let mut stretched = pos.clone();
        stretched[1] = [2.0, 0.0, 0.0];
        let forces = elem.elastic_forces(&stretched);
        let mag = (forces[1][0] * forces[1][0]
            + forces[1][1] * forces[1][1]
            + forces[1][2] * forces[1][2])
            .sqrt();
        assert!(
            mag > 1e-6,
            "stretched element should have nonzero force, got {mag}"
        );
    }

    #[test]
    fn test_corot_fem4_newton_third_law() {
        let pos = unit_tet_pos();
        let elem = CorotFemElement4::new([0, 1, 2, 3], &pos, 1000.0, 0.3);
        let mut stretched = pos.clone();
        stretched[1] = [1.5, 0.0, 0.0];
        let forces = elem.elastic_forces(&stretched);
        let mut sum = [0.0_f64; 3];
        for k in 0..4 {
            for d in 0..3 {
                sum[d] += forces[k][d];
            }
        }
        for d in 0..3 {
            assert!(
                sum[d].abs() < 1e-8,
                "Newton 3rd law: force sum[{d}] = {} should be 0",
                sum[d]
            );
        }
    }

    #[test]
    fn test_corot_fem4_stiffness_matrix_size() {
        let pos = unit_tet_pos();
        let elem = CorotFemElement4::new([0, 1, 2, 3], &pos, 1000.0, 0.3);
        let k = elem.stiffness_matrix(&pos);
        assert_eq!(
            k.len(),
            144,
            "stiffness matrix should be 12×12 = 144 entries"
        );
    }

    #[test]
    fn test_corot_fem4_strain_energy_zero_at_rest() {
        let pos = unit_tet_pos();
        let elem = CorotFemElement4::new([0, 1, 2, 3], &pos, 1000.0, 0.3);
        let e = elem.strain_energy(&pos);
        assert!(
            e.abs() < 1e-8,
            "strain energy at rest should be ~0, got {e}"
        );
    }

    #[test]
    fn test_corot_fem4_strain_energy_positive_stretched() {
        let pos = unit_tet_pos();
        let elem = CorotFemElement4::new([0, 1, 2, 3], &pos, 1000.0, 0.3);
        let mut stretched = pos.clone();
        stretched[1] = [2.0, 0.0, 0.0];
        let e = elem.strain_energy(&stretched);
        assert!(
            e > 0.0,
            "strain energy should be positive when stretched, got {e}"
        );
    }

    // ── FemSoftBodyVerlet ────────────────────────────────────────────────

    fn make_verlet_body(dt: f64) -> FemSoftBodyVerlet {
        let pos = unit_tet_pos();
        let masses = vec![1.0; 4];
        let elem = CorotFemElement4::new([0, 1, 2, 3], &pos, 1000.0, 0.3);
        FemSoftBodyVerlet::new(pos, masses, vec![elem], [0.0, -9.81, 0.0], dt, 0.1, 0.0)
    }

    #[test]
    fn test_verlet_body_initial_ke_zero() {
        let body = make_verlet_body(0.001);
        assert!(
            body.kinetic_energy().abs() < 1e-15,
            "initial KE should be 0"
        );
    }

    #[test]
    fn test_verlet_body_step_moves_nodes() {
        // Velocity Verlet: first step accumulates velocity, second step moves positions.
        let mut body = make_verlet_body(0.001);
        let pos_before = body.positions.clone();
        body.step(); // computes accel from gravity, updates vel
        body.step(); // now positions move under nonzero velocity
        let moved = body
            .positions
            .iter()
            .zip(pos_before.iter())
            .any(|(a, b)| (0..3).any(|d| (a[d] - b[d]).abs() > 1e-15));
        assert!(moved, "nodes should move after two steps");
    }

    #[test]
    fn test_verlet_body_pinned_node_stays_put() {
        let mut body = make_verlet_body(0.001);
        body.pin_node(0);
        let pos0_before = body.positions[0];
        for _ in 0..10 {
            body.step();
        }
        for d in 0..3 {
            assert!(
                (body.positions[0][d] - pos0_before[d]).abs() < 1e-14,
                "pinned node 0 should not move in dim {d}"
            );
        }
    }

    #[test]
    fn test_verlet_body_deformed_has_positive_strain_energy() {
        let pos = unit_tet_pos();
        let masses = vec![1.0; 4];
        let elem = CorotFemElement4::new([0, 1, 2, 3], &pos, 1000.0, 0.3);
        let mut body = FemSoftBodyVerlet::new(pos, masses, vec![elem], [0.0; 3], 0.001, 0.0, 0.0);
        body.positions[1] = [2.0, 0.0, 0.0];
        assert!(
            body.strain_energy() > 0.0,
            "deformed body should have positive strain energy"
        );
    }

    #[test]
    fn test_verlet_body_rayleigh_reduces_energy() {
        // With strong damping, energy should decrease over time
        let pos = unit_tet_pos();
        let masses = vec![1.0; 4];
        let elem = CorotFemElement4::new([0, 1, 2, 3], &pos, 1000.0, 0.3);
        let mut body = FemSoftBodyVerlet::new(pos, masses, vec![elem], [0.0; 3], 0.001, 5.0, 0.0);
        // Start with deformed configuration
        body.positions[1] = [1.5, 0.0, 0.0];
        let e0 = body.mechanical_energy();
        for _ in 0..200 {
            body.step();
        }
        let e1 = body.mechanical_energy();
        assert!(
            e1 < e0,
            "Rayleigh damping should reduce mechanical energy: e0={e0}, e1={e1}"
        );
    }

    // ── assemble_internal_forces ─────────────────────────────────────────

    #[test]
    fn test_assemble_internal_forces_at_rest() {
        let pos = unit_tet_pos();
        let elem = CorotFemElement4::new([0, 1, 2, 3], &pos, 1000.0, 0.3);
        let forces = assemble_internal_forces(&pos, &[elem]);
        for i in 0..4 {
            let mag = (forces[i][0] * forces[i][0]
                + forces[i][1] * forces[i][1]
                + forces[i][2] * forces[i][2])
                .sqrt();
            assert!(
                mag < 1e-8,
                "internal force at rest should be ~0 for node {i}, got {mag}"
            );
        }
    }

    #[test]
    fn test_assemble_internal_forces_sum_zero() {
        let pos = unit_tet_pos();
        let elem = CorotFemElement4::new([0, 1, 2, 3], &pos, 1000.0, 0.3);
        let mut stretched = pos.clone();
        stretched[1] = [1.8, 0.0, 0.0];
        let forces = assemble_internal_forces(&stretched, &[elem]);
        let mut sum = [0.0_f64; 3];
        for f in &forces {
            for d in 0..3 {
                sum[d] += f[d];
            }
        }
        for d in 0..3 {
            assert!(
                sum[d].abs() < 1e-8,
                "internal force sum[{d}] = {} should be 0",
                sum[d]
            );
        }
    }
}
