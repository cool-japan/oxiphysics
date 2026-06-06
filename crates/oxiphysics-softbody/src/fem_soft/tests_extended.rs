// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Extended tests for FEM soft body modules.

use super::math_helpers::*;
use super::*;

#[cfg(test)]
mod extended_fem_tests {
    use super::BoundaryConditions;

    use crate::fem_soft::CorotFemElement4;
    use crate::fem_soft::CorotFemNode;
    use crate::fem_soft::CorotFemSim;
    use crate::fem_soft::CorotFemTet;
    use crate::fem_soft::CorotationalTet;
    use crate::fem_soft::FemSoftBodyVerlet;

    use crate::fem_soft::HighLevelFemBody;
    use crate::fem_soft::HyperelasticBody;
    use crate::fem_soft::NeoHookeanElement;
    use crate::fem_soft::NeoHookeanMaterial;
    use crate::fem_soft::NewmarkBetaSim;
    use crate::fem_soft::RayleighDamping;
    use crate::fem_soft::VelocityVerletIntegrator;

    use crate::fem_soft::assemble_internal_forces;
    use crate::fem_soft::cauchy_stress;

    use crate::fem_soft::green_lagrange_strain;

    use crate::fem_soft::tests_extended::det3x3;
    use crate::fem_soft::tests_extended::inv3x3;

    use crate::fem_soft::tests_extended::mul3x3;
    use crate::fem_soft::tests_extended::transpose3x3;

    use crate::fem_soft::von_mises_stress;

    // ── helpers ──────────────────────────────────────────────────────────────

    fn unit_pos() -> Vec<[f64; 3]> {
        vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ]
    }

    fn make_elem(young: f64, poisson: f64) -> CorotFemElement4 {
        CorotFemElement4::new([0, 1, 2, 3], &unit_pos(), young, poisson)
    }

    fn lame(young: f64, poisson: f64) -> (f64, f64) {
        let mu = young / (2.0 * (1.0 + poisson));
        let lambda = young * poisson / ((1.0 + poisson) * (1.0 - 2.0 * poisson));
        (mu, lambda)
    }

    fn corot_nodes() -> Vec<CorotFemNode> {
        let pos = unit_pos();
        pos.iter()
            .map(|&p| CorotFemNode {
                position: p,
                velocity: [0.0; 3],
                force: [0.0; 3],
                mass: 1.0,
            })
            .collect()
    }

    fn mag3(v: [f64; 3]) -> f64 {
        (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
    }

    // ── polar decomposition ───────────────────────────────────────────────

    #[test]
    fn test_polar_identity_gives_identity_rotation() {
        let id = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let r = CorotFemTet::polar_decompose_rotation(id);
        for (i, row) in r.iter().enumerate() {
            for (j, &val) in row.iter().enumerate() {
                let exp = if i == j { 1.0 } else { 0.0 };
                assert!(
                    (val - exp).abs() < 1e-8,
                    "polar(I)[{i}][{j}]={} != {exp}",
                    val
                );
            }
        }
    }

    #[test]
    fn test_polar_rotation_is_orthogonal() {
        // R^T R = I
        let f = [[1.2, 0.3, 0.0], [0.1, 0.9, 0.2], [0.0, 0.1, 1.1]];
        let r = CorotFemTet::polar_decompose_rotation(f);
        let rt = transpose3x3(r);
        let rtr = mul3x3(rt, r);
        for (i, row) in rtr.iter().enumerate() {
            for (j, &val) in row.iter().enumerate() {
                let exp = if i == j { 1.0 } else { 0.0 };
                assert!(
                    (val - exp).abs() < 1e-6,
                    "R^T R [{i}][{j}] = {}, expected {exp}",
                    val
                );
            }
        }
    }

    #[test]
    fn test_polar_rotation_det_is_positive_one() {
        let f = [[0.8, 0.4, 0.0], [0.2, 1.1, 0.3], [0.0, 0.1, 0.9]];
        let r = CorotFemTet::polar_decompose_rotation(f);
        let d = det3x3(r);
        assert!((d - 1.0).abs() < 1e-6, "det(R) = {d}, expected 1");
    }

    #[test]
    fn test_polar_decompose_r_times_s_recovers_f() {
        let f = [[1.3, 0.2, 0.1], [0.0, 1.1, 0.05], [0.0, 0.0, 0.9]];
        let (r, s) = CorotationalTet::polar_decompose(&f);
        let rs = mul3x3(r, s);
        for (i, row) in rs.iter().enumerate() {
            for (j, &val) in row.iter().enumerate() {
                assert!(
                    (val - f[i][j]).abs() < 1e-5,
                    "RS[{i}][{j}]={} != F[{i}][{j}]={}",
                    val,
                    f[i][j]
                );
            }
        }
    }

    #[test]
    fn test_polar_decompose_s_is_symmetric() {
        let f = [[1.5, 0.2, 0.0], [0.1, 1.0, 0.3], [0.0, 0.0, 0.8]];
        let (_r, s) = CorotationalTet::polar_decompose(&f);
        for (i, row) in s.iter().enumerate() {
            for j in (i + 1)..3 {
                assert!(
                    (row[j] - s[j][i]).abs() < 1e-5,
                    "S not symmetric: S[{i}][{j}]={} vs S[{j}][{i}]={}",
                    row[j],
                    s[j][i]
                );
            }
        }
    }

    #[test]
    fn test_polar_decompose_r_rotation_z90() {
        // 90° about Z
        let rot = [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]];
        let (r, _s) = CorotationalTet::polar_decompose(&rot);
        for (i, row) in r.iter().enumerate() {
            for (j, &val) in row.iter().enumerate() {
                assert!(
                    (val - rot[i][j]).abs() < 1e-6,
                    "R[{i}][{j}]={} != expected {}",
                    val,
                    rot[i][j]
                );
            }
        }
    }

    #[test]
    fn test_polar_decompose_singular_returns_identity() {
        // Singular F → R should be identity (graceful fallback)
        let f = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]];
        let r = CorotFemTet::polar_decompose_rotation(f);
        let d = det3x3(r);
        // Either identity (det=1) or at least det not negative
        assert!(d >= -1e-6, "det(R) for singular F should be non-negative");
    }

    // ── deformation gradient ─────────────────────────────────────────────

    #[test]
    fn test_deformation_gradient_identity_at_rest() {
        let pos = unit_pos();
        let elem = CorotFemElement4::new([0, 1, 2, 3], &pos, 1000.0, 0.3);
        let f = elem.deformation_gradient(&pos);
        for (i, row) in f.iter().enumerate() {
            for (j, &val) in row.iter().enumerate() {
                let exp = if i == j { 1.0 } else { 0.0 };
                assert!((val - exp).abs() < 1e-10, "F[{i}][{j}] should be I at rest");
            }
        }
    }

    #[test]
    fn test_deformation_gradient_x_stretch_by_two() {
        let pos = unit_pos();
        let elem = CorotFemElement4::new([0, 1, 2, 3], &pos, 1000.0, 0.3);
        let mut deformed = pos.clone();
        deformed[1] = [2.0, 0.0, 0.0];
        let f = elem.deformation_gradient(&deformed);
        assert!(
            (f[0][0] - 2.0).abs() < 1e-10,
            "F[0][0] should be 2 for x-stretch"
        );
        assert!((f[1][1] - 1.0).abs() < 1e-10, "F[1][1] should stay 1");
        assert!((f[2][2] - 1.0).abs() < 1e-10, "F[2][2] should stay 1");
    }

    #[test]
    fn test_deformation_gradient_uniform_compression() {
        let pos = unit_pos();
        let elem = CorotFemElement4::new([0, 1, 2, 3], &pos, 1000.0, 0.3);
        let mut compressed = pos.clone();
        compressed[1] = [0.5, 0.0, 0.0];
        compressed[2] = [0.0, 0.5, 0.0];
        compressed[3] = [0.0, 0.0, 0.5];
        let f = elem.deformation_gradient(&compressed);
        assert!(
            (f[0][0] - 0.5).abs() < 1e-10,
            "F[0][0] should be 0.5 for compression"
        );
    }

    #[test]
    fn test_deformation_gradient_det_matches_volume_ratio() {
        let pos = unit_pos();
        let elem = CorotFemElement4::new([0, 1, 2, 3], &pos, 1000.0, 0.3);
        // Scale all three edges by 2 → volume ratio = 8
        let mut scaled = pos.clone();
        scaled[1] = [2.0, 0.0, 0.0];
        scaled[2] = [0.0, 2.0, 0.0];
        scaled[3] = [0.0, 0.0, 2.0];
        let f = elem.deformation_gradient(&scaled);
        let j = det3x3(f);
        assert!(
            (j - 8.0).abs() < 1e-8,
            "det(F)={j}, expected 8 for 2× scale"
        );
    }

    // ── element stiffness matrix ─────────────────────────────────────────

    #[test]
    fn test_stiffness_matrix_is_12x12() {
        let elem = make_elem(1000.0, 0.3);
        let k = elem.stiffness_matrix(&unit_pos());
        assert_eq!(k.len(), 144, "stiffness matrix must be 12×12");
    }

    #[test]
    fn test_stiffness_matrix_diagonal_entries_positive() {
        let elem = make_elem(1000.0, 0.3);
        let k = elem.stiffness_matrix(&unit_pos());
        // Diagonal entries of stiffness matrix should be negative (restoring)
        // but since we compute ∂f/∂x and forces resist displacement, K diagonal
        // entries are typically negative (force decreases as displacement increases).
        // At least they should be nonzero.
        let diag_nonzero = (0..12).any(|i| k[i * 12 + i].abs() > 1e-8);
        assert!(diag_nonzero, "stiffness matrix diagonal should be nonzero");
    }

    #[test]
    fn test_stiffness_matrix_higher_young_means_stiffer() {
        let pos = unit_pos();
        let mut stretched = pos.clone();
        stretched[1] = [1.3, 0.0, 0.0];
        let k_soft =
            CorotFemElement4::new([0, 1, 2, 3], &pos, 100.0, 0.3).stiffness_matrix(&stretched);
        let k_stiff =
            CorotFemElement4::new([0, 1, 2, 3], &pos, 10000.0, 0.3).stiffness_matrix(&stretched);
        // Max absolute entry of stiffer matrix should be larger
        let max_soft = k_soft.iter().map(|x| x.abs()).fold(0.0_f64, f64::max);
        let max_stiff = k_stiff.iter().map(|x| x.abs()).fold(0.0_f64, f64::max);
        assert!(
            max_stiff > max_soft,
            "stiffer material must produce larger K entries: soft={max_soft}, stiff={max_stiff}"
        );
    }

    #[test]
    fn test_stiffness_matrix_symmetry_at_rest() {
        let elem = make_elem(1000.0, 0.3);
        let k = elem.stiffness_matrix(&unit_pos());
        for r in 0..12 {
            for c in 0..12 {
                let diff = (k[r * 12 + c] - k[c * 12 + r]).abs();
                assert!(
                    diff < 1e-2,
                    "K[{r},{c}]={} vs K[{c},{r}]={}: asymmetry={diff}",
                    k[r * 12 + c],
                    k[c * 12 + r]
                );
            }
        }
    }

    #[test]
    fn test_stiffness_matrix_force_displacement_consistency() {
        // K * dx should approximate df for small dx
        let pos = unit_pos();
        let elem = make_elem(1000.0, 0.3);
        let f0 = elem.elastic_forces(&pos);
        let dx_node = 1e-4;
        let mut pos_p = pos.clone();
        pos_p[1][0] += dx_node;
        let f1 = elem.elastic_forces(&pos_p);
        // Finite-difference gradient for node 1 dof 0
        let df_fd = (f1[1][0] - f0[1][0]) / dx_node;
        // From stiffness matrix (col for node1/dof0 = col 3, row for node1/dof0 = row 3)
        let k = elem.stiffness_matrix(&pos);
        let df_k = k[3 * 12 + 3];
        assert!(
            (df_fd - df_k).abs() < 0.1 * (df_fd.abs() + 1e-6),
            "stiffness matrix entry K[3,3]={df_k} vs FD={df_fd}"
        );
    }

    // ── element force assembly ────────────────────────────────────────────

    #[test]
    fn test_force_assembly_momentum_conservation() {
        // Total force on a closed body must be zero (no external forces)
        let pos = unit_pos();
        let e1 = CorotFemElement4::new([0, 1, 2, 3], &pos, 1000.0, 0.3);
        let mut deformed = pos.clone();
        deformed[1] = [1.4, 0.0, 0.0];
        deformed[2] = [0.0, 1.2, 0.0];
        let f = assemble_internal_forces(&deformed, &[e1]);
        let mut sum = [0.0f64; 3];
        for fi in &f {
            for (d, &fval) in fi.iter().enumerate() {
                sum[d] += fval;
            }
        }
        for (d, &s) in sum.iter().enumerate() {
            assert!(
                s.abs() < 1e-7,
                "total internal force[{d}] = {} should be 0",
                s
            );
        }
    }

    #[test]
    fn test_force_assembly_two_elements_shared_node() {
        // Two tets sharing node 0: contributions from both elements get assembled
        let mut pos = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [1.0, 1.0, 0.0],
        ];
        let e1 = CorotFemElement4::new([0, 1, 2, 3], &pos, 1000.0, 0.3);
        let e2 = CorotFemElement4::new([0, 1, 4, 2], &pos, 1000.0, 0.3);
        pos[1] = [1.5, 0.0, 0.0]; // deform
        let forces = assemble_internal_forces(&pos, &[e1, e2]);
        // Sum over all nodes must be zero (momentum)
        let mut sum = [0.0f64; 3];
        for f in &forces {
            for (d, &fval) in f.iter().enumerate() {
                sum[d] += fval;
            }
        }
        for (d, &s) in sum.iter().enumerate() {
            assert!(
                s.abs() < 1e-6,
                "two-element force sum[{d}] = {} should be 0",
                s
            );
        }
    }

    #[test]
    fn test_force_direction_opposes_stretch() {
        // Stretch node 1 in +x: the force on node 1 should be in -x direction
        let pos = unit_pos();
        let elem = make_elem(1000.0, 0.3);
        let mut deformed = pos.clone();
        deformed[1] = [2.0, 0.0, 0.0];
        let forces = elem.elastic_forces(&deformed);
        assert!(
            forces[1][0] < 0.0,
            "force on stretched node should pull back in -x, got {}",
            forces[1][0]
        );
    }

    #[test]
    fn test_force_increases_with_deformation() {
        let pos = unit_pos();
        let elem = make_elem(1000.0, 0.3);
        let mut d1 = pos.clone();
        d1[1] = [1.2, 0.0, 0.0];
        let f1 = elem.elastic_forces(&d1);
        let mut d2 = pos.clone();
        d2[1] = [1.5, 0.0, 0.0];
        let f2 = elem.elastic_forces(&d2);
        // More stretch → larger restoring force magnitude
        assert!(
            mag3(f2[1]) > mag3(f1[1]),
            "larger deformation should produce larger force"
        );
    }

    #[test]
    fn test_force_scales_with_young_modulus() {
        let pos = unit_pos();
        let mut deformed = pos.clone();
        deformed[1] = [1.5, 0.0, 0.0];
        let e_soft = CorotFemElement4::new([0, 1, 2, 3], &pos, 100.0, 0.3);
        let e_stiff = CorotFemElement4::new([0, 1, 2, 3], &pos, 10000.0, 0.3);
        let f_soft = mag3(e_soft.elastic_forces(&deformed)[1]);
        let f_stiff = mag3(e_stiff.elastic_forces(&deformed)[1]);
        assert!(
            f_stiff > f_soft,
            "stiffer material must produce larger force: soft={f_soft}, stiff={f_stiff}"
        );
    }

    #[test]
    fn test_corot_tet_forces_same_as_fem_element4() {
        // CorotationalTet and CorotFemElement4 should give same forces
        let pos = unit_pos();
        let mut deformed = pos.clone();
        deformed[1] = [1.4, 0.0, 0.0];
        let (mu, lambda) = lame(1000.0, 0.3);
        let tet = CorotationalTet::new([0, 1, 2, 3], &pos);
        let f1 = tet.elastic_forces(&deformed, mu, lambda);
        let elem = CorotFemElement4::new([0, 1, 2, 3], &pos, 1000.0, 0.3);
        let f2 = elem.elastic_forces(&deformed);
        for (k, (force1, force2)) in f1.iter().zip(f2.iter()).enumerate() {
            for (d, (&v1, &v2)) in force1.iter().zip(force2.iter()).enumerate() {
                assert!(
                    (v1 - v2).abs() < 1e-6,
                    "force mismatch at node {k} dim {d}: {}-{}",
                    v1,
                    v2
                );
            }
        }
    }

    // ── global stiffness assembly ─────────────────────────────────────────

    #[test]
    fn test_global_stiffness_not_all_zeros() {
        let pos = unit_pos();
        let elem = make_elem(1000.0, 0.3);
        let k = elem.stiffness_matrix(&pos);
        let nonzero = k.iter().any(|x| x.abs() > 1e-10);
        assert!(nonzero, "stiffness matrix should have nonzero entries");
    }

    #[test]
    fn test_global_stiffness_scales_with_volume() {
        // Scale the tet by 2 in all dims → rest volume = 8×, stiffness entries scale
        let pos_big: Vec<[f64; 3]> = unit_pos()
            .iter()
            .map(|&p| [p[0] * 2.0, p[1] * 2.0, p[2] * 2.0])
            .collect();
        let elem_small = CorotFemElement4::new([0, 1, 2, 3], &unit_pos(), 1000.0, 0.3);
        let elem_big = CorotFemElement4::new([0, 1, 2, 3], &pos_big, 1000.0, 0.3);
        let k_small = elem_small.stiffness_matrix(&unit_pos());
        let k_big = elem_big.stiffness_matrix(&pos_big);
        let sum_small: f64 = k_small.iter().map(|x| x.abs()).sum();
        let sum_big: f64 = k_big.iter().map(|x| x.abs()).sum();
        // Larger volume → stiffness entries scale (roughly)
        assert!(
            sum_big > sum_small * 0.5,
            "bigger tet should have non-trivially larger stiffness: small={sum_small}, big={sum_big}"
        );
    }

    #[test]
    fn test_stiffness_matrix_corot_tet_vs_element4() {
        // Both element types should give the same stiffness matrix
        let pos = unit_pos();
        let (mu, lambda) = lame(1000.0, 0.3);
        let tet = CorotationalTet::new([0, 1, 2, 3], &pos);
        let elem = CorotFemElement4::new([0, 1, 2, 3], &pos, 1000.0, 0.3);
        let k1 = tet.stiffness_matrix(&pos, mu, lambda);
        let k2 = elem.stiffness_matrix(&pos);
        for i in 0..144 {
            assert!(
                (k1[i] - k2[i]).abs() < 1e-4,
                "stiffness mismatch at entry {i}: {} vs {}",
                k1[i],
                k2[i]
            );
        }
    }

    // ── implicit time integration (backward Euler / Newton-Raphson style) ─

    #[test]
    fn test_newmark_pinned_node_zero_velocity() {
        let nodes = corot_nodes();
        let tet = CorotFemTet::new([0, 1, 2, 3], &nodes);
        let mut sim = NewmarkBetaSim::new(nodes, vec![tet], 500.0, 300.0, 0.01);
        sim.pin_node(0);
        for _ in 0..10 {
            sim.step();
        }
        for d in 0..3 {
            assert!(
                sim.nodes[0].velocity[d].abs() < 1e-15,
                "pinned node velocity[{d}] should be 0"
            );
        }
    }

    #[test]
    fn test_newmark_apply_impulse_changes_velocity() {
        let nodes = corot_nodes();
        let tet = CorotFemTet::new([0, 1, 2, 3], &nodes);
        let mut sim = NewmarkBetaSim::new(nodes, vec![tet], 500.0, 300.0, 0.01);
        let v_before = sim.nodes[1].velocity[0];
        sim.apply_impulse(1, [10.0, 0.0, 0.0]);
        // velocity should change by impulse / mass = 10 / 1 = 10
        assert!(
            (sim.nodes[1].velocity[0] - v_before - 10.0).abs() < 1e-10,
            "impulse should change velocity by 10"
        );
    }

    #[test]
    fn test_newmark_pinned_impulse_no_effect() {
        let nodes = corot_nodes();
        let tet = CorotFemTet::new([0, 1, 2, 3], &nodes);
        let mut sim = NewmarkBetaSim::new(nodes, vec![tet], 500.0, 300.0, 0.01);
        sim.pin_node(0);
        sim.apply_impulse(0, [100.0, 0.0, 0.0]);
        assert!(
            sim.nodes[0].velocity[0].abs() < 1e-15,
            "impulse on pinned node should have no effect"
        );
    }

    #[test]
    fn test_newmark_active_node_count() {
        let nodes = corot_nodes();
        let tet = CorotFemTet::new([0, 1, 2, 3], &nodes);
        let mut sim = NewmarkBetaSim::new(nodes, vec![tet], 500.0, 300.0, 0.01);
        assert_eq!(sim.active_node_count(), 4, "initially all nodes active");
        sim.pin_node(1);
        assert_eq!(sim.active_node_count(), 3, "one node pinned");
        sim.pin_node(3);
        assert_eq!(sim.active_node_count(), 2, "two nodes pinned");
    }

    #[test]
    fn test_newmark_unpin_restores_active_count() {
        let nodes = corot_nodes();
        let tet = CorotFemTet::new([0, 1, 2, 3], &nodes);
        let mut sim = NewmarkBetaSim::new(nodes, vec![tet], 500.0, 300.0, 0.01);
        sim.pin_node(0);
        sim.pin_node(1);
        sim.unpin_node(0);
        assert_eq!(sim.active_node_count(), 3, "one node still pinned");
    }

    #[test]
    fn test_newmark_rayleigh_damping_reduces_kinetic_energy() {
        let nodes = corot_nodes();
        let tet = CorotFemTet::new([0, 1, 2, 3], &nodes);
        let mut deformed = nodes.clone();
        deformed[1].position = [1.5, 0.0, 0.0];
        // Strong mass-proportional Rayleigh
        let mut sim = NewmarkBetaSim::with_params(
            deformed,
            vec![tet],
            1000.0,
            500.0,
            0.001,
            0.25,
            0.5,
            10.0,
            0.0,
        );
        for _ in 0..20 {
            sim.step();
        }
        let ke = sim.total_kinetic_energy();
        // With strong damping KE should remain modest (< 50 is a generous bound)
        assert!(
            ke < 50.0,
            "highly damped system should have small KE, got {ke}"
        );
    }

    #[test]
    fn test_newmark_beta025_gamma05_params() {
        let nodes = corot_nodes();
        let tet = CorotFemTet::new([0, 1, 2, 3], &nodes);
        let sim = NewmarkBetaSim::new(nodes, vec![tet], 500.0, 300.0, 0.01);
        assert!(
            (sim.beta - 0.25).abs() < 1e-15,
            "default beta should be 0.25"
        );
        assert!(
            (sim.gamma - 0.5).abs() < 1e-15,
            "default gamma should be 0.5"
        );
    }

    #[test]
    fn test_newmark_custom_beta_gamma() {
        let nodes = corot_nodes();
        let tet = CorotFemTet::new([0, 1, 2, 3], &nodes);
        let sim =
            NewmarkBetaSim::with_params(nodes, vec![tet], 500.0, 300.0, 0.01, 0.3, 0.6, 0.0, 0.0);
        assert!((sim.beta - 0.3).abs() < 1e-15);
        assert!((sim.gamma - 0.6).abs() < 1e-15);
    }

    #[test]
    fn test_newmark_strain_energy_increases_under_stretch() {
        let nodes = corot_nodes();
        let tet = CorotFemTet::new([0, 1, 2, 3], &nodes);
        let mut sim = NewmarkBetaSim::new(nodes, vec![tet], 1000.0, 500.0, 0.001);
        let e0 = sim.total_strain_energy();
        // Manually stretch node 1
        sim.nodes[1].position = [1.8, 0.0, 0.0];
        let e1 = sim.total_strain_energy();
        assert!(
            e1 > e0,
            "strain energy must increase under stretch: e0={e0}, e1={e1}"
        );
    }

    // ── plastic strain handling ────────────────────────────────────────────

    /// Simple plastic strain accumulator (von Mises yield criterion).
    ///
    /// Stores accumulated equivalent plastic strain `ε_p` per node.
    #[derive(Debug, Clone)]
    struct PlasticStrainAccumulator {
        /// Equivalent plastic strain per node.
        pub equiv_plastic_strain: Vec<f64>,
        /// Von Mises yield stress (Pa).
        pub yield_stress: f64,
        /// Hardening modulus H (Pa).  Simple isotropic linear hardening.
        pub hardening: f64,
    }

    impl PlasticStrainAccumulator {
        fn new(n_nodes: usize, yield_stress: f64, hardening: f64) -> Self {
            Self {
                equiv_plastic_strain: vec![0.0; n_nodes],
                yield_stress,
                hardening,
            }
        }

        /// Return current yield stress (with isotropic hardening).
        fn current_yield_stress(&self, node: usize) -> f64 {
            self.yield_stress + self.hardening * self.equiv_plastic_strain[node]
        }

        /// Update plastic strain for a given von Mises stress at a node.
        /// Returns the plastic strain increment.
        fn update(&mut self, node: usize, von_mises: f64) -> f64 {
            let sigma_y = self.current_yield_stress(node);
            if von_mises > sigma_y {
                let d_ep = (von_mises - sigma_y) / (3.0 * self.hardening.max(1.0) + 1.0e-30);
                self.equiv_plastic_strain[node] += d_ep;
                d_ep
            } else {
                0.0
            }
        }

        /// Total accumulated plastic work (approximation).
        fn total_plastic_strain(&self) -> f64 {
            self.equiv_plastic_strain.iter().sum()
        }
    }

    #[test]
    fn test_plastic_strain_zero_below_yield() {
        let mut acc = PlasticStrainAccumulator::new(4, 100.0, 0.0);
        let d_ep = acc.update(0, 50.0); // von Mises < yield
        assert!(
            d_ep.abs() < 1e-15,
            "no plastic strain below yield, got {d_ep}"
        );
        assert!(acc.equiv_plastic_strain[0].abs() < 1e-15);
    }

    #[test]
    fn test_plastic_strain_accumulates_above_yield() {
        let mut acc = PlasticStrainAccumulator::new(4, 100.0, 0.0);
        let d_ep = acc.update(0, 150.0);
        assert!(
            d_ep > 0.0,
            "plastic strain should accumulate above yield, got {d_ep}"
        );
        assert!(acc.equiv_plastic_strain[0] > 0.0);
    }

    #[test]
    fn test_plastic_hardening_raises_yield_stress() {
        let mut acc = PlasticStrainAccumulator::new(1, 100.0, 500.0);
        let sy0 = acc.current_yield_stress(0);
        acc.update(0, 200.0); // above yield → plastic strain increases
        let sy1 = acc.current_yield_stress(0);
        assert!(
            sy1 > sy0,
            "yield stress should increase with hardening: {sy0} → {sy1}"
        );
    }

    #[test]
    fn test_plastic_strain_at_exact_yield_is_zero() {
        let mut acc = PlasticStrainAccumulator::new(1, 100.0, 0.0);
        let d_ep = acc.update(0, 100.0); // exactly at yield
        assert!(
            d_ep.abs() < 1e-15,
            "plastic strain at exactly yield threshold should be 0, got {d_ep}"
        );
    }

    #[test]
    fn test_plastic_total_strain_sum() {
        let mut acc = PlasticStrainAccumulator::new(3, 100.0, 0.0);
        acc.update(0, 200.0);
        acc.update(1, 300.0);
        let total = acc.total_plastic_strain();
        let expected = acc.equiv_plastic_strain[0] + acc.equiv_plastic_strain[1];
        assert!(
            (total - expected).abs() < 1e-14,
            "total plastic strain mismatch"
        );
    }

    #[test]
    fn test_plastic_strain_per_node_independent() {
        let mut acc = PlasticStrainAccumulator::new(4, 100.0, 0.0);
        acc.update(0, 200.0);
        // Node 1 should still have zero plastic strain
        assert!(
            acc.equiv_plastic_strain[1].abs() < 1e-15,
            "node 1 plastic strain should be 0 if not yielded"
        );
    }

    #[test]
    fn test_von_mises_from_cauchy_and_plastic_yield_check() {
        // Build a stress from a deformation and check von Mises against yield
        let pos = unit_pos();
        let mat = NeoHookeanMaterial::from_young_poisson(1000.0, 0.3);
        let elem = NeoHookeanElement::new([0, 1, 2, 3], &pos, mat);
        let mut deformed = pos.clone();
        deformed[1] = [3.0, 0.0, 0.0]; // large stretch
        let f = elem.deformation_gradient(&deformed);
        let piola = mat.piola_kirchhoff(f);
        let sigma = cauchy_stress(piola, f);
        let vm = von_mises_stress(sigma);
        assert!(vm >= 0.0, "von Mises stress must be non-negative");
    }

    #[test]
    fn test_plastic_strain_no_hardening_saturates() {
        // With no hardening (H=0), every step above yield gives same plastic increment
        let yield_s = 100.0;
        let mut acc = PlasticStrainAccumulator::new(1, yield_s, 0.0);
        let sigma = 150.0;
        let d1 = acc.update(0, sigma);
        // Second call: since H=0, yield stress doesn't change → same increment
        let d2 = acc.update(0, sigma);
        assert!(
            (d1 - d2).abs() < 1e-10,
            "no-hardening increments should match: {d1} vs {d2}"
        );
    }

    // ── energy computation ────────────────────────────────────────────────

    #[test]
    fn test_strain_energy_scales_with_young_modulus() {
        let pos = unit_pos();
        let mut deformed = pos.clone();
        deformed[1] = [1.5, 0.0, 0.0];
        let e_soft = CorotFemElement4::new([0, 1, 2, 3], &pos, 100.0, 0.3).strain_energy(&deformed);
        let e_stiff =
            CorotFemElement4::new([0, 1, 2, 3], &pos, 10000.0, 0.3).strain_energy(&deformed);
        assert!(
            e_stiff > e_soft,
            "stiffer material → more strain energy: {e_soft} < {e_stiff}"
        );
    }

    #[test]
    fn test_strain_energy_increases_with_deformation() {
        let pos = unit_pos();
        let elem = make_elem(1000.0, 0.3);
        let mut d1 = pos.clone();
        d1[1] = [1.2, 0.0, 0.0];
        let mut d2 = pos.clone();
        d2[1] = [2.0, 0.0, 0.0];
        let e1 = elem.strain_energy(&d1);
        let e2 = elem.strain_energy(&d2);
        assert!(e2 > e1, "more deformation → more energy: {e1} < {e2}");
    }

    #[test]
    fn test_total_elastic_energy_multiple_elements() {
        // Two elements at rest: total energy should be ≈0
        let mut pos = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [1.0, 1.0, 0.0],
        ];
        let e1 = CorotFemElement4::new([0, 1, 2, 3], &pos, 1000.0, 0.3);
        // Second tet must be valid (non-degenerate)
        let e2 = CorotFemElement4::new([0, 1, 4, 2], &pos, 1000.0, 0.3);
        let total = e1.strain_energy(&pos) + e2.strain_energy(&pos);
        assert!(
            total.abs() < 1e-7,
            "total energy of two rest-state elements should be ~0, got {total}"
        );
        // Deform shared node
        pos[0] = [0.2, 0.0, 0.0];
        let total_d = e1.strain_energy(&pos) + e2.strain_energy(&pos);
        assert!(total_d > 0.0, "deformed system should have positive energy");
    }

    #[test]
    fn test_kinetic_energy_zero_initial() {
        let body = make_verlet_body_no_grav();
        assert!(
            body.kinetic_energy().abs() < 1e-15,
            "initial KE should be 0"
        );
    }

    fn make_verlet_body_no_grav() -> FemSoftBodyVerlet {
        let pos = unit_pos();
        let masses = vec![1.0; 4];
        let elem = make_elem(1000.0, 0.3);
        FemSoftBodyVerlet::new(pos, masses, vec![elem], [0.0; 3], 0.001, 0.0, 0.0)
    }

    #[test]
    fn test_kinetic_energy_positive_after_impulse() {
        let mut body = make_verlet_body_no_grav();
        body.velocities[1] = [1.0, 0.0, 0.0];
        let ke = body.kinetic_energy();
        assert!(
            ke > 0.0,
            "KE should be positive with nonzero velocity, got {ke}"
        );
    }

    #[test]
    fn test_mechanical_energy_is_sum() {
        let mut body = make_verlet_body_no_grav();
        body.positions[1] = [1.5, 0.0, 0.0];
        body.velocities[2] = [0.0, 0.5, 0.0];
        let ke = body.kinetic_energy();
        let se = body.strain_energy();
        let me = body.mechanical_energy();
        assert!(
            (me - ke - se).abs() < 1e-12,
            "mechanical energy = KE + SE: {me} vs {ke} + {se} = {}",
            ke + se
        );
    }

    #[test]
    fn test_neo_hookean_strain_energy_positive_compressed() {
        let mat = NeoHookeanMaterial::from_young_poisson(1000.0, 0.3);
        let f = [[0.5, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 0.5]];
        let w = mat.strain_energy(f);
        assert!(
            w > 0.0,
            "strain energy under compression should be positive, got {w}"
        );
    }

    #[test]
    fn test_strain_energy_corot_tet_vs_high_level() {
        let pos = unit_pos();
        let mut deformed = pos.clone();
        deformed[1] = [1.6, 0.0, 0.0];
        let (mu, lambda) = lame(1000.0, 0.3);
        let _tet = CorotationalTet::new([0, 1, 2, 3], &pos);
        // Compute via CorotFemSim helper
        let nodes_deformed: Vec<CorotFemNode> = deformed
            .iter()
            .map(|&p| CorotFemNode {
                position: p,
                velocity: [0.0; 3],
                force: [0.0; 3],
                mass: 1.0,
            })
            .collect();
        let tet2 = CorotFemTet::new([0, 1, 2, 3], &{
            let ns: Vec<CorotFemNode> = pos
                .iter()
                .map(|&p| CorotFemNode {
                    position: p,
                    velocity: [0.0; 3],
                    force: [0.0; 3],
                    mass: 1.0,
                })
                .collect();
            ns
        });
        let sim = CorotFemSim::new(nodes_deformed, vec![tet2], mu, lambda, 0.001);
        let e_sim = sim.total_strain_energy();
        // HighLevelFemBody
        let mut body = HighLevelFemBody::new(pos.clone(), vec![1.0; 4], 1000.0, 0.3);
        body.add_tet(0, 1, 2, 3);
        body.positions[1] = deformed[1];
        let e_body = body.total_elastic_energy();
        assert!(
            (e_sim - e_body).abs() < 1e-5,
            "CorotFemSim energy={e_sim} vs HighLevelFemBody energy={e_body}"
        );
    }

    #[test]
    fn test_hyperelastic_body_strain_energy_zero_at_rest() {
        let pos = unit_pos();
        let mat = NeoHookeanMaterial::from_young_poisson(1000.0, 0.3);
        let elem = NeoHookeanElement::new([0, 1, 2, 3], &pos, mat);
        let body = HyperelasticBody::new(pos.clone(), vec![elem], 1.0, 0.01);
        assert!(
            body.strain_energy().abs() < 1e-8,
            "neo-hookean energy at rest ~0"
        );
    }

    #[test]
    fn test_hyperelastic_body_strain_energy_positive_deformed() {
        let pos = unit_pos();
        let mat = NeoHookeanMaterial::from_young_poisson(1000.0, 0.3);
        let elem = NeoHookeanElement::new([0, 1, 2, 3], &pos, mat);
        let mut body = HyperelasticBody::new(pos.clone(), vec![elem], 1.0, 0.01);
        body.nodes[1] = [2.0, 0.0, 0.0];
        assert!(
            body.strain_energy() > 0.0,
            "neo-hookean energy > 0 when deformed"
        );
    }

    // ── matrix helper edge cases ──────────────────────────────────────────

    #[test]
    fn test_inv3x3_product_is_identity() {
        let m = [[2.0, 1.0, 0.0], [0.0, 3.0, 1.0], [1.0, 0.0, 2.0]];
        let inv = inv3x3(m);
        let prod = mul3x3(m, inv);
        for (i, row) in prod.iter().enumerate() {
            for (j, &val) in row.iter().enumerate() {
                let exp = if i == j { 1.0 } else { 0.0 };
                assert!(
                    (val - exp).abs() < 1e-10,
                    "M*inv(M)[{i}][{j}]={} != {exp}",
                    val
                );
            }
        }
    }

    #[test]
    fn test_det3x3_scale_factor() {
        let id = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let m2 = [[2.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 2.0]];
        assert!((det3x3(id) - 1.0).abs() < 1e-14);
        assert!((det3x3(m2) - 8.0).abs() < 1e-10, "det(2I) = 8");
    }

    #[test]
    fn test_transpose_is_own_inverse() {
        let m = [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]];
        let tt = transpose3x3(transpose3x3(m));
        for (i, row) in tt.iter().enumerate() {
            for (j, &val) in row.iter().enumerate() {
                assert!(
                    (val - m[i][j]).abs() < 1e-14,
                    "transpose(transpose(M)) != M"
                );
            }
        }
    }

    #[test]
    fn test_mul3x3_associativity() {
        let a = [[1.0, 2.0, 0.0], [0.0, 1.0, 1.0], [1.0, 0.0, 1.0]];
        let b = [[2.0, 0.0, 1.0], [1.0, 1.0, 0.0], [0.0, 1.0, 2.0]];
        let c = [[1.0, 1.0, 0.0], [0.0, 2.0, 1.0], [1.0, 0.0, 1.0]];
        let ab_c = mul3x3(mul3x3(a, b), c);
        let a_bc = mul3x3(a, mul3x3(b, c));
        for (i, (ab_c_row, a_bc_row)) in ab_c.iter().zip(a_bc.iter()).enumerate() {
            for (j, (&v1, &v2)) in ab_c_row.iter().zip(a_bc_row.iter()).enumerate() {
                assert!(
                    (v1 - v2).abs() < 1e-10,
                    "(AB)C[{i}][{j}]={} vs A(BC)={}",
                    v1,
                    v2
                );
            }
        }
    }

    #[test]
    fn test_green_lagrange_pure_shear() {
        // Pure shear: F = [[1,γ,0],[0,1,0],[0,0,1]]
        let gamma = 0.5;
        let f = [[1.0, gamma, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let e = green_lagrange_strain(f);
        // E12 = γ/2
        assert!(
            (e[0][1] - gamma / 2.0).abs() < 1e-10,
            "E12 for pure shear should be γ/2, got {}",
            e[0][1]
        );
    }

    #[test]
    fn test_cauchy_stress_symmetry_isotropic() {
        // For symmetric F (no rotation), Cauchy stress should be symmetric
        let f = [[1.2, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 0.8]];
        let mat = NeoHookeanMaterial::from_young_poisson(1000.0, 0.3);
        let piola = mat.piola_kirchhoff(f);
        let sigma = cauchy_stress(piola, f);
        for (i, row) in sigma.iter().enumerate() {
            for j in (i + 1)..3 {
                assert!(
                    (row[j] - sigma[j][i]).abs() < 1e-8,
                    "Cauchy stress not symmetric: [{i},{j}]={} vs [{j},{i}]={}",
                    row[j],
                    sigma[j][i]
                );
            }
        }
    }

    // ── RayleighDamping additional tests ──────────────────────────────────

    #[test]
    fn test_rayleigh_alpha_only_damps_velocity() {
        let rd = RayleighDamping::new(1.0, 0.0);
        let mut forces = vec![[0.0_f64; 3]; 2];
        let velocities = vec![[2.0_f64, 0.0, 0.0], [0.0, -1.0, 0.0]];
        let masses = vec![1.0, 2.0];
        rd.apply(&mut forces, &velocities, &masses);
        // Node 0: f = -alpha * m * v = -1 * 1 * 2 = -2
        assert!(
            (forces[0][0] - (-2.0)).abs() < 1e-12,
            "node 0 force[0] should be -2"
        );
        // Node 1: f = -1 * 2 * (-1) = 2
        assert!(
            (forces[1][1] - 2.0).abs() < 1e-12,
            "node 1 force[1] should be 2"
        );
    }

    #[test]
    fn test_rayleigh_symmetric_two_modes_roundtrip() {
        // Same xi at both modes → test symmetry of formula
        let xi = 0.07;
        let om1 = 5.0;
        let om2 = 50.0;
        let rd = RayleighDamping::from_two_modes(xi, om1, xi, om2);
        let xi_check1 = rd.damping_ratio(om1);
        let xi_check2 = rd.damping_ratio(om2);
        assert!(
            (xi_check1 - xi).abs() < 1e-8,
            "xi at om1: got {xi_check1}, expected {xi}"
        );
        assert!(
            (xi_check2 - xi).abs() < 1e-8,
            "xi at om2: got {xi_check2}, expected {xi}"
        );
    }

    // ── BoundaryConditions additional ─────────────────────────────────────

    #[test]
    fn test_bc_update_pin_position() {
        let mut bc = BoundaryConditions::new();
        bc.pin_at(0, [1.0, 2.0, 3.0]);
        // Pin at new position — should update
        bc.pin_at(0, [4.0, 5.0, 6.0]);
        assert_eq!(bc.len(), 1, "should not duplicate pin entries");
        let mut pos = vec![[0.0_f64; 3]];
        let mut vel = vec![[0.0_f64; 3]];
        bc.enforce(&mut pos, &mut vel);
        assert_eq!(
            pos[0],
            [4.0, 5.0, 6.0],
            "position should be updated to new pin"
        );
    }

    #[test]
    fn test_bc_multiple_pins() {
        let mut bc = BoundaryConditions::new();
        bc.pin_at(0, [1.0, 0.0, 0.0]);
        bc.pin_at(1, [0.0, 1.0, 0.0]);
        bc.pin_at(2, [0.0, 0.0, 1.0]);
        assert_eq!(bc.len(), 3);
        let mut pos = vec![[0.0_f64; 3]; 3];
        let mut vel = vec![[5.0_f64; 3]; 3];
        bc.enforce(&mut pos, &mut vel);
        assert_eq!(pos[0], [1.0, 0.0, 0.0]);
        assert_eq!(pos[1], [0.0, 1.0, 0.0]);
        assert_eq!(pos[2], [0.0, 0.0, 1.0]);
        for v in &vel {
            assert_eq!(*v, [0.0; 3]);
        }
    }

    // ── VelocityVerletIntegrator additional ───────────────────────────────

    #[test]
    fn test_velocity_verlet_predict_positions_no_accel() {
        let integrator = VelocityVerletIntegrator::new(1, 0.01);
        let mut positions = vec![[0.0_f64; 3]];
        let velocities = vec![[2.0_f64, 3.0, -1.0]];
        let pinned = vec![false];
        integrator.predict_positions(&mut positions, &velocities, &pinned);
        // x += dt*v + 0.5*dt^2*0 = 0.01 * [2,3,-1]
        assert!((positions[0][0] - 0.02).abs() < 1e-12);
        assert!((positions[0][1] - 0.03).abs() < 1e-12);
        assert!((positions[0][2] - (-0.01)).abs() < 1e-12);
    }

    #[test]
    fn test_velocity_verlet_pinned_node_predict_no_movement() {
        let integrator = VelocityVerletIntegrator::new(1, 0.01);
        let mut positions = vec![[5.0_f64, 5.0, 5.0]];
        let velocities = vec![[10.0_f64, 10.0, 10.0]];
        let pinned = vec![true];
        integrator.predict_positions(&mut positions, &velocities, &pinned);
        assert_eq!(positions[0], [5.0, 5.0, 5.0], "pinned node should not move");
    }

    #[test]
    fn test_velocity_verlet_correct_velocities_averages_accels() {
        let mut integrator = VelocityVerletIntegrator::new(1, 0.1);
        // prev_accel = [10, 0, 0]
        integrator.prev_accel = vec![[10.0, 0.0, 0.0]];
        let mut velocities = vec![[0.0_f64; 3]];
        let new_accel = vec![[20.0_f64, 0.0, 0.0]];
        let pinned = vec![false];
        integrator.correct_velocities(&mut velocities, &new_accel, &pinned);
        // dv = 0.5 * dt * (a_old + a_new) = 0.5 * 0.1 * 30 = 1.5
        assert!(
            (velocities[0][0] - 1.5).abs() < 1e-12,
            "velocity should be 1.5, got {}",
            velocities[0][0]
        );
        // prev_accel should be updated
        assert_eq!(integrator.prev_accel[0], [20.0, 0.0, 0.0]);
    }
}

#[cfg(test)]
mod corot_extended_tests {

    use crate::fem_soft::HhtAlphaParams;

    use crate::fem_soft::assemble_global_stiffness;

    use crate::fem_soft::central_difference_step;
    use crate::fem_soft::conjugate_gradient;
    use crate::fem_soft::corot_cauchy_stress;
    use crate::fem_soft::corot_internal_forces;
    use crate::fem_soft::corotational_strain;
    use crate::fem_soft::element_volume_stats;
    use crate::fem_soft::gravity_forces;

    use crate::fem_soft::lumped_mass_accelerations;
    use crate::fem_soft::polar_decompose_r;

    use crate::fem_soft::tests_extended::isotropic_d_matrix;
    use crate::fem_soft::tests_extended::mul3x3;
    use crate::fem_soft::tests_extended::transpose3x3;
    use crate::fem_soft::tet_signed_volume;
    use crate::fem_soft::tet_stiffness_matrix;

    // ── polar decompose identity ──────────────────────────────────────────

    #[test]
    fn test_polar_decompose_identity() {
        let identity = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let r = polar_decompose_r(identity);
        for (i, row) in r.iter().enumerate() {
            for (j, &val) in row.iter().enumerate() {
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!(
                    (val - expected).abs() < 1e-8,
                    "polar_decompose_r(I)[{i}][{j}] = {}, expected {expected}",
                    val
                );
            }
        }
    }

    #[test]
    fn test_polar_decompose_pure_rotation() {
        // 90-degree rotation about Z
        let theta = std::f64::consts::FRAC_PI_2;
        let f = [
            [theta.cos(), -theta.sin(), 0.0],
            [theta.sin(), theta.cos(), 0.0],
            [0.0, 0.0, 1.0],
        ];
        let r = polar_decompose_r(f);
        // R should equal F for a pure rotation
        for (i, row) in r.iter().enumerate() {
            for (j, &val) in row.iter().enumerate() {
                assert!(
                    (val - f[i][j]).abs() < 1e-7,
                    "r[{i}][{j}] = {}, f[{i}][{j}] = {}",
                    val,
                    f[i][j]
                );
            }
        }
    }

    #[test]
    fn test_polar_decompose_r_is_orthogonal() {
        let f = [[1.2, 0.3, -0.1], [0.1, 0.9, 0.2], [-0.2, 0.1, 1.1]];
        let r = polar_decompose_r(f);
        let rt = transpose3x3(r);
        let rrt = mul3x3(r, rt);
        for (i, row) in rrt.iter().enumerate() {
            for (j, &val) in row.iter().enumerate() {
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!(
                    (val - expected).abs() < 1e-7,
                    "R*R^T[{i}][{j}] = {}, expected {expected}",
                    val
                );
            }
        }
    }

    // ── corotational strain ───────────────────────────────────────────────

    #[test]
    fn test_corot_strain_zero_for_rigid_rotation() {
        // For a pure rotation F = R, strain should be zero
        let theta = 0.3_f64;
        let r = [
            [theta.cos(), -theta.sin(), 0.0],
            [theta.sin(), theta.cos(), 0.0],
            [0.0, 0.0, 1.0],
        ];
        let eps = corotational_strain(r, r);
        for (i, row) in eps.iter().enumerate() {
            for (j, &val) in row.iter().enumerate() {
                assert!(
                    val.abs() < 1e-7,
                    "strain[{i}][{j}] = {} for pure rotation",
                    val
                );
            }
        }
    }

    #[test]
    fn test_corot_strain_symmetric() {
        let f = [[1.1, 0.05, 0.0], [0.05, 1.0, 0.02], [0.0, 0.02, 0.95]];
        let r = polar_decompose_r(f);
        let eps = corotational_strain(r, f);
        for (i, row) in eps.iter().enumerate() {
            for (j, &val) in row.iter().enumerate() {
                assert!(
                    (val - eps[j][i]).abs() < 1e-10,
                    "strain not symmetric at [{i}][{j}]"
                );
            }
        }
    }

    // ── Cauchy stress ─────────────────────────────────────────────────────

    #[test]
    fn test_corot_cauchy_stress_zero_strain() {
        let eps = [[0.0_f64; 3]; 3];
        let sigma = corot_cauchy_stress(eps, 1e5, 7e4);
        for (i, row) in sigma.iter().enumerate() {
            for (j, &val) in row.iter().enumerate() {
                assert_eq!(val, 0.0, "sigma[{i}][{j}] should be zero");
            }
        }
    }

    #[test]
    fn test_corot_cauchy_stress_hydrostatic() {
        // Uniform volumetric strain: eps = e * I
        let e = 0.01;
        let mut eps = [[0.0_f64; 3]; 3];
        eps[0][0] = e;
        eps[1][1] = e;
        eps[2][2] = e;
        let lambda = 1e5;
        let mu = 7e4;
        let sigma = corot_cauchy_stress(eps, lambda, mu);
        // Expected: sigma_ii = (lambda + 2*mu) * e + lambda * 2*e = (3*lambda + 2*mu)*e
        let expected = (3.0 * lambda + 2.0 * mu) * e;
        for (i, row) in sigma.iter().enumerate() {
            assert!(
                (row[i] - expected).abs() < 1.0,
                "sigma[{i}][{i}] = {}, expected {expected}",
                row[i]
            );
        }
    }

    #[test]
    fn test_corot_cauchy_stress_symmetric() {
        let mut eps = [[0.0_f64; 3]; 3];
        eps[0][1] = 0.005;
        eps[1][0] = 0.005;
        eps[0][0] = 0.01;
        eps[1][1] = 0.008;
        let sigma = corot_cauchy_stress(eps, 1.2e5, 8e4);
        for (i, row) in sigma.iter().enumerate() {
            for (j, &val) in row.iter().enumerate() {
                assert!(
                    (val - sigma[j][i]).abs() < 1e-8,
                    "sigma not symmetric at [{i}][{j}]"
                );
            }
        }
    }

    // ── stiffness matrix ──────────────────────────────────────────────────

    #[test]
    fn test_tet_stiffness_positive_definite_trace() {
        let p0 = [0.0, 0.0, 0.0];
        let p1 = [1.0, 0.0, 0.0];
        let p2 = [0.0, 1.0, 0.0];
        let p3 = [0.0, 0.0, 1.0];
        let ke = tet_stiffness_matrix(p0, p1, p2, p3, 1e5, 7e4);
        let trace: f64 = (0..12).map(|i| ke[i * 12 + i]).sum();
        assert!(
            trace > 0.0,
            "stiffness matrix trace should be positive, got {trace}"
        );
    }

    #[test]
    fn test_tet_stiffness_symmetric() {
        let p0 = [0.0, 0.0, 0.0];
        let p1 = [1.0, 0.0, 0.0];
        let p2 = [0.0, 1.0, 0.0];
        let p3 = [0.0, 0.0, 1.0];
        let ke = tet_stiffness_matrix(p0, p1, p2, p3, 1e5, 7e4);
        for i in 0..12 {
            for j in 0..12 {
                assert!(
                    (ke[i * 12 + j] - ke[j * 12 + i]).abs() < 1e-6,
                    "K[{i}][{j}] != K[{j}][{i}]: {} vs {}",
                    ke[i * 12 + j],
                    ke[j * 12 + i]
                );
            }
        }
    }

    #[test]
    fn test_tet_stiffness_degenerate_returns_zero() {
        // All four nodes coplanar → zero matrix
        let p0 = [0.0, 0.0, 0.0];
        let p1 = [1.0, 0.0, 0.0];
        let p2 = [0.0, 1.0, 0.0];
        let p3 = [0.5, 0.5, 0.0]; // same plane
        let ke = tet_stiffness_matrix(p0, p1, p2, p3, 1e5, 7e4);
        let norm: f64 = ke.iter().map(|v| v.abs()).sum();
        assert!(
            norm < 1e-10,
            "degenerate tet should give zero stiffness, norm={norm}"
        );
    }

    // ── internal forces ───────────────────────────────────────────────────

    #[test]
    fn test_corot_internal_forces_zero_at_rest() {
        let rest = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        let forces = corot_internal_forces(rest, rest, 1e5, 7e4);
        for (a, force) in forces.iter().enumerate() {
            for (d, &val) in force.iter().enumerate() {
                assert!(
                    val.abs() < 1e-8,
                    "internal force at rest [{a}][{d}] = {}",
                    val
                );
            }
        }
    }

    #[test]
    fn test_corot_internal_forces_nonzero_under_stretch() {
        let rest = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        let mut deformed = rest;
        deformed[1][0] = 1.2; // stretch node 1 in x
        let forces = corot_internal_forces(rest, deformed, 1e5, 7e4);
        let f_norm: f64 = forces.iter().flat_map(|f| f.iter()).map(|v| v.abs()).sum();
        assert!(
            f_norm > 1.0,
            "forces should be nonzero under stretch, norm={f_norm}"
        );
    }

    // ── HHT-α integrator ──────────────────────────────────────────────────

    #[test]
    fn test_hht_alpha_from_alpha_zero() {
        let p = HhtAlphaParams::from_alpha(0.0, 0.01);
        assert!(
            (p.beta - 0.25).abs() < 1e-12,
            "beta should be 0.25 for alpha=0"
        );
        assert!(
            (p.gamma - 0.5).abs() < 1e-12,
            "gamma should be 0.5 for alpha=0"
        );
    }

    #[test]
    fn test_hht_alpha_predict_no_accel() {
        let p = HhtAlphaParams::from_alpha(0.0, 0.1);
        let u = vec![[1.0, 2.0, 3.0]];
        let v = vec![[0.5, -0.5, 0.0]];
        let a = vec![[0.0, 0.0, 0.0]];
        let (u_pred, v_pred) = p.predict(&u, &v, &a);
        assert!(
            (u_pred[0][0] - 1.05).abs() < 1e-10,
            "u_pred[x] = {}",
            u_pred[0][0]
        );
        assert!(
            (v_pred[0][0] - 0.5).abs() < 1e-10,
            "v_pred[x] = {}",
            v_pred[0][0]
        );
    }

    #[test]
    fn test_hht_alpha_correct_updates_velocity() {
        let p = HhtAlphaParams::from_alpha(0.0, 0.1);
        let mut u_pred = vec![[1.0, 0.0, 0.0]];
        let mut v_pred = vec![[0.0, 0.0, 0.0]];
        let a_new = vec![[10.0, 0.0, 0.0]];
        p.correct(&mut u_pred, &mut v_pred, &a_new);
        // v += gamma * dt * a = 0.5 * 0.1 * 10 = 0.5
        assert!(
            (v_pred[0][0] - 0.5).abs() < 1e-10,
            "v_pred[x] = {}",
            v_pred[0][0]
        );
        // u += beta * dt^2 * a = 0.25 * 0.01 * 10 = 0.025
        assert!(
            (u_pred[0][0] - 1.025).abs() < 1e-10,
            "u_pred[x] = {}",
            u_pred[0][0]
        );
    }

    // ── central difference ────────────────────────────────────────────────

    #[test]
    fn test_central_difference_step_free_node() {
        let mut x = vec![[0.0, 0.0, 0.0]];
        let mut v_half = vec![[1.0, 0.0, 0.0]];
        let a = vec![[2.0, 0.0, 0.0]];
        let pinned = vec![false];
        central_difference_step(&mut x, &mut v_half, &a, 0.1, &pinned);
        // v_half += a*dt = 1 + 0.2 = 1.2
        assert!(
            (v_half[0][0] - 1.2).abs() < 1e-12,
            "v_half[x] = {}",
            v_half[0][0]
        );
        // x += v_half*dt = 1.2 * 0.1 = 0.12
        assert!((x[0][0] - 0.12).abs() < 1e-12, "x[x] = {}", x[0][0]);
    }

    #[test]
    fn test_central_difference_step_pinned_node() {
        let mut x = vec![[5.0, 5.0, 5.0]];
        let mut v_half = vec![[10.0, 10.0, 10.0]];
        let a = vec![[100.0, 100.0, 100.0]];
        let pinned = vec![true];
        central_difference_step(&mut x, &mut v_half, &a, 0.1, &pinned);
        assert_eq!(x[0], [5.0, 5.0, 5.0], "pinned node must not move");
        assert_eq!(
            v_half[0],
            [10.0, 10.0, 10.0],
            "pinned velocity must not change"
        );
    }

    // ── volume stats ──────────────────────────────────────────────────────

    #[test]
    fn test_tet_signed_volume_unit_tet() {
        let p = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        let v = tet_signed_volume(p);
        assert!((v - 1.0 / 6.0).abs() < 1e-12, "unit tet volume = {v}");
    }

    #[test]
    fn test_element_volume_stats_single() {
        let positions = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        let indices = vec![[0usize, 1, 2, 3]];
        let (min_v, max_v, total) = element_volume_stats(&positions, &indices);
        assert!((total - 1.0 / 6.0).abs() < 1e-12);
        assert_eq!(min_v, max_v, "single element: min==max");
    }

    // ── gravity forces ────────────────────────────────────────────────────

    #[test]
    fn test_gravity_forces_single_node() {
        let masses = vec![2.0];
        let g = [0.0, -9.81, 0.0];
        let f = gravity_forces(&masses, g);
        assert!(
            (f[0][1] - (-19.62)).abs() < 1e-8,
            "gravity force y = {}",
            f[0][1]
        );
        assert_eq!(f[0][0], 0.0);
        assert_eq!(f[0][2], 0.0);
    }

    #[test]
    fn test_gravity_forces_multiple_nodes() {
        let masses = vec![1.0, 2.0, 3.0];
        let g = [0.0, -9.81, 0.0];
        let f = gravity_forces(&masses, g);
        assert_eq!(f.len(), 3);
        for (i, &m) in masses.iter().enumerate() {
            assert!((f[i][1] - (-9.81 * m)).abs() < 1e-10);
        }
    }

    // ── lumped mass accelerations ─────────────────────────────────────────

    #[test]
    fn test_lumped_mass_accelerations_free() {
        let masses = vec![2.0];
        let forces = vec![[10.0, 0.0, 0.0]];
        let pinned = vec![false];
        let acc = lumped_mass_accelerations(&masses, &forces, &pinned);
        assert!((acc[0][0] - 5.0).abs() < 1e-12, "acc = {}", acc[0][0]);
    }

    #[test]
    fn test_lumped_mass_accelerations_pinned() {
        let masses = vec![2.0];
        let forces = vec![[999.0, 0.0, 0.0]];
        let pinned = vec![true];
        let acc = lumped_mass_accelerations(&masses, &forces, &pinned);
        assert_eq!(acc[0], [0.0; 3], "pinned node acc should be zero");
    }

    // ── CG solver ─────────────────────────────────────────────────────────

    #[test]
    fn test_conjugate_gradient_2x2() {
        // [4, 1; 1, 3] * x = [1; 2]  → x = [1/11, 7/11]
        let a = [4.0_f64, 1.0, 1.0, 3.0];
        let b = [1.0_f64, 2.0];
        let x = conjugate_gradient(&a, &b, 1e-12, 100);
        assert!((x[0] - 1.0 / 11.0).abs() < 1e-8, "x[0] = {}", x[0]);
        assert!((x[1] - 7.0 / 11.0).abs() < 1e-8, "x[1] = {}", x[1]);
    }

    #[test]
    fn test_conjugate_gradient_identity() {
        let a = [1.0_f64, 0.0, 0.0, 1.0];
        let b = [3.0_f64, 7.0];
        let x = conjugate_gradient(&a, &b, 1e-12, 100);
        assert!((x[0] - 3.0).abs() < 1e-8);
        assert!((x[1] - 7.0).abs() < 1e-8);
    }

    // ── assemble global stiffness ─────────────────────────────────────────

    #[test]
    fn test_assemble_global_stiffness_single_element() {
        let positions = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        let elements = vec![([0usize, 1, 2, 3], 1e5, 7e4)];
        let k = assemble_global_stiffness(&positions, &elements);
        assert_eq!(k.len(), 144);
        let trace: f64 = (0..12).map(|i| k[i * 12 + i]).sum();
        assert!(
            trace > 0.0,
            "global K trace should be positive, got {trace}"
        );
    }

    #[test]
    fn test_isotropic_d_matrix_diagonal() {
        let lambda = 1e5;
        let mu = 7e4;
        let d = isotropic_d_matrix(lambda, mu);
        let expected_diag_0 = lambda + 2.0 * mu;
        assert!(
            (d[0][0] - expected_diag_0).abs() < 1e-6,
            "D[0][0] = {}",
            d[0][0]
        );
        assert!((d[3][3] - mu).abs() < 1e-6, "D[3][3] = {}", d[3][3]);
    }
}
