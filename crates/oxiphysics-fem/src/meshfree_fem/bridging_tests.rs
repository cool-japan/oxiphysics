// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Tests of the bridging-domain coupling: pairing, energy shares, symmetric
//! blending and the penalty ties.

use super::functions::{Point2, generate_regular_nodes, mls_shape_2d};
use super::types::{
    CouplingInterface, CouplingMethod, CouplingNode, KernelType, MeshfreeNode, SupportDomain,
};
use crate::Error;

/// Meshfree nodes on the line `y = 0.5` at the given abscissae.
fn line_cloud(xs: &[f64], support_radius: f64) -> SupportDomain {
    let nodes = xs
        .iter()
        .map(|&x| MeshfreeNode {
            pos: [x, 0.5],
            support_radius,
            value: 0.0,
        })
        .collect();
    SupportDomain::new(nodes, 1.0)
}

fn count_nonzero(k: &[f64]) -> usize {
    k.iter().filter(|v| **v != 0.0).count()
}

/// A symmetric `n x n` matrix with all entries non-zero.
fn test_matrix(n: usize, seed: f64) -> Vec<f64> {
    let mut k = vec![0.0; n * n];
    for i in 0..n {
        for j in i..n {
            let v = 1.0 + seed + 0.37 * i as f64 + 0.11 * j as f64 + 0.05 * (i * j) as f64;
            k[i * n + j] = v;
            k[j * n + i] = v;
        }
    }
    k
}

#[test]
fn three_coincident_pairs_give_twenty_four_penalty_entries() {
    // Overlap strip [0.1, 0.9]; meshfree nodes at 0.2, 0.5, 0.8 (+ two outside the strip);
    // FEM nodes coincide with the first three and one extra node sits at 0.65.
    let domain = line_cloud(&[0.0, 0.2, 0.5, 0.8, 1.0], 0.2);
    let fem_nodes: Vec<Point2> = vec![[0.0, 0.5], [0.2, 0.5], [0.5, 0.5], [0.8, 0.5], [0.65, 0.5]];
    let interface = super::setup_bridging_domain(&domain, &fem_nodes, 0.1, 0.9);

    let pairs: Vec<&CouplingNode> = interface
        .nodes
        .iter()
        .filter(|n| n.meshfree_idx.is_some() && n.fem_idx.is_some())
        .collect();
    assert_eq!(pairs.len(), 3, "{:?}", interface.nodes);
    assert_eq!(
        pairs
            .iter()
            .map(|n| (n.meshfree_idx, n.fem_idx))
            .collect::<Vec<_>>(),
        vec![(Some(1), Some(1)), (Some(2), Some(2)), (Some(3), Some(3))]
    );
    // The extra FEM node has no meshfree partner but stays in the interface.
    let unpaired: Vec<&CouplingNode> = interface
        .nodes
        .iter()
        .filter(|n| n.meshfree_idx.is_none())
        .collect();
    assert_eq!(unpaired.len(), 1);
    assert_eq!(unpaired[0].fem_idx, Some(4));
    // Nodes outside the strip are not part of the interface.
    assert_eq!(interface.num_coupling_nodes(), 4);

    // Each pair ties two components with four entries each: 3 pairs x 2 DOFs x 4 = 24.
    let k = interface.compute_penalty_coupling(10, 10);
    assert_eq!(k.len(), 20 * 20);
    assert_eq!(count_nonzero(&k), 24);
    // The coupling matrix of a rigid translation of both models is zero.
    let ones: Vec<f64> = vec![1.0; 20];
    for row in 0..20 {
        let f: f64 = (0..20).map(|col| k[row * 20 + col] * ones[col]).sum();
        assert!(f.abs() < 1e-9, "row {row}: {f}");
    }
}

#[test]
fn pairing_is_one_to_one_and_prefers_the_nearest_node() {
    // Two meshfree nodes compete for the same FEM node: the nearer one wins and the
    // other falls back to the second FEM node, or stays unpaired.
    let domain = line_cloud(&[0.30, 0.34], 1.0);
    let fem_nodes: Vec<Point2> = vec![[0.31, 0.5], [0.70, 0.5]];
    let interface = super::setup_bridging_domain(&domain, &fem_nodes, 0.0, 1.0);
    let by_meshfree = |i: usize| {
        interface
            .nodes
            .iter()
            .find(|n| n.meshfree_idx == Some(i))
            .map(|n| n.fem_idx)
    };
    // Node 0 (distance 0.01) wins FEM node 0; node 1 then takes FEM node 1, which lies at
    // distance 0.36 < 0.5 (half the support radius).
    assert_eq!(by_meshfree(0), Some(Some(0)));
    assert_eq!(by_meshfree(1), Some(Some(1)));
    // No FEM node is used twice.
    let mut used: Vec<usize> = interface.nodes.iter().filter_map(|n| n.fem_idx).collect();
    used.sort_unstable();
    used.dedup();
    assert_eq!(
        used.len(),
        interface
            .nodes
            .iter()
            .filter(|n| n.fem_idx.is_some())
            .count()
    );
    // A node farther than half the support radius from every FEM node stays unpaired.
    let far = super::setup_bridging_domain(&line_cloud(&[0.3], 0.2), &[[0.45, 0.5]], 0.0, 1.0);
    assert_eq!(far.nodes.len(), 2);
    assert!(
        far.nodes
            .iter()
            .all(|n| n.meshfree_idx.is_none() || n.fem_idx.is_none())
    );
}

#[test]
fn blend_weight_is_the_meshfree_share_for_both_kinds_of_nodes() {
    let domain = line_cloud(&[0.2, 0.5, 0.8], 0.2);
    let fem_nodes: Vec<Point2> = vec![[0.2, 0.5], [0.65, 0.5]];
    let interface = super::setup_bridging_domain(&domain, &fem_nodes, 0.1, 0.9);
    for node in &interface.nodes {
        let expected = ((node.pos[0] - 0.1) / 0.8).clamp(0.0, 1.0);
        assert!(
            (node.blend_weight - expected).abs() < 1e-14,
            "node at {}: {} vs {expected}",
            node.pos[0],
            node.blend_weight
        );
    }
    // The FEM-only node at 0.65 is 0.6875 meshfree (the old code stored 1 − w = 0.3125).
    let fem_only = interface
        .nodes
        .iter()
        .find(|n| n.meshfree_idx.is_none())
        .expect("unpaired FEM node");
    assert!((fem_only.blend_weight - 0.6875).abs() < 1e-14);
}

#[test]
fn compute_bridging_weights_follows_position_and_stored_origin() {
    let mut interface = CouplingInterface::new(CouplingMethod::BridgingDomain, 1.0e6);
    interface.overlap_x_min = 1.0;
    for (i, &x) in [0.5, 1.0, 1.25, 1.5, 2.0, 3.0].iter().enumerate() {
        interface.add_node(CouplingNode {
            pos: [x, 0.0],
            meshfree_idx: Some(i),
            fem_idx: Some(i),
            blend_weight: 0.5,
        });
    }
    interface.compute_bridging_weights(1.0);
    let w: Vec<f64> = interface.nodes.iter().map(|n| n.blend_weight).collect();
    assert_eq!(w, vec![0.0, 0.0, 0.25, 0.5, 1.0, 1.0]);
    // A different width rescales the ramp; a non-positive width is a sharp step at x_min.
    interface.compute_bridging_weights(4.0);
    let w: Vec<f64> = interface.nodes.iter().map(|n| n.blend_weight).collect();
    assert_eq!(w, vec![0.0, 0.0, 0.0625, 0.125, 0.25, 0.5]);
    interface.compute_bridging_weights(0.0);
    let w: Vec<f64> = interface.nodes.iter().map(|n| n.blend_weight).collect();
    assert_eq!(w, vec![0.0, 0.0, 1.0, 1.0, 1.0, 1.0]);
}

#[test]
fn blended_matrices_add_up_to_the_original_stiffness() {
    // Identical K in both models on coincident nodes: the two energy shares are
    // complementary, so K_mf + K_fem = K and the energy is neither lost nor doubled.
    let mut interface = CouplingInterface::new(CouplingMethod::BridgingDomain, 1.0e6);
    interface.overlap_x_min = 0.1;
    for i in 0..5 {
        interface.add_node(CouplingNode {
            pos: [0.25 * i as f64, 0.5],
            meshfree_idx: Some(i),
            fem_idx: Some(i),
            blend_weight: 0.0,
        });
    }
    interface.compute_bridging_weights(1.0);
    let n = 10;
    let k0 = test_matrix(n, 0.0);
    let mut k_mf = k0.clone();
    let mut k_fem = k0.clone();
    interface
        .apply_bridging_blend(&mut k_mf, n, &mut k_fem, n)
        .expect("valid interface");
    let mut changed = false;
    for idx in 0..n * n {
        assert!(
            (k_mf[idx] + k_fem[idx] - k0[idx]).abs() < 1e-12 * k0[idx].abs(),
            "entry {idx}: {} + {} != {}",
            k_mf[idx],
            k_fem[idx],
            k0[idx]
        );
        changed |= (k_mf[idx] - k0[idx]).abs() > 1e-9;
    }
    assert!(
        changed,
        "a non-trivial blend must change the meshfree matrix"
    );
    // Both blended matrices stay symmetric (scaling only the rows would not).
    for i in 0..n {
        for j in 0..n {
            assert!((k_mf[i * n + j] - k_mf[j * n + i]).abs() < 1e-14);
            assert!((k_fem[i * n + j] - k_fem[j * n + i]).abs() < 1e-14);
        }
    }
}

#[test]
fn uniform_share_scales_the_whole_block_and_other_nodes_are_untouched() {
    // Interface nodes 1 and 2 of 4 carry the share 0.3; the entries between nodes
    // outside the interface keep their value.
    let mut interface = CouplingInterface::new(CouplingMethod::BridgingDomain, 1.0e6);
    for i in [1usize, 2] {
        interface.add_node(CouplingNode {
            pos: [i as f64, 0.0],
            meshfree_idx: Some(i),
            fem_idx: Some(i),
            blend_weight: 0.3,
        });
    }
    let n = 8;
    let k0 = test_matrix(n, 1.0);
    let mut k_mf = k0.clone();
    let mut k_fem = k0.clone();
    interface
        .apply_bridging_blend(&mut k_mf, n, &mut k_fem, n)
        .expect("valid interface");
    for row in 0..n {
        for col in 0..n {
            let in_mf_block = (2..6).contains(&row) && (2..6).contains(&col);
            let outside = !(2..6).contains(&row) && !(2..6).contains(&col);
            let idx = row * n + col;
            if in_mf_block {
                assert!(
                    (k_mf[idx] - 0.3 * k0[idx]).abs() < 1e-13,
                    "mf ({row},{col})"
                );
                assert!(
                    (k_fem[idx] - 0.7 * k0[idx]).abs() < 1e-13,
                    "fem ({row},{col})"
                );
            } else if outside {
                assert_eq!(k_mf[idx], k0[idx], "mf ({row},{col})");
                assert_eq!(k_fem[idx], k0[idx], "fem ({row},{col})");
            } else {
                // Mixed block: the average of the node shares (1 and 0.3 / 1 and 0.7).
                assert!(
                    (k_mf[idx] - 0.65 * k0[idx]).abs() < 1e-13,
                    "mf ({row},{col})"
                );
                assert!(
                    (k_fem[idx] - 0.85 * k0[idx]).abs() < 1e-13,
                    "fem ({row},{col})"
                );
            }
        }
    }
}

#[test]
fn blend_rejects_inconsistent_interfaces_and_leaves_the_matrices_alone() {
    let n = 6;
    let k0 = test_matrix(n, 0.5);
    let node = |mf: Option<usize>, fem: Option<usize>, w: f64| CouplingNode {
        pos: [0.0, 0.0],
        meshfree_idx: mf,
        fem_idx: fem,
        blend_weight: w,
    };
    let cases: Vec<(&str, Vec<CouplingNode>)> = vec![
        ("meshfree node out of range", vec![node(Some(3), None, 0.5)]),
        ("FEM node out of range", vec![node(None, Some(3), 0.5)]),
        (
            "meshfree node twice",
            vec![node(Some(1), None, 0.5), node(Some(1), Some(0), 0.2)],
        ),
        (
            "FEM node twice",
            vec![node(None, Some(2), 0.5), node(Some(0), Some(2), 0.2)],
        ),
        ("weight above one", vec![node(Some(0), Some(0), 1.5)]),
        ("weight below zero", vec![node(Some(0), Some(0), -0.1)]),
        ("NaN weight", vec![node(Some(0), Some(0), f64::NAN)]),
    ];
    for (label, nodes) in cases {
        let mut interface = CouplingInterface::new(CouplingMethod::BridgingDomain, 1.0);
        for cn in nodes {
            interface.add_node(cn);
        }
        let (mut a, mut b) = (k0.clone(), k0.clone());
        let result = interface.apply_bridging_blend(&mut a, n, &mut b, n);
        assert!(
            matches!(result, Err(Error::InvalidInput(_))),
            "{label}: {result:?}"
        );
        assert_eq!(a, k0, "{label}: meshfree matrix modified");
        assert_eq!(b, k0, "{label}: FEM matrix modified");
    }
    // Wrong sizes and odd DOF counts.
    let interface = CouplingInterface::new(CouplingMethod::BridgingDomain, 1.0);
    let (mut a, mut b) = (k0.clone(), k0.clone());
    assert!(
        interface
            .apply_bridging_blend(&mut a, n, &mut b[..10], n)
            .is_err()
    );
    assert!(
        interface
            .apply_bridging_blend(&mut a[..35], n, &mut b, n)
            .is_err()
    );
    let mut odd = vec![0.0; 25];
    let mut odd2 = vec![0.0; 25];
    assert!(
        interface
            .apply_bridging_blend(&mut odd, 5, &mut odd2, 5)
            .is_err()
    );
    // An empty interface changes nothing.
    interface
        .apply_bridging_blend(&mut a, n, &mut b, n)
        .expect("empty interface is valid");
    assert_eq!(a, k0);
}

#[test]
fn bridging_interface_of_a_setup_is_accepted_by_the_blend() {
    // End to end: setup -> blend with the model sizes of the two clouds.
    let nodes = generate_regular_nodes(0.0, 0.0, 1.0, 1.0, 6, 3, 1.5);
    let domain = SupportDomain::new(nodes, 1.0);
    let fem_nodes: Vec<Point2> = (0..4)
        .flat_map(|i| (0..3).map(move |j| [0.4 + 0.2 * i as f64, 0.5 * j as f64]))
        .collect();
    let interface = super::setup_bridging_domain(&domain, &fem_nodes, 0.4, 0.8);
    assert!(
        interface
            .nodes
            .iter()
            .any(|n| n.meshfree_idx.is_some() && n.fem_idx.is_some())
    );
    let n_mf = 2 * domain.node_count();
    let n_fem = 2 * fem_nodes.len();
    let mut k_mf = test_matrix(n_mf, 0.0);
    let mut k_fem = test_matrix(n_fem, 0.0);
    let (k_mf0, k_fem0) = (k_mf.clone(), k_fem.clone());
    interface
        .apply_bridging_blend(&mut k_mf, n_mf, &mut k_fem, n_fem)
        .expect("interface built by setup_bridging_domain");
    assert_ne!(k_mf, k_mf0);
    assert_ne!(k_fem, k_fem0);
}

#[test]
fn field_penalty_ties_the_mls_field_not_the_nodal_parameter() {
    let nodes = generate_regular_nodes(0.0, 0.0, 1.0, 1.0, 5, 5, 2.5);
    let domain = SupportDomain::new(nodes, 1.0);
    let ktype = KernelType::CubicSpline;
    let n_mf_nodes = domain.node_count();
    let ndof_mf = 2 * n_mf_nodes;
    // One pair: meshfree node 12 (the centre) with FEM node 0.
    let centre = domain.nodes[12].pos;
    let mut interface = CouplingInterface::new(CouplingMethod::BridgingDomain, 1.0e4);
    interface.add_node(CouplingNode {
        pos: centre,
        meshfree_idx: Some(12),
        fem_idx: Some(0),
        blend_weight: 0.5,
    });
    let ndof_fem = 4;
    let ndof = ndof_mf + ndof_fem;
    let k_field = interface
        .compute_field_penalty_coupling(&domain, ktype, ndof_mf, ndof_fem)
        .expect("shape functions exist at the centre");
    let k_nodal = interface.compute_penalty_coupling(ndof_mf, ndof_fem);

    // Parameters û: an arbitrary non-nodal pattern; the FEM value equals the MLS field.
    let shape = mls_shape_2d(centre, &domain, ktype).expect("shape functions");
    let parameters: Vec<f64> = (0..ndof_mf)
        .map(|i| 0.01 * ((i * 7 % 11) as f64) - 0.03)
        .collect();
    let field_x: f64 = shape
        .neighbours
        .iter()
        .zip(shape.phi.iter())
        .map(|(&j, &phi)| phi * parameters[2 * j])
        .sum();
    let field_y: f64 = shape
        .neighbours
        .iter()
        .zip(shape.phi.iter())
        .map(|(&j, &phi)| phi * parameters[2 * j + 1])
        .sum();
    let mut state = parameters.clone();
    state.extend_from_slice(&[field_x, field_y, 0.0, 0.0]);
    let apply = |k: &[f64]| -> Vec<f64> {
        (0..ndof)
            .map(|r| (0..ndof).map(|c| k[r * ndof + c] * state[c]).sum())
            .collect()
    };
    let max_abs = |v: &[f64]| v.iter().fold(0.0_f64, |m, x| m.max(x.abs()));
    // The continuity condition u_h(x_c) = u_fem holds, so the field tie exerts no force ...
    assert!(
        max_abs(&apply(&k_field)) < 1e-9,
        "field tie force {}",
        max_abs(&apply(&k_field))
    );
    // ... whereas the nodal tie (û_12 − u_fem) is violated and pushes back with p · mismatch.
    let mismatch = (parameters[24] - field_x).abs();
    assert!(
        mismatch > 1e-4,
        "test setup: the parameter equals the field ({mismatch})"
    );
    assert!(max_abs(&apply(&k_nodal)) > 1.0e3 * mismatch);

    // The field tie is symmetric and positive semi-definite (p a aᵀ).
    for r in 0..ndof {
        for c in 0..ndof {
            assert!((k_field[r * ndof + c] - k_field[c * ndof + r]).abs() < 1e-9);
        }
    }
    let mut probe = parameters.clone();
    probe.extend_from_slice(&[0.2, -0.1, 0.0, 0.0]);
    let k_probe: Vec<f64> = (0..ndof)
        .map(|r| (0..ndof).map(|c| k_field[r * ndof + c] * probe[c]).sum())
        .collect();
    let energy: f64 = probe.iter().zip(k_probe.iter()).map(|(v, kv)| v * kv).sum();
    assert!(
        energy > 0.0,
        "a violated continuity condition stores penalty energy: {energy}"
    );
}

#[test]
fn field_penalty_validates_its_input() {
    let nodes = generate_regular_nodes(0.0, 0.0, 1.0, 1.0, 5, 5, 2.5);
    let domain = SupportDomain::new(nodes, 1.0);
    let mut interface = CouplingInterface::new(CouplingMethod::BridgingDomain, 1.0e4);
    interface.add_node(CouplingNode {
        pos: domain.nodes[12].pos,
        meshfree_idx: Some(12),
        fem_idx: Some(3),
        blend_weight: 0.5,
    });
    // Wrong meshfree DOF count.
    assert!(
        interface
            .compute_field_penalty_coupling(&domain, KernelType::CubicSpline, 10, 8)
            .is_err()
    );
    // FEM node 3 needs DOFs 6 and 7.
    assert!(
        interface
            .compute_field_penalty_coupling(&domain, KernelType::CubicSpline, 50, 6)
            .is_err()
    );
    // Coupling node outside the cloud: no MLS shape functions.
    let mut outside = CouplingInterface::new(CouplingMethod::BridgingDomain, 1.0e4);
    outside.add_node(CouplingNode {
        pos: [50.0, 50.0],
        meshfree_idx: Some(0),
        fem_idx: Some(0),
        blend_weight: 0.5,
    });
    assert!(matches!(
        outside.compute_field_penalty_coupling(&domain, KernelType::CubicSpline, 50, 2),
        Err(Error::InvalidInput(_))
    ));
    // Unpaired coupling nodes carry no tie.
    let mut unpaired = CouplingInterface::new(CouplingMethod::BridgingDomain, 1.0e4);
    unpaired.add_node(CouplingNode {
        pos: [0.5, 0.5],
        meshfree_idx: Some(12),
        fem_idx: None,
        blend_weight: 0.5,
    });
    let k = unpaired
        .compute_field_penalty_coupling(&domain, KernelType::CubicSpline, 50, 2)
        .expect("no pairs");
    assert!(k.iter().all(|v| *v == 0.0));
}

#[test]
fn setup_handles_empty_and_degenerate_overlaps() {
    let domain = line_cloud(&[0.2, 0.5], 0.2);
    // No node in the strip.
    let none = super::setup_bridging_domain(&domain, &[[0.3, 0.5]], 0.6, 0.9);
    assert_eq!(none.num_coupling_nodes(), 0);
    // Reversed bounds select nothing.
    let reversed = super::setup_bridging_domain(&domain, &[[0.3, 0.5]], 0.9, 0.1);
    assert_eq!(reversed.num_coupling_nodes(), 0);
    // A zero-width strip is a sharp interface; a node on its edge has the meshfree share 0.
    let sharp = super::setup_bridging_domain(&domain, &[[0.2, 0.5]], 0.2, 0.2);
    assert_eq!(sharp.num_coupling_nodes(), 1);
    assert_eq!(sharp.nodes[0].blend_weight, 0.0);
    assert_eq!(sharp.method, CouplingMethod::BridgingDomain);
    assert_eq!(sharp.overlap_x_min, 0.2);
}
