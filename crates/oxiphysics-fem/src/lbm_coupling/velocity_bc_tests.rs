// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Tests of the velocity boundary conditions that
//! [`FemLbmCoupling::transfer_velocity_to_lbm`] derives from the IB points.

use super::*;

// ------------------------------------------------------------------
// Velocity boundary conditions imposed on the LBM grid
// ------------------------------------------------------------------

/// Coupling whose IB points sit at `points` and move with `velocities`.
fn moving_points_coupling(
    points: &[[f64; 3]],
    velocities: &[[f64; 3]],
    spacing: f64,
) -> FemLbmCoupling {
    let mut coupling = FemLbmCoupling::new(FemLbmCouplingConfig::default());
    coupling
        .set_interface_nodes((0..points.len()).collect(), points.to_vec(), spacing)
        .expect("valid interface");
    coupling
        .update_fem_state(points, velocities)
        .expect("matching lengths");
    coupling
}

/// Kernel weight `phi(x/h - i) phi(y/h - j) phi(z/h - k)` of the node
/// `[i, j, k]` for a point at `pos`, from the public 1-D kernel.
fn node_weight(pos: [f64; 3], h: f64, node: [usize; 3]) -> f64 {
    (0..3)
        .map(|d| delta_function(pos[d] / h - node[d] as f64, 1.0))
        .product()
}

/// Velocity imposed at `node`, if the node is in the list.
fn bc_at(bcs: &LbmVelocityBcList, node: [usize; 3]) -> Option<[f64; 3]> {
    bcs.iter().find(|(n, _)| *n == node).map(|(_, v)| *v)
}

fn assert_vec_close(got: [f64; 3], want: [f64; 3], tol: f64, what: &str) {
    for (g, w) in got.iter().zip(want.iter()) {
        assert!(
            (g - w).abs() < tol,
            "{what}: got {got:?}, expected {want:?}"
        );
    }
}

/// Deterministic pseudo-random numbers in `[0, 1)` (64-bit LCG).
fn lcg_unit(state: &mut u64) -> f64 {
    *state = state
        .wrapping_mul(6_364_136_223_846_793_005)
        .wrapping_add(1_442_695_040_888_963_407);
    (*state >> 11) as f64 / (1_u64 << 53) as f64
}

/// A single IB point imposes exactly its own velocity on every node of its
/// kernel support.  The old code multiplied the kernel weight in twice and
/// divided by one factor only, i.e. it imposed `omega * u`.
#[test]
fn test_transfer_velocity_single_point_imposes_its_velocity_on_every_stencil_node() {
    let dims = [12, 12, 12];
    let velocity = [0.7, -0.2, 1.3];
    let coupling = moving_points_coupling(&[[5.3, 4.8, 6.1]], &[velocity], 1.0);
    let bcs = coupling
        .transfer_velocity_to_lbm(dims)
        .expect("initialised");

    // |x/h - i| < 3/2 along every axis: nodes 4..=6 x 4..=6 x 5..=7, each
    // exactly once and in ascending row-major order
    let mut expected_nodes = Vec::new();
    for i in 4..=6 {
        for j in 4..=6 {
            for k in 5..=7 {
                expected_nodes.push([i, j, k]);
            }
        }
    }
    let nodes: Vec<[usize; 3]> = bcs.iter().map(|(node, _)| *node).collect();
    assert_eq!(nodes, expected_nodes);
    for (node, bc) in &bcs {
        assert_vec_close(*bc, velocity, 1e-12, &format!("node {node:?}"));
    }
}

/// For a point exactly on a grid node the old code imposed `(2/3)^3 u` on
/// that node and only `u / 216` on the corners of the stencil.
#[test]
fn test_transfer_velocity_point_on_a_node_is_not_attenuated() {
    let dims = [12, 12, 12];
    let velocity = [-1.5, 0.25, 0.5];
    let coupling = moving_points_coupling(&[[6.0, 6.0, 6.0]], &[velocity], 1.0);
    let bcs = coupling
        .transfer_velocity_to_lbm(dims)
        .expect("initialised");
    assert_eq!(bcs.len(), 27);
    for node in [[6, 6, 6], [5, 5, 5], [7, 7, 7], [5, 7, 6]] {
        let bc = bc_at(&bcs, node).expect("node in the support");
        assert_vec_close(bc, velocity, 1e-12, &format!("node {node:?}"));
    }
}

#[test]
fn test_transfer_velocity_blends_overlapping_points_by_kernel_weight() {
    let dims = [12, 12, 12];
    let points = [[5.0, 5.0, 5.0], [6.0, 5.0, 5.0]];
    let velocities = [[1.0, 0.0, 0.0], [0.0, 2.0, -1.0]];
    let coupling = moving_points_coupling(&points, &velocities, 1.0);
    let bcs = coupling
        .transfer_velocity_to_lbm(dims)
        .expect("initialised");
    // x-support of the first point: 4..=6, of the second: 5..=7
    assert_eq!(bcs.len(), 4 * 3 * 3);

    let mut shared = 0;
    for (node, bc) in &bcs {
        let w0 = node_weight(points[0], 1.0, *node);
        let w1 = node_weight(points[1], 1.0, *node);
        let mut expected = [0.0; 3];
        for ((e, a), b) in expected
            .iter_mut()
            .zip(velocities[0].iter())
            .zip(velocities[1].iter())
        {
            *e = (w0 * a + w1 * b) / (w0 + w1);
        }
        assert_vec_close(*bc, expected, 1e-12, &format!("node {node:?}"));
        if w0 > 0.0 && w1 > 0.0 {
            shared += 1;
        }
    }
    assert_eq!(shared, 2 * 3 * 3);

    // On the axis through both points the weights are phi(0) : phi(1) = 4 : 1
    // under the first point and 1 : 4 under the second one.
    let under_first = bc_at(&bcs, [5, 5, 5]).expect("node in the support");
    assert_vec_close(
        under_first,
        [0.8, 0.4, -0.2],
        1e-12,
        "under the first point",
    );
    let under_second = bc_at(&bcs, [6, 5, 5]).expect("node in the support");
    assert_vec_close(
        under_second,
        [0.2, 1.6, -0.8],
        1e-12,
        "under the second point",
    );
    // Nodes that only one point reaches keep that point's velocity.
    let first_only = bc_at(&bcs, [4, 5, 5]).expect("node in the support");
    assert_vec_close(first_only, velocities[0], 1e-12, "first point only");
    let second_only = bc_at(&bcs, [7, 5, 5]).expect("node in the support");
    assert_vec_close(second_only, velocities[1], 1e-12, "second point only");
}

#[test]
fn test_transfer_velocity_is_independent_of_grid_spacing() {
    let dims = [10, 10, 10];
    let points = [[4.3, 5.1, 3.8], [5.0, 5.6, 4.2]];
    let velocities = [[0.3, -0.1, 0.5], [-0.4, 0.2, 0.1]];
    let reference = moving_points_coupling(&points, &velocities, 1.0)
        .transfer_velocity_to_lbm(dims)
        .expect("initialised");
    assert!(!reference.is_empty());
    for h in [0.25, 3.0] {
        let scaled: Vec<[f64; 3]> = points
            .iter()
            .map(|p| [p[0] * h, p[1] * h, p[2] * h])
            .collect();
        let bcs = moving_points_coupling(&scaled, &velocities, h)
            .transfer_velocity_to_lbm(dims)
            .expect("initialised");
        assert_eq!(bcs.len(), reference.len(), "h = {h}");
        for ((node_ref, bc_ref), (node, bc)) in reference.iter().zip(bcs.iter()) {
            assert_eq!(node_ref, node, "h = {h}");
            assert_vec_close(*bc, *bc_ref, 1e-10, &format!("h = {h}, node {node:?}"));
        }
    }
}

#[test]
fn test_transfer_velocity_clips_the_stencil_at_the_grid_boundary() {
    let dims = [6, 6, 6];
    let velocity = [0.5, 0.25, -0.75];

    // Low faces: x = 0.2 -> nodes 0..=1 (node -1 does not exist), y = 0 ->
    // nodes 0..=1, z = 3 -> nodes 2..=4.
    let low = moving_points_coupling(&[[0.2, 0.0, 3.0]], &[velocity], 1.0);
    let bcs = low.transfer_velocity_to_lbm(dims).expect("initialised");
    assert_eq!(bcs.len(), 2 * 2 * 3);
    for (node, bc) in &bcs {
        assert!(
            node[0] <= 1 && node[1] <= 1 && (2..=4).contains(&node[2]),
            "{node:?}"
        );
        assert_vec_close(*bc, velocity, 1e-12, &format!("node {node:?}"));
    }

    // High faces: x = 5.9 -> node 5 only (6 and 7 do not exist), y = 5 ->
    // nodes 4..=5, z = 2.5 -> nodes 2..=3 (the weights at +-1.5 vanish).
    let high = moving_points_coupling(&[[5.9, 5.0, 2.5]], &[velocity], 1.0);
    let bcs = high.transfer_velocity_to_lbm(dims).expect("initialised");
    let nodes: Vec<[usize; 3]> = bcs.iter().map(|(node, _)| *node).collect();
    assert_eq!(nodes, [[5, 4, 2], [5, 4, 3], [5, 5, 2], [5, 5, 3]]);
    for (node, bc) in &bcs {
        assert_vec_close(*bc, velocity, 1e-12, &format!("node {node:?}"));
    }
}

#[test]
fn test_transfer_velocity_ignores_points_outside_the_grid_and_empty_grids() {
    let velocity = [1.0, 0.0, 0.0];
    let outside = moving_points_coupling(
        &[
            [20.0, 3.0, 3.0],
            [-4.0, 3.0, 3.0],
            [3.0, 3.0, -1.6],
            [f64::NAN, 3.0, 3.0],
        ],
        &[velocity; 4],
        1.0,
    );
    let bcs = outside
        .transfer_velocity_to_lbm([6, 6, 6])
        .expect("initialised");
    assert!(bcs.is_empty());

    let inside = moving_points_coupling(&[[3.0, 3.0, 3.0]], &[velocity], 1.0);
    for dims in [[0, 6, 6], [6, 0, 6], [6, 6, 0]] {
        let bcs = inside.transfer_velocity_to_lbm(dims).expect("initialised");
        assert!(bcs.is_empty(), "dims {dims:?}");
    }
}

#[test]
fn test_transfer_velocity_rejects_invalid_grid_spacing() {
    let mut coupling = moving_points_coupling(&[[3.0, 3.0, 3.0]], &[[1.0, 0.0, 0.0]], 1.0);
    for bad in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        coupling.lbm_grid_spacing = bad;
        assert!(
            matches!(
                coupling.transfer_velocity_to_lbm([8, 8, 8]),
                Err(CouplingError::InvalidParameter(_))
            ),
            "spacing {bad}"
        );
    }
}

/// Scattered points: the list is strictly ascending in the node index, a
/// repeated call gives bit-identical output, and every entry equals the
/// brute-force kernel-weighted average over all points.
#[test]
fn test_transfer_velocity_scattered_points_match_brute_force_average() {
    let dims = [14, 14, 14];
    let mut state = 0x2545_f491_4f6c_dd1d_u64;
    let mut points = Vec::new();
    let mut velocities = Vec::new();
    for _ in 0..60 {
        points.push([
            3.0 + 8.0 * lcg_unit(&mut state),
            3.0 + 8.0 * lcg_unit(&mut state),
            3.0 + 8.0 * lcg_unit(&mut state),
        ]);
        velocities.push([
            lcg_unit(&mut state) - 0.5,
            lcg_unit(&mut state) - 0.5,
            lcg_unit(&mut state) - 0.5,
        ]);
    }
    let coupling = moving_points_coupling(&points, &velocities, 1.0);
    let first = coupling
        .transfer_velocity_to_lbm(dims)
        .expect("initialised");
    let second = coupling
        .transfer_velocity_to_lbm(dims)
        .expect("initialised");
    assert_eq!(first, second);
    assert!(
        first.windows(2).all(|pair| pair[0].0 < pair[1].0),
        "nodes must be unique and ascending"
    );

    let mut checked = 0;
    for i in 0..dims[0] {
        for j in 0..dims[1] {
            for k in 0..dims[2] {
                let node = [i, j, k];
                let weights: Vec<f64> = points.iter().map(|p| node_weight(*p, 1.0, node)).collect();
                let total: f64 = weights.iter().sum();
                let imposed = bc_at(&first, node);
                if total <= 1e-12 {
                    // nodes outside every support are not reported
                    assert!(total > 0.0 || imposed.is_none(), "node {node:?}");
                    continue;
                }
                let mut expected = [0.0; 3];
                for (w, v) in weights.iter().zip(velocities.iter()) {
                    for (e, v_d) in expected.iter_mut().zip(v.iter()) {
                        *e += w * v_d / total;
                    }
                }
                let bc = imposed.expect("node in the support of some point");
                assert_vec_close(bc, expected, 1e-12, &format!("node {node:?}"));
                checked += 1;
            }
        }
    }
    assert!(checked > 500, "only {checked} nodes checked");
}

#[test]
fn test_kernel_stencil_weights_form_a_partition_of_unity() {
    let dims = [16, 16, 16];
    let h = 0.5;
    let mut state = 7_u64;
    for _ in 0..200 {
        let pos = [
            h * (4.0 + 8.0 * lcg_unit(&mut state)),
            h * (4.0 + 8.0 * lcg_unit(&mut state)),
            h * (4.0 + 8.0 * lcg_unit(&mut state)),
        ];
        let stencil = kernel_stencil(pos, h, dims);
        let total: f64 = stencil.iter().map(|(_, weight)| weight).sum();
        assert!((total - 1.0).abs() < 1e-12, "sum {total} at {pos:?}");
        assert!((8..=27).contains(&stencil.len()), "{} nodes", stencil.len());
        assert!(stencil.iter().all(|(_, weight)| *weight > 0.0));
    }
}

#[test]
fn test_kernel_axis_weights_edge_cases() {
    // no nodes: empty axis, non-finite or absurdly distant coordinates
    assert!(kernel_axis_weights(3.0, 0).is_empty());
    assert!(kernel_axis_weights(f64::INFINITY, 8).is_empty());
    assert!(kernel_axis_weights(f64::NAN, 8).is_empty());
    assert!(kernel_axis_weights(1.0e300, 8).is_empty());
    assert!(kernel_axis_weights(-1.0e300, 8).is_empty());

    // the kernel vanishes at a distance of 3/2 ...
    assert!(kernel_axis_weights(-1.5, 8).is_empty());
    assert!(kernel_axis_weights(8.5, 8).is_empty());
    // ... and is small but positive just inside the support
    let near_edge = kernel_axis_weights(-1.4, 8);
    assert_eq!(near_edge.len(), 1);
    assert_eq!(near_edge[0].0, 0);
    assert!(near_edge[0].1 > 0.0 && near_edge[0].1 < 0.02);

    // a one-node axis keeps its own node: phi(0) = 2/3
    let single = kernel_axis_weights(0.0, 1);
    assert_eq!(single.len(), 1);
    assert_eq!(single[0].0, 0);
    assert!((single[0].1 - 2.0 / 3.0).abs() < 1e-15);

    // nodes are ascending and the weights positive
    let weights = kernel_axis_weights(3.3, 8);
    assert_eq!(
        weights.iter().map(|(node, _)| *node).collect::<Vec<_>>(),
        [2, 3, 4]
    );
    assert!(weights.iter().all(|(_, weight)| *weight > 0.0));
}

#[test]
fn test_coupling_step_returns_the_structure_velocity_as_boundary_condition() {
    let velocity = [1.0, -2.0, 0.5];
    let mut coupling = moving_points_coupling(&[[3.0, 3.0, 3.0]], &[velocity], 1.0);
    let dims = [8, 8, 8];
    let total = 512;
    let density = vec![1.0; total];
    let fluid_velocity = vec![[0.0; 3]; total];
    let bcs = coupling
        .coupling_step(&density, &fluid_velocity, 1.0, dims)
        .expect("fields match the grid");
    assert_eq!(bcs.len(), 27);
    for (node, bc) in &bcs {
        assert_vec_close(*bc, velocity, 1e-12, &format!("node {node:?}"));
    }
}
