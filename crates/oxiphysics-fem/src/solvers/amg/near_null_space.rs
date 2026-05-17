// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Near-null-space vectors (rigid body modes) for smoothed-aggregation AMG.

#![allow(dead_code)]

/// Compute rigid body modes for a 3D elasticity aggregate.
///
/// For scalar problems (`n_dof_per_node = 1`), returns a single vector of ones.
///
/// For 3D elasticity (`n_dof_per_node = 3`), returns up to 6 vectors:
/// - 3 translations (constant displacement fields)
/// - 3 rotations about the centroid of the nodes
///
/// The vectors are orthonormalized via modified Gram-Schmidt before return.
/// Linearly dependent vectors (norm < 1e-14 after projection) are removed.
///
/// # Arguments
/// * `node_coords` — 3D coordinates for each node in the aggregate
/// * `n_dof_per_node` — 1 for scalar problems, 3 for 3D elasticity
pub fn rigid_body_modes_3d(node_coords: &[[f64; 3]], n_dof_per_node: usize) -> Vec<Vec<f64>> {
    let k = node_coords.len();
    if k == 0 {
        return Vec::new();
    }

    let dof_count = k * n_dof_per_node;

    if n_dof_per_node == 1 {
        // Scalar: constant-one vector
        return vec![vec![1.0; dof_count]];
    }

    // Compute centroid
    let mut cx = 0.0f64;
    let mut cy = 0.0f64;
    let mut cz = 0.0f64;
    for coord in node_coords {
        cx += coord[0];
        cy += coord[1];
        cz += coord[2];
    }
    let inv_k = 1.0 / k as f64;
    cx *= inv_k;
    cy *= inv_k;
    cz *= inv_k;

    // Build 6 rigid body mode vectors
    let mut vecs: Vec<Vec<f64>> = Vec::with_capacity(6);

    // 3 Translations
    // e1: (1,0,0,...,1,0,0)
    let mut e1 = vec![0.0f64; dof_count];
    let mut e2 = vec![0.0f64; dof_count];
    let mut e3 = vec![0.0f64; dof_count];
    for j in 0..k {
        e1[3 * j] = 1.0;
        e2[3 * j + 1] = 1.0;
        e3[3 * j + 2] = 1.0;
    }
    vecs.push(e1);
    vecs.push(e2);
    vecs.push(e3);

    // 3 Rotations about centroid
    // Rotation about x-axis: ω_x for node j = (0, -(z_j-cz), y_j-cy)
    let mut rx = vec![0.0f64; dof_count];
    // Rotation about y-axis: ω_y for node j = (z_j-cz, 0, -(x_j-cx))
    let mut ry = vec![0.0f64; dof_count];
    // Rotation about z-axis: ω_z for node j = (-(y_j-cy), x_j-cx, 0)
    let mut rz = vec![0.0f64; dof_count];

    for j in 0..k {
        let xr = node_coords[j][0] - cx;
        let yr = node_coords[j][1] - cy;
        let zr = node_coords[j][2] - cz;

        // Rotation about x: (0, -z, y)
        rx[3 * j] = 0.0;
        rx[3 * j + 1] = -zr;
        rx[3 * j + 2] = yr;

        // Rotation about y: (z, 0, -x)
        ry[3 * j] = zr;
        ry[3 * j + 1] = 0.0;
        ry[3 * j + 2] = -xr;

        // Rotation about z: (-y, x, 0)
        rz[3 * j] = -yr;
        rz[3 * j + 1] = xr;
        rz[3 * j + 2] = 0.0;
    }
    vecs.push(rx);
    vecs.push(ry);
    vecs.push(rz);

    // Orthonormalize via modified Gram-Schmidt
    gram_schmidt(&mut vecs);
    vecs
}

/// Modified Gram-Schmidt orthonormalization in-place.
///
/// Normalizes each vector after projecting out previously processed vectors.
/// Removes linearly dependent vectors (post-projection norm < 1e-14).
pub fn gram_schmidt(vecs: &mut Vec<Vec<f64>>) {
    let mut i = 0;
    while i < vecs.len() {
        // Compute norm of vecs[i]
        let norm_sq: f64 = vecs[i].iter().map(|v| v * v).sum();
        let norm = norm_sq.sqrt();

        if norm < 1e-14 {
            // Remove linearly dependent vector
            vecs.remove(i);
            continue;
        }

        // Normalize
        let inv_norm = 1.0 / norm;
        for v in vecs[i].iter_mut() {
            *v *= inv_norm;
        }

        // Project out from all subsequent vectors
        // Clone the current vector so we can mutate the rest of the slice
        let qi: Vec<f64> = vecs[i].clone();
        let qi_len = qi.len();
        for vj in vecs[(i + 1)..].iter_mut() {
            let dot: f64 = vj.iter().zip(qi.iter()).map(|(a, b)| a * b).sum();
            for k in 0..qi_len {
                vj[k] -= dot * qi[k];
            }
        }

        i += 1;
    }
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    fn dot(a: &[f64], b: &[f64]) -> f64 {
        a.iter().zip(b.iter()).map(|(x, y)| x * y).sum()
    }

    #[test]
    fn test_rbm_6_vectors_3d() {
        // 4 nodes at different 3D positions → should produce 6 orthonormal RBMs
        let coords: Vec<[f64; 3]> = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];

        let vecs = rigid_body_modes_3d(&coords, 3);

        assert_eq!(vecs.len(), 6, "Expected 6 RBMs, got {}", vecs.len());

        // Each pair must be orthogonal
        for i in 0..6 {
            // Self-dot must be ≈ 1.0 (unit vectors)
            let self_dot = dot(&vecs[i], &vecs[i]);
            assert!(
                (self_dot - 1.0).abs() < 1e-12,
                "RBM {i} is not unit: ||v||² = {self_dot}"
            );

            for j in (i + 1)..6 {
                let cross_dot = dot(&vecs[i], &vecs[j]);
                assert!(
                    cross_dot.abs() < 1e-12,
                    "RBMs {i} and {j} are not orthogonal: dot = {cross_dot}"
                );
            }
        }
    }

    #[test]
    fn test_rbm_scalar() {
        let coords: Vec<[f64; 3]> = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
        let vecs = rigid_body_modes_3d(&coords, 1);
        assert_eq!(vecs.len(), 1, "Scalar should produce 1 vector");
        assert_eq!(
            vecs[0].len(),
            2,
            "Scalar vector should have n_nodes entries"
        );
    }

    #[test]
    fn test_gram_schmidt_orthonormal() {
        // 3 random-ish vectors of length 4
        let mut vecs: Vec<Vec<f64>> = vec![
            vec![1.0, 2.0, 3.0, 4.0],
            vec![2.0, 3.0, 4.0, 5.0],
            vec![3.0, 1.0, 4.0, 1.0],
        ];

        gram_schmidt(&mut vecs);

        // After GS, must be orthonormal
        let n_vecs = vecs.len();
        for i in 0..n_vecs {
            let norm_sq = dot(&vecs[i], &vecs[i]);
            assert!(
                (norm_sq - 1.0).abs() < 1e-12,
                "Vector {i} not unit after GS: ||v||² = {norm_sq}"
            );
            for j in (i + 1)..n_vecs {
                let d = dot(&vecs[i], &vecs[j]);
                assert!(
                    d.abs() < 1e-12,
                    "Vectors {i} and {j} not orthogonal after GS: dot = {d}"
                );
            }
        }
    }

    #[test]
    fn test_gram_schmidt_removes_dependent() {
        // Two identical vectors — one should be removed
        let mut vecs: Vec<Vec<f64>> = vec![
            vec![1.0, 0.0, 0.0],
            vec![1.0, 0.0, 0.0], // duplicate
            vec![0.0, 1.0, 0.0],
        ];
        gram_schmidt(&mut vecs);
        // Should produce 2 orthonormal vectors
        assert_eq!(
            vecs.len(),
            2,
            "Expected 2 independent vectors, got {}",
            vecs.len()
        );
    }
}
