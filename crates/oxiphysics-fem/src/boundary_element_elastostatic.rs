// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Influence matrices of the elastostatic boundary element method.
//!
//! The collocation equations of 3-D linear elasticity with constant elements
//! are `H u = G t`.  `G` collects the Kelvin displacement kernel `U_ij`, `H` the
//! Kelvin traction kernel `T_ij`; the singular self terms are handled
//! analytically (`G`) and by the rigid-body condition (`H`).

use std::f64::consts::PI;

use super::{BoundaryMesh, DenseMatrix, kelvin_displacement_3d, kelvin_traction_3d, normalize3};

/// Self-influence block of the elastostatic single-layer matrix `G` of a
/// constant element collocated at its centroid.
///
/// The element is replaced by the disk of equal area (radius `R = sqrt(A/pi)`)
/// and the Kelvin displacement kernel `U_ij` is integrated over it in closed
/// form.  In the local frame (two tangents `t1`, `t2` and the normal `n`; the
/// vector `x - y` lies in the element plane, so `r_n = 0`):
///
/// * tangential: `G_tt = R (7 - 8 nu) / (16 G (1 - nu))`,
/// * normal: `G_nn = R (3 - 4 nu) / (8 G (1 - nu))`,
/// * no tangential-tangential or tangential-normal coupling.
///
/// Rotated to global axes the block is `G_tt (I - n n^T) + G_nn n n^T`, so no
/// tangent basis has to be constructed.
pub(super) fn elastostatic_self_g_block(
    area: f64,
    normal: &[f64; 3],
    shear_modulus: f64,
    nu: f64,
) -> [[f64; 3]; 3] {
    let radius = (area / PI).sqrt();
    let denominator = 16.0 * shear_modulus * (1.0 - nu);
    let tangential = radius * (7.0 - 8.0 * nu) / denominator;
    let normal_part = 2.0 * radius * (3.0 - 4.0 * nu) / denominator;
    let mut n = *normal;
    normalize3(&mut n);
    let mut block = [[0.0; 3]; 3];
    for (i, row) in block.iter_mut().enumerate() {
        for (j, entry) in row.iter_mut().enumerate() {
            let identity = if i == j { tangential } else { 0.0 };
            *entry = identity + (normal_part - tangential) * n[i] * n[j];
        }
    }
    block
}

/// Assemble H and G matrices for 3-D elastostatics BEM (constant elements).
///
/// Each element carries 3 DOFs (u_x, u_y, u_z), so H and G are `(3n x 3n)` and
/// the collocation equations read `H u = G t` (`u`, `t` = displacement and
/// traction at the element centroids).
///
/// * Off-diagonal blocks use the Kelvin kernels `U_ij` and `T_ij` at the element
///   centroid times the element area; `T_ij` is evaluated with the field point
///   (the element that carries the normal) first.
/// * The self block of `G` is the closed-form disk integral of `U_ij`
///   (see `elastostatic_self_g_block`).
/// * The diagonal blocks of `H` follow from the rigid-body condition: a
///   translation `u = e_k` of every element produces no traction, so for every
///   pair of directions `H[3i+a, 3i+b] = -sum_{j != i} H[3i+a, 3j+b]`.  This
///   fills the complete 3x3 diagonal block (free term plus the principal-value
///   integral), not just its diagonal.
pub fn assemble_elastostatic_bem(
    mesh: &BoundaryMesh,
    shear_modulus: f64,
    nu: f64,
) -> (DenseMatrix, DenseMatrix) {
    let n = mesh.num_elements();
    let dim = 3 * n;
    let mut h_mat = DenseMatrix::zeros(dim, dim);
    let mut g_mat = DenseMatrix::zeros(dim, dim);

    for ie in 0..n {
        let xi = mesh.element_centroid(ie);
        for je in 0..n {
            let aj = mesh.elements[je].area;
            if ie == je {
                let block =
                    elastostatic_self_g_block(aj, &mesh.elements[je].normal, shear_modulus, nu);
                for (di, block_row) in block.iter().enumerate() {
                    for (dj, &value) in block_row.iter().enumerate() {
                        g_mat.set(ie * 3 + di, je * 3 + dj, value);
                    }
                }
                continue;
            }
            let xj = mesh.element_centroid(je);
            let nj = mesh.elements[je].normal;
            for di in 0..3_usize {
                for dj in 0..3_usize {
                    let u_ij = kelvin_displacement_3d(&xi, &xj, di, dj, shear_modulus, nu);
                    // The normal `nj` belongs to the integration (field) point
                    // `xj`, which therefore comes first.
                    let t_ij = kelvin_traction_3d(&xj, &xi, &nj, di, dj, nu);
                    g_mat.set(ie * 3 + di, je * 3 + dj, u_ij * aj);
                    h_mat.set(ie * 3 + di, je * 3 + dj, t_ij * aj);
                }
            }
        }
    }

    // Rigid body motion for elastostatics: H * {e_k} = 0 for each direction k, so
    // the whole 3x3 diagonal block of H is minus the sum of the off-diagonal ones.
    for i in 0..n {
        for di in 0..3_usize {
            for dk in 0..3_usize {
                let mut off_sum = 0.0;
                for j in 0..n {
                    if j != i {
                        off_sum += h_mat.get(i * 3 + di, j * 3 + dk);
                    }
                }
                h_mat.set(i * 3 + di, i * 3 + dk, -off_sum);
            }
        }
    }

    (h_mat, g_mat)
}
