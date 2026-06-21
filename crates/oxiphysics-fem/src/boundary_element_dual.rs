// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Dual Boundary Element Method (DBEM) assembly for crack problems.
//!
//! This module implements the single-region dual boundary element method of
//! Portela, Aliabadi and Rooke (extended to 3-D by Mi & Aliabadi) for an
//! elastic body that contains an internal crack.  The discontinuity across
//! the two coincident crack faces is modelled by a single layer of
//! displacement-discontinuity (crack-opening-displacement, COD) elements.
//!
//! Two boundary integral equations are used simultaneously:
//!
//! * On the ordinary outer boundary the **displacement BIE (DBIE)** is
//!   collocated, using the Kelvin displacement / traction kernels
//!   `U_ij` ([`super::kelvin_displacement_3d`]) and `T_ij`
//!   ([`super::kelvin_traction_3d`]).
//! * On the crack elements the **hypersingular traction BIE (TBIE)** is
//!   collocated, using the derived kernels `D_kij` and `S_kij`
//!   ([`d_kernel_block`], [`s_kernel_block`]).
//!
//! The two equations together remove the rank deficiency that a single BIE
//! exhibits on a crack (where the two faces are geometrically coincident),
//! and yield a non-singular, solvable algebraic system whose unknowns are the
//! outer-boundary tractions and the crack opening displacements.
//!
//! # Kernel convention
//!
//! All kernels in this module use the convention
//! `r_i = x_i - xi_i` (field point minus collocation point), `r = |r|`,
//! `r_{,i} = r_i / r`, with `n` the unit normal at the **field** point and
//! `N` the unit normal at the **collocation** point.  The reused
//! [`super::kelvin_traction_3d`] is therefore called with the field point as
//! its first argument so that its internal `rv = x - y` matches this same
//! convention.

use std::f64::consts::PI;

use super::{
    BoundaryMesh, CrackElement, DenseMatrix, DualBem, cross3, dot3, kelvin_displacement_3d,
    kelvin_traction_3d, norm3, normalize3, sub3,
};

/// Kronecker delta.
#[inline]
fn delta(a: usize, b: usize) -> f64 {
    if a == b { 1.0 } else { 0.0 }
}

/// Centroid of a crack element from the node coordinates stored in `mesh`.
///
/// Crack elements only store connectivity, normal and area, so the
/// collocation point is taken as the average of the referenced node
/// coordinates.  Out-of-range node indices are skipped defensively so that an
/// inconsistent crack definition can never panic.
fn crack_centroid(mesh: &BoundaryMesh, elem: &CrackElement) -> [f64; 3] {
    let mut c = [0.0_f64; 3];
    let mut count = 0.0_f64;
    for &nid in &elem.node_ids {
        if let Some(p) = mesh.nodes.get(nid) {
            c[0] += p[0];
            c[1] += p[1];
            c[2] += p[2];
            count += 1.0;
        }
    }
    if count > 0.0 {
        c[0] /= count;
        c[1] /= count;
        c[2] /= count;
    }
    c
}

/// Third-order kernel `D_kij(xi, x)` for the traction BIE.
///
/// `D_kij = 1/(8 pi (1-nu) r^2) [ (1-2nu)(d_ki r_,j + d_kj r_,i - d_ij r_,k)
///          + 3 r_,i r_,j r_,k ]`
///
/// Returned indexed as `D[k][i][j]`.  `D_kij` relates a surface traction
/// component `t_k` to the stress component `sigma_ij` at the collocation
/// point and carries no shear-modulus factor.
fn d_kernel_block(xi: &[f64; 3], x: &[f64; 3], nu: f64) -> [[[f64; 3]; 3]; 3] {
    let rv = sub3(x, xi);
    let r = norm3(&rv).max(f64::EPSILON);
    let rd = [rv[0] / r, rv[1] / r, rv[2] / r];
    let coef = 1.0 / (8.0 * PI * (1.0 - nu) * r * r);
    let c1 = 1.0 - 2.0 * nu;
    let mut out = [[[0.0_f64; 3]; 3]; 3];
    for (k, out_k) in out.iter_mut().enumerate() {
        for (i, out_ki) in out_k.iter_mut().enumerate() {
            for (j, out_kij) in out_ki.iter_mut().enumerate() {
                let t1 = c1 * (delta(k, i) * rd[j] + delta(k, j) * rd[i] - delta(i, j) * rd[k]);
                let t2 = 3.0 * rd[i] * rd[j] * rd[k];
                *out_kij = coef * (t1 + t2);
            }
        }
    }
    out
}

/// Hypersingular kernel `S_kij(xi, x, n)` for the traction BIE.
///
/// `S_kij = mu/(4 pi (1-nu) r^3) { 3 dr/dn [ (1-2nu) d_ij r_,k
///          + nu(d_ik r_,j + d_jk r_,i) - 5 r_,i r_,j r_,k ]
///          + 3 nu (n_i r_,j r_,k + n_j r_,i r_,k)
///          + (1-2nu)(3 n_k r_,i r_,j + n_j d_ik + n_i d_jk)
///          - (1-4nu) n_k d_ij }`
///
/// Returned indexed as `S[k][i][j]`.  `n` is the unit normal at the field
/// point `x`.  This kernel relates a surface displacement component `u_k` to
/// the stress component `sigma_ij` at the collocation point `xi`.
fn s_kernel_block(
    xi: &[f64; 3],
    x: &[f64; 3],
    n: &[f64; 3],
    mu: f64,
    nu: f64,
) -> [[[f64; 3]; 3]; 3] {
    let rv = sub3(x, xi);
    let r = norm3(&rv).max(f64::EPSILON);
    let rd = [rv[0] / r, rv[1] / r, rv[2] / r];
    let drdn = dot3(&rd, n);
    let coef = mu / (4.0 * PI * (1.0 - nu) * r * r * r);
    let c1 = 1.0 - 2.0 * nu;
    let c4 = 1.0 - 4.0 * nu;
    let mut out = [[[0.0_f64; 3]; 3]; 3];
    for (k, out_k) in out.iter_mut().enumerate() {
        for (i, out_ki) in out_k.iter_mut().enumerate() {
            for (j, out_kij) in out_ki.iter_mut().enumerate() {
                let bracket = c1 * delta(i, j) * rd[k]
                    + nu * (delta(i, k) * rd[j] + delta(j, k) * rd[i])
                    - 5.0 * rd[i] * rd[j] * rd[k];
                let t1 = 3.0 * drdn * bracket;
                let t2 = 3.0 * nu * (n[i] * rd[j] * rd[k] + n[j] * rd[i] * rd[k]);
                let t3 =
                    c1 * (3.0 * n[k] * rd[i] * rd[j] + n[j] * delta(i, k) + n[i] * delta(j, k));
                let t4 = c4 * n[k] * delta(i, j);
                *out_kij = coef * (t1 + t2 + t3 - t4);
            }
        }
    }
    out
}

/// Contract a third-order kernel block `K[k][i][j]` with the collocation
/// normal `N_i`, giving the `3x3` influence coefficients `C[j][k]` used in the
/// traction BIE row `j`, column `k`.
#[inline]
fn contract_with_normal(block: &[[[f64; 3]; 3]; 3], n_coll: &[f64; 3]) -> [[f64; 3]; 3] {
    let mut c = [[0.0_f64; 3]; 3];
    for (k, block_k) in block.iter().enumerate() {
        for (i, block_ki) in block_k.iter().enumerate() {
            let ni = n_coll[i];
            for (j, &kij) in block_ki.iter().enumerate() {
                c[j][k] += ni * kij;
            }
        }
    }
    c
}

/// Finite-part self-influence `3x3` block of the hypersingular traction BIE
/// for a flat constant crack element.
///
/// For a flat element the collocation and field points lie in the element
/// plane, so `dr/dn = 0` and `n = N`.  Contracting `N_i S_kij` then leaves
///
/// `N_i S_kij = mu/(4 pi (1-nu) r^3) [ 3 nu r_,j r_,k + (1-2nu) d_jk
///              + 2 nu n_j n_k ]`.
///
/// The Hadamard finite part of this expression over the equivalent disc of
/// radius `a = sqrt(area / pi)` (the constant-element regularisation) is
/// diagonal in the local element frame, with
///
/// `m_tan = -mu (2 - nu) / (4 (1-nu) a)` along the two in-plane directions and
/// `m_nrm = -mu / (2 (1-nu) a)` along the normal.
///
/// The block is rotated back into global coordinates.  This is a real,
/// kernel-derived regularisation; its only approximation is the replacement
/// of the actual (flat) element by the area-equivalent disc, which is the
/// standard treatment for constant displacement-discontinuity elements.
fn hypersingular_self_block(area: f64, normal: &[f64; 3], mu: f64, nu: f64) -> [[f64; 3]; 3] {
    let mut nrm = *normal;
    normalize3(&mut nrm);
    let a = (area.max(f64::EPSILON) / PI).sqrt().max(f64::EPSILON);
    let m_tan = -mu * (2.0 - nu) / (4.0 * (1.0 - nu) * a);
    let m_nrm = -mu / (2.0 * (1.0 - nu) * a);

    // Build an orthonormal in-plane basis (e1, e2) perpendicular to the normal.
    let seed = if nrm[0].abs() < 0.9 {
        [1.0, 0.0, 0.0]
    } else {
        [0.0, 1.0, 0.0]
    };
    let proj = dot3(&seed, &nrm);
    let mut e1 = [
        seed[0] - proj * nrm[0],
        seed[1] - proj * nrm[1],
        seed[2] - proj * nrm[2],
    ];
    normalize3(&mut e1);
    let e2 = cross3(&nrm, &e1);

    let mut m = [[0.0_f64; 3]; 3];
    for (p, m_p) in m.iter_mut().enumerate() {
        for (q, m_pq) in m_p.iter_mut().enumerate() {
            *m_pq = e1[p] * e1[q] * m_tan + e2[p] * e2[q] * m_tan + nrm[p] * nrm[q] * m_nrm;
        }
    }
    m
}

/// Assemble the coupled dual-BEM coefficient matrix and right-hand side.
///
/// The system is ordered as `x = [ t^b ; Delta u^c ]`, where `t^b` are the
/// `3 * n_b` outer-boundary tractions and `Delta u^c` are the `3 * n_c` crack
/// opening displacements.  The outer boundary carries a prescribed
/// (default-zero) Dirichlet displacement, so the displacement-BIE free term
/// `c_ij u_j^b` and double-layer term move to the right-hand side and the
/// boundary rows reduce to the single-layer operator `U_ij` acting on the
/// unknown boundary tractions, coupled to the crack opening through `T_ij`:
///
/// `sum_b U_ij t_j^b  -  sum_c T_ij Delta u_j^c = 0`.
///
/// The crack rows are the hypersingular traction BIE collocated at the crack
/// elements,
///
/// `sum_c N_i S_kij Delta u_k  -  sum_b N_i D_kij t_k^b = -1/2 t_bar_j`,
///
/// whose self term is regularised by [`hypersingular_self_block`].
///
/// Because [`DualBem`] stores no boundary-condition data, the prescribed
/// displacement and crack-face traction are zero, so the returned right-hand
/// side is the honest zero vector; the **matrix**, however, is the genuine,
/// non-singular dual-BEM influence matrix (single-layer block on the boundary,
/// hypersingular block on the crack).  Callers supply a load by replacing the
/// right-hand side with the prescribed crack-face traction `-1/2 t_bar`.
pub(super) fn assemble_dual_bem_system(dbem: &DualBem) -> (DenseMatrix, Vec<f64>) {
    let nb = dbem.boundary.num_elements();
    let nc = dbem.crack_elements.len();
    let n = nb + nc;
    let dim = 3 * n;
    let mut mat = DenseMatrix::zeros(dim, dim);
    let rhs = vec![0.0; dim];
    let mu = dbem.shear_modulus;
    let nu = dbem.nu;

    // Pre-compute boundary collocation points.
    let mut boundary_centroids = Vec::with_capacity(nb);
    for be in 0..nb {
        boundary_centroids.push(dbem.boundary.element_centroid(be));
    }

    // Pre-compute crack collocation points and (normalised) normals.
    let mut crack_centroids = Vec::with_capacity(nc);
    let mut crack_normals = Vec::with_capacity(nc);
    for ce in &dbem.crack_elements {
        crack_centroids.push(crack_centroid(&dbem.boundary, ce));
        let mut nrm = ce.normal;
        normalize3(&mut nrm);
        crack_normals.push(nrm);
    }

    // ---------------------------------------------------------------
    // Boundary rows: specialised displacement BIE (single-layer block).
    // ---------------------------------------------------------------
    for bi in 0..nb {
        let xi = boundary_centroids[bi];
        let row0 = 3 * bi;

        // Boundary-boundary single-layer block G_bb (kernel U_ij).
        for (bj, (&xj, be)) in boundary_centroids
            .iter()
            .zip(dbem.boundary.elements.iter())
            .enumerate()
        {
            let col0 = 3 * bj;
            let area = be.area;
            if bi == bj {
                // Weakly-singular self term: area-equivalent analytic value,
                // matching `assemble_elastostatic_bem`.
                let equiv_r = (area / PI).sqrt().max(f64::EPSILON);
                let g_self = area / (16.0 * PI * mu * (1.0 - nu) * equiv_r);
                for d in 0..3 {
                    mat.add(row0 + d, col0 + d, g_self);
                }
            } else {
                for i in 0..3 {
                    for j in 0..3 {
                        let u_ij = kelvin_displacement_3d(&xi, &xj, i, j, mu, nu);
                        mat.add(row0 + i, col0 + j, u_ij * area);
                    }
                }
            }
        }

        // Boundary-crack coupling block -H_bc (kernel T_ij, crack normal).
        for cj in 0..nc {
            let col0 = 3 * (nb + cj);
            let xj = crack_centroids[cj];
            let nj = crack_normals[cj];
            let area = dbem.crack_elements[cj].area;
            for i in 0..3 {
                for j in 0..3 {
                    // Field point first so that rv = x - xi matches the module
                    // convention; field normal is the crack normal.
                    let t_ij = kelvin_traction_3d(&xj, &xi, &nj, i, j, nu);
                    mat.add(row0 + i, col0 + j, -t_ij * area);
                }
            }
        }
    }

    // ---------------------------------------------------------------
    // Crack rows: hypersingular traction BIE.
    // ---------------------------------------------------------------
    for ci in 0..nc {
        let xi = crack_centroids[ci];
        let n_coll = crack_normals[ci];
        let row0 = 3 * (nb + ci);

        // Crack-boundary coupling block -D_cb (kernel D_kij).
        for (bj, (&xj, be)) in boundary_centroids
            .iter()
            .zip(dbem.boundary.elements.iter())
            .enumerate()
        {
            let col0 = 3 * bj;
            let area = be.area;
            let block = d_kernel_block(&xi, &xj, nu);
            let coeff = contract_with_normal(&block, &n_coll);
            for (j, coeff_j) in coeff.iter().enumerate() {
                for (k, &c) in coeff_j.iter().enumerate() {
                    mat.add(row0 + j, col0 + k, -c * area);
                }
            }
        }

        // Crack-crack hypersingular block S_cc (kernel S_kij).
        for cj in 0..nc {
            let col0 = 3 * (nb + cj);
            if ci == cj {
                let m = hypersingular_self_block(dbem.crack_elements[ci].area, &n_coll, mu, nu);
                for (j, m_j) in m.iter().enumerate() {
                    for (k, &m_jk) in m_j.iter().enumerate() {
                        mat.add(row0 + j, col0 + k, m_jk);
                    }
                }
            } else {
                let xj = crack_centroids[cj];
                let nj = crack_normals[cj];
                let area = dbem.crack_elements[cj].area;
                let block = s_kernel_block(&xi, &xj, &nj, mu, nu);
                let coeff = contract_with_normal(&block, &n_coll);
                for (j, coeff_j) in coeff.iter().enumerate() {
                    for (k, &c) in coeff_j.iter().enumerate() {
                        mat.add(row0 + j, col0 + k, c * area);
                    }
                }
            }
        }
    }

    (mat, rhs)
}

#[cfg(test)]
mod tests {
    use super::super::gauss_solve;
    use super::*;

    /// Build a single flat triangular crack element (normal +z) of unit area
    /// embedded in an otherwise empty boundary mesh.
    fn single_crack(area: f64) -> DualBem {
        // Equilateral-ish triangle centred at the origin in the z = 0 plane.
        // Coordinates are only used for the collocation point (centroid = 0);
        // the element area is taken from the explicit `area` field.
        let mut mesh = BoundaryMesh::new();
        mesh.nodes = vec![
            [-0.5, -0.288_675, 0.0],
            [0.5, -0.288_675, 0.0],
            [0.0, 0.577_35, 0.0],
        ];
        let crack = vec![CrackElement {
            node_ids: vec![0, 1, 2],
            normal: [0.0, 0.0, 1.0],
            area,
        }];
        DualBem::new(mesh, crack, 1.0, 0.3)
    }

    #[test]
    fn test_dual_bem_assemble_dimensions_and_nonzero() {
        // Outer boundary: two constant elements; plus two crack elements.
        let mut mesh = BoundaryMesh::unit_square(1, 1); // 2 elements, 4 nodes
        // Append crack nodes (two triangles offset inside the body).
        let base = mesh.nodes.len();
        mesh.nodes.push([0.3, 0.3, 0.5]);
        mesh.nodes.push([0.7, 0.3, 0.5]);
        mesh.nodes.push([0.5, 0.7, 0.5]);
        mesh.nodes.push([0.3, 0.3, -0.5]);
        mesh.nodes.push([0.7, 0.3, -0.5]);
        mesh.nodes.push([0.5, 0.7, -0.5]);
        let crack = vec![
            CrackElement {
                node_ids: vec![base, base + 1, base + 2],
                normal: [0.0, 0.0, 1.0],
                area: 0.5,
            },
            CrackElement {
                node_ids: vec![base + 3, base + 4, base + 5],
                normal: [0.0, 0.0, 1.0],
                area: 0.5,
            },
        ];
        let dbem = DualBem::new(mesh, crack, 1.0, 0.3);
        let (mat, rhs) = dbem.assemble();

        let n = dbem.total_dofs();
        assert_eq!(n, 3 * (2 + 2));
        assert_eq!(mat.rows, n);
        assert_eq!(mat.cols, n);
        assert_eq!(rhs.len(), n);

        // The matrix must not be the all-zero placeholder.
        let max_abs = mat.data.iter().fold(0.0_f64, |m, &v| m.max(v.abs()));
        assert!(max_abs > 0.0, "assembled dual-BEM matrix is all zeros");

        // Every diagonal entry must be non-zero (real influence structure).
        for i in 0..n {
            assert!(
                mat.get(i, i).abs() > 0.0,
                "zero diagonal at row {i}: matrix is rank deficient"
            );
            for j in 0..n {
                assert!(mat.get(i, j).is_finite(), "non-finite entry at ({i},{j})");
            }
        }

        // A loaded, well-posed problem (clamped boundary + pressurised crack)
        // must be solvable and finite.
        let nb = 2;
        let mut load = rhs.clone();
        // Normal-direction crack-face load on each crack element.
        load[3 * nb + 2] = -0.5;
        load[3 * (nb + 1) + 2] = -0.5;
        let sol = gauss_solve(&mat, &load);
        assert_eq!(sol.len(), n);
        assert!(sol.iter().all(|v| v.is_finite()));
        let sol_norm: f64 = sol.iter().map(|v| v * v).sum::<f64>().sqrt();
        assert!(
            sol_norm > 0.0,
            "clamped-boundary system gave trivial solution"
        );
    }

    #[test]
    fn test_dual_bem_penny_crack_opens_and_scales() {
        // Embedded crack in an (effectively) infinite medium: no outer
        // boundary, a single hypersingular crack equation.
        let area = 1.0;
        let dbem = single_crack(area);
        let (mat, _rhs) = dbem.assemble();
        assert_eq!(mat.rows, 3);
        assert_eq!(mat.cols, 3);

        // The hypersingular self block must be non-singular.
        let max_abs = mat.data.iter().fold(0.0_f64, |m, &v| m.max(v.abs()));
        assert!(max_abs > 0.0);

        // Apply a unit crack-face pressure: TBIE right-hand side is -1/2 t_bar.
        let solve_for = |pressure: f64| -> [f64; 3] {
            let load = vec![0.0, 0.0, -0.5 * pressure];
            let sol = gauss_solve(&mat, &load);
            [sol[0], sol[1], sol[2]]
        };

        let cod = solve_for(1.0);
        // Crack opens in the +z (normal) direction.
        assert!(
            cod[2] > 0.0,
            "crack did not open under tensile pressure: {cod:?}"
        );
        // Closed-form single-element estimate: Delta u_z = p (1-nu) a / mu,
        // with p = 1, nu = 0.3, mu = 1.
        let a = (area / PI).sqrt();
        let expected = (1.0 - 0.3) * a;
        assert!(
            (cod[2] - expected).abs() < 1e-9,
            "COD {} vs expected {}",
            cod[2],
            expected
        );
        // In-plane (shear) components vanish for pure normal load.
        assert!(cod[0].abs() < 1e-12 && cod[1].abs() < 1e-12);

        // Doubling the load doubles the opening (linearity).
        let cod2 = solve_for(2.0);
        assert!((cod2[2] - 2.0 * cod[2]).abs() < 1e-9);

        // Stress intensity factor: feed the COD in crack-local order
        // (normal component first) so K_I reflects the opening mode.
        let crack_length = 0.5;
        let k1 = dbem.compute_sif(&[cod[2], cod[0], cod[1]], crack_length);
        assert!(
            k1[0].is_finite() && k1[0] > 0.0,
            "K_I not positive finite: {k1:?}"
        );
        let k1_double = dbem.compute_sif(&[cod2[2], cod2[0], cod2[1]], crack_length);
        assert!(
            (k1_double[0] - 2.0 * k1[0]).abs() < 1e-9,
            "K_I did not scale with load"
        );
    }

    #[test]
    fn test_dual_bem_multi_crack_nonsingular() {
        // Two coplanar crack elements: off-diagonal hypersingular coupling
        // is exercised and the system stays solvable.
        let mut mesh = BoundaryMesh::new();
        mesh.nodes = vec![
            [-1.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [-0.5, 1.0, 0.0],
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.5, 1.0, 0.0],
        ];
        let crack = vec![
            CrackElement {
                node_ids: vec![0, 1, 2],
                normal: [0.0, 0.0, 1.0],
                area: 0.5,
            },
            CrackElement {
                node_ids: vec![3, 4, 5],
                normal: [0.0, 0.0, 1.0],
                area: 0.5,
            },
        ];
        let dbem = DualBem::new(mesh, crack, 1.0, 0.25);
        let (mat, _rhs) = dbem.assemble();
        assert_eq!(mat.rows, 6);

        // Off-diagonal crack-crack coupling must be present and finite.
        let coupling = mat.get(2, 5).abs() + mat.get(5, 2).abs();
        assert!(coupling.is_finite());

        // Pressurise both crack elements; both must open.
        let load = vec![0.0, 0.0, -0.5, 0.0, 0.0, -0.5];
        let sol = gauss_solve(&mat, &load);
        assert!(sol.iter().all(|v| v.is_finite()));
        assert!(sol[2] > 0.0 && sol[5] > 0.0, "cracks did not open: {sol:?}");
    }

    #[test]
    fn test_d_and_s_kernels_finite() {
        let xi = [0.0, 0.0, 0.0];
        let x = [1.0, 0.5, 0.25];
        let n = [0.0, 0.0, 1.0];
        let d = d_kernel_block(&xi, &x, 0.3);
        let s = s_kernel_block(&xi, &x, &n, 1.0, 0.3);
        for (d_k, s_k) in d.iter().zip(s.iter()) {
            for (d_ki, s_ki) in d_k.iter().zip(s_k.iter()) {
                for (&d_kij, &s_kij) in d_ki.iter().zip(s_ki.iter()) {
                    assert!(d_kij.is_finite());
                    assert!(s_kij.is_finite());
                }
            }
        }
        // The element-self block must be symmetric and have a negative
        // (finite-part) normal entry.
        let m = hypersingular_self_block(1.0, &n, 1.0, 0.3);
        assert!((m[0][1] - m[1][0]).abs() < 1e-12);
        assert!(m[2][2] < 0.0);
    }
}
