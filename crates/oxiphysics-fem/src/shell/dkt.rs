// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Discrete Kirchhoff triangle (DKT) plate-bending element.
//!
//! This is the 9-DOF thin-plate triangle of Batoz, Bathe and Ho (1980),
//! *A study of three-node triangular plate bending elements*, Int. J. Numer.
//! Meth. Engng 15, 1771-1812.  It is the bending part of
//! [`FlatShellElement`](super::FlatShellElement).
//!
//! # Formulation
//!
//! The plate deflection `w` is never interpolated inside the element.  Instead
//! the *slope* field `φ = (∂w/∂x, ∂w/∂y)` is interpolated with the quadratic
//! shape functions `N_1 … N_6` of the six-node triangle (three vertices, three
//! mid-side nodes), and the Kirchhoff condition "the rotation equals the slope
//! of the deflection" is imposed discretely at the six nodes:
//!
//! * at the vertices `φ` is the nodal slope, taken from the rotation DOFs;
//! * along every side `w` is a cubic Hermite polynomial of the arc length
//!   (end values and end tangential slopes), which fixes the *tangential* slope
//!   at the mid-side node,
//!   `φ_s = 3 (w_j − w_i) / (2 l) − (φ_s,i + φ_s,j) / 4`;
//! * the *normal* slope varies linearly along the side, so at the mid-side node
//!   `φ_n = (φ_n,i + φ_n,j) / 2`.
//!
//! The mid-side slopes are therefore linear functions of the nine DOFs
//! `(w, θx, θy)` per vertex, and so is the whole slope field.  The curvatures
//! `κ = (∂φ_x/∂x, ∂φ_y/∂y, ∂φ_x/∂y + ∂φ_y/∂x)` follow by differentiating the
//! quadratic interpolation, and `K = ∫ Bᵀ D_b B dA` is integrated exactly with
//! the three-point Hammer rule (`B` is linear in the area coordinates).
//!
//! # Conventions
//!
//! Degrees of freedom per vertex are `(w, θx, θy) = (w, ∂w/∂y, −∂w/∂x)`, the
//! right-handed rotations about the x- and y-axes, matching
//! [`RectPlateElement`](super::RectPlateElement) and the bending DOFs of
//! [`FlatShellElement`](super::FlatShellElement).  The element is orientation
//! independent: clockwise and counter-clockwise vertex orderings give the same
//! physical element.
//!
//! The element has exactly three zero-energy modes (the rigid motions
//! `w = c₀ + c₁ x + c₂ y`), reproduces constant curvature states exactly (patch
//! test) and converges with `O(h²)` in the deflection.

/// Row `a` holds the dependence of one slope component (`∂w/∂x` or `∂w/∂y`) at
/// node `a` of the six-node interpolation on the nine element DOFs.
type SlopeRows = [[f64; 9]; 6];

/// Three-point Hammer rule on the reference triangle `(ξ, η)`; every point
/// carries one third of the element area.  Exact for polynomials of degree 2.
const HAMMER_POINTS: [(f64, f64); 3] = [
    (1.0 / 6.0, 1.0 / 6.0),
    (2.0 / 3.0, 1.0 / 6.0),
    (1.0 / 6.0, 2.0 / 3.0),
];

/// The three sides as `(start, end)` vertex pairs.  Side `k` carries the
/// mid-side node `3 + k`, the node opposite vertex `k`.
const SIDES: [(usize, usize); 3] = [(1, 2), (2, 0), (0, 1)];

/// Twice the signed area of the triangle, or `None` for a degenerate (zero
/// area, collapsed or non-finite) element.
fn twice_signed_area(nodes: &[[f64; 2]; 3]) -> Option<f64> {
    let two_a = (nodes[1][0] - nodes[0][0]) * (nodes[2][1] - nodes[0][1])
        - (nodes[2][0] - nodes[0][0]) * (nodes[1][1] - nodes[0][1]);
    let longest_side_sq = SIDES
        .iter()
        .map(|&(i, j)| {
            let dx = nodes[j][0] - nodes[i][0];
            let dy = nodes[j][1] - nodes[i][1];
            dx * dx + dy * dy
        })
        .fold(0.0_f64, f64::max);
    if two_a.is_finite() && two_a.abs() > 1.0e-14 * longest_side_sq {
        Some(two_a)
    } else {
        None
    }
}

/// Derivatives `(∂N_a/∂ξ, ∂N_a/∂η)` of the six quadratic shape functions
/// `N = (L₁(2L₁−1), L₂(2L₂−1), L₃(2L₃−1), 4L₂L₃, 4L₃L₁, 4L₁L₂)` with the area
/// coordinates `L₁ = 1−ξ−η`, `L₂ = ξ`, `L₃ = η`.
fn shape_gradients(xi: f64, eta: f64) -> ([f64; 6], [f64; 6]) {
    let l1 = 1.0 - xi - eta;
    let d_xi = [
        1.0 - 4.0 * l1,
        4.0 * xi - 1.0,
        0.0,
        4.0 * eta,
        -4.0 * eta,
        4.0 * (1.0 - 2.0 * xi - eta),
    ];
    let d_eta = [
        1.0 - 4.0 * l1,
        0.0,
        4.0 * eta - 1.0,
        4.0 * xi,
        4.0 * (1.0 - xi - 2.0 * eta),
        -4.0 * xi,
    ];
    (d_xi, d_eta)
}

/// Linear maps from the nine element DOFs to the slope components
/// `(∂w/∂x, ∂w/∂y)` at the six nodes, with the discrete Kirchhoff constraints
/// applied at the mid-side nodes.  `None` if a side has zero length.
fn slope_map(nodes: &[[f64; 2]; 3]) -> Option<(SlopeRows, SlopeRows)> {
    let mut slope_x: SlopeRows = [[0.0; 9]; 6];
    let mut slope_y: SlopeRows = [[0.0; 9]; 6];
    for vertex in 0..3 {
        // θy = −∂w/∂x and θx = ∂w/∂y.
        slope_x[vertex][3 * vertex + 2] = -1.0;
        slope_y[vertex][3 * vertex + 1] = 1.0;
    }
    for (side, &(i, j)) in SIDES.iter().enumerate() {
        let dx = nodes[j][0] - nodes[i][0];
        let dy = nodes[j][1] - nodes[i][1];
        let length = dx.hypot(dy);
        if length.is_nan() || length <= 0.0 {
            return None;
        }
        // Unit tangent s (i -> j) and unit normal n; the result is invariant
        // under the sign of either.
        let (sx, sy) = (dx / length, dy / length);
        let (nx, ny) = (-sy, sx);
        let mut tangential = [0.0_f64; 9];
        let mut normal = [0.0_f64; 9];
        for col in 0..9 {
            let sum_x = slope_x[i][col] + slope_x[j][col];
            let sum_y = slope_y[i][col] + slope_y[j][col];
            // φ_s = 3 (w_j − w_i) / (2 l) − (φ_s,i + φ_s,j) / 4
            tangential[col] = -0.25 * (sx * sum_x + sy * sum_y);
            // φ_n = (φ_n,i + φ_n,j) / 2
            normal[col] = 0.5 * (nx * sum_x + ny * sum_y);
        }
        tangential[3 * j] += 1.5 / length;
        tangential[3 * i] -= 1.5 / length;
        let node = 3 + side;
        for col in 0..9 {
            slope_x[node][col] = sx * tangential[col] + nx * normal[col];
            slope_y[node][col] = sy * tangential[col] + ny * normal[col];
        }
    }
    Some((slope_x, slope_y))
}

/// Curvature-DOF matrix `B (3 × 9)` at the reference point `(ξ, η)` for a
/// non-degenerate triangle with the given slope maps.
fn curvature_rows(
    nodes: &[[f64; 2]; 3],
    two_a: f64,
    slope_x: &SlopeRows,
    slope_y: &SlopeRows,
    xi: f64,
    eta: f64,
) -> [[f64; 9]; 3] {
    let (d_xi, d_eta) = shape_gradients(xi, eta);
    let (dx21, dy21) = (nodes[1][0] - nodes[0][0], nodes[1][1] - nodes[0][1]);
    let (dx31, dy31) = (nodes[2][0] - nodes[0][0], nodes[2][1] - nodes[0][1]);
    let mut b = [[0.0_f64; 9]; 3];
    for node in 0..6 {
        // Physical gradient of N_a from the inverse Jacobian of the affine map.
        let dn_dx = (dy31 * d_xi[node] - dy21 * d_eta[node]) / two_a;
        let dn_dy = (-dx31 * d_xi[node] + dx21 * d_eta[node]) / two_a;
        for col in 0..9 {
            b[0][col] += dn_dx * slope_x[node][col];
            b[1][col] += dn_dy * slope_y[node][col];
            b[2][col] += dn_dy * slope_x[node][col] + dn_dx * slope_y[node][col];
        }
    }
    b
}

/// Curvature-displacement matrix `B (3 × 9)` of the DKT element at the point
/// with area-coordinate parameters `(ξ, η)` (`L₂ = ξ`, `L₃ = η`, `L₁ = 1−ξ−η`).
///
/// The curvature vector is `κ = B d = (∂²w/∂x², ∂²w/∂y², 2 ∂²w/∂x∂y)` for the
/// nodal DOF vector `d = (w₁, θx₁, θy₁, w₂, θx₂, θy₂, w₃, θx₃, θy₃)` with
/// `(θx, θy) = (∂w/∂y, −∂w/∂x)`.  Moments follow as `M = D_b κ`.
///
/// Returns `None` for a degenerate triangle.
pub fn dkt_curvature_matrix(nodes: &[[f64; 2]; 3], xi: f64, eta: f64) -> Option<[[f64; 9]; 3]> {
    let two_a = twice_signed_area(nodes)?;
    let (slope_x, slope_y) = slope_map(nodes)?;
    Some(curvature_rows(nodes, two_a, &slope_x, &slope_y, xi, eta))
}

/// 9 × 9 bending stiffness `K = ∫ Bᵀ D_b B dA` of the DKT element.
///
/// `nodes` are the three vertex positions in the plane; `d_b` is the plate
/// constitutive matrix relating the curvature vector
/// `(∂²w/∂x², ∂²w/∂y², 2 ∂²w/∂x∂y)` to the moments, e.g.
/// `D [[1, ν, 0], [ν, 1, 0], [0, 0, (1 − ν)/2]]` with the flexural rigidity
/// `D = E t³ / (12 (1 − ν²))`.
///
/// The matrix is symmetric, positive semi-definite and has rank six: its three
/// null vectors are the rigid motions of the plane.  DOF order per vertex is
/// `(w, θx, θy) = (w, ∂w/∂y, −∂w/∂x)`.  A degenerate triangle gives the zero
/// matrix.
pub fn dkt_bending_stiffness(nodes: &[[f64; 2]; 3], d_b: &[[f64; 3]; 3]) -> [[f64; 9]; 9] {
    let mut k = [[0.0_f64; 9]; 9];
    let (Some(two_a), Some((slope_x, slope_y))) = (twice_signed_area(nodes), slope_map(nodes))
    else {
        return k;
    };
    let weight = two_a.abs() / 6.0;
    for &(xi, eta) in &HAMMER_POINTS {
        let b = curvature_rows(nodes, two_a, &slope_x, &slope_y, xi, eta);
        let mut db = [[0.0_f64; 9]; 3];
        for (db_row, d_row) in db.iter_mut().zip(d_b.iter()) {
            for (col, entry) in db_row.iter_mut().enumerate() {
                *entry = d_row
                    .iter()
                    .zip(b.iter())
                    .map(|(d_value, b_row)| d_value * b_row[col])
                    .sum();
            }
        }
        for (i, k_row) in k.iter_mut().enumerate() {
            for (j, entry) in k_row.iter_mut().enumerate() {
                *entry += weight * (0..3).map(|r| b[r][i] * db[r][j]).sum::<f64>();
            }
        }
    }
    // Remove the rounding asymmetry of the triple products.
    let mut symmetric = [[0.0_f64; 9]; 9];
    for (i, symmetric_row) in symmetric.iter_mut().enumerate() {
        for (j, entry) in symmetric_row.iter_mut().enumerate() {
            *entry = 0.5 * (k[i][j] + k[j][i]);
        }
    }
    symmetric
}

/// Equivalent nodal load vector of a uniform transverse pressure `q`
/// (positive along `+w`) on the triangle, in DKT DOF order.
///
/// The deflection of the DKT element is not interpolated inside the element, so
/// no consistent load exists; the standard choice (Batoz et al., 1980) lumps the
/// resultant `q A` equally to the three deflection DOFs and applies no nodal
/// moments.
pub fn dkt_uniform_pressure_load(nodes: &[[f64; 2]; 3], q: f64) -> [f64; 9] {
    let mut f = [0.0_f64; 9];
    if let Some(two_a) = twice_signed_area(nodes) {
        let share = q * 0.5 * two_a.abs() / 3.0;
        for vertex in 0..3 {
            f[3 * vertex] = share;
        }
    }
    f
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::eigen::jacobi_eigen_dense;
    use crate::solvers::CholeskyDense;

    /// `D_b` for `E t³ / (12 (1 − ν²)) = 1`-scaled unit plates.
    fn plate_matrix(e: f64, nu: f64, t: f64) -> [[f64; 3]; 3] {
        let d = e * t * t * t / (12.0 * (1.0 - nu * nu));
        [
            [d, d * nu, 0.0],
            [d * nu, d, 0.0],
            [0.0, 0.0, d * (1.0 - nu) / 2.0],
        ]
    }

    /// Triangles of different shape and orientation (one is clockwise).
    fn sample_triangles() -> Vec<[[f64; 2]; 3]> {
        vec![
            [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]],
            [[0.0, 0.0], [1.0, 0.2], [0.3, 0.9]],
            [[0.2, -0.1], [2.5, 0.4], [1.1, 0.35]],
            [[0.0, 0.0], [0.0, 1.0], [1.0, 0.0]],
            [[-0.5, 0.2], [0.4, -0.3], [0.1, 0.8]],
        ]
    }

    fn quad_form(k: &[[f64; 9]; 9], d: &[f64; 9]) -> f64 {
        k.iter()
            .zip(d.iter())
            .map(|(row, di)| {
                0.5 * di
                    * row
                        .iter()
                        .zip(d.iter())
                        .map(|(kij, dj)| kij * dj)
                        .sum::<f64>()
            })
            .sum()
    }

    fn triangle_area(nodes: &[[f64; 2]; 3]) -> f64 {
        0.5 * ((nodes[1][0] - nodes[0][0]) * (nodes[2][1] - nodes[0][1])
            - (nodes[2][0] - nodes[0][0]) * (nodes[1][1] - nodes[0][1]))
            .abs()
    }

    /// Nodal DOFs of `w = c0 + c1 x + c2 y`.
    fn rigid_dofs(nodes: &[[f64; 2]; 3], c: [f64; 3]) -> [f64; 9] {
        let mut d = [0.0; 9];
        for (v, p) in nodes.iter().enumerate() {
            d[3 * v] = c[0] + c[1] * p[0] + c[2] * p[1];
            d[3 * v + 1] = c[2];
            d[3 * v + 2] = -c[1];
        }
        d
    }

    /// Nodal DOFs of `w = (a x² + b y²)/2 + c x y`.
    fn quadratic_dofs(nodes: &[[f64; 2]; 3], a: f64, b: f64, c: f64) -> [f64; 9] {
        let mut d = [0.0; 9];
        for (v, p) in nodes.iter().enumerate() {
            d[3 * v] = 0.5 * (a * p[0] * p[0] + b * p[1] * p[1]) + c * p[0] * p[1];
            d[3 * v + 1] = b * p[1] + c * p[0];
            d[3 * v + 2] = -(a * p[0] + c * p[1]);
        }
        d
    }

    #[test]
    fn test_dkt_stiffness_is_symmetric_with_exactly_three_zero_modes() {
        let d_b = plate_matrix(210e9, 0.3, 0.02);
        for nodes in sample_triangles() {
            let k = dkt_bending_stiffness(&nodes, &d_b);
            let scale = k.iter().flatten().fold(0.0_f64, |m, v| m.max(v.abs()));
            for (i, row) in k.iter().enumerate() {
                for (j, &entry) in row.iter().enumerate() {
                    assert_eq!(entry, k[j][i], "asymmetry at ({i},{j})");
                }
            }
            let rows: Vec<Vec<f64>> = k.iter().map(|r| r.to_vec()).collect();
            let (eigenvalues, _) = jacobi_eigen_dense(&rows, 9);
            let lambda_max = eigenvalues[8];
            assert!(
                eigenvalues[..3].iter().all(|l| l.abs() < 1e-9 * lambda_max),
                "three rigid modes expected: {eigenvalues:?}"
            );
            assert!(
                eigenvalues[3..].iter().all(|&l| l > 1e-4 * lambda_max),
                "rank must be six: {eigenvalues:?} (scale {scale})"
            );
        }
    }

    #[test]
    fn test_dkt_rigid_motions_carry_no_force() {
        let d_b = plate_matrix(1.0, 0.3, 1.0);
        for nodes in sample_triangles() {
            let k = dkt_bending_stiffness(&nodes, &d_b);
            let k_norm = k.iter().flatten().fold(0.0_f64, |m, v| m.max(v.abs()));
            for c in [
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [0.0, 0.0, 1.0],
                [0.7, -1.3, 0.4],
            ] {
                let d = rigid_dofs(&nodes, c);
                for row in &k {
                    let force: f64 = row.iter().zip(d.iter()).map(|(a, b)| a * b).sum();
                    assert!(force.abs() < 1e-12 * k_norm, "mode {c:?}: force {force}");
                }
            }
        }
    }

    #[test]
    fn test_dkt_constant_curvature_patch_test() {
        let (e, nu, t) = (2.0, 0.25, 0.3);
        let d_b = plate_matrix(e, nu, t);
        for nodes in sample_triangles() {
            let k = dkt_bending_stiffness(&nodes, &d_b);
            for (a, b, c) in [
                (1.0, 0.0, 0.0),
                (0.0, 1.0, 0.0),
                (0.0, 0.0, 1.0),
                (0.7, -0.4, 0.25),
            ] {
                let d = quadratic_dofs(&nodes, a, b, c);
                // Exact curvature vector (w,xx, w,yy, 2 w,xy) = (a, b, 2c).
                let kappa = [a, b, 2.0 * c];
                let mut energy = 0.0;
                for i in 0..3 {
                    for j in 0..3 {
                        energy += 0.5 * kappa[i] * d_b[i][j] * kappa[j];
                    }
                }
                let exact = energy * triangle_area(&nodes);
                let fe = quad_form(&k, &d);
                assert!(
                    (fe - exact).abs() < 1e-12 * exact.max(1e-30),
                    "{nodes:?} ({a},{b},{c}): {fe} vs {exact}"
                );
                // Pointwise: B d equals the exact constant curvature everywhere.
                for (xi, eta) in [(0.2, 0.3), (1.0 / 6.0, 2.0 / 3.0), (0.0, 0.0), (0.5, 0.5)] {
                    let bm = dkt_curvature_matrix(&nodes, xi, eta).expect("valid triangle");
                    for r in 0..3 {
                        let got: f64 = (0..9).map(|c9| bm[r][c9] * d[c9]).sum();
                        assert!(
                            (got - kappa[r]).abs() < 1e-11,
                            "curvature {r} at ({xi},{eta}): {got} vs {}",
                            kappa[r]
                        );
                    }
                }
            }
        }
    }

    #[test]
    fn test_dkt_is_independent_of_vertex_orientation_and_numbering() {
        let d_b = plate_matrix(1.0, 0.3, 1.0);
        let nodes = [[0.0, 0.0], [1.0, 0.2], [0.3, 0.9]];
        let k = dkt_bending_stiffness(&nodes, &d_b);
        // Same physical triangle with vertices listed (0, 2, 1) and (1, 2, 0).
        for perm in [[0usize, 2, 1], [1, 2, 0], [2, 1, 0]] {
            let permuted = [nodes[perm[0]], nodes[perm[1]], nodes[perm[2]]];
            let kp = dkt_bending_stiffness(&permuted, &d_b);
            for a in 0..9 {
                for b in 0..9 {
                    let (va, ca) = (perm[a / 3], a % 3);
                    let (vb, cb) = (perm[b / 3], b % 3);
                    let reference = k[3 * va + ca][3 * vb + cb];
                    assert!(
                        (kp[a][b] - reference).abs() < 1e-12 * (1.0 + reference.abs()),
                        "perm {perm:?} entry ({a},{b}): {} vs {reference}",
                        kp[a][b]
                    );
                }
            }
        }
    }

    #[test]
    fn test_dkt_energy_is_invariant_under_in_plane_rotation() {
        let d_b = plate_matrix(1.0, 0.3, 1.0);
        let nodes = [[0.0, 0.0], [1.0, 0.2], [0.3, 0.9]];
        let alpha: f64 = 0.7;
        let (s, c) = alpha.sin_cos();
        let rotate = |p: [f64; 2]| [c * p[0] - s * p[1], s * p[0] + c * p[1]];
        let rotated = [rotate(nodes[0]), rotate(nodes[1]), rotate(nodes[2])];
        let k = dkt_bending_stiffness(&nodes, &d_b);
        let kr = dkt_bending_stiffness(&rotated, &d_b);
        // A generic displacement state; (θx, θy) rotates as a vector.
        let d: [f64; 9] = [0.3, -0.2, 0.5, 0.1, 0.4, -0.3, -0.25, 0.15, 0.05];
        let mut dr = d;
        for v in 0..3 {
            let [tx, ty] = rotate([d[3 * v + 1], d[3 * v + 2]]);
            dr[3 * v + 1] = tx;
            dr[3 * v + 2] = ty;
        }
        let u = quad_form(&k, &d);
        let ur = quad_form(&kr, &dr);
        assert!((u - ur).abs() < 1e-12 * u.abs(), "{u} vs {ur}");
    }

    #[test]
    fn test_dkt_degenerate_triangle_gives_zero_stiffness() {
        let d_b = plate_matrix(1.0, 0.3, 1.0);
        for nodes in [
            [[0.0, 0.0], [1.0, 0.0], [2.0, 0.0]],
            [[0.0, 0.0], [0.0, 0.0], [1.0, 1.0]],
            [[1.0, 1.0], [1.0, 1.0], [1.0, 1.0]],
        ] {
            let k = dkt_bending_stiffness(&nodes, &d_b);
            assert!(k.iter().flatten().all(|v| *v == 0.0));
            assert!(dkt_curvature_matrix(&nodes, 0.2, 0.2).is_none());
            assert!(
                dkt_uniform_pressure_load(&nodes, 1.0)
                    .iter()
                    .all(|v| *v == 0.0)
            );
        }
    }

    #[test]
    fn test_dkt_pressure_load_resultant_and_distribution() {
        let nodes = [[0.0, 0.0], [2.0, 0.0], [0.0, 3.0]];
        let f = dkt_uniform_pressure_load(&nodes, 5.0);
        assert!(
            (f[0] - 5.0).abs() < 1e-14 && (f[3] - 5.0).abs() < 1e-14 && (f[6] - 5.0).abs() < 1e-14
        );
        assert_eq!(
            (f[1], f[2], f[4], f[5], f[7], f[8]),
            (0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
        );
        let flipped = dkt_uniform_pressure_load(&[nodes[0], nodes[2], nodes[1]], 5.0);
        assert!((flipped[0] - 5.0).abs() < 1e-14);
    }

    /// `w_c D / (q a⁴)` at the centre of a uniformly loaded square plate of
    /// side `a = 1`, from a quarter-plate model (symmetry planes `x = 0`,
    /// `y = 0`) of `cells × cells` squares, each cut into two triangles with
    /// alternating diagonals, simply supported or clamped on the outer edges.
    fn quarter_plate_centre_deflection(cells: usize, clamped: bool) -> f64 {
        let d_b = plate_matrix(1.0, 0.3, 1.0);
        let flexural_rigidity = d_b[0][0];
        let h = 0.5 / cells as f64;
        let nn = cells + 1;
        let node_id = |i: usize, j: usize| j * nn + i;
        let position = |n: usize| [(n % nn) as f64 * h, (n / nn) as f64 * h];
        let ndof = 3 * nn * nn;
        let mut k = vec![0.0_f64; ndof * ndof];
        let mut f = vec![0.0_f64; ndof];
        for j in 0..cells {
            for i in 0..cells {
                let (n00, n10, n11, n01) = (
                    node_id(i, j),
                    node_id(i + 1, j),
                    node_id(i + 1, j + 1),
                    node_id(i, j + 1),
                );
                let triangles = if (i + j).is_multiple_of(2) {
                    [[n00, n10, n11], [n00, n11, n01]]
                } else {
                    [[n00, n10, n01], [n10, n11, n01]]
                };
                for tri in triangles {
                    let coords = [position(tri[0]), position(tri[1]), position(tri[2])];
                    let ke = dkt_bending_stiffness(&coords, &d_b);
                    let fe = dkt_uniform_pressure_load(&coords, 1.0);
                    let dofs: Vec<usize> = tri
                        .iter()
                        .flat_map(|&n| [3 * n, 3 * n + 1, 3 * n + 2])
                        .collect();
                    for (r, &dr) in dofs.iter().enumerate() {
                        f[dr] += fe[r];
                        for (c, &dc) in dofs.iter().enumerate() {
                            k[dr * ndof + dc] += ke[r][c];
                        }
                    }
                }
            }
        }
        let mut fixed = vec![false; ndof];
        for j in 0..nn {
            for i in 0..nn {
                let n = node_id(i, j);
                if i == 0 {
                    fixed[3 * n + 2] = true; // symmetry: ∂w/∂x = 0
                }
                if j == 0 {
                    fixed[3 * n + 1] = true; // symmetry: ∂w/∂y = 0
                }
                if i == cells || j == cells {
                    fixed[3 * n] = true;
                    if clamped {
                        fixed[3 * n + 1] = true;
                        fixed[3 * n + 2] = true;
                    } else {
                        if i == cells {
                            fixed[3 * n + 1] = true; // edge along y: θx = 0
                        }
                        if j == cells {
                            fixed[3 * n + 2] = true; // edge along x: θy = 0
                        }
                    }
                }
            }
        }
        let free: Vec<usize> = (0..ndof).filter(|&d| !fixed[d]).collect();
        let nf = free.len();
        let mut k_ff = vec![0.0_f64; nf * nf];
        let mut f_f = vec![0.0_f64; nf];
        for (r, &dr) in free.iter().enumerate() {
            f_f[r] = f[dr];
            for (c, &dc) in free.iter().enumerate() {
                k_ff[r * nf + c] = k[dr * ndof + dc];
            }
        }
        let chol = CholeskyDense::new(nf);
        let factor = chol.factorize(&k_ff).expect("constrained plate is SPD");
        let u = chol.solve(&factor, &f_f);
        let centre = 3 * node_id(0, 0);
        let position_in_free = free
            .iter()
            .position(|&d| d == centre)
            .expect("centre deflection is a free DOF");
        u[position_in_free] * flexural_rigidity
    }

    #[test]
    fn test_dkt_clamped_square_plate_matches_reference_deflection() {
        // Clamped, uniform load: w_max = 0.00126 q a⁴ / D (Timoshenko).
        let w8 = quarter_plate_centre_deflection(8, true);
        assert!(
            (w8 - 0.00126).abs() < 0.03 * 0.00126,
            "8x8 quarter-plate centre deflection {w8} vs 0.00126"
        );
        // O(h²) convergence towards the converged value 0.0012653.
        let exact = 0.001_265_3;
        let err2 = (quarter_plate_centre_deflection(2, true) - exact).abs();
        let err4 = (quarter_plate_centre_deflection(4, true) - exact).abs();
        let err8 = (w8 - exact).abs();
        // A 4 x 4 quarter model is an 8 x 8 mesh of the whole plate (2 x 8 x 8
        // = 128 triangles): already within 3 % of the converged deflection.
        assert!(
            err4 < 0.03 * exact,
            "8x8 full-plate centre deflection error {err4} vs {exact}"
        );
        assert!(
            err4 < 0.5 * err2 && err8 < 0.5 * err4,
            "{err2} {err4} {err8}"
        );
    }

    #[test]
    fn test_dkt_simply_supported_square_plate_matches_navier() {
        // Simply supported, uniform load: w_max = 0.00406 q a⁴ / D (Navier).
        let w8 = quarter_plate_centre_deflection(8, false);
        assert!(
            (w8 - 0.00406).abs() < 0.03 * 0.00406,
            "8x8 quarter-plate centre deflection {w8} vs 0.00406"
        );
        let exact = 0.004_062_35;
        let err4 = (quarter_plate_centre_deflection(4, false) - exact).abs();
        let err8 = (w8 - exact).abs();
        assert!(err8 < 0.5 * err4, "no convergence: {err4} -> {err8}");
        // The 4 x 4 quarter model is an 8 x 8 mesh of the whole plate.
        assert!(
            err4 < 0.03 * exact,
            "8x8 full-plate centre deflection error {err4} vs {exact}"
        );
    }
}
