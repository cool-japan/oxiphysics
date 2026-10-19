//! Regression tests for the covariance eigen-solver behind
//! [`estimate_normals`] (findings F2268 / F2269).
//!
//! The previous solver started its power iteration at `[1, 0, 0]` and shifted
//! the matrix by the absolute amount `max_eigenvalue + 1`.  The shift dwarfs
//! the covariance of any cloud whose spacing is below a metre, so the iteration
//! barely moved from its start vector, and an axis-aligned start vector was
//! never left at all.

use super::functions::smallest_eigenvector_3x3;
use super::{PointCloud, estimate_normals};

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn norm(a: [f64; 3]) -> f64 {
    dot(a, a).sqrt()
}

fn unit(a: [f64; 3]) -> [f64; 3] {
    let n = norm(a);
    [a[0] / n, a[1] / n, a[2] / n]
}

/// Rotation matrix about the (not necessarily unit) `axis` by `angle` radians
/// (Rodrigues' formula).
fn rotation(axis: [f64; 3], angle: f64) -> [[f64; 3]; 3] {
    let [x, y, z] = unit(axis);
    let (s, c) = angle.sin_cos();
    let t = 1.0 - c;
    [
        [t * x * x + c, t * x * y - s * z, t * x * z + s * y],
        [t * x * y + s * z, t * y * y + c, t * y * z - s * x],
        [t * x * z - s * y, t * y * z + s * x, t * z * z + c],
    ]
}

fn apply(m: [[f64; 3]; 3], v: [f64; 3]) -> [f64; 3] {
    [dot(m[0], v), dot(m[1], v), dot(m[2], v)]
}

fn column(m: [[f64; 3]; 3], c: usize) -> [f64; 3] {
    [m[0][c], m[1][c], m[2][c]]
}

/// `R · diag(l) · Rᵀ`, a symmetric matrix with eigenvalues `l` and the columns
/// of `R` as eigenvectors.
fn symmetric_from_eigen(r: [[f64; 3]; 3], l: [f64; 3]) -> [[f64; 3]; 3] {
    let mut m = [[0.0; 3]; 3];
    for (i, row) in m.iter_mut().enumerate() {
        for (j, cell) in row.iter_mut().enumerate() {
            *cell = (0..3).map(|k| r[i][k] * l[k] * r[j][k]).sum();
        }
    }
    m
}

/// Regular `nu x nv` grid of points `origin + i*step_u*u + j*step_v*v`.
fn grid_cloud(
    origin: [f64; 3],
    (u, step_u, nu): ([f64; 3], f64, usize),
    (v, step_v, nv): ([f64; 3], f64, usize),
) -> PointCloud {
    let mut cloud = PointCloud::new();
    for i in 0..nu {
        for j in 0..nv {
            let a = i as f64 * step_u;
            let b = j as f64 * step_v;
            cloud.positions.push([
                origin[0] + a * u[0] + b * v[0],
                origin[1] + a * u[1] + b * v[1],
                origin[2] + a * u[2] + b * v[2],
            ]);
        }
    }
    cloud
}

// -- the eigen-solver itself ---------------------------------------------------

#[test]
fn smallest_eigenvector_of_diag_1_4_0_is_z() {
    // The case from the audit: an axis-aligned start vector stayed on x.
    let v = smallest_eigenvector_3x3([[1.0, 0.0, 0.0], [0.0, 4.0, 0.0], [0.0, 0.0, 0.0]]);
    assert!((norm(v) - 1.0).abs() < 1e-12, "not a unit vector: {v:?}");
    assert!((v[2].abs() - 1.0).abs() < 1e-12, "expected +-z, got {v:?}");
    assert!(v[0].abs() < 1e-12 && v[1].abs() < 1e-12, "got {v:?}");
}

#[test]
fn smallest_eigenvector_follows_the_smallest_eigenvalue_on_every_axis() {
    // (0, 1, 4) assigned to the axes in all six orders.
    let permutations: [[usize; 3]; 6] = [
        [0, 1, 2],
        [0, 2, 1],
        [1, 0, 2],
        [1, 2, 0],
        [2, 0, 1],
        [2, 1, 0],
    ];
    let values = [0.0, 1.0, 4.0];
    for perm in permutations {
        let mut diag = [0.0; 3];
        for (rank, &axis) in perm.iter().enumerate() {
            diag[axis] = values[rank];
        }
        let m = [
            [diag[0], 0.0, 0.0],
            [0.0, diag[1], 0.0],
            [0.0, 0.0, diag[2]],
        ];
        let v = smallest_eigenvector_3x3(m);
        // `perm[0]` is the axis that carries the eigenvalue 0.
        assert!(
            (v[perm[0]].abs() - 1.0).abs() < 1e-12,
            "diag {diag:?}: expected the axis {} but got {v:?}",
            perm[0]
        );
    }
}

#[test]
fn smallest_eigenvector_is_independent_of_the_matrix_scale() {
    let r = rotation([1.0, 2.0, 3.0], 30.0_f64.to_radians());
    let expected = column(r, 0); // eigenvector of the smallest eigenvalue 0.5
    for exponent in [-12, -8, -4, 0, 4, 12] {
        let scale = 10.0_f64.powi(exponent);
        let m = symmetric_from_eigen(r, [0.5 * scale, 2.0 * scale, 3.0 * scale]);
        let v = smallest_eigenvector_3x3(m);
        let alignment = dot(v, expected).abs();
        assert!(
            alignment > 1.0 - 1e-9,
            "scale 1e{exponent}: |v . expected| = {alignment}, v = {v:?}"
        );
    }
}

#[test]
fn smallest_eigenvector_of_a_tiny_nearly_isotropic_covariance() {
    // eigenvalues 1.00e-4, 1.01e-4 and 1.30e-4 differ by ~1 %: the old shifted
    // power iteration left the start vector (1, 1, 1)/sqrt(3) almost untouched.
    let r = rotation([0.3, -1.0, 0.7], 41.0_f64.to_radians());
    let m = symmetric_from_eigen(r, [1.30e-4, 1.01e-4, 1.00e-4]);
    let v = smallest_eigenvector_3x3(m);
    assert!(dot(v, column(r, 2)).abs() > 1.0 - 1e-9, "got {v:?}");
}

#[test]
fn smallest_eigenvector_of_the_zero_matrix_is_still_a_unit_vector() {
    let v = smallest_eigenvector_3x3([[0.0; 3]; 3]);
    assert!((norm(v) - 1.0).abs() < 1e-12, "got {v:?}");
}

// -- estimate_normals ----------------------------------------------------------

#[test]
fn normals_of_an_axis_aligned_patch_are_z_at_every_length_scale() {
    for spacing in [1e-3, 1.0, 1e3] {
        let cloud = grid_cloud(
            [5.0 * spacing, -3.0 * spacing, 7.0 * spacing],
            ([1.0, 0.0, 0.0], spacing, 9),
            ([0.0, 1.0, 0.0], spacing, 9),
        );
        let normals = estimate_normals(&cloud, 8);
        assert_eq!(normals.len(), 81);
        for (i, n) in normals.iter().enumerate() {
            assert!(
                (n[2].abs() - 1.0).abs() < 1e-9,
                "spacing {spacing}: normal {i} is {n:?}"
            );
        }
    }
}

#[test]
fn normals_of_an_anisotropic_patch_are_z() {
    // x-spread 1e-2 and y-spread 4e-2 give a covariance of order 1e-4..1e-3
    // with an exactly zero z eigenvalue.
    let cloud = grid_cloud(
        [0.0, 0.0, 2.0],
        ([1.0, 0.0, 0.0], 1e-2, 12),
        ([0.0, 1.0, 0.0], 4e-2, 12),
    );
    // 20 neighbours reach the adjacent rows, so the neighbourhoods are not
    // collinear (a collinear neighbourhood has no unique plane normal).
    let normals = estimate_normals(&cloud, 20);
    for (i, n) in normals.iter().enumerate() {
        assert!((n[2].abs() - 1.0).abs() < 1e-9, "normal {i} is {n:?}");
    }
}

#[test]
fn normals_of_a_rotated_patch_match_the_analytic_normal() {
    let r = rotation([1.0, 1.0, 0.0], 30.0_f64.to_radians());
    let u = apply(r, [1.0, 0.0, 0.0]);
    let v = apply(r, [0.0, 1.0, 0.0]);
    let n_true = apply(r, [0.0, 0.0, 1.0]);
    // Fixture check: the axis lies in the xy-plane, so the plane normal is the
    // z axis tilted by exactly 30 degrees (far from axis-aligned).
    assert!((n_true[2] - 30.0_f64.to_radians().cos()).abs() < 1e-12);
    for spacing in [1e-3, 1.0, 1e3] {
        let cloud = grid_cloud(
            [spacing, 2.0 * spacing, -spacing],
            (u, spacing, 9),
            (v, spacing, 9),
        );
        let normals = estimate_normals(&cloud, 8);
        for (i, n) in normals.iter().enumerate() {
            assert!((norm(*n) - 1.0).abs() < 1e-9, "normal {i} not unit: {n:?}");
            let alignment = dot(*n, n_true).abs();
            assert!(
                alignment > 0.9999,
                "spacing {spacing}: normal {i} = {n:?}, |n . n_true| = {alignment}"
            );
        }
    }
}

#[test]
fn normals_of_a_unit_sphere_point_outwards() {
    // Fibonacci lattice: deterministic, near-uniform samples of the unit sphere.
    let count = 600;
    let golden_angle = std::f64::consts::PI * (3.0 - 5.0_f64.sqrt());
    let mut cloud = PointCloud::new();
    for i in 0..count {
        let y = 1.0 - 2.0 * (i as f64 + 0.5) / count as f64;
        let ring = (1.0 - y * y).sqrt();
        let phi = golden_angle * i as f64;
        cloud
            .positions
            .push([ring * phi.cos(), y, ring * phi.sin()]);
    }
    let normals = estimate_normals(&cloud, 10);
    for (p, n) in cloud.positions.iter().zip(&normals) {
        assert!((norm(*n) - 1.0).abs() < 1e-9, "not a unit normal: {n:?}");
        let outward = dot(*n, *p);
        assert!(outward > 0.95, "normal {n:?} at {p:?}: n . p = {outward}");
    }
}
