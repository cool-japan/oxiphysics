//! Regression tests for PCA based normal estimation.
//!
//! The normal of a neighbourhood is the eigenvector of the smallest eigenvalue
//! of its covariance matrix. An earlier shifted power iteration started from the
//! x axis and returned an in-plane vector for neighbourhoods lying in an
//! axis-aligned plane.

use super::types::{NormalEstimation, PointCloud};

fn norm(a: [f64; 3]) -> f64 {
    (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt()
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn unit(a: [f64; 3]) -> [f64; 3] {
    let length = norm(a);
    [a[0] / length, a[1] / length, a[2] / length]
}

/// Distance between `normal` and `expected` up to the sign of the normal.
fn normal_error(normal: [f64; 3], expected: [f64; 3]) -> f64 {
    let same = norm([
        normal[0] - expected[0],
        normal[1] - expected[1],
        normal[2] - expected[2],
    ]);
    let flipped = norm([
        normal[0] + expected[0],
        normal[1] + expected[1],
        normal[2] + expected[2],
    ]);
    same.min(flipped)
}

fn mat_vec(m: &[[f64; 3]; 3], v: [f64; 3]) -> [f64; 3] {
    [
        m[0][0] * v[0] + m[0][1] * v[1] + m[0][2] * v[2],
        m[1][0] * v[0] + m[1][1] * v[1] + m[1][2] * v[2],
        m[2][0] * v[0] + m[2][1] * v[1] + m[2][2] * v[2],
    ]
}

/// Grid of `nx * ny` points mapped through `place`.
fn grid_cloud(nx: usize, ny: usize, place: impl Fn(f64, f64) -> [f64; 3]) -> PointCloud {
    let mut cloud = PointCloud::new();
    for i in 0..nx {
        for j in 0..ny {
            cloud.add_point(place(i as f64, j as f64));
        }
    }
    cloud
}

/// Both public entry points, evaluated with a neighbourhood covering the cloud.
fn estimated_normals(cloud: &PointCloud) -> Vec<Vec<[f64; 3]>> {
    let k = cloud.points.len();
    let from_estimator = NormalEstimation::new(k).estimate(cloud);
    let mut pca_cloud = PointCloud::from_points(cloud.points.clone());
    pca_cloud.estimate_normals_pca(k);
    vec![from_estimator, pca_cloud.normals]
}

fn assert_all_normals(cloud: &PointCloud, expected: [f64; 3], tolerance: f64) {
    for normals in estimated_normals(cloud) {
        assert_eq!(normals.len(), cloud.points.len());
        for normal in normals {
            assert!(
                (norm(normal) - 1.0).abs() < 1e-9,
                "normal {normal:?} is not a unit vector"
            );
            let error = normal_error(normal, expected);
            assert!(
                error < tolerance,
                "normal {normal:?} differs from {expected:?} by {error}"
            );
        }
    }
}

#[test]
fn square_grid_in_the_xz_plane_has_the_y_axis_as_normal() {
    // Covariance diag(a, 0, a): the old iteration collapsed to the zero vector.
    let cloud = grid_cloud(5, 5, |x, z| [x, 0.0, z]);
    assert_all_normals(&cloud, [0.0, 1.0, 0.0], 1e-9);
}

#[test]
fn wide_grid_in_the_xz_plane_has_the_y_axis_as_normal() {
    // Covariance diag(a, 0, c) with a > c.
    let cloud = grid_cloud(5, 3, |x, z| [x, 0.0, z]);
    assert_all_normals(&cloud, [0.0, 1.0, 0.0], 1e-9);
}

#[test]
fn grid_in_the_xy_plane_with_wider_y_spread_has_the_z_axis_as_normal() {
    // Covariance diag(a, b, 0) with b > a: the old iteration stayed on the x axis.
    let cloud = grid_cloud(3, 6, |x, y| [x, y, 0.0]);
    assert_all_normals(&cloud, [0.0, 0.0, 1.0], 1e-9);
}

#[test]
fn grid_in_the_yz_plane_has_the_x_axis_as_normal() {
    let cloud = grid_cloud(4, 7, |y, z| [2.0, y, z]);
    assert_all_normals(&cloud, [1.0, 0.0, 0.0], 1e-9);
}

/// Grid points on the plane through `origin` with unit normal `normal`, spaced
/// `spacing` apart in the first in-plane direction and `0.55 * spacing` apart in
/// the second one so that the in-plane covariance eigenvalues differ.
fn tilted_plane_cloud(normal: [f64; 3], origin: [f64; 3], spacing: f64) -> PointCloud {
    let helper = if normal[2].abs() < 0.9 {
        [0.0, 0.0, 1.0]
    } else {
        [1.0, 0.0, 0.0]
    };
    let u = unit(cross(normal, helper));
    let v = cross(normal, u);
    grid_cloud(7, 5, |i, j| {
        let s = (i - 3.0) * spacing;
        let t = (j - 2.0) * 0.55 * spacing;
        [
            origin[0] + s * u[0] + t * v[0],
            origin[1] + s * u[1] + t * v[1],
            origin[2] + s * u[2] + t * v[2],
        ]
    })
}

#[test]
fn tilted_plane_normal_is_recovered_to_nine_digits() {
    let normal = unit([1.0, 2.0, 3.0]);
    let cloud = tilted_plane_cloud(normal, [0.3, -0.2, 0.9], 0.31);
    assert_all_normals(&cloud, normal, 1e-9);
}

#[test]
fn normals_of_a_plane_through_an_axis_direction_are_recovered() {
    // A plane whose normal is exactly a diagonal of the xy plane.
    let normal = unit([1.0, 1.0, 0.0]);
    let cloud = tilted_plane_cloud(normal, [1.0, 1.0, 1.0], 0.5);
    assert_all_normals(&cloud, normal, 1e-9);
}

#[test]
fn normal_estimation_does_not_depend_on_the_scale_of_the_cloud() {
    let normal = unit([-2.0, 1.0, 0.5]);
    for spacing in [1e-7, 1e-3, 1.0, 1e3, 1e6] {
        let cloud = tilted_plane_cloud(normal, [spacing, -2.0 * spacing, 0.5 * spacing], spacing);
        assert_all_normals(&cloud, normal, 1e-8);
    }
}

#[test]
fn local_neighbourhoods_of_a_plane_agree_with_the_plane_normal() {
    let normal = unit([0.0, 0.3, 1.0]);
    let cloud = tilted_plane_cloud(normal, [0.0; 3], 1.0);
    let normals = NormalEstimation::new(8).estimate(&cloud);
    assert_eq!(normals.len(), cloud.points.len());
    for n in normals {
        assert!(normal_error(n, normal) < 1e-9, "local normal {n:?}");
    }
}

#[test]
fn smallest_eigenvector_of_diagonal_matrices_is_a_coordinate_axis() {
    let cases = [
        (
            [[5.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 2.0]],
            [0.0, 1.0, 0.0],
        ),
        (
            [[1.0, 0.0, 0.0], [0.0, 3.0, 0.0], [0.0, 0.0, 0.5]],
            [0.0, 0.0, 1.0],
        ),
        (
            [[0.1, 0.0, 0.0], [0.0, 3.0, 0.0], [0.0, 0.0, 4.0]],
            [1.0, 0.0, 0.0],
        ),
    ];
    for (matrix, expected) in cases {
        let found = NormalEstimation::smallest_eigenvector_3x3_sym(&matrix);
        assert!(
            normal_error(found, expected) < 1e-12,
            "{matrix:?} gave {found:?}, expected {expected:?}"
        );
    }
}

#[test]
fn smallest_eigenvector_of_a_rotated_matrix_is_the_rotated_axis() {
    // R is a rotation by 0.7 rad about the axis (1, 1, 1) / sqrt(3) (Rodrigues).
    let axis = unit([1.0, 1.0, 1.0]);
    let (sin, cos) = 0.7_f64.sin_cos();
    let k = [
        [0.0, -axis[2], axis[1]],
        [axis[2], 0.0, -axis[0]],
        [-axis[1], axis[0], 0.0],
    ];
    let mut rotation = [[0.0; 3]; 3];
    for (i, row) in rotation.iter_mut().enumerate() {
        for (j, entry) in row.iter_mut().enumerate() {
            let identity = if i == j { 1.0 } else { 0.0 };
            let k_squared: f64 = (0..3).map(|m| k[i][m] * k[m][j]).sum();
            *entry = identity + sin * k[i][j] + (1.0 - cos) * k_squared;
        }
    }
    let eigenvalues = [4.0, 1.0, 0.25];
    let mut matrix = [[0.0; 3]; 3];
    for (i, row) in matrix.iter_mut().enumerate() {
        for (j, entry) in row.iter_mut().enumerate() {
            *entry = (0..3)
                .map(|m| rotation[i][m] * eigenvalues[m] * rotation[j][m])
                .sum();
        }
    }
    let expected = [rotation[0][2], rotation[1][2], rotation[2][2]];
    let found = NormalEstimation::smallest_eigenvector_3x3_sym(&matrix);
    assert!(
        normal_error(found, expected) < 1e-10,
        "gave {found:?}, expected {expected:?}"
    );
}

#[test]
fn smallest_eigenvector_of_a_repeated_eigenvalue_matrix_is_still_an_eigenvector() {
    let matrix = [[2.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 2.0]];
    let found = NormalEstimation::smallest_eigenvector_3x3_sym(&matrix);
    assert!((norm(found) - 1.0).abs() < 1e-12);
    let image = mat_vec(&matrix, found);
    assert!(
        norm([
            image[0] - 2.0 * found[0],
            image[1] - 2.0 * found[1],
            image[2] - 2.0 * found[2]
        ]) < 1e-12
    );
}

#[test]
fn smallest_eigenvector_of_degenerate_matrices_falls_back_to_the_z_axis() {
    let zero = [[0.0; 3]; 3];
    assert_eq!(
        NormalEstimation::smallest_eigenvector_3x3_sym(&zero),
        [0.0, 0.0, 1.0]
    );
    let not_a_number = [[f64::NAN, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    assert_eq!(
        NormalEstimation::smallest_eigenvector_3x3_sym(&not_a_number),
        [0.0, 0.0, 1.0]
    );
    let infinite = [[f64::INFINITY, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    assert_eq!(
        NormalEstimation::smallest_eigenvector_3x3_sym(&infinite),
        [0.0, 0.0, 1.0]
    );
}
