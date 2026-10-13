//! Regression tests for the inertia tensor of a hull.
//!
//! The tensor of a uniform solid is taken about its centre of mass. An earlier
//! version took it about the mean of the vertices and approximated the integral
//! of every tetrahedron by its corner values, which drops the cross terms.

use super::{ConvexHull3DVec, hull_inertia_tensor};
use oxiphysics_core::math::{Mat3, Vec3};

const TOLERANCE: f64 = 1e-12;

fn assert_tensor_close(actual: &[[f64; 3]; 3], expected: &[[f64; 3]; 3], tolerance: f64) {
    for (i, (actual_row, expected_row)) in actual.iter().zip(expected).enumerate() {
        for (j, (a, e)) in actual_row.iter().zip(expected_row).enumerate() {
            assert!(
                (a - e).abs() <= tolerance * (1.0 + e.abs()),
                "I[{i}][{j}] = {a}, expected {e}"
            );
        }
    }
}

fn diagonal(a: f64, b: f64, c: f64) -> [[f64; 3]; 3] {
    [[a, 0.0, 0.0], [0.0, b, 0.0], [0.0, 0.0, c]]
}

/// Eight corners of the axis-aligned box `[lo, hi]`.
fn box_corners(lo: Vec3, hi: Vec3) -> Vec<Vec3> {
    let mut corners = Vec::new();
    for x in [lo.x, hi.x] {
        for y in [lo.y, hi.y] {
            for z in [lo.z, hi.z] {
                corners.push(Vec3::new(x, y, z));
            }
        }
    }
    corners
}

fn hull_of(points: &[Vec3]) -> ConvexHull3DVec {
    ConvexHull3DVec::build(points).expect("the points span a volume")
}

#[test]
fn unit_cube_of_unit_mass_has_inertia_one_sixth() {
    let corners = box_corners(Vec3::zeros(), Vec3::new(1.0, 1.0, 1.0));
    let inertia = hull_inertia_tensor(&hull_of(&corners), 1.0);
    assert_tensor_close(
        &inertia,
        &diagonal(1.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0),
        TOLERANCE,
    );
}

#[test]
fn box_inertia_matches_the_closed_form() {
    // I_xx = m (b^2 + c^2) / 12 and cyclic, for an off-centre box.
    let (a, b, c, m) = (2.0, 3.0, 0.5, 4.5);
    let lo = Vec3::new(1.0, -4.0, 2.5);
    let corners = box_corners(lo, lo + Vec3::new(a, b, c));
    let inertia = hull_inertia_tensor(&hull_of(&corners), m);
    let expected = diagonal(
        m * (b * b + c * c) / 12.0,
        m * (a * a + c * c) / 12.0,
        m * (a * a + b * b) / 12.0,
    );
    assert_tensor_close(&inertia, &expected, TOLERANCE);
}

#[test]
fn square_pyramid_is_taken_about_its_centre_of_mass() {
    // Base side a = 2 in the plane z = 0, apex at height h = 2. The centre of mass
    // lies at h / 4 = 0.5, whereas the vertex mean lies at 0.4.
    // I_zz = m a^2 / 10 and I_xx = I_yy = m (a^2 / 20 + 3 h^2 / 80).
    let points = [
        Vec3::new(-1.0, -1.0, 0.0),
        Vec3::new(1.0, -1.0, 0.0),
        Vec3::new(1.0, 1.0, 0.0),
        Vec3::new(-1.0, 1.0, 0.0),
        Vec3::new(0.0, 0.0, 2.0),
    ];
    let m = 3.0;
    let inertia = hull_inertia_tensor(&hull_of(&points), m);
    let (a, h) = (2.0, 2.0);
    let transverse = m * (a * a / 20.0 + 3.0 * h * h / 80.0);
    let axial = m * a * a / 10.0;
    assert_tensor_close(
        &inertia,
        &diagonal(transverse, transverse, axial),
        TOLERANCE,
    );
}

/// Inertia of the uniform tetrahedron from the closed form of Tonon (2004).
///
/// The second moments about the origin are `V / 10 (sum_{i <= j} v_ip v_jp)` for
/// `x_p^2` and `V / 20 (sum_i 2 v_ip v_iq + sum_{i != j} v_ip v_jq)` for
/// `x_p x_q`; the parallel-axis theorem then moves the tensor to the centre of
/// mass.
fn tetrahedron_inertia_about_centre_of_mass(v: &[Vec3; 4], mass: f64) -> [[f64; 3]; 3] {
    let volume = (v[1] - v[0])
        .dot(&(v[2] - v[0]).cross(&(v[3] - v[0])))
        .abs()
        / 6.0;
    let density = mass / volume;
    let second_moment = |p: usize, q: usize| -> f64 {
        if p == q {
            let pair_sum: f64 = v
                .iter()
                .enumerate()
                .flat_map(|(i, vi)| v.iter().skip(i).map(move |vj| vi[p] * vj[p]))
                .sum();
            density * volume / 10.0 * pair_sum
        } else {
            let cross_sum: f64 = v
                .iter()
                .enumerate()
                .flat_map(|(i, vi)| {
                    v.iter().enumerate().map(move |(j, vj)| {
                        if i == j {
                            2.0 * vi[p] * vi[q]
                        } else {
                            vi[p] * vj[q]
                        }
                    })
                })
                .sum();
            density * volume / 20.0 * cross_sum
        }
    };
    let trace = second_moment(0, 0) + second_moment(1, 1) + second_moment(2, 2);
    let about_origin = Mat3::from_fn(|p, q| {
        let delta = if p == q { trace } else { 0.0 };
        delta - second_moment(p, q)
    });
    let centre = (v[0] + v[1] + v[2] + v[3]) / 4.0;
    let shift = Mat3::identity() * centre.dot(&centre) - centre * centre.transpose();
    let about_centre = about_origin - shift * mass;
    std::array::from_fn(|i| std::array::from_fn(|j| about_centre[(i, j)]))
}

#[test]
fn off_centre_tetrahedron_matches_tonons_closed_form() {
    let vertices = [
        Vec3::new(0.3, -0.2, 0.1),
        Vec3::new(1.7, 0.4, 0.2),
        Vec3::new(0.5, 1.9, -0.3),
        Vec3::new(0.2, 0.6, 1.5),
    ];
    let mass = 2.5;
    let expected = tetrahedron_inertia_about_centre_of_mass(&vertices, mass);
    let inertia = hull_inertia_tensor(&hull_of(&vertices), mass);
    assert_tensor_close(&inertia, &expected, TOLERANCE);
    // The cross terms do not vanish for this skewed tetrahedron.
    assert!(expected[0][1].abs() > 1e-3 && expected[1][2].abs() > 1e-3);
}

#[test]
fn far_from_the_origin_tetrahedron_keeps_its_inertia() {
    let base = [
        Vec3::new(0.3, -0.2, 0.1),
        Vec3::new(1.7, 0.4, 0.2),
        Vec3::new(0.5, 1.9, -0.3),
        Vec3::new(0.2, 0.6, 1.5),
    ];
    let shift = Vec3::new(1.0e3, -2.0e3, 5.0e2);
    let moved = base.map(|v| v + shift);
    let expected = hull_inertia_tensor(&hull_of(&base), 1.0);
    let inertia = hull_inertia_tensor(&hull_of(&moved), 1.0);
    assert_tensor_close(&inertia, &expected, 1e-9);
}

#[test]
fn inertia_scales_linearly_with_mass() {
    let corners = box_corners(Vec3::new(-1.0, -2.0, -0.5), Vec3::new(1.0, 2.0, 0.5));
    let hull = hull_of(&corners);
    let unit = hull_inertia_tensor(&hull, 1.0);
    let heavy = hull_inertia_tensor(&hull, 7.0);
    let expected: [[f64; 3]; 3] =
        std::array::from_fn(|i| std::array::from_fn(|j| 7.0 * unit[i][j]));
    assert_tensor_close(&heavy, &expected, TOLERANCE);
}

#[test]
fn rotated_hull_has_the_rotated_inertia_tensor() {
    let lo = Vec3::new(-1.0, -2.0, -0.5);
    let hi = Vec3::new(1.0, 2.0, 0.5);
    let angle = 0.7_f64;
    let (s, c) = angle.sin_cos();
    let rotation = Mat3::new(c, -s, 0.0, s, c, 0.0, 0.0, 0.0, 1.0);
    let corners = box_corners(lo, hi);
    let rotated: Vec<Vec3> = corners.iter().map(|p| rotation * p).collect();
    let base = hull_inertia_tensor(&hull_of(&corners), 2.0);
    let base_matrix = Mat3::from_fn(|i, j| base[i][j]);
    let expected_matrix = rotation * base_matrix * rotation.transpose();
    let expected: [[f64; 3]; 3] =
        std::array::from_fn(|i| std::array::from_fn(|j| expected_matrix[(i, j)]));
    let inertia = hull_inertia_tensor(&hull_of(&rotated), 2.0);
    assert_tensor_close(&inertia, &expected, 1e-11);
    // The rotated box has a non-zero product of inertia.
    assert!(inertia[0][1].abs() > 1e-3);
}

#[test]
fn inertia_tensor_is_symmetric_and_positive_definite() {
    let points = [
        Vec3::new(0.0, 0.0, 0.0),
        Vec3::new(2.0, 0.1, 0.0),
        Vec3::new(0.3, 1.5, 0.2),
        Vec3::new(0.4, 0.3, 1.1),
        Vec3::new(1.2, 1.1, 1.3),
        Vec3::new(-0.5, 0.8, 0.6),
    ];
    let inertia = hull_inertia_tensor(&hull_of(&points), 1.3);
    for (i, row) in inertia.iter().enumerate() {
        for (j, value) in row.iter().enumerate() {
            assert_eq!(*value, inertia[j][i], "I[{i}][{j}] is not symmetric");
        }
    }
    let matrix = Mat3::from_fn(|i, j| inertia[i][j]);
    let eigenvalues = matrix.symmetric_eigenvalues();
    assert!(
        eigenvalues.iter().all(|&lambda| lambda > 0.0),
        "eigenvalues {eigenvalues:?}"
    );
    // Principal moments satisfy the triangle inequality of a rigid body.
    let (a, b, c) = (eigenvalues[0], eigenvalues[1], eigenvalues[2]);
    assert!(a + b >= c - 1e-12 && a + c >= b - 1e-12 && b + c >= a - 1e-12);
}

#[test]
fn inertia_does_not_depend_on_the_winding_of_the_faces() {
    let corners = box_corners(Vec3::new(0.0, 0.0, 0.0), Vec3::new(1.0, 2.0, 3.0));
    let hull = hull_of(&corners);
    let reference = hull_inertia_tensor(&hull, 1.0);

    let inverted = ConvexHull3DVec {
        vertices: hull.vertices.clone(),
        faces: hull.faces.iter().map(|f| [f[0], f[2], f[1]]).collect(),
    };
    assert_tensor_close(&hull_inertia_tensor(&inverted, 1.0), &reference, TOLERANCE);

    let mixed = ConvexHull3DVec {
        vertices: hull.vertices.clone(),
        faces: hull
            .faces
            .iter()
            .enumerate()
            .map(|(k, f)| if k % 2 == 0 { *f } else { [f[0], f[2], f[1]] })
            .collect(),
    };
    assert_tensor_close(&hull_inertia_tensor(&mixed, 1.0), &reference, TOLERANCE);
}

#[test]
fn hulls_without_volume_have_no_inertia() {
    let zero = [[0.0; 3]; 3];
    let empty = ConvexHull3DVec {
        vertices: Vec::new(),
        faces: Vec::new(),
    };
    assert_eq!(hull_inertia_tensor(&empty, 1.0), zero);

    let flat = ConvexHull3DVec {
        vertices: vec![
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
            Vec3::new(1.0, 1.0, 0.0),
        ],
        faces: vec![[0, 1, 2], [1, 3, 2], [0, 2, 1], [1, 2, 3]],
    };
    assert_eq!(hull_inertia_tensor(&flat, 1.0), zero);
}

#[test]
fn small_hulls_are_not_mistaken_for_degenerate_ones() {
    // A cube of 1 micrometre side has a volume of 1e-18, far below the former
    // absolute cut-off of 1e-12.
    let side = 1.0e-6;
    let corners = box_corners(Vec3::zeros(), Vec3::new(side, side, side));
    let inertia = hull_inertia_tensor(&hull_of(&corners), 1.0);
    let expected = side * side / 6.0;
    assert_tensor_close(&inertia, &diagonal(expected, expected, expected), 1e-9);
}
