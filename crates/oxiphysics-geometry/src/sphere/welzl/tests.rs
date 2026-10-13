//! Regression tests for the circumspheres and the minimum enclosing sphere.

use super::*;
use crate::sphere::{minimum_enclosing_sphere, random_points_on_sphere};

/// Deterministic uniform samples in `[0, 1)` from a 64-bit LCG.
struct Lcg(u64);

impl Lcg {
    fn next_unit(&mut self) -> f64 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        (self.0 >> 11) as f64 / (1u64 << 53) as f64
    }

    fn point_in_box(&mut self, half_extent: f64) -> [f64; 3] {
        [
            (self.next_unit() * 2.0 - 1.0) * half_extent,
            (self.next_unit() * 2.0 - 1.0) * half_extent,
            (self.next_unit() * 2.0 - 1.0) * half_extent,
        ]
    }
}

fn assert_vec_close(actual: [f64; 3], expected: [f64; 3], tolerance: f64) {
    let distance = sphere_dist3(actual, expected);
    assert!(
        distance <= tolerance,
        "{actual:?} differs from {expected:?} by {distance} (tolerance {tolerance})"
    );
}

fn assert_encloses(points: &[[f64; 3]], center: [f64; 3], radius: f64, tolerance: f64) {
    for p in points {
        let distance = sphere_dist3(*p, center);
        assert!(
            distance <= radius + tolerance,
            "point {p:?} is {distance} from the centre {center:?}, radius {radius}"
        );
    }
}

/// Smallest sphere through every point of `support` (one to four affinely
/// independent points), from the normal equations `2 E E^T lambda = diag(E E^T)`
/// of the edge vectors `E`. It shares no formula with the implementation under
/// test. Returns `None` for affinely dependent supports.
fn reference_support_sphere(support: &[[f64; 3]]) -> Option<([f64; 3], f64)> {
    let (first, rest) = support.split_first()?;
    let edges: Vec<[f64; 3]> = rest.iter().map(|q| sub3sph(*q, *first)).collect();
    let k = edges.len();
    let mut matrix: Vec<Vec<f64>> = edges
        .iter()
        .map(|ei| {
            let mut row: Vec<f64> = edges.iter().map(|ej| 2.0 * dot3sph(*ei, *ej)).collect();
            row.push(dot3sph(*ei, *ei));
            row
        })
        .collect();
    for col in 0..k {
        let pivot_row =
            (col..k).max_by(|&a, &b| matrix[a][col].abs().total_cmp(&matrix[b][col].abs()))?;
        matrix.swap(col, pivot_row);
        let pivot = matrix[col][col];
        if pivot.abs() < 1e-12 {
            return None;
        }
        let pivot_values: Vec<f64> = matrix[col].iter().map(|value| value / pivot).collect();
        for (row, values) in matrix.iter_mut().enumerate() {
            if row == col {
                values.clone_from(&pivot_values);
            } else {
                let factor = values[col];
                for (value, reduced) in values.iter_mut().zip(&pivot_values) {
                    *value -= factor * reduced;
                }
            }
        }
    }
    let mut center = *first;
    for (row, edge) in matrix.iter().zip(&edges) {
        let lambda = row[k];
        for (c, e) in center.iter_mut().zip(edge) {
            *c += lambda * e;
        }
    }
    Some((center, sphere_dist3(*first, center)))
}

/// Exact minimum enclosing sphere of a small point set: the smallest sphere
/// among all support-set spheres of one to four points that encloses every point.
fn brute_force_minimum_sphere(points: &[[f64; 3]]) -> ([f64; 3], f64) {
    let mut best: Option<([f64; 3], f64)> = None;
    for subset in 1_u32..(1 << points.len()) {
        if subset.count_ones() > 4 {
            continue;
        }
        let support: Vec<[f64; 3]> = points
            .iter()
            .enumerate()
            .filter(|(index, _)| subset & (1 << index) != 0)
            .map(|(_, p)| *p)
            .collect();
        let Some((center, radius)) = reference_support_sphere(&support) else {
            continue;
        };
        let encloses_all = points
            .iter()
            .all(|p| sphere_dist3(*p, center) <= radius * (1.0 + 1e-9) + 1e-12);
        if encloses_all && best.is_none_or(|(_, current)| radius < current) {
            best = Some((center, radius));
        }
    }
    best.unwrap_or(([0.0; 3], 0.0))
}

// ── circumsphere_3 ───────────────────────────────────────────────────────────

#[test]
fn circumsphere_3_right_triangle_centre_is_the_hypotenuse_midpoint() {
    let (center, radius) = circumsphere_3([0.0; 3], [2.0, 0.0, 0.0], [0.0, 2.0, 0.0]);
    assert_vec_close(center, [1.0, 1.0, 0.0], 1e-12);
    assert!((radius - 2.0_f64.sqrt()).abs() < 1e-12, "radius {radius}");
}

#[test]
fn circumsphere_3_tilted_triangle_is_equidistant_and_coplanar() {
    let a = [1.0, 2.0, 3.0];
    let b = [4.0, 2.0, -1.0];
    let c = [0.0, 5.0, 2.0];
    let (center, radius) = circumsphere_3(a, b, c);
    for p in [a, b, c] {
        assert!(
            (sphere_dist3(p, center) - radius).abs() < 1e-12,
            "vertex {p:?} is not on the sphere"
        );
    }
    let normal = cross3sph(sub3sph(b, a), sub3sph(c, a));
    assert!(
        dot3sph(sub3sph(center, a), normal).abs() < 1e-10,
        "the centre must lie in the plane of the triangle"
    );
}

#[test]
fn circumsphere_3_is_independent_of_the_vertex_order() {
    let a = [1.0, 2.0, 3.0];
    let b = [4.0, 2.0, -1.0];
    let c = [0.0, 5.0, 2.0];
    let (reference, reference_radius) = circumsphere_3(a, b, c);
    for (p, q, r) in [(a, c, b), (b, a, c), (b, c, a), (c, a, b), (c, b, a)] {
        let (center, radius) = circumsphere_3(p, q, r);
        assert_vec_close(center, reference, 1e-12);
        assert!((radius - reference_radius).abs() < 1e-12);
    }
}

#[test]
fn circumsphere_3_obtuse_triangle_matches_the_analytic_circumcentre() {
    // Centre (5, y0, 0) with 25 + y0^2 = (0.5 - y0)^2, i.e. y0 = -24.75.
    let (center, radius) = circumsphere_3([0.0; 3], [10.0, 0.0, 0.0], [5.0, 0.5, 0.0]);
    assert_vec_close(center, [5.0, -24.75, 0.0], 1e-10);
    assert!((radius - 25.25).abs() < 1e-10, "radius {radius}");
}

#[test]
fn circumsphere_3_thin_needle_triangle_stays_accurate() {
    // Apex at the origin, base of half-width h far away: by symmetry the centre
    // is (cx, 0, 0) with cx^2 = (10 - cx)^2 + h^2.
    let h = 1e-6;
    let (center, radius) = circumsphere_3([0.0; 3], [10.0, h, 0.0], [10.0, -h, 0.0]);
    let cx = (100.0 + h * h) / 20.0;
    assert_vec_close(center, [cx, 0.0, 0.0], 1e-9);
    assert!((radius - cx).abs() < 1e-9, "radius {radius}");
}

#[test]
fn circumsphere_3_collinear_points_fall_back_to_the_enclosing_diametral_sphere() {
    let a = [0.0; 3];
    let b = [1.0, 0.0, 0.0];
    let c = [3.0, 0.0, 0.0];
    for (p, q, r) in [(a, b, c), (b, c, a), (c, a, b), (c, b, a)] {
        let (center, radius) = circumsphere_3(p, q, r);
        assert_vec_close(center, [1.5, 0.0, 0.0], 1e-12);
        assert!((radius - 1.5).abs() < 1e-12, "radius {radius}");
    }
}

#[test]
fn circumsphere_3_coincident_points_collapse_to_a_point_sphere() {
    let p = [1.0, 1.0, 1.0];
    let (center, radius) = circumsphere_3(p, p, p);
    assert_vec_close(center, p, 1e-15);
    assert_eq!(radius, 0.0);

    let (center, radius) = circumsphere_3([0.0; 3], [0.0; 3], [2.0, 0.0, 0.0]);
    assert_vec_close(center, [1.0, 0.0, 0.0], 1e-15);
    assert!((radius - 1.0).abs() < 1e-15);
}

// ── circumsphere_4 ───────────────────────────────────────────────────────────

#[test]
fn circumsphere_4_regular_tetrahedron_recovers_the_exact_sphere() {
    let center = [1.0, -2.0, 3.0];
    let radius = 2.0;
    let s = radius / 3.0_f64.sqrt();
    let corners = [
        [1.0, 1.0, 1.0],
        [1.0, -1.0, -1.0],
        [-1.0, 1.0, -1.0],
        [-1.0, -1.0, 1.0],
    ]
    .map(|v| {
        [
            center[0] + s * v[0],
            center[1] + s * v[1],
            center[2] + s * v[2],
        ]
    });
    let (c, r) = circumsphere_4(corners[0], corners[1], corners[2], corners[3]);
    assert_vec_close(c, center, 1e-12);
    assert!((r - radius).abs() < 1e-12, "radius {r}");
}

#[test]
fn circumsphere_4_unit_corner_tetrahedron_has_the_cube_centre() {
    let (center, radius) =
        circumsphere_4([0.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]);
    assert_vec_close(center, [0.5, 0.5, 0.5], 1e-12);
    assert!(
        (radius - 3.0_f64.sqrt() / 2.0).abs() < 1e-12,
        "radius {radius}"
    );
}

#[test]
fn circumsphere_4_skewed_tetrahedron_is_equidistant_for_every_vertex_order() {
    let corners = [
        [0.0, 0.0, 0.0],
        [3.0, 0.0, 0.0],
        [0.5, 2.0, 0.0],
        [1.0, 1.0, 4.0],
    ];
    let (reference, reference_radius) =
        circumsphere_4(corners[0], corners[1], corners[2], corners[3]);
    for corner in &corners {
        assert!(
            (sphere_dist3(*corner, reference) - reference_radius).abs() < 1e-12,
            "corner {corner:?} is not on the sphere"
        );
    }
    let mut checked = 0;
    for (i, a) in corners.iter().enumerate() {
        for (j, b) in corners.iter().enumerate().filter(|(j, _)| *j != i) {
            for (k, c) in corners
                .iter()
                .enumerate()
                .filter(|(k, _)| *k != i && *k != j)
            {
                for (_, d) in corners
                    .iter()
                    .enumerate()
                    .filter(|(l, _)| *l != i && *l != j && *l != k)
                {
                    let (center, radius) = circumsphere_4(*a, *b, *c, *d);
                    assert_vec_close(center, reference, 1e-12);
                    assert!((radius - reference_radius).abs() < 1e-12);
                    checked += 1;
                }
            }
        }
    }
    assert_eq!(checked, 24, "every vertex order must be visited");
}

#[test]
fn circumsphere_4_needle_tetrahedron_uses_a_well_conditioned_corner() {
    // Apex at the origin, tiny base far away. Seen from the apex the three edge
    // vectors are nearly parallel, so Cramer's rule is only accurate at a base
    // corner. By symmetry the centre is (cx, 0, 0) with cx^2 = (10 - cx)^2 + h^2.
    let h = 1e-5;
    let (center, radius) =
        circumsphere_4([0.0; 3], [10.0, h, 0.0], [10.0, -h, 0.0], [10.0, 0.0, h]);
    let cx = (100.0 + h * h) / 20.0;
    assert_vec_close(center, [cx, 0.0, 0.0], 1e-9);
    assert!((radius - cx).abs() < 1e-9, "radius {radius}");
}

#[test]
fn circumsphere_4_concyclic_coplanar_points_use_the_circle_centre() {
    let (center, radius) =
        circumsphere_4([0.0; 3], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0]);
    assert_vec_close(center, [0.5, 0.5, 0.0], 1e-12);
    assert!(
        (radius - 0.5 * 2.0_f64.sqrt()).abs() < 1e-12,
        "radius {radius}"
    );
}

#[test]
fn circumsphere_4_non_concyclic_coplanar_points_return_an_enclosing_sphere() {
    // The hypotenuse of the right triangle (0,0), (4,0), (0,3) is a diameter and
    // the fourth point lies inside the resulting sphere.
    let corners = [
        [0.0, 0.0, 0.0],
        [4.0, 0.0, 0.0],
        [0.0, 3.0, 0.0],
        [1.0, 1.0, 0.0],
    ];
    let (center, radius) = circumsphere_4(corners[0], corners[1], corners[2], corners[3]);
    assert_vec_close(center, [2.0, 1.5, 0.0], 1e-12);
    assert!((radius - 2.5).abs() < 1e-12, "radius {radius}");
    assert_encloses(&corners, center, radius, 1e-12);
}

// ── minimum_enclosing_sphere ─────────────────────────────────────────────────

#[test]
fn minimum_enclosing_sphere_of_a_right_triangle_is_its_hypotenuse_sphere() {
    let points = [[0.0; 3], [2.0, 0.0, 0.0], [0.0, 2.0, 0.0]];
    let (center, radius) = minimum_enclosing_sphere(&points);
    assert_vec_close(center, [1.0, 1.0, 0.0], 1e-12);
    assert!((radius - 2.0_f64.sqrt()).abs() < 1e-12, "radius {radius}");
}

#[test]
fn minimum_enclosing_sphere_of_an_obtuse_triangle_uses_the_longest_edge() {
    let points = [[0.0; 3], [10.0, 0.0, 0.0], [5.0, 0.5, 0.0]];
    let (center, radius) = minimum_enclosing_sphere(&points);
    assert_vec_close(center, [5.0, 0.0, 0.0], 1e-12);
    assert!((radius - 5.0).abs() < 1e-12, "radius {radius}");
}

#[test]
fn minimum_enclosing_sphere_of_nothing_and_of_one_point() {
    assert_eq!(minimum_enclosing_sphere(&[]), ([0.0; 3], 0.0));
    let (center, radius) = minimum_enclosing_sphere(&[[1.0, -2.0, 3.0]]);
    assert_eq!(center, [1.0, -2.0, 3.0]);
    assert_eq!(radius, 0.0);
}

#[test]
fn minimum_enclosing_sphere_of_identical_points_is_a_point() {
    let points = vec![[1.0, 2.0, 3.0]; 1000];
    let (center, radius) = minimum_enclosing_sphere(&points);
    assert_vec_close(center, [1.0, 2.0, 3.0], 1e-15);
    assert_eq!(radius, 0.0);
}

#[test]
fn minimum_enclosing_sphere_of_collinear_points_spans_the_extremes() {
    let points: Vec<[f64; 3]> = (0..100).map(|i| [0.1 * i as f64, 0.0, 0.0]).collect();
    let (center, radius) = minimum_enclosing_sphere(&points);
    assert_vec_close(center, [4.95, 0.0, 0.0], 1e-9);
    assert!((radius - 4.95).abs() < 1e-9, "radius {radius}");
}

#[test]
fn minimum_enclosing_sphere_of_a_coplanar_grid_is_centred_on_the_grid() {
    let mut points = Vec::new();
    for i in 0..5 {
        for j in 0..5 {
            points.push([i as f64, j as f64, 1.0]);
        }
    }
    let (center, radius) = minimum_enclosing_sphere(&points);
    assert_vec_close(center, [2.0, 2.0, 1.0], 1e-9);
    assert!((radius - 8.0_f64.sqrt()).abs() < 1e-9, "radius {radius}");
}

#[test]
fn minimum_enclosing_sphere_1000_seeded_points_are_all_enclosed_tightly() {
    let mut rng = Lcg(2024);
    let points: Vec<[f64; 3]> = (0..1000).map(|_| rng.point_in_box(5.0)).collect();
    let (center, radius) = minimum_enclosing_sphere(&points);
    assert_encloses(&points, center, radius, 1e-9);

    // The sphere is tight: at least two input points sit on its surface.
    let on_surface = points
        .iter()
        .filter(|p| (sphere_dist3(**p, center) - radius).abs() < 1e-9)
        .count();
    assert!(on_surface >= 2, "only {on_surface} points on the surface");

    // The radius is at least half the diameter and, by Jung's theorem, at most
    // sqrt(3/8) times the diameter.
    let mut diameter = 0.0_f64;
    for (i, p) in points.iter().enumerate() {
        for q in &points[i + 1..] {
            diameter = diameter.max(sphere_dist3(*p, *q));
        }
    }
    assert!(
        radius >= 0.5 * diameter - 1e-9,
        "radius {radius}, diameter {diameter}"
    );
    assert!(
        radius <= (3.0_f64 / 8.0).sqrt() * diameter + 1e-9,
        "radius {radius}, diameter {diameter}"
    );
}

#[test]
fn minimum_enclosing_sphere_matches_brute_force_on_small_random_sets() {
    for seed in 0..300_u64 {
        let mut rng = Lcg(seed.wrapping_mul(7919).wrapping_add(13));
        let count = 2 + (seed % 8) as usize;
        let points: Vec<[f64; 3]> = (0..count).map(|_| rng.point_in_box(3.0)).collect();
        let (center, radius) = minimum_enclosing_sphere(&points);
        let (reference_center, reference_radius) = brute_force_minimum_sphere(&points);
        assert!(
            (radius - reference_radius).abs() < 1e-9 * reference_radius.max(1.0),
            "seed {seed}: radius {radius} vs brute force {reference_radius}"
        );
        assert_vec_close(center, reference_center, 1e-8);
    }
}

#[test]
fn minimum_enclosing_sphere_does_not_depend_on_the_input_order() {
    let mut rng = Lcg(5);
    let points: Vec<[f64; 3]> = (0..400).map(|_| rng.point_in_box(2.0)).collect();
    let mut reversed = points.clone();
    reversed.reverse();
    let (center, radius) = minimum_enclosing_sphere(&points);
    let (reversed_center, reversed_radius) = minimum_enclosing_sphere(&reversed);
    assert_vec_close(center, reversed_center, 1e-9);
    assert!((radius - reversed_radius).abs() < 1e-9);
}

#[test]
fn minimum_enclosing_sphere_of_cospherical_points_is_that_sphere() {
    let points = random_points_on_sphere(3.0, 5000, 7);
    let (center, radius) = minimum_enclosing_sphere(&points);
    assert_vec_close(center, [0.0; 3], 1e-9);
    assert!((radius - 3.0).abs() < 1e-9, "radius {radius}");
}

#[test]
fn minimum_enclosing_sphere_handles_200_000_points_without_deep_recursion() {
    let mut rng = Lcg(99);
    let points: Vec<[f64; 3]> = (0..200_000).map(|_| rng.point_in_box(10.0)).collect();
    let (center, radius) = minimum_enclosing_sphere(&points);
    assert_encloses(&points, center, radius, 1e-9);
    // The cube [-10, 10]^3 has circumradius 10 sqrt(3); dense random points get
    // close to it but never beyond.
    assert!(radius <= 10.0 * 3.0_f64.sqrt() + 1e-9, "radius {radius}");
    assert!(radius > 17.0, "radius {radius}");

    let spherical = random_points_on_sphere(1.0, 200_000, 5);
    let (center, radius) = minimum_enclosing_sphere(&spherical);
    assert_vec_close(center, [0.0; 3], 1e-9);
    assert!((radius - 1.0).abs() < 1e-9, "radius {radius}");
}

// ── helpers ──────────────────────────────────────────────────────────────────

#[test]
fn shuffle_points_is_a_deterministic_permutation() {
    let original: Vec<[f64; 3]> = (0..100).map(|i| [i as f64, 0.0, 0.0]).collect();
    let mut first = original.clone();
    let mut second = original.clone();
    shuffle_points(&mut first, 17);
    shuffle_points(&mut second, 17);
    assert_eq!(first, second);
    assert_ne!(first, original);

    let mut sorted = first.clone();
    sorted.sort_by(|a, b| a[0].total_cmp(&b[0]));
    assert_eq!(sorted, original);
}

#[test]
fn enclosing_sphere_of_candidates_handles_empty_and_single_inputs() {
    assert_eq!(enclosing_sphere_of_candidates(&[]), ([0.0; 3], 0.0));
    assert_eq!(
        enclosing_sphere_of_candidates(&[[2.0, 3.0, 4.0]]),
        ([2.0, 3.0, 4.0], 0.0)
    );
}
