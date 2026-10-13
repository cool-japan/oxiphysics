//! Robustness regression tests for the QuickHull builder.
//!
//! The hull of a point cloud must be a closed, consistently wound surface that
//! encloses every input point, whatever the order, the scale or the symmetry of
//! the input.  An earlier version oriented new faces against an arbitrary
//! vertex of the surviving hull; on symmetric inputs such as a sorted
//! icosphere this produced open surfaces with hundreds of bogus faces, and the
//! seed selection rejected valid clouds depending on which side of the first
//! triangle the fourth point lay.

use super::{ConvexHull3DVec, is_point_set_convex};
use oxiphysics_core::math::{Vec3, quat_from_axis_angle};
use std::collections::HashMap;

/// Small deterministic linear congruential generator for reproducible clouds.
struct Lcg(u64);

impl Lcg {
    fn new(seed: u64) -> Self {
        Self(seed)
    }

    /// Uniform sample in `[0, 1)`.
    fn unit(&mut self) -> f64 {
        self.0 = self
            .0
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        (self.0 >> 33) as f64 / (1u64 << 31) as f64
    }

    /// Uniform sample in `[lo, hi)`.
    fn range(&mut self, lo: f64, hi: f64) -> f64 {
        lo + (hi - lo) * self.unit()
    }

    /// Uniform index in `0..n` (`n > 0`).
    fn below(&mut self, n: usize) -> usize {
        ((self.unit() * n as f64) as usize).min(n - 1)
    }
}

fn shuffled(points: &[Vec3], seed: u64) -> Vec<Vec3> {
    let mut rng = Lcg::new(seed);
    let mut out = points.to_vec();
    for i in (1..out.len()).rev() {
        let j = rng.below(i + 1);
        out.swap(i, j);
    }
    out
}

fn sorted(points: &[Vec3]) -> Vec<Vec3> {
    let mut out = points.to_vec();
    out.sort_by(|a, b| {
        a.x.total_cmp(&b.x)
            .then(a.y.total_cmp(&b.y))
            .then(a.z.total_cmp(&b.z))
    });
    out
}

/// Icosahedron split once and projected on the sphere: 42 points in convex
/// position with many exactly symmetric configurations. Returns the vertices
/// and the 80 outward-facing triangles.
fn icosphere(radius: f64) -> (Vec<Vec3>, Vec<[usize; 3]>) {
    let phi = (1.0 + 5.0_f64.sqrt()) / 2.0;
    let raw = [
        [-1.0, phi, 0.0],
        [1.0, phi, 0.0],
        [-1.0, -phi, 0.0],
        [1.0, -phi, 0.0],
        [0.0, -1.0, phi],
        [0.0, 1.0, phi],
        [0.0, -1.0, -phi],
        [0.0, 1.0, -phi],
        [phi, 0.0, -1.0],
        [phi, 0.0, 1.0],
        [-phi, 0.0, -1.0],
        [-phi, 0.0, 1.0],
    ];
    let mut vertices: Vec<Vec3> = raw
        .iter()
        .map(|v| Vec3::new(v[0], v[1], v[2]).normalize() * radius)
        .collect();
    let base: [[usize; 3]; 20] = [
        [0, 11, 5],
        [0, 5, 1],
        [0, 1, 7],
        [0, 7, 10],
        [0, 10, 11],
        [1, 5, 9],
        [5, 11, 4],
        [11, 10, 2],
        [10, 7, 6],
        [7, 1, 8],
        [3, 9, 4],
        [3, 4, 2],
        [3, 2, 6],
        [3, 6, 8],
        [3, 8, 9],
        [4, 9, 5],
        [2, 4, 11],
        [6, 2, 10],
        [8, 6, 7],
        [9, 8, 1],
    ];
    let mut midpoints: HashMap<(usize, usize), usize> = HashMap::new();
    let mut midpoint = |vertices: &mut Vec<Vec3>, a: usize, b: usize| -> usize {
        let key = (a.min(b), a.max(b));
        if let Some(&index) = midpoints.get(&key) {
            return index;
        }
        let m = ((vertices[a] + vertices[b]) * 0.5).normalize() * radius;
        vertices.push(m);
        midpoints.insert(key, vertices.len() - 1);
        vertices.len() - 1
    };
    let mut triangles = Vec::new();
    for tri in &base {
        let ab = midpoint(&mut vertices, tri[0], tri[1]);
        let bc = midpoint(&mut vertices, tri[1], tri[2]);
        let ca = midpoint(&mut vertices, tri[2], tri[0]);
        triangles.push([tri[0], ab, ca]);
        triangles.push([tri[1], bc, ab]);
        triangles.push([tri[2], ca, bc]);
        triangles.push([ab, bc, ca]);
    }
    (vertices, triangles)
}

/// Volume enclosed by a closed outward-facing triangle mesh.
fn mesh_volume(vertices: &[Vec3], triangles: &[[usize; 3]]) -> f64 {
    triangles
        .iter()
        .map(|t| vertices[t[0]].dot(&vertices[t[1]].cross(&vertices[t[2]])) / 6.0)
        .sum()
}

fn diagonal(points: &[Vec3]) -> f64 {
    let mut lo = points[0];
    let mut hi = points[0];
    for p in points {
        lo = lo.inf(p);
        hi = hi.sup(p);
    }
    (hi - lo).norm()
}

/// Asserts that `hull` is a closed, outward-wound triangulated sphere that
/// encloses all of `points`.
fn assert_valid_hull(points: &[Vec3], hull: &ConvexHull3DVec) {
    // Every directed edge occurs once and is matched by its reverse.
    let mut directed: HashMap<(usize, usize), usize> = HashMap::new();
    for f in &hull.faces {
        for k in 0..3 {
            *directed.entry((f[k], f[(k + 1) % 3])).or_insert(0) += 1;
        }
    }
    for (&(a, b), &count) in &directed {
        assert_eq!(count, 1, "directed edge ({a}, {b}) is used {count} times");
        assert_eq!(
            directed.get(&(b, a)),
            Some(&1),
            "edge ({a}, {b}) has no reverse partner"
        );
    }
    assert_eq!(hull.euler_characteristic(), 2, "not a topological sphere");

    // All input points lie behind every face plane; the faces look outward.
    let tolerance = 1e-9 * diagonal(points);
    let centroid = hull.vertices.iter().sum::<Vec3>() / hull.vertices.len() as f64;
    for f in &hull.faces {
        let (a, b, c) = (
            hull.vertices[f[0]],
            hull.vertices[f[1]],
            hull.vertices[f[2]],
        );
        let normal = (b - a).cross(&(c - a));
        let length = normal.norm();
        assert!(length > 0.0, "zero-area face {f:?}");
        let normal = normal / length;
        assert!(normal.dot(&(centroid - a)) < 0.0, "face {f:?} looks inward");
        for p in points {
            assert!(
                normal.dot(&(p - a)) <= tolerance,
                "point {p:?} lies outside face {f:?}"
            );
        }
    }
    assert!(hull.volume() > 0.0);
}

fn same_vertex_set(a: &ConvexHull3DVec, b: &ConvexHull3DVec) -> bool {
    let (a, b) = (sorted(&a.vertices), sorted(&b.vertices));
    a.len() == b.len() && a.iter().zip(&b).all(|(p, q)| (p - q).norm() < 1e-12)
}

fn random_cloud(rng: &mut Lcg, kind: usize, n: usize) -> Vec<Vec3> {
    (0..n)
        .map(|_| match kind {
            0 => Vec3::new(
                rng.range(-1.0, 1.0),
                rng.range(-1.0, 1.0),
                rng.range(-1.0, 1.0),
            ),
            1 => {
                let v = Vec3::new(
                    rng.range(-1.0, 1.0),
                    rng.range(-1.0, 1.0),
                    rng.range(-1.0, 1.0),
                );
                // Points on a sphere: all of them are hull vertices.
                v.try_normalize(1e-9).unwrap_or_else(Vec3::x) * 2.0
            }
            _ => Vec3::new(
                rng.range(-3.0, 3.0),
                rng.range(-1.0, 1.0),
                rng.range(0.0, 1e-2),
            ),
        })
        .collect()
}

#[test]
fn sorted_icosphere_builds_the_exact_hull() {
    let (vertices, triangles) = icosphere(2.0);
    let expected_volume = mesh_volume(&vertices, &triangles);
    assert!(expected_volume > 29.0 && expected_volume < 30.0);
    for points in [
        vertices.clone(),
        sorted(&vertices),
        sorted(&vertices).into_iter().rev().collect(),
    ] {
        let hull = ConvexHull3DVec::build(&points).expect("icosphere is a 3D point set");
        assert_valid_hull(&points, &hull);
        assert_eq!(hull.n_vertices(), 42);
        assert_eq!(hull.n_faces(), 80);
        assert!((hull.volume() - expected_volume).abs() < 1e-9 * expected_volume);
        assert!(is_point_set_convex(&points));
    }
}

#[test]
fn hull_does_not_depend_on_the_point_order() {
    let (icosphere_points, _) = icosphere(1.5);
    let mut rng = Lcg::new(12);
    let ball = random_cloud(&mut rng, 0, 80);
    for points in [icosphere_points, ball] {
        let reference = ConvexHull3DVec::build(&points).expect("reference hull");
        let mut orders = vec![sorted(&points), sorted(&points).into_iter().rev().collect()];
        orders.extend((0..6).map(|seed| shuffled(&points, seed)));
        for order in orders {
            let hull = ConvexHull3DVec::build(&order).expect("hull of a permutation");
            assert_valid_hull(&order, &hull);
            assert!(same_vertex_set(&reference, &hull));
            assert!((hull.volume() - reference.volume()).abs() < 1e-9 * reference.volume());
        }
    }
}

#[test]
fn pyramid_builds_with_the_apex_on_either_side_of_the_base() {
    let base = [
        Vec3::new(0.0, 0.0, 0.0),
        Vec3::new(1.0, 0.0, 0.0),
        Vec3::new(1.0, 1.0, 0.0),
        Vec3::new(0.0, 1.0, 0.0),
    ];
    for height in [1.0, -1.0, 0.25, -3.0] {
        let apex = Vec3::new(0.5, 0.5, height);
        let mut points = base.to_vec();
        points.push(apex);
        let mut orders = vec![points.clone(), sorted(&points)];
        orders.push(std::iter::once(apex).chain(base).collect());
        orders.extend((0..4).map(|seed| shuffled(&points, seed)));
        for order in orders {
            let hull = ConvexHull3DVec::build(&order)
                .unwrap_or_else(|| panic!("pyramid of height {height} must build"));
            assert_valid_hull(&order, &hull);
            assert_eq!(hull.n_vertices(), 5);
            assert_eq!(hull.n_faces(), 6);
            let expected = height.abs() / 3.0;
            assert!((hull.volume() - expected).abs() < 1e-12, "height {height}");
        }
    }
}

#[test]
fn construction_is_independent_of_the_scale_of_the_input() {
    let mut rng = Lcg::new(77);
    let unit_cloud = random_cloud(&mut rng, 0, 60);
    let reference = ConvexHull3DVec::build(&unit_cloud).expect("unit-scale hull");
    for scale in [1e-6, 1e-4, 1e-2, 1.0, 1e3, 1e6] {
        let points: Vec<Vec3> = unit_cloud.iter().map(|p| p * scale).collect();
        let hull = ConvexHull3DVec::build(&points)
            .unwrap_or_else(|| panic!("scale {scale} must build a hull"));
        assert_valid_hull(&points, &hull);
        assert_eq!(hull.n_vertices(), reference.n_vertices(), "scale {scale}");
        let expected = reference.volume() * scale.powi(3);
        assert!(
            (hull.volume() - expected).abs() < 1e-9 * expected,
            "scale {scale}: volume {} vs {expected}",
            hull.volume()
        );
    }
}

#[test]
fn degenerate_and_invalid_input_has_no_hull() {
    let p = |x: f64, y: f64, z: f64| Vec3::new(x, y, z);
    // Too few points.
    assert!(ConvexHull3DVec::build(&[]).is_none());
    assert!(
        ConvexHull3DVec::build(&[p(0.0, 0.0, 0.0), p(1.0, 0.0, 0.0), p(0.0, 1.0, 0.0)]).is_none()
    );
    // All points identical.
    assert!(ConvexHull3DVec::build(&[p(1.0, 2.0, 3.0); 6]).is_none());
    // Collinear points.
    let line: Vec<Vec3> = (0..8)
        .map(|i| p(i as f64, 2.0 * i as f64, -(i as f64)))
        .collect();
    assert!(ConvexHull3DVec::build(&line).is_none());
    // Coplanar points, in a plane that is not axis aligned.
    let mut rng = Lcg::new(5);
    let plane: Vec<Vec3> = (0..20)
        .map(|_| {
            let (u, v) = (rng.range(-1.0, 1.0), rng.range(-1.0, 1.0));
            p(u + v, u - v, 0.5 * u + 2.0 * v)
        })
        .collect();
    assert!(ConvexHull3DVec::build(&plane).is_none());
    // Non-finite coordinates are rejected instead of panicking.
    let cube = [
        p(0.0, 0.0, 0.0),
        p(1.0, 0.0, 0.0),
        p(0.0, 1.0, 0.0),
        p(0.0, 0.0, 1.0),
    ];
    for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        let mut points = cube.to_vec();
        points.push(p(0.2, bad, 0.2));
        assert!(ConvexHull3DVec::build(&points).is_none(), "{bad}");
    }
    // The finite subset itself is fine.
    assert!(ConvexHull3DVec::build(&cube).is_some());
}

#[test]
fn lattice_with_many_coplanar_and_collinear_points_is_enclosed_exactly() {
    let mut lattice = Vec::new();
    for x in 0..4 {
        for y in 0..4 {
            for z in 0..4 {
                lattice.push(Vec3::new(f64::from(x), f64::from(y), f64::from(z)));
            }
        }
    }
    let mut orders = vec![lattice.clone(), sorted(&lattice)];
    orders.extend((0..4).map(|seed| shuffled(&lattice, seed)));
    for order in orders {
        let hull = ConvexHull3DVec::build(&order).expect("lattice hull");
        assert_valid_hull(&order, &hull);
        assert!((hull.volume() - 27.0).abs() < 1e-9);
        assert!(order.iter().all(|p| hull.contains_point(*p)));
        for corner in [
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(3.0, 0.0, 0.0),
            Vec3::new(3.0, 3.0, 3.0),
            Vec3::new(0.0, 3.0, 3.0),
        ] {
            assert!(hull.vertices.iter().any(|v| (v - corner).norm() < 1e-12));
        }
    }
}

#[test]
fn random_clouds_of_every_flavour_are_enclosed() {
    let mut rng = Lcg::new(2024);
    for kind in 0..3 {
        for n in [5, 12, 40, 150] {
            for _ in 0..4 {
                let points = random_cloud(&mut rng, kind, n);
                let hull = ConvexHull3DVec::build(&points)
                    .unwrap_or_else(|| panic!("kind {kind}, n {n} must build"));
                assert_valid_hull(&points, &hull);
                if kind == 1 {
                    // Points on a sphere are all extreme.
                    assert_eq!(hull.n_vertices(), n);
                }
            }
        }
    }
}

#[test]
fn rotated_symmetric_clouds_keep_their_hull() {
    // A rotated cube with extra points on its faces, edges and in its interior.
    let rotation = quat_from_axis_angle(&Vec3::new(1.0, 2.0, 3.0).normalize(), 0.7);
    let mut cloud = Vec::new();
    for x in [-1.0, 1.0] {
        for y in [-1.0, 1.0] {
            for z in [-1.0, 1.0] {
                cloud.push(Vec3::new(x, y, z));
            }
        }
    }
    cloud.extend([
        Vec3::new(0.0, 0.0, 0.0),
        Vec3::new(1.0, 0.0, 0.0),
        Vec3::new(0.3, -0.2, 1.0),
        Vec3::new(1.0, 1.0, 0.0),
        Vec3::new(0.0, -1.0, -1.0),
    ]);
    let cloud: Vec<Vec3> = cloud.iter().map(|p| rotation * p).collect();
    for order in [cloud.clone(), sorted(&cloud), shuffled(&cloud, 3)] {
        let hull = ConvexHull3DVec::build(&order).expect("rotated cube hull");
        assert_valid_hull(&order, &hull);
        assert!((hull.volume() - 8.0).abs() < 1e-9);
        for corner in cloud.iter().take(8) {
            assert!(hull.vertices.iter().any(|v| (v - corner).norm() < 1e-12));
        }
    }
}
