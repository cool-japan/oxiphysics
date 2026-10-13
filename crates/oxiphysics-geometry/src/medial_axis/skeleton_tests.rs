// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Regression tests for the 3-D medial-axis constructions
//! ([`super::MedialAxis::from_points`], [`super::CurveSkeleton`] and
//! [`super::laplacian_contraction`]) on shapes with a known medial axis.

use super::*;
use std::f64::consts::PI;

/// Deterministic xorshift64* generator for reproducible jitter.
struct XorShift(u64);

impl XorShift {
    fn next_f64(&mut self) -> f64 {
        let mut x = self.0;
        x ^= x >> 12;
        x ^= x << 25;
        x ^= x >> 27;
        self.0 = x;
        (x.wrapping_mul(0x2545_F491_4F6C_DD1D) >> 11) as f64 / (1u64 << 53) as f64
    }

    /// Uniform in `[-amp, amp]`.
    fn jitter(&mut self, amp: f64) -> f64 {
        (2.0 * self.next_f64() - 1.0) * amp
    }
}

/// Fibonacci-sphere directions restricted to `z >= 0` (if `upper`) or
/// `z <= 0`.
fn hemisphere_dirs(count: usize, upper: bool) -> Vec<[f64; 3]> {
    let total = 2 * count;
    let golden = PI * (3.0 - 5.0_f64.sqrt());
    (0..total)
        .filter_map(|i| {
            let z = 1.0 - (i as f64 + 0.5) * 2.0 / total as f64;
            let r = (1.0 - z * z).sqrt();
            let th = golden * i as f64;
            let d = [r * th.cos(), r * th.sin(), z];
            if (upper && z > 0.0) || (!upper && z < 0.0) {
                Some(d)
            } else {
                None
            }
        })
        .collect()
}

/// Jittered surface sample of a capsule of radius `r` around the segment
/// `(0,0,0)–(0,0,len)`. Its medial axis is exactly that segment with
/// constant radius `r`.
fn capsule(r: f64, len: f64, spacing: f64, seed: u64) -> Vec<[f64; 3]> {
    let mut rng = XorShift(seed);
    let mut pts = Vec::new();
    let around = ((2.0 * PI * r) / spacing).round().max(6.0) as usize;
    let rings = (len / spacing).round().max(1.0) as usize;
    let golden = PI * (3.0 - 5.0_f64.sqrt());
    for k in 0..=rings {
        let z = len * k as f64 / rings as f64;
        for j in 0..around {
            let th = 2.0 * PI * j as f64 / around as f64 + golden * k as f64;
            let rr = r * (1.0 + rng.jitter(2e-3));
            pts.push([rr * th.cos(), rr * th.sin(), z + rng.jitter(2e-3 * len)]);
        }
    }
    let cap_count = ((2.0 * PI * r * r) / (spacing * spacing)).round() as usize;
    for (upper, cz) in [(true, len), (false, 0.0)] {
        for d in hemisphere_dirs(cap_count, upper) {
            let rr = r * (1.0 + rng.jitter(2e-3));
            pts.push([d[0] * rr, d[1] * rr, cz + d[2] * rr]);
        }
    }
    pts
}

/// Distance from `p` to the segment `a–b`.
fn seg_dist(p: [f64; 3], a: [f64; 3], b: [f64; 3]) -> f64 {
    point_to_segment_dist(p, a, b)
}

/// Fibonacci sample of a sphere of radius `r` centred at `c`.
fn sphere(c: [f64; 3], r: f64, count: usize) -> Vec<[f64; 3]> {
    let golden = PI * (3.0 - 5.0_f64.sqrt());
    (0..count)
        .map(|i| {
            let z = 1.0 - (i as f64 + 0.5) * 2.0 / count as f64;
            let rho = (1.0 - z * z).sqrt();
            let th = golden * i as f64;
            [
                c[0] + r * rho * th.cos(),
                c[1] + r * rho * th.sin(),
                c[2] + r * z,
            ]
        })
        .collect()
}

/// Jittered sample of a torus with major radius `big_r` and tube radius
/// `r` around the z axis. Its medial axis is the core circle.
fn torus(big_r: f64, r: f64, nu: usize, nv: usize, seed: u64) -> Vec<[f64; 3]> {
    let mut rng = XorShift(seed);
    let mut pts = Vec::with_capacity(nu * nv);
    let golden = PI * (3.0 - 5.0_f64.sqrt());
    for i in 0..nu {
        let u = 2.0 * PI * i as f64 / nu as f64;
        for j in 0..nv {
            let v = 2.0 * PI * j as f64 / nv as f64 + golden * i as f64;
            let rr = r * (1.0 + rng.jitter(2e-3));
            let uu = u + rng.jitter(1e-3);
            let rho = big_r + rr * v.cos();
            pts.push([rho * uu.cos(), rho * uu.sin(), rr * v.sin()]);
        }
    }
    pts
}

/// Surface sample of the union of two capsules of radius `r` around the
/// segments `(-3,0,0)–(3,0,0)` and `(0,0,0)–(0,0,3)`: a T-shaped solid whose
/// curve skeleton has a junction at the origin.
fn t_shape(r: f64, spacing: f64) -> Vec<[f64; 3]> {
    let bar: Vec<[f64; 3]> = capsule(r, 6.0, spacing, 11)
        .into_iter()
        .map(|p| [p[2] - 3.0, p[1], p[0]])
        .collect();
    let stem = capsule(r, 3.0, spacing, 13);
    let in_bar = |p: [f64; 3]| seg_dist(p, [-3.0, 0.0, 0.0], [3.0, 0.0, 0.0]) < r * 0.999;
    let in_stem = |p: [f64; 3]| seg_dist(p, [0.0, 0.0, 0.0], [0.0, 0.0, 3.0]) < r * 0.999;
    bar.into_iter()
        .filter(|&p| !in_stem(p))
        .chain(stem.into_iter().filter(|&p| !in_bar(p)))
        .collect()
}

const CAPSULE_A: [f64; 3] = [0.0, 0.0, 0.0];
const CAPSULE_B: [f64; 3] = [0.0, 0.0, 4.0];

/// Assert that a skeleton of the standard capsule lies on its axis segment,
/// has the capsule radius, covers the segment and is a single branch.
fn assert_capsule_skeleton(ma: &MedialAxis, what: &str, max_dev: f64, radius_tol: f64) {
    assert!(!ma.points.is_empty(), "{what}: empty skeleton");
    for p in &ma.points {
        let d = seg_dist(p.position, CAPSULE_A, CAPSULE_B);
        assert!(
            d <= max_dev,
            "{what}: point {:?} is {d} from the axis",
            p.position
        );
        assert!(
            (p.radius - 1.0).abs() <= radius_tol,
            "{what}: radius {} at {:?}",
            p.radius,
            p.position
        );
        assert!((p.feature_size - p.radius).abs() < 1e-12);
    }
    let zmin = ma
        .points
        .iter()
        .map(|p| p.position[2])
        .fold(f64::INFINITY, f64::min);
    let zmax = ma
        .points
        .iter()
        .map(|p| p.position[2])
        .fold(f64::NEG_INFINITY, f64::max);
    assert!(
        zmin <= 0.4 && zmax >= 3.6,
        "{what}: covers only z in [{zmin}, {zmax}]"
    );
    assert_eq!(ma.branch_count(), 1, "{what}: branches {:?}", ma.branches);
    assert!(
        ma.junctions.is_empty(),
        "{what}: junctions {:?}",
        ma.junctions
    );
    assert!(ma.is_tubular(), "{what}: not tubular");
    let len = ma.total_length();
    assert!((3.4..=5.5).contains(&len), "{what}: length {len}");
    // Every interior point of the branch has a unit tangent along the axis.
    let branch = &ma.branches[0];
    let axial = branch[1..branch.len() - 1]
        .iter()
        .filter(|&&i| ma.points[i].tangent[2].abs() > 0.7)
        .count();
    assert!(
        axial * 10 >= (branch.len() - 2) * 7,
        "{what}: tangents are not axial"
    );
}

#[test]
fn voronoi_skeleton_of_capsule_is_its_axis_segment() {
    let pts = capsule(1.0, 4.0, 0.25, 7);
    let ma = MedialAxis::from_points(&pts, SkeletonMethod::Voronoi).expect("voronoi");
    assert_capsule_skeleton(&ma, "voronoi", 0.15, 0.12);
}

#[test]
fn ball_rolling_skeleton_of_capsule_is_its_axis_segment() {
    let pts = capsule(1.0, 4.0, 0.25, 7);
    let ma = MedialAxis::from_points(&pts, SkeletonMethod::BallRolling).expect("ball rolling");
    assert_capsule_skeleton(&ma, "ball rolling", 0.05, 0.05);
}

#[test]
fn thinning_skeleton_of_capsule_is_its_axis_segment() {
    let pts = capsule(1.0, 4.0, 0.25, 7);
    let ma = MedialAxis::from_points(&pts, SkeletonMethod::Thinning).expect("thinning");
    // Voxel size is 6 / 64 ≈ 0.094; allow two voxels of deviation.
    assert_capsule_skeleton(&ma, "thinning", 0.2, 0.2);
}

#[test]
fn medial_surface_of_capsule_samples_the_axis_balls() {
    let pts = capsule(1.0, 4.0, 0.25, 7);
    let ma = MedialAxis::from_points(&pts, SkeletonMethod::MedialSurface).expect("surface");
    assert!(ma.points.len() > 100, "only {} poles", ma.points.len());
    assert!(ma.branches.is_empty() && ma.junctions.is_empty());
    for p in &ma.points {
        assert!(
            seg_dist(p.position, CAPSULE_A, CAPSULE_B) < 0.2,
            "{:?}",
            p.position
        );
        assert!((p.radius - 1.0).abs() < 0.2, "radius {}", p.radius);
    }
    // The union of the medial balls reconstructs the capsule: every sample
    // is within 5 % of the radius of some ball's boundary.
    for s in &pts {
        let best = ma
            .points
            .iter()
            .map(|p| (dist3(*s, p.position) - p.radius).abs())
            .fold(f64::INFINITY, f64::min);
        assert!(best < 0.05, "sample {s:?} not on any medial ball ({best})");
    }
}

#[test]
fn every_method_reduces_a_sphere_to_its_centre() {
    let c = [1.0, -2.0, 0.5];
    let pts = sphere(c, 2.0, 500);
    for method in [
        SkeletonMethod::Voronoi,
        SkeletonMethod::MedialSurface,
        SkeletonMethod::BallRolling,
        SkeletonMethod::Thinning,
    ] {
        let ma = MedialAxis::from_points(&pts, method).expect("sphere");
        let tol = if method == SkeletonMethod::Thinning {
            0.1
        } else {
            1e-3
        };
        assert_eq!(ma.points.len(), 1, "{method:?}: {} points", ma.points.len());
        let p = &ma.points[0];
        assert!(
            dist3(p.position, c) < tol,
            "{method:?}: centre {:?}",
            p.position
        );
        assert!(
            (p.radius - 2.0).abs() < tol,
            "{method:?}: radius {}",
            p.radius
        );
        assert!(ma.junctions.is_empty());
        assert_eq!(ma.total_length(), 0.0);
    }
}

#[test]
fn skeleton_of_t_shape_has_a_junction_at_the_crossing() {
    let pts = t_shape(0.6, 0.15);
    for method in [
        SkeletonMethod::Voronoi,
        SkeletonMethod::BallRolling,
        SkeletonMethod::Thinning,
    ] {
        let ma = MedialAxis::from_points(&pts, method).expect("t shape");
        assert!(!ma.junctions.is_empty(), "{method:?}: no junction");
        assert!(
            ma.branch_count() >= 3,
            "{method:?}: {} branches",
            ma.branch_count()
        );
        let near_origin = ma
            .junctions
            .iter()
            .any(|&j| dist3(ma.points[j].position, [0.0; 3]) < 0.6);
        assert!(near_origin, "{method:?}: junctions {:?}", ma.junctions);
        let ext = |k: usize, f: fn(f64, f64) -> f64, init: f64| {
            ma.points.iter().map(|p| p.position[k]).fold(init, f)
        };
        assert!(
            ext(0, f64::max, f64::NEG_INFINITY) > 2.5,
            "{method:?}: +x arm"
        );
        assert!(ext(0, f64::min, f64::INFINITY) < -2.5, "{method:?}: -x arm");
        assert!(
            ext(2, f64::max, f64::NEG_INFINITY) > 2.5,
            "{method:?}: +z arm"
        );
        for p in &ma.points {
            let d = seg_dist(p.position, [-3.0, 0.0, 0.0], [3.0, 0.0, 0.0]).min(seg_dist(
                p.position,
                [0.0; 3],
                [0.0, 0.0, 3.0],
            ));
            assert!(
                d < 0.3,
                "{method:?}: point {:?} off the skeleton by {d}",
                p.position
            );
        }
        // The skeleton graph of a T is a tree (χ = 1).
        assert_eq!(euler_characteristic(&ma), 1, "{method:?}");
    }
}

#[test]
fn skeleton_of_torus_is_one_closed_loop() {
    let pts = torus(2.0, 0.5, 90, 22, 5);
    for method in [
        SkeletonMethod::Voronoi,
        SkeletonMethod::BallRolling,
        SkeletonMethod::Thinning,
    ] {
        let ma = MedialAxis::from_points(&pts, method).expect("torus");
        assert!(
            ma.junctions.is_empty(),
            "{method:?}: junctions {:?}",
            ma.junctions
        );
        assert_eq!(
            ma.branch_count(),
            1,
            "{method:?}: branches {:?}",
            ma.branches
        );
        let b = &ma.branches[0];
        assert_eq!(b.first(), b.last(), "{method:?}: loop is not closed");
        // A single cycle has χ = V − E = 0.
        assert_eq!(euler_characteristic(&ma), 0, "{method:?}");
        for p in &ma.points {
            let rho = (p.position[0].powi(2) + p.position[1].powi(2)).sqrt();
            let d = ((rho - 2.0).powi(2) + p.position[2].powi(2)).sqrt();
            assert!(
                d < 0.12,
                "{method:?}: point {:?} is {d} off the core circle",
                p.position
            );
            assert!(
                (p.radius - 0.5).abs() < 0.1,
                "{method:?}: radius {}",
                p.radius
            );
        }
        let len = ma.total_length();
        let circle = 2.0 * PI * 2.0;
        assert!(
            (len - circle).abs() < 0.15 * circle,
            "{method:?}: loop length {len} vs {circle}"
        );
    }
}

#[test]
fn medial_axis_rejects_bad_input_and_config() {
    let pts = capsule(1.0, 4.0, 0.25, 7);
    let mut bad = pts.clone();
    bad[3][1] = f64::NAN;
    assert!(matches!(
        MedialAxis::from_points(&bad, SkeletonMethod::Voronoi),
        Err(crate::error::Error::DegenerateGeometry { .. })
    ));
    let configs = [
        MedialAxisConfig {
            min_radius_factor: -1.0,
            ..MedialAxisConfig::default()
        },
        MedialAxisConfig {
            prune_scale: f64::NAN,
            ..MedialAxisConfig::default()
        },
        MedialAxisConfig {
            loop_factor: 0.0,
            ..MedialAxisConfig::default()
        },
        MedialAxisConfig {
            thinning_resolution: 2,
            ..MedialAxisConfig::default()
        },
    ];
    for cfg in configs {
        assert!(matches!(
            MedialAxis::from_points_with_config(&pts, SkeletonMethod::Voronoi, &cfg),
            Err(crate::error::Error::OutOfRange { .. })
        ));
    }
    // A resolution filter larger than any medial ball leaves nothing.
    let strict = MedialAxisConfig {
        min_radius_factor: 1e6,
        ..MedialAxisConfig::default()
    };
    assert!(matches!(
        MedialAxis::from_points_with_config(&pts, SkeletonMethod::Voronoi, &strict),
        Err(crate::error::Error::DegenerateGeometry { .. })
    ));
}

#[test]
fn duplicate_samples_do_not_change_the_skeleton() {
    let pts = capsule(1.0, 4.0, 0.25, 7);
    let mut doubled = pts.clone();
    doubled.extend(pts.iter().take(50).copied());
    let a = MedialAxis::from_points(&pts, SkeletonMethod::BallRolling).expect("a");
    let b = MedialAxis::from_points(&doubled, SkeletonMethod::BallRolling).expect("b");
    assert_capsule_skeleton(&b, "ball rolling with duplicates", 0.05, 0.05);
    assert!(hausdorff_distance(&a, &b) < 0.3);
}

#[test]
fn prune_scale_controls_twig_removal() {
    // With twig pruning disabled (scale 0) the T keeps at least the junction;
    // with an enormous scale every twig hanging off a junction is swallowed
    // by the junction ball, leaving a junction-free skeleton.
    let pts = t_shape(0.6, 0.15);
    let loose = MedialAxisConfig {
        prune_scale: 0.0,
        ..MedialAxisConfig::default()
    };
    let tight = MedialAxisConfig {
        prune_scale: 100.0,
        ..MedialAxisConfig::default()
    };
    let a = MedialAxis::from_points_with_config(&pts, SkeletonMethod::Voronoi, &loose).expect("a");
    let b = MedialAxis::from_points_with_config(&pts, SkeletonMethod::Voronoi, &tight).expect("b");
    assert!(!a.junctions.is_empty());
    assert!(b.junctions.is_empty(), "junctions {:?}", b.junctions);
    assert!(b.total_length() < a.total_length());
}

#[test]
fn laplacian_contraction_moves_the_capsule_onto_its_axis() {
    let pts = capsule(1.0, 4.0, 0.25, 7);
    let mean_axis_distance = |c: &[[f64; 3]]| {
        c.iter()
            .map(|p| (p[0] * p[0] + p[1] * p[1]).sqrt())
            .sum::<f64>()
            / c.len() as f64
    };
    let mut previous = mean_axis_distance(&pts);
    for steps in [1, 2, 4, 8] {
        let c = laplacian_contraction(&pts, steps).expect("contract");
        assert_eq!(c.len(), pts.len());
        let d = mean_axis_distance(&c);
        assert!(d < previous, "{steps} steps: {d} >= {previous}");
        previous = d;
    }
    assert!(
        previous < 0.05,
        "8 steps leave mean axis distance {previous}"
    );
    let c = laplacian_contraction(&pts, 8).expect("contract");
    for p in &c {
        assert!(
            seg_dist(*p, [0.0, 0.0, -0.15], [0.0, 0.0, 4.15]) < 0.1,
            "{p:?}"
        );
    }
    assert!(matches!(
        laplacian_contraction(&pts, 0),
        Err(crate::error::Error::OutOfRange { .. })
    ));
}

#[test]
fn curve_skeleton_of_capsule_is_a_path_along_the_axis() {
    let pts = capsule(1.0, 4.0, 0.25, 7);
    for steps in [6, 8] {
        let cs = CurveSkeleton::from_surface_points(&pts, steps).expect("curve skeleton");
        assert!(cs.node_count() >= 4, "{steps}: {} nodes", cs.node_count());
        assert_eq!(cs.leaf_nodes().len(), 2, "{steps}");
        assert!(cs.junction_nodes().is_empty(), "{steps}");
        assert_eq!(cs.edge_count(), cs.node_count() - 1, "{steps}: not a path");
        for nd in &cs.nodes {
            // Nodes sit at the centres of the fitted medial balls: on the
            // axis segment (end nodes at the hemisphere centres).
            let d = seg_dist(nd.position, [0.0, 0.0, -0.25], [0.0, 0.0, 4.25]);
            assert!(d < 0.05, "{steps}: node {:?}", nd.position);
            assert!(
                (nd.radius - 1.0).abs() < 0.08,
                "{steps}: radius {}",
                nd.radius
            );
        }
        let len = cs.total_length();
        assert!((3.5..=6.0).contains(&len), "{steps}: length {len}");
    }
}

#[test]
fn curve_skeleton_of_torus_is_a_cycle() {
    let pts = torus(2.0, 0.5, 72, 18, 9);
    let cs = CurveSkeleton::from_surface_points(&pts, 8).expect("torus skeleton");
    assert!(cs.leaf_nodes().is_empty(), "leaves {:?}", cs.leaf_nodes());
    assert!(cs.junction_nodes().is_empty());
    assert_eq!(cs.edge_count(), cs.node_count(), "a cycle has V = E");
    for nd in &cs.nodes {
        let rho = (nd.position[0].powi(2) + nd.position[1].powi(2)).sqrt();
        let d = ((rho - 2.0).powi(2) + nd.position[2].powi(2)).sqrt();
        assert!(
            d < 0.05,
            "node {:?} is {d} off the core circle",
            nd.position
        );
        assert!((nd.radius - 0.5).abs() < 0.05, "radius {}", nd.radius);
    }
    let circle = 2.0 * PI * 2.0;
    assert!((cs.total_length() - circle).abs() < 0.1 * circle);
}

#[test]
fn curve_skeleton_of_t_shape_has_one_junction_and_three_arms() {
    let pts = t_shape(0.6, 0.15);
    let cs = CurveSkeleton::from_surface_points(&pts, 10).expect("t skeleton");
    let junctions = cs.junction_nodes();
    assert_eq!(junctions.len(), 1, "junctions {junctions:?}");
    assert!(dist3(cs.nodes[junctions[0]].position, [0.0; 3]) < 0.5);
    let leaves = cs.leaf_nodes();
    assert_eq!(leaves.len(), 3, "leaves {leaves:?}");
    // One leaf per arm.
    let tips: Vec<[f64; 3]> = leaves.iter().map(|&l| cs.nodes[l].position).collect();
    assert!(tips.iter().any(|p| p[0] > 2.3));
    assert!(tips.iter().any(|p| p[0] < -2.3));
    assert!(tips.iter().any(|p| p[2] > 2.3));
    assert_eq!(cs.edge_count(), cs.node_count() - 1, "a T is a tree");
}
