// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Regression tests for the conservative-advancement shape casts.
//!
//! The oblique-approach tests pin the time of impact to its analytic value: the
//! original `d / |vel|` advancement shrank the gap by only `1 - cos(theta)` per
//! step and ran out of iterations before reaching a contact at 60 to 85 degrees
//! from the normal. The box, capsule and convex-hull casts are also checked
//! against independent exact oracles (slab intersection, a 13-axis
//! separating-axis interval test, feature-pair distances).

use super::advancement::{Probe, advance, separation_direction, unit_vector};
use super::box_queries::{aabb_separation, box_corners, segment_box_distance};
use super::hull_distance::{HullDistance, hull_distance};
use super::*;

// ─────────────────────────────────────────────────────────────────────────────
// Helpers
// ─────────────────────────────────────────────────────────────────────────────

fn close(a: f64, b: f64, tol: f64) -> bool {
    (a - b).abs() <= tol
}

fn vec_close(a: [f64; 3], b: [f64; 3], tol: f64) -> bool {
    len3(sub3(a, b)) <= tol
}

fn default_cfg() -> ShapeCastConfig {
    ShapeCastConfig::default()
}

fn cfg_until(max_toc: f64) -> ShapeCastConfig {
    ShapeCastConfig {
        max_toc,
        ..ShapeCastConfig::default()
    }
}

/// A mover that reaches `contact_pos` at time `contact_time` while travelling at
/// `angle_deg` from the -y direction (in the xy plane) with the given `speed`.
fn oblique_motion(
    contact_pos: [f64; 3],
    angle_deg: f64,
    speed: f64,
    contact_time: f64,
) -> LinearMotion {
    let angle = angle_deg.to_radians();
    let vel = [speed * angle.sin(), -speed * angle.cos(), 0.0];
    LinearMotion::new(sub3(contact_pos, scale3(vel, contact_time)), vel)
}

/// Deterministic xorshift generator for the randomised oracle tests.
struct Rng(u64);

impl Rng {
    fn new(seed: u64) -> Self {
        Rng(seed.max(1))
    }

    fn next_u64(&mut self) -> u64 {
        let mut x = self.0;
        x ^= x >> 12;
        x ^= x << 25;
        x ^= x >> 27;
        self.0 = x;
        x.wrapping_mul(0x2545_F491_4F6C_DD1D)
    }

    /// Uniform in `[0, 1)`.
    fn unit(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
    }

    fn range(&mut self, lo: f64, hi: f64) -> f64 {
        lo + (hi - lo) * self.unit()
    }

    fn vec(&mut self, lo: f64, hi: f64) -> [f64; 3] {
        [self.range(lo, hi), self.range(lo, hi), self.range(lo, hi)]
    }
}

/// Exact first-contact time of two axis-aligned boxes (the first one moving) by
/// slab intersection: the time interval during which the boxes overlap on every
/// axis.
fn slab_toi(
    center_a: [f64; 3],
    half_a: [f64; 3],
    vel: [f64; 3],
    center_b: [f64; 3],
    half_b: [f64; 3],
    max_t: f64,
) -> Option<f64> {
    let mut t_in = f64::NEG_INFINITY;
    let mut t_out = f64::INFINITY;
    for k in 0..3 {
        let reach = half_a[k] + half_b[k];
        let offset = center_a[k] - center_b[k];
        if vel[k].abs() < 1e-300 {
            if offset.abs() > reach {
                return None;
            }
        } else {
            let t0 = (-reach - offset) / vel[k];
            let t1 = (reach - offset) / vel[k];
            t_in = t_in.max(t0.min(t1));
            t_out = t_out.min(t0.max(t1));
        }
    }
    if t_in > t_out || t_out < 0.0 || t_in > max_t {
        None
    } else {
        Some(t_in.max(0.0))
    }
}

/// Entry and exit time of a box translating by `vel` against a static triangle,
/// from the 13-axis separating-axis test (box axes, triangle normal and the nine
/// box-axis x triangle-edge cross products), each axis giving an interval of
/// times during which the projections overlap (Ericson, RTCD 5.5.8).
///
/// The shapes overlap during `[t_in, t_out]` if `t_in <= t_out`; a negative
/// duration measures how far from touching they stay. `None` means the motion is
/// parallel to a separating axis that never lets them overlap.
fn sat_box_triangle_interval(
    center: [f64; 3],
    half: [f64; 3],
    vel: [f64; 3],
    tri: [[f64; 3]; 3],
) -> Option<(f64, f64)> {
    let edges = [
        sub3(tri[1], tri[0]),
        sub3(tri[2], tri[1]),
        sub3(tri[0], tri[2]),
    ];
    let box_axes = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    let mut axes: Vec<[f64; 3]> = box_axes.to_vec();
    axes.push(cross3(edges[0], edges[1]));
    for axis in box_axes {
        for edge in edges {
            axes.push(cross3(axis, edge));
        }
    }
    let mut t_in = f64::NEG_INFINITY;
    let mut t_out = f64::INFINITY;
    for axis in axes {
        let norm = len3(axis);
        if norm < 1e-12 {
            continue;
        }
        let l = scale3(axis, 1.0 / norm);
        let c = dot3(center, l);
        let r = half[0] * l[0].abs() + half[1] * l[1].abs() + half[2] * l[2].abs();
        let proj = [dot3(tri[0], l), dot3(tri[1], l), dot3(tri[2], l)];
        let tri_min = proj.iter().copied().fold(f64::INFINITY, f64::min);
        let tri_max = proj.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        let u = dot3(vel, l);
        if u.abs() < 1e-12 {
            if c - r > tri_max || c + r < tri_min {
                return None;
            }
        } else {
            let t0 = (tri_min - r - c) / u;
            let t1 = (tri_max + r - c) / u;
            t_in = t_in.max(t0.min(t1));
            t_out = t_out.min(t0.max(t1));
        }
    }
    Some((t_in, t_out))
}

/// Exact distance between two disjoint triangles: the minimum over all
/// vertex-face and edge-edge feature pairs.
fn triangle_triangle_distance(a: [[f64; 3]; 3], b: [[f64; 3]; 3]) -> f64 {
    let mut best = f64::INFINITY;
    for p in a {
        let (q, _) = crate_point_tri(p, b[0], b[1], b[2]);
        best = best.min(dist3(p, q));
    }
    for p in b {
        let (q, _) = crate_point_tri(p, a[0], a[1], a[2]);
        best = best.min(dist3(p, q));
    }
    for i in 0..3 {
        for j in 0..3 {
            let (d, _, _) = seg_seg_dist(a[i], a[(i + 1) % 3], b[j], b[(j + 1) % 3]);
            best = best.min(d);
        }
    }
    best
}

fn translated(points: &[[f64; 3]], offset: [f64; 3]) -> Vec<[f64; 3]> {
    points.iter().map(|&p| add3(p, offset)).collect()
}

// ─────────────────────────────────────────────────────────────────────────────
// Advancement driver
// ─────────────────────────────────────────────────────────────────────────────

/// A point 3 units above the plane y = 0 with radius 0.5 approaches it at any
/// angle: the gap is linear in time, so the Newton step lands on the contact
/// after one step (two probes: the start and the contact).
#[test]
fn advance_reaches_a_planar_contact_in_one_step_at_any_angle() {
    for angle_deg in [0.0_f64, 60.0, 80.0, 85.0, 89.0] {
        let angle = angle_deg.to_radians();
        let vel = [10.0 * angle.sin(), -10.0 * angle.cos(), 0.0];
        let mut probes = 0usize;
        let result = advance(vel, &cfg_until(100.0), |t| {
            probes += 1;
            let y = 3.0 + vel[1] * t;
            Some(Probe {
                gap: y - 0.5,
                direction: Some([0.0, 1.0, 0.0]),
                fallback: None,
                witness_a: [0.0, y, 0.0],
                witness_b: [0.0; 3],
            })
        });
        let expected = 2.5 / (10.0 * angle.cos());
        assert!(result.hit, "angle {angle_deg}");
        assert!(
            close(result.toi, expected, 1e-9 * expected.max(1.0)),
            "angle {angle_deg}: toi {} expected {expected}",
            result.toi
        );
        assert_eq!(probes, 2, "angle {angle_deg}");
    }
}

#[test]
fn advance_misses_when_the_gap_is_not_closing() {
    let probe = |_t: f64| {
        Some(Probe {
            gap: 1.0,
            direction: Some([0.0, 1.0, 0.0]),
            fallback: None,
            witness_a: [0.0, 1.0, 0.0],
            witness_b: [0.0; 3],
        })
    };
    // Moving away from the plane.
    assert!(!advance([1.0, 0.5, 0.0], &cfg_until(10.0), probe).hit);
    // Sliding parallel to the plane.
    assert!(!advance([1.0, 0.0, 0.0], &cfg_until(10.0), probe).hit);
    // Not moving at all.
    assert!(!advance([0.0; 3], &cfg_until(10.0), probe).hit);
}

#[test]
fn advance_ends_without_a_hit_when_nothing_can_be_measured() {
    let none = advance([0.0, -1.0, 0.0], &cfg_until(1.0), |_t| None);
    assert!(!none.hit);
    for gap in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        let result = advance([0.0, -1.0, 0.0], &cfg_until(1.0), |_t| {
            Some(Probe {
                gap,
                direction: None,
                fallback: None,
                witness_a: [0.0; 3],
                witness_b: [0.0; 3],
            })
        });
        assert!(!result.hit, "gap {gap}");
    }
}

#[test]
fn advance_keeps_stepping_on_the_last_direction_when_the_direction_is_lost() {
    // The direction is only known at t = 0; afterwards the closest features are
    // reported as coincident although the gap is still above the tolerance.
    let vel = [0.0, -1.0, 0.0];
    let result = advance(vel, &cfg_until(10.0), |t| {
        let y = 2.0 + vel[1] * t;
        Some(Probe {
            gap: y,
            direction: if t == 0.0 {
                Some([0.0, 1.0, 0.0])
            } else {
                None
            },
            fallback: Some([1.0, 0.0, 0.0]),
            witness_a: [0.0, y, 0.0],
            witness_b: [0.0; 3],
        })
    });
    assert!(result.hit);
    assert!(close(result.toi, 2.0, 1e-9), "toi={}", result.toi);
    // The earlier direction wins over the probe's fallback.
    assert!(vec_close(result.normal_at_toi, [0.0, 1.0, 0.0], 1e-12));
}

#[test]
fn advance_respects_target_distance_and_tolerance() {
    let vel = [0.0, -1.0, 0.0];
    let probe_at = |t: f64| {
        let y = 5.0 + vel[1] * t;
        Some(Probe {
            gap: y,
            direction: Some([0.0, 1.0, 0.0]),
            fallback: None,
            witness_a: [0.0, y, 0.0],
            witness_b: [0.0; 3],
        })
    };
    let target = ShapeCastConfig {
        max_toc: 10.0,
        target_distance: 0.75,
        toi_tolerance: 0.0,
    };
    assert!(close(advance(vel, &target, probe_at).toi, 4.25, 1e-12));
    // Already within the target separation: a hit at t = 0, not a penetration.
    let near = advance(
        vel,
        &ShapeCastConfig {
            max_toc: 10.0,
            target_distance: 6.0,
            toi_tolerance: 0.0,
        },
        probe_at,
    );
    assert!(near.hit && near.toi == 0.0 && !near.is_penetrating);
    // A hit contact beyond max_toc is a miss, a contact exactly at max_toc a hit.
    assert!(!advance(vel, &cfg_until(4.9), probe_at).hit);
    assert!(advance(vel, &cfg_until(5.0), probe_at).hit);
}

#[test]
fn advance_flags_initial_overlap_as_penetrating() {
    let result = advance([0.0, -1.0, 0.0], &cfg_until(1.0), |_t| {
        Some(Probe {
            gap: -0.25,
            direction: None,
            fallback: Some([0.0, 0.0, 1.0]),
            witness_a: [0.0; 3],
            witness_b: [0.0; 3],
        })
    });
    assert!(result.hit && result.is_penetrating);
    assert_eq!(result.toi, 0.0);
    assert!(vec_close(result.normal_at_toi, [0.0, 0.0, 1.0], 1e-12));
}

#[test]
fn separation_direction_ignores_rounding_noise() {
    let far = [1.0e3, 2.0e3, 3.0e3];
    assert!(separation_direction(far, far).is_none());
    assert!(separation_direction(add3(far, [1e-9, 0.0, 0.0]), far).is_none());
    let n = separation_direction(add3(far, [0.0, 1e-3, 0.0]), far);
    assert!(n.is_some_and(|n| vec_close(n, [0.0, 1.0, 0.0], 1e-12)));
    assert!(unit_vector([0.0; 3]).is_none());
    assert!(unit_vector([f64::NAN, 0.0, 0.0]).is_none());
}

// ─────────────────────────────────────────────────────────────────────────────
// Sphere casts
// ─────────────────────────────────────────────────────────────────────────────

#[test]
fn sphere_vs_box_oblique_approach_hits_at_the_analytic_time() {
    // Sphere of radius 0.5 reaches the top face (y = 1) of a large box when its
    // centre is at y = 1.5; the approach is planar, so the time is exact.
    for angle_deg in [60.0, 80.0, 85.0] {
        let motion = oblique_motion([0.0, 1.5, 0.0], angle_deg, 10.0, 0.4);
        let result = SphereCast::vs_box(
            &motion,
            0.5,
            [0.0; 3],
            [5.0, 1.0, 5.0],
            &ShapeCastConfig::default(),
        );
        assert!(result.hit && !result.is_penetrating, "angle {angle_deg}");
        assert!(
            close(result.toi, 0.4, 1e-9),
            "angle {angle_deg}: {}",
            result.toi
        );
        assert!(vec_close(result.normal_at_toi, [0.0, 1.0, 0.0], 1e-9));
        assert!(vec_close(result.witness_b, [0.0, 1.0, 0.0], 1e-6));
        assert!(vec_close(result.witness_a, [0.0, 1.0, 0.0], 1e-6));
    }
}

#[test]
fn sphere_vs_box_centre_inside_uses_the_nearest_face_normal() {
    let motion = LinearMotion::new([0.9, 0.0, 0.0], [0.0, -1.0, 0.0]);
    let result = SphereCast::vs_box(
        &motion,
        0.25,
        [0.0; 3],
        [1.0, 1.0, 1.0],
        &ShapeCastConfig::default(),
    );
    assert!(result.hit && result.is_penetrating);
    assert_eq!(result.toi, 0.0);
    assert!(vec_close(result.normal_at_toi, [1.0, 0.0, 0.0], 1e-12));
}

#[test]
fn sphere_vs_triangle_oblique_approach_hits_at_the_analytic_time() {
    // Large triangle in the plane y = 0 containing the origin.
    let (v0, v1, v2) = ([-50.0, 0.0, -50.0], [50.0, 0.0, -50.0], [0.0, 0.0, 50.0]);
    for angle_deg in [60.0, 80.0, 85.0] {
        let motion = oblique_motion([0.0, 0.5, 0.0], angle_deg, 10.0, 0.4);
        let result = SphereCast::vs_triangle(&motion, 0.5, v0, v1, v2, &ShapeCastConfig::default());
        assert!(result.hit, "angle {angle_deg}");
        assert!(
            close(result.toi, 0.4, 1e-9),
            "angle {angle_deg}: {}",
            result.toi
        );
        assert!(vec_close(result.normal_at_toi, [0.0, 1.0, 0.0], 1e-9));
        assert!(vec_close(result.witness_b, [0.0; 3], 1e-6));
    }
}

#[test]
fn sphere_vs_triangle_resting_on_the_plane_uses_the_side_it_came_from() {
    // Centre exactly in the triangle plane: the normal is the triangle normal,
    // oriented against the motion.
    let (v0, v1, v2) = ([-5.0, 0.0, -5.0], [5.0, 0.0, -5.0], [0.0, 0.0, 5.0]);
    let from_above = LinearMotion::new([0.0; 3], [0.0, -1.0, 0.0]);
    let up = SphereCast::vs_triangle(&from_above, 0.1, v0, v1, v2, &ShapeCastConfig::default());
    assert!(up.hit && up.is_penetrating);
    assert!(vec_close(up.normal_at_toi, [0.0, 1.0, 0.0], 1e-12));
    let from_below = LinearMotion::new([0.0; 3], [0.0, 1.0, 0.0]);
    let down = SphereCast::vs_triangle(&from_below, 0.1, v0, v1, v2, &ShapeCastConfig::default());
    assert!(down.hit && down.is_penetrating);
    assert!(vec_close(down.normal_at_toi, [0.0, -1.0, 0.0], 1e-12));
}

#[test]
fn mesh_bvh_cast_keeps_a_hit_exactly_at_max_toc() {
    // The sphere touches the triangle at t = 1 = max_toc: a strict `<` against the
    // initial "no hit" record (toi = max_toc) used to discard it.
    let mesh = MeshBvhCast::new(
        vec![[-1.0, -1.0, 0.0], [1.0, -1.0, 0.0], [0.0, 1.0, 0.0]],
        vec![0, 1, 2],
    );
    let motion = LinearMotion::new([0.0, 0.0, -1.25], [0.0, 0.0, 1.0]);
    let result = mesh.sphere_cast(&motion, 0.25, &ShapeCastConfig::default());
    assert!(result.hit);
    assert_eq!(result.toi, 1.0);
    assert!(vec_close(result.normal_at_toi, [0.0, 0.0, -1.0], 1e-12));
}

// ─────────────────────────────────────────────────────────────────────────────
// Capsule casts
// ─────────────────────────────────────────────────────────────────────────────

#[test]
fn capsule_vs_capsule_oblique_approach_hits_at_the_analytic_time() {
    // Capsule A (axis along x, radius 0.25) sweeps towards the static capsule B
    // (axis along z through the origin, radius 0.25). The axes are 0.5 apart when
    // A's centre is at (0, 0.5, 0), reached at t = 0.4 at every angle.
    for angle_deg in [60.0, 80.0, 85.0] {
        let motion = oblique_motion([0.0, 0.5, 0.0], angle_deg, 10.0, 0.4);
        let result = CapsuleCast::vs_capsule(
            &motion,
            [-0.5, 0.0, 0.0],
            [0.5, 0.0, 0.0],
            0.25,
            [0.0, 0.0, -0.5],
            [0.0, 0.0, 0.5],
            0.25,
            &ShapeCastConfig::default(),
        );
        assert!(result.hit, "angle {angle_deg}");
        assert!(
            close(result.toi, 0.4, 1e-6),
            "angle {angle_deg}: {}",
            result.toi
        );
        assert!(vec_close(result.normal_at_toi, [0.0, 1.0, 0.0], 1e-6));
        assert!(vec_close(result.witness_a, [0.0, 0.5, 0.0], 1e-5));
        assert!(vec_close(result.witness_b, [0.0; 3], 1e-5));
    }
}

#[test]
fn capsule_vs_capsule_hits_with_the_default_config_when_the_gap_never_reaches_exactly_zero() {
    // With the default target distance of 0 the old loop needed `d <= 0`, which a
    // geometric shrink by 0.5 per step (60 degrees) never reaches within 64 steps.
    let motion = oblique_motion([0.0, 0.5, 0.0], 60.0, 1.0, 0.8);
    let result = CapsuleCast::vs_capsule(
        &motion,
        [-0.5, 0.0, 0.0],
        [0.5, 0.0, 0.0],
        0.25,
        [0.0, 0.0, -0.5],
        [0.0, 0.0, 0.5],
        0.25,
        &ShapeCastConfig::default(),
    );
    assert!(result.hit);
    assert!(close(result.toi, 0.8, 1e-6), "toi={}", result.toi);
}

#[test]
fn capsule_vs_box_detects_a_mid_segment_hit_on_a_narrow_box() {
    // A 10 unit long horizontal capsule falls onto a box that is only 0.2 wide in
    // x. Only the middle of the capsule hits; both end points stay about 5 units
    // away from the box in x, so an end-point-only distance never sees it.
    let motion = LinearMotion::new([0.0, 10.0, 0.0], [0.0, -10.0, 0.0]);
    let result = CapsuleCast::vs_box(
        &motion,
        [-5.0, 0.0, 0.0],
        [5.0, 0.0, 0.0],
        0.5,
        [0.0; 3],
        [0.1, 0.5, 0.5],
        &default_cfg(),
    );
    assert!(result.hit && !result.is_penetrating);
    // The axis touches y = 0.5 + 0.5 = 1.0 at t = 0.9.
    assert!(close(result.toi, 0.9, 1e-9), "toi={}", result.toi);
    assert!(vec_close(result.normal_at_toi, [0.0, 1.0, 0.0], 1e-9));
    // The axis lies flat over the 0.2 wide top face: any point of it within the
    // face is a valid witness, and the box point is the one straight below it.
    assert!(result.witness_a[0].abs() <= 0.1 + 1e-9);
    assert!(close(result.witness_a[1], 1.0, 1e-9) && close(result.witness_a[2], 0.0, 1e-9));
    assert!(vec_close(
        result.witness_b,
        [result.witness_a[0], 0.5, 0.0],
        1e-9
    ));
}

#[test]
fn capsule_vs_box_end_cap_hit_still_works() {
    // Vertical capsule (axis y in [-1, 1], radius 0.5) lands with its lower cap.
    let motion = LinearMotion::new([0.0, 10.0, 0.0], [0.0, -10.0, 0.0]);
    let result = CapsuleCast::vs_box(
        &motion,
        [0.0, -1.0, 0.0],
        [0.0, 1.0, 0.0],
        0.5,
        [0.0; 3],
        [1.0, 1.0, 1.0],
        &default_cfg(),
    );
    // The lower axis end point reaches y = 1.5 at t = (10 - 1 - 1.5) / 10.
    assert!(result.hit);
    assert!(close(result.toi, 0.75, 1e-9), "toi={}", result.toi);
    assert!(vec_close(result.normal_at_toi, [0.0, 1.0, 0.0], 1e-9));
    assert!(vec_close(result.witness_a, [0.0, 1.5, 0.0], 1e-9));
    assert!(vec_close(result.witness_b, [0.0, 1.0, 0.0], 1e-9));
}

#[test]
fn capsule_vs_box_edge_contact_has_a_diagonal_normal() {
    // Vertical capsule (radius 0.5) slides diagonally into a vertical box edge at
    // (1, y, 1): the contact normal is the horizontal diagonal and the witnesses
    // lie on the capsule axis and the edge.
    let motion = LinearMotion::new([3.0, 0.0, 3.0], [-2.0, 0.0, -2.0]);
    let result = CapsuleCast::vs_box(
        &motion,
        [0.0, -0.5, 0.0],
        [0.0, 0.5, 0.0],
        0.5,
        [0.0; 3],
        [1.0, 1.0, 1.0],
        &default_cfg(),
    );
    // The horizontal distance to the edge is sqrt(2) * (2 - 2t) = 0.5.
    let expected = 1.0 - 0.5 / (2.0 * 2.0_f64.sqrt());
    assert!(result.hit);
    assert!(close(result.toi, expected, 1e-9), "toi={}", result.toi);
    let diagonal = 0.5_f64.sqrt();
    assert!(vec_close(
        result.normal_at_toi,
        [diagonal, 0.0, diagonal],
        1e-9
    ));
    assert!(close(result.witness_b[0], 1.0, 1e-9) && close(result.witness_b[2], 1.0, 1e-9));
    let pos = motion.position_at(result.toi);
    assert!(close(result.witness_a[0], pos[0], 1e-9));
    assert!(close(result.witness_a[2], pos[2], 1e-9));
}

#[test]
fn capsule_vs_box_oblique_approach_hits_at_the_analytic_time() {
    // Horizontal capsule (axis along z, radius 0.25) comes down at 60 to 85 degrees
    // from the vertical onto the top face (y = 0.5) of a wide box; its axis is
    // 0.75 above the face at t = 0.4.
    for angle_deg in [60.0, 80.0, 85.0] {
        let motion = oblique_motion([0.0, 0.75, 0.0], angle_deg, 10.0, 0.4);
        let result = CapsuleCast::vs_box(
            &motion,
            [0.0, 0.0, -0.5],
            [0.0, 0.0, 0.5],
            0.25,
            [0.0; 3],
            [5.0, 0.5, 5.0],
            &ShapeCastConfig::default(),
        );
        assert!(result.hit, "angle {angle_deg}");
        assert!(
            close(result.toi, 0.4, 1e-9),
            "angle {angle_deg}: {}",
            result.toi
        );
        assert!(vec_close(result.normal_at_toi, [0.0, 1.0, 0.0], 1e-9));
    }
}

#[test]
fn capsule_vs_box_misses_and_penetration() {
    // Passes beside the box.
    let aside = LinearMotion::new([5.0, 10.0, 0.0], [0.0, -10.0, 0.0]);
    let miss = CapsuleCast::vs_box(
        &aside,
        [-0.5, 0.0, 0.0],
        [0.5, 0.0, 0.0],
        0.5,
        [0.0; 3],
        [1.0, 1.0, 1.0],
        &default_cfg(),
    );
    assert!(!miss.hit);
    // Starts with the axis inside the box: initial penetration, exit normal.
    let inside = LinearMotion::new([0.0, 0.8, 0.0], [0.0, -1.0, 0.0]);
    let overlap = CapsuleCast::vs_box(
        &inside,
        [-0.5, 0.0, 0.0],
        [0.5, 0.0, 0.0],
        0.25,
        [0.0; 3],
        [1.0, 1.0, 1.0],
        &default_cfg(),
    );
    assert!(overlap.hit && overlap.is_penetrating);
    assert!(vec_close(overlap.normal_at_toi, [0.0, 1.0, 0.0], 1e-12));
}

// ─────────────────────────────────────────────────────────────────────────────
// Segment-box distance
// ─────────────────────────────────────────────────────────────────────────────

/// Minimum distance from a segment to a box found by dense sampling followed by
/// repeated local refinement: an independent (inexact) oracle.
fn sampled_segment_box_distance(p0: [f64; 3], p1: [f64; 3], lo: [f64; 3], hi: [f64; 3]) -> f64 {
    let dist_at = |s: f64| {
        let p = add3(p0, scale3(sub3(p1, p0), s));
        let q = [
            p[0].clamp(lo[0], hi[0]),
            p[1].clamp(lo[1], hi[1]),
            p[2].clamp(lo[2], hi[2]),
        ];
        dist3(p, q)
    };
    let samples = 4000;
    let mut best_s = 0.0;
    let mut best = f64::INFINITY;
    for i in 0..=samples {
        let s = i as f64 / samples as f64;
        let d = dist_at(s);
        if d < best {
            best = d;
            best_s = s;
        }
    }
    let mut half_width = 1.0 / samples as f64;
    for _ in 0..4 {
        let (lo_s, hi_s) = (
            (best_s - half_width).max(0.0),
            (best_s + half_width).min(1.0),
        );
        for i in 0..=samples {
            let s = lo_s + (hi_s - lo_s) * i as f64 / samples as f64;
            let d = dist_at(s);
            if d < best {
                best = d;
                best_s = s;
            }
        }
        half_width /= samples as f64 / 2.0;
    }
    best
}

#[test]
fn segment_box_distance_matches_dense_sampling_on_random_inputs() {
    let mut rng = Rng::new(0x5EED_0001);
    for case in 0..400 {
        let center = rng.vec(-1.0, 1.0);
        let half = [
            rng.range(0.05, 1.0),
            rng.range(0.05, 1.0),
            rng.range(0.05, 1.0),
        ];
        let lo = sub3(center, half);
        let hi = add3(center, half);
        let p0 = rng.vec(-4.0, 4.0);
        let p1 = rng.vec(-4.0, 4.0);
        let exact = segment_box_distance(p0, p1, lo, hi);
        let sampled = sampled_segment_box_distance(p0, p1, lo, hi);
        assert!(
            exact.distance <= sampled + 1e-12,
            "case {case}: exact {} above sampled {sampled}",
            exact.distance
        );
        assert!(
            sampled - exact.distance <= 1e-7,
            "case {case}: exact {} far below sampled {sampled}",
            exact.distance
        );
        // The reported pair is consistent: on the segment and on the box.
        let d = sub3(p1, p0);
        let s = if dot3(d, d) > 0.0 {
            dot3(sub3(exact.segment_point, p0), d) / dot3(d, d)
        } else {
            0.0
        };
        assert!((-1e-12..=1.0 + 1e-12).contains(&s), "case {case}: s={s}");
        assert!(vec_close(
            exact.segment_point,
            add3(p0, scale3(d, s.clamp(0.0, 1.0))),
            1e-9
        ));
        for k in 0..3 {
            assert!(exact.box_point[k] >= lo[k] - 1e-12 && exact.box_point[k] <= hi[k] + 1e-12);
        }
        assert!(close(
            exact.distance,
            dist3(exact.segment_point, exact.box_point),
            1e-12
        ));
    }
}

#[test]
fn segment_box_distance_is_zero_when_the_segment_crosses_the_box() {
    let lo = [-1.0; 3];
    let hi = [1.0; 3];
    // Passes through the box without an end point inside it.
    let through = segment_box_distance([-3.0, 0.2, 0.1], [3.0, 0.2, 0.1], lo, hi);
    assert_eq!(through.distance, 0.0);
    assert!(through.segment_point.iter().all(|c| c.abs() <= 1.0));
    // Both end points outside; the segment x + y = 2 only touches the box edge
    // (1, 1, z) at (1, 1, 0).
    let touching = segment_box_distance([2.0, 0.0, 0.0], [0.0, 2.0, 0.0], lo, hi);
    assert!(close(touching.distance, 0.0, 1e-12));
    assert!(vec_close(touching.box_point, [1.0, 1.0, 0.0], 1e-12));
    // The segment x + y = 2.5 passes the edge (1, 1, z) at distance 0.5 / sqrt(2).
    let cutting = segment_box_distance([2.0, 0.5, 0.0], [0.5, 2.0, 0.0], lo, hi);
    assert!(close(cutting.distance, 0.125_f64.sqrt(), 1e-12));
    assert!(vec_close(cutting.segment_point, [1.25, 1.25, 0.0], 1e-12));
    assert!(vec_close(cutting.box_point, [1.0, 1.0, 0.0], 1e-12));
}

#[test]
fn segment_box_distance_degenerate_and_axis_parallel_segments() {
    let lo = [-1.0; 3];
    let hi = [1.0; 3];
    // A segment that is a single point.
    let point = segment_box_distance([3.0, 4.0, 0.0], [3.0, 4.0, 0.0], lo, hi);
    assert!(close(
        point.distance,
        (2.0_f64 * 2.0 + 3.0 * 3.0).sqrt(),
        1e-12
    ));
    assert!(vec_close(point.box_point, [1.0, 1.0, 0.0], 1e-12));
    // Parallel to a face at height 2: distance 1 everywhere along the overlap.
    let parallel = segment_box_distance([-5.0, 2.0, 0.0], [5.0, 2.0, 0.0], lo, hi);
    assert!(close(parallel.distance, 1.0, 1e-12));
    assert!((-1.0..=1.0).contains(&parallel.segment_point[0]));
    // Segment running along an edge line outside the box.
    let edge = segment_box_distance([-5.0, 2.0, 2.0], [5.0, 2.0, 2.0], lo, hi);
    assert!(close(edge.distance, 2.0_f64.sqrt(), 1e-12));
}

// ─────────────────────────────────────────────────────────────────────────────
// Box-box distance
// ─────────────────────────────────────────────────────────────────────────────

#[test]
fn aabb_separation_face_edge_corner_and_overlap() {
    let unit_lo = [-1.0; 3];
    let unit_hi = [1.0; 3];
    // Face to face: separated along y only; x and z overlap in [0, 1].
    let face = aabb_separation([0.0, 3.0, 0.0], [2.0, 5.0, 2.0], unit_lo, unit_hi);
    assert!(close(face.distance, 2.0, 1e-12));
    assert!(
        face.direction
            .is_some_and(|n| vec_close(n, [0.0, 1.0, 0.0], 1e-12))
    );
    assert!(close(face.point_a[1], 3.0, 1e-12) && close(face.point_b[1], 1.0, 1e-12));
    assert!(close(face.point_a[0], 0.5, 1e-12) && close(face.point_b[0], 0.5, 1e-12));
    // Edge to edge: separated along x and y by 3 and 4.
    let edge = aabb_separation([4.0, -8.0, -0.5], [6.0, -5.0, 0.5], unit_lo, unit_hi);
    assert!(close(edge.distance, 5.0, 1e-12));
    assert!(
        edge.direction
            .is_some_and(|n| vec_close(n, [0.6, -0.8, 0.0], 1e-12))
    );
    // Corner to corner.
    let corner = aabb_separation([2.0, 2.0, 3.0], [4.0, 4.0, 4.0], unit_lo, unit_hi);
    assert!(close(corner.distance, (1.0_f64 + 1.0 + 4.0).sqrt(), 1e-12));
    assert!(vec_close(corner.point_a, [2.0, 2.0, 3.0], 1e-12));
    assert!(vec_close(corner.point_b, [1.0, 1.0, 1.0], 1e-12));
    // Overlap: no separation, no direction.
    let overlap = aabb_separation([0.5, 0.5, 0.5], [2.0, 2.0, 2.0], unit_lo, unit_hi);
    assert_eq!(overlap.distance, 0.0);
    assert!(overlap.direction.is_none());
    // Touching faces.
    let touching = aabb_separation([1.0, -1.0, -1.0], [3.0, 1.0, 1.0], unit_lo, unit_hi);
    assert_eq!(touching.distance, 0.0);
}

// ─────────────────────────────────────────────────────────────────────────────
// Box casts
// ─────────────────────────────────────────────────────────────────────────────

#[test]
fn box_vs_box_oblique_approach_hits_at_the_analytic_time() {
    // Box A (half 0.5) lands on the top face (y = 1) of box B (half 1) when its
    // centre is at (0, 1.5, 0), at t = 0.4 for every angle.
    for angle_deg in [60.0, 80.0, 85.0] {
        let motion = oblique_motion([0.0, 1.5, 0.0], angle_deg, 10.0, 0.4);
        let result = BoxCast::vs_box(
            &motion,
            [0.5; 3],
            [0.0; 3],
            [1.0; 3],
            &ShapeCastConfig::default(),
        );
        assert!(result.hit && !result.is_penetrating, "angle {angle_deg}");
        assert!(
            close(result.toi, 0.4, 1e-6),
            "angle {angle_deg}: {}",
            result.toi
        );
        assert!(vec_close(result.normal_at_toi, [0.0, 1.0, 0.0], 1e-6));
        assert!(close(result.witness_a[1], 1.0, 1e-6) && close(result.witness_b[1], 1.0, 1e-6));
        let analytic = slab_toi(motion.pos, [0.5; 3], motion.vel, [0.0; 3], [1.0; 3], 1.0);
        assert!(analytic.is_some_and(|t| close(t, 0.4, 1e-12)));
    }
}

#[test]
fn box_vs_box_corner_contact_reports_the_diagonal_normal() {
    // Head-on corner approach: A's lower corner meets B's upper corner (1, 1, .)
    // when A's centre is at (1.5, 1.5, 0).
    let motion = LinearMotion::new([4.0, 4.0, 0.0], [-3.0, -3.0, 0.0]);
    let result = BoxCast::vs_box(&motion, [0.5; 3], [0.0; 3], [1.0; 3], &default_cfg());
    assert!(result.hit);
    assert!(close(result.toi, 2.5 / 3.0, 1e-6), "toi={}", result.toi);
    let diagonal = 0.5_f64.sqrt();
    assert!(vec_close(
        result.normal_at_toi,
        [diagonal, diagonal, 0.0],
        1e-6
    ));
    assert!(close(result.witness_b[0], 1.0, 1e-6) && close(result.witness_b[1], 1.0, 1e-6));
}

#[test]
fn box_vs_box_initial_overlap_uses_the_axis_of_least_penetration() {
    let motion = LinearMotion::new([0.75, 0.0, 0.0], [0.0, -1.0, 0.0]);
    let result = BoxCast::vs_box(&motion, [0.5; 3], [0.0; 3], [0.5; 3], &default_cfg());
    assert!(result.hit && result.is_penetrating);
    assert_eq!(result.toi, 0.0);
    assert!(vec_close(result.normal_at_toi, [1.0, 0.0, 0.0], 1e-12));
}

#[test]
fn box_vs_box_matches_the_slab_oracle_on_random_inputs() {
    let mut rng = Rng::new(0x5EED_0002);
    let cfg = ShapeCastConfig {
        max_toc: 1.0,
        target_distance: 0.0,
        toi_tolerance: 1e-9,
    };
    let mut hits = 0;
    for case in 0..500 {
        let half_a = [
            rng.range(0.1, 1.0),
            rng.range(0.1, 1.0),
            rng.range(0.1, 1.0),
        ];
        let half_b = [
            rng.range(0.1, 1.0),
            rng.range(0.1, 1.0),
            rng.range(0.1, 1.0),
        ];
        let start = rng.vec(-6.0, 6.0);
        let vel = scale3(rng.vec(-1.0, 1.0), 8.0);
        let center_b = rng.vec(-1.0, 1.0);
        let motion = LinearMotion::new(start, vel);
        let expected = slab_toi(start, half_a, vel, center_b, half_b, 1.0);
        let result = BoxCast::vs_box(&motion, half_a, center_b, half_b, &cfg);
        match expected {
            Some(t) => {
                // Skip contacts that are within the tolerance band of the end of the sweep.
                if t > 1.0 - 1e-6 {
                    continue;
                }
                hits += 1;
                assert!(result.hit, "case {case}: missed a contact at {t}");
                assert!(
                    close(result.toi, t, 1e-6 * (1.0 + t)),
                    "case {case}: toi {} expected {t}",
                    result.toi
                );
            }
            None => {
                // Skip near-grazing no-hit cases.
                let near = slab_toi(start, half_a, vel, center_b, half_b, 1.0 + 1e-6);
                if near.is_some() {
                    continue;
                }
                assert!(!result.hit, "case {case}: false hit at {}", result.toi);
            }
        }
    }
    assert!(
        hits > 30,
        "the random scenes should contain contacts, got {hits}"
    );
}

#[test]
fn box_vs_triangle_lands_on_the_middle_of_a_large_triangle() {
    // A box of half extent 0.1 dropped from y = 5 at 10 units per step onto the
    // centre of a 200 x 200 triangle in the plane y = 0. Only the triangle's
    // far-away vertices were measured before, so the box overshot it.
    let tri = [
        [-100.0, 0.0, -100.0],
        [100.0, 0.0, -100.0],
        [0.0, 0.0, 100.0],
    ];
    let motion = LinearMotion::new([0.0, 5.0, 0.0], [0.0, -10.0, 0.0]);
    let result = BoxCast::vs_triangle(&motion, [0.1; 3], tri[0], tri[1], tri[2], &default_cfg());
    assert!(result.hit && !result.is_penetrating);
    assert!(close(result.toi, 0.49, 1e-9), "toi={}", result.toi);
    assert!(vec_close(result.normal_at_toi, [0.0, 1.0, 0.0], 1e-9));
    // The witnesses are the facing points of the box bottom and of the triangle.
    assert!(close(result.witness_a[1], 0.0, 1e-9));
    assert!(close(result.witness_b[1], 0.0, 1e-9));
    assert!(result.witness_a[0].abs() <= 0.1 + 1e-9 && result.witness_a[2].abs() <= 0.1 + 1e-9);
    let (nearest, _) = crate_point_tri(result.witness_b, tri[0], tri[1], tri[2]);
    assert!(vec_close(nearest, result.witness_b, 1e-9));
}

#[test]
fn box_vs_triangle_oblique_approach_hits_at_the_analytic_time() {
    let tri = [
        [-100.0, 0.0, -100.0],
        [100.0, 0.0, -100.0],
        [0.0, 0.0, 100.0],
    ];
    for angle_deg in [60.0, 80.0, 85.0] {
        // The box centre is 0.1 above the plane at the contact, t = 0.4.
        let motion = oblique_motion([0.0, 0.1, 0.0], angle_deg, 10.0, 0.4);
        let result = BoxCast::vs_triangle(
            &motion,
            [0.1; 3],
            tri[0],
            tri[1],
            tri[2],
            &ShapeCastConfig::default(),
        );
        assert!(result.hit, "angle {angle_deg}");
        assert!(
            close(result.toi, 0.4, 1e-6),
            "angle {angle_deg}: {}",
            result.toi
        );
        assert!(vec_close(result.normal_at_toi, [0.0, 1.0, 0.0], 1e-6));
    }
}

#[test]
fn box_vs_triangle_edge_contact_has_the_in_plane_normal() {
    // Triangle in the plane y = 0 with the hypotenuse x + z = 1. A box of half
    // extent 0.1 slides horizontally towards it: its vertical edge meets the
    // hypotenuse when x_c + z_c = 1.2, i.e. at t = 2.4.
    let (v0, v1, v2) = ([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]);
    let motion = LinearMotion::new([3.0, 0.0, 3.0], [-1.0, 0.0, -1.0]);
    let result = BoxCast::vs_triangle(&motion, [0.1; 3], v0, v1, v2, &cfg_until(5.0));
    assert!(result.hit);
    assert!(close(result.toi, 2.4, 1e-6), "toi={}", result.toi);
    let diagonal = 0.5_f64.sqrt();
    assert!(vec_close(
        result.normal_at_toi,
        [diagonal, 0.0, diagonal],
        1e-6
    ));
    assert!(close(result.witness_b[1], 0.0, 1e-9));
    assert!(close(result.witness_b[0] + result.witness_b[2], 1.0, 1e-6));
}

#[test]
fn box_vs_triangle_passes_beside_and_above_without_a_hit() {
    let tri = [[-1.0, 0.0, -1.0], [1.0, 0.0, -1.0], [0.0, 0.0, 1.0]];
    // Falls well outside the triangle's footprint.
    let beside = LinearMotion::new([5.0, 5.0, 0.0], [0.0, -10.0, 0.0]);
    assert!(!BoxCast::vs_triangle(&beside, [0.1; 3], tri[0], tri[1], tri[2], &default_cfg()).hit);
    // Moving away from the triangle.
    let away = LinearMotion::new([0.0, 5.0, 0.0], [0.0, 10.0, 0.0]);
    assert!(!BoxCast::vs_triangle(&away, [0.1; 3], tri[0], tri[1], tri[2], &default_cfg()).hit);
    // Too slow to arrive within max_toc.
    let slow = LinearMotion::new([0.0, 5.0, 0.0], [0.0, -1.0, 0.0]);
    assert!(!BoxCast::vs_triangle(&slow, [0.1; 3], tri[0], tri[1], tri[2], &default_cfg()).hit);
}

#[test]
fn box_vs_triangle_matches_the_separating_axis_oracle_on_random_inputs() {
    let mut rng = Rng::new(0x5EED_0003);
    let cfg = ShapeCastConfig {
        max_toc: 1.0,
        target_distance: 0.0,
        toi_tolerance: 1e-9,
    };
    let mut hits = 0;
    let mut misses = 0;
    for case in 0..600 {
        let half = [
            rng.range(0.05, 0.5),
            rng.range(0.05, 0.5),
            rng.range(0.05, 0.5),
        ];
        let tri = [rng.vec(-2.0, 2.0), rng.vec(-2.0, 2.0), rng.vec(-2.0, 2.0)];
        // Start 3 to 6 units away and sweep roughly towards the origin.
        let direction = unit_vector(rng.vec(-1.0, 1.0)).unwrap_or([0.0, 1.0, 0.0]);
        let start = scale3(direction, rng.range(3.0, 6.0));
        let aim = scale3(rng.vec(-1.0, 1.0), 1.5);
        let vel = scale3(sub3(aim, start), rng.range(0.6, 1.4));
        let motion = LinearMotion::new(start, vel);
        let result = BoxCast::vs_triangle(&motion, half, tri[0], tri[1], tri[2], &cfg);
        let Some((t_in, t_out)) = sat_box_triangle_interval(start, half, vel, tri) else {
            assert!(!result.hit, "case {case}: false hit at {}", result.toi);
            misses += 1;
            continue;
        };
        // Skip grazing contacts and contacts at the very ends of the sweep window,
        // where the exact oracle and the tolerance-based cast may legitimately differ.
        let margin = 1e-5;
        if (t_in - t_out).abs() < margin || (t_in - 1.0).abs() < margin || t_out.abs() < margin {
            continue;
        }
        if t_in <= t_out && t_out >= 0.0 && t_in <= 1.0 {
            let expected = t_in.max(0.0);
            hits += 1;
            assert!(result.hit, "case {case}: missed a contact at {expected}");
            assert!(
                close(result.toi, expected, 1e-6 * (1.0 + expected)),
                "case {case}: toi {} expected {expected}",
                result.toi
            );
        } else {
            misses += 1;
            assert!(!result.hit, "case {case}: false hit at {}", result.toi);
        }
    }
    assert!(hits > 20, "expected random contacts, got {hits}");
    assert!(misses > 20, "expected random misses, got {misses}");
}

// ─────────────────────────────────────────────────────────────────────────────
// Convex hull distance
// ─────────────────────────────────────────────────────────────────────────────

fn assert_valid_witnesses(
    nearest: &HullDistance,
    on_a: impl Fn([f64; 3]) -> f64,
    on_b: impl Fn([f64; 3]) -> f64,
    label: &str,
) {
    assert!(
        on_a(nearest.point_a) <= 1e-9,
        "{label}: witness A off its hull by {}",
        on_a(nearest.point_a)
    );
    assert!(
        on_b(nearest.point_b) <= 1e-9,
        "{label}: witness B off its hull by {}",
        on_b(nearest.point_b)
    );
    assert!(
        close(
            nearest.distance,
            dist3(nearest.point_a, nearest.point_b),
            1e-9 * (1.0 + nearest.distance)
        ),
        "{label}: witnesses do not realise the distance"
    );
}

#[test]
fn hull_distance_of_boxes_matches_the_closed_form() {
    let mut rng = Rng::new(0x5EED_0004);
    let mut separated = 0;
    let mut overlapping = 0;
    for case in 0..500 {
        let half_a = [
            rng.range(0.1, 1.0),
            rng.range(0.1, 1.0),
            rng.range(0.1, 1.0),
        ];
        let half_b = [
            rng.range(0.1, 1.0),
            rng.range(0.1, 1.0),
            rng.range(0.1, 1.0),
        ];
        let center_a = rng.vec(-2.5, 2.5);
        let center_b = rng.vec(-2.5, 2.5);
        let corners_a = box_corners(half_a);
        let corners_b = translated(&box_corners(half_b), center_b);
        let expected = aabb_separation(
            sub3(center_a, half_a),
            add3(center_a, half_a),
            sub3(center_b, half_b),
            add3(center_b, half_b),
        );
        let got = hull_distance(&corners_a, center_a, &corners_b);
        let Some(got) = got else {
            panic!("case {case}: no result for non-empty hulls");
        };
        assert!(
            close(got.distance, expected.distance, 1e-9),
            "case {case}: distance {} expected {}",
            got.distance,
            expected.distance
        );
        if expected.distance > 0.0 {
            separated += 1;
            assert!(!got.intersecting, "case {case}");
        } else {
            overlapping += 1;
            assert!(got.intersecting, "case {case}");
        }
        let (lo_a, hi_a) = (sub3(center_a, half_a), add3(center_a, half_a));
        let (lo_b, hi_b) = (sub3(center_b, half_b), add3(center_b, half_b));
        let outside = |lo: [f64; 3], hi: [f64; 3]| {
            move |p: [f64; 3]| {
                dist3(
                    p,
                    [
                        p[0].clamp(lo[0], hi[0]),
                        p[1].clamp(lo[1], hi[1]),
                        p[2].clamp(lo[2], hi[2]),
                    ],
                )
            }
        };
        assert_valid_witnesses(
            &got,
            outside(lo_a, hi_a),
            outside(lo_b, hi_b),
            &format!("case {case}"),
        );
    }
    assert!(
        separated > 100 && overlapping > 20,
        "{separated} / {overlapping}"
    );
}

#[test]
fn hull_distance_point_against_box_and_degenerate_hulls() {
    let corners = box_corners([1.0; 3]);
    // Point against a box: distance to the clamped point.
    for p in [
        [3.0, 0.5, -0.5],
        [2.0, 2.0, 2.0],
        [0.0, 0.0, 5.0],
        [0.25, 0.25, 0.25],
    ] {
        let q = [
            p[0].clamp(-1.0, 1.0),
            p[1].clamp(-1.0, 1.0),
            p[2].clamp(-1.0, 1.0),
        ];
        let got = hull_distance(&[p], [0.0; 3], &corners);
        assert!(
            got.is_some_and(|g| close(g.distance, dist3(p, q), 1e-12)),
            "point {p:?}"
        );
    }
    // Two single points.
    let points = hull_distance(&[[0.0, 0.0, 0.0]], [3.0, 4.0, 0.0], &[[0.0; 3]]);
    assert!(points.is_some_and(|g| close(g.distance, 5.0, 1e-12) && !g.intersecting));
    // Segment against segment: skew segments one unit apart, then crossing ones.
    let skew = hull_distance(
        &[[-1.0, 0.0, 0.0], [1.0, 0.0, 0.0]],
        [0.0, 0.0, 1.0],
        &[[0.0, -1.0, 0.0], [0.0, 1.0, 0.0]],
    );
    assert!(skew.is_some_and(|g| close(g.distance, 1.0, 1e-12)));
    let crossing = hull_distance(
        &[[-1.0, 0.0, 0.0], [1.0, 0.0, 0.0]],
        [0.0; 3],
        &[[0.0, -1.0, 0.0], [0.0, 1.0, 0.0]],
    );
    assert!(crossing.is_some_and(|g| g.intersecting && g.distance == 0.0));
    // Collinear overlapping segments intersect; collinear disjoint ones do not.
    let overlapping = hull_distance(
        &[[0.0; 3], [2.0, 0.0, 0.0]],
        [1.0, 0.0, 0.0],
        &[[0.0; 3], [2.0, 0.0, 0.0]],
    );
    assert!(overlapping.is_some_and(|g| g.intersecting));
    let disjoint = hull_distance(
        &[[0.0; 3], [1.0, 0.0, 0.0]],
        [3.0, 0.0, 0.0],
        &[[0.0; 3], [1.0, 0.0, 0.0]],
    );
    assert!(disjoint.is_some_and(|g| close(g.distance, 2.0, 1e-12)));
    // A flat square against a point above it and a point beside it.
    let square = [
        [-1.0, 0.0, -1.0],
        [1.0, 0.0, -1.0],
        [1.0, 0.0, 1.0],
        [-1.0, 0.0, 1.0],
    ];
    let above = hull_distance(&[[0.5, 2.0, 0.5]], [0.0; 3], &square);
    assert!(above.is_some_and(|g| close(g.distance, 2.0, 1e-12)));
    let beside = hull_distance(&[[4.0, 3.0, 0.0]], [0.0; 3], &square);
    assert!(beside.is_some_and(|g| close(g.distance, 18.0_f64.sqrt(), 1e-12)));
    // Empty hulls have no distance.
    assert!(hull_distance(&[], [0.0; 3], &corners).is_none());
    assert!(hull_distance(&corners, [0.0; 3], &[]).is_none());
}

#[test]
fn hull_distance_of_intersecting_and_touching_boxes() {
    let corners = box_corners([1.0; 3]);
    let overlapping = hull_distance(&corners, [0.5, 0.25, 0.0], &corners);
    assert!(overlapping.is_some_and(|g| g.intersecting && g.distance == 0.0));
    let touching = hull_distance(&corners, [2.0, 0.0, 0.0], &corners);
    assert!(touching.is_some_and(|g| g.distance < 1e-12));
    let nested = hull_distance(&box_corners([0.25; 3]), [0.0; 3], &corners);
    assert!(nested.is_some_and(|g| g.intersecting));
}

#[test]
fn hull_distance_of_triangles_matches_the_feature_pair_oracle() {
    let mut rng = Rng::new(0x5EED_0005);
    for case in 0..500 {
        // Separated by the plane x = 0, so the feature-pair minimum is the distance.
        let mut side = |lo: f64, hi: f64| -> [[f64; 3]; 3] {
            std::array::from_fn(|_| {
                [
                    rng.range(lo, hi),
                    rng.range(-2.0, 2.0),
                    rng.range(-2.0, 2.0),
                ]
            })
        };
        let a = side(-3.0, -0.05);
        let b = side(0.05, 3.0);
        let expected = triangle_triangle_distance(a, b);
        let got = hull_distance(&a, [0.0; 3], &b);
        let Some(got) = got else {
            panic!("case {case}: no result");
        };
        assert!(
            close(got.distance, expected, 1e-9),
            "case {case}: distance {} expected {expected}",
            got.distance
        );
        assert!(!got.intersecting, "case {case}");
        let on_triangle = |t: [[f64; 3]; 3]| {
            move |p: [f64; 3]| {
                let (q, _) = crate_point_tri(p, t[0], t[1], t[2]);
                dist3(p, q)
            }
        };
        assert_valid_witnesses(
            &got,
            on_triangle(a),
            on_triangle(b),
            &format!("case {case}"),
        );
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Convex casts
// ─────────────────────────────────────────────────────────────────────────────

#[test]
fn convex_cast_cube_stops_on_the_face_of_a_large_cube() {
    // A unit cube falls from y = 10 onto the top face (y = 5) of a 10 x 10 x 10
    // cube. The closest *vertex* pair is about 7.8 units apart, which made the
    // vertex-pair "GJK" step far past the face; the true gap is 4.5 units.
    let shape_a = box_corners([0.5; 3]);
    let shape_b = box_corners([5.0; 3]);
    let motion = LinearMotion::new([0.0, 10.0, 0.0], [0.0, -10.0, 0.0]);
    let result = ConvexCast::cast(&motion, &shape_a, &shape_b, &default_cfg());
    assert!(result.hit && !result.is_penetrating);
    assert!(close(result.toi, 0.45, 1e-9), "toi={}", result.toi);
    assert!(vec_close(result.normal_at_toi, [0.0, 1.0, 0.0], 1e-9));
    assert!(close(result.witness_a[1], 5.0, 1e-9) && close(result.witness_b[1], 5.0, 1e-9));
    // The cube stops on the face plane: its centre is half an edge above it.
    let stop = motion.position_at(result.toi);
    assert!(close(stop[1], 5.5, 1e-8));
}

#[test]
fn convex_cast_with_empty_shapes_is_a_miss_instead_of_a_panic() {
    let cube = box_corners([0.5; 3]);
    let motion = LinearMotion::new([0.0, 5.0, 0.0], [0.0, -1.0, 0.0]);
    assert!(!ConvexCast::cast(&motion, &[], &cube, &default_cfg()).hit);
    assert!(!ConvexCast::cast(&motion, &cube, &[], &default_cfg()).hit);
    assert!(!ConvexCast::cast(&motion, &[], &[], &default_cfg()).hit);
}

#[test]
fn convex_cast_oblique_approach_hits_at_the_analytic_time() {
    for angle_deg in [60.0, 80.0, 85.0] {
        let motion = oblique_motion([0.0, 1.5, 0.0], angle_deg, 10.0, 0.4);
        let result = ConvexCast::cast(
            &motion,
            &box_corners([0.5; 3]),
            &box_corners([1.0; 3]),
            &ShapeCastConfig::default(),
        );
        assert!(result.hit, "angle {angle_deg}");
        assert!(
            close(result.toi, 0.4, 1e-6),
            "angle {angle_deg}: {}",
            result.toi
        );
        assert!(vec_close(result.normal_at_toi, [0.0, 1.0, 0.0], 1e-6));
    }
}

#[test]
fn convex_cast_curved_gap_with_a_target_distance_matches_the_quadratic_root() {
    // Point against point: the gap |p0 + t v| is a curved (non-planar) function of
    // time. With a target separation of 0.5 the contact is the smaller root of
    // (t - 3)^2 + 0.3^2 = 0.5^2, i.e. t = 3 - 0.4 = 2.6, with normal (-0.8, 0.6, 0).
    let motion = LinearMotion::new([-3.0, 0.3, 0.0], [1.0, 0.0, 0.0]);
    let cfg = ShapeCastConfig {
        max_toc: 10.0,
        target_distance: 0.5,
        toi_tolerance: 1e-12,
    };
    let result = ConvexCast::cast(&motion, &[[0.0; 3]], &[[0.0; 3]], &cfg);
    assert!(result.hit && !result.is_penetrating);
    assert!(close(result.toi, 2.6, 1e-9), "toi={}", result.toi);
    assert!(vec_close(result.normal_at_toi, [-0.8, 0.6, 0.0], 1e-9));
    assert!(vec_close(result.witness_a, [-0.4, 0.3, 0.0], 1e-9));
    assert!(vec_close(result.witness_b, [0.0; 3], 1e-12));
    // Passing at 0.51 > 0.5 misses, passing at 0.49 < 0.5 hits.
    let wide = LinearMotion::new([-5.0, 0.51, 0.0], [1.0, 0.0, 0.0]);
    assert!(!ConvexCast::cast(&wide, &[[0.0; 3]], &[[0.0; 3]], &cfg).hit);
    let narrow = LinearMotion::new([-5.0, 0.49, 0.0], [1.0, 0.0, 0.0]);
    let grazing = ConvexCast::cast(&narrow, &[[0.0; 3]], &[[0.0; 3]], &cfg);
    assert!(grazing.hit);
    let expected = 5.0 - (0.25_f64 - 0.49 * 0.49).sqrt();
    assert!(
        close(grazing.toi, expected, 1e-6),
        "toi={} expected {expected}",
        grazing.toi
    );
}

#[test]
fn convex_cast_ignores_separating_and_parallel_motion() {
    let cube = box_corners([0.5; 3]);
    let away = LinearMotion::new([0.0, 3.0, 0.0], [0.0, 1.0, 0.0]);
    assert!(!ConvexCast::cast(&away, &cube, &cube, &cfg_until(100.0)).hit);
    let sliding = LinearMotion::new([0.0, 3.0, 0.0], [1.0, 0.0, 0.0]);
    assert!(!ConvexCast::cast(&sliding, &cube, &cube, &cfg_until(100.0)).hit);
    let still = LinearMotion::new([0.0, 3.0, 0.0], [0.0; 3]);
    assert!(!ConvexCast::cast(&still, &cube, &cube, &cfg_until(100.0)).hit);
}

#[test]
fn convex_cast_initial_overlap_is_penetrating_and_uses_the_reverse_motion() {
    let cube = box_corners([0.5; 3]);
    let motion = LinearMotion::new([0.25, 0.0, 0.0], [0.0, -2.0, 0.0]);
    let result = ConvexCast::cast(&motion, &cube, &cube, &default_cfg());
    assert!(result.hit && result.is_penetrating);
    assert_eq!(result.toi, 0.0);
    assert!(vec_close(result.normal_at_toi, [0.0, 1.0, 0.0], 1e-12));
}

#[test]
fn convex_cast_matches_the_slab_oracle_on_random_boxes() {
    let mut rng = Rng::new(0x5EED_0006);
    let cfg = ShapeCastConfig {
        max_toc: 1.0,
        target_distance: 0.0,
        toi_tolerance: 1e-9,
    };
    let mut hits = 0;
    for case in 0..300 {
        let half_a = [
            rng.range(0.1, 1.0),
            rng.range(0.1, 1.0),
            rng.range(0.1, 1.0),
        ];
        let half_b = [
            rng.range(0.1, 1.0),
            rng.range(0.1, 1.0),
            rng.range(0.1, 1.0),
        ];
        let start = rng.vec(-6.0, 6.0);
        let vel = scale3(rng.vec(-1.0, 1.0), 8.0);
        let center_b = rng.vec(-1.0, 1.0);
        let motion = LinearMotion::new(start, vel);
        let expected = slab_toi(start, half_a, vel, center_b, half_b, 1.0);
        let result = ConvexCast::cast(
            &motion,
            &box_corners(half_a),
            &translated(&box_corners(half_b), center_b),
            &cfg,
        );
        match expected {
            Some(t) if t <= 1.0 - 1e-6 => {
                hits += 1;
                assert!(result.hit, "case {case}: missed a contact at {t}");
                assert!(
                    close(result.toi, t, 1e-6 * (1.0 + t)),
                    "case {case}: toi {} expected {t}",
                    result.toi
                );
            }
            Some(_) => {}
            None => {
                if slab_toi(start, half_a, vel, center_b, half_b, 1.0 + 1e-6).is_none() {
                    assert!(!result.hit, "case {case}: false hit at {}", result.toi);
                }
            }
        }
    }
    assert!(hits > 20, "expected random contacts, got {hits}");
}

// ─────────────────────────────────────────────────────────────────────────────
// Shape cast configuration
// ─────────────────────────────────────────────────────────────────────────────

#[test]
fn target_distance_stops_the_cast_early_by_the_requested_margin() {
    let motion = LinearMotion::new([0.0, 5.0, 0.0], [0.0, -10.0, 0.0]);
    let base = BoxCast::vs_box(&motion, [0.5; 3], [0.0; 3], [1.0; 3], &default_cfg());
    // Contact when the centre is at y = 1.5: t = 0.35.
    assert!(close(base.toi, 0.35, 1e-9), "toi={}", base.toi);
    let skin = ShapeCastConfig {
        target_distance: 0.2,
        ..ShapeCastConfig::default()
    };
    let early = BoxCast::vs_box(&motion, [0.5; 3], [0.0; 3], [1.0; 3], &skin);
    assert!(close(early.toi, 0.33, 1e-6), "toi={}", early.toi);
    assert!(close(early.witness_a[1] - early.witness_b[1], 0.2, 1e-6));
}
