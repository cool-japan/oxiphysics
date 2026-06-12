// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! One-shot full box-box contact manifold generation.
//!
//! Generates up to four contact points for a box-box overlap in a single pass
//! using separating-axis selection plus Sutherland-Hodgman face clipping. The
//! design goal is the *no-rock guarantee*: a box resting flat on the ground
//! produces four symmetric corner contacts with a single shared face normal, so
//! a sequential-impulse solver applies zero net torque and the box does not rock.
//!
//! All returned contacts carry a normal pointing FROM body B TOWARD body A
//! (the crate convention). Contact points straddle the interface symmetrically
//! (`point_a = m + n_ba * d/2`, `point_b = m - n_ba * d/2`) so swapping the
//! dispatch order (A<->B) negates the normal and swaps the points, leaving the
//! manifold geometrically identical.

use oxiphysics_core::Transform;
use oxiphysics_core::math::Vec3;
use oxiphysics_geometry::BoxShape;

use super::convex_manifold::{Polyhedron, polyhedron_manifold_from_normal};
use super::types::NarrowPhaseResult;
use crate::contact_generation::{add3, cross3, dot3, norm3, normalize3, scale3, sub3};
use crate::types::{CollisionPair, Contact, ContactManifold};

/// Reciprocal of the spatial grid used to quantise contact sort keys. The grid
/// (1e-4 m) is far above floating-point noise yet far below the separation of
/// distinct box features, so the contact ordering is stable under tiny nudges.
const SORT_QUANTUM_INV: f64 = 1.0e4;

/// Build the three world-space axes (columns of the rotation matrix) of a body.
///
/// Column `i` is the world direction of the body's local axis `i`. Returned
/// explicitly (rather than via an indexed loop) to keep clippy free of the
/// `needless_range_loop` lint.
fn world_axes(t: &Transform) -> [[f64; 3]; 3] {
    let r = t.rotation.to_rotation_matrix();
    let m = r.matrix();
    [
        [m[(0, 0)], m[(1, 0)], m[(2, 0)]],
        [m[(0, 1)], m[(1, 1)], m[(2, 1)]],
        [m[(0, 2)], m[(1, 2)], m[(2, 2)]],
    ]
}

/// Half-width of a box projected onto a (not necessarily unit) world direction `l`.
///
/// `half[i]` is the half-extent along local axis `i`, `axes[i]` the world
/// direction of that axis. The projection radius is `sum_i half[i] * |axes[i] . l|`.
fn projection_radius(half: [f64; 3], axes: [[f64; 3]; 3], l: [f64; 3]) -> f64 {
    half.iter()
        .zip(axes.iter())
        .map(|(&h, ax)| h * dot3(*ax, l).abs())
        .sum()
}

/// Which separating axis won the SAT minimisation.
///
/// The face variants are payload-free: the face-face manifold is built by the
/// shared convex helper, which re-selects the reference/incident faces from the
/// world normals, so the SAT axis index is no longer needed downstream. The
/// edge variant keeps both local edge-axis indices because the edge-edge witness
/// path restricts its segment search to those parallel edges.
enum BoxSatAxis {
    /// Face normal of box A won the minimisation.
    FaceA,
    /// Face normal of box B won the minimisation.
    FaceB,
    /// Cross product of A's local edge `i` and B's local edge `j`.
    EdgeEdge(usize, usize),
}

/// Result of the detailed box-box separating-axis test.
struct BoxSatResult {
    /// Penetration depth along `normal` (clamped non-negative).
    depth: f64,
    /// Unit minimum-translation normal, pointing FROM B TOWARD A.
    normal: [f64; 3],
    /// The axis that produced the minimum overlap.
    axis: BoxSatAxis,
}

/// Full separating-axis test for two oriented boxes in world space.
///
/// Uses projection-radius SAT (not the single-point `box_box_sat`). The 15
/// candidate axes are the 3 face normals of A, the 3 of B, and the 9 edge-edge
/// cross products. Face axes are evaluated first; an edge axis only wins if it
/// beats the best face overlap by more than `edge_bias`, which prevents a
/// resting box from selecting a numerically-equal edge axis and rocking.
fn box_box_sat_detailed(
    center_a: [f64; 3],
    axes_a: [[f64; 3]; 3],
    half_a: [f64; 3],
    center_b: [f64; 3],
    axes_b: [[f64; 3]; 3],
    half_b: [f64; 3],
) -> Option<BoxSatResult> {
    let center_diff = sub3(center_a, center_b);

    // Overlap along a unit axis `l`; `None` if the boxes are separated on it.
    let overlap_on = |l: [f64; 3]| -> Option<f64> {
        let radius_a = projection_radius(half_a, axes_a, l);
        let radius_b = projection_radius(half_b, axes_b, l);
        let center_dist = dot3(center_diff, l);
        let separation = center_dist.abs() - (radius_a + radius_b);
        if separation > 0.0 {
            None
        } else {
            Some(-separation)
        }
    };

    let mut min_overlap = f64::INFINITY;
    let mut best_axis = BoxSatAxis::FaceA;
    let mut best_dir = [0.0_f64, 0.0, 0.0];

    // Face axes of A first.
    for ax in axes_a.iter() {
        match overlap_on(*ax) {
            None => return None,
            Some(overlap) => {
                if overlap < min_overlap {
                    min_overlap = overlap;
                    best_axis = BoxSatAxis::FaceA;
                    best_dir = *ax;
                }
            }
        }
    }

    // Face axes of B next.
    for ax in axes_b.iter() {
        match overlap_on(*ax) {
            None => return None,
            Some(overlap) => {
                if overlap < min_overlap {
                    min_overlap = overlap;
                    best_axis = BoxSatAxis::FaceB;
                    best_dir = *ax;
                }
            }
        }
    }

    // Edge-bias: a resting box must prefer the face axis over a numerically
    // equal edge axis. Without this margin, floating-point ties on a flat
    // resting contact can hand the manifold to a spurious edge axis, which
    // collapses the four corner contacts to one and lets the box rock.
    let scale = half_a
        .iter()
        .chain(half_b.iter())
        .cloned()
        .fold(0.0_f64, f64::max)
        .max(1.0);
    let edge_bias = 1e-4 * scale;

    // Edge-edge axes: cross of every A edge direction with every B edge direction.
    for (i, ax_a) in axes_a.iter().enumerate() {
        for (j, ax_b) in axes_b.iter().enumerate() {
            let axis = cross3(*ax_a, *ax_b);
            // Parallel edges produce a degenerate (near-zero) cross; skip them.
            if norm3(axis) < 1e-6 {
                continue;
            }
            let l = normalize3(axis);
            match overlap_on(l) {
                None => return None,
                Some(overlap) => {
                    if overlap < min_overlap - edge_bias {
                        min_overlap = overlap;
                        best_axis = BoxSatAxis::EdgeEdge(i, j);
                        best_dir = l;
                    }
                }
            }
        }
    }

    // Orient the winning normal so it points FROM B TOWARD A.
    let mut n = best_dir;
    if dot3(n, center_diff) < 0.0 {
        n = [-n[0], -n[1], -n[2]];
    }
    let n = normalize3(n);

    Some(BoxSatResult {
        depth: min_overlap.max(0.0),
        normal: n,
        axis: best_axis,
    })
}

/// World-space endpoints of the box edges parallel to local axis `axis`.
///
/// A box has exactly four edges along each local axis; an edge of `edge_list`
/// runs parallel to axis `k` when its two local vertices differ only in
/// coordinate `k`. The endpoints are transformed into world space.
fn world_edges_along_axis(
    verts: &[[f64; 3]; 8],
    transform: &Transform,
    axis: usize,
) -> Vec<([f64; 3], [f64; 3])> {
    let mut edges = Vec::with_capacity(4);
    for (i0, i1) in BoxShape::edge_list() {
        let local_dir = sub3(verts[i1], verts[i0]);
        // The dominant non-zero component identifies the parallel local axis.
        let mut dominant = 0;
        let mut best = local_dir[0].abs();
        for (k, &component) in local_dir.iter().enumerate().skip(1) {
            if component.abs() > best {
                best = component.abs();
                dominant = k;
            }
        }
        if dominant == axis {
            let p0 = transform.transform_point(&Vec3::from(verts[i0]));
            let p1 = transform.transform_point(&Vec3::from(verts[i1]));
            edges.push(([p0.x, p0.y, p0.z], [p1.x, p1.y, p1.z]));
        }
    }
    edges
}

/// Closest points between two segments `[p1,q1]` and `[p2,q2]` (Ericson).
///
/// Returns `(closest_on_seg1, closest_on_seg2, distance)`. Handles degenerate
/// (point-like) segments by clamping the parametric coordinates to `[0,1]`.
fn closest_segment_segment(
    p1: [f64; 3],
    q1: [f64; 3],
    p2: [f64; 3],
    q2: [f64; 3],
) -> ([f64; 3], [f64; 3], f64) {
    let d1 = sub3(q1, p1);
    let d2 = sub3(q2, p2);
    let r = sub3(p1, p2);
    let a = dot3(d1, d1);
    let e = dot3(d2, d2);
    let f = dot3(d2, r);
    let eps = 1e-12;

    let s;
    let t;

    if a <= eps && e <= eps {
        // Both segments degenerate to points.
        s = 0.0;
        t = 0.0;
    } else if a <= eps {
        // First segment degenerate.
        s = 0.0;
        t = (f / e).clamp(0.0, 1.0);
    } else {
        let c = dot3(d1, r);
        if e <= eps {
            // Second segment degenerate.
            t = 0.0;
            s = (-c / a).clamp(0.0, 1.0);
        } else {
            let b = dot3(d1, d2);
            let denom = a * e - b * b;
            let s0 = if denom != 0.0 {
                ((b * f - c * e) / denom).clamp(0.0, 1.0)
            } else {
                0.0
            };
            let t0 = (b * s0 + f) / e;
            if t0 < 0.0 {
                t = 0.0;
                s = (-c / a).clamp(0.0, 1.0);
            } else if t0 > 1.0 {
                t = 1.0;
                s = ((b - c) / a).clamp(0.0, 1.0);
            } else {
                s = s0;
                t = t0;
            }
        }
    }

    let c1 = add3(p1, scale3(d1, s));
    let c2 = add3(p2, scale3(d2, t));
    let dist = norm3(sub3(c1, c2));
    (c1, c2, dist)
}

/// Reduce a contact set to at most four representative points.
///
/// Selection is deterministic with lowest-index tie-breaking:
/// (1) deepest point, (2) farthest from #1, (3) maximum-area triangle with
/// {#1,#2}, (4) maximises area against triangle {#1,#2,#3}. With four or fewer
/// inputs the set is returned unchanged.
fn reduce_contacts_to_4(contacts: Vec<Contact>) -> Vec<Contact> {
    if contacts.len() <= 4 {
        return contacts;
    }

    let pts: Vec<Vec3> = contacts.iter().map(|c| c.point()).collect();
    let n = contacts.len();

    // (1) Deepest contact (lowest index breaks ties).
    let mut i0 = 0;
    let mut best_depth = f64::NEG_INFINITY;
    for (i, c) in contacts.iter().enumerate() {
        if c.depth > best_depth {
            best_depth = c.depth;
            i0 = i;
        }
    }

    // (2) Farthest from #1.
    let mut i1 = i0;
    let mut best_d = f64::NEG_INFINITY;
    for i in 0..n {
        if i == i0 {
            continue;
        }
        let d = (pts[i] - pts[i0]).norm();
        if d > best_d {
            best_d = d;
            i1 = i;
        }
    }

    // (3) Largest triangle area with {#1, #2}.
    let edge = pts[i1] - pts[i0];
    let mut i2 = i0;
    let mut best_area = f64::NEG_INFINITY;
    for i in 0..n {
        if i == i0 || i == i1 {
            continue;
        }
        let area = (edge.cross(&(pts[i] - pts[i0]))).norm();
        if area > best_area {
            best_area = area;
            i2 = i;
        }
    }

    // (4) Maximise total triangle area against {#1, #2, #3}.
    let mut i3 = i0;
    let mut best_area3 = f64::NEG_INFINITY;
    for i in 0..n {
        if i == i0 || i == i1 || i == i2 {
            continue;
        }
        let p = pts[i];
        let a012 = (pts[i1] - pts[i0]).cross(&(p - pts[i0])).norm();
        let a123 = (pts[i2] - pts[i1]).cross(&(p - pts[i1])).norm();
        let a203 = (pts[i0] - pts[i2]).cross(&(p - pts[i2])).norm();
        let total = a012 + a123 + a203;
        if total > best_area3 {
            best_area3 = total;
            i3 = i;
        }
    }

    vec![
        contacts[i0].clone(),
        contacts[i1].clone(),
        contacts[i2].clone(),
        contacts[i3].clone(),
    ]
}

/// One-shot full box-box contact manifold (up to 4 points). Normal points B->A.
pub fn box_box_manifold(
    box_a: &BoxShape,
    transform_a: &Transform,
    box_b: &BoxShape,
    transform_b: &Transform,
    pair: CollisionPair,
) -> NarrowPhaseResult {
    let center_a = [
        transform_a.position.x,
        transform_a.position.y,
        transform_a.position.z,
    ];
    let axes_a = world_axes(transform_a);
    let half_a = [
        box_a.half_extents.x,
        box_a.half_extents.y,
        box_a.half_extents.z,
    ];

    let center_b = [
        transform_b.position.x,
        transform_b.position.y,
        transform_b.position.z,
    ];
    let axes_b = world_axes(transform_b);
    let half_b = [
        box_b.half_extents.x,
        box_b.half_extents.y,
        box_b.half_extents.z,
    ];

    let sat = match box_box_sat_detailed(center_a, axes_a, half_a, center_b, axes_b, half_b) {
        Some(s) => s,
        None => return NarrowPhaseResult::separated(),
    };

    let n_ba = Vec3::from(sat.normal);
    let mut contacts: Vec<Contact> = Vec::new();

    match sat.axis {
        BoxSatAxis::FaceA | BoxSatAxis::FaceB => {
            // Face-face overlap: delegate to the shared convex face-clip helper.
            // Two oriented boxes become world-space polygon polyhedra and the
            // generic Sutherland-Hodgman reference/incident clip produces the
            // contacts — the same algorithm the box-specific path used to inline,
            // now shared with the general convex-convex manifold so both stay in
            // lock-step. The helper picks ref/incident faces by `most_parallel`
            // dot, which for a box agrees with the SAT-axis face selection, and
            // applies the same reduce-to-4 + stable spatial sort.
            let poly_a = Polyhedron::from_box(box_a, transform_a);
            let poly_b = Polyhedron::from_box(box_b, transform_b);
            return NarrowPhaseResult::contact(polyhedron_manifold_from_normal(
                &poly_a,
                &poly_b,
                sat.normal,
                sat.depth,
                pair,
            ));
        }
        BoxSatAxis::EdgeEdge(edge_axis_a, edge_axis_b) => {
            let verts_a = box_a.vertex_list();
            let verts_b = box_b.vertex_list();

            // The SAT minimum came from cross(A.axis[edge_axis_a],
            // B.axis[edge_axis_b]); restrict the search to the four edges of each
            // box that run parallel to those local axes (reads both stored axis
            // indices and avoids the full 12x12 sweep).
            let edges_a = world_edges_along_axis(&verts_a, transform_a, edge_axis_a);
            let edges_b = world_edges_along_axis(&verts_b, transform_b, edge_axis_b);

            let mut best_c1 = [0.0_f64, 0.0, 0.0];
            let mut best_c2 = [0.0_f64, 0.0, 0.0];
            let mut best_dist = f64::INFINITY;
            let mut best_aligned_dist = f64::INFINITY;
            let mut have_aligned = false;
            let mut have_any = false;

            for (a0, a1) in &edges_a {
                for (b0, b1) in &edges_b {
                    let (c1, c2, dist) = closest_segment_segment(*a0, *a1, *b0, *b1);
                    let dir = sub3(c1, c2);
                    let alignment = if norm3(dir) > 1e-9 {
                        dot3(normalize3(dir), sat.normal).abs()
                    } else {
                        0.0
                    };

                    // Globally closest pair (fallback if nothing aligns well).
                    if dist < best_dist {
                        best_dist = dist;
                        if !have_aligned {
                            best_c1 = c1;
                            best_c2 = c2;
                        }
                        have_any = true;
                    }

                    // Prefer pairs whose closest-approach direction matches the
                    // SAT normal; among those keep the closest.
                    if alignment > 0.5 && dist < best_aligned_dist {
                        best_aligned_dist = dist;
                        best_c1 = c1;
                        best_c2 = c2;
                        have_aligned = true;
                    }
                }
            }

            if have_aligned || have_any {
                // c1 lies on A's edge -> point_a; c2 on B's edge -> point_b.
                let pa = Vec3::from(best_c1);
                let pb = Vec3::from(best_c2);
                let depth = sat.depth.max(0.0);
                contacts.push(Contact::new(pa, pb, n_ba, depth));
            }
        }
    }

    let contacts = reduce_contacts_to_4(contacts);
    let mut contacts = contacts;
    // feature-id continuity: plain `Contact` carries no feature id; warm-start
    // continuity is provided by this stable spatial ordering plus the manifold
    // cache's local-position matching. (RichContact's feature_a/feature_b is a
    // future upgrade path for explicit feature-id continuity.)
    //
    // The keys are quantised to a coarse grid before comparison so that
    // sub-grid floating-point noise (e.g. two co-linear corners whose clipped
    // x differs only in the last bits) cannot flip the ordering, while a real
    // positional change (orders of magnitude larger than the grid) still sorts
    // deterministically. This keeps the i-th contact identity stable from frame
    // to frame even under a tiny nudge. See test_feature_id_stability_under_nudge.
    let q = |v: f64| (v * SORT_QUANTUM_INV).round();
    contacts.sort_by(|c1, c2| {
        q(c1.point_b.x)
            .partial_cmp(&q(c2.point_b.x))
            .unwrap_or(std::cmp::Ordering::Equal)
            .then(
                q(c1.point_b.z)
                    .partial_cmp(&q(c2.point_b.z))
                    .unwrap_or(std::cmp::Ordering::Equal),
            )
            .then(
                q(c1.point_b.y)
                    .partial_cmp(&q(c2.point_b.y))
                    .unwrap_or(std::cmp::Ordering::Equal),
            )
    });

    let mut manifold = ContactManifold::new(pair);
    for c in contacts {
        manifold.add_contact(c);
    }
    NarrowPhaseResult::contact(manifold)
}

#[cfg(test)]
mod tests {
    use super::*;
    use oxiphysics_core::Transform;
    use std::f64::consts::FRAC_PI_4;

    const IDENTITY_AXES: [[f64; 3]; 3] = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];

    #[test]
    fn test_sat_axisymmetric_box_on_box_picks_face_y() {
        let half = [0.5, 0.5, 0.5];
        let center_b = [0.0, 0.0, 0.0];
        let center_a = [0.0, 0.99, 0.0];
        let sat =
            box_box_sat_detailed(center_a, IDENTITY_AXES, half, center_b, IDENTITY_AXES, half)
                .expect("overlap");
        // A face axis must win (not an edge axis); the +y normal confirms it is
        // specifically the y face that produced the minimum overlap.
        assert!(matches!(sat.axis, BoxSatAxis::FaceA | BoxSatAxis::FaceB));
        assert!(sat.normal[1] > 0.9);
        assert!((sat.depth - 0.01).abs() < 1e-3);
    }

    #[test]
    fn test_sat_separated_returns_none() {
        let half = [0.5, 0.5, 0.5];
        let center_a = [0.0, 100.0, 0.0];
        let center_b = [0.0, 0.0, 0.0];
        let result =
            box_box_sat_detailed(center_a, IDENTITY_AXES, half, center_b, IDENTITY_AXES, half);
        assert!(result.is_none());
    }

    #[test]
    fn test_sat_edge_axis_selected() {
        // Box A axis-aligned at the origin; box B rotated 45 deg about the body
        // diagonal (1,1,1) and pushed out along (0.75,0.75,0.75) so the closest
        // features are an edge of A and an edge of B. The minimum separating axis
        // is then a genuine edge-edge cross (verified empirically as EdgeEdge).
        let half = [0.5, 0.5, 0.5];
        let ta = Transform::from_position(Vec3::new(0.0, 0.0, 0.0));
        let tb = Transform::from_axis_angle(
            Vec3::new(0.75, 0.75, 0.75),
            Vec3::new(1.0, 1.0, 1.0),
            FRAC_PI_4,
        );
        let center_a = [ta.position.x, ta.position.y, ta.position.z];
        let center_b = [tb.position.x, tb.position.y, tb.position.z];
        let axes_a = world_axes(&ta);
        let axes_b = world_axes(&tb);
        let sat =
            box_box_sat_detailed(center_a, axes_a, half, center_b, axes_b, half).expect("overlap");
        // STRICT path: the diagonal corner-to-edge configuration selects an
        // edge-edge axis with a unit normal and positive depth.
        assert!(matches!(sat.axis, BoxSatAxis::EdgeEdge(_, _)));
        assert!(sat.depth > 0.0);
        assert!((norm3(sat.normal) - 1.0).abs() < 1e-6);
    }

    #[test]
    fn test_box_on_ground_4_points() {
        let box_a = BoxShape::new(Vec3::new(0.5, 0.5, 0.5));
        let box_b = BoxShape::new(Vec3::new(10.0, 0.5, 10.0));
        let tb = Transform::from_position(Vec3::new(0.0, -0.5, 0.0));
        let ta = Transform::from_position(Vec3::new(0.0, 0.5 - 0.005, 0.0));

        let result = box_box_manifold(&box_a, &ta, &box_b, &tb, CollisionPair::new(0, 1));
        let manifold = result.manifold.expect("overlapping -> manifold");
        assert_eq!(manifold.contacts.len(), 4);

        let corners = [
            [-0.5, 0.0, -0.5],
            [0.5, 0.0, -0.5],
            [0.5, 0.0, 0.5],
            [-0.5, 0.0, 0.5],
        ];
        for corner in &corners {
            let found = manifold.contacts.iter().any(|c| {
                (c.point_b.x - corner[0]).abs() < 1e-2
                    && (c.point_b.y - corner[1]).abs() < 1e-2
                    && (c.point_b.z - corner[2]).abs() < 1e-2
            });
            assert!(found, "expected a contact near corner {corner:?}");
        }

        for c in &manifold.contacts {
            assert!(c.normal.y > 0.999);
            assert!(c.normal.x.abs() < 1e-3);
            assert!(c.normal.z.abs() < 1e-3);
            assert!((c.depth - 0.005).abs() < 1e-3);
        }
    }

    #[test]
    fn test_box_on_ground_no_rocking() {
        let box_a = BoxShape::new(Vec3::new(0.5, 0.5, 0.5));
        let box_b = BoxShape::new(Vec3::new(10.0, 0.5, 10.0));
        let tb = Transform::from_position(Vec3::new(0.0, -0.5, 0.0));
        let ta = Transform::from_position(Vec3::new(0.0, 0.5 - 0.005, 0.0));

        let result = box_box_manifold(&box_a, &ta, &box_b, &tb, CollisionPair::new(0, 1));
        let manifold = result.manifold.expect("overlapping -> manifold");
        assert_eq!(manifold.contacts.len(), 4);

        let cross3 = |a: [f64; 3], b: [f64; 3]| {
            [
                a[1] * b[2] - a[2] * b[1],
                a[2] * b[0] - a[0] * b[2],
                a[0] * b[1] - a[1] * b[0],
            ]
        };
        let dot3 = |a: [f64; 3], b: [f64; 3]| a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
        let mul_diag = |d: [f64; 3], v: [f64; 3]| [d[0] * v[0], d[1] * v[1], d[2] * v[2]];

        let inv_mass_a = 1.0_f64;
        let inv_inertia_a = [6.0_f64, 6.0, 6.0];
        let c_a = [ta.position.x, ta.position.y, ta.position.z];
        let mut v_a = [0.0_f64, -1.0, 0.0];
        let mut w_a = [0.0_f64, 0.0, 0.0];

        // Accumulated normal impulse per contact. Sequential impulse clamps the
        // TOTAL accumulated impulse to be non-negative (Catto/Box2D form), not
        // each per-iteration increment. Per-increment clamping admits a spurious
        // fixed point where the body translates up and spins while every contact
        // shows non-negative normal velocity; accumulation lets a later pass take
        // back earlier over-impulse and drives a symmetric resting box to w = 0.
        let n_contacts = manifold.contacts.len();
        let mut acc = vec![0.0_f64; n_contacts];

        let run_pass = |v_a: &mut [f64; 3], w_a: &mut [f64; 3], acc: &mut [f64]| {
            for (i, c) in manifold.contacts.iter().enumerate() {
                let r_a = [
                    c.point_a.x - c_a[0],
                    c.point_a.y - c_a[1],
                    c.point_a.z - c_a[2],
                ];
                let n = [c.normal.x, c.normal.y, c.normal.z];
                let wxr = cross3(*w_a, r_a);
                let v_point_a = [v_a[0] + wxr[0], v_a[1] + wxr[1], v_a[2] + wxr[2]];
                let vn = dot3(v_point_a, n);
                let rxn = cross3(r_a, n);
                let denom = inv_mass_a + dot3(n, cross3(mul_diag(inv_inertia_a, rxn), r_a));
                let d_lambda = -vn / denom;
                let old = acc[i];
                let new = (old + d_lambda).max(0.0);
                let applied = new - old;
                acc[i] = new;
                v_a[0] += n[0] * inv_mass_a * applied;
                v_a[1] += n[1] * inv_mass_a * applied;
                v_a[2] += n[2] * inv_mass_a * applied;
                let dl = cross3(r_a, [n[0] * applied, n[1] * applied, n[2] * applied]);
                let dw = mul_diag(inv_inertia_a, dl);
                w_a[0] += dw[0];
                w_a[1] += dw[1];
                w_a[2] += dw[2];
            }
        };

        for _ in 0..10 {
            run_pass(&mut v_a, &mut w_a, &mut acc);
        }
        let w_norm = (w_a[0] * w_a[0] + w_a[1] * w_a[1] + w_a[2] * w_a[2]).sqrt();
        assert!(w_norm < 1e-6, "resting box should not spin: |w| = {w_norm}");

        for _ in 0..30 {
            // Re-drive downward each step; warm-start the accumulator (as a real
            // solver does) so the contact keeps holding the box flat.
            v_a = [0.0, -1.0, 0.0];
            for _ in 0..10 {
                run_pass(&mut v_a, &mut w_a, &mut acc);
            }
            let wn = (w_a[0] * w_a[0] + w_a[1] * w_a[1] + w_a[2] * w_a[2]).sqrt();
            assert!(wn < 1e-3, "box rocked during settling: |w| = {wn}");
        }
    }

    #[test]
    fn test_two_stacked_boxes_stable() {
        let ground = BoxShape::new(Vec3::new(10.0, 0.5, 10.0));
        let cube_a = BoxShape::new(Vec3::new(0.5, 0.5, 0.5));
        let cube_b = BoxShape::new(Vec3::new(0.5, 0.5, 0.5));

        let t_ground = Transform::from_position(Vec3::new(0.0, -0.5, 0.0));
        let t_a = Transform::from_position(Vec3::new(0.0, 0.5 - 0.005, 0.0));
        let t_b = Transform::from_position(Vec3::new(0.0, 1.5 - 0.01, 0.0));

        // ground (A=lower) vs cube A (B=upper): normal B->A points down.
        let m_ground_a =
            box_box_manifold(&ground, &t_ground, &cube_a, &t_a, CollisionPair::new(0, 1))
                .manifold
                .expect("ground-A overlap");
        assert_eq!(m_ground_a.contacts.len(), 4);
        // cube A (lower) vs cube B (upper).
        let m_a_b = box_box_manifold(&cube_a, &t_a, &cube_b, &t_b, CollisionPair::new(1, 2))
            .manifold
            .expect("A-B overlap");
        assert_eq!(m_a_b.contacts.len(), 4);

        let cross3 = |a: [f64; 3], b: [f64; 3]| {
            [
                a[1] * b[2] - a[2] * b[1],
                a[2] * b[0] - a[0] * b[2],
                a[0] * b[1] - a[1] * b[0],
            ]
        };
        let dot3 = |a: [f64; 3], b: [f64; 3]| a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
        let mul_diag = |d: [f64; 3], v: [f64; 3]| [d[0] * v[0], d[1] * v[1], d[2] * v[2]];

        // Bodies: index 0 = ground (static), 1 = cube A, 2 = cube B.
        let inv_mass = [0.0_f64, 1.0, 1.0];
        let inv_inertia = [[0.0_f64, 0.0, 0.0], [6.0, 6.0, 6.0], [6.0, 6.0, 6.0]];
        let centers = [
            [
                t_ground.position.x,
                t_ground.position.y,
                t_ground.position.z,
            ],
            [t_a.position.x, t_a.position.y, t_a.position.z],
            [t_b.position.x, t_b.position.y, t_b.position.z],
        ];
        let mut vel = [[0.0_f64; 3]; 3];
        let mut omega = [[0.0_f64; 3]; 3];

        let b_start = centers[2];

        // Accumulated normal impulse per contact, per manifold (warm-started
        // across steps). Sequential impulse clamps the TOTAL accumulated impulse
        // non-negative, which is what keeps a stack quiet: per-increment clamping
        // admits spurious spin/drift fixed points.
        let mut acc_ground = vec![0.0_f64; m_ground_a.contacts.len()];
        let mut acc_ab = vec![0.0_f64; m_a_b.contacts.len()];

        // Solve one contact between bodies `ia` (A side) and `ib` (B side),
        // accumulating its normal impulse in `acc`.
        let solve_contact = |vel: &mut [[f64; 3]; 3],
                             omega: &mut [[f64; 3]; 3],
                             acc: &mut f64,
                             ia: usize,
                             ib: usize,
                             point_a: [f64; 3],
                             point_b: [f64; 3],
                             n: [f64; 3]| {
            let r_a = [
                point_a[0] - centers[ia][0],
                point_a[1] - centers[ia][1],
                point_a[2] - centers[ia][2],
            ];
            let r_b = [
                point_b[0] - centers[ib][0],
                point_b[1] - centers[ib][1],
                point_b[2] - centers[ib][2],
            ];
            let wxr_a = cross3(omega[ia], r_a);
            let wxr_b = cross3(omega[ib], r_b);
            let vp_a = [
                vel[ia][0] + wxr_a[0],
                vel[ia][1] + wxr_a[1],
                vel[ia][2] + wxr_a[2],
            ];
            let vp_b = [
                vel[ib][0] + wxr_b[0],
                vel[ib][1] + wxr_b[1],
                vel[ib][2] + wxr_b[2],
            ];
            let v_rel = [vp_a[0] - vp_b[0], vp_a[1] - vp_b[1], vp_a[2] - vp_b[2]];
            let vn = dot3(v_rel, n);
            let rxn_a = cross3(r_a, n);
            let rxn_b = cross3(r_b, n);
            let ang_a = dot3(n, cross3(mul_diag(inv_inertia[ia], rxn_a), r_a));
            let ang_b = dot3(n, cross3(mul_diag(inv_inertia[ib], rxn_b), r_b));
            let denom = inv_mass[ia] + inv_mass[ib] + ang_a + ang_b;
            if denom <= 0.0 {
                return;
            }
            let d_lambda = -vn / denom;
            let old = *acc;
            let new = (old + d_lambda).max(0.0);
            let applied = new - old;
            *acc = new;
            let imp = [n[0] * applied, n[1] * applied, n[2] * applied];
            vel[ia][0] += imp[0] * inv_mass[ia];
            vel[ia][1] += imp[1] * inv_mass[ia];
            vel[ia][2] += imp[2] * inv_mass[ia];
            vel[ib][0] -= imp[0] * inv_mass[ib];
            vel[ib][1] -= imp[1] * inv_mass[ib];
            vel[ib][2] -= imp[2] * inv_mass[ib];
            let dl_a = mul_diag(inv_inertia[ia], cross3(r_a, imp));
            let dl_b = mul_diag(inv_inertia[ib], cross3(r_b, imp));
            omega[ia][0] += dl_a[0];
            omega[ia][1] += dl_a[1];
            omega[ia][2] += dl_a[2];
            omega[ib][0] -= dl_b[0];
            omega[ib][1] -= dl_b[1];
            omega[ib][2] -= dl_b[2];
        };

        for _ in 0..60 {
            // Gravity drive: both cubes pushed down each step.
            vel[1][1] = -1.0;
            vel[2][1] = -1.0;
            // Several Gauss-Seidel sweeps over both manifolds per step, with the
            // accumulators warm-started across steps (stack stays at rest).
            for _ in 0..12 {
                // ground(0)-A(1): A side = ground (0), B side = cube A (1).
                for (i, c) in m_ground_a.contacts.iter().enumerate() {
                    solve_contact(
                        &mut vel,
                        &mut omega,
                        &mut acc_ground[i],
                        0,
                        1,
                        [c.point_a.x, c.point_a.y, c.point_a.z],
                        [c.point_b.x, c.point_b.y, c.point_b.z],
                        [c.normal.x, c.normal.y, c.normal.z],
                    );
                }
                // A(1)-B(2): A side = cube A (1), B side = cube B (2).
                for (i, c) in m_a_b.contacts.iter().enumerate() {
                    solve_contact(
                        &mut vel,
                        &mut omega,
                        &mut acc_ab[i],
                        1,
                        2,
                        [c.point_a.x, c.point_a.y, c.point_a.z],
                        [c.point_b.x, c.point_b.y, c.point_b.z],
                        [c.normal.x, c.normal.y, c.normal.z],
                    );
                }
            }
        }

        // Top box must not drift horizontally nor spin.
        let dx = (vel[2][0]).abs();
        let dz = (vel[2][2]).abs();
        assert!(dx < 1e-3, "top box drifted in x: {dx}");
        assert!(dz < 1e-3, "top box drifted in z: {dz}");
        let _ = b_start;
        let wb =
            (omega[2][0] * omega[2][0] + omega[2][1] * omega[2][1] + omega[2][2] * omega[2][2])
                .sqrt();
        assert!(wb < 1e-3, "top box spun: |w_b| = {wb}");
    }

    #[test]
    fn test_feature_id_stability_under_nudge() {
        let box_a = BoxShape::new(Vec3::new(0.5, 0.5, 0.5));
        let ground = BoxShape::new(Vec3::new(10.0, 0.5, 10.0));
        let t_ground = Transform::from_position(Vec3::new(0.0, -0.5, 0.0));

        let ta1 = Transform::from_position(Vec3::new(0.0, 0.5 - 0.005, 0.0));
        let ta2 = Transform::from_position(Vec3::new(0.001, 0.5 - 0.005, 0.0));

        let m1 = box_box_manifold(&box_a, &ta1, &ground, &t_ground, CollisionPair::new(0, 1))
            .manifold
            .expect("m1");
        let m2 = box_box_manifold(&box_a, &ta2, &ground, &t_ground, CollisionPair::new(0, 1))
            .manifold
            .expect("m2");

        assert_eq!(m1.contacts.len(), 4);
        assert_eq!(m2.contacts.len(), 4);

        for i in 0..m1.contacts.len() {
            let p1 = m1.contacts[i].point_b;
            let p2 = m2.contacts[i].point_b;
            assert!((p2.x - p1.x - 0.001).abs() < 5e-4);
            assert!((p2.y - p1.y).abs() < 5e-4);
            assert!((p2.z - p1.z).abs() < 5e-4);
        }
    }

    #[test]
    fn test_normal_sign_both_orders() {
        let box_a = BoxShape::new(Vec3::new(0.5, 0.5, 0.5));
        let box_b = BoxShape::new(Vec3::new(0.5, 0.5, 0.5));
        let tb = Transform::from_position(Vec3::new(0.0, 0.0, 0.0));
        let ta = Transform::from_position(Vec3::new(0.0, 0.99, 0.0));

        let m1 = box_box_manifold(&box_a, &ta, &box_b, &tb, CollisionPair::new(0, 1))
            .manifold
            .expect("m1");
        let m2 = box_box_manifold(&box_b, &tb, &box_a, &ta, CollisionPair::new(1, 0))
            .manifold
            .expect("m2");

        assert!(!m1.contacts.is_empty());
        assert!(!m2.contacts.is_empty());
        assert_eq!(m1.contacts.len(), m2.contacts.len());

        for c1 in &m1.contacts {
            let p1 = c1.point();
            let mut best = 0usize;
            let mut best_d = f64::INFINITY;
            for (k, c2) in m2.contacts.iter().enumerate() {
                let d = (p1 - c2.point()).norm();
                if d < best_d {
                    best_d = d;
                    best = k;
                }
            }
            let c2 = &m2.contacts[best];
            assert!((c2.normal - (-c1.normal)).norm() < 1e-6);
            assert!((c1.depth - c2.depth).abs() < 1e-6);
            assert!((c1.point_a - c2.point_b).norm() < 1e-6);
            assert!((c1.point_b - c2.point_a).norm() < 1e-6);
        }
    }
}
