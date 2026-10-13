// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Generic convex-convex narrow phase: GJK distance and EPA penetration depth.
//!
//! Every convex collider is a convex *core* inflated by a rounding margin
//! (see [`Core`]).  GJK is run on the cores to obtain either their separation
//! distance and closest points, or - when the cores overlap - a simplex that
//! seeds the expanding polytope algorithm (EPA), which then yields the exact
//! minimum translation vector.  The margins are applied analytically
//! afterwards, so round shapes (spheres, capsules) never degrade the
//! conditioning of the iteration.
//!
//! The Minkowski difference of the two cores is `M = A - B`; a point of `M`
//! is stored together with the two shape points that produced it so that
//! witness points on both shapes can be recovered from barycentric weights.

use super::geom::{Vec3, cross, dot, length, length_sq, madd, neg, normalized, scale, sub};
use super::narrowphase::Contact;
use super::shapes::Core;

/// Maximum number of GJK iterations.
const GJK_MAX_ITERATIONS: usize = 64;

/// Maximum number of EPA polytope expansions.
const EPA_MAX_ITERATIONS: usize = 128;

/// Squared distance below which the origin is considered to lie on the simplex.
const ORIGIN_EPS_SQ: f64 = 1e-24;

/// Relative progress threshold that terminates GJK.
const GJK_RELATIVE_TOLERANCE: f64 = 1e-12;

/// Relative tolerance of the EPA termination test.
const EPA_RELATIVE_TOLERANCE: f64 = 1e-9;

/// A vertex of the Minkowski difference with its generating shape points.
#[derive(Debug, Clone, Copy)]
struct SupportPoint {
    /// `a - b`.
    w: Vec3,
    /// Support point on shape A.
    a: Vec3,
    /// Support point on shape B.
    b: Vec3,
}

/// Support point of `A - B` in direction `dir`.
fn minkowski_support(core_a: &Core, core_b: &Core, dir: Vec3) -> SupportPoint {
    let a = core_a.support(dir);
    let b = core_b.support(neg(dir));
    SupportPoint { w: sub(a, b), a, b }
}

/// Closest point to the origin on triangle `abc`, with barycentric weights
/// (Ericson, *Real-Time Collision Detection*, 5.1.5 specialised to `p = 0`).
fn closest_on_triangle(a: Vec3, b: Vec3, c: Vec3) -> (Vec3, [f64; 3]) {
    let ab = sub(b, a);
    let ac = sub(c, a);
    let ap = neg(a);
    let d1 = dot(ab, ap);
    let d2 = dot(ac, ap);
    if d1 <= 0.0 && d2 <= 0.0 {
        return (a, [1.0, 0.0, 0.0]);
    }
    let bp = neg(b);
    let d3 = dot(ab, bp);
    let d4 = dot(ac, bp);
    if d3 >= 0.0 && d4 <= d3 {
        return (b, [0.0, 1.0, 0.0]);
    }
    let vc = d1 * d4 - d3 * d2;
    if vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0 {
        let denom = d1 - d3;
        let v = if denom.abs() > 0.0 { d1 / denom } else { 0.0 };
        return (madd(a, ab, v), [1.0 - v, v, 0.0]);
    }
    let cp = neg(c);
    let d5 = dot(ab, cp);
    let d6 = dot(ac, cp);
    if d6 >= 0.0 && d5 <= d6 {
        return (c, [0.0, 0.0, 1.0]);
    }
    let vb = d5 * d2 - d1 * d6;
    if vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0 {
        let denom = d2 - d6;
        let w = if denom.abs() > 0.0 { d2 / denom } else { 0.0 };
        return (madd(a, ac, w), [1.0 - w, 0.0, w]);
    }
    let va = d3 * d6 - d5 * d4;
    if va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0 {
        let denom = (d4 - d3) + (d5 - d6);
        let w = if denom.abs() > 0.0 {
            (d4 - d3) / denom
        } else {
            0.0
        };
        return (madd(b, sub(c, b), w), [0.0, 1.0 - w, w]);
    }
    let sum = va + vb + vc;
    if sum.abs() <= 0.0 {
        return (a, [1.0, 0.0, 0.0]);
    }
    let denom = 1.0 / sum;
    let v = vb * denom;
    let w = vc * denom;
    (madd(madd(a, ab, v), ac, w), [1.0 - v - w, v, w])
}

/// Result of reducing a simplex towards the origin.
struct Reduction {
    /// Closest point of the simplex to the origin.
    point: Vec3,
    /// Retained vertices as `(index into the old simplex, barycentric weight)`.
    kept: Vec<(usize, f64)>,
    /// `true` when the origin lies inside a (non-degenerate) tetrahedron.
    inside: bool,
}

/// Find the point of the simplex closest to the origin and the minimal
/// sub-simplex that carries it.
fn reduce_simplex(simplex: &[SupportPoint]) -> Reduction {
    match simplex.len() {
        0 => Reduction {
            point: [0.0; 3],
            kept: Vec::new(),
            inside: false,
        },
        1 => Reduction {
            point: simplex[0].w,
            kept: vec![(0, 1.0)],
            inside: false,
        },
        2 => {
            let a = simplex[0].w;
            let b = simplex[1].w;
            let ab = sub(b, a);
            let len_sq = dot(ab, ab);
            if len_sq <= 1e-300 {
                return Reduction {
                    point: a,
                    kept: vec![(0, 1.0)],
                    inside: false,
                };
            }
            let t = -dot(a, ab) / len_sq;
            if t <= 0.0 {
                Reduction {
                    point: a,
                    kept: vec![(0, 1.0)],
                    inside: false,
                }
            } else if t >= 1.0 {
                Reduction {
                    point: b,
                    kept: vec![(1, 1.0)],
                    inside: false,
                }
            } else {
                Reduction {
                    point: madd(a, ab, t),
                    kept: vec![(0, 1.0 - t), (1, t)],
                    inside: false,
                }
            }
        }
        3 => {
            let (point, w) = closest_on_triangle(simplex[0].w, simplex[1].w, simplex[2].w);
            let kept = w
                .iter()
                .enumerate()
                .filter(|(_, weight)| **weight > 0.0)
                .map(|(i, weight)| (i, *weight))
                .collect();
            Reduction {
                point,
                kept,
                inside: false,
            }
        }
        _ => reduce_tetrahedron(simplex),
    }
}

/// [`reduce_simplex`] for four vertices.
fn reduce_tetrahedron(simplex: &[SupportPoint]) -> Reduction {
    const FACES: [([usize; 3], usize); 4] = [
        ([0, 1, 2], 3),
        ([0, 3, 1], 2),
        ([0, 2, 3], 1),
        ([1, 3, 2], 0),
    ];
    let v: Vec<Vec3> = simplex.iter().map(|s| s.w).collect();
    let extent = v
        .iter()
        .fold(0.0_f64, |m, p| {
            m.max(p[0].abs()).max(p[1].abs()).max(p[2].abs())
        })
        .max(1e-300);
    let volume = dot(sub(v[1], v[0]), cross(sub(v[2], v[0]), sub(v[3], v[0])));
    let degenerate = volume.abs() <= 1e-12 * extent * extent * extent;

    let mut best: Option<(f64, Vec3, [usize; 3], [f64; 3])> = None;
    let mut outside_any = false;
    for (face, opposite) in FACES {
        let [i, j, k] = face;
        let normal = cross(sub(v[j], v[i]), sub(v[k], v[i]));
        let outside = if degenerate {
            true
        } else if let Some(mut n) = normalized(normal) {
            if dot(n, sub(v[opposite], v[i])) > 0.0 {
                n = neg(n);
            }
            dot(n, neg(v[i])) > 1e-12 * extent
        } else {
            true
        };
        if !outside {
            continue;
        }
        outside_any = true;
        let (point, weights) = closest_on_triangle(v[i], v[j], v[k]);
        let d2 = dot(point, point);
        if best.as_ref().is_none_or(|(bd, _, _, _)| d2 < *bd) {
            best = Some((d2, point, face, weights));
        }
    }
    match best {
        Some((_, point, face, weights)) if outside_any => {
            let kept = face
                .iter()
                .zip(weights.iter())
                .filter(|(_, w)| **w > 0.0)
                .map(|(i, w)| (*i, *w))
                .collect();
            Reduction {
                point,
                kept,
                inside: false,
            }
        }
        _ => Reduction {
            point: [0.0; 3],
            kept: vec![(0, 0.25), (1, 0.25), (2, 0.25), (3, 0.25)],
            inside: true,
        },
    }
}

/// Outcome of the GJK distance query.
enum GjkOutcome {
    /// The cores are disjoint.
    Separated {
        /// Distance between the cores.
        distance: f64,
        /// Closest point on core A.
        closest_a: Vec3,
        /// Closest point on core B.
        closest_b: Vec3,
    },
    /// The cores touch or overlap; the simplex encloses (or contains) the origin.
    Overlapping(Vec<SupportPoint>),
}

/// Weighted combination of the simplex points on shape A and shape B.
fn witness_points(simplex: &[SupportPoint], weights: &[f64]) -> (Vec3, Vec3) {
    let mut pa = [0.0; 3];
    let mut pb = [0.0; 3];
    for (s, w) in simplex.iter().zip(weights.iter()) {
        pa = madd(pa, s.a, *w);
        pb = madd(pb, s.b, *w);
    }
    (pa, pb)
}

/// Run GJK on the two cores.
fn gjk(core_a: &Core, core_b: &Core) -> GjkOutcome {
    let mut dir = sub(core_a.center(), core_b.center());
    if length_sq(dir) < ORIGIN_EPS_SQ {
        dir = [1.0, 0.0, 0.0];
    }
    let mut simplex = vec![minkowski_support(core_a, core_b, dir)];
    let mut weights = vec![1.0];
    let mut v = simplex[0].w;
    for _ in 0..GJK_MAX_ITERATIONS {
        let v_sq = length_sq(v);
        if v_sq < ORIGIN_EPS_SQ {
            return GjkOutcome::Overlapping(simplex);
        }
        let w = minkowski_support(core_a, core_b, neg(v));
        // `v` is the point of the current simplex closest to the origin and
        // `w` the extreme point of M along -v; if `w` does not lie measurably
        // beyond `v` the iteration has converged on the closest point of M.
        let progress = v_sq - dot(v, w.w);
        if progress <= GJK_RELATIVE_TOLERANCE * v_sq {
            break;
        }
        simplex.push(w);
        let reduction = reduce_simplex(&simplex);
        if reduction.inside {
            return GjkOutcome::Overlapping(simplex);
        }
        let reduced: Vec<SupportPoint> = reduction.kept.iter().map(|(i, _)| simplex[*i]).collect();
        weights = reduction.kept.iter().map(|(_, w)| *w).collect();
        simplex = reduced;
        v = reduction.point;
    }
    let (closest_a, closest_b) = witness_points(&simplex, &weights);
    GjkOutcome::Separated {
        distance: length(v),
        closest_a,
        closest_b,
    }
}

/// Grow a degenerate simplex (fewer than four affinely independent vertices)
/// into a tetrahedron of Minkowski-difference points.
///
/// Returns `None` when the Minkowski difference is flat and no such
/// tetrahedron exists.
fn complete_tetrahedron(
    core_a: &Core,
    core_b: &Core,
    simplex: &mut Vec<SupportPoint>,
) -> Option<()> {
    const AXES: [Vec3; 6] = [
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
        [-1.0, 0.0, 0.0],
        [0.0, -1.0, 0.0],
        [0.0, 0.0, -1.0],
    ];
    let extent =
        |pts: &[SupportPoint]| -> f64 { pts.iter().fold(1e-12_f64, |m, p| m.max(length(p.w))) };
    // Drop duplicate vertices first.
    let mut unique: Vec<SupportPoint> = Vec::with_capacity(4);
    for s in simplex.iter() {
        let scale_ref = extent(simplex);
        if unique
            .iter()
            .all(|u| length(sub(u.w, s.w)) > 1e-10 * scale_ref)
        {
            unique.push(*s);
        }
    }
    *simplex = unique;
    if simplex.is_empty() {
        simplex.push(minkowski_support(core_a, core_b, AXES[0]));
    }
    while simplex.len() < 4 {
        let scale_ref = extent(simplex).max(1.0);
        let candidates: Vec<Vec3> = match simplex.len() {
            1 => AXES.to_vec(),
            2 => {
                let edge = sub(simplex[1].w, simplex[0].w);
                let mut dirs = Vec::with_capacity(12);
                for axis in AXES {
                    let d = cross(edge, axis);
                    if length_sq(d) > 1e-20 {
                        dirs.push(d);
                    }
                }
                dirs
            }
            _ => {
                let n = cross(
                    sub(simplex[1].w, simplex[0].w),
                    sub(simplex[2].w, simplex[0].w),
                );
                vec![n, neg(n)]
            }
        };
        let mut added = false;
        for dir in candidates {
            let s = minkowski_support(core_a, core_b, dir);
            let independent = match simplex.len() {
                1 => length(sub(s.w, simplex[0].w)) > 1e-9 * scale_ref,
                2 => {
                    let e = sub(simplex[1].w, simplex[0].w);
                    length(cross(e, sub(s.w, simplex[0].w))) > 1e-9 * scale_ref * scale_ref
                }
                _ => {
                    let e1 = sub(simplex[1].w, simplex[0].w);
                    let e2 = sub(simplex[2].w, simplex[0].w);
                    dot(cross(e1, e2), sub(s.w, simplex[0].w)).abs()
                        > 1e-9 * scale_ref * scale_ref * scale_ref
                }
            };
            if independent {
                simplex.push(s);
                added = true;
                break;
            }
        }
        if !added {
            return None;
        }
    }
    Some(())
}

/// A triangular face of the EPA polytope.
#[derive(Debug, Clone, Copy)]
struct Face {
    /// Vertex indices, wound so that `normal` is the outward normal.
    idx: [usize; 3],
    /// Outward unit normal.
    normal: Vec3,
    /// Signed distance of the face plane from the origin along `normal`.
    dist: f64,
}

/// Build a face from three polytope vertices, orienting its normal away from
/// `interior`.
fn make_face(verts: &[SupportPoint], idx: [usize; 3], interior: Vec3) -> Option<Face> {
    let [i, j, k] = idx;
    let a = verts[i].w;
    let edge_1 = sub(verts[j].w, a);
    let edge_2 = sub(verts[k].w, a);
    let area_vector = cross(edge_1, edge_2);
    // A sliver whose area is lost in rounding has no reliable normal.
    if length_sq(area_vector) <= 1e-24 * length_sq(edge_1) * length_sq(edge_2) {
        return None;
    }
    let normal = normalized(area_vector)?;
    if dot(normal, sub(a, interior)) >= 0.0 {
        Some(Face {
            idx,
            normal,
            dist: dot(normal, a),
        })
    } else {
        let flipped = neg(normal);
        Some(Face {
            idx: [i, k, j],
            normal: flipped,
            dist: dot(flipped, a),
        })
    }
}

/// Penetration of two overlapping cores as found by EPA.
struct EpaResult {
    /// Outward normal of `M` at the boundary point closest to the origin.
    normal: Vec3,
    /// Distance from the origin to that boundary point.
    depth: f64,
    /// Witness point on core A.
    point_a: Vec3,
    /// Witness point on core B.
    point_b: Vec3,
}

/// Expanding polytope algorithm seeded with a tetrahedron around the origin.
fn epa(core_a: &Core, core_b: &Core, start: &[SupportPoint]) -> Option<EpaResult> {
    let mut verts: Vec<SupportPoint> = start.to_vec();
    let interior = {
        let mut c = [0.0; 3];
        for s in &verts {
            c = madd(c, s.w, 1.0 / verts.len() as f64);
        }
        c
    };
    let mut faces: Vec<Face> = Vec::with_capacity(32);
    for idx in [[0, 1, 2], [0, 3, 1], [0, 2, 3], [1, 3, 2]] {
        faces.push(make_face(&verts, idx, interior)?);
    }
    let mut best = faces[0];
    for _ in 0..EPA_MAX_ITERATIONS {
        best = *faces.iter().min_by(|a, b| a.dist.total_cmp(&b.dist))?;
        let w = minkowski_support(core_a, core_b, best.normal);
        let reach = dot(best.normal, w.w);
        if reach - best.dist <= EPA_RELATIVE_TOLERANCE * (1.0 + reach.abs()) {
            break;
        }
        // Remove the faces visible from the new point and stitch the horizon.
        let visible: Vec<bool> = faces
            .iter()
            .map(|f| dot(f.normal, sub(w.w, verts[f.idx[0]].w)) > 1e-14)
            .collect();
        if !visible.iter().any(|v| *v) {
            break;
        }
        let mut horizon: Vec<(usize, usize)> = Vec::new();
        for (face, is_visible) in faces.iter().zip(visible.iter()) {
            if !*is_visible {
                continue;
            }
            for e in 0..3 {
                let from = face.idx[e];
                let to = face.idx[(e + 1) % 3];
                let shared = faces
                    .iter()
                    .zip(visible.iter())
                    .any(|(other, other_visible)| {
                        *other_visible
                            && (0..3).any(|f| other.idx[f] == to && other.idx[(f + 1) % 3] == from)
                    });
                if !shared {
                    horizon.push((from, to));
                }
            }
        }
        let new_index = verts.len();
        verts.push(w);
        let mut kept: Vec<Face> = faces
            .iter()
            .zip(visible.iter())
            .filter(|(_, v)| !**v)
            .map(|(f, _)| *f)
            .collect();
        for (from, to) in horizon {
            if let Some(face) = make_face(&verts, [from, to, new_index], interior) {
                kept.push(face);
            }
        }
        if kept.len() < 4 {
            break;
        }
        faces = kept;
    }
    let p = scale(best.normal, best.dist);
    let a = verts[best.idx[0]];
    let b = verts[best.idx[1]];
    let c = verts[best.idx[2]];
    let weights = barycentric(p, a.w, b.w, c.w);
    let (point_a, point_b) = witness_points(&[a, b, c], &weights);
    Some(EpaResult {
        normal: best.normal,
        depth: best.dist.max(0.0),
        point_a,
        point_b,
    })
}

/// Barycentric coordinates of `p` (assumed to lie on the plane of `abc`),
/// clamped to the triangle.
fn barycentric(p: Vec3, a: Vec3, b: Vec3, c: Vec3) -> [f64; 3] {
    let v0 = sub(b, a);
    let v1 = sub(c, a);
    let v2 = sub(p, a);
    let d00 = dot(v0, v0);
    let d01 = dot(v0, v1);
    let d11 = dot(v1, v1);
    let d20 = dot(v2, v0);
    let d21 = dot(v2, v1);
    let denom = d00 * d11 - d01 * d01;
    if denom.abs() <= 1e-300 {
        return [1.0 / 3.0; 3];
    }
    let v = (d11 * d20 - d01 * d21) / denom;
    let w = (d00 * d21 - d01 * d20) / denom;
    let u = 1.0 - v - w;
    let clamped = [u.max(0.0), v.max(0.0), w.max(0.0)];
    let sum = clamped[0] + clamped[1] + clamped[2];
    if sum <= 0.0 {
        [1.0 / 3.0; 3]
    } else {
        [clamped[0] / sum, clamped[1] / sum, clamped[2] / sum]
    }
}

/// Contact between two convex shapes given as `core + margin`.
///
/// The returned normal points from B toward A; `point_a` / `point_b` lie on the
/// (rounded) surfaces of the two shapes.  `None` when the shapes are disjoint.
pub(super) fn convex_contact(
    core_a: &Core,
    margin_a: f64,
    core_b: &Core,
    margin_b: f64,
) -> Option<Contact> {
    let reach = margin_a + margin_b;
    match gjk(core_a, core_b) {
        GjkOutcome::Separated {
            distance,
            closest_a,
            closest_b,
        } => {
            if distance >= reach {
                return None;
            }
            // `closest_a - closest_b` points from B toward A.
            let normal = normalized(sub(closest_a, closest_b))?;
            Some(Contact {
                normal,
                depth: reach - distance,
                point_a: madd(closest_a, normal, -margin_a),
                point_b: madd(closest_b, normal, margin_b),
            })
        }
        GjkOutcome::Overlapping(mut simplex) => {
            complete_tetrahedron(core_a, core_b, &mut simplex)?;
            let result = epa(core_a, core_b, &simplex)?;
            // EPA's normal is the outward normal of M = A - B, i.e. the
            // direction in which B must move relative to A; the contact
            // normal (B toward A) is its opposite.
            let normal = neg(result.normal);
            Some(Contact {
                normal,
                depth: result.depth + reach,
                point_a: madd(result.point_a, normal, -margin_a),
                point_b: madd(result.point_b, normal, margin_b),
            })
        }
    }
}
