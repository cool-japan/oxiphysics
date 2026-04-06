//! Self-contained narrow-phase collision detection using plain f64 arrays.
//!
//! This module provides standalone functions that operate on raw geometric
//! data without any external dependencies on nalgebra or shape abstractions.

// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

#[allow(unused_imports)]
use super::functions::*;

// ── Vector helpers ────────────────────────────────────────────────────────────

#[allow(dead_code)]
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[allow(dead_code)]
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

#[allow(dead_code)]
fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

#[allow(dead_code)]
fn add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

#[allow(dead_code)]
fn scale(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

#[allow(dead_code)]
fn len(a: [f64; 3]) -> f64 {
    (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt()
}

#[allow(dead_code)]
fn len_sq(a: [f64; 3]) -> f64 {
    a[0] * a[0] + a[1] * a[1] + a[2] * a[2]
}

#[allow(dead_code)]
fn normalize(a: [f64; 3]) -> [f64; 3] {
    let l = len(a);
    if l > 1e-10 {
        scale(a, 1.0 / l)
    } else {
        [0.0, 1.0, 0.0]
    }
}

#[allow(dead_code)]
fn neg(a: [f64; 3]) -> [f64; 3] {
    [-a[0], -a[1], -a[2]]
}

// ── Contact types ─────────────────────────────────────────────────────────────

/// A single contact point between two shapes.
#[allow(dead_code)]
#[derive(Debug, Clone, PartialEq)]
pub struct Contact {
    /// Witness point on shape A in world space.
    pub point_a: [f64; 3],
    /// Witness point on shape B in world space.
    pub point_b: [f64; 3],
    /// Contact normal pointing from B toward A (unit length).
    pub normal: [f64; 3],
    /// Penetration depth (positive means overlapping).
    pub depth: f64,
}

/// Result of a narrow-phase collision query.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub enum NarrowPhaseResult {
    /// The shapes are separated; no contact.
    NoContact,
    /// A single contact point was found.
    SingleContact(Contact),
    /// Multiple contact points were found (contact manifold).
    Manifold(Vec<Contact>),
}

impl NarrowPhaseResult {
    /// Returns `true` if there is at least one contact.
    #[allow(dead_code)]
    pub fn has_contact(&self) -> bool {
        !matches!(self, NarrowPhaseResult::NoContact)
    }

    /// Returns the first contact if any.
    #[allow(dead_code)]
    pub fn first_contact(&self) -> Option<&Contact> {
        match self {
            NarrowPhaseResult::NoContact => None,
            NarrowPhaseResult::SingleContact(c) => Some(c),
            NarrowPhaseResult::Manifold(v) => v.first(),
        }
    }
}

// ── Sphere vs Sphere ──────────────────────────────────────────────────────────

/// Sphere–sphere narrow-phase test.
///
/// Returns `NoContact` when the spheres are separated, `SingleContact`
/// when they overlap.
#[allow(dead_code)]
pub fn sphere_sphere(ca: [f64; 3], ra: f64, cb: [f64; 3], rb: f64) -> NarrowPhaseResult {
    let d = sub(cb, ca);
    let dist = len(d);
    let depth = ra + rb - dist;
    if depth < 0.0 {
        return NarrowPhaseResult::NoContact;
    }
    let normal = if dist > 1e-10 {
        scale(d, 1.0 / dist)
    } else {
        [0.0, 1.0, 0.0]
    };
    let point_a = add(ca, scale(normal, ra));
    let point_b = sub(cb, scale(normal, rb));
    NarrowPhaseResult::SingleContact(Contact {
        point_a,
        point_b,
        normal,
        depth,
    })
}

// ── Sphere vs Box ─────────────────────────────────────────────────────────────

/// Sphere–axis-aligned box narrow-phase test.
///
/// The box is defined by its `box_center` and `half_extents`.
/// Returns `NoContact` when the sphere does not reach the box.
#[allow(dead_code)]
pub fn sphere_box(
    sphere_center: [f64; 3],
    radius: f64,
    box_center: [f64; 3],
    half_extents: [f64; 3],
) -> NarrowPhaseResult {
    // Closest point on box to sphere center (work in box-local space).
    let local = sub(sphere_center, box_center);
    let clamped = [
        local[0].clamp(-half_extents[0], half_extents[0]),
        local[1].clamp(-half_extents[1], half_extents[1]),
        local[2].clamp(-half_extents[2], half_extents[2]),
    ];
    let diff = sub(local, clamped);
    let dist_sq = len_sq(diff);

    if dist_sq >= radius * radius {
        return NarrowPhaseResult::NoContact;
    }

    let dist = dist_sq.sqrt();
    let (normal, depth) = if dist > 1e-10 {
        // Sphere centre is outside the box.
        (normalize(diff), radius - dist)
    } else {
        // Sphere centre is inside the box; push out through nearest face.
        let dx = half_extents[0] - local[0].abs();
        let dy = half_extents[1] - local[1].abs();
        let dz = half_extents[2] - local[2].abs();
        let n = if dx <= dy && dx <= dz {
            [local[0].signum(), 0.0, 0.0]
        } else if dy <= dz {
            [0.0, local[1].signum(), 0.0]
        } else {
            [0.0, 0.0, local[2].signum()]
        };
        let pen = dx.min(dy).min(dz) + radius;
        (n, pen)
    };

    // Convert closest point back to world space.
    let closest_world = add(box_center, clamped);
    let point_a = sub(sphere_center, scale(normal, radius));
    let point_b = closest_world;

    NarrowPhaseResult::SingleContact(Contact {
        point_a,
        point_b,
        normal,
        depth,
    })
}

// ── Box vs Box (SAT) ──────────────────────────────────────────────────────────

/// OBB–OBB SAT collision test using 15 separating axes.
///
/// `center_a`/`center_b` are world-space box centres.
/// `half_a`/`half_b` are half-extents.
/// `rot_a`/`rot_b` are rotation matrices whose columns are the box's local axes.
///
/// Returns `NoContact` if any axis separates the boxes, otherwise the contact
/// with minimum penetration depth.
#[allow(dead_code)]
pub fn box_box_sat(
    center_a: [f64; 3],
    half_a: [f64; 3],
    rot_a: [[f64; 3]; 3],
    center_b: [f64; 3],
    half_b: [f64; 3],
    rot_b: [[f64; 3]; 3],
) -> NarrowPhaseResult {
    // Translation in world space.
    let t_world = sub(center_b, center_a);

    // Express t in A's frame.
    let t = [
        dot(rot_a[0], t_world),
        dot(rot_a[1], t_world),
        dot(rot_a[2], t_world),
    ];

    // R[i][j] = dot(rot_a[i], rot_b[j])
    let mut r = [[0.0f64; 3]; 3];
    let mut abs_r = [[0.0f64; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            r[i][j] = dot(rot_a[i], rot_b[j]);
            abs_r[i][j] = r[i][j].abs() + 1e-6;
        }
    }

    let mut min_depth = f64::INFINITY;
    let mut best_axis = [0.0f64, 1.0, 0.0];

    // ── A's face normals (axes 0..2) ─────────────────────────────────────────
    for i in 0..3 {
        let ra = half_a[i];
        let rb = half_b[0] * abs_r[i][0] + half_b[1] * abs_r[i][1] + half_b[2] * abs_r[i][2];
        let depth = (ra + rb) - t[i].abs();
        if depth < 0.0 {
            return NarrowPhaseResult::NoContact;
        }
        if depth < min_depth {
            min_depth = depth;
            let sign = if t[i] >= 0.0 { 1.0 } else { -1.0 };
            best_axis = rot_a[i];
            best_axis[0] *= sign;
            best_axis[1] *= sign;
            best_axis[2] *= sign;
        }
    }

    // ── B's face normals (axes 3..5) ─────────────────────────────────────────
    for j in 0..3 {
        let ra = half_a[0] * abs_r[0][j] + half_a[1] * abs_r[1][j] + half_a[2] * abs_r[2][j];
        let rb = half_b[j];
        let sep = t[0] * r[0][j] + t[1] * r[1][j] + t[2] * r[2][j];
        let depth = (ra + rb) - sep.abs();
        if depth < 0.0 {
            return NarrowPhaseResult::NoContact;
        }
        if depth < min_depth {
            min_depth = depth;
            let sign = if sep >= 0.0 { 1.0 } else { -1.0 };
            best_axis = rot_b[j];
            best_axis[0] *= sign;
            best_axis[1] *= sign;
            best_axis[2] *= sign;
        }
    }

    // ── Edge-edge cross products (axes 6..14) ─────────────────────────────────
    let edge_pairs: [(usize, usize); 9] = [
        (0, 0),
        (0, 1),
        (0, 2),
        (1, 0),
        (1, 1),
        (1, 2),
        (2, 0),
        (2, 1),
        (2, 2),
    ];

    for (i, j) in edge_pairs {
        let i1 = (i + 1) % 3;
        let i2 = (i + 2) % 3;
        let j1 = (j + 1) % 3;
        let j2 = (j + 2) % 3;

        let ra = half_a[i1] * abs_r[i2][j] + half_a[i2] * abs_r[i1][j];
        let rb = half_b[j1] * abs_r[i][j2] + half_b[j2] * abs_r[i][j1];
        let sep = t[i2] * r[i1][j] - t[i1] * r[i2][j];
        let depth = (ra + rb) - sep.abs();
        if depth < 0.0 {
            return NarrowPhaseResult::NoContact;
        }

        // Build the cross-product axis in world space.
        let axis = cross(rot_a[i], rot_b[j]);
        let axis_len = len(axis);
        if axis_len > 1e-6 && depth < min_depth {
            min_depth = depth;
            let axis_n = scale(axis, 1.0 / axis_len);
            // Orient toward B.
            let sign = if dot(axis_n, t_world) >= 0.0 {
                1.0
            } else {
                -1.0
            };
            best_axis = scale(axis_n, sign);
        }
    }

    let point_a = add(center_a, scale(best_axis, min_depth * 0.5));
    let point_b = sub(center_b, scale(best_axis, min_depth * 0.5));

    NarrowPhaseResult::SingleContact(Contact {
        point_a,
        point_b,
        normal: best_axis,
        depth: min_depth,
    })
}

// ── Capsule vs Capsule ────────────────────────────────────────────────────────

/// Capsule–capsule narrow-phase test.
///
/// Each capsule is defined by two segment endpoints `(a0, a1)` and radii.
/// Finds the closest points between the two segments, then performs a
/// sphere–sphere test at those points.
#[allow(dead_code)]
pub fn capsule_capsule(
    a0: [f64; 3],
    a1: [f64; 3],
    ra: f64,
    b0: [f64; 3],
    b1: [f64; 3],
    rb: f64,
) -> NarrowPhaseResult {
    let (_, _, ca, cb) = segment_segment_closest(a0, a1, b0, b1);
    sphere_sphere(ca, ra, cb, rb)
}

// ── Closest point utilities ───────────────────────────────────────────────────

/// Closest point on segment `(a, b)` to point `p`.
///
/// Returns the closest point.
#[allow(dead_code)]
pub fn closest_point_on_segment(p: [f64; 3], a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    let ab = sub(b, a);
    let len_sq_ab = len_sq(ab);
    if len_sq_ab < 1e-20 {
        return a;
    }
    let t = (dot(sub(p, a), ab) / len_sq_ab).clamp(0.0, 1.0);
    add(a, scale(ab, t))
}

/// Closest point on triangle `(a, b, c)` to point `p`.
///
/// Uses the Ericson barycentric projection method.
#[allow(dead_code)]
pub fn closest_point_on_triangle(p: [f64; 3], a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> [f64; 3] {
    let ab = sub(b, a);
    let ac = sub(c, a);
    let ap = sub(p, a);

    let d1 = dot(ab, ap);
    let d2 = dot(ac, ap);
    if d1 <= 0.0 && d2 <= 0.0 {
        return a;
    }

    let bp = sub(p, b);
    let d3 = dot(ab, bp);
    let d4 = dot(ac, bp);
    if d3 >= 0.0 && d4 <= d3 {
        return b;
    }

    let vc = d1 * d4 - d3 * d2;
    if vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0 {
        let v = d1 / (d1 - d3);
        return add(a, scale(ab, v));
    }

    let cp = sub(p, c);
    let d5 = dot(ab, cp);
    let d6 = dot(ac, cp);
    if d6 >= 0.0 && d5 <= d6 {
        return c;
    }

    let vb = d5 * d2 - d1 * d6;
    if vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0 {
        let w = d2 / (d2 - d6);
        return add(a, scale(ac, w));
    }

    let va = d3 * d6 - d5 * d4;
    if va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0 {
        let w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        return add(b, scale(sub(c, b), w));
    }

    let denom = 1.0 / (va + vb + vc);
    let v = vb * denom;
    let w = vc * denom;
    add(add(a, scale(ab, v)), scale(ac, w))
}

/// Closest points between two segments `(a0, a1)` and `(b0, b1)`.
///
/// Returns `(t_a, t_b, point_on_a, point_on_b)`.
#[allow(dead_code)]
pub fn segment_segment_closest(
    a0: [f64; 3],
    a1: [f64; 3],
    b0: [f64; 3],
    b1: [f64; 3],
) -> (f64, f64, [f64; 3], [f64; 3]) {
    let d1 = sub(a1, a0);
    let d2 = sub(b1, b0);
    let r = sub(a0, b0);

    let a = dot(d1, d1);
    let e = dot(d2, d2);
    let f = dot(d2, r);

    let (s, t) = if a <= 1e-10 && e <= 1e-10 {
        (0.0_f64, 0.0_f64)
    } else if a <= 1e-10 {
        (0.0, (f / e).clamp(0.0, 1.0))
    } else {
        let c = dot(d1, r);
        if e <= 1e-10 {
            ((-c / a).clamp(0.0, 1.0), 0.0)
        } else {
            let b = dot(d1, d2);
            let denom = a * e - b * b;
            let s0 = if denom.abs() > 1e-10 {
                ((b * f - c * e) / denom).clamp(0.0, 1.0)
            } else {
                0.0
            };
            let t_raw = (b * s0 + f) / e;
            if t_raw < 0.0 {
                ((-c / a).clamp(0.0, 1.0), 0.0)
            } else if t_raw > 1.0 {
                (((b - c) / a).clamp(0.0, 1.0), 1.0)
            } else {
                (s0, t_raw)
            }
        }
    };

    let pa = add(a0, scale(d1, s));
    let pb = add(b0, scale(d2, t));
    (s, t, pa, pb)
}

// ── GJK ───────────────────────────────────────────────────────────────────────

/// A GJK simplex (up to 4 points).
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct GjkSimplex {
    pts: [[f64; 3]; 4],
    count: usize,
}

impl GjkSimplex {
    /// Create an empty simplex.
    #[allow(dead_code)]
    pub fn new() -> Self {
        Self {
            pts: [[0.0; 3]; 4],
            count: 0,
        }
    }

    /// Add a point to the simplex.
    #[allow(dead_code)]
    pub fn push(&mut self, p: [f64; 3]) {
        if self.count < 4 {
            self.pts[self.count] = p;
            self.count += 1;
        }
    }

    /// Number of points in the simplex.
    #[allow(dead_code)]
    pub fn len(&self) -> usize {
        self.count
    }

    /// Returns `true` if the simplex has no points.
    #[allow(dead_code)]
    pub fn is_empty(&self) -> bool {
        self.count == 0
    }
}

impl Default for GjkSimplex {
    fn default() -> Self {
        Self::new()
    }
}

/// Support function for the Minkowski difference of two convex shapes.
///
/// `support_a(d)` returns the support point of shape A in direction `d`.
/// `support_b(d)` returns the support point of shape B in direction `d`.
/// The Minkowski difference support is `support_a(d) - support_b(-d)`.
fn minkowski_diff_support<A, B>(dir: [f64; 3], support_a: &A, support_b: &B) -> [f64; 3]
where
    A: Fn([f64; 3]) -> [f64; 3],
    B: Fn([f64; 3]) -> [f64; 3],
{
    let pa = support_a(dir);
    let pb = support_b(neg(dir));
    sub(pa, pb)
}

/// GJK intersection test for two convex shapes.
///
/// `support_a` and `support_b` are support-function closures.
/// Returns `true` if the shapes intersect (origin inside Minkowski difference).
#[allow(dead_code)]
pub fn gjk_intersect<A, B>(support_a: A, support_b: B) -> bool
where
    A: Fn([f64; 3]) -> [f64; 3],
    B: Fn([f64; 3]) -> [f64; 3],
{
    let mut simplex = GjkSimplex::new();
    let mut dir = [1.0f64, 0.0, 0.0];

    let first = minkowski_diff_support(dir, &support_a, &support_b);
    simplex.push(first);
    dir = neg(first);

    if len_sq(dir) < 1e-20 {
        // Origin is on the first support point — shapes touch.
        return true;
    }

    for _ in 0..64 {
        let a = minkowski_diff_support(dir, &support_a, &support_b);
        if dot(a, dir) < 0.0 {
            // New point doesn't pass origin in direction `dir`.
            return false;
        }
        simplex.push(a);

        match do_simplex_gjk(&mut simplex, &mut dir) {
            GjkStep::ContainsOrigin => return true,
            GjkStep::Continue => {}
        }
    }
    // Conservative: assume intersecting after max iterations.
    true
}

enum GjkStep {
    ContainsOrigin,
    Continue,
}

/// Process simplex and update direction for GJK.
fn do_simplex_gjk(simplex: &mut GjkSimplex, dir: &mut [f64; 3]) -> GjkStep {
    match simplex.count {
        2 => {
            // Line case.
            let b = simplex.pts[0];
            let a = simplex.pts[1];
            let ab = sub(b, a);
            let ao = neg(a);
            if dot(ab, ao) > 0.0 {
                *dir = cross(cross(ab, ao), ab);
            } else {
                simplex.pts[0] = a;
                simplex.count = 1;
                *dir = ao;
            }
            GjkStep::Continue
        }
        3 => {
            // Triangle case.
            let c = simplex.pts[0];
            let b = simplex.pts[1];
            let a = simplex.pts[2];
            let ab = sub(b, a);
            let ac = sub(c, a);
            let ao = neg(a);
            let abc = cross(ab, ac);

            if dot(cross(abc, ac), ao) > 0.0 {
                if dot(ac, ao) > 0.0 {
                    simplex.pts[0] = c;
                    simplex.pts[1] = a;
                    simplex.count = 2;
                    *dir = cross(cross(ac, ao), ac);
                } else {
                    simplex.pts[0] = b;
                    simplex.pts[1] = a;
                    simplex.count = 2;
                    return do_simplex_gjk(simplex, dir);
                }
            } else if dot(cross(ab, abc), ao) > 0.0 {
                simplex.pts[0] = b;
                simplex.pts[1] = a;
                simplex.count = 2;
                return do_simplex_gjk(simplex, dir);
            } else {
                // Origin is above/below triangle plane.
                if dot(abc, ao) > 0.0 {
                    *dir = abc;
                } else {
                    // Flip winding.
                    simplex.pts[0] = b;
                    simplex.pts[1] = c;
                    // pts[2] stays as a
                    *dir = neg(abc);
                }
            }
            GjkStep::Continue
        }
        4 => {
            // Tetrahedron case — check if origin is inside.
            let d = simplex.pts[0];
            let c = simplex.pts[1];
            let b = simplex.pts[2];
            let a = simplex.pts[3];
            let ab = sub(b, a);
            let ac = sub(c, a);
            let ad = sub(d, a);
            let ao = neg(a);

            let abc = cross(ab, ac);
            let acd = cross(ac, ad);
            let adb = cross(ad, ab);

            if dot(abc, ao) > 0.0 {
                // Remove d, keep triangle abc.
                simplex.pts[0] = c;
                simplex.pts[1] = b;
                simplex.pts[2] = a;
                simplex.count = 3;
                return do_simplex_gjk(simplex, dir);
            }
            if dot(acd, ao) > 0.0 {
                simplex.pts[0] = d;
                simplex.pts[1] = c;
                simplex.pts[2] = a;
                simplex.count = 3;
                return do_simplex_gjk(simplex, dir);
            }
            if dot(adb, ao) > 0.0 {
                simplex.pts[0] = b;
                simplex.pts[1] = d;
                simplex.pts[2] = a;
                simplex.count = 3;
                return do_simplex_gjk(simplex, dir);
            }
            GjkStep::ContainsOrigin
        }
        _ => GjkStep::Continue,
    }
}

// ── Contact manifold clip (box-box reference face) ────────────────────────────

/// Clip a polygon against one side-plane of a reference face.
///
/// `poly` is the current polygon vertices.
/// `plane_normal` and `plane_offset` define the clip plane (`n·x = offset`).
/// Points with `n·x > offset` are kept (or interpolated at the boundary).
#[allow(dead_code)]
fn clip_polygon_by_plane(
    poly: &[[f64; 3]],
    plane_normal: [f64; 3],
    plane_offset: f64,
) -> Vec<[f64; 3]> {
    let mut out = Vec::new();
    if poly.is_empty() {
        return out;
    }
    let n = poly.len();
    for i in 0..n {
        let curr = poly[i];
        let next = poly[(i + 1) % n];
        let dc = dot(plane_normal, curr) - plane_offset;
        let dn = dot(plane_normal, next) - plane_offset;

        if dc >= 0.0 {
            out.push(curr);
        }
        // Edge crosses boundary.
        if (dc < 0.0) != (dn < 0.0) {
            let t = dc / (dc - dn);
            let interp = add(curr, scale(sub(next, curr), t));
            out.push(interp);
        }
    }
    out
}

/// Build a box-box contact manifold by clipping the incident face polygon
/// against the reference face's four side planes.
///
/// `ref_center` is the world-space centre of the reference face.
/// `ref_u` and `ref_v` are the two tangent axes of the reference face.
/// `ref_half_u` and `ref_half_v` are the face's half-extents along those axes.
/// `ref_normal` is the face normal (pointing outward from reference box).
/// `incident_verts` are the four corners of the incident face.
/// Returns a vector of contact points (those below the reference face plane).
#[allow(dead_code)]
pub fn clip_incident_face(
    ref_center: [f64; 3],
    ref_u: [f64; 3],
    ref_v: [f64; 3],
    ref_half_u: f64,
    ref_half_v: f64,
    ref_normal: [f64; 3],
    incident_verts: &[[f64; 3]; 4],
) -> Vec<Contact> {
    // Start with incident face polygon.
    let mut poly: Vec<[f64; 3]> = incident_verts.to_vec();

    // Clip by the four edge planes of the reference face.
    poly = clip_polygon_by_plane(&poly, ref_u, dot(ref_u, ref_center) + ref_half_u);
    poly = clip_polygon_by_plane(&poly, neg(ref_u), -dot(ref_u, ref_center) + ref_half_u);
    poly = clip_polygon_by_plane(&poly, ref_v, dot(ref_v, ref_center) + ref_half_v);
    poly = clip_polygon_by_plane(&poly, neg(ref_v), -dot(ref_v, ref_center) + ref_half_v);

    // Keep points below (or on) the reference face plane.
    let ref_plane_offset = dot(ref_normal, ref_center);
    poly.retain(|&p| dot(ref_normal, p) <= ref_plane_offset + 1e-6);

    poly.into_iter()
        .map(|pt| {
            let proj = add(
                ref_center,
                add(
                    scale(ref_u, dot(ref_u, sub(pt, ref_center))),
                    scale(ref_v, dot(ref_v, sub(pt, ref_center))),
                ),
            );
            let depth = ref_plane_offset - dot(ref_normal, pt);
            Contact {
                point_a: proj,
                point_b: pt,
                normal: ref_normal,
                depth: depth.max(0.0),
            }
        })
        .collect()
}

// ── Triangle-AABB intersection ─────────────────────────────────────────────────

/// Triangle-AABB intersection test using the Separating Axis Theorem.
///
/// `tri_verts` – the three triangle vertices in world space.
/// `box_center` – centre of the AABB.
/// `half_extents` – half-extents of the AABB along x, y, z.
///
/// Returns `true` if the triangle overlaps the box.
#[allow(dead_code)]
pub fn triangle_aabb_intersect(
    tri_verts: &[[f64; 3]; 3],
    box_center: [f64; 3],
    half_extents: [f64; 3],
) -> bool {
    // Translate triangle to box-centred space.
    let v: [[f64; 3]; 3] = [
        sub(tri_verts[0], box_center),
        sub(tri_verts[1], box_center),
        sub(tri_verts[2], box_center),
    ];

    // Test 3 face-normal axes of the AABB.
    for a in 0..3 {
        let mn = v[0][a].min(v[1][a]).min(v[2][a]);
        let mx = v[0][a].max(v[1][a]).max(v[2][a]);
        if mn > half_extents[a] || mx < -half_extents[a] {
            return false;
        }
    }

    // Test triangle normal axis.
    let e0 = sub(v[1], v[0]);
    let e1 = sub(v[2], v[1]);
    let n_tri = cross(e0, e1);
    let d = dot(n_tri, v[0]);
    // Radius of AABB projected onto triangle normal.
    let r = half_extents[0] * n_tri[0].abs()
        + half_extents[1] * n_tri[1].abs()
        + half_extents[2] * n_tri[2].abs();
    if d.abs() > r {
        return false;
    }

    // Test 9 cross-product axes (edge × box-axis).
    let edges = [e0, e1, sub(v[0], v[2])];
    let box_axes: [[f64; 3]; 3] = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];

    for &edge in &edges {
        for &axis in &box_axes {
            let sep_axis = cross(edge, axis);
            if len_sq(sep_axis) < 1e-20 {
                continue;
            }
            let p: [f64; 3] = [
                dot(sep_axis, v[0]),
                dot(sep_axis, v[1]),
                dot(sep_axis, v[2]),
            ];
            let p_min = p[0].min(p[1]).min(p[2]);
            let p_max = p[0].max(p[1]).max(p[2]);
            let r_box = half_extents[0] * sep_axis[0].abs()
                + half_extents[1] * sep_axis[1].abs()
                + half_extents[2] * sep_axis[2].abs();
            if p_min > r_box || p_max < -r_box {
                return false;
            }
        }
    }

    true
}

// ── OBB-Capsule intersection ───────────────────────────────────────────────────

/// OBB–capsule narrow-phase test.
///
/// The capsule is defined by its segment `(c0, c1)` and radius `cr`.
/// The OBB is defined by `center`, `half_extents`, and rotation matrix `rot`.
/// Returns `NoContact` or a `SingleContact`.
#[allow(dead_code)]
pub fn obb_capsule(
    obb_center: [f64; 3],
    obb_half: [f64; 3],
    obb_rot: [[f64; 3]; 3],
    cap_a: [f64; 3],
    cap_b: [f64; 3],
    cap_r: f64,
) -> NarrowPhaseResult {
    // Transform capsule segment endpoints into OBB-local space.
    let local_a = world_to_obb_local(sub(cap_a, obb_center), obb_rot);
    let local_b = world_to_obb_local(sub(cap_b, obb_center), obb_rot);

    // Find closest point on the segment to the AABB (OBB in local = AABB).
    let cp = closest_point_segment_aabb(local_a, local_b, obb_half);

    // Sphere test at closest point.
    let diff = sub(local_a, cp);
    // The actual closest point on the segment to the AABB surface needs finding.
    // We use: closest segment point to the AABB closest point.
    let cp_seg = closest_point_on_segment(cp, local_a, local_b);
    let d = sub(cp_seg, cp);
    let dist = len(d);
    let depth = cap_r - dist;
    if depth < 0.0 {
        return NarrowPhaseResult::NoContact;
    }

    let _ = diff; // silence unused warning
    let normal_local = if dist > 1e-10 {
        normalize(d)
    } else {
        [0.0, 1.0, 0.0]
    };
    // Transform normal back to world space.
    let normal_world = obb_local_to_world(normal_local, obb_rot);

    let pa = add(
        obb_center,
        obb_local_to_world(add(cp, scale(normal_local, depth * 0.5)), obb_rot),
    );
    let pb = add(
        obb_center,
        obb_local_to_world(
            sub(cp_seg, scale(normal_local, cap_r - depth * 0.5)),
            obb_rot,
        ),
    );

    NarrowPhaseResult::SingleContact(Contact {
        point_a: pa,
        point_b: pb,
        normal: normal_world,
        depth,
    })
}

/// Closest point on segment `(a, b)` to a box `[-half, +half]` in 3D.
fn closest_point_segment_aabb(a: [f64; 3], b: [f64; 3], half: [f64; 3]) -> [f64; 3] {
    // Clamp parametric t to find closest segment point to AABB.
    // We find the closest point on the box to the segment endpoints and midpoint.
    let candidates = [a, b, add(scale(a, 0.5), scale(b, 0.5))];
    let mut best = [0.0f64; 3];
    let mut best_dist = f64::INFINITY;
    for &p in &candidates {
        let clamped = [
            p[0].clamp(-half[0], half[0]),
            p[1].clamp(-half[1], half[1]),
            p[2].clamp(-half[2], half[2]),
        ];
        let d = len_sq(sub(p, clamped));
        if d < best_dist {
            best_dist = d;
            best = clamped;
        }
    }
    best
}

fn world_to_obb_local(v: [f64; 3], rot: [[f64; 3]; 3]) -> [f64; 3] {
    [dot(rot[0], v), dot(rot[1], v), dot(rot[2], v)]
}

fn obb_local_to_world(v: [f64; 3], rot: [[f64; 3]; 3]) -> [f64; 3] {
    [
        rot[0][0] * v[0] + rot[1][0] * v[1] + rot[2][0] * v[2],
        rot[0][1] * v[0] + rot[1][1] * v[1] + rot[2][1] * v[2],
        rot[0][2] * v[0] + rot[1][2] * v[1] + rot[2][2] * v[2],
    ]
}

// ── AABB-Capsule intersection ──────────────────────────────────────────────────

/// AABB–capsule narrow-phase test.
///
/// Finds the closest point on the capsule segment to the AABB, then performs
/// a sphere–box test at that point.
#[allow(dead_code)]
pub fn aabb_capsule(
    box_center: [f64; 3],
    half_extents: [f64; 3],
    cap_a: [f64; 3],
    cap_b: [f64; 3],
    cap_r: f64,
) -> NarrowPhaseResult {
    // Closest point on segment to box centre.
    let closest_on_seg = closest_point_on_segment(box_center, cap_a, cap_b);
    // Now test the sphere at that point against the AABB.
    sphere_box(closest_on_seg, cap_r, box_center, half_extents)
}

// ── Heightfield-Sphere ────────────────────────────────────────────────────────

/// A simple heightfield defined on a uniform 2D grid.
///
/// The heightfield stores heights at grid points (x, z) with spacing `cell_size`.
/// Origin is at `(origin_x, origin_z)`.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct Heightfield {
    /// Heights at grid points (row-major: row = z index, col = x index).
    pub heights: Vec<f64>,
    /// Number of columns (x direction).
    pub n_x: usize,
    /// Number of rows (z direction).
    pub n_z: usize,
    /// Cell size (same in x and z).
    pub cell_size: f64,
    /// World-space x origin.
    pub origin_x: f64,
    /// World-space z origin.
    pub origin_z: f64,
}

impl Heightfield {
    /// Create a flat heightfield at `y = 0`.
    #[allow(dead_code)]
    pub fn flat(n_x: usize, n_z: usize, cell_size: f64) -> Self {
        Self {
            heights: vec![0.0; n_x * n_z],
            n_x,
            n_z,
            cell_size,
            origin_x: 0.0,
            origin_z: 0.0,
        }
    }

    /// Get the height at grid cell (ix, iz).
    #[allow(dead_code)]
    pub fn height_at(&self, ix: usize, iz: usize) -> f64 {
        self.heights[iz * self.n_x + ix]
    }

    /// Bilinearly interpolated height at world-space (x, z).
    #[allow(dead_code)]
    pub fn height_world(&self, x: f64, z: f64) -> f64 {
        let lx = (x - self.origin_x) / self.cell_size;
        let lz = (z - self.origin_z) / self.cell_size;
        let ix = (lx as usize).min(self.n_x - 2);
        let iz = (lz as usize).min(self.n_z - 2);
        let fx = (lx - ix as f64).clamp(0.0, 1.0);
        let fz = (lz - iz as f64).clamp(0.0, 1.0);

        let h00 = self.height_at(ix, iz);
        let h10 = self.height_at(ix + 1, iz);
        let h01 = self.height_at(ix, iz + 1);
        let h11 = self.height_at(ix + 1, iz + 1);

        h00 * (1.0 - fx) * (1.0 - fz)
            + h10 * fx * (1.0 - fz)
            + h01 * (1.0 - fx) * fz
            + h11 * fx * fz
    }

    /// Surface normal at world-space (x, z) via central differences.
    #[allow(dead_code)]
    pub fn normal_world(&self, x: f64, z: f64) -> [f64; 3] {
        let h = self.cell_size * 0.5;
        let dydx = (self.height_world(x + h, z) - self.height_world(x - h, z)) / (2.0 * h);
        let dydz = (self.height_world(x, z + h) - self.height_world(x, z - h)) / (2.0 * h);
        normalize([-dydx, 1.0, -dydz])
    }
}

/// Heightfield–sphere narrow-phase test.
///
/// Finds the heightfield height directly beneath the sphere centre and checks
/// for penetration.
#[allow(dead_code)]
pub fn heightfield_sphere(
    hf: &Heightfield,
    sphere_center: [f64; 3],
    radius: f64,
) -> NarrowPhaseResult {
    let x = sphere_center[0];
    let z = sphere_center[2];
    let y_surf = hf.height_world(x, z);
    let depth = radius - (sphere_center[1] - y_surf);

    if depth < 0.0 {
        return NarrowPhaseResult::NoContact;
    }

    let normal = hf.normal_world(x, z);
    let point_b = [x, y_surf, z];
    let point_a = sub(sphere_center, scale(normal, radius));

    NarrowPhaseResult::SingleContact(Contact {
        point_a,
        point_b,
        normal,
        depth,
    })
}

// ── Heightfield-Box ───────────────────────────────────────────────────────────

/// Heightfield–AABB narrow-phase test.
///
/// Iterates over the grid cells beneath the box AABB footprint and tests each
/// cell's representative height against the box bottom face.
#[allow(dead_code)]
pub fn heightfield_box(
    hf: &Heightfield,
    box_center: [f64; 3],
    half_extents: [f64; 3],
) -> NarrowPhaseResult {
    let x_min = box_center[0] - half_extents[0];
    let x_max = box_center[0] + half_extents[0];
    let z_min = box_center[2] - half_extents[2];
    let z_max = box_center[2] + half_extents[2];

    let mut deepest: Option<Contact> = None;
    let mut max_depth = 0.0f64;

    // Sample at 4 corners and centre of the box footprint.
    let samples = [
        [x_min, z_min],
        [x_max, z_min],
        [x_min, z_max],
        [x_max, z_max],
        [box_center[0], box_center[2]],
    ];

    for &[sx, sz] in &samples {
        let y_surf = hf.height_world(sx, sz);
        let box_bottom = box_center[1] - half_extents[1];
        let depth = y_surf - box_bottom;
        if depth > max_depth {
            max_depth = depth;
            let normal = hf.normal_world(sx, sz);
            deepest = Some(Contact {
                point_a: [sx, box_bottom, sz],
                point_b: [sx, y_surf, sz],
                normal,
                depth,
            });
        }
    }

    match deepest {
        Some(c) => NarrowPhaseResult::SingleContact(c),
        None => NarrowPhaseResult::NoContact,
    }
}

// ── Triangle–sphere narrow phase ──────────────────────────────────────────────

/// Triangle–sphere narrow-phase test.
///
/// Finds the closest point on the triangle to the sphere centre and performs
/// a point-sphere overlap check.
#[allow(dead_code)]
pub fn triangle_sphere(
    tri_a: [f64; 3],
    tri_b: [f64; 3],
    tri_c: [f64; 3],
    sphere_center: [f64; 3],
    radius: f64,
) -> NarrowPhaseResult {
    let cp = closest_point_on_triangle(sphere_center, tri_a, tri_b, tri_c);
    let d = sub(sphere_center, cp);
    let dist = len(d);
    let depth = radius - dist;

    if depth < 0.0 {
        return NarrowPhaseResult::NoContact;
    }

    let normal = if dist > 1e-10 {
        scale(d, 1.0 / dist)
    } else {
        // Degenerate: use triangle normal.
        let ab = sub(tri_b, tri_a);
        let ac = sub(tri_c, tri_a);
        normalize(cross(ab, ac))
    };

    let point_a = sub(sphere_center, scale(normal, radius));
    NarrowPhaseResult::SingleContact(Contact {
        point_a,
        point_b: cp,
        normal,
        depth,
    })
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    fn approx_eq(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }

    fn vec_approx_eq(a: [f64; 3], b: [f64; 3], tol: f64) -> bool {
        approx_eq(a[0], b[0], tol) && approx_eq(a[1], b[1], tol) && approx_eq(a[2], b[2], tol)
    }

    // ── sphere_sphere ────────────────────────────────────────────────────────

    #[test]
    fn test_sphere_sphere_touching() {
        let ca = [0.0f64, 0.0, 0.0];
        let cb = [2.0f64, 0.0, 0.0];
        match sphere_sphere(ca, 1.0, cb, 1.0) {
            NarrowPhaseResult::SingleContact(c) => {
                assert!(approx_eq(c.depth, 0.0, 1e-10), "depth={}", c.depth);
            }
            _ => panic!("expected SingleContact for touching spheres"),
        }
    }

    #[test]
    fn test_sphere_sphere_overlapping() {
        let ca = [0.0f64, 0.0, 0.0];
        let cb = [1.5f64, 0.0, 0.0];
        match sphere_sphere(ca, 1.0, cb, 1.0) {
            NarrowPhaseResult::SingleContact(c) => {
                assert!(approx_eq(c.depth, 0.5, 1e-10), "depth={}", c.depth);
                let nlen = len(c.normal);
                assert!(approx_eq(nlen, 1.0, 1e-10), "normal length={}", nlen);
            }
            _ => panic!("expected SingleContact for overlapping spheres"),
        }
    }

    #[test]
    fn test_sphere_sphere_separated() {
        let ca = [0.0f64, 0.0, 0.0];
        let cb = [5.0f64, 0.0, 0.0];
        assert!(
            matches!(
                sphere_sphere(ca, 1.0, cb, 1.0),
                NarrowPhaseResult::NoContact
            ),
            "separated spheres must return NoContact"
        );
    }

    // ── sphere_box ───────────────────────────────────────────────────────────

    #[test]
    fn test_sphere_box_hit() {
        // Sphere center at [1.3, 0, 0], radius 0.5; box center at [0,0,0] half_extents [1,1,1].
        // Closest point on box surface = [1, 0, 0]; dist = 0.3 < radius → contact.
        let sc = [1.3f64, 0.0, 0.0];
        let bc = [0.0f64, 0.0, 0.0];
        let he = [1.0f64, 1.0, 1.0];
        match sphere_box(sc, 0.5, bc, he) {
            NarrowPhaseResult::SingleContact(c) => {
                assert!(c.depth > 0.0, "depth={}", c.depth);
                assert!(approx_eq(len(c.normal), 1.0, 1e-10));
            }
            _ => panic!("expected contact"),
        }
    }

    #[test]
    fn test_sphere_box_miss() {
        let sc = [5.0f64, 0.0, 0.0];
        let bc = [0.0f64, 0.0, 0.0];
        let he = [1.0f64, 1.0, 1.0];
        assert!(
            matches!(sphere_box(sc, 0.5, bc, he), NarrowPhaseResult::NoContact),
            "far sphere must return NoContact"
        );
    }

    // ── box_box_sat ──────────────────────────────────────────────────────────

    #[test]
    fn test_box_box_axis_aligned_overlapping() {
        let identity = [[1.0f64, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let ca = [0.0f64, 0.0, 0.0];
        let cb = [1.5f64, 0.0, 0.0];
        let half = [1.0f64, 1.0, 1.0];
        match box_box_sat(ca, half, identity, cb, half, identity) {
            NarrowPhaseResult::SingleContact(c) => {
                assert!(c.depth > 0.0, "depth={}", c.depth);
                assert!(
                    approx_eq(c.depth, 0.5, 1e-5),
                    "expected depth≈0.5, got {}",
                    c.depth
                );
            }
            _ => panic!("expected contact for overlapping boxes"),
        }
    }

    #[test]
    fn test_box_box_axis_aligned_separated() {
        let identity = [[1.0f64, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let ca = [0.0f64, 0.0, 0.0];
        let cb = [5.0f64, 0.0, 0.0];
        let half = [1.0f64, 1.0, 1.0];
        assert!(
            matches!(
                box_box_sat(ca, half, identity, cb, half, identity),
                NarrowPhaseResult::NoContact
            ),
            "separated boxes must return NoContact"
        );
    }

    // ── capsule_capsule ──────────────────────────────────────────────────────

    #[test]
    fn test_capsule_capsule_overlap() {
        let a0 = [0.0f64, -1.0, 0.0];
        let a1 = [0.0f64, 1.0, 0.0];
        let b0 = [0.8f64, -1.0, 0.0];
        let b1 = [0.8f64, 1.0, 0.0];
        match capsule_capsule(a0, a1, 0.5, b0, b1, 0.5) {
            NarrowPhaseResult::SingleContact(c) => {
                assert!(c.depth > 0.0, "depth={}", c.depth);
                assert!(approx_eq(len(c.normal), 1.0, 1e-10));
            }
            _ => panic!("expected contact for overlapping capsules"),
        }
    }

    #[test]
    fn test_capsule_capsule_separated() {
        let a0 = [0.0f64, 0.0, 0.0];
        let a1 = [1.0f64, 0.0, 0.0];
        let b0 = [10.0f64, 0.0, 0.0];
        let b1 = [11.0f64, 0.0, 0.0];
        assert!(
            matches!(
                capsule_capsule(a0, a1, 0.3, b0, b1, 0.3),
                NarrowPhaseResult::NoContact
            ),
            "far capsules must return NoContact"
        );
    }

    // ── closest_point_on_segment ─────────────────────────────────────────────

    #[test]
    fn test_closest_point_on_segment_midpoint() {
        let p = [1.0f64, 1.0, 0.0];
        let a = [0.0f64, 0.0, 0.0];
        let b = [2.0f64, 0.0, 0.0];
        let cp = closest_point_on_segment(p, a, b);
        assert!(vec_approx_eq(cp, [1.0, 0.0, 0.0], 1e-10), "cp={:?}", cp);
    }

    #[test]
    fn test_closest_point_on_segment_endpoint() {
        let p = [-1.0f64, 0.0, 0.0];
        let a = [0.0f64, 0.0, 0.0];
        let b = [2.0f64, 0.0, 0.0];
        let cp = closest_point_on_segment(p, a, b);
        assert!(vec_approx_eq(cp, [0.0, 0.0, 0.0], 1e-10), "cp={:?}", cp);
    }

    #[test]
    fn test_closest_point_on_triangle_vertex() {
        // Point behind vertex a — should clamp to a.
        let p = [-1.0f64, -1.0, 0.0];
        let a = [0.0f64, 0.0, 0.0];
        let b = [2.0f64, 0.0, 0.0];
        let c = [1.0f64, 2.0, 0.0];
        let cp = closest_point_on_triangle(p, a, b, c);
        assert!(vec_approx_eq(cp, a, 1e-10), "cp={:?}", cp);
    }

    // ── gjk_intersect ────────────────────────────────────────────────────────

    #[test]
    fn test_gjk_sphere_pair_overlap() {
        // Two unit spheres at [0,0,0] and [1,0,0]: they overlap.
        let support_a = |d: [f64; 3]| -> [f64; 3] {
            let n = normalize(d);
            scale(n, 1.0) // sphere A at origin, radius 1
        };
        let support_b = |d: [f64; 3]| -> [f64; 3] {
            let n = normalize(d);
            add([1.0, 0.0, 0.0], scale(n, 1.0)) // sphere B at [1,0,0], radius 1
        };
        assert!(
            gjk_intersect(support_a, support_b),
            "overlapping spheres must intersect"
        );
    }

    #[test]
    fn test_gjk_sphere_pair_separated() {
        // Two unit spheres at [0,0,0] and [5,0,0]: separated.
        let support_a = |d: [f64; 3]| -> [f64; 3] { scale(normalize(d), 1.0) };
        let support_b =
            |d: [f64; 3]| -> [f64; 3] { add([5.0, 0.0, 0.0], scale(normalize(d), 1.0)) };
        assert!(
            !gjk_intersect(support_a, support_b),
            "separated spheres must not intersect"
        );
    }

    // ── NarrowPhaseResult helpers ────────────────────────────────────────────

    #[test]
    fn test_narrow_phase_result_has_contact() {
        assert!(!NarrowPhaseResult::NoContact.has_contact());
        let c = Contact {
            point_a: [0.0; 3],
            point_b: [0.0; 3],
            normal: [0.0, 1.0, 0.0],
            depth: 0.1,
        };
        assert!(NarrowPhaseResult::SingleContact(c.clone()).has_contact());
        assert!(NarrowPhaseResult::Manifold(vec![c]).has_contact());
    }

    // ── triangle_aabb_intersect ─────────────────────────────────────────────

    #[test]
    fn test_triangle_aabb_intersect_overlap() {
        // Triangle in xy-plane centred at origin, box also centred at origin.
        let tri: [[f64; 3]; 3] = [[-0.5, 0.0, -0.5], [0.5, 0.0, -0.5], [0.0, 0.0, 0.5]];
        let box_c = [0.0f64, 0.0, 0.0];
        let half = [1.0f64, 1.0, 1.0];
        assert!(
            triangle_aabb_intersect(&tri, box_c, half),
            "Triangle inside box must intersect"
        );
    }

    #[test]
    fn test_triangle_aabb_intersect_separated() {
        // Triangle far away from the box.
        let tri: [[f64; 3]; 3] = [[10.0, 0.0, 0.0], [11.0, 0.0, 0.0], [10.5, 1.0, 0.0]];
        let box_c = [0.0f64, 0.0, 0.0];
        let half = [1.0f64, 1.0, 1.0];
        assert!(
            !triangle_aabb_intersect(&tri, box_c, half),
            "Separated triangle must not intersect box"
        );
    }

    #[test]
    fn test_triangle_aabb_intersect_edge_touch() {
        // Triangle with one vertex just touching the box face.
        let tri: [[f64; 3]; 3] = [
            [1.0, 0.0, 0.0], // on +x face
            [2.0, 0.5, 0.0],
            [2.0, -0.5, 0.0],
        ];
        let box_c = [0.0f64, 0.0, 0.0];
        let half = [1.0f64, 1.0, 1.0];
        // The vertex at x=1 is exactly on the face; should intersect.
        assert!(
            triangle_aabb_intersect(&tri, box_c, half),
            "Triangle touching box face must intersect"
        );
    }

    // ── aabb_capsule ────────────────────────────────────────────────────────

    #[test]
    fn test_aabb_capsule_overlap() {
        // Capsule along y-axis from (0,-2,0) to (0,2,0) with radius 0.5.
        // Box at origin with half-extents [1,1,1] → capsule passes through → contact.
        let box_c = [0.0f64, 0.0, 0.0];
        let half = [1.0f64, 1.0, 1.0];
        let cap_a = [0.0f64, -2.0, 0.0];
        let cap_b = [0.0f64, 2.0, 0.0];
        match aabb_capsule(box_c, half, cap_a, cap_b, 0.5) {
            NarrowPhaseResult::SingleContact(c) => {
                assert!(c.depth >= 0.0, "depth={}", c.depth);
            }
            _ => panic!("expected contact for overlapping AABB-capsule"),
        }
    }

    #[test]
    fn test_aabb_capsule_separated() {
        let box_c = [0.0f64, 0.0, 0.0];
        let half = [0.5f64, 0.5, 0.5];
        let cap_a = [5.0f64, 0.0, 0.0];
        let cap_b = [6.0f64, 0.0, 0.0];
        assert!(
            matches!(
                aabb_capsule(box_c, half, cap_a, cap_b, 0.3),
                NarrowPhaseResult::NoContact
            ),
            "Far capsule must return NoContact"
        );
    }

    // ── obb_capsule ─────────────────────────────────────────────────────────

    #[test]
    fn test_obb_capsule_axis_aligned_overlap() {
        let identity = [[1.0f64, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let obb_c = [0.0f64, 0.0, 0.0];
        let obb_h = [1.0f64, 1.0, 1.0];
        let cap_a = [0.0f64, -2.0, 0.0];
        let cap_b = [0.0f64, 2.0, 0.0];
        let result = obb_capsule(obb_c, obb_h, identity, cap_a, cap_b, 0.5);
        assert!(result.has_contact(), "OBB-capsule should have contact");
    }

    // ── heightfield_sphere ──────────────────────────────────────────────────

    #[test]
    fn test_heightfield_sphere_above() {
        let hf = Heightfield::flat(5, 5, 1.0);
        // Sphere far above the flat heightfield → no contact.
        let sphere_c = [2.0f64, 5.0, 2.0];
        assert!(
            matches!(
                heightfield_sphere(&hf, sphere_c, 0.5),
                NarrowPhaseResult::NoContact
            ),
            "Sphere high above heightfield must return NoContact"
        );
    }

    #[test]
    fn test_heightfield_sphere_below() {
        let hf = Heightfield::flat(5, 5, 1.0);
        // Sphere centre at y=0.3 with radius=0.5: bottom is at y=-0.2 → penetration = 0.5-0.3 = 0.2.
        let sphere_c = [2.0f64, 0.3, 2.0];
        match heightfield_sphere(&hf, sphere_c, 0.5) {
            NarrowPhaseResult::SingleContact(c) => {
                assert!(
                    approx_eq(c.depth, 0.2, 1e-6),
                    "Expected depth 0.2, got {}",
                    c.depth
                );
            }
            _ => panic!("expected contact for sphere below heightfield"),
        }
    }

    #[test]
    fn test_heightfield_sphere_contact_normal() {
        let hf = Heightfield::flat(5, 5, 1.0);
        let sphere_c = [2.0f64, 0.3, 2.0];
        match heightfield_sphere(&hf, sphere_c, 0.5) {
            NarrowPhaseResult::SingleContact(c) => {
                // Normal should point up for a flat surface.
                assert!(
                    approx_eq(c.normal[1], 1.0, 1e-6),
                    "Normal should be (0,1,0) for flat surface: {:?}",
                    c.normal
                );
            }
            _ => panic!("expected contact"),
        }
    }

    // ── heightfield_box ─────────────────────────────────────────────────────

    #[test]
    fn test_heightfield_box_above() {
        let hf = Heightfield::flat(5, 5, 1.0);
        let box_c = [2.0f64, 5.0, 2.0];
        let half = [0.5f64, 0.5, 0.5];
        assert!(
            matches!(
                heightfield_box(&hf, box_c, half),
                NarrowPhaseResult::NoContact
            ),
            "Box far above heightfield must return NoContact"
        );
    }

    #[test]
    fn test_heightfield_box_penetrating() {
        let hf = Heightfield::flat(5, 5, 1.0);
        // Box bottom at y = 0.0 - 0.5 = -0.5 → penetrates flat surface at y=0.
        let box_c = [2.0f64, 0.0, 2.0];
        let half = [0.5f64, 0.5, 0.5];
        match heightfield_box(&hf, box_c, half) {
            NarrowPhaseResult::SingleContact(c) => {
                assert!(
                    c.depth > 0.0,
                    "Penetrating box should have positive depth: {}",
                    c.depth
                );
            }
            _ => panic!("expected contact for penetrating box"),
        }
    }

    // ── triangle_sphere ──────────────────────────────────────────────────────

    #[test]
    fn test_triangle_sphere_contact() {
        // Sphere centred 0.3 above the triangle → contact with depth = radius - 0.3.
        let tri_a = [-1.0f64, 0.0, -1.0];
        let tri_b = [1.0f64, 0.0, -1.0];
        let tri_c = [0.0f64, 0.0, 1.0];
        let sphere_c = [0.0f64, 0.3, 0.0];
        let radius = 0.5;
        match triangle_sphere(tri_a, tri_b, tri_c, sphere_c, radius) {
            NarrowPhaseResult::SingleContact(c) => {
                assert!(approx_eq(c.depth, 0.2, 1e-6), "depth={}", c.depth);
                assert!(
                    approx_eq(len(c.normal), 1.0, 1e-10),
                    "normal not unit: {:?}",
                    c.normal
                );
            }
            _ => panic!("expected contact"),
        }
    }

    #[test]
    fn test_triangle_sphere_no_contact() {
        let tri_a = [-1.0f64, 0.0, -1.0];
        let tri_b = [1.0f64, 0.0, -1.0];
        let tri_c = [0.0f64, 0.0, 1.0];
        let sphere_c = [0.0f64, 5.0, 0.0];
        assert!(
            matches!(
                triangle_sphere(tri_a, tri_b, tri_c, sphere_c, 0.5),
                NarrowPhaseResult::NoContact
            ),
            "Sphere far above triangle must return NoContact"
        );
    }

    #[test]
    fn test_heightfield_bilinear_flat() {
        let hf = Heightfield::flat(4, 4, 1.0);
        // For a flat heightfield, height at any point should be 0.
        assert!(approx_eq(hf.height_world(1.5, 1.5), 0.0, 1e-12));
        assert!(approx_eq(hf.height_world(0.1, 2.9), 0.0, 1e-12));
    }

    #[test]
    fn test_heightfield_normal_flat_is_up() {
        let hf = Heightfield::flat(4, 4, 1.0);
        let n = hf.normal_world(1.5, 1.5);
        assert!(
            vec_approx_eq(n, [0.0, 1.0, 0.0], 1e-6),
            "Flat HF normal: {:?}",
            n
        );
    }
}
