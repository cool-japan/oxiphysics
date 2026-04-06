//! Contact manifold clipping for narrow-phase collision detection.
//!
//! Implements Sutherland-Hodgman polygon clipping in 2D and 3D, box-box
//! manifold construction via reference/incident face identification and
//! polygon clipping, and manifold reduction to a bounded contact set.

#![allow(clippy::needless_range_loop)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0
#![allow(dead_code)]

#[allow(unused_imports)]
use super::functions::*;

// ─── ContactPoint ────────────────────────────────────────────────────────────

/// A single contact point produced by narrow-phase collision detection.
#[derive(Debug, Clone)]
pub struct ContactPoint {
    /// World-space position of the contact.
    pub point: [f64; 3],
    /// Contact normal (unit vector, from body B toward body A).
    pub normal: [f64; 3],
    /// Penetration depth (positive when bodies overlap).
    pub depth: f64,
}

// ─── ContactManifold ─────────────────────────────────────────────────────────

/// A set of contact points sharing a common contact normal.
#[derive(Debug, Clone)]
pub struct ContactManifold {
    /// Individual contact points in this manifold.
    pub contacts: Vec<ContactPoint>,
    /// Shared contact normal for the manifold.
    pub normal: [f64; 3],
}

impl ContactManifold {
    /// Create an empty manifold with the given normal.
    pub fn new(normal: [f64; 3]) -> Self {
        Self {
            contacts: Vec::new(),
            normal,
        }
    }
}

// ─── Math helpers ─────────────────────────────────────────────────────────────

#[inline]
fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[inline]
fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

#[inline]
fn add3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

#[inline]
fn scale3(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

#[inline]
fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

#[inline]
fn norm3(a: [f64; 3]) -> f64 {
    dot3(a, a).sqrt()
}

#[inline]
fn normalize3(a: [f64; 3]) -> [f64; 3] {
    let n = norm3(a);
    if n < 1e-15 { a } else { scale3(a, 1.0 / n) }
}

// ─── Sutherland-Hodgman 2-D ───────────────────────────────────────────────────

/// Clip a 2-D `subject` polygon against a convex `clip` polygon using the
/// Sutherland–Hodgman algorithm.
///
/// Both polygons are given as ordered vertex lists (clockwise or
/// counter-clockwise, but both in the same winding order). An empty `clip`
/// polygon returns an empty result.
pub fn sutherland_hodgman_2d(subject: &[[f64; 2]], clip: &[[f64; 2]]) -> Vec<[f64; 2]> {
    if clip.is_empty() || subject.is_empty() {
        return Vec::new();
    }

    let mut output: Vec<[f64; 2]> = subject.to_vec();

    let n = clip.len();
    for i in 0..n {
        if output.is_empty() {
            return Vec::new();
        }
        let input = output.clone();
        output.clear();

        let edge_start = clip[i];
        let edge_end = clip[(i + 1) % n];

        // Edge direction and its left-hand normal (inside half-space normal).
        let ex = edge_end[0] - edge_start[0];
        let ey = edge_end[1] - edge_start[1];

        // inside(p): positive when p is on the left (inside) of the directed edge.
        let inside = |p: [f64; 2]| -> bool {
            let dx = p[0] - edge_start[0];
            let dy = p[1] - edge_start[1];
            // cross product (edge × dp); positive → left side
            ex * dy - ey * dx >= 0.0
        };

        // Intersection of segment S→E with the current clip edge.
        // We solve for t in: S + t*(E-S) = edge_start + u*(edge_end-edge_start)
        // Using 2-D cross products:
        //   (E-S) × (edge_end-edge_start) * t = (edge_start-S) × (edge_end-edge_start)
        let intersect = |s: [f64; 2], e: [f64; 2]| -> [f64; 2] {
            let dx = e[0] - s[0]; // seg dir
            let dy = e[1] - s[1];
            // 2-D cross: seg_dir × clip_dir
            let denom = dx * ey - dy * ex;
            if denom.abs() < 1e-15 {
                return s;
            }
            // (clip_start - S) × clip_dir
            let t = ((edge_start[0] - s[0]) * ey - (edge_start[1] - s[1]) * ex) / denom;
            [s[0] + t * dx, s[1] + t * dy]
        };

        let k = input.len();
        for j in 0..k {
            let current = input[j];
            let previous = input[(j + k - 1) % k];

            if inside(current) {
                if !inside(previous) {
                    output.push(intersect(previous, current));
                }
                output.push(current);
            } else if inside(previous) {
                output.push(intersect(previous, current));
            }
        }
    }

    output
}

// ─── Half-space clip (3-D) ────────────────────────────────────────────────────

/// Clip a 3-D polygon (`poly`) against a single half-space defined by
/// `plane_normal · x ≥ plane_d`.
///
/// Vertices on the positive side (inside) are kept; vertices on the negative
/// side are removed, with new vertices inserted at the intersection with the
/// plane boundary.
pub fn clip_polygon_by_halfspace(
    poly: &[[f64; 3]],
    plane_normal: [f64; 3],
    plane_d: f64,
) -> Vec<[f64; 3]> {
    if poly.is_empty() {
        return Vec::new();
    }

    let mut output = Vec::with_capacity(poly.len() + 1);
    let n = poly.len();

    for i in 0..n {
        let current = poly[i];
        let previous = poly[(i + n - 1) % n];

        let d_cur = dot3(plane_normal, current) - plane_d;
        let d_prev = dot3(plane_normal, previous) - plane_d;

        if d_cur >= 0.0 {
            // current is inside
            if d_prev < 0.0 {
                // previous was outside → insert intersection
                let t = d_prev / (d_prev - d_cur);
                output.push(add3(previous, scale3(sub3(current, previous), t)));
            }
            output.push(current);
        } else if d_prev >= 0.0 {
            // current is outside, previous was inside → insert intersection
            let t = d_prev / (d_prev - d_cur);
            output.push(add3(previous, scale3(sub3(current, previous), t)));
        }
    }

    output
}

// ─── Plane projection helpers ────────────────────────────────────────────────

/// Project a 3-D point `pt` onto the 2-D coordinate system spanned by `u`
/// and `v` at `plane_origin`.
pub fn project_point_to_plane(
    pt: [f64; 3],
    plane_origin: [f64; 3],
    u: [f64; 3],
    v: [f64; 3],
) -> [f64; 2] {
    let d = sub3(pt, plane_origin);
    [dot3(d, u), dot3(d, v)]
}

/// Reconstruct a 3-D point from its 2-D projection `pt2d` on the plane
/// defined by `plane_origin`, `u`, and `v`.
pub fn unproject_point_from_plane(
    pt2d: [f64; 2],
    plane_origin: [f64; 3],
    u: [f64; 3],
    v: [f64; 3],
) -> [f64; 3] {
    add3(plane_origin, add3(scale3(u, pt2d[0]), scale3(v, pt2d[1])))
}

// ─── Box manifold construction ────────────────────────────────────────────────

/// Extract the 3×3 rotation sub-matrix from a column-major 4×4 homogeneous
/// transform as three column vectors (the local X, Y, Z axes in world space).
fn rotation_axes(t: [[f64; 4]; 4]) -> [[f64; 3]; 3] {
    [
        [t[0][0], t[1][0], t[2][0]],
        [t[0][1], t[1][1], t[2][1]],
        [t[0][2], t[1][2], t[2][2]],
    ]
}

/// Extract the translation vector from a 4×4 homogeneous transform.
fn translation(t: [[f64; 4]; 4]) -> [f64; 3] {
    [t[0][3], t[1][3], t[2][3]]
}

/// Transform a local-space point by a 4×4 homogeneous matrix.
fn transform_point(m: [[f64; 4]; 4], p: [f64; 3]) -> [f64; 3] {
    [
        m[0][0] * p[0] + m[0][1] * p[1] + m[0][2] * p[2] + m[0][3],
        m[1][0] * p[0] + m[1][1] * p[1] + m[1][2] * p[2] + m[1][3],
        m[2][0] * p[0] + m[2][1] * p[1] + m[2][2] * p[2] + m[2][3],
    ]
}

/// Return the face of a box (given its world axes and half-extents) whose
/// outward normal is most aligned with `dir`.  Returns the face normal and
/// the four world-space corner vertices.
fn best_face(
    center: [f64; 3],
    axes: [[f64; 3]; 3],
    half: [f64; 3],
    dir: [f64; 3],
) -> ([f64; 3], Vec<[f64; 3]>) {
    // Find the axis most aligned with `dir`
    let mut best_dot = -f64::INFINITY;
    let mut best_axis = 0usize;
    let mut best_sign = 1.0_f64;
    for i in 0..3 {
        let d = dot3(axes[i], dir);
        if d.abs() > best_dot {
            best_dot = d.abs();
            best_axis = i;
            best_sign = d.signum();
        }
    }

    let face_normal = scale3(axes[best_axis], best_sign);

    // The two tangent axes
    let u_axis = best_axis;
    let v_axis = (best_axis + 1) % 3;
    let w_axis = (best_axis + 2) % 3;

    // Face center
    let fc = add3(center, scale3(face_normal, half[best_axis]));

    // Four corners
    let hu = scale3(axes[u_axis], half[u_axis]);
    let hv = scale3(axes[w_axis], half[w_axis]);
    let _ = v_axis; // suppress unused warning
    let corners = vec![
        add3(add3(fc, hu), hv),
        add3(sub3(fc, hu), hv),
        sub3(sub3(fc, hu), hv),
        add3(sub3(fc, hv), hu),
    ];

    (face_normal, corners)
}

/// Build a contact manifold for a box–box collision by:
///
/// 1. Identifying the reference face (on box A, most facing the contact normal)
///    and the incident face (on box B, most facing −normal).
/// 2. Clipping the incident face polygon against the four side planes of the
///    reference face.
/// 3. Projecting clipped points onto the reference plane and keeping those at
///    or below it (i.e., penetrating).
///
/// The transform matrices are **row-major** 4×4 homogeneous matrices where
/// `transform[row][col]`.  Translation is stored in the last column
/// (`transform[row][3]`).
#[allow(clippy::too_many_arguments)]
pub fn build_box_manifold(
    half_extents_a: [f64; 3],
    transform_a: [[f64; 4]; 4],
    half_extents_b: [f64; 3],
    transform_b: [[f64; 4]; 4],
    contact_normal: [f64; 3],
    penetration_depth: f64,
) -> ContactManifold {
    let center_a = translation(transform_a);
    let center_b = translation(transform_b);
    let axes_a = rotation_axes(transform_a);
    let axes_b = rotation_axes(transform_b);

    // Reference face: box A face whose normal is most aligned with contact_normal
    let (ref_normal, _ref_corners) = best_face(center_a, axes_a, half_extents_a, contact_normal);

    // Incident face: box B face whose normal is most aligned with -contact_normal
    let neg_normal = scale3(contact_normal, -1.0);
    let (_inc_normal, inc_corners) = best_face(center_b, axes_b, half_extents_b, neg_normal);

    // Build 2-D coordinate system on the reference plane
    // Choose a tangent vector perpendicular to ref_normal
    let tangent_u = {
        let candidate = if ref_normal[0].abs() < 0.9 {
            [1.0_f64, 0.0, 0.0]
        } else {
            [0.0_f64, 1.0, 0.0]
        };
        normalize3(cross3(candidate, ref_normal))
    };
    let tangent_v = cross3(ref_normal, tangent_u);

    // Project incident face corners onto 2-D
    let inc_2d: Vec<[f64; 2]> = inc_corners
        .iter()
        .map(|&p| project_point_to_plane(p, center_a, tangent_u, tangent_v))
        .collect();

    // Build 2-D reference face clipping polygon (a quad in 2-D)
    // We compute the half-extents of A projected onto the reference plane.
    // For simplicity we use the reference face corners projected into 2-D.
    let ref_face_normal_world = ref_normal;
    // Recompute reference face corners for clipping
    let best_axis_idx = {
        let mut best = 0usize;
        let mut best_d = -f64::INFINITY;
        for i in 0..3 {
            let d = dot3(axes_a[i], ref_face_normal_world).abs();
            if d > best_d {
                best_d = d;
                best = i;
            }
        }
        best
    };
    let u_idx = (best_axis_idx + 1) % 3;
    let w_idx = (best_axis_idx + 2) % 3;

    // Reference face in 2-D (clipping polygon)
    let hu = half_extents_a[u_idx];
    let hv = half_extents_a[w_idx];
    let ref_clip_2d: Vec<[f64; 2]> = vec![[hu, hv], [-hu, hv], [-hu, -hv], [hu, -hv]];

    // Sutherland-Hodgman clip in 2-D
    let clipped_2d = sutherland_hodgman_2d(&inc_2d, &ref_clip_2d);

    // Unproject back to 3-D and keep points below (or on) the reference plane
    let ref_plane_d = dot3(ref_normal, center_a) + half_extents_a[best_axis_idx];

    let mut manifold = ContactManifold::new(contact_normal);
    for pt2d in &clipped_2d {
        let pt3d = unproject_point_from_plane(*pt2d, center_a, tangent_u, tangent_v);
        let signed_dist = dot3(ref_normal, pt3d) - ref_plane_d;
        if signed_dist <= 1e-4 {
            // depth is how far below the reference plane the point is
            let depth = (-signed_dist).max(0.0).max(penetration_depth * 0.01);
            manifold.contacts.push(ContactPoint {
                point: pt3d,
                normal: contact_normal,
                depth,
            });
        }
    }

    manifold
}

// ─── Manifold reduction ───────────────────────────────────────────────────────

/// Reduce a manifold to at most `max_contacts` points using a farthest-point
/// heuristic:
///
/// 1. Keep the contact with the greatest penetration depth.
/// 2. Iteratively add the contact that is farthest (in 3-D) from the already
///    selected set.
pub fn reduce_manifold(manifold: &mut ContactManifold, max_contacts: usize) {
    if manifold.contacts.len() <= max_contacts {
        return;
    }

    let contacts = &manifold.contacts;
    let n = contacts.len();
    let mut selected: Vec<usize> = Vec::with_capacity(max_contacts);

    // Seed with the deepest point
    let deepest = contacts
        .iter()
        .enumerate()
        .max_by(|(_, a), (_, b)| {
            a.depth
                .partial_cmp(&b.depth)
                .unwrap_or(std::cmp::Ordering::Equal)
        })
        .map(|(i, _)| i)
        .unwrap_or(0);
    selected.push(deepest);

    while selected.len() < max_contacts {
        let mut best_idx = 0usize;
        let mut best_dist = -1.0_f64;

        for i in 0..n {
            if selected.contains(&i) {
                continue;
            }
            // Minimum distance from this candidate to the already-selected set
            let min_d = selected
                .iter()
                .map(|&j| {
                    let d = sub3(contacts[i].point, contacts[j].point);
                    dot3(d, d)
                })
                .fold(f64::INFINITY, f64::min);

            if min_d > best_dist {
                best_dist = min_d;
                best_idx = i;
            }
        }

        selected.push(best_idx);
    }

    // Retain only the selected contacts
    let kept: Vec<ContactPoint> = selected.iter().map(|&i| contacts[i].clone()).collect();
    manifold.contacts = kept;
}

// ─── Warm-start impulse data ──────────────────────────────────────────────────

/// Warm-start impulse data stored per contact point across simulation frames.
///
/// These values are transferred from the previous frame's solver to seed the
/// current frame's iterative solver (warm-starting), which dramatically
/// improves convergence speed and stability.
#[derive(Debug, Clone, Default)]
pub struct WarmStartData {
    /// Accumulated normal impulse from the previous frame (N·s).
    pub normal_impulse: f64,
    /// Accumulated tangential impulse (friction) in the contact plane (N·s).
    /// Stores impulses along both tangent axes.
    pub tangent_impulse: [f64; 2],
    /// Number of consecutive simulation frames this contact has been active.
    pub age: u32,
}

impl WarmStartData {
    /// Create new zero-initialised warm-start data.
    pub fn new() -> Self {
        Self::default()
    }

    /// Scale all impulses by `factor` (used for time-step ratio correction).
    pub fn scale(&mut self, factor: f64) {
        self.normal_impulse *= factor;
        self.tangent_impulse[0] *= factor;
        self.tangent_impulse[1] *= factor;
    }

    /// Reset impulses to zero (contact broken or new contact).
    pub fn reset(&mut self) {
        self.normal_impulse = 0.0;
        self.tangent_impulse = [0.0, 0.0];
        self.age = 0;
    }

    /// Increment age counter (called once per frame when contact persists).
    pub fn tick(&mut self) {
        self.age = self.age.saturating_add(1);
    }

    /// Whether this contact is considered "warm" (has valid cached impulses).
    pub fn is_warm(&self) -> bool {
        self.age > 0
    }
}

// ─── WarmContactPoint ─────────────────────────────────────────────────────────

/// A contact point with embedded warm-start data for solver efficiency.
#[derive(Debug, Clone)]
pub struct WarmContactPoint {
    /// Geometric contact data.
    pub contact: ContactPoint,
    /// Warm-start impulse data from the previous frame.
    pub warm: WarmStartData,
}

impl WarmContactPoint {
    /// Create a new `WarmContactPoint` with zero warm-start data.
    pub fn new(contact: ContactPoint) -> Self {
        Self {
            contact,
            warm: WarmStartData::new(),
        }
    }

    /// Create a new point, copying warm-start data from a previous contact.
    pub fn with_warm(contact: ContactPoint, warm: WarmStartData) -> Self {
        Self { contact, warm }
    }
}

// ─── 4-point manifold builder ────────────────────────────────────────────────

/// Reduce a set of contact points to exactly four representatives that
/// maximise the covered area in the contact plane.
///
/// This is the standard algorithm used in production rigid-body engines:
/// 1. Keep the deepest contact point.
/// 2. Keep the point farthest from point 1.
/// 3. Keep the point that maximises the triangle area (p1, p2, candidate).
/// 4. Keep the point that maximises the quadrilateral area (p1…p3, candidate).
///
/// If fewer than 4 points are provided they are returned unchanged.
pub fn build_4point_manifold(points: &[ContactPoint]) -> Vec<ContactPoint> {
    if points.len() <= 4 {
        return points.to_vec();
    }

    // Step 1: deepest penetration.
    let i0 = points
        .iter()
        .enumerate()
        .max_by(|(_, a), (_, b)| {
            a.depth
                .partial_cmp(&b.depth)
                .unwrap_or(std::cmp::Ordering::Equal)
        })
        .map(|(i, _)| i)
        .unwrap_or(0);

    // Step 2: farthest from p0.
    let p0 = points[i0].point;
    let i1 = points
        .iter()
        .enumerate()
        .filter(|(i, _)| *i != i0)
        .max_by(|(_, a), (_, b)| {
            dist_sq3(a.point, p0)
                .partial_cmp(&dist_sq3(b.point, p0))
                .unwrap_or(std::cmp::Ordering::Equal)
        })
        .map(|(i, _)| i)
        .unwrap_or(if i0 == 0 { 1 } else { 0 });

    // Step 3: maximise triangle area with p0, p1.
    let p1 = points[i1].point;
    let i2 = points
        .iter()
        .enumerate()
        .filter(|(i, _)| *i != i0 && *i != i1)
        .max_by(|(_, a), (_, b)| {
            triangle_area_sq(p0, p1, a.point)
                .partial_cmp(&triangle_area_sq(p0, p1, b.point))
                .unwrap_or(std::cmp::Ordering::Equal)
        })
        .map(|(i, _)| i);

    let i2 = match i2 {
        Some(idx) => idx,
        None => return vec![points[i0].clone(), points[i1].clone()],
    };

    // Step 4: maximise quadrilateral area.
    let p2 = points[i2].point;
    let i3 = points
        .iter()
        .enumerate()
        .filter(|(i, _)| *i != i0 && *i != i1 && *i != i2)
        .max_by(|(_, a), (_, b)| {
            quad_area_sq(p0, p1, p2, a.point)
                .partial_cmp(&quad_area_sq(p0, p1, p2, b.point))
                .unwrap_or(std::cmp::Ordering::Equal)
        })
        .map(|(i, _)| i);

    let mut result = vec![points[i0].clone(), points[i1].clone(), points[i2].clone()];
    if let Some(idx) = i3 {
        result.push(points[idx].clone());
    }
    result
}

// ─── Frame-to-frame contact persistence ──────────────────────────────────────

/// Threshold (squared) for matching contact points across frames.
const MATCH_DIST_SQ_MANIFOLD: f64 = 5e-4; // (2.2 cm)²

/// Transfer warm-start data from a previous-frame manifold to a newly computed
/// manifold for the same body pair.
///
/// For each new contact point the algorithm searches the old manifold for a
/// spatially close match (within `MATCH_DIST_SQ_MANIFOLD`). When a match is
/// found the old [`WarmStartData`] is transferred; the age counter is
/// incremented. Unmatched new points receive zero-initialised warm data.
pub fn transfer_warm_start(
    new_contacts: &[ContactPoint],
    old_warm: &[WarmContactPoint],
) -> Vec<WarmContactPoint> {
    new_contacts
        .iter()
        .map(|new_pt| {
            // Find the closest point in the old manifold.
            let best = old_warm.iter().min_by(|a, b| {
                let da = dist_sq3(a.contact.point, new_pt.point);
                let db = dist_sq3(b.contact.point, new_pt.point);
                da.partial_cmp(&db).unwrap_or(std::cmp::Ordering::Equal)
            });

            if let Some(old_pt) = best {
                let d_sq = dist_sq3(old_pt.contact.point, new_pt.point);
                if d_sq < MATCH_DIST_SQ_MANIFOLD {
                    let mut warm = old_pt.warm.clone();
                    warm.tick();
                    return WarmContactPoint::with_warm(new_pt.clone(), warm);
                }
            }

            WarmContactPoint::new(new_pt.clone())
        })
        .collect()
}

// ─── WarmManifold ─────────────────────────────────────────────────────────────

/// A contact manifold with per-point warm-start impulse data.
///
/// Wraps a [`ContactManifold`] and adds the frame-to-frame persistence state
/// needed to warm-start the iterative constraint solver.
#[derive(Debug, Clone)]
pub struct WarmManifold {
    /// Contact normal for this manifold.
    pub normal: [f64; 3],
    /// Contact points with warm-start data.
    pub points: Vec<WarmContactPoint>,
}

impl WarmManifold {
    /// Create an empty `WarmManifold`.
    pub fn new(normal: [f64; 3]) -> Self {
        Self {
            normal,
            points: Vec::new(),
        }
    }

    /// Build a `WarmManifold` from a [`ContactManifold`], transferring
    /// warm-start data from the previous frame (`prev`).
    ///
    /// Reduces the contact set to at most 4 representative points first.
    pub fn from_manifold(manifold: &ContactManifold, prev: Option<&WarmManifold>) -> Self {
        let reduced = build_4point_manifold(&manifold.contacts);

        let points = if let Some(old) = prev {
            transfer_warm_start(&reduced, &old.points)
        } else {
            reduced
                .iter()
                .map(|c| WarmContactPoint::new(c.clone()))
                .collect()
        };

        Self {
            normal: manifold.normal,
            points,
        }
    }

    /// Number of contact points.
    pub fn len(&self) -> usize {
        self.points.len()
    }

    /// Whether this manifold has no contact points.
    pub fn is_empty(&self) -> bool {
        self.points.is_empty()
    }

    /// Iterate over the warm contact points.
    pub fn iter(&self) -> std::slice::Iter<'_, WarmContactPoint> {
        self.points.iter()
    }

    /// Mutable access to warm contact points.
    pub fn iter_mut(&mut self) -> std::slice::IterMut<'_, WarmContactPoint> {
        self.points.iter_mut()
    }

    /// Apply an impulse scale factor to all warm-start data.
    ///
    /// Used when the simulation timestep changes between frames.
    pub fn scale_impulses(&mut self, factor: f64) {
        for pt in &mut self.points {
            pt.warm.scale(factor);
        }
    }

    /// Reset all warm-start data (e.g., after a large time gap or teleport).
    pub fn reset_warm_start(&mut self) {
        for pt in &mut self.points {
            pt.warm.reset();
        }
    }

    /// Average contact normal recomputed from all points.
    pub fn recompute_normal(&mut self) {
        if self.points.is_empty() {
            return;
        }
        let mut avg = [0.0f64; 3];
        for pt in &self.points {
            avg = add3(avg, pt.contact.normal);
        }
        let n = self.points.len() as f64;
        self.normal = normalize3(scale3(avg, 1.0 / n));
    }
}

// ─── Geometry helpers for 4-point manifold ───────────────────────────────────

#[inline]
fn dist_sq3(a: [f64; 3], b: [f64; 3]) -> f64 {
    let d = sub3(a, b);
    dot3(d, d)
}

#[inline]
fn triangle_area_sq(a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> f64 {
    let e1 = sub3(b, a);
    let e2 = sub3(c, a);
    let cx = cross3(e1, e2);
    dot3(cx, cx)
}

#[inline]
fn quad_area_sq(a: [f64; 3], b: [f64; 3], c: [f64; 3], d: [f64; 3]) -> f64 {
    triangle_area_sq(a, b, c) + triangle_area_sq(a, c, d)
}

// ─── Persistent patch tracking ───────────────────────────────────────────────

/// A contact patch that persists across frames, tracking history of contacts.
///
/// Useful for integrating contact forces and detecting sustained contact.
#[derive(Debug, Clone)]
pub struct PersistentPatch {
    /// Shared contact normal.
    pub normal: [f64; 3],
    /// Warm contact points (max 4).
    pub points: Vec<WarmContactPoint>,
    /// Number of consecutive frames this patch has been active.
    pub frame_count: u32,
    /// Whether this patch was active in the most recent frame.
    pub active: bool,
}

impl PersistentPatch {
    /// Create a new empty persistent patch.
    pub fn new(normal: [f64; 3]) -> Self {
        Self {
            normal,
            points: Vec::new(),
            frame_count: 0,
            active: true,
        }
    }

    /// Update the patch from a new [`ContactManifold`].
    ///
    /// Transfers warm-start data from matching previous points, increments
    /// `frame_count`, and reduces to at most 4 points.
    pub fn update(&mut self, manifold: &ContactManifold) {
        let reduced = build_4point_manifold(&manifold.contacts);
        let new_warm = transfer_warm_start(&reduced, &self.points);
        self.points = new_warm;
        self.normal = manifold.normal;
        self.frame_count = self.frame_count.saturating_add(1);
        self.active = true;
    }

    /// Mark the patch as inactive (no collision detected this frame).
    pub fn deactivate(&mut self) {
        self.active = false;
    }

    /// Whether the patch has been active for at least `n` consecutive frames.
    pub fn is_sustained(&self, n: u32) -> bool {
        self.active && self.frame_count >= n
    }

    /// Maximum penetration depth across all contact points.
    pub fn max_depth(&self) -> f64 {
        self.points
            .iter()
            .map(|p| p.contact.depth)
            .fold(0.0f64, f64::max)
    }

    /// Average contact position across all points.
    pub fn average_position(&self) -> Option<[f64; 3]> {
        if self.points.is_empty() {
            return None;
        }
        let n = self.points.len() as f64;
        let mut avg = [0.0f64; 3];
        for p in &self.points {
            avg = add3(avg, p.contact.point);
        }
        Some(scale3(avg, 1.0 / n))
    }
}

// ─── Contact normal smoothing ────────────────────────────────────────────────

/// Smoothly blend a new contact normal with a previous one.
///
/// Interpolates between `prev_normal` and `new_normal` by `alpha` ∈ \[0,1\]
/// and re-normalises. `alpha = 1.0` returns `new_normal` unchanged.
pub fn smooth_contact_normal(prev_normal: [f64; 3], new_normal: [f64; 3], alpha: f64) -> [f64; 3] {
    let blended = add3(scale3(prev_normal, 1.0 - alpha), scale3(new_normal, alpha));
    normalize3(blended)
}

/// Apply exponential moving average smoothing to contact normals in a manifold.
///
/// Iterates over all contact points and blends their normals with the manifold
/// normal using factor `alpha`.
pub fn smooth_manifold_normals(manifold: &mut ContactManifold, alpha: f64) {
    let base = manifold.normal;
    for pt in &mut manifold.contacts {
        pt.normal = smooth_contact_normal(base, pt.normal, alpha);
    }
}

// ─── Speculative contacts ─────────────────────────────────────────────────────

/// A speculative contact used for predictive collision response.
///
/// Generated before penetration occurs, based on the distance between shapes
/// and their relative velocities.
#[derive(Debug, Clone)]
pub struct SpeculativeContactPoint {
    /// World-space position of the speculative contact.
    pub point: [f64; 3],
    /// Contact normal (from B toward A).
    pub normal: [f64; 3],
    /// Signed separation distance (negative = overlapping, positive = gap).
    pub separation: f64,
    /// Relative velocity along the normal (closing speed).
    pub closing_speed: f64,
}

impl SpeculativeContactPoint {
    /// Create a speculative contact point.
    pub fn new(point: [f64; 3], normal: [f64; 3], separation: f64, closing_speed: f64) -> Self {
        Self {
            point,
            normal,
            separation,
            closing_speed,
        }
    }

    /// Whether the speculative contact will produce a collision within `dt`.
    pub fn will_collide(&self, dt: f64) -> bool {
        self.separation - self.closing_speed * dt < 0.0
    }

    /// Predicted penetration depth at end of `dt`.
    pub fn predicted_depth(&self, dt: f64) -> f64 {
        -(self.separation - self.closing_speed * dt)
    }
}

/// Generate speculative contacts for a contact manifold given relative velocity.
///
/// Returns contacts where the closing speed times `dt` exceeds the separation.
pub fn generate_speculative_contacts(
    manifold: &ContactManifold,
    normal: [f64; 3],
    closing_speed: f64,
    separation: f64,
    dt: f64,
) -> Vec<SpeculativeContactPoint> {
    manifold
        .contacts
        .iter()
        .filter_map(|pt| {
            let spec = SpeculativeContactPoint::new(pt.point, normal, separation, closing_speed);
            if spec.will_collide(dt) {
                Some(spec)
            } else {
                None
            }
        })
        .collect()
}

// ─── Manifold frame rotation ──────────────────────────────────────────────────

/// Rotate all contact points in a manifold by a 3×3 rotation matrix.
///
/// The rotation matrix is given in row-major format `[[row0], [row1], [row2]]`.
pub fn rotate_manifold_points(manifold: &mut ContactManifold, rot: [[f64; 3]; 3]) {
    for pt in &mut manifold.contacts {
        pt.point = mat3_mul_vec3(rot, pt.point);
        pt.normal = mat3_mul_vec3(rot, pt.normal);
    }
    manifold.normal = mat3_mul_vec3(rot, manifold.normal);
}

/// Multiply a 3×3 rotation matrix by a 3-vector.
#[inline]
fn mat3_mul_vec3(m: [[f64; 3]; 3], v: [f64; 3]) -> [f64; 3] {
    [
        m[0][0] * v[0] + m[0][1] * v[1] + m[0][2] * v[2],
        m[1][0] * v[0] + m[1][1] * v[1] + m[1][2] * v[2],
        m[2][0] * v[0] + m[2][1] * v[1] + m[2][2] * v[2],
    ]
}

// ─── Friction basis vectors ───────────────────────────────────────────────────

/// Compute an orthonormal friction basis `(t1, t2)` from a contact normal.
///
/// Returns two vectors tangent to the contact plane such that
/// `{normal, t1, t2}` form a right-handed frame.
pub fn friction_basis(normal: [f64; 3]) -> ([f64; 3], [f64; 3]) {
    // Pick a vector not parallel to normal
    let candidate = if normal[0].abs() < 0.9 {
        [1.0f64, 0.0, 0.0]
    } else {
        [0.0f64, 1.0, 0.0]
    };
    let t1 = normalize3(cross3(normal, candidate));
    let t2 = cross3(normal, t1);
    (t1, t2)
}

/// Decompose a velocity vector into normal and tangential components.
///
/// Returns `(v_normal, v_tangential)` where `v_normal` is along the contact
/// normal and `v_tangential` is in the contact plane.
pub fn decompose_velocity(vel: [f64; 3], normal: [f64; 3]) -> ([f64; 3], [f64; 3]) {
    let vn = dot3(vel, normal);
    let v_normal = scale3(normal, vn);
    let v_tangential = sub3(vel, v_normal);
    (v_normal, v_tangential)
}

/// Compute the friction impulse magnitude using Coulomb friction.
///
/// Given the normal impulse magnitude `lambda_n` and friction coefficient `mu`,
/// returns the maximum tangential impulse magnitude.
pub fn coulomb_friction_limit(lambda_n: f64, mu: f64) -> f64 {
    (lambda_n * mu).max(0.0)
}

// ─── Manifold quality metrics ─────────────────────────────────────────────────

/// Compute the quality score of a manifold based on point spread and depths.
///
/// Returns a value in `[0, 1]` where 1.0 means well-distributed, deep contacts.
/// Uses the ratio of the actual spread to the maximum possible spread.
pub fn manifold_quality(manifold: &ContactManifold) -> f64 {
    let n = manifold.contacts.len();
    if n == 0 {
        return 0.0;
    }

    // Average depth score
    let avg_depth = manifold.contacts.iter().map(|p| p.depth).sum::<f64>() / n as f64;
    let depth_score = (avg_depth / (avg_depth + 1.0)).min(1.0);

    if n < 2 {
        return depth_score * 0.5;
    }

    // Spread score: max pairwise distance among contact points
    let mut max_dist_sq = 0.0f64;
    for i in 0..n {
        for j in (i + 1)..n {
            let d = sub3(manifold.contacts[i].point, manifold.contacts[j].point);
            let dsq = dot3(d, d);
            if dsq > max_dist_sq {
                max_dist_sq = dsq;
            }
        }
    }
    let spread_score = (max_dist_sq.sqrt() / (max_dist_sq.sqrt() + 1.0)).min(1.0);

    (depth_score + spread_score) * 0.5
}

// ─── Tests ────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // ── Sutherland-Hodgman 2D ──────────────────────────────────────────────────

    #[test]
    fn test_sh_clip_square_against_x_ge_0() {
        // Subject: unit square centred at origin
        let subject = [[-1.0, -1.0], [1.0, -1.0], [1.0, 1.0], [-1.0, 1.0]];
        // Clip polygon for x ≥ 0: a large right half-plane represented as a
        // convex polygon (a tall vertical strip)
        let clip = [[0.0, -2.0], [2.0, -2.0], [2.0, 2.0], [0.0, 2.0]];

        let result = sutherland_hodgman_2d(&subject, &clip);

        assert!(!result.is_empty(), "clipped polygon must not be empty");
        // All result vertices must have x ≥ 0
        for v in &result {
            assert!(
                v[0] >= -1e-10,
                "vertex {:?} has negative x after clipping against x≥0",
                v
            );
        }
        // The result should have 4 vertices (the right half of the square)
        assert_eq!(
            result.len(),
            4,
            "right half of square should have 4 vertices, got {:?}",
            result
        );
        // Area check: the right half of the unit square has area 2.0
        let area = polygon_area_2d(&result);
        assert!(
            (area - 2.0).abs() < 1e-10,
            "expected area 2.0, got {}",
            area
        );
    }

    #[test]
    fn test_sh_no_clipping_needed() {
        // Subject fully inside the clip polygon
        let subject = [[0.1, 0.1], [0.9, 0.1], [0.9, 0.9], [0.1, 0.9]];
        let clip = [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
        let result = sutherland_hodgman_2d(&subject, &clip);
        assert_eq!(result.len(), 4);
    }

    #[test]
    fn test_sh_fully_outside() {
        let subject = [[-3.0, -3.0], [-2.0, -3.0], [-2.0, -2.0], [-3.0, -2.0]];
        let clip = [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
        let result = sutherland_hodgman_2d(&subject, &clip);
        assert!(
            result.is_empty(),
            "fully outside polygon should be clipped to empty"
        );
    }

    // ── clip_polygon_by_halfspace ──────────────────────────────────────────────

    #[test]
    fn test_clip_polygon_above_plane() {
        // A quad entirely above z = 0 plane (normal = [0,0,1], d = 0)
        let quad = [
            [1.0, 1.0, 1.0],
            [-1.0, 1.0, 1.0],
            [-1.0, -1.0, 1.0],
            [1.0, -1.0, 1.0],
        ];
        let result = clip_polygon_by_halfspace(&quad, [0.0, 0.0, 1.0], 0.0);
        // All 4 vertices are above (inside) the plane so the polygon is unchanged.
        assert_eq!(
            result.len(),
            4,
            "quad fully above plane must be unchanged, got {:?}",
            result
        );
    }

    #[test]
    fn test_clip_polygon_straddles_plane() {
        // A quad that straddles z = 0: two vertices above, two below.
        let quad = [
            [1.0, 0.0, 1.0],
            [-1.0, 0.0, 1.0],
            [-1.0, 0.0, -1.0],
            [1.0, 0.0, -1.0],
        ];
        let result = clip_polygon_by_halfspace(&quad, [0.0, 0.0, 1.0], 0.0);
        // Should produce a 4-vertex polygon (the top half)
        assert!(
            result.len() >= 3,
            "straddling quad must produce at least a triangle, got {:?}",
            result
        );
        for v in &result {
            assert!(
                v[2] >= -1e-10,
                "all clipped vertices must have z≥0, got {:?}",
                v
            );
        }
    }

    #[test]
    fn test_clip_polygon_below_plane() {
        // All vertices below z=0 → empty result
        let quad = [
            [1.0, 0.0, -1.0],
            [-1.0, 0.0, -1.0],
            [-1.0, 0.0, -2.0],
            [1.0, 0.0, -2.0],
        ];
        let result = clip_polygon_by_halfspace(&quad, [0.0, 0.0, 1.0], 0.0);
        assert!(
            result.is_empty(),
            "quad fully below plane must be clipped to empty"
        );
    }

    // ── reduce_manifold ───────────────────────────────────────────────────────

    #[test]
    fn test_reduce_manifold_8_to_4() {
        // Create 8 contacts spread around a unit circle in the XY plane.
        let mut manifold = ContactManifold::new([0.0, 0.0, 1.0]);
        for i in 0..8 {
            let angle = (i as f64) * std::f64::consts::TAU / 8.0;
            manifold.contacts.push(ContactPoint {
                point: [angle.cos(), angle.sin(), 0.0],
                normal: [0.0, 0.0, 1.0],
                depth: 0.1 + (i as f64) * 0.01, // slightly varying depths
            });
        }

        reduce_manifold(&mut manifold, 4);

        assert_eq!(
            manifold.contacts.len(),
            4,
            "manifold should be reduced to 4 contacts"
        );
    }

    #[test]
    fn test_reduce_manifold_already_small() {
        let mut manifold = ContactManifold::new([0.0, 1.0, 0.0]);
        manifold.contacts.push(ContactPoint {
            point: [0.0, 0.0, 0.0],
            normal: [0.0, 1.0, 0.0],
            depth: 0.05,
        });
        reduce_manifold(&mut manifold, 4);
        assert_eq!(manifold.contacts.len(), 1);
    }

    // ── project/unproject round-trip ─────────────────────────────────────────

    #[test]
    fn test_project_unproject_roundtrip() {
        let origin = [1.0, 2.0, 3.0];
        let u = normalize3([1.0, 0.0, 0.0]);
        let v = normalize3([0.0, 1.0, 0.0]);
        let pt = [2.5, 3.5, 3.0]; // lies on the plane (z == origin.z)
        let pt2d = project_point_to_plane(pt, origin, u, v);
        let pt3d = unproject_point_from_plane(pt2d, origin, u, v);
        for i in 0..3 {
            assert!(
                (pt3d[i] - pt[i]).abs() < 1e-12,
                "round-trip mismatch at component {}: expected {}, got {}",
                i,
                pt[i],
                pt3d[i]
            );
        }
    }

    // ── Helper ────────────────────────────────────────────────────────────────

    /// Signed area of a 2-D polygon (positive for CCW winding).
    fn polygon_area_2d(poly: &[[f64; 2]]) -> f64 {
        let n = poly.len();
        let mut area = 0.0_f64;
        for i in 0..n {
            let j = (i + 1) % n;
            area += poly[i][0] * poly[j][1];
            area -= poly[j][0] * poly[i][1];
        }
        (area / 2.0).abs()
    }

    // ── WarmStartData ─────────────────────────────────────────────────────────

    #[test]
    fn test_warm_start_data_default_zero() {
        let w = WarmStartData::new();
        assert_eq!(w.normal_impulse, 0.0);
        assert_eq!(w.tangent_impulse, [0.0, 0.0]);
        assert_eq!(w.age, 0);
        assert!(!w.is_warm());
    }

    #[test]
    fn test_warm_start_data_tick_increments_age() {
        let mut w = WarmStartData::new();
        w.tick();
        assert_eq!(w.age, 1);
        assert!(w.is_warm());
        w.tick();
        assert_eq!(w.age, 2);
    }

    #[test]
    fn test_warm_start_data_scale() {
        let mut w = WarmStartData::new();
        w.normal_impulse = 10.0;
        w.tangent_impulse = [4.0, -2.0];
        w.scale(0.5);
        assert!((w.normal_impulse - 5.0).abs() < 1e-12);
        assert!((w.tangent_impulse[0] - 2.0).abs() < 1e-12);
        assert!((w.tangent_impulse[1] - (-1.0)).abs() < 1e-12);
    }

    #[test]
    fn test_warm_start_data_reset() {
        let mut w = WarmStartData::new();
        w.normal_impulse = 5.0;
        w.tangent_impulse = [1.0, 2.0];
        w.age = 10;
        w.reset();
        assert_eq!(w.normal_impulse, 0.0);
        assert_eq!(w.tangent_impulse, [0.0, 0.0]);
        assert_eq!(w.age, 0);
    }

    // ── build_4point_manifold ─────────────────────────────────────────────────

    fn make_cp(x: f64, z: f64, depth: f64) -> ContactPoint {
        ContactPoint {
            point: [x, 0.0, z],
            normal: [0.0, 1.0, 0.0],
            depth,
        }
    }

    #[test]
    fn test_build_4point_manifold_reduces_to_4() {
        // 8 points spread on a 2×2 grid, depths vary.
        let pts: Vec<ContactPoint> = vec![
            make_cp(-1.0, -1.0, 0.1),
            make_cp(1.0, -1.0, 0.2),
            make_cp(1.0, 1.0, 0.15),
            make_cp(-1.0, 1.0, 0.05),
            make_cp(0.0, 0.0, 0.25), // deepest — seed
            make_cp(0.5, -0.5, 0.08),
            make_cp(-0.5, 0.5, 0.12),
            make_cp(0.0, 1.0, 0.07),
        ];
        let result = build_4point_manifold(&pts);
        assert_eq!(result.len(), 4, "should reduce to exactly 4 points");
    }

    #[test]
    fn test_build_4point_manifold_fewer_than_4_unchanged() {
        let pts: Vec<ContactPoint> = vec![make_cp(0.0, 0.0, 0.1), make_cp(1.0, 0.0, 0.1)];
        let result = build_4point_manifold(&pts);
        assert_eq!(result.len(), 2, "fewer than 4 points returned unchanged");
    }

    #[test]
    fn test_build_4point_manifold_deepest_included() {
        let pts: Vec<ContactPoint> = vec![
            make_cp(-1.0, 0.0, 0.01),
            make_cp(0.0, 0.0, 0.99), // deepest
            make_cp(1.0, 0.0, 0.01),
            make_cp(0.0, 1.0, 0.01),
            make_cp(0.0, -1.0, 0.01),
        ];
        let result = build_4point_manifold(&pts);
        // The deepest point ([0,0,0], depth=0.99) must be in the result.
        let has_deepest = result.iter().any(|p| (p.depth - 0.99).abs() < 1e-9);
        assert!(
            has_deepest,
            "deepest point must be included in 4-point manifold"
        );
    }

    // ── transfer_warm_start ───────────────────────────────────────────────────

    #[test]
    fn test_transfer_warm_start_matching_point() {
        // Old manifold: one warm contact at origin.
        let old_pt = {
            let mut cp = WarmContactPoint::new(ContactPoint {
                point: [0.0, 0.0, 0.0],
                normal: [0.0, 1.0, 0.0],
                depth: 0.05,
            });
            cp.warm.normal_impulse = 42.0;
            cp.warm.tick();
            cp
        };
        let old_warm = vec![old_pt];

        // New contact very close to the old one.
        let new_contacts = vec![ContactPoint {
            point: [0.001, 0.0, 0.001], // within threshold
            normal: [0.0, 1.0, 0.0],
            depth: 0.04,
        }];

        let result = transfer_warm_start(&new_contacts, &old_warm);
        assert_eq!(result.len(), 1);
        assert!(
            result[0].warm.is_warm(),
            "warm-start should transfer to close contact"
        );
        assert!(
            (result[0].warm.normal_impulse - 42.0).abs() < 1e-9,
            "normal_impulse should be transferred"
        );
        assert_eq!(result[0].warm.age, 2, "age should be incremented");
    }

    #[test]
    fn test_transfer_warm_start_distant_point_gets_zero() {
        let old_pt = {
            let mut cp = WarmContactPoint::new(ContactPoint {
                point: [0.0, 0.0, 0.0],
                normal: [0.0, 1.0, 0.0],
                depth: 0.05,
            });
            cp.warm.normal_impulse = 10.0;
            cp.warm.tick();
            cp
        };
        let old_warm = vec![old_pt];

        // New contact far away.
        let new_contacts = vec![ContactPoint {
            point: [100.0, 0.0, 0.0],
            normal: [0.0, 1.0, 0.0],
            depth: 0.04,
        }];

        let result = transfer_warm_start(&new_contacts, &old_warm);
        assert_eq!(result.len(), 1);
        assert!(
            !result[0].warm.is_warm(),
            "distant contact should get no warm-start"
        );
        assert_eq!(result[0].warm.normal_impulse, 0.0);
    }

    // ── WarmManifold ──────────────────────────────────────────────────────────

    #[test]
    fn test_warm_manifold_from_manifold_no_prev() {
        let mut m = ContactManifold::new([0.0, 1.0, 0.0]);
        for i in 0..4 {
            m.contacts.push(make_cp(i as f64, 0.0, 0.1));
        }
        let wm = WarmManifold::from_manifold(&m, None);
        assert_eq!(wm.len(), 4);
        for pt in wm.iter() {
            assert!(!pt.warm.is_warm(), "first frame: no warm data");
        }
    }

    #[test]
    fn test_warm_manifold_scale_impulses() {
        let mut wm = WarmManifold::new([0.0, 1.0, 0.0]);
        let mut wcp = WarmContactPoint::new(make_cp(0.0, 0.0, 0.1));
        wcp.warm.normal_impulse = 8.0;
        wm.points.push(wcp);

        wm.scale_impulses(0.25);
        assert!((wm.points[0].warm.normal_impulse - 2.0).abs() < 1e-12);
    }

    #[test]
    fn test_warm_manifold_reset_warm_start() {
        let mut wm = WarmManifold::new([0.0, 1.0, 0.0]);
        let mut wcp = WarmContactPoint::new(make_cp(0.0, 0.0, 0.1));
        wcp.warm.normal_impulse = 5.0;
        wcp.warm.age = 3;
        wm.points.push(wcp);

        wm.reset_warm_start();
        assert_eq!(wm.points[0].warm.normal_impulse, 0.0);
        assert_eq!(wm.points[0].warm.age, 0);
    }

    #[test]
    fn test_warm_manifold_recompute_normal() {
        let mut wm = WarmManifold::new([1.0, 0.0, 0.0]);
        let normals = [[0.0, 1.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.8, 0.6]];
        for n in normals {
            let mut wcp = WarmContactPoint::new(ContactPoint {
                point: [0.0, 0.0, 0.0],
                normal: n,
                depth: 0.1,
            });
            wcp.warm = WarmStartData::new();
            wm.points.push(wcp);
        }
        wm.recompute_normal();
        let nl = (wm.normal[0] * wm.normal[0]
            + wm.normal[1] * wm.normal[1]
            + wm.normal[2] * wm.normal[2])
            .sqrt();
        assert!(
            (nl - 1.0).abs() < 1e-9,
            "recomputed normal should be unit length"
        );
    }

    // ── PersistentPatch ───────────────────────────────────────────────────────

    #[test]
    fn test_persistent_patch_new_empty() {
        let patch = PersistentPatch::new([0.0, 1.0, 0.0]);
        assert!(patch.points.is_empty());
        assert_eq!(patch.frame_count, 0);
        assert!(patch.active);
    }

    #[test]
    fn test_persistent_patch_update_increments_frame() {
        let mut patch = PersistentPatch::new([0.0, 1.0, 0.0]);
        let mut m = ContactManifold::new([0.0, 1.0, 0.0]);
        m.contacts.push(make_cp(0.0, 0.0, 0.1));
        patch.update(&m);
        assert_eq!(patch.frame_count, 1);
    }

    #[test]
    fn test_persistent_patch_update_twice() {
        let mut patch = PersistentPatch::new([0.0, 1.0, 0.0]);
        let mut m = ContactManifold::new([0.0, 1.0, 0.0]);
        m.contacts.push(make_cp(0.0, 0.0, 0.1));
        patch.update(&m);
        patch.update(&m);
        assert_eq!(patch.frame_count, 2);
    }

    #[test]
    fn test_persistent_patch_deactivate() {
        let mut patch = PersistentPatch::new([0.0, 1.0, 0.0]);
        patch.deactivate();
        assert!(!patch.active);
    }

    #[test]
    fn test_persistent_patch_is_sustained() {
        let mut patch = PersistentPatch::new([0.0, 1.0, 0.0]);
        let mut m = ContactManifold::new([0.0, 1.0, 0.0]);
        m.contacts.push(make_cp(0.0, 0.0, 0.1));
        patch.update(&m);
        patch.update(&m);
        assert!(patch.is_sustained(2), "should be sustained after 2 frames");
        assert!(
            !patch.is_sustained(3),
            "should not be sustained for 3 frames yet"
        );
    }

    #[test]
    fn test_persistent_patch_max_depth() {
        let mut patch = PersistentPatch::new([0.0, 1.0, 0.0]);
        let mut m = ContactManifold::new([0.0, 1.0, 0.0]);
        m.contacts.push(make_cp(0.0, 0.0, 0.05));
        m.contacts.push(make_cp(1.0, 0.0, 0.10));
        patch.update(&m);
        assert!(
            (patch.max_depth() - 0.10).abs() < 1e-10,
            "max_depth should be 0.10"
        );
    }

    #[test]
    fn test_persistent_patch_average_position() {
        let mut patch = PersistentPatch::new([0.0, 1.0, 0.0]);
        let mut m = ContactManifold::new([0.0, 1.0, 0.0]);
        m.contacts.push(make_cp(0.0, 0.0, 0.1));
        m.contacts.push(make_cp(2.0, 0.0, 0.1));
        patch.update(&m);
        let avg = patch.average_position().unwrap();
        assert!(
            (avg[0] - 1.0).abs() < 0.01,
            "average x should be ~1.0, got {}",
            avg[0]
        );
    }

    // ── smooth_contact_normal ─────────────────────────────────────────────────

    #[test]
    fn test_smooth_contact_normal_alpha_one() {
        let prev = [1.0f64, 0.0, 0.0];
        let new = [0.0, 1.0, 0.0];
        let result = smooth_contact_normal(prev, new, 1.0);
        assert!((result[0]).abs() < 1e-10, "alpha=1 should give new normal");
        assert!((result[1] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_smooth_contact_normal_alpha_zero() {
        let prev = [1.0f64, 0.0, 0.0];
        let new = [0.0, 1.0, 0.0];
        let result = smooth_contact_normal(prev, new, 0.0);
        assert!(
            (result[0] - 1.0).abs() < 1e-10,
            "alpha=0 should give prev normal"
        );
        assert!((result[1]).abs() < 1e-10);
    }

    #[test]
    fn test_smooth_contact_normal_midpoint_unit_length() {
        let prev = [1.0f64, 0.0, 0.0];
        let new = [0.0, 1.0, 0.0];
        let result = smooth_contact_normal(prev, new, 0.5);
        let l = (result[0] * result[0] + result[1] * result[1] + result[2] * result[2]).sqrt();
        assert!(
            (l - 1.0).abs() < 1e-10,
            "smoothed normal should be unit length, got {l}"
        );
    }

    #[test]
    fn test_smooth_manifold_normals() {
        let mut m = ContactManifold::new([0.0, 1.0, 0.0]);
        m.contacts.push(ContactPoint {
            point: [0.0, 0.0, 0.0],
            normal: [1.0, 0.0, 0.0],
            depth: 0.1,
        });
        smooth_manifold_normals(&mut m, 0.5);
        // The contact point normal should be blended
        let n = m.contacts[0].normal;
        let l = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
        assert!(
            (l - 1.0).abs() < 1e-10,
            "smoothed normal should be unit length"
        );
    }

    // ── SpeculativeContactPoint ───────────────────────────────────────────────

    #[test]
    fn test_speculative_contact_will_collide() {
        // separation=0.1, closing_speed=1.0, dt=0.2 → after dt: 0.1 - 0.2 = -0.1 < 0 → collide
        let spec = SpeculativeContactPoint::new([0.0; 3], [0.0, 1.0, 0.0], 0.1, 1.0);
        assert!(spec.will_collide(0.2));
    }

    #[test]
    fn test_speculative_contact_will_not_collide() {
        // separation=1.0, closing_speed=0.1, dt=0.1 → after dt: 1.0 - 0.01 > 0 → no collide
        let spec = SpeculativeContactPoint::new([0.0; 3], [0.0, 1.0, 0.0], 1.0, 0.1);
        assert!(!spec.will_collide(0.1));
    }

    #[test]
    fn test_speculative_contact_predicted_depth() {
        let spec = SpeculativeContactPoint::new([0.0; 3], [0.0, 1.0, 0.0], 0.1, 1.0);
        let depth = spec.predicted_depth(0.2);
        assert!((depth - 0.1).abs() < 1e-10, "Expected 0.1, got {depth}");
    }

    #[test]
    fn test_generate_speculative_contacts() {
        let mut m = ContactManifold::new([0.0, 1.0, 0.0]);
        m.contacts.push(ContactPoint {
            point: [0.0, 0.0, 0.0],
            normal: [0.0, 1.0, 0.0],
            depth: 0.05,
        });
        m.contacts.push(ContactPoint {
            point: [1.0, 0.0, 0.0],
            normal: [0.0, 1.0, 0.0],
            depth: 0.05,
        });
        let speculative = generate_speculative_contacts(&m, [0.0, 1.0, 0.0], 5.0, 0.1, 0.1);
        // 0.1 - 5.0*0.1 = -0.4 < 0 → both should collide
        assert_eq!(
            speculative.len(),
            2,
            "Both speculative contacts should be generated"
        );
    }

    #[test]
    fn test_generate_speculative_contacts_none_collide() {
        let mut m = ContactManifold::new([0.0, 1.0, 0.0]);
        m.contacts.push(ContactPoint {
            point: [0.0, 0.0, 0.0],
            normal: [0.0, 1.0, 0.0],
            depth: 0.05,
        });
        // separation huge, slow closing speed, small dt → won't collide
        let speculative = generate_speculative_contacts(&m, [0.0, 1.0, 0.0], 0.001, 100.0, 0.0001);
        assert_eq!(
            speculative.len(),
            0,
            "No speculative contacts should be generated"
        );
    }

    // ── rotate_manifold_points ────────────────────────────────────────────────

    #[test]
    fn test_rotate_manifold_points_identity() {
        let identity = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let mut m = ContactManifold::new([0.0, 1.0, 0.0]);
        m.contacts.push(ContactPoint {
            point: [1.0, 2.0, 3.0],
            normal: [0.0, 1.0, 0.0],
            depth: 0.1,
        });
        rotate_manifold_points(&mut m, identity);
        assert!((m.contacts[0].point[0] - 1.0).abs() < 1e-10);
        assert!((m.contacts[0].point[1] - 2.0).abs() < 1e-10);
        assert!((m.contacts[0].point[2] - 3.0).abs() < 1e-10);
    }

    #[test]
    fn test_rotate_manifold_points_90_about_z() {
        // 90° rotation about Z: (1,0,0) → (0,1,0)
        let rot_z_90 = [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]];
        let mut m = ContactManifold::new([1.0, 0.0, 0.0]);
        m.contacts.push(ContactPoint {
            point: [1.0, 0.0, 0.0],
            normal: [1.0, 0.0, 0.0],
            depth: 0.1,
        });
        rotate_manifold_points(&mut m, rot_z_90);
        assert!(
            (m.contacts[0].point[0]).abs() < 1e-10,
            "x should be 0 after rotation"
        );
        assert!(
            (m.contacts[0].point[1] - 1.0).abs() < 1e-10,
            "y should be 1 after rotation"
        );
    }

    // ── friction_basis ────────────────────────────────────────────────────────

    #[test]
    fn test_friction_basis_orthogonal_to_normal() {
        let normals: [[f64; 3]; 4] = [
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            normalize3([1.0, 1.0, 1.0]),
        ];
        for n in normals {
            let (t1, t2) = friction_basis(n);
            assert!(
                dot3(n, t1).abs() < 1e-9,
                "t1 should be perpendicular to normal"
            );
            assert!(
                dot3(n, t2).abs() < 1e-9,
                "t2 should be perpendicular to normal"
            );
            assert!(
                dot3(t1, t2).abs() < 1e-9,
                "t1 and t2 should be perpendicular"
            );
            let l1 = (t1[0] * t1[0] + t1[1] * t1[1] + t1[2] * t1[2]).sqrt();
            let l2 = (t2[0] * t2[0] + t2[1] * t2[1] + t2[2] * t2[2]).sqrt();
            assert!((l1 - 1.0).abs() < 1e-9, "t1 should be unit length");
            assert!((l2 - 1.0).abs() < 1e-9, "t2 should be unit length");
        }
    }

    // ── decompose_velocity ────────────────────────────────────────────────────

    #[test]
    fn test_decompose_velocity_normal_component() {
        let vel = [0.0, -3.0, 0.0];
        let normal = [0.0, 1.0, 0.0];
        let (vn, vt) = decompose_velocity(vel, normal);
        assert!(
            (vn[1] - (-3.0)).abs() < 1e-10,
            "normal component should be [0,-3,0]"
        );
        assert!(
            vt[0].abs() < 1e-10 && vt[1].abs() < 1e-10 && vt[2].abs() < 1e-10,
            "tangential component should be zero for pure normal velocity"
        );
    }

    #[test]
    fn test_decompose_velocity_tangential_component() {
        let vel = [2.0, 0.0, 0.0];
        let normal = [0.0, 1.0, 0.0];
        let (vn, vt) = decompose_velocity(vel, normal);
        assert!(
            vn[0].abs() < 1e-10 && vn[1].abs() < 1e-10,
            "normal component should be zero"
        );
        assert!(
            (vt[0] - 2.0).abs() < 1e-10,
            "tangential component should be [2,0,0]"
        );
    }

    // ── coulomb_friction_limit ────────────────────────────────────────────────

    #[test]
    fn test_coulomb_friction_limit_basic() {
        let limit = coulomb_friction_limit(10.0, 0.5);
        assert!((limit - 5.0).abs() < 1e-10, "Expected 5.0, got {limit}");
    }

    #[test]
    fn test_coulomb_friction_limit_negative_normal() {
        // Negative normal impulse should clamp to zero
        let limit = coulomb_friction_limit(-1.0, 0.5);
        assert_eq!(
            limit, 0.0,
            "Negative normal impulse should give 0 friction limit"
        );
    }

    // ── manifold_quality ──────────────────────────────────────────────────────

    #[test]
    fn test_manifold_quality_empty() {
        let m = ContactManifold::new([0.0, 1.0, 0.0]);
        let q = manifold_quality(&m);
        assert_eq!(q, 0.0, "Empty manifold should have quality 0");
    }

    #[test]
    fn test_manifold_quality_single_point() {
        let mut m = ContactManifold::new([0.0, 1.0, 0.0]);
        m.contacts.push(ContactPoint {
            point: [0.0; 3],
            normal: [0.0, 1.0, 0.0],
            depth: 0.1,
        });
        let q = manifold_quality(&m);
        assert!(q > 0.0 && q <= 1.0, "Quality should be in (0,1], got {q}");
    }

    #[test]
    fn test_manifold_quality_four_wide_points() {
        let mut m = ContactManifold::new([0.0, 1.0, 0.0]);
        for &x in &[-1.0f64, 1.0, 0.0, 0.0] {
            m.contacts.push(ContactPoint {
                point: [x, 0.0, 0.0],
                normal: [0.0, 1.0, 0.0],
                depth: 0.5,
            });
        }
        let q = manifold_quality(&m);
        assert!(q > 0.0 && q <= 1.0, "Quality should be in (0,1], got {q}");
    }
}
