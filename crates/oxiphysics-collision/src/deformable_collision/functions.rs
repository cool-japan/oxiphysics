//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::{
    ConstraintParams, DeformBvh, DeformCcdResult, DeformCollisionConfig, DeformConstraint,
    DeformContact, DeformContactType, DeformManifold, DeformableMesh, PenaltyParams, RigidBox,
    RigidPlane, RigidSphere, TriFace,
};

/// Subtract two 3-vectors.
#[inline]
pub(super) fn v3_sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
/// Add two 3-vectors.
#[inline]
pub(super) fn v3_add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}
/// Scale a 3-vector.
#[inline]
pub(super) fn v3_scale(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}
/// Dot product of two 3-vectors.
#[inline]
pub(super) fn v3_dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
/// Cross product of two 3-vectors.
#[inline]
pub(super) fn v3_cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
/// Squared length of a 3-vector.
#[inline]
pub(super) fn v3_len_sq(a: [f64; 3]) -> f64 {
    v3_dot(a, a)
}
/// Length of a 3-vector.
#[inline]
pub(super) fn v3_len(a: [f64; 3]) -> f64 {
    v3_len_sq(a).sqrt()
}
/// Normalise a 3-vector. Returns zero vector for near-zero input.
#[inline]
pub(super) fn v3_norm(a: [f64; 3]) -> [f64; 3] {
    let l = v3_len(a);
    if l < 1e-30 {
        [0.0; 3]
    } else {
        v3_scale(a, 1.0 / l)
    }
}
/// Component-wise minimum.
#[inline]
pub(super) fn v3_min(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0].min(b[0]), a[1].min(b[1]), a[2].min(b[2])]
}
/// Component-wise maximum.
#[inline]
pub(super) fn v3_max(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0].max(b[0]), a[1].max(b[1]), a[2].max(b[2])]
}
/// Negate a 3-vector.
#[inline]
pub(super) fn v3_neg(a: [f64; 3]) -> [f64; 3] {
    [-a[0], -a[1], -a[2]]
}
/// Linear interpolation between two 3-vectors.
#[inline]
pub(super) fn v3_lerp(a: [f64; 3], b: [f64; 3], t: f64) -> [f64; 3] {
    v3_add(v3_scale(a, 1.0 - t), v3_scale(b, t))
}
/// Distance squared between two points.
#[inline]
pub(super) fn v3_dist_sq(a: [f64; 3], b: [f64; 3]) -> f64 {
    v3_len_sq(v3_sub(b, a))
}
/// Default tolerance for geometric tests.
pub(super) const EPSILON: f64 = 1e-10;
/// Maximum iterations for CCD bisection.
pub(super) const CCD_MAX_ITER: usize = 64;
/// Default spatial hash cell size.
pub(super) const DEFAULT_CELL_SIZE: f64 = 1.0;
/// Compute the closest point on triangle (a,b,c) to point p.
/// Returns (closest_point, barycentric_coords, squared_distance).
pub fn closest_point_on_triangle(
    p: [f64; 3],
    a: [f64; 3],
    b: [f64; 3],
    c: [f64; 3],
) -> ([f64; 3], [f64; 3], f64) {
    let ab = v3_sub(b, a);
    let ac = v3_sub(c, a);
    let ap = v3_sub(p, a);
    let d1 = v3_dot(ab, ap);
    let d2 = v3_dot(ac, ap);
    if d1 <= 0.0 && d2 <= 0.0 {
        return (a, [1.0, 0.0, 0.0], v3_dist_sq(p, a));
    }
    let bp = v3_sub(p, b);
    let d3 = v3_dot(ab, bp);
    let d4 = v3_dot(ac, bp);
    if d3 >= 0.0 && d4 <= d3 {
        return (b, [0.0, 1.0, 0.0], v3_dist_sq(p, b));
    }
    let vc = d1 * d4 - d3 * d2;
    if vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0 {
        let v = d1 / (d1 - d3);
        let pt = v3_add(a, v3_scale(ab, v));
        return (pt, [1.0 - v, v, 0.0], v3_dist_sq(p, pt));
    }
    let cp = v3_sub(p, c);
    let d5 = v3_dot(ab, cp);
    let d6 = v3_dot(ac, cp);
    if d6 >= 0.0 && d5 <= d6 {
        return (c, [0.0, 0.0, 1.0], v3_dist_sq(p, c));
    }
    let vb = d5 * d2 - d1 * d6;
    if vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0 {
        let w = d2 / (d2 - d6);
        let pt = v3_add(a, v3_scale(ac, w));
        return (pt, [1.0 - w, 0.0, w], v3_dist_sq(p, pt));
    }
    let va = d3 * d6 - d5 * d4;
    if va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0 {
        let w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        let pt = v3_add(b, v3_scale(v3_sub(c, b), w));
        return (pt, [0.0, 1.0 - w, w], v3_dist_sq(p, pt));
    }
    let denom = 1.0 / (va + vb + vc);
    let v = vb * denom;
    let w = vc * denom;
    let pt = v3_add(a, v3_add(v3_scale(ab, v), v3_scale(ac, w)));
    (pt, [1.0 - v - w, v, w], v3_dist_sq(p, pt))
}
/// Test vertex-face proximity. Returns a contact if the vertex is within
/// `thickness` of the triangle face.
pub fn vertex_face_contact(
    vertex: [f64; 3],
    a: [f64; 3],
    b: [f64; 3],
    c: [f64; 3],
    thickness: f64,
    feature_a: usize,
    feature_b: usize,
) -> Option<DeformContact> {
    let (closest, bary, dist_sq) = closest_point_on_triangle(vertex, a, b, c);
    let thresh_sq = thickness * thickness;
    if dist_sq > thresh_sq || dist_sq < EPSILON * EPSILON {
        return None;
    }
    let dist = dist_sq.sqrt();
    let normal = v3_norm(v3_sub(vertex, closest));
    let depth = thickness - dist;
    if depth <= 0.0 {
        return None;
    }
    Some(DeformContact::new(
        closest,
        normal,
        depth,
        bary,
        DeformContactType::VertexFace,
        feature_a,
        feature_b,
    ))
}
/// Closest points between two line segments (p0,p1) and (q0,q1).
/// Returns (s, t, closest_on_seg1, closest_on_seg2, dist_sq).
pub fn closest_points_edge_edge(
    p0: [f64; 3],
    p1: [f64; 3],
    q0: [f64; 3],
    q1: [f64; 3],
) -> (f64, f64, [f64; 3], [f64; 3], f64) {
    let d1 = v3_sub(p1, p0);
    let d2 = v3_sub(q1, q0);
    let r = v3_sub(p0, q0);
    let a = v3_dot(d1, d1);
    let e = v3_dot(d2, d2);
    let f = v3_dot(d2, r);
    if a < EPSILON && e < EPSILON {
        return (0.0, 0.0, p0, q0, v3_dist_sq(p0, q0));
    }
    let (mut s, mut t);
    if a < EPSILON {
        s = 0.0;
        t = (f / e).clamp(0.0, 1.0);
    } else {
        let c = v3_dot(d1, r);
        if e < EPSILON {
            t = 0.0;
            s = (-c / a).clamp(0.0, 1.0);
        } else {
            let b = v3_dot(d1, d2);
            let denom = a * e - b * b;
            if denom.abs() > EPSILON {
                s = ((b * f - c * e) / denom).clamp(0.0, 1.0);
            } else {
                s = 0.0;
            }
            t = (b * s + f) / e;
            if t < 0.0 {
                t = 0.0;
                s = (-c / a).clamp(0.0, 1.0);
            } else if t > 1.0 {
                t = 1.0;
                s = ((b - c) / a).clamp(0.0, 1.0);
            }
        }
    }
    let cp = v3_add(p0, v3_scale(d1, s));
    let cq = v3_add(q0, v3_scale(d2, t));
    (s, t, cp, cq, v3_dist_sq(cp, cq))
}
/// Test edge-edge proximity. Returns a contact if the edges are within
/// `thickness` of each other.
pub fn edge_edge_contact(
    p0: [f64; 3],
    p1: [f64; 3],
    q0: [f64; 3],
    q1: [f64; 3],
    thickness: f64,
    feature_a: usize,
    feature_b: usize,
) -> Option<DeformContact> {
    let (_s, _t, cp, cq, dist_sq) = closest_points_edge_edge(p0, p1, q0, q1);
    let thresh_sq = thickness * thickness;
    if dist_sq > thresh_sq || dist_sq < EPSILON * EPSILON {
        return None;
    }
    let dist = dist_sq.sqrt();
    let normal = v3_norm(v3_sub(cp, cq));
    let depth = thickness - dist;
    if depth <= 0.0 {
        return None;
    }
    let mid = v3_scale(v3_add(cp, cq), 0.5);
    Some(DeformContact::new(
        mid,
        normal,
        depth,
        [0.0; 3],
        DeformContactType::EdgeEdge,
        feature_a,
        feature_b,
    ))
}
/// Detect all vertex-face contacts between two deformable meshes.
pub fn detect_vertex_face_contacts(
    mesh_a: &DeformableMesh,
    mesh_b: &DeformableMesh,
    bvh_a: &DeformBvh,
    bvh_b: &DeformBvh,
    config: &DeformCollisionConfig,
) -> Vec<DeformContact> {
    let mut contacts = Vec::new();
    let pairs = bvh_a.query_overlap_pairs_margin(bvh_b, config.thickness);
    for (fi_a, fi_b) in pairs {
        let fb = &mesh_b.faces[fi_b];
        let ta = mesh_b.positions[fb.v0];
        let tb = mesh_b.positions[fb.v1];
        let tc = mesh_b.positions[fb.v2];
        let fa = &mesh_a.faces[fi_a];
        for &vi in &[fa.v0, fa.v1, fa.v2] {
            let p = mesh_a.positions[vi];
            if let Some(c) = vertex_face_contact(p, ta, tb, tc, config.thickness, vi, fi_b) {
                contacts.push(c);
            }
        }
        let fa_ref = &mesh_a.faces[fi_a];
        let ua = mesh_a.positions[fa_ref.v0];
        let ub = mesh_a.positions[fa_ref.v1];
        let uc = mesh_a.positions[fa_ref.v2];
        for &vi in &[fb.v0, fb.v1, fb.v2] {
            let p = mesh_b.positions[vi];
            if let Some(mut c) = vertex_face_contact(p, ua, ub, uc, config.thickness, vi, fi_a) {
                c.normal = v3_neg(c.normal);
                contacts.push(c);
            }
        }
    }
    contacts
}
/// Detect all edge-edge contacts between two deformable meshes.
pub fn detect_edge_edge_contacts(
    mesh_a: &DeformableMesh,
    mesh_b: &DeformableMesh,
    bvh_a: &DeformBvh,
    bvh_b: &DeformBvh,
    config: &DeformCollisionConfig,
) -> Vec<DeformContact> {
    let mut contacts = Vec::new();
    if !config.edge_edge {
        return contacts;
    }
    let pairs = bvh_a.query_overlap_pairs_margin(bvh_b, config.thickness);
    for (fi_a, fi_b) in pairs {
        let fa = &mesh_a.faces[fi_a];
        let fb = &mesh_b.faces[fi_b];
        let edges_a = fa.edges();
        let edges_b = fb.edges();
        for (ea_idx, &(ea0, ea1)) in edges_a.iter().enumerate() {
            for (eb_idx, &(eb0, eb1)) in edges_b.iter().enumerate() {
                let p0 = mesh_a.positions[ea0];
                let p1 = mesh_a.positions[ea1];
                let q0 = mesh_b.positions[eb0];
                let q1 = mesh_b.positions[eb1];
                if let Some(c) = edge_edge_contact(
                    p0,
                    p1,
                    q0,
                    q1,
                    config.thickness,
                    fi_a * 3 + ea_idx,
                    fi_b * 3 + eb_idx,
                ) {
                    contacts.push(c);
                }
            }
        }
    }
    contacts
}
/// Full deformable-deformable collision detection.
pub fn detect_deformable_deformable(
    mesh_a: &DeformableMesh,
    mesh_b: &DeformableMesh,
    bvh_a: &DeformBvh,
    bvh_b: &DeformBvh,
    config: &DeformCollisionConfig,
) -> Vec<DeformContact> {
    let mut contacts = detect_vertex_face_contacts(mesh_a, mesh_b, bvh_a, bvh_b, config);
    let ee = detect_edge_edge_contacts(mesh_a, mesh_b, bvh_a, bvh_b, config);
    contacts.extend(ee);
    contacts
}
/// Detect contacts between a deformable mesh and a rigid sphere.
pub fn detect_deformable_sphere(
    mesh: &DeformableMesh,
    sphere: &RigidSphere,
    thickness: f64,
) -> Vec<DeformContact> {
    let mut contacts = Vec::new();
    let total_r = sphere.radius + thickness;
    let total_r_sq = total_r * total_r;
    for (vi, &pos) in mesh.positions.iter().enumerate() {
        let dist_sq = v3_dist_sq(pos, sphere.centre);
        if dist_sq < total_r_sq && dist_sq > EPSILON * EPSILON {
            let dist = dist_sq.sqrt();
            let normal = v3_norm(v3_sub(pos, sphere.centre));
            let depth = total_r - dist;
            if depth > 0.0 {
                let contact_pos = v3_add(sphere.centre, v3_scale(normal, sphere.radius));
                contacts.push(DeformContact::new(
                    contact_pos,
                    normal,
                    depth,
                    [0.0; 3],
                    DeformContactType::VertexRigid,
                    vi,
                    0,
                ));
            }
        }
    }
    contacts
}
/// Detect contacts between a deformable mesh and a rigid plane.
pub fn detect_deformable_plane(
    mesh: &DeformableMesh,
    plane: &RigidPlane,
    thickness: f64,
) -> Vec<DeformContact> {
    let mut contacts = Vec::new();
    for (vi, &pos) in mesh.positions.iter().enumerate() {
        let signed_dist = v3_dot(pos, plane.normal) - plane.offset;
        if signed_dist < thickness {
            let depth = thickness - signed_dist;
            let contact_pos = v3_sub(pos, v3_scale(plane.normal, signed_dist));
            contacts.push(DeformContact::new(
                contact_pos,
                plane.normal,
                depth,
                [0.0; 3],
                DeformContactType::VertexRigid,
                vi,
                0,
            ));
        }
    }
    contacts
}
/// Closest point on an AABB to a point.
pub(super) fn closest_point_on_aabb(p: [f64; 3], bmin: [f64; 3], bmax: [f64; 3]) -> [f64; 3] {
    [
        p[0].clamp(bmin[0], bmax[0]),
        p[1].clamp(bmin[1], bmax[1]),
        p[2].clamp(bmin[2], bmax[2]),
    ]
}
/// Detect contacts between a deformable mesh and a rigid box.
pub fn detect_deformable_box(
    mesh: &DeformableMesh,
    rigid_box: &RigidBox,
    thickness: f64,
) -> Vec<DeformContact> {
    let mut contacts = Vec::new();
    for (vi, &pos) in mesh.positions.iter().enumerate() {
        let cp = closest_point_on_aabb(pos, rigid_box.min, rigid_box.max);
        let _dist_sq = v3_dist_sq(pos, cp);
        let inside = pos[0] >= rigid_box.min[0]
            && pos[0] <= rigid_box.max[0]
            && pos[1] >= rigid_box.min[1]
            && pos[1] <= rigid_box.max[1]
            && pos[2] >= rigid_box.min[2]
            && pos[2] <= rigid_box.max[2];
        if inside {
            let mut best_dist = f64::MAX;
            let mut best_normal = [0.0, 1.0, 0.0];
            let axes: [[f64; 3]; 6] = [
                [1.0, 0.0, 0.0],
                [-1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [0.0, -1.0, 0.0],
                [0.0, 0.0, 1.0],
                [0.0, 0.0, -1.0],
            ];
            let face_dists = [
                rigid_box.max[0] - pos[0],
                pos[0] - rigid_box.min[0],
                rigid_box.max[1] - pos[1],
                pos[1] - rigid_box.min[1],
                rigid_box.max[2] - pos[2],
                pos[2] - rigid_box.min[2],
            ];
            for (i, &fd) in face_dists.iter().enumerate() {
                if fd < best_dist {
                    best_dist = fd;
                    best_normal = axes[i];
                }
            }
            contacts.push(DeformContact::new(
                pos,
                best_normal,
                best_dist + thickness,
                [0.0; 3],
                DeformContactType::VertexRigid,
                vi,
                0,
            ));
        } else if _dist_sq < thickness * thickness {
            let dist = _dist_sq.sqrt();
            let normal = v3_norm(v3_sub(pos, cp));
            let depth = thickness - dist;
            if depth > 0.0 {
                contacts.push(DeformContact::new(
                    cp,
                    normal,
                    depth,
                    [0.0; 3],
                    DeformContactType::VertexRigid,
                    vi,
                    0,
                ));
            }
        }
    }
    contacts
}
/// Build a vertex-to-face adjacency map. Returns for each vertex the set of
/// face indices that share it.
pub fn build_adjacency(mesh: &DeformableMesh) -> Vec<Vec<usize>> {
    let mut adj = vec![Vec::new(); mesh.num_vertices()];
    for (fi, f) in mesh.faces.iter().enumerate() {
        adj[f.v0].push(fi);
        adj[f.v1].push(fi);
        adj[f.v2].push(fi);
    }
    adj
}
/// Check if two faces share a vertex.
pub(super) fn faces_share_vertex(fa: &TriFace, fb: &TriFace) -> bool {
    let va = [fa.v0, fa.v1, fa.v2];
    let vb = [fb.v0, fb.v1, fb.v2];
    for &a in &va {
        for &b in &vb {
            if a == b {
                return true;
            }
        }
    }
    false
}
/// Detect self-collision within a single deformable mesh using BVH.
pub fn detect_self_collision(
    mesh: &DeformableMesh,
    bvh: &DeformBvh,
    config: &DeformCollisionConfig,
) -> Vec<DeformContact> {
    let mut contacts = Vec::new();
    if !config.self_collision {
        return contacts;
    }
    let pairs = bvh.query_self_overlap();
    let adjacency = if config.adjacency_skip {
        Some(build_adjacency(mesh))
    } else {
        None
    };
    for (fi_a, fi_b) in pairs {
        let fa = &mesh.faces[fi_a];
        let fb = &mesh.faces[fi_b];
        if config.adjacency_skip && faces_share_vertex(fa, fb) {
            continue;
        }
        let tb0 = mesh.positions[fb.v0];
        let tb1 = mesh.positions[fb.v1];
        let tb2 = mesh.positions[fb.v2];
        for &vi in &[fa.v0, fa.v1, fa.v2] {
            if config.adjacency_skip
                && let Some(ref adj) = adjacency
                && adj[vi].contains(&fi_b)
            {
                continue;
            }
            let p = mesh.positions[vi];
            if let Some(c) = vertex_face_contact(p, tb0, tb1, tb2, config.thickness, vi, fi_b) {
                contacts.push(c);
            }
        }
        let ta0 = mesh.positions[fa.v0];
        let ta1 = mesh.positions[fa.v1];
        let ta2 = mesh.positions[fa.v2];
        for &vi in &[fb.v0, fb.v1, fb.v2] {
            if config.adjacency_skip
                && let Some(ref adj) = adjacency
                && adj[vi].contains(&fi_a)
            {
                continue;
            }
            let p = mesh.positions[vi];
            if let Some(c) = vertex_face_contact(p, ta0, ta1, ta2, config.thickness, vi, fi_a) {
                contacts.push(c);
            }
        }
        if config.edge_edge {
            let edges_a = fa.edges();
            let edges_b = fb.edges();
            for (ea_idx, &(ea0, ea1)) in edges_a.iter().enumerate() {
                for (eb_idx, &(eb0, eb1)) in edges_b.iter().enumerate() {
                    if ea0 == eb0 || ea0 == eb1 || ea1 == eb0 || ea1 == eb1 {
                        continue;
                    }
                    let p0 = mesh.positions[ea0];
                    let p1 = mesh.positions[ea1];
                    let q0 = mesh.positions[eb0];
                    let q1 = mesh.positions[eb1];
                    if let Some(c) = edge_edge_contact(
                        p0,
                        p1,
                        q0,
                        q1,
                        config.thickness,
                        fi_a * 3 + ea_idx,
                        fi_b * 3 + eb_idx,
                    ) {
                        contacts.push(c);
                    }
                }
            }
        }
    }
    contacts
}
/// Number of uniform sub-intervals for the coarse CCD pre-scan.
///
/// A naive monotone bisection fails when the distance-vs-time curve dips below
/// `thickness` mid-path while both endpoints stay above it (e.g. a thin shell
/// the primitive passes straight through). The coarse pre-scan locates the
/// earliest sub-interval that brackets such a crossing so the subsequent
/// bisection converges to the true earliest time of impact.
pub(super) const CCD_PRESCAN_SUBDIVISIONS: usize = 128;
/// Coarse uniform pre-scan followed by bisection refinement.
///
/// `dist` maps a time `t` in `[0, 1]` to the (non-squared) distance between the
/// two moving primitives. Returns the earliest time in `[0, 1]` whose distance
/// is strictly below `thickness`, or `None` if the primitives never come within
/// `thickness` over the interval.
///
/// The caller is responsible for the `t == 0` early-out; this routine scans
/// sub-interval right endpoints `t = k/N` for `k` in `1..=N` and refines the
/// first bracket that dips below `thickness`. If no sample dips below, it also
/// refines around the interior minimum sample to catch a thin dip that falls
/// entirely between two coarse samples.
pub(super) fn ccd_prescan_bisect(thickness: f64, dist: impl Fn(f64) -> f64) -> Option<f64> {
    let n = CCD_PRESCAN_SUBDIVISIONS;
    let inv_n = 1.0 / n as f64;
    // Coarse uniform pre-scan. Track the minimum sample for the fallback. The
    // `t = 0` sample is the caller's early-out responsibility but is still
    // included here so an interior minimum is detected relative to it.
    let mut min_dist = dist(0.0);
    let mut kmin = 0usize;
    for k in 1..=n {
        let t = k as f64 * inv_n;
        let d = dist(t);
        if d < min_dist {
            min_dist = d;
            kmin = k;
        }
        if d < thickness {
            // First sub-interval whose right endpoint dips below thickness:
            // [lo, hi] brackets a crossing because the left endpoint sample at
            // (k-1)/N is, by construction of this left-to-right scan, still
            // >= thickness.
            let lo = (k - 1) as f64 * inv_n;
            let hi = t;
            return Some(refine_crossing(lo, hi, thickness, &dist));
        }
    }
    // No coarse sample dipped below thickness. Guard against a thin dip that
    // lives entirely inside a sub-interval around the global minimum.
    if (1..=n - 1).contains(&kmin) {
        let lo = (kmin - 1) as f64 * inv_n;
        let hi = (kmin + 1) as f64 * inv_n;
        if let Some(t) = refine_minimum(lo, hi, thickness, &dist) {
            return Some(t);
        }
    }
    None
}
/// Bisect `[lo, hi]` where `dist(lo) >= thickness` and `dist(hi) < thickness`,
/// converging on the earliest time whose distance is below `thickness`.
fn refine_crossing(mut lo: f64, mut hi: f64, thickness: f64, dist: impl Fn(f64) -> f64) -> f64 {
    let mut best_t = hi;
    for _ in 0..CCD_MAX_ITER {
        let mid = (lo + hi) * 0.5;
        if dist(mid) < thickness {
            hi = mid;
            best_t = mid;
        } else {
            lo = mid;
        }
        if (hi - lo) < 1e-8 {
            break;
        }
    }
    best_t
}
/// Narrow `[lo, hi]` toward an interior minimum looking for any sub-thickness
/// time. Returns the earliest such time, or `None` if the dip never reaches
/// below `thickness`. Endpoints are not assumed to bracket a crossing.
fn refine_minimum(
    mut lo: f64,
    mut hi: f64,
    thickness: f64,
    dist: impl Fn(f64) -> f64,
) -> Option<f64> {
    let mut found: Option<f64> = None;
    let mut d_lo = dist(lo);
    let mut d_hi = dist(hi);
    for _ in 0..CCD_MAX_ITER {
        let mid = (lo + hi) * 0.5;
        let d_mid = dist(mid);
        if d_mid < thickness {
            // Sub-thickness time located. Both `lo` and `hi` are still
            // >= thickness here (the initial endpoints came from above-thickness
            // pre-scan samples and every endpoint update moved to an
            // above-thickness `mid`), so `[lo, mid]` brackets a crossing.
            // Refine it for the earliest sub-thickness time.
            found = Some(refine_crossing(lo, mid, thickness, &dist));
            break;
        }
        // Move the endpoint with the larger distance inward toward the smaller,
        // walking the bracket toward the minimum.
        if d_lo >= d_hi {
            lo = mid;
            d_lo = d_mid;
        } else {
            hi = mid;
            d_hi = d_mid;
        }
        if (hi - lo) < 1e-8 {
            break;
        }
    }
    found
}
/// Vertex-face CCD using bisection: does the moving vertex cross the moving
/// triangle during \[0, 1\]?
pub fn ccd_vertex_face(
    v0: [f64; 3],
    v1: [f64; 3],
    a0: [f64; 3],
    a1: [f64; 3],
    b0: [f64; 3],
    b1: [f64; 3],
    c0: [f64; 3],
    c1: [f64; 3],
    thickness: f64,
) -> Option<f64> {
    let (_, _, dist0) = closest_point_on_triangle(v0, a0, b0, c0);
    if dist0.sqrt() < thickness {
        return Some(0.0);
    }
    // Coarse uniform pre-scan + bisection: a thin shell the vertex passes
    // straight through has both endpoints above `thickness` while the true
    // minimum is mid-path, so a naive monotone bisection would walk away from
    // it. The pre-scan brackets the earliest crossing before refining.
    ccd_prescan_bisect(thickness, |t| {
        let vm = v3_lerp(v0, v1, t);
        let am = v3_lerp(a0, a1, t);
        let bm = v3_lerp(b0, b1, t);
        let cm = v3_lerp(c0, c1, t);
        let (_, _, dsq) = closest_point_on_triangle(vm, am, bm, cm);
        dsq.sqrt()
    })
}
/// Edge-edge CCD using bisection.
pub fn ccd_edge_edge(
    p0_start: [f64; 3],
    p0_end: [f64; 3],
    p1_start: [f64; 3],
    p1_end: [f64; 3],
    q0_start: [f64; 3],
    q0_end: [f64; 3],
    q1_start: [f64; 3],
    q1_end: [f64; 3],
    thickness: f64,
) -> Option<f64> {
    let (_, _, _, _, dist0) = closest_points_edge_edge(p0_start, p1_start, q0_start, q1_start);
    if dist0.sqrt() < thickness {
        return Some(0.0);
    }
    // Coarse uniform pre-scan + bisection (see `ccd_vertex_face`): the edges may
    // sweep through each other with both endpoints above `thickness`, so bracket
    // the earliest crossing before refining.
    ccd_prescan_bisect(thickness, |t| {
        let pm0 = v3_lerp(p0_start, p0_end, t);
        let pm1 = v3_lerp(p1_start, p1_end, t);
        let qm0 = v3_lerp(q0_start, q0_end, t);
        let qm1 = v3_lerp(q1_start, q1_end, t);
        let (_, _, _, _, dsq) = closest_points_edge_edge(pm0, pm1, qm0, qm1);
        dsq.sqrt()
    })
}
/// Vertex-face CCD via analytic cubic solver, falling back to bisection on miss
/// for numerical robustness. Arg order matches `ccd_vertex_face`.
pub fn ccd_vertex_face_analytic(
    v0: [f64; 3],
    v1: [f64; 3],
    a0: [f64; 3],
    a1: [f64; 3],
    b0: [f64; 3],
    b1: [f64; 3],
    c0: [f64; 3],
    c1: [f64; 3],
    thickness: f64,
) -> Option<f64> {
    let vv = v3_sub(v1, v0);
    let av = v3_sub(a1, a0);
    let bv = v3_sub(b1, b0);
    let cv = v3_sub(c1, c0);
    if let Some(t) =
        super::cubic_toi::cubic_toi_vertex_face(a0, b0, c0, v0, av, bv, cv, vv, thickness)
    {
        return Some(t);
    }
    ccd_vertex_face(v0, v1, a0, a1, b0, b1, c0, c1, thickness)
}

/// Edge-edge CCD via analytic cubic solver, falling back to bisection on miss.
/// Edge 1 = p0->p1, Edge 2 = q0->q1. Arg order matches `ccd_edge_edge`.
pub fn ccd_edge_edge_analytic(
    p0_start: [f64; 3],
    p0_end: [f64; 3],
    p1_start: [f64; 3],
    p1_end: [f64; 3],
    q0_start: [f64; 3],
    q0_end: [f64; 3],
    q1_start: [f64; 3],
    q1_end: [f64; 3],
    thickness: f64,
) -> Option<f64> {
    let v_p0 = v3_sub(p0_end, p0_start);
    let v_p1 = v3_sub(p1_end, p1_start);
    let v_q0 = v3_sub(q0_end, q0_start);
    let v_q1 = v3_sub(q1_end, q1_start);
    if let Some(t) = super::cubic_toi::cubic_toi_edge_edge(
        p0_start, p1_start, q0_start, q1_start, v_p0, v_p1, v_q0, v_q1, thickness,
    ) {
        return Some(t);
    }
    ccd_edge_edge(
        p0_start, p0_end, p1_start, p1_end, q0_start, q0_end, q1_start, q1_end, thickness,
    )
}
/// Run CCD over an entire deformable mesh pair. Returns the earliest contact.
pub fn ccd_deformable_deformable(
    mesh_a_prev: &DeformableMesh,
    mesh_a_curr: &DeformableMesh,
    mesh_b_prev: &DeformableMesh,
    mesh_b_curr: &DeformableMesh,
    bvh_a: &DeformBvh,
    bvh_b: &DeformBvh,
    thickness: f64,
) -> Option<DeformCcdResult> {
    let pairs = bvh_a.query_overlap_pairs(bvh_b);
    let mut earliest: Option<DeformCcdResult> = None;
    for (fi_a, fi_b) in pairs {
        let fa = &mesh_a_curr.faces[fi_a];
        let fb = &mesh_b_curr.faces[fi_b];
        for &vi in &[fa.v0, fa.v1, fa.v2] {
            let vp = mesh_a_prev.positions[vi];
            let vc = mesh_a_curr.positions[vi];
            let a_p = mesh_b_prev.positions[fb.v0];
            let a_c = mesh_b_curr.positions[fb.v0];
            let b_p = mesh_b_prev.positions[fb.v1];
            let b_c = mesh_b_curr.positions[fb.v1];
            let c_p = mesh_b_prev.positions[fb.v2];
            let c_c = mesh_b_curr.positions[fb.v2];
            if let Some(toi) = ccd_vertex_face(vp, vc, a_p, a_c, b_p, b_c, c_p, c_c, thickness) {
                let update = match &earliest {
                    Some(e) => toi < e.toi,
                    None => true,
                };
                if update {
                    let pos = v3_lerp(vp, vc, toi);
                    let at = v3_lerp(a_p, a_c, toi);
                    let bt = v3_lerp(b_p, b_c, toi);
                    let ct = v3_lerp(c_p, c_c, toi);
                    let n = v3_norm(v3_cross(v3_sub(bt, at), v3_sub(ct, at)));
                    earliest = Some(DeformCcdResult {
                        toi,
                        normal: n,
                        position: pos,
                        feature_a: vi,
                        feature_b: fi_b,
                    });
                }
            }
        }
    }
    earliest
}
/// Compute penalty force for a single deformable contact.
/// Returns force to apply to the vertex (feature_a).
pub fn compute_penalty_force(
    contact: &DeformContact,
    velocity: [f64; 3],
    params: &PenaltyParams,
) -> [f64; 3] {
    let vn = v3_dot(velocity, contact.normal);
    let fn_mag = params.stiffness * contact.depth - params.damping * vn;
    let fn_mag = fn_mag.max(0.0);
    let f_normal = v3_scale(contact.normal, fn_mag);
    let vt = v3_sub(velocity, v3_scale(contact.normal, vn));
    let vt_len = v3_len(vt);
    let f_friction = if vt_len > EPSILON {
        let max_friction = params.friction * fn_mag;
        let friction_mag = max_friction.min(params.damping * vt_len);
        v3_scale(v3_norm(vt), -friction_mag)
    } else {
        [0.0; 3]
    };
    v3_add(f_normal, f_friction)
}
/// Apply penalty forces for all contacts to a deformable mesh.
/// Modifies velocities in-place.
pub fn apply_penalty_forces(
    mesh: &mut DeformableMesh,
    contacts: &[DeformContact],
    params: &PenaltyParams,
    dt: f64,
) {
    for c in contacts {
        let vi = c.feature_a;
        if vi >= mesh.num_vertices() {
            continue;
        }
        let vel = mesh.velocities[vi];
        let force = compute_penalty_force(c, vel, params);
        let inv_m = mesh.inv_mass[vi];
        if inv_m > 0.0 {
            let dv = v3_scale(force, inv_m * dt);
            mesh.velocities[vi] = v3_add(vel, dv);
        }
    }
}
/// Generate position constraints from contacts.
pub fn generate_constraints(contacts: &[DeformContact]) -> Vec<DeformConstraint> {
    contacts
        .iter()
        .map(|c| DeformConstraint {
            vertex: c.feature_a,
            target: v3_add(c.position, v3_scale(c.normal, c.depth)),
            normal: c.normal,
            depth: c.depth,
        })
        .collect()
}
/// Solve position constraints (Gauss-Seidel projection).
pub fn solve_constraints(
    mesh: &mut DeformableMesh,
    constraints: &[DeformConstraint],
    params: &ConstraintParams,
) {
    for _ in 0..params.iterations {
        for c in constraints {
            let vi = c.vertex;
            if vi >= mesh.num_vertices() || mesh.inv_mass[vi] <= 0.0 {
                continue;
            }
            let pos = mesh.positions[vi];
            let penetration = v3_dot(v3_sub(c.target, pos), c.normal);
            if penetration <= 0.0 {
                continue;
            }
            let correction = v3_scale(c.normal, penetration * params.relaxation);
            mesh.positions[vi] = v3_add(pos, correction);
        }
    }
}
/// Apply friction during constraint solve.
pub fn apply_constraint_friction(
    mesh: &mut DeformableMesh,
    constraints: &[DeformConstraint],
    params: &ConstraintParams,
    dt: f64,
) {
    if dt < EPSILON {
        return;
    }
    let inv_dt = 1.0 / dt;
    for c in constraints {
        let vi = c.vertex;
        if vi >= mesh.num_vertices() || mesh.inv_mass[vi] <= 0.0 {
            continue;
        }
        let vel = mesh.velocities[vi];
        let vn = v3_dot(vel, c.normal);
        let vt = v3_sub(vel, v3_scale(c.normal, vn));
        let vt_len = v3_len(vt);
        if vt_len < EPSILON {
            continue;
        }
        let max_friction = params.friction * (vn.abs() + c.depth * inv_dt);
        let friction_mag = max_friction.min(vt_len);
        let friction_dir = v3_norm(vt);
        mesh.velocities[vi] = v3_sub(vel, v3_scale(friction_dir, friction_mag));
    }
}
/// Full pipeline: detect contacts and generate manifold.
pub fn generate_contact_manifold(
    mesh_a: &DeformableMesh,
    mesh_b: &DeformableMesh,
    bvh_a: &DeformBvh,
    bvh_b: &DeformBvh,
    config: &DeformCollisionConfig,
    body_a_id: usize,
    body_b_id: usize,
    max_contacts: usize,
) -> DeformManifold {
    let raw_contacts = detect_deformable_deformable(mesh_a, mesh_b, bvh_a, bvh_b, config);
    let mut manifold = DeformManifold::new(body_a_id, body_b_id, max_contacts);
    for c in raw_contacts {
        manifold.add_contact(c);
    }
    manifold
}
/// Generate manifold for self-collision.
pub fn generate_self_collision_manifold(
    mesh: &DeformableMesh,
    bvh: &DeformBvh,
    config: &DeformCollisionConfig,
    body_id: usize,
    max_contacts: usize,
) -> DeformManifold {
    let raw_contacts = detect_self_collision(mesh, bvh, config);
    let mut manifold = DeformManifold::new(body_id, body_id, max_contacts);
    for c in raw_contacts {
        manifold.add_contact(c);
    }
    manifold
}
/// Create a simple quad mesh (two triangles) for testing.
#[cfg(test)]
pub(super) fn make_quad_mesh(origin: [f64; 3], size: f64) -> DeformableMesh {
    let p0 = origin;
    let p1 = v3_add(origin, [size, 0.0, 0.0]);
    let p2 = v3_add(origin, [size, 0.0, size]);
    let p3 = v3_add(origin, [0.0, 0.0, size]);
    let positions = vec![p0, p1, p2, p3];
    let faces = vec![TriFace::new(0, 1, 2), TriFace::new(0, 2, 3)];
    let inv_mass = vec![1.0; 4];
    DeformableMesh::new(positions, faces, inv_mass)
}
/// Create a tetrahedron mesh for testing.
#[cfg(test)]
pub(super) fn make_tetra_mesh(centre: [f64; 3], scale: f64) -> DeformableMesh {
    let positions = vec![
        v3_add(centre, v3_scale([1.0, 1.0, 1.0], scale)),
        v3_add(centre, v3_scale([1.0, -1.0, -1.0], scale)),
        v3_add(centre, v3_scale([-1.0, 1.0, -1.0], scale)),
        v3_add(centre, v3_scale([-1.0, -1.0, 1.0], scale)),
    ];
    let faces = vec![
        TriFace::new(0, 1, 2),
        TriFace::new(0, 1, 3),
        TriFace::new(0, 2, 3),
        TriFace::new(1, 2, 3),
    ];
    let inv_mass = vec![1.0; 4];
    DeformableMesh::new(positions, faces, inv_mass)
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::deformable_collision::DeformAabb;
    use crate::deformable_collision::SpatialHash;
    #[test]
    fn test_v3_sub() {
        let a = [3.0, 2.0, 1.0];
        let b = [1.0, 1.0, 1.0];
        let r = v3_sub(a, b);
        assert!((r[0] - 2.0).abs() < 1e-12);
        assert!((r[1] - 1.0).abs() < 1e-12);
        assert!((r[2]).abs() < 1e-12);
    }
    #[test]
    fn test_v3_add() {
        let r = v3_add([1.0, 2.0, 3.0], [4.0, 5.0, 6.0]);
        assert!((r[0] - 5.0).abs() < 1e-12);
        assert!((r[1] - 7.0).abs() < 1e-12);
        assert!((r[2] - 9.0).abs() < 1e-12);
    }
    #[test]
    fn test_v3_cross() {
        let x = [1.0, 0.0, 0.0];
        let y = [0.0, 1.0, 0.0];
        let z = v3_cross(x, y);
        assert!((z[0]).abs() < 1e-12);
        assert!((z[1]).abs() < 1e-12);
        assert!((z[2] - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_v3_norm_unit() {
        let a = [3.0, 0.0, 0.0];
        let n = v3_norm(a);
        assert!((n[0] - 1.0).abs() < 1e-12);
        assert!((v3_len(n) - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_v3_norm_zero() {
        let n = v3_norm([0.0; 3]);
        assert!(v3_len(n).abs() < 1e-12);
    }
    #[test]
    fn test_v3_lerp() {
        let a = [0.0, 0.0, 0.0];
        let b = [10.0, 20.0, 30.0];
        let m = v3_lerp(a, b, 0.5);
        assert!((m[0] - 5.0).abs() < 1e-12);
        assert!((m[1] - 10.0).abs() < 1e-12);
        assert!((m[2] - 15.0).abs() < 1e-12);
    }
    #[test]
    fn test_aabb_overlap() {
        let a = DeformAabb::new([0.0; 3], [2.0; 3]);
        let b = DeformAabb::new([1.0; 3], [3.0; 3]);
        assert!(a.overlaps(&b));
        assert!(b.overlaps(&a));
    }
    #[test]
    fn test_aabb_no_overlap() {
        let a = DeformAabb::new([0.0; 3], [1.0; 3]);
        let b = DeformAabb::new([2.0; 3], [3.0; 3]);
        assert!(!a.overlaps(&b));
    }
    #[test]
    fn test_aabb_expand_point() {
        let mut aabb = DeformAabb::empty();
        aabb.expand_point([1.0, 2.0, 3.0]);
        aabb.expand_point([-1.0, -2.0, -3.0]);
        assert!((aabb.min[0] - (-1.0)).abs() < 1e-12);
        assert!((aabb.max[1] - 2.0).abs() < 1e-12);
    }
    #[test]
    fn test_aabb_surface_area() {
        let a = DeformAabb::new([0.0; 3], [1.0, 2.0, 3.0]);
        let sa = a.surface_area();
        assert!((sa - 22.0).abs() < 1e-12);
    }
    #[test]
    fn test_aabb_from_triangle() {
        let aabb = DeformAabb::from_triangle([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
        assert!(aabb.min[0].abs() < 1e-12);
        assert!((aabb.max[0] - 1.0).abs() < 1e-12);
        assert!((aabb.max[1] - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_aabb_union() {
        let a = DeformAabb::new([0.0; 3], [1.0; 3]);
        let b = DeformAabb::new([-1.0; 3], [2.0; 3]);
        let u = a.union(&b);
        assert!((u.min[0] - (-1.0)).abs() < 1e-12);
        assert!((u.max[0] - 2.0).abs() < 1e-12);
    }
    #[test]
    fn test_aabb_centre() {
        let a = DeformAabb::new([0.0; 3], [2.0; 3]);
        let c = a.centre();
        assert!((c[0] - 1.0).abs() < 1e-12);
        assert!((c[1] - 1.0).abs() < 1e-12);
        assert!((c[2] - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_aabb_expand_margin() {
        let mut a = DeformAabb::new([0.0; 3], [1.0; 3]);
        a.expand_margin(0.5);
        assert!((a.min[0] - (-0.5)).abs() < 1e-12);
        assert!((a.max[0] - 1.5).abs() < 1e-12);
    }
    #[test]
    fn test_tri_face_edges() {
        let f = TriFace::new(0, 2, 1);
        let edges = f.edges();
        assert_eq!(edges[0], (0, 2));
        assert_eq!(edges[1], (1, 2));
        assert_eq!(edges[2], (0, 1));
    }
    #[test]
    fn test_mesh_creation() {
        let mesh = make_quad_mesh([0.0; 3], 1.0);
        assert_eq!(mesh.num_vertices(), 4);
        assert_eq!(mesh.num_faces(), 2);
    }
    #[test]
    fn test_mesh_face_normal() {
        let mesh = make_quad_mesh([0.0; 3], 1.0);
        let n = mesh.face_normal_unit(0);
        assert!(n[1].abs() > 0.9);
    }
    #[test]
    fn test_mesh_bounding_box() {
        let mesh = make_quad_mesh([0.0; 3], 2.0);
        let bb = mesh.bounding_box();
        assert!(bb.min[0].abs() < 1e-12);
        assert!((bb.max[0] - 2.0).abs() < 1e-12);
    }
    #[test]
    fn test_mesh_collect_edges() {
        let mesh = make_quad_mesh([0.0; 3], 1.0);
        let edges = mesh.collect_edges();
        assert_eq!(edges.len(), 5);
    }
    #[test]
    fn test_mesh_snapshot() {
        let mut mesh = make_quad_mesh([0.0; 3], 1.0);
        mesh.positions[0] = [99.0, 99.0, 99.0];
        mesh.snapshot_positions();
        assert!((mesh.prev_positions[0][0] - 99.0).abs() < 1e-12);
    }
    #[test]
    fn test_closest_point_inside() {
        let a = [0.0, 0.0, 0.0];
        let b = [2.0, 0.0, 0.0];
        let c = [0.0, 2.0, 0.0];
        let p = [0.5, 0.5, 1.0];
        let (cp, bary, _dsq) = closest_point_on_triangle(p, a, b, c);
        assert!(cp[2].abs() < 1e-12);
        assert!(bary[0] > 0.0 && bary[1] > 0.0 && bary[2] > 0.0);
    }
    #[test]
    fn test_closest_point_vertex() {
        let a = [0.0, 0.0, 0.0];
        let b = [1.0, 0.0, 0.0];
        let c = [0.0, 1.0, 0.0];
        let p = [-1.0, -1.0, 0.0];
        let (_cp, bary, _dsq) = closest_point_on_triangle(p, a, b, c);
        assert!((bary[0] - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_edge_edge_parallel() {
        let p0 = [0.0, 0.0, 0.0];
        let p1 = [1.0, 0.0, 0.0];
        let q0 = [0.0, 1.0, 0.0];
        let q1 = [1.0, 1.0, 0.0];
        let (_, _, _, _, dsq) = closest_points_edge_edge(p0, p1, q0, q1);
        assert!((dsq - 1.0).abs() < 1e-8);
    }
    #[test]
    fn test_edge_edge_crossing() {
        let p0 = [0.0, 0.0, -1.0];
        let p1 = [0.0, 0.0, 1.0];
        let q0 = [-1.0, 0.5, 0.0];
        let q1 = [1.0, 0.5, 0.0];
        let (_, _, _, _, dsq) = closest_points_edge_edge(p0, p1, q0, q1);
        assert!((dsq - 0.25).abs() < 1e-8);
    }
    #[test]
    fn test_bvh_build() {
        let mesh = make_quad_mesh([0.0; 3], 1.0);
        let bvh = DeformBvh::build(&mesh);
        assert!(!bvh.nodes.is_empty());
    }
    #[test]
    fn test_bvh_refit() {
        let mut mesh = make_quad_mesh([0.0; 3], 1.0);
        let mut bvh = DeformBvh::build(&mesh);
        mesh.positions[1] = [5.0, 0.0, 0.0];
        bvh.refit(&mesh);
        assert!(bvh.nodes[0].aabb.max[0] >= 5.0);
    }
    #[test]
    fn test_bvh_overlap_pairs() {
        let mesh_a = make_quad_mesh([0.0; 3], 1.0);
        let mesh_b = make_quad_mesh([0.5, 0.0, 0.5], 1.0);
        let bvh_a = DeformBvh::build(&mesh_a);
        let bvh_b = DeformBvh::build(&mesh_b);
        let pairs = bvh_a.query_overlap_pairs(&bvh_b);
        assert!(!pairs.is_empty());
    }
    #[test]
    fn test_bvh_no_overlap() {
        let mesh_a = make_quad_mesh([0.0; 3], 1.0);
        let mesh_b = make_quad_mesh([100.0, 100.0, 100.0], 1.0);
        let bvh_a = DeformBvh::build(&mesh_a);
        let bvh_b = DeformBvh::build(&mesh_b);
        let pairs = bvh_a.query_overlap_pairs(&bvh_b);
        assert!(pairs.is_empty());
    }
    #[test]
    fn test_self_overlap_tetra() {
        let mesh = make_tetra_mesh([0.0; 3], 1.0);
        let bvh = DeformBvh::build(&mesh);
        let pairs = bvh.query_self_overlap();
        for (a, b) in &pairs {
            assert!(a < b);
        }
    }
    #[test]
    fn test_self_collision_adjacency_skip() {
        let mesh = make_tetra_mesh([0.0; 3], 0.5);
        let bvh = DeformBvh::build(&mesh);
        let config = DeformCollisionConfig {
            thickness: 10.0,
            adjacency_skip: true,
            self_collision: true,
            edge_edge: false,
        };
        let contacts = detect_self_collision(&mesh, &bvh, &config);
        assert!(contacts.is_empty());
    }
    #[test]
    fn test_spatial_hash_insert_query() {
        let mesh = make_quad_mesh([0.0; 3], 1.0);
        let mut sh = SpatialHash::new(0.5);
        sh.build_from_mesh(&mesh);
        assert!(sh.num_cells() > 0);
        let aabb = DeformAabb::new([0.0; 3], [0.5; 3]);
        let result = sh.query_aabb(&aabb);
        assert!(!result.is_empty());
    }
    #[test]
    fn test_spatial_hash_self_pairs() {
        let mesh = make_quad_mesh([0.0; 3], 1.0);
        let mut sh = SpatialHash::new(2.0);
        sh.build_from_mesh(&mesh);
        let pairs = sh.query_self_pairs();
        assert!(!pairs.is_empty());
    }
    #[test]
    fn test_spatial_hash_clear() {
        let mesh = make_quad_mesh([0.0; 3], 1.0);
        let mut sh = SpatialHash::new(0.5);
        sh.build_from_mesh(&mesh);
        assert!(sh.num_cells() > 0);
        sh.clear();
        assert_eq!(sh.num_cells(), 0);
    }
    #[test]
    fn test_vertex_face_contact_hit() {
        let a = [0.0, 0.0, 0.0];
        let b = [2.0, 0.0, 0.0];
        let c = [0.0, 2.0, 0.0];
        let p = [0.5, 0.5, 0.005];
        let contact = vertex_face_contact(p, a, b, c, 0.01, 0, 0);
        assert!(contact.is_some());
        let ct = contact.unwrap();
        assert!(ct.depth > 0.0);
    }
    #[test]
    fn test_vertex_face_contact_miss() {
        let a = [0.0, 0.0, 0.0];
        let b = [2.0, 0.0, 0.0];
        let c = [0.0, 2.0, 0.0];
        let p = [0.5, 0.5, 5.0];
        let contact = vertex_face_contact(p, a, b, c, 0.01, 0, 0);
        assert!(contact.is_none());
    }
    #[test]
    fn test_edge_edge_contact_hit() {
        let p0 = [0.0, 0.0, -1.0];
        let p1 = [0.0, 0.0, 1.0];
        let q0 = [-1.0, 0.005, 0.0];
        let q1 = [1.0, 0.005, 0.0];
        let contact = edge_edge_contact(p0, p1, q0, q1, 0.01, 0, 0);
        assert!(contact.is_some());
    }
    #[test]
    fn test_edge_edge_contact_miss() {
        let p0 = [0.0, 0.0, -1.0];
        let p1 = [0.0, 0.0, 1.0];
        let q0 = [-1.0, 5.0, 0.0];
        let q1 = [1.0, 5.0, 0.0];
        let contact = edge_edge_contact(p0, p1, q0, q1, 0.01, 0, 0);
        assert!(contact.is_none());
    }
    #[test]
    fn test_deformable_plane_contact() {
        let mesh = make_quad_mesh([0.0, 0.01, 0.0], 1.0);
        let plane = RigidPlane {
            normal: [0.0, 1.0, 0.0],
            offset: 0.0,
        };
        let contacts = detect_deformable_plane(&mesh, &plane, 0.05);
        assert!(!contacts.is_empty());
        for c in &contacts {
            assert!(c.depth > 0.0);
        }
    }
    #[test]
    fn test_deformable_sphere_contact() {
        let mesh = make_quad_mesh([0.0; 3], 0.1);
        let sphere = RigidSphere {
            centre: [0.05, 0.5, 0.05],
            radius: 0.5,
        };
        let contacts = detect_deformable_sphere(&mesh, &sphere, 0.05);
        assert!(!contacts.is_empty());
    }
    #[test]
    fn test_deformable_sphere_no_contact() {
        let mesh = make_quad_mesh([0.0; 3], 1.0);
        let sphere = RigidSphere {
            centre: [0.5, 100.0, 0.5],
            radius: 1.0,
        };
        let contacts = detect_deformable_sphere(&mesh, &sphere, 0.01);
        assert!(contacts.is_empty());
    }
    #[test]
    fn test_deformable_box_inside() {
        let mut mesh = make_quad_mesh([0.0; 3], 0.1);
        mesh.positions[0] = [0.5, 0.5, 0.5];
        let rigid_box = RigidBox {
            min: [0.0; 3],
            max: [1.0; 3],
        };
        let contacts = detect_deformable_box(&mesh, &rigid_box, 0.01);
        assert!(!contacts.is_empty());
    }
    #[test]
    fn test_ccd_vertex_face_hit() {
        let v0 = [0.5, 1.0, 0.5];
        let v1 = [0.5, -1.0, 0.5];
        let a0 = [0.0, 0.0, 0.0];
        let a1 = [0.0, 0.0, 0.0];
        let b0 = [2.0, 0.0, 0.0];
        let b1 = [2.0, 0.0, 0.0];
        let c0 = [0.0, 0.0, 2.0];
        let c1 = [0.0, 0.0, 2.0];
        let toi = ccd_vertex_face(v0, v1, a0, a1, b0, b1, c0, c1, 0.01);
        assert!(toi.is_some());
        assert!(toi.unwrap() > 0.0 && toi.unwrap() < 1.0);
    }
    #[test]
    fn test_ccd_vertex_face_miss() {
        let v0 = [0.5, 1.0, 0.5];
        let v1 = [0.5, 2.0, 0.5];
        let a0 = [0.0, 0.0, 0.0];
        let a1 = [0.0, 0.0, 0.0];
        let b0 = [2.0, 0.0, 0.0];
        let b1 = [2.0, 0.0, 0.0];
        let c0 = [0.0, 0.0, 2.0];
        let c1 = [0.0, 0.0, 2.0];
        let toi = ccd_vertex_face(v0, v1, a0, a1, b0, b1, c0, c1, 0.01);
        assert!(toi.is_none());
    }
    #[test]
    fn test_ccd_edge_edge_hit() {
        let p0s = [0.0, 0.0, -1.0];
        let p0e = [0.0, 0.0, -1.0];
        let p1s = [0.0, 0.0, 1.0];
        let p1e = [0.0, 0.0, 1.0];
        let q0s = [-1.0, 1.0, 0.0];
        let q0e = [-1.0, -1.0, 0.0];
        let q1s = [1.0, 1.0, 0.0];
        let q1e = [1.0, -1.0, 0.0];
        let toi = ccd_edge_edge(p0s, p0e, p1s, p1e, q0s, q0e, q1s, q1e, 0.01);
        assert!(toi.is_some());
    }
    #[test]
    fn test_penalty_force_pushes_out() {
        let contact = DeformContact::new(
            [0.0; 3],
            [0.0, 1.0, 0.0],
            0.1,
            [0.0; 3],
            DeformContactType::VertexFace,
            0,
            0,
        );
        let vel = [0.0; 3];
        let params = PenaltyParams::default();
        let force = compute_penalty_force(&contact, vel, &params);
        assert!(force[1] > 0.0);
    }
    #[test]
    fn test_penalty_force_damping() {
        let contact = DeformContact::new(
            [0.0; 3],
            [0.0, 1.0, 0.0],
            0.1,
            [0.0; 3],
            DeformContactType::VertexFace,
            0,
            0,
        );
        let vel = [0.0, 10.0, 0.0];
        let params = PenaltyParams::default();
        let force_still = compute_penalty_force(&contact, [0.0; 3], &params);
        let force_moving = compute_penalty_force(&contact, vel, &params);
        assert!(force_moving[1] < force_still[1]);
    }
    #[test]
    fn test_apply_penalty_forces() {
        let mut mesh = make_quad_mesh([0.0; 3], 1.0);
        let contact = DeformContact::new(
            [0.0; 3],
            [0.0, 1.0, 0.0],
            0.1,
            [0.0; 3],
            DeformContactType::VertexFace,
            0,
            0,
        );
        let params = PenaltyParams::default();
        apply_penalty_forces(&mut mesh, &[contact], &params, 0.01);
        assert!(mesh.velocities[0][1] > 0.0);
    }
    #[test]
    fn test_constraint_solve() {
        let mut mesh = make_quad_mesh([0.0, -0.1, 0.0], 1.0);
        let constraint = DeformConstraint {
            vertex: 0,
            target: [0.0, 0.0, 0.0],
            normal: [0.0, 1.0, 0.0],
            depth: 0.1,
        };
        let params = ConstraintParams::default();
        solve_constraints(&mut mesh, &[constraint], &params);
        assert!(mesh.positions[0][1] > -0.1);
    }
    #[test]
    fn test_constraint_pinned_vertex() {
        let mut mesh = make_quad_mesh([0.0, -0.1, 0.0], 1.0);
        mesh.inv_mass[0] = 0.0;
        let constraint = DeformConstraint {
            vertex: 0,
            target: [0.0, 0.0, 0.0],
            normal: [0.0, 1.0, 0.0],
            depth: 0.1,
        };
        let params = ConstraintParams::default();
        let original_y = mesh.positions[0][1];
        solve_constraints(&mut mesh, &[constraint], &params);
        assert!((mesh.positions[0][1] - original_y).abs() < 1e-12);
    }
    #[test]
    fn test_manifold_add_reduce() {
        let mut m = DeformManifold::new(0, 1, 3);
        for i in 0..5 {
            m.add_contact(DeformContact::new(
                [0.0; 3],
                [0.0, 1.0, 0.0],
                i as f64 * 0.1,
                [0.0; 3],
                DeformContactType::VertexFace,
                i,
                0,
            ));
        }
        assert!(m.len() <= 3);
        assert!(m.max_depth() >= 0.3);
    }
    #[test]
    fn test_manifold_average_normal() {
        let mut m = DeformManifold::new(0, 1, 10);
        m.add_contact(DeformContact::new(
            [0.0; 3],
            [0.0, 1.0, 0.0],
            0.1,
            [0.0; 3],
            DeformContactType::VertexFace,
            0,
            0,
        ));
        m.add_contact(DeformContact::new(
            [0.0; 3],
            [0.0, 1.0, 0.0],
            0.1,
            [0.0; 3],
            DeformContactType::VertexFace,
            1,
            0,
        ));
        let avg = m.average_normal();
        assert!((avg[1] - 1.0).abs() < 1e-6);
    }
    #[test]
    fn test_manifold_is_empty() {
        let m = DeformManifold::new(0, 1, 10);
        assert!(m.is_empty());
    }
    #[test]
    fn test_generate_constraints() {
        let contacts = vec![DeformContact::new(
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            0.05,
            [0.0; 3],
            DeformContactType::VertexFace,
            0,
            0,
        )];
        let constraints = generate_constraints(&contacts);
        assert_eq!(constraints.len(), 1);
        assert_eq!(constraints[0].vertex, 0);
    }
    #[test]
    fn test_build_adjacency() {
        let mesh = make_tetra_mesh([0.0; 3], 1.0);
        let adj = build_adjacency(&mesh);
        assert_eq!(adj.len(), 4);
        for v in &adj {
            assert_eq!(v.len(), 3);
        }
    }
    #[test]
    fn test_faces_share_vertex_true() {
        let a = TriFace::new(0, 1, 2);
        let b = TriFace::new(0, 3, 4);
        assert!(faces_share_vertex(&a, &b));
    }
    #[test]
    fn test_faces_share_vertex_false() {
        let a = TriFace::new(0, 1, 2);
        let b = TriFace::new(3, 4, 5);
        assert!(!faces_share_vertex(&a, &b));
    }
    #[test]
    fn test_contact_type() {
        let c = DeformContact::new(
            [0.0; 3],
            [0.0, 1.0, 0.0],
            0.1,
            [0.0; 3],
            DeformContactType::EdgeEdge,
            0,
            0,
        );
        assert_eq!(c.contact_type, DeformContactType::EdgeEdge);
    }
    #[test]
    fn test_full_pipeline_overlapping_meshes() {
        let mesh_a = make_quad_mesh([0.0; 3], 1.0);
        let mesh_b = make_quad_mesh([0.0, 0.005, 0.0], 1.0);
        let bvh_a = DeformBvh::build(&mesh_a);
        let bvh_b = DeformBvh::build(&mesh_b);
        let config = DeformCollisionConfig {
            thickness: 0.02,
            edge_edge: false,
            self_collision: false,
            adjacency_skip: false,
        };
        let manifold =
            generate_contact_manifold(&mesh_a, &mesh_b, &bvh_a, &bvh_b, &config, 0, 1, 32);
        assert!(!manifold.is_empty());
    }
    #[test]
    fn test_full_pipeline_separated_meshes() {
        let mesh_a = make_quad_mesh([0.0; 3], 1.0);
        let mesh_b = make_quad_mesh([0.0, 100.0, 0.0], 1.0);
        let bvh_a = DeformBvh::build(&mesh_a);
        let bvh_b = DeformBvh::build(&mesh_b);
        let config = DeformCollisionConfig::default();
        let manifold =
            generate_contact_manifold(&mesh_a, &mesh_b, &bvh_a, &bvh_b, &config, 0, 1, 32);
        assert!(manifold.is_empty());
    }
    #[test]
    fn test_friction_response() {
        let mut mesh = make_quad_mesh([0.0; 3], 1.0);
        mesh.velocities[0] = [1.0, 0.0, 0.0];
        let constraint = DeformConstraint {
            vertex: 0,
            target: [0.0, 0.05, 0.0],
            normal: [0.0, 1.0, 0.0],
            depth: 0.05,
        };
        let params = ConstraintParams {
            iterations: 1,
            relaxation: 1.0,
            friction: 0.5,
        };
        apply_constraint_friction(&mut mesh, &[constraint], &params, 0.01);
        assert!(mesh.velocities[0][0] < 1.0);
    }
}
