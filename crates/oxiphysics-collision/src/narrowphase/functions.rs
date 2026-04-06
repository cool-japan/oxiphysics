//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions_2::ray_vs_capsule;
use super::types::{
    CompoundShape, ContactFeature, ContactFilter, FeatureContact, NarrowPhaseContact,
    PointQueryResult, RayCastResult, ShapeKind, TriangleMesh,
};
pub(super) fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
pub(super) fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
pub(super) fn add3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}
pub(super) fn scale3(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}
pub(super) fn len3(a: [f64; 3]) -> f64 {
    dot3(a, a).sqrt()
}
pub(super) fn normalize3(a: [f64; 3]) -> [f64; 3] {
    let l = len3(a);
    if l > 1e-12 {
        scale3(a, 1.0 / l)
    } else {
        [0.0, 1.0, 0.0]
    }
}
#[allow(dead_code)]
pub(super) fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
/// Top-level routing function: compute a contact between two world-space shapes.
///
/// Picks the fastest specialized algorithm when the shape combination is known,
/// and returns `None` if the shapes are separated.
#[allow(dead_code)]
pub fn shape_shape_contact(a: &ShapeKind, b: &ShapeKind) -> Option<NarrowPhaseContact> {
    match (a, b) {
        (
            ShapeKind::Sphere {
                center: ca,
                radius: ra,
            },
            ShapeKind::Sphere {
                center: cb,
                radius: rb,
            },
        ) => sphere_vs_sphere(*ca, *ra, *cb, *rb),
        (ShapeKind::Sphere { center, radius }, ShapeKind::Plane { normal, offset }) => {
            sphere_vs_plane(*center, *radius, *normal, *offset)
        }
        (ShapeKind::Plane { normal, offset }, ShapeKind::Sphere { center, radius }) => {
            sphere_vs_plane(*center, *radius, *normal, *offset).map(|c| c.flipped())
        }
        (
            ShapeKind::Sphere { center, radius },
            ShapeKind::Box {
                center: bc,
                half_extents,
            },
        ) => sphere_vs_box(*center, *radius, *bc, *half_extents),
        (
            ShapeKind::Box {
                center: bc,
                half_extents,
            },
            ShapeKind::Sphere { center, radius },
        ) => sphere_vs_box(*center, *radius, *bc, *half_extents).map(|c| c.flipped()),
        (
            ShapeKind::Capsule {
                p0: a0,
                p1: a1,
                radius: ra,
            },
            ShapeKind::Capsule {
                p0: b0,
                p1: b1,
                radius: rb,
            },
        ) => capsule_vs_capsule(*a0, *a1, *ra, *b0, *b1, *rb),
        (ShapeKind::Sphere { center, radius }, ShapeKind::Capsule { p0, p1, radius: rc }) => {
            sphere_vs_capsule(*center, *radius, *p0, *p1, *rc)
        }
        (ShapeKind::Capsule { p0, p1, radius: rc }, ShapeKind::Sphere { center, radius }) => {
            sphere_vs_capsule(*center, *radius, *p0, *p1, *rc).map(|c| c.flipped())
        }
        _ => gjk_fallback(a, b),
    }
}
pub(super) fn sphere_vs_sphere(
    ca: [f64; 3],
    ra: f64,
    cb: [f64; 3],
    rb: f64,
) -> Option<NarrowPhaseContact> {
    let diff = sub3(cb, ca);
    let dist = len3(diff);
    let depth = ra + rb - dist;
    if depth < 0.0 {
        return None;
    }
    let normal = if dist > 1e-12 {
        scale3(diff, 1.0 / dist)
    } else {
        [0.0, 1.0, 0.0]
    };
    Some(NarrowPhaseContact {
        normal,
        depth,
        point_a: add3(ca, scale3(normal, ra)),
        point_b: sub3(cb, scale3(normal, rb)),
    })
}
pub(super) fn sphere_vs_plane(
    center: [f64; 3],
    radius: f64,
    plane_normal: [f64; 3],
    plane_offset: f64,
) -> Option<NarrowPhaseContact> {
    let dist = dot3(center, plane_normal) - plane_offset;
    let depth = radius - dist;
    if depth < 0.0 {
        return None;
    }
    let point_a = sub3(center, scale3(plane_normal, radius));
    let point_b = sub3(center, scale3(plane_normal, dist));
    Some(NarrowPhaseContact {
        normal: plane_normal,
        depth,
        point_a,
        point_b,
    })
}
pub(super) fn sphere_vs_box(
    sphere_center: [f64; 3],
    radius: f64,
    box_center: [f64; 3],
    half_extents: [f64; 3],
) -> Option<NarrowPhaseContact> {
    let local = sub3(sphere_center, box_center);
    let closest = [
        local[0].clamp(-half_extents[0], half_extents[0]),
        local[1].clamp(-half_extents[1], half_extents[1]),
        local[2].clamp(-half_extents[2], half_extents[2]),
    ];
    let diff = sub3(local, closest);
    let dist_sq = dot3(diff, diff);
    let r2 = radius * radius;
    if dist_sq > r2 {
        return None;
    }
    let dist = dist_sq.sqrt();
    let (normal, depth) = if dist < 1e-12 {
        let overlaps = [
            (
                half_extents[0] - local[0].abs(),
                if local[0] >= 0.0 {
                    [1.0, 0.0, 0.0]
                } else {
                    [-1.0, 0.0, 0.0]
                },
            ),
            (
                half_extents[1] - local[1].abs(),
                if local[1] >= 0.0 {
                    [0.0, 1.0, 0.0]
                } else {
                    [0.0, -1.0, 0.0]
                },
            ),
            (
                half_extents[2] - local[2].abs(),
                if local[2] >= 0.0 {
                    [0.0, 0.0, 1.0]
                } else {
                    [0.0, 0.0, -1.0]
                },
            ),
        ];
        let (pen, axis) = overlaps
            .iter()
            .copied()
            .min_by(|(a, _), (b, _)| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal))
            .unwrap_or((0.0, [0.0; 3]));
        (axis, pen + radius)
    } else {
        (normalize3(diff), radius - dist)
    };
    let point_b_world = add3(box_center, closest);
    let point_a = sub3(sphere_center, scale3(normal, radius));
    Some(NarrowPhaseContact {
        normal,
        depth,
        point_a,
        point_b: point_b_world,
    })
}
pub(super) fn capsule_vs_capsule(
    a0: [f64; 3],
    a1: [f64; 3],
    ra: f64,
    b0: [f64; 3],
    b1: [f64; 3],
    rb: f64,
) -> Option<NarrowPhaseContact> {
    let (pa, pb) = closest_segment_segment(a0, a1, b0, b1);
    let diff = sub3(pb, pa);
    let dist = len3(diff);
    let depth = ra + rb - dist;
    if depth < 0.0 {
        return None;
    }
    let normal = if dist > 1e-12 {
        scale3(diff, 1.0 / dist)
    } else {
        [0.0, 1.0, 0.0]
    };
    Some(NarrowPhaseContact {
        normal,
        depth,
        point_a: add3(pa, scale3(normal, ra)),
        point_b: sub3(pb, scale3(normal, rb)),
    })
}
pub(super) fn sphere_vs_capsule(
    sphere_center: [f64; 3],
    rs: f64,
    cap_p0: [f64; 3],
    cap_p1: [f64; 3],
    rc: f64,
) -> Option<NarrowPhaseContact> {
    let closest = closest_point_on_segment(sphere_center, cap_p0, cap_p1);
    let diff = sub3(sphere_center, closest);
    let dist = len3(diff);
    let depth = rs + rc - dist;
    if depth < 0.0 {
        return None;
    }
    let normal = if dist > 1e-12 {
        scale3(diff, 1.0 / dist)
    } else {
        [0.0, 1.0, 0.0]
    };
    Some(NarrowPhaseContact {
        normal,
        depth,
        point_a: sub3(sphere_center, scale3(normal, rs)),
        point_b: add3(closest, scale3(normal, rc)),
    })
}
/// Very simple GJK-style fallback: check AABB overlap as a quick rejection,
/// then test support points to see whether the shapes overlap.
pub(super) fn gjk_fallback(a: &ShapeKind, b: &ShapeKind) -> Option<NarrowPhaseContact> {
    let (a_min, a_max) = a.aabb();
    let (b_min, b_max) = b.aabb();
    for i in 0..3 {
        if a_max[i] < b_min[i] || b_max[i] < a_min[i] {
            return None;
        }
    }
    let axes: [[f64; 3]; 3] = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    let mut min_depth = f64::INFINITY;
    let mut best_normal = [0.0_f64, 1.0, 0.0];
    for axis in axes {
        let sa = dot3(a.support(axis), axis) - dot3(a.support(scale3(axis, -1.0)), axis);
        let sb = dot3(b.support(axis), axis) - dot3(b.support(scale3(axis, -1.0)), axis);
        let _ = (sa, sb);
        let proj_a_max = dot3(a.support(axis), axis);
        let proj_a_min = dot3(a.support(scale3(axis, -1.0)), scale3(axis, -1.0));
        let proj_b_max = dot3(b.support(axis), axis);
        let proj_b_min = dot3(b.support(scale3(axis, -1.0)), scale3(axis, -1.0));
        let overlap = proj_a_max.min(proj_b_max) - (-proj_a_min).max(-proj_b_min);
        if overlap < 0.0 {
            return None;
        }
        if overlap < min_depth {
            min_depth = overlap;
            best_normal = axis;
        }
    }
    let ca = a.support(best_normal);
    let cb = b.support(scale3(best_normal, -1.0));
    Some(NarrowPhaseContact {
        normal: best_normal,
        depth: min_depth,
        point_a: ca,
        point_b: cb,
    })
}
/// Closest point on segment \[p, q\] to point `v`.
pub(super) fn closest_point_on_segment(v: [f64; 3], p: [f64; 3], q: [f64; 3]) -> [f64; 3] {
    let pq = sub3(q, p);
    let pv = sub3(v, p);
    let denom = dot3(pq, pq);
    if denom < 1e-24 {
        return p;
    }
    let t = (dot3(pv, pq) / denom).clamp(0.0, 1.0);
    add3(p, scale3(pq, t))
}
/// Closest pair of points between two segments.
pub(super) fn closest_segment_segment(
    a0: [f64; 3],
    a1: [f64; 3],
    b0: [f64; 3],
    b1: [f64; 3],
) -> ([f64; 3], [f64; 3]) {
    let d1 = sub3(a1, a0);
    let d2 = sub3(b1, b0);
    let r = sub3(a0, b0);
    let a = dot3(d1, d1);
    let e = dot3(d2, d2);
    let f = dot3(d2, r);
    let (s, t) = if a <= 1e-12 && e <= 1e-12 {
        (0.0, 0.0)
    } else if a <= 1e-12 {
        (0.0, (f / e).clamp(0.0, 1.0))
    } else {
        let c = dot3(d1, r);
        if e <= 1e-12 {
            ((-c / a).clamp(0.0, 1.0), 0.0)
        } else {
            let b = dot3(d1, d2);
            let denom = a * e - b * b;
            let s0 = if denom > 1e-12 {
                ((b * f - c * e) / denom).clamp(0.0, 1.0)
            } else {
                0.0
            };
            let t0 = (b * s0 + f) / e;
            if t0 < 0.0 {
                ((-c / a).clamp(0.0, 1.0), 0.0)
            } else if t0 > 1.0 {
                (((b - c) / a).clamp(0.0, 1.0), 1.0)
            } else {
                (s0, t0)
            }
        }
    };
    (add3(a0, scale3(d1, s)), add3(b0, scale3(d2, t)))
}
/// Dispatcher that also records the contact feature.
///
/// Extends [`shape_shape_contact`] with feature tagging for faster warm-starting
/// and incremental manifold updates.
#[allow(dead_code)]
pub fn shape_shape_feature_contact(a: &ShapeKind, b: &ShapeKind) -> Option<FeatureContact> {
    let contact = shape_shape_contact(a, b)?;
    let feature = match (a, b) {
        (ShapeKind::Sphere { .. }, ShapeKind::Sphere { .. })
        | (ShapeKind::Sphere { .. }, ShapeKind::Capsule { .. })
        | (ShapeKind::Capsule { .. }, ShapeKind::Sphere { .. }) => {
            ContactFeature::VertexFace { vertex: 0, face: 0 }
        }
        (ShapeKind::Sphere { .. }, ShapeKind::Plane { .. })
        | (ShapeKind::Plane { .. }, ShapeKind::Sphere { .. }) => {
            ContactFeature::VertexFace { vertex: 0, face: 0 }
        }
        (ShapeKind::Sphere { .. }, ShapeKind::Box { .. })
        | (ShapeKind::Box { .. }, ShapeKind::Sphere { .. }) => {
            ContactFeature::VertexFace { vertex: 0, face: 0 }
        }
        (ShapeKind::Capsule { .. }, ShapeKind::Capsule { .. }) => ContactFeature::EdgeEdge {
            edge_a: 0,
            edge_b: 0,
        },
        _ => ContactFeature::Unknown,
    };
    Some(FeatureContact { contact, feature })
}
/// Route a convex shape against every triangle of a concave mesh and collect contacts.
///
/// Returns all contacts whose depth is ≥ `min_depth`.
#[allow(dead_code)]
pub fn convex_vs_mesh(
    convex: &ShapeKind,
    mesh: &TriangleMesh,
    min_depth: f64,
) -> Vec<NarrowPhaseContact> {
    let mut contacts = Vec::new();
    for i in 0..mesh.tri_count() {
        let [v0, v1, v2] = mesh.triangle(i);
        let centroid = scale3(add3(add3(v0, v1), v2), 1.0 / 3.0);
        let r0 = len3(sub3(v0, centroid));
        let r1 = len3(sub3(v1, centroid));
        let r2 = len3(sub3(v2, centroid));
        let tri_radius = r0.max(r1).max(r2);
        let tri_sphere = ShapeKind::Sphere {
            center: centroid,
            radius: tri_radius,
        };
        if let Some(c) = shape_shape_contact(convex, &tri_sphere)
            && c.depth >= min_depth
        {
            contacts.push(c);
        }
    }
    contacts
}
/// Run narrowphase between a compound shape and a single shape.
///
/// Returns one contact per colliding child–shape pair (shallowest first).
#[allow(dead_code)]
pub fn compound_vs_shape(
    compound: &CompoundShape,
    other: &ShapeKind,
    filter: &ContactFilter,
) -> Vec<NarrowPhaseContact> {
    let mut contacts: Vec<NarrowPhaseContact> = compound
        .children
        .iter()
        .filter_map(|child| {
            let c = shape_shape_contact(child, other)?;
            filter.apply(c)
        })
        .collect();
    contacts.sort_by(|a, b| {
        a.depth
            .partial_cmp(&b.depth)
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    contacts
}
/// Run narrowphase between two compound shapes.
///
/// Tests every child pair.  Returns all contacts passing the filter.
#[allow(dead_code)]
pub fn compound_vs_compound(
    a: &CompoundShape,
    b: &CompoundShape,
    filter: &ContactFilter,
) -> Vec<NarrowPhaseContact> {
    let mut contacts = Vec::new();
    for ca in &a.children {
        for cb in &b.children {
            if let Some(c) = shape_shape_contact(ca, cb)
                && let Some(fc) = filter.apply(c)
            {
                contacts.push(fc);
            }
        }
    }
    contacts
}
/// Test whether `point` is inside `shape` and return distance information.
#[allow(dead_code)]
pub fn point_query(point: [f64; 3], shape: &ShapeKind) -> PointQueryResult {
    match shape {
        ShapeKind::Sphere { center, radius } => {
            let diff = sub3(point, *center);
            let dist = len3(diff);
            let normal = if dist > 1e-12 {
                scale3(diff, 1.0 / dist)
            } else {
                [0.0, 1.0, 0.0]
            };
            let signed_dist = dist - radius;
            let closest = add3(*center, scale3(normal, *radius));
            PointQueryResult {
                is_inside: signed_dist < 0.0,
                closest_surface_point: closest,
                signed_distance: signed_dist,
                normal,
            }
        }
        ShapeKind::Box {
            center,
            half_extents,
        } => {
            let local = sub3(point, *center);
            let clamped = [
                local[0].clamp(-half_extents[0], half_extents[0]),
                local[1].clamp(-half_extents[1], half_extents[1]),
                local[2].clamp(-half_extents[2], half_extents[2]),
            ];
            let on_surface = add3(*center, clamped);
            let diff_to_point = sub3(point, on_surface);
            let dist = len3(diff_to_point);
            let is_inside = local[0].abs() <= half_extents[0]
                && local[1].abs() <= half_extents[1]
                && local[2].abs() <= half_extents[2];
            let signed_dist = if is_inside { -dist } else { dist };
            let normal = if dist > 1e-12 {
                scale3(diff_to_point, 1.0 / dist)
            } else {
                let overlaps = [
                    (
                        half_extents[0] - local[0].abs(),
                        if local[0] >= 0.0 {
                            [1.0_f64, 0.0, 0.0]
                        } else {
                            [-1.0, 0.0, 0.0]
                        },
                    ),
                    (
                        half_extents[1] - local[1].abs(),
                        if local[1] >= 0.0 {
                            [0.0, 1.0, 0.0]
                        } else {
                            [0.0, -1.0, 0.0]
                        },
                    ),
                    (
                        half_extents[2] - local[2].abs(),
                        if local[2] >= 0.0 {
                            [0.0, 0.0, 1.0]
                        } else {
                            [0.0, 0.0, -1.0]
                        },
                    ),
                ];
                overlaps
                    .iter()
                    .copied()
                    .min_by(|(a, _), (b, _)| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal))
                    .expect("overlaps is a fixed 3-element array")
                    .1
            };
            PointQueryResult {
                is_inside,
                closest_surface_point: on_surface,
                signed_distance: signed_dist,
                normal,
            }
        }
        ShapeKind::Capsule { p0, p1, radius } => {
            let seg = sub3(*p1, *p0);
            let pq = sub3(point, *p0);
            let denom = dot3(seg, seg);
            let t = if denom < 1e-24 {
                0.0
            } else {
                (dot3(pq, seg) / denom).clamp(0.0, 1.0)
            };
            let closest_axis = add3(*p0, scale3(seg, t));
            let diff = sub3(point, closest_axis);
            let dist = len3(diff);
            let normal = if dist > 1e-12 {
                scale3(diff, 1.0 / dist)
            } else {
                [0.0, 1.0, 0.0]
            };
            let signed_dist = dist - radius;
            let closest = add3(closest_axis, scale3(normal, *radius));
            PointQueryResult {
                is_inside: signed_dist < 0.0,
                closest_surface_point: closest,
                signed_distance: signed_dist,
                normal,
            }
        }
        ShapeKind::Plane { normal, offset } => {
            let n = normalize3(*normal);
            let dist = dot3(point, n) - offset;
            let closest = sub3(point, scale3(n, dist));
            PointQueryResult {
                is_inside: dist < 0.0,
                closest_surface_point: closest,
                signed_distance: dist,
                normal: n,
            }
        }
        ShapeKind::Convex { vertices } => {
            if vertices.is_empty() {
                return PointQueryResult {
                    is_inside: false,
                    closest_surface_point: point,
                    signed_distance: f64::INFINITY,
                    normal: [0.0, 1.0, 0.0],
                };
            }
            let closest_v = vertices
                .iter()
                .min_by(|a, b| {
                    len3(sub3(point, **a))
                        .partial_cmp(&len3(sub3(point, **b)))
                        .unwrap_or(std::cmp::Ordering::Equal)
                })
                .copied()
                .expect("vertices is non-empty (checked above)");
            let diff = sub3(point, closest_v);
            let dist = len3(diff);
            let normal = if dist > 1e-12 {
                scale3(diff, 1.0 / dist)
            } else {
                [0.0, 1.0, 0.0]
            };
            PointQueryResult {
                is_inside: false,
                closest_surface_point: closest_v,
                signed_distance: dist,
                normal,
            }
        }
    }
}
/// Cast a ray from `origin` in direction `dir` (not necessarily unit-length)
/// against `shape`.
///
/// Returns `Some(result)` if the ray hits within `[0, max_toi]`.
#[allow(dead_code)]
pub fn ray_cast(
    origin: [f64; 3],
    dir: [f64; 3],
    max_toi: f64,
    shape: &ShapeKind,
) -> Option<RayCastResult> {
    match shape {
        ShapeKind::Sphere { center, radius } => {
            ray_vs_sphere(origin, dir, max_toi, *center, *radius)
        }
        ShapeKind::Plane { normal, offset } => ray_vs_plane(origin, dir, max_toi, *normal, *offset),
        ShapeKind::Box {
            center,
            half_extents,
        } => ray_vs_aabb(origin, dir, max_toi, *center, *half_extents),
        ShapeKind::Capsule { p0, p1, radius } => {
            ray_vs_capsule(origin, dir, max_toi, *p0, *p1, *radius)
        }
        _ => {
            let br = shape.bounding_radius();
            let (mn, mx) = shape.aabb();
            let c = scale3(add3(mn, mx), 0.5);
            ray_vs_sphere(origin, dir, max_toi, c, br)
        }
    }
}
pub(super) fn ray_vs_sphere(
    o: [f64; 3],
    d: [f64; 3],
    max_toi: f64,
    center: [f64; 3],
    radius: f64,
) -> Option<RayCastResult> {
    let oc = sub3(o, center);
    let a = dot3(d, d);
    let b = 2.0 * dot3(oc, d);
    let c = dot3(oc, oc) - radius * radius;
    let disc = b * b - 4.0 * a * c;
    if disc < 0.0 {
        return None;
    }
    let sqrt_disc = disc.sqrt();
    let t = (-b - sqrt_disc) / (2.0 * a);
    let t = if t < 0.0 {
        (-b + sqrt_disc) / (2.0 * a)
    } else {
        t
    };
    if t < 0.0 || t > max_toi {
        return None;
    }
    let hit = add3(o, scale3(d, t));
    let normal = normalize3(sub3(hit, center));
    Some(RayCastResult {
        hit: true,
        toi: t,
        hit_point: hit,
        normal,
    })
}
pub(super) fn ray_vs_plane(
    o: [f64; 3],
    d: [f64; 3],
    max_toi: f64,
    normal: [f64; 3],
    offset: f64,
) -> Option<RayCastResult> {
    let n = normalize3(normal);
    let denom = dot3(d, n);
    if denom.abs() < 1e-12 {
        return None;
    }
    let t = (offset - dot3(o, n)) / denom;
    if t < 0.0 || t > max_toi {
        return None;
    }
    let hit = add3(o, scale3(d, t));
    Some(RayCastResult {
        hit: true,
        toi: t,
        hit_point: hit,
        normal: n,
    })
}
pub(super) fn ray_vs_aabb(
    o: [f64; 3],
    d: [f64; 3],
    max_toi: f64,
    center: [f64; 3],
    half: [f64; 3],
) -> Option<RayCastResult> {
    let mn = sub3(center, half);
    let mx = add3(center, half);
    let mut t_min = 0.0_f64;
    let mut t_max = max_toi;
    let mut hit_normal = [0.0_f64; 3];
    for i in 0..3 {
        if d[i].abs() < 1e-12 {
            if o[i] < mn[i] || o[i] > mx[i] {
                return None;
            }
        } else {
            let inv = 1.0 / d[i];
            let mut t1 = (mn[i] - o[i]) * inv;
            let mut t2 = (mx[i] - o[i]) * inv;
            let mut n = [-1.0_f64, 0.0, 0.0];
            if t1 > t2 {
                std::mem::swap(&mut t1, &mut t2);
                n[i] = 1.0;
            } else {
                n = [0.0; 3];
                n[i] = -1.0;
            }
            if t1 > t_min {
                t_min = t1;
                hit_normal = n;
            }
            if t2 < t_max {
                t_max = t2;
            }
            if t_min > t_max {
                return None;
            }
        }
    }
    if t_min < 0.0 {
        return None;
    }
    let hit = add3(o, scale3(d, t_min));
    Some(RayCastResult {
        hit: true,
        toi: t_min,
        hit_point: hit,
        normal: hit_normal,
    })
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::BatchNarrowPhase;
    #[test]
    fn sphere_aabb_correct() {
        let s = ShapeKind::Sphere {
            center: [1.0, 2.0, 3.0],
            radius: 2.0,
        };
        let (mn, mx) = s.aabb();
        assert!((mn[0] - (-1.0)).abs() < 1e-12);
        assert!((mx[0] - 3.0).abs() < 1e-12);
    }
    #[test]
    fn box_aabb_correct() {
        let b = ShapeKind::Box {
            center: [0.0; 3],
            half_extents: [1.0, 2.0, 3.0],
        };
        let (mn, mx) = b.aabb();
        assert!((mn[0] - (-1.0)).abs() < 1e-12);
        assert!((mx[1] - 2.0).abs() < 1e-12);
        assert!((mx[2] - 3.0).abs() < 1e-12);
    }
    #[test]
    fn capsule_aabb_correct() {
        let c = ShapeKind::Capsule {
            p0: [0.0, -1.0, 0.0],
            p1: [0.0, 1.0, 0.0],
            radius: 0.5,
        };
        let (mn, mx) = c.aabb();
        assert!((mn[1] - (-1.5)).abs() < 1e-12);
        assert!((mx[1] - 1.5).abs() < 1e-12);
    }
    #[test]
    fn sphere_sphere_overlapping_contact() {
        let a = ShapeKind::Sphere {
            center: [0.0; 3],
            radius: 1.0,
        };
        let b = ShapeKind::Sphere {
            center: [1.5, 0.0, 0.0],
            radius: 1.0,
        };
        let c = shape_shape_contact(&a, &b).expect("should collide");
        assert!(c.depth > 0.0, "depth={}", c.depth);
        assert!((c.depth - 0.5).abs() < 1e-10, "depth={}", c.depth);
    }
    #[test]
    fn sphere_sphere_separated_none() {
        let a = ShapeKind::Sphere {
            center: [0.0; 3],
            radius: 0.5,
        };
        let b = ShapeKind::Sphere {
            center: [5.0, 0.0, 0.0],
            radius: 0.5,
        };
        assert!(shape_shape_contact(&a, &b).is_none());
    }
    #[test]
    fn sphere_plane_touching() {
        let s = ShapeKind::Sphere {
            center: [0.0, 1.0, 0.0],
            radius: 1.0,
        };
        let p = ShapeKind::Plane {
            normal: [0.0, 1.0, 0.0],
            offset: 0.0,
        };
        let c = shape_shape_contact(&s, &p).expect("touching");
        assert!(c.depth >= 0.0, "depth={}", c.depth);
    }
    #[test]
    fn sphere_plane_above_no_contact() {
        let s = ShapeKind::Sphere {
            center: [0.0, 2.0, 0.0],
            radius: 0.5,
        };
        let p = ShapeKind::Plane {
            normal: [0.0, 1.0, 0.0],
            offset: 0.0,
        };
        assert!(shape_shape_contact(&s, &p).is_none());
    }
    #[test]
    fn sphere_box_overlapping() {
        let s = ShapeKind::Sphere {
            center: [0.9, 0.0, 0.0],
            radius: 0.5,
        };
        let b = ShapeKind::Box {
            center: [0.0; 3],
            half_extents: [1.0; 3],
        };
        let c = shape_shape_contact(&s, &b).expect("should overlap");
        assert!(c.depth > 0.0, "depth={}", c.depth);
    }
    #[test]
    fn sphere_box_separated() {
        let s = ShapeKind::Sphere {
            center: [3.0, 0.0, 0.0],
            radius: 0.5,
        };
        let b = ShapeKind::Box {
            center: [0.0; 3],
            half_extents: [1.0; 3],
        };
        assert!(shape_shape_contact(&s, &b).is_none());
    }
    #[test]
    fn capsule_capsule_overlap() {
        let a = ShapeKind::Capsule {
            p0: [0.0, 0.0, 0.0],
            p1: [0.0, 2.0, 0.0],
            radius: 0.5,
        };
        let b = ShapeKind::Capsule {
            p0: [0.8, 0.0, 0.0],
            p1: [0.8, 2.0, 0.0],
            radius: 0.5,
        };
        let c = shape_shape_contact(&a, &b).expect("should overlap");
        assert!(c.depth > 0.0, "depth={}", c.depth);
    }
    #[test]
    fn sphere_capsule_overlap() {
        let s = ShapeKind::Sphere {
            center: [0.0, 1.0, 0.0],
            radius: 0.5,
        };
        let c_shape = ShapeKind::Capsule {
            p0: [0.3, 0.0, 0.0],
            p1: [0.3, 2.0, 0.0],
            radius: 0.5,
        };
        let c = shape_shape_contact(&s, &c_shape).expect("should overlap");
        assert!(c.depth > 0.0);
    }
    #[test]
    fn contact_flipped_swaps_normal_and_points() {
        let c = NarrowPhaseContact {
            normal: [1.0, 0.0, 0.0],
            depth: 0.1,
            point_a: [1.0, 0.0, 0.0],
            point_b: [2.0, 0.0, 0.0],
        };
        let f = c.flipped();
        assert!((f.normal[0] - (-1.0)).abs() < 1e-12);
        assert_eq!(f.point_a, c.point_b);
        assert_eq!(f.point_b, c.point_a);
    }
    #[test]
    fn contact_midpoint() {
        let c = NarrowPhaseContact {
            normal: [0.0, 1.0, 0.0],
            depth: 0.2,
            point_a: [0.0, 0.0, 0.0],
            point_b: [2.0, 0.0, 0.0],
        };
        let mid = c.midpoint();
        assert!((mid[0] - 1.0).abs() < 1e-12);
    }
    #[test]
    fn filter_discards_shallow() {
        let f = ContactFilter {
            min_depth: 0.05,
            ..Default::default()
        };
        let c = NarrowPhaseContact {
            normal: [0.0, 1.0, 0.0],
            depth: 0.01,
            point_a: [0.0; 3],
            point_b: [0.0; 3],
        };
        assert!(f.apply(c).is_none());
    }
    #[test]
    fn filter_clamps_deep() {
        let f = ContactFilter {
            max_depth: 0.5,
            ..Default::default()
        };
        let c = NarrowPhaseContact {
            normal: [0.0, 1.0, 0.0],
            depth: 2.0,
            point_a: [0.0; 3],
            point_b: [0.0; 3],
        };
        let out = f.apply(c).unwrap();
        assert!((out.depth - 0.5).abs() < 1e-12);
    }
    #[test]
    fn filter_flips_normal() {
        let f = ContactFilter {
            flip_normal: true,
            ..Default::default()
        };
        let c = NarrowPhaseContact {
            normal: [1.0, 0.0, 0.0],
            depth: 0.1,
            point_a: [0.0; 3],
            point_b: [0.0; 3],
        };
        let out = f.apply(c).unwrap();
        assert!((out.normal[0] - (-1.0)).abs() < 1e-12);
    }
    #[test]
    fn batch_run_returns_correct_count() {
        let shapes = vec![
            ShapeKind::Sphere {
                center: [0.0; 3],
                radius: 1.0,
            },
            ShapeKind::Sphere {
                center: [1.5, 0.0, 0.0],
                radius: 1.0,
            },
            ShapeKind::Sphere {
                center: [10.0, 0.0, 0.0],
                radius: 0.5,
            },
        ];
        let pairs = vec![(0, 1), (0, 2), (1, 2)];
        let batch = BatchNarrowPhase::new();
        let results = batch.run(&shapes, &pairs);
        assert_eq!(results.len(), 3);
        assert!(results[0].is_some(), "0-1 should collide");
        assert!(results[1].is_none(), "0-2 too far");
        assert!(results[2].is_none(), "1-2 too far");
    }
    #[test]
    fn batch_run_compact_only_contacts() {
        let shapes = vec![
            ShapeKind::Sphere {
                center: [0.0; 3],
                radius: 1.0,
            },
            ShapeKind::Sphere {
                center: [1.5, 0.0, 0.0],
                radius: 1.0,
            },
            ShapeKind::Sphere {
                center: [10.0, 0.0, 0.0],
                radius: 0.5,
            },
        ];
        let pairs = vec![(0, 1), (0, 2)];
        let batch = BatchNarrowPhase::new();
        let results = batch.run_compact(&shapes, &pairs);
        assert_eq!(results.len(), 1);
        assert_eq!(results[0].0, (0, 1));
    }
    #[test]
    fn batch_run_out_of_bounds_index_returns_none() {
        let shapes = vec![ShapeKind::Sphere {
            center: [0.0; 3],
            radius: 1.0,
        }];
        let pairs = vec![(0, 99)];
        let batch = BatchNarrowPhase::new();
        let results = batch.run(&shapes, &pairs);
        assert!(results[0].is_none());
    }
    #[test]
    fn sphere_support_is_center_plus_radius() {
        let s = ShapeKind::Sphere {
            center: [1.0, 0.0, 0.0],
            radius: 2.0,
        };
        let sup = s.support([1.0, 0.0, 0.0]);
        assert!((sup[0] - 3.0).abs() < 1e-10, "sup={sup:?}");
    }
    #[test]
    fn box_support_corner() {
        let b = ShapeKind::Box {
            center: [0.0; 3],
            half_extents: [1.0; 3],
        };
        let sup = b.support([1.0, 1.0, 1.0]);
        assert!((sup[0] - 1.0).abs() < 1e-12);
        assert!((sup[1] - 1.0).abs() < 1e-12);
        assert!((sup[2] - 1.0).abs() < 1e-12);
    }
    #[test]
    fn convex_bounding_radius_is_max_vertex_len() {
        let v = vec![[3.0, 0.0, 0.0], [0.0, 4.0, 0.0], [1.0, 1.0, 1.0]];
        let s = ShapeKind::Convex { vertices: v };
        let r = s.bounding_radius();
        assert!((r - 4.0).abs() < 1e-12, "r={r}");
    }
}
