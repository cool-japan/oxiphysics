//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

pub use super::specialized::*;

use super::functions::{
    add3, dot3, len3, normalize3, point_query, ray_cast, ray_vs_sphere, scale3, sub3,
};
use super::types::{PointQueryResult, RayCastResult, SegmentCastResult, ShapeKind};

pub(super) fn ray_vs_capsule(
    o: [f64; 3],
    d: [f64; 3],
    max_toi: f64,
    p0: [f64; 3],
    p1: [f64; 3],
    radius: f64,
) -> Option<RayCastResult> {
    let seg = sub3(p1, p0);
    let seg_len_sq = dot3(seg, seg);
    if seg_len_sq < 1e-24 {
        return ray_vs_sphere(o, d, max_toi, p0, radius);
    }
    let seg_n = scale3(seg, 1.0 / seg_len_sq.sqrt());
    let dp = sub3(o, p0);
    let dp_perp = sub3(dp, scale3(seg_n, dot3(dp, seg_n)));
    let d_perp = sub3(d, scale3(seg_n, dot3(d, seg_n)));
    let a = dot3(d_perp, d_perp);
    let b = 2.0 * dot3(dp_perp, d_perp);
    let c = dot3(dp_perp, dp_perp) - radius * radius;
    let mut best_t = f64::INFINITY;
    let mut best_n = [0.0_f64, 1.0, 0.0];
    if a > 1e-24 {
        let disc = b * b - 4.0 * a * c;
        if disc >= 0.0 {
            let sqrt_d = disc.sqrt();
            for sign in &[-1.0_f64, 1.0] {
                let t = (-b + sign * sqrt_d) / (2.0 * a);
                if t >= 0.0 && t <= max_toi {
                    let hit = add3(o, scale3(d, t));
                    let proj = dot3(sub3(hit, p0), seg_n);
                    let seg_len = seg_len_sq.sqrt();
                    if proj >= 0.0 && proj <= seg_len && t < best_t {
                        best_t = t;
                        let axis_pt = add3(p0, scale3(seg_n, proj));
                        best_n = normalize3(sub3(hit, axis_pt));
                    }
                }
            }
        }
    }
    for &cap_center in &[p0, p1] {
        if let Some(r) = ray_vs_sphere(o, d, max_toi, cap_center, radius)
            && r.toi < best_t
        {
            best_t = r.toi;
            best_n = r.normal;
        }
    }
    if best_t <= max_toi {
        let hit = add3(o, scale3(d, best_t));
        Some(RayCastResult {
            hit: true,
            toi: best_t,
            hit_point: hit,
            normal: best_n,
        })
    } else {
        None
    }
}
/// Cast a line segment from `start` to `end` against `shape`.
///
/// Returns `Some` if the segment intersects the shape.
pub fn segment_cast(
    start: [f64; 3],
    end: [f64; 3],
    shape: &ShapeKind,
) -> Option<SegmentCastResult> {
    let dir = sub3(end, start);
    let len = len3(dir);
    if len < 1e-24 {
        let pq = point_query(start, shape);
        return if pq.is_inside {
            Some(SegmentCastResult {
                hit: true,
                t: 0.0,
                hit_point: start,
                normal: pq.normal,
            })
        } else {
            None
        };
    }
    let result = ray_cast(start, dir, 1.0, shape)?;
    Some(SegmentCastResult {
        hit: result.hit,
        t: result.toi,
        hit_point: result.hit_point,
        normal: result.normal,
    })
}
/// Cast a ray against multiple shapes and return all hits, sorted by TOI.
pub fn ray_cast_batch(
    origin: [f64; 3],
    dir: [f64; 3],
    max_toi: f64,
    shapes: &[ShapeKind],
) -> Vec<(usize, RayCastResult)> {
    let mut hits: Vec<(usize, RayCastResult)> = shapes
        .iter()
        .enumerate()
        .filter_map(|(i, s)| ray_cast(origin, dir, max_toi, s).map(|r| (i, r)))
        .collect();
    hits.sort_by(|(_, a), (_, b)| {
        a.toi
            .partial_cmp(&b.toi)
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    hits
}
/// Query a point against multiple shapes and return all shapes that contain it.
pub fn point_query_batch(point: [f64; 3], shapes: &[ShapeKind]) -> Vec<(usize, PointQueryResult)> {
    shapes
        .iter()
        .enumerate()
        .filter_map(|(i, s)| {
            let r = point_query(point, s);
            if r.is_inside { Some((i, r)) } else { None }
        })
        .collect()
}
#[cfg(test)]
mod extra_tests {
    use super::super::*;
    use crate::narrowphase::TriangleMesh;
    use crate::narrowphase::types::{
        CompoundShape, ContactFeature, ContactFilter, FeatureContact, NarrowPhaseContact,
    };
    #[test]
    fn feature_contact_sphere_sphere_is_vertex_face() {
        let a = ShapeKind::Sphere {
            center: [0.0; 3],
            radius: 1.0,
        };
        let b = ShapeKind::Sphere {
            center: [1.5, 0.0, 0.0],
            radius: 1.0,
        };
        let fc = shape_shape_feature_contact(&a, &b).unwrap();
        assert!(matches!(fc.feature, ContactFeature::VertexFace { .. }));
    }
    #[test]
    fn feature_contact_capsule_capsule_is_edge_edge() {
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
        let fc = shape_shape_feature_contact(&a, &b).unwrap();
        assert!(matches!(fc.feature, ContactFeature::EdgeEdge { .. }));
    }
    #[test]
    fn feature_contact_no_overlap_returns_none() {
        let a = ShapeKind::Sphere {
            center: [0.0; 3],
            radius: 0.5,
        };
        let b = ShapeKind::Sphere {
            center: [10.0, 0.0, 0.0],
            radius: 0.5,
        };
        assert!(shape_shape_feature_contact(&a, &b).is_none());
    }
    #[test]
    fn feature_contact_from_plain() {
        let c = NarrowPhaseContact {
            normal: [1.0, 0.0, 0.0],
            depth: 0.1,
            point_a: [0.0; 3],
            point_b: [0.0; 3],
        };
        let fc = FeatureContact::from_plain(c);
        assert_eq!(fc.feature, ContactFeature::Unknown);
    }
    #[test]
    fn triangle_mesh_tri_count() {
        let mesh = TriangleMesh::new(vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [1.0, 0.0, 1.0],
            [0.0, 1.0, 1.0],
        ]);
        assert_eq!(mesh.tri_count(), 2);
    }
    #[test]
    fn triangle_mesh_get_triangle() {
        let mesh = TriangleMesh::new(vec![[1.0, 0.0, 0.0], [2.0, 0.0, 0.0], [0.0, 1.0, 0.0]]);
        let tri = mesh.triangle(0);
        assert_eq!(tri[0], [1.0, 0.0, 0.0]);
    }
    #[test]
    fn triangle_mesh_face_normal_z_axis() {
        let mesh = TriangleMesh::new(vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]);
        let n = mesh.face_normal(0);
        assert!(n[2].abs() > 1e-9, "normal z should be non-zero: {n:?}");
    }
    #[test]
    fn convex_vs_mesh_sphere_hits_close_triangle() {
        let mesh = TriangleMesh::new(vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]);
        let sphere = ShapeKind::Sphere {
            center: [0.3, 0.3, 0.2],
            radius: 1.0,
        };
        let contacts = convex_vs_mesh(&sphere, &mesh, 0.0);
        assert!(
            !contacts.is_empty(),
            "sphere overlapping triangle mesh centroid should produce contacts"
        );
    }
    #[test]
    fn convex_vs_mesh_sphere_misses_far_mesh() {
        let mesh = TriangleMesh::new(vec![
            [10.0, 10.0, 0.0],
            [11.0, 10.0, 0.0],
            [10.0, 11.0, 0.0],
        ]);
        let sphere = ShapeKind::Sphere {
            center: [0.0, 0.0, 0.0],
            radius: 0.5,
        };
        let contacts = convex_vs_mesh(&sphere, &mesh, 0.0);
        assert!(contacts.is_empty(), "sphere far from mesh should miss");
    }
    #[test]
    fn compound_aabb_covers_all_children() {
        let c = CompoundShape::new(vec![
            ShapeKind::Sphere {
                center: [0.0, 0.0, 0.0],
                radius: 1.0,
            },
            ShapeKind::Sphere {
                center: [5.0, 0.0, 0.0],
                radius: 1.0,
            },
        ]);
        let (mn, mx) = c.aabb();
        assert!(mn[0] <= -1.0, "mn[0]={}", mn[0]);
        assert!(mx[0] >= 6.0, "mx[0]={}", mx[0]);
    }
    #[test]
    fn compound_vs_shape_hit() {
        let compound = CompoundShape::new(vec![
            ShapeKind::Sphere {
                center: [0.0, 0.0, 0.0],
                radius: 1.0,
            },
            ShapeKind::Sphere {
                center: [3.0, 0.0, 0.0],
                radius: 1.0,
            },
        ]);
        let other = ShapeKind::Sphere {
            center: [0.5, 0.0, 0.0],
            radius: 0.8,
        };
        let contacts = compound_vs_shape(&compound, &other, &ContactFilter::default());
        assert!(
            !contacts.is_empty(),
            "compound should collide with nearby sphere"
        );
    }
    #[test]
    fn compound_vs_shape_miss() {
        let compound = CompoundShape::new(vec![ShapeKind::Sphere {
            center: [0.0, 0.0, 0.0],
            radius: 0.5,
        }]);
        let other = ShapeKind::Sphere {
            center: [10.0, 0.0, 0.0],
            radius: 0.5,
        };
        let contacts = compound_vs_shape(&compound, &other, &ContactFilter::default());
        assert!(contacts.is_empty(), "compound far from other should miss");
    }
    #[test]
    fn compound_vs_compound_both_hit() {
        let a = CompoundShape::new(vec![ShapeKind::Sphere {
            center: [0.0, 0.0, 0.0],
            radius: 1.0,
        }]);
        let b = CompoundShape::new(vec![ShapeKind::Sphere {
            center: [1.0, 0.0, 0.0],
            radius: 1.0,
        }]);
        let contacts = compound_vs_compound(&a, &b, &ContactFilter::default());
        assert!(!contacts.is_empty());
    }
    #[test]
    fn compound_vs_compound_miss() {
        let a = CompoundShape::new(vec![ShapeKind::Sphere {
            center: [0.0, 0.0, 0.0],
            radius: 0.3,
        }]);
        let b = CompoundShape::new(vec![ShapeKind::Sphere {
            center: [10.0, 0.0, 0.0],
            radius: 0.3,
        }]);
        let contacts = compound_vs_compound(&a, &b, &ContactFilter::default());
        assert!(contacts.is_empty());
    }
    #[test]
    fn compound_vs_shape_results_sorted_by_depth() {
        let compound = CompoundShape::new(vec![
            ShapeKind::Sphere {
                center: [0.0, 0.0, 0.0],
                radius: 2.0,
            },
            ShapeKind::Sphere {
                center: [1.8, 0.0, 0.0],
                radius: 0.5,
            },
        ]);
        let other = ShapeKind::Sphere {
            center: [2.0, 0.0, 0.0],
            radius: 0.5,
        };
        let contacts = compound_vs_shape(&compound, &other, &ContactFilter::default());
        for w in contacts.windows(2) {
            assert!(
                w[0].depth <= w[1].depth,
                "contacts should be sorted by depth"
            );
        }
    }
    #[test]
    fn point_inside_sphere() {
        let s = ShapeKind::Sphere {
            center: [0.0; 3],
            radius: 2.0,
        };
        let r = point_query([0.5, 0.0, 0.0], &s);
        assert!(r.is_inside, "point inside sphere should be inside");
        assert!(
            r.signed_distance < 0.0,
            "signed dist inside: {}",
            r.signed_distance
        );
    }
    #[test]
    fn point_outside_sphere() {
        let s = ShapeKind::Sphere {
            center: [0.0; 3],
            radius: 1.0,
        };
        let r = point_query([3.0, 0.0, 0.0], &s);
        assert!(!r.is_inside);
        assert!(
            r.signed_distance > 0.0,
            "signed dist outside: {}",
            r.signed_distance
        );
    }
    #[test]
    fn point_on_sphere_surface() {
        let s = ShapeKind::Sphere {
            center: [0.0; 3],
            radius: 1.0,
        };
        let r = point_query([1.0, 0.0, 0.0], &s);
        assert!(
            r.signed_distance.abs() < 1e-10,
            "surface dist: {}",
            r.signed_distance
        );
    }
    #[test]
    fn point_inside_box() {
        let b = ShapeKind::Box {
            center: [0.0; 3],
            half_extents: [1.0; 3],
        };
        let r = point_query([0.0, 0.0, 0.0], &b);
        assert!(r.is_inside);
    }
    #[test]
    fn point_outside_box() {
        let b = ShapeKind::Box {
            center: [0.0; 3],
            half_extents: [1.0; 3],
        };
        let r = point_query([3.0, 0.0, 0.0], &b);
        assert!(!r.is_inside);
        assert!(r.signed_distance > 0.0);
    }
    #[test]
    fn point_inside_capsule() {
        let c = ShapeKind::Capsule {
            p0: [0.0, 0.0, 0.0],
            p1: [0.0, 2.0, 0.0],
            radius: 1.0,
        };
        let r = point_query([0.0, 1.0, 0.0], &c);
        assert!(r.is_inside, "centre of capsule should be inside");
    }
    #[test]
    fn point_outside_capsule() {
        let c = ShapeKind::Capsule {
            p0: [0.0, 0.0, 0.0],
            p1: [0.0, 2.0, 0.0],
            radius: 0.5,
        };
        let r = point_query([5.0, 1.0, 0.0], &c);
        assert!(!r.is_inside);
    }
    #[test]
    fn point_below_plane_is_inside() {
        let p = ShapeKind::Plane {
            normal: [0.0, 1.0, 0.0],
            offset: 0.0,
        };
        let r = point_query([0.0, -1.0, 0.0], &p);
        assert!(r.is_inside, "below plane should be inside");
    }
    #[test]
    fn point_above_plane_is_outside() {
        let p = ShapeKind::Plane {
            normal: [0.0, 1.0, 0.0],
            offset: 0.0,
        };
        let r = point_query([0.0, 1.0, 0.0], &p);
        assert!(!r.is_inside);
    }
    #[test]
    fn ray_hits_sphere() {
        let s = ShapeKind::Sphere {
            center: [5.0, 0.0, 0.0],
            radius: 1.0,
        };
        let r = ray_cast([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], 10.0, &s);
        assert!(r.is_some(), "ray should hit sphere");
        let res = r.unwrap();
        assert!(res.toi > 0.0 && res.toi < 10.0, "toi={}", res.toi);
    }
    #[test]
    fn ray_misses_sphere() {
        let s = ShapeKind::Sphere {
            center: [5.0, 0.0, 0.0],
            radius: 0.3,
        };
        let r = ray_cast([0.0, 0.0, 0.0], [0.0, 1.0, 0.0], 10.0, &s);
        assert!(r.is_none(), "perpendicular ray should miss sphere");
    }
    #[test]
    fn ray_hits_plane() {
        let p = ShapeKind::Plane {
            normal: [0.0, 1.0, 0.0],
            offset: 0.0,
        };
        let r = ray_cast([0.0, 5.0, 0.0], [0.0, -1.0, 0.0], 10.0, &p);
        assert!(r.is_some(), "downward ray should hit XZ plane");
        let res = r.unwrap();
        assert!((res.toi - 5.0).abs() < 1e-6, "toi={}", res.toi);
    }
    #[test]
    fn ray_misses_plane_going_away() {
        let p = ShapeKind::Plane {
            normal: [0.0, 1.0, 0.0],
            offset: 0.0,
        };
        let r = ray_cast([0.0, 1.0, 0.0], [0.0, 1.0, 0.0], 10.0, &p);
        assert!(r.is_none(), "ray going away from plane should miss");
    }
    #[test]
    fn ray_hits_box() {
        let b = ShapeKind::Box {
            center: [5.0, 0.0, 0.0],
            half_extents: [1.0; 3],
        };
        let r = ray_cast([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], 10.0, &b);
        assert!(r.is_some(), "ray should hit box");
    }
    #[test]
    fn ray_misses_box() {
        let b = ShapeKind::Box {
            center: [0.0, 0.0, 0.0],
            half_extents: [1.0; 3],
        };
        let r = ray_cast([5.0, 5.0, 5.0], [1.0, 0.0, 0.0], 10.0, &b);
        assert!(r.is_none(), "ray offset and parallel should miss box");
    }
    #[test]
    fn ray_hits_capsule() {
        let c = ShapeKind::Capsule {
            p0: [5.0, -1.0, 0.0],
            p1: [5.0, 1.0, 0.0],
            radius: 1.0,
        };
        let r = ray_cast([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], 10.0, &c);
        assert!(r.is_some(), "ray should hit capsule");
    }
    #[test]
    fn segment_cast_hits_sphere() {
        let s = ShapeKind::Sphere {
            center: [0.5, 0.0, 0.0],
            radius: 0.3,
        };
        let r = segment_cast([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], &s);
        assert!(r.is_some(), "segment should hit sphere");
    }
    #[test]
    fn segment_cast_misses_sphere() {
        let s = ShapeKind::Sphere {
            center: [5.0, 0.0, 0.0],
            radius: 0.3,
        };
        let r = segment_cast([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], &s);
        assert!(r.is_none(), "short segment should miss distant sphere");
    }
    #[test]
    fn segment_cast_degenerate_point_inside() {
        let s = ShapeKind::Sphere {
            center: [0.0; 3],
            radius: 2.0,
        };
        let r = segment_cast([0.0, 0.0, 0.0], [0.0, 0.0, 0.0], &s);
        assert!(r.is_some(), "point inside sphere should count as hit");
    }
    #[test]
    fn segment_cast_degenerate_point_outside() {
        let s = ShapeKind::Sphere {
            center: [10.0, 0.0, 0.0],
            radius: 0.5,
        };
        let r = segment_cast([0.0, 0.0, 0.0], [0.0, 0.0, 0.0], &s);
        assert!(r.is_none(), "degenerate point outside sphere should miss");
    }
    #[test]
    fn ray_cast_batch_sorted_by_toi() {
        let shapes = vec![
            ShapeKind::Sphere {
                center: [8.0, 0.0, 0.0],
                radius: 1.0,
            },
            ShapeKind::Sphere {
                center: [3.0, 0.0, 0.0],
                radius: 1.0,
            },
            ShapeKind::Sphere {
                center: [5.0, 0.0, 0.0],
                radius: 1.0,
            },
        ];
        let hits = ray_cast_batch([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], 20.0, &shapes);
        assert_eq!(hits.len(), 3, "should hit all three spheres");
        for w in hits.windows(2) {
            assert!(w[0].1.toi <= w[1].1.toi, "hits should be sorted by toi");
        }
    }
    #[test]
    fn ray_cast_batch_empty_shapes() {
        let hits = ray_cast_batch([0.0; 3], [1.0, 0.0, 0.0], 10.0, &[]);
        assert!(hits.is_empty());
    }
    #[test]
    fn point_query_batch_finds_containing_shapes() {
        let shapes = vec![
            ShapeKind::Sphere {
                center: [0.0; 3],
                radius: 2.0,
            },
            ShapeKind::Sphere {
                center: [0.0; 3],
                radius: 5.0,
            },
            ShapeKind::Sphere {
                center: [10.0, 0.0, 0.0],
                radius: 1.0,
            },
        ];
        let hits = point_query_batch([0.0, 0.0, 0.0], &shapes);
        assert_eq!(
            hits.len(),
            2,
            "point at origin should be inside the first two spheres"
        );
    }
    #[test]
    fn point_query_batch_no_hits() {
        let shapes = vec![ShapeKind::Sphere {
            center: [10.0, 0.0, 0.0],
            radius: 0.5,
        }];
        let hits = point_query_batch([0.0, 0.0, 0.0], &shapes);
        assert!(hits.is_empty());
    }
    #[test]
    fn contact_feature_eq() {
        assert_eq!(ContactFeature::Unknown, ContactFeature::Unknown);
        assert_ne!(
            ContactFeature::FaceFace {
                face_a: 0,
                face_b: 1
            },
            ContactFeature::FaceFace {
                face_a: 0,
                face_b: 2
            },
        );
    }
}
