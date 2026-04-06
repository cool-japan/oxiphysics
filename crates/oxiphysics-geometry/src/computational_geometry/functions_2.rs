//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;

#[cfg(test)]
mod tests {
    use super::*;
    use crate::computational_geometry::HalfEdgeMesh;
    use crate::computational_geometry::Line2D;
    use crate::computational_geometry::LineArrangement;
    use crate::computational_geometry::Obstacle2D;
    use crate::computational_geometry::Point2;
    use crate::computational_geometry::Segment2D;
    use crate::computational_geometry::SlabPointLocator;
    use crate::computational_geometry::VisibilityGraph;
    #[test]
    fn test_point2_add_sub() {
        let a = Point2::new(1.0, 2.0);
        let b = Point2::new(3.0, 4.0);
        let s = a.add(b);
        assert!((s.x - 4.0).abs() < 1e-12);
        assert!((s.y - 6.0).abs() < 1e-12);
        let d = b.sub(a);
        assert!((d.x - 2.0).abs() < 1e-12);
    }
    #[test]
    fn test_point2_cross() {
        let a = Point2::new(1.0, 0.0);
        let b = Point2::new(0.0, 1.0);
        assert!((a.cross(b) - 1.0).abs() < 1e-12);
        assert!((b.cross(a) + 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_point2_dist() {
        let a = Point2::new(0.0, 0.0);
        let b = Point2::new(3.0, 4.0);
        assert!((a.dist(b) - 5.0).abs() < 1e-10);
    }
    #[test]
    fn test_cross2_ccw() {
        let a = Point2::new(0.0, 0.0);
        let b = Point2::new(1.0, 0.0);
        let c = Point2::new(0.5, 1.0);
        assert!(Point2::cross2(a, b, c) > 0.0);
    }
    #[test]
    fn test_cross3_orthogonal() {
        let x: Point3 = [1.0, 0.0, 0.0];
        let y: Point3 = [0.0, 1.0, 0.0];
        let z = cross3(x, y);
        assert!((z[0]).abs() < 1e-12);
        assert!((z[1]).abs() < 1e-12);
        assert!((z[2] - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_normalize3() {
        let v: Point3 = [3.0, 4.0, 0.0];
        let n = normalize3(v);
        assert!((mag3(n) - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_dot3() {
        let a: Point3 = [1.0, 0.0, 0.0];
        let b: Point3 = [0.0, 1.0, 0.0];
        assert!((dot3(a, b)).abs() < 1e-12);
        let c: Point3 = [2.0, 3.0, 4.0];
        assert!((dot3(c, c) - 29.0).abs() < 1e-12);
    }
    #[test]
    fn test_halfedge_add_triangle() {
        let mut mesh = HalfEdgeMesh::new();
        let v0 = mesh.add_vertex([0.0, 0.0, 0.0]);
        let v1 = mesh.add_vertex([1.0, 0.0, 0.0]);
        let v2 = mesh.add_vertex([0.0, 1.0, 0.0]);
        let fid = mesh.add_triangle(v0, v1, v2);
        assert_eq!(mesh.num_faces(), 1);
        let verts = mesh.face_vertices(fid);
        assert_eq!(verts.len(), 3);
    }
    #[test]
    fn test_halfedge_face_normal_up() {
        let mut mesh = HalfEdgeMesh::new();
        let v0 = mesh.add_vertex([0.0, 0.0, 0.0]);
        let v1 = mesh.add_vertex([1.0, 0.0, 0.0]);
        let v2 = mesh.add_vertex([0.0, 1.0, 0.0]);
        let fid = mesh.add_triangle(v0, v1, v2);
        let n = mesh.face_normal(fid);
        assert!(n[2] > 0.9);
    }
    #[test]
    fn test_halfedge_face_centroid() {
        let mut mesh = HalfEdgeMesh::new();
        let v0 = mesh.add_vertex([0.0, 0.0, 0.0]);
        let v1 = mesh.add_vertex([3.0, 0.0, 0.0]);
        let v2 = mesh.add_vertex([0.0, 3.0, 0.0]);
        let fid = mesh.add_triangle(v0, v1, v2);
        let c = mesh.face_centroid(fid);
        assert!((c[0] - 1.0).abs() < 1e-10);
        assert!((c[1] - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_halfedge_twin_links() {
        let mut mesh = HalfEdgeMesh::new();
        let v0 = mesh.add_vertex([0.0, 0.0, 0.0]);
        let v1 = mesh.add_vertex([1.0, 0.0, 0.0]);
        let v2 = mesh.add_vertex([0.0, 1.0, 0.0]);
        let v3 = mesh.add_vertex([1.0, 1.0, 0.0]);
        mesh.add_triangle(v0, v1, v2);
        mesh.add_triangle(v1, v3, v2);
        mesh.build_twin_links();
        let has_twin = mesh.half_edges.iter().any(|he| he.twin != usize::MAX);
        assert!(has_twin);
    }
    #[test]
    fn test_convex_hull_3d_tetrahedron() {
        let pts: Vec<Point3> = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        let hull = convex_hull_3d(&pts).expect("Hull should succeed");
        assert_eq!(hull.faces.len(), 4);
    }
    #[test]
    fn test_convex_hull_3d_cube_faces() {
        let pts: Vec<Point3> = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [1.0, 0.0, 1.0],
            [0.0, 1.0, 1.0],
            [1.0, 1.0, 1.0],
        ];
        let hull = convex_hull_3d(&pts).expect("Cube hull should succeed");
        assert!(hull.faces.len() >= 12);
    }
    #[test]
    fn test_convex_hull_3d_volume_tetrahedron() {
        let pts: Vec<Point3> = vec![
            [0.0, 0.0, 0.0],
            [6.0, 0.0, 0.0],
            [0.0, 6.0, 0.0],
            [0.0, 0.0, 6.0],
        ];
        let hull = convex_hull_3d(&pts).expect("Hull should succeed");
        let vol = hull.volume();
        assert!((vol - 36.0).abs() < 1.0, "Volume was {vol}");
    }
    #[test]
    fn test_convex_hull_3d_surface_area_tetrahedron() {
        let pts: Vec<Point3> = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        let hull = convex_hull_3d(&pts).expect("Hull should succeed");
        let sa = hull.surface_area();
        assert!(sa > 2.0 && sa < 3.0, "Surface area was {sa}");
    }
    #[test]
    fn test_sutherland_hodgman_unit_squares_overlap() {
        let sq: Vec<Point2> = vec![
            Point2::new(0.0, 0.0),
            Point2::new(1.0, 0.0),
            Point2::new(1.0, 1.0),
            Point2::new(0.0, 1.0),
        ];
        let result = sutherland_hodgman(&sq, &sq);
        assert!(result.len() >= 4);
        let area = polygon_area(&result).abs();
        assert!((area - 1.0).abs() < 1e-9);
    }
    #[test]
    fn test_sutherland_hodgman_half_overlap() {
        let subject = vec![
            Point2::new(0.0, 0.0),
            Point2::new(2.0, 0.0),
            Point2::new(2.0, 1.0),
            Point2::new(0.0, 1.0),
        ];
        let clip = vec![
            Point2::new(0.0, 0.0),
            Point2::new(1.0, 0.0),
            Point2::new(1.0, 1.0),
            Point2::new(0.0, 1.0),
        ];
        let result = sutherland_hodgman(&subject, &clip);
        let area = polygon_area(&result).abs();
        assert!((area - 1.0).abs() < 1e-9);
    }
    #[test]
    fn test_sutherland_hodgman_no_overlap() {
        let sq1 = vec![
            Point2::new(0.0, 0.0),
            Point2::new(1.0, 0.0),
            Point2::new(1.0, 1.0),
            Point2::new(0.0, 1.0),
        ];
        let sq2 = vec![
            Point2::new(2.0, 0.0),
            Point2::new(3.0, 0.0),
            Point2::new(3.0, 1.0),
            Point2::new(2.0, 1.0),
        ];
        let result = sutherland_hodgman(&sq1, &sq2);
        assert!(result.is_empty());
    }
    #[test]
    fn test_polygon_area_square() {
        let sq = vec![
            Point2::new(0.0, 0.0),
            Point2::new(2.0, 0.0),
            Point2::new(2.0, 2.0),
            Point2::new(0.0, 2.0),
        ];
        assert!((polygon_area(&sq) - 4.0).abs() < 1e-10);
    }
    #[test]
    fn test_minkowski_sum_two_unit_squares() {
        let sq = vec![
            Point2::new(0.0, 0.0),
            Point2::new(1.0, 0.0),
            Point2::new(1.0, 1.0),
            Point2::new(0.0, 1.0),
        ];
        let result = minkowski_sum(&sq, &sq);
        assert!(result.len() >= 4);
        let area = polygon_area(&result).abs();
        assert!(area > 3.5 && area < 4.5, "Area was {area}");
    }
    #[test]
    fn test_minkowski_sum_nonempty() {
        let tri = vec![
            Point2::new(0.0, 0.0),
            Point2::new(1.0, 0.0),
            Point2::new(0.0, 1.0),
        ];
        let sq = vec![
            Point2::new(0.0, 0.0),
            Point2::new(1.0, 0.0),
            Point2::new(1.0, 1.0),
            Point2::new(0.0, 1.0),
        ];
        let result = minkowski_sum(&tri, &sq);
        assert!(!result.is_empty());
    }
    #[test]
    fn test_line_arrangement_three_lines() {
        let lines = vec![
            Line2D::through(Point2::new(0.0, 0.0), Point2::new(1.0, 0.0)),
            Line2D::through(Point2::new(0.0, 1.0), Point2::new(1.0, 2.0)),
            Line2D::through(Point2::new(0.0, 3.0), Point2::new(3.0, 0.0)),
        ];
        let arr = LineArrangement::build(&lines);
        assert_eq!(arr.num_vertices(), 3);
    }
    #[test]
    fn test_line_arrangement_two_lines_one_vertex() {
        let lines = vec![
            Line2D::through(Point2::new(0.0, 0.0), Point2::new(1.0, 0.0)),
            Line2D::through(Point2::new(0.5, -1.0), Point2::new(0.5, 1.0)),
        ];
        let arr = LineArrangement::build(&lines);
        assert_eq!(arr.num_vertices(), 1);
        let v = &arr.vertices[0];
        assert!((v.point.x - 0.5).abs() < 1e-9);
        assert!((v.point.y).abs() < 1e-9);
    }
    #[test]
    fn test_line_intersect_2d_basic() {
        let l1 = Line2D {
            a: 1.0,
            b: 0.0,
            c: 2.0,
        };
        let l2 = Line2D {
            a: 0.0,
            b: 1.0,
            c: 3.0,
        };
        let pt = line_intersect_2d(&l1, &l2).expect("Should intersect");
        assert!((pt.x - 2.0).abs() < 1e-10);
        assert!((pt.y - 3.0).abs() < 1e-10);
    }
    #[test]
    fn test_line_intersect_parallel() {
        let l1 = Line2D {
            a: 1.0,
            b: 0.0,
            c: 1.0,
        };
        let l2 = Line2D {
            a: 1.0,
            b: 0.0,
            c: 2.0,
        };
        assert!(line_intersect_2d(&l1, &l2).is_none());
    }
    #[test]
    fn test_slab_locator_build() {
        let segs = vec![
            Segment2D::new(Point2::new(0.0, 1.0), Point2::new(4.0, 1.0)),
            Segment2D::new(Point2::new(0.0, 2.0), Point2::new(4.0, 2.0)),
        ];
        let locator = SlabPointLocator::build(&segs, &[None, None]);
        assert!(!locator.slab_xs.is_empty());
    }
    #[test]
    fn test_slab_locator_query_above_first() {
        let segs = vec![Segment2D::new(Point2::new(0.0, 1.0), Point2::new(4.0, 1.0))];
        let locator = SlabPointLocator::build(&segs, &[None]);
        let result = locator.locate(Point2::new(2.0, 0.5));
        assert_eq!(result, Some(0));
    }
    #[test]
    fn test_slab_locator_query_above_all() {
        let segs = vec![Segment2D::new(Point2::new(0.0, 1.0), Point2::new(4.0, 1.0))];
        let locator = SlabPointLocator::build(&segs, &[None]);
        let result = locator.locate(Point2::new(2.0, 5.0));
        assert!(result.is_none());
    }
    #[test]
    fn test_visibility_graph_no_obstacles() {
        let vg = VisibilityGraph::build(&[]);
        assert_eq!(vg.num_nodes(), 0);
    }
    #[test]
    fn test_visibility_graph_square_obstacle() {
        let obs = Obstacle2D::new(vec![
            Point2::new(0.0, 0.0),
            Point2::new(1.0, 0.0),
            Point2::new(1.0, 1.0),
            Point2::new(0.0, 1.0),
        ]);
        let vg = VisibilityGraph::build(&[obs]);
        assert_eq!(vg.num_nodes(), 4);
    }
    #[test]
    fn test_visibility_graph_adjacent_vertices_visible() {
        let obs = Obstacle2D::new(vec![
            Point2::new(0.0, 0.0),
            Point2::new(1.0, 0.0),
            Point2::new(1.0, 1.0),
            Point2::new(0.0, 1.0),
        ]);
        let vg = VisibilityGraph::build(&[obs]);
        let v0_sees_v1 = vg.edges[0].contains(&1);
        assert!(v0_sees_v1);
    }
    #[test]
    fn test_art_gallery_triangle() {
        let tri = vec![
            Point2::new(0.0, 0.0),
            Point2::new(4.0, 0.0),
            Point2::new(2.0, 4.0),
        ];
        let result = art_gallery_greedy(&tri, 3);
        assert!(!result.guards.is_empty());
        assert!(result.covered.iter().all(|&c| c));
    }
    #[test]
    fn test_art_gallery_square_one_guard() {
        let sq = vec![
            Point2::new(0.0, 0.0),
            Point2::new(4.0, 0.0),
            Point2::new(4.0, 4.0),
            Point2::new(0.0, 4.0),
        ];
        let result = art_gallery_greedy(&sq, 4);
        assert!(result.covered.iter().all(|&c| c));
    }
    #[test]
    fn test_check_full_coverage_with_guard() {
        let tri = vec![
            Point2::new(0.0, 0.0),
            Point2::new(1.0, 0.0),
            Point2::new(0.0, 1.0),
        ];
        assert!(check_full_coverage(&tri, &[0]));
    }
    #[test]
    fn test_delaunay_four_points() {
        let pts = vec![
            Point2::new(0.0, 0.0),
            Point2::new(1.0, 0.0),
            Point2::new(1.0, 1.0),
            Point2::new(0.0, 1.0),
        ];
        let tris = delaunay_2d(&pts);
        assert_eq!(tris.len(), 2);
    }
    #[test]
    fn test_delaunay_triangle_count() {
        let pts: Vec<Point2> = (0..6)
            .map(|i| {
                let angle = i as f64 * std::f64::consts::TAU / 6.0;
                Point2::new(angle.cos(), angle.sin())
            })
            .collect();
        let tris = delaunay_2d(&pts);
        assert!(!tris.is_empty());
    }
    #[test]
    fn test_circumcircle_2d_basic() {
        let a = Point2::new(0.0, 0.0);
        let b = Point2::new(1.0, 0.0);
        let c = Point2::new(0.5, 1.0);
        let result = circumcircle_2d(a, b, c);
        assert!(result.is_some());
        let (center, r2) = result.unwrap();
        let r2a = center.dist_sq(a);
        let r2b = center.dist_sq(b);
        assert!((r2a - r2b).abs() < 1e-9);
        assert!((r2a - r2).abs() < 1e-9);
    }
    #[test]
    fn test_in_circumcircle_2d() {
        let a = Point2::new(0.0, 0.0);
        let b = Point2::new(2.0, 0.0);
        let c = Point2::new(1.0, 2.0);
        let inside = Point2::new(1.0, 0.5);
        let outside = Point2::new(5.0, 5.0);
        assert!(in_circumcircle_2d(a, b, c, inside));
        assert!(!in_circumcircle_2d(a, b, c, outside));
    }
    #[test]
    fn test_voronoi_from_delaunay_four_points() {
        let pts = vec![
            Point2::new(0.0, 0.0),
            Point2::new(2.0, 0.0),
            Point2::new(2.0, 2.0),
            Point2::new(0.0, 2.0),
        ];
        let tris = delaunay_2d(&pts);
        let cells = voronoi_from_delaunay(&pts, &tris);
        assert_eq!(cells.len(), 4);
        for cell in &cells {
            assert!(!cell.circumcenters.is_empty());
        }
    }
    #[test]
    fn test_ccw_positive() {
        let a = Point2::new(0.0, 0.0);
        let b = Point2::new(1.0, 0.0);
        let c = Point2::new(0.5, 1.0);
        assert!(ccw(a, b, c));
        assert!(!ccw(a, c, b));
    }
    #[test]
    fn test_collinear() {
        let a = Point2::new(0.0, 0.0);
        let b = Point2::new(1.0, 1.0);
        let c = Point2::new(2.0, 2.0);
        assert!(collinear(a, b, c));
        let d = Point2::new(1.0, 1.5);
        assert!(!collinear(a, b, d));
    }
    #[test]
    fn test_point_in_triangle() {
        let a = Point2::new(0.0, 0.0);
        let b = Point2::new(3.0, 0.0);
        let c = Point2::new(0.0, 3.0);
        assert!(point_in_triangle(Point2::new(0.5, 0.5), a, b, c));
        assert!(!point_in_triangle(Point2::new(3.0, 3.0), a, b, c));
    }
    #[test]
    fn test_point_in_convex_polygon() {
        let sq = vec![
            Point2::new(0.0, 0.0),
            Point2::new(4.0, 0.0),
            Point2::new(4.0, 4.0),
            Point2::new(0.0, 4.0),
        ];
        assert!(point_in_convex_polygon(Point2::new(2.0, 2.0), &sq));
        assert!(!point_in_convex_polygon(Point2::new(5.0, 2.0), &sq));
    }
    #[test]
    fn test_convex_hull_2d_square() {
        let pts = vec![
            Point2::new(0.0, 0.0),
            Point2::new(1.0, 0.0),
            Point2::new(1.0, 1.0),
            Point2::new(0.0, 1.0),
            Point2::new(0.5, 0.5),
        ];
        let hull = convex_hull_2d(&pts);
        assert_eq!(hull.len(), 4);
    }
    #[test]
    fn test_convex_hull_2d_triangle() {
        let pts = vec![
            Point2::new(0.0, 0.0),
            Point2::new(2.0, 0.0),
            Point2::new(1.0, 2.0),
            Point2::new(1.0, 0.5),
        ];
        let hull = convex_hull_2d(&pts);
        assert_eq!(hull.len(), 3);
    }
    #[test]
    fn test_segments_properly_intersect_yes() {
        let p = Point2::new(0.0, 0.0);
        let q = Point2::new(2.0, 2.0);
        let a = Point2::new(0.0, 2.0);
        let b = Point2::new(2.0, 0.0);
        assert!(segments_properly_intersect(p, q, a, b));
    }
    #[test]
    fn test_segments_properly_intersect_no() {
        let p = Point2::new(0.0, 0.0);
        let q = Point2::new(1.0, 0.0);
        let a = Point2::new(2.0, 0.0);
        let b = Point2::new(3.0, 0.0);
        assert!(!segments_properly_intersect(p, q, a, b));
    }
    #[test]
    fn test_segment2d_y_at() {
        let s = Segment2D::new(Point2::new(0.0, 0.0), Point2::new(4.0, 4.0));
        assert!((s.y_at(2.0).unwrap() - 2.0).abs() < 1e-10);
        assert!(s.y_at(5.0).is_none());
    }
    #[test]
    fn test_line2d_through_signed_dist() {
        let l = Line2D::through(Point2::new(0.0, 2.0), Point2::new(1.0, 2.0));
        let above = Point2::new(0.5, 5.0);
        let below = Point2::new(0.5, 0.0);
        assert!(l.signed_dist(above) * l.signed_dist(below) < 0.0);
    }
    #[test]
    fn test_convex_hull_3d_insufficient_points() {
        let pts: Vec<Point3> = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        assert!(convex_hull_3d(&pts).is_none());
    }
    #[test]
    fn test_obstacle_edges_count() {
        let obs = Obstacle2D::new(vec![
            Point2::new(0.0, 0.0),
            Point2::new(1.0, 0.0),
            Point2::new(1.0, 1.0),
            Point2::new(0.0, 1.0),
        ]);
        assert_eq!(obs.edges().len(), 4);
    }
    #[test]
    fn test_euler_face_count() {
        let lines = vec![
            Line2D::through(Point2::new(0.0, 0.0), Point2::new(1.0, 0.5)),
            Line2D::through(Point2::new(0.0, 1.0), Point2::new(1.0, 0.0)),
            Line2D::through(Point2::new(0.0, 0.5), Point2::new(1.0, 1.0)),
        ];
        let arr = LineArrangement::build(&lines);
        assert_eq!(arr.euler_face_count(), 7);
    }
    #[test]
    fn test_halfedge_vertex_faces() {
        let mut mesh = HalfEdgeMesh::new();
        let v0 = mesh.add_vertex([0.0, 0.0, 0.0]);
        let v1 = mesh.add_vertex([1.0, 0.0, 0.0]);
        let v2 = mesh.add_vertex([0.0, 1.0, 0.0]);
        mesh.add_triangle(v0, v1, v2);
        mesh.build_twin_links();
        let faces = mesh.vertex_faces(v0);
        assert!(!faces.is_empty());
    }
    #[test]
    fn test_convex_hull_2d_circle_approx() {
        let n = 8usize;
        let pts: Vec<Point2> = (0..n)
            .map(|i| {
                let t = i as f64 * std::f64::consts::TAU / n as f64;
                Point2::new(t.cos(), t.sin())
            })
            .collect();
        let hull = convex_hull_2d(&pts);
        assert_eq!(hull.len(), n);
    }
    #[test]
    fn test_sutherland_hodgman_triangle_clip() {
        let subject = vec![
            Point2::new(-1.0, -1.0),
            Point2::new(3.0, -1.0),
            Point2::new(1.0, 3.0),
        ];
        let clip = vec![
            Point2::new(0.0, 0.0),
            Point2::new(1.0, 0.0),
            Point2::new(1.0, 1.0),
            Point2::new(0.0, 1.0),
        ];
        let result = sutherland_hodgman(&subject, &clip);
        assert!(!result.is_empty());
        let area = polygon_area(&result).abs();
        assert!(area > 0.9 && area <= 1.0 + 1e-9, "area={area}");
    }
    #[test]
    fn test_voronoi_circumcenters_equidistant() {
        let pts = vec![
            Point2::new(0.0, 0.0),
            Point2::new(2.0, 0.0),
            Point2::new(1.0, 2.0),
        ];
        let tris = delaunay_2d(&pts);
        let cells = voronoi_from_delaunay(&pts, &tris);
        for t in &tris {
            if let Some((center, r2)) = circumcircle_2d(pts[t.a], pts[t.b], pts[t.c]) {
                let ra = center.dist_sq(pts[t.a]);
                let rb = center.dist_sq(pts[t.b]);
                assert!((ra - r2).abs() < 1e-8);
                assert!((rb - r2).abs() < 1e-8);
            }
        }
        assert_eq!(cells.len(), 3);
    }
}
