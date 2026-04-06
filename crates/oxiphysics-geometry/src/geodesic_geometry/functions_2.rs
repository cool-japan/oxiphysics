//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::needless_range_loop)]
#[allow(unused_imports)]
use super::functions::*;
use std::cmp::Ordering;
use std::collections::BinaryHeap;

use super::functions::{angle_at_vertex, dist3};
use super::types::{DistNode, GeodesicMesh};

/// Find all vertices within a geodesic radius of a source vertex.
///
/// Returns pairs of `(vertex_index, distance)` sorted by distance.
pub fn geodesic_ball(mesh: &GeodesicMesh, source: usize, radius: f64) -> Vec<(usize, f64)> {
    let adj = mesh.vertex_adjacency();
    let n = mesh.num_vertices();
    let mut dist = vec![f64::INFINITY; n];
    let mut heap = BinaryHeap::new();
    dist[source] = 0.0;
    heap.push(DistNode {
        vertex: source,
        dist: 0.0,
    });
    while let Some(DistNode { vertex: u, dist: d }) = heap.pop() {
        if d > dist[u] {
            continue;
        }
        if d > radius {
            continue;
        }
        for &v in &adj[u] {
            let w = dist3(mesh.vertices[u], mesh.vertices[v]);
            let new_dist = dist[u] + w;
            if new_dist < dist[v] && new_dist <= radius {
                dist[v] = new_dist;
                heap.push(DistNode {
                    vertex: v,
                    dist: new_dist,
                });
            }
        }
    }
    let mut result: Vec<(usize, f64)> = dist
        .iter()
        .enumerate()
        .filter(|&(_, d)| d.is_finite() && *d <= radius)
        .map(|(i, &d)| (i, d))
        .collect();
    result.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(Ordering::Equal));
    result
}
/// Compute discrete Gaussian curvature at a vertex using the angle defect.
///
/// For an interior vertex: `K = (2*pi - sum_of_angles) / voronoi_area`.
pub fn gaussian_curvature_at_vertex(mesh: &GeodesicMesh, vi: usize) -> f64 {
    let v2f = mesh.vertex_to_faces();
    let mut angle_sum = 0.0;
    for &fi in &v2f[vi] {
        let f = mesh.faces[fi];
        let local_idx = f.iter().position(|&v| v == vi).expect("element must exist");
        let a = mesh.vertices[f[local_idx]];
        let b = mesh.vertices[f[(local_idx + 1) % 3]];
        let c = mesh.vertices[f[(local_idx + 2) % 3]];
        angle_sum += angle_at_vertex(b, a, c);
    }
    let defect = 2.0 * std::f64::consts::PI - angle_sum;
    let area = mesh.voronoi_area(vi);
    defect / area
}
/// Compute discrete Gaussian curvature at all vertices.
pub fn gaussian_curvature_field(mesh: &GeodesicMesh) -> Vec<f64> {
    (0..mesh.num_vertices())
        .map(|vi| gaussian_curvature_at_vertex(mesh, vi))
        .collect()
}
/// Verify the Gauss-Bonnet theorem: sum of angle defects = 2*pi*chi
/// where chi is the Euler characteristic.
///
/// Returns the total angle defect (should be approximately 2*pi*chi for a
/// closed mesh).
pub fn total_angle_defect(mesh: &GeodesicMesh) -> f64 {
    let v2f = mesh.vertex_to_faces();
    let mut total = 0.0;
    for vi in 0..mesh.num_vertices() {
        let mut angle_sum = 0.0;
        for &fi in &v2f[vi] {
            let f = mesh.faces[fi];
            let local_idx = f.iter().position(|&v| v == vi).expect("element must exist");
            let a = mesh.vertices[f[local_idx]];
            let b = mesh.vertices[f[(local_idx + 1) % 3]];
            let c = mesh.vertices[f[(local_idx + 2) % 3]];
            angle_sum += angle_at_vertex(b, a, c);
        }
        total += 2.0 * std::f64::consts::PI - angle_sum;
    }
    total
}
/// Find the approximate geodesic midpoint between two vertices.
///
/// Traces the shortest path and returns the vertex closest to the midpoint
/// of the path.
pub fn geodesic_midpoint(mesh: &GeodesicMesh, v0: usize, v1: usize) -> usize {
    let dist = dijkstra_geodesic(mesh, &[v0]);
    let path = trace_geodesic_path(mesh, &dist, v0, v1);
    if path.is_empty() {
        return v0;
    }
    let total_len = geodesic_path_length(mesh, &path);
    let half_len = total_len * 0.5;
    let mut accum = 0.0;
    for i in 0..path.len() - 1 {
        let seg_len = dist3(mesh.vertices[path[i]], mesh.vertices[path[i + 1]]);
        if accum + seg_len >= half_len {
            let d_to_i = half_len - accum;
            let d_to_next = (accum + seg_len) - half_len;
            return if d_to_i <= d_to_next {
                path[i]
            } else {
                path[i + 1]
            };
        }
        accum += seg_len;
    }
    *path.last().expect("collection should not be empty")
}
/// Compute the all-pairs geodesic distance matrix for a set of vertices.
///
/// Returns an `n x n` matrix (flattened row-major) where `n = vertices.len()`.
pub fn geodesic_distance_matrix(mesh: &GeodesicMesh, vertices: &[usize]) -> Vec<f64> {
    let n = vertices.len();
    let mut matrix = vec![0.0; n * n];
    for (i, &vi) in vertices.iter().enumerate() {
        let dist = dijkstra_geodesic(mesh, &[vi]);
        for (j, &vj) in vertices.iter().enumerate() {
            matrix[i * n + j] = dist[vj];
        }
    }
    matrix
}
/// Build a simple tetrahedron mesh for testing.
pub(super) fn build_tetrahedron() -> GeodesicMesh {
    let vertices = vec![
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.5, 0.866, 0.0],
        [0.5, 0.289, 0.816],
    ];
    let faces = vec![[0, 1, 2], [0, 1, 3], [1, 2, 3], [0, 2, 3]];
    GeodesicMesh::new(vertices, faces)
}
/// Build a simple flat quad (2 triangles) for testing.
pub(super) fn build_flat_quad() -> GeodesicMesh {
    let vertices = vec![
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [1.0, 1.0, 0.0],
        [0.0, 1.0, 0.0],
    ];
    let faces = vec![[0, 1, 2], [0, 2, 3]];
    GeodesicMesh::new(vertices, faces)
}
/// Build a simple planar grid mesh.
pub(super) fn build_grid_mesh(nx: usize, ny: usize, spacing: f64) -> GeodesicMesh {
    let mut vertices = Vec::new();
    for iy in 0..ny {
        for ix in 0..nx {
            vertices.push([ix as f64 * spacing, iy as f64 * spacing, 0.0]);
        }
    }
    let mut faces = Vec::new();
    for iy in 0..(ny - 1) {
        for ix in 0..(nx - 1) {
            let v00 = iy * nx + ix;
            let v10 = iy * nx + ix + 1;
            let v01 = (iy + 1) * nx + ix;
            let v11 = (iy + 1) * nx + ix + 1;
            faces.push([v00, v10, v11]);
            faces.push([v00, v11, v01]);
        }
    }
    GeodesicMesh::new(vertices, faces)
}
/// Build an icosphere-like mesh for testing (octahedron subdivision).
pub(super) fn build_octahedron() -> GeodesicMesh {
    let vertices = vec![
        [1.0, 0.0, 0.0],
        [-1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, -1.0, 0.0],
        [0.0, 0.0, 1.0],
        [0.0, 0.0, -1.0],
    ];
    let faces = vec![
        [0, 2, 4],
        [2, 1, 4],
        [1, 3, 4],
        [3, 0, 4],
        [0, 5, 2],
        [2, 5, 1],
        [1, 5, 3],
        [3, 5, 0],
    ];
    GeodesicMesh::new(vertices, faces)
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::geodesic_geometry::HeatMethodParams;
    #[test]
    fn test_dijkstra_source_distance_zero() {
        let mesh = build_tetrahedron();
        let dist = dijkstra_geodesic(&mesh, &[0]);
        assert!((dist[0] - 0.0).abs() < 1e-10, "Source distance should be 0");
    }
    #[test]
    fn test_dijkstra_triangle_inequality() {
        let mesh = build_tetrahedron();
        let dist = dijkstra_geodesic(&mesh, &[0]);
        let d01 = dist[1];
        let d02 = dist[2];
        let d12 = dist3(mesh.vertices[1], mesh.vertices[2]);
        assert!(d02 <= d01 + d12 + 1e-10, "Triangle inequality violated");
    }
    #[test]
    fn test_dijkstra_symmetry() {
        let mesh = build_tetrahedron();
        let dist_from_0 = dijkstra_geodesic(&mesh, &[0]);
        let dist_from_1 = dijkstra_geodesic(&mesh, &[1]);
        assert!(
            (dist_from_0[1] - dist_from_1[0]).abs() < 1e-10,
            "Geodesic distance should be symmetric"
        );
    }
    #[test]
    fn test_dijkstra_flat_quad_diagonal() {
        let mesh = build_flat_quad();
        let dist = dijkstra_geodesic(&mesh, &[0]);
        let d02 = dist[2];
        assert!(
            d02 < 1.5,
            "Distance across diagonal should be close to sqrt(2), got {:.6}",
            d02
        );
    }
    #[test]
    fn test_dijkstra_multi_source() {
        let mesh = build_grid_mesh(5, 5, 1.0);
        let dist = dijkstra_geodesic(&mesh, &[0, 24]);
        for d in &dist {
            assert!(d.is_finite(), "All vertices should be reachable");
        }
        assert!(dist[12] < 5.0, "Center should be close to a source");
    }
    #[test]
    fn test_fmm_source_distance_zero() {
        let mesh = build_tetrahedron();
        let dist = fast_marching_geodesic(&mesh, &[0]);
        assert!(
            (dist[0] - 0.0).abs() < 1e-10,
            "FMM source distance should be 0"
        );
    }
    #[test]
    fn test_fmm_all_finite() {
        let mesh = build_tetrahedron();
        let dist = fast_marching_geodesic(&mesh, &[0]);
        for (i, &d) in dist.iter().enumerate() {
            assert!(
                d.is_finite(),
                "FMM vertex {} should have finite distance",
                i
            );
        }
    }
    #[test]
    fn test_fmm_no_greater_than_dijkstra() {
        let mesh = build_grid_mesh(4, 4, 1.0);
        let dijk = dijkstra_geodesic(&mesh, &[0]);
        let fmm = fast_marching_geodesic(&mesh, &[0]);
        for (i, (&dj, &df)) in dijk.iter().zip(fmm.iter()).enumerate() {
            assert!(
                df <= dj + 0.1,
                "FMM distance at {} ({:.6}) should not greatly exceed Dijkstra ({:.6})",
                i,
                df,
                dj
            );
        }
    }
    #[test]
    fn test_fmm_flat_grid_accuracy() {
        let mesh = build_grid_mesh(5, 5, 1.0);
        let dist = fast_marching_geodesic(&mesh, &[0]);
        let v_40 = 4;
        assert!(
            (dist[v_40] - 4.0).abs() < 0.5,
            "FMM distance on flat grid to (4,0) should be ~4.0, got {:.6}",
            dist[v_40]
        );
    }
    #[test]
    fn test_heat_method_source_zero() {
        let mesh = build_grid_mesh(5, 5, 1.0);
        let params = HeatMethodParams::default();
        let dist = heat_method_geodesic(&mesh, &[0], &params);
        let min_d = dist.iter().copied().fold(f64::INFINITY, f64::min);
        assert!(
            (dist[0] - min_d).abs() < 1e-6,
            "Source should have minimum distance"
        );
    }
    #[test]
    fn test_heat_method_monotone_on_line() {
        let mesh = build_grid_mesh(5, 1, 1.0);
        let params = HeatMethodParams::default();
        let dist = heat_method_geodesic(&mesh, &[0], &params);
        for i in 1..5 {
            assert!(
                dist[i] >= dist[i - 1] - 1e-6,
                "Heat method should give monotone distances on a line: d[{}]={:.6} < d[{}]={:.6}",
                i,
                dist[i],
                i - 1,
                dist[i - 1]
            );
        }
    }
    #[test]
    fn test_voronoi_single_site() {
        let mesh = build_tetrahedron();
        let (labels, cells) = geodesic_voronoi(&mesh, &[0]);
        assert_eq!(cells.len(), 1);
        for &l in &labels {
            assert_eq!(l, 0, "Single site => all vertices belong to it");
        }
        assert_eq!(cells[0].vertices.len(), 4);
    }
    #[test]
    fn test_voronoi_two_sites() {
        let mesh = build_grid_mesh(5, 5, 1.0);
        let (labels, cells) = geodesic_voronoi(&mesh, &[0, 24]);
        assert_eq!(cells.len(), 2);
        assert!(!cells[0].vertices.is_empty());
        assert!(!cells[1].vertices.is_empty());
        assert_eq!(labels[0], 0);
        assert_eq!(labels[24], 1);
    }
    #[test]
    fn test_voronoi_covers_all_vertices() {
        let mesh = build_grid_mesh(4, 4, 1.0);
        let (_labels, cells) = geodesic_voronoi(&mesh, &[0, 5, 10]);
        let total: usize = cells.iter().map(|c| c.vertices.len()).sum();
        assert_eq!(total, 16, "All 16 vertices should be assigned");
    }
    #[test]
    fn test_exp_map_zero_vector() {
        let mesh = build_tetrahedron();
        let result = discrete_exp_map(&mesh, 0, [0.0, 0.0], 100);
        assert_eq!(
            result.vertex, 0,
            "Zero tangent vector should stay at source"
        );
    }
    #[test]
    fn test_exp_map_returns_valid_face() {
        let mesh = build_tetrahedron();
        let result = discrete_exp_map(&mesh, 0, [0.1, 0.05], 100);
        assert!(
            result.face < mesh.num_faces(),
            "Result face should be valid"
        );
    }
    #[test]
    fn test_log_map_same_vertex() {
        let mesh = build_tetrahedron();
        let result = discrete_log_map(&mesh, 0, 0);
        assert!(
            (result.distance - 0.0).abs() < 1e-10,
            "Same vertex should have zero distance"
        );
        assert!((result.coords[0]).abs() < 1e-10);
        assert!((result.coords[1]).abs() < 1e-10);
    }
    #[test]
    fn test_log_map_distance_matches_dijkstra() {
        let mesh = build_tetrahedron();
        let result = discrete_log_map(&mesh, 0, 2);
        let dist = dijkstra_geodesic(&mesh, &[0]);
        assert!(
            (result.distance - dist[2]).abs() < 1e-10,
            "Log map distance should match Dijkstra"
        );
    }
    #[test]
    fn test_path_source_to_self() {
        let mesh = build_tetrahedron();
        let dist = dijkstra_geodesic(&mesh, &[0]);
        let path = trace_geodesic_path(&mesh, &dist, 0, 0);
        assert_eq!(path.len(), 1);
        assert_eq!(path[0], 0);
    }
    #[test]
    fn test_path_endpoints() {
        let mesh = build_grid_mesh(5, 5, 1.0);
        let dist = dijkstra_geodesic(&mesh, &[0]);
        let path = trace_geodesic_path(&mesh, &dist, 0, 24);
        assert_eq!(*path.first().unwrap(), 0, "Path should start at source");
        assert_eq!(*path.last().unwrap(), 24, "Path should end at target");
    }
    #[test]
    fn test_path_length_consistent() {
        let mesh = build_grid_mesh(4, 4, 1.0);
        let dist = dijkstra_geodesic(&mesh, &[0]);
        let path = trace_geodesic_path(&mesh, &dist, 0, 15);
        let length = geodesic_path_length(&mesh, &path);
        let dijk_dist = dist[15];
        assert!(
            length >= dijk_dist - 1e-6,
            "Path length {:.6} should be >= Dijkstra distance {:.6}",
            length,
            dijk_dist
        );
    }
    #[test]
    fn test_parallel_transport_preserves_length() {
        let mesh = build_grid_mesh(4, 4, 1.0);
        let path = vec![0, 1, 2, 3];
        let initial = [0.0, 1.0, 0.0];
        let transported = parallel_transport(&mesh, &path, initial);
        let initial_len = len3(transported[0]);
        for (i, &v) in transported.iter().enumerate() {
            let l = len3(v);
            assert!(
                (l - initial_len).abs() < 0.3,
                "Transport should approximately preserve length at step {}: got {:.6}, expected {:.6}",
                i,
                l,
                initial_len
            );
        }
    }
    #[test]
    fn test_parallel_transport_stays_tangent() {
        let mesh = build_octahedron();
        let path = vec![0, 2, 1];
        let initial = [0.0, 0.0, 1.0];
        let transported = parallel_transport(&mesh, &path, initial);
        for (i, &vi) in path.iter().enumerate() {
            let n = mesh.vertex_normal(vi);
            let dot = dot3(transported[i], n).abs();
            assert!(
                dot < 0.5,
                "Transported vector should be approximately tangent at step {}: normal dot = {:.6}",
                i,
                dot
            );
        }
    }
    #[test]
    fn test_fps_returns_requested_count() {
        let mesh = build_grid_mesh(5, 5, 1.0);
        let samples = farthest_point_sampling(&mesh, 0, 5);
        assert_eq!(samples.len(), 5, "FPS should return 5 samples");
    }
    #[test]
    fn test_fps_first_is_initial() {
        let mesh = build_grid_mesh(5, 5, 1.0);
        let samples = farthest_point_sampling(&mesh, 7, 3);
        assert_eq!(samples[0], 7, "First sample should be the initial vertex");
    }
    #[test]
    fn test_fps_unique_samples() {
        let mesh = build_grid_mesh(5, 5, 1.0);
        let samples = farthest_point_sampling(&mesh, 0, 10);
        let mut unique = samples.clone();
        unique.sort();
        unique.dedup();
        assert_eq!(unique.len(), samples.len(), "FPS samples should be unique");
    }
    #[test]
    fn test_fps_fmm_returns_requested_count() {
        let mesh = build_grid_mesh(4, 4, 1.0);
        let samples = farthest_point_sampling_fmm(&mesh, 0, 4);
        assert_eq!(samples.len(), 4, "FPS-FMM should return 4 samples");
    }
    #[test]
    fn test_centroid_single_point() {
        let mesh = build_tetrahedron();
        let c = geodesic_centroid(&mesh, &[2]);
        assert_eq!(c, 2, "Centroid of single point should be that point");
    }
    #[test]
    fn test_centroid_symmetric_points() {
        let mesh = build_grid_mesh(3, 3, 1.0);
        let corners = [0, 2, 6, 8];
        let c = geodesic_centroid(&mesh, &corners);
        assert!(c < mesh.num_vertices());
    }
    #[test]
    fn test_isolines_on_grid() {
        let mesh = build_grid_mesh(5, 5, 1.0);
        let dist = dijkstra_geodesic(&mesh, &[0]);
        let iso = geodesic_isolines(&mesh, &dist, 2.0);
        assert!(
            !iso.is_empty(),
            "Should find isoline points at distance 2.0 on a 5x5 grid"
        );
    }
    #[test]
    fn test_geodesic_ball_contains_source() {
        let mesh = build_grid_mesh(5, 5, 1.0);
        let ball = geodesic_ball(&mesh, 12, 1.5);
        assert!(
            ball.iter().any(|&(v, _)| v == 12),
            "Ball should contain the source vertex"
        );
    }
    #[test]
    fn test_geodesic_ball_respects_radius() {
        let mesh = build_grid_mesh(5, 5, 1.0);
        let ball = geodesic_ball(&mesh, 0, 2.0);
        for &(_v, d) in &ball {
            assert!(
                d <= 2.0 + 1e-10,
                "All vertices in ball should be within radius, got {:.6}",
                d
            );
        }
    }
    #[test]
    fn test_gauss_bonnet_tetrahedron() {
        let mesh = build_tetrahedron();
        let defect = total_angle_defect(&mesh);
        let expected = 4.0 * std::f64::consts::PI;
        assert!(
            (defect - expected).abs() < 1.0,
            "Total angle defect should be ~4*pi for tetrahedron, got {:.6}",
            defect
        );
    }
    #[test]
    fn test_gaussian_curvature_field_length() {
        let mesh = build_octahedron();
        let field = gaussian_curvature_field(&mesh);
        assert_eq!(field.len(), mesh.num_vertices());
    }
    #[test]
    fn test_midpoint_same_vertex() {
        let mesh = build_tetrahedron();
        let mid = geodesic_midpoint(&mesh, 0, 0);
        assert_eq!(mid, 0);
    }
    #[test]
    fn test_midpoint_on_path() {
        let mesh = build_grid_mesh(5, 5, 1.0);
        let mid = geodesic_midpoint(&mesh, 0, 4);
        assert_eq!(mid, 2, "Midpoint of 0..4 on bottom row should be 2");
    }
    #[test]
    fn test_distance_matrix_diagonal_zero() {
        let mesh = build_tetrahedron();
        let verts = vec![0, 1, 2, 3];
        let mat = geodesic_distance_matrix(&mesh, &verts);
        for i in 0..4 {
            assert!(
                (mat[i * 4 + i] - 0.0).abs() < 1e-10,
                "Diagonal should be zero"
            );
        }
    }
    #[test]
    fn test_distance_matrix_symmetric() {
        let mesh = build_tetrahedron();
        let verts = vec![0, 1, 2, 3];
        let mat = geodesic_distance_matrix(&mesh, &verts);
        for i in 0..4 {
            for j in 0..4 {
                assert!(
                    (mat[i * 4 + j] - mat[j * 4 + i]).abs() < 1e-10,
                    "Distance matrix should be symmetric"
                );
            }
        }
    }
    #[test]
    fn test_mean_edge_length_unit_grid() {
        let mesh = build_grid_mesh(3, 3, 1.0);
        let h = mean_edge_length(&mesh);
        assert!(
            h > 0.9 && h < 1.5,
            "Mean edge length should be reasonable, got {:.6}",
            h
        );
    }
    #[test]
    fn test_tetrahedron_vertex_count() {
        let mesh = build_tetrahedron();
        assert_eq!(mesh.num_vertices(), 4);
        assert_eq!(mesh.num_faces(), 4);
    }
    #[test]
    fn test_grid_mesh_dimensions() {
        let mesh = build_grid_mesh(5, 5, 1.0);
        assert_eq!(mesh.num_vertices(), 25);
        assert_eq!(mesh.num_faces(), 32);
    }
    #[test]
    fn test_mesh_total_area_flat_quad() {
        let mesh = build_flat_quad();
        let area = mesh.total_area();
        assert!(
            (area - 1.0).abs() < 1e-10,
            "Unit quad should have area 1.0, got {:.6}",
            area
        );
    }
}
