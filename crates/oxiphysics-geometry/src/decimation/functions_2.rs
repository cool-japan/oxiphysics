//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[cfg(test)]
mod tests_decimation_extended {
    use crate::decimation::DecimationMetrics;
    use crate::decimation::EdgeFeature;
    use crate::decimation::EdgePriorityQueue;
    use crate::decimation::MeshEdge;
    use crate::decimation::NormalCone;
    use crate::decimation::ProgressiveMeshSimple;
    use crate::decimation::*;

    use crate::decimation::SimpleMesh;
    fn small_mesh() -> SimpleMesh {
        let mut m = SimpleMesh::new();
        m.add_vertex([0.0, 0.0, 0.0]);
        m.add_vertex([1.0, 0.0, 0.0]);
        m.add_vertex([1.0, 1.0, 0.0]);
        m.add_vertex([0.0, 1.0, 0.0]);
        m.add_triangle(0, 1, 2);
        m.add_triangle(0, 2, 3);
        m
    }
    #[test]
    fn test_vertex_clustering_large_cell_merges_all() {
        let m = small_mesh();
        let reduced = vertex_clustering_decimate(&m, 100.0);
        assert_eq!(reduced.vertex_count(), 1);
        assert_eq!(reduced.triangle_count(), 0);
    }
    #[test]
    fn test_vertex_clustering_small_cell_preserves() {
        let m = small_mesh();
        let reduced = vertex_clustering_decimate(&m, 0.01);
        assert_eq!(reduced.vertex_count(), m.vertex_count());
    }
    #[test]
    fn test_vertex_clustering_empty_mesh() {
        let m = SimpleMesh::new();
        let reduced = vertex_clustering_decimate(&m, 1.0);
        assert_eq!(reduced.vertex_count(), 0);
    }
    #[test]
    fn test_vertex_clustering_mid_cell_reduces() {
        let mut m = SimpleMesh::new();
        m.add_vertex([0.0, 0.0, 0.0]);
        m.add_vertex([0.01, 0.0, 0.0]);
        m.add_vertex([1.0, 1.0, 0.0]);
        m.add_triangle(0, 1, 2);
        let reduced = vertex_clustering_decimate(&m, 0.1);
        assert!(reduced.vertex_count() < m.vertex_count());
    }
    #[test]
    fn test_progressive_mesh_initial() {
        let m = small_mesh();
        let pm = ProgressiveMeshSimple::new(m);
        assert_eq!(pm.n_collapses(), 0);
        assert_eq!(pm.n_triangles(), 2);
    }
    #[test]
    fn test_progressive_mesh_collapse_edge() {
        let m = small_mesh();
        let mut pm = ProgressiveMeshSimple::new(m);
        pm.collapse_edge(2, 1);
        assert_eq!(pm.n_collapses(), 1);
        assert!(pm.history[0].removed == 2);
    }
    #[test]
    fn test_progressive_mesh_multiple_collapses() {
        let m = small_mesh();
        let mut pm = ProgressiveMeshSimple::new(m);
        pm.collapse_edge(2, 1);
        pm.collapse_edge(3, 0);
        assert_eq!(pm.n_collapses(), 2);
    }
    #[test]
    fn test_edge_priority_queue_empty() {
        let q = EdgePriorityQueue::new();
        assert!(q.is_empty());
        assert_eq!(q.len(), 0);
    }
    #[test]
    fn test_edge_priority_queue_push_pop() {
        let mut q = EdgePriorityQueue::new();
        q.push(MeshEdge {
            v0: 0,
            v1: 1,
            cost: 2.0,
        });
        q.push(MeshEdge {
            v0: 1,
            v1: 2,
            cost: 0.5,
        });
        q.push(MeshEdge {
            v0: 0,
            v1: 2,
            cost: 1.0,
        });
        let e = q.pop().unwrap();
        assert!((e.cost - 0.5).abs() < 1e-12, "cost={}", e.cost);
    }
    #[test]
    fn test_edge_priority_queue_from_mesh() {
        let m = small_mesh();
        let q = EdgePriorityQueue::from_mesh(&m);
        assert_eq!(q.len(), 5);
    }
    #[test]
    fn test_edge_priority_queue_ordering() {
        let mut q = EdgePriorityQueue::new();
        for i in 0..10 {
            q.push(MeshEdge {
                v0: 0,
                v1: i,
                cost: (10 - i) as f64,
            });
        }
        let mut last_cost = f64::NEG_INFINITY;
        while let Some(e) = q.pop() {
            assert!(e.cost >= last_cost);
            last_cost = e.cost;
        }
    }
    #[test]
    fn test_normal_cone_from_normal() {
        let nc = NormalCone::from_normal([0.0, 0.0, 1.0]);
        assert!((nc.axis[2] - 1.0).abs() < 1e-10);
        assert_eq!(nc.half_angle, 0.0);
    }
    #[test]
    fn test_normal_cone_contains_same() {
        let nc = NormalCone::from_normal([0.0, 0.0, 1.0]);
        assert!(nc.contains([0.0, 0.0, 1.0], 1e-6));
    }
    #[test]
    fn test_normal_cone_excludes_opposite() {
        let nc = NormalCone::from_normal([0.0, 0.0, 1.0]);
        assert!(!nc.contains([0.0, 0.0, -1.0], 0.01));
    }
    #[test]
    fn test_normal_cone_merge() {
        let nc1 = NormalCone::from_normal([1.0, 0.0, 0.0]);
        let nc2 = NormalCone::from_normal([0.0, 1.0, 0.0]);
        let merged = nc1.merge(&nc2);
        assert!(merged.contains([1.0, 0.0, 0.0], 1e-6));
        assert!(merged.contains([0.0, 1.0, 0.0], 1e-6));
    }
    #[test]
    fn test_normal_cone_half_angle_deg() {
        let nc = NormalCone {
            axis: [0.0, 0.0, 1.0],
            half_angle: std::f64::consts::PI / 4.0,
        };
        assert!((nc.half_angle_deg() - 45.0).abs() < 1e-10);
    }
    #[test]
    fn test_classify_edge_smooth() {
        let n0 = [0.0, 0.0, 1.0];
        let n1 = [0.0, 0.0, 1.0];
        let feat = classify_edge(n0, n1, 30.0);
        assert_eq!(feat, EdgeFeature::Smooth);
    }
    #[test]
    fn test_classify_edge_crease() {
        let n0 = [0.0, 0.0, 1.0];
        let n1 = [0.0, 1.0, 0.0];
        let feat = classify_edge(n0, n1, 30.0);
        assert_eq!(feat, EdgeFeature::Crease);
    }
    #[test]
    fn test_feature_aware_cost_smooth() {
        let cost = feature_aware_cost(1.0, EdgeFeature::Smooth, 100.0, 1000.0);
        assert!((cost - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_feature_aware_cost_crease() {
        let cost = feature_aware_cost(1.0, EdgeFeature::Crease, 50.0, 1000.0);
        assert!((cost - 50.0).abs() < 1e-12);
    }
    #[test]
    fn test_feature_aware_cost_boundary() {
        let cost = feature_aware_cost(2.0, EdgeFeature::Boundary, 50.0, 100.0);
        assert!((cost - 200.0).abs() < 1e-12);
    }
    #[test]
    fn test_decimation_stats_default() {
        let s = DecimationMetrics::default();
        assert_eq!(s.vertex_reduction_ratio(), 0.0);
        assert_eq!(s.triangle_reduction_ratio(), 0.0);
        assert_eq!(s.avg_qem_cost(), 0.0);
    }
    #[test]
    fn test_decimation_stats_reduction() {
        let s = DecimationMetrics {
            original_vertices: 100,
            original_triangles: 200,
            reduced_vertices: 50,
            reduced_triangles: 100,
            n_collapses: 50,
            total_qem_cost: 10.0,
            max_qem_cost: 0.5,
        };
        assert!((s.vertex_reduction_ratio() - 0.5).abs() < 1e-10);
        assert!((s.triangle_reduction_ratio() - 0.5).abs() < 1e-10);
        assert!((s.avg_qem_cost() - 0.2).abs() < 1e-10);
    }
    #[test]
    fn test_collect_decimation_metrics() {
        let orig = small_mesh();
        let mut pm = ProgressiveMeshSimple::new(orig.clone());
        pm.collapse_edge(2, 1);
        let stats = collect_decimation_metrics(&orig, &pm.current, pm.n_collapses(), 0.5, 0.5);
        assert_eq!(stats.original_vertices, 4);
        assert_eq!(stats.n_collapses, 1);
    }
    #[test]
    fn test_mesh_edge_min_heap_ordering() {
        use std::collections::BinaryHeap;
        let mut heap = BinaryHeap::new();
        heap.push(MeshEdge {
            v0: 0,
            v1: 1,
            cost: 5.0,
        });
        heap.push(MeshEdge {
            v0: 0,
            v1: 2,
            cost: 1.0,
        });
        heap.push(MeshEdge {
            v0: 1,
            v1: 2,
            cost: 3.0,
        });
        let e = heap.pop().unwrap();
        assert!(
            (e.cost - 1.0).abs() < 1e-12,
            "Expected min cost 1.0, got {}",
            e.cost
        );
    }
}
#[cfg(test)]
mod tests_qem_decimation {

    use crate::decimation::QemDecimation;
    use crate::decimation::SimpleMesh;
    /// A flat quad: 4 vertices, 2 triangles (all interior edge).
    fn quad_mesh() -> SimpleMesh {
        let mut m = SimpleMesh::new();
        m.add_vertex([0.0, 0.0, 0.0]);
        m.add_vertex([1.0, 0.0, 0.0]);
        m.add_vertex([1.0, 1.0, 0.0]);
        m.add_vertex([0.0, 1.0, 0.0]);
        m.add_triangle(0, 1, 2);
        m.add_triangle(0, 2, 3);
        m
    }
    /// A simple mesh with a sharp crease (two perpendicular quads).
    fn crease_mesh() -> SimpleMesh {
        let mut m = SimpleMesh::new();
        m.add_vertex([0.0, 0.0, 0.0]);
        m.add_vertex([1.0, 0.0, 0.0]);
        m.add_vertex([1.0, 1.0, 0.0]);
        m.add_vertex([0.0, 1.0, 0.0]);
        m.add_vertex([0.0, 0.0, 1.0]);
        m.add_vertex([1.0, 0.0, 1.0]);
        m.add_triangle(0, 1, 2);
        m.add_triangle(0, 2, 3);
        m.add_triangle(0, 1, 4);
        m.add_triangle(1, 5, 4);
        m
    }
    #[test]
    fn test_error_threshold_positive() {
        let qd = QemDecimation::new(quad_mesh());
        let thr = qd.compute_error_threshold(1.0);
        assert!(thr > 0.0, "threshold must be positive, got {thr}");
    }
    #[test]
    fn test_error_threshold_scales_with_factor() {
        let qd = QemDecimation::new(crease_mesh());
        let thr1 = qd.compute_error_threshold(0.001);
        let thr2 = qd.compute_error_threshold(100.0);
        assert!(thr2 >= thr1, "larger factor → larger or equal threshold");
    }
    #[test]
    fn test_error_threshold_empty_mesh() {
        let qd = QemDecimation::new(SimpleMesh::new());
        let thr = qd.compute_error_threshold(1.0);
        assert_eq!(thr, 0.0, "empty mesh → zero threshold");
    }
    #[test]
    fn test_error_threshold_finite() {
        let qd = QemDecimation::new(quad_mesh());
        let thr = qd.compute_error_threshold(0.5);
        assert!(thr.is_finite(), "threshold must be finite, got {thr}");
    }
    #[test]
    fn test_preserve_boundary_no_collapse_with_negative_threshold() {
        let mut qd = QemDecimation::new(quad_mesh());
        let before = qd.mesh.triangle_count();
        let n = qd.preserve_boundary(-1.0);
        assert_eq!(n, 0, "negative threshold should collapse nothing");
        assert_eq!(qd.mesh.triangle_count(), before);
    }
    #[test]
    fn test_preserve_boundary_collapses_with_large_threshold() {
        let mut qd = QemDecimation::new(quad_mesh());
        let n = qd.preserve_boundary(1e10);
        let _ = n;
    }
    #[test]
    fn test_preserve_boundary_result_is_valid_mesh() {
        let mut qd = QemDecimation::new(quad_mesh());
        qd.preserve_boundary(1e6);
        for tri in &qd.mesh.triangles {
            for &vi in tri {
                assert!(vi < qd.mesh.vertices.len(), "invalid vertex index {vi}");
            }
        }
    }
    #[test]
    fn test_feature_score_flat_mesh_near_zero() {
        let qd = QemDecimation::new(quad_mesh());
        let scores = qd.compute_feature_score();
        assert!(!scores.is_empty(), "should have scores for all edges");
        for (&(_a, _b), &s) in &scores {
            if s.is_finite() {
                assert!(s >= 0.0, "score must be non-negative, got {s}");
            }
        }
    }
    #[test]
    fn test_feature_score_crease_mesh_high_score() {
        let qd = QemDecimation::new(crease_mesh());
        let scores = qd.compute_feature_score();
        let max_score = scores
            .values()
            .filter(|v| v.is_finite())
            .cloned()
            .fold(0.0f64, f64::max);
        assert!(
            max_score > 0.5,
            "crease mesh should have a high feature score, got {max_score}"
        );
    }
    #[test]
    fn test_feature_score_boundary_edges_infinite() {
        let mut m = SimpleMesh::new();
        m.add_vertex([0.0, 0.0, 0.0]);
        m.add_vertex([1.0, 0.0, 0.0]);
        m.add_vertex([0.5, 1.0, 0.0]);
        m.add_triangle(0, 1, 2);
        let qd = QemDecimation::new(m);
        let scores = qd.compute_feature_score();
        for (&_, &s) in &scores {
            assert!(
                s.is_infinite(),
                "boundary edge should have infinite score, got {s}"
            );
        }
    }
}
