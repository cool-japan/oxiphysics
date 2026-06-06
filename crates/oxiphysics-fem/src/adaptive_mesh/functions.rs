//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::{AdaptiveMesh2D, ConformingConstraint, Element2D};

pub(super) fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
pub(super) fn add3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}
pub(super) fn scale3(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}
pub(super) fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
pub(super) fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
pub(super) fn norm3(a: [f64; 3]) -> f64 {
    dot3(a, a).sqrt()
}
pub(super) fn midpoint3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    scale3(add3(a, b), 0.5)
}
pub(super) fn clamp_f64(v: f64, lo: f64, hi: f64) -> f64 {
    if v < lo {
        lo
    } else if v > hi {
        hi
    } else {
        v
    }
}
#[cfg(test)]
mod tests {

    use crate::adaptive_mesh::*;
    /// Build a regular tetrahedron with unit edge length.
    fn regular_tet() -> [[f64; 3]; 4] {
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.5, f64::sqrt(3.0) / 2.0, 0.0],
            [0.5, f64::sqrt(3.0) / 6.0, f64::sqrt(6.0) / 3.0],
        ]
    }
    fn flat_tet() -> [[f64; 3]; 4] {
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.5, 0.001, 0.0],
            [0.5, 0.0005, 0.0001],
        ]
    }
    fn simple_mesh() -> AdaptiveMesh {
        let nodes = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.5, 1.0, 0.0],
            [0.5, 0.333, 1.0],
            [0.5, 0.333, -1.0],
        ];
        let tets = vec![[0, 1, 2, 3], [0, 1, 2, 4]];
        AdaptiveMesh::new(nodes, tets)
    }
    #[test]
    fn quality_regular_tet_close_to_one() {
        let tet = regular_tet();
        let q = MeshQualityMetric::element_quality(tet);
        assert!(q > 0.9, "regular tet quality should be near 1, got {}", q);
    }
    #[test]
    fn quality_flat_tet_low() {
        let tet = flat_tet();
        let q = MeshQualityMetric::element_quality(tet);
        assert!(q < 0.3, "flat tet quality should be low, got {}", q);
    }
    #[test]
    fn aspect_ratio_regular_tet_near_one() {
        let tet = regular_tet();
        let ar = MeshQualityMetric::aspect_ratio(tet);
        assert!(
            ar < 2.0,
            "regular tet aspect ratio should be near 1, got {}",
            ar
        );
    }
    #[test]
    fn aspect_ratio_flat_tet_large() {
        let tet = flat_tet();
        let ar = MeshQualityMetric::aspect_ratio(tet);
        assert!(
            ar > 5.0,
            "flat tet aspect ratio should be large, got {}",
            ar
        );
    }
    #[test]
    fn jacobian_det_regular_tet_positive() {
        let tet = regular_tet();
        let j = MeshQualityMetric::jacobian_determinant(tet);
        assert!(j > 0.0, "regular tet jacobian should be positive");
    }
    #[test]
    fn jacobian_det_inverted_tet_negative() {
        let mut tet = regular_tet();
        tet.swap(1, 2);
        let j = MeshQualityMetric::jacobian_determinant(tet);
        assert!(j < 0.0, "inverted tet jacobian should be negative");
    }
    #[test]
    fn min_dihedral_angle_regular_tet_positive() {
        let tet = regular_tet();
        let a = MeshQualityMetric::min_dihedral_angle(tet);
        assert!(a > 0.0, "min dihedral angle should be positive");
    }
    #[test]
    fn volume_unit_tet_correct() {
        let tet = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        let v = MeshQualityMetric::volume(tet);
        assert!(
            (v - 1.0 / 6.0).abs() < 1e-12,
            "volume should be 1/6, got {}",
            v
        );
    }
    #[test]
    fn quality_degenerate_tet_is_zero() {
        let tet = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [3.0, 0.0, 0.0],
        ];
        let q = MeshQualityMetric::element_quality(tet);
        assert!(q < 1e-12, "coplanar tet quality should be 0");
    }
    #[test]
    fn refinement_criterion_error_based_debug() {
        let c = RefinementCriterion::ErrorBased { threshold: 0.1 };
        let s = format!("{:?}", c);
        assert!(s.contains("ErrorBased"));
    }
    #[test]
    fn refinement_criterion_uniform_debug() {
        let c = RefinementCriterion::Uniform;
        let s = format!("{:?}", c);
        assert!(s.contains("Uniform"));
    }
    #[test]
    fn queue_push_pop_order() {
        let mut q = RefinementQueue::new();
        q.push(0, 1.0);
        q.push(1, 3.0);
        q.push(2, 2.0);
        assert_eq!(q.pop(), Some(1));
        assert_eq!(q.pop(), Some(2));
        assert_eq!(q.pop(), Some(0));
        assert_eq!(q.pop(), None);
    }
    #[test]
    fn queue_len_correct() {
        let mut q = RefinementQueue::new();
        assert_eq!(q.len(), 0);
        q.push(0, 1.0);
        q.push(1, 2.0);
        assert_eq!(q.len(), 2);
    }
    #[test]
    fn queue_is_empty() {
        let mut q = RefinementQueue::new();
        assert!(q.is_empty());
        q.push(0, 1.0);
        assert!(!q.is_empty());
    }
    #[test]
    fn queue_mark_all_above_threshold() {
        let mut q = RefinementQueue::new();
        let errors = [0.05, 0.2, 0.01, 0.15];
        q.mark_all_above_threshold(&errors, 0.1);
        assert_eq!(q.len(), 2);
    }
    #[test]
    fn queue_default_is_empty() {
        let q = RefinementQueue::default();
        assert!(q.is_empty());
    }
    #[test]
    fn edge_split_adds_midpoint() {
        let mut nodes = vec![[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let mid_idx = EdgeSplit::split_edge(&mut nodes, 0, 1);
        assert_eq!(mid_idx, 2);
        assert!((nodes[2][0] - 1.0).abs() < 1e-12);
    }
    #[test]
    fn split_longest_edge_finds_longest() {
        let nodes = vec![
            [0.0, 0.0, 0.0],
            [100.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
        ];
        let tet = [0, 1, 2, 3];
        let (a, b) = EdgeSplit::split_longest_edge(&tet, &nodes);
        assert!(
            a == 1 || b == 1,
            "longest edge should involve node 1, got ({}, {})",
            a,
            b
        );
    }
    #[test]
    fn edge_collapse_midpoint_correct() {
        let mp = EdgeCollapse::midpoint([0.0, 0.0, 0.0], [2.0, 4.0, 6.0]);
        assert!((mp[0] - 1.0).abs() < 1e-12);
        assert!((mp[1] - 2.0).abs() < 1e-12);
        assert!((mp[2] - 3.0).abs() < 1e-12);
    }
    #[test]
    fn edge_collapse_candidate_finds_shortest() {
        let nodes = vec![
            [0.0, 0.0, 0.0],
            [0.1, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 1.0],
        ];
        let tet = [0, 1, 2, 3];
        let (a, b) = EdgeCollapse::candidate(&tet, &nodes);
        assert!(
            (a == 0 && b == 1) || (a == 1 && b == 0),
            "shortest edge should be 0-1"
        );
    }
    #[test]
    fn edge_collapse_removes_degenerate_tets() {
        let mut nodes = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.5, 1.0, 0.0],
            [0.5, 0.5, 1.0],
        ];
        let mut tets = vec![[0, 1, 2, 3]];
        EdgeCollapse::collapse(&mut nodes, &mut tets, 0, 1);
        assert_eq!(tets.len(), 0);
    }
    #[test]
    fn optimize_node_moves_to_centroid() {
        let node = [0.0, 0.0, 0.0];
        let neighbors = [[2.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 2.0]];
        let opt = NodeSmoothing::optimize_node(node, &neighbors);
        assert!((opt[0] - 2.0 / 3.0).abs() < 1e-12);
    }
    #[test]
    fn optimize_node_no_neighbors_returns_same() {
        let node = [1.0, 2.0, 3.0];
        let opt = NodeSmoothing::optimize_node(node, &[]);
        assert_eq!(opt, node);
    }
    #[test]
    fn laplacian_smooth_moves_interior_node() {
        let mut nodes = vec![
            [0.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [0.0, 2.0, 0.0],
            [1.0, 1.0, 0.0],
        ];
        let connectivity = vec![vec![], vec![], vec![], vec![0, 1, 2]];
        let boundary = vec![true, true, true, false];
        NodeSmoothing::laplacian_smooth(&mut nodes, &connectivity, 5, &boundary);
        assert!((nodes[3][0] - 2.0 / 3.0).abs() < 1e-12);
        assert!((nodes[3][1] - 2.0 / 3.0).abs() < 1e-12);
    }
    #[test]
    fn weighted_laplacian_respects_weights() {
        let mut nodes = vec![[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [0.0, 0.0, 0.0]];
        let connectivity = vec![vec![], vec![], vec![0, 1]];
        let weights = vec![vec![], vec![], vec![1.0, 3.0]];
        let boundary = vec![true, true, false];
        NodeSmoothing::weighted_laplacian(&mut nodes, &weights, &connectivity, 1, &boundary);
        assert!((nodes[2][0] - 1.5).abs() < 1e-12);
    }
    #[test]
    fn adaptive_mesh_new_creates_correct_levels() {
        let mesh = simple_mesh();
        assert_eq!(mesh.refinement_level.len(), mesh.tets.len());
        assert!(mesh.refinement_level.iter().all(|&l| l == 0));
    }
    #[test]
    fn refine_element_increases_tet_count() {
        let mut mesh = simple_mesh();
        let n_before = mesh.tets.len();
        mesh.refine_element(0);
        assert!(mesh.tets.len() > n_before);
    }
    #[test]
    fn refine_element_increases_node_count() {
        let mut mesh = simple_mesh();
        let n_before = mesh.nodes.len();
        mesh.refine_element(0);
        assert!(mesh.nodes.len() > n_before);
    }
    #[test]
    fn min_quality_returns_finite() {
        let mesh = simple_mesh();
        let q = mesh.min_quality();
        assert!(q.is_finite());
    }
    #[test]
    fn mean_quality_returns_finite() {
        let mesh = simple_mesh();
        let q = mesh.mean_quality();
        assert!(q.is_finite());
        assert!((0.0..=1.0).contains(&q));
    }
    #[test]
    fn quality_histogram_bins_count_all_tets() {
        let mesh = simple_mesh();
        let n_bins = 5;
        let hist = mesh.quality_histogram(n_bins);
        assert_eq!(hist.len(), n_bins);
        let total: usize = hist.iter().map(|(_, c)| c).sum();
        assert_eq!(total, mesh.tets.len());
    }
    #[test]
    fn refine_to_criterion_uniform_doubles_tets() {
        let mut mesh = simple_mesh();
        let n_before = mesh.tets.len();
        let errors = vec![0.0; n_before];
        mesh.refine_to_criterion(&RefinementCriterion::Uniform, &errors);
        assert!(mesh.tets.len() >= n_before * 2);
    }
    #[test]
    fn refine_to_criterion_error_based_selects_correct() {
        let mut mesh = simple_mesh();
        let errors = vec![0.5, 0.01];
        let n_before = mesh.tets.len();
        mesh.refine_to_criterion(&RefinementCriterion::ErrorBased { threshold: 0.1 }, &errors);
        assert!(mesh.tets.len() > n_before);
    }
    #[test]
    fn smooth_runs_without_panic() {
        let mut mesh = simple_mesh();
        mesh.smooth(3);
        assert!(mesh.nodes.iter().all(|n| n.iter().all(|x| x.is_finite())));
    }
    #[test]
    fn octree_insert_and_query_sphere_finds_element() {
        let mut root = OctreeNode::new([0.0, 0.0, 0.0], 10.0);
        root.insert(0, [-1.0, -1.0, -1.0], [1.0, 1.0, 1.0], 4);
        let results = root.query_sphere([0.0, 0.0, 0.0], 2.0);
        assert!(!results.is_empty(), "should find inserted element");
        assert!(results.contains(&0));
    }
    #[test]
    fn octree_query_sphere_misses_far_element() {
        let mut root = OctreeNode::new([0.0, 0.0, 0.0], 10.0);
        root.insert(0, [8.0, 8.0, 8.0], [9.0, 9.0, 9.0], 4);
        let results = root.query_sphere([0.0, 0.0, 0.0], 1.0);
        assert!(
            !results.contains(&0),
            "far element should not appear in near sphere query"
        );
    }
    #[test]
    fn octree_query_box_finds_element() {
        let mut root = OctreeNode::new([0.0, 0.0, 0.0], 10.0);
        root.insert(5, [-0.5, -0.5, -0.5], [0.5, 0.5, 0.5], 4);
        let results = root.query_box([-1.0, -1.0, -1.0], [1.0, 1.0, 1.0]);
        assert!(results.contains(&5));
    }
    #[test]
    fn octree_query_box_misses_non_overlapping() {
        let mut root = OctreeNode::new([0.0, 0.0, 0.0], 10.0);
        root.insert(3, [5.0, 5.0, 5.0], [6.0, 6.0, 6.0], 4);
        let results = root.query_box([-1.0, -1.0, -1.0], [1.0, 1.0, 1.0]);
        assert!(!results.contains(&3));
    }
    #[test]
    fn octree_multiple_elements() {
        let mut root = OctreeNode::new([0.0, 0.0, 0.0], 20.0);
        for i in 0..8 {
            let x = (i as f64) * 2.0;
            root.insert(i, [x, 0.0, 0.0], [x + 1.0, 1.0, 1.0], 3);
        }
        let results = root.query_box([0.0, 0.0, 0.0], [15.0, 1.0, 1.0]);
        assert!(!results.is_empty());
    }
    #[test]
    fn red_refine_produces_8_tets() {
        let nodes = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        let tets = vec![[0, 1, 2, 3]];
        let mut lr = LocalRefinement::new(nodes, tets);
        let children = lr.red_refine(0);
        assert_eq!(lr.tets.len(), 8, "red refine should produce exactly 8 tets");
        for &c in &children {
            assert!(c < lr.tets.len(), "child index {} out of bounds", c);
        }
    }
    #[test]
    fn red_refine_inserts_6_midpoint_nodes() {
        let nodes = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        let tets = vec![[0, 1, 2, 3]];
        let mut lr = LocalRefinement::new(nodes, tets);
        lr.red_refine(0);
        assert_eq!(lr.nodes.len(), 10, "4 original + 6 midpoints = 10 nodes");
    }
    #[test]
    fn green_refine_produces_2_tets() {
        let nodes = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        let tets = vec![[0, 1, 2, 3]];
        let mut lr = LocalRefinement::new(nodes, tets);
        let children = lr.green_refine(0, 0);
        assert_eq!(lr.tets.len(), 2, "green refine should produce 2 tets");
        for &c in &children {
            assert!(c < lr.tets.len());
        }
    }
    #[test]
    fn green_refine_inserts_1_midpoint_node() {
        let nodes = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        let tets = vec![[0, 1, 2, 3]];
        let mut lr = LocalRefinement::new(nodes, tets);
        lr.green_refine(0, 0);
        assert_eq!(lr.nodes.len(), 5, "4 original + 1 midpoint = 5 nodes");
    }
    #[test]
    fn green_refine_all_edges_valid() {
        for edge_idx in 0..6 {
            let nodes = vec![
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [0.0, 0.0, 1.0],
            ];
            let tets = vec![[0, 1, 2, 3]];
            let mut lr = LocalRefinement::new(nodes, tets);
            lr.green_refine(0, edge_idx);
            assert_eq!(lr.tets.len(), 2, "edge {} should produce 2 tets", edge_idx);
        }
    }
}
/// Compute a residual-based error indicator for element `elem`.
///
/// Uses the L2 norm of the nodal-value variation over the element as a
/// proxy for the residual error.
///
/// # Arguments
/// * `values`    – solution values at each node (indexed globally)
/// * `positions` – 2-D node positions
/// * `elem`      – element whose error is estimated
pub fn error_indicator_residual(values: &[f64], _positions: &[[f64; 2]], elem: &Element2D) -> f64 {
    if elem.nodes.is_empty() {
        return 0.0;
    }
    let node_vals: Vec<f64> = elem
        .nodes
        .iter()
        .filter(|&&ni| ni < values.len())
        .map(|&ni| values[ni])
        .collect();
    if node_vals.is_empty() {
        return 0.0;
    }
    let mean = node_vals.iter().sum::<f64>() / node_vals.len() as f64;
    let variance = node_vals.iter().map(|v| (v - mean).powi(2)).sum::<f64>();
    variance.sqrt()
}
/// Dörfler (bulk) marking strategy.
///
/// Marks the smallest set of elements whose cumulative error exceeds
/// `theta` times the total error.
///
/// # Arguments
/// * `errors` – per-element error indicators
/// * `theta`  – fraction of total error to capture (0 < θ ≤ 1)
///
/// Returns the sorted indices of marked elements.
pub fn mark_elements_doerfler(errors: &[f64], theta: f64) -> Vec<usize> {
    if errors.is_empty() {
        return Vec::new();
    }
    let total: f64 = errors.iter().sum();
    if total <= 0.0 {
        return Vec::new();
    }
    let target = theta * total;
    let mut order: Vec<usize> = (0..errors.len()).collect();
    order.sort_by(|&a, &b| {
        errors[b]
            .partial_cmp(&errors[a])
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    let mut cumulative = 0.0;
    let mut marked = Vec::new();
    for &i in &order {
        if cumulative >= target {
            break;
        }
        marked.push(i);
        cumulative += errors[i];
    }
    marked.sort_unstable();
    marked
}
/// Refine all marked elements in `mesh`.
///
/// # Arguments
/// * `mesh`   – mutable adaptive mesh
/// * `marked` – indices of elements to refine
pub fn refine_mesh(mesh: &mut AdaptiveMesh2D, marked: &[usize]) {
    let mut to_refine: Vec<usize> = marked.to_vec();
    to_refine.sort_unstable();
    for &idx in &to_refine {
        mesh.refine_element(idx);
    }
}
/// Coarsen leaf elements whose error falls below `coarsen_threshold`.
///
/// # Arguments
/// * `mesh`               – mutable adaptive mesh
/// * `coarsen_threshold`  – elements with error below this are coarsened
/// * `errors`             – per-element error indicators (indexed by leaf order)
pub fn coarsen_mesh(mesh: &mut AdaptiveMesh2D, coarsen_threshold: f64, errors: &[f64]) {
    let leaves = mesh.leaf_elements();
    let coarsen_ids: Vec<usize> = leaves
        .iter()
        .zip(errors.iter())
        .filter(|&(_, &e)| e < coarsen_threshold)
        .map(|(&id, _)| id)
        .collect();
    for id in coarsen_ids {
        if id < mesh.elements.len() {
            mesh.elements[id].is_leaf = false;
        }
    }
}
/// Detect hanging (non-conforming) nodes in a 2-D adaptive mesh.
///
/// A hanging node is one that lies on the interior of an edge of a
/// coarser neighbour element.  This simplified implementation detects
/// nodes that are not shared by all their edge-adjacent elements.
///
/// Returns a list of [`ConformingConstraint`]s.
pub fn hanging_nodes(mesh: &AdaptiveMesh2D) -> Vec<ConformingConstraint> {
    let mut edge_to_elems: std::collections::HashMap<(usize, usize), Vec<usize>> =
        std::collections::HashMap::new();
    for (ei, elem) in mesh.elements.iter().enumerate() {
        if !elem.is_leaf {
            continue;
        }
        let n = elem.nodes.len();
        for k in 0..n {
            let a = elem.nodes[k];
            let b = elem.nodes[(k + 1) % n];
            let key = if a < b { (a, b) } else { (b, a) };
            edge_to_elems.entry(key).or_default().push(ei);
        }
    }
    let mut constraints = Vec::new();
    for edge in edge_to_elems.keys() {
        let mid_x = (mesh.nodes[edge.0][0] + mesh.nodes[edge.1][0]) / 2.0;
        let mid_y = (mesh.nodes[edge.0][1] + mesh.nodes[edge.1][1]) / 2.0;
        for (ni, &npos) in mesh.nodes.iter().enumerate() {
            if ni == edge.0 || ni == edge.1 {
                continue;
            }
            if (npos[0] - mid_x).abs() < 1e-14 && (npos[1] - mid_y).abs() < 1e-14 {
                constraints.push(ConformingConstraint {
                    master_node: edge.0,
                    slave_node: ni,
                    weight: 0.5,
                });
                constraints.push(ConformingConstraint {
                    master_node: edge.1,
                    slave_node: ni,
                    weight: 0.5,
                });
            }
        }
    }
    constraints
}
/// Compute the 2-D element quality (aspect ratio) for a quad/triangle element.
///
/// Returns a value in (0, 1] where 1 is optimal.
///
/// # Arguments
/// * `mesh`    – adaptive mesh
/// * `elem_id` – index of the element to evaluate
pub fn mesh_quality_2d(mesh: &AdaptiveMesh2D, elem_id: usize) -> f64 {
    if elem_id >= mesh.elements.len() {
        return 0.0;
    }
    let elem = &mesh.elements[elem_id];
    if elem.nodes.len() < 3 {
        return 0.0;
    }
    let pts: Vec<[f64; 2]> = elem
        .nodes
        .iter()
        .filter(|&&ni| ni < mesh.nodes.len())
        .map(|&ni| mesh.nodes[ni])
        .collect();
    if pts.len() < 3 {
        return 0.0;
    }
    let n = pts.len();
    let mut min_edge = f64::INFINITY;
    let mut max_edge = 0.0_f64;
    for k in 0..n {
        let a = pts[k];
        let b = pts[(k + 1) % n];
        let d = ((a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2)).sqrt();
        if d < min_edge {
            min_edge = d;
        }
        if d > max_edge {
            max_edge = d;
        }
    }
    if max_edge < 1e-30 {
        return 0.0;
    }
    (min_edge / max_edge).clamp(0.0, 1.0)
}
/// Bisect a 1-D element `elem_id` by inserting its midpoint.
///
/// The element `[nodes[elem_id[0\\]], nodes[elem_id[1\]]]` is replaced by two
/// child elements sharing the new midpoint node.
///
/// # Arguments
/// * `nodes`   – mutable list of 1-D node coordinates
/// * `elements`– mutable list of 1-D elements (each is `[left_node, right_node]`)
/// * `elem_id` – index of the element to bisect
pub fn bisection_refinement_1d(
    nodes: &mut Vec<f64>,
    elements: &mut Vec<[usize; 2]>,
    elem_id: usize,
) {
    if elem_id >= elements.len() {
        return;
    }
    let [a, b] = elements[elem_id];
    if a >= nodes.len() || b >= nodes.len() {
        return;
    }
    let mid = (nodes[a] + nodes[b]) / 2.0;
    let mid_idx = nodes.len();
    nodes.push(mid);
    elements[elem_id] = [a, mid_idx];
    elements.push([mid_idx, b]);
}
#[cfg(test)]
mod tests_2d {
    use super::*;

    #[test]
    fn test_new_unit_square_creates_elements() {
        let mesh = AdaptiveMesh2D::new_unit_square(2);
        assert_eq!(mesh.element_count(), 4);
    }
    #[test]
    fn test_new_unit_square_node_count() {
        let mesh = AdaptiveMesh2D::new_unit_square(2);
        assert_eq!(mesh.node_count(), 9);
    }
    #[test]
    fn test_leaf_elements_all_initially() {
        let mesh = AdaptiveMesh2D::new_unit_square(3);
        assert_eq!(mesh.leaf_elements().len(), mesh.element_count());
    }
    #[test]
    fn test_refine_element_increases_element_count() {
        let mut mesh = AdaptiveMesh2D::new_unit_square(1);
        let before = mesh.element_count();
        mesh.refine_element(0);
        assert!(mesh.element_count() > before);
    }
    #[test]
    fn test_refine_element_increases_node_count() {
        let mut mesh = AdaptiveMesh2D::new_unit_square(1);
        let before = mesh.node_count();
        mesh.refine_element(0);
        assert!(mesh.node_count() > before);
    }
    #[test]
    fn test_refine_marks_parent_non_leaf() {
        let mut mesh = AdaptiveMesh2D::new_unit_square(1);
        mesh.refine_element(0);
        assert!(!mesh.elements[0].is_leaf);
    }
    #[test]
    fn test_leaf_elements_after_refinement() {
        let mut mesh = AdaptiveMesh2D::new_unit_square(1);
        mesh.refine_element(0);
        let leaves = mesh.leaf_elements();
        assert_eq!(leaves.len(), 4);
    }
    #[test]
    fn test_doerfler_marking_theta_1_marks_all() {
        let errors = vec![1.0, 2.0, 3.0, 4.0];
        let marked = mark_elements_doerfler(&errors, 1.0);
        assert_eq!(marked.len(), errors.len());
    }
    #[test]
    fn test_doerfler_marking_theta_0_marks_none() {
        let errors = vec![1.0, 2.0, 3.0];
        let marked = mark_elements_doerfler(&errors, 0.0);
        assert!(marked.is_empty());
    }
    #[test]
    fn test_doerfler_marking_sorted_output() {
        let errors = vec![1.0, 4.0, 2.0, 3.0];
        let marked = mark_elements_doerfler(&errors, 0.4);
        assert!(marked.contains(&1));
    }
    #[test]
    fn test_error_indicator_non_negative() {
        let mesh = AdaptiveMesh2D::new_unit_square(2);
        let values = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0];
        for elem in &mesh.elements {
            let e = error_indicator_residual(&values, &mesh.nodes, elem);
            assert!(e >= 0.0);
        }
    }
    #[test]
    fn test_error_indicator_zero_uniform_values() {
        let mesh = AdaptiveMesh2D::new_unit_square(2);
        let values = vec![5.0; 9];
        for elem in &mesh.elements {
            let e = error_indicator_residual(&values, &mesh.nodes, elem);
            assert!(e < 1e-12);
        }
    }
    #[test]
    fn test_mesh_quality_unit_square_element() {
        let mesh = AdaptiveMesh2D::new_unit_square(1);
        let q = mesh_quality_2d(&mesh, 0);
        assert!(q > 0.0 && q <= 1.0 + 1e-12);
    }
    #[test]
    fn test_mesh_quality_out_of_bounds() {
        let mesh = AdaptiveMesh2D::new_unit_square(1);
        assert_eq!(mesh_quality_2d(&mesh, 999), 0.0);
    }
    #[test]
    fn test_bisection_1d_increases_node_count() {
        let mut nodes = vec![0.0, 1.0];
        let mut elements = vec![[0usize, 1usize]];
        bisection_refinement_1d(&mut nodes, &mut elements, 0);
        assert_eq!(nodes.len(), 3);
    }
    #[test]
    fn test_bisection_1d_inserts_midpoint() {
        let mut nodes = vec![0.0, 2.0];
        let mut elements = vec![[0usize, 1usize]];
        bisection_refinement_1d(&mut nodes, &mut elements, 0);
        assert!((nodes[2] - 1.0).abs() < 1e-14);
    }
    #[test]
    fn test_bisection_1d_increases_element_count() {
        let mut nodes = vec![0.0, 1.0];
        let mut elements = vec![[0usize, 1usize]];
        bisection_refinement_1d(&mut nodes, &mut elements, 0);
        assert_eq!(elements.len(), 2);
    }
    #[test]
    fn test_bisection_1d_children_cover_parent() {
        let mut nodes = vec![0.0, 4.0];
        let mut elements = vec![[0usize, 1usize]];
        bisection_refinement_1d(&mut nodes, &mut elements, 0);
        let left = elements[0];
        let right = elements[1];
        assert!((nodes[left[0]] - 0.0).abs() < 1e-14);
        assert!((nodes[left[1]] - nodes[right[0]]).abs() < 1e-14);
        assert!((nodes[right[1]] - 4.0).abs() < 1e-14);
    }
    #[test]
    fn test_bisection_1d_out_of_bounds() {
        let mut nodes = vec![0.0, 1.0];
        let mut elements = vec![[0usize, 1usize]];
        bisection_refinement_1d(&mut nodes, &mut elements, 99);
        assert_eq!(nodes.len(), 2);
    }
    #[test]
    fn test_refine_mesh_marks_elements() {
        let mut mesh = AdaptiveMesh2D::new_unit_square(2);
        let before = mesh.element_count();
        let marked = vec![0];
        refine_mesh(&mut mesh, &marked);
        assert!(mesh.element_count() > before);
    }
    #[test]
    fn test_conforming_constraints_initial_mesh_none() {
        let mesh = AdaptiveMesh2D::new_unit_square(1);
        let constraints = hanging_nodes(&mesh);
        assert!(constraints.is_empty());
    }
    #[test]
    fn test_coarsen_mesh_reduces_leaves() {
        let mut mesh = AdaptiveMesh2D::new_unit_square(2);
        let errors = vec![0.001, 0.002, 0.001, 0.001];
        let before_leaves = mesh.leaf_elements().len();
        coarsen_mesh(&mut mesh, 0.01, &errors);
        assert!(mesh.leaf_elements().len() <= before_leaves);
    }
    #[test]
    fn test_element2d_fields() {
        let elem = Element2D {
            id: 5,
            nodes: vec![0, 1, 2, 3],
            level: 2,
            is_leaf: true,
        };
        assert_eq!(elem.id, 5);
        assert_eq!(elem.level, 2);
        assert!(elem.is_leaf);
    }
    #[test]
    fn test_conforming_constraint_fields() {
        let c = ConformingConstraint {
            master_node: 0,
            slave_node: 4,
            weight: 0.5,
        };
        assert_eq!(c.master_node, 0);
        assert_eq!(c.slave_node, 4);
        assert!((c.weight - 0.5).abs() < 1e-14);
    }
}
