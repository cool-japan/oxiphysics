//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::*;
/// Hanging-node constraint: a slave node is linearly interpolated from
/// a master node with a given weight.
pub struct ConformingConstraint {
    /// Index of the master node.
    pub master_node: usize,
    /// Index of the slave (hanging) node.
    pub slave_node: usize,
    /// Interpolation weight (0 < w ≤ 1).
    pub weight: f64,
}
/// An adaptive tetrahedral mesh supporting refinement, coarsening, and smoothing.
pub struct AdaptiveMesh {
    /// Node positions.
    pub nodes: Vec<[f64; 3]>,
    /// Tetrahedra, each defined by 4 node indices.
    pub tets: Vec<[usize; 4]>,
    /// Per-element refinement level (0 = coarsest).
    pub refinement_level: Vec<u8>,
}
impl AdaptiveMesh {
    /// Create a new `AdaptiveMesh` from node positions and tetrahedral connectivity.
    pub fn new(nodes: Vec<[f64; 3]>, tets: Vec<[usize; 4]>) -> Self {
        let n_tets = tets.len();
        Self {
            nodes,
            tets,
            refinement_level: vec![0; n_tets],
        }
    }
    /// Refine element `idx` by bisecting its longest edge.
    ///
    /// Inserts a new midpoint node and replaces the tetrahedron with two children.
    pub fn refine_element(&mut self, idx: usize) {
        if idx >= self.tets.len() {
            return;
        }
        let tet = self.tets[idx];
        let (ei, ej) = EdgeSplit::split_longest_edge(&tet, &self.nodes);
        let mid_idx = EdgeSplit::split_edge(&mut self.nodes, ei, ej);
        let others: Vec<usize> = tet
            .iter()
            .copied()
            .filter(|&n| n != ei && n != ej)
            .collect();
        let old_level = self.refinement_level[idx];
        self.tets.remove(idx);
        self.refinement_level.remove(idx);
        if others.len() >= 2 {
            let o0 = others[0];
            let o1 = others[1];
            self.tets.push([ei, mid_idx, o0, o1]);
            self.refinement_level.push(old_level + 1);
            self.tets.push([ej, mid_idx, o0, o1]);
            self.refinement_level.push(old_level + 1);
        }
    }
    /// Coarsen element `idx` by collapsing its shortest edge (if quality is maintained).
    pub fn coarsen_element(&mut self, idx: usize) {
        if idx >= self.tets.len() {
            return;
        }
        let tet = self.tets[idx];
        let (ci, cj) = EdgeCollapse::candidate(&tet, &self.nodes);
        if EdgeCollapse::is_valid(&self.nodes, &self.tets, ci, cj) {
            EdgeCollapse::collapse(&mut self.nodes, &mut self.tets, ci, cj);
            if self.refinement_level.len() > self.tets.len() {
                self.refinement_level.truncate(self.tets.len());
            }
            while self.refinement_level.len() < self.tets.len() {
                self.refinement_level.push(0);
            }
        }
    }
    /// Apply `n_iter` iterations of Laplacian smoothing to interior nodes.
    pub fn smooth(&mut self, n_iter: usize) {
        let n_nodes = self.nodes.len();
        let mut connectivity: Vec<std::collections::HashSet<usize>> =
            vec![std::collections::HashSet::new(); n_nodes];
        for tet in &self.tets {
            for &i in tet {
                for &j in tet {
                    if i != j {
                        connectivity[i].insert(j);
                    }
                }
            }
        }
        let conn_vec: Vec<Vec<usize>> = connectivity
            .into_iter()
            .map(|s| s.into_iter().collect())
            .collect();
        let boundary = vec![false; n_nodes];
        NodeSmoothing::laplacian_smooth(&mut self.nodes, &conn_vec, n_iter, &boundary);
    }
    /// Compute a quality histogram.
    ///
    /// Returns a vector of `(bin_lower_bound, count)` pairs for `n_bins` equal-width bins
    /// spanning \[0, 1\].
    pub fn quality_histogram(&self, n_bins: usize) -> Vec<(f64, usize)> {
        let n = n_bins.max(1);
        let mut counts = vec![0usize; n];
        for tet in &self.tets {
            let pts = [
                self.nodes[tet[0]],
                self.nodes[tet[1]],
                self.nodes[tet[2]],
                self.nodes[tet[3]],
            ];
            let q = MeshQualityMetric::element_quality(pts);
            let bin = ((q * n as f64).floor() as usize).min(n - 1);
            counts[bin] += 1;
        }
        (0..n).map(|i| (i as f64 / n as f64, counts[i])).collect()
    }
    /// Return the minimum element quality across all tetrahedra.
    pub fn min_quality(&self) -> f64 {
        self.tets
            .iter()
            .map(|tet| {
                let pts = [
                    self.nodes[tet[0]],
                    self.nodes[tet[1]],
                    self.nodes[tet[2]],
                    self.nodes[tet[3]],
                ];
                MeshQualityMetric::element_quality(pts)
            })
            .fold(f64::INFINITY, f64::min)
    }
    /// Return the mean element quality across all tetrahedra.
    pub fn mean_quality(&self) -> f64 {
        if self.tets.is_empty() {
            return 0.0;
        }
        let sum: f64 = self
            .tets
            .iter()
            .map(|tet| {
                let pts = [
                    self.nodes[tet[0]],
                    self.nodes[tet[1]],
                    self.nodes[tet[2]],
                    self.nodes[tet[3]],
                ];
                MeshQualityMetric::element_quality(pts)
            })
            .sum();
        sum / self.tets.len() as f64
    }
    /// Refine elements according to a `RefinementCriterion` and per-element error indicators.
    pub fn refine_to_criterion(&mut self, criterion: &RefinementCriterion, errors: &[f64]) {
        let mut to_refine: Vec<usize> = Vec::new();
        match criterion {
            RefinementCriterion::ErrorBased { threshold } => {
                for (i, &err) in errors.iter().enumerate() {
                    if err > *threshold {
                        to_refine.push(i);
                    }
                }
            }
            RefinementCriterion::GradientBased { threshold } => {
                for (i, &grad) in errors.iter().enumerate() {
                    if grad > *threshold {
                        to_refine.push(i);
                    }
                }
            }
            RefinementCriterion::GeometricBased { target_size } => {
                for (i, tet) in self.tets.iter().enumerate() {
                    let (ei, ej) = EdgeSplit::split_longest_edge(tet, &self.nodes);
                    let len = norm3(sub3(self.nodes[ei], self.nodes[ej]));
                    if len > *target_size {
                        to_refine.push(i);
                    }
                }
            }
            RefinementCriterion::Uniform => {
                to_refine = (0..self.tets.len()).collect();
            }
        }
        to_refine.sort_unstable_by(|a, b| b.cmp(a));
        to_refine.dedup();
        for idx in to_refine {
            if idx < self.tets.len() {
                self.refine_element(idx);
            }
        }
    }
}
/// A 2-D quad or triangle element in an adaptive mesh.
pub struct Element2D {
    /// Unique element identifier.
    pub id: usize,
    /// Indices into the parent mesh node list.
    pub nodes: Vec<usize>,
    /// Refinement level (0 = coarsest).
    pub level: usize,
    /// `true` if this element has no children (active element).
    pub is_leaf: bool,
}
/// Mesh smoothing by repositioning interior nodes.
pub struct NodeSmoothing;
impl NodeSmoothing {
    /// Apply `iter` iterations of Laplacian smoothing.
    ///
    /// Each node is moved to the average of its neighbours.
    /// Boundary nodes (listed in `boundary`) are fixed.
    pub fn laplacian_smooth(
        nodes: &mut [[f64; 3]],
        connectivity: &[Vec<usize>],
        iter: usize,
        boundary: &[bool],
    ) {
        for _ in 0..iter {
            let old: Vec<[f64; 3]> = nodes.to_vec();
            for (idx, neighbors) in connectivity.iter().enumerate() {
                if boundary.get(idx).copied().unwrap_or(false) {
                    continue;
                }
                if neighbors.is_empty() {
                    continue;
                }
                let mut avg = [0.0f64; 3];
                for &nb in neighbors {
                    avg = add3(avg, old[nb]);
                }
                nodes[idx] = scale3(avg, 1.0 / neighbors.len() as f64);
            }
        }
    }
    /// Weighted Laplacian smoothing: each neighbour is weighted by `weights[i][k]`.
    pub fn weighted_laplacian(
        nodes: &mut [[f64; 3]],
        weights: &[Vec<f64>],
        connectivity: &[Vec<usize>],
        iter: usize,
        boundary: &[bool],
    ) {
        for _ in 0..iter {
            let old: Vec<[f64; 3]> = nodes.to_vec();
            for (idx, neighbors) in connectivity.iter().enumerate() {
                if boundary.get(idx).copied().unwrap_or(false) {
                    continue;
                }
                if neighbors.is_empty() {
                    continue;
                }
                let w_row = &weights[idx];
                let mut sum = [0.0f64; 3];
                let mut total_w = 0.0;
                for (k, &nb) in neighbors.iter().enumerate() {
                    let w = w_row.get(k).copied().unwrap_or(1.0);
                    sum = add3(sum, scale3(old[nb], w));
                    total_w += w;
                }
                if total_w > 1e-15 {
                    nodes[idx] = scale3(sum, 1.0 / total_w);
                }
            }
        }
    }
    /// Optimise a single node by moving it to the centroid of its neighbours.
    pub fn optimize_node(node: [f64; 3], neighbors: &[[f64; 3]]) -> [f64; 3] {
        if neighbors.is_empty() {
            return node;
        }
        let mut sum = [0.0f64; 3];
        for &nb in neighbors {
            sum = add3(sum, nb);
        }
        scale3(sum, 1.0 / neighbors.len() as f64)
    }
}
/// Operations for collapsing an edge (merging two nodes into one).
pub struct EdgeCollapse;
impl EdgeCollapse {
    /// Return a `(node_i, node_j)` collapse candidate pair from a tet.
    ///
    /// Picks the shortest edge as the collapse candidate.
    pub fn candidate(tet: &[usize; 4], nodes: &[[f64; 3]]) -> (usize, usize) {
        let pairs = [
            (tet[0], tet[1]),
            (tet[0], tet[2]),
            (tet[0], tet[3]),
            (tet[1], tet[2]),
            (tet[1], tet[3]),
            (tet[2], tet[3]),
        ];
        let mut best = pairs[0];
        let mut best_len = norm3(sub3(nodes[pairs[0].0], nodes[pairs[0].1]));
        for &(i, j) in pairs.iter().skip(1) {
            let l = norm3(sub3(nodes[i], nodes[j]));
            if l < best_len {
                best_len = l;
                best = (i, j);
            }
        }
        best
    }
    /// Check whether collapsing edge (i, j) preserves mesh validity.
    ///
    /// Returns `true` when the collapse is deemed safe (no element inversion expected).
    pub fn is_valid(nodes: &[[f64; 3]], tets: &[[usize; 4]], i: usize, j: usize) -> bool {
        let mp = Self::midpoint(nodes[i], nodes[j]);
        for tet in tets {
            if !tet.contains(&i) && !tet.contains(&j) {
                continue;
            }
            let mut new_tet = [[0.0f64; 3]; 4];
            for (k, &n) in tet.iter().enumerate() {
                if n == i || n == j {
                    new_tet[k] = mp;
                } else {
                    new_tet[k] = nodes[n];
                }
            }
            if MeshQualityMetric::jacobian_determinant(new_tet) < 0.0 {
                return false;
            }
        }
        true
    }
    /// Perform the edge collapse: node j is merged into node i (midpoint).
    ///
    /// Returns the updated tet list with degenerate tets (collapsed to a line/point) removed.
    pub fn collapse(nodes: &mut [[f64; 3]], tets: &mut Vec<[usize; 4]>, i: usize, j: usize) {
        let mp = Self::midpoint(nodes[i], nodes[j]);
        nodes[i] = mp;
        for tet in tets.iter_mut() {
            for n in tet.iter_mut() {
                if *n == j {
                    *n = i;
                }
            }
        }
        tets.retain(|t| {
            t[0] != t[1]
                && t[0] != t[2]
                && t[0] != t[3]
                && t[1] != t[2]
                && t[1] != t[3]
                && t[2] != t[3]
        });
    }
    /// Compute the midpoint between two nodes.
    pub fn midpoint(ni: [f64; 3], nj: [f64; 3]) -> [f64; 3] {
        midpoint3(ni, nj)
    }
}
/// Operations for splitting an edge (inserting a new midpoint node).
pub struct EdgeSplit;
impl EdgeSplit {
    /// Split edge (i, j) by inserting a midpoint node.
    ///
    /// Returns the index of the newly inserted node.
    pub fn split_edge(nodes: &mut Vec<[f64; 3]>, i: usize, j: usize) -> usize {
        let mp = midpoint3(nodes[i], nodes[j]);
        nodes.push(mp);
        nodes.len() - 1
    }
    /// Return the indices `(a, b)` of the longest edge in the tetrahedron.
    pub fn split_longest_edge(tet: &[usize; 4], nodes: &[[f64; 3]]) -> (usize, usize) {
        let pairs = [
            (tet[0], tet[1]),
            (tet[0], tet[2]),
            (tet[0], tet[3]),
            (tet[1], tet[2]),
            (tet[1], tet[3]),
            (tet[2], tet[3]),
        ];
        let mut best = pairs[0];
        let mut best_len = norm3(sub3(nodes[pairs[0].0], nodes[pairs[0].1]));
        for &(i, j) in pairs.iter().skip(1) {
            let l = norm3(sub3(nodes[i], nodes[j]));
            if l > best_len {
                best_len = l;
                best = (i, j);
            }
        }
        best
    }
}
/// 2-D adaptive mesh with quad/triangle elements.
///
/// Supports hierarchical refinement via quad-tree–style element splitting.
pub struct AdaptiveMesh2D {
    /// Node positions in the 2-D plane.
    pub nodes: Vec<[f64; 2]>,
    /// All elements (both active leaves and refined parents).
    pub elements: Vec<Element2D>,
    /// Maximum allowed refinement level.
    pub max_level: usize,
}
impl AdaptiveMesh2D {
    /// Create an `n_init × n_init` quad mesh on the unit square \[0,1\]².
    ///
    /// # Arguments
    /// * `n_init` – number of divisions along each axis (must be ≥ 1)
    pub fn new_unit_square(n_init: usize) -> Self {
        let n = n_init.max(1);
        let mut nodes = Vec::new();
        for j in 0..=n {
            for i in 0..=n {
                nodes.push([i as f64 / n as f64, j as f64 / n as f64]);
            }
        }
        let mut elements = Vec::new();
        let mut id = 0;
        let stride = n + 1;
        for j in 0..n {
            for i in 0..n {
                let bl = j * stride + i;
                let br = bl + 1;
                let tl = bl + stride;
                let tr = tl + 1;
                elements.push(Element2D {
                    id,
                    nodes: vec![bl, br, tr, tl],
                    level: 0,
                    is_leaf: true,
                });
                id += 1;
            }
        }
        Self {
            nodes,
            elements,
            max_level: 8,
        }
    }
    /// Total number of nodes.
    pub fn node_count(&self) -> usize {
        self.nodes.len()
    }
    /// Total number of elements (leaf and non-leaf).
    pub fn element_count(&self) -> usize {
        self.elements.len()
    }
    /// Indices of all leaf elements (active elements).
    pub fn leaf_elements(&self) -> Vec<usize> {
        self.elements
            .iter()
            .enumerate()
            .filter(|(_, e)| e.is_leaf)
            .map(|(i, _)| i)
            .collect()
    }
    /// Refine element `elem_id` by splitting it into four quad children.
    ///
    /// Marks the parent as non-leaf and inserts four children.
    /// Does nothing if the element is not a leaf or is already at `max_level`.
    pub fn refine_element(&mut self, elem_id: usize) {
        if elem_id >= self.elements.len() {
            return;
        }
        if !self.elements[elem_id].is_leaf {
            return;
        }
        if self.elements[elem_id].level >= self.max_level {
            return;
        }
        let parent_nodes = self.elements[elem_id].nodes.clone();
        if parent_nodes.len() < 4 {
            return;
        }
        let p: Vec<[f64; 2]> = parent_nodes.iter().map(|&ni| self.nodes[ni]).collect();
        let cx = (p[0][0] + p[1][0] + p[2][0] + p[3][0]) / 4.0;
        let cy = (p[0][1] + p[1][1] + p[2][1] + p[3][1]) / 4.0;
        let cn = self.nodes.len();
        self.nodes.push([cx, cy]);
        let mut en = [0usize; 4];
        for k in 0..4 {
            let kn = (k + 1) % 4;
            let mx = (p[k][0] + p[kn][0]) / 2.0;
            let my = (p[k][1] + p[kn][1]) / 2.0;
            en[k] = self.nodes.len();
            self.nodes.push([mx, my]);
        }
        let lev = self.elements[elem_id].level;
        let next_id = self.elements.len();
        let children = [
            vec![parent_nodes[0], en[0], cn, en[3]],
            vec![en[0], parent_nodes[1], en[1], cn],
            vec![cn, en[1], parent_nodes[2], en[2]],
            vec![en[3], cn, en[2], parent_nodes[3]],
        ];
        self.elements[elem_id].is_leaf = false;
        for (k, child_nodes) in children.iter().enumerate() {
            self.elements.push(Element2D {
                id: next_id + k,
                nodes: child_nodes.clone(),
                level: lev + 1,
                is_leaf: true,
            });
        }
    }
}
/// Collection of quality metrics for tetrahedral elements.
pub struct MeshQualityMetric;
impl MeshQualityMetric {
    /// Compute the aspect ratio of a tetrahedron (longest edge / shortest altitude).
    ///
    /// Returns a value ≥ 1; perfect regular tetrahedron → `√6 / (2√(2/3)) ≈ 1`.
    /// Larger values indicate more distorted elements.
    pub fn aspect_ratio(tet: [[f64; 3]; 4]) -> f64 {
        let edges = Self::all_edges(tet);
        let max_edge = edges.iter().cloned().fold(0.0_f64, f64::max);
        let vol = Self::volume(tet).abs();
        if vol < 1e-30 {
            return f64::INFINITY;
        }
        let face_area = Self::total_face_area(tet);
        let inradius = 3.0 * vol / face_area;
        if inradius < 1e-30 {
            return f64::INFINITY;
        }
        max_edge / (3.0 * inradius)
    }
    /// Compute the minimum dihedral angle (radians) across all 6 edges of a tetrahedron.
    pub fn min_dihedral_angle(tet: [[f64; 3]; 4]) -> f64 {
        const EDGE_PAIRS: [(usize, usize, usize, usize); 6] = [
            (0, 1, 2, 3),
            (0, 2, 1, 3),
            (0, 3, 1, 2),
            (1, 2, 0, 3),
            (1, 3, 0, 2),
            (2, 3, 0, 1),
        ];
        let mut min_angle = f64::INFINITY;
        for (a, b, c, d) in EDGE_PAIRS {
            let angle = Self::dihedral_angle(tet[a], tet[b], tet[c], tet[d]);
            if angle < min_angle {
                min_angle = angle;
            }
        }
        min_angle
    }
    /// Compute the Jacobian determinant of the linear mapping from the reference tet.
    pub fn jacobian_determinant(tet: [[f64; 3]; 4]) -> f64 {
        let e1 = sub3(tet[1], tet[0]);
        let e2 = sub3(tet[2], tet[0]);
        let e3 = sub3(tet[3], tet[0]);
        e1[0] * (e2[1] * e3[2] - e2[2] * e3[1]) - e1[1] * (e2[0] * e3[2] - e2[2] * e3[0])
            + e1[2] * (e2[0] * e3[1] - e2[1] * e3[0])
    }
    /// Composite element quality score in \[0, 1\] (1 = perfect regular tet).
    ///
    /// Uses the mean-ratio metric: η = 12·(3V)^(2/3) / (Σ|eᵢ|²).
    pub fn element_quality(tet: [[f64; 3]; 4]) -> f64 {
        let vol = Self::volume(tet).abs();
        if vol < 1e-30 {
            return 0.0;
        }
        let edges = Self::all_edges(tet);
        let sum_sq: f64 = edges.iter().map(|e| e * e).sum();
        if sum_sq < 1e-30 {
            return 0.0;
        }
        let numerator = 12.0 * (3.0 * vol).powf(2.0 / 3.0);
        (numerator / sum_sq).clamp(0.0, 1.0)
    }
    /// Signed volume of a tetrahedron (positive when nodes are in CCW order).
    pub fn volume(tet: [[f64; 3]; 4]) -> f64 {
        let e1 = sub3(tet[1], tet[0]);
        let e2 = sub3(tet[2], tet[0]);
        let e3 = sub3(tet[3], tet[0]);
        dot3(cross3(e1, e2), e3) / 6.0
    }
    fn all_edges(tet: [[f64; 3]; 4]) -> [f64; 6] {
        [
            norm3(sub3(tet[1], tet[0])),
            norm3(sub3(tet[2], tet[0])),
            norm3(sub3(tet[3], tet[0])),
            norm3(sub3(tet[2], tet[1])),
            norm3(sub3(tet[3], tet[1])),
            norm3(sub3(tet[3], tet[2])),
        ]
    }
    fn total_face_area(tet: [[f64; 3]; 4]) -> f64 {
        let faces: [[usize; 3]; 4] = [[1, 2, 3], [0, 2, 3], [0, 1, 3], [0, 1, 2]];
        let mut total = 0.0;
        for f in &faces {
            let ab = sub3(tet[f[1]], tet[f[0]]);
            let ac = sub3(tet[f[2]], tet[f[0]]);
            total += norm3(cross3(ab, ac)) * 0.5;
        }
        total
    }
    /// Dihedral angle along the edge (a,b) shared by the triangles (a,b,c) and (a,b,d).
    fn dihedral_angle(a: [f64; 3], b: [f64; 3], c: [f64; 3], d: [f64; 3]) -> f64 {
        let ab = sub3(b, a);
        let ac = sub3(c, a);
        let ad = sub3(d, a);
        let n1 = cross3(ab, ac);
        let n2 = cross3(ab, ad);
        let cos_angle = dot3(n1, n2) / (norm3(n1) * norm3(n2)).max(1e-30);
        cos_angle.clamp(-1.0, 1.0).acos()
    }
}
/// Red-green refinement for conforming tetrahedral meshes.
///
/// - **Red refinement**: bisect all 6 edges of a tetrahedron into 8 children.
/// - **Green refinement**: bisect one edge of a neighbouring tet to restore conformity.
pub struct LocalRefinement {
    /// Node positions (shared with the mesh).
    pub nodes: Vec<[f64; 3]>,
    /// Tetrahedral connectivity.
    pub tets: Vec<[usize; 4]>,
}
impl LocalRefinement {
    /// Create a new `LocalRefinement` instance.
    pub fn new(nodes: Vec<[f64; 3]>, tets: Vec<[usize; 4]>) -> Self {
        Self { nodes, tets }
    }
    /// Red-refine tetrahedron `tet_idx`: bisect all 6 edges and return the 8 child tet indices.
    ///
    /// The original tetrahedron is removed from `self.tets` and replaced by 8 children.
    /// Returns the indices of the 8 new tetrahedra in `self.tets`.
    pub fn red_refine(&mut self, tet_idx: usize) -> [usize; 8] {
        if tet_idx >= self.tets.len() {
            return [0; 8];
        }
        let tet = self.tets[tet_idx];
        let m01 = self.get_or_insert_midpoint(tet[0], tet[1]);
        let m02 = self.get_or_insert_midpoint(tet[0], tet[2]);
        let m03 = self.get_or_insert_midpoint(tet[0], tet[3]);
        let m12 = self.get_or_insert_midpoint(tet[1], tet[2]);
        let m13 = self.get_or_insert_midpoint(tet[1], tet[3]);
        let m23 = self.get_or_insert_midpoint(tet[2], tet[3]);
        self.tets.remove(tet_idx);
        let children: [[usize; 4]; 8] = [
            [tet[0], m01, m02, m03],
            [tet[1], m01, m12, m13],
            [tet[2], m02, m12, m23],
            [tet[3], m03, m13, m23],
            [m01, m02, m12, m23],
            [m01, m12, m13, m23],
            [m01, m02, m03, m23],
            [m01, m03, m13, m23],
        ];
        let base = self.tets.len();
        for child in &children {
            self.tets.push(*child);
        }
        [
            base,
            base + 1,
            base + 2,
            base + 3,
            base + 4,
            base + 5,
            base + 6,
            base + 7,
        ]
    }
    /// Green-refine tetrahedron `tet_idx` along the given `edge_idx` (0–5).
    ///
    /// Bisects the single specified edge and returns the indices of the 2 new tetrahedra.
    pub fn green_refine(&mut self, tet_idx: usize, edge_idx: usize) -> [usize; 2] {
        if tet_idx >= self.tets.len() {
            return [0; 2];
        }
        let tet = self.tets[tet_idx];
        const EDGE_TABLE: [(usize, usize); 6] = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)];
        let (ea, eb) = EDGE_TABLE[edge_idx.min(5)];
        let mid = self.get_or_insert_midpoint(tet[ea], tet[eb]);
        let others: Vec<usize> = tet
            .iter()
            .enumerate()
            .filter(|&(k, _)| k != ea && k != eb)
            .map(|(_, &n)| n)
            .collect();
        self.tets.remove(tet_idx);
        let o0 = others.first().copied().unwrap_or(0);
        let o1 = others.get(1).copied().unwrap_or(0);
        let base = self.tets.len();
        self.tets.push([tet[ea], mid, o0, o1]);
        self.tets.push([tet[eb], mid, o0, o1]);
        [base, base + 1]
    }
    fn get_or_insert_midpoint(&mut self, i: usize, j: usize) -> usize {
        let mp = midpoint3(self.nodes[i], self.nodes[j]);
        for (k, &n) in self.nodes.iter().enumerate() {
            if (n[0] - mp[0]).abs() < 1e-14
                && (n[1] - mp[1]).abs() < 1e-14
                && (n[2] - mp[2]).abs() < 1e-14
            {
                return k;
            }
        }
        self.nodes.push(mp);
        self.nodes.len() - 1
    }
}
/// Priority queue of element indices awaiting refinement.
pub struct RefinementQueue {
    /// Internal storage of `(element_index, priority)` pairs.
    pub(super) heap: Vec<(usize, f64)>,
}
impl RefinementQueue {
    /// Create an empty refinement queue.
    pub fn new() -> Self {
        Self { heap: Vec::new() }
    }
    /// Push an element index with the given priority (higher = refined first).
    pub fn push(&mut self, elem_idx: usize, priority: f64) {
        self.heap.push((elem_idx, priority));
        let n = self.heap.len();
        let mut i = n - 1;
        while i > 0 && self.heap[i].1 > self.heap[i - 1].1 {
            self.heap.swap(i, i - 1);
            i -= 1;
        }
    }
    /// Pop the highest-priority element index.
    pub fn pop(&mut self) -> Option<usize> {
        if self.heap.is_empty() {
            None
        } else {
            Some(self.heap.remove(0).0)
        }
    }
    /// Mark all elements whose error exceeds `tol` by pushing them with error as priority.
    pub fn mark_all_above_threshold(&mut self, errors: &[f64], tol: f64) {
        for (idx, &err) in errors.iter().enumerate() {
            if err > tol {
                self.push(idx, err);
            }
        }
    }
    /// Number of elements in the queue.
    pub fn len(&self) -> usize {
        self.heap.len()
    }
    /// Returns `true` when no elements are queued.
    pub fn is_empty(&self) -> bool {
        self.heap.is_empty()
    }
}
/// Criterion used to decide which elements should be refined.
#[derive(Debug, Clone)]
pub enum RefinementCriterion {
    /// Refine elements whose error indicator exceeds `threshold`.
    ErrorBased {
        /// Error threshold above which an element is refined.
        threshold: f64,
    },
    /// Refine elements where the solution gradient magnitude exceeds `threshold`.
    GradientBased {
        /// Gradient magnitude threshold.
        threshold: f64,
    },
    /// Refine elements whose longest edge exceeds `target_size`.
    GeometricBased {
        /// Target edge length for refinement.
        target_size: f64,
    },
    /// Refine all elements uniformly (one level).
    Uniform,
}
/// A node in an octree spatial index over tetrahedral elements.
pub struct OctreeNode {
    /// Center of this octree cell.
    pub center: [f64; 3],
    /// Half-extent (half-width) of this cubic cell.
    pub half: f64,
    /// Eight child octants (None = leaf or empty).
    pub children: [Option<Box<OctreeNode>>; 8],
    /// Element indices stored at this node.
    pub elements: Vec<usize>,
}
impl OctreeNode {
    /// Create a new empty leaf octree node.
    pub fn new(center: [f64; 3], half: f64) -> Self {
        Self {
            center,
            half,
            children: [None, None, None, None, None, None, None, None],
            elements: Vec::new(),
        }
    }
    /// Insert element `elem` whose axis-aligned bounding box spans \[`bbox_min`, `bbox_max`\].
    pub fn insert(
        &mut self,
        elem: usize,
        bbox_min: [f64; 3],
        bbox_max: [f64; 3],
        max_depth: usize,
    ) {
        if max_depth == 0 || !self.aabb_overlaps_cell(bbox_min, bbox_max) {
            self.elements.push(elem);
            return;
        }
        let child_half = self.half * 0.5;
        let mut inserted = false;
        for oct in 0..8usize {
            let child_center = self.child_center(oct);
            let child_min = [
                child_center[0] - child_half,
                child_center[1] - child_half,
                child_center[2] - child_half,
            ];
            let child_max = [
                child_center[0] + child_half,
                child_center[1] + child_half,
                child_center[2] + child_half,
            ];
            if bbox_max[0] >= child_min[0]
                && bbox_min[0] <= child_max[0]
                && bbox_max[1] >= child_min[1]
                && bbox_min[1] <= child_max[1]
                && bbox_max[2] >= child_min[2]
                && bbox_min[2] <= child_max[2]
            {
                if self.children[oct].is_none() {
                    self.children[oct] = Some(Box::new(OctreeNode::new(child_center, child_half)));
                }
                self.children[oct]
                    .as_mut()
                    .expect("child node was just created")
                    .insert(elem, bbox_min, bbox_max, max_depth - 1);
                inserted = true;
            }
        }
        if !inserted {
            self.elements.push(elem);
        }
    }
    /// Query all elements whose bounding box may overlap a sphere at `center` with `radius`.
    pub fn query_sphere(&self, center: [f64; 3], radius: f64) -> Vec<usize> {
        let mut result = Vec::new();
        if !self.sphere_overlaps_cell(center, radius) {
            return result;
        }
        result.extend_from_slice(&self.elements);
        for c in self.children.iter().flatten() {
            result.extend(c.query_sphere(center, radius));
        }
        result
    }
    /// Query all elements whose bounding box may overlap the given AABB \[`min`, `max`\].
    pub fn query_box(&self, min: [f64; 3], max: [f64; 3]) -> Vec<usize> {
        let mut result = Vec::new();
        if !self.aabb_overlaps_cell(min, max) {
            return result;
        }
        result.extend_from_slice(&self.elements);
        for c in self.children.iter().flatten() {
            result.extend(c.query_box(min, max));
        }
        result
    }
    fn child_center(&self, oct: usize) -> [f64; 3] {
        let h = self.half * 0.5;
        let dx = if oct & 1 != 0 { h } else { -h };
        let dy = if oct & 2 != 0 { h } else { -h };
        let dz = if oct & 4 != 0 { h } else { -h };
        [
            self.center[0] + dx,
            self.center[1] + dy,
            self.center[2] + dz,
        ]
    }
    fn aabb_overlaps_cell(&self, min: [f64; 3], max: [f64; 3]) -> bool {
        let cmin = [
            self.center[0] - self.half,
            self.center[1] - self.half,
            self.center[2] - self.half,
        ];
        let cmax = [
            self.center[0] + self.half,
            self.center[1] + self.half,
            self.center[2] + self.half,
        ];
        max[0] >= cmin[0]
            && min[0] <= cmax[0]
            && max[1] >= cmin[1]
            && min[1] <= cmax[1]
            && max[2] >= cmin[2]
            && min[2] <= cmax[2]
    }
    fn sphere_overlaps_cell(&self, center: [f64; 3], radius: f64) -> bool {
        let cmin = [
            self.center[0] - self.half,
            self.center[1] - self.half,
            self.center[2] - self.half,
        ];
        let cmax = [
            self.center[0] + self.half,
            self.center[1] + self.half,
            self.center[2] + self.half,
        ];
        let dx = clamp_f64(center[0], cmin[0], cmax[0]) - center[0];
        let dy = clamp_f64(center[1], cmin[1], cmax[1]) - center[1];
        let dz = clamp_f64(center[2], cmin[2], cmax[2]) - center[2];
        dx * dx + dy * dy + dz * dz <= radius * radius
    }
}
