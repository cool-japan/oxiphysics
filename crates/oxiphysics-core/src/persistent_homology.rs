// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Persistent homology and topological data analysis.
//!
//! Provides simplicial complexes, filtered complexes, persistence diagrams,
//! Vietoris-Rips filtration, persistence images, and barcode statistics.

use std::collections::HashMap;

// ─────────────────────────────────────────────────────────────────────────────
// Simplex
// ─────────────────────────────────────────────────────────────────────────────

/// A simplex defined by its vertices, dimension, and filtration value.
///
/// A `k`-simplex has `k+1` vertices (e.g., 0-simplex = vertex, 1-simplex = edge).
#[derive(Debug, Clone, PartialEq)]
pub struct Simplex {
    /// Sorted vertex indices forming this simplex.
    pub vertices: Vec<usize>,
    /// Homological dimension (0 = vertex, 1 = edge, 2 = triangle, …).
    pub dimension: usize,
    /// Filtration value at which this simplex enters the complex.
    pub filtration: f64,
}

impl Simplex {
    /// Create a new simplex from a list of vertices and filtration value.
    ///
    /// The dimension is inferred as `vertices.len() - 1`.
    /// Vertices are sorted in ascending order.
    pub fn new(mut vertices: Vec<usize>, filtration: f64) -> Self {
        vertices.sort_unstable();
        let dimension = vertices.len().saturating_sub(1);
        Self {
            vertices,
            dimension,
            filtration,
        }
    }

    /// Return the faces (boundary simplices) of this simplex.
    ///
    /// A `k`-simplex has `k+1` faces of dimension `k-1`, each formed by
    /// removing one vertex.
    pub fn faces(&self) -> Vec<Simplex> {
        if self.dimension == 0 {
            return vec![];
        }
        (0..self.vertices.len())
            .map(|i| {
                let verts: Vec<usize> = self
                    .vertices
                    .iter()
                    .enumerate()
                    .filter(|(j, _)| *j != i)
                    .map(|(_, &v)| v)
                    .collect();
                Simplex::new(verts, self.filtration)
            })
            .collect()
    }

    /// Return `true` if this simplex contains all vertices of `other`.
    pub fn contains(&self, other: &Simplex) -> bool {
        other.vertices.iter().all(|v| self.vertices.contains(v))
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// SimplicialComplex
// ─────────────────────────────────────────────────────────────────────────────

/// A finite simplicial complex with boundary operators and topological invariants.
#[derive(Debug, Clone)]
pub struct SimplicialComplex {
    /// All simplices grouped by dimension.
    pub simplices: Vec<Vec<Simplex>>,
    /// Maximum dimension present in the complex.
    pub max_dimension: usize,
}

impl SimplicialComplex {
    /// Create an empty simplicial complex.
    pub fn new() -> Self {
        Self {
            simplices: vec![],
            max_dimension: 0,
        }
    }

    /// Add a simplex (and all its faces recursively) to the complex.
    pub fn add_simplex(&mut self, s: Simplex) {
        // First ensure faces are present
        let faces = s.faces();
        for face in faces {
            self.add_simplex(face);
        }
        let dim = s.dimension;
        while self.simplices.len() <= dim {
            self.simplices.push(vec![]);
        }
        if dim > self.max_dimension {
            self.max_dimension = dim;
        }
        // Avoid duplicates
        if !self.simplices[dim].iter().any(|x| x.vertices == s.vertices) {
            self.simplices[dim].push(s);
        }
    }

    /// Return the number of simplices of dimension `k`.
    pub fn count(&self, k: usize) -> usize {
        self.simplices.get(k).map_or(0, |v| v.len())
    }

    /// Compute the Euler characteristic χ = Σ (-1)^k * f_k.
    ///
    /// Where `f_k` is the number of `k`-simplices.
    pub fn euler_characteristic(&self) -> i64 {
        self.simplices
            .iter()
            .enumerate()
            .map(|(k, s)| {
                if k % 2 == 0 {
                    s.len() as i64
                } else {
                    -(s.len() as i64)
                }
            })
            .sum()
    }

    /// Compute the boundary matrix for dimension `k`.
    ///
    /// Entry `[i][j]` is 1 if the `j`-th `(k-1)`-simplex is a face of the
    /// `i`-th `k`-simplex (with sign ignored for simplicity), else 0.
    pub fn boundary_matrix(&self, k: usize) -> Vec<Vec<i32>> {
        if k == 0 || k > self.max_dimension {
            return vec![];
        }
        let k_simplices = &self.simplices[k];
        let km1_simplices = &self.simplices[k - 1];
        let rows = k_simplices.len();
        let cols = km1_simplices.len();
        let mut mat = vec![vec![0i32; cols]; rows];
        for (i, sigma) in k_simplices.iter().enumerate() {
            for (j, tau) in km1_simplices.iter().enumerate() {
                if sigma.contains(tau) {
                    mat[i][j] = 1;
                }
            }
        }
        mat
    }
}

impl Default for SimplicialComplex {
    fn default() -> Self {
        Self::new()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// FilteredComplex
// ─────────────────────────────────────────────────────────────────────────────

/// A simplicial complex with a filtration ordering.
///
/// Simplices are ordered by their filtration values, and then by dimension.
/// The boundary matrix algorithm operates on this ordering.
#[derive(Debug, Clone)]
pub struct FilteredComplex {
    /// All simplices in filtration order (sorted by `filtration` then `dimension`).
    pub ordered_simplices: Vec<Simplex>,
}

impl FilteredComplex {
    /// Create a filtered complex from a list of simplices.
    ///
    /// Simplices are sorted by filtration value then dimension.
    pub fn new(mut simplices: Vec<Simplex>) -> Self {
        simplices.sort_by(|a, b| {
            a.filtration
                .partial_cmp(&b.filtration)
                .unwrap_or(std::cmp::Ordering::Equal)
                .then(a.dimension.cmp(&b.dimension))
        });
        Self {
            ordered_simplices: simplices,
        }
    }

    /// Return the number of simplices in this complex.
    pub fn len(&self) -> usize {
        self.ordered_simplices.len()
    }

    /// Return `true` if the complex is empty.
    pub fn is_empty(&self) -> bool {
        self.ordered_simplices.is_empty()
    }

    /// Compute the boundary matrix column for simplex at index `i`.
    ///
    /// Returns a sorted vector of indices `j < i` such that
    /// `ordered_simplices[j]` is a codimension-1 face of `ordered_simplices[i]`.
    pub fn boundary_column(&self, i: usize) -> Vec<usize> {
        let sigma = &self.ordered_simplices[i];
        if sigma.dimension == 0 {
            return vec![];
        }
        let faces = sigma.faces();
        let mut col = Vec::new();
        for face in &faces {
            if let Some(j) = self.ordered_simplices[..i]
                .iter()
                .rposition(|s| s.vertices == face.vertices)
            {
                col.push(j);
            }
        }
        col.sort_unstable();
        col
    }

    /// Run the standard persistence algorithm (column reduction).
    ///
    /// Returns the pivot pairs `(i, j)` where column `j` has pivot row `i`,
    /// meaning simplex `j` dies when simplex `i` is born.
    pub fn reduce(&self) -> Vec<(usize, usize)> {
        let n = self.ordered_simplices.len();
        // Reduced columns stored as sorted Vec<usize>
        let mut columns: Vec<Vec<usize>> = (0..n).map(|i| self.boundary_column(i)).collect();
        // low[j] = lowest set bit index in column j, or None
        let low = |col: &Vec<usize>| col.last().copied();
        let mut pivot_col: HashMap<usize, usize> = HashMap::new();
        let mut pairs = Vec::new();

        for j in 0..n {
            loop {
                if columns[j].is_empty() {
                    break;
                }
                let l = low(&columns[j]).expect("column j is non-empty");
                if let Some(&k) = pivot_col.get(&l) {
                    // Add column k to column j (XOR-style over Z/2Z)
                    let col_k = columns[k].clone();
                    let col_j = std::mem::take(&mut columns[j]);
                    columns[j] = xor_sorted_vecs(col_j, col_k);
                } else {
                    pivot_col.insert(l, j);
                    pairs.push((l, j));
                    break;
                }
            }
        }
        pairs
    }
}

/// Symmetric difference of two sorted vectors (XOR over Z/2Z).
fn xor_sorted_vecs(mut a: Vec<usize>, mut b: Vec<usize>) -> Vec<usize> {
    a.sort_unstable();
    b.sort_unstable();
    let mut result = Vec::new();
    let (mut i, mut j) = (0, 0);
    while i < a.len() && j < b.len() {
        match a[i].cmp(&b[j]) {
            std::cmp::Ordering::Less => {
                result.push(a[i]);
                i += 1;
            }
            std::cmp::Ordering::Greater => {
                result.push(b[j]);
                j += 1;
            }
            std::cmp::Ordering::Equal => {
                // Cancel (both drop out)
                i += 1;
                j += 1;
            }
        }
    }
    result.extend_from_slice(&a[i..]);
    result.extend_from_slice(&b[j..]);
    result
}

// ─────────────────────────────────────────────────────────────────────────────
// PersistenceDiagram
// ─────────────────────────────────────────────────────────────────────────────

/// A persistence diagram: collection of birth-death pairs per dimension.
#[derive(Debug, Clone)]
pub struct PersistenceDiagram {
    /// Birth-death pairs `(birth, death)` per homological dimension.
    pub pairs: Vec<(f64, f64, usize)>,
    /// Essential classes (never destroyed): `(birth, dimension)`.
    pub essential: Vec<(f64, usize)>,
}

impl PersistenceDiagram {
    /// Create a `PersistenceDiagram` from a `FilteredComplex` using standard reduction.
    pub fn from_filtered_complex(fc: &FilteredComplex) -> Self {
        let pairs_idx = fc.reduce();
        let n = fc.ordered_simplices.len();
        let mut paired = vec![false; n];
        let mut pairs = Vec::new();
        let mut essential = Vec::new();

        for (birth_idx, death_idx) in &pairs_idx {
            let b = fc.ordered_simplices[*birth_idx].filtration;
            let d = fc.ordered_simplices[*death_idx].filtration;
            let dim = fc.ordered_simplices[*birth_idx].dimension;
            if (d - b).abs() > 1e-15 {
                pairs.push((b, d, dim));
            }
            paired[*birth_idx] = true;
            paired[*death_idx] = true;
        }
        // Unpaired simplices → essential classes
        for (i, s) in fc.ordered_simplices.iter().enumerate() {
            if !paired[i] {
                essential.push((s.filtration, s.dimension));
            }
        }
        Self { pairs, essential }
    }

    /// Return Betti numbers up to dimension `max_dim` at filtration value `t`.
    ///
    /// β_k(t) = number of classes born at or before `t` that are still alive at `t`.
    pub fn betti_numbers(&self, t: f64, max_dim: usize) -> Vec<usize> {
        let mut betti = vec![0usize; max_dim + 1];
        for &(b, d, dim) in &self.pairs {
            if dim <= max_dim && b <= t && d > t {
                betti[dim] += 1;
            }
        }
        for &(b, dim) in &self.essential {
            if dim <= max_dim && b <= t {
                betti[dim] += 1;
            }
        }
        betti
    }

    /// Return total persistence (sum of lifetimes for all finite pairs).
    pub fn total_persistence(&self) -> f64 {
        self.pairs.iter().map(|(b, d, _)| d - b).sum()
    }

    /// Return all finite pairs for dimension `k`.
    pub fn pairs_for_dim(&self, k: usize) -> Vec<(f64, f64)> {
        self.pairs
            .iter()
            .filter(|&&(_, _, d)| d == k)
            .map(|&(b, d, _)| (b, d))
            .collect()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// RipsFiltration
// ─────────────────────────────────────────────────────────────────────────────

/// Vietoris-Rips filtration built from a point cloud or distance matrix.
#[derive(Debug, Clone)]
pub struct RipsFiltration {
    /// Distance matrix (symmetric, zero diagonal).
    pub distances: Vec<Vec<f64>>,
    /// Number of points.
    pub n_points: usize,
    /// Maximum filtration radius.
    pub max_radius: f64,
    /// Maximum simplex dimension to build.
    pub max_dim: usize,
}

impl RipsFiltration {
    /// Build a Rips filtration from a 2D point cloud (Euclidean distances).
    ///
    /// `max_radius` controls how far to grow the filtration.
    /// `max_dim` controls the maximum simplex dimension built.
    pub fn from_point_cloud(points: &[[f64; 2]], max_radius: f64, max_dim: usize) -> Self {
        let n = points.len();
        let mut distances = vec![vec![0.0_f64; n]; n];
        for i in 0..n {
            for j in (i + 1)..n {
                let dx = points[i][0] - points[j][0];
                let dy = points[i][1] - points[j][1];
                let d = (dx * dx + dy * dy).sqrt();
                distances[i][j] = d;
                distances[j][i] = d;
            }
        }
        Self {
            distances,
            n_points: n,
            max_radius,
            max_dim,
        }
    }

    /// Build a Rips filtration from a precomputed distance matrix.
    pub fn from_distance_matrix(distances: Vec<Vec<f64>>, max_radius: f64, max_dim: usize) -> Self {
        let n = distances.len();
        Self {
            distances,
            n_points: n,
            max_radius,
            max_dim,
        }
    }

    /// Build a `FilteredComplex` from this filtration.
    pub fn build_complex(&self) -> FilteredComplex {
        let mut simplices = Vec::new();
        let n = self.n_points;

        // 0-simplices (vertices) — born at filtration 0
        for i in 0..n {
            simplices.push(Simplex::new(vec![i], 0.0));
        }

        // 1-simplices (edges)
        for i in 0..n {
            for j in (i + 1)..n {
                let d = self.distances[i][j];
                if d <= self.max_radius {
                    simplices.push(Simplex::new(vec![i, j], d));
                }
            }
        }

        // 2-simplices (triangles) if max_dim >= 2
        if self.max_dim >= 2 {
            for i in 0..n {
                for j in (i + 1)..n {
                    for k in (j + 1)..n {
                        let d = self.distances[i][j]
                            .max(self.distances[i][k])
                            .max(self.distances[j][k]);
                        if d <= self.max_radius {
                            simplices.push(Simplex::new(vec![i, j, k], d));
                        }
                    }
                }
            }
        }

        // 3-simplices (tetrahedra) if max_dim >= 3
        if self.max_dim >= 3 && n >= 4 {
            for i in 0..n {
                for j in (i + 1)..n {
                    for k in (j + 1)..n {
                        for l in (k + 1)..n {
                            let d = self.distances[i][j]
                                .max(self.distances[i][k])
                                .max(self.distances[i][l])
                                .max(self.distances[j][k])
                                .max(self.distances[j][l])
                                .max(self.distances[k][l]);
                            if d <= self.max_radius {
                                simplices.push(Simplex::new(vec![i, j, k, l], d));
                            }
                        }
                    }
                }
            }
        }

        FilteredComplex::new(simplices)
    }

    /// Compute persistence diagram for this Rips filtration.
    pub fn persistence_diagram(&self) -> PersistenceDiagram {
        let fc = self.build_complex();
        PersistenceDiagram::from_filtered_complex(&fc)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// PersistenceImage
// ─────────────────────────────────────────────────────────────────────────────

/// Persistence image: a pixelated representation of a persistence diagram.
///
/// Points are represented in birth-persistence coordinates, then convolved
/// with a Gaussian kernel and discretized onto a pixel grid.
#[derive(Debug, Clone)]
pub struct PersistenceImage {
    /// Pixel grid (row-major, rows = persistence axis, cols = birth axis).
    pub pixels: Vec<Vec<f64>>,
    /// Number of grid rows (persistence resolution).
    pub n_rows: usize,
    /// Number of grid cols (birth resolution).
    pub n_cols: usize,
    /// Minimum birth value.
    pub birth_min: f64,
    /// Maximum birth value.
    pub birth_max: f64,
    /// Maximum persistence value.
    pub pers_max: f64,
    /// Gaussian bandwidth (sigma).
    pub sigma: f64,
}

impl PersistenceImage {
    /// Build a persistence image from a persistence diagram.
    ///
    /// `n_rows` and `n_cols` control resolution.
    /// `sigma` is the Gaussian kernel bandwidth.
    pub fn from_diagram(
        diagram: &PersistenceDiagram,
        n_rows: usize,
        n_cols: usize,
        sigma: f64,
    ) -> Self {
        // Collect finite pairs in birth-persistence coords
        let pts: Vec<(f64, f64)> = diagram
            .pairs
            .iter()
            .filter(|&&(b, d, _)| d.is_finite() && d > b)
            .map(|&(b, d, _)| (b, d - b))
            .collect();

        if pts.is_empty() {
            return Self {
                pixels: vec![vec![0.0; n_cols]; n_rows],
                n_rows,
                n_cols,
                birth_min: 0.0,
                birth_max: 1.0,
                pers_max: 1.0,
                sigma,
            };
        }

        let birth_min = pts.iter().map(|(b, _)| *b).fold(f64::INFINITY, f64::min);
        let birth_max = pts
            .iter()
            .map(|(b, _)| *b)
            .fold(f64::NEG_INFINITY, f64::max);
        let pers_max = pts
            .iter()
            .map(|(_, p)| *p)
            .fold(f64::NEG_INFINITY, f64::max);

        let birth_range = (birth_max - birth_min).max(1e-10);
        let pers_range = pers_max.max(1e-10);

        let mut pixels = vec![vec![0.0_f64; n_cols]; n_rows];
        let two_sigma_sq = 2.0 * sigma * sigma;

        for (bx, py) in &pts {
            // Ramp weight: persistence / max_persistence
            let weight = py / pers_range;
            for (row, pixel_row) in pixels.iter_mut().enumerate() {
                // row 0 = high persistence
                let p_grid =
                    pers_range * (n_rows - 1 - row) as f64 / (n_rows as f64 - 1.0).max(1.0);
                for (col, pixel) in pixel_row.iter_mut().enumerate() {
                    let b_grid =
                        birth_min + birth_range * col as f64 / (n_cols as f64 - 1.0).max(1.0);
                    let db = bx - b_grid;
                    let dp = py - p_grid;
                    let gauss = (-(db * db + dp * dp) / two_sigma_sq).exp();
                    *pixel += weight * gauss;
                }
            }
        }

        Self {
            pixels,
            n_rows,
            n_cols,
            birth_min,
            birth_max,
            pers_max,
            sigma,
        }
    }

    /// Flatten the pixel grid to a 1D feature vector (row-major).
    pub fn to_vector(&self) -> Vec<f64> {
        self.pixels
            .iter()
            .flat_map(|row| row.iter().copied())
            .collect()
    }

    /// Return the maximum pixel value.
    pub fn max_value(&self) -> f64 {
        self.pixels
            .iter()
            .flat_map(|row| row.iter())
            .cloned()
            .fold(f64::NEG_INFINITY, f64::max)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// BarcodeStatistics
// ─────────────────────────────────────────────────────────────────────────────

/// Summary statistics computed from a persistence barcode.
#[derive(Debug, Clone)]
pub struct BarcodeStatistics {
    /// Sum of finite bar lengths.
    pub total_persistence: f64,
    /// Mean finite bar length.
    pub mean_persistence: f64,
    /// Median finite bar length.
    pub median_persistence: f64,
    /// Maximum finite bar length.
    pub max_persistence: f64,
    /// Number of finite bars.
    pub n_finite: usize,
    /// Number of essential (infinite) classes.
    pub n_essential: usize,
    /// Bottleneck distance approximation to the empty diagram.
    ///
    /// Equals the half-persistence of the most persistent finite bar.
    pub bottleneck_to_empty: f64,
}

impl BarcodeStatistics {
    /// Compute statistics from a `PersistenceDiagram`.
    pub fn from_diagram(diagram: &PersistenceDiagram) -> Self {
        let mut lifetimes: Vec<f64> = diagram
            .pairs
            .iter()
            .filter(|&&(b, d, _)| d.is_finite() && d > b)
            .map(|&(b, d, _)| d - b)
            .collect();
        lifetimes.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));

        let n_finite = lifetimes.len();
        let total_persistence = lifetimes.iter().sum::<f64>();
        let mean_persistence = if n_finite > 0 {
            total_persistence / n_finite as f64
        } else {
            0.0
        };
        let median_persistence = if n_finite > 0 {
            if n_finite % 2 == 1 {
                lifetimes[n_finite / 2]
            } else {
                (lifetimes[n_finite / 2 - 1] + lifetimes[n_finite / 2]) / 2.0
            }
        } else {
            0.0
        };
        let max_persistence = lifetimes.last().copied().unwrap_or(0.0);
        let bottleneck_to_empty = max_persistence / 2.0;
        let n_essential = diagram.essential.len();

        Self {
            total_persistence,
            mean_persistence,
            median_persistence,
            max_persistence,
            n_finite,
            n_essential,
            bottleneck_to_empty,
        }
    }

    /// Approximate bottleneck distance between two diagrams.
    ///
    /// Uses greedy matching by L∞ distance in birth-persistence coordinates.
    pub fn bottleneck_distance(a: &PersistenceDiagram, b: &PersistenceDiagram) -> f64 {
        let pts_a: Vec<(f64, f64)> = a
            .pairs
            .iter()
            .filter(|&&(b_, d, _)| d.is_finite() && d > b_)
            .map(|&(b_, d, _)| (b_, d - b_))
            .collect();
        let pts_b: Vec<(f64, f64)> = b
            .pairs
            .iter()
            .filter(|&&(b_, d, _)| d.is_finite() && d > b_)
            .map(|&(b_, d, _)| (b_, d - b_))
            .collect();

        let linf = |pa: (f64, f64), pb: (f64, f64)| -> f64 {
            (pa.0 - pb.0).abs().max((pa.1 - pb.1).abs())
        };
        let diag_dist = |p: (f64, f64)| -> f64 { p.1 / 2.0 };

        let mut max_dist = 0.0_f64;
        let mut used = vec![false; pts_b.len()];

        for &pa in &pts_a {
            let best = pts_b
                .iter()
                .enumerate()
                .filter(|(idx, _)| !used[*idx])
                .map(|(idx, &pb)| (idx, linf(pa, pb)))
                .min_by(|x, y| x.1.partial_cmp(&y.1).unwrap_or(std::cmp::Ordering::Equal));
            match best {
                Some((idx, d)) if d < diag_dist(pa) => {
                    used[idx] = true;
                    max_dist = max_dist.max(d);
                }
                _ => {
                    max_dist = max_dist.max(diag_dist(pa));
                }
            }
        }
        for (idx, &pb) in pts_b.iter().enumerate() {
            if !used[idx] {
                max_dist = max_dist.max(diag_dist(pb));
            }
        }
        max_dist
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // ── Simplex ───────────────────────────────────────────────────────────────

    #[test]
    fn test_simplex_vertex() {
        let s = Simplex::new(vec![3], 0.5);
        assert_eq!(s.dimension, 0);
        assert_eq!(s.vertices, vec![3]);
        assert!((s.filtration - 0.5_f64).abs() < 1e-12);
    }

    #[test]
    fn test_simplex_edge_dimension() {
        let s = Simplex::new(vec![0, 1], 1.0);
        assert_eq!(s.dimension, 1);
    }

    #[test]
    fn test_simplex_triangle_dimension() {
        let s = Simplex::new(vec![0, 1, 2], 2.0);
        assert_eq!(s.dimension, 2);
    }

    #[test]
    fn test_simplex_vertices_sorted() {
        let s = Simplex::new(vec![3, 1, 0, 2], 0.0);
        assert_eq!(s.vertices, vec![0, 1, 2, 3]);
    }

    #[test]
    fn test_simplex_vertex_has_no_faces() {
        let s = Simplex::new(vec![5], 0.0);
        assert!(s.faces().is_empty());
    }

    #[test]
    fn test_simplex_edge_has_two_faces() {
        let s = Simplex::new(vec![0, 1], 1.0);
        let faces = s.faces();
        assert_eq!(faces.len(), 2);
        assert!(faces.iter().all(|f| f.dimension == 0));
    }

    #[test]
    fn test_simplex_triangle_has_three_faces() {
        let s = Simplex::new(vec![0, 1, 2], 1.5);
        let faces = s.faces();
        assert_eq!(faces.len(), 3);
        assert!(faces.iter().all(|f| f.dimension == 1));
    }

    #[test]
    fn test_simplex_contains() {
        let edge = Simplex::new(vec![0, 1], 1.0);
        let v0 = Simplex::new(vec![0], 0.0);
        let v2 = Simplex::new(vec![2], 0.0);
        assert!(edge.contains(&v0));
        assert!(!edge.contains(&v2));
    }

    // ── SimplicialComplex ─────────────────────────────────────────────────────

    #[test]
    fn test_simplicial_complex_empty() {
        let sc = SimplicialComplex::new();
        assert_eq!(sc.count(0), 0);
        assert_eq!(sc.euler_characteristic(), 0);
    }

    #[test]
    fn test_simplicial_complex_add_triangle() {
        let mut sc = SimplicialComplex::new();
        sc.add_simplex(Simplex::new(vec![0, 1, 2], 1.0));
        assert_eq!(sc.count(0), 3); // 3 vertices
        assert_eq!(sc.count(1), 3); // 3 edges
        assert_eq!(sc.count(2), 1); // 1 triangle
    }

    #[test]
    fn test_simplicial_complex_euler_triangle() {
        let mut sc = SimplicialComplex::new();
        sc.add_simplex(Simplex::new(vec![0, 1, 2], 1.0));
        // χ = V - E + F = 3 - 3 + 1 = 1
        assert_eq!(sc.euler_characteristic(), 1);
    }

    #[test]
    fn test_simplicial_complex_no_duplicates() {
        let mut sc = SimplicialComplex::new();
        sc.add_simplex(Simplex::new(vec![0, 1], 1.0));
        sc.add_simplex(Simplex::new(vec![0, 1], 2.0)); // duplicate
        assert_eq!(sc.count(1), 1);
    }

    #[test]
    fn test_simplicial_complex_boundary_matrix_edge() {
        let mut sc = SimplicialComplex::new();
        sc.add_simplex(Simplex::new(vec![0, 1], 1.0));
        sc.add_simplex(Simplex::new(vec![1, 2], 1.0));
        let mat = sc.boundary_matrix(1);
        assert_eq!(mat.len(), 2); // 2 edges
    }

    #[test]
    fn test_simplicial_complex_default() {
        let sc = SimplicialComplex::default();
        assert_eq!(sc.count(0), 0);
    }

    // ── FilteredComplex ───────────────────────────────────────────────────────

    #[test]
    fn test_filtered_complex_ordering() {
        let simplices = vec![
            Simplex::new(vec![0, 1], 2.0),
            Simplex::new(vec![0], 0.0),
            Simplex::new(vec![1], 0.0),
        ];
        let fc = FilteredComplex::new(simplices);
        // Vertices first (dim 0, filt 0) then edge (dim 1, filt 2)
        assert_eq!(fc.ordered_simplices[0].dimension, 0);
        assert_eq!(fc.ordered_simplices[2].dimension, 1);
    }

    #[test]
    fn test_filtered_complex_len() {
        let simplices = vec![
            Simplex::new(vec![0], 0.0),
            Simplex::new(vec![1], 0.0),
            Simplex::new(vec![0, 1], 1.0),
        ];
        let fc = FilteredComplex::new(simplices);
        assert_eq!(fc.len(), 3);
        assert!(!fc.is_empty());
    }

    #[test]
    fn test_filtered_complex_empty() {
        let fc = FilteredComplex::new(vec![]);
        assert_eq!(fc.len(), 0);
        assert!(fc.is_empty());
    }

    #[test]
    fn test_filtered_complex_boundary_column_vertex() {
        let simplices = vec![Simplex::new(vec![0], 0.0)];
        let fc = FilteredComplex::new(simplices);
        assert!(fc.boundary_column(0).is_empty());
    }

    #[test]
    fn test_filtered_complex_reduce_two_points() {
        let simplices = vec![
            Simplex::new(vec![0], 0.0),
            Simplex::new(vec![1], 0.0),
            Simplex::new(vec![0, 1], 1.0),
        ];
        let fc = FilteredComplex::new(simplices);
        let pairs = fc.reduce();
        // The edge should pair with one of the vertices
        assert!(!pairs.is_empty());
    }

    // ── PersistenceDiagram ────────────────────────────────────────────────────

    #[test]
    fn test_persistence_diagram_two_points() {
        let simplices = vec![
            Simplex::new(vec![0], 0.0),
            Simplex::new(vec![1], 0.0),
            Simplex::new(vec![0, 1], 1.5),
        ];
        let fc = FilteredComplex::new(simplices);
        let diag = PersistenceDiagram::from_filtered_complex(&fc);
        // H0: one component dies when edge is added
        assert!(!diag.pairs.is_empty() || !diag.essential.is_empty());
    }

    #[test]
    fn test_persistence_diagram_betti_numbers() {
        let simplices = vec![
            Simplex::new(vec![0], 0.0),
            Simplex::new(vec![1], 0.0),
            Simplex::new(vec![0, 1], 1.0),
        ];
        let fc = FilteredComplex::new(simplices);
        let diag = PersistenceDiagram::from_filtered_complex(&fc);
        // At t=0.5: 2 components, at t=2.0: 1 component
        let b0 = diag.betti_numbers(2.0, 1);
        assert_eq!(b0[0], 1);
    }

    #[test]
    fn test_persistence_diagram_total_persistence() {
        let simplices = vec![
            Simplex::new(vec![0], 0.0),
            Simplex::new(vec![1], 0.0),
            Simplex::new(vec![0, 1], 1.0),
        ];
        let fc = FilteredComplex::new(simplices);
        let diag = PersistenceDiagram::from_filtered_complex(&fc);
        assert!(diag.total_persistence() >= 0.0);
    }

    // ── RipsFiltration ────────────────────────────────────────────────────────

    #[test]
    fn test_rips_from_point_cloud_distances() {
        let pts = [[0.0_f64, 0.0_f64], [1.0, 0.0], [0.0, 1.0]];
        let rips = RipsFiltration::from_point_cloud(&pts, 2.0, 2);
        assert_eq!(rips.n_points, 3);
        assert!((rips.distances[0][1] - 1.0_f64).abs() < 1e-12);
    }

    #[test]
    fn test_rips_build_complex_vertices() {
        let pts = [[0.0_f64, 0.0_f64], [1.0, 0.0], [2.0, 0.0]];
        let rips = RipsFiltration::from_point_cloud(&pts, 5.0, 1);
        let fc = rips.build_complex();
        let n_vertices = fc
            .ordered_simplices
            .iter()
            .filter(|s| s.dimension == 0)
            .count();
        assert_eq!(n_vertices, 3);
    }

    #[test]
    fn test_rips_build_complex_edges_within_radius() {
        let pts = [[0.0_f64, 0.0_f64], [1.0, 0.0], [10.0, 0.0]];
        let rips = RipsFiltration::from_point_cloud(&pts, 1.5, 1);
        let fc = rips.build_complex();
        let n_edges = fc
            .ordered_simplices
            .iter()
            .filter(|s| s.dimension == 1)
            .count();
        // Only edge [0,1] has dist 1.0 ≤ 1.5
        assert_eq!(n_edges, 1);
    }

    #[test]
    fn test_rips_from_distance_matrix() {
        let d = vec![
            vec![0.0_f64, 1.0, 2.0],
            vec![1.0, 0.0, 1.5],
            vec![2.0, 1.5, 0.0],
        ];
        let rips = RipsFiltration::from_distance_matrix(d, 2.0, 2);
        assert_eq!(rips.n_points, 3);
    }

    #[test]
    fn test_rips_persistence_diagram_single_point() {
        let pts = [[0.0_f64, 0.0_f64]];
        let rips = RipsFiltration::from_point_cloud(&pts, 1.0, 0);
        let diag = rips.persistence_diagram();
        // Single point: one essential H0 class
        assert!(!diag.essential.is_empty());
    }

    #[test]
    fn test_rips_triangle_has_triangles() {
        let pts = [[0.0_f64, 0.0_f64], [1.0, 0.0], [0.5, 1.0]];
        let rips = RipsFiltration::from_point_cloud(&pts, 3.0, 2);
        let fc = rips.build_complex();
        let n_tri = fc
            .ordered_simplices
            .iter()
            .filter(|s| s.dimension == 2)
            .count();
        assert_eq!(n_tri, 1);
    }

    // ── PersistenceImage ──────────────────────────────────────────────────────

    #[test]
    fn test_persistence_image_shape() {
        let diag = PersistenceDiagram {
            pairs: vec![(0.0, 1.0, 0), (0.5, 2.0, 0)],
            essential: vec![],
        };
        let img = PersistenceImage::from_diagram(&diag, 10, 10, 0.1);
        assert_eq!(img.pixels.len(), 10);
        assert_eq!(img.pixels[0].len(), 10);
    }

    #[test]
    fn test_persistence_image_non_negative() {
        let diag = PersistenceDiagram {
            pairs: vec![(0.0, 1.0, 0), (0.5, 2.0, 1)],
            essential: vec![],
        };
        let img = PersistenceImage::from_diagram(&diag, 8, 8, 0.2);
        for row in &img.pixels {
            for &v in row {
                assert!(v >= 0.0);
            }
        }
    }

    #[test]
    fn test_persistence_image_empty_diagram() {
        let diag = PersistenceDiagram {
            pairs: vec![],
            essential: vec![],
        };
        let img = PersistenceImage::from_diagram(&diag, 5, 5, 0.1);
        let total: f64 = img.to_vector().iter().sum();
        assert!(total.abs() < 1e-12);
    }

    #[test]
    fn test_persistence_image_to_vector_length() {
        let diag = PersistenceDiagram {
            pairs: vec![(0.0, 1.0, 0)],
            essential: vec![],
        };
        let img = PersistenceImage::from_diagram(&diag, 6, 7, 0.1);
        assert_eq!(img.to_vector().len(), 42);
    }

    #[test]
    fn test_persistence_image_max_value_positive() {
        let diag = PersistenceDiagram {
            pairs: vec![(0.0, 1.0, 0)],
            essential: vec![],
        };
        let img = PersistenceImage::from_diagram(&diag, 5, 5, 0.5);
        assert!(img.max_value() > 0.0);
    }

    // ── BarcodeStatistics ─────────────────────────────────────────────────────

    #[test]
    fn test_barcode_stats_empty() {
        let diag = PersistenceDiagram {
            pairs: vec![],
            essential: vec![],
        };
        let stats = BarcodeStatistics::from_diagram(&diag);
        assert_eq!(stats.n_finite, 0);
        assert!((stats.total_persistence).abs() < 1e-12);
    }

    #[test]
    fn test_barcode_stats_single_bar() {
        let diag = PersistenceDiagram {
            pairs: vec![(0.0, 2.0, 0)],
            essential: vec![],
        };
        let stats = BarcodeStatistics::from_diagram(&diag);
        assert_eq!(stats.n_finite, 1);
        assert!((stats.total_persistence - 2.0_f64).abs() < 1e-12);
        assert!((stats.mean_persistence - 2.0_f64).abs() < 1e-12);
        assert!((stats.max_persistence - 2.0_f64).abs() < 1e-12);
    }

    #[test]
    fn test_barcode_stats_median_odd() {
        let diag = PersistenceDiagram {
            pairs: vec![(0.0, 1.0, 0), (0.0, 2.0, 0), (0.0, 3.0, 0)],
            essential: vec![],
        };
        let stats = BarcodeStatistics::from_diagram(&diag);
        assert!((stats.median_persistence - 2.0_f64).abs() < 1e-12);
    }

    #[test]
    fn test_barcode_stats_median_even() {
        let diag = PersistenceDiagram {
            pairs: vec![(0.0, 1.0, 0), (0.0, 3.0, 0)],
            essential: vec![],
        };
        let stats = BarcodeStatistics::from_diagram(&diag);
        assert!((stats.median_persistence - 2.0_f64).abs() < 1e-12);
    }

    #[test]
    fn test_barcode_stats_essential_counted() {
        let diag = PersistenceDiagram {
            pairs: vec![],
            essential: vec![(0.0, 0), (0.0, 1)],
        };
        let stats = BarcodeStatistics::from_diagram(&diag);
        assert_eq!(stats.n_essential, 2);
    }

    #[test]
    fn test_barcode_bottleneck_identical() {
        let diag = PersistenceDiagram {
            pairs: vec![(0.0, 1.0, 0), (1.0, 3.0, 1)],
            essential: vec![],
        };
        let d = BarcodeStatistics::bottleneck_distance(&diag, &diag);
        assert!(d < 1e-12, "identical diagrams have bottleneck distance 0");
    }

    #[test]
    fn test_barcode_bottleneck_nonneg() {
        let a = PersistenceDiagram {
            pairs: vec![(0.0, 1.0, 0)],
            essential: vec![],
        };
        let b = PersistenceDiagram {
            pairs: vec![(0.1, 1.2, 0)],
            essential: vec![],
        };
        assert!(BarcodeStatistics::bottleneck_distance(&a, &b) >= 0.0);
    }

    #[test]
    fn test_barcode_bottleneck_empty_vs_nonempty() {
        let empty = PersistenceDiagram {
            pairs: vec![],
            essential: vec![],
        };
        let nonempty = PersistenceDiagram {
            pairs: vec![(0.0, 2.0, 0)],
            essential: vec![],
        };
        let d = BarcodeStatistics::bottleneck_distance(&empty, &nonempty);
        // Should equal diag_dist of (0, 2) = 1.0
        assert!((d - 1.0_f64).abs() < 1e-12);
    }

    #[test]
    fn test_xor_sorted_vecs_cancel() {
        let a = vec![0, 1, 2];
        let b = vec![1, 2, 3];
        let result = xor_sorted_vecs(a, b);
        assert_eq!(result, vec![0, 3]);
    }

    #[test]
    fn test_xor_sorted_vecs_empty() {
        let result = xor_sorted_vecs(vec![], vec![1, 2]);
        assert_eq!(result, vec![1, 2]);
    }
}
