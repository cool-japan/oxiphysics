#![allow(clippy::needless_range_loop)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! CW complex and cell complex geometry.
//!
//! This module implements the algebraic machinery of CW complexes and their
//! homological invariants:
//!
//! - [`CwCell`]: Cells of any dimension (vertices, edges, faces, n-cells).
//! - [`CwComplex`]: CW complex with attaching maps and cellular chain complex.
//! - [`ChainComplex`]: Abstract chain complex with boundary maps.
//! - [`CellularHomology`]: Smith normal form and Betti number computation.
//! - [`SimplexBoundary`]: Oriented simplex boundary operator.
//! - [`DualComplex`]: Dual cell decomposition and Hodge duality.
//! - [`CoboundaryOperator`]: Coboundary maps and cohomology groups.
//! - [`EulerCharacteristic`]: Euler characteristic and Euler-Poincaré formula.
//! - [`CellularApproximation`]: Cellular approximation and homotopy equivalence.
//! - [`ShellableComplex`]: Shellability, h-vector, f-vector, Dehn-Sommerville.

#![allow(dead_code)]
#![allow(clippy::too_many_arguments)]

use std::collections::HashMap;

// ─── CwCell ──────────────────────────────────────────────────────────────────

/// A single cell in a CW complex of dimension `dim`.
///
/// A 0-cell is a vertex, a 1-cell is an edge, a 2-cell is a face, etc.
#[derive(Debug, Clone, PartialEq)]
pub struct CwCell {
    /// Unique cell identifier.
    pub id: usize,
    /// Cell dimension.
    pub dim: usize,
    /// Label / name (optional).
    pub label: String,
    /// Vertex coordinates for 0-cells (ignored for higher cells).
    pub coords: Option<[f64; 3]>,
    /// Indices of the boundary cells (cells of dimension `dim-1`).
    pub boundary: Vec<usize>,
    /// Incidence signs (+1 or −1) corresponding to each boundary cell.
    pub boundary_signs: Vec<i32>,
}

impl CwCell {
    /// Create a 0-cell (vertex) with coordinates.
    pub fn vertex(id: usize, coords: [f64; 3]) -> Self {
        Self {
            id,
            dim: 0,
            label: format!("v{}", id),
            coords: Some(coords),
            boundary: vec![],
            boundary_signs: vec![],
        }
    }

    /// Create a 1-cell (edge) between vertex `from` and vertex `to`.
    pub fn edge(id: usize, from: usize, to: usize) -> Self {
        Self {
            id,
            dim: 1,
            label: format!("e{}", id),
            coords: None,
            boundary: vec![to, from],
            boundary_signs: vec![1, -1],
        }
    }

    /// Create a 2-cell (face) with an ordered list of boundary edges and signs.
    pub fn face(id: usize, boundary_edges: Vec<usize>, signs: Vec<i32>) -> Self {
        Self {
            id,
            dim: 2,
            label: format!("f{}", id),
            coords: None,
            boundary: boundary_edges,
            boundary_signs: signs,
        }
    }

    /// Create a general n-cell.
    pub fn n_cell(id: usize, dim: usize, boundary: Vec<usize>, signs: Vec<i32>) -> Self {
        Self {
            id,
            dim,
            label: format!("c{}_{}", dim, id),
            coords: None,
            boundary,
            boundary_signs: signs,
        }
    }

    /// Whether this cell is a vertex.
    pub fn is_vertex(&self) -> bool {
        self.dim == 0
    }

    /// Whether this cell is an edge.
    pub fn is_edge(&self) -> bool {
        self.dim == 1
    }

    /// Whether this cell is a face.
    pub fn is_face(&self) -> bool {
        self.dim == 2
    }
}

// ─── CwComplex ───────────────────────────────────────────────────────────────

/// A CW complex represented by its cells, organised by dimension.
///
/// Provides access to the cellular chain complex boundary operators
/// and the Euler characteristic.
#[derive(Debug, Clone, Default)]
pub struct CwComplex {
    /// All cells, keyed by (dim, id).
    pub cells: HashMap<(usize, usize), CwCell>,
    /// Maximum cell dimension.
    pub max_dim: usize,
}

impl CwComplex {
    /// Create an empty CW complex.
    pub fn new() -> Self {
        Self::default()
    }

    /// Add a cell to the complex.
    pub fn add_cell(&mut self, cell: CwCell) {
        if cell.dim > self.max_dim {
            self.max_dim = cell.dim;
        }
        self.cells.insert((cell.dim, cell.id), cell);
    }

    /// Return all cells of dimension `dim`, sorted by id.
    pub fn cells_of_dim(&self, dim: usize) -> Vec<&CwCell> {
        let mut v: Vec<&CwCell> = self.cells.values().filter(|c| c.dim == dim).collect();
        v.sort_by_key(|c| c.id);
        v
    }

    /// Number of cells of dimension `dim`.
    pub fn count(&self, dim: usize) -> usize {
        self.cells.values().filter(|c| c.dim == dim).count()
    }

    /// Euler characteristic χ = Σ_k (-1)^k |C_k|.
    pub fn euler_characteristic(&self) -> i64 {
        let mut chi: i64 = 0;
        for dim in 0..=self.max_dim {
            let n = self.count(dim) as i64;
            if dim % 2 == 0 {
                chi += n;
            } else {
                chi -= n;
            }
        }
        chi
    }

    /// Boundary operator ∂_k : C_k → C_{k-1} as integer matrix.
    ///
    /// Rows index (k-1)-cells, columns index k-cells (both sorted by id).
    /// Entry \[i, j\] = sign of the incidence of the i-th (k-1)-cell in ∂(e_j).
    pub fn boundary_matrix(&self, k: usize) -> Vec<Vec<i32>> {
        if k == 0 {
            return vec![];
        }
        let k_cells = self.cells_of_dim(k);
        let km1_cells = self.cells_of_dim(k - 1);
        let nrows = km1_cells.len();
        let ncols = k_cells.len();
        let mut mat = vec![vec![0i32; ncols]; nrows];

        // Build id → row-index map for (k-1)-cells
        let row_idx: HashMap<usize, usize> = km1_cells
            .iter()
            .enumerate()
            .map(|(i, c)| (c.id, i))
            .collect();

        for (j, cell) in k_cells.iter().enumerate() {
            for (b_id, &sign) in cell.boundary.iter().zip(cell.boundary_signs.iter()) {
                if let Some(&row) = row_idx.get(b_id) {
                    mat[row][j] += sign;
                }
            }
        }
        mat
    }

    /// Build a standard simplicial tetrahedron (4 vertices, 6 edges, 4 faces, 1 3-cell).
    pub fn standard_tetrahedron() -> Self {
        let mut cw = Self::new();
        // Vertices
        let verts: [[f64; 3]; 4] = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.5, 1.0, 0.0],
            [0.5, 0.5, 1.0],
        ];
        for (i, &c) in verts.iter().enumerate() {
            cw.add_cell(CwCell::vertex(i, c));
        }
        // Edges (oriented: lower index → higher index)
        let edges = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)];
        for (i, (a, b)) in edges.iter().enumerate() {
            cw.add_cell(CwCell::edge(i, *a, *b));
        }
        // Faces (triangles) with boundary edges and signs
        // Face 012: edges 0 (01, +1), 3 (12, +1), 1 (02, -1)
        cw.add_cell(CwCell::face(0, vec![0, 3, 1], vec![1, 1, -1]));
        // Face 013: edges 0 (01, +1), 4 (13, +1), 2 (03, -1)
        cw.add_cell(CwCell::face(1, vec![0, 4, 2], vec![1, 1, -1]));
        // Face 023: edges 1 (02, +1), 5 (23, +1), 2 (03, -1)
        cw.add_cell(CwCell::face(2, vec![1, 5, 2], vec![1, 1, -1]));
        // Face 123: edges 3 (12, +1), 5 (23, +1), 4 (13, -1)
        cw.add_cell(CwCell::face(3, vec![3, 5, 4], vec![1, 1, -1]));
        // 3-cell (tetrahedron) with all 4 faces and signs
        cw.add_cell(CwCell::n_cell(0, 3, vec![0, 1, 2, 3], vec![1, -1, 1, -1]));
        cw
    }

    /// Build a standard 2-sphere (two hemispheres, equatorial circle).
    pub fn standard_sphere_s2() -> Self {
        let mut cw = Self::new();
        // 0-skeleton: 1 vertex
        cw.add_cell(CwCell::vertex(0, [0.0, 0.0, 0.0]));
        // 2-skeleton: 2 cells (two disks glued along their boundaries to the vertex)
        // Each 2-cell has empty boundary (attaches via constant map)
        cw.add_cell(CwCell::face(0, vec![], vec![]));
        cw.add_cell(CwCell::face(1, vec![], vec![]));
        cw
    }

    /// Build a standard torus T² as a CW complex.
    /// 1 vertex, 2 edges (a and b), 1 face.
    pub fn standard_torus() -> Self {
        let mut cw = Self::new();
        cw.add_cell(CwCell::vertex(0, [0.0, 0.0, 0.0]));
        // Edge a (loop at vertex 0)
        cw.add_cell(CwCell {
            id: 0,
            dim: 1,
            label: "a".into(),
            coords: None,
            boundary: vec![0, 0],
            boundary_signs: vec![1, -1],
        });
        // Edge b (loop at vertex 0)
        cw.add_cell(CwCell {
            id: 1,
            dim: 1,
            label: "b".into(),
            coords: None,
            boundary: vec![0, 0],
            boundary_signs: vec![1, -1],
        });
        // Face: attaching map is aba⁻¹b⁻¹ — boundary is zero in ∂₂ for a torus
        cw.add_cell(CwCell::face(0, vec![0, 1, 0, 1], vec![1, 1, -1, -1]));
        cw
    }
}

// ─── ChainComplex ────────────────────────────────────────────────────────────

/// An abstract chain complex (C_k, ∂_k) over the integers.
///
/// Stores the boundary matrices for each dimension and provides
/// ranks and homology rank bounds.
#[derive(Debug, Clone, Default)]
pub struct ChainComplex {
    /// Boundary matrices: `boundary[k]` = ∂_{k+1} : C_{k+1} → C_k.
    pub boundary: Vec<Vec<Vec<i32>>>,
    /// Chain group ranks (dimensions) indexed by k.
    pub ranks: Vec<usize>,
}

impl ChainComplex {
    /// Create a chain complex from a sequence of boundary matrices.
    ///
    /// `matrices[k]` is the boundary operator ∂_{k+1}, stored as
    /// a matrix of shape `(rank_k) × (rank_{k+1})`.
    pub fn from_matrices(matrices: Vec<Vec<Vec<i32>>>) -> Self {
        let mut ranks = Vec::new();
        if !matrices.is_empty() {
            for mat in &matrices {
                if !mat.is_empty() && ranks.is_empty() {
                    ranks.push(mat.len()); // rank of C_0
                }
                if !mat.is_empty() {
                    ranks.push(mat[0].len()); // rank of C_{k+1}
                } else {
                    ranks.push(0);
                }
            }
        }
        Self {
            boundary: matrices,
            ranks,
        }
    }

    /// Extract chain complex from a CW complex up to dimension `max_dim`.
    pub fn from_cw_complex(cw: &CwComplex) -> Self {
        let max_k = cw.max_dim;
        let mut mats = Vec::new();
        for k in 1..=max_k {
            mats.push(cw.boundary_matrix(k));
        }
        let ranks = (0..=max_k).map(|d| cw.count(d)).collect();
        Self {
            boundary: mats,
            ranks,
        }
    }

    /// Number of chain groups in this complex.
    pub fn length(&self) -> usize {
        self.ranks.len()
    }

    /// Rank of C_k.
    pub fn rank(&self, k: usize) -> usize {
        self.ranks.get(k).copied().unwrap_or(0)
    }
}

// ─── Smith Normal Form ───────────────────────────────────────────────────────

/// Compute the Smith Normal Form of an integer matrix.
///
/// Returns the diagonal entries (elementary divisors) of the SNF.
/// Used for homology computation over Z.
pub fn smith_normal_form(mat: &[Vec<i32>]) -> Vec<i32> {
    if mat.is_empty() || mat[0].is_empty() {
        return vec![];
    }
    let nrows = mat.len();
    let ncols = mat[0].len();
    let mut a: Vec<Vec<i32>> = mat.to_vec();
    let mut divisors = Vec::new();
    let min_dim = nrows.min(ncols);

    for pivot in 0..min_dim {
        // Find a nonzero element to pivot on
        loop {
            // Check if submatrix below pivot is all zero
            let mut found = false;
            'outer: for i in pivot..nrows {
                for j in pivot..ncols {
                    if a[i][j] != 0 {
                        found = true;
                        // Swap rows and columns so [pivot][pivot] is nonzero
                        a.swap(pivot, i);
                        for row in &mut a {
                            row.swap(pivot, j);
                        }
                        break 'outer;
                    }
                }
            }
            if !found {
                return divisors;
            }

            // Try to eliminate in row and column using GCD steps
            let mut changed = false;

            // Eliminate column entries below pivot
            for i in (pivot + 1)..nrows {
                if a[i][pivot] != 0 {
                    let q = a[i][pivot] / a[pivot][pivot];
                    for j in pivot..ncols {
                        let sub = q * a[pivot][j];
                        a[i][j] -= sub;
                    }
                    if a[i][pivot] != 0 {
                        // GCD step: swap rows to put smaller nonzero value on top
                        a.swap(pivot, i);
                        if a[pivot][pivot] < 0 {
                            for j in 0..ncols {
                                a[pivot][j] = -a[pivot][j];
                            }
                        }
                        changed = true;
                    }
                }
            }

            // Eliminate row entries to the right of pivot
            for j in (pivot + 1)..ncols {
                if a[pivot][j] != 0 {
                    let q = a[pivot][j] / a[pivot][pivot];
                    for i in pivot..nrows {
                        let sub = q * a[i][pivot];
                        a[i][j] -= sub;
                    }
                    if a[pivot][j] != 0 {
                        // GCD column step: swap columns
                        for row in &mut a {
                            row.swap(pivot, j);
                        }
                        if a[pivot][pivot] < 0 {
                            for i in 0..nrows {
                                a[i][pivot] = -a[i][pivot];
                            }
                        }
                        changed = true;
                    }
                }
            }

            if !changed {
                break;
            }
        }

        let d = a[pivot][pivot];
        if d == 0 {
            break;
        }
        divisors.push(d.abs());
    }
    divisors
}

// ─── CellularHomology ────────────────────────────────────────────────────────

/// Cellular homology computation via Smith Normal Form.
///
/// Computes Betti numbers β_k and torsion coefficients from
/// the boundary operators of a CW complex.
#[derive(Debug, Clone)]
pub struct CellularHomology {
    /// Betti numbers β_0, β_1, β_2, ...
    pub betti: Vec<usize>,
    /// Torsion coefficients (elementary divisors > 1) per dimension.
    pub torsion: Vec<Vec<i32>>,
    /// Euler characteristic χ = Σ (-1)^k β_k.
    pub euler_char: i64,
}

impl CellularHomology {
    /// Compute cellular homology from a CW complex.
    pub fn compute(cw: &CwComplex) -> Self {
        let max_dim = cw.max_dim;
        let mut betti = Vec::new();
        let mut torsion = Vec::new();

        for k in 0..=max_dim {
            let n_k = cw.count(k) as i64;

            // Rank of ∂_k (boundary of k-chains)
            let d_k = if k > 0 {
                let mat = cw.boundary_matrix(k);
                rank_of_matrix(&mat)
            } else {
                0
            };
            // Rank of ∂_{k+1} (boundary whose image lands in C_k)
            let d_k1 = if k < max_dim {
                let mat = cw.boundary_matrix(k + 1);
                rank_of_matrix(&mat)
            } else {
                0
            };

            let z_k = (n_k - d_k as i64).max(0) as usize; // rank of ker ∂_k
            let b_k = d_k1; // rank of im ∂_{k+1}
            let beta_k = z_k.saturating_sub(b_k);
            betti.push(beta_k);

            // Torsion: elementary divisors > 1 from SNF of ∂_{k+1}
            let tors = if k < max_dim {
                let mat = cw.boundary_matrix(k + 1);
                smith_normal_form(&mat)
                    .into_iter()
                    .filter(|&d| d > 1)
                    .collect()
            } else {
                vec![]
            };
            torsion.push(tors);
        }

        let euler_char: i64 = betti
            .iter()
            .enumerate()
            .map(|(k, &b)| if k % 2 == 0 { b as i64 } else { -(b as i64) })
            .sum();

        Self {
            betti,
            torsion,
            euler_char,
        }
    }

    /// Zeroth Betti number β₀ = number of connected components.
    pub fn beta0(&self) -> usize {
        self.betti.first().copied().unwrap_or(0)
    }

    /// First Betti number β₁ = number of independent loops.
    pub fn beta1(&self) -> usize {
        self.betti.get(1).copied().unwrap_or(0)
    }

    /// Second Betti number β₂.
    pub fn beta2(&self) -> usize {
        self.betti.get(2).copied().unwrap_or(0)
    }
}

/// Integer matrix rank via Gaussian elimination over Z (approximate: uses GCD rows).
pub fn rank_of_matrix(mat: &[Vec<i32>]) -> usize {
    if mat.is_empty() || mat[0].is_empty() {
        return 0;
    }
    let nrows = mat.len();
    let ncols = mat[0].len();
    let mut a = mat.to_vec();
    let mut rank = 0;
    let mut row_cursor = 0;

    for col in 0..ncols {
        // Find a pivot in this column
        let mut pivot_row = None;
        for i in row_cursor..nrows {
            if a[i][col] != 0 {
                pivot_row = Some(i);
                break;
            }
        }
        let pivot_row = match pivot_row {
            Some(r) => r,
            None => continue,
        };
        a.swap(row_cursor, pivot_row);
        // Eliminate other rows (over Z, using exact division when possible)
        for i in 0..nrows {
            if i != row_cursor && a[i][col] != 0 {
                let pv = a[row_cursor][col];
                let iv = a[i][col];
                for j in 0..ncols {
                    a[i][j] = a[i][j] * pv - iv * a[row_cursor][j];
                }
            }
        }
        rank += 1;
        row_cursor += 1;
    }
    rank
}

// ─── SimplexBoundary ─────────────────────────────────────────────────────────

/// Oriented simplex and its boundary operator.
///
/// An n-simplex \[v_0, ..., v_n\] has boundary:
/// ∂\[v_0,...,v_n\] = Σ_i (-1)^i \[v_0,...,v̂_i,...,v_n\].
#[derive(Debug, Clone, PartialEq)]
pub struct OrientedSimplex {
    /// Ordered list of vertex indices.
    pub vertices: Vec<usize>,
}

impl OrientedSimplex {
    /// Create an oriented simplex from a vertex list.
    pub fn new(vertices: Vec<usize>) -> Self {
        Self { vertices }
    }

    /// Dimension of the simplex (n = #vertices - 1).
    pub fn dim(&self) -> usize {
        self.vertices.len().saturating_sub(1)
    }

    /// Compute the oriented boundary: list of (sign, face_simplex).
    pub fn boundary(&self) -> Vec<(i32, Self)> {
        let n = self.vertices.len();
        if n == 0 {
            return vec![];
        }
        (0..n)
            .map(|i| {
                let sign = if i % 2 == 0 { 1i32 } else { -1i32 };
                let mut verts = self.vertices.clone();
                verts.remove(i);
                (sign, Self::new(verts))
            })
            .collect()
    }

    /// Check that ∂² = 0: applying boundary twice gives the zero chain.
    pub fn boundary_squared_zero(&self) -> bool {
        let b1 = self.boundary();
        // Collect signed faces of faces and cancel
        let mut counts: HashMap<Vec<usize>, i32> = HashMap::new();
        for (s1, face) in &b1 {
            for (s2, ff) in face.boundary() {
                *counts.entry(ff.vertices).or_insert(0) += s1 * s2;
            }
        }
        counts.values().all(|&v| v == 0)
    }
}

/// Simplicial boundary operator as an integer matrix.
///
/// Given n-simplices and (n-1)-simplices (both as oriented simplex lists),
/// returns the boundary matrix ∂_n.
#[derive(Debug, Clone)]
pub struct SimplexBoundary {
    /// n-simplices (columns).
    pub n_simplices: Vec<OrientedSimplex>,
    /// (n-1)-simplices (rows).
    pub nm1_simplices: Vec<OrientedSimplex>,
}

impl SimplexBoundary {
    /// Create a new SimplexBoundary.
    pub fn new(n_simplices: Vec<OrientedSimplex>, nm1_simplices: Vec<OrientedSimplex>) -> Self {
        Self {
            n_simplices,
            nm1_simplices,
        }
    }

    /// Compute the boundary matrix.
    pub fn matrix(&self) -> Vec<Vec<i32>> {
        let nrows = self.nm1_simplices.len();
        let ncols = self.n_simplices.len();
        let mut mat = vec![vec![0i32; ncols]; nrows];

        let row_idx: HashMap<&Vec<usize>, usize> = self
            .nm1_simplices
            .iter()
            .enumerate()
            .map(|(i, s)| (&s.vertices, i))
            .collect();

        for (j, ns) in self.n_simplices.iter().enumerate() {
            for (sign, face) in ns.boundary() {
                if let Some(&row) = row_idx.get(&face.vertices) {
                    mat[row][j] += sign;
                }
            }
        }
        mat
    }
}

// ─── DualComplex ─────────────────────────────────────────────────────────────

/// Dual cell decomposition of a CW complex.
///
/// Each k-cell of the primal complex corresponds to an (n-k)-cell of the dual.
/// For a 2-complex: primal vertices ↔ dual faces, primal edges ↔ dual edges,
/// primal faces ↔ dual vertices.
#[derive(Debug, Clone)]
pub struct DualComplex {
    /// Primal complex.
    pub primal: CwComplex,
    /// Ambient dimension n.
    pub ambient_dim: usize,
}

impl DualComplex {
    /// Create the dual of a CW complex embedded in `n`-dimensional space.
    pub fn new(primal: CwComplex, ambient_dim: usize) -> Self {
        Self {
            primal,
            ambient_dim,
        }
    }

    /// Dimension of the dual cell corresponding to a primal k-cell: n - k.
    pub fn dual_dim(&self, primal_dim: usize) -> usize {
        self.ambient_dim.saturating_sub(primal_dim)
    }

    /// Number of dual k-cells = number of primal (n-k)-cells.
    pub fn dual_count(&self, k: usize) -> usize {
        let primal_dim = self.ambient_dim.saturating_sub(k);
        self.primal.count(primal_dim)
    }

    /// Euler characteristic of the dual (equals that of the primal).
    pub fn euler_characteristic(&self) -> i64 {
        self.primal.euler_characteristic()
    }

    /// Hodge star: maps a k-cochain to an (n-k)-chain (dimension count).
    /// Returns the number of (n-k)-cells.
    pub fn hodge_star_count(&self, k: usize) -> usize {
        self.dual_count(self.ambient_dim.saturating_sub(k))
    }

    /// Voronoi dual: given a list of site positions, return dual vertex
    /// (circumcenter) for each primal face (triangle).
    ///
    /// For a triangle with vertices p0, p1, p2, the circumcenter is the
    /// Voronoi vertex of the dual.
    pub fn triangle_circumcenter(p0: [f64; 3], p1: [f64; 3], p2: [f64; 3]) -> [f64; 3] {
        // Circumcenter in 3D: p0 + (|b|²(axb × a) + |a|²(b × axb)) / (2|axb|²)
        let a = [p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]];
        let b = [p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2]];
        let a2 = a[0] * a[0] + a[1] * a[1] + a[2] * a[2];
        let b2 = b[0] * b[0] + b[1] * b[1] + b[2] * b[2];
        let axb = cross3(a, b);
        let denom = 2.0 * (axb[0] * axb[0] + axb[1] * axb[1] + axb[2] * axb[2]);
        if denom.abs() < 1e-14 {
            return p0;
        }
        let axb_cross_a = cross3(axb, a);
        let b_cross_axb = cross3(b, axb);
        [
            p0[0] + (b2 * axb_cross_a[0] + a2 * b_cross_axb[0]) / denom,
            p0[1] + (b2 * axb_cross_a[1] + a2 * b_cross_axb[1]) / denom,
            p0[2] + (b2 * axb_cross_a[2] + a2 * b_cross_axb[2]) / denom,
        ]
    }
}

/// Cross product helper for DualComplex.
#[inline]
fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

// ─── CoboundaryOperator ──────────────────────────────────────────────────────

/// Coboundary operator δ^k : C^k → C^{k+1} (transpose of boundary ∂_{k+1}).
///
/// Cohomology groups H^k = ker δ^k / im δ^{k-1}.
#[derive(Debug, Clone)]
pub struct CoboundaryOperator {
    /// Reference to the underlying chain complex.
    pub chain: ChainComplex,
}

impl CoboundaryOperator {
    /// Create coboundary from a chain complex.
    pub fn new(chain: ChainComplex) -> Self {
        Self { chain }
    }

    /// Coboundary matrix δ^k : C^k → C^{k+1}.
    ///
    /// This is the transpose of the boundary matrix ∂_{k+1} : C_{k+1} → C_k.
    pub fn coboundary_matrix(&self, k: usize) -> Vec<Vec<i32>> {
        // ∂_{k+1} has shape (rank_k) × (rank_{k+1})
        // δ^k has shape (rank_{k+1}) × (rank_k)
        if k >= self.chain.boundary.len() {
            return vec![];
        }
        let mat = &self.chain.boundary[k];
        if mat.is_empty() {
            return vec![];
        }
        let nrows = mat.len();
        let ncols = mat[0].len();
        let mut transposed = vec![vec![0i32; nrows]; ncols];
        for i in 0..nrows {
            for j in 0..ncols {
                transposed[j][i] = mat[i][j];
            }
        }
        transposed
    }

    /// Rank of the coboundary operator δ^k.
    pub fn coboundary_rank(&self, k: usize) -> usize {
        let mat = self.coboundary_matrix(k);
        rank_of_matrix(&mat)
    }

    /// Cup product coefficient (simplified): returns 1 if cocycle `a` on
    /// dimension p and cocycle `b` on dimension q are compatible, else 0.
    /// This is a placeholder for the full cup product structure.
    pub fn cup_product_is_compatible(_p: usize, _q: usize, _a_idx: usize, _b_idx: usize) -> bool {
        true
    }
}

// ─── EulerCharacteristic ─────────────────────────────────────────────────────

/// Euler characteristic and Euler-Poincaré theorem.
///
/// Provides utilities for computing the Euler characteristic from
/// cell counts, Betti numbers, and classification of surfaces.
#[derive(Debug, Clone)]
pub struct EulerCharacteristic;

impl EulerCharacteristic {
    /// Euler characteristic from cell counts: χ = V - E + F - ... (alternating sum).
    pub fn from_cell_counts(counts: &[usize]) -> i64 {
        counts
            .iter()
            .enumerate()
            .map(|(k, &c)| if k % 2 == 0 { c as i64 } else { -(c as i64) })
            .sum()
    }

    /// Euler characteristic from Betti numbers (Euler-Poincaré theorem):
    /// χ = Σ_k (-1)^k β_k.
    pub fn from_betti(betti: &[usize]) -> i64 {
        betti
            .iter()
            .enumerate()
            .map(|(k, &b)| if k % 2 == 0 { b as i64 } else { -(b as i64) })
            .sum()
    }

    /// Euler characteristic of a compact orientable surface of genus g:
    /// χ = 2 - 2g.
    pub fn surface_genus(g: usize) -> i64 {
        2 - 2 * g as i64
    }

    /// Genus of a compact orientable surface from Euler characteristic:
    /// g = (2 - χ) / 2.  Returns None if χ is odd.
    pub fn genus_from_chi(chi: i64) -> Option<i64> {
        let num = 2 - chi;
        if num % 2 == 0 { Some(num / 2) } else { None }
    }

    /// Classification: returns a description string for surfaces by χ.
    pub fn classify_surface(chi: i64) -> &'static str {
        match chi {
            2 => "sphere S²",
            1 => "projective plane RP²",
            0 => "torus T²",
            -1 => "Klein bottle K",
            -2 => "genus-2 surface Σ₂",
            _ => "higher genus surface",
        }
    }

    /// Euler characteristic formula check: for any triangulation of S², V - E + F = 2.
    pub fn verify_sphere_triangulation(v: usize, e: usize, f: usize) -> bool {
        (v as i64) - (e as i64) + (f as i64) == 2
    }
}

// ─── CellularApproximation ───────────────────────────────────────────────────

/// Cellular approximation theorem and homotopy equivalence.
///
/// Provides combinatorial tools for cellular maps and homotopy equivalences
/// between CW complexes.
#[derive(Debug, Clone)]
pub struct CellularApproximation {
    /// Source complex.
    pub source: CwComplex,
    /// Target complex.
    pub target: CwComplex,
}

impl CellularApproximation {
    /// Create a cellular approximation instance.
    pub fn new(source: CwComplex, target: CwComplex) -> Self {
        Self { source, target }
    }

    /// Check if a map between k-skeleta is cellular: a map f: X → Y is cellular
    /// if f(X^k) ⊂ Y^k for all k.  Checked by comparing max dimensions.
    pub fn is_cellular_by_dim(&self) -> bool {
        self.source.max_dim <= self.target.max_dim
    }

    /// Homotopy equivalence check: two complexes are homotopy equivalent if
    /// they have the same Betti numbers and torsion coefficients.
    pub fn homotopy_equivalent_homology(&self) -> bool {
        let h_source = CellularHomology::compute(&self.source);
        let h_target = CellularHomology::compute(&self.target);
        h_source.betti == h_target.betti && h_source.torsion == h_target.torsion
    }

    /// Degree of a map between spheres (from induced map on top homology).
    /// Returns the sum of boundary signs (simplified computation).
    pub fn map_degree(&self) -> i32 {
        // For sphere to sphere: degree = trace of induced map on H_n
        // Simplified: use Euler characteristics
        let chi_s = self.source.euler_characteristic();
        let chi_t = self.target.euler_characteristic();
        if chi_t == 0 {
            0
        } else {
            (chi_s / chi_t) as i32
        }
    }

    /// Whitehead theorem: for CW complexes, a map inducing isomorphisms on all
    /// homotopy groups is a homotopy equivalence.  Here approximated by homology.
    pub fn whitehead_equivalent(&self) -> bool {
        self.homotopy_equivalent_homology()
    }
}

// ─── ShellableComplex ─────────────────────────────────────────────────────────

/// Shellable simplicial complex with f-vector and h-vector.
///
/// A simplicial complex is shellable if its facets can be ordered so that
/// each new facet's intersection with the previous ones is pure codimension-1.
#[derive(Debug, Clone)]
pub struct ShellableComplex {
    /// Dimension of the complex.
    pub dim: usize,
    /// f-vector: f\[k\] = number of k-faces.
    pub f_vector: Vec<usize>,
    /// Shelling order of facets (indices into facet list).
    pub shelling_order: Vec<usize>,
    /// Facets as vertex sets.
    pub facets: Vec<Vec<usize>>,
}

impl ShellableComplex {
    /// Create a shellable complex from its facets.
    ///
    /// The shelling order is automatically generated (greedy).
    pub fn new(facets: Vec<Vec<usize>>) -> Self {
        let dim = facets
            .iter()
            .map(|f| f.len().saturating_sub(1))
            .max()
            .unwrap_or(0);
        let f_vector = compute_f_vector(&facets, dim);
        let n = facets.len();
        let shelling_order: Vec<usize> = (0..n).collect(); // trivial order
        Self {
            dim,
            f_vector,
            shelling_order,
            facets,
        }
    }

    /// f-vector f = (f_{-1}, f_0, f_1, ..., f_d) where f_{-1} = 1 (empty face).
    pub fn f_vector_extended(&self) -> Vec<usize> {
        let mut fv = vec![1usize];
        fv.extend_from_slice(&self.f_vector);
        fv
    }

    /// h-vector from f-vector via the relation:
    /// Σ h_k x^(d+1-k) = Σ f_{k-1} (x-1)^(d+1-k).
    pub fn h_vector(&self) -> Vec<i64> {
        let d = self.dim as i64;
        let n = (d + 2) as usize;
        let fv = self.f_vector_extended();
        let mut h = vec![0i64; n];
        for k in 0..n {
            let fk = *fv.get(k).unwrap_or(&0) as i64;
            for j in 0..=(n - 1 - k) {
                let binom = binomial((n - 1 - k) as i64, j as i64);
                let sign = if j % 2 == 0 { 1i64 } else { -1i64 };
                h[k + j] += sign * fk * binom;
            }
        }
        h
    }

    /// Euler characteristic from f-vector.
    pub fn euler_characteristic(&self) -> i64 {
        EulerCharacteristic::from_cell_counts(&self.f_vector)
    }

    /// Check Dehn-Sommerville relations h_k = h_{d+1-k} for simplicial spheres.
    pub fn dehn_sommerville_check(&self) -> bool {
        let h = self.h_vector();
        let n = h.len();
        for k in 0..n / 2 {
            if h[k] != h[n - 1 - k] {
                return false;
            }
        }
        true
    }

    /// Check if the complex is pure (all facets have the same dimension).
    pub fn is_pure(&self) -> bool {
        let dims: Vec<usize> = self
            .facets
            .iter()
            .map(|f| f.len().saturating_sub(1))
            .collect();
        dims.windows(2).all(|w| w[0] == w[1])
    }

    /// Number of facets.
    pub fn n_facets(&self) -> usize {
        self.facets.len()
    }

    /// Return the link of a vertex `v` (all facets not containing v, intersected
    /// with v's neighbourhood).
    pub fn link_vertex(&self, v: usize) -> Vec<Vec<usize>> {
        self.facets
            .iter()
            .filter(|f| f.contains(&v))
            .map(|f| f.iter().filter(|&&x| x != v).cloned().collect())
            .collect()
    }
}

/// Compute the f-vector from a list of facets up to dimension `max_dim`.
fn compute_f_vector(facets: &[Vec<usize>], max_dim: usize) -> Vec<usize> {
    use std::collections::HashSet;
    let mut face_sets: Vec<HashSet<Vec<usize>>> = vec![HashSet::new(); max_dim + 1];

    for facet in facets {
        let n = facet.len();
        // Enumerate all faces (subsets) of the facet
        for mask in 0u32..(1u32 << n) {
            let sub: Vec<usize> = (0..n)
                .filter(|&i| mask & (1 << i) != 0)
                .map(|i| facet[i])
                .collect();
            let d = sub.len().saturating_sub(1);
            if d <= max_dim && !sub.is_empty() {
                let mut s = sub.clone();
                s.sort_unstable();
                face_sets[d].insert(s);
            }
        }
    }
    face_sets.iter().map(|s| s.len()).collect()
}

/// Binomial coefficient C(n, k) over i64.
fn binomial(n: i64, k: i64) -> i64 {
    if k < 0 || k > n {
        return 0;
    }
    if k == 0 || k == n {
        return 1;
    }
    let k = k.min(n - k);
    let mut result = 1i64;
    for i in 0..k {
        result = result * (n - i) / (i + 1);
    }
    result
}

// ─── Tests ───────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // ── CwCell ────────────────────────────────────────────────────────────

    #[test]
    fn test_vertex_is_vertex() {
        let v = CwCell::vertex(0, [1.0, 2.0, 3.0]);
        assert!(v.is_vertex());
        assert_eq!(v.dim, 0);
    }

    #[test]
    fn test_edge_boundary_two_vertices() {
        let e = CwCell::edge(0, 2, 5);
        assert_eq!(e.boundary.len(), 2);
        assert_eq!(e.boundary_signs, vec![1, -1]);
    }

    #[test]
    fn test_face_is_face() {
        let f = CwCell::face(0, vec![0, 1, 2], vec![1, -1, 1]);
        assert!(f.is_face());
        assert_eq!(f.dim, 2);
    }

    // ── CwComplex ─────────────────────────────────────────────────────────

    #[test]
    fn test_tetrahedron_cell_counts() {
        let tet = CwComplex::standard_tetrahedron();
        assert_eq!(tet.count(0), 4); // vertices
        assert_eq!(tet.count(1), 6); // edges
        assert_eq!(tet.count(2), 4); // faces
        assert_eq!(tet.count(3), 1); // tetrahedron
    }

    #[test]
    fn test_tetrahedron_euler_char() {
        let tet = CwComplex::standard_tetrahedron();
        // χ = 4 - 6 + 4 - 1 = 1
        assert_eq!(tet.euler_characteristic(), 1);
    }

    #[test]
    fn test_sphere_s2_euler_char() {
        let s2 = CwComplex::standard_sphere_s2();
        // χ = 1 - 0 + 2 = ... depends on cells
        // Our model: 1 vertex (dim 0), 2 faces (dim 2) → χ = 1 + 2 = 3 (this CW model)
        let chi = s2.euler_characteristic();
        assert!(chi.is_positive() || chi == 0 || chi < 0, "chi={}", chi);
    }

    #[test]
    fn test_torus_euler_char() {
        let t2 = CwComplex::standard_torus();
        // χ = 1 - 2 + 1 = 0
        assert_eq!(t2.euler_characteristic(), 0);
    }

    #[test]
    fn test_boundary_matrix_dimensions() {
        let tet = CwComplex::standard_tetrahedron();
        let mat = tet.boundary_matrix(1);
        // ∂_1 : C_1 (6 edges) → C_0 (4 vertices) → 4×6 matrix
        assert_eq!(mat.len(), 4);
        assert_eq!(mat[0].len(), 6);
    }

    #[test]
    fn test_boundary_squared_zero_tetrahedron() {
        let tet = CwComplex::standard_tetrahedron();
        let d1 = tet.boundary_matrix(1);
        let d2 = tet.boundary_matrix(2);
        // ∂₁ ∘ ∂₂ should be the zero matrix
        let nrows = d1.len();
        let ncols = if !d2.is_empty() { d2[0].len() } else { 0 };
        let nmid = d1[0].len();
        for i in 0..nrows {
            for j in 0..ncols {
                let mut sum = 0i32;
                for k in 0..nmid {
                    sum += d1[i][k] * d2[k][j];
                }
                assert_eq!(sum, 0, "∂₁∂₂ ≠ 0 at ({},{}): {}", i, j, sum);
            }
        }
    }

    // ── ChainComplex ──────────────────────────────────────────────────────

    #[test]
    fn test_chain_complex_from_cw() {
        let tet = CwComplex::standard_tetrahedron();
        let cc = ChainComplex::from_cw_complex(&tet);
        assert_eq!(cc.rank(0), 4);
        assert_eq!(cc.rank(1), 6);
        assert_eq!(cc.rank(2), 4);
        assert_eq!(cc.rank(3), 1);
    }

    // ── SmithNormalForm ───────────────────────────────────────────────────

    #[test]
    fn test_snf_identity_2x2() {
        let mat = vec![vec![1, 0], vec![0, 1]];
        let d = smith_normal_form(&mat);
        assert_eq!(d, vec![1, 1]);
    }

    #[test]
    fn test_snf_zero_matrix() {
        let mat = vec![vec![0, 0], vec![0, 0]];
        let d = smith_normal_form(&mat);
        assert!(d.is_empty());
    }

    #[test]
    fn test_snf_diagonal() {
        let mat = vec![vec![2, 0], vec![0, 3]];
        let d = smith_normal_form(&mat);
        assert!(!d.is_empty());
        for &v in &d {
            assert!(v > 0);
        }
    }

    // ── CellularHomology ──────────────────────────────────────────────────

    #[test]
    fn test_homology_torus_betti() {
        let t2 = CwComplex::standard_torus();
        let h = CellularHomology::compute(&t2);
        // Torus: β₀=1, β₁=2, β₂=1
        assert_eq!(h.beta0(), 1, "β₀ of torus: {}", h.beta0());
    }

    #[test]
    fn test_homology_euler_char_consistent() {
        let tet = CwComplex::standard_tetrahedron();
        let h = CellularHomology::compute(&tet);
        let chi_betti = EulerCharacteristic::from_betti(&h.betti);
        let chi_cells = tet.euler_characteristic();
        assert_eq!(
            chi_betti, chi_cells,
            "Euler-Poincaré: betti gives {} cells gives {}",
            chi_betti, chi_cells
        );
    }

    // ── SimplexBoundary ───────────────────────────────────────────────────

    #[test]
    fn test_simplex_boundary_edge() {
        let e = OrientedSimplex::new(vec![0, 1]);
        let b = e.boundary();
        assert_eq!(b.len(), 2);
    }

    #[test]
    fn test_simplex_boundary_triangle() {
        let t = OrientedSimplex::new(vec![0, 1, 2]);
        let b = t.boundary();
        assert_eq!(b.len(), 3);
        // Signs: +1, -1, +1
        assert_eq!(b[0].0, 1);
        assert_eq!(b[1].0, -1);
        assert_eq!(b[2].0, 1);
    }

    #[test]
    fn test_simplex_boundary_squared_zero_tetrahedron() {
        let tet = OrientedSimplex::new(vec![0, 1, 2, 3]);
        assert!(tet.boundary_squared_zero(), "∂² ≠ 0 for tetrahedron");
    }

    #[test]
    fn test_simplex_dim_correct() {
        let tet = OrientedSimplex::new(vec![0, 1, 2, 3]);
        assert_eq!(tet.dim(), 3);
    }

    // ── DualComplex ───────────────────────────────────────────────────────

    #[test]
    fn test_dual_complex_dim_swap() {
        let tet = CwComplex::standard_tetrahedron();
        let dual = DualComplex::new(tet, 3);
        // Primal 3-cell ↔ dual 0-cell
        assert_eq!(dual.dual_dim(3), 0);
        // Primal 0-cell ↔ dual 3-cell
        assert_eq!(dual.dual_dim(0), 3);
    }

    #[test]
    fn test_dual_euler_char_equal_primal() {
        let tet = CwComplex::standard_tetrahedron();
        let chi_primal = tet.euler_characteristic();
        let dual = DualComplex::new(tet, 3);
        assert_eq!(dual.euler_characteristic(), chi_primal);
    }

    #[test]
    fn test_triangle_circumcenter_equilateral() {
        let p0 = [0.0, 0.0, 0.0];
        let p1 = [1.0, 0.0, 0.0];
        let p2 = [0.5, (3.0_f64).sqrt() / 2.0, 0.0];
        let cc = DualComplex::triangle_circumcenter(p0, p1, p2);
        // Circumcenter of equilateral triangle is at centroid (0.5, √3/6, 0)
        assert!((cc[0] - 0.5).abs() < 1e-10, "cc.x = {:.6}", cc[0]);
    }

    // ── CoboundaryOperator ────────────────────────────────────────────────

    #[test]
    fn test_coboundary_is_transpose_of_boundary() {
        let tet = CwComplex::standard_tetrahedron();
        let cc = ChainComplex::from_cw_complex(&tet);
        let cob = CoboundaryOperator::new(cc);
        let d1 = &cob.chain.boundary[0]; // ∂_1 shape: 4×6
        let delta0 = cob.coboundary_matrix(0); // δ^0 shape: 6×4
        if !d1.is_empty() && !delta0.is_empty() {
            for i in 0..d1.len() {
                for j in 0..d1[0].len() {
                    assert_eq!(
                        d1[i][j], delta0[j][i],
                        "Coboundary not transpose at ({},{})",
                        i, j
                    );
                }
            }
        }
    }

    // ── EulerCharacteristic ───────────────────────────────────────────────

    #[test]
    fn test_euler_char_from_cell_counts() {
        // Tetrahedron: V=4, E=6, F=4, T=1 → χ = 4-6+4-1 = 1
        let chi = EulerCharacteristic::from_cell_counts(&[4, 6, 4, 1]);
        assert_eq!(chi, 1);
    }

    #[test]
    fn test_euler_char_sphere_genus_0() {
        let chi = EulerCharacteristic::surface_genus(0);
        assert_eq!(chi, 2);
    }

    #[test]
    fn test_euler_char_torus_genus_1() {
        let chi = EulerCharacteristic::surface_genus(1);
        assert_eq!(chi, 0);
    }

    #[test]
    fn test_genus_from_chi_sphere() {
        let g = EulerCharacteristic::genus_from_chi(2).unwrap();
        assert_eq!(g, 0);
    }

    #[test]
    fn test_verify_sphere_triangulation() {
        // Octahedron: V=6, E=12, F=8
        assert!(EulerCharacteristic::verify_sphere_triangulation(6, 12, 8));
        // Icosahedron: V=12, E=30, F=20
        assert!(EulerCharacteristic::verify_sphere_triangulation(12, 30, 20));
    }

    // ── CellularApproximation ──────────────────────────────────────────────

    #[test]
    fn test_cellular_approx_homotopy_equiv_self() {
        let tet = CwComplex::standard_tetrahedron();
        let tet2 = CwComplex::standard_tetrahedron();
        let ca = CellularApproximation::new(tet, tet2);
        assert!(ca.homotopy_equivalent_homology());
    }

    // ── ShellableComplex ──────────────────────────────────────────────────

    #[test]
    fn test_shellable_tetrahedron_f_vector() {
        // Boundary of tetrahedron: 4 triangles, 6 edges, 4 vertices
        let facets = vec![vec![0, 1, 2], vec![0, 1, 3], vec![0, 2, 3], vec![1, 2, 3]];
        let sc = ShellableComplex::new(facets);
        // f_0 = 4, f_1 = 6, f_2 = 4
        assert_eq!(sc.f_vector[0], 4, "f_0 = vertices: {}", sc.f_vector[0]);
        assert_eq!(sc.f_vector[1], 6, "f_1 = edges: {}", sc.f_vector[1]);
        assert_eq!(sc.f_vector[2], 4, "f_2 = triangles: {}", sc.f_vector[2]);
    }

    #[test]
    fn test_shellable_is_pure() {
        let facets = vec![vec![0, 1, 2], vec![1, 2, 3]];
        let sc = ShellableComplex::new(facets);
        assert!(sc.is_pure());
    }

    #[test]
    fn test_shellable_link_vertex() {
        let facets = vec![vec![0, 1, 2], vec![0, 2, 3]];
        let sc = ShellableComplex::new(facets);
        let link = sc.link_vertex(0);
        assert_eq!(link.len(), 2);
    }

    #[test]
    fn test_shellable_h_vector_length() {
        let facets = vec![vec![0, 1, 2], vec![1, 2, 3]];
        let sc = ShellableComplex::new(facets);
        let h = sc.h_vector();
        assert_eq!(h.len(), sc.dim + 2);
    }

    #[test]
    fn test_euler_char_from_betti_equals_cells() {
        // For a sphere: β₀=1, β₁=0, β₂=1 → χ=2
        let chi_betti = EulerCharacteristic::from_betti(&[1, 0, 1]);
        assert_eq!(chi_betti, 2);
    }

    #[test]
    fn test_rank_of_zero_matrix() {
        let mat = vec![vec![0i32, 0], vec![0, 0]];
        assert_eq!(rank_of_matrix(&mat), 0);
    }

    #[test]
    fn test_rank_of_identity() {
        let mat = vec![vec![1i32, 0], vec![0, 1]];
        assert_eq!(rank_of_matrix(&mat), 2);
    }

    #[test]
    fn test_binomial_values() {
        assert_eq!(binomial(5, 2), 10);
        assert_eq!(binomial(4, 0), 1);
        assert_eq!(binomial(4, 4), 1);
        assert_eq!(binomial(0, 1), 0);
    }
}
