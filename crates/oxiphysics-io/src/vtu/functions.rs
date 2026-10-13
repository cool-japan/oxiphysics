//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::{VtuAnnotation, VtuCellQuality, VtuWriter};
use crate::small_sym_eigen::symmetric_eigen;

/// VTK cell type: vertex (1 point).
pub const VTK_VERTEX: u8 = 1;
/// VTK cell type: poly vertex (variable points).
pub const VTK_POLY_VERTEX: u8 = 2;
/// VTK cell type: line (2 points).
pub const VTK_LINE: u8 = 3;
/// VTK cell type: triangle (3 points).
pub const VTK_TRIANGLE: u8 = 5;
/// VTK cell type: triangle strip (variable points).
pub const VTK_TRIANGLE_STRIP: u8 = 6;
/// VTK cell type: polygon (variable points).
pub const VTK_POLYGON: u8 = 7;
/// VTK cell type: quad (4 points).
pub const VTK_QUAD: u8 = 9;
/// VTK cell type: tetrahedron (4 points).
pub const VTK_TET: u8 = 10;
/// VTK cell type: voxel (8 points, axis-aligned hex).
pub const VTK_VOXEL: u8 = 11;
/// VTK cell type: hexahedron (8 points).
pub const VTK_HEX: u8 = 12;
/// VTK cell type: wedge / pentahedron (6 points).
pub const VTK_WEDGE: u8 = 13;
/// VTK cell type: pyramid (5 points).
pub const VTK_PYRAMID: u8 = 14;
/// Return the expected number of nodes for a given VTK cell type.
/// Returns 0 for variable-size cells (polygon, poly-vertex, etc.).
pub fn cell_type_num_nodes(cell_type: u8) -> usize {
    match cell_type {
        VTK_VERTEX => 1,
        VTK_POLY_VERTEX => 0,
        VTK_LINE => 2,
        VTK_TRIANGLE => 3,
        VTK_TRIANGLE_STRIP => 0,
        VTK_POLYGON => 0,
        VTK_QUAD => 4,
        VTK_TET => 4,
        VTK_VOXEL => 8,
        VTK_HEX => 8,
        VTK_WEDGE => 6,
        VTK_PYRAMID => 5,
        _ => 0,
    }
}
/// VTK cell type: quadratic triangle (6 points).
pub const VTK_QUADRATIC_TRIANGLE: u8 = 22;
/// VTK cell type: quadratic quad (8 points).
pub const VTK_QUADRATIC_QUAD: u8 = 23;
/// VTK cell type: quadratic tetrahedron (10 points).
pub const VTK_QUADRATIC_TET: u8 = 24;
/// VTK cell type: quadratic hexahedron (20 points).
pub const VTK_QUADRATIC_HEX: u8 = 25;
/// Return the number of nodes for common higher-order VTK cell types.
pub fn higher_order_cell_num_nodes(cell_type: u8) -> usize {
    match cell_type {
        VTK_QUADRATIC_TRIANGLE => 6,
        VTK_QUADRATIC_QUAD => 8,
        VTK_QUADRATIC_TET => 10,
        VTK_QUADRATIC_HEX => 20,
        other => cell_type_num_nodes(other),
    }
}
/// Local node pairs of the edges of a VTK line cell.
const LINE_EDGES: [(usize, usize); 1] = [(0, 1)];
/// Local node pairs of the 3 edges of a VTK triangle.
const TRIANGLE_EDGES: [(usize, usize); 3] = [(0, 1), (1, 2), (2, 0)];
/// Local node pairs of the 4 edges of a VTK quadrilateral.
const QUAD_EDGES: [(usize, usize); 4] = [(0, 1), (1, 2), (2, 3), (3, 0)];
/// Local node pairs of the 6 edges of a VTK tetrahedron.
const TET_EDGES: [(usize, usize); 6] = [(0, 1), (1, 2), (2, 0), (0, 3), (1, 3), (2, 3)];
/// Local node pairs of the 8 edges of a VTK pyramid (quadrilateral base
/// `0..4`, apex `4`).
const PYRAMID_EDGES: [(usize, usize); 8] = [
    (0, 1),
    (1, 2),
    (2, 3),
    (3, 0),
    (0, 4),
    (1, 4),
    (2, 4),
    (3, 4),
];
/// Local node pairs of the 9 edges of a VTK wedge (triangles `0..3` and
/// `3..6`, joined by the edges `i -> i + 3`).
const WEDGE_EDGES: [(usize, usize); 9] = [
    (0, 1),
    (1, 2),
    (2, 0),
    (3, 4),
    (4, 5),
    (5, 3),
    (0, 3),
    (1, 4),
    (2, 5),
];
/// Local node pairs of the 12 edges of a VTK hexahedron.
const HEX_EDGES: [(usize, usize); 12] = [
    (0, 1),
    (1, 2),
    (2, 3),
    (3, 0),
    (4, 5),
    (5, 6),
    (6, 7),
    (7, 4),
    (0, 4),
    (1, 5),
    (2, 6),
    (3, 7),
];
/// For every corner of a VTK hexahedron the three neighbouring nodes whose
/// edge vectors form a right-handed triple (`det = +1` for a unit cube).
const HEX_CORNER_NEIGHBOURS: [[usize; 3]; 8] = [
    [1, 3, 4],
    [2, 0, 5],
    [3, 1, 6],
    [0, 2, 7],
    [7, 5, 0],
    [4, 6, 1],
    [5, 7, 2],
    [6, 4, 3],
];
/// Right-handed corner triples of the four base corners of a VTK pyramid; the
/// third neighbour is always the apex.
const PYRAMID_BASE_CORNER_NEIGHBOURS: [[usize; 3]; 4] =
    [[1, 3, 4], [2, 0, 4], [3, 1, 4], [0, 2, 4]];
/// Corner triples of a VTK wedge.  VTK orders the base triangle `(0, 1, 2)` so
/// that its right-hand normal points *away* from the opposite triangle
/// `(3, 4, 5)`; the triples are arranged such that a valid (non-inverted) wedge
/// has positive corner determinants.
const WEDGE_CORNER_NEIGHBOURS: [[usize; 3]; 6] = [
    [2, 1, 3],
    [0, 2, 4],
    [1, 0, 5],
    [4, 5, 0],
    [5, 3, 1],
    [3, 4, 2],
];
/// Permutation that turns the node order of a `VTK_VOXEL` (x varies fastest)
/// into the node order of a `VTK_HEX`.
const VOXEL_TO_HEX: [usize; 8] = [0, 1, 3, 2, 4, 5, 7, 6];

/// The linear cell shapes whose quality can be measured.
#[derive(Clone, Copy)]
enum QualityShape {
    Line,
    Triangle,
    Quad,
    Tet,
    Hex,
    Wedge,
    Pyramid,
}

impl QualityShape {
    /// Shape used to measure a VTK cell type.  Quadratic cells are measured on
    /// their corner nodes, which come first in VTK ordering.
    fn from_cell_type(cell_type: u8) -> Option<Self> {
        match cell_type {
            VTK_LINE => Some(Self::Line),
            VTK_TRIANGLE | VTK_QUADRATIC_TRIANGLE => Some(Self::Triangle),
            VTK_QUAD | VTK_QUADRATIC_QUAD => Some(Self::Quad),
            VTK_TET | VTK_QUADRATIC_TET => Some(Self::Tet),
            VTK_HEX | VTK_VOXEL | VTK_QUADRATIC_HEX => Some(Self::Hex),
            VTK_WEDGE => Some(Self::Wedge),
            VTK_PYRAMID => Some(Self::Pyramid),
            _ => None,
        }
    }

    /// Number of corner nodes of the shape.
    fn corner_count(self) -> usize {
        match self {
            Self::Line => 2,
            Self::Triangle => 3,
            Self::Quad | Self::Tet => 4,
            Self::Pyramid => 5,
            Self::Wedge => 6,
            Self::Hex => 8,
        }
    }

    /// Edge table of the shape.
    fn edges(self) -> &'static [(usize, usize)] {
        match self {
            Self::Line => &LINE_EDGES,
            Self::Triangle => &TRIANGLE_EDGES,
            Self::Quad => &QUAD_EDGES,
            Self::Tet => &TET_EDGES,
            Self::Pyramid => &PYRAMID_EDGES,
            Self::Wedge => &WEDGE_EDGES,
            Self::Hex => &HEX_EDGES,
        }
    }
}

#[inline]
fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[inline]
fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// Jacobian of the corner spanned by the edge vectors `a`, `b`, `c`,
/// normalised by the edge lengths: `det(a, b, c) / (|a| |b| |c|)`.  It is `+1`
/// for three orthogonal edges in right-handed order and `0` for a degenerate
/// corner (a zero-length edge).
fn scaled_corner_jacobian(a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> f64 {
    let (la, lb, lc) = (norm3(a), norm3(b), norm3(c));
    if la == 0.0 || lb == 0.0 || lc == 0.0 {
        return 0.0;
    }
    dot3(a, cross3(b, c)) / (la * lb * lc)
}

/// Minimum of the scaled corner Jacobians over a table of corner triples:
/// `triples[k]` lists the three neighbours of node `k` whose edge vectors from
/// node `k` form a right-handed triple.
fn min_corner_jacobian(p: &[[f64; 3]], triples: &[[usize; 3]]) -> f64 {
    let mut worst = f64::INFINITY;
    for (corner, &[a, b, c]) in triples.iter().enumerate() {
        let origin = p[corner];
        let value =
            scaled_corner_jacobian(sub3(p[a], origin), sub3(p[b], origin), sub3(p[c], origin));
        worst = worst.min(value);
    }
    worst
}

/// Scaled Jacobian of a tetrahedron (Verdict definition):
/// `sqrt(2) * det(p1 - p0, p2 - p0, p3 - p0)` divided by the largest product of
/// the three edge lengths meeting at a vertex.  A regular tetrahedron gives
/// `1`, an inverted one a negative value.
fn tet_scaled_jacobian(p: &[[f64; 3]]) -> f64 {
    let det = dot3(sub3(p[1], p[0]), cross3(sub3(p[2], p[0]), sub3(p[3], p[0])));
    let length = |i: usize, j: usize| norm3(sub3(p[i], p[j]));
    let products = [
        length(0, 1) * length(0, 2) * length(0, 3),
        length(1, 0) * length(1, 2) * length(1, 3),
        length(2, 0) * length(2, 1) * length(2, 3),
        length(3, 0) * length(3, 1) * length(3, 2),
    ];
    let largest = products.iter().copied().fold(0.0_f64, f64::max);
    if largest > 0.0 {
        std::f64::consts::SQRT_2 * det / largest
    } else {
        0.0
    }
}

/// Scaled Jacobian of a planar polygon (triangle or quadrilateral): the
/// minimum over the corners of the sine of the signed corner angle, measured
/// against the element normal (Newell's normal of the vertex loop).  A
/// reflex corner or an inverted element yields a negative value, a degenerate
/// element `0`.
fn polygon_scaled_jacobian(p: &[[f64; 3]]) -> f64 {
    let n = p.len();
    let mut normal = [0.0_f64; 3];
    for i in 0..n {
        let c = cross3(sub3(p[i], p[0]), sub3(p[(i + 1) % n], p[0]));
        for (acc, v) in normal.iter_mut().zip(c) {
            *acc += v;
        }
    }
    let len = norm3(normal);
    if len == 0.0 {
        return 0.0;
    }
    let unit = [normal[0] / len, normal[1] / len, normal[2] / len];
    let mut worst = f64::INFINITY;
    for i in 0..n {
        let next = sub3(p[(i + 1) % n], p[i]);
        let prev = sub3(p[(i + n - 1) % n], p[i]);
        let (ln, lp) = (norm3(next), norm3(prev));
        let value = if ln == 0.0 || lp == 0.0 {
            0.0
        } else {
            dot3(cross3(next, prev), unit) / (ln * lp)
        };
        worst = worst.min(value);
    }
    worst
}

/// Quality record of a cell that cannot be evaluated: unsupported cell type,
/// wrong number of nodes for its type, a connectivity entry that does not
/// refer to an existing point, or a non-finite coordinate.
fn unmeasurable_cell_quality(cell_index: usize, cell_type: u8) -> VtuCellQuality {
    VtuCellQuality {
        cell_index,
        cell_type,
        aspect_ratio: 0.0,
        min_edge: 0.0,
        max_edge: 0.0,
        scaled_jacobian: 0.0,
    }
}

/// Measure one cell; see [`compute_cell_quality`].
fn cell_quality(
    cell_index: usize,
    cell_type: u8,
    node_ids: &[usize],
    points: &[[f64; 3]],
) -> VtuCellQuality {
    let Some(shape) = QualityShape::from_cell_type(cell_type) else {
        return unmeasurable_cell_quality(cell_index, cell_type);
    };
    if node_ids.len() != higher_order_cell_num_nodes(cell_type) {
        return unmeasurable_cell_quality(cell_index, cell_type);
    }
    let mut corners: Vec<[f64; 3]> = Vec::with_capacity(shape.corner_count());
    for &node in &node_ids[..shape.corner_count()] {
        match points.get(node) {
            Some(&position) => corners.push(position),
            None => return unmeasurable_cell_quality(cell_index, cell_type),
        }
    }
    if corners
        .iter()
        .flatten()
        .any(|coordinate| !coordinate.is_finite())
    {
        return unmeasurable_cell_quality(cell_index, cell_type);
    }
    if cell_type == VTK_VOXEL {
        let voxel = corners.clone();
        for (slot, &source) in corners.iter_mut().zip(VOXEL_TO_HEX.iter()) {
            *slot = voxel[source];
        }
    }

    let mut min_edge = f64::INFINITY;
    let mut max_edge = 0.0_f64;
    for &(a, b) in shape.edges() {
        let length = norm3(sub3(corners[a], corners[b]));
        min_edge = min_edge.min(length);
        max_edge = max_edge.max(length);
    }
    let aspect_ratio = if max_edge > 0.0 {
        min_edge / max_edge
    } else {
        0.0
    };

    let scaled_jacobian = match shape {
        QualityShape::Line => {
            if max_edge > 0.0 {
                1.0
            } else {
                0.0
            }
        }
        QualityShape::Triangle => 2.0 / 3.0_f64.sqrt() * polygon_scaled_jacobian(&corners),
        QualityShape::Quad => polygon_scaled_jacobian(&corners),
        QualityShape::Tet => tet_scaled_jacobian(&corners),
        QualityShape::Hex => min_corner_jacobian(&corners, &HEX_CORNER_NEIGHBOURS),
        QualityShape::Wedge => {
            2.0 / 3.0_f64.sqrt() * min_corner_jacobian(&corners, &WEDGE_CORNER_NEIGHBOURS)
        }
        QualityShape::Pyramid => {
            std::f64::consts::SQRT_2
                * min_corner_jacobian(&corners, &PYRAMID_BASE_CORNER_NEIGHBOURS)
        }
    };

    VtuCellQuality {
        cell_index,
        cell_type,
        aspect_ratio,
        min_edge,
        max_edge,
        scaled_jacobian,
    }
}

/// Compute per-cell quality metrics for all cells in a `VtuWriter`.
///
/// The metrics are measured on the *true edges* of every cell shape (1 for a
/// line, 3 triangle, 4 quadrilateral, 6 tetrahedron, 8 pyramid, 9 wedge, 12
/// hexahedron/voxel), not on the cyclic node order:
///
/// * `min_edge` / `max_edge` are the shortest and longest edge,
///   `aspect_ratio = min_edge / max_edge` (`1` for an equilateral cell).
/// * `scaled_jacobian` follows the Verdict definitions: the minimum over the
///   element corners of the Jacobian determinant normalised by the lengths of
///   the edges spanning the corner, rescaled so that the ideal element of each
///   shape scores `1` (unit cube, unit square, regular tetrahedron,
///   equilateral triangle, equilateral prism, pyramid with eight equal
///   edges).  Degenerate elements score `0` and inverted or non-convex ones a
///   negative value.  Quadrilaterals and triangles are measured against their
///   own element normal, so they may sit anywhere in 3-D space.
///
/// VTK node ordering is assumed: tetrahedron and pyramid base `(0, 1, 2)` have
/// their right-hand normal pointing towards the remaining node(s), the wedge
/// base `(0, 1, 2)` points away from `(3, 4, 5)`, and `VTK_VOXEL` is first
/// permuted to hexahedron order.  Quadratic triangles, quadrilaterals,
/// tetrahedra and hexahedra are measured on their corner nodes.
///
/// A line has no Jacobian; it reports `1` when it has a non-zero length.
///
/// A cell that cannot be measured (unsupported type such as a vertex or a
/// polygon, a node count that does not match the type, a connectivity entry
/// outside the point list, or a NaN/infinite coordinate) reports all-zero
/// metrics.
pub fn compute_cell_quality(writer: &VtuWriter) -> Vec<VtuCellQuality> {
    let points = writer.points();
    writer
        .cell_types()
        .iter()
        .zip(writer.cells_per_cell().iter())
        .enumerate()
        .map(|(cell_index, (&cell_type, node_ids))| {
            cell_quality(cell_index, cell_type, node_ids, points)
        })
        .collect()
}
#[inline]
pub(super) fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
#[inline]
pub(super) fn norm3(v: [f64; 3]) -> f64 {
    (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
}
/// Convert a structured grid (Nx × Ny × Nz) to VTU unstructured HEX8 mesh.
///
/// `coords` is a flat array of `(Nx+1) × (Ny+1) × (Nz+1)` node positions
/// in (x, y, z) row-major order (z-fastest).
pub fn structured_to_vtu(nx: usize, ny: usize, nz: usize, coords: &[[f64; 3]]) -> VtuWriter {
    let mut writer = VtuWriter::new();
    writer.add_points(coords);
    let idx =
        |i: usize, j: usize, k: usize| -> usize { i * (ny + 1) * (nz + 1) + j * (nz + 1) + k };
    for i in 0..nx {
        for j in 0..ny {
            for k in 0..nz {
                let n = [
                    idx(i, j, k),
                    idx(i + 1, j, k),
                    idx(i + 1, j + 1, k),
                    idx(i, j + 1, k),
                    idx(i, j, k + 1),
                    idx(i + 1, j, k + 1),
                    idx(i + 1, j + 1, k + 1),
                    idx(i, j + 1, k + 1),
                ];
                writer.add_cell(n.to_vec(), VTK_HEX);
            }
        }
    }
    writer
}
/// Generate a uniform Cartesian structured grid as VTU.
///
/// Domain: \[x0, x1\] × \[y0, y1\] × \[z0, z1\] with `nx`×`ny`×`nz` cells.
pub fn uniform_grid_to_vtu(
    nx: usize,
    ny: usize,
    nz: usize,
    x0: f64,
    x1: f64,
    y0: f64,
    y1: f64,
    z0: f64,
    z1: f64,
) -> VtuWriter {
    let dx = (x1 - x0) / nx as f64;
    let dy = (y1 - y0) / ny as f64;
    let dz = (z1 - z0) / nz as f64;
    let mut coords = Vec::with_capacity((nx + 1) * (ny + 1) * (nz + 1));
    for i in 0..=(nx) {
        for j in 0..=(ny) {
            for k in 0..=(nz) {
                coords.push([x0 + i as f64 * dx, y0 + j as f64 * dy, z0 + k as f64 * dz]);
            }
        }
    }
    structured_to_vtu(nx, ny, nz, &coords)
}
/// Interpolate point (node) scalar data to cell-centre values by averaging.
///
/// For each cell the value is the arithmetic mean of its node values.
/// Returns a vec with one value per cell.
pub fn point_to_cell_data(writer: &VtuWriter, point_data: &[f64]) -> Vec<f64> {
    let cells = writer.cells_per_cell();
    let mut result = Vec::with_capacity(cells.len());
    for node_ids in cells {
        if node_ids.is_empty() {
            result.push(0.0);
            continue;
        }
        let sum: f64 = node_ids
            .iter()
            .map(|&ni| *point_data.get(ni).unwrap_or(&0.0))
            .sum();
        result.push(sum / node_ids.len() as f64);
    }
    result
}
/// Interpolate cell-centre scalar data to node values by averaging over
/// all cells that share each node (unweighted).
///
/// Returns a vec with one value per point.
pub fn cell_to_point_data(writer: &VtuWriter, cell_data: &[f64]) -> Vec<f64> {
    let n_pts = writer.num_points();
    let mut sums = vec![0.0_f64; n_pts];
    let mut counts = vec![0usize; n_pts];
    let cells = writer.cells_per_cell();
    for (ci, node_ids) in cells.iter().enumerate() {
        let cv = *cell_data.get(ci).unwrap_or(&0.0);
        for &ni in node_ids {
            if ni < n_pts {
                sums[ni] += cv;
                counts[ni] += 1;
            }
        }
    }
    sums.iter()
        .zip(counts.iter())
        .map(|(s, &c)| if c > 0 { s / c as f64 } else { 0.0 })
        .collect()
}
/// Estimate the gradient of a scalar field at each node with a weighted
/// least-squares fit over the node's neighbours (all nodes that share a cell
/// with it).
///
/// For node `i`, neighbour `j`, `d_j = x_j - x_i` and `ds_j = s_j - s_i` the
/// gradient `g` minimises `sum_j w_j (g . d_j - ds_j)^2` with `w_j = 1 / |d_j|^2`.
/// The `3 x 3` normal equations `A g = b` (`A = sum w d d^T`,
/// `b = sum w ds d`) are solved with the pseudo-inverse of the symmetric matrix
/// `A` obtained from its eigen-decomposition (eigenvalues below `1e-12` times
/// the largest one are dropped).  Linear fields are therefore reproduced
/// exactly on any mesh, and planar or collinear neighbourhoods yield the
/// minimum-norm gradient (zero component normal to the plane) instead of a
/// singular system.
///
/// Nodes that have no neighbour, no scalar value or only coincident neighbours
/// get a zero gradient.  Cost is `O(n_nodes * k)` for `k` neighbours per node.
///
/// Returns interleaved `[gx, gy, gz, gx, gy, gz, ...]`.
pub fn estimate_gradient(writer: &VtuWriter, scalars: &[f64]) -> Vec<f64> {
    let n_pts = writer.num_points();
    let pts = writer.points();
    let mut adj: Vec<Vec<usize>> = vec![Vec::new(); n_pts];
    for node_ids in writer.cells_per_cell() {
        for &ni in node_ids {
            if ni >= n_pts {
                continue;
            }
            for &nj in node_ids {
                if ni != nj && nj < n_pts {
                    adj[ni].push(nj);
                }
            }
        }
    }
    for neighbours in &mut adj {
        neighbours.sort_unstable();
        neighbours.dedup();
    }
    let mut grad = vec![0.0_f64; n_pts * 3];
    for (i, neighbours) in adj.iter().enumerate() {
        let Some(&si) = scalars.get(i) else {
            continue;
        };
        let pi = pts[i];
        let mut normal = [[0.0_f64; 3]; 3];
        let mut rhs = [0.0_f64; 3];
        for &j in neighbours {
            let Some(&sj) = scalars.get(j) else {
                continue;
            };
            let d = sub3(pts[j], pi);
            let d2 = dot3(d, d);
            if d2 < 1e-30 {
                continue;
            }
            let w = 1.0 / d2;
            let ds = sj - si;
            for (r, row) in normal.iter_mut().enumerate() {
                for (c, entry) in row.iter_mut().enumerate() {
                    *entry += w * d[r] * d[c];
                }
                rhs[r] += w * ds * d[r];
            }
        }
        let g = solve_symmetric_pseudo_inverse(normal, rhs);
        grad[i * 3..i * 3 + 3].copy_from_slice(&g);
    }
    grad
}

/// Minimum-norm least-squares solution of `a x = b` for a symmetric positive
/// semi-definite `3 x 3` matrix `a`: the pseudo-inverse is built from the
/// eigen-decomposition, dropping eigenvalues below `1e-12 * lambda_max`.
fn solve_symmetric_pseudo_inverse(a: [[f64; 3]; 3], b: [f64; 3]) -> [f64; 3] {
    let eigen = symmetric_eigen(a);
    let largest = eigen.values[0];
    let mut x = [0.0_f64; 3];
    if largest.is_nan() || largest <= 0.0 {
        return x;
    }
    let cutoff = 1e-12 * largest;
    for (value, vector) in eigen.values.iter().zip(eigen.vectors.iter()) {
        if *value > cutoff {
            let coefficient = dot3(*vector, b) / value;
            for (xi, vi) in x.iter_mut().zip(vector.iter()) {
                *xi += coefficient * vi;
            }
        }
    }
    x
}
/// Write a VTU file with metadata annotations in the XML header comment.
pub fn write_vtu_with_annotations(writer: &VtuWriter, annotations: &[VtuAnnotation]) -> String {
    let xml = writer.to_xml();
    let mut header = String::from("<!-- VTU Annotations\n");
    for ann in annotations {
        header.push_str(&format!("  {}: {}\n", ann.key, ann.value));
    }
    header.push_str("-->\n");
    format!("{}{}", header, xml)
}
/// Compute the axis-aligned bounding box of all points.
///
/// Returns `(min, max)` or `None` if there are no points.
pub fn compute_bounding_box(points: &[[f64; 3]]) -> Option<([f64; 3], [f64; 3])> {
    if points.is_empty() {
        return None;
    }
    let mut min = points[0];
    let mut max = points[0];
    for &p in &points[1..] {
        for k in 0..3 {
            if p[k] < min[k] {
                min[k] = p[k];
            }
            if p[k] > max[k] {
                max[k] = p[k];
            }
        }
    }
    Some((min, max))
}
/// Compute the diagonal length of a bounding box.
pub fn bbox_diagonal(min: [f64; 3], max: [f64; 3]) -> f64 {
    let dx = max[0] - min[0];
    let dy = max[1] - min[1];
    let dz = max[2] - min[2];
    (dx * dx + dy * dy + dz * dz).sqrt()
}
/// Compute the volume (Lx × Ly × Lz) of a bounding box.
pub fn bbox_volume(min: [f64; 3], max: [f64; 3]) -> f64 {
    (max[0] - min[0]).abs() * (max[1] - min[1]).abs() * (max[2] - min[2]).abs()
}
/// Filter cells whose average scalar value (computed from point data) is
/// above a threshold.  Returns the point indices of qualifying cells.
pub fn threshold_cells(writer: &VtuWriter, point_data: &[f64], threshold: f64) -> Vec<usize> {
    let cell_vals = point_to_cell_data(writer, point_data);
    cell_vals
        .iter()
        .enumerate()
        .filter(|&(_, &v)| v > threshold)
        .map(|(i, _)| i)
        .collect()
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::vtu::PvdWriter;
    use crate::vtu::PvtuWriter;
    use crate::vtu::VtuReader;
    use crate::vtu::VtuTimeSeries;
    use crate::vtu::VtuWriterLegacy;
    use oxiphysics_core::math::Vec3;
    #[test]
    fn test_vtu_write() {
        let path = std::env::temp_dir().join("oxiphy_test.vtu");
        let pts = vec![
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
            Vec3::new(0.0, 0.0, 1.0),
        ];
        let cells = vec![[0, 1, 2, 3]];
        let pvals = vec![100.0, 200.0, 300.0, 400.0];
        let cvals = vec![0.166];
        let pdata = vec![("temperature", pvals.as_slice())];
        let cdata = vec![("volume", cvals.as_slice())];
        VtuWriterLegacy::write(path.to_str().unwrap_or(""), &pts, &cells, &pdata, &cdata).unwrap();
        let content = std::fs::read_to_string(&path).unwrap();
        assert!(content.contains("<VTKFile"));
        assert!(content.contains("UnstructuredGrid"));
        assert!(content.contains("temperature"));
        std::fs::remove_file(&path).ok();
    }
    #[test]
    fn test_vtu_builder_empty_grid_xml() {
        let writer = VtuWriter::new();
        let xml = writer.to_xml();
        assert!(xml.contains("<VTKFile"), "missing VTKFile tag: {xml}");
        assert!(xml.contains("NumberOfPoints=\"0\""));
        assert!(xml.contains("NumberOfCells=\"0\""));
    }
    #[test]
    fn test_vtu_builder_single_point() {
        let mut writer = VtuWriter::new();
        writer.add_point([1.0, 2.0, 3.0]);
        let xml = writer.to_xml();
        assert!(xml.contains("NumberOfPoints=\"1\""));
        assert!(
            xml.contains("1 2 3") || xml.contains("1.0 2.0 3.0"),
            "point coords missing: {xml}"
        );
    }
    #[test]
    fn test_vtu_builder_triangle_cell() {
        let mut writer = VtuWriter::new();
        writer.add_point([0.0, 0.0, 0.0]);
        writer.add_point([1.0, 0.0, 0.0]);
        writer.add_point([0.0, 1.0, 0.0]);
        writer.add_cell(vec![0, 1, 2], VTK_TRIANGLE);
        let xml = writer.to_xml();
        assert!(xml.contains("NumberOfCells=\"1\""));
        assert!(xml.contains('5'), "triangle cell type 5 missing: {xml}");
    }
    #[test]
    fn test_vtu_builder_point_data_included() {
        let mut writer = VtuWriter::new();
        writer.add_point([0.0, 0.0, 0.0]);
        writer.add_point([1.0, 0.0, 0.0]);
        writer.add_point_data_scalar("pressure", &[101.0, 202.0]);
        let xml = writer.to_xml();
        assert!(xml.contains("pressure"), "scalar field name missing: {xml}");
        assert!(xml.contains("101"), "scalar value missing: {xml}");
    }
    #[test]
    fn test_vtu_builder_vector_data_included() {
        let mut writer = VtuWriter::new();
        writer.add_point([0.0, 0.0, 0.0]);
        writer.add_point_data_vector("velocity", &[[1.0, 0.0, 0.0]]);
        let xml = writer.to_xml();
        assert!(xml.contains("velocity"), "vector field name missing: {xml}");
        assert!(xml.contains("NumberOfComponents=\"3\""));
    }
    #[test]
    fn test_vtu_reader_parses_points() {
        let mut writer = VtuWriter::new();
        writer.add_point([1.5, 2.5, 3.5]);
        writer.add_point([4.0, 5.0, 6.0]);
        let xml = writer.to_xml();
        let reader = VtuReader::from_xml(&xml).unwrap();
        assert_eq!(reader.points().len(), 2, "expected 2 points");
        assert!((reader.points()[0][0] - 1.5).abs() < 1e-6);
        assert!((reader.points()[1][2] - 6.0).abs() < 1e-6);
    }
    #[test]
    fn test_pvd_format() {
        let pvd = PvdWriter::write(&[(0.0, "frame_0.vtu"), (0.1, "frame_1.vtu")]);
        assert!(pvd.contains("Collection"), "PVD missing Collection: {pvd}");
        assert!(pvd.contains("frame_0.vtu"));
        assert!(pvd.contains("timestep=\"0\"") || pvd.contains("timestep=\"0.0\""));
    }
    #[test]
    fn test_vtu_write_to_file() {
        let path = std::env::temp_dir().join("oxiphy_builder_test.vtu");
        let mut writer = VtuWriter::new();
        writer.add_point([0.0, 0.0, 0.0]);
        writer.add_point([1.0, 0.0, 0.0]);
        writer.add_point([0.0, 1.0, 0.0]);
        writer.add_cell(vec![0, 1, 2], VTK_TRIANGLE);
        writer.write_to_file(path.to_str().unwrap_or("")).unwrap();
        let content = std::fs::read_to_string(&path).unwrap();
        assert!(content.contains("<VTKFile"));
        std::fs::remove_file(&path).ok();
    }
    #[test]
    fn test_cell_type_num_nodes() {
        assert_eq!(cell_type_num_nodes(VTK_VERTEX), 1);
        assert_eq!(cell_type_num_nodes(VTK_LINE), 2);
        assert_eq!(cell_type_num_nodes(VTK_TRIANGLE), 3);
        assert_eq!(cell_type_num_nodes(VTK_QUAD), 4);
        assert_eq!(cell_type_num_nodes(VTK_TET), 4);
        assert_eq!(cell_type_num_nodes(VTK_HEX), 8);
        assert_eq!(cell_type_num_nodes(VTK_WEDGE), 6);
        assert_eq!(cell_type_num_nodes(VTK_PYRAMID), 5);
        assert_eq!(cell_type_num_nodes(VTK_POLYGON), 0);
        assert_eq!(cell_type_num_nodes(255), 0);
    }
    #[test]
    fn test_add_points_bulk() {
        let mut writer = VtuWriter::new();
        let first = writer.add_points(&[[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]);
        assert_eq!(first, 0);
        assert_eq!(writer.num_points(), 3);
    }
    #[test]
    fn test_add_triangle_shorthand() {
        let mut writer = VtuWriter::new();
        writer.add_points(&[[0.0; 3]; 3]);
        writer.add_triangle(0, 1, 2);
        assert_eq!(writer.num_cells(), 1);
        let xml = writer.to_xml();
        assert!(xml.contains('5'));
    }
    #[test]
    fn test_add_quad_shorthand() {
        let mut writer = VtuWriter::new();
        writer.add_points(&[[0.0; 3]; 4]);
        writer.add_quad(0, 1, 2, 3);
        assert_eq!(writer.num_cells(), 1);
        let xml = writer.to_xml();
        assert!(xml.contains('9'));
    }
    #[test]
    fn test_add_tet_shorthand() {
        let mut writer = VtuWriter::new();
        writer.add_points(&[[0.0; 3]; 4]);
        writer.add_tet(0, 1, 2, 3);
        assert_eq!(writer.num_cells(), 1);
    }
    #[test]
    fn test_add_hex_shorthand() {
        let mut writer = VtuWriter::new();
        writer.add_points(&[[0.0; 3]; 8]);
        writer.add_hex(0, 1, 2, 3, 4, 5, 6, 7);
        assert_eq!(writer.num_cells(), 1);
        let xml = writer.to_xml();
        assert!(xml.contains("12"));
    }
    #[test]
    fn test_cell_data_scalar_in_xml() {
        let mut writer = VtuWriter::new();
        writer.add_points(&[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]);
        writer.add_triangle(0, 1, 2);
        writer.add_cell_data_scalar("volume", &[0.5]);
        let xml = writer.to_xml();
        assert!(
            xml.contains("<CellData>"),
            "CellData section missing: {xml}"
        );
        assert!(xml.contains("volume"), "cell scalar name missing: {xml}");
        assert!(xml.contains("0.5"), "cell scalar value missing: {xml}");
    }
    #[test]
    fn test_cell_data_vector_in_xml() {
        let mut writer = VtuWriter::new();
        writer.add_points(&[[0.0; 3]; 3]);
        writer.add_triangle(0, 1, 2);
        writer.add_cell_data_vector("normal", &[[0.0, 0.0, 1.0]]);
        let xml = writer.to_xml();
        assert!(xml.contains("normal"));
        assert!(xml.contains("<CellData>"));
    }
    #[test]
    fn test_bounding_box_empty() {
        let writer = VtuWriter::new();
        assert!(writer.bounding_box().is_none());
    }
    #[test]
    fn test_bounding_box_multiple_points() {
        let mut writer = VtuWriter::new();
        writer.add_points(&[[1.0, -2.0, 3.0], [-1.0, 5.0, 0.0], [0.0, 0.0, 10.0]]);
        let (min, max) = writer.bounding_box().unwrap();
        assert!((min[0] - (-1.0)).abs() < 1e-12);
        assert!((min[1] - (-2.0)).abs() < 1e-12);
        assert!((min[2] - 0.0).abs() < 1e-12);
        assert!((max[0] - 1.0).abs() < 1e-12);
        assert!((max[1] - 5.0).abs() < 1e-12);
        assert!((max[2] - 10.0).abs() < 1e-12);
    }
    #[test]
    fn test_validate_valid_writer() {
        let mut writer = VtuWriter::new();
        writer.add_points(&[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]);
        writer.add_triangle(0, 1, 2);
        writer.add_point_data_scalar("temp", &[1.0, 2.0, 3.0]);
        writer.add_cell_data_scalar("vol", &[0.5]);
        assert!(writer.validate().is_ok());
    }
    #[test]
    fn test_validate_bad_connectivity() {
        let mut writer = VtuWriter::new();
        writer.add_point([0.0, 0.0, 0.0]);
        writer.add_cell(vec![0, 99], VTK_LINE);
        assert!(writer.validate().is_err());
    }
    #[test]
    fn test_validate_mismatched_point_data() {
        let mut writer = VtuWriter::new();
        writer.add_points(&[[0.0; 3]; 3]);
        writer.add_point_data_scalar("temp", &[1.0, 2.0]);
        assert!(writer.validate().is_err());
    }
    #[test]
    fn test_validate_mismatched_cell_data() {
        let mut writer = VtuWriter::new();
        writer.add_points(&[[0.0; 3]; 3]);
        writer.add_triangle(0, 1, 2);
        writer.add_cell_data_scalar("vol", &[0.5, 1.0]);
        assert!(writer.validate().is_err());
    }
    #[test]
    fn test_reader_extracts_num_cells() {
        let mut writer = VtuWriter::new();
        writer.add_points(&[[0.0; 3]; 4]);
        writer.add_triangle(0, 1, 2);
        writer.add_cell(vec![2, 3], VTK_LINE);
        let xml = writer.to_xml();
        let reader = VtuReader::from_xml(&xml).unwrap();
        assert_eq!(reader.num_cells(), 2);
    }
    #[test]
    fn test_reader_cell_types() {
        let mut writer = VtuWriter::new();
        writer.add_points(&[[0.0; 3]; 4]);
        writer.add_triangle(0, 1, 2);
        writer.add_cell(vec![2, 3], VTK_LINE);
        let xml = writer.to_xml();
        let reader = VtuReader::from_xml(&xml).unwrap();
        assert_eq!(reader.cell_types().len(), 2);
        assert_eq!(reader.cell_types()[0], VTK_TRIANGLE);
        assert_eq!(reader.cell_types()[1], VTK_LINE);
    }
    #[test]
    fn test_reader_point_data_names() {
        let mut writer = VtuWriter::new();
        writer.add_points(&[[0.0; 3]; 2]);
        writer.add_point_data_scalar("pressure", &[1.0, 2.0]);
        writer.add_point_data_scalar("temperature", &[300.0, 310.0]);
        let xml = writer.to_xml();
        let reader = VtuReader::from_xml(&xml).unwrap();
        let names = reader.point_data_names();
        assert!(names.iter().any(|n| n == "pressure"));
        assert!(names.iter().any(|n| n == "temperature"));
    }
    #[test]
    fn test_reader_cell_data_names() {
        let mut writer = VtuWriter::new();
        writer.add_points(&[[0.0; 3]; 3]);
        writer.add_triangle(0, 1, 2);
        writer.add_cell_data_scalar("volume", &[0.5]);
        let xml = writer.to_xml();
        let reader = VtuReader::from_xml(&xml).unwrap();
        let names = reader.cell_data_names();
        assert!(names.iter().any(|n| n == "volume"));
    }
    #[test]
    fn test_reader_bounding_box() {
        let mut writer = VtuWriter::new();
        writer.add_points(&[[-1.0, -2.0, -3.0], [4.0, 5.0, 6.0]]);
        let xml = writer.to_xml();
        let reader = VtuReader::from_xml(&xml).unwrap();
        let (min, max) = reader.bounding_box().unwrap();
        assert!((min[0] - (-1.0)).abs() < 1e-6);
        assert!((max[2] - 6.0).abs() < 1e-6);
    }
    #[test]
    fn test_pvd_parallel() {
        let pvd = PvdWriter::write_parallel(&[
            (0.0, vec!["frame_0_0.vtu", "frame_0_1.vtu"]),
            (0.1, vec!["frame_1_0.vtu", "frame_1_1.vtu"]),
        ]);
        assert!(pvd.contains("frame_0_0.vtu"));
        assert!(pvd.contains("frame_0_1.vtu"));
        assert!(pvd.contains("part=\"0\""));
        assert!(pvd.contains("part=\"1\""));
    }
    #[test]
    fn test_pvtu_writer() {
        let pvtu = PvtuWriter::write(
            &["piece_0.vtu", "piece_1.vtu"],
            &["temperature"],
            &["volume"],
        );
        assert!(pvtu.contains("PUnstructuredGrid"));
        assert!(pvtu.contains("piece_0.vtu"));
        assert!(pvtu.contains("piece_1.vtu"));
        assert!(pvtu.contains("temperature"));
        assert!(pvtu.contains("volume"));
        assert!(pvtu.contains("PPoints"));
    }
    #[test]
    fn test_roundtrip_write_read_with_cell_data() {
        let mut writer = VtuWriter::new();
        writer.add_points(&[
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.5, 0.5, 1.0],
        ]);
        writer.add_tet(0, 1, 2, 3);
        writer.add_point_data_scalar("density", &[1.0, 1.1, 1.2, 1.3]);
        writer.add_cell_data_scalar("quality", &[0.95]);
        assert!(writer.validate().is_ok());
        let xml = writer.to_xml();
        let reader = VtuReader::from_xml(&xml).unwrap();
        assert_eq!(reader.points().len(), 4);
        assert_eq!(reader.num_cells(), 1);
        assert!(reader.point_data_names().iter().any(|n| n == "density"));
        assert!(reader.cell_data_names().iter().any(|n| n == "quality"));
    }
    #[test]
    fn test_default_trait() {
        let w = VtuWriter::default();
        assert_eq!(w.num_points(), 0);
        assert_eq!(w.num_cells(), 0);
    }
    #[test]
    fn test_multiple_cells_mixed_types() {
        let mut writer = VtuWriter::new();
        writer.add_points(&[[0.0; 3]; 10]);
        writer.add_triangle(0, 1, 2);
        writer.add_quad(3, 4, 5, 6);
        writer.add_cell(vec![7, 8], VTK_LINE);
        writer.add_cell(vec![9], VTK_VERTEX);
        assert_eq!(writer.num_cells(), 4);
        let xml = writer.to_xml();
        assert!(xml.contains("NumberOfCells=\"4\""));
    }
    #[test]
    fn test_time_series_num_frames() {
        let ts = VtuTimeSeries::new("/tmp", "test_ts");
        assert_eq!(ts.num_frames(), 0);
        assert!(ts.times().is_empty());
    }
}
