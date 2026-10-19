// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! V-HACD: Voxel-based Hierarchical Approximate Convex Decomposition.
//!
//! The decomposition works on a voxel grid, as the original V-HACD does:
//!
//! 1. The input triangle soup is rasterized into a uniform voxel grid
//!    (SAT-based conservative surface voxelisation followed by an interior
//!    flood fill).
//! 2. The occupied voxels fall into 6-connected components (BFS).
//! 3. A *cluster* is a set of voxels.  Its convex hull is the hull of the
//!    corner points of its voxels, i.e. the hull of the union of the voxel
//!    cubes.  The *concavity* of a cluster is the fraction of its hull volume
//!    that is empty space: `(hull_volume - voxel_volume) / hull_volume`, except
//!    that empty space within one voxel of the cluster does not count.  The
//!    voxels of a convex solid with a tilted face form a staircase, so the hull
//!    of the voxels always leaves a thin empty layer along its faces; counting
//!    it would make every tilted convex solid look concave (a rotated cube
//!    reads 10 to 20 % concave), and no split could remove it.  Features
//!    narrower than two voxels, such as a slot or a pinhole, vanish in the same
//!    way, which is the resolution limit of a voxel method.
//! 4. Clusters are split best-first, the one with the largest concavity first,
//!    until every cluster is below the concavity threshold, no split makes
//!    progress, or the maximum number of parts is reached.  A cluster is cut by
//!    an axis-aligned plane: about eight evenly spaced candidate planes per axis
//!    and the planes at the largest steps of the cross-section area are tried.
//!    A candidate costs the empty hull space that remains in the two children,
//!    as a fraction of the hull volume of the parent, plus a penalty for
//!    unbalanced children; the cheapest plane wins.  The penalty makes a ring
//!    break at a plane through the middle of its hole instead of being peeled
//!    one thin slice at a time, which removes nothing of the hole.
//! 5. Every cluster becomes one [`ConvexPart`]: the vertices and triangles of
//!    its hull, the exact hull volume, the centroid of the hull solid and the
//!    concavity of the cluster.
//!
//! The voxel hull covers the part of the input it stands for, but it is
//! larger than the input by up to one voxel on each side, so the volumes of the
//! parts converge to the volume of the input as the resolution grows.

use oxiphysics_core::math::Vec3;
use std::collections::BTreeSet;

use crate::convex_decomposition::{ConvexPart, ConvexParts};
use crate::quickhull::ConvexHull3DVec;
use crate::voxel_grid::VoxelGrid;

/// Number of evenly spaced split planes tried per axis.
const EVEN_SPLIT_CANDIDATES: usize = 8;
/// Number of split planes per axis taken at the largest steps of the
/// cross-section area.
const PROFILE_STEP_CANDIDATES: usize = 2;
/// Children whose hull volume is not at least this fraction below the volume of
/// the parent hull make no progress.
const MIN_RELATIVE_VOLUME_DROP: f64 = 1e-6;
/// Weight of the imbalance of the children, `|n_a - n_b| / n` in voxels, in the
/// cost of a candidate split.  A balanced cut is preferred unless an unbalanced
/// one removes at least this much of the hull volume per unit of imbalance.
const BALANCE_WEIGHT: f64 = 0.2;
/// Scores closer than this count as equal when choosing a split.
const SCORE_TOLERANCE: f64 = 1e-12;

// ---------------------------------------------------------------------------
// Public types
// ---------------------------------------------------------------------------

/// Configuration for V-HACD voxel decomposition.
#[derive(Debug, Clone)]
pub struct VHacdConfig {
    /// Grid resolution along the longest axis (other axes scale proportionally).
    ///
    /// Higher values give better fidelity at the cost of memory and compute.
    /// `0` is treated as `1`.  Default: `64`.
    pub resolution: usize,
    /// Maximum concavity of a cluster that is accepted as convex.
    ///
    /// The concavity is the fraction of the volume of the cluster's convex
    /// hull that is empty space, a number in `[0, 1)`.  Empty space within one
    /// voxel of the cluster does not count: the voxel grid leaves such a thin
    /// layer along the hull even of a convex solid (see the module
    /// documentation).  A value of `0.05` accepts a cluster when at most 5 % of
    /// its hull is empty space farther than a voxel from the cluster.
    /// Default: `0.05`.
    pub concavity_threshold: f64,
    /// Maximum number of nested splits of a connected component.  Default: `8`.
    pub max_depth: u32,
    /// A cluster with fewer than twice this many voxels is not split any
    /// further and becomes a leaf part.  Default: `32`.
    pub min_voxels_per_cluster: usize,
    /// Maximum number of output parts.  Default: `64`.
    ///
    /// When the voxelised input has more connected components than this, only
    /// the largest ones are kept.
    pub max_parts: usize,
}

impl Default for VHacdConfig {
    fn default() -> Self {
        Self {
            resolution: 64,
            concavity_threshold: 0.05,
            max_depth: 8,
            min_voxels_per_cluster: 32,
            max_parts: 64,
        }
    }
}

/// Error variants for V-HACD decomposition.
#[derive(Debug, thiserror::Error)]
pub enum VHacdError {
    /// The supplied triangle list was empty.
    #[error("input triangle soup is empty")]
    EmptyInput,
}

// ---------------------------------------------------------------------------
// Voxel clusters
// ---------------------------------------------------------------------------

/// Voxel coordinates `[ix, iy, iz]`.
type VoxelCoord = [usize; 3];

/// A set of voxels together with the convex hull of their corner points.
struct Cluster {
    /// Coordinates of the voxels of the cluster.
    voxels: Vec<VoxelCoord>,
    /// Convex hull of the corner points of the voxels.
    hull: ConvexHull3DVec,
    /// Volume of `hull`.
    hull_volume: f64,
    /// Fraction of the hull volume that is empty space more than one voxel away
    /// from the voxels; `0` until [`Cluster::measure`] has been called.
    concavity: f64,
    /// Number of splits that led to this cluster.
    depth: u32,
    /// Set when no split of the cluster makes progress.
    frozen: bool,
}

impl Cluster {
    /// Build the cluster of `voxels` and its hull, `None` if the hull cannot be
    /// built (non-finite grid coordinates).  The concavity is not measured yet.
    fn new(voxels: Vec<VoxelCoord>, depth: u32, grid: &VoxelGrid) -> Option<Self> {
        let points = corner_points(&voxels, grid);
        let hull = ConvexHull3DVec::build(&points)?;
        let hull_volume = hull.volume();
        if !(hull_volume > 0.0 && hull_volume.is_finite()) {
            return None;
        }
        Some(Self {
            voxels,
            hull,
            hull_volume,
            concavity: 0.0,
            depth,
            frozen: false,
        })
    }

    /// Measure the concavity; the candidates of a split are only measured once
    /// the best of them is chosen.
    fn measure(mut self, grid: &VoxelGrid) -> Self {
        let empty = empty_hull_volume(&self.hull, &self.voxels, grid);
        self.concavity = (empty / self.hull_volume).clamp(0.0, 1.0);
        self
    }
}

/// Grow the set of `true` cells of `cells`, a box of `size` cells, by one cell
/// along `axis` in both directions.
fn dilate_along(cells: &mut [bool], size: [usize; 3], axis: usize) {
    let stride = [1, size[0], size[0] * size[1]][axis];
    let length = size[axis];
    let source = cells.to_vec();
    for (i, cell) in cells.iter_mut().enumerate() {
        let coordinate = (i / stride) % length;
        let before = coordinate > 0 && source[i - stride];
        let after = coordinate + 1 < length && source[i + stride];
        *cell = source[i] || before || after;
    }
}

/// Volume of the empty space of the hull of `voxels` that is farther than one
/// voxel from every voxel (not even touching one at an edge or a corner).
///
/// The cells of the box around the cluster are marked near an occupied cell by
/// growing the occupied cells by one cell along each axis in turn.  A cell
/// counts when its centre lies inside the hull and it is not near an occupied
/// cell.  The part of a row of cell centres that lies inside the hull is an
/// interval, found by clipping the row against the face planes of the hull; the
/// rows run along the longest side of the box, which keeps their number small.
fn empty_hull_volume(hull: &ConvexHull3DVec, voxels: &[VoxelCoord], grid: &VoxelGrid) -> f64 {
    let mut low = [usize::MAX; 3];
    let mut high = [0usize; 3];
    for v in voxels {
        for axis in 0..3 {
            low[axis] = low[axis].min(v[axis]);
            high[axis] = high[axis].max(v[axis]);
        }
    }
    if voxels.is_empty() {
        return 0.0;
    }
    let size = [
        high[0] - low[0] + 1,
        high[1] - low[1] + 1,
        high[2] - low[2] + 1,
    ];
    let strides = [1, size[0], size[0] * size[1]];
    let mut near = vec![false; size[0] * size[1] * size[2]];
    for v in voxels {
        let local = [v[0] - low[0], v[1] - low[1], v[2] - low[2]];
        near[local[0] * strides[0] + local[1] * strides[1] + local[2] * strides[2]] = true;
    }
    for axis in 0..3 {
        dilate_along(&mut near, size, axis);
    }

    // Face planes `normal . x <= offset` of the hull, normals pointing out.
    let vertex_count = hull.vertices.len().max(1) as f64;
    let inside = hull.vertices.iter().fold(Vec3::zeros(), |sum, v| sum + v) / vertex_count;
    let planes: Vec<([f64; 3], f64)> = hull
        .faces
        .iter()
        .filter_map(|face| {
            let a = hull.vertices[face[0]];
            let normal = (hull.vertices[face[1]] - a).cross(&(hull.vertices[face[2]] - a));
            let length = normal.norm();
            (length > 0.0).then(|| {
                let sign = if normal.dot(&(inside - a)) > 0.0 {
                    -1.0
                } else {
                    1.0
                };
                let unit = normal * (sign / length);
                ([unit.x, unit.y, unit.z], unit.dot(&a))
            })
        })
        .collect();

    let mut axes = [0usize, 1, 2];
    axes.sort_by_key(|&axis| std::cmp::Reverse(size[axis]));
    let [a, b, c] = axes;
    let centre = |axis: usize, local: usize| {
        grid.origin[axis] + ((low[axis] + local) as f64 + 0.5) * grid.step
    };
    let mut empty_cells = 0usize;
    for local_b in 0..size[b] {
        let centre_b = centre(b, local_b);
        for local_c in 0..size[c] {
            let centre_c = centre(c, local_c);
            // Interval of the hull along the row, in world coordinates.
            let mut lower = f64::NEG_INFINITY;
            let mut upper = f64::INFINITY;
            let mut outside = false;
            for (normal, offset) in &planes {
                let rest = offset - normal[b] * centre_b - normal[c] * centre_c;
                if normal[a] > 1e-12 {
                    upper = upper.min(rest / normal[a]);
                } else if normal[a] < -1e-12 {
                    lower = lower.max(rest / normal[a]);
                } else if rest < 0.0 {
                    outside = true;
                    break;
                }
            }
            if outside {
                continue;
            }
            let to_local = |x: f64| (x - grid.origin[a]) / grid.step - 0.5 - low[a] as f64;
            let first = to_local(lower).ceil().max(0.0) as usize;
            let last = (to_local(upper).floor().min(size[a] as f64 - 1.0)).max(-1.0);
            if last < first as f64 {
                continue;
            }
            let base = local_b * strides[b] + local_c * strides[c];
            empty_cells += (first..=last as usize)
                .filter(|&local_a| !near[base + local_a * strides[a]])
                .count();
        }
    }
    empty_cells as f64 * grid.voxel_volume()
}

/// Axes ordered by decreasing grid size.
fn axes_by_size(dims: [usize; 3]) -> [usize; 3] {
    let mut axes = [0usize, 1, 2];
    axes.sort_by_key(|&a| std::cmp::Reverse(dims[a]));
    axes
}

/// Keep only the first and the last point of every line parallel to `axis`
/// (points that share the other two coordinates): a point strictly between two
/// others on such a line is not a vertex of the convex hull.
fn keep_line_endpoints(points: &mut Vec<[usize; 3]>, axis: usize) {
    let (p, q) = match axis {
        0 => (1, 2),
        1 => (0, 2),
        _ => (0, 1),
    };
    points.sort_unstable_by_key(|pt| (pt[p], pt[q], pt[axis]));
    points.dedup();
    let mut kept: Vec<[usize; 3]> = Vec::with_capacity(points.len());
    let mut start = 0;
    while start < points.len() {
        let mut end = start;
        while end + 1 < points.len()
            && points[end + 1][p] == points[start][p]
            && points[end + 1][q] == points[start][q]
        {
            end += 1;
        }
        kept.push(points[start]);
        if end > start {
            kept.push(points[end]);
        }
        start = end + 1;
    }
    *points = kept;
}

/// Corner points of the voxels that span the convex hull of all their corners.
///
/// The hull of the union of the voxel cubes is the hull of the corners of the
/// surface voxels, and every hull vertex is an end point of the corner set
/// along each of the three axes.  The points are found without enumerating the
/// eight corners of every voxel: first the end points of the corner set along
/// the longest axis are read off the occupied voxel rows, then points that lie
/// between two others along the other two axes are dropped.  The hull of the
/// returned points equals the hull of all corners.
fn corner_points(voxels: &[VoxelCoord], grid: &VoxelGrid) -> Vec<Vec3> {
    let [a, b, c] = axes_by_size(grid.dims);
    let (db, dc) = (grid.dims[b], grid.dims[c]);

    // First and last occupied voxel of every row along axis `a`.
    let mut row_first = vec![usize::MAX; db * dc];
    let mut row_last = vec![0usize; db * dc];
    for v in voxels {
        let row = v[b] * dc + v[c];
        row_first[row] = row_first[row].min(v[a]);
        row_last[row] = row_last[row].max(v[a]);
    }

    // Corner lattice points `(along a, along b, along c)`: the corner set on
    // the lattice line `(lb, lc)` comes from the up to four rows around it.
    let mut points: Vec<[usize; 3]> = Vec::new();
    for lb in 0..=db {
        for lc in 0..=dc {
            let mut first = usize::MAX;
            let mut last = 0usize;
            let mut any = false;
            for ib in [lb.wrapping_sub(1), lb] {
                for ic in [lc.wrapping_sub(1), lc] {
                    if ib < db && ic < dc {
                        let row = ib * dc + ic;
                        if row_first[row] != usize::MAX {
                            any = true;
                            first = first.min(row_first[row]);
                            last = last.max(row_last[row] + 1);
                        }
                    }
                }
            }
            if any {
                points.push([first, lb, lc]);
                points.push([last, lb, lc]);
            }
        }
    }
    keep_line_endpoints(&mut points, 1);
    keep_line_endpoints(&mut points, 2);

    points
        .iter()
        .map(|p| {
            let mut lattice = [0usize; 3];
            lattice[a] = p[0];
            lattice[b] = p[1];
            lattice[c] = p[2];
            Vec3::new(
                grid.origin[0] + lattice[0] as f64 * grid.step,
                grid.origin[1] + lattice[1] as f64 * grid.step,
                grid.origin[2] + lattice[2] as f64 * grid.step,
            )
        })
        .collect()
}

/// Candidate split positions along one axis for a cluster that occupies the
/// voxel slices `lo..=hi`, with `counts[i]` voxels in slice `i`.
///
/// A position `s` cuts between the slices `s - 1` and `s`.  The candidates are
/// evenly spaced positions and the positions at the largest steps of the slice
/// occupancy, where the cross-section of the cluster changes abruptly (the
/// wall of a notch).
fn split_positions(counts: &[usize], lo: usize, hi: usize) -> Vec<usize> {
    if hi <= lo {
        return Vec::new();
    }
    let extent = hi - lo + 1;
    let mut positions: BTreeSet<usize> = BTreeSet::new();
    for k in 1..=EVEN_SPLIT_CANDIDATES {
        let offset = ((extent * k) as f64 / (EVEN_SPLIT_CANDIDATES + 1) as f64).round() as usize;
        positions.insert((lo + offset).clamp(lo + 1, hi));
    }
    let mut steps: Vec<(usize, usize)> = (lo + 1..=hi)
        .map(|i| (counts[i].abs_diff(counts[i - 1]), i))
        .filter(|&(step, _)| step > 0)
        .collect();
    steps.sort_by(|x, y| y.0.cmp(&x.0).then(x.1.cmp(&y.1)));
    for &(_, i) in steps.iter().take(PROFILE_STEP_CANDIDATES) {
        positions.insert(i);
    }
    positions.into_iter().collect()
}

// ---------------------------------------------------------------------------
// VHacdVoxel
// ---------------------------------------------------------------------------

/// Voxel-based approximate convex decomposer.
///
/// Accepts a triangle soup (each element is `[x0,y0,z0, x1,y1,z1, x2,y2,z2]`)
/// and decomposes it into approximately-convex `ConvexPart`s stored in a
/// `ConvexParts` result.
pub struct VHacdVoxel {
    /// Algorithm configuration.
    pub config: VHacdConfig,
}

impl VHacdVoxel {
    /// Create a new decomposer with the given configuration.
    pub fn new(config: VHacdConfig) -> Self {
        Self { config }
    }

    /// Decompose a triangle soup into approximately-convex parts.
    ///
    /// The hulls of the parts cover the voxelised input.  Each part carries the
    /// real volume of its hull, the centroid of the hull solid and the
    /// concavity of its voxel cluster (see the module documentation).
    ///
    /// # Errors
    /// Returns [`VHacdError::EmptyInput`] when `tris` is empty.
    pub fn decompose(&self, tris: &[[f64; 9]]) -> Result<ConvexParts, VHacdError> {
        if tris.is_empty() {
            return Err(VHacdError::EmptyInput);
        }

        // ------------------------------------------------------------------
        // Step 1: bounding box and grid construction
        // ------------------------------------------------------------------
        let (bmin, bmax) = bounding_box_of_tris(tris);
        let extents = [bmax[0] - bmin[0], bmax[1] - bmin[1], bmax[2] - bmin[2]];
        let max_extent = extents.iter().cloned().fold(f64::NEG_INFINITY, f64::max);

        let step = if max_extent > 0.0 {
            max_extent / self.config.resolution.max(1) as f64
        } else {
            1.0
        };

        let dims = [
            ((extents[0] / step).ceil() as usize + 1).max(1),
            ((extents[1] / step).ceil() as usize + 1).max(1),
            ((extents[2] / step).ceil() as usize + 1).max(1),
        ];

        // ------------------------------------------------------------------
        // Step 2: rasterise triangles (surface + interior flood-fill)
        // ------------------------------------------------------------------
        let mut grid = VoxelGrid::new(bmin, step, dims);
        for tri in tris {
            rasterize_triangle_into_grid(tri, &mut grid);
        }
        flood_fill_interior(&mut grid);

        // ------------------------------------------------------------------
        // Step 3: BFS connected components, largest first so that the part
        // limit drops the smallest ones
        // ------------------------------------------------------------------
        let mut components = grid.connected_components();
        components.sort_by_key(|component| std::cmp::Reverse(component.len()));
        let mut clusters: Vec<Cluster> = components
            .iter()
            .filter_map(|component| {
                let voxels: Vec<VoxelCoord> = component
                    .iter()
                    .map(|&idx| {
                        let (ix, iy, iz) = grid_coords(idx, &grid);
                        [ix, iy, iz]
                    })
                    .collect();
                Cluster::new(voxels, 0, &grid)
            })
            .take(self.config.max_parts)
            .map(|cluster| cluster.measure(&grid))
            .collect();

        // ------------------------------------------------------------------
        // Step 4: best-first splitting
        // ------------------------------------------------------------------
        self.split_best_first(&mut clusters, &grid);

        let parts = clusters
            .iter()
            .map(|cluster| hull_to_convex_part(&cluster.hull, cluster.concavity))
            .collect();
        Ok(ConvexParts { parts })
    }

    // -----------------------------------------------------------------------
    // Private helpers
    // -----------------------------------------------------------------------

    /// Whether `cluster` is worth splitting.
    fn can_split(&self, cluster: &Cluster) -> bool {
        !cluster.frozen
            && cluster.concavity > self.config.concavity_threshold
            && cluster.depth < self.config.max_depth
            && cluster.voxels.len() >= 2 * self.config.min_voxels_per_cluster.max(1)
    }

    /// Split the cluster of largest concavity until all are acceptable, no
    /// split makes progress, or the part limit is reached.
    fn split_best_first(&self, clusters: &mut Vec<Cluster>, grid: &VoxelGrid) {
        while clusters.len() < self.config.max_parts {
            let next = clusters
                .iter()
                .enumerate()
                .filter(|(_, cluster)| self.can_split(cluster))
                .max_by(|(_, a), (_, b)| a.concavity.total_cmp(&b.concavity))
                .map(|(index, _)| index);
            let Some(index) = next else { break };
            match self.best_split(&clusters[index], grid) {
                Some((first, second)) => {
                    clusters[index] = first;
                    clusters.push(second);
                }
                None => clusters[index].frozen = true,
            }
        }
    }

    /// Best axis-aligned split of `parent`: the two children that leave the
    /// least empty hull space, with a penalty for unbalanced children (ties:
    /// the smaller total hull volume), or `None` if no candidate plane gives
    /// two children that both have a smaller hull than the parent.  The
    /// children are returned measured.
    fn best_split(&self, parent: &Cluster, grid: &VoxelGrid) -> Option<(Cluster, Cluster)> {
        let mut lo = [usize::MAX; 3];
        let mut hi = [0usize; 3];
        let mut counts: [Vec<usize>; 3] = [
            vec![0; grid.dims[0]],
            vec![0; grid.dims[1]],
            vec![0; grid.dims[2]],
        ];
        for v in &parent.voxels {
            for axis in 0..3 {
                lo[axis] = lo[axis].min(v[axis]);
                hi[axis] = hi[axis].max(v[axis]);
                counts[axis][v[axis]] += 1;
            }
        }

        let largest_child_hull = parent.hull_volume * (1.0 - MIN_RELATIVE_VOLUME_DROP);
        let voxel_count = parent.voxels.len() as f64;
        let parent_voxel_volume = voxel_count * grid.voxel_volume();
        let mut best: Option<(f64, f64, Cluster, Cluster)> = None;
        for axis in 0..3 {
            for position in split_positions(&counts[axis], lo[axis], hi[axis]) {
                let (below, above): (Vec<VoxelCoord>, Vec<VoxelCoord>) = parent
                    .voxels
                    .iter()
                    .copied()
                    .partition(|v| v[axis] < position);
                let (Some(first), Some(second)) = (
                    Cluster::new(below, parent.depth + 1, grid),
                    Cluster::new(above, parent.depth + 1, grid),
                ) else {
                    continue;
                };
                if !(first.hull_volume < largest_child_hull
                    && second.hull_volume < largest_child_hull)
                {
                    continue;
                }
                // The empty hull space left in the children is the sum of the
                // hull volumes minus the volume of the voxels, which the split
                // does not change.
                let total_volume = first.hull_volume + second.hull_volume;
                let imbalance =
                    first.voxels.len().abs_diff(second.voxels.len()) as f64 / voxel_count;
                let score = (total_volume - parent_voxel_volume) / parent.hull_volume
                    + BALANCE_WEIGHT * imbalance;
                let better = match &best {
                    None => true,
                    Some((best_score, best_volume, _, _)) => {
                        score < best_score - SCORE_TOLERANCE
                            || (score <= best_score + SCORE_TOLERANCE
                                && total_volume < *best_volume)
                    }
                };
                if better {
                    best = Some((score, total_volume, first, second));
                }
            }
        }
        best.map(|(_, _, first, second)| (first.measure(grid), second.measure(grid)))
    }
}

// ---------------------------------------------------------------------------
// Standalone helper functions
// ---------------------------------------------------------------------------

/// Decompose linear voxel index into `(ix, iy, iz)`.
#[inline]
fn grid_coords(idx: usize, grid: &VoxelGrid) -> (usize, usize, usize) {
    let iz = idx / (grid.dims[0] * grid.dims[1]);
    let rem = idx % (grid.dims[0] * grid.dims[1]);
    let iy = rem / grid.dims[0];
    let ix = rem % grid.dims[0];
    (ix, iy, iz)
}

/// Compute the bounding box of a triangle soup, with a small padding.
fn bounding_box_of_tris(tris: &[[f64; 9]]) -> ([f64; 3], [f64; 3]) {
    let mut bmin = [f64::INFINITY; 3];
    let mut bmax = [f64::NEG_INFINITY; 3];
    for tri in tris {
        for v in 0..3 {
            for ax in 0..3 {
                let c = tri[v * 3 + ax];
                if c < bmin[ax] {
                    bmin[ax] = c;
                }
                if c > bmax[ax] {
                    bmax[ax] = c;
                }
            }
        }
    }
    let pad = 1e-6;
    for ax in 0..3 {
        bmin[ax] -= pad;
        bmax[ax] += pad;
    }
    (bmin, bmax)
}

/// Conservative SAT-based triangle rasterisation: marks every voxel whose AABB
/// intersects the triangle surface.
fn rasterize_triangle_into_grid(tri: &[f64; 9], grid: &mut VoxelGrid) {
    let v0 = [tri[0], tri[1], tri[2]];
    let v1 = [tri[3], tri[4], tri[5]];
    let v2 = [tri[6], tri[7], tri[8]];

    // Triangle AABB
    let tri_min = [
        v0[0].min(v1[0]).min(v2[0]),
        v0[1].min(v1[1]).min(v2[1]),
        v0[2].min(v1[2]).min(v2[2]),
    ];
    let tri_max = [
        v0[0].max(v1[0]).max(v2[0]),
        v0[1].max(v1[1]).max(v2[1]),
        v0[2].max(v1[2]).max(v2[2]),
    ];

    let ix_min = ((tri_min[0] - grid.origin[0]) / grid.step).floor().max(0.0) as usize;
    let iy_min = ((tri_min[1] - grid.origin[1]) / grid.step).floor().max(0.0) as usize;
    let iz_min = ((tri_min[2] - grid.origin[2]) / grid.step).floor().max(0.0) as usize;
    let ix_max = (((tri_max[0] - grid.origin[0]) / grid.step).ceil() as usize)
        .min(grid.dims[0].saturating_sub(1));
    let iy_max = (((tri_max[1] - grid.origin[1]) / grid.step).ceil() as usize)
        .min(grid.dims[1].saturating_sub(1));
    let iz_max = (((tri_max[2] - grid.origin[2]) / grid.step).ceil() as usize)
        .min(grid.dims[2].saturating_sub(1));

    for iz in iz_min..=iz_max {
        for iy in iy_min..=iy_max {
            for ix in ix_min..=ix_max {
                let vox_min = [
                    grid.origin[0] + ix as f64 * grid.step,
                    grid.origin[1] + iy as f64 * grid.step,
                    grid.origin[2] + iz as f64 * grid.step,
                ];
                let vox_max = [
                    vox_min[0] + grid.step,
                    vox_min[1] + grid.step,
                    vox_min[2] + grid.step,
                ];
                if aabb_triangle_intersect(&vox_min, &vox_max, &v0, &v1, &v2) {
                    grid.mark(ix, iy, iz);
                }
            }
        }
    }
}

/// Interior flood-fill: BFS from all grid-boundary voxels that are **empty**,
/// collecting the exterior region.  Any occupied voxel or one adjacent to the
/// exterior is treated as a surface; everything else is interior and marked.
///
/// Implementation: invert the grid, BFS from corner, mark all reachable empty
/// voxels as "exterior"; remaining unmarked-and-unoccupied voxels are interior
/// and get marked occupied.
fn flood_fill_interior(grid: &mut VoxelGrid) {
    use bitvec::prelude::*;
    use std::collections::VecDeque;

    let [nx, ny, nz] = grid.dims;
    let n = nx * ny * nz;
    // `exterior` tracks empty voxels reachable from the boundary
    let mut exterior: BitVec = bitvec![0; n];
    let mut queue: VecDeque<usize> = VecDeque::new();

    // Seed from all boundary voxels that are currently empty
    let seed_if_empty =
        |idx: usize, occupied: &BitVec, exterior: &mut BitVec, queue: &mut VecDeque<usize>| {
            if !occupied[idx] && !exterior[idx] {
                exterior.set(idx, true);
                queue.push_back(idx);
            }
        };

    // ±z faces
    for iz in [0usize, nz.saturating_sub(1)] {
        for iy in 0..ny {
            for ix in 0..nx {
                let idx = grid.index(ix, iy, iz);
                seed_if_empty(idx, &grid.occupied, &mut exterior, &mut queue);
            }
        }
    }
    // ±y faces
    for iy in [0usize, ny.saturating_sub(1)] {
        for iz in 0..nz {
            for ix in 0..nx {
                let idx = grid.index(ix, iy, iz);
                seed_if_empty(idx, &grid.occupied, &mut exterior, &mut queue);
            }
        }
    }
    // ±x faces
    for ix in [0usize, nx.saturating_sub(1)] {
        for iz in 0..nz {
            for iy in 0..ny {
                let idx = grid.index(ix, iy, iz);
                seed_if_empty(idx, &grid.occupied, &mut exterior, &mut queue);
            }
        }
    }

    // BFS to flood exterior
    let neighbors: [(i64, i64, i64); 6] = [
        (1, 0, 0),
        (-1, 0, 0),
        (0, 1, 0),
        (0, -1, 0),
        (0, 0, 1),
        (0, 0, -1),
    ];
    while let Some(idx) = queue.pop_front() {
        let iz = idx / (nx * ny);
        let rem = idx % (nx * ny);
        let iy = rem / nx;
        let ix = rem % nx;
        for (dx, dy, dz) in neighbors {
            let nxi = ix as i64 + dx;
            let nyi = iy as i64 + dy;
            let nzi = iz as i64 + dz;
            if nxi < 0
                || nyi < 0
                || nzi < 0
                || nxi >= nx as i64
                || nyi >= ny as i64
                || nzi >= nz as i64
            {
                continue;
            }
            let nidx = grid.index(nxi as usize, nyi as usize, nzi as usize);
            if !grid.occupied[nidx] && !exterior[nidx] {
                exterior.set(nidx, true);
                queue.push_back(nidx);
            }
        }
    }

    // Any voxel that is neither occupied nor exterior is interior → mark it
    for i in 0..n {
        if !grid.occupied[i] && !exterior[i] {
            grid.occupied.set(i, true);
        }
    }
}

/// Separating-Axis Theorem (SAT) test between an AABB and a triangle.
/// Returns `true` if they overlap.
fn aabb_triangle_intersect(
    aabb_min: &[f64; 3],
    aabb_max: &[f64; 3],
    v0: &[f64; 3],
    v1: &[f64; 3],
    v2: &[f64; 3],
) -> bool {
    let center = [
        (aabb_min[0] + aabb_max[0]) * 0.5,
        (aabb_min[1] + aabb_max[1]) * 0.5,
        (aabb_min[2] + aabb_max[2]) * 0.5,
    ];
    let half = [
        (aabb_max[0] - aabb_min[0]) * 0.5,
        (aabb_max[1] - aabb_min[1]) * 0.5,
        (aabb_max[2] - aabb_min[2]) * 0.5,
    ];

    // Translate so AABB is centred at origin
    let tv0 = [v0[0] - center[0], v0[1] - center[1], v0[2] - center[2]];
    let tv1 = [v1[0] - center[0], v1[1] - center[1], v1[2] - center[2]];
    let tv2 = [v2[0] - center[0], v2[1] - center[1], v2[2] - center[2]];

    // Test 3 AABB face normals (coordinate axes)
    for ax in 0..3 {
        let p0 = tv0[ax];
        let p1 = tv1[ax];
        let p2 = tv2[ax];
        let mn = p0.min(p1).min(p2);
        let mx = p0.max(p1).max(p2);
        if mn > half[ax] || mx < -half[ax] {
            return false;
        }
    }

    // Triangle edges
    let e0 = [tv1[0] - tv0[0], tv1[1] - tv0[1], tv1[2] - tv0[2]];
    let e1 = [tv2[0] - tv0[0], tv2[1] - tv0[1], tv2[2] - tv0[2]];
    let e2 = [tv2[0] - tv1[0], tv2[1] - tv1[1], tv2[2] - tv1[2]];

    // Triangle face normal
    let normal = cross3(&e0, &e1);
    if !test_axis_sep(&normal, &half, &tv0, &tv1, &tv2) {
        return false;
    }

    // 9 cross products: edge × AABB axis
    let edges = [e0, e1, e2];
    let aabb_axes = [[1.0_f64, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    for e in &edges {
        for a in &aabb_axes {
            let ax = cross3(e, a);
            let len_sq = ax[0] * ax[0] + ax[1] * ax[1] + ax[2] * ax[2];
            if len_sq < 1e-30 {
                continue; // degenerate axis, skip
            }
            if !test_axis_sep(&ax, &half, &tv0, &tv1, &tv2) {
                return false;
            }
        }
    }

    true
}

#[inline]
fn cross3(a: &[f64; 3], b: &[f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

#[inline]
fn dot3(a: &[f64; 3], b: &[f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// Returns `true` when `axis` does **not** separate the AABB (at origin, half-extents
/// `half`) from the triangle `(tv0, tv1, tv2)`.
#[inline]
fn test_axis_sep(
    axis: &[f64; 3],
    half: &[f64; 3],
    tv0: &[f64; 3],
    tv1: &[f64; 3],
    tv2: &[f64; 3],
) -> bool {
    let r = half[0] * axis[0].abs() + half[1] * axis[1].abs() + half[2] * axis[2].abs();
    let p0 = dot3(axis, tv0);
    let p1 = dot3(axis, tv1);
    let p2 = dot3(axis, tv2);
    let mn = p0.min(p1).min(p2);
    let mx = p0.max(p1).max(p2);
    !(mn > r || mx < -r)
}

/// Convert a `ConvexHull3DVec` into a `ConvexPart`.
///
/// The volume and the centroid come from the hull itself, by the exact
/// signed-tetrahedron decomposition from an interior point, not from the voxel
/// count, so they are correct for any resolution.  `concavity` is the
/// concavity of the voxel cluster the hull belongs to.
fn hull_to_convex_part(hull: &ConvexHull3DVec, concavity: f64) -> ConvexPart {
    let vertices: Vec<[f64; 3]> = hull.vertices.iter().map(|v| [v.x, v.y, v.z]).collect();
    let indices: Vec<[u32; 3]> = hull
        .faces
        .iter()
        .map(|f| {
            let index = |i: usize| u32::try_from(i).unwrap_or(u32::MAX);
            [index(f[0]), index(f[1]), index(f[2])]
        })
        .collect();
    let count = hull.vertices.len().max(1) as f64;
    let interior = hull.vertices.iter().fold(Vec3::zeros(), |sum, v| sum + v) / count;
    let mut six_volume = 0.0;
    let mut weighted = Vec3::zeros();
    for f in &hull.faces {
        let a = hull.vertices[f[0]] - interior;
        let b = hull.vertices[f[1]] - interior;
        let c = hull.vertices[f[2]] - interior;
        let det = a.dot(&b.cross(&c));
        six_volume += det;
        weighted += (a + b + c) * (det / 4.0);
    }
    let volume = six_volume.abs() / 6.0;
    let centroid = if six_volume != 0.0 {
        interior + weighted / six_volume
    } else {
        interior
    };
    ConvexPart {
        vertices,
        concavity,
        indices,
        volume,
        centroid: [centroid.x, centroid.y, centroid.z],
    }
}

// ---------------------------------------------------------------------------
// Unit tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    fn box_tris(min: [f64; 3], max: [f64; 3]) -> Vec<[f64; 9]> {
        let [x0, y0, z0] = min;
        let [x1, y1, z1] = max;
        vec![
            // -X face
            [x0, y0, z0, x0, y1, z0, x0, y1, z1],
            [x0, y0, z0, x0, y1, z1, x0, y0, z1],
            // +X face
            [x1, y0, z0, x1, y1, z1, x1, y1, z0],
            [x1, y0, z0, x1, y0, z1, x1, y1, z1],
            // -Y face
            [x0, y0, z0, x1, y0, z1, x1, y0, z0],
            [x0, y0, z0, x0, y0, z1, x1, y0, z1],
            // +Y face
            [x0, y1, z0, x1, y1, z0, x1, y1, z1],
            [x0, y1, z0, x1, y1, z1, x0, y1, z1],
            // -Z face
            [x0, y0, z0, x1, y0, z0, x1, y1, z0],
            [x0, y0, z0, x1, y1, z0, x0, y1, z0],
            // +Z face
            [x0, y0, z1, x1, y1, z1, x1, y0, z1],
            [x0, y0, z1, x0, y1, z1, x1, y1, z1],
        ]
    }

    /// L-shaped solid of volume 3: the boxes `[0,2] x [0,1] x [0,1]` and
    /// `[0,1] x [1,2] x [0,1]`.  The hull of the L-shaped profile has area 3.5,
    /// so the hull of the solid has volume 3.5.
    fn l_shape_tris() -> Vec<[f64; 9]> {
        let mut tris = box_tris([0.0, 0.0, 0.0], [2.0, 1.0, 1.0]);
        tris.extend(box_tris([0.0, 1.0, 0.0], [1.0, 2.0, 1.0]));
        tris
    }

    /// Deterministic pseudo-random numbers in `[0, 1)` (64-bit LCG).
    fn lcg(seed: u64) -> impl FnMut() -> f64 {
        let mut state = seed;
        move || {
            state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            (state >> 11) as f64 / (1u64 << 53) as f64
        }
    }

    /// Whether `p` lies inside the hull stored in `part` (outward triangles).
    fn part_contains(part: &ConvexPart, p: [f64; 3], tolerance: f64) -> bool {
        part.indices.iter().all(|f| {
            let a = Vec3::from(part.vertices[f[0] as usize]);
            let b = Vec3::from(part.vertices[f[1] as usize]);
            let c = Vec3::from(part.vertices[f[2] as usize]);
            let normal = (b - a).cross(&(c - a)).normalize();
            normal.dot(&(Vec3::from(p) - a)) <= tolerance
        })
    }

    fn assert_vertices_covered(tris: &[[f64; 9]], parts: &ConvexParts) {
        for tri in tris {
            for corner in 0..3 {
                let p = [tri[3 * corner], tri[3 * corner + 1], tri[3 * corner + 2]];
                assert!(
                    parts.parts.iter().any(|part| part_contains(part, p, 1e-9)),
                    "vertex {p:?} is not covered by any part"
                );
            }
        }
    }

    #[test]
    fn test_unit_box_one_part() {
        let tris = box_tris([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let cfg = VHacdConfig {
            resolution: 32,
            ..Default::default()
        };
        let decomposer = VHacdVoxel::new(cfg);
        let result = decomposer
            .decompose(&tris)
            .expect("decompose should succeed");
        assert_eq!(result.parts.len(), 1, "unit box should give 1 part");
        let part = &result.parts[0];
        assert!(
            (part.volume - 1.0).abs() / 1.0 <= 0.10,
            "volume within 10% of 1.0, got {}",
            part.volume
        );
        assert!(
            part.concavity < 1e-9,
            "a box is convex, got {}",
            part.concavity
        );
        assert_eq!(part.vertices.len(), 8);
        assert_eq!(part.indices.len(), 12);
        for axis in 0..3 {
            assert!(
                (part.centroid[axis] - 0.5).abs() < 0.02,
                "centroid {:?}",
                part.centroid
            );
        }
    }

    #[test]
    fn test_empty_input_error() {
        let decomposer = VHacdVoxel::new(VHacdConfig::default());
        assert!(
            decomposer.decompose(&[]).is_err(),
            "empty input should return error"
        );
    }

    #[test]
    fn corner_points_keep_the_hull_of_all_voxel_corners() {
        let mut next = lcg(0x0cd0_0b10);
        let grid = VoxelGrid::new([-1.5, 0.25, 3.0], 0.5, [11, 8, 9]);
        let mut voxels: Vec<VoxelCoord> = Vec::new();
        for ix in 0..11usize {
            for iy in 0..8usize {
                for iz in 0..9usize {
                    let (dx, dy, dz) = (ix as f64 - 5.0, iy as f64 - 3.5, iz as f64 - 4.0);
                    // A ball with a bite taken out, plus random speckles.
                    let ball = dx * dx + dy * dy + dz * dz <= 12.0;
                    let bite = dx > 2.0 && dy > 0.0;
                    if (ball && !bite) || next() < 0.05 {
                        voxels.push([ix, iy, iz]);
                    }
                }
            }
        }
        let mut all_corners: Vec<Vec3> = Vec::new();
        for v in &voxels {
            for dx in 0..2 {
                for dy in 0..2 {
                    for dz in 0..2 {
                        all_corners.push(Vec3::new(
                            grid.origin[0] + (v[0] + dx) as f64 * grid.step,
                            grid.origin[1] + (v[1] + dy) as f64 * grid.step,
                            grid.origin[2] + (v[2] + dz) as f64 * grid.step,
                        ));
                    }
                }
            }
        }
        let reduced = corner_points(&voxels, &grid);
        assert!(
            reduced.len() < all_corners.len() / 4,
            "{} points",
            reduced.len()
        );
        let full_hull = ConvexHull3DVec::build(&all_corners).expect("full hull");
        let reduced_hull = ConvexHull3DVec::build(&reduced).expect("reduced hull");
        let (full, small) = (full_hull.volume(), reduced_hull.volume());
        assert!((full - small).abs() <= 1e-9 * full, "{full} vs {small}");
        assert!(
            all_corners.iter().all(|&p| reduced_hull.contains_point(p)),
            "every corner lies inside the hull of the reduced points"
        );
    }

    #[test]
    fn split_positions_cover_even_spacing_and_the_largest_steps() {
        // 10 slices of 6 voxels followed by 10 slices of 3: one big step at 10.
        let mut counts = vec![6usize; 10];
        counts.extend(vec![3usize; 10]);
        let positions = split_positions(&counts, 0, 19);
        assert!(
            positions.windows(2).all(|w| w[0] < w[1]),
            "sorted and unique"
        );
        assert!(positions.iter().all(|&p| (1..=19).contains(&p)));
        assert!(positions.contains(&10), "the step position is a candidate");
        assert!(positions.len() >= EVEN_SPLIT_CANDIDATES);
        assert!(split_positions(&counts, 5, 5).is_empty());
    }

    #[test]
    fn concavity_of_an_l_shape_is_measured_against_its_real_hull() {
        // No splitting allowed: the single part is the hull of the whole L.
        let tris = l_shape_tris();
        let cfg = VHacdConfig {
            max_depth: 0,
            ..Default::default()
        };
        let result = VHacdVoxel::new(cfg)
            .decompose(&tris)
            .expect("decompose should succeed");
        assert_eq!(result.parts.len(), 1);
        let part = &result.parts[0];
        // The hull of the solid is 3.5; the bounding box would be 4.0.
        assert!(
            (part.volume - 3.5).abs() <= 0.35,
            "hull volume near 3.5, got {}",
            part.volume
        );
        // The empty space of the hull is 1 - 3/3.5 = 0.143 of it; the space
        // within a voxel of the L is not counted, which leaves most of it.
        assert!(
            part.concavity > 0.09 && part.concavity < 0.143,
            "concavity of the L between 0.09 and 0.143, got {}",
            part.concavity
        );
    }

    #[test]
    fn l_shape_is_split_into_convex_parts_with_the_right_total_volume() {
        let tris = l_shape_tris();
        let cfg = VHacdConfig::default();
        let threshold = cfg.concavity_threshold;
        let result = VHacdVoxel::new(cfg)
            .decompose(&tris)
            .expect("decompose should succeed");
        assert!(
            result.parts.len() >= 2,
            "an L is not convex, got {} part(s)",
            result.parts.len()
        );
        let total: f64 = result.parts.iter().map(|p| p.volume).sum();
        assert!(
            (total - 3.0).abs() / 3.0 <= 0.10,
            "total volume within 10% of 3.0, got {total}"
        );
        for part in &result.parts {
            assert!(
                part.concavity <= threshold + 1e-12,
                "every leaf is below the threshold, got {}",
                part.concavity
            );
            assert!(part.volume > 0.0);
        }
        assert_vertices_covered(&tris, &result);
    }

    #[test]
    fn part_limit_keeps_the_hulls_covering_the_input() {
        let tris = l_shape_tris();
        let cfg = VHacdConfig {
            max_parts: 2,
            ..Default::default()
        };
        let result = VHacdVoxel::new(cfg)
            .decompose(&tris)
            .expect("decompose should succeed");
        assert!(result.parts.len() <= 2);
        assert_vertices_covered(&tris, &result);
    }

    #[test]
    fn smallest_components_are_dropped_first_when_over_the_part_limit() {
        let mut tris = box_tris([0.0, 0.0, 0.0], [2.0, 2.0, 2.0]);
        tris.extend(box_tris([3.0, 0.0, 0.0], [4.0, 1.0, 1.0]));
        tris.extend(box_tris([5.0, 0.0, 0.0], [5.5, 0.5, 0.5]));
        let cfg = VHacdConfig {
            max_parts: 2,
            ..Default::default()
        };
        let result = VHacdVoxel::new(cfg)
            .decompose(&tris)
            .expect("decompose should succeed");
        assert_eq!(result.parts.len(), 2);
        assert!(
            result.parts.iter().all(|p| p.centroid[0] < 4.5),
            "the small box is gone"
        );
        let mut volumes: Vec<f64> = result.parts.iter().map(|p| p.volume).collect();
        volumes.sort_by(f64::total_cmp);
        assert!(volumes[0] > 0.9 && volumes[1] > 7.5, "volumes {volumes:?}");
    }

    #[test]
    fn flat_input_gives_one_thin_part() {
        let tris = vec![[0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0]];
        let result = VHacdVoxel::new(VHacdConfig::default())
            .decompose(&tris)
            .expect("a flat triangle is valid input");
        assert_eq!(result.parts.len(), 1);
        assert!(result.parts[0].volume > 0.0);
    }

    #[test]
    fn zero_resolution_is_treated_as_one() {
        let tris = box_tris([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let coarse = |resolution: usize| {
            let cfg = VHacdConfig {
                resolution,
                ..Default::default()
            };
            VHacdVoxel::new(cfg)
                .decompose(&tris)
                .expect("decompose should succeed")
        };
        let zero = coarse(0);
        let one = coarse(1);
        assert_eq!(zero.parts.len(), 1, "a part, not an empty result");
        assert_eq!(one.parts.len(), 1);
        assert_eq!(zero.parts[0].volume, one.parts[0].volume);
        assert!(zero.parts[0].volume > 0.0);
    }

    #[test]
    fn hull_part_reports_the_centroid_of_the_solid() {
        // Square pyramid with base [0,1]^2 at z = 0 and apex (0.5, 0.5, 1):
        // volume 1/3, centroid at a quarter of the height (the average of the
        // vertices would be at a fifth).
        let points = [
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(1.0, 1.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
            Vec3::new(0.5, 0.5, 1.0),
        ];
        let hull = ConvexHull3DVec::build(&points).expect("pyramid hull");
        let part = hull_to_convex_part(&hull, 0.25);
        assert!((part.volume - 1.0 / 3.0).abs() < 1e-12);
        assert!((part.centroid[0] - 0.5).abs() < 1e-12);
        assert!((part.centroid[1] - 0.5).abs() < 1e-12);
        assert!((part.centroid[2] - 0.25).abs() < 1e-12);
        assert_eq!(part.concavity, 0.25);
        assert_eq!(part.indices.len(), hull.faces.len());
    }

    // -----------------------------------------------------------------------
    // Convex solids that are not aligned with the grid, and rings
    // -----------------------------------------------------------------------

    /// Rotate `p` by `angles` about the x, y and z axes (in that order).
    fn rotate(p: [f64; 3], angles: [f64; 3]) -> [f64; 3] {
        let (sx, cx) = angles[0].sin_cos();
        let (sy, cy) = angles[1].sin_cos();
        let (sz, cz) = angles[2].sin_cos();
        let (y1, z1) = (cx * p[1] - sx * p[2], sx * p[1] + cx * p[2]);
        let (x2, z2) = (cy * p[0] + sy * z1, -sy * p[0] + cy * z1);
        [cz * x2 - sz * y1, sz * x2 + cz * y1, z2]
    }

    fn rotated(tris: &[[f64; 9]], angles: [f64; 3]) -> Vec<[f64; 9]> {
        tris.iter()
            .map(|t| {
                let a = rotate([t[0], t[1], t[2]], angles);
                let b = rotate([t[3], t[4], t[5]], angles);
                let c = rotate([t[6], t[7], t[8]], angles);
                [a[0], a[1], a[2], b[0], b[1], b[2], c[0], c[1], c[2]]
            })
            .collect()
    }

    /// Regular octahedron `|x| + |y| + |z| <= 1` of volume 4/3.
    fn octahedron_tris() -> Vec<[f64; 9]> {
        let v = [
            [1.0, 0.0, 0.0],
            [-1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, -1.0, 0.0],
            [0.0, 0.0, 1.0],
            [0.0, 0.0, -1.0],
        ];
        let faces = [
            [0, 2, 4],
            [2, 1, 4],
            [1, 3, 4],
            [3, 0, 4],
            [2, 0, 5],
            [1, 2, 5],
            [3, 1, 5],
            [0, 3, 5],
        ];
        faces
            .iter()
            .map(|f| {
                let (a, b, c) = (v[f[0]], v[f[1]], v[f[2]]);
                [a[0], a[1], a[2], b[0], b[1], b[2], c[0], c[1], c[2]]
            })
            .collect()
    }

    /// Torus of revolution about the z axis with `nu` x `nv` quads.
    fn torus_tris(major: f64, minor: f64, nu: usize, nv: usize) -> Vec<[f64; 9]> {
        let point = |i: usize, j: usize| {
            let u = std::f64::consts::TAU * (i % nu) as f64 / nu as f64;
            let v = std::f64::consts::TAU * (j % nv) as f64 / nv as f64;
            let ring = major + minor * v.cos();
            [ring * u.cos(), ring * u.sin(), minor * v.sin()]
        };
        let mut tris = Vec::new();
        for i in 0..nu {
            for j in 0..nv {
                let (a, b) = (point(i, j), point(i + 1, j));
                let (c, d) = (point(i + 1, j + 1), point(i, j + 1));
                tris.push([a[0], a[1], a[2], b[0], b[1], b[2], c[0], c[1], c[2]]);
                tris.push([a[0], a[1], a[2], c[0], c[1], c[2], d[0], d[1], d[2]]);
            }
        }
        tris
    }

    fn assert_single_convex_part(tris: &[[f64; 9]], resolution: usize, volume: f64) {
        let cfg = VHacdConfig {
            resolution,
            ..Default::default()
        };
        let result = VHacdVoxel::new(cfg)
            .decompose(tris)
            .expect("decompose should succeed");
        assert_eq!(result.parts.len(), 1, "a convex solid is one part");
        let part = &result.parts[0];
        assert!(part.concavity < 1e-9, "concavity {}", part.concavity);
        assert!(
            part.vertices.len() > 8,
            "the hull follows the solid, not its bounding box: {} vertices",
            part.vertices.len()
        );
        // The hull of the voxels covers the solid and exceeds it by at most a
        // voxel on each side.
        assert!(
            part.volume >= volume && part.volume <= 1.8 * volume,
            "hull volume {} for a solid of {volume}",
            part.volume
        );
        assert_vertices_covered(tris, &result);
    }

    #[test]
    fn tilted_convex_solids_are_one_part() {
        // The voxels of a solid whose faces are not aligned with the grid form
        // staircases, so the hull of the voxels is far from filled: the plain
        // fraction of empty hull volume is about 17 % for the rotated cube at
        // this resolution.  That layer is not concavity and no split removes
        // it, so the solid must stay in one piece.
        let cube = box_tris([-0.5, -0.5, -0.5], [0.5, 0.5, 0.5]);
        assert_single_convex_part(&rotated(&cube, [0.5, 0.3, 0.7]), 32, 1.0);
        assert_single_convex_part(&octahedron_tris(), 32, 4.0 / 3.0);
        assert_single_convex_part(&rotated(&octahedron_tris(), [0.2, 0.9, 0.4]), 32, 4.0 / 3.0);
    }

    #[test]
    fn rotated_l_shape_is_still_split_and_covered() {
        let tris = rotated(&l_shape_tris(), [0.5, 0.3, 0.7]);
        let cfg = VHacdConfig {
            resolution: 32,
            ..Default::default()
        };
        let result = VHacdVoxel::new(cfg)
            .decompose(&tris)
            .expect("decompose should succeed");
        assert!(
            result.parts.len() >= 2,
            "an L is not convex in any orientation, got {} part(s)",
            result.parts.len()
        );
        assert!(
            result.parts.len() <= 8,
            "{} parts for two boxes",
            result.parts.len()
        );
        assert!(result.parts.iter().all(|p| p.volume > 0.0));
        assert_vertices_covered(&tris, &result);
    }

    #[test]
    fn ring_is_cut_around_its_hole() {
        // The hull of a ring covers its hole; the parts must not.
        let tris = torus_tris(1.0, 0.35, 24, 12);
        let cfg = VHacdConfig {
            resolution: 24,
            ..Default::default()
        };
        let result = VHacdVoxel::new(cfg)
            .decompose(&tris)
            .expect("decompose should succeed");
        assert!(
            result.parts.len() >= 4,
            "a ring needs several parts, got {}",
            result.parts.len()
        );
        assert!(
            result
                .parts
                .iter()
                .all(|part| !part_contains(part, [0.0, 0.0, 0.0], 0.0)),
            "the centre of the hole must stay free"
        );
        assert_vertices_covered(&tris, &result);
    }

    /// Voxels `(ix, iy, iz)` of `grid` for which `occupied` holds.
    fn voxels_where(
        grid: &VoxelGrid,
        mut occupied: impl FnMut(usize, usize, usize) -> bool,
    ) -> Vec<VoxelCoord> {
        let mut voxels = Vec::new();
        for ix in 0..grid.dims[0] {
            for iy in 0..grid.dims[1] {
                for iz in 0..grid.dims[2] {
                    if occupied(ix, iy, iz) {
                        voxels.push([ix, iy, iz]);
                    }
                }
            }
        }
        voxels
    }

    /// Volume of the cells inside the hull that are not within one voxel of an
    /// occupied cell, counted one cell at a time with the signed distance of the
    /// hull.  A cell whose centre is within `1e-9` of the hull surface counts
    /// when `inclusive` is set and does not otherwise.
    fn brute_force_empty_hull_volume(
        hull: &ConvexHull3DVec,
        voxels: &[VoxelCoord],
        grid: &VoxelGrid,
        inclusive: bool,
    ) -> f64 {
        let occupied: std::collections::HashSet<VoxelCoord> = voxels.iter().copied().collect();
        let near_occupied = |ix: usize, iy: usize, iz: usize| {
            (ix.saturating_sub(1)..=ix + 1).any(|x| {
                (iy.saturating_sub(1)..=iy + 1)
                    .any(|y| (iz.saturating_sub(1)..=iz + 1).any(|z| occupied.contains(&[x, y, z])))
            })
        };
        let margin = if inclusive { 1e-9 } else { -1e-9 };
        let mut cells = 0usize;
        for ix in 0..grid.dims[0] {
            for iy in 0..grid.dims[1] {
                for iz in 0..grid.dims[2] {
                    if near_occupied(ix, iy, iz) {
                        continue;
                    }
                    let [x, y, z] = grid.voxel_center(ix, iy, iz);
                    if hull.signed_distance(Vec3::new(x, y, z)) <= margin {
                        cells += 1;
                    }
                }
            }
        }
        cells as f64 * grid.voxel_volume()
    }

    #[test]
    fn empty_hull_volume_agrees_with_a_brute_force_count() {
        let grid = VoxelGrid::new([-0.5, 0.25, 2.0], 0.5, [14, 14, 9]);
        // L-shaped slab, a ring, a tilted convex wedge and speckles.
        let l_shape = voxels_where(&grid, |ix, iy, iz| {
            iz < 5 && ix < 12 && iy < 12 && (ix < 4 || iy < 4)
        });
        let ring = voxels_where(&grid, |ix, iy, iz| {
            let (dx, dy) = (ix as f64 - 6.5, iy as f64 - 6.5);
            let r2 = dx * dx + dy * dy;
            iz < 6 && (9.0..=30.0).contains(&r2)
        });
        let wedge = voxels_where(&grid, |ix, iy, iz| iz < 6 && 2 * ix + 3 * iy <= 30);
        let mut next = lcg(0x00de_e9a9);
        let speckles = voxels_where(&grid, |_, _, _| next() < 0.03);
        for (name, voxels) in [
            ("L", l_shape),
            ("ring", ring),
            ("wedge", wedge),
            ("speckles", speckles),
        ] {
            let points = corner_points(&voxels, &grid);
            let hull = ConvexHull3DVec::build(&points).expect("hull");
            let measured = empty_hull_volume(&hull, &voxels, &grid);
            let lower = brute_force_empty_hull_volume(&hull, &voxels, &grid, false);
            let upper = brute_force_empty_hull_volume(&hull, &voxels, &grid, true);
            assert!(
                lower - 1e-9 <= measured && measured <= upper + 1e-9,
                "{name}: empty hull volume {measured} not in [{lower}, {upper}]"
            );
            if name == "wedge" {
                assert_eq!(measured, 0.0, "a convex staircase leaves no deep gap");
            } else {
                assert!(
                    measured > 0.0,
                    "{name} has empty space away from its voxels"
                );
            }
        }
    }

    #[test]
    fn thin_washer_is_cut_around_its_hole() {
        // A washer three voxels thick: its hull is a disc, and the hole is
        // empty space far from any voxel although it is only a voxel and a
        // half from the top and the bottom of the hull.
        let (outer, inner, half_thickness) = (1.0_f64, 0.5_f64, 0.05_f64);
        let segments = 32;
        let point =
            |radius: f64, angle: f64, z: f64| [radius * angle.cos(), radius * angle.sin(), z];
        let mut tris: Vec<[f64; 9]> = Vec::new();
        for k in 0..segments {
            let a0 = std::f64::consts::TAU * k as f64 / segments as f64;
            let a1 = std::f64::consts::TAU * (k + 1) as f64 / segments as f64;
            let mut quad = |p: [f64; 3], q: [f64; 3], r: [f64; 3], t: [f64; 3]| {
                tris.push([p[0], p[1], p[2], q[0], q[1], q[2], r[0], r[1], r[2]]);
                tris.push([p[0], p[1], p[2], r[0], r[1], r[2], t[0], t[1], t[2]]);
            };
            let h = half_thickness;
            quad(
                point(inner, a0, h),
                point(outer, a0, h),
                point(outer, a1, h),
                point(inner, a1, h),
            );
            quad(
                point(inner, a0, -h),
                point(inner, a1, -h),
                point(outer, a1, -h),
                point(outer, a0, -h),
            );
            quad(
                point(outer, a0, -h),
                point(outer, a1, -h),
                point(outer, a1, h),
                point(outer, a0, h),
            );
            quad(
                point(inner, a0, -h),
                point(inner, a0, h),
                point(inner, a1, h),
                point(inner, a1, -h),
            );
        }
        let cfg = VHacdConfig {
            resolution: 64,
            ..Default::default()
        };
        let result = VHacdVoxel::new(cfg)
            .decompose(&tris)
            .expect("decompose should succeed");
        assert!(
            result.parts.len() >= 2,
            "a washer is not convex, got {} part(s)",
            result.parts.len()
        );
        assert!(
            result
                .parts
                .iter()
                .all(|part| !part_contains(part, [0.0, 0.0, 0.0], 0.0)),
            "the centre of the hole must stay free"
        );
        assert_vertices_covered(&tris, &result);
    }
}
