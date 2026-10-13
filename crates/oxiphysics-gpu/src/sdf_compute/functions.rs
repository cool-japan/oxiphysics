//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use rayon::prelude::*;

/// `fast_sweeping_update` now lives in the [`eikonal`](super::eikonal) module; it is
/// re-exported here so that the path it was first published under keeps working.
pub use super::eikonal::fast_sweeping_update;
use super::eikonal::sdf_redistance;
use super::mc_tables::{MC_CORNER_OFFSETS, MC_EDGE_CORNERS, MC_EDGE_TABLE, MC_TRI_TABLE};
use super::types::{SdfGrid, SphereTraceResult, Triangle};

/// Compute the union of two SDFs: `min(a, b)` pointwise.
pub fn union_sdf(a: &SdfGrid, b: &SdfGrid) -> SdfGrid {
    assert_eq!(a.nx, b.nx, "union_sdf: nx mismatch");
    assert_eq!(a.ny, b.ny, "union_sdf: ny mismatch");
    assert_eq!(a.nz, b.nz, "union_sdf: nz mismatch");
    let values: Vec<f64> = a
        .values
        .par_iter()
        .zip(b.values.par_iter())
        .map(|(&av, &bv)| av.min(bv))
        .collect();
    SdfGrid {
        nx: a.nx,
        ny: a.ny,
        nz: a.nz,
        dx: a.dx,
        origin: a.origin,
        values,
    }
}
/// Compute the intersection of two SDFs: `max(a, b)` pointwise.
pub fn intersection_sdf(a: &SdfGrid, b: &SdfGrid) -> SdfGrid {
    assert_eq!(a.nx, b.nx, "intersection_sdf: nx mismatch");
    assert_eq!(a.ny, b.ny, "intersection_sdf: ny mismatch");
    assert_eq!(a.nz, b.nz, "intersection_sdf: nz mismatch");
    let values: Vec<f64> = a
        .values
        .par_iter()
        .zip(b.values.par_iter())
        .map(|(&av, &bv)| av.max(bv))
        .collect();
    SdfGrid {
        nx: a.nx,
        ny: a.ny,
        nz: a.nz,
        dx: a.dx,
        origin: a.origin,
        values,
    }
}
/// Compute the subtraction of SDF b from SDF a: `max(a, -b)` pointwise.
///
/// The result is the region inside `a` but outside `b`.
pub fn subtraction_sdf(a: &SdfGrid, b: &SdfGrid) -> SdfGrid {
    assert_eq!(a.nx, b.nx, "subtraction_sdf: nx mismatch");
    assert_eq!(a.ny, b.ny, "subtraction_sdf: ny mismatch");
    assert_eq!(a.nz, b.nz, "subtraction_sdf: nz mismatch");
    let values: Vec<f64> = a
        .values
        .par_iter()
        .zip(b.values.par_iter())
        .map(|(&av, &bv)| av.max(-bv))
        .collect();
    SdfGrid {
        nx: a.nx,
        ny: a.ny,
        nz: a.nz,
        dx: a.dx,
        origin: a.origin,
        values,
    }
}
/// Compute the shell SDF: `|sdf(p)| - thickness/2`.
pub fn shell_sdf(grid: &SdfGrid, thickness: f64) -> SdfGrid {
    let half = thickness / 2.0;
    let values: Vec<f64> = grid.values.par_iter().map(|&v| v.abs() - half).collect();
    SdfGrid {
        nx: grid.nx,
        ny: grid.ny,
        nz: grid.nz,
        dx: grid.dx,
        origin: grid.origin,
        values,
    }
}
/// Smooth union of two SDFs using polynomial smoothing.
///
/// `k` controls the smoothing radius. Larger `k` means sharper transition.
pub fn smooth_union_sdf(a: &SdfGrid, b: &SdfGrid, k: f64) -> SdfGrid {
    assert_eq!(a.nx, b.nx);
    assert_eq!(a.ny, b.ny);
    assert_eq!(a.nz, b.nz);
    let values: Vec<f64> = a
        .values
        .par_iter()
        .zip(b.values.par_iter())
        .map(|(&av, &bv)| {
            let h = (0.5 + 0.5 * (bv - av) / k).clamp(0.0, 1.0);
            bv * (1.0 - h) + av * h - k * h * (1.0 - h)
        })
        .collect();
    SdfGrid {
        nx: a.nx,
        ny: a.ny,
        nz: a.nz,
        dx: a.dx,
        origin: a.origin,
        values,
    }
}
/// Perform sphere tracing (ray marching) against an SDF grid.
///
/// Traces a ray from `origin` in `direction` (must be unit-length)
/// until the SDF value is below `surface_threshold` or `max_t` is reached.
pub fn sphere_trace(
    grid: &SdfGrid,
    ray_origin: [f64; 3],
    ray_direction: [f64; 3],
    max_t: f64,
    max_iterations: usize,
    surface_threshold: f64,
) -> SphereTraceResult {
    let mut t = 0.0;
    let mut pos = ray_origin;
    for iter in 0..max_iterations {
        let sdf_val = match grid.sample(pos) {
            Some(v) => v,
            None => {
                return SphereTraceResult {
                    hit: false,
                    position: pos,
                    t,
                    iterations: iter,
                };
            }
        };
        if sdf_val < surface_threshold {
            return SphereTraceResult {
                hit: true,
                position: pos,
                t,
                iterations: iter,
            };
        }
        t += sdf_val;
        if t > max_t {
            return SphereTraceResult {
                hit: false,
                position: pos,
                t,
                iterations: iter,
            };
        }
        pos = [
            ray_origin[0] + ray_direction[0] * t,
            ray_origin[1] + ray_direction[1] * t,
            ray_origin[2] + ray_direction[2] * t,
        ];
    }
    SphereTraceResult {
        hit: false,
        position: pos,
        t,
        iterations: max_iterations,
    }
}
/// Convert a triangle mesh to an SDF by computing the distance from
/// each grid cell to the nearest triangle.
///
/// This is a brute-force O(N*M) approach where N is grid cells and
/// M is triangles. For large meshes, use spatial acceleration.
///
/// `vertices` - vertex positions
/// `triangles` - triangle indices (3 per triangle)
pub fn mesh_to_sdf(grid: &mut SdfGrid, vertices: &[[f64; 3]], triangles: &[[usize; 3]]) {
    let ny = grid.ny;
    let nz = grid.nz;
    let dx = grid.dx;
    let origin = grid.origin;
    grid.values.par_iter_mut().enumerate().for_each(|(idx, v)| {
        let i = idx / (ny * nz);
        let j = (idx / nz) % ny;
        let k = idx % nz;
        let p = [
            origin[0] + (i as f64 + 0.5) * dx,
            origin[1] + (j as f64 + 0.5) * dx,
            origin[2] + (k as f64 + 0.5) * dx,
        ];
        let mut min_dist = f64::MAX;
        for tri in triangles {
            let a = vertices[tri[0]];
            let b = vertices[tri[1]];
            let c = vertices[tri[2]];
            let dist = point_triangle_distance(&p, &a, &b, &c);
            if dist < min_dist {
                min_dist = dist;
            }
        }
        *v = min_dist;
    });
}
/// Compute the distance from a point to a triangle.
pub(super) fn point_triangle_distance(
    p: &[f64; 3],
    a: &[f64; 3],
    b: &[f64; 3],
    c: &[f64; 3],
) -> f64 {
    let ab = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
    let ac = [c[0] - a[0], c[1] - a[1], c[2] - a[2]];
    let ap = [p[0] - a[0], p[1] - a[1], p[2] - a[2]];
    let d1 = dot3(&ab, &ap);
    let d2 = dot3(&ac, &ap);
    if d1 <= 0.0 && d2 <= 0.0 {
        return dist3(p, a);
    }
    let bp = [p[0] - b[0], p[1] - b[1], p[2] - b[2]];
    let d3 = dot3(&ab, &bp);
    let d4 = dot3(&ac, &bp);
    if d3 >= 0.0 && d4 <= d3 {
        return dist3(p, b);
    }
    let vc = d1 * d4 - d3 * d2;
    if vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0 {
        let v = d1 / (d1 - d3);
        let proj = [a[0] + ab[0] * v, a[1] + ab[1] * v, a[2] + ab[2] * v];
        return dist3(p, &proj);
    }
    let cp = [p[0] - c[0], p[1] - c[1], p[2] - c[2]];
    let d5 = dot3(&ab, &cp);
    let d6 = dot3(&ac, &cp);
    if d6 >= 0.0 && d5 <= d6 {
        return dist3(p, c);
    }
    let vb = d5 * d2 - d1 * d6;
    if vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0 {
        let w = d2 / (d2 - d6);
        let proj = [a[0] + ac[0] * w, a[1] + ac[1] * w, a[2] + ac[2] * w];
        return dist3(p, &proj);
    }
    let va = d3 * d6 - d5 * d4;
    if va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0 {
        let w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        let bc = [c[0] - b[0], c[1] - b[1], c[2] - b[2]];
        let proj = [b[0] + bc[0] * w, b[1] + bc[1] * w, b[2] + bc[2] * w];
        return dist3(p, &proj);
    }
    let denom = 1.0 / (va + vb + vc);
    let v = vb * denom;
    let w = vc * denom;
    let proj = [
        a[0] + ab[0] * v + ac[0] * w,
        a[1] + ab[1] * v + ac[1] * w,
        a[2] + ab[2] * v + ac[2] * w,
    ];
    dist3(p, &proj)
}
pub(super) fn dot3(a: &[f64; 3], b: &[f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
pub(super) fn dist3(a: &[f64; 3], b: &[f64; 3]) -> f64 {
    let dx = a[0] - b[0];
    let dy = a[1] - b[1];
    let dz = a[2] - b[2];
    (dx * dx + dy * dy + dz * dz).sqrt()
}
/// Evaluate an SDF at a batch of query points.
///
/// Returns the SDF value at each query point, or `f64::MAX` if outside the grid.
pub fn evaluate_sdf_batch(grid: &SdfGrid, points: &[[f64; 3]]) -> Vec<f64> {
    points
        .par_iter()
        .map(|&p| grid.sample(p).unwrap_or(f64::MAX))
        .collect()
}
/// Compute the SDF gradient at a batch of query points.
pub fn gradient_sdf_batch(grid: &SdfGrid, points: &[[f64; 3]]) -> Vec<[f64; 3]> {
    points
        .par_iter()
        .map(|&p| grid.gradient_at_point(p).unwrap_or([0.0; 3]))
        .collect()
}
/// Count the number of cells where the SDF is negative (inside the surface).
pub fn count_interior_cells(grid: &SdfGrid) -> usize {
    grid.values.par_iter().filter(|&&v| v < 0.0).count()
}
/// Compute the approximate volume enclosed by the zero level-set.
///
/// Simply counts interior cells and multiplies by cell volume.
pub fn approximate_volume(grid: &SdfGrid) -> f64 {
    let count = count_interior_cells(grid);
    count as f64 * grid.dx * grid.dx * grid.dx
}
/// Interpolate the position of an edge intersection given two corner SDF values.
///
/// `p1`/`p2` are the end points of the edge and `val1`/`val2` the field values
/// there.  The result lies on the edge (the interpolation parameter is clamped to
/// `[0, 1]`); if the parameter is not finite (NaN field values) the midpoint is
/// returned.
#[inline]
pub(super) fn interpolate_vertex(
    p1: [f64; 3],
    p2: [f64; 3],
    val1: f64,
    val2: f64,
    iso: f64,
) -> [f64; 3] {
    let t = (iso - val1) / (val2 - val1);
    let t = if t.is_finite() {
        t.clamp(0.0, 1.0)
    } else {
        0.5
    };
    [
        p1[0] + t * (p2[0] - p1[0]),
        p1[1] + t * (p2[1] - p1[1]),
        p1[2] + t * (p2[2] - p1[2]),
    ]
}
/// Triangles generated by one marching-cubes cell (at most five).
pub(super) struct CellTriangles {
    /// Triangle vertices in world space; only the first `count` are valid.
    pub(super) triangles: [[[f64; 3]; 3]; 5],
    /// Number of valid triangles.
    pub(super) count: usize,
}
/// Run marching cubes on the single cell whose lowest corner is grid point
/// `(i, j, k)`.
///
/// Corner `c` of the cell is the grid point `(i, j, k) + MC_CORNER_OFFSETS[c]`
/// and is inside when its value is below `isovalue`.  Every crossed edge yields
/// one vertex, linearly interpolated *from its end point with the lower flat grid
/// index*: the two to four cells that share an edge therefore compute bit-for-bit
/// identical vertices, which is what makes the extracted mesh watertight.  The
/// caller guarantees that grid point `(i + 1, j + 1, k + 1)` exists.
pub(super) fn polygonize_cell(
    grid: &SdfGrid,
    i: usize,
    j: usize,
    k: usize,
    isovalue: f64,
) -> CellTriangles {
    let mut cell = CellTriangles {
        triangles: [[[0.0; 3]; 3]; 5],
        count: 0,
    };
    let mut values = [0.0_f64; 8];
    let mut case = 0_usize;
    for (corner, offset) in MC_CORNER_OFFSETS.iter().enumerate() {
        let value = grid.get(i + offset[0], j + offset[1], k + offset[2]);
        values[corner] = value;
        if value < isovalue {
            case |= 1 << corner;
        }
    }
    let crossed = MC_EDGE_TABLE[case];
    if crossed == 0 {
        return cell;
    }
    let corner_position = |corner: usize| -> [f64; 3] {
        let offset = MC_CORNER_OFFSETS[corner];
        grid.world_pos(i + offset[0], j + offset[1], k + offset[2])
    };
    let mut vertices = [[0.0_f64; 3]; 12];
    for (edge, corners) in MC_EDGE_CORNERS.iter().enumerate() {
        if crossed & (1 << edge) == 0 {
            continue;
        }
        let [a, b] = *corners;
        // Lexicographic order of the corner offsets equals the flat-index order.
        let (lower, upper) = if MC_CORNER_OFFSETS[a] <= MC_CORNER_OFFSETS[b] {
            (a, b)
        } else {
            (b, a)
        };
        vertices[edge] = interpolate_vertex(
            corner_position(lower),
            corner_position(upper),
            values[lower],
            values[upper],
            isovalue,
        );
    }
    for edges in MC_TRI_TABLE[case].chunks_exact(3) {
        let (Ok(e0), Ok(e1), Ok(e2)) = (
            usize::try_from(edges[0]),
            usize::try_from(edges[1]),
            usize::try_from(edges[2]),
        ) else {
            break;
        };
        if cell.count == cell.triangles.len() {
            break;
        }
        // The table winds into the solid; reverse it so the triangle faces outwards.
        cell.triangles[cell.count] = [vertices[e0], vertices[e2], vertices[e1]];
        cell.count += 1;
    }
    cell
}
/// Extract a triangle mesh from the SDF grid at the given isovalue
/// using the marching cubes algorithm (Lorensen & Cline, with the standard
/// 256-case triangle table).
///
/// Returns a list of triangles. For a signed distance field, `isovalue = 0.0`
/// extracts the zero level-set (the surface).  A grid point is inside when its
/// value is below `isovalue`, and the triangles are wound so that their
/// right-handed normals point towards increasing field values, i.e. out of the
/// solid.  Vertices are the linearly interpolated crossings of the grid edges and
/// are shared bit-for-bit by neighbouring cells, so a closed surface that the
/// grid resolves is watertight: every undirected mesh edge belongs to exactly
/// two triangles.  (Under-resolved features can hit the topological ambiguity of
/// the original algorithm, which may leave small holes.)
///
/// Cells are processed in parallel; the output order is deterministic (cells in
/// ascending `i`, then `j`, then `k`).
pub fn marching_cubes(grid: &SdfGrid, isovalue: f64) -> Vec<Triangle> {
    if grid.nx < 2 || grid.ny < 2 || grid.nz < 2 {
        return Vec::new();
    }
    let ny = grid.ny;
    let nz = grid.nz;
    let slabs: Vec<Vec<Triangle>> = (0..grid.nx - 1)
        .into_par_iter()
        .map(|i| {
            let mut slab = Vec::new();
            for j in 0..ny - 1 {
                for k in 0..nz - 1 {
                    let cell = polygonize_cell(grid, i, j, k, isovalue);
                    slab.extend(cell.triangles[..cell.count].iter().map(|&v| Triangle { v }));
                }
            }
            slab
        })
        .collect();
    slabs.into_iter().flatten().collect()
}
/// Count the number of triangles in the extracted mesh.
pub fn mesh_triangle_count(grid: &SdfGrid, isovalue: f64) -> usize {
    marching_cubes(grid, isovalue).len()
}
/// Apply a 3×3×3 box-filter (Gaussian-like) convolution to an SDF grid.
///
/// Each cell is replaced by the weighted average of its 3×3×3 neighbourhood.
/// `sigma` controls the Gaussian kernel width (in cells).
pub fn sdf_gaussian_blur(grid: &SdfGrid, sigma: f64) -> SdfGrid {
    let nx = grid.nx;
    let ny = grid.ny;
    let nz = grid.nz;
    let mut kernel = [[[0.0f64; 3]; 3]; 3];
    let mut kernel_sum = 0.0;
    let s2 = 2.0 * sigma * sigma;
    for di in -1i32..=1 {
        for dj in -1i32..=1 {
            for dk in -1i32..=1 {
                let r2 = (di * di + dj * dj + dk * dk) as f64;
                let w = (-r2 / s2).exp();
                kernel[(di + 1) as usize][(dj + 1) as usize][(dk + 1) as usize] = w;
                kernel_sum += w;
            }
        }
    }
    let mut out = SdfGrid::new(nx, ny, nz, grid.dx, grid.origin);
    for i in 0..nx {
        for j in 0..ny {
            for k in 0..nz {
                let mut acc = 0.0;
                let mut wt = 0.0;
                for di in -1i32..=1 {
                    for dj in -1i32..=1 {
                        for dk in -1i32..=1 {
                            let ni = i as i32 + di;
                            let nj = j as i32 + dj;
                            let nk = k as i32 + dk;
                            if ni >= 0
                                && ni < nx as i32
                                && nj >= 0
                                && nj < ny as i32
                                && nk >= 0
                                && nk < nz as i32
                            {
                                let w =
                                    kernel[(di + 1) as usize][(dj + 1) as usize][(dk + 1) as usize];
                                acc += w * grid.get(ni as usize, nj as usize, nk as usize);
                                wt += w;
                            }
                        }
                    }
                }
                let v = if wt > 1e-15 {
                    acc / wt
                } else {
                    grid.get(i, j, k)
                };
                out.set(i, j, k, v);
            }
        }
    }
    let _ = kernel_sum;
    out
}
/// `true` when the grid dimensions are consistent with its storage, it is
/// non-empty and its cell size is positive and finite.
pub(super) fn grid_is_usable(grid: &SdfGrid) -> bool {
    let cells = grid
        .nx
        .checked_mul(grid.ny)
        .and_then(|plane| plane.checked_mul(grid.nz));
    cells == Some(grid.values.len())
        && !grid.values.is_empty()
        && grid.dx.is_finite()
        && grid.dx > 0.0
}
/// Copy of `grid` with its dimensions, spacing and origin but new `values`.
fn with_values(grid: &SdfGrid, values: Vec<f64>) -> SdfGrid {
    SdfGrid {
        nx: grid.nx,
        ny: grid.ny,
        nz: grid.nz,
        dx: grid.dx,
        origin: grid.origin,
        values,
    }
}
/// One explicit step `dst = src + coeff * D(src)` where `D` is the sum over the
/// three axes of the second difference `v₊ − 2v + v₋`.
///
/// `coeff` therefore already contains the `1/dx²` of the Laplacian.  Along an axis
/// on which a cell lacks a neighbour on either side the second difference is
/// taken as zero (a linear extrapolation boundary condition).  Rows of constant
/// `i` are processed in parallel.
fn laplacian_step(src: &[f64], dst: &mut [f64], dims: [usize; 3], coeff: f64) {
    let [nx, ny, nz] = dims;
    let stride_i = ny * nz;
    dst.par_chunks_mut(stride_i)
        .enumerate()
        .for_each(|(i, slab)| {
            for j in 0..ny {
                for k in 0..nz {
                    let idx = i * stride_i + j * nz + k;
                    let v = src[idx];
                    let mut second = 0.0;
                    if i > 0 && i + 1 < nx {
                        second += src[idx + stride_i] - 2.0 * v + src[idx - stride_i];
                    }
                    if j > 0 && j + 1 < ny {
                        second += src[idx + nz] - 2.0 * v + src[idx - nz];
                    }
                    if k > 0 && k + 1 < nz {
                        second += src[idx + 1] - 2.0 * v + src[idx - 1];
                    }
                    slab[j * nz + k] = v + coeff * second;
                }
            }
        });
}
/// Laplacian sharpening (unsharp mask) of the SDF.
///
/// Returns `grid - amount * Δgrid`, where `Δ` is the 7-point Laplacian
/// (`Σ (v₊ − 2v + v₋) / dx²`, zero along axes on which a cell has no two
/// neighbours).  Subtracting the Laplacian is anti-diffusion: for `amount > 0`
/// local extrema are exaggerated (a negative spike becomes more negative) and
/// grid-scale noise is amplified by up to `1 + 12·amount/dx²`, while a negative
/// `amount` smooths exactly like [`sdf_laplacian_smooth`] without sub-stepping.
/// A grid with inconsistent dimensions or a non-positive cell size is returned
/// unchanged.
pub fn sdf_laplacian_sharpen(grid: &SdfGrid, amount: f64) -> SdfGrid {
    if !grid_is_usable(grid) {
        return with_values(grid, grid.values.clone());
    }
    let mut values = vec![0.0; grid.values.len()];
    let coeff = -amount / (grid.dx * grid.dx);
    laplacian_step(
        &grid.values,
        &mut values,
        [grid.nx, grid.ny, grid.nz],
        coeff,
    );
    with_values(grid, values)
}
/// SDF dilation: offset the surface outward by `offset` units.
///
/// Subtracts `offset` from all SDF values, which moves the zero level set
/// outward along the gradient; a negative offset erodes.  Because only the values
/// are shifted, the result stays an exact distance field on the outside of the
/// shape but not wherever the offset reaches the medial axis; run
/// [`sdf_redistance`] to restore the unit-gradient
/// property.
pub fn sdf_dilate(grid: &SdfGrid, offset: f64) -> SdfGrid {
    let values: Vec<f64> = grid.values.par_iter().map(|&v| v - offset).collect();
    SdfGrid {
        nx: grid.nx,
        ny: grid.ny,
        nz: grid.nz,
        dx: grid.dx,
        origin: grid.origin,
        values,
    }
}
/// SDF erosion: offset the surface inward by `offset` units.
///
/// Equivalent to `sdf_dilate(grid, -offset)`.
pub fn sdf_erode(grid: &SdfGrid, offset: f64) -> SdfGrid {
    sdf_dilate(grid, -offset)
}
/// SDF morphological opening: erode, re-initialise the distance field, dilate.
///
/// The solid is eroded by `offset`, so every part thinner than `2 * offset`
/// vanishes; the eroded field is then re-initialised with
/// [`sdf_redistance`], because a uniform shift keeps the
/// values that belonged to the removed parts and dilating that field would simply
/// grow them back; finally the remainder is dilated by `offset`.  The result
/// removes protrusions and islands narrower than `2 * offset`, keeps larger
/// features (rounding their convex corners) and is a subset of the input, up to
/// the first-order accuracy of the re-initialisation (a fraction of a cell).  A
/// shape that is entirely removed by the erosion yields an empty field (`+f64::MAX`
/// everywhere).  A negative `offset` exchanges the two operations, i.e. computes
/// [`sdf_close`] for `-offset`.
pub fn sdf_open(grid: &SdfGrid, offset: f64) -> SdfGrid {
    let eroded = sdf_erode(grid, offset);
    let redistanced = sdf_redistance(&eroded);
    sdf_dilate(&redistanced, offset)
}
/// SDF morphological closing: dilate, re-initialise the distance field, erode.
///
/// The solid is dilated by `offset`, so gaps and holes narrower than
/// `2 * offset` are filled and nearby parts merge; the dilated field is
/// re-initialised with [`sdf_redistance`] so that the
/// following erosion starts from the true distance to the merged shape instead of
/// undoing the dilation exactly; finally the result is eroded by `offset`.  The
/// result fills small holes and gaps, keeps the rest of the boundary and is a
/// superset of the input, up to the first-order accuracy of the re-initialisation
/// (a fraction of a cell).  A negative `offset` exchanges the two operations, i.e.
/// computes [`sdf_open`] for `-offset`.
pub fn sdf_close(grid: &SdfGrid, offset: f64) -> SdfGrid {
    let dilated = sdf_dilate(grid, offset);
    let redistanced = sdf_redistance(&dilated);
    sdf_erode(&redistanced, offset)
}
/// Signed distance field offset: produce an iso-surface at `offset`.
///
/// The new surface is the locus of points where the original SDF = `offset`.
pub fn sdf_offset_surface(grid: &SdfGrid, offset: f64) -> SdfGrid {
    let values: Vec<f64> = grid.values.par_iter().map(|&v| v - offset).collect();
    SdfGrid {
        nx: grid.nx,
        ny: grid.ny,
        nz: grid.nz,
        dx: grid.dx,
        origin: grid.origin,
        values,
    }
}
/// Iterative Laplacian smoothing of the SDF.
///
/// Advances the heat equation `∂φ/∂t = Δφ` by `dt` pseudo-time units in each of
/// `n_iterations` iterations (small `dt` → mild smoothing), with the explicit
/// update `v ← v + dt_sub·Δv`.  The explicit 3-D scheme is only stable for
/// `dt_sub / dx² ≤ 1/6`, so every iteration is split into
/// `m = ⌈6·dt/dx²⌉` equal sub-steps of `dt_sub = dt/m`; at that bound each update
/// is a convex combination of a cell and its six neighbours, hence smoothing never
/// creates new extrema (`max|φ|` does not increase) and noise decays instead of
/// growing.  The cost is `n_iterations · m` sweeps, i.e. proportional to
/// `dt/dx²`.  Along an axis on which a cell lacks neighbours on both sides the
/// second difference is taken as zero, so linear fields are preserved.
///
/// A grid with inconsistent dimensions, a non-positive or non-finite `dt` or
/// cell size, or `n_iterations == 0` is returned unchanged.
pub fn sdf_laplacian_smooth(grid: &SdfGrid, n_iterations: usize, dt: f64) -> SdfGrid {
    let ratio = dt / (grid.dx * grid.dx);
    if n_iterations == 0
        || !grid_is_usable(grid)
        || !dt.is_finite()
        || dt <= 0.0
        || !ratio.is_finite()
        || ratio <= 0.0
    {
        return with_values(grid, grid.values.clone());
    }
    let mut substeps = (6.0 * ratio).ceil().max(1.0) as usize;
    while ratio / (substeps as f64) * 6.0 > 1.0 {
        substeps += 1;
    }
    let coeff = ratio / substeps as f64;
    let dims = [grid.nx, grid.ny, grid.nz];
    let mut current = grid.values.clone();
    let mut next = vec![0.0; current.len()];
    for _ in 0..n_iterations.saturating_mul(substeps) {
        laplacian_step(&current, &mut next, dims, coeff);
        std::mem::swap(&mut current, &mut next);
    }
    with_values(grid, current)
}
/// Mean-curvature smoothing of the SDF.
///
/// For a signed distance field `|∇φ| = 1`, so the level-set mean-curvature flow
/// `∂φ/∂t = |∇φ|·∇·(∇φ/|∇φ|)` reduces to the heat equation `∂φ/∂t = Δφ` (the
/// Laplacian of a distance function is the sum of the principal curvatures of its
/// level sets).  This function therefore advances that flow by `step` pseudo-time
/// units with one stable [`sdf_laplacian_smooth`] iteration.  The approximation
/// holds only while the field is close to a distance function and the smoothed
/// field is not renormalised; call [`sdf_redistance`]
/// afterwards when an exact unit-gradient field is needed.
pub fn sdf_mean_curvature_smooth(grid: &SdfGrid, step: f64) -> SdfGrid {
    sdf_laplacian_smooth(grid, 1, step)
}
#[cfg(test)]
mod tests {
    use super::*;
    fn make_sphere_grid(nx: usize, dx: f64, center: [f64; 3], radius: f64) -> SdfGrid {
        let origin = [0.0; 3];
        let mut g = SdfGrid::new(nx, nx, nx, dx, origin);
        g.compute_sphere_sdf(center, radius);
        g
    }
    #[test]
    fn test_sphere_center_is_negative_radius() {
        let n = 21usize;
        let dx = 0.1;
        let radius = 0.4;
        let mid = (n / 2) as f64 + 0.5;
        let center = [mid * dx, mid * dx, mid * dx];
        let g = make_sphere_grid(n, dx, center, radius);
        let c = n / 2;
        let sdf_val = g.get(c, c, c);
        assert!(
            (sdf_val - (-radius)).abs() < dx,
            "centre value {sdf_val} should be close to -{radius}"
        );
    }
    #[test]
    fn test_box_far_outside_positive() {
        let n = 11usize;
        let dx = 0.1;
        let mut g = SdfGrid::new(n, n, n, dx, [0.0; 3]);
        let box_center = [0.55, 0.55, 0.55];
        let half_extents = [0.1, 0.1, 0.1];
        g.compute_box_sdf(box_center, half_extents);
        let v = g.get(0, 0, 0);
        assert!(
            v > 0.0,
            "far-outside cell should have positive SDF, got {v}"
        );
    }
    #[test]
    fn test_gradient_on_sphere_surface_is_unit() {
        let n = 41usize;
        let dx = 0.05;
        let radius = 0.5;
        let mid = (n / 2) as f64 + 0.5;
        let center = [mid * dx, mid * dx, mid * dx];
        let g = make_sphere_grid(n, dx, center, radius);
        let c = n / 2;
        let surface_i = c + (radius / dx) as usize;
        let grad = g.gradient_at(surface_i, c, c);
        let mag = (grad[0].powi(2) + grad[1].powi(2) + grad[2].powi(2)).sqrt();
        assert!(
            (mag - 1.0).abs() < 0.1,
            "gradient magnitude should be close to 1.0, got {mag}"
        );
        assert!(grad[0] > 0.5, "gradient should point outward, got {grad:?}");
    }
    #[test]
    fn test_union_sdf_inside_either() {
        let n = 21usize;
        let dx = 0.1;
        let radius = 0.3;
        let origin = [0.0; 3];
        let c_float = (n / 2) as f64 + 0.5;
        let center_a = [(c_float - 3.0) * dx, c_float * dx, c_float * dx];
        let center_b = [(c_float + 3.0) * dx, c_float * dx, c_float * dx];
        let mut ga = SdfGrid::new(n, n, n, dx, origin);
        ga.compute_sphere_sdf(center_a, radius);
        let mut gb = SdfGrid::new(n, n, n, dx, origin);
        gb.compute_sphere_sdf(center_b, radius);
        let u = union_sdf(&ga, &gb);
        let cy = n / 2;
        let cz = n / 2;
        let ia = n / 2 - 3;
        assert!(
            u.get(ia, cy, cz) < 0.0,
            "inside sphere A should be negative in union"
        );
        let ib = n / 2 + 3;
        assert!(
            u.get(ib, cy, cz) < 0.0,
            "inside sphere B should be negative in union"
        );
    }
    #[test]
    fn test_total_cells() {
        let g = SdfGrid::new(4, 5, 6, 0.1, [0.0; 3]);
        assert_eq!(g.total_cells(), 4 * 5 * 6);
    }
    /// Subtracting a large sphere from a small sphere should produce
    /// positive values everywhere.
    #[test]
    fn test_subtraction_sdf() {
        let n = 11usize;
        let dx = 0.2;
        let origin = [0.0; 3];
        let center = [1.1, 1.1, 1.1];
        let mut small = SdfGrid::new(n, n, n, dx, origin);
        small.compute_sphere_sdf(center, 0.3);
        let mut large = SdfGrid::new(n, n, n, dx, origin);
        large.compute_sphere_sdf(center, 0.5);
        let result = subtraction_sdf(&small, &large);
        let c = n / 2;
        assert!(
            result.get(c, c, c) > 0.0,
            "subtraction centre should be positive, got {}",
            result.get(c, c, c)
        );
    }
    /// Intersection of two overlapping spheres should be smaller than either.
    #[test]
    fn test_intersection_sdf() {
        let n = 21usize;
        let dx = 0.1;
        let origin = [0.0; 3];
        let radius = 0.5;
        let c = (n / 2) as f64 + 0.5;
        let center_a = [(c - 1.0) * dx, c * dx, c * dx];
        let center_b = [(c + 1.0) * dx, c * dx, c * dx];
        let mut ga = SdfGrid::new(n, n, n, dx, origin);
        ga.compute_sphere_sdf(center_a, radius);
        let mut gb = SdfGrid::new(n, n, n, dx, origin);
        gb.compute_sphere_sdf(center_b, radius);
        let inter = intersection_sdf(&ga, &gb);
        let mid = n / 2;
        let val = inter.get(mid, mid, mid);
        assert!(
            val < 0.0,
            "midpoint of intersection should be inside, got {val}"
        );
    }
    /// Smooth union should produce smaller values than min at the seam.
    #[test]
    fn test_smooth_union() {
        let n = 11usize;
        let dx = 0.2;
        let origin = [0.0; 3];
        let radius = 0.3;
        let c = (n / 2) as f64 + 0.5;
        let center_a = [(c - 2.0) * dx, c * dx, c * dx];
        let center_b = [(c + 2.0) * dx, c * dx, c * dx];
        let mut ga = SdfGrid::new(n, n, n, dx, origin);
        ga.compute_sphere_sdf(center_a, radius);
        let mut gb = SdfGrid::new(n, n, n, dx, origin);
        gb.compute_sphere_sdf(center_b, radius);
        let su = smooth_union_sdf(&ga, &gb, 0.5);
        let u = union_sdf(&ga, &gb);
        let mid = n / 2;
        assert!(
            su.get(mid, mid, mid) <= u.get(mid, mid, mid) + 0.1,
            "smooth union should not be much larger than union"
        );
    }
    /// Sample at a cell centre should match the cell value.
    #[test]
    fn test_sample_at_cell_center() {
        let n = 11usize;
        let dx = 0.1;
        let mut g = SdfGrid::new(n, n, n, dx, [0.0; 3]);
        g.compute_sphere_sdf([0.55, 0.55, 0.55], 0.3);
        let c = n / 2;
        let pos = g.world_pos(c, c, c);
        let sampled = g.sample(pos);
        assert!(sampled.is_some());
        let cell_val = g.get(c, c, c);
        assert!(
            (sampled.unwrap() - cell_val).abs() < 0.05,
            "sampled = {:?}, cell = {cell_val}",
            sampled
        );
    }
    /// Sample outside the grid should return None.
    #[test]
    fn test_sample_outside_grid() {
        let g = SdfGrid::new(5, 5, 5, 0.1, [0.0; 3]);
        assert!(g.sample([-1.0, -1.0, -1.0]).is_none());
    }
    /// Ray pointing toward a sphere should hit.
    #[test]
    fn test_sphere_trace_hit() {
        let n = 41usize;
        let dx = 0.05;
        let center = [1.0, 1.0, 1.0];
        let radius = 0.3;
        let mut g = SdfGrid::new(n, n, n, dx, [0.0; 3]);
        g.compute_sphere_sdf(center, radius);
        let result = sphere_trace(&g, [0.1, 1.0, 1.0], [1.0, 0.0, 0.0], 5.0, 100, dx * 0.5);
        assert!(result.hit, "ray should hit the sphere");
        assert!(result.t > 0.0, "t should be positive");
    }
    /// Ray pointing away from a sphere should miss.
    #[test]
    fn test_sphere_trace_miss() {
        let n = 21usize;
        let dx = 0.1;
        let center = [1.0, 1.0, 1.0];
        let radius = 0.3;
        let mut g = SdfGrid::new(n, n, n, dx, [0.0; 3]);
        g.compute_sphere_sdf(center, radius);
        let result = sphere_trace(&g, [0.1, 1.0, 1.0], [-1.0, 0.0, 0.0], 5.0, 100, dx * 0.5);
        assert!(!result.hit, "ray should miss the sphere");
    }
    /// Distance from a point to a triangle vertex.
    #[test]
    fn test_point_triangle_distance_at_vertex() {
        let a = [0.0, 0.0, 0.0];
        let b = [1.0, 0.0, 0.0];
        let c = [0.0, 1.0, 0.0];
        let p = [-1.0, 0.0, 0.0];
        let d = point_triangle_distance(&p, &a, &b, &c);
        assert!((d - 1.0).abs() < 1e-10, "distance = {d}, expected 1.0");
    }
    /// Distance from a point directly above the triangle center.
    #[test]
    fn test_point_triangle_distance_above() {
        let a = [0.0, 0.0, 0.0];
        let b = [1.0, 0.0, 0.0];
        let c = [0.0, 1.0, 0.0];
        let p = [0.2, 0.2, 1.0];
        let d = point_triangle_distance(&p, &a, &b, &c);
        assert!((d - 1.0).abs() < 1e-10, "distance = {d}, expected 1.0");
    }
    /// Batch evaluation should return one value per point.
    #[test]
    fn test_evaluate_sdf_batch() {
        let n = 11usize;
        let dx = 0.2;
        let mut g = SdfGrid::new(n, n, n, dx, [0.0; 3]);
        g.compute_sphere_sdf([1.1, 1.1, 1.1], 0.5);
        let points = vec![[1.1, 1.1, 1.1], [0.1, 0.1, 0.1]];
        let vals = evaluate_sdf_batch(&g, &points);
        assert_eq!(vals.len(), 2);
        assert!(vals[0] < 0.0, "centre should be negative, got {}", vals[0]);
    }
    /// Volume of a sphere should be approximately 4/3 * pi * r^3.
    #[test]
    fn test_approximate_volume_sphere() {
        let n = 41usize;
        let dx = 0.05;
        let radius = 0.5;
        let center = [1.0, 1.0, 1.0];
        let mut g = SdfGrid::new(n, n, n, dx, [0.0; 3]);
        g.compute_sphere_sdf(center, radius);
        let vol = approximate_volume(&g);
        let expected = 4.0 / 3.0 * std::f64::consts::PI * radius.powi(3);
        assert!(
            (vol - expected).abs() / expected < 0.2,
            "volume = {vol}, expected ~{expected}"
        );
    }
    /// Cylinder SDF: centre should be inside.
    #[test]
    fn test_cylinder_sdf_center_inside() {
        let n = 21usize;
        let dx = 0.1;
        let mut g = SdfGrid::new(n, n, n, dx, [0.0; 3]);
        g.compute_cylinder_sdf([1.0, 1.0], 0.3, 1.2);
        let c = n / 2;
        let val = g.get(10, 10, c);
        assert!(val < 0.0, "cylinder centre should be inside, got {val}");
    }
    /// Torus SDF: ring centre should be inside.
    #[test]
    fn test_torus_sdf() {
        let n = 21usize;
        let dx = 0.1;
        let center = [1.0, 1.0, 1.0];
        let mut g = SdfGrid::new(n, n, n, dx, [0.0; 3]);
        g.compute_torus_sdf(center, 0.4, 0.1);
        let ring_x = ((center[0] + 0.4) / dx - 0.5) as usize;
        let ring_y = (center[1] / dx - 0.5) as usize;
        let ring_z = (center[2] / dx - 0.5) as usize;
        if ring_x < n && ring_y < n && ring_z < n {
            let val = g.get(ring_x, ring_y, ring_z);
            assert!(val < 0.1, "ring point should be near surface, got {val}");
        }
    }
}
