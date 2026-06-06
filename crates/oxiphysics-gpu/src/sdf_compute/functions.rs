//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use rayon::prelude::*;

use super::types::{SdfGrid, SphereTraceResult, Triangle};

/// One pass of the fast sweeping method.
pub fn fast_sweeping_update(grid: &mut SdfGrid) {
    let nx = grid.nx;
    let ny = grid.ny;
    let nz = grid.nz;
    let dx = grid.dx;
    for i in 0..nx {
        for j in 0..ny {
            for k in 0..nz {
                let cur = grid.get(i, j, k);
                let mut best = cur;
                if i > 0 {
                    let candidate = grid.get(i - 1, j, k) + dx;
                    if candidate < best {
                        best = candidate;
                    }
                }
                if i + 1 < nx {
                    let candidate = grid.get(i + 1, j, k) + dx;
                    if candidate < best {
                        best = candidate;
                    }
                }
                if j > 0 {
                    let candidate = grid.get(i, j - 1, k) + dx;
                    if candidate < best {
                        best = candidate;
                    }
                }
                if j + 1 < ny {
                    let candidate = grid.get(i, j + 1, k) + dx;
                    if candidate < best {
                        best = candidate;
                    }
                }
                if k > 0 {
                    let candidate = grid.get(i, j, k - 1) + dx;
                    if candidate < best {
                        best = candidate;
                    }
                }
                if k + 1 < nz {
                    let candidate = grid.get(i, j, k + 1) + dx;
                    if candidate < best {
                        best = candidate;
                    }
                }
                if best < cur {
                    grid.set(i, j, k, best);
                }
            }
        }
    }
}
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
/// Marching cubes edge table (12 edges per cube).
/// Each entry encodes which edges are intersected for a given vertex mask.
#[rustfmt::skip]
pub(super) const MC_EDGE_TABLE: [u16; 256] = [
    0x000, 0x109, 0x203, 0x30a, 0x406, 0x50f, 0x605, 0x70c, 0x80c, 0x905, 0xa0f, 0xb06,
    0xc0a, 0xd03, 0xe09, 0xf00, 0x190, 0x099, 0x393, 0x29a, 0x596, 0x49f, 0x795, 0x69c,
    0x99c, 0x895, 0xb9f, 0xa96, 0xd9a, 0xc93, 0xf99, 0xe90, 0x230, 0x339, 0x033, 0x13a,
    0x636, 0x73f, 0x435, 0x53c, 0xa3c, 0xb35, 0x83f, 0x936, 0xe3a, 0xf33, 0xc39, 0xd30,
    0x3a0, 0x2a9, 0x1a3, 0x0aa, 0x7a6, 0x6af, 0x5a5, 0x4ac, 0xbac, 0xaa5, 0x9af, 0x8a6,
    0xfaa, 0xea3, 0xda9, 0xca0, 0x460, 0x569, 0x663, 0x76a, 0x066, 0x16f, 0x265, 0x36c,
    0xc6c, 0xd65, 0xe6f, 0xf66, 0x86a, 0x963, 0xa69, 0xb60, 0x5f0, 0x4f9, 0x7f3, 0x6fa,
    0x1f6, 0x0ff, 0x3f5, 0x2fc, 0xdfc, 0xcf5, 0xfff, 0xef6, 0x9fa, 0x8f3, 0xbf9, 0xaf0,
    0x650, 0x759, 0x453, 0x55a, 0x256, 0x35f, 0x055, 0x15c, 0xe5c, 0xf55, 0xc5f, 0xd56,
    0xa5a, 0xb53, 0x859, 0x950, 0x7c0, 0x6c9, 0x5c3, 0x4ca, 0x3c6, 0x2cf, 0x1c5, 0x0cc,
    0xfcc, 0xec5, 0xdcf, 0xcc6, 0xbca, 0xac3, 0x9c9, 0x8c0, 0x8c0, 0x9c9, 0xac3, 0xbca,
    0xcc6, 0xdcf, 0xec5, 0xfcc, 0x0cc, 0x1c5, 0x2cf, 0x3c6, 0x4ca, 0x5c3, 0x6c9, 0x7c0,
    0x950, 0x859, 0xb53, 0xa5a, 0xd56, 0xc5f, 0xf55, 0xe5c, 0x15c, 0x055, 0x35f, 0x256,
    0x55a, 0x453, 0x759, 0x650, 0xaf0, 0xbf9, 0x8f3, 0x9fa, 0xef6, 0xfff, 0xcf5, 0xdfc,
    0x2fc, 0x3f5, 0x0ff, 0x1f6, 0x6fa, 0x7f3, 0x4f9, 0x5f0, 0xb60, 0xa69, 0x963, 0x86a,
    0xf66, 0xe6f, 0xd65, 0xc6c, 0x36c, 0x265, 0x16f, 0x066, 0x76a, 0x663, 0x569, 0x460,
    0xca0, 0xda9, 0xea3, 0xfaa, 0x8a6, 0x9af, 0xaa5, 0xbac, 0x4ac, 0x5a5, 0x6af, 0x7a6,
    0x0aa, 0x1a3, 0x2a9, 0x3a0, 0xd30, 0xc39, 0xf33, 0xe3a, 0x936, 0x835, 0xb3f, 0xa36,
    0x53c, 0x435, 0x73f, 0x636, 0x13a, 0x033, 0x339, 0x230, 0xe90, 0xf99, 0xc93, 0xd9a,
    0xa96, 0xb9f, 0x895, 0x99c, 0x69c, 0x795, 0x49f, 0x596, 0x29a, 0x393, 0x099, 0x190,
    0xf00, 0xe09, 0xd03, 0xc0a, 0xb06, 0xa0f, 0x905, 0x80c, 0x70c, 0x605, 0x50f, 0x406,
    0x30a, 0x203, 0x109, 0x000,
];
/// Interpolate the position of an edge intersection given two corner SDF values.
#[inline]
pub(super) fn interpolate_vertex(
    p1: [f64; 3],
    p2: [f64; 3],
    val1: f64,
    val2: f64,
    iso: f64,
) -> [f64; 3] {
    if (val2 - val1).abs() < 1e-15 {
        return p1;
    }
    let t = (iso - val1) / (val2 - val1);
    [
        p1[0] + t * (p2[0] - p1[0]),
        p1[1] + t * (p2[1] - p1[1]),
        p1[2] + t * (p2[2] - p1[2]),
    ]
}
/// Extract a triangle mesh from the SDF grid at the given isovalue
/// using the marching cubes algorithm (simplified 3-case variant).
///
/// Returns a list of triangles. For a signed distance field, `isovalue = 0.0`
/// extracts the zero level-set (the surface).
pub fn marching_cubes(grid: &SdfGrid, isovalue: f64) -> Vec<Triangle> {
    let nx = grid.nx;
    let ny = grid.ny;
    let nz = grid.nz;
    if nx < 2 || ny < 2 || nz < 2 {
        return Vec::new();
    }
    let mut triangles = Vec::new();
    for i in 0..nx - 1 {
        for j in 0..ny - 1 {
            for k in 0..nz - 1 {
                let corners: [[f64; 3]; 8] = [
                    grid.world_pos(i, j, k),
                    grid.world_pos(i + 1, j, k),
                    grid.world_pos(i + 1, j + 1, k),
                    grid.world_pos(i, j + 1, k),
                    grid.world_pos(i, j, k + 1),
                    grid.world_pos(i + 1, j, k + 1),
                    grid.world_pos(i + 1, j + 1, k + 1),
                    grid.world_pos(i, j + 1, k + 1),
                ];
                let vals: [f64; 8] = [
                    grid.get(i, j, k),
                    grid.get(i + 1, j, k),
                    grid.get(i + 1, j + 1, k),
                    grid.get(i, j + 1, k),
                    grid.get(i, j, k + 1),
                    grid.get(i + 1, j, k + 1),
                    grid.get(i + 1, j + 1, k + 1),
                    grid.get(i, j + 1, k + 1),
                ];
                let mut cube_idx = 0u8;
                for (c, &v) in vals.iter().enumerate() {
                    if v < isovalue {
                        cube_idx |= 1 << c;
                    }
                }
                let edge_flags = MC_EDGE_TABLE[cube_idx as usize];
                if edge_flags == 0 {
                    continue;
                }
                let mut verts = [[0.0f64; 3]; 12];
                if edge_flags & 0x001 != 0 {
                    verts[0] =
                        interpolate_vertex(corners[0], corners[1], vals[0], vals[1], isovalue);
                }
                if edge_flags & 0x002 != 0 {
                    verts[1] =
                        interpolate_vertex(corners[1], corners[2], vals[1], vals[2], isovalue);
                }
                if edge_flags & 0x004 != 0 {
                    verts[2] =
                        interpolate_vertex(corners[2], corners[3], vals[2], vals[3], isovalue);
                }
                if edge_flags & 0x008 != 0 {
                    verts[3] =
                        interpolate_vertex(corners[3], corners[0], vals[3], vals[0], isovalue);
                }
                if edge_flags & 0x010 != 0 {
                    verts[4] =
                        interpolate_vertex(corners[4], corners[5], vals[4], vals[5], isovalue);
                }
                if edge_flags & 0x020 != 0 {
                    verts[5] =
                        interpolate_vertex(corners[5], corners[6], vals[5], vals[6], isovalue);
                }
                if edge_flags & 0x040 != 0 {
                    verts[6] =
                        interpolate_vertex(corners[6], corners[7], vals[6], vals[7], isovalue);
                }
                if edge_flags & 0x080 != 0 {
                    verts[7] =
                        interpolate_vertex(corners[7], corners[4], vals[7], vals[4], isovalue);
                }
                if edge_flags & 0x100 != 0 {
                    verts[8] =
                        interpolate_vertex(corners[0], corners[4], vals[0], vals[4], isovalue);
                }
                if edge_flags & 0x200 != 0 {
                    verts[9] =
                        interpolate_vertex(corners[1], corners[5], vals[1], vals[5], isovalue);
                }
                if edge_flags & 0x400 != 0 {
                    verts[10] =
                        interpolate_vertex(corners[2], corners[6], vals[2], vals[6], isovalue);
                }
                if edge_flags & 0x800 != 0 {
                    verts[11] =
                        interpolate_vertex(corners[3], corners[7], vals[3], vals[7], isovalue);
                }
                let active: Vec<[f64; 3]> = (0..12)
                    .filter(|&e| edge_flags & (1 << e) != 0)
                    .map(|e| verts[e])
                    .collect();
                if active.len() >= 3 {
                    for tri_idx in 1..active.len() - 1 {
                        triangles.push(Triangle {
                            v: [active[0], active[tri_idx], active[tri_idx + 1]],
                        });
                    }
                }
            }
        }
    }
    triangles
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
/// Laplacian (sharpening) convolution on the SDF.
///
/// Returns grid + `amount` * Laplacian(grid), which sharpens features.
pub fn sdf_laplacian_sharpen(grid: &SdfGrid, amount: f64) -> SdfGrid {
    let nx = grid.nx;
    let ny = grid.ny;
    let nz = grid.nz;
    let inv_dx2 = 1.0 / (grid.dx * grid.dx);
    let mut out = SdfGrid::new(nx, ny, nz, grid.dx, grid.origin);
    for i in 0..nx {
        for j in 0..ny {
            for k in 0..nz {
                let v = grid.get(i, j, k);
                let lx = if i > 0 && i + 1 < nx {
                    (grid.get(i + 1, j, k) - 2.0 * v + grid.get(i - 1, j, k)) * inv_dx2
                } else {
                    0.0
                };
                let ly = if j > 0 && j + 1 < ny {
                    (grid.get(i, j + 1, k) - 2.0 * v + grid.get(i, j - 1, k)) * inv_dx2
                } else {
                    0.0
                };
                let lz = if k > 0 && k + 1 < nz {
                    (grid.get(i, j, k + 1) - 2.0 * v + grid.get(i, j, k - 1)) * inv_dx2
                } else {
                    0.0
                };
                out.set(i, j, k, v + amount * (lx + ly + lz));
            }
        }
    }
    out
}
/// SDF dilation: offset the surface outward by `offset` units.
///
/// Simply subtracts `offset` from all SDF values.
/// Positive offset = dilation (expand solid), negative = erosion.
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
/// SDF morphological opening: erode then dilate.
///
/// Removes small protrusions (blobs smaller than `offset`).
pub fn sdf_open(grid: &SdfGrid, offset: f64) -> SdfGrid {
    let eroded = sdf_erode(grid, offset);
    sdf_dilate(&eroded, offset)
}
/// SDF morphological closing: dilate then erode.
///
/// Fills small holes (gaps smaller than `offset`).
pub fn sdf_close(grid: &SdfGrid, offset: f64) -> SdfGrid {
    let dilated = sdf_dilate(grid, offset);
    sdf_erode(&dilated, offset)
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
/// Applies `n_iterations` of a simple diffusion smoother.
/// `dt` is the pseudo-time step (small values → mild smoothing).
pub fn sdf_laplacian_smooth(grid: &SdfGrid, n_iterations: usize, dt: f64) -> SdfGrid {
    let mut current = SdfGrid {
        nx: grid.nx,
        ny: grid.ny,
        nz: grid.nz,
        dx: grid.dx,
        origin: grid.origin,
        values: grid.values.clone(),
    };
    for _ in 0..n_iterations {
        let sharpened = sdf_laplacian_sharpen(&current, -dt);
        current = sharpened;
    }
    current
}
/// Mean-curvature smoothing of the SDF (simplified).
///
/// Uses the divergence of the gradient (Laplacian as proxy for mean curvature).
/// `step` is the smoothing step size.
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
