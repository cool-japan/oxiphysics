//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)
use rayon::prelude::*;

use super::types::{DistanceQuery, GpuSdfGrid, SdfGrid, SdfShape};

/// Bilateral filter for SDF: preserves sharp features while smoothing noise.
///
/// Uses a spatial Gaussian kernel weighted by the range (SDF value) difference.
pub fn sdf_bilateral_filter(grid: &SdfGrid, sigma_s: f64, sigma_r: f64) -> SdfGrid {
    let nx = grid.nx;
    let ny = grid.ny;
    let nz = grid.nz;
    let mut out = SdfGrid::new(nx, ny, nz, grid.dx, grid.origin);
    let s2 = 2.0 * sigma_s * sigma_s;
    let r2 = 2.0 * sigma_r * sigma_r;
    for i in 0..nx {
        for j in 0..ny {
            for k in 0..nz {
                let v0 = grid.get(i, j, k);
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
                                let vn = grid.get(ni as usize, nj as usize, nk as usize);
                                let dist2 = (di * di + dj * dj + dk * dk) as f64;
                                let w_s = (-dist2 / s2).exp();
                                let w_r = (-(v0 - vn) * (v0 - vn) / r2).exp();
                                let w = w_s * w_r;
                                acc += w * vn;
                                wt += w;
                            }
                        }
                    }
                }
                out.set(i, j, k, if wt > 1e-15 { acc / wt } else { v0 });
            }
        }
    }
    out
}
/// Query the distance field at a point and compute surface normal + closest point.
pub fn query_distance_field(grid: &SdfGrid, pos: [f64; 3]) -> Option<DistanceQuery> {
    let dist = grid.sample(pos)?;
    let grad = grid.gradient_at_point(pos).unwrap_or([0.0; 3]);
    let grad_len = (grad[0] * grad[0] + grad[1] * grad[1] + grad[2] * grad[2]).sqrt();
    let normal = if grad_len > 1e-15 {
        [grad[0] / grad_len, grad[1] / grad_len, grad[2] / grad_len]
    } else {
        [0.0, 0.0, 1.0]
    };
    let closest_point = [
        pos[0] - dist * normal[0],
        pos[1] - dist * normal[1],
        pos[2] - dist * normal[2],
    ];
    Some(DistanceQuery {
        distance: dist,
        normal,
        closest_point,
        is_inside: dist < 0.0,
    })
}
/// Batch query multiple points and return DistanceQuery results.
pub fn query_distance_field_batch(
    grid: &SdfGrid,
    points: &[[f64; 3]],
) -> Vec<Option<DistanceQuery>> {
    points
        .par_iter()
        .map(|&p| query_distance_field(grid, p))
        .collect()
}
/// Find the zero-crossing (surface) along a 1-D ray by bisection.
///
/// Returns the parameter `t` such that `grid.sample(origin + t * direction) ≈ 0`,
/// or `None` if no crossing found in `[t_min, t_max]`.
pub fn find_zero_crossing(
    grid: &SdfGrid,
    origin: [f64; 3],
    direction: [f64; 3],
    t_min: f64,
    t_max: f64,
    n_bisect: usize,
) -> Option<f64> {
    let sample_at = |t: f64| -> Option<f64> {
        let p = [
            origin[0] + t * direction[0],
            origin[1] + t * direction[1],
            origin[2] + t * direction[2],
        ];
        grid.sample(p)
    };
    let v_min = sample_at(t_min)?;
    let v_max = sample_at(t_max)?;
    if v_min * v_max > 0.0 {
        return None;
    }
    let mut lo = t_min;
    let mut hi = t_max;
    let mut v_lo = v_min;
    for _ in 0..n_bisect {
        let mid = (lo + hi) * 0.5;
        let v_mid = sample_at(mid)?;
        if v_lo * v_mid <= 0.0 {
            hi = mid;
        } else {
            lo = mid;
            v_lo = v_mid;
        }
    }
    Some((lo + hi) * 0.5)
}
/// Compute the projected area of a surface onto the xy-plane.
///
/// Counts grid cells with SDF < 0 in the bottom z-slice.
pub fn projected_area_xy(grid: &SdfGrid) -> f64 {
    let ny = grid.ny;
    let nz = grid.nz;
    let nx = grid.nx;
    let mut count = 0usize;
    for i in 0..nx {
        for j in 0..ny {
            let occupied = (0..nz).any(|k| grid.get(i, j, k) < 0.0);
            if occupied {
                count += 1;
            }
        }
    }
    count as f64 * grid.dx * grid.dx
}
#[cfg(test)]
mod tests_new_sdf {
    use super::super::functions::*;
    use super::*;

    fn sphere_grid(n: usize, dx: f64, radius: f64) -> SdfGrid {
        let center = [(n as f64 * 0.5) * dx; 3];
        let mut g = SdfGrid::new(n, n, n, dx, [0.0; 3]);
        g.compute_sphere_sdf(center, radius);
        g
    }
    #[test]
    fn test_marching_cubes_sphere_produces_triangles() {
        let g = sphere_grid(15, 0.1, 0.5);
        let tris = marching_cubes(&g, 0.0);
        assert!(
            !tris.is_empty(),
            "marching cubes on sphere should produce triangles"
        );
    }
    #[test]
    fn test_marching_cubes_all_positive_no_triangles() {
        let mut g = SdfGrid::new(5, 5, 5, 0.1, [0.0; 3]);
        g.values.iter_mut().for_each(|v| *v = 1.0);
        let tris = marching_cubes(&g, 0.0);
        assert!(tris.is_empty(), "no triangles when all SDF > 0");
    }
    #[test]
    fn test_marching_cubes_triangle_count_increases_with_resolution() {
        let g_lo = sphere_grid(8, 0.2, 0.5);
        let g_hi = sphere_grid(20, 0.1, 0.5);
        let n_lo = mesh_triangle_count(&g_lo, 0.0);
        let n_hi = mesh_triangle_count(&g_hi, 0.0);
        assert!(
            n_hi >= n_lo,
            "finer grid should produce at least as many triangles: lo={n_lo}, hi={n_hi}"
        );
    }
    #[test]
    fn test_marching_cubes_small_grid() {
        let mut g = SdfGrid::new(2, 2, 2, 0.5, [0.0; 3]);
        g.set(0, 0, 0, -0.1);
        g.values.iter_mut().skip(1).for_each(|v| *v = 0.5);
        let tris = marching_cubes(&g, 0.0);
        let _ = tris;
    }
    #[test]
    fn test_gaussian_blur_preserves_size() {
        let g = sphere_grid(10, 0.1, 0.4);
        let blurred = sdf_gaussian_blur(&g, 1.0);
        assert_eq!(blurred.nx, g.nx);
        assert_eq!(blurred.ny, g.ny);
        assert_eq!(blurred.nz, g.nz);
    }
    #[test]
    fn test_gaussian_blur_reduces_extremes() {
        let g = sphere_grid(15, 0.1, 0.5);
        let (lo_before, hi_before) = g
            .values
            .iter()
            .fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), &v| {
                (lo.min(v), hi.max(v))
            });
        let blurred = sdf_gaussian_blur(&g, 1.5);
        let (lo_after, hi_after) = blurred
            .values
            .iter()
            .fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), &v| {
                (lo.min(v), hi.max(v))
            });
        assert!(
            lo_after >= lo_before - 1e-6,
            "blur should raise minimum: {lo_before} → {lo_after}"
        );
        assert!(
            hi_after <= hi_before + 1e-6,
            "blur should lower maximum: {hi_before} → {hi_after}"
        );
    }
    #[test]
    fn test_laplacian_sharpen_size() {
        let g = sphere_grid(8, 0.1, 0.3);
        let sharp = sdf_laplacian_sharpen(&g, 0.001);
        assert_eq!(sharp.values.len(), g.values.len());
    }
    #[test]
    fn test_sdf_dilate_expands() {
        let g = sphere_grid(15, 0.1, 0.3);
        let dilated = sdf_dilate(&g, 0.1);
        let count_orig = g.values.iter().filter(|&&v| v < 0.0).count();
        let count_dil = dilated.values.iter().filter(|&&v| v < 0.0).count();
        assert!(count_dil >= count_orig, "dilation should expand interior");
    }
    #[test]
    fn test_sdf_erode_shrinks() {
        let g = sphere_grid(15, 0.1, 0.3);
        let eroded = sdf_erode(&g, 0.05);
        let count_orig = g.values.iter().filter(|&&v| v < 0.0).count();
        let count_er = eroded.values.iter().filter(|&&v| v < 0.0).count();
        assert!(count_er <= count_orig, "erosion should shrink interior");
    }
    #[test]
    fn test_sdf_open_leq_original() {
        let g = sphere_grid(15, 0.1, 0.3);
        let opened = sdf_open(&g, 0.05);
        let count_orig = g.values.iter().filter(|&&v| v < 0.0).count();
        let count_open = opened.values.iter().filter(|&&v| v < 0.0).count();
        assert!(
            count_open <= count_orig + 5,
            "open should not significantly expand"
        );
    }
    #[test]
    fn test_sdf_close_geq_original() {
        let g = sphere_grid(15, 0.1, 0.3);
        let closed = sdf_close(&g, 0.05);
        let count_orig = g.values.iter().filter(|&&v| v < 0.0).count();
        let count_close = closed.values.iter().filter(|&&v| v < 0.0).count();
        assert!(
            count_close >= count_orig - 5,
            "close should not significantly shrink"
        );
    }
    #[test]
    fn test_sdf_offset_surface() {
        let g = sphere_grid(15, 0.1, 0.3);
        let offset = sdf_offset_surface(&g, 0.05);
        for (&orig, &off) in g.values.iter().zip(offset.values.iter()) {
            assert!((off - (orig - 0.05)).abs() < 1e-12);
        }
    }
    #[test]
    fn test_laplacian_smooth_preserves_size() {
        let g = sphere_grid(8, 0.1, 0.3);
        let smoothed = sdf_laplacian_smooth(&g, 3, 0.01);
        assert_eq!(smoothed.values.len(), g.values.len());
    }
    #[test]
    fn test_mean_curvature_smooth_size() {
        let g = sphere_grid(8, 0.1, 0.3);
        let smoothed = sdf_mean_curvature_smooth(&g, 0.001);
        assert_eq!(smoothed.values.len(), g.values.len());
    }
    #[test]
    fn test_bilateral_filter_size() {
        let g = sphere_grid(8, 0.1, 0.3);
        let filtered = sdf_bilateral_filter(&g, 1.5, 0.1);
        assert_eq!(filtered.values.len(), g.values.len());
    }
    #[test]
    fn test_bilateral_filter_preserves_sign() {
        let n = 15usize;
        let dx = 0.1;
        let center = [(n as f64 * 0.5) * dx; 3];
        let mut g = SdfGrid::new(n, n, n, dx, [0.0; 3]);
        g.compute_sphere_sdf(center, 0.5);
        let c = n / 2;
        assert!(g.get(c, c, c) < 0.0, "center should be inside");
        let filtered = sdf_bilateral_filter(&g, 1.0, 0.2);
        assert!(
            filtered.get(c, c, c) < 0.0,
            "center should remain inside after filter"
        );
    }
    #[test]
    fn test_query_distance_field_inside() {
        let n = 21usize;
        let dx = 0.1;
        let center = [(n as f64 * 0.5) * dx; 3];
        let mut g = SdfGrid::new(n, n, n, dx, [0.0; 3]);
        g.compute_sphere_sdf(center, 0.5);
        let q = query_distance_field(&g, center).expect("should return query");
        assert!(q.is_inside, "center should be inside");
        assert!(q.distance < 0.0, "distance at center should be negative");
    }
    #[test]
    fn test_query_distance_field_outside() {
        let n = 21usize;
        let dx = 0.1;
        let center = [(n as f64 * 0.5) * dx; 3];
        let mut g = SdfGrid::new(n, n, n, dx, [0.0; 3]);
        g.compute_sphere_sdf(center, 0.3);
        let far = [center[0] + 0.8, center[1], center[2]];
        if let Some(q) = query_distance_field(&g, far)
            && q.distance.is_finite()
        {
            assert!(!q.is_inside, "far point should be outside");
        }
    }
    #[test]
    fn test_query_batch() {
        let g = sphere_grid(15, 0.1, 0.4);
        let center = [(15_f64 * 0.5) * 0.1; 3];
        let pts = vec![center, [0.0, 0.0, 0.0]];
        let results = query_distance_field_batch(&g, &pts);
        assert_eq!(results.len(), 2);
    }
    #[test]
    fn test_find_zero_crossing() {
        let n = 31usize;
        let dx = 0.05;
        let center = [(n as f64 * 0.5) * dx; 3];
        let radius = 0.4;
        let mut g = SdfGrid::new(n, n, n, dx, [0.0; 3]);
        g.compute_sphere_sdf(center, radius);
        let origin = [0.1, center[1], center[2]];
        let direction = [1.0, 0.0, 0.0];
        let t = find_zero_crossing(&g, origin, direction, 0.0, 1.0, 20);
        assert!(t.is_some(), "should find zero crossing");
        let t_val = t.unwrap();
        let expected_t = center[0] - radius - origin[0];
        assert!(
            (t_val - expected_t).abs() < 0.1,
            "t_val={t_val}, expected≈{expected_t}"
        );
    }
    #[test]
    fn test_projected_area_sphere() {
        let n = 21usize;
        let dx = 0.1;
        let center = [(n as f64 * 0.5) * dx; 3];
        let radius = 0.4;
        let mut g = SdfGrid::new(n, n, n, dx, [0.0; 3]);
        g.compute_sphere_sdf(center, radius);
        let area = projected_area_xy(&g);
        let expected_area = std::f64::consts::PI * radius * radius;
        assert!(
            (area - expected_area).abs() / expected_area < 0.3,
            "projected area {area} vs expected {expected_area}"
        );
    }
}
/// Generate a [`GpuSdfGrid`] by evaluating `shape` at every cell centre.
///
/// * `origin`  — world-space corner of the grid.
/// * `extent`  — total size in each dimension.
/// * `res`     — number of cells in each dimension.
pub fn generate_sdf_grid(
    shape: &SdfShape,
    origin: [f64; 3],
    extent: [f64; 3],
    res: usize,
) -> GpuSdfGrid {
    let cell_size = extent[0] / res as f64;
    let mut grid = GpuSdfGrid::new(res, res, res, origin, cell_size);
    for ix in 0..res {
        for iy in 0..res {
            for iz in 0..res {
                let p = grid.cell_center(ix, iy, iz);
                let idx = grid.index(ix, iy, iz);
                grid.data[idx] = shape.signed_distance(p);
            }
        }
    }
    grid
}
/// Simplified surface extraction from a [`GpuSdfGrid`].
///
/// For each cell of the grid, the function checks whether any of the 12 cube
/// edges cross the iso-surface.  For each crossing edge, it linearly
/// interpolates the crossing point and emits a degenerate triangle
/// (three copies of the crossing point) as a placeholder.  This is *not* a
/// full marching-cubes implementation — it counts crossings and returns one
/// "triangle" per crossing — but it exercises the grid access pattern and
/// confirms that crossings are correctly detected.
///
/// Returns a `Vec<[[f64;3\];3]>` of triangles (three vertices each).
pub fn march_surface(grid: &GpuSdfGrid, iso: f64) -> Vec<[[f64; 3]; 3]> {
    let mut triangles: Vec<[[f64; 3]; 3]> = Vec::new();
    if grid.nx < 2 || grid.ny < 2 || grid.nz < 2 {
        return triangles;
    }
    for ix in 0..grid.nx - 1 {
        for iy in 0..grid.ny - 1 {
            for iz in 0..grid.nz - 1 {
                let corners: [[usize; 3]; 8] = [
                    [ix, iy, iz],
                    [ix + 1, iy, iz],
                    [ix + 1, iy + 1, iz],
                    [ix, iy + 1, iz],
                    [ix, iy, iz + 1],
                    [ix + 1, iy, iz + 1],
                    [ix + 1, iy + 1, iz + 1],
                    [ix, iy + 1, iz + 1],
                ];
                let edges: [[usize; 2]; 12] = [
                    [0, 1],
                    [1, 2],
                    [2, 3],
                    [3, 0],
                    [4, 5],
                    [5, 6],
                    [6, 7],
                    [7, 4],
                    [0, 4],
                    [1, 5],
                    [2, 6],
                    [3, 7],
                ];
                for edge in &edges {
                    let [ia, ib] = *edge;
                    let ca = corners[ia];
                    let cb = corners[ib];
                    let da = grid.get(ca[0], ca[1], ca[2]);
                    let db = grid.get(cb[0], cb[1], cb[2]);
                    if (da < iso) != (db < iso) {
                        let t = if (db - da).abs() > 1e-12 {
                            (iso - da) / (db - da)
                        } else {
                            0.5
                        };
                        let pa = grid.cell_center(ca[0], ca[1], ca[2]);
                        let pb = grid.cell_center(cb[0], cb[1], cb[2]);
                        let pt = [
                            pa[0] + t * (pb[0] - pa[0]),
                            pa[1] + t * (pb[1] - pa[1]),
                            pa[2] + t * (pb[2] - pa[2]),
                        ];
                        triangles.push([pt, pt, pt]);
                    }
                }
            }
        }
    }
    triangles
}
#[cfg(test)]
mod gpu_sdf_tests {

    use crate::sdf_compute::SdfCombine;
    use crate::sdf_compute::SdfShape;
    use crate::sdf_compute::generate_sdf_grid;
    use crate::sdf_compute::march_surface;
    #[test]
    fn test_sphere_sdf_at_center() {
        let r = 1.5;
        let center = [0.0, 0.0, 0.0];
        let shape = SdfShape::Sphere { center, r };
        let d = shape.signed_distance(center);
        assert!(
            (d - (-r)).abs() < 1e-12,
            "SDF at center should be -r, got {d}"
        );
    }
    #[test]
    fn test_sphere_sdf_outside() {
        let r = 1.0;
        let shape = SdfShape::Sphere {
            center: [0.0; 3],
            r,
        };
        let d = shape.signed_distance([3.0, 0.0, 0.0]);
        assert!((d - 2.0).abs() < 1e-12, "SDF outside sphere, got {d}");
    }
    #[test]
    fn test_box_sdf_outside() {
        let shape = SdfShape::Box3 {
            center: [0.0; 3],
            half: [1.0, 1.0, 1.0],
        };
        let d = shape.signed_distance([3.0, 0.0, 0.0]);
        assert!(d > 0.0, "SDF outside box should be positive, got {d}");
    }
    #[test]
    fn test_smooth_union_between_two_spheres() {
        let sa = SdfShape::Sphere {
            center: [-0.5, 0.0, 0.0],
            r: 1.0,
        };
        let sb = SdfShape::Sphere {
            center: [0.5, 0.0, 0.0],
            r: 1.0,
        };
        let combo = SdfCombine::SmoothUnion(sa, sb, 0.5);
        let d = combo.signed_distance([0.0, 0.0, 0.0]);
        assert!(d < 0.0, "smooth union midpoint should be inside, got {d}");
    }
    #[test]
    fn test_hard_union_and_intersection() {
        let sa = SdfShape::Sphere {
            center: [0.0; 3],
            r: 2.0,
        };
        let sb = SdfShape::Sphere {
            center: [0.0; 3],
            r: 1.0,
        };
        let union = SdfCombine::Union(sa.clone(), sb.clone());
        let inter = SdfCombine::Intersection(sa, sb);
        let p = [0.0, 0.0, 0.0];
        assert!((union.signed_distance(p) - (-2.0)).abs() < 1e-12);
        assert!((inter.signed_distance(p) - (-1.0)).abs() < 1e-12);
    }
    #[test]
    fn test_generate_sdf_grid_sphere_center() {
        let r = 0.4;
        let shape = SdfShape::Sphere {
            center: [0.5, 0.5, 0.5],
            r,
        };
        let grid = generate_sdf_grid(&shape, [0.0; 3], [1.0, 1.0, 1.0], 11);
        let mid = 5;
        let d = grid.get(mid, mid, mid);
        assert!(
            (d - (-r)).abs() < 0.1,
            "grid at sphere centre ≈ -r, got {d}"
        );
    }
    #[test]
    fn test_grid_gradient_points_away_from_sphere() {
        let r = 0.3;
        let center = [0.5, 0.5, 0.5];
        let shape = SdfShape::Sphere { center, r };
        let grid = generate_sdf_grid(&shape, [0.0; 3], [1.0, 1.0, 1.0], 21);
        let p = [center[0] + r + 0.05, center[1], center[2]];
        let grad = grid.gradient_at(p);
        assert!(
            grad[0] > 0.0,
            "gradient x should be positive, got {:?}",
            grad
        );
    }
    #[test]
    fn test_march_surface_finds_crossings() {
        let r = 0.3;
        let shape = SdfShape::Sphere {
            center: [0.5, 0.5, 0.5],
            r,
        };
        let grid = generate_sdf_grid(&shape, [0.0; 3], [1.0, 1.0, 1.0], 11);
        let tris = march_surface(&grid, 0.0);
        assert!(
            !tris.is_empty(),
            "marching cubes should find iso-surface crossings for a sphere"
        );
    }
    #[test]
    fn test_capsule_sdf() {
        let shape = SdfShape::Capsule {
            a: [0.0, 0.0, 0.0],
            b: [1.0, 0.0, 0.0],
            r: 0.5,
        };
        let d = shape.signed_distance([0.5, 0.0, 0.0]);
        assert!(d < 0.0, "midpoint inside capsule, got {d}");
        let d_far = shape.signed_distance([5.0, 0.0, 0.0]);
        assert!(d_far > 0.0, "far point outside capsule, got {d_far}");
    }
}
