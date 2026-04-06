//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::needless_range_loop)]
use crate::shape::RayHit;
use oxiphysics_core::math::{Real, Vec3};

use super::types::{HeightField, HeightfieldRayTraversal, HeightfieldRaycast};

/// Compute the area of a triangle from 3 vertices as arrays.
#[allow(dead_code)]
pub(super) fn tri_area(a: &[Real; 3], b: &[Real; 3], c: &[Real; 3]) -> Real {
    let ab = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
    let ac = [c[0] - a[0], c[1] - a[1], c[2] - a[2]];
    let cross = [
        ab[1] * ac[2] - ab[2] * ac[1],
        ab[2] * ac[0] - ab[0] * ac[2],
        ab[0] * ac[1] - ab[1] * ac[0],
    ];
    0.5 * (cross[0] * cross[0] + cross[1] * cross[1] + cross[2] * cross[2]).sqrt()
}
/// Intersect a ray (in XZ plane) with an axis-aligned rectangle \[0, w\] x \[0, h\].
/// Returns (t_enter, t_exit) or None if no intersection.
#[allow(dead_code)]
pub(super) fn ray_aabb_xz(
    ox: Real,
    oz: Real,
    dx: Real,
    dz: Real,
    w: Real,
    h: Real,
    max_toi: Real,
) -> Option<(Real, Real)> {
    let (mut tmin, mut tmax) = (Real::NEG_INFINITY, Real::INFINITY);
    if dx.abs() < 1e-12 {
        if ox < 0.0 || ox > w {
            return None;
        }
    } else {
        let t1 = -ox / dx;
        let t2 = (w - ox) / dx;
        let (ta, tb) = if t1 < t2 { (t1, t2) } else { (t2, t1) };
        tmin = tmin.max(ta);
        tmax = tmax.min(tb);
    }
    if dz.abs() < 1e-12 {
        if oz < 0.0 || oz > h {
            return None;
        }
    } else {
        let t1 = -oz / dz;
        let t2 = (h - oz) / dz;
        let (ta, tb) = if t1 < t2 { (t1, t2) } else { (t2, t1) };
        tmin = tmin.max(ta);
        tmax = tmax.min(tb);
    }
    if tmin > tmax || tmax < 0.0 || tmin > max_toi {
        return None;
    }
    Some((tmin, tmax))
}
/// Möller–Trumbore ray-triangle intersection.
/// Returns a RayHit if the ray hits the triangle within \[0, max_toi\].
pub(super) fn ray_triangle(
    origin: &Vec3,
    direction: &Vec3,
    max_toi: Real,
    v0: &Vec3,
    v1: &Vec3,
    v2: &Vec3,
) -> Option<RayHit> {
    let edge1 = v1 - v0;
    let edge2 = v2 - v0;
    let h = direction.cross(&edge2);
    let det = edge1.dot(&h);
    if det.abs() < 1e-10 {
        return None;
    }
    let inv_det = 1.0 / det;
    let s = origin - v0;
    let u = inv_det * s.dot(&h);
    if !(0.0..=1.0).contains(&u) {
        return None;
    }
    let q = s.cross(&edge1);
    let v = inv_det * direction.dot(&q);
    if v < 0.0 || u + v > 1.0 {
        return None;
    }
    let t = inv_det * edge2.dot(&q);
    if t < 0.0 || t > max_toi {
        return None;
    }
    let point = origin + direction * t;
    let normal = edge1.cross(&edge2).normalize();
    Some(RayHit {
        point,
        normal,
        toi: t,
    })
}
/// DDA ray traversal over the heightfield cells.
///
/// For each cell along the ray (in the XZ plane), tests the two triangles
/// using bilinear interpolation and returns the closest hit.
#[allow(dead_code)]
pub fn heightfield_ray_traverse(
    hf: &HeightField,
    ray_origin: [f64; 3],
    ray_dir: [f64; 3],
    max_t: f64,
) -> Option<HeightfieldRaycast> {
    if hf.rows < 2 || hf.cols < 2 {
        return None;
    }
    use oxiphysics_core::math::Vec3;
    let grid_w = (hf.cols - 1) as f64 * hf.scale_x;
    let grid_h = (hf.rows - 1) as f64 * hf.scale_z;
    let (t_enter, t_exit) = ray_aabb_xz(
        ray_origin[0],
        ray_origin[2],
        ray_dir[0],
        ray_dir[2],
        grid_w,
        grid_h,
        max_t,
    )?;
    let eps = 1e-8;
    let step = hf.scale_x.min(hf.scale_z) * 0.5;
    let mut t = t_enter.max(0.0);
    let mut best: Option<HeightfieldRaycast> = None;
    while t <= t_exit.min(max_t) {
        let px = ray_origin[0] + ray_dir[0] * t;
        let pz = ray_origin[2] + ray_dir[2] * t;
        let col = ((px / hf.scale_x).floor() as isize)
            .max(0)
            .min((hf.cols - 2) as isize) as usize;
        let row = ((pz / hf.scale_z).floor() as isize)
            .max(0)
            .min((hf.rows - 2) as isize) as usize;
        let v00 = Vec3::new(
            col as f64 * hf.scale_x,
            hf.height_at(row, col),
            row as f64 * hf.scale_z,
        );
        let v10 = Vec3::new(
            (col + 1) as f64 * hf.scale_x,
            hf.height_at(row, col + 1),
            row as f64 * hf.scale_z,
        );
        let v01 = Vec3::new(
            col as f64 * hf.scale_x,
            hf.height_at(row + 1, col),
            (row + 1) as f64 * hf.scale_z,
        );
        let v11 = Vec3::new(
            (col + 1) as f64 * hf.scale_x,
            hf.height_at(row + 1, col + 1),
            (row + 1) as f64 * hf.scale_z,
        );
        let o = Vec3::new(ray_origin[0], ray_origin[1], ray_origin[2]);
        let d = Vec3::new(ray_dir[0], ray_dir[1], ray_dir[2]);
        for (a, b, c) in [(&v00, &v10, &v11), (&v00, &v11, &v01)] {
            if let Some(hit) = ray_triangle(&o, &d, max_t, a, b, c)
                && best
                    .as_ref()
                    .is_none_or(|bh: &HeightfieldRaycast| hit.toi < bh.t)
            {
                let dot = hit.normal.dot(&d);
                let n = if dot > 0.0 {
                    [-hit.normal.x, -hit.normal.y, -hit.normal.z]
                } else {
                    [hit.normal.x, hit.normal.y, hit.normal.z]
                };
                best = Some(HeightfieldRaycast {
                    cell_ix: col,
                    cell_iz: row,
                    t: hit.toi,
                    normal: n,
                });
            }
        }
        if best.is_some() {
            return best;
        }
        t += step + eps;
    }
    best
}
/// Laplacian smoothing with a blending factor alpha.
///
/// Each interior vertex is blended: `h_new = (1-alpha)*h + alpha*avg_neighbors`.
#[allow(dead_code)]
pub fn heightfield_smooth(hf: &mut HeightField, iterations: usize, alpha: f64) {
    for _ in 0..iterations {
        let mut new_h = hf.heights.clone();
        for row in 1..(hf.rows - 1) {
            for col in 1..(hf.cols - 1) {
                let avg = (hf.height_at(row - 1, col)
                    + hf.height_at(row + 1, col)
                    + hf.height_at(row, col - 1)
                    + hf.height_at(row, col + 1))
                    / 4.0;
                let cur = hf.height_at(row, col);
                new_h[row * hf.cols + col] = (1.0 - alpha) * cur + alpha * avg;
            }
        }
        hf.heights = new_h;
    }
}
/// Simple thermal erosion: for each cell, compute height excess over the
/// average of lower neighbors and redistribute `rate` fraction to them.
#[allow(dead_code)]
pub fn heightfield_erode(hf: &mut HeightField, rate: f64, iterations: usize) {
    for _ in 0..iterations {
        let mut delta = vec![0.0f64; hf.rows * hf.cols];
        for row in 0..hf.rows {
            for col in 0..hf.cols {
                let h = hf.height_at(row, col);
                let neighbors: [(isize, isize); 4] = [
                    (row as isize - 1, col as isize),
                    (row as isize + 1, col as isize),
                    (row as isize, col as isize - 1),
                    (row as isize, col as isize + 1),
                ];
                let mut lower_count = 0usize;
                let mut total_diff = 0.0f64;
                for (nr, nc) in neighbors {
                    if nr >= 0 && nr < hf.rows as isize && nc >= 0 && nc < hf.cols as isize {
                        let nh = hf.height_at(nr as usize, nc as usize);
                        if nh < h {
                            total_diff += h - nh;
                            lower_count += 1;
                        }
                    }
                }
                if lower_count > 0 && total_diff > 0.0 {
                    let transfer = rate * total_diff / lower_count as f64;
                    for (nr, nc) in neighbors {
                        if nr >= 0 && nr < hf.rows as isize && nc >= 0 && nc < hf.cols as isize {
                            let nh = hf.height_at(nr as usize, nc as usize);
                            if nh < h {
                                let frac = (h - nh) / total_diff;
                                delta[nr as usize * hf.cols + nc as usize] += transfer * frac;
                                delta[row * hf.cols + col] -= transfer * frac;
                            }
                        }
                    }
                }
            }
        }
        for i in 0..hf.heights.len() {
            hf.heights[i] += delta[i];
        }
    }
}
/// Compute per-vertex normals for all grid vertices via finite differences.
///
/// Returns a flat `Vec<[f64;3]>` in row-major order (same as `compute_all_normals`).
#[allow(dead_code)]
pub fn heightfield_compute_normals(hf: &HeightField) -> Vec<[f64; 3]> {
    hf.compute_all_normals()
}
/// Enumerate all triangles in the heightfield as a flat list of vertex triples.
#[allow(dead_code)]
pub fn heightfield_to_triangle_list(hf: &HeightField) -> Vec<[[f64; 3]; 3]> {
    let mut tris =
        Vec::with_capacity((hf.rows.saturating_sub(1)) * (hf.cols.saturating_sub(1)) * 2);
    for row in 0..(hf.rows.saturating_sub(1)) {
        for col in 0..(hf.cols.saturating_sub(1)) {
            let v00 = [
                col as f64 * hf.scale_x,
                hf.height_at(row, col),
                row as f64 * hf.scale_z,
            ];
            let v10 = [
                (col + 1) as f64 * hf.scale_x,
                hf.height_at(row, col + 1),
                row as f64 * hf.scale_z,
            ];
            let v01 = [
                col as f64 * hf.scale_x,
                hf.height_at(row + 1, col),
                (row + 1) as f64 * hf.scale_z,
            ];
            let v11 = [
                (col + 1) as f64 * hf.scale_x,
                hf.height_at(row + 1, col + 1),
                (row + 1) as f64 * hf.scale_z,
            ];
            tris.push([v00, v10, v11]);
            tris.push([v00, v11, v01]);
        }
    }
    tris
}
/// Return `(min_height, max_height)` over all grid vertices.
#[allow(dead_code)]
pub fn heightfield_min_max(hf: &HeightField) -> (f64, f64) {
    let mut min = f64::INFINITY;
    let mut max = f64::NEG_INFINITY;
    for &h in &hf.heights {
        if h < min {
            min = h;
        }
        if h > max {
            max = h;
        }
    }
    (min, max)
}
/// Compute the axis-aligned bounding box of the heightfield.
///
/// Returns `(min_corner, max_corner)` as `[f64; 3]` arrays.
#[allow(dead_code)]
pub fn heightfield_aabb(hf: &HeightField) -> ([f64; 3], [f64; 3]) {
    let (min_h, max_h) = heightfield_min_max(hf);
    let min = [0.0, min_h, 0.0];
    let max = [
        hf.scale_x * (hf.cols.saturating_sub(1)) as f64,
        max_h,
        hf.scale_z * (hf.rows.saturating_sub(1)) as f64,
    ];
    (min, max)
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::HeightField;

    #[test]
    fn test_heightfield_bounding_box() {
        let heights = vec![0.0, 1.0, 2.0, 3.0, 1.0, 0.5, 1.5, 2.5, 0.0];
        let hf = HeightField::new(heights, 3, 3, 1.0, 1.0);
        let (bb_min, bb_max) = hf.aabb();
        assert!((bb_min[1] - 0.0).abs() < 1e-10);
        assert!((bb_max[1] - 3.0).abs() < 1e-10);
    }
    #[test]
    fn test_heightfield_height_at() {
        let heights = vec![1.0, 2.0, 3.0, 4.0];
        let hf = HeightField::new(heights, 2, 2, 1.0, 1.0);
        assert!((hf.height_at(0, 1) - 2.0).abs() < 1e-10);
        assert!((hf.height_at(1, 0) - 3.0).abs() < 1e-10);
    }
    #[test]
    fn test_heightfield_ray_cast() {
        let heights = vec![0.0, 0.0, 0.0, 0.0];
        let hf = HeightField::new(heights, 2, 2, 1.0, 1.0);
        let origin = Vec3::new(0.5, 5.0, 0.5);
        let dir = Vec3::new(0.0, -1.0, 0.0);
        let hit = hf.ray_cast(&origin, &dir, 100.0);
        assert!(hit.is_some(), "expected a hit on flat terrain");
        let hit = hit.unwrap();
        assert!(
            (hit.toi - 5.0).abs() < 1e-6,
            "expected toi=5, got {}",
            hit.toi
        );
        assert!(
            (hit.point.y).abs() < 1e-6,
            "expected hit at y=0, got y={}",
            hit.point.y
        );
    }
    #[test]
    fn test_heightfield_volume_flat() {
        let heights = vec![2.0, 2.0, 2.0, 2.0];
        let hf = HeightField::new(heights, 2, 2, 1.0, 1.0);
        assert!(
            (hf.volume() - 2.0).abs() < 1e-10,
            "expected volume=2.0, got {}",
            hf.volume()
        );
    }
    #[test]
    fn test_heightfield_volume_positive() {
        let heights = vec![0.0, 1.0, 2.0, 3.0, 1.0, 0.5, 1.5, 2.5, 0.0];
        let hf = HeightField::new(heights, 3, 3, 1.0, 1.0);
        assert!(
            hf.volume() > 0.0,
            "expected positive volume for non-flat terrain"
        );
    }
    #[test]
    fn test_flat_field_const_height() {
        let hf = HeightField::from_fn(4, 4, 1.0, |_, _| 5.0);
        assert!((hf.min_height() - 5.0).abs() < 1e-10);
        assert!((hf.max_height() - 5.0).abs() < 1e-10);
        for row in 0..4 {
            for col in 0..4 {
                assert!((hf.height_at(row, col) - 5.0).abs() < 1e-10);
            }
        }
    }
    #[test]
    fn test_interpolation_corners() {
        let hf = HeightField::new(vec![0.0, 1.0, 2.0, 3.0], 2, 2, 1.0, 1.0);
        assert!((hf.height_at_uv(0.0, 0.0) - 0.0).abs() < 1e-10);
        assert!((hf.height_at_uv(1.0, 0.0) - 1.0).abs() < 1e-10);
        assert!((hf.height_at_uv(0.0, 1.0) - 2.0).abs() < 1e-10);
        assert!((hf.height_at_uv(1.0, 1.0) - 3.0).abs() < 1e-10);
    }
    #[test]
    fn test_interpolation_center() {
        let hf = HeightField::new(vec![0.0, 1.0, 2.0, 3.0], 2, 2, 1.0, 1.0);
        let center = hf.height_at_uv(0.5, 0.5);
        assert!((center - 1.5).abs() < 1e-10, "expected 1.5, got {}", center);
    }
    #[test]
    fn test_tessellation_tri_count() {
        let hf = HeightField::from_fn(3, 4, 1.0, |_, _| 0.0);
        let (verts, tris) = hf.to_triangle_mesh();
        assert_eq!(verts.len(), 12);
        assert_eq!(tris.len(), 12);
    }
    #[test]
    fn test_tessellation_vertex_count() {
        let hf = HeightField::from_fn(5, 5, 0.5, |c, r| (c + r) as Real);
        let (verts, tris) = hf.to_triangle_mesh();
        assert_eq!(verts.len(), 25);
        assert_eq!(tris.len(), 32);
    }
    #[test]
    fn test_normal_flat_field_points_up() {
        let hf = HeightField::from_fn(4, 4, 1.0, |_, _| 0.0);
        for row in 0..4 {
            for col in 0..4 {
                let n = hf.normal_at_grid(col, row);
                assert!((n[0]).abs() < 1e-10, "nx should be 0");
                assert!((n[1] - 1.0).abs() < 1e-10, "ny should be 1");
                assert!((n[2]).abs() < 1e-10, "nz should be 0");
            }
        }
    }
    #[test]
    fn test_smoothing_reduces_roughness() {
        let mut hf =
            HeightField::from_fn(5, 5, 1.0, |c, r| if c == 2 && r == 2 { 10.0 } else { 0.0 });
        let range_before = hf.max_height() - hf.min_height();
        hf.smooth(3);
        let range_after = hf.max_height() - hf.min_height();
        assert!(
            range_after < range_before,
            "smoothing should reduce height range: {} >= {}",
            range_after,
            range_before
        );
    }
    #[test]
    fn test_surface_area_flat() {
        let hf = HeightField::from_fn(3, 3, 1.0, |_, _| 0.0);
        let area = hf.surface_area();
        assert!(
            (area - 4.0).abs() < 1e-10,
            "expected area=4.0, got {}",
            area
        );
    }
    #[test]
    fn test_ray_cast_grid_flat() {
        let hf = HeightField::from_fn(4, 4, 1.0, |_, _| 0.0);
        let origin = [1.5, 10.0, 1.5];
        let dir = [0.0, -1.0, 0.0];
        let hit = hf.ray_cast_grid(origin, dir, 100.0);
        assert!(hit.is_some(), "should hit flat terrain");
        let (toi, normal) = hit.unwrap();
        assert!((toi - 10.0).abs() < 1e-4, "toi={} expected 10", toi);
        assert!((normal[1] - 1.0).abs() < 1e-4, "normal should point up");
    }
    #[test]
    fn test_ray_cast_grid_miss() {
        let hf = HeightField::from_fn(4, 4, 1.0, |_, _| 0.0);
        let origin = [0.0, 5.0, 1.5];
        let dir = [1.0, 0.0, 0.0];
        let hit = hf.ray_cast_grid(origin, dir, 100.0);
        assert!(hit.is_none(), "horizontal ray should miss flat terrain");
    }
    #[test]
    fn test_from_fn_creates_correct_grid() {
        let hf = HeightField::from_fn(3, 2, 2.0, |c, r| (c * 10 + r) as Real);
        assert_eq!(hf.cols, 3);
        assert_eq!(hf.rows, 2);
        assert!((hf.scale_x - 2.0).abs() < 1e-10);
        assert!((hf.height_at(0, 0) - 0.0).abs() < 1e-10);
        assert!((hf.height_at(0, 2) - 20.0).abs() < 1e-10);
        assert!((hf.height_at(1, 0) - 1.0).abs() < 1e-10);
        assert!((hf.height_at(1, 1) - 11.0).abs() < 1e-10);
    }
    #[test]
    fn test_min_max_height() {
        let hf = HeightField::new(vec![-3.0, 1.0, 5.0, 2.0], 2, 2, 1.0, 1.0);
        assert!((hf.min_height() - (-3.0)).abs() < 1e-10);
        assert!((hf.max_height() - 5.0).abs() < 1e-10);
    }
    #[test]
    fn test_lod_downsample_reduces_size() {
        let hf = HeightField::from_fn(8, 8, 1.0, |c, r| (c + r) as Real);
        let lod = hf.lod_downsample(2);
        assert_eq!(lod.cols, 4);
        assert_eq!(lod.rows, 4);
        assert!((lod.scale_x - 2.0).abs() < 1e-10);
    }
    #[test]
    fn test_lod_factor_one_unchanged() {
        let hf = HeightField::from_fn(4, 4, 1.0, |c, r| (c * r) as Real);
        let lod = hf.lod_downsample(1);
        assert_eq!(lod.rows, hf.rows);
        assert_eq!(lod.cols, hf.cols);
    }
    #[test]
    fn test_lod_pyramid_has_multiple_levels() {
        let hf = HeightField::from_fn(16, 16, 1.0, |_, _| 0.0);
        let pyramid = hf.lod_pyramid(4);
        assert!(pyramid.len() >= 2, "should produce at least 2 levels");
        assert!(pyramid[1].cols <= pyramid[0].cols);
    }
    #[test]
    fn test_normal_at_world_flat_is_up() {
        let hf = HeightField::from_fn(4, 4, 1.0, |_, _| 0.0);
        let n = hf.normal_at_world(1.5, 1.5);
        assert!(
            (n[1] - 1.0).abs() < 1e-6,
            "flat terrain normal should be [0,1,0]"
        );
    }
    #[test]
    fn test_serialize_deserialize_roundtrip() {
        let hf = HeightField::from_fn(4, 4, 2.0, |c, r| (c + r) as Real);
        let data = hf.serialize();
        let hf2 = HeightField::deserialize(&data).expect("deserialize failed");
        assert_eq!(hf2.rows, hf.rows);
        assert_eq!(hf2.cols, hf.cols);
        assert!((hf2.scale_x - hf.scale_x).abs() < 1e-10);
        for (a, b) in hf.heights.iter().zip(hf2.heights.iter()) {
            assert!((a - b).abs() < 1e-10);
        }
    }
    #[test]
    fn test_deserialize_bad_data_returns_none() {
        let data = vec![2.0, 2.0, 1.0];
        assert!(HeightField::deserialize(&data).is_none());
    }
    #[test]
    fn test_compute_all_normals_count() {
        let hf = HeightField::from_fn(3, 4, 1.0, |_, _| 0.0);
        let normals = hf.compute_all_normals();
        assert_eq!(normals.len(), 3 * 4);
    }
    #[test]
    fn test_compute_all_normals_unit_length() {
        let hf = HeightField::from_fn(4, 4, 1.0, |c, r| (c as Real).sin() + (r as Real).cos());
        for n in hf.compute_all_normals() {
            let len = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
            assert!((len - 1.0).abs() < 1e-9, "normal not unit length: {len}");
        }
    }
    #[test]
    fn test_height_at_world_flat() {
        let hf = HeightField::from_fn(4, 4, 1.0, |_, _| 3.0);
        let h = hf.height_at_world(1.5, 2.0);
        assert!(
            (h - 3.0).abs() < 1e-9,
            "flat field should return 3.0 everywhere"
        );
    }
    #[test]
    fn test_height_at_world_clamps_to_extents() {
        let hf = HeightField::from_fn(3, 3, 1.0, |c, _r| c as Real);
        let h = hf.height_at_world(100.0, 1.0);
        assert!(h.is_finite());
    }
    #[test]
    fn test_ray_cast_bounded_flat_hit() {
        let hf = HeightField::from_fn(8, 8, 1.0, |_, _| 0.0);
        let origin = [3.5, 5.0, 3.5];
        let dir = [0.0, -1.0, 0.0];
        let hit = hf.ray_cast_bounded(origin, dir, 20.0, 1000);
        assert!(hit.is_some(), "should hit flat terrain");
        let (toi, _) = hit.unwrap();
        assert!((toi - 5.0).abs() < 0.5, "toi should be near 5, got {toi}");
    }
    #[test]
    fn test_ray_cast_bounded_miss_max_steps() {
        let hf = HeightField::from_fn(4, 4, 1.0, |_, _| 0.0);
        let origin = [0.5, 10.0, 0.5];
        let dir = [1.0, 0.0, 0.0];
        let hit = hf.ray_cast_bounded(origin, dir, 1000.0, 1);
        let _ = hit;
    }
    #[test]
    fn test_heightfield_min_max_uniform() {
        let hf = HeightField::from_fn(4, 4, 1.0, |_, _| 7.0);
        let (mn, mx) = heightfield_min_max(&hf);
        assert!((mn - 7.0).abs() < 1e-10, "min={mn}");
        assert!((mx - 7.0).abs() < 1e-10, "max={mx}");
    }
    #[test]
    fn test_heightfield_min_max_varying() {
        let hf = HeightField::new(vec![-1.0, 3.0, 5.0, 2.0], 2, 2, 1.0, 1.0);
        let (mn, mx) = heightfield_min_max(&hf);
        assert!((mn - (-1.0)).abs() < 1e-10);
        assert!((mx - 5.0).abs() < 1e-10);
    }
    #[test]
    fn test_heightfield_aabb_contains_all_vertices() {
        let hf = HeightField::new(vec![-2.0, 4.0, 1.0, 3.0], 2, 2, 2.0, 3.0);
        let (bb_min, bb_max) = heightfield_aabb(&hf);
        for row in 0..hf.rows {
            for col in 0..hf.cols {
                let x = col as f64 * hf.scale_x;
                let y = hf.height_at(row, col);
                let z = row as f64 * hf.scale_z;
                assert!(
                    x >= bb_min[0] && x <= bb_max[0],
                    "x={x} not in [{},{}]",
                    bb_min[0],
                    bb_max[0]
                );
                assert!(
                    y >= bb_min[1] && y <= bb_max[1],
                    "y={y} not in [{},{}]",
                    bb_min[1],
                    bb_max[1]
                );
                assert!(
                    z >= bb_min[2] && z <= bb_max[2],
                    "z={z} not in [{},{}]",
                    bb_min[2],
                    bb_max[2]
                );
            }
        }
    }
    #[test]
    fn test_heightfield_smooth_reduces_variance() {
        let mut hf =
            HeightField::from_fn(5, 5, 1.0, |c, r| if c == 2 && r == 2 { 100.0 } else { 0.0 });
        let var_before: f64 = {
            let mean = hf.heights.iter().sum::<f64>() / hf.heights.len() as f64;
            hf.heights.iter().map(|&h| (h - mean).powi(2)).sum::<f64>() / hf.heights.len() as f64
        };
        heightfield_smooth(&mut hf, 3, 1.0);
        let var_after: f64 = {
            let mean = hf.heights.iter().sum::<f64>() / hf.heights.len() as f64;
            hf.heights.iter().map(|&h| (h - mean).powi(2)).sum::<f64>() / hf.heights.len() as f64
        };
        assert!(
            var_after < var_before,
            "smoothing should reduce variance: {var_after} >= {var_before}"
        );
    }
    #[test]
    fn test_heightfield_compute_normals_nonzero() {
        let hf = HeightField::from_fn(4, 4, 1.0, |c, r| (c + r) as f64 * 0.5);
        let normals = heightfield_compute_normals(&hf);
        assert_eq!(normals.len(), 4 * 4);
        for n in &normals {
            let len = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
            assert!(len > 0.5, "normal should be non-zero, len={len}");
        }
    }
    #[test]
    fn test_heightfield_to_triangle_list_count() {
        let hf = HeightField::from_fn(3, 4, 1.0, |_, _| 0.0);
        let tris = heightfield_to_triangle_list(&hf);
        assert_eq!(tris.len(), 12, "expected 12 triangles, got {}", tris.len());
    }
    #[test]
    fn test_heightfield_ray_traverse_hits_flat() {
        let hf = HeightField::from_fn(4, 4, 1.0, |_, _| 0.0);
        let origin = [1.5, 5.0, 1.5];
        let dir = [0.0, -1.0, 0.0];
        let hit = heightfield_ray_traverse(&hf, origin, dir, 100.0);
        assert!(hit.is_some(), "should hit flat terrain");
        let h = hit.unwrap();
        assert!((h.t - 5.0).abs() < 0.1, "t≈5, got {}", h.t);
    }
    #[test]
    fn test_heightfield_erode_lowers_peaks() {
        let mut hf =
            HeightField::from_fn(5, 5, 1.0, |c, r| if c == 2 && r == 2 { 10.0 } else { 0.0 });
        let peak_before = hf.height_at(2, 2);
        heightfield_erode(&mut hf, 0.5, 5);
        let peak_after = hf.height_at(2, 2);
        assert!(peak_after < peak_before, "erosion should lower the peak");
    }
}
/// Height-field ray cast using DDA grid traversal.
///
/// Traverses grid cells in order along the ray projection, testing the two
/// triangles of each cell against the ray.  Returns the world-space hit
/// position `[x, y, z]` of the first (nearest) intersection, or `None` on
/// a miss.
///
/// For purely vertical rays (no X or Z component), a single-cell look-up is
/// performed at the ray's (x, z) position.
#[allow(dead_code)]
pub fn heightfield_ray_intersect(
    hf: &HeightField,
    ray_origin: [f64; 3],
    ray_dir: [f64; 3],
) -> Option<[f64; 3]> {
    let max_t = 1e15_f64;
    let o = Vec3::new(ray_origin[0], ray_origin[1], ray_origin[2]);
    let d = Vec3::new(ray_dir[0], ray_dir[1], ray_dir[2]);
    let horizontal_len = (ray_dir[0] * ray_dir[0] + ray_dir[2] * ray_dir[2]).sqrt();
    if horizontal_len < 1e-12 {
        let px = ray_origin[0];
        let pz = ray_origin[2];
        let grid_w = (hf.cols - 1) as f64 * hf.scale_x;
        let grid_h = (hf.rows - 1) as f64 * hf.scale_z;
        if px < 0.0 || px > grid_w || pz < 0.0 || pz > grid_h {
            return None;
        }
        let col = ((px / hf.scale_x).floor() as usize).min(hf.cols - 2);
        let row = ((pz / hf.scale_z).floor() as usize).min(hf.rows - 2);
        let v00 = Vec3::new(
            col as f64 * hf.scale_x,
            hf.height_at(row, col),
            row as f64 * hf.scale_z,
        );
        let v10 = Vec3::new(
            (col + 1) as f64 * hf.scale_x,
            hf.height_at(row, col + 1),
            row as f64 * hf.scale_z,
        );
        let v01 = Vec3::new(
            col as f64 * hf.scale_x,
            hf.height_at(row + 1, col),
            (row + 1) as f64 * hf.scale_z,
        );
        let v11 = Vec3::new(
            (col + 1) as f64 * hf.scale_x,
            hf.height_at(row + 1, col + 1),
            (row + 1) as f64 * hf.scale_z,
        );
        for (a, b, c) in [(&v00, &v10, &v11), (&v00, &v11, &v01)] {
            if let Some(hit) = ray_triangle(&o, &d, max_t, a, b, c) {
                return Some([hit.point.x, hit.point.y, hit.point.z]);
            }
        }
        return None;
    }
    let mut traversal = HeightfieldRayTraversal::new(hf, ray_origin, ray_dir, max_t)?;
    while let Some((col, row, _)) = traversal.next_cell() {
        if col >= hf.cols - 1 || row >= hf.rows - 1 {
            continue;
        }
        let v00 = Vec3::new(
            col as f64 * hf.scale_x,
            hf.height_at(row, col),
            row as f64 * hf.scale_z,
        );
        let v10 = Vec3::new(
            (col + 1) as f64 * hf.scale_x,
            hf.height_at(row, col + 1),
            row as f64 * hf.scale_z,
        );
        let v01 = Vec3::new(
            col as f64 * hf.scale_x,
            hf.height_at(row + 1, col),
            (row + 1) as f64 * hf.scale_z,
        );
        let v11 = Vec3::new(
            (col + 1) as f64 * hf.scale_x,
            hf.height_at(row + 1, col + 1),
            (row + 1) as f64 * hf.scale_z,
        );
        for (a, b, c) in [(&v00, &v10, &v11), (&v00, &v11, &v01)] {
            if let Some(hit) = ray_triangle(&o, &d, max_t, a, b, c) {
                return Some([hit.point.x, hit.point.y, hit.point.z]);
            }
        }
    }
    None
}
/// Generate a triangle mesh from a height grid.
///
/// Returns a flat list of triangle indices `[i0, i1, i2]` into the vertex
/// array produced by `to_triangle_mesh`.  Each quad cell is split into two
/// triangles.  The total count is `(rows-1) * (cols-1) * 2`.
#[allow(dead_code)]
pub fn heightfield_tessellate(hf: &HeightField) -> Vec<[usize; 3]> {
    if hf.rows < 2 || hf.cols < 2 {
        return vec![];
    }
    let mut tris = Vec::with_capacity((hf.rows - 1) * (hf.cols - 1) * 2);
    for row in 0..(hf.rows - 1) {
        for col in 0..(hf.cols - 1) {
            let i00 = row * hf.cols + col;
            let i10 = row * hf.cols + (col + 1);
            let i01 = (row + 1) * hf.cols + col;
            let i11 = (row + 1) * hf.cols + (col + 1);
            tris.push([i00, i10, i11]);
            tris.push([i00, i11, i01]);
        }
    }
    tris
}
/// Apply Gaussian blur to the height field, returning a new smoothed `HeightField`.
///
/// `sigma` controls the width of the kernel (in grid cells).  A 1-D separable
/// Gaussian kernel of radius `ceil(3*sigma)` is applied first in the column
/// direction then in the row direction.  Boundary cells are handled by
/// clamping the kernel window to the grid extent.
#[allow(dead_code)]
pub fn heightfield_smooth_gaussian(hf: &HeightField, sigma: f64) -> HeightField {
    if sigma < 1e-10 {
        return hf.clone();
    }
    let radius = (3.0 * sigma).ceil() as isize;
    let kernel_len = (2 * radius + 1) as usize;
    let mut kernel = Vec::with_capacity(kernel_len);
    for k in -radius..=radius {
        let w = (-(k as f64 * k as f64) / (2.0 * sigma * sigma)).exp();
        kernel.push(w);
    }
    let kernel_sum: f64 = kernel.iter().sum();
    let kernel: Vec<f64> = kernel.iter().map(|&w| w / kernel_sum).collect();
    let mut tmp = vec![0.0f64; hf.rows * hf.cols];
    for row in 0..hf.rows {
        for col in 0..hf.cols {
            let mut acc = 0.0f64;
            let mut wsum = 0.0f64;
            for (ki, k) in (-radius..=radius).enumerate() {
                let c = (col as isize + k).clamp(0, hf.cols as isize - 1) as usize;
                acc += kernel[ki] * hf.height_at(row, c);
                wsum += kernel[ki];
            }
            tmp[row * hf.cols + col] = if wsum > 1e-15 {
                acc / wsum
            } else {
                hf.height_at(row, col)
            };
        }
    }
    let mut out = vec![0.0f64; hf.rows * hf.cols];
    for row in 0..hf.rows {
        for col in 0..hf.cols {
            let mut acc = 0.0f64;
            let mut wsum = 0.0f64;
            for (ki, k) in (-radius..=radius).enumerate() {
                let r = (row as isize + k).clamp(0, hf.rows as isize - 1) as usize;
                acc += kernel[ki] * tmp[r * hf.cols + col];
                wsum += kernel[ki];
            }
            out[row * hf.cols + col] = if wsum > 1e-15 {
                acc / wsum
            } else {
                tmp[row * hf.cols + col]
            };
        }
    }
    HeightField::new(out, hf.rows, hf.cols, hf.scale_x, hf.scale_z)
}
/// Compute the surface normal at grid position `(i, j)` using central differences.
///
/// `i` is the column index (X direction) and `j` is the row index (Z direction).
/// Returns a unit normal `[nx, ny, nz]`.
#[allow(dead_code)]
pub fn heightfield_normal_at(hf: &HeightField, i: usize, j: usize) -> [f64; 3] {
    hf.normal_at_grid(i, j)
}
/// Bilinear interpolation of height at arbitrary world-space position `(x, z)`.
///
/// Clamps to the grid extents.
#[allow(dead_code)]
pub fn heightfield_height_at_xy(hf: &HeightField, x: f64, z: f64) -> f64 {
    hf.height_at_world(x, z)
}
#[cfg(test)]
mod tests_extended {

    use crate::HeightField;
    use crate::heightfield::HeightfieldRayTraversal;
    use crate::heightfield::heightfield_height_at_xy;
    use crate::heightfield::heightfield_normal_at;
    use crate::heightfield::heightfield_ray_intersect;
    use crate::heightfield::heightfield_smooth_gaussian;
    use crate::heightfield::heightfield_tessellate;
    #[test]
    fn test_ray_traversal_visits_cells() {
        let hf = HeightField::from_fn(5, 5, 1.0, |_, _| 0.0);
        let origin = [0.1, 10.0, 0.1];
        let dir = [1.0, 0.0, 0.0];
        let traversal = HeightfieldRayTraversal::new(&hf, origin, dir, 100.0);
        assert!(traversal.is_some(), "should create traversal");
        let mut t = traversal.unwrap();
        let mut count = 0;
        while t.next_cell().is_some() {
            count += 1;
        }
        assert!(count >= 1, "should visit at least 1 cell");
    }
    #[test]
    fn test_ray_traversal_outside_grid_returns_none() {
        let hf = HeightField::from_fn(4, 4, 1.0, |_, _| 0.0);
        let origin = [100.0, 10.0, 100.0];
        let dir = [1.0, 0.0, 0.0];
        let traversal = HeightfieldRayTraversal::new(&hf, origin, dir, 50.0);
        assert!(traversal.is_none(), "ray outside grid should return None");
    }
    #[test]
    fn test_ray_intersect_flat_plane() {
        let hf = HeightField::from_fn(4, 4, 1.0, |_, _| 0.0);
        let origin = [1.5, 5.0, 1.5];
        let dir = [0.0, -1.0, 0.0];
        let hit = heightfield_ray_intersect(&hf, origin, dir);
        assert!(hit.is_some(), "should intersect flat terrain");
        let p = hit.unwrap();
        assert!(p[1].abs() < 0.01, "hit y should be ~0, got {}", p[1]);
    }
    #[test]
    fn test_ray_intersect_misses_upward_ray() {
        let hf = HeightField::from_fn(4, 4, 1.0, |_, _| 0.0);
        let origin = [1.5, 5.0, 1.5];
        let dir = [0.0, 1.0, 0.0];
        let hit = heightfield_ray_intersect(&hf, origin, dir);
        assert!(hit.is_none(), "upward ray should miss terrain");
    }
    #[test]
    fn test_tessellate_triangle_count() {
        let hf = HeightField::from_fn(3, 4, 1.0, |_, _| 0.0);
        let tris = heightfield_tessellate(&hf);
        assert_eq!(tris.len(), 12, "expected 12 triangles, got {}", tris.len());
    }
    #[test]
    fn test_tessellate_5x5_grid() {
        let hf = HeightField::from_fn(5, 5, 1.0, |_, _| 0.0);
        let tris = heightfield_tessellate(&hf);
        assert_eq!(tris.len(), 32, "expected 32 triangles, got {}", tris.len());
    }
    #[test]
    fn test_tessellate_single_cell_2x2() {
        let hf = HeightField::from_fn(2, 2, 1.0, |_, _| 0.0);
        let tris = heightfield_tessellate(&hf);
        assert_eq!(tris.len(), 2);
    }
    #[test]
    fn test_tessellate_indices_in_range() {
        let hf = HeightField::from_fn(4, 4, 1.0, |c, r| (c + r) as f64);
        let tris = heightfield_tessellate(&hf);
        let max_idx = hf.rows * hf.cols;
        for tri in &tris {
            for &idx in tri.iter() {
                assert!(idx < max_idx, "index {idx} out of range {max_idx}");
            }
        }
    }
    #[test]
    fn test_gaussian_smooth_flat_stays_flat() {
        let hf = HeightField::from_fn(6, 6, 1.0, |_, _| 3.0);
        let smoothed = heightfield_smooth_gaussian(&hf, 1.0);
        for &h in &smoothed.heights {
            assert!(
                (h - 3.0).abs() < 1e-6,
                "flat field should stay flat after Gaussian smooth, got {h}"
            );
        }
    }
    #[test]
    fn test_gaussian_smooth_reduces_spike() {
        let mut hf = HeightField::from_fn(7, 7, 1.0, |_, _| 0.0);
        hf.heights[3 * 7 + 3] = 100.0;
        let smoothed = heightfield_smooth_gaussian(&hf, 1.5);
        assert!(
            smoothed.height_at(3, 3) < 100.0,
            "spike should be reduced by Gaussian smoothing"
        );
        assert!(
            smoothed.height_at(3, 4) > 0.0,
            "neighbors should gain height after smoothing"
        );
    }
    #[test]
    fn test_gaussian_smooth_zero_sigma_is_identity() {
        let hf = HeightField::from_fn(4, 4, 1.0, |c, r| (c * r) as f64);
        let smoothed = heightfield_smooth_gaussian(&hf, 0.0);
        for (a, b) in hf.heights.iter().zip(smoothed.heights.iter()) {
            assert!(
                (a - b).abs() < 1e-10,
                "zero sigma should return identical heights"
            );
        }
    }
    #[test]
    fn test_normal_at_flat_plane_is_up() {
        let hf = HeightField::from_fn(5, 5, 1.0, |_, _| 0.0);
        for j in 0..hf.rows {
            for i in 0..hf.cols {
                let n = heightfield_normal_at(&hf, i, j);
                assert!(
                    n[0].abs() < 1e-10,
                    "nx should be 0 at ({i},{j}), got {}",
                    n[0]
                );
                assert!(
                    (n[1] - 1.0).abs() < 1e-10,
                    "ny should be 1 at ({i},{j}), got {}",
                    n[1]
                );
                assert!(
                    n[2].abs() < 1e-10,
                    "nz should be 0 at ({i},{j}), got {}",
                    n[2]
                );
            }
        }
    }
    #[test]
    fn test_normal_at_tilted_plane_has_nonzero_components() {
        let hf = HeightField::from_fn(4, 4, 1.0, |c, _| c as f64);
        let n = heightfield_normal_at(&hf, 2, 2);
        assert!(
            n[0].abs() > 0.1,
            "tilted plane should have nonzero nx, got {}",
            n[0]
        );
        let len = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
        assert!(
            (len - 1.0).abs() < 1e-9,
            "normal should be unit length, got {len}"
        );
    }
    #[test]
    fn test_bilinear_interpolation_center() {
        let hf = HeightField::new(vec![0.0, 1.0, 2.0, 3.0], 2, 2, 1.0, 1.0);
        let h = heightfield_height_at_xy(&hf, 0.5, 0.5);
        assert!(
            (h - 1.5).abs() < 1e-9,
            "bilinear center should be 1.5, got {h}"
        );
    }
    #[test]
    fn test_bilinear_interpolation_corner() {
        let hf = HeightField::new(vec![0.0, 1.0, 2.0, 3.0], 2, 2, 1.0, 1.0);
        let h = heightfield_height_at_xy(&hf, 0.0, 0.0);
        assert!((h - 0.0).abs() < 1e-9, "corner (0,0) should be 0, got {h}");
    }
    #[test]
    fn test_bilinear_interpolation_clamps() {
        let hf = HeightField::from_fn(3, 3, 1.0, |_, _| 7.0);
        let h = heightfield_height_at_xy(&hf, 999.0, 999.0);
        assert!(
            (h - 7.0).abs() < 1e-9,
            "clamped height should be 7, got {h}"
        );
    }
}
#[cfg(test)]
mod tests_convenience {

    use crate::HeightField;

    #[test]
    fn test_aabb_flat_field() {
        let hf = HeightField::from_fn(4, 4, 2.0, |_, _| 5.0);
        let (mn, mx) = hf.aabb();
        assert!((mn[0] - 0.0).abs() < 1e-9);
        assert!((mn[1] - 5.0).abs() < 1e-9);
        assert!((mn[2] - 0.0).abs() < 1e-9);
        assert!((mx[0] - 6.0).abs() < 1e-9);
        assert!((mx[1] - 5.0).abs() < 1e-9);
        assert!((mx[2] - 6.0).abs() < 1e-9);
    }
    #[test]
    fn test_aabb_varying_heights() {
        let mut heights = vec![0.0f64; 3 * 3];
        heights[4] = 10.0;
        let hf = HeightField::new(heights, 3, 3, 1.0, 1.0);
        let (mn, mx) = hf.aabb();
        assert!((mn[1] - 0.0).abs() < 1e-9, "min height should be 0");
        assert!((mx[1] - 10.0).abs() < 1e-9, "max height should be 10");
    }
    #[test]
    fn test_aabb_contains_all_vertices() {
        let hf = HeightField::from_fn(5, 5, 1.0, |col, row| (col + row) as f64 * 0.5);
        let (mn, mx) = hf.aabb();
        for row in 0..hf.rows {
            for col in 0..hf.cols {
                let x = col as f64 * hf.scale_x;
                let y = hf.height_at(row, col);
                let z = row as f64 * hf.scale_z;
                assert!(x >= mn[0] && x <= mx[0], "x out of aabb");
                assert!(y >= mn[1] && y <= mx[1], "y out of aabb");
                assert!(z >= mn[2] && z <= mx[2], "z out of aabb");
            }
        }
    }
    #[test]
    fn test_height_at_xz_flat_field() {
        let hf = HeightField::from_fn(4, 4, 1.0, |_, _| 3.0);
        assert!((hf.height_at_xz(1.5, 1.5) - 3.0).abs() < 1e-9);
    }
    #[test]
    fn test_height_at_xz_interpolates() {
        let hf = HeightField::from_fn(4, 4, 1.0, |col, _row| col as f64);
        let h = hf.height_at_xz(0.5, 0.0);
        assert!((h - 0.5).abs() < 0.01, "expected ~0.5, got {h}");
    }
    #[test]
    fn test_normals_count() {
        let hf = HeightField::from_fn(3, 4, 1.0, |_, _| 0.0);
        let normals = hf.normals();
        assert_eq!(normals.len(), 3 * 4);
    }
    #[test]
    fn test_normals_unit_length() {
        let hf = HeightField::from_fn(4, 4, 1.0, |col, row| (col + row) as f64 * 0.3);
        for n in hf.normals() {
            let len = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
            assert!((len - 1.0).abs() < 1e-9, "normal not unit: len={len}");
        }
    }
    #[test]
    fn test_normals_flat_field_point_up() {
        let hf = HeightField::from_fn(4, 4, 1.0, |_, _| 0.0);
        for n in hf.normals() {
            assert!(
                (n[1] - 1.0).abs() < 1e-9,
                "flat field normal should point up"
            );
        }
    }
    #[test]
    fn test_tessellate_vertex_count() {
        let hf = HeightField::from_fn(3, 4, 1.0, |_, _| 0.0);
        let (verts, _) = hf.tessellate();
        assert_eq!(verts.len(), 3 * 4);
    }
    #[test]
    fn test_tessellate_triangle_count() {
        let hf = HeightField::from_fn(3, 4, 1.0, |_, _| 0.0);
        let (_, tris) = hf.tessellate();
        assert_eq!(tris.len(), 2 * 3 * 2);
    }
    #[test]
    fn test_tessellate_indices_in_range() {
        let hf = HeightField::from_fn(4, 4, 1.0, |col, row| (col as f64 + row as f64) * 0.1);
        let (verts, tris) = hf.tessellate();
        for tri in &tris {
            for &idx in tri {
                assert!(
                    idx < verts.len(),
                    "index {idx} out of range {}",
                    verts.len()
                );
            }
        }
    }
    #[test]
    fn test_tessellate_vertex_positions() {
        let hf = HeightField::from_fn(3, 3, 2.0, |col, row| (col + row) as f64);
        let (verts, _) = hf.tessellate();
        let v = verts[5];
        assert!((v[0] - 4.0).abs() < 1e-9, "x should be col*scale_x = 4");
        assert!((v[2] - 2.0).abs() < 1e-9, "z should be row*scale_z = 2");
        assert!((v[1] - 3.0).abs() < 1e-9, "y = col+row = 3");
    }
    #[test]
    fn test_ray_intersect_flat_hit() {
        let hf = HeightField::from_fn(5, 5, 1.0, |_, _| 0.0);
        let origin = [2.0, 10.0, 2.0];
        let dir = [0.0, -1.0, 0.0];
        let result = hf.ray_intersect(origin, dir);
        assert!(result.is_some(), "should hit flat terrain");
        let (t, _normal) = result.unwrap();
        assert!(t > 0.0 && t < 20.0, "t should be positive, got {t}");
    }
    #[test]
    fn test_ray_intersect_miss_above_terrain() {
        let hf = HeightField::from_fn(5, 5, 1.0, |_, _| 0.0);
        let origin = [2.0, 1.0, 2.0];
        let dir = [0.0, 1.0, 0.0];
        let result = hf.ray_intersect(origin, dir);
        assert!(result.is_none(), "upward ray should miss terrain");
    }
    #[test]
    fn test_ray_intersect_normal_is_unit() {
        let hf = HeightField::from_fn(5, 5, 1.0, |_, _| 0.0);
        let origin = [2.0, 10.0, 2.0];
        let dir = [0.0, -1.0, 0.0];
        if let Some((_t, n)) = hf.ray_intersect(origin, dir) {
            let len = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
            assert!((len - 1.0).abs() < 1e-6, "normal should be unit, len={len}");
        }
    }
}
#[cfg(test)]
mod tests_algorithms {

    use crate::HeightField;

    #[test]
    fn test_slope_at_flat_is_zero() {
        let hf = HeightField::from_fn(4, 4, 1.0, |_, _| 5.0);
        for row in 0..hf.rows {
            for col in 0..hf.cols {
                let s = hf.slope_at(col, row);
                assert!(s.abs() < 1e-10, "flat field slope should be 0, got {s}");
            }
        }
    }
    #[test]
    fn test_slope_at_tilted_plane() {
        let hf = HeightField::from_fn(5, 5, 1.0, |col, _| col as f64);
        let s = hf.slope_at(2, 2);
        assert!((s - 1.0).abs() < 1e-9, "slope should be 1.0, got {s}");
    }
    #[test]
    fn test_slope_map_count() {
        let hf = HeightField::from_fn(3, 4, 1.0, |_, _| 0.0);
        let sm = hf.slope_map();
        assert_eq!(sm.len(), 3 * 4);
    }
    #[test]
    fn test_curvature_at_flat_is_zero() {
        let hf = HeightField::from_fn(5, 5, 1.0, |_, _| 3.0);
        for row in 1..hf.rows - 1 {
            for col in 1..hf.cols - 1 {
                let k = hf.curvature_at(col, row);
                assert!(k.abs() < 1e-10, "flat curvature should be 0, got {k}");
            }
        }
    }
    #[test]
    fn test_curvature_map_count() {
        let hf = HeightField::from_fn(4, 5, 1.0, |_, _| 0.0);
        let cm = hf.curvature_map();
        assert_eq!(cm.len(), 4 * 5);
    }
    #[test]
    fn test_curvature_at_peak_is_negative() {
        let center = 2usize;
        let hf = HeightField::from_fn(5, 5, 1.0, |col, row| {
            let dc = col as f64 - center as f64;
            let dr = row as f64 - center as f64;
            -(dc * dc + dr * dr)
        });
        let k = hf.curvature_at(center, center);
        assert!(k < 0.0, "peak curvature should be negative, got {k}");
    }
    #[test]
    fn test_closest_vertex_flat_field() {
        let hf = HeightField::from_fn(4, 4, 1.0, |_, _| 0.0);
        let q = [2.0, 100.0, 1.0];
        let (pt, row, col) = hf.closest_vertex(q);
        assert_eq!(row, 1, "expected row=1, got {row}");
        assert_eq!(col, 2, "expected col=2, got {col}");
        assert!((pt[0] - 2.0).abs() < 1e-9);
        assert!((pt[2] - 1.0).abs() < 1e-9);
    }
    #[test]
    fn test_resample_preserves_extents() {
        let hf = HeightField::from_fn(4, 4, 1.0, |col, row| (col + row) as f64 * 0.5);
        let resampled = hf.resample(7, 7);
        let orig_max_x = (hf.cols - 1) as f64 * hf.scale_x;
        let new_max_x = (resampled.cols - 1) as f64 * resampled.scale_x;
        assert!((orig_max_x - new_max_x).abs() < 1e-9, "x extent mismatch");
        assert_eq!(resampled.rows, 7);
        assert_eq!(resampled.cols, 7);
    }
    #[test]
    fn test_resample_flat_field_stays_flat() {
        let hf = HeightField::from_fn(3, 3, 1.0, |_, _| 7.0);
        let r = hf.resample(5, 5);
        for &h in &r.heights {
            assert!(
                (h - 7.0).abs() < 1e-9,
                "flat field should stay flat, got {h}"
            );
        }
    }
    #[test]
    fn test_hydraulic_erode_lowers_peaks() {
        let mut hf = HeightField::from_fn(
            5,
            5,
            1.0,
            |col, row| {
                if col == 2 && row == 2 { 10.0 } else { 0.0 }
            },
        );
        let peak_before = hf.height_at(2, 2);
        hf.hydraulic_erode(0.3, 3);
        let peak_after = hf.height_at(2, 2);
        assert!(
            peak_after < peak_before,
            "hydraulic erosion should lower peak: {peak_after} >= {peak_before}"
        );
    }
    #[test]
    fn test_hydraulic_erode_flat_no_change() {
        let mut hf = HeightField::from_fn(4, 4, 1.0, |_, _| 5.0);
        let before: Vec<f64> = hf.heights.clone();
        hf.hydraulic_erode(0.3, 5);
        for (a, b) in before.iter().zip(hf.heights.iter()) {
            assert!((a - b).abs() < 1e-10, "flat field should not change");
        }
    }
    #[test]
    fn test_flow_accumulation_length() {
        let hf = HeightField::from_fn(4, 4, 1.0, |col, row| (col + row) as f64);
        let fa = hf.flow_accumulation();
        assert_eq!(fa.len(), 4 * 4, "flow accumulation map size mismatch");
    }
    #[test]
    fn test_flow_accumulation_all_positive() {
        let hf = HeightField::from_fn(4, 4, 1.0, |col, row| (col + row) as f64);
        let fa = hf.flow_accumulation();
        for &v in &fa {
            assert!(v >= 1.0, "each cell should accumulate at least 1, got {v}");
        }
    }
    #[test]
    fn test_clamp_heights() {
        let mut hf = HeightField::new(vec![-5.0, 3.0, 7.0, 10.0], 2, 2, 1.0, 1.0);
        hf.clamp_heights(0.0, 5.0);
        assert!((hf.heights[0] - 0.0).abs() < 1e-10, "should clamp -5 to 0");
        assert!((hf.heights[1] - 3.0).abs() < 1e-10, "3 unchanged");
        assert!((hf.heights[2] - 5.0).abs() < 1e-10, "7 clamped to 5");
        assert!((hf.heights[3] - 5.0).abs() < 1e-10, "10 clamped to 5");
    }
    #[test]
    fn test_scale_heights() {
        let mut hf = HeightField::new(vec![1.0, 2.0, 3.0, 4.0], 2, 2, 1.0, 1.0);
        hf.scale_heights(3.0);
        assert!((hf.heights[0] - 3.0).abs() < 1e-10);
        assert!((hf.heights[3] - 12.0).abs() < 1e-10);
    }
    #[test]
    fn test_offset_heights() {
        let mut hf = HeightField::new(vec![0.0, 1.0, 2.0, 3.0], 2, 2, 1.0, 1.0);
        hf.offset_heights(10.0);
        assert!((hf.heights[0] - 10.0).abs() < 1e-10);
        assert!((hf.heights[3] - 13.0).abs() < 1e-10);
    }
    #[test]
    fn test_normalize_heights() {
        let mut hf = HeightField::new(vec![0.0, 5.0, 10.0, 20.0], 2, 2, 1.0, 1.0);
        hf.normalize_heights();
        assert!((hf.min_height() - 0.0).abs() < 1e-10, "min should be 0");
        assert!((hf.max_height() - 1.0).abs() < 1e-10, "max should be 1");
    }
    #[test]
    fn test_normalize_uniform_heights() {
        let mut hf = HeightField::from_fn(3, 3, 1.0, |_, _| 7.0);
        hf.normalize_heights();
        for &h in &hf.heights {
            assert!(
                (h - 0.0).abs() < 1e-10,
                "uniform field should normalize to 0"
            );
        }
    }
    #[test]
    fn test_invert_heights() {
        let mut hf = HeightField::new(vec![0.0, 5.0, 10.0, 20.0], 2, 2, 1.0, 1.0);
        hf.invert_heights();
        assert!(
            (hf.min_height() - 0.0).abs() < 1e-10,
            "min should still be 0"
        );
        assert!(
            (hf.max_height() - 20.0).abs() < 1e-10,
            "max should still be 20"
        );
        assert!((hf.heights[0] - 20.0).abs() < 1e-10);
        assert!((hf.heights[3] - 0.0).abs() < 1e-10);
    }
    #[test]
    fn test_mean_height_flat() {
        let hf = HeightField::from_fn(4, 4, 1.0, |_, _| 3.5);
        assert!((hf.mean_height() - 3.5).abs() < 1e-10);
    }
    #[test]
    fn test_mean_height_varying() {
        let hf = HeightField::new(vec![0.0, 1.0, 2.0, 3.0], 2, 2, 1.0, 1.0);
        assert!((hf.mean_height() - 1.5).abs() < 1e-10);
    }
    #[test]
    fn test_height_variance_flat_is_zero() {
        let hf = HeightField::from_fn(4, 4, 1.0, |_, _| 5.0);
        assert!(
            hf.height_variance() < 1e-10,
            "flat field variance should be 0"
        );
    }
    #[test]
    fn test_height_variance_positive_for_varied() {
        let hf = HeightField::new(vec![0.0, 10.0, 0.0, 10.0], 2, 2, 1.0, 1.0);
        assert!(
            hf.height_variance() > 0.0,
            "varied field should have positive variance"
        );
    }
    #[test]
    fn test_count_peaks_flat() {
        let hf = HeightField::from_fn(4, 4, 1.0, |_, _| 5.0);
        let peaks = hf.count_peaks();
        assert_eq!(peaks, 0, "flat field should have no peaks");
    }
    #[test]
    fn test_count_peaks_single_spike() {
        let hf = HeightField::from_fn(
            5,
            5,
            1.0,
            |col, row| {
                if col == 2 && row == 2 { 10.0 } else { 0.0 }
            },
        );
        let peaks = hf.count_peaks();
        assert_eq!(peaks, 1, "single spike should give 1 peak, got {peaks}");
    }
    #[test]
    fn test_count_peaks_two_spikes() {
        let hf = HeightField::from_fn(7, 7, 1.0, |col, row| {
            if (col == 1 && row == 1) || (col == 5 && row == 5) {
                10.0
            } else {
                0.0
            }
        });
        let peaks = hf.count_peaks();
        assert_eq!(peaks, 2, "two spikes should give 2 peaks, got {peaks}");
    }
}
