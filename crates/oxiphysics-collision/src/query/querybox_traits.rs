// # QueryBox - Trait Implementations
//
// This module contains trait implementations for `QueryBox`.
//
// ## Implemented Traits
//
// - `ShapeQuery`
//
// 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::*;
use super::types::{PointQueryResult, QueryBox, RayCastResult};

impl ShapeQuery for QueryBox {
    fn ray_cast(
        &self,
        ray_origin: [f64; 3],
        ray_dir: [f64; 3],
        max_toi: f64,
    ) -> Option<RayCastResult> {
        let lo = sub(self.center, self.half_extents);
        let hi = add(self.center, self.half_extents);
        let mut t_min = 0.0_f64;
        let mut t_max = max_toi;
        let mut hit_axis = 0usize;
        let mut hit_sign = 1.0_f64;
        for i in 0..3 {
            let d = ray_dir[i];
            let o = ray_origin[i];
            if d.abs() < 1e-15 {
                if o < lo[i] || o > hi[i] {
                    return None;
                }
            } else {
                let inv_d = 1.0 / d;
                let mut t1 = (lo[i] - o) * inv_d;
                let mut t2 = (hi[i] - o) * inv_d;
                let sign = if t1 <= t2 { -1.0_f64 } else { 1.0_f64 };
                if t1 > t2 {
                    core::mem::swap(&mut t1, &mut t2);
                }
                if t1 > t_min {
                    t_min = t1;
                    hit_axis = i;
                    hit_sign = sign;
                }
                if t2 < t_max {
                    t_max = t2;
                }
                if t_min > t_max {
                    return None;
                }
            }
        }
        if t_min < 0.0 || t_min > max_toi {
            return None;
        }
        let point = add(ray_origin, scale(ray_dir, t_min));
        let mut normal = [0.0; 3];
        normal[hit_axis] = hit_sign;
        Some(RayCastResult {
            toi: t_min,
            normal,
            point,
        })
    }
    fn point_query(&self, point: [f64; 3]) -> PointQueryResult {
        let lo = sub(self.center, self.half_extents);
        let hi = add(self.center, self.half_extents);
        let clamped = [
            clamp(point[0], lo[0], hi[0]),
            clamp(point[1], lo[1], hi[1]),
            clamp(point[2], lo[2], hi[2]),
        ];
        let inside = point[0] >= lo[0]
            && point[0] <= hi[0]
            && point[1] >= lo[1]
            && point[1] <= hi[1]
            && point[2] >= lo[2]
            && point[2] <= hi[2];
        if inside {
            let depths = [
                hi[0] - point[0],
                point[0] - lo[0],
                hi[1] - point[1],
                point[1] - lo[1],
                hi[2] - point[2],
                point[2] - lo[2],
            ];
            let min_depth = depths.iter().cloned().fold(f64::INFINITY, f64::min);
            PointQueryResult {
                point: clamped,
                distance: -min_depth,
                inside: true,
            }
        } else {
            let dist = norm(sub(point, clamped));
            PointQueryResult {
                point: clamped,
                distance: dist,
                inside: false,
            }
        }
    }
    fn aabb(&self) -> ([f64; 3], [f64; 3]) {
        let lo = sub(self.center, self.half_extents);
        let hi = add(self.center, self.half_extents);
        (lo, hi)
    }
}
