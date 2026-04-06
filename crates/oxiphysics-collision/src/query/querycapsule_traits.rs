// # QueryCapsule - Trait Implementations
//
// This module contains trait implementations for `QueryCapsule`.
//
// ## Implemented Traits
//
// - `ShapeQuery`
//
// 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::*;
use super::types::{PointQueryResult, QueryCapsule, QuerySphere, RayCastResult};

impl ShapeQuery for QueryCapsule {
    fn ray_cast(
        &self,
        ray_origin: [f64; 3],
        ray_dir: [f64; 3],
        max_toi: f64,
    ) -> Option<RayCastResult> {
        let ab = sub(self.p1, self.p0);
        let ao = sub(ray_origin, self.p0);
        let r = self.radius;
        let d = ray_dir;
        let ab_len = norm(ab);
        if ab_len < 1e-15 {
            let sphere = QuerySphere {
                center: self.p0,
                radius: r,
            };
            return sphere.ray_cast(ray_origin, ray_dir, max_toi);
        }
        let axis = scale(ab, 1.0 / ab_len);
        let d_par = dot(d, axis);
        let d_perp = sub(d, scale(axis, d_par));
        let ao_par = dot(ao, axis);
        let ao_perp = sub(ao, scale(axis, ao_par));
        let a = dot(d_perp, d_perp);
        let b = dot(ao_perp, d_perp);
        let c_coef = dot(ao_perp, ao_perp) - r * r;
        let mut best_toi = f64::INFINITY;
        let mut best_point = [0.0; 3];
        let mut best_normal = [0.0; 3];
        if a > 1e-30 {
            let disc = b * b - a * c_coef;
            if disc >= 0.0 {
                let sqrt_disc = disc.sqrt();
                for &sign in &[-1.0_f64, 1.0] {
                    let t = (-b + sign * sqrt_disc) / a;
                    if t < 0.0 || t > max_toi {
                        continue;
                    }
                    let hit = add(ray_origin, scale(d, t));
                    let hit_par = dot(sub(hit, self.p0), axis);
                    if hit_par < 0.0 || hit_par > ab_len {
                        continue;
                    }
                    if t < best_toi {
                        best_toi = t;
                        best_point = hit;
                        let proj = add(self.p0, scale(axis, hit_par));
                        best_normal = normalize(sub(hit, proj));
                    }
                }
            }
        }
        for &cap_center in &[self.p0, self.p1] {
            let sphere = QuerySphere {
                center: cap_center,
                radius: r,
            };
            if let Some(res) = sphere.ray_cast(ray_origin, d, max_toi) {
                let hit_par = dot(sub(res.point, self.p0), axis);
                let valid = if cap_center == self.p0 {
                    hit_par <= 0.0
                } else {
                    hit_par >= ab_len
                };
                if valid && res.toi < best_toi {
                    best_toi = res.toi;
                    best_point = res.point;
                    best_normal = res.normal;
                }
            }
        }
        if best_toi.is_finite() && best_toi <= max_toi {
            Some(RayCastResult {
                toi: best_toi,
                normal: best_normal,
                point: best_point,
            })
        } else {
            None
        }
    }
    fn point_query(&self, point: [f64; 3]) -> PointQueryResult {
        let (seg_pt, _t) = self.closest_point_on_segment(point);
        let diff = sub(point, seg_pt);
        let dist_from_axis = norm(diff);
        let signed_dist = dist_from_axis - self.radius;
        let inside = signed_dist < 0.0;
        let closest = if dist_from_axis < 1e-15 {
            add(seg_pt, [self.radius, 0.0, 0.0])
        } else {
            add(seg_pt, scale(normalize(diff), self.radius))
        };
        PointQueryResult {
            point: closest,
            distance: signed_dist,
            inside,
        }
    }
    fn aabb(&self) -> ([f64; 3], [f64; 3]) {
        let r = self.radius;
        let lo = [
            self.p0[0].min(self.p1[0]) - r,
            self.p0[1].min(self.p1[1]) - r,
            self.p0[2].min(self.p1[2]) - r,
        ];
        let hi = [
            self.p0[0].max(self.p1[0]) + r,
            self.p0[1].max(self.p1[1]) + r,
            self.p0[2].max(self.p1[2]) + r,
        ];
        (lo, hi)
    }
}
