// # QuerySphere - Trait Implementations
//
// This module contains trait implementations for `QuerySphere`.
//
// ## Implemented Traits
//
// - `ShapeQuery`
//
// 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::*;
use super::types::{PointQueryResult, QuerySphere, RayCastResult};

impl ShapeQuery for QuerySphere {
    fn ray_cast(
        &self,
        ray_origin: [f64; 3],
        ray_dir: [f64; 3],
        max_toi: f64,
    ) -> Option<RayCastResult> {
        let m = sub(ray_origin, self.center);
        let b = dot(m, ray_dir);
        let c = dot(m, m) - self.radius * self.radius;
        let disc = b * b - c;
        if disc < 0.0 {
            return None;
        }
        let sqrt_disc = disc.sqrt();
        let t0 = -b - sqrt_disc;
        let t1 = -b + sqrt_disc;
        let toi = if t0 >= 0.0 {
            t0
        } else if t1 >= 0.0 {
            t1
        } else {
            return None;
        };
        if toi > max_toi {
            return None;
        }
        let point = add(ray_origin, scale(ray_dir, toi));
        let normal = normalize(sub(point, self.center));
        Some(RayCastResult { toi, normal, point })
    }
    fn point_query(&self, point: [f64; 3]) -> PointQueryResult {
        let d = sub(point, self.center);
        let dist_from_center = norm(d);
        let signed_dist = dist_from_center - self.radius;
        let inside = signed_dist < 0.0;
        let closest = if dist_from_center < 1e-15 {
            add(self.center, [self.radius, 0.0, 0.0])
        } else {
            add(self.center, scale(normalize(d), self.radius))
        };
        PointQueryResult {
            point: closest,
            distance: signed_dist,
            inside,
        }
    }
    fn aabb(&self) -> ([f64; 3], [f64; 3]) {
        let r = self.radius;
        let c = self.center;
        (
            [c[0] - r, c[1] - r, c[2] - r],
            [c[0] + r, c[1] + r, c[2] + r],
        )
    }
}
