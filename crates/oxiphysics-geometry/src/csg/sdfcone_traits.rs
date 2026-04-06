//! # SdfCone - Trait Implementations
//!
//! This module contains trait implementations for `SdfCone`.
//!
//! ## Implemented Traits
//!
//! - `ImplicitSurface`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::ImplicitSurface;
#[allow(unused_imports)]
use super::functions::*;
use super::functions::{length, normalize};
use super::types::SdfCone;

impl ImplicitSurface for SdfCone {
    fn sdf(&self, p: [f64; 3]) -> f64 {
        let dx = p[0] - self.apex[0];
        let dy = p[1] - self.apex[1];
        let dz = p[2] - self.apex[2];
        let r = (dx * dx + dz * dz).sqrt();
        let h = dy;
        let sin_a = self.half_angle.sin();
        let cos_a = self.half_angle.cos();
        let q = [(r * r + h * h).sqrt(), 0.0_f64, 0.0];
        let c = [
            sin_a * q[0] - cos_a * q[1],
            cos_a * q[0] + sin_a * q[1],
            0.0,
        ];
        let c_clamped = [c[0].max(0.0), c[1].clamp(0.0, self.height), 0.0];
        let cone_sdf = c[1].max(0.0)
            * (-1.0_f64)
                .min(1.0)
                .min(if c[0] < 0.0 && c[1] < 0.0 { -1.0 } else { 1.0 });
        let _ = (c_clamped, cone_sdf);
        let cone_r_at_h = h.max(0.0).min(self.height) * self.half_angle.tan();
        let excess_r = r - cone_r_at_h;
        let h_dist = h - self.height;
        if h < 0.0 {
            length([dx, dy, dz])
        } else if h > self.height {
            let base_r = self.height * self.half_angle.tan();
            let dr = (r - base_r).max(0.0);
            (dr * dr + h_dist * h_dist).sqrt()
        } else {
            let side_normal_r = cos_a;
            let side_normal_h = -sin_a;
            excess_r * side_normal_r + (h - h * 1.0) * side_normal_h
        }
    }
    fn gradient(&self, p: [f64; 3]) -> [f64; 3] {
        const EPS: f64 = 1e-5;
        let gx = self.sdf([p[0] + EPS, p[1], p[2]]) - self.sdf([p[0] - EPS, p[1], p[2]]);
        let gy = self.sdf([p[0], p[1] + EPS, p[2]]) - self.sdf([p[0], p[1] - EPS, p[2]]);
        let gz = self.sdf([p[0], p[1], p[2] + EPS]) - self.sdf([p[0], p[1], p[2] - EPS]);
        normalize([gx, gy, gz])
    }
}
