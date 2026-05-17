// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Sphere-tracing ray marcher for SDF scenes.

use super::helpers::{add3, normalize3, scale3};
use super::operators::sdf_normal;

// ─────────────────────────────────────────────────────────────────────────────
// Ray Marching (sphere tracing)
// ─────────────────────────────────────────────────────────────────────────────

/// Result of a ray-march intersection test.
#[derive(Debug, Clone)]
pub struct RayMarchHit {
    /// Distance along the ray \[same units as SDF\].
    pub t: f64,
    /// Hit point.
    pub point: [f64; 3],
    /// Outward surface normal at the hit point.
    pub normal: [f64; 3],
    /// SDF value at the hit point (should be near zero).
    pub sdf_value: f64,
    /// Number of iterations used.
    pub steps: usize,
}

/// Sphere-tracing ray marcher for SDF scenes.
#[derive(Debug, Clone)]
pub struct RayMarcher {
    /// Maximum ray-march distance.
    pub t_max: f64,
    /// Surface hit tolerance.
    pub tolerance: f64,
    /// Maximum sphere-tracing iterations.
    pub max_steps: usize,
    /// Step scale factor (< 1 for over-relaxation safety).
    pub step_scale: f64,
}

impl RayMarcher {
    /// Construct a ray marcher with default parameters.
    pub fn new() -> Self {
        Self {
            t_max: 100.0,
            tolerance: 1e-4,
            max_steps: 256,
            step_scale: 0.95,
        }
    }

    /// Construct a ray marcher with custom parameters.
    pub fn with_params(t_max: f64, tolerance: f64, max_steps: usize, step_scale: f64) -> Self {
        Self {
            t_max,
            tolerance,
            max_steps,
            step_scale,
        }
    }

    /// March a ray from `origin` in direction `dir` (should be unit length)
    /// through the SDF `f`.
    ///
    /// Returns `Some(RayMarchHit)` if a surface is found, `None` otherwise.
    pub fn march<F>(&self, f: &F, origin: [f64; 3], dir: [f64; 3]) -> Option<RayMarchHit>
    where
        F: Fn([f64; 3]) -> f64,
    {
        let d = normalize3(dir);
        let mut t = 0.0;
        for step in 0..self.max_steps {
            let p = add3(origin, scale3(d, t));
            let sdf = f(p);
            if sdf.abs() < self.tolerance {
                let normal = sdf_normal(f, p, 1e-4);
                return Some(RayMarchHit {
                    t,
                    point: p,
                    normal,
                    sdf_value: sdf,
                    steps: step + 1,
                });
            }
            if t > self.t_max {
                break;
            }
            t += sdf.abs() * self.step_scale;
        }
        None
    }

    /// Cast a shadow ray: returns true if the path from `p` toward `light_dir`
    /// is unobstructed within distance `max_dist`.
    pub fn shadow<F>(&self, f: &F, p: [f64; 3], light_dir: [f64; 3], max_dist: f64) -> f64
    where
        F: Fn([f64; 3]) -> f64,
    {
        let d = normalize3(light_dir);
        let mut t = 0.01; // offset to avoid self-intersection
        let mut shadow = 1.0_f64;
        for _ in 0..self.max_steps {
            if t >= max_dist {
                break;
            }
            let q = add3(p, scale3(d, t));
            let h = f(q);
            if h < self.tolerance {
                return 0.0;
            }
            shadow = shadow.min(8.0 * h / t);
            t += h;
        }
        shadow.clamp(0.0, 1.0)
    }

    /// Compute ambient occlusion at surface point `p` with normal `n`.
    pub fn ambient_occlusion<F>(&self, f: &F, p: [f64; 3], n: [f64; 3], step: f64) -> f64
    where
        F: Fn([f64; 3]) -> f64,
    {
        let mut occ = 0.0;
        let mut scale = 1.0;
        for i in 0..5 {
            let dist = (i + 1) as f64 * step;
            let q = add3(p, scale3(n, dist));
            occ += scale * (dist - f(q));
            scale *= 0.5;
        }
        1.0 - occ.clamp(0.0, 1.0)
    }
}

impl Default for RayMarcher {
    fn default() -> Self {
        Self::new()
    }
}
