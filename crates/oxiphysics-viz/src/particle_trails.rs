// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Particle trail, pathline, and streamline visualization.
//!
//! Provides data structures and algorithms for rendering motion history of
//! particles, integrating pathlines through velocity fields, and computing
//! streamlines for steady flows.  All geometry is returned as plain Rust data
//! structures — no GPU dependency.

#![allow(dead_code)]

use std::collections::HashMap;

// ─────────────────────────────────────────────────────────────────────────────
// TrailPoint
// ─────────────────────────────────────────────────────────────────────────────

/// A single sample recorded along a particle trail.
#[derive(Debug, Clone, PartialEq)]
pub struct TrailPoint {
    /// World-space position \[x, y, z\].
    pub position: [f64; 3],
    /// Simulation time at which this sample was recorded.
    pub time: f64,
    /// RGBA color of this sample (components in \[0, 1\]).
    pub color: [f64; 4],
    /// Rendered point size (pixels or NDC-relative units).
    pub size: f64,
}

impl TrailPoint {
    /// Construct a new `TrailPoint`.
    pub fn new(position: [f64; 3], time: f64, color: [f64; 4], size: f64) -> Self {
        Self {
            position,
            time,
            color,
            size,
        }
    }

    /// Construct with default white color and unit size.
    pub fn simple(position: [f64; 3], time: f64) -> Self {
        Self::new(position, time, [1.0, 1.0, 1.0, 1.0], 1.0)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// ParticleTrail
// ─────────────────────────────────────────────────────────────────────────────

/// A time-ordered sequence of [`TrailPoint`]s representing a particle's history.
///
/// Optionally limits the total number of stored points (`max_length`) and can
/// apply an alpha-fade effect based on sample age.
#[derive(Debug, Clone)]
pub struct ParticleTrail {
    /// Stored trail samples, oldest first.
    pub points: Vec<TrailPoint>,
    /// Maximum number of points to keep. `0` means unlimited.
    pub max_length: usize,
    /// When `true`, the alpha channel of each point is multiplied by its
    /// relative age (newest = 1.0, oldest = 0.0).
    pub fade_alpha: bool,
}

impl ParticleTrail {
    /// Create an empty trail.
    ///
    /// - `max_length`: maximum number of recorded points (0 = unlimited).
    /// - `fade_alpha`: whether to fade older points toward transparency.
    pub fn new(max_length: usize, fade_alpha: bool) -> Self {
        Self {
            points: Vec::new(),
            max_length,
            fade_alpha,
        }
    }

    /// Append a new point to the trail.
    ///
    /// If `max_length > 0` and the trail is at capacity, the oldest point is
    /// removed before inserting the new one.
    pub fn add_point(&mut self, point: TrailPoint) {
        if self.max_length > 0 && self.points.len() >= self.max_length {
            self.points.remove(0);
        }
        self.points.push(point);

        if self.fade_alpha {
            self.apply_fade();
        }
    }

    /// Remove all points older than `current_time - max_age`.
    pub fn trim_old(&mut self, current_time: f64, max_age: f64) {
        let cutoff = current_time - max_age;
        self.points.retain(|p| p.time >= cutoff);
        if self.fade_alpha {
            self.apply_fade();
        }
    }

    /// Re-apply the alpha fade based on relative age (newest = alpha 1.0).
    fn apply_fade(&mut self) {
        let n = self.points.len();
        if n == 0 {
            return;
        }
        for (i, pt) in self.points.iter_mut().enumerate() {
            let t = if n > 1 {
                i as f64 / (n - 1) as f64
            } else {
                1.0
            };
            pt.color[3] = t;
        }
    }

    /// Return the number of points in the trail.
    pub fn len(&self) -> usize {
        self.points.len()
    }

    /// Return `true` if the trail contains no points.
    pub fn is_empty(&self) -> bool {
        self.points.is_empty()
    }

    /// Clear all points.
    pub fn clear(&mut self) {
        self.points.clear();
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// TrailStyle
// ─────────────────────────────────────────────────────────────────────────────

/// Controls how a trail is colored and sized.
#[derive(Debug, Clone)]
pub struct TrailStyle {
    /// Which attribute drives the color of each segment.
    pub color_by: ColorBy,
    /// Which attribute drives the width of each segment.
    pub width_by: WidthBy,
    /// Alpha multiplier applied as a function of relative age (1.0 = no extra
    /// decay, 0.0 = fully transparent at the tail).
    pub alpha_decay: f64,
    /// Base RGBA color used when `color_by == ColorBy::Constant`.
    pub base_color: [f64; 4],
    /// Base width used when `width_by == WidthBy::Constant`.
    pub base_width: f64,
}

/// Attribute used to determine per-sample color.
#[derive(Debug, Clone, PartialEq)]
pub enum ColorBy {
    /// Color by relative sample age (0 = oldest, 1 = newest).
    Age,
    /// Color by instantaneous speed (requires velocity information).
    Speed,
    /// Color by particle identifier (hash-based).
    Id,
    /// Constant color for all samples.
    Constant,
}

/// Attribute used to determine per-sample width.
#[derive(Debug, Clone, PartialEq)]
pub enum WidthBy {
    /// Constant width.
    Constant,
    /// Width proportional to particle speed.
    Speed,
    /// Width proportional to relative age (thicker = newer).
    Age,
}

impl Default for TrailStyle {
    fn default() -> Self {
        Self {
            color_by: ColorBy::Age,
            width_by: WidthBy::Constant,
            alpha_decay: 1.0,
            base_color: [1.0, 1.0, 1.0, 1.0],
            base_width: 1.0,
        }
    }
}

impl TrailStyle {
    /// Compute the effective alpha for a sample with relative age `t ∈ [0, 1]`
    /// (0 = oldest, 1 = newest).
    pub fn alpha_at(&self, t: f64) -> f64 {
        t.powf(self.alpha_decay.max(0.0))
    }

    /// Compute the effective width for a sample with relative age `t` and
    /// speed `speed`.
    pub fn width_at(&self, t: f64, speed: f64) -> f64 {
        match self.width_by {
            WidthBy::Constant => self.base_width,
            WidthBy::Speed => self.base_width * speed.max(0.0),
            WidthBy::Age => self.base_width * t,
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// TrailRenderer
// ─────────────────────────────────────────────────────────────────────────────

/// Manages a collection of particle trails indexed by particle ID.
#[derive(Debug, Default)]
pub struct TrailRenderer {
    /// Per-particle trail data.
    pub trails: HashMap<usize, ParticleTrail>,
    /// Default maximum trail length for newly created trails.
    pub default_max_length: usize,
    /// Default fade-alpha setting for newly created trails.
    pub default_fade_alpha: bool,
}

impl TrailRenderer {
    /// Create a new `TrailRenderer` with specified defaults.
    pub fn new(default_max_length: usize, default_fade_alpha: bool) -> Self {
        Self {
            trails: HashMap::new(),
            default_max_length,
            default_fade_alpha,
        }
    }

    /// Record a new position sample for particle `particle_id`.
    ///
    /// Creates a trail entry if this is the first sample for the given ID.
    pub fn update(&mut self, particle_id: usize, pos: [f64; 3], time: f64) {
        let trail = self.trails.entry(particle_id).or_insert_with(|| {
            ParticleTrail::new(self.default_max_length, self.default_fade_alpha)
        });
        trail.add_point(TrailPoint::simple(pos, time));
    }

    /// Update a particle with an explicit color and size.
    pub fn update_colored(
        &mut self,
        particle_id: usize,
        pos: [f64; 3],
        time: f64,
        color: [f64; 4],
        size: f64,
    ) {
        let trail = self.trails.entry(particle_id).or_insert_with(|| {
            ParticleTrail::new(self.default_max_length, self.default_fade_alpha)
        });
        trail.add_point(TrailPoint::new(pos, time, color, size));
    }

    /// Return a reference to the trail for the given particle, or `None` if
    /// no data has been recorded yet.
    pub fn get_trail(&self, id: usize) -> Option<&ParticleTrail> {
        self.trails.get(&id)
    }

    /// Return a mutable reference to the trail for the given particle.
    pub fn get_trail_mut(&mut self, id: usize) -> Option<&mut ParticleTrail> {
        self.trails.get_mut(&id)
    }

    /// Remove the trail for a particle.
    pub fn remove_trail(&mut self, id: usize) {
        self.trails.remove(&id);
    }

    /// Trim all trails, discarding samples older than `current_time - max_age`.
    pub fn trim_all(&mut self, current_time: f64, max_age: f64) {
        for trail in self.trails.values_mut() {
            trail.trim_old(current_time, max_age);
        }
    }

    /// Return the total number of tracked particles.
    pub fn particle_count(&self) -> usize {
        self.trails.len()
    }

    /// Clear all trail data.
    pub fn clear(&mut self) {
        self.trails.clear();
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// PathlineIntegrator
// ─────────────────────────────────────────────────────────────────────────────

/// Integrates pathlines through a time-dependent (or steady) velocity field
/// using the classical 4th-order Runge–Kutta method.
pub struct PathlineIntegrator<F>
where
    F: Fn([f64; 3], f64) -> [f64; 3],
{
    /// Velocity field callback: `f(pos, t) -> velocity`.
    pub velocity_field: F,
    /// Integration time step size.
    pub step_size: f64,
}

impl<F> PathlineIntegrator<F>
where
    F: Fn([f64; 3], f64) -> [f64; 3],
{
    /// Create a new integrator with the given velocity field and step size.
    pub fn new(velocity_field: F, step_size: f64) -> Self {
        Self {
            velocity_field,
            step_size,
        }
    }

    /// Integrate from `start` at time `t0` to time `t_end` using RK4.
    ///
    /// Returns a `Vec` of positions including the starting point.
    /// If `t_end <= t0` or `step_size <= 0`, returns only the start point.
    pub fn integrate_rk4(&self, start: [f64; 3], t0: f64, t_end: f64) -> Vec<[f64; 3]> {
        if self.step_size <= 0.0 || t_end <= t0 {
            return vec![start];
        }

        let mut positions = vec![start];
        let mut pos = start;
        let mut t = t0;
        let h = self.step_size;

        while t < t_end - 1e-14 {
            let dt = h.min(t_end - t);
            pos = rk4_step(&self.velocity_field, pos, t, dt);
            t += dt;
            positions.push(pos);
        }

        positions
    }
}

/// Single RK4 step for a velocity field `f(pos, t)`.
fn rk4_step<F>(f: &F, pos: [f64; 3], t: f64, h: f64) -> [f64; 3]
where
    F: Fn([f64; 3], f64) -> [f64; 3],
{
    let k1 = f(pos, t);
    let k2 = f(vec3_add(pos, vec3_scale(k1, 0.5 * h)), t + 0.5 * h);
    let k3 = f(vec3_add(pos, vec3_scale(k2, 0.5 * h)), t + 0.5 * h);
    let k4 = f(vec3_add(pos, vec3_scale(k3, h)), t + h);

    [
        pos[0] + h / 6.0 * (k1[0] + 2.0 * k2[0] + 2.0 * k3[0] + k4[0]),
        pos[1] + h / 6.0 * (k1[1] + 2.0 * k2[1] + 2.0 * k3[1] + k4[1]),
        pos[2] + h / 6.0 * (k1[2] + 2.0 * k2[2] + 2.0 * k3[2] + k4[2]),
    ]
}

// ─────────────────────────────────────────────────────────────────────────────
// StreamlineSeeder
// ─────────────────────────────────────────────────────────────────────────────

/// Seeds streamlines from a set of starting positions and integrates them
/// through a steady velocity field.
pub struct StreamlineSeeder<F>
where
    F: Fn([f64; 3]) -> [f64; 3],
{
    /// Starting positions for streamlines.
    pub seed_positions: Vec<[f64; 3]>,
    /// Total arc-length (in time units) to integrate each streamline.
    pub integration_length: f64,
    /// Step size for the RK4 integrator.
    pub step_size: f64,
    /// Steady velocity field: `f(pos) -> velocity`.
    pub velocity_field: F,
}

impl<F> StreamlineSeeder<F>
where
    F: Fn([f64; 3]) -> [f64; 3],
{
    /// Create a new seeder.
    pub fn new(
        seed_positions: Vec<[f64; 3]>,
        integration_length: f64,
        step_size: f64,
        velocity_field: F,
    ) -> Self {
        Self {
            seed_positions,
            integration_length,
            step_size,
            velocity_field,
        }
    }

    /// Compute a streamline for each seed position.
    ///
    /// Returns one `Vec<[f64;3]>` per seed, containing the integrated path.
    pub fn compute_streamlines(&self) -> Vec<Vec<[f64; 3]>> {
        self.seed_positions
            .iter()
            .map(|&seed| {
                // Wrap the steady field as time-independent for the RK4 integrator.
                let integrator = PathlineIntegrator::new(
                    |p: [f64; 3], _t: f64| (self.velocity_field)(p),
                    self.step_size,
                );
                integrator.integrate_rk4(seed, 0.0, self.integration_length)
            })
            .collect()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// RibbonTrail
// ─────────────────────────────────────────────────────────────────────────────

/// A ribbon trail formed by two offset parallel curves around a particle path.
///
/// Each pair of consecutive center-line points is expanded by `width` in the
/// direction perpendicular to the velocity, producing a flat ribbon in 3-D.
#[derive(Debug, Clone)]
pub struct RibbonTrail {
    /// Center-line samples (oldest first).
    pub center_points: Vec<[f64; 3]>,
    /// Half-width of the ribbon (in world-space units).
    pub half_width: f64,
    /// Up-vector hint used to compute the ribbon's lateral direction when the
    /// velocity is nearly parallel to it.
    pub up_hint: [f64; 3],
}

impl RibbonTrail {
    /// Create an empty ribbon trail.
    pub fn new(half_width: f64) -> Self {
        Self {
            center_points: Vec::new(),
            half_width,
            up_hint: [0.0, 1.0, 0.0],
        }
    }

    /// Append a new center-line point.
    pub fn add_point(&mut self, pos: [f64; 3]) {
        self.center_points.push(pos);
    }

    /// Compute the ribbon geometry as a flat list of quad vertices.
    ///
    /// Returns a `Vec` of `([f64;3], [f64;3])` pairs, where each pair is the
    /// left and right edge positions for one center-line sample.  The length
    /// equals `center_points.len()`.
    ///
    /// When only one point exists, both edge vertices coincide at that point.
    pub fn compute_ribbon_geometry(&self) -> Vec<([f64; 3], [f64; 3])> {
        let n = self.center_points.len();
        if n == 0 {
            return Vec::new();
        }
        if n == 1 {
            let p = self.center_points[0];
            return vec![(p, p)];
        }

        let mut result = Vec::with_capacity(n);

        for i in 0..n {
            // Estimate tangent via finite differences.
            let tangent = if i == 0 {
                vec3_sub(self.center_points[1], self.center_points[0])
            } else if i == n - 1 {
                vec3_sub(self.center_points[n - 1], self.center_points[n - 2])
            } else {
                vec3_sub(self.center_points[i + 1], self.center_points[i - 1])
            };

            let tangent_norm = vec3_normalize(tangent);
            // Compute a lateral direction perpendicular to the tangent.
            let lateral = vec3_normalize(vec3_cross(tangent_norm, self.up_hint));
            // If tangent is parallel to up_hint, use an alternative.
            let lateral = if vec3_length(lateral) < 1e-8 {
                vec3_normalize(vec3_cross(tangent_norm, [1.0, 0.0, 0.0]))
            } else {
                lateral
            };

            let offset = vec3_scale(lateral, self.half_width);
            let p = self.center_points[i];
            let left = vec3_add(p, offset);
            let right = vec3_sub(p, offset);
            result.push((left, right));
        }

        result
    }

    /// Return the number of center-line points.
    pub fn len(&self) -> usize {
        self.center_points.len()
    }

    /// Return `true` if there are no center-line points.
    pub fn is_empty(&self) -> bool {
        self.center_points.is_empty()
    }

    /// Clear all center-line points.
    pub fn clear(&mut self) {
        self.center_points.clear();
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Vector helpers (no nalgebra — only [f64; 3])
// ─────────────────────────────────────────────────────────────────────────────

#[inline]
fn vec3_add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

#[inline]
fn vec3_sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

#[inline]
fn vec3_scale(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

#[inline]
fn vec3_length(a: [f64; 3]) -> f64 {
    (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt()
}

#[inline]
fn vec3_normalize(a: [f64; 3]) -> [f64; 3] {
    let len = vec3_length(a);
    if len < 1e-14 {
        a
    } else {
        vec3_scale(a, 1.0 / len)
    }
}

#[inline]
fn vec3_cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

#[inline]
fn vec3_dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[inline]
fn vec3_dist(a: [f64; 3], b: [f64; 3]) -> f64 {
    vec3_length(vec3_sub(a, b))
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // ── TrailPoint ────────────────────────────────────────────────────────────

    #[test]
    fn test_trail_point_new() {
        let pt = TrailPoint::new([1.0, 2.0, 3.0], 0.5, [1.0, 0.0, 0.0, 1.0], 2.0);
        assert_eq!(pt.position, [1.0, 2.0, 3.0]);
        assert_eq!(pt.time, 0.5);
        assert_eq!(pt.color[0], 1.0);
        assert_eq!(pt.size, 2.0);
    }

    #[test]
    fn test_trail_point_simple_defaults() {
        let pt = TrailPoint::simple([0.0, 0.0, 0.0], 1.0);
        assert_eq!(pt.color, [1.0, 1.0, 1.0, 1.0]);
        assert_eq!(pt.size, 1.0);
    }

    // ── ParticleTrail ─────────────────────────────────────────────────────────

    #[test]
    fn test_particle_trail_starts_empty() {
        let trail = ParticleTrail::new(10, false);
        assert!(trail.is_empty());
        assert_eq!(trail.len(), 0);
    }

    #[test]
    fn test_particle_trail_add_point() {
        let mut trail = ParticleTrail::new(10, false);
        trail.add_point(TrailPoint::simple([0.0, 0.0, 0.0], 0.0));
        assert_eq!(trail.len(), 1);
    }

    #[test]
    fn test_particle_trail_max_length_enforced() {
        let mut trail = ParticleTrail::new(3, false);
        for i in 0..5 {
            trail.add_point(TrailPoint::simple([i as f64, 0.0, 0.0], i as f64));
        }
        assert_eq!(trail.len(), 3, "trail should be capped at max_length=3");
        // Oldest points should have been dropped; last point should be t=4.
        assert_eq!(trail.points.last().unwrap().time, 4.0);
    }

    #[test]
    fn test_particle_trail_unlimited_when_max_zero() {
        let mut trail = ParticleTrail::new(0, false);
        for i in 0..20 {
            trail.add_point(TrailPoint::simple([i as f64, 0.0, 0.0], i as f64));
        }
        assert_eq!(trail.len(), 20);
    }

    #[test]
    fn test_particle_trail_fade_alpha_applied() {
        let mut trail = ParticleTrail::new(0, true);
        trail.add_point(TrailPoint::simple([0.0, 0.0, 0.0], 0.0));
        trail.add_point(TrailPoint::simple([1.0, 0.0, 0.0], 1.0));
        // With 2 points: oldest gets alpha 0.0, newest gets 1.0.
        assert!((trail.points[0].color[3] - 0.0).abs() < 1e-10);
        assert!((trail.points[1].color[3] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_particle_trail_trim_old_removes_stale_points() {
        let mut trail = ParticleTrail::new(0, false);
        for i in 0..10 {
            trail.add_point(TrailPoint::simple([i as f64, 0.0, 0.0], i as f64));
        }
        trail.trim_old(9.0, 3.0); // keep t >= 6
        for pt in &trail.points {
            assert!(
                pt.time >= 6.0,
                "stale point t={} should be removed",
                pt.time
            );
        }
    }

    #[test]
    fn test_particle_trail_trim_old_keeps_recent() {
        let mut trail = ParticleTrail::new(0, false);
        trail.add_point(TrailPoint::simple([0.0, 0.0, 0.0], 10.0));
        trail.trim_old(11.0, 5.0); // keep t >= 6 → 10 is kept
        assert_eq!(trail.len(), 1);
    }

    #[test]
    fn test_particle_trail_clear() {
        let mut trail = ParticleTrail::new(0, false);
        trail.add_point(TrailPoint::simple([0.0, 0.0, 0.0], 0.0));
        trail.clear();
        assert!(trail.is_empty());
    }

    // ── TrailRenderer ─────────────────────────────────────────────────────────

    #[test]
    fn test_trail_renderer_update_creates_trail() {
        let mut renderer = TrailRenderer::new(10, false);
        renderer.update(0, [1.0, 2.0, 3.0], 0.0);
        assert!(renderer.get_trail(0).is_some());
    }

    #[test]
    fn test_trail_renderer_update_appends_points() {
        let mut renderer = TrailRenderer::new(10, false);
        renderer.update(0, [0.0, 0.0, 0.0], 0.0);
        renderer.update(0, [1.0, 0.0, 0.0], 1.0);
        assert_eq!(renderer.get_trail(0).unwrap().len(), 2);
    }

    #[test]
    fn test_trail_renderer_multiple_particles() {
        let mut renderer = TrailRenderer::new(10, false);
        renderer.update(0, [0.0, 0.0, 0.0], 0.0);
        renderer.update(1, [1.0, 0.0, 0.0], 0.0);
        renderer.update(2, [2.0, 0.0, 0.0], 0.0);
        assert_eq!(renderer.particle_count(), 3);
    }

    #[test]
    fn test_trail_renderer_get_trail_missing_returns_none() {
        let renderer = TrailRenderer::new(10, false);
        assert!(renderer.get_trail(99).is_none());
    }

    #[test]
    fn test_trail_renderer_remove_trail() {
        let mut renderer = TrailRenderer::new(10, false);
        renderer.update(0, [0.0, 0.0, 0.0], 0.0);
        renderer.remove_trail(0);
        assert!(renderer.get_trail(0).is_none());
    }

    #[test]
    fn test_trail_renderer_trim_all() {
        let mut renderer = TrailRenderer::new(0, false);
        for i in 0..5 {
            renderer.update(0, [i as f64, 0.0, 0.0], i as f64);
        }
        renderer.trim_all(4.0, 2.0); // keep t >= 2
        let trail = renderer.get_trail(0).unwrap();
        for pt in &trail.points {
            assert!(pt.time >= 2.0);
        }
    }

    #[test]
    fn test_trail_renderer_clear() {
        let mut renderer = TrailRenderer::new(10, false);
        renderer.update(0, [0.0, 0.0, 0.0], 0.0);
        renderer.clear();
        assert_eq!(renderer.particle_count(), 0);
    }

    #[test]
    fn test_trail_renderer_update_colored() {
        let mut renderer = TrailRenderer::new(10, false);
        renderer.update_colored(5, [0.0, 0.0, 0.0], 1.0, [1.0, 0.0, 0.0, 1.0], 3.0);
        let trail = renderer.get_trail(5).unwrap();
        let pt = &trail.points[0];
        assert_eq!(pt.color[0], 1.0);
        assert_eq!(pt.size, 3.0);
    }

    // ── TrailStyle ────────────────────────────────────────────────────────────

    #[test]
    fn test_trail_style_alpha_at_endpoints() {
        let style = TrailStyle::default();
        assert!((style.alpha_at(0.0) - 0.0).abs() < 1e-10);
        assert!((style.alpha_at(1.0) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_trail_style_width_constant() {
        let style = TrailStyle {
            width_by: WidthBy::Constant,
            base_width: 2.5,
            ..Default::default()
        };
        assert!((style.width_at(0.5, 10.0) - 2.5).abs() < 1e-10);
    }

    #[test]
    fn test_trail_style_width_by_age() {
        let style = TrailStyle {
            width_by: WidthBy::Age,
            base_width: 2.0,
            ..Default::default()
        };
        let w = style.width_at(0.5, 0.0);
        assert!((w - 1.0).abs() < 1e-10);
    }

    // ── PathlineIntegrator (RK4) ───────────────────────────────────────────────

    #[test]
    fn test_rk4_uniform_field_x() {
        // v = (1, 0, 0) → integrate from 0 to 1 with h=0.1 → final x = 1.0
        let integrator = PathlineIntegrator::new(|_p, _t| [1.0, 0.0, 0.0], 0.1);
        let path = integrator.integrate_rk4([0.0, 0.0, 0.0], 0.0, 1.0);
        let last = path.last().unwrap();
        assert!(
            (last[0] - 1.0).abs() < 1e-10,
            "x should be 1.0: {}",
            last[0]
        );
        assert!(last[1].abs() < 1e-10);
        assert!(last[2].abs() < 1e-10);
    }

    #[test]
    fn test_rk4_returns_start_when_t_end_le_t0() {
        let integrator = PathlineIntegrator::new(|_p, _t| [1.0, 0.0, 0.0], 0.1);
        let path = integrator.integrate_rk4([1.0, 2.0, 3.0], 5.0, 5.0);
        assert_eq!(path.len(), 1);
        assert_eq!(path[0], [1.0, 2.0, 3.0]);
    }

    #[test]
    fn test_rk4_step_count() {
        // 10 steps of 0.1 from t=0 to t=1 → 11 points (start + 10 steps).
        let integrator = PathlineIntegrator::new(|_p, _t| [0.0, 1.0, 0.0], 0.1);
        let path = integrator.integrate_rk4([0.0, 0.0, 0.0], 0.0, 1.0);
        assert_eq!(path.len(), 11);
    }

    #[test]
    fn test_rk4_circular_field_radius_preserved() {
        // v = (-y, x, 0): circular motion; radius should be preserved.
        let integrator = PathlineIntegrator::new(|p: [f64; 3], _t| [-p[1], p[0], 0.0], 0.01);
        let path = integrator.integrate_rk4([1.0, 0.0, 0.0], 0.0, 1.0);
        let r_start = vec3_length(path[0]);
        let r_end = vec3_length(*path.last().unwrap());
        assert!(
            (r_end - r_start).abs() < 1e-4,
            "RK4 should preserve radius: start={r_start}, end={r_end}"
        );
    }

    #[test]
    fn test_rk4_negative_step_returns_start() {
        let integrator = PathlineIntegrator::new(|_p, _t| [1.0, 0.0, 0.0], -0.1);
        let path = integrator.integrate_rk4([0.0, 0.0, 0.0], 0.0, 1.0);
        assert_eq!(path.len(), 1);
    }

    #[test]
    fn test_rk4_zero_velocity_field() {
        // v = (0,0,0): particle stays at start.
        let integrator = PathlineIntegrator::new(|_p, _t| [0.0, 0.0, 0.0], 0.1);
        let path = integrator.integrate_rk4([3.0, 4.0, 5.0], 0.0, 1.0);
        for pt in &path {
            assert!((pt[0] - 3.0).abs() < 1e-14);
            assert!((pt[1] - 4.0).abs() < 1e-14);
            assert!((pt[2] - 5.0).abs() < 1e-14);
        }
    }

    // ── StreamlineSeeder ──────────────────────────────────────────────────────

    #[test]
    fn test_streamline_seeder_count() {
        let seeds = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let seeder = StreamlineSeeder::new(seeds, 1.0, 0.1, |_p| [1.0, 0.0, 0.0]);
        let lines = seeder.compute_streamlines();
        assert_eq!(lines.len(), 3, "one streamline per seed");
    }

    #[test]
    fn test_streamline_seeder_uniform_field_straight() {
        let seeds = vec![[0.0, 0.0, 0.0]];
        let seeder = StreamlineSeeder::new(seeds, 1.0, 0.1, |_p| [1.0, 0.0, 0.0]);
        let lines = seeder.compute_streamlines();
        let path = &lines[0];
        // Should move purely in x.
        for pt in path {
            assert!(pt[1].abs() < 1e-10);
            assert!(pt[2].abs() < 1e-10);
        }
    }

    #[test]
    fn test_streamline_seeder_empty_seeds() {
        let seeder = StreamlineSeeder::new(vec![], 1.0, 0.1, |_p| [1.0, 0.0, 0.0]);
        let lines = seeder.compute_streamlines();
        assert!(lines.is_empty());
    }

    #[test]
    fn test_streamline_seeder_integration_length_controls_path_length() {
        let seeds = vec![[0.0, 0.0, 0.0]];
        let s_short = StreamlineSeeder::new(seeds.clone(), 0.5, 0.1, |_p| [1.0, 0.0, 0.0]);
        let s_long = StreamlineSeeder::new(seeds, 2.0, 0.1, |_p| [1.0, 0.0, 0.0]);
        let short_line = &s_short.compute_streamlines()[0];
        let long_line = &s_long.compute_streamlines()[0];
        assert!(
            long_line.len() > short_line.len(),
            "longer integration should produce more points"
        );
    }

    // ── RibbonTrail ───────────────────────────────────────────────────────────

    #[test]
    fn test_ribbon_trail_empty() {
        let ribbon = RibbonTrail::new(0.5);
        assert!(ribbon.is_empty());
        assert_eq!(ribbon.compute_ribbon_geometry().len(), 0);
    }

    #[test]
    fn test_ribbon_trail_single_point() {
        let mut ribbon = RibbonTrail::new(0.5);
        ribbon.add_point([1.0, 2.0, 3.0]);
        let geom = ribbon.compute_ribbon_geometry();
        assert_eq!(geom.len(), 1);
        // Both edges coincide at the single point.
        assert_eq!(geom[0].0, geom[0].1);
    }

    #[test]
    fn test_ribbon_trail_two_points_produce_two_edge_pairs() {
        let mut ribbon = RibbonTrail::new(1.0);
        ribbon.add_point([0.0, 0.0, 0.0]);
        ribbon.add_point([1.0, 0.0, 0.0]);
        let geom = ribbon.compute_ribbon_geometry();
        assert_eq!(geom.len(), 2);
    }

    #[test]
    fn test_ribbon_trail_width_is_correct() {
        // Motion along +x, ribbon up = +y → edges should be offset in ±z.
        let mut ribbon = RibbonTrail::new(1.0);
        ribbon.up_hint = [0.0, 1.0, 0.0];
        ribbon.add_point([0.0, 0.0, 0.0]);
        ribbon.add_point([1.0, 0.0, 0.0]);
        let geom = ribbon.compute_ribbon_geometry();
        let (left, right) = geom[0];
        let dist = vec3_dist(left, right);
        // Expected distance = 2 * half_width = 2.0
        assert!(
            (dist - 2.0).abs() < 1e-8,
            "ribbon width should be 2*half_width=2.0, got {dist}"
        );
    }

    #[test]
    fn test_ribbon_trail_clear() {
        let mut ribbon = RibbonTrail::new(0.5);
        ribbon.add_point([0.0, 0.0, 0.0]);
        ribbon.clear();
        assert!(ribbon.is_empty());
    }

    #[test]
    fn test_ribbon_trail_edges_symmetric_around_center() {
        // For a straight path along +x with up = +y, edges should be symmetric
        // in z around the center line (y=0 plane).
        let mut ribbon = RibbonTrail::new(0.5);
        ribbon.up_hint = [0.0, 1.0, 0.0];
        for i in 0..4 {
            ribbon.add_point([i as f64, 0.0, 0.0]);
        }
        let geom = ribbon.compute_ribbon_geometry();
        for (left, right) in &geom {
            let mid_z = (left[2] + right[2]) / 2.0;
            assert!(
                mid_z.abs() < 1e-8,
                "ribbon edges should be symmetric about center: mid_z={mid_z}"
            );
        }
    }

    // ── Vector helpers ─────────────────────────────────────────────────────────

    #[test]
    fn test_vec3_add() {
        let r = vec3_add([1.0, 2.0, 3.0], [4.0, 5.0, 6.0]);
        assert_eq!(r, [5.0, 7.0, 9.0]);
    }

    #[test]
    fn test_vec3_sub() {
        let r = vec3_sub([5.0, 7.0, 9.0], [1.0, 2.0, 3.0]);
        assert_eq!(r, [4.0, 5.0, 6.0]);
    }

    #[test]
    fn test_vec3_normalize_unit_vector() {
        let n = vec3_normalize([3.0, 0.0, 4.0]);
        let len = vec3_length(n);
        assert!((len - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_vec3_cross_orthogonal() {
        let a = [1.0, 0.0, 0.0];
        let b = [0.0, 1.0, 0.0];
        let c = vec3_cross(a, b);
        assert!((c[2] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_vec3_dot() {
        let d = vec3_dot([1.0, 2.0, 3.0], [4.0, 5.0, 6.0]);
        assert!((d - 32.0).abs() < 1e-10);
    }
}
