#![allow(clippy::needless_range_loop, clippy::ptr_arg)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Path and trajectory visualization.
//!
//! Provides data structures and algorithms for rendering particle and object
//! paths, including tube/ribbon geometry generation, arc-length
//! parameterization, Frenet frames, curvature, path smoothing, and
//! aerospace flight-path utilities.  No GPU or windowing dependency.

#![allow(dead_code)]
#![allow(clippy::too_many_arguments)]

use std::collections::HashMap;

// ── Vector helpers ────────────────────────────────────────────────────────────

fn vec3_add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

fn vec3_sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn vec3_scale(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

fn vec3_dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn vec3_cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn vec3_length(a: [f64; 3]) -> f64 {
    (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt()
}

fn vec3_normalize(a: [f64; 3]) -> [f64; 3] {
    let len = vec3_length(a);
    if len < 1e-14 {
        return [0.0, 0.0, 1.0];
    }
    [a[0] / len, a[1] / len, a[2] / len]
}

fn vec3_lerp(a: [f64; 3], b: [f64; 3], t: f64) -> [f64; 3] {
    vec3_add(vec3_scale(a, 1.0 - t), vec3_scale(b, t))
}

// ── PathStyle ─────────────────────────────────────────────────────────────────

/// Visual style for rendering a path.
#[derive(Debug, Clone, PartialEq)]
pub enum PathStyle {
    /// Continuous solid line.
    Solid,
    /// Dashed line pattern.
    Dashed,
    /// Dotted line pattern.
    Dotted,
    /// Color gradient along the path.
    Gradient,
    /// Flat ribbon following the path.
    Ribbon,
    /// Volumetric tube around the path.
    Tube,
    /// Arrows indicating direction of travel.
    Arrow,
}

// ── PathColor ─────────────────────────────────────────────────────────────────

/// Color scheme for a rendered path.
#[derive(Debug, Clone, PartialEq)]
pub enum PathColor {
    /// Uniform RGBA color (components in \[0, 1\]).
    Uniform([f64; 4]),
    /// Color mapped from speed magnitude.
    VelocityMapped,
    /// Color mapped from normalized simulation time.
    TimeMapped,
    /// Color mapped from a per-point scalar field.
    ScalarMapped(Vec<f64>),
    /// Full rainbow (hue sweeps 0 → 360° along the path).
    Rainbow,
}

// ── PathSegment ───────────────────────────────────────────────────────────────

/// A single rendered line segment of a path.
#[derive(Debug, Clone, PartialEq)]
pub struct PathSegment {
    /// World-space start position.
    pub start: [f64; 3],
    /// World-space end position.
    pub end: [f64; 3],
    /// Rendered line/ribbon width.
    pub width: f64,
    /// RGBA color at this segment.
    pub color: [f64; 4],
    /// Whether to render an arrow head at `end`.
    pub arrow_at_end: bool,
}

impl PathSegment {
    /// Construct a new `PathSegment`.
    pub fn new(
        start: [f64; 3],
        end: [f64; 3],
        width: f64,
        color: [f64; 4],
        arrow_at_end: bool,
    ) -> Self {
        Self {
            start,
            end,
            width,
            color,
            arrow_at_end,
        }
    }

    /// Length of this segment.
    pub fn length(&self) -> f64 {
        vec3_length(vec3_sub(self.end, self.start))
    }
}

// ── PathData ──────────────────────────────────────────────────────────────────

/// Raw path data: a polyline with optional attributes.
#[derive(Debug, Clone)]
pub struct PathData {
    /// Ordered sequence of world-space points.
    pub points: Vec<[f64; 3]>,
    /// Simulation timestamps, one per point.
    pub timestamps: Vec<f64>,
    /// Optional per-point velocity vectors.
    pub velocities: Option<Vec<[f64; 3]>>,
    /// Optional per-point scalar field values (e.g. temperature).
    pub scalars: Option<Vec<f64>>,
    /// Visual style.
    pub style: PathStyle,
    /// Color scheme.
    pub color: PathColor,
}

impl PathData {
    /// Construct a new `PathData` from points and timestamps.
    pub fn new(
        points: Vec<[f64; 3]>,
        timestamps: Vec<f64>,
        style: PathStyle,
        color: PathColor,
    ) -> Self {
        Self {
            points,
            timestamps,
            velocities: None,
            scalars: None,
            style,
            color,
        }
    }

    /// Total arc length of the path.
    pub fn arc_length(&self) -> f64 {
        self.points
            .windows(2)
            .map(|w| vec3_length(vec3_sub(w[1], w[0])))
            .sum()
    }

    /// Number of points in the path.
    pub fn len(&self) -> usize {
        self.points.len()
    }

    /// Return `true` if the path has no points.
    pub fn is_empty(&self) -> bool {
        self.points.is_empty()
    }
}

// ── Arc-length parameterization ───────────────────────────────────────────────

/// Compute cumulative arc-length values for `points`.
///
/// Returns a `Vec` of the same length as `points` where the first entry is
/// `0.0` and each subsequent entry is the cumulative distance from the start.
pub fn arc_length_parameterize(points: &[[f64; 3]]) -> Vec<f64> {
    let mut s = Vec::with_capacity(points.len());
    s.push(0.0);
    for w in points.windows(2) {
        let prev = *s.last().expect("collection should not be empty");
        s.push(prev + vec3_length(vec3_sub(w[1], w[0])));
    }
    s
}

// ── Path smoothing ────────────────────────────────────────────────────────────

/// Smooth `points` with a Gaussian-weighted window of half-width `window`.
///
/// Points near the ends are treated with a smaller effective window to avoid
/// boundary artefacts.  `window = 0` returns a clone of the input unchanged.
pub fn smooth_path(points: &[[f64; 3]], window: usize) -> Vec<[f64; 3]> {
    if window == 0 || points.len() < 3 {
        return points.to_vec();
    }
    let n = points.len();
    let sigma = window as f64;
    let mut out = Vec::with_capacity(n);
    for i in 0..n {
        let lo = i.saturating_sub(window);
        let hi = (i + window + 1).min(n);
        let mut sum = [0.0f64; 3];
        let mut wt_total = 0.0f64;
        for j in lo..hi {
            let dx = (i as f64) - (j as f64);
            let wt = (-0.5 * (dx / sigma).powi(2)).exp();
            sum = vec3_add(sum, vec3_scale(points[j], wt));
            wt_total += wt;
        }
        out.push(vec3_scale(sum, 1.0 / wt_total));
    }
    out
}

// ── Frenet frame ──────────────────────────────────────────────────────────────

/// Compute the Frenet–Serret frame at a point.
///
/// Returns `[T, N, B]` (tangent, normal, binormal) as a 3×3 row-major matrix.
/// - `p` is unused (reserved for future curvature-at-point use).
/// - `tangent` — unit tangent vector (need not be pre-normalized).
/// - `up` — reference up vector used to compute N when the natural normal is degenerate.
pub fn frenet_frame(p: &[f64; 3], tangent: &[f64; 3], up: &[f64; 3]) -> [[f64; 3]; 3] {
    let _ = p; // future use
    let t = vec3_normalize(*tangent);
    // binormal = T × up, then normal = B × T
    let mut b = vec3_normalize(vec3_cross(t, *up));
    if vec3_length(b) < 1e-10 {
        // tangent is parallel to up — choose a different reference
        let alt_up = if t[0].abs() < 0.9 {
            [1.0, 0.0, 0.0]
        } else {
            [0.0, 1.0, 0.0]
        };
        b = vec3_normalize(vec3_cross(t, alt_up));
    }
    let normal = vec3_normalize(vec3_cross(b, t));
    [t, normal, b]
}

// ── Path curvature ────────────────────────────────────────────────────────────

/// Estimate the curvature of a polyline at index `i` using finite differences.
///
/// Returns `0.0` at the endpoints or when the segment length is near zero.
pub fn path_curvature(points: &[[f64; 3]], i: usize) -> f64 {
    if i == 0 || i + 1 >= points.len() {
        return 0.0;
    }
    let pm = points[i - 1];
    let p = points[i];
    let pp = points[i + 1];
    let d1 = vec3_sub(p, pm);
    let d2 = vec3_sub(pp, p);
    let l1 = vec3_length(d1);
    let l2 = vec3_length(d2);
    if l1 < 1e-14 || l2 < 1e-14 {
        return 0.0;
    }
    let t1 = vec3_normalize(d1);
    let t2 = vec3_normalize(d2);
    let dt = vec3_sub(t2, t1);
    let ds = 0.5 * (l1 + l2);
    vec3_length(dt) / ds
}

// ── PathColormap ──────────────────────────────────────────────────────────────

/// Maps scalar values to RGBA colors using named colormaps.
#[derive(Debug, Clone)]
pub struct PathColormap;

impl PathColormap {
    /// Map `t ∈ [0, 1]` using the cool-warm diverging colormap.
    ///
    /// Blue (cool) at `t=0`, white at `t=0.5`, red (warm) at `t=1`.
    pub fn cool_warm(t: f64) -> [f64; 4] {
        let t = t.clamp(0.0, 1.0);
        let r = t;
        let g = 1.0 - (2.0 * (t - 0.5)).abs();
        let b = 1.0 - t;
        [r, g, b, 1.0]
    }

    /// Map `t ∈ [0, 1]` using the classic jet colormap.
    ///
    /// Blue → cyan → green → yellow → red.
    pub fn jet(t: f64) -> [f64; 4] {
        let t = t.clamp(0.0, 1.0);
        let r = (1.5 - (4.0 * t - 3.0).abs()).clamp(0.0, 1.0);
        let g = (1.5 - (4.0 * t - 2.0).abs()).clamp(0.0, 1.0);
        let b = (1.5 - (4.0 * t - 1.0).abs()).clamp(0.0, 1.0);
        [r, g, b, 1.0]
    }

    /// Map `t ∈ [0, 1]` to a full rainbow hue (HSV with S=V=1).
    pub fn rainbow(t: f64) -> [f64; 4] {
        let t = t.clamp(0.0, 1.0);
        let h = t * 360.0;
        let (r, g, b) = hsv_to_rgb(h, 1.0, 1.0);
        [r, g, b, 1.0]
    }
}

fn hsv_to_rgb(h: f64, s: f64, v: f64) -> (f64, f64, f64) {
    let h = h.rem_euclid(360.0);
    let hi = (h / 60.0).floor() as i32;
    let f = h / 60.0 - hi as f64;
    let p = v * (1.0 - s);
    let q = v * (1.0 - f * s);
    let t = v * (1.0 - (1.0 - f) * s);
    match hi {
        0 => (v, t, p),
        1 => (q, v, p),
        2 => (p, v, t),
        3 => (p, q, v),
        4 => (t, p, v),
        _ => (v, p, q),
    }
}

// ── PathIntersection ──────────────────────────────────────────────────────────

/// Represents an intersection point between two paths.
#[derive(Debug, Clone, PartialEq)]
pub struct PathIntersection {
    /// Approximate 3D world position of the crossing.
    pub position: [f64; 3],
    /// Index in path A just before the crossing.
    pub idx_a: usize,
    /// Index in path B just before the crossing.
    pub idx_b: usize,
    /// Distance between the two paths at the closest approach.
    pub distance: f64,
}

/// Find approximate intersections (close approaches) between two 3-D polylines.
///
/// Returns all pairs of segments whose midpoints are within `threshold` of
/// each other.  For exact 3-D intersections `threshold` can be small (e.g. `0.1`).
pub fn find_path_intersections(
    path_a: &[[f64; 3]],
    path_b: &[[f64; 3]],
    threshold: f64,
) -> Vec<PathIntersection> {
    let mut result = Vec::new();
    for (i, w_a) in path_a.windows(2).enumerate() {
        let mid_a = vec3_scale(vec3_add(w_a[0], w_a[1]), 0.5);
        for (j, w_b) in path_b.windows(2).enumerate() {
            let mid_b = vec3_scale(vec3_add(w_b[0], w_b[1]), 0.5);
            let dist = vec3_length(vec3_sub(mid_a, mid_b));
            if dist <= threshold {
                result.push(PathIntersection {
                    position: vec3_scale(vec3_add(mid_a, mid_b), 0.5),
                    idx_a: i,
                    idx_b: j,
                    distance: dist,
                });
            }
        }
    }
    result
}

// ── PathRenderer ─────────────────────────────────────────────────────────────

/// Converts `PathData` into renderable geometry.
#[derive(Debug, Clone)]
pub struct PathRenderer {
    /// Default visual style.
    pub style: PathStyle,
    /// Default color scheme.
    pub color: PathColor,
}

impl PathRenderer {
    /// Create a new `PathRenderer` with the given style and color.
    pub fn new(style: PathStyle, color: PathColor) -> Self {
        Self { style, color }
    }

    /// Generate one `PathSegment` per consecutive pair of points in `path`.
    pub fn generate_segments(path: &PathData) -> Vec<PathSegment> {
        let n = path.points.len();
        if n < 2 {
            return Vec::new();
        }
        let arc = arc_length_parameterize(&path.points);
        let total = *arc.last().unwrap_or(&1.0);
        let mut segs = Vec::with_capacity(n - 1);
        for i in 0..n - 1 {
            let t = if total > 1e-14 {
                arc[i] / total
            } else {
                i as f64 / (n - 1) as f64
            };
            let color = Self::resolve_color(&path.color, t, path, i);
            segs.push(PathSegment {
                start: path.points[i],
                end: path.points[i + 1],
                width: 1.0,
                color,
                arrow_at_end: false,
            });
        }
        segs
    }

    fn resolve_color(scheme: &PathColor, t: f64, path: &PathData, idx: usize) -> [f64; 4] {
        match scheme {
            PathColor::Uniform(c) => *c,
            PathColor::TimeMapped => PathColormap::cool_warm(t),
            PathColor::Rainbow => PathColormap::rainbow(t),
            PathColor::VelocityMapped => {
                if let Some(vels) = &path.velocities {
                    if let Some(&v) = vels.get(idx) {
                        let speed = vec3_length(v);
                        PathColormap::jet((speed / 10.0).clamp(0.0, 1.0))
                    } else {
                        [1.0, 1.0, 1.0, 1.0]
                    }
                } else {
                    [1.0, 1.0, 1.0, 1.0]
                }
            }
            PathColor::ScalarMapped(scalars) => {
                let v = scalars.get(idx).copied().unwrap_or(0.0);
                PathColormap::jet(v.clamp(0.0, 1.0))
            }
        }
    }

    /// Generate a tube mesh around `path`.
    ///
    /// Returns `(vertices, triangles)` where each triangle is a triple of vertex indices.
    pub fn generate_tube(
        path: &PathData,
        radius: f64,
        n_sides: usize,
    ) -> (Vec<[f64; 3]>, Vec<[usize; 3]>) {
        let n_sides = n_sides.max(3);
        let n = path.points.len();
        if n < 2 {
            return (Vec::new(), Vec::new());
        }
        let up = [0.0f64, 1.0, 0.0];
        let mut verts: Vec<[f64; 3]> = Vec::new();
        for i in 0..n {
            let tangent = if i + 1 < n {
                vec3_sub(path.points[i + 1], path.points[i])
            } else {
                vec3_sub(path.points[i], path.points[i - 1])
            };
            let frame = frenet_frame(&path.points[i], &tangent, &up);
            let normal = frame[1]; // N
            let binormal = frame[2]; // B
            for s in 0..n_sides {
                let angle = 2.0 * std::f64::consts::PI * (s as f64) / (n_sides as f64);
                let offset = vec3_add(
                    vec3_scale(normal, radius * angle.cos()),
                    vec3_scale(binormal, radius * angle.sin()),
                );
                verts.push(vec3_add(path.points[i], offset));
            }
        }
        let mut tris: Vec<[usize; 3]> = Vec::new();
        for ring in 0..n - 1 {
            for s in 0..n_sides {
                let s_next = (s + 1) % n_sides;
                let a = ring * n_sides + s;
                let b = ring * n_sides + s_next;
                let c = (ring + 1) * n_sides + s;
                let d = (ring + 1) * n_sides + s_next;
                tris.push([a, b, c]);
                tris.push([b, d, c]);
            }
        }
        (verts, tris)
    }

    /// Generate a ribbon mesh along `path`.
    ///
    /// `normal` is the ribbon's face normal (need not be pre-normalized).
    /// Returns `(vertices, triangles)`.
    pub fn generate_ribbon(
        path: &PathData,
        width: f64,
        normal: [f64; 3],
    ) -> (Vec<[f64; 3]>, Vec<[usize; 3]>) {
        let n = path.points.len();
        if n < 2 {
            return (Vec::new(), Vec::new());
        }
        let half_w = width * 0.5;
        let mut verts: Vec<[f64; 3]> = Vec::new();
        for i in 0..n {
            let tangent = if i + 1 < n {
                vec3_sub(path.points[i + 1], path.points[i])
            } else {
                vec3_sub(path.points[i], path.points[i - 1])
            };
            let t = vec3_normalize(tangent);
            let side = vec3_normalize(vec3_cross(t, normal));
            verts.push(vec3_add(path.points[i], vec3_scale(side, half_w)));
            verts.push(vec3_add(path.points[i], vec3_scale(side, -half_w)));
        }
        let mut tris: Vec<[usize; 3]> = Vec::new();
        for i in 0..n - 1 {
            let a = i * 2;
            let b = i * 2 + 1;
            let c = (i + 1) * 2;
            let d = (i + 1) * 2 + 1;
            tris.push([a, b, c]);
            tris.push([b, d, c]);
        }
        (verts, tris)
    }

    /// Add arrow-head markers (via `arrow_at_end = true`) to segments spaced
    /// approximately every `1.0 / density` of the total path length.
    pub fn add_arrow_heads(segments: &mut Vec<PathSegment>, size: f64) {
        let _ = size;
        if segments.is_empty() {
            return;
        }
        // Mark last segment and every ~10th segment.
        let n = segments.len();
        let step = (n / 5).max(1);
        for i in (0..n).step_by(step) {
            segments[i].arrow_at_end = true;
        }
        segments[n - 1].arrow_at_end = true;
    }
}

// ── TrajectorySet ─────────────────────────────────────────────────────────────

/// A collection of named trajectories (one per particle or object).
#[derive(Debug, Clone, Default)]
pub struct TrajectorySet {
    /// Map from object ID to its path.
    pub trajectories: HashMap<usize, PathData>,
}

impl TrajectorySet {
    /// Create an empty `TrajectorySet`.
    pub fn new() -> Self {
        Self::default()
    }

    /// Add or replace the trajectory for `id`.
    pub fn add_trajectory(&mut self, id: usize, path: PathData) {
        self.trajectories.insert(id, path);
    }

    /// Render all trajectories using `renderer`, returning a flat list of segments.
    pub fn render_all(&self, _renderer: &PathRenderer) -> Vec<PathSegment> {
        let mut all = Vec::new();
        for path in self.trajectories.values() {
            all.extend(PathRenderer::generate_segments(path));
        }
        all
    }

    /// Return a new `TrajectorySet` keeping only the sub-path within `[t_min, t_max]`.
    pub fn filter_by_time(&self, t_min: f64, t_max: f64) -> TrajectorySet {
        let mut out = TrajectorySet::new();
        for (&id, path) in &self.trajectories {
            let filtered: Vec<usize> = path
                .timestamps
                .iter()
                .enumerate()
                .filter(|&(_, &t)| t >= t_min && t <= t_max)
                .map(|(i, _)| i)
                .collect();
            if filtered.is_empty() {
                continue;
            }
            let points = filtered.iter().map(|&i| path.points[i]).collect();
            let timestamps = filtered.iter().map(|&i| path.timestamps[i]).collect();
            let velocities = path
                .velocities
                .as_ref()
                .map(|v| filtered.iter().map(|&i| v[i]).collect());
            let scalars = path
                .scalars
                .as_ref()
                .map(|s| filtered.iter().map(|&i| s[i]).collect());
            out.add_trajectory(
                id,
                PathData {
                    points,
                    timestamps,
                    velocities,
                    scalars,
                    style: path.style.clone(),
                    color: path.color.clone(),
                },
            );
        }
        out
    }

    /// Total arc length across all trajectories.
    pub fn total_length(&self) -> f64 {
        self.trajectories.values().map(|p| p.arc_length()).sum()
    }

    /// Speed statistics (min, mean, max) computed from velocities across all trajectories.
    ///
    /// Returns `(0, 0, 0)` when no velocity data is available.
    pub fn speed_statistics(&self) -> (f64, f64, f64) {
        let speeds: Vec<f64> = self
            .trajectories
            .values()
            .filter_map(|p| p.velocities.as_ref())
            .flat_map(|vels| vels.iter().map(|&v| vec3_length(v)))
            .collect();
        if speeds.is_empty() {
            return (0.0, 0.0, 0.0);
        }
        let min_s = speeds.iter().cloned().fold(f64::INFINITY, f64::min);
        let max_s = speeds.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let mean_s = speeds.iter().sum::<f64>() / speeds.len() as f64;
        (min_s, mean_s, max_s)
    }
}

// ── FlightPath ────────────────────────────────────────────────────────────────

/// Aerospace flight-path visualization with altitude-based coloring.
#[derive(Debug, Clone)]
pub struct FlightPath {
    /// Underlying path data (x = longitude-like, y = altitude, z = latitude-like).
    pub path: PathData,
    /// Reference sea-level altitude (in simulation units).
    pub sea_level: f64,
    /// Maximum altitude for color scaling.
    pub max_altitude: f64,
}

impl FlightPath {
    /// Create a new `FlightPath`.
    pub fn new(path: PathData, sea_level: f64, max_altitude: f64) -> Self {
        Self {
            path,
            sea_level,
            max_altitude,
        }
    }

    /// Generate segments where color encodes altitude (y component).
    pub fn altitude_colored_segments(&self) -> Vec<PathSegment> {
        let n = self.path.points.len();
        if n < 2 {
            return Vec::new();
        }
        let range = (self.max_altitude - self.sea_level).max(1.0);
        (0..n - 1)
            .map(|i| {
                let alt = self.path.points[i][1];
                let t = ((alt - self.sea_level) / range).clamp(0.0, 1.0);
                let color = PathColormap::cool_warm(t);
                PathSegment::new(
                    self.path.points[i],
                    self.path.points[i + 1],
                    1.0,
                    color,
                    false,
                )
            })
            .collect()
    }

    /// Maximum altitude reached along the flight path.
    pub fn peak_altitude(&self) -> f64 {
        self.path
            .points
            .iter()
            .map(|p| p[1])
            .fold(f64::NEG_INFINITY, f64::max)
    }

    /// Minimum altitude (closest to ground) along the flight path.
    pub fn min_altitude(&self) -> f64 {
        self.path
            .points
            .iter()
            .map(|p| p[1])
            .fold(f64::INFINITY, f64::min)
    }
}

// ── ParticleTrails ────────────────────────────────────────────────────────────

/// Fade-out history trails for a set of particles.
///
/// Each particle accumulates a fixed-length rolling window of past positions.
/// Alpha fades from `1.0` (newest) to `0.0` (oldest).
#[derive(Debug, Clone)]
pub struct ParticleTrails {
    /// Map from particle ID to its position history.
    pub trails: HashMap<usize, Vec<[f64; 3]>>,
    /// Maximum number of history points per particle.
    pub max_history: usize,
    /// Base RGBA color (alpha is modulated by age).
    pub base_color: [f64; 4],
}

impl ParticleTrails {
    /// Create a new `ParticleTrails` container.
    pub fn new(max_history: usize, base_color: [f64; 4]) -> Self {
        Self {
            trails: HashMap::new(),
            max_history: max_history.max(1),
            base_color,
        }
    }

    /// Append a new position for particle `id`, pruning oldest entries if needed.
    pub fn update(&mut self, id: usize, position: [f64; 3]) {
        let hist = self.trails.entry(id).or_default();
        hist.push(position);
        if hist.len() > self.max_history {
            hist.remove(0);
        }
    }

    /// Generate faded `PathSegment`s for particle `id`.
    ///
    /// Returns an empty `Vec` if the particle has fewer than two history points.
    pub fn segments_for(&self, id: usize) -> Vec<PathSegment> {
        let hist = match self.trails.get(&id) {
            Some(h) if h.len() >= 2 => h,
            _ => return Vec::new(),
        };
        let n = hist.len();
        (0..n - 1)
            .map(|i| {
                // newest entry is at the back → highest alpha
                let alpha = (i + 1) as f64 / n as f64;
                let color = [
                    self.base_color[0],
                    self.base_color[1],
                    self.base_color[2],
                    self.base_color[3] * alpha,
                ];
                PathSegment::new(hist[i], hist[i + 1], 1.0, color, false)
            })
            .collect()
    }

    /// Total number of tracked particles.
    pub fn n_particles(&self) -> usize {
        self.trails.len()
    }

    /// Clear all history.
    pub fn clear(&mut self) {
        self.trails.clear();
    }
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    fn straight_path(n: usize) -> PathData {
        let points: Vec<[f64; 3]> = (0..n).map(|i| [i as f64, 0.0, 0.0]).collect();
        let timestamps: Vec<f64> = (0..n).map(|i| i as f64 * 0.1).collect();
        PathData::new(
            points,
            timestamps,
            PathStyle::Solid,
            PathColor::Uniform([1.0, 1.0, 1.0, 1.0]),
        )
    }

    // -- arc_length_parameterize --

    #[test]
    fn test_arc_length_straight_line() {
        let pts: Vec<[f64; 3]> = (0..5).map(|i| [i as f64, 0.0, 0.0]).collect();
        let s = arc_length_parameterize(&pts);
        assert_eq!(s.len(), 5);
        assert!((s[0]).abs() < 1e-12);
        assert!((s[4] - 4.0).abs() < 1e-10);
    }

    #[test]
    fn test_arc_length_single_point() {
        let pts = vec![[0.0, 0.0, 0.0]];
        let s = arc_length_parameterize(&pts);
        assert_eq!(s.len(), 1);
        assert!(s[0].abs() < 1e-12);
    }

    #[test]
    fn test_arc_length_monotone() {
        let pts: Vec<[f64; 3]> = (0..10).map(|i| [i as f64, (i as f64).sin(), 0.0]).collect();
        let s = arc_length_parameterize(&pts);
        for w in s.windows(2) {
            assert!(w[1] >= w[0]);
        }
    }

    // -- smooth_path --

    #[test]
    fn test_smooth_path_zero_window_identity() {
        let pts: Vec<[f64; 3]> = vec![[0.0, 0.0, 0.0], [1.0, 1.0, 0.0], [2.0, 0.0, 0.0]];
        let smoothed = smooth_path(&pts, 0);
        assert_eq!(smoothed.len(), pts.len());
        for (a, b) in pts.iter().zip(smoothed.iter()) {
            assert!((a[0] - b[0]).abs() < 1e-10);
        }
    }

    #[test]
    fn test_smooth_path_output_length_preserved() {
        let pts: Vec<[f64; 3]> = (0..10).map(|i| [i as f64, 0.0, 0.0]).collect();
        let smoothed = smooth_path(&pts, 2);
        assert_eq!(smoothed.len(), pts.len());
    }

    #[test]
    fn test_smooth_path_reduces_spike() {
        // Insert a sharp spike at index 5; smoothing should reduce it.
        let mut pts: Vec<[f64; 3]> = (0..11).map(|i| [i as f64, 0.0, 0.0]).collect();
        pts[5][1] = 100.0; // spike
        let smoothed = smooth_path(&pts, 2);
        // After smoothing, index 5 y-value should be significantly less than 100.
        assert!(smoothed[5][1] < 50.0);
    }

    // -- frenet_frame --

    #[test]
    fn test_frenet_frame_tangent_normalized() {
        let p = [0.0, 0.0, 0.0];
        let t = [2.0, 0.0, 0.0];
        let up = [0.0, 1.0, 0.0];
        let frame = frenet_frame(&p, &t, &up);
        let tlen = vec3_length(frame[0]);
        assert!((tlen - 1.0).abs() < 1e-10, "tangent not unit: {tlen}");
    }

    #[test]
    fn test_frenet_frame_orthogonality() {
        let p = [0.0, 0.0, 0.0];
        let t = [1.0, 1.0, 0.0];
        let up = [0.0, 0.0, 1.0];
        let frame = frenet_frame(&p, &t, &up);
        let dot_tn = vec3_dot(frame[0], frame[1]);
        let dot_tb = vec3_dot(frame[0], frame[2]);
        assert!(dot_tn.abs() < 1e-10, "T·N not zero: {dot_tn}");
        assert!(dot_tb.abs() < 1e-10, "T·B not zero: {dot_tb}");
    }

    // -- path_curvature --

    #[test]
    fn test_curvature_straight_line_is_zero() {
        let pts: Vec<[f64; 3]> = (0..5).map(|i| [i as f64, 0.0, 0.0]).collect();
        for i in 1..4 {
            let k = path_curvature(&pts, i);
            assert!(
                k.abs() < 1e-10,
                "curvature of straight line should be 0, got {k}"
            );
        }
    }

    #[test]
    fn test_curvature_endpoint_is_zero() {
        let pts = vec![[0.0, 0.0, 0.0], [1.0, 1.0, 0.0], [2.0, 0.0, 0.0]];
        assert!(path_curvature(&pts, 0).abs() < 1e-10);
        assert!(path_curvature(&pts, 2).abs() < 1e-10);
    }

    #[test]
    fn test_curvature_corner_is_positive() {
        // 90-degree corner: curvature > 0.
        let pts = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0]];
        let k = path_curvature(&pts, 1);
        assert!(k > 0.0, "corner should have positive curvature, got {k}");
    }

    // -- PathColormap --

    #[test]
    fn test_cool_warm_endpoints() {
        let cold = PathColormap::cool_warm(0.0);
        let warm = PathColormap::cool_warm(1.0);
        assert!(cold[2] > cold[0], "t=0 should be more blue than red");
        assert!(warm[0] > warm[2], "t=1 should be more red than blue");
    }

    #[test]
    fn test_jet_midpoint_green() {
        let c = PathColormap::jet(0.5);
        assert!(c[1] > c[0], "jet midpoint g > r");
        assert!(c[1] > c[2], "jet midpoint g > b");
    }

    #[test]
    fn test_rainbow_endpoints_distinct() {
        // t=0 → red (hue 0°), t=0.5 → cyan (hue 180°) — clearly distinct.
        let c0 = PathColormap::rainbow(0.0);
        let c1 = PathColormap::rainbow(0.5);
        let diff: f64 = (0..3).map(|i| (c0[i] - c1[i]).abs()).sum();
        assert!(
            diff > 0.5,
            "rainbow t=0 and t=0.5 should differ, diff={diff}"
        );
    }

    // -- PathSegment --

    #[test]
    fn test_path_segment_length() {
        let seg = PathSegment::new(
            [0.0, 0.0, 0.0],
            [3.0, 4.0, 0.0],
            1.0,
            [1.0, 1.0, 1.0, 1.0],
            false,
        );
        assert!((seg.length() - 5.0).abs() < 1e-10);
    }

    // -- PathData --

    #[test]
    fn test_path_data_arc_length() {
        let d = straight_path(6); // 5 unit segments
        assert!((d.arc_length() - 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_path_data_len() {
        let d = straight_path(8);
        assert_eq!(d.len(), 8);
    }

    #[test]
    fn test_path_data_is_empty() {
        let d = PathData::new(
            vec![],
            vec![],
            PathStyle::Solid,
            PathColor::Uniform([1.0, 1.0, 1.0, 1.0]),
        );
        assert!(d.is_empty());
    }

    // -- PathRenderer --

    #[test]
    fn test_generate_segments_count() {
        let d = straight_path(5);
        let segs = PathRenderer::generate_segments(&d);
        assert_eq!(segs.len(), 4); // n-1 segments for n points
    }

    #[test]
    fn test_generate_segments_empty_path() {
        let d = straight_path(1);
        let segs = PathRenderer::generate_segments(&d);
        assert!(segs.is_empty());
    }

    #[test]
    fn test_generate_tube_vertex_count() {
        let d = straight_path(4);
        let (verts, _) = PathRenderer::generate_tube(&d, 0.1, 6);
        assert_eq!(verts.len(), 4 * 6); // n_points * n_sides
    }

    #[test]
    fn test_generate_tube_triangle_count() {
        let d = straight_path(4);
        let (_, tris) = PathRenderer::generate_tube(&d, 0.1, 4);
        // (n_points-1) * n_sides * 2 triangles
        assert_eq!(tris.len(), 3 * 4 * 2);
    }

    #[test]
    fn test_generate_ribbon_vertex_count() {
        let d = straight_path(5);
        let (verts, _) = PathRenderer::generate_ribbon(&d, 0.2, [0.0, 1.0, 0.0]);
        assert_eq!(verts.len(), 5 * 2); // 2 verts per cross-section
    }

    #[test]
    fn test_add_arrow_heads_last_segment() {
        let d = straight_path(6);
        let mut segs = PathRenderer::generate_segments(&d);
        PathRenderer::add_arrow_heads(&mut segs, 0.1);
        assert!(segs.last().unwrap().arrow_at_end);
    }

    // -- TrajectorySet --

    #[test]
    fn test_trajectory_set_total_length() {
        let mut ts = TrajectorySet::new();
        ts.add_trajectory(0, straight_path(6)); // arc = 5
        ts.add_trajectory(1, straight_path(6)); // arc = 5
        assert!((ts.total_length() - 10.0).abs() < 1e-10);
    }

    #[test]
    fn test_trajectory_set_filter_by_time() {
        let mut ts = TrajectorySet::new();
        ts.add_trajectory(0, straight_path(10)); // timestamps 0.0…0.9
        let filtered = ts.filter_by_time(0.2, 0.5);
        let path = filtered.trajectories.get(&0).unwrap();
        // timestamps in [0.2, 0.5] → indices 2,3,4,5
        assert_eq!(path.points.len(), 4);
    }

    #[test]
    fn test_trajectory_set_speed_statistics_no_vel() {
        let mut ts = TrajectorySet::new();
        ts.add_trajectory(0, straight_path(5));
        let (mn, mean, mx) = ts.speed_statistics();
        assert!((mn + mean + mx).abs() < 1e-10); // all zero
    }

    #[test]
    fn test_trajectory_set_speed_statistics_with_vel() {
        let mut ts = TrajectorySet::new();
        let mut p = straight_path(4);
        p.velocities = Some(vec![
            [1.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [3.0, 0.0, 0.0],
            [4.0, 0.0, 0.0],
        ]);
        ts.add_trajectory(0, p);
        let (mn, _mean, mx) = ts.speed_statistics();
        assert!((mn - 1.0).abs() < 1e-10);
        assert!((mx - 4.0).abs() < 1e-10);
    }

    // -- FlightPath --

    #[test]
    fn test_flight_path_peak_altitude() {
        let pts = vec![[0.0, 100.0, 0.0], [1.0, 500.0, 0.0], [2.0, 200.0, 0.0]];
        let ts = vec![0.0, 1.0, 2.0];
        let pd = PathData::new(
            pts,
            ts,
            PathStyle::Solid,
            PathColor::Uniform([1.0, 1.0, 1.0, 1.0]),
        );
        let fp = FlightPath::new(pd, 0.0, 1000.0);
        assert!((fp.peak_altitude() - 500.0).abs() < 1e-10);
    }

    #[test]
    fn test_flight_path_altitude_segments_count() {
        let pts: Vec<[f64; 3]> = (0..5).map(|i| [i as f64, i as f64 * 100.0, 0.0]).collect();
        let ts: Vec<f64> = (0..5).map(|i| i as f64).collect();
        let pd = PathData::new(
            pts,
            ts,
            PathStyle::Solid,
            PathColor::Uniform([1.0, 1.0, 1.0, 1.0]),
        );
        let fp = FlightPath::new(pd, 0.0, 500.0);
        let segs = fp.altitude_colored_segments();
        assert_eq!(segs.len(), 4);
    }

    // -- ParticleTrails --

    #[test]
    fn test_particle_trails_max_history() {
        let mut trails = ParticleTrails::new(3, [1.0, 1.0, 1.0, 1.0]);
        for i in 0..10 {
            trails.update(0, [i as f64, 0.0, 0.0]);
        }
        assert!(trails.trails[&0].len() <= 3);
    }

    #[test]
    fn test_particle_trails_segments_alpha_fade() {
        let mut trails = ParticleTrails::new(5, [1.0, 0.0, 0.0, 1.0]);
        for i in 0..5 {
            trails.update(0, [i as f64, 0.0, 0.0]);
        }
        let segs = trails.segments_for(0);
        assert!(!segs.is_empty());
        // Earlier segments should have lower alpha than later ones.
        let first_alpha = segs[0].color[3];
        let last_alpha = segs.last().unwrap().color[3];
        assert!(last_alpha >= first_alpha);
    }

    #[test]
    fn test_particle_trails_clear() {
        let mut trails = ParticleTrails::new(5, [1.0, 1.0, 1.0, 1.0]);
        trails.update(0, [0.0, 0.0, 0.0]);
        trails.clear();
        assert_eq!(trails.n_particles(), 0);
    }

    // -- PathIntersection --

    #[test]
    fn test_find_intersections_same_path() {
        let pts: Vec<[f64; 3]> = (0..5).map(|i| [i as f64, 0.0, 0.0]).collect();
        let hits = find_path_intersections(&pts, &pts, 0.01);
        // Diagonal (same segment), every segment matches itself.
        assert!(!hits.is_empty());
    }

    #[test]
    fn test_find_intersections_no_overlap() {
        let a: Vec<[f64; 3]> = (0..3).map(|i| [i as f64, 0.0, 0.0]).collect();
        let b: Vec<[f64; 3]> = (0..3).map(|i| [i as f64, 100.0, 0.0]).collect();
        let hits = find_path_intersections(&a, &b, 1.0);
        assert!(hits.is_empty());
    }
}
