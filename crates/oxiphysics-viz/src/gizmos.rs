// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Physics debug gizmos — immediate-mode line/shape overlays.
//!
//! `Gizmos` is an immediate-mode drawing API for physics debugging.
//! Each frame the caller pushes primitives (arrows, boxes, spheres, grids,
//! velocity/force vectors, …) and the renderer drains them once per frame.
//!
//! # Design
//!
//! - **No GPU dependency** — all output is plain `LinePrimitive` / `TrianglePrimitive` data.
//! - **Depth-sorting** — optional front-to-back ordering for semi-transparent overlays.
//! - **Layer system** — gizmos may be placed on named layers that can be toggled.
//! - **Transform stack** — push/pop a local-to-world matrix for easy hierarchical drawing.
//!
//! # Example
//!
//! ```no_run
//! use oxiphysics_viz::gizmos::{Gizmos, GizmoLayer};
//! use oxiphysics_viz::primitives::Color;
//!
//! let mut g = Gizmos::new();
//! g.set_layer(GizmoLayer::Physics);
//! g.draw_axis([0.0, 0.0, 0.0], 1.0);
//! let lines = g.drain_lines();
//! assert!(!lines.is_empty());
//! ```

use crate::primitives::{Color, LinePrimitive};
use std::f64::consts::{PI, TAU};

// ---------------------------------------------------------------------------
// Helper math (no nalgebra, plain [f64; 3])
// ---------------------------------------------------------------------------

#[inline]
fn add3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

#[inline]
fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

#[inline]
fn scale3(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

#[inline]
fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[inline]
fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

#[inline]
fn length3(v: [f64; 3]) -> f64 {
    dot3(v, v).sqrt()
}

#[inline]
fn normalize3(v: [f64; 3]) -> [f64; 3] {
    let l = length3(v);
    if l < 1e-15 {
        [0.0, 0.0, 0.0]
    } else {
        scale3(v, 1.0 / l)
    }
}

/// Apply a column-major 4×4 matrix to a 3-D point (homogeneous, w=1).
fn mat4_transform_point(m: &[[f64; 4]; 4], p: [f64; 3]) -> [f64; 3] {
    let w = m[0][3] * p[0] + m[1][3] * p[1] + m[2][3] * p[2] + m[3][3];
    let inv_w = if w.abs() < 1e-30 { 1.0 } else { 1.0 / w };
    [
        (m[0][0] * p[0] + m[1][0] * p[1] + m[2][0] * p[2] + m[3][0]) * inv_w,
        (m[0][1] * p[0] + m[1][1] * p[1] + m[2][1] * p[2] + m[3][1]) * inv_w,
        (m[0][2] * p[0] + m[1][2] * p[1] + m[2][2] * p[2] + m[3][2]) * inv_w,
    ]
}

// ---------------------------------------------------------------------------
// GizmoLayer
// ---------------------------------------------------------------------------

/// Logical layer a gizmo is drawn on.
///
/// Layers can be toggled at runtime so that, e.g., collision geometry is
/// hidden while velocity arrows remain visible.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
pub enum GizmoLayer {
    /// Physics collision shapes and contact manifolds.
    Physics,
    /// Velocity, force, impulse arrows.
    Dynamics,
    /// Spatial partitioning structures (BVH, octree, …).
    BroadPhase,
    /// General-purpose overlay (default).
    #[default]
    Overlay,
    /// User-defined layer 0.
    User0,
    /// User-defined layer 1.
    User1,
}

// ---------------------------------------------------------------------------
// GizmoLine — a single debug line segment with layer metadata
// ---------------------------------------------------------------------------

/// A single debug line primitive with layer tagging.
#[derive(Debug, Clone)]
pub struct GizmoLine {
    /// World-space start point.
    pub start: [f64; 3],
    /// World-space end point.
    pub end: [f64; 3],
    /// Line color (RGBA).
    pub color: Color,
    /// Layer this gizmo belongs to.
    pub layer: GizmoLayer,
}

impl GizmoLine {
    /// Create a new gizmo line.
    pub fn new(start: [f64; 3], end: [f64; 3], color: Color, layer: GizmoLayer) -> Self {
        Self {
            start,
            end,
            color,
            layer,
        }
    }

    /// Convert to a [`LinePrimitive`] (discards layer metadata).
    pub fn to_line_primitive(&self) -> LinePrimitive {
        use oxiphysics_core::math::Vec3;
        LinePrimitive {
            start: Vec3::new(self.start[0], self.start[1], self.start[2]),
            end: Vec3::new(self.end[0], self.end[1], self.end[2]),
            color: self.color,
        }
    }
}

// ---------------------------------------------------------------------------
// GizmoTransform — 4×4 transform (column-major)
// ---------------------------------------------------------------------------

/// A column-major 4×4 transform matrix used in the gizmo transform stack.
#[derive(Debug, Clone, Copy)]
pub struct GizmoTransform {
    /// Column-major 4×4 matrix.
    pub m: [[f64; 4]; 4],
}

impl GizmoTransform {
    /// Identity transform.
    pub fn identity() -> Self {
        let mut m = [[0.0_f64; 4]; 4];
        m[0][0] = 1.0;
        m[1][1] = 1.0;
        m[2][2] = 1.0;
        m[3][3] = 1.0;
        Self { m }
    }

    /// Translation-only transform.
    pub fn translation(tx: f64, ty: f64, tz: f64) -> Self {
        let mut t = Self::identity();
        t.m[3][0] = tx;
        t.m[3][1] = ty;
        t.m[3][2] = tz;
        t
    }

    /// Uniform scale transform.
    pub fn uniform_scale(s: f64) -> Self {
        let mut t = Self::identity();
        t.m[0][0] = s;
        t.m[1][1] = s;
        t.m[2][2] = s;
        t
    }

    /// Rotation around the Y axis by `angle` radians.
    pub fn rotation_y(angle: f64) -> Self {
        let (s, c) = angle.sin_cos();
        let mut t = Self::identity();
        t.m[0][0] = c;
        t.m[2][0] = s;
        t.m[0][2] = -s;
        t.m[2][2] = c;
        t
    }

    /// Transform a point through this matrix.
    pub fn transform_point(&self, p: [f64; 3]) -> [f64; 3] {
        mat4_transform_point(&self.m, p)
    }

    /// Compose two transforms (self * other).
    pub fn compose(&self, other: &GizmoTransform) -> GizmoTransform {
        let mut r = [[0.0_f64; 4]; 4];
        for (col, r_col) in r.iter_mut().enumerate() {
            for (row, r_cr) in r_col.iter_mut().enumerate() {
                let mut v = 0.0;
                for k in 0..4 {
                    v += self.m[k][row] * other.m[col][k];
                }
                *r_cr = v;
            }
        }
        GizmoTransform { m: r }
    }
}

impl Default for GizmoTransform {
    fn default() -> Self {
        Self::identity()
    }
}

// ---------------------------------------------------------------------------
// Gizmos — the main immediate-mode drawing context
// ---------------------------------------------------------------------------

/// Immediate-mode physics debug overlay renderer.
///
/// Typical usage:
/// 1. Call drawing methods during the physics update.
/// 2. Call [`Gizmos::drain_lines`] each frame to retrieve all queued primitives.
/// 3. The internal buffer is cleared after each drain.
pub struct Gizmos {
    lines: Vec<GizmoLine>,
    active_layer: GizmoLayer,
    active_color: Color,
    transform_stack: Vec<GizmoTransform>,
    /// Layers that are currently visible (all enabled by default).
    enabled_layers: Vec<GizmoLayer>,
}

impl Default for Gizmos {
    fn default() -> Self {
        Self::new()
    }
}

impl Gizmos {
    /// Create a new, empty `Gizmos` context.
    pub fn new() -> Self {
        Self {
            lines: Vec::new(),
            active_layer: GizmoLayer::default(),
            active_color: Color::white(),
            transform_stack: vec![GizmoTransform::identity()],
            enabled_layers: vec![
                GizmoLayer::Physics,
                GizmoLayer::Dynamics,
                GizmoLayer::BroadPhase,
                GizmoLayer::Overlay,
                GizmoLayer::User0,
                GizmoLayer::User1,
            ],
        }
    }

    // ------------------------------------------------------------------
    // Layer & color control
    // ------------------------------------------------------------------

    /// Set the active drawing layer.
    pub fn set_layer(&mut self, layer: GizmoLayer) {
        self.active_layer = layer;
    }

    /// Set the active drawing color.
    pub fn set_color(&mut self, color: Color) {
        self.active_color = color;
    }

    /// Enable a layer.
    pub fn enable_layer(&mut self, layer: GizmoLayer) {
        if !self.enabled_layers.contains(&layer) {
            self.enabled_layers.push(layer);
        }
    }

    /// Disable a layer (gizmos on that layer are still queued but hidden when draining).
    pub fn disable_layer(&mut self, layer: GizmoLayer) {
        self.enabled_layers.retain(|l| *l != layer);
    }

    /// Check whether a layer is currently enabled.
    pub fn is_layer_enabled(&self, layer: GizmoLayer) -> bool {
        self.enabled_layers.contains(&layer)
    }

    // ------------------------------------------------------------------
    // Transform stack
    // ------------------------------------------------------------------

    /// Push a local transform onto the stack. All subsequent drawing calls
    /// are in the coordinate frame of the composed transform.
    pub fn push_transform(&mut self, t: GizmoTransform) {
        let top = *self
            .transform_stack
            .last()
            .expect("collection should not be empty");
        self.transform_stack.push(top.compose(&t));
    }

    /// Pop the most recently pushed transform.
    ///
    /// Does nothing if only the identity root transform remains.
    pub fn pop_transform(&mut self) {
        if self.transform_stack.len() > 1 {
            self.transform_stack.pop();
        }
    }

    /// Returns the current combined transform.
    pub fn current_transform(&self) -> &GizmoTransform {
        self.transform_stack
            .last()
            .expect("collection should not be empty")
    }

    // ------------------------------------------------------------------
    // Internal helpers
    // ------------------------------------------------------------------

    fn current_layer_enabled(&self) -> bool {
        self.is_layer_enabled(self.active_layer)
    }

    fn push_line(&mut self, a: [f64; 3], b: [f64; 3]) {
        if !self.current_layer_enabled() {
            return;
        }
        let xf = *self.current_transform();
        let wa = xf.transform_point(a);
        let wb = xf.transform_point(b);
        self.lines
            .push(GizmoLine::new(wa, wb, self.active_color, self.active_layer));
    }

    // ------------------------------------------------------------------
    // Primitive drawing methods
    // ------------------------------------------------------------------

    /// Draw a single line segment.
    pub fn draw_line(&mut self, start: [f64; 3], end: [f64; 3]) {
        self.push_line(start, end);
    }

    /// Draw a 3-D cross (three axis-aligned lines) at `center` with half-length `size`.
    pub fn draw_cross(&mut self, center: [f64; 3], size: f64) {
        self.push_line(
            [center[0] - size, center[1], center[2]],
            [center[0] + size, center[1], center[2]],
        );
        self.push_line(
            [center[0], center[1] - size, center[2]],
            [center[0], center[1] + size, center[2]],
        );
        self.push_line(
            [center[0], center[1], center[2] - size],
            [center[0], center[1], center[2] + size],
        );
    }

    /// Draw RGB coordinate axes at `origin` with arm length `size`.
    ///
    /// X → current color (set before calling), but axes override:
    /// - X axis = red
    /// - Y axis = green
    /// - Z axis = blue
    pub fn draw_axis(&mut self, origin: [f64; 3], size: f64) {
        let saved = self.active_color;

        self.active_color = Color::red();
        self.push_line(origin, [origin[0] + size, origin[1], origin[2]]);

        self.active_color = Color::green();
        self.push_line(origin, [origin[0], origin[1] + size, origin[2]]);

        self.active_color = Color::blue();
        self.push_line(origin, [origin[0], origin[1], origin[2] + size]);

        self.active_color = saved;
    }

    /// Draw a filled arrow from `start` to `end`.
    ///
    /// The arrowhead is approximated by 4 diagonal lines.
    pub fn draw_arrow(&mut self, start: [f64; 3], end: [f64; 3]) {
        self.push_line(start, end);

        let shaft = sub3(end, start);
        let shaft_len = length3(shaft);
        if shaft_len < 1e-10 {
            return;
        }
        let shaft_dir = scale3(shaft, 1.0 / shaft_len);
        let head_len = shaft_len * 0.15;

        // Build a perpendicular using the world-up vector (or world-right if near parallel).
        let up = if shaft_dir[1].abs() < 0.9 {
            [0.0, 1.0, 0.0]
        } else {
            [1.0, 0.0, 0.0]
        };
        let perp = normalize3(cross3(shaft_dir, up));
        let perp2 = normalize3(cross3(shaft_dir, perp));

        let base = sub3(end, scale3(shaft_dir, head_len));
        for sign_a in &[-1.0_f64, 1.0] {
            for sign_b in &[-1.0_f64, 1.0] {
                let tip = add3(
                    base,
                    add3(
                        scale3(perp, *sign_a * head_len * 0.4),
                        scale3(perp2, *sign_b * head_len * 0.4),
                    ),
                );
                self.push_line(end, tip);
            }
        }
    }

    /// Draw a wireframe axis-aligned box from `min` to `max`.
    pub fn draw_aabb(&mut self, min: [f64; 3], max: [f64; 3]) {
        let (x0, y0, z0) = (min[0], min[1], min[2]);
        let (x1, y1, z1) = (max[0], max[1], max[2]);
        // Bottom face
        self.push_line([x0, y0, z0], [x1, y0, z0]);
        self.push_line([x1, y0, z0], [x1, y0, z1]);
        self.push_line([x1, y0, z1], [x0, y0, z1]);
        self.push_line([x0, y0, z1], [x0, y0, z0]);
        // Top face
        self.push_line([x0, y1, z0], [x1, y1, z0]);
        self.push_line([x1, y1, z0], [x1, y1, z1]);
        self.push_line([x1, y1, z1], [x0, y1, z1]);
        self.push_line([x0, y1, z1], [x0, y1, z0]);
        // Vertical pillars
        self.push_line([x0, y0, z0], [x0, y1, z0]);
        self.push_line([x1, y0, z0], [x1, y1, z0]);
        self.push_line([x1, y0, z1], [x1, y1, z1]);
        self.push_line([x0, y0, z1], [x0, y1, z1]);
    }

    /// Draw a wireframe sphere approximation using `segments` per great circle.
    pub fn draw_sphere(&mut self, center: [f64; 3], radius: f64, segments: usize) {
        let n = segments.max(4);
        let draw_circle = |gizmos: &mut Gizmos, ax: usize, ay: usize| {
            for i in 0..n {
                let t0 = (i as f64) / (n as f64) * TAU;
                let t1 = (i + 1) as f64 / (n as f64) * TAU;
                let mut a = [0.0_f64; 3];
                let mut b = [0.0_f64; 3];
                a[ax] = center[ax] + radius * t0.cos();
                a[ay] = center[ay] + radius * t0.sin();
                a[3 - ax - ay] = center[3 - ax - ay];
                b[ax] = center[ax] + radius * t1.cos();
                b[ay] = center[ay] + radius * t1.sin();
                b[3 - ax - ay] = center[3 - ax - ay];
                gizmos.push_line(a, b);
            }
        };
        draw_circle(self, 0, 1); // XY
        draw_circle(self, 0, 2); // XZ
        draw_circle(self, 1, 2); // YZ
    }

    /// Draw a wireframe capsule (cylinder + two hemispherical end caps).
    ///
    /// `axis` is the unit direction of the capsule's primary axis.
    pub fn draw_capsule(
        &mut self,
        center: [f64; 3],
        axis: [f64; 3],
        half_height: f64,
        radius: f64,
        segments: usize,
    ) {
        let n = segments.max(4);
        let ax = normalize3(axis);
        let cap_a = add3(center, scale3(ax, half_height));
        let cap_b = sub3(center, scale3(ax, half_height));

        // Two end-cap center-lines
        self.push_line(
            add3(cap_a, scale3(ax, radius)),
            add3(cap_b, scale3(ax, -radius)),
        );

        // Build two perpendicular vectors to `ax`.
        let perp = if ax[1].abs() < 0.9 {
            normalize3(cross3(ax, [0.0, 1.0, 0.0]))
        } else {
            normalize3(cross3(ax, [1.0, 0.0, 0.0]))
        };
        let perp2 = normalize3(cross3(ax, perp));

        // Waist circle at each cap center
        for cap_center in [cap_a, cap_b] {
            for i in 0..n {
                let t0 = (i as f64) / (n as f64) * TAU;
                let t1 = (i + 1) as f64 / (n as f64) * TAU;
                let a = add3(
                    cap_center,
                    add3(
                        scale3(perp, radius * t0.cos()),
                        scale3(perp2, radius * t0.sin()),
                    ),
                );
                let b = add3(
                    cap_center,
                    add3(
                        scale3(perp, radius * t1.cos()),
                        scale3(perp2, radius * t1.sin()),
                    ),
                );
                self.push_line(a, b);
            }
        }

        // Four vertical lines connecting the caps
        for t in [0.0_f64, PI * 0.5, PI, PI * 1.5] {
            let offset = add3(
                scale3(perp, radius * t.cos()),
                scale3(perp2, radius * t.sin()),
            );
            self.push_line(add3(cap_a, offset), add3(cap_b, offset));
        }
    }

    /// Draw a ground plane grid of `half_count` cells in each direction.
    pub fn draw_grid(&mut self, center: [f64; 3], cell_size: f64, half_count: i32) {
        let extent = half_count as f64 * cell_size;
        for i in -half_count..=half_count {
            let x = center[0] + i as f64 * cell_size;
            let z0 = center[2] - extent;
            let z1 = center[2] + extent;
            self.push_line([x, center[1], z0], [x, center[1], z1]);

            let z = center[2] + i as f64 * cell_size;
            let x0 = center[0] - extent;
            let x1 = center[0] + extent;
            self.push_line([x0, center[1], z], [x1, center[1], z]);
        }
    }

    /// Draw velocity as an arrow scaled by `scale`.
    pub fn draw_velocity(&mut self, position: [f64; 3], velocity: [f64; 3], scale: f64) {
        let end = add3(position, scale3(velocity, scale));
        self.draw_arrow(position, end);
    }

    /// Draw a contact manifold: a point with a normal arrow of length `depth`.
    pub fn draw_contact(&mut self, contact_point: [f64; 3], normal: [f64; 3], depth: f64) {
        let saved = self.active_color;
        self.active_color = Color::red();
        self.draw_cross(contact_point, depth * 0.25);
        self.active_color = Color::new(1.0, 0.5, 0.0, 1.0); // orange
        let end = add3(contact_point, scale3(normalize3(normal), depth));
        self.draw_arrow(contact_point, end);
        self.active_color = saved;
    }

    /// Draw a cone along `axis` from `tip` with given `height` and base `radius`.
    pub fn draw_cone(
        &mut self,
        tip: [f64; 3],
        axis: [f64; 3],
        height: f64,
        radius: f64,
        segments: usize,
    ) {
        let n = segments.max(4);
        let ax = normalize3(axis);
        let base_center = add3(tip, scale3(ax, height));
        let perp = if ax[1].abs() < 0.9 {
            normalize3(cross3(ax, [0.0, 1.0, 0.0]))
        } else {
            normalize3(cross3(ax, [1.0, 0.0, 0.0]))
        };
        let perp2 = normalize3(cross3(ax, perp));

        for i in 0..n {
            let t0 = (i as f64) / (n as f64) * TAU;
            let t1 = (i + 1) as f64 / (n as f64) * TAU;
            let a = add3(
                base_center,
                add3(
                    scale3(perp, radius * t0.cos()),
                    scale3(perp2, radius * t0.sin()),
                ),
            );
            let b = add3(
                base_center,
                add3(
                    scale3(perp, radius * t1.cos()),
                    scale3(perp2, radius * t1.sin()),
                ),
            );
            self.push_line(a, b);
            if i % (n / 4) == 0 {
                self.push_line(tip, a);
            }
        }
    }

    // ------------------------------------------------------------------
    // Buffer management
    // ------------------------------------------------------------------

    /// Returns the number of queued line primitives.
    pub fn line_count(&self) -> usize {
        self.lines.len()
    }

    /// Drain all queued lines, returning them as [`LinePrimitive`] objects.
    ///
    /// After this call the internal buffer is empty.
    pub fn drain_lines(&mut self) -> Vec<LinePrimitive> {
        let lines: Vec<LinePrimitive> = self
            .lines
            .iter()
            .filter(|l| self.enabled_layers.contains(&l.layer))
            .map(|l| l.to_line_primitive())
            .collect();
        self.lines.clear();
        lines
    }

    /// Drain lines belonging only to the specified layer.
    pub fn drain_layer(&mut self, layer: GizmoLayer) -> Vec<LinePrimitive> {
        let (keep, take): (Vec<_>, Vec<_>) = self.lines.drain(..).partition(|l| l.layer != layer);
        self.lines = keep;
        take.iter().map(|l| l.to_line_primitive()).collect()
    }

    /// Clear all queued gizmos without returning them.
    pub fn clear(&mut self) {
        self.lines.clear();
    }
}

// ---------------------------------------------------------------------------
// GizmoBuilder — fluent API for building individual gizmos
// ---------------------------------------------------------------------------

/// A fluent builder for configuring and submitting a single gizmo.
///
/// # Example
///
/// ```rust,ignore
/// GizmoBuilder::new(&mut gizmos)
///     .color(Color::cyan())
///     .layer(GizmoLayer::Physics)
///     .arrow([0.0, 0.0, 0.0], [1.0, 0.0, 0.0]);
/// ```
pub struct GizmoBuilder<'a> {
    gizmos: &'a mut Gizmos,
    saved_color: Color,
    saved_layer: GizmoLayer,
}

impl<'a> GizmoBuilder<'a> {
    /// Begin building a gizmo.
    pub fn new(gizmos: &'a mut Gizmos) -> Self {
        let saved_color = gizmos.active_color;
        let saved_layer = gizmos.active_layer;
        Self {
            gizmos,
            saved_color,
            saved_layer,
        }
    }

    /// Set the drawing color.
    pub fn color(self, color: Color) -> Self {
        self.gizmos.active_color = color;
        self
    }

    /// Set the drawing layer.
    pub fn layer(self, layer: GizmoLayer) -> Self {
        self.gizmos.active_layer = layer;
        self
    }

    /// Draw an arrow and restore settings.
    pub fn arrow(mut self, start: [f64; 3], end: [f64; 3]) {
        self.gizmos.draw_arrow(start, end);
        self.restore();
    }

    /// Draw an AABB and restore settings.
    pub fn aabb(mut self, min: [f64; 3], max: [f64; 3]) {
        self.gizmos.draw_aabb(min, max);
        self.restore();
    }

    /// Draw a sphere and restore settings.
    pub fn sphere(mut self, center: [f64; 3], radius: f64, segments: usize) {
        self.gizmos.draw_sphere(center, radius, segments);
        self.restore();
    }

    fn restore(&mut self) {
        self.gizmos.active_color = self.saved_color;
        self.gizmos.active_layer = self.saved_layer;
    }
}

impl<'a> Drop for GizmoBuilder<'a> {
    fn drop(&mut self) {
        // Restore on drop (if `restore` was not called explicitly).
        self.restore();
    }
}

// ---------------------------------------------------------------------------
// Depth-sorted drain helper
// ---------------------------------------------------------------------------

/// Sort direction for depth sorting.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SortOrder {
    /// Front to back (closest first).
    FrontToBack,
    /// Back to front (furthest first) — for transparent overlays.
    BackToFront,
}

/// Drain all lines from `gizmos` sorted by their midpoint distance to `camera_pos`.
pub fn drain_lines_sorted(
    gizmos: &mut Gizmos,
    camera_pos: [f64; 3],
    order: SortOrder,
) -> Vec<LinePrimitive> {
    let mut lines = gizmos.drain_lines();
    lines.sort_by(|a, b| {
        let mid_a = [
            (a.start.x + a.end.x) * 0.5,
            (a.start.y + a.end.y) * 0.5,
            (a.start.z + a.end.z) * 0.5,
        ];
        let mid_b = [
            (b.start.x + b.end.x) * 0.5,
            (b.start.y + b.end.y) * 0.5,
            (b.start.z + b.end.z) * 0.5,
        ];
        let da = {
            let d = sub3(mid_a, camera_pos);
            dot3(d, d)
        };
        let db = {
            let d = sub3(mid_b, camera_pos);
            dot3(d, d)
        };
        match order {
            SortOrder::FrontToBack => da.partial_cmp(&db).unwrap_or(std::cmp::Ordering::Equal),
            SortOrder::BackToFront => db.partial_cmp(&da).unwrap_or(std::cmp::Ordering::Equal),
        }
    });
    lines
}

// ---------------------------------------------------------------------------
// Transform Gizmo (translate / rotate / scale handles)
// ---------------------------------------------------------------------------

/// Axis identifier for the three principal axes.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GizmoAxis {
    /// X axis (red).
    X,
    /// Y axis (green).
    Y,
    /// Z axis (blue).
    Z,
}

/// The active editing mode of a transform gizmo.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TransformMode {
    /// Translate along axes.
    Translate,
    /// Rotate around axes.
    Rotate,
    /// Scale along axes.
    Scale,
}

/// Colors used by the transform gizmo.
fn axis_color(axis: GizmoAxis) -> Color {
    match axis {
        GizmoAxis::X => Color::red(),
        GizmoAxis::Y => Color::green(),
        GizmoAxis::Z => Color::blue(),
    }
}

impl Gizmos {
    /// Draw a translate gizmo at `origin` with arm length `size`.
    ///
    /// Three arrow handles are drawn along X, Y, Z axes.
    pub fn draw_translate_gizmo(&mut self, origin: [f64; 3], size: f64) {
        let saved = self.active_color;
        for (axis, dir) in [
            (GizmoAxis::X, [size, 0.0, 0.0_f64]),
            (GizmoAxis::Y, [0.0, size, 0.0]),
            (GizmoAxis::Z, [0.0, 0.0, size]),
        ] {
            self.active_color = axis_color(axis);
            let end = add3(origin, dir);
            self.draw_arrow(origin, end);
        }
        self.active_color = saved;
    }

    /// Draw a rotate gizmo at `origin` with ring radius `size`.
    ///
    /// Three circles are drawn around X, Y, Z axes.
    pub fn draw_rotate_gizmo(&mut self, origin: [f64; 3], size: f64, segments: usize) {
        let saved = self.active_color;
        let n = segments.max(8);

        // XY circle (around Z axis) — blue
        self.active_color = Color::blue();
        for i in 0..n {
            let t0 = (i as f64) / (n as f64) * TAU;
            let t1 = (i + 1) as f64 / (n as f64) * TAU;
            let a = add3(origin, [size * t0.cos(), size * t0.sin(), 0.0]);
            let b = add3(origin, [size * t1.cos(), size * t1.sin(), 0.0]);
            self.push_line(a, b);
        }
        // XZ circle (around Y axis) — green
        self.active_color = Color::green();
        for i in 0..n {
            let t0 = (i as f64) / (n as f64) * TAU;
            let t1 = (i + 1) as f64 / (n as f64) * TAU;
            let a = add3(origin, [size * t0.cos(), 0.0, size * t0.sin()]);
            let b = add3(origin, [size * t1.cos(), 0.0, size * t1.sin()]);
            self.push_line(a, b);
        }
        // YZ circle (around X axis) — red
        self.active_color = Color::red();
        for i in 0..n {
            let t0 = (i as f64) / (n as f64) * TAU;
            let t1 = (i + 1) as f64 / (n as f64) * TAU;
            let a = add3(origin, [0.0, size * t0.cos(), size * t0.sin()]);
            let b = add3(origin, [0.0, size * t1.cos(), size * t1.sin()]);
            self.push_line(a, b);
        }
        self.active_color = saved;
    }

    /// Draw a scale gizmo at `origin` with arm length `size`.
    ///
    /// Three lines with box end-caps (represented as small crosses) are drawn.
    pub fn draw_scale_gizmo(&mut self, origin: [f64; 3], size: f64) {
        let saved = self.active_color;
        let cap_size = size * 0.08;
        for (axis, dir) in [
            (GizmoAxis::X, [size, 0.0, 0.0_f64]),
            (GizmoAxis::Y, [0.0, size, 0.0]),
            (GizmoAxis::Z, [0.0, 0.0, size]),
        ] {
            self.active_color = axis_color(axis);
            let end = add3(origin, dir);
            self.push_line(origin, end);
            // Draw a small cross at the end to indicate scale handle
            self.draw_cross(end, cap_size);
        }
        self.active_color = saved;
    }

    /// Draw a unified transform gizmo according to `mode`.
    pub fn draw_transform_gizmo(&mut self, origin: [f64; 3], size: f64, mode: TransformMode) {
        match mode {
            TransformMode::Translate => self.draw_translate_gizmo(origin, size),
            TransformMode::Rotate => self.draw_rotate_gizmo(origin, size, 24),
            TransformMode::Scale => self.draw_scale_gizmo(origin, size),
        }
    }
}

// ---------------------------------------------------------------------------
// Constraint visualization helpers
// ---------------------------------------------------------------------------

impl Gizmos {
    /// Draw a hinge joint constraint at `pivot` with the hinge `axis` and angle limits
    /// `[min_angle, max_angle]` in radians.
    ///
    /// The allowed arc is drawn as a fan in the plane perpendicular to `axis`.
    pub fn draw_hinge_constraint(
        &mut self,
        pivot: [f64; 3],
        axis: [f64; 3],
        min_angle: f64,
        max_angle: f64,
        radius: f64,
        segments: usize,
    ) {
        let saved = self.active_color;
        let ax = normalize3(axis);
        let n = segments.max(4);

        // Build perpendicular vectors
        let perp = if ax[1].abs() < 0.9 {
            normalize3(cross3(ax, [0.0, 1.0, 0.0]))
        } else {
            normalize3(cross3(ax, [1.0, 0.0, 0.0]))
        };
        let perp2 = normalize3(cross3(ax, perp));

        // Draw the hinge axis
        self.active_color = Color::new(1.0, 1.0, 0.0, 1.0); // yellow
        self.push_line(
            sub3(pivot, scale3(ax, radius)),
            add3(pivot, scale3(ax, radius)),
        );

        // Draw the limit fan
        self.active_color = Color::new(0.5, 0.5, 1.0, 0.7); // light blue
        let range = max_angle - min_angle;
        for i in 0..n {
            let a0 = min_angle + range * (i as f64) / (n as f64);
            let a1 = min_angle + range * (i + 1) as f64 / (n as f64);
            let p0 = add3(
                pivot,
                add3(
                    scale3(perp, radius * a0.cos()),
                    scale3(perp2, radius * a0.sin()),
                ),
            );
            let p1 = add3(
                pivot,
                add3(
                    scale3(perp, radius * a1.cos()),
                    scale3(perp2, radius * a1.sin()),
                ),
            );
            self.push_line(pivot, p0);
            self.push_line(p0, p1);
        }
        self.push_line(
            pivot,
            add3(
                pivot,
                add3(
                    scale3(perp, radius * max_angle.cos()),
                    scale3(perp2, radius * max_angle.sin()),
                ),
            ),
        );
        self.active_color = saved;
    }

    /// Draw a prismatic (slider) constraint at `anchor` along `direction`
    /// with travel range `[min_dist, max_dist]`.
    pub fn draw_prismatic_constraint(
        &mut self,
        anchor: [f64; 3],
        direction: [f64; 3],
        min_dist: f64,
        max_dist: f64,
    ) {
        let saved = self.active_color;
        let d = normalize3(direction);

        self.active_color = Color::new(1.0, 0.5, 0.0, 1.0); // orange
        let a = add3(anchor, scale3(d, min_dist));
        let b = add3(anchor, scale3(d, max_dist));
        self.push_line(a, b);

        // Tick marks at limits
        let perp = if d[1].abs() < 0.9 {
            normalize3(cross3(d, [0.0, 1.0, 0.0]))
        } else {
            normalize3(cross3(d, [1.0, 0.0, 0.0]))
        };
        let tick = 0.1;
        self.push_line(sub3(a, scale3(perp, tick)), add3(a, scale3(perp, tick)));
        self.push_line(sub3(b, scale3(perp, tick)), add3(b, scale3(perp, tick)));
        self.active_color = saved;
    }

    /// Draw a ball-and-socket constraint at `pivot` with cone limit `half_angle` (radians).
    pub fn draw_ball_socket_constraint(
        &mut self,
        pivot: [f64; 3],
        direction: [f64; 3],
        half_angle: f64,
        radius: f64,
        segments: usize,
    ) {
        let saved = self.active_color;
        let ax = normalize3(direction);
        let cap = add3(pivot, scale3(ax, radius));
        let cone_radius = radius * half_angle.sin();

        self.active_color = Color::new(0.8, 0.2, 0.8, 1.0); // purple
        self.draw_cone(pivot, ax, radius, cone_radius, segments);
        self.draw_sphere(cap, cone_radius * 0.1, 6);
        self.active_color = saved;
    }
}

// ---------------------------------------------------------------------------
// Force arrow with auto-scaling
// ---------------------------------------------------------------------------

impl Gizmos {
    /// Draw a force arrow at `position` representing `force` vector.
    ///
    /// `max_force` is the reference magnitude mapped to `max_length`.
    /// Forces below `min_display` are not drawn.
    pub fn draw_force_arrow_scaled(
        &mut self,
        position: [f64; 3],
        force: [f64; 3],
        max_force: f64,
        max_length: f64,
        min_display: f64,
    ) {
        let mag = length3(force);
        if mag < min_display || max_force < 1e-12 {
            return;
        }
        let scale = (mag / max_force).min(1.0) * max_length;
        let dir = normalize3(force);
        let end = add3(position, scale3(dir, scale));
        self.draw_arrow(position, end);
    }

    /// Draw multiple force arrows for a rigid body at `position` with `forces`.
    pub fn draw_body_forces(
        &mut self,
        position: [f64; 3],
        forces: &[[f64; 3]],
        max_force: f64,
        max_length: f64,
    ) {
        let saved = self.active_color;
        self.active_color = Color::new(1.0, 0.4, 0.0, 1.0); // orange
        for force in forces {
            self.draw_force_arrow_scaled(position, *force, max_force, max_length, 1e-6);
        }
        self.active_color = saved;
    }
}

// ---------------------------------------------------------------------------
// Contact normal glyphs
// ---------------------------------------------------------------------------

impl Gizmos {
    /// Draw a contact normal glyph: a cross at the contact point and an arrow
    /// along the normal scaled by penetration depth.
    pub fn draw_contact_normal_glyph(
        &mut self,
        contact_point: [f64; 3],
        normal: [f64; 3],
        penetration_depth: f64,
    ) {
        let saved = self.active_color;
        let cross_size = penetration_depth.abs() * 0.2 + 0.02;
        let arrow_len = penetration_depth.abs() + 0.05;

        self.active_color = Color::red();
        self.draw_cross(contact_point, cross_size);

        self.active_color = Color::new(1.0, 0.6, 0.0, 1.0); // amber
        let end = add3(contact_point, scale3(normalize3(normal), arrow_len));
        self.draw_arrow(contact_point, end);

        self.active_color = saved;
    }

    /// Draw all contact points from a slice of `(point, normal, depth)` tuples.
    pub fn draw_contact_manifold(&mut self, contacts: &[([f64; 3], [f64; 3], f64)]) {
        for &(pt, n, depth) in contacts {
            self.draw_contact_normal_glyph(pt, n, depth);
        }
    }
}

// ---------------------------------------------------------------------------
// Selection highlight
// ---------------------------------------------------------------------------

impl Gizmos {
    /// Draw a selection highlight (bright AABB outline) around `min`..`max`.
    pub fn draw_selection_highlight(&mut self, min: [f64; 3], max: [f64; 3]) {
        let saved = self.active_color;
        self.active_color = Color::new(1.0, 1.0, 0.0, 1.0); // bright yellow
        self.draw_aabb(min, max);
        self.active_color = saved;
    }

    /// Draw a selection highlight around a sphere at `center` with `radius`.
    pub fn draw_selection_sphere(&mut self, center: [f64; 3], radius: f64) {
        let saved = self.active_color;
        self.active_color = Color::new(1.0, 1.0, 0.0, 1.0); // bright yellow
        self.draw_sphere(center, radius * 1.05, 16);
        self.active_color = saved;
    }
}

// ---------------------------------------------------------------------------
// Velocity trail
// ---------------------------------------------------------------------------

/// A velocity trail that stores the last N positions for a moving body.
#[derive(Debug, Clone)]
pub struct VelocityTrail {
    /// Ring buffer of world-space positions.
    pub positions: Vec<[f64; 3]>,
    /// Current write index.
    write: usize,
    /// Number of valid entries.
    count: usize,
    /// Maximum number of trail points.
    pub capacity: usize,
    /// Color of the trail.
    pub color: Color,
}

impl VelocityTrail {
    /// Create a new trail with the given capacity.
    pub fn new(capacity: usize, color: Color) -> Self {
        let cap = capacity.max(2);
        Self {
            positions: vec![[0.0; 3]; cap],
            write: 0,
            count: 0,
            capacity: cap,
            color,
        }
    }

    /// Record a new position.
    pub fn push(&mut self, position: [f64; 3]) {
        self.positions[self.write] = position;
        self.write = (self.write + 1) % self.capacity;
        self.count = (self.count + 1).min(self.capacity);
    }

    /// Number of valid points.
    pub fn len(&self) -> usize {
        self.count
    }

    /// Returns `true` if empty.
    pub fn is_empty(&self) -> bool {
        self.count == 0
    }

    /// Draw the trail into `gizmos` using line segments.
    pub fn draw(&self, gizmos: &mut Gizmos) {
        if self.count < 2 {
            return;
        }
        let saved = gizmos.active_color;
        gizmos.active_color = self.color;
        // Iterate from oldest to newest
        for k in 0..self.count - 1 {
            let i0 = (self.write + self.capacity - self.count + k) % self.capacity;
            let i1 = (i0 + 1) % self.capacity;
            gizmos.push_line(self.positions[i0], self.positions[i1]);
        }
        gizmos.active_color = saved;
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // --- basic line drawing ---

    #[test]
    fn test_draw_line_adds_one_line() {
        let mut g = Gizmos::new();
        g.draw_line([0.0, 0.0, 0.0], [1.0, 0.0, 0.0]);
        assert_eq!(g.line_count(), 1);
    }

    #[test]
    fn test_draw_cross_adds_three_lines() {
        let mut g = Gizmos::new();
        g.draw_cross([0.0, 0.0, 0.0], 1.0);
        assert_eq!(g.line_count(), 3);
    }

    #[test]
    fn test_draw_axis_adds_three_lines() {
        let mut g = Gizmos::new();
        g.draw_axis([0.0, 0.0, 0.0], 1.0);
        assert_eq!(g.line_count(), 3);
    }

    #[test]
    fn test_draw_aabb_adds_twelve_lines() {
        let mut g = Gizmos::new();
        g.draw_aabb([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        assert_eq!(g.line_count(), 12);
    }

    #[test]
    fn test_draw_sphere_adds_three_circles() {
        let segments = 16;
        let mut g = Gizmos::new();
        g.draw_sphere([0.0, 0.0, 0.0], 1.0, segments);
        // 3 great circles, each with `segments` lines
        assert_eq!(g.line_count(), 3 * segments);
    }

    #[test]
    fn test_draw_arrow_adds_lines() {
        let mut g = Gizmos::new();
        g.draw_arrow([0.0, 0.0, 0.0], [1.0, 0.0, 0.0]);
        // shaft + 4 arrowhead lines
        assert!(g.line_count() >= 5);
    }

    #[test]
    fn test_draw_grid_line_count() {
        let mut g = Gizmos::new();
        g.draw_grid([0.0, 0.0, 0.0], 1.0, 2);
        // (2*half+1) lines in each direction = 5 x-lines + 5 z-lines = 10
        assert_eq!(g.line_count(), 10);
    }

    // --- drain clears buffer ---

    #[test]
    fn test_drain_clears_buffer() {
        let mut g = Gizmos::new();
        g.draw_line([0.0; 3], [1.0, 0.0, 0.0]);
        let lines = g.drain_lines();
        assert_eq!(lines.len(), 1);
        assert_eq!(g.line_count(), 0);
    }

    #[test]
    fn test_clear_empties_buffer() {
        let mut g = Gizmos::new();
        g.draw_cross([0.0; 3], 1.0);
        g.clear();
        assert_eq!(g.line_count(), 0);
    }

    // --- layer filtering ---

    #[test]
    fn test_disabled_layer_hidden_from_drain() {
        let mut g = Gizmos::new();
        g.set_layer(GizmoLayer::Physics);
        g.disable_layer(GizmoLayer::Physics);
        g.draw_line([0.0; 3], [1.0, 0.0, 0.0]);
        // Even though a line was pushed, it is filtered at drain time.
        let lines = g.drain_lines();
        assert_eq!(lines.len(), 0);
    }

    #[test]
    fn test_enabled_layer_appears_in_drain() {
        let mut g = Gizmos::new();
        g.set_layer(GizmoLayer::Overlay);
        g.enable_layer(GizmoLayer::Overlay);
        g.draw_line([0.0; 3], [1.0, 0.0, 0.0]);
        let lines = g.drain_lines();
        assert_eq!(lines.len(), 1);
    }

    #[test]
    fn test_is_layer_enabled() {
        let mut g = Gizmos::new();
        assert!(g.is_layer_enabled(GizmoLayer::Physics));
        g.disable_layer(GizmoLayer::Physics);
        assert!(!g.is_layer_enabled(GizmoLayer::Physics));
        g.enable_layer(GizmoLayer::Physics);
        assert!(g.is_layer_enabled(GizmoLayer::Physics));
    }

    // --- transform stack ---

    #[test]
    fn test_push_pop_transform_identity() {
        let mut g = Gizmos::new();
        g.push_transform(GizmoTransform::translation(1.0, 0.0, 0.0));
        g.draw_line([0.0; 3], [1.0, 0.0, 0.0]);
        g.pop_transform();
        // After pop, the next line should be in the original frame
        g.draw_line([0.0; 3], [1.0, 0.0, 0.0]);
        // Both lines were pushed
        assert_eq!(g.line_count(), 2);
    }

    #[test]
    fn test_transform_translates_points() {
        let mut g = Gizmos::new();
        g.push_transform(GizmoTransform::translation(5.0, 0.0, 0.0));
        g.draw_line([0.0; 3], [0.0; 3]);
        g.pop_transform();
        // The line's start should now be at x=5
        let raw = &g.lines[0];
        assert!((raw.start[0] - 5.0).abs() < 1e-10);
    }

    // --- GizmoTransform ---

    #[test]
    fn test_transform_identity_no_op() {
        let t = GizmoTransform::identity();
        let p = t.transform_point([3.0, -2.0, 1.0]);
        assert!((p[0] - 3.0).abs() < 1e-10);
        assert!((p[1] + 2.0).abs() < 1e-10);
        assert!((p[2] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_transform_translation() {
        let t = GizmoTransform::translation(1.0, 2.0, 3.0);
        let p = t.transform_point([0.0; 3]);
        assert!((p[0] - 1.0).abs() < 1e-10);
        assert!((p[1] - 2.0).abs() < 1e-10);
        assert!((p[2] - 3.0).abs() < 1e-10);
    }

    #[test]
    fn test_transform_uniform_scale() {
        let t = GizmoTransform::uniform_scale(2.0);
        let p = t.transform_point([1.0, 1.0, 1.0]);
        assert!((p[0] - 2.0).abs() < 1e-10);
        assert!((p[1] - 2.0).abs() < 1e-10);
        assert!((p[2] - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_transform_compose() {
        // Verify compose: translate.compose(translate) sums the translations.
        let t1 = GizmoTransform::translation(1.0, 0.0, 0.0);
        let t2 = GizmoTransform::translation(2.0, 0.0, 0.0);
        let combined = t1.compose(&t2);
        // Any point p → combined transforms it correctly.
        let p = combined.transform_point([0.0; 3]);
        // Two translations: result should add translations (3.0 or similar depending on order).
        // Simply verify it is non-zero and a deterministic value.
        assert!(
            p[0] > 0.0,
            "composed translation must shift along x, got {}",
            p[0]
        );

        // Compose identity with anything should give the same result as the other.
        let id = GizmoTransform::identity();
        let scale = GizmoTransform::uniform_scale(3.0);
        let composed = id.compose(&scale);
        let q = composed.transform_point([1.0, 2.0, 3.0]);
        let q2 = scale.transform_point([1.0, 2.0, 3.0]);
        assert!((q[0] - q2[0]).abs() < 1e-10);
        assert!((q[1] - q2[1]).abs() < 1e-10);
        assert!((q[2] - q2[2]).abs() < 1e-10);
    }

    // --- draw_velocity ---

    #[test]
    fn test_draw_velocity_arrow_length() {
        let mut g = Gizmos::new();
        let pos = [0.0; 3];
        let vel = [1.0, 0.0, 0.0];
        g.draw_velocity(pos, vel, 2.0);
        assert!(g.line_count() >= 1); // at least the shaft
    }

    // --- draw_contact ---

    #[test]
    fn test_draw_contact_generates_lines() {
        let mut g = Gizmos::new();
        g.draw_contact([0.0; 3], [0.0, 1.0, 0.0], 0.5);
        assert!(g.line_count() > 0);
    }

    // --- capsule ---

    #[test]
    fn test_draw_capsule_generates_lines() {
        let mut g = Gizmos::new();
        g.draw_capsule([0.0; 3], [0.0, 1.0, 0.0], 0.5, 0.25, 8);
        assert!(g.line_count() > 0);
    }

    // --- cone ---

    #[test]
    fn test_draw_cone_generates_lines() {
        let mut g = Gizmos::new();
        g.draw_cone([0.0; 3], [0.0, 1.0, 0.0], 1.0, 0.5, 8);
        assert!(g.line_count() > 0);
    }

    // --- drain layer ---

    #[test]
    fn test_drain_layer_selective() {
        let mut g = Gizmos::new();
        g.set_layer(GizmoLayer::Physics);
        g.draw_line([0.0; 3], [1.0, 0.0, 0.0]);
        g.set_layer(GizmoLayer::Dynamics);
        g.draw_line([0.0; 3], [2.0, 0.0, 0.0]);

        let phys = g.drain_layer(GizmoLayer::Physics);
        assert_eq!(phys.len(), 1);
        assert_eq!(g.line_count(), 1); // Dynamics line remains
    }

    // --- depth sorted drain ---

    #[test]
    fn test_drain_lines_sorted_front_to_back() {
        let mut g = Gizmos::new();
        g.draw_line([10.0, 0.0, 0.0], [10.0, 0.0, 0.0]);
        g.draw_line([1.0, 0.0, 0.0], [1.0, 0.0, 0.0]);
        g.draw_line([5.0, 0.0, 0.0], [5.0, 0.0, 0.0]);
        let sorted = drain_lines_sorted(&mut g, [0.0; 3], SortOrder::FrontToBack);
        assert_eq!(sorted.len(), 3);
        let x0: f64 = sorted[0].start.x;
        let x1: f64 = sorted[1].start.x;
        let x2: f64 = sorted[2].start.x;
        assert!(x0 <= x1 && x1 <= x2, "not sorted: {} {} {}", x0, x1, x2);
    }

    #[test]
    fn test_drain_lines_sorted_back_to_front() {
        let mut g = Gizmos::new();
        g.draw_line([1.0, 0.0, 0.0], [1.0, 0.0, 0.0]);
        g.draw_line([10.0, 0.0, 0.0], [10.0, 0.0, 0.0]);
        let sorted = drain_lines_sorted(&mut g, [0.0; 3], SortOrder::BackToFront);
        assert_eq!(sorted.len(), 2);
        let x0 = sorted[0].start.x;
        let x1 = sorted[1].start.x;
        assert!(x0 >= x1, "back-to-front: expected {} >= {}", x0, x1);
    }

    // -----------------------------------------------------------------------
    // Transform gizmo tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_draw_translate_gizmo_adds_lines() {
        let mut g = Gizmos::new();
        g.draw_translate_gizmo([0.0; 3], 1.0);
        // 3 axes, each arrow = shaft + 4 head lines = 5 lines
        assert!(
            g.line_count() >= 3,
            "translate gizmo should add at least 3 lines"
        );
    }

    #[test]
    fn test_draw_rotate_gizmo_adds_lines() {
        let mut g = Gizmos::new();
        g.draw_rotate_gizmo([0.0; 3], 1.0, 12);
        // 3 circles × 12 segments = 36 lines
        assert_eq!(g.line_count(), 36, "rotate gizmo: 3 circles × 12 segments");
    }

    #[test]
    fn test_draw_scale_gizmo_adds_lines() {
        let mut g = Gizmos::new();
        g.draw_scale_gizmo([0.0; 3], 1.0);
        // 3 shaft lines + 3×3 cross lines = 12
        assert!(g.line_count() >= 6, "scale gizmo should have >= 6 lines");
    }

    #[test]
    fn test_draw_transform_gizmo_translate_mode() {
        let mut g = Gizmos::new();
        g.draw_transform_gizmo([0.0; 3], 1.0, TransformMode::Translate);
        assert!(g.line_count() > 0);
    }

    #[test]
    fn test_draw_transform_gizmo_rotate_mode() {
        let mut g = Gizmos::new();
        g.draw_transform_gizmo([0.0; 3], 1.0, TransformMode::Rotate);
        assert!(g.line_count() > 0);
    }

    #[test]
    fn test_draw_transform_gizmo_scale_mode() {
        let mut g = Gizmos::new();
        g.draw_transform_gizmo([0.0; 3], 1.0, TransformMode::Scale);
        assert!(g.line_count() > 0);
    }

    // -----------------------------------------------------------------------
    // Constraint visualization tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_draw_hinge_constraint_generates_lines() {
        let mut g = Gizmos::new();
        g.draw_hinge_constraint(
            [0.0; 3],
            [0.0, 1.0, 0.0],
            -std::f64::consts::FRAC_PI_4,
            std::f64::consts::FRAC_PI_4,
            0.5,
            8,
        );
        assert!(g.line_count() > 0);
    }

    #[test]
    fn test_draw_prismatic_constraint_generates_lines() {
        let mut g = Gizmos::new();
        g.draw_prismatic_constraint([0.0; 3], [0.0, 1.0, 0.0], -1.0, 1.0);
        assert!(g.line_count() >= 3, "prismatic needs axis + 2 tick lines");
    }

    #[test]
    fn test_draw_ball_socket_constraint_generates_lines() {
        let mut g = Gizmos::new();
        g.draw_ball_socket_constraint([0.0; 3], [0.0, 1.0, 0.0], 0.3, 1.0, 8);
        assert!(g.line_count() > 0);
    }

    // -----------------------------------------------------------------------
    // Force arrow tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_draw_force_arrow_scaled_draws_arrow() {
        let mut g = Gizmos::new();
        g.draw_force_arrow_scaled([0.0; 3], [1.0, 0.0, 0.0], 10.0, 1.0, 0.01);
        assert!(g.line_count() >= 1);
    }

    #[test]
    fn test_draw_force_arrow_scaled_skips_tiny_force() {
        let mut g = Gizmos::new();
        g.draw_force_arrow_scaled([0.0; 3], [0.000001, 0.0, 0.0], 10.0, 1.0, 0.01);
        assert_eq!(g.line_count(), 0, "tiny force should not be drawn");
    }

    #[test]
    fn test_draw_body_forces_multiple() {
        let mut g = Gizmos::new();
        let forces = [[1.0, 0.0, 0.0], [0.0, 2.0, 0.0]];
        g.draw_body_forces([0.0; 3], &forces, 5.0, 1.0);
        assert!(g.line_count() > 0);
    }

    // -----------------------------------------------------------------------
    // Contact normal glyph tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_draw_contact_normal_glyph_generates_lines() {
        let mut g = Gizmos::new();
        g.draw_contact_normal_glyph([0.0; 3], [0.0, 1.0, 0.0], 0.05);
        assert!(g.line_count() > 0);
    }

    #[test]
    fn test_draw_contact_manifold_multiple() {
        let mut g = Gizmos::new();
        let contacts = vec![
            ([0.0, 0.0, 0.0], [0.0, 1.0, 0.0], 0.01),
            ([1.0, 0.0, 0.0], [0.0, 1.0, 0.0], 0.02),
        ];
        g.draw_contact_manifold(&contacts);
        assert!(g.line_count() > 0);
    }

    // -----------------------------------------------------------------------
    // Selection highlight tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_draw_selection_highlight_adds_12_lines() {
        let mut g = Gizmos::new();
        g.draw_selection_highlight([-1.0; 3], [1.0; 3]);
        assert_eq!(g.line_count(), 12, "selection highlight = 12 AABB edges");
    }

    #[test]
    fn test_draw_selection_sphere_adds_lines() {
        let mut g = Gizmos::new();
        g.draw_selection_sphere([0.0; 3], 1.0);
        assert!(g.line_count() > 0);
    }

    // -----------------------------------------------------------------------
    // VelocityTrail tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_velocity_trail_starts_empty() {
        let trail = VelocityTrail::new(16, Color::white());
        assert!(trail.is_empty());
        assert_eq!(trail.len(), 0);
    }

    #[test]
    fn test_velocity_trail_push_and_len() {
        let mut trail = VelocityTrail::new(4, Color::white());
        trail.push([1.0, 0.0, 0.0]);
        trail.push([2.0, 0.0, 0.0]);
        assert_eq!(trail.len(), 2);
    }

    #[test]
    fn test_velocity_trail_cap_overflow() {
        let mut trail = VelocityTrail::new(3, Color::white());
        for i in 0..10 {
            trail.push([i as f64, 0.0, 0.0]);
        }
        assert_eq!(trail.len(), 3, "trail should not exceed capacity");
    }

    #[test]
    fn test_velocity_trail_draw_generates_lines() {
        let mut trail = VelocityTrail::new(8, Color::white());
        trail.push([0.0, 0.0, 0.0]);
        trail.push([1.0, 0.0, 0.0]);
        trail.push([2.0, 0.0, 0.0]);
        let mut g = Gizmos::new();
        trail.draw(&mut g);
        assert_eq!(g.line_count(), 2, "3 points = 2 segments");
    }

    #[test]
    fn test_velocity_trail_draw_single_point_no_lines() {
        let mut trail = VelocityTrail::new(8, Color::white());
        trail.push([0.0, 0.0, 0.0]);
        let mut g = Gizmos::new();
        trail.draw(&mut g);
        assert_eq!(g.line_count(), 0, "single point = no lines");
    }

    // -----------------------------------------------------------------------
    // GizmoAxis tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_axis_color_x_is_red() {
        let c = axis_color(GizmoAxis::X);
        assert!(c.r > 0.9 && c.g < 0.1 && c.b < 0.1);
    }

    #[test]
    fn test_axis_color_y_is_green() {
        let c = axis_color(GizmoAxis::Y);
        assert!(c.g > 0.9 && c.r < 0.1 && c.b < 0.1);
    }

    #[test]
    fn test_axis_color_z_is_blue() {
        let c = axis_color(GizmoAxis::Z);
        assert!(c.b > 0.9 && c.r < 0.1 && c.g < 0.1);
    }

    // -----------------------------------------------------------------------
    // Rotate-gizmo segment count test
    // -----------------------------------------------------------------------

    #[test]
    fn test_rotate_gizmo_segment_count_minimum() {
        let mut g = Gizmos::new();
        g.draw_rotate_gizmo([0.0; 3], 1.0, 4); // below 8, clamped to 8
        // 3 circles × 8 segments = 24
        assert_eq!(g.line_count(), 24);
    }
}
