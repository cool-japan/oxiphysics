// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Accumulated debug draw commands for a single frame.

/// RGBA color with components in `[0, 1]`.
#[derive(Debug, Clone, Copy)]
pub struct Color {
    /// Red channel.
    pub r: f32,
    /// Green channel.
    pub g: f32,
    /// Blue channel.
    pub b: f32,
    /// Alpha channel.
    pub a: f32,
}

impl Color {
    /// Opaque red.
    pub const RED: Color = Color {
        r: 1.0,
        g: 0.0,
        b: 0.0,
        a: 1.0,
    };
    /// Opaque green.
    pub const GREEN: Color = Color {
        r: 0.0,
        g: 1.0,
        b: 0.0,
        a: 1.0,
    };
    /// Opaque blue.
    pub const BLUE: Color = Color {
        r: 0.0,
        g: 0.0,
        b: 1.0,
        a: 1.0,
    };
    /// Opaque white.
    pub const WHITE: Color = Color {
        r: 1.0,
        g: 1.0,
        b: 1.0,
        a: 1.0,
    };
    /// Opaque yellow.
    pub const YELLOW: Color = Color {
        r: 1.0,
        g: 1.0,
        b: 0.0,
        a: 1.0,
    };

    /// Construct a color from RGBA components.
    pub fn rgba(r: f32, g: f32, b: f32, a: f32) -> Self {
        Self { r, g, b, a }
    }

    /// Return a copy of this color with a different alpha value.
    pub fn with_alpha(self, a: f32) -> Self {
        Self { a, ..self }
    }
}

/// A single debug draw command.
#[derive(Debug, Clone)]
pub enum DebugCommand {
    /// A line segment.
    Line {
        /// Start position.
        start: [f32; 3],
        /// End position.
        end: [f32; 3],
        /// Render color.
        color: Color,
        /// Line thickness in world units.
        thickness: f32,
    },
    /// A sphere outline.
    Sphere {
        /// Centre position.
        center: [f32; 3],
        /// Sphere radius.
        radius: f32,
        /// Render color.
        color: Color,
    },
    /// An axis-aligned box outline.
    Box {
        /// Centre position.
        center: [f32; 3],
        /// Half-extents on each axis.
        half_extents: [f32; 3],
        /// Render color.
        color: Color,
    },
    /// A directional arrow.
    Arrow {
        /// Arrow origin.
        origin: [f32; 3],
        /// Unit direction vector.
        direction: [f32; 3],
        /// Length of the arrow.
        length: f32,
        /// Render color.
        color: Color,
    },
    /// Screen-space text label.
    Text {
        /// World-space anchor position.
        position: [f32; 3],
        /// Label string.
        text: String,
        /// Render color.
        color: Color,
        /// Font size in points.
        size: f32,
    },
    /// A single point.
    Point {
        /// Position.
        position: [f32; 3],
        /// Render color.
        color: Color,
        /// Point size in pixels.
        size: f32,
    },
}

/// Accumulated debug draw buffer; cleared once per frame.
pub struct DebugOverlay {
    commands: Vec<DebugCommand>,
    /// Maximum number of commands this overlay will accept.
    pub max_commands: usize,
}

impl DebugOverlay {
    /// Create a new overlay with the given command capacity limit.
    pub fn new(max_commands: usize) -> Self {
        Self {
            commands: Vec::new(),
            max_commands,
        }
    }

    /// Remove all pending commands.
    pub fn clear(&mut self) {
        self.commands.clear();
    }

    /// Add a line segment.
    pub fn line(&mut self, start: [f32; 3], end: [f32; 3], color: Color) {
        if self.is_full() {
            return;
        }
        self.commands.push(DebugCommand::Line {
            start,
            end,
            color,
            thickness: 1.0,
        });
    }

    /// Add a sphere outline.
    pub fn sphere(&mut self, center: [f32; 3], radius: f32, color: Color) {
        if self.is_full() {
            return;
        }
        self.commands.push(DebugCommand::Sphere {
            center,
            radius,
            color,
        });
    }

    /// Add the 12 wireframe edges of an axis-aligned bounding box.
    pub fn aabb(&mut self, min: [f32; 3], max: [f32; 3], color: Color) {
        let [x0, y0, z0] = min;
        let [x1, y1, z1] = max;

        // 8 corners of the AABB
        let corners: [[f32; 3]; 8] = [
            [x0, y0, z0],
            [x1, y0, z0],
            [x1, y1, z0],
            [x0, y1, z0],
            [x0, y0, z1],
            [x1, y0, z1],
            [x1, y1, z1],
            [x0, y1, z1],
        ];

        // 12 edges (pairs of corner indices)
        let edges: [(usize, usize); 12] = [
            // bottom face
            (0, 1),
            (1, 2),
            (2, 3),
            (3, 0),
            // top face
            (4, 5),
            (5, 6),
            (6, 7),
            (7, 4),
            // verticals
            (0, 4),
            (1, 5),
            (2, 6),
            (3, 7),
        ];

        for (a, b) in edges {
            self.line(corners[a], corners[b], color);
        }
    }

    /// Add a directional arrow.
    pub fn arrow(&mut self, origin: [f32; 3], dir: [f32; 3], len: f32, color: Color) {
        if self.is_full() {
            return;
        }
        self.commands.push(DebugCommand::Arrow {
            origin,
            direction: dir,
            length: len,
            color,
        });
    }

    /// Add a point marker.
    pub fn point(&mut self, pos: [f32; 3], color: Color) {
        if self.is_full() {
            return;
        }
        self.commands.push(DebugCommand::Point {
            position: pos,
            color,
            size: 4.0,
        });
    }

    /// Add a text label at a world-space position.
    pub fn text(&mut self, pos: [f32; 3], text: impl Into<String>, color: Color) {
        if self.is_full() {
            return;
        }
        self.commands.push(DebugCommand::Text {
            position: pos,
            text: text.into(),
            color,
            size: 12.0,
        });
    }

    /// Borrow the slice of pending commands.
    pub fn commands(&self) -> &[DebugCommand] {
        &self.commands
    }

    /// Number of pending commands.
    pub fn len(&self) -> usize {
        self.commands.len()
    }

    /// Returns `true` if there are no pending commands.
    pub fn is_empty(&self) -> bool {
        self.commands.is_empty()
    }

    /// Returns `true` when the command buffer has reached its capacity limit.
    pub fn is_full(&self) -> bool {
        self.commands.len() >= self.max_commands
    }
}

// ─── ContactDebug ─────────────────────────────────────────────────────────────

/// Visualizes contact points, normals, and penetration depths.
pub struct ContactDebug<'a> {
    overlay: &'a mut DebugOverlay,
    /// Color used for contact normal arrows.
    pub normal_color: Color,
    /// Color used for contact point markers.
    pub point_color: Color,
    /// Scale factor applied to the penetration depth when drawing the normal arrow.
    pub depth_scale: f32,
}

impl<'a> ContactDebug<'a> {
    /// Create a new `ContactDebug` that draws into `overlay`.
    pub fn new(overlay: &'a mut DebugOverlay) -> Self {
        Self {
            overlay,
            normal_color: Color::RED,
            point_color: Color::YELLOW,
            depth_scale: 5.0,
        }
    }

    /// Draw a contact: a point marker at `position` and an arrow along `normal`
    /// scaled by `penetration_depth * depth_scale`.
    pub fn draw_contact(&mut self, position: [f32; 3], normal: [f32; 3], penetration_depth: f32) {
        self.overlay.point(position, self.point_color);
        let len = penetration_depth * self.depth_scale;
        self.overlay.arrow(position, normal, len, self.normal_color);
    }

    /// Draw a batch of contacts.
    ///
    /// Each entry is `(position, normal, penetration_depth)`.
    pub fn draw_contacts(&mut self, contacts: &[([f32; 3], [f32; 3], f32)]) {
        for &(pos, nor, depth) in contacts {
            self.draw_contact(pos, nor, depth);
        }
    }
}

// ─── ForceDebug ───────────────────────────────────────────────────────────────

/// Draws force and torque vectors on rigid bodies.
pub struct ForceDebug<'a> {
    overlay: &'a mut DebugOverlay,
    /// Color for force arrows.
    pub force_color: Color,
    /// Color for torque arrows.
    pub torque_color: Color,
    /// World-unit scale: arrow length = force_magnitude × scale.
    pub force_scale: f32,
    /// World-unit scale for torques.
    pub torque_scale: f32,
}

impl<'a> ForceDebug<'a> {
    /// Create a new `ForceDebug` that draws into `overlay`.
    pub fn new(overlay: &'a mut DebugOverlay) -> Self {
        Self {
            overlay,
            force_color: Color::rgba(1.0, 0.5, 0.0, 1.0),
            torque_color: Color::rgba(0.0, 0.5, 1.0, 1.0),
            force_scale: 0.01,
            torque_scale: 0.05,
        }
    }

    /// Draw a force vector applied at `body_center`.
    pub fn draw_force(&mut self, body_center: [f32; 3], force: [f32; 3]) {
        let mag = vec3_len(force);
        if mag < 1e-10 {
            return;
        }
        let dir = vec3_scale(force, 1.0 / mag);
        let len = mag * self.force_scale;
        self.overlay.arrow(body_center, dir, len, self.force_color);
    }

    /// Draw a torque vector on a body at `body_center`.
    pub fn draw_torque(&mut self, body_center: [f32; 3], torque: [f32; 3]) {
        let mag = vec3_len(torque);
        if mag < 1e-10 {
            return;
        }
        let dir = vec3_scale(torque, 1.0 / mag);
        let len = mag * self.torque_scale;
        self.overlay.arrow(body_center, dir, len, self.torque_color);
    }
}

/// Euclidean length of a 3-vector.
#[inline]
fn vec3_len(v: [f32; 3]) -> f32 {
    (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
}

/// Scale a 3-vector by scalar `s`.
#[inline]
fn vec3_scale(v: [f32; 3], s: f32) -> [f32; 3] {
    [v[0] * s, v[1] * s, v[2] * s]
}

// ─── PerformanceOverlay ───────────────────────────────────────────────────────

/// Tracks per-frame timing statistics for an on-screen performance overlay.
///
/// Call [`PerformanceOverlay::record_frame`] once per rendered frame to keep a
/// rolling window of frame durations.
pub struct PerformanceOverlay {
    /// Rolling history of frame durations in seconds.
    frame_times: std::collections::VecDeque<f32>,
    /// Maximum number of samples stored.
    pub history_size: usize,
    /// Last recorded physics step time in seconds.
    pub last_physics_step_ms: f32,
}

impl PerformanceOverlay {
    /// Create a new overlay with the given rolling-window `history_size`.
    pub fn new(history_size: usize) -> Self {
        Self {
            frame_times: std::collections::VecDeque::with_capacity(history_size),
            history_size: history_size.max(1),
            last_physics_step_ms: 0.0,
        }
    }

    /// Record a frame that took `dt` seconds to render.
    pub fn record_frame(&mut self, dt: f32) {
        if self.frame_times.len() >= self.history_size {
            self.frame_times.pop_front();
        }
        self.frame_times.push_back(dt.max(0.0));
    }

    /// Record the last physics step time in milliseconds.
    pub fn record_physics_step(&mut self, step_ms: f32) {
        self.last_physics_step_ms = step_ms;
    }

    /// Average frames per second over the recorded history.
    pub fn fps(&self) -> f32 {
        if self.frame_times.is_empty() {
            return 0.0;
        }
        let avg_dt = self.avg_frame_time_ms() / 1000.0;
        if avg_dt < 1e-10 { 0.0 } else { 1.0 / avg_dt }
    }

    /// Average frame time in milliseconds over the recorded history.
    pub fn avg_frame_time_ms(&self) -> f32 {
        if self.frame_times.is_empty() {
            return 0.0;
        }
        let sum: f32 = self.frame_times.iter().sum();
        (sum / self.frame_times.len() as f32) * 1000.0
    }

    /// Maximum frame time in milliseconds over the recorded history.
    pub fn max_frame_time_ms(&self) -> f32 {
        self.frame_times.iter().cloned().fold(0.0_f32, f32::max) * 1000.0
    }

    /// Minimum frame time in milliseconds over the recorded history.
    pub fn min_frame_time_ms(&self) -> f32 {
        self.frame_times.iter().cloned().fold(f32::MAX, f32::min) * 1000.0
    }

    /// Number of samples currently recorded.
    pub fn sample_count(&self) -> usize {
        self.frame_times.len()
    }

    /// Returns `true` when no samples have been recorded yet.
    pub fn is_empty(&self) -> bool {
        self.frame_times.is_empty()
    }

    /// Emit text debug commands into `overlay` at `screen_pos`.
    ///
    /// Three lines are emitted: FPS, ms/frame, and physics step time.
    pub fn draw(&self, overlay: &mut DebugOverlay, screen_pos: [f32; 3]) {
        let fps_line = format!("FPS: {:.1}", self.fps());
        let ms_line = format!("Frame: {:.2} ms", self.avg_frame_time_ms());
        let phys_line = format!("Physics: {:.2} ms", self.last_physics_step_ms);

        let dy = 0.03_f32; // vertical offset per line (NDC or world units)
        overlay.text(screen_pos, fps_line, Color::WHITE);
        overlay.text(
            [screen_pos[0], screen_pos[1] - dy, screen_pos[2]],
            ms_line,
            Color::WHITE,
        );
        overlay.text(
            [screen_pos[0], screen_pos[1] - 2.0 * dy, screen_pos[2]],
            phys_line,
            Color::WHITE,
        );
    }
}

// ─── TextPrimitive ────────────────────────────────────────────────────────

/// A simple text rendering primitive for debug overlay.
///
/// Stores a label with position, color, font size, and optional background.
#[derive(Debug, Clone)]
pub struct TextPrimitive {
    /// Screen-space or world-space anchor position.
    pub position: [f32; 3],
    /// The text string.
    pub text: String,
    /// Text color.
    pub color: Color,
    /// Font size in points.
    pub font_size: f32,
    /// Optional background color (None = transparent).
    pub background: Option<Color>,
}

impl TextPrimitive {
    /// Create a basic text primitive.
    pub fn new(position: [f32; 3], text: impl Into<String>, color: Color, font_size: f32) -> Self {
        Self {
            position,
            text: text.into(),
            color,
            font_size,
            background: None,
        }
    }

    /// Set a background color.
    pub fn with_background(mut self, bg: Color) -> Self {
        self.background = Some(bg);
        self
    }

    /// Approximate width in characters (for layout purposes).
    pub fn char_count(&self) -> usize {
        self.text.len()
    }
}

/// A column of text lines for debug display.
#[derive(Debug, Clone)]
pub struct TextColumn {
    /// All text lines in this column.
    pub lines: Vec<TextPrimitive>,
    /// Vertical spacing between lines (world units or NDC).
    pub line_spacing: f32,
}

impl TextColumn {
    /// Create a new text column.
    pub fn new(line_spacing: f32) -> Self {
        Self {
            lines: Vec::new(),
            line_spacing,
        }
    }

    /// Add a line of text.
    pub fn add_line(&mut self, text: impl Into<String>, color: Color, font_size: f32) {
        let y_offset = self.lines.len() as f32 * self.line_spacing;
        self.lines.push(TextPrimitive::new(
            [0.0, -y_offset, 0.0],
            text,
            color,
            font_size,
        ));
    }

    /// Emit all lines into a DebugOverlay at the given base position.
    pub fn draw(&self, overlay: &mut DebugOverlay, base_pos: [f32; 3]) {
        for line in &self.lines {
            let pos = [
                base_pos[0] + line.position[0],
                base_pos[1] + line.position[1],
                base_pos[2] + line.position[2],
            ];
            overlay.text(pos, &line.text, line.color);
        }
    }

    /// Number of lines.
    pub fn line_count(&self) -> usize {
        self.lines.len()
    }
}

// ─── OverlayStatistics ────────────────────────────────────────────────────

/// Collects and displays simulation statistics.
#[derive(Debug, Clone)]
pub struct OverlayStatistics {
    /// Label-value pairs.
    pub entries: Vec<(String, String)>,
}

impl OverlayStatistics {
    /// Create empty stats.
    pub fn new() -> Self {
        Self {
            entries: Vec::new(),
        }
    }

    /// Add a string entry.
    pub fn add(&mut self, label: impl Into<String>, value: impl Into<String>) {
        self.entries.push((label.into(), value.into()));
    }

    /// Add a numeric entry with formatting.
    pub fn add_f64(&mut self, label: impl Into<String>, value: f64, decimals: usize) {
        self.entries
            .push((label.into(), format!("{:.prec$}", value, prec = decimals)));
    }

    /// Draw all entries as text commands.
    pub fn draw(&self, overlay: &mut DebugOverlay, base_pos: [f32; 3], line_height: f32) {
        for (i, (label, value)) in self.entries.iter().enumerate() {
            let y = base_pos[1] - i as f32 * line_height;
            overlay.text(
                [base_pos[0], y, base_pos[2]],
                format!("{}: {}", label, value),
                Color::WHITE,
            );
        }
    }

    /// Number of entries.
    pub fn count(&self) -> usize {
        self.entries.len()
    }
}

impl Default for OverlayStatistics {
    fn default() -> Self {
        Self::new()
    }
}

// ─── ConstraintViz ────────────────────────────────────────────────────────

/// Visualizes constraints (springs, distance constraints) as lines.
pub struct ConstraintViz<'a> {
    overlay: &'a mut DebugOverlay,
    /// Color for satisfied constraints.
    pub satisfied_color: Color,
    /// Color for violated constraints.
    pub violated_color: Color,
    /// Threshold for violation detection (relative error).
    pub violation_threshold: f32,
}

impl<'a> ConstraintViz<'a> {
    /// Create a new constraint visualizer.
    pub fn new(overlay: &'a mut DebugOverlay) -> Self {
        Self {
            overlay,
            satisfied_color: Color::GREEN,
            violated_color: Color::RED,
            violation_threshold: 0.05,
        }
    }

    /// Draw a distance constraint between two points.
    ///
    /// `current_dist` and `rest_dist` determine whether the constraint is
    /// satisfied or violated.
    pub fn draw_distance(&mut self, p1: [f32; 3], p2: [f32; 3], current_dist: f32, rest_dist: f32) {
        let error = ((current_dist - rest_dist) / rest_dist.max(1e-6)).abs();
        let color = if error > self.violation_threshold {
            self.violated_color
        } else {
            self.satisfied_color
        };
        self.overlay.line(p1, p2, color);
    }

    /// Draw a batch of distance constraints.
    ///
    /// Each entry is `(point_a, point_b, current_distance, rest_distance)`.
    pub fn draw_batch(&mut self, constraints: &[([f32; 3], [f32; 3], f32, f32)]) {
        for &(p1, p2, cd, rd) in constraints {
            self.draw_distance(p1, p2, cd, rd);
        }
    }
}

// ─── VelocityArrowViz ─────────────────────────────────────────────────────

/// Visualizes velocity vectors as arrows on bodies.
pub struct VelocityArrowViz<'a> {
    overlay: &'a mut DebugOverlay,
    /// Scale factor for velocity magnitude to arrow length.
    pub scale: f32,
    /// Color for velocity arrows.
    pub color: Color,
}

impl<'a> VelocityArrowViz<'a> {
    /// Create a new velocity arrow visualizer.
    pub fn new(overlay: &'a mut DebugOverlay, scale: f32) -> Self {
        Self {
            overlay,
            scale,
            color: Color::BLUE,
        }
    }

    /// Draw a velocity arrow at `position` for velocity `velocity`.
    pub fn draw(&mut self, position: [f32; 3], velocity: [f32; 3]) {
        let mag = vec3_len(velocity);
        if mag < 1e-10 {
            return;
        }
        let dir = vec3_scale(velocity, 1.0 / mag);
        let len = mag * self.scale;
        self.overlay.arrow(position, dir, len, self.color);
    }

    /// Draw velocity arrows for a batch of (position, velocity) pairs.
    pub fn draw_batch(&mut self, bodies: &[([f32; 3], [f32; 3])]) {
        for &(pos, vel) in bodies {
            self.draw(pos, vel);
        }
    }
}

// ─── ForceDiagram ─────────────────────────────────────────────────────────

/// Draws a force diagram showing all forces acting on a body.
pub struct ForceDiagram<'a> {
    overlay: &'a mut DebugOverlay,
    /// Scale factor for forces.
    pub scale: f32,
}

impl<'a> ForceDiagram<'a> {
    /// Create a new force diagram visualizer.
    pub fn new(overlay: &'a mut DebugOverlay, scale: f32) -> Self {
        Self { overlay, scale }
    }

    /// Draw a force vector with a label.
    pub fn draw_force(&mut self, body_pos: [f32; 3], force: [f32; 3], color: Color, label: &str) {
        let mag = vec3_len(force);
        if mag < 1e-10 {
            return;
        }
        let dir = vec3_scale(force, 1.0 / mag);
        let len = mag * self.scale;
        self.overlay.arrow(body_pos, dir, len, color);
        // Place label at the tip of the arrow
        let tip = [
            body_pos[0] + dir[0] * len,
            body_pos[1] + dir[1] * len,
            body_pos[2] + dir[2] * len,
        ];
        self.overlay.text(tip, label, color);
    }

    /// Draw gravity, normal, and friction forces in one call.
    pub fn draw_free_body(
        &mut self,
        body_pos: [f32; 3],
        gravity: [f32; 3],
        normal: [f32; 3],
        friction: [f32; 3],
    ) {
        self.draw_force(body_pos, gravity, Color::YELLOW, "W");
        self.draw_force(body_pos, normal, Color::GREEN, "N");
        self.draw_force(body_pos, friction, Color::RED, "f");
    }
}

// ─── Grid ─────────────────────────────────────────────────────────────────

impl DebugOverlay {
    /// Draw a ground-plane grid of `n x n` cells centred at the origin.
    ///
    /// The grid lies in the XZ plane at `y = height`.
    pub fn grid(&mut self, n: usize, cell_size: f32, height: f32, color: Color) {
        let half = (n as f32 * cell_size) * 0.5;
        for i in 0..=n {
            let t = -half + i as f32 * cell_size;
            // X-parallel line
            self.line([t, height, -half], [t, height, half], color);
            // Z-parallel line
            self.line([-half, height, t], [half, height, t], color);
        }
    }

    /// Draw coordinate axes at the origin.
    pub fn axes(&mut self, length: f32) {
        self.arrow([0.0; 3], [1.0, 0.0, 0.0], length, Color::RED);
        self.arrow([0.0; 3], [0.0, 1.0, 0.0], length, Color::GREEN);
        self.arrow([0.0; 3], [0.0, 0.0, 1.0], length, Color::BLUE);
    }

    // -----------------------------------------------------------------------
    // BVH tree debug visualization
    // -----------------------------------------------------------------------

    /// Draw the nodes of a Bounding Volume Hierarchy (BVH) as axis-aligned
    /// box outlines.
    ///
    /// Each entry in `bvh_nodes` is a `(min, max)` pair describing one BVH
    /// node's axis-aligned bounding box in world space.  The `level` parameter
    /// selects a colour tint: level 0 is bright green, deeper levels fade to
    /// yellow and red to indicate tree depth.
    ///
    /// At most `max_nodes` entries are drawn to avoid flooding the overlay
    /// command buffer.
    pub fn draw_bvh_tree(
        &mut self,
        bvh_nodes: &[([f32; 3], [f32; 3])],
        levels: &[usize],
        max_nodes: usize,
    ) {
        let n = bvh_nodes.len().min(levels.len()).min(max_nodes);
        for i in 0..n {
            let (min, max) = bvh_nodes[i];
            let level = levels[i];
            // Colour: level 0 = green, deeper = more red.
            let t = (level as f32 * 0.15).clamp(0.0, 1.0);
            let color = Color::rgba(t, 1.0 - t * 0.5, 0.0, 1.0);
            self.aabb(min, max, color);
        }
    }

    // -----------------------------------------------------------------------
    // Velocity field arrows
    // -----------------------------------------------------------------------

    /// Draw velocity arrows for a set of sample points in a vector field.
    ///
    /// `positions` and `velocities` must have the same length.  Each arrow is
    /// scaled by `arrow_scale` so that physically meaningful magnitudes produce
    /// clearly visible arrows.  Arrows shorter than `min_mag` (in world units,
    /// before scaling) are skipped to avoid visual clutter near stagnation
    /// points.
    pub fn draw_velocity_field(
        &mut self,
        positions: &[[f32; 3]],
        velocities: &[[f32; 3]],
        arrow_scale: f32,
        min_mag: f32,
        color: Color,
    ) {
        let n = positions.len().min(velocities.len());
        for i in 0..n {
            let v = velocities[i];
            let mag = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
            if mag < min_mag {
                continue;
            }
            let dir = [v[0] / mag, v[1] / mag, v[2] / mag];
            let len = mag * arrow_scale;
            self.arrow(positions[i], dir, len, color);
        }
    }

    // -----------------------------------------------------------------------
    // Stress tensor glyphs
    // -----------------------------------------------------------------------

    /// Draw stress tensor glyphs as sets of principal-axis arrows.
    ///
    /// Each glyph consists of three arrows along the principal stress
    /// directions, colour-coded by sign: tensile (positive) stresses are
    /// drawn in `tensile_color`, compressive (negative) in
    /// `compressive_color`.  Arrow length is proportional to the eigenvalue
    /// magnitude times `scale`.
    ///
    /// `positions` is the world-space centre of each glyph.
    /// `principal_stresses` contains three principal stress values
    /// `[σ₁, σ₂, σ₃]` for each position.
    /// `principal_axes` contains the corresponding unit eigenvectors
    /// `[[e1x, e1y, e1z\], [e2x, e2y, e2z], [e3x, e3y, e3z]]` per glyph.
    pub fn draw_stress_tensor_glyphs(
        &mut self,
        positions: &[[f32; 3]],
        principal_stresses: &[[f32; 3]],
        principal_axes: &[[[f32; 3]; 3]],
        scale: f32,
        tensile_color: Color,
        compressive_color: Color,
    ) {
        let n = positions
            .len()
            .min(principal_stresses.len())
            .min(principal_axes.len());
        for i in 0..n {
            let pos = positions[i];
            let stresses = principal_stresses[i];
            let axes = principal_axes[i];
            for k in 0..3 {
                let sigma = stresses[k];
                if sigma.abs() < 1e-10 {
                    continue;
                }
                let color = if sigma >= 0.0 {
                    tensile_color
                } else {
                    compressive_color
                };
                let dir = axes[k];
                let len = sigma.abs() * scale;
                // Draw in both directions for a symmetric glyph.
                self.arrow(pos, dir, len, color);
                let neg_dir = [-dir[0], -dir[1], -dir[2]];
                self.arrow(pos, neg_dir, len, color);
            }
        }
    }

    /// Draw a circle (approximated by line segments) in the XZ plane.
    pub fn circle_xz(&mut self, center: [f32; 3], radius: f32, segments: usize, color: Color) {
        let segs = segments.max(3);
        let step = 2.0 * std::f32::consts::PI / segs as f32;
        for i in 0..segs {
            let a0 = i as f32 * step;
            let a1 = (i + 1) as f32 * step;
            let p0 = [
                center[0] + radius * a0.cos(),
                center[1],
                center[2] + radius * a0.sin(),
            ];
            let p1 = [
                center[0] + radius * a1.cos(),
                center[1],
                center[2] + radius * a1.sin(),
            ];
            self.line(p0, p1, color);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_overlay_clear() {
        let mut overlay = DebugOverlay::new(100);
        overlay.line([0.0; 3], [1.0, 0.0, 0.0], Color::RED);
        overlay.sphere([0.0; 3], 1.0, Color::GREEN);
        assert_eq!(overlay.len(), 2);
        overlay.clear();
        assert!(overlay.is_empty());
    }

    #[test]
    fn test_overlay_aabb_expands_to_lines() {
        let mut overlay = DebugOverlay::new(100);
        overlay.aabb([-1.0; 3], [1.0; 3], Color::WHITE);
        assert_eq!(
            overlay.len(),
            12,
            "aabb() must emit exactly 12 line commands"
        );
        for cmd in overlay.commands() {
            assert!(matches!(cmd, DebugCommand::Line { .. }));
        }
    }

    #[test]
    fn test_overlay_max_commands() {
        let mut overlay = DebugOverlay::new(3);
        for _ in 0..10 {
            overlay.line([0.0; 3], [1.0, 0.0, 0.0], Color::BLUE);
        }
        assert_eq!(overlay.len(), 3, "overlay must not exceed max_commands");
        assert!(overlay.is_full());
    }

    #[test]
    fn test_color_constants() {
        assert!((Color::RED.r - 1.0).abs() < 1e-6);
        assert!(Color::RED.g.abs() < 1e-6);
        assert!(Color::RED.b.abs() < 1e-6);

        assert!((Color::GREEN.g - 1.0).abs() < 1e-6);
        assert!(Color::GREEN.r.abs() < 1e-6);

        assert!((Color::BLUE.b - 1.0).abs() < 1e-6);
        assert!((Color::WHITE.r - 1.0).abs() < 1e-6);
        assert!((Color::WHITE.g - 1.0).abs() < 1e-6);
        assert!((Color::WHITE.b - 1.0).abs() < 1e-6);

        assert!((Color::YELLOW.r - 1.0).abs() < 1e-6);
        assert!((Color::YELLOW.g - 1.0).abs() < 1e-6);
        assert!(Color::YELLOW.b.abs() < 1e-6);
    }

    // ── ContactDebug ─────────────────────────────────────────────────────────

    #[test]
    fn test_contact_debug_emits_point_and_arrow() {
        let mut overlay = DebugOverlay::new(100);
        {
            let mut cd = ContactDebug::new(&mut overlay);
            cd.draw_contact([0.0; 3], [0.0, 1.0, 0.0], 0.1);
        }
        // Should emit 1 point + 1 arrow = 2 commands.
        assert_eq!(overlay.len(), 2);
        assert!(matches!(overlay.commands()[0], DebugCommand::Point { .. }));
        assert!(matches!(overlay.commands()[1], DebugCommand::Arrow { .. }));
    }

    #[test]
    fn test_contact_debug_batch() {
        let mut overlay = DebugOverlay::new(100);
        {
            let mut cd = ContactDebug::new(&mut overlay);
            let contacts = vec![
                ([0.0; 3], [0.0, 1.0, 0.0], 0.05_f32),
                ([1.0, 0.0, 0.0], [1.0, 0.0, 0.0], 0.1_f32),
            ];
            cd.draw_contacts(&contacts);
        }
        // 2 contacts × 2 commands each = 4.
        assert_eq!(overlay.len(), 4);
    }

    // ── ForceDebug ────────────────────────────────────────────────────────────

    #[test]
    fn test_force_debug_emits_arrow_for_nonzero_force() {
        let mut overlay = DebugOverlay::new(100);
        {
            let mut fd = ForceDebug::new(&mut overlay);
            fd.draw_force([0.0; 3], [10.0, 0.0, 0.0]);
        }
        assert_eq!(overlay.len(), 1);
        assert!(matches!(overlay.commands()[0], DebugCommand::Arrow { .. }));
    }

    #[test]
    fn test_force_debug_zero_force_no_command() {
        let mut overlay = DebugOverlay::new(100);
        {
            let mut fd = ForceDebug::new(&mut overlay);
            fd.draw_force([0.0; 3], [0.0, 0.0, 0.0]);
        }
        assert!(overlay.is_empty(), "zero force should not emit any command");
    }

    #[test]
    fn test_torque_debug_emits_arrow() {
        let mut overlay = DebugOverlay::new(100);
        {
            let mut fd = ForceDebug::new(&mut overlay);
            fd.draw_torque([0.0; 3], [0.0, 5.0, 0.0]);
        }
        assert_eq!(overlay.len(), 1);
        assert!(matches!(overlay.commands()[0], DebugCommand::Arrow { .. }));
    }

    // ── PerformanceOverlay ────────────────────────────────────────────────────

    #[test]
    fn test_perf_overlay_fps_constant_framerate() {
        let mut perf = PerformanceOverlay::new(10);
        for _ in 0..10 {
            perf.record_frame(1.0 / 60.0); // 60 fps
        }
        let fps = perf.fps();
        assert!((fps - 60.0).abs() < 0.5, "expected ~60 fps, got {fps}");
    }

    #[test]
    fn test_perf_overlay_empty_fps_zero() {
        let perf = PerformanceOverlay::new(10);
        assert!(perf.is_empty());
        assert_eq!(perf.fps(), 0.0);
    }

    #[test]
    fn test_perf_overlay_rolling_window() {
        let mut perf = PerformanceOverlay::new(3);
        for i in 0..10 {
            perf.record_frame(i as f32 * 0.01 + 0.001);
        }
        assert_eq!(
            perf.sample_count(),
            3,
            "rolling window should cap at history_size"
        );
    }

    #[test]
    fn test_perf_overlay_draw_emits_three_text_commands() {
        let mut perf = PerformanceOverlay::new(10);
        perf.record_frame(1.0 / 30.0);
        perf.record_physics_step(2.5);
        let mut overlay = DebugOverlay::new(100);
        perf.draw(&mut overlay, [0.0; 3]);
        assert_eq!(overlay.len(), 3, "draw() should emit 3 Text commands");
        for cmd in overlay.commands() {
            assert!(matches!(cmd, DebugCommand::Text { .. }));
        }
    }

    #[test]
    fn test_perf_overlay_max_min_frame_time() {
        let mut perf = PerformanceOverlay::new(10);
        perf.record_frame(0.010); // 10 ms
        perf.record_frame(0.020); // 20 ms
        perf.record_frame(0.015); // 15 ms
        assert!((perf.max_frame_time_ms() - 20.0).abs() < 0.1);
        assert!((perf.min_frame_time_ms() - 10.0).abs() < 0.1);
    }

    // ── TextPrimitive ─────────────────────────────────────────────────────

    #[test]
    fn test_text_primitive_char_count() {
        let tp = TextPrimitive::new([0.0; 3], "Hello", Color::WHITE, 12.0);
        assert_eq!(tp.char_count(), 5);
    }

    #[test]
    fn test_text_primitive_with_background() {
        let tp =
            TextPrimitive::new([0.0; 3], "Test", Color::WHITE, 12.0).with_background(Color::RED);
        assert!(tp.background.is_some());
    }

    // ── TextColumn ────────────────────────────────────────────────────────

    #[test]
    fn test_text_column_add_lines() {
        let mut col = TextColumn::new(0.03);
        col.add_line("Line 1", Color::WHITE, 12.0);
        col.add_line("Line 2", Color::WHITE, 12.0);
        assert_eq!(col.line_count(), 2);
    }

    #[test]
    fn test_text_column_draw() {
        let mut col = TextColumn::new(0.03);
        col.add_line("FPS: 60", Color::GREEN, 12.0);
        col.add_line("Objects: 100", Color::WHITE, 12.0);
        let mut overlay = DebugOverlay::new(100);
        col.draw(&mut overlay, [0.0; 3]);
        assert_eq!(overlay.len(), 2);
    }

    // ── OverlayStatistics ─────────────────────────────────────────────────

    #[test]
    fn test_overlay_statistics_add() {
        let mut stats = OverlayStatistics::new();
        stats.add("Bodies", "42");
        stats.add_f64("Energy", 123.456, 2);
        assert_eq!(stats.count(), 2);
    }

    #[test]
    fn test_overlay_statistics_draw() {
        let mut stats = OverlayStatistics::new();
        stats.add("Test", "Value");
        let mut overlay = DebugOverlay::new(100);
        stats.draw(&mut overlay, [0.0; 3], 0.03);
        assert_eq!(overlay.len(), 1);
    }

    #[test]
    fn test_overlay_statistics_default() {
        let stats: OverlayStatistics = Default::default();
        assert_eq!(stats.count(), 0);
    }

    // ── ConstraintViz ─────────────────────────────────────────────────────

    #[test]
    fn test_constraint_viz_satisfied() {
        let mut overlay = DebugOverlay::new(100);
        {
            let mut cv = ConstraintViz::new(&mut overlay);
            cv.draw_distance([0.0; 3], [1.0, 0.0, 0.0], 1.0, 1.0);
        }
        assert_eq!(overlay.len(), 1);
        // Should be green (satisfied)
        if let DebugCommand::Line { color, .. } = &overlay.commands()[0] {
            assert!((color.g - 1.0).abs() < 1e-6);
        }
    }

    #[test]
    fn test_constraint_viz_violated() {
        let mut overlay = DebugOverlay::new(100);
        {
            let mut cv = ConstraintViz::new(&mut overlay);
            cv.draw_distance([0.0; 3], [1.0, 0.0, 0.0], 2.0, 1.0); // 100% error
        }
        assert_eq!(overlay.len(), 1);
        // Should be red (violated)
        if let DebugCommand::Line { color, .. } = &overlay.commands()[0] {
            assert!((color.r - 1.0).abs() < 1e-6);
        }
    }

    #[test]
    fn test_constraint_viz_batch() {
        let mut overlay = DebugOverlay::new(100);
        {
            let mut cv = ConstraintViz::new(&mut overlay);
            let constraints = vec![
                ([0.0; 3], [1.0, 0.0, 0.0], 1.0_f32, 1.0_f32),
                ([0.0; 3], [2.0, 0.0, 0.0], 2.0, 1.0),
            ];
            cv.draw_batch(&constraints);
        }
        assert_eq!(overlay.len(), 2);
    }

    // ── VelocityArrowViz ──────────────────────────────────────────────────

    #[test]
    fn test_velocity_arrow_viz() {
        let mut overlay = DebugOverlay::new(100);
        {
            let mut vv = VelocityArrowViz::new(&mut overlay, 0.1);
            vv.draw([0.0; 3], [10.0, 0.0, 0.0]);
        }
        assert_eq!(overlay.len(), 1);
    }

    #[test]
    fn test_velocity_arrow_zero_velocity() {
        let mut overlay = DebugOverlay::new(100);
        {
            let mut vv = VelocityArrowViz::new(&mut overlay, 0.1);
            vv.draw([0.0; 3], [0.0, 0.0, 0.0]);
        }
        assert!(overlay.is_empty(), "Zero velocity should produce no arrow");
    }

    #[test]
    fn test_velocity_arrow_batch() {
        let mut overlay = DebugOverlay::new(100);
        {
            let mut vv = VelocityArrowViz::new(&mut overlay, 0.1);
            vv.draw_batch(&[
                ([0.0; 3], [1.0, 0.0, 0.0]),
                ([1.0, 0.0, 0.0], [0.0, 1.0, 0.0]),
            ]);
        }
        assert_eq!(overlay.len(), 2);
    }

    // ── ForceDiagram ──────────────────────────────────────────────────────

    #[test]
    fn test_force_diagram_draw_force() {
        let mut overlay = DebugOverlay::new(100);
        {
            let mut fd = ForceDiagram::new(&mut overlay, 0.01);
            fd.draw_force([0.0; 3], [0.0, -9810.0, 0.0], Color::YELLOW, "W");
        }
        // 1 arrow + 1 text = 2
        assert_eq!(overlay.len(), 2);
    }

    #[test]
    fn test_force_diagram_free_body() {
        let mut overlay = DebugOverlay::new(100);
        {
            let mut fd = ForceDiagram::new(&mut overlay, 0.01);
            fd.draw_free_body(
                [0.0; 3],
                [0.0, -9810.0, 0.0],
                [0.0, 9810.0, 0.0],
                [100.0, 0.0, 0.0],
            );
        }
        // 3 forces * 2 commands (arrow + text) = 6
        assert_eq!(overlay.len(), 6);
    }

    // ── Grid and Axes ─────────────────────────────────────────────────────

    #[test]
    fn test_grid_draws_lines() {
        let mut overlay = DebugOverlay::new(1000);
        overlay.grid(4, 1.0, 0.0, Color::WHITE);
        // 4x4 grid → (4+1)*2 = 10 lines
        assert_eq!(overlay.len(), 10);
    }

    #[test]
    fn test_axes_draws_three_arrows() {
        let mut overlay = DebugOverlay::new(100);
        overlay.axes(1.0);
        assert_eq!(overlay.len(), 3);
    }

    #[test]
    fn test_circle_xz_draws_segments() {
        let mut overlay = DebugOverlay::new(1000);
        overlay.circle_xz([0.0; 3], 1.0, 16, Color::WHITE);
        assert_eq!(overlay.len(), 16);
    }

    // ── Color ─────────────────────────────────────────────────────────────

    #[test]
    fn test_color_with_alpha() {
        let c = Color::RED.with_alpha(0.5);
        assert!((c.a - 0.5).abs() < 1e-6);
        assert!((c.r - 1.0).abs() < 1e-6);
    }

    // ── draw_bvh_tree ──────────────────────────────────────────────────────

    #[test]
    fn test_draw_bvh_tree_emits_boxes() {
        let mut overlay = DebugOverlay::new(10000);
        let nodes = vec![([-1.0_f32; 3], [1.0_f32; 3]), ([-0.5_f32; 3], [0.5_f32; 3])];
        let levels = vec![0usize, 1];
        overlay.draw_bvh_tree(&nodes, &levels, 10);
        // 2 AABBs × 12 lines each = 24 lines
        assert_eq!(
            overlay.len(),
            24,
            "BVH tree with 2 nodes should emit 24 lines"
        );
    }

    #[test]
    fn test_draw_bvh_tree_max_nodes_limit() {
        let mut overlay = DebugOverlay::new(10000);
        let nodes: Vec<_> = (0..10).map(|_| ([-1.0_f32; 3], [1.0_f32; 3])).collect();
        let levels: Vec<_> = (0..10).collect();
        overlay.draw_bvh_tree(&nodes, &levels, 3);
        // Only 3 AABBs drawn (max_nodes=3) → 36 lines
        assert_eq!(overlay.len(), 36, "max_nodes=3 should draw exactly 3 boxes");
    }

    #[test]
    fn test_draw_bvh_tree_empty_no_commands() {
        let mut overlay = DebugOverlay::new(1000);
        overlay.draw_bvh_tree(&[], &[], 100);
        assert!(overlay.is_empty());
    }

    // ── draw_velocity_field ────────────────────────────────────────────────

    #[test]
    fn test_velocity_field_emits_arrows() {
        let mut overlay = DebugOverlay::new(1000);
        let pos = vec![[0.0_f32; 3], [1.0, 0.0, 0.0]];
        let vel = vec![[1.0_f32, 0.0, 0.0], [0.0, 1.0, 0.0]];
        overlay.draw_velocity_field(&pos, &vel, 1.0, 0.01, Color::BLUE);
        assert_eq!(
            overlay.len(),
            2,
            "Two velocity arrows should produce 2 commands"
        );
    }

    #[test]
    fn test_velocity_field_skips_low_magnitude() {
        let mut overlay = DebugOverlay::new(1000);
        let pos = vec![[0.0_f32; 3]];
        let vel = vec![[0.0_f32, 0.0, 0.0001]]; // below min_mag=0.01
        overlay.draw_velocity_field(&pos, &vel, 1.0, 0.01, Color::RED);
        assert!(
            overlay.is_empty(),
            "Sub-threshold velocity should be skipped"
        );
    }

    #[test]
    fn test_velocity_field_empty_input() {
        let mut overlay = DebugOverlay::new(1000);
        overlay.draw_velocity_field(&[], &[], 1.0, 0.0, Color::WHITE);
        assert!(overlay.is_empty());
    }

    // ── draw_stress_tensor_glyphs ──────────────────────────────────────────

    #[test]
    fn test_stress_glyphs_emit_arrows() {
        let mut overlay = DebugOverlay::new(10000);
        let positions = vec![[0.0_f32; 3]];
        let stresses = vec![[1.0_f32, -0.5, 0.0]]; // σ3=0 skipped
        let axes = vec![[[1.0_f32, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]];
        overlay.draw_stress_tensor_glyphs(
            &positions,
            &stresses,
            &axes,
            1.0,
            Color::RED,
            Color::BLUE,
        );
        // σ1=1 → 2 arrows, σ2=-0.5 → 2 arrows, σ3=0 → skipped = 4 arrows
        assert_eq!(
            overlay.len(),
            4,
            "2 non-zero principal stresses × 2 directions = 4 arrows"
        );
    }

    #[test]
    fn test_stress_glyphs_zero_stresses_no_output() {
        let mut overlay = DebugOverlay::new(1000);
        let positions = vec![[0.0_f32; 3]];
        let stresses = vec![[0.0_f32; 3]];
        let axes = vec![[[1.0_f32, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]];
        overlay.draw_stress_tensor_glyphs(
            &positions,
            &stresses,
            &axes,
            1.0,
            Color::RED,
            Color::BLUE,
        );
        assert!(
            overlay.is_empty(),
            "All-zero stresses should produce no glyphs"
        );
    }

    #[test]
    fn test_stress_glyphs_empty_input() {
        let mut overlay = DebugOverlay::new(1000);
        overlay.draw_stress_tensor_glyphs(&[], &[], &[], 1.0, Color::WHITE, Color::WHITE);
        assert!(overlay.is_empty());
    }
}
