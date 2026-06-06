//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::*;
use std::f64::consts::PI;

/// A leader line connecting a point of interest to its annotation label.
#[derive(Debug, Clone)]
pub struct LeaderLine {
    /// Start point (the feature being annotated).
    pub start: [f64; 3],
    /// End point (where the label sits).
    pub end: [f64; 3],
    /// Optional intermediate waypoints for bent leader lines.
    pub waypoints: Vec<[f64; 3]>,
    /// Arrow head size at the start.
    pub arrow_size: f64,
    /// Whether to show an arrowhead.
    pub show_arrow: bool,
}
impl LeaderLine {
    /// Create a straight leader line.
    pub fn new(start: [f64; 3], end: [f64; 3]) -> Self {
        Self {
            start,
            end,
            waypoints: Vec::new(),
            arrow_size: 0.08,
            show_arrow: true,
        }
    }
    /// Create a bent leader line with one intermediate waypoint.
    pub fn with_bend(start: [f64; 3], bend: [f64; 3], end: [f64; 3]) -> Self {
        Self {
            start,
            end,
            waypoints: vec![bend],
            arrow_size: 0.08,
            show_arrow: true,
        }
    }
    /// Total length of the leader line path.
    pub fn length(&self) -> f64 {
        let mut total = 0.0;
        let mut prev = self.start;
        for &wp in &self.waypoints {
            total += dist3(prev, wp);
            prev = wp;
        }
        total += dist3(prev, self.end);
        total
    }
    /// Number of line segments in the leader.
    pub fn segment_count(&self) -> usize {
        self.waypoints.len() + 1
    }
    /// Get all points in the leader path (start, waypoints, end).
    pub fn path_points(&self) -> Vec<[f64; 3]> {
        let mut pts = Vec::with_capacity(self.waypoints.len() + 2);
        pts.push(self.start);
        pts.extend_from_slice(&self.waypoints);
        pts.push(self.end);
        pts
    }
    /// Direction of the first segment (from start).
    pub fn initial_direction(&self) -> [f64; 3] {
        let next = if self.waypoints.is_empty() {
            self.end
        } else {
            self.waypoints[0]
        };
        normalize3(sub3(next, self.start))
    }
}
/// A complete area measurement annotation.
#[derive(Debug, Clone)]
pub struct AreaMeasurement {
    /// Polygon vertices.
    pub points: Vec<[f64; 3]>,
    /// Precision (decimal places).
    pub precision: usize,
    /// Unit label.
    pub unit: String,
}
impl AreaMeasurement {
    /// Create an area measurement for a polygon.
    pub fn new(points: Vec<[f64; 3]>, unit: &str) -> Self {
        Self {
            points,
            precision: 3,
            unit: unit.to_string(),
        }
    }
    /// Measured area.
    pub fn area(&self) -> f64 {
        measure_polygon_area(&self.points)
    }
    /// Perimeter of the polygon.
    pub fn perimeter(&self) -> f64 {
        measure_polygon_perimeter(&self.points)
    }
    /// Centroid of the polygon.
    pub fn centroid(&self) -> [f64; 3] {
        compute_centroid(&self.points)
    }
    /// Formatted display string.
    pub fn display(&self) -> String {
        format!(
            "{:.prec$} {}\u{00b2}",
            self.area(),
            self.unit,
            prec = self.precision
        )
    }
}
/// An arc annotation showing an angle between two rays emanating from a vertex.
#[derive(Debug, Clone)]
pub struct AngleArc {
    /// The vertex point where the two rays meet.
    pub vertex: [f64; 3],
    /// Direction of the first ray (normalized).
    pub ray_a: [f64; 3],
    /// Direction of the second ray (normalized).
    pub ray_b: [f64; 3],
    /// Radius of the arc from the vertex.
    pub radius: f64,
    /// Number of segments to approximate the arc.
    pub segments: usize,
}
impl AngleArc {
    /// Create a new angle arc annotation.
    pub fn new(vertex: [f64; 3], ray_a: [f64; 3], ray_b: [f64; 3], radius: f64) -> Self {
        Self {
            vertex,
            ray_a: normalize3(ray_a),
            ray_b: normalize3(ray_b),
            radius,
            segments: 16,
        }
    }
    /// Compute the angle in radians between the two rays.
    pub fn angle_rad(&self) -> f64 {
        let d = dot3(self.ray_a, self.ray_b).clamp(-1.0, 1.0);
        d.acos()
    }
    /// Compute the angle in degrees.
    pub fn angle_deg(&self) -> f64 {
        self.angle_rad().to_degrees()
    }
    /// Generate the arc polyline points using spherical interpolation.
    pub fn arc_points(&self) -> Vec<[f64; 3]> {
        let angle = self.angle_rad();
        if angle.abs() < 1e-12 {
            return vec![add3(self.vertex, scale3(self.ray_a, self.radius))];
        }
        let n = self.segments.max(2);
        let sin_angle = angle.sin();
        let mut pts = Vec::with_capacity(n + 1);
        for i in 0..=n {
            let t = i as f64 / n as f64;
            let theta = t * angle;
            let w_a = ((1.0 - t) * angle).sin() / sin_angle;
            let w_b = (t * angle).sin() / sin_angle;
            let dir = add3(scale3(self.ray_a, w_a), scale3(self.ray_b, w_b));
            let _ = theta;
            pts.push(add3(self.vertex, scale3(dir, self.radius)));
        }
        pts
    }
    /// Midpoint direction of the arc.
    pub fn mid_direction(&self) -> [f64; 3] {
        normalize3(add3(self.ray_a, self.ray_b))
    }
    /// Position at the midpoint of the arc (for label placement).
    pub fn label_position(&self) -> [f64; 3] {
        add3(self.vertex, scale3(self.mid_direction(), self.radius * 1.3))
    }
    /// Display label for the angle.
    pub fn display_label(&self) -> String {
        format!("{:.1}\u{00b0}", self.angle_deg())
    }
    /// Arc length.
    pub fn arc_length(&self) -> f64 {
        self.radius * self.angle_rad()
    }
}
/// Style configuration for annotations.
#[derive(Debug, Clone)]
pub struct AnnotationStyle {
    /// Line width for dimension and leader lines.
    pub line_width: f64,
    /// Arrow head length.
    pub arrow_length: f64,
    /// Arrow head width.
    pub arrow_width: f64,
    /// Font size for labels.
    pub font_size: f64,
    /// Color as RGBA \[0..1\].
    pub color: [f64; 4],
    /// Text background color (RGBA).
    pub text_bg_color: [f64; 4],
    /// Whether to show text background.
    pub show_text_bg: bool,
}
impl AnnotationStyle {
    /// Create an annotation style with a custom color.
    pub fn with_color(color: [f64; 4]) -> Self {
        Self {
            color,
            ..Self::default()
        }
    }
    /// Compute the arrow head area (approximate triangle).
    pub fn arrow_area(&self) -> f64 {
        0.5 * self.arrow_length * self.arrow_width
    }
    /// Create a style for dimensioning.
    pub fn dimension_style() -> Self {
        Self {
            line_width: 0.8,
            arrow_length: 0.12,
            arrow_width: 0.04,
            font_size: 10.0,
            color: [0.0, 0.0, 0.0, 1.0],
            text_bg_color: [1.0, 1.0, 1.0, 0.9],
            show_text_bg: true,
        }
    }
    /// Create a style for callouts.
    pub fn callout_style() -> Self {
        Self {
            line_width: 1.2,
            arrow_length: 0.08,
            arrow_width: 0.04,
            font_size: 11.0,
            color: [0.2, 0.2, 0.8, 1.0],
            text_bg_color: [0.95, 0.95, 1.0, 0.95],
            show_text_bg: true,
        }
    }
}
/// A complete distance measurement annotation.
#[derive(Debug, Clone)]
pub struct DistanceMeasurement {
    /// Dimension line.
    pub dimension: DimensionLine,
    /// Precision (decimal places).
    pub precision: usize,
    /// Unit label.
    pub unit: String,
}
impl DistanceMeasurement {
    /// Create a distance measurement between two points.
    pub fn new(start: [f64; 3], end: [f64; 3], offset: f64, unit: &str) -> Self {
        Self {
            dimension: DimensionLine::new(start, end, offset),
            precision: 3,
            unit: unit.to_string(),
        }
    }
    /// Measured distance.
    pub fn distance(&self) -> f64 {
        self.dimension.length()
    }
    /// Formatted display string.
    pub fn display(&self) -> String {
        let val = self.distance();
        format!("{:.prec$} {}", val, self.unit, prec = self.precision)
    }
}
/// A complete angle measurement annotation.
#[derive(Debug, Clone)]
pub struct AngleMeasurement {
    /// The angle arc.
    pub arc: AngleArc,
    /// Precision (decimal places).
    pub precision: usize,
}
impl AngleMeasurement {
    /// Create an angle measurement.
    pub fn new(vertex: [f64; 3], ray_a: [f64; 3], ray_b: [f64; 3], radius: f64) -> Self {
        Self {
            arc: AngleArc::new(vertex, ray_a, ray_b, radius),
            precision: 1,
        }
    }
    /// Measured angle in degrees.
    pub fn angle_deg(&self) -> f64 {
        self.arc.angle_deg()
    }
    /// Formatted display string.
    pub fn display(&self) -> String {
        format!("{:.prec$}\u{00b0}", self.angle_deg(), prec = self.precision)
    }
}
/// Anchor point for text or label placement.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum TextAnchor {
    /// Top-left corner.
    TopLeft,
    /// Top-center.
    TopCenter,
    /// Top-right corner.
    TopRight,
    /// Middle-left.
    MiddleLeft,
    /// Center.
    #[default]
    Center,
    /// Middle-right.
    MiddleRight,
    /// Bottom-left corner.
    BottomLeft,
    /// Bottom-center.
    BottomCenter,
    /// Bottom-right corner.
    BottomRight,
}
impl TextAnchor {
    /// Compute the offset from the anchor to the center of a box with given
    /// `width` and `height`.
    pub fn offset_to_center(&self, width: f64, height: f64) -> [f64; 2] {
        let hw = width * 0.5;
        let hh = height * 0.5;
        match self {
            Self::TopLeft => [hw, -hh],
            Self::TopCenter => [0.0, -hh],
            Self::TopRight => [-hw, -hh],
            Self::MiddleLeft => [hw, 0.0],
            Self::Center => [0.0, 0.0],
            Self::MiddleRight => [-hw, 0.0],
            Self::BottomLeft => [hw, hh],
            Self::BottomCenter => [0.0, hh],
            Self::BottomRight => [-hw, hh],
        }
    }
}
/// A scale bar indicating real-world distance on a drawing.
#[derive(Debug, Clone)]
pub struct ScaleBar {
    /// Bottom-left corner position.
    pub position: [f64; 3],
    /// Length in world units.
    pub length: f64,
    /// Number of divisions.
    pub divisions: usize,
    /// Height of the bar.
    pub bar_height: f64,
    /// Unit label (e.g., "m", "mm").
    pub unit: String,
    /// Scale factor (e.g., 1:100).
    pub scale_factor: f64,
}
impl ScaleBar {
    /// Create a new scale bar.
    pub fn new(position: [f64; 3], length: f64, divisions: usize, unit: &str) -> Self {
        Self {
            position,
            length,
            divisions: divisions.max(1),
            bar_height: length * 0.05,
            unit: unit.to_string(),
            scale_factor: 1.0,
        }
    }
    /// Division width.
    pub fn division_width(&self) -> f64 {
        self.length / self.divisions as f64
    }
    /// Tick positions along the scale bar.
    pub fn tick_positions(&self) -> Vec<[f64; 3]> {
        let mut ticks = Vec::with_capacity(self.divisions + 1);
        for i in 0..=self.divisions {
            let x = self.position[0] + i as f64 * self.division_width();
            ticks.push([x, self.position[1], self.position[2]]);
        }
        ticks
    }
    /// Tick labels (real-world distance values).
    pub fn tick_labels(&self) -> Vec<String> {
        (0..=self.divisions)
            .map(|i| {
                let val = i as f64 * self.division_width() * self.scale_factor;
                format!("{:.1} {}", val, self.unit)
            })
            .collect()
    }
    /// Bounding box of the scale bar.
    pub fn bounding_box(&self) -> BoundingBox2D {
        BoundingBox2D::new(
            [self.position[0], self.position[1]],
            [
                self.position[0] + self.length,
                self.position[1] + self.bar_height,
            ],
        )
    }
    /// End position (right end of the scale bar).
    pub fn end_position(&self) -> [f64; 3] {
        [
            self.position[0] + self.length,
            self.position[1],
            self.position[2],
        ]
    }
    /// Total represented real-world length.
    pub fn real_world_length(&self) -> f64 {
        self.length * self.scale_factor
    }
}
/// Result of a label placement attempt.
#[derive(Debug, Clone)]
pub struct PlacedLabel {
    /// Position of the label center.
    pub position: [f64; 2],
    /// Bounding box of the label.
    pub bounds: BoundingBox2D,
    /// The original (ideal) position before adjustment.
    pub ideal_position: [f64; 2],
    /// Whether the label was displaced from its ideal position.
    pub displaced: bool,
}
/// A polyline annotation (chain of connected line segments).
#[derive(Debug, Clone)]
pub struct PolylineAnnotation {
    /// Points of the polyline.
    pub points: Vec<[f64; 3]>,
    /// Whether the polyline is closed.
    pub closed: bool,
    /// Label placed at the centroid.
    pub label: String,
}
impl PolylineAnnotation {
    /// Create a new polyline annotation.
    pub fn new(points: Vec<[f64; 3]>, closed: bool) -> Self {
        Self {
            points,
            closed,
            label: String::new(),
        }
    }
    /// Total length of the polyline.
    pub fn length(&self) -> f64 {
        if self.points.len() < 2 {
            return 0.0;
        }
        let mut len = 0.0;
        for i in 0..self.points.len() - 1 {
            len += dist3(self.points[i], self.points[i + 1]);
        }
        if self.closed && self.points.len() >= 2 {
            len += dist3(
                *self.points.last().expect("collection should not be empty"),
                self.points[0],
            );
        }
        len
    }
    /// Number of segments.
    pub fn segment_count(&self) -> usize {
        if self.points.len() < 2 {
            return 0;
        }
        let base = self.points.len() - 1;
        if self.closed { base + 1 } else { base }
    }
    /// Centroid of the polyline points.
    pub fn centroid(&self) -> [f64; 3] {
        compute_centroid(&self.points)
    }
}
/// A container for annotations that manages placement and collision avoidance.
#[derive(Debug, Clone)]
pub struct AnnotationLayer {
    /// Dimension lines.
    pub dimensions: Vec<DimensionLine>,
    /// Leader lines.
    pub leaders: Vec<LeaderLine>,
    /// Callout boxes.
    pub callouts: Vec<CalloutBox>,
    /// Angle arcs.
    pub angle_arcs: Vec<AngleArc>,
    /// Scale bars.
    pub scale_bars: Vec<ScaleBar>,
    /// Annotation anchors.
    pub anchors: Vec<AnnotationAnchor>,
}
impl AnnotationLayer {
    /// Create an empty annotation layer.
    pub fn new() -> Self {
        Self {
            dimensions: Vec::new(),
            leaders: Vec::new(),
            callouts: Vec::new(),
            angle_arcs: Vec::new(),
            scale_bars: Vec::new(),
            anchors: Vec::new(),
        }
    }
    /// Add a dimension line.
    pub fn add_dimension(&mut self, dim: DimensionLine) {
        self.dimensions.push(dim);
    }
    /// Add a leader line.
    pub fn add_leader(&mut self, leader: LeaderLine) {
        self.leaders.push(leader);
    }
    /// Add a callout box.
    pub fn add_callout(&mut self, callout: CalloutBox) {
        self.callouts.push(callout);
    }
    /// Add an angle arc.
    pub fn add_angle_arc(&mut self, arc: AngleArc) {
        self.angle_arcs.push(arc);
    }
    /// Add a scale bar.
    pub fn add_scale_bar(&mut self, bar: ScaleBar) {
        self.scale_bars.push(bar);
    }
    /// Add an annotation anchor.
    pub fn add_anchor(&mut self, anchor: AnnotationAnchor) {
        self.anchors.push(anchor);
    }
    /// Total number of annotations.
    pub fn annotation_count(&self) -> usize {
        self.dimensions.len()
            + self.leaders.len()
            + self.callouts.len()
            + self.angle_arcs.len()
            + self.scale_bars.len()
    }
    /// Clear all annotations.
    pub fn clear(&mut self) {
        self.dimensions.clear();
        self.leaders.clear();
        self.callouts.clear();
        self.angle_arcs.clear();
        self.scale_bars.clear();
        self.anchors.clear();
    }
    /// Collect all callout bounding boxes for collision detection.
    pub fn callout_bounds(&self) -> Vec<BoundingBox2D> {
        self.callouts.iter().map(|c| c.bounding_box()).collect()
    }
    /// Update all anchors with new feature positions.
    pub fn update_anchors(&mut self, new_positions: &[[f64; 3]]) {
        for (anchor, &pos) in self.anchors.iter_mut().zip(new_positions.iter()) {
            anchor.update(pos);
        }
    }
}
/// A small coordinate axes widget for orientation reference.
#[derive(Debug, Clone)]
pub struct CoordinateAxesWidget {
    /// Origin position.
    pub origin: [f64; 3],
    /// Length of each axis arm.
    pub arm_length: f64,
    /// Arrow head size.
    pub arrow_size: f64,
    /// Labels for each axis.
    pub labels: [String; 3],
    /// Whether to show axis labels.
    pub show_labels: bool,
}
impl CoordinateAxesWidget {
    /// Create a default XYZ axes widget at the given origin.
    pub fn new(origin: [f64; 3], arm_length: f64) -> Self {
        Self {
            origin,
            arm_length,
            arrow_size: arm_length * 0.15,
            labels: ["X".to_string(), "Y".to_string(), "Z".to_string()],
            show_labels: true,
        }
    }
    /// Tip positions of the three axes.
    pub fn axis_tips(&self) -> [[f64; 3]; 3] {
        [
            add3(self.origin, [self.arm_length, 0.0, 0.0]),
            add3(self.origin, [0.0, self.arm_length, 0.0]),
            add3(self.origin, [0.0, 0.0, self.arm_length]),
        ]
    }
    /// Label positions (slightly beyond the tips).
    pub fn label_positions(&self) -> [[f64; 3]; 3] {
        let ext = self.arm_length * 1.15;
        [
            add3(self.origin, [ext, 0.0, 0.0]),
            add3(self.origin, [0.0, ext, 0.0]),
            add3(self.origin, [0.0, 0.0, ext]),
        ]
    }
    /// Line segments for the axes: returns 3 pairs (origin, tip).
    pub fn axis_lines(&self) -> [([f64; 3], [f64; 3]); 3] {
        let tips = self.axis_tips();
        [
            (self.origin, tips[0]),
            (self.origin, tips[1]),
            (self.origin, tips[2]),
        ]
    }
    /// Set custom labels.
    pub fn with_labels(mut self, x: &str, y: &str, z: &str) -> Self {
        self.labels = [x.to_string(), y.to_string(), z.to_string()];
        self
    }
}
/// A grid overlay annotation for reference.
#[derive(Debug, Clone)]
pub struct GridAnnotation {
    /// Origin (bottom-left corner).
    pub origin: [f64; 3],
    /// Width along X.
    pub width: f64,
    /// Height along Y.
    pub height: f64,
    /// Number of divisions along X.
    pub divisions_x: usize,
    /// Number of divisions along Y.
    pub divisions_y: usize,
    /// Whether to show labels at grid intersections.
    pub show_labels: bool,
}
impl GridAnnotation {
    /// Create a new grid annotation.
    pub fn new(
        origin: [f64; 3],
        width: f64,
        height: f64,
        divisions_x: usize,
        divisions_y: usize,
    ) -> Self {
        Self {
            origin,
            width,
            height,
            divisions_x: divisions_x.max(1),
            divisions_y: divisions_y.max(1),
            show_labels: false,
        }
    }
    /// Cell width.
    pub fn cell_width(&self) -> f64 {
        self.width / self.divisions_x as f64
    }
    /// Cell height.
    pub fn cell_height(&self) -> f64 {
        self.height / self.divisions_y as f64
    }
    /// Number of grid intersection points.
    pub fn intersection_count(&self) -> usize {
        (self.divisions_x + 1) * (self.divisions_y + 1)
    }
    /// Number of line segments in the grid.
    pub fn line_count(&self) -> usize {
        (self.divisions_x + 1) + (self.divisions_y + 1)
    }
    /// Generate horizontal line endpoints.
    pub fn horizontal_lines(&self) -> Vec<([f64; 3], [f64; 3])> {
        let mut lines = Vec::with_capacity(self.divisions_y + 1);
        for j in 0..=self.divisions_y {
            let y = self.origin[1] + j as f64 * self.cell_height();
            lines.push((
                [self.origin[0], y, self.origin[2]],
                [self.origin[0] + self.width, y, self.origin[2]],
            ));
        }
        lines
    }
    /// Generate vertical line endpoints.
    pub fn vertical_lines(&self) -> Vec<([f64; 3], [f64; 3])> {
        let mut lines = Vec::with_capacity(self.divisions_x + 1);
        for i in 0..=self.divisions_x {
            let x = self.origin[0] + i as f64 * self.cell_width();
            lines.push((
                [x, self.origin[1], self.origin[2]],
                [x, self.origin[1] + self.height, self.origin[2]],
            ));
        }
        lines
    }
    /// Total number of cells.
    pub fn cell_count(&self) -> usize {
        self.divisions_x * self.divisions_y
    }
}
/// A dimension line indicating a measured distance between two points.
///
/// Consists of two extension lines, a dimension line with arrows, and a text
/// label showing the measured value.
#[derive(Debug, Clone)]
pub struct DimensionLine {
    /// Start point of the dimension.
    pub start: [f64; 3],
    /// End point of the dimension.
    pub end: [f64; 3],
    /// Offset distance for the dimension line from the measured points.
    pub offset: f64,
    /// Optional label text (if empty, auto-generated from length).
    pub label: String,
    /// Arrow head size.
    pub arrow_size: f64,
    /// Extension line overshoot beyond the dimension line.
    pub extension_overshoot: f64,
}
impl DimensionLine {
    /// Create a new dimension line with default styling.
    pub fn new(start: [f64; 3], end: [f64; 3], offset: f64) -> Self {
        Self {
            start,
            end,
            offset,
            label: String::new(),
            arrow_size: 0.1,
            extension_overshoot: 0.05,
        }
    }
    /// Create a dimension line with a custom label.
    pub fn with_label(start: [f64; 3], end: [f64; 3], offset: f64, label: &str) -> Self {
        let mut dim = Self::new(start, end, offset);
        dim.label = label.to_string();
        dim
    }
    /// Measured length between start and end.
    pub fn length(&self) -> f64 {
        dist3(self.start, self.end)
    }
    /// Direction vector from start to end (normalized).
    pub fn direction(&self) -> [f64; 3] {
        normalize3(sub3(self.end, self.start))
    }
    /// Midpoint of the dimension line.
    pub fn midpoint(&self) -> [f64; 3] {
        midpoint3(self.start, self.end)
    }
    /// Label text: returns custom label if set, otherwise formatted length.
    pub fn display_label(&self) -> String {
        if self.label.is_empty() {
            format!("{:.3}", self.length())
        } else {
            self.label.clone()
        }
    }
    /// Compute the offset direction perpendicular to the dimension line.
    ///
    /// Uses the Z-axis as a reference to determine the perpendicular; falls
    /// back to Y-axis if the dimension is parallel to Z.
    pub fn offset_direction(&self) -> [f64; 3] {
        let dir = self.direction();
        let up = if dir[2].abs() > 0.99 {
            [0.0, 1.0, 0.0]
        } else {
            [0.0, 0.0, 1.0]
        };
        normalize3(cross3(dir, up))
    }
    /// Offset start point (where the dimension line actually sits).
    pub fn offset_start(&self) -> [f64; 3] {
        add3(self.start, scale3(self.offset_direction(), self.offset))
    }
    /// Offset end point.
    pub fn offset_end(&self) -> [f64; 3] {
        add3(self.end, scale3(self.offset_direction(), self.offset))
    }
    /// Generate extension line endpoints: returns `[(ext_start_a, ext_end_a), (ext_start_b, ext_end_b)]`.
    pub fn extension_lines(&self) -> [([f64; 3], [f64; 3]); 2] {
        let off_dir = scale3(
            self.offset_direction(),
            self.offset + self.extension_overshoot,
        );
        [
            (self.start, add3(self.start, off_dir)),
            (self.end, add3(self.end, off_dir)),
        ]
    }
    /// Generate arrow tip points at both ends of the dimension line.
    ///
    /// Returns `(tip_start, tip_end)`.
    pub fn arrow_tips(&self) -> ([f64; 3], [f64; 3]) {
        (self.offset_start(), self.offset_end())
    }
    /// Compute the bounding box of the entire dimension annotation in 2D (XY plane).
    pub fn bounding_box_2d(&self) -> BoundingBox2D {
        let os = self.offset_start();
        let oe = self.offset_end();
        let ext = self.extension_lines();
        let mut min_x = self.start[0].min(self.end[0]).min(os[0]).min(oe[0]);
        let mut max_x = self.start[0].max(self.end[0]).max(os[0]).max(oe[0]);
        let mut min_y = self.start[1].min(self.end[1]).min(os[1]).min(oe[1]);
        let mut max_y = self.start[1].max(self.end[1]).max(os[1]).max(oe[1]);
        for (a, b) in &ext {
            min_x = min_x.min(a[0]).min(b[0]);
            max_x = max_x.max(a[0]).max(b[0]);
            min_y = min_y.min(a[1]).min(b[1]);
            max_y = max_y.max(a[1]).max(b[1]);
        }
        BoundingBox2D::new([min_x, min_y], [max_x, max_y])
    }
}
/// Strategy for anchoring an annotation to a feature.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AnchorStrategy {
    /// Fixed position — the annotation does not move.
    Fixed,
    /// Follow the feature point.
    FollowPoint,
    /// Maintain a fixed offset from the feature point.
    FixedOffset,
    /// Snap to the nearest grid intersection.
    SnapToGrid,
}
/// A compass rose annotation for orientation reference.
#[derive(Debug, Clone)]
pub struct CompassRose {
    /// Center position.
    pub center: [f64; 3],
    /// Radius of the compass.
    pub radius: f64,
    /// North direction angle (radians from +Y, CCW).
    pub north_angle: f64,
    /// Number of cardinal directions to show (4, 8, or 16).
    pub cardinal_count: usize,
}
impl CompassRose {
    /// Create a compass rose with north pointing up.
    pub fn new(center: [f64; 3], radius: f64) -> Self {
        Self {
            center,
            radius,
            north_angle: 0.0,
            cardinal_count: 8,
        }
    }
    /// Get the endpoint for a cardinal direction at a given angle.
    pub fn direction_point(&self, angle_rad: f64) -> [f64; 3] {
        let total_angle = angle_rad + self.north_angle;
        [
            self.center[0] + self.radius * total_angle.sin(),
            self.center[1] + self.radius * total_angle.cos(),
            self.center[2],
        ]
    }
    /// Generate line endpoints for all cardinal directions.
    pub fn direction_lines(&self) -> Vec<([f64; 3], [f64; 3])> {
        let n = self.cardinal_count.max(4);
        let mut lines = Vec::with_capacity(n);
        for i in 0..n {
            let angle = 2.0 * PI * i as f64 / n as f64;
            let tip = self.direction_point(angle);
            lines.push((self.center, tip));
        }
        lines
    }
    /// Labels for the cardinal directions (N, NE, E, SE, S, SW, W, NW).
    pub fn direction_labels(&self) -> Vec<(&'static str, [f64; 3])> {
        let labels_8 = ["N", "NE", "E", "SE", "S", "SW", "W", "NW"];
        let n = self.cardinal_count.clamp(4, 8);
        let step = 8 / n;
        let label_r = self.radius * 1.2;
        let mut result = Vec::with_capacity(n);
        for i in 0..n {
            let idx = i * step;
            let angle = 2.0 * PI * i as f64 / n as f64 + self.north_angle;
            let pos = [
                self.center[0] + label_r * angle.sin(),
                self.center[1] + label_r * angle.cos(),
                self.center[2],
            ];
            result.push((labels_8[idx % 8], pos));
        }
        result
    }
}
/// Axis-aligned bounding box in 2D for text layout and collision detection.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BoundingBox2D {
    /// Minimum corner \[x, y\].
    pub min: [f64; 2],
    /// Maximum corner \[x, y\].
    pub max: [f64; 2],
}
impl BoundingBox2D {
    /// Create a new bounding box.
    pub fn new(min: [f64; 2], max: [f64; 2]) -> Self {
        Self { min, max }
    }
    /// Create from center, width, and height.
    pub fn from_center(center: [f64; 2], width: f64, height: f64) -> Self {
        let hw = width * 0.5;
        let hh = height * 0.5;
        Self {
            min: [center[0] - hw, center[1] - hh],
            max: [center[0] + hw, center[1] + hh],
        }
    }
    /// Width of the bounding box.
    pub fn width(&self) -> f64 {
        self.max[0] - self.min[0]
    }
    /// Height of the bounding box.
    pub fn height(&self) -> f64 {
        self.max[1] - self.min[1]
    }
    /// Area of the bounding box.
    pub fn area(&self) -> f64 {
        self.width() * self.height()
    }
    /// Center point.
    pub fn center(&self) -> [f64; 2] {
        [
            (self.min[0] + self.max[0]) * 0.5,
            (self.min[1] + self.max[1]) * 0.5,
        ]
    }
    /// Check if this box overlaps another box.
    pub fn overlaps(&self, other: &BoundingBox2D) -> bool {
        self.min[0] < other.max[0]
            && self.max[0] > other.min[0]
            && self.min[1] < other.max[1]
            && self.max[1] > other.min[1]
    }
    /// Check if a 2D point is inside this box.
    pub fn contains_point(&self, p: [f64; 2]) -> bool {
        p[0] >= self.min[0] && p[0] <= self.max[0] && p[1] >= self.min[1] && p[1] <= self.max[1]
    }
    /// Expand the box by a margin on all sides.
    pub fn expand(&self, margin: f64) -> Self {
        Self {
            min: [self.min[0] - margin, self.min[1] - margin],
            max: [self.max[0] + margin, self.max[1] + margin],
        }
    }
    /// Merge with another bounding box, returning the union.
    pub fn union(&self, other: &BoundingBox2D) -> Self {
        Self {
            min: [self.min[0].min(other.min[0]), self.min[1].min(other.min[1])],
            max: [self.max[0].max(other.max[0]), self.max[1].max(other.max[1])],
        }
    }
    /// Compute the overlap area between two bounding boxes (0 if no overlap).
    pub fn overlap_area(&self, other: &BoundingBox2D) -> f64 {
        let ox = (self.max[0].min(other.max[0]) - self.min[0].max(other.min[0])).max(0.0);
        let oy = (self.max[1].min(other.max[1]) - self.min[1].max(other.min[1])).max(0.0);
        ox * oy
    }
    /// Perimeter of the bounding box.
    pub fn perimeter(&self) -> f64 {
        2.0 * (self.width() + self.height())
    }
}
/// An annotation anchor that ties a label position to a feature position.
#[derive(Debug, Clone)]
pub struct AnnotationAnchor {
    /// The feature point being annotated.
    pub feature_point: [f64; 3],
    /// The current label position.
    pub label_position: [f64; 3],
    /// Anchoring strategy.
    pub strategy: AnchorStrategy,
    /// Fixed offset (used by FixedOffset strategy).
    pub offset: [f64; 3],
    /// Grid spacing (used by SnapToGrid strategy).
    pub grid_spacing: f64,
}
impl AnnotationAnchor {
    /// Create a new annotation anchor.
    pub fn new(
        feature_point: [f64; 3],
        label_position: [f64; 3],
        strategy: AnchorStrategy,
    ) -> Self {
        let offset = sub3(label_position, feature_point);
        Self {
            feature_point,
            label_position,
            strategy,
            offset,
            grid_spacing: 1.0,
        }
    }
    /// Update the label position based on a new feature point location.
    pub fn update(&mut self, new_feature_point: [f64; 3]) {
        self.feature_point = new_feature_point;
        match self.strategy {
            AnchorStrategy::Fixed => {}
            AnchorStrategy::FollowPoint => {
                self.label_position = new_feature_point;
            }
            AnchorStrategy::FixedOffset => {
                self.label_position = add3(new_feature_point, self.offset);
            }
            AnchorStrategy::SnapToGrid => {
                let raw = add3(new_feature_point, self.offset);
                self.label_position = snap_to_grid(raw, self.grid_spacing);
            }
        }
    }
    /// Distance from the feature point to the label.
    pub fn leader_length(&self) -> f64 {
        dist3(self.feature_point, self.label_position)
    }
    /// Direction from feature to label (normalized).
    pub fn leader_direction(&self) -> [f64; 3] {
        normalize3(sub3(self.label_position, self.feature_point))
    }
}
/// A callout box with text content, positioned at a given center.
#[derive(Debug, Clone)]
pub struct CalloutBox {
    /// Center position of the callout box.
    pub center: [f64; 3],
    /// Width of the callout box.
    pub width: f64,
    /// Height of the callout box.
    pub height: f64,
    /// Text content.
    pub text: String,
    /// Corner radius for rounded corners (0 = sharp).
    pub corner_radius: f64,
    /// Padding inside the box.
    pub padding: f64,
    /// Text anchor within the box.
    pub anchor: TextAnchor,
}
impl CalloutBox {
    /// Create a new callout box.
    pub fn new(center: [f64; 3], width: f64, height: f64, text: &str) -> Self {
        Self {
            center,
            width,
            height,
            text: text.to_string(),
            corner_radius: 0.0,
            padding: 0.05,
            anchor: TextAnchor::Center,
        }
    }
    /// Create a callout with rounded corners.
    pub fn rounded(center: [f64; 3], width: f64, height: f64, text: &str, radius: f64) -> Self {
        let mut cb = Self::new(center, width, height, text);
        cb.corner_radius = radius;
        cb
    }
    /// Area of the callout box.
    pub fn area(&self) -> f64 {
        self.width * self.height
    }
    /// Perimeter of the callout box (approximate if rounded).
    pub fn perimeter(&self) -> f64 {
        if self.corner_radius > 0.0 {
            let r = self
                .corner_radius
                .min(self.width * 0.5)
                .min(self.height * 0.5);
            let straight_w = self.width - 2.0 * r;
            let straight_h = self.height - 2.0 * r;
            2.0 * (straight_w + straight_h) + 2.0 * PI * r
        } else {
            2.0 * (self.width + self.height)
        }
    }
    /// Bounding box in 2D.
    pub fn bounding_box(&self) -> BoundingBox2D {
        BoundingBox2D::from_center([self.center[0], self.center[1]], self.width, self.height)
    }
    /// Available text area (after padding).
    pub fn text_area(&self) -> f64 {
        let w = (self.width - 2.0 * self.padding).max(0.0);
        let h = (self.height - 2.0 * self.padding).max(0.0);
        w * h
    }
    /// Get the four corner points of the box in the XY plane.
    pub fn corners(&self) -> [[f64; 3]; 4] {
        let hw = self.width * 0.5;
        let hh = self.height * 0.5;
        let [cx, cy, cz] = self.center;
        [
            [cx - hw, cy - hh, cz],
            [cx + hw, cy - hh, cz],
            [cx + hw, cy + hh, cz],
            [cx - hw, cy + hh, cz],
        ]
    }
    /// Check whether a 2D point (XY) lies inside the callout.
    pub fn contains_point_2d(&self, p: [f64; 2]) -> bool {
        self.bounding_box().contains_point(p)
    }
}
