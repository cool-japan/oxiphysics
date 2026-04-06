//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::too_many_arguments)]
use std::collections::HashMap;

use super::functions::Point3;
#[allow(unused_imports)]
use super::functions::*;

/// A topological noise candidate (a low-persistence feature to be removed).
#[derive(Debug, Clone)]
pub struct NoisePair {
    /// Birth value.
    pub birth: f64,
    /// Death value.
    pub death: f64,
    /// Dimension.
    pub dimension: usize,
    /// Whether this pair is classified as noise (persistence ≤ threshold).
    pub is_noise: bool,
}
impl NoisePair {
    /// Persistence of this pair.
    pub fn persistence(&self) -> f64 {
        (self.death - self.birth).max(0.0)
    }
}
/// Tracks a collection of features over time.
#[derive(Debug, Clone, Default)]
pub struct TopologicalFeatureTracker {
    /// All tracked features indexed by id.
    pub features: HashMap<usize, TrackedFeature>,
    /// Next feature id to assign.
    pub(super) next_id: usize,
}
impl TopologicalFeatureTracker {
    /// Create an empty tracker.
    pub fn new() -> Self {
        Self {
            features: HashMap::new(),
            next_id: 0,
        }
    }
    /// Register a new feature. Returns the assigned id.
    pub fn register_feature(&mut self, dimension: usize) -> usize {
        let id = self.next_id;
        self.next_id += 1;
        self.features.insert(id, TrackedFeature::new(id, dimension));
        id
    }
    /// Record a birth-death pair at a given time for feature `id`.
    pub fn record(&mut self, id: usize, time: f64, pair: BirthDeathPair) {
        if let Some(feat) = self.features.get_mut(&id) {
            feat.record(time, pair);
        }
    }
    /// Number of tracked features.
    pub fn feature_count(&self) -> usize {
        self.features.len()
    }
}
/// A line segment connecting two 3-D points.
#[derive(Debug, Clone)]
pub struct LineSegment {
    /// Start point of the segment.
    pub start: Point3,
    /// End point of the segment.
    pub end: Point3,
    /// Display color.
    pub color: Color4,
    /// Line width hint (in pixels or NDC units).
    pub width: f32,
}
impl LineSegment {
    /// Construct a new line segment.
    pub fn new(start: Point3, end: Point3, color: Color4) -> Self {
        Self {
            start,
            end,
            color,
            width: 1.0,
        }
    }
    /// Construct with explicit width.
    pub fn with_width(start: Point3, end: Point3, color: Color4, width: f32) -> Self {
        Self {
            start,
            end,
            color,
            width,
        }
    }
    /// Euclidean length of the segment.
    pub fn length(&self) -> f64 {
        let dx = self.end[0] - self.start[0];
        let dy = self.end[1] - self.start[1];
        let dz = self.end[2] - self.start[2];
        (dx * dx + dy * dy + dz * dz).sqrt()
    }
}
/// Type of a Morse critical point.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CriticalPointType {
    /// Local minimum (index 0).
    Minimum,
    /// Saddle of index 1.
    Saddle1,
    /// Saddle of index 2.
    Saddle2,
    /// Local maximum (index n for n-dimensional manifold).
    Maximum,
}
impl CriticalPointType {
    /// Return the Morse index.
    pub fn morse_index(&self) -> usize {
        match self {
            Self::Minimum => 0,
            Self::Saddle1 => 1,
            Self::Saddle2 => 2,
            Self::Maximum => 3,
        }
    }
    /// Standard color coding: min=blue, saddle1=cyan, saddle2=orange, max=red.
    pub fn default_color(&self) -> Color4 {
        match self {
            Self::Minimum => Color4::blue(),
            Self::Saddle1 => Color4::cyan(),
            Self::Saddle2 => Color4::orange(),
            Self::Maximum => Color4::red(),
        }
    }
}
/// Visualization for a Jacobi set (critical set of a bivariate map).
#[derive(Debug, Clone, Default)]
pub struct JacobiSetViz {
    /// All Jacobi segments.
    pub segments: Vec<JacobiSegment>,
}
impl JacobiSetViz {
    /// Create an empty Jacobi set visualization.
    pub fn new() -> Self {
        Self {
            segments: Vec::new(),
        }
    }
    /// Add a segment.
    pub fn add_segment(&mut self, seg: JacobiSegment) {
        self.segments.push(seg);
    }
    /// Build all line segments.
    pub fn build_lines(&self) -> Vec<LineSegment> {
        self.segments.iter().map(|s| s.to_line()).collect()
    }
    /// Compute the Jacobi set from a bivariate field on a 2-D grid.
    ///
    /// A sample point lies on the Jacobi set when the 2×2 Jacobian has rank < 2,
    /// detected by |det(J)| < `threshold`.
    pub fn from_bivariate_field_2d(
        field: impl Fn(f64, f64) -> [f64; 2],
        x0: f64,
        x1: f64,
        nx: usize,
        y0: f64,
        y1: f64,
        ny: usize,
        det_threshold: f64,
    ) -> Self {
        let dx = (x1 - x0) / (nx as f64 - 1.0).max(1.0);
        let dy = (y1 - y0) / (ny as f64 - 1.0).max(1.0);
        let mut segments = Vec::new();
        for j in 0..ny.saturating_sub(1) {
            for i in 0..nx.saturating_sub(1) {
                let x = x0 + i as f64 * dx;
                let y = y0 + j as f64 * dy;
                let [f0, g0] = field(x, y);
                let [f1, _] = field(x + dx, y);
                let [_, g1] = field(x, y + dy);
                let dfdx = (f1 - f0) / dx;
                let dgdy = (g1 - g0) / dy;
                let dfdy = 0.0_f64;
                let dgdx = 0.0_f64;
                let det = dfdx * dgdy - dfdy * dgdx;
                if det.abs() < det_threshold {
                    let start = [x, y, 0.0];
                    let end = [x + dx, y + dy, 0.0];
                    segments.push(JacobiSegment::new(start, end, 1));
                }
            }
        }
        Self { segments }
    }
    /// Number of segments.
    pub fn segment_count(&self) -> usize {
        self.segments.len()
    }
}
/// Visualization for a contour tree.
#[derive(Debug, Clone, Default)]
pub struct ContourTreeViz {
    /// Vertex scalar values (indexed by vertex id).
    pub vertices: Vec<f64>,
    /// Vertex positions.
    pub positions: Vec<Point3>,
    /// Edges.
    pub edges: Vec<ContourTreeEdge>,
}
impl ContourTreeViz {
    /// Create an empty contour tree.
    pub fn new() -> Self {
        Self {
            vertices: Vec::new(),
            positions: Vec::new(),
            edges: Vec::new(),
        }
    }
    /// Add a vertex and return its index.
    pub fn add_vertex(&mut self, value: f64, position: Point3) -> usize {
        let idx = self.vertices.len();
        self.vertices.push(value);
        self.positions.push(position);
        idx
    }
    /// Add an edge.
    pub fn add_edge(&mut self, lower: usize, upper: usize, color: Color4) {
        self.edges.push(ContourTreeEdge {
            lower,
            upper,
            color,
        });
    }
    /// Build line segments from edges.
    pub fn build_lines(&self) -> Vec<LineSegment> {
        self.edges
            .iter()
            .filter_map(|e| {
                let p0 = self.positions.get(e.lower)?;
                let p1 = self.positions.get(e.upper)?;
                Some(LineSegment::new(*p0, *p1, e.color))
            })
            .collect()
    }
    /// Number of leaf nodes (degree 1 vertices).
    pub fn leaf_count(&self) -> usize {
        let n = self.vertices.len();
        let mut degree = vec![0usize; n];
        for e in &self.edges {
            if e.lower < n {
                degree[e.lower] += 1;
            }
            if e.upper < n {
                degree[e.upper] += 1;
            }
        }
        degree.iter().filter(|&&d| d == 1).count()
    }
}
/// An edge in the contour tree.
#[derive(Debug, Clone)]
pub struct ContourTreeEdge {
    /// Lower vertex (by scalar value).
    pub lower: usize,
    /// Upper vertex (by scalar value).
    pub upper: usize,
    /// Color for this edge.
    pub color: Color4,
}
/// An arc in a Reeb graph connecting two nodes.
#[derive(Debug, Clone)]
pub struct ReebArc {
    /// Source node id.
    pub source: usize,
    /// Target node id.
    pub target: usize,
    /// Sampled points along the arc for rendering.
    pub polyline: Vec<Point3>,
    /// Display color.
    pub color: Color4,
}
/// Visualization data for a Reeb graph.
#[derive(Debug, Clone, Default)]
pub struct ReebGraphViz {
    /// All nodes.
    pub nodes: Vec<ReebNode>,
    /// All arcs.
    pub arcs: Vec<ReebArc>,
}
impl ReebGraphViz {
    /// Create an empty Reeb graph visualization.
    pub fn new() -> Self {
        Self {
            nodes: Vec::new(),
            arcs: Vec::new(),
        }
    }
    /// Add a node and return its id.
    pub fn add_node(&mut self, value: f64, position: Point3, kind: CriticalPointType) -> usize {
        let id = self.nodes.len();
        self.nodes.push(ReebNode {
            id,
            value,
            position,
            kind,
        });
        id
    }
    /// Add an arc between two nodes.
    pub fn add_arc(&mut self, source: usize, target: usize, polyline: Vec<Point3>, color: Color4) {
        self.arcs.push(ReebArc {
            source,
            target,
            polyline,
            color,
        });
    }
    /// Build line segments for all arcs.
    pub fn build_lines(&self) -> Vec<LineSegment> {
        let mut lines = Vec::new();
        for arc in &self.arcs {
            let pts = &arc.polyline;
            for i in 0..pts.len().saturating_sub(1) {
                lines.push(LineSegment::new(pts[i], pts[i + 1], arc.color));
            }
        }
        lines
    }
    /// Build node glyphs.
    pub fn build_node_glyphs(&self, radius: f64) -> Vec<PointGlyph> {
        self.nodes
            .iter()
            .map(|n| {
                PointGlyph::labeled(
                    n.position,
                    n.kind.default_color(),
                    radius,
                    format!("{:.6}", n.value),
                )
            })
            .collect()
    }
    /// Look up a node by id.
    pub fn node(&self, id: usize) -> Option<&ReebNode> {
        self.nodes.get(id)
    }
}
/// Visualization of Euler characteristic over a sequence of threshold values.
#[derive(Debug, Clone, Default)]
pub struct EulerCharacteristicViz {
    /// Threshold values (x-axis of the curve).
    pub thresholds: Vec<f64>,
    /// Euler characteristic at each threshold.
    pub chi_values: Vec<i64>,
}
impl EulerCharacteristicViz {
    /// Construct from parallel threshold/chi vectors.
    pub fn new(thresholds: Vec<f64>, chi_values: Vec<i64>) -> Self {
        assert_eq!(thresholds.len(), chi_values.len(), "lengths must match");
        Self {
            thresholds,
            chi_values,
        }
    }
    /// Build line segments tracing the Euler characteristic curve in 2-D.
    ///
    /// The threshold maps to world X, chi maps to world Y, Z is fixed at 0.
    /// The curve is scaled so that `scale_y` world units span the full chi range.
    pub fn build_curve_lines(&self, scale_y: f64) -> Vec<LineSegment> {
        if self.thresholds.len() < 2 {
            return Vec::new();
        }
        let chi_min = *self.chi_values.iter().min().unwrap_or(&0) as f64;
        let chi_max = *self.chi_values.iter().max().unwrap_or(&0) as f64;
        let chi_range = (chi_max - chi_min).max(1.0);
        let t_min = self
            .thresholds
            .iter()
            .cloned()
            .fold(f64::INFINITY, f64::min);
        let t_max = self
            .thresholds
            .iter()
            .cloned()
            .fold(f64::NEG_INFINITY, f64::max);
        let t_range = (t_max - t_min).max(1e-14);
        let mut lines = Vec::new();
        let n = self.thresholds.len();
        for i in 0..n - 1 {
            let x0 = (self.thresholds[i] - t_min) / t_range;
            let y0 = (self.chi_values[i] as f64 - chi_min) / chi_range * scale_y;
            let x1 = (self.thresholds[i + 1] - t_min) / t_range;
            let y1 = (self.chi_values[i + 1] as f64 - chi_min) / chi_range * scale_y;
            let color = if self.chi_values[i] >= 0 {
                Color4::green()
            } else {
                Color4::red()
            };
            lines.push(LineSegment::new([x0, y0, 0.0], [x1, y1, 0.0], color));
        }
        lines
    }
    /// Return the minimum Euler characteristic.
    pub fn chi_min(&self) -> Option<i64> {
        self.chi_values.iter().cloned().min()
    }
    /// Return the maximum Euler characteristic.
    pub fn chi_max(&self) -> Option<i64> {
        self.chi_values.iter().cloned().max()
    }
}
/// A horizontal bar in a persistence barcode.
#[derive(Debug, Clone)]
pub struct BarcodeBar {
    /// Birth value (left endpoint).
    pub birth: f64,
    /// Death value (right endpoint; may be infinity for essential features).
    pub death: f64,
    /// Vertical position (lane index).
    pub lane: usize,
    /// Homology dimension.
    pub dimension: usize,
}
impl BarcodeBar {
    /// Build a [`LineSegment`] for this bar.
    ///
    /// `max_death` clamps infinite bars.
    pub fn to_line(&self, max_death: f64, lane_height: f64) -> LineSegment {
        let d = if self.death.is_infinite() {
            max_death
        } else {
            self.death
        };
        let y = self.lane as f64 * lane_height;
        let color = match self.dimension {
            0 => Color4::blue(),
            1 => Color4::green(),
            2 => Color4::red(),
            _ => Color4::white(),
        };
        LineSegment::with_width([self.birth, y, 0.0], [d, y, 0.0], color, 2.0)
    }
}
/// A single birth/death pair from a persistence diagram.
#[derive(Debug, Clone, PartialEq)]
pub struct BirthDeathPair {
    /// Filtration value at which the feature is born.
    pub birth: f64,
    /// Filtration value at which the feature dies (`f64::INFINITY` if essential).
    pub death: f64,
    /// Homology dimension.
    pub dimension: usize,
}
impl BirthDeathPair {
    /// Create a new birth-death pair.
    pub fn new(birth: f64, death: f64, dimension: usize) -> Self {
        Self {
            birth,
            death,
            dimension,
        }
    }
    /// Persistence = death − birth.
    pub fn persistence(&self) -> f64 {
        if self.death.is_infinite() {
            f64::INFINITY
        } else {
            self.death - self.birth
        }
    }
    /// Whether this pair represents an essential class.
    pub fn is_essential(&self) -> bool {
        self.death.is_infinite()
    }
    /// Return the position of this pair in the birth-death plane as `[birth, death]`.
    pub fn diagram_point(&self) -> [f64; 2] {
        [self.birth, self.death]
    }
}
/// Visualization of all unstable manifolds in a scalar field.
#[derive(Debug, Clone, Default)]
pub struct UnstableManifoldViz {
    /// All arcs.
    pub arcs: Vec<UnstableManifoldArc>,
    /// Color for unstable manifold lines.
    pub color: Color4,
}
impl UnstableManifoldViz {
    /// Create an empty visualization with a given color.
    pub fn new(color: Color4) -> Self {
        Self {
            arcs: Vec::new(),
            color,
        }
    }
    /// Add an arc.
    pub fn add_arc(&mut self, arc: UnstableManifoldArc) {
        self.arcs.push(arc);
    }
    /// Build all arc lines.
    pub fn build_lines(&self) -> Vec<LineSegment> {
        self.arcs
            .iter()
            .flat_map(|a| a.to_lines(self.color))
            .collect()
    }
    /// Total number of sample points across all arcs.
    pub fn total_points(&self) -> usize {
        self.arcs.iter().map(|a| a.points.len()).sum()
    }
    /// Trace a single gradient ascent arc from a seed point.
    ///
    /// Uses a fixed step-size Euler integration of the gradient field.
    pub fn trace_arc(
        gradient: impl Fn(f64, f64, f64) -> [f64; 3],
        seed: Point3,
        step_size: f64,
        max_steps: usize,
    ) -> Vec<Point3> {
        let mut pts = vec![seed];
        let mut pos = seed;
        for _ in 0..max_steps {
            let g = gradient(pos[0], pos[1], pos[2]);
            let norm = (g[0] * g[0] + g[1] * g[1] + g[2] * g[2]).sqrt();
            if norm < 1e-14 {
                break;
            }
            let dir = [g[0] / norm, g[1] / norm, g[2] / norm];
            pos = [
                pos[0] + dir[0] * step_size,
                pos[1] + dir[1] * step_size,
                pos[2] + dir[2] * step_size,
            ];
            pts.push(pos);
        }
        pts
    }
}
/// Visualization for a persistence barcode.
#[derive(Debug, Clone, Default)]
pub struct PersistenceBarcodeViz {
    /// All bars.
    pub bars: Vec<BarcodeBar>,
    /// Height between barcode lanes.
    pub lane_height: f64,
}
impl PersistenceBarcodeViz {
    /// Construct from birth-death pairs, assigning lanes sequentially.
    pub fn from_pairs(pairs: &[BirthDeathPair], lane_height: f64) -> Self {
        let bars = pairs
            .iter()
            .enumerate()
            .map(|(i, p)| BarcodeBar {
                birth: p.birth,
                death: p.death,
                lane: i,
                dimension: p.dimension,
            })
            .collect();
        Self { bars, lane_height }
    }
    /// Build all bar line segments.
    pub fn build_lines(&self, max_death: f64) -> Vec<LineSegment> {
        self.bars
            .iter()
            .map(|b| b.to_line(max_death, self.lane_height))
            .collect()
    }
    /// Number of bars.
    pub fn bar_count(&self) -> usize {
        self.bars.len()
    }
    /// Number of essential (infinite) bars.
    pub fn essential_count(&self) -> usize {
        self.bars.iter().filter(|b| b.death.is_infinite()).count()
    }
}
/// Visualization of a gradient vector field on a regular grid.
#[derive(Debug, Clone, Default)]
pub struct GradientVectorFieldViz {
    /// All gradient arrows.
    pub arrows: Vec<GradientArrow>,
    /// Scale factor applied to arrow lengths for display.
    pub arrow_scale: f64,
}
impl GradientVectorFieldViz {
    /// Create an empty gradient field visualization.
    pub fn new(arrow_scale: f64) -> Self {
        Self {
            arrows: Vec::new(),
            arrow_scale,
        }
    }
    /// Sample a scalar field on a regular 2-D grid and compute finite-difference gradients.
    ///
    /// The grid spans `[x0, x1] × [y0, y1]` with `nx × ny` samples.
    /// `field` maps `(x, y)` → scalar value.
    pub fn from_scalar_field_2d(
        field: impl Fn(f64, f64) -> f64,
        x0: f64,
        x1: f64,
        nx: usize,
        y0: f64,
        y1: f64,
        ny: usize,
        colormap: impl Fn(f64) -> Color4,
    ) -> Self {
        let dx = (x1 - x0) / (nx as f64 - 1.0).max(1.0);
        let dy = (y1 - y0) / (ny as f64 - 1.0).max(1.0);
        let mut arrows = Vec::with_capacity(nx * ny);
        for j in 0..ny {
            for i in 0..nx {
                let x = x0 + i as f64 * dx;
                let y = y0 + j as f64 * dy;
                let fx_p = field(x + 0.5 * dx, y);
                let fx_m = field(x - 0.5 * dx, y);
                let fy_p = field(x, y + 0.5 * dy);
                let fy_m = field(x, y - 0.5 * dy);
                let gx = (fx_p - fx_m) / dx;
                let gy = (fy_p - fy_m) / dy;
                let mag = (gx * gx + gy * gy).sqrt();
                let color = colormap(mag);
                arrows.push(GradientArrow {
                    origin: [x, y, 0.0],
                    direction: [gx, gy, 0.0],
                    magnitude: mag,
                    color,
                });
            }
        }
        Self {
            arrows,
            arrow_scale: 1.0,
        }
    }
    /// Build all arrows as line segments.
    pub fn build_lines(&self) -> Vec<LineSegment> {
        self.arrows
            .iter()
            .map(|a| a.to_line(self.arrow_scale))
            .collect()
    }
    /// Maximum gradient magnitude across all arrows.
    pub fn max_magnitude(&self) -> f64 {
        self.arrows
            .iter()
            .map(|a| a.magnitude)
            .fold(0.0_f64, f64::max)
    }
}
/// A pair cancellation step in topological simplification.
#[derive(Debug, Clone)]
pub struct CancellationStep {
    /// Index of the first (lower) critical point being cancelled.
    pub cp_index_a: usize,
    /// Index of the second (higher) critical point being cancelled.
    pub cp_index_b: usize,
    /// Persistence of this cancellation.
    pub persistence: f64,
    /// Whether this step has been applied.
    pub applied: bool,
}
/// A critical point in a scalar field (Morse theory).
#[derive(Debug, Clone)]
pub struct CriticalPoint {
    /// World-space position of the critical point.
    pub position: Point3,
    /// Scalar field value at this point.
    pub value: f64,
    /// Type of the critical point.
    pub kind: CriticalPointType,
    /// Index within the parent dataset.
    pub index: usize,
}
impl CriticalPoint {
    /// Construct a new critical point.
    pub fn new(position: Point3, value: f64, kind: CriticalPointType, index: usize) -> Self {
        Self {
            position,
            value,
            kind,
            index,
        }
    }
    /// Build a [`PointGlyph`] for this critical point.
    pub fn to_glyph(&self, radius: f64) -> PointGlyph {
        PointGlyph::labeled(
            self.position,
            self.kind.default_color(),
            radius,
            format!("{:.6}", self.value),
        )
    }
}
/// A topological feature tracked over multiple time steps.
#[derive(Debug, Clone)]
pub struct TrackedFeature {
    /// Unique feature identifier (stable across time steps).
    pub id: usize,
    /// Homology dimension.
    pub dimension: usize,
    /// Birth-death pairs indexed by time step.
    pub time_series: Vec<(f64, BirthDeathPair)>,
}
impl TrackedFeature {
    /// Create a new tracked feature.
    pub fn new(id: usize, dimension: usize) -> Self {
        Self {
            id,
            dimension,
            time_series: Vec::new(),
        }
    }
    /// Record the feature's birth-death pair at a given time.
    pub fn record(&mut self, time: f64, pair: BirthDeathPair) {
        self.time_series.push((time, pair));
    }
    /// Mean persistence over all recorded steps.
    pub fn mean_persistence(&self) -> f64 {
        if self.time_series.is_empty() {
            return 0.0;
        }
        let sum: f64 = self
            .time_series
            .iter()
            .map(|(_, p)| {
                if p.persistence().is_finite() {
                    p.persistence()
                } else {
                    0.0
                }
            })
            .sum();
        sum / self.time_series.len() as f64
    }
    /// Build a line tracing mean persistence over time in 2-D.
    pub fn build_persistence_trace(&self) -> Vec<LineSegment> {
        if self.time_series.len() < 2 {
            return Vec::new();
        }
        let color = match self.dimension {
            0 => Color4::blue(),
            1 => Color4::green(),
            2 => Color4::red(),
            _ => Color4::white(),
        };
        let mut lines = Vec::new();
        for i in 0..self.time_series.len() - 1 {
            let (t0, p0) = &self.time_series[i];
            let (t1, p1) = &self.time_series[i + 1];
            let y0 = if p0.persistence().is_finite() {
                p0.persistence()
            } else {
                0.0
            };
            let y1 = if p1.persistence().is_finite() {
                p1.persistence()
            } else {
                0.0
            };
            lines.push(LineSegment::new([*t0, y0, 0.0], [*t1, y1, 0.0], color));
        }
        lines
    }
}
/// A segment of the Jacobi set between two scalar fields.
#[derive(Debug, Clone)]
pub struct JacobiSegment {
    /// Start position.
    pub start: Point3,
    /// End position.
    pub end: Point3,
    /// Jacobi set segment type:
    /// 0 = regular, 1 = fold (rank drops), 2 = cusp.
    pub segment_type: u8,
}
impl JacobiSegment {
    /// Create a new Jacobi segment.
    pub fn new(start: Point3, end: Point3, segment_type: u8) -> Self {
        Self {
            start,
            end,
            segment_type,
        }
    }
    /// Build a [`LineSegment`] for rendering.
    pub fn to_line(&self) -> LineSegment {
        let color = match self.segment_type {
            0 => Color4::cyan(),
            1 => Color4::magenta(),
            _ => Color4::yellow(),
        };
        LineSegment::new(self.start, self.end, color)
    }
}
/// A node in a Reeb graph.
#[derive(Debug, Clone)]
pub struct ReebNode {
    /// Node identifier.
    pub id: usize,
    /// Scalar value at this node.
    pub value: f64,
    /// World-space position for rendering.
    pub position: Point3,
    /// Type of the associated critical point.
    pub kind: CriticalPointType,
}
/// A simple RGBA color with f32 components in \[0, 1\].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Color4 {
    /// Red channel.
    pub r: f32,
    /// Green channel.
    pub g: f32,
    /// Blue channel.
    pub b: f32,
    /// Alpha channel.
    pub a: f32,
}
impl Color4 {
    /// Create a new color.
    pub fn new(r: f32, g: f32, b: f32, a: f32) -> Self {
        Self { r, g, b, a }
    }
    /// Fully opaque red.
    pub fn red() -> Self {
        Self::new(1.0, 0.0, 0.0, 1.0)
    }
    /// Fully opaque green.
    pub fn green() -> Self {
        Self::new(0.0, 1.0, 0.0, 1.0)
    }
    /// Fully opaque blue.
    pub fn blue() -> Self {
        Self::new(0.0, 0.0, 1.0, 1.0)
    }
    /// Fully opaque white.
    pub fn white() -> Self {
        Self::new(1.0, 1.0, 1.0, 1.0)
    }
    /// Fully opaque yellow.
    pub fn yellow() -> Self {
        Self::new(1.0, 1.0, 0.0, 1.0)
    }
    /// Fully opaque cyan.
    pub fn cyan() -> Self {
        Self::new(0.0, 1.0, 1.0, 1.0)
    }
    /// Fully opaque magenta.
    pub fn magenta() -> Self {
        Self::new(1.0, 0.0, 1.0, 1.0)
    }
    /// Fully opaque orange.
    pub fn orange() -> Self {
        Self::new(1.0, 0.5, 0.0, 1.0)
    }
    /// Linear interpolation between two colors.
    pub fn lerp(a: Self, b: Self, t: f32) -> Self {
        let t = t.clamp(0.0, 1.0);
        Self::new(
            a.r + (b.r - a.r) * t,
            a.g + (b.g - a.g) * t,
            a.b + (b.b - a.b) * t,
            a.a + (b.a - a.a) * t,
        )
    }
}
/// A colored point glyph for rendering a feature at a 3-D location.
#[derive(Debug, Clone)]
pub struct PointGlyph {
    /// World-space position.
    pub position: Point3,
    /// Display color.
    pub color: Color4,
    /// Glyph radius in world units.
    pub radius: f64,
    /// Optional label string.
    pub label: Option<String>,
}
impl PointGlyph {
    /// Create a point glyph at the given position.
    pub fn new(position: Point3, color: Color4, radius: f64) -> Self {
        Self {
            position,
            color,
            radius,
            label: None,
        }
    }
    /// Create a labeled glyph.
    pub fn labeled(position: Point3, color: Color4, radius: f64, label: impl Into<String>) -> Self {
        Self {
            position,
            color,
            radius,
            label: Some(label.into()),
        }
    }
}
/// Visualization for topological noise removal via persistence thresholding.
#[derive(Debug, Clone, Default)]
pub struct TopologicalNoiseRemovalViz {
    /// All birth-death pairs.
    pub pairs: Vec<NoisePair>,
    /// Current noise threshold.
    pub threshold: f64,
}
impl TopologicalNoiseRemovalViz {
    /// Construct with a set of pairs and an initial threshold.
    pub fn new(pairs: Vec<BirthDeathPair>, threshold: f64) -> Self {
        let noise_pairs = pairs
            .into_iter()
            .map(|p| {
                let persistence = if p.death.is_infinite() {
                    f64::INFINITY
                } else {
                    p.death - p.birth
                };
                let is_noise = persistence.is_finite() && persistence <= threshold;
                NoisePair {
                    birth: p.birth,
                    death: if p.death.is_infinite() {
                        p.birth + 1e6
                    } else {
                        p.death
                    },
                    dimension: p.dimension,
                    is_noise,
                }
            })
            .collect();
        Self {
            pairs: noise_pairs,
            threshold,
        }
    }
    /// Apply a new threshold, updating `is_noise` flags.
    pub fn update_threshold(&mut self, threshold: f64) {
        self.threshold = threshold;
        for p in &mut self.pairs {
            p.is_noise = p.persistence() <= threshold;
        }
    }
    /// Count how many pairs are classified as noise.
    pub fn noise_count(&self) -> usize {
        self.pairs.iter().filter(|p| p.is_noise).count()
    }
    /// Count how many pairs are considered significant.
    pub fn signal_count(&self) -> usize {
        self.pairs.iter().filter(|p| !p.is_noise).count()
    }
    /// Build glyphs for noise pairs (shown in dim grey) and signal pairs (colored by dimension).
    pub fn build_diagram_glyphs(&self, max_val: f64) -> Vec<PointGlyph> {
        self.pairs
            .iter()
            .map(|p| {
                let death_vis = p.death.min(max_val);
                let pos = [p.birth, death_vis, 0.0];
                let color = if p.is_noise {
                    Color4::new(0.4, 0.4, 0.4, 0.5)
                } else {
                    match p.dimension {
                        0 => Color4::blue(),
                        1 => Color4::green(),
                        2 => Color4::red(),
                        _ => Color4::white(),
                    }
                };
                PointGlyph::new(pos, color, if p.is_noise { 0.01 } else { 0.02 })
            })
            .collect()
    }
}
/// Visualization geometry for a persistence diagram.
///
/// Produces scatter points and the diagonal line for the birth-death plane.
#[derive(Debug, Clone, Default)]
pub struct PersistenceDiagramViz {
    /// The persistence pairs.
    pub pairs: Vec<BirthDeathPair>,
    /// Color assigned to dimension 0 pairs.
    pub color_dim0: Color4,
    /// Color assigned to dimension 1 pairs.
    pub color_dim1: Color4,
    /// Color assigned to dimension 2 pairs.
    pub color_dim2: Color4,
    /// Point radius in diagram space.
    pub point_radius: f64,
    /// Whether to draw the diagonal line.
    pub show_diagonal: bool,
}
impl PersistenceDiagramViz {
    /// Construct a new visualization with default colors.
    pub fn new(pairs: Vec<BirthDeathPair>) -> Self {
        Self {
            pairs,
            color_dim0: Color4::blue(),
            color_dim1: Color4::green(),
            color_dim2: Color4::red(),
            point_radius: 0.02,
            show_diagonal: true,
        }
    }
    /// Get the color for a given dimension.
    pub fn color_for_dim(&self, dim: usize) -> Color4 {
        match dim {
            0 => self.color_dim0,
            1 => self.color_dim1,
            2 => self.color_dim2,
            _ => Color4::white(),
        }
    }
    /// Build diagram glyphs (one per pair) in a 2-D layout embedded into 3-D.
    ///
    /// The birth axis maps to world X and the death axis maps to world Y.
    /// Infinite death values are clamped to `max_death`.
    pub fn build_glyphs(&self, max_death: f64) -> Vec<PointGlyph> {
        self.pairs
            .iter()
            .map(|p| {
                let death_clamped = if p.death.is_infinite() {
                    max_death
                } else {
                    p.death
                };
                let pos = [p.birth, death_clamped, 0.0];
                PointGlyph::new(pos, self.color_for_dim(p.dimension), self.point_radius)
            })
            .collect()
    }
    /// Build the diagonal line segment from (min_val, min_val) to (max_val, max_val).
    pub fn build_diagonal(&self, min_val: f64, max_val: f64) -> Option<LineSegment> {
        if !self.show_diagonal {
            return None;
        }
        Some(LineSegment::new(
            [min_val, min_val, 0.0],
            [max_val, max_val, 0.0],
            Color4::new(0.5, 0.5, 0.5, 0.8),
        ))
    }
    /// Compute the bottleneck distance between this diagram and `other`
    /// (approximation: max matched persistence difference).
    pub fn bottleneck_distance_approx(&self, other: &Self) -> f64 {
        let mut dist = 0.0_f64;
        for dim in 0..3 {
            let mut a: Vec<&BirthDeathPair> =
                self.pairs.iter().filter(|p| p.dimension == dim).collect();
            let mut b: Vec<&BirthDeathPair> =
                other.pairs.iter().filter(|p| p.dimension == dim).collect();
            a.sort_by(|x, y| {
                x.birth
                    .partial_cmp(&y.birth)
                    .unwrap_or(std::cmp::Ordering::Equal)
            });
            b.sort_by(|x, y| {
                x.birth
                    .partial_cmp(&y.birth)
                    .unwrap_or(std::cmp::Ordering::Equal)
            });
            for (pa, pb) in a.iter().zip(b.iter()) {
                let dp = (pa.persistence() - pb.persistence()).abs();
                if dp > dist {
                    dist = dp;
                }
            }
        }
        dist
    }
}
/// A fiber surface is a level set in a bivariate scalar field.
/// This visualization stores sampled contour points on the fiber.
#[derive(Debug, Clone, Default)]
pub struct FiberSurfaceViz {
    /// Sample points on the fiber surface.
    pub points: Vec<Point3>,
    /// Color of the fiber surface.
    pub color: Color4,
    /// Target value pair `(f, g)` defining the fiber.
    pub target: [f64; 2],
}
impl FiberSurfaceViz {
    /// Construct a new fiber surface visualization.
    pub fn new(target: [f64; 2], color: Color4) -> Self {
        Self {
            points: Vec::new(),
            color,
            target,
        }
    }
    /// Sample a fiber surface from a bivariate field on a 3-D grid.
    ///
    /// `field` maps `(x, y, z)` to `[f, g]`. The fiber is where `f ≈ target[0]`
    /// AND `g ≈ target[1]`. Points within `tolerance` of both targets are kept.
    pub fn sample_from_field(
        field: impl Fn(f64, f64, f64) -> [f64; 2],
        bounds: [[f64; 2]; 3],
        steps: [usize; 3],
        target: [f64; 2],
        tolerance: f64,
        color: Color4,
    ) -> Self {
        let [bx, by, bz] = bounds;
        let [nx, ny, nz] = steps;
        let dx = (bx[1] - bx[0]) / (nx as f64 - 1.0).max(1.0);
        let dy = (by[1] - by[0]) / (ny as f64 - 1.0).max(1.0);
        let dz = (bz[1] - bz[0]) / (nz as f64 - 1.0).max(1.0);
        let mut points = Vec::new();
        for k in 0..nz {
            for j in 0..ny {
                for i in 0..nx {
                    let x = bx[0] + i as f64 * dx;
                    let y = by[0] + j as f64 * dy;
                    let z = bz[0] + k as f64 * dz;
                    let [f, g] = field(x, y, z);
                    if (f - target[0]).abs() < tolerance && (g - target[1]).abs() < tolerance {
                        points.push([x, y, z]);
                    }
                }
            }
        }
        Self {
            points,
            color,
            target,
        }
    }
    /// Number of sample points.
    pub fn point_count(&self) -> usize {
        self.points.len()
    }
}
/// A single unstable manifold arc (gradient flow line from saddle to max).
#[derive(Debug, Clone)]
pub struct UnstableManifoldArc {
    /// Sampled points along the arc.
    pub points: Vec<Point3>,
    /// Source saddle index.
    pub saddle_index: usize,
    /// Target maximum index.
    pub max_index: usize,
}
impl UnstableManifoldArc {
    /// Build line segments for this arc.
    pub fn to_lines(&self, color: Color4) -> Vec<LineSegment> {
        let mut lines = Vec::new();
        for i in 0..self.points.len().saturating_sub(1) {
            lines.push(LineSegment::new(self.points[i], self.points[i + 1], color));
        }
        lines
    }
}
/// Visualization data for a Morse complex.
#[derive(Debug, Clone)]
pub struct MorseComplexViz {
    /// All critical points.
    pub critical_points: Vec<CriticalPoint>,
    /// Integral curves connecting saddles to minima (ascending manifolds).
    pub ascending_arcs: Vec<LineSegment>,
    /// Integral curves connecting saddles to maxima (descending manifolds).
    pub descending_arcs: Vec<LineSegment>,
}
impl MorseComplexViz {
    /// Create an empty Morse complex visualization.
    pub fn new() -> Self {
        Self {
            critical_points: Vec::new(),
            ascending_arcs: Vec::new(),
            descending_arcs: Vec::new(),
        }
    }
    /// Add a critical point.
    pub fn add_critical_point(&mut self, cp: CriticalPoint) {
        self.critical_points.push(cp);
    }
    /// Add an ascending arc (saddle → minimum).
    pub fn add_ascending_arc(&mut self, arc: LineSegment) {
        self.ascending_arcs.push(arc);
    }
    /// Add a descending arc (saddle → maximum).
    pub fn add_descending_arc(&mut self, arc: LineSegment) {
        self.descending_arcs.push(arc);
    }
    /// Count critical points of the given type.
    pub fn count_of_type(&self, kind: CriticalPointType) -> usize {
        self.critical_points
            .iter()
            .filter(|cp| cp.kind == kind)
            .count()
    }
    /// Compute the Euler characteristic: χ = #min − #saddle1 + #saddle2 − #max (3-D Morse).
    pub fn euler_characteristic(&self) -> i64 {
        let mins = self.count_of_type(CriticalPointType::Minimum) as i64;
        let s1 = self.count_of_type(CriticalPointType::Saddle1) as i64;
        let s2 = self.count_of_type(CriticalPointType::Saddle2) as i64;
        let maxs = self.count_of_type(CriticalPointType::Maximum) as i64;
        mins - s1 + s2 - maxs
    }
    /// Build all glyphs for rendering.
    pub fn build_glyphs(&self, radius: f64) -> Vec<PointGlyph> {
        render_critical_points(&self.critical_points, radius)
    }
    /// Collect all arcs into one flat list (ascending first, then descending).
    pub fn all_arcs(&self) -> Vec<&LineSegment> {
        self.ascending_arcs
            .iter()
            .chain(self.descending_arcs.iter())
            .collect()
    }
}
/// A single gradient arrow glyph at a sample point.
#[derive(Debug, Clone)]
pub struct GradientArrow {
    /// Origin of the gradient vector.
    pub origin: Point3,
    /// Gradient direction (not necessarily normalized).
    pub direction: Point3,
    /// Magnitude of the gradient.
    pub magnitude: f64,
    /// Display color.
    pub color: Color4,
}
impl GradientArrow {
    /// Build a [`LineSegment`] representing this arrow scaled by `scale`.
    pub fn to_line(&self, scale: f64) -> LineSegment {
        let end = [
            self.origin[0] + self.direction[0] * scale,
            self.origin[1] + self.direction[1] * scale,
            self.origin[2] + self.direction[2] * scale,
        ];
        LineSegment::new(self.origin, end, self.color)
    }
}
/// Visualization of topology simplification by pair cancellation.
#[derive(Debug, Clone, Default)]
pub struct TopologySimplificationViz {
    /// Original critical points.
    pub original_cps: Vec<CriticalPoint>,
    /// Planned cancellations sorted by persistence.
    pub cancellations: Vec<CancellationStep>,
    /// Current simplification threshold.
    pub threshold: f64,
}
impl TopologySimplificationViz {
    /// Construct a new simplification visualization.
    pub fn new(original_cps: Vec<CriticalPoint>) -> Self {
        Self {
            original_cps,
            cancellations: Vec::new(),
            threshold: 0.0,
        }
    }
    /// Add a cancellation pair.
    pub fn add_cancellation(&mut self, cp_a: usize, cp_b: usize, persistence: f64) {
        self.cancellations.push(CancellationStep {
            cp_index_a: cp_a,
            cp_index_b: cp_b,
            persistence,
            applied: false,
        });
        self.cancellations.sort_by(|a, b| {
            a.persistence
                .partial_cmp(&b.persistence)
                .unwrap_or(std::cmp::Ordering::Equal)
        });
    }
    /// Apply all cancellations with persistence ≤ threshold.
    pub fn apply_threshold(&mut self, threshold: f64) {
        self.threshold = threshold;
        for c in &mut self.cancellations {
            c.applied = c.persistence <= threshold;
        }
    }
    /// Return the indices of active (non-cancelled) critical points.
    pub fn active_indices(&self) -> Vec<usize> {
        let mut cancelled = std::collections::HashSet::new();
        for c in &self.cancellations {
            if c.applied {
                cancelled.insert(c.cp_index_a);
                cancelled.insert(c.cp_index_b);
            }
        }
        (0..self.original_cps.len())
            .filter(|i| !cancelled.contains(i))
            .collect()
    }
    /// Build glyphs for active critical points only.
    pub fn build_active_glyphs(&self, radius: f64) -> Vec<PointGlyph> {
        self.active_indices()
            .iter()
            .map(|&i| self.original_cps[i].to_glyph(radius))
            .collect()
    }
    /// Number of cancelled pairs.
    pub fn cancelled_count(&self) -> usize {
        self.cancellations.iter().filter(|c| c.applied).count()
    }
}
