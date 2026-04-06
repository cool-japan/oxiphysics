// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Scientific plotting data structures and SVG output.
//!
//! Provides composable building blocks for 2-D scientific plots: axes,
//! series, layouts, and an SVG renderer.  All types are pure data — no
//! windowing or GPU dependency is required.

#![allow(dead_code)]
#![allow(clippy::too_many_arguments)]

// ---------------------------------------------------------------------------
// ColorScheme
// ---------------------------------------------------------------------------

/// Named color scheme for rendering plot elements.
#[derive(Debug, Clone, PartialEq)]
pub enum ColorScheme {
    /// Black-and-white grayscale.
    Grayscale,
    /// Viridis perceptual colormap (purple → green → yellow).
    Viridis,
    /// Jet colormap (blue → cyan → green → yellow → red).
    Jet,
    /// Custom RGB color given as `(red, green, blue)` in `[0, 1]`.
    Custom(f64, f64, f64),
}

impl ColorScheme {
    /// Return an SVG-compatible `rgb(r,g,b)` string for the scheme.
    ///
    /// For named colormaps, returns the color at the midpoint of the map.
    pub fn to_svg_color(&self) -> String {
        match self {
            ColorScheme::Grayscale => "rgb(128,128,128)".into(),
            ColorScheme::Viridis => "rgb(68,1,84)".into(),
            ColorScheme::Jet => "rgb(0,0,255)".into(),
            ColorScheme::Custom(r, g, b) => {
                let ri = (r.clamp(0.0, 1.0) * 255.0) as u8;
                let gi = (g.clamp(0.0, 1.0) * 255.0) as u8;
                let bi = (b.clamp(0.0, 1.0) * 255.0) as u8;
                format!("rgb({ri},{gi},{bi})")
            }
        }
    }

    /// Return a suggested color string for a series index.
    ///
    /// Cycles through a short palette of distinct colors.
    pub fn palette_color(idx: usize) -> String {
        const PALETTE: &[&str] = &[
            "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f",
            "#bcbd22", "#17becf",
        ];
        PALETTE[idx % PALETTE.len()].to_string()
    }
}

// ---------------------------------------------------------------------------
// PlotType
// ---------------------------------------------------------------------------

/// The visual representation type for a data series.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum PlotType {
    /// Connected line plot.
    Line,
    /// Scatter plot (unconnected markers).
    Scatter,
    /// Vertical bar chart.
    Bar,
    /// Frequency histogram.
    Histogram,
    /// Contour / isolines plot.
    Contour,
    /// 3-D surface (rendered as a flat projection in SVG).
    Surface,
    /// Vector / quiver field.
    VectorField,
    /// Line plot with error bars.
    ErrorBar,
}

// ---------------------------------------------------------------------------
// LineStyle
// ---------------------------------------------------------------------------

/// Stroke style for line primitives.
#[derive(Debug, Clone, PartialEq)]
pub enum LineStyle {
    /// Solid line.
    Solid,
    /// Dashed line.
    Dashed,
    /// Dotted line.
    Dotted,
    /// Dash-dot pattern.
    DashDot,
}

impl LineStyle {
    /// Return the SVG `stroke-dasharray` attribute value for this style.
    pub fn svg_dash_array(&self) -> Option<&'static str> {
        match self {
            LineStyle::Solid => None,
            LineStyle::Dashed => Some("8,4"),
            LineStyle::Dotted => Some("2,4"),
            LineStyle::DashDot => Some("8,4,2,4"),
        }
    }
}

// ---------------------------------------------------------------------------
// MarkerShape
// ---------------------------------------------------------------------------

/// Shape used for data-point markers.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum MarkerShape {
    /// No marker.
    None,
    /// Circle marker.
    Circle,
    /// Square marker.
    Square,
    /// Triangle (pointing up).
    Triangle,
    /// Diamond marker.
    Diamond,
    /// Cross / plus marker.
    Cross,
}

// ---------------------------------------------------------------------------
// DataPoint2D
// ---------------------------------------------------------------------------

/// A 2-D data point with an optional error bar.
#[derive(Debug, Clone, PartialEq)]
pub struct DataPoint2D {
    /// Horizontal coordinate.
    pub x: f64,
    /// Vertical coordinate.
    pub y: f64,
    /// Optional symmetric error (used by [`PlotType::ErrorBar`]).
    pub error: Option<f64>,
}

impl DataPoint2D {
    /// Construct a simple point without an error bar.
    pub fn new(x: f64, y: f64) -> Self {
        Self { x, y, error: None }
    }

    /// Construct a point with a symmetric error bar of half-width `err`.
    pub fn with_error(x: f64, y: f64, err: f64) -> Self {
        Self {
            x,
            y,
            error: Some(err),
        }
    }
}

// ---------------------------------------------------------------------------
// PlotSeries
// ---------------------------------------------------------------------------

/// A single data series in a scientific plot.
#[derive(Debug, Clone)]
pub struct PlotSeries {
    /// Human-readable label for this series (used in the legend).
    pub label: String,
    /// The data points.
    pub data: Vec<DataPoint2D>,
    /// Visual representation type.
    pub plot_type: PlotType,
    /// Color scheme for this series.
    pub color: ColorScheme,
    /// Stroke style for line-based plot types.
    pub line_style: LineStyle,
    /// Marker shape.
    pub marker: MarkerShape,
    /// Stroke width in SVG user units.
    pub line_width: f64,
    /// Marker radius in SVG user units.
    pub marker_size: f64,
}

impl PlotSeries {
    /// Construct a new series with default styling.
    pub fn new(label: impl Into<String>, data: Vec<DataPoint2D>, plot_type: PlotType) -> Self {
        Self {
            label: label.into(),
            data,
            plot_type,
            color: ColorScheme::Custom(0.12, 0.47, 0.71),
            line_style: LineStyle::Solid,
            marker: MarkerShape::None,
            line_width: 1.5,
            marker_size: 4.0,
        }
    }

    /// Set the color scheme (builder pattern).
    pub fn with_color(mut self, color: ColorScheme) -> Self {
        self.color = color;
        self
    }

    /// Set the line style (builder pattern).
    pub fn with_line_style(mut self, style: LineStyle) -> Self {
        self.line_style = style;
        self
    }

    /// Set the marker shape (builder pattern).
    pub fn with_marker(mut self, marker: MarkerShape) -> Self {
        self.marker = marker;
        self
    }

    /// Return the x-range `(min, max)` of the data, or `None` if the data is empty.
    pub fn x_range(&self) -> Option<(f64, f64)> {
        if self.data.is_empty() {
            return None;
        }
        let min = self.data.iter().map(|p| p.x).fold(f64::INFINITY, f64::min);
        let max = self
            .data
            .iter()
            .map(|p| p.x)
            .fold(f64::NEG_INFINITY, f64::max);
        Some((min, max))
    }

    /// Return the y-range `(min, max)` of the data, or `None` if the data is empty.
    pub fn y_range(&self) -> Option<(f64, f64)> {
        if self.data.is_empty() {
            return None;
        }
        let min = self.data.iter().map(|p| p.y).fold(f64::INFINITY, f64::min);
        let max = self
            .data
            .iter()
            .map(|p| p.y)
            .fold(f64::NEG_INFINITY, f64::max);
        Some((min, max))
    }
}

// ---------------------------------------------------------------------------
// AxisScale
// ---------------------------------------------------------------------------

/// Scale type for an axis.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum AxisScale {
    /// Linear scale.
    Linear,
    /// Logarithmic (base-10) scale.
    Log,
}

// ---------------------------------------------------------------------------
// Axis
// ---------------------------------------------------------------------------

/// Configuration for a single plot axis.
#[derive(Debug, Clone)]
pub struct Axis {
    /// Axis label text.
    pub label: String,
    /// Explicit data range `(min, max)`.  `None` means auto-range.
    pub range: Option<(f64, f64)>,
    /// Scale type.
    pub scale: AxisScale,
    /// Tick spacing.  `None` means auto-spaced ticks.
    pub tick_spacing: Option<f64>,
    /// Number of grid lines to draw along this axis.
    pub grid_lines: usize,
}

impl Default for Axis {
    fn default() -> Self {
        Self {
            label: String::new(),
            range: None,
            scale: AxisScale::Linear,
            tick_spacing: None,
            grid_lines: 5,
        }
    }
}

impl Axis {
    /// Create an axis with the given label and automatic range.
    pub fn new(label: impl Into<String>) -> Self {
        Self {
            label: label.into(),
            ..Default::default()
        }
    }

    /// Set the explicit range (builder pattern).
    pub fn with_range(mut self, min: f64, max: f64) -> Self {
        self.range = Some((min, max));
        self
    }

    /// Set the scale type (builder pattern).
    pub fn with_scale(mut self, scale: AxisScale) -> Self {
        self.scale = scale;
        self
    }

    /// Set the tick spacing (builder pattern).
    pub fn with_tick_spacing(mut self, spacing: f64) -> Self {
        self.tick_spacing = Some(spacing);
        self
    }

    /// Compute a sensible range for this axis from a set of (min, max) hints.
    ///
    /// Returns the explicit range if set, otherwise pads the data range by 5 %.
    pub fn effective_range(&self, data_min: f64, data_max: f64) -> (f64, f64) {
        if let Some(r) = self.range {
            return r;
        }
        if (data_max - data_min).abs() < 1e-15 {
            return (data_min - 1.0, data_max + 1.0);
        }
        let pad = (data_max - data_min) * 0.05;
        (data_min - pad, data_max + pad)
    }

    /// Return tick positions for this axis given its effective range.
    ///
    /// Produces at most `grid_lines + 1` ticks.
    pub fn tick_positions(&self, lo: f64, hi: f64) -> Vec<f64> {
        let n = (self.grid_lines + 1).max(2);
        let step = (hi - lo) / (n - 1) as f64;
        (0..n).map(|i| lo + i as f64 * step).collect()
    }
}

// ---------------------------------------------------------------------------
// LegendPosition
// ---------------------------------------------------------------------------

/// Where to place the plot legend.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum LegendPosition {
    /// Top-left corner of the plot area.
    TopLeft,
    /// Top-right corner of the plot area.
    TopRight,
    /// Bottom-left corner of the plot area.
    BottomLeft,
    /// Bottom-right corner of the plot area.
    BottomRight,
    /// Do not render a legend.
    Hidden,
}

// ---------------------------------------------------------------------------
// PlotMargins
// ---------------------------------------------------------------------------

/// Margin (padding) around the plot area in SVG user units.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PlotMargins {
    /// Top margin.
    pub top: f64,
    /// Right margin.
    pub right: f64,
    /// Bottom margin.
    pub bottom: f64,
    /// Left margin.
    pub left: f64,
}

impl Default for PlotMargins {
    fn default() -> Self {
        Self {
            top: 40.0,
            right: 30.0,
            bottom: 60.0,
            left: 70.0,
        }
    }
}

// ---------------------------------------------------------------------------
// PlotLayout
// ---------------------------------------------------------------------------

/// Grid layout for a multi-panel figure.
#[derive(Debug, Clone)]
pub struct PlotLayout {
    /// Number of rows in the subplot grid.
    pub rows: usize,
    /// Number of columns in the subplot grid.
    pub cols: usize,
    /// Margins around each panel.
    pub margins: PlotMargins,
    /// Legend position for all panels.
    pub legend_position: LegendPosition,
    /// Horizontal spacing between panels in SVG user units.
    pub h_spacing: f64,
    /// Vertical spacing between panels in SVG user units.
    pub v_spacing: f64,
}

impl Default for PlotLayout {
    fn default() -> Self {
        Self {
            rows: 1,
            cols: 1,
            margins: PlotMargins::default(),
            legend_position: LegendPosition::TopRight,
            h_spacing: 20.0,
            v_spacing: 20.0,
        }
    }
}

impl PlotLayout {
    /// Construct a layout with `rows × cols` panels.
    pub fn new(rows: usize, cols: usize) -> Self {
        Self {
            rows: rows.max(1),
            cols: cols.max(1),
            ..Default::default()
        }
    }

    /// Total number of panels.
    pub fn n_panels(&self) -> usize {
        self.rows * self.cols
    }
}

// ---------------------------------------------------------------------------
// ScientificPlot
// ---------------------------------------------------------------------------

/// A complete scientific plot with title, axes, series, and layout.
#[derive(Debug, Clone)]
pub struct ScientificPlot {
    /// Plot title.
    pub title: String,
    /// Horizontal (x) axis configuration.
    pub x_axis: Axis,
    /// Vertical (y) axis configuration.
    pub y_axis: Axis,
    /// Ordered list of data series.
    pub series: Vec<PlotSeries>,
    /// Layout (subplot grid, margins, legend).
    pub layout: PlotLayout,
    /// Overall SVG canvas width in user units.
    pub width: f64,
    /// Overall SVG canvas height in user units.
    pub height: f64,
}

impl Default for ScientificPlot {
    fn default() -> Self {
        Self::new("", 640.0, 480.0)
    }
}

impl ScientificPlot {
    /// Construct an empty plot with the given title and canvas size.
    pub fn new(title: impl Into<String>, width: f64, height: f64) -> Self {
        Self {
            title: title.into(),
            x_axis: Axis::new("x"),
            y_axis: Axis::new("y"),
            series: Vec::new(),
            layout: PlotLayout::default(),
            width: width.max(1.0),
            height: height.max(1.0),
        }
    }

    /// Append a data series to the plot.
    pub fn add_series(&mut self, s: PlotSeries) {
        self.series.push(s);
    }

    /// Replace the x-axis configuration.
    pub fn set_x_axis(&mut self, axis: Axis) {
        self.x_axis = axis;
    }

    /// Replace the y-axis configuration.
    pub fn set_y_axis(&mut self, axis: Axis) {
        self.y_axis = axis;
    }

    /// Compute the bounding box of all data across all series.
    ///
    /// Returns `None` if there are no data points.
    pub fn data_bounds(&self) -> Option<(f64, f64, f64, f64)> {
        let (mut xmin, mut xmax) = (f64::INFINITY, f64::NEG_INFINITY);
        let (mut ymin, mut ymax) = (f64::INFINITY, f64::NEG_INFINITY);
        let mut any = false;

        for s in &self.series {
            for p in &s.data {
                if p.x < xmin {
                    xmin = p.x;
                }
                if p.x > xmax {
                    xmax = p.x;
                }
                if p.y < ymin {
                    ymin = p.y;
                }
                if p.y > ymax {
                    ymax = p.y;
                }
                any = true;
            }
        }
        if any {
            Some((xmin, xmax, ymin, ymax))
        } else {
            None
        }
    }

    /// Render the plot to an SVG string.
    ///
    /// Produces a self-contained SVG document with a title, axis labels,
    /// tick marks, grid lines, and all series.  Error bars are rendered
    /// for series of type [`PlotType::ErrorBar`].
    pub fn to_svg_string(&self) -> String {
        let m = &self.layout.margins;
        let plot_w = self.width - m.left - m.right;
        let plot_h = self.height - m.top - m.bottom;

        // Determine axis ranges
        let (x_data_min, x_data_max, y_data_min, y_data_max) =
            self.data_bounds().unwrap_or((0.0, 1.0, 0.0, 1.0));
        let (xlo, xhi) = self.x_axis.effective_range(x_data_min, x_data_max);
        let (ylo, yhi) = self.y_axis.effective_range(y_data_min, y_data_max);

        let x_range = (xhi - xlo).max(1e-15);
        let y_range = (yhi - ylo).max(1e-15);

        // Coordinate mappers: data → SVG pixel
        let map_x = |x: f64| m.left + (x - xlo) / x_range * plot_w;
        let map_y = |y: f64| m.top + plot_h - (y - ylo) / y_range * plot_h;

        let mut svg = String::with_capacity(8192);

        // --- SVG header ---
        svg.push_str(&format!(
            r#"<svg xmlns="http://www.w3.org/2000/svg" width="{}" height="{}">"#,
            self.width, self.height
        ));
        svg.push('\n');
        svg.push_str(r#"<rect width="100%" height="100%" fill="white"/>"#);
        svg.push('\n');

        // --- Title ---
        if !self.title.is_empty() {
            svg.push_str(&format!(
                r#"<text x="{}" y="20" text-anchor="middle" font-size="16" font-weight="bold">{}</text>"#,
                self.width / 2.0,
                escape_xml(&self.title)
            ));
            svg.push('\n');
        }

        // --- Plot area border ---
        svg.push_str(&format!(
            "<rect x=\"{}\" y=\"{}\" width=\"{}\" height=\"{}\" fill=\"none\" stroke=\"#333\" stroke-width=\"1\"/>",
            m.left, m.top, plot_w, plot_h
        ));
        svg.push('\n');

        // --- Grid lines ---
        let x_ticks = self.x_axis.tick_positions(xlo, xhi);
        let y_ticks = self.y_axis.tick_positions(ylo, yhi);

        for &xt in &x_ticks {
            let px = map_x(xt);
            svg.push_str(&format!(
                "<line x1=\"{px:.1}\" y1=\"{:.1}\" x2=\"{px:.1}\" y2=\"{:.1}\" stroke=\"#ddd\" stroke-width=\"0.5\"/>",
                m.top, m.top + plot_h
            ));
            svg.push('\n');
            // Tick label
            svg.push_str(&format!(
                r#"<text x="{px:.1}" y="{:.1}" text-anchor="middle" font-size="10">{:.3}</text>"#,
                m.top + plot_h + 15.0,
                xt
            ));
            svg.push('\n');
        }

        for &yt in &y_ticks {
            let py = map_y(yt);
            svg.push_str(&format!(
                "<line x1=\"{:.1}\" y1=\"{py:.1}\" x2=\"{:.1}\" y2=\"{py:.1}\" stroke=\"#ddd\" stroke-width=\"0.5\"/>",
                m.left, m.left + plot_w
            ));
            svg.push('\n');
            svg.push_str(&format!(
                r#"<text x="{:.1}" y="{py:.1}" text-anchor="end" dominant-baseline="middle" font-size="10">{:.3}</text>"#,
                m.left - 5.0,
                yt
            ));
            svg.push('\n');
        }

        // --- Axis labels ---
        svg.push_str(&format!(
            r#"<text x="{:.1}" y="{:.1}" text-anchor="middle" font-size="12">{}</text>"#,
            m.left + plot_w / 2.0,
            m.top + plot_h + 45.0,
            escape_xml(&self.x_axis.label)
        ));
        svg.push('\n');

        svg.push_str(&format!(
            r#"<text x="-{:.1}" y="{:.1}" text-anchor="middle" font-size="12" transform="rotate(-90)">{}</text>"#,
            m.top + plot_h / 2.0,
            15.0,
            escape_xml(&self.y_axis.label)
        ));
        svg.push('\n');

        // --- Series ---
        for (idx, series) in self.series.iter().enumerate() {
            let color = if series.color == ColorScheme::Grayscale
                || series.color == ColorScheme::Viridis
                || series.color == ColorScheme::Jet
            {
                series.color.to_svg_color()
            } else {
                // Use palette if the user hasn't set a custom color explicitly
                // (default color is ~#1f77b4 which is palette index 0)
                ColorScheme::palette_color(idx)
            };

            match series.plot_type {
                PlotType::Line | PlotType::ErrorBar => {
                    // Draw the line
                    if series.data.len() >= 2 {
                        let mut path = String::from("M");
                        for (i, p) in series.data.iter().enumerate() {
                            let px = map_x(p.x);
                            let py = map_y(p.y);
                            if i == 0 {
                                path.push_str(&format!("{px:.1},{py:.1}"));
                            } else {
                                path.push_str(&format!(" L{px:.1},{py:.1}"));
                            }
                        }
                        let dash = series
                            .line_style
                            .svg_dash_array()
                            .map(|d| format!(r#" stroke-dasharray="{d}""#))
                            .unwrap_or_default();
                        svg.push_str(&format!(
                            r#"<path d="{path}" fill="none" stroke="{color}" stroke-width="{:.1}"{dash}/>"#,
                            series.line_width
                        ));
                        svg.push('\n');
                    }
                    // Draw error bars
                    if series.plot_type == PlotType::ErrorBar {
                        for p in &series.data {
                            if let Some(err) = p.error {
                                let px = map_x(p.x);
                                let py_lo = map_y(p.y - err);
                                let py_hi = map_y(p.y + err);
                                svg.push_str(&format!(
                                    r#"<line x1="{px:.1}" y1="{py_lo:.1}" x2="{px:.1}" y2="{py_hi:.1}" stroke="{color}" stroke-width="1"/>"#
                                ));
                                svg.push('\n');
                                // Cap lines
                                for py_cap in [py_lo, py_hi] {
                                    svg.push_str(&format!(
                                        r#"<line x1="{:.1}" y1="{py_cap:.1}" x2="{:.1}" y2="{py_cap:.1}" stroke="{color}" stroke-width="1"/>"#,
                                        px - 4.0, px + 4.0
                                    ));
                                    svg.push('\n');
                                }
                            }
                        }
                    }
                }
                PlotType::Scatter => {
                    for p in &series.data {
                        let px = map_x(p.x);
                        let py = map_y(p.y);
                        svg.push_str(&render_marker(
                            series.marker,
                            px,
                            py,
                            series.marker_size,
                            &color,
                        ));
                        svg.push('\n');
                    }
                }
                PlotType::Bar => {
                    let bar_w = if series.data.len() > 1 {
                        let step = (series.data[1].x - series.data[0].x).abs();
                        map_x(xlo + step) - map_x(xlo)
                    } else {
                        20.0
                    } * 0.8;
                    let baseline = map_y(0.0_f64.max(ylo));
                    for p in &series.data {
                        let px = map_x(p.x) - bar_w / 2.0;
                        let py = map_y(p.y);
                        let bar_h = (baseline - py).abs();
                        let bar_top = py.min(baseline);
                        svg.push_str(&format!(
                            r#"<rect x="{px:.1}" y="{bar_top:.1}" width="{bar_w:.1}" height="{bar_h:.1}" fill="{color}" opacity="0.8"/>"#
                        ));
                        svg.push('\n');
                    }
                }
                PlotType::Histogram => {
                    // Re-use bar rendering
                    for p in &series.data {
                        let px = map_x(p.x) - 5.0;
                        let py = map_y(p.y);
                        let bar_h = map_y(ylo) - py;
                        svg.push_str(&format!(
                            r#"<rect x="{px:.1}" y="{py:.1}" width="10" height="{bar_h:.1}" fill="{color}" opacity="0.7"/>"#
                        ));
                        svg.push('\n');
                    }
                }
                PlotType::Contour | PlotType::Surface | PlotType::VectorField => {
                    // For complex plot types, emit a placeholder text annotation.
                    svg.push_str(&format!(
                        "<text x=\"{:.1}\" y=\"{:.1}\" font-size=\"10\" fill=\"#888\">[{:?} data, n={}]</text>",
                        m.left + 5.0,
                        m.top + 14.0 * (idx as f64 + 1.0),
                        series.plot_type,
                        series.data.len()
                    ));
                    svg.push('\n');
                }
            }

            // Marker overlay for Line plots
            if series.plot_type == PlotType::Line && series.marker != MarkerShape::None {
                for p in &series.data {
                    let px = map_x(p.x);
                    let py = map_y(p.y);
                    svg.push_str(&render_marker(
                        series.marker,
                        px,
                        py,
                        series.marker_size,
                        &color,
                    ));
                    svg.push('\n');
                }
            }
        }

        // --- Legend ---
        if self.layout.legend_position != LegendPosition::Hidden && !self.series.is_empty() {
            let leg_x = match self.layout.legend_position {
                LegendPosition::TopLeft | LegendPosition::BottomLeft => m.left + 8.0,
                _ => m.left + plot_w - 140.0,
            };
            let leg_y = match self.layout.legend_position {
                LegendPosition::TopLeft | LegendPosition::TopRight => m.top + 8.0,
                _ => m.top + plot_h - (self.series.len() as f64 * 18.0 + 10.0),
            };
            let leg_h = self.series.len() as f64 * 18.0 + 10.0;
            svg.push_str(&format!(
                "<rect x=\"{leg_x:.1}\" y=\"{leg_y:.1}\" width=\"135\" height=\"{leg_h:.1}\" fill=\"white\" stroke=\"#aaa\" stroke-width=\"0.5\" rx=\"3\"/>"
            ));
            svg.push('\n');
            for (i, s) in self.series.iter().enumerate() {
                let row_y = leg_y + 8.0 + i as f64 * 18.0 + 9.0;
                let c = ColorScheme::palette_color(i);
                svg.push_str(&format!(
                    r#"<line x1="{:.1}" y1="{row_y:.1}" x2="{:.1}" y2="{row_y:.1}" stroke="{c}" stroke-width="2"/>"#,
                    leg_x + 5.0,
                    leg_x + 25.0
                ));
                svg.push('\n');
                svg.push_str(&format!(
                    r#"<text x="{:.1}" y="{:.1}" font-size="11" dominant-baseline="middle">{}</text>"#,
                    leg_x + 30.0,
                    row_y,
                    escape_xml(&s.label)
                ));
                svg.push('\n');
            }
        }

        svg.push_str("</svg>\n");
        svg
    }
}

// ---------------------------------------------------------------------------
// Internal helpers
// ---------------------------------------------------------------------------

/// Render an SVG shape for a data-point marker.
fn render_marker(shape: MarkerShape, cx: f64, cy: f64, r: f64, color: &str) -> String {
    match shape {
        MarkerShape::None => String::new(),
        MarkerShape::Circle => {
            format!(r#"<circle cx="{cx:.1}" cy="{cy:.1}" r="{r:.1}" fill="{color}"/>"#)
        }
        MarkerShape::Square => {
            let s = r * 1.5;
            format!(
                r#"<rect x="{:.1}" y="{:.1}" width="{:.1}" height="{:.1}" fill="{color}"/>"#,
                cx - s / 2.0,
                cy - s / 2.0,
                s,
                s
            )
        }
        MarkerShape::Triangle => {
            let pts = format!(
                "{:.1},{:.1} {:.1},{:.1} {:.1},{:.1}",
                cx,
                cy - r,
                cx - r,
                cy + r * 0.577,
                cx + r,
                cy + r * 0.577
            );
            format!(r#"<polygon points="{pts}" fill="{color}"/>"#)
        }
        MarkerShape::Diamond => {
            let pts = format!(
                "{:.1},{:.1} {:.1},{:.1} {:.1},{:.1} {:.1},{:.1}",
                cx,
                cy - r,
                cx + r,
                cy,
                cx,
                cy + r,
                cx - r,
                cy
            );
            format!(r#"<polygon points="{pts}" fill="{color}"/>"#)
        }
        MarkerShape::Cross => {
            format!(
                r#"<line x1="{:.1}" y1="{cy:.1}" x2="{:.1}" y2="{cy:.1}" stroke="{color}" stroke-width="2"/><line x1="{cx:.1}" y1="{:.1}" x2="{cx:.1}" y2="{:.1}" stroke="{color}" stroke-width="2"/>"#,
                cx - r,
                cx + r,
                cy - r,
                cy + r
            )
        }
    }
}

/// Escape characters that are special in XML/SVG text nodes.
fn escape_xml(s: &str) -> String {
    s.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
        .replace('"', "&quot;")
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    fn line_series(n: usize) -> PlotSeries {
        let data = (0..n)
            .map(|i| DataPoint2D::new(i as f64, (i as f64).sin()))
            .collect();
        PlotSeries::new("sin", data, PlotType::Line)
    }

    fn scatter_series() -> PlotSeries {
        let data = vec![
            DataPoint2D::new(0.0, 0.0),
            DataPoint2D::new(1.0, 2.0),
            DataPoint2D::new(2.0, 1.0),
        ];
        PlotSeries::new("scatter", data, PlotType::Scatter)
    }

    // ---- ColorScheme tests ----

    #[test]
    fn test_color_scheme_custom_svg() {
        let c = ColorScheme::Custom(1.0, 0.0, 0.0);
        assert_eq!(c.to_svg_color(), "rgb(255,0,0)");
    }

    #[test]
    fn test_color_scheme_grayscale_svg() {
        let svg = ColorScheme::Grayscale.to_svg_color();
        assert!(svg.starts_with("rgb("));
    }

    #[test]
    fn test_color_scheme_custom_clamp() {
        let c = ColorScheme::Custom(2.0, -1.0, 0.5);
        let svg = c.to_svg_color();
        // r should be clamped to 255, g to 0
        assert!(svg.contains("255"));
        assert!(svg.contains(",0,"));
    }

    #[test]
    fn test_palette_color_cycles() {
        let c0 = ColorScheme::palette_color(0);
        let c10 = ColorScheme::palette_color(10); // wraps back to index 0
        assert_eq!(c0, c10);
    }

    // ---- LineStyle tests ----

    #[test]
    fn test_line_style_solid_has_no_dash() {
        assert!(LineStyle::Solid.svg_dash_array().is_none());
    }

    #[test]
    fn test_line_style_dashed_has_dash() {
        assert!(LineStyle::Dashed.svg_dash_array().is_some());
    }

    #[test]
    fn test_line_style_dotted_dash() {
        let d = LineStyle::Dotted.svg_dash_array().unwrap();
        assert!(d.contains("2"));
    }

    // ---- DataPoint2D tests ----

    #[test]
    fn test_data_point_no_error() {
        let p = DataPoint2D::new(1.0, 2.0);
        assert_eq!(p.x, 1.0);
        assert_eq!(p.y, 2.0);
        assert!(p.error.is_none());
    }

    #[test]
    fn test_data_point_with_error() {
        let p = DataPoint2D::with_error(1.0, 2.0, 0.5);
        assert_eq!(p.error, Some(0.5));
    }

    // ---- PlotSeries tests ----

    #[test]
    fn test_series_x_range() {
        let s = line_series(5);
        let (lo, hi) = s.x_range().unwrap();
        assert_eq!(lo, 0.0);
        assert_eq!(hi, 4.0);
    }

    #[test]
    fn test_series_y_range_nonempty() {
        let s = line_series(10);
        let (lo, hi) = s.y_range().unwrap();
        assert!(lo <= hi);
    }

    #[test]
    fn test_series_empty_range_none() {
        let s = PlotSeries::new("empty", vec![], PlotType::Line);
        assert!(s.x_range().is_none());
        assert!(s.y_range().is_none());
    }

    #[test]
    fn test_series_builder_color() {
        let s = line_series(3).with_color(ColorScheme::Viridis);
        assert_eq!(s.color, ColorScheme::Viridis);
    }

    #[test]
    fn test_series_builder_line_style() {
        let s = line_series(3).with_line_style(LineStyle::Dashed);
        assert_eq!(s.line_style, LineStyle::Dashed);
    }

    #[test]
    fn test_series_builder_marker() {
        let s = scatter_series().with_marker(MarkerShape::Circle);
        assert_eq!(s.marker, MarkerShape::Circle);
    }

    // ---- Axis tests ----

    #[test]
    fn test_axis_default_auto_range() {
        let ax = Axis::new("x");
        assert!(ax.range.is_none());
    }

    #[test]
    fn test_axis_explicit_range() {
        let ax = Axis::new("x").with_range(-5.0, 5.0);
        assert_eq!(ax.range, Some((-5.0, 5.0)));
    }

    #[test]
    fn test_axis_effective_range_respects_explicit() {
        let ax = Axis::new("x").with_range(0.0, 10.0);
        let (lo, hi) = ax.effective_range(-100.0, 100.0);
        assert_eq!(lo, 0.0);
        assert_eq!(hi, 10.0);
    }

    #[test]
    fn test_axis_effective_range_auto_pads() {
        let ax = Axis::default();
        let (lo, hi) = ax.effective_range(0.0, 10.0);
        assert!(lo < 0.0);
        assert!(hi > 10.0);
    }

    #[test]
    fn test_axis_effective_range_degenerate() {
        let ax = Axis::default();
        let (lo, hi) = ax.effective_range(5.0, 5.0);
        assert!(lo < 5.0 && hi > 5.0);
    }

    #[test]
    fn test_axis_tick_positions_count() {
        let ax = Axis {
            grid_lines: 4,
            ..Default::default()
        };
        let ticks = ax.tick_positions(0.0, 1.0);
        assert_eq!(ticks.len(), 5); // grid_lines + 1
    }

    #[test]
    fn test_axis_tick_positions_include_endpoints() {
        let ax = Axis::default();
        let ticks = ax.tick_positions(0.0, 1.0);
        assert!((ticks[0] - 0.0).abs() < 1e-12);
        assert!((ticks.last().unwrap() - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_axis_scale_log() {
        let ax = Axis::new("y").with_scale(AxisScale::Log);
        assert_eq!(ax.scale, AxisScale::Log);
    }

    // ---- PlotLayout tests ----

    #[test]
    fn test_layout_n_panels() {
        let layout = PlotLayout::new(2, 3);
        assert_eq!(layout.n_panels(), 6);
    }

    #[test]
    fn test_layout_min_one_panel() {
        let layout = PlotLayout::new(0, 0);
        assert_eq!(layout.n_panels(), 1);
    }

    // ---- ScientificPlot construction ----

    #[test]
    fn test_plot_add_series() {
        let mut plot = ScientificPlot::new("test", 640.0, 480.0);
        plot.add_series(line_series(5));
        assert_eq!(plot.series.len(), 1);
    }

    #[test]
    fn test_plot_set_axes() {
        let mut plot = ScientificPlot::new("test", 640.0, 480.0);
        plot.set_x_axis(Axis::new("time (s)").with_range(0.0, 10.0));
        plot.set_y_axis(Axis::new("amplitude").with_range(-1.0, 1.0));
        assert_eq!(plot.x_axis.label, "time (s)");
        assert_eq!(plot.y_axis.range, Some((-1.0, 1.0)));
    }

    #[test]
    fn test_plot_data_bounds_empty() {
        let plot = ScientificPlot::default();
        assert!(plot.data_bounds().is_none());
    }

    #[test]
    fn test_plot_data_bounds_with_series() {
        let mut plot = ScientificPlot::default();
        plot.add_series(line_series(5));
        let (xmin, xmax, _, _) = plot.data_bounds().unwrap();
        assert_eq!(xmin, 0.0);
        assert_eq!(xmax, 4.0);
    }

    // ---- SVG output tests ----

    #[test]
    fn test_svg_contains_svg_tag() {
        let plot = ScientificPlot::new("Title", 640.0, 480.0);
        let svg = plot.to_svg_string();
        assert!(svg.contains("<svg"), "SVG output must open with <svg");
        assert!(svg.contains("</svg>"), "SVG output must close with </svg>");
    }

    #[test]
    fn test_svg_contains_title() {
        let mut plot = ScientificPlot::new("My Plot", 640.0, 480.0);
        plot.add_series(line_series(5));
        let svg = plot.to_svg_string();
        assert!(svg.contains("My Plot"), "SVG should contain the title text");
    }

    #[test]
    fn test_svg_contains_axis_label() {
        let mut plot = ScientificPlot::new("P", 640.0, 480.0);
        plot.set_x_axis(Axis::new("Distance [m]"));
        plot.add_series(line_series(3));
        let svg = plot.to_svg_string();
        assert!(
            svg.contains("Distance [m]"),
            "SVG should contain x-axis label"
        );
    }

    #[test]
    fn test_svg_line_series_contains_path() {
        let mut plot = ScientificPlot::new("P", 640.0, 480.0);
        plot.add_series(line_series(5));
        let svg = plot.to_svg_string();
        assert!(
            svg.contains("<path"),
            "Line series should produce <path> element"
        );
    }

    #[test]
    fn test_svg_scatter_contains_circle() {
        let mut plot = ScientificPlot::new("P", 640.0, 480.0);
        let s = scatter_series().with_marker(MarkerShape::Circle);
        plot.add_series(s);
        let svg = plot.to_svg_string();
        assert!(
            svg.contains("<circle"),
            "Scatter series should produce circles"
        );
    }

    #[test]
    fn test_svg_bar_series_contains_rect() {
        let mut plot = ScientificPlot::new("P", 640.0, 480.0);
        let data = vec![
            DataPoint2D::new(1.0, 3.0),
            DataPoint2D::new(2.0, 5.0),
            DataPoint2D::new(3.0, 2.0),
        ];
        plot.add_series(PlotSeries::new("bars", data, PlotType::Bar));
        let svg = plot.to_svg_string();
        assert!(
            svg.contains("<rect"),
            "Bar chart should produce <rect> elements"
        );
    }

    #[test]
    fn test_svg_error_bar_present() {
        let mut plot = ScientificPlot::new("P", 640.0, 480.0);
        let data = vec![
            DataPoint2D::with_error(0.0, 1.0, 0.2),
            DataPoint2D::with_error(1.0, 2.0, 0.3),
        ];
        plot.add_series(PlotSeries::new("eb", data, PlotType::ErrorBar));
        let svg = plot.to_svg_string();
        assert!(
            svg.contains("<line"),
            "Error bar should contain line elements"
        );
    }

    #[test]
    fn test_svg_no_title_when_empty() {
        let plot = ScientificPlot::new("", 640.0, 480.0);
        let svg = plot.to_svg_string();
        // No title text element should be present
        assert!(
            !svg.contains("font-weight=\"bold\""),
            "empty title should not render title element"
        );
    }

    #[test]
    fn test_svg_legend_present() {
        let mut plot = ScientificPlot::new("P", 640.0, 480.0);
        plot.add_series(line_series(3));
        let svg = plot.to_svg_string();
        // Legend contains the series label
        assert!(
            svg.contains("sin"),
            "SVG legend should contain series label"
        );
    }

    #[test]
    fn test_svg_hidden_legend() {
        let mut plot = ScientificPlot::new("P", 640.0, 480.0);
        plot.layout.legend_position = LegendPosition::Hidden;
        plot.add_series(line_series(3));
        let svg = plot.to_svg_string();
        // The legend rectangle should be absent
        assert!(
            !svg.contains("stroke=\"#aaa\""),
            "hidden legend should not render border"
        );
    }

    #[test]
    fn test_xml_escape() {
        let s = escape_xml("a<b>&\"c");
        assert_eq!(s, "a&lt;b&gt;&amp;&quot;c");
    }

    #[test]
    fn test_svg_size_attribute() {
        let plot = ScientificPlot::new("T", 800.0, 600.0);
        let svg = plot.to_svg_string();
        assert!(svg.contains("width=\"800\""));
        assert!(svg.contains("height=\"600\""));
    }

    #[test]
    fn test_svg_multiple_series() {
        let mut plot = ScientificPlot::new("Multi", 640.0, 480.0);
        plot.add_series(line_series(5));
        plot.add_series(scatter_series().with_marker(MarkerShape::Circle));
        let svg = plot.to_svg_string();
        assert!(svg.contains("<path"), "first series renders line");
        assert!(svg.contains("<circle"), "second series renders circles");
    }

    #[test]
    fn test_render_marker_none_empty() {
        let s = render_marker(MarkerShape::None, 10.0, 10.0, 4.0, "red");
        assert!(s.is_empty());
    }

    #[test]
    fn test_render_marker_diamond_polygon() {
        let s = render_marker(MarkerShape::Diamond, 10.0, 10.0, 4.0, "blue");
        assert!(s.contains("<polygon"));
    }

    #[test]
    fn test_render_marker_cross_lines() {
        let s = render_marker(MarkerShape::Cross, 10.0, 10.0, 4.0, "green");
        assert!(s.contains("<line"));
    }
}
