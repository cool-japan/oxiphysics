//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::*;
use crate::primitives::Color;

/// A detailed legend for a colormap, including raster and tick labels.
pub struct ColormapLegend {
    /// Colormap to display.
    pub colormap: Colormap,
    /// Minimum value mapped to t=0.
    pub vmin: f64,
    /// Maximum value mapped to t=1.
    pub vmax: f64,
    /// Label text for the legend.
    pub label: String,
    /// Number of tick marks.
    pub n_ticks: usize,
}
impl ColormapLegend {
    /// Create a new legend.
    pub fn new(colormap: Colormap, vmin: f64, vmax: f64, label: &str) -> Self {
        Self {
            colormap,
            vmin,
            vmax,
            label: label.to_owned(),
            n_ticks: 5,
        }
    }
    /// Set the number of tick marks.
    pub fn with_n_ticks(mut self, n: usize) -> Self {
        self.n_ticks = n.max(2);
        self
    }
    /// Generate a flat 2D raster `[width * height]` of RGBA `[u8;4]` pixels.
    ///
    /// The raster is a horizontal gradient strip `width` pixels wide and
    /// `height` rows tall.
    pub fn render_raster(&self, width: usize, height: usize) -> Vec<[u8; 4]> {
        if width == 0 || height == 0 {
            return Vec::new();
        }
        let mut pixels = Vec::with_capacity(width * height);
        for _row in 0..height {
            for col in 0..width {
                let t = col as f64 / (width - 1).max(1) as f64;
                let value = self.vmin + t * (self.vmax - self.vmin);
                let color = map_scalar(value, self.vmin, self.vmax, self.colormap);
                pixels.push([
                    (color.r.clamp(0.0, 1.0) * 255.0).round() as u8,
                    (color.g.clamp(0.0, 1.0) * 255.0).round() as u8,
                    (color.b.clamp(0.0, 1.0) * 255.0).round() as u8,
                    255,
                ]);
            }
        }
        pixels
    }
    /// Generate formatted tick labels.
    pub fn tick_labels(&self) -> Vec<String> {
        let n = self.n_ticks.max(2);
        (0..n)
            .map(|i| {
                let t = i as f64 / (n - 1) as f64;
                let v = self.vmin + t * (self.vmax - self.vmin);
                let range = (self.vmax - self.vmin).abs();
                if range >= 100.0 {
                    format!("{:.0}", v)
                } else if range >= 1.0 {
                    format!("{:.2}", v)
                } else {
                    format!("{:.4}", v)
                }
            })
            .collect()
    }
    /// Compute the pixel x-coordinate for a given value.
    ///
    /// Returns `None` if the value is outside `[vmin, vmax]`.
    pub fn value_to_pixel_x(&self, value: f64, width: usize) -> Option<usize> {
        let range = self.vmax - self.vmin;
        if range.abs() < 1e-30 {
            return Some(width / 2);
        }
        let t = (value - self.vmin) / range;
        if !(0.0..=1.0).contains(&t) {
            return None;
        }
        Some((t * (width - 1) as f64).round() as usize)
    }
}
/// Interpolation mode for colormap sampling.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum InterpolationMode {
    /// Linear interpolation between stops.
    Linear,
    /// Nearest-neighbor (step to the nearest stop).
    Nearest,
    /// Step function (uses the left stop color until the next stop).
    Step,
}
/// A color stop used by `CustomColormapBuilder` and `SplineColormap`.
#[derive(Debug, Clone)]
pub struct ColorStop {
    /// Position in \[0, 1\].
    pub t: f64,
    /// RGBA color at this position.
    pub color: Color,
}
impl ColorStop {
    fn s_t(&self) -> f64 {
        self.t
    }
}
/// Builder for constructing colormaps from user-defined color stops.
///
/// Stops are stored sorted by `t`.  When fewer than two stops are provided,
/// [`build`](CustomColormapBuilder::build) will still return a valid (constant)
/// colormap that clamps to the single stop color.
#[derive(Debug, Clone, Default)]
pub struct CustomColormapBuilder {
    pub(super) stops: Vec<ColorStop>,
    /// Interpolation mode applied when sampling the built colormap.
    pub interpolation: LinearInterp,
}
impl CustomColormapBuilder {
    /// Create a new builder with no stops.
    pub fn new() -> Self {
        Self::default()
    }
    /// Add a color stop at position `t` (clamped to \[0, 1\]).
    pub fn add_stop(mut self, t: f64, color: Color) -> Self {
        let t = t.clamp(0.0, 1.0);
        self.stops.push(ColorStop { t, color });
        self.stops
            .sort_by(|a, b| a.t.partial_cmp(&b.t).unwrap_or(std::cmp::Ordering::Equal));
        self
    }
    /// Set the interpolation mode.
    pub fn with_interpolation(mut self, mode: LinearInterp) -> Self {
        self.interpolation = mode;
        self
    }
    /// Sample the colormap at `t` in \[0, 1\].
    ///
    /// Returns a linearly (or smooth-step) interpolated `Color`.
    pub fn sample(&self, t: f64) -> Color {
        let t = t.clamp(0.0, 1.0);
        if self.stops.is_empty() {
            return Color {
                r: 0.0,
                g: 0.0,
                b: 0.0,
                a: 1.0,
            };
        }
        if self.stops.len() == 1 {
            return self.stops[0].color;
        }
        let first = &self.stops[0];
        let last = &self.stops[self.stops.len() - 1];
        if t <= first.t {
            return first.color;
        }
        if t >= last.t {
            return last.color;
        }
        let idx = self.stops.partition_point(|s| s.s_t() <= t) - 1;
        let lo = &self.stops[idx];
        let hi = &self.stops[(idx + 1).min(self.stops.len() - 1)];
        let span = hi.t - lo.t;
        let raw_alpha = if span < 1e-12 { 0.0 } else { (t - lo.t) / span };
        let alpha = match self.interpolation {
            LinearInterp::Linear => raw_alpha,
            LinearInterp::SmoothStep => raw_alpha * raw_alpha * (3.0 - 2.0 * raw_alpha),
        };
        lerp_color(lo.color, hi.color, alpha as f32)
    }
    /// Build a list of `n` evenly-spaced colors (for use as a discrete palette).
    pub fn build(&self, n: usize) -> Vec<Color> {
        let n = n.max(1);
        (0..n)
            .map(|i| {
                let t = if n == 1 {
                    0.5
                } else {
                    i as f64 / (n - 1) as f64
                };
                self.sample(t)
            })
            .collect()
    }
}
/// A colorbar for legend visualization.
pub struct Colorbar {
    /// Colormap to use.
    pub colormap: Colormap,
    /// Minimum value.
    pub vmin: f64,
    /// Maximum value.
    pub vmax: f64,
    /// Label for the colorbar.
    pub label: String,
    /// Number of tick marks.
    pub n_ticks: usize,
}
impl Colorbar {
    /// Create a new colorbar.
    pub fn new(colormap: Colormap, vmin: f64, vmax: f64, label: &str) -> Self {
        Self {
            colormap,
            vmin,
            vmax,
            label: label.to_owned(),
            n_ticks: 5,
        }
    }
    /// Set the number of ticks.
    pub fn with_ticks(mut self, n: usize) -> Self {
        self.n_ticks = n;
        self
    }
    /// Generate the tick values.
    pub fn tick_values(&self) -> Vec<f64> {
        if self.n_ticks < 2 {
            return vec![self.vmin];
        }
        let step = (self.vmax - self.vmin) / (self.n_ticks - 1) as f64;
        (0..self.n_ticks)
            .map(|i| self.vmin + step * i as f64)
            .collect()
    }
    /// Generate a raster of colors for the colorbar.
    ///
    /// Returns `n_pixels` colors sampled uniformly from vmin to vmax.
    pub fn raster(&self, n_pixels: usize) -> Vec<Color> {
        if n_pixels == 0 {
            return Vec::new();
        }
        (0..n_pixels)
            .map(|i| {
                let t = if n_pixels == 1 {
                    0.5
                } else {
                    i as f64 / (n_pixels - 1) as f64
                };
                let value = self.vmin + t * (self.vmax - self.vmin);
                map_scalar(value, self.vmin, self.vmax, self.colormap)
            })
            .collect()
    }
}
/// Scientific colormap producing RGBA `[u8; 4]` output.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ColorMap {
    /// Viridis (dark purple -> yellow).
    Viridis,
    /// Plasma (dark blue -> yellow).
    Plasma,
    /// Inferno (black -> yellow-white).
    Inferno,
    /// Magma (black -> white-yellow).
    Magma,
    /// Turbo (Google rainbow).
    Turbo,
    /// Red-Blue diverging.
    RdBu,
    /// Seismic diverging (blue-white-red).
    Seismic,
    /// Brown-Blue-Green diverging.
    BrBG,
}
impl ColorMap {
    /// Sample the colormap at parameter `t in [0, 1]`.
    ///
    /// Returns an RGBA `[u8; 4]` array (alpha is always 255).
    pub fn sample(self, t: f64) -> [u8; 4] {
        let t_clamped = t.clamp(0.0, 1.0) as f32;
        let col = match self {
            ColorMap::Viridis => viridis(t_clamped),
            ColorMap::Plasma => plasma(t_clamped),
            ColorMap::Inferno => inferno(t_clamped),
            ColorMap::Magma => magma(t_clamped),
            ColorMap::Turbo => turbo(t_clamped),
            ColorMap::RdBu => rdbu(t_clamped),
            ColorMap::Seismic => seismic(t_clamped),
            ColorMap::BrBG => brbg(t_clamped),
        };
        [
            (col.r.clamp(0.0, 1.0) * 255.0).round() as u8,
            (col.g.clamp(0.0, 1.0) * 255.0).round() as u8,
            (col.b.clamp(0.0, 1.0) * 255.0).round() as u8,
            255,
        ]
    }
}
impl ColorMap {
    /// Sample the Cividis perceptually-uniform colormap.
    pub fn sample_cividis(t: f64) -> [u8; 4] {
        let col = cividis_fn(t.clamp(0.0, 1.0) as f32);
        [
            (col.r.clamp(0.0, 1.0) * 255.0).round() as u8,
            (col.g.clamp(0.0, 1.0) * 255.0).round() as u8,
            (col.b.clamp(0.0, 1.0) * 255.0).round() as u8,
            255,
        ]
    }
    /// Sample the IsoRainbow isoluminant colormap.
    pub fn sample_iso_rainbow(t: f64) -> [u8; 4] {
        let col = iso_rainbow_fn(t.clamp(0.0, 1.0) as f32);
        [
            (col.r.clamp(0.0, 1.0) * 255.0).round() as u8,
            (col.g.clamp(0.0, 1.0) * 255.0).round() as u8,
            (col.b.clamp(0.0, 1.0) * 255.0).round() as u8,
            255,
        ]
    }
}
/// A colormap that wraps smoothly at `t = 0` / `t = 1`.
///
/// Useful for periodic data such as angles, phases, or times of day.
#[derive(Debug, Clone)]
pub struct CyclicColormap {
    /// The cyclic color style variant to use.
    pub style: CyclicStyle,
    /// Number of full cycles across \[0, 1\].  Default is 1.
    pub cycles: f64,
}
impl CyclicColormap {
    /// Create a new cyclic colormap with a given style and one cycle.
    pub fn new(style: CyclicStyle) -> Self {
        Self { style, cycles: 1.0 }
    }
    /// Set the number of full cycles across \[0, 1\].
    pub fn with_cycles(mut self, cycles: f64) -> Self {
        self.cycles = cycles.max(0.01);
        self
    }
    /// Sample the colormap at `t` in \[0, 1\].
    pub fn sample(&self, t: f64) -> Color {
        let phase = (t * self.cycles).fract();
        let phase = if phase < 0.0 { phase + 1.0 } else { phase };
        match self.style {
            CyclicStyle::HsvWheel => {
                let (r, g, b) = hsv_to_rgb(phase, 1.0, 1.0);
                Color {
                    r: r as f32,
                    g: g as f32,
                    b: b as f32,
                    a: 1.0,
                }
            }
            CyclicStyle::Twilight => sample_twilight(phase as f32),
            CyclicStyle::Phase => sample_phase(phase as f32),
        }
    }
    /// Build a discrete palette of `n` colors.
    pub fn build(&self, n: usize) -> Vec<Color> {
        let n = n.max(1);
        (0..n).map(|i| self.sample(i as f64 / n as f64)).collect()
    }
}
/// A perceptual-ish colormap that sweeps through HSV hue space.
///
/// Unlike [`CyclicStyle::HsvWheel`] which always cycles fully, `HsvColormap`
/// can sweep a *partial* hue range and supports configurable saturation and
/// value (brightness).
#[derive(Debug, Clone)]
pub struct HsvColormap {
    /// Starting hue in \[0, 1\] (0 = red, 1/3 = green, 2/3 = blue).
    pub hue_start: f64,
    /// Ending hue in \[0, 1\].
    pub hue_end: f64,
    /// Saturation (0 = grey, 1 = fully saturated).
    pub saturation: f64,
    /// Value / brightness (0 = black, 1 = full brightness).
    pub value: f64,
}
impl HsvColormap {
    /// Full rainbow from red → violet.
    pub fn rainbow() -> Self {
        Self {
            hue_start: 0.0,
            hue_end: 5.0 / 6.0,
            saturation: 1.0,
            value: 1.0,
        }
    }
    /// Create a custom HSV colormap.
    pub fn new(hue_start: f64, hue_end: f64, saturation: f64, value: f64) -> Self {
        Self {
            hue_start: hue_start.clamp(0.0, 1.0),
            hue_end: hue_end.clamp(0.0, 1.0),
            saturation: saturation.clamp(0.0, 1.0),
            value: value.clamp(0.0, 1.0),
        }
    }
    /// Sample at `t` in \[0, 1\].
    pub fn sample(&self, t: f64) -> Color {
        let t = t.clamp(0.0, 1.0);
        let hue = self.hue_start + (self.hue_end - self.hue_start) * t;
        let (r, g, b) = hsv_to_rgb(hue, self.saturation, self.value);
        Color {
            r: r as f32,
            g: g as f32,
            b: b as f32,
            a: 1.0,
        }
    }
    /// Build a discrete palette of `n` evenly-spaced colors.
    pub fn build(&self, n: usize) -> Vec<Color> {
        let n = n.max(1);
        (0..n)
            .map(|i| {
                let t = if n == 1 {
                    0.5
                } else {
                    i as f64 / (n - 1) as f64
                };
                self.sample(t)
            })
            .collect()
    }
}
/// Linear interpolation strategy marker used by `CustomColormapBuilder`.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub enum LinearInterp {
    /// Simple linear blend between adjacent stops.
    #[default]
    Linear,
    /// Smooth-step (cubic Hermite) blending.
    SmoothStep,
}
/// A categorical (qualitative) colormap for discrete categories.
///
/// Uses a fixed palette of distinct colors.
pub struct CategoricalColormap {
    /// The palette of distinct colors.
    pub colors: Vec<Color>,
}
impl CategoricalColormap {
    /// Create a categorical colormap with a default 10-color palette.
    pub fn default_palette() -> Self {
        Self {
            colors: vec![
                Color::new(0.122, 0.467, 0.706, 1.0),
                Color::new(1.000, 0.498, 0.055, 1.0),
                Color::new(0.173, 0.627, 0.173, 1.0),
                Color::new(0.839, 0.153, 0.157, 1.0),
                Color::new(0.580, 0.404, 0.741, 1.0),
                Color::new(0.549, 0.337, 0.294, 1.0),
                Color::new(0.890, 0.467, 0.761, 1.0),
                Color::new(0.498, 0.498, 0.498, 1.0),
                Color::new(0.737, 0.741, 0.133, 1.0),
                Color::new(0.090, 0.745, 0.812, 1.0),
            ],
        }
    }
    /// Create a categorical colormap from custom colors.
    pub fn from_colors(colors: Vec<Color>) -> Self {
        Self { colors }
    }
    /// Map a category index to a color.
    ///
    /// Wraps around if the index exceeds the palette size.
    pub fn map_category(&self, index: usize) -> Color {
        if self.colors.is_empty() {
            return Color::new(0.0, 0.0, 0.0, 1.0);
        }
        self.colors[index % self.colors.len()]
    }
    /// Map a slice of category indices to colors.
    pub fn map_categories(&self, indices: &[usize]) -> Vec<Color> {
        indices.iter().map(|&i| self.map_category(i)).collect()
    }
    /// Number of distinct colors in the palette.
    pub fn palette_size(&self) -> usize {
        self.colors.len()
    }
}
/// Style of cyclic (wrapping) colormap.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CyclicStyle {
    /// Full HSV hue wheel: hue rotates 0 → 360° and wraps.
    HsvWheel,
    /// Twilight-inspired: blue → white → red → black → blue.
    Twilight,
    /// Phase: cyan → yellow → cyan (good for angle data).
    Phase,
}
/// A colormap that uses Catmull-Rom spline interpolation between color stops.
///
/// For fewer than four stops the implementation falls back to linear blending
/// so the API remains consistent regardless of the number of stops.
#[derive(Debug, Clone)]
pub struct SplineColormap {
    pub(super) stops: Vec<ColorStop>,
}
impl SplineColormap {
    /// Create a spline colormap from a list of `(t, color)` pairs.
    pub fn new(stops: impl IntoIterator<Item = (f64, Color)>) -> Self {
        let mut v: Vec<ColorStop> = stops
            .into_iter()
            .map(|(t, color)| ColorStop {
                t: t.clamp(0.0, 1.0),
                color,
            })
            .collect();
        v.sort_by(|a, b| a.t.partial_cmp(&b.t).unwrap_or(std::cmp::Ordering::Equal));
        Self { stops: v }
    }
    /// Sample the spline colormap at `t` in \[0, 1\].
    pub fn sample(&self, t: f64) -> Color {
        let t = t.clamp(0.0, 1.0);
        let n = self.stops.len();
        if n == 0 {
            return Color {
                r: 0.0,
                g: 0.0,
                b: 0.0,
                a: 1.0,
            };
        }
        if n == 1 {
            return self.stops[0].color;
        }
        if n < 4 {
            let idx = self.stops.partition_point(|s| s.t <= t).saturating_sub(1);
            let lo = &self.stops[idx];
            let hi = &self.stops[(idx + 1).min(n - 1)];
            let span = hi.t - lo.t;
            let alpha = if span < 1e-12 { 0.0 } else { (t - lo.t) / span };
            return lerp_color(lo.color, hi.color, alpha as f32);
        }
        let idx = self
            .stops
            .partition_point(|s| s.t <= t)
            .saturating_sub(1)
            .min(n - 2);
        let i0 = if idx == 0 { 0 } else { idx - 1 };
        let i1 = idx;
        let i2 = (idx + 1).min(n - 1);
        let i3 = (idx + 2).min(n - 1);
        let span = self.stops[i2].t - self.stops[i1].t;
        let u = if span < 1e-12 {
            0.0
        } else {
            (t - self.stops[i1].t) / span
        };
        let u = u.clamp(0.0, 1.0) as f32;
        catmull_rom_color(
            self.stops[i0].color,
            self.stops[i1].color,
            self.stops[i2].color,
            self.stops[i3].color,
            u,
        )
    }
    /// Build a discrete palette of `n` colors.
    pub fn build(&self, n: usize) -> Vec<Color> {
        let n = n.max(1);
        (0..n)
            .map(|i| {
                let t = if n == 1 {
                    0.5
                } else {
                    i as f64 / (n - 1) as f64
                };
                self.sample(t)
            })
            .collect()
    }
}
/// Available colormaps for scalar visualization.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Colormap {
    /// Jet colormap (blue -> cyan -> green -> yellow -> red).
    Jet,
    /// Viridis colormap (dark purple -> teal -> yellow).
    Viridis,
    /// Cool-warm diverging colormap (blue -> white -> red).
    CoolWarm,
    /// Grayscale (black -> white).
    Grayscale,
    /// Plasma colormap.
    Plasma,
    /// Inferno colormap.
    Inferno,
    /// Magma colormap.
    Magma,
    /// Turbo colormap.
    Turbo,
    /// Red-Blue diverging colormap.
    RdBu,
    /// Seismic diverging colormap.
    Seismic,
    /// Brown-Blue-Green diverging colormap.
    BrBG,
    /// Cividis colormap (perceptually uniform, deuteranopia-friendly).
    Cividis,
    /// Isoluminant rainbow colormap.
    IsoRainbow,
}
impl Colormap {
    /// Returns whether this colormap is perceptually uniform.
    pub fn is_perceptually_uniform(self) -> bool {
        matches!(
            self,
            Colormap::Viridis
                | Colormap::Plasma
                | Colormap::Inferno
                | Colormap::Magma
                | Colormap::Turbo
                | Colormap::Cividis
        )
    }
    /// Returns whether this colormap is diverging.
    pub fn is_diverging(self) -> bool {
        matches!(
            self,
            Colormap::CoolWarm | Colormap::RdBu | Colormap::Seismic | Colormap::BrBG
        )
    }
    /// Compute the CIELab L* lightness values for `n` evenly-spaced samples of
    /// this colormap.
    ///
    /// Each returned value is the L* component (0 = black, 100 = white) of the
    /// colormap colour at the corresponding normalised position in `[0, 1]`.
    /// Returns an empty vector when `n == 0`.
    pub fn compute_perceptual_lightness(self, n: usize) -> Vec<f64> {
        if n == 0 {
            return Vec::new();
        }
        (0..n)
            .map(|i| {
                let t = if n == 1 {
                    0.5
                } else {
                    i as f64 / (n - 1) as f64
                };
                let c = map_scalar(t, 0.0, 1.0, self);
                let (l, _, _) = rgb_to_lab(c.r as f64, c.g as f64, c.b as f64);
                l
            })
            .collect()
    }
    /// Compute the CIEDE2000 color difference between each consecutive pair of
    /// `n` evenly-spaced colormap samples.
    ///
    /// Returns a vector of length `n - 1`.  A perceptually uniform colormap
    /// should have approximately equal values across all steps.  Returns an
    /// empty vector when `n < 2`.
    pub fn compute_color_difference(self, n: usize) -> Vec<f64> {
        if n < 2 {
            return Vec::new();
        }
        let samples: Vec<Color> = (0..n)
            .map(|i| {
                let t = i as f64 / (n - 1) as f64;
                map_scalar(t, 0.0, 1.0, self)
            })
            .collect();
        samples.windows(2).map(|w| ciede2000(w[0], w[1])).collect()
    }
    /// Generate a categorical palette of `n` visually distinct colors.
    ///
    /// Colors are selected by spacing hues uniformly around the HSV color
    /// wheel at maximum saturation and a mid-range value, ensuring sufficient
    /// visual separation for up to ~20 categories.  For larger `n` the palette
    /// wraps with adjusted lightness levels to maintain distinctiveness.
    pub fn generate_categorical_palette(n: usize) -> Vec<Color> {
        if n == 0 {
            return Vec::new();
        }
        let mut out = Vec::with_capacity(n);
        for i in 0..n {
            let pass = i / n.max(1);
            let step = if pass == 0 { n } else { n / 2 };
            let step = step.max(1);
            let h = ((i * 360 / step) % 360) as f32;
            let s = 0.85_f32;
            let v = if pass == 0 { 0.90_f32 } else { 0.65_f32 };
            let (r, g, b) = hsv_to_rgb_f32(h, s, v);
            out.push(Color::new(r, g, b, 1.0));
        }
        out
    }
}
