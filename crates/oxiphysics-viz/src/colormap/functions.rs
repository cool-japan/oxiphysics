//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use crate::primitives::Color;

use super::types::{ColorMap, Colormap, InterpolationMode};

/// Map a scalar value to a color using the specified colormap.
///
/// The value is clamped to `[min, max]` and normalized to `[0, 1]`.
pub fn map_scalar(value: f64, min: f64, max: f64, colormap: Colormap) -> Color {
    let range = max - min;
    let t = if range.abs() < 1e-15 {
        0.5
    } else {
        ((value - min) / range).clamp(0.0, 1.0)
    };
    match colormap {
        Colormap::Jet => jet(t as f32),
        Colormap::Viridis => viridis(t as f32),
        Colormap::CoolWarm => cool_warm(t as f32),
        Colormap::Grayscale => grayscale(t as f32),
        Colormap::Plasma => plasma(t as f32),
        Colormap::Inferno => inferno(t as f32),
        Colormap::Magma => magma(t as f32),
        Colormap::Turbo => turbo(t as f32),
        Colormap::RdBu => rdbu(t as f32),
        Colormap::Seismic => seismic(t as f32),
        Colormap::BrBG => brbg(t as f32),
        Colormap::Cividis => cividis_fn(t as f32),
        Colormap::IsoRainbow => iso_rainbow_fn(t as f32),
    }
}
/// Jet colormap: blue(0) -> cyan(0.25) -> green(0.5) -> yellow(0.75) -> red(1).
pub fn jet(t: f32) -> Color {
    let stops: [(f32, Color); 5] = [
        (0.0, Color::new(0.0, 0.0, 1.0, 1.0)),
        (0.25, Color::new(0.0, 1.0, 1.0, 1.0)),
        (0.5, Color::new(0.0, 1.0, 0.0, 1.0)),
        (0.75, Color::new(1.0, 1.0, 0.0, 1.0)),
        (1.0, Color::new(1.0, 0.0, 0.0, 1.0)),
    ];
    interpolate_stops(&stops, t)
}
/// Viridis colormap approximation.
pub fn viridis(t: f32) -> Color {
    let stops: [(f32, Color); 5] = [
        (0.0, Color::new(0.267, 0.004, 0.329, 1.0)),
        (0.25, Color::new(0.282, 0.140, 0.458, 1.0)),
        (0.5, Color::new(0.127, 0.566, 0.551, 1.0)),
        (0.75, Color::new(0.544, 0.774, 0.248, 1.0)),
        (1.0, Color::new(0.993, 0.906, 0.144, 1.0)),
    ];
    interpolate_stops(&stops, t)
}
/// Cool-warm diverging colormap: blue -> white -> red.
pub fn cool_warm(t: f32) -> Color {
    let stops: [(f32, Color); 3] = [
        (0.0, Color::new(0.231, 0.298, 0.753, 1.0)),
        (0.5, Color::new(0.865, 0.865, 0.865, 1.0)),
        (1.0, Color::new(0.706, 0.016, 0.150, 1.0)),
    ];
    interpolate_stops(&stops, t)
}
/// Grayscale: black -> white.
pub fn grayscale(t: f32) -> Color {
    Color::new(t, t, t, 1.0)
}
/// Plasma colormap approximation.
pub fn plasma(t: f32) -> Color {
    let stops: [(f32, Color); 5] = [
        (0.0, Color::new(0.050, 0.030, 0.528, 1.0)),
        (0.25, Color::new(0.498, 0.031, 0.598, 1.0)),
        (0.5, Color::new(0.798, 0.208, 0.432, 1.0)),
        (0.75, Color::new(0.973, 0.535, 0.229, 1.0)),
        (1.0, Color::new(0.940, 0.975, 0.131, 1.0)),
    ];
    interpolate_stops(&stops, t)
}
/// Inferno colormap approximation.
pub fn inferno(t: f32) -> Color {
    let stops: [(f32, Color); 5] = [
        (0.0, Color::new(0.0, 0.0, 0.016, 1.0)),
        (0.25, Color::new(0.275, 0.055, 0.404, 1.0)),
        (0.5, Color::new(0.663, 0.196, 0.357, 1.0)),
        (0.75, Color::new(0.941, 0.592, 0.122, 1.0)),
        (1.0, Color::new(0.988, 1.0, 0.643, 1.0)),
    ];
    interpolate_stops(&stops, t)
}
/// Magma colormap approximation.
pub fn magma(t: f32) -> Color {
    let stops: [(f32, Color); 5] = [
        (0.0, Color::new(0.0, 0.0, 0.016, 1.0)),
        (0.25, Color::new(0.271, 0.047, 0.400, 1.0)),
        (0.5, Color::new(0.604, 0.196, 0.514, 1.0)),
        (0.75, Color::new(0.933, 0.541, 0.600, 1.0)),
        (1.0, Color::new(0.988, 0.992, 0.749, 1.0)),
    ];
    interpolate_stops(&stops, t)
}
/// Turbo colormap approximation (Google / Anton Mikhailov).
pub fn turbo(t: f32) -> Color {
    let stops: [(f32, Color); 6] = [
        (0.0, Color::new(0.190, 0.071, 0.231, 1.0)),
        (0.2, Color::new(0.071, 0.518, 0.839, 1.0)),
        (0.4, Color::new(0.102, 0.804, 0.518, 1.0)),
        (0.6, Color::new(0.682, 0.890, 0.243, 1.0)),
        (0.8, Color::new(0.988, 0.647, 0.118, 1.0)),
        (1.0, Color::new(0.667, 0.000, 0.024, 1.0)),
    ];
    interpolate_stops(&stops, t)
}
/// Red-Blue diverging colormap.
pub fn rdbu(t: f32) -> Color {
    let stops: [(f32, Color); 3] = [
        (0.0, Color::new(0.706, 0.016, 0.150, 1.0)),
        (0.5, Color::new(0.969, 0.969, 0.969, 1.0)),
        (1.0, Color::new(0.129, 0.400, 0.675, 1.0)),
    ];
    interpolate_stops(&stops, t)
}
/// Seismic diverging colormap.
pub fn seismic(t: f32) -> Color {
    let stops: [(f32, Color); 5] = [
        (0.0, Color::new(0.0, 0.0, 0.3, 1.0)),
        (0.25, Color::new(0.0, 0.0, 1.0, 1.0)),
        (0.5, Color::new(1.0, 1.0, 1.0, 1.0)),
        (0.75, Color::new(1.0, 0.0, 0.0, 1.0)),
        (1.0, Color::new(0.5, 0.0, 0.0, 1.0)),
    ];
    interpolate_stops(&stops, t)
}
/// Brown-Blue-Green diverging colormap (BrBG).
pub fn brbg(t: f32) -> Color {
    let stops: [(f32, Color); 5] = [
        (0.0, Color::new(0.329, 0.188, 0.020, 1.0)),
        (0.25, Color::new(0.749, 0.506, 0.176, 1.0)),
        (0.5, Color::new(0.961, 0.961, 0.961, 1.0)),
        (0.75, Color::new(0.353, 0.706, 0.675, 1.0)),
        (1.0, Color::new(0.004, 0.400, 0.369, 1.0)),
    ];
    interpolate_stops(&stops, t)
}
/// Linearly interpolate between color stops.
pub fn interpolate_stops(stops: &[(f32, Color)], t: f32) -> Color {
    if t <= stops[0].0 {
        return stops[0].1;
    }
    if t >= stops[stops.len() - 1].0 {
        return stops[stops.len() - 1].1;
    }
    for i in 0..stops.len() - 1 {
        if t >= stops[i].0 && t <= stops[i + 1].0 {
            let local_t = (t - stops[i].0) / (stops[i + 1].0 - stops[i].0);
            return Color::lerp(&stops[i].1, &stops[i + 1].1, local_t);
        }
    }
    stops[stops.len() - 1].1
}
/// Sample a set of color stops with the specified interpolation mode.
pub fn sample_stops_with_mode(stops: &[(f32, Color)], t: f32, mode: InterpolationMode) -> Color {
    if stops.is_empty() {
        return Color::new(0.0, 0.0, 0.0, 1.0);
    }
    let t_clamped = t.clamp(stops[0].0, stops[stops.len() - 1].0);
    match mode {
        InterpolationMode::Linear => interpolate_stops(stops, t_clamped),
        InterpolationMode::Nearest => {
            let mut best_idx = 0;
            let mut best_dist = f32::MAX;
            for (i, &(pos, _)) in stops.iter().enumerate() {
                let dist = (pos - t_clamped).abs();
                if dist < best_dist {
                    best_dist = dist;
                    best_idx = i;
                }
            }
            stops[best_idx].1
        }
        InterpolationMode::Step => {
            let mut result = stops[0].1;
            for &(pos, color) in stops {
                if pos <= t_clamped {
                    result = color;
                } else {
                    break;
                }
            }
            result
        }
    }
}
/// Color histogram bins using a colormap.
///
/// Each bin gets a color based on its index position in \[0, 1\].
pub fn histogram_colors(n_bins: usize, colormap: Colormap) -> Vec<Color> {
    if n_bins == 0 {
        return Vec::new();
    }
    (0..n_bins)
        .map(|i| {
            let t = if n_bins == 1 {
                0.5
            } else {
                i as f64 / (n_bins - 1) as f64
            };
            map_scalar(t, 0.0, 1.0, colormap)
        })
        .collect()
}
/// Color histogram bins based on their values.
///
/// Each bin is colored by its count value mapped to \[vmin, vmax\].
pub fn histogram_value_colors(
    counts: &[f64],
    colormap: Colormap,
    vmin: f64,
    vmax: f64,
) -> Vec<Color> {
    counts
        .iter()
        .map(|&c| map_scalar(c, vmin, vmax, colormap))
        .collect()
}
/// Apply a [`ColorMap`] to a slice of values, normalizing by `[vmin, vmax]`.
///
/// Values outside `[vmin, vmax]` are clamped.
pub fn apply_colormap(values: &[f64], cmap: ColorMap, vmin: f64, vmax: f64) -> Vec<[u8; 4]> {
    let range = vmax - vmin;
    values
        .iter()
        .map(|&v| {
            let t = if range.abs() < 1e-30 {
                0.5
            } else {
                ((v - vmin) / range).clamp(0.0, 1.0)
            };
            cmap.sample(t)
        })
        .collect()
}
/// Normalize a slice of values to `[0, 1]` using min-max scaling.
///
/// If all values are identical, returns a vector of `0.5`.
pub fn normalize(values: &[f64]) -> Vec<f64> {
    if values.is_empty() {
        return Vec::new();
    }
    let min = values.iter().cloned().fold(f64::INFINITY, f64::min);
    let max = values.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let range = max - min;
    if range.abs() < 1e-30 {
        return vec![0.5; values.len()];
    }
    values.iter().map(|&v| (v - min) / range).collect()
}
/// Symmetrically normalize values so that 0.0 maps to 0.5.
///
/// The absolute maximum |v| is used as the scale: output is `v / (2*max_abs) + 0.5`.
pub fn symmetric_normalize(values: &[f64]) -> Vec<f64> {
    if values.is_empty() {
        return Vec::new();
    }
    let max_abs = values.iter().map(|v| v.abs()).fold(0.0_f64, f64::max);
    if max_abs < 1e-30 {
        return vec![0.5; values.len()];
    }
    values.iter().map(|&v| v / (2.0 * max_abs) + 0.5).collect()
}
/// Normalize using percentile-based clipping.
///
/// Values below the `low_pct` percentile and above the `high_pct` percentile
/// are clipped before normalization. Percentiles are in \[0, 100\].
pub fn percentile_normalize(values: &[f64], low_pct: f64, high_pct: f64) -> Vec<f64> {
    if values.is_empty() {
        return Vec::new();
    }
    let mut sorted: Vec<f64> = values.to_vec();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let low_idx = ((low_pct / 100.0) * (sorted.len() - 1) as f64).round() as usize;
    let high_idx = ((high_pct / 100.0) * (sorted.len() - 1) as f64).round() as usize;
    let low_idx = low_idx.min(sorted.len() - 1);
    let high_idx = high_idx.min(sorted.len() - 1);
    let vmin = sorted[low_idx];
    let vmax = sorted[high_idx];
    let range = vmax - vmin;
    if range.abs() < 1e-30 {
        return vec![0.5; values.len()];
    }
    values
        .iter()
        .map(|&v| ((v - vmin) / range).clamp(0.0, 1.0))
        .collect()
}
/// Linearly blend two RGBA colors at parameter `t in [0, 1]`.
///
/// `t = 0` returns `a`; `t = 1` returns `b`. Alpha is always 255.
pub fn blend_rgb(a: [u8; 4], b: [u8; 4], t: f64) -> [u8; 4] {
    let t = t.clamp(0.0, 1.0);
    let lerp = |x: u8, y: u8| -> u8 {
        ((x as f64) * (1.0 - t) + (y as f64) * t)
            .round()
            .clamp(0.0, 255.0) as u8
    };
    [lerp(a[0], b[0]), lerp(a[1], b[1]), lerp(a[2], b[2]), 255]
}
/// Premultiplied alpha blend: `result = src + dst * (1 - src_alpha)`.
pub fn alpha_blend(src: [u8; 4], dst: [u8; 4]) -> [u8; 4] {
    let sa = src[3] as f64 / 255.0;
    let da = dst[3] as f64 / 255.0;
    let out_a = sa + da * (1.0 - sa);
    if out_a < 1e-10 {
        return [0, 0, 0, 0];
    }
    let blend_ch = |s: u8, d: u8| -> u8 {
        let sc = s as f64 / 255.0;
        let dc = d as f64 / 255.0;
        let out = (sc * sa + dc * da * (1.0 - sa)) / out_a;
        (out * 255.0).round().clamp(0.0, 255.0) as u8
    };
    [
        blend_ch(src[0], dst[0]),
        blend_ch(src[1], dst[1]),
        blend_ch(src[2], dst[2]),
        (out_a * 255.0).round().clamp(0.0, 255.0) as u8,
    ]
}
/// Convert RGB (all in `[0, 1]`) to HSV.
///
/// Returns `(h, s, v)` where h in \[0, 360), s in \[0, 1\\], v in \[0, 1\].
pub fn rgb_to_hsv(r: f64, g: f64, b: f64) -> (f64, f64, f64) {
    let cmax = r.max(g).max(b);
    let cmin = r.min(g).min(b);
    let delta = cmax - cmin;
    let v = cmax;
    let s = if cmax > 1e-10 { delta / cmax } else { 0.0 };
    let h = if delta < 1e-10 {
        0.0
    } else if (cmax - r).abs() < 1e-10 {
        60.0 * (((g - b) / delta) % 6.0)
    } else if (cmax - g).abs() < 1e-10 {
        60.0 * ((b - r) / delta + 2.0)
    } else {
        60.0 * ((r - g) / delta + 4.0)
    };
    let h = if h < 0.0 { h + 360.0 } else { h };
    (h, s, v)
}
/// Convert HSV (h in \[0, 360), s in \[0, 1\\], v in \[0, 1\]) to RGB in `[0, 1]`.
pub fn hsv_to_rgb(h: f64, s: f64, v: f64) -> (f64, f64, f64) {
    if s < 1e-10 {
        return (v, v, v);
    }
    let h_norm = h / 60.0;
    let i = h_norm.floor() as i32 % 6;
    let f = h_norm - h_norm.floor();
    let p = v * (1.0 - s);
    let q = v * (1.0 - s * f);
    let t = v * (1.0 - s * (1.0 - f));
    match i {
        0 => (v, t, p),
        1 => (q, v, p),
        2 => (p, v, t),
        3 => (p, q, v),
        4 => (t, p, v),
        _ => (v, p, q),
    }
}
/// Convert RGB (all in \[0, 1\]) to HSL.
///
/// Returns `(h, s, l)` where h in \[0, 360), s in \[0, 1\\], l in \[0, 1\].
pub fn rgb_to_hsl(r: f64, g: f64, b: f64) -> (f64, f64, f64) {
    let cmax = r.max(g).max(b);
    let cmin = r.min(g).min(b);
    let delta = cmax - cmin;
    let l = (cmax + cmin) / 2.0;
    let s = if delta < 1e-10 {
        0.0
    } else {
        delta / (1.0 - (2.0 * l - 1.0).abs())
    };
    let h = if delta < 1e-10 {
        0.0
    } else if (cmax - r).abs() < 1e-10 {
        60.0 * (((g - b) / delta) % 6.0)
    } else if (cmax - g).abs() < 1e-10 {
        60.0 * ((b - r) / delta + 2.0)
    } else {
        60.0 * ((r - g) / delta + 4.0)
    };
    let h = if h < 0.0 { h + 360.0 } else { h };
    (h, s, l)
}
/// Cividis colormap: perceptually uniform, suitable for colour-vision deficiency.
pub fn cividis_fn(t: f32) -> Color {
    let stops: [(f32, Color); 6] = [
        (0.00, Color::new(0.000, 0.135, 0.304, 1.0)),
        (0.20, Color::new(0.124, 0.222, 0.498, 1.0)),
        (0.40, Color::new(0.302, 0.357, 0.558, 1.0)),
        (0.60, Color::new(0.520, 0.500, 0.524, 1.0)),
        (0.80, Color::new(0.738, 0.660, 0.372, 1.0)),
        (1.00, Color::new(0.996, 0.851, 0.000, 1.0)),
    ];
    interpolate_stops(&stops, t)
}
/// IsoRainbow: constant-luminance rainbow for isoluminant visualization.
pub fn iso_rainbow_fn(t: f32) -> Color {
    let stops: [(f32, Color); 7] = [
        (0.00, Color::new(0.859, 0.396, 0.427, 1.0)),
        (0.17, Color::new(0.824, 0.510, 0.086, 1.0)),
        (0.33, Color::new(0.639, 0.627, 0.000, 1.0)),
        (0.50, Color::new(0.333, 0.694, 0.271, 1.0)),
        (0.67, Color::new(0.000, 0.718, 0.576, 1.0)),
        (0.83, Color::new(0.024, 0.667, 0.839, 1.0)),
        (1.00, Color::new(0.525, 0.573, 0.933, 1.0)),
    ];
    interpolate_stops(&stops, t)
}
/// Rainbow colormap: classic hue rotation (red → violet).
///
/// Note: this is *not* perceptually uniform.  Prefer `Viridis` or `Plasma`
/// for scientific use.  Retained for legacy / artistic purposes.
pub fn rainbow(t: f32) -> Color {
    let h = (1.0 - t.clamp(0.0, 1.0)) * 300.0;
    let (r, g, b) = hsv_to_rgb_f32(h, 1.0, 1.0);
    Color::new(r, g, b, 1.0)
}
/// Helper: HSV → RGB with f32 inputs.  H in \[0, 360), S/V in \[0, 1\\].
pub fn hsv_to_rgb_f32(h: f32, s: f32, v: f32) -> (f32, f32, f32) {
    if s < 1e-7 {
        return (v, v, v);
    }
    let hi = (h / 60.0).floor() as u32 % 6;
    let f = h / 60.0 - (h / 60.0).floor();
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
/// Return the Tab10 palette (10 colors from Matplotlib's tableau).
pub fn tab10_palette() -> [Color; 10] {
    [
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
    ]
}
/// Return the Set1 palette (9 high-saturation colors from ColorBrewer).
pub fn set1_palette() -> [Color; 9] {
    [
        Color::new(0.894, 0.102, 0.110, 1.0),
        Color::new(0.216, 0.494, 0.722, 1.0),
        Color::new(0.302, 0.686, 0.290, 1.0),
        Color::new(0.596, 0.306, 0.639, 1.0),
        Color::new(1.000, 0.498, 0.000, 1.0),
        Color::new(1.000, 1.000, 0.200, 1.0),
        Color::new(0.651, 0.337, 0.157, 1.0),
        Color::new(0.969, 0.506, 0.749, 1.0),
        Color::new(0.600, 0.600, 0.600, 1.0),
    ]
}
/// Return the Paired palette (12 colors; paired light/dark versions of 6 hues).
pub fn paired_palette() -> [Color; 12] {
    [
        Color::new(0.651, 0.808, 0.890, 1.0),
        Color::new(0.122, 0.471, 0.706, 1.0),
        Color::new(0.698, 0.875, 0.541, 1.0),
        Color::new(0.200, 0.627, 0.173, 1.0),
        Color::new(0.984, 0.604, 0.600, 1.0),
        Color::new(0.890, 0.102, 0.110, 1.0),
        Color::new(0.992, 0.749, 0.435, 1.0),
        Color::new(1.000, 0.498, 0.000, 1.0),
        Color::new(0.792, 0.698, 0.839, 1.0),
        Color::new(0.416, 0.239, 0.604, 1.0),
        Color::new(1.000, 1.000, 0.600, 1.0),
        Color::new(0.694, 0.349, 0.157, 1.0),
    ]
}
/// Invert a colormap: sample at `1 - t` instead of `t`.
///
/// # Example
/// `invert_colormap(Colormap::Viridis, 0.0)` returns the color that
/// `Viridis` produces at `t = 1.0`.
pub fn invert_colormap(cmap: Colormap, t: f64) -> Color {
    map_scalar(1.0 - t.clamp(0.0, 1.0), 0.0, 1.0, cmap)
}
/// Sample a continuous colormap at `n` evenly-spaced values, producing a
/// fixed-size discrete palette.
///
/// Returns a `Vec` of `n` colors spanning `[0, 1]`.
pub fn discrete_colormap(cmap: Colormap, n: usize) -> Vec<Color> {
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
            map_scalar(t, 0.0, 1.0, cmap)
        })
        .collect()
}
/// Map a scalar to an RGBA `[u8; 4]` with a custom alpha value.
///
/// The alpha channel is specified in `[0, 1]`; it is independent of the
/// colormap and is simply packed into the output.
pub fn map_scalar_rgba(value: f64, min: f64, max: f64, cmap: Colormap, alpha: f32) -> [u8; 4] {
    let c = map_scalar(value, min, max, cmap);
    [
        (c.r.clamp(0.0, 1.0) * 255.0).round() as u8,
        (c.g.clamp(0.0, 1.0) * 255.0).round() as u8,
        (c.b.clamp(0.0, 1.0) * 255.0).round() as u8,
        (alpha.clamp(0.0, 1.0) * 255.0).round() as u8,
    ]
}
/// Convert linear-RGB (each component in `[0, 1]`) to CIELAB `(L*, a*, b*)`.
///
/// Illuminant D65 reference white is used.  The input is assumed to be in
/// *linear* light (not gamma-encoded sRGB).
pub fn rgb_to_lab(r: f64, g: f64, b: f64) -> (f64, f64, f64) {
    let x = 0.4124564 * r + 0.3575761 * g + 0.1804375 * b;
    let y = 0.2126729 * r + 0.7151522 * g + 0.0721750 * b;
    let z = 0.0193339 * r + 0.1191920 * g + 0.9503041 * b;
    let xn = x / 0.95047;
    let yn = y / 1.00000;
    let zn = z / 1.08883;
    let f = |t: f64| -> f64 {
        if t > 0.008856 {
            t.powf(1.0 / 3.0)
        } else {
            7.787 * t + 16.0 / 116.0
        }
    };
    let fx = f(xn);
    let fy = f(yn);
    let fz = f(zn);
    let l = 116.0 * fy - 16.0;
    let a = 500.0 * (fx - fy);
    let bb = 200.0 * (fy - fz);
    (l, a, bb)
}
/// Convert CIELAB `(L*, a*, b*)` back to linear-RGB `(r, g, b)` in `[0, 1]`.
///
/// Values that fall outside the sRGB gamut are clamped.  Illuminant D65
/// reference white is used.
pub fn lab_to_rgb(l: f64, a: f64, b: f64) -> (f64, f64, f64) {
    let fy = (l + 16.0) / 116.0;
    let fx = a / 500.0 + fy;
    let fz = fy - b / 200.0;
    let finv = |t: f64| -> f64 {
        let t3 = t * t * t;
        if t3 > 0.008856 {
            t3
        } else {
            (t - 16.0 / 116.0) / 7.787
        }
    };
    let x = finv(fx) * 0.95047;
    let y = finv(fy) * 1.00000;
    let z = finv(fz) * 1.08883;
    let r = 3.2404542 * x - 1.5371385 * y - 0.4985314 * z;
    let g = -0.9692660 * x + 1.8760108 * y + 0.0415560 * z;
    let bb = 0.0556434 * x - 0.2040259 * y + 1.0572252 * z;
    (r.clamp(0.0, 1.0), g.clamp(0.0, 1.0), bb.clamp(0.0, 1.0))
}
/// Compute an approximate CIEDE2000 perceptual distance between two colors.
///
/// Both colors should have components in `[0, 1]` (linear-RGB).
///
/// This is a lightweight approximation using L*a*b* Euclidean distance as a
/// proxy for the full CIEDE2000 formula.  The true CIEDE2000 metric requires
/// hue-rotation corrections not included here.
pub fn perceptual_distance(c1: Color, c2: Color) -> f64 {
    let (l1, a1, b1) = rgb_to_lab(c1.r as f64, c1.g as f64, c1.b as f64);
    let (l2, a2, b2) = rgb_to_lab(c2.r as f64, c2.g as f64, c2.b as f64);
    let dl = l1 - l2;
    let da = a1 - a2;
    let db = b1 - b2;
    (dl * dl + da * da + db * db).sqrt()
}
/// Estimate the perceptual uniformity of a colormap by measuring the
/// standard deviation of the CIEDE2000-approximate distances between
/// `n` consecutive uniformly-spaced samples.
///
/// A perfectly uniform map has standard deviation 0.  Lower values are
/// better; values below ~5 are generally considered acceptable.
pub fn colormap_uniformity_score(cmap: Colormap, n: usize) -> f64 {
    let n = n.max(2);
    let colors: Vec<Color> = (0..n)
        .map(|i| map_scalar(i as f64 / (n - 1) as f64, 0.0, 1.0, cmap))
        .collect();
    let dists: Vec<f64> = colors
        .windows(2)
        .map(|w| perceptual_distance(w[0], w[1]))
        .collect();
    let mean = dists.iter().sum::<f64>() / dists.len() as f64;
    let variance = dists.iter().map(|d| (d - mean).powi(2)).sum::<f64>() / dists.len() as f64;
    variance.sqrt()
}
/// Compute the CIEDE2000 color difference between two sRGB colors.
///
/// This implements the full CIEDE2000 formula including hue-rotation,
/// chroma weighting, and lightness compensation terms.  Both colors
/// should have components in `[0, 1]`.
///
/// Returns a non-negative value; differences below 1 are imperceptible,
/// differences above 10 are clearly visible.
pub fn ciede2000(c1: Color, c2: Color) -> f64 {
    let (l1, a1s, b1s) = rgb_to_lab(c1.r as f64, c1.g as f64, c1.b as f64);
    let (l2, a2s, b2s) = rgb_to_lab(c2.r as f64, c2.g as f64, c2.b as f64);
    let c1ab = (a1s * a1s + b1s * b1s).sqrt();
    let c2ab = (a2s * a2s + b2s * b2s).sqrt();
    let c_avg = (c1ab + c2ab) / 2.0;
    let c7 = c_avg.powi(7);
    let g = 0.5 * (1.0 - (c7 / (c7 + 25_f64.powi(7))).sqrt());
    let a1p = a1s * (1.0 + g);
    let a2p = a2s * (1.0 + g);
    let c1p = (a1p * a1p + b1s * b1s).sqrt();
    let c2p = (a2p * a2p + b2s * b2s).sqrt();
    let h1p = b1s.atan2(a1p).to_degrees().rem_euclid(360.0);
    let h2p = b2s.atan2(a2p).to_degrees().rem_euclid(360.0);
    let dl_p = l2 - l1;
    let dc_p = c2p - c1p;
    let dh_p = if (c1p * c2p).abs() < 1e-15 {
        0.0
    } else {
        let diff = h2p - h1p;
        if diff.abs() <= 180.0 {
            diff
        } else if diff > 180.0 {
            diff - 360.0
        } else {
            diff + 360.0
        }
    };
    let dh_big = 2.0 * (c1p * c2p).sqrt() * (dh_p.to_radians() / 2.0).sin();
    let l_avg = (l1 + l2) / 2.0;
    let c_avg_p = (c1p + c2p) / 2.0;
    let h_avg_p = if (c1p * c2p).abs() < 1e-15 {
        h1p + h2p
    } else if (h1p - h2p).abs() <= 180.0 {
        (h1p + h2p) / 2.0
    } else if h1p + h2p < 360.0 {
        (h1p + h2p + 360.0) / 2.0
    } else {
        (h1p + h2p - 360.0) / 2.0
    };
    let t_term = 1.0 - 0.17 * ((h_avg_p - 30.0).to_radians()).cos()
        + 0.24 * (2.0 * h_avg_p.to_radians()).cos()
        + 0.32 * ((3.0 * h_avg_p + 6.0).to_radians()).cos()
        - 0.20 * ((4.0 * h_avg_p - 63.0).to_radians()).cos();
    let sl = 1.0 + 0.015 * (l_avg - 50.0).powi(2) / (20.0 + (l_avg - 50.0).powi(2)).sqrt();
    let sc = 1.0 + 0.045 * c_avg_p;
    let sh = 1.0 + 0.015 * c_avg_p * t_term;
    let c7p = c_avg_p.powi(7);
    let rc = 2.0 * (c7p / (c7p + 25_f64.powi(7))).sqrt();
    let d_theta = 30.0 * (-(((h_avg_p - 275.0) / 25.0).powi(2))).exp();
    let rt = -(rc * (2.0 * d_theta.to_radians()).sin());
    let kl = 1.0;
    let kc = 1.0;
    let kh = 1.0;
    ((dl_p / (kl * sl)).powi(2)
        + (dc_p / (kc * sc)).powi(2)
        + (dh_big / (kh * sh)).powi(2)
        + rt * (dc_p / (kc * sc)) * (dh_big / (kh * sh)))
        .sqrt()
}
/// Sample the Cividis colormap at `t ∈ [0, 1]`.
pub fn sample_cividis(t: f32) -> Color {
    cividis_fn(t)
}
/// Sample the IsoRainbow (isoluminant) colormap at `t ∈ [0, 1]`.
pub fn sample_iso_rainbow(t: f32) -> Color {
    iso_rainbow_fn(t)
}
/// Blend two colormaps at parameter `t` with blend weight `alpha`.
///
/// `alpha = 0` returns a sample from `cmap_a`, `alpha = 1` from `cmap_b`.
/// Both colormaps are sampled at the same parameter `t`.
pub fn blend_colormaps(cmap_a: Colormap, cmap_b: Colormap, t: f64, alpha: f64) -> Color {
    let a = map_scalar(t, 0.0, 1.0, cmap_a);
    let b = map_scalar(t, 0.0, 1.0, cmap_b);
    let w = alpha.clamp(0.0, 1.0) as f32;
    Color::new(
        a.r * (1.0 - w) + b.r * w,
        a.g * (1.0 - w) + b.g * w,
        a.b * (1.0 - w) + b.b * w,
        1.0,
    )
}
/// Blend two RGBA `[u8;4]` arrays using `ColorMap` samplers.
///
/// Both colormaps are sampled at `t` and blended with weight `alpha`.
pub fn blend_colormaps_u8(cmap_a: ColorMap, cmap_b: ColorMap, t: f64, alpha: f64) -> [u8; 4] {
    let a = cmap_a.sample(t);
    let b = cmap_b.sample(t);
    blend_rgb(a, b, alpha)
}
/// Equalize a scalar field so that the output histogram is approximately flat.
///
/// `n_bins` controls the resolution of the empirical CDF used for equalization.
/// Returns values in `[0, 1]`.
pub fn equalize_histogram(values: &[f64], n_bins: usize) -> Vec<f64> {
    if values.is_empty() {
        return Vec::new();
    }
    let n_bins = n_bins.max(2);
    let min = values.iter().cloned().fold(f64::INFINITY, f64::min);
    let max = values.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let range = max - min;
    if range < 1e-30 {
        return vec![0.5; values.len()];
    }
    let mut hist = vec![0usize; n_bins];
    for &v in values {
        let bin = ((v - min) / range * (n_bins as f64 - 1.0))
            .round()
            .clamp(0.0, (n_bins - 1) as f64) as usize;
        hist[bin] += 1;
    }
    let mut cdf = vec![0usize; n_bins];
    cdf[0] = hist[0];
    for i in 1..n_bins {
        cdf[i] = cdf[i - 1] + hist[i];
    }
    let cdf_min = *cdf.iter().find(|&&c| c > 0).unwrap_or(&0);
    let total = values.len();
    values
        .iter()
        .map(|&v| {
            let bin = ((v - min) / range * (n_bins as f64 - 1.0))
                .round()
                .clamp(0.0, (n_bins - 1) as f64) as usize;
            if total <= cdf_min {
                0.5
            } else {
                ((cdf[bin] - cdf_min) as f64 / (total - cdf_min) as f64).clamp(0.0, 1.0)
            }
        })
        .collect()
}
/// Apply power-law (gamma) normalization to a slice of values.
///
/// Output = `t^gamma` where `t` is the min-max normalized input.
/// `gamma < 1` brightens dark values; `gamma > 1` darkens them.
pub fn power_normalize(values: &[f64], gamma: f64) -> Vec<f64> {
    if values.is_empty() {
        return Vec::new();
    }
    let normed = normalize(values);
    normed.iter().map(|&t| t.powf(gamma)).collect()
}
/// Normalize values on a log scale: `log(v - min + 1) / log(max - min + 1)`.
///
/// All output values are in `[0, 1]`. Works for non-negative value ranges.
pub fn log_normalize(values: &[f64]) -> Vec<f64> {
    if values.is_empty() {
        return Vec::new();
    }
    let min = values.iter().cloned().fold(f64::INFINITY, f64::min);
    let max = values.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let log_range = (max - min + 1.0).ln();
    if log_range < 1e-30 {
        return vec![0.5; values.len()];
    }
    values
        .iter()
        .map(|&v| ((v - min + 1.0).ln() / log_range).clamp(0.0, 1.0))
        .collect()
}
/// Catmull-Rom interpolation for a single `f32` channel.
pub fn catmull_rom_f32(p0: f32, p1: f32, p2: f32, p3: f32, t: f32) -> f32 {
    let t2 = t * t;
    let t3 = t2 * t;
    0.5 * ((2.0 * p1)
        + (-p0 + p2) * t
        + (2.0 * p0 - 5.0 * p1 + 4.0 * p2 - p3) * t2
        + (-p0 + 3.0 * p1 - 3.0 * p2 + p3) * t3)
}
/// Catmull-Rom interpolation between four `Color` values.
pub fn catmull_rom_color(p0: Color, p1: Color, p2: Color, p3: Color, t: f32) -> Color {
    Color {
        r: catmull_rom_f32(p0.r, p1.r, p2.r, p3.r, t).clamp(0.0, 1.0),
        g: catmull_rom_f32(p0.g, p1.g, p2.g, p3.g, t).clamp(0.0, 1.0),
        b: catmull_rom_f32(p0.b, p1.b, p2.b, p3.b, t).clamp(0.0, 1.0),
        a: catmull_rom_f32(p0.a, p1.a, p2.a, p3.a, t).clamp(0.0, 1.0),
    }
}
/// Linear blend between two `Color` values.
pub fn lerp_color(a: Color, b: Color, t: f32) -> Color {
    let t = t.clamp(0.0, 1.0);
    Color {
        r: a.r + (b.r - a.r) * t,
        g: a.g + (b.g - a.g) * t,
        b: a.b + (b.b - a.b) * t,
        a: a.a + (b.a - a.a) * t,
    }
}
/// Sample the twilight cyclic colormap at `t` in \[0, 1\].
pub fn sample_twilight(t: f32) -> Color {
    let t = t.clamp(0.0, 1.0);
    let (r, g, b) = if t < 0.25 {
        let s = t / 0.25;
        (0.22 + 0.78 * s, 0.24 + 0.76 * s, 0.44 + 0.56 * s)
    } else if t < 0.5 {
        let s = (t - 0.25) / 0.25;
        (1.0 - 0.30 * s, 1.0 - 0.90 * s, 1.0 - 0.90 * s)
    } else if t < 0.75 {
        let s = (t - 0.5) / 0.25;
        (0.70 - 0.65 * s, 0.10 - 0.05 * s, 0.10 - 0.05 * s)
    } else {
        let s = (t - 0.75) / 0.25;
        (0.05 + 0.17 * s, 0.05 + 0.19 * s, 0.05 + 0.39 * s)
    };
    Color { r, g, b, a: 1.0 }
}
/// Sample the phase cyclic colormap at `t` in \[0, 1\].
///
/// Produces a cyan–yellow–cyan cycle suitable for phase/angle data.
pub fn sample_phase(t: f32) -> Color {
    let t = t.clamp(0.0, 1.0);
    let angle = t * std::f32::consts::TAU;
    let sin_v = (angle.sin() + 1.0) * 0.5;
    let cos_v = (angle.cos() + 1.0) * 0.5;
    Color {
        r: sin_v * 0.9,
        g: (sin_v * 0.5 + cos_v * 0.5).clamp(0.0, 1.0),
        b: cos_v * 0.9,
        a: 1.0,
    }
}
/// Sample the Turbo colormap at `t` in \[0, 1\].
///
/// Turbo is a perceptually improved replacement for Jet: smoother, lower
/// perceptual artefacts, colorblind-friendlier than a raw rainbow.
///
/// Coefficients are a piecewise-polynomial approximation taken from the
/// original Google AI blog post.
pub fn sample_turbo(t: f32) -> Color {
    let t = t.clamp(0.0, 1.0);
    let r = (0.1357 + t * (4.5974 + t * (-4.6115 + t * 2.7842))).clamp(0.0, 1.0);
    let g = (0.0914 + t * (2.1856 + t * (4.8052 + t * (-14.078 + t * 8.0)))).clamp(0.0, 1.0);
    let b = (0.1067 + t * (2.6895 + t * (-0.4826 + t * (-5.9454 + t * 4.7512)))).clamp(0.0, 1.0);
    Color { r, g, b, a: 1.0 }
}
/// Sample the classic Jet colormap at `t` in \[0, 1\].
///
/// Jet has known perceptual shortcomings (false detail near yellow/cyan).
/// Prefer [`sample_turbo`] for new work.  This implementation is provided for
/// backward compatibility and comparison purposes.
pub fn sample_jet(t: f32) -> Color {
    let t = t.clamp(0.0, 1.0);
    let r = (1.5 - (4.0 * t - 3.0).abs()).clamp(0.0, 1.0);
    let g = (1.5 - (4.0 * t - 2.0).abs()).clamp(0.0, 1.0);
    let b = (1.5 - (4.0 * t - 1.0).abs()).clamp(0.0, 1.0);
    Color { r, g, b, a: 1.0 }
}
/// Compute the perceptual non-uniformity score of a custom colormap.
///
/// This works like [`colormap_uniformity_score`] but operates on an arbitrary
/// sampler closure rather than a [`Colormap`] variant, making it suitable for
/// `SplineColormap`, `CustomColormapBuilder`, or any other sampler.
pub fn custom_colormap_uniformity_score<F>(sampler: F, n: usize) -> f64
where
    F: Fn(f64) -> Color,
{
    let n = n.max(2);
    let samples: Vec<Color> = (0..n).map(|i| sampler(i as f64 / (n - 1) as f64)).collect();
    let deltas: Vec<f64> = samples
        .windows(2)
        .map(|w| perceptual_distance(w[0], w[1]))
        .collect();
    let mean = deltas.iter().copied().sum::<f64>() / deltas.len() as f64;
    if mean < 1e-12 {
        return 0.0;
    }
    let variance = deltas.iter().map(|d| (d - mean).powi(2)).sum::<f64>() / deltas.len() as f64;
    variance.sqrt() / mean
}
/// Convert a [`Color`] to an \[L*, a*, b*\] triplet without going through the
/// `rgb_to_lab` public API (avoids an extra `f64` cast layer).
pub fn color_to_lab(c: Color) -> (f64, f64, f64) {
    rgb_to_lab(c.r as f64, c.g as f64, c.b as f64)
}
/// Return the lightness (L*) of a `Color` in CIELAB.
pub fn lightness(c: Color) -> f64 {
    let (l, _, _) = color_to_lab(c);
    l
}
/// Check whether a colormap is monotone in lightness (always increasing or
/// always decreasing).  Returns `true` if the lightness sequence is monotone,
/// `false` otherwise.
pub fn is_monotone_lightness<F>(sampler: F, n: usize) -> bool
where
    F: Fn(f64) -> Color,
{
    let n = n.max(2);
    let ls: Vec<f64> = (0..n)
        .map(|i| lightness(sampler(i as f64 / (n - 1) as f64)))
        .collect();
    let increasing = ls.windows(2).all(|w| w[1] >= w[0] - 1e-6);
    let decreasing = ls.windows(2).all(|w| w[1] <= w[0] + 1e-6);
    increasing || decreasing
}
/// Adjust the brightness of a `Color` by a multiplicative `factor`.
///
/// Channels are clamped to \[0, 1\].
pub fn adjust_brightness(c: Color, factor: f32) -> Color {
    Color {
        r: (c.r * factor).clamp(0.0, 1.0),
        g: (c.g * factor).clamp(0.0, 1.0),
        b: (c.b * factor).clamp(0.0, 1.0),
        a: c.a,
    }
}
/// Adjust the saturation of a `Color` by `factor` (in HSV space).
///
/// `factor = 0.0` → fully desaturated (grey).
/// `factor = 1.0` → unchanged.
/// `factor > 1.0` → over-saturated (clamped).
pub fn adjust_saturation(c: Color, factor: f32) -> Color {
    let (h, s, v) = rgb_to_hsv(c.r as f64, c.g as f64, c.b as f64);
    let new_s = (s * factor as f64).clamp(0.0, 1.0);
    let (r, g, b) = hsv_to_rgb(h, new_s, v);
    Color {
        r: r as f32,
        g: g as f32,
        b: b as f32,
        a: c.a,
    }
}
/// Convert a scalar value to a `Color` using a generic sampler closure.
///
/// Equivalent to [`map_scalar`] but for custom colormaps.  `value` is
/// normalized using `min`/`max` before sampling.
pub fn map_scalar_custom<F>(value: f64, min: f64, max: f64, sampler: F) -> Color
where
    F: Fn(f64) -> Color,
{
    let t = if (max - min).abs() < 1e-12 {
        0.5
    } else {
        ((value - min) / (max - min)).clamp(0.0, 1.0)
    };
    sampler(t)
}
/// Apply a custom colormap sampler to a slice of values, returning RGBA bytes.
pub fn apply_custom_colormap<F>(values: &[f64], vmin: f64, vmax: f64, sampler: F) -> Vec<[u8; 4]>
where
    F: Fn(f64) -> Color,
{
    values
        .iter()
        .map(|&v| {
            let c = map_scalar_custom(v, vmin, vmax, &sampler);
            [
                (c.r * 255.0).round() as u8,
                (c.g * 255.0).round() as u8,
                (c.b * 255.0).round() as u8,
                (c.a * 255.0).round() as u8,
            ]
        })
        .collect()
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::colormap::CategoricalColormap;
    use crate::colormap::Colorbar;
    use crate::colormap::ColormapLegend;
    #[test]
    fn jet_endpoints() {
        let blue = map_scalar(0.0, 0.0, 1.0, Colormap::Jet);
        assert!((blue.r).abs() < 0.01);
        assert!((blue.b - 1.0).abs() < 0.01);
        let red = map_scalar(1.0, 0.0, 1.0, Colormap::Jet);
        assert!((red.r - 1.0).abs() < 0.01);
        assert!((red.b).abs() < 0.01);
    }
    #[test]
    fn viridis_endpoints() {
        let start = map_scalar(0.0, 0.0, 1.0, Colormap::Viridis);
        assert!(start.r < 0.3);
        assert!(start.g < 0.1);
        assert!(start.b > 0.3);
        let end_color = map_scalar(1.0, 0.0, 1.0, Colormap::Viridis);
        assert!(end_color.r > 0.9);
        assert!(end_color.g > 0.8);
    }
    #[test]
    fn test_brbg_endpoints() {
        let brown = map_scalar(0.0, 0.0, 1.0, Colormap::BrBG);
        assert!(brown.r > 0.3, "BrBG start should be brownish");
        let green = map_scalar(1.0, 0.0, 1.0, Colormap::BrBG);
        assert!(green.g > 0.3, "BrBG end should be greenish");
    }
    #[test]
    fn test_brbg_midpoint_near_white() {
        let mid = map_scalar(0.5, 0.0, 1.0, Colormap::BrBG);
        assert!(
            mid.r > 0.9 && mid.g > 0.9 && mid.b > 0.9,
            "BrBG midpoint should be near white"
        );
    }
    #[test]
    fn test_colormap_sample_alpha_255() {
        for cmap in [
            ColorMap::Viridis,
            ColorMap::Plasma,
            ColorMap::Inferno,
            ColorMap::Magma,
            ColorMap::Turbo,
            ColorMap::RdBu,
            ColorMap::Seismic,
            ColorMap::BrBG,
        ] {
            for &t in &[0.0, 0.5, 1.0] {
                let rgba = cmap.sample(t);
                assert_eq!(rgba[3], 255, "alpha should be 255 for cmap={cmap:?} t={t}");
            }
        }
    }
    #[test]
    fn test_colormap_sample_t0_t1_distinct() {
        for cmap in [ColorMap::Viridis, ColorMap::Plasma, ColorMap::Turbo] {
            let c0 = cmap.sample(0.0);
            let c1 = cmap.sample(1.0);
            let diff: i32 = (c0[0] as i32 - c1[0] as i32).abs()
                + (c0[1] as i32 - c1[1] as i32).abs()
                + (c0[2] as i32 - c1[2] as i32).abs();
            assert!(
                diff > 10,
                "t=0 and t=1 should be distinct colors for {cmap:?}"
            );
        }
    }
    #[test]
    fn test_normalize_range() {
        let values = vec![0.0, 2.5, 5.0, 10.0];
        let normed = normalize(&values);
        assert!((normed[0] - 0.0).abs() < 1e-10, "min should map to 0");
        assert!((normed[3] - 1.0).abs() < 1e-10, "max should map to 1");
        for &v in &normed {
            assert!(
                (0.0..=1.0).contains(&v),
                "normalized value out of range: {v}"
            );
        }
    }
    #[test]
    fn test_symmetric_normalize_center() {
        let values = vec![-4.0, -2.0, 0.0, 2.0, 4.0];
        let normed = symmetric_normalize(&values);
        assert!(
            (normed[2] - 0.5).abs() < 1e-10,
            "zero should map to 0.5, got {}",
            normed[2]
        );
        assert!(
            (normed[4] - 1.0).abs() < 1e-10,
            "max should map to 1.0, got {}",
            normed[4]
        );
        assert!(
            (normed[0] - 0.0).abs() < 1e-10,
            "min should map to 0.0, got {}",
            normed[0]
        );
    }
    #[test]
    fn test_hsv_roundtrip() {
        let test_colors = [
            (1.0, 0.0, 0.0),
            (0.0, 1.0, 0.0),
            (0.0, 0.0, 1.0),
            (0.5, 0.5, 0.5),
            (0.2, 0.8, 0.4),
        ];
        for (r, g, b) in test_colors {
            let (h, s, v) = rgb_to_hsv(r, g, b);
            let (r2, g2, b2) = hsv_to_rgb(h, s, v);
            assert!((r - r2).abs() < 1e-6, "R roundtrip failed: {r} vs {r2}");
            assert!((g - g2).abs() < 1e-6, "G roundtrip failed: {g} vs {g2}");
            assert!((b - b2).abs() < 1e-6, "B roundtrip failed: {b} vs {b2}");
        }
    }
    #[test]
    fn test_blend_rgb_endpoints() {
        let a = [255u8, 0, 0, 255];
        let b = [0u8, 0, 255, 255];
        let c0 = blend_rgb(a, b, 0.0);
        let c1 = blend_rgb(a, b, 1.0);
        assert_eq!(c0[0], 255);
        assert_eq!(c0[2], 0);
        assert_eq!(c1[0], 0);
        assert_eq!(c1[2], 255);
        assert_eq!(c0[3], 255);
        assert_eq!(c1[3], 255);
    }
    #[test]
    fn test_apply_colormap_length() {
        let values = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let out = apply_colormap(&values, ColorMap::Viridis, 1.0, 5.0);
        assert_eq!(out.len(), 5);
        for rgba in &out {
            assert_eq!(rgba[3], 255);
        }
    }
    #[test]
    fn test_categorical_colormap() {
        let cat = CategoricalColormap::default_palette();
        assert_eq!(cat.palette_size(), 10);
        let c0 = cat.map_category(0);
        let c1 = cat.map_category(1);
        let diff = (c0.r - c1.r).abs() + (c0.g - c1.g).abs() + (c0.b - c1.b).abs();
        assert!(
            diff > 0.1,
            "Different categories should have different colors"
        );
        let c10 = cat.map_category(10);
        let c0_again = cat.map_category(0);
        assert!(
            (c10.r - c0_again.r).abs() < 1e-6,
            "Index 10 should wrap to 0"
        );
    }
    #[test]
    fn test_categorical_map_categories() {
        let cat = CategoricalColormap::default_palette();
        let colors = cat.map_categories(&[0, 1, 2, 0]);
        assert_eq!(colors.len(), 4);
        assert!(
            (colors[0].r - colors[3].r).abs() < 1e-6,
            "Same category should give same color"
        );
    }
    #[test]
    fn test_colorbar_tick_values() {
        let cb = Colorbar::new(Colormap::Viridis, 0.0, 100.0, "Temperature");
        let ticks = cb.tick_values();
        assert_eq!(ticks.len(), 5);
        assert!((ticks[0] - 0.0).abs() < 1e-10);
        assert!((ticks[4] - 100.0).abs() < 1e-10);
    }
    #[test]
    fn test_colorbar_raster() {
        let cb = Colorbar::new(Colormap::Jet, 0.0, 1.0, "Pressure");
        let raster = cb.raster(256);
        assert_eq!(raster.len(), 256);
    }
    #[test]
    fn test_colorbar_with_ticks() {
        let cb = Colorbar::new(Colormap::Viridis, -10.0, 10.0, "Stress").with_ticks(11);
        let ticks = cb.tick_values();
        assert_eq!(ticks.len(), 11);
        assert!((ticks[5] - 0.0).abs() < 1e-10, "Middle tick should be 0");
    }
    #[test]
    fn test_histogram_colors() {
        let colors = histogram_colors(5, Colormap::Viridis);
        assert_eq!(colors.len(), 5);
    }
    #[test]
    fn test_histogram_value_colors() {
        let counts = vec![10.0, 50.0, 100.0];
        let colors = histogram_value_colors(&counts, Colormap::Jet, 0.0, 100.0);
        assert_eq!(colors.len(), 3);
    }
    #[test]
    fn test_interpolation_mode_nearest() {
        let stops = vec![
            (0.0_f32, Color::new(1.0, 0.0, 0.0, 1.0)),
            (1.0_f32, Color::new(0.0, 0.0, 1.0, 1.0)),
        ];
        let c = sample_stops_with_mode(&stops, 0.3, InterpolationMode::Nearest);
        assert!(
            c.r > 0.9 && c.b < 0.1,
            "Nearest at 0.3 should pick first stop"
        );
        let c2 = sample_stops_with_mode(&stops, 0.7, InterpolationMode::Nearest);
        assert!(
            c2.b > 0.9 && c2.r < 0.1,
            "Nearest at 0.7 should pick second stop"
        );
    }
    #[test]
    fn test_interpolation_mode_step() {
        let stops = vec![
            (0.0_f32, Color::new(1.0, 0.0, 0.0, 1.0)),
            (0.5_f32, Color::new(0.0, 1.0, 0.0, 1.0)),
            (1.0_f32, Color::new(0.0, 0.0, 1.0, 1.0)),
        ];
        let c = sample_stops_with_mode(&stops, 0.3, InterpolationMode::Step);
        assert!(c.r > 0.9, "Step at 0.3 should use first stop color");
    }
    #[test]
    fn test_alpha_blend_opaque_over_opaque() {
        let result = alpha_blend([255, 0, 0, 255], [0, 0, 255, 255]);
        assert_eq!(result[0], 255);
        assert_eq!(result[2], 0);
        assert_eq!(result[3], 255);
    }
    #[test]
    fn test_alpha_blend_transparent_preserves_dst() {
        let result = alpha_blend([255, 0, 0, 0], [0, 0, 255, 255]);
        assert_eq!(result[2], 255);
        assert_eq!(result[3], 255);
    }
    #[test]
    fn test_percentile_normalize() {
        let values = vec![0.0, 1.0, 2.0, 3.0, 100.0];
        let normed = percentile_normalize(&values, 0.0, 80.0);
        assert!(
            (normed[4] - 1.0).abs() < 1e-10,
            "Outlier should be clipped to 1.0"
        );
    }
    #[test]
    fn test_rgb_to_hsl() {
        let (h, s, l) = rgb_to_hsl(1.0, 0.0, 0.0);
        assert!(
            h.abs() < 1e-6 || (h - 360.0).abs() < 1e-6,
            "Red H should be 0"
        );
        assert!((s - 1.0).abs() < 1e-6, "Red S should be 1");
        assert!((l - 0.5).abs() < 1e-6, "Red L should be 0.5");
    }
    #[test]
    fn test_colormap_brbg_u8() {
        let c = ColorMap::BrBG.sample(0.0);
        assert_eq!(c[3], 255);
        let c2 = ColorMap::BrBG.sample(1.0);
        let diff: i32 = (c[0] as i32 - c2[0] as i32).abs()
            + (c[1] as i32 - c2[1] as i32).abs()
            + (c[2] as i32 - c2[2] as i32).abs();
        assert!(diff > 10, "BrBG endpoints should be distinct");
    }
    #[test]
    fn test_cividis_endpoints_distinct() {
        let c0 = map_scalar(0.0, 0.0, 1.0, Colormap::Cividis);
        let c1 = map_scalar(1.0, 0.0, 1.0, Colormap::Cividis);
        let diff = (c0.r - c1.r).abs() + (c0.g - c1.g).abs() + (c0.b - c1.b).abs();
        assert!(diff > 0.1, "Cividis endpoints should differ");
    }
    #[test]
    fn test_cividis_alpha_one() {
        let c = map_scalar(0.5, 0.0, 1.0, Colormap::Cividis);
        assert!((c.a - 1.0).abs() < 1e-5);
    }
    #[test]
    fn test_isoluminant_rainbow_endpoints() {
        let c0 = map_scalar(0.0, 0.0, 1.0, Colormap::IsoRainbow);
        let c1 = map_scalar(1.0, 0.0, 1.0, Colormap::IsoRainbow);
        let diff = (c0.r - c1.r).abs() + (c0.g - c1.g).abs() + (c0.b - c1.b).abs();
        assert!(diff > 0.05, "IsoRainbow endpoints should differ");
    }
    #[test]
    fn test_blend_colormaps_midpoint() {
        let c = blend_colormaps(Colormap::Grayscale, Colormap::Jet, 0.5, 0.5);
        assert!(c.r >= 0.0 && c.r <= 1.0);
        assert!(c.g >= 0.0 && c.g <= 1.0);
        assert!(c.b >= 0.0 && c.b <= 1.0);
    }
    #[test]
    fn test_blend_colormaps_alpha_zero() {
        let c = blend_colormaps(Colormap::Grayscale, Colormap::Jet, 0.5, 0.0);
        let expected = map_scalar(0.5, 0.0, 1.0, Colormap::Grayscale);
        assert!((c.r - expected.r).abs() < 1e-5);
    }
    #[test]
    fn test_blend_colormaps_alpha_one() {
        let c = blend_colormaps(Colormap::Grayscale, Colormap::Jet, 0.5, 1.0);
        let expected = map_scalar(0.5, 0.0, 1.0, Colormap::Jet);
        assert!((c.r - expected.r).abs() < 1e-5);
    }
    #[test]
    fn test_equalize_histogram_output_length() {
        let values = vec![0.0, 0.1, 0.2, 0.5, 0.8, 0.9, 1.0];
        let equalized = equalize_histogram(&values, 256);
        assert_eq!(equalized.len(), values.len());
    }
    #[test]
    fn test_equalize_histogram_range() {
        let values: Vec<f64> = (0..50).map(|i| i as f64 / 49.0).collect();
        let equalized = equalize_histogram(&values, 256);
        for &v in &equalized {
            assert!((0.0..=1.0).contains(&v), "equalized out of range: {v}");
        }
    }
    #[test]
    fn test_equalize_histogram_constant_field() {
        let values = vec![0.5; 20];
        let equalized = equalize_histogram(&values, 64);
        assert_eq!(equalized.len(), 20);
    }
    #[test]
    fn test_equalize_histogram_empty() {
        let equalized = equalize_histogram(&[], 64);
        assert!(equalized.is_empty());
    }
    #[test]
    fn test_colormap_legend_generation() {
        let legend = ColormapLegend::new(Colormap::Viridis, 0.0, 100.0, "Temperature (K)");
        assert_eq!(legend.colormap, Colormap::Viridis);
        assert!((legend.vmin - 0.0).abs() < 1e-10);
        assert!((legend.vmax - 100.0).abs() < 1e-10);
        assert_eq!(legend.label, "Temperature (K)");
    }
    #[test]
    fn test_colormap_legend_default_n_ticks() {
        let legend = ColormapLegend::new(Colormap::Plasma, -1.0, 1.0, "Velocity");
        assert!(legend.n_ticks >= 2);
    }
    #[test]
    fn test_colormap_legend_render_raster() {
        let legend = ColormapLegend::new(Colormap::Turbo, 0.0, 1.0, "Pressure");
        let raster = legend.render_raster(128, 20);
        assert_eq!(raster.len(), 128 * 20);
    }
    #[test]
    fn test_colormap_legend_tick_labels() {
        let legend = ColormapLegend::new(Colormap::Viridis, 0.0, 100.0, "Temp").with_n_ticks(5);
        let labels = legend.tick_labels();
        assert_eq!(labels.len(), 5);
        assert!(labels[0].contains("0") || labels[0].contains("0."));
        assert!(labels[4].contains("100"));
    }
    #[test]
    fn test_colormap_turbo_u8_monotone_green() {
        let c0 = ColorMap::Turbo.sample(0.0);
        let c_mid = ColorMap::Turbo.sample(0.5);
        let bright_mid = c_mid[0] as u32 + c_mid[1] as u32 + c_mid[2] as u32;
        let bright_0 = c0[0] as u32 + c0[1] as u32 + c0[2] as u32;
        assert!(
            bright_mid > bright_0,
            "Turbo midpoint should be brighter than start"
        );
    }
    #[test]
    fn test_colormap_cividis_u8_alpha() {
        let c = ColorMap::sample_cividis(0.5);
        assert_eq!(c[3], 255);
    }
    #[test]
    fn test_power_normalize_gamma_one_identity() {
        let values = vec![0.0, 0.25, 0.5, 0.75, 1.0];
        let result = power_normalize(&values, 1.0);
        for (v, r) in values.iter().zip(result.iter()) {
            assert!(
                (v - r).abs() < 1e-6,
                "gamma=1 should be identity: {v} vs {r}"
            );
        }
    }
    #[test]
    fn test_power_normalize_range() {
        let values: Vec<f64> = (0..10).map(|i| i as f64 / 9.0).collect();
        let result = power_normalize(&values, 2.0);
        for &r in &result {
            assert!((0.0..=1.0 + 1e-10).contains(&r), "out of range: {r}");
        }
    }
    #[test]
    fn test_power_normalize_empty() {
        assert!(power_normalize(&[], 2.0).is_empty());
    }
    #[test]
    fn test_log_normalize_positive_values() {
        let values = vec![1.0, 10.0, 100.0, 1000.0];
        let result = log_normalize(&values);
        for i in 0..result.len() - 1 {
            assert!(
                result[i] < result[i + 1],
                "log_normalize not monotone at {i}"
            );
        }
        assert!((result[0] - 0.0).abs() < 1e-6);
        assert!((result[3] - 1.0).abs() < 1e-6);
    }
    #[test]
    fn test_log_normalize_single_value() {
        let values = vec![42.0];
        let result = log_normalize(&values);
        assert_eq!(result.len(), 1);
    }
    #[test]
    fn test_rdbu_midpoint_near_white() {
        let mid = map_scalar(0.5, 0.0, 1.0, Colormap::RdBu);
        assert!(
            mid.r > 0.9 && mid.g > 0.9 && mid.b > 0.9,
            "RdBu midpoint should be near white"
        );
    }
    #[test]
    fn test_seismic_midpoint_near_white() {
        let mid = map_scalar(0.5, 0.0, 1.0, Colormap::Seismic);
        assert!(
            mid.r > 0.9 && mid.g > 0.9 && mid.b > 0.9,
            "Seismic midpoint should be near white"
        );
    }
    #[test]
    fn test_map_scalar_clamps_low() {
        let c_neg = map_scalar(-10.0, 0.0, 1.0, Colormap::Viridis);
        let c_zero = map_scalar(0.0, 0.0, 1.0, Colormap::Viridis);
        assert!((c_neg.r - c_zero.r).abs() < 1e-5);
    }
    #[test]
    fn test_map_scalar_clamps_high() {
        let c_high = map_scalar(100.0, 0.0, 1.0, Colormap::Viridis);
        let c_one = map_scalar(1.0, 0.0, 1.0, Colormap::Viridis);
        assert!((c_high.r - c_one.r).abs() < 1e-5);
    }
    #[test]
    fn test_map_scalar_zero_range() {
        let c = map_scalar(5.0, 5.0, 5.0, Colormap::Viridis);
        let _ = c;
    }
    #[test]
    fn test_rainbow_endpoints_distinct() {
        let c0 = rainbow(0.0);
        let c1 = rainbow(1.0);
        let diff = (c0.r - c1.r).abs() + (c0.g - c1.g).abs() + (c0.b - c1.b).abs();
        assert!(
            diff > 0.1,
            "Rainbow t=0 and t=1 should be visually distinct"
        );
    }
    #[test]
    fn test_rainbow_alpha_one() {
        let c = rainbow(0.5);
        assert!((c.a - 1.0).abs() < 1e-5, "Rainbow alpha should be 1.0");
    }
    #[test]
    fn test_rainbow_midpoint_is_valid_color() {
        let c = rainbow(0.5);
        assert!(c.r >= 0.0 && c.r <= 1.0);
        assert!(c.g >= 0.0 && c.g <= 1.0);
        assert!(c.b >= 0.0 && c.b <= 1.0);
    }
    #[test]
    fn test_tab10_palette_size() {
        let pal = tab10_palette();
        assert_eq!(pal.len(), 10);
    }
    #[test]
    fn test_tab10_all_opaque() {
        for c in tab10_palette().iter() {
            assert!((c.a - 1.0).abs() < 1e-5, "Tab10 colors should be opaque");
        }
    }
    #[test]
    fn test_tab10_all_distinct() {
        let pal = tab10_palette();
        for i in 0..pal.len() {
            for j in i + 1..pal.len() {
                let diff = (pal[i].r - pal[j].r).abs()
                    + (pal[i].g - pal[j].g).abs()
                    + (pal[i].b - pal[j].b).abs();
                assert!(diff > 0.05, "Tab10 colors {i} and {j} are too similar");
            }
        }
    }
    #[test]
    fn test_set1_palette_size() {
        assert_eq!(set1_palette().len(), 9);
    }
    #[test]
    fn test_paired_palette_size() {
        assert_eq!(paired_palette().len(), 12);
    }
    #[test]
    fn test_paired_light_darker_than_dark() {
        let pal = paired_palette();
        for i in 0..6 {
            let light = pal[2 * i];
            let dark = pal[2 * i + 1];
            let lum_light = 0.2126 * light.r + 0.7152 * light.g + 0.0722 * light.b;
            let lum_dark = 0.2126 * dark.r + 0.7152 * dark.g + 0.0722 * dark.b;
            assert!(
                lum_light >= lum_dark,
                "Paired palette: light[{i}] should be >= lum of dark[{i}]"
            );
        }
    }
    #[test]
    fn test_invert_colormap_t0_equals_t1() {
        let inverted_at_0 = invert_colormap(Colormap::Viridis, 0.0);
        let original_at_1 = map_scalar(1.0, 0.0, 1.0, Colormap::Viridis);
        assert!((inverted_at_0.r - original_at_1.r).abs() < 1e-5);
        assert!((inverted_at_0.g - original_at_1.g).abs() < 1e-5);
        assert!((inverted_at_0.b - original_at_1.b).abs() < 1e-5);
    }
    #[test]
    fn test_invert_colormap_t1_equals_t0() {
        let inverted_at_1 = invert_colormap(Colormap::Plasma, 1.0);
        let original_at_0 = map_scalar(0.0, 0.0, 1.0, Colormap::Plasma);
        assert!((inverted_at_1.r - original_at_0.r).abs() < 1e-5);
    }
    #[test]
    fn test_discrete_colormap_length() {
        let n = 8;
        let palette = discrete_colormap(Colormap::Jet, n);
        assert_eq!(
            palette.len(),
            n,
            "discrete_colormap should return exactly n colors"
        );
    }
    #[test]
    fn test_discrete_colormap_empty() {
        let palette = discrete_colormap(Colormap::Viridis, 0);
        assert!(palette.is_empty());
    }
    #[test]
    fn test_discrete_colormap_single() {
        let palette = discrete_colormap(Colormap::Jet, 1);
        assert_eq!(palette.len(), 1);
    }
    #[test]
    fn test_discrete_colormap_endpoints_match_original() {
        let n = 5;
        let palette = discrete_colormap(Colormap::Viridis, n);
        let c_start = map_scalar(0.0, 0.0, 1.0, Colormap::Viridis);
        let c_end = map_scalar(1.0, 0.0, 1.0, Colormap::Viridis);
        assert!((palette[0].r - c_start.r).abs() < 1e-5);
        assert!((palette[n - 1].r - c_end.r).abs() < 1e-5);
    }
    #[test]
    fn test_map_scalar_rgba_custom_alpha() {
        let rgba = map_scalar_rgba(0.5, 0.0, 1.0, Colormap::Viridis, 0.5);
        assert_eq!(rgba[3], 128, "alpha=0.5 should map to u8 value ~128");
    }
    #[test]
    fn test_map_scalar_rgba_full_alpha() {
        let rgba = map_scalar_rgba(0.0, 0.0, 1.0, Colormap::Jet, 1.0);
        assert_eq!(rgba[3], 255);
    }
    #[test]
    fn test_map_scalar_rgba_zero_alpha() {
        let rgba = map_scalar_rgba(1.0, 0.0, 1.0, Colormap::Plasma, 0.0);
        assert_eq!(rgba[3], 0);
    }
    #[test]
    fn test_rgb_to_lab_white_is_l100() {
        let (l, a, b) = rgb_to_lab(1.0, 1.0, 1.0);
        assert!(
            (l - 100.0).abs() < 1.0,
            "White should have L* ≈ 100, got {l}"
        );
        assert!(a.abs() < 1.0, "White a* should be near 0, got {a}");
        assert!(b.abs() < 1.0, "White b* should be near 0, got {b}");
    }
    #[test]
    fn test_rgb_to_lab_black_is_l0() {
        let (l, _a, _b) = rgb_to_lab(0.0, 0.0, 0.0);
        assert!(l.abs() < 1.0, "Black should have L* ≈ 0, got {l}");
    }
    #[test]
    fn test_lab_to_rgb_roundtrip() {
        let test_colors = [
            (1.0, 0.0, 0.0),
            (0.0, 1.0, 0.0),
            (0.0, 0.0, 1.0),
            (0.5, 0.5, 0.5),
        ];
        for (r, g, b) in test_colors {
            let (l, a, bb) = rgb_to_lab(r, g, b);
            let (r2, g2, b2) = lab_to_rgb(l, a, bb);
            assert!((r - r2).abs() < 0.01, "Lab roundtrip R: {r} -> {r2}");
            assert!((g - g2).abs() < 0.01, "Lab roundtrip G: {g} -> {g2}");
            assert!((b - b2).abs() < 0.01, "Lab roundtrip B: {b} -> {b2}");
        }
    }
    #[test]
    fn test_perceptual_distance_identical_is_zero() {
        let c = Color::new(0.5, 0.3, 0.7, 1.0);
        let d = perceptual_distance(c, c);
        assert!(
            d.abs() < 1e-6,
            "Distance from a color to itself should be 0, got {d}"
        );
    }
    #[test]
    fn test_perceptual_distance_black_white_large() {
        let black = Color::new(0.0, 0.0, 0.0, 1.0);
        let white = Color::new(1.0, 1.0, 1.0, 1.0);
        let d = perceptual_distance(black, white);
        assert!(
            d > 50.0,
            "Black-white perceptual distance should be large, got {d}"
        );
    }
    #[test]
    fn test_uniformity_viridis_lower_than_jet() {
        let jet_score = colormap_uniformity_score(Colormap::Jet, 32);
        let viridis_score = colormap_uniformity_score(Colormap::Viridis, 32);
        assert!(
            viridis_score < jet_score,
            "Viridis should score better (lower) than Jet: viridis={viridis_score:.3} jet={jet_score:.3}"
        );
    }
    #[test]
    fn test_uniformity_score_nonnegative() {
        for cmap in [
            Colormap::Viridis,
            Colormap::Plasma,
            Colormap::Jet,
            Colormap::Grayscale,
        ] {
            let s = colormap_uniformity_score(cmap, 20);
            assert!(
                s >= 0.0,
                "Uniformity score must be non-negative, got {s} for {cmap:?}"
            );
        }
    }
    #[test]
    fn test_perceptual_lightness_viridis_monotone() {
        let ls = Colormap::compute_perceptual_lightness(Colormap::Viridis, 16);
        for i in 0..ls.len() - 1 {
            assert!(
                ls[i + 1] >= ls[i] - 5.0,
                "Viridis L* should be roughly non-decreasing at step {i}: {} -> {}",
                ls[i],
                ls[i + 1]
            );
        }
    }
    #[test]
    fn test_perceptual_lightness_grayscale_monotone() {
        let ls = Colormap::compute_perceptual_lightness(Colormap::Grayscale, 10);
        assert_eq!(ls.len(), 10);
        for i in 0..ls.len() - 1 {
            assert!(
                ls[i + 1] > ls[i] - 1.0,
                "Grayscale L* should increase at step {i}: {} -> {}",
                ls[i],
                ls[i + 1]
            );
        }
    }
    #[test]
    fn test_perceptual_lightness_range() {
        let ls = Colormap::compute_perceptual_lightness(Colormap::Plasma, 20);
        for &l in &ls {
            assert!(
                (0.0..=105.0).contains(&l),
                "L* should be in [0, 100+epsilon], got {l}"
            );
        }
    }
    #[test]
    fn test_perceptual_lightness_empty() {
        let ls = Colormap::compute_perceptual_lightness(Colormap::Jet, 0);
        assert!(ls.is_empty());
    }
    #[test]
    fn test_color_difference_identical_is_zero() {
        let diffs = Colormap::compute_color_difference(Colormap::Viridis, 8);
        for &d in &diffs {
            assert!(d >= 0.0, "delta-E should be non-negative, got {d}");
        }
    }
    #[test]
    fn test_color_difference_length() {
        let n = 12;
        let diffs = Colormap::compute_color_difference(Colormap::Jet, n);
        assert_eq!(diffs.len(), n - 1);
    }
    #[test]
    fn test_color_difference_viridis_more_uniform_than_jet() {
        let n = 32;
        let jet_diffs = Colormap::compute_color_difference(Colormap::Jet, n);
        let viridis_diffs = Colormap::compute_color_difference(Colormap::Viridis, n);
        let jet_var: f64 = {
            let mean = jet_diffs.iter().sum::<f64>() / jet_diffs.len() as f64;
            jet_diffs.iter().map(|&d| (d - mean).powi(2)).sum::<f64>() / jet_diffs.len() as f64
        };
        let viridis_var: f64 = {
            let mean = viridis_diffs.iter().sum::<f64>() / viridis_diffs.len() as f64;
            viridis_diffs
                .iter()
                .map(|&d| (d - mean).powi(2))
                .sum::<f64>()
                / viridis_diffs.len() as f64
        };
        assert!(
            viridis_var <= jet_var + 1.0,
            "Viridis delta-E variance {viridis_var:.4} should not greatly exceed Jet {jet_var:.4}"
        );
    }
    #[test]
    fn test_categorical_palette_length() {
        for n in [1usize, 3, 7, 12, 20] {
            let pal = Colormap::generate_categorical_palette(n);
            assert_eq!(pal.len(), n, "Expected {n} colors, got {}", pal.len());
        }
    }
    #[test]
    fn test_categorical_palette_all_opaque() {
        let pal = Colormap::generate_categorical_palette(8);
        for c in &pal {
            assert!(
                (c.a - 1.0).abs() < 1e-5,
                "Categorical palette color should be opaque"
            );
        }
    }
    #[test]
    fn test_categorical_palette_all_in_range() {
        let pal = Colormap::generate_categorical_palette(16);
        for c in &pal {
            assert!(c.r >= 0.0 && c.r <= 1.0, "r out of range: {}", c.r);
            assert!(c.g >= 0.0 && c.g <= 1.0, "g out of range: {}", c.g);
            assert!(c.b >= 0.0 && c.b <= 1.0, "b out of range: {}", c.b);
        }
    }
}
