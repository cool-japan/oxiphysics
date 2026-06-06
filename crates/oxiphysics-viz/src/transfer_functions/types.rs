//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

/// Opacity profile type.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum OpacityProfileKind {
    /// Linearly ramps from 0 to 1.
    LinearRamp,
    /// Constant opacity.
    Constant(f64),
    /// Gaussian bump centered at `center` with standard deviation `sigma`.
    Gaussian {
        /// Center position in \[0, 1\].
        center: f64,
        /// Standard deviation of the Gaussian.
        sigma: f64,
    },
    /// Step function: 0 below `threshold`, 1 at or above.
    Step {
        /// Scalar position at which opacity jumps from 0 to 1.
        threshold: f64,
    },
    /// Tent / triangle profile peaked at `center`.
    Tent {
        /// Position of the tent peak in \[0, 1\].
        center: f64,
        /// Half-width of the tent base.
        half_width: f64,
    },
}
/// A piecewise-linear RGBA transfer function.
///
/// The transfer function is defined by a sorted list of [`ControlPoint`]s.
/// Lookup uses linear interpolation between adjacent control points.
///
/// # Example
///
/// ```no_run
/// use oxiphysics_viz::transfer_functions::{TransferFunction, Rgba};
///
/// let mut tf = TransferFunction::new();
/// tf.add_point(0.0, Rgba::new(0.0, 0.0, 1.0, 0.0)); // transparent blue at cold
/// tf.add_point(1.0, Rgba::new(1.0, 0.0, 0.0, 1.0)); // opaque red at hot
/// let color = tf.sample(0.5);
/// assert!((color.r - 0.5).abs() < 1e-6);
/// ```
#[derive(Debug, Clone)]
pub struct TransferFunction {
    /// Sorted list of control points.
    pub(super) points: Vec<ControlPoint>,
}
impl TransferFunction {
    /// Create an empty transfer function.
    pub fn new() -> Self {
        Self { points: Vec::new() }
    }
    /// Add or replace a control point at `value`.
    ///
    /// The internal list is kept sorted.
    pub fn add_point(&mut self, value: f64, color: Rgba) {
        let cp = ControlPoint::new(value, color);
        let pos = self.points.partition_point(|p| p.value < cp.value);
        if pos < self.points.len() && (self.points[pos].value - cp.value).abs() < 1e-12 {
            self.points[pos] = cp;
        } else {
            self.points.insert(pos, cp);
        }
    }
    /// Remove the control point nearest to `value` (within `tol`).
    ///
    /// Returns `true` if a point was removed.
    pub fn remove_near(&mut self, value: f64, tol: f64) -> bool {
        if let Some(idx) = self
            .points
            .iter()
            .position(|p| (p.value - value).abs() <= tol)
        {
            self.points.remove(idx);
            true
        } else {
            false
        }
    }
    /// Number of control points.
    pub fn point_count(&self) -> usize {
        self.points.len()
    }
    /// Read-only access to the control points.
    pub fn points(&self) -> &[ControlPoint] {
        &self.points
    }
    /// Sample the transfer function at `t` (clamped to \[0, 1\]).
    ///
    /// Returns transparent black if no control points are defined.
    pub fn sample(&self, t: f64) -> Rgba {
        let t = t.clamp(0.0, 1.0);
        if self.points.is_empty() {
            return Rgba::transparent();
        }
        if self.points.len() == 1 {
            return self.points[0].color;
        }
        if t <= self.points[0].value {
            return self.points[0].color;
        }
        if t >= self.points[self.points.len() - 1].value {
            return self.points[self.points.len() - 1].color;
        }
        let hi = self.points.partition_point(|p| p.value <= t);
        let lo = hi.saturating_sub(1);
        let lo = lo.min(self.points.len() - 2);
        let hi = lo + 1;
        let a = &self.points[lo];
        let b = &self.points[hi];
        let span = b.value - a.value;
        let local_t = if span < 1e-15 {
            0.5
        } else {
            (t - a.value) / span
        };
        Rgba::lerp(a.color, b.color, local_t)
    }
    /// Build a lookup table of `n` evenly-spaced RGBA samples.
    pub fn build_lut(&self, n: usize) -> Vec<Rgba> {
        (0..n)
            .map(|i| {
                let t = i as f64 / (n.saturating_sub(1).max(1)) as f64;
                self.sample(t)
            })
            .collect()
    }
    /// Hot gas / temperature: transparent dark-blue → opaque yellow-white.
    pub fn temperature() -> Self {
        let mut tf = Self::new();
        tf.add_point(0.0, Rgba::new(0.0, 0.0, 0.2, 0.0));
        tf.add_point(0.2, Rgba::new(0.0, 0.0, 0.8, 0.05));
        tf.add_point(0.45, Rgba::new(0.8, 0.0, 0.0, 0.3));
        tf.add_point(0.7, Rgba::new(1.0, 0.5, 0.0, 0.6));
        tf.add_point(0.9, Rgba::new(1.0, 1.0, 0.0, 0.85));
        tf.add_point(1.0, Rgba::new(1.0, 1.0, 1.0, 1.0));
        tf
    }
    /// Pressure: blue (low pressure) → white (ambient) → red (high pressure).
    pub fn pressure() -> Self {
        let mut tf = Self::new();
        tf.add_point(0.0, Rgba::new(0.0, 0.0, 1.0, 0.8));
        tf.add_point(0.5, Rgba::new(1.0, 1.0, 1.0, 0.05));
        tf.add_point(1.0, Rgba::new(1.0, 0.0, 0.0, 0.8));
        tf
    }
    /// Velocity magnitude: transparent at zero, opaque green/yellow at high speed.
    pub fn velocity_magnitude() -> Self {
        let mut tf = Self::new();
        tf.add_point(0.0, Rgba::new(0.0, 0.0, 1.0, 0.0));
        tf.add_point(0.3, Rgba::new(0.0, 1.0, 0.5, 0.3));
        tf.add_point(0.7, Rgba::new(1.0, 1.0, 0.0, 0.7));
        tf.add_point(1.0, Rgba::new(1.0, 0.2, 0.0, 0.9));
        tf
    }
    /// Density / mass: transparent at low density, opaque smoke-gray at high.
    pub fn density() -> Self {
        let mut tf = Self::new();
        tf.add_point(0.0, Rgba::new(0.2, 0.2, 0.2, 0.0));
        tf.add_point(0.5, Rgba::new(0.5, 0.5, 0.5, 0.4));
        tf.add_point(1.0, Rgba::new(0.8, 0.8, 0.8, 0.9));
        tf
    }
}
/// A single control point in a piecewise-linear transfer function.
#[derive(Debug, Clone, Copy)]
pub struct ControlPoint {
    /// Normalized scalar value in \[0, 1\].
    pub value: f64,
    /// RGBA color at this control point.
    pub color: Rgba,
}
impl ControlPoint {
    /// Construct a control point at `value` with the given RGBA color.
    pub fn new(value: f64, color: Rgba) -> Self {
        Self {
            value: value.clamp(0.0, 1.0),
            color,
        }
    }
}
/// Represents the state of an interactive transfer function editor.
///
/// Tracks which control point is selected, the current scalar range, and
/// provides convenience methods for mouse-based editing workflows.
#[derive(Debug, Clone)]
pub struct TransferFunctionEditor {
    /// The transfer function being edited.
    pub tf: TransferFunction,
    /// Scalar minimum for domain mapping.
    pub scalar_min: f64,
    /// Scalar maximum for domain mapping.
    pub scalar_max: f64,
    /// Index of the currently selected control point, if any.
    pub selected: Option<usize>,
    /// Whether the editor is in opacity-only mode (color kept fixed).
    pub opacity_only: bool,
}
impl TransferFunctionEditor {
    /// Create an editor for the given transfer function.
    pub fn new(tf: TransferFunction, scalar_min: f64, scalar_max: f64) -> Self {
        Self {
            tf,
            scalar_min,
            scalar_max,
            selected: None,
            opacity_only: false,
        }
    }
    /// Map a raw scalar value to \[0, 1\] for TF lookup.
    pub fn normalize(&self, scalar: f64) -> f64 {
        let range = self.scalar_max - self.scalar_min;
        if range.abs() < 1e-30 {
            return 0.5;
        }
        ((scalar - self.scalar_min) / range).clamp(0.0, 1.0)
    }
    /// Sample the transfer function at a raw scalar value.
    pub fn sample_scalar(&self, scalar: f64) -> Rgba {
        self.tf.sample(self.normalize(scalar))
    }
    /// Select the control point nearest to the given normalized position `t`.
    ///
    /// Selects the closest point within `tolerance`.
    pub fn select_nearest(&mut self, t: f64, tolerance: f64) {
        let best = self
            .tf
            .points()
            .iter()
            .enumerate()
            .min_by(|(_, a), (_, b)| {
                let da = (a.value - t).abs();
                let db = (b.value - t).abs();
                da.partial_cmp(&db).unwrap_or(std::cmp::Ordering::Equal)
            });
        if let Some((idx, cp)) = best {
            if (cp.value - t).abs() <= tolerance {
                self.selected = Some(idx);
            } else {
                self.selected = None;
            }
        }
    }
    /// Move the selected control point to position `new_t`.
    ///
    /// Does nothing if no point is selected.
    pub fn move_selected(&mut self, new_t: f64) {
        if let Some(idx) = self.selected {
            let color = self.tf.points()[idx].color;
            let old_val = self.tf.points()[idx].value;
            let _ = self.tf.remove_near(old_val, 1e-10);
            self.tf.add_point(new_t, color);
            let new_idx = self
                .tf
                .points()
                .iter()
                .position(|p| (p.value - new_t).abs() < 1e-12);
            self.selected = new_idx;
        }
    }
    /// Add a new control point at position `t` sampled from the current TF color at that position.
    pub fn add_at(&mut self, t: f64) {
        let color = self.tf.sample(t);
        self.tf.add_point(t, color);
    }
    /// Remove the currently selected control point.
    ///
    /// Does nothing if no point is selected.
    pub fn remove_selected(&mut self) {
        if let Some(idx) = self.selected.take()
            && idx < self.tf.points().len()
        {
            let val = self.tf.points()[idx].value;
            self.tf.remove_near(val, 1e-12);
        }
    }
}
/// A pre-integration table that maps `(s_back, s_front)` pairs to RGBA.
///
/// Pre-integration avoids aliasing artifacts by analytically integrating
/// the transfer function over each ray segment.  The table is indexed by
/// a discretized `n × n` grid over `[0, 1]²`.
#[derive(Debug, Clone)]
pub struct PreintegrationTable {
    /// Number of entries per axis.
    pub n: usize,
    /// Flat `n × n` table indexed `[s_back * n + s_front]`.
    pub table: Vec<Rgba>,
}
impl PreintegrationTable {
    /// Build the pre-integration table from a [`TransferFunction`].
    ///
    /// The table is computed by numerically integrating the TF over segments.
    /// `n` is the resolution per axis; `sub_steps` controls integration accuracy.
    pub fn build(tf: &TransferFunction, n: usize, sub_steps: usize) -> Self {
        let n = n.max(2);
        let sub = sub_steps.max(4);
        let mut table = Vec::with_capacity(n * n);
        let sub_f = sub as f64;
        for ib in 0..n {
            for _if in 0..n {
                let sb = ib as f64 / (n - 1) as f64;
                let sf = _if as f64 / (n - 1) as f64;
                let lo = sb.min(sf);
                let hi = sb.max(sf);
                if (hi - lo).abs() < 1e-12 {
                    table.push(tf.sample(lo));
                    continue;
                }
                let mut acc_rgb = [0.0_f64; 3];
                let mut acc_alpha = 0.0_f64;
                for s in 0..sub {
                    let t0 = lo + (hi - lo) * s as f64 / sub_f;
                    let t1 = lo + (hi - lo) * (s + 1) as f64 / sub_f;
                    let t_mid = (t0 + t1) / 2.0;
                    let dt = t1 - t0;
                    let c = tf.sample(t_mid);
                    let a = 1.0 - (-c.a * dt * 100.0).exp();
                    let trans = 1.0 - acc_alpha;
                    acc_rgb[0] += trans * a * c.r;
                    acc_rgb[1] += trans * a * c.g;
                    acc_rgb[2] += trans * a * c.b;
                    acc_alpha += trans * a;
                    if acc_alpha >= 0.99 {
                        break;
                    }
                }
                table.push(Rgba::new(
                    acc_rgb[0].clamp(0.0, 1.0),
                    acc_rgb[1].clamp(0.0, 1.0),
                    acc_rgb[2].clamp(0.0, 1.0),
                    acc_alpha.clamp(0.0, 1.0),
                ));
            }
        }
        Self { n, table }
    }
    /// Look up the pre-integrated color for the segment `[s_back, s_front]`.
    pub fn sample(&self, s_back: f64, s_front: f64) -> Rgba {
        let s_back = s_back.clamp(0.0, 1.0);
        let s_front = s_front.clamp(0.0, 1.0);
        let ib = ((s_back * (self.n - 1) as f64).round() as usize).min(self.n - 1);
        let ifr = ((s_front * (self.n - 1) as f64).round() as usize).min(self.n - 1);
        self.table[ib * self.n + ifr]
    }
}
/// A transfer function that decouples color (from a [`GradientColormap`])
/// and opacity (from an [`OpacityProfile`]).
///
/// This design is common in advanced volume rendering UIs where artists control
/// color hue and transparency independently.
#[derive(Debug, Clone)]
pub struct TwoPartTransferFunction {
    /// Color source.
    pub colormap: GradientColormap,
    /// Opacity source.
    pub opacity: OpacityProfile,
    /// Global opacity scale applied on top of the profile.
    pub global_opacity: f64,
}
impl TwoPartTransferFunction {
    /// Construct a two-part TF.
    pub fn new(colormap: GradientColormap, opacity: OpacityProfile) -> Self {
        Self {
            colormap,
            opacity,
            global_opacity: 1.0,
        }
    }
    /// Set the global opacity multiplier.
    pub fn with_global_opacity(mut self, alpha: f64) -> Self {
        self.global_opacity = alpha.clamp(0.0, 1.0);
        self
    }
    /// Sample at `t` ∈ \[0, 1\], returning RGBA.
    pub fn sample(&self, t: f64) -> Rgba {
        let [r, g, b] = self.colormap.sample(t);
        let a = (self.opacity.evaluate(t) * self.global_opacity).clamp(0.0, 1.0);
        Rgba::new(r, g, b, a)
    }
    /// Build a full lookup table of `n` entries.
    pub fn build_lut(&self, n: usize) -> Vec<Rgba> {
        (0..n)
            .map(|i| {
                let t = i as f64 / (n.saturating_sub(1).max(1)) as f64;
                self.sample(t)
            })
            .collect()
    }
    /// Convert to a plain [`TransferFunction`] by sampling at `resolution` points.
    pub fn to_transfer_function(&self, resolution: usize) -> TransferFunction {
        let mut tf = TransferFunction::new();
        let n = resolution.max(2);
        for i in 0..n {
            let t = i as f64 / (n - 1) as f64;
            tf.add_point(t, self.sample(t));
        }
        tf
    }
    /// Pre-built X-ray-style: gray colormap, opacity rises steeply for dense material.
    pub fn xray_style() -> Self {
        let mut cm = GradientColormap::new();
        cm.add_stop(ColorStop::new(0.0, 0.0, 0.0, 0.0));
        cm.add_stop(ColorStop::new(1.0, 1.0, 1.0, 1.0));
        let op = OpacityProfile::new(OpacityProfileKind::Gaussian {
            center: 0.75,
            sigma: 0.15,
        });
        Self::new(cm, op).with_global_opacity(0.9)
    }
    /// Pre-built fire-style: cool → hot colors with sigmoid opacity ramp.
    pub fn fire_style() -> Self {
        let mut cm = GradientColormap::new();
        cm.add_stop(ColorStop::new(0.0, 0.0, 0.0, 0.0));
        cm.add_stop(ColorStop::new(0.25, 0.5, 0.0, 0.0));
        cm.add_stop(ColorStop::new(0.6, 1.0, 0.4, 0.0));
        cm.add_stop(ColorStop::new(0.85, 1.0, 1.0, 0.2));
        cm.add_stop(ColorStop::new(1.0, 1.0, 1.0, 1.0));
        let op = OpacityProfile::new(OpacityProfileKind::LinearRamp).with_scale(0.85);
        Self::new(cm, op)
    }
}
/// A lookup table that assigns material classes to normalized scalar values.
///
/// Useful for CT/MRI volume rendering where distinct materials (bone, tissue,
/// air) occupy identifiable Hounsfield-unit ranges.
#[derive(Debug, Clone, Default)]
pub struct ClassificationTable {
    pub(super) classes: Vec<MaterialClass>,
}
impl ClassificationTable {
    /// Create an empty classification table.
    pub fn new() -> Self {
        Self {
            classes: Vec::new(),
        }
    }
    /// Register a material class.
    pub fn add_class(&mut self, class: MaterialClass) {
        self.classes.push(class);
    }
    /// Find the first material class containing `t`, if any.
    pub fn classify(&self, t: f64) -> Option<&MaterialClass> {
        self.classes.iter().find(|c| c.contains(t))
    }
    /// Classify `t` and return an RGBA color.
    ///
    /// Returns transparent black if no class matches.
    pub fn sample(&self, t: f64) -> Rgba {
        if let Some(cls) = self.classify(t) {
            let mut c = cls.color;
            c.a = (c.a * cls.opacity_scale).clamp(0.0, 1.0);
            c
        } else {
            Rgba::transparent()
        }
    }
    /// Convert the classification table into a [`TransferFunction`] with
    /// `resolution` evenly-spaced samples.
    pub fn to_transfer_function(&self, resolution: usize) -> TransferFunction {
        let mut tf = TransferFunction::new();
        let n = resolution.max(2);
        for i in 0..n {
            let t = i as f64 / (n - 1) as f64;
            tf.add_point(t, self.sample(t));
        }
        tf
    }
    /// Pre-built CT classification for 4 tissues (normalized to \[0, 1\]).
    ///
    /// Approximate Hounsfield normalization:
    /// - Air: \[0.0, 0.1\] → transparent gray
    /// - Soft tissue: \[0.3, 0.6\] → semi-transparent pink
    /// - Bone: \[0.7, 1.0\] → opaque white/yellow
    pub fn ct_tissues() -> Self {
        let mut ct = Self::new();
        ct.add_class(MaterialClass::new(
            "air",
            0.0,
            0.1,
            Rgba::new(0.5, 0.5, 0.5, 0.01),
            0.0,
        ));
        ct.add_class(MaterialClass::new(
            "soft_tissue",
            0.3,
            0.6,
            Rgba::new(0.9, 0.6, 0.55, 0.4),
            0.6,
        ));
        ct.add_class(MaterialClass::new(
            "bone",
            0.7,
            1.0,
            Rgba::new(0.95, 0.95, 0.8, 1.0),
            1.0,
        ));
        ct
    }
}
/// Quality presets for the ray-marching step size.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RenderQuality {
    /// Coarse — fast preview, `step = 0.05`.
    Preview,
    /// Medium — balanced quality/performance, `step = 0.01`.
    Medium,
    /// High — production quality, `step = 0.002`.
    High,
}
impl RenderQuality {
    /// Step size in normalized \[0, 1\] space for this quality level.
    pub fn step_size(&self) -> f64 {
        match self {
            RenderQuality::Preview => 0.05,
            RenderQuality::Medium => 0.01,
            RenderQuality::High => 0.002,
        }
    }
}
/// A standalone piecewise-linear opacity map with explicit knots.
///
/// Similar to [`OpacityProfile`] but with user-specified knots at arbitrary
/// positions, allowing fine-grained manual opacity editing.
#[derive(Debug, Clone)]
pub struct PiecewiseLinearOpacityMap {
    /// Sorted list of `(position, opacity)` knots.
    pub(super) knots: Vec<(f64, f64)>,
}
impl PiecewiseLinearOpacityMap {
    /// Create an empty opacity map.
    pub fn new() -> Self {
        Self { knots: Vec::new() }
    }
    /// Add or replace a knot at `position`.
    pub fn add_knot(&mut self, position: f64, opacity: f64) {
        let pos = position.clamp(0.0, 1.0);
        let op = opacity.clamp(0.0, 1.0);
        let idx = self.knots.partition_point(|&(p, _)| p < pos);
        if idx < self.knots.len() && (self.knots[idx].0 - pos).abs() < 1e-12 {
            self.knots[idx] = (pos, op);
        } else {
            self.knots.insert(idx, (pos, op));
        }
    }
    /// Number of knots.
    pub fn knot_count(&self) -> usize {
        self.knots.len()
    }
    /// Evaluate opacity at `t` in `[0, 1]`.
    pub fn evaluate(&self, t: f64) -> f64 {
        let t = t.clamp(0.0, 1.0);
        if self.knots.is_empty() {
            return 0.0;
        }
        if self.knots.len() == 1 {
            return self.knots[0].1;
        }
        if t <= self.knots[0].0 {
            return self.knots[0].1;
        }
        if t >= self.knots[self.knots.len() - 1].0 {
            return self.knots[self.knots.len() - 1].1;
        }
        let hi = self.knots.partition_point(|&(p, _)| p <= t);
        let lo = hi.saturating_sub(1).min(self.knots.len() - 2);
        let (p0, o0) = self.knots[lo];
        let (p1, o1) = self.knots[lo + 1];
        let span = p1 - p0;
        let local = if span < 1e-15 { 0.5 } else { (t - p0) / span };
        o0 + (o1 - o0) * local
    }
    /// Combine with a [`GradientColormap`] to produce a full RGBA TF.
    pub fn to_transfer_function(&self, colormap: &GradientColormap) -> TransferFunction {
        let mut tf = TransferFunction::new();
        for &(pos, op) in &self.knots {
            let [r, g, b] = colormap.sample(pos);
            tf.add_point(pos, Rgba::new(r, g, b, op));
        }
        tf
    }
    /// Build a sampled lookup table of opacity values.
    pub fn build_lut(&self, n: usize) -> Vec<f64> {
        (0..n)
            .map(|i| {
                let t = i as f64 / (n.saturating_sub(1).max(1)) as f64;
                self.evaluate(t)
            })
            .collect()
    }
}
/// A single color stop in a gradient colormap.
#[derive(Debug, Clone, Copy)]
pub struct ColorStop {
    /// Position of this stop in \[0, 1\].
    pub position: f64,
    /// RGB color at this stop (alpha ignored; set by opacity profile).
    pub rgb: [f64; 3],
}
impl ColorStop {
    /// Construct a color stop.
    pub fn new(position: f64, r: f64, g: f64, b: f64) -> Self {
        Self {
            position: position.clamp(0.0, 1.0),
            rgb: [r, g, b],
        }
    }
}
/// A material class assigned to a scalar interval.
#[derive(Debug, Clone)]
pub struct MaterialClass {
    /// Name of this material class (e.g., `"bone"`, `"soft_tissue"`).
    pub name: String,
    /// Normalized scalar range: `[lo, hi]` in \[0, 1\].
    pub range: [f64; 2],
    /// Base color for this class.
    pub color: Rgba,
    /// Opacity multiplier for this class.
    pub opacity_scale: f64,
}
impl MaterialClass {
    /// Construct a material class.
    pub fn new(name: impl Into<String>, lo: f64, hi: f64, color: Rgba, opacity_scale: f64) -> Self {
        Self {
            name: name.into(),
            range: [lo.clamp(0.0, 1.0), hi.clamp(0.0, 1.0)],
            color,
            opacity_scale: opacity_scale.clamp(0.0, 1.0),
        }
    }
    /// Returns `true` if `t` falls within this class's range.
    pub fn contains(&self, t: f64) -> bool {
        t >= self.range[0] && t <= self.range[1]
    }
}
/// A standalone 1-D opacity (alpha) profile.
#[derive(Debug, Clone, Copy)]
pub struct OpacityProfile {
    /// The profile shape.
    pub kind: OpacityProfileKind,
    /// Global opacity scale factor applied after profile evaluation.
    pub scale: f64,
}
impl OpacityProfile {
    /// Create a new opacity profile.
    pub fn new(kind: OpacityProfileKind) -> Self {
        Self { kind, scale: 1.0 }
    }
    /// Apply a global opacity scale.
    pub fn with_scale(mut self, scale: f64) -> Self {
        self.scale = scale.clamp(0.0, 1.0);
        self
    }
    /// Evaluate opacity at normalized value `t` in \[0, 1\].
    pub fn evaluate(&self, t: f64) -> f64 {
        let t = t.clamp(0.0, 1.0);
        let raw = match self.kind {
            OpacityProfileKind::LinearRamp => t,
            OpacityProfileKind::Constant(c) => c.clamp(0.0, 1.0),
            OpacityProfileKind::Gaussian { center, sigma } => {
                if sigma < 1e-30 {
                    return 0.0;
                }
                let x = (t - center) / sigma;
                (-0.5 * x * x).exp()
            }
            OpacityProfileKind::Step { threshold } => {
                if t >= threshold {
                    1.0
                } else {
                    0.0
                }
            }
            OpacityProfileKind::Tent { center, half_width } => {
                if half_width < 1e-30 {
                    return 0.0;
                }
                let d = (t - center).abs() / half_width;
                (1.0 - d).clamp(0.0, 1.0)
            }
        };
        (raw * self.scale).clamp(0.0, 1.0)
    }
}
/// A multi-stop gradient colormap that maps scalar values to RGB.
///
/// Unlike [`TransferFunction`], `GradientColormap` deals only with color,
/// not opacity.  Combine with an [`OpacityProfile`] for full volume rendering.
#[derive(Debug, Clone)]
pub struct GradientColormap {
    pub(super) stops: Vec<ColorStop>,
}
impl GradientColormap {
    /// Create an empty gradient colormap.
    pub fn new() -> Self {
        Self { stops: Vec::new() }
    }
    /// Add or replace a color stop at `stop.position`.
    pub fn add_stop(&mut self, stop: ColorStop) {
        let pos = self.stops.partition_point(|s| s.position < stop.position);
        if pos < self.stops.len() && (self.stops[pos].position - stop.position).abs() < 1e-12 {
            self.stops[pos] = stop;
        } else {
            self.stops.insert(pos, stop);
        }
    }
    /// Number of stops.
    pub fn stop_count(&self) -> usize {
        self.stops.len()
    }
    /// Sample the gradient at position `t` in \[0, 1\].
    ///
    /// Returns black if no stops are defined.
    pub fn sample(&self, t: f64) -> [f64; 3] {
        let t = t.clamp(0.0, 1.0);
        if self.stops.is_empty() {
            return [0.0; 3];
        }
        if self.stops.len() == 1 {
            return self.stops[0].rgb;
        }
        if t <= self.stops[0].position {
            return self.stops[0].rgb;
        }
        if t >= self.stops[self.stops.len() - 1].position {
            return self.stops[self.stops.len() - 1].rgb;
        }
        let hi = self.stops.partition_point(|s| s.position <= t);
        let lo = hi.saturating_sub(1).min(self.stops.len() - 2);
        let hi = lo + 1;
        let a = &self.stops[lo];
        let b = &self.stops[hi];
        let span = b.position - a.position;
        let local_t = if span < 1e-15 {
            0.5
        } else {
            (t - a.position) / span
        };
        [
            a.rgb[0] + (b.rgb[0] - a.rgb[0]) * local_t,
            a.rgb[1] + (b.rgb[1] - a.rgb[1]) * local_t,
            a.rgb[2] + (b.rgb[2] - a.rgb[2]) * local_t,
        ]
    }
    /// Sample and combine with an opacity profile to produce a full RGBA value.
    pub fn sample_with_opacity(&self, t: f64, opacity: &OpacityProfile) -> Rgba {
        let [r, g, b] = self.sample(t);
        Rgba::new(r, g, b, opacity.evaluate(t))
    }
    /// Build a lookup table of `n` evenly-spaced RGB samples.
    pub fn build_rgb_lut(&self, n: usize) -> Vec<[f64; 3]> {
        (0..n)
            .map(|i| {
                let t = i as f64 / (n.saturating_sub(1).max(1)) as f64;
                self.sample(t)
            })
            .collect()
    }
    /// Diverging blue–white–red colormap (useful for signed quantities).
    pub fn diverging_bwr() -> Self {
        let mut g = Self::new();
        g.add_stop(ColorStop::new(0.0, 0.0, 0.2, 0.8));
        g.add_stop(ColorStop::new(0.5, 1.0, 1.0, 1.0));
        g.add_stop(ColorStop::new(1.0, 0.8, 0.0, 0.0));
        g
    }
    /// Sequential blue-green-yellow colormap (similar to YlGnBu reversed).
    pub fn sequential_bgy() -> Self {
        let mut g = Self::new();
        g.add_stop(ColorStop::new(0.0, 0.0, 0.0, 0.5));
        g.add_stop(ColorStop::new(0.4, 0.0, 0.6, 0.4));
        g.add_stop(ColorStop::new(1.0, 0.9, 0.9, 0.1));
        g
    }
    /// Inferno-inspired: black → dark purple → orange → yellow.
    pub fn inferno_approx() -> Self {
        let mut g = Self::new();
        g.add_stop(ColorStop::new(0.0, 0.0, 0.0, 0.0));
        g.add_stop(ColorStop::new(0.25, 0.2, 0.0, 0.4));
        g.add_stop(ColorStop::new(0.5, 0.6, 0.1, 0.2));
        g.add_stop(ColorStop::new(0.75, 0.95, 0.55, 0.05));
        g.add_stop(ColorStop::new(1.0, 0.99, 0.99, 0.65));
        g
    }
    /// Rainbow colormap: violet → blue → cyan → green → yellow → red.
    pub fn rainbow() -> Self {
        let mut g = Self::new();
        g.add_stop(ColorStop::new(0.0, 0.5, 0.0, 1.0));
        g.add_stop(ColorStop::new(0.166, 0.0, 0.0, 1.0));
        g.add_stop(ColorStop::new(0.333, 0.0, 1.0, 1.0));
        g.add_stop(ColorStop::new(0.5, 0.0, 1.0, 0.0));
        g.add_stop(ColorStop::new(0.666, 1.0, 1.0, 0.0));
        g.add_stop(ColorStop::new(0.833, 1.0, 0.5, 0.0));
        g.add_stop(ColorStop::new(1.0, 1.0, 0.0, 0.0));
        g
    }
}
/// A simple fixed-bucket histogram for scalar data analysis.
///
/// Used to automatically place control points at data-dense or interesting
/// regions of the value range.
#[derive(Debug, Clone)]
pub struct HistogramAnalyzer {
    /// Histogram bucket counts.
    pub(super) counts: Vec<u64>,
    /// Minimum scalar value of the data.
    pub data_min: f64,
    /// Maximum scalar value of the data.
    pub data_max: f64,
}
impl HistogramAnalyzer {
    /// Build a histogram from raw scalar data with `num_buckets` bins.
    ///
    /// Returns `None` if `data` is empty or all values are identical.
    pub fn from_data(data: &[f64], num_buckets: usize) -> Option<Self> {
        if data.is_empty() || num_buckets == 0 {
            return None;
        }
        let data_min = data.iter().cloned().fold(f64::INFINITY, f64::min);
        let data_max = data.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        if (data_max - data_min).abs() < 1e-30 {
            return None;
        }
        let mut counts = vec![0u64; num_buckets];
        let range = data_max - data_min;
        for &v in data {
            let t = ((v - data_min) / range).clamp(0.0, 1.0 - 1e-12);
            let idx = (t * num_buckets as f64) as usize;
            counts[idx.min(num_buckets - 1)] += 1;
        }
        Some(Self {
            counts,
            data_min,
            data_max,
        })
    }
    /// Number of buckets.
    pub fn num_buckets(&self) -> usize {
        self.counts.len()
    }
    /// Raw count in bucket `i`.
    pub fn count(&self, i: usize) -> u64 {
        self.counts.get(i).copied().unwrap_or(0)
    }
    /// Normalized frequency (fraction of total) in bucket `i`.
    pub fn frequency(&self, i: usize) -> f64 {
        let total: u64 = self.counts.iter().sum();
        if total == 0 {
            return 0.0;
        }
        self.count(i) as f64 / total as f64
    }
    /// Normalized center value (in \[0, 1\]) of bucket `i`.
    pub fn bucket_center_normalized(&self, i: usize) -> f64 {
        let n = self.counts.len();
        if n == 0 {
            return 0.0;
        }
        (i as f64 + 0.5) / n as f64
    }
    /// Raw scalar center of bucket `i`.
    pub fn bucket_center_raw(&self, i: usize) -> f64 {
        self.data_min + self.bucket_center_normalized(i) * (self.data_max - self.data_min)
    }
    /// Find the index of the bucket with the maximum count.
    pub fn peak_bucket(&self) -> usize {
        self.counts
            .iter()
            .enumerate()
            .max_by_key(|&(_, c)| c)
            .map(|(i, _)| i)
            .unwrap_or(0)
    }
    /// Compute the cumulative distribution at bucket `i` (inclusive).
    pub fn cdf(&self, i: usize) -> f64 {
        let total: u64 = self.counts.iter().sum();
        if total == 0 {
            return 0.0;
        }
        let cum: u64 = self.counts[..=i.min(self.counts.len() - 1)].iter().sum();
        cum as f64 / total as f64
    }
    /// Find the normalized value at which `fraction` of the data lies below it.
    ///
    /// E.g. `percentile(0.95)` returns the 95th-percentile position in \[0, 1\].
    pub fn percentile(&self, fraction: f64) -> f64 {
        let fraction = fraction.clamp(0.0, 1.0);
        for i in 0..self.counts.len() {
            if self.cdf(i) >= fraction {
                return self.bucket_center_normalized(i);
            }
        }
        1.0
    }
    /// Suggest an automatic [`TransferFunction`] based on the histogram.
    ///
    /// The algorithm places sparse/transparent control points at low-frequency
    /// regions and opaque control points near the data peak.
    pub fn auto_transfer_function(&self, colormap: &GradientColormap) -> TransferFunction {
        let mut tf = TransferFunction::new();
        let peak = self.peak_bucket();
        let n = self.counts.len();
        for i in 0..n {
            let t = self.bucket_center_normalized(i);
            let freq = self.frequency(i);
            let dist_from_peak = (i as f64 - peak as f64).abs() / n as f64;
            let opacity = (freq * 10.0).min(1.0) * (1.0 - dist_from_peak * 2.0).max(0.0);
            let [r, g, b] = colormap.sample(t);
            tf.add_point(t, Rgba::new(r, g, b, opacity));
        }
        tf
    }
}
/// Parameters controlling a CPU direct-volume-rendering pass.
#[derive(Debug, Clone)]
pub struct VolumeRenderingParams {
    /// Transfer function to apply at each sample.
    pub transfer_function: TransferFunction,
    /// Ray-march step size (in normalized \[0, 1\] space).
    pub step_size: f64,
    /// Early termination threshold: stop accumulating when opacity exceeds this.
    pub early_termination_threshold: f64,
    /// Number of ray samples (derived from `step_size` and unit-cube path length).
    pub max_samples: usize,
    /// Enable gradient-magnitude-based opacity shading.
    pub gradient_shading: bool,
    /// Ambient lighting coefficient for gradient shading.
    pub ambient: f64,
    /// Diffuse lighting coefficient for gradient shading.
    pub diffuse: f64,
    /// Empty-space skipping threshold: skip voxels with opacity below this.
    pub empty_space_threshold: f64,
}
impl VolumeRenderingParams {
    /// Construct params with a specific quality preset.
    pub fn with_quality(quality: RenderQuality) -> Self {
        let step = quality.step_size();
        let max_samples = (1.0 / step * 1.74) as usize;
        Self {
            step_size: step,
            max_samples: max_samples.max(8),
            ..Self::default()
        }
    }
    /// Apply gradient shading settings.
    pub fn with_lighting(mut self, ambient: f64, diffuse: f64) -> Self {
        self.ambient = ambient.clamp(0.0, 1.0);
        self.diffuse = diffuse.clamp(0.0, 1.0);
        self
    }
    /// Set the empty-space skipping threshold.
    pub fn with_empty_space_threshold(mut self, threshold: f64) -> Self {
        self.empty_space_threshold = threshold.clamp(0.0, 1.0);
        self
    }
}
/// A 64-bit floating-point RGBA color.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Rgba {
    /// Red channel in \[0, 1\].
    pub r: f64,
    /// Green channel in \[0, 1\].
    pub g: f64,
    /// Blue channel in \[0, 1\].
    pub b: f64,
    /// Alpha (opacity) channel in \[0, 1\].
    pub a: f64,
}
impl Rgba {
    /// Construct an RGBA color from components.
    pub fn new(r: f64, g: f64, b: f64, a: f64) -> Self {
        Self { r, g, b, a }
    }
    /// Opaque black.
    pub fn black() -> Self {
        Self::new(0.0, 0.0, 0.0, 1.0)
    }
    /// Opaque white.
    pub fn white() -> Self {
        Self::new(1.0, 1.0, 1.0, 1.0)
    }
    /// Transparent black.
    pub fn transparent() -> Self {
        Self::new(0.0, 0.0, 0.0, 0.0)
    }
    /// Linear interpolation between two colors.
    pub fn lerp(a: Rgba, b: Rgba, t: f64) -> Rgba {
        let t = t.clamp(0.0, 1.0);
        Rgba {
            r: a.r + (b.r - a.r) * t,
            g: a.g + (b.g - a.g) * t,
            b: a.b + (b.b - a.b) * t,
            a: a.a + (b.a - a.a) * t,
        }
    }
    /// Multiply color channels by a scalar (alpha unchanged).
    pub fn scale_rgb(self, s: f64) -> Self {
        Rgba::new(
            (self.r * s).clamp(0.0, 1.0),
            (self.g * s).clamp(0.0, 1.0),
            (self.b * s).clamp(0.0, 1.0),
            self.a,
        )
    }
    /// Pre-multiply alpha: RGB = RGB * alpha.
    pub fn premultiply_alpha(self) -> Self {
        Rgba::new(self.r * self.a, self.g * self.a, self.b * self.a, self.a)
    }
    /// Convert to an `[f32; 4]` array (RGBA order).
    pub fn to_f32_array(self) -> [f32; 4] {
        [self.r as f32, self.g as f32, self.b as f32, self.a as f32]
    }
    /// Convert to an `[u8; 4]` array (RGBA order, gamma-linear to 0..=255).
    pub fn to_u8_array(self) -> [u8; 4] {
        fn to_u8(v: f64) -> u8 {
            (v.clamp(0.0, 1.0) * 255.0).round() as u8
        }
        [to_u8(self.r), to_u8(self.g), to_u8(self.b), to_u8(self.a)]
    }
}
/// A 2-D transfer function keyed on `(density, gradient_magnitude)`, both
/// normalized to `[0, 1]`.
///
/// The table is stored as a flat `rows × cols` grid.  Bilinear interpolation
/// is used at lookup time.
#[derive(Debug, Clone)]
pub struct MultiDimensionalTransferFunction {
    /// Number of columns (density axis).
    pub cols: usize,
    /// Number of rows (gradient magnitude axis).
    pub rows: usize,
    /// Flat row-major array of `rows * cols` RGBA entries.
    pub table: Vec<Rgba>,
}
impl MultiDimensionalTransferFunction {
    /// Create a new 2-D TF filled with transparent black.
    pub fn new(cols: usize, rows: usize) -> Self {
        let cols = cols.max(1);
        let rows = rows.max(1);
        Self {
            cols,
            rows,
            table: vec![Rgba::transparent(); cols * rows],
        }
    }
    /// Set the RGBA value at `(density_idx, grad_idx)`.
    pub fn set(&mut self, density_idx: usize, grad_idx: usize, color: Rgba) {
        let i = grad_idx.min(self.rows - 1) * self.cols + density_idx.min(self.cols - 1);
        self.table[i] = color;
    }
    /// Sample the 2-D TF with bilinear interpolation.
    ///
    /// Both `density` and `grad_mag` should be in `[0, 1]`.
    pub fn sample(&self, density: f64, grad_mag: f64) -> Rgba {
        let density = density.clamp(0.0, 1.0);
        let grad_mag = grad_mag.clamp(0.0, 1.0);
        let fc = density * (self.cols - 1) as f64;
        let fr = grad_mag * (self.rows - 1) as f64;
        let c0 = fc.floor() as usize;
        let r0 = fr.floor() as usize;
        let c1 = (c0 + 1).min(self.cols - 1);
        let r1 = (r0 + 1).min(self.rows - 1);
        let tc = fc - c0 as f64;
        let tr = fr - r0 as f64;
        let idx = |r: usize, c: usize| r * self.cols + c;
        let c00 = self.table[idx(r0, c0)];
        let c01 = self.table[idx(r0, c1)];
        let c10 = self.table[idx(r1, c0)];
        let c11 = self.table[idx(r1, c1)];
        let lo = Rgba::lerp(c00, c01, tc);
        let hi = Rgba::lerp(c10, c11, tc);
        Rgba::lerp(lo, hi, tr)
    }
    /// Convert to a flat 1-D [`TransferFunction`] by averaging over all gradient
    /// magnitude rows.
    pub fn flatten_to_1d(&self) -> TransferFunction {
        let mut tf = TransferFunction::new();
        for ci in 0..self.cols {
            let t = ci as f64 / (self.cols - 1).max(1) as f64;
            let mut acc = Rgba::new(0.0, 0.0, 0.0, 0.0);
            for ri in 0..self.rows {
                let c = self.table[ri * self.cols + ci];
                acc.r += c.r;
                acc.g += c.g;
                acc.b += c.b;
                acc.a += c.a;
            }
            let n = self.rows as f64;
            tf.add_point(t, Rgba::new(acc.r / n, acc.g / n, acc.b / n, acc.a / n));
        }
        tf
    }
}
/// Computes opacity by combining a base scalar TF with gradient-magnitude boosting.
///
/// At each sample, a 1-D base opacity (from a [`PiecewiseLinearOpacityMap`]) is
/// multiplied by a factor that increases near boundaries (high gradient magnitude).
#[derive(Debug, Clone)]
pub struct GradientMagnitudeOpacity {
    /// Base opacity map (applied to normalized scalar value).
    pub base_map: PiecewiseLinearOpacityMap,
    /// Gradient sharpness: higher values give steeper boundary enhancement.
    pub sharpness: f64,
    /// Minimum opacity floor (prevents fully transparent regions).
    pub opacity_floor: f64,
}
impl GradientMagnitudeOpacity {
    /// Create with default settings.
    pub fn new(base_map: PiecewiseLinearOpacityMap) -> Self {
        Self {
            base_map,
            sharpness: 3.0,
            opacity_floor: 0.0,
        }
    }
    /// Set the sharpness parameter.
    pub fn with_sharpness(mut self, s: f64) -> Self {
        self.sharpness = s.max(0.0);
        self
    }
    /// Set the minimum opacity floor.
    pub fn with_floor(mut self, floor: f64) -> Self {
        self.opacity_floor = floor.clamp(0.0, 1.0);
        self
    }
    /// Evaluate opacity at scalar `t` with gradient magnitude `g` (both in \[0,1\]).
    pub fn evaluate(&self, t: f64, g: f64) -> f64 {
        let base = self.base_map.evaluate(t);
        let boost = (g * self.sharpness).clamp(0.0, 1.0);
        let op = (base + boost * (1.0 - base)).clamp(0.0, 1.0);
        op.max(self.opacity_floor)
    }
}
