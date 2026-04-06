//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::needless_range_loop, clippy::ptr_arg)]
#[allow(unused_imports)]
use super::functions::*;
/// Pressure field colourmap over a mesh.
#[allow(dead_code)]
pub struct PressureColormap {
    /// Minimum pressure for clamping.
    pub p_min: f64,
    /// Maximum pressure for clamping.
    pub p_max: f64,
    /// Number of iso-contour levels.
    pub contour_levels: usize,
    /// Colormap name (stored for reference).
    pub colormap_name: String,
}
impl PressureColormap {
    /// Create a pressure colormap.
    pub fn new(p_min: f64, p_max: f64) -> Self {
        Self {
            p_min,
            p_max,
            contour_levels: 10,
            colormap_name: "viridis".to_string(),
        }
    }
    /// Map a pressure value to an RGBA colour using a viridis-like palette.
    pub fn map(&self, pressure: f64) -> [f32; 4] {
        let t =
            ((pressure - self.p_min) / (self.p_max - self.p_min).max(1e-12)).clamp(0.0, 1.0) as f32;
        let r = (t * t * 1.5).min(1.0);
        let g = (t * 1.2 - 0.1).clamp(0.0, 1.0);
        let b = ((1.0 - t) * 0.9).clamp(0.0, 1.0);
        [r, g, b, 1.0]
    }
    /// Build coloured pressure samples for mesh vertices.
    pub fn apply(&self, positions: &[[f64; 3]], pressures: &[f64]) -> Vec<PressureSample> {
        positions
            .iter()
            .zip(pressures.iter())
            .map(|(&pos, &p)| PressureSample {
                position: pos,
                pressure: p,
                color: self.map(p),
            })
            .collect()
    }
    /// Compute iso-contour threshold values.
    pub fn iso_contours(&self) -> Vec<f64> {
        let n = self.contour_levels;
        (0..=n)
            .map(|i| self.p_min + (self.p_max - self.p_min) * i as f64 / n as f64)
            .collect()
    }
}
/// Extracts the free surface from a signed-distance level-set φ field.
#[allow(dead_code)]
pub struct FreeSurfaceExtractor {
    /// Grid resolution per axis.
    pub resolution: usize,
    /// Grid domain size.
    pub domain_size: f64,
}
impl FreeSurfaceExtractor {
    /// Create a free-surface extractor.
    pub fn new(resolution: usize, domain_size: f64) -> Self {
        Self {
            resolution,
            domain_size,
        }
    }
    /// Extract triangles from a signed-distance field φ (zero = free surface).
    ///
    /// Uses a simplified marching-cubes-like approach: for each edge
    /// where φ changes sign, insert a vertex at the crossing.
    pub fn extract(&self, phi: &[f64]) -> Vec<FreeSurfaceTriangle> {
        let n = self.resolution;
        let h = self.domain_size / n as f64;
        let idx = |x: usize, y: usize, z: usize| x * n * n + y * n + z;
        let mut tris = Vec::new();
        for ix in 0..n.saturating_sub(1) {
            for iy in 0..n.saturating_sub(1) {
                for iz in 0..n.saturating_sub(1) {
                    let v000 = phi[idx(ix, iy, iz)];
                    let v100 = phi[idx(ix + 1, iy, iz)];
                    if (v000 < 0.0) != (v100 < 0.0) {
                        let t = -v000 / (v100 - v000 + 1e-15);
                        let px = ix as f64 * h + t * h;
                        let py = iy as f64 * h;
                        let pz = iz as f64 * h;
                        let p = [px, py, pz];
                        let dp = h * 0.5;
                        let verts = [p, [px, py + dp, pz], [px, py, pz + dp]];
                        let e1 = sub3(verts[1], verts[0]);
                        let e2 = sub3(verts[2], verts[0]);
                        let normal = norm3(cross3(e1, e2));
                        tris.push(FreeSurfaceTriangle { verts, normal });
                    }
                }
            }
        }
        tris
    }
}
/// Flow statistics: mean, RMS, turbulence intensity, and Reynolds stress.
#[allow(dead_code)]
#[derive(Debug, Clone, Default)]
pub struct FlowStatistics {
    /// Mean velocity components `[u_mean, v_mean, w_mean]`.
    pub mean: [f64; 3],
    /// RMS (root-mean-square) of each velocity component.
    pub rms: [f64; 3],
    /// Turbulence intensity TI = sqrt(2k/3) / |U_mean|.
    pub turbulence_intensity: f64,
    /// Reynolds stress tensor components `[u'u', v'v', w'w', u'v', u'w', v'w']`.
    pub reynolds_stress: [f64; 6],
}
impl FlowStatistics {
    /// Compute statistics from a set of instantaneous velocity samples.
    pub fn compute(samples: &[[f64; 3]]) -> Self {
        let n = samples.len();
        if n == 0 {
            return Self::default();
        }
        let nf = n as f64;
        let mut mean = [0.0_f64; 3];
        for s in samples {
            for k in 0..3 {
                mean[k] += s[k];
            }
        }
        for k in 0..3 {
            mean[k] /= nf;
        }
        let mut rms = [0.0_f64; 3];
        let mut rs = [0.0_f64; 6];
        for s in samples {
            let u = s[0] - mean[0];
            let v = s[1] - mean[1];
            let w = s[2] - mean[2];
            rms[0] += u * u;
            rms[1] += v * v;
            rms[2] += w * w;
            rs[0] += u * u;
            rs[1] += v * v;
            rs[2] += w * w;
            rs[3] += u * v;
            rs[4] += u * w;
            rs[5] += v * w;
        }
        for k in 0..3 {
            rms[k] = (rms[k] / nf).sqrt();
        }
        for k in 0..6 {
            rs[k] /= nf;
        }
        let tke = 0.5 * (rs[0] + rs[1] + rs[2]);
        let u_mag = len3(mean);
        let ti = if u_mag > 1e-12 {
            (2.0 * tke / 3.0).sqrt() / u_mag
        } else {
            0.0
        };
        Self {
            mean,
            rms,
            turbulence_intensity: ti,
            reynolds_stress: rs,
        }
    }
    /// Turbulent kinetic energy TKE = 0.5 * (u'² + v'² + w'²).
    pub fn tke(&self) -> f64 {
        0.5 * (self.reynolds_stress[0] + self.reynolds_stress[1] + self.reynolds_stress[2])
    }
}
/// A simple foam renderer that determines which particles should spawn foam.
#[allow(dead_code)]
pub struct FoamRenderer {
    /// Velocity threshold: particles above this speed spawn foam.
    pub velocity_threshold: f64,
    /// Curvature threshold for foam generation.
    pub curvature_threshold: f64,
}
impl FoamRenderer {
    /// Create a foam renderer with given thresholds.
    pub fn new(velocity_threshold: f64, curvature_threshold: f64) -> Self {
        Self {
            velocity_threshold,
            curvature_threshold,
        }
    }
    /// Decide whether a particle with given speed should spawn foam.
    pub fn should_foam(&self, speed: f64) -> bool {
        speed >= self.velocity_threshold
    }
    /// Count foam particles from a velocity field.
    pub fn count_foam_particles(&self, speeds: &[f64]) -> usize {
        speeds.iter().filter(|&&s| self.should_foam(s)).count()
    }
    /// Compute foam intensity as fraction of particles that foam.
    pub fn foam_intensity(&self, speeds: &[f64]) -> f64 {
        let n = speeds.len();
        if n == 0 {
            return 0.0;
        }
        self.count_foam_particles(speeds) as f64 / n as f64
    }
}
/// A coloured pressure sample on a mesh vertex.
#[allow(dead_code)]
pub struct PressureSample {
    /// 3D position.
    pub position: [f64; 3],
    /// Pressure value.
    pub pressure: f64,
    /// Mapped colour.
    pub color: [f32; 4],
}
/// A single rendered bubble.
#[allow(dead_code)]
pub struct Bubble {
    /// World-space centre.
    pub center: [f64; 3],
    /// Bubble radius.
    pub radius: f64,
    /// Refraction index.
    pub ior: f64,
    /// Opacity (0 = fully transparent).
    pub opacity: f32,
    /// Caustic intensity hint.
    pub caustic_intensity: f32,
}
impl Bubble {
    /// Create a bubble.
    pub fn new(center: [f64; 3], radius: f64) -> Self {
        Self {
            center,
            radius,
            ior: 1.0003,
            opacity: 0.15,
            caustic_intensity: 0.3,
        }
    }
    /// Compute refracted ray direction (Snell's law, simplified).
    ///
    /// `incident` is the unit incoming ray; `normal` is the surface normal.
    /// Returns the refracted unit vector, or `None` on total internal reflection.
    pub fn refract(&self, incident: [f64; 3], normal: [f64; 3], n1: f64) -> Option<[f64; 3]> {
        let n = n1 / self.ior;
        let cos_i = -dot3(incident, normal);
        let sin2_t = n * n * (1.0 - cos_i * cos_i);
        if sin2_t > 1.0 {
            return None;
        }
        let cos_t = (1.0 - sin2_t).sqrt();
        let refracted = add3(scale3(incident, n), scale3(normal, n * cos_i - cos_t));
        Some(norm3(refracted))
    }
}
/// Renders a streakline: positions of particles that have passed through
/// a fixed injection point over time.
#[allow(dead_code)]
pub struct StreaklineRenderer {
    /// Maximum length of the streak (number of positions).
    pub max_length: usize,
    /// Stored positions.
    pub positions: Vec<[f64; 3]>,
    /// Color gradient (interpolated along streak).
    pub color_head: [f32; 4],
    /// Tail color.
    pub color_tail: [f32; 4],
}
impl StreaklineRenderer {
    /// Create a new streakline renderer.
    pub fn new(max_length: usize) -> Self {
        Self {
            max_length,
            positions: Vec::new(),
            color_head: [1.0, 1.0, 1.0, 1.0],
            color_tail: [0.2, 0.2, 1.0, 0.3],
        }
    }
    /// Emit a new particle position from the injection point.
    pub fn emit(&mut self, pos: [f64; 3]) {
        self.positions.push(pos);
        while self.positions.len() > self.max_length {
            self.positions.remove(0);
        }
    }
    /// Compute the color at position index `i` along the streak.
    pub fn color_at(&self, i: usize) -> [f32; 4] {
        let n = self.positions.len().max(1);
        let t = i as f32 / (n - 1).max(1) as f32;
        lerp_color(self.color_tail, self.color_head, t)
    }
}
/// A triangle of the extracted free surface.
#[allow(dead_code)]
pub struct FreeSurfaceTriangle {
    /// Vertex positions.
    pub verts: [[f64; 3]; 3],
    /// Normal (outward from fluid).
    pub normal: [f64; 3],
}
/// A single velocity glyph at a sample point.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct VelocityGlyph {
    /// Sample origin.
    pub origin: [f64; 3],
    /// Tip of the glyph (for Arrow/Hedgehog style).
    pub tip: [f64; 3],
    /// Colour based on magnitude (RGBA, 0..1).
    pub color: [f32; 4],
    /// Glyph display style.
    pub style: GlyphStyle,
    /// Velocity magnitude at this glyph.
    pub magnitude: f64,
}
impl VelocityGlyph {
    /// Build a glyph from a velocity vector, scale, and reference magnitude for colour mapping.
    pub fn from_velocity(
        origin: [f64; 3],
        velocity: [f64; 3],
        scale: f64,
        ref_magnitude: f64,
        style: GlyphStyle,
    ) -> Self {
        let mag = len3(velocity);
        let tip = add3(origin, scale3(velocity, scale));
        let t = if ref_magnitude > 1e-12 {
            (mag / ref_magnitude).clamp(0.0, 1.0)
        } else {
            0.5
        };
        let color = [t as f32, 0.0, (1.0 - t) as f32, 1.0];
        Self {
            origin,
            tip,
            color,
            style,
            magnitude: mag,
        }
    }
    /// Build a set of glyphs from a `FluidField`.
    pub fn build_from_field(
        field: &FluidField,
        scale: f64,
        style: GlyphStyle,
    ) -> Vec<VelocityGlyph> {
        let max_mag = field
            .vector
            .iter()
            .map(|v| len3(*v))
            .fold(0.0_f64, f64::max);
        let ref_mag = if max_mag < 1e-12 { 1.0 } else { max_mag };
        let [nx, ny, nz] = field.dims;
        let mut glyphs = Vec::with_capacity(nx * ny * nz);
        for iz in 0..nz {
            for iy in 0..ny {
                for ix in 0..nx {
                    let idx = field.idx(ix, iy, iz);
                    let pos = field.position(ix, iy, iz);
                    let vel = field.vector[idx];
                    glyphs.push(VelocityGlyph::from_velocity(
                        pos, vel, scale, ref_mag, style,
                    ));
                }
            }
        }
        glyphs
    }
}
/// Fluid simulation statistics panel.
#[allow(dead_code)]
pub struct FluidStats {
    /// Reynolds number.
    pub reynolds: f64,
    /// Mach number.
    pub mach: f64,
    /// CFL number.
    pub cfl: f64,
    /// Maximum velocity magnitude.
    pub max_velocity: f64,
    /// Total fluid mass.
    pub total_mass: f64,
    /// Time step.
    pub dt: f64,
}
impl FluidStats {
    /// Create fluid stats.
    pub fn new() -> Self {
        Self {
            reynolds: 0.0,
            mach: 0.0,
            cfl: 0.0,
            max_velocity: 0.0,
            total_mass: 0.0,
            dt: 0.0,
        }
    }
    /// Compute Reynolds number: Re = ρ U L / μ.
    pub fn compute_reynolds(&mut self, density: f64, velocity: f64, length: f64, viscosity: f64) {
        self.reynolds = density * velocity * length / viscosity.max(1e-12);
    }
    /// Compute Mach number: Ma = U / c_sound.
    pub fn compute_mach(&mut self, velocity: f64, sound_speed: f64) {
        self.mach = velocity / sound_speed.max(1e-12);
    }
    /// Compute CFL number: CFL = U Δt / Δx.
    pub fn compute_cfl(&mut self, velocity: f64, dt: f64, dx: f64) {
        self.cfl = velocity * dt / dx.max(1e-12);
        self.dt = dt;
    }
    /// Format the stats as a multi-line string.
    pub fn format(&self) -> String {
        format!(
            "Re={:.1} Ma={:.3} CFL={:.3} Umax={:.3} m={:.3} dt={:.6}",
            self.reynolds, self.mach, self.cfl, self.max_velocity, self.total_mass, self.dt
        )
    }
}
/// A single computed streamline.
#[allow(dead_code)]
pub struct Streamline3D {
    /// Sequence of points along the streamline.
    pub points: Vec<[f64; 3]>,
    /// Velocity magnitudes at each point.
    pub magnitudes: Vec<f64>,
}
/// Velocity arrow field renderer.
///
/// Converts a vector field sample into renderable arrow glyphs.
#[allow(dead_code)]
pub struct VelocityArrows {
    /// Scale factor: length = magnitude * scale.
    pub scale: f64,
    /// Minimum arrow magnitude to display.
    pub min_magnitude: f64,
    /// Arrow colour at maximum magnitude.
    pub color_max: [f32; 4],
    /// Arrow colour at zero magnitude.
    pub color_min: [f32; 4],
    /// Reference magnitude for colour mapping.
    pub ref_magnitude: f64,
}
impl VelocityArrows {
    /// Create a velocity arrow renderer.
    pub fn new(scale: f64, ref_magnitude: f64) -> Self {
        Self {
            scale,
            min_magnitude: 0.0,
            color_max: [1.0, 0.0, 0.0, 1.0],
            color_min: [0.0, 0.0, 1.0, 1.0],
            ref_magnitude,
        }
    }
    /// Build an arrow from a velocity sample.
    pub fn make_arrow(&self, origin: [f64; 3], velocity: [f64; 3]) -> Option<VelocityArrow> {
        let mag = len3(velocity);
        if mag < self.min_magnitude {
            return None;
        }
        let tip = add3(origin, scale3(velocity, self.scale));
        let t = (mag / self.ref_magnitude).clamp(0.0, 1.0) as f32;
        let color = lerp_color(self.color_min, self.color_max, t);
        Some(VelocityArrow { origin, tip, color })
    }
    /// Build arrows for a grid of velocity samples.
    pub fn build_field(&self, origins: &[[f64; 3]], velocities: &[[f64; 3]]) -> Vec<VelocityArrow> {
        origins
            .iter()
            .zip(velocities.iter())
            .filter_map(|(&o, &v)| self.make_arrow(o, v))
            .collect()
    }
}
/// Line Integral Convolution renderer for 2D velocity field visualization.
#[allow(dead_code)]
pub struct LicRenderer {
    /// Image width in pixels.
    pub width: usize,
    /// Image height in pixels.
    pub height: usize,
    /// Integration length (number of steps each direction).
    pub steps: usize,
    /// Step size in normalized image coordinates.
    pub step_size: f64,
}
impl LicRenderer {
    /// Create a new LIC renderer.
    pub fn new(width: usize, height: usize, steps: usize, step_size: f64) -> Self {
        Self {
            width,
            height,
            steps,
            step_size,
        }
    }
    /// Compute the LIC image from velocity field (vx, vy) and noise texture.
    ///
    /// Returns a flat image of `width * height` intensity values.
    pub fn compute(&self, vx: &[f64], vy: &[f64], noise: &[f64]) -> Vec<f64> {
        let w = self.width;
        let h = self.height;
        let mut out = vec![0.0f64; w * h];
        for py in 0..h {
            for px in 0..w {
                let mut sum = 0.0;
                let mut count = 0.0;
                let mut cx = px as f64 / w as f64;
                let mut cy = py as f64 / h as f64;
                for step in 0..(self.steps * 2 + 1) {
                    let sign = if step <= self.steps { 1.0 } else { -1.0 };
                    let ix = (cx * w as f64).clamp(0.0, (w - 1) as f64) as usize;
                    let iy = (cy * h as f64).clamp(0.0, (h - 1) as f64) as usize;
                    let idx = iy * w + ix;
                    if idx < noise.len() {
                        sum += noise[idx];
                        count += 1.0;
                    }
                    if idx < vx.len() {
                        let v_len = (vx[idx] * vx[idx] + vy[idx] * vy[idx]).sqrt().max(1e-9);
                        cx += sign * self.step_size * vx[idx] / v_len;
                        cy += sign * self.step_size * vy[idx] / v_len;
                    }
                    cx = cx.clamp(0.0, 1.0);
                    cy = cy.clamp(0.0, 1.0);
                }
                out[py * w + px] = if count > 0.0 { sum / count } else { 0.0 };
            }
        }
        out
    }
}
/// Water surface wave visualization.
#[allow(dead_code)]
pub struct WaveViz {
    /// Normal map resolution.
    pub normal_map_res: usize,
    /// Fresnel base reflectance (F0).
    pub fresnel_f0: f64,
    /// Foam threshold (wave height fraction).
    pub foam_threshold: f64,
    /// Computed normal map (nx, ny, nz per pixel).
    pub normal_map: Vec<[f32; 3]>,
    /// Foam mask (0 = no foam, 1 = full foam).
    pub foam_mask: Vec<f32>,
}
impl WaveViz {
    /// Create a wave visualization.
    pub fn new(normal_map_res: usize) -> Self {
        let n2 = normal_map_res * normal_map_res;
        Self {
            normal_map_res,
            fresnel_f0: 0.02,
            foam_threshold: 0.7,
            normal_map: vec![[0.0, 1.0, 0.0]; n2],
            foam_mask: vec![0.0; n2],
        }
    }
    /// Compute Schlick Fresnel reflectance.
    pub fn fresnel(&self, cos_theta: f64) -> f64 {
        self.fresnel_f0 + (1.0 - self.fresnel_f0) * (1.0 - cos_theta).powi(5)
    }
    /// Update the normal map from a height field.
    pub fn update_normals(&mut self, heights: &[f64], cell_size: f64) {
        let n = self.normal_map_res;
        for i in 0..n {
            for j in 0..n {
                let hi_plus = heights[((i + 1).min(n - 1)) * n + j];
                let hi_minus = heights[i.saturating_sub(1) * n + j];
                let hj_plus = heights[i * n + (j + 1).min(n - 1)];
                let hj_minus = heights[i * n + j.saturating_sub(1)];
                let dx = (hi_plus - hi_minus) / (2.0 * cell_size);
                let dz = (hj_plus - hj_minus) / (2.0 * cell_size);
                let normal = norm3([-dx, 1.0, -dz]);
                self.normal_map[i * n + j] = [normal[0] as f32, normal[1] as f32, normal[2] as f32];
            }
        }
    }
    /// Update the foam mask based on wave height relative to the mean.
    pub fn update_foam(&mut self, heights: &[f64]) {
        let mean: f64 = heights.iter().sum::<f64>() / heights.len().max(1) as f64;
        let max_h = heights.iter().cloned().fold(f64::MIN, f64::max);
        let range = (max_h - mean).max(1e-9);
        for (i, &h) in heights.iter().enumerate() {
            let t = ((h - mean) / range).clamp(0.0, 1.0);
            self.foam_mask[i] = if t > self.foam_threshold {
                (t - self.foam_threshold) as f32 / (1.0 - self.foam_threshold as f32)
            } else {
                0.0
            };
        }
    }
}
/// Style for velocity field glyphs.
#[allow(dead_code)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GlyphStyle {
    /// Simple arrow glyphs.
    Arrow,
    /// Hedgehog-style lines (no arrowhead).
    Hedgehog,
    /// Line Integral Convolution texture-style glyphs.
    LIC,
}
/// A structured 2-D or 3-D scalar/vector field on a regular grid.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct FluidField {
    /// Number of cells along x, y, z.
    pub dims: [usize; 3],
    /// Physical cell spacing along each axis.
    pub spacing: [f64; 3],
    /// Flat scalar values stored in x-major order (size = dims\[0\]*dims\[1\]*dims\[2\]).
    pub scalar: Vec<f64>,
    /// Flat vector values stored as `[vx, vy, vz]` per cell.
    pub vector: Vec<[f64; 3]>,
}
impl FluidField {
    /// Create a zero-initialised field with the given dimensions and spacing.
    pub fn new(dims: [usize; 3], spacing: [f64; 3]) -> Self {
        let n = dims[0] * dims[1] * dims[2];
        Self {
            dims,
            spacing,
            scalar: vec![0.0; n],
            vector: vec![[0.0; 3]; n],
        }
    }
    /// Linear index for cell `(ix, iy, iz)`.
    pub fn idx(&self, ix: usize, iy: usize, iz: usize) -> usize {
        iz * self.dims[0] * self.dims[1] + iy * self.dims[0] + ix
    }
    /// Physical position of cell centre `(ix, iy, iz)`.
    pub fn position(&self, ix: usize, iy: usize, iz: usize) -> [f64; 3] {
        [
            (ix as f64 + 0.5) * self.spacing[0],
            (iy as f64 + 0.5) * self.spacing[1],
            (iz as f64 + 0.5) * self.spacing[2],
        ]
    }
    /// Total number of cells.
    pub fn num_cells(&self) -> usize {
        self.dims[0] * self.dims[1] * self.dims[2]
    }
    /// Scalar range `(min, max)` across all cells.
    pub fn scalar_range(&self) -> (f64, f64) {
        let mn = self.scalar.iter().cloned().fold(f64::INFINITY, f64::min);
        let mx = self
            .scalar
            .iter()
            .cloned()
            .fold(f64::NEG_INFINITY, f64::max);
        (mn, mx)
    }
}
/// A single splat (Gaussian footprint) for a fluid particle.
#[allow(dead_code)]
pub struct Splat {
    /// Screen-space position.
    pub position: [f32; 2],
    /// Splat radius (pixels).
    pub radius: f32,
    /// Colour and opacity.
    pub color: [f32; 4],
    /// Depth value for depth-sorted blending.
    pub depth: f32,
}
/// A 2D slice extracted from a 3D scalar field.
#[allow(dead_code)]
pub struct VolumeSlice {
    /// Grid resolution of the slice (pixels per axis).
    pub resolution: usize,
    /// Axis-aligned plane: 0 = YZ, 1 = XZ, 2 = XY.
    pub axis: usize,
    /// Slice position along the axis (world units).
    pub position: f64,
    /// Extracted pixel values.
    pub pixels: Vec<f64>,
    /// Pixel size (world units).
    pub pixel_size: f64,
}
impl VolumeSlice {
    /// Create a volume slice extractor.
    pub fn new(resolution: usize, axis: usize, position: f64, pixel_size: f64) -> Self {
        Self {
            resolution,
            axis,
            position,
            pixels: vec![0.0; resolution * resolution],
            pixel_size,
        }
    }
    /// Extract the slice from a 3D scalar field using bilinear interpolation.
    pub fn extract(&mut self, field: &[f64], field_res: usize, field_size: f64) {
        let h = field_size / field_res as f64;
        let slice_idx = (self.position / h).clamp(0.0, (field_res - 1) as f64) as usize;
        let n = self.resolution;
        for i in 0..n {
            for j in 0..n {
                let (fx, fy, fz) = match self.axis {
                    0 => (slice_idx, i, j),
                    1 => (i, slice_idx, j),
                    _ => (i, j, slice_idx),
                };
                let fx = fx.min(field_res - 1);
                let fy = fy.min(field_res - 1);
                let fz = fz.min(field_res - 1);
                self.pixels[i * n + j] = field[fx * field_res * field_res + fy * field_res + fz];
            }
        }
    }
    /// Bilinear interpolation within the extracted slice.
    pub fn sample_bilinear(&self, u: f64, v: f64) -> f64 {
        let n = self.resolution;
        let fx = (u * (n - 1) as f64).clamp(0.0, (n - 1) as f64);
        let fy = (v * (n - 1) as f64).clamp(0.0, (n - 1) as f64);
        let ix = fx as usize;
        let iy = fy as usize;
        let tx = fx - ix as f64;
        let ty = fy - iy as f64;
        let ix1 = (ix + 1).min(n - 1);
        let iy1 = (iy + 1).min(n - 1);
        let v00 = self.pixels[ix * n + iy];
        let v10 = self.pixels[ix1 * n + iy];
        let v01 = self.pixels[ix * n + iy1];
        let v11 = self.pixels[ix1 * n + iy1];
        (1.0 - tx) * (1.0 - ty) * v00
            + tx * (1.0 - ty) * v10
            + (1.0 - tx) * ty * v01
            + tx * ty * v11
    }
}
/// A single velocity arrow glyph.
#[allow(dead_code)]
pub struct VelocityArrow {
    /// Arrow base position.
    pub origin: [f64; 3],
    /// Arrow tip position.
    pub tip: [f64; 3],
    /// Arrow colour (RGBA, 0..1).
    pub color: [f32; 4],
}
/// LBM lattice visualisation: density field, velocity magnitude, streamlines.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct LbmViz {
    /// Grid dimensions `[nx, ny, nz]`.
    pub dims: [usize; 3],
    /// Density field (one value per cell).
    pub density: Vec<f64>,
    /// Velocity field (three components per cell).
    pub velocity: Vec<[f64; 3]>,
}
impl LbmViz {
    /// Create a zero-initialised LBM visualisation grid.
    pub fn new(dims: [usize; 3]) -> Self {
        let n = dims[0] * dims[1] * dims[2];
        Self {
            dims,
            density: vec![1.0; n],
            velocity: vec![[0.0; 3]; n],
        }
    }
    /// Linear index for cell `(ix, iy, iz)`.
    pub fn idx(&self, ix: usize, iy: usize, iz: usize) -> usize {
        iz * self.dims[0] * self.dims[1] + iy * self.dims[0] + ix
    }
    /// Velocity magnitude field.
    pub fn velocity_magnitude(&self) -> Vec<f64> {
        self.velocity.iter().map(|v| len3(*v)).collect()
    }
    /// Density range `(min, max)`.
    pub fn density_range(&self) -> (f64, f64) {
        let mn = self.density.iter().cloned().fold(f64::INFINITY, f64::min);
        let mx = self
            .density
            .iter()
            .cloned()
            .fold(f64::NEG_INFINITY, f64::max);
        (mn, mx)
    }
    /// Maximum velocity magnitude across all cells.
    pub fn max_velocity_magnitude(&self) -> f64 {
        self.velocity
            .iter()
            .map(|v| len3(*v))
            .fold(0.0_f64, f64::max)
    }
    /// Generate velocity arrow glyphs for a 2-D XY slice at `iz`.
    pub fn slice_glyphs(&self, iz: usize, scale: f64) -> Vec<VelocityGlyph> {
        let [nx, ny, _] = self.dims;
        let mut glyphs = Vec::new();
        let max_mag = self.max_velocity_magnitude().max(1e-12);
        for iy in 0..ny {
            for ix in 0..nx {
                let idx = self.idx(ix, iy, iz);
                let origin = [ix as f64, iy as f64, iz as f64];
                glyphs.push(VelocityGlyph::from_velocity(
                    origin,
                    self.velocity[idx],
                    scale,
                    max_mag,
                    GlyphStyle::Arrow,
                ));
            }
        }
        glyphs
    }
}
/// Export a sequence of fluid frames for animation.
#[allow(dead_code)]
#[derive(Debug, Default)]
pub struct AnimationExport {
    /// Accumulated frames (each frame is a flat list of glyph origins).
    pub frames: Vec<Vec<[f64; 3]>>,
    /// Simulation time for each frame.
    pub times: Vec<f64>,
    /// Target frames-per-second for playback.
    pub fps: f64,
}
impl AnimationExport {
    /// Create an empty `AnimationExport`.
    pub fn new(fps: f64) -> Self {
        Self {
            frames: Vec::new(),
            times: Vec::new(),
            fps,
        }
    }
    /// Push a new frame (list of particle/glyph positions) at simulation time `t`.
    pub fn push_frame(&mut self, positions: Vec<[f64; 3]>, t: f64) {
        self.frames.push(positions);
        self.times.push(t);
    }
    /// Number of frames accumulated.
    pub fn num_frames(&self) -> usize {
        self.frames.len()
    }
    /// Duration of the animation in seconds.
    pub fn duration(&self) -> f64 {
        if self.times.is_empty() {
            0.0
        } else {
            self.times[self.times.len() - 1] - self.times[0]
        }
    }
    /// Playback time in seconds for frame `i` based on fps.
    pub fn playback_time(&self, frame_index: usize) -> f64 {
        frame_index as f64 / self.fps.max(1e-10)
    }
    /// Build a simple CSV export string (frame index, sim time, particle count).
    pub fn to_csv(&self) -> String {
        let mut out = String::from("frame,sim_time,particle_count\n");
        for (i, (frame, &t)) in self.frames.iter().zip(self.times.iter()).enumerate() {
            out.push_str(&format!("{},{:.6},{}\n", i, t, frame.len()));
        }
        out
    }
}
/// A triangle produced by marching cubes.
#[allow(dead_code)]
pub struct IsoTriangle {
    /// Vertex positions.
    pub vertices: [[f64; 3]; 3],
    /// Per-vertex normals.
    pub normals: [[f64; 3]; 3],
    /// Iso-level value at which this triangle was extracted.
    pub level: f64,
}

impl IsoTriangle {
    /// Compute the face normal from the triangle vertices using the cross product.
    /// Returns a unit-length normal vector.
    pub fn face_normal(&self) -> [f64; 3] {
        let [a, b, c] = self.vertices;
        let ab = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
        let ac = [c[0] - a[0], c[1] - a[1], c[2] - a[2]];
        let nx = ab[1] * ac[2] - ab[2] * ac[1];
        let ny = ab[2] * ac[0] - ab[0] * ac[2];
        let nz = ab[0] * ac[1] - ab[1] * ac[0];
        let len = (nx * nx + ny * ny + nz * nz).sqrt();
        if len < 1e-15 {
            [0.0, 0.0, 0.0]
        } else {
            [nx / len, ny / len, nz / len]
        }
    }

    /// Compute the centroid (average of the three vertices).
    pub fn centroid(&self) -> [f64; 3] {
        let [a, b, c] = self.vertices;
        [
            (a[0] + b[0] + c[0]) / 3.0,
            (a[1] + b[1] + c[1]) / 3.0,
            (a[2] + b[2] + c[2]) / 3.0,
        ]
    }
}
/// Streamline renderer using Runge-Kutta integration.
#[allow(dead_code)]
pub struct StreamlineRenderer {
    /// Maximum number of steps per streamline.
    pub max_steps: usize,
    /// Integration step size (world units).
    pub step_size: f64,
    /// Minimum velocity magnitude to continue integration.
    pub min_velocity: f64,
    /// Tube radius for mesh generation.
    pub tube_radius: f64,
}
impl StreamlineRenderer {
    /// Create a streamline renderer.
    pub fn new(max_steps: usize, step_size: f64) -> Self {
        Self {
            max_steps,
            step_size,
            min_velocity: 1e-6,
            tube_radius: 0.01,
        }
    }
    /// Trace a streamline from `seed` using RK4 integration.
    ///
    /// `velocity_fn` maps a 3D position to a velocity vector.
    pub fn trace<F>(&self, seed: [f64; 3], velocity_fn: F) -> Streamline3D
    where
        F: Fn([f64; 3]) -> [f64; 3],
    {
        let mut points = vec![seed];
        let mut magnitudes = Vec::new();
        let mut pos = seed;
        for _ in 0..self.max_steps {
            let v = velocity_fn(pos);
            let mag = len3(v);
            magnitudes.push(mag);
            if mag < self.min_velocity {
                break;
            }
            let k1 = v;
            let k2 = velocity_fn(add3(pos, scale3(k1, self.step_size * 0.5)));
            let k3 = velocity_fn(add3(pos, scale3(k2, self.step_size * 0.5)));
            let k4 = velocity_fn(add3(pos, scale3(k3, self.step_size)));
            let step = scale3(
                add3(add3(k1, scale3(k2, 2.0)), add3(scale3(k3, 2.0), k4)),
                self.step_size / 6.0,
            );
            pos = add3(pos, step);
            points.push(pos);
        }
        Streamline3D { points, magnitudes }
    }
    /// Trace multiple streamlines from a seed array.
    pub fn trace_many<F>(&self, seeds: &[[f64; 3]], velocity_fn: F) -> Vec<Streamline3D>
    where
        F: Fn([f64; 3]) -> [f64; 3],
    {
        seeds.iter().map(|&s| self.trace(s, &velocity_fn)).collect()
    }
}
/// Isosurface visualizer using marching cubes.
#[allow(dead_code)]
pub struct IsosurfaceViz {
    /// Isosurface threshold value.
    pub threshold: f64,
    /// Grid resolution (cells per axis).
    pub resolution: usize,
    /// Cell size (world units).
    pub cell_size: f64,
    /// Whether to compute ambient occlusion.
    pub ambient_occlusion: bool,
}
impl IsosurfaceViz {
    /// Create an isosurface visualizer.
    pub fn new(threshold: f64, resolution: usize, domain_size: f64) -> Self {
        Self {
            threshold,
            resolution,
            cell_size: domain_size / resolution as f64,
            ambient_occlusion: false,
        }
    }
    /// Interpolate the edge crossing between two scalar values.
    pub fn edge_interp(&self, p0: [f64; 3], v0: f64, p1: [f64; 3], v1: f64) -> [f64; 3] {
        let dv = v1 - v0;
        if dv.abs() < 1e-12 {
            return p0;
        }
        let t = (self.threshold - v0) / dv;
        add3(p0, scale3(sub3(p1, p0), t))
    }
    /// Compute a vertex normal by finite differences of the scalar field.
    pub fn compute_normal(&self, scalars: &[f64], ix: usize, iy: usize, iz: usize) -> [f64; 3] {
        let n = self.resolution;
        let idx = |x: usize, y: usize, z: usize| x * n * n + y * n + z;
        let h = self.cell_size;
        let dx = if ix + 1 < n && ix > 0 {
            scalars[idx(ix + 1, iy, iz)] - scalars[idx(ix - 1, iy, iz)]
        } else {
            0.0
        };
        let dy = if iy + 1 < n && iy > 0 {
            scalars[idx(ix, iy + 1, iz)] - scalars[idx(ix, iy - 1, iz)]
        } else {
            0.0
        };
        let dz = if iz + 1 < n && iz > 0 {
            scalars[idx(ix, iy, iz + 1)] - scalars[idx(ix, iy, iz - 1)]
        } else {
            0.0
        };
        norm3([dx / (2.0 * h), dy / (2.0 * h), dz / (2.0 * h)])
    }
    /// Extract a simplified isosurface: returns all edge-interpolated crossing points.
    pub fn extract_crossings(&self, scalars: &[f64]) -> Vec<[f64; 3]> {
        let n = self.resolution;
        let mut crossings = Vec::new();
        let idx = |x: usize, y: usize, z: usize| x * n * n + y * n + z;
        for ix in 0..n.saturating_sub(1) {
            for iy in 0..n.saturating_sub(1) {
                for iz in 0..n.saturating_sub(1) {
                    let v0 = scalars[idx(ix, iy, iz)];
                    let v1 = scalars[idx(ix + 1, iy, iz)];
                    if (v0 < self.threshold) != (v1 < self.threshold) {
                        let p0 = [
                            ix as f64 * self.cell_size,
                            iy as f64 * self.cell_size,
                            iz as f64 * self.cell_size,
                        ];
                        let p1 = [
                            (ix + 1) as f64 * self.cell_size,
                            iy as f64 * self.cell_size,
                            iz as f64 * self.cell_size,
                        ];
                        crossings.push(self.edge_interp(p0, v0, p1, v1));
                    }
                }
            }
        }
        crossings
    }
}
/// A single metaball in a metaball field.
#[allow(dead_code)]
pub struct Metaball {
    /// World-space center.
    pub center: [f64; 3],
    /// Effective radius.
    pub radius: f64,
    /// Strength (scalar multiplier).
    pub strength: f64,
}
/// Vorticity field visualizer.
#[allow(dead_code)]
pub struct VorticityViz {
    /// Minimum vorticity magnitude to display.
    pub min_vorticity: f64,
    /// Lambda2 threshold for vortex core extraction.
    pub lambda2_threshold: f64,
    /// Colormap range maximum.
    pub max_vorticity: f64,
}
impl VorticityViz {
    /// Create a vorticity visualizer.
    pub fn new(max_vorticity: f64) -> Self {
        Self {
            min_vorticity: 0.0,
            lambda2_threshold: -0.1,
            max_vorticity,
        }
    }
    /// Map vorticity magnitude to a colour (red = high, blue = low).
    pub fn color(&self, magnitude: f64) -> [f32; 4] {
        let t = (magnitude / self.max_vorticity.max(1e-12)).clamp(0.0, 1.0) as f32;
        [t, 0.0, 1.0 - t, 1.0]
    }
    /// Compute vorticity from a velocity gradient (Jacobian matrix, row-major 3×3).
    ///
    /// Returns the vorticity vector `(∂w/∂y - ∂v/∂z, ∂u/∂z - ∂w/∂x, ∂v/∂x - ∂u/∂y)`.
    pub fn vorticity_from_jacobian(&self, j: [f64; 9]) -> [f64; 3] {
        [j[7] - j[5], j[2] - j[6], j[3] - j[1]]
    }
    /// Extract vortex core lines by testing lambda2 criterion.
    ///
    /// Returns indices of cells where lambda2 < threshold.
    pub fn vortex_core_cells(&self, lambda2_values: &[f64]) -> Vec<usize> {
        lambda2_values
            .iter()
            .enumerate()
            .filter(|&(_, &v)| v < self.lambda2_threshold)
            .map(|(i, _)| i)
            .collect()
    }
}
/// Bubble renderer for particle-based fluids.
#[allow(dead_code)]
pub struct BubbleRenderer {
    /// Bubbles to render.
    pub bubbles: Vec<Bubble>,
    /// Maximum number of bubbles.
    pub max_bubbles: usize,
    /// Water IOR.
    pub water_ior: f64,
}
impl BubbleRenderer {
    /// Create a bubble renderer.
    pub fn new(max_bubbles: usize) -> Self {
        Self {
            bubbles: Vec::new(),
            max_bubbles,
            water_ior: 1.333,
        }
    }
    /// Add a bubble, respecting the maximum count.
    pub fn add_bubble(&mut self, b: Bubble) {
        if self.bubbles.len() < self.max_bubbles {
            self.bubbles.push(b);
        }
    }
    /// Remove bubbles that have risen above `surface_y`.
    pub fn prune_above_surface(&mut self, surface_y: f64) {
        self.bubbles.retain(|b| b.center[1] < surface_y);
    }
    /// Count visible bubbles (opacity > threshold).
    pub fn visible_count(&self, min_opacity: f32) -> usize {
        self.bubbles
            .iter()
            .filter(|b| b.opacity >= min_opacity)
            .count()
    }
}
/// A metaball implicit surface field.
///
/// Used for particle-based fluid surface reconstruction.
#[allow(dead_code)]
pub struct MetaballField {
    /// List of metaballs.
    pub balls: Vec<Metaball>,
    /// Iso-threshold for surface extraction.
    pub threshold: f64,
}
impl MetaballField {
    /// Create a new metaball field.
    pub fn new(balls: Vec<Metaball>, threshold: f64) -> Self {
        Self { balls, threshold }
    }
    /// Evaluate the metaball scalar field at `point`.
    ///
    /// Field = sum_i strength_i / (||p - c_i||² + epsilon)
    pub fn evaluate(&self, point: [f64; 3]) -> f64 {
        self.balls.iter().fold(0.0, |acc, b| {
            let d2 = {
                let dx = point[0] - b.center[0];
                let dy = point[1] - b.center[1];
                let dz = point[2] - b.center[2];
                dx * dx + dy * dy + dz * dz
            };
            acc + b.strength / (d2 + 1e-6)
        })
    }
    /// Extract iso-surface crossings along the X axis for a 1D scan.
    pub fn scan_x(&self, y: f64, z: f64, x_min: f64, x_max: f64, steps: usize) -> Vec<f64> {
        let mut crossings = Vec::new();
        let dx = (x_max - x_min) / steps as f64;
        let mut prev = self.evaluate([x_min, y, z]);
        for i in 1..=steps {
            let x = x_min + i as f64 * dx;
            let v = self.evaluate([x, y, z]);
            if (prev < self.threshold) != (v < self.threshold) {
                let t = (self.threshold - prev) / (v - prev + 1e-15);
                crossings.push(x_min + (i as f64 - 1.0 + t) * dx);
            }
            prev = v;
        }
        crossings
    }
}
/// A pathline integrator that tracks particle trajectories over time.
#[allow(dead_code)]
pub struct PathlineIntegrator {
    /// Maximum number of history steps.
    pub max_steps: usize,
    /// Integration step size.
    pub dt: f64,
    /// Particle paths (one `Vec<[f64;3]>` per seed).
    pub paths: Vec<Vec<[f64; 3]>>,
}
impl PathlineIntegrator {
    /// Create a new pathline integrator.
    pub fn new(max_steps: usize, dt: f64) -> Self {
        Self {
            max_steps,
            dt,
            paths: Vec::new(),
        }
    }
    /// Add a seed particle.
    pub fn add_seed(&mut self, pos: [f64; 3]) {
        self.paths.push(vec![pos]);
    }
    /// Advance all paths by one step using the given velocity field.
    pub fn step<F: Fn([f64; 3]) -> [f64; 3]>(&mut self, velocity_fn: F) {
        let dt = self.dt;
        for path in &mut self.paths {
            if path.len() >= self.max_steps {
                continue;
            }
            if let Some(&last) = path.last() {
                let v = velocity_fn(last);
                let next = [
                    last[0] + v[0] * dt,
                    last[1] + v[1] * dt,
                    last[2] + v[2] * dt,
                ];
                path.push(next);
            }
        }
    }
    /// Get the total length of all paths.
    pub fn total_points(&self) -> usize {
        self.paths.iter().map(|p| p.len()).sum()
    }
}
/// Particle splatting renderer (SPH/LBM particles).
#[allow(dead_code)]
pub struct ParticleSplat {
    /// Gaussian splat sigma (fraction of radius).
    pub sigma: f64,
    /// Minimum opacity for a particle to be rendered.
    pub min_opacity: f32,
    /// Blending mode name (stored for reference).
    pub blend_mode: String,
    /// Splat scale (pixels per world unit).
    pub scale: f64,
}
impl ParticleSplat {
    /// Create a particle splatting renderer.
    pub fn new(scale: f64) -> Self {
        Self {
            sigma: 0.3,
            min_opacity: 0.01,
            blend_mode: "additive".to_string(),
            scale,
        }
    }
    /// Compute the Gaussian weight at distance `d` from the splat centre.
    pub fn gaussian_weight(&self, d: f64, radius: f64) -> f64 {
        let s = radius * self.sigma;
        (-0.5 * (d / s).powi(2)).exp()
    }
    /// Build a splat for a particle at world position `pos` projected to screen `screen_pos`.
    pub fn build_splat(
        &self,
        screen_pos: [f32; 2],
        world_radius: f64,
        depth: f32,
        color: [f32; 4],
    ) -> Option<Splat> {
        let pixel_radius = (world_radius * self.scale) as f32;
        if color[3] < self.min_opacity {
            return None;
        }
        Some(Splat {
            position: screen_pos,
            radius: pixel_radius,
            color,
            depth,
        })
    }
    /// Sort splats back-to-front for correct alpha blending.
    pub fn depth_sort(splats: &mut Vec<Splat>) {
        splats.sort_by(|a, b| {
            b.depth
                .partial_cmp(&a.depth)
                .unwrap_or(std::cmp::Ordering::Equal)
        });
    }
}
/// Visualise a free surface from a VOF or level-set scalar field.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct FreeSurfaceViz {
    /// Level-set or VOF field (positive = fluid, negative = air).
    pub field: FluidField,
    /// Iso-value at which the surface is extracted (typically 0.0 or 0.5).
    pub iso_value: f64,
}
impl FreeSurfaceViz {
    /// Create a `FreeSurfaceViz` from a `FluidField`.
    pub fn new(field: FluidField, iso_value: f64) -> Self {
        Self { field, iso_value }
    }
    /// Extract crossing points (simplified marching-squares/cubes edge points).
    ///
    /// Returns pairs of positions where the iso-surface crosses grid edges.
    pub fn extract_crossings(&self) -> Vec<[f64; 3]> {
        let [nx, ny, nz] = self.field.dims;
        let mut crossings = Vec::new();
        for iz in 0..nz {
            for iy in 0..ny {
                for ix in 0..nx.saturating_sub(1) {
                    let i0 = self.field.idx(ix, iy, iz);
                    let i1 = self.field.idx(ix + 1, iy, iz);
                    let v0 = self.field.scalar[i0] - self.iso_value;
                    let v1 = self.field.scalar[i1] - self.iso_value;
                    if v0 * v1 < 0.0 {
                        let t = v0 / (v0 - v1);
                        let p0 = self.field.position(ix, iy, iz);
                        let p1 = self.field.position(ix + 1, iy, iz);
                        crossings.push([
                            p0[0] + t * (p1[0] - p0[0]),
                            p0[1] + t * (p1[1] - p0[1]),
                            p0[2] + t * (p1[2] - p0[2]),
                        ]);
                    }
                }
            }
        }
        crossings
    }
    /// Return the fraction of cells with scalar > iso_value (fluid fraction).
    pub fn fluid_fraction(&self) -> f64 {
        let total = self.field.scalar.len();
        if total == 0 {
            return 0.0;
        }
        let fluid = self
            .field
            .scalar
            .iter()
            .filter(|&&v| v > self.iso_value)
            .count();
        fluid as f64 / total as f64
    }
}
