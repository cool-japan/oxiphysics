//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::*;
/// Configuration for point-sprite particle rendering.
#[derive(Debug, Clone)]
pub struct PointSpriteConfig {
    /// Particle radius in world-space units.
    pub size_world: f32,
    /// Whether to stretch the sprite along the velocity direction.
    pub use_velocity_stretch: bool,
    /// Stretch multiplier applied to the velocity magnitude.
    pub stretch_factor: f32,
}
/// A single trail point.
#[derive(Debug, Clone)]
pub struct TrailPoint {
    /// Position.
    pub position: [f32; 3],
    /// Age of this point (seconds since creation).
    pub age: f32,
    /// Width at this point.
    pub width: f32,
    /// Color at this point.
    pub color: [f32; 4],
}
/// A velocity glyph arrow with base, tip, and color.
#[derive(Debug, Clone)]
pub struct VelocityGlyph {
    /// Arrow base (particle position).
    pub base: [f32; 3],
    /// Arrow tip (base + velocity * scale).
    pub tip: [f32; 3],
    /// RGBA color (typically mapped from speed).
    pub color: [f32; 4],
    /// Velocity magnitude.
    pub magnitude: f32,
}
/// A particle trail that stores a history of positions.
#[derive(Debug, Clone)]
pub struct ParticleTrail {
    /// Trail points (most recent last).
    pub points: Vec<TrailPoint>,
    /// Maximum number of points to keep.
    pub max_points: usize,
    /// Maximum age before a point is removed (seconds).
    pub max_age: f32,
}
impl ParticleTrail {
    /// Create a new particle trail.
    pub fn new(max_points: usize, max_age: f32) -> Self {
        Self {
            points: Vec::with_capacity(max_points),
            max_points,
            max_age,
        }
    }
    /// Add a new point to the trail.
    pub fn add_point(&mut self, position: [f32; 3], width: f32, color: [f32; 4]) {
        if self.points.len() >= self.max_points {
            self.points.remove(0);
        }
        self.points.push(TrailPoint {
            position,
            age: 0.0,
            width,
            color,
        });
    }
    /// Age all points by `dt` seconds and remove those that exceed `max_age`.
    pub fn update(&mut self, dt: f32) {
        for p in &mut self.points {
            p.age += dt;
        }
        self.points.retain(|p| p.age < self.max_age);
    }
    /// Generate line strip vertex data: `[px, py, pz, cr, cg, cb, ca]` per vertex.
    pub fn to_line_strip_data(&self) -> Vec<f32> {
        let mut out = Vec::with_capacity(self.points.len() * 7);
        for p in &self.points {
            out.extend_from_slice(&p.position);
            let alpha_factor = 1.0 - (p.age / self.max_age).clamp(0.0, 1.0);
            out.push(p.color[0]);
            out.push(p.color[1]);
            out.push(p.color[2]);
            out.push(p.color[3] * alpha_factor);
        }
        out
    }
    /// Number of active trail points.
    pub fn point_count(&self) -> usize {
        self.points.len()
    }
}
/// CPU-side buffer holding up to `max_particles` particle instances.
#[derive(Debug, Clone)]
pub struct ParticleBuffer {
    /// Stored particle instances.
    pub instances: Vec<ParticleInstance>,
    /// Maximum number of particles this buffer can hold.
    pub max_particles: usize,
}
impl ParticleBuffer {
    /// Create a new empty buffer with the given capacity.
    pub fn new(max_particles: usize) -> Self {
        Self {
            instances: Vec::with_capacity(max_particles),
            max_particles,
        }
    }
    /// Add a particle at `pos` with the given `radius` and `color`.
    ///
    /// Returns `false` if the buffer is already full (particle is not added).
    pub fn add(&mut self, pos: [f32; 3], radius: f32, color: [f32; 4]) -> bool {
        if self.instances.len() >= self.max_particles {
            return false;
        }
        self.instances.push(ParticleInstance {
            position: pos,
            radius,
            color,
            velocity: [0.0; 3],
        });
        true
    }
    /// Remove all particles from the buffer.
    pub fn clear(&mut self) {
        self.instances.clear();
    }
    /// Return the number of particles currently in the buffer.
    pub fn count(&self) -> usize {
        self.instances.len()
    }
    /// Produce interleaved vertex data: `[px, py, pz, radius, cr, cg, cb, ca]` per particle.
    pub fn to_vertex_data(&self) -> Vec<f32> {
        let mut out = Vec::with_capacity(self.instances.len() * 8);
        for p in &self.instances {
            out.push(p.position[0]);
            out.push(p.position[1]);
            out.push(p.position[2]);
            out.push(p.radius);
            out.push(p.color[0]);
            out.push(p.color[1]);
            out.push(p.color[2]);
            out.push(p.color[3]);
        }
        out
    }
    /// Build a buffer from a slice of positions, assigning the same `radius` and `color`
    /// to every particle.
    pub fn from_positions(positions: &[[f32; 3]], radius: f32, color: [f32; 4]) -> Self {
        let mut buf = Self::new(positions.len());
        for &pos in positions {
            buf.add(pos, radius, color);
        }
        buf
    }
}
/// A single velocity arrow glyph.
#[derive(Debug, Clone)]
pub struct VelocityArrow {
    /// Tail of the arrow.
    pub origin: [f32; 3],
    /// Unit direction of the arrow.
    pub direction: [f32; 3],
    /// Length of the arrow (original velocity magnitude × scale).
    pub magnitude: f32,
    /// RGBA color.
    pub color: [f32; 4],
}
/// A single particle instance ready to be uploaded to a GPU buffer.
#[derive(Debug, Clone, PartialEq)]
pub struct ParticleInstance {
    /// World-space position \[x, y, z\].
    pub position: [f32; 3],
    /// Sphere radius in world units.
    pub radius: f32,
    /// RGBA color, each component in \[0, 1\].
    pub color: [f32; 4],
    /// World-space velocity \[vx, vy, vz\].
    pub velocity: [f32; 3],
}
/// Buffer of velocity arrows.
#[derive(Debug, Clone)]
pub struct ArrowBuffer {
    /// Stored arrows.
    pub arrows: Vec<VelocityArrow>,
}
impl ArrowBuffer {
    /// Create an empty arrow buffer.
    pub fn new() -> Self {
        Self { arrows: Vec::new() }
    }
    /// Add an arrow at `pos` for velocity `vel`.
    ///
    /// The arrow length is `|vel| * scale`.  Color is taken from `colorizer` using
    /// the velocity magnitude.
    pub fn add(&mut self, pos: [f32; 3], vel: [f32; 3], scale: f32, colorizer: &ScalarColorizer) {
        let mag = (vel[0] * vel[0] + vel[1] * vel[1] + vel[2] * vel[2]).sqrt();
        let direction = if mag > f32::EPSILON {
            [vel[0] / mag, vel[1] / mag, vel[2] / mag]
        } else {
            [0.0, 0.0, 0.0]
        };
        let color = colorizer.colorize(mag);
        self.arrows.push(VelocityArrow {
            origin: pos,
            direction,
            magnitude: mag * scale,
            color,
        });
    }
    /// Produce line vertex data: each arrow becomes two points `[x0,y0,z0, x1,y1,z1]`.
    pub fn to_line_data(&self) -> Vec<f32> {
        let mut out = Vec::with_capacity(self.arrows.len() * 6);
        for a in &self.arrows {
            let tip = [
                a.origin[0] + a.direction[0] * a.magnitude,
                a.origin[1] + a.direction[1] * a.magnitude,
                a.origin[2] + a.direction[2] * a.magnitude,
            ];
            out.extend_from_slice(&a.origin);
            out.extend_from_slice(&tip);
        }
        out
    }
    /// Return the maximum arrow magnitude, or 0.0 if the buffer is empty.
    pub fn max_magnitude(&self) -> f32 {
        self.arrows
            .iter()
            .map(|a| a.magnitude)
            .fold(0.0_f32, f32::max)
    }
}
/// A single metaball (implicit sphere blob).
#[derive(Debug, Clone, Copy)]
pub struct Metaball {
    /// Centre position.
    pub center: [f32; 3],
    /// Strength coefficient.
    pub strength: f32,
    /// Falloff radius.
    pub radius: f32,
}
/// Maps scalar values from a given range to colors using a `ColorMap`.
#[derive(Debug, Clone)]
pub struct ScalarColorizer {
    /// Minimum value of the input range (maps to t = 0).
    pub min_val: f32,
    /// Maximum value of the input range (maps to t = 1).
    pub max_val: f32,
    /// Color map used for the mapping.
    pub color_map: ColorMap,
}
impl ScalarColorizer {
    /// Create a new colorizer for the interval `[min_val, max_val]`.
    pub fn new(min_val: f32, max_val: f32, color_map: ColorMap) -> Self {
        Self {
            min_val,
            max_val,
            color_map,
        }
    }
    /// Normalize `value` to \[0, 1\] relative to `[min_val, max_val]` and return the
    /// corresponding color.  If the range has zero width, clamps to 0.5.
    pub fn colorize(&self, value: f32) -> [f32; 4] {
        let range = self.max_val - self.min_val;
        let t = if range.abs() < f32::EPSILON {
            0.5
        } else {
            ((value - self.min_val) / range).clamp(0.0, 1.0)
        };
        self.color_map.map(t)
    }
    /// Colorize a slice of scalar values.
    pub fn colorize_array(&self, values: &[f32]) -> Vec<[f32; 4]> {
        values.iter().map(|&v| self.colorize(v)).collect()
    }
}
/// A billboarded particle quad that always faces the camera.
#[derive(Debug, Clone)]
pub struct BillboardParticle {
    /// World-space centre position.
    pub position: [f32; 3],
    /// Half-size of the billboard quad.
    pub half_size: f32,
    /// RGBA color.
    pub color: [f32; 4],
    /// Texture index (for sprite sheets).
    pub texture_id: u32,
}
/// Renderer that manages a list of point sprites and generates vertex data.
#[derive(Debug, Clone)]
pub struct PointSpriteRenderer {
    /// Configuration for all sprites.
    pub config: PointSpriteConfig,
    /// Sprite list.
    pub sprites: Vec<PointSprite>,
}
impl PointSpriteRenderer {
    /// Create a new empty renderer.
    pub fn new(config: PointSpriteConfig) -> Self {
        Self {
            config,
            sprites: Vec::new(),
        }
    }
    /// Add a sprite at `pos` with velocity `vel` and `color`.
    ///
    /// The size is computed from `config.size_world` plus optional velocity stretch.
    pub fn add(&mut self, pos: [f32; 3], vel: [f32; 3], color: [f32; 4]) {
        let speed = (vel[0] * vel[0] + vel[1] * vel[1] + vel[2] * vel[2]).sqrt();
        let size = if self.config.use_velocity_stretch {
            self.config.size_world + speed * self.config.stretch_factor
        } else {
            self.config.size_world
        };
        self.sprites.push(PointSprite {
            position: pos,
            velocity: vel,
            color,
            size,
        });
    }
    /// Number of sprites.
    pub fn count(&self) -> usize {
        self.sprites.len()
    }
    /// Clear all sprites.
    pub fn clear(&mut self) {
        self.sprites.clear();
    }
    /// Generate flat vertex data per sprite:
    /// `[px, py, pz, vx, vy, vz, cr, cg, cb, ca, size]` (11 floats/sprite).
    pub fn to_vertex_data(&self) -> Vec<f32> {
        let mut out = Vec::with_capacity(self.sprites.len() * 11);
        for s in &self.sprites {
            out.extend_from_slice(&s.position);
            out.extend_from_slice(&s.velocity);
            out.extend_from_slice(&s.color);
            out.push(s.size);
        }
        out
    }
    /// Sort sprites back-to-front for alpha blending.
    pub fn sort_back_to_front(&mut self, camera_pos: [f32; 3]) {
        self.sprites.sort_by(|a, b| {
            let da = {
                let d = [
                    a.position[0] - camera_pos[0],
                    a.position[1] - camera_pos[1],
                    a.position[2] - camera_pos[2],
                ];
                d[0] * d[0] + d[1] * d[1] + d[2] * d[2]
            };
            let db = {
                let d = [
                    b.position[0] - camera_pos[0],
                    b.position[1] - camera_pos[1],
                    b.position[2] - camera_pos[2],
                ];
                d[0] * d[0] + d[1] * d[1] + d[2] * d[2]
            };
            db.partial_cmp(&da).unwrap_or(std::cmp::Ordering::Equal)
        });
    }
}
/// Rendering statistics collected during a particle draw pass.
#[derive(Debug, Clone, Copy, Default)]
pub struct ParticleRenderStats {
    /// Number of particles that passed frustum culling.
    pub n_visible: usize,
    /// Number of particles culled.
    pub n_culled: usize,
    /// Number of draw calls issued.
    pub draw_calls: usize,
}
/// Gaussian-splatting renderer configuration.
#[derive(Debug, Clone, Copy)]
pub struct SplatRenderer {
    /// Standard deviation of the Gaussian splat in world units.
    pub sigma: f64,
    /// Radius beyond which the splat contribution is zero.
    pub cutoff_radius: f64,
}
impl SplatRenderer {
    /// Compute the Gaussian weight for a pixel at `center_uv` when splatting
    /// a particle at world position `pos` with radius `radius`.
    ///
    /// `pixel_size` is the world-space size of one pixel.
    /// Returns a value in `[0, 1]` representing the Gaussian contribution.
    pub fn splat_to_buffer(
        &self,
        _pos: [f64; 3],
        _radius: f64,
        center_uv: [f64; 2],
        pixel_size: f64,
    ) -> f64 {
        let dx = (center_uv[0] - 0.5) * pixel_size;
        let dy = (center_uv[1] - 0.5) * pixel_size;
        let r2 = dx * dx + dy * dy;
        let two_sigma2 = 2.0 * self.sigma * self.sigma;
        if two_sigma2 < 1e-30 {
            return 0.0;
        }
        (-r2 / two_sigma2).exp()
    }
}
/// A collection of metaballs that together form an implicit density field.
#[derive(Debug, Clone)]
pub struct MetaballField {
    /// All contributing blobs.
    pub blobs: Vec<Metaball>,
}
impl MetaballField {
    /// Create an empty metaball field.
    pub fn new() -> Self {
        Self { blobs: Vec::new() }
    }
    /// Add a blob at `center` with given `strength` and `radius`.
    pub fn add(&mut self, center: [f32; 3], strength: f32, radius: f32) {
        self.blobs.push(Metaball {
            center,
            strength,
            radius,
        });
    }
    /// Evaluate the combined density field at world-space point `p`.
    pub fn evaluate(&self, p: [f32; 3]) -> f32 {
        self.blobs
            .iter()
            .map(|b| {
                let dx = p[0] - b.center[0];
                let dy = p[1] - b.center[1];
                let dz = p[2] - b.center[2];
                let r2 = dx * dx + dy * dy + dz * dz;
                let r = r2.sqrt();
                if r >= b.radius {
                    0.0
                } else {
                    let t = 1.0 - r / b.radius;
                    b.strength * t * t
                }
            })
            .sum()
    }
    /// Extract an isosurface at `threshold` using a simple grid march.
    ///
    /// `grid_dims` is `[nx, ny, nz]`.  Returns `Vec<[[f32; 3\]; 3]>` triangles.
    pub fn extract_isosurface(
        &self,
        box_min: [f32; 3],
        box_max: [f32; 3],
        grid_dims: [usize; 3],
        threshold: f32,
    ) -> Vec<[[f32; 3]; 3]> {
        let [nx, ny, nz] = grid_dims;
        let sx = (box_max[0] - box_min[0]) / nx as f32;
        let sy = (box_max[1] - box_min[1]) / ny as f32;
        let sz = (box_max[2] - box_min[2]) / nz as f32;
        let mut tris = Vec::new();
        for ix in 0..nx.saturating_sub(1) {
            for iy in 0..ny.saturating_sub(1) {
                for iz in 0..nz.saturating_sub(1) {
                    let cx = box_min[0] + (ix as f32 + 0.5) * sx;
                    let cy = box_min[1] + (iy as f32 + 0.5) * sy;
                    let cz = box_min[2] + (iz as f32 + 0.5) * sz;
                    let mut inside = 0u8;
                    let corners = [
                        [cx, cy, cz],
                        [cx + sx, cy, cz],
                        [cx + sx, cy + sy, cz],
                        [cx, cy + sy, cz],
                        [cx, cy, cz + sz],
                        [cx + sx, cy, cz + sz],
                        [cx + sx, cy + sy, cz + sz],
                        [cx, cy + sy, cz + sz],
                    ];
                    for (k, &c) in corners.iter().enumerate() {
                        if self.evaluate(c) >= threshold {
                            inside |= 1 << k;
                        }
                    }
                    if inside == 0 || inside == 0xFF {
                        continue;
                    }
                    let h = sx * 0.3;
                    tris.push([[cx - h, cy, cz], [cx + h, cy, cz], [cx, cy + h, cz]]);
                }
            }
        }
        tris
    }
}
impl MetaballField {
    /// Add a particle (blob) at `pos` with radius `r` and unit strength.
    pub fn add_particle(&mut self, pos: [f32; 3], r: f32) {
        self.add(pos, 1.0, r);
    }
    /// Sample the metaball density field at point `p` (alias for `evaluate`).
    pub fn sample_field(&self, p: [f32; 3]) -> f64 {
        self.evaluate(p) as f64
    }
}
/// Sorting buffer that maintains a depth-sorted index list for back-to-front
/// alpha-blended rendering.
#[derive(Debug, Clone, Default)]
pub struct ParticleSortBuffer {
    /// Sorted indices into a particle array.
    pub indices: Vec<usize>,
    /// Corresponding view-space depths (positive = in front of camera).
    pub depths: Vec<f64>,
}
impl ParticleSortBuffer {
    /// Create an empty sort buffer.
    pub fn new() -> Self {
        Self {
            indices: Vec::new(),
            depths: Vec::new(),
        }
    }
    /// Rebuild the index/depth arrays and sort back-to-front.
    ///
    /// `camera_pos` is the world-space camera position.
    pub fn sort_back_to_front(&mut self, camera_pos: [f64; 3], positions: &[[f64; 3]]) {
        let n = positions.len();
        self.indices = (0..n).collect();
        self.depths = positions
            .iter()
            .map(|p| {
                let dx = p[0] - camera_pos[0];
                let dy = p[1] - camera_pos[1];
                let dz = p[2] - camera_pos[2];
                dx * dx + dy * dy + dz * dz
            })
            .collect();
        self.indices.sort_by(|&a, &b| {
            self.depths[b]
                .partial_cmp(&self.depths[a])
                .unwrap_or(std::cmp::Ordering::Equal)
        });
    }
}
/// A billboard particle with f64 position and radius for high-precision use.
#[derive(Debug, Clone, PartialEq)]
pub struct BillboardSprite {
    /// World-space position \[x, y, z\] in f64.
    pub pos: [f64; 3],
    /// Sphere radius in world units.
    pub radius: f64,
    /// RGBA color, each component in \[0, 1\].
    pub color: [f32; 4],
    /// Per-particle opacity override.
    pub opacity: f32,
}
/// A sprite entry holding position, velocity, color and computed screen size.
#[derive(Debug, Clone)]
pub struct PointSprite {
    /// World-space position.
    pub position: [f32; 3],
    /// Velocity (used if `use_velocity_stretch` is enabled).
    pub velocity: [f32; 3],
    /// RGBA color.
    pub color: [f32; 4],
    /// Computed world-space size for this sprite.
    pub size: f32,
}
/// Pre-defined color maps for scalar visualization.
#[derive(Debug, Clone, PartialEq)]
pub enum ColorMap {
    /// Perceptually uniform blue→green→yellow palette.
    Viridis,
    /// Perceptually uniform purple→orange→yellow palette.
    Plasma,
    /// Black→red→yellow→white thermal ramp.
    Hot,
    /// Cyan→magenta linear blend.
    Cool,
    /// Linear grey ramp.
    Grayscale,
}
impl ColorMap {
    /// Map scalar parameter `t` ∈ \[0, 1\] to an RGBA color.
    pub fn map(&self, t: f32) -> [f32; 4] {
        match self {
            ColorMap::Viridis => {
                let stops: [[f32; 4]; 5] = [
                    [0.0, 0.0, 0.3, 1.0],
                    [0.0, 0.25, 0.8, 0.4],
                    [0.2, 0.5, 0.7, 0.5],
                    [0.5, 0.75, 0.8, 0.3],
                    [0.99, 1.0, 0.9, 0.14],
                ];
                piecewise_lerp(&stops, t)
            }
            ColorMap::Plasma => {
                let stops: [[f32; 4]; 5] = [
                    [0.05, 0.03, 0.53, 1.0],
                    [0.46, 0.00, 0.66, 1.0],
                    [0.80, 0.10, 0.46, 1.0],
                    [0.97, 0.46, 0.10, 1.0],
                    [0.94, 0.98, 0.13, 1.0],
                ];
                piecewise_lerp(&stops, t)
            }
            ColorMap::Hot => {
                let stops: [[f32; 4]; 4] = [
                    [0.0, 0.0, 0.0, 1.0],
                    [1.0, 0.0, 0.0, 1.0],
                    [1.0, 1.0, 0.0, 1.0],
                    [1.0, 1.0, 1.0, 1.0],
                ];
                piecewise_lerp(&stops, t)
            }
            ColorMap::Cool => {
                let t = t.clamp(0.0, 1.0);
                [t, 1.0 - t, 1.0, 1.0]
            }
            ColorMap::Grayscale => {
                let t = t.clamp(0.0, 1.0);
                [t, t, t, 1.0]
            }
        }
    }
}
/// SPH (Smoothed Particle Hydrodynamics) density field evaluator.
pub struct SphDensityField {
    /// Smoothing length h.
    pub h: f32,
}
impl SphDensityField {
    /// Create a new SPH density field with smoothing length `h`.
    pub fn new(h: f32) -> Self {
        Self { h }
    }
    /// Cubic spline kernel W(r, h).
    fn kernel(&self, r: f32) -> f32 {
        let q = r / self.h;
        let sigma = 1.0 / (std::f32::consts::PI * self.h * self.h * self.h);
        if q < 1.0 {
            sigma * (1.0 - 1.5 * q * q + 0.75 * q * q * q)
        } else if q < 2.0 {
            sigma * 0.25 * (2.0 - q).powi(3)
        } else {
            0.0
        }
    }
    /// Evaluate the SPH density at point `p` given a list of `(position, mass)` particles.
    pub fn density(&self, p: [f32; 3], particles: &[([f32; 3], f32)]) -> f32 {
        particles
            .iter()
            .map(|(pos, mass)| {
                let dx = p[0] - pos[0];
                let dy = p[1] - pos[1];
                let dz = p[2] - pos[2];
                let r = (dx * dx + dy * dy + dz * dz).sqrt();
                mass * self.kernel(r)
            })
            .sum()
    }
    /// Generate a 2D slice of density values at z=0 on a uniform grid.
    ///
    /// Returns `[nx * ny]` density values.
    pub fn slice_z0(
        &self,
        particles: &[([f32; 3], f32)],
        xrange: [f32; 2],
        yrange: [f32; 2],
        nx: usize,
        ny: usize,
    ) -> Vec<f32> {
        let mut grid = Vec::with_capacity(nx * ny);
        for iy in 0..ny {
            for ix in 0..nx {
                let x = xrange[0] + (ix as f32 / (nx - 1).max(1) as f32) * (xrange[1] - xrange[0]);
                let y = yrange[0] + (iy as f32 / (ny - 1).max(1) as f32) * (yrange[1] - yrange[0]);
                grid.push(self.density([x, y, 0.0], particles));
            }
        }
        grid
    }
}
/// Level of detail system for particle rendering.
///
/// `thresholds[i]` is the camera distance above which LOD level `i+1` is used.
pub struct LodSystem {
    /// Distance thresholds in ascending order.
    pub thresholds: Vec<f32>,
}
impl LodSystem {
    /// Create a LOD system with the given distance thresholds.
    pub fn new(thresholds: Vec<f32>) -> Self {
        Self { thresholds }
    }
    /// Return the LOD level (0 = highest detail) for the given camera distance.
    pub fn level_for_distance(&self, distance: f32) -> usize {
        for (i, &t) in self.thresholds.iter().enumerate() {
            if distance < t {
                return i;
            }
        }
        self.thresholds.len().saturating_sub(1)
    }
    /// Return an appropriate render radius scaled by LOD level.
    pub fn render_radius(&self, base_radius: f32, camera_distance: f32) -> f32 {
        let level = self.level_for_distance(camera_distance);
        let scale = 1.0 / (1.0 + level as f32 * 0.5);
        (base_radius * scale) / camera_distance.max(f32::EPSILON)
    }
    /// Return the fraction of particles to render at a given LOD level.
    ///
    /// Level 0 → 1.0 (all), level 1 → 0.5, level 2 → 0.25, etc.
    pub fn particle_fraction(&self, level: usize) -> f32 {
        1.0 / (1u32 << level.min(8)) as f32
    }
}
