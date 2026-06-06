//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::{cross, mip_ray};

/// Accumulates emission-absorption along a ray.
///
/// Each step contributes an emission color and an absorption coefficient.
pub struct EmissionAbsorptionAccumulator {
    /// Accumulated color `[r, g, b]`.
    pub color: [f32; 3],
    /// Accumulated transmittance (starts at 1, decreases toward 0).
    pub transmittance: f32,
}
impl EmissionAbsorptionAccumulator {
    /// Create a new accumulator (transparent, black).
    pub fn new() -> Self {
        Self {
            color: [0.0; 3],
            transmittance: 1.0,
        }
    }
    /// Add one step: `ds` is the step length, `emission` is the emitted radiance,
    /// `absorption` is the extinction coefficient.
    pub fn step(&mut self, ds: f32, emission: [f32; 3], absorption: f32) {
        let t = (-absorption * ds).exp();
        let contrib = self.transmittance * (1.0 - t);
        self.color[0] += contrib * emission[0];
        self.color[1] += contrib * emission[1];
        self.color[2] += contrib * emission[2];
        self.transmittance *= t;
    }
    /// Return `[r, g, b, alpha]` where `alpha = 1 - transmittance`.
    pub fn result(&self) -> [f32; 4] {
        [
            self.color[0],
            self.color[1],
            self.color[2],
            1.0 - self.transmittance,
        ]
    }
}
/// A triangle vertex produced by isosurface extraction.
#[derive(Debug, Clone, Copy)]
pub struct IsoVertex {
    /// Vertex position in world space.
    pub position: Vec3,
    /// Surface normal at the vertex.
    pub normal: Vec3,
}
/// A ray defined by an origin and a (unit) direction.
#[derive(Debug, Clone, Copy)]
pub struct Ray {
    /// Ray origin in world space.
    pub origin: Vec3,
    /// Unit direction of the ray.
    pub direction: Vec3,
}
impl Ray {
    /// Create a new ray from an origin and direction.
    pub fn new(origin: Vec3, direction: Vec3) -> Self {
        Self { origin, direction }
    }
    /// Return the point at parameter `t` along the ray.
    pub fn at(&self, t: f32) -> Vec3 {
        self.origin.add(&self.direction.scale(t))
    }
}
/// Ray marcher with adaptive step size based on local density gradient.
///
/// In regions of high gradient (boundaries), the step size is reduced for
/// better accuracy.  In homogeneous regions, larger steps are used.
pub struct AdaptiveRayMarcher {
    /// Base (minimum) step size.
    pub base_step: f32,
    /// Maximum step size (used in empty regions).
    pub max_step: f32,
    /// Early termination opacity threshold.
    pub opacity_threshold: f32,
    /// Gradient threshold: below this, use max_step; above, use base_step.
    pub gradient_threshold: f32,
}
impl AdaptiveRayMarcher {
    /// Create an adaptive ray marcher.
    pub fn new(base_step: f32, max_step: f32) -> Self {
        Self {
            base_step,
            max_step,
            opacity_threshold: 0.99,
            gradient_threshold: 0.05,
        }
    }
    /// March a ray through the volume with adaptive step size.
    ///
    /// Returns `[r, g, b, alpha]`.
    pub fn march(&self, ray: &Ray, volume: &Volume, tf: &TransferFunction) -> [f32; 4] {
        let Some((t_enter, t_exit)) = volume.aabb_intersect(ray) else {
            return [0.0; 4];
        };
        let mut t = t_enter;
        let mut acc = [0.0f32; 4];
        while t < t_exit {
            let p = ray.at(t);
            let density = volume.sample_trilinear(p);
            let sample = tf.evaluate(density);
            let p_fwd = ray.at(t + self.base_step);
            let d_fwd = volume.sample_trilinear(p_fwd);
            let gradient = (d_fwd - density).abs();
            let step = if gradient < self.gradient_threshold {
                self.max_step
            } else {
                self.base_step
            };
            let sa = sample[3];
            let oma = 1.0 - acc[3];
            acc[0] += oma * sample[0] * sa;
            acc[1] += oma * sample[1] * sa;
            acc[2] += oma * sample[2] * sa;
            acc[3] += oma * sa;
            if acc[3] >= self.opacity_threshold {
                break;
            }
            t += step;
        }
        acc
    }
}
/// Maps a scalar density value to an RGBA colour via piecewise-linear interpolation
/// over a sorted list of control points.
pub struct TransferFunction {
    /// Sorted by density (first element of the tuple).
    pub control_points: Vec<(f32, [f32; 4])>,
}
impl TransferFunction {
    /// Create an empty transfer function with no control points.
    pub fn new() -> Self {
        Self {
            control_points: Vec::new(),
        }
    }
    /// Add a control point.  The list is kept sorted by density.
    pub fn add_control_point(&mut self, density: f32, color: [f32; 3], opacity: f32) {
        let rgba = [color[0], color[1], color[2], opacity];
        self.control_points.push((density, rgba));
        self.control_points
            .sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
    }
    /// Evaluate the transfer function at `density` by linear interpolation.
    pub fn evaluate(&self, density: f32) -> [f32; 4] {
        if self.control_points.is_empty() {
            return [0.0; 4];
        }
        if density <= self.control_points[0].0 {
            return self.control_points[0].1;
        }
        let last = self
            .control_points
            .last()
            .expect("collection should not be empty");
        if density >= last.0 {
            return last.1;
        }
        let pos = self.control_points.partition_point(|cp| cp.0 <= density);
        let (d0, c0) = self.control_points[pos - 1];
        let (d1, c1) = self.control_points[pos];
        let t = (density - d0) / (d1 - d0);
        [
            c0[0] + t * (c1[0] - c0[0]),
            c0[1] + t * (c1[1] - c0[1]),
            c0[2] + t * (c1[2] - c0[2]),
            c0[3] + t * (c1[3] - c0[3]),
        ]
    }
    /// Preset transfer function for plasma / fire.
    pub fn plasma_transfer_function() -> Self {
        let mut tf = Self::new();
        tf.add_control_point(0.00, [0.0, 0.0, 0.0], 0.0);
        tf.add_control_point(0.10, [0.2, 0.0, 0.4], 0.02);
        tf.add_control_point(0.30, [0.8, 0.1, 0.0], 0.15);
        tf.add_control_point(0.55, [1.0, 0.5, 0.0], 0.40);
        tf.add_control_point(0.75, [1.0, 0.9, 0.1], 0.70);
        tf.add_control_point(1.00, [1.0, 1.0, 1.0], 1.00);
        tf
    }
    /// Preset transfer function for smoke.
    pub fn smoke_transfer_function() -> Self {
        let mut tf = Self::new();
        tf.add_control_point(0.00, [0.0, 0.0, 0.0], 0.0);
        tf.add_control_point(0.05, [0.5, 0.5, 0.5], 0.05);
        tf.add_control_point(0.30, [0.4, 0.4, 0.4], 0.25);
        tf.add_control_point(0.70, [0.2, 0.2, 0.2], 0.60);
        tf.add_control_point(1.00, [0.1, 0.1, 0.1], 0.90);
        tf
    }
}
/// An occupancy grid for fast empty-space skipping.
///
/// The volume is subdivided into bricks of size `brick_size³`.  A brick is
/// marked occupied if its maximum density exceeds a threshold.
pub struct OccupancyGrid {
    /// Number of bricks along each axis.
    pub nx: usize,
    /// Number of bricks along y axis.
    pub ny: usize,
    /// Number of bricks along z axis.
    pub nz: usize,
    /// Brick size in voxels.
    pub brick_size: usize,
    /// Flat array: `true` if the brick has density above threshold.
    pub occupied: Vec<bool>,
}
impl OccupancyGrid {
    /// Build an occupancy grid from a volume.
    pub fn build(volume: &Volume, brick_size: usize, threshold: f32) -> Self {
        let bs = brick_size.max(1);
        let nx = volume.nx.div_ceil(bs);
        let ny = volume.ny.div_ceil(bs);
        let nz = volume.nz.div_ceil(bs);
        let mut occupied = vec![false; nx * ny * nz];
        for bx in 0..nx {
            for by in 0..ny {
                for bz in 0..nz {
                    let x0 = bx * bs;
                    let y0 = by * bs;
                    let z0 = bz * bs;
                    let mut any_above = false;
                    'outer: for ix in x0..((x0 + bs).min(volume.nx)) {
                        for iy in y0..((y0 + bs).min(volume.ny)) {
                            for iz in z0..((z0 + bs).min(volume.nz)) {
                                if volume.get_value(ix, iy, iz) >= threshold {
                                    any_above = true;
                                    break 'outer;
                                }
                            }
                        }
                    }
                    occupied[bz * ny * nx + by * nx + bx] = any_above;
                }
            }
        }
        Self {
            nx,
            ny,
            nz,
            brick_size: bs,
            occupied,
        }
    }
    /// Check whether the world-space position `p` falls in an occupied brick.
    ///
    /// `volume` is used to map world-space to voxel coordinates.
    pub fn is_occupied_at(&self, volume: &Volume, p: Vec3) -> bool {
        let fx = (p.x - volume.origin.x) / volume.size.x * volume.nx as f32;
        let fy = (p.y - volume.origin.y) / volume.size.y * volume.ny as f32;
        let fz = (p.z - volume.origin.z) / volume.size.z * volume.nz as f32;
        if fx < 0.0
            || fy < 0.0
            || fz < 0.0
            || fx >= volume.nx as f32
            || fy >= volume.ny as f32
            || fz >= volume.nz as f32
        {
            return false;
        }
        let bx = (fx as usize) / self.brick_size;
        let by = (fy as usize) / self.brick_size;
        let bz = (fz as usize) / self.brick_size;
        let bx = bx.min(self.nx - 1);
        let by = by.min(self.ny - 1);
        let bz = bz.min(self.nz - 1);
        self.occupied[bz * self.ny * self.nx + by * self.nx + bx]
    }
    /// Fraction of bricks that are occupied.
    pub fn occupancy_ratio(&self) -> f32 {
        let total = self.occupied.len();
        if total == 0 {
            return 0.0;
        }
        let occ = self.occupied.iter().filter(|&&b| b).count();
        occ as f32 / total as f32
    }
}
/// A 3-component f32 vector used throughout the volume renderer.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Vec3 {
    /// X component.
    pub x: f32,
    /// Y component.
    pub y: f32,
    /// Z component.
    pub z: f32,
}
impl Vec3 {
    /// Create a new vector with the given components.
    pub fn new(x: f32, y: f32, z: f32) -> Self {
        Self { x, y, z }
    }
    /// Return the zero vector.
    pub fn zero() -> Self {
        Self {
            x: 0.0,
            y: 0.0,
            z: 0.0,
        }
    }
    /// Component-wise addition.
    pub fn add(&self, o: &Vec3) -> Vec3 {
        Vec3::new(self.x + o.x, self.y + o.y, self.z + o.z)
    }
    /// Component-wise subtraction.
    pub fn sub(&self, o: &Vec3) -> Vec3 {
        Vec3::new(self.x - o.x, self.y - o.y, self.z - o.z)
    }
    /// Uniform scaling.
    pub fn scale(&self, s: f32) -> Vec3 {
        Vec3::new(self.x * s, self.y * s, self.z * s)
    }
    /// Dot product.
    pub fn dot(&self, o: &Vec3) -> f32 {
        self.x * o.x + self.y * o.y + self.z * o.z
    }
    /// Euclidean length.
    pub fn length(&self) -> f32 {
        self.dot(self).sqrt()
    }
    /// Return a unit-length copy, or the zero vector if length is near zero.
    pub fn normalize(&self) -> Vec3 {
        let len = self.length();
        if len > 1e-12 {
            self.scale(1.0 / len)
        } else {
            Vec3::zero()
        }
    }
}
/// High-level renderer that ties together Camera, RayMarcher and TransferFunction.
pub struct VolumeRenderer {
    /// The camera used to generate primary rays.
    pub camera: Camera,
    /// The ray marcher that steps through the volume.
    pub marcher: RayMarcher,
    /// Transfer function mapping density to colour and opacity.
    pub transfer_function: TransferFunction,
    /// Background colour `[r, g, b]`.
    pub background: [f32; 3],
}
impl VolumeRenderer {
    /// Create a new volume renderer with default black background.
    pub fn new(camera: Camera, marcher: RayMarcher, tf: TransferFunction) -> Self {
        Self {
            camera,
            marcher,
            transfer_function: tf,
            background: [0.0, 0.0, 0.0],
        }
    }
    /// Render the volume into a `width * height` pixel buffer.
    /// Each element is `[r, g, b, a]` in row-major order.
    pub fn render(&self, volume: &Volume, width: usize, height: usize) -> Vec<[f32; 4]> {
        let mut pixels = Vec::with_capacity(width * height);
        let wf = width as f32;
        let hf = height as f32;
        for py in 0..height {
            for px in 0..width {
                let ray = self.camera.generate_ray(px as f32, py as f32, wf, hf);
                let rgba = self.marcher.march(&ray, volume, &self.transfer_function);
                let one_minus_a = 1.0 - rgba[3];
                let out = [
                    rgba[0] + one_minus_a * self.background[0],
                    rgba[1] + one_minus_a * self.background[1],
                    rgba[2] + one_minus_a * self.background[2],
                    rgba[3],
                ];
                pixels.push(out);
            }
        }
        pixels
    }
}
/// A MIP renderer that generates images from a volume.
pub struct MipRenderer {
    /// Camera generating rays.
    pub camera: Camera,
    /// Step size for marching.
    pub step_size: f32,
}
impl MipRenderer {
    /// Create a new MIP renderer.
    pub fn new(camera: Camera, step_size: f32) -> Self {
        Self { camera, step_size }
    }
    /// Render using MIP mode.
    pub fn render(
        &self,
        volume: &Volume,
        tf: &TransferFunction,
        width: usize,
        height: usize,
    ) -> Vec<[f32; 4]> {
        let mut pixels = Vec::with_capacity(width * height);
        for py in 0..height {
            for px in 0..width {
                let ray =
                    self.camera
                        .generate_ray(px as f32, py as f32, width as f32, height as f32);
                pixels.push(mip_ray(&ray, volume, tf, self.step_size));
            }
        }
        pixels
    }
}
/// A direct volume rendering accumulator that combines color and density
/// via emission-absorption integration.
///
/// This models each voxel as both emitting (providing color) and absorbing
/// (reducing transmittance) light along the ray.
pub struct DvrAccumulator {
    /// Accumulated RGBA color.
    pub color: [f32; 4],
    /// Remaining transmittance (starts at 1.0).
    pub transmittance: f32,
    /// Total number of integration steps taken.
    pub steps: usize,
}
impl DvrAccumulator {
    /// Create a fresh accumulator.
    pub fn new() -> Self {
        Self {
            color: [0.0; 4],
            transmittance: 1.0,
            steps: 0,
        }
    }
    /// Integrate one step of emission-absorption with step size `ds`.
    ///
    /// `rgba` is the sample from the transfer function; `absorption_scale`
    /// multiplies the optical depth.
    pub fn integrate_step(&mut self, rgba: [f32; 4], ds: f32, absorption_scale: f32) {
        let optical_depth = rgba[3] * absorption_scale * ds;
        let t = (-optical_depth).exp();
        let contrib = self.transmittance * (1.0 - t);
        self.color[0] += contrib * rgba[0];
        self.color[1] += contrib * rgba[1];
        self.color[2] += contrib * rgba[2];
        self.color[3] = 1.0 - (self.color[3].max(0.0) + contrib).clamp(0.0, 1.0);
        self.transmittance *= t;
        self.steps += 1;
    }
    /// Finalize: return `[r, g, b, 1.0 - transmittance]`.
    pub fn result(&self) -> [f32; 4] {
        [
            self.color[0].clamp(0.0, 1.0),
            self.color[1].clamp(0.0, 1.0),
            self.color[2].clamp(0.0, 1.0),
            (1.0 - self.transmittance).clamp(0.0, 1.0),
        ]
    }
    /// Returns `true` if transmittance has dropped below `threshold`.
    pub fn is_terminated(&self, threshold: f32) -> bool {
        self.transmittance < threshold
    }
}
/// Performs front-to-back compositing along a ray through a `Volume`.
pub struct RayMarcher {
    /// Distance between successive sample points along the ray.
    pub step_size: f32,
    /// Terminate early when accumulated opacity exceeds this threshold.
    pub max_opacity: f32,
}
impl RayMarcher {
    /// Create a new ray marcher with the given step size.
    pub fn new(step_size: f32) -> Self {
        Self {
            step_size,
            max_opacity: 0.99,
        }
    }
    /// March `ray` through `volume`, evaluate `tf` at each sample, and composite.
    /// Returns `[r, g, b, alpha]`.
    pub fn march(&self, ray: &Ray, volume: &Volume, tf: &TransferFunction) -> [f32; 4] {
        let Some((t_enter, t_exit)) = volume.aabb_intersect(ray) else {
            return [0.0, 0.0, 0.0, 0.0];
        };
        let mut t = t_enter;
        let mut acc_r = 0.0_f32;
        let mut acc_g = 0.0_f32;
        let mut acc_b = 0.0_f32;
        let mut acc_a = 0.0_f32;
        while t < t_exit {
            let p = ray.at(t);
            let density = volume.sample_trilinear(p);
            let sample = tf.evaluate(density);
            let sr = sample[0];
            let sg = sample[1];
            let sb = sample[2];
            let sa = sample[3];
            let one_minus_a = 1.0 - acc_a;
            acc_r += one_minus_a * sr * sa;
            acc_g += one_minus_a * sg * sa;
            acc_b += one_minus_a * sb * sa;
            acc_a += one_minus_a * sa;
            if acc_a >= self.max_opacity {
                break;
            }
            t += self.step_size;
        }
        [acc_r, acc_g, acc_b, acc_a]
    }
}
/// A ray marcher that uses an [`OccupancyGrid`] to skip empty bricks.
pub struct EmptySpaceSkippingMarcher {
    /// Step size inside occupied bricks.
    pub fine_step: f32,
    /// Step size in empty bricks (coarser).
    pub coarse_step: f32,
    /// Opacity threshold for early termination.
    pub opacity_threshold: f32,
}
impl EmptySpaceSkippingMarcher {
    /// Create with the given step sizes.
    pub fn new(fine_step: f32, coarse_step: f32) -> Self {
        Self {
            fine_step,
            coarse_step,
            opacity_threshold: 0.99,
        }
    }
    /// March the ray, skipping empty bricks with larger steps.
    pub fn march(
        &self,
        ray: &Ray,
        volume: &Volume,
        tf: &TransferFunction,
        occupancy: &OccupancyGrid,
    ) -> ([f32; 4], usize) {
        let Some((t_enter, t_exit)) = volume.aabb_intersect(ray) else {
            return ([0.0; 4], 0);
        };
        let mut t = t_enter;
        let mut acc = [0.0f32; 4];
        let mut steps = 0usize;
        while t < t_exit {
            let p = ray.at(t);
            let step = if occupancy.is_occupied_at(volume, p) {
                self.fine_step
            } else {
                self.coarse_step
            };
            if occupancy.is_occupied_at(volume, p) {
                let density = volume.sample_trilinear(p);
                let sample = tf.evaluate(density);
                let sa = sample[3];
                let oma = 1.0 - acc[3];
                acc[0] += oma * sample[0] * sa;
                acc[1] += oma * sample[1] * sa;
                acc[2] += oma * sample[2] * sa;
                acc[3] += oma * sa;
                if acc[3] >= self.opacity_threshold {
                    break;
                }
            }
            steps += 1;
            t += step;
        }
        (acc, steps)
    }
}
/// Ray marcher with configurable early termination and step statistics.
pub struct EarlyTerminationRayMarcher {
    /// Step size for marching.
    pub step_size: f32,
    /// Stop when accumulated opacity exceeds this threshold.
    pub opacity_threshold: f32,
}
impl EarlyTerminationRayMarcher {
    /// March a ray through the volume and return `([r,g,b,a], step_count)`.
    pub fn march_with_stats(
        &self,
        ray: &Ray,
        volume: &Volume,
        tf: &TransferFunction,
    ) -> ([f32; 4], usize) {
        let Some((t_enter, t_exit)) = volume.aabb_intersect(ray) else {
            return ([0.0, 0.0, 0.0, 0.0], 0);
        };
        let mut t = t_enter;
        let mut acc = [0.0f32; 4];
        let mut steps = 0usize;
        while t < t_exit {
            let p = ray.at(t);
            let density = volume.sample_trilinear(p);
            let sample = tf.evaluate(density);
            let sa = sample[3];
            let oma = 1.0 - acc[3];
            acc[0] += oma * sample[0] * sa;
            acc[1] += oma * sample[1] * sa;
            acc[2] += oma * sample[2] * sa;
            acc[3] += oma * sa;
            steps += 1;
            if acc[3] >= self.opacity_threshold {
                break;
            }
            t += self.step_size;
        }
        (acc, steps)
    }
}
/// A perspective camera that generates rays for volume rendering.
pub struct Camera {
    /// Camera position in world space.
    pub position: Vec3,
    /// Look-at target point.
    pub target: Vec3,
    /// Up direction hint.
    pub up: Vec3,
    /// Vertical field-of-view in radians.
    pub fov_y: f32,
}
impl Camera {
    /// Create a new perspective camera.
    pub fn new(pos: Vec3, target: Vec3, up: Vec3, fov_y: f32) -> Self {
        Self {
            position: pos,
            target,
            up,
            fov_y,
        }
    }
    /// Generate a ray through pixel (px, py) for an image of size (width, height).
    /// px and py are in `[0, width)` and `[0, height)` respectively (top-left origin).
    pub fn generate_ray(&self, px: f32, py: f32, width: f32, height: f32) -> Ray {
        let forward = self.target.sub(&self.position).normalize();
        let right = cross(&forward, &self.up).normalize();
        let up = cross(&right, &forward).normalize();
        let aspect = width / height;
        let half_h = (self.fov_y * 0.5).tan();
        let half_w = aspect * half_h;
        let ndc_x = (2.0 * (px + 0.5) / width - 1.0) * half_w;
        let ndc_y = (1.0 - 2.0 * (py + 0.5) / height) * half_h;
        let dir = forward
            .add(&right.scale(ndc_x))
            .add(&up.scale(ndc_y))
            .normalize();
        Ray::new(self.position, dir)
    }
}
/// Transfer function mapping scalar density values to RGBA via linear
/// interpolation over a list of control points (f64 density, \[f32;4\] colour).
#[derive(Debug, Clone)]
pub struct TransferFunctionF64 {
    /// Control points sorted by density value (ascending).
    pub control_points: Vec<(f64, [f32; 4])>,
}
impl TransferFunctionF64 {
    /// Create an empty transfer function.
    pub fn new() -> Self {
        Self {
            control_points: Vec::new(),
        }
    }
    /// Add a control point; list is kept sorted.
    pub fn add_point(&mut self, density: f64, color: [f32; 4]) {
        self.control_points.push((density, color));
        self.control_points
            .sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
    }
    /// Sample the transfer function at `value` by linear interpolation.
    pub fn sample(&self, value: f64) -> [f32; 4] {
        if self.control_points.is_empty() {
            return [0.0; 4];
        }
        if value <= self.control_points[0].0 {
            return self.control_points[0].1;
        }
        let last = self
            .control_points
            .last()
            .expect("collection should not be empty");
        if value >= last.0 {
            return last.1;
        }
        let pos = self.control_points.partition_point(|cp| cp.0 <= value);
        let (d0, c0) = self.control_points[pos - 1];
        let (d1, c1) = self.control_points[pos];
        let t = ((value - d0) / (d1 - d0)) as f32;
        [
            c0[0] + t * (c1[0] - c0[0]),
            c0[1] + t * (c1[1] - c0[1]),
            c0[2] + t * (c1[2] - c0[2]),
            c0[3] + t * (c1[3] - c0[3]),
        ]
    }
}
/// A brick-based volume that subdivides a volume into smaller bricks for
/// cache-friendly access and occlusion culling.
pub struct BrickVolume {
    /// Full volume dimensions along x.
    pub nx: usize,
    /// Full volume dimensions along y.
    pub ny: usize,
    /// Full volume dimensions along z.
    pub nz: usize,
    /// Brick size (all bricks are `brick_size` voxels on a side).
    pub brick_size: usize,
    /// Number of bricks along x.
    pub n_bricks_x: usize,
    /// Number of bricks along y.
    pub n_bricks_y: usize,
    /// Number of bricks along z.
    pub n_bricks_z: usize,
    /// Flat storage: brick-major, then voxel-major within each brick.
    pub data: Vec<f32>,
    /// World-space origin.
    pub origin: Vec3,
    /// World-space size.
    pub size: Vec3,
}
impl BrickVolume {
    /// Create a new zero-filled brick volume.
    pub fn new(
        nx: usize,
        ny: usize,
        nz: usize,
        brick_size: usize,
        origin: Vec3,
        size: Vec3,
    ) -> Self {
        let bs = brick_size.max(1);
        let nbx = nx.div_ceil(bs);
        let nby = ny.div_ceil(bs);
        let nbz = nz.div_ceil(bs);
        let total = nbx * nby * nbz * bs * bs * bs;
        Self {
            nx,
            ny,
            nz,
            brick_size: bs,
            n_bricks_x: nbx,
            n_bricks_y: nby,
            n_bricks_z: nbz,
            data: vec![0.0; total],
            origin,
            size,
        }
    }
    /// Flat index into the brick-layout data array.
    fn brick_index(&self, ix: usize, iy: usize, iz: usize) -> usize {
        let bs = self.brick_size;
        let bx = ix / bs;
        let by = iy / bs;
        let bz = iz / bs;
        let lx = ix % bs;
        let ly = iy % bs;
        let lz = iz % bs;
        let brick_flat = bz * self.n_bricks_y * self.n_bricks_x + by * self.n_bricks_x + bx;
        brick_flat * bs * bs * bs + lz * bs * bs + ly * bs + lx
    }
    /// Set a voxel value.
    pub fn set_value(&mut self, ix: usize, iy: usize, iz: usize, v: f32) {
        let idx = self.brick_index(ix, iy, iz);
        if idx < self.data.len() {
            self.data[idx] = v;
        }
    }
    /// Get a voxel value.
    pub fn get_value(&self, ix: usize, iy: usize, iz: usize) -> f32 {
        let idx = self.brick_index(ix, iy, iz);
        if idx < self.data.len() {
            self.data[idx]
        } else {
            0.0
        }
    }
    /// Fill all voxels with a constant value.
    pub fn fill(&mut self, v: f32) {
        for x in &mut self.data {
            *x = v;
        }
    }
    /// Total number of bricks.
    pub fn total_bricks(&self) -> usize {
        self.n_bricks_x * self.n_bricks_y * self.n_bricks_z
    }
    /// Convert a BrickVolume to a regular Volume.
    pub fn to_volume(&self) -> Volume {
        let mut vol = Volume::new(self.nx, self.ny, self.nz, self.origin, self.size);
        for ix in 0..self.nx {
            for iy in 0..self.ny {
                for iz in 0..self.nz {
                    vol.set_value(ix, iy, iz, self.get_value(ix, iy, iz));
                }
            }
        }
        vol
    }
}
/// A 3-D scalar density field stored as a flat array (row-major z-fastest).
pub struct Volume {
    /// Number of cells along x.
    pub nx: usize,
    /// Number of cells along y.
    pub ny: usize,
    /// Number of cells along z.
    pub nz: usize,
    /// Density value for each (ix, iy, iz) cell.
    pub values: Vec<f32>,
    /// World-space position of the (0,0,0) corner.
    pub origin: Vec3,
    /// World-space extents of the entire grid.
    pub size: Vec3,
}
impl Volume {
    /// Create a new zero-filled volume with the given dimensions.
    pub fn new(nx: usize, ny: usize, nz: usize, origin: Vec3, size: Vec3) -> Self {
        Self {
            nx,
            ny,
            nz,
            values: vec![0.0; nx * ny * nz],
            origin,
            size,
        }
    }
    /// Flat index: ix * ny * nz + iy * nz + iz
    pub fn idx(&self, ix: usize, iy: usize, iz: usize) -> usize {
        ix * self.ny * self.nz + iy * self.nz + iz
    }
    /// Set the density value at the given voxel indices.
    pub fn set_value(&mut self, ix: usize, iy: usize, iz: usize, v: f32) {
        let i = self.idx(ix, iy, iz);
        self.values[i] = v;
    }
    /// Get the density value at the given voxel indices.
    pub fn get_value(&self, ix: usize, iy: usize, iz: usize) -> f32 {
        let i = self.idx(ix, iy, iz);
        self.values[i]
    }
    /// Trilinear interpolation of density at world-space point `p`.
    /// Returns 0 for any point outside the volume bounding box.
    pub fn sample_trilinear(&self, p: Vec3) -> f32 {
        let fx = (p.x - self.origin.x) / self.size.x * self.nx as f32;
        let fy = (p.y - self.origin.y) / self.size.y * self.ny as f32;
        let fz = (p.z - self.origin.z) / self.size.z * self.nz as f32;
        if fx < 0.0
            || fy < 0.0
            || fz < 0.0
            || fx >= self.nx as f32
            || fy >= self.ny as f32
            || fz >= self.nz as f32
        {
            return 0.0;
        }
        let ix0 = fx.floor() as usize;
        let iy0 = fy.floor() as usize;
        let iz0 = fz.floor() as usize;
        let ix1 = (ix0 + 1).min(self.nx - 1);
        let iy1 = (iy0 + 1).min(self.ny - 1);
        let iz1 = (iz0 + 1).min(self.nz - 1);
        let tx = fx - ix0 as f32;
        let ty = fy - iy0 as f32;
        let tz = fz - iz0 as f32;
        let c000 = self.get_value(ix0, iy0, iz0);
        let c001 = self.get_value(ix0, iy0, iz1);
        let c010 = self.get_value(ix0, iy1, iz0);
        let c011 = self.get_value(ix0, iy1, iz1);
        let c100 = self.get_value(ix1, iy0, iz0);
        let c101 = self.get_value(ix1, iy0, iz1);
        let c110 = self.get_value(ix1, iy1, iz0);
        let c111 = self.get_value(ix1, iy1, iz1);
        let c00 = c000 * (1.0 - tz) + c001 * tz;
        let c01 = c010 * (1.0 - tz) + c011 * tz;
        let c10 = c100 * (1.0 - tz) + c101 * tz;
        let c11 = c110 * (1.0 - tz) + c111 * tz;
        let c0 = c00 * (1.0 - ty) + c01 * ty;
        let c1 = c10 * (1.0 - ty) + c11 * ty;
        c0 * (1.0 - tx) + c1 * tx
    }
    /// Axis-aligned bounding-box slab test.
    /// Returns `Some((t_enter, t_exit))` when the ray intersects the box, or `None`.
    pub fn aabb_intersect(&self, ray: &Ray) -> Option<(f32, f32)> {
        let min = self.origin;
        let max = Vec3::new(
            self.origin.x + self.size.x,
            self.origin.y + self.size.y,
            self.origin.z + self.size.z,
        );
        let inv = |d: f32| {
            if d.abs() > 1e-12 {
                1.0 / d
            } else {
                f32::INFINITY
            }
        };
        let inv_dx = inv(ray.direction.x);
        let inv_dy = inv(ray.direction.y);
        let inv_dz = inv(ray.direction.z);
        let tx1 = (min.x - ray.origin.x) * inv_dx;
        let tx2 = (max.x - ray.origin.x) * inv_dx;
        let ty1 = (min.y - ray.origin.y) * inv_dy;
        let ty2 = (max.y - ray.origin.y) * inv_dy;
        let tz1 = (min.z - ray.origin.z) * inv_dz;
        let tz2 = (max.z - ray.origin.z) * inv_dz;
        let t_enter = tx1.min(tx2).max(ty1.min(ty2)).max(tz1.min(tz2));
        let t_exit = tx1.max(tx2).min(ty1.max(ty2)).min(tz1.max(tz2));
        if t_exit > t_enter && t_exit > 0.0 {
            Some((t_enter.max(0.0), t_exit))
        } else {
            None
        }
    }
}
impl Volume {
    /// Compute minimum and maximum density values in the volume.
    pub fn min_max(&self) -> (f32, f32) {
        let mut min_v = f32::MAX;
        let mut max_v = f32::MIN;
        for &v in &self.values {
            if v < min_v {
                min_v = v;
            }
            if v > max_v {
                max_v = v;
            }
        }
        (min_v, max_v)
    }
    /// Compute the average density.
    pub fn average_density(&self) -> f32 {
        if self.values.is_empty() {
            return 0.0;
        }
        let sum: f32 = self.values.iter().sum();
        sum / self.values.len() as f32
    }
    /// Total number of voxels.
    pub fn voxel_count(&self) -> usize {
        self.nx * self.ny * self.nz
    }
    /// World-space cell size along each axis.
    pub fn cell_size(&self) -> Vec3 {
        Vec3::new(
            self.size.x / self.nx as f32,
            self.size.y / self.ny as f32,
            self.size.z / self.nz as f32,
        )
    }
    /// Set all values to a constant.
    pub fn fill(&mut self, value: f32) {
        for v in &mut self.values {
            *v = value;
        }
    }
    /// Add a sphere of constant density to the volume.
    pub fn add_sphere(&mut self, center: Vec3, radius: f32, density: f32) {
        let r2 = radius * radius;
        for ix in 0..self.nx {
            for iy in 0..self.ny {
                for iz in 0..self.nz {
                    let x = self.origin.x + (ix as f32 + 0.5) * self.size.x / self.nx as f32;
                    let y = self.origin.y + (iy as f32 + 0.5) * self.size.y / self.ny as f32;
                    let z = self.origin.z + (iz as f32 + 0.5) * self.size.z / self.nz as f32;
                    let dx = x - center.x;
                    let dy = y - center.y;
                    let dz = z - center.z;
                    if dx * dx + dy * dy + dz * dz <= r2 {
                        let idx = self.idx(ix, iy, iz);
                        self.values[idx] += density;
                    }
                }
            }
        }
    }
}
impl Volume {
    /// Estimate the gradient at voxel `(ix, iy, iz)` using central differences.
    ///
    /// Returns `[gx, gy, gz]` in world-space units.
    pub fn gradient(&self, ix: usize, iy: usize, iz: usize) -> [f32; 3] {
        let cell = self.cell_size();
        let ix0 = ix.saturating_sub(1);
        let ix1 = (ix + 1).min(self.nx - 1);
        let iy0 = iy.saturating_sub(1);
        let iy1 = (iy + 1).min(self.ny - 1);
        let iz0 = iz.saturating_sub(1);
        let iz1 = (iz + 1).min(self.nz - 1);
        let dx =
            (self.get_value(ix1, iy, iz) - self.get_value(ix0, iy, iz)) / (2.0 * cell.x).max(1e-12);
        let dy =
            (self.get_value(ix, iy1, iz) - self.get_value(ix, iy0, iz)) / (2.0 * cell.y).max(1e-12);
        let dz =
            (self.get_value(ix, iy, iz1) - self.get_value(ix, iy, iz0)) / (2.0 * cell.z).max(1e-12);
        [dx, dy, dz]
    }
    /// Magnitude of the gradient at voxel `(ix, iy, iz)`.
    pub fn gradient_magnitude(&self, ix: usize, iy: usize, iz: usize) -> f32 {
        let g = self.gradient(ix, iy, iz);
        (g[0] * g[0] + g[1] * g[1] + g[2] * g[2]).sqrt()
    }
    /// Normalized surface normal at voxel `(ix, iy, iz)`.
    ///
    /// Returns the zero vector if the gradient magnitude is negligible.
    pub fn normal_at(&self, ix: usize, iy: usize, iz: usize) -> [f32; 3] {
        let g = self.gradient(ix, iy, iz);
        let mag = (g[0] * g[0] + g[1] * g[1] + g[2] * g[2]).sqrt();
        if mag < 1e-12 {
            [0.0; 3]
        } else {
            [g[0] / mag, g[1] / mag, g[2] / mag]
        }
    }
}
impl Volume {
    /// Compute a histogram of density values with `n_bins` bins.
    ///
    /// Returns a `Vec`u32` of length `n_bins` containing bin counts.
    pub fn histogram(&self, n_bins: usize) -> Vec<u32> {
        let n = n_bins.max(1);
        let (min_v, max_v) = self.min_max();
        let range = (max_v - min_v).max(1e-10);
        let mut bins = vec![0u32; n];
        for &v in &self.values {
            let t = ((v - min_v) / range).clamp(0.0, 1.0 - 1e-10);
            let idx = (t * n as f32) as usize;
            bins[idx.min(n - 1)] += 1;
        }
        bins
    }
    /// Normalize the volume values so that the maximum is 1.0.
    pub fn normalize(&mut self) {
        let (_, max_v) = self.min_max();
        if max_v < 1e-12 {
            return;
        }
        for v in &mut self.values {
            *v /= max_v;
        }
    }
    /// Clamp all density values to `\[lo, hi\]`.
    pub fn clamp_values(&mut self, lo: f32, hi: f32) {
        for v in &mut self.values {
            *v = v.clamp(lo, hi);
        }
    }
}
/// A 3-D scalar field stored as a flat `Vec`f64`, voxel-centred layout.
///
/// Index ordering: `ix * ny * nz + iy * nz + iz` (x-major).
#[derive(Debug, Clone)]
pub struct VolumeGrid {
    /// Data array: `nx * ny * nz` values.
    pub data: Vec<f64>,
    /// Grid dimension along x.
    pub nx: usize,
    /// Grid dimension along y.
    pub ny: usize,
    /// Grid dimension along z.
    pub nz: usize,
    /// World-space voxel size (assumed isotropic).
    pub voxel_size: f64,
}
impl VolumeGrid {
    /// Create a new zero-filled volume.
    pub fn new(nx: usize, ny: usize, nz: usize, voxel_size: f64) -> Self {
        Self {
            data: vec![0.0_f64; nx * ny * nz],
            nx,
            ny,
            nz,
            voxel_size,
        }
    }
    #[inline]
    pub(super) fn idx(&self, ix: usize, iy: usize, iz: usize) -> usize {
        ix * self.ny * self.nz + iy * self.nz + iz
    }
    /// Trilinear interpolation at world-space position `pos` (origin = (0,0,0)).
    ///
    /// Returns 0 for positions outside the grid.
    pub fn sample_trilinear(&self, pos: [f64; 3]) -> f64 {
        let fx = pos[0] / self.voxel_size;
        let fy = pos[1] / self.voxel_size;
        let fz = pos[2] / self.voxel_size;
        if fx < 0.0
            || fy < 0.0
            || fz < 0.0
            || fx >= self.nx as f64
            || fy >= self.ny as f64
            || fz >= self.nz as f64
        {
            return 0.0;
        }
        let ix0 = fx.floor() as usize;
        let iy0 = fy.floor() as usize;
        let iz0 = fz.floor() as usize;
        let ix1 = (ix0 + 1).min(self.nx - 1);
        let iy1 = (iy0 + 1).min(self.ny - 1);
        let iz1 = (iz0 + 1).min(self.nz - 1);
        let tx = fx - ix0 as f64;
        let ty = fy - iy0 as f64;
        let tz = fz - iz0 as f64;
        let c000 = self.data[self.idx(ix0, iy0, iz0)];
        let c001 = self.data[self.idx(ix0, iy0, iz1)];
        let c010 = self.data[self.idx(ix0, iy1, iz0)];
        let c011 = self.data[self.idx(ix0, iy1, iz1)];
        let c100 = self.data[self.idx(ix1, iy0, iz0)];
        let c101 = self.data[self.idx(ix1, iy0, iz1)];
        let c110 = self.data[self.idx(ix1, iy1, iz0)];
        let c111 = self.data[self.idx(ix1, iy1, iz1)];
        let c00 = c000 * (1.0 - tz) + c001 * tz;
        let c01 = c010 * (1.0 - tz) + c011 * tz;
        let c10 = c100 * (1.0 - tz) + c101 * tz;
        let c11 = c110 * (1.0 - tz) + c111 * tz;
        let c0 = c00 * (1.0 - ty) + c01 * ty;
        let c1 = c10 * (1.0 - ty) + c11 * ty;
        c0 * (1.0 - tx) + c1 * tx
    }
    /// Estimate the gradient using central differences.
    ///
    /// Uses a step of `voxel_size` in each direction.
    pub fn gradient_at(&self, pos: [f64; 3]) -> [f64; 3] {
        let h = self.voxel_size;
        let gx = (self.sample_trilinear([pos[0] + h, pos[1], pos[2]])
            - self.sample_trilinear([pos[0] - h, pos[1], pos[2]]))
            / (2.0 * h);
        let gy = (self.sample_trilinear([pos[0], pos[1] + h, pos[2]])
            - self.sample_trilinear([pos[0], pos[1] - h, pos[2]]))
            / (2.0 * h);
        let gz = (self.sample_trilinear([pos[0], pos[1], pos[2] + h])
            - self.sample_trilinear([pos[0], pos[1], pos[2] - h]))
            / (2.0 * h);
        [gx, gy, gz]
    }
}
/// Configuration for the `ray_march` function.
#[derive(Debug, Clone, Copy)]
pub struct RayMarchSettings {
    /// Step size along the ray in world units.
    pub step_size: f64,
    /// Maximum number of steps.
    pub max_steps: usize,
    /// Absorption coefficient (Beer-Lambert).
    pub absorption: f64,
    /// Scattering coefficient (unused in the simple integrator, provided for extension).
    pub scattering: f64,
}
