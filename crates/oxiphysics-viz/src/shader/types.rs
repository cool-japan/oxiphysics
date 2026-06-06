//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use std::collections::HashMap;
use std::f64::consts::PI;

use super::functions::{dot3, identity4, look_at, normalize3, normalize3_shader, sh9_basis};

/// A typed uniform buffer that stores named uniform values with ordering.
pub struct UniformBuffer {
    /// Uniform values keyed by name.
    pub values: HashMap<String, UniformValue>,
    /// Insertion order (for deterministic layout).
    pub order: Vec<String>,
    /// Binding point index.
    pub binding: u32,
}
impl UniformBuffer {
    /// Create a new empty uniform buffer at the given binding point.
    pub fn new(binding: u32) -> Self {
        Self {
            values: HashMap::new(),
            order: Vec::new(),
            binding,
        }
    }
    /// Insert or replace a uniform value.
    pub fn set(&mut self, name: &str, value: UniformValue) {
        let key = name.to_owned();
        if !self.values.contains_key(&key) {
            self.order.push(key.clone());
        }
        self.values.insert(key, value);
    }
    /// Look up a uniform value by name.
    pub fn get(&self, name: &str) -> Option<&UniformValue> {
        self.values.get(name)
    }
    /// Number of uniforms in this buffer.
    pub fn len(&self) -> usize {
        self.values.len()
    }
    /// Whether the buffer is empty.
    pub fn is_empty(&self) -> bool {
        self.values.is_empty()
    }
    /// Remove a uniform by name. Returns true if it was present.
    pub fn remove(&mut self, name: &str) -> bool {
        if self.values.remove(name).is_some() {
            self.order.retain(|n| n != name);
            true
        } else {
            false
        }
    }
    /// Clear all uniforms.
    pub fn clear(&mut self) {
        self.values.clear();
        self.order.clear();
    }
}
/// Orthographic projection parameters for a directional shadow map.
pub struct ShadowMapSetup {
    /// Light direction (normalized).
    pub light_dir: [f64; 3],
    /// Orthographic half-width.
    pub half_width: f64,
    /// Orthographic half-height.
    pub half_height: f64,
    /// Near clip plane distance.
    pub near: f64,
    /// Far clip plane distance.
    pub far: f64,
    /// Shadow map resolution (pixels).
    pub resolution: u32,
    /// Shadow bias to prevent acne.
    pub bias: f64,
}
impl ShadowMapSetup {
    /// Create a shadow map with reasonable defaults for a directional light.
    pub fn directional(light_dir: [f64; 3], extent: f64) -> Self {
        Self {
            light_dir: normalize3(light_dir),
            half_width: extent,
            half_height: extent,
            near: 0.1,
            far: extent * 4.0,
            resolution: 2048,
            bias: 0.005,
        }
    }
    /// Compute the orthographic projection matrix (column-major 4x4).
    pub fn ortho_matrix(&self) -> [[f64; 4]; 4] {
        let l = -self.half_width;
        let r = self.half_width;
        let b = -self.half_height;
        let t = self.half_height;
        let n = self.near;
        let f = self.far;
        let mut m = [[0.0_f64; 4]; 4];
        m[0][0] = 2.0 / (r - l);
        m[1][1] = 2.0 / (t - b);
        m[2][2] = -2.0 / (f - n);
        m[3][0] = -(r + l) / (r - l);
        m[3][1] = -(t + b) / (t - b);
        m[3][2] = -(f + n) / (f - n);
        m[3][3] = 1.0;
        m
    }
    /// Compute a simple look-at view matrix for the shadow camera.
    ///
    /// `center` is the world position the shadow camera looks at.
    pub fn view_matrix(&self, center: [f64; 3]) -> [[f64; 4]; 4] {
        let ld = self.light_dir;
        let eye = [
            center[0] - ld[0] * self.far * 0.5,
            center[1] - ld[1] * self.far * 0.5,
            center[2] - ld[2] * self.far * 0.5,
        ];
        look_at(eye, center, [0.0, 1.0, 0.0])
    }
    /// Texel size in world space.
    pub fn texel_size(&self) -> f64 {
        (2.0 * self.half_width) / (self.resolution as f64)
    }
}
/// A skeleton — a list of joints.
#[derive(Debug, Clone, Default)]
pub struct Skeleton {
    /// All joints, indexed by their position in this array.
    pub joints: Vec<Joint>,
}
impl Skeleton {
    /// Create an empty skeleton.
    pub fn new() -> Self {
        Self { joints: Vec::new() }
    }
    /// Add a joint and return its index.
    pub fn add_joint(&mut self, joint: Joint) -> usize {
        let idx = self.joints.len();
        self.joints.push(joint);
        idx
    }
    /// Number of joints.
    pub fn joint_count(&self) -> usize {
        self.joints.len()
    }
    /// Lookup a joint by name.
    pub fn find_joint(&self, name: &str) -> Option<usize> {
        self.joints.iter().position(|j| j.name == name)
    }
}
/// A joint (bone) in a skeleton, described by its inverse bind-pose matrix.
#[derive(Debug, Clone)]
pub struct Joint {
    /// Name of the joint.
    pub name: String,
    /// Index of the parent joint (`None` for root).
    pub parent: Option<usize>,
    /// Inverse bind-pose matrix (column-major 4x4).
    pub inv_bind: [[f64; 4]; 4],
}
impl Joint {
    /// Create a new joint with an identity inverse bind-pose matrix.
    pub fn new(name: &str, parent: Option<usize>) -> Self {
        Self {
            name: name.to_owned(),
            parent,
            inv_bind: identity4(),
        }
    }
}
/// Simple BSSRDF (bidirectional scattering-surface reflectance distribution
/// function) approximation for translucent materials.
///
/// Uses the Jensen 2001 dipole model scaled by pre-integrated profile weights.
#[derive(Debug, Clone)]
pub struct BssrdfApprox {
    /// Scattering coefficients per channel (σ_s).
    pub sigma_s: [f64; 3],
    /// Absorption coefficients per channel (σ_a).
    pub sigma_a: [f64; 3],
    /// Index of refraction.
    pub ior: f64,
}
impl BssrdfApprox {
    /// Reduced scattering coefficient: σ'_s = σ_s (1 - g) (Henyey-Greenstein g=0).
    pub fn sigma_s_prime(&self) -> [f64; 3] {
        self.sigma_s
    }
    /// Effective transport coefficient: σ_tr = sqrt(3 σ_a σ'_t).
    pub fn sigma_tr(&self) -> [f64; 3] {
        std::array::from_fn(|i| {
            let sigma_t_prime = self.sigma_s[i] + self.sigma_a[i];
            (3.0 * self.sigma_a[i] * sigma_t_prime).sqrt()
        })
    }
    /// Diffuse mean free path: l = 1 / σ_tr.
    pub fn mean_free_path(&self) -> [f64; 3] {
        self.sigma_tr()
            .map(|s| if s > 1e-30 { 1.0 / s } else { f64::MAX })
    }
    /// Evaluate the diffusion profile R(r) per channel.
    pub fn diffusion_profile(&self, r: f64) -> [f64; 3] {
        let sigma_tr = self.sigma_tr();
        let sp = self.sigma_s_prime();
        std::array::from_fn(|i| {
            if r < 1e-12 {
                return sp[i];
            }
            let sigma_t = sp[i] + self.sigma_a[i];
            let d = 1.0 / (3.0 * sigma_t.max(1e-30));
            (-sigma_tr[i] * r).exp() / (4.0 * PI * d * r)
        })
    }
    /// Total transmittance (fraction of light exiting after traversing distance d).
    pub fn transmittance(&self, distance: f64) -> [f64; 3] {
        let sigma_tr = self.sigma_tr();
        sigma_tr.map(|s| (-s * distance).exp())
    }
}
/// Options for generating a GLSL shader preamble.
#[derive(Debug, Clone)]
pub struct GlslOptions {
    /// GLSL version string, e.g. `"330 core"`.
    pub version: String,
    /// Whether to emit precision qualifiers (useful for ES targets).
    pub precision_mediump: bool,
    /// Extra `#define` macros (name -> optional value).
    pub defines: Vec<(String, Option<String>)>,
    /// Extra `#extension` directives (extension name -> behavior).
    pub extensions: Vec<(String, String)>,
}
impl GlslOptions {
    /// Default options for GLSL 330 core (desktop OpenGL).
    pub fn glsl_330() -> Self {
        Self {
            version: "330 core".to_owned(),
            precision_mediump: false,
            defines: Vec::new(),
            extensions: Vec::new(),
        }
    }
    /// Default options for GLSL ES 300 (WebGL 2 / mobile).
    pub fn glsl_es_300() -> Self {
        Self {
            version: "300 es".to_owned(),
            precision_mediump: true,
            defines: Vec::new(),
            extensions: Vec::new(),
        }
    }
    /// Add a `#define` without a value.
    pub fn define_flag(mut self, name: &str) -> Self {
        self.defines.push((name.to_owned(), None));
        self
    }
    /// Add a `#define` with a value.
    pub fn define(mut self, name: &str, value: &str) -> Self {
        self.defines.push((name.to_owned(), Some(value.to_owned())));
        self
    }
    /// Add a `#extension` directive.
    pub fn extension(mut self, name: &str, behavior: &str) -> Self {
        self.extensions.push((name.to_owned(), behavior.to_owned()));
        self
    }
    /// Generate the preamble string.
    pub fn preamble(&self) -> String {
        let mut out = format!("#version {}\n", self.version);
        for (ext, beh) in &self.extensions {
            out.push_str(&format!("#extension {ext} : {beh}\n"));
        }
        if self.precision_mediump {
            out.push_str("precision mediump float;\n");
        }
        for (name, val) in &self.defines {
            if let Some(v) = val {
                out.push_str(&format!("#define {name} {v}\n"));
            } else {
                out.push_str(&format!("#define {name}\n"));
            }
        }
        out
    }
}
/// A very simple image-based lighting probe encoded as 9 spherical harmonics
/// coefficients (L1 irradiance).
///
/// This represents ambient lighting that varies smoothly with surface normal.
#[derive(Debug, Clone)]
pub struct IblProbe {
    /// 9 SH coefficients for the R, G, B channels (L0 + L1 + L2, band 0..2).
    ///
    /// Layout: L00, L1-1, L10, L11, L2-2, L2-1, L20, L21, L22.
    pub sh9: [[f64; 3]; 9],
    /// Overall intensity multiplier.
    pub intensity: f64,
}
impl IblProbe {
    /// Create a probe with uniform grey illumination.
    pub fn uniform_grey(value: f64) -> Self {
        let mut sh9 = [[0.0_f64; 3]; 9];
        sh9[0] = [value; 3];
        Self {
            sh9,
            intensity: 1.0,
        }
    }
    /// Evaluate the irradiance in direction `dir` (does not need to be normalised).
    ///
    /// Uses a simplified 4-band (L0+L1) evaluation; higher L2 bands are zeroed
    /// if the caller only fills sh9\[0..4\].
    pub fn evaluate(&self, dir: [f64; 3]) -> [f64; 3] {
        let d = normalize3_shader(dir);
        let sh = sh9_basis(d);
        let mut result = [0.0_f64; 3];
        for (sh_i, sh9_i) in sh.iter().zip(self.sh9.iter()) {
            for (res_c, s9_c) in result.iter_mut().zip(sh9_i.iter()) {
                *res_c += sh_i * s9_c;
            }
        }
        [
            (result[0] * self.intensity).max(0.0),
            (result[1] * self.intensity).max(0.0),
            (result[2] * self.intensity).max(0.0),
        ]
    }
}
/// A pose: one local transform matrix per joint.
#[derive(Debug, Clone)]
pub struct Pose {
    /// Local joint transforms, one per joint.
    pub local_transforms: Vec<[[f64; 4]; 4]>,
}
impl Pose {
    /// Create a bind-pose (all identity matrices) for a skeleton of `n` joints.
    pub fn bind_pose(n: usize) -> Self {
        Self {
            local_transforms: vec![identity4(); n],
        }
    }
    /// Number of joints in this pose.
    pub fn joint_count(&self) -> usize {
        self.local_transforms.len()
    }
}
/// A named render pass that owns a collection of materials.
pub struct RenderPass {
    /// Human-readable name for this pass.
    pub name: String,
    /// Background clear colour (RGBA, linear).
    pub clear_color: [f64; 4],
    /// Materials submitted to this pass.
    pub materials: Vec<Material>,
}
impl RenderPass {
    /// Create an empty render pass with a black clear colour.
    pub fn new(name: &str) -> Self {
        Self {
            name: name.to_owned(),
            clear_color: [0.0, 0.0, 0.0, 1.0],
            materials: Vec::new(),
        }
    }
    /// Set the clear colour.
    pub fn set_clear_color(&mut self, r: f64, g: f64, b: f64, a: f64) {
        self.clear_color = [r, g, b, a];
    }
    /// Append a material to the pass.
    pub fn add_material(&mut self, m: Material) {
        self.materials.push(m);
    }
    /// Return the number of materials currently in the pass.
    pub fn material_count(&self) -> usize {
        self.materials.len()
    }
}
/// Per-vertex skinning data: joint indices and blend weights for up to 4 joints.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SkinningVertex {
    /// Index of up to 4 influencing joints.
    pub joint_indices: [u32; 4],
    /// Blend weights corresponding to each joint (should sum to ~1.0).
    pub weights: [f32; 4],
}
impl SkinningVertex {
    /// Create a vertex influenced by a single joint with full weight.
    pub fn single_joint(joint: u32) -> Self {
        Self {
            joint_indices: [joint, 0, 0, 0],
            weights: [1.0, 0.0, 0.0, 0.0],
        }
    }
    /// Create a vertex with two influences.
    pub fn two_joints(j0: u32, w0: f32, j1: u32, w1: f32) -> Self {
        Self {
            joint_indices: [j0, j1, 0, 0],
            weights: [w0, w1, 0.0, 0.0],
        }
    }
    /// Normalise the weights so they sum to 1.0.
    pub fn normalized(mut self) -> Self {
        let total: f32 = self.weights.iter().sum();
        if total > 1e-6 {
            for w in &mut self.weights {
                *w /= total;
            }
        }
        self
    }
    /// Sum of the weights (should be 1.0 after normalization).
    pub fn weight_sum(&self) -> f32 {
        self.weights.iter().sum()
    }
}
/// A simple PCF shadow sampler working on a floating-point depth map.
#[derive(Debug, Clone)]
pub struct PcfSampler {
    /// The depth map (row-major: `depths[row * width + col]`).
    pub depths: Vec<f32>,
    /// Width of the shadow map in texels.
    pub width: usize,
    /// Height of the shadow map in texels.
    pub height: usize,
    /// Shadow bias to reduce self-shadowing.
    pub bias: f32,
}
impl PcfSampler {
    /// Create a new PCF sampler from a depth map.
    pub fn new(depths: Vec<f32>, width: usize, height: usize, bias: f32) -> Self {
        Self {
            depths,
            width,
            height,
            bias,
        }
    }
    /// Create a sampler filled with a constant depth value.
    pub fn constant(depth: f32, width: usize, height: usize) -> Self {
        Self::new(vec![depth; width * height], width, height, 0.001)
    }
    /// Fetch the stored depth at texel `(u, v)` (both in \[0, 1\]).
    pub fn fetch_depth(&self, u: f64, v: f64) -> f32 {
        let col = ((u * self.width as f64) as usize).min(self.width.saturating_sub(1));
        let row = ((v * self.height as f64) as usize).min(self.height.saturating_sub(1));
        self.depths[row * self.width + col]
    }
    /// Sample PCF shadow factor at `(u, v)` for a fragment at shadow-space
    /// depth `receiver_depth`.
    ///
    /// `radius` controls the kernel radius in texels.
    /// Returns a value in \[0, 1\] where 0 = fully shadowed, 1 = fully lit.
    pub fn sample_pcf(&self, u: f64, v: f64, receiver_depth: f32, radius: usize) -> f64 {
        let r = radius as isize;
        let mut lit = 0.0_f64;
        let mut total = 0.0_f64;
        for dv in -r..=r {
            for du in -r..=r {
                let uu = u + du as f64 / self.width as f64;
                let vv = v + dv as f64 / self.height as f64;
                let uu = uu.clamp(0.0, 1.0 - 1e-6);
                let vv = vv.clamp(0.0, 1.0 - 1e-6);
                let d = self.fetch_depth(uu, vv);
                if receiver_depth - self.bias <= d {
                    lit += 1.0;
                }
                total += 1.0;
            }
        }
        if total > 0.0 { lit / total } else { 1.0 }
    }
    /// Sample with a 3x3 kernel (radius = 1).
    pub fn sample_pcf_3x3(&self, u: f64, v: f64, receiver_depth: f32) -> f64 {
        self.sample_pcf(u, v, receiver_depth, 1)
    }
}
/// Parameters for a particle sprite (billboard).
#[derive(Debug, Clone, Copy)]
pub struct ParticleSpriteParams {
    /// Particle radius in world space.
    pub radius: f64,
    /// Base colour.
    pub color: [f64; 4],
    /// Soft particle blend distance.
    pub soft_distance: f64,
    /// Sprite texture (0 = circle, 1 = gaussian, 2 = ring).
    pub sprite_type: u8,
    /// Glow factor (additive brightness).
    pub glow: f64,
}
impl ParticleSpriteParams {
    /// Default circular particle.
    pub fn circle(radius: f64, color: [f64; 4]) -> Self {
        Self {
            radius,
            color,
            soft_distance: 0.1,
            sprite_type: 0,
            glow: 0.0,
        }
    }
    /// Gaussian soft particle (for fire/smoke).
    pub fn gaussian(radius: f64, color: [f64; 4]) -> Self {
        Self {
            radius,
            color,
            soft_distance: 0.2,
            sprite_type: 1,
            glow: 0.3,
        }
    }
    /// Evaluate sprite alpha at local UV coordinate (u,v in \[-1,1\]).
    pub fn alpha(&self, u: f64, v: f64) -> f64 {
        let r2 = u * u + v * v;
        let r = r2.sqrt();
        match self.sprite_type {
            0 => {
                if r <= 1.0 {
                    self.color[3]
                } else {
                    0.0
                }
            }
            1 => {
                let g = (-3.0 * r2).exp();
                g * self.color[3]
            }
            2 => {
                let ring = 1.0 - (2.0 * r - 1.0).abs();
                ring.max(0.0) * self.color[3]
            }
            _ => {
                if r <= 1.0 {
                    self.color[3]
                } else {
                    0.0
                }
            }
        }
    }
    /// Soft particle depth blend: fades out near geometry.
    pub fn soft_blend(&self, particle_depth: f64, scene_depth: f64) -> f64 {
        let diff = (scene_depth - particle_depth).max(0.0);
        (diff / self.soft_distance.max(1e-12)).min(1.0)
    }
    /// Final colour including glow.
    pub fn final_color(&self, alpha: f64) -> [f64; 4] {
        let scale = 1.0 + self.glow;
        [
            (self.color[0] * scale).min(1.0),
            (self.color[1] * scale).min(1.0),
            (self.color[2] * scale).min(1.0),
            alpha,
        ]
    }
}
/// Volumetric fog parameters.
#[derive(Debug, Clone, Copy)]
pub struct VolumetricFogParams {
    /// Fog density (extinction coefficient, m⁻¹).
    pub density: f64,
    /// Scattering albedo \[0, 1\].
    pub scatter_albedo: f64,
    /// Fog colour (RGB).
    pub color: [f64; 3],
    /// Height falloff (1/m): fog is denser at lower altitudes.
    pub height_falloff: f64,
    /// Base height (fog is at maximum below this).
    pub base_height: f64,
    /// Phase function anisotropy g (Henyey-Greenstein, -1..1).
    pub phase_g: f64,
}
impl VolumetricFogParams {
    /// Uniform thin fog.
    pub fn thin_fog() -> Self {
        Self {
            density: 0.01,
            scatter_albedo: 0.9,
            color: [0.85, 0.88, 0.92],
            height_falloff: 0.5,
            base_height: 0.0,
            phase_g: 0.0,
        }
    }
    /// Dense ground fog.
    pub fn ground_fog() -> Self {
        Self {
            density: 0.1,
            scatter_albedo: 0.95,
            color: [0.9, 0.92, 0.95],
            height_falloff: 2.0,
            base_height: 0.0,
            phase_g: 0.2,
        }
    }
    /// Beer-Lambert transmittance over distance d at height y.
    pub fn transmittance(&self, distance: f64, height: f64) -> f64 {
        let h_factor = if self.height_falloff > 0.0 {
            (-self.height_falloff * (height - self.base_height).max(0.0)).exp()
        } else {
            1.0
        };
        let sigma_e = self.density * h_factor;
        (-sigma_e * distance).exp()
    }
    /// Henyey-Greenstein phase function.
    pub fn phase(&self, cos_theta: f64) -> f64 {
        let g = self.phase_g;
        let denom = 1.0 + g * g - 2.0 * g * cos_theta;
        (1.0 - g * g) / (4.0 * PI * denom.powf(1.5))
    }
    /// In-scattered radiance (single-scattering estimate).
    ///
    /// `n_steps` = ray march steps, `ray_length` = total march distance.
    pub fn inscatter(
        &self,
        ray_origin: [f64; 3],
        ray_dir: [f64; 3],
        ray_length: f64,
        light_dir: [f64; 3],
        light_color: [f64; 3],
        n_steps: usize,
    ) -> [f64; 3] {
        if n_steps == 0 {
            return [0.0; 3];
        }
        let step = ray_length / n_steps as f64;
        let cos_theta = dot3(ray_dir, light_dir);
        let ph = self.phase(cos_theta);
        let mut accum = [0.0f64; 3];
        let mut transmittance_acc = 1.0f64;
        for i in 0..n_steps {
            let t = (i as f64 + 0.5) * step;
            let pos = [
                ray_origin[0] + ray_dir[0] * t,
                ray_origin[1] + ray_dir[1] * t,
                ray_origin[2] + ray_dir[2] * t,
            ];
            let height = pos[1];
            let h_factor = if self.height_falloff > 0.0 {
                (-self.height_falloff * (height - self.base_height).max(0.0)).exp()
            } else {
                1.0
            };
            let sigma_e = self.density * h_factor;
            let sigma_s = sigma_e * self.scatter_albedo;
            let tr_step = (-sigma_e * step).exp();
            for ch in 0..3 {
                accum[ch] +=
                    transmittance_acc * sigma_s * ph * light_color[ch] * self.color[ch] * step;
            }
            transmittance_acc *= tr_step;
        }
        accum
    }
    /// Fog blend: lerp(scene_color, fog_color, 1 - transmittance).
    pub fn blend_scene(&self, scene_color: [f64; 3], distance: f64, height: f64) -> [f64; 3] {
        let t = self.transmittance(distance, height);
        std::array::from_fn(|i| t * scene_color[i] + (1.0 - t) * self.color[i])
    }
}
/// Source code and metadata for a single shader.
pub struct ShaderSource {
    /// Pipeline stage this shader targets.
    pub stage: ShaderStage,
    /// GLSL/WGSL/SPIR-V source text.
    pub source: String,
    /// Name of the entry-point function.
    pub entry_point: String,
}
impl ShaderSource {
    /// Create a vertex shader with the conventional `main` entry point.
    pub fn vertex(source: &str) -> Self {
        Self {
            stage: ShaderStage::Vertex,
            source: source.to_owned(),
            entry_point: "main".to_owned(),
        }
    }
    /// Create a fragment shader with the conventional `main` entry point.
    pub fn fragment(source: &str) -> Self {
        Self {
            stage: ShaderStage::Fragment,
            source: source.to_owned(),
            entry_point: "main".to_owned(),
        }
    }
    /// Create a compute shader with the conventional `main` entry point.
    pub fn compute(source: &str) -> Self {
        Self {
            stage: ShaderStage::Compute,
            source: source.to_owned(),
            entry_point: "main".to_owned(),
        }
    }
}
/// Phong shading parameters.
pub struct PhongMaterial {
    /// Ambient reflectance coefficient (RGB).
    pub ambient: [f64; 3],
    /// Diffuse reflectance coefficient (RGB).
    pub diffuse: [f64; 3],
    /// Specular reflectance coefficient (RGB).
    pub specular: [f64; 3],
    /// Specular shininess exponent.
    pub shininess: f64,
}
impl PhongMaterial {
    /// A default mid-gray Phong material.
    pub fn default_gray() -> Self {
        Self {
            ambient: [0.1, 0.1, 0.1],
            diffuse: [0.6, 0.6, 0.6],
            specular: [0.3, 0.3, 0.3],
            shininess: 32.0,
        }
    }
    /// Convert to a generic [`Material`] carrying Phong shader stubs and
    /// the material parameters as uniforms.
    pub fn to_material(&self) -> Material {
        const VS: &str = r#"
// Phong vertex shader (stub)
void main() {
    gl_Position = u_mvp * a_position;
    v_normal = mat3(u_normal_matrix) * a_normal;
    v_world_pos = (u_model * a_position).xyz;
}
"#;
        const FS: &str = r#"
// Phong fragment shader (stub)
void main() {
    vec3 n = normalize(v_normal);
    vec3 l = normalize(u_light_dir);
    float ndotl = max(dot(n, l), 0.0);
    vec3 r = reflect(-l, n);
    vec3 v = normalize(u_view_dir);
    float spec = pow(max(dot(r, v), 0.0), u_shininess);
    vec3 color = u_ambient + u_diffuse * ndotl + u_specular * spec;
    frag_color = vec4(color, 1.0);
}
"#;
        let mut mat = Material::new(
            "phong",
            ShaderSource::vertex(VS),
            ShaderSource::fragment(FS),
        );
        mat.set_uniform("u_ambient", UniformValue::Vec3(self.ambient));
        mat.set_uniform("u_diffuse", UniformValue::Vec3(self.diffuse));
        mat.set_uniform("u_specular", UniformValue::Vec3(self.specular));
        mat.set_uniform("u_shininess", UniformValue::Float(self.shininess));
        mat
    }
    /// Compute Phong irradiance at a surface point.
    pub fn phong_lighting(
        mat: &PhongMaterial,
        normal: [f64; 3],
        light_dir: [f64; 3],
        view_dir: [f64; 3],
    ) -> [f64; 3] {
        let n = normalize3(normal);
        let l = normalize3(light_dir);
        let v = normalize3(view_dir);
        let n_dot_l = dot3(n, l).max(0.0);
        let r = [
            2.0 * n_dot_l * n[0] - l[0],
            2.0 * n_dot_l * n[1] - l[1],
            2.0 * n_dot_l * n[2] - l[2],
        ];
        let r_dot_v = dot3(r, v).max(0.0);
        let spec = r_dot_v.powf(mat.shininess);
        [
            mat.ambient[0] + mat.diffuse[0] * n_dot_l + mat.specular[0] * spec,
            mat.ambient[1] + mat.diffuse[1] * n_dot_l + mat.specular[1] * spec,
            mat.ambient[2] + mat.diffuse[2] * n_dot_l + mat.specular[2] * spec,
        ]
    }
}
/// A named render material that pairs vertex/fragment shaders with uniforms.
pub struct Material {
    /// Human-readable name.
    pub name: String,
    /// Vertex shader source.
    pub vertex_shader: ShaderSource,
    /// Fragment shader source.
    pub fragment_shader: ShaderSource,
    /// Per-material uniform values keyed by name.
    pub uniforms: HashMap<String, UniformValue>,
}
impl Material {
    /// Construct a new material from the given shaders with an empty uniform map.
    pub fn new(name: &str, vs: ShaderSource, fs: ShaderSource) -> Self {
        Self {
            name: name.to_owned(),
            vertex_shader: vs,
            fragment_shader: fs,
            uniforms: HashMap::new(),
        }
    }
    /// Insert or replace a uniform value.
    pub fn set_uniform(&mut self, name: &str, value: UniformValue) {
        self.uniforms.insert(name.to_owned(), value);
    }
    /// Look up a uniform value by name.
    pub fn get_uniform(&self, name: &str) -> Option<&UniformValue> {
        self.uniforms.get(name)
    }
    /// Number of uniforms.
    pub fn uniform_count(&self) -> usize {
        self.uniforms.len()
    }
}
/// Which pipeline stage a shader belongs to.
pub enum ShaderStage {
    /// Vertex processing stage.
    Vertex,
    /// Fragment/pixel shading stage.
    Fragment,
    /// General-purpose compute stage.
    Compute,
}
/// A typed uniform value that can be uploaded to a shader.
pub enum UniformValue {
    /// Single 64-bit float.
    Float(f64),
    /// Three-component vector.
    Vec3([f64; 3]),
    /// Four-component vector.
    Vec4([f64; 4]),
    /// 4x4 column-major matrix.
    Mat4([[f64; 4]; 4]),
    /// 32-bit signed integer.
    Int(i32),
    /// Boolean flag.
    Bool(bool),
    /// Two-component vector.
    Vec2([f64; 2]),
}
/// Parameters for a subsurface scattering (SSS) material.
///
/// Uses the dipole diffusion approximation for skin/wax/marble.
#[derive(Debug, Clone)]
pub struct SssParams {
    /// Scatter mean free path length (per RGB channel, in world units).
    pub mfp: [f64; 3],
    /// Albedo at each wavelength (fraction scattered vs absorbed).
    pub albedo: [f64; 3],
    /// Index of refraction (n ≈ 1.4 for skin).
    pub ior: f64,
    /// Subsurface tint colour (multiplied into diffuse term).
    pub tint: [f64; 3],
    /// Blend weight between SSS and standard diffuse.
    pub sss_weight: f64,
}
impl SssParams {
    /// Realistic skin parameters.
    pub fn skin() -> Self {
        Self {
            mfp: [3.67e-3, 1.37e-3, 6.8e-4],
            albedo: [0.44, 0.22, 0.13],
            ior: 1.40,
            tint: [1.0, 0.85, 0.70],
            sss_weight: 0.8,
        }
    }
    /// Marble parameters.
    pub fn marble() -> Self {
        Self {
            mfp: [2.19e-3, 2.62e-3, 3.00e-3],
            albedo: [0.83, 0.79, 0.75],
            ior: 1.55,
            tint: [0.95, 0.93, 0.90],
            sss_weight: 0.5,
        }
    }
    /// Fresnel reflectance (Schlick approximation, scalar).
    pub fn fresnel_schlick_scalar(&self, cos_theta: f64) -> f64 {
        let r0 = ((self.ior - 1.0) / (self.ior + 1.0)).powi(2);
        r0 + (1.0 - r0) * (1.0 - cos_theta).powi(5)
    }
    /// Dipole diffuse reflectance R(r) at surface distance r.
    /// Simplified single-channel formula: R(r) = A / (4π) * (e^(-σtr*r)/r^2).
    pub fn dipole_reflectance(&self, r: f64, channel: usize) -> f64 {
        if r < 1e-12 {
            return self.albedo[channel.min(2)];
        }
        let mfp = self.mfp[channel.min(2)].max(1e-12);
        let sigma_tr = (3.0 / (mfp * mfp)).sqrt();
        let a = self.albedo[channel.min(2)];
        a / (4.0 * PI) * (-sigma_tr * r).exp() / (r * r)
    }
    /// Evaluate SSS lighting for a surface point.
    ///
    /// `n_dot_l` is the dot product of surface normal and light direction.
    /// `r` is the distance from the shading point to the light sample.
    pub fn evaluate_sss(&self, n_dot_l: f64, r: f64, light_color: [f64; 3]) -> [f64; 3] {
        let diffuse_weight = n_dot_l.max(0.0);
        std::array::from_fn(|ch| {
            let rr = self.dipole_reflectance(r, ch);
            let standard = self.tint[ch] * diffuse_weight * light_color[ch];
            let sss = rr * PI * light_color[ch] * self.tint[ch];
            (1.0 - self.sss_weight) * standard + self.sss_weight * sss
        })
    }
}
/// Physically-based rendering (PBR) material parameters (metallic-roughness workflow).
pub struct PbrMaterial {
    /// Albedo / base colour (linear RGB).
    pub base_color: [f64; 3],
    /// Metallic factor in \[0, 1\].
    pub metallic: f64,
    /// Perceptual roughness in \[0, 1\].
    pub roughness: f64,
    /// Emissive radiance (linear RGB).
    pub emissive: [f64; 3],
    /// Ambient occlusion factor in \[0, 1\].
    pub ao: f64,
}

impl Default for PbrMaterial {
    /// Sensible default: opaque mid-gray dielectric, roughness 0.5.
    fn default() -> Self {
        Self {
            base_color: [0.5, 0.5, 0.5],
            metallic: 0.0,
            roughness: 0.5,
            emissive: [0.0, 0.0, 0.0],
            ao: 1.0,
        }
    }
}

impl PbrMaterial {
    /// Create a metallic material.
    pub fn metallic(base_color: [f64; 3], roughness: f64) -> Self {
        Self {
            base_color,
            metallic: 1.0,
            roughness,
            emissive: [0.0; 3],
            ao: 1.0,
        }
    }
    /// Create a dielectric (non-metallic) material.
    pub fn dielectric(base_color: [f64; 3], roughness: f64) -> Self {
        Self {
            base_color,
            metallic: 0.0,
            roughness,
            emissive: [0.0; 3],
            ao: 1.0,
        }
    }
    /// Compute F0 (reflectance at normal incidence).
    ///
    /// For dielectrics, F0 is typically 0.04.
    /// For metals, F0 is the base color.
    pub fn f0(&self) -> [f64; 3] {
        let dielectric_f0 = 0.04;
        [
            dielectric_f0 * (1.0 - self.metallic) + self.base_color[0] * self.metallic,
            dielectric_f0 * (1.0 - self.metallic) + self.base_color[1] * self.metallic,
            dielectric_f0 * (1.0 - self.metallic) + self.base_color[2] * self.metallic,
        ]
    }
    /// GGX (Trowbridge-Reitz) normal distribution function.
    ///
    /// `D(h) = alpha^2 / (pi * ((n.h)^2 * (alpha^2 - 1) + 1)^2)`
    pub fn ggx_distribution(n_dot_h: f64, roughness: f64) -> f64 {
        let alpha = roughness * roughness;
        let alpha2 = alpha * alpha;
        let denom = n_dot_h * n_dot_h * (alpha2 - 1.0) + 1.0;
        alpha2 / (PI * denom * denom)
    }
    /// Schlick approximation to the Fresnel reflectance.
    ///
    /// `F = F0 + (1 - F0) * (1 - cos_theta)^5`
    pub fn schlick_fresnel(cos_theta: f64, f0: [f64; 3]) -> [f64; 3] {
        let t = (1.0 - cos_theta).powi(5);
        [
            f0[0] + (1.0 - f0[0]) * t,
            f0[1] + (1.0 - f0[1]) * t,
            f0[2] + (1.0 - f0[2]) * t,
        ]
    }
    /// Schlick-GGX geometry function for a single direction.
    ///
    /// `G1(n, v) = (n.v) / ((n.v) * (1 - k) + k)`
    ///
    /// where `k = (roughness + 1)^2 / 8` for direct lighting.
    pub fn geometry_schlick_ggx(n_dot_v: f64, roughness: f64) -> f64 {
        let r = roughness + 1.0;
        let k = r * r / 8.0;
        let denom = n_dot_v * (1.0 - k) + k;
        if denom.abs() < 1e-30 {
            return 0.0;
        }
        n_dot_v / denom
    }
    /// Smith's geometry function combining view and light directions.
    ///
    /// `G = G1(n, v) * G1(n, l)`
    pub fn geometry_smith(n_dot_v: f64, n_dot_l: f64, roughness: f64) -> f64 {
        Self::geometry_schlick_ggx(n_dot_v.max(0.0), roughness)
            * Self::geometry_schlick_ggx(n_dot_l.max(0.0), roughness)
    }
    /// Convert PBR parameters to a [`Material`] with uniforms.
    pub fn to_material(&self) -> Material {
        let vs = ShaderSource::vertex("void main() { /* PBR VS stub */ }");
        let fs = ShaderSource::fragment("void main() { /* PBR FS stub */ }");
        let mut mat = Material::new("pbr", vs, fs);
        mat.set_uniform("u_base_color", UniformValue::Vec3(self.base_color));
        mat.set_uniform("u_metallic", UniformValue::Float(self.metallic));
        mat.set_uniform("u_roughness", UniformValue::Float(self.roughness));
        mat.set_uniform("u_emissive", UniformValue::Vec3(self.emissive));
        mat.set_uniform("u_ao", UniformValue::Float(self.ao));
        mat
    }
}
/// A morph target: displacement vectors per vertex.
#[derive(Debug, Clone)]
pub struct MorphTarget {
    /// Name for this target (e.g. "smile", "blink_left").
    pub name: String,
    /// Displacement vectors per vertex.
    pub displacements: Vec<[f32; 3]>,
}
impl MorphTarget {
    /// Create a morph target.
    pub fn new(name: &str, displacements: Vec<[f32; 3]>) -> Self {
        Self {
            name: name.to_owned(),
            displacements,
        }
    }
    /// Number of vertices in this morph target.
    pub fn vertex_count(&self) -> usize {
        self.displacements.len()
    }
}
/// Anisotropic specular BRDF parameters (Ashikhmin-Shirley model).
#[derive(Debug, Clone, Copy)]
pub struct AnisotropicBrdf {
    /// Roughness along tangent direction.
    pub alpha_x: f64,
    /// Roughness along bitangent direction.
    pub alpha_y: f64,
    /// Specular reflectance at normal incidence (F0).
    pub f0: [f64; 3],
}
impl AnisotropicBrdf {
    /// Isotropic version (alpha_x = alpha_y).
    pub fn isotropic(roughness: f64, f0: [f64; 3]) -> Self {
        Self {
            alpha_x: roughness,
            alpha_y: roughness,
            f0,
        }
    }
    /// Fresnel (Schlick) per channel.
    pub fn fresnel(&self, cos_theta: f64) -> [f64; 3] {
        self.f0.map(|f| f + (1.0 - f) * (1.0 - cos_theta).powi(5))
    }
    /// GGX NDF for anisotropic surfaces.
    ///
    /// h = half vector, t = tangent, b = bitangent, n = normal.
    pub fn ndf_ggx_aniso(
        &self,
        h: [f64; 3],
        tangent: [f64; 3],
        bitangent: [f64; 3],
        normal: [f64; 3],
    ) -> f64 {
        let ax = self.alpha_x.max(1e-4);
        let ay = self.alpha_y.max(1e-4);
        let h_dot_t = dot3(h, tangent);
        let h_dot_b = dot3(h, bitangent);
        let h_dot_n = dot3(h, normal).max(0.0);
        let denom = (h_dot_t / ax).powi(2) + (h_dot_b / ay).powi(2) + h_dot_n * h_dot_n;
        1.0 / (PI * ax * ay * denom * denom)
    }
    /// Evaluate anisotropic BRDF value.
    ///
    /// Returns RGB specular contribution.
    pub fn evaluate(
        &self,
        view_dir: [f64; 3],
        light_dir: [f64; 3],
        normal: [f64; 3],
        tangent: [f64; 3],
        bitangent: [f64; 3],
    ) -> [f64; 3] {
        let h = normalize3([
            view_dir[0] + light_dir[0],
            view_dir[1] + light_dir[1],
            view_dir[2] + light_dir[2],
        ]);
        let n_dot_l = dot3(normal, light_dir).max(0.0);
        let n_dot_v = dot3(normal, view_dir).max(1e-4);
        let v_dot_h = dot3(view_dir, h).max(0.0);
        let d = self.ndf_ggx_aniso(h, tangent, bitangent, normal);
        let fresnel = self.fresnel(v_dot_h);
        let g = (4.0_f64 * n_dot_l * n_dot_v).max(1e-6);
        fresnel.map(|f| d * f / g * n_dot_l)
    }
    /// Hemispherical directional reflectance (approximate energy check).
    pub fn directional_albedo(&self, n_dot_v: f64) -> [f64; 3] {
        self.fresnel(n_dot_v.max(0.0))
    }
}
/// Thin-film iridescence parameters.
#[derive(Debug, Clone, Copy)]
pub struct IridescenceParams {
    /// Film thickness in nanometres.
    pub thickness_nm: f64,
    /// Index of refraction of the thin film.
    pub ior_film: f64,
    /// Index of refraction of the substrate.
    pub ior_substrate: f64,
    /// Iridescence strength \[0, 1\].
    pub strength: f64,
}
impl IridescenceParams {
    /// Typical soap bubble.
    pub fn soap_bubble() -> Self {
        Self {
            thickness_nm: 300.0,
            ior_film: 1.34,
            ior_substrate: 1.0,
            strength: 1.0,
        }
    }
    /// Oil film on water.
    pub fn oil_film() -> Self {
        Self {
            thickness_nm: 200.0,
            ior_film: 1.47,
            ior_substrate: 1.33,
            strength: 0.8,
        }
    }
    /// Optical path difference for a given angle of incidence.
    pub fn opd(&self, cos_theta_film: f64) -> f64 {
        2.0 * self.ior_film * self.thickness_nm * cos_theta_film
    }
    /// Evaluate iridescent color for given wavelengths.
    ///
    /// `wavelengths_nm` should be \[R, G, B\] wavelengths (e.g. \[700, 546, 435\]).
    pub fn evaluate_rgb(&self, cos_theta: f64, wavelengths_nm: [f64; 3]) -> [f64; 3] {
        let sin_t = (1.0 - cos_theta * cos_theta).sqrt() / self.ior_film;
        let sin_t_clamped = sin_t.clamp(-1.0, 1.0);
        let cos_theta_film = (1.0 - sin_t_clamped * sin_t_clamped).sqrt().max(0.0);
        wavelengths_nm.map(|lambda| {
            let opd = self.opd(cos_theta_film);
            let phase = 2.0 * PI * opd / lambda.max(1.0);
            let interference = ((phase.cos() + 1.0) * 0.5).powi(2);
            self.strength * interference + (1.0 - self.strength)
        })
    }
    /// Fresnel reflectance at film interface (simplified, unpolarised).
    pub fn fresnel_film(&self, cos_theta: f64) -> f64 {
        let n1 = 1.0;
        let n2 = self.ior_film;
        let r0 = ((n1 - n2) / (n1 + n2)).powi(2);
        r0 + (1.0 - r0) * (1.0 - cos_theta).powi(5)
    }
}
