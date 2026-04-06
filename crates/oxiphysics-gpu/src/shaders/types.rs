//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions_2::*;
use std::collections::HashMap;

#[allow(unused_imports)]
use super::functions::*;
use super::functions::{
    BOUNDARY_ENFORCE_WGSL, BROADPHASE_SORT_SHADER, INTEGRATE_WGSL, LBM_BGK_D2Q9_WGSL,
    LBM_STREAMING_SHADER, RIGID_INTEGRATE_SHADER, SPH_DENSITY_WGSL, SPH_FORCE_WGSL,
};

/// Reflection data extracted from a (mock) SPIR-V module or WGSL source.
#[derive(Debug, Clone)]
pub struct SpirVModule {
    /// Entry point function names.
    pub entry_points: Vec<String>,
    /// Number of bindings.
    pub binding_count: usize,
    /// Workgroup size \[x, y, z\].
    pub workgroup_size: [u32; 3],
    /// Raw SPIR-V bytes.
    pub spirv_bytes: Vec<u8>,
}
impl SpirVModule {
    /// Build a SpirV module from WGSL source (mock reflection).
    pub fn from_wgsl(source: &str) -> Self {
        let mut entry_points = Vec::new();
        for line in source.lines() {
            let trimmed = line.trim();
            if let Some(pos) = trimmed.find("fn ") {
                let rest = &trimmed[pos + 3..];
                let name: String = rest
                    .chars()
                    .take_while(|c| c.is_alphanumeric() || *c == '_')
                    .collect();
                if !name.is_empty() {
                    entry_points.push(name);
                }
            }
        }
        let binding_count = source.matches("@binding(").count();
        let workgroup_size = parse_workgroup_size(source);
        let spirv_bytes = mock_compile_to_spirv(
            source,
            entry_points.first().map(|s| s.as_str()).unwrap_or("main"),
        );
        Self {
            entry_points,
            binding_count,
            workgroup_size,
            spirv_bytes,
        }
    }
    /// Check whether this module has a specific entry point.
    pub fn has_entry_point(&self, name: &str) -> bool {
        self.entry_points.iter().any(|e| e == name)
    }
    /// Size of the SPIR-V binary in bytes.
    pub fn byte_size(&self) -> usize {
        self.spirv_bytes.len()
    }
}
/// A push constant range that can be set per draw/dispatch call.
#[derive(Debug, Clone, Copy)]
pub struct PushConstantRange {
    /// Byte offset into the push constant block.
    pub offset: u32,
    /// Size in bytes.
    pub size: u32,
    /// Shader stage.
    pub stage: ShaderStage,
}
impl PushConstantRange {
    /// Create a new push constant range.
    pub fn new(offset: u32, size: u32, stage: ShaderStage) -> Self {
        Self {
            offset,
            size,
            stage,
        }
    }
    /// Check if this range fits within the standard 128-byte GPU limit.
    pub fn fits_standard_limit(&self) -> bool {
        self.offset + self.size <= 128
    }
}
/// Texture address mode.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AddressMode {
    /// Clamp coordinates to \[0, 1\].
    ClampToEdge,
    /// Wrap coordinates (repeat).
    Repeat,
    /// Mirror-repeat coordinates.
    MirrorRepeat,
}
/// A layout describing all bindings in a descriptor set group.
#[derive(Debug, Clone)]
pub struct DescriptorSetLayout {
    /// Group index (set number).
    pub group: u32,
    /// All binding entries.
    pub bindings: Vec<DescriptorBinding>,
}
impl DescriptorSetLayout {
    /// Create a new empty layout for the given group.
    pub fn new(group: u32) -> Self {
        Self {
            group,
            bindings: Vec::new(),
        }
    }
    /// Add a storage buffer binding.
    pub fn add_storage_buffer(&mut self, binding: u32, stage: ShaderStage, read_only: bool) {
        self.bindings.push(DescriptorBinding {
            binding,
            descriptor_type: DescriptorType::StorageBuffer,
            stage,
            read_only,
        });
    }
    /// Add a uniform buffer binding.
    pub fn add_uniform_buffer(&mut self, binding: u32, stage: ShaderStage) {
        self.bindings.push(DescriptorBinding {
            binding,
            descriptor_type: DescriptorType::UniformBuffer,
            stage,
            read_only: true,
        });
    }
    /// Add a combined image sampler binding.
    pub fn add_sampler(&mut self, binding: u32, stage: ShaderStage) {
        self.bindings.push(DescriptorBinding {
            binding,
            descriptor_type: DescriptorType::CombinedImageSampler,
            stage,
            read_only: true,
        });
    }
    /// Number of bindings in this layout.
    pub fn len(&self) -> usize {
        self.bindings.len()
    }
    /// Whether the layout has no bindings.
    pub fn is_empty(&self) -> bool {
        self.bindings.is_empty()
    }
}
/// A storage binding entry in a bind group.
#[derive(Debug, Clone)]
pub struct StorageBinding {
    /// Name in WGSL.
    pub name: String,
    /// Binding index.
    pub binding: u32,
    /// Whether read-only.
    pub read_only: bool,
}
/// Manages shader source files and detects changes for hot-reloading.
#[derive(Debug, Default)]
pub struct ShaderHotReloadManager {
    /// Map from shader name to current source.
    pub(super) sources: HashMap<String, String>,
    /// Map from shader name to source hash (for change detection).
    pub(super) hashes: HashMap<String, u64>,
}
impl ShaderHotReloadManager {
    /// Create a new hot-reload manager.
    pub fn new() -> Self {
        Self::default()
    }
    /// Start watching a shader file with initial source.
    pub fn watch(&mut self, name: &str, source: &str) {
        let hash = simple_hash(source);
        self.sources.insert(name.to_string(), source.to_string());
        self.hashes.insert(name.to_string(), hash);
    }
    /// Stop watching a shader.
    pub fn unwatch(&mut self, name: &str) {
        self.sources.remove(name);
        self.hashes.remove(name);
    }
    /// Update a shader's source. Returns `true` if the source changed.
    pub fn update(&mut self, name: &str, new_source: &str) -> bool {
        let new_hash = simple_hash(new_source);
        if let Some(old_hash) = self.hashes.get(name)
            && *old_hash == new_hash
        {
            return false;
        }
        self.sources
            .insert(name.to_string(), new_source.to_string());
        self.hashes.insert(name.to_string(), new_hash);
        true
    }
    /// Check if a shader is being watched.
    pub fn is_watched(&self, name: &str) -> bool {
        self.sources.contains_key(name)
    }
    /// Get the current source for a watched shader.
    pub fn get_source(&self, name: &str) -> Option<&str> {
        self.sources.get(name).map(|s| s.as_str())
    }
    /// Return names of all watched shaders.
    pub fn watched_names(&self) -> Vec<&str> {
        self.sources.keys().map(|s| s.as_str()).collect()
    }
    /// Number of watched shaders.
    pub fn len(&self) -> usize {
        self.sources.len()
    }
    /// Whether no shaders are watched.
    pub fn is_empty(&self) -> bool {
        self.sources.is_empty()
    }
}
/// Type of a GPU descriptor.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DescriptorType {
    /// Uniform buffer.
    UniformBuffer,
    /// Storage buffer.
    StorageBuffer,
    /// Combined image sampler.
    CombinedImageSampler,
    /// Storage image.
    StorageImage,
}
/// Registry mapping shader names to their [`ShaderMetadata`].
#[derive(Debug, Default)]
pub struct ShaderMetaRegistry {
    pub(super) entries: HashMap<String, ShaderMetadata>,
}
impl ShaderMetaRegistry {
    /// Create an empty registry.
    pub fn new() -> Self {
        Self::default()
    }
    /// Register a shader under `name`.
    pub fn register(&mut self, name: &str, meta: ShaderMetadata) {
        self.entries.insert(name.to_string(), meta);
    }
    /// Look up a shader by name.
    pub fn lookup(&self, name: &str) -> Option<&ShaderMetadata> {
        self.entries.get(name)
    }
    /// Return all registered shader names.
    pub fn all_names(&self) -> Vec<&str> {
        self.entries.keys().map(|s| s.as_str()).collect()
    }
    /// Number of registered shaders.
    pub fn len(&self) -> usize {
        self.entries.len()
    }
    /// True if no shaders are registered.
    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }
}
/// A specialization constant for shader compilation.
#[derive(Debug, Clone)]
pub struct SpecializationConstant {
    /// Name of the constant in the shader source.
    pub name: String,
    /// Default value as a string.
    pub default_value: String,
    /// Description of what the constant controls.
    pub description: String,
}
impl SpecializationConstant {
    /// Create a new specialization constant.
    pub fn new(name: &str, default_value: &str, description: &str) -> Self {
        Self {
            name: name.to_string(),
            default_value: default_value.to_string(),
            description: description.to_string(),
        }
    }
}
/// A parameterised WGSL shader template.
///
/// Placeholders of the form `{KEY}` in the template string are replaced by
/// the corresponding values supplied via [`ShaderTemplate::instantiate`].
#[derive(Debug, Clone)]
pub struct ShaderTemplate {
    /// The raw template source with `{KEY}` placeholders.
    pub template: String,
}
impl ShaderTemplate {
    /// Create a new shader template from a source string.
    pub fn new(template: impl Into<String>) -> Self {
        Self {
            template: template.into(),
        }
    }
    /// Instantiate the template by replacing all `{KEY}` placeholders.
    ///
    /// # Arguments
    /// * `params` – map from placeholder name (without braces) to replacement string.
    ///
    /// # Returns
    /// A new `String` with all recognised placeholders replaced.  Unknown
    /// placeholders are left as-is.
    pub fn instantiate(&self, params: &HashMap<&str, &str>) -> String {
        let mut result = self.template.clone();
        for (key, value) in params {
            let placeholder = format!("{{{}}}", key);
            result = result.replace(&placeholder, value);
        }
        result
    }
    /// Return a list of placeholder names found in the template.
    ///
    /// Only returns names that look like valid template placeholders:
    /// `{IDENTIFIER}` where IDENTIFIER is uppercase letters, digits, and
    /// underscores, starting with an uppercase letter.
    pub fn placeholders(&self) -> Vec<String> {
        let mut result = Vec::new();
        let chars: Vec<char> = self.template.chars().collect();
        let mut i = 0;
        while i < chars.len() {
            if chars[i] == '{' {
                if i + 1 < chars.len() && chars[i + 1].is_ascii_uppercase() {
                    let start = i + 1;
                    let mut end = start;
                    while end < chars.len() && chars[end] != '}' {
                        end += 1;
                    }
                    if end < chars.len() {
                        let name: String = chars[start..end].iter().collect();
                        let is_valid = name
                            .chars()
                            .all(|c| c.is_ascii_uppercase() || c.is_ascii_digit() || c == '_');
                        if is_valid && !result.contains(&name) {
                            result.push(name);
                        }
                        i = end + 1;
                    } else {
                        i += 1;
                    }
                } else {
                    i += 1;
                }
            } else {
                i += 1;
            }
        }
        result
    }
    /// Check if all placeholders are provided in the params map.
    pub fn all_placeholders_provided(&self, params: &HashMap<&str, &str>) -> bool {
        for p in self.placeholders() {
            if !params.contains_key(p.as_str()) {
                return false;
            }
        }
        true
    }
}
/// Texture format for attachments.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TextureFormat {
    /// 8-bit RGBA unorm.
    Rgba8Unorm,
    /// 8-bit RGBA sRGB.
    Rgba8Srgb,
    /// 16-bit float RGBA (HDR).
    Rgba16Float,
    /// 32-bit float single channel.
    R32Float,
    /// 32-bit depth.
    Depth32Float,
    /// 24-bit depth + 8-bit stencil.
    Depth24PlusStencil8,
}
/// A flexible bind group layout builder.
#[derive(Debug, Clone, Default)]
pub struct BindGroupLayout {
    pub(super) uniforms: Vec<UniformBinding>,
    pub(super) storages: Vec<StorageBinding>,
}
impl BindGroupLayout {
    /// Create an empty layout.
    pub fn new() -> Self {
        Self::default()
    }
    /// Add a uniform buffer binding.
    pub fn add_uniform(&mut self, name: &str, _group: u32, binding: u32, size_bytes: u32) {
        self.uniforms.push(UniformBinding {
            name: name.to_string(),
            binding,
            size_bytes,
        });
    }
    /// Add a storage buffer binding.
    pub fn add_storage(&mut self, name: &str, _group: u32, binding: u32, read_only: bool) {
        self.storages.push(StorageBinding {
            name: name.to_string(),
            binding,
            read_only,
        });
    }
    /// Total number of bindings.
    pub fn binding_count(&self) -> usize {
        self.uniforms.len() + self.storages.len()
    }
    /// Whether this layout has no bindings.
    pub fn is_empty(&self) -> bool {
        self.uniforms.is_empty() && self.storages.is_empty()
    }
    /// Generate a WGSL snippet declaring all bindings.
    pub fn to_wgsl_snippet(&self) -> String {
        let mut out = String::new();
        for u in &self.uniforms {
            out.push_str(&format!(
                "@group(0) @binding({}) var<uniform> {}: {};\n",
                u.binding,
                u.name.to_lowercase(),
                u.name
            ));
        }
        for s in &self.storages {
            let access = if s.read_only { "read" } else { "read_write" };
            out.push_str(&format!(
                "@group(0) @binding({}) var<storage, {}> {}: array<f32>;\n",
                s.binding,
                access,
                s.name.to_lowercase()
            ));
        }
        out
    }
}
/// Description of a depth/stencil attachment.
#[derive(Debug, Clone)]
pub struct DepthAttachmentDesc {
    /// Format of the depth attachment.
    pub format: TextureFormat,
    /// Load operation for depth.
    pub load_op: LoadOp,
    /// Store operation for depth.
    pub store_op: StoreOp,
    /// Clear depth value.
    pub clear_depth: f32,
}
/// A pipeline for compiling WGSL shaders with preprocessing steps.
///
/// Steps: resolve includes -> apply specialization -> validate -> cache.
pub struct ShaderCompilationPipeline {
    pub(super) includes: HashMap<String, String>,
    pub(super) cache: ShaderCache,
}
impl ShaderCompilationPipeline {
    /// Create a new compilation pipeline.
    pub fn new() -> Self {
        Self {
            includes: HashMap::new(),
            cache: ShaderCache::new(),
        }
    }
    /// Register an include file for `#include` resolution.
    pub fn add_include(&mut self, name: &str, source: &str) {
        self.includes.insert(name.to_string(), source.to_string());
    }
    /// Compile a shader source through the pipeline.
    ///
    /// Returns the processed source string. Caches the result.
    pub fn compile(
        &mut self,
        name: &str,
        source: &str,
        spec_map: Option<&SpecializationMap>,
    ) -> Result<String, String> {
        if let Some(cached) = self.cache.entries.get(name) {
            return Ok(cached.clone());
        }
        let includes_ref: HashMap<&str, &str> = self
            .includes
            .iter()
            .map(|(k, v)| (k.as_str(), v.as_str()))
            .collect();
        let resolved = resolve_includes(source, &includes_ref);
        let specialized = if let Some(sm) = spec_map {
            sm.apply(&resolved)
        } else {
            resolved
        };
        if !validate_wgsl_structure(&specialized) {
            return Err(format!("shader '{}' failed structural validation", name));
        }
        self.cache
            .entries
            .insert(name.to_string(), specialized.clone());
        Ok(specialized)
    }
    /// Return the number of cached shaders.
    pub fn cache_size(&self) -> usize {
        self.cache.len()
    }
    /// Clear the compilation cache.
    pub fn clear_cache(&mut self) {
        self.cache.clear();
    }
}
/// A registry of named compute shader descriptors.
#[derive(Debug, Default)]
pub struct ShaderRegistry {
    pub(super) shaders: HashMap<String, ComputeShaderDesc>,
}
impl ShaderRegistry {
    /// Create a new empty registry.
    pub fn new() -> Self {
        Self::default()
    }
    /// Register a shader under the given name.
    pub fn register(&mut self, name: impl Into<String>, desc: ComputeShaderDesc) {
        self.shaders.insert(name.into(), desc);
    }
    /// Retrieve a shader descriptor by name.
    pub fn get(&self, name: &str) -> Option<&ComputeShaderDesc> {
        self.shaders.get(name)
    }
    /// Return the number of registered shaders.
    pub fn len(&self) -> usize {
        self.shaders.len()
    }
    /// Return true if the registry is empty.
    pub fn is_empty(&self) -> bool {
        self.shaders.is_empty()
    }
    /// Return an iterator over registered shader names.
    pub fn names(&self) -> impl Iterator<Item = &str> {
        self.shaders.keys().map(|s| s.as_str())
    }
    /// Remove a shader from the registry.
    pub fn unregister(&mut self, name: &str) -> Option<ComputeShaderDesc> {
        self.shaders.remove(name)
    }
    /// Check if a shader is registered.
    pub fn contains(&self, name: &str) -> bool {
        self.shaders.contains_key(name)
    }
    /// Create a registry pre-populated with the built-in shaders.
    pub fn with_builtins() -> Self {
        let mut reg = Self::new();
        reg.register(
            "sph_density",
            ComputeShaderDesc::new("main", [64, 1, 1], SPH_DENSITY_WGSL),
        );
        reg.register(
            "sph_force",
            ComputeShaderDesc::new("main", [64, 1, 1], SPH_FORCE_WGSL),
        );
        reg.register(
            "integrate",
            ComputeShaderDesc::new("main", [64, 1, 1], INTEGRATE_WGSL),
        );
        reg.register(
            "lbm_bgk_d2q9",
            ComputeShaderDesc::new("main", [64, 1, 1], LBM_BGK_D2Q9_WGSL),
        );
        reg.register(
            "lbm_streaming",
            ComputeShaderDesc::new("main", [64, 1, 1], LBM_STREAMING_SHADER),
        );
        reg.register(
            "rigid_integrate",
            ComputeShaderDesc::new("main", [64, 1, 1], RIGID_INTEGRATE_SHADER),
        );
        reg.register(
            "broadphase_sort",
            ComputeShaderDesc::new("main", [64, 1, 1], BROADPHASE_SORT_SHADER),
        );
        reg.register(
            "boundary_enforce",
            ComputeShaderDesc::new("main", [64, 1, 1], BOUNDARY_ENFORCE_WGSL),
        );
        reg
    }
}
/// Descriptor for a WGSL compute shader.
#[derive(Debug, Clone)]
pub struct ComputeShaderDesc {
    /// Name of the entry-point function (e.g. `"main"`).
    pub entry_point: String,
    /// Workgroup size `[x, y, z]`.
    pub workgroup_size: [u32; 3],
    /// Full WGSL source string.
    pub source: String,
}
impl ComputeShaderDesc {
    /// Create a new compute shader descriptor.
    pub fn new(
        entry_point: impl Into<String>,
        workgroup_size: [u32; 3],
        source: impl Into<String>,
    ) -> Self {
        Self {
            entry_point: entry_point.into(),
            workgroup_size,
            source: source.into(),
        }
    }
    /// Number of threads per workgroup.
    pub fn threads_per_workgroup(&self) -> u32 {
        self.workgroup_size[0] * self.workgroup_size[1] * self.workgroup_size[2]
    }
    /// Count the number of binding annotations in the shader source.
    pub fn binding_count(&self) -> usize {
        self.source.matches("@binding(").count()
    }
}
/// Sampler descriptor for creating GPU samplers.
#[derive(Debug, Clone)]
pub struct SamplerDesc {
    /// Minification filter.
    pub filter_min: FilterMode,
    /// Magnification filter.
    pub filter_mag: FilterMode,
    /// Address mode for all axes.
    pub address_mode: AddressMode,
    /// Anisotropy level (1 = disabled).
    pub anisotropy: u32,
    /// LOD bias.
    pub lod_bias: f32,
    /// Maximum LOD level.
    pub lod_max: f32,
}
impl SamplerDesc {
    /// Create a linear sampler with clamp-to-edge.
    pub fn linear() -> Self {
        Self {
            filter_min: FilterMode::Linear,
            filter_mag: FilterMode::Linear,
            address_mode: AddressMode::ClampToEdge,
            anisotropy: 1,
            lod_bias: 0.0,
            lod_max: 16.0,
        }
    }
    /// Create a nearest-neighbor sampler.
    pub fn nearest() -> Self {
        Self {
            filter_min: FilterMode::Nearest,
            filter_mag: FilterMode::Nearest,
            address_mode: AddressMode::ClampToEdge,
            anisotropy: 1,
            lod_bias: 0.0,
            lod_max: 0.0,
        }
    }
    /// Create an anisotropic sampler with repeat addressing.
    pub fn anisotropic(max_anisotropy: u32) -> Self {
        Self {
            filter_min: FilterMode::Linear,
            filter_mag: FilterMode::Linear,
            address_mode: AddressMode::Repeat,
            anisotropy: max_anisotropy,
            lod_bias: 0.0,
            lod_max: 16.0,
        }
    }
}
/// Description of a full render pass.
#[derive(Debug, Clone)]
pub struct RenderPassDesc {
    /// Color attachments.
    pub color_attachments: Vec<ColorAttachmentDesc>,
    /// Optional depth/stencil attachment.
    pub depth_attachment: Option<DepthAttachmentDesc>,
    /// Name of this render pass (for debugging).
    pub name: String,
}
impl RenderPassDesc {
    /// Create a simple render pass with one RGBA8 color attachment, no depth.
    pub fn new_simple_color() -> Self {
        Self {
            color_attachments: vec![ColorAttachmentDesc {
                format: TextureFormat::Rgba8Unorm,
                load_op: LoadOp::Clear,
                store_op: StoreOp::Store,
                clear_color: [0.0, 0.0, 0.0, 1.0],
            }],
            depth_attachment: None,
            name: "SimpleColor".to_string(),
        }
    }
    /// Create a render pass with an HDR color attachment and depth buffer.
    pub fn new_with_depth() -> Self {
        Self {
            color_attachments: vec![ColorAttachmentDesc {
                format: TextureFormat::Rgba16Float,
                load_op: LoadOp::Clear,
                store_op: StoreOp::Store,
                clear_color: [0.0, 0.0, 0.0, 1.0],
            }],
            depth_attachment: Some(DepthAttachmentDesc {
                format: TextureFormat::Depth32Float,
                load_op: LoadOp::Clear,
                store_op: StoreOp::Store,
                clear_depth: 1.0,
            }),
            name: "ColorDepth".to_string(),
        }
    }
    /// Number of total attachments (color + optional depth).
    pub fn total_attachment_count(&self) -> usize {
        self.color_attachments.len()
            + if self.depth_attachment.is_some() {
                1
            } else {
                0
            }
    }
}
/// High-level category for a GPU compute shader.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub enum ShaderVariant {
    /// General physics integration / force kernels.
    Physics,
    /// Collision detection and response.
    Collision,
    /// Smoothed particle hydrodynamics.
    Sph,
    /// Lattice Boltzmann method.
    Lbm,
    /// Rigid body dynamics.
    RigidBody,
    /// Neural network inference.
    NeuralInference,
}
/// Metadata describing a single compute shader.
#[derive(Debug, Clone)]
pub struct ShaderMetadata {
    /// The high-level shader category.
    pub variant: ShaderVariant,
    /// Name of the entry-point function (e.g. `"main"`).
    pub entry_point: String,
    /// Workgroup size `[x, y, z]`.
    pub workgroup_size: [u32; 3],
    /// Number of bind groups required by the shader.
    pub bind_group_count: u32,
}
impl ShaderMetadata {
    /// Create a new metadata record.
    pub fn new(
        variant: ShaderVariant,
        entry_point: impl Into<String>,
        workgroup_size: [u32; 3],
        bind_group_count: u32,
    ) -> Self {
        Self {
            variant,
            entry_point: entry_point.into(),
            workgroup_size,
            bind_group_count,
        }
    }
    /// Total number of threads per workgroup.
    pub fn threads_per_workgroup(&self) -> u32 {
        self.workgroup_size[0] * self.workgroup_size[1] * self.workgroup_size[2]
    }
}
/// A simple cache for compiled shader sources.
///
/// Caches shader source strings by a composite key of name + specialization.
#[derive(Debug, Default)]
pub struct ShaderCache {
    pub(super) entries: HashMap<String, String>,
}
impl ShaderCache {
    /// Create a new empty cache.
    pub fn new() -> Self {
        Self::default()
    }
    /// Get a cached shader source by key, or compute and cache it.
    pub fn get_or_insert(&mut self, key: &str, compute: impl FnOnce() -> String) -> &str {
        self.entries.entry(key.to_string()).or_insert_with(compute)
    }
    /// Check if a shader is cached.
    pub fn contains(&self, key: &str) -> bool {
        self.entries.contains_key(key)
    }
    /// Return the number of cached entries.
    pub fn len(&self) -> usize {
        self.entries.len()
    }
    /// Check if the cache is empty.
    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }
    /// Clear the cache.
    pub fn clear(&mut self) {
        self.entries.clear();
    }
    /// Remove a specific entry.
    pub fn remove(&mut self, key: &str) -> Option<String> {
        self.entries.remove(key)
    }
}
/// Shader stage flags.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ShaderStage {
    /// Vertex shader stage.
    Vertex,
    /// Fragment shader stage.
    Fragment,
    /// Compute shader stage.
    Compute,
    /// All stages.
    All,
}
/// A parameterised shader template that replaces `#define KEY` tokens.
///
/// Instantiation replaces occurrences of each key found in `defines`
/// with the corresponding value.  The key must appear as a whole word
/// (surrounded by non-identifier characters or at start/end of line).
#[derive(Debug, Clone)]
pub struct ShaderTemplateV2 {
    /// Raw template source text.
    pub source: String,
    /// Map from token name to replacement value.
    pub defines: HashMap<String, String>,
}
impl ShaderTemplateV2 {
    /// Create a new template with the given source and defines.
    pub fn new(source: impl Into<String>, defines: HashMap<String, String>) -> Self {
        Self {
            source: source.into(),
            defines,
        }
    }
    /// Instantiate the template by performing all `#define` substitutions.
    ///
    /// Each entry in `defines` replaces all occurrences of the key
    /// in the source with the corresponding value.
    pub fn instantiate(&self) -> String {
        let mut result = self.source.clone();
        for (key, value) in &self.defines {
            result = result.replace(key.as_str(), value.as_str());
        }
        result
    }
}
/// Description of a color attachment in a render pass.
#[derive(Debug, Clone)]
pub struct ColorAttachmentDesc {
    /// Texture format of the attachment.
    pub format: TextureFormat,
    /// Load operation.
    pub load_op: LoadOp,
    /// Store operation.
    pub store_op: StoreOp,
    /// Clear color \[r, g, b, a\].
    pub clear_color: [f32; 4],
}
/// A uniform binding entry in a bind group.
#[derive(Debug, Clone)]
pub struct UniformBinding {
    /// Name in WGSL.
    pub name: String,
    /// Binding index.
    pub binding: u32,
    /// Size in bytes.
    pub size_bytes: u32,
}
/// Store operation for an attachment.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StoreOp {
    /// Store the rendered content.
    Store,
    /// Discard the content after rendering.
    Discard,
}
/// A single binding entry in a descriptor set layout.
#[derive(Debug, Clone)]
pub struct DescriptorBinding {
    /// Binding index.
    pub binding: u32,
    /// Type of the descriptor.
    pub descriptor_type: DescriptorType,
    /// Shader stage that uses this binding.
    pub stage: ShaderStage,
    /// Whether a storage buffer is read-only.
    pub read_only: bool,
}
/// Texture filtering mode.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FilterMode {
    /// Nearest neighbor (point) filtering.
    Nearest,
    /// Bilinear/trilinear filtering.
    Linear,
}
/// Load operation for an attachment.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LoadOp {
    /// Load the existing content.
    Load,
    /// Clear to a specified value.
    Clear,
    /// Content doesn't matter.
    DontCare,
}
/// A collection of specialization constants for a shader.
#[derive(Debug, Clone, Default)]
pub struct SpecializationMap {
    pub(super) constants: Vec<SpecializationConstant>,
    pub(super) overrides: HashMap<String, String>,
}
impl SpecializationMap {
    /// Create a new empty specialization map.
    pub fn new() -> Self {
        Self::default()
    }
    /// Define a specialization constant with a default value.
    pub fn define(&mut self, name: &str, default_value: &str, description: &str) {
        self.constants.push(SpecializationConstant::new(
            name,
            default_value,
            description,
        ));
    }
    /// Override a constant's value.
    pub fn set(&mut self, name: &str, value: &str) {
        self.overrides.insert(name.to_string(), value.to_string());
    }
    /// Get the effective value for a constant (override or default).
    pub fn get(&self, name: &str) -> Option<&str> {
        if let Some(v) = self.overrides.get(name) {
            return Some(v.as_str());
        }
        for c in &self.constants {
            if c.name == name {
                return Some(c.default_value.as_str());
            }
        }
        None
    }
    /// Return the number of defined constants.
    pub fn len(&self) -> usize {
        self.constants.len()
    }
    /// Check if there are no constants defined.
    pub fn is_empty(&self) -> bool {
        self.constants.is_empty()
    }
    /// Apply specialization constants to a shader source by replacing
    /// `const {NAME} = {DEFAULT};` patterns with overridden values.
    pub fn apply(&self, source: &str) -> String {
        let mut result = source.to_string();
        for c in &self.constants {
            let value = self
                .overrides
                .get(&c.name)
                .map(|s| s.as_str())
                .unwrap_or(&c.default_value);
            let old = format!("const {} = {};", c.name, c.default_value);
            let new = format!("const {} = {};", c.name, value);
            result = result.replace(&old, &new);
        }
        result
    }
}
/// Description of a uniform buffer binding.
#[derive(Debug, Clone)]
pub struct UniformBufferDesc {
    /// Name of the uniform block.
    pub name: String,
    /// Bind group index.
    pub group: u32,
    /// Binding slot within the group.
    pub binding: u32,
    /// Size of the uniform buffer in bytes.
    pub size_bytes: u32,
}
impl UniformBufferDesc {
    /// Create a new uniform buffer description.
    pub fn new(name: &str, group: u32, binding: u32, size_bytes: u32) -> Self {
        Self {
            name: name.to_string(),
            group,
            binding,
            size_bytes,
        }
    }
    /// Generate the WGSL binding annotation string.
    pub fn wgsl_annotation(&self) -> String {
        format!(
            "@group({}) @binding({}) var<uniform> {}: {};",
            self.group,
            self.binding,
            self.name.to_lowercase(),
            self.name
        )
    }
}
/// Cache for compiled shader bytecode (mock compiled blobs).
///
/// Supports a maximum total byte budget; when the budget is exceeded,
/// `evict_oldest` removes the oldest-inserted entry.
#[derive(Debug)]
pub struct BytecodeShaderCache {
    /// Map from shader name to compiled bytecode.
    pub cache: HashMap<String, Vec<u8>>,
    /// Ordered insertion keys for LRU-like eviction.
    pub(super) insertion_order: Vec<String>,
    /// Maximum total bytes allowed.
    pub max_size: usize,
}
impl BytecodeShaderCache {
    /// Create a new cache with the given byte budget.
    pub fn new(max_size: usize) -> Self {
        Self {
            cache: HashMap::new(),
            insertion_order: Vec::new(),
            max_size,
        }
    }
    /// Insert (or replace) a shader's compiled bytecode.
    ///
    /// If the total size would exceed `max_size`, `evict_oldest` is called
    /// before inserting.
    pub fn insert(&mut self, name: &str, bytecode: Vec<u8>) {
        if self.cache.contains_key(name) {
            self.insertion_order.retain(|k| k != name);
        }
        self.cache.insert(name.to_string(), bytecode);
        self.insertion_order.push(name.to_string());
        while self.total_bytes() > self.max_size && !self.insertion_order.is_empty() {
            self.evict_oldest();
        }
    }
    /// Retrieve compiled bytecode by shader name.
    pub fn get(&self, name: &str) -> Option<&Vec<u8>> {
        self.cache.get(name)
    }
    /// Evict the oldest cached entry.  No-op if cache is empty.
    pub fn evict_oldest(&mut self) {
        if let Some(oldest) = self.insertion_order.first().cloned() {
            self.insertion_order.remove(0);
            self.cache.remove(&oldest);
        }
    }
    /// Total bytes currently stored across all cached entries.
    pub fn total_bytes(&self) -> usize {
        self.cache.values().map(|v| v.len()).sum()
    }
    /// Number of cached entries.
    pub fn len(&self) -> usize {
        self.cache.len()
    }
    /// True if the cache is empty.
    pub fn is_empty(&self) -> bool {
        self.cache.is_empty()
    }
}
