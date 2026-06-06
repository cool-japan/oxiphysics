//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)
use super::functions::*;
use std::collections::HashMap;

/// A capacity-bounded LRU cache for compiled shader variants.
///
/// When the cache exceeds `capacity` entries the least-recently-accessed
/// variant is evicted.
pub struct VariantCache {
    /// Maximum number of variants to hold.
    pub capacity: usize,
    /// Stored variants in insertion order (we use swap-remove for eviction).
    pub(super) entries: Vec<(ShaderKey, CompiledVariant, u64)>,
    /// Monotonic access counter.
    pub(super) clock: u64,
}
impl VariantCache {
    /// Create a new empty cache with the given capacity.
    pub fn new(capacity: usize) -> Self {
        Self {
            capacity: capacity.max(1),
            entries: Vec::new(),
            clock: 0,
        }
    }
    /// Insert a compiled variant.  Evicts the LRU entry if at capacity.
    pub fn insert(&mut self, variant: CompiledVariant) {
        if let Some(pos) = self.entries.iter().position(|(k, _, _)| *k == variant.key) {
            self.clock += 1;
            self.entries[pos] = (variant.key.clone(), variant, self.clock);
            return;
        }
        if self.entries.len() >= self.capacity {
            let lru_pos = self
                .entries
                .iter()
                .enumerate()
                .min_by_key(|(_, (_, _, t))| *t)
                .map(|(i, _)| i)
                .unwrap_or(0);
            self.entries.swap_remove(lru_pos);
        }
        self.clock += 1;
        self.entries
            .push((variant.key.clone(), variant, self.clock));
    }
    /// Look up a variant by key.  Bumps its access time on hit.
    pub fn get(&mut self, key: &ShaderKey) -> Option<&CompiledVariant> {
        if let Some(pos) = self.entries.iter().position(|(k, _, _)| k == key) {
            self.clock += 1;
            self.entries[pos].2 = self.clock;
            return Some(&self.entries[pos].1);
        }
        None
    }
    /// Number of variants currently cached.
    pub fn len(&self) -> usize {
        self.entries.len()
    }
    /// True when the cache is empty.
    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }
    /// Remove all entries.
    pub fn clear(&mut self) {
        self.entries.clear();
        self.clock = 0;
    }
    /// Total binary bytes cached.
    pub fn total_binary_bytes(&self) -> usize {
        self.entries
            .iter()
            .map(|(_, v, _)| v.binary_size_bytes)
            .sum()
    }
}
/// A compiled, instantiated shader variant ready for dispatch.
#[derive(Debug, Clone)]
pub struct CompiledVariant {
    /// The key that identifies this variant.
    pub key: ShaderKey,
    /// The fully-substituted WGSL source text.
    pub wgsl: String,
    /// The effective workgroup size after substitution.
    pub workgroup_size: [u32; 3],
    /// Byte size of the compiled binary mock (always `wgsl.len()` bytes here).
    pub binary_size_bytes: usize,
}
impl CompiledVariant {
    fn new(key: ShaderKey, wgsl: String, workgroup_size: [u32; 3]) -> Self {
        let binary_size_bytes = wgsl.len();
        Self {
            key,
            wgsl,
            workgroup_size,
            binary_size_bytes,
        }
    }
}
/// A registry of named [`VariantProfile`]s with dependency resolution.
#[derive(Debug, Default)]
pub struct VariantProfileRegistry {
    pub(super) profiles: HashMap<String, VariantProfile>,
}
impl VariantProfileRegistry {
    /// Create an empty registry.
    pub fn new() -> Self {
        Self::default()
    }
    /// Register a profile.  Overwrites any existing profile with the same name.
    pub fn register(&mut self, profile: VariantProfile) {
        self.profiles.insert(profile.name.clone(), profile);
    }
    /// Resolve a profile by name, merging inherited defines.
    ///
    /// Returns `None` if the profile (or any base) is not found.
    pub fn resolve(&self, name: &str) -> Option<HashMap<String, String>> {
        let profile = self.profiles.get(name)?;
        let mut merged = if let Some(base) = &profile.base {
            self.resolve(base)?
        } else {
            HashMap::new()
        };
        for (k, v) in &profile.defines {
            merged.insert(k.clone(), v.clone());
        }
        Some(merged)
    }
    /// Return a sorted list of all registered profile names.
    pub fn profile_names(&self) -> Vec<&str> {
        let mut names: Vec<&str> = self.profiles.keys().map(|s| s.as_str()).collect();
        names.sort_unstable();
        names
    }
}
/// Tracks which shaders depend on (include) which other shaders.
///
/// When a shader is modified, all shaders that depend on it (transitively)
/// also need to be recompiled.
#[derive(Debug, Default)]
pub struct ShaderDependencyGraph {
    /// `depends_on[A]` = set of shaders A includes / depends on.
    pub(super) depends_on: HashMap<String, Vec<String>>,
}
impl ShaderDependencyGraph {
    /// Create an empty dependency graph.
    pub fn new() -> Self {
        Self::default()
    }
    /// Record that shader `dependent` depends on `dependency`.
    pub fn add_dependency(&mut self, dependent: impl Into<String>, dependency: impl Into<String>) {
        self.depends_on
            .entry(dependent.into())
            .or_default()
            .push(dependency.into());
    }
    /// Return all shaders that directly depend on `dependency`.
    pub fn direct_dependents(&self, dependency: &str) -> Vec<&str> {
        self.depends_on
            .iter()
            .filter(|(_, deps)| deps.iter().any(|d| d == dependency))
            .map(|(name, _)| name.as_str())
            .collect()
    }
    /// Return all shaders that transitively depend on `changed_shader`
    /// (BFS / DFS over the dependency graph).
    pub fn transitive_dependents(&self, changed_shader: &str) -> Vec<String> {
        let mut visited = std::collections::HashSet::new();
        let mut queue = std::collections::VecDeque::new();
        queue.push_back(changed_shader.to_string());
        while let Some(current) = queue.pop_front() {
            for (name, deps) in &self.depends_on {
                if deps.contains(&current) && !visited.contains(name.as_str()) {
                    visited.insert(name.clone());
                    queue.push_back(name.clone());
                }
            }
        }
        let mut result: Vec<String> = visited.into_iter().collect();
        result.sort_unstable();
        result
    }
    /// Return the direct dependencies of `shader`.
    pub fn direct_dependencies(&self, shader: &str) -> &[String] {
        self.depends_on
            .get(shader)
            .map(Vec::as_slice)
            .unwrap_or(&[])
    }
    /// Return all shader names that have any recorded dependencies.
    pub fn shaders_with_deps(&self) -> Vec<&str> {
        let mut names: Vec<&str> = self.depends_on.keys().map(|s| s.as_str()).collect();
        names.sort_unstable();
        names
    }
}
/// Uniquely identifies a compiled shader variant.
///
/// Two keys are equal when the base shader name and all defines match.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct ShaderKey {
    /// The base shader name (e.g. `"sph_density"`).
    pub name: String,
    /// Sorted define pairs `(KEY, VALUE)` baked into this variant.
    pub defines: Vec<(String, String)>,
}
impl ShaderKey {
    /// Create a key from a name and an unsorted define map.
    pub fn new(name: impl Into<String>, defines: &HashMap<String, String>) -> Self {
        let mut sorted: Vec<(String, String)> = defines
            .iter()
            .map(|(k, v)| (k.clone(), v.clone()))
            .collect();
        sorted.sort_by(|a, b| a.0.cmp(&b.0));
        Self {
            name: name.into(),
            defines: sorted,
        }
    }
    /// Create a key with no defines.
    pub fn bare(name: impl Into<String>) -> Self {
        Self {
            name: name.into(),
            defines: Vec::new(),
        }
    }
    /// Return a deterministic string fingerprint suitable for file names.
    pub fn fingerprint(&self) -> String {
        if self.defines.is_empty() {
            return self.name.clone();
        }
        let suffix: String = self
            .defines
            .iter()
            .map(|(k, v)| format!("{k}_{v}"))
            .collect::<Vec<_>>()
            .join("__");
        format!("{}__{}", self.name, suffix)
    }
}
/// A named set of compile-time defines forming a "variant profile".
///
/// Variant profiles can inherit from a base profile, inheriting its defines
/// while optionally overriding individual values.
#[derive(Debug, Clone)]
pub struct VariantProfile {
    /// Unique name for this profile.
    pub name: String,
    /// The defines in this profile.
    pub defines: HashMap<String, String>,
    /// Optional base profile name (inherited defines are merged first).
    pub base: Option<String>,
}
impl VariantProfile {
    /// Create a new variant profile with no base.
    pub fn new(name: impl Into<String>) -> Self {
        Self {
            name: name.into(),
            defines: HashMap::new(),
            base: None,
        }
    }
    /// Create a variant profile that inherits from `base`.
    pub fn with_base(name: impl Into<String>, base: impl Into<String>) -> Self {
        Self {
            name: name.into(),
            defines: HashMap::new(),
            base: Some(base.into()),
        }
    }
    /// Set a define on this profile.
    pub fn set(mut self, key: impl Into<String>, value: impl Into<String>) -> Self {
        self.defines.insert(key.into(), value.into());
        self
    }
}
/// Tracks shader source modification times to support hot reload.
///
/// In production this would watch file system events; here timestamps
/// are set manually for testing.
#[derive(Debug, Default)]
pub struct HotReloadTracker {
    /// `name → last_modified` (monotonic counter, not wall clock).
    pub(super) timestamps: HashMap<String, u64>,
    /// `name → last_compiled` counter.
    pub(super) compiled_at: HashMap<String, u64>,
}
impl HotReloadTracker {
    /// Create a new tracker.
    pub fn new() -> Self {
        Self::default()
    }
    /// Record a modification for `name` at the given logical time.
    pub fn touch(&mut self, name: impl Into<String>, time: u64) {
        self.timestamps.insert(name.into(), time);
    }
    /// Record that `name` was compiled at the given logical time.
    pub fn record_compile(&mut self, name: impl Into<String>, time: u64) {
        self.compiled_at.insert(name.into(), time);
    }
    /// Returns `true` when `name` has been modified since it was last compiled.
    pub fn needs_recompile(&self, name: &str) -> bool {
        let modified = self.timestamps.get(name).copied().unwrap_or(0);
        let compiled = self.compiled_at.get(name).copied().unwrap_or(0);
        modified > compiled
    }
    /// Return all shader names that need recompilation.
    pub fn stale_shaders(&self) -> Vec<&str> {
        self.timestamps
            .keys()
            .filter(|n| self.needs_recompile(n))
            .map(|s| s.as_str())
            .collect()
    }
}
impl HotReloadTracker {
    /// Touch multiple shaders at once with the same logical timestamp.
    pub fn touch_batch(&mut self, names: &[&str], time: u64) {
        for &name in names {
            self.touch(name, time);
        }
    }
    /// Record that all currently stale shaders have been recompiled at `time`.
    pub fn flush_stale(&mut self, time: u64) {
        let stale: Vec<String> = self
            .timestamps
            .keys()
            .filter(|n| self.needs_recompile(n))
            .cloned()
            .collect();
        for name in stale {
            self.compiled_at.insert(name, time);
        }
    }
    /// Return the set of shader names that have been modified but never compiled.
    pub fn never_compiled(&self) -> Vec<&str> {
        self.timestamps
            .keys()
            .filter(|n| !self.compiled_at.contains_key(n.as_str()))
            .map(|s| s.as_str())
            .collect()
    }
}
/// Central registry for shader sources and compiled variants.
///
/// Call [`ShaderRegistry::register`] to add a [`ShaderSource`], then
/// [`ShaderRegistry::get_or_compile`] to obtain a compiled variant,
/// which is transparently cached in the internal [`VariantCache`].
pub struct ShaderRegistry {
    /// All registered shader sources, keyed by name.
    pub(super) sources: HashMap<String, ShaderSource>,
    /// Compiled variant cache.
    pub(super) cache: VariantCache,
    /// Number of cache hits since creation.
    pub cache_hits: u64,
    /// Number of compilations triggered since creation.
    pub compilations: u64,
}
impl ShaderRegistry {
    /// Create a new empty registry with the given cache capacity.
    pub fn new(cache_capacity: usize) -> Self {
        Self {
            sources: HashMap::new(),
            cache: VariantCache::new(cache_capacity),
            cache_hits: 0,
            compilations: 0,
        }
    }
    /// Register a shader source.  Overwrites any existing source with the
    /// same name and clears cached variants for that name.
    pub fn register(&mut self, source: ShaderSource) {
        let name = source.name.clone();
        self.sources.insert(name.clone(), source);
        self.cache.entries.retain(|(k, _, _)| k.name != name);
    }
    /// Get or compile a variant.
    ///
    /// Returns `Err` when the source name is unknown or when a required
    /// placeholder is not supplied in `defines`.
    pub fn get_or_compile(
        &mut self,
        name: &str,
        defines: &HashMap<String, String>,
    ) -> Result<&CompiledVariant, RegistryError> {
        let key = ShaderKey::new(name, defines);
        if self.cache.get(&key).is_some() {
            self.cache_hits += 1;
            return Ok(self
                .cache
                .get(&key)
                .expect("cache entry must exist after insertion"));
        }
        let source = self
            .sources
            .get(name)
            .ok_or_else(|| RegistryError::UnknownShader(name.to_string()))?;
        for ph in &source.placeholders {
            if !defines.contains_key(ph.as_str()) {
                return Err(RegistryError::MissingDefine {
                    shader: name.to_string(),
                    define: ph.clone(),
                });
            }
        }
        let wgsl = source.instantiate(defines);
        let workgroup_size = source.workgroup_size;
        let variant = CompiledVariant::new(key.clone(), wgsl, workgroup_size);
        self.compilations += 1;
        self.cache.insert(variant);
        Ok(self
            .cache
            .get(&key)
            .expect("cache entry must exist after insertion"))
    }
    /// Return all registered shader names.
    pub fn shader_names(&self) -> Vec<&str> {
        self.sources.keys().map(|s| s.as_str()).collect()
    }
    /// Return the number of entries currently in the variant cache.
    pub fn cached_count(&self) -> usize {
        self.cache.len()
    }
    /// Invalidate a single shader by name (clears its cached variants).
    pub fn invalidate(&mut self, name: &str) {
        self.cache.entries.retain(|(k, _, _)| k.name != name);
    }
    /// Invalidate all cached variants.
    pub fn invalidate_all(&mut self) {
        self.cache.clear();
    }
}
impl ShaderRegistry {
    /// Compile a specific shader variant with the given options.
    ///
    /// Unlike [`get_or_compile`](ShaderRegistry::get_or_compile) this method
    /// always forces a fresh compilation (not served from the cache) and
    /// respects the `max_source_bytes` limit in `opts`.
    pub fn compile_variant(
        &mut self,
        name: &str,
        defines: &HashMap<String, String>,
        opts: &ShaderCompileOptions,
    ) -> Result<CompiledVariant, RegistryError> {
        let source = self
            .sources
            .get(name)
            .ok_or_else(|| RegistryError::UnknownShader(name.to_string()))?;
        let mut merged = defines.clone();
        for (k, v) in &opts.extra_defines {
            merged.entry(k.clone()).or_insert_with(|| v.clone());
        }
        for ph in &source.placeholders {
            if !merged.contains_key(ph.as_str()) {
                return Err(RegistryError::MissingDefine {
                    shader: name.to_string(),
                    define: ph.clone(),
                });
            }
        }
        let wgsl = source.instantiate(&merged);
        if opts.max_source_bytes > 0 && wgsl.len() > opts.max_source_bytes {
            return Err(RegistryError::SourceTooLarge {
                shader: name.to_string(),
                size: wgsl.len(),
                limit: opts.max_source_bytes,
            });
        }
        let workgroup_size = opts.workgroup_size;
        let key = ShaderKey::new(name, &merged);
        self.compilations += 1;
        Ok(CompiledVariant::new(key, wgsl, workgroup_size))
    }
}
impl ShaderRegistry {
    /// Check a [`HotReloadTracker`] and invalidate any stale shader variants.
    ///
    /// Returns the list of shader names that were invalidated.
    pub fn apply_hot_reload(&mut self, tracker: &HotReloadTracker) -> Vec<String> {
        let stale: Vec<String> = tracker
            .stale_shaders()
            .into_iter()
            .map(|s| s.to_string())
            .collect();
        for name in &stale {
            self.invalidate(name);
        }
        stale
    }
    /// Return all shader names registered in this registry.
    pub fn registered_count(&self) -> usize {
        self.sources.len()
    }
    /// Return the size (in WGSL source bytes) of the named shader, or 0.
    pub fn source_bytes(&self, name: &str) -> usize {
        self.sources.get(name).map(|s| s.wgsl.len()).unwrap_or(0)
    }
}
/// A raw WGSL shader source with associated metadata.
#[derive(Debug, Clone)]
pub struct ShaderSource {
    /// Human-readable identifier (must be unique in the registry).
    pub name: String,
    /// Raw WGSL text (may contain `{{PLACEHOLDER}}` tokens).
    pub wgsl: String,
    /// Suggested workgroup size `[x, y, z]`.
    pub workgroup_size: [u32; 3],
    /// List of placeholder names this source accepts.
    pub placeholders: Vec<String>,
}
impl ShaderSource {
    /// Create a new shader source.
    pub fn new(name: impl Into<String>, wgsl: impl Into<String>, workgroup_size: [u32; 3]) -> Self {
        let wgsl_str: String = wgsl.into();
        let placeholders = collect_placeholders(&wgsl_str);
        Self {
            name: name.into(),
            wgsl: wgsl_str,
            workgroup_size,
            placeholders,
        }
    }
    /// Instantiate this source by substituting `{{KEY}}` tokens with values.
    ///
    /// Unknown keys are left untouched.  Returns the instantiated WGSL string.
    pub fn instantiate(&self, defines: &HashMap<String, String>) -> String {
        let mut out = self.wgsl.clone();
        for (k, v) in defines {
            let token = format!("{{{{{}}}}}", k);
            out = out.replace(&token, v);
        }
        out
    }
    /// Return the total thread count per workgroup.
    pub fn threads_per_group(&self) -> u32 {
        self.workgroup_size[0] * self.workgroup_size[1] * self.workgroup_size[2]
    }
}
/// A typed specialization constant (analogous to Vulkan spec constants).
#[derive(Debug, Clone, PartialEq)]
pub enum SpecConstValue {
    /// Integer constant.
    Int(i64),
    /// Unsigned integer constant.
    Uint(u64),
    /// Floating-point constant.
    Float(f64),
    /// Boolean constant.
    Bool(bool),
}
impl SpecConstValue {
    /// Render the value as a WGSL literal string.
    pub fn to_wgsl(&self) -> String {
        match self {
            SpecConstValue::Int(v) => v.to_string(),
            SpecConstValue::Uint(v) => format!("{v}u"),
            SpecConstValue::Float(v) => format!("{v}"),
            SpecConstValue::Bool(v) => v.to_string(),
        }
    }
}
/// A named specialization constant bound to a shader.
#[derive(Debug, Clone, PartialEq)]
pub struct SpecializationConstant {
    /// Name used in the WGSL source as `{{SPEC_NAME}}`.
    pub name: String,
    /// Default value.
    pub default_value: SpecConstValue,
    /// Optional override value set at pipeline creation time.
    pub override_value: Option<SpecConstValue>,
}
impl SpecializationConstant {
    /// Create a new specialization constant with a default value.
    pub fn new(name: impl Into<String>, default: SpecConstValue) -> Self {
        Self {
            name: name.into(),
            default_value: default,
            override_value: None,
        }
    }
    /// Set an override value.
    pub fn with_override(mut self, value: SpecConstValue) -> Self {
        self.override_value = Some(value);
        self
    }
    /// Return the effective value (override if set, else default).
    pub fn effective_value(&self) -> &SpecConstValue {
        self.override_value.as_ref().unwrap_or(&self.default_value)
    }
}
/// A simple pipeline cache mapping [`PipelineCacheKey`]s to byte blobs
/// (mock: stores a label string for each entry).
#[derive(Debug, Default)]
pub struct PipelineCache {
    pub(super) entries: HashMap<PipelineCacheKey, String>,
    /// Number of cache hits.
    pub hits: u64,
    /// Number of cache misses.
    pub misses: u64,
}
impl PipelineCache {
    /// Create an empty pipeline cache.
    pub fn new() -> Self {
        Self::default()
    }
    /// Insert a pipeline entry.
    pub fn insert(&mut self, key: PipelineCacheKey, label: impl Into<String>) {
        self.entries.insert(key, label.into());
    }
    /// Look up a pipeline.  Returns `Some(&label)` on hit, `None` on miss.
    pub fn get(&mut self, key: &PipelineCacheKey) -> Option<&str> {
        if let Some(v) = self.entries.get(key) {
            self.hits += 1;
            Some(v.as_str())
        } else {
            self.misses += 1;
            None
        }
    }
    /// Number of entries in the cache.
    pub fn len(&self) -> usize {
        self.entries.len()
    }
    /// Returns `true` when the cache is empty.
    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }
    /// Clear all entries.
    pub fn clear(&mut self) {
        self.entries.clear();
    }
    /// Hit rate in `[0.0, 1.0]`.
    pub fn hit_rate(&self) -> f64 {
        let total = self.hits + self.misses;
        if total == 0 {
            return 0.0;
        }
        self.hits as f64 / total as f64
    }
}
/// Describes a GPU compute pipeline built from a shader variant.
#[derive(Debug, Clone)]
pub struct PipelineDescriptor {
    /// The shader variant key used by this pipeline.
    pub key: ShaderKey,
    /// Bind group layouts (number of bind groups required).
    pub bind_group_count: u32,
    /// Push-constant block size in bytes (0 = none).
    pub push_constant_bytes: u32,
    /// Human-readable label for debugging.
    pub label: String,
}
impl PipelineDescriptor {
    /// Create a new pipeline descriptor.
    pub fn new(
        key: ShaderKey,
        bind_group_count: u32,
        push_constant_bytes: u32,
        label: impl Into<String>,
    ) -> Self {
        Self {
            key,
            bind_group_count,
            push_constant_bytes,
            label: label.into(),
        }
    }
    /// Validate the descriptor (mock: checks workgroup product).
    pub fn validate(&self) -> Result<(), RegistryError> {
        if self.bind_group_count == 0 {
            return Err(RegistryError::UnknownShader(
                "pipeline has no bind groups".into(),
            ));
        }
        Ok(())
    }
}
/// Options controlling how a shader variant is compiled.
#[derive(Debug, Clone, PartialEq)]
pub struct ShaderCompileOptions {
    /// Workgroup sizes to try when specialising the shader (x, y, z).
    pub workgroup_size: [u32; 3],
    /// Additional `#define`-style substitutions beyond the registry defaults.
    pub extra_defines: HashMap<String, String>,
    /// Optional byte-size budget.  Compilation fails if the instantiated WGSL
    /// exceeds this limit (0 = unlimited).
    pub max_source_bytes: usize,
}
impl ShaderCompileOptions {
    /// Construct default options.
    pub fn new() -> Self {
        Self {
            workgroup_size: [64, 1, 1],
            extra_defines: HashMap::new(),
            max_source_bytes: 0,
        }
    }
}
/// A set of specialization constants for a shader.
#[derive(Debug, Clone, Default)]
pub struct SpecConstSet {
    /// The constants in this set, keyed by name.
    pub constants: HashMap<String, SpecializationConstant>,
}
impl SpecConstSet {
    /// Create an empty set.
    pub fn new() -> Self {
        Self::default()
    }
    /// Add a constant.
    pub fn add(&mut self, constant: SpecializationConstant) {
        self.constants.insert(constant.name.clone(), constant);
    }
    /// Build a `HashMap<String, String>` suitable for passing to
    /// [`ShaderSource::instantiate`].
    pub fn to_defines(&self) -> HashMap<String, String> {
        self.constants
            .iter()
            .map(|(name, c)| (name.clone(), c.effective_value().to_wgsl()))
            .collect()
    }
    /// Return `true` if a constant named `name` is in this set.
    pub fn has(&self, name: &str) -> bool {
        self.constants.contains_key(name)
    }
    /// Return the effective WGSL string for constant `name`, or `None`.
    pub fn get_wgsl(&self, name: &str) -> Option<String> {
        self.constants
            .get(name)
            .map(|c| c.effective_value().to_wgsl())
    }
}
/// Errors returned by the shader registry.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum RegistryError {
    /// No source registered under the given name.
    UnknownShader(String),
    /// A required `{{DEFINE}}` was not supplied.
    MissingDefine {
        /// Shader that requires the missing define.
        shader: String,
        /// The missing define name.
        define: String,
    },
    /// The instantiated WGSL exceeds the permitted byte limit.
    SourceTooLarge {
        /// Shader that is too large.
        shader: String,
        /// Actual size in bytes.
        size: usize,
        /// Permitted byte limit.
        limit: usize,
    },
}
/// Composite key for a GPU pipeline cache entry.
///
/// Uniquely identifies a compiled pipeline based on the shader variant and
/// pipeline parameters.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct PipelineCacheKey {
    /// Shader variant key.
    pub shader_key: ShaderKey,
    /// Workgroup size.
    pub workgroup_size: [u32; 3],
    /// Push constant size in bytes.
    pub push_constant_bytes: u32,
    /// Bind group layout signature (e.g. a hash of the layout).
    pub layout_hash: u64,
}
impl PipelineCacheKey {
    /// Create a new cache key.
    pub fn new(
        shader_key: ShaderKey,
        workgroup_size: [u32; 3],
        push_constant_bytes: u32,
        layout_hash: u64,
    ) -> Self {
        Self {
            shader_key,
            workgroup_size,
            push_constant_bytes,
            layout_hash,
        }
    }
    /// Compute an FNV-1a hash of the entire key.
    pub fn hash_key(&self) -> u64 {
        let repr = format!(
            "{}_{:?}_{}_{}",
            self.shader_key.fingerprint(),
            self.workgroup_size,
            self.push_constant_bytes,
            self.layout_hash,
        );
        compute_cache_key(&repr)
    }
}
