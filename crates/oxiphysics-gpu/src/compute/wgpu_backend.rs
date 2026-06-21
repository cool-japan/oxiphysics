// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! WebGPU (wgpu) compute backend for the OxiPhysics GPU acceleration layer.
//!
//! This module provides [`WgpuBackend`] which implements `ComputeBackend` using
//! the `wgpu` crate for cross-platform GPU compute (Vulkan, Metal, DX12, WebGPU).
//!
//! ## Feature flag
//!
//! The stub [`WgpuBackend`] in this module is **always** compiled: it stores
//! buffers as CPU-side `Vec<f64>` shadows and never touches a GPU, so the crate
//! builds without the `wgpu` dependency.  The real on-device backend
//! `real::WgpuBackendReal` and its `wgpu` / `pollster` / `bytemuck` dependencies
//! are gated behind the `wgpu-backend` Cargo feature:
//!
//! ```toml
//! [dependencies]
//! oxiphysics-gpu = { features = ["wgpu-backend"] }
//! ```
//!
//! ## Enabling the dependency
//!
//! To activate the wgpu backend, add `wgpu` to the crate's `Cargo.toml`:
//!
//! ```toml
//! [features]
//! wgpu-backend = ["wgpu"]
//!
//! [dependencies]
//! wgpu = { version = "0.20", optional = true }
//! ```
//!
//! ## Architecture
//!
//! ```text
//!  WgpuBackend
//!   ├── wgpu::Device / wgpu::Queue          ← GPU device & command queue
//!   ├── Vec<WgpuBufferEntry>                 ← Registered GPU buffers
//!   │     ├── wgpu::Buffer (device memory)
//!   │     └── size, usage flags
//!   └── ShaderRegistry                       ← Compiled WGSL compute shaders
//!
//!  Compute pipeline:
//!    write_buffer → [upload via staging] → dispatch(kernel) → [readback via staging] → read_buffer
//! ```
//!
//! ## Usage (when feature is enabled)
//!
//! ```ignore
//! use oxiphysics_gpu::compute::wgpu_backend::WgpuBackend;
//! use oxiphysics_gpu::compute::ComputeBackend;
//!
//! let backend = WgpuBackend::new_async().await?;
//! let handle = backend.create_buffer(1024);
//! backend.write_buffer(handle, &vec![1.0_f64; 128]);
//! // ... dispatch kernel ...
//! let data = backend.read_buffer(handle);
//! ```

// ── BufferHandle (re-used from parent module) ─────────────────────────────────

/// Opaque handle to a GPU buffer allocated by a `ComputeBackend`.
///
/// This type mirrors the one in the parent `compute` module so that
/// [`WgpuBackend`] can implement the same `ComputeBackend` trait.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct WgpuBufferHandle(pub usize);

// ── WgpuDeviceInfo ────────────────────────────────────────────────────────────

/// Information about the GPU device selected by the wgpu adapter.
#[derive(Debug, Clone, Default)]
pub struct WgpuDeviceInfo {
    /// Human-readable device name (e.g. `"NVIDIA GeForce RTX 4090"`).
    pub name: String,
    /// Backend API in use: `"Vulkan"`, `"Metal"`, `"Dx12"`, `"WebGpu"`, or `"None"`.
    pub backend: String,
    /// Driver version string (if available).
    pub driver_version: String,
    /// Total VRAM in bytes (0 if not reported by the adapter).
    pub vram_bytes: u64,
    /// Whether the device supports 64-bit floating-point storage.
    pub supports_f64: bool,
    /// Maximum workgroup size (x, y, z).
    pub max_workgroup_size: [u32; 3],
}

// ── WgpuBackend ───────────────────────────────────────────────────────────────

/// WebGPU compute backend.
///
/// This struct is **always** a CPU-emulation backend: it stores buffers as
/// CPU-side `Vec<f64>` shadows and its [`dispatch`](Self::dispatch) is a no-op
/// (it does not execute WGSL).  [`try_new`](Self::try_new) therefore returns
/// `Err(WgpuInitError::NotAvailable)`; use [`new_stub`](Self::new_stub) to get
/// an instance for exercising the buffer API.
///
/// The real on-device GPU backend is `real::WgpuBackendReal`, compiled under the
/// `wgpu-backend` feature; it owns the actual `wgpu::Device` / `Queue` and runs
/// WGSL kernels via `dispatch_wgsl`.
#[derive(Debug)]
pub struct WgpuBackend {
    /// Device info (populated at initialisation).
    pub device_info: WgpuDeviceInfo,
    /// Allocated CPU-side buffers (mirrors GPU allocations).
    ///
    /// In the stub implementation these are plain `Vec<f64>` acting as
    /// stand-ins for actual `wgpu::Buffer` objects.  A full implementation
    /// wraps `wgpu::Buffer` behind `Arc<Mutex<…>>` to allow async reads.
    buffers: Vec<WgpuBufferEntry>,
    /// Whether the backend is operational.
    available: bool,
}

/// Internal buffer entry storing metadata and a CPU-side shadow copy.
#[derive(Debug, Clone)]
struct WgpuBufferEntry {
    /// Byte capacity of the GPU buffer (8 × `len` for f64 arrays).
    capacity: usize,
    /// CPU-side shadow for upload/download (avoids wgpu dep in stub).
    shadow: Vec<f64>,
}

impl WgpuBackend {
    /// Attempt to create a GPU-backed `WgpuBackend`.
    ///
    /// This stub type never owns a real GPU device, so it **always** returns
    /// `Err(WgpuInitError::NotAvailable)` — an honest "not available", never a
    /// fabricated success.  Construct a CPU-emulation instance explicitly with
    /// [`new_stub`](Self::new_stub) when you need to exercise the buffer API, or
    /// use the real on-device backend `real::WgpuBackendReal` (the `wgpu-backend`
    /// feature) for actual GPU compute.
    pub fn try_new() -> Result<Self, WgpuInitError> {
        // The real GPU adapter/device setup lives in `real::WgpuBackendReal`
        // (feature-gated): it calls `wgpu::Instance::request_adapter` and
        // `adapter.request_device`.  This stub owns no device, so it is honestly
        // unavailable rather than pretending to have initialised a GPU.
        Err(WgpuInitError::NotAvailable)
    }

    /// Create a stub backend for testing that stores data in CPU memory.
    ///
    /// This is equivalent to what `try_new` would return on a headless system
    /// but without returning an error — useful for unit testing backend logic.
    pub fn new_stub() -> Self {
        Self {
            device_info: WgpuDeviceInfo {
                name: "CPU stub".to_string(),
                backend: "None".to_string(),
                ..Default::default()
            },
            buffers: Vec::new(),
            available: false,
        }
    }

    /// Return `true` if a real GPU device is available.
    pub fn is_available(&self) -> bool {
        self.available
    }

    /// Return device information for diagnostics.
    pub fn device_info(&self) -> &WgpuDeviceInfo {
        &self.device_info
    }

    // ── Buffer management ────────────────────────────────────────────────────

    /// Allocate a GPU buffer that can hold `len` `f64` values.
    ///
    /// Returns a [`WgpuBufferHandle`] that can be passed to [`Self::write_buffer`]
    /// and [`Self::read_buffer`].
    ///
    /// In the stub implementation, a CPU-side shadow `Vec<f64>` is allocated.
    /// In the full wgpu implementation, `wgpu::Device::create_buffer` is called
    /// with `STORAGE | COPY_SRC | COPY_DST` usage flags.
    pub fn create_buffer(&mut self, len: usize) -> WgpuBufferHandle {
        let handle = WgpuBufferHandle(self.buffers.len());
        self.buffers.push(WgpuBufferEntry {
            capacity: len,
            shadow: vec![0.0; len],
        });
        handle
    }

    /// Upload `data` to the GPU buffer at `handle`.
    ///
    /// In the stub, data is copied into the CPU-side shadow.
    /// In the full implementation, `queue.write_buffer` is used.
    pub fn write_buffer(&mut self, handle: WgpuBufferHandle, data: &[f64]) {
        if let Some(entry) = self.buffers.get_mut(handle.0) {
            let len = data.len().min(entry.capacity);
            entry.shadow[..len].copy_from_slice(&data[..len]);
        }
    }

    /// Download data from the GPU buffer at `handle`.
    ///
    /// In the stub, data is read from the CPU-side shadow.
    /// In the full implementation, a staging buffer is created, the command
    /// `encoder.copy_buffer_to_buffer` is executed, and the staging buffer is
    /// mapped for reading.
    pub fn read_buffer(&self, handle: WgpuBufferHandle) -> Vec<f64> {
        self.buffers
            .get(handle.0)
            .map(|e| e.shadow.clone())
            .unwrap_or_default()
    }

    // ── Dispatch ─────────────────────────────────────────────────────────────

    /// Dispatch a compute kernel with `work_groups_x` workgroups.
    ///
    /// **This stub does not execute any kernel.**  It is a no-op pass-through:
    /// every buffer is left exactly as written.  It exists only so that code
    /// targeting the backend interface compiles and runs without a GPU — callers
    /// must therefore **not** treat the buffers as "computed" after a stub
    /// dispatch (doing so would turn unchanged inputs into a fabricated result).
    /// The real implementation `real::WgpuBackendReal::dispatch_wgsl` looks a
    /// `ComputePipeline` up from the shader cache and calls
    /// `encoder.dispatch_workgroups`.
    ///
    /// # Arguments
    ///
    /// * `kernel_name` — name of the WGSL shader entry point
    /// * `buffers`     — input/output buffer handles
    /// * `work_groups_x` — number of workgroups in the X dimension
    pub fn dispatch(
        &mut self,
        kernel_name: &str,
        buffers: &[WgpuBufferHandle],
        work_groups_x: u32,
    ) {
        // No-op: this stub owns no GPU device and cannot execute WGSL, so the
        // buffers are left unchanged.  Callers must not assume they were
        // modified — the real dispatch is `real::WgpuBackendReal::dispatch_wgsl`.
        let _ = (kernel_name, buffers, work_groups_x);
    }

    // ── WGSL shader registry ──────────────────────────────────────────────────

    /// Register a WGSL compute shader source under `name`.
    ///
    /// **This stub discards the source** — it neither stores nor compiles it,
    /// because the stub never dispatches WGSL.  The real implementation
    /// `real::WgpuBackendReal` compiles the WGSL lazily on first dispatch via
    /// `device.create_shader_module` and caches the resulting pipeline.
    pub fn register_shader(&mut self, name: &str, wgsl_source: &str) {
        // No-op: the stub does not execute WGSL, so there is nothing to compile
        // or store.  The real pipeline cache lives in `real::WgpuBackendReal`.
        let _ = (name, wgsl_source);
    }
}

// ── Built-in WGSL kernels ─────────────────────────────────────────────────────

/// WGSL source for a parallel prefix sum (exclusive scan) kernel.
///
/// This is the Blelloch algorithm adapted for WGSL with a workgroup of 256 threads.
pub const WGSL_PARALLEL_SCAN: &str = r#"
// Exclusive parallel prefix sum (Blelloch up-sweep / down-sweep)
// Workgroup size: 256 threads
// Binding 0: input buffer (read)
// Binding 1: output buffer (write)
// Binding 2: uniform { n: u32, pass: u32 }

@group(0) @binding(0) var<storage, read> input:  array<f32>;
@group(0) @binding(1) var<storage, read_write> output: array<f32>;

struct Params { n: u32, pass: u32 }
@group(0) @binding(2) var<uniform> params: Params;

var<workgroup> shared: array<f32, 256>;

@compute @workgroup_size(256)
fn exclusive_scan(@builtin(global_invocation_id) gid: vec3<u32>,
                  @builtin(local_invocation_id) lid: vec3<u32>) {
    let n = params.n;
    let i = gid.x;

    // Load
    shared[lid.x] = select(0.0, input[i], i < n);
    workgroupBarrier();

    // Up-sweep (reduce)
    var stride: u32 = 1u;
    loop {
        if stride >= 256u { break; }
        if lid.x % (stride * 2u) == (stride * 2u - 1u) {
            shared[lid.x] += shared[lid.x - stride];
        }
        workgroupBarrier();
        stride = stride * 2u;
    }

    // Down-sweep
    if lid.x == 255u { shared[255] = 0.0; }
    workgroupBarrier();
    stride = 128u;
    loop {
        if stride == 0u { break; }
        if lid.x % (stride * 2u) == (stride * 2u - 1u) {
            let tmp = shared[lid.x - stride];
            shared[lid.x - stride] = shared[lid.x];
            shared[lid.x] += tmp;
        }
        workgroupBarrier();
        stride = stride / 2u;
    }

    // Store
    if i < n { output[i] = shared[lid.x]; }
}
"#;

/// WGSL source for a simple SPH density kernel.
///
/// Computes particle density via a cubic-spline kernel with radius `h`.
pub const WGSL_SPH_DENSITY: &str = r#"
// SPH density kernel — W_spline3 smoothing
// Binding 0: positions array (x0,y0,z0, x1,y1,z1, ...)
// Binding 1: densities output (one per particle)
// Binding 2: uniform { n: u32, h: f32, mass: f32 }

struct SphParams { n: u32, h: f32, mass: f32 }
@group(0) @binding(0) var<storage, read>       positions: array<f32>;
@group(0) @binding(1) var<storage, read_write> densities: array<f32>;
@group(0) @binding(2) var<uniform>             params:    SphParams;

fn w_spline3(r: f32, h: f32) -> f32 {
    let q = r / h;
    let sigma = 3.0 / (2.0 * 3.14159265358979 * h * h * h);
    if q < 1.0 {
        return sigma * (2.0/3.0 - q*q + 0.5*q*q*q);
    } else if q < 2.0 {
        let t = 2.0 - q;
        return sigma * (1.0/6.0) * t*t*t;
    } else {
        return 0.0;
    }
}

@compute @workgroup_size(64)
fn sph_density(@builtin(global_invocation_id) gid: vec3<u32>) {
    let i = gid.x;
    let n = params.n;
    if i >= n { return; }

    let xi = vec3<f32>(positions[i*3u], positions[i*3u+1u], positions[i*3u+2u]);
    var density: f32 = 0.0;

    for (var j: u32 = 0u; j < n; j++) {
        let xj = vec3<f32>(positions[j*3u], positions[j*3u+1u], positions[j*3u+2u]);
        let r = length(xi - xj);
        density += params.mass * w_spline3(r, params.h);
    }

    densities[i] = density;
}
"#;

/// WGSL source for parallel BVH ray traversal.
///
/// Traverses a linearized BVH (LBVH) to find ray–box intersections.
/// This is a stub; real traversal requires the full BVH node buffer layout.
pub const WGSL_BVH_TRAVERSAL: &str = r#"
// Parallel BVH ray traversal stub
// Each thread handles one ray; BVH nodes are in binding 0.

struct Ray { origin: vec3<f32>, dir: vec3<f32>, t_max: f32 }
struct BvhNode { lo: vec3<f32>, hi: vec3<f32>, left: u32, right: u32, is_leaf: u32, prim: u32 }
struct HitResult { hit: u32, t: f32, prim: u32 }

@group(0) @binding(0) var<storage, read>       nodes:   array<BvhNode>;
@group(0) @binding(1) var<storage, read>        rays:    array<Ray>;
@group(0) @binding(2) var<storage, read_write> results: array<HitResult>;
@group(0) @binding(3) var<uniform>             num_rays: u32;

fn ray_aabb(ray: Ray, lo: vec3<f32>, hi: vec3<f32>) -> f32 {
    let inv_dir = 1.0 / ray.dir;
    let t0 = (lo - ray.origin) * inv_dir;
    let t1 = (hi - ray.origin) * inv_dir;
    let t_min = max(max(min(t0.x, t1.x), min(t0.y, t1.y)), min(t0.z, t1.z));
    let t_max_box = min(min(max(t0.x, t1.x), max(t0.y, t1.y)), max(t0.z, t1.z));
    if t_max_box < t_min || t_min > ray.t_max { return -1.0; }
    return t_min;
}

@compute @workgroup_size(64)
fn bvh_traverse(@builtin(global_invocation_id) gid: vec3<u32>) {
    let rid = gid.x;
    if rid >= num_rays { return; }
    let ray = rays[rid];
    results[rid] = HitResult(0u, ray.t_max, 0xFFFFFFFFu);

    // Iterative DFS stack (max depth 32)
    var stack: array<u32, 32>;
    var sp: i32 = 0;
    stack[0] = 0u;

    loop {
        if sp < 0 { break; }
        let node_idx = stack[sp]; sp--;
        let node = nodes[node_idx];

        let t = ray_aabb(ray, node.lo, node.hi);
        if t < 0.0 { continue; }

        if node.is_leaf != 0u {
            if t < results[rid].t {
                results[rid] = HitResult(1u, t, node.prim);
            }
        } else {
            if sp < 30 { sp++; stack[sp] = node.left; }
            if sp < 30 { sp++; stack[sp] = node.right; }
        }
    }
}
"#;

// ── WgpuInitError ─────────────────────────────────────────────────────────────

/// Error returned when the wgpu backend cannot be initialised.
#[derive(Debug, Clone, PartialEq)]
pub enum WgpuInitError {
    /// No compatible GPU adapter was found.
    NoAdapter,
    /// The `wgpu-backend` feature is not enabled; this is a stub build.
    NotAvailable,
    /// The device request failed (e.g. out of memory).
    DeviceRequestFailed(String),
    /// A required GPU feature is disabled or not supported.
    FeatureDisabled,
    /// Device creation failed with the given error string.
    DeviceRequest(String),
    /// A buffer handle is out of range.
    InvalidHandle(usize),
    /// A mutex was poisoned (should not occur in practice).
    PoisonedLock,
}

impl std::fmt::Display for WgpuInitError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            WgpuInitError::NoAdapter => write!(f, "No compatible GPU adapter found"),
            WgpuInitError::NotAvailable => write!(f, "wgpu-backend feature not enabled"),
            WgpuInitError::DeviceRequestFailed(s) => write!(f, "Device request failed: {s}"),
            WgpuInitError::FeatureDisabled => write!(f, "Required GPU feature is not available"),
            WgpuInitError::DeviceRequest(s) => write!(f, "Device request error: {s}"),
            WgpuInitError::InvalidHandle(h) => write!(f, "Invalid buffer handle: {h}"),
            WgpuInitError::PoisonedLock => write!(f, "Internal mutex was poisoned"),
        }
    }
}

impl std::error::Error for WgpuInitError {}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn try_new_returns_not_available_in_stub_build() {
        let result = WgpuBackend::try_new();
        assert!(matches!(result, Err(WgpuInitError::NotAvailable)));
    }

    #[test]
    fn stub_backend_write_read_roundtrip() {
        let mut backend = WgpuBackend::new_stub();
        let handle = backend.create_buffer(4);
        let data = vec![1.0_f64, 2.0, 3.0, 4.0];
        backend.write_buffer(handle, &data);
        let out = backend.read_buffer(handle);
        assert_eq!(out, data);
    }

    #[test]
    fn stub_dispatch_is_noop() {
        let mut backend = WgpuBackend::new_stub();
        let h = backend.create_buffer(8);
        let before = backend.read_buffer(h);
        backend.dispatch("sph_density", &[h], 1);
        let after = backend.read_buffer(h);
        assert_eq!(before, after, "stub dispatch should not modify buffers");
    }

    #[test]
    fn wgsl_kernels_are_non_empty() {
        assert!(!WGSL_PARALLEL_SCAN.is_empty());
        assert!(!WGSL_SPH_DENSITY.is_empty());
        assert!(!WGSL_BVH_TRAVERSAL.is_empty());
    }

    #[test]
    fn device_info_stub_has_name() {
        let backend = WgpuBackend::new_stub();
        assert!(!backend.device_info().name.is_empty());
    }

    #[test]
    fn wgpu_init_error_display() {
        assert!(!WgpuInitError::NotAvailable.to_string().is_empty());
        assert!(!WgpuInitError::NoAdapter.to_string().is_empty());
        assert!(!WgpuInitError::FeatureDisabled.to_string().is_empty());
        assert!(
            !WgpuInitError::DeviceRequest("oom".into())
                .to_string()
                .is_empty()
        );
        assert!(!WgpuInitError::InvalidHandle(7).to_string().is_empty());
        assert!(!WgpuInitError::PoisonedLock.to_string().is_empty());
    }
}

// ── Real wgpu backend (feature-gated) ─────────────────────────────────────────

/// Real wgpu compute backend, enabled only with the `wgpu-backend` feature.
///
/// Provides GPU buffer management, WGSL shader dispatch, and CPU-side readback
/// using `wgpu` 29's cross-platform Vulkan / Metal / DX12 backends.
///
/// # Thread safety
///
/// `wgpu::Device` and `wgpu::Queue` are `Send + Sync`.  The shader cache is
/// protected by a `Mutex`, making `WgpuBackendReal` safe to share across
/// threads (though individual dispatches are synchronous on the calling thread).
///
/// # Usage
///
/// ```ignore
/// // With the wgpu-backend feature enabled:
/// use oxiphysics_gpu::compute::wgpu_backend::real::WgpuBackendReal;
///
/// let mut backend = WgpuBackendReal::try_new()?;
/// let h = backend.create_buffer_f64(128);
/// backend.write_buffer_f64(h, &vec![1.0_f64; 128]);
/// backend.dispatch_wgsl(
///     WGSL_SPH_DENSITY, "sph_density",
///     &[(h, wgpu::BufferBindingType::Storage { read_only: false })],
///     [2, 1, 1],
/// )?;
/// let out = backend.read_buffer_f64(h);
/// ```
#[cfg(feature = "wgpu-backend")]
pub mod real {
    use super::{WgpuBufferHandle, WgpuDeviceInfo, WgpuInitError};
    use std::collections::HashMap;
    use std::hash::{DefaultHasher, Hash, Hasher};
    use std::sync::{Arc, Mutex};

    // ── Internal shader-cache entry ──────────────────────────────────────────

    struct ShaderCacheEntry {
        pipeline: Arc<wgpu::ComputePipeline>,
    }

    // ── WgpuBackendReal ──────────────────────────────────────────────────────

    /// Real GPU compute backend backed by `wgpu` 29.
    ///
    /// Obtain an instance via [`WgpuBackendReal::try_new`] (synchronous,
    /// blocks the thread) or [`WgpuBackendReal::try_new_async`] from within
    /// an async context.
    pub struct WgpuBackendReal {
        device: Arc<wgpu::Device>,
        queue: Arc<wgpu::Queue>,
        /// Device information (name, backend, driver).
        pub device_info: WgpuDeviceInfo,
        /// Allocated GPU buffers, indexed by `WgpuBufferHandle.0`.
        buffers: Vec<Option<Arc<wgpu::Buffer>>>,
        /// Byte size of each buffer (parallel to `buffers`).
        buffer_sizes: Vec<u64>,
        /// Compiled pipeline cache, keyed by a hash of WGSL source + entry point.
        shader_cache: Mutex<HashMap<u64, ShaderCacheEntry>>,
    }

    impl WgpuBackendReal {
        // ── Construction ─────────────────────────────────────────────────────

        /// Create a real GPU backend, blocking the calling thread.
        ///
        /// Returns `Err` if no compatible GPU adapter is found or if device
        /// creation fails.  Prefer [`try_new_async`](Self::try_new_async) from
        /// within an `async` context.
        pub fn try_new() -> Result<Self, WgpuInitError> {
            pollster::block_on(Self::try_new_async())
        }

        /// Create a real GPU backend asynchronously.
        ///
        /// This is the preferred entry point from `async` contexts (tokio,
        /// wasm-bindgen-futures, etc.).
        pub async fn try_new_async() -> Result<Self, WgpuInitError> {
            let instance =
                wgpu::Instance::new(wgpu::InstanceDescriptor::new_without_display_handle());

            let adapter = instance
                .request_adapter(&wgpu::RequestAdapterOptions {
                    power_preference: wgpu::PowerPreference::HighPerformance,
                    compatible_surface: None,
                    force_fallback_adapter: false,
                })
                .await
                .map_err(|_| WgpuInitError::NoAdapter)?;

            let info = adapter.get_info();

            let desc = wgpu::DeviceDescriptor {
                label: Some("oxiphysics-wgpu"),
                required_features: wgpu::Features::empty(),
                required_limits: adapter.limits(),
                ..Default::default()
            };

            let (device, queue) = adapter
                .request_device(&desc)
                .await
                .map_err(|e| WgpuInitError::DeviceRequest(e.to_string()))?;

            let device_info = WgpuDeviceInfo {
                name: info.name.clone(),
                backend: format!("{:?}", info.backend),
                driver_version: info.driver_info.clone(),
                // VRAM is not exposed by wgpu's AdapterInfo; use 0 as sentinel.
                vram_bytes: 0,
                // GPU-native f64 requires a device extension not in the base profile.
                supports_f64: false,
                // Conservative defaults matching most desktop GPU limits.
                max_workgroup_size: [256, 256, 64],
            };

            Ok(Self {
                device: Arc::new(device),
                queue: Arc::new(queue),
                device_info,
                buffers: Vec::new(),
                buffer_sizes: Vec::new(),
                shader_cache: Mutex::new(HashMap::new()),
            })
        }

        /// Return `true` — this struct always wraps a real GPU device.
        pub fn is_available(&self) -> bool {
            true
        }

        // ── Buffer management ─────────────────────────────────────────────────

        /// Allocate a GPU storage buffer of `size_bytes` bytes.
        ///
        /// The buffer is created with `STORAGE | COPY_SRC | COPY_DST` usage
        /// flags so that it can be used as a shader binding and for staged
        /// CPU read/write.
        pub fn create_buffer_storage(&mut self, size_bytes: u64) -> WgpuBufferHandle {
            let handle = WgpuBufferHandle(self.buffers.len());
            let buf = self.device.create_buffer(&wgpu::BufferDescriptor {
                label: None,
                size: size_bytes,
                usage: wgpu::BufferUsages::STORAGE
                    | wgpu::BufferUsages::COPY_SRC
                    | wgpu::BufferUsages::COPY_DST,
                mapped_at_creation: false,
            });
            self.buffers.push(Some(Arc::new(buf)));
            self.buffer_sizes.push(size_bytes);
            handle
        }

        /// Allocate a GPU buffer sized for `len` `f64` values.
        ///
        /// Internally the data is stored as `f32` on the GPU (8 bytes per
        /// element to maintain the same stride).
        pub fn create_buffer_f64(&mut self, len: usize) -> WgpuBufferHandle {
            // We store f64 values packed as two f32s to preserve stride; or
            // simply allocate 8 bytes per element and use the f32 path with
            // two floats per logical element. For simplicity, the current
            // implementation casts f64→f32 on write and f32→f64 on read, so
            // we only need 4 bytes per element on the GPU.
            self.create_buffer_storage((len * 4) as u64)
        }

        /// Upload `data` to the GPU buffer at `handle`, casting `f64` → `f32`.
        ///
        /// # Panics
        ///
        /// Does nothing (silently returns) if `handle` is out of range.
        pub fn write_buffer_f64(&self, handle: WgpuBufferHandle, data: &[f64]) {
            if let Some(Some(buf)) = self.buffers.get(handle.0) {
                let f32_data: Vec<f32> = data.iter().map(|&v| v as f32).collect();
                self.queue
                    .write_buffer(buf, 0, bytemuck::cast_slice(&f32_data));
            }
        }

        /// Download data from the GPU buffer at `handle`, casting `f32` → `f64`.
        ///
        /// This blocks the calling thread until the GPU has finished all
        /// outstanding work and the readback mapping is complete.
        ///
        /// Returns an empty `Vec` if the handle is invalid or the readback fails.
        pub fn read_buffer_f64(&self, handle: WgpuBufferHandle) -> Vec<f64> {
            let buf = match self.buffers.get(handle.0).and_then(|b| b.as_ref()) {
                Some(b) => b.clone(),
                None => return Vec::new(),
            };
            let size = self.buffer_sizes[handle.0];

            // Create a CPU-visible staging buffer for the readback.
            let staging = self.device.create_buffer(&wgpu::BufferDescriptor {
                label: Some("oxiphysics_staging_readback"),
                size,
                usage: wgpu::BufferUsages::MAP_READ | wgpu::BufferUsages::COPY_DST,
                mapped_at_creation: false,
            });

            // Record and submit the copy command.
            let mut encoder = self
                .device
                .create_command_encoder(&wgpu::CommandEncoderDescriptor { label: None });
            encoder.copy_buffer_to_buffer(&buf, 0, &staging, 0, size);
            self.queue.submit(std::iter::once(encoder.finish()));

            // Map the staging buffer for reading.
            let slice = staging.slice(..);
            let (tx, rx) = std::sync::mpsc::channel();
            slice.map_async(wgpu::MapMode::Read, move |result| {
                let _ = tx.send(result);
            });

            // Block until the GPU has completed and the mapping is ready.
            if let Err(_e) = self.device.poll(wgpu::PollType::Wait {
                submission_index: None,
                timeout: None,
            }) {
                return Vec::new();
            }

            // Check that the mapping succeeded.
            if rx.recv().ok().and_then(|r| r.ok()).is_none() {
                return Vec::new();
            }

            let mapped = slice.get_mapped_range();
            let f32_data: &[f32] = bytemuck::cast_slice(&mapped);
            let result: Vec<f64> = f32_data.iter().map(|&v| v as f64).collect();
            drop(mapped);
            staging.unmap();
            result
        }

        // ── Dispatch ──────────────────────────────────────────────────────────

        /// Upload raw bytes to the GPU buffer at `handle`.
        ///
        /// The byte slice must fit within the buffer's allocated size.
        /// Does nothing (silently returns) if `handle` is out of range.
        pub fn queue_write_buffer_raw(&self, handle: &WgpuBufferHandle, data: &[u8]) {
            if let Some(Some(buf)) = self.buffers.get(handle.0) {
                self.queue.write_buffer(buf, 0, data);
            }
        }

        /// Upload `f32` data directly to the GPU buffer at `handle` (no f64→f32 cast).
        ///
        /// Does nothing (silently returns) if `handle` is out of range.
        pub fn queue_write_buffer_f32(&self, handle: &WgpuBufferHandle, data: &[f32]) {
            if let Some(Some(buf)) = self.buffers.get(handle.0) {
                self.queue.write_buffer(buf, 0, bytemuck::cast_slice(data));
            }
        }

        /// Download raw `f32` values from the GPU buffer at `handle`.
        ///
        /// Returns an empty `Vec` if the handle is invalid or the readback fails.
        pub fn read_buffer_f32(&self, handle: WgpuBufferHandle) -> Vec<f32> {
            let buf = match self.buffers.get(handle.0).and_then(|b| b.as_ref()) {
                Some(b) => b.clone(),
                None => return Vec::new(),
            };
            let size = self.buffer_sizes[handle.0];

            let staging = self.device.create_buffer(&wgpu::BufferDescriptor {
                label: Some("oxiphysics_staging_readback_f32"),
                size,
                usage: wgpu::BufferUsages::MAP_READ | wgpu::BufferUsages::COPY_DST,
                mapped_at_creation: false,
            });

            let mut encoder = self
                .device
                .create_command_encoder(&wgpu::CommandEncoderDescriptor { label: None });
            encoder.copy_buffer_to_buffer(&buf, 0, &staging, 0, size);
            self.queue.submit(std::iter::once(encoder.finish()));

            let slice = staging.slice(..);
            let (tx, rx) = std::sync::mpsc::channel();
            slice.map_async(wgpu::MapMode::Read, move |result| {
                let _ = tx.send(result);
            });

            if let Err(_e) = self.device.poll(wgpu::PollType::Wait {
                submission_index: None,
                timeout: None,
            }) {
                return Vec::new();
            }

            if rx.recv().ok().and_then(|r| r.ok()).is_none() {
                return Vec::new();
            }

            let mapped = slice.get_mapped_range();
            let result: Vec<f32> = bytemuck::cast_slice::<u8, f32>(&mapped).to_vec();
            drop(mapped);
            staging.unmap();
            result
        }

        /// Download raw `u32` values from the GPU buffer at `handle`.
        ///
        /// Mirrors [`read_buffer_f32`](Self::read_buffer_f32) but reinterprets
        /// the mapped bytes as `u32`.  This is required for integer readbacks
        /// (scan, compaction, histogram) where a float round-trip would corrupt
        /// NaN-pattern bit values.
        ///
        /// Returns an empty `Vec` if the handle is invalid or the readback fails.
        pub fn read_buffer_u32(&self, handle: WgpuBufferHandle) -> Vec<u32> {
            let buf = match self.buffers.get(handle.0).and_then(|b| b.as_ref()) {
                Some(b) => b.clone(),
                None => return Vec::new(),
            };
            let size = self.buffer_sizes[handle.0];

            let staging = self.device.create_buffer(&wgpu::BufferDescriptor {
                label: Some("oxiphysics_staging_readback_u32"),
                size,
                usage: wgpu::BufferUsages::MAP_READ | wgpu::BufferUsages::COPY_DST,
                mapped_at_creation: false,
            });

            let mut encoder = self
                .device
                .create_command_encoder(&wgpu::CommandEncoderDescriptor { label: None });
            encoder.copy_buffer_to_buffer(&buf, 0, &staging, 0, size);
            self.queue.submit(std::iter::once(encoder.finish()));

            let slice = staging.slice(..);
            let (tx, rx) = std::sync::mpsc::channel();
            slice.map_async(wgpu::MapMode::Read, move |result| {
                let _ = tx.send(result);
            });

            if let Err(_e) = self.device.poll(wgpu::PollType::Wait {
                submission_index: None,
                timeout: None,
            }) {
                return Vec::new();
            }

            if rx.recv().ok().and_then(|r| r.ok()).is_none() {
                return Vec::new();
            }

            let mapped = slice.get_mapped_range();
            let result: Vec<u32> = bytemuck::cast_slice::<u8, u32>(&mapped).to_vec();
            drop(mapped);
            staging.unmap();
            result
        }

        // ── Dispatch ──────────────────────────────────────────────────────────

        /// Compute the 3-D workgroup dispatch counts for `n_items` elements.
        ///
        /// Returns `[0, 1, 1]` for `n_items == 0` (no-op dispatch).
        pub fn dispatch_count_for(n_items: usize, workgroup_size: u32) -> [u32; 3] {
            crate::compute::timestamp::dispatch_count_for(n_items, workgroup_size)
        }

        /// Compile and dispatch a WGSL compute shader.
        ///
        /// The pipeline is compiled lazily and cached by a hash of
        /// `(wgsl_src, entry_point)`, so repeated calls with the same shader
        /// do not recompile.
        ///
        /// # Parameters
        ///
        /// * `wgsl_src`    — WGSL shader source code.
        /// * `entry_point` — Name of the `@compute` entry point function.
        /// * `buffers`     — Ordered list of `(handle, binding_type)` pairs.
        ///   Binding index in the WGSL shader corresponds to the position in
        ///   this slice (binding 0 = `buffers[0]`, etc.).
        /// * `workgroups`  — `[x, y, z]` dispatch counts.
        ///
        /// # Errors
        ///
        /// Returns `Err(WgpuInitError::InvalidHandle)` if any buffer handle is
        /// out of range.  Returns `Err(WgpuInitError::PoisonedLock)` if the
        /// shader-cache mutex is poisoned (should not occur in practice).
        pub fn dispatch_wgsl(
            &self,
            wgsl_src: &str,
            entry_point: &str,
            buffers: &[(WgpuBufferHandle, wgpu::BufferBindingType)],
            workgroups: [u32; 3],
        ) -> Result<(), WgpuInitError> {
            // Hash the shader source + entry point to key the pipeline cache.
            let mut hasher = DefaultHasher::new();
            wgsl_src.hash(&mut hasher);
            entry_point.hash(&mut hasher);
            let key = hasher.finish();

            // Obtain or compile the pipeline.
            let pipeline: Arc<wgpu::ComputePipeline> = {
                let mut cache = self.shader_cache.lock().unwrap_or_else(|e| e.into_inner());

                if let Some(entry) = cache.get(&key) {
                    entry.pipeline.clone()
                } else {
                    let module = self
                        .device
                        .create_shader_module(wgpu::ShaderModuleDescriptor {
                            label: Some(entry_point),
                            source: wgpu::ShaderSource::Wgsl(wgsl_src.into()),
                        });
                    let pipeline = Arc::new(self.device.create_compute_pipeline(
                        &wgpu::ComputePipelineDescriptor {
                            label: Some(entry_point),
                            layout: None,
                            module: &module,
                            entry_point: Some(entry_point),
                            compilation_options: wgpu::PipelineCompilationOptions::default(),
                            cache: None,
                        },
                    ));
                    cache.insert(
                        key,
                        ShaderCacheEntry {
                            pipeline: pipeline.clone(),
                        },
                    );
                    pipeline
                }
            };

            // Derive the bind-group layout from the compiled pipeline.
            let bg_layout = pipeline.get_bind_group_layout(0);

            // Build the bind-group entries.
            let mut entries: Vec<wgpu::BindGroupEntry> = Vec::with_capacity(buffers.len());
            for (i, (handle, _binding_type)) in buffers.iter().enumerate() {
                let buf = self
                    .buffers
                    .get(handle.0)
                    .and_then(|b| b.as_ref())
                    .ok_or(WgpuInitError::InvalidHandle(handle.0))?;
                entries.push(wgpu::BindGroupEntry {
                    binding: i as u32,
                    resource: buf.as_entire_binding(),
                });
            }

            let bind_group = self.device.create_bind_group(&wgpu::BindGroupDescriptor {
                label: None,
                layout: &bg_layout,
                entries: &entries,
            });

            // Record and submit the compute pass.
            let mut encoder = self
                .device
                .create_command_encoder(&wgpu::CommandEncoderDescriptor { label: None });
            {
                let mut pass = encoder.begin_compute_pass(&wgpu::ComputePassDescriptor {
                    label: None,
                    timestamp_writes: None,
                });
                pass.set_pipeline(&pipeline);
                pass.set_bind_group(0, &bind_group, &[]);
                pass.dispatch_workgroups(workgroups[0], workgroups[1], workgroups[2]);
            }
            self.queue.submit(std::iter::once(encoder.finish()));

            // Block until the GPU has finished (synchronous dispatch).
            self.device
                .poll(wgpu::PollType::Wait {
                    submission_index: None,
                    timeout: None,
                })
                .map_err(|_| WgpuInitError::DeviceRequest("poll failed".into()))?;

            Ok(())
        }
    }

    // ── Feature-gated tests ───────────────────────────────────────────────────

    #[cfg(test)]
    mod tests {
        use super::*;

        /// Helper: attempt to create a real backend, returning `None` if no GPU
        /// is available (e.g. in headless CI).
        fn try_backend() -> Option<WgpuBackendReal> {
            WgpuBackendReal::try_new().ok()
        }

        #[test]
        fn real_backend_try_new_succeeds_or_gracefully_fails() {
            // This test always passes: it either succeeds (GPU present) or
            // returns None (headless / CI environment).
            match WgpuBackendReal::try_new() {
                Ok(b) => {
                    assert!(b.is_available());
                    assert!(!b.device_info.backend.is_empty());
                }
                Err(e) => {
                    // NoAdapter is the expected error in headless CI.
                    eprintln!("No GPU adapter available: {e}");
                }
            }
        }

        #[test]
        fn real_backend_create_and_write_buffer() {
            let Some(mut backend) = try_backend() else {
                return;
            };
            let data = vec![1.0_f64, 2.0, 3.0, 4.0];
            let handle = backend.create_buffer_f64(data.len());
            backend.write_buffer_f64(handle, &data);
            // write_buffer_f64 is fire-and-forget; we just verify no panic.
            assert!(handle.0 < backend.buffers.len());
        }

        #[test]
        fn real_backend_buffer_roundtrip() {
            let Some(mut backend) = try_backend() else {
                return;
            };
            let data = vec![1.0_f64, 2.0, 3.0, 4.0];
            let handle = backend.create_buffer_f64(data.len());
            backend.write_buffer_f64(handle, &data);
            let out = backend.read_buffer_f64(handle);
            // f64→f32→f64 loses precision; check within f32 rounding.
            assert_eq!(out.len(), data.len());
            for (&expected, &got) in data.iter().zip(out.iter()) {
                assert!(
                    (expected as f32 - got as f32).abs() < 1e-5,
                    "roundtrip mismatch: expected {expected}, got {got}"
                );
            }
        }

        #[test]
        fn real_backend_dispatch_scale_shader() {
            let Some(mut backend) = try_backend() else {
                return;
            };
            use super::super::WgpuBackend;

            // A simple WGSL shader that multiplies each f32 element by 2.
            const SCALE_BY_TWO: &str = r#"
@group(0) @binding(0) var<storage, read>       input_buf:  array<f32>;
@group(0) @binding(1) var<storage, read_write> output_buf: array<f32>;

@compute @workgroup_size(64)
fn scale_by_two(@builtin(global_invocation_id) gid: vec3<u32>) {
    let i = gid.x;
    if i < arrayLength(&input_buf) {
        output_buf[i] = input_buf[i] * 2.0;
    }
}
"#;
            let n: usize = 4;
            let input_data: Vec<f32> = (1..=n as u32).map(|x| x as f32).collect();
            let in_handle = backend.create_buffer_storage((n * 4) as u64);
            let out_handle = backend.create_buffer_storage((n * 4) as u64);

            backend.queue.write_buffer(
                backend.buffers[in_handle.0].as_ref().unwrap(),
                0,
                bytemuck::cast_slice(&input_data),
            );

            // Dispatch: 1 workgroup of 64 threads covers n=4 elements.
            let result = backend.dispatch_wgsl(
                SCALE_BY_TWO,
                "scale_by_two",
                &[
                    (
                        in_handle,
                        wgpu::BufferBindingType::Storage { read_only: true },
                    ),
                    (
                        out_handle,
                        wgpu::BufferBindingType::Storage { read_only: false },
                    ),
                ],
                [1, 1, 1],
            );
            assert!(result.is_ok(), "dispatch_wgsl failed: {:?}", result.err());

            // Readback via staging and verify.
            let out = backend.read_buffer_f64(out_handle);
            assert_eq!(out.len(), n);
            for (i, &v) in out.iter().enumerate() {
                let expected = (i + 1) as f64 * 2.0;
                assert!(
                    (v - expected).abs() < 0.01,
                    "element {i}: expected {expected}, got {v}"
                );
            }

            // Regression guard: stub backend still works.
            let mut stub = WgpuBackend::new_stub();
            let h = stub.create_buffer(4);
            let _ = stub.read_buffer(h);
        }

        #[test]
        fn dispatch_count_for_zero_items() {
            assert_eq!(WgpuBackendReal::dispatch_count_for(0, 64), [0, 1, 1]);
        }

        #[test]
        fn dispatch_count_for_65_items() {
            assert_eq!(WgpuBackendReal::dispatch_count_for(65, 64), [2, 1, 1]);
        }

        #[test]
        fn dispatch_count_for_exact_workgroup() {
            assert_eq!(WgpuBackendReal::dispatch_count_for(256, 64), [4, 1, 1]);
        }
    }
}
