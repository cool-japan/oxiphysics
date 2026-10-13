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
                  @builtin(local_invocation_id) lid: vec3<u32>,
                  @builtin(num_workgroups) nwg: vec3<u32>) {
    let n = params.n;
    // Dispatches wider than one dimension allows are folded into y.
    let i = gid.x + gid.y * nwg.x * 256u;

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
fn sph_density(@builtin(global_invocation_id) gid: vec3<u32>,
               @builtin(num_workgroups) nwg: vec3<u32>) {
    // Dispatches wider than one dimension allows are folded into y.
    let i = gid.x + gid.y * nwg.x * 64u;
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
fn bvh_traverse(@builtin(global_invocation_id) gid: vec3<u32>,
                @builtin(num_workgroups) nwg: vec3<u32>) {
    // Dispatches wider than one dimension allows are folded into y.
    let rid = gid.x + gid.y * nwg.x * 64u;
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
    /// Mapping a GPU buffer for readback failed: the map request was rejected,
    /// its callback never ran, or the mapped range could not be obtained.
    MapFailed(String),
    /// Waiting for the GPU to finish its work (a device poll) failed.
    PollFailed(String),
    /// The GPU rejected an operation.
    ///
    /// Carries the description of the validation, out-of-memory or internal
    /// error that wgpu reported - for example invalid WGSL, a binding that does
    /// not match the shader, a dispatch beyond the device limits, or an invalid
    /// queue write.  Such errors are returned instead of aborting the caller.
    Validation(String),
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
            WgpuInitError::MapFailed(s) => write!(f, "Buffer readback mapping failed: {s}"),
            WgpuInitError::PollFailed(s) => write!(f, "Waiting for the GPU failed: {s}"),
            WgpuInitError::Validation(s) => write!(f, "GPU error: {s}"),
        }
    }
}

impl std::error::Error for WgpuInitError {}

impl From<crate::Error> for WgpuInitError {
    /// A crate-level error (for example a dispatch beyond the workgroup limits)
    /// reaches the backend as a [`WgpuInitError::Validation`] with its message.
    fn from(error: crate::Error) -> Self {
        WgpuInitError::Validation(error.to_string())
    }
}

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

    #[test]
    fn wgpu_init_error_new_variants_display_their_message() {
        let map = WgpuInitError::MapFailed("range is not mapped".into()).to_string();
        assert!(map.contains("range is not mapped"), "{map}");
        let poll = WgpuInitError::PollFailed("timed out".into()).to_string();
        assert!(poll.contains("timed out"), "{poll}");
        let validation = WgpuInitError::Validation("bad binding".into()).to_string();
        assert!(validation.contains("bad binding"), "{validation}");
    }

    #[test]
    fn crate_error_converts_to_validation_error() {
        let converted: WgpuInitError = crate::Error::DispatchLimitExceeded {
            dispatch_size: 70_000,
            limit: 65_535,
        }
        .into();
        match converted {
            WgpuInitError::Validation(message) => {
                assert!(
                    message.contains("70000") && message.contains("65535"),
                    "{message}"
                );
            }
            other => panic!("expected Validation, got {other:?}"),
        }
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
/// let out = backend.read_buffer_f64(h)?;
/// ```
///
/// # Errors instead of aborts
///
/// wgpu reports invalid shaders, binding mismatches, oversized dispatches and
/// similar mistakes as validation errors, which abort the process by default.
/// This backend installs its own handler for such errors and wraps shader
/// compilation, dispatch and readback in error scopes, so they come back as
/// [`WgpuInitError::Validation`].  An error raised by an operation that has no
/// result of its own (for example `queue_write_buffer_*` with a bad size, or a
/// failed allocation) is kept and returned by the next
/// [`dispatch_wgsl`](real::WgpuBackendReal::dispatch_wgsl) or `read_buffer_*`
/// call (see [`take_pending_error`](real::WgpuBackendReal::take_pending_error)).
#[cfg(feature = "wgpu-backend")]
pub mod real {
    use super::{WgpuBufferHandle, WgpuDeviceInfo, WgpuInitError};
    use std::collections::HashMap;
    use std::hash::{DefaultHasher, Hash, Hasher};
    use std::sync::{Arc, Mutex, PoisonError};
    use std::time::Duration;

    /// How long a readback waits for the map callback after the device poll
    /// returned (the callback normally ran before the poll returns).
    const MAP_CALLBACK_TIMEOUT: Duration = Duration::from_secs(30);

    /// One-line description of a wgpu error, prefixed with its class.
    fn describe_error(error: &wgpu::Error) -> String {
        match error {
            wgpu::Error::OutOfMemory { .. } => "out of memory".to_string(),
            wgpu::Error::Validation { description, .. } => {
                format!("validation error: {description}")
            }
            wgpu::Error::Internal { description, .. } => {
                format!("internal error: {description}")
            }
        }
    }

    /// Keep `message` in `slot` unless an earlier error is already waiting: the
    /// first error is the root cause, later ones are usually its consequences.
    fn record_error(slot: &Mutex<Option<String>>, message: String) {
        let mut guard = slot.lock().unwrap_or_else(PoisonError::into_inner);
        guard.get_or_insert(message);
    }

    /// Out-of-memory, internal and validation error scopes pushed around a
    /// piece of GPU work.
    ///
    /// The fields are declared innermost scope first: scopes must be popped in
    /// reverse order of their creation, and fields drop in declaration order,
    /// so a value dropped on an early return still unwinds the stack correctly.
    struct ErrorScopes {
        validation: wgpu::ErrorScopeGuard,
        internal: wgpu::ErrorScopeGuard,
        out_of_memory: wgpu::ErrorScopeGuard,
    }

    impl ErrorScopes {
        fn push(device: &wgpu::Device) -> Self {
            let out_of_memory = device.push_error_scope(wgpu::ErrorFilter::OutOfMemory);
            let internal = device.push_error_scope(wgpu::ErrorFilter::Internal);
            let validation = device.push_error_scope(wgpu::ErrorFilter::Validation);
            Self {
                validation,
                internal,
                out_of_memory,
            }
        }

        /// Pop all scopes, innermost first, and report the first error caught.
        fn finish(self) -> Result<(), WgpuInitError> {
            let validation = pollster::block_on(self.validation.pop());
            let internal = pollster::block_on(self.internal.pop());
            let out_of_memory = pollster::block_on(self.out_of_memory.pop());
            match validation.or(internal).or(out_of_memory) {
                Some(error) => Err(WgpuInitError::Validation(describe_error(&error))),
                None => Ok(()),
            }
        }
    }

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
        /// Byte size of each buffer (parallel to `buffers`; `0` for a freed slot).
        buffer_sizes: Vec<u64>,
        /// Indices of destroyed buffer slots, reused by `create_buffer_storage`
        /// so that the slot table does not grow with every create/destroy cycle.
        free_slots: Vec<usize>,
        /// Compiled pipeline cache, keyed by a hash of WGSL source + entry point.
        shader_cache: Mutex<HashMap<u64, ShaderCacheEntry>>,
        /// First error wgpu reported outside any error scope, kept until the next
        /// fallible call returns it.  Shared with the device's uncaptured-error
        /// handler.
        pending_error: Arc<Mutex<Option<String>>>,
        /// `max_compute_workgroups_per_dimension` of the device.
        max_workgroups_per_dimension: u32,
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
                    apply_limit_buckets: false,
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

            // wgpu aborts the process on errors nobody handles; route them into
            // `pending_error` so that the next fallible call can return them.
            let pending_error: Arc<Mutex<Option<String>>> = Arc::new(Mutex::new(None));
            let error_sink = Arc::clone(&pending_error);
            device.on_uncaptured_error(Arc::new(move |error: wgpu::Error| {
                record_error(&error_sink, describe_error(&error));
            }));
            let max_workgroups_per_dimension = device.limits().max_compute_workgroups_per_dimension;

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
                free_slots: Vec::new(),
                shader_cache: Mutex::new(HashMap::new()),
                pending_error,
                max_workgroups_per_dimension,
            })
        }

        /// `max_compute_workgroups_per_dimension` of the device: the largest
        /// workgroup count `dispatch_wgsl` accepts along one dimension.
        pub fn max_workgroups_per_dimension(&self) -> u32 {
            self.max_workgroups_per_dimension
        }

        /// Take the first error wgpu reported outside an error scope since the
        /// last call, if any.
        ///
        /// Operations without a result of their own - `queue_write_buffer_*`,
        /// buffer allocation - cannot return a failure.  Their errors (a write
        /// larger than the buffer, an allocation beyond the device limits, ...)
        /// are kept instead of aborting the process and are returned by the next
        /// `dispatch_wgsl` / `read_buffer_*` call; call this to look at them
        /// sooner.
        pub fn take_pending_error(&self) -> Option<String> {
            self.pending_error
                .lock()
                .unwrap_or_else(PoisonError::into_inner)
                .take()
        }

        /// Fail with the pending uncaptured error, if there is one.
        fn check_pending_error(&self) -> Result<(), WgpuInitError> {
            match self.take_pending_error() {
                Some(message) => Err(WgpuInitError::Validation(message)),
                None => Ok(()),
            }
        }

        /// Run `body` inside error scopes.
        ///
        /// An error wgpu raises while `body` runs is returned as
        /// [`WgpuInitError::Validation`] instead of reaching the uncaptured-error
        /// handler.  An error returned by `body` itself takes precedence.
        fn scoped<T>(
            &self,
            body: impl FnOnce() -> Result<T, WgpuInitError>,
        ) -> Result<T, WgpuInitError> {
            let scopes = ErrorScopes::push(&self.device);
            let outcome = body();
            let scope_result = scopes.finish();
            let value = outcome?;
            scope_result?;
            Ok(value)
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
        ///
        /// The size is rounded up to a multiple of 4 bytes with a minimum of
        /// 4: wgpu cannot copy, map or bind a buffer whose size is not 4-byte
        /// aligned or that is empty, and would report a validation error for
        /// every later use.  [`buffer_size`](Self::buffer_size) returns the
        /// rounded size.  A request beyond the device's buffer-size limit is
        /// reported by the next `dispatch_wgsl` / `read_buffer_*` call (see
        /// [`take_pending_error`](Self::take_pending_error)).
        ///
        /// A slot freed by [`destroy_buffer`](Self::destroy_buffer) is reused
        /// before the slot table grows, so a handle value can reappear after its
        /// buffer was destroyed.  Every copy of a destroyed handle is therefore
        /// invalid and must be dropped by the caller.
        pub fn create_buffer_storage(&mut self, size_bytes: u64) -> WgpuBufferHandle {
            let size_bytes = size_bytes.max(4).saturating_add(3) & !3;
            let buf = Arc::new(self.device.create_buffer(&wgpu::BufferDescriptor {
                label: None,
                size: size_bytes,
                usage: wgpu::BufferUsages::STORAGE
                    | wgpu::BufferUsages::COPY_SRC
                    | wgpu::BufferUsages::COPY_DST,
                mapped_at_creation: false,
            }));
            if let Some(slot) = self.free_slots.pop()
                && let (Some(entry), Some(size)) =
                    (self.buffers.get_mut(slot), self.buffer_sizes.get_mut(slot))
            {
                *entry = Some(buf);
                *size = size_bytes;
                return WgpuBufferHandle(slot);
            }
            let handle = WgpuBufferHandle(self.buffers.len());
            self.buffers.push(Some(buf));
            self.buffer_sizes.push(size_bytes);
            handle
        }

        /// Free the GPU buffer behind `handle` immediately.
        ///
        /// The native resource is released with `wgpu::Buffer::destroy` (not left
        /// to garbage collection), the slot is cleared and its index goes onto a
        /// free-list that [`create_buffer_storage`](Self::create_buffer_storage)
        /// reuses.  Returns `true` if a live buffer was destroyed and `false` for
        /// an unknown or already destroyed handle (idempotent).
        ///
        /// After this call every copy of `handle` is invalid: using one in a
        /// dispatch yields `WgpuInitError::InvalidHandle` until the slot is
        /// handed out again by a later allocation.
        pub fn destroy_buffer(&mut self, handle: WgpuBufferHandle) -> bool {
            let Some(slot) = self.buffers.get_mut(handle.0) else {
                return false;
            };
            let Some(buf) = slot.take() else {
                return false;
            };
            buf.destroy();
            if let Some(size) = self.buffer_sizes.get_mut(handle.0) {
                *size = 0;
            }
            self.free_slots.push(handle.0);
            true
        }

        /// Byte size of the live buffer behind `handle`, or `None` if the handle
        /// is unknown or the buffer was destroyed.
        pub fn buffer_size(&self, handle: WgpuBufferHandle) -> Option<u64> {
            match self.buffers.get(handle.0) {
                Some(Some(_)) => self.buffer_sizes.get(handle.0).copied(),
                _ => None,
            }
        }

        /// Make sure the buffer behind `*handle` holds at least `min_bytes` bytes.
        ///
        /// If the buffer is already large enough nothing happens and `false` is
        /// returned.  Otherwise a new buffer is allocated with geometric growth
        /// (`max(min_bytes, 2 × old size)`, rounded up to a multiple of 4 and at
        /// least 16 bytes), the old buffer is destroyed, `*handle` is updated to
        /// the new buffer and `true` is returned.  A destroyed or unknown handle
        /// is simply replaced by a fresh allocation.
        ///
        /// The previous contents are **not** copied, so this is meant for
        /// per-call scratch buffers (rays, results, staging) that are rewritten
        /// before every dispatch.
        pub fn ensure_buffer(&mut self, handle: &mut WgpuBufferHandle, min_bytes: u64) -> bool {
            let current = self.buffer_size(*handle);
            if let Some(size) = current
                && size >= min_bytes
            {
                return false;
            }
            let grown = match current {
                Some(size) => min_bytes.max(size.saturating_mul(2)),
                None => min_bytes,
            };
            let new_size = grown.max(16).div_ceil(4) * 4;
            let replacement = self.create_buffer_storage(new_size);
            if current.is_some() {
                self.destroy_buffer(*handle);
            }
            *handle = replacement;
            true
        }

        /// Number of buffers currently alive on this backend (destroyed slots
        /// excluded).  Constant in steady state when callers recycle their
        /// scratch buffers instead of allocating per call.
        pub fn live_buffer_count(&self) -> usize {
            self.buffers.iter().filter(|b| b.is_some()).count()
        }

        /// Total bytes of all live buffers (destroyed slots count as zero).
        pub fn live_buffer_bytes(&self) -> u64 {
            self.buffers
                .iter()
                .zip(self.buffer_sizes.iter())
                .filter(|(b, _)| b.is_some())
                .map(|(_, &size)| size)
                .sum()
        }

        /// Number of compute pipelines compiled and cached so far (one per distinct
        /// WGSL source and entry point).
        ///
        /// A long-lived workload compiles its pipelines on the first dispatch and
        /// this count then stays constant; a count that keeps growing means the
        /// workload recreates the backend or feeds new shader text on every call.
        pub fn compiled_pipeline_count(&self) -> usize {
            self.shader_cache
                .lock()
                .unwrap_or_else(|e| e.into_inner())
                .len()
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
        /// # Errors
        ///
        /// See [`read_buffer_f32`](Self::read_buffer_f32).
        pub fn read_buffer_f64(&self, handle: WgpuBufferHandle) -> Result<Vec<f64>, WgpuInitError> {
            Ok(self
                .read_buffer_pod::<f32>(handle)?
                .into_iter()
                .map(f64::from)
                .collect())
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
        /// The result holds the whole buffer (see
        /// [`buffer_size`](Self::buffer_size)); callers that use a prefix of it
        /// must check the length they need.
        ///
        /// # Errors
        ///
        /// * [`WgpuInitError::InvalidHandle`] if `handle` is unknown or destroyed.
        /// * [`WgpuInitError::MapFailed`] if the staging buffer could not be mapped.
        /// * [`WgpuInitError::PollFailed`] if waiting for the GPU failed.
        /// * [`WgpuInitError::Validation`] if wgpu rejected the copy, or an
        ///   earlier operation left an error behind (see
        ///   [`take_pending_error`](Self::take_pending_error)).
        ///
        /// A failed readback is never reported as an empty or zero-filled `Vec`.
        pub fn read_buffer_f32(&self, handle: WgpuBufferHandle) -> Result<Vec<f32>, WgpuInitError> {
            self.read_buffer_pod::<f32>(handle)
        }

        /// Download raw `u32` values from the GPU buffer at `handle`.
        ///
        /// Mirrors [`read_buffer_f32`](Self::read_buffer_f32) but reinterprets
        /// the mapped bytes as `u32`.  This is required for integer readbacks
        /// (scan, compaction, histogram) where a float round-trip would corrupt
        /// NaN-pattern bit values.
        ///
        /// # Errors
        ///
        /// See [`read_buffer_f32`](Self::read_buffer_f32).
        pub fn read_buffer_u32(&self, handle: WgpuBufferHandle) -> Result<Vec<u32>, WgpuInitError> {
            self.read_buffer_pod::<u32>(handle)
        }

        /// Copy the whole buffer behind `handle` into a staging buffer, map it
        /// and return its contents as a `Vec<T>`.
        fn read_buffer_pod<T: bytemuck::Pod>(
            &self,
            handle: WgpuBufferHandle,
        ) -> Result<Vec<T>, WgpuInitError> {
            self.check_pending_error()?;
            let buffer = self
                .buffers
                .get(handle.0)
                .and_then(|b| b.as_ref())
                .ok_or(WgpuInitError::InvalidHandle(handle.0))?;
            let size = self
                .buffer_sizes
                .get(handle.0)
                .copied()
                .ok_or(WgpuInitError::InvalidHandle(handle.0))?;
            let element = std::mem::size_of::<T>() as u64;
            if element == 0 || size % element != 0 {
                return Err(WgpuInitError::MapFailed(format!(
                    "a buffer of {size} bytes does not hold a whole number of \
                     {element}-byte elements"
                )));
            }

            // Stage the copy and request the mapping inside error scopes, so a
            // rejected copy or map comes back as an error instead of an abort.
            let (sender, receiver) = std::sync::mpsc::channel();
            let staging = self.scoped(|| {
                let staging = self.device.create_buffer(&wgpu::BufferDescriptor {
                    label: Some("oxiphysics_staging_readback"),
                    size,
                    usage: wgpu::BufferUsages::MAP_READ | wgpu::BufferUsages::COPY_DST,
                    mapped_at_creation: false,
                });
                let mut encoder = self
                    .device
                    .create_command_encoder(&wgpu::CommandEncoderDescriptor { label: None });
                encoder.copy_buffer_to_buffer(buffer, 0, &staging, 0, size);
                self.queue.submit(std::iter::once(encoder.finish()));
                staging
                    .slice(..)
                    .map_async(wgpu::MapMode::Read, move |result| {
                        let _ = sender.send(result);
                    });
                Ok(staging)
            })?;

            // Block until the GPU has completed and the mapping is ready.
            self.device
                .poll(wgpu::PollType::Wait {
                    submission_index: None,
                    timeout: None,
                })
                .map_err(|e| WgpuInitError::PollFailed(e.to_string()))?;
            match receiver.recv_timeout(MAP_CALLBACK_TIMEOUT) {
                Ok(Ok(())) => {}
                Ok(Err(e)) => return Err(WgpuInitError::MapFailed(e.to_string())),
                Err(e) => {
                    return Err(WgpuInitError::MapFailed(format!(
                        "the map callback did not run: {e}"
                    )));
                }
            }

            let mapped = staging
                .slice(..)
                .get_mapped_range()
                .map_err(|e| WgpuInitError::MapFailed(e.to_string()))?;
            // Copy through a byte view of a freshly allocated, correctly aligned
            // `Vec<T>`: the mapped pointer's alignment is not guaranteed.
            let mut out = vec![<T as bytemuck::Zeroable>::zeroed(); (size / element) as usize];
            bytemuck::cast_slice_mut::<T, u8>(&mut out).copy_from_slice(&mapped);
            drop(mapped);
            staging.unmap();
            Ok(out)
        }

        // ── Dispatch ──────────────────────────────────────────────────────────

        /// Compute the 3-D workgroup dispatch counts for `n_items` elements.
        ///
        /// Returns `[0, 1, 1]` for `n_items == 0` (no-op dispatch).  Inputs that
        /// need more than 65535 workgroups are folded into the y dimension; the
        /// kernel must then linearise its indices with `num_workgroups` - see
        /// [`dispatch_count_for`](crate::compute::timestamp::dispatch_count_for).
        ///
        /// # Errors
        ///
        /// [`crate::Error::DispatchLimitExceeded`] if `n_items` cannot be
        /// covered even by the folded grid.
        pub fn dispatch_count_for(n_items: usize, workgroup_size: u32) -> crate::Result<[u32; 3]> {
            crate::compute::timestamp::dispatch_count_for(n_items, workgroup_size)
        }

        /// Compile `entry_point` of `wgsl_src`, or fetch it from the pipeline cache.
        ///
        /// A shader that fails to compile is reported as an error and is **not**
        /// cached, so a later corrected call is not shadowed by a broken pipeline.
        fn compiled_pipeline(
            &self,
            wgsl_src: &str,
            entry_point: &str,
        ) -> Result<Arc<wgpu::ComputePipeline>, WgpuInitError> {
            // Hash the shader source + entry point to key the pipeline cache.
            let mut hasher = DefaultHasher::new();
            wgsl_src.hash(&mut hasher);
            entry_point.hash(&mut hasher);
            let key = hasher.finish();

            let mut cache = self
                .shader_cache
                .lock()
                .unwrap_or_else(PoisonError::into_inner);
            if let Some(entry) = cache.get(&key) {
                return Ok(entry.pipeline.clone());
            }
            let pipeline = self.scoped(|| {
                let module = self
                    .device
                    .create_shader_module(wgpu::ShaderModuleDescriptor {
                        label: Some(entry_point),
                        source: wgpu::ShaderSource::Wgsl(wgsl_src.into()),
                    });
                Ok(Arc::new(self.device.create_compute_pipeline(
                    &wgpu::ComputePipelineDescriptor {
                        label: Some(entry_point),
                        layout: None,
                        module: &module,
                        entry_point: Some(entry_point),
                        compilation_options: wgpu::PipelineCompilationOptions::default(),
                        cache: None,
                    },
                )))
            })?;
            cache.insert(
                key,
                ShaderCacheEntry {
                    pipeline: pipeline.clone(),
                },
            );
            Ok(pipeline)
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
        /// * `workgroups`  — `[x, y, z]` dispatch counts, each at most
        ///   [`max_workgroups_per_dimension`](Self::max_workgroups_per_dimension).
        ///
        /// # Errors
        ///
        /// * [`WgpuInitError::InvalidHandle`] if any buffer handle is unknown or
        ///   destroyed.
        /// * [`WgpuInitError::Validation`] if the shader does not compile, the
        ///   bindings do not match the shader, a dimension of `workgroups` exceeds
        ///   the device limit, or an earlier operation left an error behind (see
        ///   [`take_pending_error`](Self::take_pending_error)).  wgpu would abort
        ///   the process for each of these; they are returned instead.
        /// * [`WgpuInitError::PollFailed`] if waiting for the GPU failed.
        pub fn dispatch_wgsl(
            &self,
            wgsl_src: &str,
            entry_point: &str,
            buffers: &[(WgpuBufferHandle, wgpu::BufferBindingType)],
            workgroups: [u32; 3],
        ) -> Result<(), WgpuInitError> {
            self.check_pending_error()?;

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

            let limit = self.max_workgroups_per_dimension;
            if workgroups.iter().any(|&count| count > limit) {
                return Err(WgpuInitError::Validation(format!(
                    "a dispatch of {}x{}x{} workgroups exceeds the device limit of {limit} \
                     per dimension",
                    workgroups[0], workgroups[1], workgroups[2]
                )));
            }

            let pipeline = self.compiled_pipeline(wgsl_src, entry_point)?;

            // Record and submit the compute pass.
            self.scoped(|| {
                // Derive the bind-group layout from the compiled pipeline.
                let bg_layout = pipeline.get_bind_group_layout(0);
                let bind_group = self.device.create_bind_group(&wgpu::BindGroupDescriptor {
                    label: None,
                    layout: &bg_layout,
                    entries: &entries,
                });
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
                Ok(())
            })?;

            // Block until the GPU has finished (synchronous dispatch).
            self.device
                .poll(wgpu::PollType::Wait {
                    submission_index: None,
                    timeout: None,
                })
                .map_err(|e| WgpuInitError::PollFailed(e.to_string()))?;

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
            let out = backend.read_buffer_f64(handle).expect("readback");
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
            let out = backend.read_buffer_f64(out_handle).expect("readback");
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

        /// Adds one to every element (entry `add_one`) or doubles it (`double`).
        const ARITHMETIC_WGSL: &str = r#"
@group(0) @binding(0) var<storage, read_write> data: array<f32>;

@compute @workgroup_size(64)
fn add_one(@builtin(global_invocation_id) gid: vec3<u32>) {
    if gid.x < arrayLength(&data) {
        data[gid.x] = data[gid.x] + 1.0;
    }
}

@compute @workgroup_size(64)
fn double(@builtin(global_invocation_id) gid: vec3<u32>) {
    if gid.x < arrayLength(&data) {
        data[gid.x] = data[gid.x] * 2.0;
    }
}
"#;

        fn read_write() -> wgpu::BufferBindingType {
            wgpu::BufferBindingType::Storage { read_only: false }
        }

        #[test]
        fn destroy_buffer_frees_the_slot_and_is_idempotent() {
            let Some(mut backend) = try_backend() else {
                return;
            };
            let a = backend.create_buffer_storage(64);
            let b = backend.create_buffer_storage(128);
            assert_eq!(backend.live_buffer_count(), 2);
            assert_eq!(backend.live_buffer_bytes(), 192);
            assert_eq!(backend.buffer_size(a), Some(64));

            assert!(backend.destroy_buffer(a), "a live buffer is destroyed");
            assert!(!backend.destroy_buffer(a), "destroying twice is a no-op");
            assert!(
                !backend.destroy_buffer(WgpuBufferHandle(99)),
                "an unknown handle is not destroyed"
            );
            assert_eq!(backend.live_buffer_count(), 1);
            assert_eq!(backend.live_buffer_bytes(), 128);
            assert_eq!(backend.buffer_size(a), None);
            assert_eq!(backend.buffer_size(b), Some(128));

            // The destroyed handle is rejected instead of reaching the GPU.
            let result =
                backend.dispatch_wgsl(ARITHMETIC_WGSL, "add_one", &[(a, read_write())], [1, 1, 1]);
            assert!(matches!(result, Err(WgpuInitError::InvalidHandle(_))));
            assert!(matches!(
                backend.read_buffer_f32(a),
                Err(WgpuInitError::InvalidHandle(_))
            ));
            // The surviving buffer is unaffected.
            let result =
                backend.dispatch_wgsl(ARITHMETIC_WGSL, "add_one", &[(b, read_write())], [1, 1, 1]);
            assert!(result.is_ok());
        }

        #[test]
        fn create_buffer_storage_reuses_destroyed_slots() {
            let Some(mut backend) = try_backend() else {
                return;
            };
            let a = backend.create_buffer_storage(16);
            let b = backend.create_buffer_storage(32);
            let c = backend.create_buffer_storage(48);
            assert_eq!((a.0, b.0, c.0), (0, 1, 2));
            backend.destroy_buffer(b);
            let d = backend.create_buffer_storage(80);
            assert_eq!(d.0, b.0, "the freed slot must be reused");
            assert_eq!(backend.buffer_size(d), Some(80));
            assert_eq!(backend.live_buffer_count(), 3);

            let e = backend.create_buffer_storage(16);
            assert_eq!(e.0, 3, "the table grows only when no slot is free");

            // 1000 create/destroy cycles must not grow the slot table or leak.
            for _ in 0..1000 {
                let h = backend.create_buffer_storage(256);
                assert!(backend.destroy_buffer(h));
            }
            assert_eq!(backend.live_buffer_count(), 4);
            assert_eq!(backend.live_buffer_bytes(), 16 + 80 + 48 + 16);
            let f = backend.create_buffer_storage(16);
            assert_eq!(f.0, 4, "churn used a single extra slot");
        }

        #[test]
        fn ensure_buffer_grows_geometrically_and_keeps_the_buffer_count() {
            let Some(mut backend) = try_backend() else {
                return;
            };
            let mut handle = backend.create_buffer_storage(64);
            let count = backend.live_buffer_count();
            assert!(
                !backend.ensure_buffer(&mut handle, 64),
                "exactly large enough"
            );
            assert!(
                !backend.ensure_buffer(&mut handle, 10),
                "larger than needed"
            );
            assert_eq!(backend.buffer_size(handle), Some(64));

            let before = handle;
            assert!(backend.ensure_buffer(&mut handle, 65));
            assert_eq!(backend.buffer_size(handle), Some(128), "capacity doubles");
            assert_eq!(
                backend.buffer_size(before),
                None,
                "the old buffer is destroyed"
            );
            assert_eq!(backend.live_buffer_count(), count);

            assert!(backend.ensure_buffer(&mut handle, 1000));
            assert_eq!(backend.buffer_size(handle), Some(1000));
            assert!(backend.ensure_buffer(&mut handle, 1001));
            assert_eq!(backend.buffer_size(handle), Some(2000));
            assert_eq!(backend.live_buffer_count(), count);
            assert_eq!(backend.live_buffer_bytes(), 2000);

            // A destroyed handle is replaced by a fresh allocation of at least 16 bytes
            // rounded up to a multiple of 4.
            let gone = backend.create_buffer_storage(8);
            backend.destroy_buffer(gone);
            let mut replaced = gone;
            assert!(backend.ensure_buffer(&mut replaced, 3));
            assert_eq!(backend.buffer_size(replaced), Some(16));
            assert!(backend.ensure_buffer(&mut replaced, 17));
            assert_eq!(backend.buffer_size(replaced), Some(32));

            // The grown buffer is fully usable.
            backend.queue_write_buffer_f32(&handle, &[1.0, 2.0, 3.0]);
            let out = backend.read_buffer_f32(handle).expect("readback");
            assert_eq!(out.len(), 500);
            assert_eq!(&out[..3], &[1.0, 2.0, 3.0]);
        }

        #[test]
        fn pipelines_are_compiled_once_per_shader_and_entry_point() {
            let Some(mut backend) = try_backend() else {
                return;
            };
            assert_eq!(backend.compiled_pipeline_count(), 0);
            let data = backend.create_buffer_storage(16);
            backend.queue_write_buffer_f32(&data, &[1.0, 2.0, 3.0, 4.0]);
            for _ in 0..3 {
                backend
                    .dispatch_wgsl(
                        ARITHMETIC_WGSL,
                        "add_one",
                        &[(data, read_write())],
                        [1, 1, 1],
                    )
                    .expect("add_one dispatch");
            }
            assert_eq!(
                backend.compiled_pipeline_count(),
                1,
                "repeat dispatches reuse the pipeline"
            );
            backend
                .dispatch_wgsl(
                    ARITHMETIC_WGSL,
                    "double",
                    &[(data, read_write())],
                    [1, 1, 1],
                )
                .expect("double dispatch");
            assert_eq!(
                backend.compiled_pipeline_count(),
                2,
                "another entry point is another pipeline"
            );
            // ((x + 1) + 1 + 1) * 2
            assert_eq!(
                backend.read_buffer_f32(data).expect("readback"),
                vec![8.0, 10.0, 12.0, 14.0]
            );
        }

        #[test]
        fn dispatch_count_for_zero_items() {
            assert_eq!(
                WgpuBackendReal::dispatch_count_for(0, 64).unwrap(),
                [0, 1, 1]
            );
        }

        #[test]
        fn dispatch_count_for_65_items() {
            assert_eq!(
                WgpuBackendReal::dispatch_count_for(65, 64).unwrap(),
                [2, 1, 1]
            );
        }

        #[test]
        fn dispatch_count_for_exact_workgroup() {
            assert_eq!(
                WgpuBackendReal::dispatch_count_for(256, 64).unwrap(),
                [4, 1, 1]
            );
        }

        #[test]
        fn dispatch_count_for_folds_large_inputs_and_reports_impossible_ones() {
            let [x, y, z] = WgpuBackendReal::dispatch_count_for(5_000_000, 64).unwrap();
            assert!(x <= 65_535 && y > 1 && y <= 65_535 && z == 1);
            assert!(WgpuBackendReal::dispatch_count_for(usize::MAX, 64).is_err());
        }

        // ── Errors instead of aborts (F1890, F1893) ───────────────────────────

        /// Writes `data[i] = data[i] + 1` with the folded-dispatch index formula
        /// (`idx = gid.x + gid.y * nwg.x * 64`).
        const FOLDED_INCREMENT_WGSL: &str = r#"
@group(0) @binding(0) var<storage, read_write> data: array<u32>;

@compute @workgroup_size(64)
fn increment(
    @builtin(global_invocation_id) gid: vec3<u32>,
    @builtin(num_workgroups) nwg: vec3<u32>,
) {
    let idx = gid.x + gid.y * nwg.x * 64u;
    if idx < arrayLength(&data) {
        data[idx] = data[idx] + 1u;
    }
}
"#;

        #[test]
        fn invalid_wgsl_is_an_error_and_is_not_cached() {
            let Some(mut backend) = try_backend() else {
                return;
            };
            let data = backend.create_buffer_storage(16);
            let broken = "@compute @workgroup_size(1)\nfn main() { let x: u32 = ; }\n";
            let result = backend.dispatch_wgsl(broken, "main", &[(data, read_write())], [1, 1, 1]);
            assert!(
                matches!(result, Err(WgpuInitError::Validation(_))),
                "invalid WGSL must be reported, got {result:?}"
            );
            assert_eq!(
                backend.compiled_pipeline_count(),
                0,
                "a pipeline that failed to compile must not be cached"
            );
            // The same call keeps failing (it is not served from a broken cache entry) ...
            let again = backend.dispatch_wgsl(broken, "main", &[(data, read_write())], [1, 1, 1]);
            assert!(matches!(again, Err(WgpuInitError::Validation(_))));
            // ... and the backend still works afterwards.
            backend.queue_write_buffer_f32(&data, &[1.0, 2.0, 3.0, 4.0]);
            backend
                .dispatch_wgsl(
                    ARITHMETIC_WGSL,
                    "double",
                    &[(data, read_write())],
                    [1, 1, 1],
                )
                .expect("valid dispatch after a failed one");
            assert_eq!(
                backend.read_buffer_f32(data).expect("readback"),
                vec![2.0, 4.0, 6.0, 8.0]
            );
        }

        #[test]
        fn unknown_entry_point_is_an_error() {
            let Some(mut backend) = try_backend() else {
                return;
            };
            let data = backend.create_buffer_storage(16);
            let result = backend.dispatch_wgsl(
                ARITHMETIC_WGSL,
                "no_such_kernel",
                &[(data, read_write())],
                [1, 1, 1],
            );
            assert!(
                matches!(result, Err(WgpuInitError::Validation(_))),
                "{result:?}"
            );
            assert_eq!(backend.compiled_pipeline_count(), 0);
        }

        #[test]
        fn binding_count_mismatch_is_an_error() {
            let Some(mut backend) = try_backend() else {
                return;
            };
            let a = backend.create_buffer_storage(16);
            let b = backend.create_buffer_storage(16);
            // The shader declares one binding; passing two (or none) is a mismatch.
            let too_many = backend.dispatch_wgsl(
                ARITHMETIC_WGSL,
                "add_one",
                &[(a, read_write()), (b, read_write())],
                [1, 1, 1],
            );
            assert!(
                matches!(too_many, Err(WgpuInitError::Validation(_))),
                "{too_many:?}"
            );
            let none = backend.dispatch_wgsl(ARITHMETIC_WGSL, "add_one", &[], [1, 1, 1]);
            assert!(
                matches!(none, Err(WgpuInitError::Validation(_))),
                "{none:?}"
            );
            // A correct call right after still succeeds.
            backend
                .dispatch_wgsl(ARITHMETIC_WGSL, "add_one", &[(a, read_write())], [1, 1, 1])
                .expect("matching bindings");
        }

        #[test]
        fn dispatch_beyond_the_workgroup_limit_is_an_error() {
            let Some(mut backend) = try_backend() else {
                return;
            };
            let data = backend.create_buffer_storage(16);
            let limit = backend.max_workgroups_per_dimension();
            assert!(
                limit >= 65_535,
                "WebGPU guarantees at least 65535 workgroups"
            );
            for shape in [[limit + 1, 1, 1], [1, limit + 1, 1], [1, 1, limit + 1]] {
                let result = backend.dispatch_wgsl(
                    ARITHMETIC_WGSL,
                    "add_one",
                    &[(data, read_write())],
                    shape,
                );
                match result {
                    Err(WgpuInitError::Validation(message)) => {
                        assert!(message.contains("exceeds"), "{message}");
                    }
                    other => panic!("shape {shape:?}: expected a validation error, got {other:?}"),
                }
            }
            // Zero-sized dispatches are legal no-ops.
            backend
                .dispatch_wgsl(
                    ARITHMETIC_WGSL,
                    "add_one",
                    &[(data, read_write())],
                    [0, 1, 1],
                )
                .expect("empty dispatch");
        }

        #[test]
        fn reading_an_unknown_handle_is_an_error_not_an_empty_vec() {
            let Some(backend) = try_backend() else {
                return;
            };
            let missing = WgpuBufferHandle(41);
            assert!(matches!(
                backend.read_buffer_f32(missing),
                Err(WgpuInitError::InvalidHandle(41))
            ));
            assert!(matches!(
                backend.read_buffer_u32(missing),
                Err(WgpuInitError::InvalidHandle(41))
            ));
            assert!(matches!(
                backend.read_buffer_f64(missing),
                Err(WgpuInitError::InvalidHandle(41))
            ));
        }

        #[test]
        fn u32_readback_preserves_every_bit_pattern() {
            let Some(mut backend) = try_backend() else {
                return;
            };
            let words: Vec<u32> = vec![0, 1, 0x7FC0_0001, 0xFFFF_FFFF, 0x8000_0000, 0xDEAD_BEEF];
            let buffer = backend.create_buffer_storage((words.len() * 4) as u64);
            backend.queue_write_buffer_raw(&buffer, bytemuck::cast_slice(&words));
            assert_eq!(backend.read_buffer_u32(buffer).expect("readback"), words);
            // The same bytes read as f32 keep their length.
            assert_eq!(
                backend.read_buffer_f32(buffer).expect("readback").len(),
                words.len()
            );
        }

        #[test]
        fn buffer_sizes_are_rounded_up_to_a_multiple_of_four() {
            let Some(mut backend) = try_backend() else {
                return;
            };
            for (requested, expected) in [(0u64, 4u64), (1, 4), (4, 4), (5, 8), (8, 8), (15, 16)] {
                let handle = backend.create_buffer_storage(requested);
                assert_eq!(
                    backend.buffer_size(handle),
                    Some(expected),
                    "requested {requested} bytes"
                );
            }
            // An "empty" buffer is a valid dispatch target and readback source.
            let empty = backend.create_buffer_storage(0);
            backend
                .dispatch_wgsl(
                    ARITHMETIC_WGSL,
                    "add_one",
                    &[(empty, read_write())],
                    [1, 1, 1],
                )
                .expect("dispatch on a minimum-sized buffer");
            assert_eq!(backend.read_buffer_f32(empty).expect("readback"), vec![1.0]);
        }

        #[test]
        fn oversized_queue_write_is_reported_by_the_next_call_instead_of_aborting() {
            let Some(mut backend) = try_backend() else {
                return;
            };
            let data = backend.create_buffer_storage(16);
            // 32 bytes into a 16-byte buffer: a validation error that aborts the
            // process under wgpu's default handler.
            backend.queue_write_buffer_raw(&data, &[0u8; 32]);
            let result = backend.dispatch_wgsl(
                ARITHMETIC_WGSL,
                "add_one",
                &[(data, read_write())],
                [1, 1, 1],
            );
            assert!(
                matches!(result, Err(WgpuInitError::Validation(_))),
                "the failed write must surface, got {result:?}"
            );
            // The error was consumed: the backend is usable again.
            backend
                .dispatch_wgsl(
                    ARITHMETIC_WGSL,
                    "add_one",
                    &[(data, read_write())],
                    [1, 1, 1],
                )
                .expect("dispatch after the reported error");
        }

        #[test]
        fn take_pending_error_returns_and_clears_the_first_error() {
            let Some(mut backend) = try_backend() else {
                return;
            };
            assert_eq!(backend.take_pending_error(), None);
            let data = backend.create_buffer_storage(16);
            backend.queue_write_buffer_raw(&data, &[0u8; 64]);
            let first = backend.take_pending_error();
            assert!(first.is_some(), "the oversized write must be recorded");
            assert_eq!(
                backend.take_pending_error(),
                None,
                "taking clears the error"
            );
        }

        #[test]
        fn unaligned_queue_write_does_not_abort() {
            let Some(mut backend) = try_backend() else {
                return;
            };
            let data = backend.create_buffer_storage(16);
            backend.queue_write_buffer_raw(&data, &[1u8, 2, 3]);
            assert!(matches!(
                backend.read_buffer_u32(data),
                Err(WgpuInitError::Validation(_))
            ));
            assert!(backend.read_buffer_u32(data).is_ok());
        }

        #[test]
        fn buffer_larger_than_the_device_limit_is_reported_not_aborted() {
            let Some(mut backend) = try_backend() else {
                return;
            };
            let handle = backend.create_buffer_storage(u64::MAX / 4);
            assert!(matches!(
                backend.read_buffer_u32(handle),
                Err(WgpuInitError::Validation(_))
            ));
            backend.destroy_buffer(handle);
            let fine = backend.create_buffer_storage(16);
            assert!(backend.read_buffer_u32(fine).is_ok());
        }

        // ── Folded dispatches (F1880) ─────────────────────────────────────────

        /// Run `increment` over `n` zeroed words with the given dispatch shape
        /// and return the buffer.
        fn run_increment(backend: &mut WgpuBackendReal, n: usize, shape: [u32; 3]) -> Vec<u32> {
            let data = backend.create_buffer_storage((n * 4) as u64);
            backend.queue_write_buffer_raw(&data, &vec![0u8; n * 4]);
            backend
                .dispatch_wgsl(
                    FOLDED_INCREMENT_WGSL,
                    "increment",
                    &[(data, read_write())],
                    shape,
                )
                .expect("increment dispatch");
            let out = backend.read_buffer_u32(data).expect("readback");
            backend.destroy_buffer(data);
            out
        }

        #[test]
        fn folded_dispatch_touches_every_element_exactly_once() {
            let Some(mut backend) = try_backend() else {
                return;
            };
            // A tiny artificial limit forces a genuine 2-D fold on 3000 elements.
            let shape = crate::compute::timestamp::dispatch_count_with_limit(3000, 64, 7).unwrap();
            assert_eq!(shape, [7, 7, 1], "47 workgroups folded into a 7 x 7 grid");
            let out = run_increment(&mut backend, 3000, shape);
            assert_eq!(out.len(), 3000);
            assert!(
                out.iter().all(|&v| v == 1),
                "every element must be processed exactly once"
            );
        }

        #[test]
        fn dispatch_of_five_million_elements_is_folded_and_complete() {
            let Some(mut backend) = try_backend() else {
                return;
            };
            let n = 5_000_000usize;
            let shape = WgpuBackendReal::dispatch_count_for(n, 64).unwrap();
            assert!(
                shape[1] > 1,
                "5e6 / 64 = 78125 workgroups exceeds one dimension"
            );
            let out = run_increment(&mut backend, n, shape);
            assert_eq!(out.len(), n);
            let wrong = out.iter().filter(|&&v| v != 1).count();
            assert_eq!(wrong, 0, "{wrong} of {n} elements were skipped or repeated");
        }
    }
}
