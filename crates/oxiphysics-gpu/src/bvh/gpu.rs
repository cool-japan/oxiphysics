// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! GPU-accelerated BVH traversal using a persistent `WgpuBackendReal` instance.

use super::cpu::{Bvh, flatten, ray_aabb_t};
use super::types::{FlatBvhNode, GpuRay};

/// WGSL source for the BVH traversal kernel.
#[cfg(feature = "wgpu-backend")]
const BVH_TRAVERSAL_WGSL: &str = include_str!("../shaders/bvh_traversal.wgsl");

// ============================================================================
// BvhGpuState — per-BVH GPU resources (allocated once at construction)
// ============================================================================

/// Most rays traced by one dispatch.
///
/// Larger batches are split: one workgroup handles 64 rays and a dispatch may use
/// at most 65 535 workgroups per dimension, and the ray buffer stays far below the
/// device's storage-binding limit.
#[cfg(feature = "wgpu-backend")]
const MAX_RAYS_PER_DISPATCH: usize = 1 << 20;

/// Rays the ray and result buffers hold before the first growth.
#[cfg(feature = "wgpu-backend")]
const INITIAL_RAY_CAPACITY: usize = 256;

/// The device and every GPU buffer of one BVH traverser.
///
/// The node, primitive-AABB, primitive-index and object-ID buffers are uploaded
/// once at construction.  The 16-byte parameter buffer is created once and
/// rewritten for every dispatch.  The ray and result buffers are grow-only: they
/// are replaced (old buffer destroyed, capacity at least doubled) only when a batch
/// holds more rays than any earlier one, so the number of live buffers never
/// changes after construction and a traversal allocates nothing in steady state.
#[cfg(feature = "wgpu-backend")]
pub(crate) struct BvhGpuDevice {
    /// The long-lived wgpu device, queue and pipeline cache.
    pub(crate) backend: crate::compute::wgpu_backend::real::WgpuBackendReal,
    /// BVH nodes: 8 words each (`min.xyz`, `max.xyz`, `left_first`, `count`).
    pub(crate) nodes_buf: crate::compute::WgpuBufferHandle,
    /// Primitive AABB buffer: 6 × f32 per primitive [min_xyz, max_xyz].
    pub(crate) prim_aabbs_buf: crate::compute::WgpuBufferHandle,
    /// Primitive-index buffer: one u32 per entry in the flat prim-index array.
    pub(crate) prim_indices_buf: crate::compute::WgpuBufferHandle,
    /// Object-ID buffer: one i32 per primitive.
    pub(crate) object_ids_buf: crate::compute::WgpuBufferHandle,
    /// `[n_nodes, n_rays, n_prims, 0]`, rewritten for every dispatch.
    pub(crate) params_buf: crate::compute::WgpuBufferHandle,
    /// Ray buffer: 8 × f32 per ray, grow-only.
    pub(crate) rays_buf: crate::compute::WgpuBufferHandle,
    /// Hit buffer: one i32 per ray, grow-only.
    pub(crate) results_buf: crate::compute::WgpuBufferHandle,
    /// Host staging for the ray upload, reused across calls.
    pub(crate) ray_words: Vec<f32>,
}

/// Persistent GPU resources for a BVH.
///
/// `BvhGpuState` keeps the [`BvhGpuDevice`] (the `WgpuBackendReal` and all of its
/// buffers) behind a `Mutex`, because `traverse_rays` only takes `&self`.
///
/// # Thread safety
///
/// `WgpuBackendReal` is `Send + Sync` (device/queue are `Arc`-wrapped and the
/// shader cache uses `Mutex`).  The additional `Mutex` serialises dispatches and the
/// resizing of the grow-only buffers, so callers who share a `BvhGpuTraverser`
/// across threads can do so safely.  `BvhGpuTraverser` is never placed inside a
/// `rayon::ParallelIterator` closure in the current code-base; the lock overhead is
/// negligible compared to GPU dispatch latency.
#[cfg(feature = "wgpu-backend")]
pub(crate) struct BvhGpuState {
    /// Device and buffers; locked per dispatch.
    pub(crate) device: std::sync::Mutex<BvhGpuDevice>,
    /// Number of dispatches completed (observability / reuse test).
    pub(crate) dispatch_count: std::sync::atomic::AtomicU64,
    /// Monotonically increasing ID assigned at construction (reuse test).
    pub(crate) creation_id: u64,
}

// ============================================================================
// BvhTraverserInner — enum over CPU-only and GPU variants
// ============================================================================

pub(crate) enum BvhTraverserInner {
    /// CPU-only fallback.
    Cpu,
    /// GPU backend with pre-uploaded primitive buffers.
    #[cfg(feature = "wgpu-backend")]
    Gpu(Box<BvhGpuState>),
}

// ============================================================================
// BvhGpuTraverser
// ============================================================================

/// GPU-accelerated BVH ray traversal.
///
/// Encodes a flat BVH into GPU-resident buffers and dispatches the
/// `bvh_traversal.wgsl` kernel.  Falls back to CPU traversal when no GPU
/// adapter is available.
///
/// # Usage
///
/// ```rust
/// use oxiphysics_gpu::bvh::{Aabb, Bvh, BvhPrimitive, BvhGpuTraverser, GpuRay};
///
/// let prims = vec![
///     BvhPrimitive::new(Aabb::new([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]), 0),
/// ];
/// let bvh = Bvh::build(prims);
/// let traverser = BvhGpuTraverser::new(&bvh);
/// let rays = vec![GpuRay::new([0.5, 0.5, -1.0], [0.0, 0.0, 1.0], 100.0)];
/// let hits = traverser.traverse_rays(&rays);
/// assert_eq!(hits.len(), 1);
/// ```
pub struct BvhGpuTraverser {
    /// Flat BVH nodes (CPU copy kept for fallback).
    pub(crate) flat_nodes: Vec<FlatBvhNode>,
    /// Primitive indices (CPU copy).
    pub(crate) prim_indices: Vec<usize>,
    /// Primitives (CPU copy).
    pub(crate) primitives: Vec<super::types::BvhPrimitive>,
    /// GPU resources (gated behind feature).
    pub(crate) inner: BvhTraverserInner,
}

/// Monotonically increasing counter for creation_id assignment.
#[cfg(feature = "wgpu-backend")]
static CREATION_COUNTER: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(1);

impl BvhGpuTraverser {
    /// Create a traverser from a BVH.
    ///
    /// Flattens the BVH and uploads the node + primitive buffers to the GPU.
    /// Falls back to CPU traversal if no GPU adapter is found.
    pub fn new(bvh: &Bvh) -> Self {
        let (flat_nodes, prim_indices) = flatten(bvh);
        let primitives = bvh.primitives.clone();

        #[cfg(feature = "wgpu-backend")]
        {
            use crate::compute::wgpu_backend::real::WgpuBackendReal;
            use std::sync::atomic::Ordering;

            if let Ok(mut backend) = WgpuBackendReal::try_new() {
                // ── Upload the BVH once; every later traverse_rays call reuses it ──

                // nodes: 8 words per node [min.xyz, max.xyz, left_first, count]
                let node_words: Vec<u32> = flat_nodes
                    .iter()
                    .flat_map(|n| {
                        [
                            n.aabb.min[0].to_bits(),
                            n.aabb.min[1].to_bits(),
                            n.aabb.min[2].to_bits(),
                            n.aabb.max[0].to_bits(),
                            n.aabb.max[1].to_bits(),
                            n.aabb.max[2].to_bits(),
                            n.left_first,
                            n.count,
                        ]
                    })
                    .collect();
                let nodes_buf =
                    backend.create_buffer_storage((node_words.len() * 4).max(16) as u64);
                backend.queue_write_buffer_raw(&nodes_buf, bytemuck::cast_slice(&node_words));

                // prim_aabbs: 6 × f32 per primitive
                let prim_aabb_f32s: Vec<f32> = primitives
                    .iter()
                    .flat_map(|p| {
                        [
                            p.aabb.min[0],
                            p.aabb.min[1],
                            p.aabb.min[2],
                            p.aabb.max[0],
                            p.aabb.max[1],
                            p.aabb.max[2],
                        ]
                    })
                    .collect();
                let prim_aabbs_buf =
                    backend.create_buffer_storage((prim_aabb_f32s.len() * 4).max(16) as u64);
                backend.queue_write_buffer_f32(&prim_aabbs_buf, &prim_aabb_f32s);

                // prim_indices: one u32 per entry
                let prim_u32s: Vec<u32> = prim_indices.iter().map(|&i| i as u32).collect();
                let prim_indices_buf =
                    backend.create_buffer_storage((prim_u32s.len() * 4).max(16) as u64);
                backend.queue_write_buffer_raw(&prim_indices_buf, bytemuck::cast_slice(&prim_u32s));

                // object_ids: one i32 per primitive
                let obj_ids: Vec<i32> = primitives.iter().map(|p| p.object_id as i32).collect();
                let object_ids_buf =
                    backend.create_buffer_storage((obj_ids.len() * 4).max(16) as u64);
                backend.queue_write_buffer_raw(&object_ids_buf, bytemuck::cast_slice(&obj_ids));

                // Per-call buffers, created once: the parameters are rewritten and
                // the ray/result buffers grow on demand.
                let params_buf = backend.create_buffer_storage(16);
                let rays_buf = backend.create_buffer_storage((INITIAL_RAY_CAPACITY * 32) as u64);
                let results_buf = backend.create_buffer_storage((INITIAL_RAY_CAPACITY * 4) as u64);

                let creation_id = CREATION_COUNTER.fetch_add(1, Ordering::Relaxed);

                return Self {
                    flat_nodes,
                    prim_indices,
                    primitives,
                    inner: BvhTraverserInner::Gpu(Box::new(BvhGpuState {
                        device: std::sync::Mutex::new(BvhGpuDevice {
                            backend,
                            nodes_buf,
                            prim_aabbs_buf,
                            prim_indices_buf,
                            object_ids_buf,
                            params_buf,
                            rays_buf,
                            results_buf,
                            ray_words: Vec::new(),
                        }),
                        dispatch_count: std::sync::atomic::AtomicU64::new(0),
                        creation_id,
                    })),
                };
            }
        }

        Self {
            flat_nodes,
            prim_indices,
            primitives,
            inner: BvhTraverserInner::Cpu,
        }
    }

    /// Create a CPU-only traverser (useful for testing without a GPU).
    pub fn new_cpu(bvh: &Bvh) -> Self {
        let (flat_nodes, prim_indices) = flatten(bvh);
        Self {
            flat_nodes,
            prim_indices,
            primitives: bvh.primitives.clone(),
            inner: BvhTraverserInner::Cpu,
        }
    }

    /// Returns `true` if using a real GPU backend.
    pub fn is_gpu(&self) -> bool {
        match &self.inner {
            BvhTraverserInner::Cpu => false,
            #[cfg(feature = "wgpu-backend")]
            BvhTraverserInner::Gpu(_) => true,
        }
    }

    /// Returns the current dispatch count (GPU variant only; always 0 for CPU).
    #[cfg(feature = "wgpu-backend")]
    pub fn dispatch_count(&self) -> u64 {
        match &self.inner {
            BvhTraverserInner::Cpu => 0,
            BvhTraverserInner::Gpu(state) => state
                .dispatch_count
                .load(std::sync::atomic::Ordering::Relaxed),
        }
    }

    /// Returns the creation_id of the underlying GPU state (for reuse tests).
    #[cfg(feature = "wgpu-backend")]
    pub fn creation_id(&self) -> Option<u64> {
        match &self.inner {
            BvhTraverserInner::Cpu => None,
            BvhTraverserInner::Gpu(state) => Some(state.creation_id),
        }
    }

    /// Number of GPU buffers currently alive on the traverser's device, or `None`
    /// for a CPU-only traverser.
    ///
    /// The BVH is uploaded once and the ray/result buffers are reused and grown in
    /// place, so this stays constant for the lifetime of the traverser.
    #[cfg(feature = "wgpu-backend")]
    pub fn gpu_buffer_count(&self) -> Option<usize> {
        match &self.inner {
            BvhTraverserInner::Cpu => None,
            BvhTraverserInner::Gpu(state) => state
                .device
                .lock()
                .ok()
                .map(|device| device.backend.live_buffer_count()),
        }
    }

    /// Total size in bytes of the live GPU buffers of the traverser's device, or
    /// `None` for a CPU-only traverser.
    ///
    /// Grows only when a call traces more rays than any earlier call.
    #[cfg(feature = "wgpu-backend")]
    pub fn gpu_buffer_bytes(&self) -> Option<u64> {
        match &self.inner {
            BvhTraverserInner::Cpu => None,
            BvhTraverserInner::Gpu(state) => state
                .device
                .lock()
                .ok()
                .map(|device| device.backend.live_buffer_bytes()),
        }
    }

    /// Traverse the BVH for each ray.
    ///
    /// Returns a `Vec<i32>` of length `rays.len()`.  Each element is either:
    /// - the `object_id` of the first hit leaf's primitive, or
    /// - `-1` if no intersection was found.
    pub fn traverse_rays(&self, rays: &[GpuRay]) -> Vec<i32> {
        if rays.is_empty() {
            return Vec::new();
        }
        match &self.inner {
            BvhTraverserInner::Cpu => self.traverse_rays_cpu(rays),
            // An empty BVH has nothing to hit; the kernel returns before writing any
            // result in that case, so answer on the host.
            #[cfg(feature = "wgpu-backend")]
            BvhTraverserInner::Gpu(_) if self.flat_nodes.is_empty() => self.traverse_rays_cpu(rays),
            #[cfg(feature = "wgpu-backend")]
            BvhTraverserInner::Gpu(state) => self
                .traverse_rays_gpu(state, rays)
                .unwrap_or_else(|_| self.traverse_rays_cpu(rays)),
        }
    }

    // ── CPU traversal ─────────────────────────────────────────────────────────

    fn traverse_rays_cpu(&self, rays: &[GpuRay]) -> Vec<i32> {
        rays.iter()
            .map(|ray| self.traverse_single_cpu(ray))
            .collect()
    }

    fn traverse_single_cpu(&self, ray: &GpuRay) -> i32 {
        if self.flat_nodes.is_empty() {
            return -1;
        }
        let inv_dir = [
            1.0 / ray.direction[0],
            1.0 / ray.direction[1],
            1.0 / ray.direction[2],
        ];
        let origin = ray.origin;
        let max_t = ray.max_t;
        let mut best_hit: i32 = -1;
        let mut best_t = max_t;

        let mut stack = Vec::with_capacity(64);
        stack.push(0usize);

        while let Some(idx) = stack.pop() {
            let node = &self.flat_nodes[idx];
            // Slab test
            if ray_aabb_t(origin, inv_dir, &node.aabb).is_none() {
                continue;
            }
            if node.count > 0 {
                // Leaf: check each primitive
                let start = node.left_first as usize;
                let end = (start + node.count as usize).min(self.prim_indices.len());
                for &pi in &self.prim_indices[start..end] {
                    if pi >= self.primitives.len() {
                        continue;
                    }
                    if let Some((t_near, _)) =
                        ray_aabb_t(origin, inv_dir, &self.primitives[pi].aabb)
                        && t_near >= 0.0
                        && t_near < best_t
                    {
                        best_t = t_near;
                        best_hit = self.primitives[pi].object_id as i32;
                    }
                }
            } else {
                let right = node.left_first as usize;
                let left = idx + 1;
                if right < self.flat_nodes.len() {
                    stack.push(right);
                }
                if left < self.flat_nodes.len() && left != right {
                    stack.push(left);
                }
            }
        }
        best_hit
    }

    // ── GPU traversal ─────────────────────────────────────────────────────────

    /// Trace `rays` on the GPU in batches of at most [`MAX_RAYS_PER_DISPATCH`].
    ///
    /// Returns `Err` when the device lock is poisoned, a dispatch fails or a
    /// read-back comes back short; `traverse_rays` then answers on the host.
    #[cfg(feature = "wgpu-backend")]
    fn traverse_rays_gpu(
        &self,
        state: &BvhGpuState,
        rays: &[GpuRay],
    ) -> Result<Vec<i32>, crate::GpuError> {
        use std::sync::atomic::Ordering;

        let too_large = |what: &str| {
            crate::GpuError::ShaderDispatch(format!("{what} count does not fit a 32-bit index"))
        };
        let n_nodes = u32::try_from(self.flat_nodes.len()).map_err(|_| too_large("BVH node"))?;
        let n_prims =
            u32::try_from(self.prim_indices.len()).map_err(|_| too_large("BVH primitive"))?;

        let mut device = state.device.lock().map_err(|_| {
            crate::GpuError::ShaderDispatch(
                "BVH GPU state lock poisoned by an earlier panic".into(),
            )
        })?;
        let mut hits = Vec::with_capacity(rays.len());
        for batch in rays.chunks(MAX_RAYS_PER_DISPATCH) {
            hits.extend(device.trace_batch(batch, n_nodes, n_prims)?);
            state.dispatch_count.fetch_add(1, Ordering::Relaxed);
        }
        Ok(hits)
    }
}

#[cfg(feature = "wgpu-backend")]
impl BvhGpuDevice {
    /// Upload `batch` (at most [`MAX_RAYS_PER_DISPATCH`] rays), dispatch the
    /// traversal kernel and read the hit of every ray back.
    ///
    /// The ray and result buffers are grown in place when the batch is larger than
    /// any earlier one; nothing else is allocated on the device.
    fn trace_batch(
        &mut self,
        batch: &[GpuRay],
        n_nodes: u32,
        n_prims: u32,
    ) -> Result<Vec<i32>, crate::GpuError> {
        let n_rays = batch.len();
        let n_rays_u32 = u32::try_from(n_rays).map_err(|_| {
            crate::GpuError::ShaderDispatch("ray batch does not fit a 32-bit index".into())
        })?;

        // Rays: 8 f32 per ray [origin.xyz, direction.xyz, max_t, pad].
        self.ray_words.clear();
        for ray in batch {
            self.ray_words.extend_from_slice(&[
                ray.origin[0],
                ray.origin[1],
                ray.origin[2],
                ray.direction[0],
                ray.direction[1],
                ray.direction[2],
                ray.max_t,
                0.0,
            ]);
        }
        self.backend
            .ensure_buffer(&mut self.rays_buf, (self.ray_words.len() * 4) as u64);
        self.backend
            .ensure_buffer(&mut self.results_buf, (n_rays * 4) as u64);
        self.backend
            .queue_write_buffer_f32(&self.rays_buf, &self.ray_words);
        self.backend.queue_write_buffer_raw(
            &self.params_buf,
            bytemuck::cast_slice(&[n_nodes, n_rays_u32, n_prims, 0_u32]),
        );

        let read_only = wgpu::BufferBindingType::Storage { read_only: true };
        let read_write = wgpu::BufferBindingType::Storage { read_only: false };
        let workgroups = crate::compute::timestamp::dispatch_count_for(n_rays, 64)
            .map_err(|e| crate::GpuError::ShaderDispatch(e.to_string()))?;
        self.backend
            .dispatch_wgsl(
                BVH_TRAVERSAL_WGSL,
                "main",
                &[
                    (self.params_buf, read_only),
                    (self.nodes_buf, read_only),
                    (self.rays_buf, read_only),
                    (self.results_buf, read_write),
                    (self.prim_indices_buf, read_only),
                    (self.object_ids_buf, read_only),
                    (self.prim_aabbs_buf, read_only),
                ],
                workgroups,
            )
            .map_err(|e| crate::GpuError::ShaderDispatch(e.to_string()))?;

        // The grown result buffer may be larger than this batch: keep the first n_rays.
        let raw = self
            .backend
            .read_buffer_u32(self.results_buf)
            .map_err(|e| crate::GpuError::ReadBack(e.to_string()))?;
        if raw.len() < n_rays {
            return Err(crate::GpuError::ReadBack(format!(
                "expected {n_rays} BVH hits, the read-back returned {}",
                raw.len()
            )));
        }
        Ok(raw[..n_rays].iter().map(|&word| word as i32).collect())
    }
}
