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

/// Persistent GPU resources for a BVH.
///
/// `BvhGpuState` owns a `WgpuBackendReal` and three primitive buffers that are
/// uploaded once at construction time.  Per-call allocations are limited to the
/// rays buffer and the results buffer, which are re-allocated on demand when the
/// ray count increases.
///
/// # Thread safety
///
/// `WgpuBackendReal` is `Send + Sync` (device/queue are `Arc`-wrapped and the
/// shader cache uses `Mutex`).  We wrap the backend in an additional `Mutex` so
/// that callers who share a `BvhGpuTraverser` across threads can do so safely
/// without re-entrant dispatch issues.  `BvhGpuTraverser` is never placed inside
/// a `rayon::ParallelIterator` closure in the current code-base (verified by the
/// send-bound audit in Step 0 of the block spec); nevertheless we default to
/// `Mutex` so the type is unconditionally `Send + Sync`.
#[cfg(feature = "wgpu-backend")]
pub(crate) struct BvhGpuState {
    /// Backend stored once; locked per dispatch.
    ///
    /// Using `Mutex<WgpuBackendReal>` rather than `RefCell` so that
    /// `BvhGpuTraverser` is `Send + Sync` regardless of call-site threading.
    /// The audit found no parallel call sites today, but the Mutex overhead is
    /// negligible compared to GPU dispatch latency.
    pub(crate) backend: std::sync::Mutex<crate::compute::wgpu_backend::real::WgpuBackendReal>,
    /// Primitive AABB buffer: 6 × f32 per primitive [min_xyz, max_xyz].
    pub(crate) prim_aabbs_buf: crate::compute::WgpuBufferHandle,
    /// Primitive-index buffer: one u32 per entry in the flat prim-index array.
    pub(crate) prim_indices_buf: crate::compute::WgpuBufferHandle,
    /// Object-ID buffer: one i32 per primitive.
    pub(crate) object_ids_buf: crate::compute::WgpuBufferHandle,
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
                // ── Upload primitive data (once, reused across all traverse_rays calls) ──

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

                let creation_id = CREATION_COUNTER.fetch_add(1, Ordering::Relaxed);

                return Self {
                    flat_nodes,
                    prim_indices,
                    primitives,
                    inner: BvhTraverserInner::Gpu(Box::new(BvhGpuState {
                        backend: std::sync::Mutex::new(backend),
                        prim_aabbs_buf,
                        prim_indices_buf,
                        object_ids_buf,
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

    #[cfg(feature = "wgpu-backend")]
    fn traverse_rays_gpu(
        &self,
        state: &BvhGpuState,
        rays: &[GpuRay],
    ) -> Result<Vec<i32>, crate::GpuError> {
        use std::sync::atomic::Ordering;

        let n_rays = rays.len() as u32;
        let n_nodes = self.flat_nodes.len() as u32;
        let n_prims = self.prim_indices.len() as u32;

        // Lock the backend for this dispatch.
        let mut backend = state
            .backend
            .lock()
            .expect("BvhGpuState backend lock poisoned");

        // ── Encode params [n_nodes, n_rays, n_prims, 0] (per-call) ──────────
        let params_data: [u32; 4] = [n_nodes, n_rays, n_prims, 0];
        let params_buf = backend.create_buffer_storage(16);
        backend.queue_write_buffer_raw(&params_buf, bytemuck::cast_slice(&params_data));

        // ── Encode BVH nodes (per-call, same data each time) ─────────────────
        // Layout: [min_x, min_y, min_z, max_x, max_y, max_z, left_first_bits, count_bits]
        let node_f32s: Vec<f32> = self
            .flat_nodes
            .iter()
            .flat_map(|n| {
                [
                    n.aabb.min[0],
                    n.aabb.min[1],
                    n.aabb.min[2],
                    n.aabb.max[0],
                    n.aabb.max[1],
                    n.aabb.max[2],
                    f32::from_bits(n.left_first),
                    f32::from_bits(n.count),
                ]
            })
            .collect();
        let nodes_buf = backend.create_buffer_storage((node_f32s.len() * 4).max(16) as u64);
        backend.queue_write_buffer_f32(&nodes_buf, &node_f32s);

        // ── Encode rays: 8 f32 per ray ────────────────────────────────────────
        let ray_f32s: Vec<f32> = rays
            .iter()
            .flat_map(|r| {
                [
                    r.origin[0],
                    r.origin[1],
                    r.origin[2],
                    r.direction[0],
                    r.direction[1],
                    r.direction[2],
                    r.max_t,
                    0.0_f32, // pad
                ]
            })
            .collect();
        let rays_buf = backend.create_buffer_storage((ray_f32s.len() * 4).max(16) as u64);
        backend.queue_write_buffer_f32(&rays_buf, &ray_f32s);

        // ── Results buffer ────────────────────────────────────────────────────
        let results_buf = backend.create_buffer_storage((n_rays as usize * 4).max(16) as u64);

        // ── Dispatch ──────────────────────────────────────────────────────────
        let workgroups_x = n_rays.div_ceil(64);
        backend
            .dispatch_wgsl(
                BVH_TRAVERSAL_WGSL,
                "main",
                &[
                    (
                        params_buf,
                        wgpu::BufferBindingType::Storage { read_only: true },
                    ),
                    (
                        nodes_buf,
                        wgpu::BufferBindingType::Storage { read_only: true },
                    ),
                    (
                        rays_buf,
                        wgpu::BufferBindingType::Storage { read_only: true },
                    ),
                    (
                        results_buf,
                        wgpu::BufferBindingType::Storage { read_only: false },
                    ),
                    (
                        state.prim_indices_buf,
                        wgpu::BufferBindingType::Storage { read_only: true },
                    ),
                    (
                        state.object_ids_buf,
                        wgpu::BufferBindingType::Storage { read_only: true },
                    ),
                    (
                        state.prim_aabbs_buf,
                        wgpu::BufferBindingType::Storage { read_only: true },
                    ),
                ],
                [workgroups_x, 1, 1],
            )
            .map_err(|e| crate::GpuError::ShaderDispatch(e.to_string()))?;

        // ── Read results back ─────────────────────────────────────────────────
        let raw = backend.read_buffer_f32(results_buf);
        let hits: Vec<i32> = raw
            .iter()
            .take(n_rays as usize)
            .map(|&f| f32::to_bits(f) as i32)
            .collect();

        // Pad to n_rays if buffer was short.
        let mut result = hits;
        result.resize(n_rays as usize, -1);

        // Increment observability counter.
        state.dispatch_count.fetch_add(1, Ordering::Relaxed);

        Ok(result)
    }
}
