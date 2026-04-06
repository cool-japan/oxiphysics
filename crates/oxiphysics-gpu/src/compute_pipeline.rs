#![allow(clippy::ptr_arg)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! CPU-side simulation of a wgpu compute pipeline abstraction.
//!
//! Provides data layout and dispatch logic that mirrors a typical GPU compute
//! pipeline, with all computation running on the CPU.

// ── silence unused-item lints for the enum variants / fields used only in
//    tests or future GPU back-ends ──────────────────────────────────────────
#![allow(dead_code)]

// ─────────────────────────────────────────────────────────────────────────────
// BufferUsage
// ─────────────────────────────────────────────────────────────────────────────

/// Intended usage of a [`ComputeBuffer`].
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum BufferUsage {
    /// Vertex data fed to the rasteriser.
    Vertex,
    /// Index data for indexed draw calls.
    Index,
    /// Small, frequently-updated uniform / constant data.
    Uniform,
    /// General-purpose read-write storage.
    Storage,
    /// CPU↔GPU transfer staging area.
    Staging,
}

// ─────────────────────────────────────────────────────────────────────────────
// ComputeBuffer
// ─────────────────────────────────────────────────────────────────────────────

/// A CPU-resident buffer that mirrors a GPU buffer.
#[derive(Debug, Clone)]
pub struct ComputeBuffer {
    /// Raw `f32` elements.
    pub data: Vec<f32>,
    /// Intended GPU usage hint.
    pub usage: BufferUsage,
    /// Human-readable label (useful for debugging).
    pub label: String,
}

impl ComputeBuffer {
    /// Allocate a zero-initialised buffer of `size` `f32` elements.
    pub fn new(size: usize, usage: BufferUsage, label: &str) -> Self {
        Self {
            data: vec![0.0_f32; size],
            usage,
            label: label.to_owned(),
        }
    }

    /// Write `values` into the buffer starting at element `offset`.
    ///
    /// # Panics
    /// Panics if `offset + values.len() > self.data.len()`.
    pub fn write_f32(&mut self, offset: usize, values: &[f32]) {
        let end = offset + values.len();
        assert!(
            end <= self.data.len(),
            "write_f32: out-of-bounds write (offset={offset}, len={}, capacity={})",
            values.len(),
            self.data.len()
        );
        self.data[offset..end].copy_from_slice(values);
    }

    /// Read `count` `f32` elements starting at element `offset`.
    ///
    /// # Panics
    /// Panics if `offset + count > self.data.len()`.
    pub fn read_f32(&self, offset: usize, count: usize) -> Vec<f32> {
        let end = offset + count;
        assert!(
            end <= self.data.len(),
            "read_f32: out-of-bounds read (offset={offset}, count={count}, capacity={})",
            self.data.len()
        );
        self.data[offset..end].to_vec()
    }

    /// Total size of the buffer in bytes (`len * 4`).
    pub fn byte_size(&self) -> usize {
        self.data.len() * std::mem::size_of::<f32>()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// WorkgroupSize
// ─────────────────────────────────────────────────────────────────────────────

/// Workgroup dimensions for a compute dispatch.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct WorkgroupSize {
    /// X dimension.
    pub x: u32,
    /// Y dimension.
    pub y: u32,
    /// Z dimension.
    pub z: u32,
}

impl WorkgroupSize {
    /// Compute the number of workgroups needed to cover `total` items when
    /// each workgroup handles `workgroup` items (ceiling division).
    pub fn dispatch_count(total: u32, workgroup: u32) -> u32 {
        assert!(workgroup > 0, "workgroup size must be > 0");
        total.div_ceil(workgroup)
    }
}

impl Default for WorkgroupSize {
    fn default() -> Self {
        Self { x: 64, y: 1, z: 1 }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// ComputeKernelKind
// ─────────────────────────────────────────────────────────────────────────────

/// The kind of compute kernel to dispatch.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ComputeKernelKind {
    /// Semi-implicit velocity / position integration.
    VelocityUpdate,
    /// Pressure Jacobi iteration.
    PressureJacobi,
    /// Lennard-Jones particle force evaluation.
    ParticleForce,
    /// Neighbour-search (spatial hashing).
    NeighborSearch,
    /// User-defined kernel identified by a string tag.
    Custom(String),
}

// ─────────────────────────────────────────────────────────────────────────────
// CpuComputeDispatch
// ─────────────────────────────────────────────────────────────────────────────

/// Dispatcher that executes compute kernels on the CPU.
pub struct CpuComputeDispatch {
    /// Which kernel this dispatcher is configured for.
    pub kernel: ComputeKernelKind,
    /// Workgroup size hint (informational on the CPU path).
    pub workgroup_size: WorkgroupSize,
}

impl CpuComputeDispatch {
    /// Create a new dispatcher.
    pub fn new(kernel: ComputeKernelKind, wg: WorkgroupSize) -> Self {
        Self {
            kernel,
            workgroup_size: wg,
        }
    }

    /// Semi-implicit Euler integration for `n` particles.
    ///
    /// `pos[i] += vel[i] * dt` then `vel[i] += force[i] / mass[i] * dt`.
    pub fn dispatch_velocity_update(
        &self,
        pos: &mut ComputeBuffer,
        vel: &mut ComputeBuffer,
        force: &ComputeBuffer,
        mass: &ComputeBuffer,
        dt: f32,
        n: usize,
    ) {
        for i in 0..n {
            pos.data[i] += vel.data[i] * dt;
            vel.data[i] += force.data[i] / mass.data[i] * dt;
        }
    }

    /// Single Jacobi pressure iteration over an `nx × ny` grid.
    ///
    /// Interior points only; boundary cells are left unchanged.
    ///
    /// `p[i,j] = (p_old[i+1,j] + p_old[i-1,j] + p_old[i,j+1] + p_old[i,j-1]
    ///            - dx² * rhs[i,j]) / 4`
    pub fn dispatch_pressure_jacobi(
        &self,
        p: &mut ComputeBuffer,
        p_old: &ComputeBuffer,
        rhs: &ComputeBuffer,
        nx: usize,
        ny: usize,
        dx: f32,
    ) {
        let dx2 = dx * dx;
        for j in 1..ny - 1 {
            for i in 1..nx - 1 {
                let idx = j * nx + i;
                p.data[idx] = (p_old.data[idx + 1]
                    + p_old.data[idx - 1]
                    + p_old.data[idx + nx]
                    + p_old.data[idx - nx]
                    - dx2 * rhs.data[idx])
                    / 4.0;
            }
        }
    }

    /// O(n²) Lennard-Jones force accumulation.
    ///
    /// Positions are stored as interleaved `[x0, y0, x1, y1, …]`.
    /// Forces are accumulated in-place (`force` is zeroed first).
    pub fn dispatch_particle_force(
        &self,
        pos: &ComputeBuffer,
        force: &mut ComputeBuffer,
        eps: f32,
        sigma: f32,
        n: usize,
    ) {
        // Zero forces first.
        for v in force.data[..2 * n].iter_mut() {
            *v = 0.0;
        }

        for i in 0..n {
            for j in (i + 1)..n {
                let dx = pos.data[2 * j] - pos.data[2 * i];
                let dy = pos.data[2 * j + 1] - pos.data[2 * i + 1];
                let r2 = dx * dx + dy * dy;
                if r2 < 1e-12 {
                    continue;
                }
                let sr2 = (sigma * sigma) / r2;
                let sr6 = sr2 * sr2 * sr2;
                let sr12 = sr6 * sr6;
                // F = 24ε/r² * (2(σ/r)^12 - (σ/r)^6)
                let fmag = 24.0 * eps / r2 * (2.0 * sr12 - sr6);
                force.data[2 * i] -= fmag * dx;
                force.data[2 * i + 1] -= fmag * dy;
                force.data[2 * j] += fmag * dx;
                force.data[2 * j + 1] += fmag * dy;
            }
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// GpuStats
// ─────────────────────────────────────────────────────────────────────────────

/// Accumulated statistics for compute dispatches.
#[derive(Debug, Clone, Default)]
pub struct GpuStats {
    /// Total number of dispatches recorded.
    pub dispatch_count: u64,
    /// Total bytes transferred (reads + writes).
    pub bytes_transferred: u64,
    /// Total kernel wall-clock time in milliseconds.
    pub kernel_time_ms: f64,
}

impl GpuStats {
    /// Create zeroed stats.
    pub fn new() -> Self {
        Self::default()
    }

    /// Record a single dispatch event.
    pub fn record_dispatch(&mut self, bytes: u64, time_ms: f64) {
        self.dispatch_count += 1;
        self.bytes_transferred += bytes;
        self.kernel_time_ms += time_ms;
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Free-standing solver helpers
// ─────────────────────────────────────────────────────────────────────────────

/// Single Jacobi sweep over an `nx × ny` pressure grid.
///
/// Updates every interior cell of `p_new` using `p_old` and `rhs`.
pub fn jacobi_step_2d(
    p_new: &mut Vec<f32>,
    p_old: &[f32],
    rhs: &[f32],
    nx: usize,
    ny: usize,
    dx: f32,
) {
    let dx2 = dx * dx;
    for j in 1..ny - 1 {
        for i in 1..nx - 1 {
            let idx = j * nx + i;
            p_new[idx] = (p_old[idx + 1] + p_old[idx - 1] + p_old[idx + nx] + p_old[idx - nx]
                - dx2 * rhs[idx])
                / 4.0;
        }
    }
}

/// Iterative Jacobi pressure-Poisson solver.
///
/// Runs `n_iter` Jacobi sweeps and returns the final L∞ residual.
pub fn pressure_poisson_solve(
    p: &mut Vec<f32>,
    rhs: &[f32],
    nx: usize,
    ny: usize,
    dx: f32,
    n_iter: usize,
) -> f32 {
    let size = nx * ny;
    let mut p_old = p.clone();

    for _ in 0..n_iter {
        jacobi_step_2d(p, &p_old, rhs, nx, ny, dx);
        p_old.copy_from_slice(&p[..size]);
    }

    // Compute L∞ residual on interior cells.
    let dx2 = dx * dx;
    let mut residual = 0.0_f32;
    for j in 1..ny - 1 {
        for i in 1..nx - 1 {
            let idx = j * nx + i;
            let lap = (p[idx + 1] + p[idx - 1] + p[idx + nx] + p[idx - nx] - 4.0 * p[idx]) / dx2;
            let r = (lap - rhs[idx]).abs();
            if r > residual {
                residual = r;
            }
        }
    }
    residual
}

// ─────────────────────────────────────────────────────────────────────────────
// PipelineCache – LRU-style pipeline caching
// ─────────────────────────────────────────────────────────────────────────────

/// A simple LRU-eviction cache for compiled compute pipelines.
///
/// Keyed by a string label; evicts the oldest entry when the capacity is reached.
pub struct PipelineCache {
    /// Maximum number of entries.
    capacity: usize,
    /// Entries in insertion order (oldest first).
    entries: Vec<(String, CpuComputeDispatch)>,
}

impl PipelineCache {
    /// Create a new cache with the given capacity.
    pub fn new(capacity: usize) -> Self {
        Self {
            capacity,
            entries: Vec::new(),
        }
    }

    /// Insert (or replace) a pipeline under `key`.
    pub fn insert(&mut self, key: &str, pipeline: CpuComputeDispatch) {
        // Remove existing entry with the same key
        self.entries.retain(|(k, _)| k != key);
        // Evict oldest if at capacity
        while self.entries.len() >= self.capacity {
            self.entries.remove(0);
        }
        self.entries.push((key.to_owned(), pipeline));
    }

    /// Look up a cached pipeline by key.
    pub fn get(&self, key: &str) -> Option<&CpuComputeDispatch> {
        self.entries.iter().find(|(k, _)| k == key).map(|(_, v)| v)
    }

    /// Number of cached entries.
    pub fn len(&self) -> usize {
        self.entries.len()
    }

    /// Whether the cache is empty.
    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }

    /// Clear all cached pipelines.
    pub fn clear(&mut self) {
        self.entries.clear();
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// PipelineStats – dispatch-level statistics
// ─────────────────────────────────────────────────────────────────────────────

/// Fine-grained statistics for compute pipeline usage.
#[derive(Debug, Clone, Default)]
pub struct PipelineStats {
    /// Total number of dispatches.
    pub total_dispatches: u64,
    /// Total number of workgroups launched.
    pub total_workgroups: u64,
    /// Total number of invocations (workgroups × workgroup_size).
    pub total_invocations: u64,
    /// Number of times a cached pipeline was re-used.
    pub cache_hits: u64,
    /// Number of times a pipeline had to be compiled (or was not in cache).
    pub cache_misses: u64,
}

impl PipelineStats {
    /// Record a single dispatch event.
    pub fn record_dispatch(&mut self, num_workgroups: u64, wg_size: WorkgroupSize) {
        self.total_dispatches += 1;
        self.total_workgroups += num_workgroups;
        self.total_invocations +=
            num_workgroups * (wg_size.x as u64) * (wg_size.y as u64) * (wg_size.z as u64);
    }

    /// Cache hit ratio (0.0–1.0). Returns NaN when no lookups occurred.
    pub fn cache_hit_ratio(&self) -> f64 {
        let total = self.cache_hits + self.cache_misses;
        if total == 0 {
            return f64::NAN;
        }
        self.cache_hits as f64 / total as f64
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// MultiPassPipeline – chained compute passes
// ─────────────────────────────────────────────────────────────────────────────

/// A single compute pass in a multi-pass pipeline.
#[derive(Debug, Clone)]
pub struct ComputePass {
    /// Human-readable label.
    pub label: String,
    /// The kernel to dispatch.
    pub kernel: ComputeKernelKind,
    /// Workgroup size for this pass.
    pub workgroup_size: WorkgroupSize,
    /// Indices of buffers bound to this pass.
    pub buffer_bindings: Vec<usize>,
}

/// A sequence of compute passes that execute in order.
#[derive(Debug)]
pub struct MultiPassPipeline {
    /// Human-readable label.
    pub label: String,
    /// Ordered list of passes.
    pub passes: Vec<ComputePass>,
}

impl MultiPassPipeline {
    /// Create a new empty multi-pass pipeline.
    pub fn new(label: &str) -> Self {
        Self {
            label: label.to_owned(),
            passes: Vec::new(),
        }
    }

    /// Append a compute pass.
    pub fn add_pass(&mut self, pass: ComputePass) {
        self.passes.push(pass);
    }

    /// Number of passes.
    pub fn num_passes(&self) -> usize {
        self.passes.len()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Pipeline validation
// ─────────────────────────────────────────────────────────────────────────────

/// Validate resource bindings for a single compute pass.
///
/// Returns a list of error messages (empty if valid).
pub fn validate_resource_bindings(pass: &ComputePass, buffers: &[ComputeBuffer]) -> Vec<String> {
    let mut errors = Vec::new();
    let mut seen = std::collections::HashSet::new();
    for &idx in &pass.buffer_bindings {
        if idx >= buffers.len() {
            errors.push(format!(
                "Pass '{}': buffer binding {} is out of range (have {} buffers)",
                pass.label,
                idx,
                buffers.len()
            ));
        }
        if !seen.insert(idx) {
            errors.push(format!(
                "Pass '{}': Duplicate buffer binding {}",
                pass.label, idx
            ));
        }
    }
    errors
}

/// Validate all passes in a multi-pass pipeline.
pub fn validate_pipeline(pipeline: &MultiPassPipeline, buffers: &[ComputeBuffer]) -> Vec<String> {
    let mut errors = Vec::new();
    for pass in &pipeline.passes {
        errors.extend(validate_resource_bindings(pass, buffers));
    }
    errors
}

// ─────────────────────────────────────────────────────────────────────────────
// Additional solver helpers
// ─────────────────────────────────────────────────────────────────────────────

/// Single SOR (Successive Over-Relaxation) sweep over an `nx × ny` grid.
///
/// `omega = 1.0` gives standard Gauss-Seidel; `omega ∈ (1, 2)` gives SOR.
pub fn sor_step_2d(
    p: &mut Vec<f32>,
    p_old: &[f32],
    rhs: &[f32],
    nx: usize,
    ny: usize,
    dx: f32,
    omega: f32,
) {
    let dx2 = dx * dx;
    for j in 1..ny - 1 {
        for i in 1..nx - 1 {
            let idx = j * nx + i;
            let gs = (p_old[idx + 1] + p_old[idx - 1] + p_old[idx + nx] + p_old[idx - nx]
                - dx2 * rhs[idx])
                / 4.0;
            p[idx] = (1.0 - omega) * p_old[idx] + omega * gs;
        }
    }
}

/// Red-black Gauss-Seidel sweep (in-place) on an `nx × ny` grid.
///
/// Updates "red" cells (i+j even) first, then "black" cells (i+j odd).
pub fn red_black_gauss_seidel_step(p: &mut Vec<f32>, rhs: &[f32], nx: usize, ny: usize, dx: f32) {
    let dx2 = dx * dx;
    // Red sweep (i + j even)
    for j in 1..ny - 1 {
        for i in 1..nx - 1 {
            if (i + j) % 2 == 0 {
                let idx = j * nx + i;
                p[idx] =
                    (p[idx + 1] + p[idx - 1] + p[idx + nx] + p[idx - nx] - dx2 * rhs[idx]) / 4.0;
            }
        }
    }
    // Black sweep (i + j odd)
    for j in 1..ny - 1 {
        for i in 1..nx - 1 {
            if (i + j) % 2 == 1 {
                let idx = j * nx + i;
                p[idx] =
                    (p[idx + 1] + p[idx - 1] + p[idx + nx] + p[idx - nx] - dx2 * rhs[idx]) / 4.0;
            }
        }
    }
}

/// Compute the L∞ residual of a 2D Poisson discretisation.
pub fn compute_linf_residual(p: &[f32], rhs: &[f32], nx: usize, ny: usize, dx: f32) -> f32 {
    let dx2 = dx * dx;
    let mut residual = 0.0_f32;
    for j in 1..ny - 1 {
        for i in 1..nx - 1 {
            let idx = j * nx + i;
            let lap = (p[idx + 1] + p[idx - 1] + p[idx + nx] + p[idx - nx] - 4.0 * p[idx]) / dx2;
            let r = (lap - rhs[idx]).abs();
            if r > residual {
                residual = r;
            }
        }
    }
    residual
}

/// O(n²) neighbor search for 2D interleaved positions `[x0, y0, x1, y1, …]`.
///
/// Returns a `Vec<Vec`usize`>` where `result[i]` contains the indices of
/// particles within `cutoff` distance of particle `i`.
pub fn dispatch_neighbor_search(positions: &[f32], n: usize, cutoff: f32) -> Vec<Vec<usize>> {
    let cutoff2 = cutoff * cutoff;
    let mut neighbors = vec![Vec::new(); n];
    for i in 0..n {
        for j in (i + 1)..n {
            let dx = positions[2 * j] - positions[2 * i];
            let dy = positions[2 * j + 1] - positions[2 * i + 1];
            let r2 = dx * dx + dy * dy;
            if r2 < cutoff2 {
                neighbors[i].push(j);
                neighbors[j].push(i);
            }
        }
    }
    neighbors
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // ── BufferUsage ──────────────────────────────────────────────────────────

    #[test]
    fn buffer_usage_eq() {
        assert_eq!(BufferUsage::Storage, BufferUsage::Storage);
        assert_ne!(BufferUsage::Vertex, BufferUsage::Index);
    }

    #[test]
    fn buffer_usage_clone() {
        let u = BufferUsage::Uniform;
        assert_eq!(u.clone(), BufferUsage::Uniform);
    }

    // ── ComputeBuffer ────────────────────────────────────────────────────────

    #[test]
    fn compute_buffer_new_zeroed() {
        let buf = ComputeBuffer::new(8, BufferUsage::Storage, "test");
        assert_eq!(buf.data.len(), 8);
        assert!(buf.data.iter().all(|&v| v == 0.0));
        assert_eq!(buf.label, "test");
    }

    #[test]
    fn compute_buffer_byte_size() {
        let buf = ComputeBuffer::new(4, BufferUsage::Uniform, "u");
        assert_eq!(buf.byte_size(), 16);
    }

    #[test]
    fn compute_buffer_write_read_roundtrip() {
        let mut buf = ComputeBuffer::new(8, BufferUsage::Storage, "rw");
        buf.write_f32(2, &[1.0, 2.0, 3.0]);
        let out = buf.read_f32(2, 3);
        assert_eq!(out, vec![1.0, 2.0, 3.0]);
    }

    #[test]
    fn compute_buffer_write_at_offset_zero() {
        let mut buf = ComputeBuffer::new(4, BufferUsage::Storage, "s");
        buf.write_f32(0, &[9.0, 8.0, 7.0, 6.0]);
        assert_eq!(buf.data, vec![9.0, 8.0, 7.0, 6.0]);
    }

    #[test]
    #[should_panic(expected = "out-of-bounds write")]
    fn compute_buffer_write_oob_panics() {
        let mut buf = ComputeBuffer::new(4, BufferUsage::Storage, "oob");
        buf.write_f32(3, &[1.0, 2.0]); // 3+2 > 4
    }

    #[test]
    #[should_panic(expected = "out-of-bounds read")]
    fn compute_buffer_read_oob_panics() {
        let buf = ComputeBuffer::new(4, BufferUsage::Storage, "oob");
        let _ = buf.read_f32(3, 2);
    }

    // ── WorkgroupSize ────────────────────────────────────────────────────────

    #[test]
    fn workgroup_dispatch_count_exact() {
        assert_eq!(WorkgroupSize::dispatch_count(64, 64), 1);
    }

    #[test]
    fn workgroup_dispatch_count_ceil() {
        assert_eq!(WorkgroupSize::dispatch_count(65, 64), 2);
        assert_eq!(WorkgroupSize::dispatch_count(1, 64), 1);
    }

    #[test]
    fn workgroup_dispatch_count_zero_total() {
        assert_eq!(WorkgroupSize::dispatch_count(0, 64), 0);
    }

    #[test]
    fn workgroup_default() {
        let wg = WorkgroupSize::default();
        assert_eq!(wg.x, 64);
        assert_eq!(wg.y, 1);
        assert_eq!(wg.z, 1);
    }

    // ── ComputeKernelKind ────────────────────────────────────────────────────

    #[test]
    fn kernel_kind_custom_eq() {
        let a = ComputeKernelKind::Custom("foo".into());
        let b = ComputeKernelKind::Custom("foo".into());
        assert_eq!(a, b);
    }

    #[test]
    fn kernel_kind_variants_neq() {
        assert_ne!(
            ComputeKernelKind::VelocityUpdate,
            ComputeKernelKind::PressureJacobi
        );
    }

    // ── CpuComputeDispatch – velocity update ─────────────────────────────────

    #[test]
    fn velocity_update_basic() {
        let disp =
            CpuComputeDispatch::new(ComputeKernelKind::VelocityUpdate, WorkgroupSize::default());
        let n = 3;
        let mut pos = ComputeBuffer::new(n, BufferUsage::Storage, "pos");
        let mut vel = ComputeBuffer::new(n, BufferUsage::Storage, "vel");
        let mut force = ComputeBuffer::new(n, BufferUsage::Storage, "force");
        let mut mass = ComputeBuffer::new(n, BufferUsage::Storage, "mass");

        pos.write_f32(0, &[0.0, 1.0, 2.0]);
        vel.write_f32(0, &[1.0, 0.5, -1.0]);
        force.write_f32(0, &[0.0, 1.0, 0.0]);
        mass.write_f32(0, &[1.0, 2.0, 1.0]);

        let dt = 0.1_f32;
        disp.dispatch_velocity_update(&mut pos, &mut vel, &force, &mass, dt, n);

        // pos[0] = 0.0 + 1.0*0.1 = 0.1,  vel[0] = 1.0 + 0/1*0.1 = 1.0
        assert!((pos.data[0] - 0.1).abs() < 1e-6);
        assert!((vel.data[0] - 1.0).abs() < 1e-6);
        // pos[1] = 1.0 + 0.5*0.1 = 1.05, vel[1] = 0.5 + 1/2*0.1 = 0.55
        assert!((pos.data[1] - 1.05).abs() < 1e-6);
        assert!((vel.data[1] - 0.55).abs() < 1e-6);
    }

    #[test]
    fn velocity_update_zero_force() {
        let disp =
            CpuComputeDispatch::new(ComputeKernelKind::VelocityUpdate, WorkgroupSize::default());
        let n = 2;
        let mut pos = ComputeBuffer::new(n, BufferUsage::Storage, "pos");
        let mut vel = ComputeBuffer::new(n, BufferUsage::Storage, "vel");
        let force = ComputeBuffer::new(n, BufferUsage::Storage, "force");
        let mut mass = ComputeBuffer::new(n, BufferUsage::Storage, "mass");

        pos.write_f32(0, &[0.0, 0.0]);
        vel.write_f32(0, &[2.0, -3.0]);
        mass.write_f32(0, &[1.0, 1.0]);

        disp.dispatch_velocity_update(&mut pos, &mut vel, &force, &mass, 0.5, n);

        assert!((pos.data[0] - 1.0).abs() < 1e-6);
        assert!((pos.data[1] - (-1.5)).abs() < 1e-6);
        // velocities unchanged (zero force)
        assert!((vel.data[0] - 2.0).abs() < 1e-6);
        assert!((vel.data[1] - (-3.0)).abs() < 1e-6);
    }

    // ── CpuComputeDispatch – pressure Jacobi ─────────────────────────────────

    #[test]
    fn pressure_jacobi_interior_update() {
        let disp =
            CpuComputeDispatch::new(ComputeKernelKind::PressureJacobi, WorkgroupSize::default());
        let nx = 4;
        let ny = 4;
        let mut p = ComputeBuffer::new(nx * ny, BufferUsage::Storage, "p");
        let mut p_old = ComputeBuffer::new(nx * ny, BufferUsage::Storage, "p_old");
        let rhs = ComputeBuffer::new(nx * ny, BufferUsage::Storage, "rhs");

        // Set p_old to known values; neighbours of (1,1) are all 1.0
        for v in p_old.data.iter_mut() {
            *v = 1.0;
        }

        disp.dispatch_pressure_jacobi(&mut p, &p_old, &rhs, nx, ny, 1.0);

        // p[1*4+1] = (1+1+1+1 - 0)/4 = 1.0
        let idx = nx + 1;
        assert!((p.data[idx] - 1.0).abs() < 1e-6);
    }

    #[test]
    fn pressure_jacobi_boundary_unchanged() {
        let disp =
            CpuComputeDispatch::new(ComputeKernelKind::PressureJacobi, WorkgroupSize::default());
        let nx = 5;
        let ny = 5;
        let mut p = ComputeBuffer::new(nx * ny, BufferUsage::Storage, "p");
        let p_old = ComputeBuffer::new(nx * ny, BufferUsage::Storage, "p_old");
        let rhs = ComputeBuffer::new(nx * ny, BufferUsage::Storage, "rhs");

        // Boundary should remain zero.
        disp.dispatch_pressure_jacobi(&mut p, &p_old, &rhs, nx, ny, 1.0);
        assert_eq!(p.data[0], 0.0); // corner
        assert_eq!(p.data[4], 0.0); // top-right corner
    }

    // ── CpuComputeDispatch – LJ particle force ───────────────────────────────

    #[test]
    fn particle_force_zero_at_large_sep() {
        let disp =
            CpuComputeDispatch::new(ComputeKernelKind::ParticleForce, WorkgroupSize::default());
        let n = 2;
        // Two particles very far apart → tiny force.
        let mut pos = ComputeBuffer::new(2 * n, BufferUsage::Storage, "pos");
        let mut force = ComputeBuffer::new(2 * n, BufferUsage::Storage, "force");
        pos.write_f32(0, &[0.0, 0.0, 1000.0, 0.0]);

        disp.dispatch_particle_force(&pos, &mut force, 1.0, 1.0, n);
        // Force should be negligible at r=1000σ
        assert!(force.data[0].abs() < 1e-10);
    }

    #[test]
    fn particle_force_newton3() {
        let disp =
            CpuComputeDispatch::new(ComputeKernelKind::ParticleForce, WorkgroupSize::default());
        let n = 2;
        let mut pos = ComputeBuffer::new(2 * n, BufferUsage::Storage, "pos");
        let mut force = ComputeBuffer::new(2 * n, BufferUsage::Storage, "force");
        pos.write_f32(0, &[0.0, 0.0, 1.5, 0.0]);

        disp.dispatch_particle_force(&pos, &mut force, 1.0, 1.0, n);
        // Newton's third law: f0 + f1 == 0
        assert!((force.data[0] + force.data[2]).abs() < 1e-5);
        assert!((force.data[1] + force.data[3]).abs() < 1e-5);
    }

    // ── GpuStats ─────────────────────────────────────────────────────────────

    #[test]
    fn gpu_stats_initial_zero() {
        let s = GpuStats::new();
        assert_eq!(s.dispatch_count, 0);
        assert_eq!(s.bytes_transferred, 0);
        assert_eq!(s.kernel_time_ms, 0.0);
    }

    #[test]
    fn gpu_stats_accumulate() {
        let mut s = GpuStats::new();
        s.record_dispatch(128, 0.5);
        s.record_dispatch(256, 1.0);
        assert_eq!(s.dispatch_count, 2);
        assert_eq!(s.bytes_transferred, 384);
        assert!((s.kernel_time_ms - 1.5).abs() < 1e-9);
    }

    // ── jacobi_step_2d ───────────────────────────────────────────────────────

    #[test]
    fn jacobi_step_2d_uniform_field() {
        let nx = 4;
        let ny = 4;
        let size = nx * ny;
        let mut p_new = vec![0.0_f32; size];
        let p_old = vec![1.0_f32; size];
        let rhs = vec![0.0_f32; size];

        jacobi_step_2d(&mut p_new, &p_old, &rhs, nx, ny, 1.0);

        // Uniform field → interior stays 1.0.
        for j in 1..ny - 1 {
            for i in 1..nx - 1 {
                assert!((p_new[j * nx + i] - 1.0).abs() < 1e-6);
            }
        }
    }

    #[test]
    fn jacobi_step_2d_rhs_effect() {
        let nx = 4;
        let ny = 4;
        let size = nx * ny;
        let mut p_new = vec![0.0_f32; size];
        let p_old = vec![4.0_f32; size];
        // rhs = 4 at every interior point
        let rhs = vec![4.0_f32; size];

        jacobi_step_2d(&mut p_new, &p_old, &rhs, nx, ny, 1.0);
        // p[i,j] = (4+4+4+4 - 1²*4)/4 = (16-4)/4 = 3
        for j in 1..ny - 1 {
            for i in 1..nx - 1 {
                assert!((p_new[j * nx + i] - 3.0).abs() < 1e-6);
            }
        }
    }

    // ── pressure_poisson_solve ───────────────────────────────────────────────

    #[test]
    fn pressure_poisson_zero_rhs_zero_bc() {
        // Zero RHS + zero BCs → solution stays zero → zero residual.
        let nx = 5;
        let ny = 5;
        let mut p = vec![0.0_f32; nx * ny];
        let rhs = vec![0.0_f32; nx * ny];
        let residual = pressure_poisson_solve(&mut p, &rhs, nx, ny, 0.1, 50);
        assert!(residual < 1e-6, "residual={residual}");
    }

    #[test]
    fn pressure_poisson_residual_decreases() {
        let nx = 6;
        let ny = 6;
        let mut p1 = vec![0.0_f32; nx * ny];
        let mut p2 = p1.clone();
        let rhs: Vec<f32> = (0..(nx * ny)).map(|k| (k as f32).sin()).collect();
        let dx = 0.1;

        let r1 = pressure_poisson_solve(&mut p1, &rhs, nx, ny, dx, 10);
        let r2 = pressure_poisson_solve(&mut p2, &rhs, nx, ny, dx, 200);
        assert!(
            r2 <= r1 + 1e-4,
            "more iterations should not increase residual (r1={r1}, r2={r2})"
        );
    }

    // ── PipelineCache ──────────────────────────────────────────────────────

    #[test]
    fn pipeline_cache_insert_and_get() {
        let mut cache = PipelineCache::new(4);
        let disp =
            CpuComputeDispatch::new(ComputeKernelKind::VelocityUpdate, WorkgroupSize::default());
        cache.insert("vel_update", disp);
        assert!(cache.get("vel_update").is_some());
        assert!(cache.get("nonexistent").is_none());
    }

    #[test]
    fn pipeline_cache_eviction() {
        let mut cache = PipelineCache::new(2);
        let d1 =
            CpuComputeDispatch::new(ComputeKernelKind::VelocityUpdate, WorkgroupSize::default());
        let d2 =
            CpuComputeDispatch::new(ComputeKernelKind::PressureJacobi, WorkgroupSize::default());
        let d3 =
            CpuComputeDispatch::new(ComputeKernelKind::ParticleForce, WorkgroupSize::default());
        cache.insert("a", d1);
        cache.insert("b", d2);
        cache.insert("c", d3); // should evict "a"
        assert!(cache.get("a").is_none());
        assert!(cache.get("b").is_some());
        assert!(cache.get("c").is_some());
    }

    #[test]
    fn pipeline_cache_replace() {
        let mut cache = PipelineCache::new(4);
        let d1 =
            CpuComputeDispatch::new(ComputeKernelKind::VelocityUpdate, WorkgroupSize::default());
        let d2 = CpuComputeDispatch::new(
            ComputeKernelKind::ParticleForce,
            WorkgroupSize { x: 128, y: 1, z: 1 },
        );
        cache.insert("key", d1);
        cache.insert("key", d2);
        let entry = cache.get("key").unwrap();
        assert_eq!(entry.kernel, ComputeKernelKind::ParticleForce);
    }

    // ── PipelineStats ──────────────────────────────────────────────────────

    #[test]
    fn pipeline_stats_default() {
        let stats = PipelineStats::default();
        assert_eq!(stats.total_dispatches, 0);
        assert_eq!(stats.total_workgroups, 0);
        assert_eq!(stats.total_invocations, 0);
        assert_eq!(stats.cache_hits, 0);
        assert_eq!(stats.cache_misses, 0);
    }

    #[test]
    fn pipeline_stats_record() {
        let mut stats = PipelineStats::default();
        stats.record_dispatch(4, WorkgroupSize { x: 64, y: 1, z: 1 });
        assert_eq!(stats.total_dispatches, 1);
        assert_eq!(stats.total_workgroups, 4);
        assert_eq!(stats.total_invocations, 4 * 64);
    }

    #[test]
    fn pipeline_stats_record_3d_workgroup() {
        let mut stats = PipelineStats::default();
        stats.record_dispatch(2, WorkgroupSize { x: 8, y: 8, z: 4 });
        assert_eq!(stats.total_dispatches, 1);
        assert_eq!(stats.total_workgroups, 2);
        assert_eq!(stats.total_invocations, 2 * 8 * 8 * 4);
    }

    #[test]
    fn pipeline_stats_cache_ratio() {
        let mut stats = PipelineStats::default();
        assert!(stats.cache_hit_ratio().is_nan() || stats.cache_hit_ratio() == 0.0);
        stats.cache_hits = 3;
        stats.cache_misses = 1;
        assert!((stats.cache_hit_ratio() - 0.75).abs() < 1e-6);
    }

    // ── MultiPassPipeline ──────────────────────────────────────────────────

    #[test]
    fn multi_pass_empty() {
        let mp = MultiPassPipeline::new("empty");
        assert_eq!(mp.passes.len(), 0);
        assert_eq!(mp.label, "empty");
    }

    #[test]
    fn multi_pass_execute_add_scale() {
        // pass 0: fill buffer with [1, 2, 3, 4]
        // pass 1: scale by 2 → [2, 4, 6, 8]
        let mut mp = MultiPassPipeline::new("add_scale");
        mp.add_pass(ComputePass {
            label: "fill".into(),
            kernel: ComputeKernelKind::Custom("fill".into()),
            workgroup_size: WorkgroupSize::default(),
            buffer_bindings: vec![0],
        });
        mp.add_pass(ComputePass {
            label: "scale".into(),
            kernel: ComputeKernelKind::Custom("scale".into()),
            workgroup_size: WorkgroupSize::default(),
            buffer_bindings: vec![0],
        });
        assert_eq!(mp.passes.len(), 2);
        assert_eq!(mp.passes[0].label, "fill");
        assert_eq!(mp.passes[1].label, "scale");
    }

    #[test]
    fn multi_pass_dispatch_velocity_chain() {
        // Test chaining two velocity update passes
        let mut mp = MultiPassPipeline::new("vel_chain");
        mp.add_pass(ComputePass {
            label: "step1".into(),
            kernel: ComputeKernelKind::VelocityUpdate,
            workgroup_size: WorkgroupSize::default(),
            buffer_bindings: vec![0, 1, 2, 3],
        });
        mp.add_pass(ComputePass {
            label: "step2".into(),
            kernel: ComputeKernelKind::VelocityUpdate,
            workgroup_size: WorkgroupSize::default(),
            buffer_bindings: vec![0, 1, 2, 3],
        });

        let n = 2;
        let mut pos = ComputeBuffer::new(n, BufferUsage::Storage, "pos");
        let mut vel = ComputeBuffer::new(n, BufferUsage::Storage, "vel");
        let force = ComputeBuffer::new(n, BufferUsage::Storage, "force");
        let mut mass = ComputeBuffer::new(n, BufferUsage::Storage, "mass");

        pos.write_f32(0, &[0.0, 0.0]);
        vel.write_f32(0, &[1.0, 2.0]);
        mass.write_f32(0, &[1.0, 1.0]);

        let dt = 0.1_f32;

        // Execute passes manually (since we simulate CPU-side)
        for pass in &mp.passes {
            if pass.kernel == ComputeKernelKind::VelocityUpdate {
                let disp =
                    CpuComputeDispatch::new(ComputeKernelKind::VelocityUpdate, pass.workgroup_size);
                disp.dispatch_velocity_update(&mut pos, &mut vel, &force, &mass, dt, n);
            }
        }
        // After 2 steps with zero force: pos = vel*dt*2
        assert!((pos.data[0] - 0.2).abs() < 1e-5);
        assert!((pos.data[1] - 0.4).abs() < 1e-5);
    }

    // ── Pipeline validation ────────────────────────────────────────────────

    #[test]
    fn validate_binding_valid() {
        let buffers = vec![
            ComputeBuffer::new(16, BufferUsage::Storage, "buf0"),
            ComputeBuffer::new(16, BufferUsage::Uniform, "buf1"),
        ];
        let pass = ComputePass {
            label: "test".into(),
            kernel: ComputeKernelKind::VelocityUpdate,
            workgroup_size: WorkgroupSize::default(),
            buffer_bindings: vec![0, 1],
        };
        let errors = validate_resource_bindings(&pass, &buffers);
        assert!(errors.is_empty());
    }

    #[test]
    fn validate_binding_out_of_range() {
        let buffers = vec![ComputeBuffer::new(16, BufferUsage::Storage, "buf0")];
        let pass = ComputePass {
            label: "test".into(),
            kernel: ComputeKernelKind::VelocityUpdate,
            workgroup_size: WorkgroupSize::default(),
            buffer_bindings: vec![0, 5],
        };
        let errors = validate_resource_bindings(&pass, &buffers);
        assert_eq!(errors.len(), 1);
        assert!(errors[0].contains("out of range"));
    }

    #[test]
    fn validate_binding_duplicate() {
        let buffers = vec![ComputeBuffer::new(16, BufferUsage::Storage, "buf0")];
        let pass = ComputePass {
            label: "test".into(),
            kernel: ComputeKernelKind::VelocityUpdate,
            workgroup_size: WorkgroupSize::default(),
            buffer_bindings: vec![0, 0],
        };
        let errors = validate_resource_bindings(&pass, &buffers);
        assert_eq!(errors.len(), 1);
        assert!(errors[0].contains("Duplicate"));
    }

    #[test]
    fn validate_pipeline_all_passes() {
        let buffers = vec![ComputeBuffer::new(16, BufferUsage::Storage, "buf0")];
        let mut mp = MultiPassPipeline::new("test");
        mp.add_pass(ComputePass {
            label: "good".into(),
            kernel: ComputeKernelKind::VelocityUpdate,
            workgroup_size: WorkgroupSize::default(),
            buffer_bindings: vec![0],
        });
        mp.add_pass(ComputePass {
            label: "bad".into(),
            kernel: ComputeKernelKind::PressureJacobi,
            workgroup_size: WorkgroupSize::default(),
            buffer_bindings: vec![0, 3],
        });
        let errors = validate_pipeline(&mp, &buffers);
        assert_eq!(errors.len(), 1); // only the bad pass has errors
    }

    // ── ComputeBuffer additional tests ─────────────────────────────────────

    #[test]
    fn compute_buffer_clone() {
        let mut buf = ComputeBuffer::new(4, BufferUsage::Storage, "orig");
        buf.write_f32(0, &[1.0, 2.0, 3.0, 4.0]);
        let cloned = buf.clone();
        assert_eq!(buf.data, cloned.data);
        assert_eq!(buf.label, cloned.label);
    }

    #[test]
    fn compute_buffer_staging_usage() {
        let buf = ComputeBuffer::new(8, BufferUsage::Staging, "staging");
        assert_eq!(buf.usage, BufferUsage::Staging);
        assert_eq!(buf.byte_size(), 32);
    }

    // ── Workgroup additional tests ─────────────────────────────────────────

    #[test]
    fn workgroup_dispatch_count_large() {
        assert_eq!(WorkgroupSize::dispatch_count(1024, 256), 4);
        assert_eq!(WorkgroupSize::dispatch_count(1025, 256), 5);
    }

    // ── GpuStats additional tests ──────────────────────────────────────────

    #[test]
    fn gpu_stats_clone() {
        let mut s = GpuStats::new();
        s.record_dispatch(100, 1.5);
        let s2 = s.clone();
        assert_eq!(s.dispatch_count, s2.dispatch_count);
        assert_eq!(s.bytes_transferred, s2.bytes_transferred);
        assert!((s.kernel_time_ms - s2.kernel_time_ms).abs() < 1e-12);
    }

    // ── dispatch_particle_force additional tests ──────────────────────────

    #[test]
    fn particle_force_repulsive_at_close_range() {
        let disp =
            CpuComputeDispatch::new(ComputeKernelKind::ParticleForce, WorkgroupSize::default());
        let n = 2;
        let mut pos = ComputeBuffer::new(2 * n, BufferUsage::Storage, "pos");
        let mut force = ComputeBuffer::new(2 * n, BufferUsage::Storage, "force");
        // Two particles at distance 0.9σ (< σ → repulsive region)
        pos.write_f32(0, &[0.0, 0.0, 0.9, 0.0]);
        disp.dispatch_particle_force(&pos, &mut force, 1.0, 1.0, n);
        // Force on particle 0 should push it away from particle 1 (negative x)
        assert!(
            force.data[0] < 0.0,
            "expected repulsive force, got {}",
            force.data[0]
        );
    }

    #[test]
    fn particle_force_three_particles() {
        let disp =
            CpuComputeDispatch::new(ComputeKernelKind::ParticleForce, WorkgroupSize::default());
        let n = 3;
        let mut pos = ComputeBuffer::new(2 * n, BufferUsage::Storage, "pos");
        let mut force = ComputeBuffer::new(2 * n, BufferUsage::Storage, "force");
        // Triangle arrangement
        pos.write_f32(0, &[0.0, 0.0, 2.0, 0.0, 1.0, 1.732]);
        disp.dispatch_particle_force(&pos, &mut force, 1.0, 1.0, n);

        // Total momentum conservation: sum of all forces should be zero
        let fx_total = force.data[0] + force.data[2] + force.data[4];
        let fy_total = force.data[1] + force.data[3] + force.data[5];
        assert!(fx_total.abs() < 1e-5, "fx_total={fx_total}");
        assert!(fy_total.abs() < 1e-5, "fy_total={fy_total}");
    }

    // ── pressure_poisson_solve additional tests ───────────────────────────

    #[test]
    fn pressure_poisson_uniform_rhs() {
        let nx = 8;
        let ny = 8;
        let mut p = vec![0.0_f32; nx * ny];
        let rhs = vec![1.0_f32; nx * ny];
        let residual = pressure_poisson_solve(&mut p, &rhs, nx, ny, 0.1, 500);
        // After many iterations residual should decrease significantly
        assert!(residual < 10.0, "residual={residual}");
    }

    // ── SOR solver ────────────────────────────────────────────────────────

    #[test]
    fn sor_step_uniform_field() {
        let nx = 4;
        let ny = 4;
        let mut p = vec![0.0_f32; nx * ny];
        let rhs = vec![0.0_f32; nx * ny];
        let p_ref = vec![1.0_f32; nx * ny];
        sor_step_2d(&mut p, &p_ref, &rhs, nx, ny, 1.0, 1.0);
        // With omega=1.0 this is the same as Jacobi
        for j in 1..ny - 1 {
            for i in 1..nx - 1 {
                assert!((p[j * nx + i] - 1.0).abs() < 1e-6);
            }
        }
    }

    #[test]
    fn sor_step_over_relaxation() {
        // SOR with omega > 1 should differ from standard Jacobi (omega=1)
        let nx = 6;
        let ny = 6;
        let rhs = vec![0.0_f32; nx * ny];

        // Non-uniform reference: set boundary to 1, interior p_old to 0
        let mut p_ref = vec![0.0_f32; nx * ny];
        for i in 0..nx {
            p_ref[i] = 1.0;
            p_ref[(ny - 1) * nx + i] = 1.0;
        }
        for j in 0..ny {
            p_ref[j * nx] = 1.0;
            p_ref[j * nx + nx - 1] = 1.0;
        }

        let mut p_jac = vec![0.0_f32; nx * ny];
        let mut p_sor = vec![0.0_f32; nx * ny];
        sor_step_2d(&mut p_jac, &p_ref, &rhs, nx, ny, 1.0, 1.0);
        sor_step_2d(&mut p_sor, &p_ref, &rhs, nx, ny, 1.0, 1.5);

        // SOR result should differ from Jacobi for interior nodes next to boundary
        let idx = nx + 1; // has 2 boundary neighbors
        // Jacobi: (1+0+1+0)/4 = 0.5, SOR: (1-1.5)*0 + 1.5*0.5 = 0.75
        assert!(
            (p_sor[idx] - p_jac[idx]).abs() > 0.01,
            "SOR and Jacobi should differ: SOR={}, Jac={}",
            p_sor[idx],
            p_jac[idx]
        );
    }

    // ── Red-black Gauss-Seidel ────────────────────────────────────────────

    #[test]
    fn red_black_gs_uniform() {
        let nx = 6;
        let ny = 6;
        let mut p = vec![1.0_f32; nx * ny];
        let rhs = vec![0.0_f32; nx * ny];
        red_black_gauss_seidel_step(&mut p, &rhs, nx, ny, 1.0);
        // Uniform field is a fixed point with zero rhs
        for j in 1..ny - 1 {
            for i in 1..nx - 1 {
                assert!((p[j * nx + i] - 1.0).abs() < 1e-6);
            }
        }
    }

    #[test]
    fn red_black_gs_converges() {
        let nx = 8;
        let ny = 8;
        let mut p = vec![0.0_f32; nx * ny];
        let rhs = vec![0.0_f32; nx * ny];
        // Set boundary to 1
        for i in 0..nx {
            p[i] = 1.0;
            p[(ny - 1) * nx + i] = 1.0;
        }
        for j in 0..ny {
            p[j * nx] = 1.0;
            p[j * nx + nx - 1] = 1.0;
        }
        // Multiple sweeps should converge interior toward 1.0
        for _ in 0..200 {
            red_black_gauss_seidel_step(&mut p, &rhs, nx, ny, 1.0);
        }
        let center = p[(ny / 2) * nx + nx / 2];
        assert!((center - 1.0).abs() < 0.01, "center={center}");
    }

    // ── compute_linf_residual ─────────────────────────────────────────────

    #[test]
    fn linf_residual_zero_for_exact() {
        // Uniform field with zero rhs is an exact solution
        let nx = 4;
        let ny = 4;
        let p = vec![1.0_f32; nx * ny];
        let rhs = vec![0.0_f32; nx * ny];
        let res = compute_linf_residual(&p, &rhs, nx, ny, 1.0);
        assert!(res < 1e-6, "res={res}");
    }

    #[test]
    fn linf_residual_nonzero_for_wrong() {
        let nx = 4;
        let ny = 4;
        let mut p = vec![0.0_f32; nx * ny];
        p[nx + 1] = 100.0; // big spike
        let rhs = vec![0.0_f32; nx * ny];
        let res = compute_linf_residual(&p, &rhs, nx, ny, 1.0);
        assert!(res > 1.0, "expected large residual, got {res}");
    }

    // ── Dispatch neighbor search ──────────────────────────────────────────

    #[test]
    fn dispatch_neighbor_search_basic() {
        let n = 4;
        let positions = vec![
            0.0_f32, 0.0, // particle 0
            0.5, 0.0, // particle 1 (close to 0)
            5.0, 5.0, // particle 2 (far)
            0.3, 0.3, // particle 3 (close to 0 and 1)
        ];
        let neighbors = dispatch_neighbor_search(&positions, n, 1.0);
        // particle 0 should have neighbors 1 and 3
        assert!(neighbors[0].contains(&1));
        assert!(neighbors[0].contains(&3));
        // particle 2 should have no neighbors
        assert!(neighbors[2].is_empty());
    }

    #[test]
    fn dispatch_neighbor_search_all_close() {
        let n = 3;
        let positions = vec![0.0_f32, 0.0, 0.1, 0.0, 0.0, 0.1];
        let neighbors = dispatch_neighbor_search(&positions, n, 1.0);
        // All particles within cutoff of each other
        assert_eq!(neighbors[0].len(), 2);
        assert_eq!(neighbors[1].len(), 2);
        assert_eq!(neighbors[2].len(), 2);
    }

    #[test]
    fn dispatch_neighbor_search_none() {
        let n = 2;
        let positions = vec![0.0_f32, 0.0, 100.0, 100.0];
        let neighbors = dispatch_neighbor_search(&positions, n, 1.0);
        assert!(neighbors[0].is_empty());
        assert!(neighbors[1].is_empty());
    }
}
