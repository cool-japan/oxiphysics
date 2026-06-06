//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use std::cell::RefCell;
use std::collections::HashMap;

/// Records a sequence of kernel dispatches for batched execution.
///
/// Compute passes accumulate dispatch commands and execute them in order.
pub struct ComputePass {
    /// Recorded dispatch commands: (kernel_name, work_size).
    pub(super) commands: Vec<(String, usize)>,
}
impl ComputePass {
    /// Create a new empty compute pass.
    pub fn new() -> Self {
        Self {
            commands: Vec::new(),
        }
    }
    /// Record a dispatch command.
    pub fn dispatch(&mut self, kernel_name: &str, work_size: usize) {
        self.commands.push((kernel_name.to_string(), work_size));
    }
    /// Return the number of recorded commands.
    pub fn num_commands(&self) -> usize {
        self.commands.len()
    }
    /// Return the recorded commands.
    pub fn commands(&self) -> &[(String, usize)] {
        &self.commands
    }
    /// Clear all recorded commands.
    pub fn clear(&mut self) {
        self.commands.clear();
    }
    /// Total work items across all recorded dispatches.
    pub fn total_work_items(&self) -> usize {
        self.commands.iter().map(|(_, ws)| ws).sum()
    }
}
/// Describes how a buffer is used in a compute pass.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BufferUsage {
    /// Buffer is read-only (storage, read).
    ReadOnly,
    /// Buffer is write-only (storage, read_write but only written).
    WriteOnly,
    /// Buffer is read-write.
    ReadWrite,
    /// Buffer is a uniform (small, read-only parameters).
    Uniform,
}
/// A single GPU command entry.
#[derive(Debug, Clone)]
pub enum GpuCommand {
    /// Copy from one buffer to another.
    CopyBuffer {
        /// Source buffer identifier.
        src: BufferId,
        /// Destination buffer identifier.
        dst: BufferId,
        /// Number of bytes to copy.
        size: usize,
    },
    /// Dispatch a compute kernel.
    DispatchCompute {
        /// Name of the compute kernel.
        kernel_name: String,
        /// Workgroup counts for each dimension.
        workgroups: [u32; 3],
    },
    /// Insert a pipeline barrier.
    Barrier(PipelineBarrier),
    /// Set a push constant value.
    PushConstant {
        /// Push constant name.
        name: String,
        /// Push constant value.
        value: f64,
    },
}
/// Tracks buffer lifecycle (creation, writes, reads) for debugging.
pub struct ResourceLifecycle {
    pub(super) events: Vec<ResourceEvent>,
}
impl ResourceLifecycle {
    /// Create a new lifecycle tracker.
    pub fn new() -> Self {
        Self { events: Vec::new() }
    }
    /// Record a creation event.
    pub fn record_create(&mut self, id: BufferId, size: usize) {
        self.events.push(ResourceEvent::Created(id, size));
    }
    /// Record a write event.
    pub fn record_write(&mut self, id: BufferId) {
        self.events.push(ResourceEvent::Written(id));
    }
    /// Record a read event.
    pub fn record_read(&mut self, id: BufferId) {
        self.events.push(ResourceEvent::Read(id));
    }
    /// Record a destroy event.
    pub fn record_destroy(&mut self, id: BufferId) {
        self.events.push(ResourceEvent::Destroyed(id));
    }
    /// Return all events.
    pub fn events(&self) -> &[ResourceEvent] {
        &self.events
    }
    /// Return the number of events recorded.
    pub fn len(&self) -> usize {
        self.events.len()
    }
    /// Check if no events have been recorded.
    pub fn is_empty(&self) -> bool {
        self.events.is_empty()
    }
    /// Clear all events.
    pub fn clear(&mut self) {
        self.events.clear();
    }
    /// Count events of a specific type for a given buffer.
    pub fn count_writes(&self, id: BufferId) -> usize {
        self.events
            .iter()
            .filter(|e| matches!(e, ResourceEvent::Written(bid) if * bid == id))
            .count()
    }
    /// Count reads for a given buffer.
    pub fn count_reads(&self, id: BufferId) -> usize {
        self.events
            .iter()
            .filter(|e| matches!(e, ResourceEvent::Read(bid) if * bid == id))
            .count()
    }
}
/// Specifies the type of pipeline barrier needed between passes.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum PipelineBarrier {
    /// Ensure all writes to storage buffers are visible before reading.
    StorageReadAfterWrite,
    /// Ensure all writes to uniform buffers are visible.
    UniformReadAfterWrite,
    /// Full barrier (all types).
    Full,
    /// No barrier needed.
    None,
}
/// GPU occupancy model (simplified).
///
/// Models occupancy as the ratio of active warps to maximum concurrent warps.
#[derive(Debug, Clone)]
pub struct OccupancyModel {
    /// Total number of compute units (SMs / CUs).
    pub compute_units: u32,
    /// Maximum warps per compute unit.
    pub max_warps_per_cu: u32,
    /// Warp size (threads per warp, typically 32 for NVIDIA or 64 for AMD).
    pub warp_size: u32,
    /// Shared memory per compute unit (bytes).
    pub shared_mem_per_cu: u32,
    /// Registers per compute unit.
    pub registers_per_cu: u32,
}
impl OccupancyModel {
    /// Create a model resembling a mid-range discrete GPU.
    pub fn mid_range() -> Self {
        Self {
            compute_units: 32,
            max_warps_per_cu: 32,
            warp_size: 32,
            shared_mem_per_cu: 48 * 1024,
            registers_per_cu: 65536,
        }
    }
    /// Estimate the theoretical occupancy (0.0–1.0) for a kernel.
    ///
    /// Occupancy is limited by:
    /// 1. Workgroup size (must not exceed warp_size * max_warps_per_cu).
    /// 2. Shared memory usage.
    /// 3. Register usage.
    pub fn estimate_occupancy(
        &self,
        workgroup_size: u32,
        shared_mem_bytes: u32,
        registers_per_thread: u32,
    ) -> f64 {
        let warps_per_wg = workgroup_size.div_ceil(self.warp_size);
        let max_wg_by_warps = self.max_warps_per_cu / warps_per_wg.max(1);
        let max_wg_by_smem = self
            .shared_mem_per_cu
            .checked_div(shared_mem_bytes)
            .unwrap_or(u32::MAX);
        let regs_per_wg = registers_per_thread * workgroup_size;
        let max_wg_by_regs = self
            .registers_per_cu
            .checked_div(regs_per_wg)
            .unwrap_or(u32::MAX);
        let active_wg = max_wg_by_warps.min(max_wg_by_smem).min(max_wg_by_regs);
        let active_warps = (active_wg * warps_per_wg).min(self.max_warps_per_cu);
        (active_warps as f64 / self.max_warps_per_cu as f64).clamp(0.0, 1.0)
    }
    /// Total theoretical peak throughput in GFLOP/s (mock model).
    ///
    /// Assumes 2 FP32 ops per clock per SIMD unit.
    pub fn peak_gflops(&self, clock_mhz: f64) -> f64 {
        let simd_width = self.warp_size as f64;
        2.0 * simd_width * self.compute_units as f64 * clock_mhz * 1e6 / 1e9
    }
}
/// A recorded sequence of GPU commands (mock encoder).
///
/// Records commands for later submission; models wgpu-style recording.
pub struct GpuCommandEncoder {
    pub(super) label: String,
    pub(super) commands: Vec<GpuCommand>,
}
impl GpuCommandEncoder {
    /// Create a new command encoder with a debug label.
    pub fn new(label: impl Into<String>) -> Self {
        Self {
            label: label.into(),
            commands: Vec::new(),
        }
    }
    /// Record a buffer copy command.
    pub fn copy_buffer(&mut self, src: BufferId, dst: BufferId, size: usize) {
        self.commands
            .push(GpuCommand::CopyBuffer { src, dst, size });
    }
    /// Record a compute dispatch.
    pub fn dispatch_compute(&mut self, kernel_name: &str, workgroups: [u32; 3]) {
        self.commands.push(GpuCommand::DispatchCompute {
            kernel_name: kernel_name.to_string(),
            workgroups,
        });
    }
    /// Insert a pipeline barrier.
    pub fn insert_barrier(&mut self, barrier: PipelineBarrier) {
        self.commands.push(GpuCommand::Barrier(barrier));
    }
    /// Set a named push constant.
    pub fn push_constant(&mut self, name: &str, value: f64) {
        self.commands.push(GpuCommand::PushConstant {
            name: name.to_string(),
            value,
        });
    }
    /// Return the label of this encoder.
    pub fn label(&self) -> &str {
        &self.label
    }
    /// Number of recorded commands.
    pub fn command_count(&self) -> usize {
        self.commands.len()
    }
    /// Return the recorded commands.
    pub fn commands(&self) -> &[GpuCommand] {
        &self.commands
    }
    /// Reset the encoder (clear recorded commands).
    pub fn reset(&mut self) {
        self.commands.clear();
    }
    /// "Submit" the recorded commands: replay them on the dispatcher.
    ///
    /// For copy commands, data is transferred between buffers.
    /// Other commands are noted but not executed (they are mock-only).
    pub fn submit(&self, dispatcher: &mut ComputeDispatcher) -> Result<(), GpuError> {
        for cmd in &self.commands {
            if let GpuCommand::CopyBuffer { src, dst, .. } = cmd {
                dispatcher.copy_buffer(*src, *dst)?;
            }
        }
        Ok(())
    }
}
/// Manages GPU buffers and dispatches parallel map/reduce operations.
///
/// This is a CPU-side simulation of GPU compute dispatch that uses a thread
/// pool metaphor (sequential execution) for testing without a real GPU.
pub struct ComputeDispatcher {
    pub(super) buffers: HashMap<BufferId, GpuBuffer>,
    pub(super) next_id: u32,
}
impl ComputeDispatcher {
    /// Create a new dispatcher with no buffers.
    pub fn new() -> Self {
        Self {
            buffers: HashMap::new(),
            next_id: 0,
        }
    }
    /// Allocate a new buffer of `size` f64 elements, optionally pre-loaded
    /// with `initial_data`.  Returns the new buffer's [`BufferId`].
    pub fn create_buffer(&mut self, size: usize, initial_data: Option<&[f64]>) -> BufferId {
        let id = BufferId(self.next_id);
        self.next_id += 1;
        let buf = match initial_data {
            Some(data) => {
                let mut b = GpuBuffer::new(size);
                let copy_len = data.len().min(size);
                b.data[..copy_len].copy_from_slice(&data[..copy_len]);
                b
            }
            None => GpuBuffer::new(size),
        };
        self.buffers.insert(id, buf);
        id
    }
    /// Write `data` into the buffer identified by `id`.
    ///
    /// # Errors
    /// Returns [`GpuError::InvalidBuffer`] if `id` is not registered.
    pub fn write_buffer(&mut self, id: BufferId, data: &[f64]) -> Result<(), GpuError> {
        match self.buffers.get_mut(&id) {
            Some(buf) => {
                buf.data = data.to_vec();
                buf.size = data.len();
                Ok(())
            }
            None => Err(GpuError::InvalidBuffer(id)),
        }
    }
    /// Read the contents of the buffer identified by `id`.
    ///
    /// # Errors
    /// Returns [`GpuError::InvalidBuffer`] if `id` is not registered.
    pub fn read_buffer(&self, id: BufferId) -> Result<Vec<f64>, GpuError> {
        self.buffers
            .get(&id)
            .map(|b| b.data.clone())
            .ok_or(GpuError::InvalidBuffer(id))
    }
    /// Return the number of buffers currently managed.
    pub fn num_buffers(&self) -> usize {
        self.buffers.len()
    }
    /// Check if a buffer exists.
    pub fn has_buffer(&self, id: BufferId) -> bool {
        self.buffers.contains_key(&id)
    }
    /// Return the size of a buffer.
    pub fn buffer_size(&self, id: BufferId) -> Result<usize, GpuError> {
        self.buffers
            .get(&id)
            .map(|b| b.size)
            .ok_or(GpuError::InvalidBuffer(id))
    }
    /// Destroy (remove) a buffer.
    pub fn destroy_buffer(&mut self, id: BufferId) -> Result<(), GpuError> {
        self.buffers
            .remove(&id)
            .map(|_| ())
            .ok_or(GpuError::InvalidBuffer(id))
    }
    /// Copy data from one buffer to another.
    pub fn copy_buffer(&mut self, src: BufferId, dst: BufferId) -> Result<(), GpuError> {
        let src_data = self
            .buffers
            .get(&src)
            .ok_or(GpuError::InvalidBuffer(src))?
            .data
            .clone();
        let dst_buf = self
            .buffers
            .get_mut(&dst)
            .ok_or(GpuError::InvalidBuffer(dst))?;
        if src_data.len() != dst_buf.size {
            return Err(GpuError::SizeMismatch {
                expected: dst_buf.size,
                got: src_data.len(),
            });
        }
        dst_buf.data = src_data;
        Ok(())
    }
    /// Dispatch a parallel map: `out[i] = f(in[i])` for every element.
    ///
    /// # Errors
    /// Returns an error if either buffer is invalid or sizes differ.
    pub fn dispatch_map(
        &mut self,
        buf_in: BufferId,
        buf_out: BufferId,
        f: impl Fn(f64) -> f64,
    ) -> Result<(), GpuError> {
        let input = self
            .buffers
            .get(&buf_in)
            .ok_or(GpuError::InvalidBuffer(buf_in))?
            .data
            .clone();
        let out_buf = self
            .buffers
            .get_mut(&buf_out)
            .ok_or(GpuError::InvalidBuffer(buf_out))?;
        if input.len() != out_buf.size {
            return Err(GpuError::SizeMismatch {
                expected: out_buf.size,
                got: input.len(),
            });
        }
        out_buf.data = input.iter().map(|&x| f(x)).collect();
        Ok(())
    }
    /// Dispatch a parallel map with index: `out[i] = f(i, in[i])`.
    pub fn dispatch_map_indexed(
        &mut self,
        buf_in: BufferId,
        buf_out: BufferId,
        f: impl Fn(usize, f64) -> f64,
    ) -> Result<(), GpuError> {
        let input = self
            .buffers
            .get(&buf_in)
            .ok_or(GpuError::InvalidBuffer(buf_in))?
            .data
            .clone();
        let out_buf = self
            .buffers
            .get_mut(&buf_out)
            .ok_or(GpuError::InvalidBuffer(buf_out))?;
        if input.len() != out_buf.size {
            return Err(GpuError::SizeMismatch {
                expected: out_buf.size,
                got: input.len(),
            });
        }
        out_buf.data = input.iter().enumerate().map(|(i, &x)| f(i, x)).collect();
        Ok(())
    }
    /// Dispatch a zip-map: `out[i] = f(a[i], b[i])`.
    pub fn dispatch_zip_map(
        &mut self,
        buf_a: BufferId,
        buf_b: BufferId,
        buf_out: BufferId,
        f: impl Fn(f64, f64) -> f64,
    ) -> Result<(), GpuError> {
        let a_data = self
            .buffers
            .get(&buf_a)
            .ok_or(GpuError::InvalidBuffer(buf_a))?
            .data
            .clone();
        let b_data = self
            .buffers
            .get(&buf_b)
            .ok_or(GpuError::InvalidBuffer(buf_b))?
            .data
            .clone();
        if a_data.len() != b_data.len() {
            return Err(GpuError::SizeMismatch {
                expected: a_data.len(),
                got: b_data.len(),
            });
        }
        let out_buf = self
            .buffers
            .get_mut(&buf_out)
            .ok_or(GpuError::InvalidBuffer(buf_out))?;
        if a_data.len() != out_buf.size {
            return Err(GpuError::SizeMismatch {
                expected: out_buf.size,
                got: a_data.len(),
            });
        }
        out_buf.data = a_data
            .iter()
            .zip(b_data.iter())
            .map(|(&a, &b)| f(a, b))
            .collect();
        Ok(())
    }
    /// Dispatch a parallel reduce: folds all elements in `buf` using `f`.
    ///
    /// Mimics a GPU tree-reduction (sequential here for correctness).
    ///
    /// # Errors
    /// Returns [`GpuError::InvalidBuffer`] or [`GpuError::EmptyBuffer`].
    pub fn dispatch_reduce(
        &self,
        buf: BufferId,
        f: impl Fn(f64, f64) -> f64,
    ) -> Result<f64, GpuError> {
        let data = self.buffers.get(&buf).ok_or(GpuError::InvalidBuffer(buf))?;
        let mut iter = data.data.iter().copied();
        let first = iter.next().ok_or(GpuError::EmptyBuffer)?;
        Ok(iter.fold(first, f))
    }
    /// Dispatch a mock SPH density kernel.
    ///
    /// Computes a simplified SPH density estimate for each particle:
    ///
    /// ```text
    /// rho_i = sum_j m_j * W(|r_i - r_j|, h)
    /// ```
    ///
    /// where `W(r, h) = max(0, 1 - (r/h)^2)` (simplified poly-6 mock).
    ///
    /// Buffer layout (flat f64):
    /// * `pos_buf` — `[x0, y0, z0, x1, y1, z1, ...]`
    /// * `mass_buf` — `[m0, m1, ...]`
    /// * `out_density_buf` — written with `[rho0, rho1, ...]`
    ///
    /// # Errors
    /// Returns an error if any buffer id is invalid.
    pub fn dispatch_sph_density(
        &mut self,
        pos_buf: BufferId,
        mass_buf: BufferId,
        h: f64,
        out_density_buf: BufferId,
    ) -> Result<(), GpuError> {
        let positions = self
            .buffers
            .get(&pos_buf)
            .ok_or(GpuError::InvalidBuffer(pos_buf))?
            .data
            .clone();
        let masses = self
            .buffers
            .get(&mass_buf)
            .ok_or(GpuError::InvalidBuffer(mass_buf))?
            .data
            .clone();
        let n = positions.len() / 3;
        let h2 = h * h;
        let mut densities = vec![0.0f64; n];
        for i in 0..n {
            let xi = positions[i * 3];
            let yi = positions[i * 3 + 1];
            let zi = positions[i * 3 + 2];
            let mut rho = 0.0;
            for j in 0..n {
                let dx = xi - positions[j * 3];
                let dy = yi - positions[j * 3 + 1];
                let dz = zi - positions[j * 3 + 2];
                let r2 = dx * dx + dy * dy + dz * dz;
                if r2 < h2 {
                    let q = 1.0 - r2 / h2;
                    rho += masses[j] * q * q;
                }
            }
            densities[i] = rho;
        }
        let out_buf = self
            .buffers
            .get_mut(&out_density_buf)
            .ok_or(GpuError::InvalidBuffer(out_density_buf))?;
        out_buf.data = densities;
        out_buf.size = n;
        Ok(())
    }
    /// Dispatch a tree-based parallel reduction on a buffer.
    ///
    /// Simulates a GPU tree reduction: repeatedly halves the active range,
    /// summing adjacent elements until one value remains.
    ///
    /// Returns the reduced value (identity `0.0` for an empty buffer).
    pub fn dispatch_reduction_tree(&self, buf: BufferId) -> Result<f64, GpuError> {
        let data = self
            .buffers
            .get(&buf)
            .ok_or(GpuError::InvalidBuffer(buf))?
            .data
            .clone();
        if data.is_empty() {
            return Ok(0.0);
        }
        let mut work = data;
        let mut len = work.len();
        while len > 1 {
            let half = len / 2;
            for i in 0..half {
                work[i] = work[i * 2] + work[i * 2 + 1];
            }
            if len % 2 == 1 {
                work[half] = work[len - 1];
                len = half + 1;
            } else {
                len = half;
            }
        }
        Ok(work[0])
    }
    /// Dispatch an inclusive prefix scan (cumulative sum) on a buffer.
    ///
    /// Writes `out[i] = sum(in[0..=i])` into `out_buf`.
    ///
    /// Uses a sequential Hillis-Steele-style scan for correctness.
    pub fn dispatch_inclusive_scan(
        &mut self,
        buf_in: BufferId,
        buf_out: BufferId,
    ) -> Result<(), GpuError> {
        let data = self
            .buffers
            .get(&buf_in)
            .ok_or(GpuError::InvalidBuffer(buf_in))?
            .data
            .clone();
        let n = data.len();
        let mut result = data;
        for i in 1..n {
            result[i] += result[i - 1];
        }
        let out = self
            .buffers
            .get_mut(&buf_out)
            .ok_or(GpuError::InvalidBuffer(buf_out))?;
        out.data = result;
        out.size = n;
        Ok(())
    }
    /// Dispatch a 2-bit radix sort on a buffer of non-negative f64 values.
    ///
    /// Values are cast to u64 (bit-cast) and sorted by their bits using
    /// counting sort passes with 2-bit digits.  32 passes cover all 64 bits.
    /// For non-negative IEEE 754 doubles the bit pattern order matches numeric
    /// order.  Returns the sorted data as a new `Vec`f64` (input unchanged).
    pub fn dispatch_radix_sort(&self, buf: BufferId) -> Result<Vec<f64>, GpuError> {
        let data = self
            .buffers
            .get(&buf)
            .ok_or(GpuError::InvalidBuffer(buf))?
            .data
            .clone();
        let n = data.len();
        if n == 0 {
            return Ok(Vec::new());
        }
        let mut keys: Vec<u64> = data.iter().map(|&v| v.to_bits()).collect();
        for pass in 0..32usize {
            let shift = pass * 2;
            let mut counts = [0usize; 4];
            for &k in &keys {
                counts[((k >> shift) & 0x3) as usize] += 1;
            }
            let mut starts = [0usize; 4];
            for i in 1..4 {
                starts[i] = starts[i - 1] + counts[i - 1];
            }
            let mut out = vec![0u64; n];
            let mut pos = starts;
            for &k in &keys {
                let digit = ((k >> shift) & 0x3) as usize;
                out[pos[digit]] = k;
                pos[digit] += 1;
            }
            keys = out;
        }
        Ok(keys.iter().map(|&bits| f64::from_bits(bits)).collect())
    }
}
/// Opaque handle to a GPU/CPU buffer (usize-indexed, used by ComputeBackend).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct BufferHandle(pub usize);
/// Specification for a GPU compute kernel dispatch.
#[derive(Debug, Clone)]
pub struct KernelSpec {
    /// Human-readable kernel name.
    pub name: String,
    /// Number of threads per workgroup `\[X, Y, Z\]`.
    pub workgroup_size: [u32; 3],
    /// Ordered list of buffer bindings for this kernel.
    pub buffer_bindings: Vec<BufferId>,
}
impl KernelSpec {
    /// Create a new kernel spec with a 1-D workgroup.
    pub fn new(name: impl Into<String>, workgroup_x: u32, buffer_bindings: Vec<BufferId>) -> Self {
        Self {
            name: name.into(),
            workgroup_size: [workgroup_x, 1, 1],
            buffer_bindings,
        }
    }
    /// Create a kernel spec with a 3-D workgroup size.
    pub fn with_workgroup_3d(
        name: impl Into<String>,
        workgroup_size: [u32; 3],
        buffer_bindings: Vec<BufferId>,
    ) -> Self {
        Self {
            name: name.into(),
            workgroup_size,
            buffer_bindings,
        }
    }
    /// Compute the number of workgroups needed for `total_items` in the X dimension.
    pub fn num_workgroups_x(&self, total_items: u32) -> u32 {
        total_items.div_ceil(self.workgroup_size[0])
    }
    /// Total threads per workgroup.
    pub fn threads_per_workgroup(&self) -> u32 {
        self.workgroup_size[0] * self.workgroup_size[1] * self.workgroup_size[2]
    }
}
/// A CPU-resident buffer that mimics a GPU storage buffer.
#[derive(Debug, Clone)]
pub struct GpuBuffer {
    /// Buffer contents as f64.
    pub data: Vec<f64>,
    /// Declared capacity of the buffer (may differ from `data.len()` if
    /// the buffer was created with a fixed size but partially written).
    pub size: usize,
}
impl GpuBuffer {
    /// Create a zero-filled buffer of `size` elements.
    pub fn new(size: usize) -> Self {
        Self {
            data: vec![0.0; size],
            size,
        }
    }
    /// Create a buffer pre-loaded with `initial_data`.
    pub fn from_data(initial_data: Vec<f64>) -> Self {
        let size = initial_data.len();
        Self {
            data: initial_data,
            size,
        }
    }
    /// Fill the buffer with a constant value.
    pub fn fill(&mut self, value: f64) {
        for v in &mut self.data {
            *v = value;
        }
    }
    /// Clear the buffer (set all elements to 0).
    pub fn clear(&mut self) {
        self.fill(0.0);
    }
    /// Get a slice of the buffer data.
    pub fn as_slice(&self) -> &[f64] {
        &self.data
    }
    /// Get a mutable slice of the buffer data.
    pub fn as_mut_slice(&mut self) -> &mut [f64] {
        &mut self.data
    }
    /// Number of bytes the buffer would occupy on GPU (f64 = 8 bytes each).
    pub fn byte_size(&self) -> usize {
        self.size * std::mem::size_of::<f64>()
    }
}
/// Errors that can occur during GPU (or CPU-fallback) dispatch.
#[derive(Debug, Clone, PartialEq)]
pub enum GpuError {
    /// The specified buffer id is not registered.
    InvalidBuffer(BufferId),
    /// Input and output buffers must have the same length.
    SizeMismatch {
        /// Expected buffer size.
        expected: usize,
        /// Actual buffer size received.
        got: usize,
    },
    /// The reduction was attempted on an empty buffer.
    EmptyBuffer,
    /// A kernel or operation was not found.
    NotFound(String),
}
/// CPU fallback compute backend.
///
/// Stores buffers as `Vec`f64` in memory and dispatches kernels on the CPU.
pub struct CpuBackend {
    pub(super) buffers: RefCell<Vec<Vec<f64>>>,
}
impl CpuBackend {
    /// Create a new CPU backend.
    pub fn new() -> Self {
        Self {
            buffers: RefCell::new(Vec::new()),
        }
    }
    /// Return the number of buffers currently allocated.
    pub fn num_buffers(&self) -> usize {
        self.buffers.borrow().len()
    }
    /// Return the total number of f64 elements across all buffers.
    pub fn total_elements(&self) -> usize {
        self.buffers.borrow().iter().map(|b| b.len()).sum()
    }
}
/// A single resource event for lifecycle tracking.
#[derive(Debug, Clone)]
pub enum ResourceEvent {
    /// Buffer was created.
    Created(BufferId, usize),
    /// Buffer was written to.
    Written(BufferId),
    /// Buffer was read from.
    Read(BufferId),
    /// Buffer was destroyed.
    Destroyed(BufferId),
}
/// A record of divergent branches observed in a kernel.
#[derive(Debug, Clone, Default)]
pub struct WarpDivergenceRecord {
    /// Number of branch instructions encountered.
    pub total_branches: u64,
    /// Number of branches where threads diverged (not all took same path).
    pub divergent_branches: u64,
}
impl WarpDivergenceRecord {
    /// Compute the divergence rate (0.0 = no divergence, 1.0 = fully divergent).
    pub fn divergence_rate(&self) -> f64 {
        if self.total_branches == 0 {
            0.0
        } else {
            self.divergent_branches as f64 / self.total_branches as f64
        }
    }
    /// Estimated performance penalty from divergence (relative slowdown factor).
    ///
    /// A simple model: penalty = 1 + divergence_rate * (warp_size - 1) / warp_size.
    pub fn performance_penalty(&self, warp_size: u32) -> f64 {
        let rate = self.divergence_rate();
        1.0 + rate * (warp_size as f64 - 1.0) / warp_size as f64
    }
}
/// A mock GPU timeline semaphore for synchronising multi-pass GPU work.
///
/// On real GPU APIs (Vulkan, D3D12), timeline semaphores allow the CPU to
/// wait for a specific GPU progress point.  This mock records signal and
/// wait operations for testing.
pub struct TimelineSemaphore {
    /// Current value of the semaphore counter.
    pub value: u64,
    /// History of signalled values.
    pub(super) signal_history: Vec<u64>,
    /// History of wait requests.
    pub(super) wait_history: Vec<u64>,
}
impl TimelineSemaphore {
    /// Create a new semaphore starting at value 0.
    pub fn new() -> Self {
        Self {
            value: 0,
            signal_history: Vec::new(),
            wait_history: Vec::new(),
        }
    }
    /// Signal the semaphore to `new_value`.  Value must be monotonically increasing.
    pub fn signal(&mut self, new_value: u64) {
        assert!(
            new_value > self.value,
            "semaphore values must increase monotonically"
        );
        self.value = new_value;
        self.signal_history.push(new_value);
    }
    /// Record a wait-for request.  In a mock, this checks `wait_value <= current`.
    ///
    /// Returns `true` if the semaphore has already reached `wait_value`.
    pub fn wait(&mut self, wait_value: u64) -> bool {
        self.wait_history.push(wait_value);
        self.value >= wait_value
    }
    /// Return the current semaphore value.
    pub fn current_value(&self) -> u64 {
        self.value
    }
    /// Number of times the semaphore has been signalled.
    pub fn signal_count(&self) -> usize {
        self.signal_history.len()
    }
}
/// GPU memory bandwidth model.
///
/// Estimates the effective bandwidth and the roofline-model bound for a kernel.
#[derive(Debug, Clone)]
pub struct MemoryBandwidthModel {
    /// Peak memory bandwidth in GB/s.
    pub peak_bandwidth_gbs: f64,
    /// Peak compute throughput in GFLOP/s.
    pub peak_compute_gflops: f64,
}
impl MemoryBandwidthModel {
    /// Create a model for a mid-range discrete GPU.
    pub fn mid_range() -> Self {
        Self {
            peak_bandwidth_gbs: 480.0,
            peak_compute_gflops: 10000.0,
        }
    }
    /// Compute the arithmetic intensity (FLOP/byte) of a kernel.
    ///
    /// `flops` – total floating-point operations.
    /// `bytes_accessed` – total bytes read/written.
    pub fn arithmetic_intensity(flops: f64, bytes_accessed: f64) -> f64 {
        if bytes_accessed < 1e-30 {
            f64::INFINITY
        } else {
            flops / bytes_accessed
        }
    }
    /// Roofline performance estimate (GFLOP/s) given arithmetic intensity.
    pub fn roofline_performance(&self, arithmetic_intensity: f64) -> f64 {
        let bw_bound = arithmetic_intensity * self.peak_bandwidth_gbs;
        bw_bound.min(self.peak_compute_gflops)
    }
    /// Estimated kernel execution time in milliseconds.
    ///
    /// `flops` – total FLOPs in the kernel.
    /// `bytes_accessed` – total bytes transferred.
    pub fn estimated_runtime_ms(&self, flops: f64, bytes_accessed: f64) -> f64 {
        let intensity = Self::arithmetic_intensity(flops, bytes_accessed);
        let perf_gflops = self.roofline_performance(intensity);
        if perf_gflops < 1e-30 {
            return f64::INFINITY;
        }
        (flops / (perf_gflops * 1e9)) * 1e3
    }
    /// Whether the kernel is bandwidth-bound or compute-bound.
    ///
    /// Returns `true` if bandwidth-bound (arithmetic intensity below the ridge point).
    pub fn is_bandwidth_bound(&self, arithmetic_intensity: f64) -> bool {
        let ridge_point = self.peak_compute_gflops / self.peak_bandwidth_gbs;
        arithmetic_intensity < ridge_point
    }
}
/// A binding entry associating a buffer with a binding index and usage.
#[derive(Debug, Clone, Copy)]
pub struct BufferBinding {
    /// Binding index in the shader (e.g. @binding(0)).
    pub binding: u32,
    /// The buffer to bind.
    pub buffer_id: BufferId,
    /// Usage of the buffer in this binding.
    pub usage: BufferUsage,
}
impl BufferBinding {
    /// Create a new buffer binding.
    pub fn new(binding: u32, buffer_id: BufferId, usage: BufferUsage) -> Self {
        Self {
            binding,
            buffer_id,
            usage,
        }
    }
    /// Shorthand for a read-only binding.
    pub fn read(binding: u32, buffer_id: BufferId) -> Self {
        Self::new(binding, buffer_id, BufferUsage::ReadOnly)
    }
    /// Shorthand for a write-only binding.
    pub fn write(binding: u32, buffer_id: BufferId) -> Self {
        Self::new(binding, buffer_id, BufferUsage::WriteOnly)
    }
    /// Shorthand for a read-write binding.
    pub fn read_write(binding: u32, buffer_id: BufferId) -> Self {
        Self::new(binding, buffer_id, BufferUsage::ReadWrite)
    }
    /// Shorthand for a uniform binding.
    pub fn uniform(binding: u32, buffer_id: BufferId) -> Self {
        Self::new(binding, buffer_id, BufferUsage::Uniform)
    }
}
/// Newtype handle for GPU buffers in the dispatcher model.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct BufferId(pub u32);
