//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::{
    BufferBinding, BufferHandle, BufferId, BufferUsage, PipelineBarrier, WarpDivergenceRecord,
};

/// Trait for a compute backend (GPU or CPU fallback).
pub trait ComputeBackend {
    /// Human-readable name of this backend.
    fn name(&self) -> &str;
    /// Allocate a buffer that can hold `size` f64 elements.
    fn create_buffer(&self, size: usize) -> BufferHandle;
    /// Write `data` into the buffer referenced by `handle`.
    fn write_buffer(&self, handle: BufferHandle, data: &[f64]);
    /// Read the full contents of the buffer referenced by `handle`.
    fn read_buffer(&self, handle: BufferHandle) -> Vec<f64>;
    /// Dispatch a compute kernel over `work_size` work items.
    fn dispatch(&self, kernel: &dyn ComputeKernel, work_size: usize);
}
/// Trait for a compute kernel that can be dispatched on a backend.
pub trait ComputeKernel {
    /// Human-readable name of this kernel.
    fn name(&self) -> &str;
    /// Execute the kernel over `work_size` work items.
    ///
    /// * `inputs`  – read-only input slices.
    /// * `outputs` – mutable output vectors (pre-allocated by the caller).
    /// * `work_size` – number of logical work items.
    fn execute(&self, inputs: &[&[f64]], outputs: &mut [Vec<f64>], work_size: usize);
}
/// Compute the number of workgroups needed for a 1D dispatch.
pub fn compute_num_workgroups(total_items: u32, workgroup_size: u32) -> u32 {
    total_items.div_ceil(workgroup_size)
}
/// Compute workgroup counts for a 3D dispatch.
pub fn compute_num_workgroups_3d(total: [u32; 3], workgroup_size: [u32; 3]) -> [u32; 3] {
    [
        total[0].div_ceil(workgroup_size[0]),
        total[1].div_ceil(workgroup_size[1]),
        total[2].div_ceil(workgroup_size[2]),
    ]
}
/// Determine the required pipeline barrier between two kernel passes.
///
/// If the output buffers of pass A overlap with the input buffers of pass B,
/// a read-after-write barrier is needed.
pub fn required_barrier(
    pass_a_outputs: &[BufferId],
    pass_b_inputs: &[BufferId],
) -> PipelineBarrier {
    let overlap = pass_a_outputs.iter().any(|out| pass_b_inputs.contains(out));
    if overlap {
        PipelineBarrier::StorageReadAfterWrite
    } else {
        PipelineBarrier::None
    }
}
/// Detect whether any buffers in a pass alias the same storage.
///
/// Two bindings alias if they reference the same `BufferId` with incompatible
/// usages (e.g., one is write and the other is read in the same pass).
pub fn detect_aliasing(bindings: &[BufferBinding]) -> Vec<(u32, u32)> {
    let mut conflicts = Vec::new();
    for i in 0..bindings.len() {
        for j in (i + 1)..bindings.len() {
            if bindings[i].buffer_id == bindings[j].buffer_id {
                let write_i = matches!(
                    bindings[i].usage,
                    BufferUsage::WriteOnly | BufferUsage::ReadWrite
                );
                let read_j = matches!(
                    bindings[j].usage,
                    BufferUsage::ReadOnly | BufferUsage::ReadWrite
                );
                let write_j = matches!(
                    bindings[j].usage,
                    BufferUsage::WriteOnly | BufferUsage::ReadWrite
                );
                let read_i = matches!(
                    bindings[i].usage,
                    BufferUsage::ReadOnly | BufferUsage::ReadWrite
                );
                if write_i && read_j || write_j && read_i {
                    conflicts.push((bindings[i].binding, bindings[j].binding));
                }
            }
        }
    }
    conflicts
}
/// Simulate warp divergence by analysing a boolean predicate over work items.
///
/// Groups work items into warps and checks if all threads take the same branch.
/// Returns a [`WarpDivergenceRecord`].
pub fn analyse_warp_divergence(predicates: &[bool], warp_size: usize) -> WarpDivergenceRecord {
    if predicates.is_empty() || warp_size == 0 {
        return WarpDivergenceRecord::default();
    }
    let mut total = 0u64;
    let mut divergent = 0u64;
    let n_warps = predicates.len().div_ceil(warp_size);
    for w in 0..n_warps {
        let start = w * warp_size;
        let end = (start + warp_size).min(predicates.len());
        let slice = &predicates[start..end];
        total += 1;
        let all_true = slice.iter().all(|&v| v);
        let all_false = slice.iter().all(|&v| !v);
        if !all_true && !all_false {
            divergent += 1;
        }
    }
    WarpDivergenceRecord {
        total_branches: total,
        divergent_branches: divergent,
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::CpuBackend;
    use crate::compute::ComputeDispatcher;
    use crate::compute::ComputePass;
    use crate::compute::GpuBuffer;
    use crate::compute::GpuCommand;
    use crate::compute::GpuCommandEncoder;
    use crate::compute::GpuError;
    use crate::compute::KernelSpec;
    use crate::compute::MemoryBandwidthModel;
    use crate::compute::OccupancyModel;
    use crate::compute::ResourceLifecycle;
    use crate::compute::TimelineSemaphore;
    #[test]
    fn cpu_backend_buffer_roundtrip() {
        let backend = CpuBackend::new();
        let buf = backend.create_buffer(4);
        backend.write_buffer(buf, &[1.0, 2.0, 3.0, 4.0]);
        let data = backend.read_buffer(buf);
        assert_eq!(data, vec![1.0, 2.0, 3.0, 4.0]);
    }
    #[test]
    fn dispatcher_buffer_write_read_roundtrip() {
        let mut d = ComputeDispatcher::new();
        let id = d.create_buffer(5, None);
        d.write_buffer(id, &[1.0, 2.0, 3.0, 4.0, 5.0]).unwrap();
        let out = d.read_buffer(id).unwrap();
        assert_eq!(out, vec![1.0, 2.0, 3.0, 4.0, 5.0]);
    }
    #[test]
    fn dispatcher_buffer_initial_data() {
        let mut d = ComputeDispatcher::new();
        let id = d.create_buffer(3, Some(&[10.0, 20.0, 30.0]));
        let out = d.read_buffer(id).unwrap();
        assert_eq!(out, vec![10.0, 20.0, 30.0]);
    }
    #[test]
    fn dispatcher_invalid_buffer_read_errors() {
        let d = ComputeDispatcher::new();
        let bad_id = BufferId(99);
        assert_eq!(d.read_buffer(bad_id), Err(GpuError::InvalidBuffer(bad_id)));
    }
    #[test]
    fn dispatch_map_identity() {
        let mut d = ComputeDispatcher::new();
        let src = d.create_buffer(4, Some(&[1.0, 2.0, 3.0, 4.0]));
        let dst = d.create_buffer(4, None);
        d.dispatch_map(src, dst, |x| x).unwrap();
        assert_eq!(d.read_buffer(dst).unwrap(), vec![1.0, 2.0, 3.0, 4.0]);
    }
    #[test]
    fn dispatch_map_scale_by_two() {
        let mut d = ComputeDispatcher::new();
        let src = d.create_buffer(3, Some(&[1.0, 2.0, 3.0]));
        let dst = d.create_buffer(3, None);
        d.dispatch_map(src, dst, |x| x * 2.0).unwrap();
        assert_eq!(d.read_buffer(dst).unwrap(), vec![2.0, 4.0, 6.0]);
    }
    #[test]
    fn dispatch_reduce_sum() {
        let mut d = ComputeDispatcher::new();
        let id = d.create_buffer(5, Some(&[1.0, 2.0, 3.0, 4.0, 5.0]));
        let sum = d.dispatch_reduce(id, |a, b| a + b).unwrap();
        assert!((sum - 15.0).abs() < 1e-12);
    }
    #[test]
    fn dispatch_reduce_max() {
        let mut d = ComputeDispatcher::new();
        let id = d.create_buffer(5, Some(&[3.0, 1.0, 7.0, 2.0, 5.0]));
        let max = d.dispatch_reduce(id, f64::max).unwrap();
        assert!((max - 7.0).abs() < 1e-12);
    }
    #[test]
    fn dispatch_reduce_empty_errors() {
        let mut d = ComputeDispatcher::new();
        let id = d.create_buffer(0, None);
        assert_eq!(
            d.dispatch_reduce(id, |a, b| a + b),
            Err(GpuError::EmptyBuffer)
        );
    }
    #[test]
    fn sph_density_single_particle_self_contribution_positive() {
        let mut d = ComputeDispatcher::new();
        let pos = d.create_buffer(3, Some(&[0.0, 0.0, 0.0]));
        let mass = d.create_buffer(1, Some(&[1.0]));
        let out = d.create_buffer(1, None);
        d.dispatch_sph_density(pos, mass, 1.0, out).unwrap();
        let density = d.read_buffer(out).unwrap();
        assert_eq!(density.len(), 1);
        assert!((density[0] - 1.0).abs() < 1e-12);
    }
    #[test]
    fn sph_density_two_particles_within_kernel_positive() {
        let mut d = ComputeDispatcher::new();
        let pos = d.create_buffer(6, Some(&[0.0, 0.0, 0.0, 0.5, 0.0, 0.0]));
        let mass = d.create_buffer(2, Some(&[1.0, 1.0]));
        let out = d.create_buffer(2, None);
        d.dispatch_sph_density(pos, mass, 2.0, out).unwrap();
        let density = d.read_buffer(out).unwrap();
        assert_eq!(density.len(), 2);
        assert!(
            density[0] > 0.0,
            "density[0] should be positive: {}",
            density[0]
        );
        assert!(
            density[1] > 0.0,
            "density[1] should be positive: {}",
            density[1]
        );
    }
    #[test]
    fn sph_density_particles_outside_kernel_zero_cross_contribution() {
        let mut d = ComputeDispatcher::new();
        let pos = d.create_buffer(6, Some(&[0.0, 0.0, 0.0, 100.0, 0.0, 0.0]));
        let mass = d.create_buffer(2, Some(&[1.0, 1.0]));
        let out = d.create_buffer(2, None);
        d.dispatch_sph_density(pos, mass, 1.0, out).unwrap();
        let density = d.read_buffer(out).unwrap();
        assert!((density[0] - 1.0).abs() < 1e-12);
        assert!((density[1] - 1.0).abs() < 1e-12);
    }
    #[test]
    fn kernel_spec_creation() {
        let b0 = BufferId(0);
        let b1 = BufferId(1);
        let spec = KernelSpec::new("sph_density", 64, vec![b0, b1]);
        assert_eq!(spec.name, "sph_density");
        assert_eq!(spec.workgroup_size, [64, 1, 1]);
        assert_eq!(spec.buffer_bindings.len(), 2);
    }
    #[test]
    fn gpu_buffer_new_zeros() {
        let buf = GpuBuffer::new(8);
        assert_eq!(buf.size, 8);
        assert!(buf.data.iter().all(|&v| v == 0.0));
    }
    #[test]
    fn test_buffer_binding_shorthands() {
        let id = BufferId(5);
        let br = BufferBinding::read(0, id);
        assert_eq!(br.usage, BufferUsage::ReadOnly);
        let bw = BufferBinding::write(1, id);
        assert_eq!(bw.usage, BufferUsage::WriteOnly);
        let brw = BufferBinding::read_write(2, id);
        assert_eq!(brw.usage, BufferUsage::ReadWrite);
        let bu = BufferBinding::uniform(3, id);
        assert_eq!(bu.usage, BufferUsage::Uniform);
    }
    #[test]
    fn test_kernel_spec_3d_workgroup() {
        let spec = KernelSpec::with_workgroup_3d("test", [8, 8, 4], vec![]);
        assert_eq!(spec.workgroup_size, [8, 8, 4]);
        assert_eq!(spec.threads_per_workgroup(), 256);
    }
    #[test]
    fn test_kernel_spec_num_workgroups() {
        let spec = KernelSpec::new("test", 64, vec![]);
        assert_eq!(spec.num_workgroups_x(100), 2);
        assert_eq!(spec.num_workgroups_x(64), 1);
        assert_eq!(spec.num_workgroups_x(65), 2);
    }
    #[test]
    fn test_gpu_buffer_fill_and_clear() {
        let mut buf = GpuBuffer::new(5);
        buf.fill(42.0);
        assert!(buf.data.iter().all(|&v| (v - 42.0).abs() < 1e-12));
        buf.clear();
        assert!(buf.data.iter().all(|&v| v == 0.0));
    }
    #[test]
    fn test_gpu_buffer_byte_size() {
        let buf = GpuBuffer::new(10);
        assert_eq!(buf.byte_size(), 80);
    }
    #[test]
    fn test_gpu_buffer_as_slice() {
        let buf = GpuBuffer::from_data(vec![1.0, 2.0, 3.0]);
        assert_eq!(buf.as_slice(), &[1.0, 2.0, 3.0]);
    }
    #[test]
    fn test_cpu_backend_num_buffers() {
        let backend = CpuBackend::new();
        assert_eq!(backend.num_buffers(), 0);
        backend.create_buffer(10);
        assert_eq!(backend.num_buffers(), 1);
        backend.create_buffer(5);
        assert_eq!(backend.num_buffers(), 2);
    }
    #[test]
    fn test_cpu_backend_total_elements() {
        let backend = CpuBackend::new();
        backend.create_buffer(10);
        backend.create_buffer(5);
        assert_eq!(backend.total_elements(), 15);
    }
    #[test]
    fn test_dispatcher_num_buffers() {
        let mut d = ComputeDispatcher::new();
        assert_eq!(d.num_buffers(), 0);
        d.create_buffer(5, None);
        assert_eq!(d.num_buffers(), 1);
    }
    #[test]
    fn test_dispatcher_has_buffer() {
        let mut d = ComputeDispatcher::new();
        let id = d.create_buffer(5, None);
        assert!(d.has_buffer(id));
        assert!(!d.has_buffer(BufferId(999)));
    }
    #[test]
    fn test_dispatcher_buffer_size() {
        let mut d = ComputeDispatcher::new();
        let id = d.create_buffer(7, None);
        assert_eq!(d.buffer_size(id).unwrap(), 7);
    }
    #[test]
    fn test_dispatcher_destroy_buffer() {
        let mut d = ComputeDispatcher::new();
        let id = d.create_buffer(5, None);
        assert!(d.has_buffer(id));
        d.destroy_buffer(id).unwrap();
        assert!(!d.has_buffer(id));
    }
    #[test]
    fn test_dispatcher_destroy_invalid_buffer_errors() {
        let mut d = ComputeDispatcher::new();
        assert_eq!(
            d.destroy_buffer(BufferId(42)),
            Err(GpuError::InvalidBuffer(BufferId(42)))
        );
    }
    #[test]
    fn test_dispatcher_copy_buffer() {
        let mut d = ComputeDispatcher::new();
        let src = d.create_buffer(3, Some(&[1.0, 2.0, 3.0]));
        let dst = d.create_buffer(3, None);
        d.copy_buffer(src, dst).unwrap();
        assert_eq!(d.read_buffer(dst).unwrap(), vec![1.0, 2.0, 3.0]);
    }
    #[test]
    fn test_dispatcher_copy_buffer_size_mismatch() {
        let mut d = ComputeDispatcher::new();
        let src = d.create_buffer(3, Some(&[1.0, 2.0, 3.0]));
        let dst = d.create_buffer(5, None);
        assert!(d.copy_buffer(src, dst).is_err());
    }
    #[test]
    fn test_dispatch_map_indexed() {
        let mut d = ComputeDispatcher::new();
        let src = d.create_buffer(4, Some(&[10.0, 20.0, 30.0, 40.0]));
        let dst = d.create_buffer(4, None);
        d.dispatch_map_indexed(src, dst, |i, x| x + i as f64)
            .unwrap();
        assert_eq!(d.read_buffer(dst).unwrap(), vec![10.0, 21.0, 32.0, 43.0]);
    }
    #[test]
    fn test_dispatch_zip_map() {
        let mut d = ComputeDispatcher::new();
        let a = d.create_buffer(3, Some(&[1.0, 2.0, 3.0]));
        let b = d.create_buffer(3, Some(&[10.0, 20.0, 30.0]));
        let out = d.create_buffer(3, None);
        d.dispatch_zip_map(a, b, out, |x, y| x + y).unwrap();
        assert_eq!(d.read_buffer(out).unwrap(), vec![11.0, 22.0, 33.0]);
    }
    #[test]
    fn test_compute_pass_recording() {
        let mut pass = ComputePass::new();
        assert_eq!(pass.num_commands(), 0);
        pass.dispatch("density", 1000);
        pass.dispatch("force", 1000);
        pass.dispatch("integrate", 1000);
        assert_eq!(pass.num_commands(), 3);
        assert_eq!(pass.total_work_items(), 3000);
        assert_eq!(pass.commands()[0].0, "density");
        assert_eq!(pass.commands()[1].1, 1000);
    }
    #[test]
    fn test_compute_pass_clear() {
        let mut pass = ComputePass::new();
        pass.dispatch("test", 100);
        assert_eq!(pass.num_commands(), 1);
        pass.clear();
        assert_eq!(pass.num_commands(), 0);
    }
    #[test]
    fn test_resource_lifecycle_tracking() {
        let mut lifecycle = ResourceLifecycle::new();
        assert!(lifecycle.is_empty());
        let id = BufferId(0);
        lifecycle.record_create(id, 100);
        lifecycle.record_write(id);
        lifecycle.record_write(id);
        lifecycle.record_read(id);
        assert_eq!(lifecycle.len(), 4);
        assert_eq!(lifecycle.count_writes(id), 2);
        assert_eq!(lifecycle.count_reads(id), 1);
    }
    #[test]
    fn test_resource_lifecycle_clear() {
        let mut lifecycle = ResourceLifecycle::new();
        lifecycle.record_create(BufferId(0), 10);
        lifecycle.clear();
        assert!(lifecycle.is_empty());
    }
    #[test]
    fn test_compute_num_workgroups() {
        assert_eq!(compute_num_workgroups(100, 64), 2);
        assert_eq!(compute_num_workgroups(64, 64), 1);
        assert_eq!(compute_num_workgroups(1, 64), 1);
    }
    #[test]
    fn test_compute_num_workgroups_3d() {
        let wg = compute_num_workgroups_3d([100, 100, 100], [8, 8, 8]);
        assert_eq!(wg, [13, 13, 13]);
    }
    #[test]
    fn test_gpu_error_display() {
        let e = GpuError::InvalidBuffer(BufferId(5));
        assert!(format!("{e}").contains("5"));
        let e2 = GpuError::SizeMismatch {
            expected: 10,
            got: 5,
        };
        assert!(format!("{e2}").contains("10"));
        let e3 = GpuError::EmptyBuffer;
        assert!(format!("{e3}").contains("empty"));
        let e4 = GpuError::NotFound("test".to_string());
        assert!(format!("{e4}").contains("test"));
    }
    #[test]
    fn test_command_encoder_basic() {
        let mut enc = GpuCommandEncoder::new("test_pass");
        assert_eq!(enc.label(), "test_pass");
        assert_eq!(enc.command_count(), 0);
        enc.dispatch_compute("density", [64, 1, 1]);
        enc.dispatch_compute("force", [64, 1, 1]);
        enc.insert_barrier(PipelineBarrier::StorageReadAfterWrite);
        assert_eq!(enc.command_count(), 3);
    }
    #[test]
    fn test_command_encoder_reset() {
        let mut enc = GpuCommandEncoder::new("enc");
        enc.dispatch_compute("k", [1, 1, 1]);
        enc.reset();
        assert_eq!(enc.command_count(), 0);
    }
    #[test]
    fn test_command_encoder_submit_copies() {
        let mut enc = GpuCommandEncoder::new("enc");
        let mut d = ComputeDispatcher::new();
        let src = d.create_buffer(3, Some(&[1.0, 2.0, 3.0]));
        let dst = d.create_buffer(3, None);
        enc.copy_buffer(src, dst, 3);
        enc.submit(&mut d).unwrap();
        assert_eq!(d.read_buffer(dst).unwrap(), vec![1.0, 2.0, 3.0]);
    }
    #[test]
    fn test_command_encoder_push_constant() {
        let mut enc = GpuCommandEncoder::new("enc");
        enc.push_constant("dt", 0.001);
        assert_eq!(enc.command_count(), 1);
        match &enc.commands()[0] {
            GpuCommand::PushConstant { name, value } => {
                assert_eq!(name, "dt");
                assert!((value - 0.001).abs() < 1e-15);
            }
            _ => panic!("expected PushConstant"),
        }
    }
    #[test]
    fn test_required_barrier_overlap() {
        let a_out = vec![BufferId(0), BufferId(1)];
        let b_in = vec![BufferId(1), BufferId(2)];
        let barrier = required_barrier(&a_out, &b_in);
        assert_eq!(barrier, PipelineBarrier::StorageReadAfterWrite);
    }
    #[test]
    fn test_required_barrier_no_overlap() {
        let a_out = vec![BufferId(0)];
        let b_in = vec![BufferId(5)];
        let barrier = required_barrier(&a_out, &b_in);
        assert_eq!(barrier, PipelineBarrier::None);
    }
    #[test]
    fn test_detect_aliasing_conflict() {
        let bindings = vec![
            BufferBinding::write(0, BufferId(10)),
            BufferBinding::read(1, BufferId(10)),
        ];
        let conflicts = detect_aliasing(&bindings);
        assert!(!conflicts.is_empty(), "should detect aliasing conflict");
    }
    #[test]
    fn test_detect_aliasing_no_conflict() {
        let bindings = vec![
            BufferBinding::read(0, BufferId(10)),
            BufferBinding::read(1, BufferId(11)),
        ];
        let conflicts = detect_aliasing(&bindings);
        assert!(conflicts.is_empty(), "no conflict expected");
    }
    #[test]
    fn test_detect_aliasing_same_buffer_two_reads() {
        let bindings = vec![
            BufferBinding::read(0, BufferId(5)),
            BufferBinding::read(1, BufferId(5)),
        ];
        let conflicts = detect_aliasing(&bindings);
        assert!(conflicts.is_empty());
    }
    #[test]
    fn test_timeline_semaphore_signal_and_wait() {
        let mut sem = TimelineSemaphore::new();
        assert_eq!(sem.current_value(), 0);
        sem.signal(1);
        assert_eq!(sem.current_value(), 1);
        assert!(sem.wait(1));
        assert!(!sem.wait(2));
        sem.signal(3);
        assert!(sem.wait(3));
        assert_eq!(sem.signal_count(), 2);
    }
    #[test]
    fn test_timeline_semaphore_default() {
        let sem = TimelineSemaphore::default();
        assert_eq!(sem.current_value(), 0);
    }
    #[test]
    fn test_occupancy_full_when_unconstrained() {
        let model = OccupancyModel::mid_range();
        let occ = model.estimate_occupancy(64, 0, 32);
        assert!(
            occ > 0.5,
            "occupancy should be high for small workgroup, got {occ}"
        );
    }
    #[test]
    fn test_occupancy_limited_by_shared_memory() {
        let model = OccupancyModel::mid_range();
        let occ = model.estimate_occupancy(64, model.shared_mem_per_cu, 1);
        let occ_limited = model.estimate_occupancy(64, model.shared_mem_per_cu / 2, 1);
        assert!(
            occ <= occ_limited,
            "more smem usage should give lower or equal occupancy"
        );
    }
    #[test]
    fn test_occupancy_bounded_to_one() {
        let model = OccupancyModel::mid_range();
        let occ = model.estimate_occupancy(1, 0, 0);
        assert!((0.0..=1.0).contains(&occ));
    }
    #[test]
    fn test_peak_gflops_positive() {
        let model = OccupancyModel::mid_range();
        let gflops = model.peak_gflops(1500.0);
        assert!(gflops > 0.0);
    }
    #[test]
    fn test_warp_divergence_none() {
        let predicates = vec![true; 32];
        let rec = analyse_warp_divergence(&predicates, 32);
        assert_eq!(rec.divergent_branches, 0);
        assert!((rec.divergence_rate()).abs() < 1e-12);
    }
    #[test]
    fn test_warp_divergence_full() {
        let predicates: Vec<bool> = (0..32).map(|i| i % 2 == 0).collect();
        let rec = analyse_warp_divergence(&predicates, 32);
        assert_eq!(rec.divergent_branches, 1);
        assert!((rec.divergence_rate() - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_warp_divergence_penalty() {
        let rec = WarpDivergenceRecord {
            total_branches: 10,
            divergent_branches: 5,
        };
        let penalty = rec.performance_penalty(32);
        assert!(
            penalty > 1.0 && penalty < 2.0,
            "penalty should be > 1, got {penalty}"
        );
    }
    #[test]
    fn test_warp_divergence_empty() {
        let rec = analyse_warp_divergence(&[], 32);
        assert_eq!(rec.total_branches, 0);
        assert!((rec.divergence_rate()).abs() < 1e-12);
    }
    #[test]
    fn test_memory_bandwidth_arithmetic_intensity() {
        let intensity = MemoryBandwidthModel::arithmetic_intensity(1000.0, 100.0);
        assert!((intensity - 10.0).abs() < 1e-12);
    }
    #[test]
    fn test_memory_bandwidth_zero_bytes() {
        let intensity = MemoryBandwidthModel::arithmetic_intensity(100.0, 0.0);
        assert!(intensity.is_infinite());
    }
    #[test]
    fn test_roofline_bandwidth_bound() {
        let model = MemoryBandwidthModel::mid_range();
        let perf = model.roofline_performance(0.1);
        let expected = 0.1 * model.peak_bandwidth_gbs;
        assert!(
            (perf - expected).abs() < 1e-6,
            "bandwidth-bound perf mismatch"
        );
    }
    #[test]
    fn test_roofline_compute_bound() {
        let model = MemoryBandwidthModel::mid_range();
        let perf = model.roofline_performance(1e9);
        assert!((perf - model.peak_compute_gflops).abs() < 1e-6);
    }
    #[test]
    fn test_is_bandwidth_bound() {
        let model = MemoryBandwidthModel::mid_range();
        let ridge = model.peak_compute_gflops / model.peak_bandwidth_gbs;
        assert!(model.is_bandwidth_bound(ridge * 0.5));
        assert!(!model.is_bandwidth_bound(ridge * 2.0));
    }
    #[test]
    fn test_estimated_runtime_ms_positive() {
        let model = MemoryBandwidthModel::mid_range();
        let t = model.estimated_runtime_ms(1e12, 1e9);
        assert!(t > 0.0 && t.is_finite());
    }
    #[test]
    fn test_reduction_tree_sum() {
        let mut d = ComputeDispatcher::new();
        let buf = d.create_buffer(4, Some(&[1.0, 2.0, 3.0, 4.0]));
        let result = d.dispatch_reduction_tree(buf).unwrap();
        assert!(
            (result - 10.0).abs() < 1e-12,
            "sum should be 10, got {result}"
        );
    }
    #[test]
    fn test_reduction_tree_empty() {
        let mut d = ComputeDispatcher::new();
        let buf = d.create_buffer(0, Some(&[]));
        let result = d.dispatch_reduction_tree(buf).unwrap();
        assert_eq!(result, 0.0);
    }
    #[test]
    fn test_reduction_tree_single_element() {
        let mut d = ComputeDispatcher::new();
        let buf = d.create_buffer(1, Some(&[42.0]));
        let result = d.dispatch_reduction_tree(buf).unwrap();
        assert!((result - 42.0).abs() < 1e-12);
    }
    #[test]
    fn test_reduction_tree_power_of_two() {
        let data: Vec<f64> = (1..=8).map(|x| x as f64).collect();
        let mut d = ComputeDispatcher::new();
        let buf = d.create_buffer(8, Some(&data));
        let result = d.dispatch_reduction_tree(buf).unwrap();
        assert!((result - 36.0).abs() < 1e-12, "1+2+…+8=36, got {result}");
    }
    #[test]
    fn test_inclusive_scan_basic() {
        let mut d = ComputeDispatcher::new();
        let buf_in = d.create_buffer(4, Some(&[1.0, 2.0, 3.0, 4.0]));
        let buf_out = d.create_buffer(4, None);
        d.dispatch_inclusive_scan(buf_in, buf_out).unwrap();
        let result = d.read_buffer(buf_out).unwrap();
        let expected = [1.0, 3.0, 6.0, 10.0];
        for (a, b) in result.iter().zip(expected.iter()) {
            assert!((a - b).abs() < 1e-12, "mismatch: {a} vs {b}");
        }
    }
    #[test]
    fn test_inclusive_scan_single() {
        let mut d = ComputeDispatcher::new();
        let buf_in = d.create_buffer(1, Some(&[7.0]));
        let buf_out = d.create_buffer(1, None);
        d.dispatch_inclusive_scan(buf_in, buf_out).unwrap();
        let result = d.read_buffer(buf_out).unwrap();
        assert!((result[0] - 7.0).abs() < 1e-12);
    }
    #[test]
    fn test_radix_sort_basic() {
        let data = vec![5.0, 1.0, 3.0, 2.0, 4.0];
        let mut d = ComputeDispatcher::new();
        let buf = d.create_buffer(5, Some(&data));
        let sorted = d.dispatch_radix_sort(buf).unwrap();
        for w in sorted.windows(2) {
            assert!(w[0] <= w[1], "not sorted: {} > {}", w[0], w[1]);
        }
    }
    #[test]
    fn test_radix_sort_empty() {
        let mut d = ComputeDispatcher::new();
        let buf = d.create_buffer(0, Some(&[]));
        let sorted = d.dispatch_radix_sort(buf).unwrap();
        assert!(sorted.is_empty());
    }
    #[test]
    fn test_radix_sort_already_sorted() {
        let data = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let mut d = ComputeDispatcher::new();
        let buf = d.create_buffer(5, Some(&data));
        let sorted = d.dispatch_radix_sort(buf).unwrap();
        assert_eq!(sorted, data);
    }
    #[test]
    fn test_radix_sort_length_preserved() {
        let data: Vec<f64> = (0..16).map(|i| (16 - i) as f64).collect();
        let mut d = ComputeDispatcher::new();
        let buf = d.create_buffer(16, Some(&data));
        let sorted = d.dispatch_radix_sort(buf).unwrap();
        assert_eq!(sorted.len(), 16);
    }
}
