//! # CpuBackend - Trait Implementations
//!
//! This module contains trait implementations for `CpuBackend`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//! - `ComputeBackend`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)
use super::functions::{ComputeBackend, ComputeKernel};
use super::types::{BufferHandle, CpuBackend};

impl Default for CpuBackend {
    fn default() -> Self {
        Self::new()
    }
}

impl ComputeBackend for CpuBackend {
    fn name(&self) -> &str {
        "CPU Fallback"
    }
    fn create_buffer(&self, size: usize) -> BufferHandle {
        let mut bufs = self.buffers.borrow_mut();
        let id = bufs.len();
        bufs.push(vec![0.0; size]);
        BufferHandle(id)
    }
    fn write_buffer(&self, handle: BufferHandle, data: &[f64]) {
        let mut bufs = self.buffers.borrow_mut();
        bufs[handle.0] = data.to_vec();
    }
    fn read_buffer(&self, handle: BufferHandle) -> Vec<f64> {
        let bufs = self.buffers.borrow();
        bufs[handle.0].clone()
    }
    fn dispatch(&self, kernel: &dyn ComputeKernel, work_size: usize) {
        let mut outputs: Vec<Vec<f64>> = Vec::new();
        kernel.execute(&[], &mut outputs, work_size);
    }
}
