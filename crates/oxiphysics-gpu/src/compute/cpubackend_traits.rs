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
use crate::{Error, Result};

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
    fn write_buffer(&self, handle: BufferHandle, data: &[f64]) -> Result<()> {
        let mut bufs = self.buffers.borrow_mut();
        let slot = bufs
            .get_mut(handle.0)
            .ok_or(Error::InvalidBufferHandle(handle.0))?;
        slot.clear();
        slot.extend_from_slice(data);
        Ok(())
    }
    fn read_buffer(&self, handle: BufferHandle) -> Result<Vec<f64>> {
        self.buffers
            .borrow()
            .get(handle.0)
            .cloned()
            .ok_or(Error::InvalidBufferHandle(handle.0))
    }
    fn dispatch(
        &self,
        kernel: &dyn ComputeKernel,
        inputs: &[BufferHandle],
        outputs: &[BufferHandle],
        work_size: usize,
    ) -> Result<()> {
        for (position, handle) in outputs.iter().enumerate() {
            if outputs[..position].contains(handle) {
                return Err(Error::general(format!(
                    "buffer {} is bound to more than one output of kernel '{}'",
                    handle.0,
                    kernel.name()
                )));
            }
        }
        // Snapshot the inputs and the pre-sized outputs, then release the borrow
        // before the kernel runs so that it may itself use this backend.
        let (input_data, mut output_data) = {
            let bufs = self.buffers.borrow();
            let fetch = |handle: &BufferHandle| {
                bufs.get(handle.0)
                    .cloned()
                    .ok_or(Error::InvalidBufferHandle(handle.0))
            };
            let input_data = inputs.iter().map(&fetch).collect::<Result<Vec<_>>>()?;
            let output_data = outputs.iter().map(&fetch).collect::<Result<Vec<_>>>()?;
            (input_data, output_data)
        };
        let input_views: Vec<&[f64]> = input_data.iter().map(Vec::as_slice).collect();
        kernel.execute(&input_views, &mut output_data, work_size);
        let mut bufs = self.buffers.borrow_mut();
        for (handle, data) in outputs.iter().zip(output_data) {
            let slot = bufs
                .get_mut(handle.0)
                .ok_or(Error::InvalidBufferHandle(handle.0))?;
            *slot = data;
        }
        Ok(())
    }
}
