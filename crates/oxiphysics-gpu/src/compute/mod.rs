//! Auto-generated module structure

pub mod buffer_pool;
pub mod computedispatcher_traits;
pub mod computepass_traits;
pub mod cpubackend_traits;
pub mod functions;
pub mod gpuerror_traits;
pub mod resourcelifecycle_traits;
pub mod timelinesemaphore_traits;
pub mod timestamp;
pub mod types;
pub mod wgpu_backend;

pub use buffer_pool::BufferPool;
pub use timestamp::{ComputeDispatchTimer, dispatch_count_for};
pub use wgpu_backend::{WgpuBackend, WgpuBufferHandle, WgpuDeviceInfo, WgpuInitError};

pub mod cuda_backend;
pub use cuda_backend::{CudaBackend, CudaBufferHandle, CudaDeviceInfo, CudaInitError};

// Re-export all types
pub use functions::*;
pub use types::*;
