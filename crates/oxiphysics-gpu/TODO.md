# oxiphysics-gpu TODO

Last updated: 2026-04-06 / v0.1.0

## Phase 1: Foundation
- [x] Define core types and traits (ComputeBackend, ComputeKernel, BufferHandle)
- [x] Implement basic error handling
- [x] Add unit tests

## Phase 2: Core Implementation (CPU backend)
- [x] CpuBackend: full CPU-fallback implementation of ComputeBackend
- [x] Kernel dispatch utilities: dispatch_count, aligned_size, linear_index_3d
- [x] DispatchTimer profiling
- [x] ParticleSystem: position/velocity buffers, neighbor queries
- [x] BVH spatial acceleration (bvh module)
- [x] Cell-list neighbor search (cell_list module)
- [x] SDF compute (sdf_compute module)
- [x] Parallel sort (parallel_sort module)
- [x] Grid reduction (grid_reduce module)
- [x] Flux compute (flux_compute module)
- [x] Sparse GPU operations — CPU side (sparse_gpu module)
- [x] Compute pipeline management (compute_pipeline, pipeline modules)
- [x] Shader registry stubs (shader_registry, shaders modules)
- [x] Neural compute kernels — CPU (neural_compute module)
- [x] Integration tests (2,748 public items, 2,740 tests, 0 stubs)
- [x] Performance benchmarks (basic)

## Phase 3: GPU backends (planned)
- [ ] wgpu backend (v0.2.0)
  - [ ] wgpu device/adapter initialization
  - [ ] Buffer upload/download via wgpu
  - [ ] WGSL compute shaders for particle kernels
  - [ ] wgpu-based BVH traversal
- [ ] CUDA backend via cudarc (v0.3.0)
  - [ ] cudarc device context
  - [ ] CUDA kernel launch wrappers
  - [ ] Unified memory support
- [ ] Benchmark: CPU vs wgpu vs CUDA
- [ ] Extended examples (GPU-accelerated SPH, LBM)
