# oxiphysics-gpu

**Status: [Stable]** — v0.1.3 (2026-06-06)

[![Tests](https://img.shields.io/badge/tests-2811-green)](https://github.com/cool-japan/oxiphysics)
[![docs.rs](https://img.shields.io/docsrs/oxiphysics-gpu)](https://docs.rs/oxiphysics-gpu)

GPU-accelerated compute abstractions for the [OxiPhysics](https://github.com/cool-japan/oxiphysics) engine.

> **Note:** The wgpu compute backend is enabled by default on desktop (the `wgpu-backend` feature) with WGSL kernels for SPH density, LBM D3Q19 BGK, and BVH traversal, each parity-tested against the CPU reference. A transparent Rayon CPU fallback is used when no GPU adapter is available, and for WASM / headless / minimal builds via `--no-default-features`. An optional CUDA backend (the `cuda-backend` feature, cudarc dynamic-loading, runtime CUDA ≥ 12.0) is available pending RTX-class hardware verification.

## Features

- **wgpu compute backend** (default `wgpu-backend` feature): real `wgpu` Instance/Adapter/Device/Queue; WGSL compute kernels for SPH density, LBM D3Q19 BGK (streaming + collision), and BVH traversal, each with GPU buffer readback and a CPU-parity test against the reference path
- **Transparent CPU fallback**: a Rayon CPU path is used when no GPU adapter is available, and for WASM / headless / minimal builds via `--no-default-features` (or the `cpu-only` feature)
- **Optional CUDA backend** (`cuda-backend` feature, cudarc dynamic-loading, runtime CUDA ≥ 12.0): non-default, pending RTX-class hardware verification
- **Compute abstraction**: `ComputeBackend` trait + `CpuBackend` implementation
- **Kernel dispatch**: `ComputeKernel`, `BufferHandle`, `DispatchTimer`; `dispatch_count`, `aligned_size`, `linear_index_3d` utilities
- **Particle system**: `ParticleSystem` — position/velocity buffers, neighbor queries
- **Spatial acceleration**: BVH (bounding volume hierarchy), cell-list neighbor search, SDF compute (`sdf_compute`)
- **Parallel primitives**: parallel sort (`parallel_sort`), grid reduction (`grid_reduce`), flux compute (`flux_compute`)
- **Sparse GPU**: sparse matrix/vector operations on CPU (`sparse_gpu`)
- **Pipeline**: compute pipeline management (`compute_pipeline`, `pipeline`), shader registry (`shader_registry`, `shaders`)
- **Neural compute**: neural network forward-pass compute kernels (`neural_compute`)
- 2,748 public API items, 2,811 tests — 0 stubs

## Key Exports

```rust
use oxiphysics_gpu::{
    ComputeBackend, CpuBackend, ComputeKernel, BufferHandle,
    dispatch_count, aligned_size, linear_index_3d,
    DispatchTimer, ParticleSystem,
};
```

## Roadmap

| Version | Feature |
|---------|---------|
| 0.1.x | CPU backend + wgpu compute backend (SPH / LBM D3Q19 BGK / BVH WGSL kernels, CPU-parity tested) — `wgpu-backend` default on desktop |
| 0.1.x | Optional CUDA backend via cudarc (`cuda-backend` feature) — available, pending RTX-class hardware verification |

## License

Apache-2.0 — Part of the [OxiPhysics](https://github.com/cool-japan/oxiphysics) project by COOLJAPAN OU (Team Kitasan)
