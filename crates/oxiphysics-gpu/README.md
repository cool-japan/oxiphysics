# oxiphysics-gpu

**Status: [Partial]** — v0.1.0 (2026-04-06)

[![Tests](https://img.shields.io/badge/tests-2740-yellow)](https://github.com/cool-japan/oxiphysics)
[![docs.rs](https://img.shields.io/docsrs/oxiphysics-gpu)](https://docs.rs/oxiphysics-gpu)

GPU-accelerated compute abstractions for the [OxiPhysics](https://github.com/cool-japan/oxiphysics) engine.

> **Note:** v0.1.0 ships a **CPU backend only**. No wgpu or CUDA dependencies are present in this release.
> GPU dispatch (wgpu/CUDA) is planned for v0.2.0.

## Features

- **Compute abstraction**: `ComputeBackend` trait + `CpuBackend` implementation
- **Kernel dispatch**: `ComputeKernel`, `BufferHandle`, `DispatchTimer`; `dispatch_count`, `aligned_size`, `linear_index_3d` utilities
- **Particle system**: `ParticleSystem` — position/velocity buffers, neighbor queries
- **Spatial acceleration**: BVH (bounding volume hierarchy), cell-list neighbor search, SDF compute (`sdf_compute`)
- **Parallel primitives**: parallel sort (`parallel_sort`), grid reduction (`grid_reduce`), flux compute (`flux_compute`)
- **Sparse GPU**: sparse matrix/vector operations on CPU (`sparse_gpu`)
- **Pipeline**: compute pipeline management (`compute_pipeline`, `pipeline`), shader registry (`shader_registry`, `shaders`)
- **Neural compute**: neural network forward-pass compute kernels (`neural_compute`)
- 2,748 public API items, 2,740 tests — 0 stubs

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
| 0.1.0 | CPU backend (current) |
| 0.2.0 | wgpu backend (planned) |
| 0.3.0 | CUDA backend via cudarc (planned) |

## License

Apache-2.0 — Part of the [OxiPhysics](https://github.com/cool-japan/oxiphysics) project by COOLJAPAN OU (Team Kitasan)
