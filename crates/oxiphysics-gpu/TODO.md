# oxiphysics-gpu TODO

Last updated: 2026-06-06 / v0.1.3

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
- [x] Integration tests (2,748 public items, 2,811 tests, 0 stubs)
- [x] Performance benchmarks (basic)

## Phase 3: GPU backends (planned)
- [x] wgpu backend skeleton (v0.2.0 — `wgpu_backend` module with `WgpuBackend` struct, feature-gated stub)
  - [x] wgpu device/adapter initialization (stub — see `WgpuBackend::try_new` TODO comment)
  - [x] Buffer upload/download via wgpu (stub with CPU shadow + `write_buffer`/`read_buffer`)
  - [x] WGSL compute shaders for particle kernels (`WGSL_SPH_DENSITY`, `WGSL_PARALLEL_SCAN`)
  - [x] wgpu-based BVH traversal (`WGSL_BVH_TRAVERSAL`)
- [x] CUDA backend skeleton via cudarc (v0.3.0 — `cuda_backend` module with `CudaBackend` struct, `CudaBufferHandle`, `CudaDeviceInfo`, `CudaInitError`, PTX kernel stubs)
  - [x] cudarc device context (stub with `try_new(ordinal)` and TODO comments for real impl)
  - [x] CUDA kernel launch wrappers (`launch`, `register_kernel`, `compile_and_register`)
  - [x] Unified memory support (`alloc_unified` — stub; real impl via `cudarc::alloc_zeros_unified`)
- [x] Benchmark: CPU vs wgpu vs CUDA (`gpu_bench` — `GpuBenchHarness`, SPH density + LBM + scan timing across backends)
- [x] Extended examples (GPU-accelerated SPH, LBM) (`sph_gpu` — WCSPH CPU+wgpu, `lbm_gpu` — D3Q19 BGK lid-driven cavity)

## Phase 4: wgpu Backend Activation (v0.2.0)

> **Goal:** Promote `wgpu-backend` from compiling-stub to a real compute path. Today `WgpuBackend::try_new` holds a TODO comment and operates on CPU `Vec<f64>` shadows; the `dispatch()` function does not dispatch to GPU. Phase 4 turns the scaffold into a working pipeline with measurable speedup over the Rayon CPU baseline.

### 4.1 Real device initialization
- [x] `WgpuBackend::try_new` — real `Instance / Adapter / (Device, Queue)` via `pollster::block_on` (planned 2026-04-24; bundles 4.1–4.1d)
  - **Goal:** `WgpuBackend::try_new()` returns a backend wrapping a real `wgpu::Device` on machines with any adapter, and `WgpuInitError::NoAdapter` on headless CI. `available` flag reflects reality.
  - **Design:** `wgpu::Instance::new(Backends::all())`, `pollster::block_on(instance.request_adapter(HighPerformance))`, `adapter.request_device(...)`, populate `WgpuDeviceInfo` from `adapter.get_info()` + `adapter.limits()`. Remove existing `available: false` short-circuit. `#[cfg(feature = "wgpu-backend")]` retained.
  - **Prerequisites:** Add `wgpu`, `pollster`, `bytemuck` to root `Cargo.toml` `[workspace.dependencies]` (grep-confirmed absent 2026-04-24), then reference via `*.workspace = true` in `crates/oxiphysics-gpu/Cargo.toml` under `wgpu-backend` feature gate.
  - **Files:** `src/compute/wgpu_backend.rs` (struct rewrite, `try_new`), `Cargo.toml` (workspace dep references), root `Cargo.toml` (workspace dep additions)
  - **Tests:** `unit::try_new_with_adapter_succeeds` (skip not fail on headless), `unit::no_adapter_error_typed`, `unit::device_info_populated`
- [x] Capture `AdapterInfo` into `WgpuDeviceInfo` (vendor, device name, backend, limits) (planned 2026-04-24; part of 4.1)
- [x] Remove the `available: false` short-circuit; set `available = true` only on successful adapter request (planned 2026-04-24; part of 4.1)
- [x] Typed `WgpuInitError::NoAdapter` fallback if no adapter found on headless CI (planned 2026-04-24; part of 4.1)

### 4.2 Real buffer pipeline
- [x] Real `wgpu::Buffer` handles (STORAGE | COPY_SRC | COPY_DST; MAP_READ for readback) (planned 2026-04-24; bundles 4.2–4.2e)
  - **Goal:** `WgpuBufferHandle` wraps a real `wgpu::Buffer`; `write_buffer`/`read_buffer` move data via `queue.write_buffer` / staging-buffer `map_async`. f64↔f32 conversion at API boundary. Reusable buffer pool.
  - **Design:** `WgpuBufferHandle { id: BufferId, buffer: Arc<wgpu::Buffer>, size_bytes: u64, usage: BufferUsages, scalar: Scalar }`. `write_buffer` → f64→f32 projection + `queue.write_buffer`. `read_buffer` → staging `MAP_READ|COPY_DST`, `copy_buffer_to_buffer`, `map_async(Read)`, `device.poll(Wait)`. Buffer pool: size-bucketed freelist keyed on `(next_pow2(size_bytes), usage)`, capped at 64/bucket.
  - **Files:** `src/compute/wgpu_backend.rs` (buffer methods), `src/compute/buffer_pool.rs` (NEW, ~250 LoC)
  - **Tests:** `unit::buffer_round_trip_f64`, `unit::buffer_round_trip_f32_via_f64_api` (rel err ≤ 1e-6), `unit::pool_reuses_buckets`; all gated on `try_new().is_ok()`
- [x] `write_buffer` — `queue.write_buffer` into GPU memory; no CPU shadow (planned 2026-04-24; part of 4.2)
- [x] `read_buffer` — staging-buffer copy + `map_async(MapMode::Read)` with proper async fence (planned 2026-04-24; part of 4.2)
- [x] `f32` as on-device scalar; transparently convert from `f64` in the API layer (planned 2026-04-24; part of 4.2)
- [x] Buffer pool / free-list to avoid per-frame allocation churn (planned 2026-04-24; part of 4.2)

### 4.3 Real compute dispatch
- [x] `dispatch` — compile WGSL via `device.create_shader_module`, cached `ComputePipeline` (by shader-hash), bind groups (planned 2026-04-24; bundles 4.3–4.3c)
  - **Goal:** `dispatch(kernel_src, bind_groups, workgroups)` compiles WGSL once, caches `ComputePipeline` by fnv1a hash of source+entry, encodes compute pass, optionally records timestamps.
  - **Design:** `ShaderCache { map: HashMap<u64, Arc<wgpu::ComputePipeline>> }`. `BindGroupSpec { buffers: Vec<(u32, WgpuBufferHandle, BindingType)> }`. `dispatch_count_for(n_items, workgroup_size) -> [u32; 3]` = `(n + wg - 1) / wg` in x-dim. Timestamp queries feature-gated on `wgpu::Features::TIMESTAMP_QUERY` → `QuerySet` size 2, resolve → `DispatchTimer::elapsed_ns()`.
  - **Files:** `src/compute/wgpu_backend.rs` (dispatch, shader cache, bind group builder), `src/compute/timestamp.rs` (NEW, ~180 LoC)
  - **Tests:** `unit::dispatch_count_edge_cases`, `unit::shader_cache_hit`, `integration::timestamp_nonzero` (skip if TIMESTAMP_QUERY unsupported)
- [x] Workgroup size tuning; expose `dispatch_count_for(n_items, workgroup_size)` (planned 2026-04-24; part of 4.3)
- [x] Real fence / timestamp queries via `QuerySet` (feature: `TIMESTAMP_QUERY`) for `DispatchTimer` (planned 2026-04-24; part of 4.3)

### 4.4 Kernel activation (end-to-end tests)
- [x] `WGSL_SPH_DENSITY` — bind particles, dispatch, read back density; smoke test with 4 particles (2026-04-24; `tests/wgpu_kernels.rs`)
  - **Goal:** Each existing WGSL constant wired through real dispatch and validated against CPU reference.
  - **Design:** SPH density: bind positions/mass/h (baked-in constants to avoid uniform-upload API gap), dispatch, non-zero density verified. Parallel scan: N=256 copy-kernel dispatch verifies round-trip. BVH traversal & LBM full-parity: deferred to Phase 5 (need raw-bytes uniform upload path for params).
  - **Files:** `src/sph_gpu.rs`, `src/lbm_gpu.rs`, `src/compute/wgpu_backend.rs` (dispatch wiring), `tests/wgpu_kernels.rs` (NEW — all tests gated on `try_new().is_ok()`, skip-not-fail on headless CI)
  - **Tests:** `test_wgpu_sph_density_dispatch_smoke`, `test_wgpu_parallel_scan_parity`, `test_wgpu_buffer_round_trip`, `test_wgpu_shader_cache_hit`, `test_wgpu_backend_is_available`, `test_wgpu_dispatch_count_for`; all skip-not-fail on headless CI
- [x] `WGSL_PARALLEL_SCAN` — copy-kernel dispatch parity N=256 (2026-04-24; part of 4.4; full Blelloch deferred — `pass` keyword reserved in wgpu 29)
- [x] `WGSL_BVH_TRAVERSAL` — full traversal parity — completed 2026-05-11
  - **Goal:** `BvhGpuTraverser::traverse_rays` runs end-to-end on `WgpuBackendReal`, reusing the device/queue/pipeline across calls, with a CPU-parity test on a 10⁵-leaf BVH proving hit-index equality.
  - **Design:**
    0. **Send-bound audit (pre-design)**: subagent runs `rg -n 'BvhGpuTraverser|par_iter|rayon::scope|spawn' crates/oxiphysics-*/src/ --no-heading` to find every call site. If `BvhGpuTraverser` is *ever* moved into a `rayon::ParallelIterator` closure, use `std::sync::Mutex<WgpuBackendReal>` (or `parking_lot::Mutex` if already a workspace dep) instead of `RefCell` so the type stays `Send + Sync` unconditionally. Document the decision in a comment on the field.
    1. **Refactor `BvhGpuTraverser` to single-init**: move `prim_aabbs_buf`, `prim_indices_buf`, `object_ids_buf` into `BvhGpuState` (uploaded at construction). Stop re-allocating `WgpuBackendReal` inside `traverse_rays_gpu` — use the one stored in `BvhGpuState`. Use the interior-mutability primitive chosen in step 0 (`Mutex<WgpuBackendReal>` by default).
    2. **Per-call rays buffer**: only the rays buffer and results buffer are re-allocated/resized per call.
    3. **Dispatch**: keep using `BVH_TRAVERSAL_WGSL`; call `dispatch_wgsl` with the cached pipeline.
    4. **Typed uniform helper (optional polish)**: add `WgpuBackendReal::create_buffer_uniform_raw(&mut self, data: &[u8]) -> WgpuBufferHandle`.
    5. **Splitrs `bvh.rs` (2245 lines)**: split into `crates/oxiphysics-gpu/src/bvh/{mod.rs, types.rs, cpu.rs, gpu.rs}` using `splitrs`.
    6. **Observability counter**: add `pub(crate) dispatch_count: AtomicU64` to `BvhGpuState`.
  - **Files:** `crates/oxiphysics-gpu/src/bvh.rs` → split into `bvh/{mod.rs,types.rs,cpu.rs,gpu.rs}`; `compute/wgpu_backend.rs`; `tests/wgpu_kernels.rs`; `TODO.md`; root `TODO.md`
  - **Tests:** `test_bvh_gpu_parity_10e5_leaves`, `test_bvh_gpu_traverser_reuses_state` (dispatch_count-based, no timing), `test_bvh_gpu_traverser_send_across_threads` (if parallel call sites found)
  - **Risk:** Mutex poisoning under panic; splitrs layout differences — verify with `cargo check -p oxiphysics-gpu` after split.
- [x] LBM D3Q19 BGK step — streaming + collision kernels, lid-driven cavity smoke test — completed 2026-05-11
  - **Goal:** `LbmSimulation::step()` runs the D3Q19 BGK collision+streaming on `WgpuBackendReal` using `shaders/lbm_bgk_d3q19.wgsl`, with a 32³ lid-driven cavity smoke test asserting mean-velocity / mean-density agreement within 1e-3 over 500 steps.
  - **Design:**
    0. **WGSL indexing verification (pre-design)**: read `crates/oxiphysics-gpu/src/shaders/lbm_bgk_d3q19.wgsl` end-to-end; document exact `params` buffer layout, `f_in`/`f_out` indexing formula, streaming offsets, boundary conventions. Write findings to `/tmp/lbm_wgsl_layout.md`. Do NOT proceed until verified.
    0b. **Round-trip marshalling unit test (pre-dispatch)**: `test_lbm_soa_to_gpu_buffer_roundtrip` — 4×4×4 SoA with unique values; flatten→unflatten; assert reconstruction. CPU-only, no GPU required.
    1. **Replace stub `WgpuBackend` with `WgpuBackendReal`** in `LbmSimulation`; add `Option<LbmGpuState>` field.
    2. **`LbmGpuState`**: owns `WgpuBackendReal`, `f_in_buf`, `f_out_buf`, `params_buf`, ping-pong bool.
    3. **Upload**: SoA `Vec<Vec<f64>>` → flat `Vec<f32>` using exact indexing from step 0.
    4. **Dispatch**: `dispatch_wgsl` with workgroup `[ceil(nx/8), ceil(ny/8), ceil(nz/8)]`.
    5. **Boundary conditions**: verify `omega_bits` `bitcast<f32>` round-trip matches WGSL.
    6. **Ping-pong**: swap `f_in`/`f_out` handles each step (no data copy).
    7. **Readback**: only when macroscopic moments queried.
  - **Files:** `crates/oxiphysics-gpu/src/lbm_gpu.rs`; `shaders/lbm_bgk_d3q19.wgsl` (verify omega_bits); `tests/wgpu_kernels.rs`; `TODO.md`; root `TODO.md`
  - **Tests:** `test_lbm_soa_to_gpu_buffer_roundtrip` (CPU-only, gating), `test_lbm_d3q19_lid_cavity_gpu_vs_cpu` (GPU, skip-with-print if unavailable), `test_lbm_d3q19_gpu_resident_stepping`
  - **Risk:** SoA→flat marshalling off-by-one; omega_bits packing; ping-pong buffer confusion.

### 4.5 Benchmarks & CI
- [x] `gpu_bench` harness updated — `cpu_vs_wgpu_comparison` method added; CPU inclusive-scan vs wgpu copy-dispatch (2026-04-24; `src/gpu_bench.rs`)
  - **Goal:** `GpuBenchHarness::cpu_vs_wgpu_comparison(n)` benchmarks CPU scan vs wgpu dispatch; returns CPU-only when no adapter available.
  - **Files:** `src/gpu_bench.rs` (new method), `Cargo.toml` (default feature promoted)
- [x] Regression test: wgpu SPH ≥ 5× CPU at N = 10⁵ (env-gated: set OXIPHYSICS_RTX_BENCH=1) (planned 2026-05-14, landed 2026-05-14)
  - **Note (2026-05-14):** `tests/wgpu_sph_speedup.rs` smoke test always passes; RTX assertion activates with env var.
- [x] `wgpu-backend` promoted to `default = ["wgpu-backend"]` (2026-04-24; `Cargo.toml`)

### 4.6 Scope NOT in Phase 4 (parked)
- CUDA backend activation (`cuda_backend`) stays at skeleton for v0.3.0; Phase 4 is wgpu-only.
- Custom compute fences / `wgpu::CommandEncoder` orchestration across multi-pass pipelines — move to Phase 5 if needed.
- WebGPU in browser (distinct from desktop wgpu) — handled by `oxiphysics-wasm` Phase 7 demos.
