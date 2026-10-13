# oxiphysics-gpu TODO

**Version:** 0.1.3 | **Updated:** 2026-06-11 | **Status:** v0.1.x complete — production baseline; forward roadmap active
**Scale:** 63,794 SLoC | 2,811 tests

> **Constraints:** Pure-Rust default features (`wgpu-backend` is the default GPU path). The optional `cuda-backend` feature via cudarc is the approved FFI pattern and stays non-default. Benchmarks/CI are local scripts or `xtask` subcommands only — no new workflow yamls. Every promoted kernel lands Phase-22-style: CPU-parity test + env-gated speedup test (`OXIPHYSICS_RTX_BENCH` pattern), skip-not-fail on headless machines.

## Completed (v0.1.0 – v0.1.3)

### Phase 1: Foundation
- [x] Define core types and traits (ComputeBackend, ComputeKernel, BufferHandle)
- [x] Implement basic error handling
- [x] Add unit tests

### Phase 2: Core Implementation (CPU backend)
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

### Phase 3: GPU backends (shipped — wgpu real; CUDA optional/pending hardware)
- [x] wgpu backend (`wgpu_backend` module — real `WgpuBackendReal` with `wgpu` Instance/Adapter/Device/Queue; the `WgpuBackend` struct remains as the CPU-shadow fallback)
  - [x] wgpu device/adapter initialization (real — `WgpuBackendReal::try_new` / `try_new_async` via `pollster::block_on`)
  - [x] Buffer upload/download via wgpu (real — `queue.write_buffer` + staging-buffer `map_async` readback)
  - [x] WGSL compute shaders for particle kernels (`WGSL_SPH_DENSITY`, `WGSL_PARALLEL_SCAN`)
  - [x] wgpu-based BVH traversal (`WGSL_BVH_TRAVERSAL`)
- [x] CUDA backend via cudarc (`cuda_backend` module with `CudaBackend`, `CudaBufferHandle`, `CudaDeviceInfo`, `CudaInitError`) — behind the optional, non-default `cuda-backend` feature; pending RTX-class hardware verification
  - [x] cudarc device context (`try_new(ordinal)`)
  - [x] CUDA kernel launch wrappers (`launch`, `register_kernel`, `compile_and_register`)
  - [x] Unified memory support (`alloc_unified` via `cudarc::alloc_zeros_unified`)
- [x] Benchmark: CPU vs wgpu vs CUDA (`gpu_bench` — `GpuBenchHarness`, SPH density + LBM + scan timing across backends)
- [x] Extended examples (GPU-accelerated SPH, LBM) (`sph_gpu` — WCSPH CPU+wgpu, `lbm_gpu` — D3Q19 BGK lid-driven cavity)

### Phase 4: wgpu Backend Activation (shipped — root TODO Phase 22 / KF-4, 2026-05-14)

> **Status:** Complete. The `wgpu-backend` feature is a real compute path via `WgpuBackendReal` (real `wgpu::Device`/`Queue`, `queue.write_buffer`, staging readback, cached `ComputePipeline` dispatch) and is enabled by default on desktop. SPH density, LBM D3Q19 BGK, and BVH traversal kernels run end-to-end with CPU-parity tests. The ≥5× speedup regression test is env-gated (`OXIPHYSICS_RTX_BENCH=1`) and requires RTX-class hardware; no measured speedup is claimed here.

#### 4.1 Real device initialization
- [x] `WgpuBackend::try_new` — real `Instance / Adapter / (Device, Queue)` via `pollster::block_on` (planned 2026-04-24; bundles 4.1–4.1d)
  - **Goal:** `WgpuBackend::try_new()` returns a backend wrapping a real `wgpu::Device` on machines with any adapter, and `WgpuInitError::NoAdapter` on headless CI. `available` flag reflects reality.
  - **Design:** `wgpu::Instance::new(Backends::all())`, `pollster::block_on(instance.request_adapter(HighPerformance))`, `adapter.request_device(...)`, populate `WgpuDeviceInfo` from `adapter.get_info()` + `adapter.limits()`. Remove existing `available: false` short-circuit. `#[cfg(feature = "wgpu-backend")]` retained.
  - **Prerequisites:** Add `wgpu`, `pollster`, `bytemuck` to root `Cargo.toml` `[workspace.dependencies]` (grep-confirmed absent 2026-04-24), then reference via `*.workspace = true` in `crates/oxiphysics-gpu/Cargo.toml` under `wgpu-backend` feature gate.
  - **Files:** `src/compute/wgpu_backend.rs` (struct rewrite, `try_new`), `Cargo.toml` (workspace dep references), root `Cargo.toml` (workspace dep additions)
  - **Tests:** `unit::try_new_with_adapter_succeeds` (skip not fail on headless), `unit::no_adapter_error_typed`, `unit::device_info_populated`
- [x] Capture `AdapterInfo` into `WgpuDeviceInfo` (vendor, device name, backend, limits) (planned 2026-04-24; part of 4.1)
- [x] Remove the `available: false` short-circuit; set `available = true` only on successful adapter request (planned 2026-04-24; part of 4.1)
- [x] Typed `WgpuInitError::NoAdapter` fallback if no adapter found on headless CI (planned 2026-04-24; part of 4.1)

#### 4.2 Real buffer pipeline
- [x] Real `wgpu::Buffer` handles (STORAGE | COPY_SRC | COPY_DST; MAP_READ for readback) (planned 2026-04-24; bundles 4.2–4.2e)
  - **Goal:** `WgpuBufferHandle` wraps a real `wgpu::Buffer`; `write_buffer`/`read_buffer` move data via `queue.write_buffer` / staging-buffer `map_async`. f64↔f32 conversion at API boundary. Reusable buffer pool.
  - **Design:** `WgpuBufferHandle { id: BufferId, buffer: Arc<wgpu::Buffer>, size_bytes: u64, usage: BufferUsages, scalar: Scalar }`. `write_buffer` → f64→f32 projection + `queue.write_buffer`. `read_buffer` → staging `MAP_READ|COPY_DST`, `copy_buffer_to_buffer`, `map_async(Read)`, `device.poll(Wait)`. Buffer pool: size-bucketed freelist keyed on `(next_pow2(size_bytes), usage)`, capped at 64/bucket.
  - **Files:** `src/compute/wgpu_backend.rs` (buffer methods), `src/compute/buffer_pool.rs` (NEW, ~250 LoC)
  - **Tests:** `unit::buffer_round_trip_f64`, `unit::buffer_round_trip_f32_via_f64_api` (rel err ≤ 1e-6), `unit::pool_reuses_buckets`; all gated on `try_new().is_ok()`
- [x] `write_buffer` — `queue.write_buffer` into GPU memory; no CPU shadow (planned 2026-04-24; part of 4.2)
- [x] `read_buffer` — staging-buffer copy + `map_async(MapMode::Read)` with proper async fence (planned 2026-04-24; part of 4.2)
- [x] `f32` as on-device scalar; transparently convert from `f64` in the API layer (planned 2026-04-24; part of 4.2)
- [x] Buffer pool / free-list to avoid per-frame allocation churn (planned 2026-04-24; part of 4.2)

#### 4.3 Real compute dispatch
- [x] `dispatch` — compile WGSL via `device.create_shader_module`, cached `ComputePipeline` (by shader-hash), bind groups (planned 2026-04-24; bundles 4.3–4.3c)
  - **Goal:** `dispatch(kernel_src, bind_groups, workgroups)` compiles WGSL once, caches `ComputePipeline` by fnv1a hash of source+entry, encodes compute pass, optionally records timestamps.
  - **Design:** `ShaderCache { map: HashMap<u64, Arc<wgpu::ComputePipeline>> }`. `BindGroupSpec { buffers: Vec<(u32, WgpuBufferHandle, BindingType)> }`. `dispatch_count_for(n_items, workgroup_size) -> [u32; 3]` = `(n + wg - 1) / wg` in x-dim. Timestamp queries feature-gated on `wgpu::Features::TIMESTAMP_QUERY` → `QuerySet` size 2, resolve → `DispatchTimer::elapsed_ns()`.
  - **Files:** `src/compute/wgpu_backend.rs` (dispatch, shader cache, bind group builder), `src/compute/timestamp.rs` (NEW, ~180 LoC)
  - **Tests:** `unit::dispatch_count_edge_cases`, `unit::shader_cache_hit`, `integration::timestamp_nonzero` (skip if TIMESTAMP_QUERY unsupported)
- [x] Workgroup size tuning; expose `dispatch_count_for(n_items, workgroup_size)` (planned 2026-04-24; part of 4.3)
- [x] Real fence / timestamp queries via `QuerySet` (feature: `TIMESTAMP_QUERY`) for `DispatchTimer` (planned 2026-04-24; part of 4.3)

#### 4.4 Kernel activation (end-to-end tests)
- [x] `WGSL_SPH_DENSITY` — bind particles, dispatch, read back density; smoke test with 4 particles (2026-04-24; `tests/wgpu_kernels.rs`)
  - **Goal:** Each existing WGSL constant wired through real dispatch and validated against CPU reference.
  - **Design:** SPH density: bind positions/mass/h (baked-in constants to avoid uniform-upload API gap), dispatch, non-zero density verified. Parallel scan: N=256 copy-kernel dispatch verifies round-trip. BVH traversal & LBM full-parity: deferred to Phase 5 (need raw-bytes uniform upload path for params).
  - **Files:** `src/sph_gpu.rs`, `src/lbm_gpu.rs`, `src/compute/wgpu_backend.rs` (dispatch wiring), `tests/wgpu_kernels.rs` (NEW — all tests gated on `try_new().is_ok()`, skip-not-fail on headless CI)
  - **Tests:** `test_wgpu_sph_density_dispatch_smoke`, `test_wgpu_parallel_scan_parity`, `test_wgpu_buffer_round_trip`, `test_wgpu_shader_cache_hit`, `test_wgpu_backend_is_available`, `test_wgpu_dispatch_count_for`; all skip-not-fail on headless CI
- [x] `WGSL_PARALLEL_SCAN` — copy-kernel dispatch parity N=256 (2026-04-24; part of 4.4; full Blelloch deferred — `pass` keyword reserved in wgpu 29)
- [x] `WGSL_BVH_TRAVERSAL` — full traversal parity — completed 2026-05-11
  - **Goal:** `BvhGpuTraverser::traverse_rays` runs end-to-end on `WgpuBackendReal`, reusing the device/queue/pipeline across calls, with a CPU-parity test on a 10⁵-leaf BVH proving hit-index equality.
  - **Design:** single-init `BvhGpuState` (prim AABBs/indices/object-ids uploaded at construction), per-call rays/results buffers only, cached `BVH_TRAVERSAL_WGSL` pipeline via `dispatch_wgsl`, `Mutex<WgpuBackendReal>` for Send+Sync, `bvh.rs` (2,245 lines) split via splitrs into `bvh/{mod,types,cpu,gpu}.rs`, `dispatch_count: AtomicU64` observability counter.
  - **Tests:** `test_bvh_gpu_parity_10e5_leaves`, `test_bvh_gpu_traverser_reuses_state` (dispatch_count-based, no timing), `test_bvh_gpu_traverser_send_across_threads`
- [x] LBM D3Q19 BGK step — streaming + collision kernels, lid-driven cavity smoke test — completed 2026-05-11
  - **Goal:** `LbmSimulation::step()` runs the D3Q19 BGK collision+streaming on `WgpuBackendReal` using `shaders/lbm_bgk_d3q19.wgsl`, with a 32³ lid-driven cavity smoke test asserting mean-velocity / mean-density agreement within 1e-3 over 500 steps.
  - **Design:** WGSL indexing verified end-to-end first; `LbmGpuState` owns backend + `f_in`/`f_out`/`params` buffers; SoA→flat f32 upload; workgroup `[ceil(nx/8), ceil(ny/8), ceil(nz/8)]`; `omega_bits` bitcast round-trip verified; ping-pong handle swap per step; readback only when macroscopic moments queried.
  - **Tests:** `test_lbm_soa_to_gpu_buffer_roundtrip` (CPU-only, gating), `test_lbm_d3q19_lid_cavity_gpu_vs_cpu` (GPU, skip-with-print if unavailable), `test_lbm_d3q19_gpu_resident_stepping`

#### 4.5 Benchmarks & CI
- [x] `gpu_bench` harness updated — `cpu_vs_wgpu_comparison` method added; CPU inclusive-scan vs wgpu copy-dispatch (2026-04-24; `src/gpu_bench.rs`)
  - **Goal:** `GpuBenchHarness::cpu_vs_wgpu_comparison(n)` benchmarks CPU scan vs wgpu dispatch; returns CPU-only when no adapter available.
  - **Files:** `src/gpu_bench.rs` (new method), `Cargo.toml` (default feature promoted)
- [x] Regression test: wgpu SPH ≥ 5× CPU at N = 10⁵ (env-gated: set OXIPHYSICS_RTX_BENCH=1) (planned 2026-05-14, landed 2026-05-14)
  - **Note (2026-05-14):** `tests/wgpu_sph_speedup.rs` smoke test always passes; RTX assertion activates with env var.
- [x] `wgpu-backend` promoted to `default = ["wgpu-backend"]` (2026-04-24; `Cargo.toml`)

#### 4.6 Scope NOT in Phase 4 (parked at ship time; now scheduled in the forward roadmap)
- CUDA backend activation (`cuda_backend`) stays at skeleton for v0.3.0; Phase 4 is wgpu-only. → now an explicit `[~]` item under v0.3.0 below (hardware-gated).
- Custom compute fences / `wgpu::CommandEncoder` orchestration across multi-pass pipelines — move to Phase 5 if needed. → now under Deferred / research track.
- WebGPU in browser (distinct from desktop wgpu) — handled by `oxiphysics-wasm` Phase 7 demos. → tracked in the oxiphysics-wasm TODO (v0.3.0 "WebGPU compute in browser").

### Known in-code limitation markers (audited 2026-06-11)

Correcting the record: the often-quoted "59 TODOs" figure is a miscount. The audited taxonomy is **3 literal `// TODO` blocks + ~38 "(CPU mock)" module banners + 7 undispatched WGSL templates + 3 CUDA stub-PTX notes**. The CPU mocks are honest, tested reference implementations — not missing code. Full inventory below; forward roadmap items reference these as "(closes marker N)".

1. `src/compute/wgpu_backend.rs:133,245,271` — legacy `WgpuBackend` stub TODOs — obsolete: `WgpuBackendReal` in the same file already does all three
2. `src/compute/cuda_backend.rs:183,234,287` — "(stub — replace with actual nvcc-compiled PTX)" x3
3. `src/gpu_bench.rs:129` — "Try CUDA — this is a stub that always reports unavailable"
4. `src/gpu_rigid.rs` — "GPU-accelerated rigid body batch simulation (CPU mock)" — SAP broadphase + sequential-impulse on CPU
5. `src/gpu_cloth.rs` — cloth springs/bending/XPBD/collision "(CPU mock implementation)"
6. `src/gpu_sort.rs` + `src/parallel_sort.rs` — bitonic/radix sort "(CPU simulation)"; `parallel_sort.rs:1022` "Simulates the GPU per-block histogram + scatter"
7. `src/gpu_reduction.rs` + `src/grid_reduce.rs` — scan/reduce/histogram/compact "warp-level primitives (simulated on CPU)"; `grid_reduce.rs:256` "In a real GPU kernel this would use atomicAdd"
8. `src/gpu_sparse_solver.rs` — CSR SpMV + CG/PCG "(CPU mock)"
9. `src/gpu_fem_assembly.rs`, `src/gpu_md_solver.rs`, `src/gpu_fluid_euler.rs` (MAC grid), `src/gpu_lbm.rs` (D2Q9) — CPU mock backends
10. `src/gpu_sdf.rs`, `src/gpu_voxel.rs`, `src/gpu_mesh_processing.rs`, `src/gpu_ray_tracing.rs`, `src/gpu_thermal.rs`, `src/gpu_particle_system.rs` — CPU mocks
11. `src/scheduler.rs:494–528` — "Simulated async compute queue", placeholder outputs
12. `src/memory.rs:36` — "CPU-side mock for a typed GPU buffer backed by Vec<u8>"
13. `src/shaders/functions.rs` — WGSL templates SPH_DENSITY / INTEGRATE / LBM_BGK_D2Q9 / CELL_LIST / SDF_COMPUTE / SPH_FORCE / BOUNDARY_ENFORCE registered in shader_registry but never dispatched end-to-end
14. `src/compute/wgpu_backend.rs` `WGSL_PARALLEL_SCAN` — copy-kernel only; full Blelloch deferred ("`pass` keyword reserved in wgpu 29")
15. `src/kernels/{rigid,md_force}/types.rs` — "CPU-side mock of the GPU contact-batch dispatch kernel"

## v0.2.0 — Close the CPU-mock gap on the hot path

Sequencing (root TODO Phase 26): reduction lib → radix sort → LBVH → SPH pipeline → cloth → contact solver. Exit criteria: zero literal TODOs in `src/compute/`; 250k-particle SPH dam break ≥60 steps/s on M-series; conformance suite v1 enumerates every dispatched kernel.

### Foundation primitives
- [x] GPU reduction-primitive library (scan/reduce/compact/histogram, real WGSL) — (closes markers 7, 14) (planned 2026-06-11)
  - **Goal:** exclusive scan of 16M u32 <5 ms on M-series; parity vs `gpu_reduction.rs`/`grid_reduce.rs` CPU mocks; replaces copy-kernel `WGSL_PARALLEL_SCAN` (rename identifiers to dodge the `pass` reserved word).
  - **Design:** workgroup Blelloch + device-level scan-of-block-sums; new `src/kernels_wgsl/` dir of `include_str!` shaders dispatched via `WgpuBackendReal::dispatch_wgsl`. Foundation for everything below.
  - **Files:** NEW src/kernels_wgsl/{scan,reduce,compact,histogram}.wgsl (include_str!), NEW src/gpu_primitives.rs (wrappers: exclusive_scan_u32, reduce_sum_f32/reduce_max_f32, compact_u32, histogram_u32 via WgpuBackendReal); gpu_reduction.rs + grid_reduce.rs re-documented as CPU reference implementations
  - **Tests:** parity vs CPU mocks at sizes {1, 255, 256, 257, 65k, 1M, 16M} incl. non-power-of-two (skip_if_no_gpu! pattern); f32 sum tolerance ∝ n·ε; histogram exact; compact order-stable; perf targets env-gated (OXIPHYSICS_GPU_BENCH=1), not hard asserts
  - **Risk:** wgpu 29 WGSL reserved words (`pass`) — avoid in identifiers; Blelloch block-sums recursion for >256² elements; naga validation at test time
- [x] GPU radix sort (u32/u64 key, 8-bit digit, 4/8 passes) — (closes marker 6) (planned 2026-06-12) (shipped 2026-06-12; real WGSL histogram+scatter dispatched on Metal, stable per-tile decoupled scatter, GPU exclusive-scan offsets)
  - **Goal:** 10M keys sorted <15 ms; parity vs `parallel_sort.rs` CPU radix. (dep: reduction lib)
  - **Design:** per-pass digit histogram + exclusive scan + scatter, built directly on the reduction lib (the structure `parallel_sort.rs:1022` already simulates on CPU).
  - **Files:** `src/kernels_wgsl/{radix_histogram,radix_scatter}.wgsl` (NEW); new `src/gpu_radix.rs` driver (`radix_sort_u32_gpu`, `radix_sort_pairs_gpu`); `parallel_sort.rs` stays as the parity reference.
  - **Tests:** radix parity vs CPU `radix_sort_u32` at {255,256,257,1024,65536,1_000_000} incl. duplicate-heavy + reverse-sorted; 16M env-gated behind `OXIPHYSICS_GPU_BENCH`; pairs sort keeps payload aligned to keys.
  - **Risk:** wgpu 29 reserved word `pass` (avoid in identifiers); ping-pong key+payload buffers across 4 LSD passes; stable scatter = global_digit_offset[digit] + local_rank_within_digit.
- [x] Radix-sort LBVH build on GPU (Morton + Karras topology) (planned 2026-06-12) (shipped 2026-06-12 as HYBRID: GPU radix sort of (morton,id) pairs → CPU compute_bvh_from_sorted hierarchy; bit-identical to CPU reference, 100k-leaf topology+ray parity on Metal)
  - **Goal:** rebuild 1M-leaf BVH <8 ms; hit-parity vs `bvh/cpu.rs` builder on the existing 10^5-leaf traversal test. (dep: radix sort)
  - **Files:** new `src/gpu_lbvh.rs` (`gpu_lbvh_build`); `src/kernels_wgsl/lbvh.wgsl` (Karras 2012 internal-node range/split); reuses `gpu_radix::radix_sort_pairs_gpu` for (morton,id) sort and `bvh::compute_bvh_from_sorted` as the CPU reference.
  - **Tests:** LBVH 1M-leaf (or 100k to match existing budget) build → topology or ray-hit parity vs CPU `compute_bvh_from_sorted`; empty + single-element edge cases (no panic).
  - **Risk:** Karras longest-common-prefix (δ) with index-augmented tie-break for duplicate Morton codes; fall back to GPU-sort + CPU `compute_bvh_from_sorted` hierarchy if on-GPU topology parity is flaky.
- [ ] On-GPU Karras hierarchy kernel (replace the hybrid's CPU `compute_bvh_from_sorted` step) — (follow-on to the shipped hybrid LBVH)
  - **Goal:** build the BVH internal-node hierarchy fully on the GPU via the Karras 2012 parallel longest-common-prefix (LCP) radix tree, eliminating the CPU hierarchy step; ray-hit parity vs the hybrid on the 100k-leaf scene.
  - **Files:** new `src/kernels_wgsl/lbvh.wgsl` (Karras internal-node range + split with index-augmented δ tie-break for duplicate Morton codes); `src/gpu_lbvh.rs` (on-GPU hierarchy path behind the same `gpu_lbvh_build` entry, hybrid kept as fallback).

### Solver pipelines on GPU
- [x] SPH full pipeline on GPU — (promotes 4 of the 7 marker-13 templates) (planned 2026-06-12) (shipped 2026-06-13: real GPU-resident WCSPH pipeline on Metal — `sph_cell_list`→GPU radix-sort pairs→`histogram`+`exclusive_scan`→`sph_density`→`sph_force`→`sph_integrate`→`sph_boundary`, all dispatched through `WgpuBackendReal`; legacy CPU-shadow `WgpuBackend` stub dropped from `sph_gpu.rs`; cell-list-neighbor parity vs CPU brute force = exact, dam-break density L2 = 2e-6 over 20 steps; 5 new WGSL kernels + 10 `functions.rs` registry templates now naga-validated in `tests/wgsl_validation.rs`)
  - **Goal:** 250k-particle dam break ≥60 steps/s on M-series, 1e-4 L2 parity vs CPU path. — L2 parity DONE (2e-6, exceeds the 1e-4 target); steps/s perf target split out below.
  - **Design:** promote shader_registry templates (CELL_LIST, SPH_FORCE, INTEGRATE, BOUNDARY_ENFORCE) from validated-only to dispatched kernels alongside the already-real density kernel. (dep: reduction lib for cell-list compaction) — DONE: implemented as dedicated `src/kernels_wgsl/sph_*.wgsl` kernels (the `functions.rs` templates were repaired to valid WGSL and are now naga-validated; the dispatched pipeline uses the new `kernels_wgsl` constants).
  - **Files:** `src/sph_gpu.rs` (pipeline orchestration), `src/shaders/functions.rs` (templates repaired to valid WGSL), `src/kernels_wgsl/sph_{cell_list,density,force,integrate,boundary}.wgsl` (NEW dispatched kernels), `src/kernels_wgsl/mod.rs` (constants), `tests/wgpu_kernels.rs` + `tests/wgsl_validation.rs` (NEW).
  - **Tests:** `test_cell_list_correctness` (GPU cell-list neighbor sets vs CPU brute force) + `test_sph_gpu_dam_break_parity` (~200-particle 2-D dam break, 20 steps, density L2 + momentum bound) in `tests/wgpu_kernels.rs`; 7 naga compile checks in `tests/wgsl_validation.rs`.
- [ ] SPH GPU 250k @ ≥60 steps/s perf target — (follow-on to the shipped SPH pipeline; correctness already at 2e-6 L2)
  - **Goal:** 250k-particle dam break ≥60 steps/s on M-series. The sort+scan fold is DONE (F1971/F1972): the spatial-hash rebuild now runs on the one resident backend of `SphGpuState` (`src/sph_gpu/resident.rs`: `sph_cell_list` → multi-pass GPU radix sort with a multi-level GPU scan → `sph_cell_histogram` → GPU bucket scan) with every buffer preallocated at construction, no key/id readback, no per-step device creation and every pipeline compiled once (tests: buffer count constant over 1000 steps, pipeline count fixed after the first step). Remaining: measure the throughput at 250k particles and add an env-gated (`OXIPHYSICS_RTX_BENCH`) steps/s assertion; the per-step upload/readback of the particle state is the next cost to remove (v0.3.0 residency).
- [ ] XPBD cloth kernels — (closes marker 5)
  - **Goal:** 65k-vertex cloth at 120 substeps/s; parity vs `gpu_cloth.rs` CPU mock within 1e-4 max vertex error over 100 steps.
  - **Design:** graph-colored constraint batches, one dispatch per color.
  - **Files:** `src/kernels_wgsl/xpbd_cloth.wgsl` (NEW); `gpu_cloth.rs` becomes the CPU parity reference.
- [ ] GPU contact solver (graph-colored PGS) — (closes marker 4 and the rigid half of 15)
  - **Goal:** 50k contacts × 8 iterations <4 ms; post-solve velocity parity vs `gpu_rigid.rs` sequential impulse within 1e-3.
  - **Design:** greedy coloring on CPU (later GPU), per-color dispatch over `kernels/rigid` contact batches.
  - **Files:** `src/kernels/rigid/` (contact-batch dispatch goes real), `src/kernels_wgsl/pgs_contact.wgsl` (NEW); `gpu_rigid.rs` stays as the parity reference.

### Profiling, cleanup & conformance
- [ ] wgpu timestamp-query profiling integrated with umbrella profiler — verified partially exists (`src/compute/timestamp.rs` GPU timestamp pair + CPU fallback; Phase 4.3 shipped QuerySet path). Remaining: per-kernel scopes exported into `oxiphysics::profiler::FrameReport` / folded stacks.
  - **Goal:** flamegraph shows named GPU kernel spans.
  - **Files:** `src/compute/timestamp.rs` (named scopes), umbrella `profiler.rs` `FrameReport` ingestion (pairs with the umbrella chrome-tracing item).
- [ ] Retire legacy `WgpuBackend` CPU-shadow stub — (closes marker 1)
  - **Goal:** zero literal TODO comments in `src/compute/`; one entry point (`WgpuBackendReal` + explicit `CpuBackend` fallback); the 3 TODO blocks deleted or rewritten as doc pointers.
  - **Files:** `src/compute/wgpu_backend.rs:133,245,271` (the three obsolete TODO sites).
- [ ] Kernel parity conformance suite v1
  - **Goal:** every dispatched WGSL kernel enumerated in a const registry with a seeded parity test vs CPU reference + tolerance table (extends `tests/wgpu_kernels.rs`/`tests/lbm_kernels.rs`); skip-not-fail on headless.
  - **Files:** kernel registry const in `src/shaders/` (next to `shader_registry`), `tests/wgpu_kernels.rs`, `tests/lbm_kernels.rs`.

## v0.3.0 — Residency, accuracy, second backend

Exit criteria: 10k-body rigid scene fully resident (≤2 readbacks/s); DS f64-emulated reductions opt-in; autotuner ≥10% win on at least one adapter family.

- [ ] Full-pipeline GPU residency (broadphase→narrowphase→solver, no per-step readback)
  - **Goal:** 10k-body rigid scene with ≤2 buffer readbacks/second (readback only on user query).
  - **Design:** persistent bind groups, ping-pong state, indirect dispatch for variable pair counts; LBM already proved the resident-stepping pattern (`test_lbm_d3q19_gpu_resident_stepping`).
  - **Files:** `src/compute/wgpu_backend.rs` (indirect dispatch, persistent bind groups); orchestration across the v0.2.0 kernel set.
- [ ] Double-single (DS) f64-emulation WGSL arithmetic (wgpu has no f64)
  - **Goal:** DS dot product of 1M values within 1e-12 rel of CPU f64; opt-in per kernel for constraint RHS and reductions.
  - **Design:** two-f32 Dekker/Knuth add/mul library as a WGSL include prefix.
  - **Files:** `src/kernels_wgsl/ds_math.wgsl` (NEW include prefix) + opt-in flag on dispatch entry points.
- [ ] Async compute + double-buffered transfer overlap — (closes marker 11)
  - **Goal:** ≥20% step-time reduction at 1M particles from upload/compute overlap.
  - **Design:** replace `scheduler.rs` simulated queue with real submissions + the existing frame-graph types.
  - **Files:** `src/scheduler.rs:494–528` (simulated queue replaced with real submissions).
- [ ] Kernel autotuner
  - **Goal:** autotuned workgroup size ({32..256} search, cached per `AdapterInfo` in a JSON cache file) beats fixed-64 by ≥10% on at least one adapter family.
  - **Files:** tuning search in `src/compute/wgpu_backend.rs`; cache keyed on `AdapterInfo` (JSON under the user cache dir — never a hardcoded absolute path).
- [~] CUDA backend activation (feature-gated `cuda-backend`, non-default — cudarc is the approved FFI pattern; default stays pure Rust; needs RTX hardware — env-gated tests) — (closes markers 2, 3)
  - **Goal:** SPH density + scan parity tests pass on RTX hardware (env-gated like `OXIPHYSICS_RTX_BENCH`); the 3 stub-PTX kernels replaced with real compiled kernels; `gpu_bench.rs:129` CUDA stub wired.
  - **Verified:** plumbing exists (device/buffers/launch/unified-memory under cudarc), pending hardware.
  - **Ready when:** RTX-class hardware is available for the env-gated verification runs.

## v1.0 — Production & validation

- [ ] Multi-GPU device groups (dep: residency)
  - **Goal:** 1.6x scaling on 2 adapters at 4M SPH particles via spatial domain decomposition + halo exchange.
  - **Files:** device-group layer in `src/compute/` (NEW); feeds the root Phase 30 10M-particle single-workstation goal.
- [ ] Conformance suite across ALL kernels + honest module taxonomy — (closes remaining markers 8, 9, 10, 12 + residuals of 13/15)
  - **Goal:** zero modules whose docs claim "GPU" while CPU-only; each of the ~38 CPU-mock modules either promoted to a dispatched kernel or re-documented as "CPU reference implementation"; adapter capability matrix doc.
  - **Scope detail:** marker-13 residual templates (LBM_BGK_D2Q9, SDF_COMPUTE) get an explicit dispatched-or-reference verdict; `kernels/md_force` (marker 15) likewise.

## Deferred / research track

- [~] Custom compute fences / `wgpu::CommandEncoder` orchestration across multi-pass pipelines (parked since Phase 4.6) — Ready when: v0.3.0 full-pipeline residency or async-compute items need explicit multi-pass frame graphs beyond single-dispatch encode.
- [~] WebGPU in browser (distinct from desktop wgpu) — Ready when: the v0.2.0 SPH GPU pipeline lands here; implementation is owned by `oxiphysics-wasm` (its v0.3.0 "WebGPU compute in browser" item compiles this crate's wgpu backend to wasm32).

<!-- AUDIT-2026-10-12 BEGIN -->
## Audit 2026-10-12 — verified findings (Phase 32 remediation)

Verified: 195 (critical 6, high 44, medium 98, low 47); 36 further raw findings pending second-round verification. IDs refer to the audit ledger; `[x]` = fixed.

### Critical (6)
- [ ] **F1794** (critical, fabrication, medium/medium) `fluid_gpu.rs:499` — LBMGpuKernels: fake velocities, no-op streaming, identity bounce-back
- [ ] **F1799** (critical, fabrication, medium/medium) `fluid_sim_gpu/functions.rs:402` — FLIP G2P: x-only transfer, zero FLIP delta, advects by dx
- [ ] **F1831** (critical, bug, easy/trivial) `gpu_md_solver.rs:167` — LJ pair force sign inverted: repulsion becomes attraction
- [x] **F1882** (critical, fabrication, hard/medium) `lbm_gpu.rs:192` — LBM GPU kernel ignores body force, solid walls and lid velocity
- [x] **F1883** (critical, bug, hard/large) `lbm_gpu.rs:376` — Default GPU LBM step drops walls, body force and lid
- [x] **F1885** (critical, bug, hard/medium) `lbm_gpu.rs:423` — LbmSimulation GPU path ignores walls, Guo forcing and lid velocity

### High (44)
- [ ] **F1752** (high, bug, easy/small) `cell_list/types.rs:672` — CellList::update_incremental overwrites next cell's first atom
- [ ] **F1753** (high, bug, medium/small) `collision_gpu/functions.rs:116` — GJK never tests tetrahedron containment; distances are wrong
- [ ] **F1755** (high, bug, easy/small) `collision_gpu/types.rs:547` — Persistent manifolds never expire and accumulate duplicate points
- [ ] **F1756** (high, bug, easy/trivial) `collision_gpu/types.rs:637` — Contact cache zips all pairs with filtered contacts: wrong keys
- [ ] **F1759** (high, fabrication, medium/small) `compute/cpubackend_traits.rs:38` — CpuBackend::dispatch discards all kernel work
- [ ] **F1768** (high, bug, medium/medium) `compute/timestamp.rs:123` — 1-D dispatch exceeds 65535 workgroup limit for large N
- [ ] **F1778** (high, api, medium/medium) `compute/wgpu_backend.rs:706` — read_buffer_* return empty Vec on failure; callers return garbage
- [ ] **F1781** (high, panic, medium/small) `compute/wgpu_backend.rs:931` — No wgpu error scopes: validation errors panic instead of returning Err
- [ ] **F1789** (high, fabrication, medium/medium) `deformable_gpu/types.rs:580` — FEMGpuKernel forces are three edge springs, not tetrahedral FEM
- [ ] **F1790** (high, fabrication, hard/medium) `deformable_gpu/types.rs:1051` — CorotationalElement forces are edge springs; nonzero at rest pose
- [ ] **F1798** (high, bug, easy/trivial) `fluid_sim_gpu/functions.rs:86` — SPH pressure term is acceleration but integrator divides by density
- [ ] **F1803** (high, fabrication, medium/small) `fluid_sim_gpu/types.rs:641` — FLIP pressure solve is a no-op: fluid flags never set
- [ ] **F1804** (high, bug, medium/small) `flux_compute.rs:682` — Roe Euler flux wave strengths ignore energy jump and gamma
- [ ] **F1805** (high, bug, easy/trivial) `flux_compute.rs:872` — Characteristic reconstruction swaps w+/w- bindings, flipping u1 sign
- [ ] **F1807** (high, bug, medium/small) `gpu_bench.rs:159` — CPU-only paths silently run GPU; LbmGpuSolver::new_cpu reads stale state
- [ ] **F1823** (high, bug, easy/small) `gpu_fluid_euler.rs:255` — Euler pressure projection removes only dt times the divergence
- [ ] **F1838** (high, fabrication, medium/medium) `gpu_nn.rs:234` — BackpropGpu fabricates gradients; Conv1D silently runs as Dense
- [ ] **F1840** (high, perf, easy/trivial) `gpu_particles.rs:482` — ParticlePool::recycle_dead is O(capacity × free) on every tick
- [x] **F1841** (high, perf, medium/medium) `gpu_primitives.rs:21` — Each primitive call creates a fresh wgpu Instance/Adapter/Device
- [x] **F1843** (high, perf, medium/medium) `gpu_radix.rs:165` — Radix sort/scan create a fresh wgpu device per call and pass
- [x] **F1849** (high, bug, easy/trivial) `gpu_sdf.rs:337` — sdf_smooth_union operands swapped; returns max not min
- [x] **F1852** (high, bug, easy/small) `gpu_sparse_solver.rs:61` — add_entry/finalize builds corrupt CSR for unordered or repeated input
- [ ] **F1854** (high, bug, easy/trivial) `gpu_sph_density.rs:199` — CPU artificial viscosity sign reversed vs WGSL; approaching particles attract
- [x] **F1857** (high, fabrication, hard/medium) `gpu_sph_pressure.rs:251` — pcisph_gpu_correction never predicts positions or corrects pressure
- [x] **F1858** (high, bug, easy/trivial) `gpu_sph_solver.rs:86` — Wendland C2 uses support-2h form truncated at h
- [ ] **F1859** (high, bug, easy/trivial) `gpu_sph_solver.rs:94` — Wendland kernel uses the support-2h form but truncates at r=h
- [ ] **F1860** (high, bug, medium/medium) `gpu_sph_solver.rs:131` — gpu_sph_step: unit mass, fixed rho0, 5.5x-off kernel, anti-viscosity
- [x] **F1861** (high, bug, easy/small) `gpu_sph_solver.rs:243` — SPH viscosity uses negative dW/dr as Laplacian: anti-diffusive
- [x] **F1874** (high, bug, easy/small) `kernels/md_force/functions.rs:321` — Harmonic angle forces have reversed (anti-restoring) sign
- [x] **F1875** (high, bug, medium/medium) `kernels/rigid/functions.rs:65` — Rigid integrators use inverse mass as inverse inertia
- [ ] **F1876** (high, bug, medium/small) `kernels/rigid/types.rs:152` — ContactBatchProcessor applies inverse mass twice to contact impulses
- [x] **F1884** (high, perf, easy/small) `lbm_gpu.rs:399` — LbmSimulation retries GPU init and eprintln!s on every step
- [ ] **F1888** (high, bug, medium/small) `lbm_gpu.rs:500` — CPU LBM Guo forcing lacks source term and injects ωF/2 per step
- [x] **F1890** (high, fabrication, medium/small) `lbm_gpu.rs:577` — Lid velocity BC writes only solid cells; set_lid_velocity no-op
- [x] **F1896** (high, bug, medium/medium) `lbm_gpu.rs:1117` — LBM lid and walls ignored on default GPU path; CPU lid dead
- [ ] **F1908** (high, bug, easy/small) `parallel/functions.rs:273` — Parallel scans write through as_ptr() of immutable Vec (UB)
- [ ] **F1922** (high, fabrication, medium/medium) `pipeline/functions.rs:8` — PhysicsPipeline reports fabricated pair/constraint counts; narrowphase no-op
- [ ] **F1930** (high, bug, easy/trivial) `raytracing/functions.rs:352` — path_trace_step uses wrong 2cosθ weight; panics on empty materials
- [x] **F1941** (high, fabrication, medium/medium) `sdf_compute/functions.rs:10` — fast_sweeping_update is a one-direction L1 chamfer, not Eikonal
- [x] **F1945** (high, fabrication, medium/medium) `sdf_compute/functions.rs:487` — marching_cubes fans edge vertices without a triangle table
- [x] **F1947** (high, fabrication, medium/medium) `sdf_compute/functions.rs:625` — sdf_open/sdf_close are exact identity operations
- [x] **F1948** (high, bug, easy/small) `sdf_compute/functions.rs:664` — Laplacian smooth runs backward heat equation; sharpen actually smooths
- [x] **F1971** (high, bug, medium/medium) `sph_gpu.rs:711` — GPU buffers never freed; SPH and BVH leak VRAM every call
- [x] **F1972** (high, perf, medium/medium) `sph_gpu.rs:747` — Each SPH GPU step creates three wgpu devices and recompiles pipelines

### Medium (98)
- [ ] **F1743** (medium, api, easy/small) `crates/oxiphysics-gpu/Cargo.toml:19` — cpu-only and headless-ci features are no-ops despite README claim
- [ ] **F1748** (medium, perf, easy/small) `bvh/gpu.rs:264` — CPU BVH traversal lacks closest-hit culling; serial over rays
- [ ] **F1749** (medium, bug, medium/small) `cell_list/types.rs:246` — build_parallel hardcodes cell_size=1 and origin 0
- [ ] **F1750** (medium, bug, easy/small) `cell_list/types.rs:420` — for_each_pair silently misses pairs when cutoff > cell_size
- [ ] **F1751** (medium, bug, easy/small) `cell_list/types.rs:429` — GpuCellList::for_each_pair misses pairs when cutoff > cell_size
- [ ] **F1761** (medium, stub, hard/medium) `compute/cuda_backend.rs:183` — PTX_* kernel constants are header-only stubs advertised as shipped
- [ ] **F1762** (medium, bug, medium/small) `compute/cuda_backend.rs:579` — CUDA handle index desyncs from real buffer after alloc failure
- [ ] **F1763** (medium, bug, medium/small) `compute/cuda_backend.rs:585` — CUDA buffer handles desync from device buffers after a failed alloc
- [ ] **F1764** (medium, api, medium/small) `compute/cuda_backend.rs:830` — CUDA launch/write/read silently swallow errors; raw-pointer split
- [ ] **F1765** (medium, panic, easy/trivial) `compute/cuda_backend.rs:881` — query_device_info lacks catch_unwind; panics without CUDA driver
- [ ] **F1771** (medium, api, easy/small) `compute/types.rs:753` — compute::GpuError lacks std::error::Error; clashes with crate::GpuError
- [ ] **F1773** (medium, test, easy/small) `compute/wgpu_backend.rs:274` — WGSL_PARALLEL_SCAN uses reserved words; shipped shaders not naga-tested
- [ ] **F1779** (medium, bug, easy/small) `compute/wgpu_backend.rs:740` — Readback failures return empty Vec; consumers fabricate zeros
- [ ] **F1780** (medium, panic, easy/medium) `compute/wgpu_backend.rs:764` — Buffer writes lack size/handle checks; wgpu validation panics
- [ ] **F1782** (medium, api, medium/small) `compute/wgpu_backend.rs:962` — dispatch_wgsl: invalid WGSL or bindings panic; binding type ignored
- [ ] **F1783** (medium, perf, medium/medium) `compute/wgpu_backend.rs:993` — Every dispatch/readback is a blocking round-trip with fresh allocations
- [ ] **F1785** (medium, panic, easy/trivial) `compute_pipeline.rs:366` — PipelineCache::new(0) then insert panics; 'LRU' is FIFO
- [ ] **F1791** (medium, bug, medium/medium) `deformable_gpu/types.rs:1176` — XPBD volume uses compliance as stiffness; shape match no rotation
- [ ] **F1792** (medium, fabrication, medium/medium) `deformable_gpu/types.rs:1212` — XPBD shape matching lacks rotation; volume compliance acts as stiffness
- [ ] **F1793** (medium, api, medium/medium) `error.rs:9` — Crate Error enum is unused, not non_exhaustive, no From impls
- [ ] **F1795** (medium, perf, medium/small) `fluid_gpu.rs:663` — SPHGpuKernels density, pressure and viscosity loops are O(n^2) serial
- [ ] **F1796** (medium, bug, medium/small) `fluid_gpu.rs:756` — SPH viscosity uses W, not Laplacian; NS viscosity unused
- [ ] **F1797** (medium, perf, medium/small) `fluid_sim_gpu/functions.rs:34` — SPH density/forces are O(n^2); 'parallel' variant is a serial copy
- [ ] **F1802** (medium, bug, medium/small) `fluid_sim_gpu/types.rs:424` — LbmD2Q9 inlet/outlet applied once, then overwritten by periodic streaming
- [ ] **F1806** (medium, bug, easy/trivial) `gpu_bench.rs:125` — Benchmark wgpu detection uses always-failing stub
- [ ] **F1809** (medium, bug, easy/small) `gpu_bench.rs:538` — Benchmarks gate wgpu on the stub and mislabel GPU time as CPU
- [ ] **F1811** (medium, bug, medium/small) `gpu_cloth.rs:520` — Cloth 'XPBD' lacks λ accumulation; gravity hardcoded
- [ ] **F1812** (medium, fabrication, easy/small) `gpu_collision_detection.rs:292` — batch_gjk_dispatch is a bounding-sphere heuristic, not GJK
- [ ] **F1813** (medium, fabrication, easy/small) `gpu_collision_ext.rs:287` — gpu_gjk_distance uses vertex-pair minimum, wrong for convex hulls
- [ ] **F1818** (medium, bug, medium/small) `gpu_fluid.rs:118` — MacGrid projection drops ρ/Δt; viscosity/density unused
- [ ] **F1819** (medium, fabrication, medium/small) `gpu_fluid.rs:178` — FluidSimGpu ignores viscosity/density; projection removes only dt*div
- [ ] **F1820** (medium, test, easy/medium) `gpu_fluid.rs:529` — About 300 tests contain no assertion; several discard the computed result
- [ ] **F1821** (medium, bug, medium/small) `gpu_fluid_euler.rs:125` — vorticity_at stencil wrong; confinement drops fz, no dt
- [ ] **F1822** (medium, bug, medium/small) `gpu_fluid_euler.rs:190` — Euler projection Poisson RHS lacks 1/dt; per-sweep p clone
- [ ] **F1824** (medium, bug, medium/small) `gpu_fluid_euler.rs:342` — w-velocity advection ignores y-velocity; back-trace pinned to row
- [ ] **F1825** (medium, bug, easy/small) `gpu_fluid_euler.rs:400` — Vorticity confinement drops fz, omits dt*dx, underflows for nz==0
- [ ] **F1826** (medium, bug, medium/small) `gpu_fluid_euler.rs:496` — w-velocity advection ignores y motion (fixed j row)
- [ ] **F1828** (medium, bug, medium/small) `gpu_lbm.rs:291` — NoSlip walls BGK-collided after bounce-back; 'ZouHe' is equilibrium
- [ ] **F1830** (medium, perf, medium/small) `gpu_md_solver.rs:148` — MD forces and energy are O(n^2) all-pairs with sqrt before cutoff
- [ ] **F1832** (medium, fabrication, medium/small) `gpu_mesh_processing.rs:208` — gpu_loop_subdivision is midpoint subdivision, not Loop
- [ ] **F1833** (medium, perf, medium/medium) `gpu_mesh_processing.rs:326` — Mesh decimation is O(T²): full rescan and remap per collapse
- [ ] **F1834** (medium, perf, easy/small) `gpu_mesh_processing.rs:423` — gpu_weld_vertices is O(V^2) with sqrt per pair
- [ ] **F1835** (medium, missing-feature, hard/medium) `gpu_neural_solver.rs:164` — PINN solver has no training, derivatives, or symmetry-breaking initialisation
- [ ] **F1839** (medium, bug, medium/small) `gpu_particles.rs:275` — Turbulence curl noise differentiates uncorrelated hash white noise
- [ ] **F1847** (medium, bug, easy/small) `gpu_rigid.rs:185` — apply_impulse uses body-space inverse inertia on world-space torque
- [ ] **F1848** (medium, stub, medium/medium) `gpu_rigid.rs:366` — Sequential-impulse solver is linear-only with hardcoded restitution
- [ ] **F1851** (medium, bug, easy/small) `gpu_sort.rs:253` — counting_sort silently clamps values above max_val and can over-allocate
- [ ] **F1855** (medium, perf, medium/medium) `gpu_sph_density.rs:216` — gpu_neighbor_list documented as cell list but is O(N²) brute force
- [ ] **F1863** (medium, fabrication, easy/small) `gpu_thermal.rs:231` — Neumann thermal boundary condition is a silent no-op
- [ ] **F1864** (medium, bug, easy/small) `gpu_voxel.rs:107` — Voxelizer marks voxel 0 for triangles at negative coordinates
- [ ] **F1869** (medium, perf, easy/small) `kernels/broadphase.rs:223` — UniformGridGpu ignores dims: unbounded cell loop and O(P*N) find
- [ ] **F1870** (medium, perf, easy/small) `kernels/broadphase.rs:255` — UniformGridGpu::query_pairs linear-searches AABBs for every candidate pair
- [ ] **F1873** (medium, stub, hard/medium) `kernels/md_force/functions.rs:215` — pppm_mesh_energy_estimate returns Σq²/N, no Green's function
- [ ] **F1877** (medium, perf, medium/medium) `kernels/sph.rs:163` — Public CPU SPH kernels are all-pairs O(N²) despite in-file NeighborList
- [ ] **F1878** (medium, bug, easy/small) `kernels/sph.rs:357` — NeighborList query unclamped while build clamps: particles lose neighbours
- [ ] **F1880** (medium, fabrication, medium/medium) `kernels/sph.rs:757` — surface_tension_force claims CSF curvature but never computes it
- [ ] **F1881** (medium, bug, easy/small) `kernels/sph.rs:997` — Pressure 'force' returns acceleration but wcsph_euler_step divides by mass
- [ ] **F1886** (medium, panic, easy/small) `lbm_gpu.rs:441` — LBM step expect() panics; GPU init is retried and logged every step
- [ ] **F1887** (medium, bug, easy/small) `lbm_gpu.rs:463` — Failed GPU readback silently zeroes all LBM populations
- [ ] **F1892** (medium, bug, medium/small) `lbm_gpu.rs:771` — LbmGpuSolver GPU/CPU paths differ; CPU variant reads stale data
- [ ] **F1893** (medium, api, easy/small) `lbm_gpu.rs:780` — LBM/SPH configs never validated; invalid omega silently replaced
- [ ] **F1894** (medium, bug, easy/small) `lbm_gpu.rs:829` — LbmGpuSolver::new_cpu and CPU fallback still run the GPU path
- [ ] **F1900** (medium, bug, easy/small) `memory.rs:233` — Buffer pool read/write trust stale caller handle offsets
- [ ] **F1901** (medium, fabrication, easy/trivial) `memory.rs:794` — TransferQueue::copy_from_gpu returns zeros instead of pool data
- [ ] **F1903** (medium, bug, easy/trivial) `neural_compute/types.rs:899` — GELU uses sqrt(2/sqrt(pi)) instead of sqrt(2/pi)
- [ ] **F1910** (medium, panic, easy/trivial) `parallel_sort.rs:341` — nth_element recursion depth O(n) on duplicate keys: stack overflow
- [ ] **F1911** (medium, panic, easy/trivial) `parallel_sort.rs:371` — nth_element degrades to O(n²) with n-deep recursion on equal keys
- [ ] **F1912** (medium, bug, easy/small) `parallel_sort.rs:1167` — histogram_bucket_sort skips sorting when value range < f64::EPSILON
- [ ] **F1915** (medium, bug, easy/small) `particle_system/types.rs:69` — GridParticleCollision nondeterministic (HashMap order), misses pairs
- [ ] **F1918** (medium, bug, easy/trivial) `particle_system/types.rs:794` — retire_expired condition is contradictory; always returns 0
- [ ] **F1919** (medium, stub, easy/small) `particle_system/types.rs:957` — Cone emitter discards half_angle; behaves as Point
- [ ] **F1920** (medium, fabrication, medium/small) `path_tracer.rs:478` — PathTracerScene lights never contribute; material index panics
- [ ] **F1921** (medium, bug, easy/small) `path_tracer.rs:707` — render_pass NaN/inf for width/height 1 or samples_per_pass 0; serial
- [ ] **F1923** (medium, bug, easy/small) `pipeline/types.rs:664` — GpuMemoryPool::free accepts double/overlapping frees, re-hands memory
- [ ] **F1924** (medium, bug, easy/small) `ray_marching.rs:410` — SoftShadow::compute fully shadows points on the surface
- [ ] **F1932** (medium, bug, easy/trivial) `raytracing/types.rs:132` — PathState::russian_roulette is deterministic and biased
- [ ] **F1934** (medium, panic, medium/medium) `raytracing/types.rs:659` — Bvh::build panics on NaN centroid; sweep SAH re-sorts 4x per node
- [ ] **F1936** (medium, fabrication, medium/medium) `raytracing/types.rs:1016` — RenderConfig spp/AO/soft-shadow/DoF/max_depth never read
- [ ] **F1938** (medium, panic, easy/trivial) `scheduler.rs:254` — TaskGraph::critical_path panics on NaN estimated_ms
- [ ] **F1939** (medium, bug, easy/small) `scheduler.rs:345` — TaskScheduler barriers stored but never enforced
- [ ] **F1940** (medium, bug, easy/small) `scheduler.rs:362` — TaskScheduler ignores barriers registered via add_barrier
- [ ] **F1942** (medium, bug, medium/medium) `sdf_compute/functions.rs:227` — mesh_to_sdf yields unsigned distance; unchecked index panic
- [ ] **F1943** (medium, bug, medium/medium) `sdf_compute/functions.rs:251` — mesh_to_sdf writes unsigned distances into a signed grid
- [ ] **F1944** (medium, bug, easy/small) `sdf_compute/functions.rs:367` — MC_EDGE_TABLE entries 213-215 wrong; no triangle table
- [ ] **F1949** (medium, bug, easy/small) `sdf_compute/functions_2.rs:387` — generate_sdf_grid ignores extent[1] and extent[2]
- [ ] **F1950** (medium, stub, medium/small) `sdf_compute/functions_2.rs:411` — march_surface emits degenerate placeholder triangles
- [ ] **F1951** (medium, fabrication, medium/medium) `sdf_compute/functions_2.rs:462` — SDF surface extractors emit degenerate or mis-triangulated meshes
- [ ] **F1952** (medium, bug, easy/small) `sdf_compute/types.rs:426` — SdfGrid::gradient_at_point substitutes 0.0 for out-of-grid samples
- [ ] **F1953** (medium, panic, easy/trivial) `shader_registry/types.rs:122` — VariantProfileRegistry::resolve recurses forever on cyclic bases
- [ ] **F1954** (medium, bug, easy/small) `shader_registry/types.rs:408` — CompiledVariant.workgroup_size ignores substituted WGSL
- [ ] **F1957** (medium, fabrication, medium/medium) `shaders/functions.rs:462` — SpirVModule::from_wgsl returns fabricated SPIR-V bytes
- [ ] **F1959** (medium, bug, easy/small) `shaders/types.rs:473` — ShaderCompilationPipeline cache keyed by name only
- [ ] **F1960** (medium, panic, easy/small) `sparse_gpu.rs:96` — SparseTriplet::to_csr panics on out-of-range row; add asserts
- [ ] **F1962** (medium, bug, easy/medium) `sparse_gpu.rs:464` — CG/BiCGSTAB/PCG: no validation/breakdown checks, Vec alloc per iteration
- [ ] **F1963** (medium, bug, easy/small) `sparse_gpu.rs:516` — BiCGSTAB omega=0 yields NaN; breakdowns look like convergence
- [ ] **F1967** (medium, bug, easy/small) `sph_gpu.rs:322` — Default SPH particle mass rho0*(2h)^3 gives ~8x rest density
- [ ] **F1968** (medium, perf, medium/medium) `sph_gpu.rs:403` — SPH CPU fallback is serial O(N^2) all-pairs
- [ ] **F1969** (medium, perf, medium/medium) `sph_gpu.rs:409` — CPU SPH fallback is single-threaded O(N²)

### Low (47) — tracked in the external ledger
F1744, F1745, F1758, F1760, F1766, F1767, F1769, F1770, F1772, F1774, F1775, F1776, F1777, F1787, F1788, F1808, F1814, F1817, F1827, F1845, F1846, F1850, F1853, F1856, F1866, F1867, F1871, F1872, F1879, F1889, F1891, F1897, F1898, F1899, F1902, F1905, F1914, F1916, F1917, F1925, F1927, F1928, F1937, F1955, F1956, F1970, F1973

<!-- AUDIT-2026-10-12 END -->
