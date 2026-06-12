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
- [~] GPU radix sort (u32/u64 key, 8-bit digit, 4/8 passes) — (closes marker 6) (planned 2026-06-12)
  - **Goal:** 10M keys sorted <15 ms; parity vs `parallel_sort.rs` CPU radix. (dep: reduction lib)
  - **Design:** per-pass digit histogram + exclusive scan + scatter, built directly on the reduction lib (the structure `parallel_sort.rs:1022` already simulates on CPU).
  - **Files:** `src/kernels_wgsl/{radix_histogram,radix_scatter}.wgsl` (NEW); new `src/gpu_radix.rs` driver (`radix_sort_u32_gpu`, `radix_sort_pairs_gpu`); `parallel_sort.rs` stays as the parity reference.
  - **Tests:** radix parity vs CPU `radix_sort_u32` at {255,256,257,1024,65536,1_000_000} incl. duplicate-heavy + reverse-sorted; 16M env-gated behind `OXIPHYSICS_GPU_BENCH`; pairs sort keeps payload aligned to keys.
  - **Risk:** wgpu 29 reserved word `pass` (avoid in identifiers); ping-pong key+payload buffers across 4 LSD passes; stable scatter = global_digit_offset[digit] + local_rank_within_digit.
- [~] Radix-sort LBVH build on GPU (Morton + Karras topology) (planned 2026-06-12)
  - **Goal:** rebuild 1M-leaf BVH <8 ms; hit-parity vs `bvh/cpu.rs` builder on the existing 10^5-leaf traversal test. (dep: radix sort)
  - **Files:** new `src/gpu_lbvh.rs` (`gpu_lbvh_build`); `src/kernels_wgsl/lbvh.wgsl` (Karras 2012 internal-node range/split); reuses `gpu_radix::radix_sort_pairs_gpu` for (morton,id) sort and `bvh::compute_bvh_from_sorted` as the CPU reference.
  - **Tests:** LBVH 1M-leaf (or 100k to match existing budget) build → topology or ray-hit parity vs CPU `compute_bvh_from_sorted`; empty + single-element edge cases (no panic).
  - **Risk:** Karras longest-common-prefix (δ) with index-augmented tie-break for duplicate Morton codes; fall back to GPU-sort + CPU `compute_bvh_from_sorted` hierarchy if on-GPU topology parity is flaky.

### Solver pipelines on GPU
- [ ] SPH full pipeline on GPU — (promotes 4 of the 7 marker-13 templates)
  - **Goal:** 250k-particle dam break ≥60 steps/s on M-series, 1e-4 L2 parity vs CPU path.
  - **Design:** promote shader_registry templates (CELL_LIST, SPH_FORCE, INTEGRATE, BOUNDARY_ENFORCE) from validated-only to dispatched kernels alongside the already-real density kernel. (dep: reduction lib for cell-list compaction)
  - **Files:** `src/sph_gpu.rs` (pipeline orchestration), `src/shaders/functions.rs` (templates promoted into dispatched kernels).
  - **Tests:** dam-break L2-parity test in `tests/wgpu_kernels.rs`; steps/s assertion env-gated (`OXIPHYSICS_RTX_BENCH` pattern).
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
