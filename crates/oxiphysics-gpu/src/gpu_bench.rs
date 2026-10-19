// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Performance benchmarks comparing CPU, wgpu, and CUDA backends.
//!
//! This module measures throughput for the core compute kernels (SPH density,
//! LBM collision, parallel scan) across available backends and reports
//! wall-clock timing and effective GFLOP/s estimates.
//!
//! # Quick usage
//!
//! ```
//! use oxiphysics_gpu::gpu_bench::{GpuBenchHarness, BackendKind};
//!
//! let mut h = GpuBenchHarness::new();
//!
//! // Benchmark SPH density summation for 256 particles
//! let reports = h.bench_sph_density(256);
//! for r in &reports {
//!     println!("{}", r);
//! }
//!
//! // Compare available backends
//! let available = GpuBenchHarness::available_backends();
//! assert!(available.contains(&BackendKind::Cpu));
//! ```
//!
//! # Backend labels are measured, not assumed
//!
//! A [`BackendKind::Cpu`] report always comes from a simulation built with
//! [`ComputePreference::CpuOnly`], which never creates a wgpu device, so the
//! timing really is the CPU reference even on a host with a GPU.  A
//! [`BackendKind::Wgpu`] report is produced only from a simulation built with
//! [`ComputePreference::GpuRequired`] whose every timed step ran on the GPU; on
//! hosts without a usable adapter the report is simply absent.

use std::time::{Duration, Instant};

use crate::GpuError;
use crate::compute_preference::ComputePreference;
use crate::lbm_gpu::{LbmConfig, LbmSimulation};
use crate::sph_gpu::{SphConfig, SphSimulation};

// ── BackendKind ───────────────────────────────────────────────────────────────

/// Identifies a compute backend for benchmark reporting.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum BackendKind {
    /// CPU (Rayon-parallel fallback).
    Cpu,
    /// wgpu (WebGPU compute shaders).
    Wgpu,
    /// CUDA via cudarc.
    Cuda,
}

impl std::fmt::Display for BackendKind {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Cpu => write!(f, "CPU"),
            Self::Wgpu => write!(f, "wgpu"),
            Self::Cuda => write!(f, "CUDA"),
        }
    }
}

// ── GpuBenchReport ────────────────────────────────────────────────────────────

/// Result of a single GPU benchmark run.
#[derive(Debug, Clone)]
pub struct GpuBenchReport {
    /// Kernel / benchmark name.
    pub name: String,
    /// Which backend was measured.
    pub backend: BackendKind,
    /// Problem size (particles, cells, …).
    pub n: usize,
    /// Number of timed iterations.
    pub iterations: u32,
    /// Total wall-clock time.
    pub total: Duration,
    /// Mean time per iteration.
    pub mean: Duration,
    /// Estimated throughput (MFLOP/s or Mparticles/s depending on kernel).
    pub mflops: Option<f64>,
}

impl std::fmt::Display for GpuBenchReport {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "[{:<5} {:>20}] n={:>6} mean={:.3}µs",
            self.backend,
            self.name,
            self.n,
            self.mean.as_secs_f64() * 1e6
        )?;
        if let Some(mf) = self.mflops {
            write!(f, "  {:.1} MFLOPs", mf)?;
        }
        Ok(())
    }
}

// ── Benchmark problem builders ────────────────────────────────────────────────

/// SPH benchmark configuration for `n` particles in the `[-10, 10]³` box
/// (no gravity: a pure density/pressure workload).
fn sph_bench_config(n: usize) -> SphConfig {
    SphConfig {
        n_particles: n,
        smoothing_h: 0.1,
        rest_density: 1000.0,
        gravity: 0.0,
        domain_min: [-10.; 3],
        domain_max: [10.; 3],
        ..SphConfig::default()
    }
}

/// Scatter the `n` particles of `sim` on the regular benchmark grid.
fn seed_sph_grid(sim: &mut SphSimulation, n: usize) {
    let side = (n as f64).cbrt().ceil() as usize + 1;
    for idx in 0..n {
        sim.state.pos_x[idx] = (idx % side) as f64 * 0.1 - 5.0;
        sim.state.pos_y[idx] = ((idx / side) % side) as f64 * 0.1;
        sim.state.pos_z[idx] = (idx / (side * side)) as f64 * 0.1;
    }
}

/// SPH benchmark simulation pinned to the CPU reference.
///
/// The simulation is built with [`ComputePreference::CpuOnly`], so no GPU
/// adapter is requested and every step really runs on the CPU — the reports
/// built from it may honestly carry [`BackendKind::Cpu`].
fn build_cpu_sph(cfg: &SphConfig, n: usize) -> SphSimulation {
    let mut sim = SphSimulation::new(SphConfig {
        compute: ComputePreference::CpuOnly,
        ..cfg.clone()
    });
    seed_sph_grid(&mut sim, n);
    sim
}

/// SPH benchmark simulation that must run on the GPU, or `None` when no GPU
/// path is available (no adapter, wgpu feature off, unsupported size).
fn build_gpu_sph(cfg: &SphConfig, n: usize) -> Option<SphSimulation> {
    let mut sim = SphSimulation::try_new(SphConfig {
        compute: ComputePreference::GpuRequired,
        ..cfg.clone()
    })
    .ok()?;
    seed_sph_grid(&mut sim, n);
    Some(sim)
}

/// LBM benchmark configuration pinned to the CPU reference.
fn lbm_bench_config(nx: usize, ny: usize, nz: usize) -> LbmConfig {
    LbmConfig {
        nx,
        ny,
        nz,
        tau: 0.6,
        compute: ComputePreference::CpuOnly,
        ..LbmConfig::default()
    }
}

/// Time `iterations` SPH steps after `warmup` untimed ones.
///
/// Fails with the step error when a [`ComputePreference::GpuRequired`]
/// simulation cannot run a step on the GPU, so such a run is never reported.
fn time_sph_steps(
    sim: &mut SphSimulation,
    warmup: u32,
    iterations: u32,
) -> Result<Duration, GpuError> {
    let dt = 1.0 / 60.0;
    for _ in 0..warmup {
        sim.try_step(dt)?;
    }
    let t0 = Instant::now();
    for _ in 0..iterations {
        sim.try_step(dt)?;
    }
    Ok(t0.elapsed())
}

// ── GpuBenchHarness ───────────────────────────────────────────────────────────

/// Timing harness for GPU/CPU backend comparison benchmarks.
pub struct GpuBenchHarness {
    /// Warm-up iterations (not timed).
    pub warmup: u32,
    /// Timed iterations.
    pub iterations: u32,
    /// Collected reports.
    pub reports: Vec<GpuBenchReport>,
}

impl GpuBenchHarness {
    /// Create a harness with 2 warm-up and 5 timed iterations.
    pub fn new() -> Self {
        Self {
            warmup: 2,
            iterations: 5,
            reports: Vec::new(),
        }
    }

    /// Return which backends are available in this build.
    ///
    /// CPU is always available.  wgpu and CUDA depend on feature flags and
    /// device availability (they appear in the list only when initialisation
    /// succeeds).
    pub fn available_backends() -> Vec<BackendKind> {
        let mut out = vec![BackendKind::Cpu];

        // Try wgpu — succeeds when a real adapter and device can be created.
        // (The `WgpuBackend` CPU-shadow type never owns a device and always
        // reports itself unavailable, so it cannot be used for detection.)
        #[cfg(feature = "wgpu-backend")]
        {
            if crate::compute::wgpu_backend::real::WgpuBackendReal::try_new().is_ok() {
                out.push(BackendKind::Wgpu);
            }
        }

        // Try CUDA — this is a stub that always reports unavailable.
        if crate::compute::cuda_backend::CudaBackend::try_new(0).is_ok() {
            out.push(BackendKind::Cuda);
        }

        out
    }

    // ── SPH density benchmark ─────────────────────────────────────────────────

    /// Benchmark SPH density summation for `n` particles on all available backends.
    ///
    /// Each particle's density is recomputed from scratch each call to avoid
    /// caching effects.  FLOPs estimated as 10 × N² (distance + kernel eval).
    pub fn bench_sph_density(&mut self, n: usize) -> Vec<GpuBenchReport> {
        let cfg = sph_bench_config(n);
        let mut out = Vec::new();
        // FLOPs: 10 per particle pair (distance + kernel evaluation).
        let flops = 10.0 * n as f64 * n as f64;

        // ── CPU path ──────────────────────────────────────────────────────────
        // Pinned to the CPU reference: no GPU backend is ever initialised.
        let mut cpu_sim = build_cpu_sph(&cfg, n);
        if let Ok(total) = time_sph_steps(&mut cpu_sim, self.warmup, self.iterations) {
            let mflops = flops / (total.as_secs_f64() / self.iterations as f64) / 1e6;
            let r = GpuBenchReport {
                name: "sph_density".to_string(),
                backend: BackendKind::Cpu,
                n,
                iterations: self.iterations,
                total,
                mean: total / self.iterations,
                mflops: Some(mflops),
            };
            out.push(r.clone());
            self.reports.push(r);
        }

        // ── wgpu path (only when every step really ran on the GPU) ────────────
        if let Some(mut gpu_sim) = build_gpu_sph(&cfg, n)
            && let Ok(total) = time_sph_steps(&mut gpu_sim, self.warmup, self.iterations)
        {
            let mflops = flops / (total.as_secs_f64() / self.iterations as f64) / 1e6;
            let r = GpuBenchReport {
                name: "sph_density".to_string(),
                backend: BackendKind::Wgpu,
                n,
                iterations: self.iterations,
                total,
                mean: total / self.iterations,
                mflops: Some(mflops),
            };
            out.push(r.clone());
            self.reports.push(r);
        }

        out
    }

    // ── LBM benchmark ─────────────────────────────────────────────────────────

    /// Benchmark one LBM BGK step on an `nx × ny × nz` domain.
    ///
    /// FLOPs estimated as 120 × nc (19 distribution reads + BGK + streaming).
    pub fn bench_lbm_step(&mut self, nx: usize, ny: usize, nz: usize) -> Vec<GpuBenchReport> {
        let cfg = lbm_bench_config(nx, ny, nz);
        let nc = nx * ny * nz;
        let mut out = Vec::new();

        // CPU path
        {
            // `lbm_bench_config` pins the CPU reference: this report is labelled
            // `BackendKind::Cpu`, so the default wgpu path must not run on the
            // GPU behind its back.
            let mut sim = LbmSimulation::new(cfg.clone());
            sim.set_lid_velocity(0.1, 0.0, 0.0);

            for _ in 0..self.warmup {
                sim.step();
            }
            let t0 = Instant::now();
            for _ in 0..self.iterations {
                sim.step();
            }
            let total = t0.elapsed();

            let flops = 120.0 * nc as f64;
            let mflops = flops / (total.as_secs_f64() / self.iterations as f64) / 1e6;

            let r = GpuBenchReport {
                name: format!("lbm_bgk_{}x{}x{}", nx, ny, nz),
                backend: BackendKind::Cpu,
                n: nc,
                iterations: self.iterations,
                total,
                mean: total / self.iterations,
                mflops: Some(mflops),
            };
            out.push(r.clone());
            self.reports.push(r);
        }

        out
    }

    // ── Particle scan benchmark ───────────────────────────────────────────────

    /// Benchmark parallel prefix scan on `n` f64 elements (CPU Rayon scan).
    ///
    /// FLOPs = 2n (N adds in up-sweep + N adds in down-sweep).
    pub fn bench_parallel_scan(&mut self, n: usize) -> GpuBenchReport {
        let data: Vec<f64> = (0..n).map(|i| i as f64 + 1.0).collect();

        for _ in 0..self.warmup {
            let _ = inclusive_scan_cpu(&data);
        }
        let t0 = Instant::now();
        let mut result = Vec::new();
        for _ in 0..self.iterations {
            result = inclusive_scan_cpu(&data);
        }
        let total = t0.elapsed();
        let _ = result;

        let flops = 2.0 * n as f64;
        let mflops = flops / (total.as_secs_f64() / self.iterations as f64) / 1e6;

        let r = GpuBenchReport {
            name: "parallel_scan".to_string(),
            backend: BackendKind::Cpu,
            n,
            iterations: self.iterations,
            total,
            mean: total / self.iterations,
            mflops: Some(mflops),
        };
        self.reports.push(r.clone());
        r
    }

    // ── Full suite ────────────────────────────────────────────────────────────

    /// Run the complete GPU benchmark suite and return a formatted summary.
    ///
    /// ```
    /// use oxiphysics_gpu::gpu_bench::GpuBenchHarness;
    /// let mut h = GpuBenchHarness::new();
    /// let summary = h.run_full_suite();
    /// assert!(!summary.is_empty());
    /// ```
    pub fn run_full_suite(&mut self) -> String {
        self.bench_sph_density(64);
        self.bench_sph_density(256);
        self.bench_lbm_step(8, 8, 8);
        self.bench_lbm_step(16, 16, 4);
        self.bench_parallel_scan(1024);
        self.bench_parallel_scan(65536);

        let mut out = format!("{} benchmarks\n", self.reports.len());
        for r in &self.reports {
            out.push_str(&format!("  {}\n", r));
        }
        out
    }

    /// Benchmark CPU inclusive scan vs wgpu copy dispatch for `n` f64 elements.
    ///
    /// Both sides operate on the same data (a ramp of 0.0..n).  The wgpu side
    /// dispatches a copy shader (since f32 on-device means scan parity is a
    /// different test).  Returns a `Vec` with one CPU report, and optionally one
    /// wgpu report if an adapter is available.
    ///
    /// If no GPU adapter is present, only the CPU report is returned (no panic).
    ///
    /// ```
    /// use oxiphysics_gpu::gpu_bench::GpuBenchHarness;
    /// let mut h = GpuBenchHarness::new();
    /// let reports = h.cpu_vs_wgpu_comparison(1000);
    /// assert!(!reports.is_empty());
    /// assert_eq!(reports[0].name, "cpu_copy_scan");
    /// ```
    pub fn cpu_vs_wgpu_comparison(&mut self, n: usize) -> Vec<GpuBenchReport> {
        let mut out = Vec::new();

        // ── CPU path: inclusive scan ──────────────────────────────────────────
        let data: Vec<f64> = (0..n).map(|i| i as f64).collect();
        for _ in 0..self.warmup {
            let _ = inclusive_scan_cpu(&data);
        }
        let t0 = std::time::Instant::now();
        for _ in 0..self.iterations {
            let _ = inclusive_scan_cpu(&data);
        }
        let total_cpu = t0.elapsed();
        let mean_cpu = total_cpu / self.iterations;
        let flops = 2.0 * n as f64;
        let mflops_cpu = flops / (total_cpu.as_secs_f64() / self.iterations as f64) / 1e6;

        let cpu_report = GpuBenchReport {
            name: "cpu_copy_scan".to_string(),
            backend: BackendKind::Cpu,
            n,
            iterations: self.iterations,
            total: total_cpu,
            mean: mean_cpu,
            mflops: Some(mflops_cpu),
        };
        out.push(cpu_report.clone());
        self.reports.push(cpu_report);

        // ── wgpu path (feature-gated) ─────────────────────────────────────────
        #[cfg(feature = "wgpu-backend")]
        {
            use crate::compute::wgpu_backend::real::WgpuBackendReal;

            let backend_result = WgpuBackendReal::try_new();
            if let Ok(mut backend) = backend_result
                && let Ok(workgroups) = WgpuBackendReal::dispatch_count_for(n, 64)
            {
                const COPY_WGSL: &str = r#"
@group(0) @binding(0) var<storage, read>       in_buf:  array<f32>;
@group(0) @binding(1) var<storage, read_write> out_buf: array<f32>;

@compute @workgroup_size(64)
fn main(
    @builtin(global_invocation_id) gid: vec3<u32>,
    @builtin(num_workgroups) nwg: vec3<u32>,
) {
    let i = gid.x + gid.y * nwg.x * 64u;
    if (i < arrayLength(&in_buf)) {
        out_buf[i] = in_buf[i];
    }
}
"#;
                let in_buf = backend.create_buffer_f64(n);
                let out_buf = backend.create_buffer_f64(n);
                backend.write_buffer_f64(in_buf, &data);

                // Warm-up
                for _ in 0..self.warmup {
                    let _ = backend.dispatch_wgsl(
                        COPY_WGSL,
                        "main",
                        &[
                            (in_buf, wgpu::BufferBindingType::Storage { read_only: true }),
                            (
                                out_buf,
                                wgpu::BufferBindingType::Storage { read_only: false },
                            ),
                        ],
                        workgroups,
                    );
                }

                let t0 = std::time::Instant::now();
                for _ in 0..self.iterations {
                    let _ = backend.dispatch_wgsl(
                        COPY_WGSL,
                        "main",
                        &[
                            (in_buf, wgpu::BufferBindingType::Storage { read_only: true }),
                            (
                                out_buf,
                                wgpu::BufferBindingType::Storage { read_only: false },
                            ),
                        ],
                        workgroups,
                    );
                }
                let total_wgpu = t0.elapsed();
                let mean_wgpu = total_wgpu / self.iterations;
                let mflops_wgpu = flops / (total_wgpu.as_secs_f64() / self.iterations as f64) / 1e6;

                let wgpu_report = GpuBenchReport {
                    name: "wgpu_copy_dispatch".to_string(),
                    backend: BackendKind::Wgpu,
                    n,
                    iterations: self.iterations,
                    total: total_wgpu,
                    mean: mean_wgpu,
                    mflops: Some(mflops_wgpu),
                };
                out.push(wgpu_report.clone());
                self.reports.push(wgpu_report);
            }
        }

        out
    }

    /// Benchmark the SPH density kernel on CPU and wgpu backends side-by-side.
    ///
    /// Builds an SPH simulation with `n` particles arranged in a uniform grid
    /// inside the domain `[-10, 10]³`. Runs the full `SphSimulation::step`
    /// (density + pressure + accel + integrate) on both backends and returns
    /// timing reports.
    ///
    /// If no wgpu adapter is available, only the CPU report is returned.
    ///
    /// # Example
    /// ```
    /// let mut h = oxiphysics_gpu::gpu_bench::GpuBenchHarness::new();
    /// let reports = h.cpu_vs_wgpu_sph(64);
    /// assert!(!reports.is_empty());
    /// ```
    pub fn cpu_vs_wgpu_sph(&mut self, n: usize) -> Vec<GpuBenchReport> {
        let cfg = sph_bench_config(n);
        let mut out = Vec::new();

        // ── CPU path ──────────────────────────────────────────────────────────
        // Pinned to the CPU reference: no GPU backend is ever initialised.
        let mut cpu_sim = build_cpu_sph(&cfg, n);
        if let Ok(total) = time_sph_steps(&mut cpu_sim, self.warmup, self.iterations) {
            let r = GpuBenchReport {
                name: "sph_density_cpu".to_string(),
                backend: BackendKind::Cpu,
                n,
                iterations: self.iterations,
                total,
                mean: total / self.iterations,
                mflops: None,
            };
            out.push(r.clone());
            self.reports.push(r);
        }

        // ── wgpu path (only when every step really ran on the GPU) ────────────
        if let Some(mut gpu_sim) = build_gpu_sph(&cfg, n)
            && let Ok(total) = time_sph_steps(&mut gpu_sim, self.warmup, self.iterations)
        {
            let r = GpuBenchReport {
                name: "sph_density_wgpu".to_string(),
                backend: BackendKind::Wgpu,
                n,
                iterations: self.iterations,
                total,
                mean: total / self.iterations,
                mflops: None,
            };
            out.push(r.clone());
            self.reports.push(r);
        }

        out
    }

    /// Run SPH density on CPU and (optionally) CUDA; return timing reports.
    ///
    /// The CPU path runs the same `SphSimulation::step` loop as
    /// [`Self::cpu_vs_wgpu_sph`] but is tagged with `"cuda_sph_density_cpu"`.
    ///
    /// Under the `cuda-backend` feature, a second report is added when a CUDA
    /// device is available at runtime.  If no CUDA driver is present (e.g. on
    /// macOS) only the CPU report is returned — no panic.
    ///
    /// # Example
    /// ```
    /// let mut h = oxiphysics_gpu::gpu_bench::GpuBenchHarness::new();
    /// let reports = h.cpu_vs_cuda_sph(64);
    /// assert!(!reports.is_empty());
    /// assert!(reports[0].name.contains("sph_density"));
    /// assert!(reports[0].mean > std::time::Duration::ZERO);
    /// ```
    pub fn cpu_vs_cuda_sph(&mut self, n: usize) -> Vec<GpuBenchReport> {
        let cfg = sph_bench_config(n);
        let mut out = Vec::new();

        // ── CPU path ──────────────────────────────────────────────────────────
        // Pinned to the CPU reference: no GPU backend is ever initialised.
        let mut cpu_sim = build_cpu_sph(&cfg, n);
        if let Ok(total) = time_sph_steps(&mut cpu_sim, self.warmup, self.iterations) {
            let r = GpuBenchReport {
                name: "cuda_sph_density_cpu".to_string(),
                backend: BackendKind::Cpu,
                n,
                iterations: self.iterations,
                total,
                mean: total / self.iterations,
                mflops: None,
            };
            out.push(r.clone());
            self.reports.push(r);
        }

        // ── CUDA path (feature + runtime gated) ───────────────────────────────
        #[cfg(feature = "cuda-backend")]
        {
            use crate::compute::cuda_backend::{CUDA_SPH_DENSITY_SRC, CudaBackend};

            if let Ok(mut backend) = CudaBackend::try_new(0) {
                // Compile and register the SPH density kernel.
                let compiled =
                    backend.compile_and_register("sph_density_kernel", CUDA_SPH_DENSITY_SRC);
                if compiled.is_ok() {
                    // Build position buffer (n × 3 doubles, interleaved xyz).
                    let side = (n as f64).cbrt().ceil() as usize + 1;
                    let mut positions = vec![0.0_f64; n * 3];
                    for idx in 0..n {
                        positions[3 * idx] = (idx % side) as f64 * 0.1 - 5.0;
                        positions[3 * idx + 1] = ((idx / side) % side) as f64 * 0.1;
                        positions[3 * idx + 2] = (idx / (side * side)) as f64 * 0.1;
                    }

                    let pos_buf = backend.create_buffer(n * 3);
                    let den_buf = backend.create_buffer(n);
                    backend.write_buffer(pos_buf, &positions);

                    let block_x: u32 = 256;
                    let grid_x = (n as u32).div_ceil(block_x);

                    // The kernel signature is:
                    //   sph_density_kernel(const double*, double*, int, double, double)
                    // so we forward (n_particles, smoothing_h, particle_mass)
                    // as scalar arguments after the two buffer arguments.
                    let n_i32 = [n as i32];
                    let scalars_f64 = [
                        cfg.smoothing_h,
                        // Use a unit particle mass if the config left it as
                        // zero (cpu_vs_cuda_sph does not run SphSimulation::new
                        // for the GPU path, so the default 0.0 mass is fine to
                        // override for a smoke benchmark).
                        if cfg.particle_mass > 0.0 {
                            cfg.particle_mass
                        } else {
                            1.0
                        },
                    ];

                    // Warm-up
                    for _ in 0..self.warmup {
                        backend.launch_with_scalars(
                            "sph_density_kernel",
                            &[pos_buf, den_buf],
                            &n_i32,
                            &scalars_f64,
                            grid_x,
                            block_x,
                        );
                        backend.synchronize();
                    }

                    let t0 = Instant::now();
                    for _ in 0..self.iterations {
                        backend.launch_with_scalars(
                            "sph_density_kernel",
                            &[pos_buf, den_buf],
                            &n_i32,
                            &scalars_f64,
                            grid_x,
                            block_x,
                        );
                        backend.synchronize();
                    }
                    let total = t0.elapsed();

                    let r = GpuBenchReport {
                        name: "cuda_sph_density_gpu".to_string(),
                        backend: BackendKind::Cuda,
                        n,
                        iterations: self.iterations,
                        total,
                        mean: total / self.iterations,
                        mflops: None,
                    };
                    out.push(r.clone());
                    self.reports.push(r);
                }
            }
        }

        out
    }

    /// Print a comparison table for all collected reports.
    pub fn print_comparison(&self) {
        println!("\n{:=<75}", "");
        println!(
            "{:<5} {:<22} {:>8} {:>12} {:>10}",
            "Back", "Kernel", "N", "Mean (µs)", "MFLOPs"
        );
        println!("{:=<75}", "");
        for r in &self.reports {
            let mf = r.mflops.map_or("—".to_string(), |m| format!("{:.1}", m));
            println!(
                "{:<5} {:<22} {:>8} {:>12.3} {:>10}",
                r.backend,
                r.name,
                r.n,
                r.mean.as_secs_f64() * 1e6,
                mf
            );
        }
        println!("{:=<75}", "");
    }
}

impl Default for GpuBenchHarness {
    fn default() -> Self {
        Self::new()
    }
}

// ── SpeedupReport ─────────────────────────────────────────────────────────────

/// Computed speedup between a CPU and a wgpu benchmark report.
#[derive(Debug, Clone)]
pub struct SpeedupReport {
    /// Mean time for the baseline (CPU) backend.
    pub cpu_mean: Duration,
    /// Mean time for the accelerated (wgpu) backend, if available.
    pub wgpu_mean: Option<Duration>,
    /// Speedup ratio = cpu_mean / wgpu_mean, if wgpu was measured.
    pub speedup: Option<f64>,
}

/// Compute a speedup ratio from a pair of bench reports.
///
/// Expects `reports[0]` to be the CPU report and `reports[1]` (if present)
/// to be the wgpu report.  Returns `SpeedupReport { speedup: None }` if
/// only one report is present (GPU unavailable).
///
/// # Example
/// ```
/// use oxiphysics_gpu::gpu_bench::{GpuBenchHarness, compute_speedup};
/// let mut h = GpuBenchHarness::new();
/// let reports = h.cpu_vs_wgpu_sph(64);
/// let sr = compute_speedup(&reports);
/// assert!(sr.cpu_mean.as_secs_f64() > 0.0);
/// ```
pub fn compute_speedup(reports: &[GpuBenchReport]) -> SpeedupReport {
    let cpu_mean = reports.first().map(|r| r.mean).unwrap_or(Duration::ZERO);
    let wgpu_mean = reports.get(1).map(|r| r.mean);
    let speedup = wgpu_mean.map(|wm| {
        if wm.as_secs_f64() > 0.0 {
            cpu_mean.as_secs_f64() / wm.as_secs_f64()
        } else {
            f64::INFINITY
        }
    });
    SpeedupReport {
        cpu_mean,
        wgpu_mean,
        speedup,
    }
}

// ── CudaSpeedupReport ─────────────────────────────────────────────────────────

/// Computed speedup between a CPU and a CUDA benchmark report.
#[derive(Debug, Clone)]
pub struct CudaSpeedupReport {
    /// Mean time for the baseline (CPU) backend.
    pub cpu_mean: Duration,
    /// Mean time for the CUDA backend, if available.
    pub cuda_mean: Option<Duration>,
    /// Speedup ratio = cpu_mean / cuda_mean, if CUDA was measured.
    pub speedup: Option<f64>,
}

/// Compute a speedup ratio from a pair of bench reports produced by
/// [`GpuBenchHarness::cpu_vs_cuda_sph`].
///
/// Expects `reports[0]` to be the CPU report and `reports[1]` (if present)
/// to be the CUDA report.  Returns `CudaSpeedupReport { speedup: None }` if
/// only one report is present (CUDA unavailable).
///
/// # Example
/// ```
/// use oxiphysics_gpu::gpu_bench::{GpuBenchHarness, compute_cuda_speedup};
/// let mut h = GpuBenchHarness::new();
/// let reports = h.cpu_vs_cuda_sph(64);
/// let sr = compute_cuda_speedup(&reports);
/// assert!(sr.cpu_mean.as_secs_f64() >= 0.0);
/// ```
pub fn compute_cuda_speedup(reports: &[GpuBenchReport]) -> CudaSpeedupReport {
    let cpu_mean = reports.first().map(|r| r.mean).unwrap_or(Duration::ZERO);
    let cuda_mean = reports.get(1).map(|r| r.mean);
    let speedup = cuda_mean.map(|cm| {
        if cm.as_secs_f64() > 0.0 {
            cpu_mean.as_secs_f64() / cm.as_secs_f64()
        } else {
            f64::INFINITY
        }
    });
    CudaSpeedupReport {
        cpu_mean,
        cuda_mean,
        speedup,
    }
}

// ── CPU helpers ───────────────────────────────────────────────────────────────

/// Sequential inclusive prefix scan (Σ) on `f64` elements.
///
/// Returns a `Vec<f64>` where `out[i] = Σ_{j≤i} data[j]`.
pub fn inclusive_scan_cpu(data: &[f64]) -> Vec<f64> {
    let mut out = Vec::with_capacity(data.len());
    let mut acc = 0.0_f64;
    for &v in data {
        acc += v;
        out.push(acc);
    }
    out
}

// ── tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_available_backends_has_cpu() {
        let b = GpuBenchHarness::available_backends();
        assert!(b.contains(&BackendKind::Cpu));
    }

    /// wgpu is listed exactly when a real adapter can be created: the old
    /// detection went through the always-failing CPU-shadow type.
    #[test]
    fn test_available_backends_reports_wgpu_when_an_adapter_exists() {
        let listed = GpuBenchHarness::available_backends().contains(&BackendKind::Wgpu);
        #[cfg(feature = "wgpu-backend")]
        assert_eq!(
            listed,
            crate::compute::wgpu_backend::real::WgpuBackendReal::try_new().is_ok()
        );
        #[cfg(not(feature = "wgpu-backend"))]
        assert!(!listed);
    }

    #[test]
    fn test_inclusive_scan() {
        let data = vec![1.0, 2.0, 3.0, 4.0];
        let out = inclusive_scan_cpu(&data);
        assert_eq!(out, vec![1.0, 3.0, 6.0, 10.0]);
    }

    #[test]
    fn test_bench_sph_density_returns_at_least_cpu() {
        let mut h = GpuBenchHarness {
            warmup: 0,
            iterations: 1,
            reports: Vec::new(),
        };
        let reports = h.bench_sph_density(8);
        assert!(!reports.is_empty());
        assert_eq!(reports[0].backend, BackendKind::Cpu);
    }

    #[test]
    fn test_bench_lbm_step() {
        let mut h = GpuBenchHarness {
            warmup: 0,
            iterations: 1,
            reports: Vec::new(),
        };
        let reports = h.bench_lbm_step(4, 4, 4);
        assert_eq!(reports.len(), 1);
        assert_eq!(reports[0].n, 64);
    }

    #[test]
    fn test_bench_parallel_scan() {
        let mut h = GpuBenchHarness {
            warmup: 0,
            iterations: 1,
            reports: Vec::new(),
        };
        let r = h.bench_parallel_scan(100);
        assert_eq!(r.n, 100);
        assert!(r.mflops.is_some());
    }

    #[test]
    fn test_run_full_suite() {
        let mut h = GpuBenchHarness {
            warmup: 0,
            iterations: 1,
            reports: Vec::new(),
        };
        let summary = h.run_full_suite();
        assert!(summary.contains("benchmarks"));
    }

    #[test]
    fn test_backend_display() {
        assert_eq!(format!("{}", BackendKind::Cpu), "CPU");
        assert_eq!(format!("{}", BackendKind::Wgpu), "wgpu");
        assert_eq!(format!("{}", BackendKind::Cuda), "CUDA");
    }

    /// The simulations behind every `BackendKind::Cpu` report must be pinned to
    /// the CPU reference: on a host with a GPU adapter the default
    /// configuration would silently run on the GPU.
    #[test]
    fn test_cpu_bench_simulations_never_touch_the_gpu() {
        let cfg = sph_bench_config(8);
        let mut sph = build_cpu_sph(&cfg, 8);
        assert_eq!(sph.config.compute, ComputePreference::CpuOnly);
        assert!(!sph.has_gpu());
        time_sph_steps(&mut sph, 1, 2).expect("CPU steps cannot fail");
        assert!(!sph.has_gpu());
        assert!(sph.gpu_error().is_none(), "the GPU must not even be probed");

        let mut lbm = LbmSimulation::new(lbm_bench_config(4, 4, 4));
        lbm.set_lid_velocity(0.1, 0.0, 0.0);
        for _ in 0..3 {
            lbm.step();
        }
        assert!(!lbm.has_gpu());
        assert!(lbm.gpu_error().is_none(), "the GPU must not even be probed");
        assert_eq!(lbm.step_count, 3);
    }

    /// A GPU benchmark simulation either really has a GPU pipeline and runs its
    /// steps on it, or does not exist (no CPU duplicate under a GPU label).
    #[test]
    fn test_gpu_bench_simulation_is_gpu_or_absent() {
        let cfg = sph_bench_config(8);
        if let Some(mut sim) = build_gpu_sph(&cfg, 8) {
            assert_eq!(sim.config.compute, ComputePreference::GpuRequired);
            assert!(sim.has_gpu());
            time_sph_steps(&mut sim, 0, 2).expect("the GPU takes every step");
            assert!(sim.has_gpu());
        }
    }

    /// Every wgpu report must come from a GPU-required run, and the CPU report
    /// always leads.
    #[test]
    fn test_sph_reports_label_backends_truthfully() {
        let mut h = GpuBenchHarness {
            warmup: 0,
            iterations: 1,
            reports: Vec::new(),
        };
        let reports = h.cpu_vs_wgpu_sph(8);
        assert_eq!(reports[0].backend, BackendKind::Cpu);
        assert_eq!(reports[0].name, "sph_density_cpu");
        assert!(reports.len() <= 2);
        if let Some(gpu) = reports.get(1) {
            assert_eq!(gpu.backend, BackendKind::Wgpu);
            assert_eq!(gpu.name, "sph_density_wgpu");
            assert!(build_gpu_sph(&sph_bench_config(8), 8).is_some());
        }
    }

    #[test]
    fn test_seed_sph_grid_layout() {
        let cfg = sph_bench_config(8);
        let sim = build_cpu_sph(&cfg, 8);
        // side = ceil(cbrt(8)) + 1 = 3: x runs fastest, then y, then z.
        let expected = |idx: usize| {
            [
                (idx % 3) as f64 * 0.1 - 5.0,
                ((idx / 3) % 3) as f64 * 0.1,
                (idx / 9) as f64 * 0.1,
            ]
        };
        for idx in 0..8 {
            let want = expected(idx);
            assert_eq!(sim.state.pos_x[idx], want[0]);
            assert_eq!(sim.state.pos_y[idx], want[1]);
            assert_eq!(sim.state.pos_z[idx], want[2]);
        }
    }
}
