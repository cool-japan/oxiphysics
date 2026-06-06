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

use std::time::{Duration, Instant};

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

        // Try wgpu — succeeds when the GPU driver is present.
        if crate::compute::WgpuBackend::try_new().is_ok() {
            out.push(BackendKind::Wgpu);
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
        let cfg = SphConfig {
            n_particles: n,
            smoothing_h: 0.1,
            rest_density: 1000.0,
            gravity: 0.0, // no gravity — pure density bench
            domain_min: [-10.; 3],
            domain_max: [10.; 3],
            ..SphConfig::default()
        };

        let mut out = Vec::new();

        // ── CPU path ──────────────────────────────────────────────────────────
        {
            // Build a CPU-only sim (no GPU backend will be tried)
            let mut sim = SphSimulation::new(cfg.clone());
            // Scatter particles in a regular grid
            let side = (n as f64).cbrt().ceil() as usize + 1;
            for (idx, i) in (0..n).enumerate() {
                let x = (idx % side) as f64 * 0.1 - 5.0;
                let y = ((idx / side) % side) as f64 * 0.1;
                let z = (idx / (side * side)) as f64 * 0.1;
                sim.state.pos_x[i] = x;
                sim.state.pos_y[i] = y;
                sim.state.pos_z[i] = z;
            }

            // Warm-up
            for _ in 0..self.warmup {
                sim.step(1.0 / 60.0);
            }

            let t0 = Instant::now();
            for _ in 0..self.iterations {
                sim.step(1.0 / 60.0);
            }
            let total = t0.elapsed();

            let flops = 10.0 * n as f64 * n as f64;
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

        // ── wgpu path (if available) ──────────────────────────────────────────
        if crate::compute::WgpuBackend::try_new().is_ok() {
            let mut sim = SphSimulation::new(cfg.clone());
            let side = (n as f64).cbrt().ceil() as usize + 1;
            for (idx, i) in (0..n).enumerate() {
                sim.state.pos_x[i] = (idx % side) as f64 * 0.1 - 5.0;
                sim.state.pos_y[i] = ((idx / side) % side) as f64 * 0.1;
                sim.state.pos_z[i] = (idx / (side * side)) as f64 * 0.1;
            }

            for _ in 0..self.warmup {
                sim.step(1.0 / 60.0);
            }
            let t0 = Instant::now();
            for _ in 0..self.iterations {
                sim.step(1.0 / 60.0);
            }
            let total = t0.elapsed();

            let backend = if sim.has_gpu() {
                BackendKind::Wgpu
            } else {
                BackendKind::Cpu
            };
            let flops = 10.0 * n as f64 * n as f64;
            let mflops = flops / (total.as_secs_f64() / self.iterations as f64) / 1e6;

            let r = GpuBenchReport {
                name: "sph_density".to_string(),
                backend,
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
        let cfg = LbmConfig {
            nx,
            ny,
            nz,
            tau: 0.6,
            rho0: 1.0,
            force_x: 0.0,
            force_y: 0.0,
            force_z: 0.0,
        };
        let nc = nx * ny * nz;
        let mut out = Vec::new();

        // CPU path
        {
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
            if let Ok(mut backend) = backend_result {
                const COPY_WGSL: &str = r#"
@group(0) @binding(0) var<storage, read>       in_buf:  array<f32>;
@group(0) @binding(1) var<storage, read_write> out_buf: array<f32>;

@compute @workgroup_size(64)
fn main(@builtin(global_invocation_id) gid: vec3<u32>) {
    let i = gid.x;
    if (i < arrayLength(&in_buf)) {
        out_buf[i] = in_buf[i];
    }
}
"#;
                let in_buf = backend.create_buffer_f64(n);
                let out_buf = backend.create_buffer_f64(n);
                backend.write_buffer_f64(in_buf, &data);

                let workgroups = WgpuBackendReal::dispatch_count_for(n, 64);

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
        let cfg = SphConfig {
            n_particles: n,
            smoothing_h: 0.1,
            rest_density: 1000.0,
            gravity: 0.0,
            domain_min: [-10.; 3],
            domain_max: [10.; 3],
            ..SphConfig::default()
        };

        let mut out = Vec::new();

        // ── CPU path ──────────────────────────────────────────────────────────
        {
            let mut sim = SphSimulation::new(cfg.clone());
            let side = (n as f64).cbrt().ceil() as usize + 1;
            for (idx, i) in (0..n).enumerate() {
                let x = (idx % side) as f64 * 0.1 - 5.0;
                let y = ((idx / side) % side) as f64 * 0.1;
                let z = (idx / (side * side)) as f64 * 0.1;
                sim.state.pos_x[i] = x;
                sim.state.pos_y[i] = y;
                sim.state.pos_z[i] = z;
            }
            for _ in 0..self.warmup {
                sim.step(1.0 / 60.0);
            }
            let t0 = Instant::now();
            for _ in 0..self.iterations {
                sim.step(1.0 / 60.0);
            }
            let total = t0.elapsed();

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

        // ── wgpu path (runtime-gated) ─────────────────────────────────────────
        if crate::compute::WgpuBackend::try_new().is_ok() {
            let mut sim = SphSimulation::new(cfg.clone());
            let side = (n as f64).cbrt().ceil() as usize + 1;
            for (idx, i) in (0..n).enumerate() {
                sim.state.pos_x[i] = (idx % side) as f64 * 0.1 - 5.0;
                sim.state.pos_y[i] = ((idx / side) % side) as f64 * 0.1;
                sim.state.pos_z[i] = (idx / (side * side)) as f64 * 0.1;
            }

            for _ in 0..self.warmup {
                sim.step(1.0 / 60.0);
            }
            let t0 = Instant::now();
            for _ in 0..self.iterations {
                sim.step(1.0 / 60.0);
            }
            let total = t0.elapsed();

            let backend = if sim.has_gpu() {
                BackendKind::Wgpu
            } else {
                BackendKind::Cpu
            };

            let r = GpuBenchReport {
                name: "sph_density_wgpu".to_string(),
                backend,
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
        let cfg = crate::sph_gpu::SphConfig {
            n_particles: n,
            smoothing_h: 0.1,
            rest_density: 1000.0,
            gravity: 0.0,
            domain_min: [-10.; 3],
            domain_max: [10.; 3],
            ..crate::sph_gpu::SphConfig::default()
        };
        let mut out = Vec::new();

        // ── CPU path ──────────────────────────────────────────────────────────
        {
            let mut sim = crate::sph_gpu::SphSimulation::new(cfg.clone());
            let side = (n as f64).cbrt().ceil() as usize + 1;
            for idx in 0..n {
                let x = (idx % side) as f64 * 0.1 - 5.0;
                let y = ((idx / side) % side) as f64 * 0.1;
                let z = (idx / (side * side)) as f64 * 0.1;
                sim.state.pos_x[idx] = x;
                sim.state.pos_y[idx] = y;
                sim.state.pos_z[idx] = z;
            }
            for _ in 0..self.warmup {
                sim.step(1.0 / 60.0);
            }
            let t0 = Instant::now();
            for _ in 0..self.iterations {
                sim.step(1.0 / 60.0);
            }
            let total = t0.elapsed();

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
}
