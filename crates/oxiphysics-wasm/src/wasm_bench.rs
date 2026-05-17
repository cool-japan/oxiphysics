// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! WASM physics engine performance benchmarks.
//!
//! Measures simulation throughput for the three main simulation paths:
//! - Rigid-body stepping (`WasmPhysicsEngine::step`)
//! - Particle-system stepping (`WasmParticleSystem::step`)
//! - SPH fluid stepping (`WasmSphSim::step`)
//!
//! ## Platform notes
//!
//! On `wasm32-unknown-unknown` `std::time::Instant` is not available.
//! All timing-dependent code in this module is gated behind
//! `#[cfg(not(target_arch = "wasm32"))]`.  On WASM targets only the
//! JS-friendly `WasmBenchResult` struct and its methods are available.
//!
//! ## Usage (native)
//!
//! ```no_run
//! use oxiphysics_wasm::wasm_bench::run_wasm_bench_suite;
//! let report = run_wasm_bench_suite(true);
//! println!("{}", report);
//! ```

use wasm_bindgen::prelude::*;

use crate::engine::{WasmPhysicsEngine, WasmSphSim};
#[cfg(target_arch = "wasm32")]
use crate::particle_system::{WasmParticleConfig, WasmParticleSystem};

// ============================================================================
// WasmBenchResult — JS-friendly flat result (always available)
// ============================================================================

/// A flat benchmark result suitable for JavaScript consumption.
///
/// Timing is expressed in seconds (f64) so that JS-side `performance.now()`
/// measurements can be directly stored here.
#[wasm_bindgen]
#[derive(Debug, Clone)]
pub struct WasmBenchResult {
    /// Human-readable benchmark name.
    pub n: u32,
    /// Number of simulation steps measured.
    pub iterations: u32,
    /// Total wall-clock time in seconds.
    pub total_secs: f64,
    /// Mean time per step in seconds.
    pub mean_secs: f64,
    /// Steps per second.
    pub steps_per_sec: f64,
    /// Simulated entities per second.
    pub throughput: f64,
    /// Optional MFLOPS estimate (`-1.0` = unknown).
    pub mflops: f64,
    name: String,
}

#[wasm_bindgen]
impl WasmBenchResult {
    /// Create from raw timing values.
    pub fn new(name: String, n: u32, iterations: u32, total_secs: f64, mflops: f64) -> Self {
        let mean_secs = total_secs / (iterations.max(1) as f64);
        let steps_per_sec = iterations as f64 / total_secs.max(1e-12);
        let throughput = steps_per_sec * n as f64;
        Self {
            name,
            n,
            iterations,
            total_secs,
            mean_secs,
            steps_per_sec,
            throughput,
            mflops,
        }
    }

    /// Return the benchmark name.
    pub fn name(&self) -> String {
        self.name.clone()
    }

    /// Serialize to a JSON string.
    pub fn to_json(&self) -> String {
        format!(
            r#"{{"name":"{}","n":{},"iterations":{},"totalSecs":{},"meanSecs":{},"stepsPerSec":{},"throughput":{},"mflops":{}}}"#,
            self.name,
            self.n,
            self.iterations,
            self.total_secs,
            self.mean_secs,
            self.steps_per_sec,
            self.throughput,
            self.mflops,
        )
    }
}

// ============================================================================
// RapierReference — comparison table (always available, no timing dependency)
// ============================================================================

/// Reference throughput numbers from Rapier 0.18 WASM benchmarks (2024).
#[wasm_bindgen]
pub struct RapierReference;

#[wasm_bindgen]
impl RapierReference {
    /// Bodies per second for Rapier's "falling cubes" benchmark (WASM).
    pub fn falling_cubes_throughput() -> f64 {
        50_000.0
    }

    /// Particles per second for a basic verlet particle system.
    pub fn particle_throughput() -> f64 {
        5_000_000.0
    }

    /// SPH particles per second for a 1000-particle sim.
    pub fn sph_throughput() -> f64 {
        30_000.0
    }
}

// ============================================================================
// Timing-dependent code — native only
// ============================================================================

/// Native-only benchmark internals using `std::time::Instant`.
#[cfg(not(target_arch = "wasm32"))]
pub mod native {
    use std::time::{Duration, Instant};

    use super::{WasmBenchResult, WasmPhysicsEngine, WasmSphSim};
    use crate::particle_system::{WasmParticleConfig, WasmParticleSystem};

    /// Timing report for a single benchmark run with `std::time` types.
    #[derive(Debug, Clone)]
    pub struct WasmBenchReport {
        /// Human-readable benchmark name.
        pub name: String,
        /// Problem size.
        pub n: usize,
        /// Number of measured iterations.
        pub iterations: u32,
        /// Total wall-clock time.
        pub total: Duration,
        /// Mean time per step.
        pub mean: Duration,
        /// Steps per second.
        pub steps_per_sec: f64,
        /// Simulated entities per second.
        pub throughput: f64,
        /// Optional MFLOPS estimate.
        pub mflops: Option<f64>,
    }

    impl WasmBenchReport {
        pub(super) fn new(
            name: &str,
            n: usize,
            iterations: u32,
            total: Duration,
            mflops: Option<f64>,
        ) -> Self {
            let mean = total / iterations.max(1);
            let steps_per_sec = iterations as f64 / total.as_secs_f64().max(1e-12);
            let throughput = steps_per_sec * n as f64;
            Self {
                name: name.to_string(),
                n,
                iterations,
                total,
                mean,
                steps_per_sec,
                throughput,
                mflops,
            }
        }

        /// Convert to the JS-friendly flat result.
        pub fn to_wasm_result(&self) -> WasmBenchResult {
            WasmBenchResult::new(
                self.name.clone(),
                self.n as u32,
                self.iterations,
                self.total.as_secs_f64(),
                self.mflops.unwrap_or(-1.0),
            )
        }
    }

    impl std::fmt::Display for WasmBenchReport {
        fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
            write!(
                f,
                "{:40}  n={:6}  iters={:4}  mean={:>9}  throughput={:.2e}/s",
                self.name,
                self.n,
                self.iterations,
                format_duration(self.mean),
                self.throughput,
            )?;
            if let Some(mf) = self.mflops {
                write!(f, "  {:.1} MFLOP/s", mf)?;
            }
            Ok(())
        }
    }

    fn format_duration(d: Duration) -> String {
        let ns = d.as_nanos();
        if ns < 1_000 {
            format!("{} ns", ns)
        } else if ns < 1_000_000 {
            format!("{:.1} µs", ns as f64 / 1_000.0)
        } else if ns < 1_000_000_000 {
            format!("{:.1} ms", ns as f64 / 1_000_000.0)
        } else {
            format!("{:.2} s", d.as_secs_f64())
        }
    }

    /// Generic timing harness.
    pub struct WasmBenchHarness {
        pub warmup: u32,
        pub iterations: u32,
        reports: Vec<WasmBenchReport>,
    }

    impl WasmBenchHarness {
        pub fn new(warmup: u32, iterations: u32) -> Self {
            Self {
                warmup,
                iterations,
                reports: Vec::new(),
            }
        }

        pub fn run<F>(&mut self, name: &str, n: usize, mflops: Option<f64>, mut f: F)
        where
            F: FnMut(),
        {
            for _ in 0..self.warmup {
                f();
            }
            let t0 = Instant::now();
            for _ in 0..self.iterations {
                f();
            }
            let total = t0.elapsed();
            self.reports.push(WasmBenchReport::new(
                name,
                n,
                self.iterations,
                total,
                mflops,
            ));
        }

        pub fn reports(&self) -> &[WasmBenchReport] {
            &self.reports
        }

        pub fn print_summary(&self) -> String {
            let mut out = String::from("=== WASM Physics Benchmark Suite ===\n");
            for r in &self.reports {
                out.push_str(&format!("{}\n", r));
            }
            out
        }
    }

    /// Benchmark rigid-body stepping with `n` falling spheres.
    pub fn bench_rigid_bodies(n: usize, iterations: u32) -> WasmBenchReport {
        let mut engine = WasmPhysicsEngine::new(0.0, -9.81, 0.0);
        let ground = engine.add_static_body(0.0, 0.0, 0.0);
        engine.add_plane_collider(ground, 0.0, 1.0, 0.0, 0.0);

        let cols = (n as f64).cbrt().ceil() as usize;
        let mut added = 0;
        'outer: for ix in 0..cols {
            for iy in 0..cols {
                for iz in 0..cols {
                    if added >= n {
                        break 'outer;
                    }
                    let x = ix as f64 * 1.1;
                    let y = 1.0 + iy as f64 * 1.1;
                    let z = iz as f64 * 1.1;
                    let b = engine.add_dynamic_body(1.0, x, y, z);
                    engine.add_sphere_collider(b, 0.5);
                    added += 1;
                }
            }
        }

        let dt = 1.0 / 60.0;
        for _ in 0..5 {
            engine.step(dt);
        }

        let t0 = Instant::now();
        for _ in 0..iterations {
            engine.step(dt);
        }
        let total = t0.elapsed();

        let mflops = (n as f64 * 30.0 * iterations as f64) / (total.as_secs_f64() * 1e6);
        WasmBenchReport::new(
            &format!("rigid_bodies(n={})", n),
            n,
            iterations,
            total,
            Some(mflops),
        )
    }

    /// Benchmark particle-system stepping with `n` particles.
    pub fn bench_particle_system(n: usize, iterations: u32) -> WasmBenchReport {
        let config = WasmParticleConfig {
            max_particles: n,
            lifetime_range: [5.0, 10.0],
            gravity: [0.0, -9.81, 0.0],
            ..WasmParticleConfig::default_earth()
        };
        let mut sys = WasmParticleSystem::new(config);
        sys.spawn(n);

        let dt = 1.0 / 60.0;
        for _ in 0..5 {
            sys.step(dt);
        }

        let t0 = Instant::now();
        for _ in 0..iterations {
            sys.step(dt);
        }
        let total = t0.elapsed();

        let mflops = (n as f64 * 15.0 * iterations as f64) / (total.as_secs_f64() * 1e6);
        WasmBenchReport::new(
            &format!("particle_system(n={})", n),
            n,
            iterations,
            total,
            Some(mflops),
        )
    }

    /// Benchmark SPH fluid stepping with `n` particles.
    pub fn bench_sph_fluid(n: usize, iterations: u32) -> WasmBenchReport {
        let mut sim = WasmSphSim::new();
        sim.configure(0.1, 1000.0, 200.0, 0.01, 0.02);
        sim.set_gravity(0.0, -9.81, 0.0);

        let cols = (n as f64).cbrt().ceil() as usize;
        let spacing = 0.08;
        let mut added = 0;
        'outer: for ix in 0..cols {
            for iy in 0..cols {
                for iz in 0..cols {
                    if added >= n {
                        break 'outer;
                    }
                    sim.add_particle(
                        ix as f64 * spacing,
                        iy as f64 * spacing,
                        iz as f64 * spacing,
                    );
                    added += 1;
                }
            }
        }

        let dt = 1e-4;
        for _ in 0..3 {
            sim.step(dt);
        }

        let t0 = Instant::now();
        for _ in 0..iterations {
            sim.step(dt);
        }
        let total = t0.elapsed();

        let mflops = (n as f64 * 200.0 * iterations as f64) / (total.as_secs_f64() * 1e6);
        WasmBenchReport::new(
            &format!("sph_fluid(n={})", n),
            n,
            iterations,
            total,
            Some(mflops),
        )
    }

    /// Run the complete WASM benchmark suite.
    pub fn run_wasm_bench_suite(verbose: bool) -> String {
        let sizes: &[(usize, u32)] = &[(100, 200), (500, 100), (1000, 50)];
        let sph_sizes: &[(usize, u32)] = &[(100, 200), (500, 50), (1000, 20)];

        let mut all_reports: Vec<WasmBenchReport> = Vec::new();
        let harness = WasmBenchHarness::new(3, 1);

        if verbose {
            eprintln!("[wasm_bench] Running rigid body benchmarks...");
        }
        for &(n, iters) in sizes {
            all_reports.push(bench_rigid_bodies(n, iters));
        }

        if verbose {
            eprintln!("[wasm_bench] Running particle system benchmarks...");
        }
        for &(n, iters) in sizes {
            all_reports.push(bench_particle_system(n, iters));
        }

        if verbose {
            eprintln!("[wasm_bench] Running SPH fluid benchmarks...");
        }
        for &(n, iters) in sph_sizes {
            all_reports.push(bench_sph_fluid(n, iters));
        }

        if verbose {
            eprintln!("[wasm_bench] Running WebGL bridge extraction benchmark...");
        }
        let webgl_report = {
            let mut engine = WasmPhysicsEngine::new(0.0, -9.81, 0.0);
            for i in 0..500 {
                let b = engine.add_dynamic_body(1.0, i as f64 * 0.1, 0.0, 0.0);
                engine.add_sphere_collider(b, 0.5);
            }
            let iters = 1000u32;
            for _ in 0..10 {
                let _ = crate::webgl_bridge::WebGlBridge::extract_render_frame(&engine);
            }
            let t0 = Instant::now();
            for _ in 0..iters {
                let _ = crate::webgl_bridge::WebGlBridge::extract_render_frame(&engine);
            }
            let total = t0.elapsed();
            WasmBenchReport::new(
                "webgl_bridge::extract_render_frame(n=500)",
                500,
                iters,
                total,
                None,
            )
        };
        all_reports.push(webgl_report);
        let _ = harness.reports();

        let mut out = String::from("=== OxiPhysics WASM Benchmark Suite ===\n");
        out.push_str("  warmup: 3 iters  |  platform: native (WASM target in browser)\n\n");
        for r in &all_reports {
            out.push_str(&format!("{}\n", r));
        }
        out.push('\n');

        // Comparison table
        out.push_str("=== OxiPhysics vs Rapier WASM Reference ===\n");
        out.push_str(&format!(
            "{:<32} {:>16} {:>18} {:>8}\n",
            "Benchmark", "OxiPhysics/s", "Rapier-ref/s", "Ratio"
        ));
        out.push_str(&"-".repeat(80));
        out.push('\n');

        let references: &[(&str, f64)] = &[
            ("rigid", super::RapierReference::falling_cubes_throughput()),
            ("particle", super::RapierReference::particle_throughput()),
            ("sph", super::RapierReference::sph_throughput()),
        ];
        for r in &all_reports {
            let ref_val = references
                .iter()
                .find(|(k, _)| r.name.contains(k))
                .map(|(_, v)| *v)
                .unwrap_or(0.0);
            let ratio = if ref_val > 0.0 {
                r.throughput / ref_val
            } else {
                0.0
            };
            out.push_str(&format!(
                "{:<32} {:>16.2e} {:>18.2e} {:>7.2}×\n",
                &r.name[..r.name.len().min(32)],
                r.throughput,
                ref_val,
                ratio
            ));
        }
        out
    }
}

// ============================================================================
// Public re-exports for convenience
// ============================================================================

#[cfg(not(target_arch = "wasm32"))]
pub use native::{
    WasmBenchHarness, WasmBenchReport, bench_particle_system, bench_rigid_bodies, bench_sph_fluid,
    run_wasm_bench_suite,
};

// ============================================================================
// WASM-facing entry points (JS-callable stubs that return WasmBenchResult)
// ============================================================================

/// Run a quick rigid-body benchmark on the WASM/JS side.
///
/// On native targets this delegates to `native::bench_rigid_bodies`.
/// On WASM this runs the sim but measures no time (returns 0.0 duration);
/// the caller should wrap it with `performance.now()` on the JS side.
#[wasm_bindgen]
pub fn wasm_bench_rigid_bodies(n: u32, iterations: u32) -> WasmBenchResult {
    #[cfg(not(target_arch = "wasm32"))]
    {
        native::bench_rigid_bodies(n as usize, iterations).to_wasm_result()
    }
    #[cfg(target_arch = "wasm32")]
    {
        // Run the simulation without native timing; JS side should measure with performance.now()
        let mut engine = WasmPhysicsEngine::new(0.0, -9.81, 0.0);
        let ground = engine.add_static_body(0.0, 0.0, 0.0);
        engine.add_plane_collider(ground, 0.0, 1.0, 0.0, 0.0);
        for i in 0..(n as usize) {
            let b = engine.add_dynamic_body(1.0, (i % 10) as f64, (i / 10) as f64 + 1.0, 0.0);
            engine.add_sphere_collider(b, 0.5);
        }
        let dt = 1.0 / 60.0;
        for _ in 0..iterations {
            engine.step(dt);
        }
        WasmBenchResult::new(format!("rigid_bodies(n={})", n), n, iterations, 0.0, -1.0)
    }
}

/// Run a quick particle-system benchmark.
#[wasm_bindgen]
pub fn wasm_bench_particle_system(n: u32, iterations: u32) -> WasmBenchResult {
    #[cfg(not(target_arch = "wasm32"))]
    {
        native::bench_particle_system(n as usize, iterations).to_wasm_result()
    }
    #[cfg(target_arch = "wasm32")]
    {
        let config = WasmParticleConfig {
            max_particles: n as usize,
            lifetime_range: [5.0, 10.0],
            gravity: [0.0, -9.81, 0.0],
            ..WasmParticleConfig::default_earth()
        };
        let mut sys = WasmParticleSystem::new(config);
        sys.spawn(n as usize);
        let dt = 1.0 / 60.0;
        for _ in 0..iterations {
            sys.step(dt);
        }
        WasmBenchResult::new(
            format!("particle_system(n={})", n),
            n,
            iterations,
            0.0,
            -1.0,
        )
    }
}

/// Run a quick SPH fluid benchmark.
#[wasm_bindgen]
pub fn wasm_bench_sph_fluid(n: u32, iterations: u32) -> WasmBenchResult {
    #[cfg(not(target_arch = "wasm32"))]
    {
        native::bench_sph_fluid(n as usize, iterations).to_wasm_result()
    }
    #[cfg(target_arch = "wasm32")]
    {
        let mut sim = WasmSphSim::new();
        sim.configure(0.1, 1000.0, 200.0, 0.01, 0.02);
        sim.set_gravity(0.0, -9.81, 0.0);
        for i in 0..(n as usize) {
            let x = (i % 10) as f64 * 0.08;
            let y = (i / 10) as f64 * 0.08;
            sim.add_particle(x, y, 0.0);
        }
        let dt = 1e-4;
        for _ in 0..iterations {
            sim.step(dt);
        }
        WasmBenchResult::new(format!("sph_fluid(n={})", n), n, iterations, 0.0, -1.0)
    }
}

// ============================================================================
// Unit tests
// ============================================================================

#[cfg(test)]
#[cfg(not(target_arch = "wasm32"))]
mod tests {
    use super::native::*;

    #[test]
    fn bench_rigid_bodies_runs() {
        let r = bench_rigid_bodies(10, 5);
        assert!(r.steps_per_sec > 0.0);
        assert!(r.throughput > 0.0);
        assert!(r.mflops.unwrap() > 0.0);
    }

    #[test]
    fn bench_particle_system_runs() {
        let r = bench_particle_system(50, 5);
        assert!(r.steps_per_sec > 0.0);
    }

    #[test]
    fn bench_sph_fluid_runs() {
        let r = bench_sph_fluid(20, 3);
        assert!(r.steps_per_sec > 0.0);
    }

    #[test]
    fn harness_records_reports() {
        let mut h = WasmBenchHarness::new(1, 3);
        let mut counter = 0u32;
        h.run("test", 100, Some(1.0), || {
            counter += 1;
        });
        assert_eq!(h.reports().len(), 1);
        assert_eq!(counter, 4); // 1 warmup + 3 measured
        assert!(h.reports()[0].steps_per_sec > 0.0);
    }
}
