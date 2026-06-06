// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Benchmark suite for the OxiPhysics materials constitutive kernels.
//!
//! Measures batch throughput (MFLOP/s or Mpoints/s) for each kernel in
//! [`crate::simd_paths`], covering:
//!
//! - Linear isotropic elasticity
//! - Von Mises plasticity
//! - Neo-Hookean hyperelasticity
//! - Viscoplasticity (Perzyna)
//! - Miner's rule fatigue damage
//! - Thermo-elastic coupling
//! - J2 return-mapping
//!
//! ## Usage
//!
//! ```
//! use oxiphysics_materials::materials_bench::{MatBenchHarness, bench_elastic};
//!
//! let mut h = MatBenchHarness::new();
//!
//! let r = h.run("elastic_n1000", || bench_elastic(1000));
//! println!("{}", r);
//! ```

use std::time::{Duration, Instant};

use crate::simd_paths::{
    SoaMaterialPoints, elastic_stress_batch, elastic_stress_batch_x4, miner_damage_batch,
    neo_hookean_stress_batch, return_mapping_batch, thermal_expansion_stress_batch,
    viscoplastic_rate_batch, von_mises_yield_batch,
};

// ── MatBenchReport ────────────────────────────────────────────────────────────

/// Result of one material batch kernel benchmark run.
#[derive(Debug, Clone)]
pub struct MatBenchReport {
    /// Human-readable kernel name.
    pub name: String,
    /// Number of timed iterations.
    pub iterations: u32,
    /// Number of material integration points.
    pub n: usize,
    /// Total wall-clock time for all iterations.
    pub total: Duration,
    /// Mean time per iteration.
    pub mean: Duration,
    /// Estimated throughput in MFLOP/s (kernel-specific estimate).
    pub mflops: Option<f64>,
    /// Throughput in Mpoints/s (1e6 points/second).
    pub mpoints_per_sec: f64,
}

impl std::fmt::Display for MatBenchReport {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "[{:<28}] n={:>7} mean={:.2}µs  {:.1} Mpts/s",
            self.name,
            self.n,
            self.mean.as_secs_f64() * 1e6,
            self.mpoints_per_sec
        )?;
        if let Some(mf) = self.mflops {
            write!(f, "  ({:.1} MFLOPs)", mf)?;
        }
        Ok(())
    }
}

// ── MatBenchHarness ───────────────────────────────────────────────────────────

/// Timing harness for material constitutive kernel benchmarks.
pub struct MatBenchHarness {
    /// Number of warm-up iterations (not timed).
    pub warmup: u32,
    /// Number of timed iterations.
    pub iterations: u32,
    /// Collected reports.
    pub reports: Vec<MatBenchReport>,
}

impl MatBenchHarness {
    /// Create a harness with 2 warm-up and 5 timed iterations.
    pub fn new() -> Self {
        Self {
            warmup: 2,
            iterations: 5,
            reports: Vec::new(),
        }
    }

    /// Create a faster harness (1 warm-up, 3 timed).
    pub fn fast() -> Self {
        Self {
            warmup: 1,
            iterations: 3,
            reports: Vec::new(),
        }
    }

    /// Run closure `f` and record a `MatBenchReport`.
    pub fn run<F>(&mut self, name: &str, mut f: F) -> MatBenchReport
    where
        F: FnMut() -> MatBenchResult,
    {
        let mut last = MatBenchResult::default();
        for _ in 0..self.warmup {
            last = f();
        }

        let t0 = Instant::now();
        for _ in 0..self.iterations {
            last = f();
        }
        let total = t0.elapsed();
        let mean = total / self.iterations;

        let sec = mean.as_secs_f64().max(1e-12);
        let mpoints = last.n as f64 / sec / 1e6;

        let r = MatBenchReport {
            name: name.to_string(),
            iterations: self.iterations,
            n: last.n,
            total,
            mean,
            mflops: last.mflops,
            mpoints_per_sec: mpoints,
        };
        self.reports.push(r.clone());
        r
    }

    /// Print a summary table of all collected reports.
    pub fn print_summary(&self) {
        println!("\n{:=<80}", "");
        println!(
            "{:<28} {:>8} {:>12} {:>12} {:>10}",
            "Kernel", "N", "Mean (µs)", "Mpts/s", "MFLOPs"
        );
        println!("{:=<80}", "");
        for r in &self.reports {
            let mf = r.mflops.map_or("—".to_string(), |m| format!("{:.1}", m));
            println!(
                "{:<28} {:>8} {:>12.3} {:>12.1} {:>10}",
                r.name,
                r.n,
                r.mean.as_secs_f64() * 1e6,
                r.mpoints_per_sec,
                mf
            );
        }
        println!("{:=<80}", "");
    }
}

impl Default for MatBenchHarness {
    fn default() -> Self {
        Self::new()
    }
}

// ── MatBenchResult ────────────────────────────────────────────────────────────

/// Intermediate output from a benchmark closure.
#[derive(Debug, Default, Clone)]
pub struct MatBenchResult {
    /// Problem size.
    pub n: usize,
    /// Optional MFLOP estimate.
    pub mflops: Option<f64>,
}

// ── Kernel benchmarks ─────────────────────────────────────────────────────────

/// Benchmark linear isotropic elasticity (scalar loop).
///
/// FLOPs ≈ 18 × n  (6 stress components, 3 additions each).
pub fn bench_elastic(n: usize) -> MatBenchResult {
    let mut pts = SoaMaterialPoints::new(n);
    let e = 210e9_f64;
    let nu = 0.3_f64;
    // Uniform uniaxial strain
    pts.strain_xx.fill(1e-3);
    elastic_stress_batch(&mut pts, e, nu);
    let flops = 18.0 * n as f64;
    MatBenchResult {
        n,
        mflops: Some(flops / 1e6),
    }
}

/// Benchmark linear isotropic elasticity (x4 unrolled loop).
///
/// Same physics as [`bench_elastic`], but using the unrolled 4-wide variant.
pub fn bench_elastic_x4(n: usize) -> MatBenchResult {
    let mut pts = SoaMaterialPoints::new(n);
    let e = 210e9_f64;
    let nu = 0.3_f64;
    pts.strain_xx.fill(1e-3);
    elastic_stress_batch_x4(&mut pts, e, nu);
    let flops = 18.0 * n as f64;
    MatBenchResult {
        n,
        mflops: Some(flops / 1e6),
    }
}

/// Benchmark von Mises yield surface evaluation.
///
/// FLOPs ≈ 10 × n  (deviatoric + Frobenius + sqrt + compare).
pub fn bench_von_mises(n: usize) -> MatBenchResult {
    let mut pts = SoaMaterialPoints::new(n);
    // Apply a stress state slightly above yield
    pts.stress_xx.fill(300e6);
    pts.stress_yy.fill(-50e6);
    let _f = von_mises_yield_batch(&pts);
    let flops = 10.0 * n as f64;
    MatBenchResult {
        n,
        mflops: Some(flops / 1e6),
    }
}

/// Benchmark Neo-Hookean hyperelastic Cauchy stress.
///
/// FLOPs ≈ 50 × n  (det F, left Cauchy-Green, stress assembly).
pub fn bench_neo_hookean(n: usize) -> MatBenchResult {
    let mut pts = SoaMaterialPoints::new(n);
    let e = 210e9_f64;
    let nu = 0.3_f64;
    // Perturb F slightly from identity
    for i in 0..n {
        pts.f11[i] = 1.001;
        pts.f22[i] = 0.999;
    }
    neo_hookean_stress_batch(&mut pts, e, nu);
    let flops = 50.0 * n as f64;
    MatBenchResult {
        n,
        mflops: Some(flops / 1e6),
    }
}

/// Benchmark Perzyna viscoplastic flow rate.
///
/// FLOPs ≈ 12 × n.
pub fn bench_viscoplastic(n: usize) -> MatBenchResult {
    let mut pts = SoaMaterialPoints::new(n);
    pts.stress_xx.fill(400e6); // overstressed
    // Compute yield values first
    let yield_vals = von_mises_yield_batch(&pts);
    let _rates = viscoplastic_rate_batch(&pts, &yield_vals, 1.0, 1e6);
    let flops = 12.0 * n as f64;
    MatBenchResult {
        n,
        mflops: Some(flops / 1e6),
    }
}

/// Benchmark Miner's rule fatigue damage accumulation (Basquin relation).
///
/// FLOPs ≈ 5 × n.
pub fn bench_fatigue(n: usize) -> MatBenchResult {
    let mut pts = SoaMaterialPoints::new(n);
    pts.stress_xx.fill(300e6);
    // miner_damage_batch(pts, sigma_f_prime, b_exponent, n_applied)
    miner_damage_batch(&mut pts, 500e6, -0.1, 1000.0);
    let flops = 5.0 * n as f64;
    MatBenchResult {
        n,
        mflops: Some(flops / 1e6),
    }
}

/// Benchmark thermo-elastic stress correction.
///
/// FLOPs ≈ 8 × n.
pub fn bench_thermal_expansion(n: usize) -> MatBenchResult {
    let mut pts = SoaMaterialPoints::new(n);
    for i in 0..n {
        pts.temperature[i] = 400.0;
    } // 400 K (hot)
    thermal_expansion_stress_batch(&mut pts, 210e9, 0.3, 12e-6, 293.15);
    let flops = 8.0 * n as f64;
    MatBenchResult {
        n,
        mflops: Some(flops / 1e6),
    }
}

/// Benchmark J2 return-mapping (elasto-plastic radial return).
///
/// FLOPs ≈ 30 × n.
pub fn bench_return_mapping(n: usize) -> MatBenchResult {
    let mut pts = SoaMaterialPoints::new(n);
    pts.stress_xx.fill(500e6); // overstressed; will be corrected to yield surface
    // return_mapping_batch(pts, e, nu, h) — h = isotropic hardening modulus
    return_mapping_batch(&mut pts, 210e9, 0.3, 1e9);
    let flops = 30.0 * n as f64;
    MatBenchResult {
        n,
        mflops: Some(flops / 1e6),
    }
}

// ── Full benchmark suite ──────────────────────────────────────────────────────

/// Run the complete materials benchmark suite and return a formatted summary.
///
/// ```
/// use oxiphysics_materials::materials_bench::run_mat_bench_suite;
/// let s = run_mat_bench_suite(false);
/// assert!(!s.is_empty());
/// ```
pub fn run_mat_bench_suite(verbose: bool) -> String {
    let mut h = MatBenchHarness::fast();

    let sizes = [1_000usize, 10_000, 100_000];

    for &n in &sizes {
        h.run(&format!("elastic_n{}", n), || bench_elastic(n));
        h.run(&format!("elastic_x4_n{}", n), || bench_elastic_x4(n));
        h.run(&format!("von_mises_n{}", n), || bench_von_mises(n));
        h.run(&format!("neo_hookean_n{}", n), || bench_neo_hookean(n));
        h.run(&format!("viscoplastic_n{}", n), || bench_viscoplastic(n));
        h.run(&format!("fatigue_n{}", n), || bench_fatigue(n));
        h.run(&format!("thermal_exp_n{}", n), || {
            bench_thermal_expansion(n)
        });
        h.run(&format!("return_map_n{}", n), || bench_return_mapping(n));
    }

    let mut out = String::new();
    if verbose {
        for r in &h.reports {
            out.push_str(&format!("{}\n", r));
        }
    } else {
        out = format!("{} material benchmarks completed", h.reports.len());
    }
    out
}

// ── tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_bench_elastic_runs() {
        let r = bench_elastic(100);
        assert_eq!(r.n, 100);
        assert!(r.mflops.is_some());
    }

    #[test]
    fn test_bench_elastic_x4_runs() {
        let r = bench_elastic_x4(100);
        assert_eq!(r.n, 100);
    }

    #[test]
    fn test_bench_von_mises_runs() {
        let r = bench_von_mises(100);
        assert_eq!(r.n, 100);
    }

    #[test]
    fn test_bench_neo_hookean_runs() {
        let r = bench_neo_hookean(100);
        assert_eq!(r.n, 100);
    }

    #[test]
    fn test_bench_viscoplastic_runs() {
        let r = bench_viscoplastic(100);
        assert_eq!(r.n, 100);
    }

    #[test]
    fn test_bench_fatigue_runs() {
        let r = bench_fatigue(100);
        assert_eq!(r.n, 100);
    }

    #[test]
    fn test_bench_thermal_expansion_runs() {
        let r = bench_thermal_expansion(100);
        assert_eq!(r.n, 100);
    }

    #[test]
    fn test_bench_return_mapping_runs() {
        let r = bench_return_mapping(100);
        assert_eq!(r.n, 100);
    }

    #[test]
    fn test_harness_collects_reports() {
        let mut h = MatBenchHarness {
            warmup: 0,
            iterations: 2,
            reports: Vec::new(),
        };
        h.run("el_100", || bench_elastic(100));
        h.run("vm_100", || bench_von_mises(100));
        assert_eq!(h.reports.len(), 2);
    }

    #[test]
    fn test_run_mat_bench_suite_compact() {
        let s = run_mat_bench_suite(false);
        assert!(s.contains("material benchmarks"), "output: {}", s);
    }

    #[test]
    fn test_run_mat_bench_suite_verbose() {
        let s = run_mat_bench_suite(true);
        assert!(s.contains("elastic"), "output: {}", s);
    }
}
