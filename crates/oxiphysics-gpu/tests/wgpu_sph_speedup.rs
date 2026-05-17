// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

use oxiphysics_gpu::gpu_bench::{GpuBenchHarness, compute_speedup};

/// Smoke test: `cpu_vs_wgpu_sph` returns at least one report and CPU time is finite.
/// Does NOT assert any speedup ratio — hardware-agnostic.
#[test]
fn wgpu_sph_speedup_smoke() {
    let mut h = GpuBenchHarness {
        warmup: 0,
        iterations: 1,
        reports: Vec::new(),
    };
    let reports = h.cpu_vs_wgpu_sph(512);

    assert!(!reports.is_empty(), "should return at least one report");
    assert!(
        reports[0].name.contains("sph"),
        "report name should mention sph, got: {}",
        reports[0].name
    );
    assert!(
        reports[0].mean.as_secs_f64() > 0.0,
        "CPU mean time should be positive, got {:?}",
        reports[0].mean
    );

    let speedup = compute_speedup(&reports);
    assert!(
        speedup.cpu_mean.as_secs_f64() > 0.0,
        "cpu_mean should be positive"
    );

    // Print observed speedup for informational purposes
    if let Some(s) = speedup.speedup {
        eprintln!("Observed SPH speedup (wgpu/CPU) at n=512: {:.2}x", s);
    } else {
        eprintln!("No wgpu adapter available; only CPU path measured");
    }
}

/// RTX-class assertion: wgpu SPH should be >= 5x faster than CPU at N=100_000.
///
/// This test is ENV-GATED: set `OXIPHYSICS_RTX_BENCH=1` to activate.
/// Without the env var the test exits immediately to keep CI fast.
/// The full bench and assertion only run on RTX-class hardware.
#[test]
fn wgpu_sph_speedup_at_100k() {
    let rtx_mode = std::env::var("OXIPHYSICS_RTX_BENCH").is_ok();
    if !rtx_mode {
        eprintln!("OXIPHYSICS_RTX_BENCH not set; skipping N=100_000 benchmark");
        return;
    }

    let mut h = GpuBenchHarness {
        warmup: 1,
        iterations: 3,
        reports: Vec::new(),
    };
    let reports = h.cpu_vs_wgpu_sph(100_000);

    let speedup = compute_speedup(&reports);

    if let Some(s) = speedup.speedup {
        eprintln!("SPH speedup at N=100_000: {:.2}x", s);
        assert!(
            s >= 5.0,
            "Expected wgpu SPH >= 5x CPU at N=100_000 (OXIPHYSICS_RTX_BENCH mode), got {:.2}x",
            s
        );
    } else {
        panic!("OXIPHYSICS_RTX_BENCH is set but no wgpu adapter was found");
    }
}
