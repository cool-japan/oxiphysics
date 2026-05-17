// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Integration tests for the CUDA compute backend.
//!
//! These tests compile and run on any platform regardless of whether the
//! `cuda-backend` feature is enabled or CUDA hardware is present.  Tests that
//! require an RTX GPU are gated behind the `OXIPHYSICS_CUDA_BENCH=1` env var.

use std::time::Duration;

use oxiphysics_gpu::{
    compute::cuda_backend::CudaBackend,
    gpu_bench::{GpuBenchHarness, compute_cuda_speedup},
};

// ── Test 1: stub round-trip ───────────────────────────────────────────────────

/// Hardware-agnostic smoke test.
///
/// Verifies that:
/// - `CudaBackend::try_new(0)` does **not** panic.  It may return `Err` when
///   the `cuda-backend` feature is off or no CUDA driver is present, or `Ok`
///   when a real CUDA device is available.
/// - `cpu_vs_cuda_sph(512)` returns at least one report tagged "sph_density".
/// - The first report's `mean` is strictly positive.
/// - `compute_cuda_speedup` does not panic on the output.
#[test]
fn cuda_stub_roundtrip() {
    // try_new must not panic; its outcome depends on the build features and
    // the host hardware, so both `Ok` and `Err` are acceptable.
    let _ = CudaBackend::try_new(0);

    let mut h = GpuBenchHarness {
        warmup: 0,
        iterations: 1,
        reports: Vec::new(),
    };

    let reports = h.cpu_vs_cuda_sph(512);

    assert!(
        !reports.is_empty(),
        "cpu_vs_cuda_sph must return at least one report"
    );
    assert!(
        reports[0].name.contains("sph_density"),
        "report name should contain 'sph_density', got: {}",
        reports[0].name
    );
    assert!(
        reports[0].mean > Duration::ZERO,
        "CPU mean time must be positive, got {:?}",
        reports[0].mean
    );

    let sr = compute_cuda_speedup(&reports);
    assert!(
        sr.cpu_mean > Duration::ZERO,
        "cpu_mean should be positive in the speedup report"
    );
    // The CUDA report is only present when a real device is available.
    assert!(
        sr.cuda_mean.is_none() || sr.cuda_mean.unwrap() > Duration::ZERO,
        "cuda_mean, if present, must be positive"
    );
}

// ── Test 2: RTX hardware speedup (env-gated) ─────────────────────────────────

/// Hardware speedup assertion.
///
/// Only runs when `OXIPHYSICS_CUDA_BENCH=1` is set.
/// Asserts ≥ 5× speedup at N = 100 000 when a CUDA device is present.
///
/// Skip in CI unless you have an RTX-class runner:
/// ```text
/// OXIPHYSICS_CUDA_BENCH=1 cargo nextest run -p oxiphysics-gpu
///   --features cuda-backend cuda_sph_speedup_at_100k
/// ```
#[test]
fn cuda_sph_speedup_at_100k() {
    if std::env::var("OXIPHYSICS_CUDA_BENCH").as_deref() != Ok("1") {
        // Skip when the env gate is not set.
        return;
    }

    let mut h = GpuBenchHarness {
        warmup: 2,
        iterations: 5,
        reports: Vec::new(),
    };

    let reports = h.cpu_vs_cuda_sph(100_000);
    assert!(
        reports.len() >= 2,
        "expected CPU + CUDA reports, got {} (is CUDA hardware present?)",
        reports.len()
    );

    let sr = compute_cuda_speedup(&reports);
    let speedup = sr
        .speedup
        .expect("CUDA report present but no speedup computed");
    assert!(
        speedup >= 5.0,
        "expected ≥ 5× CUDA speedup over CPU at N=100k, got {speedup:.2}×"
    );
}
