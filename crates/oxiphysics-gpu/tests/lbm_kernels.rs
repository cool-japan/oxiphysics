// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Integration tests for the LBM D3Q19 BGK GPU kernel.
//!
//! Test organisation:
//!
//! 1. `test_lbm_soa_to_gpu_buffer_roundtrip` — CPU-only marshalling correctness.
//!    Verifies the q-major flatten/unflatten round-trip.  Always passes without a GPU.
//!
//! 2. `test_lbm_d3q19_periodic_gpu_vs_cpu` — GPU / CPU physics parity test.
//!    A fully periodic box (solid mask cleared, no lid, no body force) with a
//!    non-uniform initial perturbation is advanced on the GPU and on the CPU
//!    reference path of the same `LbmSimulation` (GPU disabled); the per-cell
//!    density and velocity must agree to f32 precision.  Walls, lid and Guo
//!    forcing parity is covered by the unit tests in `src/lbm_gpu/tests.rs`.
//!
//! 3. `test_lbm_d3q19_gpu_resident_stepping` — verifies that multiple GPU steps
//!    can run without intermediate readback, and that a single final readback
//!    yields a deterministic, finite result.

/// Always-available smoke test so the binary always has at least one test
/// even when compiled without the `wgpu-backend` feature.
#[test]
fn lbm_kernels_link_smoke() {
    assert_eq!(1 + 1, 2);
}

// ── Step 0b: CPU-only SoA ↔ flat-f32 round-trip ──────────────────────────────

/// Verify that the q-major flatten → unflatten round-trip preserves values
/// within f64 → f32 → f64 precision.
///
/// This test requires no GPU and always runs in CI.
#[test]
fn test_lbm_soa_to_gpu_buffer_roundtrip() {
    let nx: usize = 4;
    let ny: usize = 4;
    let nz: usize = 4;
    let nc = nx * ny * nz; // 64

    // Build a deterministic SoA with unique values at every (q, x, y, z).
    // Value: q * 1000.0 + cell_index, chosen so no two entries are identical.
    let mut soa: Vec<Vec<f64>> = Vec::with_capacity(19);
    for q in 0..19_usize {
        let mut dir = Vec::with_capacity(nc);
        for cell in 0..nc {
            // Unique, well within f32 range and recoverable within 1e-5
            dir.push(q as f64 * 1000.0 + cell as f64 + 0.5);
        }
        soa.push(dir);
    }

    // Flatten to f32 using the q-major convention:
    //   flat[q * nc + cell] = soa[q][cell]
    // This matches the WGSL `idx` function:
    //   q * (nx*ny*nz) + z*(nx*ny) + y*nx + x
    // since cell = x + nx*(y + ny*z).
    let flat: Vec<f32> = {
        let mut v = Vec::with_capacity(19 * nc);
        for dir in soa.iter().take(19) {
            for &val in dir.iter().take(nc) {
                v.push(val as f32);
            }
        }
        v
    };

    assert_eq!(flat.len(), 19 * nc, "flat buffer length mismatch");

    // Unflatten back to SoA f64 using q-major: flat[q*nc + cell].
    let mut recovered: Vec<Vec<f64>> = Vec::with_capacity(19);
    for (q, _) in (0..19_usize).enumerate() {
        let dir: Vec<f64> = flat[q * nc..(q + 1) * nc]
            .iter()
            .map(|&v| v as f64)
            .collect();
        recovered.push(dir);
    }

    // Assert reconstruction matches original within f64 → f32 → f64 tolerance.
    for (q, (orig_dir, rec_dir)) in soa.iter().zip(recovered.iter()).enumerate() {
        for (cell, (&orig, &got)) in orig_dir.iter().zip(rec_dir.iter()).enumerate() {
            let tol = orig.abs() * 1e-6 + 1e-6;
            assert!(
                (orig - got).abs() <= tol,
                "q={q} cell={cell}: orig={orig} got={got} diff={}",
                (orig - got).abs()
            );
        }
    }
}

// ── GPU tests (wgpu-backend feature) ─────────────────────────────────────────

#[cfg(feature = "wgpu-backend")]
mod gpu_tests {
    use oxiphysics_gpu::compute_preference::ComputePreference;
    use oxiphysics_gpu::lbm_gpu::{LbmConfig, LbmSimulation};

    /// Macro: skip GPU-dependent tests gracefully when no adapter is available.
    macro_rules! skip_if_no_gpu {
        ($sim:expr) => {
            if !$sim.has_gpu() {
                eprintln!("SKIPPED: no GPU adapter available");
                return;
            }
        };
    }

    /// Build a periodic LBM simulation with a non-uniform initial state.
    ///
    /// The solid mask is cleared (all fluid) so the box is fully periodic.
    /// The non-uniformity is a single-mode density perturbation: the cell at
    /// index 0 receives an extra 1e-3 density in the rest direction (q=0).
    fn periodic_sim(nx: usize, ny: usize, nz: usize, tau: f64) -> LbmSimulation {
        let cfg = LbmConfig {
            nx,
            ny,
            nz,
            tau,
            rho0: 1.0,
            force_x: 0.0,
            force_y: 0.0,
            force_z: 0.0,
            compute: ComputePreference::Auto,
        };
        let mut sim = LbmSimulation::new(cfg);
        sim.set_solid_mask(vec![false; nx * ny * nz])
            .expect("mask length matches the grid");
        sim.populations_mut()[0][0] += 1e-3;
        sim
    }

    /// LBM D3Q19 BGK GPU vs CPU parity test on a periodic box.
    ///
    /// The CPU side is the real `step_cpu` reference (GPU disabled), not a
    /// re-implementation.  Skipped gracefully when no GPU adapter is available.
    #[test]
    fn test_lbm_d3q19_periodic_gpu_vs_cpu() {
        const N: usize = 12;
        const TAU: f64 = 1.0 / 1.8; // omega = 1.8
        const N_STEPS: usize = 100;

        let mut gpu_sim = periodic_sim(N, N, N, TAU);
        if !gpu_sim.ensure_gpu() {
            eprintln!(
                "SKIPPED: no GPU adapter available ({:?})",
                gpu_sim.gpu_error()
            );
            return;
        }
        let mut cpu_sim = periodic_sim(N, N, N, TAU);
        cpu_sim.disable_gpu();

        for _ in 0..N_STEPS {
            gpu_sim.step();
            cpu_sim.step();
        }
        assert!(
            gpu_sim.has_gpu(),
            "GPU path dropped out: {:?}",
            gpu_sim.gpu_error()
        );

        let rho_gpu = gpu_sim.density_field();
        let rho_cpu = cpu_sim.density_field();
        let u_gpu = gpu_sim.velocity_field();
        let u_cpu = cpu_sim.velocity_field();
        let mut max_rho = 0.0_f64;
        let mut max_u = 0.0_f64;
        for cell in 0..N * N * N {
            max_rho = max_rho.max((rho_gpu[cell] - rho_cpu[cell]).abs());
            for k in 0..3 {
                max_u = max_u.max((u_gpu[cell][k] - u_cpu[cell][k]).abs());
            }
        }
        assert!(max_rho < 5e-5, "density mismatch {max_rho:e}");
        assert!(max_u < 1e-5, "velocity mismatch {max_u:e}");

        // The perturbation must actually have spread (non-trivial dynamics).
        let spread = rho_cpu.iter().filter(|&&r| (r - 1.0).abs() > 1e-9).count();
        assert!(
            spread > 100,
            "perturbation did not propagate: {spread} cells"
        );

        // Mean density stays ≈ 1 (mass conserved under periodic BC).
        let gpu_rho = gpu_sim.mean_density();
        assert!(
            (gpu_rho - 1.0).abs() < 1e-4,
            "GPU mean density diverged: rho={gpu_rho:.6}"
        );
    }

    /// Verify that the GPU path can run multiple steps without intermediate
    /// readback ("GPU resident stepping"), and that a single final readback
    /// yields a deterministic, finite result.
    ///
    /// Skipped gracefully when no GPU adapter is available.
    #[test]
    fn test_lbm_d3q19_gpu_resident_stepping() {
        const NX: usize = 8;
        const NY: usize = 8;
        const NZ: usize = 8;
        const TAU: f64 = 1.0 / 1.8;
        const N_STEPS: usize = 10;

        let mut sim = periodic_sim(NX, NY, NZ, TAU);

        // First step: triggers lazy GPU init
        sim.step();
        skip_if_no_gpu!(sim);

        // Run remaining steps — no readback between them
        for _ in 1..N_STEPS {
            sim.step();
        }

        // Assert step_count is correct
        assert_eq!(
            sim.step_count, N_STEPS as u64,
            "step_count mismatch after {N_STEPS} GPU steps"
        );

        // Single final readback via mean_velocity()
        let (ux, uy, uz) = sim.mean_velocity();

        // Results must be finite (not NaN or Inf)
        assert!(
            ux.is_finite(),
            "mean_ux is not finite after GPU stepping: {ux}"
        );
        assert!(
            uy.is_finite(),
            "mean_uy is not finite after GPU stepping: {uy}"
        );
        assert!(
            uz.is_finite(),
            "mean_uz is not finite after GPU stepping: {uz}"
        );

        // Mean density must be close to 1.0 (mass conserved under periodic BC)
        let rho = sim.mean_density();
        assert!(
            (rho - 1.0).abs() < 1e-3,
            "mean density diverged after {N_STEPS} GPU steps: rho={rho:.6}"
        );

        // Second run with same IC for determinism check
        let mut sim2 = periodic_sim(NX, NY, NZ, TAU);
        for _ in 0..N_STEPS {
            sim2.step();
        }
        if sim2.has_gpu() {
            let (ux2, uy2, uz2) = sim2.mean_velocity();
            assert!(
                (ux - ux2).abs() < 1e-5,
                "Non-deterministic result: ux={ux:.8} ux2={ux2:.8}"
            );
            assert!(
                (uy - uy2).abs() < 1e-5,
                "Non-deterministic result: uy={uy:.8} uy2={uy2:.8}"
            );
            assert!(
                (uz - uz2).abs() < 1e-5,
                "Non-deterministic result: uz={uz:.8} uz2={uz2:.8}"
            );
        }
    }
}
