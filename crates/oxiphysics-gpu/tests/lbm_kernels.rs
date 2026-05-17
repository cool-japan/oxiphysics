// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Integration tests for the LBM D3Q19 BGK GPU kernel.
//!
//! Test organisation:
//!
//! 1. `test_lbm_soa_to_gpu_buffer_roundtrip` — CPU-only marshalling correctness.
//!    Verifies the q-major flatten/unflatten round-trip.  Always passes without a GPU.
//!
//! 2. `test_lbm_d3q19_lid_cavity_gpu_vs_cpu` — GPU / CPU physics parity test.
//!    Uses **periodic-only** boundary conditions (no solid mask, no lid, no body
//!    force) so that both the WGSL shader and the CPU reference compute identical
//!    physics.  A non-uniform initial perturbation in density/velocity makes the
//!    simulation non-trivial; after 500 steps the mean velocities must agree
//!    within 1e-3.
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
    /// The solid mask is cleared (all fluid) so that the GPU WGSL shader
    /// (which uses periodic boundary conditions on all faces) and the CPU
    /// reference path compute identical physics.
    ///
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
        };
        let mut sim = LbmSimulation::new(cfg);
        // Clear all solid cells so physics is fully periodic.
        for v in sim.solid.iter_mut() {
            *v = false;
        }
        // Apply a tiny density perturbation so the system is non-trivial.
        sim.f[0][0] += 1e-3;
        sim
    }

    /// Run a reference CPU simulation in **periodic** mode:
    /// no solid mask, no lid, no body force.
    ///
    /// We re-implement the CPU step here as a minimal periodic BGK kernel that
    /// matches the WGSL shader exactly, to avoid the CPU `step_cpu()` applying
    /// bounce-back walls.
    fn run_cpu_periodic(
        nx: usize,
        ny: usize,
        nz: usize,
        tau: f64,
        n_steps: usize,
    ) -> (f64, f64, f64) {
        use oxiphysics_gpu::lbm_gpu::{D3Q19_EX, D3Q19_EY, D3Q19_EZ, D3Q19_W};

        let nc = nx * ny * nz;
        let omega = 1.0 / tau;

        // Initialise at rest with tiny density perturbation at cell 0.
        let mut f: Vec<Vec<f64>> = (0..19)
            .map(|q| {
                let mut dir = vec![D3Q19_W[q]; nc];
                if q == 0 {
                    dir[0] += 1e-3;
                }
                dir
            })
            .collect();
        let mut f_tmp: Vec<Vec<f64>> = vec![vec![0.0; nc]; 19];

        let nx_i = nx as i32;
        let ny_i = ny as i32;
        let nz_i = nz as i32;

        for _ in 0..n_steps {
            // BGK collision (in-place)
            for z in 0..nz {
                for y in 0..ny {
                    for x in 0..nx {
                        let cell = x + nx * (y + ny * z);
                        let mut rho = 0.0_f64;
                        let mut ux = 0.0_f64;
                        let mut uy = 0.0_f64;
                        let mut uz = 0.0_f64;
                        for i in 0..19 {
                            let fi = f[i][cell];
                            rho += fi;
                            ux += fi * D3Q19_EX[i] as f64;
                            uy += fi * D3Q19_EY[i] as f64;
                            uz += fi * D3Q19_EZ[i] as f64;
                        }
                        if rho < 1e-10 {
                            continue;
                        }
                        ux /= rho;
                        uy /= rho;
                        uz /= rho;
                        let u2 = ux * ux + uy * uy + uz * uz;
                        for i in 0..19 {
                            let eu = D3Q19_EX[i] as f64 * ux
                                + D3Q19_EY[i] as f64 * uy
                                + D3Q19_EZ[i] as f64 * uz;
                            let feq =
                                D3Q19_W[i] * rho * (1.0 + 3.0 * eu + 4.5 * eu * eu - 1.5 * u2);
                            f[i][cell] += omega * (feq - f[i][cell]);
                        }
                    }
                }
            }

            // Periodic streaming (pull-style, matching WGSL)
            for dir in f_tmp.iter_mut().take(19) {
                for v in dir.iter_mut().take(nc) {
                    *v = 0.0;
                }
            }
            for z in 0..nz {
                for y in 0..ny {
                    for x in 0..nx {
                        let cell = x + nx * (y + ny * z);
                        for i in 0..19 {
                            // Pull from upstream (periodic wrap)
                            let sx = ((x as i32 - D3Q19_EX[i] + nx_i).rem_euclid(nx_i)) as usize;
                            let sy = ((y as i32 - D3Q19_EY[i] + ny_i).rem_euclid(ny_i)) as usize;
                            let sz = ((z as i32 - D3Q19_EZ[i] + nz_i).rem_euclid(nz_i)) as usize;
                            f_tmp[i][cell] = f[i][sx + nx * (sy + ny * sz)];
                        }
                    }
                }
            }
            std::mem::swap(&mut f, &mut f_tmp);
        }

        // Compute mean velocity over all cells.
        // Build per-cell (rho, ux, uy, uz) by accumulating across directions.
        let mut cell_rho = vec![0.0_f64; nc];
        let mut cell_ux = vec![0.0_f64; nc];
        let mut cell_uy = vec![0.0_f64; nc];
        let mut cell_uz = vec![0.0_f64; nc];
        for (i, dir) in f.iter().enumerate().take(19) {
            let ex = D3Q19_EX[i] as f64;
            let ey = D3Q19_EY[i] as f64;
            let ez = D3Q19_EZ[i] as f64;
            for (fi, (rho_c, (ux_c, (uy_c, uz_c)))) in dir.iter().zip(
                cell_rho.iter_mut().zip(
                    cell_ux
                        .iter_mut()
                        .zip(cell_uy.iter_mut().zip(cell_uz.iter_mut())),
                ),
            ) {
                *rho_c += fi;
                *ux_c += fi * ex;
                *uy_c += fi * ey;
                *uz_c += fi * ez;
            }
        }
        let mut sum_ux = 0.0_f64;
        let mut sum_uy = 0.0_f64;
        let mut sum_uz = 0.0_f64;
        for (rho, (ux, (uy, uz))) in cell_rho
            .iter()
            .zip(cell_ux.iter().zip(cell_uy.iter().zip(cell_uz.iter())))
        {
            if *rho > 1e-10 {
                sum_ux += ux / rho;
                sum_uy += uy / rho;
                sum_uz += uz / rho;
            }
        }
        let n = nc as f64;
        (sum_ux / n, sum_uy / n, sum_uz / n)
    }

    /// LBM D3Q19 BGK GPU vs CPU parity test.
    ///
    /// Both paths use periodic boundary conditions (no solid mask, no lid, no
    /// body force) so they compute identical physics.  After 500 steps the
    /// mean velocities must agree within 1e-3 (GPU uses f32 internally).
    ///
    /// Skipped gracefully when no GPU adapter is available.
    #[ignore = "heavy: GPU physics parity test (>300s), run manually in isolated environment"]
    #[test]
    fn test_lbm_d3q19_lid_cavity_gpu_vs_cpu() {
        const NX: usize = 32;
        const NY: usize = 32;
        const NZ: usize = 32;
        const TAU: f64 = 1.0 / 1.8; // omega = 1.8
        const N_STEPS: usize = 500;

        // GPU path via LbmSimulation (periodic)
        let mut gpu_sim = periodic_sim(NX, NY, NZ, TAU);

        // Force lazy GPU init by taking one step
        gpu_sim.step();

        // Check whether a real GPU was acquired
        if !gpu_sim.has_gpu() {
            eprintln!("SKIPPED: no GPU adapter available for test_lbm_d3q19_lid_cavity_gpu_vs_cpu");
            return;
        }

        // Run remaining steps on GPU
        for _ in 1..N_STEPS {
            gpu_sim.step();
        }

        let (gpu_ux, gpu_uy, gpu_uz) = gpu_sim.mean_velocity();

        // CPU reference using the periodic mini-kernel
        let (cpu_ux, cpu_uy, cpu_uz) = run_cpu_periodic(NX, NY, NZ, TAU, N_STEPS);

        let diff_ux = (gpu_ux - cpu_ux).abs();
        let diff_uy = (gpu_uy - cpu_uy).abs();
        let diff_uz = (gpu_uz - cpu_uz).abs();

        assert!(
            diff_ux < 1e-3,
            "mean_ux mismatch: gpu={gpu_ux:.6} cpu={cpu_ux:.6} diff={diff_ux:.6e}"
        );
        assert!(
            diff_uy < 1e-3,
            "mean_uy mismatch: gpu={gpu_uy:.6} cpu={cpu_uy:.6} diff={diff_uy:.6e}"
        );
        assert!(
            diff_uz < 1e-3,
            "mean_uz mismatch: gpu={gpu_uz:.6} cpu={cpu_uz:.6} diff={diff_uz:.6e}"
        );

        // Mean density should still be approximately 1.0 (mass conserved)
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
