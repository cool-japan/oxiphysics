// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Unit tests for the D3Q19 LBM (CPU reference, GPU slot handling and
//! CPU-vs-GPU parity).

use super::*;

// ── helpers ───────────────────────────────────────────────────────────────────

fn cfg(nx: usize, ny: usize, nz: usize, tau: f64) -> LbmConfig {
    LbmConfig {
        nx,
        ny,
        nz,
        tau,
        ..LbmConfig::default()
    }
}

/// Simulation pinned to the CPU reference path.
fn cpu_sim(config: LbmConfig) -> LbmSimulation {
    let mut sim = LbmSimulation::new(config);
    sim.disable_gpu();
    sim
}

fn total_mass(sim: &mut LbmSimulation) -> f64 {
    sim.populations()
        .iter()
        .map(|dir| dir.iter().sum::<f64>())
        .sum()
}

/// Mean x-velocity over the fluid cells of the row `y` (all x, z).
fn row_mean_ux(sim: &mut LbmSimulation, y: usize) -> f64 {
    let (nx, nz) = (sim.config.nx, sim.config.nz);
    let mut sum = 0.0;
    let mut count = 0usize;
    for z in 0..nz {
        for x in 0..nx {
            let idx = x + nx * (y + sim.config.ny * z);
            if sim.solid()[idx] {
                continue;
            }
            sum += sim.cell_macro(x, y, z).1[0];
            count += 1;
        }
    }
    sum / count.max(1) as f64
}

/// Turn the default duct into a plane channel: clear the z-faces so only the
/// y = 0 and y = ny − 1 walls remain (periodic in x and z).
fn open_z_faces(sim: &mut LbmSimulation) {
    let (nx, ny, nz) = (sim.config.nx, sim.config.ny, sim.config.nz);
    for y in 1..ny - 1 {
        for x in 0..nx {
            sim.set_solid(x, y, 0, false).expect("in range");
            sim.set_solid(x, y, nz - 1, false).expect("in range");
        }
    }
}

// ── lattice constants ─────────────────────────────────────────────────────────

#[test]
fn test_d3q19_weights_sum() {
    let sum: f64 = D3Q19_W.iter().sum();
    assert!(
        (sum - 1.0).abs() < 1e-12,
        "weights should sum to 1, got {sum}"
    );
}

#[test]
fn test_d3q19_opposite_index() {
    for i in 0..19 {
        let j = D3Q19_OPP[i];
        assert_eq!(D3Q19_EX[i], -D3Q19_EX[j]);
        assert_eq!(D3Q19_EY[i], -D3Q19_EY[j]);
        assert_eq!(D3Q19_EZ[i], -D3Q19_EZ[j]);
    }
}

#[test]
fn test_lbm_construction() {
    let sim = LbmSimulation::new(cfg(4, 4, 4, 0.6));
    assert_eq!(sim.config.n_cells(), 64);
    assert_eq!(sim.f.len(), 19);
    assert_eq!(sim.f[0].len(), 64);
    // Default walls: y and z faces solid, x periodic.
    assert!(sim.solid()[0]);
    let interior = 1 + 4 * (1 + 4);
    assert!(!sim.solid()[interior]);
}

#[test]
fn test_lbm_viscosity() {
    let nu = cfg(16, 16, 16, 0.6).viscosity();
    // ν = (1/3)(τ − 0.5) = (1/3)(0.1) = 1/30
    assert!((nu - 1.0 / 30.0).abs() < 1e-10, "nu={nu}");
}

// ── WGSL kernel ───────────────────────────────────────────────────────────────

/// The kernel must parse and validate under naga, use a portable workgroup
/// size, and reference all four bindings (the auto-derived bind-group layout
/// drops unused bindings, which would make the 4-entry bind group invalid).
#[test]
fn test_lbm_wgsl_kernel_validates() {
    let module = match naga::front::wgsl::parse_str(LBM_BGK_D3Q19_WGSL) {
        Ok(m) => m,
        Err(e) => panic!(
            "LBM WGSL parse failed:\n{}",
            e.emit_to_string(LBM_BGK_D3Q19_WGSL)
        ),
    };
    let mut validator = naga::valid::Validator::new(
        naga::valid::ValidationFlags::all(),
        naga::valid::Capabilities::all(),
    );
    let info = match validator.validate(&module) {
        Ok(info) => info,
        Err(e) => panic!("LBM WGSL validation failed: {e:?}"),
    };
    let (ep_index, ep) = module
        .entry_points
        .iter()
        .enumerate()
        .find(|(_, ep)| ep.name == "main")
        .expect("entry point `main`");
    assert_eq!(ep.workgroup_size, [4, 4, 4]);
    #[cfg(feature = "wgpu-backend")]
    assert_eq!(ep.workgroup_size, [LBM_WORKGROUP_EDGE; 3]);
    let ep_info = info.get_entry_point(ep_index);
    let mut used_bindings: Vec<u32> = module
        .global_variables
        .iter()
        .filter(|(handle, _)| !ep_info[*handle].is_empty())
        .filter_map(|(_, var)| var.binding.as_ref().map(|b| b.binding))
        .collect();
    used_bindings.sort_unstable();
    assert_eq!(used_bindings, vec![0, 1, 2, 3]);
}

// ── CPU reference physics ─────────────────────────────────────────────────────

#[test]
fn test_lbm_mass_conservation_cpu_exact() {
    let mut sim = cpu_sim(cfg(6, 6, 6, 0.6));
    sim.set_lid_velocity(0.05, 0.0, 0.0);
    let before = total_mass(&mut sim);
    for _ in 0..50 {
        sim.step();
    }
    let after = total_mass(&mut sim);
    assert!(
        (before - after).abs() < 1e-12 * before,
        "mass drift: before={before} after={after}"
    );
}

#[test]
fn test_lbm_mass_conservation_default_path() {
    let mut sim = LbmSimulation::new(cfg(6, 6, 6, 0.6));
    let before = total_mass(&mut sim);
    for _ in 0..10 {
        sim.step();
    }
    let after = total_mass(&mut sim);
    assert!(
        (before - after).abs() / before < 1e-5,
        "mass should be conserved: before={before:.8} after={after:.8}"
    );
}

/// A zero lid velocity leaves a fluid at rest exactly at rest.
#[test]
fn test_lbm_rest_state_is_stationary() {
    let mut sim = cpu_sim(cfg(5, 5, 5, 0.7));
    let initial = sim.f.clone();
    for _ in 0..20 {
        sim.step();
    }
    for (a, b) in initial.iter().zip(sim.f.iter()) {
        for (x, y) in a.iter().zip(b.iter()) {
            assert!((x - y).abs() < 1e-15, "rest state drifted: {x} vs {y}");
        }
    }
}

/// Regression for the dead lid BC: the lid must drive the fluid below it.
#[test]
fn test_lbm_lid_driven_velocity() {
    let mut sim = cpu_sim(cfg(6, 8, 6, 0.6));
    sim.set_lid_velocity(0.05, 0.0, 0.0);
    for _ in 0..200 {
        sim.step();
    }
    let ny = sim.config.ny;
    let top = row_mean_ux(&mut sim, ny - 2);
    let bottom = row_mean_ux(&mut sim, 1);
    assert!(
        top > 1e-3,
        "row below lid should move with +x lid, got {top}"
    );
    assert!(top < 0.05, "fluid cannot outrun the lid: {top}");
    assert!(
        top > bottom,
        "velocity must decay away from the lid: top={top} bottom={bottom}"
    );
    // Profile is monotone in y (Couette-like, x periodic).
    let mut prev = f64::NEG_INFINITY;
    for y in 1..ny - 1 {
        let u = row_mean_ux(&mut sim, y);
        assert!(u > prev, "non-monotone lid profile at y={y}: {u} <= {prev}");
        prev = u;
    }
    let (ux, uy, uz) = sim.mean_velocity();
    assert!(ux > 0.0, "mean ux={ux}");
    assert!(uy.abs() < 1e-3 * ux.abs().max(1e-12) + 1e-12, "uy={uy}");
    assert!(uz.abs() < 1e-12, "uz={uz} must vanish by symmetry");
}

/// In a fully periodic box the Guo scheme adds exactly F of momentum per
/// cell and step; the physical mean velocity after k steps is (k + ½)F/ρ.
/// The previous velocity-shift-only forcing added ωF/2 per step instead.
#[test]
fn test_lbm_guo_forcing_momentum_exact() {
    let mut config = cfg(4, 4, 4, 0.8);
    config.force_x = 1e-5;
    config.force_z = -2e-6;
    let mut sim = cpu_sim(config);
    sim.set_solid_mask(vec![false; 64]).expect("mask length");
    let steps = 10;
    for _ in 0..steps {
        sim.step();
    }
    let (ux, uy, uz) = sim.mean_velocity();
    let k = steps as f64 + 0.5;
    assert!(
        (ux - k * 1e-5).abs() < 1e-14,
        "ux={ux} expected {}",
        k * 1e-5
    );
    assert!(uy.abs() < 1e-15, "uy={uy}");
    assert!(
        (uz + k * 2e-6).abs() < 1e-14,
        "uz={uz} expected {}",
        -k * 2e-6
    );
}

/// Guo forcing adds exactly `F` of momentum per cell and step: in a fully
/// periodic box the summed lattice momentum `Σ_cells Σ_i f_i e_i` equals
/// `steps · F · cells` (the previous velocity-shift-only forcing injected
/// `ωF/2` per step instead).
#[test]
fn test_lbm_guo_forcing_adds_exactly_f_momentum_per_cell_and_step() {
    let mut config = cfg(4, 4, 4, 0.8);
    config.force_x = 1e-5;
    let mut sim = cpu_sim(config);
    sim.set_solid_mask(vec![false; 64]).expect("mask length");
    let steps = 25;
    for _ in 0..steps {
        sim.step();
    }
    let mut momentum = [0.0_f64; 3];
    for (i, dir) in sim.populations().iter().enumerate() {
        let total: f64 = dir.iter().sum();
        momentum[0] += D3Q19_EX[i] as f64 * total;
        momentum[1] += D3Q19_EY[i] as f64 * total;
        momentum[2] += D3Q19_EZ[i] as f64 * total;
    }
    let expected = steps as f64 * 1e-5 * 64.0;
    assert!(
        (momentum[0] - expected).abs() < 1e-12,
        "x momentum {} expected {expected}",
        momentum[0]
    );
    assert!(momentum[1].abs() < 1e-12 && momentum[2].abs() < 1e-12);
}

#[test]
fn test_lbm_body_force_accelerates_flow() {
    let mut config = cfg(4, 6, 6, 0.55);
    config.force_x = 1e-5;
    let mut sim = cpu_sim(config);
    for _ in 0..100 {
        sim.step();
    }
    let (ux, _, _) = sim.mean_velocity();
    assert!(
        ux > 1e-4,
        "body force in +X should produce ux > 0, got {ux}"
    );
}

/// Plane Poiseuille flow between halfway bounce-back walls at y = 0.5 and
/// y = ny − 1.5.  At τ = ½ + √(3/16) BGK with halfway bounce-back places the
/// wall exactly, so the steady profile is the analytic parabola
/// u(y) = F/(2ρν) (y − ½)(ny − 3/2 − y).
#[test]
fn test_lbm_poiseuille_matches_analytic_profile() {
    let tau = 0.5 + (3.0_f64 / 16.0).sqrt();
    let force = 1e-5;
    let mut config = cfg(4, 10, 4, tau);
    config.force_x = force;
    let nu = config.viscosity();
    let mut sim = cpu_sim(config);
    open_z_faces(&mut sim);
    for _ in 0..3000 {
        sim.step();
    }
    let ny = sim.config.ny;
    let y_hi = ny as f64 - 1.5;
    let u_max = force / (2.0 * nu) * (y_hi - 0.5).powi(2) / 4.0;
    for y in 1..ny - 1 {
        let yf = y as f64;
        let expected = force / (2.0 * nu) * (yf - 0.5) * (y_hi - yf);
        let got = row_mean_ux(&mut sim, y);
        assert!(
            (got - expected).abs() < 1e-9 * u_max,
            "y={y}: got {got:.9e}, analytic {expected:.9e}"
        );
    }
}

// ── geometry setters ──────────────────────────────────────────────────────────

#[test]
fn test_lbm_set_solid_bounds_and_mask_length() {
    let mut sim = LbmSimulation::new(cfg(4, 4, 4, 0.6));
    assert!(sim.set_solid(4, 0, 0, true).is_err());
    assert!(sim.set_solid(0, 0, 4, true).is_err());
    assert!(sim.set_solid_mask(vec![false; 10]).is_err());
    let idx = 1 + 4 * (1 + 4);
    assert!(!sim.solid()[idx]);
    sim.set_solid(1, 1, 1, true).expect("in range");
    assert!(sim.solid()[idx]);
    sim.set_solid_mask(vec![false; 64]).expect("mask length");
    assert!(sim.solid().iter().all(|&s| !s));
}

#[test]
fn test_lbm_cell_macro_out_of_range_is_zero() {
    let mut sim = cpu_sim(cfg(4, 4, 4, 0.6));
    assert_eq!(sim.cell_macro(4, 0, 0), (0.0, [0.0; 3]));
    let (rho, u) = sim.cell_macro(1, 1, 1);
    assert!((rho - 1.0).abs() < 1e-12);
    assert_eq!(u, [0.0; 3]);
}

// ── GPU slot handling ─────────────────────────────────────────────────────────

#[cfg(feature = "wgpu-backend")]
#[test]
fn test_lbm_gpu_init_attempted_once() {
    let mut sim = LbmSimulation::new(cfg(4, 4, 4, 0.6));
    for _ in 0..100 {
        sim.step();
    }
    assert_eq!(sim.gpu_init_attempts, 1);
    assert_eq!(sim.step_count, 100);
    // Either the GPU is active, or its failure is remembered.
    assert!(sim.has_gpu() != sim.gpu_error().is_some());
}

#[cfg(feature = "wgpu-backend")]
#[test]
fn test_lbm_unavailable_gpu_is_not_retried() {
    let mut sim = LbmSimulation::new(cfg(4, 4, 4, 0.6));
    sim.gpu = GpuSlot::Unavailable(GpuError::BackendInit("no adapter (test)".to_string()));
    for _ in 0..100 {
        sim.step();
    }
    assert_eq!(
        sim.gpu_init_attempts, 0,
        "failed GPU init must not be retried"
    );
    assert_eq!(sim.step_count, 100);
    assert!(!sim.has_gpu());
    assert!(matches!(sim.gpu_error(), Some(GpuError::BackendInit(_))));
    assert!(!sim.ensure_gpu());
    assert_eq!(sim.gpu_init_attempts, 0);
}

#[cfg(feature = "wgpu-backend")]
#[test]
fn test_lbm_disabled_gpu_never_initialises() {
    let mut sim = cpu_sim(cfg(4, 4, 4, 0.6));
    for _ in 0..100 {
        sim.step();
    }
    assert_eq!(sim.gpu_init_attempts, 0);
    assert!(!sim.has_gpu());
    assert!(sim.gpu_error().is_none());
}

#[cfg(feature = "wgpu-backend")]
#[test]
fn test_lbm_unflatten_rejects_short_readback() {
    assert!(matches!(
        unflatten_f32_to_soa(&[], 8),
        Err(GpuError::ReadBack(_))
    ));
    let flat: Vec<f32> = (0..19 * 3).map(|v| v as f32).collect();
    let soa = unflatten_f32_to_soa(&flat, 3).expect("exact length");
    assert_eq!(soa.len(), 19);
    assert_eq!(soa[2], vec![6.0, 7.0, 8.0]);
    assert_eq!(flatten_soa_to_f32(&soa, 3), flat);
}

// ── CPU vs GPU parity (skipped without an adapter) ────────────────────────────

/// Build a (GPU, CPU) pair from the same configuration, or `None` when no GPU
/// adapter is available.
#[cfg(feature = "wgpu-backend")]
fn gpu_cpu_pair(config: LbmConfig) -> Option<(LbmSimulation, LbmSimulation)> {
    let mut gpu = LbmSimulation::new(config.clone());
    if !gpu.ensure_gpu() {
        eprintln!("SKIPPED: no GPU adapter ({:?})", gpu.gpu_error());
        return None;
    }
    Some((gpu, cpu_sim(config)))
}

/// Max absolute difference of the density and velocity fields.
///
/// Observed GPU(f32)-vs-CPU(f64) differences after 100–200 steps are
/// ≈ 5e-6 in ρ and ≈ 2e-7 in u; the tests allow ~10× that.
#[cfg(feature = "wgpu-backend")]
fn field_diff(a: &mut LbmSimulation, b: &mut LbmSimulation) -> (f64, f64) {
    let rho_diff = a
        .density_field()
        .iter()
        .zip(b.density_field().iter())
        .map(|(x, y)| (x - y).abs())
        .fold(0.0_f64, f64::max);
    let u_diff = a
        .velocity_field()
        .iter()
        .zip(b.velocity_field().iter())
        .flat_map(|(x, y)| (0..3).map(move |k| (x[k] - y[k]).abs()))
        .fold(0.0_f64, f64::max);
    (rho_diff, u_diff)
}

#[cfg(feature = "wgpu-backend")]
fn step_both(gpu: &mut LbmSimulation, cpu: &mut LbmSimulation, n: usize) {
    for _ in 0..n {
        gpu.step();
        cpu.step();
    }
    assert!(gpu.has_gpu(), "GPU path dropped out: {:?}", gpu.gpu_error());
}

#[cfg(feature = "wgpu-backend")]
#[test]
fn test_lbm_gpu_cpu_parity_poiseuille() {
    let mut config = cfg(8, 12, 8, 0.8);
    config.force_x = 1e-5;
    let Some((mut gpu, mut cpu)) = gpu_cpu_pair(config) else {
        return;
    };
    step_both(&mut gpu, &mut cpu, 200);
    let (rho_diff, u_diff) = field_diff(&mut gpu, &mut cpu);
    let (ux, _, _) = cpu.mean_velocity();
    assert!(ux > 1e-4, "channel flow must develop, ux={ux}");
    assert!(rho_diff < 5e-5, "density mismatch {rho_diff:e}");
    assert!(u_diff < 2e-6, "velocity mismatch {u_diff:e}");
}

#[cfg(feature = "wgpu-backend")]
#[test]
fn test_lbm_gpu_cpu_parity_lid_cavity() {
    let Some((mut gpu, mut cpu)) = gpu_cpu_pair(cfg(12, 12, 12, 0.6)) else {
        return;
    };
    gpu.set_lid_velocity(0.05, 0.0, 0.0);
    cpu.set_lid_velocity(0.05, 0.0, 0.0);
    step_both(&mut gpu, &mut cpu, 150);
    let (rho_diff, u_diff) = field_diff(&mut gpu, &mut cpu);
    let ny = gpu.config.ny;
    let top = row_mean_ux(&mut gpu, ny - 2);
    assert!(top > 1e-3, "GPU lid must drive the fluid, got {top}");
    assert!(rho_diff < 5e-5, "density mismatch {rho_diff:e}");
    assert!(u_diff < 2e-6, "velocity mismatch {u_diff:e}");
    // Tangential lid: mass is conserved on the GPU as well.
    let rho = gpu.mean_density();
    assert!((rho - 1.0).abs() < 1e-4, "GPU mean density {rho}");
}

/// Geometry, forcing, lid and population changes made between GPU steps must
/// reach the GPU (solid mask / params / populations re-upload).
#[cfg(feature = "wgpu-backend")]
#[test]
fn test_lbm_gpu_cpu_parity_mid_run_changes() {
    let Some((mut gpu, mut cpu)) = gpu_cpu_pair(cfg(10, 10, 10, 0.7)) else {
        return;
    };
    step_both(&mut gpu, &mut cpu, 20);

    // Obstacle + body force + lid switched on mid-run.
    for sim in [&mut gpu, &mut cpu] {
        for y in 3..6 {
            sim.set_solid(5, y, 5, true).expect("in range");
        }
        sim.config.force_x = 2e-5;
        sim.set_lid_velocity(0.03, 0.0, 0.01);
    }
    step_both(&mut gpu, &mut cpu, 40);

    // Perturb populations through the GPU-aware accessor.
    for sim in [&mut gpu, &mut cpu] {
        let cell = 2 + 10 * (4 + 10 * 4);
        sim.populations_mut()[1][cell] += 1e-3;
        sim.config.tau = 0.9;
    }
    step_both(&mut gpu, &mut cpu, 40);

    let (rho_diff, u_diff) = field_diff(&mut gpu, &mut cpu);
    assert!(rho_diff < 5e-5, "density mismatch {rho_diff:e}");
    assert!(u_diff < 2e-6, "velocity mismatch {u_diff:e}");
    let (ux, _, _) = gpu.mean_velocity();
    assert!(ux > 1e-4, "forced GPU flow must develop, ux={ux}");
}

#[cfg(feature = "wgpu-backend")]
#[test]
fn test_lbm_gpu_resident_stepping_deterministic() {
    let mut config = cfg(8, 8, 8, 0.6);
    config.force_x = 1e-5;
    let mut a = LbmSimulation::new(config.clone());
    let mut b = LbmSimulation::new(config);
    for _ in 0..10 {
        a.step();
        b.step();
    }
    if !a.has_gpu() || !b.has_gpu() {
        eprintln!("SKIPPED: no GPU adapter");
        return;
    }
    let fa = a.populations().to_vec();
    let fb = b.populations().to_vec();
    assert_eq!(fa, fb, "GPU stepping must be deterministic");
    assert!(fa.iter().flatten().all(|v| v.is_finite()));
}

/// Destroy GPU buffers of a ready simulation to provoke dispatch failures.
#[cfg(feature = "wgpu-backend")]
fn sabotage_gpu(sim: &mut LbmSimulation, also_input: bool) {
    let GpuSlot::Ready(state) = &mut sim.gpu else {
        panic!("GPU slot must be ready");
    };
    let solid = state.solid_buf;
    assert!(state.backend.destroy_buffer(solid));
    if also_input {
        let input = state.input_buf();
        assert!(state.backend.destroy_buffer(input));
    }
}

/// A dispatch failure mid-run hands the latest GPU populations to the CPU
/// reference (no `expect`, no lost state) and is remembered, not retried.
#[cfg(feature = "wgpu-backend")]
#[test]
fn test_lbm_gpu_dispatch_failure_falls_back_with_state() {
    let mut config = cfg(8, 8, 8, 0.7);
    config.force_x = 1e-5;
    let Some((mut gpu, mut cpu)) = gpu_cpu_pair(config) else {
        return;
    };
    gpu.set_lid_velocity(0.04, 0.0, 0.0);
    cpu.set_lid_velocity(0.04, 0.0, 0.0);
    step_both(&mut gpu, &mut cpu, 20);

    sabotage_gpu(&mut gpu, false);
    for _ in 0..20 {
        gpu.step();
        cpu.step();
    }
    assert!(!gpu.has_gpu());
    assert!(
        matches!(gpu.gpu_error(), Some(GpuError::ShaderDispatch(_))),
        "unexpected error state: {:?}",
        gpu.gpu_error()
    );
    assert_eq!(gpu.step_count, 40);
    assert_eq!(gpu.gpu_init_attempts, 1);
    let (rho_diff, u_diff) = field_diff(&mut gpu, &mut cpu);
    assert!(rho_diff < 5e-5, "density mismatch {rho_diff:e}");
    assert!(u_diff < 2e-6, "velocity mismatch {u_diff:e}");
}

/// When the GPU populations cannot be recovered `try_step` reports
/// `ReadBack` without stepping; `step` then continues on the CPU.
#[cfg(feature = "wgpu-backend")]
#[test]
fn test_lbm_gpu_unrecoverable_failure_reports_readback() {
    let Some((mut gpu, _cpu)) = gpu_cpu_pair(cfg(6, 6, 6, 0.7)) else {
        return;
    };
    for _ in 0..5 {
        gpu.step();
    }
    assert!(gpu.has_gpu());
    sabotage_gpu(&mut gpu, true);

    let result = gpu.try_step();
    assert!(
        matches!(result, Err(GpuError::ReadBack(_))),
        "expected ReadBack error, got {result:?}"
    );
    assert_eq!(gpu.step_count, 5, "no step may be taken on failure");
    assert!(matches!(gpu.gpu_error(), Some(GpuError::ReadBack(_))));
    assert!(!gpu.has_gpu());

    gpu.step();
    assert_eq!(gpu.step_count, 6);
    assert!(gpu.try_step().is_ok());
    assert_eq!(gpu.step_count, 7);
    assert!(gpu.populations().iter().flatten().all(|v| v.is_finite()));
}

// ── ComputePreference ─────────────────────────────────────────────────────────

/// Whether this host exposes a wgpu adapter at all.
fn adapter_present() -> bool {
    #[cfg(feature = "wgpu-backend")]
    {
        WgpuBackendReal::try_new().is_ok()
    }
    #[cfg(not(feature = "wgpu-backend"))]
    {
        false
    }
}

#[test]
fn test_lbm_config_defaults_to_auto_preference() {
    assert_eq!(LbmConfig::default().compute, ComputePreference::Auto);
}

/// `CpuOnly` in the configuration never initialises the GPU - not even on a
/// host that has an adapter - and runs the exact CPU reference.
#[test]
fn test_lbm_cpu_only_config_never_initialises_the_gpu() {
    let mut config = cfg(6, 6, 6, 0.7);
    config.force_x = 1e-5;
    config.compute = ComputePreference::CpuOnly;
    let mut sim = LbmSimulation::new(config.clone());
    assert!(!sim.ensure_gpu());
    sim.set_lid_velocity(0.03, 0.0, 0.0);
    let mut reference = cpu_sim(config);
    reference.set_lid_velocity(0.03, 0.0, 0.0);
    for _ in 0..30 {
        sim.step();
        reference.step();
    }
    assert!(!sim.has_gpu());
    assert!(sim.gpu_error().is_none());
    #[cfg(feature = "wgpu-backend")]
    assert_eq!(sim.gpu_init_attempts, 0);
    assert_eq!(sim.step_count, 30);
    assert_eq!(sim.f, reference.f, "CpuOnly must be the CPU reference");
}

/// `CpuOnly` selected on a live simulation switches the GPU off at the next
/// step, before any initialisation is attempted.
#[test]
fn test_lbm_cpu_only_selected_after_construction() {
    let mut sim = LbmSimulation::new(cfg(4, 4, 4, 0.6));
    sim.config.compute = ComputePreference::CpuOnly;
    for _ in 0..10 {
        sim.step();
    }
    assert!(!sim.has_gpu());
    assert!(sim.gpu_error().is_none());
    #[cfg(feature = "wgpu-backend")]
    assert_eq!(sim.gpu_init_attempts, 0);
    assert_eq!(sim.step_count, 10);
}

#[test]
fn test_lbm_disable_gpu_sets_the_cpu_only_preference() {
    let mut sim = LbmSimulation::new(cfg(4, 4, 4, 0.6));
    assert_eq!(sim.config.compute, ComputePreference::Auto);
    sim.disable_gpu();
    assert_eq!(sim.config.compute, ComputePreference::CpuOnly);
    assert!(!sim.has_gpu());
}

/// `GpuRequired` either takes the step on the GPU or reports a typed init
/// error and leaves the simulation untouched - never a silent CPU step.
#[test]
fn test_lbm_gpu_required_runs_on_the_gpu_or_reports_an_error() {
    let mut config = cfg(6, 6, 6, 0.7);
    config.compute = ComputePreference::GpuRequired;
    let mut sim = LbmSimulation::new(config);
    let before = sim.f.clone();
    match sim.try_step() {
        Ok(()) => {
            assert!(adapter_present());
            assert!(sim.has_gpu());
            assert_eq!(sim.step_count, 1);
        }
        Err(err) => {
            assert!(matches!(err, GpuError::BackendInit(_)), "{err:?}");
            assert_eq!(sim.step_count, 0);
            assert_eq!(sim.f, before, "no CPU step may be taken");
            sim.step();
            assert_eq!(
                sim.step_count, 0,
                "the infallible step must not fall back either"
            );
            assert_eq!(sim.f, before);
        }
    }
}

/// A required GPU that is disabled explicitly cannot take steps: the
/// contradiction is reported, not resolved by running on the CPU.
#[cfg(feature = "wgpu-backend")]
#[test]
fn test_lbm_gpu_required_with_disabled_gpu_is_an_error() {
    let mut sim = cpu_sim(cfg(4, 4, 4, 0.6));
    sim.config.compute = ComputePreference::GpuRequired;
    let before = sim.f.clone();
    let err = sim.try_step().expect_err("GPU required but disabled");
    assert!(matches!(err, GpuError::BackendInit(_)), "{err:?}");
    assert_eq!(sim.step_count, 0);
    sim.step();
    assert_eq!(sim.step_count, 0);
    assert_eq!(sim.f, before);
}

/// A dispatch failure of a required GPU is reported with its original error;
/// the populations are recovered but no CPU step replaces the failed one.
#[cfg(feature = "wgpu-backend")]
#[test]
fn test_lbm_gpu_required_dispatch_failure_is_not_hidden() {
    let mut config = cfg(8, 8, 8, 0.7);
    config.force_x = 1e-5;
    config.compute = ComputePreference::GpuRequired;
    let mut sim = LbmSimulation::new(config);
    if !sim.ensure_gpu() {
        eprintln!("SKIPPED: no GPU adapter ({:?})", sim.gpu_error());
        return;
    }
    for _ in 0..5 {
        sim.try_step().expect("healthy GPU step");
    }
    sabotage_gpu(&mut sim, false);
    let err = sim
        .try_step()
        .expect_err("the sabotaged dispatch must fail");
    assert!(matches!(err, GpuError::ShaderDispatch(_)), "{err:?}");
    assert_eq!(sim.step_count, 5, "the failed step must not be counted");
    assert!(!sim.has_gpu());
    // The recovered populations are the valid state after five GPU steps.
    assert!(sim.populations().iter().flatten().all(|v| v.is_finite()));
    let (ux, _, _) = sim.mean_velocity();
    assert!(
        ux > 1e-5,
        "the five healthy steps must be preserved, ux={ux}"
    );
    // Neither entry point falls back to the CPU afterwards.
    assert!(sim.try_step().is_err());
    sim.step();
    assert_eq!(sim.step_count, 5);
}

// ── LbmGpuSolver ──────────────────────────────────────────────────────────────

/// The CPU variant must observe the state it steps: with a body force the
/// velocity read back through the solver grows (on a GPU host the old variant
/// stepped on the GPU but read the stale CPU copy forever).
#[test]
fn test_lbm_gpu_solver_cpu_variant_reads_fresh_state() {
    let mut solver = LbmGpuSolver::new_cpu(6, 8, 6, 1.5);
    assert!(!solver.is_gpu());
    assert!(solver.gpu_error().is_none());
    solver.simulation_mut().config.force_x = 1e-4;
    let cell = 2 + 6 * (3 + 8 * 3);
    let rest = solver.read_velocity()[cell];
    for _ in 0..50 {
        solver.step().expect("CPU step");
    }
    let moved = solver.read_velocity()[cell];
    assert!(
        moved[0] > 1e-3,
        "interior velocity must grow under the body force, got {moved:?}"
    );
    assert!(moved[0] > 10.0 * rest[0].abs().max(1e-6));
    assert!(moved[1].abs() < 1e-6 && moved[2].abs() < 1e-6);
    let again = solver.read_velocity()[cell];
    assert_eq!(moved, again, "repeated reads observe the same state");
    assert!(!solver.is_gpu());
    assert_eq!(
        solver.simulation().config.compute,
        ComputePreference::CpuOnly
    );
}

/// `LbmGpuSolver::new` keeps the `Auto` preference while it runs on a GPU and
/// pins `CpuOnly` when no GPU can be initialised (the reason stays available),
/// so the configuration always states the device that actually runs.
#[test]
fn test_lbm_gpu_solver_new_pins_cpu_only_when_no_gpu_is_available() {
    let mut solver = LbmGpuSolver::new(4, 4, 4, 1.5);
    if solver.is_gpu() {
        assert_eq!(solver.simulation().config.compute, ComputePreference::Auto);
        assert!(solver.gpu_error().is_none());
    } else {
        assert_eq!(
            solver.simulation().config.compute,
            ComputePreference::CpuOnly
        );
        #[cfg(feature = "wgpu-backend")]
        assert!(
            solver.gpu_error().is_some(),
            "the initialisation failure must stay observable"
        );
        #[cfg(not(feature = "wgpu-backend"))]
        assert!(solver.gpu_error().is_none());
    }
    // Either way the solver steps and reports a consistent field.
    solver.step().expect("step");
    assert_eq!(solver.read_density().len(), 64);
}

#[test]
fn test_lbm_gpu_solver_construction() {
    let solver = LbmGpuSolver::new_cpu(8, 8, 8, 1.5);
    assert_eq!(solver.nx, 8);
    assert_eq!(solver.ny, 8);
    assert_eq!(solver.nz, 8);
    assert!((solver.omega - 1.5).abs() < 1e-6);
    assert!(!solver.is_gpu());
    assert!((solver.simulation().config.tau - 1.0 / 1.5).abs() < 1e-6);
}

#[test]
fn test_lbm_gpu_solver_density_init() {
    let mut solver = LbmGpuSolver::new_cpu(4, 4, 4, 1.5);
    let rho = solver.read_density();
    assert_eq!(rho.len(), 64);
    for &r in &rho {
        assert!((r - 1.0).abs() < 1e-4, "rho={r}");
    }
}

/// Mass conservation on the CPU variant under a moving lid and a body force.
///
/// A fluid at rest conserves mass trivially (and so would a solver that only
/// ever returns its stale initial state), so the flow must demonstrably develop
/// while the total mass stays constant.
#[test]
fn test_lbm_gpu_solver_step_conserves_mass_cpu() {
    let mut solver = LbmGpuSolver::new_cpu(8, 8, 8, 1.5);
    solver.set_lid_velocity(0.1, 0.0, 0.0);
    solver.simulation_mut().config.force_x = 1e-4;
    let total_mass = |rho: &[f32]| rho.iter().map(|&r| f64::from(r)).sum::<f64>();
    let peak_speed = |vel: &[[f32; 3]]| {
        vel.iter()
            .map(|u| f64::from(u[0] * u[0] + u[1] * u[1] + u[2] * u[2]).sqrt())
            .fold(0.0_f64, f64::max)
    };
    let mass_before = total_mass(&solver.read_density());
    // At rest apart from the Guo half-force correction F / 2.
    let speed_before = peak_speed(&solver.read_velocity());
    for _ in 0..10 {
        solver.step().expect("step failed");
    }
    assert_eq!(solver.step_count, 10);
    let speed_after = peak_speed(&solver.read_velocity());
    assert!(
        speed_after > 1e-3 && speed_after > 10.0 * speed_before,
        "the lid and the body force must set the fluid in motion: \
         peak speed {speed_before} -> {speed_after}"
    );
    let mass_after = total_mass(&solver.read_density());
    let rel_err = (mass_before - mass_after).abs() / mass_before;
    assert!(rel_err < 1e-5, "mass not conserved: rel_err={rel_err:e}");
}

/// Lid-driven cavity through the solver API (GPU when available): mass is
/// conserved and the lid actually drives the fluid on whichever path runs.
#[test]
fn test_lbm_gpu_lid_driven_cavity() {
    let (nx, ny, nz) = (16_u32, 16_u32, 16_u32);
    let mut solver = LbmGpuSolver::new(nx, ny, nz, 1.5);
    solver.set_lid_velocity(0.1, 0.0, 0.0);
    for _ in 0..100 {
        solver.step().expect("LBM step failed");
    }
    let rho = solver.read_density();
    let total_rho: f32 = rho.iter().sum();
    let n = (nx * ny * nz) as f32;
    let rel_err = (total_rho - n).abs() / n;
    assert!(rel_err < 1e-4, "density not conserved: rel_err={rel_err:e}");

    let vel = solver.read_velocity();
    let (nxu, nyu) = (nx as usize, ny as usize);
    let below_lid = 8 + nxu * ((nyu - 2) + nyu * 8);
    assert!(
        vel[below_lid][0] > 1e-3,
        "cell below lid should move in +x, got {:?}",
        vel[below_lid]
    );
}
