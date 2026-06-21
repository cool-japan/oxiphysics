use super::*;

// --- WasmSphConfig ---

#[test]
fn test_sph_config_default_valid() {
    let cfg = WasmSphConfig::default();
    assert!(cfg.validate().is_ok());
}

#[test]
fn test_sph_config_invalid_radius() {
    let cfg = WasmSphConfig {
        particle_radius: -0.01,
        ..Default::default()
    };
    assert!(cfg.validate().is_err());
}

#[test]
fn test_sph_config_smoothing_too_small() {
    let cfg = WasmSphConfig {
        particle_radius: 0.1,
        smoothing_length: 0.05,
        ..Default::default()
    };
    assert!(cfg.validate().is_err());
}

#[test]
fn test_sph_kernel_wendland_at_zero() {
    let cfg = WasmSphConfig::default();
    let w = cfg.kernel_wendland(0.0);
    assert!(w > 0.0);
}

#[test]
fn test_sph_kernel_wendland_beyond_h_is_zero_like() {
    let cfg = WasmSphConfig::default();
    let w = cfg.kernel_wendland(cfg.smoothing_length * 2.0);
    assert!(w >= 0.0);
}

#[test]
fn test_sph_kernel_cubic_decays() {
    let cfg = WasmSphConfig::default();
    let w0 = cfg.kernel_cubic(0.0);
    let w1 = cfg.kernel_cubic(cfg.smoothing_length * 0.5);
    assert!(w0 > w1);
}

// --- WasmSphParticle ---

#[test]
fn test_sph_particle_kinetic_energy() {
    let mut p = WasmSphParticle::new(0, [0.0; 3], 1.0);
    p.velocity = [1.0, 0.0, 0.0];
    assert!((p.kinetic_energy() - 0.5).abs() < 1e-10);
}

#[test]
fn test_sph_particle_speed() {
    let mut p = WasmSphParticle::new(1, [0.0; 3], 1.0);
    p.velocity = [3.0, 4.0, 0.0];
    assert!((p.speed() - 5.0).abs() < 1e-10);
}

#[test]
fn test_sph_particle_boundary() {
    let p = WasmSphParticle::boundary(0, [0.0; 3], 1.0);
    assert!(p.is_boundary);
}

// --- WasmSphSimulation ---

#[test]
fn test_sph_sim_add_particle() {
    let mut sim = WasmSphSimulation::new(WasmSphConfig::default());
    let idx = sim.add_particle([0.0, 1.0, 0.0], 0.1);
    assert_eq!(idx, 0);
    assert_eq!(sim.particles.len(), 1);
}

#[test]
fn test_sph_sim_add_particles_flat() {
    let mut sim = WasmSphSimulation::new(WasmSphConfig::default());
    let positions = vec![0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0];
    sim.add_particles(&positions, 0.1);
    assert_eq!(sim.particles.len(), 3);
}

#[test]
fn test_sph_sim_step_advances_time() {
    let mut sim = WasmSphSimulation::new(WasmSphConfig::default());
    sim.add_particle([0.0, 1.0, 0.0], 0.1);
    sim.step(0.01);
    assert!((sim.time - 0.01).abs() < 1e-15);
    assert_eq!(sim.step_count, 1);
}

#[test]
fn test_sph_sim_gravity_accelerates_particle() {
    let mut sim = WasmSphSimulation::new(WasmSphConfig::default());
    sim.add_particle([0.0, 1.0, 0.0], 0.1);
    sim.step(0.1);
    let vy = sim.particles[0].velocity[1];
    assert!(vy < 0.0, "gravity should accelerate particle downward");
}

#[test]
fn test_sph_sim_get_positions_length() {
    let mut sim = WasmSphSimulation::new(WasmSphConfig::default());
    sim.add_particle([0.0; 3], 0.1);
    sim.add_particle([1.0, 0.0, 0.0], 0.1);
    let pos = sim.get_positions();
    assert_eq!(pos.len(), 6);
}

#[test]
fn test_sph_sim_fluid_count_excludes_boundary() {
    let mut sim = WasmSphSimulation::new(WasmSphConfig::default());
    sim.add_particle([0.0; 3], 0.1);
    sim.particles
        .push(WasmSphParticle::boundary(99, [1.0; 3], 0.1));
    assert_eq!(sim.fluid_count(), 1);
}

#[test]
fn test_sph_sim_total_kinetic_energy() {
    let mut sim = WasmSphSimulation::new(WasmSphConfig::default());
    let idx = sim.add_particle([0.0; 3], 1.0);
    sim.particles[idx].velocity = [1.0, 0.0, 0.0];
    let ke = sim.total_kinetic_energy();
    assert!((ke - 0.5).abs() < 1e-10);
}

// --- L1: SPH density summation tests ---

#[test]
fn test_sph_density_single_particle_self_contribution() {
    use std::f64::consts::PI;
    let h = 0.2_f64;
    let m = 0.5_f64;
    let expected = m * 8.0 / (PI * h.powi(3));
    let densities = compute_sph_densities(&[WasmSphParticle::new(0, [0.0; 3], m)], h);
    assert!(
        (densities[0] - expected).abs() < 1e-12,
        "single-particle density {:.6} != expected {:.6}",
        densities[0],
        expected
    );
}

#[test]
fn test_sph_density_cubic_lattice_symmetry() {
    let h = 0.2_f64;
    let m = 0.1_f64;
    let s = h / 2.0;
    let positions: Vec<[f64; 3]> = vec![
        [0.0, 0.0, 0.0],
        [s, 0.0, 0.0],
        [0.0, s, 0.0],
        [s, s, 0.0],
        [0.0, 0.0, s],
        [s, 0.0, s],
        [0.0, s, s],
        [s, s, s],
    ];
    let particles: Vec<WasmSphParticle> = positions
        .iter()
        .enumerate()
        .map(|(i, &pos)| WasmSphParticle::new(i as u64, pos, m))
        .collect();
    let densities = compute_sph_densities(&particles, h);
    let rho0 = densities[0];
    for (i, &rho) in densities.iter().enumerate() {
        assert!(
            (rho - rho0).abs() < 1e-12,
            "particle {i} density {rho:.6} differs from particle 0 density {rho0:.6}"
        );
    }
    let expected = m
        * (cubic_spline_compact(0.0, h)
            + 3.0 * cubic_spline_compact(s, h)
            + 3.0 * cubic_spline_compact(s * std::f64::consts::SQRT_2, h)
            + cubic_spline_compact(s * 3.0_f64.sqrt(), h));
    assert!(
        (rho0 - expected).abs() < 1e-10,
        "cubic-lattice density {rho0:.6} != expected {expected:.6}"
    );
    assert!(rho0 > 0.0 && rho0.is_finite());
}

#[test]
fn test_sph_density_step_not_stub() {
    let mut sim = WasmSphSimulation::new(WasmSphConfig::default());
    let h = sim.config.smoothing_length;
    sim.add_particle([0.0, 0.0, 0.0], 0.1);
    sim.add_particle([h * 0.3, 0.0, 0.0], 0.1);
    sim.step(0.001);
    let rho0 = sim.particles[0].density;
    let rho1 = sim.particles[1].density;
    assert!(
        rho0 > 0.0 && rho0 < sim.config.rest_density * 0.5,
        "density {rho0:.4} should be well below rest_density after real SPH summation"
    );
    assert!(
        (rho0 - rho1).abs() < 1e-8,
        "symmetric pair should have equal densities: {rho0:.6} vs {rho1:.6}"
    );
}

// --- WasmLbmConfig ---

#[test]
fn test_lbm_config_d2q9_valid() {
    let cfg = WasmLbmConfig::d2q9(64, 64, 0.1);
    assert!(cfg.validate().is_ok());
}

#[test]
fn test_lbm_config_d3q19_valid() {
    let cfg = WasmLbmConfig::d3q19(32, 32, 32, 0.1);
    assert!(cfg.validate().is_ok());
}

#[test]
fn test_lbm_config_invalid_grid() {
    let cfg = WasmLbmConfig {
        nx: 0,
        ..Default::default()
    };
    assert!(cfg.validate().is_err());
}

#[test]
fn test_lbm_config_invalid_tau() {
    let cfg = WasmLbmConfig {
        tau: 0.4,
        ..Default::default()
    };
    assert!(cfg.validate().is_err());
}

#[test]
fn test_lbm_config_reynolds_number() {
    let cfg = WasmLbmConfig::d2q9(64, 64, 0.1);
    let re = cfg.reynolds_number(0.1, 10.0);
    assert!((re - 10.0).abs() < 1e-10);
}

// --- WasmLbmSimulation ---

#[test]
fn test_lbm_sim_node_count() {
    let cfg = WasmLbmConfig::d2q9(8, 16, 0.1);
    let sim = WasmLbmSimulation::new(cfg);
    assert_eq!(sim.node_count(), 128);
}

#[test]
fn test_lbm_sim_step_increments_count() {
    let mut sim = WasmLbmSimulation::new(WasmLbmConfig::default());
    sim.step();
    assert_eq!(sim.step_count, 1);
}

#[test]
fn test_lbm_sim_set_boundary() {
    let mut sim = WasmLbmSimulation::new(WasmLbmConfig::d2q9(8, 8, 0.1));
    sim.set_boundary(0, 0, 0, true);
    assert!(sim.boundary[0]);
}

#[test]
fn test_lbm_sim_pressure_field_length() {
    let sim = WasmLbmSimulation::new(WasmLbmConfig::d2q9(4, 4, 0.1));
    let p = sim.get_pressure_field();
    assert_eq!(p.len(), 16);
}

// --- L2: LBM D3Q19 BGK facade tests ---

#[test]
fn test_lbm_d3q19_step_no_panic() {
    let mut sim = WasmLbmSimulation::new(WasmLbmConfig::d3q19(4, 4, 4, 0.1));
    sim.step();
    assert_eq!(sim.step_count, 1);
}

#[test]
fn test_lbm_d3q19_density_nonzero_after_init() {
    let sim = WasmLbmSimulation::new(WasmLbmConfig::d3q19(4, 4, 4, 0.1));
    let mean = sim.mean_density();
    assert!(
        (mean - 1.0).abs() < 1e-10,
        "initial mean density should be 1.0, got {mean}"
    );
    assert!(sim.density.iter().all(|&rho| rho > 0.0));
}

#[test]
fn test_lbm_d3q19_density_conservation() {
    let mut sim = WasmLbmSimulation::new(WasmLbmConfig::d3q19(8, 8, 4, 0.6));
    let total_init: f64 = sim.density.iter().sum();
    for _ in 0..10 {
        sim.step();
    }
    let total_final: f64 = sim.density.iter().sum();
    let rel_err = (total_final - total_init).abs() / total_init;
    assert!(
        rel_err < 1e-10,
        "total density changed by relative error {rel_err:.2e} (should be < 1e-10)"
    );
}

// --- WasmFluidStats ---

#[test]
fn test_fluid_stats_from_sph_empty() {
    let sim = WasmSphSimulation::new(WasmSphConfig::default());
    let stats = WasmFluidStats::from_sph(&sim, 0.01, 0.05);
    assert_eq!(stats.active_count, 0);
}

#[test]
fn test_fluid_stats_cfl_ok() {
    let mut stats = WasmFluidStats::new();
    stats.cfl = 0.5;
    assert!(stats.is_cfl_ok());
    stats.cfl = 1.1;
    assert!(!stats.is_cfl_ok());
}

// --- WasmMultiphaseConfig ---

#[test]
fn test_multiphase_density_ratio() {
    let cfg = WasmMultiphaseConfig::default();
    let ratio = cfg.density_ratio();
    assert!(ratio > 100.0);
}

#[test]
fn test_multiphase_capillary_length() {
    let cfg = WasmMultiphaseConfig::default();
    let l = cfg.capillary_length();
    assert!(l > 0.001 && l < 0.01);
}

// --- WasmFluidCoupling ---

#[test]
fn test_coupling_buoyancy() {
    let c = WasmFluidCoupling::sphere_in_water();
    let f = c.buoyancy(0.001);
    assert!((f - 9.81).abs() < 0.01);
}

#[test]
fn test_coupling_stokes_drag() {
    let c = WasmFluidCoupling::sphere_in_water();
    let f = c.stokes_drag(0.01, 0.1);
    assert!(f > 0.0);
}

#[test]
fn test_coupling_reynolds() {
    let c = WasmFluidCoupling::sphere_in_water();
    let re = c.reynolds(1.0, 0.01);
    assert!((re - 10000.0).abs() < 1.0);
}

// --- WasmParticleEmitter ---

#[test]
fn test_emitter_emit_zero_dt() {
    let mut em = WasmParticleEmitter::new([0.0; 3], [0.0, 1.0, 0.0], 100.0, 1.0, 5.0, 0.01);
    let particles = em.emit(0.0);
    assert_eq!(particles.len(), 0);
}

#[test]
fn test_emitter_emit_one_second() {
    let mut em = WasmParticleEmitter::new([0.0; 3], [0.0, 1.0, 0.0], 10.0, 1.0, 5.0, 0.01);
    let particles = em.emit(1.0);
    assert_eq!(particles.len(), 10);
}

#[test]
fn test_emitter_total_emitted() {
    let mut em = WasmParticleEmitter::new([0.0; 3], [1.0, 0.0, 0.0], 5.0, 2.0, 10.0, 0.1);
    em.emit(1.0);
    assert_eq!(em.total_emitted(), 5);
}

#[test]
fn test_emitter_reset() {
    let mut em = WasmParticleEmitter::new([0.0; 3], [1.0, 0.0, 0.0], 10.0, 1.0, 5.0, 0.01);
    em.emit(1.0);
    em.reset();
    assert_eq!(em.total_emitted(), 0);
}

// --- WasmFlowAnalyzer ---

#[test]
fn test_flow_analyzer_add_seeds() {
    let mut fa = WasmFlowAnalyzer::new([0.0, 1.0, 0.0, 1.0, 0.0, 1.0]);
    fa.add_seed(StreamlineSeed::new([0.5, 0.5, 0.5], 1.0, 0.1));
    assert_eq!(fa.seeds.len(), 1);
}

#[test]
fn test_flow_analyzer_integrate_streamlines() {
    let mut fa = WasmFlowAnalyzer::new([0.0, 10.0, 0.0, 10.0, 0.0, 10.0]);
    fa.add_seed(StreamlineSeed::new([1.0, 5.0, 5.0], 2.0, 0.5));
    fa.integrate_streamlines([1.0, 0.0, 0.0]);
    assert_eq!(fa.streamline_count(), 1);
    assert!(fa.streamlines[0].len() >= 2);
}

#[test]
fn test_flow_analyzer_ftle_requires_velocity_field() {
    // Without a configured velocity field the FTLE must report an honest
    // error rather than fabricating a field.
    let mut fa = WasmFlowAnalyzer::new([0.0, 1.0, 0.0, 1.0, 0.0, 1.0]);
    assert!(fa.compute_ftle(4, 4, 1, 1.0).is_err());
}

#[test]
fn test_flow_analyzer_ftle_length() {
    let mut fa = WasmFlowAnalyzer::new([0.0, 1.0, 0.0, 1.0, 0.0, 1.0]);
    fa.set_uniform_velocity_field([1.0, 0.0, 0.0]);
    let ftle = fa
        .compute_ftle(4, 4, 1, 1.0)
        .expect("uniform field FTLE should succeed");
    assert_eq!(ftle.len(), 16);
}

#[test]
fn test_flow_analyzer_ftle_uniform_flow_is_zero() {
    // A spatially-uniform translation flow has flow-map Jacobian = I, so
    // λ_max(∇Fᵀ∇F) = 1 and the FTLE is identically zero everywhere — NOT
    // the old |sin·cos| pattern.
    let mut fa = WasmFlowAnalyzer::new([0.0, 1.0, 0.0, 1.0, 0.0, 0.0]);
    fa.set_uniform_velocity_field([0.7, -0.3, 0.0]);
    let ftle = fa
        .compute_ftle(5, 5, 1, 2.0)
        .expect("uniform field FTLE should succeed")
        .to_vec();
    assert_eq!(ftle.len(), 25);
    for (i, &v) in ftle.iter().enumerate() {
        assert!(
            v.abs() < 1e-6,
            "uniform-flow FTLE at node {i} should be ~0, got {v}"
        );
    }
}

#[test]
fn test_flow_analyzer_ftle_saddle_matches_strain_rate() {
    // Linear saddle field u = (a·x, −a·y, 0). The analytic flow map is
    // F(x,y) = (x·e^{a t}, y·e^{−a t}), giving λ_max(∇Fᵀ∇F) = e^{2 a t}
    // and FTLE = (1/t)·½·ln(e^{2 a t}) = a everywhere.
    let a = 0.8_f64;
    let t = 1.5_f64;
    let mut fa = WasmFlowAnalyzer::new([-1.0, 1.0, -1.0, 1.0, 0.0, 0.0]);
    fa.set_linear_velocity_field([[a, 0.0, 0.0], [0.0, -a, 0.0], [0.0, 0.0, 0.0]], [0.0; 3]);
    let ftle = fa
        .compute_ftle(6, 6, 1, t)
        .expect("saddle field FTLE should succeed")
        .to_vec();
    assert_eq!(ftle.len(), 36);
    for (i, &v) in ftle.iter().enumerate() {
        assert!(
            (v - a).abs() < 1e-3,
            "saddle FTLE at node {i} should be a={a}, got {v}"
        );
    }
    // Sanity: this is decidedly not the discarded |sin·cos| pattern, whose
    // values vary node-to-node between 0 and ~1.
    let max = ftle.iter().cloned().fold(f64::MIN, f64::max);
    let min = ftle.iter().cloned().fold(f64::MAX, f64::min);
    assert!(
        (max - min).abs() < 1e-3,
        "saddle FTLE should be spatially uniform; spread = {}",
        max - min
    );
}

#[test]
fn test_flow_analyzer_ftle_rotation_is_zero() {
    // Solid-body rotation is volume- and length-preserving (∇F is a pure
    // rotation), so the FTLE is zero everywhere.
    let omega = 1.2_f64;
    let mut fa = WasmFlowAnalyzer::new([-1.0, 1.0, -1.0, 1.0, 0.0, 0.0]);
    fa.set_linear_velocity_field(
        [[0.0, -omega, 0.0], [omega, 0.0, 0.0], [0.0, 0.0, 0.0]],
        [0.0; 3],
    );
    let ftle = fa
        .compute_ftle(5, 5, 1, 1.0)
        .expect("rotation field FTLE should succeed")
        .to_vec();
    for (i, &v) in ftle.iter().enumerate() {
        assert!(
            v.abs() < 1e-3,
            "rotation FTLE at node {i} should be ~0, got {v}"
        );
    }
}

#[test]
fn test_flow_analyzer_ftle_invalid_time_errors() {
    let mut fa = WasmFlowAnalyzer::new([0.0, 1.0, 0.0, 1.0, 0.0, 0.0]);
    fa.set_uniform_velocity_field([1.0, 0.0, 0.0]);
    assert!(fa.compute_ftle(4, 4, 1, 0.0).is_err());
    assert!(fa.compute_ftle(4, 4, 1, f64::NAN).is_err());
    assert!(fa.compute_ftle(0, 0, 0, 1.0).is_err());
}

#[test]
fn test_flow_analyzer_clear() {
    let mut fa = WasmFlowAnalyzer::new([0.0, 1.0, 0.0, 1.0, 0.0, 1.0]);
    fa.add_seed(StreamlineSeed::new([0.5; 3], 1.0, 0.1));
    fa.integrate_streamlines([1.0, 0.0, 0.0]);
    fa.clear_results();
    assert_eq!(fa.streamline_count(), 0);
}

#[test]
fn test_flow_analyzer_total_points() {
    let mut fa = WasmFlowAnalyzer::new([0.0, 100.0, 0.0, 100.0, 0.0, 100.0]);
    fa.add_seed(StreamlineSeed::new([0.0, 50.0, 50.0], 10.0, 1.0));
    fa.add_seed(StreamlineSeed::new([0.0, 25.0, 50.0], 5.0, 1.0));
    fa.integrate_streamlines([1.0, 0.0, 0.0]);
    assert!(fa.total_points() > 0);
}

#[test]
fn test_strouhal_viv() {
    let st = WasmFluidCoupling::strouhal_viv(0.2, 0.01, 1.0);
    assert!((st - 0.002).abs() < 1e-10);
}

// --- H4: RK4 circular streamline ---

#[test]
fn rk4_circular_streamline() {
    let fa = WasmFlowAnalyzer::new([-10.0, 10.0, -10.0, 10.0, -10.0, 10.0]);
    let h = 0.01_f64;
    let steps = 100_usize;
    let mut pos = [1.0_f64, 0.0, 0.0];
    let vel_fn = |p: [f64; 3]| -> [f64; 3] { [-p[1], p[0], 0.0] };
    for _ in 0..steps {
        pos = fa.rk4_step(pos, h, &vel_fn);
    }
    let radius = (pos[0] * pos[0] + pos[1] * pos[1]).sqrt();
    assert!(
        (radius - 1.0).abs() < 0.01,
        "RK4 radius error too large: radius = {radius}"
    );
}
