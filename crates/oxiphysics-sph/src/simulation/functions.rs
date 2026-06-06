//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[cfg(test)]
use crate::simulation::types::*;
#[cfg(test)]
use crate::simulation::types_sim::*;
#[cfg(test)]
use oxiphysics_core::math::Vec3;

#[cfg(test)]
mod tests {
    use super::*;
    use crate::boundary_sph::BoundaryPlane;
    use crate::kernel::CubicSplineKernel;
    use crate::particle::SphParticle;

    #[test]
    fn dam_break_particles_fall_and_spread() {
        let params = SphSimulationParams {
            gravity: Vec3::new(0.0, -9.81, 0.0),
            smoothing_length: 0.08,
            sound_speed: 50.0,
            rest_density: 1000.0,
            viscosity: 0.01,
            solver_type: SolverType::Wcsph,
            min_dt: 1e-5,
            max_dt: 0.001,
            ..Default::default()
        };
        let mut sim = SphSimulation::new(params);
        sim.boundaries
            .add_plane(BoundaryPlane::new(Vec3::zeros(), Vec3::new(0.0, 1.0, 0.0)));
        let spacing: f64 = 0.04;
        let mass = 1000.0 * spacing.powi(3);
        for i in 0..4 {
            for j in 0..4 {
                for k in 0..4 {
                    let pos = Vec3::new(
                        0.1 + i as f64 * spacing,
                        0.2 + j as f64 * spacing,
                        0.1 + k as f64 * spacing,
                    );
                    sim.particles
                        .add_particle(&SphParticle::new(pos, Vec3::zeros(), mass));
                }
            }
        }
        let initial_y: Vec<f64> = sim.particles.positions.iter().map(|p| p.y).collect();
        let kernel = CubicSplineKernel;
        for _ in 0..20 {
            sim.step(0.0005, &kernel);
        }
        let final_y: Vec<f64> = sim.particles.positions.iter().map(|p| p.y).collect();
        let avg_initial: f64 = initial_y.iter().sum::<f64>() / initial_y.len() as f64;
        let avg_final: f64 = final_y.iter().sum::<f64>() / final_y.len() as f64;
        assert!(
            avg_final < avg_initial,
            "Particles should fall: avg_y went from {avg_initial} to {avg_final}"
        );
    }
    #[test]
    fn sim_config_default_water() {
        let cfg = SphSimConfig::default_water();
        assert!((cfg.rest_density - 1000.0).abs() < 1e-10);
        assert!(cfg.dt > 0.0);
        assert!(cfg.kernel_radius > 0.0);
    }
    #[test]
    fn sim_config_cfl_dt_positive() {
        let cfg = SphSimConfig::default_water();
        let dt = cfg.cfl_dt();
        assert!(dt > 0.0, "CFL dt must be positive");
        assert!(dt < 1.0, "CFL dt should be small");
    }
    #[test]
    fn sim_config_sound_speed_positive() {
        let cfg = SphSimConfig::default_water();
        let c0 = cfg.estimated_sound_speed();
        assert!(c0 > 0.0);
    }
    #[test]
    fn sim_state_advance() {
        let mut s = SphSimState::new();
        assert_eq!(s.step, 0);
        assert!((s.time).abs() < 1e-14);
        s.advance(0.001);
        assert_eq!(s.step, 1);
        assert!((s.time - 0.001).abs() < 1e-14);
    }
    #[test]
    fn sph_sim_new_creates_particles() {
        let cfg = SphSimConfig::default_water();
        let sim = SphSim::new(cfg, 8);
        assert_eq!(sim.len(), 8);
        assert!(!sim.is_empty());
    }
    #[test]
    fn sph_sim_step_advances_time() {
        let cfg = SphSimConfig::default_water();
        let mut sim = SphSim::new(cfg, 8);
        sim.step();
        assert_eq!(sim.state.step, 1);
        assert!(sim.state.time > 0.0);
    }
    #[test]
    fn sph_sim_density_sum_positive() {
        let cfg = SphSimConfig::default_water();
        let mut sim = SphSim::new(cfg, 27);
        sim.compute_density_sum();
        for &d in &sim.densities {
            assert!(d > 0.0, "density must be positive after summation");
        }
    }
    #[test]
    fn sph_sim_kinetic_energy_zero_initially() {
        let cfg = SphSimConfig::default_water();
        let sim = SphSim::new(cfg, 8);
        assert!((sim.kinetic_energy()).abs() < 1e-14);
    }
    #[test]
    fn sph_sim_gravity_increases_ke() {
        let cfg = SphSimConfig::default_water();
        let mut sim = SphSim::new(cfg, 8);
        let ke_before = sim.kinetic_energy();
        for _ in 0..5 {
            sim.step();
        }
        let ke_after = sim.kinetic_energy();
        assert!(
            ke_after > ke_before,
            "Gravity should increase KE: before={ke_before}, after={ke_after}"
        );
    }
    #[test]
    fn sph_sim_cubic_kernel_positive_at_zero() {
        let w = SphSim::cubic_w(0.0, 0.1);
        assert!(w > 0.0, "kernel at r=0 must be positive");
    }
    #[test]
    fn sph_sim_cubic_kernel_zero_outside_support() {
        let w = SphSim::cubic_w(0.25, 0.1);
        assert!((w).abs() < 1e-14, "kernel outside support must be zero");
    }
}
#[cfg(test)]
mod tests_extended {
    use super::*;

    #[test]
    fn sph_phase_water_density() {
        let p = SphPhase::water();
        assert!((p.rest_density - 1000.0).abs() < 1e-10);
    }
    #[test]
    fn sph_phase_tait_pressure_at_rest() {
        let p = SphPhase::water();
        let pressure = p.tait_pressure(1000.0);
        assert!(
            pressure.abs() < 1e-6,
            "P at rest density should be 0, got {pressure}"
        );
    }
    #[test]
    fn sph_phase_tait_pressure_compressed() {
        let p = SphPhase::water();
        let pressure = p.tait_pressure(1100.0);
        assert!(
            pressure > 0.0,
            "compressed water should have positive pressure"
        );
    }
    #[test]
    fn sph_phase_oil_has_higher_viscosity_than_air() {
        let oil = SphPhase::oil();
        let air = SphPhase::air();
        assert!(oil.viscosity > air.viscosity);
    }
    #[test]
    fn adaptive_timestep_select_clamps_to_range() {
        let mut ats = AdaptiveTimestep::new(0.001, 1e-6, 0.01);
        let dt = ats.select(1000.0, 0.1, 100.0, 1e-6);
        assert!(dt >= 1e-6);
        assert!(dt <= 0.01);
    }
    #[test]
    fn adaptive_timestep_mean_dt_monotone() {
        let mut ats = AdaptiveTimestep::new(0.001, 1e-6, 0.01);
        for _ in 0..5 {
            ats.select(10.0, 0.1, 100.0, 1e-4);
        }
        assert_eq!(ats.steps_taken(), 5);
        let mean = ats.mean_dt();
        assert!(mean > 0.0 && mean.is_finite());
    }
    #[test]
    fn adaptive_timestep_cfl_limit_positive() {
        let ats = AdaptiveTimestep::new(0.001, 1e-6, 0.01);
        let dt_cfl = ats.cfl_limit(0.0, 0.1, 100.0);
        assert!(dt_cfl > 0.0);
    }
    #[test]
    fn adaptive_timestep_viscous_limit_large_for_small_nu() {
        let ats = AdaptiveTimestep::new(0.001, 1e-8, 1.0);
        let dt_v = ats.viscous_limit(0.1, 1e-6);
        assert!(dt_v > 0.0 && dt_v.is_finite());
    }
    #[test]
    fn simulation_stats_initial_state() {
        let s = SimulationStats::new();
        assert_eq!(s.num_steps, 0);
        assert_eq!(s.total_time, 0.0);
        assert!(s.min_dt.is_infinite());
    }
    #[test]
    fn simulation_stats_record_step() {
        let mut s = SimulationStats::new();
        s.record_step(0.001, 1.5, 0.3);
        assert_eq!(s.num_steps, 1);
        assert!((s.total_time - 0.001).abs() < 1e-15);
        assert!((s.mean_dt() - 0.001).abs() < 1e-15);
        assert!((s.last_kinetic_energy - 1.5).abs() < 1e-14);
        assert!((s.last_max_speed - 0.3).abs() < 1e-14);
    }
    #[test]
    fn simulation_stats_min_max_dt() {
        let mut s = SimulationStats::new();
        s.record_step(0.01, 0.0, 0.0);
        s.record_step(0.001, 0.0, 0.0);
        s.record_step(0.005, 0.0, 0.0);
        assert!((s.min_dt - 0.001).abs() < 1e-15);
        assert!((s.max_dt - 0.01).abs() < 1e-15);
    }
    #[test]
    fn sim_snapshot_from_sim() {
        let cfg = SphSimConfig::default_water();
        let sim = SphSim::new(cfg, 8);
        let snap = SimSnapshot::from_sim(&sim);
        assert_eq!(snap.len(), 8);
        assert!(!snap.is_empty());
        assert_eq!(snap.step, 0);
    }
    #[test]
    fn sim_snapshot_centre_of_mass_empty() {
        let snap = SimSnapshot {
            time: 0.0,
            step: 0,
            positions: vec![],
            velocities: vec![],
            densities: vec![],
        };
        let com = snap.centre_of_mass();
        assert_eq!(com, [0.0; 3]);
    }
    #[test]
    fn sim_snapshot_centre_of_mass_single() {
        let snap = SimSnapshot {
            time: 0.0,
            step: 0,
            positions: vec![[1.0, 2.0, 3.0]],
            velocities: vec![[0.0; 3]],
            densities: vec![1000.0],
        };
        let com = snap.centre_of_mass();
        assert!((com[0] - 1.0).abs() < 1e-14);
        assert!((com[1] - 2.0).abs() < 1e-14);
        assert!((com[2] - 3.0).abs() < 1e-14);
    }
    #[test]
    fn sph_runner_runs_n_steps() {
        let cfg = SphSimConfig::default_water();
        let sim = SphSim::new(cfg, 8);
        let mut runner = SphRunner::new(sim, 0);
        runner.run_steps(5);
        assert_eq!(runner.stats.num_steps, 5);
    }
    #[test]
    fn sph_runner_snapshots_collected() {
        let cfg = SphSimConfig::default_water();
        let sim = SphSim::new(cfg, 8);
        let mut runner = SphRunner::new(sim, 2);
        runner.run_steps(6);
        assert_eq!(runner.num_snapshots(), 3);
    }
    #[test]
    fn sph_runner_stats_min_max_dt() {
        let cfg = SphSimConfig::default_water();
        let sim = SphSim::new(cfg, 4);
        let mut runner = SphRunner::new(sim, 0);
        runner.run_steps(3);
        assert!(runner.stats.min_dt.is_finite());
        assert!(runner.stats.max_dt.is_finite());
        assert!(runner.stats.max_dt >= runner.stats.min_dt);
    }
    #[test]
    fn multi_phase_add_phase_and_particle() {
        let mut mp = MultiPhaseParticleSet::new();
        let w_id = mp.add_phase(SphPhase::water());
        let o_id = mp.add_phase(SphPhase::oil());
        mp.add_particle([0.0; 3], [0.0; 3], 0.001, w_id);
        mp.add_particle([1.0, 0.0, 0.0], [0.0; 3], 0.001, o_id);
        assert_eq!(mp.len(), 2);
        assert!(!mp.is_empty());
        let counts = mp.count_by_phase();
        assert_eq!(counts[w_id], 1);
        assert_eq!(counts[o_id], 1);
    }
    #[test]
    fn multi_phase_update_pressures_at_rest() {
        let mut mp = MultiPhaseParticleSet::new();
        let w_id = mp.add_phase(SphPhase::water());
        mp.add_particle([0.0; 3], [0.0; 3], 1.0, w_id);
        mp.densities[0] = 1000.0;
        mp.update_pressures();
        assert!(
            mp.pressures[0].abs() < 1e-3,
            "P at rest density ~0, got {}",
            mp.pressures[0]
        );
    }
    #[test]
    fn multi_phase_kinetic_energy_stationary() {
        let mut mp = MultiPhaseParticleSet::new();
        let w_id = mp.add_phase(SphPhase::water());
        mp.add_particle([0.0; 3], [0.0; 3], 1.0, w_id);
        assert!((mp.kinetic_energy()).abs() < 1e-14);
    }
    #[test]
    fn multi_phase_kinetic_energy_moving() {
        let mut mp = MultiPhaseParticleSet::new();
        let w_id = mp.add_phase(SphPhase::water());
        mp.add_particle([0.0; 3], [2.0, 0.0, 0.0], 1.0, w_id);
        assert!((mp.kinetic_energy() - 2.0).abs() < 1e-14);
    }
    #[test]
    fn multi_phase_centre_of_mass_symmetric() {
        let mut mp = MultiPhaseParticleSet::new();
        let w_id = mp.add_phase(SphPhase::water());
        mp.add_particle([-1.0, 0.0, 0.0], [0.0; 3], 1.0, w_id);
        mp.add_particle([1.0, 0.0, 0.0], [0.0; 3], 1.0, w_id);
        let com = mp.centre_of_mass();
        assert!(com[0].abs() < 1e-14, "com.x should be 0, got {}", com[0]);
    }
    #[test]
    fn multi_phase_default_is_empty() {
        let mp = MultiPhaseParticleSet::default();
        assert!(mp.is_empty());
    }
}
#[cfg(test)]
mod tests_periodic {
    use super::*;
    use crate::particle::SphParticle;

    /// PBC wrap: positions outside \[0, L) must be folded back.
    #[test]
    fn pbc_wraps_positions_into_box() {
        let cfg = SphSimConfig::default_water();
        let mut inner = SphSim::new(cfg, 1);
        inner.positions[0] = [1.5, -0.1, 2.3];
        let box_size = [1.0, 1.0, 1.0];
        let mut sim = PeriodicSphSim::new(inner, box_size);
        sim.apply_pbc();
        for k in 0..3 {
            assert!(
                sim.inner.positions[0][k] >= 0.0,
                "pos[{k}] = {} must be >= 0",
                sim.inner.positions[0][k]
            );
            assert!(
                sim.inner.positions[0][k] < 1.0,
                "pos[{k}] = {} must be < 1",
                sim.inner.positions[0][k]
            );
        }
    }
    /// PBC wrap: positions already inside box must not change.
    #[test]
    fn pbc_does_not_move_particles_already_in_box() {
        let cfg = SphSimConfig::default_water();
        let mut inner = SphSim::new(cfg, 1);
        inner.positions[0] = [0.3, 0.7, 0.5];
        let box_size = [1.0, 1.0, 1.0];
        let mut sim = PeriodicSphSim::new(inner, box_size);
        sim.apply_pbc();
        assert!((sim.inner.positions[0][0] - 0.3).abs() < 1e-14);
        assert!((sim.inner.positions[0][1] - 0.7).abs() < 1e-14);
        assert!((sim.inner.positions[0][2] - 0.5).abs() < 1e-14);
    }
    /// PBC wrap: particle exactly at L wraps to 0.
    #[test]
    fn pbc_wraps_exactly_at_boundary() {
        let cfg = SphSimConfig::default_water();
        let mut inner = SphSim::new(cfg, 1);
        inner.positions[0] = [1.0, 1.0, 1.0];
        let mut sim = PeriodicSphSim::new(inner, [1.0; 3]);
        sim.apply_pbc();
        for k in 0..3 {
            assert!(
                (sim.inner.positions[0][k]).abs() < 1e-14,
                "pos[{k}]={} expected 0",
                sim.inner.positions[0][k]
            );
        }
    }
    /// Momentum conservation: no external forces → total momentum constant.
    ///
    /// Uses a gravity-free configuration so pressure + viscosity interactions
    /// conserve momentum (symmetric kernel summation).
    #[test]
    fn momentum_conserved_without_gravity() {
        let mut cfg = SphSimConfig::default_water();
        cfg.gravity = [0.0; 3];
        let inner = SphSim::new(cfg, 8);
        let mut sim = PeriodicSphSim::new(inner, [1.0; 3]);
        let mom_before = sim.total_momentum();
        for _ in 0..5 {
            sim.step();
        }
        let mom_after = sim.total_momentum();
        for k in 0..3 {
            let delta = (mom_after[k] - mom_before[k]).abs();
            assert!(
                delta < 1e-6,
                "momentum[{k}] changed by {delta}: {mom_before:?} → {mom_after:?}"
            );
        }
    }
    /// Energy with gravity: total mechanical energy should change (gravity does work).
    #[test]
    fn total_energy_changes_under_gravity() {
        let cfg = SphSimConfig::default_water();
        let inner = SphSim::new(cfg, 8);
        let mut sim = PeriodicSphSim::new(inner, [2.0; 3]);
        let e_before = sim.total_energy();
        for _ in 0..10 {
            sim.step();
        }
        let e_after = sim.total_energy();
        assert!(
            (e_after - e_before).abs() > 0.0,
            "energy must change: {e_before} → {e_after}"
        );
    }
    /// Total momentum is finite after steps with gravity.
    #[test]
    fn total_momentum_finite_after_steps() {
        let cfg = SphSimConfig::default_water();
        let inner = SphSim::new(cfg, 8);
        let mut sim = PeriodicSphSim::new(inner, [2.0; 3]);
        for _ in 0..5 {
            sim.step();
        }
        let mom = sim.total_momentum();
        for (k, &m) in mom.iter().enumerate() {
            assert!(m.is_finite(), "momentum[{k}] not finite: {}", m);
        }
    }
    /// Surface tension step: total energy must remain finite.
    #[test]
    fn step_with_surface_tension_energy_finite() {
        let cfg = SphSimConfig::default_water();
        let inner = SphSim::new(cfg, 8);
        let mut sim = PeriodicSphSim::new(inner, [2.0; 3]);
        for _ in 0..3 {
            sim.step_with_surface_tension();
        }
        let e = sim.total_energy();
        assert!(
            e.is_finite(),
            "energy must be finite after surface tension steps: {e}"
        );
    }
    /// Surface tension forces are zero when sigma = 0.
    #[test]
    fn surface_tension_zero_when_sigma_is_zero() {
        let mut cfg = SphSimConfig::default_water();
        cfg.surface_tension_coeff = 0.0;
        let mut sim = SphSim::new(cfg, 8);
        sim.compute_density_sum();
        let forces_before: Vec<[f64; 3]> = sim.forces.clone();
        sim.compute_surface_tension_forces();
        for (fb, fa) in forces_before.iter().zip(sim.forces.iter()) {
            for k in 0..3 {
                assert!((fb[k] - fa[k]).abs() < 1e-30, "force changed with sigma=0");
            }
        }
    }
    /// `SphSimulation::apply_pbc` wraps positions on its own particle set.
    #[test]
    fn sph_simulation_apply_pbc_wraps() {
        let params = SphSimulationParams::default();
        let mut sim = SphSimulation::new(params);
        sim.particles.add_particle(&SphParticle::new(
            Vec3::new(1.5, -0.2, 3.1),
            Vec3::zeros(),
            1.0,
        ));
        sim.apply_pbc([1.0, 1.0, 1.0]);
        let p = &sim.particles.positions[0];
        assert!(p.x >= 0.0 && p.x < 1.0, "x={}", p.x);
        assert!(p.y >= 0.0 && p.y < 1.0, "y={}", p.y);
        assert!(p.z >= 0.0 && p.z < 1.0, "z={}", p.z);
    }
    /// `SphSimulation::total_momentum` returns finite values.
    #[test]
    fn sph_simulation_total_momentum_finite() {
        let params = SphSimulationParams::default();
        let mut sim = SphSimulation::new(params);
        sim.particles.add_particle(&SphParticle::new(
            Vec3::new(0.1, 0.2, 0.3),
            Vec3::new(1.0, -0.5, 0.3),
            2.0,
        ));
        let mom = sim.total_momentum();
        for (k, &m) in mom.iter().enumerate() {
            assert!(m.is_finite(), "momentum[{k}] not finite");
        }
    }
    /// `SphSimulation::total_energy` returns finite value.
    #[test]
    fn sph_simulation_total_energy_finite() {
        let params = SphSimulationParams::default();
        let mut sim = SphSimulation::new(params);
        sim.particles.add_particle(&SphParticle::new(
            Vec3::new(0.1, 1.0, 0.3),
            Vec3::new(0.5, 0.2, 0.0),
            1.0,
        ));
        let e = sim.total_energy();
        assert!(e.is_finite(), "energy must be finite: {e}");
    }
    /// `PeriodicSphSim::len` and `is_empty` work correctly.
    #[test]
    fn periodic_sim_len_and_is_empty() {
        let cfg = SphSimConfig::default_water();
        let inner = SphSim::new(cfg, 4);
        let sim = PeriodicSphSim::new(inner, [1.0; 3]);
        assert_eq!(sim.len(), 4);
        assert!(!sim.is_empty());
    }
}
#[cfg(test)]
mod tests_wcsph_sim {
    use super::*;
    use crate::kernel::CubicSplineKernel;
    use crate::particle::SphParticle;

    fn make_sim_with_particles(n: usize) -> WcSphSim {
        let mut sim = WcSphSimBuilder::water().with_h(0.15).build();
        let spacing = 0.05_f64;
        let side = (n as f64).cbrt().ceil() as usize;
        let mass = sim.rho0 * spacing.powi(3);
        let mut count = 0;
        'outer: for i in 0..side {
            for j in 0..side {
                for k in 0..side {
                    if count >= n {
                        break 'outer;
                    }
                    sim.particles.add_raw(
                        [i as f64 * spacing, j as f64 * spacing, k as f64 * spacing],
                        [0.0; 3],
                        mass,
                        0.0,
                        false,
                    );
                    count += 1;
                }
            }
        }
        sim
    }
    #[test]
    fn wcsph_sim_new_empty() {
        let sim = WcSphSimBuilder::water().build();
        assert!(sim.is_empty());
        assert_eq!(sim.len(), 0);
    }
    #[test]
    fn wcsph_sim_builder_sets_h() {
        let sim = WcSphSimBuilder::water().with_h(0.2).build();
        assert!((sim.h - 0.2).abs() < 1e-14);
    }
    #[test]
    fn wcsph_sim_builder_sets_rho0() {
        let sim = WcSphSimBuilder::water().with_rho0(800.0).build();
        assert!((sim.rho0 - 800.0).abs() < 1e-14);
    }
    #[test]
    fn wcsph_sim_builder_sets_mu() {
        let sim = WcSphSimBuilder::water().with_mu(0.1).build();
        assert!((sim.mu - 0.1).abs() < 1e-14);
    }
    #[test]
    fn wcsph_sim_tait_b_positive() {
        let sim = WcSphSimBuilder::water().build();
        assert!(sim.tait_b() > 0.0);
    }
    #[test]
    fn wcsph_sim_tait_b_formula() {
        let sim = WcSphSim::new(0.1, 1000.0, 100.0, 7.0, 1e-3, 0.0, [0.0; 3]);
        let expected = 1000.0 * 100.0 * 100.0 / 7.0;
        assert!((sim.tait_b() - expected).abs() < 1e-6);
    }
    #[test]
    fn wcsph_sim_density_step_positive() {
        let mut sim = make_sim_with_particles(8);
        sim.density_step();
        for &d in &sim.particles.densities {
            assert!(d > 0.0, "density must be positive after density step");
        }
    }
    #[test]
    fn wcsph_sim_pressure_step_at_rest_near_zero() {
        let mut sim = make_sim_with_particles(8);
        sim.density_step();
        for d in &mut sim.particles.densities {
            *d = sim.rho0;
        }
        sim.pressure_step();
        for &p in &sim.particles.pressures {
            assert!(p.abs() < 1.0, "P at rest density should be ~0, got {p}");
        }
    }
    #[test]
    fn wcsph_sim_force_step_no_nan() {
        let mut sim = make_sim_with_particles(8);
        sim.density_step();
        sim.pressure_step();
        sim.force_step();
        for a in &sim.particles.accelerations {
            assert!(
                a[0].is_finite() && a[1].is_finite() && a[2].is_finite(),
                "acceleration must be finite"
            );
        }
    }
    #[test]
    fn wcsph_sim_gravity_in_acceleration() {
        let mut sim = make_sim_with_particles(4);
        sim.density_step();
        sim.pressure_step();
        sim.force_step();
        let ay_sum: f64 = sim.particles.accelerations.iter().map(|a| a[1]).sum();
        assert!(
            ay_sum < 0.0,
            "net y-acceleration should be negative under gravity, got {ay_sum}"
        );
    }
    #[test]
    fn wcsph_sim_step_advances_time() {
        let mut sim = make_sim_with_particles(8);
        let dt = 1e-4;
        sim.step(dt);
        assert!((sim.time - dt).abs() < 1e-15, "time should advance by dt");
    }
    #[test]
    fn wcsph_sim_step_records_energy() {
        let mut sim = make_sim_with_particles(8);
        sim.step(1e-4);
        assert_eq!(sim.kinetic_energy_history.len(), 1);
        assert_eq!(sim.potential_energy_history.len(), 1);
    }
    #[test]
    fn wcsph_sim_gravity_increases_kinetic_energy() {
        let mut sim = make_sim_with_particles(8);
        let ke0 = sim.kinetic_energy();
        for _ in 0..10 {
            sim.step(1e-4);
        }
        let ke1 = sim.kinetic_energy();
        assert!(ke1 > ke0, "gravity must increase KE: {ke0} → {ke1}");
    }
    #[test]
    fn wcsph_sim_total_energy_finite() {
        let mut sim = make_sim_with_particles(8);
        sim.step(1e-4);
        assert!(sim.total_energy().is_finite());
    }
    #[test]
    fn wcsph_sim_total_momentum_finite() {
        let mut sim = make_sim_with_particles(8);
        sim.step(1e-4);
        let mom = sim.total_momentum();
        for (k, &m) in mom.iter().enumerate() {
            assert!(m.is_finite(), "momentum[{k}] must be finite");
        }
    }
    #[test]
    fn wcsph_sim_momentum_zero_gravity_conserved() {
        let mut sim = WcSphSimBuilder::water()
            .with_gravity([0.0; 3])
            .with_h(0.15)
            .build();
        let spacing = 0.05_f64;
        let mass = 1000.0 * spacing.powi(3);
        for i in 0..4 {
            for j in 0..4 {
                sim.particles.add_raw(
                    [i as f64 * spacing, j as f64 * spacing, 0.0],
                    [0.0; 3],
                    mass,
                    0.0,
                    false,
                );
            }
        }
        let mom_before = sim.total_momentum();
        for _ in 0..5 {
            sim.step(1e-5);
        }
        let mom_after = sim.total_momentum();
        for k in 0..3 {
            let delta = (mom_after[k] - mom_before[k]).abs();
            assert!(
                delta < 1e-6,
                "momentum[{k}] should be conserved (no gravity), delta={delta}"
            );
        }
    }
    #[test]
    fn wcsph_sim_leapfrog_step_finite_energy() {
        let mut sim = make_sim_with_particles(8);
        sim.density_step();
        sim.pressure_step();
        sim.force_step();
        sim.step_leapfrog(1e-4);
        assert!(sim.total_energy().is_finite());
    }
    #[test]
    fn wcsph_sim_verlet_step_advances_time() {
        let mut sim = make_sim_with_particles(8);
        sim.density_step();
        sim.pressure_step();
        sim.force_step();
        let dt = 1e-4;
        sim.step_verlet(dt);
        assert!((sim.time - dt).abs() < 1e-15);
    }
    #[test]
    fn wcsph_sim_cfl_dt_within_bounds() {
        let sim = make_sim_with_particles(4);
        let dt = sim.cfl_dt(0.4, 1e-6, 0.01);
        assert!((1e-6..=0.01).contains(&dt), "dt={dt} out of [1e-6, 0.01]");
    }
    #[test]
    fn wcsph_sim_run_adaptive_reaches_target_time() {
        let mut sim = make_sim_with_particles(4);
        let target = 1e-3;
        sim.run_adaptive(target, 0.4, 1e-7, 1e-3);
        assert!(
            sim.time >= target * 0.999,
            "sim.time = {} should be ≥ {target}",
            sim.time
        );
    }
    #[test]
    fn wcsph_sim_surface_tension_with_sigma_finite() {
        let mut sim = WcSphSimBuilder::water()
            .with_sigma(0.0728)
            .with_h(0.15)
            .build();
        let spacing = 0.05_f64;
        let mass = 1000.0 * spacing.powi(3);
        for i in 0..4 {
            sim.particles
                .add_raw([i as f64 * spacing, 0.0, 0.0], [0.0; 3], mass, 0.0, false);
        }
        sim.step(1e-4);
        assert!(
            sim.total_energy().is_finite(),
            "energy must remain finite with surface tension"
        );
    }
    #[test]
    fn energy_tracker_new_is_empty() {
        let et = EnergyTracker::new();
        assert!(et.is_empty());
        assert_eq!(et.len(), 0);
    }
    #[test]
    fn energy_tracker_record_and_len() {
        let mut et = EnergyTracker::new();
        et.record(0.0, 1.0, 2.0);
        et.record(0.1, 2.0, 1.5);
        assert_eq!(et.len(), 2);
        assert!(!et.is_empty());
    }
    #[test]
    fn energy_tracker_total_energies_correct() {
        let mut et = EnergyTracker::new();
        et.record(0.0, 3.0, 4.0);
        et.record(0.1, 1.0, 2.0);
        let totals = et.total_energies();
        assert!((totals[0] - 7.0).abs() < 1e-14);
        assert!((totals[1] - 3.0).abs() < 1e-14);
    }
    #[test]
    fn energy_tracker_mean_kinetic() {
        let mut et = EnergyTracker::new();
        et.record(0.0, 2.0, 0.0);
        et.record(0.1, 4.0, 0.0);
        assert!((et.mean_kinetic() - 3.0).abs() < 1e-14);
    }
    #[test]
    fn energy_tracker_max_kinetic() {
        let mut et = EnergyTracker::new();
        et.record(0.0, 2.0, 0.0);
        et.record(0.1, 7.0, 0.0);
        et.record(0.2, 3.0, 0.0);
        assert!((et.max_kinetic() - 7.0).abs() < 1e-14);
    }
    #[test]
    fn cfl_controller_standard_dt_within_bounds() {
        let ctrl = CflController::standard(1e-6, 0.01);
        let dt = ctrl.combined_dt(0.1, 100.0, 0.0, 1e-6);
        assert!((1e-6..=0.01).contains(&dt));
    }
    #[test]
    fn cfl_controller_cfl_dt_formula() {
        let ctrl = CflController::standard(1e-8, 1.0);
        let dt = ctrl.cfl_dt(0.1, 100.0, 0.0);
        assert!((dt - 4e-4).abs() < 1e-14, "dt = {dt}");
    }
    #[test]
    fn cfl_controller_viscous_dt_formula() {
        let ctrl = CflController::standard(1e-8, 1.0);
        let dt = ctrl.viscous_dt(0.1, 1e-4);
        assert!((dt - 12.5).abs() < 1e-10, "dt = {dt}");
    }
    #[test]
    fn cfl_controller_combined_min_of_two() {
        let ctrl = CflController {
            cfl: 0.4,
            viscous_factor: 0.125,
            dt_min: 1e-10,
            dt_max: 1.0,
        };
        let dt_cfl = ctrl.cfl_dt(0.1, 100.0, 0.0);
        let dt_visc = ctrl.viscous_dt(0.1, 1e-4);
        let dt_combined = ctrl.combined_dt(0.1, 100.0, 0.0, 1e-4);
        assert!((dt_combined - dt_cfl.min(dt_visc)).abs() < 1e-14);
    }
    #[test]
    fn cfl_controller_zero_nu_returns_cfl_limited() {
        let ctrl = CflController::standard(1e-8, 1.0);
        let dt = ctrl.combined_dt(0.1, 100.0, 0.0, 0.0);
        let expected_cfl = ctrl.cfl_dt(0.1, 100.0, 0.0);
        let expected = expected_cfl.clamp(1e-8, 1.0);
        assert!((dt - expected).abs() < 1e-14, "dt={dt} expected={expected}");
    }
    #[test]
    fn sph_sim_config_default_water_gravity_negative_y() {
        let cfg = SphSimConfig::default_water();
        assert!(cfg.gravity[1] < 0.0, "gravity must point down");
    }
    #[test]
    fn sph_sim_config_cfl_dt_inversely_proportional_to_sound_speed() {
        let mut cfg = SphSimConfig::default_water();
        let dt1 = cfg.cfl_dt();
        cfg.gravity = [0.0, -9.81 * 4.0, 0.0];
        let dt2 = cfg.cfl_dt();
        assert!(dt2 < dt1, "larger gravity → larger c0 → smaller CFL dt");
    }
    #[test]
    fn sph_sim_step_multiple_times_increases_time() {
        let cfg = SphSimConfig::default_water();
        let mut sim = SphSim::new(cfg, 8);
        let t0 = sim.state.time;
        for _ in 0..10 {
            sim.step();
        }
        assert!(sim.state.time > t0);
        assert_eq!(sim.state.step, 10);
    }
    #[test]
    fn sph_sim_pressure_forces_symmetric() {
        let mut cfg = SphSimConfig::default_water();
        cfg.gravity = [0.0; 3];
        let mut sim = SphSim::new(cfg, 0);
        let m = 1.0;
        sim.positions.push([0.0, 0.0, 0.0]);
        sim.positions.push([0.08, 0.0, 0.0]);
        sim.velocities.push([0.0; 3]);
        sim.velocities.push([0.0; 3]);
        sim.densities.push(1000.0);
        sim.densities.push(1000.0);
        sim.pressures.push(0.0);
        sim.pressures.push(0.0);
        sim.forces.push([0.0; 3]);
        sim.forces.push([0.0; 3]);
        sim.masses.push(m);
        sim.masses.push(m);
        sim.compute_density_sum();
        sim.compute_pressure_forces();
        let fx_total = sim.forces[0][0] + sim.forces[1][0];
        assert!(fx_total.abs() < 1e-6, "Total force not zero: {fx_total}");
    }
    #[test]
    fn sph_sim_zero_particles_step_ok() {
        let cfg = SphSimConfig::default_water();
        let mut sim = SphSim::new(cfg, 0);
        sim.step();
        assert_eq!(sim.state.step, 1);
    }
    #[test]
    fn sph_sim_positions_in_unit_cube() {
        let cfg = SphSimConfig::default_water();
        let sim = SphSim::new(cfg, 27);
        for pos in &sim.positions {
            for &c in pos {
                assert!(
                    (0.0..=1.0).contains(&c),
                    "position component out of [0,1]: {c}"
                );
            }
        }
    }
    #[test]
    fn adaptive_timestep_cfl_limit_decreases_with_high_velocity() {
        let mut ctrl = AdaptiveTimestep::new(0.001, 1e-6, 0.01);
        let dt_slow = ctrl.select(0.1, 0.1, 100.0, 1e-4);
        let mut ctrl2 = AdaptiveTimestep::new(0.001, 1e-6, 0.01);
        let dt_fast = ctrl2.select(1000.0, 0.1, 100.0, 1e-4);
        assert!(dt_fast < dt_slow, "Higher velocity should reduce dt");
    }
    #[test]
    fn adaptive_timestep_records_history() {
        let mut ctrl = AdaptiveTimestep::new(0.001, 1e-8, 0.01);
        ctrl.select(0.1, 0.1, 100.0, 1e-4);
        ctrl.select(0.2, 0.1, 100.0, 1e-4);
        ctrl.select(0.5, 0.1, 100.0, 1e-4);
        assert_eq!(ctrl.steps_taken(), 3);
    }
    #[test]
    fn adaptive_timestep_mean_dt_positive() {
        let mut ctrl = AdaptiveTimestep::new(0.001, 1e-8, 0.01);
        for i in 1..=5 {
            ctrl.select(i as f64 * 0.1, 0.1, 100.0, 1e-4);
        }
        assert!(ctrl.mean_dt() > 0.0);
    }
    #[test]
    fn simulation_stats_record_and_mean_dt() {
        let mut stats = SimulationStats::new();
        stats.record_step(0.001, 10.0, 1.5);
        stats.record_step(0.002, 12.0, 2.0);
        assert_eq!(stats.num_steps, 2);
        assert!((stats.total_time - 0.003).abs() < 1e-14);
        assert!((stats.mean_dt() - 0.0015).abs() < 1e-14);
    }
    #[test]
    fn simulation_stats_min_max_dt() {
        let mut stats = SimulationStats::new();
        stats.record_step(0.005, 0.0, 0.0);
        stats.record_step(0.001, 0.0, 0.0);
        stats.record_step(0.003, 0.0, 0.0);
        assert!((stats.min_dt - 0.001).abs() < 1e-14);
        assert!((stats.max_dt - 0.005).abs() < 1e-14);
    }
    #[test]
    fn simulation_stats_zero_steps_mean_dt_zero() {
        let stats = SimulationStats::new();
        assert_eq!(stats.mean_dt(), 0.0);
    }
    #[test]
    fn sim_snapshot_from_sim_captures_positions() {
        let cfg = SphSimConfig::default_water();
        let sim = SphSim::new(cfg, 4);
        let snap = SimSnapshot::from_sim(&sim);
        assert_eq!(snap.len(), 4);
        assert_eq!(snap.positions.len(), 4);
        assert_eq!(snap.densities.len(), 4);
    }
    #[test]
    fn sim_snapshot_centre_of_mass_equal_mass() {
        let cfg = SphSimConfig::default_water();
        let mut sim = SphSim::new(cfg, 0);
        sim.positions = vec![[0.0; 3], [2.0, 0.0, 0.0]];
        sim.velocities = vec![[0.0; 3]; 2];
        sim.densities = vec![0.0; 2];
        sim.pressures = vec![0.0; 2];
        sim.forces = vec![[0.0; 3]; 2];
        sim.masses = vec![1.0; 2];
        let snap = SimSnapshot::from_sim(&sim);
        let com = snap.centre_of_mass();
        assert!((com[0] - 1.0).abs() < 1e-14);
    }
    #[test]
    fn sph_runner_snapshots_count() {
        let cfg = SphSimConfig::default_water();
        let sim = SphSim::new(cfg, 4);
        let mut runner = SphRunner::new(sim, 2);
        runner.run_steps(6);
        assert_eq!(runner.num_snapshots(), 3);
    }
    #[test]
    fn sph_runner_stats_steps_match() {
        let cfg = SphSimConfig::default_water();
        let sim = SphSim::new(cfg, 4);
        let mut runner = SphRunner::new(sim, 0);
        runner.run_steps(5);
        assert_eq!(runner.stats.num_steps, 5);
    }
    #[test]
    fn multi_phase_particle_set_add_phases_and_particles() {
        let mut mps = MultiPhaseParticleSet::new();
        let water_id = mps.add_phase(SphPhase::water());
        let air_id = mps.add_phase(SphPhase::air());
        mps.add_particle([0.0; 3], [0.0; 3], 0.001, water_id);
        mps.add_particle([1.0, 0.0, 0.0], [0.0; 3], 0.001, air_id);
        assert_eq!(mps.len(), 2);
        assert_eq!(mps.phase_ids[0], water_id);
        assert_eq!(mps.phase_ids[1], air_id);
    }
    #[test]
    fn sph_phase_tait_pressure_at_rest_density_is_zero() {
        let phase = SphPhase::water();
        let p = phase.tait_pressure(phase.rest_density);
        assert!(
            p.abs() < 1e-6,
            "Tait pressure at rest density should be ~0, got {p}"
        );
    }
    #[test]
    fn sph_phase_tait_pressure_compressed_positive() {
        let phase = SphPhase::water();
        let p = phase.tait_pressure(phase.rest_density * 1.1);
        assert!(p > 0.0, "Compressed fluid should have positive pressure");
    }
    #[test]
    fn sph_simulation_params_default_gravity_downward() {
        let params = SphSimulationParams::default();
        assert!(params.gravity.y < 0.0, "Default gravity should be downward");
    }
    #[test]
    fn sph_simulation_wcsph_run_short() {
        let params = SphSimulationParams {
            smoothing_length: 0.1,
            min_dt: 1e-4,
            max_dt: 0.001,
            ..Default::default()
        };
        let mut sim = SphSimulation::new(params);
        let spacing = 0.05_f64;
        let mass = 1000.0 * spacing.powi(3);
        for i in 0..2 {
            for j in 0..2 {
                let pos = Vec3::new(0.1 + i as f64 * spacing, 0.5 + j as f64 * spacing, 0.1);
                sim.particles
                    .add_particle(&SphParticle::new(pos, Vec3::zeros(), mass));
            }
        }
        let kernel = CubicSplineKernel;
        sim.run(0.005, &kernel);
        assert!(sim.time > 0.0, "Simulation time should advance");
    }
    #[test]
    fn sph_sim_cubic_gradient_zero_at_origin() {
        let dw = SphSim::cubic_grad_w(0.0, 0.1);
        assert!((dw).abs() < 1e-14, "gradient at r=0 must be zero, got {dw}");
    }
    #[test]
    fn sph_sim_cubic_gradient_zero_outside_support() {
        let dw = SphSim::cubic_grad_w(0.25, 0.1);
        assert!((dw).abs() < 1e-14, "gradient outside support must be 0");
    }
    #[test]
    fn sph_sim_dist2_correct() {
        let a = [0.0_f64; 3];
        let b = [3.0, 4.0, 0.0];
        let d2 = SphSim::dist2(&a, &b);
        assert!((d2 - 25.0).abs() < 1e-14);
    }
    #[test]
    fn courant_number_stationary_fluid_equals_c0_dt_over_h() {
        let params = SphSimulationParams {
            smoothing_length: 0.1,
            sound_speed: 100.0,
            min_dt: 1e-6,
            max_dt: 0.01,
            ..SphSimulationParams::default()
        };
        let mut sim = SphSimulation::new(params);
        let mass = 1000.0 * 0.05_f64.powi(3);
        sim.particles.add_particle(&SphParticle::new(
            Vec3::new(0.5, 0.5, 0.5),
            Vec3::zeros(),
            mass,
        ));
        let kernel = CubicSplineKernel;
        let _ = kernel;
        let dt = 0.0004;
        let c = sim.compute_courant_number(dt);
        assert!(
            (c - 0.4).abs() < 1e-12,
            "Stationary Courant number should be c0*dt/h=0.4, got {c}"
        );
    }
    #[test]
    fn courant_number_moving_particle_exceeds_stationary() {
        let params = SphSimulationParams {
            smoothing_length: 0.1,
            sound_speed: 100.0,
            min_dt: 1e-6,
            max_dt: 0.01,
            ..SphSimulationParams::default()
        };
        let mut sim = SphSimulation::new(params);
        let mass = 0.001;
        let mut p = SphParticle::new(Vec3::new(0.5, 0.5, 0.5), Vec3::new(10.0, 0.0, 0.0), mass);
        p.velocity = Vec3::new(10.0, 0.0, 0.0);
        sim.particles.add_particle(&p);
        let dt = 0.0004;
        let c_moving = sim.compute_courant_number(dt);
        assert!(
            c_moving > 0.4,
            "Moving particle should give C > 0.4, got {c_moving}"
        );
    }
    #[test]
    fn adaptive_timestep_stationary_gives_cfl_dt() {
        let params = SphSimulationParams {
            smoothing_length: 0.1,
            sound_speed: 100.0,
            min_dt: 1e-6,
            max_dt: 1.0,
            ..SphSimulationParams::default()
        };
        let mut sim = SphSimulation::new(params);
        let mass = 0.001;
        sim.particles.add_particle(&SphParticle::new(
            Vec3::new(0.5, 0.5, 0.5),
            Vec3::zeros(),
            mass,
        ));
        let dt = sim.adaptive_timestep(0.4);
        assert!(
            (dt - 0.0004).abs() < 1e-12,
            "Stationary adaptive dt should be 0.0004, got {dt}"
        );
    }
    #[test]
    fn adaptive_timestep_clamped_to_max_dt() {
        let params = SphSimulationParams {
            smoothing_length: 0.1,
            sound_speed: 1.0,
            min_dt: 1e-6,
            max_dt: 0.005,
            ..SphSimulationParams::default()
        };
        let sim = SphSimulation::new(params);
        let dt = sim.adaptive_timestep(0.4);
        assert!(
            (dt - 0.005).abs() < 1e-14,
            "dt should be clamped to max_dt=0.005, got {dt}"
        );
    }
    #[test]
    fn adaptive_timestep_clamped_to_min_dt() {
        let params = SphSimulationParams {
            smoothing_length: 0.1,
            sound_speed: 1e8,
            min_dt: 1e-3,
            max_dt: 1.0,
            ..SphSimulationParams::default()
        };
        let mut sim = SphSimulation::new(params);
        sim.particles.add_particle(&SphParticle::new(
            Vec3::new(0.5, 0.5, 0.5),
            Vec3::zeros(),
            0.001,
        ));
        let dt = sim.adaptive_timestep(0.4);
        assert!(
            (dt - 1e-3).abs() < 1e-14,
            "dt should be clamped to min_dt=1e-3, got {dt}"
        );
    }
    #[test]
    fn mechanical_energy_empty_simulation_is_zero() {
        let params = SphSimulationParams::default();
        let sim = SphSimulation::new(params);
        assert_eq!(sim.compute_mechanical_energy(), 0.0);
    }
    #[test]
    fn mechanical_energy_positive_for_elevated_stationary_particle() {
        let params = SphSimulationParams::default();
        let mut sim = SphSimulation::new(params);
        let mass = 1.0;
        sim.particles.add_particle(&SphParticle::new(
            Vec3::new(0.0, 1.0, 0.0),
            Vec3::zeros(),
            mass,
        ));
        let e = sim.compute_mechanical_energy();
        assert!(
            e > 0.0,
            "Elevated stationary particle should have positive mechanical energy"
        );
    }
    #[test]
    fn mechanical_energy_equals_total_energy() {
        let params = SphSimulationParams::default();
        let mut sim = SphSimulation::new(params);
        let mass = 2.0;
        sim.particles.add_particle(&SphParticle::new(
            Vec3::new(0.0, 2.0, 0.0),
            Vec3::new(1.0, 0.0, 0.0),
            mass,
        ));
        assert_eq!(sim.compute_mechanical_energy(), sim.total_energy());
    }
}
#[cfg(test)]
mod tests_simulation_ext {
    use super::*;

    #[test]
    fn ensemble_run_produces_correct_replica_count() {
        let cfg = SphSimConfig::default_water();
        let mut ens = EnsembleSimulation::new(cfg, 4, 3, 2);
        ens.run(0.0);
        assert_eq!(ens.final_ke.len(), 3);
    }
    #[test]
    fn ensemble_mean_ke_is_finite_and_non_negative() {
        let cfg = SphSimConfig::default_water();
        let mut ens = EnsembleSimulation::new(cfg, 4, 5, 3);
        ens.run(0.001);
        assert!(ens.mean_ke() >= 0.0);
        assert!(ens.mean_ke().is_finite());
    }
    #[test]
    fn ensemble_all_finite_after_run() {
        let cfg = SphSimConfig::default_water();
        let mut ens = EnsembleSimulation::new(cfg, 4, 3, 2);
        ens.run(0.0);
        assert!(ens.all_finite());
    }
    #[test]
    fn ensemble_std_ke_non_negative() {
        let cfg = SphSimConfig::default_water();
        let mut ens = EnsembleSimulation::new(cfg, 4, 4, 2);
        ens.run(0.001);
        assert!(ens.std_ke() >= 0.0);
    }
    #[test]
    fn ensemble_no_runs_mean_ke_is_zero() {
        let cfg = SphSimConfig::default_water();
        let ens = EnsembleSimulation::new(cfg, 4, 3, 2);
        assert_eq!(ens.mean_ke(), 0.0);
    }
    #[test]
    fn two_phase_total_particles_sum() {
        let cfg = SphSimConfig::default_water();
        let a = SphSim::new(cfg.clone(), 4);
        let b = SphSim::new(cfg, 6);
        let sim = TwoPhaseSimState::new(a, b, 0.07);
        assert_eq!(sim.total_particles(), 10);
    }
    #[test]
    fn two_phase_step_advances_time() {
        let cfg = SphSimConfig::default_water();
        let dt = cfg.dt;
        let a = SphSim::new(cfg.clone(), 2);
        let b = SphSim::new(cfg, 2);
        let mut sim = TwoPhaseSimState::new(a, b, 0.07);
        sim.step();
        assert!((sim.time - dt).abs() < 1e-14);
    }
    #[test]
    fn two_phase_total_ke_finite() {
        let cfg = SphSimConfig::default_water();
        let a = SphSim::new(cfg.clone(), 4);
        let b = SphSim::new(cfg, 4);
        let mut sim = TwoPhaseSimState::new(a, b, 0.07);
        for _ in 0..3 {
            sim.step();
        }
        assert!(sim.total_ke().is_finite());
    }
    #[test]
    fn two_phase_total_momentum_finite() {
        let cfg = SphSimConfig::default_water();
        let a = SphSim::new(cfg.clone(), 4);
        let b = SphSim::new(cfg, 4);
        let mut sim = TwoPhaseSimState::new(a, b, 0.07);
        sim.step();
        let mom = sim.total_momentum();
        for &m in &mom {
            assert!(m.is_finite());
        }
    }
    #[test]
    fn harness_total_steps_correct() {
        let cfg = SphSimConfig::default_water();
        let mut sim = SphSim::new(cfg, 4);
        let mut harness = BenchmarkHarness::new(2, 5);
        harness.run(&mut sim);
        assert_eq!(harness.total_steps(), 7);
    }
    #[test]
    fn harness_multiple_runs_accumulate() {
        let cfg = SphSimConfig::default_water();
        let mut sim = SphSim::new(cfg, 2);
        let mut harness = BenchmarkHarness::new(1, 2);
        harness.run(&mut sim);
        harness.run(&mut sim);
        assert_eq!(harness.runs, 2);
        assert_eq!(harness.total_steps(), 6);
    }
    #[test]
    fn apply_periodic_bc_wraps_position() {
        let params = SphSimulationParams::default();
        let mut sim = SphSimulation::new(params);
        sim.add_particle_raw([1.5, 0.5, 0.5], [0.0; 3], 0.001);
        sim.apply_periodic_bc([1.0; 3]);
        assert!(sim.particles.positions[0].x < 1.0, "x should wrap");
        assert!(
            sim.particles.positions[0].x >= 0.0,
            "x must not be negative"
        );
    }
    #[test]
    fn total_mass_sums_correctly() {
        let params = SphSimulationParams::default();
        let mut sim = SphSimulation::new(params);
        sim.add_particle_raw([0.0; 3], [0.0; 3], 1.0);
        sim.add_particle_raw([0.0; 3], [0.0; 3], 2.0);
        assert!((sim.total_mass() - 3.0).abs() < 1e-14);
    }
    #[test]
    fn center_of_mass_single_particle_at_position() {
        let params = SphSimulationParams::default();
        let mut sim = SphSimulation::new(params);
        sim.add_particle_raw([3.0, 4.0, 5.0], [0.0; 3], 1.0);
        let com = sim.center_of_mass();
        assert!((com[0] - 3.0).abs() < 1e-14);
        assert!((com[1] - 4.0).abs() < 1e-14);
    }
    #[test]
    fn center_of_mass_empty_is_zero() {
        let params = SphSimulationParams::default();
        let sim = SphSimulation::new(params);
        let com = sim.center_of_mass();
        assert_eq!(com, [0.0; 3]);
    }
    #[test]
    fn kinetic_energy_zero_for_stationary_particles() {
        let params = SphSimulationParams::default();
        let mut sim = SphSimulation::new(params);
        sim.add_particle_raw([0.5, 0.5, 0.5], [0.0; 3], 1.0);
        assert_eq!(sim.kinetic_energy_sim(), 0.0);
    }
    #[test]
    fn kinetic_energy_nonzero_for_moving_particle() {
        let params = SphSimulationParams::default();
        let mut sim = SphSimulation::new(params);
        sim.add_particle_raw([0.5, 0.5, 0.5], [1.0, 0.0, 0.0], 2.0);
        assert!((sim.kinetic_energy_sim() - 1.0).abs() < 1e-14);
    }
    #[test]
    fn total_momentum_stationary_is_zero() {
        let params = SphSimulationParams::default();
        let mut sim = SphSimulation::new(params);
        sim.add_particle_raw([0.5, 0.5, 0.5], [0.0; 3], 1.0);
        let mom = sim.total_momentum();
        for &m in &mom {
            assert_eq!(m, 0.0);
        }
    }
    #[test]
    fn bounding_box_single_particle() {
        let params = SphSimulationParams::default();
        let mut sim = SphSimulation::new(params);
        sim.add_particle_raw([2.0, 3.0, 4.0], [0.0; 3], 1.0);
        let (lo, hi) = sim.bounding_box();
        assert!((lo[0] - 2.0).abs() < 1e-14);
        assert!((hi[1] - 3.0).abs() < 1e-14);
    }
    #[test]
    fn bounding_box_empty_is_zeros() {
        let params = SphSimulationParams::default();
        let sim = SphSimulation::new(params);
        let (lo, hi) = sim.bounding_box();
        assert_eq!(lo, [0.0; 3]);
        assert_eq!(hi, [0.0; 3]);
    }
    #[test]
    fn mean_density_empty_is_zero() {
        let params = SphSimulationParams::default();
        let sim = SphSimulation::new(params);
        assert_eq!(sim.mean_density(), 0.0);
    }
    #[test]
    fn potential_energy_at_rest_is_finite() {
        let cfg = SphSimConfig::default_water();
        let sim = SphSim::new(cfg, 4);
        assert!(sim.potential_energy().is_finite());
    }
    #[test]
    fn total_mechanical_energy_finite_after_steps() {
        let cfg = SphSimConfig::default_water();
        let mut sim = SphSim::new(cfg, 4);
        for _ in 0..3 {
            sim.step();
        }
        assert!(sim.total_mechanical_energy().is_finite());
    }
    #[test]
    fn linear_momentum_all_zero_initially() {
        let cfg = SphSimConfig::default_water();
        let sim = SphSim::new(cfg, 4);
        let p = sim.linear_momentum();
        for &pk in &p {
            assert_eq!(pk, 0.0);
        }
    }
    #[test]
    fn flat_positions_correct_length() {
        let cfg = SphSimConfig::default_water();
        let sim = SphSim::new(cfg, 5);
        let flat = sim.flat_positions();
        assert_eq!(flat.len(), 15, "5 particles × 3 coordinates = 15");
    }
    #[test]
    fn max_inter_particle_distance_positive_for_multiple_particles() {
        let cfg = SphSimConfig::default_water();
        let sim = SphSim::new(cfg, 4);
        assert!(sim.max_inter_particle_distance() > 0.0);
    }
    #[test]
    fn apply_floor_reflection_bounces_below_floor() {
        let cfg = SphSimConfig::default_water();
        let mut sim = SphSim::new(cfg, 0);
        sim.positions = vec![[0.5, -0.1, 0.5]];
        sim.velocities = vec![[0.0, -1.0, 0.0]];
        sim.masses = vec![1.0];
        sim.densities = vec![0.0];
        sim.pressures = vec![0.0];
        sim.forces = vec![[0.0; 3]];
        sim.apply_floor_reflection(0.0, 1.0);
        assert!(sim.positions[0][1] >= 0.0);
        assert!(sim.velocities[0][1] >= 0.0);
    }
    #[test]
    fn scale_masses_changes_all_masses() {
        let cfg = SphSimConfig::default_water();
        let mut sim = SphSim::new(cfg, 4);
        let original_total: f64 = sim.masses.iter().sum();
        sim.scale_masses(2.0);
        let new_total: f64 = sim.masses.iter().sum();
        assert!((new_total - 2.0 * original_total).abs() < 1e-10);
    }
    #[test]
    fn mean_particle_spacing_decreases_with_more_particles() {
        let cfg = SphSimConfig::default_water();
        let sim_coarse = SphSim::new(cfg.clone(), 8);
        let sim_fine = SphSim::new(cfg, 64);
        assert!(
            sim_fine.mean_particle_spacing() < sim_coarse.mean_particle_spacing(),
            "Finer grid should have smaller mean spacing"
        );
    }
    #[test]
    fn mean_pressure_initially_zero() {
        let cfg = SphSimConfig::default_water();
        let sim = SphSim::new(cfg, 4);
        assert_eq!(sim.mean_pressure(), 0.0);
    }
    #[test]
    fn count_negative_pressure_zero_initially() {
        let cfg = SphSimConfig::default_water();
        let sim = SphSim::new(cfg, 4);
        assert_eq!(sim.count_negative_pressure(), 0);
    }
    #[test]
    fn rescale_velocities_multiplies_by_factor() {
        let params = SphSimulationParams::default();
        let mut sim = SphSimulation::new(params);
        sim.add_particle_raw([0.5, 0.5, 0.5], [2.0, 0.0, 0.0], 1.0);
        sim.rescale_velocities(0.5);
        assert!((sim.particles.velocities[0].x - 1.0).abs() < 1e-14);
    }
    #[test]
    fn count_fast_particles_zero_when_stationary() {
        let params = SphSimulationParams::default();
        let mut sim = SphSimulation::new(params);
        sim.add_particle_raw([0.5, 0.5, 0.5], [0.0; 3], 1.0);
        assert_eq!(sim.count_fast_particles(0.1), 0);
    }
    #[test]
    fn count_fast_particles_detects_fast_ones() {
        let params = SphSimulationParams::default();
        let mut sim = SphSimulation::new(params);
        sim.add_particle_raw([0.5, 0.5, 0.5], [100.0, 0.0, 0.0], 1.0);
        assert_eq!(sim.count_fast_particles(50.0), 1);
    }
    #[test]
    fn density_std_zero_for_uniform_density() {
        let params = SphSimulationParams::default();
        let mut sim = SphSimulation::new(params);
        for _ in 0..4 {
            sim.add_particle_raw([0.1, 0.1, 0.1], [0.0; 3], 0.001);
        }
        for d in &mut sim.particles.densities {
            *d = 1000.0;
        }
        assert!(sim.density_std() < 1e-10);
    }
}
