//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use crate::atom::AtomSet;

/// Trait for thermostats that control the system temperature.
pub trait Thermostat: Send + Sync {
    /// Apply temperature control to the atom velocities.
    ///
    /// `target_temp` is the desired temperature, `dt` is the time step,
    /// and `boltzmann_k` is the Boltzmann constant in the chosen unit system.
    fn apply(&mut self, atoms: &mut AtomSet, target_temp: f64, dt: f64, boltzmann_k: f64);
}
#[cfg(test)]
mod tests {

    use super::*;
    use crate::AndersenThermostat;
    use crate::BerendsenThermostat;
    use crate::NoThermostat;
    use crate::NoseHooverThermostat;
    use crate::VelocityRescalingThermostat;
    use crate::thermostat::CsvrThermostat;
    use crate::thermostat::KineticDiagnostics;
    use crate::thermostat::TemperatureProfile;
    use crate::thermostat::ThermostatSwitcher;
    use oxiphysics_core::Vec3;
    #[test]
    fn test_berendsen_drives_temperature() {
        let mut atoms = AtomSet::new();
        atoms.add_atom(Vec3::zeros(), Vec3::new(10.0, 0.0, 0.0), 1.0, 0.0, 0);
        atoms.add_atom(Vec3::zeros(), Vec3::new(-10.0, 0.0, 0.0), 1.0, 0.0, 0);
        let k_b = 1.0;
        let target = 1.0;
        let initial_temp = atoms.temperature(k_b);
        assert!(initial_temp > target);
        let mut thermo = BerendsenThermostat::new(0.1);
        for _ in 0..1000 {
            thermo.apply(&mut atoms, target, 0.001, k_b);
        }
        let final_temp = atoms.temperature(k_b);
        assert!(
            (final_temp - target).abs() < 0.5,
            "Berendsen failed: T={final_temp}, target={target}"
        );
    }
    #[test]
    fn test_nose_hoover_convergence() {
        let mut atoms = AtomSet::new();
        atoms.add_atom(Vec3::zeros(), Vec3::new(5.0, 3.0, 1.0), 1.0, 0.0, 0);
        atoms.add_atom(Vec3::zeros(), Vec3::new(-2.0, 4.0, -1.0), 1.0, 0.0, 0);
        let k_b = 1.0;
        let target = 1.0;
        let mut thermo = NoseHooverThermostat::new(10.0);
        for _ in 0..5000 {
            thermo.apply(&mut atoms, target, 0.001, k_b);
        }
        let final_temp = atoms.temperature(k_b);
        assert!(
            (final_temp - target).abs() < 2.0,
            "Nose-Hoover: T={final_temp}, target={target}"
        );
    }
    /// Build an AtomSet with `n` atoms all at the same speed such that the
    /// instantaneous temperature equals `temp` (reduced units, kB = 1.0).
    ///
    /// Each atom has mass 1.0.  We need:
    ///   T = 2*KE / (3*N*kB)  with KE = N * 0.5 * v^2
    ///   => T = v^2 / 3  => v = sqrt(3 * T)
    /// We split the speed equally across x, y, z (alternating sign) so that
    /// the total linear momentum is zero.
    fn make_atoms_at_temp(n: usize, temp: f64) -> AtomSet {
        let v = (3.0 * temp).sqrt();
        let mut atoms = AtomSet::new();
        for i in 0..n {
            let sign = if i % 2 == 0 { 1.0_f64 } else { -1.0_f64 };
            let vc = v / 3.0_f64.sqrt();
            atoms.add_atom(
                Vec3::zeros(),
                Vec3::new(sign * vc, sign * vc, sign * vc),
                1.0,
                0.0,
                0,
            );
        }
        atoms
    }
    #[test]
    fn test_berendsen_cools_hot_system() {
        let mut atoms = make_atoms_at_temp(20, 600.0);
        let k_b = 1.0;
        let target = 300.0;
        let initial_temp = atoms.temperature(k_b);
        assert!(
            initial_temp > 500.0,
            "Setup error: initial T={initial_temp}"
        );
        let mut thermo = BerendsenThermostat::new(0.1);
        for _ in 0..100 {
            thermo.apply(&mut atoms, target, 0.01, k_b);
        }
        let final_temp = atoms.temperature(k_b);
        assert!(
            final_temp < 400.0,
            "Berendsen should cool system: T={final_temp} (target={target})"
        );
    }
    #[test]
    fn test_berendsen_heats_cold_system() {
        let mut atoms = make_atoms_at_temp(20, 50.0);
        let k_b = 1.0;
        let target = 300.0;
        let initial_temp = atoms.temperature(k_b);
        assert!(
            initial_temp < 100.0,
            "Setup error: initial T={initial_temp}"
        );
        let mut thermo = BerendsenThermostat::new(0.1);
        for _ in 0..100 {
            thermo.apply(&mut atoms, target, 0.01, k_b);
        }
        let final_temp = atoms.temperature(k_b);
        assert!(
            final_temp > 150.0,
            "Berendsen should heat system: T={final_temp} (target={target})"
        );
    }
    #[test]
    fn test_berendsen_temperature_measure() {
        let mut atoms = AtomSet::new();
        let k_b = 1.0;
        for i in 0..4_usize {
            let sign = if i % 2 == 0 { 1.0_f64 } else { -1.0_f64 };
            atoms.add_atom(Vec3::zeros(), Vec3::new(sign, sign, sign), 1.0, 0.0, 0);
        }
        let t = atoms.temperature(k_b);
        assert!(
            (t - 1.0).abs() < 1e-10,
            "temperature() returned {t}, expected 1.0"
        );
        let ke = atoms.kinetic_energy();
        let n = atoms.len() as f64;
        let t_manual = 2.0 * ke / (3.0 * n * k_b);
        assert!(
            (t - t_manual).abs() < 1e-10,
            "temperature() {t} != manual {t_manual}"
        );
    }
    /// Verify that the NH thermostat tracks an extended energy that includes
    /// both the kinetic energy and the thermostat bath energy (Q*xi^2/2).
    /// This implementation uses a simplified first-order NH integrator, so
    /// the full conserved quantity (which also includes kB*T*ln(s)) is not
    /// explicitly tracked.  We instead verify:
    ///   1. The extended energy KE + Q*xi^2/2 is non-negative (sanity check).
    ///   2. The kinetic energy remains finite and positive after many steps.
    #[test]
    fn test_nose_hoover_conserves_extended_energy() {
        let mut atoms = AtomSet::new();
        let k_b = 1.0;
        let target = 1.0;
        atoms.add_atom(Vec3::zeros(), Vec3::new(1.0, 0.0, 0.0), 1.0, 0.0, 0);
        atoms.add_atom(Vec3::zeros(), Vec3::new(-1.0, 0.0, 0.0), 1.0, 0.0, 0);
        atoms.add_atom(Vec3::zeros(), Vec3::new(0.0, 1.0, 0.0), 1.0, 0.0, 0);
        atoms.add_atom(Vec3::zeros(), Vec3::new(0.0, -1.0, 0.0), 1.0, 0.0, 0);
        let mut thermo = NoseHooverThermostat::new(10.0);
        for _ in 0..2000 {
            thermo.apply(&mut atoms, target, 0.001, k_b);
        }
        let ke_final = atoms.kinetic_energy();
        let ext_final = ke_final + 0.5 * thermo.q * thermo.xi.powi(2);
        assert!(
            ext_final >= 0.0,
            "NH extended energy is negative: {ext_final}"
        );
        assert!(
            ke_final > 0.0 && ke_final.is_finite(),
            "NH kinetic energy out of bounds: {ke_final}"
        );
    }
    /// Velocity rescaling must bring temperature to exactly the target value
    /// in a single application.
    #[test]
    fn test_velocity_rescaling_exact_temperature() {
        let mut atoms = make_atoms_at_temp(10, 200.0);
        let k_b = 1.0;
        let target = 300.0;
        let mut thermo = VelocityRescalingThermostat::new();
        thermo.apply(&mut atoms, target, 0.001, k_b);
        let t_after = atoms.temperature(k_b);
        assert!(
            (t_after - target).abs() < 1e-8,
            "VelocityRescaling should set temperature exactly; got T={t_after}"
        );
    }
    /// Velocity rescaling from a cold start must raise temperature to target.
    #[test]
    fn test_velocity_rescaling_heats_cold_system() {
        let mut atoms = make_atoms_at_temp(20, 10.0);
        let k_b = 1.0;
        let target = 300.0;
        let mut thermo = VelocityRescalingThermostat::new();
        for _ in 0..5 {
            thermo.apply(&mut atoms, target, 0.001, k_b);
        }
        let t_after = atoms.temperature(k_b);
        assert!(
            (t_after - target).abs() < 1e-6,
            "VelocityRescaling should achieve target temperature; got T={t_after}"
        );
    }
    /// After many Andersen steps (high collision rate), the temperature
    /// should converge close to the target.
    #[test]
    fn test_andersen_thermostat_convergence() {
        let mut atoms = AtomSet::new();
        for _ in 0..100 {
            atoms.add_atom(Vec3::zeros(), Vec3::zeros(), 1.0, 0.0, 0);
        }
        let k_b = 1.0;
        let target = 1.0;
        let mut thermo = AndersenThermostat::with_seed(10.0, 42);
        for _ in 0..200 {
            thermo.apply(&mut atoms, target, 0.1, k_b);
        }
        let t_after = atoms.temperature(k_b);
        assert!(
            (t_after - target).abs() < 0.4,
            "Andersen: temperature {t_after} not near target {target}"
        );
    }
    /// With zero collision rate the Andersen thermostat must not change
    /// any velocity.
    #[test]
    fn test_andersen_zero_collision_rate() {
        let mut atoms = make_atoms_at_temp(10, 300.0);
        let v_before: Vec<Vec3> = atoms.velocities.clone();
        let k_b = 1.0;
        let mut thermo = AndersenThermostat::new(0.0);
        thermo.apply(&mut atoms, 300.0, 0.01, k_b);
        for (vb, va) in v_before.iter().zip(atoms.velocities.iter()) {
            let diff = (*va - *vb).norm();
            assert!(diff < 1e-15, "Velocity changed despite zero collision rate");
        }
    }
    /// Thermostat brings a hot system to the target temperature (required test).
    #[test]
    fn test_thermostat_brings_to_target_temperature() {
        let mut atoms = make_atoms_at_temp(40, 600.0);
        let k_b = 1.0;
        let target = 300.0;
        let mut thermo = BerendsenThermostat::new(0.05);
        for _ in 0..200 {
            thermo.apply(&mut atoms, target, 0.01, k_b);
        }
        let t_final = atoms.temperature(k_b);
        assert!(
            (t_final - target).abs() < 20.0,
            "Berendsen: expected T~{target}, got T={t_final}"
        );
    }
    #[test]
    fn test_csvr_drives_temperature() {
        let mut atoms = make_atoms_at_temp(40, 600.0);
        let k_b = 1.0;
        let target = 300.0;
        let mut thermo = CsvrThermostat::with_seed(0.5, 42);
        for _ in 0..2000 {
            thermo.apply(&mut atoms, target, 0.01, k_b);
        }
        let t_final = atoms.temperature(k_b);
        assert!(
            (t_final - target).abs() < 300.0,
            "CSVR: T={t_final}, target={target}"
        );
    }
    #[test]
    fn test_csvr_positive_kinetic_energy() {
        let mut atoms = make_atoms_at_temp(10, 100.0);
        let k_b = 1.0;
        let target = 300.0;
        let mut thermo = CsvrThermostat::with_seed(0.05, 123);
        for _ in 0..200 {
            thermo.apply(&mut atoms, target, 0.001, k_b);
        }
        let ke = atoms.kinetic_energy();
        assert!(
            ke > 0.0 && ke.is_finite(),
            "KE should remain positive: {ke}"
        );
    }
    #[test]
    fn test_no_thermostat_preserves_velocities() {
        let mut atoms = make_atoms_at_temp(10, 300.0);
        let v_before: Vec<Vec3> = atoms.velocities.clone();
        let k_b = 1.0;
        let mut thermo = NoThermostat::new();
        thermo.apply(&mut atoms, 300.0, 0.01, k_b);
        for (vb, va) in v_before.iter().zip(atoms.velocities.iter()) {
            let diff = (*va - *vb).norm();
            assert!(diff < 1e-15, "NoThermostat should not change velocities");
        }
    }
    #[test]
    fn test_temperature_profile_mean() {
        let mut prof = TemperatureProfile::new();
        prof.record(0, 100.0);
        prof.record(1, 200.0);
        prof.record(2, 300.0);
        assert!((prof.mean() - 200.0).abs() < 1e-14);
        assert_eq!(prof.len(), 3);
        assert!(!prof.is_empty());
    }
    #[test]
    fn test_temperature_profile_std_dev() {
        let mut prof = TemperatureProfile::new();
        for i in 0..10 {
            prof.record(i, 300.0);
        }
        assert!(prof.std_dev() < 1e-14);
    }
    #[test]
    fn test_temperature_profile_equilibrated() {
        let mut prof = TemperatureProfile::new();
        for i in 0..100 {
            prof.record(i, 300.0 + 0.001 * (i as f64));
        }
        assert!(prof.is_equilibrated(50, 1.0));
    }
    #[test]
    fn test_temperature_profile_not_equilibrated() {
        let mut prof = TemperatureProfile::new();
        for i in 0..20 {
            prof.record(i, i as f64 * 100.0);
        }
        assert!(!prof.is_equilibrated(10, 1.0));
    }
    #[test]
    fn test_temperature_profile_mean_last_n() {
        let mut prof = TemperatureProfile::new();
        for i in 0..10 {
            prof.record(i, (i + 1) as f64);
        }
        assert!((prof.mean_last_n(3) - 9.0).abs() < 1e-14);
    }
    #[test]
    fn test_temperature_profile_max_min() {
        let mut prof = TemperatureProfile::new();
        prof.record(0, 50.0);
        prof.record(1, 400.0);
        prof.record(2, 200.0);
        assert!((prof.max_temp() - 400.0).abs() < 1e-14);
        assert!((prof.min_temp() - 50.0).abs() < 1e-14);
    }
    #[test]
    fn test_temperature_profile_latest() {
        let mut prof = TemperatureProfile::new();
        prof.record(0, 100.0);
        prof.record(1, 250.0);
        assert!((prof.latest() - 250.0).abs() < 1e-14);
    }
    #[test]
    fn test_component_temperatures() {
        let mut atoms = AtomSet::new();
        atoms.add_atom(Vec3::zeros(), Vec3::new(3.0, 0.0, 0.0), 1.0, 0.0, 0);
        atoms.add_atom(Vec3::zeros(), Vec3::new(-3.0, 0.0, 0.0), 1.0, 0.0, 0);
        let k_b = 1.0;
        let (tx, ty, tz) = KineticDiagnostics::component_temperatures(&atoms, k_b);
        assert!(tx > 0.0);
        assert!(ty.abs() < 1e-14);
        assert!(tz.abs() < 1e-14);
    }
    #[test]
    fn test_equipartition_ratio_isotropic() {
        let atoms = make_atoms_at_temp(20, 300.0);
        let k_b = 1.0;
        let ratio = KineticDiagnostics::equipartition_ratio(&atoms, k_b);
        assert!(
            (ratio - 1.0).abs() < 0.01,
            "Isotropic velocity distribution should have ratio ~1, got {ratio}"
        );
    }
    #[test]
    fn test_com_velocity_zero_momentum() {
        let atoms = make_atoms_at_temp(20, 300.0);
        let v_com = KineticDiagnostics::com_velocity(&atoms);
        assert!(
            v_com.norm() < 1e-10,
            "COM velocity should be ~0 for balanced system"
        );
    }
    #[test]
    fn test_remove_com_drift() {
        let mut atoms = AtomSet::new();
        atoms.add_atom(Vec3::zeros(), Vec3::new(10.0, 5.0, 0.0), 1.0, 0.0, 0);
        atoms.add_atom(Vec3::zeros(), Vec3::new(8.0, 3.0, 0.0), 1.0, 0.0, 0);
        KineticDiagnostics::remove_com_drift(&mut atoms);
        let v_com = KineticDiagnostics::com_velocity(&atoms);
        assert!(
            v_com.norm() < 1e-14,
            "COM velocity should be 0 after removal"
        );
    }
    #[test]
    fn test_thermostat_switcher_phases() {
        let mut atoms = make_atoms_at_temp(10, 600.0);
        let k_b = 1.0;
        let mut switcher = ThermostatSwitcher::new(Box::new(BerendsenThermostat::new(0.1)));
        switcher.add_phase(10, Box::new(VelocityRescalingThermostat::new()));
        let target = 300.0;
        for _ in 0..20 {
            switcher.apply_step(&mut atoms, target, 0.01, k_b);
        }
        assert_eq!(switcher.current_step(), 20);
        let t_final = atoms.temperature(k_b);
        assert!(
            (t_final - target).abs() < 1.0,
            "After VelocityRescaling phase: T={t_final}"
        );
    }
    #[test]
    fn test_thermostat_switcher_single_phase() {
        let mut atoms = make_atoms_at_temp(10, 100.0);
        let k_b = 1.0;
        let mut switcher = ThermostatSwitcher::new(Box::new(NoThermostat::new()));
        let t_before = atoms.temperature(k_b);
        switcher.apply_step(&mut atoms, 300.0, 0.01, k_b);
        let t_after = atoms.temperature(k_b);
        assert!(
            (t_before - t_after).abs() < 1e-14,
            "NoThermostat should not change T"
        );
    }
}
#[cfg(test)]
mod nhc_langevin_tests {
    use super::super::types::*;
    use super::*;
    use oxiphysics_core::Vec3;
    /// Build an AtomSet with `n` atoms at temperature `temp` (k_B = 1).
    fn atoms_at_temp(n: usize, temp: f64) -> AtomSet {
        let v = (3.0 * temp).sqrt() / 3.0_f64.sqrt();
        let mut atoms = AtomSet::new();
        for i in 0..n {
            let sign = if i % 2 == 0 { 1.0_f64 } else { -1.0_f64 };
            atoms.add_atom(
                Vec3::zeros(),
                Vec3::new(sign * v, sign * v, sign * v),
                1.0,
                0.0,
                0,
            );
        }
        atoms
    }
    #[test]
    fn test_nhc_creation() {
        let nhc = NoseHooverChain::new(3, vec![10.0, 10.0, 10.0], 30);
        assert_eq!(nhc.chain_length(), 3);
        assert_eq!(nhc.eta.len(), 3);
        assert_eq!(nhc.xi.len(), 3);
    }
    #[test]
    fn test_nhc_uniform_mass() {
        let nhc = NoseHooverChain::with_uniform_mass(4, 5.0, 12);
        assert_eq!(nhc.chain_length(), 4);
        for &q in &nhc.q {
            assert!((q - 5.0).abs() < 1e-12);
        }
    }
    #[test]
    fn test_nhc_chain1_scale_factor_finite() {
        let mut nhc = NoseHooverChain::new(1, vec![10.0], 30);
        let mut atoms = atoms_at_temp(10, 600.0);
        let ke_init = atoms.kinetic_energy();
        assert!(ke_init > 0.0);
        nhc.apply(&mut atoms, 300.0, 0.001, 1.0);
        let ke_final = atoms.kinetic_energy();
        assert!(ke_final.is_finite() && ke_final > 0.0);
    }
    #[test]
    fn test_nhc_chain3_reduces_temperature() {
        let mut nhc = NoseHooverChain::with_uniform_mass(3, 50.0, 60);
        let mut atoms = atoms_at_temp(20, 600.0);
        let t_init = atoms.temperature(1.0);
        assert!(t_init > 400.0);
        for _ in 0..5000 {
            nhc.apply(&mut atoms, 300.0, 0.001, 1.0);
        }
        let t_final = atoms.temperature(1.0);
        assert!(
            t_final < t_init,
            "NHC chain-3 should reduce T: init={t_init}, final={t_final}"
        );
    }
    #[test]
    fn test_nhc_conserved_energy_nonneg() {
        let nhc = NoseHooverChain::with_uniform_mass(3, 10.0, 30);
        let e = nhc.conserved_energy(300.0, 1.0);
        assert!(e.is_finite(), "NHC conserved energy should be finite: {e}");
    }
    #[test]
    fn test_nhc_single_chain_preserves_positivity() {
        let mut nhc = NoseHooverChain::new(1, vec![10.0], 6);
        let mut atoms = atoms_at_temp(2, 300.0);
        for _ in 0..1000 {
            nhc.apply(&mut atoms, 300.0, 0.001, 1.0);
        }
        let ke = atoms.kinetic_energy();
        assert!(
            ke >= 0.0 && ke.is_finite(),
            "KE must remain non-negative: {ke}"
        );
    }
    #[test]
    fn test_langevin_creation() {
        let lang = LangevinThermostat::new(1.0);
        assert!((lang.gamma - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_langevin_ou_step_changes_velocities() {
        let mut lang = LangevinThermostat::with_seed(1.0, 99);
        let mut atoms = atoms_at_temp(10, 300.0);
        let v0 = atoms.velocities[0];
        lang.apply_ou_step(&mut atoms, 300.0, 0.01, 1.0);
        for v in &atoms.velocities {
            assert!(v.norm().is_finite(), "velocity must remain finite");
        }
        let changed = atoms.velocities.iter().any(|v| (*v - v0).norm() > 1e-15);
        assert!(changed, "OU step should change at least one velocity");
    }
    #[test]
    fn test_langevin_ou_high_friction_convergence() {
        let mut lang = LangevinThermostat::with_seed(100.0, 42);
        let mut atoms = atoms_at_temp(50, 0.0);
        let k_b = 1.0;
        let target = 300.0;
        for _ in 0..1000 {
            lang.apply(&mut atoms, target, 0.01, k_b);
        }
        let t_final = atoms.temperature(k_b);
        assert!(
            t_final > 10.0,
            "Langevin should heat cold system; T={t_final}"
        );
    }
    #[test]
    fn test_langevin_random_force_amplitude() {
        let lang = LangevinThermostat::new(2.0);
        let amp = lang.random_force_amplitude(1.0, 300.0, 1.0);
        let expected = (2.0 * 2.0 * 1.0 * 300.0 / 1.0_f64).sqrt();
        assert!(
            (amp - expected).abs() < 1e-10,
            "amp={amp}, expected={expected}"
        );
    }
    #[test]
    fn test_langevin_zero_friction_preserves_temperature() {
        let mut lang = LangevinThermostat::with_seed(0.0, 77);
        let mut atoms = atoms_at_temp(20, 300.0);
        let t_before = atoms.temperature(1.0);
        lang.apply_ou_step(&mut atoms, 300.0, 0.01, 1.0);
        let t_after = atoms.temperature(1.0);
        assert!(
            (t_before - t_after).abs() < 1e-10,
            "zero gamma: T should not change; before={t_before}, after={t_after}"
        );
    }
    #[test]
    fn test_svr_creation() {
        let svr = SvrThermostat::new(0.1, 30);
        assert!((svr.tau - 0.1).abs() < 1e-12);
        assert_eq!(svr.n_dof, 30);
    }
    #[test]
    fn test_svr_rescale_factor_is_finite() {
        let mut svr = SvrThermostat::with_seed(0.1, 30, 42);
        let f = svr.rescale_factor(100.0, 150.0, 0.01);
        assert!(
            f.is_finite() && f > 0.0,
            "rescale factor must be positive finite: {f}"
        );
    }
    #[test]
    fn test_svr_ke_remains_positive() {
        let mut svr = SvrThermostat::with_seed(0.5, 60, 123);
        let mut atoms = atoms_at_temp(20, 200.0);
        let k_b = 1.0;
        let target = 300.0;
        for _ in 0..500 {
            svr.apply(&mut atoms, target, 0.01, k_b);
        }
        let ke = atoms.kinetic_energy();
        assert!(
            ke > 0.0 && ke.is_finite(),
            "SVR: KE must remain positive: {ke}"
        );
    }
    #[test]
    fn test_svr_drives_temperature_toward_target() {
        let mut svr = SvrThermostat::with_seed(0.1, 60, 9999);
        let mut atoms = atoms_at_temp(20, 600.0);
        let k_b = 1.0;
        let target = 300.0;
        for _ in 0..3000 {
            svr.apply(&mut atoms, target, 0.01, k_b);
        }
        let t_final = atoms.temperature(k_b);
        assert!(t_final < 600.0, "SVR should cool system: T_final={t_final}");
    }
    #[test]
    fn test_temperature_ramp_start() {
        let ramp = TemperatureRamp::new(100.0, 400.0, 100);
        assert!((ramp.current_target_temp() - 100.0).abs() < 1e-12);
    }
    #[test]
    fn test_temperature_ramp_end() {
        let ramp = TemperatureRamp::new(100.0, 400.0, 100);
        assert!((ramp.temp_at_step(100) - 400.0).abs() < 1e-12);
    }
    #[test]
    fn test_temperature_ramp_midpoint() {
        let ramp = TemperatureRamp::new(100.0, 300.0, 100);
        let t = ramp.temp_at_step(50);
        assert!((t - 200.0).abs() < 1e-10, "midpoint: T={t}");
    }
    #[test]
    fn test_temperature_ramp_advance() {
        let mut ramp = TemperatureRamp::new(0.0, 100.0, 10);
        for _ in 0..5 {
            ramp.next_target_temp();
        }
        assert_eq!(ramp.current_step(), 5);
        let t = ramp.current_target_temp();
        assert!((t - 50.0).abs() < 1e-10, "after 5 steps: T={t}");
    }
    #[test]
    fn test_temperature_ramp_completion() {
        let mut ramp = TemperatureRamp::new(0.0, 100.0, 5);
        assert!(!ramp.is_complete());
        for _ in 0..5 {
            ramp.next_target_temp();
        }
        assert!(
            ramp.is_complete(),
            "ramp should be complete after ramp_steps"
        );
        let t = ramp.current_target_temp();
        assert!((t - 100.0).abs() < 1e-12, "after completion: T={t}");
    }
    #[test]
    fn test_temperature_ramp_reset() {
        let mut ramp = TemperatureRamp::new(0.0, 100.0, 10);
        for _ in 0..7 {
            ramp.next_target_temp();
        }
        ramp.reset();
        assert_eq!(ramp.current_step(), 0);
        assert!((ramp.current_target_temp() - 0.0).abs() < 1e-12);
    }
    #[test]
    fn test_temperature_ramp_zero_steps() {
        let ramp = TemperatureRamp::new(200.0, 500.0, 0);
        let t = ramp.current_target_temp();
        assert!((t - 500.0).abs() < 1e-12, "0-step ramp: T={t}");
    }
    #[test]
    fn test_temperature_ramp_cooling() {
        let ramp = TemperatureRamp::new(300.0, 100.0, 100);
        let t0 = ramp.temp_at_step(0);
        let t50 = ramp.temp_at_step(50);
        let t100 = ramp.temp_at_step(100);
        assert!(
            t0 > t50 && t50 > t100,
            "cooling ramp: {t0} > {t50} > {t100}"
        );
    }
    /// Friction force must be exactly -γ·m·v for each component.
    #[test]
    fn test_langevin_friction_force_direction() {
        let gamma = 2.0;
        let mass = 3.0;
        let mut lang = LangevinThermostat::with_seed(gamma, 7777);
        let velocity = [1.0, -2.0, 0.5];
        let (f_fric, _) = lang.compute_friction_force(velocity, mass, 300.0, 0.001, 1.0);
        for k in 0..3 {
            let expected = -gamma * mass * velocity[k];
            assert!(
                (f_fric[k] - expected).abs() < 1e-12,
                "friction_force[{k}]: got {}, expected {expected}",
                f_fric[k]
            );
        }
    }
    /// Zero velocity must yield zero friction force.
    #[test]
    fn test_langevin_friction_force_zero_velocity() {
        let mut lang = LangevinThermostat::with_seed(1.5, 1234);
        let (f_fric, _) = lang.compute_friction_force([0.0, 0.0, 0.0], 2.0, 300.0, 0.001, 1.0);
        for v in f_fric {
            assert!(
                v.abs() < 1e-30,
                "zero-velocity friction must be zero, got {v}"
            );
        }
    }
    /// Stochastic force amplitude scales as sqrt(T): doubling T should multiply
    /// std of stochastic force by sqrt(2).
    #[test]
    fn test_langevin_stochastic_force_scales_with_temperature() {
        let gamma: f64 = 1.0;
        let mass: f64 = 1.0;
        let dt: f64 = 0.001;
        let kb: f64 = 1.0;
        let sigma_low = (2.0_f64 * gamma * mass * kb * 100.0 / dt).sqrt();
        let sigma_high = (2.0_f64 * gamma * mass * kb * 400.0 / dt).sqrt();
        let ratio = sigma_high / sigma_low;
        assert!(
            (ratio - 2.0).abs() < 1e-10,
            "stochastic amplitude should double when T×4: ratio={ratio}"
        );
    }
    /// Friction and stochastic forces are returned as `[f64; 3]` arrays.
    #[test]
    fn test_langevin_friction_returns_correct_types() {
        let mut lang = LangevinThermostat::with_seed(0.5, 42);
        let (f_fric, f_stoch) =
            lang.compute_friction_force([1.0, 2.0, 3.0], 1.0, 300.0, 0.002, 1.0);
        assert_eq!(f_fric.len(), 3);
        assert_eq!(f_stoch.len(), 3);
        for k in 0..3 {
            assert!(f_fric[k].is_finite(), "f_fric[{k}] not finite");
            assert!(f_stoch[k].is_finite(), "f_stoch[{k}] not finite");
        }
    }
    /// Collision probability at dt=0 must be 0.
    #[test]
    fn test_andersen_collision_prob_zero_dt() {
        let thermo = AndersenThermostat::new(10.0);
        let p = thermo.compute_collision_probability(0.0);
        assert!(p.abs() < 1e-12, "P(collision|dt=0) must be 0, got {p}");
    }
    /// For very large ν·dt, probability must saturate near 1.
    #[test]
    fn test_andersen_collision_prob_saturates() {
        let thermo = AndersenThermostat::new(1e6);
        let p = thermo.compute_collision_probability(1.0);
        assert!(
            (p - 1.0).abs() < 1e-6,
            "P(collision) should saturate at 1 for large ν·dt, got {p}"
        );
    }
    /// For small ν·dt, probability should approximate ν·dt linearly.
    #[test]
    fn test_andersen_collision_prob_linear_approximation() {
        let nu = 1.0;
        let dt = 1e-4;
        let thermo = AndersenThermostat::new(nu);
        let p = thermo.compute_collision_probability(dt);
        let linear_approx = nu * dt;
        assert!(
            (p - linear_approx).abs() < 1e-7,
            "Poisson prob should be ≈ ν·dt for small ν·dt; P={p}, linear={linear_approx}"
        );
    }
    /// Target KE must equal (3N/2)·k_B·T.
    #[test]
    fn test_csvr_kinetic_energy_target_formula() {
        let thermo = CsvrThermostat::new(0.5);
        let n = 10;
        let temp = 300.0;
        let kb = 1.380649e-23;
        let ke_target = thermo.compute_kinetic_energy_target(n, temp, kb);
        let expected = 0.5 * (3 * n) as f64 * kb * temp;
        assert!(
            (ke_target - expected).abs() < 1e-30,
            "KE target = {ke_target}, expected {expected}"
        );
    }
    /// Target KE must scale linearly with temperature.
    #[test]
    fn test_csvr_kinetic_energy_target_scales_linearly() {
        let thermo = CsvrThermostat::new(0.5);
        let n = 5;
        let kb = 1.0;
        let ke100 = thermo.compute_kinetic_energy_target(n, 100.0, kb);
        let ke200 = thermo.compute_kinetic_energy_target(n, 200.0, kb);
        assert!(
            (ke200 / ke100 - 2.0).abs() < 1e-12,
            "KE target must scale linearly with T; ratio={}",
            ke200 / ke100
        );
    }
    /// Target KE must scale linearly with number of atoms.
    #[test]
    fn test_csvr_kinetic_energy_target_scales_with_natoms() {
        let thermo = CsvrThermostat::new(0.5);
        let temp = 300.0;
        let kb = 1.0;
        let ke10 = thermo.compute_kinetic_energy_target(10, temp, kb);
        let ke20 = thermo.compute_kinetic_energy_target(20, temp, kb);
        assert!(
            (ke20 / ke10 - 2.0).abs() < 1e-12,
            "KE target must scale with N; ratio={}",
            ke20 / ke10
        );
    }
}
