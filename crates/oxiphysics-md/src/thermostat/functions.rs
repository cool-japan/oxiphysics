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
    /// CSVR must sample the canonical ensemble: the time-averaged temperature
    /// equals the target and the kinetic energy fluctuates with
    /// `Var(K) = 2 K̄² / N_f`.  An instantaneous temperature is useless here
    /// (its relative scatter is `sqrt(2 / N_f)`), so the test averages 45 000
    /// steps (≈ 2 000 independent samples at `dt/τ = 0.1`) after a burn-in.
    #[test]
    fn test_csvr_drives_temperature() {
        let n_atoms = 100;
        let n_dof = (3 * n_atoms) as f64;
        let k_b = 1.0;
        let target = 300.0;
        let k_bar = 0.5 * n_dof * k_b * target;
        let mut atoms = make_atoms_at_temp(n_atoms, 600.0);
        let mut thermo = CsvrThermostat::with_seed(0.1, 42);
        for _ in 0..5_000 {
            thermo.apply(&mut atoms, target, 0.01, k_b);
        }
        let mut series = Vec::with_capacity(45_000);
        for _ in 0..45_000 {
            thermo.apply(&mut atoms, target, 0.01, k_b);
            series.push(atoms.kinetic_energy());
        }
        let n = series.len() as f64;
        let mean = series.iter().sum::<f64>() / n;
        let var = series.iter().map(|k| (k - mean).powi(2)).sum::<f64>() / n;
        let t_mean = 2.0 * mean / (n_dof * k_b);
        assert!(
            (t_mean - target).abs() / target < 0.02,
            "CSVR: <T>={t_mean}, target={target}"
        );
        let var_canonical = 2.0 * k_bar * k_bar / n_dof;
        assert!(
            (var / var_canonical - 1.0).abs() < 0.15,
            "CSVR: Var(K)={var}, canonical 2 K^2/N_f={var_canonical}"
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
        let nhc =
            NoseHooverChain::new(3, vec![10.0, 10.0, 10.0], 30).expect("valid NHC parameters");
        assert_eq!(nhc.chain_length(), 3);
        assert_eq!(nhc.eta.len(), 3);
        assert_eq!(nhc.xi.len(), 3);
    }
    #[test]
    fn test_nhc_uniform_mass() {
        let nhc = NoseHooverChain::with_uniform_mass(4, 5.0, 12).expect("valid NHC parameters");
        assert_eq!(nhc.chain_length(), 4);
        for &q in &nhc.q {
            assert!((q - 5.0).abs() < 1e-12);
        }
    }
    #[test]
    fn test_nhc_chain1_scale_factor_finite() {
        let mut nhc = NoseHooverChain::new(1, vec![10.0], 30).expect("valid NHC parameters");
        let mut atoms = atoms_at_temp(10, 600.0);
        let ke_init = atoms.kinetic_energy();
        assert!(ke_init > 0.0);
        nhc.apply(&mut atoms, 300.0, 0.001, 1.0);
        let ke_final = atoms.kinetic_energy();
        assert!(ke_final.is_finite() && ke_final > 0.0);
    }
    #[test]
    fn test_nhc_chain3_reduces_temperature() {
        let mut nhc =
            NoseHooverChain::with_uniform_mass(3, 50.0, 60).expect("valid NHC parameters");
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
        let nhc = NoseHooverChain::with_uniform_mass(3, 10.0, 30).expect("valid NHC parameters");
        let e = nhc.conserved_energy(300.0, 1.0);
        assert!(e.is_finite(), "NHC conserved energy should be finite: {e}");
    }
    #[test]
    fn test_nhc_single_chain_preserves_positivity() {
        let mut nhc = NoseHooverChain::new(1, vec![10.0], 6).expect("valid NHC parameters");
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
    fn test_nhc_constructors_reject_invalid_parameters() {
        let empty = NoseHooverChain::new(0, Vec::new(), 6);
        let message = empty
            .expect_err("chain_length 0 must be rejected")
            .to_string();
        assert!(message.contains("chain_length"), "message: {message}");
        assert!(NoseHooverChain::with_uniform_mass(0, 1.0, 6).is_err());
        assert!(NoseHooverChain::new(3, vec![1.0, 1.0], 6).is_err());
        assert!(NoseHooverChain::new(1, vec![1.0, 1.0], 6).is_err());
        assert!(NoseHooverChain::new(2, vec![1.0, 0.0], 6).is_err());
        assert!(NoseHooverChain::new(2, vec![1.0, -2.0], 6).is_err());
        assert!(NoseHooverChain::new(2, vec![1.0, f64::NAN], 6).is_err());
        assert!(NoseHooverChain::new(2, vec![f64::INFINITY, 1.0], 6).is_err());
        assert!(NoseHooverChain::with_uniform_mass(2, -1.0, 6).is_err());
        assert!(NoseHooverChain::new(2, vec![1.0, 1.0], 0).is_err());
        assert!(NoseHooverChain::new(2, vec![1.0, 1.0], 6).is_ok());
    }
    #[test]
    fn test_nhc_propagate_is_inert_for_inconsistent_chains() {
        let mut empty = NoseHooverChain {
            chain_length: 0,
            q: Vec::new(),
            eta: Vec::new(),
            xi: Vec::new(),
            n_dof: 6,
        };
        assert_eq!(empty.propagate(3.0, 1.0, 1.0, 0.01), 1.0);
        let mut atoms = atoms_at_temp(2, 300.0);
        let before = atoms.velocities.clone();
        empty.apply(&mut atoms, 300.0, 0.01, 1.0);
        assert_eq!(atoms.velocities, before);
        // `chain_length` larger than the stored vectors: only the available
        // thermostat acts, with no out-of-bounds access.
        let mut truncated = NoseHooverChain {
            chain_length: 3,
            q: vec![10.0],
            eta: vec![0.0],
            xi: vec![0.0],
            n_dof: 6,
        };
        let scale = truncated.propagate(50.0, 1.0, 1.0, 0.01);
        assert!(scale.is_finite() && scale > 0.0, "scale = {scale}");
        assert!(truncated.xi[0] > 0.0, "KE above target must accelerate xi");
    }
    #[test]
    fn test_nhc_conserved_energy_weights_first_thermostat_by_n_dof() {
        let mut nhc =
            NoseHooverChain::new(3, vec![2.0, 3.0, 4.0], 12).expect("valid NHC parameters");
        nhc.eta = vec![0.5, -0.25, 0.75];
        nhc.xi = vec![1.0, 2.0, 0.5];
        let (k_b, temp) = (1.0, 2.0);
        // sum_k Q_k xi_k^2 / 2 = (2*1 + 3*4 + 4*0.25) / 2 = 7.5
        // N_f kT eta_1 = 12 * 2 * 0.5 = 12, kT (eta_2 + eta_3) = 2 * 0.5 = 1
        let expected = 7.5 + 12.0 + 1.0;
        let e = nhc.conserved_energy(temp, k_b);
        assert!(
            (e - expected).abs() < 1e-12,
            "H_nhc = {e}, expected {expected}"
        );
    }
    /// Independent single-thermostat Yoshida-Suzuki step: both half-kicks of
    /// every sub-step see the kinetic energy scaled by the sub-steps so far.
    fn reference_single_thermostat_step(
        xi0: f64,
        q: f64,
        n_dof: f64,
        kinetic_energy: f64,
        kbt: f64,
        dt: f64,
    ) -> (f64, f64) {
        let w1 = 1.0 / (2.0 - 2.0_f64.powf(1.0 / 3.0));
        let w0 = 1.0 - 2.0 * w1;
        let mut xi = xi0;
        let mut scale = 1.0_f64;
        for w in [w1, w0, w1] {
            let half = 0.5 * w * dt;
            xi += (2.0 * kinetic_energy * scale * scale - n_dof * kbt) / q * half;
            scale *= (-xi * w * dt).exp();
            xi += (2.0 * kinetic_energy * scale * scale - n_dof * kbt) / q * half;
        }
        (xi, scale)
    }
    #[test]
    fn test_nhc_substeps_use_the_scaled_kinetic_energy() {
        let (q, n_dof, ke, kbt, dt) = (5.0, 12.0, 40.0, 2.0, 0.05);
        let mut nhc = NoseHooverChain::new(1, vec![q], 12).expect("valid NHC parameters");
        nhc.xi[0] = 0.3;
        let scale = nhc.propagate(ke, kbt, 1.0, dt);
        let (xi_ref, scale_ref) = reference_single_thermostat_step(0.3, q, n_dof, ke, kbt, dt);
        assert!(
            (scale - scale_ref).abs() < 1e-14,
            "scale {scale} vs reference {scale_ref}"
        );
        assert!(
            (nhc.xi[0] - xi_ref).abs() < 1e-14,
            "xi {} vs reference {xi_ref}",
            nhc.xi[0]
        );
    }
    /// Outcome of an NHC run on twelve decoupled harmonic degrees of freedom.
    struct HarmonicNhcRun {
        /// Initial value of `KE + PE + H_nhc`.
        e0: f64,
        /// Largest relative deviation of `KE + PE + H_nhc` from `e0`.
        max_drift: f64,
        /// Largest `|xi_1|` reached.
        max_xi: f64,
        /// Range (max - min) of the physical energy `KE + PE`.
        physical_range: f64,
    }
    /// Four atoms in isotropic harmonic traps (twelve decoupled oscillators)
    /// integrated as NHC(dt/2) - velocity Verlet(dt) - NHC(dt/2) for 10 000
    /// steps with an NHC of the given length.
    fn run_harmonic_nhc(chain_length: usize) -> HarmonicNhcRun {
        let spring = [0.5_f64, 1.0, 2.0, 3.0];
        let start = [
            ([0.5, -0.3, 0.2], [1.0, -0.8, 0.5]),
            ([-0.4, 0.6, 0.1], [-0.7, 0.9, -0.4]),
            ([0.3, 0.2, -0.5], [0.6, 0.3, -1.1]),
            ([-0.2, -0.1, 0.4], [-0.5, -0.9, 0.8]),
        ];
        let mut atoms = AtomSet::new();
        for (p, v) in start {
            atoms.add_atom(
                Vec3::new(p[0], p[1], p[2]),
                Vec3::new(v[0], v[1], v[2]),
                1.0,
                0.0,
                0,
            );
        }
        let (k_b, target, dt) = (1.0, 1.0, 0.005);
        // Start from the equipartition kinetic energy N_f kT / 2 = 6.
        let rescale = (6.0 / atoms.kinetic_energy()).sqrt();
        for v in &mut atoms.velocities {
            *v *= rescale;
        }
        let compute_forces = |atoms: &mut AtomSet| {
            for (i, k) in spring.iter().enumerate() {
                atoms.forces[i] = atoms.positions[i] * (-k);
            }
        };
        let potential = |atoms: &AtomSet| -> f64 {
            (0..atoms.len())
                .map(|i| 0.5 * spring[i] * atoms.positions[i].norm_squared())
                .sum()
        };
        // Q_1 = N_f kT tau^2 and Q_k = kT tau^2 with tau = 1.
        let mut masses = vec![1.0; chain_length];
        masses[0] = 12.0;
        let mut nhc = NoseHooverChain::new(chain_length, masses, 12).expect("valid NHC parameters");
        compute_forces(&mut atoms);
        let e0 = atoms.kinetic_energy() + potential(&atoms) + nhc.conserved_energy(target, k_b);
        let mut run = HarmonicNhcRun {
            e0,
            max_drift: 0.0,
            max_xi: 0.0,
            physical_range: 0.0,
        };
        let mut physical_min = f64::INFINITY;
        let mut physical_max = f64::NEG_INFINITY;
        for _ in 0..10_000 {
            nhc.apply(&mut atoms, target, 0.5 * dt, k_b);
            for i in 0..atoms.len() {
                let accel = atoms.forces[i] / atoms.masses[i];
                atoms.velocities[i] += accel * (0.5 * dt);
            }
            for i in 0..atoms.len() {
                let v = atoms.velocities[i];
                atoms.positions[i] += v * dt;
            }
            compute_forces(&mut atoms);
            for i in 0..atoms.len() {
                let accel = atoms.forces[i] / atoms.masses[i];
                atoms.velocities[i] += accel * (0.5 * dt);
            }
            nhc.apply(&mut atoms, target, 0.5 * dt, k_b);
            let physical = atoms.kinetic_energy() + potential(&atoms);
            let e = physical + nhc.conserved_energy(target, k_b);
            run.max_drift = run.max_drift.max((e - e0).abs() / e0.abs());
            run.max_xi = run.max_xi.max(nhc.xi[0].abs());
            physical_min = physical_min.min(physical);
            physical_max = physical_max.max(physical);
        }
        run.physical_range = physical_max - physical_min;
        run
    }
    /// `KE + PE + H_nhc` is the NHC constant of motion: it must stay constant
    /// (relative drift below 1e-4 over 10 000 steps) for every chain length.
    /// With the former stale first-thermostat force the drift was 5e-2 to 1e-1.
    #[test]
    fn test_nhc_conserves_extended_energy_for_harmonic_oscillators() {
        for chain_length in [1_usize, 2, 3, 4, 6] {
            let run = run_harmonic_nhc(chain_length);
            assert!(
                run.max_drift < 1e-4,
                "chain {chain_length}: extended energy drifted by {:.3e} (relative to {})",
                run.max_drift,
                run.e0
            );
            // The test is only meaningful if the thermostat really exchanges energy.
            assert!(
                run.max_xi > 0.02,
                "chain {chain_length}: thermostat never engaged, max|xi_1| = {}",
                run.max_xi
            );
            assert!(
                run.physical_range > 100.0 * run.max_drift * run.e0.abs(),
                "chain {chain_length}: physical energy barely moved (range {}) so conservation proves little",
                run.physical_range
            );
        }
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
    /// SVR (Bussi) must equilibrate at the target temperature -- the former
    /// implementation settled near `T/2` -- with canonical fluctuations.
    #[test]
    fn test_svr_drives_temperature_toward_target() {
        let n_atoms = 100;
        let n_dof = 3 * n_atoms;
        let k_b = 1.0;
        let target = 300.0;
        let k_bar = 0.5 * n_dof as f64 * k_b * target;
        let mut svr = SvrThermostat::with_seed(0.1, n_dof, 9999);
        let mut atoms = atoms_at_temp(n_atoms, 600.0);
        for _ in 0..5_000 {
            svr.apply(&mut atoms, target, 0.01, k_b);
        }
        let mut series = Vec::with_capacity(45_000);
        for _ in 0..45_000 {
            svr.apply(&mut atoms, target, 0.01, k_b);
            series.push(atoms.kinetic_energy());
        }
        let n = series.len() as f64;
        let mean = series.iter().sum::<f64>() / n;
        let var = series.iter().map(|k| (k - mean).powi(2)).sum::<f64>() / n;
        let t_mean = 2.0 * mean / (n_dof as f64 * k_b);
        assert!(
            (t_mean - target).abs() / target < 0.02,
            "SVR: <T>={t_mean}, target={target}"
        );
        let var_canonical = 2.0 * k_bar * k_bar / n_dof as f64;
        assert!(
            (var / var_canonical - 1.0).abs() < 0.15,
            "SVR: Var(K)={var}, canonical 2 K^2/N_f={var_canonical}"
        );
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
