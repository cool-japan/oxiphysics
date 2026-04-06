//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use std::f64::consts::PI;

/// Elementary charge (C).
pub(super) const ELEM_CHARGE: f64 = 1.602_176_634e-19;
/// Vacuum permittivity (F/m).
pub(super) const EPS0: f64 = 8.854_187_817e-12;
/// Avogadro's number (mol⁻¹).
pub(super) const AVOGADRO: f64 = 6.022_140_76e23;
/// Boltzmann constant (J/K).
pub(super) const KB: f64 = 1.380_649e-23;
/// Faraday constant (C/mol).
pub(super) const FARADAY: f64 = 96_485.332_123;
/// Conversion factor: Joules per particle → kJ/mol.
pub(super) const J_TO_KJMOL: f64 = AVOGADRO / 1000.0;
/// Gas constant R = k_B * N_A (J mol⁻¹ K⁻¹).
pub(super) const R_GAS: f64 = 8.314_462_618;
/// kBT at 298.15 K in kJ/mol.
pub(super) const KBT_298: f64 = R_GAS * 298.15 / 1000.0;
/// Compute the Born solvation free energy ΔG_Born (kJ/mol).
///
/// ```text
/// ΔG = − (q² N_A) / (8π ε₀ r) · (1/ε_in − 1/ε_out)
/// ```
///
/// # Arguments
/// * `q`       – ionic charge (elementary charge units)
/// * `r`       – Born radius (m)
/// * `eps_in`  – inner dielectric constant
/// * `eps_out` – outer (solvent) dielectric constant
pub fn born_energy(q: f64, r: f64, eps_in: f64, eps_out: f64) -> f64 {
    if r == 0.0 || eps_in == 0.0 || eps_out == 0.0 {
        return 0.0;
    }
    let q_coulomb = q * ELEM_CHARGE;
    let prefactor = (q_coulomb * q_coulomb) / (8.0 * PI * EPS0 * r);
    -prefactor * (1.0 / eps_in - 1.0 / eps_out) * J_TO_KJMOL
}
/// Compute the Debye screening length λ_D (m).
///
/// ```text
/// λ_D = sqrt( ε₀ ε_r k_B T / (2 N_A e² I) )
/// ```
///
/// # Arguments
/// * `eps`            – relative dielectric constant
/// * `ionic_strength` – ionic strength I (mol/m³)
/// * `t`              – temperature (K)
pub fn debye_length(eps: f64, ionic_strength: f64, t: f64) -> f64 {
    if ionic_strength == 0.0 {
        return f64::INFINITY;
    }
    let numerator = EPS0 * eps * KB * t;
    let denominator = 2.0 * AVOGADRO * ELEM_CHARGE * ELEM_CHARGE * ionic_strength;
    (numerator / denominator).sqrt()
}
/// Compute the solvent-accessible surface area (SASA) for a sphere.
///
/// SASA = 4π (r_solute + r_probe)²
///
/// # Arguments
/// * `radius`       – solute radius (Å or any unit)
/// * `probe_radius` – probe radius (typically 1.4 Å for water)
pub fn solvent_accessible_area(radius: f64, probe_radius: f64) -> f64 {
    let r_eff = radius + probe_radius;
    4.0 * PI * r_eff * r_eff
}
/// Compute the Lennard-Jones dispersion contribution to solvation (kJ/mol).
///
/// Uses a simple London-type approximation: G_disp = −ε_disp · (σ/r)^6
/// integrated over a spherical cavity of radius r.
///
/// # Arguments
/// * `eps_disp` – dispersion energy scale (kJ/mol)
/// * `sigma`    – collision diameter (Å)
/// * `r`        – cavity radius (Å)
pub fn dispersion_energy(eps_disp: f64, sigma: f64, r: f64) -> f64 {
    if r <= 0.0 {
        return 0.0;
    }
    -eps_disp * (sigma / r).powi(6)
}
/// Compute the repulsion energy (kJ/mol) at contact.
///
/// G_rep = ε_rep · (σ/r)^12
///
/// # Arguments
/// * `eps_rep` – repulsion scale (kJ/mol)
/// * `sigma`   – collision diameter (Å)
/// * `r`       – contact distance (Å)
pub fn repulsion_energy(eps_rep: f64, sigma: f64, r: f64) -> f64 {
    if r <= 0.0 {
        return f64::INFINITY;
    }
    eps_rep * (sigma / r).powi(12)
}
#[cfg(test)]
mod tests {
    use super::super::types::*;
    use super::*;
    #[test]
    fn test_born_energy_negative_for_ion_transfer() {
        let dg = born_energy(1.0, 1.9e-10, 1.0, 78.5);
        assert!(dg < 0.0);
    }
    #[test]
    fn test_born_energy_zero_radius() {
        assert_eq!(born_energy(1.0, 0.0, 1.0, 78.5), 0.0);
    }
    #[test]
    fn test_born_energy_equal_dielectrics_zero() {
        let dg = born_energy(1.0, 1e-10, 78.5, 78.5);
        assert!(dg.abs() < 1e-6);
    }
    #[test]
    fn test_born_energy_divalent_quadruple() {
        let dg1 = born_energy(1.0, 1e-10, 1.0, 78.5);
        let dg2 = born_energy(2.0, 1e-10, 1.0, 78.5);
        assert!((dg2 / dg1 - 4.0).abs() < 1e-6);
    }
    #[test]
    fn test_debye_length_zero_ionic_inf() {
        assert!(debye_length(78.5, 0.0, 298.15).is_infinite());
    }
    #[test]
    fn test_debye_length_physiological() {
        let ld = debye_length(78.5, 150.0, 298.15);
        assert!(ld > 5e-10 && ld < 2e-9, "Debye ≈ 0.8 nm, got {ld:.2e}");
    }
    #[test]
    fn test_sasa_sphere() {
        let sasa = solvent_accessible_area(1.0, 1.4);
        let expected = 4.0 * PI * 2.4 * 2.4;
        assert!((sasa - expected).abs() < 1e-8);
    }
    #[test]
    fn test_dispersion_energy_negative() {
        let e = dispersion_energy(1.0, 1.0, 1.0);
        assert!(e <= 0.0);
    }
    #[test]
    fn test_dispersion_energy_zero_r() {
        assert_eq!(dispersion_energy(1.0, 1.0, 0.0), 0.0);
    }
    #[test]
    fn test_repulsion_energy_positive() {
        assert!(repulsion_energy(1.0, 1.0, 1.0) > 0.0);
    }
    #[test]
    fn test_repulsion_energy_zero_r_infinite() {
        assert!(repulsion_energy(1.0, 1.0, 0.0).is_infinite());
    }
    #[test]
    fn test_born_solvation_negative_for_water() {
        let b = BornSolvation::new(1.0, 1.9e-10, 78.5, 298.15);
        assert!(b.solvation_energy() < 0.0);
    }
    #[test]
    fn test_born_solvation_self_energy_vacuum_positive() {
        let b = BornSolvation::new(1.0, 1.9e-10, 78.5, 298.15);
        assert!(b.self_energy_vacuum() > 0.0);
    }
    #[test]
    fn test_born_solvation_self_energy_solvent_less_than_vacuum() {
        let b = BornSolvation::new(1.0, 1.9e-10, 78.5, 298.15);
        assert!(b.self_energy_solvent() < b.self_energy_vacuum());
    }
    #[test]
    fn test_born_solvation_scales_quadratically_with_charge() {
        let b1 = BornSolvation::new(1.0, 1.9e-10, 78.5, 298.15);
        let b2 = BornSolvation::new(2.0, 1.9e-10, 78.5, 298.15);
        let ratio = b2.solvation_energy() / b1.solvation_energy();
        assert!((ratio - 4.0).abs() < 1e-5);
    }
    #[test]
    fn test_born_solvation_entropy_zero() {
        let b = BornSolvation::new(1.0, 1.9e-10, 78.5, 298.15);
        assert_eq!(b.entropy_contribution(), 0.0);
    }
    #[test]
    fn test_born_solvation_scaled_energy() {
        let b = BornSolvation::new(1.0, 1.9e-10, 78.5, 298.15);
        let e = b.solvation_energy();
        assert!((b.scaled_energy(2.0) - 2.0 * e).abs() < 1e-10 * e.abs());
    }
    #[test]
    fn test_gb_single_atom_approximates_born() {
        let charges = vec![1.0_f64];
        let positions = vec![0.0_f64, 0.0, 0.0];
        let intrinsic_radii = vec![1.9e-10_f64];
        let scaling_factors = vec![0.85_f64];
        let gb = GeneralizedBorn::new(charges, positions, intrinsic_radii, scaling_factors, 78.5);
        let e_gb = gb.solvation_energy();
        let e_born = born_energy(1.0, 1.9e-10, 1.0, 78.5);
        assert!(e_gb < 0.0);
        assert!((e_gb / e_born - 1.0).abs() < 0.5);
    }
    #[test]
    fn test_gb_two_atoms_negative_energy() {
        let charges = vec![0.5_f64, -0.5];
        let positions = vec![0.0_f64, 0.0, 0.0, 3e-10, 0.0, 0.0];
        let radii = vec![1.5e-10_f64, 1.5e-10];
        let scales = vec![0.85_f64, 0.85];
        let gb = GeneralizedBorn::new(charges, positions, radii, scales, 78.5);
        let e = gb.solvation_energy();
        assert!(e.is_finite());
    }
    #[test]
    fn test_gb_function_rij_zero_equals_sqrt_alpha() {
        let alpha_i = 2e-10_f64;
        let alpha_j = 2e-10_f64;
        let fgb = GeneralizedBorn::gb_function(0.0, alpha_i, alpha_j);
        assert!((fgb - (alpha_i * alpha_j).sqrt()).abs() < 1e-20);
    }
    #[test]
    fn test_gb_effective_born_radius_positive() {
        let charges = vec![1.0_f64];
        let positions = vec![0.0_f64, 0.0, 0.0];
        let radii = vec![1.9e-10_f64];
        let scales = vec![0.85_f64];
        let gb = GeneralizedBorn::new(charges, positions, radii, scales, 78.5);
        assert!(gb.effective_born_radius(0) > 0.0);
    }
    #[test]
    fn test_gb_distance_formula() {
        let charges = vec![0.0_f64, 0.0];
        let positions = vec![0.0_f64, 0.0, 0.0, 3.0, 4.0, 0.0];
        let radii = vec![1.5e-10_f64, 1.5e-10];
        let scales = vec![0.85_f64, 0.85];
        let gb = GeneralizedBorn::new(charges, positions, radii, scales, 78.5);
        assert!((gb.distance(0, 1) - 5.0).abs() < 1e-12);
    }
    #[test]
    fn test_gb_salt_correction_zero_kappa() {
        let charges = vec![1.0_f64];
        let positions = vec![0.0_f64, 0.0, 0.0];
        let radii = vec![1.9e-10_f64];
        let scales = vec![0.85_f64];
        let gb = GeneralizedBorn::new(charges, positions, radii, scales, 78.5);
        assert_eq!(gb.salt_correction(0.0), 0.0);
    }
    #[test]
    fn test_solvation_shell_total_hydration() {
        let s = SolvationShell::new(6, 12, 3.5, 2.5);
        assert_eq!(s.total_hydration_number(), 18);
    }
    #[test]
    fn test_solvation_shell_residence_time_positive() {
        let s = SolvationShell::new(6, 12, 3.5, 2.5);
        assert!(s.compute_residence_time() > 0.0);
    }
    #[test]
    fn test_solvation_shell_volume_1_positive() {
        let s = SolvationShell::new(6, 12, 3.5, 2.5);
        assert!(s.shell_volume(1) > 0.0);
    }
    #[test]
    fn test_solvation_shell_volume_invalid_zero() {
        let s = SolvationShell::new(6, 12, 3.5, 2.5);
        assert_eq!(s.shell_volume(5), 0.0);
    }
    #[test]
    fn test_solvation_shell_rdf_peak() {
        let mut s = SolvationShell::new(6, 12, 3.5, 2.5);
        s.set_rdf(vec![1.0, 2.0, 3.0], vec![0.5, 2.8, 1.0]);
        assert!((s.rdf_peak() - 2.8).abs() < 1e-12);
    }
    #[test]
    fn test_solvation_shell_rdf_first_peak_position() {
        let mut s = SolvationShell::new(6, 12, 3.5, 2.5);
        s.set_rdf(vec![1.0, 2.0, 3.0], vec![0.5, 2.8, 1.0]);
        let pos = s.rdf_first_peak_position(2.0);
        assert!((pos - 2.0).abs() < 1e-12);
    }
    #[test]
    fn test_solvation_shell_coordination_number_rdf() {
        let mut s = SolvationShell::new(6, 12, 3.5, 2.5);
        let r: Vec<f64> = (0..20).map(|i| 1.0 + i as f64 * 0.1).collect();
        let gr = vec![1.0_f64; r.len()];
        s.set_rdf(r, gr);
        let rho = 0.033_f64;
        let n = s.coordination_number_from_rdf(rho, 1.0, 3.0);
        assert!(n > 0.0);
    }
    #[test]
    fn test_solvation_shell_cumulative_increasing() {
        let mut s = SolvationShell::new(6, 12, 3.5, 2.5);
        let r: Vec<f64> = (0..10).map(|i| 1.0 + i as f64 * 0.2).collect();
        let gr = vec![1.0_f64; r.len()];
        s.set_rdf(r, gr);
        let n_r = s.cumulative_coordination(0.033);
        assert!(n_r[n_r.len() - 1] >= n_r[0]);
    }
    #[test]
    fn test_solvation_shell_first_boundary_returns_finite() {
        let mut s = SolvationShell::new(6, 12, 3.5, 2.5);
        let r: Vec<f64> = (0..30).map(|i| 1.0 + i as f64 * 0.2).collect();
        let gr: Vec<f64> = r
            .iter()
            .enumerate()
            .map(|(i, _)| {
                if i < 5 {
                    i as f64 * 0.5
                } else {
                    (14 - i as i32).max(0) as f64 * 0.3
                }
            })
            .collect();
        s.set_rdf(r, gr);
        let rb = s.first_shell_boundary();
        assert!(rb.is_finite());
    }
    #[test]
    fn test_hydration_born_negative() {
        let h = HydrationFreeEnergy::new(1.0, 1.9e-10, 78.5, 0.0, 2.0, 0.0, 0.0, 2.0);
        assert!(h.born_contribution() < 0.0);
    }
    #[test]
    fn test_hydration_cavity_positive() {
        let h = HydrationFreeEnergy::new(0.0, 1.9e-10, 78.5, 0.02, 3.0, 0.0, 0.0, 3.0);
        assert!(h.cavity_contribution() > 0.0);
    }
    #[test]
    fn test_hydration_total_not_nan() {
        let h = HydrationFreeEnergy::new(1.0, 1.9e-10, 78.5, 0.02, 3.0, 1.0, 1.0, 3.0);
        assert!(!h.total_solvation_energy().is_nan());
    }
    #[test]
    fn test_hydration_dispersion_negative() {
        let h = HydrationFreeEnergy::new(0.0, 1.9e-10, 1.0, 0.0, 3.0, 1.0, 0.0, 3.0);
        assert!(h.dispersion_contribution() <= 0.0);
    }
    #[test]
    fn test_hydration_repulsion_positive() {
        let h = HydrationFreeEnergy::new(0.0, 1.9e-10, 1.0, 0.0, 3.0, 0.0, 1.0, 3.0);
        assert!(h.repulsion_contribution() >= 0.0);
    }
    #[test]
    fn test_pb_ionic_strength_nacl() {
        let mut pb = PoissonBoltzmann::new(78.5, 298.15);
        pb.add_ion(100.0, 1.0);
        pb.add_ion(100.0, -1.0);
        assert!((pb.ionic_strength() - 100.0).abs() < 1e-10);
    }
    #[test]
    fn test_pb_debye_length_positive() {
        let mut pb = PoissonBoltzmann::new(78.5, 298.15);
        pb.add_ion(100.0, 1.0);
        pb.add_ion(100.0, -1.0);
        assert!(pb.debye_length() > 0.0 && pb.debye_length().is_finite());
    }
    #[test]
    fn test_pb_debye_huckel_activity_negative() {
        let mut pb = PoissonBoltzmann::new(78.5, 298.15);
        pb.add_ion(100.0, 1.0);
        pb.add_ion(100.0, -1.0);
        let lngamma = pb.debye_huckel_activity(1.0, 1.0);
        assert!(lngamma <= 0.0);
    }
    #[test]
    fn test_pb_dh_potential_decays_with_r() {
        let mut pb = PoissonBoltzmann::new(78.5, 298.15);
        pb.add_ion(100.0, 1.0);
        pb.add_ion(100.0, -1.0);
        let q = ELEM_CHARGE;
        let phi1 = pb.dh_potential(q, 1e-10);
        let phi2 = pb.dh_potential(q, 5e-10);
        assert!(phi1 > phi2);
    }
    #[test]
    fn test_pb_electrostatic_free_energy_negative_for_ion() {
        let mut pb = PoissonBoltzmann::new(78.5, 298.15);
        pb.add_ion(100.0, 1.0);
        pb.add_ion(100.0, -1.0);
        let dg = pb.electrostatic_free_energy(1.0, 2e-10);
        assert!(dg <= 0.0);
    }
    #[test]
    fn test_pb_empty_ions_infinite_debye() {
        let pb = PoissonBoltzmann::new(78.5, 298.15);
        assert!(pb.debye_length().is_infinite());
    }
    #[test]
    fn test_widom_excess_chemical_potential_zero_insertions() {
        let s = SolvationMD::new(298.15, 1.0, 3.0);
        assert_eq!(s.excess_chemical_potential(), 0.0);
    }
    #[test]
    fn test_widom_lj_insertion_accumulates() {
        let mut s = SolvationMD::new(298.15, 1.0, 3.0);
        s.lj_insertion(3.0);
        assert_eq!(s.n_insertions(), 1);
    }
    #[test]
    fn test_widom_excess_chemical_potential_finite() {
        let mut s = SolvationMD::new(298.15, 0.5, 3.0);
        for i in 1..=10 {
            s.lj_insertion((i as f64) * 0.5 + 2.5);
        }
        let mu_ex = s.excess_chemical_potential();
        assert!(mu_ex.is_finite());
    }
    #[test]
    fn test_widom_standard_error_finite() {
        let mut s = SolvationMD::new(298.15, 0.5, 3.0);
        for i in 1..=5 {
            s.lj_insertion((i as f64) * 0.5 + 2.5);
        }
        assert!(s.standard_error().is_finite());
    }
    #[test]
    fn test_widom_n_insertions() {
        let mut s = SolvationMD::new(298.15, 1.0, 3.0);
        s.add_insertion(0.0);
        s.add_insertion(1.0);
        assert_eq!(s.n_insertions(), 2);
    }
    #[test]
    fn test_excess_mu_no_samples_infinity() {
        let ecp = ExcessChemicalPotential::new(298.15);
        assert!(ecp.estimate().is_infinite());
    }
    #[test]
    fn test_excess_mu_zero_energy_zero_mu() {
        let mut ecp = ExcessChemicalPotential::new(298.15);
        for _ in 0..100 {
            ecp.insert(0.0);
        }
        assert!(ecp.estimate().abs() < 1e-8);
    }
    #[test]
    fn test_excess_mu_negative_insertion_positive_mu() {
        let mut ecp = ExcessChemicalPotential::new(298.15);
        for _ in 0..50 {
            ecp.insert(-10.0);
        }
        assert!(ecp.estimate() < 0.0);
    }
    #[test]
    fn test_excess_mu_n_samples() {
        let mut ecp = ExcessChemicalPotential::new(298.15);
        ecp.insert(1.0);
        ecp.insert(2.0);
        assert_eq!(ecp.n_samples(), 2);
    }
    #[test]
    fn test_excess_mu_mean_insertion_energy() {
        let mut ecp = ExcessChemicalPotential::new(298.15);
        ecp.insert(4.0);
        ecp.insert(6.0);
        assert!((ecp.mean_insertion_energy() - 5.0).abs() < 1e-10);
    }
    #[test]
    fn test_excess_mu_block_average_converges() {
        let mut ecp = ExcessChemicalPotential::new(298.15);
        for _ in 0..100 {
            ecp.insert(1.0);
        }
        let (mu_ba, _std) = ecp.block_average(10);
        assert!((mu_ba - ecp.estimate()).abs() < 1.0);
    }
    #[test]
    fn test_excess_mu_standard_error_positive() {
        let mut ecp = ExcessChemicalPotential::new(298.15);
        for i in 0..10 {
            ecp.insert(i as f64 * 0.1);
        }
        assert!(ecp.standard_error() >= 0.0);
    }
    #[test]
    fn test_ti_integrate_empty_zero() {
        let ti = ThermodynamicIntegrationSolvation::new(298.15);
        assert_eq!(ti.integrate(), 0.0);
    }
    #[test]
    fn test_ti_integrate_constant_dhdl() {
        let mut ti = ThermodynamicIntegrationSolvation::new(298.15);
        ti.add_point(0.0, 10.0, 0.1);
        ti.add_point(1.0, 10.0, 0.1);
        assert!((ti.integrate() - 10.0).abs() < 1e-10);
    }
    #[test]
    fn test_ti_integrate_linear_dhdl() {
        let mut ti = ThermodynamicIntegrationSolvation::new(298.15);
        for i in 0..=10 {
            let lam = i as f64 * 0.1;
            ti.add_point(lam, lam, 0.01);
        }
        assert!((ti.integrate() - 0.5).abs() < 1e-4);
    }
    #[test]
    fn test_ti_uncertainty_positive() {
        let mut ti = ThermodynamicIntegrationSolvation::new(298.15);
        ti.add_point(0.0, 5.0, 0.1);
        ti.add_point(1.0, 5.0, 0.1);
        assert!(ti.uncertainty() > 0.0);
    }
    #[test]
    fn test_ti_simpson_matches_trapz_for_constant() {
        let mut ti = ThermodynamicIntegrationSolvation::new(298.15);
        ti.add_point(0.0, 5.0, 0.1);
        ti.add_point(0.5, 5.0, 0.1);
        ti.add_point(1.0, 5.0, 0.1);
        assert!((ti.integrate_simpson() - 5.0).abs() < 1e-10);
        assert!((ti.integrate() - 5.0).abs() < 1e-10);
    }
    #[test]
    fn test_ti_partial_integration() {
        let mut ti = ThermodynamicIntegrationSolvation::new(298.15);
        ti.add_point(0.0, 10.0, 0.1);
        ti.add_point(0.5, 10.0, 0.1);
        ti.add_point(1.0, 10.0, 0.1);
        let full = ti.integrate();
        let half = ti.partial_free_energy(0.0, 0.5);
        assert!(half < full && half > 0.0);
    }
    #[test]
    fn test_cosmo_factor_water() {
        let cs = ContinuumSolvation::new(78.5, 3.0, 0.02, 1.0);
        let f = cs.cosmo_factor();
        assert!(f > 0.9 && f < 1.0);
    }
    #[test]
    fn test_cosmo_factor_vacuum() {
        let cs = ContinuumSolvation::new(1.0, 3.0, 0.02, 0.0);
        assert!((cs.cosmo_factor()).abs() < 1e-10);
    }
    #[test]
    fn test_cosmo_electrostatic_negative_for_charged() {
        let cs = ContinuumSolvation::new(78.5, 2.0, 0.0, 1.0);
        assert!(cs.electrostatic_energy() <= 0.0);
    }
    #[test]
    fn test_cosmo_cavity_energy_positive() {
        let cs = ContinuumSolvation::new(78.5, 3.0, 0.02, 0.0);
        assert!(cs.cavity_energy() > 0.0);
    }
    #[test]
    fn test_cosmo_effective_dielectric_bulk_far() {
        let cs = ContinuumSolvation::new(78.5, 3.0, 0.0, 0.0);
        let eps_eff = cs.effective_dielectric(100.0);
        assert!(eps_eff > 70.0);
    }
    #[test]
    fn test_cosmo_effective_dielectric_inside() {
        let cs = ContinuumSolvation::new(78.5, 3.0, 0.0, 0.0);
        let eps_eff = cs.effective_dielectric(-10.0);
        assert!(eps_eff < 5.0);
    }
    #[test]
    fn test_hydrophobic_energy_positive() {
        let he = HydrophobicEffect::new(100.0, 0.0245);
        assert!((he.compute_hydrophobic_energy() - 2.45).abs() < 1e-6);
    }
    #[test]
    fn test_hydrophobic_buried_fraction() {
        let he = HydrophobicEffect::new(200.0, 0.02);
        let full = he.buried_energy(1.0);
        let half = he.buried_energy(0.5);
        assert!((half - full * 0.5).abs() < 1e-12);
    }
    #[test]
    fn test_hydrophobic_transfer_free_energy() {
        let he = HydrophobicEffect::new(100.0, 0.02).with_transfer(0.01, 50.0);
        assert!((he.transfer_free_energy() - 0.5).abs() < 1e-10);
    }
    #[test]
    fn test_hydrophobic_total_nonpolar() {
        let he = HydrophobicEffect::new(100.0, 0.02).with_transfer(0.01, 50.0);
        let expected = he.compute_hydrophobic_energy() + he.transfer_free_energy();
        assert!((he.total_nonpolar_energy() - expected).abs() < 1e-10);
    }
    #[test]
    fn test_hydrophobic_sphere_sasa() {
        let sasa = HydrophobicEffect::sphere_sasa(2.0, 1.4);
        let expected = solvent_accessible_area(2.0, 1.4);
        assert!((sasa - expected).abs() < 1e-10);
    }
    #[test]
    fn test_hydrophobic_sum_sasa() {
        let radii = vec![1.5_f64, 2.0, 2.5];
        let sum = HydrophobicEffect::sum_sasa(&radii, 1.4);
        let expected: f64 = radii.iter().map(|&r| solvent_accessible_area(r, 1.4)).sum();
        assert!((sum - expected).abs() < 1e-8);
    }
    #[test]
    fn test_solvation_fe_compute_born() {
        let sfe = SolvationFreeEnergy::new(1.0, 78.5, 1.9e-10);
        let dg = sfe.compute_born(0.0, 0.0, 0.0);
        assert!((dg - sfe.total_born_energy()).abs() < 1e-12);
    }
    #[test]
    fn test_solvation_fe_negative() {
        let sfe = SolvationFreeEnergy::new(1.0, 78.5, 1.9e-10);
        assert!(sfe.total_born_energy() < 0.0);
    }
    #[test]
    fn test_transfer_total_ddg() {
        let mut tfe = TransferFreeEnergy::new("water", "octanol");
        tfe.add_fragment("methyl", -1.0);
        tfe.add_fragment("oh", 2.0);
        assert!((tfe.total_ddg() - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_transfer_correlation_perfect() {
        let mut tfe = TransferFreeEnergy::new("water", "octanol");
        tfe.add_fragment("a", 1.0);
        tfe.add_fragment("b", 2.0);
        tfe.add_fragment("c", 3.0);
        let r = tfe.correlation_with_logp(&[1.0, 2.0, 3.0]);
        assert!((r - 1.0).abs() < 1e-6);
    }
    #[test]
    fn test_transfer_mismatched_length_zero() {
        let mut tfe = TransferFreeEnergy::new("water", "octanol");
        tfe.add_fragment("a", 1.0);
        assert_eq!(tfe.correlation_with_logp(&[1.0, 2.0]), 0.0);
    }
    #[test]
    #[allow(clippy::assertions_on_constants)]
    fn test_constants_avogadro_positive() {
        assert!(AVOGADRO > 6e23);
    }
    #[test]
    #[allow(clippy::assertions_on_constants)]
    fn test_constants_kb_positive() {
        assert!(KB > 0.0);
    }
    #[test]
    #[allow(clippy::assertions_on_constants)]
    fn test_kbt_298_reasonable() {
        assert!(KBT_298 > 2.0 && KBT_298 < 3.0);
    }
}
