//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::GalvanicSeriesEntry;

/// Faraday constant \[C/mol\]
pub const FARADAY: f64 = 96_485.0;
/// Universal gas constant \[J/(mol·K)\]
pub const GAS_CONSTANT: f64 = 8.314_462_618;
/// Nernst equation: equilibrium potential from concentration ratio.
///
/// `E = E0 − (RT / (n·F)) · ln(C_red / C_ox)`
///
/// # Arguments
/// * `e0`         — standard electrode potential \[V\]
/// * `n_electrons`— electrons transferred per reaction event
/// * `temp`       — temperature \[K\]
/// * `c_red`      — concentration of reduced species \[mol/m³\]
/// * `c_ox`       — concentration of oxidised species \[mol/m³\]
pub fn nernst_potential(e0: f64, n_electrons: u32, temp: f64, c_red: f64, c_ox: f64) -> f64 {
    let rt_over_nf = GAS_CONSTANT * temp / ((n_electrons as f64) * FARADAY);
    e0 - rt_over_nf * (c_red / c_ox).ln()
}
/// Equilibrium cell voltage from two half-cell potentials.
///
/// `E_cell = E_cathode − E_anode`
pub fn cell_voltage(e_cathode: f64, e_anode: f64) -> f64 {
    e_cathode - e_anode
}
/// Tafel analysis utilities.
///
/// Extracts electrochemical kinetic parameters from Tafel plot data.
/// A Tafel plot shows log|i| vs. overpotential η; the slope gives
/// the Tafel slope b = 2.303·RT/(α·F).
#[allow(dead_code)]
pub fn tafel_slope(alpha: f64, temp_k: f64) -> f64 {
    2.303 * GAS_CONSTANT * temp_k / (alpha * FARADAY)
}
/// Exchange current density from Tafel extrapolation.
///
/// From the linear Tafel region: `log i = log i0 + η/b`
/// At η = 0: `i0 = 10^(log_i_at_eta0)` where log_i_at_eta0 is the
/// y-intercept of the Tafel line.
///
/// # Arguments
/// * `eta_1`, `log_i_1` — point 1 on the Tafel line (V, log10(A/m²))
/// * `eta_2`, `log_i_2` — point 2 on the Tafel line (V, log10(A/m²))
#[allow(dead_code)]
pub fn exchange_current_from_tafel(eta_1: f64, log_i_1: f64, eta_2: f64, log_i_2: f64) -> f64 {
    if (eta_2 - eta_1).abs() < 1e-12 {
        return 0.0;
    }
    let slope = (log_i_2 - log_i_1) / (eta_2 - eta_1);
    let log_i0 = log_i_1 - slope * eta_1;
    10.0_f64.powf(log_i0)
}
/// Corrosion penetration rate \[mm/year\] from corrosion current density.
///
/// Uses the Faraday equation: `CPR = (i_corr * M) / (n * F * ρ) * K`
/// where K converts to mm/year (SI input units: A/m², g/mol, g/cm³).
///
/// # Arguments
/// * `i_corr`     — corrosion current density \[A/m²\]
/// * `molar_mass` — molar mass \[g/mol\]
/// * `n_electrons`— valence
/// * `density`    — metal density \[g/cm³\]
#[allow(dead_code)]
pub fn corrosion_penetration_rate(
    i_corr: f64,
    molar_mass: f64,
    n_electrons: u32,
    density: f64,
) -> f64 {
    let i_ua_cm2 = i_corr * 0.1;
    0.00327 * i_ua_cm2 * molar_mass / (n_electrons as f64 * density)
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::electrochemistry::BatteryCell;
    use crate::electrochemistry::BatteryModel;
    use crate::electrochemistry::ButlerVolmerKinetics;
    use crate::electrochemistry::ConcentrationPolarisation;
    use crate::electrochemistry::CorrosionRate;
    use crate::electrochemistry::DiffusionLayer;

    use crate::electrochemistry::ElectrodeKinetics;
    use crate::electrochemistry::ElectrolyteConductivity;
    use crate::electrochemistry::Electroplating;
    use crate::electrochemistry::EvansDiagram;
    use crate::electrochemistry::FuelCellStack;

    use crate::electrochemistry::GalvanicCouple;
    use crate::electrochemistry::ImpedanceModel;
    use crate::electrochemistry::LiIonDegradation;

    use crate::electrochemistry::SeiLayer;

    pub(super) const F_OVER_RT: f64 = 38.924;
    pub(super) const F_CONST: f64 = 96_485.0;
    #[test]
    fn test_bv_zero_eta_gives_zero_current() {
        let ek = ElectrodeKinetics::new(10.0, 0.5, 0.5, 298.15);
        let j = ek.current_density(0.0, F_OVER_RT);
        assert!(j.abs() < 1.0e-10, "BV current at η=0 should be 0, got {j}");
    }
    #[test]
    fn test_bv_positive_eta_positive_current() {
        let ek = ElectrodeKinetics::new(10.0, 0.5, 0.5, 298.15);
        let j = ek.current_density(0.1, F_OVER_RT);
        assert!(
            j > 0.0,
            "BV current at positive η should be positive, got {j}"
        );
    }
    #[test]
    fn test_limiting_current_proportional_to_diffusivity() {
        let dl1 = DiffusionLayer::new(1.0e-4, 1.0e-9, 1000.0);
        let dl2 = DiffusionLayer::new(1.0e-4, 2.0e-9, 1000.0);
        let jl1 = dl1.limiting_current(2, F_CONST, 1.0);
        let jl2 = dl2.limiting_current(2, F_CONST, 1.0);
        let ratio = jl2 / jl1;
        assert!(
            (ratio - 2.0).abs() < 1.0e-10,
            "Limiting current should scale linearly with D, got ratio {ratio}"
        );
    }
    #[test]
    fn test_battery_discharge_reduces_soc() {
        let mut cell = BatteryCell::new(10.0, 3.7, 0.05);
        assert!((cell.soc - 1.0).abs() < 1.0e-12);
        cell.discharge(10.0, 3600.0);
        assert!(
            cell.soc < 1.0,
            "SOC should decrease after discharge, got {}",
            cell.soc
        );
        assert!(
            cell.soc.abs() < 1.0e-10,
            "SOC should be ~0 after full discharge"
        );
    }
    #[test]
    fn test_terminal_voltage_at_zero_current_equals_ocv() {
        let cell = BatteryCell::new(10.0, 3.7, 0.05);
        let ocv = cell.open_circuit_voltage();
        let vt = cell.terminal_voltage(0.0);
        assert!(
            (vt - ocv).abs() < 1.0e-12,
            "Terminal voltage at I=0 should equal OCV: {ocv}, got {vt}"
        );
    }
    #[test]
    fn test_bvk_zero_eta() {
        let bv = ButlerVolmerKinetics::new(5.0, 0.5, 0.5, 298.15);
        let j = bv.current_density(0.0);
        assert!(j.abs() < 1.0e-10, "got {j}");
    }
    #[test]
    fn test_bvk_positive_eta_anodic() {
        let bv = ButlerVolmerKinetics::new(5.0, 0.5, 0.5, 298.15);
        let j = bv.current_density(0.1);
        assert!(
            j > 0.0,
            "anodic overpotential should give positive j, got {j}"
        );
    }
    #[test]
    fn test_bvk_negative_eta_cathodic() {
        let bv = ButlerVolmerKinetics::new(5.0, 0.5, 0.5, 298.15);
        let j = bv.current_density(-0.1);
        assert!(
            j < 0.0,
            "cathodic overpotential should give negative j, got {j}"
        );
    }
    #[test]
    fn test_bvk_f_over_rt_at_298() {
        let bv = ButlerVolmerKinetics::new(1.0, 0.5, 0.5, 298.15);
        let q = bv.f_over_rt();
        assert!((q - 38.924).abs() < 0.01, "got {q}");
    }
    #[test]
    fn test_bvk_charge_transfer_resistance_symmetric() {
        let bv = ButlerVolmerKinetics::new(10.0, 0.5, 0.5, 298.15);
        let rct = bv.charge_transfer_resistance();
        let expected = GAS_CONSTANT * 298.15 / (10.0 * FARADAY);
        assert!(
            (rct - expected).abs() < 1.0e-12,
            "got {rct}, expected {expected}"
        );
    }
    #[test]
    fn test_bvk_tafel_slope_anodic() {
        let bv = ButlerVolmerKinetics::new(1.0, 0.5, 0.5, 298.15);
        let b = bv.tafel_slope_anodic();
        assert!(b > 0.10 && b < 0.14, "got {b}");
    }
    #[test]
    fn test_nernst_equal_concentrations() {
        let e = nernst_potential(0.34, 2, 298.15, 1.0, 1.0);
        assert!((e - 0.34).abs() < 1.0e-12, "got {e}");
    }
    #[test]
    fn test_nernst_shifts_with_concentration() {
        let e_eq = nernst_potential(0.34, 2, 298.15, 1.0, 1.0);
        let e_ox = nernst_potential(0.34, 2, 298.15, 0.01, 1.0);
        assert!(e_ox > e_eq, "e_ox={e_ox} should be > e_eq={e_eq}");
    }
    #[test]
    fn test_cell_voltage_difference() {
        let v = cell_voltage(1.23, 0.0);
        assert!((v - 1.23).abs() < 1.0e-12, "got {v}");
    }
    #[test]
    fn test_battery_model_initial_state() {
        let bat = BatteryModel::new(50.0, 0.01, 0.05, 2000.0);
        assert!((bat.soc - 1.0).abs() < 1.0e-12);
        assert_eq!(bat.u_rc, 0.0);
    }
    #[test]
    fn test_battery_model_terminal_voltage_no_load() {
        let bat = BatteryModel::new(50.0, 0.01, 0.05, 2000.0);
        let ocv = bat.open_circuit_voltage();
        let vt = bat.terminal_voltage(0.0);
        assert!((vt - ocv).abs() < 1.0e-12, "got {vt}");
    }
    #[test]
    fn test_battery_model_step_discharges_soc() {
        let mut bat = BatteryModel::new(1.0, 0.01, 0.05, 1000.0);
        let soc_init = bat.soc;
        bat.step(1.0, 100.0);
        assert!(bat.soc < soc_init, "SOC should decrease");
    }
    #[test]
    fn test_battery_model_power_positive_on_discharge() {
        let bat = BatteryModel::new(10.0, 0.01, 0.05, 500.0);
        let p = bat.power(5.0);
        assert!(p > 0.0, "power should be positive on discharge, got {p}");
    }
    #[test]
    fn test_battery_model_depletion_flag() {
        let mut bat = BatteryModel::new(1.0, 0.01, 0.05, 1000.0);
        bat.step(1.0, 3600.0);
        assert!(
            bat.is_depleted(),
            "battery should be depleted after full discharge"
        );
    }
    #[test]
    fn test_corrosion_rate_iron() {
        let cr = CorrosionRate::new(1.0e-3, 55.85, 2, 7.87);
        let rate = cr.mm_per_year();
        assert!(rate > 0.0, "corrosion rate should be positive, got {rate}");
    }
    #[test]
    fn test_corrosion_rate_higher_current_faster() {
        let cr1 = CorrosionRate::new(1.0e-3, 55.85, 2, 7.87);
        let cr2 = CorrosionRate::new(2.0e-3, 55.85, 2, 7.87);
        assert!(cr2.mm_per_year() > cr1.mm_per_year());
    }
    #[test]
    fn test_stern_geary_positive() {
        let i_corr = CorrosionRate::from_polarisation_resistance(1000.0, 0.12, 0.12);
        assert!(i_corr > 0.0, "i_corr should be positive, got {i_corr}");
    }
    #[test]
    fn test_fuel_cell_voltage_at_zero_current() {
        let fc = FuelCellStack::new(1.23, 1.0, 0.07, 1.5e-5, 1.5e4, 3e-5, 8e-5);
        let v = fc.cell_voltage(0.0);
        assert!((v - 1.23).abs() < 0.01, "got {v}");
    }
    #[test]
    fn test_fuel_cell_voltage_decreases_with_current() {
        let fc = FuelCellStack::new(1.23, 1.0, 0.07, 1.5e-5, 1.5e4, 3e-5, 8e-5);
        let v1 = fc.cell_voltage(100.0);
        let v2 = fc.cell_voltage(1000.0);
        assert!(
            v2 < v1,
            "voltage should decrease with current density: v1={v1}, v2={v2}"
        );
    }
    #[test]
    fn test_fuel_cell_voltage_zero_at_limiting() {
        let fc = FuelCellStack::new(1.23, 1.0, 0.07, 1.5e-5, 1.5e4, 3e-5, 8e-5);
        let v = fc.cell_voltage(1.5e4);
        assert_eq!(v, 0.0, "voltage at limiting current should be 0");
    }
    #[test]
    fn test_fuel_cell_power_density_positive() {
        let fc = FuelCellStack::new(1.23, 1.0, 0.07, 1.5e-5, 1.5e4, 3e-5, 8e-5);
        let p = fc.power_density(500.0);
        assert!(p > 0.0, "got {p}");
    }
    #[test]
    fn test_fuel_cell_efficiency_less_than_unity() {
        let fc = FuelCellStack::new(1.23, 1.0, 0.07, 1.5e-5, 1.5e4, 3e-5, 8e-5);
        let eff = fc.efficiency(500.0);
        assert!(
            eff > 0.0 && eff <= 1.0,
            "efficiency should be in (0,1], got {eff}"
        );
    }
    #[test]
    fn test_impedance_z_real_dc_limit() {
        let imp = ImpedanceModel::new(1.0, 10.0, 1e-3, 0.01);
        let z_r = imp.z_real(0.0);
        assert!(
            (z_r - 11.0).abs() < 1e-6,
            "DC Z_real = R_s + R_ct, got {z_r}"
        );
    }
    #[test]
    fn test_impedance_z_real_high_frequency() {
        let imp = ImpedanceModel::new(1.0, 10.0, 1e-3, 0.01);
        let z_r = imp.z_real(1e8);
        assert!(
            z_r < 11.0,
            "high-ω Z_real should approach R_s = 1.0, got {z_r}"
        );
    }
    #[test]
    fn test_impedance_magnitude_positive() {
        let imp = ImpedanceModel::new(1.0, 10.0, 1e-3, 0.01);
        let z = imp.z_magnitude(100.0);
        assert!(z > 0.0, "impedance magnitude should be positive, got {z}");
    }
    #[test]
    fn test_impedance_nyquist_point_real_part_positive() {
        let imp = ImpedanceModel::new(1.0, 10.0, 1e-3, 0.01);
        let (zr, _zi) = imp.nyquist_point(100.0);
        assert!(zr > 0.0, "Z' should be positive, got {zr}");
    }
    #[test]
    fn test_impedance_angular_frequency() {
        let omega = ImpedanceModel::angular_frequency(1.0);
        assert!((omega - 2.0 * std::f64::consts::PI).abs() < 1e-10);
    }
    #[test]
    fn test_liion_initial_soh_unity() {
        let bat = LiIonDegradation::new(50.0, 0.05);
        assert!((bat.soh_capacity() - 1.0).abs() < 1e-12);
        assert!((bat.soh_resistance() - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_liion_capacity_fade() {
        let mut bat = LiIonDegradation::new(50.0, 0.05);
        bat.update_cycle_fade(100.0, 0.01);
        assert!(
            bat.soh_capacity() < 1.0,
            "capacity should decrease after cycling"
        );
    }
    #[test]
    fn test_liion_calendar_aging_increases_resistance() {
        let mut bat = LiIonDegradation::new(50.0, 0.05);
        bat.update_calendar_aging(365.0, 0.01);
        assert!(
            bat.resistance > bat.resistance_init,
            "resistance should grow with age"
        );
    }
    #[test]
    fn test_liion_not_end_of_life_initially() {
        let bat = LiIonDegradation::new(50.0, 0.05);
        assert!(!bat.is_end_of_life(), "fresh battery should not be at EOL");
    }
    #[test]
    fn test_liion_end_of_life_after_heavy_cycling() {
        let mut bat = LiIonDegradation::new(50.0, 0.05);
        bat.update_cycle_fade(10000.0, 0.025);
        assert!(
            bat.is_end_of_life(),
            "heavily cycled battery should reach EOL"
        );
    }
    #[test]
    fn test_liion_remaining_cycles_positive() {
        let bat = LiIonDegradation::new(50.0, 0.05);
        let rem = bat.estimated_remaining_cycles(0.001);
        assert!(rem.is_none() || rem.unwrap() > 0);
    }
    #[test]
    fn test_electrolyte_molar_conductivity_dilute() {
        let ec = ElectrolyteConductivity::new(0.01499, 9.39e-3, 298.15);
        let lambda = ec.molar_conductivity(0.001);
        assert!(
            lambda > 0.0 && lambda <= ec.lambda_0,
            "Λ should be in (0, Λ₀]"
        );
    }
    #[test]
    fn test_electrolyte_specific_conductivity_increases_with_conc() {
        let ec = ElectrolyteConductivity::new(0.01499, 0.001, 298.15);
        let k1 = ec.specific_conductivity(0.01);
        let k2 = ec.specific_conductivity(0.1);
        assert!(
            k2 > k1,
            "conductivity should increase with concentration in dilute range"
        );
    }
    #[test]
    fn test_electrolyte_resistivity_finite() {
        let ec = ElectrolyteConductivity::new(0.01499, 0.001, 298.15);
        let rho = ec.resistivity(0.01);
        assert!(
            rho > 0.0 && rho.is_finite(),
            "resistivity should be finite, got {rho}"
        );
    }
    #[test]
    fn test_electrolyte_temperature_correction_positive() {
        let ec = ElectrolyteConductivity::new(0.01499, 0.001, 350.0);
        let kappa_corr = ec.temperature_corrected_conductivity(0.1, 298.15, 15000.0);
        assert!(
            kappa_corr > 0.0,
            "temperature-corrected conductivity should be positive"
        );
    }
    #[test]
    fn test_electrolyte_transference_number_sum_unity() {
        let t_plus = ElectrolyteConductivity::transference_number_cation(76.3, 76.3);
        assert!(
            (t_plus - 0.5).abs() < 1e-10,
            "equal λ → t+ = 0.5, got {t_plus}"
        );
    }
    #[test]
    fn test_electroplating_mass_deposited_positive() {
        let ep = Electroplating::new(58.69, 2, 0.95);
        let m = ep.mass_deposited_g(1.0, 3600.0);
        assert!(m > 0.0, "mass deposited should be positive, got {m}");
    }
    #[test]
    fn test_electroplating_mass_proportional_to_time() {
        let ep = Electroplating::new(58.69, 2, 1.0);
        let m1 = ep.mass_deposited_g(1.0, 1000.0);
        let m2 = ep.mass_deposited_g(1.0, 2000.0);
        assert!(
            (m2 / m1 - 2.0).abs() < 1e-10,
            "mass ∝ time, got ratio {}",
            m2 / m1
        );
    }
    #[test]
    fn test_electroplating_mass_proportional_to_current() {
        let ep = Electroplating::new(58.69, 2, 1.0);
        let m1 = ep.mass_deposited_g(1.0, 1000.0);
        let m2 = ep.mass_deposited_g(3.0, 1000.0);
        assert!(
            (m2 / m1 - 3.0).abs() < 1e-10,
            "mass ∝ current, got ratio {}",
            m2 / m1
        );
    }
    #[test]
    fn test_electroplating_faraday_known_value() {
        let ep = Electroplating::new(63.55, 2, 1.0);
        let m = ep.mass_deposited_g(1.0, 1.0);
        assert!((m - 63.55 / (2.0 * 96485.0)).abs() < 1e-8, "got {m}");
    }
    #[test]
    fn test_electroplating_thickness_positive() {
        let ep = Electroplating::new(58.69, 2, 0.95);
        let th = ep.thickness_um(1.0, 3600.0, 0.01, 8.9);
        assert!(th > 0.0, "thickness should be positive, got {th}");
    }
    #[test]
    fn test_electroplating_time_for_mass_positive() {
        let ep = Electroplating::new(63.55, 2, 1.0);
        let t = ep.time_for_mass(1.0, 1.0);
        assert!(t > 0.0, "time should be positive, got {t}");
    }
    #[test]
    fn test_electroplating_plating_rate_positive() {
        let ep = Electroplating::new(58.69, 2, 0.95);
        let rate = ep.plating_rate_um_per_min(500.0, 8.9);
        assert!(rate > 0.0, "plating rate should be positive, got {rate}");
    }
    #[test]
    fn test_evans_diagram_anodic_current_increases_with_potential() {
        let evans = EvansDiagram::new(1e-3, -0.44, 0.5, 1e-5, 0.0, 0.5, 298.15);
        let i1 = evans.anodic_current(-0.44);
        let i2 = evans.anodic_current(-0.30);
        assert!(
            i2 > i1,
            "anodic current should increase with potential: {i1} vs {i2}"
        );
    }
    #[test]
    fn test_evans_diagram_cathodic_current_decreases_with_potential() {
        let evans = EvansDiagram::new(1e-3, -0.44, 0.5, 1e-5, 0.0, 0.5, 298.15);
        let i1 = evans.cathodic_current(-0.1);
        let i2 = evans.cathodic_current(-0.5);
        assert!(
            i2 > i1,
            "cathodic current should increase as potential decreases: {i1} vs {i2}"
        );
    }
    #[test]
    fn test_evans_diagram_corrosion_potential_between_standard_potentials() {
        let evans = EvansDiagram::new(1e-3, -0.44, 0.5, 1e-5, 0.0, 0.5, 298.15);
        let e_corr = evans.corrosion_potential(-0.6, 0.1, 1e-6);
        assert!(e_corr.is_some(), "corrosion potential should be found");
        let e = e_corr.unwrap();
        assert!(
            e > -0.6 && e < 0.1,
            "E_corr should be between electrode potentials, got {e}"
        );
    }
    #[test]
    fn test_evans_diagram_corrosion_current_positive() {
        let evans = EvansDiagram::new(1e-3, -0.44, 0.5, 1e-5, 0.0, 0.5, 298.15);
        let i_corr = evans.corrosion_current(-0.6, 0.1, 1e-6);
        assert!(i_corr.is_some(), "corrosion current should be computed");
        assert!(
            i_corr.unwrap() > 0.0,
            "corrosion current should be positive"
        );
    }
    #[test]
    fn test_galvanic_couple_anodic_current_increases() {
        let couple = GalvanicCouple::new(-1.1, 0.08, 1e-4, -0.44, 0.12, 1e-5);
        let i1 = couple.anodic_current_1(-1.1);
        let i2 = couple.anodic_current_1(-0.9);
        assert!(i2 > i1, "anodic Zn current should increase as E increases");
    }
    #[test]
    fn test_galvanic_couple_potential_between_corrosion_potentials() {
        let couple = GalvanicCouple::new(-1.1, 0.08, 1e-4, -0.44, 0.12, 1e-5);
        let e_c = couple.couple_potential(1e-6);
        if let Some(e) = e_c {
            assert!(
                e > -1.3 && e < -0.3,
                "couple potential should be between Zn and Fe corrosion potentials, got {e}"
            );
        }
    }
    #[test]
    fn test_concentration_polarisation_surface_conc_at_zero_current() {
        let cp = ConcentrationPolarisation::new(1000.0, 1000.0, 2, 298.15);
        let c_s = cp.surface_concentration(0.0);
        assert!(
            (c_s - 1000.0).abs() < 1e-10,
            "at zero current, C_s = C_bulk, got {c_s}"
        );
    }
    #[test]
    fn test_concentration_polarisation_surface_conc_decreases_with_current() {
        let cp = ConcentrationPolarisation::new(1000.0, 500.0, 2, 298.15);
        let c_s1 = cp.surface_concentration(100.0);
        let c_s2 = cp.surface_concentration(400.0);
        assert!(c_s2 < c_s1, "C_s decreases with current: {c_s1} vs {c_s2}");
    }
    #[test]
    fn test_concentration_polarisation_overpotential_zero_at_zero_current() {
        let cp = ConcentrationPolarisation::new(1000.0, 1000.0, 2, 298.15);
        let eta = cp.overpotential(0.0);
        assert!(
            eta.abs() < 1e-10,
            "no overpotential at zero current, got {eta}"
        );
    }
    #[test]
    fn test_concentration_polarisation_overpotential_increases_with_current() {
        let cp = ConcentrationPolarisation::new(1000.0, 500.0, 2, 298.15);
        let eta1 = cp.overpotential(100.0);
        let eta2 = cp.overpotential(400.0);
        assert!(
            eta2 > eta1,
            "overpotential increases with current: {eta1} vs {eta2}"
        );
    }
    #[test]
    fn test_concentration_polarisation_current_fraction() {
        let cp = ConcentrationPolarisation::new(1000.0, 500.0, 2, 298.15);
        let frac = cp.current_fraction(250.0);
        assert!((frac - 0.5).abs() < 1e-10, "fraction = j/j_L, got {frac}");
    }
    #[test]
    fn test_sei_layer_initial_resistance() {
        let sei = SeiLayer::new(0.005, 1e-5, 298.15, 30000.0);
        let r = sei.resistance(0.0, 298.15);
        assert!((r - 0.005).abs() < 1e-10, "initial R_SEI = R0, got {r}");
    }
    #[test]
    fn test_sei_layer_resistance_grows_with_time() {
        let sei = SeiLayer::new(0.005, 1e-5, 298.15, 30000.0);
        let r1 = sei.resistance(1e6, 298.15);
        let r2 = sei.resistance(4e6, 298.15);
        assert!(
            r2 > r1,
            "SEI resistance should grow with time: {r1} vs {r2}"
        );
    }
    #[test]
    fn test_sei_layer_resistance_greater_at_higher_temp() {
        let sei = SeiLayer::new(0.005, 1e-5, 298.15, 30000.0);
        let r_low = sei.resistance(1e6, 298.15);
        let r_high = sei.resistance(1e6, 333.15);
        assert!(
            r_high > r_low,
            "higher T → faster SEI growth: {r_low} vs {r_high}"
        );
    }
    #[test]
    fn test_sei_layer_time_to_double_positive() {
        let sei = SeiLayer::new(0.01, 1e-4, 298.15, 30000.0);
        let t = sei.time_to_double_resistance(298.15);
        assert!(t > 0.0, "time to double should be positive, got {t}");
    }
    #[test]
    fn test_tafel_slope_standard_value() {
        let b = tafel_slope(0.5, 298.15);
        assert!((b - 0.1183).abs() < 0.002, "Tafel slope = {b}");
    }
    #[test]
    fn test_tafel_slope_increases_with_temperature() {
        let b1 = tafel_slope(0.5, 298.15);
        let b2 = tafel_slope(0.5, 350.0);
        assert!(b2 > b1, "Tafel slope increases with T: {b1} vs {b2}");
    }
    #[test]
    fn test_exchange_current_from_tafel_extrapolation() {
        let log_i0 = (1.0e-3_f64).log10();
        let log_i1 = log_i0 + 0.1 / 0.12;
        let i0 = exchange_current_from_tafel(0.0, log_i0, 0.1, log_i1);
        assert!((i0 - 1.0e-3).abs() / 1.0e-3 < 0.01, "i0 = {i0}");
    }
    #[test]
    fn test_corrosion_penetration_rate_positive() {
        let cpr = corrosion_penetration_rate(1.0e-3, 55.85, 2, 7.87);
        assert!(cpr > 0.0, "CPR should be positive, got {cpr}");
    }
    #[test]
    fn test_corrosion_penetration_rate_proportional_to_current() {
        let cpr1 = corrosion_penetration_rate(1.0e-3, 55.85, 2, 7.87);
        let cpr2 = corrosion_penetration_rate(2.0e-3, 55.85, 2, 7.87);
        assert!(
            (cpr2 / cpr1 - 2.0).abs() < 1e-10,
            "CPR ∝ i_corr, got ratio {}",
            cpr2 / cpr1
        );
    }
}
/// Symmetric Butler-Volmer: current density with equal transfer coefficients.
///
/// `j = 2 * j0 * sinh(α·F·η / RT)`  for αa = αc = α.
///
/// # Arguments
/// * `j0`    — exchange current density \[A/m²\]
/// * `alpha` — transfer coefficient (0 < α < 1, typically 0.5)
/// * `eta`   — overpotential \[V\]
/// * `temp`  — temperature \[K\]
#[allow(dead_code)]
pub fn symmetric_bv_current(j0: f64, alpha: f64, eta: f64, temp: f64) -> f64 {
    let q = FARADAY / (GAS_CONSTANT * temp);
    2.0 * j0 * (alpha * q * eta).sinh()
}
/// Tafel (high-overpotential) approximation to Butler-Volmer.
///
/// For η >> RT/(αF): `j ≈ j0 * exp(α·F·η / RT)` (anodic branch)
///
/// # Arguments
/// * `j0`    — exchange current density \[A/m²\]
/// * `alpha` — transfer coefficient
/// * `eta`   — overpotential \[V\] (positive = anodic)
/// * `temp`  — temperature \[K\]
#[allow(dead_code)]
pub fn tafel_current_anodic(j0: f64, alpha: f64, eta: f64, temp: f64) -> f64 {
    let q = FARADAY / (GAS_CONSTANT * temp);
    j0 * (alpha * q * eta).exp()
}
/// Linearised Butler-Volmer for small overpotentials (η << RT/F).
///
/// `j ≈ j0 * (αa + αc) * F * η / RT`
///
/// # Arguments
/// * `j0`      — exchange current density \[A/m²\]
/// * `alpha_a`, `alpha_c` — transfer coefficients
/// * `eta`     — overpotential \[V\]
/// * `temp`    — temperature \[K\]
#[allow(dead_code)]
pub fn linearised_bv_current(j0: f64, alpha_a: f64, alpha_c: f64, eta: f64, temp: f64) -> f64 {
    let q = FARADAY / (GAS_CONSTANT * temp);
    j0 * (alpha_a + alpha_c) * q * eta
}
/// Nernst equation with activities (dimensionless).
///
/// `E = E0 - (RT/(n·F)) * ln(a_red / a_ox)`
///
/// where `a_red` and `a_ox` are the activities of reduced and oxidised species.
///
/// # Arguments
/// * `e0`      — standard potential \[V\]
/// * `n`       — electrons transferred
/// * `temp`    — temperature \[K\]
/// * `a_red`   — activity of reduced form (dimensionless)
/// * `a_ox`    — activity of oxidised form
#[allow(dead_code)]
pub fn nernst_with_activity(e0: f64, n: u32, temp: f64, a_red: f64, a_ox: f64) -> f64 {
    let rt_over_nf = GAS_CONSTANT * temp / (n as f64 * FARADAY);
    e0 - rt_over_nf * (a_red / a_ox).ln()
}
/// Nernst potential shift \[V\] due to a tenfold change in concentration ratio.
///
/// `ΔE = -(RT / n·F) * ln(10)`  per decade change.
#[allow(dead_code)]
pub fn nernst_shift_per_decade(n: u32, temp: f64) -> f64 {
    let rt_over_nf = GAS_CONSTANT * temp / (n as f64 * FARADAY);
    rt_over_nf * 10.0_f64.ln()
}
/// Levich equation for a rotating disk electrode (RDE).
///
/// `j_L = 0.620 * n * F * D^(2/3) * ω^(1/2) * ν^(-1/6) * C_bulk`
///
/// where:
/// - `n` — electrons transferred
/// - `D` — diffusion coefficient \[m²/s\]
/// - `ω` — angular velocity \[rad/s\]
/// - `ν` — kinematic viscosity \[m²/s\]
/// - `C_bulk` — bulk concentration \[mol/m³\]
///
/// Returns the limiting current density \[A/m²\].
#[allow(dead_code)]
pub fn levich_limiting_current(
    n_electrons: u32,
    diffusivity: f64,
    angular_velocity: f64,
    kinematic_viscosity: f64,
    bulk_concentration: f64,
) -> f64 {
    0.620
        * (n_electrons as f64)
        * FARADAY
        * diffusivity.powf(2.0 / 3.0)
        * angular_velocity.sqrt()
        * kinematic_viscosity.powf(-1.0 / 6.0)
        * bulk_concentration
}
/// Koutecky-Levich analysis: combined kinetic and diffusion control.
///
/// `1/j = 1/j_k + 1/j_L`
///
/// where `j_k` is the kinetic current density and `j_L` is the Levich limiting current.
/// Returns the total measured current density \[A/m²\].
#[allow(dead_code)]
pub fn koutecky_levich_current(j_kinetic: f64, j_levich: f64) -> f64 {
    if j_kinetic.abs() < f64::EPSILON || j_levich.abs() < f64::EPSILON {
        return 0.0;
    }
    1.0 / (1.0 / j_kinetic + 1.0 / j_levich)
}
/// Standard galvanic series (selected metals, vs SHE in seawater).
#[allow(dead_code)]
pub const GALVANIC_SERIES: &[GalvanicSeriesEntry] = &[
    GalvanicSeriesEntry {
        name: "Magnesium",
        potential_v: -1.63,
    },
    GalvanicSeriesEntry {
        name: "Zinc",
        potential_v: -1.10,
    },
    GalvanicSeriesEntry {
        name: "Aluminium 2024",
        potential_v: -0.80,
    },
    GalvanicSeriesEntry {
        name: "Carbon steel",
        potential_v: -0.61,
    },
    GalvanicSeriesEntry {
        name: "Cast iron",
        potential_v: -0.58,
    },
    GalvanicSeriesEntry {
        name: "Stainless 304 (act)",
        potential_v: -0.50,
    },
    GalvanicSeriesEntry {
        name: "Lead",
        potential_v: -0.40,
    },
    GalvanicSeriesEntry {
        name: "Tin",
        potential_v: -0.32,
    },
    GalvanicSeriesEntry {
        name: "Nickel (passive)",
        potential_v: -0.10,
    },
    GalvanicSeriesEntry {
        name: "Copper",
        potential_v: 0.05,
    },
    GalvanicSeriesEntry {
        name: "Silver",
        potential_v: 0.40,
    },
    GalvanicSeriesEntry {
        name: "Gold",
        potential_v: 0.90,
    },
];
/// Galvanic current enhancement factor due to cathode-to-anode area ratio.
///
/// For a galvanic couple, a large cathodic area relative to the anode
/// accelerates anodic corrosion:
///
/// `i_anode = i_couple * (A_cathode / A_anode)`
///
/// Returns the accelerated corrosion current \[A/m²\] on the anode.
#[allow(dead_code)]
pub fn galvanic_area_ratio_current(
    i_couple_a: f64,
    area_cathode_m2: f64,
    area_anode_m2: f64,
) -> f64 {
    if area_anode_m2 < f64::EPSILON {
        return f64::INFINITY;
    }
    i_couple_a / area_anode_m2 * area_cathode_m2 / area_anode_m2
}
#[cfg(test)]
mod tests_new {

    use crate::electrochemistry::DoubleLayerCapacitance;

    use crate::electrochemistry::GALVANIC_SERIES;

    use crate::electrochemistry::PeukertModel;

    use crate::electrochemistry::galvanic_area_ratio_current;
    use crate::electrochemistry::koutecky_levich_current;
    use crate::electrochemistry::levich_limiting_current;
    use crate::electrochemistry::linearised_bv_current;
    use crate::electrochemistry::nernst_shift_per_decade;
    use crate::electrochemistry::nernst_with_activity;
    use crate::electrochemistry::symmetric_bv_current;
    use crate::electrochemistry::tafel_current_anodic;
    #[test]
    fn test_peukert_available_capacity_at_nominal_rate() {
        let p = PeukertModel::lead_acid_100ah();
        let q = p.available_capacity_ah(p.i_nominal_a);
        assert!(
            (q - p.capacity_nominal_ah).abs() < 1e-6,
            "At I_nom, Q = Q_nom, got {q}"
        );
    }
    #[test]
    fn test_peukert_capacity_decreases_at_higher_current() {
        let p = PeukertModel::lead_acid_100ah();
        let q_nom = p.available_capacity_ah(p.i_nominal_a);
        let q_high = p.available_capacity_ah(p.i_nominal_a * 5.0);
        assert!(
            q_high < q_nom,
            "Higher discharge rate → lower available capacity"
        );
    }
    #[test]
    fn test_peukert_discharge_time_positive() {
        let p = PeukertModel::lead_acid_100ah();
        let t = p.discharge_time_h(10.0);
        assert!(t > 0.0, "Discharge time should be positive, got {t}");
    }
    #[test]
    fn test_peukert_discharge_time_decreases_with_current() {
        let p = PeukertModel::lead_acid_100ah();
        let t1 = p.discharge_time_h(5.0);
        let t2 = p.discharge_time_h(20.0);
        assert!(
            t2 < t1,
            "Higher current → shorter discharge time: {t1} vs {t2}"
        );
    }
    #[test]
    fn test_peukert_c_rate_at_nominal() {
        let p = PeukertModel::new(100.0, 5.0, 1.2);
        let cr = p.c_rate(5.0);
        assert!((cr - 0.05).abs() < 1e-10, "C-rate = 5/100 = 0.05, got {cr}");
    }
    #[test]
    fn test_peukert_capacity_fade_factor_unity_at_nominal() {
        let p = PeukertModel::lead_acid_100ah();
        let f = p.capacity_fade_factor(p.i_nominal_a);
        assert!((f - 1.0).abs() < 1e-6, "No fade at nominal rate, got {f}");
    }
    #[test]
    fn test_peukert_li_ion_less_fade_than_lead_acid() {
        let la = PeukertModel::lead_acid_100ah();
        let li = PeukertModel::li_ion_10ah();
        let current_la = la.capacity_nominal_ah * 10.0;
        let current_li = li.capacity_nominal_ah * 10.0;
        let fade_la = la.capacity_fade_factor(current_la);
        let fade_li = li.capacity_fade_factor(current_li);
        assert!(
            fade_li > fade_la,
            "Li-ion retains more capacity at high C-rate"
        );
    }
    #[test]
    fn test_symmetric_bv_zero_at_zero_eta() {
        let j = symmetric_bv_current(10.0, 0.5, 0.0, 298.15);
        assert!(j.abs() < 1e-10, "symmetric BV: j=0 at η=0, got {j}");
    }
    #[test]
    fn test_symmetric_bv_odd_function_of_eta() {
        let j_pos = symmetric_bv_current(5.0, 0.5, 0.1, 298.15);
        let j_neg = symmetric_bv_current(5.0, 0.5, -0.1, 298.15);
        assert!(
            (j_pos + j_neg).abs() < 1e-10,
            "symmetric BV is odd: j(η) = -j(-η)"
        );
    }
    #[test]
    fn test_tafel_current_increases_with_eta() {
        let j1 = tafel_current_anodic(1.0, 0.5, 0.1, 298.15);
        let j2 = tafel_current_anodic(1.0, 0.5, 0.3, 298.15);
        assert!(j2 > j1, "Tafel current increases with overpotential");
    }
    #[test]
    fn test_linearised_bv_proportional_to_eta() {
        let j1 = linearised_bv_current(1.0, 0.5, 0.5, 0.001, 298.15);
        let j2 = linearised_bv_current(1.0, 0.5, 0.5, 0.002, 298.15);
        assert!(
            (j2 / j1 - 2.0).abs() < 1e-6,
            "linearised BV ∝ η, ratio = {}",
            j2 / j1
        );
    }
    #[test]
    fn test_nernst_with_activity_equal_activities() {
        let e = nernst_with_activity(0.34, 2, 298.15, 1.0, 1.0);
        assert!(
            (e - 0.34).abs() < 1e-12,
            "Equal activities → E = E0, got {e}"
        );
    }
    #[test]
    fn test_nernst_shift_per_decade_positive() {
        let shift = nernst_shift_per_decade(1, 298.15);
        assert!(
            shift > 0.05 && shift < 0.07,
            "Nernst shift ~59 mV/decade, got {shift}"
        );
    }
    #[test]
    fn test_nernst_shift_per_decade_halved_for_two_electrons() {
        let shift_1 = nernst_shift_per_decade(1, 298.15);
        let shift_2 = nernst_shift_per_decade(2, 298.15);
        assert!((shift_1 / shift_2 - 2.0).abs() < 1e-10, "Shift ∝ 1/n");
    }
    #[test]
    fn test_levich_limiting_current_positive() {
        let j_l = levich_limiting_current(1, 1e-9, 100.0, 1e-6, 1.0);
        assert!(j_l > 0.0, "Levich j_L should be positive, got {j_l}");
    }
    #[test]
    fn test_levich_increases_with_rotation_speed() {
        let j_l1 = levich_limiting_current(1, 1e-9, 100.0, 1e-6, 1.0);
        let j_l2 = levich_limiting_current(1, 1e-9, 400.0, 1e-6, 1.0);
        assert!(
            (j_l2 / j_l1 - 2.0).abs() < 0.01,
            "j_L ∝ √ω: ratio = {}",
            j_l2 / j_l1
        );
    }
    #[test]
    fn test_koutecky_levich_less_than_either_limit() {
        let j_k = 100.0;
        let j_l = 200.0;
        let j = koutecky_levich_current(j_k, j_l);
        assert!(j < j_k && j < j_l, "Combined current < min(j_k, j_L): {j}");
    }
    #[test]
    fn test_koutecky_levich_approaches_kinetic_limit_at_large_jl() {
        let j_k = 10.0;
        let j_l = 1e10;
        let j = koutecky_levich_current(j_k, j_l);
        assert!((j - j_k).abs() / j_k < 1e-5, "KL → j_k for large j_L: {j}");
    }
    #[test]
    fn test_galvanic_series_has_entries() {
        assert!(
            !GALVANIC_SERIES.is_empty(),
            "Galvanic series should have entries"
        );
    }
    #[test]
    fn test_galvanic_series_copper_more_noble_than_zinc() {
        let cu = GALVANIC_SERIES.iter().find(|e| e.name == "Copper").unwrap();
        let zn = GALVANIC_SERIES.iter().find(|e| e.name == "Zinc").unwrap();
        assert!(
            cu.potential_v > zn.potential_v,
            "Copper more noble than Zinc"
        );
    }
    #[test]
    fn test_galvanic_area_ratio_current_larger_cathode() {
        let i_couple = 1.0e-3;
        let a_c = 1.0;
        let a_a = 0.01;
        let j_anode = galvanic_area_ratio_current(i_couple, a_c, a_a);
        assert!(
            j_anode > i_couple / a_a,
            "Large cathode area accelerates anode corrosion"
        );
    }
    #[test]
    fn test_debye_length_positive() {
        let dl = DoubleLayerCapacitance::kcl_01_mol_l();
        let kappa_inv = dl.debye_length();
        assert!(
            kappa_inv > 0.0 && kappa_inv < 0.01,
            "Debye length should be nm-scale, got {kappa_inv} m"
        );
    }
    #[test]
    fn test_debye_length_decreases_with_concentration() {
        let c1 = DoubleLayerCapacitance::new(10.0, 78.5 * 8.854e-12, 298.15);
        let c2 = DoubleLayerCapacitance::new(100.0, 78.5 * 8.854e-12, 298.15);
        assert!(
            c2.debye_length() < c1.debye_length(),
            "Debye length decreases with concentration"
        );
    }
    #[test]
    fn test_double_layer_capacitance_at_pzc_positive() {
        let dl = DoubleLayerCapacitance::kcl_01_mol_l();
        let c = dl.capacitance_at_pzc();
        assert!(c > 0.0, "Capacitance at PZC should be positive, got {c}");
    }
    #[test]
    fn test_differential_capacitance_minimum_at_pzc() {
        let dl = DoubleLayerCapacitance::kcl_01_mol_l();
        let c_pzc = dl.differential_capacitance(0.0);
        let c_off = dl.differential_capacitance(0.2);
        assert!(
            c_off >= c_pzc,
            "Capacitance minimum at PZC: c_pzc={c_pzc}, c_off={c_off}"
        );
    }
    #[test]
    fn test_differential_capacitance_even_function() {
        let dl = DoubleLayerCapacitance::kcl_01_mol_l();
        let c_pos = dl.differential_capacitance(0.1);
        let c_neg = dl.differential_capacitance(-0.1);
        assert!(
            (c_pos - c_neg).abs() < 1e-20,
            "Differential capacitance is even: {c_pos} vs {c_neg}"
        );
    }
}
