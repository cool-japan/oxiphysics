// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Energy materials module: Li-ion batteries, fuel cells, solar cells,
//! thermoelectrics, supercapacitors, hydrogen storage, nuclear materials,
//! and piezoelectric energy harvesting.
//!
//! All functions use SI units unless otherwise stated.

#![allow(dead_code)]
#![allow(unused_imports)]
#![allow(clippy::too_many_arguments)]

use std::f64::consts::PI;

/// Universal gas constant (J/mol·K).
const R: f64 = 8.314;
/// Faraday constant (C/mol).
const F_FARADAY: f64 = 96_485.0;
/// Boltzmann constant (J/K).
const K_B: f64 = 1.380649e-23;
/// Elementary charge (C).
const Q_E: f64 = 1.602176634e-19;
/// Planck constant (J·s).
const H_PLANCK: f64 = 6.62607015e-34;
/// Speed of light (m/s).
const C_LIGHT: f64 = 2.99792458e8;

// ---------------------------------------------------------------------------
// Lithium-Ion Battery — Butler-Volmer Kinetics & Diffusion
// ---------------------------------------------------------------------------

/// Parameters for a Li-ion battery electrode (half-cell).
#[derive(Clone, Debug)]
pub struct LiIonElectrode {
    /// Exchange current density i₀ (A/m²).
    pub i0: f64,
    /// Transfer coefficient alpha (dimensionless, typically 0.5).
    pub alpha: f64,
    /// Temperature (K).
    pub temperature: f64,
    /// Electrode thickness (m).
    pub thickness: f64,
    /// Solid-phase diffusion coefficient Ds (m²/s).
    pub d_solid: f64,
    /// Electrolyte diffusion coefficient De (m²/s).
    pub d_electrolyte: f64,
    /// Maximum solid concentration cs_max (mol/m³).
    pub cs_max: f64,
    /// Electrode volume fraction of active material.
    pub epsilon_s: f64,
    /// Particle radius (m).
    pub particle_radius: f64,
    /// Open circuit potential function (stoichiometry -> V).
    pub ocp_a: f64,
    /// OCP slope.
    pub ocp_b: f64,
}

impl LiIonElectrode {
    /// Typical LiCoO₂ cathode defaults.
    pub fn lco_cathode() -> Self {
        Self {
            i0: 2.0,
            alpha: 0.5,
            temperature: 298.15,
            thickness: 80.0e-6,
            d_solid: 1.0e-14,
            d_electrolyte: 7.5e-10,
            cs_max: 51_410.0,
            epsilon_s: 0.5,
            particle_radius: 2.0e-6,
            ocp_a: 4.2,
            ocp_b: -0.5,
        }
    }

    /// Graphite anode defaults.
    pub fn graphite_anode() -> Self {
        Self {
            i0: 3.6,
            alpha: 0.5,
            temperature: 298.15,
            thickness: 100.0e-6,
            d_solid: 3.9e-14,
            d_electrolyte: 7.5e-10,
            cs_max: 26_390.0,
            epsilon_s: 0.471,
            particle_radius: 5.0e-6,
            ocp_a: 0.2,
            ocp_b: 0.1,
        }
    }

    /// Butler-Volmer current density (A/m²) as function of overpotential eta (V).
    pub fn butler_volmer(&self, eta: f64) -> f64 {
        let f_rt = F_FARADAY / (R * self.temperature);
        self.i0 * ((self.alpha * f_rt * eta).exp() - ((self.alpha - 1.0) * f_rt * eta).exp())
    }

    /// Linearised BV (low overpotential approximation): i ≈ i₀ * F/(R*T) * eta.
    pub fn linear_butler_volmer(&self, eta: f64) -> f64 {
        self.i0 * F_FARADAY / (R * self.temperature) * eta
    }

    /// Open circuit potential (OCP) as function of state-of-charge x (0..1).
    pub fn ocp(&self, x: f64) -> f64 {
        self.ocp_a + self.ocp_b * x
    }

    /// Solid diffusion time constant (s): tau = R_p^2 / (15 * Ds).
    pub fn diffusion_time_constant(&self) -> f64 {
        self.particle_radius * self.particle_radius / (15.0 * self.d_solid)
    }

    /// Effective reaction surface area per unit volume (m⁻¹).
    pub fn specific_interfacial_area(&self) -> f64 {
        3.0 * self.epsilon_s / self.particle_radius
    }

    /// Capacity per unit area (A·s/m²) at 1C rate.
    pub fn areal_capacity(&self) -> f64 {
        self.cs_max * F_FARADAY * self.epsilon_s * self.thickness
    }

    /// State of charge from surface concentration cs_surf / cs_max.
    pub fn state_of_charge(&self, cs_surf: f64) -> f64 {
        (cs_surf / self.cs_max).clamp(0.0, 1.0)
    }

    /// Tafel slope (V/decade) for the electrode reaction.
    pub fn tafel_slope(&self) -> f64 {
        R * self.temperature / (self.alpha * F_FARADAY) * 10.0f64.ln()
    }
}

/// Battery cell model combining anode and cathode.
#[derive(Clone, Debug)]
pub struct LiIonCell {
    /// Cathode electrode.
    pub cathode: LiIonElectrode,
    /// Anode electrode.
    pub anode: LiIonElectrode,
    /// Electrolyte resistance (Ohm·m²).
    pub r_electrolyte: f64,
    /// Contact resistance (Ohm·m²).
    pub r_contact: f64,
    /// Nominal capacity (A·h).
    pub capacity_ah: f64,
}

impl LiIonCell {
    /// Standard 18650 LCO/graphite cell.
    pub fn nmc18650() -> Self {
        Self {
            cathode: LiIonElectrode::lco_cathode(),
            anode: LiIonElectrode::graphite_anode(),
            r_electrolyte: 5.0e-5,
            r_contact: 5.0e-6,
            capacity_ah: 2.5,
        }
    }

    /// Cell open circuit voltage.
    pub fn ocv(&self, soc_cathode: f64, soc_anode: f64) -> f64 {
        self.cathode.ocp(soc_cathode) - self.anode.ocp(soc_anode)
    }

    /// Terminal voltage during discharge at current density j (A/m²).
    pub fn terminal_voltage(&self, j: f64, soc: f64) -> f64 {
        let ocv = self.cathode.ocp(soc) - self.anode.ocp(1.0 - soc);
        let eta_c = self.cathode.butler_volmer(j / self.cathode.i0);
        let eta_a = self.anode.butler_volmer(-j / self.anode.i0);
        let v_ir = j * (self.r_electrolyte + self.r_contact);
        ocv - eta_c + eta_a - v_ir
    }

    /// C-rate given current I (A) and capacity.
    pub fn c_rate(&self, current_a: f64) -> f64 {
        current_a / self.capacity_ah
    }

    /// Ragone plot: specific energy (Wh/kg) vs specific power (W/kg).
    pub fn ragone_point(&self, c_rate: f64, cell_mass_kg: f64) -> (f64, f64) {
        let current = c_rate * self.capacity_ah;
        let voltage = 3.6; // Average discharge voltage
        let power = current * voltage;
        let energy = self.capacity_ah * voltage;
        (energy / cell_mass_kg, power / cell_mass_kg)
    }
}

// ---------------------------------------------------------------------------
// Fuel Cell — Nafion Proton Conductivity and GDL Transport
// ---------------------------------------------------------------------------

/// Nafion membrane properties for PEM fuel cells.
#[derive(Clone, Debug)]
pub struct NafionMembrane {
    /// Membrane thickness (m).
    pub thickness: f64,
    /// Water content lambda (H₂O molecules per SO₃⁻ group).
    pub water_content: f64,
    /// Temperature (K).
    pub temperature: f64,
    /// Equivalent weight (g/mol SO₃⁻).
    pub eq_weight: f64,
    /// Dry density (kg/m³).
    pub dry_density: f64,
}

impl NafionMembrane {
    /// Nafion 117 default parameters.
    pub fn nafion117() -> Self {
        Self {
            thickness: 183.0e-6,
            water_content: 14.0,
            temperature: 353.15,
            eq_weight: 1100.0,
            dry_density: 2100.0,
        }
    }

    /// Springer proton conductivity (S/m) as function of temperature and water content.
    pub fn proton_conductivity(&self) -> f64 {
        let lambda = self.water_content;
        let t = self.temperature;
        // Springer et al. (1991)
        (0.005139 * lambda - 0.00326) * (1268.0 * (1.0 / 303.15 - 1.0 / t)).exp() * 1000.0 // S/m conversion from S/cm
    }

    /// Membrane resistance (Ohm·m²).
    pub fn area_specific_resistance(&self) -> f64 {
        self.thickness / self.proton_conductivity()
    }

    /// Water diffusivity in Nafion (m²/s) — Motupally correlation.
    pub fn water_diffusivity(&self) -> f64 {
        let lambda = self.water_content;
        let t = self.temperature;
        let d_w = if lambda < 3.0 {
            1e-10 * (0.87 * (3.0 - lambda) + 2.95 * lambda)
        } else if lambda < 17.0 {
            1e-10 * (2.95 + 1.00 * (lambda - 3.0))
        } else {
            1e-10 * (1.54 + 0.198 * lambda)
        };
        d_w * (2416.0 * (1.0 / 303.15 - 1.0 / t)).exp()
    }

    /// Sorption isotherm: equilibrium lambda vs. relative humidity a_w.
    pub fn sorption_isotherm(&self, a_w: f64) -> f64 {
        // Springer model
        0.043 + 17.81 * a_w - 39.85 * a_w * a_w + 36.0 * a_w * a_w * a_w
    }

    /// Electroosmotic drag coefficient (H₂O/H⁺) at given lambda.
    pub fn electroosmotic_drag(&self) -> f64 {
        2.5 * self.water_content / 22.0
    }
}

/// Gas diffusion layer (GDL) transport model.
#[derive(Clone, Debug)]
pub struct GasDiffusionLayer {
    /// Porosity (0..1).
    pub porosity: f64,
    /// Tortuosity.
    pub tortuosity: f64,
    /// Thickness (m).
    pub thickness: f64,
    /// Contact angle (rad) — for liquid water flooding.
    pub contact_angle: f64,
    /// Permeability (m²).
    pub permeability: f64,
}

impl GasDiffusionLayer {
    /// Toray carbon paper TGP-H-120 defaults.
    pub fn toray_tgp_h120() -> Self {
        Self {
            porosity: 0.78,
            tortuosity: 1.5,
            thickness: 370.0e-6,
            contact_angle: 110.0f64.to_radians(),
            permeability: 7.5e-12,
        }
    }

    /// Effective diffusivity via Bruggeman model: D_eff = D_bulk * eps^1.5.
    pub fn effective_diffusivity(&self, d_bulk: f64) -> f64 {
        d_bulk * self.porosity.powf(1.5)
    }

    /// Knudsen diffusivity for small pores (m²/s).
    ///
    /// `d_pore`: mean pore diameter (m). `m_molar`: molar mass (kg/mol).
    pub fn knudsen_diffusivity(&self, d_pore: f64, m_molar: f64, temperature: f64) -> f64 {
        (d_pore / 3.0) * (8.0 * R * temperature / (PI * m_molar)).sqrt()
    }

    /// Darcy velocity (m/s) under pressure gradient dp_dx (Pa/m).
    pub fn darcy_velocity(&self, viscosity: f64, dp_dx: f64) -> f64 {
        self.permeability / viscosity * dp_dx
    }

    /// Capillary pressure (Pa) for a pore of diameter d_pore.
    pub fn capillary_pressure(&self, d_pore: f64, surface_tension: f64) -> f64 {
        4.0 * surface_tension * self.contact_angle.cos().abs() / d_pore
    }
}

/// Polarisation curve model for a PEM fuel cell.
#[derive(Clone, Debug)]
pub struct PemFuelCell {
    /// Nafion membrane.
    pub membrane: NafionMembrane,
    /// Cathode GDL.
    pub cathode_gdl: GasDiffusionLayer,
    /// Anode GDL.
    pub anode_gdl: GasDiffusionLayer,
    /// Cathode exchange current density (A/m²).
    pub i0_c: f64,
    /// Anode exchange current density (A/m²).
    pub i0_a: f64,
    /// Cathode transfer coefficient.
    pub alpha_c: f64,
    /// Nernst open circuit voltage (V) at standard conditions.
    pub e_nernst: f64,
    /// Temperature (K).
    pub temperature: f64,
}

impl PemFuelCell {
    /// Standard hydrogen/air PEM fuel cell.
    pub fn standard() -> Self {
        Self {
            membrane: NafionMembrane::nafion117(),
            cathode_gdl: GasDiffusionLayer::toray_tgp_h120(),
            anode_gdl: GasDiffusionLayer::toray_tgp_h120(),
            i0_c: 1.0e-4,
            i0_a: 10.0,
            alpha_c: 0.5,
            e_nernst: 1.229,
            temperature: 353.15,
        }
    }

    /// Activation overpotential at cathode (V) for current density j (A/m²).
    pub fn cathode_activation_overpotential(&self, j: f64) -> f64 {
        R * self.temperature / (self.alpha_c * F_FARADAY) * (j / self.i0_c).abs().ln()
    }

    /// Ohmic overpotential (V) for current density j (A/m²).
    pub fn ohmic_overpotential(&self, j: f64) -> f64 {
        j * self.membrane.area_specific_resistance()
    }

    /// Limiting current density (A/m²) for oxygen transport.
    pub fn limiting_current_density(&self, c_bulk: f64, d_eff: f64) -> f64 {
        4.0 * F_FARADAY * d_eff * c_bulk / self.cathode_gdl.thickness
    }

    /// Cell voltage (V) at current density j.
    pub fn cell_voltage(&self, j: f64) -> f64 {
        let eta_act = self.cathode_activation_overpotential(j);
        let eta_ohm = self.ohmic_overpotential(j);
        (self.e_nernst - eta_act - eta_ohm).max(0.0)
    }

    /// Power density (W/m²) at current density j.
    pub fn power_density(&self, j: f64) -> f64 {
        j * self.cell_voltage(j)
    }
}

// ---------------------------------------------------------------------------
// Solar Cell — Shockley-Queisser Limit & p-n Junction
// ---------------------------------------------------------------------------

/// Solar cell parameters (single-junction diode model).
#[derive(Clone, Debug)]
pub struct SolarCell {
    /// Short-circuit current density Jsc (A/m²).
    pub j_sc: f64,
    /// Reverse saturation current density J0 (A/m²).
    pub j0: f64,
    /// Ideality factor n (1..2).
    pub ideality: f64,
    /// Series resistance (Ohm·m²).
    pub r_series: f64,
    /// Shunt resistance (Ohm·m²).
    pub r_shunt: f64,
    /// Temperature (K).
    pub temperature: f64,
    /// Bandgap energy Eg (eV).
    pub bandgap_ev: f64,
}

impl SolarCell {
    /// Create a new solar cell with given parameters.
    ///
    /// # Arguments
    /// * `j_sc` – Short-circuit current density (A/m²)
    /// * `j0` – Reverse saturation current density (A/m²)
    /// * `ideality` – Ideality factor (1–2)
    /// * `temperature` – Temperature (K)
    /// * `r_series` – Series resistance (Ω·m²)
    pub fn new(j_sc: f64, j0: f64, ideality: f64, temperature: f64, r_series: f64) -> Self {
        Self {
            j_sc,
            j0,
            ideality,
            r_series,
            r_shunt: 1.0e4,
            temperature,
            bandgap_ev: 1.12,
        }
    }

    /// Typical silicon solar cell (AM1.5 illumination).
    pub fn silicon() -> Self {
        Self {
            j_sc: 400.0,
            j0: 1.0e-10,
            ideality: 1.0,
            r_series: 1.0e-5,
            r_shunt: 1.0e3,
            temperature: 298.15,
            bandgap_ev: 1.12,
        }
    }

    /// GaAs solar cell.
    pub fn gaas() -> Self {
        Self {
            j_sc: 290.0,
            j0: 1.0e-15,
            ideality: 1.0,
            r_series: 5.0e-6,
            r_shunt: 5.0e3,
            temperature: 298.15,
            bandgap_ev: 1.42,
        }
    }

    /// Thermal voltage Vt = kT/q (V).
    pub fn thermal_voltage(&self) -> f64 {
        K_B * self.temperature / Q_E
    }

    /// Diode current density (A/m²) at voltage V (including parasitic resistances).
    pub fn current_density(&self, v: f64) -> f64 {
        let vt = self.thermal_voltage() * self.ideality;
        // Implicit equation; approximate solution ignoring r_series coupling
        let j_diode = self.j0 * ((v / vt).exp() - 1.0);
        let j_shunt = v / self.r_shunt;
        self.j_sc - j_diode - j_shunt
    }

    /// Open circuit voltage Voc (V).
    pub fn voc(&self) -> f64 {
        let vt = self.thermal_voltage() * self.ideality;
        vt * (self.j_sc / self.j0 + 1.0).ln()
    }

    /// Estimate fill factor (FF) using Ro = Jsc * Voc / J0-normalised formula.
    pub fn fill_factor(&self) -> f64 {
        let voc_norm = self.voc() / (self.thermal_voltage() * self.ideality);
        // Empirical formula (Green, 1982)
        (voc_norm - (voc_norm + 0.72).ln()) / (voc_norm + 1.0)
    }

    /// Maximum power point voltage Vmpp (V).
    pub fn v_mpp(&self) -> f64 {
        let voc = self.voc();
        let vt = self.thermal_voltage() * self.ideality;
        // Approximate: Vmpp ≈ Voc - Vt * ln(Vmpp / Vt + 1)
        // Iterate
        let mut v = voc * 0.8;
        for _ in 0..20 {
            v = voc - vt * (v / vt + 1.0).ln();
        }
        v
    }

    /// Maximum power density (W/m²).
    pub fn p_max(&self) -> f64 {
        let vmpp = self.v_mpp();
        let jmpp = self.current_density(vmpp);
        vmpp * jmpp.max(0.0)
    }

    /// Power conversion efficiency at AM1.5 (1000 W/m² irradiance).
    pub fn efficiency(&self, irradiance: f64) -> f64 {
        self.p_max() / irradiance
    }

    /// Shockley-Queisser detailed balance limit efficiency for bandgap Eg (eV).
    ///
    /// Uses a simplified analytical approximation.
    pub fn sq_limit(bandgap_ev: f64) -> f64 {
        // Approximate SQ limit (rough numerical estimate)
        let eg = bandgap_ev;
        if !(0.5..=4.0).contains(&eg) {
            return 0.0;
        }
        // Empirical curve fit

        0.31 * (1.0 - ((eg - 1.34) / 0.8).powi(2)).max(0.0)
    }

    /// Temperature coefficient of Voc (V/K) — typically negative.
    pub fn temp_coeff_voc(&self) -> f64 {
        -K_B / Q_E * (self.j_sc / self.j0).ln() / self.temperature
    }
}

// ---------------------------------------------------------------------------
// Thermoelectrics — Seebeck, Peltier, Thomson, ZT Figure of Merit
// ---------------------------------------------------------------------------

/// Thermoelectric material properties.
#[derive(Clone, Debug)]
pub struct Thermoelectric {
    /// Seebeck coefficient S (V/K).
    pub seebeck: f64,
    /// Electrical conductivity sigma (S/m).
    pub electrical_conductivity: f64,
    /// Thermal conductivity kappa (W/m·K).
    pub thermal_conductivity: f64,
    /// Operating temperature T (K).
    pub temperature: f64,
}

impl Thermoelectric {
    /// Bismuth telluride Bi₂Te₃ at room temperature.
    pub fn bi2te3() -> Self {
        Self {
            seebeck: 200.0e-6,
            electrical_conductivity: 1.0e5,
            thermal_conductivity: 1.0,
            temperature: 300.0,
        }
    }

    /// PbTe at high temperature.
    pub fn pbte_high_temp() -> Self {
        Self {
            seebeck: 300.0e-6,
            electrical_conductivity: 7.0e4,
            thermal_conductivity: 1.5,
            temperature: 750.0,
        }
    }

    /// Dimensionless figure of merit ZT = S² * sigma * T / kappa.
    pub fn zt(&self) -> f64 {
        self.seebeck * self.seebeck * self.electrical_conductivity * self.temperature
            / self.thermal_conductivity
    }

    /// Power factor PF = S² * sigma (W/m·K²).
    pub fn power_factor(&self) -> f64 {
        self.seebeck * self.seebeck * self.electrical_conductivity
    }

    /// Carnot efficiency limit (dimensionless).
    pub fn carnot_efficiency(&self, t_cold: f64) -> f64 {
        1.0 - t_cold / self.temperature
    }

    /// Maximum thermoelectric efficiency (generator mode).
    pub fn max_efficiency(&self, t_hot: f64, t_cold: f64) -> f64 {
        let zt_avg = self.zt(); // Using T_hot as reference
        let sqrt_term = (1.0 + zt_avg).sqrt();
        let delta_t = t_hot - t_cold;
        (delta_t / t_hot) * (sqrt_term - 1.0) / (sqrt_term + t_cold / t_hot)
    }

    /// Peltier coefficient Pi = S * T (W/A = V).
    pub fn peltier_coefficient(&self) -> f64 {
        self.seebeck * self.temperature
    }

    /// Thomson coefficient tau = T * dS/dT (V/K).
    ///
    /// Approximated using finite difference of S over a small dT.
    pub fn thomson_coefficient(&self, ds_dt: f64) -> f64 {
        self.temperature * ds_dt
    }

    /// Maximum COP for Peltier cooler.
    pub fn max_cop_cooling(&self, t_hot: f64, t_cold: f64) -> f64 {
        let zt_avg = self.zt();
        let sqrt_term = (1.0 + zt_avg).sqrt();
        let delta_t = t_hot - t_cold;
        (t_cold * (sqrt_term - t_hot / t_cold)) / (delta_t * (sqrt_term + 1.0))
    }

    /// Open circuit voltage for temperature difference dt.
    pub fn seebeck_voltage(&self, dt: f64) -> f64 {
        self.seebeck * dt
    }

    /// Thermal resistance (K·m/W) per unit cross-sectional area and length.
    pub fn thermal_resistance(&self, length: f64) -> f64 {
        length / self.thermal_conductivity
    }
}

// ---------------------------------------------------------------------------
// Supercapacitor — EDL and Pseudocapacitance
// ---------------------------------------------------------------------------

/// Electrochemical double layer (EDL) supercapacitor model.
#[derive(Clone, Debug)]
pub struct Supercapacitor {
    /// Specific capacitance per unit area (F/m²).
    pub c_specific: f64,
    /// Electrode surface area (m²).
    pub electrode_area: f64,
    /// Electrolyte resistance (Ohm).
    pub r_esr: f64,
    /// Leakage resistance (Ohm).
    pub r_leakage: f64,
    /// Maximum voltage (V).
    pub v_max: f64,
    /// Pseudocapacitance fraction (0..1).
    pub pseudo_fraction: f64,
    /// Temperature (K).
    pub temperature: f64,
}

impl Supercapacitor {
    /// Activated carbon EDL supercapacitor (500 F).
    pub fn activated_carbon_500f() -> Self {
        Self {
            c_specific: 0.2,
            electrode_area: 2500.0,
            r_esr: 1.0e-3,
            r_leakage: 1.0e5,
            v_max: 2.7,
            pseudo_fraction: 0.0,
            temperature: 298.15,
        }
    }

    /// RuO₂ pseudocapacitor.
    pub fn ruo2_pseudocap() -> Self {
        Self {
            c_specific: 0.8,
            electrode_area: 500.0,
            r_esr: 5.0e-3,
            r_leakage: 5.0e4,
            v_max: 1.2,
            pseudo_fraction: 0.6,
            temperature: 298.15,
        }
    }

    /// Total capacitance (F).
    pub fn capacitance(&self) -> f64 {
        self.c_specific * self.electrode_area * (1.0 + self.pseudo_fraction)
    }

    /// Stored energy (J) at voltage V.
    pub fn energy_stored(&self, v: f64) -> f64 {
        0.5 * self.capacitance() * v * v
    }

    /// Maximum stored energy (J).
    pub fn max_energy(&self) -> f64 {
        self.energy_stored(self.v_max)
    }

    /// Peak power (W) limited by ESR.
    pub fn peak_power(&self) -> f64 {
        self.v_max * self.v_max / (4.0 * self.r_esr)
    }

    /// Specific energy (Wh/kg) assuming electrode mass `m_kg`.
    pub fn specific_energy_wh_kg(&self, m_kg: f64) -> f64 {
        self.max_energy() / (m_kg * 3600.0)
    }

    /// RC time constant (s).
    pub fn time_constant(&self) -> f64 {
        self.r_esr * self.capacitance()
    }

    /// Voltage decay under constant leakage: V(t) = V0 * exp(-t / (R_leak * C)).
    pub fn voltage_decay(&self, v0: f64, t: f64) -> f64 {
        let tau_leak = self.r_leakage * self.capacitance();
        v0 * (-t / tau_leak).exp()
    }

    /// Charge stored (C) at voltage V.
    pub fn charge_stored(&self, v: f64) -> f64 {
        self.capacitance() * v
    }

    /// Helmholtz EDL capacitance per area (F/m²) from Stern model.
    pub fn helmholtz_capacitance(epsilon_r: f64, d_debye: f64) -> f64 {
        8.854187817e-12 * epsilon_r / d_debye
    }

    /// Debye screening length (m).
    pub fn debye_length(c_mol_l: f64, z: f64, temperature: f64) -> f64 {
        let c_mol_m3 = c_mol_l * 1000.0;
        let eps = 8.854e-12 * 80.0; // Water dielectric constant
        (eps * R * temperature / (2.0 * F_FARADAY * F_FARADAY * c_mol_m3 * z * z)).sqrt()
    }
}

// ---------------------------------------------------------------------------
// Hydrogen Storage — Metal Hydride Kinetics
// ---------------------------------------------------------------------------

/// Metal hydride hydrogen storage model.
#[derive(Clone, Debug)]
pub struct MetalHydride {
    /// Maximum hydrogen storage capacity (wt%).
    pub max_capacity_wt: f64,
    /// Activation energy for absorption Ea (J/mol).
    pub ea_absorption: f64,
    /// Activation energy for desorption Ea (J/mol).
    pub ea_desorption: f64,
    /// Pre-exponential factor ka (1/s).
    pub k0_absorption: f64,
    /// Pre-exponential factor kd (1/s).
    pub k0_desorption: f64,
    /// Equilibrium plateau pressure at 298 K (Pa).
    pub p_eq_298: f64,
    /// Enthalpy of hydride formation Delta_H (J/mol H₂).
    pub delta_h: f64,
    /// Entropy of hydride formation Delta_S (J/mol H₂/K).
    pub delta_s: f64,
    /// Temperature (K).
    pub temperature: f64,
}

impl MetalHydride {
    /// LaNi₅ hydride defaults.
    pub fn lani5() -> Self {
        Self {
            max_capacity_wt: 1.4,
            ea_absorption: 21_000.0,
            ea_desorption: 35_000.0,
            k0_absorption: 500.0,
            k0_desorption: 200.0,
            p_eq_298: 200_000.0,
            delta_h: -30_800.0,
            delta_s: -108.0,
            temperature: 298.15,
        }
    }

    /// MgH₂ high-capacity hydride.
    pub fn mgh2() -> Self {
        Self {
            max_capacity_wt: 7.6,
            ea_absorption: 60_000.0,
            ea_desorption: 80_000.0,
            k0_absorption: 100.0,
            k0_desorption: 50.0,
            p_eq_298: 0.01,
            delta_h: -75_000.0,
            delta_s: -135.0,
            temperature: 623.15,
        }
    }

    /// Van't Hoff equilibrium pressure at temperature T (Pa).
    pub fn equilibrium_pressure(&self, t: f64) -> f64 {
        let p0 = 1.0e5; // Reference pressure (Pa)
        p0 * (self.delta_h / (R * t) - self.delta_s / R).exp()
    }

    /// Absorption rate constant (1/s) at temperature T.
    pub fn absorption_rate_constant(&self, t: f64) -> f64 {
        self.k0_absorption * (-self.ea_absorption / (R * t)).exp()
    }

    /// Desorption rate constant (1/s) at temperature T.
    pub fn desorption_rate_constant(&self, t: f64) -> f64 {
        self.k0_desorption * (-self.ea_desorption / (R * t)).exp()
    }

    /// Reaction kinetics (d\[H\]/dt) (wt%/s) during absorption.
    ///
    /// `x`: current H content (wt%). `p`: applied H₂ pressure (Pa).
    pub fn absorption_rate(&self, x: f64, p: f64) -> f64 {
        let peq = self.equilibrium_pressure(self.temperature);
        let ka = self.absorption_rate_constant(self.temperature);
        if p <= peq {
            return 0.0;
        }
        ka * (self.max_capacity_wt - x) * ((p / peq).ln()).max(0.0)
    }

    /// Reaction kinetics (d\[H\]/dt) (wt%/s) during desorption.
    pub fn desorption_rate(&self, x: f64, p: f64) -> f64 {
        let peq = self.equilibrium_pressure(self.temperature);
        let kd = self.desorption_rate_constant(self.temperature);
        if p >= peq {
            return 0.0;
        }
        kd * x * ((peq / p.max(1e-10)).ln()).max(0.0)
    }

    /// Heat generated during absorption (W/kg H₂).
    pub fn absorption_heat(&self, rate_wt_per_s: f64) -> f64 {
        // Q = -Delta_H / M_H2 * rate
        rate_wt_per_s * (-self.delta_h) / 2.016 // kJ/kg approximate
    }
}

// ---------------------------------------------------------------------------
// Nuclear Materials — Radiation Damage and Swelling
// ---------------------------------------------------------------------------

/// Nuclear material radiation damage model.
#[derive(Clone, Debug)]
pub struct NuclearMaterial {
    /// Material name.
    pub name: String,
    /// Threshold displacement energy Ed (eV).
    pub e_displacement: f64,
    /// Swelling rate (%/dpa) at operating temperature.
    pub swelling_rate: f64,
    /// Irradiation hardening coefficient (MPa / sqrt(dpa)).
    pub hardening_coeff: f64,
    /// Creep compliance coefficient B (per Pa per dpa).
    pub creep_coefficient: f64,
    /// Embrittlement transition temperature shift (K/dpa).
    pub dbtt_shift_per_dpa: f64,
    /// Initial yield stress (Pa).
    pub yield_stress_0: f64,
    /// Young's modulus (Pa).
    pub modulus: f64,
}

impl NuclearMaterial {
    /// Ferritic-martensitic steel (F82H) for fusion applications.
    pub fn f82h() -> Self {
        Self {
            name: "F82H ferritic-martensitic steel".into(),
            e_displacement: 40.0,
            swelling_rate: 0.002,
            hardening_coeff: 120.0e6,
            creep_coefficient: 1.0e-6,
            dbtt_shift_per_dpa: 1.5,
            yield_stress_0: 500.0e6,
            modulus: 210.0e9,
        }
    }

    /// Austenitic stainless steel 316L.
    pub fn ss316l() -> Self {
        Self {
            name: "316L stainless steel".into(),
            e_displacement: 40.0,
            swelling_rate: 0.05,
            hardening_coeff: 200.0e6,
            creep_coefficient: 2.0e-6,
            dbtt_shift_per_dpa: 0.0,
            yield_stress_0: 270.0e6,
            modulus: 193.0e9,
        }
    }

    /// Displacement damage (dpa) from neutron fluence.
    ///
    /// Uses the NRT (Norgett-Robinson-Torrens) standard.
    pub fn dpa_from_fluence(&self, fluence: f64, sigma_d: f64) -> f64 {
        // dpa = phi * sigma_d * t (simplified)
        fluence * sigma_d
    }

    /// Radiation-induced swelling (volume fraction).
    pub fn swelling(&self, dpa: f64) -> f64 {
        // Incubation period (~1 dpa) then linear
        if dpa < 1.0 {
            return 0.0;
        }
        self.swelling_rate * (dpa - 1.0) / 100.0
    }

    /// Irradiation hardening: delta_sigma_y = k * sqrt(dpa).
    pub fn irradiation_hardening(&self, dpa: f64) -> f64 {
        self.hardening_coeff * dpa.sqrt()
    }

    /// Yield stress after irradiation.
    pub fn irradiated_yield_stress(&self, dpa: f64) -> f64 {
        self.yield_stress_0 + self.irradiation_hardening(dpa)
    }

    /// Irradiation creep strain rate (1/s) under stress sigma and flux phi (dpa/s).
    pub fn irradiation_creep_rate(&self, sigma: f64, dpas_per_s: f64) -> f64 {
        self.creep_coefficient * sigma * dpas_per_s
    }

    /// Ductile-to-brittle transition temperature shift (K).
    pub fn dbtt_shift(&self, dpa: f64) -> f64 {
        self.dbtt_shift_per_dpa * dpa
    }

    /// PKA (primary knock-on atom) energy from Lindhard theory.
    ///
    /// `e_recoil`: recoil energy (eV). Returns effective PKA energy (eV).
    pub fn pka_lindhard(&self, e_recoil: f64, atomic_mass: f64) -> f64 {
        let xi = 0.16 * atomic_mass.powf(-2.0 / 3.0);
        e_recoil / (1.0 + xi * e_recoil / self.e_displacement)
    }

    /// NRT displacement cross-section (barn, simplified).
    pub fn nrt_displacements(&self, e_pka: f64) -> f64 {
        if e_pka < self.e_displacement {
            return 0.0;
        }
        0.8 * e_pka / (2.0 * self.e_displacement)
    }
}

// ---------------------------------------------------------------------------
// Piezoelectric Energy Harvesting
// ---------------------------------------------------------------------------

/// Piezoelectric energy harvesting model (bimorph cantilever).
#[derive(Clone, Debug)]
pub struct PiezoHarvester {
    /// Piezoelectric coupling coefficient d₃₁ (C/N = m/V).
    pub d31: f64,
    /// Electromechanical coupling factor k₃₁ (dimensionless).
    pub k31: f64,
    /// Young's modulus of piezo material (Pa).
    pub e_piezo: f64,
    /// Density of piezo material (kg/m³).
    pub density: f64,
    /// Beam length L (m).
    pub length: f64,
    /// Beam width b (m).
    pub width: f64,
    /// Beam thickness h (m).
    pub thickness: f64,
    /// Permittivity at constant stress epsilon_T (F/m).
    pub permittivity: f64,
    /// Mechanical damping ratio zeta.
    pub damping: f64,
    /// Proof mass (kg).
    pub proof_mass: f64,
}

impl PiezoHarvester {
    /// PZT-5A bimorph cantilever defaults.
    pub fn pzt5a() -> Self {
        Self {
            d31: -175.0e-12,
            k31: 0.34,
            e_piezo: 61.0e9,
            density: 7750.0,
            length: 30.0e-3,
            width: 5.0e-3,
            thickness: 0.5e-3,
            permittivity: 15.0e-9,
            damping: 0.02,
            proof_mass: 1.0e-3,
        }
    }

    /// Natural (resonant) frequency f_n (Hz) for a bimorph cantilever.
    pub fn natural_frequency(&self) -> f64 {
        let ei = self.e_piezo * self.width * self.thickness.powi(3) / 12.0;
        let m = self.proof_mass + 0.23 * self.density * self.width * self.thickness * self.length;
        let k = 3.0 * ei / self.length.powi(3);
        (k / m).sqrt() / (2.0 * PI)
    }

    /// Open circuit voltage amplitude V_oc (V) at base acceleration a (m/s²) and excitation freq f.
    pub fn open_circuit_voltage(&self, a: f64, f: f64) -> f64 {
        let fn_ = self.natural_frequency();
        let omega = 2.0 * PI * f;
        let omega_n = 2.0 * PI * fn_;
        let r = omega / omega_n;
        let m = self.proof_mass;
        // Tip displacement amplitude
        let denom = ((1.0 - r * r).powi(2) + (2.0 * self.damping * r).powi(2)).sqrt();
        let x_amp =
            m * a / (self.e_piezo * self.width * self.thickness * denom / self.length.powi(2));
        // Voltage from piezoelectric coupling
        let h_pc = self.e_piezo * self.d31.abs() * 3.0 * self.thickness
            / (2.0 * self.length * self.length);
        let c_p = self.permittivity * self.width * self.length / self.thickness;
        h_pc * x_amp / c_p
    }

    /// Maximum power output (W) at resonance into matched load.
    pub fn max_power(&self, a: f64) -> f64 {
        let fn_ = self.natural_frequency();
        let voc = self.open_circuit_voltage(a, fn_);
        let c_p = self.permittivity * self.width * self.length / self.thickness;
        let omega_n = 2.0 * PI * fn_;
        let r_opt = 1.0 / (omega_n * c_p);
        voc * voc / (4.0 * r_opt)
    }

    /// Normalised power density (W/m³) for volume Vol.
    pub fn power_density(&self, a: f64) -> f64 {
        let vol = self.length * self.width * self.thickness;
        self.max_power(a) / vol
    }

    /// Effective piezoelectric coefficient (C/m²) g₃₁ = d₃₁ / epsilon_T.
    pub fn g31(&self) -> f64 {
        self.d31 / self.permittivity
    }

    /// Mechanical quality factor Qm ≈ 1 / (2 * zeta).
    pub fn quality_factor(&self) -> f64 {
        1.0 / (2.0 * self.damping)
    }

    /// Efficiency of electromechanical coupling: k²/(1+k²), in \[0,1).
    pub fn coupling_efficiency(&self) -> f64 {
        let k2 = self.k31 * self.k31;
        k2 / (1.0 + k2)
    }
}

// ---------------------------------------------------------------------------
// Electrode material capacity and voltage models
// ---------------------------------------------------------------------------

/// Capacity fading model for Li-ion batteries under cycling.
#[derive(Clone, Debug)]
pub struct CapacityFadingModel {
    /// Initial capacity (A·h).
    pub q0: f64,
    /// SEI growth rate constant k_sei (s⁻¹·⁰·⁵).
    pub k_sei: f64,
    /// Calendar aging coefficient k_cal.
    pub k_cal: f64,
    /// Cycle aging exponent.
    pub cycle_exp: f64,
}

impl CapacityFadingModel {
    /// Typical NMC 18650 cell aging model.
    pub fn nmc18650() -> Self {
        Self {
            q0: 3.0,
            k_sei: 1.5e-4,
            k_cal: 1.2e-4,
            cycle_exp: 0.5,
        }
    }

    /// Capacity after `n_cycles` under standard cycling conditions.
    pub fn capacity_after_cycles(&self, n_cycles: f64) -> f64 {
        self.q0 * (1.0 - self.k_cal * n_cycles.powf(self.cycle_exp))
    }

    /// Capacity after time `t_hours` of calendar aging.
    pub fn capacity_calendar(&self, t_hours: f64) -> f64 {
        self.q0 * (1.0 - self.k_sei * t_hours.sqrt())
    }

    /// State of health (SOH) from current capacity.
    pub fn state_of_health(&self, q_current: f64) -> f64 {
        q_current / self.q0
    }

    /// Cycle life (number of cycles to 80% SOH).
    pub fn cycle_life_80pct(&self) -> f64 {
        let fade_target = 0.20;
        (fade_target / self.k_cal).powf(1.0 / self.cycle_exp)
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    fn approx_eq(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }

    // --- LiIonElectrode ---

    #[test]
    fn test_butler_volmer_zero_eta() {
        let e = LiIonElectrode::lco_cathode();
        let i = e.butler_volmer(0.0);
        assert!(i.abs() < 1e-10, "BV at zero overpotential should be ~0");
    }

    #[test]
    fn test_butler_volmer_positive_eta() {
        let e = LiIonElectrode::lco_cathode();
        let i = e.butler_volmer(0.05);
        assert!(
            i > 0.0,
            "Positive overpotential should give positive current"
        );
    }

    #[test]
    fn test_diffusion_time_constant() {
        let e = LiIonElectrode::graphite_anode();
        let tau = e.diffusion_time_constant();
        assert!(tau > 0.0);
        // tau = R^2 / (15 * D) = (5e-6)^2 / (15 * 3.9e-14) ≈ 0.04 s
        let expected = 5.0e-6_f64.powi(2) / (15.0 * 3.9e-14);
        assert!(approx_eq(tau, expected, 1e-5));
    }

    #[test]
    fn test_areal_capacity_positive() {
        let e = LiIonElectrode::lco_cathode();
        let q = e.areal_capacity();
        assert!(q > 0.0);
    }

    #[test]
    fn test_tafel_slope_positive() {
        let e = LiIonElectrode::lco_cathode();
        assert!(e.tafel_slope() > 0.0);
    }

    // --- NafionMembrane ---

    #[test]
    fn test_nafion_conductivity_positive() {
        let m = NafionMembrane::nafion117();
        let sigma = m.proton_conductivity();
        assert!(sigma > 0.0);
    }

    #[test]
    fn test_nafion_asr_positive() {
        let m = NafionMembrane::nafion117();
        let asr = m.area_specific_resistance();
        assert!(asr > 0.0);
    }

    #[test]
    fn test_nafion_sorption_isotherm() {
        let m = NafionMembrane::nafion117();
        let lambda_0 = m.sorption_isotherm(0.0);
        let lambda_1 = m.sorption_isotherm(1.0);
        assert!(lambda_0 >= 0.0);
        assert!(lambda_1 > lambda_0);
    }

    // --- GasDiffusionLayer ---

    #[test]
    fn test_gdl_effective_diffusivity() {
        let gdl = GasDiffusionLayer::toray_tgp_h120();
        let d_eff = gdl.effective_diffusivity(2.6e-5);
        assert!(d_eff < 2.6e-5 && d_eff > 0.0);
    }

    #[test]
    fn test_gdl_darcy_velocity() {
        let gdl = GasDiffusionLayer::toray_tgp_h120();
        let v = gdl.darcy_velocity(1e-5, 1000.0);
        assert!(v > 0.0);
    }

    // --- SolarCell ---

    #[test]
    fn test_solar_cell_voc_positive() {
        let cell = SolarCell::silicon();
        let voc = cell.voc();
        assert!(voc > 0.5 && voc < 1.0);
    }

    #[test]
    fn test_solar_cell_fill_factor_range() {
        let cell = SolarCell::silicon();
        let ff = cell.fill_factor();
        assert!(ff > 0.5 && ff < 0.99);
    }

    #[test]
    fn test_solar_cell_efficiency() {
        let cell = SolarCell::silicon();
        let eta = cell.efficiency(1000.0);
        assert!(eta > 0.0 && eta < 0.5);
    }

    #[test]
    fn test_sq_limit_silicon() {
        let eta = SolarCell::sq_limit(1.12);
        assert!(eta > 0.0);
    }

    #[test]
    fn test_solar_temp_coeff_negative() {
        let cell = SolarCell::silicon();
        assert!(cell.temp_coeff_voc() < 0.0);
    }

    // --- Thermoelectric ---

    #[test]
    fn test_zt_bi2te3() {
        let te = Thermoelectric::bi2te3();
        let zt = te.zt();
        // Bi2Te3 is state-of-the-art: ZT ~ 1
        assert!(zt > 0.5 && zt < 2.0);
    }

    #[test]
    fn test_power_factor_positive() {
        let te = Thermoelectric::bi2te3();
        assert!(te.power_factor() > 0.0);
    }

    #[test]
    fn test_peltier_coefficient() {
        let te = Thermoelectric::bi2te3();
        let pi = te.peltier_coefficient();
        assert!(approx_eq(pi, te.seebeck * te.temperature, 1e-20));
    }

    #[test]
    fn test_seebeck_voltage() {
        let te = Thermoelectric::bi2te3();
        let v = te.seebeck_voltage(10.0);
        assert!(approx_eq(v, te.seebeck * 10.0, 1e-20));
    }

    #[test]
    fn test_max_efficiency_positive() {
        let te = Thermoelectric::bi2te3();
        let eta = te.max_efficiency(350.0, 250.0);
        assert!(eta > 0.0 && eta < 1.0);
    }

    // --- Supercapacitor ---

    #[test]
    fn test_supercap_capacitance() {
        let sc = Supercapacitor::activated_carbon_500f();
        let c = sc.capacitance();
        assert!(approx_eq(c, sc.c_specific * sc.electrode_area, 1.0));
    }

    #[test]
    fn test_supercap_energy_stored() {
        let sc = Supercapacitor::activated_carbon_500f();
        let e = sc.energy_stored(1.0);
        assert!(approx_eq(e, 0.5 * sc.capacitance() * 1.0, 1.0));
    }

    #[test]
    fn test_supercap_peak_power() {
        let sc = Supercapacitor::activated_carbon_500f();
        let p = sc.peak_power();
        assert!(p > 0.0);
    }

    #[test]
    fn test_supercap_voltage_decay() {
        let sc = Supercapacitor::activated_carbon_500f();
        let v0 = 2.7;
        let v1 = sc.voltage_decay(v0, 3600.0);
        assert!(v1 < v0);
    }

    #[test]
    fn test_debye_length_dilute() {
        let ld = Supercapacitor::debye_length(0.001, 1.0, 298.15);
        // ~10 nm for 1 mM
        assert!(ld > 1e-9 && ld < 1e-6);
    }

    // --- MetalHydride ---

    #[test]
    fn test_lani5_eq_pressure() {
        let mh = MetalHydride::lani5();
        let peq = mh.equilibrium_pressure(298.15);
        assert!(peq > 0.0);
    }

    #[test]
    fn test_absorption_rate_zero_below_eq() {
        let mh = MetalHydride::lani5();
        let peq = mh.equilibrium_pressure(mh.temperature);
        let rate = mh.absorption_rate(0.0, peq * 0.5);
        assert!(rate == 0.0);
    }

    #[test]
    fn test_absorption_rate_positive_above_eq() {
        let mh = MetalHydride::lani5();
        let peq = mh.equilibrium_pressure(mh.temperature);
        let rate = mh.absorption_rate(0.0, peq * 2.0);
        assert!(rate > 0.0);
    }

    // --- NuclearMaterial ---

    #[test]
    fn test_swelling_zero_at_low_dpa() {
        let nm = NuclearMaterial::f82h();
        assert_eq!(nm.swelling(0.5), 0.0);
    }

    #[test]
    fn test_swelling_positive_above_threshold() {
        let nm = NuclearMaterial::f82h();
        let sw = nm.swelling(10.0);
        assert!(sw > 0.0);
    }

    #[test]
    fn test_irradiation_hardening_increases() {
        let nm = NuclearMaterial::f82h();
        let h1 = nm.irradiation_hardening(1.0);
        let h10 = nm.irradiation_hardening(10.0);
        assert!(h10 > h1);
    }

    #[test]
    fn test_irradiated_yield_stress() {
        let nm = NuclearMaterial::ss316l();
        let sy = nm.irradiated_yield_stress(5.0);
        assert!(sy > nm.yield_stress_0);
    }

    #[test]
    fn test_nrt_displacements_zero_below_threshold() {
        let nm = NuclearMaterial::f82h();
        assert_eq!(nm.nrt_displacements(nm.e_displacement - 1.0), 0.0);
    }

    // --- PiezoHarvester ---

    #[test]
    fn test_piezo_natural_frequency() {
        let ph = PiezoHarvester::pzt5a();
        let fn_ = ph.natural_frequency();
        assert!(fn_ > 0.0 && fn_ < 10000.0);
    }

    #[test]
    fn test_piezo_quality_factor() {
        let ph = PiezoHarvester::pzt5a();
        let qm = ph.quality_factor();
        assert!(approx_eq(qm, 25.0, 1e-10));
    }

    #[test]
    fn test_piezo_coupling_efficiency() {
        let ph = PiezoHarvester::pzt5a();
        let eta = ph.coupling_efficiency();
        assert!(eta > 0.0 && eta < 1.0);
    }

    // --- CapacityFadingModel ---

    #[test]
    fn test_capacity_fading_decreases() {
        let m = CapacityFadingModel::nmc18650();
        let q100 = m.capacity_after_cycles(100.0);
        let q500 = m.capacity_after_cycles(500.0);
        assert!(q100 > q500);
    }

    #[test]
    fn test_capacity_soh_at_start() {
        let m = CapacityFadingModel::nmc18650();
        let soh = m.state_of_health(m.q0);
        assert!(approx_eq(soh, 1.0, 1e-10));
    }

    #[test]
    fn test_cycle_life_positive() {
        let m = CapacityFadingModel::nmc18650();
        let life = m.cycle_life_80pct();
        assert!(life > 0.0);
    }

    // --- LiIonCell ---

    #[test]
    fn test_cell_ocv_positive() {
        let cell = LiIonCell::nmc18650();
        let ocv = cell.ocv(0.5, 0.5);
        assert!(ocv.abs() < 5.0);
    }

    #[test]
    fn test_pem_fuel_cell_voltage() {
        let fc = PemFuelCell::standard();
        let v = fc.cell_voltage(100.0);
        assert!((0.0..=1.3).contains(&v));
    }

    #[test]
    fn test_pem_power_density() {
        let fc = PemFuelCell::standard();
        let p = fc.power_density(1000.0);
        assert!(p > 0.0);
    }
}
