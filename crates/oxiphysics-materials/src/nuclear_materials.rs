// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Nuclear materials science — radiation damage and actinide materials.
//!
//! This module provides:
//! - [`RadiationDamage`]: NRT-DPA model, Kinchin-Pease, cascade, swelling
//! - [`FuelPellet`]: Thermal conductivity, fission gas release, gap width
//! - [`ActinideMaterial`]: Lattice parameter, bulk modulus, magnetic moment
//! - [`ZircaloyClad`]: Yield strength, creep rate, oxidation, hydrogen pickup
//! - [`lindhard_electronic_stopping`]: Electronic stopping power (Lindhard model)
//! - [`orowan_strengthening`]: Radiation-induced obstacle strengthening

#![allow(dead_code)]
#![allow(clippy::too_many_arguments)]

use std::f64::consts::PI;

// ---------------------------------------------------------------------------
// Enumerations
// ---------------------------------------------------------------------------

/// Nuclear fuel composition variants.
#[derive(Debug, Clone, PartialEq)]
pub enum NuclearFuel {
    /// Uranium dioxide (standard LWR fuel).
    UO2,
    /// Mixed oxide fuel (UO₂ + PuO₂).
    MixedOxide,
    /// Uranium nitride (advanced fuel).
    UN,
    /// Uranium carbide (fast reactor fuel).
    UC,
}

/// Actinide element.
#[derive(Debug, Clone, PartialEq)]
pub enum Actinide {
    /// Uranium (Z = 92).
    Uranium,
    /// Plutonium (Z = 94).
    Plutonium,
    /// Thorium (Z = 90).
    Thorium,
    /// Neptunium (Z = 93).
    Neptunium,
    /// Americium (Z = 95).
    Americium,
}

/// Zircaloy alloy grade.
#[derive(Debug, Clone, PartialEq)]
pub enum ZircaloyGrade {
    /// Zircaloy-2 (Sn-Fe-Cr-Ni alloy, BWR cladding).
    Zircaloy2,
    /// Zircaloy-4 (Sn-Fe-Cr alloy, PWR cladding, low Ni).
    Zircaloy4,
    /// Zirlo (advanced Zr-Nb-Sn alloy).
    Zirlo,
    /// M5 (Zr-1 %Nb binary alloy).
    M5,
}

// ---------------------------------------------------------------------------
// RadiationDamage
// ---------------------------------------------------------------------------

/// Radiation damage model for crystalline metals.
///
/// Implements the NRT (Norgett-Robinson-Torrens) displacement model,
/// Kinchin-Pease approximation, cascade statistics, recombination,
/// and void swelling correlations relevant to structural materials
/// in fission and fusion reactors.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct RadiationDamage {
    /// Threshold displacement energy Ed \[eV\] for the host lattice atom.
    pub displacement_energy: f64,
    /// Primary knock-on atom (PKA) kinetic energy \[eV\].
    pub pka_energy: f64,
    /// Atomic number of the target atom.
    pub atomic_number: u32,
    /// Atomic mass of the target atom \[amu\].
    pub atomic_mass: f64,
}

impl RadiationDamage {
    /// Create a new radiation damage model.
    ///
    /// # Arguments
    /// * `displacement_energy` – Ed in eV (typical: 25–40 eV for metals).
    /// * `pka_energy` – PKA recoil energy in eV.
    /// * `atomic_number` – Z of host lattice.
    /// * `atomic_mass` – M of host lattice in amu.
    pub fn new(
        displacement_energy: f64,
        pka_energy: f64,
        atomic_number: u32,
        atomic_mass: f64,
    ) -> Self {
        Self {
            displacement_energy,
            pka_energy,
            atomic_number,
            atomic_mass,
        }
    }

    /// NRT displacements per atom (dpa) for a given neutron fluence.
    ///
    /// dpa = 0.8 · ν_NRT · fluence / (2 · Ed)
    /// where ν_NRT is the number of displacements from one PKA event.
    ///
    /// # Arguments
    /// * `fluence` – Neutron fluence \[n/m²\].
    pub fn nrt_dpa(&self, fluence: f64) -> f64 {
        let nu = self.kinchin_pease(self.pka_energy) as f64;
        0.8 * nu * fluence / (2.0 * self.displacement_energy)
    }

    /// Kinchin-Pease number of displacements from a single PKA.
    ///
    /// Returns the integer number of Frenkel pairs produced:
    /// - 0 if energy < Ed
    /// - 1 if Ed ≤ energy < 2Ed
    /// - energy/(2Ed) otherwise (integer division)
    ///
    /// # Arguments
    /// * `energy` – PKA kinetic energy in eV.
    pub fn kinchin_pease(&self, energy: f64) -> usize {
        let ed = self.displacement_energy;
        if energy < ed {
            0
        } else if energy < 2.0 * ed {
            1
        } else {
            (energy / (2.0 * ed)) as usize
        }
    }

    /// Estimate the number of atoms in a displacement cascade.
    ///
    /// Uses an empirical power law: N_cascade ≈ (E_PKA / Ed)^0.75.
    ///
    /// # Arguments
    /// * `energy` – PKA energy in eV.
    pub fn cascade_size(&self, energy: f64) -> usize {
        if energy <= self.displacement_energy {
            return 1;
        }
        let ratio = energy / self.displacement_energy;
        ratio.powf(0.75).round() as usize
    }

    /// Recombination fraction of Frenkel pairs in a cascade.
    ///
    /// Based on Brinkman's model: f_rec = 1 − 1/(1 + α·ν) where α = 0.15.
    pub fn recombination_fraction(&self) -> f64 {
        let nu = self.kinchin_pease(self.pka_energy) as f64;
        if nu < 1.0 {
            return 0.0;
        }
        let alpha = 0.15_f64;
        1.0 - 1.0 / (1.0 + alpha * nu)
    }

    /// Void swelling as a function of cumulative damage in dpa.
    ///
    /// Empirical model: ΔV/V ≈ A · (dpa − dpa_inc)^n for dpa > dpa_inc,
    /// with A = 1.5e-3 %/dpa^n, dpa_inc = 5 dpa (incubation dose), n = 1.5.
    ///
    /// Returns swelling in percent.
    ///
    /// # Arguments
    /// * `dpa` – Cumulative dose in displacements per atom.
    pub fn swelling(&self, dpa: f64) -> f64 {
        let dpa_inc = 5.0_f64;
        let a = 1.5e-3_f64;
        let n = 1.5_f64;
        if dpa <= dpa_inc {
            0.0
        } else {
            a * (dpa - dpa_inc).powf(n)
        }
    }
}

// ---------------------------------------------------------------------------
// FuelPellet
// ---------------------------------------------------------------------------

/// Nuclear fuel pellet model.
///
/// Computes thermal-mechanical properties relevant to fuel performance:
/// thermal conductivity, fission gas release, pellet-cladding gap, and
/// centerline temperature under linear heat rate.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct FuelPellet {
    /// Fuel composition.
    pub material: NuclearFuel,
    /// Burnup \[MWd/tHM\].
    pub burnup: f64,
    /// Mean fuel pellet temperature \[K\].
    pub temperature: f64,
    /// Pellet radius \[m\].
    pub radius: f64,
}

impl FuelPellet {
    /// Create a new fuel pellet.
    ///
    /// # Arguments
    /// * `material` – Fuel type.
    /// * `burnup` – Accumulated burnup in MWd/tHM.
    /// * `temperature` – Average pellet temperature in K.
    /// * `radius` – Pellet outer radius in m.
    pub fn new(material: NuclearFuel, burnup: f64, temperature: f64, radius: f64) -> Self {
        Self {
            material,
            burnup,
            temperature,
            radius,
        }
    }

    /// Thermal conductivity of the fuel \[W/(m·K)\].
    ///
    /// Uses the Lucuta correlation for UO₂:
    /// k = 1/(A + B·T) + C·T³·exp(−D/T) corrected for burnup degradation.
    /// MOX uses a 20 % reduction factor.  UN and UC use empirical values.
    pub fn thermal_conductivity(&self) -> f64 {
        let t = self.temperature;
        let base = match self.material {
            NuclearFuel::UO2 => {
                // Lucuta model: 1/(0.0452 + 2.46e-4·T) + 5.47e9·T^-2.5·exp(-16350/T)
                let a = 0.0452_f64;
                let b = 2.46e-4_f64;
                let lattice_term = 1.0 / (a + b * t);
                let electron_term = 5.47e9 * t.powf(-2.5) * (-16350.0 / t).exp();
                lattice_term + electron_term
            }
            NuclearFuel::MixedOxide => {
                let a = 0.0452_f64;
                let b = 2.46e-4_f64;
                0.80 / (a + b * t)
            }
            NuclearFuel::UN => 1.864e-2 * t.powf(0.5) + 12.0,
            NuclearFuel::UC => 20.0 * (t / 1000.0).powf(-0.2),
        };
        // Burnup degradation factor (porosity + fission products)
        let bu_factor = 1.0 - 0.02 * (self.burnup / 10_000.0).min(1.0);
        base * bu_factor
    }

    /// Fractional fission gas release (Xe + Kr) from Booth diffusion model.
    ///
    /// f_FGR = 1 − (6/π²) Σ_{n=1}^{10} exp(−n²·π²·D_eff·t/a²) / n²
    /// Simplified to: f_FGR ≈ 1 − exp(−β·Bu·T²) with empirical β.
    pub fn fission_gas_release(&self) -> f64 {
        let beta = 3.0e-12_f64;
        let raw = 1.0 - (-beta * self.burnup * self.temperature * self.temperature).exp();
        raw.clamp(0.0, 1.0)
    }

    /// Pellet-cladding gap width \[m\].
    ///
    /// Accounts for fuel swelling due to fission products and thermal expansion.
    /// Gap = initial_gap − Δr_fuel where Δr_fuel = r₀·(α·ΔT + β_sw·burnup).
    pub fn pellet_cladding_gap(&self) -> f64 {
        let initial_gap = 80e-6_f64; // 80 μm
        let alpha_th = 10e-6_f64; // thermal expansion coefficient K⁻¹
        let beta_sw = 0.5e-9_f64; // swelling coefficient per MWd/tHM
        let delta_t = (self.temperature - 300.0).max(0.0);
        let delta_r = self.radius * (alpha_th * delta_t + beta_sw * self.burnup);
        (initial_gap - delta_r).max(0.0)
    }

    /// Centerline temperature of the fuel pellet under a given linear heat rate.
    ///
    /// T_cl = T_surface + q′/(4π·k_avg) where k_avg is evaluated at T_surface.
    ///
    /// # Arguments
    /// * `q_linear` – Linear heat rate \[W/m\].
    pub fn centerline_temperature(&self, q_linear: f64) -> f64 {
        let k = self.thermal_conductivity();
        let t_surface = self.temperature;
        t_surface + q_linear / (4.0 * PI * k)
    }
}

// ---------------------------------------------------------------------------
// ActinideMaterial
// ---------------------------------------------------------------------------

/// Actinide metal or oxide material descriptor.
///
/// Provides crystallographic and electronic properties derived from
/// experimental data and relativistic DFT calculations.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct ActinideMaterial {
    /// Actinide element.
    pub element: Actinide,
    /// Formal oxidation state (e.g. +4 for UO₂).
    pub oxidation_state: i32,
    /// Temperature \[K\].
    pub temperature: f64,
}

impl ActinideMaterial {
    /// Create a new actinide material descriptor.
    ///
    /// # Arguments
    /// * `element` – Actinide element.
    /// * `oxidation_state` – Formal oxidation state.
    /// * `temperature` – Temperature in K.
    pub fn new(element: Actinide, oxidation_state: i32, temperature: f64) -> Self {
        Self {
            element,
            oxidation_state,
            temperature,
        }
    }

    /// Equilibrium lattice parameter \[Å\].
    ///
    /// Returns experimental 0 K values corrected by linear thermal expansion.
    pub fn lattice_parameter(&self) -> f64 {
        let a0 = match self.element {
            Actinide::Uranium => 2.854,   // orthorhombic α-U a-axis
            Actinide::Plutonium => 3.159, // δ-Pu fcc
            Actinide::Thorium => 5.084,   // fcc α-Th
            Actinide::Neptunium => 6.663, // α-Np orthorhombic
            Actinide::Americium => 3.468, // α-Am dhcp
        };
        let alpha = match self.element {
            Actinide::Uranium => 14e-6,
            Actinide::Plutonium => 60e-6, // anomalously large
            Actinide::Thorium => 11e-6,
            Actinide::Neptunium => 15e-6,
            Actinide::Americium => 12e-6,
        };
        a0 * (1.0 + alpha * (self.temperature - 293.0))
    }

    /// Bulk modulus \[GPa\].
    ///
    /// Returns experimental values at room temperature
    /// with a small linear temperature correction.
    pub fn bulk_modulus(&self) -> f64 {
        let b0 = match self.element {
            Actinide::Uranium => 115.0,
            Actinide::Plutonium => 35.0,
            Actinide::Thorium => 58.0,
            Actinide::Neptunium => 67.0,
            Actinide::Americium => 30.0,
        };
        let db_dt = -0.02_f64; // GPa/K
        (b0 + db_dt * (self.temperature - 300.0)).max(0.0)
    }

    /// Ordered magnetic moment \[μ_B / atom\].
    ///
    /// Returns the experimental low-temperature moment; returns 0.0 for
    /// non-magnetic phases (U, Pu in metallic form at RT).
    pub fn magnetic_moment(&self) -> f64 {
        match self.element {
            Actinide::Uranium => 0.0,   // itinerant, quenched
            Actinide::Plutonium => 0.0, // paramagnetic
            Actinide::Thorium => 0.0,   // non-magnetic
            Actinide::Neptunium => 0.3, // weak moment in α-Np
            Actinide::Americium => 0.0, // non-magnetic metal
        }
    }

    /// Ground-state electronic configuration string.
    ///
    /// Returns the abbreviated electron configuration beyond \[Rn\].
    pub fn electronic_configuration(&self) -> &str {
        match self.element {
            Actinide::Thorium => "[Rn] 6d2 7s2",
            Actinide::Uranium => "[Rn] 5f3 6d1 7s2",
            Actinide::Neptunium => "[Rn] 5f4 6d1 7s2",
            Actinide::Plutonium => "[Rn] 5f6 7s2",
            Actinide::Americium => "[Rn] 5f7 7s2",
        }
    }
}

// ---------------------------------------------------------------------------
// ZircaloyClad
// ---------------------------------------------------------------------------

/// Zircaloy cladding model for LWR fuel rods.
///
/// Computes mechanical and corrosion properties relevant to in-reactor
/// and out-of-reactor cladding performance: yield strength, creep rate,
/// oxidation kinetics, and hydrogen pickup.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct ZircaloyClad {
    /// Zircaloy alloy grade.
    pub composition: ZircaloyGrade,
    /// Cladding temperature \[K\].
    pub temperature: f64,
    /// Fast neutron fluence \[n/m²\] (E > 1 MeV).
    pub fluence: f64,
}

impl ZircaloyClad {
    /// Create a new Zircaloy cladding instance.
    ///
    /// # Arguments
    /// * `composition` – Alloy grade.
    /// * `temperature` – Temperature in K.
    /// * `fluence` – Accumulated fast fluence in n/m².
    pub fn new(composition: ZircaloyGrade, temperature: f64, fluence: f64) -> Self {
        Self {
            composition,
            temperature,
            fluence,
        }
    }

    /// 0.2 % proof (yield) strength \[MPa\].
    ///
    /// Based on the Khan-Hall correlation including radiation hardening:
    /// σ_y = (σ_y0 + Δσ_irr) · exp(−T/T_ref)^0.3
    pub fn yield_strength(&self) -> f64 {
        let sigma_y0 = match self.composition {
            ZircaloyGrade::Zircaloy2 => 380.0,
            ZircaloyGrade::Zircaloy4 => 400.0,
            ZircaloyGrade::Zirlo => 480.0,
            ZircaloyGrade::M5 => 420.0,
        };
        // Radiation hardening: Δσ ≈ A·φt^0.5
        let a_irr = 2.5e-10_f64;
        let delta_sigma = a_irr * self.fluence.sqrt();
        let t_ref = 900.0_f64;
        (sigma_y0 + delta_sigma) * (self.temperature / t_ref).powf(-0.3)
    }

    /// Steady-state in-reactor creep rate \[s⁻¹\].
    ///
    /// εdot = A · σ^n · φ · exp(−Q/(R·T))
    /// where σ = 100 MPa (reference), Q = 230 kJ/mol (thermal component).
    pub fn creep_rate(&self) -> f64 {
        let a = match self.composition {
            ZircaloyGrade::Zircaloy2 | ZircaloyGrade::Zircaloy4 => 1.0e-25,
            ZircaloyGrade::Zirlo | ZircaloyGrade::M5 => 5.0e-26,
        };
        let sigma_ref = 100.0_f64; // MPa
        let n = 3.0_f64;
        let q = 230_000.0_f64; // J/mol
        let r = 8.314_f64;
        a * sigma_ref.powf(n) * self.fluence * (-q / (r * self.temperature)).exp()
    }

    /// Zirconium oxide (ZrO₂) layer growth rate \[nm/s\].
    ///
    /// Pre-transition cubic kinetics: dr/dt = A · exp(−Q/(R·T))
    /// where A and Q depend on alloy grade.
    pub fn oxidation_rate(&self) -> f64 {
        let (a, q) = match self.composition {
            ZircaloyGrade::Zircaloy2 => (2.3e4_f64, 96_000.0_f64),
            ZircaloyGrade::Zircaloy4 => (2.0e4_f64, 96_000.0_f64),
            ZircaloyGrade::Zirlo => (1.2e4_f64, 99_000.0_f64),
            ZircaloyGrade::M5 => (1.0e4_f64, 102_000.0_f64),
        };
        let r = 8.314_f64;
        a * (-q / (r * self.temperature)).exp()
    }

    /// Cumulative hydrogen pickup in the cladding wall \[wppm\].
    ///
    /// H pickup = f_H · (δ_oxide(t) / δ_0) where f_H is the pickup fraction
    /// and δ_oxide = k_c · t^(1/3) (cubic oxidation law).
    ///
    /// # Arguments
    /// * `time` – Irradiation / corrosion time in seconds.
    pub fn hydrogen_pickup(&self, time: f64) -> f64 {
        let f_h = match self.composition {
            ZircaloyGrade::Zircaloy2 => 0.15, // BWR, higher pickup
            ZircaloyGrade::Zircaloy4 => 0.10,
            ZircaloyGrade::Zirlo => 0.08,
            ZircaloyGrade::M5 => 0.05, // lowest pickup
        };
        // Cubic oxidation: thickness proportional to t^(1/3)
        let kc = self.oxidation_rate() * 1e9; // nm/s to nm^(1/3)·s^(-1/3) scale
        let oxide_thickness = kc * time.powf(1.0 / 3.0); // nm
        // H pickup in wppm: empirical proportionality
        f_h * oxide_thickness * 0.4
    }
}

// ---------------------------------------------------------------------------
// Standalone functions
// ---------------------------------------------------------------------------

/// Lindhard electronic stopping cross section \[eV·cm²/atom\].
///
/// Calculates the low-energy electronic stopping of an ion (Z₁, A₁) in a
/// target (Z₂, A₂) using the Lindhard-Scharff-Schiøtt formula:
///
/// Se = k_L · ε^{1/2}
///
/// where ε is the reduced energy and k_L is the Lindhard factor.
///
/// # Arguments
/// * `energy` – Ion kinetic energy \[eV\].
/// * `z1` – Atomic number of the projectile.
/// * `z2` – Atomic number of the target atom.
/// * `a1` – Atomic mass of the projectile \[amu\].
/// * `a2` – Atomic mass of the target atom \[amu\].
pub fn lindhard_electronic_stopping(energy: f64, z1: f64, z2: f64, a1: f64, a2: f64) -> f64 {
    // Screening length (Thomas-Fermi, Ziegler form) [Å]
    let a_tf = 0.529 * 0.8853 / (z1.powf(2.0 / 3.0) + z2.powf(2.0 / 3.0)).sqrt();
    // Reduced mass
    let m12 = a1 * a2 / (a1 + a2);
    // Reduced energy ε (dimensionless)
    let epsilon = energy * a2 * a_tf / ((a1 + a2) * z1 * z2 * 14.4);
    // Lindhard factor k_L
    let k_l = 0.0793 * z1.powf(2.0 / 3.0) * z2.powf(1.0 / 2.0) * (a1 + a2).powf(3.0 / 2.0)
        / ((z1.powf(2.0 / 3.0) + z2.powf(2.0 / 3.0)).powf(3.0 / 4.0)
            * a1.powf(3.0 / 2.0)
            * a2.powf(1.0 / 2.0));
    let _ = m12; // mass used in extended model
    k_l * epsilon.sqrt() * a_tf
}

/// Orowan strengthening increment due to radiation-induced obstacles \[MPa\].
///
/// Computes the critical resolved shear stress increment from a dispersed
/// field of impenetrable obstacles (loops, voids, precipitates) using the
/// Orowan-Taylor relation:
///
/// Δτ = M · μ · b / (2π · L · √(1−ν))
///
/// Simplified form without Poisson's ratio:
///
/// Δτ = 0.84 · M · μ · b / L
///
/// where L = √(spacing²−(2r)²) ≈ spacing for small radii.
///
/// # Arguments
/// * `spacing` – Mean obstacle spacing \[m\] (= 1/√(N·r) for N obstacles/vol).
/// * `radius` – Obstacle radius \[m\].
/// * `shear_mod` – Shear modulus of the matrix \[Pa\].
/// * `burgers` – Burgers vector magnitude \[m\].
pub fn orowan_strengthening(spacing: f64, radius: f64, shear_mod: f64, burgers: f64) -> f64 {
    let taylor_factor = 3.06_f64; // Taylor factor M for fcc/bcc polycrystal
    // Effective obstacle spacing: bypass requires bowing between obstacles
    let l_eff = (spacing * spacing - 4.0 * radius * radius)
        .max(spacing * 0.01)
        .sqrt();
    let delta_tau = 0.84 * taylor_factor * shear_mod * burgers / l_eff;
    delta_tau * 1e-6 // Pa → MPa
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // --- RadiationDamage ---

    #[test]
    fn test_kinchin_pease_below_threshold() {
        let rd = RadiationDamage::new(40.0, 100.0, 26, 55.85);
        assert_eq!(rd.kinchin_pease(20.0), 0);
    }

    #[test]
    fn test_kinchin_pease_at_threshold() {
        let rd = RadiationDamage::new(40.0, 100.0, 26, 55.85);
        assert_eq!(rd.kinchin_pease(40.0), 1);
    }

    #[test]
    fn test_kinchin_pease_above_threshold() {
        let rd = RadiationDamage::new(40.0, 100.0, 26, 55.85);
        // 400 eV / (2*40) = 5
        assert_eq!(rd.kinchin_pease(400.0), 5);
    }

    #[test]
    fn test_nrt_dpa_positive() {
        let rd = RadiationDamage::new(40.0, 200.0, 26, 55.85);
        let dpa = rd.nrt_dpa(1e22);
        assert!(dpa > 0.0);
    }

    #[test]
    fn test_nrt_dpa_scales_with_fluence() {
        let rd = RadiationDamage::new(40.0, 200.0, 26, 55.85);
        let dpa1 = rd.nrt_dpa(1e22);
        let dpa2 = rd.nrt_dpa(2e22);
        assert!((dpa2 / dpa1 - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_cascade_size_below_threshold() {
        let rd = RadiationDamage::new(40.0, 100.0, 26, 55.85);
        assert_eq!(rd.cascade_size(30.0), 1);
    }

    #[test]
    fn test_cascade_size_increases_with_energy() {
        let rd = RadiationDamage::new(25.0, 1000.0, 26, 55.85);
        let s1 = rd.cascade_size(100.0);
        let s2 = rd.cascade_size(1000.0);
        assert!(s2 > s1);
    }

    #[test]
    fn test_recombination_fraction_range() {
        let rd = RadiationDamage::new(40.0, 400.0, 26, 55.85);
        let f = rd.recombination_fraction();
        assert!((0.0..=1.0).contains(&f));
    }

    #[test]
    fn test_swelling_below_incubation() {
        let rd = RadiationDamage::new(40.0, 100.0, 26, 55.85);
        assert_eq!(rd.swelling(3.0), 0.0);
    }

    #[test]
    fn test_swelling_above_incubation() {
        let rd = RadiationDamage::new(40.0, 100.0, 26, 55.85);
        assert!(rd.swelling(20.0) > 0.0);
    }

    #[test]
    fn test_swelling_monotone() {
        let rd = RadiationDamage::new(40.0, 100.0, 26, 55.85);
        assert!(rd.swelling(30.0) > rd.swelling(20.0));
    }

    // --- FuelPellet ---

    #[test]
    fn test_uo2_thermal_conductivity_positive() {
        let p = FuelPellet::new(NuclearFuel::UO2, 10_000.0, 800.0, 4.1e-3);
        assert!(p.thermal_conductivity() > 0.0);
    }

    #[test]
    fn test_mox_lower_conductivity_than_uo2() {
        let uo2 = FuelPellet::new(NuclearFuel::UO2, 0.0, 800.0, 4.1e-3);
        let mox = FuelPellet::new(NuclearFuel::MixedOxide, 0.0, 800.0, 4.1e-3);
        assert!(mox.thermal_conductivity() < uo2.thermal_conductivity());
    }

    #[test]
    fn test_un_thermal_conductivity_reasonable() {
        let p = FuelPellet::new(NuclearFuel::UN, 5_000.0, 900.0, 4.1e-3);
        let k = p.thermal_conductivity();
        assert!(k > 1.0 && k < 50.0);
    }

    #[test]
    fn test_fission_gas_release_at_low_temp_burnup() {
        let p = FuelPellet::new(NuclearFuel::UO2, 100.0, 500.0, 4.1e-3);
        let fgr = p.fission_gas_release();
        assert!((0.0..=1.0).contains(&fgr));
    }

    #[test]
    fn test_fission_gas_release_increases_with_burnup() {
        let p1 = FuelPellet::new(NuclearFuel::UO2, 5_000.0, 1000.0, 4.1e-3);
        let p2 = FuelPellet::new(NuclearFuel::UO2, 50_000.0, 1000.0, 4.1e-3);
        assert!(p2.fission_gas_release() >= p1.fission_gas_release());
    }

    #[test]
    fn test_pellet_cladding_gap_nonnegative() {
        let p = FuelPellet::new(NuclearFuel::UO2, 60_000.0, 1200.0, 4.1e-3);
        assert!(p.pellet_cladding_gap() >= 0.0);
    }

    #[test]
    fn test_centerline_temperature_above_surface() {
        let p = FuelPellet::new(NuclearFuel::UO2, 10_000.0, 700.0, 4.1e-3);
        let t_cl = p.centerline_temperature(20_000.0);
        assert!(t_cl > p.temperature);
    }

    // --- ActinideMaterial ---

    #[test]
    fn test_uranium_lattice_parameter_positive() {
        let am = ActinideMaterial::new(Actinide::Uranium, 0, 293.0);
        assert!(am.lattice_parameter() > 0.0);
    }

    #[test]
    fn test_plutonium_lattice_expands_with_temperature() {
        let am1 = ActinideMaterial::new(Actinide::Plutonium, 0, 300.0);
        let am2 = ActinideMaterial::new(Actinide::Plutonium, 0, 600.0);
        assert!(am2.lattice_parameter() > am1.lattice_parameter());
    }

    #[test]
    fn test_thorium_bulk_modulus_decreases_with_temperature() {
        let am1 = ActinideMaterial::new(Actinide::Thorium, 0, 300.0);
        let am2 = ActinideMaterial::new(Actinide::Thorium, 0, 800.0);
        assert!(am2.bulk_modulus() < am1.bulk_modulus());
    }

    #[test]
    fn test_neptunium_magnetic_moment_nonzero() {
        let am = ActinideMaterial::new(Actinide::Neptunium, 0, 293.0);
        assert!(am.magnetic_moment() > 0.0);
    }

    #[test]
    fn test_uranium_magnetic_moment_zero() {
        let am = ActinideMaterial::new(Actinide::Uranium, 4, 293.0);
        assert_eq!(am.magnetic_moment(), 0.0);
    }

    #[test]
    fn test_electronic_configuration_uranium() {
        let am = ActinideMaterial::new(Actinide::Uranium, 4, 293.0);
        assert!(am.electronic_configuration().contains("5f3"));
    }

    #[test]
    fn test_electronic_configuration_thorium() {
        let am = ActinideMaterial::new(Actinide::Thorium, 4, 293.0);
        assert!(am.electronic_configuration().contains("6d2"));
    }

    // --- ZircaloyClad ---

    #[test]
    fn test_yield_strength_positive() {
        let zr = ZircaloyClad::new(ZircaloyGrade::Zircaloy4, 600.0, 1e25);
        assert!(zr.yield_strength() > 0.0);
    }

    #[test]
    fn test_yield_strength_increases_with_fluence() {
        let zr1 = ZircaloyClad::new(ZircaloyGrade::Zircaloy4, 600.0, 0.0);
        let zr2 = ZircaloyClad::new(ZircaloyGrade::Zircaloy4, 600.0, 1e26);
        assert!(zr2.yield_strength() > zr1.yield_strength());
    }

    #[test]
    fn test_creep_rate_positive() {
        let zr = ZircaloyClad::new(ZircaloyGrade::Zircaloy2, 620.0, 1e25);
        assert!(zr.creep_rate() > 0.0);
    }

    #[test]
    fn test_oxidation_rate_positive() {
        let zr = ZircaloyClad::new(ZircaloyGrade::M5, 620.0, 1e25);
        assert!(zr.oxidation_rate() > 0.0);
    }

    #[test]
    fn test_m5_lower_oxidation_than_zircaloy4() {
        let zr4 = ZircaloyClad::new(ZircaloyGrade::Zircaloy4, 620.0, 1e25);
        let m5 = ZircaloyClad::new(ZircaloyGrade::M5, 620.0, 1e25);
        assert!(m5.oxidation_rate() < zr4.oxidation_rate());
    }

    #[test]
    fn test_hydrogen_pickup_increases_with_time() {
        let zr = ZircaloyClad::new(ZircaloyGrade::Zircaloy4, 620.0, 1e25);
        let h1 = zr.hydrogen_pickup(1e6);
        let h2 = zr.hydrogen_pickup(1e8);
        assert!(h2 > h1);
    }

    #[test]
    fn test_m5_lower_hydrogen_pickup_than_zircaloy2() {
        let zr2 = ZircaloyClad::new(ZircaloyGrade::Zircaloy2, 620.0, 1e25);
        let m5 = ZircaloyClad::new(ZircaloyGrade::M5, 620.0, 1e25);
        assert!(m5.hydrogen_pickup(1e7) < zr2.hydrogen_pickup(1e7));
    }

    // --- Standalone functions ---

    #[test]
    fn test_lindhard_stopping_positive() {
        // Fission fragment (Kr, Z=36, A=84) in UO2 target (U, Z=92, A=238)
        let s = lindhard_electronic_stopping(1e6, 36.0, 92.0, 84.0, 238.0);
        assert!(s > 0.0);
    }

    #[test]
    fn test_lindhard_stopping_increases_with_energy() {
        // Higher PKA energy → larger ε → larger Se (Se ∝ ε^0.5)
        let s1 = lindhard_electronic_stopping(1e5, 18.0, 26.0, 40.0, 55.85);
        let s2 = lindhard_electronic_stopping(1e6, 18.0, 26.0, 40.0, 55.85);
        assert!(s2 > s1);
    }

    #[test]
    fn test_orowan_strengthening_positive() {
        // Typical loop: spacing = 50 nm, radius = 2 nm, μ = 80 GPa, b = 0.286 nm
        let delta = orowan_strengthening(50e-9, 2e-9, 80e9, 0.286e-9);
        assert!(delta > 0.0);
    }

    #[test]
    fn test_orowan_strengthening_decreases_with_spacing() {
        let d1 = orowan_strengthening(20e-9, 1e-9, 80e9, 0.286e-9);
        let d2 = orowan_strengthening(100e-9, 1e-9, 80e9, 0.286e-9);
        assert!(d1 > d2);
    }

    #[test]
    fn test_uc_thermal_conductivity() {
        let p = FuelPellet::new(NuclearFuel::UC, 0.0, 1000.0, 4.1e-3);
        assert!(p.thermal_conductivity() > 0.0);
    }

    #[test]
    fn test_americium_electronic_config() {
        let am = ActinideMaterial::new(Actinide::Americium, 3, 300.0);
        assert!(am.electronic_configuration().contains("5f7"));
    }

    #[test]
    fn test_zirlo_yield_strength() {
        let zr = ZircaloyClad::new(ZircaloyGrade::Zirlo, 600.0, 1e25);
        assert!(zr.yield_strength() > 0.0);
    }

    #[test]
    fn test_plutonium_magnetic_moment_zero() {
        let am = ActinideMaterial::new(Actinide::Plutonium, 0, 300.0);
        assert_eq!(am.magnetic_moment(), 0.0);
    }
}
