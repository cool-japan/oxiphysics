//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use std::f64::consts::PI;

/// Energy-absorbing crash structure model (foam/tube progressive crushing).
pub struct CrashAbsorber {
    /// Designation.
    pub designation: String,
    /// Crush stress σ_c (MPa).
    pub crush_stress: f64,
    /// Trigger mechanism reduction factor (0.6–0.9).
    pub trigger_factor: f64,
    /// Cross-sectional area A (mm²).
    pub area: f64,
    /// Total crushable length L (mm).
    pub length: f64,
    /// Densification strain ε_d (typically 0.6–0.8).
    pub densification_strain: f64,
    /// Density of absorber material (kg/m³).
    pub density: f64,
}
impl CrashAbsorber {
    /// Create an aluminum foam crash absorber.
    pub fn aluminum_foam(crush_stress: f64, area: f64, length: f64) -> Self {
        CrashAbsorber {
            designation: "Aluminum foam".to_string(),
            crush_stress,
            trigger_factor: 0.8,
            area,
            length,
            densification_strain: 0.7,
            density: 400.0,
        }
    }
    /// Initial peak force with trigger mechanism (N).
    pub fn peak_force(&self) -> f64 {
        self.crush_stress * self.area * self.trigger_factor
    }
    /// Mean crush force (N).
    pub fn mean_crush_force(&self) -> f64 {
        self.crush_stress * self.area
    }
    /// Total energy absorbed (J) before densification.
    pub fn total_energy_absorbed(&self) -> f64 {
        let crush_length = self.length * self.densification_strain;
        self.mean_crush_force() * crush_length
    }
    /// Specific energy absorption SEA (J/kg).
    pub fn specific_energy_absorption(&self) -> f64 {
        let mass = self.area * self.length * 1e-9 * self.density;
        self.total_energy_absorbed() / mass.max(1e-30)
    }
    /// Crush stroke to absorb energy `e` (J).
    pub fn stroke_for_energy(&self, e: f64) -> f64 {
        e / self.mean_crush_force()
    }
}
/// Nickel-base superalloy material properties (e.g., IN718, Rene 104).
pub struct SuperalloyMaterial {
    /// Alloy designation.
    pub designation: String,
    /// 0.2% proof strength at 650°C (MPa).
    pub sigma_02_650: f64,
    /// Ultimate tensile strength at 650°C (MPa).
    pub uts_650: f64,
    /// Creep rupture life at 150 MPa, 750°C (hours).
    pub creep_rupture_life: f64,
    /// Larson-Miller constant C.
    pub larson_miller_c: f64,
    /// Volume fraction of γ' precipitates.
    pub gamma_prime_fraction: f64,
    /// Oxidation rate constant kp at 900°C (mg²/cm⁴/h).
    pub oxidation_kp: f64,
    /// Density (kg/m³).
    pub density: f64,
}
impl SuperalloyMaterial {
    /// Inconel 718 standard properties.
    pub fn in718() -> Self {
        SuperalloyMaterial {
            designation: "IN718".to_string(),
            sigma_02_650: 1000.0,
            uts_650: 1250.0,
            creep_rupture_life: 1000.0,
            larson_miller_c: 20.0,
            gamma_prime_fraction: 0.18,
            oxidation_kp: 1e-4,
            density: 8190.0,
        }
    }
    /// Larson-Miller parameter P = T * (C + log10(tr)) where T in Kelvin.
    pub fn larson_miller_parameter(&self, temp_k: f64, rupture_hours: f64) -> f64 {
        temp_k * (self.larson_miller_c + rupture_hours.log10())
    }
    /// Estimated creep rupture life at temperature and stress (simplified).
    pub fn creep_rupture_estimate(&self, temp_k: f64, stress_mpa: f64) -> f64 {
        let p_ref = self.larson_miller_parameter(923.15, self.creep_rupture_life);
        let stress_factor = (self.sigma_02_650 / stress_mpa).powi(2);
        let tr = 10.0_f64.powf(p_ref / temp_k - self.larson_miller_c);
        tr * stress_factor
    }
    /// Oxidation mass gain after time `t` (hours): Δm = sqrt(kp * t).
    pub fn oxidation_mass_gain(&self, t_hours: f64) -> f64 {
        (self.oxidation_kp * t_hours).sqrt()
    }
}
/// Ceramic matrix composite (SiC/SiC or C/SiC) material properties.
pub struct CmcMaterial {
    /// Material designation (e.g., "SiC/SiC").
    pub designation: String,
    /// Matrix-cracking stress σ_mc (MPa).
    pub matrix_cracking_stress: f64,
    /// Ultimate tensile strength UTS (MPa).
    pub uts: f64,
    /// Elastic modulus E (GPa).
    pub modulus: f64,
    /// Thermal conductivity k (W/m·K).
    pub thermal_conductivity: f64,
    /// Coefficient of thermal expansion CTE (1/K).
    pub cte: f64,
    /// Maximum use temperature (°C).
    pub max_temperature: f64,
    /// Fracture toughness KIc (MPa·√m).
    pub k_ic: f64,
}
impl CmcMaterial {
    /// SiC/SiC composite properties at room temperature.
    pub fn sic_sic() -> Self {
        CmcMaterial {
            designation: "SiC/SiC".to_string(),
            matrix_cracking_stress: 120.0,
            uts: 400.0,
            modulus: 230.0,
            thermal_conductivity: 15.0,
            cte: 4.0e-6,
            max_temperature: 1300.0,
            k_ic: 25.0,
        }
    }
    /// C/SiC composite properties at room temperature.
    pub fn c_sic() -> Self {
        CmcMaterial {
            designation: "C/SiC".to_string(),
            matrix_cracking_stress: 80.0,
            uts: 350.0,
            modulus: 85.0,
            thermal_conductivity: 10.0,
            cte: 1.5e-6,
            max_temperature: 1650.0,
            k_ic: 20.0,
        }
    }
    /// Critical crack length for brittle fracture: ac = (KIc / (Y * σ))² / π.
    pub fn critical_crack_length(&self, sigma: f64, y: f64) -> f64 {
        (self.k_ic / (y * sigma)).powi(2) / PI
    }
    /// Thermal shock resistance parameter: R = σ * (1 - ν) / (E * CTE).
    pub fn thermal_shock_resistance(&self, nu: f64) -> f64 {
        self.uts * (1.0 - nu) / (self.modulus * 1e3 * self.cte)
    }
}
/// Honeycomb sandwich panel structural properties.
pub struct HoneycombSandwich {
    /// Core cell size (mm).
    pub cell_size: f64,
    /// Core density (kg/m³).
    pub core_density: f64,
    /// Ribbon direction shear modulus GRT (MPa).
    pub g_rt: f64,
    /// Transverse direction shear modulus GLT (MPa).
    pub g_lt: f64,
    /// Face sheet thickness tf (mm).
    pub face_thickness: f64,
    /// Core thickness tc (mm).
    pub core_thickness: f64,
    /// Face sheet modulus Ef (GPa).
    pub face_modulus: f64,
}
impl HoneycombSandwich {
    /// Create an aluminum honeycomb sandwich panel.
    pub fn aluminum_honeycomb(
        cell_size: f64,
        core_density: f64,
        face_thickness: f64,
        core_thickness: f64,
    ) -> Self {
        let g_rt = 100.0 * (core_density / 100.0).powf(1.5);
        let g_lt = 50.0 * (core_density / 100.0).powf(1.5);
        HoneycombSandwich {
            cell_size,
            core_density,
            g_rt,
            g_lt,
            face_thickness,
            core_thickness,
            face_modulus: 70.0,
        }
    }
    /// Total panel thickness h (mm).
    pub fn total_thickness(&self) -> f64 {
        2.0 * self.face_thickness + self.core_thickness
    }
    /// Equivalent bending stiffness D (N·mm) per unit width.
    /// D ≈ Ef * tf * d² / 2 where d = distance between face centroids.
    pub fn flexural_rigidity(&self) -> f64 {
        let d = self.core_thickness + self.face_thickness;
        self.face_modulus * 1000.0 * self.face_thickness * d * d / 2.0
    }
    /// Face sheet wrinkling stress σ_w (MPa).
    pub fn wrinkling_stress(&self) -> f64 {
        let ec = 2.0 * self.g_lt;
        0.5 * (self.face_modulus * 1000.0 * ec * self.g_rt).powf(1.0 / 3.0)
    }
    /// Core shear stress limit τ_max (MPa) for panel of width b under load P.
    pub fn shear_stress(&self, shear_force: f64) -> f64 {
        let h = self.total_thickness();
        shear_force / (h - 2.0 * self.face_thickness)
    }
}
/// Extended titanium alloy model with creep and high-temperature behavior.
pub struct TitaniumAlloyExtended {
    /// Base designation.
    pub designation: String,
    /// Room-temperature yield strength (MPa).
    pub sigma_y_rt: f64,
    /// Elevated-temperature yield strength at 300°C (MPa).
    pub sigma_y_300: f64,
    /// Creep activation energy Q (J/mol).
    pub creep_q: f64,
    /// Creep pre-exponential constant A.
    pub creep_a: f64,
    /// Creep stress exponent n.
    pub creep_n: f64,
    /// Oxidation rate constant at 600°C kp (mg²/cm⁴/h).
    pub kp_600: f64,
    /// Young's modulus at RT (GPa).
    pub e_rt: f64,
    /// Young's modulus temperature coefficient (GPa/°C).
    pub de_dt: f64,
    /// Density (kg/m³).
    pub density: f64,
}
impl TitaniumAlloyExtended {
    /// Ti-6Al-4V elevated temperature properties.
    pub fn ti64_elevated() -> Self {
        TitaniumAlloyExtended {
            designation: "Ti-6Al-4V (elevated)".to_string(),
            sigma_y_rt: 880.0,
            sigma_y_300: 700.0,
            creep_q: 250_000.0,
            creep_a: 1.0e8,
            creep_n: 4.0,
            kp_600: 2.0e-3,
            e_rt: 114.0,
            de_dt: -0.05,
            density: 4430.0,
        }
    }
    /// Ti-6Al-2Sn-4Zr-6Mo high-temperature alloy properties.
    pub fn ti6246() -> Self {
        TitaniumAlloyExtended {
            designation: "Ti-6-2-4-6".to_string(),
            sigma_y_rt: 1100.0,
            sigma_y_300: 900.0,
            creep_q: 270_000.0,
            creep_a: 5.0e7,
            creep_n: 3.5,
            kp_600: 1.0e-3,
            e_rt: 120.0,
            de_dt: -0.055,
            density: 4540.0,
        }
    }
    /// Young's modulus at temperature T (°C).
    pub fn youngs_modulus_at(&self, t_celsius: f64) -> f64 {
        (self.e_rt + self.de_dt * t_celsius).max(50.0)
    }
    /// Interpolated yield strength at temperature T (°C), linear between RT and 300°C.
    pub fn yield_strength_at(&self, t_celsius: f64) -> f64 {
        let t = t_celsius.clamp(20.0, 300.0);
        self.sigma_y_rt + (self.sigma_y_300 - self.sigma_y_rt) * (t - 20.0) / 280.0
    }
    /// Steady-state creep rate ε̇ = A * σ^n * exp(-Q/RT) (s⁻¹).
    pub fn creep_rate(&self, sigma_mpa: f64, t_kelvin: f64) -> f64 {
        const R: f64 = 8.314;
        self.creep_a * sigma_mpa.powf(self.creep_n) * (-self.creep_q / (R * t_kelvin)).exp()
    }
    /// Time to 0.2% creep strain at constant stress and temperature.
    pub fn time_to_02pct_creep(&self, sigma_mpa: f64, t_kelvin: f64) -> f64 {
        let eps_dot = self.creep_rate(sigma_mpa, t_kelvin);
        if eps_dot < 1e-30 {
            return f64::INFINITY;
        }
        0.002 / eps_dot
    }
    /// Oxidation mass gain (mg/cm²) after `t` hours at 600°C (parabolic law).
    pub fn oxidation_gain_600(&self, t_hours: f64) -> f64 {
        (self.kp_600 * t_hours).sqrt()
    }
}
/// CFRP laminate using Classical Lamination Theory (CLT).
pub struct CfrpLaminate {
    /// Plies from bottom to top.
    pub plies: Vec<CfrpPly>,
    /// Total thickness (mm).
    pub thickness: f64,
}
impl CfrpLaminate {
    /// Create a laminate from a list of plies.
    pub fn new(plies: Vec<CfrpPly>) -> Self {
        let thickness = plies.iter().map(|p| p.thickness).sum();
        CfrpLaminate { plies, thickness }
    }
    /// Create a symmetric quasi-isotropic layup \[0/45/-45/90\]s with given ply thickness.
    pub fn quasi_isotropic(ply_thickness: f64) -> Self {
        let angles = [0.0, 45.0, -45.0, 90.0, 90.0, -45.0, 45.0, 0.0];
        let plies = angles
            .iter()
            .map(|&a| CfrpPly::im7_8552(a, ply_thickness))
            .collect();
        Self::new(plies)
    }
    /// Compute the A-matrix (in-plane stiffness, N/mm) as \[A11, A22, A12, A66\].
    pub fn a_matrix(&self) -> [f64; 4] {
        let mut a11 = 0.0f64;
        let mut a22 = 0.0f64;
        let mut a12 = 0.0f64;
        let mut a66 = 0.0f64;
        for ply in &self.plies {
            let q = ply.transformed_stiffness();
            let t = ply.thickness;
            a11 += q[0] * t;
            a22 += q[1] * t;
            a12 += q[2] * t;
            a66 += q[3] * t;
        }
        [a11, a22, a12, a66]
    }
    /// Effective in-plane longitudinal stiffness Ex (GPa).
    pub fn effective_ex(&self) -> f64 {
        let [a11, _a22, a12, _a66] = self.a_matrix();
        (a11 - a12 * a12 / a11) / self.thickness
    }
    /// First ply failure index using Tsai-Hill criterion for each ply.
    /// Returns the minimum reserve factor (>1 = no failure).
    pub fn first_ply_failure_index(&self, nx: f64, ny: f64, nxy: f64) -> f64 {
        self.plies
            .iter()
            .map(|ply| {
                let q = ply.transformed_stiffness();
                let sigma1 = q[0] * nx / ply.thickness;
                let sigma2 = q[1] * ny / ply.thickness;
                let tau12 = q[3] * nxy / ply.thickness;
                let fi = (sigma1 / ply.f1t).powi(2) - (sigma1 * sigma2) / (ply.f1t * ply.f1t)
                    + (sigma2 / ply.f2t).powi(2)
                    + (tau12 / ply.f6).powi(2);
                if fi > 0.0 {
                    1.0 / fi.sqrt()
                } else {
                    f64::INFINITY
                }
            })
            .fold(f64::INFINITY, f64::min)
    }
}
/// Aerogel thermal insulation: Knudsen effect, effective conductivity.
pub struct AerogelInsulation {
    /// Solid skeleton thermal conductivity (W/m·K).
    pub k_solid: f64,
    /// Radiative conductivity at 25°C (W/m·K).
    pub k_rad_25: f64,
    /// Porosity (0–1).
    pub porosity: f64,
    /// Mean pore diameter (nm).
    pub pore_diameter_nm: f64,
    /// Density (kg/m³).
    pub density: f64,
    /// Operating temperature (°C).
    pub temperature: f64,
    /// Gas (air) thermal conductivity at STP (W/m·K).
    pub k_gas_stp: f64,
}
impl AerogelInsulation {
    /// Silica aerogel blanket at 25°C, ambient pressure.
    pub fn silica_aerogel() -> Self {
        AerogelInsulation {
            k_solid: 0.012,
            k_rad_25: 0.003,
            porosity: 0.95,
            pore_diameter_nm: 20.0,
            density: 120.0,
            temperature: 25.0,
            k_gas_stp: 0.0257,
        }
    }
    /// Pyrogel XT (Aspen) high-temperature aerogel blanket.
    pub fn pyrogel_xt() -> Self {
        AerogelInsulation {
            k_solid: 0.015,
            k_rad_25: 0.005,
            porosity: 0.92,
            pore_diameter_nm: 15.0,
            density: 170.0,
            temperature: 25.0,
            k_gas_stp: 0.0257,
        }
    }
    /// Knudsen number Kn = λ_mfp / pore_diameter.
    /// λ_mfp (air at STP) ≈ 68 nm.
    pub fn knudsen_number(&self, pressure_pa: f64) -> f64 {
        let lambda_stp = 68.0e-9;
        let lambda = lambda_stp * 101_325.0 / pressure_pa;
        lambda / (self.pore_diameter_nm * 1.0e-9)
    }
    /// Gas conductivity reduced by Knudsen effect.
    pub fn gas_conductivity_knudsen(&self, pressure_pa: f64) -> f64 {
        let kn = self.knudsen_number(pressure_pa);
        self.k_gas_stp / (1.0 + 2.0 * kn)
    }
    /// Effective thermal conductivity at pressure and temperature.
    pub fn effective_conductivity(&self, pressure_pa: f64) -> f64 {
        let k_gas = self.gas_conductivity_knudsen(pressure_pa);
        let t_factor = ((self.temperature + 273.15) / 298.15).powi(3);
        let k_rad = self.k_rad_25 * t_factor;
        (1.0 - self.porosity) * self.k_solid + self.porosity * k_gas + k_rad
    }
    /// Vacuum thermal conductivity (≈ solid + radiative only).
    pub fn vacuum_conductivity(&self) -> f64 {
        let t_factor = ((self.temperature + 273.15) / 298.15).powi(3);
        let k_rad = self.k_rad_25 * t_factor;
        (1.0 - self.porosity) * self.k_solid + k_rad
    }
}
/// Reusable thermal protection system (TPS) tile properties.
pub struct ThermalProtection {
    /// Material designation (PICA, TUFI, RCC, FRCI).
    pub designation: String,
    /// Thermal conductivity at room temperature (W/m·K).
    pub k_room: f64,
    /// Thermal conductivity at max temperature (W/m·K).
    pub k_hot: f64,
    /// Specific heat at room temperature (J/kg·K).
    pub cp_room: f64,
    /// Emissivity ε (dimensionless, 0–1).
    pub emissivity: f64,
    /// Max use temperature (°C).
    pub max_temp: f64,
    /// Density (kg/m³).
    pub density: f64,
}
impl ThermalProtection {
    /// PICA (Phenolic Impregnated Carbon Ablator) properties.
    pub fn pica() -> Self {
        ThermalProtection {
            designation: "PICA".to_string(),
            k_room: 0.18,
            k_hot: 0.45,
            cp_room: 1700.0,
            emissivity: 0.9,
            max_temp: 1650.0,
            density: 280.0,
        }
    }
    /// Interpolated thermal conductivity at temperature T (°C).
    pub fn conductivity_at(&self, t: f64) -> f64 {
        let t_frac = (t / self.max_temp).min(1.0);
        self.k_room + (self.k_hot - self.k_room) * t_frac
    }
    /// Radiation heat rejection rate q_rad (W/m²) at surface temp T (K).
    pub fn radiation_heat_flux(&self, t_k: f64) -> f64 {
        const SIGMA: f64 = 5.67e-8;
        self.emissivity * SIGMA * t_k.powi(4)
    }
    /// Required thickness to limit back-face temperature (simplified 1D conduction).
    pub fn required_thickness(&self, q_in: f64, dt_max: f64) -> f64 {
        let k_avg = (self.k_room + self.k_hot) / 2.0;
        k_avg * dt_max / q_in
    }
}
/// Carbon-carbon (C/C) composite with oxidation protection model.
pub struct CarbonCarbonComposite {
    /// Fiber volume fraction.
    pub vf: f64,
    /// Axial (fiber direction) thermal conductivity (W/m·K).
    pub k_axial: f64,
    /// Transverse thermal conductivity (W/m·K).
    pub k_transverse: f64,
    /// Room-temperature tensile strength (MPa).
    pub strength_rt: f64,
    /// Tensile strength at 2000°C in inert atmosphere (MPa).
    pub strength_2000: f64,
    /// SiC coating thickness (μm).
    pub sic_coat_thickness: f64,
    /// Oxidation onset temperature in air (°C).
    pub oxidation_onset_temp: f64,
    /// Oxidation rate constant at 600°C (mg/cm²/h).
    pub kox_600: f64,
    /// Density (kg/m³).
    pub density: f64,
}
impl CarbonCarbonComposite {
    /// 2D woven C/C composite with SiC-SiC coating.
    pub fn cc_2d_woven() -> Self {
        CarbonCarbonComposite {
            vf: 0.50,
            k_axial: 150.0,
            k_transverse: 30.0,
            strength_rt: 250.0,
            strength_2000: 350.0,
            sic_coat_thickness: 200.0,
            oxidation_onset_temp: 450.0,
            kox_600: 0.5,
            density: 1800.0,
        }
    }
    /// 3D C/C composite for rocket nozzle throats.
    pub fn cc_3d_nozzle() -> Self {
        CarbonCarbonComposite {
            vf: 0.55,
            k_axial: 120.0,
            k_transverse: 50.0,
            strength_rt: 300.0,
            strength_2000: 400.0,
            sic_coat_thickness: 300.0,
            oxidation_onset_temp: 400.0,
            kox_600: 0.3,
            density: 1950.0,
        }
    }
    /// Effective through-thickness thermal conductivity (parallel model).
    pub fn effective_conductivity(&self) -> f64 {
        self.vf * self.k_axial + (1.0 - self.vf) * self.k_transverse
    }
    /// Oxidation recession (μm) after `t` hours at `temp` °C.
    /// Uses simplified Arrhenius scaling from 600°C reference.
    pub fn oxidation_recession_um(&self, temp_celsius: f64, t_hours: f64) -> f64 {
        if temp_celsius < self.oxidation_onset_temp {
            return 0.0;
        }
        const Q_OX: f64 = 120_000.0;
        const R_GAS: f64 = 8.314;
        let t_ref = 873.15;
        let t_k = temp_celsius + 273.15;
        let k_scaled = self.kox_600 * ((-Q_OX / R_GAS) * (1.0 / t_k - 1.0 / t_ref)).exp();
        k_scaled * t_hours
    }
    /// SiC coating protection effectiveness: fraction of surface protected.
    pub fn coating_protection_factor(&self) -> f64 {
        1.0 - (-self.sic_coat_thickness / 100.0).exp()
    }
    /// Strength retention ratio at temperature T (°C) in inert atmosphere.
    pub fn strength_retention(&self, t_celsius: f64) -> f64 {
        let ratio = (t_celsius / 2000.0).clamp(0.0, 1.0);
        let str_t = self.strength_rt * (1.0 - ratio) + self.strength_2000 * ratio;
        str_t / self.strength_rt
    }
}
/// Foam-core sandwich panel: shear, bending, and indentation failure analysis.
pub struct FoamCoreSandwich {
    /// Core material designation.
    pub core_designation: String,
    /// Core density (kg/m³).
    pub core_density: f64,
    /// Core Young's modulus Ec (MPa).
    pub ec: f64,
    /// Core shear modulus Gc (MPa).
    pub gc: f64,
    /// Core compressive strength σc (MPa).
    pub sigma_c: f64,
    /// Core shear strength τc (MPa).
    pub tau_c: f64,
    /// Face sheet thickness tf (mm).
    pub tf: f64,
    /// Face sheet Young's modulus Ef (GPa).
    pub ef: f64,
    /// Face sheet tensile strength σf (MPa).
    pub sigma_f: f64,
    /// Core thickness tc (mm).
    pub tc: f64,
}
impl FoamCoreSandwich {
    /// Rohacell 51 WF PMI foam core sandwich with CFRP faces.
    pub fn rohacell51_cfrp(tf: f64, tc: f64) -> Self {
        FoamCoreSandwich {
            core_designation: "Rohacell 51 WF".to_string(),
            core_density: 52.0,
            ec: 70.0,
            gc: 19.0,
            sigma_c: 0.9,
            tau_c: 1.3,
            tf,
            ef: 70.0,
            sigma_f: 600.0,
            tc,
        }
    }
    /// Airex C70 PVC foam with glass/epoxy faces.
    pub fn airex_c70_gfrp(tf: f64, tc: f64) -> Self {
        FoamCoreSandwich {
            core_designation: "Airex C70".to_string(),
            core_density: 75.0,
            ec: 100.0,
            gc: 35.0,
            sigma_c: 1.4,
            tau_c: 1.7,
            tf,
            ef: 20.0,
            sigma_f: 300.0,
            tc,
        }
    }
    /// Total panel thickness (mm).
    pub fn total_thickness(&self) -> f64 {
        2.0 * self.tf + self.tc
    }
    /// Bending stiffness D (N·mm) per unit width (Euler-Bernoulli, face-dominated).
    pub fn bending_stiffness(&self) -> f64 {
        let d = self.tc + self.tf;
        self.ef * 1000.0 * self.tf * d * d / 2.0
    }
    /// Core shear stiffness S = Gc * tc (N/mm per unit width).
    pub fn shear_stiffness(&self) -> f64 {
        self.gc * self.tc
    }
    /// Maximum allowable distributed load (N/mm) limited by core shear failure.
    pub fn max_load_shear(&self, span_mm: f64) -> f64 {
        2.0 * self.tau_c * self.total_thickness() / (3.0 * span_mm)
    }
    /// Indentation failure load (N) using elastic foundation model.
    pub fn indentation_load(&self, indent_radius_mm: f64) -> f64 {
        PI * indent_radius_mm * (self.ec * self.sigma_f).sqrt() * self.tf
    }
    /// Face sheet wrinkling stress (MPa) using Euler buckling of face over core.
    pub fn wrinkling_stress_foam(&self) -> f64 {
        0.5 * (self.ef * 1000.0 * self.ec * self.gc).powf(1.0 / 3.0)
    }
}
/// Nickel-base superalloy: single crystal, directional solidification, and oxidation.
pub struct NickelSuperalloyExtended {
    /// Alloy designation.
    pub designation: String,
    /// Volume fraction of γ' precipitate.
    pub gamma_prime_vf: f64,
    /// γ' lattice misfit δ (dimensionless).
    pub misfit: f64,
    /// Room-temperature tensile strength (MPa).
    pub uts_rt: f64,
    /// Tensile strength at 1000°C (MPa).
    pub uts_1000: f64,
    /// Creep rupture life at 137 MPa, 982°C (hours).
    pub creep_life_ref: f64,
    /// Larson-Miller C constant.
    pub lm_c: f64,
    /// Oxidation parabolic constant at 1000°C (mg²/cm⁴/h).
    pub kp_1000: f64,
    /// MCrAlY bond coat thickness (μm).
    pub bond_coat_thickness: f64,
    /// Density (kg/m³).
    pub density: f64,
}
impl NickelSuperalloyExtended {
    /// CMSX-4 single-crystal superalloy properties.
    pub fn cmsx4() -> Self {
        NickelSuperalloyExtended {
            designation: "CMSX-4".to_string(),
            gamma_prime_vf: 0.70,
            misfit: -0.002,
            uts_rt: 1400.0,
            uts_1000: 950.0,
            creep_life_ref: 1000.0,
            lm_c: 21.5,
            kp_1000: 5.0e-5,
            bond_coat_thickness: 100.0,
            density: 8700.0,
        }
    }
    /// Rene 142 directionally-solidified superalloy.
    pub fn rene142() -> Self {
        NickelSuperalloyExtended {
            designation: "Rene 142".to_string(),
            gamma_prime_vf: 0.65,
            misfit: -0.0015,
            uts_rt: 1350.0,
            uts_1000: 880.0,
            creep_life_ref: 800.0,
            lm_c: 21.0,
            kp_1000: 6.0e-5,
            bond_coat_thickness: 80.0,
            density: 8820.0,
        }
    }
    /// Larson-Miller creep parameter P = T*(C + log10(tr)).
    pub fn larson_miller_parameter(&self, t_kelvin: f64, tr_hours: f64) -> f64 {
        t_kelvin * (self.lm_c + tr_hours.log10())
    }
    /// Estimated creep rupture life at (T, σ) relative to reference condition.
    pub fn estimated_rupture_life(&self, t_kelvin: f64, sigma_mpa: f64) -> f64 {
        let p_ref = self.larson_miller_parameter(1255.15, self.creep_life_ref);
        let stress_ratio = self.uts_rt / sigma_mpa;
        let tr = 10.0_f64.powf(p_ref / t_kelvin - self.lm_c);
        tr * stress_ratio.powi(2)
    }
    /// Oxidation mass gain at 1000°C after `t` hours (parabolic law).
    pub fn oxidation_gain(&self, t_hours: f64) -> f64 {
        (self.kp_1000 * t_hours).sqrt()
    }
    /// γ' strengthening increment (MPa) from Orowan bypassing.
    /// Δσ = M * G * b / (λ - 2r) where λ is particle spacing.
    pub fn gamma_prime_strengthening(&self, r_nm: f64, g_gpa: f64) -> f64 {
        const M: f64 = 3.06;
        const B: f64 = 0.254e-9;
        let r = r_nm * 1.0e-9;
        let lambda = r * (2.0 * PI / (3.0 * self.gamma_prime_vf)).sqrt();
        if lambda - 2.0 * r < 1.0e-12 {
            return 0.0;
        }
        M * g_gpa * 1.0e9 * B / (lambda - 2.0 * r) / 1.0e6
    }
    /// Thermal barrier coating effectiveness: fraction of temperature drop across TBC.
    pub fn tbc_temperature_fraction(&self, k_tbc: f64, k_substrate: f64, tbc_thick_um: f64) -> f64 {
        let tbc_m = tbc_thick_um * 1.0e-6;
        let sub_m = self.bond_coat_thickness * 1.0e-6;
        let r_tbc = tbc_m / k_tbc;
        let r_sub = sub_m / k_substrate;
        r_tbc / (r_tbc + r_sub)
    }
}
/// Ablative heat shield material model (Fay-Riddell recession rate model).
pub struct AblativeHeatShield {
    /// Material designation (e.g., "AVCOAT", "SLA-561V").
    pub designation: String,
    /// Effective heat of ablation H_eff (MJ/kg).
    pub h_eff: f64,
    /// Char layer thermal conductivity kc (W/m·K).
    pub k_char: f64,
    /// Virgin material thermal conductivity kv (W/m·K).
    pub k_virgin: f64,
    /// Char density ρc (kg/m³).
    pub rho_char: f64,
    /// Virgin density ρv (kg/m³).
    pub rho_virgin: f64,
    /// Char specific heat cp_c (J/kg·K).
    pub cp_char: f64,
    /// Pyrolysis decomposition temperature (°C).
    pub pyrolysis_temp: f64,
}
impl AblativeHeatShield {
    /// AVCOAT 5026-39G/HC ablator properties (Apollo era).
    pub fn avcoat() -> Self {
        AblativeHeatShield {
            designation: "AVCOAT 5026-39G/HC".to_string(),
            h_eff: 32.0,
            k_char: 0.5,
            k_virgin: 0.3,
            rho_char: 240.0,
            rho_virgin: 520.0,
            cp_char: 1250.0,
            pyrolysis_temp: 350.0,
        }
    }
    /// Recession rate ṙ (m/s) given surface heat flux q_w (W/m²).
    /// ṙ = q_w / (ρv * H_eff * 1e6).
    pub fn recession_rate(&self, q_w: f64) -> f64 {
        q_w / (self.rho_virgin * self.h_eff * 1e6)
    }
    /// Char layer thickness after time `t` (s) at constant heat flux.
    pub fn char_thickness(&self, q_w: f64, t: f64) -> f64 {
        self.recession_rate(q_w) * t
    }
    /// Mass loss per unit area (kg/m²) for given heat flux and time.
    pub fn mass_loss(&self, q_w: f64, t: f64) -> f64 {
        self.rho_virgin * self.char_thickness(q_w, t)
    }
}
/// Ti-6Al-4V titanium alloy material model.
pub struct TitaniumAlloy {
    /// Designation.
    pub designation: String,
    /// Young's modulus E (GPa).
    pub e: f64,
    /// 0.2% yield strength (MPa).
    pub sigma_y: f64,
    /// Ultimate tensile strength UTS (MPa).
    pub uts: f64,
    /// Fatigue threshold ΔK_th (MPa·√m).
    pub delta_k_th: f64,
    /// Paris law coefficient C (m/cycle per (MPa·√m)^m).
    pub paris_c: f64,
    /// Paris law exponent m.
    pub paris_m: f64,
    /// Density (kg/m³).
    pub density: f64,
    /// α-phase volume fraction.
    pub alpha_fraction: f64,
}
impl TitaniumAlloy {
    /// Ti-6Al-4V mill-annealed standard properties.
    pub fn ti64() -> Self {
        TitaniumAlloy {
            designation: "Ti-6Al-4V".to_string(),
            e: 114.0,
            sigma_y: 880.0,
            uts: 950.0,
            delta_k_th: 3.0,
            paris_c: 1.0e-11,
            paris_m: 3.2,
            density: 4430.0,
            alpha_fraction: 0.85,
        }
    }
    /// Paris-Erdogan fatigue crack growth rate da/dN (m/cycle).
    pub fn crack_growth_rate(&self, delta_k: f64) -> f64 {
        if delta_k < self.delta_k_th {
            return 0.0;
        }
        self.paris_c * delta_k.powf(self.paris_m)
    }
    /// Number of cycles to grow crack from `a0` to `af` (m) under ΔK = Δσ * sqrt(π * a).
    pub fn fatigue_life_cycles(&self, a0: f64, af: f64, delta_sigma: f64) -> f64 {
        let steps = 1000usize;
        let da = (af - a0) / steps as f64;
        let mut n = 0.0f64;
        for i in 0..steps {
            let a = a0 + (i as f64 + 0.5) * da;
            let dk = delta_sigma * (PI * a).sqrt();
            let dadn = self.crack_growth_rate(dk);
            if dadn > 0.0 {
                n += da / dadn;
            }
        }
        n
    }
    /// SCC threshold KISCC ≈ 0.5 * KIc. Returns critical stress for SCC initiation.
    pub fn scc_critical_stress(&self, a: f64, k_ic: f64) -> f64 {
        let k_scc = 0.5 * k_ic;
        k_scc / (PI * a).sqrt()
    }
}
/// Detailed ablation char/virgin material model tracking recession front.
pub struct AblationCharModel {
    /// Virgin material density (kg/m³).
    pub rho_virgin: f64,
    /// Char density (kg/m³).
    pub rho_char: f64,
    /// Effective heat of ablation (MJ/kg).
    pub h_eff: f64,
    /// Char thermal conductivity (W/m·K).
    pub k_char: f64,
    /// Virgin thermal conductivity (W/m·K).
    pub k_virgin: f64,
    /// Virgin specific heat (J/kg·K).
    pub cp_virgin: f64,
    /// Char specific heat (J/kg·K).
    pub cp_char: f64,
    /// Pyrolysis start temperature (°C).
    pub t_pyrolysis_start: f64,
    /// Pyrolysis completion temperature (°C).
    pub t_pyrolysis_end: f64,
    /// Pyrolysis enthalpy (J/kg).
    pub h_pyrolysis: f64,
}
impl AblationCharModel {
    /// Phenolic impregnated carbon ablator (PICA) char model.
    pub fn pica_char() -> Self {
        AblationCharModel {
            rho_virgin: 270.0,
            rho_char: 180.0,
            h_eff: 28.0,
            k_char: 0.45,
            k_virgin: 0.18,
            cp_virgin: 1700.0,
            cp_char: 1250.0,
            t_pyrolysis_start: 400.0,
            t_pyrolysis_end: 700.0,
            h_pyrolysis: 2.3e6,
        }
    }
    /// Carbon phenolic ablator for ballistic re-entry.
    pub fn carbon_phenolic() -> Self {
        AblationCharModel {
            rho_virgin: 1450.0,
            rho_char: 1250.0,
            h_eff: 42.0,
            k_char: 2.5,
            k_virgin: 1.2,
            cp_virgin: 1500.0,
            cp_char: 1800.0,
            t_pyrolysis_start: 350.0,
            t_pyrolysis_end: 600.0,
            h_pyrolysis: 3.5e6,
        }
    }
    /// Charring fraction: fraction converted to char given surface temperature.
    pub fn char_fraction(&self, t_surface: f64) -> f64 {
        if t_surface <= self.t_pyrolysis_start {
            return 0.0;
        }
        if t_surface >= self.t_pyrolysis_end {
            return 1.0;
        }
        (t_surface - self.t_pyrolysis_start) / (self.t_pyrolysis_end - self.t_pyrolysis_start)
    }
    /// Effective density at a given char fraction.
    pub fn effective_density(&self, char_frac: f64) -> f64 {
        let f = char_frac.clamp(0.0, 1.0);
        self.rho_virgin * (1.0 - f) + self.rho_char * f
    }
    /// Effective thermal conductivity at a given char fraction.
    pub fn effective_conductivity(&self, char_frac: f64) -> f64 {
        let f = char_frac.clamp(0.0, 1.0);
        self.k_virgin * (1.0 - f) + self.k_char * f
    }
    /// Surface recession velocity (m/s) at given heat flux q_w (W/m²).
    pub fn recession_velocity(&self, q_w: f64) -> f64 {
        q_w / (self.rho_virgin * self.h_eff * 1.0e6)
    }
    /// Total char thickness (m) after time `t` (s) at steady heat flux.
    pub fn char_thickness_at_time(&self, q_w: f64, t: f64) -> f64 {
        self.recession_velocity(q_w) * t
    }
    /// Pyrolysis gas mass flux (kg/m²/s): ṁg = (ρv - ρc) * ṡ.
    pub fn pyrolysis_gas_flux(&self, q_w: f64) -> f64 {
        let s_dot = self.recession_velocity(q_w);
        (self.rho_virgin - self.rho_char) * s_dot
    }
    /// Backface temperature estimate using 1D steady-state conduction through char.
    pub fn backface_temperature(&self, q_w: f64, char_thick: f64, t_front: f64) -> f64 {
        if char_thick < 1e-12 {
            return t_front;
        }
        t_front - q_w * char_thick / self.k_char
    }
}
/// Full ABD matrix computation for composite laminate (CLT).
pub struct CompositeLaminateClt {
    /// Plies with angle (deg) and thickness (mm), E1/E2/G12 (GPa), nu12.
    pub plies: Vec<(f64, f64, f64, f64, f64, f64)>,
}
impl CompositeLaminateClt {
    /// Create a new CLT laminate from angle/thickness/moduli tuples.
    /// Each tuple: `(angle_deg, thickness_mm, e1_gpa, e2_gpa, g12_gpa, nu12)`.
    pub fn new(plies: Vec<(f64, f64, f64, f64, f64, f64)>) -> Self {
        CompositeLaminateClt { plies }
    }
    /// Create a cross-ply \[0/90\]ns from ply count and thickness.
    pub fn cross_ply(
        n_pairs: usize,
        ply_thickness_mm: f64,
        e1: f64,
        e2: f64,
        g12: f64,
        nu12: f64,
    ) -> Self {
        let mut plies = Vec::new();
        for _ in 0..n_pairs {
            plies.push((0.0, ply_thickness_mm, e1, e2, g12, nu12));
            plies.push((90.0, ply_thickness_mm, e1, e2, g12, nu12));
        }
        Self::new(plies)
    }
    /// Total laminate thickness (mm).
    pub fn total_thickness_mm(&self) -> f64 {
        self.plies.iter().map(|p| p.1).sum()
    }
    /// Compute A11 component of ABD matrix (N/mm).
    pub fn a11(&self) -> f64 {
        self.plies
            .iter()
            .map(|p| {
                let (angle, t, e1, e2, g12, nu12) = *p;
                let nu21 = nu12 * e2 / e1;
                let denom = 1.0 - nu12 * nu21;
                let q11 = e1 / denom;
                let q22 = e2 / denom;
                let q12 = nu12 * e2 / denom;
                let q66 = g12;
                let c = angle.to_radians().cos();
                let s = angle.to_radians().sin();
                let c2 = c * c;
                let s2 = s * s;
                let c4 = c2 * c2;
                let s4 = s2 * s2;
                let q11b = q11 * c4 + 2.0 * (q12 + 2.0 * q66) * s2 * c2 + q22 * s4;
                q11b * t
            })
            .sum()
    }
    /// Compute A22 component (N/mm).
    pub fn a22(&self) -> f64 {
        self.plies
            .iter()
            .map(|p| {
                let (angle, t, e1, e2, g12, nu12) = *p;
                let nu21 = nu12 * e2 / e1;
                let denom = 1.0 - nu12 * nu21;
                let q11 = e1 / denom;
                let q22 = e2 / denom;
                let q12 = nu12 * e2 / denom;
                let q66 = g12;
                let c = angle.to_radians().cos();
                let s = angle.to_radians().sin();
                let c2 = c * c;
                let s2 = s * s;
                let c4 = c2 * c2;
                let s4 = s2 * s2;
                let q22b = q11 * s4 + 2.0 * (q12 + 2.0 * q66) * s2 * c2 + q22 * c4;
                q22b * t
            })
            .sum()
    }
    /// Interlaminar shear stress (MPa) estimate under transverse shear V (N/mm).
    /// τ_max ≈ 1.5 * V / h for rectangular cross-section.
    pub fn interlaminar_shear_stress(&self, v_per_width: f64) -> f64 {
        let h = self.total_thickness_mm();
        1.5 * v_per_width / h
    }
    /// Tsai-Wu failure index (simplified plane stress) for the critical ply.
    /// Returns max FI over all plies (FI > 1 = failure).
    pub fn tsai_wu_index(
        &self,
        nx: f64,
        ny: f64,
        nxy: f64,
        xt: f64,
        xc: f64,
        yt: f64,
        yc: f64,
        s: f64,
    ) -> f64 {
        self.plies
            .iter()
            .map(|p| {
                let (angle, t, e1, e2, g12, nu12) = *p;
                let nu21 = nu12 * e2 / e1;
                let denom = 1.0 - nu12 * nu21;
                let q11 = e1 / denom;
                let q22 = e2 / denom;
                let q12 = nu12 * e2 / denom;
                let q66 = g12;
                let c = angle.to_radians().cos();
                let sc = angle.to_radians().sin();
                let c2 = c * c;
                let s2 = sc * sc;
                let c4 = c2 * c2;
                let s4 = s2 * s2;
                let q11b = q11 * c4 + 2.0 * (q12 + 2.0 * q66) * s2 * c2 + q22 * s4;
                let q22b = q11 * s4 + 2.0 * (q12 + 2.0 * q66) * s2 * c2 + q22 * c4;
                let q66b = (q11 + q22 - 2.0 * q12 - 2.0 * q66) * s2 * c2 + q66 * (c4 + s4);
                let sigma1 = q11b * nx / t;
                let sigma2 = q22b * ny / t;
                let tau = q66b * nxy / t;
                let f1 = 1.0 / xt - 1.0 / xc;
                let f2 = 1.0 / yt - 1.0 / yc;
                let f11 = 1.0 / (xt * xc);
                let f22 = 1.0 / (yt * yc);
                let f66 = 1.0 / (s * s);
                f1 * sigma1
                    + f2 * sigma2
                    + f11 * sigma1.powi(2)
                    + f22 * sigma2.powi(2)
                    + f66 * tau.powi(2)
            })
            .fold(f64::NEG_INFINITY, f64::max)
    }
    /// Delamination onset criterion: Ye delamination criterion for interlaminar stresses.
    /// Returns reserve factor (>1 = safe).
    pub fn ye_delamination_rf(
        &self,
        sigma_z: f64,
        tau_xz: f64,
        tau_yz: f64,
        zzt: f64,
        szx: f64,
        szy: f64,
    ) -> f64 {
        let fi = (sigma_z / zzt).powi(2) + (tau_xz / szx).powi(2) + (tau_yz / szy).powi(2);
        if fi > 0.0 {
            1.0 / fi.sqrt()
        } else {
            f64::INFINITY
        }
    }
}
/// Aircraft fatigue load spectrum (ground-air-ground cycle model).
pub struct FatigueSpectrum {
    /// Number of design service goals (flight hours).
    pub dsg_hours: f64,
    /// Flights per hour (typical utilization).
    pub flights_per_hour: f64,
    /// Maneuver load factor exceedance table: (g_level, exceedances_per_1000_hrs).
    pub exceedance_data: Vec<(f64, f64)>,
    /// Gust spectrum intensity σ_g (g).
    pub gust_sigma: f64,
}
impl FatigueSpectrum {
    /// Create a standard medium transport aircraft spectrum.
    pub fn medium_transport() -> Self {
        FatigueSpectrum {
            dsg_hours: 120_000.0,
            flights_per_hour: 0.5,
            exceedance_data: vec![
                (1.0, 20_000.0),
                (1.5, 2_000.0),
                (2.0, 200.0),
                (2.5, 20.0),
                (3.0, 2.0),
                (3.5, 0.2),
            ],
            gust_sigma: 0.3,
        }
    }
    /// Total number of flights over DSG.
    pub fn total_flights(&self) -> f64 {
        self.dsg_hours * self.flights_per_hour
    }
    /// Exceedances at given g-level over DSG (interpolated from table).
    pub fn exceedances_at(&self, g_level: f64) -> f64 {
        if self.exceedance_data.is_empty() {
            return 0.0;
        }
        if g_level <= self.exceedance_data[0].0 {
            return self.exceedance_data[0].1 * self.dsg_hours / 1000.0;
        }
        for i in 0..self.exceedance_data.len() - 1 {
            let (g1, e1) = self.exceedance_data[i];
            let (g2, e2) = self.exceedance_data[i + 1];
            if g_level >= g1 && g_level <= g2 {
                let frac = (g_level - g1) / (g2 - g1);
                let e_interp = e1 + frac * (e2 - e1);
                return e_interp * self.dsg_hours / 1000.0;
            }
        }
        0.0
    }
    /// Miner's rule cumulative damage from multiple load levels.
    /// `s_n_pairs` is `&[(stress_range, allowed_cycles_at_stress)]`.
    pub fn miners_rule_damage(&self, s_n_pairs: &[(f64, f64)]) -> f64 {
        s_n_pairs
            .iter()
            .map(|(ni, ni_allowed)| {
                if *ni_allowed > 0.0 {
                    ni / ni_allowed
                } else {
                    f64::INFINITY
                }
            })
            .sum()
    }
}
/// Space environment degradation: radiation darkening, outgassing, atomic oxygen.
pub struct SpaceEnvironmentEffects {
    /// Material designation.
    pub designation: String,
    /// Total ionizing dose (TID) threshold for property change (Gy).
    pub tid_threshold: f64,
    /// Optical absorptance change per Mrad (1 Mrad = 10 kGy).
    pub d_alpha_per_mrad: f64,
    /// Outgassing total mass loss fraction (TML, %).
    pub tml_percent: f64,
    /// Collected volatile condensable material (CVCM, %).
    pub cvcm_percent: f64,
    /// Atomic oxygen (AO) erosion yield (cm³/atom).
    pub ao_erosion_yield: f64,
    /// AO fluence at LEO per year (atoms/cm²/year).
    pub ao_fluence_leo: f64,
    /// Material initial absorptance α0.
    pub alpha_initial: f64,
    /// Material emissivity ε.
    pub emissivity: f64,
}
impl SpaceEnvironmentEffects {
    /// Kapton HN polyimide space environment properties.
    pub fn kapton_hn() -> Self {
        SpaceEnvironmentEffects {
            designation: "Kapton HN".to_string(),
            tid_threshold: 1.0e6,
            d_alpha_per_mrad: 0.01,
            tml_percent: 0.98,
            cvcm_percent: 0.01,
            ao_erosion_yield: 3.0e-24,
            ao_fluence_leo: 8.0e20,
            alpha_initial: 0.32,
            emissivity: 0.64,
        }
    }
    /// SolarBlack paint (high absorptance thermal control).
    pub fn solar_black_paint() -> Self {
        SpaceEnvironmentEffects {
            designation: "Solar Black Paint".to_string(),
            tid_threshold: 5.0e5,
            d_alpha_per_mrad: 0.005,
            tml_percent: 0.20,
            cvcm_percent: 0.02,
            ao_erosion_yield: 1.0e-25,
            ao_fluence_leo: 8.0e20,
            alpha_initial: 0.96,
            emissivity: 0.87,
        }
    }
    /// Absorptance after radiation dose `dose_gy` (Gy).
    pub fn radiation_darkened_absorptance(&self, dose_gy: f64) -> f64 {
        let dose_mrad = dose_gy / 10_000.0;
        (self.alpha_initial + self.d_alpha_per_mrad * dose_mrad).min(1.0)
    }
    /// Atomic oxygen erosion depth (μm) after `t_years` years in LEO.
    pub fn ao_erosion_depth_um(&self, t_years: f64, density_g_cm3: f64) -> f64 {
        let fluence = self.ao_fluence_leo * t_years;
        let mass_loss = self.ao_erosion_yield * fluence;
        mass_loss * density_g_cm3 * 10_000.0
    }
    /// Check if outgassing meets NASA ASTM E595 requirement (TML < 1%, CVCM < 0.1%).
    pub fn passes_outgassing(&self) -> bool {
        self.tml_percent < 1.0 && self.cvcm_percent < 0.1
    }
    /// Solar absorptance-to-emissivity ratio α/ε (thermal balance parameter).
    pub fn alpha_over_epsilon(&self) -> f64 {
        self.alpha_initial / self.emissivity
    }
}
/// Damage tolerance analysis: initial flaw, inspection interval, repair threshold.
pub struct DamageToleranceFlaw {
    /// Initial detectable flaw size a_i (m).
    pub a_initial: f64,
    /// Critical flaw size a_cr (m) at design ultimate load.
    pub a_critical: f64,
    /// Repair threshold a_r (m).
    pub a_repair: f64,
    /// Inspection interval (flight hours).
    pub inspection_interval: f64,
    /// Probability of detection at repair threshold (0–1).
    pub pod_at_repair: f64,
}
impl DamageToleranceFlaw {
    /// Create a standard damage tolerance scenario.
    pub fn new(a_initial: f64, a_critical: f64, inspection_interval: f64) -> Self {
        DamageToleranceFlaw {
            a_initial,
            a_critical,
            a_repair: a_critical / 2.0,
            inspection_interval,
            pod_at_repair: 0.9,
        }
    }
    /// Check if the flaw size `a` exceeds the repair threshold.
    pub fn requires_repair(&self, a: f64) -> bool {
        a >= self.a_repair
    }
    /// Check if the flaw size `a` exceeds the critical size.
    pub fn is_critical(&self, a: f64) -> bool {
        a >= self.a_critical
    }
    /// Crack growth life factor: ratio of critical to initial flaw size.
    pub fn life_factor(&self) -> f64 {
        self.a_critical / self.a_initial
    }
    /// Required number of inspection intervals before crack becomes critical.
    /// (Assumes linear growth rate for simplicity.)
    pub fn inspections_before_critical(&self, growth_rate_per_hour: f64) -> f64 {
        if growth_rate_per_hour <= 0.0 {
            return f64::INFINITY;
        }
        let time_to_critical = (self.a_critical - self.a_initial) / growth_rate_per_hour;
        time_to_critical / self.inspection_interval
    }
}
/// Extended CMC model: matrix cracking, fiber pullout, BN interphase.
pub struct CeramicMatrixCompositeExtended {
    /// Fiber diameter df (μm).
    pub df: f64,
    /// Fiber volume fraction Vf.
    pub vf: f64,
    /// Fiber modulus Ef (GPa).
    pub ef: f64,
    /// Matrix modulus Em (GPa).
    pub em: f64,
    /// Interface sliding stress τi (MPa).
    pub tau_i: f64,
    /// Fiber strength σfu (MPa).
    pub sigma_fu: f64,
    /// Matrix cracking stress σmc (MPa).
    pub sigma_mc: f64,
    /// BN interphase thickness (nm).
    pub bn_thickness_nm: f64,
    /// Fiber pullout length lp (mm).
    pub pullout_length: f64,
    /// Composite density (kg/m³).
    pub density: f64,
}
impl CeramicMatrixCompositeExtended {
    /// Hi-Nicalon SiC/SiC with BN interphase.
    pub fn hi_nicalon_sic_sic() -> Self {
        CeramicMatrixCompositeExtended {
            df: 14.0,
            vf: 0.40,
            ef: 270.0,
            em: 350.0,
            tau_i: 50.0,
            sigma_fu: 2000.0,
            sigma_mc: 150.0,
            bn_thickness_nm: 300.0,
            pullout_length: 0.5,
            density: 2600.0,
        }
    }
    /// Tyranno SA C/SiC.
    pub fn tyranno_sa() -> Self {
        CeramicMatrixCompositeExtended {
            df: 7.5,
            vf: 0.45,
            ef: 340.0,
            em: 350.0,
            tau_i: 30.0,
            sigma_fu: 2300.0,
            sigma_mc: 120.0,
            bn_thickness_nm: 200.0,
            pullout_length: 1.0,
            density: 2800.0,
        }
    }
    /// Composite modulus (rule of mixtures) (GPa).
    pub fn composite_modulus(&self) -> f64 {
        self.ef * self.vf + self.em * (1.0 - self.vf)
    }
    /// Matrix cracking stress using ACK model (MPa).
    /// σ_mc = (12 * τi * Ef * Vf^2 * Γm * Ec / (Em^2 * df * (1-Vf)))^(1/3)
    pub fn ack_matrix_cracking_stress(&self, gamma_m: f64) -> f64 {
        let ec = self.composite_modulus();
        let df_m = self.df * 1.0e-6;
        let num = 12.0 * self.tau_i * self.ef * self.vf.powi(2) * gamma_m * ec;
        let den = self.em.powi(2) * df_m * 1000.0 * (1.0 - self.vf);
        if den < 1.0e-30 {
            return 0.0;
        }
        (num / den).powf(1.0 / 3.0)
    }
    /// Fiber pullout energy per unit area (J/m²).
    pub fn pullout_energy(&self) -> f64 {
        let df_m = self.df * 1.0e-6;
        let lp = self.pullout_length * 1.0e-3;
        self.vf * self.tau_i * lp.powi(2) / df_m
    }
    /// Composite ultimate tensile strength (MPa) using GLS model.
    pub fn ultimate_strength(&self) -> f64 {
        self.vf * self.sigma_fu * (1.0 - self.vf).powf(0.5)
    }
    /// Thermal shock resistance parameter R' (K) = σ_mc * k / (E * CTE).
    pub fn thermal_shock_r_prime(&self, k: f64, cte: f64) -> f64 {
        let ec = self.composite_modulus();
        self.sigma_mc * k / (ec * 1.0e9 * cte)
    }
}
/// A single lamina (ply) in a composite laminate.
#[derive(Debug, Clone)]
pub struct CfrpPly {
    /// Ply angle (degrees).
    pub angle_deg: f64,
    /// Ply thickness (mm).
    pub thickness: f64,
    /// Longitudinal modulus E1 (GPa).
    pub e1: f64,
    /// Transverse modulus E2 (GPa).
    pub e2: f64,
    /// In-plane shear modulus G12 (GPa).
    pub g12: f64,
    /// Major Poisson's ratio ν12.
    pub nu12: f64,
    /// Longitudinal tensile strength F1t (MPa).
    pub f1t: f64,
    /// Transverse tensile strength F2t (MPa).
    pub f2t: f64,
    /// In-plane shear strength F6 (MPa).
    pub f6: f64,
}
impl CfrpPly {
    /// Create a standard IM7/8552 carbon fiber epoxy ply.
    pub fn im7_8552(angle_deg: f64, thickness: f64) -> Self {
        CfrpPly {
            angle_deg,
            thickness,
            e1: 161.0,
            e2: 11.4,
            g12: 5.17,
            nu12: 0.32,
            f1t: 2560.0,
            f2t: 73.0,
            f6: 90.0,
        }
    }
    /// Transformed stiffness matrix Q̄ in global coordinates (GPa).
    /// Returns \[Q11bar, Q22bar, Q12bar, Q66bar, Q16bar, Q26bar\].
    pub fn transformed_stiffness(&self) -> [f64; 6] {
        let theta = self.angle_deg.to_radians();
        let c = theta.cos();
        let s = theta.sin();
        let nu21 = self.nu12 * self.e2 / self.e1;
        let denom = 1.0 - self.nu12 * nu21;
        let q11 = self.e1 / denom;
        let q22 = self.e2 / denom;
        let q12 = self.nu12 * self.e2 / denom;
        let q66 = self.g12;
        let c2 = c * c;
        let s2 = s * s;
        let c4 = c2 * c2;
        let s4 = s2 * s2;
        let cs = c * s;
        let q11b = q11 * c4 + 2.0 * (q12 + 2.0 * q66) * s2 * c2 + q22 * s4;
        let q22b = q11 * s4 + 2.0 * (q12 + 2.0 * q66) * s2 * c2 + q22 * c4;
        let q12b = (q11 + q22 - 4.0 * q66) * s2 * c2 + q12 * (c4 + s4);
        let q66b = (q11 + q22 - 2.0 * q12 - 2.0 * q66) * s2 * c2 + q66 * (c4 + s4);
        let q16b = (q11 - q12 - 2.0 * q66) * c * s2 * c - (q22 - q12 - 2.0 * q66) * cs * s2;
        let q26b = (q11 - q12 - 2.0 * q66) * cs * s2 - (q22 - q12 - 2.0 * q66) * c2 * cs;
        [q11b, q22b, q12b, q66b, q16b, q26b]
    }
}
/// Adhesive joint model: lap shear, T-peel, climbing drum peel.
pub struct AdhesiveJoint {
    /// Adhesive designation.
    pub designation: String,
    /// Adhesive shear modulus G (MPa).
    pub g_adhesive: f64,
    /// Adhesive tensile modulus E (MPa).
    pub e_adhesive: f64,
    /// Adhesive thickness t (mm).
    pub t_adhesive: f64,
    /// Mode I fracture energy GIc (J/m²).
    pub g_ic: f64,
    /// Mode II fracture energy GIIc (J/m²).
    pub g_iic: f64,
    /// Shear strength τ_u (MPa).
    pub tau_u: f64,
    /// Tensile strength σ_u (MPa).
    pub sigma_u: f64,
    /// Overlap length l (mm).
    pub overlap_length: f64,
}
impl AdhesiveJoint {
    /// FM 300 film adhesive (aerospace structural).
    pub fn fm300(t_adhesive: f64, overlap_length: f64) -> Self {
        AdhesiveJoint {
            designation: "FM 300".to_string(),
            g_adhesive: 650.0,
            e_adhesive: 2100.0,
            t_adhesive,
            g_ic: 1500.0,
            g_iic: 2500.0,
            tau_u: 40.0,
            sigma_u: 45.0,
            overlap_length,
        }
    }
    /// EC-3448 paste adhesive.
    pub fn ec3448(t_adhesive: f64, overlap_length: f64) -> Self {
        AdhesiveJoint {
            designation: "EC-3448".to_string(),
            g_adhesive: 400.0,
            e_adhesive: 1800.0,
            t_adhesive,
            g_ic: 800.0,
            g_iic: 1500.0,
            tau_u: 25.0,
            sigma_u: 30.0,
            overlap_length,
        }
    }
    /// Volkersen shear stress distribution in single-lap joint.
    /// Returns peak shear stress at joint end (MPa).
    pub fn volkersen_peak_shear(&self, e_sub: f64, t_sub: f64, load_n_per_mm: f64) -> f64 {
        let omega = ((self.g_adhesive / self.t_adhesive) * (2.0 / (e_sub * t_sub))).sqrt();
        let l = self.overlap_length;
        let cosh_val = (omega * l).cosh();
        let sinh_val = (omega * l).sinh();
        if sinh_val < 1.0e-12 {
            return load_n_per_mm / l;
        }
        load_n_per_mm * omega * cosh_val / sinh_val
    }
    /// T-peel fracture energy (J/m²) required to propagate delamination.
    pub fn t_peel_energy(&self, peel_force_n_per_mm: f64) -> f64 {
        2.0 * peel_force_n_per_mm * 1000.0
    }
    /// Climbing drum peel mode I energy release rate G_I (J/m²).
    pub fn climbing_drum_g1(&self, m_drum: f64, r_drum: f64) -> f64 {
        m_drum / r_drum
    }
    /// Reserve factor against shear failure.
    pub fn rf_shear(&self, applied_shear: f64) -> f64 {
        self.tau_u / applied_shear.max(1.0e-15)
    }
}
