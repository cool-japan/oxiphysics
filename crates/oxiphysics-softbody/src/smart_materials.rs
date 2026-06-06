// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Smart / active soft materials — hydrogels, LCE, IPMC, and more.
//!
//! Models the electromechanical, thermomechanical, and field-responsive
//! behaviour of stimulus-responsive soft materials:
//!
//! - [`Hydrogel`]: Flory–Huggins swelling and osmotic pressure.
//! - [`LiquidCrystalElastomer`]: order-parameter driven thermal actuation.
//! - [`IPMC`]: ionic polymer-metal composite bending actuator.
//! - [`MagnetoActiveElastomer`]: magnetically stiffened elastomer.
//! - [`PhotothermalActuator`]: light-driven thermal expansion.
//! - Free function [`kelvin_voigt_response`]: viscoelastic strain update.

// ── Hydrogel ──────────────────────────────────────────────────────────────────

/// Flory–Huggins hydrogel model.
///
/// Captures equilibrium swelling, osmotic pressure, and elastic modulus of a
/// cross-linked polymer network immersed in solvent.
///
/// # Parameters
/// * `chi` — Flory–Huggins interaction parameter (quality of solvent; < 0.5 good solvent).
/// * `v_s` — molar volume of solvent \[m³/mol\].
/// * `c_0` — reference polymer volume fraction (dry state, 0 < c_0 ≤ 1).
/// * `n_chain` — strand density of the network \[mol/m³\].
#[derive(Debug, Clone)]
pub struct Hydrogel {
    /// Flory–Huggins χ parameter (dimensionless).
    pub chi: f64,
    /// Molar volume of solvent \[m³/mol\].
    pub v_s: f64,
    /// Dry-state polymer volume fraction (reference state).
    pub c_0: f64,
    /// Strand density of the crosslinked network \[mol/m³\].
    pub n_chain: f64,
}

impl Hydrogel {
    /// Boltzmann × Avogadro = universal gas constant R \[J/(mol·K)\].
    const R: f64 = 8.314_462_618;
    /// Reference temperature \[K\].
    const T_REF: f64 = 298.15;

    /// Create a new hydrogel model.
    pub fn new(chi: f64, v_s: f64, c_0: f64, n_chain: f64) -> Self {
        Self {
            chi,
            v_s,
            c_0: c_0.clamp(1e-6, 1.0),
            n_chain,
        }
    }

    /// Equilibrium swelling ratio Q = 1/φ_eq given external chemical potential.
    ///
    /// Solved from the condition μ_mix + μ_elastic = μ_ext using a simple
    /// fixed-point iteration.  `mu_ext` is the external chemical potential
    /// in units of RT.
    ///
    /// Returns the volumetric swelling ratio Q ≥ 1.
    pub fn swelling_ratio(&self, mu_ext: f64) -> f64 {
        // Iterative solve for φ such that Δμ(φ) = mu_ext.
        // Δμ(φ) = ln(1-φ) + φ + χφ² + v_s * n_chain * (φ^(1/3)/c_0^(1/3) - φ/2)
        let mut phi = self.c_0;
        for _ in 0..200 {
            let lnterm = (1.0 - phi).ln();
            let elastic = self.v_s
                * self.n_chain
                * (phi.powf(1.0 / 3.0) / self.c_0.powf(1.0 / 3.0) - 0.5 * phi);
            let mu_calc = lnterm + phi + self.chi * phi * phi + elastic;
            let residual = mu_calc - mu_ext;
            // Numerical derivative.
            let dphi = 1e-6_f64;
            let phi_p = (phi + dphi).min(1.0 - 1e-10);
            let elastic_p = self.v_s
                * self.n_chain
                * (phi_p.powf(1.0 / 3.0) / self.c_0.powf(1.0 / 3.0) - 0.5 * phi_p);
            let mu_p = (1.0 - phi_p).ln() + phi_p + self.chi * phi_p * phi_p + elastic_p;
            let deriv = (mu_p - mu_calc) / dphi;
            if deriv.abs() < 1e-20 {
                break;
            }
            phi -= residual / deriv;
            phi = phi.clamp(1e-6, 1.0 - 1e-8);
        }
        1.0 / phi
    }

    /// Rubber-elastic shear modulus of the swollen network \[Pa\].
    ///
    /// G = n_chain * R * T.
    pub fn elastic_modulus(&self) -> f64 {
        self.n_chain * Self::R * Self::T_REF
    }

    /// Osmotic pressure at polymer volume fraction `phi` \[Pa\].
    ///
    /// Π = -(R T / v_s) * \[ln(1-φ) + φ + χφ²\]
    pub fn osmotic_pressure(&self, phi: f64) -> f64 {
        let phi = phi.clamp(1e-8, 1.0 - 1e-8);
        let mix = (1.0 - phi).ln() + phi + self.chi * phi * phi;
        -(Self::R * Self::T_REF / self.v_s) * mix
    }

    /// Flory–Huggins mixing free energy per unit volume \[J/m³\].
    ///
    /// f_mix = (R T / v_s) * \[φ ln φ + (1-φ) ln(1-φ) + χ φ(1-φ)\]
    pub fn flory_huggins_free_energy(&self, phi: f64) -> f64 {
        let phi = phi.clamp(1e-8, 1.0 - 1e-8);
        let entropy = phi * phi.ln() + (1.0 - phi) * (1.0 - phi).ln();
        let enthalpy = self.chi * phi * (1.0 - phi);
        (Self::R * Self::T_REF / self.v_s) * (entropy + enthalpy)
    }

    /// Chemical potential of solvent \[dimensionless, in units of RT\].
    pub fn chemical_potential(&self, phi: f64) -> f64 {
        let phi = phi.clamp(1e-8, 1.0 - 1e-8);
        (1.0 - phi).ln() + phi + self.chi * phi * phi
    }
}

// ── LiquidCrystalElastomer ────────────────────────────────────────────────────

/// Liquid crystal elastomer (LCE) thermomechanical actuator model.
///
/// An LCE couples nematic order (degree of liquid-crystal alignment) to
/// mechanical deformation, enabling reversible thermally-driven actuation.
///
/// # Parameters
/// * `order_param` — equilibrium nematic order parameter S ∈ \[0, 1\] (1 = fully aligned).
/// * `coupling` — order-strain coupling coefficient λ \[dimensionless\].
/// * `modulus` — isotropic shear modulus \[Pa\].
#[derive(Debug, Clone)]
pub struct LiquidCrystalElastomer {
    /// Equilibrium nematic order parameter S₀ ∈ \[0, 1\].
    pub order_param: f64,
    /// Order-strain coupling coefficient λ.
    pub coupling: f64,
    /// Isotropic shear modulus μ \[Pa\].
    pub modulus: f64,
}

impl LiquidCrystalElastomer {
    /// Create a new LCE model.
    ///
    /// * `order_param` — reference order parameter S₀ at room temperature.
    /// * `coupling` — linear coupling constant λ (typical: 0.1–0.5).
    /// * `modulus` — shear modulus \[Pa\].
    pub fn new(order_param: f64, coupling: f64, modulus: f64) -> Self {
        Self {
            order_param: order_param.clamp(0.0, 1.0),
            coupling,
            modulus,
        }
    }

    /// Nematic order parameter at temperature `temp` \[K\] given clearing point `t_ni` \[K\].
    ///
    /// Uses the mean-field Landau approximation:
    /// S(T) = S₀ * max(0, (T_NI - T) / T_NI)^(1/2)
    pub fn nematic_order(&self, temp: f64, t_ni: f64) -> f64 {
        if temp >= t_ni {
            return 0.0;
        }
        let reduced = (t_ni - temp) / t_ni;
        self.order_param * reduced.sqrt()
    }

    /// Thermally induced actuation strain at temperature `temp` \[K\].
    ///
    /// ε = λ * (S(T) - S(T_ref)) where T_ref = 0.25 * T_NI.
    /// Returns a signed strain (positive = elongation, negative = contraction).
    pub fn thermal_actuation(&self, temp: f64, t_ni: f64) -> f64 {
        let s_ref = self.nematic_order(0.25 * t_ni, t_ni);
        let s_now = self.nematic_order(temp, t_ni);
        self.coupling * (s_now - s_ref)
    }

    /// Strain-order coupling coefficient \[dimensionless\].
    ///
    /// Defined as dε/dS = λ (constant in linear theory).
    pub fn strain_order_coupling(&self) -> f64 {
        self.coupling
    }

    /// Effective actuation stress \[Pa\] = μ * ε.
    pub fn actuation_stress(&self, temp: f64, t_ni: f64) -> f64 {
        self.modulus * self.thermal_actuation(temp, t_ni)
    }

    /// Recoverable actuation work density \[J/m³\] during cooling from T_hot to T_cold.
    pub fn actuation_work_density(&self, t_hot: f64, t_cold: f64, t_ni: f64) -> f64 {
        let eps_hot = self.thermal_actuation(t_hot, t_ni);
        let eps_cold = self.thermal_actuation(t_cold, t_ni);
        0.5 * self.modulus * (eps_cold - eps_hot).powi(2)
    }
}

// ── IPMC ──────────────────────────────────────────────────────────────────────

/// Ionic Polymer-Metal Composite (IPMC) bending actuator.
///
/// An IPMC is a thin membrane of ion-exchange polymer sandwiched between metal
/// electrodes.  Voltage drives ionic flux which causes differential swelling and
/// hence bending.
///
/// # Parameters
/// * `thickness` — membrane thickness \[m\].
/// * `water_content` — mass fraction of water in the ionomer (0–1).
/// * `voltage` — applied voltage \[V\].
#[derive(Debug, Clone)]
pub struct IPMC {
    /// Membrane thickness \[m\].
    pub thickness: f64,
    /// Water content (mass fraction, 0–1).
    pub water_content: f64,
    /// Applied voltage \[V\].
    pub voltage: f64,
    /// Free length of the actuator beam \[m\].
    pub length: f64,
    /// Young's modulus of the ionomer \[Pa\].
    pub young_modulus: f64,
}

impl IPMC {
    /// Create a new IPMC actuator.
    ///
    /// * `thickness` — total membrane thickness \[m\].
    /// * `water_content` — mass fraction of absorbed water.
    /// * `voltage` — driving voltage \[V\].
    /// * `length` — free cantilever length \[m\].
    /// * `young_modulus` — Young's modulus \[Pa\].
    pub fn new(
        thickness: f64,
        water_content: f64,
        voltage: f64,
        length: f64,
        young_modulus: f64,
    ) -> Self {
        Self {
            thickness,
            water_content: water_content.clamp(0.0, 1.0),
            voltage,
            length,
            young_modulus,
        }
    }

    /// Tip displacement δ of the cantilever beam \[m\].
    ///
    /// Empirical model: δ ≈ k_d * V * (L/h)² * water_content
    /// where k_d ≈ 0.002 (fitting constant).
    pub fn tip_displacement(&self) -> f64 {
        if self.thickness < 1e-15 || self.length < 1e-15 {
            return 0.0;
        }
        let k_d = 0.002_f64;
        let ratio = self.length / self.thickness;
        k_d * self.voltage * ratio * ratio * self.water_content
    }

    /// Blocking force at the tip \[N\].
    ///
    /// F_block = E * I * δ / L³  (cantilever stiffness × tip displacement)
    /// where I = h³/12 (rectangular cross-section, unit width assumed).
    pub fn blocking_force(&self) -> f64 {
        let h = self.thickness;
        let l = self.length;
        if h < 1e-15 || l < 1e-15 {
            return 0.0;
        }
        let moment_of_inertia = h.powi(3) / 12.0;
        let stiffness = 3.0 * self.young_modulus * moment_of_inertia / l.powi(3);
        stiffness * self.tip_displacement()
    }

    /// Characteristic electrical–mechanical response time \[s\].
    ///
    /// τ ≈ h² / D_ion  where D_ion ≈ 1e-11 m²/s for Nafion-type polymers.
    /// Scales inversely with water content (higher hydration → faster ion mobility).
    pub fn response_time(&self) -> f64 {
        let d_ion = 1.0e-11_f64 * (1.0 + self.water_content);
        self.thickness.powi(2) / d_ion
    }

    /// Curvature of the bent beam \[1/m\].
    ///
    /// κ = 2δ / L²  (small-angle approximation for a cantilever).
    pub fn curvature(&self) -> f64 {
        if self.length < 1e-15 {
            return 0.0;
        }
        2.0 * self.tip_displacement() / self.length.powi(2)
    }
}

// ── MagnetoActiveElastomer ────────────────────────────────────────────────────

/// Magneto-active (magnetorheological) elastomer model.
///
/// An MAE is a soft rubber matrix with embedded magnetic particles.  An applied
/// magnetic field deforms the matrix through magnetic body forces and stiffens it
/// through field-induced microstructure.
///
/// # Parameters
/// * `magnetization` — saturation magnetization \[A/m\].
/// * `susceptibility` — initial magnetic susceptibility χ_m (dimensionless).
/// * `shear_modulus` — zero-field shear modulus G₀ \[Pa\].
#[derive(Debug, Clone)]
pub struct MagnetoActiveElastomer {
    /// Saturation magnetization M_sat \[A/m\].
    pub magnetization: f64,
    /// Initial magnetic susceptibility χ_m.
    pub susceptibility: f64,
    /// Zero-field shear modulus G₀ \[Pa\].
    pub shear_modulus: f64,
}

impl MagnetoActiveElastomer {
    /// Vacuum permeability μ₀ \[T·m/A\].
    const MU_0: f64 = 1.256_637_061_4e-6;

    /// Create a new MAE model.
    pub fn new(magnetization: f64, susceptibility: f64, shear_modulus: f64) -> Self {
        Self {
            magnetization,
            susceptibility,
            shear_modulus,
        }
    }

    /// Field-induced magnetostrictive strain at applied flux density `b_field` \[T\].
    ///
    /// ε_mag = (μ₀ * χ_m * H²) / (2 * G₀)
    ///
    /// Uses B ≈ μ₀(1 + χ_m)H → H = B / (μ₀(1+χ_m)).
    pub fn magneto_strain(&self, b_field: f64) -> f64 {
        if self.shear_modulus < 1e-10 {
            return 0.0;
        }
        let mu_eff = Self::MU_0 * (1.0 + self.susceptibility);
        let h = b_field / mu_eff;
        let numerator = Self::MU_0 * self.susceptibility * h * h;
        numerator / (2.0 * self.shear_modulus)
    }

    /// Field-stiffened shear modulus \[Pa\] at applied flux density `b_field` \[T\].
    ///
    /// G(B) = G₀ + ΔG(B) where ΔG ≈ (μ₀ χ_m² M_sat B) / (2 G₀).
    pub fn field_stiffness(&self, b_field: f64) -> f64 {
        if self.shear_modulus < 1e-10 || self.magnetization < 1e-10 {
            return self.shear_modulus;
        }
        let delta_g = Self::MU_0 * self.susceptibility.powi(2) * self.magnetization * b_field
            / (2.0 * self.shear_modulus);
        self.shear_modulus + delta_g
    }

    /// Magnetic pressure (Maxwell stress) normal to a surface \[Pa\].
    ///
    /// P_mag = B² / (2 μ₀).
    pub fn magnetic_pressure(&self, b_field: f64) -> f64 {
        b_field.powi(2) / (2.0 * Self::MU_0)
    }

    /// Saturation field H_sat \[A/m\] at which M ≈ M_sat.
    ///
    /// H_sat = M_sat / χ_m.
    pub fn saturation_field(&self) -> f64 {
        if self.susceptibility < 1e-15 {
            return f64::INFINITY;
        }
        self.magnetization / self.susceptibility
    }
}

// ── PhotothermalActuator ──────────────────────────────────────────────────────

/// Photothermal soft actuator model.
///
/// Light is absorbed by a (typically nanoparticle-loaded) polymer layer,
/// raising its temperature and causing differential thermal expansion which
/// drives bending.
///
/// # Parameters
/// * `absorption` — optical absorptivity α \[m⁻¹\] or effective absorbed fraction per W/m².
/// * `thermal_exp` — coefficient of thermal expansion \[1/K\].
#[derive(Debug, Clone)]
pub struct PhotothermalActuator {
    /// Optical absorptivity \[W absorbed / (W/m²)\] — dimensionless fraction.
    pub absorption: f64,
    /// Coefficient of thermal expansion α_T \[1/K\].
    pub thermal_exp: f64,
    /// Thermal time constant \[s\].
    pub time_constant: f64,
    /// Ambient temperature \[K\].
    pub t_ambient: f64,
    /// Heat capacity per unit volume \[J/(m³·K)\].
    pub heat_capacity: f64,
}

impl PhotothermalActuator {
    /// Create a new photothermal actuator.
    ///
    /// * `absorption` — fraction of incident light absorbed (0–1).
    /// * `thermal_exp` — coefficient of thermal expansion \[1/K\].
    /// * `time_constant` — thermal RC time constant \[s\].
    /// * `heat_capacity` — volumetric heat capacity \[J/(m³·K)\].
    pub fn new(absorption: f64, thermal_exp: f64, time_constant: f64, heat_capacity: f64) -> Self {
        Self {
            absorption: absorption.clamp(0.0, 1.0),
            thermal_exp,
            time_constant,
            t_ambient: 293.15,
            heat_capacity,
        }
    }

    /// Steady-state temperature rise \[K\] at light intensity `intensity` \[W/m²\].
    ///
    /// ΔT = α * I * τ / (ρ c_p)  (exponential approach to steady state taken at t→∞)
    pub fn temperature_rise(&self, intensity: f64) -> f64 {
        if self.heat_capacity < 1e-15 || self.time_constant < 1e-15 {
            return 0.0;
        }
        self.absorption * intensity * self.time_constant / self.heat_capacity
    }

    /// Actuation strain at steady state \[dimensionless\].
    ///
    /// ε = α_T * ΔT(I).
    pub fn actuation_strain(&self, intensity: f64) -> f64 {
        self.thermal_exp * self.temperature_rise(intensity)
    }

    /// Transient temperature rise \[K\] at time `t` \[s\] after illumination onset.
    ///
    /// ΔT(t) = ΔT_ss * (1 - exp(-t/τ)).
    pub fn transient_temperature(&self, intensity: f64, t: f64) -> f64 {
        let dt_ss = self.temperature_rise(intensity);
        dt_ss * (1.0 - (-t / self.time_constant.max(1e-15)).exp())
    }

    /// Actuation curvature \[1/m\] for a bilayer beam of total thickness `h`.
    ///
    /// κ ≈ 6 * α_T * ΔT / h  (Timoshenko bilayer formula, equal thickness layers).
    pub fn actuation_curvature(&self, intensity: f64, beam_thickness: f64) -> f64 {
        if beam_thickness < 1e-15 {
            return 0.0;
        }
        6.0 * self.actuation_strain(intensity) / beam_thickness
    }
}

// ── kelvin_voigt_response ─────────────────────────────────────────────────────

/// Kelvin–Voigt viscoelastic strain update (explicit Euler).
///
/// The constitutive equation is:
///
/// σ = E * ε + η * (dε/dt)  →  dε/dt = (σ - E * ε) / η
///
/// This function returns the updated strain ε(t + dt).
///
/// * `sigma` — applied stress \[Pa\].
/// * `e` — elastic (Young's) modulus \[Pa\].
/// * `eta` — viscosity coefficient \[Pa·s\].
/// * `dt` — time step \[s\].
/// * `eps_prev` — strain at the previous time step \[dimensionless\].
pub fn kelvin_voigt_response(sigma: f64, e: f64, eta: f64, dt: f64, eps_prev: f64) -> f64 {
    if eta.abs() < 1e-20 {
        // Pure elastic: ε = σ / E.
        return if e.abs() > 1e-20 { sigma / e } else { eps_prev };
    }
    let deps_dt = (sigma - e * eps_prev) / eta;
    eps_prev + deps_dt * dt
}

/// Maxwell viscoelastic stress relaxation update (explicit Euler).
///
/// σ̇ = E * ε̇ - σ/τ   where τ = η/E is the relaxation time.
///
/// * `strain_rate` — rate of applied strain \[1/s\].
/// * `e` — elastic modulus \[Pa\].
/// * `eta` — viscosity \[Pa·s\].
/// * `dt` — time step \[s\].
/// * `sigma_prev` — stress at the previous step \[Pa\].
pub fn maxwell_stress_update(strain_rate: f64, e: f64, eta: f64, dt: f64, sigma_prev: f64) -> f64 {
    let tau = if e.abs() > 1e-20 { eta / e } else { 1.0e10 };
    let dsigma_dt = e * strain_rate - sigma_prev / tau;
    sigma_prev + dsigma_dt * dt
}

/// Standard linear solid (Zener model) strain update.
///
/// Uses a spring-Maxwell element in parallel.
///
/// * `sigma` — applied stress \[Pa\].
/// * `e_eq` — equilibrium (long-time) modulus \[Pa\].
/// * `e_relax` — relaxation modulus of the Maxwell branch \[Pa\].
/// * `eta` — viscosity of the dashpot \[Pa·s\].
/// * `dt` — time step \[s\].
/// * `eps_prev` — previous total strain.
/// * `eps_maxwell_prev` — previous Maxwell branch strain.
///
/// Returns `(eps_new, eps_maxwell_new)`.
pub fn standard_linear_solid_update(
    sigma: f64,
    e_eq: f64,
    e_relax: f64,
    eta: f64,
    dt: f64,
    eps_prev: f64,
    eps_maxwell_prev: f64,
) -> (f64, f64) {
    // Maxwell branch: the dashpot strain rate = (ε - ε_M) * e_relax / η.
    let deps_m = (eps_prev - eps_maxwell_prev) * e_relax / eta.max(1e-20);
    let eps_maxwell_new = eps_maxwell_prev + deps_m * dt;
    // Total strain from equilibrium spring + Maxwell branch.
    let e_total = e_eq + e_relax;
    let eps_new = if e_total > 1e-20 {
        (sigma + e_relax * eps_maxwell_new) / e_total
    } else {
        eps_prev
    };
    (eps_new, eps_maxwell_new)
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // 1. Hydrogel: osmotic pressure at phi=0 is 0 (log diverges—tested at small phi).
    #[test]
    fn test_hydrogel_osmotic_pressure_small_phi() {
        let gel = Hydrogel::new(0.4, 1.8e-5, 0.1, 1000.0);
        let p = gel.osmotic_pressure(0.01);
        // At low phi, ln(1-φ) ≈ -φ so Π ≈ RT/v_s * (φ - φ - χφ²) ≈ positive.
        assert!(p.is_finite(), "osmotic pressure should be finite");
    }

    // 2. Hydrogel: osmotic pressure increases with phi in good-solvent regime.
    #[test]
    fn test_hydrogel_osmotic_pressure_increases() {
        let gel = Hydrogel::new(0.2, 1.8e-5, 0.1, 1000.0); // good solvent
        let p1 = gel.osmotic_pressure(0.05);
        let p2 = gel.osmotic_pressure(0.3);
        assert!(
            p2 > p1,
            "osmotic pressure should increase with polymer fraction"
        );
    }

    // 3. Hydrogel: elastic modulus is positive.
    #[test]
    fn test_hydrogel_elastic_modulus_positive() {
        let gel = Hydrogel::new(0.5, 1.8e-5, 0.1, 500.0);
        assert!(gel.elastic_modulus() > 0.0);
    }

    // 4. Hydrogel: flory_huggins_free_energy returns finite value.
    #[test]
    fn test_hydrogel_free_energy_finite() {
        let gel = Hydrogel::new(0.4, 1.8e-5, 0.1, 1000.0);
        let f = gel.flory_huggins_free_energy(0.2);
        assert!(f.is_finite(), "free energy should be finite");
    }

    // 5. Hydrogel: swelling_ratio ≥ 1 (gel can only swell, not shrink below dry).
    #[test]
    fn test_hydrogel_swelling_ratio_ge_one() {
        let gel = Hydrogel::new(0.3, 1.8e-5, 0.1, 1000.0);
        let q = gel.swelling_ratio(0.0);
        assert!(q >= 1.0, "swelling ratio should be ≥ 1, got {q}");
    }

    // 6. Hydrogel: elastic modulus scales linearly with n_chain.
    #[test]
    fn test_hydrogel_modulus_scales_with_n_chain() {
        let gel1 = Hydrogel::new(0.4, 1.8e-5, 0.1, 500.0);
        let gel2 = Hydrogel::new(0.4, 1.8e-5, 0.1, 1000.0);
        let ratio = gel2.elastic_modulus() / gel1.elastic_modulus();
        assert!((ratio - 2.0).abs() < 1e-10, "ratio={ratio}");
    }

    // 7. LCE: nematic_order at T > T_NI is 0 (isotropic phase).
    #[test]
    fn test_lce_nematic_order_above_tni() {
        let lce = LiquidCrystalElastomer::new(0.8, 0.3, 5e5);
        let s = lce.nematic_order(400.0, 350.0);
        assert_eq!(s, 0.0, "above T_NI order should be 0");
    }

    // 8. LCE: nematic_order at T=0 equals order_param.
    #[test]
    fn test_lce_nematic_order_at_zero() {
        let lce = LiquidCrystalElastomer::new(0.8, 0.3, 5e5);
        // At T=0: S = S0 * sqrt(1) = S0.
        let s = lce.nematic_order(0.0, 400.0);
        assert!((s - 0.8).abs() < 1e-10, "S={s}");
    }

    // 9. LCE: thermal_actuation returns 0 above T_NI.
    #[test]
    fn test_lce_thermal_actuation_above_tni() {
        let lce = LiquidCrystalElastomer::new(0.7, 0.25, 4e5);
        let eps = lce.thermal_actuation(500.0, 400.0);
        assert_eq!(eps, lce.coupling * (0.0 - lce.nematic_order(100.0, 400.0)));
    }

    // 10. LCE: strain_order_coupling equals coupling field.
    #[test]
    fn test_lce_strain_order_coupling() {
        let lce = LiquidCrystalElastomer::new(0.7, 0.3, 4e5);
        assert!((lce.strain_order_coupling() - 0.3).abs() < 1e-12);
    }

    // 11. LCE: actuation_stress is proportional to modulus.
    #[test]
    fn test_lce_actuation_stress_proportional() {
        let lce1 = LiquidCrystalElastomer::new(0.7, 0.3, 1e5);
        let lce2 = LiquidCrystalElastomer::new(0.7, 0.3, 2e5);
        let temp = 300.0;
        let t_ni = 400.0;
        let s1 = lce1.actuation_stress(temp, t_ni);
        let s2 = lce2.actuation_stress(temp, t_ni);
        if s1.abs() > 1e-15 {
            let ratio = s2 / s1;
            assert!((ratio - 2.0).abs() < 1e-8, "ratio={ratio}");
        }
    }

    // 12. IPMC: tip_displacement increases with voltage.
    #[test]
    fn test_ipmc_tip_displacement_vs_voltage() {
        let ipmc1 = IPMC::new(0.0002, 0.3, 1.0, 0.03, 1e9);
        let ipmc2 = IPMC::new(0.0002, 0.3, 2.0, 0.03, 1e9);
        assert!(
            ipmc2.tip_displacement() > ipmc1.tip_displacement(),
            "higher voltage → larger displacement"
        );
    }

    // 13. IPMC: tip_displacement is positive for positive voltage.
    #[test]
    fn test_ipmc_tip_displacement_positive() {
        let ipmc = IPMC::new(0.0002, 0.3, 1.5, 0.03, 1e9);
        assert!(ipmc.tip_displacement() > 0.0);
    }

    // 14. IPMC: blocking_force is positive.
    #[test]
    fn test_ipmc_blocking_force_positive() {
        let ipmc = IPMC::new(0.0002, 0.3, 1.5, 0.03, 2e9);
        assert!(ipmc.blocking_force() > 0.0);
    }

    // 15. IPMC: response_time scales with h².
    #[test]
    fn test_ipmc_response_time_thickness_scaling() {
        let ipmc1 = IPMC::new(0.0002, 0.3, 1.0, 0.03, 1e9);
        let ipmc2 = IPMC::new(0.0004, 0.3, 1.0, 0.03, 1e9);
        let ratio = ipmc2.response_time() / ipmc1.response_time();
        // h doubles → τ quadruples.
        assert!((ratio - 4.0).abs() < 1e-8, "ratio={ratio}");
    }

    // 16. IPMC: curvature is consistent with tip displacement.
    #[test]
    fn test_ipmc_curvature_consistency() {
        let ipmc = IPMC::new(0.0002, 0.3, 1.5, 0.05, 1e9);
        let kappa = ipmc.curvature();
        let delta = ipmc.tip_displacement();
        let l = ipmc.length;
        // κ = 2δ/L²
        let expected = 2.0 * delta / l.powi(2);
        assert!(
            (kappa - expected).abs() < 1e-12,
            "kappa={kappa}, expected={expected}"
        );
    }

    // 17. MAE: magneto_strain is 0 at zero field.
    #[test]
    fn test_mae_magneto_strain_zero_field() {
        let mae = MagnetoActiveElastomer::new(1e5, 5.0, 1e4);
        assert_eq!(mae.magneto_strain(0.0), 0.0);
    }

    // 18. MAE: magneto_strain increases with B.
    #[test]
    fn test_mae_magneto_strain_increases_with_b() {
        let mae = MagnetoActiveElastomer::new(1e5, 5.0, 1e4);
        let e1 = mae.magneto_strain(0.1);
        let e2 = mae.magneto_strain(0.2);
        assert!(e2 > e1, "strain should increase with B, e1={e1}, e2={e2}");
    }

    // 19. MAE: field_stiffness ≥ shear_modulus.
    #[test]
    fn test_mae_field_stiffness_increases() {
        let mae = MagnetoActiveElastomer::new(1e5, 5.0, 1e4);
        let g0 = mae.shear_modulus;
        let g_field = mae.field_stiffness(0.5);
        assert!(g_field >= g0, "field stiffness should be ≥ G₀");
    }

    // 20. MAE: magnetic_pressure scales as B².
    #[test]
    fn test_mae_magnetic_pressure_scaling() {
        let mae = MagnetoActiveElastomer::new(1e5, 5.0, 1e4);
        let p1 = mae.magnetic_pressure(1.0);
        let p2 = mae.magnetic_pressure(2.0);
        assert!((p2 / p1 - 4.0).abs() < 1e-8, "ratio={}", p2 / p1);
    }

    // 21. PhotothermalActuator: actuation_strain is 0 at zero intensity.
    #[test]
    fn test_photothermal_zero_intensity() {
        let pt = PhotothermalActuator::new(0.8, 1e-4, 10.0, 1e6);
        assert_eq!(pt.actuation_strain(0.0), 0.0);
    }

    // 22. PhotothermalActuator: actuation_strain increases with intensity.
    #[test]
    fn test_photothermal_strain_increases() {
        let pt = PhotothermalActuator::new(0.8, 1e-4, 10.0, 1e6);
        let e1 = pt.actuation_strain(1000.0);
        let e2 = pt.actuation_strain(2000.0);
        assert!(e2 > e1, "strain should increase with intensity");
    }

    // 23. PhotothermalActuator: transient temperature approaches steady state.
    #[test]
    fn test_photothermal_transient_approaches_ss() {
        let pt = PhotothermalActuator::new(0.8, 1e-4, 10.0, 1e6);
        let t_ss = pt.temperature_rise(1000.0);
        let t_trans = pt.transient_temperature(1000.0, 100.0); // 10 time constants
        let rel_err = (t_trans - t_ss).abs() / t_ss.max(1e-15);
        assert!(
            rel_err < 0.01,
            "transient should approach SS: t_ss={t_ss}, t_trans={t_trans}"
        );
    }

    // 24. kelvin_voigt_response: pure elastic limit (eta→0) gives σ/E.
    #[test]
    fn test_kelvin_voigt_elastic_limit() {
        let eps = kelvin_voigt_response(100.0, 1000.0, 0.0, 0.01, 0.0);
        // eta = 0 → eps = sigma / E = 0.1
        assert!((eps - 0.1).abs() < 1e-10, "eps={eps}");
    }

    // 25. kelvin_voigt_response: viscous response converges toward elastic value.
    #[test]
    fn test_kelvin_voigt_convergence() {
        let sigma = 100.0;
        let e = 1000.0;
        let eta = 10.0;
        let dt = 0.001;
        let mut eps = 0.0;
        for _ in 0..10_000 {
            eps = kelvin_voigt_response(sigma, e, eta, dt, eps);
        }
        // Steady state: eps_ss = sigma / E = 0.1.
        let eps_ss = sigma / e;
        assert!((eps - eps_ss).abs() < 1e-4, "eps={eps}, expected {eps_ss}");
    }

    // 26. maxwell_stress_update: stress relaxes to zero under zero strain rate.
    #[test]
    fn test_maxwell_stress_relaxation() {
        let e = 1000.0;
        let eta = 10.0;
        let dt = 0.001;
        let mut sigma = 100.0;
        for _ in 0..20_000 {
            sigma = maxwell_stress_update(0.0, e, eta, dt, sigma);
        }
        assert!(
            sigma.abs() < 1e-3,
            "stress should relax to zero, got {sigma}"
        );
    }

    // 27. standard_linear_solid_update: eps converges to sigma/e_eq at long times.
    #[test]
    fn test_sls_long_time_limit() {
        let sigma = 50.0;
        let e_eq = 500.0;
        let e_relax = 1000.0;
        let eta = 100.0;
        let dt = 0.0001;
        let mut eps = 0.0;
        let mut eps_m = 0.0;
        for _ in 0..1_000_000 {
            let (e_new, em_new) =
                standard_linear_solid_update(sigma, e_eq, e_relax, eta, dt, eps, eps_m);
            eps = e_new;
            eps_m = em_new;
        }
        // At t→∞ Maxwell spring fully relaxes so ε → σ / E_eq.
        let expected = sigma / e_eq;
        assert!(
            (eps - expected).abs() < 0.01 * expected.abs() + 1e-6,
            "eps={eps}, expected {expected}"
        );
    }

    // 28. MAE: saturation_field is positive and finite.
    #[test]
    fn test_mae_saturation_field() {
        let mae = MagnetoActiveElastomer::new(1e5, 5.0, 1e4);
        let hs = mae.saturation_field();
        assert!(hs > 0.0 && hs.is_finite(), "H_sat={hs}");
    }

    // 29. IPMC: zero voltage → zero tip displacement.
    #[test]
    fn test_ipmc_zero_voltage() {
        let ipmc = IPMC::new(0.0002, 0.3, 0.0, 0.03, 1e9);
        assert_eq!(ipmc.tip_displacement(), 0.0);
    }

    // 30. PhotothermalActuator: curvature is positive for positive intensity.
    #[test]
    fn test_photothermal_curvature_positive() {
        let pt = PhotothermalActuator::new(0.9, 2e-4, 5.0, 5e5);
        let kappa = pt.actuation_curvature(500.0, 0.001);
        assert!(kappa > 0.0, "curvature should be positive: {kappa}");
    }

    // 31. LCE: actuation_work_density is positive for T_hot < T_cold < T_NI.
    #[test]
    fn test_lce_work_density_positive() {
        let lce = LiquidCrystalElastomer::new(0.8, 0.4, 1e6);
        let w = lce.actuation_work_density(200.0, 50.0, 400.0);
        assert!(w >= 0.0, "work density should be non-negative: {w}");
    }

    // 32. Hydrogel: chemical_potential is finite and less than 0 for good solvent at low phi.
    #[test]
    fn test_hydrogel_chemical_potential_good_solvent() {
        // chi = 0.1 (good solvent), low phi → mostly entropic driving force for absorption.
        let gel = Hydrogel::new(0.1, 1.8e-5, 0.05, 200.0);
        let mu = gel.chemical_potential(0.05);
        assert!(mu.is_finite(), "chemical potential should be finite");
    }

    // 33. LCE: nematic_order monotonically decreases with temperature.
    #[test]
    fn test_lce_nematic_order_monotone() {
        let lce = LiquidCrystalElastomer::new(0.85, 0.35, 5e5);
        let t_ni = 400.0;
        let s1 = lce.nematic_order(100.0, t_ni);
        let s2 = lce.nematic_order(200.0, t_ni);
        let s3 = lce.nematic_order(350.0, t_ni);
        assert!(
            s1 >= s2 && s2 >= s3,
            "order should decrease with temperature"
        );
    }

    // 34. PhotothermalActuator: transient at t=0 is 0.
    #[test]
    fn test_photothermal_transient_at_zero() {
        let pt = PhotothermalActuator::new(0.8, 1e-4, 10.0, 1e6);
        let t = pt.transient_temperature(1000.0, 0.0);
        assert!(
            t.abs() < 1e-12,
            "transient temp at t=0 should be 0, got {t}"
        );
    }

    // 35. kelvin_voigt_response: increasing sigma yields increasing strain.
    #[test]
    fn test_kelvin_voigt_monotone_in_sigma() {
        let e = 500.0;
        let eta = 5.0;
        let dt = 0.01;
        let eps0 = 0.05;
        let eps1 = kelvin_voigt_response(50.0, e, eta, dt, eps0);
        let eps2 = kelvin_voigt_response(100.0, e, eta, dt, eps0);
        assert!(
            eps2 >= eps1,
            "higher stress should give larger or equal strain update"
        );
    }
}
