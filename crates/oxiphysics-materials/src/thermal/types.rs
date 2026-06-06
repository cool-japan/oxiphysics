//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::R_GAS;

/// Phase change (melting/solidification) model with latent heat.
///
/// Implements the enthalpy method for solving phase change problems.
/// The apparent heat capacity is increased in the mushy zone to account
/// for the latent heat release/absorption.
pub struct PhaseChangeModel {
    /// Liquidus temperature (melting point or end of mushy zone) \[K\].
    pub t_liquidus: f64,
    /// Solidus temperature (start of mushy zone) \[K\].
    pub t_solidus: f64,
    /// Latent heat of fusion L \[J/kg\].
    pub latent_heat: f64,
    /// Solid-phase specific heat \[J/(kg·K)\].
    pub cp_solid: f64,
    /// Liquid-phase specific heat \[J/(kg·K)\].
    pub cp_liquid: f64,
    /// Material density \[kg/m³\].
    pub density: f64,
}
impl PhaseChangeModel {
    /// Create a new phase change model.
    pub fn new(
        t_solidus: f64,
        t_liquidus: f64,
        latent_heat: f64,
        cp_solid: f64,
        cp_liquid: f64,
        density: f64,
    ) -> Self {
        assert!(t_liquidus >= t_solidus, "t_liquidus must be >= t_solidus");
        Self {
            t_liquidus,
            t_solidus,
            latent_heat,
            cp_solid,
            cp_liquid,
            density,
        }
    }
    /// Liquid fraction at temperature T (linear in mushy zone).
    pub fn liquid_fraction(&self, temp: f64) -> f64 {
        if temp <= self.t_solidus {
            0.0
        } else if temp >= self.t_liquidus {
            1.0
        } else {
            (temp - self.t_solidus) / (self.t_liquidus - self.t_solidus)
        }
    }
    /// Apparent (effective) specific heat in the enthalpy method \[J/(kg·K)\].
    ///
    /// In the mushy zone, cp_eff = cp + L/(T_liq - T_sol).
    pub fn apparent_specific_heat(&self, temp: f64) -> f64 {
        let fl = self.liquid_fraction(temp);
        let cp_mix = fl * self.cp_liquid + (1.0 - fl) * self.cp_solid;
        if temp > self.t_solidus && temp < self.t_liquidus {
            let dt = (self.t_liquidus - self.t_solidus).max(f64::EPSILON);
            cp_mix + self.latent_heat / dt
        } else {
            cp_mix
        }
    }
    /// Specific enthalpy h(T) relative to T = 0.
    pub fn enthalpy(&self, temp: f64) -> f64 {
        if temp <= self.t_solidus {
            self.cp_solid * temp
        } else if temp >= self.t_liquidus {
            self.cp_solid * self.t_solidus
                + self.latent_heat
                + self.cp_liquid * (temp - self.t_liquidus)
        } else {
            let fl = self.liquid_fraction(temp);
            self.cp_solid * self.t_solidus
                + fl * self.latent_heat
                + (fl * self.cp_liquid + (1.0 - fl) * self.cp_solid) * (temp - self.t_solidus)
        }
    }
    /// Thermal conductivity (interpolated between solid and liquid).
    pub fn conductivity(&self, k_solid: f64, k_liquid: f64, temp: f64) -> f64 {
        let fl = self.liquid_fraction(temp);
        fl * k_liquid + (1.0 - fl) * k_solid
    }
}
/// Coupled thermal-stress analysis for a 1D rod under thermal loading.
///
/// Combines 1D conduction with thermal stress calculation at each node.
pub struct ThermalStressCoupling {
    /// Number of nodes.
    pub n_nodes: usize,
    /// Node temperatures \[K\].
    pub temperature: Vec<f64>,
    /// Reference temperature (stress-free state) \[K\].
    pub t_ref: f64,
    /// Thermal expansion coefficient \[1/K\].
    pub alpha: f64,
    /// Young's modulus \[Pa\].
    pub young_modulus: f64,
    /// Poisson's ratio.
    pub poisson_ratio: f64,
}
impl ThermalStressCoupling {
    /// Create a new coupled analysis with uniform initial temperature.
    pub fn new(
        n_nodes: usize,
        t_init: f64,
        t_ref: f64,
        alpha: f64,
        young_modulus: f64,
        poisson_ratio: f64,
    ) -> Self {
        Self {
            n_nodes,
            temperature: vec![t_init; n_nodes],
            t_ref,
            alpha,
            young_modulus,
            poisson_ratio,
        }
    }
    /// Thermal strain at node i: ε_th = α·(T - T_ref).
    pub fn thermal_strain(&self, node: usize) -> f64 {
        self.alpha * (self.temperature[node] - self.t_ref)
    }
    /// Thermal stress at node i for fully constrained expansion \[Pa\].
    ///
    /// σ = -E·α·(T - T_ref) / (1 - ν)  (plane-stress)
    pub fn thermal_stress(&self, node: usize) -> f64 {
        -self.young_modulus * self.alpha * (self.temperature[node] - self.t_ref)
            / (1.0 - self.poisson_ratio)
    }
    /// Maximum thermal stress magnitude across all nodes \[Pa\].
    pub fn max_thermal_stress(&self) -> f64 {
        (0..self.n_nodes)
            .map(|i| self.thermal_stress(i).abs())
            .fold(0.0_f64, f64::max)
    }
    /// Check if any node exceeds the fracture stress \[Pa\].
    pub fn has_fracture(&self, fracture_stress: f64) -> bool {
        self.max_thermal_stress() > fracture_stress
    }
}
/// Heat-affected zone (HAZ) size model for welding processes.
///
/// Based on Rosenthal's solution for a moving heat source:
/// The HAZ boundary is defined where the peak temperature equals
/// a critical temperature (e.g., Ac1 for steel).
///
/// HAZ width: w ≈ (Q / (π·ρ·cp·v)) * (1/T_crit - 1/T_melt)  (thick plate)
#[derive(Debug, Clone)]
pub struct HeatAffectedZone {
    /// Heat input Q \[J/m\] (= Power / welding speed).
    pub heat_input: f64,
    /// Material density \[kg/m³\].
    pub rho: f64,
    /// Specific heat capacity \[J/(kg·K)\].
    pub cp: f64,
    /// Thermal conductivity \[W/(m·K)\].
    pub k: f64,
    /// Welding speed v \[m/s\].
    pub welding_speed: f64,
    /// Preheat / ambient temperature T_0 \[K\].
    pub t_preheat: f64,
    /// Melting temperature T_melt \[K\].
    pub t_melt: f64,
}
impl HeatAffectedZone {
    /// Create a new HAZ model.
    pub fn new(
        heat_input: f64,
        rho: f64,
        cp: f64,
        k: f64,
        welding_speed: f64,
        t_preheat: f64,
        t_melt: f64,
    ) -> Self {
        Self {
            heat_input,
            rho,
            cp,
            k,
            welding_speed,
            t_preheat,
            t_melt,
        }
    }
    /// Thermal diffusivity α = k / (ρ·cp) \[m²/s\].
    pub fn diffusivity(&self) -> f64 {
        self.k / (self.rho * self.cp)
    }
    /// Peak temperature at distance r from weld centre (Rosenthal, thin plate).
    ///
    /// T_peak = T_0 + Q / (π·ρ·cp·t_peak) * exp(-r²/(4·α·t_peak))
    ///
    /// Simplified thick-plate peak temperature at lateral distance y from weld line:
    /// T_peak ≈ T_0 + (Q / (2·π·k)) · (v / (2·π·α)) · exp(-v·y / (2·α))
    pub fn peak_temperature_thick_plate(&self, y: f64) -> f64 {
        let alpha = self.diffusivity();
        let coeff = self.heat_input * self.welding_speed / (2.0 * std::f64::consts::PI * alpha);
        self.t_preheat
            + coeff / (2.0 * std::f64::consts::PI) * (-self.welding_speed * y / (2.0 * alpha)).exp()
    }
    /// Estimated HAZ half-width \[m\] for a given critical temperature T_crit \[K\].
    ///
    /// Uses Rosenthal thick-plate formula: the HAZ extends to where
    /// T_peak = T_crit.
    pub fn haz_halfwidth(&self, t_crit: f64) -> f64 {
        let alpha = self.diffusivity();
        if t_crit <= self.t_preheat {
            return f64::INFINITY;
        }
        let dt = t_crit - self.t_preheat;
        self.heat_input / (std::f64::consts::E * std::f64::consts::PI * self.rho * self.cp * dt)
            * (alpha / self.welding_speed).sqrt()
    }
    /// Cooling rate at the weld centreline at T_crit \[K/s\].
    ///
    /// For thick plates: dT/dt = 2·π·k·(T_crit - T_0)² / Q
    pub fn cooling_rate_centreline(&self, t_crit: f64) -> f64 {
        2.0 * std::f64::consts::PI * self.k * (t_crit - self.t_preheat).powi(2) / self.heat_input
    }
    /// Check if the HAZ would cause martensite formation given a critical
    /// cooling rate \[K/s\] for martensite.
    pub fn forms_martensite(&self, critical_cooling_rate: f64, t_crit: f64) -> bool {
        self.cooling_rate_centreline(t_crit) >= critical_cooling_rate
    }
}
/// Johnson-Cook constitutive model for thermoviscoplastic flow stress.
///
/// σ = (A + B·ε^n) · (1 + C·ln(ε̇/ε̇₀)) · (1 - T*^m)
///
/// where T* = (T - T_room) / (T_melt - T_room) is the homologous temperature.
///
/// Reference: Johnson & Cook (1983).
#[derive(Debug, Clone)]
pub struct JohnsonCookModel {
    /// Yield stress A \[Pa\].
    pub a: f64,
    /// Strain hardening coefficient B \[Pa\].
    pub b: f64,
    /// Strain hardening exponent n (dimensionless).
    pub n: f64,
    /// Strain rate sensitivity C (dimensionless).
    pub c: f64,
    /// Thermal softening exponent m (dimensionless).
    pub m: f64,
    /// Reference strain rate ε̇₀ \[1/s\].
    pub eps_dot0: f64,
    /// Room temperature T_room \[K\].
    pub t_room: f64,
    /// Melt temperature T_melt \[K\].
    pub t_melt: f64,
}
impl JohnsonCookModel {
    /// Create a new Johnson-Cook model.
    pub fn new(
        a: f64,
        b: f64,
        n: f64,
        c: f64,
        m: f64,
        eps_dot0: f64,
        t_room: f64,
        t_melt: f64,
    ) -> Self {
        Self {
            a,
            b,
            n,
            c,
            m,
            eps_dot0,
            t_room,
            t_melt,
        }
    }
    /// Typical AISI 4340 steel parameters.
    pub fn aisi_4340() -> Self {
        Self::new(792e6, 510e6, 0.26, 0.014, 1.03, 1.0, 293.0, 1793.0)
    }
    /// Homologous temperature T* ∈ \[0, 1\].
    pub fn homologous_temperature(&self, temperature: f64) -> f64 {
        if temperature <= self.t_room {
            return 0.0;
        }
        if temperature >= self.t_melt {
            return 1.0;
        }
        (temperature - self.t_room) / (self.t_melt - self.t_room)
    }
    /// Flow stress σ \[Pa\].
    ///
    /// # Arguments
    /// * `eps`     - Equivalent plastic strain (dimensionless, ≥ 0)
    /// * `eps_dot` - Equivalent plastic strain rate \[1/s\]
    /// * `temp`    - Current temperature \[K\]
    pub fn flow_stress(&self, eps: f64, eps_dot: f64, temp: f64) -> f64 {
        let term_a = self.a + self.b * eps.max(0.0).powf(self.n);
        let eps_dot_ratio = (eps_dot / self.eps_dot0).max(1.0);
        let term_b = 1.0 + self.c * eps_dot_ratio.ln();
        let t_star = self.homologous_temperature(temp);
        let term_c = 1.0 - t_star.powf(self.m);
        term_a * term_b * term_c
    }
    /// Isothermal flow stress (room temperature, reference strain rate).
    pub fn isothermal_flow_stress(&self, eps: f64) -> f64 {
        self.flow_stress(eps, self.eps_dot0, self.t_room)
    }
    /// Strain-rate sensitivity ratio: σ(ε̇) / σ(ε̇₀) at given strain and temperature.
    pub fn strain_rate_sensitivity(&self, eps_dot: f64, eps: f64, temp: f64) -> f64 {
        let sig_ref = self.flow_stress(eps, self.eps_dot0, temp);
        if sig_ref < f64::EPSILON {
            return 1.0;
        }
        self.flow_stress(eps, eps_dot, temp) / sig_ref
    }
}
/// Thermal interface (Kapitza) resistance between two materials.
///
/// Heat flux across interface: q = delta_T / R_k
/// where R_k is the thermal interface resistance (m^2*K/W).
pub struct ThermalInterfaceResistance {
    /// Interface resistance R_k (m^2*K/W).
    pub resistance: f64,
}
impl ThermalInterfaceResistance {
    /// Create a new thermal interface resistance.
    pub fn new(resistance: f64) -> Self {
        Self { resistance }
    }
    /// Typical metal-metal interface: R ~ 1e-4 m^2*K/W.
    pub fn metal_metal() -> Self {
        Self::new(1e-4)
    }
    /// Typical metal-ceramic interface: R ~ 1e-3 m^2*K/W.
    pub fn metal_ceramic() -> Self {
        Self::new(1e-3)
    }
    /// Typical semiconductor interface: R ~ 1e-8 m^2*K/W.
    pub fn semiconductor() -> Self {
        Self::new(1e-8)
    }
    /// Heat flux across the interface for a given temperature drop.
    pub fn heat_flux(&self, delta_t: f64) -> f64 {
        delta_t / self.resistance
    }
    /// Temperature drop for a given heat flux.
    pub fn temperature_drop(&self, heat_flux: f64) -> f64 {
        heat_flux * self.resistance
    }
    /// Effective conductance (1/R_k) in W/(m^2*K).
    pub fn conductance(&self) -> f64 {
        1.0 / self.resistance
    }
}
/// 3×3 thermal conductivity tensor for anisotropic materials.
///
/// The heat flux vector q = -K · ∇T where K is the conductivity tensor.
/// For orthotropic materials K is diagonal in the principal axes.
pub struct ThermalConductivityTensor {
    /// Row-major 3×3 conductivity tensor \[W/(m·K)\].
    pub k: [[f64; 3]; 3],
}
impl ThermalConductivityTensor {
    /// Create from a full 3×3 matrix.
    pub fn new(k: [[f64; 3]; 3]) -> Self {
        Self { k }
    }
    /// Create an isotropic tensor: K = k · I.
    pub fn isotropic(k: f64) -> Self {
        Self {
            k: [[k, 0.0, 0.0], [0.0, k, 0.0], [0.0, 0.0, k]],
        }
    }
    /// Create a diagonal (orthotropic) tensor from principal conductivities.
    pub fn orthotropic(kx: f64, ky: f64, kz: f64) -> Self {
        Self {
            k: [[kx, 0.0, 0.0], [0.0, ky, 0.0], [0.0, 0.0, kz]],
        }
    }
    /// Compute heat flux vector q = -K · ∇T.
    ///
    /// # Arguments
    /// * `grad_t` — temperature gradient ∇T \[K/m\]
    ///
    /// Returns heat flux \[W/m²\].
    pub fn heat_flux(&self, grad_t: [f64; 3]) -> [f64; 3] {
        let mut q = [0.0_f64; 3];
        for (q_i, k_row) in q.iter_mut().zip(self.k.iter()) {
            for (k_ij, &gt_j) in k_row.iter().zip(grad_t.iter()) {
                *q_i -= k_ij * gt_j;
            }
        }
        q
    }
    /// Effective (geometric mean) conductivity for isotropic comparison.
    pub fn effective_conductivity(&self) -> f64 {
        (self.k[0][0] * self.k[1][1] * self.k[2][2]).powf(1.0 / 3.0)
    }
    /// Check if the tensor is symmetric (required for physical validity).
    pub fn is_symmetric(&self, tol: f64) -> bool {
        for i in 0..3 {
            for j in 0..3 {
                if (self.k[i][j] - self.k[j][i]).abs() > tol {
                    return false;
                }
            }
        }
        true
    }
}
/// Ablation model for surface recession under extreme heat flux.
///
/// Models the recession of a material surface when the surface temperature
/// exceeds the ablation (sublimation or vaporization) threshold.
/// The recession rate is governed by the Arrhenius-like expression.
pub struct AblationModel {
    /// Ablation threshold temperature T_ab \[K\].
    pub t_ablation: f64,
    /// Heat of ablation (enthalpy of vaporization) L_ab \[J/kg\].
    pub heat_of_ablation: f64,
    /// Material density \[kg/m³\].
    pub density: f64,
    /// Pre-exponential factor A \[m/s\].
    pub arrhenius_a: f64,
    /// Activation temperature T_act = Ea/(R) \[K\].
    pub activation_temp: f64,
}
impl AblationModel {
    /// Create a new ablation model.
    pub fn new(
        t_ablation: f64,
        heat_of_ablation: f64,
        density: f64,
        arrhenius_a: f64,
        activation_temp: f64,
    ) -> Self {
        Self {
            t_ablation,
            heat_of_ablation,
            density,
            arrhenius_a,
            activation_temp,
        }
    }
    /// Surface recession rate \[m/s\] at temperature T.
    ///
    /// r = A * exp(-T_act / T)  for T > T_ab, else 0
    pub fn recession_rate(&self, temp: f64) -> f64 {
        if temp <= self.t_ablation {
            return 0.0;
        }
        self.arrhenius_a * (-self.activation_temp / temp).exp()
    }
    /// Mass loss rate \[kg/(m²·s)\] = density * recession_rate.
    pub fn mass_loss_rate(&self, temp: f64) -> f64 {
        self.density * self.recession_rate(temp)
    }
    /// Blowing heat flux \[W/m²\] absorbed by ablation = L_ab * mdot.
    pub fn ablative_heat_flux(&self, temp: f64) -> f64 {
        self.heat_of_ablation * self.mass_loss_rate(temp)
    }
    /// Net surface heat flux after ablative cooling \[W/m²\].
    ///
    /// q_net = q_applied - q_ablation
    pub fn net_heat_flux(&self, applied_flux: f64, temp: f64) -> f64 {
        applied_flux - self.ablative_heat_flux(temp)
    }
}
/// Enthalpy method for melting and solidification problems.
///
/// The total enthalpy H at temperature T includes sensible and latent heat:
///
/// H(T) = ρ·cp·T + ρ·L·f_l(T)
///
/// where f_l is the liquid fraction. The temperature-update is performed
/// by inverting H(T) given the updated enthalpy.
#[derive(Debug, Clone)]
pub struct EnthalpyMethod {
    /// Density \[kg/m³\].
    pub rho: f64,
    /// Specific heat capacity \[J/(kg·K)\].
    pub cp: f64,
    /// Latent heat of fusion L \[J/kg\].
    pub latent_heat: f64,
    /// Solidus temperature T_s \[K\].
    pub t_solidus: f64,
    /// Liquidus temperature T_l \[K\].
    pub t_liquidus: f64,
    /// Current temperature \[K\].
    pub temperature: f64,
    /// Current enthalpy \[J/m³\].
    pub enthalpy: f64,
}
impl EnthalpyMethod {
    /// Create a new enthalpy method model.
    pub fn new(rho: f64, cp: f64, latent_heat: f64, t_solidus: f64, t_liquidus: f64) -> Self {
        let t0 = t_solidus - 1.0;
        let h0 = rho * cp * t0;
        Self {
            rho,
            cp,
            latent_heat,
            t_solidus,
            t_liquidus,
            temperature: t0,
            enthalpy: h0,
        }
    }
    /// Liquid fraction as a linear function of temperature in mushy zone.
    pub fn liquid_fraction(&self, temp: f64) -> f64 {
        if temp <= self.t_solidus {
            0.0
        } else if temp >= self.t_liquidus {
            1.0
        } else {
            (temp - self.t_solidus) / (self.t_liquidus - self.t_solidus)
        }
    }
    /// Total enthalpy H(T) \[J/m³\].
    pub fn enthalpy_at(&self, temp: f64) -> f64 {
        self.rho * self.cp * temp + self.rho * self.latent_heat * self.liquid_fraction(temp)
    }
    /// Invert enthalpy to temperature (Newton-Raphson, max 50 iterations).
    pub fn temperature_from_enthalpy(&self, h: f64) -> f64 {
        let mut t = h / (self.rho * self.cp);
        for _ in 0..50 {
            let h_t = self.enthalpy_at(t);
            let dh_dt = self.effective_cp(t);
            let dt = (h - h_t) / dh_dt.max(f64::EPSILON);
            t += dt;
            if dt.abs() < 1e-6 {
                break;
            }
        }
        t
    }
    /// Effective (apparent) specific heat cp_eff \[J/(kg·K)\] including latent heat.
    pub fn effective_cp(&self, temp: f64) -> f64 {
        if temp > self.t_solidus && temp < self.t_liquidus {
            let dfl_dt = 1.0 / (self.t_liquidus - self.t_solidus);
            self.rho * (self.cp + self.latent_heat * dfl_dt)
        } else {
            self.rho * self.cp
        }
    }
    /// Apply a heat flux q \[W/m³\] over time step dt \[s\].
    pub fn apply_heat_flux(&mut self, q: f64, dt: f64) {
        self.enthalpy += q * dt;
        self.temperature = self.temperature_from_enthalpy(self.enthalpy);
    }
    /// Check if material is currently in the mushy zone.
    pub fn is_mushy(&self) -> bool {
        self.temperature > self.t_solidus && self.temperature < self.t_liquidus
    }
}
/// Temperature-dependent specific heat (piecewise-linear tabulation).
pub struct TemperatureDependentSpecificHeat {
    /// Temperature sample points \[K\].
    pub temps: Vec<f64>,
    /// Specific heat values at each temperature \[J/(kg*K)\].
    pub values: Vec<f64>,
}
impl TemperatureDependentSpecificHeat {
    /// Create a new table.
    pub fn new(temps: Vec<f64>, values: Vec<f64>) -> Self {
        assert_eq!(temps.len(), values.len());
        Self { temps, values }
    }
    /// Linearly interpolate specific heat at `temperature`.
    pub fn cp_at(&self, temperature: f64) -> f64 {
        let n = self.temps.len();
        if n == 0 {
            return 0.0;
        }
        if n == 1 || temperature <= self.temps[0] {
            return self.values[0];
        }
        if temperature >= self.temps[n - 1] {
            return self.values[n - 1];
        }
        let idx = self
            .temps
            .partition_point(|&t| t <= temperature)
            .saturating_sub(1);
        let t0 = self.temps[idx];
        let t1 = self.temps[idx + 1];
        let c0 = self.values[idx];
        let c1 = self.values[idx + 1];
        let frac = (temperature - t0) / (t1 - t0);
        c0 + frac * (c1 - c0)
    }
}
/// Thermal shock resistance parameters.
///
/// R  = sigma_f * (1 - nu) / (E * alpha)        -- maximum delta_T without fracture
/// R' = R * k                                     -- with heat conduction
/// R'' = sigma_f^2 * (1 - nu) / (E * alpha^2)    -- thermal stress damage resistance
pub struct ThermalShockResistance {
    /// Fracture strength sigma_f (Pa).
    pub fracture_strength: f64,
    /// Young's modulus E (Pa).
    pub young_modulus: f64,
    /// Poisson's ratio.
    pub poisson_ratio: f64,
    /// Coefficient of thermal expansion (1/K).
    pub thermal_expansion: f64,
    /// Thermal conductivity (W/(m*K)).
    pub thermal_conductivity: f64,
}
impl ThermalShockResistance {
    /// Create a new ThermalShockResistance calculator.
    pub fn new(
        fracture_strength: f64,
        young_modulus: f64,
        poisson_ratio: f64,
        thermal_expansion: f64,
        thermal_conductivity: f64,
    ) -> Self {
        Self {
            fracture_strength,
            young_modulus,
            poisson_ratio,
            thermal_expansion,
            thermal_conductivity,
        }
    }
    /// First thermal shock resistance parameter R (K).
    ///
    /// Maximum temperature change the material can withstand.
    pub fn r_parameter(&self) -> f64 {
        self.fracture_strength * (1.0 - self.poisson_ratio)
            / (self.young_modulus * self.thermal_expansion)
    }
    /// Second parameter R' = R * k (W/m).
    pub fn r_prime(&self) -> f64 {
        self.r_parameter() * self.thermal_conductivity
    }
    /// Third parameter R'' = sigma_f^2 * (1-nu) / (E * alpha^2) (Pa*K).
    pub fn r_double_prime(&self) -> f64 {
        self.fracture_strength * self.fracture_strength * (1.0 - self.poisson_ratio)
            / (self.young_modulus * self.thermal_expansion * self.thermal_expansion)
    }
    /// Check if a given temperature change will cause thermal shock failure.
    pub fn will_fail(&self, delta_t: f64) -> bool {
        delta_t.abs() > self.r_parameter()
    }
}
/// Temperature-dependent thermal conductivity (piecewise-linear tabulation).
pub struct TemperatureDependentConductivity {
    /// Temperature sample points \[K\].
    pub temps: Vec<f64>,
    /// Conductivity values at each temperature \[W/(m*K)\].
    pub conductivities: Vec<f64>,
}
impl TemperatureDependentConductivity {
    /// Create a new table.
    pub fn new(temps: Vec<f64>, conductivities: Vec<f64>) -> Self {
        assert_eq!(temps.len(), conductivities.len());
        Self {
            temps,
            conductivities,
        }
    }
    /// Linearly interpolate conductivity at `temperature`.
    pub fn conductivity_at(&self, temperature: f64) -> f64 {
        let n = self.temps.len();
        if n == 0 {
            return 0.0;
        }
        if n == 1 || temperature <= self.temps[0] {
            return self.conductivities[0];
        }
        if temperature >= self.temps[n - 1] {
            return self.conductivities[n - 1];
        }
        let idx = self
            .temps
            .partition_point(|&t| t <= temperature)
            .saturating_sub(1);
        let t0 = self.temps[idx];
        let t1 = self.temps[idx + 1];
        let k0 = self.conductivities[idx];
        let k1 = self.conductivities[idx + 1];
        let frac = (temperature - t0) / (t1 - t0);
        k0 + frac * (k1 - k0)
    }
    /// Average conductivity over \[t_low, t_high\] using the trapezoidal rule.
    pub fn mean_conductivity(&self, t_low: f64, t_high: f64) -> f64 {
        if (t_high - t_low).abs() < f64::EPSILON {
            return self.conductivity_at(t_low);
        }
        let mut pts: Vec<f64> = vec![t_low];
        for &t in &self.temps {
            if t > t_low && t < t_high {
                pts.push(t);
            }
        }
        pts.push(t_high);
        pts.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        pts.dedup_by(|a, b| (*a - *b).abs() < f64::EPSILON);
        let mut integral = 0.0;
        for i in 0..pts.len() - 1 {
            let ta = pts[i];
            let tb = pts[i + 1];
            let ka = self.conductivity_at(ta);
            let kb = self.conductivity_at(tb);
            integral += 0.5 * (ka + kb) * (tb - ta);
        }
        integral / (t_high - t_low)
    }
}
/// Debye specific heat model.
///
/// Cv = 9*N*k_B * (T/Theta_D)^3 * integral_0^{Theta_D/T} x^4*e^x/(e^x - 1)^2 dx
///
/// In the high-T limit: Cv -> 3*N*k_B (Dulong-Petit).
/// In the low-T limit: Cv ~ T^3.
pub struct DebyeModel {
    /// Debye temperature Theta_D (K).
    pub theta_d: f64,
    /// Number of atoms per formula unit (for molar heat capacity).
    pub n_atoms: f64,
}
impl DebyeModel {
    /// Create a new Debye model.
    pub fn new(theta_d: f64, n_atoms: f64) -> Self {
        Self { theta_d, n_atoms }
    }
    /// Molar heat capacity Cv (J/(mol*K)) at temperature T.
    ///
    /// Uses numerical integration with Simpson's rule.
    pub fn molar_cv(&self, temperature: f64) -> f64 {
        if temperature < 1e-10 {
            return 0.0;
        }
        let x_max = self.theta_d / temperature;
        let n_steps = 200;
        let dx = x_max / (n_steps as f64);
        let integrand = |x: f64| -> f64 {
            if x < 1e-12 {
                return 0.0;
            }
            if x > 500.0 {
                return 0.0;
            }
            let ex = x.exp();
            x.powi(4) * ex / (ex - 1.0).powi(2)
        };
        let mut sum = integrand(0.0) + integrand(x_max);
        for i in 1..n_steps {
            let x = i as f64 * dx;
            let weight = if i % 2 == 0 { 2.0 } else { 4.0 };
            sum += weight * integrand(x);
        }
        let integral = sum * dx / 3.0;
        let ratio = temperature / self.theta_d;
        9.0 * self.n_atoms * R_GAS * ratio.powi(3) * integral
    }
    /// Dulong-Petit limit: Cv = 3*n*R.
    pub fn dulong_petit(&self) -> f64 {
        3.0 * self.n_atoms * R_GAS
    }
}
/// Thermal expansion model computing volumetric and anisotropic strain.
///
/// For isotropic materials the volumetric thermal strain is simply
/// ΔV/V = 3 α ΔT, where α is the linear thermal expansion coefficient.
/// For orthotropic/anisotropic materials separate coefficients are used
/// along each principal axis.
#[derive(Debug, Clone)]
pub struct ThermalExpansion {
    /// Linear thermal expansion coefficient along x-axis \[1/K\].
    pub alpha_x: f64,
    /// Linear thermal expansion coefficient along y-axis \[1/K\].
    pub alpha_y: f64,
    /// Linear thermal expansion coefficient along z-axis \[1/K\].
    pub alpha_z: f64,
    /// Reference temperature T_ref \[K\] at which the body is stress-free.
    pub t_ref: f64,
}
impl ThermalExpansion {
    /// Create an **isotropic** thermal expansion model.
    ///
    /// # Arguments
    /// * `alpha` - Isotropic linear CTE \[1/K\]
    /// * `t_ref` - Stress-free reference temperature \[K\]
    pub fn isotropic(alpha: f64, t_ref: f64) -> Self {
        Self {
            alpha_x: alpha,
            alpha_y: alpha,
            alpha_z: alpha,
            t_ref,
        }
    }
    /// Create an **orthotropic** thermal expansion model.
    ///
    /// # Arguments
    /// * `alpha_x`, `alpha_y`, `alpha_z` - Principal CTEs \[1/K\]
    /// * `t_ref` - Stress-free reference temperature \[K\]
    pub fn orthotropic(alpha_x: f64, alpha_y: f64, alpha_z: f64, t_ref: f64) -> Self {
        Self {
            alpha_x,
            alpha_y,
            alpha_z,
            t_ref,
        }
    }
    /// Linear strain along x: ε_x = αₓ · ΔT.
    pub fn strain_x(&self, temperature: f64) -> f64 {
        self.alpha_x * (temperature - self.t_ref)
    }
    /// Linear strain along y: ε_y = αᵧ · ΔT.
    pub fn strain_y(&self, temperature: f64) -> f64 {
        self.alpha_y * (temperature - self.t_ref)
    }
    /// Linear strain along z: ε_z = α_z · ΔT.
    pub fn strain_z(&self, temperature: f64) -> f64 {
        self.alpha_z * (temperature - self.t_ref)
    }
    /// Compute the volumetric thermal strain ΔV/V = ε_x + ε_y + ε_z.
    ///
    /// For isotropic materials: ΔV/V = 3 α ΔT (exact for small strains).
    /// For orthotropic materials: ΔV/V = (α_x + α_y + α_z) ΔT.
    ///
    /// The exact large-deformation volumetric strain is:
    /// ΔV/V = (1 + ε_x)(1 + ε_y)(1 + ε_z) − 1
    /// but the linearised form is used here as it is sufficient for most
    /// engineering temperature ranges (|αΔT| << 1).
    ///
    /// # Arguments
    /// * `temperature` - Current temperature T \[K\]
    ///
    /// # Returns
    /// Volumetric strain ΔV/V \[dimensionless\]
    pub fn compute_volumetric_strain(&self, temperature: f64) -> f64 {
        let delta_t = temperature - self.t_ref;
        (self.alpha_x + self.alpha_y + self.alpha_z) * delta_t
    }
    /// Exact (large-deformation) volumetric strain J − 1.
    ///
    /// J = (1 + ε_x)(1 + ε_y)(1 + ε_z)
    pub fn compute_volumetric_strain_exact(&self, temperature: f64) -> f64 {
        let ex = 1.0 + self.strain_x(temperature);
        let ey = 1.0 + self.strain_y(temperature);
        let ez = 1.0 + self.strain_z(temperature);
        ex * ey * ez - 1.0
    }
    /// Thermal dilatation tensor as \[ε_x, ε_y, ε_z\].
    pub fn strain_tensor(&self, temperature: f64) -> [f64; 3] {
        [
            self.strain_x(temperature),
            self.strain_y(temperature),
            self.strain_z(temperature),
        ]
    }
}
/// Extended thermal shock resistance parameter calculations.
///
/// Beyond the basic R = k·σ_f / (E·α) parameter (already in ThermalShockResistance),
/// this struct adds:
/// - R'' (energy-based, including fracture toughness)
/// - Hasselman parameter R''''
/// - Thermal shock damage resistance
#[derive(Debug, Clone)]
pub struct ThermalShockParam {
    /// Fracture strength \[Pa\].
    pub sigma_f: f64,
    /// Young's modulus \[Pa\].
    pub young_modulus: f64,
    /// Poisson's ratio.
    pub nu: f64,
    /// Thermal expansion coefficient \[1/K\].
    pub alpha: f64,
    /// Thermal conductivity \[W/(m·K)\].
    pub k: f64,
    /// Fracture toughness K_Ic \[Pa·√m\].
    pub k_ic: f64,
}
impl ThermalShockParam {
    /// Create a new thermal shock parameter set.
    pub fn new(sigma_f: f64, young_modulus: f64, nu: f64, alpha: f64, k: f64, k_ic: f64) -> Self {
        Self {
            sigma_f,
            young_modulus,
            nu,
            alpha,
            k,
            k_ic,
        }
    }
    /// First thermal shock resistance parameter R \[K\]:
    /// R = σ_f · (1-ν) / (E·α)
    pub fn r_first(&self) -> f64 {
        self.sigma_f * (1.0 - self.nu) / (self.young_modulus * self.alpha)
    }
    /// Second parameter R' \[W/m\]:
    /// R' = k · σ_f · (1-ν) / (E·α)
    pub fn r_second(&self) -> f64 {
        self.k * self.r_first()
    }
    /// Third parameter R'' \[J/m\]:
    /// R'' = G_f / (E·α²) where G_f = K_Ic² / E
    pub fn r_third(&self) -> f64 {
        let g_f = self.k_ic * self.k_ic / self.young_modulus;
        g_f / (self.young_modulus * self.alpha * self.alpha)
    }
    /// Hasselman unified thermal shock parameter \[K\]:
    /// R'''' = E · G_f / (σ_f² · (1-ν))
    pub fn r_hasselman(&self) -> f64 {
        let g_f = self.k_ic * self.k_ic / self.young_modulus;
        self.young_modulus * g_f / (self.sigma_f * self.sigma_f * (1.0 - self.nu))
    }
}
/// Einstein specific heat model.
///
/// Cv = 3*N*R * (Theta_E/T)^2 * exp(Theta_E/T) / (exp(Theta_E/T) - 1)^2
pub struct EinsteinModel {
    /// Einstein temperature Theta_E (K).
    pub theta_e: f64,
    /// Number of atoms per formula unit.
    pub n_atoms: f64,
}
impl EinsteinModel {
    /// Create a new Einstein model.
    pub fn new(theta_e: f64, n_atoms: f64) -> Self {
        Self { theta_e, n_atoms }
    }
    /// Molar heat capacity Cv (J/(mol*K)) at temperature T.
    pub fn molar_cv(&self, temperature: f64) -> f64 {
        if temperature < 1e-10 {
            return 0.0;
        }
        let x = self.theta_e / temperature;
        if x > 500.0 {
            return 0.0;
        }
        let ex = x.exp();
        3.0 * self.n_atoms * R_GAS * x * x * ex / (ex - 1.0).powi(2)
    }
}
/// Fourier heat conduction in 1D via explicit finite differences.
pub struct HeatConduction1D {
    /// Number of nodes (including boundary nodes).
    pub n_nodes: usize,
    /// Domain length \[m\].
    pub length: f64,
    /// Current nodal temperatures \[K\].
    pub temperature: Vec<f64>,
    /// Material properties.
    pub material: ThermalMaterial,
}
impl HeatConduction1D {
    /// Construct a uniform-temperature rod.
    pub fn new(n_nodes: usize, length: f64, t_init: f64, material: ThermalMaterial) -> Self {
        assert!(n_nodes >= 2, "need at least 2 nodes");
        Self {
            n_nodes,
            length,
            temperature: vec![t_init; n_nodes],
            material,
        }
    }
    /// Set Dirichlet BCs at both ends.
    pub fn set_temperature_bc(&mut self, t_left: f64, t_right: f64) {
        self.temperature[0] = t_left;
        self.temperature[self.n_nodes - 1] = t_right;
    }
    /// Explicit (forward-Euler) time step.
    pub fn step_explicit(&mut self, dt: f64) -> bool {
        let dx = self.length / (self.n_nodes - 1) as f64;
        let alpha = self.material.diffusivity();
        let r = alpha * dt / (dx * dx);
        if r > 0.5 {
            return false;
        }
        let old = self.temperature.clone();
        for i in 1..self.n_nodes - 1 {
            self.temperature[i] = old[i] + r * (old[i + 1] - 2.0 * old[i] + old[i - 1]);
        }
        true
    }
    /// Explicit step with volumetric heat source Q (W/m^3).
    pub fn step_explicit_with_source(&mut self, dt: f64, source: f64) -> bool {
        let dx = self.length / (self.n_nodes - 1) as f64;
        let alpha = self.material.diffusivity();
        let r = alpha * dt / (dx * dx);
        if r > 0.5 {
            return false;
        }
        let old = self.temperature.clone();
        let rho_cp = self.material.volumetric_heat_capacity();
        for i in 1..self.n_nodes - 1 {
            self.temperature[i] =
                old[i] + r * (old[i + 1] - 2.0 * old[i] + old[i - 1]) + source * dt / rho_cp;
        }
        true
    }
    /// Maximum stable time step: dx^2 / (2 * alpha).
    pub fn critical_dt(&self) -> f64 {
        let dx = self.length / (self.n_nodes - 1) as f64;
        let alpha = self.material.diffusivity();
        (dx * dx) / (2.0 * alpha)
    }
    /// Exact steady-state temperature profile (linear between BCs).
    pub fn steady_state_temperature(&self, t_left: f64, t_right: f64) -> Vec<f64> {
        (0..self.n_nodes)
            .map(|i| {
                let frac = i as f64 / (self.n_nodes - 1) as f64;
                t_left + frac * (t_right - t_left)
            })
            .collect()
    }
    /// Heat flux at each interior node via central finite differences.
    pub fn heat_flux(&self) -> Vec<f64> {
        let dx = self.length / (self.n_nodes - 1) as f64;
        let k = self.material.thermal_conductivity;
        (1..self.n_nodes - 1)
            .map(|i| {
                let dt_dx = (self.temperature[i + 1] - self.temperature[i - 1]) / (2.0 * dx);
                -k * dt_dx
            })
            .collect()
    }
    /// Total heat stored in the rod (integral of rho*cp*T*dx).
    pub fn total_heat_content(&self) -> f64 {
        let dx = self.length / (self.n_nodes - 1) as f64;
        let rho_cp = self.material.volumetric_heat_capacity();
        self.temperature.iter().map(|&t| rho_cp * t * dx).sum()
    }
}
/// Isotropic thermal material.
pub struct ThermalMaterial {
    /// Thermal conductivity \[W/(m*K)\]
    pub thermal_conductivity: f64,
    /// Specific heat capacity \[J/(kg*K)\]
    pub specific_heat: f64,
    /// Density \[kg/m^3\]
    pub density: f64,
    /// Linear coefficient of thermal expansion \[1/K\]
    pub thermal_expansion: f64,
}
impl ThermalMaterial {
    /// Create a new `ThermalMaterial`.
    pub fn new(k: f64, cp: f64, rho: f64, alpha: f64) -> Self {
        Self {
            thermal_conductivity: k,
            specific_heat: cp,
            density: rho,
            thermal_expansion: alpha,
        }
    }
    /// Structural steel: k=50, cp=500, rho=7850, alpha=12e-6.
    pub fn steel() -> Self {
        Self::new(50.0, 500.0, 7850.0, 12e-6)
    }
    /// Aluminium alloy: k=237, cp=900, rho=2700, alpha=23e-6.
    pub fn aluminum() -> Self {
        Self::new(237.0, 900.0, 2700.0, 23e-6)
    }
    /// Copper: k=401, cp=385, rho=8960, alpha=17e-6.
    pub fn copper() -> Self {
        Self::new(401.0, 385.0, 8960.0, 17e-6)
    }
    /// Concrete: k=1.5, cp=880, rho=2300, alpha=10e-6.
    pub fn concrete() -> Self {
        Self::new(1.5, 880.0, 2300.0, 10e-6)
    }
    /// Silicon carbide (SiC): k=120, cp=750, rho=3210, alpha=4.0e-6.
    pub fn silicon_carbide() -> Self {
        Self::new(120.0, 750.0, 3210.0, 4.0e-6)
    }
    /// Alumina (Al2O3): k=30, cp=880, rho=3950, alpha=8.0e-6.
    pub fn alumina() -> Self {
        Self::new(30.0, 880.0, 3950.0, 8.0e-6)
    }
    /// Thermal diffusivity alpha = k / (rho * cp)  \[m^2/s\].
    pub fn diffusivity(&self) -> f64 {
        self.thermal_conductivity / (self.density * self.specific_heat)
    }
    /// Thermal effusivity e = sqrt(k * rho * cp)  \[J/(m^2*K*s^0.5)\].
    pub fn effusivity(&self) -> f64 {
        (self.thermal_conductivity * self.density * self.specific_heat).sqrt()
    }
    /// Free thermal expansion strain for temperature change delta_t.
    pub fn thermal_strain(&self, delta_t: f64) -> f64 {
        self.thermal_expansion * delta_t
    }
    /// Thermal stress for fully constrained expansion (plane-stress).
    ///
    /// sigma = -E * alpha * delta_T / (1 - nu)
    pub fn thermal_stress(&self, delta_t: f64, young_modulus: f64, poisson_ratio: f64) -> f64 {
        -young_modulus * self.thermal_expansion * delta_t / (1.0 - poisson_ratio)
    }
    /// Volumetric heat capacity rho * cp \[J/(m^3*K)\].
    pub fn volumetric_heat_capacity(&self) -> f64 {
        self.density * self.specific_heat
    }
    /// Thermal penetration depth at time t.
    ///
    /// delta = sqrt(4 * alpha * t)
    pub fn penetration_depth(&self, t: f64) -> f64 {
        (4.0 * self.diffusivity() * t).sqrt()
    }
    /// Compute the first thermal shock resistance parameter R \[K\].
    ///
    /// R = σ_f · (1 − ν) / (E · α)
    ///
    /// R is the maximum temperature change the material can sustain without
    /// fracture under fully constrained thermal expansion. Higher R indicates
    /// better thermal shock resistance.
    ///
    /// # Arguments
    /// * `fracture_strength` - Tensile fracture strength σ_f \[Pa\]
    /// * `young_modulus`     - Young's modulus E \[Pa\]
    /// * `poisson_ratio`     - Poisson's ratio ν
    ///
    /// # Returns
    /// Thermal shock resistance R \[K\]
    pub fn compute_thermal_shock_resistance(
        &self,
        fracture_strength: f64,
        young_modulus: f64,
        poisson_ratio: f64,
    ) -> f64 {
        fracture_strength * (1.0 - poisson_ratio) / (young_modulus * self.thermal_expansion)
    }
    /// Estimate the electronic contribution to thermal conductivity via
    /// the Wiedemann-Franz law.
    ///
    /// κ_el = L · T / ρ_e
    ///
    /// where L = 2.44 × 10⁻⁸ W·Ω/K² is the Lorenz number and ρ_e is the
    /// electrical resistivity. This gives the electron-mediated thermal
    /// conductivity, which dominates in metals.
    ///
    /// # Arguments
    /// * `temperature`           - Absolute temperature T \[K\]
    /// * `electrical_resistivity` - Electrical resistivity ρ_e \[Ω·m\]
    ///
    /// # Returns
    /// Electronic thermal conductivity κ_el \[W/(m·K)\]
    pub fn compute_wiedemann_franz(&self, temperature: f64, electrical_resistivity: f64) -> f64 {
        const LORENZ: f64 = 2.4427e-8_f64;
        LORENZ * temperature / electrical_resistivity
    }
}
/// Newtonian cooling -- lumped capacitance model.
pub struct NewtonianCooling {
    /// Mass \[kg\].
    pub mass: f64,
    /// Specific heat \[J/(kg*K)\].
    pub specific_heat: f64,
    /// Convective heat-transfer coefficient \[W/(m^2*K)\].
    pub heat_transfer_coeff: f64,
    /// Surface area \[m^2\].
    pub surface_area: f64,
    /// Current temperature \[K\].
    pub temperature: f64,
}
impl NewtonianCooling {
    /// Construct a new lumped-capacitance object.
    pub fn new(mass: f64, cp: f64, h: f64, area: f64, t_init: f64) -> Self {
        Self {
            mass,
            specific_heat: cp,
            heat_transfer_coeff: h,
            surface_area: area,
            temperature: t_init,
        }
    }
    /// Time constant tau = m * cp / (h * A) \[s\].
    pub fn time_constant(&self) -> f64 {
        self.mass * self.specific_heat / (self.heat_transfer_coeff * self.surface_area)
    }
    /// Temperature at time `time` given ambient `t_inf`.
    pub fn temperature_at_time(&self, t_inf: f64, time: f64) -> f64 {
        let tau = self.time_constant();
        t_inf + (self.temperature - t_inf) * (-time / tau).exp()
    }
    /// Advance the temperature by one Euler step of size `dt`.
    pub fn step(&mut self, t_inf: f64, dt: f64) {
        let tau = self.time_constant();
        self.temperature += -dt / tau * (self.temperature - t_inf);
    }
    /// Time needed to reach `t_target` given ambient `t_inf`.
    pub fn time_to_temperature(&self, t_inf: f64, t_target: f64) -> Option<f64> {
        let tau = self.time_constant();
        let denom = self.temperature - t_inf;
        if denom.abs() < f64::EPSILON {
            return None;
        }
        let ratio = (t_target - t_inf) / denom;
        if ratio <= 0.0 {
            return None;
        }
        Some(-tau * ratio.ln())
    }
    /// Biot number: Bi = h * L / k (characteristic length L).
    ///
    /// For a sphere, L = radius/3.
    /// Lumped capacitance is valid when Bi < 0.1.
    pub fn biot_number(&self, characteristic_length: f64, thermal_conductivity: f64) -> f64 {
        self.heat_transfer_coeff * characteristic_length / thermal_conductivity
    }
}
/// Thermal fatigue model using the Coffin-Manson relationship
/// with temperature-range dependency.
///
/// The plastic strain range Δε_p depends on the temperature range ΔT:
/// Δε_p = α · ΔT · E_correction
///
/// Coffin-Manson: N_f = C / Δε_p^β
#[derive(Debug, Clone)]
pub struct ThermalFatigue {
    /// Coffin-Manson ductility coefficient C.
    pub c_cm: f64,
    /// Coffin-Manson exponent β (typically 0.5..0.7).
    pub beta: f64,
    /// Thermal expansion coefficient \[1/K\].
    pub alpha: f64,
    /// Young's modulus \[Pa\].
    pub young_modulus: f64,
    /// Accumulated Miner's damage D.
    pub damage: f64,
}
impl ThermalFatigue {
    /// Create a new thermal fatigue model.
    pub fn new(c_cm: f64, beta: f64, alpha: f64, young_modulus: f64) -> Self {
        Self {
            c_cm,
            beta,
            alpha,
            young_modulus,
            damage: 0.0,
        }
    }
    /// Plastic strain range for a given temperature range ΔT \[K\].
    ///
    /// Δε_p = α · ΔT  (simplified, assuming full inelastic constraint)
    pub fn plastic_strain_range(&self, delta_t: f64) -> f64 {
        self.alpha * delta_t.abs()
    }
    /// Fatigue life N_f for a given temperature range ΔT.
    ///
    /// N_f = C / Δε_p^β
    pub fn fatigue_life(&self, delta_t: f64) -> f64 {
        let deps = self.plastic_strain_range(delta_t);
        if deps < f64::EPSILON {
            return f64::INFINITY;
        }
        self.c_cm / deps.powf(self.beta)
    }
    /// Accumulate damage for `n_cycles` at temperature range ΔT.
    pub fn accumulate(&mut self, delta_t: f64, n_cycles: f64) {
        let n_f = self.fatigue_life(delta_t);
        if n_f.is_finite() {
            self.damage = (self.damage + n_cycles / n_f).min(1.0);
        }
    }
    /// Check if thermal fatigue failure has occurred (D ≥ 1).
    pub fn is_failed(&self) -> bool {
        self.damage >= 1.0 - f64::EPSILON
    }
    /// Temperature range causing failure in N cycles.
    ///
    /// Invert: ΔT = Δε_p / α,  Δε_p = (C/N)^(1/β)
    pub fn critical_delta_t(&self, n_cycles: f64) -> f64 {
        if n_cycles <= 0.0 {
            return 0.0;
        }
        let deps = (self.c_cm / n_cycles).powf(1.0 / self.beta);
        deps / self.alpha
    }
}
