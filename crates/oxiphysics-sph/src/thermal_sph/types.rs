//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use std::f64::consts::PI;

use super::functions::*;
use crate::thermal_sph::types_ext::*;

/// Adiabatic compression heating model.
///
/// Accounts for the temperature increase when a fluid parcel is
/// compressed, relevant for deep-ocean and atmospheric SPH.
#[derive(Clone, Debug)]
pub struct AdiabaticCompression {
    /// Ratio of specific heats γ (Cp/Cv).
    pub gamma: f64,
    /// Reference pressure p_0 (Pa).
    pub p_ref: f64,
    /// Reference temperature T_0 (K).
    pub t_ref: f64,
}
impl AdiabaticCompression {
    /// Create an adiabatic compression model.
    pub fn new(gamma: f64, p_ref: f64, t_ref: f64) -> Self {
        Self {
            gamma,
            p_ref,
            t_ref,
        }
    }
    /// Adiabatic temperature at pressure p: T = T_0 (p/p_0)^((γ-1)/γ).
    pub fn adiabatic_temperature(&self, pressure: f64) -> f64 {
        if self.p_ref < 1e-300 {
            return self.t_ref;
        }
        let exponent = (self.gamma - 1.0) / self.gamma;
        self.t_ref * (pressure / self.p_ref).powf(exponent)
    }
    /// Temperature rate from adiabatic compression: dT/dt = (γ-1)/γ * T/p * dp/dt.
    pub fn dtemp_dt_adiabatic(&self, temp: f64, pressure: f64, dp_dt: f64) -> f64 {
        if pressure < 1e-300 {
            return 0.0;
        }
        (self.gamma - 1.0) / self.gamma * temp / pressure * dp_dt
    }
    /// Potential temperature (removes adiabatic effect): θ = T (p_0/p)^((γ-1)/γ).
    pub fn potential_temperature(&self, temp: f64, pressure: f64) -> f64 {
        if self.p_ref < 1e-300 || pressure < 1e-300 {
            return temp;
        }
        let exponent = (self.gamma - 1.0) / self.gamma;
        temp * (self.p_ref / pressure).powf(exponent)
    }
    /// Brunt-Väisälä frequency squared N² = g/θ * dθ/dz.
    pub fn brunt_vaisala_sq(&self, dtheta_dz: f64, theta: f64, g: f64) -> f64 {
        if theta.abs() < 1e-300 {
            return 0.0;
        }
        g / theta * dtheta_dz
    }
}
/// Boiling heat transfer in SPH with bubble nucleation.
///
/// Implements pool boiling correlations, bubble growth, and departure,
/// using the Jakob number to characterize the phase-change regime.
pub struct BoilingSph {
    /// Fluid saturation temperature T_sat (K).
    pub t_sat: f64,
    /// Latent heat of vaporization L (J/kg).
    pub latent_heat: f64,
    /// Liquid density ρ_l (kg/m³).
    pub rho_l: f64,
    /// Vapor density ρ_v (kg/m³).
    pub rho_v: f64,
    /// Liquid specific heat c_p (J/kg/K).
    pub cp_l: f64,
    /// Surface tension σ (N/m).
    pub surface_tension: f64,
    /// Nucleation site density N_s (sites/m²).
    pub nucleation_density: f64,
}
impl BoilingSph {
    /// Construct a boiling model.
    pub fn new(
        t_sat: f64,
        latent_heat: f64,
        rho_l: f64,
        rho_v: f64,
        cp_l: f64,
        surface_tension: f64,
        nucleation_density: f64,
    ) -> Self {
        BoilingSph {
            t_sat,
            latent_heat,
            rho_l,
            rho_v,
            cp_l,
            surface_tension,
            nucleation_density,
        }
    }
    /// Jakob number Ja = ρ_l c_p ΔT_superheat / (ρ_v L).
    pub fn jakob_number(&self, t_wall: f64) -> f64 {
        let dt_sup = (t_wall - self.t_sat).max(0.0);
        let denom = self.rho_v * self.latent_heat;
        if denom < 1e-300 {
            return 0.0;
        }
        self.rho_l * self.cp_l * dt_sup / denom
    }
    /// Critical bubble radius for nucleation: r* = 2σ / (ρ_v L ΔT / T_sat).
    pub fn critical_bubble_radius(&self, t_wall: f64) -> f64 {
        let dt_sup = (t_wall - self.t_sat).max(0.0);
        if dt_sup < 1e-300 || self.rho_v < 1e-300 || self.latent_heat < 1e-300 {
            return f64::INFINITY;
        }
        2.0 * self.surface_tension * self.t_sat / (self.rho_v * self.latent_heat * dt_sup)
    }
    /// Bubble departure frequency (Fritz correlation, Hz).
    ///
    /// f = 0.149 * sqrt(g (ρ_l - ρ_v) / σ) * g^0.5 * d_b^0.5
    pub fn departure_frequency(&self, bubble_diameter: f64) -> f64 {
        let g = 9.81_f64;
        if bubble_diameter < 1e-300 || self.surface_tension < 1e-300 {
            return 0.0;
        }
        let inner = g * (self.rho_l - self.rho_v) / self.surface_tension;
        if inner < 0.0 {
            return 0.0;
        }
        0.149 * inner.sqrt()
    }
    /// Rohsenow pool boiling heat flux correlation (W/m²).
    ///
    /// q = μ_l h_fg \[g(ρ_l-ρ_v)/σ\]^(1/2) \[c_pl ΔT_sup/(C_sf h_fg Pr^1.7)\]^3
    ///
    /// Simplified version using empirical coefficient.
    pub fn rohsenow_heat_flux(&self, t_wall: f64, mu_l: f64, prandtl: f64) -> f64 {
        let dt_sup = (t_wall - self.t_sat).max(0.0);
        let g = 9.81_f64;
        let h_fg = self.latent_heat;
        let csf = 0.013_f64;
        let n = 1.0_f64;
        let sigma = self.surface_tension;
        if sigma < 1e-300 || h_fg < 1e-300 || prandtl < 1e-300 {
            return 0.0;
        }
        let term1 = mu_l * h_fg;
        let term2 = (g * (self.rho_l - self.rho_v).max(0.0) / sigma).sqrt();
        let term3 = (self.cp_l * dt_sup / (csf * h_fg * prandtl.powf(n))).powi(3);
        term1 * term2 * term3
    }
    /// Critical heat flux (Zuber correlation) (W/m²).
    pub fn critical_heat_flux(&self) -> f64 {
        let g = 9.81_f64;
        let pi_4 = PI / 4.0;
        let sigma = self.surface_tension;
        let h_fg = self.latent_heat;
        if sigma < 1e-300 {
            return 0.0;
        }
        let term = (sigma * g * (self.rho_l - self.rho_v).max(0.0)).powf(0.25);
        pi_4 * h_fg * self.rho_v * term / self.rho_l.max(1e-300).sqrt()
    }
    /// Leidenfrost temperature (min film boiling temperature) estimate (K).
    ///
    /// Uses Berenson's correlation: T_Leid ≈ T_sat + 0.127 h_fg ρ_v / (μ_v ...).
    /// Simplified: T_Leid ≈ 1.15 T_sat (approximate).
    pub fn leidenfrost_temperature(&self) -> f64 {
        1.15 * self.t_sat
    }
}
/// Complete thermal SPH simulation state manager.
///
/// Manages a collection of thermal particles and advances them
/// using the Brookshaw heat conduction operator with optional
/// phase change, radiation, and boiling.
pub struct ThermalSphSimulation {
    /// All particles in the simulation.
    pub particles: Vec<ThermalParticle>,
    /// SPH smoothing length (m).
    pub h: f64,
    /// Current simulation time (s).
    pub time: f64,
    /// Time step (s).
    pub dt: f64,
    /// Enable phase-change model.
    pub enable_phase_change: bool,
    /// Optional phase-change parameters.
    pub phase_change: Option<PhaseChangeSph>,
}
impl ThermalSphSimulation {
    /// Construct a new simulation.
    pub fn new(h: f64, dt: f64) -> Self {
        ThermalSphSimulation {
            particles: Vec::new(),
            h,
            time: 0.0,
            dt,
            enable_phase_change: false,
            phase_change: None,
        }
    }
    /// Add a particle to the simulation.
    pub fn add_particle(&mut self, p: ThermalParticle) {
        self.particles.push(p);
    }
    /// Enable phase change with given parameters.
    pub fn with_phase_change(mut self, pc: PhaseChangeSph) -> Self {
        self.phase_change = Some(pc);
        self.enable_phase_change = true;
        self
    }
    /// Number of particles.
    pub fn particle_count(&self) -> usize {
        self.particles.len()
    }
    /// Mean temperature of all particles.
    pub fn mean_temperature(&self) -> f64 {
        if self.particles.is_empty() {
            return 0.0;
        }
        self.particles.iter().map(|p| p.temperature).sum::<f64>() / self.particles.len() as f64
    }
    /// Maximum temperature in the simulation.
    pub fn max_temperature(&self) -> f64 {
        self.particles
            .iter()
            .map(|p| p.temperature)
            .fold(f64::NEG_INFINITY, f64::max)
    }
    /// Minimum temperature in the simulation.
    pub fn min_temperature(&self) -> f64 {
        self.particles
            .iter()
            .map(|p| p.temperature)
            .fold(f64::INFINITY, f64::min)
    }
    /// Total thermal energy Σ m c_p T.
    pub fn total_thermal_energy(&self) -> f64 {
        self.particles
            .iter()
            .map(|p| p.mass * p.specific_heat * p.temperature)
            .sum()
    }
    /// Advance all particles one time step using brute-force O(N²) neighbor search.
    pub fn step(&mut self) {
        let n = self.particles.len();
        if n == 0 {
            return;
        }
        let mut d_temps = vec![0.0f64; n];
        for (i, dt_i) in d_temps.iter_mut().enumerate().take(n) {
            let mut neighbors = Vec::new();
            let pos_i = self.particles[i].pos;
            let t_i = self.particles[i].temperature;
            let h = self.h;
            for j in 0..n {
                if i == j {
                    continue;
                }
                let r_ij = sub3(pos_i, self.particles[j].pos);
                if len3(r_ij) > 2.0 * h {
                    continue;
                }
                neighbors.push((
                    r_ij,
                    self.particles[j].temperature,
                    self.particles[j].mass,
                    self.particles[j].density,
                ));
            }
            let alpha = self.particles[i].thermal_diffusivity();
            let hc = SphHeatConduction::new(alpha, h);
            *dt_i = hc.dtemp_dt(t_i, &neighbors);
        }
        for (p, &dt_rate) in self.particles.iter_mut().zip(d_temps.iter()) {
            p.temperature += self.dt * dt_rate;
        }
        self.time += self.dt;
    }
}
/// Thermo-elastic stress from SPH temperature field.
///
/// σ_thermal = -E α ΔT / (1 - 2ν) (isotropic, Voigt notation).
pub struct ThermalStress {
    /// Young's modulus E (Pa).
    pub young_modulus: f64,
    /// Thermal expansion coefficient α (1/K).
    pub alpha: f64,
    /// Poisson's ratio ν.
    pub poisson: f64,
    /// Reference temperature T_ref (K).
    pub t_ref: f64,
}
impl ThermalStress {
    /// Construct a thermo-elastic stress model.
    pub fn new(young_modulus: f64, alpha: f64, poisson: f64, t_ref: f64) -> Self {
        ThermalStress {
            young_modulus,
            alpha,
            poisson,
            t_ref,
        }
    }
    /// Isotropic thermal stress magnitude σ_th = E α ΔT / (1 - 2ν).
    pub fn thermal_stress_magnitude(&self, temp: f64) -> f64 {
        let dt = temp - self.t_ref;
        let denom = 1.0 - 2.0 * self.poisson;
        if denom.abs() < 1e-300 {
            return 0.0;
        }
        self.young_modulus * self.alpha * dt / denom
    }
    /// Voigt thermal stress vector \[σ_xx, σ_yy, σ_zz, 0, 0, 0\].
    pub fn thermal_stress_voigt(&self, temp: f64) -> [f64; 6] {
        let s = self.thermal_stress_magnitude(temp);
        [s, s, s, 0.0, 0.0, 0.0]
    }
    /// Von Mises equivalent thermal stress: σ_VM = E α |ΔT| * √(2/3) * √3.
    pub fn von_mises_thermal(&self, temp: f64) -> f64 {
        let s = self.thermal_stress_magnitude(temp).abs();
        let _ = s;
        0.0
    }
}
/// Semi-analytic Stefan problem solution for validation.
///
/// Provides the Neumann similarity solution for 1D solidification
/// of a half-space initially at melting temperature.
/// Interface position: s(t) = 2 β √(α_s t)
/// where β is the root of: β exp(β²) erf(β) = Ste / √π
/// and Ste = c_p (T_m - T_w) / L is the Stefan number.
#[derive(Clone, Debug)]
pub struct StefanProblemAnalytic {
    /// Wall temperature T_w (K) — kept below T_melt.
    pub t_wall: f64,
    /// Melting temperature T_melt (K).
    pub t_melt: f64,
    /// Thermal diffusivity of solid α_s (m²/s).
    pub alpha_solid: f64,
    /// Thermal conductivity of solid λ_s (W/m/K).
    pub lambda_solid: f64,
    /// Specific heat of solid c_p (J/kg/K).
    pub cp_solid: f64,
    /// Latent heat L (J/kg).
    pub latent_heat: f64,
    /// Density ρ (kg/m³).
    pub density: f64,
}
impl StefanProblemAnalytic {
    /// Create a new Stefan problem analytic solver.
    pub fn new(
        t_wall: f64,
        t_melt: f64,
        alpha_solid: f64,
        lambda_solid: f64,
        cp_solid: f64,
        latent_heat: f64,
        density: f64,
    ) -> Self {
        Self {
            t_wall,
            t_melt,
            alpha_solid,
            lambda_solid,
            cp_solid,
            latent_heat,
            density,
        }
    }
    /// Stefan number Ste = c_p (T_m - T_w) / L.
    pub fn stefan_number(&self) -> f64 {
        if self.latent_heat < 1e-300 {
            return f64::INFINITY;
        }
        self.cp_solid * (self.t_melt - self.t_wall) / self.latent_heat
    }
    /// Similarity parameter β (iterative Newton solve for β exp(β²) erf(β) = Ste/√π).
    pub fn similarity_parameter(&self) -> f64 {
        let ste = self.stefan_number();
        let target = ste / PI.sqrt();
        let mut beta = (ste / 2.0).sqrt().max(0.01);
        for _ in 0..30 {
            let f = beta * (beta * beta).exp() * erf_approx(beta) - target;
            let df = (beta * beta).exp() * erf_approx(beta)
                + beta * 2.0 * beta * (beta * beta).exp() * erf_approx(beta)
                + beta * (beta * beta).exp() * 2.0 / PI.sqrt() * (-(beta * beta)).exp();
            if df.abs() < 1e-300 {
                break;
            }
            let step = f / df;
            beta -= step;
            beta = beta.max(1e-10);
            if step.abs() < 1e-10 {
                break;
            }
        }
        beta
    }
    /// Interface position s(t) = 2 β √(α_s t).
    pub fn interface_position(&self, time: f64) -> f64 {
        if time < 0.0 {
            return 0.0;
        }
        let beta = self.similarity_parameter();
        2.0 * beta * (self.alpha_solid * time).sqrt()
    }
    /// Temperature profile in solid at position x and time t.
    pub fn temperature_profile(&self, x: f64, time: f64) -> f64 {
        if time < 1e-300 {
            return self.t_wall;
        }
        let beta = self.similarity_parameter();
        let xi = x / (2.0 * (self.alpha_solid * time).sqrt());
        let erf_b = erf_approx(beta);
        if erf_b.abs() < 1e-300 {
            return self.t_melt;
        }
        self.t_wall + (self.t_melt - self.t_wall) * erf_approx(xi) / erf_b
    }
}
/// Full Boussinesq approximation for SPH thermal buoyancy.
///
/// Handles variable-temperature buoyancy-driven flows, including
/// non-linear equation of state corrections and Oberbeck-Boussinesq limit.
#[derive(Clone, Debug)]
pub struct BoussinesqSph {
    /// Reference density ρ_0 (kg/m³).
    pub rho_ref: f64,
    /// Reference temperature T_0 (K).
    pub t_ref: f64,
    /// Thermal expansion coefficient β (1/K).
    pub beta: f64,
    /// Gravity vector \[gx, gy, gz\] (m/s²).
    pub gravity: [f64; 3],
    /// Enable non-linear correction (second-order β_2 term).
    pub nonlinear: bool,
    /// Second-order thermal expansion coefficient β_2 (1/K²).
    pub beta2: f64,
}
impl BoussinesqSph {
    /// Create a linear Boussinesq model.
    pub fn new(rho_ref: f64, t_ref: f64, beta: f64, gravity: [f64; 3]) -> Self {
        Self {
            rho_ref,
            t_ref,
            beta,
            gravity,
            nonlinear: false,
            beta2: 0.0,
        }
    }
    /// Create a non-linear Boussinesq model.
    pub fn nonlinear(rho_ref: f64, t_ref: f64, beta: f64, beta2: f64, gravity: [f64; 3]) -> Self {
        Self {
            rho_ref,
            t_ref,
            beta,
            gravity,
            nonlinear: true,
            beta2,
        }
    }
    /// Density variation: ρ = ρ_0 \[1 - β ΔT - β_2 ΔT²\].
    pub fn density(&self, temp: f64) -> f64 {
        let dt = temp - self.t_ref;
        if self.nonlinear {
            self.rho_ref * (1.0 - self.beta * dt - self.beta2 * dt * dt)
        } else {
            self.rho_ref * (1.0 - self.beta * dt)
        }
    }
    /// Buoyancy force per unit mass f_b = -g β ΔT (Boussinesq approx).
    pub fn buoyancy_force(&self, temp: f64) -> [f64; 3] {
        let dt = temp - self.t_ref;
        let coeff = if self.nonlinear {
            -self.beta * dt - self.beta2 * dt * dt
        } else {
            -self.beta * dt
        };
        [
            self.gravity[0] * coeff,
            self.gravity[1] * coeff,
            self.gravity[2] * coeff,
        ]
    }
    /// Rayleigh number Ra = g β ΔT H³ / (ν κ).
    pub fn rayleigh_number(&self, delta_t: f64, height: f64, nu: f64, kappa: f64) -> f64 {
        if nu < 1e-300 || kappa < 1e-300 {
            return 0.0;
        }
        let g_mag = len3(self.gravity);
        g_mag * self.beta * delta_t.abs() * height.powi(3) / (nu * kappa)
    }
    /// Convective instability criterion: Ra > Ra_crit = 1708 (horizontal layer).
    pub fn is_convectively_unstable(&self, delta_t: f64, height: f64, nu: f64, kappa: f64) -> bool {
        self.rayleigh_number(delta_t, height, nu, kappa) > 1708.0
    }
}
/// Mushy-zone momentum sink based on Carman-Kozeny porosity model.
///
/// Used in casting and welding simulations to model flow resistance
/// in the partially solidified mushy zone.
pub struct MushyZoneDrag {
    /// Carman-Kozeny constant C_k (kg/m³/s).
    pub carman_kozeny: f64,
    /// Small constant to prevent division by zero.
    pub epsilon: f64,
}
impl MushyZoneDrag {
    /// Construct a mushy-zone drag model.
    pub fn new(carman_kozeny: f64) -> Self {
        MushyZoneDrag {
            carman_kozeny,
            epsilon: 1e-6,
        }
    }
    /// Darcy drag force per unit volume: F_d = -C_k (1-f_l)² / (f_l³ + ε) * v.
    ///
    /// `liquid_fraction` f_l ∈ \[0, 1\], `velocity` fluid velocity.
    pub fn drag_force(&self, liquid_fraction: f64, velocity: [f64; 3]) -> [f64; 3] {
        let fl = liquid_fraction.clamp(0.0, 1.0);
        let numerator = (1.0 - fl).powi(2);
        let denominator = fl.powi(3) + self.epsilon;
        let coeff = -self.carman_kozeny * numerator / denominator;
        scale3(velocity, coeff)
    }
    /// Dimensionless drag coefficient D = (1-f_l)² / (f_l³ + ε).
    pub fn drag_coefficient(&self, liquid_fraction: f64) -> f64 {
        let fl = liquid_fraction.clamp(0.0, 1.0);
        (1.0 - fl).powi(2) / (fl.powi(3) + self.epsilon)
    }
    /// Critical liquid fraction below which flow is essentially blocked.
    pub fn coherency_fraction(&self) -> f64 {
        0.3
    }
}
/// SPH internal energy equation integrator.
///
/// Solves du/dt = -p/ρ² ∇·v + α∇²T + Q_chem + Q_rad where u is
/// specific internal energy.
pub struct InternalEnergyEq {
    /// Ratio of specific heats γ (for ideal gas).
    pub gamma: f64,
    /// Chemical heat source Q_chem (W/kg).
    pub q_chem: f64,
    /// Radiative heat source Q_rad (W/kg).
    pub q_rad: f64,
}
impl InternalEnergyEq {
    /// Construct an internal energy equation integrator.
    pub fn new(gamma: f64) -> Self {
        InternalEnergyEq {
            gamma,
            q_chem: 0.0,
            q_rad: 0.0,
        }
    }
    /// Set chemical heat source Q_chem (W/kg).
    pub fn with_chemical_source(mut self, q_chem: f64) -> Self {
        self.q_chem = q_chem;
        self
    }
    /// Set radiative heat source Q_rad (W/kg).
    pub fn with_radiation_source(mut self, q_rad: f64) -> Self {
        self.q_rad = q_rad;
        self
    }
    /// Pressure work term: -p/ρ² (∇·v)_sph (J/kg/s).
    pub fn pressure_work(&self, pressure: f64, density: f64, div_v: f64) -> f64 {
        if density < 1e-300 {
            return 0.0;
        }
        -pressure / (density * density) * div_v
    }
    /// SPH velocity divergence approximation.
    ///
    /// (∇·v)_i = Σ_j m_j/ρ_j (v_j - v_i) · ∇W_ij
    pub fn velocity_divergence(
        v_i: [f64; 3],
        neighbors: &[([f64; 3], [f64; 3], f64, f64)],
        h: f64,
    ) -> f64 {
        let mut div = 0.0;
        for &(r_ij, v_j, m_j, rho_j) in neighbors {
            if rho_j < 1e-300 {
                continue;
            }
            let dv = sub3(v_j, v_i);
            let gw = cubic_kernel_grad(r_ij, h);
            div += m_j / rho_j * dot3(dv, gw);
        }
        div
    }
    /// Full du/dt for a fluid particle.
    pub fn du_dt(&self, pressure: f64, density: f64, div_v: f64, conduction_rate: f64) -> f64 {
        self.pressure_work(pressure, density, div_v) + conduction_rate + self.q_chem + self.q_rad
    }
    /// Ideal gas temperature from specific internal energy: T = (γ-1) u / R.
    pub fn temperature_from_energy(&self, u: f64, r_spec: f64) -> f64 {
        if r_spec < 1e-300 {
            return 0.0;
        }
        (self.gamma - 1.0) * u / r_spec
    }
}
/// Anisotropic SPH heat conduction using full conductivity tensor.
///
/// Implements the Cleary (1998) formulation for anisotropic media.
pub struct AnisotropicHeatConduction {
    /// Conductivity tensor (W/m/K).
    pub tensor: ConductivityTensor,
    /// SPH smoothing length (m).
    pub h: f64,
}
impl AnisotropicHeatConduction {
    /// Construct with an isotropic conductivity.
    pub fn new_isotropic(lambda: f64, h: f64) -> Self {
        AnisotropicHeatConduction {
            tensor: ConductivityTensor::isotropic(lambda),
            h,
        }
    }
    /// Construct with a full conductivity tensor.
    pub fn new_tensor(tensor: ConductivityTensor, h: f64) -> Self {
        AnisotropicHeatConduction { tensor, h }
    }
    /// Compute anisotropic heat conduction rate dT/dt at particle i.
    ///
    /// Uses: dT/dt = Σ_j (m_j/ρ_j) * (λ_ij · r_ij)/|r_ij|² * (T_i - T_j) * (∇W · r_ij)
    pub fn dtemp_dt(&self, t_i: f64, rho_i: f64, neighbors: &[([f64; 3], f64, f64, f64)]) -> f64 {
        if rho_i < 1e-300 {
            return 0.0;
        }
        let mut result = 0.0;
        for &(r_ij, t_j, m_j, rho_j) in neighbors {
            let r2 = dot3(r_ij, r_ij);
            if r2 < 1e-300 || rho_j < 1e-300 {
                continue;
            }
            let grad_w = cubic_kernel_grad(r_ij, self.h);
            let lambda_r = self.tensor.apply(r_ij);
            let lambda_dot = dot3(lambda_r, grad_w);
            result += m_j / rho_j * 2.0 * (t_i - t_j) / r2 * lambda_dot;
        }
        result / rho_i
    }
}
/// Radiation in the optically thin limit for SPH.
///
/// In the optically thin approximation (τ << 1), the medium emits but
/// does not absorb radiation from other particles. Each particle radiates
/// as a blackbody into empty space.
#[derive(Clone, Debug)]
pub struct OpticallyThinRadiation {
    /// Stefan-Boltzmann constant σ (W/m²/K⁴).
    pub sigma_sb: f64,
    /// Emissivity ε ∈ \[0, 1\].
    pub emissivity: f64,
    /// Ambient radiation temperature T_rad (K).
    pub t_ambient: f64,
    /// Rosseland mean opacity κ_R (m²/kg, mass-based).
    pub opacity: f64,
}
impl OpticallyThinRadiation {
    /// Create an optically thin radiation model.
    pub fn new(emissivity: f64, t_ambient: f64, opacity: f64) -> Self {
        Self {
            sigma_sb: 5.670_374_4e-8,
            emissivity,
            t_ambient,
            opacity,
        }
    }
    /// Volumetric cooling rate per unit mass Q_rad (W/kg).
    ///
    /// Q = -4 κ σ (T⁴ - T_rad⁴) — emission minus absorbed ambient.
    pub fn volumetric_cooling(&self, temp: f64, density: f64) -> f64 {
        if density < 1e-300 {
            return 0.0;
        }
        let emission = 4.0 * self.opacity * self.sigma_sb * temp.powi(4);
        let absorption = 4.0 * self.opacity * self.sigma_sb * self.t_ambient.powi(4);
        -(emission - absorption)
    }
    /// Net radiative power emitted per unit area (W/m²).
    pub fn net_emission(&self, temp: f64) -> f64 {
        self.emissivity * self.sigma_sb * (temp.powi(4) - self.t_ambient.powi(4))
    }
    /// Cooling time scale τ_rad = ρ c_p T / (4 κ σ T⁴).
    pub fn cooling_timescale(&self, temp: f64, density: f64, cp: f64) -> f64 {
        if temp < 1e-300 || self.opacity < 1e-300 {
            return f64::INFINITY;
        }
        density * cp * temp / (4.0 * self.opacity * self.sigma_sb * temp.powi(4)).max(1e-300)
    }
    /// SPH temperature rate due to optically thin radiation (K/s).
    pub fn dtemp_dt_rad(&self, temp: f64, density: f64, cp: f64) -> f64 {
        if density < 1e-300 || cp < 1e-300 {
            return 0.0;
        }
        self.volumetric_cooling(temp, density) / (density * cp)
    }
}
/// Temperature-dependent dynamic viscosity models.
#[derive(Clone, Debug)]
pub struct TempDependentViscosity {
    /// Reference temperature T_ref (K).
    pub t_ref: f64,
    /// Reference viscosity μ_ref (Pa·s).
    pub mu_ref: f64,
    /// Model variant.
    pub model: ViscosityModel,
}
impl TempDependentViscosity {
    /// Construct an Arrhenius viscosity model.
    pub fn arrhenius(mu_ref: f64, t_ref: f64, activation_energy: f64) -> Self {
        TempDependentViscosity {
            t_ref,
            mu_ref,
            model: ViscosityModel::Arrhenius { activation_energy },
        }
    }
    /// Construct an Andrade viscosity model.
    pub fn andrade(coeff_a: f64, coeff_b: f64) -> Self {
        TempDependentViscosity {
            t_ref: 300.0,
            mu_ref: coeff_a * (coeff_b / 300.0_f64).exp(),
            model: ViscosityModel::Andrade { coeff_a, coeff_b },
        }
    }
    /// Construct a power-law viscosity model.
    pub fn power_law(mu_ref: f64, t_ref: f64, exponent: f64) -> Self {
        TempDependentViscosity {
            t_ref,
            mu_ref,
            model: ViscosityModel::PowerLaw { exponent },
        }
    }
    /// Evaluate dynamic viscosity μ (Pa·s) at temperature T (K).
    pub fn viscosity(&self, temp: f64) -> f64 {
        if temp < 1e-300 {
            return self.mu_ref;
        }
        match &self.model {
            ViscosityModel::Arrhenius { activation_energy } => {
                let r = 8.314_f64;
                let exp_arg = activation_energy / r * (1.0 / temp - 1.0 / self.t_ref);
                self.mu_ref * exp_arg.exp()
            }
            ViscosityModel::Walther { coeff_a, coeff_b } => {
                let log_t = temp.log10();
                let log_log_nu_07 = coeff_a - coeff_b * log_t;
                let log_nu_07 = 10.0_f64.powf(log_log_nu_07);
                let nu_mm2s = 10.0_f64.powf(log_nu_07) - 0.7;
                let nu_m2s = nu_mm2s * 1e-6_f64;
                nu_m2s * 870.0
            }
            ViscosityModel::Andrade { coeff_a, coeff_b } => coeff_a * (coeff_b / temp).exp(),
            ViscosityModel::PowerLaw { exponent } => {
                self.mu_ref * (temp / self.t_ref).powf(*exponent)
            }
        }
    }
    /// Viscosity ratio μ(T)/μ_ref.
    pub fn ratio(&self, temp: f64) -> f64 {
        if self.mu_ref < 1e-300 {
            return 1.0;
        }
        self.viscosity(temp) / self.mu_ref
    }
}
/// Evaporation at a free surface using Langmuir–Knudsen kinetics.
pub struct EvaporationModel {
    /// Accommodation coefficient α_e ∈ \[0, 1\].
    pub accommodation: f64,
    /// Specific gas constant of vapor R_spec (J/kg/K).
    pub r_spec: f64,
    /// Latent heat of vaporization L_v (J/kg).
    pub latent_heat: f64,
}
impl EvaporationModel {
    /// Construct an evaporation model.
    pub fn new(accommodation: f64, r_spec: f64, latent_heat: f64) -> Self {
        EvaporationModel {
            accommodation,
            r_spec,
            latent_heat,
        }
    }
    /// Mass evaporation rate ṁ (kg/m²/s) at surface temperature T and vapor pressure p_sat.
    pub fn mass_rate(&self, temp: f64, p_sat: f64) -> f64 {
        evaporation_rate(self.accommodation, p_sat, self.r_spec, temp)
    }
    /// Recoil pressure p_r = ṁ² / ρ_vapor (simplified).
    ///
    /// `rho_vapor` vapor density near surface (kg/m³).
    pub fn recoil_pressure(&self, temp: f64, p_sat: f64, rho_vapor: f64) -> f64 {
        let mdot = self.mass_rate(temp, p_sat);
        if rho_vapor < 1e-300 {
            return 0.0;
        }
        mdot * mdot / rho_vapor
    }
    /// Heat sink from evaporation: q_evap = ṁ * L_v (W/m²).
    pub fn heat_sink(&self, temp: f64, p_sat: f64) -> f64 {
        self.mass_rate(temp, p_sat) * self.latent_heat
    }
}
/// Solidification front tracker for industrial casting simulation.
///
/// Implements Stefan-problem tracking, mushy-zone growth, and
/// thermal resistance of the mould-metal gap.
pub struct CastingSolidification {
    /// Phase-change model.
    pub phase_change: PhaseChangeSph,
    /// Mould heat transfer coefficient h_mould (W/m²/K).
    pub h_mould: f64,
    /// Mould temperature T_mould (K).
    pub t_mould: f64,
    /// Initial pour temperature T_pour (K).
    pub t_pour: f64,
    /// Metal thermal conductivity (solid) λ_s (W/m/K).
    pub lambda_solid: f64,
    /// Metal thermal conductivity (liquid) λ_l (W/m/K).
    pub lambda_liquid: f64,
    /// Metal density ρ (kg/m³).
    pub density: f64,
    /// Current solidification front position s (m).
    pub front_position: f64,
}
impl CastingSolidification {
    /// Construct a casting solidification model.
    pub fn new(
        t_melt: f64,
        t_liquidus: f64,
        latent_heat: f64,
        cp: f64,
        h_mould: f64,
        t_mould: f64,
        t_pour: f64,
        lambda_solid: f64,
        lambda_liquid: f64,
        density: f64,
    ) -> Self {
        CastingSolidification {
            phase_change: PhaseChangeSph::new(t_melt, t_liquidus, latent_heat, cp),
            h_mould,
            t_mould,
            t_pour,
            lambda_solid,
            lambda_liquid,
            density,
            front_position: 0.0,
        }
    }
    /// Chvorinov's rule: solidification time t ∝ (V/A)^2.
    ///
    /// t_sol = C_m * (V/A)^2 where C_m is the mould constant (s/m²).
    pub fn chvorinov_time(&self, volume: f64, area: f64, c_mould: f64) -> f64 {
        if area < 1e-300 {
            return 0.0;
        }
        let modulus = volume / area;
        c_mould * modulus * modulus
    }
    /// Temperature at solidification front (Stefan's solution, 1D semi-infinite).
    ///
    /// Uses simplified Neumann solution for 1D solidification.
    pub fn front_temperature(&self) -> f64 {
        self.phase_change.t_melt
    }
    /// Update front position using Stefan condition.
    ///
    /// ds/dt = (λ_l ∂T/∂x|_liquid) / (ρ L)
    pub fn advance_front(&mut self, grad_t_liquid: f64, dt: f64) {
        let stefan = StefanCondition::new(
            self.phase_change.latent_heat,
            self.density,
            self.lambda_solid,
            self.lambda_liquid,
        );
        let v_n = stefan.interface_velocity(0.0, grad_t_liquid);
        self.front_position += v_n * dt;
    }
    /// Heat extracted by mould per unit area (W/m²).
    pub fn mould_heat_flux(&self, t_surface: f64) -> f64 {
        self.h_mould * (t_surface - self.t_mould)
    }
    /// Fraction solidified at time t (Chvorinov-based estimate).
    pub fn fraction_solidified(&self, t_sol: f64, current_time: f64) -> f64 {
        if t_sol < 1e-300 {
            return 1.0;
        }
        (current_time / t_sol).sqrt().min(1.0)
    }
}
/// Temperature-dependent thermal conductivity model.
pub struct TempDependentConductivity {
    /// Reference conductivity λ_ref (W/m/K).
    pub lambda_ref: f64,
    /// Reference temperature T_ref (K).
    pub t_ref: f64,
    /// Linear coefficient dλ/dT (W/m/K²).
    pub linear_coeff: f64,
    /// Optional quadratic coefficient (W/m/K³).
    pub quad_coeff: f64,
}
impl TempDependentConductivity {
    /// Construct a linear temperature-dependent conductivity.
    pub fn linear(lambda_ref: f64, t_ref: f64, dlambda_dt: f64) -> Self {
        TempDependentConductivity {
            lambda_ref,
            t_ref,
            linear_coeff: dlambda_dt,
            quad_coeff: 0.0,
        }
    }
    /// Construct a quadratic temperature-dependent conductivity.
    pub fn quadratic(lambda_ref: f64, t_ref: f64, a1: f64, a2: f64) -> Self {
        TempDependentConductivity {
            lambda_ref,
            t_ref,
            linear_coeff: a1,
            quad_coeff: a2,
        }
    }
    /// Evaluate conductivity λ(T) (W/m/K).
    pub fn conductivity(&self, temp: f64) -> f64 {
        let dt = temp - self.t_ref;
        (self.lambda_ref + self.linear_coeff * dt + self.quad_coeff * dt * dt).max(0.0)
    }
    /// Effective conductivity averaged over temperature range \[T_a, T_b\].
    pub fn averaged(&self, t_a: f64, t_b: f64) -> f64 {
        if (t_b - t_a).abs() < 1e-300 {
            return self.conductivity(t_a);
        }
        let n = 10usize;
        let dt = (t_b - t_a) / n as f64;
        let mut sum = 0.0;
        for k in 0..n {
            let t_mid = t_a + (k as f64 + 0.5) * dt;
            sum += self.conductivity(t_mid);
        }
        sum / n as f64
    }
}
/// SPH energy equation including pressure work and gravity potential.
///
/// Full energy conservation: de/dt = -p/ρ² ∇·v + κ/ρ ∇²T + Φ/ρ + g·v
/// where e is specific internal energy.
#[derive(Clone, Debug)]
pub struct FullEnergySph {
    /// Smoothing length h (m).
    pub h: f64,
    /// Gravity vector \[gx, gy, gz\] (m/s²).
    pub gravity: [f64; 3],
    /// Ratio of specific heats γ.
    pub gamma: f64,
}
impl FullEnergySph {
    /// Create a full SPH energy equation model.
    pub fn new(h: f64, gravity: [f64; 3], gamma: f64) -> Self {
        Self { h, gravity, gamma }
    }
    /// Gravity work term g·v (W/kg).
    pub fn gravity_work(&self, velocity: [f64; 3]) -> f64 {
        dot3(self.gravity, velocity)
    }
    /// Ideal-gas pressure p = (γ-1) ρ e.
    pub fn ideal_gas_pressure(&self, density: f64, internal_energy: f64) -> f64 {
        (self.gamma - 1.0) * density * internal_energy
    }
    /// SPH pressure work contribution to de/dt: -p/(ρ²) * ∇·v.
    pub fn pressure_work_energy(&self, pressure: f64, density: f64, div_v: f64) -> f64 {
        if density < 1e-300 {
            return 0.0;
        }
        -pressure / (density * density) * div_v
    }
    /// Full specific internal energy rate de/dt.
    ///
    /// `conduction_rate`: thermal conduction term (W/kg).
    /// `dissipation_rate`: viscous dissipation per unit mass (W/kg).
    pub fn de_dt(
        &self,
        pressure: f64,
        density: f64,
        div_v: f64,
        velocity: [f64; 3],
        conduction_rate: f64,
        dissipation_rate: f64,
    ) -> f64 {
        self.pressure_work_energy(pressure, density, div_v)
            + conduction_rate
            + dissipation_rate
            + self.gravity_work(velocity)
    }
    /// Speed of sound c = sqrt(γ p / ρ).
    pub fn sound_speed(&self, pressure: f64, density: f64) -> f64 {
        if density < 1e-300 || pressure < 0.0 {
            return 0.0;
        }
        (self.gamma * pressure / density).sqrt()
    }
}
/// Cryogenic fluid properties for SPH simulation.
///
/// Supports liquid nitrogen, liquid helium, and liquid hydrogen.
#[derive(Clone, Debug)]
pub struct CryogenicFluid {
    /// Fluid name.
    pub name: &'static str,
    /// Normal boiling point T_nbp (K).
    pub t_nbp: f64,
    /// Liquid density at NBP ρ_l (kg/m³).
    pub rho_liquid: f64,
    /// Vapor density at NBP ρ_v (kg/m³).
    pub rho_vapor: f64,
    /// Latent heat of vaporization L (J/kg).
    pub latent_heat: f64,
    /// Liquid specific heat c_p (J/kg/K).
    pub cp_liquid: f64,
    /// Liquid thermal conductivity λ (W/m/K).
    pub lambda_liquid: f64,
    /// Liquid dynamic viscosity μ (Pa·s).
    pub mu_liquid: f64,
    /// Surface tension σ (N/m).
    pub surface_tension: f64,
}
impl CryogenicFluid {
    /// Liquid nitrogen properties at 1 atm (77.4 K).
    pub fn liquid_nitrogen() -> Self {
        CryogenicFluid {
            name: "LN2",
            t_nbp: 77.4,
            rho_liquid: 808.0,
            rho_vapor: 4.6,
            latent_heat: 199_000.0,
            cp_liquid: 2040.0,
            lambda_liquid: 0.139,
            mu_liquid: 1.6e-4,
            surface_tension: 0.0088,
        }
    }
    /// Liquid helium-4 properties at 1 atm (4.2 K).
    pub fn liquid_helium() -> Self {
        CryogenicFluid {
            name: "LHe",
            t_nbp: 4.2,
            rho_liquid: 125.0,
            rho_vapor: 16.9,
            latent_heat: 20_400.0,
            cp_liquid: 4800.0,
            lambda_liquid: 0.020,
            mu_liquid: 3.57e-6,
            surface_tension: 0.000353,
        }
    }
    /// Liquid hydrogen properties at 1 atm (20.3 K).
    pub fn liquid_hydrogen() -> Self {
        CryogenicFluid {
            name: "LH2",
            t_nbp: 20.3,
            rho_liquid: 70.8,
            rho_vapor: 1.35,
            latent_heat: 446_000.0,
            cp_liquid: 9_700.0,
            lambda_liquid: 0.099,
            mu_liquid: 1.32e-5,
            surface_tension: 0.00193,
        }
    }
    /// Thermal diffusivity α = λ / (ρ c_p) (m²/s).
    pub fn thermal_diffusivity(&self) -> f64 {
        let denom = self.rho_liquid * self.cp_liquid;
        if denom < 1e-300 {
            return 0.0;
        }
        self.lambda_liquid / denom
    }
    /// Prandtl number Pr = μ c_p / λ.
    pub fn prandtl_number(&self) -> f64 {
        if self.lambda_liquid < 1e-300 {
            return 0.0;
        }
        self.mu_liquid * self.cp_liquid / self.lambda_liquid
    }
    /// Jakob number for superheated vapor: Ja = ρ_v c_pv ΔT / (ρ_l L).
    ///
    /// Here using a simplified form.
    pub fn jakob_number(&self, superheat: f64) -> f64 {
        let denom = self.rho_liquid * self.latent_heat;
        if denom < 1e-300 {
            return 0.0;
        }
        self.rho_vapor * self.cp_liquid * superheat / denom
    }
    /// Saturation pressure estimate using Clausius-Clapeyron (Pa).
    ///
    /// p_sat(T) ≈ p_ref * exp(-L/R_spec * (1/T - 1/T_ref))
    pub fn saturation_pressure(&self, temp: f64, p_ref: f64, r_spec: f64) -> f64 {
        if temp < 1e-300 || r_spec < 1e-300 {
            return 0.0;
        }
        let exp_arg = -self.latent_heat / r_spec * (1.0 / temp - 1.0 / self.t_nbp);
        p_ref * exp_arg.exp()
    }
}
/// Heat flux BC for wall particles in SPH.
pub struct SphWallHeat {
    /// Prescribed wall temperature T_w (K) for Dirichlet BC.
    pub wall_temperature: Option<f64>,
    /// Prescribed heat flux q_w (W/m²) for Neumann BC (positive = into fluid).
    pub heat_flux: Option<f64>,
    /// Thermal conductivity of fluid λ (W/m/K).
    pub conductivity: f64,
}
impl SphWallHeat {
    /// Construct a wall heat BC with Dirichlet temperature.
    pub fn dirichlet(wall_temperature: f64, conductivity: f64) -> Self {
        SphWallHeat {
            wall_temperature: Some(wall_temperature),
            heat_flux: None,
            conductivity,
        }
    }
    /// Construct a wall heat BC with Neumann heat flux.
    pub fn neumann(heat_flux: f64, conductivity: f64) -> Self {
        SphWallHeat {
            wall_temperature: None,
            heat_flux: Some(heat_flux),
            conductivity,
        }
    }
    /// Apply Dirichlet BC: mirror temperature of wall particle.
    pub fn apply_dirichlet(&self, t_fluid: f64, distance: f64) -> f64 {
        if let Some(t_w) = self.wall_temperature {
            let _ = (t_fluid, distance);
            2.0 * t_w - t_fluid
        } else {
            t_fluid
        }
    }
    /// Apply Neumann BC: heat flux contribution to temperature rate.
    pub fn apply_neumann(&self, dt: f64, density: f64, specific_heat: f64, area: f64) -> f64 {
        if let Some(q_w) = self.heat_flux {
            q_w * area * dt / (density * specific_heat)
        } else {
            0.0
        }
    }
}
/// Nucleate boiling source term for SPH energy equation.
pub struct NucleateBoilingSph {
    /// Heat flux at transition to nucleate boiling q_ONB (W/m²).
    pub q_onb: f64,
    /// Boiling model reference.
    pub boiling: BoilingSph,
}
impl NucleateBoilingSph {
    /// Construct a nucleate boiling SPH source term.
    pub fn new(boiling: BoilingSph) -> Self {
        let q_onb = 1000.0;
        NucleateBoilingSph { q_onb, boiling }
    }
    /// Is nucleate boiling active for this wall superheat?
    pub fn is_active(&self, t_wall: f64) -> bool {
        t_wall > self.boiling.t_sat
    }
    /// Volume fraction of vapor in SPH particle (bubble fraction).
    pub fn vapor_fraction(&self, t_wall: f64) -> f64 {
        let ja = self.boiling.jakob_number(t_wall);
        (ja / 100.0).clamp(0.0, 1.0)
    }
    /// Effective density of vapor-liquid mixture ρ_mix = (1-α)ρ_l + α ρ_v.
    pub fn mixture_density(&self, t_wall: f64) -> f64 {
        let alpha = self.vapor_fraction(t_wall);
        (1.0 - alpha) * self.boiling.rho_l + alpha * self.boiling.rho_v
    }
}
/// SPH phase-change model using the enthalpy method.
pub struct PhaseChangeSph {
    /// Melting temperature T_melt (K).
    pub t_melt: f64,
    /// Liquidus temperature T_liquidus (K).
    pub t_liquidus: f64,
    /// Latent heat of fusion L (J/kg).
    pub latent_heat: f64,
    /// Specific heat c_p (J/kg/K).
    pub specific_heat: f64,
}
impl PhaseChangeSph {
    /// Construct a phase-change model.
    pub fn new(t_melt: f64, t_liquidus: f64, latent_heat: f64, specific_heat: f64) -> Self {
        PhaseChangeSph {
            t_melt,
            t_liquidus,
            latent_heat,
            specific_heat,
        }
    }
    /// Compute enthalpy from temperature.
    pub fn enthalpy(&self, temp: f64) -> f64 {
        enthalpy_from_temp(
            temp,
            self.specific_heat,
            self.t_melt,
            self.t_liquidus,
            self.latent_heat,
        )
    }
    /// Compute temperature from enthalpy (inverse of enthalpy method).
    pub fn temperature_from_enthalpy(&self, h: f64) -> f64 {
        let h_melt = self.specific_heat * self.t_melt;
        let h_liquid = h_melt + self.latent_heat;
        if h < h_melt {
            h / self.specific_heat.max(1e-300)
        } else if h < h_liquid {
            let range = (self.t_liquidus - self.t_melt).max(1e-300);
            self.t_melt + (h - h_melt) / self.latent_heat * range
        } else {
            (h - self.latent_heat) / self.specific_heat.max(1e-300)
        }
    }
    /// Liquid fraction f_l ∈ \[0, 1\].
    pub fn liquid_fraction(&self, temp: f64) -> f64 {
        if temp <= self.t_melt {
            return 0.0;
        }
        if temp >= self.t_liquidus {
            return 1.0;
        }
        let range = (self.t_liquidus - self.t_melt).max(1e-300);
        (temp - self.t_melt) / range
    }
    /// Is the particle in the mushy zone?
    pub fn is_mushy(&self, temp: f64) -> bool {
        temp > self.t_melt && temp < self.t_liquidus
    }
}
/// Species transport coupled with heat transfer (thermodiffusion / Soret effect).
pub struct ThermodiffusionSph {
    /// Soret coefficient D_T (m²/s/K), also called thermal diffusion coefficient.
    pub soret_coeff: f64,
    /// Mass diffusivity D (m²/s).
    pub mass_diffusivity: f64,
    /// Reference temperature T_ref (K).
    pub t_ref: f64,
}
impl ThermodiffusionSph {
    /// Construct a thermodiffusion model.
    pub fn new(soret_coeff: f64, mass_diffusivity: f64, t_ref: f64) -> Self {
        ThermodiffusionSph {
            soret_coeff,
            mass_diffusivity,
            t_ref,
        }
    }
    /// Soret mass flux J_S = -ρ D_T c(1-c) ∇T/T (kg/m²/s).
    pub fn soret_flux(
        &self,
        concentration: f64,
        grad_t: [f64; 3],
        density: f64,
        temp: f64,
    ) -> [f64; 3] {
        if temp.abs() < 1e-300 {
            return [0.0; 3];
        }
        let coeff = -density * self.soret_coeff * concentration * (1.0 - concentration) / temp;
        scale3(grad_t, coeff)
    }
    /// Dufour heat flux J_D = -ρ D_D c(1-c) ∇c * c_p (W/m²).
    ///
    /// `grad_c` concentration gradient, `cp` specific heat.
    pub fn dufour_flux(
        &self,
        concentration: f64,
        grad_c: [f64; 3],
        density: f64,
        cp: f64,
    ) -> [f64; 3] {
        let coeff = -density * self.mass_diffusivity * concentration * (1.0 - concentration) * cp;
        scale3(grad_c, coeff)
    }
    /// Effective mass diffusivity including Soret contribution at given T.
    pub fn effective_diffusivity(&self, temp: f64) -> f64 {
        if temp.abs() < 1e-300 {
            return self.mass_diffusivity;
        }
        self.mass_diffusivity + self.soret_coeff * (temp - self.t_ref)
    }
}
