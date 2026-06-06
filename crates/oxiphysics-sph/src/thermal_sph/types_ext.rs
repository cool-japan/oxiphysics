//! Extended thermal SPH types: boundary conditions, radiation, convection, and more.

use std::f64::consts::PI;

use super::functions::*;
use crate::thermal_sph::types::*;

/// SPH thermal boundary condition applier.
///
/// Applies both Dirichlet (prescribed temperature) and Neumann (prescribed
/// heat flux) thermal boundary conditions in SPH using ghost particles.
#[derive(Clone, Debug)]
pub struct SphThermalBc {
    /// Boundary condition kind.
    pub bc_kind: ThermalBcKind,
    /// Smoothing length h (m).
    pub h: f64,
}
impl SphThermalBc {
    /// Create a Dirichlet BC.
    pub fn dirichlet(t_wall: f64, h: f64) -> Self {
        Self {
            bc_kind: ThermalBcKind::Dirichlet(t_wall),
            h,
        }
    }
    /// Create a Neumann BC.
    pub fn neumann(q_flux: f64, h: f64) -> Self {
        Self {
            bc_kind: ThermalBcKind::Neumann(q_flux),
            h,
        }
    }
    /// Create a Robin (convective) BC.
    pub fn robin(h_conv: f64, t_ambient: f64, h_sph: f64) -> Self {
        Self {
            bc_kind: ThermalBcKind::Robin { h_conv, t_ambient },
            h: h_sph,
        }
    }
    /// Apply the BC: returns the effective temperature for the ghost/wall particle.
    ///
    /// For Dirichlet: T_ghost = 2 T_w - T_fluid.
    /// For Neumann: T_ghost = T_fluid + q_w * d / λ (d = distance).
    /// For Robin: T_ghost = T_fluid + (h_c (T_a - T_fluid)) * d / λ.
    pub fn ghost_temperature(&self, t_fluid: f64, distance: f64, lambda: f64) -> f64 {
        match self.bc_kind {
            ThermalBcKind::Dirichlet(t_w) => 2.0 * t_w - t_fluid,
            ThermalBcKind::Neumann(q_w) => {
                if lambda < 1e-300 {
                    return t_fluid;
                }
                t_fluid + q_w * distance / lambda
            }
            ThermalBcKind::Robin { h_conv, t_ambient } => {
                let q_eff = h_conv * (t_ambient - t_fluid);
                if lambda < 1e-300 {
                    return t_fluid;
                }
                t_fluid + q_eff * distance / lambda
            }
        }
    }
    /// Heat flux from the wall to the fluid (W/m²).
    pub fn wall_heat_flux(&self, t_fluid: f64, lambda: f64, distance: f64) -> f64 {
        match self.bc_kind {
            ThermalBcKind::Dirichlet(t_w) => {
                if distance < 1e-300 {
                    return 0.0;
                }
                lambda * (t_w - t_fluid) / distance
            }
            ThermalBcKind::Neumann(q_w) => q_w,
            ThermalBcKind::Robin { h_conv, t_ambient } => h_conv * (t_ambient - t_fluid),
        }
    }
}
/// Extended temperature-dependent viscosity with cross-model blending.
///
/// Blends multiple viscosity models for different temperature ranges,
/// for example switching from Arrhenius below melting to power-law above.
#[derive(Clone, Debug)]
pub struct BlendedViscosity {
    /// Viscosity model for low temperature (T < T_blend).
    pub low_temp_model: TempDependentViscosity,
    /// Viscosity model for high temperature (T >= T_blend).
    pub high_temp_model: TempDependentViscosity,
    /// Blending temperature T_blend (K).
    pub t_blend: f64,
    /// Blending width ΔT (K).
    pub blend_width: f64,
}
impl BlendedViscosity {
    /// Create a blended viscosity model.
    pub fn new(
        low_temp_model: TempDependentViscosity,
        high_temp_model: TempDependentViscosity,
        t_blend: f64,
        blend_width: f64,
    ) -> Self {
        Self {
            low_temp_model,
            high_temp_model,
            t_blend,
            blend_width,
        }
    }
    /// Smooth blending weight w ∈ \[0, 1\] (0 = low model, 1 = high model).
    fn blend_weight(&self, temp: f64) -> f64 {
        let x = (temp - self.t_blend) / self.blend_width.max(1e-300);
        1.0 / (1.0 + (-x).exp())
    }
    /// Blended dynamic viscosity μ (Pa·s) at temperature T.
    pub fn viscosity(&self, temp: f64) -> f64 {
        let w = self.blend_weight(temp);
        let mu_low = self.low_temp_model.viscosity(temp);
        let mu_high = self.high_temp_model.viscosity(temp);
        (1.0 - w) * mu_low + w * mu_high
    }
    /// Blended viscosity derivative dμ/dT (Pa·s/K).
    pub fn viscosity_derivative(&self, temp: f64) -> f64 {
        let dt = 1e-4;
        (self.viscosity(temp + dt) - self.viscosity(temp - dt)) / (2.0 * dt)
    }
}
/// Anisotropic thermal conductivity tensor (3×3, row-major).
///
/// Stores full symmetric conductivity tensor λ_ij for anisotropic materials
/// such as composites, crystals, or geologic layers.
#[derive(Clone, Debug)]
pub struct ConductivityTensor {
    /// Row-major 3×3 tensor components λ_ij (W/m/K).
    pub lambda: [[f64; 3]; 3],
}
impl ConductivityTensor {
    /// Construct an isotropic conductivity tensor λ_ij = λ δ_ij.
    pub fn isotropic(lambda: f64) -> Self {
        let mut l = [[0.0f64; 3]; 3];
        l[0][0] = lambda;
        l[1][1] = lambda;
        l[2][2] = lambda;
        ConductivityTensor { lambda: l }
    }
    /// Construct an orthotropic tensor with principal values (λx, λy, λz).
    pub fn orthotropic(lx: f64, ly: f64, lz: f64) -> Self {
        let mut l = [[0.0f64; 3]; 3];
        l[0][0] = lx;
        l[1][1] = ly;
        l[2][2] = lz;
        ConductivityTensor { lambda: l }
    }
    /// Apply tensor to vector: result_i = Σ_j λ_ij v_j.
    pub fn apply(&self, v: [f64; 3]) -> [f64; 3] {
        let mut out = [0.0f64; 3];
        for (i, out_i) in out.iter_mut().enumerate() {
            for (j, &vj) in v.iter().enumerate() {
                *out_i += self.lambda[i][j] * vj;
            }
        }
        out
    }
    /// Effective conductivity in direction `n` (unit vector): n^T λ n.
    pub fn effective(&self, n: [f64; 3]) -> f64 {
        let lv = self.apply(n);
        dot3(n, lv)
    }
    /// Frobenius norm of the tensor.
    pub fn frobenius_norm(&self) -> f64 {
        let mut s = 0.0f64;
        for row in &self.lambda {
            for &v in row {
                s += v * v;
            }
        }
        s.sqrt()
    }
}
/// Stefan condition for solid-liquid interface velocity.
///
/// The Stefan condition states: ρ L V_n = λ_s (∂T/∂n)_s − λ_l (∂T/∂n)_l
/// where V_n is normal interface velocity, λ_s/λ_l are solid/liquid
/// conductivities, and the derivatives are normal temperature gradients.
pub struct StefanCondition {
    /// Latent heat L (J/kg).
    pub latent_heat: f64,
    /// Density at interface ρ (kg/m³).
    pub density: f64,
    /// Thermal conductivity on solid side λ_s (W/m/K).
    pub lambda_solid: f64,
    /// Thermal conductivity on liquid side λ_l (W/m/K).
    pub lambda_liquid: f64,
}
impl StefanCondition {
    /// Construct a Stefan condition model.
    pub fn new(latent_heat: f64, density: f64, lambda_solid: f64, lambda_liquid: f64) -> Self {
        StefanCondition {
            latent_heat,
            density,
            lambda_solid,
            lambda_liquid,
        }
    }
    /// Compute interface normal velocity V_n (m/s).
    ///
    /// `grad_t_solid` temperature gradient in solid (K/m, positive into solid),
    /// `grad_t_liquid` temperature gradient in liquid (K/m, positive into liquid).
    pub fn interface_velocity(&self, grad_t_solid: f64, grad_t_liquid: f64) -> f64 {
        let denom = self.density * self.latent_heat;
        if denom.abs() < 1e-300 {
            return 0.0;
        }
        (self.lambda_solid * grad_t_solid - self.lambda_liquid * grad_t_liquid) / denom
    }
    /// Heat flux imbalance at interface (W/m²).
    pub fn heat_flux_imbalance(&self, grad_t_solid: f64, grad_t_liquid: f64) -> f64 {
        self.lambda_solid * grad_t_solid - self.lambda_liquid * grad_t_liquid
    }
    /// Solidification rate (kg/m²/s) = ρ V_n.
    pub fn solidification_rate(&self, grad_t_solid: f64, grad_t_liquid: f64) -> f64 {
        self.density * self.interface_velocity(grad_t_solid, grad_t_liquid)
    }
}
/// Viscosity-temperature correlation model.
#[derive(Clone, Debug)]
pub enum ViscosityModel {
    /// Arrhenius: μ = μ_ref exp(E_a/R (1/T - 1/T_ref)).
    Arrhenius {
        /// Activation energy E_a (J/mol).
        activation_energy: f64,
    },
    /// Walther (lubricant grade): log10(log10(ν + 0.7)) = A - B log10(T).
    Walther {
        /// Coefficient A.
        coeff_a: f64,
        /// Coefficient B.
        coeff_b: f64,
    },
    /// Andrade: μ = A exp(B/T).
    Andrade {
        /// Pre-exponential factor A (Pa·s).
        coeff_a: f64,
        /// Andrade coefficient B (K).
        coeff_b: f64,
    },
    /// Power law: μ = μ_ref (T/T_ref)^n.
    PowerLaw {
        /// Exponent n (negative for liquids).
        exponent: f64,
    },
}
/// Heat pipe SPH model with effective conductivity enhancement.
///
/// Heat pipes transport heat via evaporation-condensation cycles with
/// an effective conductivity orders of magnitude higher than the wall material.
pub struct HeatPipeSph {
    /// Effective thermal conductivity λ_eff (W/m/K).
    pub lambda_eff: f64,
    /// Operating temperature range \[T_min, T_max\] (K).
    pub t_min: f64,
    /// Maximum operating temperature (K).
    pub t_max: f64,
    /// Maximum heat transport rate Q_max (W).
    pub q_max: f64,
    /// Length of heat pipe L (m).
    pub length: f64,
}
impl HeatPipeSph {
    /// Construct a heat pipe model.
    pub fn new(lambda_eff: f64, t_min: f64, t_max: f64, q_max: f64, length: f64) -> Self {
        HeatPipeSph {
            lambda_eff,
            t_min,
            t_max,
            q_max,
            length,
        }
    }
    /// Water heat pipe with groove wick (approximation).
    pub fn water_heat_pipe(length: f64) -> Self {
        HeatPipeSph {
            lambda_eff: 10_000.0,
            t_min: 320.0,
            t_max: 500.0,
            q_max: 500.0,
            length,
        }
    }
    /// Is the heat pipe operating within its temperature range?
    pub fn is_operating(&self, temp: f64) -> bool {
        temp >= self.t_min && temp <= self.t_max
    }
    /// Heat flux through the heat pipe for a given temperature difference.
    pub fn heat_flux(&self, delta_t: f64) -> f64 {
        if self.length < 1e-300 {
            return 0.0;
        }
        (self.lambda_eff * delta_t / self.length).min(self.q_max / 1.0)
    }
    /// Thermal resistance R = L / (λ_eff A) for cross-section area A (m²).
    pub fn thermal_resistance(&self, area: f64) -> f64 {
        if self.lambda_eff < 1e-300 || area < 1e-300 {
            return f64::INFINITY;
        }
        self.length / (self.lambda_eff * area)
    }
}
/// Viscous dissipation heating model for SPH.
///
/// Computes the local heat generated by viscous stresses in the flow.
/// In an incompressible Newtonian fluid: Φ = 2 μ S_ij S_ij where S_ij
/// is the strain-rate tensor.
#[derive(Clone, Debug)]
pub struct ViscousDissipationSph {
    /// Dynamic viscosity μ (Pa·s).
    pub mu: f64,
    /// SPH smoothing length (m).
    pub h: f64,
}
impl ViscousDissipationSph {
    /// Create a viscous dissipation model.
    pub fn new(mu: f64, h: f64) -> Self {
        Self { mu, h }
    }
    /// SPH estimate of velocity gradient tensor component ∂u_α/∂x_β at particle i.
    ///
    /// (∂u_α/∂x_β)_i = Σ_j m_j/ρ_j (u_α,j - u_α,i) (∂W/∂x_β)_ij
    pub fn velocity_gradient(
        v_i: [f64; 3],
        neighbors: &[([f64; 3], [f64; 3], f64, f64)],
        h: f64,
    ) -> [[f64; 3]; 3] {
        let mut grad = [[0.0f64; 3]; 3];
        for &(r_ij, v_j, m_j, rho_j) in neighbors {
            if rho_j < 1e-300 {
                continue;
            }
            let gw = cubic_kernel_grad(r_ij, h);
            let dv = sub3(v_j, v_i);
            for a in 0..3 {
                for b in 0..3 {
                    grad[a][b] += m_j / rho_j * dv[a] * gw[b];
                }
            }
        }
        grad
    }
    /// Strain-rate tensor S_ij = (∂u_i/∂x_j + ∂u_j/∂x_i) / 2.
    pub fn strain_rate(grad_v: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
        let mut s = [[0.0f64; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                s[i][j] = 0.5 * (grad_v[i][j] + grad_v[j][i]);
            }
        }
        s
    }
    /// Double contraction S:S = Σ_ij S_ij².
    pub fn strain_rate_invariant(s: [[f64; 3]; 3]) -> f64 {
        let mut inv2 = 0.0;
        for row in &s {
            for &v in row.iter() {
                inv2 += v * v;
            }
        }
        inv2
    }
    /// Dissipation rate per unit volume Φ = 2 μ S:S (W/m³).
    pub fn dissipation_rate(
        &self,
        v_i: [f64; 3],
        neighbors: &[([f64; 3], [f64; 3], f64, f64)],
    ) -> f64 {
        let grad_v = Self::velocity_gradient(v_i, neighbors, self.h);
        let s = Self::strain_rate(grad_v);
        let inv2 = Self::strain_rate_invariant(s);
        2.0 * self.mu * inv2
    }
    /// Temperature rate from dissipation dT/dt = Φ / (ρ c_p) (K/s).
    pub fn dtemp_dt(
        &self,
        v_i: [f64; 3],
        rho_i: f64,
        cp_i: f64,
        neighbors: &[([f64; 3], [f64; 3], f64, f64)],
    ) -> f64 {
        if rho_i < 1e-300 || cp_i < 1e-300 {
            return 0.0;
        }
        self.dissipation_rate(v_i, neighbors) / (rho_i * cp_i)
    }
}
/// Parameters for a thermal SPH simulation.
#[derive(Clone, Debug)]
pub struct ThermalSPHParams {
    /// Kernel smoothing length h (m).
    pub kernel_radius: f64,
    /// Equation-of-state stiffness k.
    pub eos_k: f64,
    /// Equation-of-state exponent γ.
    pub eos_gamma: f64,
    /// Dynamic viscosity μ (Pa·s).
    pub viscosity: f64,
}
impl ThermalSPHParams {
    /// Create a new [`ThermalSPHParams`].
    pub fn new(kernel_radius: f64, eos_k: f64, eos_gamma: f64, viscosity: f64) -> Self {
        Self {
            kernel_radius,
            eos_k,
            eos_gamma,
            viscosity,
        }
    }
}
/// Kind of thermal boundary condition.
#[derive(Clone, Debug, PartialEq)]
pub enum ThermalBcKind {
    /// Fixed temperature T_w (K) — Dirichlet.
    Dirichlet(f64),
    /// Fixed heat flux q_w (W/m²) — Neumann.
    Neumann(f64),
    /// Convective: h_c (T_fluid - T_ambient) — Robin.
    Robin {
        /// Convective heat transfer coefficient h_c (W/m²/K).
        h_conv: f64,
        /// Ambient temperature T_amb (K).
        t_ambient: f64,
    },
}
/// Radiation heat transfer in SPH simulations.
///
/// Implements Stefan-Boltzmann law and simplified view-factor calculations
/// for participating media with absorption and emission.
pub struct RadiationSph {
    /// Stefan-Boltzmann constant σ (W/m²/K⁴).
    pub sigma_sb: f64,
    /// Emissivity ε ∈ \[0, 1\].
    pub emissivity: f64,
    /// Absorption coefficient κ_a (1/m).
    pub absorption_coeff: f64,
    /// Scattering coefficient κ_s (1/m).
    pub scattering_coeff: f64,
}
impl RadiationSph {
    /// Construct a radiation model with default Stefan-Boltzmann constant.
    pub fn new(emissivity: f64, absorption_coeff: f64, scattering_coeff: f64) -> Self {
        RadiationSph {
            sigma_sb: 5.670_374_4e-8,
            emissivity,
            absorption_coeff,
            scattering_coeff,
        }
    }
    /// Blackbody emission power E_b = σ T^4 (W/m²).
    pub fn blackbody_emission(&self, temp: f64) -> f64 {
        self.sigma_sb * temp.powi(4)
    }
    /// Net radiation flux from surface at T_s to environment at T_env (W/m²).
    pub fn net_radiation_flux(&self, t_surface: f64, t_env: f64) -> f64 {
        self.emissivity * self.sigma_sb * (t_surface.powi(4) - t_env.powi(4))
    }
    /// Extinction coefficient κ_ext = κ_a + κ_s (1/m).
    pub fn extinction_coeff(&self) -> f64 {
        self.absorption_coeff + self.scattering_coeff
    }
    /// Optical depth τ = κ_ext * L for path length L (m).
    pub fn optical_depth(&self, path_length: f64) -> f64 {
        self.extinction_coeff() * path_length
    }
    /// Transmittance T = exp(-τ) for path length L.
    pub fn transmittance(&self, path_length: f64) -> f64 {
        (-self.optical_depth(path_length)).exp()
    }
    /// Volumetric radiation source term in SPH (W/m³).
    ///
    /// q_rad = κ_a (4 σ T^4 - G) where G is incident irradiation (W/m²).
    pub fn volumetric_source(&self, temp: f64, irradiation: f64) -> f64 {
        self.absorption_coeff * (4.0 * self.sigma_sb * temp.powi(4) - irradiation)
    }
    /// View factor between two differential surface elements (Lambertian).
    ///
    /// dF_ij = cos(θ_i) cos(θ_j) / (π r²) * dA_j
    ///
    /// `cos_theta_i` cosine of emission angle, `cos_theta_j` cosine of
    /// reception angle, `r` distance (m), `da_j` area of receiving element (m²).
    pub fn view_factor_element(cos_theta_i: f64, cos_theta_j: f64, r: f64, da_j: f64) -> f64 {
        if r < 1e-300 {
            return 0.0;
        }
        (cos_theta_i.max(0.0) * cos_theta_j.max(0.0)) / (PI * r * r) * da_j
    }
    /// Rosseland mean absorption coefficient (diffusion approximation).
    ///
    /// Used in optically thick media: q_rad = -16 σ T³ / (3 κ_R) * ∇T.
    pub fn rosseland_flux(&self, temp: f64, grad_t: [f64; 3]) -> [f64; 3] {
        if self.absorption_coeff < 1e-300 {
            return [0.0; 3];
        }
        let coeff = -16.0 * self.sigma_sb * temp.powi(3) / (3.0 * self.absorption_coeff);
        scale3(grad_t, coeff)
    }
}
/// Symmetric SPH thermal conductivity operator.
///
/// Implements the symmetric formulation that exactly conserves energy:
/// dT_i/dt = Σ_j (m_j/ρ_j) * 2 λ_ij / (λ_i + λ_j) * (T_i - T_j) / |r_ij|² * (∇W_ij · r_ij)
#[derive(Clone, Debug)]
pub struct SymmetricSphConductivity {
    /// SPH smoothing length (m).
    pub h: f64,
}
impl SymmetricSphConductivity {
    /// Create a symmetric conductivity operator.
    pub fn new(h: f64) -> Self {
        Self { h }
    }
    /// Compute symmetric inter-particle conductivity (arithmetic mean).
    ///
    /// Uses 2 λ_i λ_j / (λ_i + λ_j) (harmonic) which prevents unphysical
    /// heat flow across sharp conductivity jumps.
    pub fn inter_particle_conductivity(lambda_i: f64, lambda_j: f64) -> f64 {
        let denom = lambda_i + lambda_j;
        if denom < 1e-300 {
            return 0.0;
        }
        2.0 * lambda_i * lambda_j / denom
    }
    /// Temperature rate dT_i/dt via symmetric formulation.
    ///
    /// `lambda_i`: conductivity of particle i.
    /// `rho_i`: density of particle i.
    /// `cp_i`: specific heat of particle i.
    /// `neighbors`: (r_ij, T_j, m_j, ρ_j, λ_j) tuples.
    pub fn dtemp_dt(
        &self,
        t_i: f64,
        lambda_i: f64,
        rho_i: f64,
        cp_i: f64,
        neighbors: &[([f64; 3], f64, f64, f64, f64)],
    ) -> f64 {
        if rho_i < 1e-300 || cp_i < 1e-300 {
            return 0.0;
        }
        let mut sum = 0.0;
        for &(r_ij, t_j, m_j, rho_j, lambda_j) in neighbors {
            let r2 = dot3(r_ij, r_ij);
            if r2 < 1e-300 || rho_j < 1e-300 {
                continue;
            }
            let lambda_ij = Self::inter_particle_conductivity(lambda_i, lambda_j);
            let gw = cubic_kernel_grad(r_ij, self.h);
            let dot_rg = dot3(r_ij, gw);
            sum += m_j / rho_j * lambda_ij * 2.0 * (t_i - t_j) / r2 * dot_rg;
        }
        sum / (rho_i * cp_i)
    }
    /// Full SPH heat equation step: returns updated temperature.
    pub fn step(
        &self,
        t_i: f64,
        lambda_i: f64,
        rho_i: f64,
        cp_i: f64,
        neighbors: &[([f64; 3], f64, f64, f64, f64)],
        dt: f64,
    ) -> f64 {
        t_i + dt * self.dtemp_dt(t_i, lambda_i, rho_i, cp_i, neighbors)
    }
}
/// Rayleigh–Bénard convection driven by Boussinesq buoyancy.
pub struct RayleighBenardSph {
    /// Thermal expansion coefficient β (1/K).
    pub beta: f64,
    /// Kinematic viscosity ν (m²/s).
    pub nu: f64,
    /// Thermal diffusivity κ (m²/s).
    pub kappa: f64,
    /// Layer height H (m).
    pub height: f64,
    /// Temperature difference ΔT across layer (K).
    pub delta_t: f64,
    /// Gravitational acceleration (m/s²).
    pub g: f64,
}
impl RayleighBenardSph {
    /// Construct a Rayleigh–Bénard model.
    pub fn new(beta: f64, nu: f64, kappa: f64, height: f64, delta_t: f64, g: f64) -> Self {
        RayleighBenardSph {
            beta,
            nu,
            kappa,
            height,
            delta_t,
            g,
        }
    }
    /// Rayleigh number Ra = g β ΔT H³ / (ν κ).
    pub fn rayleigh_number(&self) -> f64 {
        if self.nu < 1e-300 || self.kappa < 1e-300 {
            return 0.0;
        }
        self.g * self.beta * self.delta_t * self.height.powi(3) / (self.nu * self.kappa)
    }
    /// Prandtl number Pr = ν / κ.
    pub fn prandtl_number(&self) -> f64 {
        if self.kappa < 1e-300 {
            return 0.0;
        }
        self.nu / self.kappa
    }
    /// Nusselt number (empirical Grossmann–Lohse scaling): Nu ≈ 0.27 Ra^{0.25}.
    pub fn nusselt_number(&self) -> f64 {
        let ra = self.rayleigh_number();
        if ra < 1.0 {
            return 1.0;
        }
        0.27 * ra.powf(0.25)
    }
    /// Boussinesq force at temperature T given reference T_ref.
    pub fn buoyancy_force(&self, temp: f64, t_ref: f64) -> [f64; 3] {
        let dt = temp - t_ref;
        boussinesq_force_sph(dt, self.beta, [0.0, 0.0, -self.g])
    }
}
/// Local Nusselt number from SPH temperature field.
pub struct ThermalConvectiveCoeff {
    /// Reference length L (m) for Nusselt definition.
    pub ref_length: f64,
    /// Fluid thermal conductivity λ (W/m/K).
    pub conductivity: f64,
    /// Wall temperature T_w (K).
    pub t_wall: f64,
    /// Far-field temperature T_inf (K).
    pub t_inf: f64,
}
impl ThermalConvectiveCoeff {
    /// Construct a convective coefficient calculator.
    pub fn new(ref_length: f64, conductivity: f64, t_wall: f64, t_inf: f64) -> Self {
        ThermalConvectiveCoeff {
            ref_length,
            conductivity,
            t_wall,
            t_inf,
        }
    }
    /// Local heat transfer coefficient h = q_w / (T_w - T_inf).
    pub fn heat_transfer_coeff(&self, q_wall: f64) -> f64 {
        let dt = self.t_wall - self.t_inf;
        if dt.abs() < 1e-300 {
            return 0.0;
        }
        q_wall / dt
    }
    /// Local Nusselt number Nu = h L / λ.
    pub fn nusselt_number(&self, q_wall: f64) -> f64 {
        let h = self.heat_transfer_coeff(q_wall);
        if self.conductivity < 1e-300 {
            return 0.0;
        }
        h * self.ref_length / self.conductivity
    }
    /// Area-averaged Nusselt number from array of wall heat fluxes.
    pub fn average_nusselt(&self, q_wall_values: &[f64]) -> f64 {
        if q_wall_values.is_empty() {
            return 0.0;
        }
        let sum: f64 = q_wall_values.iter().map(|&q| self.nusselt_number(q)).sum();
        sum / q_wall_values.len() as f64
    }
}
/// Conjugate heat transfer model at a fluid-solid interface.
///
/// Couples SPH fluid particles with solid thermal nodes via
/// a contact resistance or direct flux matching.
pub struct ConjugateHeatTransfer {
    /// Convective heat transfer coefficient h_conv (W/m²/K).
    pub h_conv: f64,
    /// Contact resistance R_c (m²K/W).
    pub contact_resistance: f64,
    /// Fluid conductivity λ_f (W/m/K).
    pub lambda_fluid: f64,
    /// Solid conductivity λ_s (W/m/K).
    pub lambda_solid: f64,
}
impl ConjugateHeatTransfer {
    /// Construct a conjugate heat transfer model (no contact resistance).
    pub fn new(h_conv: f64, lambda_fluid: f64, lambda_solid: f64) -> Self {
        ConjugateHeatTransfer {
            h_conv,
            contact_resistance: 0.0,
            lambda_fluid,
            lambda_solid,
        }
    }
    /// Construct with contact (Kapitza) resistance.
    pub fn with_contact_resistance(
        h_conv: f64,
        contact_resistance: f64,
        lambda_fluid: f64,
        lambda_solid: f64,
    ) -> Self {
        ConjugateHeatTransfer {
            h_conv,
            contact_resistance,
            lambda_fluid,
            lambda_solid,
        }
    }
    /// Effective heat transfer coefficient including contact resistance.
    ///
    /// 1/h_eff = 1/h_conv + R_c
    pub fn effective_h(&self) -> f64 {
        let denom = 1.0 / self.h_conv.max(1e-300) + self.contact_resistance;
        if denom < 1e-300 {
            return 0.0;
        }
        1.0 / denom
    }
    /// Heat flux across interface q = h_eff (T_f - T_s) (W/m²).
    pub fn interface_flux(&self, t_fluid: f64, t_solid: f64) -> f64 {
        self.effective_h() * (t_fluid - t_solid)
    }
    /// Nusselt number from Dittus-Boelter correlation for turbulent pipe flow.
    ///
    /// Nu = 0.023 Re^0.8 Pr^0.4 (heating case).
    pub fn dittus_boelter_nu(re: f64, pr: f64) -> f64 {
        0.023 * re.powf(0.8) * pr.powf(0.4)
    }
    /// Temperature at solid surface from Newton's cooling law.
    pub fn solid_surface_temp(&self, t_fluid: f64, q_flux: f64) -> f64 {
        let h = self.effective_h();
        if h < 1e-300 {
            return t_fluid;
        }
        t_fluid - q_flux / h
    }
}
/// SPH heat conduction with Cleary harmonic-mean conductivity.
///
/// The Cleary (1998) formulation uses the harmonic mean of particle
/// conductivities to prevent artificial heat flow across conductivity
/// discontinuities at material interfaces.
pub struct ClearyHeatConduction {
    /// SPH smoothing length (m).
    pub h: f64,
}
impl ClearyHeatConduction {
    /// Construct a Cleary heat conduction model.
    pub fn new(h: f64) -> Self {
        ClearyHeatConduction { h }
    }
    /// Harmonic mean of two conductivities.
    pub fn harmonic_mean(lambda_i: f64, lambda_j: f64) -> f64 {
        let denom = lambda_i + lambda_j;
        if denom < 1e-300 {
            return 0.0;
        }
        2.0 * lambda_i * lambda_j / denom
    }
    /// Temperature rate dT_i/dt using Cleary formulation.
    ///
    /// `neighbors`: (r_ij, T_j, m_j, ρ_j, λ_j) — per-neighbor data.
    /// `lambda_i` conductivity of particle i.
    /// `cp_i` specific heat of particle i.
    /// `rho_i` density of particle i.
    pub fn dtemp_dt(
        &self,
        t_i: f64,
        lambda_i: f64,
        cp_i: f64,
        rho_i: f64,
        neighbors: &[([f64; 3], f64, f64, f64, f64)],
    ) -> f64 {
        if rho_i < 1e-300 || cp_i < 1e-300 {
            return 0.0;
        }
        let mut sum = 0.0;
        for &(r_ij, t_j, m_j, rho_j, lambda_j) in neighbors {
            let r2 = dot3(r_ij, r_ij);
            if r2 < 1e-300 || rho_j < 1e-300 {
                continue;
            }
            let lambda_ij = Self::harmonic_mean(lambda_i, lambda_j);
            let grad_w = cubic_kernel_grad(r_ij, self.h);
            let dot_rg = dot3(r_ij, grad_w);
            sum += m_j / rho_j * lambda_ij * 2.0 * (t_i - t_j) / r2 * dot_rg;
        }
        sum / (rho_i * cp_i)
    }
}
/// Full SPH heat equation solver including pressure work and dissipation.
///
/// Implements the energy equation:
/// ρ c_p DT/Dt = ∇·(λ∇T) + p/ρ Dρ/Dt + Φ_diss + Q_ext
/// in SPH particle form.
#[derive(Clone, Debug)]
pub struct SphHeatEquation {
    /// Thermal conductivity operator.
    pub conduction: SymmetricSphConductivity,
    /// Include viscous dissipation heating.
    pub include_dissipation: bool,
    /// Include pressure work term.
    pub include_pressure_work: bool,
    /// External heat source Q_ext (W/kg).
    pub q_ext: f64,
}
impl SphHeatEquation {
    /// Create a full SPH heat equation solver.
    pub fn new(h: f64, q_ext: f64) -> Self {
        Self {
            conduction: SymmetricSphConductivity::new(h),
            include_dissipation: true,
            include_pressure_work: true,
            q_ext,
        }
    }
    /// Viscous dissipation rate Φ = 2 μ ε_ij ε_ij (simplified trace).
    ///
    /// Uses the SPH approximation: Φ ≈ μ |∇v + (∇v)ᵀ|²/2 ≈ μ * shear_rate².
    pub fn viscous_dissipation(mu: f64, shear_rate: f64) -> f64 {
        mu * shear_rate * shear_rate
    }
    /// Pressure work per unit mass p/(ρ²) * Dρ/Dt ≈ -p/ρ * ∇·v.
    pub fn pressure_work_term(pressure: f64, density: f64, div_v: f64) -> f64 {
        if density < 1e-300 {
            return 0.0;
        }
        -pressure / density * div_v
    }
    /// Full energy rate dT/dt for a particle.
    ///
    /// `conduction_rate`: conduction contribution (K/s).
    /// `pressure`: particle pressure (Pa).
    /// `density`: particle density (kg/m³).
    /// `cp`: specific heat (J/kg/K).
    /// `div_v`: velocity divergence (1/s).
    /// `mu`: dynamic viscosity (Pa·s).
    /// `shear_rate`: local shear rate magnitude (1/s).
    pub fn full_dtemp_dt(
        &self,
        conduction_rate: f64,
        pressure: f64,
        density: f64,
        cp: f64,
        div_v: f64,
        mu: f64,
        shear_rate: f64,
    ) -> f64 {
        let mut rate = conduction_rate + self.q_ext / cp.max(1e-300);
        if self.include_pressure_work {
            rate += Self::pressure_work_term(pressure, density, div_v) / cp.max(1e-300);
        }
        if self.include_dissipation {
            rate +=
                Self::viscous_dissipation(mu, shear_rate) / (density.max(1e-300) * cp.max(1e-300));
        }
        rate
    }
}
/// Thermal shock model for brittle fracture triggered by thermal stresses.
///
/// Computes the thermal stress intensity factor and predicts crack initiation
/// based on a critical stress criterion.
#[derive(Clone, Debug)]
pub struct ThermalShockSph {
    /// Young's modulus E (Pa).
    pub young_modulus: f64,
    /// Poisson's ratio ν.
    pub poisson: f64,
    /// Thermal expansion coefficient α (1/K).
    pub alpha_expansion: f64,
    /// Fracture toughness K_IC (Pa·√m).
    pub fracture_toughness: f64,
    /// Critical tensile stress σ_c (Pa).
    pub critical_stress: f64,
    /// Reference temperature T_ref (K).
    pub t_ref: f64,
}
impl ThermalShockSph {
    /// Create a thermal shock model.
    pub fn new(
        young_modulus: f64,
        poisson: f64,
        alpha_expansion: f64,
        fracture_toughness: f64,
        critical_stress: f64,
        t_ref: f64,
    ) -> Self {
        Self {
            young_modulus,
            poisson,
            alpha_expansion,
            fracture_toughness,
            critical_stress,
            t_ref,
        }
    }
    /// Thermal stress σ_th = -E α ΔT / (1 - ν) (plane stress biaxial).
    pub fn thermal_stress(&self, temp: f64) -> f64 {
        let dt = temp - self.t_ref;
        let denom = 1.0 - self.poisson;
        if denom.abs() < 1e-300 {
            return 0.0;
        }
        -self.young_modulus * self.alpha_expansion * dt / denom
    }
    /// Thermal stress intensity factor K_I = σ_th √(π a) for a crack of half-length a.
    pub fn stress_intensity(&self, temp: f64, crack_half_length: f64) -> f64 {
        let sigma = self.thermal_stress(temp).abs();
        sigma * (PI * crack_half_length).sqrt()
    }
    /// Is fracture initiated at this temperature for a given crack size?
    pub fn is_fractured(&self, temp: f64, crack_half_length: f64) -> bool {
        self.stress_intensity(temp, crack_half_length) >= self.fracture_toughness
    }
    /// Biot modulus M_B = (1/E α²/(1-2ν)) characterizing thermal coupling strength.
    pub fn biot_modulus(&self) -> f64 {
        let denom = self.young_modulus * self.alpha_expansion * self.alpha_expansion;
        let nu_term = 1.0 - 2.0 * self.poisson;
        if denom < 1e-300 || nu_term.abs() < 1e-300 {
            return f64::INFINITY;
        }
        nu_term / denom
    }
    /// Minimum critical temperature difference for fracture initiation.
    pub fn critical_delta_t(&self) -> f64 {
        let denom = self.young_modulus * self.alpha_expansion;
        let nu_factor = 1.0 - self.poisson;
        if denom < 1e-300 || nu_factor.abs() < 1e-300 {
            return f64::INFINITY;
        }
        self.critical_stress * nu_factor / denom
    }
    /// Damage parameter D ∈ \[0, 1\] from thermal stress (linear softening).
    pub fn damage(&self, temp: f64) -> f64 {
        let sigma = self.thermal_stress(temp).abs();
        if sigma < self.critical_stress {
            return 0.0;
        }
        ((sigma - self.critical_stress) / self.critical_stress).min(1.0)
    }
}
/// SPH particle with thermal quantities.
#[derive(Clone, Debug)]
pub struct ThermalParticle {
    /// Position (m).
    pub pos: [f64; 3],
    /// Velocity (m/s).
    pub vel: [f64; 3],
    /// Mass (kg).
    pub mass: f64,
    /// Density (kg/m³).
    pub density: f64,
    /// Pressure (Pa).
    pub pressure: f64,
    /// Temperature (K).
    pub temperature: f64,
    /// Specific heat capacity c_p (J/kg/K).
    pub specific_heat: f64,
    /// Thermal conductivity λ (W/m/K).
    pub conductivity: f64,
    /// Heat flux magnitude (W/m²).
    pub heat_flux: f64,
    /// SPH smoothing length (m).
    pub smoothing_length: f64,
}
impl ThermalParticle {
    /// Construct a thermal particle.
    pub fn new(pos: [f64; 3], vel: [f64; 3], mass: f64, density: f64, temperature: f64) -> Self {
        ThermalParticle {
            pos,
            vel,
            mass,
            density,
            pressure: 0.0,
            temperature,
            specific_heat: 4182.0,
            conductivity: 0.6,
            heat_flux: 0.0,
            smoothing_length: 0.01,
        }
    }
    /// Enthalpy per unit volume (J/m³).
    pub fn enthalpy(&self) -> f64 {
        self.density * self.specific_heat * self.temperature
    }
    /// Thermal diffusivity α = λ / (ρ c_p).
    pub fn thermal_diffusivity(&self) -> f64 {
        if self.density < 1e-300 || self.specific_heat < 1e-300 {
            return 0.0;
        }
        self.conductivity / (self.density * self.specific_heat)
    }
    /// Update temperature from enthalpy change ΔH.
    pub fn update_temperature_from_dh(&mut self, dh: f64) {
        let cp_rho = self.specific_heat * self.density;
        if cp_rho > 1e-300 {
            self.temperature += dh / cp_rho;
        }
    }
}
/// SPH heat conduction using the Brookshaw consistent Laplacian operator.
pub struct SphHeatConduction {
    /// Thermal diffusivity α (m²/s).
    pub thermal_diffusivity: f64,
    /// SPH smoothing length (m).
    pub h: f64,
}
impl SphHeatConduction {
    /// Construct a heat conduction model.
    pub fn new(thermal_diffusivity: f64, h: f64) -> Self {
        SphHeatConduction {
            thermal_diffusivity,
            h,
        }
    }
    /// Compute the temperature rate dT/dt for particle i due to conduction.
    ///
    /// `neighbors`: slice of (r_ij, T_j, m_j, ρ_j) tuples.
    pub fn dtemp_dt(&self, t_i: f64, neighbors: &[([f64; 3], f64, f64, f64)]) -> f64 {
        let r_ij_list: Vec<[f64; 3]> = neighbors.iter().map(|(r, _, _, _)| *r).collect();
        let dt_list: Vec<f64> = neighbors.iter().map(|(_, tj, _, _)| t_i - tj).collect();
        let m_list: Vec<f64> = neighbors.iter().map(|(_, _, m, _)| *m).collect();
        let rho_list: Vec<f64> = neighbors.iter().map(|(_, _, _, rho)| *rho).collect();
        let lap = brookshaw_laplacian(&r_ij_list, &dt_list, &m_list, &rho_list, self.h);
        self.thermal_diffusivity * lap
    }
    /// Explicit Euler step: T^{n+1} = T^n + dt * α ∇²T.
    pub fn step(&self, t_i: f64, neighbors: &[([f64; 3], f64, f64, f64)], dt: f64) -> f64 {
        t_i + dt * self.dtemp_dt(t_i, neighbors)
    }
}
/// Thermal expansion model: updates density and pressure.
pub struct ThermalExpansion {
    /// Thermal expansion coefficient β (1/K).
    pub beta: f64,
    /// Reference temperature T_0 (K).
    pub t_ref: f64,
    /// Reference density ρ_0 (kg/m³).
    pub rho_ref: f64,
    /// Bulk modulus K (Pa) for pressure correction.
    pub bulk_modulus: f64,
}
impl ThermalExpansion {
    /// Construct a thermal expansion model.
    pub fn new(beta: f64, t_ref: f64, rho_ref: f64, bulk_modulus: f64) -> Self {
        ThermalExpansion {
            beta,
            t_ref,
            rho_ref,
            bulk_modulus,
        }
    }
    /// Updated density: ρ = ρ_0 * (1 − β ΔT).
    pub fn density(&self, temp: f64) -> f64 {
        let dt = temp - self.t_ref;
        self.rho_ref * (1.0 - self.beta * dt)
    }
    /// Pressure change due to thermal expansion: Δp = K β ΔT.
    pub fn pressure_change(&self, temp: f64) -> f64 {
        let dt = temp - self.t_ref;
        self.bulk_modulus * self.beta * dt
    }
    /// Boussinesq buoyancy force per unit mass (upward = +z).
    pub fn boussinesq_force(&self, temp: f64, g: f64) -> [f64; 3] {
        let dt = temp - self.t_ref;
        [0.0, 0.0, self.beta * dt * g]
    }
}
/// SPH combustion model: fuel-oxidizer mixing with Arrhenius reaction rate.
pub struct CombustionSph {
    /// Pre-exponential factor A (1/s).
    pub arrhenius_a: f64,
    /// Activation energy E_a (J/mol).
    pub activation_energy: f64,
    /// Universal gas constant R (J/mol/K).
    pub gas_constant: f64,
    /// Heat of combustion Q (J/kg_fuel).
    pub heat_of_combustion: f64,
}
impl CombustionSph {
    /// Construct a combustion model.
    pub fn new(arrhenius_a: f64, activation_energy: f64, heat_of_combustion: f64) -> Self {
        CombustionSph {
            arrhenius_a,
            activation_energy,
            gas_constant: 8.314,
            heat_of_combustion,
        }
    }
    /// Arrhenius reaction rate k(T) = A exp(−E_a / (R T)).
    pub fn reaction_rate(&self, temp: f64) -> f64 {
        if temp < 1e-300 {
            return 0.0;
        }
        self.arrhenius_a * (-self.activation_energy / (self.gas_constant * temp)).exp()
    }
    /// Mass consumption rate of fuel ω_f = k(T) * Y_fuel * ρ.
    pub fn fuel_consumption_rate(&self, temp: f64, y_fuel: f64, density: f64) -> f64 {
        self.reaction_rate(temp) * y_fuel * density
    }
    /// Heat release rate q = ω_f * Q.
    pub fn heat_release_rate(&self, temp: f64, y_fuel: f64, density: f64) -> f64 {
        self.fuel_consumption_rate(temp, y_fuel, density) * self.heat_of_combustion
    }
    /// Adiabatic flame temperature (simplified): T_f = T_0 + Q * Y_f / c_p.
    pub fn adiabatic_flame_temp(&self, t_0: f64, y_fuel: f64, c_p: f64) -> f64 {
        if c_p < 1e-300 {
            return t_0;
        }
        t_0 + self.heat_of_combustion * y_fuel / c_p
    }
}
