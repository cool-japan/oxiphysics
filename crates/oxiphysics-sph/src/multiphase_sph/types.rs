//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use std::f64::consts::PI;

use super::functions::*;

/// Evaporation and condensation source terms for SPH multiphase flow.
///
/// Implements the Hertz-Knudsen mass-transfer rate at the liquid-vapour
/// interface and a simplified Rankine-Hugoniot jump across the phase boundary.
#[derive(Debug, Clone)]
pub struct PhaseTransitionSph {
    /// Latent heat of vaporisation L \[J/kg\].
    pub latent_heat: f64,
    /// Evaporation coefficient β_evap ∈ \[0, 1\].
    pub evaporation_coeff: f64,
    /// Condensation coefficient β_cond ∈ \[0, 1\].
    pub condensation_coeff: f64,
    /// Molecular weight M_w \[kg/mol\].
    pub molar_mass: f64,
    /// Saturation temperature T_sat \[K\].
    pub t_sat: f64,
    /// Saturation pressure p_sat \[Pa\].
    pub p_sat: f64,
}
impl PhaseTransitionSph {
    /// Create a new `PhaseTransitionSph` for water at 100 °C.
    pub fn water_steam() -> Self {
        Self {
            latent_heat: 2.257e6,
            evaporation_coeff: 0.1,
            condensation_coeff: 0.1,
            molar_mass: 0.018,
            t_sat: 373.15,
            p_sat: 101325.0,
        }
    }
    /// Create a new `PhaseTransitionSph` with custom parameters.
    pub fn new(
        latent_heat: f64,
        evaporation_coeff: f64,
        condensation_coeff: f64,
        molar_mass: f64,
        t_sat: f64,
        p_sat: f64,
    ) -> Self {
        Self {
            latent_heat,
            evaporation_coeff,
            condensation_coeff,
            molar_mass,
            t_sat,
            p_sat,
        }
    }
    /// Hertz-Knudsen evaporation mass flux \[kg/(m² s)\].
    ///
    /// ṁ_ev = β_e · sqrt(M / (2π R T_sat)) · (p_sat − p_v)
    /// with p_v the local vapour partial pressure.
    pub fn evaporation_flux(&self, temperature: f64, vapor_pressure: f64) -> f64 {
        let r_gas = 8.314_f64;
        let t = temperature.max(1.0_f64);
        let coeff = self.evaporation_coeff * (self.molar_mass / (2.0 * PI * r_gas * t)).sqrt();
        let dp = self.p_sat - vapor_pressure;
        coeff * dp.max(0.0)
    }
    /// Hertz-Knudsen condensation mass flux \[kg/(m² s)\].
    pub fn condensation_flux(&self, temperature: f64, vapor_pressure: f64) -> f64 {
        let r_gas = 8.314_f64;
        let t = temperature.max(1.0_f64);
        let coeff = self.condensation_coeff * (self.molar_mass / (2.0 * PI * r_gas * t)).sqrt();
        let dp = vapor_pressure - self.p_sat;
        coeff * dp.max(0.0)
    }
    /// Net mass-transfer rate (positive = evaporation, negative = condensation).
    pub fn net_mass_flux(&self, temperature: f64, vapor_pressure: f64) -> f64 {
        self.evaporation_flux(temperature, vapor_pressure)
            - self.condensation_flux(temperature, vapor_pressure)
    }
    /// Clausius-Clapeyron saturation pressure at temperature T \[K\].
    ///
    /// p_sat(T) ≈ p_ref · exp(-L M / R · (1/T - 1/T_ref))
    pub fn clausius_clapeyron_pressure(&self, temperature: f64) -> f64 {
        let r_gas = 8.314_f64;
        let exponent = -self.latent_heat * self.molar_mass / r_gas
            * (1.0 / temperature.max(1.0_f64) - 1.0 / self.t_sat);
        self.p_sat * exponent.exp()
    }
    /// Thermal energy sink due to evaporation: Q = -ṁ · L \[W/m²\].
    pub fn energy_sink(&self, mass_flux: f64) -> f64 {
        -mass_flux * self.latent_heat
    }
    /// Jakob number Ja = ρ_l c_pl ΔT / (ρ_v L).
    ///
    /// Used to characterise the superheat driving bubble growth.
    pub fn jakob_number(&self, rho_l: f64, cp_l: f64, superheat: f64, rho_v: f64) -> f64 {
        rho_l * cp_l * superheat / (rho_v * self.latent_heat).max(1e-30)
    }
}
/// Full-state SPH particle for multiphase simulations.
///
/// Carries phase label, density, pressure, velocity, surface-tension force,
/// and the smoothed colour-function value used by the CSF method.
#[derive(Debug, Clone)]
pub struct MultiphaseSphParticle {
    /// Position \[m\].
    pub position: [f64; 3],
    /// Velocity \[m/s\].
    pub velocity: [f64; 3],
    /// Particle mass \[kg\].
    pub mass: f64,
    /// SPH density estimate ρ \[kg/m³\].
    pub density: f64,
    /// Pressure p \[Pa\].
    pub pressure: f64,
    /// Phase label (0-indexed).
    pub phase_label: usize,
    /// Smoothed colour function C ∈ \[0, 1\] (1 = phase 0, 0 = phase 1).
    pub color_function: f64,
    /// Surface-tension body force f_st \[N/kg\] (acceleration).
    pub surface_tension_force: [f64; 3],
    /// Interface normal n̂ (from colour gradient).
    pub interface_normal: [f64; 3],
    /// Interface curvature κ \[1/m\].
    pub curvature: f64,
}
impl MultiphaseSphParticle {
    /// Create a resting particle belonging to `phase_label`.
    pub fn new(position: [f64; 3], phase_label: usize, mass: f64, density: f64) -> Self {
        let color_function = if phase_label == 0 { 1.0 } else { 0.0 };
        Self {
            position,
            velocity: [0.0; 3],
            mass,
            density,
            pressure: 0.0,
            phase_label,
            color_function,
            surface_tension_force: [0.0; 3],
            interface_normal: [0.0, 0.0, 1.0],
            curvature: 0.0,
        }
    }
    /// Squared distance to another particle.
    pub fn dist2_to(&self, other: &Self) -> f64 {
        let dx = self.position[0] - other.position[0];
        let dy = self.position[1] - other.position[1];
        let dz = self.position[2] - other.position[2];
        dx * dx + dy * dy + dz * dz
    }
    /// Distance to another particle.
    pub fn dist_to(&self, other: &Self) -> f64 {
        self.dist2_to(other).sqrt()
    }
    /// Kinetic energy ½ m v².
    pub fn kinetic_energy(&self) -> f64 {
        let v2: f64 = self.velocity.iter().map(|v| v * v).sum();
        0.5 * self.mass * v2
    }
    /// Returns true if this particle is at the diffuse interface (0.1 < C < 0.9).
    pub fn is_interface(&self) -> bool {
        self.color_function > 0.1 && self.color_function < 0.9
    }
    /// Linear interpolation between phase-0 property `p0` and phase-1 property `p1`.
    pub fn interpolate(&self, p0: f64, p1: f64) -> f64 {
        self.color_function * p0 + (1.0 - self.color_function) * p1
    }
}
/// Volume-of-fluid (VOF) fraction transport for SPH.
///
/// Advects the volume fraction α using SPH interpolation.
#[derive(Debug, Clone)]
pub struct VolumeFractionField {
    /// Volume fraction α at each node ∈ \[0,1\]
    pub alpha: Vec<f64>,
    /// Node positions \[m\]
    pub positions: Vec<f64>,
    /// Node velocities \[m/s\]
    pub velocities: Vec<f64>,
    /// Smoothing length h \[m\]
    pub h: f64,
    /// Number of nodes
    pub n: usize,
}
impl VolumeFractionField {
    /// Create a new 1D volume fraction field.
    pub fn new_1d(n: usize, length: f64, alpha_init: f64) -> Self {
        let h = 2.0 * length / n as f64;
        let positions = (0..n)
            .map(|i| i as f64 * length / (n as f64 - 1.0))
            .collect();
        let velocities = vec![0.0; n];
        let alpha = vec![alpha_init; n];
        Self {
            alpha,
            positions,
            velocities,
            h,
            n,
        }
    }
    /// Initialize with a step function: α = 1 for x < x_interface, α = 0 otherwise.
    pub fn step_init(&mut self, x_interface: f64) {
        for (i, &x) in self.positions.iter().enumerate() {
            self.alpha[i] = if x < x_interface { 1.0 } else { 0.0 };
        }
    }
    /// Advect the volume fraction one time step (upwind scheme).
    pub fn advect_step(&mut self, dt: f64) {
        let dx = if self.n > 1 {
            self.positions[1] - self.positions[0]
        } else {
            1.0
        };
        let alpha_old = self.alpha.clone();
        for i in 1..self.n - 1 {
            let v = self.velocities[i];
            let flux = if v >= 0.0 {
                v * alpha_old[i - 1]
            } else {
                v * alpha_old[i + 1]
            };
            self.alpha[i] = alpha_old[i] - dt / dx * (v * alpha_old[i] - flux);
            self.alpha[i] = self.alpha[i].clamp(0.0, 1.0);
        }
    }
    /// Total volume (integral of α over domain).
    pub fn total_volume(&self) -> f64 {
        let dx = if self.n > 1 {
            self.positions[1] - self.positions[0]
        } else {
            1.0
        };
        self.alpha.iter().sum::<f64>() * dx
    }
    /// Interface location from gradient of α.
    pub fn interface_location(&self) -> Option<f64> {
        for i in 0..self.n - 1 {
            if (self.alpha[i] - 0.5) * (self.alpha[i + 1] - 0.5) < 0.0 {
                let x_i = self.positions[i];
                let frac = (0.5 - self.alpha[i]) / (self.alpha[i + 1] - self.alpha[i]);
                return Some(x_i + frac * (self.positions[i + 1] - x_i));
            }
        }
        None
    }
    /// Set uniform velocity.
    pub fn set_uniform_velocity(&mut self, v: f64) {
        self.velocities.fill(v);
    }
}
/// Rayleigh-Taylor instability analysis.
///
/// A heavy fluid (ρ₁) sits on top of a light fluid (ρ₂) under gravity.
#[derive(Debug, Clone)]
pub struct RayleighTaylorInstability {
    /// Heavy fluid density ρ₁ \[kg/m³\]
    pub rho_heavy: f64,
    /// Light fluid density ρ₂ \[kg/m³\]
    pub rho_light: f64,
    /// Gravitational acceleration g \[m/s²\]
    pub g: f64,
    /// Surface tension σ \[N/m\]
    pub sigma: f64,
    /// Kinematic viscosity of heavy fluid ν₁ \[m²/s\]
    pub nu_heavy: f64,
    /// Kinematic viscosity of light fluid ν₂ \[m²/s\]
    pub nu_light: f64,
}
impl RayleighTaylorInstability {
    /// Atwood number At = (ρ₁ − ρ₂) / (ρ₁ + ρ₂).
    pub fn atwood_number(&self) -> f64 {
        (self.rho_heavy - self.rho_light) / (self.rho_heavy + self.rho_light).max(1e-30)
    }
    /// Linear growth rate for wavenumber k (inviscid, no surface tension).
    ///
    /// σ² = At · g · k
    pub fn growth_rate_inviscid(&self, k: f64) -> f64 {
        (self.atwood_number() * self.g * k).sqrt()
    }
    /// Growth rate with surface tension (modified dispersion).
    ///
    /// σ² = At · g · k − σ · k³ / (ρ₁ + ρ₂)
    pub fn growth_rate_with_tension(&self, k: f64) -> f64 {
        let at = self.atwood_number();
        let rho_sum = self.rho_heavy + self.rho_light;
        let val = at * self.g * k - self.sigma * k * k * k / rho_sum.max(1e-30);
        if val > 0.0 { val.sqrt() } else { 0.0 }
    }
    /// Critical wavenumber above which surface tension stabilizes the interface.
    ///
    /// k_c = sqrt(g·(ρ₁−ρ₂) / σ)
    pub fn critical_wavenumber(&self) -> f64 {
        let drho = self.rho_heavy - self.rho_light;
        if drho > 0.0 && self.sigma > 0.0 {
            (self.g * drho / self.sigma).sqrt()
        } else {
            f64::INFINITY
        }
    }
    /// Most dangerous wavenumber k_max = k_c / sqrt(3).
    pub fn most_dangerous_wavenumber(&self) -> f64 {
        self.critical_wavenumber() / 3.0_f64.sqrt()
    }
    /// Bubble rise velocity (Davies-Taylor, 3D):
    ///
    /// V_b = 0.23 · sqrt(g · R · At)  where R = 2π/k is the bubble radius
    pub fn bubble_rise_velocity(&self, k: f64) -> f64 {
        let at = self.atwood_number();
        let r = PI / k;
        0.23 * (self.g * r * at).sqrt()
    }
    /// Mushroom stem velocity (Dimonte):
    pub fn spike_velocity(&self, t: f64) -> f64 {
        let at = self.atwood_number();
        let k0 = self.most_dangerous_wavenumber();
        let sigma0 = self.growth_rate_inviscid(k0);
        let v_sat = (2.0 * at * self.g / k0).sqrt();
        v_sat * (1.0 - (-2.0 * sigma0 * t).exp())
    }
}
/// A full-state multiphase SPH particle with density, pressure, phase tag, and
/// colour field.
#[derive(Debug, Clone)]
pub struct MultiphaseFluidParticle {
    /// Position \[m\].
    pub position: [f64; 3],
    /// Velocity \[m/s\].
    pub velocity: [f64; 3],
    /// Particle mass \[kg\].
    pub mass: f64,
    /// SPH density estimate \[kg/m³\].
    pub density: f64,
    /// Pressure \[Pa\].
    pub pressure: f64,
    /// Phase identifier (0 = fluid A, 1 = fluid B, …).
    pub phase: u32,
    /// Colour-function value C ∈ \[0, 1\].
    pub color_field: f64,
}
impl MultiphaseFluidParticle {
    /// Create a new resting particle.
    pub fn new(position: [f64; 3], mass: f64, density: f64, phase: u32) -> Self {
        let color_field = if phase == 0 { 1.0 } else { 0.0 };
        Self {
            position,
            velocity: [0.0; 3],
            mass,
            density,
            pressure: 0.0,
            phase,
            color_field,
        }
    }
}
/// High-level multiphase SPH simulation state.
#[derive(Debug, Clone)]
pub struct MultiphaseSphSimulation {
    /// All particles (both phases)
    pub particles: Vec<ColorParticle>,
    /// Surface tension model
    pub surface_tension: CsfSurfaceTension,
    /// Gravity \[m/s²\] (signed, positive upward = negative)
    pub gravity: [f64; 2],
    /// Simulation time \[s\]
    pub time: f64,
    /// Time step \[s\]
    pub dt: f64,
    /// Smoothing length h \[m\]
    pub h: f64,
    /// Liquid density \[kg/m³\]
    pub rho_l: f64,
    /// Gas density \[kg/m³\]
    pub rho_g: f64,
    /// Step count
    pub step_count: usize,
}
impl MultiphaseSphSimulation {
    /// Create a new multiphase SPH simulation.
    pub fn new(h: f64, sigma: f64, rho_l: f64, rho_g: f64) -> Self {
        Self {
            particles: Vec::new(),
            surface_tension: CsfSurfaceTension::new(sigma, h),
            gravity: [0.0, -9.81],
            time: 0.0,
            dt: 1e-4,
            h,
            rho_l,
            rho_g,
            step_count: 0,
        }
    }
    /// Add a liquid particle.
    pub fn add_liquid_particle(&mut self, x: f64, y: f64, mass: f64) {
        self.particles.push(ColorParticle::new_liquid(x, y, mass));
    }
    /// Add a gas particle.
    pub fn add_gas_particle(&mut self, x: f64, y: f64, mass: f64) {
        self.particles.push(ColorParticle::new_gas(x, y, mass));
    }
    /// Count liquid particles.
    pub fn n_liquid(&self) -> usize {
        self.particles.iter().filter(|p| p.color > 0.5).count()
    }
    /// Count gas particles.
    pub fn n_gas(&self) -> usize {
        self.particles.iter().filter(|p| p.color <= 0.5).count()
    }
    /// Compute and update color gradients for all particles.
    pub fn update_color_gradients(&mut self) {
        compute_color_gradients(&mut self.particles, self.h);
    }
    /// Simple time step: advect particles with their velocities.
    pub fn step(&mut self) {
        self.update_color_gradients();
        let n = self.particles.len();
        for i in 0..n {
            let color = self.particles[i].color;
            let rho = self.particles[i].interpolated_density(self.rho_l, self.rho_g);
            let fg = [self.gravity[0] * rho, self.gravity[1] * rho];
            let kappa = self.particles[i].curvature;
            let normal = self.particles[i].normal;
            let phi = color - 0.5;
            let fs = self.surface_tension.csf_force(kappa, normal, phi);
            let ax = (fg[0] + fs[0]) / rho.max(1e-15);
            let ay = (fg[1] + fs[1]) / rho.max(1e-15);
            self.particles[i].vel[0] += ax * self.dt;
            self.particles[i].vel[1] += ay * self.dt;
            self.particles[i].pos[0] += self.particles[i].vel[0] * self.dt;
            self.particles[i].pos[1] += self.particles[i].vel[1] * self.dt;
        }
        self.time += self.dt;
        self.step_count += 1;
    }
    /// Total kinetic energy.
    pub fn kinetic_energy(&self) -> f64 {
        self.particles
            .iter()
            .map(|p| {
                let v2 = p.vel[0].powi(2) + p.vel[1].powi(2);
                0.5 * p.mass * v2
            })
            .sum()
    }
    /// Center of mass of liquid phase.
    pub fn liquid_center_of_mass(&self) -> [f64; 2] {
        let liquids: Vec<&ColorParticle> =
            self.particles.iter().filter(|p| p.color > 0.5).collect();
        if liquids.is_empty() {
            return [0.0; 2];
        }
        let total_mass: f64 = liquids.iter().map(|p| p.mass).sum();
        let cx: f64 = liquids.iter().map(|p| p.pos[0] * p.mass).sum::<f64>() / total_mass;
        let cy: f64 = liquids.iter().map(|p| p.pos[1] * p.mass).sum::<f64>() / total_mass;
        [cx, cy]
    }
}
/// Continuum Surface Force (CSF) method for SPH interface tracking.
///
/// Computes the smoothed colour function, its gradient (interface normal),
/// and the divergence of the normal (curvature κ = −∇·n̂).
///
/// Reference: Brackbill, Kothe & Zemach, JCP 100, 335–354 (1992).
#[derive(Debug, Clone)]
pub struct ColorFunctionSph {
    /// Surface-tension coefficient σ \[N/m\].
    pub sigma: f64,
    /// Kernel smoothing length h \[m\].
    pub h: f64,
}
impl ColorFunctionSph {
    /// Create a new `ColorFunctionSph`.
    pub fn new(sigma: f64, h: f64) -> Self {
        Self { sigma, h }
    }
    /// Compute smoothed colour function C_i via SPH kernel summation.
    ///
    /// C_i = Σ_j (m_j / ρ_j) · C_j · W(|r_ij|, h)
    pub fn smooth_color(&self, particles: &[MultiphaseSphParticle], i: usize) -> f64 {
        let xi = particles[i].position;
        let mut c_sum = 0.0_f64;
        let mut w_sum = 0.0_f64;
        for p in particles {
            let dx = xi[0] - p.position[0];
            let dy = xi[1] - p.position[1];
            let dz = xi[2] - p.position[2];
            let r = (dx * dx + dy * dy + dz * dz).sqrt();
            if r < self.h {
                let w = cubic_spline_kernel(r, self.h);
                let vol = p.mass / p.density.max(1e-30);
                c_sum += vol * p.color_function * w;
                w_sum += vol * w;
            }
        }
        if w_sum > 1e-30 {
            c_sum / w_sum
        } else {
            particles[i].color_function
        }
    }
    /// Compute colour-function gradient ∇C_i (unnormalized interface normal).
    pub fn color_gradient(&self, particles: &[MultiphaseSphParticle], i: usize) -> [f64; 3] {
        let xi = particles[i].position;
        let mut grad = [0.0_f64; 3];
        for p in particles {
            let rij = [
                xi[0] - p.position[0],
                xi[1] - p.position[1],
                xi[2] - p.position[2],
            ];
            let r = (rij[0] * rij[0] + rij[1] * rij[1] + rij[2] * rij[2]).sqrt();
            if r < self.h && r > 1e-14 {
                let dw_dr = cubic_spline_kernel_grad(r, self.h);
                let vol = p.mass / p.density.max(1e-30);
                let dc = p.color_function - particles[i].color_function;
                for d in 0..3 {
                    grad[d] += vol * dc * dw_dr * rij[d] / r;
                }
            }
        }
        grad
    }
    /// CSF surface-tension body force on particle `i`.
    ///
    /// f_st = σ κ_i ∇C_i / ρ_i
    pub fn csf_body_force(&self, particles: &[MultiphaseSphParticle], i: usize) -> [f64; 3] {
        let grad_c = self.color_gradient(particles, i);
        let grad_c_mag = (grad_c[0].powi(2) + grad_c[1].powi(2) + grad_c[2].powi(2)).sqrt();
        if grad_c_mag < 1e-14 {
            return [0.0; 3];
        }
        let n_hat: [f64; 3] = [
            grad_c[0] / grad_c_mag,
            grad_c[1] / grad_c_mag,
            grad_c[2] / grad_c_mag,
        ];
        let kappa = particles[i].curvature;
        let rho_i = particles[i].density.max(1e-30);
        let coeff = self.sigma * kappa * grad_c_mag / rho_i;
        [coeff * n_hat[0], coeff * n_hat[1], coeff * n_hat[2]]
    }
    /// Laplace pressure jump ΔP = σ κ at a spherical interface of radius R.
    pub fn laplace_pressure(&self, radius: f64) -> f64 {
        2.0 * self.sigma / radius.max(1e-30)
    }
}
/// Color function particle for multiphase SPH.
///
/// Each particle carries a color value C ∈ \[0,1\] indicating phase membership.
#[derive(Debug, Clone)]
pub struct ColorParticle {
    /// Position \[m\]
    pub pos: [f64; 2],
    /// Velocity \[m/s\]
    pub vel: [f64; 2],
    /// Color function value (0 = gas, 1 = liquid)
    pub color: f64,
    /// Mass \[kg\]
    pub mass: f64,
    /// Density \[kg/m³\]
    pub density: f64,
    /// Smoothed color function gradient
    pub color_grad: [f64; 2],
    /// Normal vector (from color gradient)
    pub normal: [f64; 2],
    /// Curvature κ
    pub curvature: f64,
    /// Pressure \[Pa\]
    pub pressure: f64,
    /// Nearest interface distance (signed)
    pub dist_to_interface: f64,
}
impl ColorParticle {
    /// Create a new liquid particle.
    pub fn new_liquid(x: f64, y: f64, mass: f64) -> Self {
        Self {
            pos: [x, y],
            vel: [0.0; 2],
            color: 1.0,
            mass,
            density: 1000.0,
            color_grad: [0.0; 2],
            normal: [0.0, 1.0],
            curvature: 0.0,
            pressure: 0.0,
            dist_to_interface: f64::MAX,
        }
    }
    /// Create a new gas particle.
    pub fn new_gas(x: f64, y: f64, mass: f64) -> Self {
        Self {
            pos: [x, y],
            vel: [0.0; 2],
            color: 0.0,
            mass,
            density: 1.225,
            color_grad: [0.0; 2],
            normal: [0.0, -1.0],
            curvature: 0.0,
            pressure: 0.0,
            dist_to_interface: f64::MAX,
        }
    }
    /// Returns true if this particle is at the interface (0.1 < C < 0.9).
    #[inline]
    pub fn is_interface(&self) -> bool {
        self.color > 0.1 && self.color < 0.9
    }
    /// Phase of this particle based on color value.
    pub fn phase(&self) -> Phase {
        if self.color > 0.5 {
            Phase::Liquid
        } else {
            Phase::Gas
        }
    }
    /// Update the normal from color gradient.
    pub fn update_normal(&mut self) {
        let norm = (self.color_grad[0].powi(2) + self.color_grad[1].powi(2)).sqrt();
        if norm > 1e-14 {
            self.normal[0] = self.color_grad[0] / norm;
            self.normal[1] = self.color_grad[1] / norm;
        }
    }
    /// Interpolated density from color value.
    pub fn interpolated_density(&self, rho_l: f64, rho_g: f64) -> f64 {
        self.color * rho_l + (1.0 - self.color) * rho_g
    }
    /// Interpolated viscosity from color value.
    pub fn interpolated_viscosity(&self, mu_l: f64, mu_g: f64) -> f64 {
        self.color * mu_l + (1.0 - self.color) * mu_g
    }
}
/// Rayleigh-Plesset equation for spherical bubble dynamics.
///
/// R·R̈ + (3/2)·Ṙ² = (p_b − p_∞) / ρ_l − 4μ_l·Ṙ/R − 2σ/(ρ_l·R)
#[derive(Debug, Clone)]
pub struct RayleighPlesset {
    /// Bubble radius R \[m\]
    pub radius: f64,
    /// Bubble wall velocity Ṙ \[m/s\]
    pub radius_dot: f64,
    /// Liquid density ρ_l \[kg/m³\]
    pub rho_l: f64,
    /// Liquid viscosity μ_l \[Pa·s\]
    pub mu_l: f64,
    /// Surface tension σ \[N/m\]
    pub sigma: f64,
    /// Gas pressure inside bubble p_b \[Pa\]
    pub p_bubble: f64,
    /// Far-field pressure p_∞ \[Pa\]
    pub p_inf: f64,
    /// Initial radius R₀ \[m\]
    pub r0: f64,
    /// Initial gas pressure p_b0 \[Pa\]
    pub p_b0: f64,
    /// Polytropic index γ
    pub gamma: f64,
}
impl RayleighPlesset {
    /// Create a new Rayleigh-Plesset bubble.
    pub fn new(r0: f64, p_inf: f64, rho_l: f64, mu_l: f64, sigma: f64) -> Self {
        let p_b0 = p_inf + 2.0 * sigma / r0;
        Self {
            radius: r0,
            radius_dot: 0.0,
            rho_l,
            mu_l,
            sigma,
            p_bubble: p_b0,
            p_inf,
            r0,
            p_b0,
            gamma: 1.4,
        }
    }
    /// Update gas pressure using adiabatic law: p_b = p_b0 · (R₀/R)^(3γ).
    pub fn update_gas_pressure(&mut self) {
        self.p_bubble = self.p_b0 * (self.r0 / self.radius.max(1e-10)).powf(3.0 * self.gamma);
    }
    /// Compute R̈ from the Rayleigh-Plesset equation.
    pub fn radius_ddot(&self) -> f64 {
        let r = self.radius.max(1e-10);
        let rho = self.rho_l;
        let p_drive = (self.p_bubble - self.p_inf) / rho;
        let viscous = 4.0 * self.mu_l * self.radius_dot / (rho * r);
        let tension = 2.0 * self.sigma / (rho * r);
        let inertial = 1.5 * self.radius_dot * self.radius_dot / r;
        (p_drive - viscous - tension) / r - inertial
    }
    /// One explicit time step (second-order Verlet).
    pub fn step(&mut self, dt: f64) {
        self.update_gas_pressure();
        let r_ddot = self.radius_ddot();
        self.radius += self.radius_dot * dt + 0.5 * r_ddot * dt * dt;
        self.radius = self.radius.max(1e-12);
        self.update_gas_pressure();
        let r_ddot_new = self.radius_ddot();
        self.radius_dot += 0.5 * (r_ddot + r_ddot_new) * dt;
    }
    /// Natural frequency of small-amplitude oscillation (Minnaert frequency).
    ///
    /// ω_M = sqrt(3γ·p_inf / (ρ_l·R₀²))
    pub fn minnaert_frequency(&self) -> f64 {
        (3.0 * self.gamma * self.p_inf / (self.rho_l * self.r0 * self.r0)).sqrt()
    }
    /// Collapse time estimate (Rayleigh collapse time).
    ///
    /// t_c ≈ 0.915 · R₀ · sqrt(ρ_l / p_inf)
    pub fn collapse_time(&self) -> f64 {
        0.915 * self.r0 * (self.rho_l / self.p_inf.max(1e-10)).sqrt()
    }
    /// Radiated pressure at distance d from bubble center.
    ///
    /// p_rad = ρ_l / d · (R·R̈ + 2·Ṙ²) · R
    pub fn radiated_pressure(&self, distance: f64) -> f64 {
        let r = self.radius;
        let rd = self.radius_dot;
        let rdd = self.radius_ddot();
        self.rho_l * r * r / distance.max(1e-10) * (rdd + 2.0 * rd * rd / r)
    }
}
/// High-level multiphase SPH simulation driver.
///
/// Manages a collection of [`MultiphaseParticle`]s belonging to registered
/// [`PhaseDesc`] phases and advances the system by calling [`MultiphaseSPH::step`].
#[derive(Debug, Clone)]
pub struct MultiphaseSPH {
    /// All particles in the simulation.
    pub particles: Vec<MultiphaseParticle>,
    /// Registered phases (indexed by `phase_id`).
    pub phases: Vec<PhaseDesc>,
    /// Kernel smoothing length h \[m\].
    pub smoothing_length: f64,
    /// Interface model choice.
    pub interface_model: InterfaceModel,
    /// Surface-tension force formulation.
    pub surface_tension_force: SurfaceTensionForce,
    /// Gravity vector \[m/s²\].
    pub gravity: [f64; 3],
    /// Time-step size \[s\].
    pub dt: f64,
}
impl MultiphaseSPH {
    /// Create a new `MultiphaseSPH` simulation.
    ///
    /// # Arguments
    /// * `smoothing_length` – kernel support radius h \[m\]
    /// * `gravity`          – body-force acceleration \[m/s²\]
    /// * `dt`               – time-step \[s\]
    pub fn new(smoothing_length: f64, gravity: [f64; 3], dt: f64) -> Self {
        Self {
            particles: Vec::new(),
            phases: Vec::new(),
            smoothing_length,
            interface_model: InterfaceModel::ColorFunction,
            surface_tension_force: SurfaceTensionForce::Csf,
            gravity,
            dt,
        }
    }
    /// Register a new phase and return its id.
    pub fn add_phase(&mut self, density: f64, viscosity: f64, surface_tension_coeff: f64) -> usize {
        let id = self.phases.len();
        self.phases.push(PhaseDesc::new(
            id,
            density,
            viscosity,
            surface_tension_coeff,
        ));
        id
    }
    /// Add a particle at `position` belonging to `phase_id`.
    pub fn add_particle(&mut self, position: [f64; 3], phase_id: usize, mass: f64) {
        self.particles
            .push(MultiphaseParticle::new(position, phase_id, mass));
    }
    /// Total number of particles.
    pub fn n_particles(&self) -> usize {
        self.particles.len()
    }
    /// Update colour functions for all particles.
    pub fn update_color_functions(&mut self) {
        let colors = compute_color_function(&self.particles, self.smoothing_length);
        for (p, c) in self.particles.iter_mut().zip(colors.iter()) {
            p.color_function = *c;
        }
    }
    /// Update interface curvatures for all particles.
    pub fn update_curvatures(&mut self) {
        let colors: Vec<f64> = self.particles.iter().map(|p| p.color_function).collect();
        let kappas = compute_curvature(&colors);
        for (p, k) in self.particles.iter_mut().zip(kappas.iter()) {
            p.curvature = *k;
        }
    }
    /// Advance the simulation by one time step.
    ///
    /// Steps performed:
    /// 1. Update colour functions.
    /// 2. Compute curvatures.
    /// 3. Apply gravity and surface-tension forces.
    /// 4. Integrate positions with symplectic Euler.
    pub fn step(&mut self) {
        self.update_color_functions();
        self.update_curvatures();
        let n = self.particles.len();
        let h = self.smoothing_length;
        let g = self.gravity;
        let dt = self.dt;
        let mut accel = vec![[0.0_f64; 3]; n];
        for (i, acc_i) in accel.iter_mut().enumerate() {
            acc_i[0] += g[0];
            acc_i[1] += g[1];
            acc_i[2] += g[2];
            if self.surface_tension_force == SurfaceTensionForce::Csf {
                let phase = self.particles[i].phase_id;
                if phase < self.phases.len() {
                    let sigma = self.phases[phase].surface_tension_coeff;
                    let rho = self.phases[phase].density;
                    let kappa = self.particles[i].curvature;
                    let mut grad_c = [0.0_f64; 3];
                    for (j, pj) in self.particles.iter().enumerate() {
                        if i == j {
                            continue;
                        }
                        let dx = self.particles[i].position[0] - pj.position[0];
                        let dy = self.particles[i].position[1] - pj.position[1];
                        let dz = self.particles[i].position[2] - pj.position[2];
                        let r = (dx * dx + dy * dy + dz * dz).sqrt();
                        if r < h && r > 1e-14 {
                            let dc = pj.color_function - self.particles[i].color_function;
                            let wg = sph_kernel_gradient_magnitude(r, h) / r;
                            grad_c[0] += pj.mass * dc * wg * dx;
                            grad_c[1] += pj.mass * dc * wg * dy;
                            grad_c[2] += pj.mass * dc * wg * dz;
                        }
                    }
                    let factor = sigma * kappa / rho;
                    acc_i[0] += factor * grad_c[0];
                    acc_i[1] += factor * grad_c[1];
                    acc_i[2] += factor * grad_c[2];
                }
            }
        }
        for (i, acc_i) in accel.iter().enumerate() {
            self.particles[i].velocity[0] += acc_i[0] * dt;
            self.particles[i].velocity[1] += acc_i[1] * dt;
            self.particles[i].velocity[2] += acc_i[2] * dt;
            self.particles[i].position[0] += self.particles[i].velocity[0] * dt;
            self.particles[i].position[1] += self.particles[i].velocity[1] * dt;
            self.particles[i].position[2] += self.particles[i].velocity[2] * dt;
        }
    }
    /// Total kinetic energy of all particles \[J\].
    pub fn kinetic_energy(&self) -> f64 {
        self.particles.iter().map(|p| p.kinetic_energy()).sum()
    }
}
/// Parameters for a two-phase immiscible SPH simulation.
#[derive(Debug, Clone)]
pub struct MultiphaseParams {
    /// Surface-tension coefficient σ \[N/m\].
    pub surface_tension_coeff: f64,
    /// Diffuse-interface half-width ε \[m\].
    pub interface_width: f64,
    /// Density ratio ρ_A / ρ_B.
    pub density_ratio: f64,
}
impl MultiphaseParams {
    /// Create a new `MultiphaseParams`.
    pub fn new(surface_tension_coeff: f64, interface_width: f64, density_ratio: f64) -> Self {
        Self {
            surface_tension_coeff,
            interface_width,
            density_ratio,
        }
    }
}
/// Phase identifier.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Phase {
    /// Liquid phase (e.g., water)
    Liquid,
    /// Gas phase (e.g., air or vapor)
    Gas,
    /// Solid phase (fixed or moving boundary)
    Solid,
}
/// Method used to model immiscible phase interfaces.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum InterfaceModel {
    /// Transport a scalar colour function C and derive normals from ∇C.
    ColorFunction,
    /// Cahn–Hilliard phase-field order parameter (diffuse interface).
    PhaseField,
    /// Sharp-interface volume-of-fluid (geometric reconstruction).
    VolumeOfFluid,
}
/// Cahn-Hilliard phase field model for diffuse interface simulation.
///
/// ∂φ/∂t = ∇·(M∇μ)
/// μ = -ε²∇²φ + f'(φ)   with  f(φ) = ¼(φ²-1)²
#[derive(Debug, Clone)]
pub struct CahnHilliardField {
    /// Order parameter φ ∈ \[-1, +1\] (−1 = gas, +1 = liquid)
    pub phi: Vec<f64>,
    /// Chemical potential μ
    pub mu: Vec<f64>,
    /// Number of nodes
    pub n: usize,
    /// Grid spacing h \[m\]
    pub h: f64,
    /// Interface width parameter ε \[m\]
    pub epsilon: f64,
    /// Mobility M \[m²·s/kg\]
    pub mobility: f64,
    /// Cahn number (ε/L)
    pub cahn_number: f64,
}
impl CahnHilliardField {
    /// Create a 1D Cahn-Hilliard field.
    pub fn new_1d(n: usize, length: f64, epsilon: f64, mobility: f64) -> Self {
        let h = length / (n as f64 - 1.0);
        let phi = (0..n)
            .map(|i| {
                let x = i as f64 * h;
                let center = length * 0.5;
                (x - center).tanh() / epsilon.max(1e-10)
            })
            .collect();
        Self {
            phi,
            mu: vec![0.0; n],
            n,
            h,
            epsilon,
            mobility,
            cahn_number: epsilon / length,
        }
    }
    /// Evaluate the bulk free energy density: f(φ) = ¼(φ²−1)².
    #[inline]
    pub fn bulk_free_energy(&self, phi: f64) -> f64 {
        0.25 * (phi * phi - 1.0).powi(2)
    }
    /// Evaluate df/dφ = φ(φ²−1).
    #[inline]
    pub fn bulk_free_energy_deriv(&self, phi: f64) -> f64 {
        phi * (phi * phi - 1.0)
    }
    /// Laplacian of φ using second-order finite difference.
    fn laplacian_phi(&self, i: usize) -> f64 {
        if i == 0 || i == self.n - 1 {
            return 0.0;
        }
        (self.phi[i + 1] - 2.0 * self.phi[i] + self.phi[i - 1]) / (self.h * self.h)
    }
    /// Compute chemical potential: μ = −ε²∇²φ + f'(φ).
    pub fn compute_mu(&mut self) {
        for i in 0..self.n {
            let lap = self.laplacian_phi(i);
            self.mu[i] =
                -self.epsilon * self.epsilon * lap + self.bulk_free_energy_deriv(self.phi[i]);
        }
    }
    /// Explicit time step: ∂φ/∂t = M·∇²μ.
    pub fn step(&mut self, dt: f64) {
        self.compute_mu();
        let mut new_phi = self.phi.clone();
        let h2 = self.h * self.h;
        for (window, (phi_mid, new_phi_mid)) in self
            .mu
            .windows(3)
            .zip(self.phi[1..].iter().zip(new_phi[1..].iter_mut()))
        {
            let lap_mu = (window[2] - 2.0 * window[1] + window[0]) / h2;
            *new_phi_mid = (phi_mid + dt * self.mobility * lap_mu).clamp(-1.5, 1.5);
        }
        self.phi = new_phi;
    }
    /// Total free energy F = ∫ \[ε²/2 |∇φ|² + f(φ)\] dx.
    pub fn total_free_energy(&self) -> f64 {
        let mut f = 0.0;
        for i in 0..self.n {
            let grad_phi = if i == 0 {
                (self.phi[1] - self.phi[0]) / self.h
            } else if i == self.n - 1 {
                (self.phi[self.n - 1] - self.phi[self.n - 2]) / self.h
            } else {
                (self.phi[i + 1] - self.phi[i - 1]) / (2.0 * self.h)
            };
            f += (0.5 * self.epsilon * self.epsilon * grad_phi * grad_phi
                + self.bulk_free_energy(self.phi[i]))
                * self.h;
        }
        f
    }
    /// Volume fraction (fraction of nodes with φ > 0).
    pub fn volume_fraction(&self) -> f64 {
        let n_pos = self.phi.iter().filter(|&&p| p > 0.0).count();
        n_pos as f64 / self.n as f64
    }
    /// Interface location (zero crossing of φ).
    pub fn interface_position(&self) -> Option<f64> {
        for i in 0..self.n - 1 {
            if self.phi[i] == 0.0 {
                return Some(i as f64 * self.h);
            }
            if self.phi[i] * self.phi[i + 1] < 0.0 {
                let x_i = i as f64 * self.h;
                let frac = -self.phi[i] / (self.phi[i + 1] - self.phi[i]);
                return Some(x_i + frac * self.h);
            }
        }
        if self.phi[self.n - 1] == 0.0 {
            return Some((self.n - 1) as f64 * self.h);
        }
        None
    }
}
/// Physical properties of a single phase.
#[derive(Debug, Clone)]
pub struct PhaseProperties {
    /// Density \[kg/m³\]
    pub density: f64,
    /// Dynamic viscosity \[Pa·s\]
    pub viscosity: f64,
    /// Phase label
    pub phase: Phase,
}
impl PhaseProperties {
    /// Create liquid-like properties.
    pub fn water() -> Self {
        Self {
            density: 1000.0,
            viscosity: 1e-3,
            phase: Phase::Liquid,
        }
    }
    /// Create air-like properties.
    pub fn air() -> Self {
        Self {
            density: 1.225,
            viscosity: 1.8e-5,
            phase: Phase::Gas,
        }
    }
    /// Density ratio with respect to another phase.
    pub fn density_ratio(&self, other: &Self) -> f64 {
        self.density / other.density.max(1e-30)
    }
    /// Atwood number At = (ρ₁ − ρ₂) / (ρ₁ + ρ₂).
    pub fn atwood_number(heavy: &Self, light: &Self) -> f64 {
        (heavy.density - light.density) / (heavy.density + light.density).max(1e-30)
    }
}
/// Continuum Surface Force (CSF) model (Brackbill et al.).
///
/// f_surface = σ · κ · n · δ
/// where δ is a smoothed Dirac delta at the interface.
#[derive(Debug, Clone)]
pub struct CsfSurfaceTension {
    /// Surface tension coefficient σ \[N/m\]
    pub sigma: f64,
    /// Interface half-width for delta approximation \[m\]
    pub epsilon: f64,
}
impl CsfSurfaceTension {
    /// Create a CSF surface tension model.
    pub fn new(sigma: f64, epsilon: f64) -> Self {
        Self { sigma, epsilon }
    }
    /// Smoothed delta function approximation at the interface.
    ///
    /// δ(φ) = (1/2ε)(1 + cos(πφ/ε)) if |φ| < ε, else 0
    pub fn delta(&self, phi: f64) -> f64 {
        if phi.abs() < self.epsilon {
            (1.0 + (PI * phi / self.epsilon).cos()) / (2.0 * self.epsilon)
        } else {
            0.0
        }
    }
    /// CSF body force at a particle with given curvature, normal, and color.
    ///
    /// f_i = σ · κ · n_i · δ(φ)
    pub fn csf_force(&self, curvature: f64, normal: [f64; 2], phi: f64) -> [f64; 2] {
        let d = self.delta(phi);
        let f = self.sigma * curvature * d;
        [f * normal[0], f * normal[1]]
    }
    /// Continuum Surface Stress (CSS) model stress tensor.
    ///
    /// T_ij = σ · (|∇C|δᵢⱼ − n_i·n_j) · |∇C|
    pub fn css_stress_tensor(&self, normal: [f64; 2], grad_c_norm: f64) -> [[f64; 2]; 2] {
        let s = self.sigma * grad_c_norm;
        [
            [
                s * (1.0 - normal[0] * normal[0]),
                -s * normal[0] * normal[1],
            ],
            [
                -s * normal[1] * normal[0],
                s * (1.0 - normal[1] * normal[1]),
            ],
        ]
    }
    /// Laplace pressure jump across a spherical interface of radius R.
    ///
    /// Δp = 2σ/R (3D) or σ/R (2D cylindrical)
    pub fn laplace_pressure_3d(&self, radius: f64) -> f64 {
        2.0 * self.sigma / radius.max(1e-30)
    }
    /// Laplace pressure for a 2-D cylindrical interface: ΔP = σ/R.
    pub fn laplace_pressure_2d(&self, radius: f64) -> f64 {
        self.sigma / radius.max(1e-30)
    }
    /// Weber number We = ρ·V²·L / σ.
    pub fn weber_number(&self, rho: f64, velocity: f64, length: f64) -> f64 {
        rho * velocity * velocity * length / self.sigma.max(1e-30)
    }
    /// Capillary number Ca = μ·V / σ.
    pub fn capillary_number(&self, viscosity: f64, velocity: f64) -> f64 {
        viscosity * velocity / self.sigma.max(1e-30)
    }
}
/// Droplet representing a liquid region.
#[derive(Debug, Clone)]
pub struct Droplet {
    /// Center of mass position \[m\]
    pub center: [f64; 2],
    /// Velocity \[m/s\]
    pub velocity: [f64; 2],
    /// Equivalent radius \[m\]
    pub radius: f64,
    /// Volume \[m³\] (2D: area \[m²\])
    pub volume: f64,
    /// Density \[kg/m³\]
    pub density: f64,
    /// Surface tension coefficient \[N/m\]
    pub sigma: f64,
    /// Unique identifier
    pub id: usize,
}
impl Droplet {
    /// Create a spherical droplet.
    pub fn new(id: usize, x: f64, y: f64, radius: f64, density: f64, sigma: f64) -> Self {
        let volume = PI * radius * radius;
        Self {
            center: [x, y],
            velocity: [0.0; 2],
            radius,
            volume,
            density,
            sigma,
            id,
        }
    }
    /// Distance between two droplet centers.
    pub fn distance_to(&self, other: &Droplet) -> f64 {
        let dx = self.center[0] - other.center[0];
        let dy = self.center[1] - other.center[1];
        (dx * dx + dy * dy).sqrt()
    }
    /// Returns true if two droplets overlap (centers closer than sum of radii).
    pub fn overlaps(&self, other: &Droplet) -> bool {
        self.distance_to(other) < self.radius + other.radius
    }
    /// Coalesce two droplets into one (volume conservation).
    ///
    /// New radius: r = sqrt((r1² + r2²)) \[2D\]
    pub fn coalesce(d1: &Droplet, d2: &Droplet) -> Droplet {
        let v_total = d1.volume + d2.volume;
        let r_new = (v_total / PI).sqrt();
        let m1 = d1.density * d1.volume;
        let m2 = d2.density * d2.volume;
        let m_total = m1 + m2;
        let cx = (m1 * d1.center[0] + m2 * d2.center[0]) / m_total.max(1e-30);
        let cy = (m1 * d1.center[1] + m2 * d2.center[1]) / m_total.max(1e-30);
        let vx = (m1 * d1.velocity[0] + m2 * d2.velocity[0]) / m_total.max(1e-30);
        let vy = (m1 * d1.velocity[1] + m2 * d2.velocity[1]) / m_total.max(1e-30);
        let mut d_new = Droplet::new(d1.id, cx, cy, r_new, d1.density, d1.sigma);
        d_new.velocity = [vx, vy];
        d_new
    }
    /// Check if droplet should break up (Weber number > critical value).
    ///
    /// We_crit ≈ 12 for spherical droplets.
    pub fn should_breakup(&self, ambient_density: f64, relative_vel: f64, we_crit: f64) -> bool {
        let we =
            ambient_density * relative_vel * relative_vel * self.radius / self.sigma.max(1e-30);
        we > we_crit
    }
    /// Break up a droplet into two child droplets (symmetric split).
    pub fn breakup(&self, split_direction: [f64; 2]) -> (Droplet, Droplet) {
        let r_child = self.radius / 2.0_f64.sqrt();
        let offset = r_child * 0.5;
        let mut d1 = Droplet::new(
            self.id * 2,
            self.center[0] + offset * split_direction[0],
            self.center[1] + offset * split_direction[1],
            r_child,
            self.density,
            self.sigma,
        );
        let mut d2 = Droplet::new(
            self.id * 2 + 1,
            self.center[0] - offset * split_direction[0],
            self.center[1] - offset * split_direction[1],
            r_child,
            self.density,
            self.sigma,
        );
        d1.velocity = self.velocity;
        d2.velocity = self.velocity;
        (d1, d2)
    }
    /// Kinetic energy of this droplet.
    pub fn kinetic_energy(&self) -> f64 {
        let mass = self.density * self.volume;
        let v2 = self.velocity[0].powi(2) + self.velocity[1].powi(2);
        0.5 * mass * v2
    }
    /// Surface energy E_s = σ · 2π · r (2D perimeter).
    pub fn surface_energy(&self) -> f64 {
        self.sigma * 2.0 * PI * self.radius
    }
}
/// Liquid-gas-solid triple junction (contact line) dynamics.
#[derive(Debug, Clone)]
pub struct TripleJunction {
    /// Contact line position \[m\]
    pub position: f64,
    /// Contact line velocity \[m/s\]
    pub velocity: f64,
    /// Contact angle model
    pub contact_model: ContactAngleModel,
    /// Current contact angle \[rad\]
    pub current_angle: f64,
    /// Friction coefficient for contact line motion
    pub friction: f64,
}
impl TripleJunction {
    /// Create a new triple junction at rest.
    pub fn new(x: f64, contact_model: ContactAngleModel) -> Self {
        let theta0 = contact_model.theta_y;
        Self {
            position: x,
            velocity: 0.0,
            contact_model,
            current_angle: theta0,
            friction: 1e-3,
        }
    }
    /// Contact line driving force: F = σ_lg · (cos θ_Y − cos θ_c).
    pub fn driving_force(&self) -> f64 {
        self.contact_model.sigma_lg * (self.contact_model.theta_y.cos() - self.current_angle.cos())
    }
    /// Update contact line position using overdamped dynamics.
    pub fn step(&mut self, dt: f64) {
        let f = self.driving_force();
        self.velocity = f / self.friction.max(1e-30);
        self.position += self.velocity * dt;
        let ca = self.contact_model.capillary_number_from_velocity();
        let ln_r = 10.0;
        self.current_angle = self.contact_model.dynamic_contact_angle(ca, ln_r);
    }
}
/// Physical properties of a single fluid phase (generalised form).
///
/// Unlike the discriminated-union [`Phase`] enum used earlier in this file,
/// this struct allows an arbitrary number of phases to be registered at
/// runtime.
#[derive(Debug, Clone)]
pub struct PhaseDesc {
    /// Unique integer identifier for this phase.
    pub id: usize,
    /// Reference density ρ₀ \[kg/m³\].
    pub density: f64,
    /// Dynamic viscosity μ \[Pa·s\].
    pub viscosity: f64,
    /// Surface-tension coefficient σ \[N/m\].
    pub surface_tension_coeff: f64,
}
impl PhaseDesc {
    /// Create a new phase descriptor.
    pub fn new(id: usize, density: f64, viscosity: f64, surface_tension_coeff: f64) -> Self {
        Self {
            id,
            density,
            viscosity,
            surface_tension_coeff,
        }
    }
    /// Capillary length √(σ / (ρ g)) \[m\].
    pub fn capillary_length(&self, g: f64) -> f64 {
        (self.surface_tension_coeff / (self.density * g)).sqrt()
    }
}
/// SPH particle carrying phase-field information.
///
/// Uses 3-D `[f64; 3]` arrays for position and velocity so that the struct
/// is usable in both 2-D (z=0) and 3-D simulations without nalgebra.
#[derive(Debug, Clone)]
pub struct MultiphaseParticle {
    /// Position \[m\].
    pub position: [f64; 3],
    /// Velocity \[m/s\].
    pub velocity: [f64; 3],
    /// Index into the phase table.
    pub phase_id: usize,
    /// Colour function value C ∈ \[0, 1\] (1 = phase 0, 0 = phase 1).
    pub color_function: f64,
    /// Interface curvature κ \[1/m\] (computed via [`compute_curvature`]).
    pub curvature: f64,
    /// Mass \[kg\].
    pub mass: f64,
}
impl MultiphaseParticle {
    /// Construct a particle at rest with the given phase id.
    pub fn new(position: [f64; 3], phase_id: usize, mass: f64) -> Self {
        let color_function = if phase_id == 0 { 1.0 } else { 0.0 };
        Self {
            position,
            velocity: [0.0; 3],
            phase_id,
            color_function,
            curvature: 0.0,
            mass,
        }
    }
    /// Squared distance to another particle.
    pub fn dist2(&self, other: &MultiphaseParticle) -> f64 {
        let dx = self.position[0] - other.position[0];
        let dy = self.position[1] - other.position[1];
        let dz = self.position[2] - other.position[2];
        dx * dx + dy * dy + dz * dz
    }
    /// Distance to another particle.
    pub fn dist(&self, other: &MultiphaseParticle) -> f64 {
        self.dist2(other).sqrt()
    }
    /// Kinetic energy of this particle.
    pub fn kinetic_energy(&self) -> f64 {
        let v2: f64 = self.velocity.iter().map(|v| v * v).sum();
        0.5 * self.mass * v2
    }
}
/// Contact angle model for liquid-solid-gas triple junction.
///
/// Young's equation: σ_sg − σ_sl = σ_lg · cos(θ_Y)
#[derive(Debug, Clone)]
pub struct ContactAngleModel {
    /// Young's contact angle θ_Y \[rad\]
    pub theta_y: f64,
    /// Dynamic advancing contact angle \[rad\]
    pub theta_adv: f64,
    /// Dynamic receding contact angle \[rad\]
    pub theta_rec: f64,
    /// Liquid-gas surface tension σ_lg \[N/m\]
    pub sigma_lg: f64,
    /// Solid-gas surface energy σ_sg \[N/m\]
    pub sigma_sg: f64,
    /// Solid-liquid surface energy σ_sl \[N/m\]
    pub sigma_sl: f64,
    /// Contact line velocity threshold for dynamic effects \[m/s\]
    pub v_th: f64,
}
impl ContactAngleModel {
    /// Create a contact angle model from interface energies.
    pub fn new(sigma_lg: f64, sigma_sg: f64, sigma_sl: f64) -> Self {
        let cos_theta = (sigma_sg - sigma_sl) / sigma_lg.max(1e-30);
        let theta_y = cos_theta.clamp(-1.0, 1.0).acos();
        let theta_adv = theta_y * 1.1;
        let theta_rec = theta_y * 0.9;
        Self {
            theta_y,
            theta_adv,
            theta_rec,
            sigma_lg,
            sigma_sg,
            sigma_sl,
            v_th: 1e-3,
        }
    }
    /// Young-Dupré work of adhesion W_a = σ_lg · (1 + cos θ_Y).
    pub fn work_of_adhesion(&self) -> f64 {
        self.sigma_lg * (1.0 + self.theta_y.cos())
    }
    /// Spreading coefficient S = σ_sg − σ_sl − σ_lg.
    pub fn spreading_coefficient(&self) -> f64 {
        self.sigma_sg - self.sigma_sl - self.sigma_lg
    }
    /// Dynamic contact angle using Cox-Voinov model.
    ///
    /// θ³(Ca) ≈ θ_Y³ + 9·Ca·ln(L/λ)
    pub fn dynamic_contact_angle(&self, ca: f64, ln_ratio: f64) -> f64 {
        let theta3 = self.theta_y.powi(3) + 9.0 * ca * ln_ratio;
        theta3.max(0.0).cbrt()
    }
    /// Capillary length l_c = sqrt(σ / (ρ · g)).
    pub fn capillary_length(&self, rho: f64, g: f64) -> f64 {
        (self.sigma_lg / (rho * g).max(1e-30)).sqrt()
    }
    /// Bond number Bo = ρ·g·R² / σ.
    pub fn bond_number(&self, rho: f64, g: f64, radius: f64) -> f64 {
        rho * g * radius * radius / self.sigma_lg.max(1e-30)
    }
    /// Returns true if the surface is hydrophilic (θ_Y < π/2).
    pub fn is_hydrophilic(&self) -> bool {
        self.theta_y < PI * 0.5
    }
    /// Returns true if the surface is superhydrophobic (θ_Y > 150°).
    pub fn is_superhydrophobic(&self) -> bool {
        self.theta_y > 150.0_f64.to_radians()
    }
}
impl ContactAngleModel {
    /// Capillary number from stored threshold.
    fn capillary_number_from_velocity(&self) -> f64 {
        self.v_th * 1e-3 / self.sigma_lg.max(1e-30)
    }
}
/// Boussinesq approximation for buoyancy-driven (density-stratified) flow.
#[derive(Debug, Clone)]
pub struct BoussinesqModel {
    /// Reference density ρ₀ \[kg/m³\]
    pub rho_ref: f64,
    /// Gravitational acceleration \[m/s²\] (positive downward)
    pub g: f64,
    /// Thermal expansion coefficient β_T \[1/K\]
    pub beta_t: f64,
    /// Solutal expansion coefficient β_c \[m³/mol\]
    pub beta_c: f64,
    /// Reference temperature T₀ \[K\]
    pub t_ref: f64,
    /// Reference concentration c₀ \[mol/m³\]
    pub c_ref: f64,
}
impl BoussinesqModel {
    /// Create a standard Boussinesq model for water.
    pub fn water(g: f64) -> Self {
        Self {
            rho_ref: 1000.0,
            g,
            beta_t: 2.07e-4,
            beta_c: 0.0,
            t_ref: 293.15,
            c_ref: 0.0,
        }
    }
    /// Buoyancy force per unit volume: f_b = −ρ₀ · g · (β_T·ΔT + β_c·Δc)
    pub fn buoyancy_force(&self, temp: f64, conc: f64) -> f64 {
        -self.rho_ref
            * self.g
            * (self.beta_t * (temp - self.t_ref) + self.beta_c * (conc - self.c_ref))
    }
    /// Effective density: ρ ≈ ρ₀ · (1 − β_T·ΔT).
    pub fn effective_density(&self, temp: f64) -> f64 {
        self.rho_ref * (1.0 - self.beta_t * (temp - self.t_ref))
    }
    /// Rayleigh number Ra = g·β_T·ΔT·L³ / (ν·α).
    pub fn rayleigh_number(&self, delta_t: f64, length: f64, nu: f64, alpha_th: f64) -> f64 {
        self.g * self.beta_t * delta_t.abs() * length.powi(3) / (nu * alpha_th).max(1e-30)
    }
    /// Prandtl number Pr = ν / α.
    pub fn prandtl_number(&self, nu: f64, alpha_th: f64) -> f64 {
        nu / alpha_th.max(1e-30)
    }
}
/// Surface-tension force model.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SurfaceTensionForce {
    /// Continuum Surface Force (Brackbill, Kothe & Zemach, 1992).
    Csf,
    /// Pairwise inter-particle surface force (Morris, 2000).
    Pairwise,
}
/// Multi-fluid SPH density summation with smoothed interface.
///
/// Uses a modified kernel summation that accounts for the density ratio at
/// the interface (Morris et al. 2000).
#[derive(Debug, Clone)]
pub struct MultiphaseDensity {
    /// Smoothing length h \[m\].
    pub h: f64,
    /// Phase-0 reference density ρ₀ \[kg/m³\].
    pub rho0: f64,
    /// Phase-1 reference density ρ₁ \[kg/m³\].
    pub rho1: f64,
}
impl MultiphaseDensity {
    /// Create a new `MultiphaseDensity`.
    pub fn new(h: f64, rho0: f64, rho1: f64) -> Self {
        Self { h, rho0, rho1 }
    }
    /// Density ratio Λ = ρ₀ / ρ₁.
    pub fn density_ratio(&self) -> f64 {
        self.rho0 / self.rho1.max(1e-30)
    }
    /// Compute SPH density for particle `i` using kernel summation.
    ///
    /// ρ_i = Σ_j m_j W(|r_ij|, h)
    pub fn compute_density(&self, particles: &[MultiphaseSphParticle], i: usize) -> f64 {
        let xi = particles[i].position;
        let mut rho = 0.0_f64;
        for p in particles {
            let dx = xi[0] - p.position[0];
            let dy = xi[1] - p.position[1];
            let dz = xi[2] - p.position[2];
            let r = (dx * dx + dy * dy + dz * dz).sqrt();
            if r < self.h {
                rho += p.mass * cubic_spline_kernel(r, self.h);
            }
        }
        rho
    }
    /// Interpolated density from colour function: ρ = C ρ₀ + (1-C) ρ₁.
    pub fn interpolated_density(&self, color: f64) -> f64 {
        color * self.rho0 + (1.0 - color) * self.rho1
    }
    /// Update densities for all particles using the summation formula.
    pub fn update_all(&self, particles: &mut [MultiphaseSphParticle]) {
        let n = particles.len();
        let densities: Vec<f64> = (0..n).map(|i| self.compute_density(particles, i)).collect();
        for (p, &rho) in particles.iter_mut().zip(densities.iter()) {
            if rho > 1e-30 {
                p.density = rho;
            }
        }
    }
}
/// Tait equation of state for each phase with density ratio.
///
/// p = K_α · ((ρ/ρ₀_α)^γ − 1)
///
/// where the stiffness K is chosen to keep Mach number < 0.1.
#[derive(Debug, Clone)]
pub struct MultiphaseEquationOfState {
    /// Reference density for phase 0 \[kg/m³\].
    pub rho0_phase0: f64,
    /// Reference density for phase 1 \[kg/m³\].
    pub rho0_phase1: f64,
    /// Tait stiffness K for phase 0 \[Pa\].
    pub k_phase0: f64,
    /// Tait stiffness K for phase 1 \[Pa\].
    pub k_phase1: f64,
    /// Polytropic exponent γ (typically 7 for water).
    pub gamma: f64,
}
impl MultiphaseEquationOfState {
    /// Create a new `MultiphaseEquationOfState`.
    pub fn new(
        rho0_phase0: f64,
        rho0_phase1: f64,
        k_phase0: f64,
        k_phase1: f64,
        gamma: f64,
    ) -> Self {
        Self {
            rho0_phase0,
            rho0_phase1,
            k_phase0,
            k_phase1,
            gamma,
        }
    }
    /// Create a water-air EOS with γ = 7.
    pub fn water_air() -> Self {
        Self::new(1000.0, 1.225, 1e5, 1e3, 7.0)
    }
    /// Pressure for a particle of the given phase and density.
    pub fn pressure(&self, density: f64, phase_label: usize) -> f64 {
        let (rho0, k) = if phase_label == 0 {
            (self.rho0_phase0, self.k_phase0)
        } else {
            (self.rho0_phase1, self.k_phase1)
        };
        k * ((density / rho0.max(1e-30)).powf(self.gamma) - 1.0)
    }
    /// Speed of sound c = sqrt(γ K / ρ₀) for the given phase.
    pub fn sound_speed(&self, phase_label: usize) -> f64 {
        let (rho0, k) = if phase_label == 0 {
            (self.rho0_phase0, self.k_phase0)
        } else {
            (self.rho0_phase1, self.k_phase1)
        };
        (self.gamma * k / rho0.max(1e-30)).sqrt()
    }
    /// Density ratio ρ₀_0 / ρ₀_1.
    pub fn density_ratio(&self) -> f64 {
        self.rho0_phase0 / self.rho0_phase1.max(1e-30)
    }
    /// Update pressures of all particles.
    pub fn update_pressures(&self, particles: &mut [MultiphaseSphParticle]) {
        for p in particles.iter_mut() {
            p.pressure = self.pressure(p.density, p.phase_label);
        }
    }
    /// Maximum time step for stability: Δt_max = C_CFL · h / c_max.
    pub fn cfl_timestep(&self, h: f64, cfl: f64) -> f64 {
        let c0 = self.sound_speed(0);
        let c1 = self.sound_speed(1);
        let c_max = c0.max(c1);
        if c_max <= 0.0 {
            f64::INFINITY
        } else {
            cfl * h / c_max
        }
    }
}
/// Rayleigh-Plesset bubble dynamics integrated within an SPH context.
///
/// Tracks a spherical bubble radius R(t) in an incompressible liquid under
/// the pressure field provided by the surrounding SPH particles.
#[derive(Debug, Clone)]
pub struct BubbleDynamicsModel {
    /// Current bubble radius R \[m\].
    pub radius: f64,
    /// Bubble-wall velocity Ṙ \[m/s\].
    pub radius_dot: f64,
    /// Initial radius R₀ \[m\].
    pub r0: f64,
    /// Liquid density ρ_l \[kg/m³\].
    pub rho_l: f64,
    /// Liquid dynamic viscosity μ_l \[Pa·s\].
    pub mu_l: f64,
    /// Surface tension σ \[N/m\].
    pub sigma: f64,
    /// Reference gas pressure p_b0 at R = R₀ \[Pa\].
    pub p_b0: f64,
    /// Far-field liquid pressure p_∞ \[Pa\].
    pub p_inf: f64,
    /// Polytropic index γ.
    pub gamma: f64,
    /// Bubble position in the SPH domain \[m\].
    pub position: [f64; 3],
}
impl BubbleDynamicsModel {
    /// Create a new `BubbleDynamicsModel` at equilibrium (p_b0 = p_∞ + 2σ/R₀).
    pub fn new(r0: f64, p_inf: f64, rho_l: f64, mu_l: f64, sigma: f64, position: [f64; 3]) -> Self {
        let p_b0 = p_inf + 2.0 * sigma / r0.max(1e-30);
        Self {
            radius: r0,
            radius_dot: 0.0,
            r0,
            rho_l,
            mu_l,
            sigma,
            p_b0,
            p_inf,
            gamma: 1.4,
            position,
        }
    }
    /// Gas pressure at radius R: p_b = p_b0 (R₀/R)^(3γ).
    pub fn gas_pressure(&self) -> f64 {
        self.p_b0 * (self.r0 / self.radius.max(1e-30)).powf(3.0 * self.gamma)
    }
    /// Rayleigh-Plesset acceleration R̈.
    pub fn radius_ddot(&self, p_ambient: f64) -> f64 {
        let r = self.radius.max(1e-30);
        let p_b = self.gas_pressure();
        let dp = (p_b - p_ambient) / self.rho_l.max(1e-30);
        let visc_term = 4.0 * self.mu_l * self.radius_dot / (self.rho_l.max(1e-30) * r);
        let surf_term = 2.0 * self.sigma / (self.rho_l.max(1e-30) * r);
        let inertial = 1.5 * self.radius_dot * self.radius_dot / r;
        (dp - visc_term - surf_term) / r - inertial
    }
    /// Advance one time step using Velocity Verlet integration.
    pub fn step(&mut self, dt: f64, p_ambient: f64) {
        let r_ddot = self.radius_ddot(p_ambient);
        self.radius += self.radius_dot * dt + 0.5 * r_ddot * dt * dt;
        self.radius = self.radius.max(1e-12);
        let r_ddot_new = self.radius_ddot(p_ambient);
        self.radius_dot += 0.5 * (r_ddot + r_ddot_new) * dt;
    }
    /// Minnaert natural frequency ω_M = sqrt(3γ p_∞ / (ρ_l R₀²)).
    pub fn minnaert_frequency(&self) -> f64 {
        (3.0 * self.gamma * self.p_inf / (self.rho_l * self.r0 * self.r0)).sqrt()
    }
    /// Rayleigh collapse time t_c ≈ 0.915 R₀ sqrt(ρ_l / p_∞).
    pub fn collapse_time(&self) -> f64 {
        0.915 * self.r0 * (self.rho_l / self.p_inf.max(1e-30)).sqrt()
    }
    /// SPH pressure at a point `x` due to the bubble, modeled as a monopole:
    ///
    /// Δp(r) ≈ ρ_l R² Ṙ / r  (far-field pressure radiation)
    pub fn radiated_pressure_at(&self, x: [f64; 3]) -> f64 {
        let dx = x[0] - self.position[0];
        let dy = x[1] - self.position[1];
        let dz = x[2] - self.position[2];
        let dist = (dx * dx + dy * dy + dz * dz).sqrt().max(1e-30);
        self.rho_l * self.radius * self.radius * self.radius_dot / dist
    }
    /// Growth-phase check: returns true if bubble is expanding (Ṙ > 0).
    pub fn is_growing(&self) -> bool {
        self.radius_dot > 0.0
    }
}
