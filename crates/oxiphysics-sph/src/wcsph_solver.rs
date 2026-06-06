// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Advanced Weakly Compressible SPH (WCSPH) solver.
//!
//! Implements a complete WCSPH pipeline including:
//!
//! - [`WcsphConfig`]: Solver configuration (EOS, kernel, compressibility)
//! - [`TaitEos`]: Tait equation of state for quasi-incompressible flow
//! - [`MorrisEos`]: Linear EOS for very weakly compressible flow
//! - [`WcsphSolverAdvanced`]: Main solver with density summation and momentum integration
//! - [`ArtificialViscosity`]: Monaghan 1992 α-β artificial viscosity
//! - [`XsphCorrection`]: XSPH position correction for particle regularization
//! - [`DeltaSph`]: Delta-SPH density diffusion (Molteni-Colagrossi)
//! - [`RiemannSph`]: SPH with Riemann solver for inter-particle fluxes
//! - [`WcsphTimeIntegrator`]: Verlet/leapfrog/predictor-corrector integration
//! - [`WcsphBoundary`]: Dynamic boundary particles, Lennard-Jones, Adami BC
//! - [`WcsphMultiPhase`]: Multi-phase WCSPH with surface tension

// ============================================================================
// Configuration
// ============================================================================

/// Kernel type selection for WCSPH.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum KernelType {
    /// Cubic spline kernel (default).
    CubicSpline,
    /// Wendland C2 kernel.
    WendlandC2,
    /// Quintic Wendland kernel.
    WendlandC4,
    /// Gaussian kernel.
    Gaussian,
}

/// Configuration for the advanced WCSPH solver.
///
/// Controls physical parameters, numerical settings, and kernel selection.
#[derive(Debug, Clone)]
pub struct WcsphConfig {
    /// Reference (rest) density ρ₀ (kg/m³).
    pub rho0: f64,
    /// Speed of sound c₀ (m/s). Typically 10× max expected velocity.
    pub c0: f64,
    /// Compressibility factor (dimensionless, typically 1-10).
    pub compressibility_factor: f64,
    /// Gamma exponent for Tait EOS (water: 7).
    pub gamma: f64,
    /// Smoothing length h (m).
    pub h: f64,
    /// Particle mass (kg).
    pub particle_mass: f64,
    /// Kernel type selection.
    pub kernel: KernelType,
}

impl WcsphConfig {
    /// Create a new WCSPH configuration with water-like parameters.
    pub fn water(rho0: f64, c0: f64, h: f64, particle_mass: f64) -> Self {
        Self {
            rho0,
            c0,
            compressibility_factor: 1.0,
            gamma: 7.0,
            h,
            particle_mass,
            kernel: KernelType::CubicSpline,
        }
    }

    /// Tait EOS stiffness parameter B = ρ₀ * c₀² / γ.
    pub fn tait_b(&self) -> f64 {
        self.rho0 * self.c0 * self.c0 / self.gamma
    }

    /// CFL-limited time step: dt = CFL * h / c0.
    pub fn cfl_dt(&self, cfl: f64) -> f64 {
        cfl * self.h / self.c0
    }
}

impl Default for WcsphConfig {
    fn default() -> Self {
        Self::water(1000.0, 14.14, 0.01, 1e-6)
    }
}

// ============================================================================
// Tait EOS
// ============================================================================

/// Tait equation of state for quasi-incompressible SPH.
///
/// p = B * ((ρ/ρ₀)^γ - 1)  where B = ρ₀ * c₀² / γ
///
/// The weak compressibility keeps density fluctuations below ~1% when
/// c₀ ≥ 10 * V_max.
#[derive(Debug, Clone)]
pub struct TaitEos {
    /// Reference density ρ₀ (kg/m³).
    pub rho0: f64,
    /// Speed of sound c₀ (m/s).
    pub c0: f64,
    /// Gamma exponent (7 for water).
    pub gamma: f64,
}

impl TaitEos {
    /// Create a new Tait EOS.
    pub fn new(rho0: f64, c0: f64, gamma: f64) -> Self {
        Self { rho0, c0, gamma }
    }

    /// Stiffness coefficient B = ρ₀ * c₀² / γ.
    pub fn b(&self) -> f64 {
        self.rho0 * self.c0 * self.c0 / self.gamma
    }

    /// Pressure from density: p = B * ((ρ/ρ₀)^γ - 1).
    pub fn pressure(&self, rho: f64) -> f64 {
        tait_pressure(rho, self.rho0, self.c0, self.gamma)
    }

    /// Density from pressure: ρ = ρ₀ * (p/B + 1)^(1/γ).
    pub fn density(&self, pressure: f64) -> f64 {
        let b = self.b();
        self.rho0 * (pressure / b + 1.0).powf(1.0 / self.gamma)
    }

    /// Sound speed at density ρ: c(ρ) = c₀ * (ρ/ρ₀)^((γ-1)/2).
    pub fn sound_speed(&self, rho: f64) -> f64 {
        tait_sound_speed(rho, self.rho0, self.c0, self.gamma)
    }

    /// Bulk modulus K = ρ * dp/dρ = γ * B * (ρ/ρ₀)^γ.
    pub fn bulk_modulus(&self, rho: f64) -> f64 {
        self.gamma * self.b() * (rho / self.rho0).powf(self.gamma)
    }
}

// ============================================================================
// Morris EOS
// ============================================================================

/// Linear (Morris) equation of state for quasi-incompressible flow.
///
/// p = c₀² * (ρ - ρ₀)
///
/// Simpler than Tait; appropriate for very weakly compressible regimes.
#[derive(Debug, Clone)]
pub struct MorrisEos {
    /// Reference density ρ₀ (kg/m³).
    pub rho0: f64,
    /// Speed of sound c₀ (m/s).
    pub c0: f64,
}

impl MorrisEos {
    /// Create a new Morris EOS.
    pub fn new(rho0: f64, c0: f64) -> Self {
        Self { rho0, c0 }
    }

    /// Pressure: p = c₀² * (ρ - ρ₀).
    pub fn pressure(&self, rho: f64) -> f64 {
        self.c0 * self.c0 * (rho - self.rho0)
    }

    /// Density: ρ = ρ₀ + p / c₀².
    pub fn density(&self, pressure: f64) -> f64 {
        self.rho0 + pressure / (self.c0 * self.c0)
    }

    /// Sound speed (constant c₀).
    pub fn sound_speed(&self) -> f64 {
        self.c0
    }
}

// ============================================================================
// WCSPH Solver (Advanced)
// ============================================================================

/// Advanced WCSPH solver combining density summation, pressure force, and viscosity.
///
/// Main simulation loop steps:
/// 1. Density summation or continuity equation
/// 2. Pressure from EOS
/// 3. Momentum equation (pressure force + viscosity + body force)
/// 4. Time integration
pub struct WcsphSolverAdvanced {
    /// Solver configuration.
    pub config: WcsphConfig,
    /// Particle positions.
    pub positions: Vec<[f64; 3]>,
    /// Particle velocities.
    pub velocities: Vec<[f64; 3]>,
    /// Particle densities.
    pub densities: Vec<f64>,
    /// Particle pressures.
    pub pressures: Vec<f64>,
    /// Current time (s).
    pub time: f64,
}

impl WcsphSolverAdvanced {
    /// Create a new WCSPH solver.
    pub fn new(config: WcsphConfig) -> Self {
        Self {
            config,
            positions: Vec::new(),
            velocities: Vec::new(),
            densities: Vec::new(),
            pressures: Vec::new(),
            time: 0.0,
        }
    }

    /// Add a particle at position with given velocity.
    pub fn add_particle(&mut self, pos: [f64; 3], vel: [f64; 3]) {
        self.positions.push(pos);
        self.velocities.push(vel);
        self.densities.push(self.config.rho0);
        self.pressures.push(0.0);
    }

    /// Number of particles.
    pub fn num_particles(&self) -> usize {
        self.positions.len()
    }

    /// Update pressures from densities using Tait EOS.
    pub fn update_pressures(&mut self) {
        let eos = TaitEos::new(self.config.rho0, self.config.c0, self.config.gamma);
        for (p, &rho) in self.pressures.iter_mut().zip(self.densities.iter()) {
            *p = eos.pressure(rho);
        }
    }

    /// Density summation at particle i using kernel W(r, h).
    ///
    /// ρᵢ = Σⱼ mⱼ * W(|xᵢ - xⱼ|, h)
    pub fn density_summation_i(&self, i: usize, kernel_fn: &dyn Fn(f64, f64) -> f64) -> f64 {
        let xi = self.positions[i];
        let h = self.config.h;
        let m = self.config.particle_mass;
        let mut rho = 0.0;
        for xj in &self.positions {
            let dx = xi[0] - xj[0];
            let dy = xi[1] - xj[1];
            let dz = xi[2] - xj[2];
            let r = (dx * dx + dy * dy + dz * dz).sqrt();
            rho += m * kernel_fn(r, h);
        }
        rho
    }

    /// Pressure force on particle i from neighbor j.
    ///
    /// f_p = -m * (pᵢ/ρᵢ² + pⱼ/ρⱼ²) * ∇W
    pub fn pressure_force_ij(
        &self,
        pi: f64,
        rhoi: f64,
        pj: f64,
        rhoj: f64,
        grad_w: [f64; 3],
        mj: f64,
    ) -> [f64; 3] {
        let coeff = -mj * (pi / (rhoi * rhoi) + pj / (rhoj * rhoj));
        [coeff * grad_w[0], coeff * grad_w[1], coeff * grad_w[2]]
    }

    /// Maximum particle speed.
    pub fn max_speed(&self) -> f64 {
        self.velocities
            .iter()
            .map(|v| (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt())
            .fold(0.0_f64, f64::max)
    }

    /// CFL time step.
    pub fn cfl_dt(&self, cfl: f64) -> f64 {
        cfl_dt_wcsph(self.max_speed(), self.config.c0, self.config.h, cfl)
    }
}

// ============================================================================
// Artificial Viscosity
// ============================================================================

/// Monaghan 1992 artificial viscosity for SPH.
///
/// Π_ij = { (-α * c̄_ij * μ_ij + β * μ_ij²) / ρ̄_ij  if v_ij · r_ij < 0
///         { 0                                          otherwise
///
/// where μ_ij = h * v_ij · r_ij / (|r_ij|² + ε*h²).
#[derive(Debug, Clone)]
pub struct ArtificialViscosity {
    /// Linear viscosity parameter α (typically 0.01–0.1).
    pub alpha: f64,
    /// Quadratic viscosity parameter β (typically 0.0–2α).
    pub beta: f64,
    /// Small stabilization parameter ε (typically 0.01).
    pub epsilon: f64,
}

impl ArtificialViscosity {
    /// Create artificial viscosity with given parameters.
    pub fn new(alpha: f64, beta: f64, epsilon: f64) -> Self {
        Self {
            alpha,
            beta,
            epsilon,
        }
    }

    /// Default Monaghan 1992 parameters: α=0.01, β=0.0, ε=0.01.
    pub fn default_params() -> Self {
        Self::new(0.01, 0.0, 0.01)
    }

    /// Compute artificial viscosity term Π_ij.
    pub fn compute(
        &self,
        v_ij: [f64; 3],
        r_ij: [f64; 3],
        r_ij_norm: f64,
        c_i: f64,
        c_j: f64,
        rho_i: f64,
        rho_j: f64,
        h: f64,
    ) -> f64 {
        let vr = v_ij[0] * r_ij[0] + v_ij[1] * r_ij[1] + v_ij[2] * r_ij[2];
        if vr >= 0.0 {
            return 0.0;
        }
        let mu = h * vr / (r_ij_norm * r_ij_norm + self.epsilon * h * h);
        let c_bar = 0.5 * (c_i + c_j);
        let rho_bar = 0.5 * (rho_i + rho_j);
        (-self.alpha * c_bar * mu + self.beta * mu * mu) / rho_bar
    }

    /// Viscous force on particle i from particle j.
    pub fn viscous_force(&self, pi_ij: f64, mj: f64, grad_w: [f64; 3]) -> [f64; 3] {
        let c = -mj * pi_ij;
        [c * grad_w[0], c * grad_w[1], c * grad_w[2]]
    }
}

// ============================================================================
// XSPH Correction
// ============================================================================

/// XSPH position correction for improved particle regularity.
///
/// Modifies the particle velocity used for position update:
/// dx_i/dt = v_i + ε_xsph * Σⱼ (mⱼ/ρ̄_ij) * (v_j - v_i) * W_ij
#[derive(Debug, Clone)]
pub struct XsphCorrection {
    /// XSPH parameter ε (typically 0.1–0.5).
    pub epsilon: f64,
}

impl XsphCorrection {
    /// Create XSPH correction with given parameter.
    pub fn new(epsilon: f64) -> Self {
        Self { epsilon }
    }

    /// Compute XSPH correction velocity for particle i from neighbor j.
    pub fn correction_ij(
        &self,
        v_i: [f64; 3],
        v_j: [f64; 3],
        rho_i: f64,
        rho_j: f64,
        mj: f64,
        w_ij: f64,
    ) -> [f64; 3] {
        let rho_bar = 0.5 * (rho_i + rho_j);
        let c = self.epsilon * mj / rho_bar * w_ij;
        [
            c * (v_j[0] - v_i[0]),
            c * (v_j[1] - v_i[1]),
            c * (v_j[2] - v_i[2]),
        ]
    }

    /// Corrected velocity for position update: v_corr = v_i + Σⱼ correction_ij.
    pub fn corrected_velocity(&self, v_i: [f64; 3], sum_correction: [f64; 3]) -> [f64; 3] {
        [
            v_i[0] + sum_correction[0],
            v_i[1] + sum_correction[1],
            v_i[2] + sum_correction[2],
        ]
    }
}

// ============================================================================
// Delta-SPH
// ============================================================================

/// Delta-SPH density diffusion term (Molteni-Colagrossi 2009).
///
/// Adds a density diffusion term to the continuity equation:
/// Dρᵢ/Dt = ... + δ * h * c₀ * Σⱼ mⱼ/ρⱼ * ψ_ij · ∇W_ij
///
/// where ψ_ij = 2(ρⱼ - ρᵢ) * r_ij / |r_ij|²
#[derive(Debug, Clone)]
pub struct DeltaSph {
    /// Delta-SPH coefficient δ (typically 0.1).
    pub delta: f64,
    /// Speed of sound c₀ (m/s).
    pub c0: f64,
    /// Smoothing length h (m).
    pub h: f64,
}

impl DeltaSph {
    /// Create a Delta-SPH diffusion term.
    pub fn new(delta: f64, c0: f64, h: f64) -> Self {
        Self { delta, c0, h }
    }

    /// Density diffusion kernel ψ_ij.
    ///
    /// ψ_ij = 2 * (ρⱼ - ρᵢ) * r_ij / |r_ij|²
    pub fn psi_ij(&self, rho_i: f64, rho_j: f64, r_ij: [f64; 3], r_ij_sq: f64) -> [f64; 3] {
        let factor = 2.0 * (rho_j - rho_i) / (r_ij_sq + 1e-30);
        [factor * r_ij[0], factor * r_ij[1], factor * r_ij[2]]
    }

    /// Density diffusion contribution from neighbor j.
    ///
    /// D_ij = δ * h * c₀ * mⱼ/ρⱼ * (ψ_ij · ∇W_ij)
    pub fn diffusion_ij(
        &self,
        rho_i: f64,
        rho_j: f64,
        mj: f64,
        r_ij: [f64; 3],
        r_ij_sq: f64,
        grad_w: [f64; 3],
    ) -> f64 {
        let psi = self.psi_ij(rho_i, rho_j, r_ij, r_ij_sq);
        let psi_dot_grad = psi[0] * grad_w[0] + psi[1] * grad_w[1] + psi[2] * grad_w[2];
        self.delta * self.h * self.c0 * mj / rho_j * psi_dot_grad
    }
}

// ============================================================================
// Riemann SPH
// ============================================================================

/// SPH with approximate Riemann solver for inter-particle fluxes.
///
/// Replaces the standard SPH pressure gradient with a Riemann-based
/// flux that provides low numerical dissipation and good shock capturing.
#[derive(Debug, Clone)]
pub struct RiemannSph {
    /// Coefficient for Riemann dissipation (0.5–1.0).
    pub beta_r: f64,
    /// Reconstruction limiter type: 0=none, 1=van Leer, 2=minmod.
    pub limiter: u8,
}

impl RiemannSph {
    /// Create a Riemann SPH solver.
    pub fn new(beta_r: f64, limiter: u8) -> Self {
        Self { beta_r, limiter }
    }

    /// Acoustic Riemann solver: compute inter-particle pressure and velocity.
    ///
    /// Returns (p_star, v_star_n) where v_star_n is the normal velocity at
    /// the Riemann interface.
    pub fn acoustic_riemann(
        &self,
        p_l: f64,
        p_r: f64,
        v_l_n: f64,
        v_r_n: f64,
        rho_l: f64,
        rho_r: f64,
        c_l: f64,
        c_r: f64,
    ) -> (f64, f64) {
        let z_l = rho_l * c_l;
        let z_r = rho_r * c_r;
        let p_star = (z_r * p_l + z_l * p_r - z_l * z_r * (v_r_n - v_l_n)) / (z_l + z_r);
        let v_star = (z_l * v_l_n + z_r * v_r_n + (p_l - p_r)) / (z_l + z_r);
        (p_star, v_star)
    }

    /// Low-dissipation Riemann solver (pressure-based reconstruction).
    pub fn low_dissipation_riemann(
        &self,
        p_i: f64,
        p_j: f64,
        u_ij_n: f64,
        c_bar: f64,
        rho_bar: f64,
    ) -> (f64, f64) {
        let dp = p_i - p_j;
        let p_star = 0.5 * (p_i + p_j) - 0.5 * self.beta_r * rho_bar * c_bar * u_ij_n;
        let u_star = 0.5 * u_ij_n - dp / (2.0 * rho_bar * c_bar + 1e-30) * self.beta_r;
        (p_star, u_star)
    }

    /// Minmod limiter for reconstruction: minmod(a, b).
    pub fn minmod(a: f64, b: f64) -> f64 {
        if a * b <= 0.0 {
            0.0
        } else if a.abs() < b.abs() {
            a
        } else {
            b
        }
    }
}

// ============================================================================
// Time Integrator
// ============================================================================

/// WCSPH time integration schemes.
///
/// Provides Verlet, leapfrog, and predictor-corrector methods with
/// CFL-based adaptive time stepping.
#[derive(Debug, Clone)]
pub struct WcsphTimeIntegrator {
    /// Current time (s).
    pub time: f64,
    /// Current time step (s).
    pub dt: f64,
    /// Maximum allowed time step (s).
    pub dt_max: f64,
    /// CFL number (typically 0.2–0.4).
    pub cfl: f64,
}

impl WcsphTimeIntegrator {
    /// Create a new time integrator.
    pub fn new(dt: f64, dt_max: f64, cfl: f64) -> Self {
        Self {
            time: 0.0,
            dt,
            dt_max,
            cfl,
        }
    }

    /// Velocity Verlet: position update.
    ///
    /// x^{n+1} = x^n + dt * v^n + 0.5 * dt² * a^n
    pub fn verlet_position(x: f64, v: f64, a: f64, dt: f64) -> f64 {
        x + dt * v + 0.5 * dt * dt * a
    }

    /// Velocity Verlet: velocity update.
    ///
    /// v^{n+1} = v^n + 0.5 * dt * (a^n + a^{n+1})
    pub fn verlet_velocity(v: f64, a_old: f64, a_new: f64, dt: f64) -> f64 {
        v + 0.5 * dt * (a_old + a_new)
    }

    /// Leapfrog integration.
    ///
    /// v^{n+1/2} = v^{n-1/2} + dt * a^n
    /// x^{n+1}   = x^n + dt * v^{n+1/2}
    pub fn leapfrog_step(x: f64, v_half: f64, a: f64, dt: f64) -> (f64, f64) {
        let v_new = v_half + dt * a;
        let x_new = x + dt * v_new;
        (x_new, v_new)
    }

    /// Predictor step: predict position and velocity.
    pub fn predict(x: f64, v: f64, a: f64, dt: f64) -> (f64, f64) {
        let v_pred = v + 0.5 * dt * a;
        let x_pred = x + dt * v_pred;
        (x_pred, v_pred)
    }

    /// Corrector step: correct velocity with new acceleration.
    pub fn correct(v_pred: f64, a_old: f64, a_new: f64, dt: f64) -> f64 {
        v_pred + 0.5 * dt * (a_new - a_old)
    }

    /// Adaptive time step based on CFL and signal speed.
    pub fn adaptive_dt(&self, v_max: f64, c0: f64, h: f64) -> f64 {
        let dt_cfl = cfl_dt_wcsph(v_max, c0, h, self.cfl);
        dt_cfl.min(self.dt_max)
    }

    /// Advance simulation time.
    pub fn advance(&mut self) {
        self.time += self.dt;
    }
}

// ============================================================================
// Boundary Conditions
// ============================================================================

/// WCSPH boundary handling: dynamic particles, Lennard-Jones, Adami BC.
#[derive(Debug, Clone)]
pub struct WcsphBoundary {
    /// Lennard-Jones repulsion distance r₀.
    pub r0: f64,
    /// Lennard-Jones stiffness D (N/kg).
    pub lennard_jones_d: f64,
    /// Lennard-Jones exponent p1 (typically 4).
    pub p1: f64,
    /// Lennard-Jones exponent p2 (typically 2).
    pub p2: f64,
    /// Adami method: use pressure extrapolation.
    pub use_adami: bool,
}

impl WcsphBoundary {
    /// Create boundary with Lennard-Jones repulsion.
    pub fn lennard_jones(r0: f64, d: f64) -> Self {
        Self {
            r0,
            lennard_jones_d: d,
            p1: 4.0,
            p2: 2.0,
            use_adami: false,
        }
    }

    /// Create boundary using Adami method.
    pub fn adami(r0: f64) -> Self {
        Self {
            r0,
            lennard_jones_d: 0.0,
            p1: 4.0,
            p2: 2.0,
            use_adami: true,
        }
    }

    /// Lennard-Jones repulsive force magnitude.
    ///
    /// f = D * ((r₀/r)^p1 - (r₀/r)^p2) / r²  if r < r₀
    pub fn lennard_jones_force(&self, r: f64) -> f64 {
        if r >= self.r0 || r < 1e-14 {
            return 0.0;
        }
        let ratio = self.r0 / r;
        self.lennard_jones_d * (ratio.powf(self.p1) - ratio.powf(self.p2)) / (r * r)
    }

    /// Adami wall pressure: extrapolated from fluid particle pressure.
    ///
    /// p_w ≈ (Σ p_f * W_fw + Σ ρ_f * (g - a_w) · r_fw * W_fw) / Σ W_fw
    pub fn adami_wall_pressure(&self, sum_p_w: f64, sum_rho_g_r_w: f64, sum_w: f64) -> f64 {
        if sum_w < 1e-14 {
            0.0
        } else {
            (sum_p_w + sum_rho_g_r_w) / sum_w
        }
    }
}

// ============================================================================
// Multi-Phase
// ============================================================================

/// Multi-phase WCSPH with density ratio and surface tension.
///
/// Handles two immiscible fluids with different densities using
/// a color function and continuum surface stress formulation.
#[derive(Debug, Clone)]
pub struct WcsphMultiPhase {
    /// Density of phase 1 (kg/m³).
    pub rho1: f64,
    /// Density of phase 2 (kg/m³).
    pub rho2: f64,
    /// Surface tension coefficient σ (N/m).
    pub sigma: f64,
    /// Interface smoothing width (m).
    pub epsilon: f64,
}

impl WcsphMultiPhase {
    /// Create a two-phase WCSPH model.
    pub fn new(rho1: f64, rho2: f64, sigma: f64, epsilon: f64) -> Self {
        Self {
            rho1,
            rho2,
            sigma,
            epsilon,
        }
    }

    /// Density ratio r = ρ₂/ρ₁.
    pub fn density_ratio(&self) -> f64 {
        self.rho2 / self.rho1
    }

    /// Color function value for phase assignment (1.0 for phase 1, 0.0 for phase 2).
    pub fn color_function(&self, is_phase1: bool) -> f64 {
        if is_phase1 { 1.0 } else { 0.0 }
    }

    /// Surface tension force via Continuum Surface Stress (CSS).
    ///
    /// f_st = σ * κ * n_hat
    /// where κ is curvature and n_hat is the unit interface normal.
    pub fn surface_tension_force(&self, curvature: f64, normal: [f64; 3]) -> [f64; 3] {
        let f = self.sigma * curvature;
        [f * normal[0], f * normal[1], f * normal[2]]
    }

    /// Smoothed interface density: interpolated from color function.
    pub fn interface_density(&self, color: f64) -> f64 {
        color * self.rho1 + (1.0 - color) * self.rho2
    }
}

// ============================================================================
// Helper Functions
// ============================================================================

/// Tait EOS pressure: p = B * ((ρ/ρ₀)^γ - 1).
///
/// B = ρ₀ * c₀² / γ.
pub fn tait_pressure(rho: f64, rho0: f64, c0: f64, gamma: f64) -> f64 {
    let b = rho0 * c0 * c0 / gamma;
    b * ((rho / rho0).powf(gamma) - 1.0)
}

/// Sound speed from Tait EOS: c(ρ) = c₀ * (ρ/ρ₀)^((γ-1)/2).
pub fn tait_sound_speed(rho: f64, rho0: f64, c0: f64, gamma: f64) -> f64 {
    c0 * (rho / rho0).powf((gamma - 1.0) / 2.0)
}

/// CFL-limited time step for WCSPH.
///
/// dt = CFL * h / (c₀ + v_max)
pub fn cfl_dt_wcsph(v_max: f64, c0: f64, h: f64, cfl: f64) -> f64 {
    let signal_speed = c0 + v_max;
    if signal_speed < 1e-14 {
        return f64::MAX;
    }
    cfl * h / signal_speed
}

/// Renormalize density to suppress density diffusion errors.
///
/// Applies Shepard filtering: ρ_i^new = ρ_i / Σⱼ (mⱼ/ρⱼ) * W_ij
pub fn renormalize_density(rho: f64, shepherd_sum: f64) -> f64 {
    if shepherd_sum.abs() < 1e-14 {
        rho
    } else {
        rho / shepherd_sum
    }
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    // WcsphConfig tests
    #[test]
    fn test_wcsph_config_tait_b() {
        let cfg = WcsphConfig::water(1000.0, 14.14, 0.01, 1e-6);
        let b = cfg.tait_b();
        assert!((b - 1000.0 * 14.14 * 14.14 / 7.0).abs() < 1e-3);
    }

    #[test]
    fn test_wcsph_config_cfl_dt() {
        let cfg = WcsphConfig::water(1000.0, 14.14, 0.01, 1e-6);
        let dt = cfg.cfl_dt(0.25);
        assert!((dt - 0.25 * 0.01 / 14.14).abs() < 1e-10);
    }

    #[test]
    fn test_wcsph_config_default() {
        let cfg = WcsphConfig::default();
        assert_eq!(cfg.gamma, 7.0);
    }

    // TaitEos tests
    #[test]
    fn test_tait_eos_pressure_at_rho0() {
        let eos = TaitEos::new(1000.0, 14.14, 7.0);
        let p = eos.pressure(1000.0);
        assert!(p.abs() < 1e-6);
    }

    #[test]
    fn test_tait_eos_pressure_positive_compression() {
        let eos = TaitEos::new(1000.0, 14.14, 7.0);
        let p = eos.pressure(1010.0);
        assert!(p > 0.0);
    }

    #[test]
    fn test_tait_eos_pressure_negative_expansion() {
        let eos = TaitEos::new(1000.0, 14.14, 7.0);
        let p = eos.pressure(990.0);
        assert!(p < 0.0);
    }

    #[test]
    fn test_tait_eos_density_roundtrip() {
        let eos = TaitEos::new(1000.0, 14.14, 7.0);
        let rho = 1005.0;
        let p = eos.pressure(rho);
        let rho_back = eos.density(p);
        assert!((rho_back - rho).abs() < 1e-6);
    }

    #[test]
    fn test_tait_eos_sound_speed_at_rho0() {
        let eos = TaitEos::new(1000.0, 14.14, 7.0);
        let c = eos.sound_speed(1000.0);
        assert!((c - 14.14).abs() < 1e-10);
    }

    #[test]
    fn test_tait_eos_bulk_modulus_positive() {
        let eos = TaitEos::new(1000.0, 14.14, 7.0);
        assert!(eos.bulk_modulus(1000.0) > 0.0);
    }

    // MorrisEos tests
    #[test]
    fn test_morris_eos_zero_pressure_at_rho0() {
        let eos = MorrisEos::new(1000.0, 14.14);
        assert_eq!(eos.pressure(1000.0), 0.0);
    }

    #[test]
    fn test_morris_eos_density_roundtrip() {
        let eos = MorrisEos::new(1000.0, 14.14);
        let p = 100.0;
        let rho = eos.density(p);
        let p_back = eos.pressure(rho);
        assert!((p_back - p).abs() < 1e-10);
    }

    #[test]
    fn test_morris_eos_sound_speed() {
        let eos = MorrisEos::new(1000.0, 14.14);
        assert_eq!(eos.sound_speed(), 14.14);
    }

    // WcsphSolverAdvanced tests
    #[test]
    fn test_wcsph_solver_add_particle() {
        let mut solver = WcsphSolverAdvanced::new(WcsphConfig::default());
        solver.add_particle([0.0, 0.0, 0.0], [1.0, 0.0, 0.0]);
        assert_eq!(solver.num_particles(), 1);
    }

    #[test]
    fn test_wcsph_solver_update_pressures() {
        let mut solver = WcsphSolverAdvanced::new(WcsphConfig::default());
        solver.add_particle([0.0, 0.0, 0.0], [0.0, 0.0, 0.0]);
        solver.update_pressures();
        // At rest density, pressure should be ~0
        assert!(solver.pressures[0].abs() < 1e-3);
    }

    #[test]
    fn test_wcsph_solver_max_speed_zero() {
        let mut solver = WcsphSolverAdvanced::new(WcsphConfig::default());
        solver.add_particle([0.0, 0.0, 0.0], [0.0, 0.0, 0.0]);
        assert_eq!(solver.max_speed(), 0.0);
    }

    #[test]
    fn test_wcsph_solver_cfl_dt() {
        let mut solver = WcsphSolverAdvanced::new(WcsphConfig::default());
        solver.add_particle([0.0, 0.0, 0.0], [1.0, 0.0, 0.0]);
        let dt = solver.cfl_dt(0.25);
        assert!(dt > 0.0);
    }

    // ArtificialViscosity tests
    #[test]
    fn test_artificial_viscosity_approaching() {
        let av = ArtificialViscosity::new(0.01, 0.0, 0.01);
        // Approaching particles: vr < 0 → Π_ij > 0 (adds dissipation via pressure-like term)
        let pi = av.compute(
            [-1.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            1.0,
            14.0,
            14.0,
            1000.0,
            1000.0,
            0.01,
        );
        assert!(pi > 0.0); // Π_ij is positive for approaching particles
    }

    #[test]
    fn test_artificial_viscosity_separating() {
        let av = ArtificialViscosity::new(0.01, 0.0, 0.01);
        // Separating particles: vr > 0 → Π = 0
        let pi = av.compute(
            [1.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            1.0,
            14.0,
            14.0,
            1000.0,
            1000.0,
            0.01,
        );
        assert_eq!(pi, 0.0);
    }

    // XsphCorrection tests
    #[test]
    fn test_xsph_zero_relative_velocity() {
        let xsph = XsphCorrection::new(0.2);
        let c = xsph.correction_ij([1.0, 0.0, 0.0], [1.0, 0.0, 0.0], 1000.0, 1000.0, 1e-6, 1.0);
        assert_eq!(c, [0.0, 0.0, 0.0]);
    }

    #[test]
    fn test_xsph_corrected_velocity_no_correction() {
        let xsph = XsphCorrection::new(0.2);
        let v = xsph.corrected_velocity([1.0, 0.0, 0.0], [0.0, 0.0, 0.0]);
        assert_eq!(v, [1.0, 0.0, 0.0]);
    }

    // DeltaSph tests
    #[test]
    fn test_delta_sph_psi_same_density() {
        let ds = DeltaSph::new(0.1, 14.14, 0.01);
        let psi = ds.psi_ij(1000.0, 1000.0, [0.1, 0.0, 0.0], 0.01);
        assert_eq!(psi, [0.0, 0.0, 0.0]);
    }

    #[test]
    fn test_delta_sph_diffusion_nonzero() {
        let ds = DeltaSph::new(0.1, 14.14, 0.01);
        let d = ds.diffusion_ij(
            1000.0,
            1010.0,
            1e-6,
            [0.01, 0.0, 0.0],
            1e-4,
            [1.0, 0.0, 0.0],
        );
        assert!(d.is_finite());
    }

    // RiemannSph tests
    #[test]
    fn test_riemann_sph_same_state() {
        let rs = RiemannSph::new(1.0, 0);
        let (p_star, v_star) =
            rs.acoustic_riemann(1000.0, 1000.0, 0.0, 0.0, 1000.0, 1000.0, 14.14, 14.14);
        assert!((p_star - 1000.0).abs() < 1e-8);
        assert!(v_star.abs() < 1e-10);
    }

    #[test]
    fn test_riemann_sph_minmod_same_sign() {
        let a = 3.0_f64;
        let b = 5.0_f64;
        assert_eq!(RiemannSph::minmod(a, b), a);
    }

    #[test]
    fn test_riemann_sph_minmod_opposite_sign() {
        assert_eq!(RiemannSph::minmod(3.0, -5.0), 0.0);
    }

    // WcsphTimeIntegrator tests
    #[test]
    fn test_verlet_position_zero_accel() {
        let x = WcsphTimeIntegrator::verlet_position(1.0, 2.0, 0.0, 0.5);
        assert!((x - 2.0).abs() < 1e-14);
    }

    #[test]
    fn test_verlet_velocity_zero_accel() {
        let v = WcsphTimeIntegrator::verlet_velocity(1.0, 0.0, 0.0, 0.5);
        assert_eq!(v, 1.0);
    }

    #[test]
    fn test_leapfrog_step() {
        let (x, v) = WcsphTimeIntegrator::leapfrog_step(0.0, 1.0, 0.0, 0.1);
        assert!((x - 0.1).abs() < 1e-14);
        assert!((v - 1.0).abs() < 1e-14);
    }

    #[test]
    fn test_time_integrator_advance() {
        let mut ti = WcsphTimeIntegrator::new(1e-3, 1e-2, 0.25);
        ti.advance();
        assert!((ti.time - 1e-3).abs() < 1e-14);
    }

    // WcsphBoundary tests
    #[test]
    fn test_lennard_jones_no_force_outside() {
        let bc = WcsphBoundary::lennard_jones(0.01, 1e-3);
        assert_eq!(bc.lennard_jones_force(0.02), 0.0);
    }

    #[test]
    fn test_lennard_jones_repulsion_inside() {
        let bc = WcsphBoundary::lennard_jones(0.01, 1e-3);
        let f = bc.lennard_jones_force(0.005);
        assert!(f > 0.0);
    }

    #[test]
    fn test_adami_pressure_zero_sum() {
        let bc = WcsphBoundary::adami(0.01);
        let p = bc.adami_wall_pressure(0.0, 0.0, 0.0);
        assert_eq!(p, 0.0);
    }

    // WcsphMultiPhase tests
    #[test]
    fn test_multiphase_density_ratio() {
        let mp = WcsphMultiPhase::new(1000.0, 1.2, 0.07, 0.01);
        assert!((mp.density_ratio() - 1.2 / 1000.0).abs() < 1e-14);
    }

    #[test]
    fn test_multiphase_interface_density() {
        let mp = WcsphMultiPhase::new(1000.0, 1.2, 0.07, 0.01);
        let rho = mp.interface_density(1.0);
        assert!((rho - 1000.0).abs() < 1e-10);
        let rho2 = mp.interface_density(0.0);
        assert!((rho2 - 1.2).abs() < 1e-10);
    }

    #[test]
    fn test_multiphase_surface_tension_force() {
        let mp = WcsphMultiPhase::new(1000.0, 1.2, 0.07, 0.01);
        let f = mp.surface_tension_force(10.0, [0.0, 1.0, 0.0]);
        assert!((f[1] - 0.7).abs() < 1e-10);
    }

    // Helper function tests
    #[test]
    fn test_tait_pressure_at_rest() {
        let p = tait_pressure(1000.0, 1000.0, 14.14, 7.0);
        assert!(p.abs() < 1e-6);
    }

    #[test]
    fn test_tait_sound_speed_at_rest() {
        let c = tait_sound_speed(1000.0, 1000.0, 14.14, 7.0);
        assert!((c - 14.14).abs() < 1e-10);
    }

    #[test]
    fn test_cfl_dt_wcsph() {
        let dt = cfl_dt_wcsph(0.0, 14.14, 0.01, 0.25);
        assert!((dt - 0.25 * 0.01 / 14.14).abs() < 1e-14);
    }

    #[test]
    fn test_renormalize_density_normal() {
        let rho = renormalize_density(1000.0, 0.95);
        assert!((rho - 1000.0 / 0.95).abs() < 1e-8);
    }

    #[test]
    fn test_renormalize_density_zero_sum() {
        let rho = renormalize_density(1000.0, 0.0);
        assert_eq!(rho, 1000.0);
    }
}
