// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Compressible SPH for high-speed flows and shocks.
//!
//! This module implements compressible Smoothed Particle Hydrodynamics (SPH)
//! including Riemann-based methods and shock-capturing formulations:
//!
//! - [`CompressibleParticle`]: State variables for a compressible SPH particle
//! - [`CompressibleSph`]: Full compressible SPH solver with ideal-gas EOS
//!   (kick-drift-kick leapfrog, optional Monaghan artificial viscosity)
//! - [`KernelDimension`]: 1-D / 2-D / 3-D normalisation of the cubic-spline kernel
//! - [`ArtificialViscosityCs`]: Monaghan-style artificial viscosity
//! - [`GodunovSph`]: Godunov-type SPH using embedded Riemann solvers
//! - [`RiemannSolverType`]: Available Riemann solver variants
//! - [`StiffenedGasEos`]: Stiffened-gas equation of state for liquids/explosives
//! - [`rankine_hugoniot`]: Rankine-Hugoniot post-shock state (ρ₂, p₂, u₂)
//!
//! # Physics Background
//!
//! Compressible SPH solves the Euler / Navier-Stokes equations in Lagrangian
//! form. The artificial viscosity of Monaghan (1992) prevents particle
//! interpenetration and captures shocks. Godunov-SPH replaces the viscosity by
//! Riemann fluxes between particle pairs, yielding sharper shock resolution.
//!
//! ## Rankine-Hugoniot Jump Conditions
//! For a normal shock of Mach number M₁:
//!  ρ₂/ρ₁ = (γ+1)M₁² / ((γ-1)M₁² + 2)
//!  p₂/p₁ = (2γM₁² - (γ-1)) / (γ+1)
//!  u₂    = u₁ + (p₂-p₁) / (ρ₁·cs₁·M₁)  (post-shock velocity in lab frame)
//!
//! # References
//! - Monaghan, J.J. (1992). Smoothed particle hydrodynamics.
//!   *Ann. Rev. Astron. Astrophys.*, 30, 543–574.
//! - Parshikov, A.N. & Medin, S.A. (2002). Smoothed particle hydrodynamics
//!   using interparticle contact algorithms. *J. Comput. Phys.*, 180, 358–382.
//! - Harlow, F.H. & Amsden, A.A. (1971). *Fluid Dynamics*. LANL Monograph.

// ============================================================================
// Math helpers
// ============================================================================

#[inline]
fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[inline]
fn len3(v: [f64; 3]) -> f64 {
    dot3(v, v).sqrt()
}

#[inline]
fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

#[inline]
fn add3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

#[inline]
fn scale3(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

// ============================================================================
// CompressibleParticle
// ============================================================================

/// State variables for a single compressible SPH particle.
///
/// Stores the full thermodynamic and kinematic state required by the
/// compressible Euler equations in Lagrangian SPH form.
#[derive(Clone, Debug)]
pub struct CompressibleParticle {
    /// Position vector **x** (m).
    pub pos: [f64; 3],
    /// Velocity vector **v** (m/s).
    pub vel: [f64; 3],
    /// Mass density ρ (kg/m³).
    pub density: f64,
    /// Thermodynamic pressure p (Pa).
    pub pressure: f64,
    /// Specific internal energy e (J/kg).
    pub energy: f64,
    /// Local sound speed c_s (m/s).
    pub sound_speed: f64,
    /// Particle mass m (kg).
    pub mass: f64,
    /// Smoothing length h (m).
    pub h: f64,
}

impl CompressibleParticle {
    /// Create a new compressible particle.
    pub fn new(
        pos: [f64; 3],
        vel: [f64; 3],
        density: f64,
        pressure: f64,
        energy: f64,
        sound_speed: f64,
        mass: f64,
        h: f64,
    ) -> Self {
        Self {
            pos,
            vel,
            density,
            pressure,
            energy,
            sound_speed,
            mass,
            h,
        }
    }

    /// Specific total energy E = e + ½|v|² (J/kg).
    pub fn total_specific_energy(&self) -> f64 {
        self.energy + 0.5 * dot3(self.vel, self.vel)
    }

    /// Total particle energy (J).
    pub fn total_energy(&self) -> f64 {
        self.mass * self.total_specific_energy()
    }

    /// Mach number M = |v| / c_s.
    pub fn mach_number(&self) -> f64 {
        if self.sound_speed.abs() < 1e-300 {
            return f64::INFINITY;
        }
        len3(self.vel) / self.sound_speed
    }
}

// ============================================================================
// CompressibleSph
// ============================================================================

/// Spatial dimension used to normalise the cubic-spline SPH kernel.
///
/// Particles always carry 3-D positions, but the kernel normalisation has to
/// match the dimensionality of the flow: particles arranged on a line need
/// the 1-D constant, a planar flow the 2-D one.  With the wrong constant the
/// summation density — and with it every pressure and force — is off by a
/// factor of order `h^(d_used − d)`.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum KernelDimension {
    /// One-dimensional flow (particles on a line).
    One,
    /// Two-dimensional (planar) flow.
    Two,
    /// Three-dimensional flow.
    #[default]
    Three,
}

impl KernelDimension {
    /// Number of spatial dimensions.
    pub fn count(self) -> usize {
        match self {
            KernelDimension::One => 1,
            KernelDimension::Two => 2,
            KernelDimension::Three => 3,
        }
    }

    /// Normalisation `σ_d(h)` of the cubic spline written as
    /// `W = σ_d (2/3 − q² + q³/2)` for `q < 1` and `σ_d (2 − q)³ / 6` for
    /// `1 ≤ q < 2`, so that `∫ W dV = 1` in `d` dimensions:
    /// `1/h` (1-D), `15/(7π h²)` (2-D) and `3/(2π h³)` (3-D).
    pub fn normalisation(self, h: f64) -> f64 {
        match self {
            KernelDimension::One => 1.0 / h,
            KernelDimension::Two => 15.0 / (7.0 * std::f64::consts::PI * h * h),
            KernelDimension::Three => 3.0 / (2.0 * std::f64::consts::PI * h * h * h),
        }
    }
}

/// Compressible SPH solver for high-speed and shock-dominated flows.
///
/// Uses an ideal-gas EOS p = (γ-1)·ρ·e and a second-order kick-drift-kick
/// leapfrog integrator.  The momentum and energy equations use the symmetric
/// SPH form
///
/// ```text
/// dv_i/dt = −Σ_j m_j (p_i/ρ_i² + p_j/ρ_j² + Π_ij) ∇_i W_ij
/// de_i/dt =  Σ_j m_j (p_i/ρ_i² + ½ Π_ij) (v_i − v_j) · ∇_i W_ij
/// ```
///
/// where `Π_ij` is the optional Monaghan artificial viscosity.  Every pair
/// interaction is antisymmetric, so linear momentum and the total energy
/// `Σ m (e + ½ v²)` are conserved by the force evaluation.
pub struct CompressibleSph {
    /// Collection of all particles.
    pub particles: Vec<CompressibleParticle>,
    /// Ratio of specific heats γ.
    pub gamma: f64,
    /// Global smoothing length h (m).
    pub h: f64,
    /// Dimensionality of the flow; selects the kernel normalisation.
    pub dim: KernelDimension,
    /// Monaghan artificial viscosity for shock capturing (`None` = inviscid).
    pub viscosity: Option<ArtificialViscosityCs>,
}

impl CompressibleSph {
    /// Create a new solver with no particles.
    ///
    /// The kernel is 3-D ([`KernelDimension::Three`]) and the flow is inviscid;
    /// use [`with_dimension`](Self::with_dimension) and
    /// [`with_viscosity`](Self::with_viscosity) to change that.
    pub fn new(gamma: f64, h: f64) -> Self {
        Self {
            particles: Vec::new(),
            gamma,
            h,
            dim: KernelDimension::Three,
            viscosity: None,
        }
    }

    /// Select the dimensionality of the flow (kernel normalisation).
    pub fn with_dimension(mut self, dim: KernelDimension) -> Self {
        self.dim = dim;
        self
    }

    /// Enable Monaghan artificial viscosity for shock capturing.
    pub fn with_viscosity(mut self, viscosity: ArtificialViscosityCs) -> Self {
        self.viscosity = Some(viscosity);
        self
    }

    /// Add a particle.
    pub fn add_particle(&mut self, p: CompressibleParticle) {
        self.particles.push(p);
    }

    /// Number of particles.
    pub fn n_particles(&self) -> usize {
        self.particles.len()
    }

    /// Ideal-gas equation of state: p = (γ - 1) · ρ · e.
    ///
    /// # Arguments
    /// * `density` – mass density ρ (kg/m³)
    /// * `energy`  – specific internal energy e (J/kg)
    pub fn ideal_gas_eos(&self, density: f64, energy: f64) -> f64 {
        (self.gamma - 1.0) * density * energy
    }

    /// Update pressure and sound speed for every particle.
    fn update_thermodynamics(&mut self) {
        let gamma = self.gamma;
        for p in &mut self.particles {
            p.pressure = (gamma - 1.0) * p.density * p.energy;
            p.sound_speed = if p.density > 0.0 && p.pressure > 0.0 {
                (gamma * p.pressure / p.density).sqrt()
            } else {
                0.0
            };
        }
    }

    /// Perform one explicit time step of size `dt` (s).
    ///
    /// Kick-drift-kick leapfrog:
    /// 1. densities (summation) and pressure / sound speed are brought up to
    ///    date for the current positions, then the accelerations and energy
    ///    rates are evaluated;
    /// 2. half-kick of velocities and specific energies, drift of the
    ///    positions with the half-step velocity;
    /// 3. density and thermodynamics are recomputed at the new positions and
    ///    the rates are re-evaluated — the velocity-dependent terms (viscosity,
    ///    compression work) use the velocity and energy predicted to the end of
    ///    the step;
    /// 4. closing half-kick from the half-step state.
    ///
    /// Nothing happens for a non-positive or non-finite `dt`.
    pub fn step(&mut self, dt: f64) {
        if !dt.is_finite() || dt <= 0.0 || self.particles.is_empty() {
            return;
        }
        let half = 0.5 * dt;

        self.density_update();
        self.update_thermodynamics();
        let (acc_start, de_start) = self.compute_rates();

        for (p, (a, de)) in self
            .particles
            .iter_mut()
            .zip(acc_start.iter().zip(de_start.iter()))
        {
            p.vel = add3(p.vel, scale3(*a, half));
            p.energy = (p.energy + half * de).max(0.0);
        }
        for p in &mut self.particles {
            p.pos = add3(p.pos, scale3(p.vel, dt));
        }

        let half_step: Vec<([f64; 3], f64)> =
            self.particles.iter().map(|p| (p.vel, p.energy)).collect();
        for (p, (a, de)) in self
            .particles
            .iter_mut()
            .zip(acc_start.iter().zip(de_start.iter()))
        {
            p.vel = add3(p.vel, scale3(*a, half));
            p.energy = (p.energy + half * de).max(0.0);
        }
        self.density_update();
        self.update_thermodynamics();
        let (acc_end, de_end) = self.compute_rates();

        for (p, ((vel_half, energy_half), (a, de))) in self
            .particles
            .iter_mut()
            .zip(half_step.iter().zip(acc_end.iter().zip(de_end.iter())))
        {
            p.vel = add3(*vel_half, scale3(*a, half));
            p.energy = (*energy_half + half * de).max(0.0);
        }
        self.update_thermodynamics();
    }

    /// Update densities using the SPH continuity equation (summation density).
    ///
    /// ρ_i = Σ_j m_j · W(|x_i - x_j|, h)
    pub fn density_update(&mut self) {
        let h = self.h;
        let dim = self.dim;
        let positions: Vec<[f64; 3]> = self.particles.iter().map(|p| p.pos).collect();
        let masses: Vec<f64> = self.particles.iter().map(|p| p.mass).collect();
        for (p_i, pos_i) in self.particles.iter_mut().zip(positions.iter()) {
            let rho: f64 = positions
                .iter()
                .zip(masses.iter())
                .map(|(pos_j, &mj)| {
                    let r = len3(sub3(*pos_i, *pos_j));
                    mj * cubic_spline_kernel(r, h, dim)
                })
                .sum::<f64>()
                .max(1e-10);
            p_i.density = rho;
        }
    }

    /// Accelerations `dv/dt` and specific-energy rates `de/dt` of all
    /// particles for the current state (positions, velocities, densities,
    /// pressures and sound speeds are taken as they are).
    ///
    /// Each unordered pair is visited once and contributes equal and opposite
    /// momentum changes, so `Σ m_i a_i = 0` up to rounding.
    fn compute_rates(&self) -> (Vec<[f64; 3]>, Vec<f64>) {
        let n = self.particles.len();
        let h = self.h;
        let mut acc = vec![[0.0_f64; 3]; n];
        let mut de_dt = vec![0.0_f64; n];
        for i in 0..n {
            let pi = &self.particles[i];
            if pi.density.abs() < 1e-300 {
                continue;
            }
            let press_i = pi.pressure / (pi.density * pi.density);
            for j in (i + 1)..n {
                let pj = &self.particles[j];
                if pj.density.abs() < 1e-300 {
                    continue;
                }
                let rij = sub3(pi.pos, pj.pos);
                let r = len3(rij);
                if r <= 1e-300 {
                    continue;
                }
                let dw = cubic_spline_kernel_grad(r, h, self.dim);
                if dw == 0.0 {
                    continue;
                }
                // ∇_i W_ij = (dW/dr) r̂_ij  (points along x_i − x_j)
                let grad_w = scale3(rij, dw / r);
                let press_j = pj.pressure / (pj.density * pj.density);
                let visc = self
                    .viscosity
                    .as_ref()
                    .map_or(0.0, |av| av.viscosity_term(pi, pj, h));
                let coeff = press_i + press_j + visc;
                acc[i] = sub3(acc[i], scale3(grad_w, pj.mass * coeff));
                acc[j] = add3(acc[j], scale3(grad_w, pi.mass * coeff));
                let v_dot_grad_w = dot3(sub3(pi.vel, pj.vel), grad_w);
                de_dt[i] += pj.mass * (press_i + 0.5 * visc) * v_dot_grad_w;
                de_dt[j] += pi.mass * (press_j + 0.5 * visc) * v_dot_grad_w;
            }
        }
        (acc, de_dt)
    }

    /// Advance velocities by `dt` with the SPH momentum equation
    /// (pressure gradient plus the optional artificial viscosity):
    ///
    /// d**v**_i/dt = -Σ_j m_j (p_i/ρ_i² + p_j/ρ_j² + Π_ij) ∇_i W_ij
    ///
    /// This is an explicit Euler kick with the densities, pressures and
    /// positions currently stored in the particles; call
    /// [`density_update`](Self::density_update) and
    /// [`compute_sound_speed`](Self::compute_sound_speed) first so they are
    /// consistent.  [`step`](Self::step) composes these building blocks into a
    /// second-order integrator.
    pub fn momentum_update(&mut self, dt: f64) {
        let (acc, _) = self.compute_rates();
        for (p, a) in self.particles.iter_mut().zip(acc.iter()) {
            p.vel = add3(p.vel, scale3(*a, dt));
        }
    }

    /// Advance specific internal energies by `dt` with the SPH energy equation:
    ///
    /// de_i/dt = Σ_j m_j (p_i/ρ_i² + ½ Π_ij) (**v**_i - **v**_j) · ∇_i W_ij
    ///
    /// Explicit Euler update with the state currently stored in the particles
    /// (see [`momentum_update`](Self::momentum_update)); energies stay
    /// non-negative.
    pub fn energy_update(&mut self, dt: f64) {
        let (_, de_dt) = self.compute_rates();
        for (p, de) in self.particles.iter_mut().zip(de_dt.iter()) {
            p.energy = (p.energy + de * dt).max(0.0);
        }
    }

    /// Recompute the pressure and sound speed c_s = √(γ·p/ρ) of every particle
    /// from its current density and specific internal energy.
    pub fn compute_sound_speed(&mut self) {
        self.update_thermodynamics();
    }

    /// Total mechanical energy E_tot = Σ m_i · (e_i + ½|v_i|²).
    pub fn total_energy(&self) -> f64 {
        self.particles.iter().map(|p| p.total_energy()).sum()
    }

    /// Total linear momentum Σ m_i **v**_i.
    pub fn total_momentum(&self) -> [f64; 3] {
        self.particles
            .iter()
            .fold([0.0; 3], |acc, p| add3(acc, scale3(p.vel, p.mass)))
    }

    /// Maximum sound speed among all particles (m/s).
    pub fn max_sound_speed(&self) -> f64 {
        self.particles
            .iter()
            .map(|p| p.sound_speed)
            .fold(0.0_f64, f64::max)
    }

    /// CFL time step: dt = CFL · h / (c_s,max + |v|_max).
    ///
    /// The estimate uses the sound speed and the particle speeds only; it does
    /// not include the signal speed of the artificial viscosity.  For strong
    /// shocks or supersonic collisions with viscosity enabled, use a CFL
    /// number of about 0.1 or smaller to keep the total-energy error small.
    pub fn dt_cfl(&self, cfl: f64) -> f64 {
        let cs_max = self.max_sound_speed();
        let v_max = self
            .particles
            .iter()
            .map(|p| len3(p.vel))
            .fold(0.0_f64, f64::max);
        let denom = cs_max + v_max;
        if denom < 1e-300 {
            return f64::INFINITY;
        }
        cfl * self.h / denom
    }
}

// ============================================================================
// SPH kernel functions (cubic spline, 1-D / 2-D / 3-D normalisation)
// ============================================================================

/// Cubic-spline SPH kernel W(r, h) with the normalisation of `dim`
/// (support radius 2h).
fn cubic_spline_kernel(r: f64, h: f64, dim: KernelDimension) -> f64 {
    if h <= 0.0 {
        return 0.0;
    }
    let q = r / h;
    let norm = dim.normalisation(h);
    if q < 1.0 {
        norm * (2.0 / 3.0 - q * q + 0.5 * q * q * q)
    } else if q < 2.0 {
        let t = 2.0 - q;
        norm * (1.0 / 6.0) * t * t * t
    } else {
        0.0
    }
}

/// Radial derivative dW/dr of the cubic-spline kernel (non-positive).
fn cubic_spline_kernel_grad(r: f64, h: f64, dim: KernelDimension) -> f64 {
    if h <= 0.0 {
        return 0.0;
    }
    let q = r / h;
    let norm = dim.normalisation(h) / h;
    if q < 1.0 {
        norm * (-2.0 * q + 1.5 * q * q)
    } else if q < 2.0 {
        let t = 2.0 - q;
        -norm * 0.5 * t * t
    } else {
        0.0
    }
}

// ============================================================================
// ArtificialViscosityCs
// ============================================================================

/// Monaghan (1992) artificial viscosity for shock capturing in compressible SPH.
///
/// The viscosity term Π_ij is added to the pressure forces between particle pairs
/// that are approaching each other (μ_ij < 0).
pub struct ArtificialViscosityCs {
    /// Linear viscosity coefficient α (typically 1.0).
    pub alpha: f64,
    /// Quadratic viscosity coefficient β (typically 2.0 for strong shocks).
    pub beta: f64,
    /// Small positive ε to prevent singularity (fraction of h²).
    pub epsilon: f64,
}

impl ArtificialViscosityCs {
    /// Create a new artificial viscosity operator.
    pub fn new(alpha: f64, beta: f64) -> Self {
        Self {
            alpha,
            beta,
            epsilon: 0.01,
        }
    }

    /// Artificial viscosity term Π_ij between particles `pi` and `pj`.
    ///
    /// Π_ij = (-α·c̄_ij·μ_ij + β·μ_ij²) / ρ̄_ij  when **v**_ij · **r**_ij < 0,
    /// and 0 otherwise.
    ///
    /// # Arguments
    /// * `pi` – particle i
    /// * `pj` – particle j
    /// * `h`  – smoothing length (m)
    pub fn viscosity_term(
        &self,
        pi: &CompressibleParticle,
        pj: &CompressibleParticle,
        h: f64,
    ) -> f64 {
        let rij = sub3(pi.pos, pj.pos);
        let vij = sub3(pi.vel, pj.vel);
        let vr = dot3(vij, rij);
        // Only apply to approaching pairs
        if vr >= 0.0 {
            return 0.0;
        }
        let r2 = dot3(rij, rij);
        let mu_ij = h * vr / (r2 + self.epsilon * h * h);
        let c_bar = 0.5 * (pi.sound_speed + pj.sound_speed);
        let rho_bar = 0.5 * (pi.density + pj.density);
        if rho_bar.abs() < 1e-300 {
            return 0.0;
        }
        (-self.alpha * c_bar * mu_ij + self.beta * mu_ij * mu_ij) / rho_bar
    }

    /// Signal velocity: v_sig = c_i + c_j - β · v_ij · r̂_ij.
    ///
    /// Used as an alternative signal-based viscosity formulation.
    pub fn signal_velocity(&self, pi: &CompressibleParticle, pj: &CompressibleParticle) -> f64 {
        let rij = sub3(pi.pos, pj.pos);
        let r = len3(rij);
        let vij = sub3(pi.vel, pj.vel);
        let w_ij = if r > 1e-300 {
            dot3(vij, scale3(rij, 1.0 / r)).min(0.0)
        } else {
            0.0
        };
        pi.sound_speed + pj.sound_speed - self.beta * w_ij
    }
}

// ============================================================================
// RiemannSolverType
// ============================================================================

/// Available Riemann solver variants for Godunov-SPH.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum RiemannSolverType {
    /// Roe approximate Riemann solver.
    Roe,
    /// HLLC (Harten-Lax-van Leer-Contact) solver.
    HllC,
    /// Exact iterative Riemann solver.
    Exact,
}

// ============================================================================
// GodunovSph
// ============================================================================

/// Godunov-type SPH solver that replaces artificial viscosity with Riemann fluxes.
///
/// Each particle-pair interaction is treated as a local 1-D Riemann problem,
/// giving sharper shock resolution and naturally satisfying the entropy condition.
pub struct GodunovSph {
    /// Which Riemann solver to use at each interface.
    pub riemann_solver: RiemannSolverType,
    /// Adiabatic index γ.
    pub gamma: f64,
}

impl GodunovSph {
    /// Create a new Godunov-SPH solver.
    pub fn new(riemann_solver: RiemannSolverType, gamma: f64) -> Self {
        Self {
            riemann_solver,
            gamma,
        }
    }

    /// Solve the 1-D Riemann problem given left/right states.
    ///
    /// Returns `(p*, u*, rho*)` — the star-region pressure, velocity, and density.
    ///
    /// # Arguments
    /// * `left`  – (ρ_L, p_L, u_L) left state
    /// * `right` – (ρ_R, p_R, u_R) right state
    pub fn solve_riemann(&self, left: (f64, f64, f64), right: (f64, f64, f64)) -> (f64, f64, f64) {
        let (rho_l, p_l, u_l) = left;
        let (rho_r, p_r, u_r) = right;
        match self.riemann_solver {
            RiemannSolverType::Roe => self.roe_solve(rho_l, p_l, u_l, rho_r, p_r, u_r),
            RiemannSolverType::HllC => self.hllc_solve(rho_l, p_l, u_l, rho_r, p_r, u_r),
            RiemannSolverType::Exact => self.exact_solve(rho_l, p_l, u_l, rho_r, p_r, u_r),
        }
    }

    /// Godunov numerical flux F* for the conservative variables.
    ///
    /// Returns `(F_rho, F_rho_u, F_E)` — the intercell fluxes.
    pub fn godunov_flux(&self, left: (f64, f64, f64), right: (f64, f64, f64)) -> (f64, f64, f64) {
        let (p_star, u_star, rho_star) = self.solve_riemann(left, right);
        let e_star = if rho_star > 0.0 && self.gamma > 1.0 {
            p_star / ((self.gamma - 1.0) * rho_star)
        } else {
            0.0
        };
        let e_tot_star = rho_star * (e_star + 0.5 * u_star * u_star);
        let f_rho = rho_star * u_star;
        let f_rho_u = rho_star * u_star * u_star + p_star;
        let f_e = u_star * (e_tot_star + p_star);
        (f_rho, f_rho_u, f_e)
    }

    // ---- Roe approximate solver (linearised) --------------------------------

    fn roe_solve(
        &self,
        rho_l: f64,
        p_l: f64,
        u_l: f64,
        rho_r: f64,
        p_r: f64,
        u_r: f64,
    ) -> (f64, f64, f64) {
        let g = self.gamma;
        let cs_l = if rho_l > 0.0 && p_l > 0.0 {
            (g * p_l / rho_l).sqrt()
        } else {
            0.0
        };
        let cs_r = if rho_r > 0.0 && p_r > 0.0 {
            (g * p_r / rho_r).sqrt()
        } else {
            0.0
        };
        // Roe-averaged velocity and enthalpy
        let sqrt_l = rho_l.sqrt();
        let sqrt_r = rho_r.sqrt();
        let denom = sqrt_l + sqrt_r;
        let u_roe = if denom > 1e-300 {
            (sqrt_l * u_l + sqrt_r * u_r) / denom
        } else {
            0.0
        };
        // Simple average for pressure and density in star state
        let u_star = 0.5 * (u_l + u_r)
            - 0.5 * (p_r - p_l) / (0.5 * (rho_l + rho_r) * (cs_l + cs_r + 1e-300));
        let p_star = 0.5 * (p_l + p_r) - 0.5 * (u_r - u_l) * 0.5 * (rho_l + rho_r) * (cs_l + cs_r);
        let p_star = p_star.max(0.0);
        let rho_star = if u_roe >= 0.0 { rho_l } else { rho_r };
        (p_star, u_star, rho_star)
    }

    // ---- HLLC solver --------------------------------------------------------

    fn hllc_solve(
        &self,
        rho_l: f64,
        p_l: f64,
        u_l: f64,
        rho_r: f64,
        p_r: f64,
        u_r: f64,
    ) -> (f64, f64, f64) {
        let g = self.gamma;
        let cs_l = if rho_l > 0.0 && p_l > 0.0 {
            (g * p_l / rho_l).sqrt()
        } else {
            0.0
        };
        let cs_r = if rho_r > 0.0 && p_r > 0.0 {
            (g * p_r / rho_r).sqrt()
        } else {
            0.0
        };
        // HLLC wave speed estimates
        let s_l = u_l - cs_l;
        let s_r = u_r + cs_r;
        // Contact wave speed
        let denom = rho_r * (s_r - u_r) - rho_l * (s_l - u_l);
        let s_star = if denom.abs() > 1e-300 {
            (p_r - p_l + rho_l * u_l * (s_l - u_l) - rho_r * u_r * (s_r - u_r)) / denom
        } else {
            0.5 * (u_l + u_r)
        };
        // HLLC pressure
        let p_star = 0.5
            * ((p_l + rho_l * (s_l - u_l) * (s_star - u_l))
                + (p_r + rho_r * (s_r - u_r) * (s_star - u_r)));
        let p_star = p_star.max(0.0);
        let rho_star = if s_star >= 0.0 { rho_l } else { rho_r };
        (p_star, s_star, rho_star)
    }

    // ---- Exact iterative solver (Newton-Raphson on p*) ----------------------

    fn exact_solve(
        &self,
        rho_l: f64,
        p_l: f64,
        u_l: f64,
        rho_r: f64,
        p_r: f64,
        u_r: f64,
    ) -> (f64, f64, f64) {
        let g = self.gamma;
        let cs_l = if rho_l > 0.0 && p_l > 0.0 {
            (g * p_l / rho_l).sqrt()
        } else {
            1e-6
        };
        let cs_r = if rho_r > 0.0 && p_r > 0.0 {
            (g * p_r / rho_r).sqrt()
        } else {
            1e-6
        };

        // Initial guess via acoustic approximation
        let mut p =
            (0.5 * (p_l + p_r) - 0.125 * (u_r - u_l) * (rho_l + rho_r) * (cs_l + cs_r)).max(1e-10);

        // Newton-Raphson iteration (up to 20 steps)
        for _ in 0..20 {
            let (fl, fl_p) = wave_function(p, p_l, rho_l, cs_l, g);
            let (fr, fr_p) = wave_function(p, p_r, rho_r, cs_r, g);
            let delta_u = u_r - u_l;
            let f = fl + fr + delta_u;
            let fp = fl_p + fr_p;
            if fp.abs() < 1e-300 {
                break;
            }
            let dp = -f / fp;
            p = (p + dp).max(1e-10);
            if dp.abs() < 1e-8 * p {
                break;
            }
        }

        let (fl, _) = wave_function(p, p_l, rho_l, cs_l, g);
        let (fr, _) = wave_function(p, p_r, rho_r, cs_r, g);
        let u_star = 0.5 * (u_l + u_r) + 0.5 * (fr - fl);

        // Post-shock density (shock side)
        let rho_star = if u_star >= 0.0 {
            shock_or_rarefaction_density(p, p_l, rho_l, g)
        } else {
            shock_or_rarefaction_density(p, p_r, rho_r, g)
        };
        (p, u_star, rho_star)
    }
}

/// Wave function f(p, p_k, ρ_k, c_k) and its derivative for the exact Riemann solver.
fn wave_function(p: f64, p_k: f64, rho_k: f64, cs_k: f64, g: f64) -> (f64, f64) {
    if p > p_k {
        // Shock branch
        let a_k = 2.0 / ((g + 1.0) * rho_k);
        let b_k = (g - 1.0) / (g + 1.0) * p_k;
        let denom = (a_k * (p + b_k)).sqrt();
        let f = (p - p_k) * (a_k / (p + b_k)).sqrt();
        let fp = if denom.abs() > 1e-300 {
            (1.0 - 0.5 * (p - p_k) / (p + b_k)) / denom
        } else {
            0.0
        };
        (f, fp)
    } else {
        // Rarefaction branch
        let exp = 2.0 * g / (g - 1.0);
        let ratio = if p_k > 1e-300 {
            (p / p_k).powf((g - 1.0) / (2.0 * g))
        } else {
            1.0
        };
        let f = 2.0 * cs_k / (g - 1.0) * (ratio - 1.0);
        let fp = 1.0 / (rho_k * cs_k) * (p / p_k).powf(-1.0 / exp);
        (f, fp)
    }
}

/// Post-shock or post-rarefaction density from star-region pressure.
fn shock_or_rarefaction_density(p_star: f64, p_k: f64, rho_k: f64, g: f64) -> f64 {
    if rho_k <= 0.0 {
        return rho_k;
    }
    if p_star > p_k {
        // Shock
        let ratio = p_star / p_k;
        rho_k * ((g + 1.0) * ratio + (g - 1.0)) / ((g - 1.0) * ratio + (g + 1.0))
    } else {
        // Rarefaction
        if p_k > 1e-300 {
            rho_k * (p_star / p_k).powf(1.0 / g)
        } else {
            rho_k
        }
    }
}

// ============================================================================
// StiffenedGasEos
// ============================================================================

/// Stiffened-gas equation of state for liquids, condensed matter, or explosives.
///
/// p = (γ - 1) · ρ · e - γ · p_∞
///
/// where p_∞ is the stiffness pressure (Pa). For p_∞ = 0 this reduces to the
/// ideal-gas EOS.
pub struct StiffenedGasEos {
    /// Ratio of specific heats γ.
    pub gamma: f64,
    /// Stiffness pressure p_∞ (Pa).
    pub p_inf: f64,
}

impl StiffenedGasEos {
    /// Create a new stiffened-gas EOS.
    pub fn new(gamma: f64, p_inf: f64) -> Self {
        Self { gamma, p_inf }
    }

    /// Pressure from density and specific internal energy.
    ///
    /// p = (γ - 1) · ρ · e - γ · p_∞
    ///
    /// # Arguments
    /// * `rho` – density ρ (kg/m³)
    /// * `e`   – specific internal energy e (J/kg)
    pub fn pressure(&self, rho: f64, e: f64) -> f64 {
        (self.gamma - 1.0) * rho * e - self.gamma * self.p_inf
    }

    /// Sound speed c_s = √(γ · (p + p_∞) / ρ).
    ///
    /// # Arguments
    /// * `rho` – density ρ (kg/m³)
    /// * `p`   – pressure p (Pa)
    pub fn sound_speed(&self, rho: f64, p: f64) -> f64 {
        if rho <= 0.0 {
            return 0.0;
        }
        let arg = self.gamma * (p + self.p_inf) / rho;
        if arg <= 0.0 {
            return 0.0;
        }
        arg.sqrt()
    }

    /// Temperature T = (p + p_∞) / (ρ · c_v · (γ - 1)).
    ///
    /// # Arguments
    /// * `rho` – density ρ (kg/m³)
    /// * `p`   – pressure p (Pa)
    pub fn temperature(&self, rho: f64, p: f64) -> f64 {
        if rho.abs() < 1e-300 || (self.gamma - 1.0).abs() < 1e-300 {
            return 0.0;
        }
        // Assume c_v = 1 J/(kg·K) for dimensionless output
        (p + self.p_inf) / (rho * (self.gamma - 1.0))
    }
}

// ============================================================================
// rankine_hugoniot
// ============================================================================

/// Post-shock state from the Rankine-Hugoniot jump conditions.
///
/// Given an upstream state (ρ₁, p₁, u₁) and shock Mach number M > 1, returns
/// the downstream state (ρ₂, p₂, u₂) in the lab frame.
///
/// Jump relations for a normal shock:
/// - ρ₂ = ρ₁ · (γ+1)·M² / ((γ-1)·M² + 2)
/// - p₂ = p₁ · (2γ·M² - (γ-1)) / (γ+1)
/// - u₂ = u₁ - cs₁ · (M - 1/M) · 2 / (γ+1)  (in the shock frame, then shifted)
///
/// # Arguments
/// * `rho1`  – upstream density ρ₁ (kg/m³)
/// * `p1`    – upstream pressure p₁ (Pa)
/// * `u1`    – upstream velocity u₁ (m/s)
/// * `mach`  – shock Mach number M ≥ 1
/// * `gamma` – adiabatic index γ
///
/// # Returns
/// `(rho2, p2, u2)` — the post-shock state.
pub fn rankine_hugoniot(rho1: f64, p1: f64, u1: f64, mach: f64, gamma: f64) -> (f64, f64, f64) {
    let m = mach.max(1.0); // enforce supersonic
    let m2 = m * m;
    let gp1 = gamma + 1.0;
    let gm1 = gamma - 1.0;

    let rho2 = rho1 * gp1 * m2 / (gm1 * m2 + 2.0);
    let p2 = p1 * (2.0 * gamma * m2 - gm1) / gp1;
    // Sound speed upstream
    let cs1 = if rho1 > 0.0 && p1 > 0.0 {
        (gamma * p1 / rho1).sqrt()
    } else {
        0.0
    };
    // Post-shock velocity (lab frame, shock moving in +x)
    let u2 = u1 - cs1 * (m - 1.0 / m) * 2.0 / gp1;

    (rho2, p2, u2)
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    // ---- CompressibleParticle -----------------------------------------------

    #[test]
    fn test_particle_total_specific_energy() {
        let p = CompressibleParticle::new(
            [0.0; 3],
            [3.0, 4.0, 0.0],
            1.0,
            1e5,
            2.5e5,
            340.0,
            0.001,
            0.01,
        );
        let ke = 0.5 * (9.0 + 16.0);
        assert!((p.total_specific_energy() - (2.5e5 + ke)).abs() < 1e-6);
    }

    #[test]
    fn test_particle_mach_number_subsonic() {
        let p = CompressibleParticle::new(
            [0.0; 3],
            [100.0, 0.0, 0.0],
            1.2,
            1e5,
            2.0e5,
            343.0,
            0.001,
            0.01,
        );
        assert!(p.mach_number() < 1.0);
    }

    #[test]
    fn test_particle_mach_number_supersonic() {
        let p = CompressibleParticle::new(
            [0.0; 3],
            [700.0, 0.0, 0.0],
            1.2,
            1e5,
            2.0e5,
            343.0,
            0.001,
            0.01,
        );
        assert!(p.mach_number() > 1.0);
    }

    // ---- CompressibleSph ----------------------------------------------------

    #[test]
    fn test_csph_ideal_gas_eos_positive() {
        let s = CompressibleSph::new(1.4, 0.1);
        let p = s.ideal_gas_eos(1.2, 2.5e5);
        assert!(p > 0.0);
    }

    #[test]
    fn test_csph_ideal_gas_eos_zero_density() {
        let s = CompressibleSph::new(1.4, 0.1);
        assert_eq!(s.ideal_gas_eos(0.0, 2.5e5), 0.0);
    }

    #[test]
    fn test_csph_n_particles() {
        let mut s = CompressibleSph::new(1.4, 0.1);
        assert_eq!(s.n_particles(), 0);
        s.add_particle(CompressibleParticle::new(
            [0.0; 3], [0.0; 3], 1.0, 1e5, 2.5e5, 340.0, 0.01, 0.1,
        ));
        assert_eq!(s.n_particles(), 1);
    }

    #[test]
    fn test_csph_density_update_positive() {
        let mut s = CompressibleSph::new(1.4, 0.5);
        s.add_particle(CompressibleParticle::new(
            [0.0; 3], [0.0; 3], 1.0, 1e5, 2.5e5, 340.0, 0.1, 0.5,
        ));
        s.add_particle(CompressibleParticle::new(
            [0.1; 3], [0.0; 3], 1.0, 1e5, 2.5e5, 340.0, 0.1, 0.5,
        ));
        s.density_update();
        for p in &s.particles {
            assert!(p.density > 0.0);
        }
    }

    #[test]
    fn test_csph_step_does_not_panic() {
        let mut s = CompressibleSph::new(1.4, 0.5);
        s.add_particle(CompressibleParticle::new(
            [0.0; 3],
            [1.0, 0.0, 0.0],
            1.0,
            1e5,
            2.5e5,
            340.0,
            0.1,
            0.5,
        ));
        s.add_particle(CompressibleParticle::new(
            [0.2; 3],
            [0.0, 0.0, 0.0],
            1.0,
            1e5,
            2.5e5,
            340.0,
            0.1,
            0.5,
        ));
        s.step(1e-5);
    }

    #[test]
    fn test_csph_compute_sound_speed_after_init() {
        let mut s = CompressibleSph::new(1.4, 0.1);
        s.add_particle(CompressibleParticle::new(
            [0.0; 3], [0.0; 3], 1.2, 1e5, 2.07e5, 0.0, 0.01, 0.1,
        ));
        s.compute_sound_speed();
        assert!(s.particles[0].sound_speed > 0.0);
    }

    #[test]
    fn test_csph_total_energy_non_negative() {
        let mut s = CompressibleSph::new(1.4, 0.1);
        s.add_particle(CompressibleParticle::new(
            [0.0; 3],
            [1.0, 0.0, 0.0],
            1.0,
            1e5,
            2.5e5,
            340.0,
            0.01,
            0.1,
        ));
        assert!(s.total_energy() >= 0.0);
    }

    #[test]
    fn test_csph_dt_cfl_positive() {
        let mut s = CompressibleSph::new(1.4, 0.1);
        s.add_particle(CompressibleParticle::new(
            [0.0; 3], [0.0; 3], 1.2, 1e5, 2.07e5, 343.0, 0.01, 0.1,
        ));
        let dt = s.dt_cfl(0.3);
        assert!(dt > 0.0 && dt.is_finite());
    }

    // ---- ArtificialViscosityCs ----------------------------------------------

    #[test]
    fn test_art_visc_approaching_positive() {
        let av = ArtificialViscosityCs::new(1.0, 2.0);
        let pi = CompressibleParticle::new(
            [0.0, 0.0, 0.0],
            [100.0, 0.0, 0.0],
            1.0,
            1e5,
            2.5e5,
            340.0,
            0.01,
            0.1,
        );
        let pj = CompressibleParticle::new(
            [0.5, 0.0, 0.0],
            [-100.0, 0.0, 0.0],
            1.0,
            1e5,
            2.5e5,
            340.0,
            0.01,
            0.1,
        );
        let pi_ij = av.viscosity_term(&pi, &pj, 0.1);
        assert!(pi_ij > 0.0, "Approaching particles must give positive Π_ij");
    }

    #[test]
    fn test_art_visc_receding_zero() {
        let av = ArtificialViscosityCs::new(1.0, 2.0);
        let pi = CompressibleParticle::new(
            [0.0, 0.0, 0.0],
            [-100.0, 0.0, 0.0],
            1.0,
            1e5,
            2.5e5,
            340.0,
            0.01,
            0.1,
        );
        let pj = CompressibleParticle::new(
            [0.5, 0.0, 0.0],
            [100.0, 0.0, 0.0],
            1.0,
            1e5,
            2.5e5,
            340.0,
            0.01,
            0.1,
        );
        assert_eq!(av.viscosity_term(&pi, &pj, 0.1), 0.0);
    }

    #[test]
    fn test_signal_velocity_positive() {
        let av = ArtificialViscosityCs::new(1.0, 2.0);
        let pi = CompressibleParticle::new(
            [0.0, 0.0, 0.0],
            [100.0, 0.0, 0.0],
            1.0,
            1e5,
            2.5e5,
            340.0,
            0.01,
            0.1,
        );
        let pj = CompressibleParticle::new(
            [0.5, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            1.0,
            1e5,
            2.5e5,
            340.0,
            0.01,
            0.1,
        );
        let vs = av.signal_velocity(&pi, &pj);
        assert!(vs > 0.0);
    }

    // ---- GodunovSph ---------------------------------------------------------

    #[test]
    fn test_godunov_roe_symmetric_state() {
        let g = GodunovSph::new(RiemannSolverType::Roe, 1.4);
        let state = (1.0, 1e5, 0.0);
        let (p_star, u_star, _) = g.solve_riemann(state, state);
        assert!((p_star - 1e5).abs() / 1e5 < 0.01);
        assert!(u_star.abs() < 1.0);
    }

    #[test]
    fn test_godunov_hllc_pressure_positive() {
        let g = GodunovSph::new(RiemannSolverType::HllC, 1.4);
        let (p_star, _, _) = g.solve_riemann((1.0, 1e5, 0.0), (0.125, 0.1e5, 0.0));
        assert!(p_star > 0.0);
    }

    #[test]
    fn test_godunov_exact_sod_pressure_between() {
        let g = GodunovSph::new(RiemannSolverType::Exact, 1.4);
        let (p_star, _, _) = g.solve_riemann((1.0, 1.0, 0.0), (0.125, 0.1, 0.0));
        // Star pressure must be between the two states
        assert!(
            p_star > 0.1 && p_star < 1.0,
            "p_star={p_star} not in (0.1, 1.0)"
        );
    }

    #[test]
    fn test_godunov_flux_mass_conservation() {
        let g = GodunovSph::new(RiemannSolverType::HllC, 1.4);
        let (f_rho, _, _) = g.godunov_flux((1.0, 1e5, 10.0), (1.0, 1e5, 10.0));
        assert!(f_rho.is_finite());
    }

    // ---- StiffenedGasEos ----------------------------------------------------

    #[test]
    fn test_stiffened_eos_ideal_gas_limit() {
        let eos = StiffenedGasEos::new(1.4, 0.0);
        let rho = 1.2;
        let e = 2.0e5;
        let p = eos.pressure(rho, e);
        let expected = (1.4 - 1.0) * rho * e;
        assert!((p - expected).abs() < 1.0);
    }

    #[test]
    fn test_stiffened_eos_sound_speed_positive() {
        let eos = StiffenedGasEos::new(7.15, 3e8);
        let cs = eos.sound_speed(1000.0, 1e5);
        assert!(cs > 0.0);
    }

    #[test]
    fn test_stiffened_eos_sound_speed_zero_rho() {
        let eos = StiffenedGasEos::new(1.4, 0.0);
        assert_eq!(eos.sound_speed(0.0, 1e5), 0.0);
    }

    #[test]
    fn test_stiffened_eos_temperature_positive() {
        let eos = StiffenedGasEos::new(1.4, 0.0);
        let t = eos.temperature(1.0, 1e5);
        assert!(t > 0.0);
    }

    #[test]
    fn test_stiffened_eos_temperature_zero_rho() {
        let eos = StiffenedGasEos::new(1.4, 0.0);
        assert_eq!(eos.temperature(0.0, 1e5), 0.0);
    }

    // ---- rankine_hugoniot ---------------------------------------------------

    #[test]
    fn test_rh_mach_one_no_jump() {
        // M = 1: no shock, ρ₂ should equal ρ₁
        let (rho2, p2, _u2) = rankine_hugoniot(1.0, 1.0, 0.0, 1.0, 1.4);
        // For M=1 the jump is minimal
        assert!((rho2 - 1.0).abs() < 1e-10, "ρ₂={rho2} at M=1");
        assert!((p2 - 1.0).abs() < 1e-10, "p₂={p2}  at M=1");
    }

    #[test]
    fn test_rh_density_increases_across_shock() {
        let (rho2, _p2, _u2) = rankine_hugoniot(1.0, 1e5, 0.0, 2.0, 1.4);
        assert!(
            rho2 > 1.0,
            "Post-shock density must exceed pre-shock density"
        );
    }

    #[test]
    fn test_rh_pressure_increases_across_shock() {
        let (_rho2, p2, _u2) = rankine_hugoniot(1.0, 1e5, 0.0, 2.0, 1.4);
        assert!(
            p2 > 1e5,
            "Post-shock pressure must exceed pre-shock pressure"
        );
    }

    #[test]
    fn test_rh_density_ratio_strong_shock_limit() {
        // In the limit M → ∞: ρ₂/ρ₁ → (γ+1)/(γ-1) for γ=1.4 → 6
        let (rho2, _p2, _u2) = rankine_hugoniot(1.0, 1.0, 0.0, 1000.0, 1.4);
        let expected_ratio = (1.4 + 1.0) / (1.4 - 1.0);
        assert!(
            (rho2 / 1.0 - expected_ratio).abs() < 0.01 * expected_ratio,
            "ρ₂/ρ₁ = {:.4}, expected ≈ {:.4}",
            rho2,
            expected_ratio
        );
    }

    #[test]
    fn test_rh_u2_finite() {
        let (_rho2, _p2, u2) = rankine_hugoniot(1.2, 1e5, 500.0, 3.0, 1.4);
        assert!(u2.is_finite());
    }
}

#[cfg(test)]
mod tests_dynamics {
    use super::*;
    use std::f64::consts::PI;

    const GAMMA: f64 = 1.4;

    /// Gas particle with the given kinematics; density, pressure and sound
    /// speed are (re)computed by the solver.
    fn gas(pos: [f64; 3], vel: [f64; 3], mass: f64, energy: f64, h: f64) -> CompressibleParticle {
        CompressibleParticle::new(pos, vel, 1.0, 0.0, energy, 0.0, mass, h)
    }

    /// `n`³ lattice of unit-mass particles with spacing `dx`.
    fn lattice_3d(n: usize, dx: f64, h: f64, energy: f64) -> CompressibleSph {
        let mut sph = CompressibleSph::new(GAMMA, h);
        for i in 0..n {
            for j in 0..n {
                for k in 0..n {
                    let pos = [i as f64 * dx, j as f64 * dx, k as f64 * dx];
                    sph.add_particle(gas(pos, [0.0; 3], dx * dx * dx, energy, h));
                }
            }
        }
        sph.density_update();
        sph.compute_sound_speed();
        sph
    }

    /// Classic Sod shock tube on `x ∈ [−0.5, 0.5]` (ρ = 1, p = 1 left of the
    /// diaphragm, ρ = 0.125, p = 0.1 right of it).
    ///
    /// All particles have the same mass `ρ_L Δx_L`: `n_left` particles with
    /// spacing `Δx_L = 0.5 / n_left` fill the left half and `n_left / 8`
    /// particles with eight times the spacing fill the right half.  Equal
    /// masses keep the summation density consistent across the density jump
    /// (variable-mass particles on a uniform lattice would need a smoothing
    /// length of several spacings to avoid an O(1) error at the contact), and
    /// the global smoothing length `h = 1.3 Δx_R` covers the coarse right side.
    fn sod_tube(n_left: usize) -> CompressibleSph {
        let n_right = n_left / 8;
        let dx_left = 0.5 / n_left as f64;
        let dx_right = 0.5 / n_right as f64;
        let mass = dx_left;
        let h = 1.3 * dx_right;
        let mut sph = CompressibleSph::new(GAMMA, h)
            .with_dimension(KernelDimension::One)
            .with_viscosity(ArtificialViscosityCs::new(1.0, 2.0));
        let left = (0..n_left).map(|i| (-0.5 + (i as f64 + 0.5) * dx_left, 1.0, 1.0));
        let right = (0..n_right).map(|i| ((i as f64 + 0.5) * dx_right, 0.125, 0.1));
        for (x, rho, p) in left.chain(right) {
            let energy = p / ((GAMMA - 1.0) * rho);
            sph.add_particle(CompressibleParticle::new(
                [x, 0.0, 0.0],
                [0.0; 3],
                rho,
                p,
                energy,
                (GAMMA * p / rho).sqrt(),
                mass,
                h,
            ));
        }
        sph.density_update();
        sph.compute_sound_speed();
        sph
    }

    /// Advance with CFL-limited steps up to `t_end`; returns the step count.
    fn evolve(sph: &mut CompressibleSph, t_end: f64, cfl: f64) -> usize {
        let mut time = 0.0;
        let mut steps = 0;
        while time < t_end {
            let dt = sph.dt_cfl(cfl).min(t_end - time);
            sph.step(dt);
            time += dt;
            steps += 1;
        }
        steps
    }

    /// Midpoint-rule integral of `f(r) * area(r) dr` over the support `[0, 2h]`.
    fn radial_integral(h: f64, area: impl Fn(f64) -> f64, f: impl Fn(f64) -> f64) -> f64 {
        let n = 40_000;
        let dr = 2.0 * h / n as f64;
        (0..n)
            .map(|i| {
                let r = (i as f64 + 0.5) * dr;
                f(r) * area(r) * dr
            })
            .sum()
    }

    // ---- kernel normalisation -------------------------------------------------

    #[test]
    fn test_kernel_dimension_counts() {
        assert_eq!(KernelDimension::One.count(), 1);
        assert_eq!(KernelDimension::Two.count(), 2);
        assert_eq!(KernelDimension::Three.count(), 3);
        assert_eq!(KernelDimension::default(), KernelDimension::Three);
    }

    #[test]
    fn test_cubic_spline_kernel_is_unit_normalised_in_each_dimension() {
        for &h in &[0.1, 1.0, 2.7] {
            let one = radial_integral(
                h,
                |_| 2.0,
                |r| cubic_spline_kernel(r, h, KernelDimension::One),
            );
            let two = radial_integral(
                h,
                |r| 2.0 * PI * r,
                |r| cubic_spline_kernel(r, h, KernelDimension::Two),
            );
            let three = radial_integral(
                h,
                |r| 4.0 * PI * r * r,
                |r| cubic_spline_kernel(r, h, KernelDimension::Three),
            );
            assert!((one - 1.0).abs() < 1.0e-6, "1-D, h = {h}: {one}");
            assert!((two - 1.0).abs() < 1.0e-6, "2-D, h = {h}: {two}");
            assert!((three - 1.0).abs() < 1.0e-6, "3-D, h = {h}: {three}");
        }
    }

    #[test]
    fn test_cubic_spline_kernel_gradient_matches_central_difference() {
        let h = 0.8;
        for dim in [
            KernelDimension::One,
            KernelDimension::Two,
            KernelDimension::Three,
        ] {
            for i in 1..32 {
                let r = i as f64 * 0.05 * h;
                let eps = 1.0e-6 * h;
                let numeric = (cubic_spline_kernel(r + eps, h, dim)
                    - cubic_spline_kernel(r - eps, h, dim))
                    / (2.0 * eps);
                let analytic = cubic_spline_kernel_grad(r, h, dim);
                assert!(analytic <= 0.0);
                assert!(
                    (analytic - numeric).abs() <= 1.0e-6 * numeric.abs().max(1.0e-9),
                    "{dim:?}, r = {r}: analytic = {analytic}, numeric = {numeric}"
                );
            }
        }
    }

    #[test]
    fn test_summation_density_reproduces_rest_density_per_dimension() {
        // Equal-mass particles on a regular 1-D, 2-D and 3-D lattice must give
        // the lattice density in the centre when the matching kernel
        // normalisation is selected (the former fixed 1-D constant was wrong
        // by a factor ~ h^(1 - d) for planar and volumetric flows).
        let dx = 0.1;
        let h = 1.3 * dx;
        let rho0 = 2.0;
        let cases: [(KernelDimension, usize, f64); 3] = [
            (KernelDimension::One, 21, rho0 * dx),
            (KernelDimension::Two, 11, rho0 * dx * dx),
            (KernelDimension::Three, 7, rho0 * dx * dx * dx),
        ];
        for (dim, n, mass) in cases {
            let mut sph = CompressibleSph::new(GAMMA, h).with_dimension(dim);
            // Number of lattice points along x, y and z (1 along unused axes).
            let counts = [
                n,
                if dim.count() >= 2 { n } else { 1 },
                if dim.count() >= 3 { n } else { 1 },
            ];
            let coord =
                |index: usize, count: usize| (index as f64 - 0.5 * (count as f64 - 1.0)) * dx;
            for i in 0..counts[0] {
                for j in 0..counts[1] {
                    for k in 0..counts[2] {
                        let pos = [
                            coord(i, counts[0]),
                            coord(j, counts[1]),
                            coord(k, counts[2]),
                        ];
                        sph.add_particle(gas(pos, [0.0; 3], mass, 1.0, h));
                    }
                }
            }
            sph.density_update();
            let middle = sph
                .particles
                .iter()
                .min_by(|a, b| len3(a.pos).total_cmp(&len3(b.pos)))
                .map(|p| p.density)
                .unwrap_or(0.0);
            assert!(
                (middle - rho0).abs() / rho0 < 0.01,
                "{dim:?}: centre density {middle}, expected {rho0}"
            );
        }
    }

    // ---- time integration -----------------------------------------------------

    #[test]
    fn test_step_drifts_an_isolated_particle_by_velocity_times_dt() {
        let mut sph = CompressibleSph::new(GAMMA, 0.5);
        sph.add_particle(gas([0.0; 3], [2.0, -1.0, 0.5], 0.1, 1.0, 0.5));
        sph.step(0.1);
        let p = &sph.particles[0];
        for (got, want) in p.pos.iter().zip([0.2, -0.1, 0.05]) {
            assert!((got - want).abs() < 1.0e-14, "pos {:?}", p.pos);
        }
        // No neighbours: no force and no energy exchange.
        assert_eq!(p.vel, [2.0, -1.0, 0.5]);
        assert!((p.energy - 1.0).abs() < 1.0e-14);
    }

    #[test]
    fn test_step_ignores_non_positive_or_non_finite_dt() {
        let mut sph = lattice_3d(2, 1.0, 1.3, 1.0);
        let before: Vec<([f64; 3], [f64; 3], f64)> = sph
            .particles
            .iter()
            .map(|p| (p.pos, p.vel, p.energy))
            .collect();
        for dt in [0.0, -0.1, f64::NAN, f64::INFINITY] {
            sph.step(dt);
        }
        for (p, (pos, vel, energy)) in sph.particles.iter().zip(before) {
            assert_eq!(p.pos, pos);
            assert_eq!(p.vel, vel);
            assert_eq!(p.energy, energy);
        }
    }

    #[test]
    fn test_step_velocity_change_scales_linearly_with_dt() {
        // Two hot particles at rest push each other apart.  For a tiny step
        // the velocity change is a·dt, so doubling dt doubles it — the old
        // code ignored dt altogether.
        let kick = |dt: f64| {
            let mut sph = CompressibleSph::new(GAMMA, 1.0);
            sph.add_particle(gas([0.0; 3], [0.0; 3], 1.0, 2.5, 1.0));
            sph.add_particle(gas([1.0, 0.0, 0.0], [0.0; 3], 1.0, 2.5, 1.0));
            sph.step(dt);
            (sph.particles[0].vel, sph.particles[1].vel)
        };
        let (v0_small, v1_small) = kick(1.0e-4);
        let (v0_large, _) = kick(2.0e-4);
        assert!(
            v0_small[0] < 0.0,
            "left particle must move to -x: {v0_small:?}"
        );
        assert!(
            v1_small[0] > 0.0,
            "right particle must move to +x: {v1_small:?}"
        );
        let ratio = v0_large[0] / v0_small[0];
        assert!((ratio - 2.0).abs() < 1.0e-3, "velocity ratio {ratio}");
        // Equal masses: exact momentum balance.
        assert!((v0_small[0] + v1_small[0]).abs() <= 1.0e-12 * v0_small[0].abs());
    }

    #[test]
    fn test_step_moves_particles_and_changes_kinematics() {
        let mut sph = lattice_3d(3, 1.0, 1.3, 1.0);
        let start: Vec<[f64; 3]> = sph.particles.iter().map(|p| p.pos).collect();
        sph.step(0.1);
        let moved = sph
            .particles
            .iter()
            .zip(start.iter())
            .map(|(p, x0)| len3(sub3(p.pos, *x0)))
            .fold(0.0_f64, f64::max);
        assert!(moved > 1.0e-4, "particles must move: {moved}");
        let vmax = sph
            .particles
            .iter()
            .map(|p| len3(p.vel))
            .fold(0.0_f64, f64::max);
        assert!(vmax > 0.0);
    }

    #[test]
    fn test_momentum_and_energy_updates_apply_dt_exactly() {
        let build = || {
            let mut sph = CompressibleSph::new(GAMMA, 1.0);
            sph.add_particle(gas([0.0; 3], [0.5, 0.0, 0.0], 1.0, 2.5, 1.0));
            sph.add_particle(gas([0.9, 0.2, 0.0], [-0.5, 0.1, 0.0], 2.0, 1.5, 1.0));
            sph.density_update();
            sph.compute_sound_speed();
            sph
        };
        // Velocity change over dt: linear in dt for a fixed state.
        let kick = |dt: f64| {
            let mut sph = build();
            let before = sph.particles[0].vel;
            sph.momentum_update(dt);
            let after = sph.particles[0].vel;
            [
                after[0] - before[0],
                after[1] - before[1],
                after[2] - before[2],
            ]
        };
        let dv_one = kick(1.0e-3);
        let dv_two = kick(2.0e-3);
        assert!(len3(dv_one) > 0.0, "momentum_update must act");
        for d in 0..3 {
            assert!((dv_two[d] - 2.0 * dv_one[d]).abs() <= 1.0e-12 * len3(dv_one));
        }
        // Specific-energy change over dt: also linear in dt.
        let heat = |dt: f64| {
            let mut sph = build();
            let before = sph.particles[1].energy;
            sph.energy_update(dt);
            sph.particles[1].energy - before
        };
        let de_one = heat(1.0e-3);
        let de_two = heat(2.0e-3);
        assert!(de_one != 0.0, "energy_update must act");
        assert!((de_two - 2.0 * de_one).abs() <= 1.0e-12 * de_one.abs());
        // Zero time step leaves the state untouched.
        let mut idle = build();
        let snapshot: Vec<([f64; 3], f64)> =
            idle.particles.iter().map(|p| (p.vel, p.energy)).collect();
        idle.momentum_update(0.0);
        idle.energy_update(0.0);
        for (p, (vel, energy)) in idle.particles.iter().zip(snapshot) {
            assert_eq!(p.vel, vel);
            assert_eq!(p.energy, energy);
        }
    }

    #[test]
    fn test_total_momentum_sums_particle_momenta() {
        let mut sph = CompressibleSph::new(GAMMA, 0.5);
        sph.add_particle(gas([0.0; 3], [1.0, 2.0, 3.0], 2.0, 1.0, 0.5));
        sph.add_particle(gas([5.0; 3], [-1.0, 0.0, 1.0], 3.0, 1.0, 0.5));
        let p = sph.total_momentum();
        for (got, want) in p.iter().zip([-1.0, 4.0, 9.0]) {
            assert!((got - want).abs() < 1.0e-14);
        }
    }

    #[test]
    fn test_pair_forces_conserve_momentum_and_energy_for_a_3d_cluster() {
        // Hot 4³ cluster with a velocity field, artificial viscosity on: the
        // pair interactions are antisymmetric, so momentum is conserved to
        // rounding and the leapfrog keeps the total energy to second order.
        let mut sph =
            lattice_3d(4, 1.0, 1.3, 1.0).with_viscosity(ArtificialViscosityCs::new(1.0, 2.0));
        for p in sph.particles.iter_mut() {
            p.vel = [
                0.2 - 0.1 * p.pos[0] + 0.05 * p.pos[1],
                -0.1 + 0.08 * p.pos[2],
                0.15 * (p.pos[0] - 1.5),
            ];
        }
        sph.density_update();
        sph.compute_sound_speed();
        let e0 = sph.total_energy();
        let mom0 = sph.total_momentum();
        assert!(len3(mom0) > 1.0);
        evolve(&mut sph, 2.0, 0.1);
        let e1 = sph.total_energy();
        let mom1 = sph.total_momentum();
        assert!(
            (e1 - e0).abs() / e0 < 2.0e-3,
            "energy drift {}",
            (e1 - e0) / e0
        );
        assert!(
            len3(sub3(mom1, mom0)) < 1.0e-12 * len3(mom0),
            "momentum {mom0:?} -> {mom1:?}"
        );
    }

    #[test]
    fn test_hot_cluster_expands_and_converts_heat_into_motion() {
        let mut sph = lattice_3d(4, 1.0, 1.3, 2.0);
        let radius = |s: &CompressibleSph| {
            let n = s.particles.len() as f64;
            let c = s.particles.iter().fold([0.0; 3], |a, p| add3(a, p.pos));
            let c = scale3(c, 1.0 / n);
            s.particles
                .iter()
                .map(|p| len3(sub3(p.pos, c)))
                .fold(0.0_f64, f64::max)
        };
        let r0 = radius(&sph);
        let u0: f64 = sph.particles.iter().map(|p| p.mass * p.energy).sum();
        evolve(&mut sph, 3.0, 0.1);
        assert!(
            radius(&sph) > 1.3 * r0,
            "cluster must expand: {} vs {r0}",
            radius(&sph)
        );
        let u1: f64 = sph.particles.iter().map(|p| p.mass * p.energy).sum();
        let ke: f64 = sph
            .particles
            .iter()
            .map(|p| 0.5 * p.mass * dot3(p.vel, p.vel))
            .sum();
        assert!(u1 < u0, "adiabatic expansion must cool the gas");
        assert!(ke > 0.1 * u0, "kinetic energy {ke}");
    }

    #[test]
    fn test_artificial_viscosity_decelerates_converging_flow_and_heats_it() {
        // Two lattice blocks approach each other along x.
        let build = |viscosity: Option<ArtificialViscosityCs>| {
            let mut sph = CompressibleSph::new(GAMMA, 1.3).with_dimension(KernelDimension::One);
            if let Some(av) = viscosity {
                sph = sph.with_viscosity(av);
            }
            for i in 0..20 {
                let x = i as f64;
                let v = if i < 10 { 1.0 } else { -1.0 };
                sph.add_particle(gas([x, 0.0, 0.0], [v, 0.0, 0.0], 1.0, 0.5, 1.3));
            }
            sph.density_update();
            sph.compute_sound_speed();
            sph
        };
        let mut plain = build(None);
        let mut damped = build(Some(ArtificialViscosityCs::new(1.0, 2.0)));
        let e0 = plain.total_energy();
        // Strong shocks need a small CFL number: dt_cfl ignores the viscous
        // signal speed of the supersonic collision.
        evolve(&mut plain, 3.0, 0.1);
        evolve(&mut damped, 3.0, 0.1);
        let kinetic = |s: &CompressibleSph| -> f64 {
            s.particles
                .iter()
                .map(|p| 0.5 * p.mass * dot3(p.vel, p.vel))
                .sum()
        };
        let thermal =
            |s: &CompressibleSph| -> f64 { s.particles.iter().map(|p| p.mass * p.energy).sum() };
        assert!(
            kinetic(&damped) < kinetic(&plain),
            "viscosity must dissipate kinetic energy"
        );
        assert!(
            thermal(&damped) > thermal(&plain),
            "dissipated energy must heat the gas"
        );
        assert!((damped.total_energy() - e0).abs() / e0 < 0.01);
        assert!((plain.total_energy() - e0).abs() / e0 < 0.01);
    }

    // ---- Sod shock tube -------------------------------------------------------

    #[test]
    fn test_sod_shock_tube_moves_particles_and_conserves_energy_and_momentum() {
        let mut sph = sod_tube(200);
        let x0: Vec<f64> = sph.particles.iter().map(|p| p.pos[0]).collect();
        let e0 = sph.total_energy();
        let steps = evolve(&mut sph, 0.15, 0.25);
        assert!(steps > 10);
        // Particles move: the gas around the diaphragm is accelerated to the
        // contact speed u* ≈ 0.927.
        let displacement = sph
            .particles
            .iter()
            .zip(x0.iter())
            .map(|(p, x)| (p.pos[0] - x).abs())
            .fold(0.0_f64, f64::max);
        assert!(displacement > 0.05, "max displacement {displacement}");
        // Total energy within 1 %, momentum conserved to rounding.
        let e1 = sph.total_energy();
        assert!(
            (e1 - e0).abs() / e0 < 0.01,
            "energy drift {}",
            (e1 - e0) / e0
        );
        let mom = sph.total_momentum();
        let scale = sph
            .particles
            .iter()
            .map(|p| p.mass * p.vel[0].abs())
            .sum::<f64>();
        assert!(len3(mom) < 1.0e-12 * scale, "momentum {mom:?}");
        // Motion stays on the tube axis.
        for p in &sph.particles {
            assert_eq!(p.pos[1], 0.0);
            assert_eq!(p.pos[2], 0.0);
        }
    }

    #[test]
    fn test_sod_shock_tube_plateaus_match_exact_riemann_solution() {
        // Exact solution (γ = 1.4): p* = 0.30313, u* = 0.92745, ρ*_L = 0.42632
        // between the rarefaction tail (x = −0.0105 at t = 0.15) and the
        // contact (x = u* t = 0.139), ρ*_R = 0.26557 between contact and shock
        // (x = 1.7522 t = 0.263).  The windows below lie inside these regions.
        let mut sph = sod_tube(200);
        evolve(&mut sph, 0.15, 0.25);
        let mean = |lo: f64, hi: f64, f: &dyn Fn(&CompressibleParticle) -> f64| -> f64 {
            let vals: Vec<f64> = sph
                .particles
                .iter()
                .filter(|p| (lo..=hi).contains(&p.pos[0]))
                .map(f)
                .collect();
            assert!(!vals.is_empty(), "no particles in [{lo}, {hi}]");
            vals.iter().sum::<f64>() / vals.len() as f64
        };
        // Star region behind the rarefaction.
        let u_left = mean(0.03, 0.12, &|p| p.vel[0]);
        let p_left = mean(0.03, 0.12, &|p| p.pressure);
        let rho_left = mean(0.03, 0.12, &|p| p.density);
        assert!((u_left - 0.92745).abs() < 0.04, "u* = {u_left}");
        assert!((p_left - 0.30313).abs() < 0.02, "p* = {p_left}");
        assert!((rho_left - 0.42632).abs() < 0.03, "rho*_L = {rho_left}");
        // Shocked gas between contact and shock.
        let u_right = mean(0.17, 0.24, &|p| p.vel[0]);
        let p_right = mean(0.17, 0.24, &|p| p.pressure);
        let rho_right = mean(0.17, 0.24, &|p| p.density);
        assert!(
            (u_right - 0.92745).abs() < 0.04,
            "u after shock = {u_right}"
        );
        assert!(
            (p_right - 0.30313).abs() < 0.02,
            "p after shock = {p_right}"
        );
        assert!((rho_right - 0.26557).abs() < 0.02, "rho*_R = {rho_right}");
        // The undisturbed left state ahead of the rarefaction head.
        let rho_far_left = mean(-0.30, -0.20, &|p| p.density);
        let u_far_left = mean(-0.30, -0.20, &|p| p.vel[0]);
        assert!(
            (rho_far_left - 1.0).abs() < 0.03,
            "left state {rho_far_left}"
        );
        assert!(u_far_left.abs() < 0.01, "left state velocity {u_far_left}");
    }
}
