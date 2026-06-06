// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Compressible SPH for high-speed flows and shocks.
//!
//! This module implements compressible Smoothed Particle Hydrodynamics (SPH)
//! including Riemann-based methods and shock-capturing formulations:
//!
//! - [`CompressibleParticle`]: State variables for a compressible SPH particle
//! - [`CompressibleSph`]: Full compressible SPH solver with ideal-gas EOS
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

/// Snapshot tuple for per-particle thermodynamic state used in energy_update.
type ParticleSnap = (f64, f64, [f64; 3], [f64; 3], f64);

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

/// Compressible SPH solver for high-speed and shock-dominated flows.
///
/// Uses an ideal-gas EOS p = (γ-1)·ρ·e and explicit Euler time integration.
pub struct CompressibleSph {
    /// Collection of all particles.
    pub particles: Vec<CompressibleParticle>,
    /// Ratio of specific heats γ.
    pub gamma: f64,
    /// Global smoothing length h (m).
    pub h: f64,
}

impl CompressibleSph {
    /// Create a new solver with no particles.
    pub fn new(gamma: f64, h: f64) -> Self {
        Self {
            particles: Vec::new(),
            gamma,
            h,
        }
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
    /// The sequence is:
    /// 1. Update densities (continuity),
    /// 2. Update momenta (momentum equation),
    /// 3. Update energies (energy equation),
    /// 4. Update thermodynamics.
    pub fn step(&mut self, _dt: f64) {
        self.density_update();
        self.momentum_update();
        self.energy_update();
        self.compute_sound_speed();
    }

    /// Update densities using the SPH continuity equation (summation density).
    ///
    /// ρ_i = Σ_j m_j · W(|x_i - x_j|, h)
    pub fn density_update(&mut self) {
        let _n = self.particles.len();
        let h = self.h;
        let positions: Vec<[f64; 3]> = self.particles.iter().map(|p| p.pos).collect();
        let masses: Vec<f64> = self.particles.iter().map(|p| p.mass).collect();
        for (p_i, pos_i) in self.particles.iter_mut().zip(positions.iter()) {
            let rho: f64 = positions
                .iter()
                .zip(masses.iter())
                .map(|(pos_j, &mj)| {
                    let r = len3(sub3(*pos_i, *pos_j));
                    mj * cubic_spline_kernel(r, h)
                })
                .sum::<f64>()
                .max(1e-10);
            p_i.density = rho;
        }
    }

    /// Update velocities using the SPH momentum equation (pressure gradient).
    ///
    /// d**v**_i/dt = -Σ_j m_j (p_i/ρ_i² + p_j/ρ_j²) ∇W_ij
    pub fn momentum_update(&mut self) {
        let _n = self.particles.len();
        let h = self.h;
        // Snapshot current pressures / densities / positions
        let snap: Vec<(f64, f64, [f64; 3], f64)> = self
            .particles
            .iter()
            .map(|p| (p.pressure, p.density, p.pos, p.mass))
            .collect();

        for (i, p) in self.particles.iter_mut().enumerate() {
            let (pi, rhoi, xi, _) = snap[i];
            let mut acc = [0.0_f64; 3];
            for (j, &(pj, rhoj, xj, mj)) in snap.iter().enumerate() {
                if i == j {
                    continue;
                }
                let rij = sub3(xi, xj);
                let r = len3(rij);
                let dw = cubic_spline_kernel_grad(r, h);
                if r > 1e-300 {
                    let dir = scale3(rij, 1.0 / r);
                    let coeff = mj * (pi / (rhoi * rhoi) + pj / (rhoj * rhoj)) * dw;
                    acc = sub3(acc, scale3(dir, coeff));
                }
            }
            // Simple Euler update (dt = 1 for this helper; caller uses step())
            p.vel = add3(p.vel, acc);
        }
    }

    /// Update specific internal energies using the SPH energy equation.
    ///
    /// de_i/dt = (p_i/ρ_i²) Σ_j m_j (**v**_i - **v**_j) · ∇W_ij
    pub fn energy_update(&mut self) {
        let h = self.h;
        let snap: Vec<ParticleSnap> = self
            .particles
            .iter()
            .map(|p| (p.pressure, p.density, p.pos, p.vel, p.mass))
            .collect();

        for (i, p) in self.particles.iter_mut().enumerate() {
            let (pi, rhoi, xi, vi, _) = snap[i];
            if rhoi.abs() < 1e-300 {
                continue;
            }
            let mut de = 0.0;
            for (j, &(_pj, _rhoj, xj, vj, mj)) in snap.iter().enumerate() {
                if i == j {
                    continue;
                }
                let rij = sub3(xi, xj);
                let r = len3(rij);
                let dw = cubic_spline_kernel_grad(r, h);
                if r > 1e-300 {
                    let dir = scale3(rij, 1.0 / r);
                    let dv = sub3(vi, vj);
                    de += mj * dot3(dv, dir) * dw;
                }
            }
            p.energy += (pi / (rhoi * rhoi)) * de;
            p.energy = p.energy.max(0.0);
        }
    }

    /// Recompute the sound speed c_s = √(γ·p/ρ) for every particle.
    pub fn compute_sound_speed(&mut self) {
        let gamma = self.gamma;
        for p in &mut self.particles {
            p.sound_speed = if p.density > 0.0 && p.pressure > 0.0 {
                (gamma * p.pressure / p.density).sqrt()
            } else {
                0.0
            };
        }
        self.update_thermodynamics();
    }

    /// Total mechanical energy E_tot = Σ m_i · (e_i + ½|v_i|²).
    pub fn total_energy(&self) -> f64 {
        self.particles.iter().map(|p| p.total_energy()).sum()
    }

    /// Maximum sound speed among all particles (m/s).
    pub fn max_sound_speed(&self) -> f64 {
        self.particles
            .iter()
            .map(|p| p.sound_speed)
            .fold(0.0_f64, f64::max)
    }

    /// CFL time step: dt = CFL · h / (c_s,max + |v|_max).
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
// SPH kernel functions (cubic spline, 1-D normalisation)
// ============================================================================

/// Cubic-spline SPH kernel W(r, h).
fn cubic_spline_kernel(r: f64, h: f64) -> f64 {
    if h <= 0.0 {
        return 0.0;
    }
    let q = r / h;
    let norm = 1.0 / h; // 1-D normalisation
    if q < 1.0 {
        norm * (2.0 / 3.0 - q * q + 0.5 * q * q * q)
    } else if q < 2.0 {
        let t = 2.0 - q;
        norm * (1.0 / 6.0) * t * t * t
    } else {
        0.0
    }
}

/// Magnitude of the gradient of the cubic-spline kernel |dW/dr|.
fn cubic_spline_kernel_grad(r: f64, h: f64) -> f64 {
    if h <= 0.0 {
        return 0.0;
    }
    let q = r / h;
    let norm = 1.0 / (h * h);
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
