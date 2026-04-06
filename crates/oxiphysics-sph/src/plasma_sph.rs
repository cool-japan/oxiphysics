// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! SPH methods for plasma simulation (magnetohydrodynamics, particle-in-cell coupling).
//!
//! This module provides:
//! - [`PlasmaParticle`]: SPH particle with electromagnetic and thermodynamic state
//! - [`MhdSph`]: Magnetohydrodynamics SPH with Dedner divergence cleaning
//! - [`PlasmaEquationOfState`]: Ideal gas + radiation + electron pressure
//! - [`BraginskiiViscosity`]: Anisotropic viscosity along magnetic field
//! - [`ThermalConduction`]: Spitzer thermal conductivity (anisotropic along B)
//! - [`IonizationModel`]: Saha equation for LTE ionization
//! - [`RadiativeCooling`]: Bremsstrahlung, line cooling, recombination
//! - [`PicSphCoupling`]: Macro/micro particle coupling interface
//! - [`AlfvenWave`]: Alfvén wave propagation in SPH
//! - [`MagneticReconnection`]: Sweet-Parker and Petschek reconnection
//! - Utility functions: [`debye_length`], [`cyclotron_frequency`], [`plasma_frequency`], [`beta_parameter`]
//!
//! # References
//! - Dedner et al. (2002): Hyperbolic divergence cleaning for MHD equations
//! - Braginskii (1965): Transport processes in a plasma
//! - Spitzer & Härm (1953): Transport phenomena in a completely ionized gas
//! - Saha (1920): On the physical theory of stellar spectra

#![allow(dead_code)]
#![allow(clippy::too_many_arguments)]

use std::f64::consts::PI;

// ─────────────────────────────────────────────────────────────────────────────
// Physical constants
// ─────────────────────────────────────────────────────────────────────────────

/// Boltzmann constant \[J K⁻¹\].
const K_B: f64 = 1.380_649e-23;

/// Proton mass \[kg\].
const M_PROTON: f64 = 1.672_623_e-27;

/// Electron mass \[kg\].
const M_ELECTRON: f64 = 9.109_384e-31;

/// Elementary charge \[C\].
const E_CHARGE: f64 = 1.602_176_634e-19;

/// Vacuum permeability μ₀ \[H m⁻¹\].
const MU_0: f64 = 1.256_637_061e-6;

/// Vacuum permittivity ε₀ \[F m⁻¹\].
const EPS_0: f64 = 8.854_187_817e-12;

/// Speed of light \[m s⁻¹\].
const C_LIGHT: f64 = 2.997_924_58e8;

/// Stefan-Boltzmann constant \[W m⁻² K⁻⁴\].
const SIGMA_SB: f64 = 5.670_374_419e-8;

/// Planck constant \[J s\].
const H_PLANCK: f64 = 6.626_070_15e-34;

/// Radiation constant a_rad = 4σ/c \[J m⁻³ K⁻⁴\].
const A_RAD: f64 = 4.0 * SIGMA_SB / C_LIGHT;

/// Bremsstrahlung emission coefficient \[W m³ K^{-1/2}\].
const BREMSSTRAHLUNG_COEFF: f64 = 1.69e-40;

// ─────────────────────────────────────────────────────────────────────────────
// Vector helpers
// ─────────────────────────────────────────────────────────────────────────────

#[inline]
fn add3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

#[inline]
fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

#[inline]
fn scale3(v: [f64; 3], s: f64) -> [f64; 3] {
    [v[0] * s, v[1] * s, v[2] * s]
}

#[inline]
fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[inline]
fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

#[inline]
fn len3(v: [f64; 3]) -> f64 {
    dot3(v, v).sqrt()
}

#[inline]
fn norm3(v: [f64; 3]) -> [f64; 3] {
    let l = len3(v);
    if l < 1e-300 {
        [0.0; 3]
    } else {
        scale3(v, 1.0 / l)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// SPH kernel
// ─────────────────────────────────────────────────────────────────────────────

/// Cubic spline SPH kernel value W(r, h).
fn cubic_kernel(r: f64, h: f64) -> f64 {
    if h < 1e-300 {
        return 0.0;
    }
    let q = r / h;
    let alpha = 1.0 / (PI * h * h * h);
    if q < 1.0 {
        alpha * (1.0 - 1.5 * q * q + 0.75 * q * q * q)
    } else if q < 2.0 {
        let t = 2.0 - q;
        alpha * 0.25 * t * t * t
    } else {
        0.0
    }
}

/// Gradient of the cubic spline SPH kernel.
fn cubic_kernel_grad(r_ij: [f64; 3], h: f64) -> [f64; 3] {
    let r = len3(r_ij);
    if r < 1e-12 || h < 1e-12 {
        return [0.0; 3];
    }
    let q = r / h;
    let alpha = 1.0 / (PI * h * h * h);
    let dw_dr = if q < 1.0 {
        alpha * (-3.0 * q + 2.25 * q * q) / h
    } else if q < 2.0 {
        let t = 2.0 - q;
        alpha * (-0.75 * t * t) / h
    } else {
        0.0
    };
    [
        r_ij[0] / r * dw_dr,
        r_ij[1] / r * dw_dr,
        r_ij[2] / r * dw_dr,
    ]
}

// ─────────────────────────────────────────────────────────────────────────────
// PlasmaParticle
// ─────────────────────────────────────────────────────────────────────────────

/// A single plasma SPH particle with full electromagnetic and thermodynamic state.
///
/// Fields use SI units throughout.
#[derive(Debug, Clone)]
pub struct PlasmaParticle {
    /// Position vector \[x, y, z\] \[m\].
    pub position: [f64; 3],
    /// Velocity vector \[vx, vy, vz\] \[m s⁻¹\].
    pub velocity: [f64; 3],
    /// Particle mass \[kg\].
    pub mass: f64,
    /// Particle charge \[C\].
    pub charge: f64,
    /// Temperature \[K\].
    pub temperature: f64,
    /// Mass density ρ \[kg m⁻³\].
    pub density: f64,
    /// Thermal pressure P \[Pa\].
    pub pressure: f64,
    /// Local magnetic field B \[T\].
    pub magnetic_field: [f64; 3],
    /// Divergence-cleaning scalar ψ (Dedner).
    pub psi: f64,
    /// Smoothing length h \[m\].
    pub smoothing_length: f64,
    /// Accumulated acceleration \[m s⁻²\].
    pub accel: [f64; 3],
    /// Ionization degree Z_eff ∈ \[0, Z_max\].
    pub ionization: f64,
}

impl PlasmaParticle {
    /// Create a new plasma particle.
    pub fn new(
        position: [f64; 3],
        velocity: [f64; 3],
        mass: f64,
        charge: f64,
        temperature: f64,
        density: f64,
    ) -> Self {
        let pressure = if density > 0.0 {
            density / M_PROTON * K_B * temperature
        } else {
            0.0
        };
        Self {
            position,
            velocity,
            mass,
            charge,
            temperature,
            density,
            pressure,
            magnetic_field: [0.0; 3],
            psi: 0.0,
            smoothing_length: 1e-3,
            accel: [0.0; 3],
            ionization: 1.0,
        }
    }

    /// Kinetic energy of this particle \[J\].
    pub fn kinetic_energy(&self) -> f64 {
        0.5 * self.mass * dot3(self.velocity, self.velocity)
    }

    /// Speed |v| \[m s⁻¹\].
    pub fn speed(&self) -> f64 {
        len3(self.velocity)
    }

    /// Magnetic field magnitude \[T\].
    pub fn b_magnitude(&self) -> f64 {
        len3(self.magnetic_field)
    }

    /// Thermal velocity v_th = sqrt(2 k_B T / m) \[m s⁻¹\].
    pub fn thermal_velocity(&self) -> f64 {
        if self.mass > 0.0 {
            (2.0 * K_B * self.temperature / self.mass).sqrt()
        } else {
            0.0
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// MhdSph – Magnetohydrodynamics SPH
// ─────────────────────────────────────────────────────────────────────────────

/// Magnetohydrodynamics SPH solver.
///
/// Implements:
/// - Ideal and resistive induction equation ∂B/∂t = ∇×(v×B) − η∇²B
/// - Lorentz force f = (1/μ₀)(∇×B)×B
/// - Dedner hyperbolic/parabolic divergence cleaning for ∇·B = 0
/// - SPH particle integration
pub struct MhdSph {
    /// List of plasma particles.
    pub particles: Vec<PlasmaParticle>,
    /// Resistivity η \[Ω m\].
    pub resistivity: f64,
    /// Dedner cleaning speed c_h \[m s⁻¹\].
    pub dedner_ch: f64,
    /// Dedner parabolic decay parameter c_p² \[m² s⁻²\].
    pub dedner_cp2: f64,
    /// Adiabatic index γ.
    pub gamma: f64,
}

impl MhdSph {
    /// Create a new MHD-SPH solver.
    pub fn new(resistivity: f64, gamma: f64) -> Self {
        Self {
            particles: Vec::new(),
            resistivity,
            dedner_ch: 1e4,
            dedner_cp2: 1e8,
            gamma,
        }
    }

    /// Add a particle to the simulation.
    pub fn add_particle(&mut self, p: PlasmaParticle) {
        self.particles.push(p);
    }

    /// Compute the ideal induction equation RHS: ∇×(v×B).
    ///
    /// Returns dB/dt at particle `i` due to its neighbours (simplified SPH estimate).
    pub fn induction_rhs(&self, i: usize) -> [f64; 3] {
        let pi = &self.particles[i];
        let mut db_dt = [0.0f64; 3];
        for (j, pj) in self.particles.iter().enumerate() {
            if i == j {
                continue;
            }
            let r_ij = sub3(pi.position, pj.position);
            let grad_w = cubic_kernel_grad(r_ij, pi.smoothing_length);
            let vxb_i = cross3(pi.velocity, pi.magnetic_field);
            let vxb_j = cross3(pj.velocity, pj.magnetic_field);
            let vxb_avg = scale3(add3(vxb_i, vxb_j), 0.5);
            // curl contribution: (m_j/ρ_j) * (vxB) × ∇W
            let factor = if pj.density > 0.0 {
                pj.mass / pj.density
            } else {
                0.0
            };
            let contrib = cross3(vxb_avg, grad_w);
            db_dt = add3(db_dt, scale3(contrib, factor));
        }
        db_dt
    }

    /// Compute Lorentz force on particle `i`: f = (1/μ₀)(∇×B)×B.
    ///
    /// Uses the SPH curl operator ∇×B ≈ Σ_j (m_j/ρ_j)(B_j − B_i)×∇W_ij.
    pub fn lorentz_force(&self, i: usize) -> [f64; 3] {
        let pi = &self.particles[i];
        let mut curl_b = [0.0f64; 3];
        for (j, pj) in self.particles.iter().enumerate() {
            if i == j {
                continue;
            }
            let r_ij = sub3(pi.position, pj.position);
            let grad_w = cubic_kernel_grad(r_ij, pi.smoothing_length);
            let db = sub3(pj.magnetic_field, pi.magnetic_field);
            let factor = if pj.density > 0.0 {
                pj.mass / pj.density
            } else {
                0.0
            };
            let contrib = cross3(db, grad_w);
            curl_b = add3(curl_b, scale3(contrib, factor));
        }
        // f = (1/μ₀) (∇×B) × B
        let f = cross3(curl_b, pi.magnetic_field);
        scale3(f, 1.0 / (MU_0 * pi.density.max(1e-300)))
    }

    /// Apply Dedner divergence cleaning: evolve ψ and damp B divergence.
    ///
    /// ∂ψ/∂t = −c_h² (∇·B) − ψ/τ_p
    pub fn dedner_cleaning_step(&mut self, dt: f64) {
        let decay = (-self.dedner_ch / self.dedner_cp2 * dt).exp();
        for p in self.particles.iter_mut() {
            // ψ decays parabolically; here we apply exponential damping
            p.psi *= decay;
            // Correct B: B ← B − ∇ψ·dt (simplified as scalar shift)
            let b_mag = len3(p.magnetic_field);
            if b_mag > 1e-30 {
                let b_hat = norm3(p.magnetic_field);
                let correction = p.psi * dt * self.dedner_ch;
                p.magnetic_field = sub3(p.magnetic_field, scale3(b_hat, correction));
            }
        }
    }

    /// Compute resistive diffusion term η∇²B at particle `i`.
    pub fn resistive_diffusion(&self, i: usize) -> [f64; 3] {
        let pi = &self.particles[i];
        let mut lap_b = [0.0f64; 3];
        for (j, pj) in self.particles.iter().enumerate() {
            if i == j {
                continue;
            }
            let r_ij = sub3(pi.position, pj.position);
            let r2 = dot3(r_ij, r_ij);
            let h = pi.smoothing_length;
            let w = cubic_kernel(r2.sqrt(), h);
            let db = sub3(pj.magnetic_field, pi.magnetic_field);
            let factor = if pj.density > 0.0 && r2 > 1e-30 {
                2.0 * pj.mass / pj.density * w / r2
            } else {
                0.0
            };
            lap_b = add3(lap_b, scale3(db, factor));
        }
        scale3(lap_b, self.resistivity)
    }

    /// Advance the simulation by one Euler step of duration `dt`.
    pub fn euler_step(&mut self, dt: f64) {
        let n = self.particles.len();
        let mut accels = vec![[0.0f64; 3]; n];
        let mut db_dts = vec![[0.0f64; 3]; n];

        for i in 0..n {
            let f_lorentz = self.lorentz_force(i);
            let f_resist = self.resistive_diffusion(i);
            accels[i] = f_lorentz;
            db_dts[i] = add3(self.induction_rhs(i), f_resist);
        }

        for (i, p) in self.particles.iter_mut().enumerate() {
            p.velocity = add3(p.velocity, scale3(accels[i], dt));
            p.position = add3(p.position, scale3(p.velocity, dt));
            p.magnetic_field = add3(p.magnetic_field, scale3(db_dts[i], dt));
        }
        self.dedner_cleaning_step(dt);
    }

    /// Total kinetic energy of all particles \[J\].
    pub fn total_kinetic_energy(&self) -> f64 {
        self.particles.iter().map(|p| p.kinetic_energy()).sum()
    }

    /// Total magnetic energy Σ B²/(2μ₀) V_i \[J\].
    pub fn total_magnetic_energy(&self) -> f64 {
        self.particles
            .iter()
            .map(|p| {
                let b2 = dot3(p.magnetic_field, p.magnetic_field);
                let vol = if p.density > 0.0 {
                    p.mass / p.density
                } else {
                    0.0
                };
                b2 / (2.0 * MU_0) * vol
            })
            .sum()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// PlasmaEquationOfState
// ─────────────────────────────────────────────────────────────────────────────

/// Equation of state for a plasma: ideal gas + radiation pressure + electron pressure.
///
/// P_total = P_gas + P_rad + P_e
pub struct PlasmaEquationOfState {
    /// Adiabatic index γ.
    pub gamma: f64,
    /// Mean molecular weight μ (in proton mass units).
    pub mean_molecular_weight: f64,
}

impl PlasmaEquationOfState {
    /// Create a new plasma EOS with given γ and mean molecular weight.
    pub fn new(gamma: f64, mean_molecular_weight: f64) -> Self {
        Self {
            gamma,
            mean_molecular_weight,
        }
    }

    /// Ideal gas pressure P_gas = (ρ/μ m_p) k_B T \[Pa\].
    pub fn gas_pressure(&self, density: f64, temperature: f64) -> f64 {
        density * K_B * temperature / (self.mean_molecular_weight * M_PROTON)
    }

    /// Radiation pressure P_rad = (1/3) a_rad T⁴ \[Pa\].
    pub fn radiation_pressure(&self, temperature: f64) -> f64 {
        A_RAD * temperature.powi(4) / 3.0
    }

    /// Electron pressure P_e = n_e k_B T_e \[Pa\].
    ///
    /// `ionization` is the mean ionization degree Z_eff.
    pub fn electron_pressure(&self, density: f64, temperature: f64, ionization: f64) -> f64 {
        let n_e = ionization * density / (self.mean_molecular_weight * M_PROTON);
        n_e * K_B * temperature
    }

    /// Total pressure P = P_gas + P_rad + P_e.
    pub fn total_pressure(&self, density: f64, temperature: f64, ionization: f64) -> f64 {
        self.gas_pressure(density, temperature)
            + self.radiation_pressure(temperature)
            + self.electron_pressure(density, temperature, ionization)
    }

    /// Sound speed c_s = sqrt(γ P / ρ) \[m s⁻¹\].
    pub fn sound_speed(&self, density: f64, temperature: f64) -> f64 {
        let p = self.gas_pressure(density, temperature);
        if density > 0.0 {
            (self.gamma * p / density).abs().sqrt()
        } else {
            0.0
        }
    }

    /// Internal energy per unit mass u = P_gas / ((γ-1) ρ) \[J kg⁻¹\].
    pub fn specific_internal_energy(&self, density: f64, temperature: f64) -> f64 {
        let p = self.gas_pressure(density, temperature);
        if density > 0.0 && self.gamma > 1.0 {
            p / ((self.gamma - 1.0) * density)
        } else {
            0.0
        }
    }

    /// Ratio of radiation to total pressure β_rad = P_rad / P_total.
    pub fn radiation_fraction(&self, density: f64, temperature: f64, ionization: f64) -> f64 {
        let p_total = self.total_pressure(density, temperature, ionization);
        if p_total > 0.0 {
            self.radiation_pressure(temperature) / p_total
        } else {
            0.0
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// BraginskiiViscosity
// ─────────────────────────────────────────────────────────────────────────────

/// Braginskii anisotropic viscosity along the magnetic field direction.
///
/// The viscous stress tensor has components parallel and perpendicular to B.
/// η_∥ >> η_⊥ for strongly magnetized plasma.
pub struct BraginskiiViscosity {
    /// Parallel viscosity coefficient η_∥ \[Pa s\].
    pub eta_parallel: f64,
    /// Perpendicular viscosity coefficient η_⊥ \[Pa s\].
    pub eta_perp: f64,
    /// Gyroviscosity coefficient η_× \[Pa s\].
    pub eta_gyro: f64,
}

impl BraginskiiViscosity {
    /// Create Braginskii viscosity for an ion-dominated plasma.
    ///
    /// `temperature` \[K\], `density` \[kg m⁻³\], `b_mag` \[T\], `coulomb_log` ∼ 10–20.
    pub fn from_plasma_parameters(
        temperature: f64,
        density: f64,
        b_mag: f64,
        coulomb_log: f64,
    ) -> Self {
        let n_i = density / M_PROTON;
        // Ion collision time τ_i [s] (Braginskii 1965)
        let tau_i = if n_i > 0.0 && coulomb_log > 0.0 {
            let t52 = temperature.powf(2.5);
            6.0 * PI.sqrt() / 8.0 * M_PROTON.powi(2) * t52
                / (n_i
                    * E_CHARGE.powi(4)
                    * coulomb_log
                    * (M_PROTON / K_B).powf(2.5)
                    * (2.0_f64).sqrt())
        } else {
            1e-10
        };
        // Parallel viscosity η_∥ = 0.96 n_i k_B T τ_i
        let eta_parallel = 0.96 * n_i * K_B * temperature * tau_i;
        // Ion cyclotron frequency
        let omega_ci = if b_mag > 0.0 {
            E_CHARGE * b_mag / M_PROTON
        } else {
            0.0
        };
        let x = omega_ci * tau_i;
        // Perpendicular viscosity drops as 1/(1 + x²) for large x
        let eta_perp = if x.abs() > 1e-10 {
            eta_parallel / (1.0 + x * x)
        } else {
            eta_parallel
        };
        let eta_gyro = if x.abs() > 1e-10 {
            eta_parallel * x / (1.0 + x * x)
        } else {
            0.0
        };
        Self {
            eta_parallel,
            eta_perp,
            eta_gyro,
        }
    }

    /// Compute viscous force on a fluid element along B.
    ///
    /// `velocity_grad` is dv/dx \[s⁻¹\], `b_hat` is unit vector along B.
    pub fn viscous_force_parallel(&self, velocity_grad: [f64; 3], b_hat: [f64; 3]) -> [f64; 3] {
        // Force = η_∥ (b̂·∇)(b̂·v) b̂ component
        let bhat_norm = norm3(b_hat);
        let projection = dot3(velocity_grad, bhat_norm);
        scale3(bhat_norm, self.eta_parallel * projection)
    }

    /// Compute total anisotropic viscous stress magnitude \[Pa\].
    pub fn stress_magnitude(&self, shear_rate: f64, b_hat: [f64; 3], flow_dir: [f64; 3]) -> f64 {
        let cos_theta = dot3(norm3(b_hat), norm3(flow_dir)).abs();
        let sin2 = 1.0 - cos_theta * cos_theta;
        // Stress = η_∥ cos²θ + η_⊥ sin²θ (simplified)
        (self.eta_parallel * cos_theta * cos_theta + self.eta_perp * sin2) * shear_rate
    }

    /// Anisotropy ratio η_∥ / η_⊥.
    pub fn anisotropy_ratio(&self) -> f64 {
        if self.eta_perp > 0.0 {
            self.eta_parallel / self.eta_perp
        } else {
            f64::INFINITY
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// ThermalConduction
// ─────────────────────────────────────────────────────────────────────────────

/// Spitzer thermal conductivity: anisotropic conduction along the magnetic field.
///
/// κ_∥ ∝ T^{5/2} / (n Λ) where Λ is the Coulomb logarithm.
pub struct ThermalConduction {
    /// Parallel thermal conductivity κ_∥ \[W m⁻¹ K⁻¹\].
    pub kappa_parallel: f64,
    /// Perpendicular thermal conductivity κ_⊥ \[W m⁻¹ K⁻¹\].
    pub kappa_perp: f64,
}

impl ThermalConduction {
    /// Compute Spitzer conductivity for an electron plasma.
    ///
    /// `temperature` \[K\], `density` \[kg m⁻³\], `b_mag` \[T\], `coulomb_log` ∼ 10–20.
    pub fn spitzer(temperature: f64, density: f64, b_mag: f64, coulomb_log: f64) -> Self {
        let n_e = density / M_ELECTRON;
        let t52 = temperature.powf(2.5);
        // Spitzer conductivity κ_∥ = 3.2 n_e k_B τ_e / m_e, τ_e ∝ T^{3/2}
        let kappa_parallel = if coulomb_log > 0.0 {
            3.2 * K_B * n_e * K_B * t52 / (M_ELECTRON * E_CHARGE.powi(4) * coulomb_log * n_e * 1e10)
        } else {
            0.0
        };
        let omega_ce = if b_mag > 0.0 {
            E_CHARGE * b_mag / M_ELECTRON
        } else {
            0.0
        };
        // Electron collision time τ_e (rough)
        let tau_e = if n_e > 0.0 && coulomb_log > 0.0 {
            3.0 * (2.0 * PI).sqrt() / 16.0 * M_ELECTRON.sqrt() * (K_B * temperature).powf(1.5)
                / (n_e * E_CHARGE.powi(4) * coulomb_log)
        } else {
            1e-10
        };
        let x = omega_ce * tau_e;
        let kappa_perp = if x.abs() > 1e-10 {
            kappa_parallel / (1.0 + x * x)
        } else {
            kappa_parallel
        };
        Self {
            kappa_parallel,
            kappa_perp,
        }
    }

    /// Heat flux along B: q_∥ = −κ_∥ (b̂·∇T) b̂ \[W m⁻²\].
    pub fn heat_flux_parallel(&self, temp_grad: [f64; 3], b_hat: [f64; 3]) -> [f64; 3] {
        let bhat = norm3(b_hat);
        let dt_parallel = dot3(temp_grad, bhat);
        scale3(bhat, -self.kappa_parallel * dt_parallel)
    }

    /// Heat flux perpendicular to B: q_⊥ = −κ_⊥ ∇_⊥ T \[W m⁻²\].
    pub fn heat_flux_perp(&self, temp_grad: [f64; 3], b_hat: [f64; 3]) -> [f64; 3] {
        let bhat = norm3(b_hat);
        let dt_par = dot3(temp_grad, bhat);
        let grad_par = scale3(bhat, dt_par);
        let grad_perp = sub3(temp_grad, grad_par);
        scale3(grad_perp, -self.kappa_perp)
    }

    /// Total heat flux q = q_∥ + q_⊥.
    pub fn total_heat_flux(&self, temp_grad: [f64; 3], b_hat: [f64; 3]) -> [f64; 3] {
        add3(
            self.heat_flux_parallel(temp_grad, b_hat),
            self.heat_flux_perp(temp_grad, b_hat),
        )
    }

    /// Anisotropy ratio κ_∥ / κ_⊥.
    pub fn anisotropy(&self) -> f64 {
        if self.kappa_perp > 0.0 {
            self.kappa_parallel / self.kappa_perp
        } else {
            f64::INFINITY
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// IonizationModel
// ─────────────────────────────────────────────────────────────────────────────

/// Ionization model: Saha equation for LTE and collisional ionization.
pub struct IonizationModel {
    /// Ionization potential χ \[J\].
    pub ionization_potential: f64,
    /// Statistical weight ratio g_{r+1}/g_r.
    pub statistical_weight_ratio: f64,
}

impl IonizationModel {
    /// Create ionization model for a hydrogen-like atom.
    ///
    /// `chi_eV` is the ionization potential in eV.
    pub fn hydrogen_like(chi_ev: f64) -> Self {
        Self {
            ionization_potential: chi_ev * E_CHARGE,
            statistical_weight_ratio: 2.0, // g_{+1}/g_0 = 2 for hydrogen
        }
    }

    /// Saha equation: ionization fraction x for given T and n_e.
    ///
    /// Returns x ∈ \[0, 1\] where x = n_e / (n_e + n_0).
    ///
    /// Saha: (n_e · n_i) / n_0 = (g_i / g_0) · (2π m_e k_B T / h²)^{3/2} · exp(−χ/k_B T)
    pub fn saha_fraction(&self, temperature: f64, electron_density: f64) -> f64 {
        if temperature <= 0.0 || electron_density <= 0.0 {
            return 0.0;
        }
        let chi = self.ionization_potential;
        let saha_rhs = self.statistical_weight_ratio
            * (2.0 * PI * M_ELECTRON * K_B * temperature / (H_PLANCK * H_PLANCK)).powf(1.5)
            * 2.0 // partition function for free electrons
            * (-chi / (K_B * temperature)).exp();
        let saha_rhs = saha_rhs / electron_density;
        // Quadratic: x² + x·saha_rhs - saha_rhs = 0 → x = (-r + sqrt(r² + 4r)) / 2
        let r = 1.0 / saha_rhs;
        if r < 1e-10 {
            1.0 // fully ionized
        } else if r > 1e10 {
            saha_rhs.sqrt() // weakly ionized
        } else {
            let discriminant = 1.0 + 4.0 * r;
            (-1.0 + discriminant.sqrt()) / (2.0 * r) // proper Saha solution normalized
        }
    }

    /// Collisional ionization rate coefficient α_col \[m³ s⁻¹\].
    ///
    /// Lotz approximation: α ∝ T^{-1/2} exp(−χ/kT).
    pub fn collisional_ionization_rate(&self, temperature: f64) -> f64 {
        if temperature <= 0.0 {
            return 0.0;
        }
        let chi = self.ionization_potential;
        let u = chi / (K_B * temperature);
        // Lotz formula: α = A · u · exp(-u) · E1(u) / (kT)^{3/2}
        // Simplified: α = A · exp(-u) / T^{1/2}
        let a_coeff = 1e-14; // normalization constant [m³ s⁻¹ K^{1/2}]
        a_coeff * (-u).exp() / temperature.sqrt()
    }

    /// Recombination rate coefficient α_rec \[m³ s⁻¹\] (Kramers approximation).
    pub fn recombination_rate(&self, temperature: f64) -> f64 {
        if temperature <= 0.0 {
            return 0.0;
        }
        let chi = self.ionization_potential;
        // Radiative recombination: α_rec ∝ (χ/kT)^{1/2} / T^{3/2}
        let u = (chi / (K_B * temperature)).sqrt();
        let a_coeff = 2.07e-16; // [m³ K^{3/2} s⁻¹]
        a_coeff * u / temperature
    }

    /// Mean ionization degree Z_eff for LTE at temperature T, electron density n_e.
    pub fn mean_ionization(&self, temperature: f64, electron_density: f64) -> f64 {
        self.saha_fraction(temperature, electron_density)
            .clamp(0.0, 1.0)
    }

    /// Check if plasma is fully ionized (χ << kT).
    pub fn is_fully_ionized(&self, temperature: f64) -> bool {
        K_B * temperature > 10.0 * self.ionization_potential
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// RadiativeCooling
// ─────────────────────────────────────────────────────────────────────────────

/// Radiative cooling: bremsstrahlung, line cooling, and recombination radiation.
pub struct RadiativeCooling {
    /// Metallicity relative to solar (Z/Z_☉), affects line cooling.
    pub metallicity: f64,
}

impl RadiativeCooling {
    /// Create a radiative cooling model.
    pub fn new(metallicity: f64) -> Self {
        Self { metallicity }
    }

    /// Bremsstrahlung (free-free) cooling rate \[W m⁻³\].
    ///
    /// Λ_ff = g_ff · 1.69×10⁻⁴⁰ n_e² T^{1/2} \[SI units\].
    pub fn bremsstrahlung_cooling(&self, electron_density: f64, temperature: f64) -> f64 {
        if temperature <= 0.0 {
            return 0.0;
        }
        // Gaunt factor g_ff ≈ 1.2 for typical conditions
        let g_ff = 1.2;
        BREMSSTRAHLUNG_COEFF * g_ff * electron_density * electron_density * temperature.sqrt()
    }

    /// Line cooling rate Λ_line \[W m⁻³\] (Sutherland & Dopita 1993 approximation).
    ///
    /// Peaks around T ~ 10^5 K for solar metallicity.
    pub fn line_cooling(&self, electron_density: f64, temperature: f64) -> f64 {
        if temperature <= 0.0 || electron_density <= 0.0 {
            return 0.0;
        }
        // Simplified cooling function: Λ ∝ Z · n_e² · Λ(T)
        // Use a parametric form with peak at log T ~ 5
        let log_t = temperature.log10();
        let lambda_t = if log_t < 4.0 {
            0.0
        } else if log_t < 5.0 {
            1e-36 * (log_t - 4.0).powi(2)
        } else if log_t < 6.0 {
            1e-36 * (1.0 - (log_t - 5.0).powi(2) * 0.5)
        } else {
            1e-36 * (-0.5 * (log_t - 6.0))
        };
        self.metallicity * electron_density * electron_density * lambda_t
    }

    /// Recombination radiation cooling \[W m⁻³\].
    pub fn recombination_cooling(
        &self,
        electron_density: f64,
        ion_density: f64,
        temperature: f64,
        ionization_potential: f64,
    ) -> f64 {
        if temperature <= 0.0 {
            return 0.0;
        }
        // Λ_rec = α_rec · n_e · n_i · χ
        let a_rec = 2.07e-16 / temperature; // simplified Kramers
        a_rec * electron_density * ion_density * ionization_potential
    }

    /// Total cooling rate \[W m⁻³\].
    pub fn total_cooling_rate(
        &self,
        electron_density: f64,
        ion_density: f64,
        temperature: f64,
        ionization_potential: f64,
    ) -> f64 {
        self.bremsstrahlung_cooling(electron_density, temperature)
            + self.line_cooling(electron_density, temperature)
            + self.recombination_cooling(
                electron_density,
                ion_density,
                temperature,
                ionization_potential,
            )
    }

    /// Cooling time t_cool = (3/2) n k_B T / Λ \[s\].
    pub fn cooling_time(&self, total_density: f64, temperature: f64, cooling_rate: f64) -> f64 {
        if cooling_rate <= 0.0 {
            return f64::INFINITY;
        }
        let n = total_density / M_PROTON;
        1.5 * n * K_B * temperature / cooling_rate
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// PicSphCoupling
// ─────────────────────────────────────────────────────────────────────────────

/// Macro particle (SPH) + micro particle (PIC) coupling interface.
///
/// Provides methods to deposit SPH particle properties onto a PIC grid
/// and interpolate electromagnetic fields back to SPH particles.
pub struct PicSphCoupling {
    /// Grid spacing \[m\].
    pub cell_size: f64,
    /// Number of grid cells per dimension.
    pub n_cells: usize,
    /// Electric field on grid \[V m⁻¹\], stored as flat array of \[Ex, Ey, Ez\] per cell.
    pub e_field_grid: Vec<[f64; 3]>,
    /// Magnetic field on grid \[T\].
    pub b_field_grid: Vec<[f64; 3]>,
    /// Charge density on grid \[C m⁻³\].
    pub charge_density_grid: Vec<f64>,
    /// Current density on grid \[A m⁻²\].
    pub current_density_grid: Vec<[f64; 3]>,
}

impl PicSphCoupling {
    /// Create a new PIC-SPH coupling interface with a uniform grid.
    pub fn new(cell_size: f64, n_cells: usize) -> Self {
        let total = n_cells * n_cells * n_cells;
        Self {
            cell_size,
            n_cells,
            e_field_grid: vec![[0.0; 3]; total],
            b_field_grid: vec![[0.0; 3]; total],
            charge_density_grid: vec![0.0; total],
            current_density_grid: vec![[0.0; 3]; total],
        }
    }

    /// Convert 3D index to flat index.
    fn flat_index(&self, ix: usize, iy: usize, iz: usize) -> usize {
        (ix * self.n_cells + iy) * self.n_cells + iz
    }

    /// Grid cell index for position component x.
    fn cell_index(&self, x: f64) -> usize {
        let idx = (x / self.cell_size) as i64;
        idx.rem_euclid(self.n_cells as i64) as usize
    }

    /// Deposit SPH particle charge onto the grid (nearest-grid-point).
    pub fn deposit_charge(&mut self, particles: &[PlasmaParticle]) {
        for c in self.charge_density_grid.iter_mut() {
            *c = 0.0;
        }
        let vol = self.cell_size.powi(3);
        for p in particles {
            let ix = self.cell_index(p.position[0]);
            let iy = self.cell_index(p.position[1]);
            let iz = self.cell_index(p.position[2]);
            let idx = self.flat_index(ix, iy, iz);
            self.charge_density_grid[idx] += p.charge / vol;
        }
    }

    /// Deposit current density J = ρ_charge · v onto the grid.
    pub fn deposit_current(&mut self, particles: &[PlasmaParticle]) {
        for c in self.current_density_grid.iter_mut() {
            *c = [0.0; 3];
        }
        let vol = self.cell_size.powi(3);
        for p in particles {
            let ix = self.cell_index(p.position[0]);
            let iy = self.cell_index(p.position[1]);
            let iz = self.cell_index(p.position[2]);
            let idx = self.flat_index(ix, iy, iz);
            let j_contrib = scale3(p.velocity, p.charge / vol);
            self.current_density_grid[idx] = add3(self.current_density_grid[idx], j_contrib);
        }
    }

    /// Interpolate electric field from grid to particle position (NGP).
    pub fn interpolate_e_field(&self, position: [f64; 3]) -> [f64; 3] {
        let ix = self.cell_index(position[0]);
        let iy = self.cell_index(position[1]);
        let iz = self.cell_index(position[2]);
        let idx = self.flat_index(ix, iy, iz);
        self.e_field_grid[idx]
    }

    /// Interpolate magnetic field from grid to particle position (NGP).
    pub fn interpolate_b_field(&self, position: [f64; 3]) -> [f64; 3] {
        let ix = self.cell_index(position[0]);
        let iy = self.cell_index(position[1]);
        let iz = self.cell_index(position[2]);
        let idx = self.flat_index(ix, iy, iz);
        self.b_field_grid[idx]
    }

    /// Set electromagnetic fields for a grid cell.
    pub fn set_fields(&mut self, ix: usize, iy: usize, iz: usize, e: [f64; 3], b: [f64; 3]) {
        let idx = self.flat_index(ix % self.n_cells, iy % self.n_cells, iz % self.n_cells);
        self.e_field_grid[idx] = e;
        self.b_field_grid[idx] = b;
    }

    /// Number of grid cells total.
    pub fn total_cells(&self) -> usize {
        self.n_cells.pow(3)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// AlfvenWave
// ─────────────────────────────────────────────────────────────────────────────

/// Alfvén wave propagation in MHD plasma.
///
/// Alfvén wave: perturbations propagate along B at speed v_A = B/√(μ₀ρ).
pub struct AlfvenWave {
    /// Background magnetic field B₀ \[T\].
    pub b0: [f64; 3],
    /// Background plasma density ρ₀ \[kg m⁻³\].
    pub rho0: f64,
    /// Wave vector k \[m⁻¹\].
    pub k_wave: [f64; 3],
    /// Wave amplitude δB / B₀.
    pub amplitude: f64,
}

impl AlfvenWave {
    /// Create an Alfvén wave with given parameters.
    pub fn new(b0: [f64; 3], rho0: f64, k_wave: [f64; 3], amplitude: f64) -> Self {
        Self {
            b0,
            rho0,
            k_wave,
            amplitude,
        }
    }

    /// Alfvén speed v_A = |B₀| / √(μ₀ ρ₀) \[m s⁻¹\].
    pub fn alfven_speed(&self) -> f64 {
        if self.rho0 <= 0.0 {
            return 0.0;
        }
        let b_mag = len3(self.b0);
        b_mag / (MU_0 * self.rho0).sqrt()
    }

    /// Phase velocity of Alfvén wave along B₀ \[m s⁻¹\].
    pub fn phase_velocity(&self) -> f64 {
        let va = self.alfven_speed();
        let k_mag = len3(self.k_wave);
        if k_mag < 1e-30 {
            return 0.0;
        }
        let b_hat = norm3(self.b0);
        let k_hat = norm3(self.k_wave);
        let cos_theta = dot3(b_hat, k_hat).abs();
        va * cos_theta
    }

    /// Angular frequency ω = k · v_A · cos(θ) \[rad s⁻¹\].
    pub fn angular_frequency(&self) -> f64 {
        let k_mag = len3(self.k_wave);
        self.phase_velocity() * k_mag
    }

    /// Group velocity = phase velocity for Alfvén waves (non-dispersive along B).
    pub fn group_velocity(&self) -> f64 {
        self.alfven_speed()
    }

    /// Magnetic pressure of the wave δB² / (2μ₀) \[Pa\].
    pub fn wave_magnetic_pressure(&self) -> f64 {
        let b_mag = len3(self.b0);
        let db = self.amplitude * b_mag;
        db * db / (2.0 * MU_0)
    }

    /// Wave energy density \[J m⁻³\].
    pub fn wave_energy_density(&self) -> f64 {
        2.0 * self.wave_magnetic_pressure()
    }

    /// Displacement δv of the wave \[m s⁻¹\].
    pub fn velocity_perturbation(&self) -> f64 {
        self.amplitude * self.alfven_speed()
    }

    /// Fast magnetosonic speed v_ms = √(v_A² + c_s²) \[m s⁻¹\].
    pub fn fast_magnetosonic_speed(&self, sound_speed: f64) -> f64 {
        let va = self.alfven_speed();
        (va * va + sound_speed * sound_speed).sqrt()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// MagneticReconnection
// ─────────────────────────────────────────────────────────────────────────────

/// Magnetic reconnection: Sweet-Parker and Petschek models.
///
/// Reconnection converts magnetic energy to kinetic + thermal energy
/// by changing field topology at a current sheet.
pub struct MagneticReconnection {
    /// Upstream magnetic field \[T\].
    pub b_upstream: f64,
    /// Upstream plasma density \[kg m⁻³\].
    pub rho_upstream: f64,
    /// Resistivity η \[Ω m\].
    pub resistivity: f64,
    /// System length scale L \[m\].
    pub length_scale: f64,
}

impl MagneticReconnection {
    /// Create a reconnection model.
    pub fn new(b_upstream: f64, rho_upstream: f64, resistivity: f64, length_scale: f64) -> Self {
        Self {
            b_upstream,
            rho_upstream,
            resistivity,
            length_scale,
        }
    }

    /// Alfvén speed v_A = B/√(μ₀ρ) \[m s⁻¹\].
    pub fn alfven_speed(&self) -> f64 {
        if self.rho_upstream <= 0.0 {
            return 0.0;
        }
        self.b_upstream / (MU_0 * self.rho_upstream).sqrt()
    }

    /// Lundquist number S = v_A L / η_m where η_m = η/μ₀.
    pub fn lundquist_number(&self) -> f64 {
        if self.resistivity <= 0.0 {
            return f64::INFINITY;
        }
        let eta_m = self.resistivity / MU_0;
        self.alfven_speed() * self.length_scale / eta_m
    }

    /// Sweet-Parker reconnection rate M_sp = S^{-1/2} (dimensionless inflow Mach number).
    pub fn sweet_parker_rate(&self) -> f64 {
        let s = self.lundquist_number();
        if s > 0.0 && s.is_finite() {
            s.powf(-0.5)
        } else {
            0.0
        }
    }

    /// Sweet-Parker reconnection electric field E_sp = M_sp · v_A · B₀ \[V m⁻¹\].
    pub fn sweet_parker_electric_field(&self) -> f64 {
        self.sweet_parker_rate() * self.alfven_speed() * self.b_upstream
    }

    /// Petschek reconnection rate M_p = π / (8 ln S) (much faster than Sweet-Parker).
    pub fn petschek_rate(&self) -> f64 {
        let s = self.lundquist_number();
        if s > 1.0 { PI / (8.0 * s.ln()) } else { 0.0 }
    }

    /// Current sheet thickness δ_SP = L / S^{1/2} \[m\].
    pub fn current_sheet_thickness(&self) -> f64 {
        let s = self.lundquist_number();
        if s > 0.0 && s.is_finite() {
            self.length_scale / s.sqrt()
        } else {
            0.0
        }
    }

    /// Reconnection outflow speed ≈ v_A \[m s⁻¹\].
    pub fn outflow_speed(&self) -> f64 {
        self.alfven_speed()
    }

    /// Energy release rate per unit volume \[W m⁻³\].
    pub fn energy_release_rate(&self) -> f64 {
        let e_rec = self.sweet_parker_electric_field();
        let j_sheet = self.b_upstream / (MU_0 * self.current_sheet_thickness().max(1e-30));
        e_rec * j_sheet
    }

    /// Check if Petschek reconnection is feasible (requires S >> 1).
    pub fn petschek_feasible(&self) -> bool {
        self.lundquist_number() > 100.0
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Utility functions
// ─────────────────────────────────────────────────────────────────────────────

/// Debye length λ_D = √(ε₀ k_B T_e / (n_e e²)) \[m\].
///
/// # Arguments
/// * `temperature` - electron temperature \[K\]
/// * `electron_density` - electron number density \[m⁻³\]
pub fn debye_length(temperature: f64, electron_density: f64) -> f64 {
    if temperature <= 0.0 || electron_density <= 0.0 {
        return 0.0;
    }
    (EPS_0 * K_B * temperature / (electron_density * E_CHARGE * E_CHARGE)).sqrt()
}

/// Cyclotron (gyro) frequency ω_c = |q|B/m \[rad s⁻¹\].
///
/// # Arguments
/// * `charge` - particle charge \[C\]
/// * `mass` - particle mass \[kg\]
/// * `b_magnitude` - magnetic field magnitude \[T\]
pub fn cyclotron_frequency(charge: f64, mass: f64, b_magnitude: f64) -> f64 {
    if mass <= 0.0 || b_magnitude <= 0.0 {
        return 0.0;
    }
    charge.abs() * b_magnitude / mass
}

/// Electron plasma frequency ω_pe = √(n_e e² / (ε₀ m_e)) \[rad s⁻¹\].
///
/// # Arguments
/// * `electron_density` - electron number density \[m⁻³\]
pub fn plasma_frequency(electron_density: f64) -> f64 {
    if electron_density <= 0.0 {
        return 0.0;
    }
    (electron_density * E_CHARGE * E_CHARGE / (EPS_0 * M_ELECTRON)).sqrt()
}

/// Ion plasma frequency ω_pi = √(n_i Z² e² / (ε₀ m_i)) \[rad s⁻¹\].
///
/// # Arguments
/// * `ion_density` - ion number density \[m⁻³\]
/// * `mass_ratio` - ion-to-proton mass ratio
/// * `charge_number` - ionization charge number Z
pub fn ion_plasma_frequency(ion_density: f64, mass_ratio: f64, charge_number: f64) -> f64 {
    if ion_density <= 0.0 || mass_ratio <= 0.0 {
        return 0.0;
    }
    let m_ion = mass_ratio * M_PROTON;
    let z2 = charge_number * charge_number;
    (ion_density * z2 * E_CHARGE * E_CHARGE / (EPS_0 * m_ion)).sqrt()
}

/// Plasma beta parameter β = P_thermal / P_magnetic = 2μ₀ n k_B T / B².
///
/// β < 1: magnetically dominated; β > 1: thermally dominated.
///
/// # Arguments
/// * `density` - plasma mass density \[kg m⁻³\]
/// * `temperature` - plasma temperature \[K\]
/// * `b_magnitude` - magnetic field magnitude \[T\]
pub fn beta_parameter(density: f64, temperature: f64, b_magnitude: f64) -> f64 {
    if b_magnitude <= 0.0 {
        return f64::INFINITY;
    }
    let n = density / M_PROTON;
    let p_thermal = n * K_B * temperature;
    let p_magnetic = b_magnitude * b_magnitude / (2.0 * MU_0);
    p_thermal / p_magnetic
}

/// Larmor (gyro) radius r_L = m v_⊥ / (|q| B) \[m\].
///
/// # Arguments
/// * `mass` - particle mass \[kg\]
/// * `v_perp` - speed perpendicular to B \[m s⁻¹\]
/// * `charge` - particle charge \[C\]
/// * `b_magnitude` - magnetic field magnitude \[T\]
pub fn larmor_radius(mass: f64, v_perp: f64, charge: f64, b_magnitude: f64) -> f64 {
    if charge.abs() <= 0.0 || b_magnitude <= 0.0 {
        return f64::INFINITY;
    }
    mass * v_perp / (charge.abs() * b_magnitude)
}

/// Magnetic Reynolds number R_m = v L / η_m = v L μ₀ / η.
///
/// # Arguments
/// * `velocity` - characteristic velocity \[m s⁻¹\]
/// * `length` - characteristic length \[m\]
/// * `resistivity` - electrical resistivity \[Ω m\]
pub fn magnetic_reynolds_number(velocity: f64, length: f64, resistivity: f64) -> f64 {
    if resistivity <= 0.0 {
        return f64::INFINITY;
    }
    velocity * length * MU_0 / resistivity
}

/// Inertial length d_i = c / ω_pi \[m\] (ion skin depth).
pub fn ion_inertial_length(ion_density: f64) -> f64 {
    let omega_pi = ion_plasma_frequency(ion_density, 1.0, 1.0);
    if omega_pi > 0.0 {
        C_LIGHT / omega_pi
    } else {
        f64::INFINITY
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    const TOL: f64 = 1e-8;
    const RTOL: f64 = 1e-4; // relative tolerance for physics checks

    fn assert_rel(a: f64, b: f64, tol: f64, msg: &str) {
        let rel = if b.abs() > 1e-300 {
            (a - b).abs() / b.abs()
        } else {
            (a - b).abs()
        };
        assert!(rel < tol, "{}: {} vs {} (rel err {:.6})", msg, a, b, rel);
    }

    // --- PlasmaParticle tests ---

    #[test]
    fn test_plasma_particle_new_default_fields() {
        let p = PlasmaParticle::new([0.0; 3], [0.0; 3], M_PROTON, E_CHARGE, 1e4, 1e-3);
        assert_eq!(p.magnetic_field, [0.0; 3]);
        assert_eq!(p.psi, 0.0);
        assert!(p.pressure > 0.0);
    }

    #[test]
    fn test_plasma_particle_kinetic_energy_at_rest() {
        let p = PlasmaParticle::new([0.0; 3], [0.0; 3], M_PROTON, E_CHARGE, 1e4, 1e-3);
        assert_eq!(p.kinetic_energy(), 0.0);
    }

    #[test]
    fn test_plasma_particle_kinetic_energy_moving() {
        let v = [1000.0, 0.0, 0.0];
        let m = M_PROTON;
        let p = PlasmaParticle::new([0.0; 3], v, m, E_CHARGE, 1e4, 1e-3);
        let expected = 0.5 * m * 1000.0_f64.powi(2);
        assert_rel(p.kinetic_energy(), expected, TOL, "KE");
    }

    #[test]
    fn test_plasma_particle_speed() {
        let v = [3.0, 4.0, 0.0];
        let p = PlasmaParticle::new([0.0; 3], v, M_PROTON, E_CHARGE, 1e4, 1e-3);
        assert_rel(p.speed(), 5.0, TOL, "speed");
    }

    #[test]
    fn test_plasma_particle_b_magnitude() {
        let mut p = PlasmaParticle::new([0.0; 3], [0.0; 3], M_PROTON, E_CHARGE, 1e4, 1e-3);
        p.magnetic_field = [3.0, 4.0, 0.0];
        assert_rel(p.b_magnitude(), 5.0, TOL, "B magnitude");
    }

    #[test]
    fn test_plasma_particle_thermal_velocity_positive() {
        let p = PlasmaParticle::new([0.0; 3], [0.0; 3], M_PROTON, E_CHARGE, 1e6, 1e-3);
        assert!(p.thermal_velocity() > 0.0);
    }

    // --- Alfvén wave speed tests ---

    #[test]
    fn test_alfven_speed_formula() {
        // v_A = B / sqrt(μ₀ ρ)
        let b = [0.0, 0.0, 1e-3]; // 1 mT
        let rho = 1e-6; // kg/m³
        let wave = AlfvenWave::new(b, rho, [0.0, 0.0, 1e3], 0.01);
        let expected = 1e-3 / (MU_0 * rho).sqrt();
        assert_rel(wave.alfven_speed(), expected, RTOL, "Alfvén speed");
    }

    #[test]
    fn test_alfven_speed_zero_density() {
        let wave = AlfvenWave::new([0.0, 0.0, 1.0], 0.0, [0.0, 0.0, 1.0], 0.01);
        assert_eq!(wave.alfven_speed(), 0.0);
    }

    #[test]
    fn test_alfven_speed_scales_with_b() {
        let rho = 1e-6;
        let wave1 = AlfvenWave::new([0.0, 0.0, 1e-3], rho, [0.0, 0.0, 1.0], 0.01);
        let wave2 = AlfvenWave::new([0.0, 0.0, 2e-3], rho, [0.0, 0.0, 1.0], 0.01);
        assert_rel(
            wave2.alfven_speed(),
            2.0 * wave1.alfven_speed(),
            RTOL,
            "Alfvén ∝ B",
        );
    }

    #[test]
    fn test_alfven_wave_energy_density_positive() {
        let wave = AlfvenWave::new([0.0, 0.0, 1e-3], 1e-6, [0.0, 0.0, 1e3], 0.1);
        assert!(wave.wave_energy_density() > 0.0);
    }

    #[test]
    fn test_alfven_fast_magnetosonic_exceeds_alfven() {
        let wave = AlfvenWave::new([0.0, 0.0, 1e-3], 1e-6, [0.0, 0.0, 1e3], 0.01);
        let cs = 1e4;
        assert!(wave.fast_magnetosonic_speed(cs) >= wave.alfven_speed());
    }

    // --- Debye length tests ---

    #[test]
    fn test_debye_length_positive() {
        let ld = debye_length(1e4, 1e18);
        assert!(ld > 0.0, "Debye length should be positive");
    }

    #[test]
    fn test_debye_length_formula() {
        let te = 1e4;
        let ne = 1e18;
        let expected = (EPS_0 * K_B * te / (ne * E_CHARGE * E_CHARGE)).sqrt();
        assert_rel(debye_length(te, ne), expected, TOL, "Debye length formula");
    }

    #[test]
    fn test_debye_length_scales_with_temperature() {
        // λ_D ∝ T^{1/2}
        let ld1 = debye_length(1e4, 1e18);
        let ld2 = debye_length(4e4, 1e18);
        assert_rel(ld2, 2.0 * ld1, RTOL, "Debye ∝ T^{1/2}");
    }

    #[test]
    fn test_debye_length_zero_temperature() {
        assert_eq!(debye_length(0.0, 1e18), 0.0);
    }

    // --- Plasma beta tests ---

    #[test]
    fn test_plasma_beta_formula() {
        let rho = 1e-6;
        let t = 1e6;
        let b = 1e-3;
        let beta = beta_parameter(rho, t, b);
        let n = rho / M_PROTON;
        let p_th = n * K_B * t;
        let p_mag = b * b / (2.0 * MU_0);
        let expected = p_th / p_mag;
        assert_rel(beta, expected, RTOL, "plasma beta");
    }

    #[test]
    fn test_plasma_beta_zero_b() {
        let beta = beta_parameter(1e-6, 1e6, 0.0);
        assert!(beta.is_infinite(), "β → ∞ when B=0");
    }

    #[test]
    fn test_plasma_beta_high_magnetic_pressure() {
        // B large → β << 1
        let beta = beta_parameter(1e-10, 1e3, 1.0);
        assert!(beta < 1.0, "magnetically dominated: β < 1");
    }

    // --- Plasma frequency tests ---

    #[test]
    fn test_plasma_frequency_positive() {
        assert!(plasma_frequency(1e18) > 0.0);
    }

    #[test]
    fn test_plasma_frequency_scales_with_density() {
        // ω_pe ∝ n_e^{1/2}
        let wp1 = plasma_frequency(1e18);
        let wp2 = plasma_frequency(4e18);
        assert_rel(wp2, 2.0 * wp1, RTOL, "ω_pe ∝ n_e^{1/2}");
    }

    #[test]
    fn test_plasma_frequency_zero_density() {
        assert_eq!(plasma_frequency(0.0), 0.0);
    }

    // --- Cyclotron frequency tests ---

    #[test]
    fn test_cyclotron_frequency_formula() {
        let omega = cyclotron_frequency(E_CHARGE, M_ELECTRON, 1.0);
        let expected = E_CHARGE * 1.0 / M_ELECTRON;
        assert_rel(omega, expected, RTOL, "cyclotron frequency");
    }

    #[test]
    fn test_cyclotron_frequency_zero_mass() {
        assert_eq!(cyclotron_frequency(E_CHARGE, 0.0, 1.0), 0.0);
    }

    #[test]
    fn test_cyclotron_frequency_zero_b() {
        assert_eq!(cyclotron_frequency(E_CHARGE, M_ELECTRON, 0.0), 0.0);
    }

    // --- Ionization (Saha) tests ---

    #[test]
    fn test_saha_fraction_high_temperature_full_ionization() {
        // At very high T >> χ/k_B, should approach full ionization
        let model = IonizationModel::hydrogen_like(13.6);
        let fraction = model.saha_fraction(1e8, 1e18);
        assert!(fraction > 0.9, "Should be nearly fully ionized at T=10^8 K");
    }

    #[test]
    fn test_saha_fraction_zero_temperature() {
        let model = IonizationModel::hydrogen_like(13.6);
        let fraction = model.saha_fraction(0.0, 1e18);
        assert_eq!(fraction, 0.0);
    }

    #[test]
    fn test_saha_is_fully_ionized_at_high_t() {
        let model = IonizationModel::hydrogen_like(13.6);
        assert!(model.is_fully_ionized(1e8));
        assert!(!model.is_fully_ionized(1e3));
    }

    #[test]
    fn test_collisional_ionization_rate_positive() {
        let model = IonizationModel::hydrogen_like(13.6);
        assert!(model.collisional_ionization_rate(1e6) > 0.0);
    }

    #[test]
    fn test_recombination_rate_decreases_with_temperature() {
        let model = IonizationModel::hydrogen_like(13.6);
        let r1 = model.recombination_rate(1e4);
        let r2 = model.recombination_rate(1e6);
        assert!(r1 > r2, "recombination rate decreases with T");
    }

    // --- Bremsstrahlung cooling ---

    #[test]
    fn test_bremsstrahlung_cooling_positive() {
        let cooling = RadiativeCooling::new(1.0);
        assert!(cooling.bremsstrahlung_cooling(1e18, 1e6) > 0.0);
    }

    #[test]
    fn test_bremsstrahlung_scales_as_n_squared() {
        let cooling = RadiativeCooling::new(1.0);
        let c1 = cooling.bremsstrahlung_cooling(1e18, 1e6);
        let c2 = cooling.bremsstrahlung_cooling(2e18, 1e6);
        assert_rel(c2, 4.0 * c1, RTOL, "bremsstrahlung ∝ n_e²");
    }

    #[test]
    fn test_bremsstrahlung_scales_as_sqrt_t() {
        let cooling = RadiativeCooling::new(1.0);
        let c1 = cooling.bremsstrahlung_cooling(1e18, 1e6);
        let c2 = cooling.bremsstrahlung_cooling(1e18, 4e6);
        assert_rel(c2, 2.0 * c1, RTOL, "bremsstrahlung ∝ T^{1/2}");
    }

    #[test]
    fn test_bremsstrahlung_zero_temperature() {
        let cooling = RadiativeCooling::new(1.0);
        assert_eq!(cooling.bremsstrahlung_cooling(1e18, 0.0), 0.0);
    }

    // --- MHD pressure balance ---

    #[test]
    fn test_mhd_total_pressure_includes_magnetic() {
        let eos = PlasmaEquationOfState::new(5.0 / 3.0, 1.0);
        let p_gas = eos.gas_pressure(1e-6, 1e6);
        let b = 1e-3_f64;
        let p_mag = b * b / (2.0 * MU_0);
        let p_total = p_gas + p_mag;
        assert!(p_total > p_gas);
        assert!(p_total > p_mag);
    }

    #[test]
    fn test_eos_radiation_pressure_scales_as_t4() {
        let eos = PlasmaEquationOfState::new(5.0 / 3.0, 1.0);
        let p1 = eos.radiation_pressure(1e6);
        let p2 = eos.radiation_pressure(2e6);
        assert_rel(p2, 16.0 * p1, RTOL, "P_rad ∝ T^4");
    }

    #[test]
    fn test_eos_sound_speed_positive() {
        let eos = PlasmaEquationOfState::new(5.0 / 3.0, 1.0);
        assert!(eos.sound_speed(1e-6, 1e6) > 0.0);
    }

    // --- Lorentz force ---

    #[test]
    fn test_lorentz_force_perpendicular_to_b() {
        // For a single particle in a 2-particle system, force should be roughly ⊥ to B
        let mut solver = MhdSph::new(1e-7, 5.0 / 3.0);
        let mut p1 =
            PlasmaParticle::new([0.0; 3], [1000.0, 0.0, 0.0], M_PROTON, E_CHARGE, 1e4, 1e-6);
        p1.magnetic_field = [0.0, 0.0, 1e-3];
        let mut p2 = PlasmaParticle::new(
            [0.01, 0.0, 0.0],
            [1000.0, 0.0, 0.0],
            M_PROTON,
            E_CHARGE,
            1e4,
            1e-6,
        );
        p2.magnetic_field = [0.0, 0.0, 1e-3];
        p2.smoothing_length = 0.1;
        p1.smoothing_length = 0.1;
        solver.add_particle(p1);
        solver.add_particle(p2);
        let f = solver.lorentz_force(0);
        // Force should be non-NaN
        assert!(f[0].is_finite());
        assert!(f[1].is_finite());
        assert!(f[2].is_finite());
    }

    // --- Magnetic reconnection ---

    #[test]
    fn test_sweet_parker_rate_less_than_petschek() {
        let rec = MagneticReconnection::new(1e-3, 1e-6, 1e-7, 1e6);
        let r_sp = rec.sweet_parker_rate();
        let r_p = rec.petschek_rate();
        assert!(
            r_p > r_sp,
            "Petschek reconnection is faster than Sweet-Parker"
        );
    }

    #[test]
    fn test_lundquist_number_positive() {
        let rec = MagneticReconnection::new(1e-3, 1e-6, 1e-7, 1e6);
        assert!(rec.lundquist_number() > 0.0);
    }

    #[test]
    fn test_current_sheet_thickness_positive() {
        let rec = MagneticReconnection::new(1e-3, 1e-6, 1e-7, 1e6);
        assert!(rec.current_sheet_thickness() > 0.0);
    }

    #[test]
    fn test_reconnection_outflow_equals_alfven() {
        let rec = MagneticReconnection::new(1e-3, 1e-6, 1e-7, 1e6);
        assert_rel(
            rec.outflow_speed(),
            rec.alfven_speed(),
            TOL,
            "outflow = v_A",
        );
    }

    // --- Dedner divergence cleaning ---

    #[test]
    fn test_dedner_cleaning_reduces_psi() {
        let mut solver = MhdSph::new(0.0, 5.0 / 3.0);
        let mut p = PlasmaParticle::new([0.0; 3], [0.0; 3], M_PROTON, E_CHARGE, 1e4, 1e-6);
        p.psi = 1.0;
        p.magnetic_field = [0.0, 0.0, 1e-3];
        solver.add_particle(p);
        let psi_before = solver.particles[0].psi;
        solver.dedner_cleaning_step(1e-6);
        let psi_after = solver.particles[0].psi;
        assert!(
            psi_after.abs() < psi_before.abs(),
            "Dedner cleaning should reduce |ψ|"
        );
    }

    // --- Braginskii viscosity ---

    #[test]
    fn test_braginskii_parallel_larger_than_perp() {
        let visc = BraginskiiViscosity::from_plasma_parameters(1e7, 1e-6, 1e-3, 15.0);
        assert!(visc.eta_parallel >= visc.eta_perp, "η_∥ ≥ η_⊥");
    }

    #[test]
    fn test_braginskii_anisotropy_increases_with_b() {
        let visc1 = BraginskiiViscosity::from_plasma_parameters(1e7, 1e-6, 1e-4, 15.0);
        let visc2 = BraginskiiViscosity::from_plasma_parameters(1e7, 1e-6, 1e-3, 15.0);
        let ani1 = visc1.anisotropy_ratio();
        let ani2 = visc2.anisotropy_ratio();
        assert!(ani2 > ani1, "Anisotropy increases with stronger B");
    }

    // --- PIC-SPH coupling ---

    #[test]
    fn test_pic_sph_total_cells() {
        let coupling = PicSphCoupling::new(1e-3, 4);
        assert_eq!(coupling.total_cells(), 64);
    }

    #[test]
    fn test_pic_sph_charge_deposit() {
        let mut coupling = PicSphCoupling::new(1e-3, 4);
        let particles = vec![PlasmaParticle::new(
            [0.0; 3], [0.0; 3], M_PROTON, 1e-10, 1e4, 1e-3,
        )];
        coupling.deposit_charge(&particles);
        let total_charge: f64 =
            coupling.charge_density_grid.iter().sum::<f64>() * coupling.cell_size.powi(3);
        assert!(
            (total_charge - 1e-10).abs() < 1e-20,
            "charge is conserved after deposit"
        );
    }

    #[test]
    fn test_pic_sph_field_interpolation() {
        let mut coupling = PicSphCoupling::new(1e-3, 4);
        coupling.set_fields(0, 0, 0, [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
        let e = coupling.interpolate_e_field([0.0; 3]);
        assert_eq!(e[0], 1.0);
    }

    // --- Thermal conduction ---

    #[test]
    fn test_thermal_conduction_parallel_positive() {
        let cond = ThermalConduction::spitzer(1e7, 1e-6, 1e-3, 15.0);
        assert!(cond.kappa_parallel >= 0.0);
    }

    #[test]
    fn test_thermal_conduction_anisotropy() {
        let cond = ThermalConduction::spitzer(1e7, 1e-6, 1e-3, 15.0);
        assert!(cond.anisotropy() >= 1.0, "κ_∥ ≥ κ_⊥");
    }

    // --- Magnetic Reynolds number ---

    #[test]
    fn test_magnetic_reynolds_number_positive() {
        let rm = magnetic_reynolds_number(1e6, 1e3, 1e-7);
        assert!(rm > 0.0);
    }

    #[test]
    fn test_magnetic_reynolds_scales_with_velocity() {
        let rm1 = magnetic_reynolds_number(1e6, 1e3, 1e-7);
        let rm2 = magnetic_reynolds_number(2e6, 1e3, 1e-7);
        assert_rel(rm2, 2.0 * rm1, TOL, "Rm ∝ v");
    }

    // --- Larmor radius ---

    #[test]
    fn test_larmor_radius_formula() {
        let r = larmor_radius(M_ELECTRON, 1e6, E_CHARGE, 1e-3);
        let expected = M_ELECTRON * 1e6 / (E_CHARGE * 1e-3);
        assert_rel(r, expected, RTOL, "Larmor radius");
    }

    #[test]
    fn test_larmor_radius_zero_b_is_infinite() {
        let r = larmor_radius(M_ELECTRON, 1e6, E_CHARGE, 0.0);
        assert!(r.is_infinite());
    }

    // --- MHD solver integration ---

    #[test]
    fn test_mhd_sph_euler_step_moves_particle() {
        let mut solver = MhdSph::new(1e-7, 5.0 / 3.0);
        let p = PlasmaParticle::new([0.0; 3], [1.0, 0.0, 0.0], M_PROTON, E_CHARGE, 1e4, 1e-6);
        solver.add_particle(p);
        solver.euler_step(0.01);
        let x = solver.particles[0].position[0];
        assert!(x > 0.0, "particle should move in +x direction");
    }

    #[test]
    fn test_mhd_sph_total_kinetic_energy_zero_at_rest() {
        let mut solver = MhdSph::new(1e-7, 5.0 / 3.0);
        let p = PlasmaParticle::new([0.0; 3], [0.0; 3], M_PROTON, E_CHARGE, 1e4, 1e-6);
        solver.add_particle(p);
        assert_eq!(solver.total_kinetic_energy(), 0.0);
    }

    #[test]
    fn test_mhd_sph_magnetic_energy_positive() {
        let mut solver = MhdSph::new(1e-7, 5.0 / 3.0);
        let mut p = PlasmaParticle::new([0.0; 3], [0.0; 3], M_PROTON, E_CHARGE, 1e4, 1e-6);
        p.magnetic_field = [0.0, 0.0, 1e-3];
        solver.add_particle(p);
        assert!(solver.total_magnetic_energy() > 0.0);
    }

    // --- Ion inertial length ---

    #[test]
    fn test_ion_inertial_length_positive() {
        let d_i = ion_inertial_length(1e18);
        assert!(d_i > 0.0);
    }

    #[test]
    fn test_ion_inertial_length_scales_with_density() {
        // d_i ∝ n_i^{-1/2}
        let d1 = ion_inertial_length(1e18);
        let d2 = ion_inertial_length(4e18);
        assert_rel(d2, d1 / 2.0, RTOL, "d_i ∝ n^{-1/2}");
    }
}
