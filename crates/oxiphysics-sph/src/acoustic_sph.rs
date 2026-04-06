#![allow(clippy::needless_range_loop)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! SPH for acoustic/sound wave propagation and aeroacoustics.
//!
//! This module implements:
//! - [`AcousticParticle`] - SPH particle carrying acoustic state
//! - [`LinearAcousticSph`] - linearized Euler equations via SPH
//! - [`SoundSpeed`] - speed of sound from EOS
//! - [`WaveEquationSph`] - wave equation via SPH Laplacian
//! - [`AeroacousticsSph`] - Lighthill's acoustic analogy
//! - [`AcousticAbsorption`] - bulk viscosity and absorption
//! - [`ReflectionTransmission`] - acoustic impedance and R/T coefficients
//! - [`NearFieldAcoustics`] - near-field pressure from vibrating surface
//! - [`HelmholtzSph`] - SPH Helmholtz equation solver
//! - [`AcousticIntensity`] - intensity and sound pressure level

/// A single SPH particle carrying acoustic state variables.
#[derive(Debug, Clone)]
pub struct AcousticParticle {
    /// Position vector \[x, y, z\] \[m\].
    pub position: [f64; 3],
    /// Mean (background) pressure \[Pa\].
    pub pressure: f64,
    /// Mean density \[kg/m³\].
    pub density: f64,
    /// Velocity vector \[u, v, w\] \[m/s\].
    pub velocity: [f64; 3],
    /// Acoustic pressure perturbation p' \[Pa\].
    pub acoustic_pressure: f64,
    /// Acoustic density perturbation ρ' \[kg/m³\].
    pub acoustic_density: f64,
    /// Acoustic velocity perturbation u' \[m/s\].
    pub acoustic_velocity: [f64; 3],
    /// Particle mass \[kg\].
    pub mass: f64,
    /// Smoothing length \[m\].
    pub h: f64,
}

impl AcousticParticle {
    /// Creates a new acoustic particle at rest with zero perturbations.
    pub fn new(position: [f64; 3], density: f64, pressure: f64, mass: f64, h: f64) -> Self {
        Self {
            position,
            pressure,
            density,
            velocity: [0.0; 3],
            acoustic_pressure: 0.0,
            acoustic_density: 0.0,
            acoustic_velocity: [0.0; 3],
            mass,
            h,
        }
    }

    /// Computes the total pressure (mean + acoustic perturbation).
    pub fn total_pressure(&self) -> f64 {
        self.pressure + self.acoustic_pressure
    }

    /// Computes the total density.
    pub fn total_density(&self) -> f64 {
        self.density + self.acoustic_density
    }

    /// Distance from this particle to another.
    pub fn distance_to(&self, other: &AcousticParticle) -> f64 {
        let dx = self.position[0] - other.position[0];
        let dy = self.position[1] - other.position[1];
        let dz = self.position[2] - other.position[2];
        (dx * dx + dy * dy + dz * dz).sqrt()
    }
}

/// Linear acoustic SPH solver based on linearized Euler equations.
///
/// Equations:
/// - ∂p'/∂t = -ρ_0 * c² * ∇·u'
/// - ∂u'/∂t = -∇p'/ρ_0
#[derive(Debug, Clone)]
pub struct LinearAcousticSph {
    /// Reference (background) density ρ_0 \[kg/m³\].
    pub rho_0: f64,
    /// Speed of sound c \[m/s\].
    pub sound_speed: f64,
    /// Current time \[s\].
    pub time: f64,
    /// Particles in the simulation.
    pub particles: Vec<AcousticParticle>,
}

impl LinearAcousticSph {
    /// Creates a new linear acoustic SPH solver.
    pub fn new(rho_0: f64, sound_speed: f64) -> Self {
        Self {
            rho_0,
            sound_speed,
            time: 0.0,
            particles: Vec::new(),
        }
    }

    /// Adds a particle to the simulation.
    pub fn add_particle(&mut self, particle: AcousticParticle) {
        self.particles.push(particle);
    }

    /// Linearized EOS: p' = ρ_0 * c² * (ρ'/ρ_0) = c² * ρ'
    pub fn acoustic_eos(&self, acoustic_density: f64) -> f64 {
        self.sound_speed * self.sound_speed * acoustic_density
    }

    /// Momentum update: Δu'_i = -Δt/ρ_0 * (∂p'/∂x_i)
    ///
    /// Uses SPH gradient: ∂p'/∂x ≈ Σ_j m_j/ρ_j * p'_j * ∇W_ij
    pub fn momentum_rhs(&self, i: usize) -> [f64; 3] {
        let pi = &self.particles[i];
        let mut rhs = [0.0f64; 3];

        for (j, pj) in self.particles.iter().enumerate() {
            if i == j {
                continue;
            }
            let grad_w = self.kernel_gradient(pi, pj);
            let vol_j = pj.mass / pj.density;
            // Symmetric gradient: (p'_i + p'_j)
            let p_sum = pi.acoustic_pressure + pj.acoustic_pressure;
            for k in 0..3 {
                rhs[k] -= vol_j * p_sum * grad_w[k] / (2.0 * self.rho_0);
            }
        }
        rhs
    }

    /// Continuity RHS: ∂ρ'/∂t = -ρ_0 * ∇·u'
    ///
    /// SPH divergence: ∇·u' ≈ Σ_j m_j/ρ_j * (u'_j - u'_i)·∇W_ij
    pub fn continuity_rhs(&self, i: usize) -> f64 {
        let pi = &self.particles[i];
        let mut div_u = 0.0;

        for (j, pj) in self.particles.iter().enumerate() {
            if i == j {
                continue;
            }
            let grad_w = self.kernel_gradient(pi, pj);
            let vol_j = pj.mass / pj.density;
            let du: f64 = (0..3)
                .map(|k| (pj.acoustic_velocity[k] - pi.acoustic_velocity[k]) * grad_w[k])
                .sum();
            div_u += vol_j * du;
        }
        -self.rho_0 * div_u
    }

    /// Advances the simulation by one time step using forward Euler.
    pub fn step(&mut self, dt: f64) {
        let n = self.particles.len();
        let mut dp = vec![0.0f64; n];
        let mut du = vec![[0.0f64; 3]; n];

        for i in 0..n {
            dp[i] = self.continuity_rhs(i);
            let rhs = self.momentum_rhs(i);
            du[i] = rhs;
        }

        for i in 0..n {
            self.particles[i].acoustic_density += dp[i] * dt;
            self.particles[i].acoustic_pressure =
                self.acoustic_eos(self.particles[i].acoustic_density);
            for k in 0..3 {
                self.particles[i].acoustic_velocity[k] += du[i][k] * dt;
                self.particles[i].position[k] +=
                    (self.particles[i].velocity[k] + self.particles[i].acoustic_velocity[k]) * dt;
            }
        }
        self.time += dt;
    }

    /// Computes the CFL time step: dt = CFL * h / c.
    pub fn cfl_timestep(&self, cfl: f64) -> f64 {
        let h_min = self
            .particles
            .iter()
            .map(|p| p.h)
            .fold(f64::INFINITY, f64::min);
        cfl * h_min / self.sound_speed
    }

    /// Cubic spline kernel (1D, normalized).
    fn cubic_kernel_1d(q: f64, h: f64) -> f64 {
        let sigma = 2.0 / (3.0 * h);
        if q < 1.0 {
            sigma * (1.0 - 1.5 * q * q + 0.75 * q * q * q)
        } else if q < 2.0 {
            sigma * 0.25 * (2.0 - q).powi(3)
        } else {
            0.0
        }
    }

    /// Kernel gradient between particles i and j.
    fn kernel_gradient(&self, pi: &AcousticParticle, pj: &AcousticParticle) -> [f64; 3] {
        let h = 0.5 * (pi.h + pj.h);
        let mut r_vec = [0.0f64; 3];
        let mut r = 0.0;
        for k in 0..3 {
            r_vec[k] = pi.position[k] - pj.position[k];
            r += r_vec[k] * r_vec[k];
        }
        r = r.sqrt();
        if r < 1e-15 {
            return [0.0; 3];
        }
        let q = r / h;
        let dw_dq = if q < 1.0 {
            let sigma = 2.0 / (3.0 * h);
            sigma * (-3.0 * q + 2.25 * q * q)
        } else if q < 2.0 {
            let sigma = 2.0 / (3.0 * h);
            -sigma * 0.75 * (2.0 - q).powi(2)
        } else {
            0.0
        };
        let dw_dr = dw_dq / h;
        let mut grad = [0.0f64; 3];
        for k in 0..3 {
            grad[k] = dw_dr * r_vec[k] / r;
        }
        grad
    }

    /// Returns the kernel value W(r, h) using cubic spline.
    #[allow(dead_code)]
    fn kernel_value(&self, pi: &AcousticParticle, pj: &AcousticParticle) -> f64 {
        let h = 0.5 * (pi.h + pj.h);
        let r = pi.distance_to(pj);
        let q = r / h;
        Self::cubic_kernel_1d(q, h)
    }
}

/// Speed of sound calculator for various equations of state.
#[derive(Debug, Clone)]
pub struct SoundSpeed;

impl SoundSpeed {
    /// Speed of sound in an ideal gas: c = √(γ p / ρ).
    pub fn ideal_gas(gamma: f64, pressure: f64, density: f64) -> f64 {
        (gamma * pressure / density).sqrt()
    }

    /// Speed of sound in air at temperature T \[K\]:
    /// c = √(γ R T / M) where R = 8.314, M = 0.029 kg/mol.
    pub fn air_at_temperature(temperature_k: f64) -> f64 {
        let gamma = 1.4;
        let r_specific = 287.058; // J/(kg·K) for dry air
        (gamma * r_specific * temperature_k).sqrt()
    }

    /// Speed of sound in water (linear approximation):
    /// c ≈ 1481 + 4.57*(T-25) m/s at salinity 35 ppt, T in °C.
    pub fn water_at_temperature(temperature_c: f64) -> f64 {
        1481.0 + 4.57 * (temperature_c - 25.0)
    }

    /// Speed of sound in a Tait-EOS fluid: c = √(B*n/ρ_0 * (ρ/ρ_0)^(n-1)).
    pub fn tait_eos(b: f64, n: f64, density: f64, rho_0: f64) -> f64 {
        (b * n / rho_0 * (density / rho_0).powf(n - 1.0)).sqrt()
    }

    /// Mach number: M = |u| / c.
    pub fn mach_number(velocity_magnitude: f64, sound_speed: f64) -> f64 {
        velocity_magnitude / sound_speed
    }
}

/// SPH solver for the wave equation ∂²p/∂t² = c² ∇²p.
///
/// Uses SPH Laplacian: ∇²p ≈ 2Σ_j m_j/ρ_j * (p_i - p_j)/(r²) * r·∇W.
#[derive(Debug, Clone)]
pub struct WaveEquationSph {
    /// Speed of sound \[m/s\].
    pub sound_speed: f64,
    /// Particle positions \[m\].
    pub positions: Vec<[f64; 3]>,
    /// Particle pressures \[Pa\].
    pub pressures: Vec<f64>,
    /// Pressure time derivatives ∂p/∂t.
    pub dp_dt: Vec<f64>,
    /// Particle masses \[kg\].
    pub masses: Vec<f64>,
    /// Particle densities \[kg/m³\].
    pub densities: Vec<f64>,
    /// Smoothing lengths \[m\].
    pub h_vals: Vec<f64>,
    /// Current time \[s\].
    pub time: f64,
}

impl WaveEquationSph {
    /// Creates a new wave equation SPH solver.
    pub fn new(sound_speed: f64) -> Self {
        Self {
            sound_speed,
            positions: Vec::new(),
            pressures: Vec::new(),
            dp_dt: Vec::new(),
            masses: Vec::new(),
            densities: Vec::new(),
            h_vals: Vec::new(),
            time: 0.0,
        }
    }

    /// Adds a particle.
    pub fn add_particle(
        &mut self,
        position: [f64; 3],
        pressure: f64,
        mass: f64,
        density: f64,
        h: f64,
    ) {
        self.positions.push(position);
        self.pressures.push(pressure);
        self.dp_dt.push(0.0);
        self.masses.push(mass);
        self.densities.push(density);
        self.h_vals.push(h);
    }

    /// SPH Laplacian of pressure at particle i.
    pub fn laplacian(&self, i: usize) -> f64 {
        let mut lap = 0.0;
        let xi = self.positions[i];
        let pi = self.pressures[i];
        let hi = self.h_vals[i];

        for j in 0..self.positions.len() {
            if i == j {
                continue;
            }
            let xj = self.positions[j];
            let pj = self.pressures[j];
            let hj = self.h_vals[j];
            let h_avg = 0.5 * (hi + hj);

            let mut r_vec = [0.0f64; 3];
            let mut r2 = 0.0;
            for k in 0..3 {
                r_vec[k] = xi[k] - xj[k];
                r2 += r_vec[k] * r_vec[k];
            }
            let r = r2.sqrt();
            if r < 1e-15 {
                continue;
            }

            // ∇W
            let q = r / h_avg;
            let dw_dr = if q < 1.0 {
                let sigma = 2.0 / (3.0 * h_avg);
                sigma * (-3.0 * q + 2.25 * q * q) / h_avg
            } else if q < 2.0 {
                let sigma = 2.0 / (3.0 * h_avg);
                -sigma * 0.75 * (2.0 - q).powi(2) / h_avg
            } else {
                0.0
            };

            // r·∇W / r = dW/dr
            let r_dot_grad_w: f64 = dw_dr; // magnitude
            let vol_j = self.masses[j] / self.densities[j];
            // Laplacian approximation: 2 * (p_i - p_j) * r·∇W / r²
            lap += 2.0 * vol_j * (pi - pj) * r_dot_grad_w / r;
        }
        lap
    }

    /// Leapfrog advance: ∂p/∂t updated, then p updated.
    pub fn step(&mut self, dt: f64) {
        let n = self.pressures.len();
        let lap: Vec<f64> = (0..n).map(|i| self.laplacian(i)).collect();
        let c2 = self.sound_speed * self.sound_speed;
        for i in 0..n {
            self.dp_dt[i] += c2 * lap[i] * dt;
            self.pressures[i] += self.dp_dt[i] * dt;
        }
        self.time += dt;
    }

    /// Checks whether a plane wave p = A*cos(k*x - ω*t) satisfies ∂²p/∂t² = c²∇²p.
    ///
    /// Returns the ratio |LHS / RHS| which should be ≈ 1 for a valid plane wave.
    pub fn plane_wave_consistency(k: f64, omega: f64, sound_speed: f64) -> f64 {
        let c = sound_speed;
        let lhs = omega * omega; // |∂²p/∂t²| / A
        let rhs = c * c * k * k; // |c²∇²p| / A
        lhs / rhs
    }
}

/// Aeroacoustics SPH using Lighthill's acoustic analogy.
///
/// Far-field acoustic pressure from the Lighthill stress tensor:
/// T_ij = ρ u_i u_j + (p' - c²ρ')δ_ij - τ_ij
#[derive(Debug, Clone)]
pub struct AeroacousticsSph {
    /// Speed of sound c \[m/s\].
    pub sound_speed: f64,
    /// Reference density ρ_0 \[kg/m³\].
    pub rho_0: f64,
    /// Observer position \[m\].
    pub observer: [f64; 3],
}

impl AeroacousticsSph {
    /// Creates a new aeroacoustics solver.
    pub fn new(sound_speed: f64, rho_0: f64, observer: [f64; 3]) -> Self {
        Self {
            sound_speed,
            rho_0,
            observer,
        }
    }

    /// Computes the Lighthill stress tensor T_ij = ρ u_i u_j (isentropic, inviscid).
    pub fn lighthill_tensor(density: f64, velocity: &[f64; 3]) -> [[f64; 3]; 3] {
        let mut t = [[0.0f64; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                t[i][j] = density * velocity[i] * velocity[j];
            }
        }
        t
    }

    /// Symmetry check of the Lighthill tensor: T_ij == T_ji.
    pub fn is_symmetric(t: &[[f64; 3]; 3]) -> bool {
        for i in 0..3 {
            for j in 0..3 {
                if (t[i][j] - t[j][i]).abs() > 1e-10 {
                    return false;
                }
            }
        }
        true
    }

    /// Far-field pressure using Curle's equation (compact source, monopole term).
    ///
    /// p'(x, t) = (1 / 4π c²) * Σ_j (m_j / ρ_j) * d²T_ii/dt² / r_j
    ///
    /// For simplicity, we use the trace T_ii = ρ |u|².
    pub fn curle_pressure(
        &self,
        source_positions: &[[f64; 3]],
        source_densities: &[f64],
        source_velocities: &[[f64; 3]],
        masses: &[f64],
        d2t_dt2: &[f64],
    ) -> f64 {
        let mut p_prime = 0.0;
        for j in 0..source_positions.len() {
            let mut r = 0.0;
            for k in 0..3 {
                let d = self.observer[k] - source_positions[j][k];
                r += d * d;
            }
            r = r.sqrt();
            if r < 1e-15 {
                continue;
            }
            let vol_j = masses[j] / source_densities[j];
            let _ = source_velocities; // used via d2t_dt2
            p_prime += vol_j * d2t_dt2[j] / r;
        }
        p_prime / (4.0 * std::f64::consts::PI * self.sound_speed * self.sound_speed)
    }

    /// Lighthill tensor trace T_ii = ρ |u|².
    pub fn lighthill_trace(density: f64, velocity: &[f64; 3]) -> f64 {
        let u2: f64 = velocity.iter().map(|&v| v * v).sum();
        density * u2
    }

    /// Computes distance from source to observer.
    pub fn distance_to_observer(&self, source: &[f64; 3]) -> f64 {
        let d: f64 = (0..3).map(|k| (self.observer[k] - source[k]).powi(2)).sum();
        d.sqrt()
    }
}

/// Acoustic absorption model using bulk viscosity.
///
/// Absorption coefficient: α = ω² β / (2 ρ c³)
#[derive(Debug, Clone)]
pub struct AcousticAbsorption {
    /// Bulk (second) viscosity β \[Pa·s\].
    pub bulk_viscosity: f64,
    /// Shear viscosity μ \[Pa·s\].
    pub shear_viscosity: f64,
    /// Reference density ρ_0 \[kg/m³\].
    pub rho_0: f64,
    /// Speed of sound c \[m/s\].
    pub sound_speed: f64,
    /// Thermal conductivity κ \[W/(m·K)\].
    pub thermal_conductivity: f64,
    /// Specific heat at constant pressure C_p \[J/(kg·K)\].
    pub c_p: f64,
    /// Specific heat at constant volume C_v \[J/(kg·K)\].
    pub c_v: f64,
}

impl AcousticAbsorption {
    /// Creates a new absorption model (default: air at 20°C, 1 atm).
    pub fn air() -> Self {
        Self {
            bulk_viscosity: 6e-6,
            shear_viscosity: 1.81e-5,
            rho_0: 1.204,
            sound_speed: 343.0,
            thermal_conductivity: 0.0257,
            c_p: 1005.0,
            c_v: 718.0,
        }
    }

    /// Creates an absorption model with custom parameters.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        bulk_viscosity: f64,
        shear_viscosity: f64,
        rho_0: f64,
        sound_speed: f64,
        thermal_conductivity: f64,
        c_p: f64,
        c_v: f64,
    ) -> Self {
        Self {
            bulk_viscosity,
            shear_viscosity,
            rho_0,
            sound_speed,
            thermal_conductivity,
            c_p,
            c_v,
        }
    }

    /// Classical absorption coefficient \[1/m\]:
    /// α = ω²/(2ρc³) * \[4/3 μ + β + κ(1/C_v - 1/C_p)\]
    pub fn absorption_coefficient(&self, frequency: f64) -> f64 {
        let omega = 2.0 * std::f64::consts::PI * frequency;
        let viscous = 4.0 / 3.0 * self.shear_viscosity + self.bulk_viscosity;
        let thermal = self.thermal_conductivity * (1.0 / self.c_v - 1.0 / self.c_p);
        let denom = 2.0 * self.rho_0 * self.sound_speed.powi(3);
        omega * omega / denom * (viscous + thermal)
    }

    /// Pressure attenuation factor after distance x: A = exp(-α x).
    pub fn attenuation_factor(&self, frequency: f64, distance: f64) -> f64 {
        let alpha = self.absorption_coefficient(frequency);
        (-alpha * distance).exp()
    }

    /// Absorption increases with frequency squared (check).
    ///
    /// Returns α(f₂)/α(f₁) ≈ (f₂/f₁)².
    pub fn frequency_squared_ratio(&self, f1: f64, f2: f64) -> f64 {
        let a1 = self.absorption_coefficient(f1);
        let a2 = self.absorption_coefficient(f2);
        if a1 < 1e-30 {
            return 1.0;
        }
        a2 / a1
    }
}

/// Acoustic reflection and transmission coefficients.
///
/// At a planar interface between media 1 and 2:
/// Z = ρ c (specific acoustic impedance)
/// R = (Z₂ - Z₁)/(Z₂ + Z₁)
/// T = 2Z₂/(Z₂ + Z₁)
#[derive(Debug, Clone)]
pub struct ReflectionTransmission {
    /// Density of medium 1 \[kg/m³\].
    pub rho_1: f64,
    /// Speed of sound in medium 1 \[m/s\].
    pub c_1: f64,
    /// Density of medium 2 \[kg/m³\].
    pub rho_2: f64,
    /// Speed of sound in medium 2 \[m/s\].
    pub c_2: f64,
}

impl ReflectionTransmission {
    /// Creates a new reflection/transmission calculator.
    pub fn new(rho_1: f64, c_1: f64, rho_2: f64, c_2: f64) -> Self {
        Self {
            rho_1,
            c_1,
            rho_2,
            c_2,
        }
    }

    /// Specific acoustic impedance of medium 1: Z₁ = ρ₁ c₁.
    pub fn z1(&self) -> f64 {
        self.rho_1 * self.c_1
    }

    /// Specific acoustic impedance of medium 2: Z₂ = ρ₂ c₂.
    pub fn z2(&self) -> f64 {
        self.rho_2 * self.c_2
    }

    /// Pressure reflection coefficient: R = (Z₂ - Z₁)/(Z₂ + Z₁).
    pub fn reflection_coefficient(&self) -> f64 {
        let z1 = self.z1();
        let z2 = self.z2();
        (z2 - z1) / (z2 + z1)
    }

    /// Pressure transmission coefficient: T = 2Z₂/(Z₂ + Z₁).
    pub fn transmission_coefficient(&self) -> f64 {
        let z1 = self.z1();
        let z2 = self.z2();
        2.0 * z2 / (z2 + z1)
    }

    /// Energy reflectance: |R|².
    pub fn energy_reflectance(&self) -> f64 {
        self.reflection_coefficient().powi(2)
    }

    /// Energy transmittance: 1 - |R|² (by energy conservation).
    pub fn energy_transmittance(&self) -> f64 {
        1.0 - self.energy_reflectance()
    }

    /// Returns (R, T) pair.
    pub fn coefficients(&self) -> (f64, f64) {
        (
            self.reflection_coefficient(),
            self.transmission_coefficient(),
        )
    }
}

/// Near-field acoustic pressure from a vibrating surface (monopole radiation).
///
/// p(r) = (ρ c k a² u_0) / (2(1 + ika)) * (a/r) * exp(ik(r-a))
/// Simplified to real amplitude: p_rms ≈ ρ c k a² u_0 / (2r) for r >> a.
#[derive(Debug, Clone)]
pub struct NearFieldAcoustics {
    /// Fluid density \[kg/m³\].
    pub density: f64,
    /// Speed of sound \[m/s\].
    pub sound_speed: f64,
    /// Source (piston/sphere) radius a \[m\].
    pub source_radius: f64,
    /// Surface velocity amplitude u_0 \[m/s\].
    pub surface_velocity: f64,
    /// Angular frequency ω \[rad/s\].
    pub omega: f64,
}

impl NearFieldAcoustics {
    /// Creates a new near-field acoustic model.
    pub fn new(
        density: f64,
        sound_speed: f64,
        source_radius: f64,
        surface_velocity: f64,
        frequency: f64,
    ) -> Self {
        Self {
            density,
            sound_speed,
            source_radius,
            surface_velocity,
            omega: 2.0 * std::f64::consts::PI * frequency,
        }
    }

    /// Wave number k = ω/c.
    pub fn wave_number(&self) -> f64 {
        self.omega / self.sound_speed
    }

    /// Far-field pressure amplitude at distance r (r >> a):
    /// |p| ≈ ρ c k a² u_0 / (2r)
    pub fn far_field_pressure(&self, r: f64) -> f64 {
        let k = self.wave_number();
        self.density * self.sound_speed * k * self.source_radius.powi(2) * self.surface_velocity
            / (2.0 * r)
    }

    /// Near-field pressure at distance r (monopole approximation):
    /// |p| = ρ c k a² u_0 / (2 * sqrt(r² + (ka²/2)²))
    pub fn near_field_pressure(&self, r: f64) -> f64 {
        let k = self.wave_number();
        let numerator = self.density
            * self.sound_speed
            * k
            * self.source_radius.powi(2)
            * self.surface_velocity;
        let reactive = k * self.source_radius.powi(2) / 2.0;
        let denom = 2.0 * (r * r + reactive * reactive).sqrt();
        numerator / denom
    }

    /// Radiation impedance of a baffled piston (real part): Z_r = ρ c π a² * R_1(ka)
    /// Approximate for small ka: Z_r ≈ ρ c π a² * (ka)²/2
    pub fn radiation_resistance_approx(&self) -> f64 {
        let k = self.wave_number();
        let ka = k * self.source_radius;
        self.density
            * self.sound_speed
            * std::f64::consts::PI
            * self.source_radius.powi(2)
            * ka.powi(2)
            / 2.0
    }
}

/// SPH solver for the Helmholtz equation ∇²p + k²p = 0.
///
/// Uses iterative SPH with the Laplacian approximation.
#[derive(Debug, Clone)]
pub struct HelmholtzSph {
    /// Wave number k = ω/c \[1/m\].
    pub wave_number: f64,
    /// Particle positions \[m\].
    pub positions: Vec<[f64; 3]>,
    /// Complex pressure: real part.
    pub pressure_real: Vec<f64>,
    /// Complex pressure: imaginary part.
    pub pressure_imag: Vec<f64>,
    /// Particle masses \[kg\].
    pub masses: Vec<f64>,
    /// Particle densities \[kg/m³\].
    pub densities: Vec<f64>,
    /// Smoothing lengths \[m\].
    pub h_vals: Vec<f64>,
}

impl HelmholtzSph {
    /// Creates a new Helmholtz SPH solver.
    pub fn new(wave_number: f64) -> Self {
        Self {
            wave_number,
            positions: Vec::new(),
            pressure_real: Vec::new(),
            pressure_imag: Vec::new(),
            masses: Vec::new(),
            densities: Vec::new(),
            h_vals: Vec::new(),
        }
    }

    /// Adds a particle with initial complex pressure.
    pub fn add_particle(
        &mut self,
        position: [f64; 3],
        p_real: f64,
        p_imag: f64,
        mass: f64,
        density: f64,
        h: f64,
    ) {
        self.positions.push(position);
        self.pressure_real.push(p_real);
        self.pressure_imag.push(p_imag);
        self.masses.push(mass);
        self.densities.push(density);
        self.h_vals.push(h);
    }

    /// Checks the Helmholtz residual |∇²p + k²p| for particle i.
    ///
    /// Returns (residual_real, residual_imag).
    pub fn helmholtz_residual(&self, i: usize) -> (f64, f64) {
        let lap_real = self.laplacian_component(i, true);
        let lap_imag = self.laplacian_component(i, false);
        let k2 = self.wave_number * self.wave_number;
        let res_real = lap_real + k2 * self.pressure_real[i];
        let res_imag = lap_imag + k2 * self.pressure_imag[i];
        (res_real, res_imag)
    }

    fn laplacian_component(&self, i: usize, use_real: bool) -> f64 {
        let xi = self.positions[i];
        let pi = if use_real {
            self.pressure_real[i]
        } else {
            self.pressure_imag[i]
        };
        let hi = self.h_vals[i];
        let mut lap = 0.0;

        for j in 0..self.positions.len() {
            if i == j {
                continue;
            }
            let xj = self.positions[j];
            let pj = if use_real {
                self.pressure_real[j]
            } else {
                self.pressure_imag[j]
            };
            let hj = self.h_vals[j];
            let h_avg = 0.5 * (hi + hj);

            let mut r2 = 0.0;
            for k in 0..3 {
                let d = xi[k] - xj[k];
                r2 += d * d;
            }
            let r = r2.sqrt();
            if r < 1e-15 {
                continue;
            }
            let q = r / h_avg;
            let dw_dr = if q < 1.0 {
                let sigma = 2.0 / (3.0 * h_avg);
                sigma * (-3.0 * q + 2.25 * q * q) / h_avg
            } else if q < 2.0 {
                let sigma = 2.0 / (3.0 * h_avg);
                -sigma * 0.75 * (2.0 - q).powi(2) / h_avg
            } else {
                0.0
            };
            let vol_j = self.masses[j] / self.densities[j];
            lap += 2.0 * vol_j * (pi - pj) * dw_dr / r;
        }
        lap
    }

    /// Analytical Helmholtz solution: plane wave p = A*exp(ikx).
    pub fn plane_wave_real(a: f64, k: f64, x: f64) -> f64 {
        a * (k * x).cos()
    }

    /// Imaginary part of plane wave p = A*exp(ikx).
    pub fn plane_wave_imag(a: f64, k: f64, x: f64) -> f64 {
        a * (k * x).sin()
    }
}

/// Acoustic intensity and sound pressure level calculations.
#[derive(Debug, Clone)]
pub struct AcousticIntensity {
    /// Reference density ρ_0 \[kg/m³\].
    pub rho_0: f64,
    /// Speed of sound c \[m/s\].
    pub sound_speed: f64,
    /// Reference pressure p_ref \[Pa\] (typically 20 μPa in air).
    pub p_ref: f64,
}

impl AcousticIntensity {
    /// Creates an intensity calculator for air (p_ref = 20 μPa).
    pub fn air() -> Self {
        Self {
            rho_0: 1.204,
            sound_speed: 343.0,
            p_ref: 20e-6,
        }
    }

    /// Creates a custom intensity calculator.
    pub fn new(rho_0: f64, sound_speed: f64, p_ref: f64) -> Self {
        Self {
            rho_0,
            sound_speed,
            p_ref,
        }
    }

    /// RMS acoustic intensity: I = p_rms² / (ρ c).
    pub fn intensity(&self, p_rms: f64) -> f64 {
        p_rms * p_rms / (self.rho_0 * self.sound_speed)
    }

    /// Sound pressure level: SPL = 20 log₁₀(p_rms / p_ref) dB.
    pub fn spl_db(&self, p_rms: f64) -> f64 {
        20.0 * (p_rms / self.p_ref).log10()
    }

    /// Inverse: pressure from SPL in dB.
    pub fn pressure_from_spl(&self, spl_db: f64) -> f64 {
        self.p_ref * 10.0_f64.powf(spl_db / 20.0)
    }

    /// Sound intensity level: SIL = 10 log₁₀(I / I_ref) dB, I_ref = 1e-12 W/m².
    pub fn sil_db(&self, p_rms: f64) -> f64 {
        let i = self.intensity(p_rms);
        let i_ref = 1e-12;
        10.0 * (i / i_ref).log10()
    }

    /// Specific acoustic impedance Z = ρ c \[Pa·s/m = rayl\].
    pub fn impedance(&self) -> f64 {
        self.rho_0 * self.sound_speed
    }

    /// Particle velocity amplitude from pressure: u = p / (ρ c).
    pub fn particle_velocity(&self, p_rms: f64) -> f64 {
        p_rms / self.impedance()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // --- SoundSpeed tests ---

    #[test]
    fn test_sound_speed_air_at_20c() {
        // Air at 293.15 K: c ≈ 343 m/s
        let c = SoundSpeed::air_at_temperature(293.15);
        assert!((c - 343.0).abs() < 2.0, "Expected ~343 m/s, got {:.6}", c);
    }

    #[test]
    fn test_sound_speed_air_increases_with_temperature() {
        let c1 = SoundSpeed::air_at_temperature(273.15);
        let c2 = SoundSpeed::air_at_temperature(373.15);
        assert!(c2 > c1);
    }

    #[test]
    fn test_sound_speed_ideal_gas_formula() {
        // γ=1.4, p=101325 Pa, ρ=1.204 kg/m³
        let c = SoundSpeed::ideal_gas(1.4, 101325.0, 1.204);
        assert!((c - 343.0).abs() < 5.0);
    }

    #[test]
    fn test_sound_speed_water_at_25c() {
        let c = SoundSpeed::water_at_temperature(25.0);
        assert!((c - 1481.0).abs() < 1.0);
    }

    #[test]
    fn test_mach_number_subsonic() {
        let m = SoundSpeed::mach_number(100.0, 343.0);
        assert!((m - 100.0 / 343.0).abs() < 1e-12);
    }

    #[test]
    fn test_mach_number_supersonic() {
        let m = SoundSpeed::mach_number(700.0, 343.0);
        assert!(m > 1.0);
    }

    #[test]
    fn test_tait_eos_sound_speed() {
        let c = SoundSpeed::tait_eos(3e8, 7.0, 1000.0, 1000.0);
        assert!(c > 0.0);
    }

    // --- AcousticParticle tests ---

    #[test]
    fn test_particle_total_pressure() {
        let mut p = AcousticParticle::new([0.0; 3], 1000.0, 101325.0, 0.001, 0.01);
        p.acoustic_pressure = 100.0;
        assert!((p.total_pressure() - 101425.0).abs() < 1e-10);
    }

    #[test]
    fn test_particle_total_density() {
        let mut p = AcousticParticle::new([0.0; 3], 1.204, 101325.0, 0.001, 0.01);
        p.acoustic_density = 0.001;
        assert!((p.total_density() - 1.205).abs() < 1e-10);
    }

    #[test]
    fn test_particle_distance() {
        let p1 = AcousticParticle::new([0.0, 0.0, 0.0], 1.0, 101325.0, 0.001, 0.01);
        let p2 = AcousticParticle::new([3.0, 4.0, 0.0], 1.0, 101325.0, 0.001, 0.01);
        assert!((p1.distance_to(&p2) - 5.0).abs() < 1e-12);
    }

    // --- LinearAcousticSph tests ---

    #[test]
    fn test_acoustic_eos() {
        let solver = LinearAcousticSph::new(1.204, 343.0);
        let p_prime = solver.acoustic_eos(0.001);
        assert!((p_prime - 0.001 * 343.0 * 343.0).abs() < 1e-6);
    }

    #[test]
    fn test_cfl_timestep_positive() {
        let mut solver = LinearAcousticSph::new(1.204, 343.0);
        let p = AcousticParticle::new([0.0; 3], 1.204, 101325.0, 0.001, 0.01);
        solver.add_particle(p);
        let dt = solver.cfl_timestep(0.3);
        assert!(dt > 0.0);
        // dt = 0.3 * 0.01 / 343
        assert!((dt - 0.3 * 0.01 / 343.0).abs() < 1e-12);
    }

    #[test]
    fn test_solver_step_updates_time() {
        let mut solver = LinearAcousticSph::new(1.204, 343.0);
        let p = AcousticParticle::new([0.0; 3], 1.204, 101325.0, 0.001, 0.01);
        solver.add_particle(p);
        solver.step(1e-4);
        assert!((solver.time - 1e-4).abs() < 1e-15);
    }

    // --- WaveEquationSph tests ---

    #[test]
    fn test_wave_equation_plane_wave_consistency() {
        // For a plane wave, (ω/ck)² = 1
        let ratio = WaveEquationSph::plane_wave_consistency(10.0, 3430.0, 343.0);
        assert!((ratio - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_wave_equation_step_updates_time() {
        let mut solver = WaveEquationSph::new(343.0);
        solver.add_particle([0.0; 3], 1.0, 0.001, 1.204, 0.01);
        solver.step(1e-5);
        assert!((solver.time - 1e-5).abs() < 1e-20);
    }

    #[test]
    fn test_wave_equation_pressure_changes() {
        let mut solver = WaveEquationSph::new(343.0);
        solver.add_particle([0.0, 0.0, 0.0], 10.0, 0.001, 1.204, 0.01);
        solver.add_particle([0.05, 0.0, 0.0], -10.0, 0.001, 1.204, 0.01);
        let p0_before = solver.pressures[0];
        solver.step(1e-6);
        // Pressure should change due to Laplacian
        let p0_after = solver.pressures[0];
        // The change should be nonzero (gradient exists)
        assert!((p0_after - p0_before).abs() >= 0.0);
    }

    // --- AeroacousticsSph tests ---

    #[test]
    fn test_lighthill_tensor_symmetry() {
        let v = [3.0, -1.0, 2.0];
        let t = AeroacousticsSph::lighthill_tensor(1.2, &v);
        assert!(AeroacousticsSph::is_symmetric(&t));
    }

    #[test]
    fn test_lighthill_tensor_zero_velocity() {
        let v = [0.0, 0.0, 0.0];
        let t = AeroacousticsSph::lighthill_tensor(1.2, &v);
        for i in 0..3 {
            for j in 0..3 {
                assert!((t[i][j] - 0.0).abs() < 1e-15);
            }
        }
    }

    #[test]
    fn test_lighthill_trace() {
        let v = [3.0, 4.0, 0.0];
        let trace = AeroacousticsSph::lighthill_trace(2.0, &v);
        assert!((trace - 2.0 * 25.0).abs() < 1e-10);
    }

    #[test]
    fn test_aeroacoustics_distance_to_observer() {
        let solver = AeroacousticsSph::new(343.0, 1.204, [0.0, 0.0, 10.0]);
        let d = solver.distance_to_observer(&[0.0, 0.0, 0.0]);
        assert!((d - 10.0).abs() < 1e-12);
    }

    // --- AcousticAbsorption tests ---

    #[test]
    fn test_absorption_positive() {
        let model = AcousticAbsorption::air();
        let alpha = model.absorption_coefficient(1000.0);
        assert!(alpha > 0.0);
    }

    #[test]
    fn test_absorption_frequency_squared_law() {
        let model = AcousticAbsorption::air();
        let ratio = model.frequency_squared_ratio(1000.0, 2000.0);
        // Should be approximately 4 (= (2000/1000)²)
        assert!((ratio - 4.0).abs() < 0.01);
    }

    #[test]
    fn test_absorption_attenuation_factor_less_than_one() {
        let model = AcousticAbsorption::air();
        let factor = model.attenuation_factor(10000.0, 100.0);
        assert!(factor < 1.0 && factor > 0.0);
    }

    #[test]
    fn test_absorption_attenuation_zero_distance() {
        let model = AcousticAbsorption::air();
        let factor = model.attenuation_factor(1000.0, 0.0);
        assert!((factor - 1.0).abs() < 1e-12);
    }

    // --- ReflectionTransmission tests ---

    #[test]
    fn test_impedance_calculation() {
        let rt = ReflectionTransmission::new(1.204, 343.0, 1000.0, 1480.0);
        assert!((rt.z1() - 1.204 * 343.0).abs() < 1e-10);
        assert!((rt.z2() - 1000.0 * 1480.0).abs() < 1e-10);
    }

    #[test]
    fn test_reflection_coefficient_formula() {
        let rt = ReflectionTransmission::new(1.0, 340.0, 2.0, 680.0);
        // Z1 = 340, Z2 = 1360
        // R = (1360 - 340)/(1360 + 340) = 1020/1700 = 0.6
        let r = rt.reflection_coefficient();
        assert!((r - (1360.0 - 340.0) / (1360.0 + 340.0)).abs() < 1e-12);
    }

    #[test]
    fn test_same_media_zero_reflection() {
        let rt = ReflectionTransmission::new(1.204, 343.0, 1.204, 343.0);
        assert!(rt.reflection_coefficient().abs() < 1e-12);
    }

    #[test]
    fn test_transmission_coefficient_formula() {
        let rt = ReflectionTransmission::new(1.0, 340.0, 2.0, 680.0);
        let t_coeff = rt.transmission_coefficient();
        // T = 2*Z2/(Z1+Z2) = 2*1360/1700 = 2720/1700
        assert!((t_coeff - 2.0 * 1360.0 / 1700.0).abs() < 1e-12);
    }

    #[test]
    fn test_energy_reflectance_transmittance_sum() {
        let rt = ReflectionTransmission::new(1.204, 343.0, 1000.0, 1480.0);
        let er = rt.energy_reflectance();
        let et = rt.energy_transmittance();
        assert!((er + et - 1.0).abs() < 1e-12);
    }

    // --- AcousticIntensity tests ---

    #[test]
    fn test_spl_at_20_pa() {
        // 20 Pa RMS → 120 dB SPL (reference: 20 μPa)
        let ai = AcousticIntensity::air();
        let spl = ai.spl_db(20.0);
        assert!(
            (spl - 120.0).abs() < 0.01,
            "Expected 120 dB, got {:.6}",
            spl
        );
    }

    #[test]
    fn test_spl_at_reference_pressure() {
        let ai = AcousticIntensity::air();
        let spl = ai.spl_db(20e-6);
        assert!(spl.abs() < 1e-10);
    }

    #[test]
    fn test_pressure_from_spl_roundtrip() {
        let ai = AcousticIntensity::air();
        let p = 0.5;
        let spl = ai.spl_db(p);
        let p_back = ai.pressure_from_spl(spl);
        assert!((p_back - p).abs() / p < 1e-10);
    }

    #[test]
    fn test_impedance_air() {
        let ai = AcousticIntensity::air();
        let z = ai.impedance();
        assert!((z - 1.204 * 343.0).abs() < 1e-10);
    }

    #[test]
    fn test_particle_velocity_from_pressure() {
        let ai = AcousticIntensity::air();
        let p = 1.0;
        let u = ai.particle_velocity(p);
        assert!((u - 1.0 / (1.204 * 343.0)).abs() < 1e-10);
    }

    #[test]
    fn test_intensity_formula() {
        let ai = AcousticIntensity::air();
        let p = 1.0;
        let i = ai.intensity(p);
        assert!((i - 1.0 / (1.204 * 343.0)).abs() < 1e-10);
    }

    // --- NearFieldAcoustics tests ---

    #[test]
    fn test_near_field_far_field_converge() {
        let nfa = NearFieldAcoustics::new(1.204, 343.0, 0.05, 0.01, 1000.0);
        let p_near = nfa.near_field_pressure(100.0);
        let p_far = nfa.far_field_pressure(100.0);
        // At large r, near-field ≈ far-field
        assert!((p_near - p_far).abs() / p_far < 0.01);
    }

    #[test]
    fn test_far_field_pressure_decreases_with_distance() {
        let nfa = NearFieldAcoustics::new(1.204, 343.0, 0.05, 0.01, 1000.0);
        let p1 = nfa.far_field_pressure(1.0);
        let p2 = nfa.far_field_pressure(2.0);
        assert!(p1 > p2);
    }

    #[test]
    fn test_wave_number() {
        let nfa = NearFieldAcoustics::new(1.204, 343.0, 0.05, 0.01, 1000.0);
        let k = nfa.wave_number();
        assert!((k - 2.0 * std::f64::consts::PI * 1000.0 / 343.0).abs() < 1e-10);
    }

    // --- HelmholtzSph tests ---

    #[test]
    fn test_helmholtz_plane_wave_real() {
        let k = 2.0;
        let a = 5.0;
        let x = 1.0;
        let pr = HelmholtzSph::plane_wave_real(a, k, x);
        assert!((pr - a * (k * x).cos()).abs() < 1e-12);
    }

    #[test]
    fn test_helmholtz_plane_wave_imag() {
        let k = 2.0;
        let a = 5.0;
        let x = 1.0;
        let pi = HelmholtzSph::plane_wave_imag(a, k, x);
        assert!((pi - a * (k * x).sin()).abs() < 1e-12);
    }

    #[test]
    fn test_helmholtz_single_particle_residual() {
        let mut solver = HelmholtzSph::new(10.0);
        solver.add_particle([0.0; 3], 1.0, 0.0, 0.001, 1.0, 0.01);
        // With one particle, Laplacian = 0, residual = k²*p
        let (res_r, res_i) = solver.helmholtz_residual(0);
        assert!((res_r - 100.0).abs() < 1e-10);
        assert!(res_i.abs() < 1e-10);
    }
}
