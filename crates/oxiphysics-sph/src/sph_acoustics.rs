#![allow(clippy::needless_range_loop)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Acoustic SPH for wave propagation.
//!
//! Implements the linearised acoustic equations in a SPH framework, including:
//!
//! - [`AcousticParticle`]: SPH particle carrying acoustic state variables
//! - [`AcousticSPH`]: Time-integration driver for acoustic SPH
//! - [`PressureWave`]: Descriptor for a harmonic pressure wave
//! - [`AbsorbingBoundary`]: Perfectly Matched Layer (PML) absorber
//! - [`AcousticSource`]: Point source (monopole / dipole / quadrupole)
//! - [`sound_speed_water`]: Empirical sound speed in sea water
//! - [`acoustic_impedance`]: Characteristic acoustic impedance Z = ρ c
//! - [`wavelength`]: Acoustic wavelength λ = c / f
//! - [`cfl_acoustic`]: CFL time-step limit for acoustic SPH

#![allow(dead_code)]

use std::f64::consts::PI;

// ============================================================================
// Math helpers (private)
// ============================================================================

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
fn len3(v: [f64; 3]) -> f64 {
    dot3(v, v).sqrt()
}

// ============================================================================
// SPH kernel (cubic spline)
// ============================================================================

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

// ============================================================================
// AcousticParticle
// ============================================================================

/// A single SPH particle carrying acoustic state variables.
pub struct AcousticParticle {
    /// Position x \[m\].
    pub position: [f64; 3],
    /// Velocity u \[m/s\].
    pub velocity: [f64; 3],
    /// Acoustic pressure p \[Pa\].
    pub pressure: f64,
    /// Density ρ \[kg/m³\].
    pub density: f64,
    /// Local speed of sound c \[m/s\].
    pub sound_speed: f64,
    /// Particle mass m \[kg\].
    pub mass: f64,
    /// Smoothing length h \[m\].
    pub smoothing_length: f64,
}

impl AcousticParticle {
    /// Create a new `AcousticParticle`.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        position: [f64; 3],
        velocity: [f64; 3],
        pressure: f64,
        density: f64,
        sound_speed: f64,
        mass: f64,
        smoothing_length: f64,
    ) -> Self {
        Self {
            position,
            velocity,
            pressure,
            density,
            sound_speed,
            mass,
            smoothing_length,
        }
    }

    /// Acoustic kinetic energy ½ ρ |u|² \[J/m³\].
    pub fn kinetic_energy_density(&self) -> f64 {
        0.5 * self.density * dot3(self.velocity, self.velocity)
    }

    /// Acoustic potential energy p² / (2 ρ c²) \[J/m³\].
    pub fn potential_energy_density(&self) -> f64 {
        let denom = 2.0 * self.density * self.sound_speed * self.sound_speed;
        if denom < 1e-300 {
            0.0
        } else {
            self.pressure * self.pressure / denom
        }
    }

    /// Acoustic intensity I = p * u \[W/m²\].
    pub fn acoustic_intensity(&self) -> [f64; 3] {
        scale3(self.velocity, self.pressure)
    }

    /// Mach number M = |u| / c.
    pub fn mach_number(&self) -> f64 {
        if self.sound_speed < 1e-300 {
            0.0
        } else {
            len3(self.velocity) / self.sound_speed
        }
    }
}

// ============================================================================
// AcousticSPH
// ============================================================================

/// Acoustic SPH solver for linear wave propagation.
///
/// Uses the linearised Euler equations in SPH form:
/// Dρ/Dt = -ρ ∇·u
/// ρ Du/Dt = -∇p
/// p = c² (ρ - ρ₀) + p₀
pub struct AcousticSPH {
    /// SPH particles.
    pub particles: Vec<AcousticParticle>,
    /// Reference (ambient) density ρ₀ \[kg/m³\].
    pub reference_density: f64,
    /// Background (ambient) pressure p₀ \[Pa\].
    pub background_pressure: f64,
    /// Current simulation time \[s\].
    pub time: f64,
}

impl AcousticSPH {
    /// Create a new `AcousticSPH` solver.
    pub fn new(reference_density: f64, background_pressure: f64) -> Self {
        Self {
            particles: Vec::new(),
            reference_density,
            background_pressure,
            time: 0.0,
        }
    }

    /// Add a particle to the simulation.
    pub fn add_particle(&mut self, p: AcousticParticle) {
        self.particles.push(p);
    }

    /// Perform one explicit (leapfrog) acoustic SPH time step of size `dt` \[s\].
    ///
    /// Computes pressure-gradient accelerations and continuity updates using
    /// the standard SPH momentum and continuity summations.
    pub fn step(&mut self, dt: f64) {
        let n = self.particles.len();
        let mut accel = vec![[0.0_f64; 3]; n];
        let mut drho_dt = vec![0.0_f64; n];

        // Accumulate SPH interactions
        for i in 0..n {
            for j in 0..n {
                if i == j {
                    continue;
                }
                let r_ij = sub3(self.particles[i].position, self.particles[j].position);
                let r = len3(r_ij);
                let h_avg =
                    0.5 * (self.particles[i].smoothing_length + self.particles[j].smoothing_length);
                if r > 2.0 * h_avg {
                    continue;
                }

                let grad_w = cubic_kernel_grad(r_ij, h_avg);
                let mj = self.particles[j].mass;
                let rho_i = self.particles[i].density;
                let rho_j = self.particles[j].density;

                // Pressure acceleration: -m_j (p_i/ρ_i² + p_j/ρ_j²) ∇W_ij
                let pi = self.particles[i].pressure;
                let pj = self.particles[j].pressure;
                let fac = if rho_i > 1e-300 && rho_j > 1e-300 {
                    mj * (pi / (rho_i * rho_i) + pj / (rho_j * rho_j))
                } else {
                    0.0
                };
                for k in 0..3 {
                    accel[i][k] -= fac * grad_w[k];
                }

                // Continuity: Dρ/Dt = Σ_j m_j (u_i - u_j) · ∇W_ij
                let du = sub3(self.particles[i].velocity, self.particles[j].velocity);
                drho_dt[i] += mj * dot3(du, grad_w);
            }
        }

        // Integrate
        for i in 0..n {
            for k in 0..3 {
                self.particles[i].velocity[k] += dt * accel[i][k];
                self.particles[i].position[k] += dt * self.particles[i].velocity[k];
            }
            self.particles[i].density += dt * drho_dt[i];
            // Linearised equation of state
            let c = self.particles[i].sound_speed;
            self.particles[i].pressure = self.background_pressure
                + c * c * (self.particles[i].density - self.reference_density);
        }

        self.time += dt;
    }

    /// Total acoustic energy in the domain \[J\].
    pub fn total_energy(&self) -> f64 {
        self.particles
            .iter()
            .map(|p| {
                let vol = if p.density > 1e-300 {
                    p.mass / p.density
                } else {
                    0.0
                };
                (p.kinetic_energy_density() + p.potential_energy_density()) * vol
            })
            .sum()
    }
}

// ============================================================================
// PressureWave
// ============================================================================

/// Descriptor for a harmonic travelling pressure wave.
pub struct PressureWave {
    /// Pressure amplitude A \[Pa\].
    pub amplitude: f64,
    /// Frequency f \[Hz\].
    pub frequency: f64,
    /// Wavelength λ \[m\].
    pub wave_length: f64,
    /// Phase velocity c \[m/s\].
    pub phase_velocity: f64,
}

impl PressureWave {
    /// Create a `PressureWave` from amplitude, frequency, and sound speed.
    pub fn new(amplitude: f64, frequency: f64, sound_speed: f64) -> Self {
        let wl = wavelength(frequency, sound_speed);
        Self {
            amplitude,
            frequency,
            wave_length: wl,
            phase_velocity: sound_speed,
        }
    }

    /// Wavenumber k = 2π / λ \[rad/m\].
    pub fn wavenumber(&self) -> f64 {
        if self.wave_length.abs() < 1e-300 {
            0.0
        } else {
            2.0 * PI / self.wave_length
        }
    }

    /// Angular frequency ω = 2π f \[rad/s\].
    pub fn angular_frequency(&self) -> f64 {
        2.0 * PI * self.frequency
    }

    /// Pressure at position x \[m\] and time t \[s\]: p(x,t) = A cos(kx - ωt).
    pub fn pressure_at(&self, x: f64, t: f64) -> f64 {
        let k = self.wavenumber();
        let omega = self.angular_frequency();
        self.amplitude * (k * x - omega * t).cos()
    }

    /// Acoustic intensity I = A² / (2 ρ c) \[W/m²\].
    pub fn intensity(&self, density: f64) -> f64 {
        if density < 1e-300 || self.phase_velocity < 1e-300 {
            return 0.0;
        }
        self.amplitude * self.amplitude / (2.0 * density * self.phase_velocity)
    }

    /// Sound pressure level SPL = 20 log10(A / p_ref) \[dB re 20 µPa\].
    pub fn sound_pressure_level(&self) -> f64 {
        let p_ref = 20e-6; // 20 µPa
        if self.amplitude < 1e-300 {
            return f64::NEG_INFINITY;
        }
        20.0 * (self.amplitude / p_ref).log10()
    }
}

// ============================================================================
// AbsorbingBoundary (PML)
// ============================================================================

/// Perfectly Matched Layer (PML) absorbing boundary condition.
///
/// Applies exponential damping to particles inside the PML region to
/// prevent spurious reflections from domain boundaries.
pub struct AbsorbingBoundary {
    /// PML absorption coefficient σ_x \[1/m\].
    pub sigma_x: f64,
    /// PML absorption coefficient σ_y \[1/m\].
    pub sigma_y: f64,
    /// PML absorption coefficient σ_z \[1/m\].
    pub sigma_z: f64,
    /// PML thickness δ \[m\].
    pub thickness: f64,
    /// Domain boundary positions \[x_min, x_max, y_min, y_max, z_min, z_max\].
    pub bounds: [f64; 6],
}

impl AbsorbingBoundary {
    /// Create a new `AbsorbingBoundary`.
    ///
    /// # Arguments
    /// * `sigma`     – maximum PML absorption coefficient \[1/m\]
    /// * `thickness` – PML layer thickness \[m\]
    /// * `bounds`    – domain extents \[x_min, x_max, y_min, y_max, z_min, z_max\]
    pub fn new(sigma: f64, thickness: f64, bounds: [f64; 6]) -> Self {
        Self {
            sigma_x: sigma,
            sigma_y: sigma,
            sigma_z: sigma,
            thickness,
            bounds,
        }
    }

    /// Compute the local PML damping factor for a given coordinate and axis index (0=x,1=y,2=z).
    fn pml_factor(&self, coord: f64, axis: usize) -> f64 {
        let (lo, hi, sigma) = match axis {
            0 => (self.bounds[0], self.bounds[1], self.sigma_x),
            1 => (self.bounds[2], self.bounds[3], self.sigma_y),
            _ => (self.bounds[4], self.bounds[5], self.sigma_z),
        };
        let d_lo = coord - lo;
        let d_hi = hi - coord;
        let dist = d_lo.min(d_hi).max(0.0);
        if dist < self.thickness {
            let depth = self.thickness - dist;
            let normalized = depth / self.thickness;
            sigma * normalized * normalized
        } else {
            0.0
        }
    }

    /// Apply PML damping to an `AcousticParticle` over time step `dt`.
    ///
    /// Multiplies velocity and pressure by exp(-σ dt) factors in each direction.
    pub fn apply_pml(&self, particle: &mut AcousticParticle, dt: f64) {
        let ax = self.pml_factor(particle.position[0], 0);
        let ay = self.pml_factor(particle.position[1], 1);
        let az = self.pml_factor(particle.position[2], 2);
        let damp_x = (-ax * dt).exp();
        let damp_y = (-ay * dt).exp();
        let damp_z = (-az * dt).exp();
        particle.velocity[0] *= damp_x;
        particle.velocity[1] *= damp_y;
        particle.velocity[2] *= damp_z;
        let mean_damp = (damp_x * damp_y * damp_z).cbrt();
        particle.pressure *= mean_damp;
    }

    /// Check whether a position is inside the PML region.
    pub fn is_in_pml(&self, pos: [f64; 3]) -> bool {
        for axis in 0..3 {
            if self.pml_factor(pos[axis], axis) > 0.0 {
                return true;
            }
        }
        false
    }
}

// ============================================================================
// AcousticSource
// ============================================================================

/// Type of acoustic source.
#[derive(Clone, Copy, PartialEq, Debug)]
pub enum SourceType {
    /// Monopole (isotropic) point source.
    Monopole,
    /// Dipole source (directional).
    Dipole,
    /// Quadrupole source.
    Quadrupole,
}

/// An acoustic point source that emits pressure waves.
pub struct AcousticSource {
    /// Position of the source \[m\].
    pub position: [f64; 3],
    /// Emission frequency \[Hz\].
    pub frequency: f64,
    /// Emission amplitude \[Pa\].
    pub amplitude: f64,
    /// Type of the acoustic source.
    pub source_type: SourceType,
}

impl AcousticSource {
    /// Create a new `AcousticSource`.
    pub fn new(
        position: [f64; 3],
        frequency: f64,
        amplitude: f64,
        source_type: SourceType,
    ) -> Self {
        Self {
            position,
            frequency,
            amplitude,
            source_type,
        }
    }

    /// Instantaneous emission pressure at time `t` \[s\].
    ///
    /// - Monopole: A sin(2π f t)
    /// - Dipole:   A sin(2π f t) * cos(2π f t)   (harmonic cross-term)
    /// - Quadrupole: A sin²(2π f t) * 2 - A        (second harmonic)
    pub fn emit(&self, t: f64) -> f64 {
        let phase = 2.0 * PI * self.frequency * t;
        match self.source_type {
            SourceType::Monopole => self.amplitude * phase.sin(),
            SourceType::Dipole => self.amplitude * phase.sin() * phase.cos(),
            SourceType::Quadrupole => self.amplitude * (2.0 * phase.sin() * phase.sin() - 1.0),
        }
    }

    /// Pressure at observer position `obs` \[m\] at time `t`, sound speed `c`.
    ///
    /// Uses free-field 1/r spherical spreading for monopole:
    /// p(r, t) = A / r * sin(ω(t - r/c))
    pub fn pressure_at(&self, obs: [f64; 3], t: f64, sound_speed: f64) -> f64 {
        let r_vec = sub3(obs, self.position);
        let r = len3(r_vec);
        if r < 1e-10 {
            return self.emit(t);
        }
        let omega = 2.0 * PI * self.frequency;
        let delay = r / sound_speed;
        match self.source_type {
            SourceType::Monopole => self.amplitude / r * (omega * (t - delay)).sin(),
            SourceType::Dipole => self.amplitude / (r * r) * (omega * (t - delay)).sin(),
            SourceType::Quadrupole => self.amplitude / (r * r * r) * (omega * (t - delay)).sin(),
        }
    }

    /// Directivity factor D for dipole (cosine of angle to axis).
    pub fn directivity(&self, obs: [f64; 3], axis: [f64; 3]) -> f64 {
        let r = sub3(obs, self.position);
        let r_len = len3(r);
        let a_len = len3(axis);
        if r_len < 1e-12 || a_len < 1e-12 {
            return 0.0;
        }
        dot3(r, axis) / (r_len * a_len)
    }
}

// ============================================================================
// Free functions
// ============================================================================

/// Speed of sound in sea water using the simplified Mackenzie equation \[m/s\].
///
/// Valid for: 2 ≤ T ≤ 30°C, 25 ≤ S ≤ 40 ppt, 0 ≤ z ≤ 8000 m depth.
///
/// # Arguments
/// * `temperature` – temperature T \[°C\]
/// * `salinity`    – salinity S \[ppt\]
pub fn sound_speed_water(temperature: f64, salinity: f64) -> f64 {
    // Mackenzie (1981) simplified formula at surface (depth = 0)
    1448.96 + 4.591 * temperature - 5.304e-2 * temperature * temperature
        + 2.374e-4 * temperature * temperature * temperature
        + 1.340 * (salinity - 35.0)
}

/// Characteristic acoustic impedance Z = ρ c \[Pa·s/m = rayl\].
///
/// # Arguments
/// * `density`     – medium density ρ \[kg/m³\]
/// * `sound_speed` – speed of sound c \[m/s\]
pub fn acoustic_impedance(density: f64, sound_speed: f64) -> f64 {
    density * sound_speed
}

/// Acoustic wavelength λ = c / f \[m\].
///
/// # Arguments
/// * `frequency`   – frequency f \[Hz\]
/// * `sound_speed` – speed of sound c \[m/s\]
pub fn wavelength(frequency: f64, sound_speed: f64) -> f64 {
    if frequency.abs() < 1e-300 {
        return f64::INFINITY;
    }
    sound_speed / frequency
}

/// CFL time-step limit for acoustic SPH: Δt_CFL = C_CFL * h / c.
///
/// # Arguments
/// * `dx`     – particle spacing / smoothing length h \[m\]
/// * `c`      – sound speed \[m/s\]
/// * `c_cfl`  – CFL safety factor (typically 0.2–0.4)
pub fn cfl_acoustic(dx: f64, c: f64) -> f64 {
    if c.abs() < 1e-300 {
        return f64::INFINITY;
    }
    0.3 * dx / c
}

/// Transmission coefficient T for a normally-incident wave at an impedance interface.
///
/// T = 2 Z2 / (Z1 + Z2)
///
/// # Arguments
/// * `z1` – impedance of medium 1 \[rayl\]
/// * `z2` – impedance of medium 2 \[rayl\]
pub fn transmission_coefficient(z1: f64, z2: f64) -> f64 {
    let denom = z1 + z2;
    if denom.abs() < 1e-300 {
        return 0.0;
    }
    2.0 * z2 / denom
}

/// Reflection coefficient R for a normally-incident wave at an impedance interface.
///
/// R = (Z2 - Z1) / (Z2 + Z1)
///
/// # Arguments
/// * `z1` – impedance of medium 1 \[rayl\]
/// * `z2` – impedance of medium 2 \[rayl\]
pub fn reflection_coefficient(z1: f64, z2: f64) -> f64 {
    let denom = z1 + z2;
    if denom.abs() < 1e-300 {
        return 0.0;
    }
    (z2 - z1) / denom
}

/// Absorption coefficient in dB/m for a plane wave with attenuation α_np \[Np/m\].
pub fn attenuation_db_per_m(alpha_np: f64) -> f64 {
    alpha_np * 20.0 / std::f64::consts::LN_10
}

/// Near-field distance of an acoustic source: d_near = A² / λ (Rayleigh distance).
///
/// # Arguments
/// * `aperture`  – source aperture A \[m\]
/// * `lambda`    – wavelength λ \[m\]
pub fn rayleigh_distance(aperture: f64, lambda: f64) -> f64 {
    if lambda.abs() < 1e-300 {
        return 0.0;
    }
    aperture * aperture / lambda
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    // --- sound_speed_water ---

    #[test]
    fn test_sound_speed_water_typical() {
        // ~1500 m/s at 15°C, 35 ppt
        let c = sound_speed_water(15.0, 35.0);
        assert!(c > 1490.0 && c < 1510.0, "c = {c}");
    }

    #[test]
    fn test_sound_speed_water_cold() {
        // Cold water ~1450 m/s
        let c = sound_speed_water(4.0, 35.0);
        assert!(c > 1440.0 && c < 1480.0, "c = {c}");
    }

    #[test]
    fn test_sound_speed_water_increases_with_temp() {
        let c1 = sound_speed_water(10.0, 35.0);
        let c2 = sound_speed_water(25.0, 35.0);
        assert!(c2 > c1);
    }

    #[test]
    fn test_sound_speed_water_increases_with_salinity() {
        let c1 = sound_speed_water(15.0, 30.0);
        let c2 = sound_speed_water(15.0, 38.0);
        assert!(c2 > c1);
    }

    // --- acoustic_impedance ---

    #[test]
    fn test_acoustic_impedance_water() {
        // Water: rho=1025, c=1500 → Z ≈ 1.54e6 rayl
        let z = acoustic_impedance(1025.0, 1500.0);
        assert!((z - 1_537_500.0).abs() < 1.0);
    }

    #[test]
    fn test_acoustic_impedance_air() {
        // Air: rho=1.2, c=343 → Z ≈ 411.6 rayl
        let z = acoustic_impedance(1.2, 343.0);
        assert!((z - 411.6).abs() < 0.1);
    }

    #[test]
    fn test_acoustic_impedance_zero_density() {
        let z = acoustic_impedance(0.0, 1500.0);
        assert!((z).abs() < 1e-12);
    }

    // --- wavelength ---

    #[test]
    fn test_wavelength_1khz_water() {
        // 1 kHz in water: λ = 1500/1000 = 1.5 m
        let lam = wavelength(1000.0, 1500.0);
        assert!((lam - 1.5).abs() < 1e-10);
    }

    #[test]
    fn test_wavelength_1khz_air() {
        // 1 kHz in air: λ = 343/1000 = 0.343 m
        let lam = wavelength(1000.0, 343.0);
        assert!((lam - 0.343).abs() < 1e-10);
    }

    #[test]
    fn test_wavelength_zero_frequency() {
        let lam = wavelength(0.0, 1500.0);
        assert!(lam.is_infinite());
    }

    #[test]
    fn test_wavelength_increases_with_speed() {
        let l1 = wavelength(1000.0, 1000.0);
        let l2 = wavelength(1000.0, 2000.0);
        assert!(l2 > l1);
    }

    // --- cfl_acoustic ---

    #[test]
    fn test_cfl_acoustic_basic() {
        let dt = cfl_acoustic(0.01, 1500.0);
        assert!((dt - 0.3 * 0.01 / 1500.0).abs() < 1e-20);
    }

    #[test]
    fn test_cfl_acoustic_positive() {
        let dt = cfl_acoustic(0.005, 343.0);
        assert!(dt > 0.0);
    }

    #[test]
    fn test_cfl_acoustic_zero_speed() {
        let dt = cfl_acoustic(0.01, 0.0);
        assert!(dt.is_infinite());
    }

    #[test]
    fn test_cfl_acoustic_smaller_dx_smaller_dt() {
        let dt1 = cfl_acoustic(0.01, 1500.0);
        let dt2 = cfl_acoustic(0.005, 1500.0);
        assert!(dt2 < dt1);
    }

    // --- AcousticParticle ---

    #[test]
    fn test_particle_new() {
        let p = AcousticParticle::new(
            [0.0; 3],
            [1.0, 0.0, 0.0],
            100.0,
            1000.0,
            1500.0,
            0.001,
            0.01,
        );
        assert!((p.sound_speed - 1500.0).abs() < 1e-10);
    }

    #[test]
    fn test_particle_kinetic_energy() {
        let p = AcousticParticle::new([0.0; 3], [1.0, 0.0, 0.0], 0.0, 1000.0, 1500.0, 0.001, 0.01);
        let ke = p.kinetic_energy_density();
        assert!((ke - 500.0).abs() < 1e-8);
    }

    #[test]
    fn test_particle_potential_energy() {
        let p = AcousticParticle::new([0.0; 3], [0.0; 3], 1000.0, 1000.0, 1500.0, 0.001, 0.01);
        let pe = p.potential_energy_density();
        assert!(pe > 0.0);
    }

    #[test]
    fn test_particle_acoustic_intensity() {
        let p = AcousticParticle::new(
            [0.0; 3],
            [2.0, 0.0, 0.0],
            500.0,
            1000.0,
            1500.0,
            0.001,
            0.01,
        );
        let i = p.acoustic_intensity();
        assert!((i[0] - 1000.0).abs() < 1e-8);
    }

    #[test]
    fn test_particle_mach_number() {
        let p = AcousticParticle::new(
            [0.0; 3],
            [150.0, 0.0, 0.0],
            0.0,
            1000.0,
            1500.0,
            0.001,
            0.01,
        );
        assert!((p.mach_number() - 0.1).abs() < 1e-10);
    }

    // --- PressureWave ---

    #[test]
    fn test_pressure_wave_wavenumber() {
        let w = PressureWave::new(100.0, 1000.0, 1500.0);
        // k = 2π f / c = 2π * 1000 / 1500
        let expected = 2.0 * PI * 1000.0 / 1500.0;
        assert!((w.wavenumber() - expected).abs() < 1e-8);
    }

    #[test]
    fn test_pressure_wave_angular_freq() {
        let w = PressureWave::new(100.0, 1000.0, 1500.0);
        assert!((w.angular_frequency() - 2000.0 * PI).abs() < 1e-6);
    }

    #[test]
    fn test_pressure_wave_pressure_at_origin() {
        let w = PressureWave::new(100.0, 0.0, 1500.0);
        // f=0 → cos(0) = 1 for any x, t
        let p = w.pressure_at(0.0, 0.0);
        assert!((p - 100.0).abs() < 1e-8);
    }

    #[test]
    fn test_pressure_wave_intensity() {
        let w = PressureWave::new(100.0, 1000.0, 1500.0);
        let i = w.intensity(1000.0);
        let expected = 100.0 * 100.0 / (2.0 * 1000.0 * 1500.0);
        assert!((i - expected).abs() < 1e-12);
    }

    #[test]
    fn test_pressure_wave_spl() {
        // A = 1 Pa → SPL = 20 log10(1 / 20e-6) ≈ 94 dB
        let w = PressureWave::new(1.0, 1000.0, 1500.0);
        let spl = w.sound_pressure_level();
        assert!((spl - 94.0).abs() < 0.1, "SPL = {spl}");
    }

    // --- AcousticSource ---

    #[test]
    fn test_source_monopole_emit() {
        let src = AcousticSource::new([0.0; 3], 1.0, 1.0, SourceType::Monopole);
        // t=0: sin(0)=0
        assert!((src.emit(0.0)).abs() < 1e-12);
        // t=0.25: sin(π/2)=1
        assert!((src.emit(0.25) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_source_dipole_emit_zero_at_t0() {
        let src = AcousticSource::new([0.0; 3], 1.0, 1.0, SourceType::Dipole);
        // sin(0)*cos(0) = 0
        assert!((src.emit(0.0)).abs() < 1e-12);
    }

    #[test]
    fn test_source_quadrupole_emit_t0() {
        let src = AcousticSource::new([0.0; 3], 1.0, 2.0, SourceType::Quadrupole);
        // 2 * sin²(0) - 1 = -1, times amplitude 2 = -2
        assert!((src.emit(0.0) - (-2.0)).abs() < 1e-10);
    }

    #[test]
    fn test_source_pressure_at_far_field() {
        let src = AcousticSource::new([0.0; 3], 1000.0, 100.0, SourceType::Monopole);
        let p = src.pressure_at([100.0, 0.0, 0.0], 0.1, 1500.0);
        // Should be of order 100/100 = 1 Pa
        assert!(p.abs() <= 2.0); // within reasonable range
    }

    #[test]
    fn test_source_directivity() {
        let src = AcousticSource::new([0.0; 3], 1000.0, 1.0, SourceType::Dipole);
        let d = src.directivity([1.0, 0.0, 0.0], [1.0, 0.0, 0.0]);
        assert!((d - 1.0).abs() < 1e-10);
    }

    // --- AbsorbingBoundary ---

    #[test]
    fn test_pml_inside_domain_no_damping() {
        let pml = AbsorbingBoundary::new(100.0, 0.1, [-1.0, 1.0, -1.0, 1.0, -1.0, 1.0]);
        // Center particle should not be in PML
        assert!(!pml.is_in_pml([0.0, 0.0, 0.0]));
    }

    #[test]
    fn test_pml_near_boundary_is_in_pml() {
        let pml = AbsorbingBoundary::new(100.0, 0.2, [-1.0, 1.0, -1.0, 1.0, -1.0, 1.0]);
        assert!(pml.is_in_pml([0.95, 0.0, 0.0]));
    }

    #[test]
    fn test_pml_apply_reduces_velocity() {
        let pml = AbsorbingBoundary::new(1000.0, 0.5, [-1.0, 1.0, -1.0, 1.0, -1.0, 1.0]);
        let mut p = AcousticParticle::new(
            [0.9, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            100.0,
            1000.0,
            1500.0,
            0.001,
            0.01,
        );
        pml.apply_pml(&mut p, 0.001);
        assert!(p.velocity[0] < 1.0);
    }

    #[test]
    fn test_pml_no_effect_interior() {
        let pml = AbsorbingBoundary::new(1000.0, 0.1, [-1.0, 1.0, -1.0, 1.0, -1.0, 1.0]);
        let mut p = AcousticParticle::new(
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            100.0,
            1000.0,
            1500.0,
            0.001,
            0.01,
        );
        pml.apply_pml(&mut p, 0.001);
        // Interior: no damping
        assert!((p.velocity[0] - 1.0).abs() < 1e-12);
    }

    // --- AcousticSPH ---

    #[test]
    fn test_acoustic_sph_new() {
        let sph = AcousticSPH::new(1000.0, 101325.0);
        assert_eq!(sph.particles.len(), 0);
        assert!((sph.time).abs() < 1e-12);
    }

    #[test]
    fn test_acoustic_sph_add_particle() {
        let mut sph = AcousticSPH::new(1000.0, 101325.0);
        let p = AcousticParticle::new([0.0; 3], [0.0; 3], 0.0, 1000.0, 1500.0, 1.0, 0.1);
        sph.add_particle(p);
        assert_eq!(sph.particles.len(), 1);
    }

    #[test]
    fn test_acoustic_sph_step_advances_time() {
        let mut sph = AcousticSPH::new(1000.0, 101325.0);
        sph.step(0.001);
        assert!((sph.time - 0.001).abs() < 1e-12);
    }

    #[test]
    fn test_acoustic_sph_single_particle_no_interaction() {
        let mut sph = AcousticSPH::new(1000.0, 0.0);
        let p = AcousticParticle::new([0.0; 3], [1.0, 0.0, 0.0], 0.0, 1000.0, 1500.0, 1.0, 0.01);
        sph.add_particle(p);
        sph.step(0.001);
        // Without neighbors, velocity unchanged; position advected
        assert!((sph.particles[0].position[0] - 0.001).abs() < 1e-12);
    }

    // --- Transmission / reflection coefficients ---

    #[test]
    fn test_transmission_coefficient_same_media() {
        // Z1 = Z2 → T = 1
        let t = transmission_coefficient(1000.0, 1000.0);
        assert!((t - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_reflection_coefficient_same_media() {
        // Z1 = Z2 → R = 0
        let r = reflection_coefficient(1000.0, 1000.0);
        assert!(r.abs() < 1e-10);
    }

    #[test]
    fn test_energy_conservation_at_interface() {
        // |R|² + Z1/Z2 * |T|² = 1
        let z1 = 400.0_f64; // air
        let z2 = 1_500_000.0_f64; // water
        let r = reflection_coefficient(z1, z2);
        let t = transmission_coefficient(z1, z2);
        let energy = r * r + (z1 / z2) * t * t;
        assert!((energy - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_rayleigh_distance() {
        // A=0.1 m, λ=0.15 m → d = 0.01/0.15 ≈ 0.0667 m
        let d = rayleigh_distance(0.1, 0.15);
        assert!((d - 0.01 / 0.15).abs() < 1e-12);
    }

    #[test]
    fn test_attenuation_db_per_m() {
        // α_np=1 Np/m → 20/ln(10) ≈ 8.686 dB/m
        let db = attenuation_db_per_m(1.0);
        assert!((db - 8.686).abs() < 0.001);
    }

    #[test]
    fn test_source_type_equality() {
        assert_eq!(SourceType::Monopole, SourceType::Monopole);
        assert_ne!(SourceType::Monopole, SourceType::Dipole);
    }
}
