// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Acoustic streaming Lattice Boltzmann Method module.
//!
//! Implements acoustic streaming phenomena on a D2Q9 lattice:
//!
//! - **Acoustic pressure field**: D2Q9 distribution functions for sound
//! - **Eckart streaming**: bulk acoustic streaming from attenuation-driven body force
//! - **Rayleigh streaming**: near-wall boundary-layer streaming
//! - **Acoustic radiation force**: Gorkov potential, primary radiation force
//! - **Standing wave patterns**: 1-D and 2-D pressure nodal structures
//! - **Traveling wave streaming**: unidirectional acoustic streaming
//! - **Piezoelectric transducer model**: pressure source at boundary
//! - **Acoustic microfluidics**: particle focusing in half-wavelength resonators
//! - **Particle manipulation**: acoustic tweezers force on spherical particles
//! - **Acoustic levitation**: gravity-balanced radiation force
//! - **Reynolds stress tensor**: acoustic forcing via Lighthill analogy
//! - **Nyborg body force**: acoustic streaming body force formulation

#![allow(dead_code)]
#![allow(clippy::too_many_arguments)]

use std::f64::consts::PI;

// ============================================================================
// Constants
// ============================================================================

/// D2Q9 lattice speed of sound squared: cs² = 1/3.
pub const CS2: f64 = 1.0 / 3.0;

/// D2Q9 lattice speed of sound: cs = 1/√3.
pub const CS: f64 = 0.577_350_269_189_625_8_f64;

/// Reference sound speed in water at 25 °C \[m/s\].
pub const C0_WATER: f64 = 1_496.0;

/// Reference density of water at 25 °C \[kg/m³\].
pub const RHO0_WATER: f64 = 997.0;

/// Reference sound speed in air at 20 °C \[m/s\].
pub const C0_AIR: f64 = 343.0;

/// Reference air density at 20 °C \[kg/m³\].
pub const RHO0_AIR: f64 = 1.2041;

/// Dynamic viscosity of water at 25 °C \[Pa·s\].
pub const MU_WATER: f64 = 8.9e-4;

/// Acoustic reference pressure (20 µPa in air).
pub const P_REF: f64 = 2.0e-5;

// ============================================================================
// 1. D2Q9 Velocity Set
// ============================================================================

/// D2Q9 discrete velocity vectors \[cx, cy\] for all 9 directions.
///
/// Direction ordering: (0,0), (1,0), (0,1), (-1,0), (0,-1), (1,1), (-1,1), (-1,-1), (1,-1).
pub const D2Q9_CX: [f64; 9] = [0.0, 1.0, 0.0, -1.0, 0.0, 1.0, -1.0, -1.0, 1.0];

/// D2Q9 cy components.
pub const D2Q9_CY: [f64; 9] = [0.0, 0.0, 1.0, 0.0, -1.0, 1.0, 1.0, -1.0, -1.0];

/// D2Q9 equilibrium weights.
pub const D2Q9_W: [f64; 9] = [
    4.0 / 9.0,
    1.0 / 9.0,
    1.0 / 9.0,
    1.0 / 9.0,
    1.0 / 9.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
];

/// D2Q9 opposite direction indices (for bounce-back).
pub const D2Q9_OPP: [usize; 9] = [0, 3, 4, 1, 2, 7, 8, 5, 6];

// ============================================================================
// 2. Equilibrium Distribution
// ============================================================================

/// Compute the D2Q9 equilibrium distribution function f_i^{eq}.
///
/// f_i^{eq} = w_i ρ \[1 + (c_i·u)/cs² + (c_i·u)²/(2cs⁴) - u²/(2cs²)\]
///
/// # Arguments
/// * `rho` – local density
/// * `ux`  – x-velocity
/// * `uy`  – y-velocity
pub fn d2q9_equilibrium(rho: f64, ux: f64, uy: f64) -> [f64; 9] {
    let u2 = ux * ux + uy * uy;
    let mut feq = [0.0_f64; 9];
    for i in 0..9 {
        let cu = D2Q9_CX[i] * ux + D2Q9_CY[i] * uy;
        feq[i] =
            D2Q9_W[i] * rho * (1.0 + cu / CS2 + cu * cu / (2.0 * CS2 * CS2) - u2 / (2.0 * CS2));
    }
    feq
}

/// Compute the D2Q9 equilibrium for acoustic (linearized) equations.
///
/// Linearized form for small perturbations around (ρ0, 0, 0):
/// f_i^{eq} = w_i \[ρ' + ρ0 (c_i·u)/cs²\]
///
/// # Arguments
/// * `rho_prime` – density perturbation
/// * `rho0`      – background density
/// * `ux`        – velocity perturbation x
/// * `uy`        – velocity perturbation y
pub fn d2q9_equilibrium_acoustic(rho_prime: f64, rho0: f64, ux: f64, uy: f64) -> [f64; 9] {
    let mut feq = [0.0_f64; 9];
    for i in 0..9 {
        let cu = D2Q9_CX[i] * ux + D2Q9_CY[i] * uy;
        feq[i] = D2Q9_W[i] * (rho_prime + rho0 * cu / CS2);
    }
    feq
}

/// Compute macroscopic density and momentum from distribution functions.
///
/// ρ = Σ_i f_i,  ρ u_x = Σ_i f_i c_{ix},  ρ u_y = Σ_i f_i c_{iy}
///
/// Returns `(rho, ux, uy)`.
pub fn d2q9_macroscopic(f: &[f64; 9]) -> (f64, f64, f64) {
    let rho: f64 = f.iter().sum();
    let rho_ux: f64 = (0..9).map(|i| f[i] * D2Q9_CX[i]).sum();
    let rho_uy: f64 = (0..9).map(|i| f[i] * D2Q9_CY[i]).sum();
    if rho.abs() < 1e-300 {
        (0.0, 0.0, 0.0)
    } else {
        (rho, rho_ux / rho, rho_uy / rho)
    }
}

// ============================================================================
// 3. BGK Collision with Acoustic Body Force
// ============================================================================

/// BGK collision operator with optional body force term (Guo forcing).
///
/// f_i^* = f_i - (f_i - f_i^{eq}) / τ + Δt F_i
///
/// where F_i is the Guo body force term:
/// F_i = w_i (1 - 1/(2τ)) \[(c_i - u)/cs² + (c_i·u)c_i/cs⁴\] · F
///
/// # Arguments
/// * `f`   – current distribution functions
/// * `tau` – relaxation time
/// * `rho` – local density
/// * `ux`  – local x-velocity
/// * `uy`  – local y-velocity
/// * `fx`  – body force x-component
/// * `fy`  – body force y-component
/// * `dt`  – time step (lattice units = 1.0 typically)
pub fn bgk_collision_with_force(
    f: &[f64; 9],
    tau: f64,
    rho: f64,
    ux: f64,
    uy: f64,
    fx: f64,
    fy: f64,
    dt: f64,
) -> [f64; 9] {
    let feq = d2q9_equilibrium(rho, ux, uy);
    let factor = D2Q9_W.map(|w| w * (1.0 - 0.5 / tau));
    let mut f_out = [0.0_f64; 9];
    for i in 0..9 {
        let cx = D2Q9_CX[i];
        let cy = D2Q9_CY[i];
        let cu = cx * ux + cy * uy;
        let force_i = factor[i]
            * ((cx - ux) / CS2 * fx
                + (cy - uy) / CS2 * fy
                + cu / (CS2 * CS2) * (cx * fx + cy * fy));
        f_out[i] = f[i] - (f[i] - feq[i]) / tau + dt * force_i;
    }
    f_out
}

/// Simple BGK collision without body force.
///
/// f_i^* = f_i - (f_i - f_i^{eq}) / τ
pub fn bgk_collision(f: &[f64; 9], tau: f64, rho: f64, ux: f64, uy: f64) -> [f64; 9] {
    bgk_collision_with_force(f, tau, rho, ux, uy, 0.0, 0.0, 1.0)
}

// ============================================================================
// 4. Streaming Step
// ============================================================================

/// Perform the D2Q9 streaming (propagation) step on a 2-D grid.
///
/// Distributes post-collision populations to neighbouring nodes.
/// Periodic boundary conditions are applied.
///
/// # Arguments
/// * `f`  – flat array of distribution functions, length 9 * nx * ny
/// * `nx` – number of x grid points
/// * `ny` – number of y grid points
///
/// Returns new distribution function array.
pub fn stream_d2q9(f: &[f64], nx: usize, ny: usize) -> Vec<f64> {
    assert_eq!(f.len(), 9 * nx * ny);
    let mut f_out = vec![0.0_f64; 9 * nx * ny];
    for iy in 0..ny {
        for ix in 0..nx {
            for i in 0..9 {
                let cx = D2Q9_CX[i] as isize;
                let cy = D2Q9_CY[i] as isize;
                let src_x = (ix as isize - cx).rem_euclid(nx as isize) as usize;
                let src_y = (iy as isize - cy).rem_euclid(ny as isize) as usize;
                let dst_idx = i + 9 * (iy * nx + ix);
                let src_idx = i + 9 * (src_y * nx + src_x);
                f_out[dst_idx] = f[src_idx];
            }
        }
    }
    f_out
}

// ============================================================================
// 5. Acoustic Pressure Field
// ============================================================================

/// Compute the acoustic pressure perturbation at every grid node.
///
/// p'(x) = cs² ρ'(x) = cs² (ρ(x) - ρ_0)
///
/// # Arguments
/// * `f`    – distribution functions, length 9 * nx * ny
/// * `nx`   – grid x size
/// * `ny`   – grid y size
/// * `rho0` – background (mean) density
///
/// Returns flat array of pressure perturbations.
pub fn acoustic_pressure_field(f: &[f64], nx: usize, ny: usize, rho0: f64) -> Vec<f64> {
    let n = nx * ny;
    assert_eq!(f.len(), 9 * n);
    let mut p = vec![0.0_f64; n];
    for node in 0..n {
        let fi: [f64; 9] = std::array::from_fn(|i| f[i + 9 * node]);
        let (rho, _, _) = d2q9_macroscopic(&fi);
        p[node] = CS2 * (rho - rho0);
    }
    p
}

/// Sound pressure level (SPL) in decibels from a pressure amplitude.
///
/// SPL = 20 log₁₀(|p'| / p_ref)
pub fn sound_pressure_level(p_amplitude: f64) -> f64 {
    20.0 * (p_amplitude.abs() / P_REF).log10()
}

/// Acoustic energy density for a plane wave.
///
/// E = p²/(2 ρ₀ c₀²)
///
/// # Arguments
/// * `p_rms`   – RMS acoustic pressure
/// * `rho0`    – fluid density
/// * `c0`      – speed of sound
pub fn acoustic_energy_density(p_rms: f64, rho0: f64, c0: f64) -> f64 {
    p_rms * p_rms / (rho0 * c0 * c0)
}

/// Compute the RMS pressure over a spatial field.
pub fn rms_pressure(p_field: &[f64]) -> f64 {
    let n = p_field.len();
    if n == 0 {
        return 0.0;
    }
    (p_field.iter().map(|&p| p * p).sum::<f64>() / n as f64).sqrt()
}

// ============================================================================
// 6. Standing Wave Patterns
// ============================================================================

/// Pressure field of a 1-D standing wave.
///
/// p(x, t) = P₀ cos(kx) cos(ωt)
///
/// # Arguments
/// * `x`   – spatial position \[m\]
/// * `t`   – time \[s\]
/// * `p0`  – pressure amplitude \[Pa\]
/// * `k`   – wavenumber \[rad/m\]
/// * `omega` – angular frequency \[rad/s\]
pub fn standing_wave_pressure_1d(x: f64, t: f64, p0: f64, k: f64, omega: f64) -> f64 {
    p0 * (k * x).cos() * (omega * t).cos()
}

/// Velocity field of a 1-D standing wave.
///
/// u(x, t) = (P₀ / ρ₀ c₀) sin(kx) sin(ωt)
pub fn standing_wave_velocity_1d(
    x: f64,
    t: f64,
    p0: f64,
    k: f64,
    omega: f64,
    rho0: f64,
    c0: f64,
) -> f64 {
    (p0 / (rho0 * c0)) * (k * x).sin() * (omega * t).sin()
}

/// Find pressure nodal positions in a standing wave resonator of length L.
///
/// For a closed-closed resonator: nodes at x = (2n-1)λ/4, n=1,2,...
/// Returns node x-coordinates for the n-th mode.
pub fn standing_wave_nodes(length: f64, n_mode: usize) -> Vec<f64> {
    let lambda = 2.0 * length / n_mode as f64;
    let n_nodes = n_mode - 1;
    (1..=n_nodes).map(|k| k as f64 * lambda / 2.0).collect()
}

/// 2-D standing wave pressure field p(x,y) = P₀ cos(kx x) cos(ky y).
pub fn standing_wave_pressure_2d(x: f64, y: f64, p0: f64, kx: f64, ky: f64) -> f64 {
    p0 * (kx * x).cos() * (ky * y).cos()
}

// ============================================================================
// 7. Traveling Wave Acoustic Streaming (Eckart Streaming)
// ============================================================================

/// Eckart streaming body force magnitude.
///
/// The acoustic streaming force per unit volume (Eckart 1948):
/// F = α_abs E / c₀ = (α_abs / ρ₀ c₀²) p_rms²
///
/// where α_abs is the acoustic absorption coefficient \[1/m\].
///
/// # Arguments
/// * `p_rms`  – RMS acoustic pressure \[Pa\]
/// * `alpha_abs` – absorption coefficient \[m⁻¹\]
/// * `rho0`   – fluid density \[kg/m³\]
/// * `c0`     – speed of sound \[m/s\]
pub fn eckart_streaming_force(p_rms: f64, alpha_abs: f64, rho0: f64, c0: f64) -> f64 {
    alpha_abs * p_rms * p_rms / (rho0 * c0 * c0)
}

/// Eckart streaming velocity (simplified, far-field).
///
/// Streaming velocity magnitude in a focused beam:
/// V_stream ≈ F · L / μ
///
/// where L is the acoustic path length and μ the dynamic viscosity.
///
/// # Arguments
/// * `f_eckart` – Eckart body force \[N/m³\]
/// * `length`   – path length \[m\]
/// * `mu`       – dynamic viscosity \[Pa·s\]
pub fn eckart_streaming_velocity(f_eckart: f64, length: f64, mu: f64) -> f64 {
    f_eckart * length / mu
}

/// Traveling wave pressure field.
///
/// p(x, t) = P₀ exp(-α_abs x) cos(k x - ω t)
pub fn traveling_wave_pressure(x: f64, t: f64, p0: f64, k: f64, omega: f64, alpha_abs: f64) -> f64 {
    p0 * (-alpha_abs * x).exp() * (k * x - omega * t).cos()
}

/// Apply Eckart streaming body force to the LBM grid.
///
/// Adds directional forcing F_x = F_eckart along the +x direction.
/// Returns updated distribution functions with the force applied.
///
/// # Arguments
/// * `f`        – current distributions, length 9 * nx * ny
/// * `nx`, `ny` – grid dimensions
/// * `tau`      – relaxation time
/// * `p_rms_field` – RMS pressure at each node (length nx * ny)
/// * `alpha_abs` – absorption coefficient
/// * `rho0`, `c0` – reference density and sound speed
/// * `mu`       – dynamic viscosity (for velocity correction)
pub fn apply_eckart_streaming(
    f: &[f64],
    nx: usize,
    ny: usize,
    tau: f64,
    p_rms_field: &[f64],
    alpha_abs: f64,
    rho0: f64,
    c0: f64,
    _mu: f64,
) -> Vec<f64> {
    let n = nx * ny;
    assert_eq!(f.len(), 9 * n);
    assert_eq!(p_rms_field.len(), n);
    let mut f_out = f.to_vec();
    for node in 0..n {
        let fi: [f64; 9] = std::array::from_fn(|i| f[i + 9 * node]);
        let (rho, ux, uy) = d2q9_macroscopic(&fi);
        let fx = eckart_streaming_force(p_rms_field[node], alpha_abs, rho0, c0);
        let f_new = bgk_collision_with_force(&fi, tau, rho, ux, uy, fx, 0.0, 1.0);
        for i in 0..9 {
            f_out[i + 9 * node] = f_new[i];
        }
    }
    f_out
}

// ============================================================================
// 8. Rayleigh Streaming
// ============================================================================

/// Rayleigh streaming velocity profile between two parallel plates.
///
/// For a standing wave between plates at y=0 and y=d, the Rayleigh streaming
/// second-order velocity (outer solution) is:
///
/// u₂(y) = -(3 U₀²)/(8 c₀) \[1 - (2y/d - 1)²\] × sin(2kx)
///
/// This function returns u₂(y) at fixed x.
///
/// # Arguments
/// * `y`    – transverse coordinate \[0, d\]
/// * `d`    – channel half-gap (full gap = 2d) \[m\]
/// * `x`    – axial position
/// * `u0`   – acoustic velocity amplitude \[m/s\]
/// * `k`    – wavenumber \[rad/m\]
/// * `c0`   – speed of sound \[m/s\]
pub fn rayleigh_streaming_velocity(y: f64, d: f64, x: f64, u0: f64, k: f64, c0: f64) -> f64 {
    let eta = 2.0 * y / d - 1.0; // normalized coordinate in [-1,1]
    -(3.0 * u0 * u0) / (8.0 * c0) * (1.0 - eta * eta) * (2.0 * k * x).sin()
}

/// Rayleigh streaming mass transport velocity (Stokes drift corrected).
///
/// u_mass = u₂ + u_Stokes
/// For a standing wave: u_mass ≈ -(3/4)(U₀²/c₀) sin(2kx) near the wall.
pub fn rayleigh_mass_transport_velocity(x: f64, u0: f64, k: f64, c0: f64) -> f64 {
    -(3.0 / 4.0) * (u0 * u0 / c0) * (2.0 * k * x).sin()
}

// ============================================================================
// 9. Acoustic Radiation Force
// ============================================================================

/// Gorkov potential for a spherical particle in an acoustic field.
///
/// U = (2π R³ / 3) \[f₁ p²/(ρ₀ c₀²) - (3/2) f₂ ρ₀ v²\]
///
/// where f₁ = 1 - ρ₀ c₀² / (ρ_p c_p²), f₂ = 2(ρ_p - ρ₀)/(2ρ_p + ρ₀).
///
/// # Arguments
/// * `radius`    – particle radius \[m\]
/// * `p_amp`     – acoustic pressure amplitude \[Pa\]
/// * `v_amp`     – acoustic velocity amplitude \[m/s\]
/// * `rho0`      – fluid density \[kg/m³\]
/// * `c0`        – fluid speed of sound \[m/s\]
/// * `rho_p`     – particle density \[kg/m³\]
/// * `c_p`       – particle speed of sound \[m/s\]
pub fn gorkov_potential(
    radius: f64,
    p_amp: f64,
    v_amp: f64,
    rho0: f64,
    c0: f64,
    rho_p: f64,
    c_p: f64,
) -> f64 {
    let f1 = 1.0 - rho0 * c0 * c0 / (rho_p * c_p * c_p);
    let f2 = 2.0 * (rho_p - rho0) / (2.0 * rho_p + rho0);
    let v_coeff = 2.0 * PI * radius.powi(3) / 3.0;
    v_coeff * (f1 * p_amp * p_amp / (rho0 * c0 * c0) - 1.5 * f2 * rho0 * v_amp * v_amp)
}

/// Primary acoustic radiation force on a sphere in a standing wave.
///
/// F_rad = -∇U = -4π k R³ E_ac Φ(f₁,f₂) sin(2kx)
///
/// where Φ = f₁/3 - f₂/2 is the acoustic contrast factor and E_ac = p₀²/(4ρ₀c₀²).
///
/// # Arguments
/// * `radius`  – particle radius \[m\]
/// * `k`       – wavenumber \[rad/m\]
/// * `x`       – particle x-position \[m\]
/// * `p0`      – standing wave pressure amplitude \[Pa\]
/// * `rho0`    – fluid density
/// * `c0`      – fluid sound speed
/// * `rho_p`   – particle density
/// * `c_p`     – particle sound speed
pub fn acoustic_radiation_force(
    radius: f64,
    k: f64,
    x: f64,
    p0: f64,
    rho0: f64,
    c0: f64,
    rho_p: f64,
    c_p: f64,
) -> f64 {
    let f1 = 1.0 - rho0 * c0 * c0 / (rho_p * c_p * c_p);
    let f2 = 2.0 * (rho_p - rho0) / (2.0 * rho_p + rho0);
    let phi = f1 / 3.0 - f2 / 2.0;
    let e_ac = p0 * p0 / (4.0 * rho0 * c0 * c0);
    -4.0 * PI * k * radius.powi(3) * e_ac * phi * (2.0 * k * x).sin()
}

/// Acoustic contrast factor Φ(f₁, f₂) for a spherical particle.
///
/// Positive Φ → particle focuses at pressure nodes (antinodes of velocity).
/// Negative Φ → particle focuses at pressure antinodes.
pub fn acoustic_contrast_factor(rho0: f64, c0: f64, rho_p: f64, c_p: f64) -> f64 {
    let f1 = 1.0 - rho0 * c0 * c0 / (rho_p * c_p * c_p);
    let f2 = 2.0 * (rho_p - rho0) / (2.0 * rho_p + rho0);
    f1 / 3.0 - f2 / 2.0
}

// ============================================================================
// 10. Piezoelectric Transducer Model
// ============================================================================

/// Piezoelectric transducer source: pressure waveform at a boundary.
///
/// Models a piezoelectric transducer driven at frequency f with amplitude A.
/// Output pressure: p_src(t) = A × sin(2π f t) × envelope(t)
///
/// # Arguments
/// * `t`         – time \[s\]
/// * `amplitude` – peak pressure amplitude \[Pa\]
/// * `frequency` – drive frequency \[Hz\]
/// * `rise_time` – rise time for Hann window \[s\] (0 → no windowing)
pub fn piezo_pressure_source(t: f64, amplitude: f64, frequency: f64, rise_time: f64) -> f64 {
    let envelope = if rise_time > 0.0 && t < rise_time {
        0.5 * (1.0 - (PI * t / rise_time).cos())
    } else {
        1.0
    };
    amplitude * (2.0 * PI * frequency * t).sin() * envelope
}

/// Impedance of a piezoelectric transducer (simplified one-port model).
///
/// Z(f) = R_s + j(ωL_s - 1/(ωC_s)) + 1/(jωC_0)
///
/// Returns (real, imaginary) parts of impedance.
///
/// # Arguments
/// * `freq`   – frequency \[Hz\]
/// * `r_s`    – series resistance \[Ω\]
/// * `l_s`    – series inductance \[H\]
/// * `c_s`    – series capacitance \[F\]
/// * `c0`     – parallel (clamped) capacitance \[F\]
pub fn piezo_impedance(freq: f64, r_s: f64, l_s: f64, c_s: f64, c0_cap: f64) -> (f64, f64) {
    let omega = 2.0 * PI * freq;
    let z_series_re = r_s;
    let z_series_im = omega * l_s - 1.0 / (omega * c_s);
    let z_parallel_im = -1.0 / (omega * c0_cap);
    // Series branch in parallel with C0
    let z_re = z_series_re;
    let z_im = z_series_im + z_parallel_im;
    (z_re, z_im)
}

/// Resonant frequency of a piezoelectric transducer (series resonance).
///
/// f_r = 1 / (2π √(L_s C_s))
pub fn piezo_resonant_frequency(l_s: f64, c_s: f64) -> f64 {
    1.0 / (2.0 * PI * (l_s * c_s).sqrt())
}

// ============================================================================
// 11. Acoustic Microfluidics — Particle Focusing
// ============================================================================

/// Equilibrium focusing position for a particle in a half-wavelength resonator.
///
/// Particles with positive acoustic contrast factor focus at the pressure node
/// located at the center of a half-wavelength resonator (x = λ/4 from each wall).
///
/// Returns the ideal equilibrium x-position.
///
/// # Arguments
/// * `channel_width` – resonator width \[m\]
pub fn half_wavelength_focus_position(channel_width: f64) -> f64 {
    channel_width / 2.0
}

/// Time to focus a particle to the pressure node.
///
/// Approximated by integrating the radiation force equation of motion:
/// m dx/dt ≈ F_rad(x) (overdamped in viscous drag).
///
/// τ_focus ~ π μ R / (2 k² R³ E_ac Φ) for small displacements.
///
/// # Arguments
/// * `radius`   – particle radius \[m\]
/// * `k`        – wavenumber \[rad/m\]
/// * `e_ac`     – acoustic energy density \[J/m³\]
/// * `phi`      – acoustic contrast factor
/// * `mu`       – fluid dynamic viscosity \[Pa·s\]
pub fn particle_focusing_time(radius: f64, k: f64, e_ac: f64, phi: f64, mu: f64) -> f64 {
    let numerator = 3.0 * mu;
    let denominator = 4.0 * k * k * radius * radius * e_ac * phi.abs();
    if denominator.abs() < 1e-300 {
        return f64::INFINITY;
    }
    numerator / denominator
}

/// Particle trajectory in 1-D acoustic focusing (overdamped dynamics).
///
/// Solves m ẋ = F_rad(x) - 6πμR ẋ + F_drag ... simplified as
/// 6πμR ẋ = F_rad(x) (Stokes drag dominates).
///
/// Returns `(t_array, x_array)`.
///
/// # Arguments
/// * `x0`        – initial position \[m\]
/// * `radius`    – particle radius \[m\]
/// * `k`         – wavenumber \[rad/m\]
/// * `p0`        – pressure amplitude \[Pa\]
/// * `rho0`      – fluid density
/// * `c0`        – speed of sound
/// * `rho_p`     – particle density
/// * `c_p`       – particle sound speed
/// * `mu`        – viscosity
/// * `t_end`     – final time \[s\]
/// * `n_steps`   – number of time steps
pub fn particle_focusing_trajectory(
    x0: f64,
    radius: f64,
    k: f64,
    p0: f64,
    rho0: f64,
    c0: f64,
    rho_p: f64,
    c_p: f64,
    mu: f64,
    t_end: f64,
    n_steps: usize,
) -> (Vec<f64>, Vec<f64>) {
    let dt = t_end / n_steps as f64;
    let drag = 6.0 * PI * mu * radius;
    let mut t = vec![0.0_f64; n_steps + 1];
    let mut x = vec![0.0_f64; n_steps + 1];
    x[0] = x0;
    for step in 0..n_steps {
        t[step + 1] = (step + 1) as f64 * dt;
        let f_rad = acoustic_radiation_force(radius, k, x[step], p0, rho0, c0, rho_p, c_p);
        let velocity = f_rad / drag;
        x[step + 1] = x[step] + dt * velocity;
    }
    (t, x)
}

// ============================================================================
// 12. Acoustic Levitation
// ============================================================================

/// Acoustic radiation force required to levitate a sphere against gravity.
///
/// Balance: F_rad = m g = (4/3) π R³ ρ_p g
///
/// Returns the required vertical radiation force \[N\].
pub fn levitation_force_required(radius: f64, rho_p: f64, g: f64) -> f64 {
    (4.0 / 3.0) * PI * radius.powi(3) * rho_p * g
}

/// Check if levitation is possible for a given acoustic configuration.
///
/// Returns true if the maximum radiation force exceeds gravity.
///
/// # Arguments
/// * `radius`    – particle radius \[m\]
/// * `rho_p`     – particle density \[kg/m³\]
/// * `g`         – gravitational acceleration \[m/s²\]
/// * `p0`        – acoustic pressure amplitude \[Pa\]
/// * `k`         – wavenumber
/// * `rho0`      – fluid density
/// * `c0`        – speed of sound
/// * `c_p`       – particle sound speed
pub fn can_levitate(
    radius: f64,
    rho_p: f64,
    g: f64,
    p0: f64,
    k: f64,
    rho0: f64,
    c0: f64,
    c_p: f64,
) -> bool {
    // Maximum radiation force (at x = λ/8 from node)
    let x_max = PI / (8.0 * k);
    let f_max = acoustic_radiation_force(radius, k, x_max, p0, rho0, c0, rho_p, c_p).abs();
    let f_grav = levitation_force_required(radius, rho_p, g);
    f_max >= f_grav
}

/// Acoustic levitation equilibrium height (simplified 1-D model).
///
/// Finds the lowest anti-node position above the reflector surface (x=0).
/// For a standing wave above a reflector: pressure node at x = λ/4.
///
/// Returns the equilibrium height \[m\].
pub fn levitation_equilibrium_height(c0: f64, frequency: f64) -> f64 {
    let wavelength = c0 / frequency;
    wavelength / 4.0
}

// ============================================================================
// 13. Reynolds Stress Tensor (Acoustic)
// ============================================================================

/// Acoustic Reynolds stress tensor component T_{ij} = ρ₀ <u_i u_j>.
///
/// For a plane wave propagating in +x direction:
/// T_{xx} = (1/2) ρ₀ U₀²,  T_{yy} = 0,  T_{xy} = 0
///
/// Returns the 2×2 tensor as \[\[Txx, Txy\\], \[Tyx, Tyy\]].
///
/// # Arguments
/// * `u0`    – acoustic velocity amplitude \[m/s\]
/// * `rho0`  – mean fluid density \[kg/m³\]
/// * `theta` – propagation angle \[rad\]
pub fn acoustic_reynolds_stress(u0: f64, rho0: f64, theta: f64) -> [[f64; 2]; 2] {
    let ux = u0 * theta.cos();
    let uy = u0 * theta.sin();
    let txx = 0.5 * rho0 * ux * ux;
    let txy = 0.5 * rho0 * ux * uy;
    let tyy = 0.5 * rho0 * uy * uy;
    [[txx, txy], [txy, tyy]]
}

/// Divergence of the Reynolds stress tensor (body force vector).
///
/// f_i = -∂T_{ij}/∂x_j  (acoustic streaming driving force)
///
/// For a Gaussian beam with amplitude P(x):
/// ∂T_{xx}/∂x = (P²/ρ₀c₀²) × α_abs × exp(-2α_abs x)
///
/// Returns (fx, fy) at position (x, y).
///
/// # Arguments
/// * `x`, `y` – position \[m\]
/// * `p0`     – peak pressure \[Pa\]
/// * `alpha_abs` – absorption coefficient \[1/m\]
/// * `rho0`, `c0` – fluid reference values
pub fn reynolds_stress_body_force(
    x: f64,
    _y: f64,
    p0: f64,
    alpha_abs: f64,
    rho0: f64,
    c0: f64,
) -> (f64, f64) {
    let p_sq = p0 * p0 * (-2.0 * alpha_abs * x).exp();
    let fx = p_sq * alpha_abs / (rho0 * c0 * c0);
    (fx, 0.0)
}

// ============================================================================
// 14. Nyborg Body Force
// ============================================================================

/// Nyborg body force for acoustic streaming (Nyborg 1958).
///
/// F = -(α_abs / ρ₀ c₀) × (P₀²/2) × exp(-2α_abs x) / c₀
///
/// Simplified steady streaming body force from the Nyborg formulation.
/// This approximates the time-averaged acoustic momentum flux gradient.
///
/// # Arguments
/// * `p0`       – acoustic pressure amplitude \[Pa\]
/// * `x`        – axial distance from source \[m\]
/// * `alpha_abs`– absorption coefficient \[1/m\]
/// * `rho0`     – fluid density \[kg/m³\]
/// * `c0`       – speed of sound \[m/s\]
pub fn nyborg_body_force(p0: f64, x: f64, alpha_abs: f64, rho0: f64, c0: f64) -> f64 {
    (alpha_abs * p0 * p0 * (-2.0 * alpha_abs * x).exp()) / (rho0 * c0 * c0)
}

/// Nyborg streaming velocity field (far-field Gaussian beam).
///
/// u_stream(x) = (α_abs P₀²)/(4 μ ρ₀ c₀²) × \[1 - exp(-2α_abs x)\]
///
/// Returns the streaming velocity at position x along the beam axis.
pub fn nyborg_streaming_velocity(
    p0: f64,
    x: f64,
    alpha_abs: f64,
    rho0: f64,
    c0: f64,
    mu: f64,
) -> f64 {
    let prefactor = alpha_abs * p0 * p0 / (4.0 * mu * rho0 * c0 * c0);
    prefactor * (1.0 - (-2.0 * alpha_abs * x).exp())
}

// ============================================================================
// 15. Full LBM Acoustic Streaming Simulation
// ============================================================================

/// State of a 2-D acoustic streaming LBM simulation.
///
/// Holds distribution functions, source parameters, and physical scales.
pub struct AcousticStreamingLbm {
    /// Grid x-size.
    pub nx: usize,
    /// Grid y-size.
    pub ny: usize,
    /// Distribution functions, length 9 * nx * ny.
    pub f: Vec<f64>,
    /// BGK relaxation time (τ).
    pub tau: f64,
    /// Background density ρ₀ (lattice units).
    pub rho0: f64,
    /// Acoustic absorption coefficient α_abs (lattice units).
    pub alpha_abs: f64,
    /// RMS pressure field (lattice units), length nx * ny.
    pub p_rms: Vec<f64>,
    /// Current simulation step.
    pub step: usize,
}

impl AcousticStreamingLbm {
    /// Create a new acoustic streaming LBM simulation.
    ///
    /// # Arguments
    /// * `nx`, `ny`   – grid dimensions
    /// * `tau`        – BGK relaxation time (> 0.5)
    /// * `rho0`       – background density
    /// * `alpha_abs`  – absorption coefficient
    pub fn new(nx: usize, ny: usize, tau: f64, rho0: f64, alpha_abs: f64) -> Self {
        assert!(tau > 0.5, "tau must be > 0.5 for stability");
        let n = nx * ny;
        let feq = d2q9_equilibrium(rho0, 0.0, 0.0);
        let mut f = vec![0.0_f64; 9 * n];
        for node in 0..n {
            for i in 0..9 {
                f[i + 9 * node] = feq[i];
            }
        }
        Self {
            nx,
            ny,
            f,
            tau,
            rho0,
            alpha_abs,
            p_rms: vec![0.0; n],
            step: 0,
        }
    }

    /// Set the RMS pressure field (precomputed from acoustic solver).
    pub fn set_rms_pressure(&mut self, p_rms: Vec<f64>) {
        assert_eq!(p_rms.len(), self.nx * self.ny);
        self.p_rms = p_rms;
    }

    /// Set a uniform RMS pressure field.
    pub fn set_uniform_pressure(&mut self, p_val: f64) {
        self.p_rms.iter_mut().for_each(|p| *p = p_val);
    }

    /// Perform one LBM step: collision + streaming with acoustic forcing.
    ///
    /// The Eckart body force is applied via the Guo scheme.
    pub fn step_with_eckart_force(&mut self, rho_phys: f64, c0_phys: f64) {
        let n = self.nx * self.ny;
        let mut f_coll = vec![0.0_f64; 9 * n];
        for node in 0..n {
            let fi: [f64; 9] = std::array::from_fn(|i| self.f[i + 9 * node]);
            let (rho, ux, uy) = d2q9_macroscopic(&fi);
            let fx = eckart_streaming_force(self.p_rms[node], self.alpha_abs, rho_phys, c0_phys);
            let f_new = bgk_collision_with_force(&fi, self.tau, rho, ux, uy, fx, 0.0, 1.0);
            for i in 0..9 {
                f_coll[i + 9 * node] = f_new[i];
            }
        }
        self.f = stream_d2q9(&f_coll, self.nx, self.ny);
        self.step += 1;
    }

    /// Extract the velocity field as (ux, uy) arrays.
    pub fn velocity_field(&self) -> (Vec<f64>, Vec<f64>) {
        let n = self.nx * self.ny;
        let mut ux_field = vec![0.0_f64; n];
        let mut uy_field = vec![0.0_f64; n];
        for node in 0..n {
            let fi: [f64; 9] = std::array::from_fn(|i| self.f[i + 9 * node]);
            let (_, ux, uy) = d2q9_macroscopic(&fi);
            ux_field[node] = ux;
            uy_field[node] = uy;
        }
        (ux_field, uy_field)
    }

    /// Compute the kinetic energy of the streaming flow.
    pub fn kinetic_energy(&self) -> f64 {
        let n = self.nx * self.ny;
        let mut ke = 0.0_f64;
        for node in 0..n {
            let fi: [f64; 9] = std::array::from_fn(|i| self.f[i + 9 * node]);
            let (rho, ux, uy) = d2q9_macroscopic(&fi);
            ke += 0.5 * rho * (ux * ux + uy * uy);
        }
        ke
    }

    /// Compute acoustic pressure field relative to background density.
    pub fn pressure_field(&self) -> Vec<f64> {
        acoustic_pressure_field(&self.f, self.nx, self.ny, self.rho0)
    }
}

// ============================================================================
// 16. Bounce-Back Boundary Condition
// ============================================================================

/// Apply bounce-back boundary condition at wall nodes.
///
/// At wall nodes (marked true in `is_wall`), populations are reflected:
/// f_i = f_{opp(i)} (no-slip wall).
///
/// Modifies distribution functions in-place.
///
/// # Arguments
/// * `f`       – distribution functions
/// * `is_wall` – boolean mask, length nx * ny
/// * `nx`, `ny` – grid dimensions
pub fn apply_bounce_back(f: &mut [f64], is_wall: &[bool], nx: usize, ny: usize) {
    let n = nx * ny;
    assert_eq!(f.len(), 9 * n);
    assert_eq!(is_wall.len(), n);
    for node in 0..n {
        if is_wall[node] {
            let fi: [f64; 9] = std::array::from_fn(|i| f[i + 9 * node]);
            for i in 0..9 {
                f[i + 9 * node] = fi[D2Q9_OPP[i]];
            }
        }
    }
}

// ============================================================================
// 17. Acoustic Tweezers
// ============================================================================

/// Acoustic tweezer gradient force on a Rayleigh particle (R ≪ λ).
///
/// F_grad = -∇U_gorkov where U is the Gorkov potential.
/// The gradient is computed numerically using central differences.
///
/// # Arguments
/// * `radius`    – particle radius \[m\]
/// * `x`, `y`   – particle position \[m\]
/// * `p_field`  – 2-D pressure amplitude field, length nx*ny
/// * `nx`, `ny` – grid dimensions
/// * `dx`, `dy` – grid spacings \[m\]
/// * `rho0`, `c0`, `rho_p`, `c_p` – material parameters
pub fn acoustic_tweezer_force(
    radius: f64,
    x: f64,
    y: f64,
    p_field: &[f64],
    nx: usize,
    ny: usize,
    dx: f64,
    dy: f64,
    rho0: f64,
    c0: f64,
    rho_p: f64,
    c_p: f64,
) -> (f64, f64) {
    // Bilinear interpolation helper
    let interpolate = |field: &[f64], px: f64, py: f64| -> f64 {
        let ix = (px / dx).floor() as isize;
        let iy = (py / dy).floor() as isize;
        let fx = px / dx - ix as f64;
        let fy = py / dy - iy as f64;
        let get = |i: isize, j: isize| -> f64 {
            let ci = i.rem_euclid(nx as isize) as usize;
            let cj = j.rem_euclid(ny as isize) as usize;
            field[cj * nx + ci]
        };
        (1.0 - fx) * (1.0 - fy) * get(ix, iy)
            + fx * (1.0 - fy) * get(ix + 1, iy)
            + (1.0 - fx) * fy * get(ix, iy + 1)
            + fx * fy * get(ix + 1, iy + 1)
    };

    let p_center = interpolate(p_field, x, y);
    let p_xp = interpolate(p_field, x + dx, y);
    let p_xm = interpolate(p_field, x - dx, y);
    let p_yp = interpolate(p_field, x, y + dy);
    let p_ym = interpolate(p_field, x, y - dy);

    let f1 = 1.0 - rho0 * c0 * c0 / (rho_p * c_p * c_p);
    let _f2 = 2.0 * (rho_p - rho0) / (2.0 * rho_p + rho0);
    let pre = 2.0 * PI * radius.powi(3) * f1 / (3.0 * rho0 * c0 * c0);

    // Gradient of p² = 2p dp/dx
    let dp_dx = (p_xp - p_xm) / (2.0 * dx);
    let dp_dy = (p_yp - p_ym) / (2.0 * dy);
    let fx = -pre * 2.0 * p_center * dp_dx;
    let fy = -pre * 2.0 * p_center * dp_dy;
    (fx, fy)
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_d2q9_equilibrium_zero_velocity() {
        // At rest: f_i = w_i * rho
        let rho = 1.2;
        let feq = d2q9_equilibrium(rho, 0.0, 0.0);
        for i in 0..9 {
            let expected = D2Q9_W[i] * rho;
            assert!(
                (feq[i] - expected).abs() < 1e-12,
                "f[{}]: {:.8} vs {:.8}",
                i,
                feq[i],
                expected
            );
        }
    }

    #[test]
    fn test_d2q9_equilibrium_mass_conservation() {
        // sum f_i^eq = rho
        let rho = 0.9;
        let feq = d2q9_equilibrium(rho, 0.1, -0.05);
        let sum: f64 = feq.iter().sum();
        assert!((sum - rho).abs() < 1e-12, "sum = {:.8}", sum);
    }

    #[test]
    fn test_d2q9_equilibrium_momentum_conservation() {
        // sum f_i^eq cx_i = rho * ux
        let rho = 1.0;
        let ux = 0.05;
        let uy = -0.03;
        let feq = d2q9_equilibrium(rho, ux, uy);
        let mom_x: f64 = (0..9).map(|i| feq[i] * D2Q9_CX[i]).sum();
        let mom_y: f64 = (0..9).map(|i| feq[i] * D2Q9_CY[i]).sum();
        assert!((mom_x - rho * ux).abs() < 1e-12, "mom_x = {:.8}", mom_x);
        assert!((mom_y - rho * uy).abs() < 1e-12, "mom_y = {:.8}", mom_y);
    }

    #[test]
    fn test_d2q9_macroscopic_roundtrip() {
        let rho = 1.05;
        let ux = 0.02;
        let uy = -0.01;
        let feq = d2q9_equilibrium(rho, ux, uy);
        let (rho_out, ux_out, uy_out) = d2q9_macroscopic(&feq);
        assert!((rho_out - rho).abs() < 1e-12);
        assert!((ux_out - ux).abs() < 1e-12);
        assert!((uy_out - uy).abs() < 1e-12);
    }

    #[test]
    fn test_bgk_collision_conserves_mass() {
        let rho = 1.0;
        let ux = 0.05;
        let uy = 0.02;
        let feq = d2q9_equilibrium(rho, ux, uy);
        // Perturb slightly
        let mut f = feq;
        f[0] += 0.01;
        f[1] -= 0.01;
        let f_out = bgk_collision(&f, 1.0, rho, ux, uy);
        let mass_in: f64 = f.iter().sum();
        let mass_out: f64 = f_out.iter().sum();
        assert!(
            (mass_in - mass_out).abs() < 1e-12,
            "mass: {:.8} vs {:.8}",
            mass_in,
            mass_out
        );
    }

    #[test]
    fn test_stream_d2q9_shape_preserved() {
        let nx = 4;
        let ny = 4;
        let n = 9 * nx * ny;
        let f: Vec<f64> = (0..n).map(|i| i as f64 * 0.1).collect();
        let f_out = stream_d2q9(&f, nx, ny);
        assert_eq!(f_out.len(), n);
        // Total mass conserved under periodic streaming
        let mass_in: f64 = f.chunks(9).map(|c| c.iter().sum::<f64>()).sum();
        let mass_out: f64 = f_out.chunks(9).map(|c| c.iter().sum::<f64>()).sum();
        assert!(
            (mass_in - mass_out).abs() < 1e-8,
            "mass: {:.8} vs {:.8}",
            mass_in,
            mass_out
        );
    }

    #[test]
    fn test_acoustic_pressure_field_zero_perturbation() {
        let nx = 4;
        let ny = 4;
        let rho0 = 1.0;
        let mut f = vec![0.0_f64; 9 * nx * ny];
        // Fill with equilibrium at rest
        let feq = d2q9_equilibrium(rho0, 0.0, 0.0);
        for node in 0..nx * ny {
            for i in 0..9 {
                f[i + 9 * node] = feq[i];
            }
        }
        let p = acoustic_pressure_field(&f, nx, ny, rho0);
        for &pi in &p {
            assert!(pi.abs() < 1e-12, "p = {:.6e}", pi);
        }
    }

    #[test]
    fn test_sound_pressure_level() {
        // SPL at reference pressure should be 0 dB
        let spl = sound_pressure_level(P_REF);
        assert!(spl.abs() < 1e-6, "SPL = {:.6}", spl);
        // Doubling pressure adds ~6 dB
        let spl2 = sound_pressure_level(2.0 * P_REF);
        assert!((spl2 - 20.0 * 2.0_f64.log10()).abs() < 1e-6);
    }

    #[test]
    fn test_standing_wave_at_node() {
        // At a pressure node x = λ/4, pressure = 0 at all times
        let omega = 2.0 * PI * 1000.0;
        let c0 = 343.0;
        let k = omega / c0;
        let x_node = PI / (2.0 * k); // = λ/4
        for t in [0.0, 1e-4, 2e-4] {
            let p = standing_wave_pressure_1d(x_node, t, 1.0, k, omega);
            assert!(p.abs() < 1e-10, "pressure at node: {:.6e}", p);
        }
    }

    #[test]
    fn test_eckart_streaming_force_scaling() {
        // Force should scale as p_rms^2
        let f1 = eckart_streaming_force(1.0, 10.0, RHO0_WATER, C0_WATER);
        let f2 = eckart_streaming_force(2.0, 10.0, RHO0_WATER, C0_WATER);
        assert!((f2 - 4.0 * f1).abs() < 1e-20, "f2/f1 = {:.6}", f2 / f1);
    }

    #[test]
    fn test_gorkov_potential_water_particle() {
        // Polystyrene in water: positive Gorkov potential
        let rho_ps = 1050.0; // kg/m3
        let c_ps = 2350.0; // m/s
        let r = 5e-6; // 5 μm radius
        let p_amp = 1e5; // 100 kPa
        let v_amp = p_amp / (RHO0_WATER * C0_WATER);
        let u = gorkov_potential(r, p_amp, v_amp, RHO0_WATER, C0_WATER, rho_ps, c_ps);
        assert!(u.is_finite(), "Gorkov potential = {:.6e}", u);
    }

    #[test]
    fn test_acoustic_contrast_factor_polystyrene() {
        // Polystyrene in water: Φ > 0, focuses at pressure nodes
        let phi = acoustic_contrast_factor(RHO0_WATER, C0_WATER, 1050.0, 2350.0);
        assert!(phi > 0.0, "contrast factor = {:.6}", phi);
    }

    #[test]
    fn test_acoustic_contrast_factor_lipid() {
        // Lipid particles: Φ < 0, focuses at pressure antinodes
        let phi = acoustic_contrast_factor(RHO0_WATER, C0_WATER, 900.0, 1450.0);
        assert!(phi < 0.0, "contrast factor = {:.6}", phi);
    }

    #[test]
    fn test_piezo_resonant_frequency() {
        // Simple LC resonance check
        let l = 1e-6; // 1 μH
        let c = 1e-9; // 1 nF
        let f_r = piezo_resonant_frequency(l, c);
        let expected = 1.0 / (2.0 * PI * (l * c).sqrt());
        assert!((f_r - expected).abs() < 1e-3);
    }

    #[test]
    fn test_piezo_pressure_source_at_zero() {
        // At t=0 with rise time > 0, envelope = 0, so source = 0
        let p = piezo_pressure_source(0.0, 1e5, 1e6, 1e-6);
        assert!(p.abs() < 1e-10, "p(0) = {:.6e}", p);
    }

    #[test]
    fn test_rayleigh_streaming_sign() {
        // At x = λ/8 (sin(2kx) > 0), streaming velocity should be negative
        let omega = 2.0 * PI * 1e6;
        let c0 = C0_WATER;
        let k = omega / c0;
        let d = 100e-6; // 100 μm channel
        let u0 = 0.1; // m/s acoustic velocity amplitude
        let x = PI / (8.0 * k); // sin(2kx) = sin(π/4) > 0
        let u = rayleigh_streaming_velocity(d / 2.0, d, x, u0, k, c0);
        assert!(u < 0.0, "Rayleigh velocity = {:.6e}", u);
    }

    #[test]
    fn test_nyborg_body_force_decay() {
        // Body force should decrease with distance
        let p0 = 1e5;
        let alpha = 0.1;
        let f1 = nyborg_body_force(p0, 0.0, alpha, RHO0_WATER, C0_WATER);
        let f2 = nyborg_body_force(p0, 1.0, alpha, RHO0_WATER, C0_WATER);
        assert!(
            f1 > f2,
            "force at x=0 ({:.6e}) should exceed x=1 ({:.6e})",
            f1,
            f2
        );
    }

    #[test]
    fn test_nyborg_streaming_velocity_saturation() {
        // Velocity should saturate (plateau) for large x
        let p0 = 1e5;
        let alpha = 10.0; // high absorption → rapid saturation
        let v_far = nyborg_streaming_velocity(p0, 1.0, alpha, RHO0_WATER, C0_WATER, MU_WATER);
        let v_very_far =
            nyborg_streaming_velocity(p0, 100.0, alpha, RHO0_WATER, C0_WATER, MU_WATER);
        assert!(
            (v_far - v_very_far).abs() / v_very_far < 1e-3,
            "v_far={:.6e}, v_very_far={:.6e}",
            v_far,
            v_very_far
        );
    }

    #[test]
    fn test_acoustic_streaming_lbm_init() {
        let sim = AcousticStreamingLbm::new(8, 8, 1.0, 1.0, 0.01);
        assert_eq!(sim.f.len(), 9 * 64);
        assert_eq!(sim.step, 0);
        // Initial kinetic energy should be near zero
        let ke = sim.kinetic_energy();
        assert!(ke < 1e-20, "KE = {:.6e}", ke);
    }

    #[test]
    fn test_acoustic_streaming_lbm_step_runs() {
        let mut sim = AcousticStreamingLbm::new(8, 8, 1.0, 1.0, 0.01);
        sim.set_uniform_pressure(0.1);
        sim.step_with_eckart_force(RHO0_WATER, C0_WATER);
        assert_eq!(sim.step, 1);
        // After one step with forcing, velocity should be non-zero somewhere
        let (ux, _uy) = sim.velocity_field();
        let max_ux = ux.iter().cloned().fold(0.0_f64, f64::max);
        assert!(max_ux > 0.0, "max ux = {:.6e}", max_ux);
    }

    #[test]
    fn test_acoustic_streaming_lbm_mass_conservation() {
        let mut sim = AcousticStreamingLbm::new(6, 6, 0.8, 1.0, 0.05);
        sim.set_uniform_pressure(0.05);
        let mass_initial: f64 = sim.f.chunks(9).map(|c| c.iter().sum::<f64>()).sum();
        for _ in 0..5 {
            sim.step_with_eckart_force(RHO0_WATER, C0_WATER);
        }
        let mass_final: f64 = sim.f.chunks(9).map(|c| c.iter().sum::<f64>()).sum();
        let rel_err = (mass_final - mass_initial).abs() / mass_initial;
        assert!(rel_err < 1e-6, "mass error = {:.6e}", rel_err);
    }

    #[test]
    fn test_bounce_back_reverses_populations() {
        let nx = 3;
        let ny = 1;
        let mut f = vec![0.0_f64; 9 * nx * ny];
        let feq = d2q9_equilibrium(1.0, 0.1, 0.0);
        f[9..18].copy_from_slice(&feq); // middle node
        let is_wall = vec![false, true, false];
        apply_bounce_back(&mut f, &is_wall, nx, ny);
        // At wall node, f[i] = f[opp(i)] (after bounce)
        for i in 0..9 {
            let val = f[i + 9];
            let expected = feq[D2Q9_OPP[i]];
            assert!(
                (val - expected).abs() < 1e-12,
                "f[{}] = {:.8} vs {:.8}",
                i,
                val,
                expected
            );
        }
    }

    #[test]
    fn test_levitation_force_spherical() {
        // F_lev = (4/3) π R³ ρ_p g
        let r = 1e-3;
        let rho_p = 2700.0;
        let g = 9.81;
        let f_lev = levitation_force_required(r, rho_p, g);
        let expected = (4.0 / 3.0) * PI * r.powi(3) * rho_p * g;
        assert!((f_lev - expected).abs() < 1e-20);
    }

    #[test]
    fn test_levitation_equilibrium_height() {
        // Height should be λ/4
        let c0 = C0_AIR;
        let freq = 40e3; // 40 kHz
        let h = levitation_equilibrium_height(c0, freq);
        let lambda = c0 / freq;
        assert!((h - lambda / 4.0).abs() < 1e-12);
    }

    #[test]
    fn test_acoustic_reynolds_stress_direction() {
        // For a wave along x (theta=0): Txx > 0, Tyy = 0
        let t = acoustic_reynolds_stress(0.1, RHO0_WATER, 0.0);
        assert!(t[0][0] > 0.0);
        assert!(t[1][1].abs() < 1e-20);
    }

    #[test]
    fn test_reynolds_stress_body_force_direction() {
        // Body force should be in +x for attenuated beam
        let (fx, fy) = reynolds_stress_body_force(0.1, 0.0, 1e5, 0.1, RHO0_WATER, C0_WATER);
        assert!(fx > 0.0, "fx = {:.6e}", fx);
        assert_eq!(fy, 0.0);
    }

    #[test]
    fn test_particle_focusing_trajectory_converges() {
        // Positive contrast factor: particle should move toward λ/4 node
        let omega = 2.0 * PI * 1e6;
        let c0 = C0_WATER;
        let k = omega / c0;
        let lambda4 = PI / (2.0 * k); // pressure node at λ/4
        let x0 = lambda4 * 1.5; // start off-center (between node and antinode)
        let radius = 5e-6;
        let p0 = 1e5;
        let rho_p = 1050.0;
        let c_p = 2350.0;
        let mu = MU_WATER;
        let t_end = 1e-3;
        let n_steps = 100;
        let (_t, x) = particle_focusing_trajectory(
            x0, radius, k, p0, RHO0_WATER, c0, rho_p, c_p, mu, t_end, n_steps,
        );
        // Particle should move (not stay stationary)
        let dx = (x.last().unwrap() - x0).abs();
        assert!(dx > 0.0, "particle did not move");
    }

    #[test]
    fn test_rms_pressure_uniform() {
        let p = vec![1.0, 1.0, 1.0, 1.0];
        assert!((rms_pressure(&p) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_acoustic_energy_density() {
        let e = acoustic_energy_density(1.0, 1.0, 1.0);
        assert!((e - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_standing_wave_nodes_count() {
        // 2nd mode: 1 node between the walls
        let nodes = standing_wave_nodes(1.0, 2);
        assert_eq!(nodes.len(), 1);
        assert!((nodes[0] - 0.5).abs() < 1e-12);
    }

    #[test]
    fn test_half_wavelength_focus_center() {
        let w = 200e-6;
        let pos = half_wavelength_focus_position(w);
        assert!((pos - w / 2.0).abs() < 1e-12);
    }

    #[test]
    fn test_d2q9_equilibrium_acoustic_sums_to_rho_prime() {
        let rho_prime = 0.1;
        let rho0 = 1.0;
        let ux = 0.02;
        let uy = -0.01;
        let feq = d2q9_equilibrium_acoustic(rho_prime, rho0, ux, uy);
        let sum: f64 = feq.iter().sum();
        assert!((sum - rho_prime).abs() < 1e-12, "sum = {:.8}", sum);
    }
}
