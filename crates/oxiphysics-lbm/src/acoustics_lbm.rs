#![allow(
    clippy::needless_range_loop,
    clippy::ptr_arg,
    clippy::too_many_arguments
)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Acoustic Lattice Boltzmann Method module.
//!
//! Implements acoustic wave simulation on a lattice, including:
//!
//! - **Linearized Euler equations** discretized on a lattice
//! - **Sound pressure level** computation (SPL in dB)
//! - **Acoustic source terms**: monopole, dipole, quadrupole
//! - **Acoustic wave propagation** on 2-D D2Q9 lattice
//! - **Perfectly Matched Layer (PML)** absorbing boundaries
//! - **Acoustic scattering** coefficient computations
//! - **Aero-acoustics** via Lighthill's equation discretized
//! - **Room acoustics**: image-source reflection / diffraction (simplified)
//! - **Acoustic streaming** body force
//! - **Acoustic radiation pressure** (Langevin / Gorkov)

#![allow(dead_code)]

use std::f64::consts::PI;

// ============================================================================
// Constants
// ============================================================================

/// Reference pressure for SPL computation (20 µPa in air).
pub const P_REF: f64 = 2.0e-5;

/// Reference sound speed in air at 20 °C \[m/s\].
pub const C0_AIR: f64 = 343.0;

/// Reference air density at 20 °C \[kg/m³\].
pub const RHO0_AIR: f64 = 1.2041;

/// D2Q9 lattice speed of sound squared: cs² = 1/3.
pub const CS2: f64 = 1.0_f64 / 3.0_f64;

/// D2Q9 lattice speed of sound.
pub const CS: f64 = 0.577_350_269_189_625_8_f64; // sqrt(1/3)

// ============================================================================
// 1. Linearized Euler Equations on Lattice
// ============================================================================

/// State vector for the linearized Euler acoustic LBM: (ρ', u_x', u_y').
///
/// ρ' is the density perturbation, u' are velocity perturbations around
/// the mean flow (u0, v0).
#[allow(dead_code)]
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AcousticState {
    /// Density perturbation ρ'.
    pub rho_prime: f64,
    /// Velocity perturbation in x-direction.
    pub ux_prime: f64,
    /// Velocity perturbation in y-direction.
    pub uy_prime: f64,
}

impl AcousticState {
    /// Create a new acoustic state with given perturbations.
    pub fn new(rho_prime: f64, ux_prime: f64, uy_prime: f64) -> Self {
        Self {
            rho_prime,
            ux_prime,
            uy_prime,
        }
    }

    /// Compute the acoustic pressure perturbation: p' = cs² * ρ'.
    pub fn pressure_perturbation(&self) -> f64 {
        CS2 * self.rho_prime
    }

    /// Compute the acoustic intensity vector I = p' * u'.
    pub fn intensity(&self) -> [f64; 2] {
        let p = self.pressure_perturbation();
        [p * self.ux_prime, p * self.uy_prime]
    }
}

/// D2Q9 equilibrium distribution for the linearized acoustic equations.
///
/// Uses the second-order expansion:
/// f_eq_i = w_i * ρ' * (1 + e_i·u/cs² + (e_i·u)²/(2cs⁴) - u²/(2cs²))
/// where u is the total velocity (mean + perturbation).
///
/// Returns array of 9 distribution function values (D2Q9 ordering).
#[allow(dead_code)]
pub fn acoustic_equilibrium(state: &AcousticState, u0: f64, v0: f64) -> [f64; 9] {
    // D2Q9 weights
    let weights = [
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
    // D2Q9 velocity vectors (ex, ey)
    let ex = [0.0, 1.0, 0.0, -1.0, 0.0, 1.0, -1.0, -1.0, 1.0];
    let ey = [0.0, 0.0, 1.0, 0.0, -1.0, 1.0, 1.0, -1.0, -1.0];

    let ux = u0 + state.ux_prime;
    let uy = v0 + state.uy_prime;
    let u2 = ux * ux + uy * uy;
    let rho = state.rho_prime;

    let mut feq = [0.0f64; 9];
    for i in 0..9 {
        let eu = ex[i] * ux + ey[i] * uy;
        feq[i] =
            weights[i] * rho * (1.0 + eu / CS2 + eu * eu / (2.0 * CS2 * CS2) - u2 / (2.0 * CS2));
    }
    feq
}

/// Perform a single BGK collision step for acoustic LBM.
///
/// `f` is the current distribution, `feq` the equilibrium, `omega` the
/// relaxation frequency (omega = 1/tau).  Returns updated distributions.
#[allow(dead_code)]
pub fn acoustic_bgk_collision(f: &[f64; 9], feq: &[f64; 9], omega: f64) -> [f64; 9] {
    let mut f_out = *f;
    for i in 0..9 {
        f_out[i] = f[i] - omega * (f[i] - feq[i]);
    }
    f_out
}

/// Extract macroscopic density perturbation and velocity perturbation from
/// distribution functions.
///
/// Returns `(rho_prime, ux_prime, uy_prime)`.
#[allow(dead_code)]
pub fn acoustic_moments(f: &[f64; 9]) -> (f64, f64, f64) {
    let ex = [0.0, 1.0, 0.0, -1.0, 0.0, 1.0, -1.0, -1.0, 1.0];
    let ey = [0.0, 0.0, 1.0, 0.0, -1.0, 1.0, 1.0, -1.0, -1.0];

    let rho: f64 = f.iter().sum();
    let ux: f64 = f.iter().enumerate().map(|(i, &fi)| fi * ex[i]).sum();
    let uy: f64 = f.iter().enumerate().map(|(i, &fi)| fi * ey[i]).sum();

    if rho.abs() > 1.0e-30 {
        (rho, ux / rho, uy / rho)
    } else {
        (rho, 0.0, 0.0)
    }
}

// ============================================================================
// 2. Sound Pressure Level
// ============================================================================

/// Compute Sound Pressure Level (SPL) in decibels from an RMS pressure.
///
/// SPL = 20 * log10(p_rms / p_ref), where p_ref = 20 µPa.
#[allow(dead_code)]
pub fn sound_pressure_level(p_rms: f64) -> f64 {
    20.0 * (p_rms.abs() / P_REF).log10()
}

/// Compute RMS pressure from a time series of pressure values.
#[allow(dead_code)]
pub fn rms_pressure(p_series: &[f64]) -> f64 {
    if p_series.is_empty() {
        return 0.0;
    }
    let mean_sq: f64 = p_series.iter().map(|&p| p * p).sum::<f64>() / p_series.len() as f64;
    mean_sq.sqrt()
}

/// Compute SPL directly from a pressure time series.
#[allow(dead_code)]
pub fn spl_from_series(p_series: &[f64]) -> f64 {
    sound_pressure_level(rms_pressure(p_series))
}

/// Weighted A-filter gain at a given frequency \[Hz\].
///
/// Uses the standard IEC 61672 A-weighting formula.
#[allow(dead_code)]
pub fn a_weighting_db(freq_hz: f64) -> f64 {
    let f2 = freq_hz * freq_hz;
    let f4 = f2 * f2;
    let numerator = 12194.0_f64.powi(2) * f4;
    let d1 = f2 + 20.6_f64.powi(2);
    let d2 = (f2 + 107.7_f64.powi(2)).sqrt() * (f2 + 737.9_f64.powi(2)).sqrt();
    let d3 = f2 + 12194.0_f64.powi(2);
    let ra = numerator / (d1 * d2 * d3);
    20.0 * ra.log10() + 2.0
}

// ============================================================================
// 3. Acoustic Source Terms
// ============================================================================

/// Monopole acoustic source term added to the density equation.
///
/// Represents a point source with volume velocity Q(t) at position (xs, ys).
/// Returns the source contribution to density at field point (x, y).
#[allow(dead_code)]
pub fn monopole_source_term(q: f64, xs: f64, ys: f64, x: f64, y: f64, c0: f64, rho0: f64) -> f64 {
    let r = ((x - xs).powi(2) + (y - ys).powi(2)).sqrt();
    if r < 1.0e-10 {
        return 0.0;
    }
    rho0 * q / (2.0 * PI * c0 * r)
}

/// Dipole acoustic source: modeled as two monopoles of opposite sign separated
/// by distance `d` along unit direction `dir`.
///
/// Returns pressure perturbation contribution at field point.
#[allow(dead_code)]
#[allow(clippy::too_many_arguments)]
pub fn dipole_source_pressure(
    force: [f64; 2],
    xs: f64,
    ys: f64,
    x: f64,
    y: f64,
    freq: f64,
    c0: f64,
    rho0: f64,
) -> f64 {
    let k = 2.0 * PI * freq / c0;
    let rx = x - xs;
    let ry = y - ys;
    let r = (rx * rx + ry * ry).sqrt();
    if r < 1.0e-10 {
        return 0.0;
    }
    // cos(theta) = F . r_hat / |F|
    let f_dot_r = force[0] * rx / r + force[1] * ry / r;
    let f_mag = (force[0] * force[0] + force[1] * force[1]).sqrt();
    if f_mag < 1.0e-30 {
        return 0.0;
    }
    // Far-field dipole: i*k / (4*pi*r) * F.r_hat * exp(-ikr)
    rho0 * k * f_dot_r / (4.0 * PI * r) * (k * r).cos()
}

/// Quadrupole source pressure from the T_ij stress tensor.
///
/// Uses the Lighthill quadrupole radiation formula in 2-D far field.
#[allow(dead_code)]
#[allow(clippy::too_many_arguments)]
pub fn quadrupole_source_pressure(
    t: [[f64; 2]; 2],
    xs: f64,
    ys: f64,
    x: f64,
    y: f64,
    freq: f64,
    c0: f64,
    rho0: f64,
) -> f64 {
    let k = 2.0 * PI * freq / c0;
    let rx = x - xs;
    let ry = y - ys;
    let r = (rx * rx + ry * ry).sqrt();
    if r < 1.0e-10 {
        return 0.0;
    }
    let rx_n = rx / r;
    let ry_n = ry / r;
    // ri * rj * T_ij / r^2 (double contraction in observer direction)
    let mut tij_rirj = 0.0;
    let r_hat = [rx_n, ry_n];
    for i in 0..2 {
        for j in 0..2 {
            tij_rirj += t[i][j] * r_hat[i] * r_hat[j];
        }
    }
    rho0 * k * k * tij_rirj / (4.0 * PI * r) * (k * r).cos()
}

/// Add a monopole source to the LBM forcing by injecting mass into f\[0\].
///
/// `source_strength` is Δρ per time step from the monopole.
/// Returns the modified distribution function.
#[allow(dead_code)]
pub fn inject_monopole_lbm(f: &[f64; 9], source_strength: f64) -> [f64; 9] {
    let mut f_out = *f;
    // Distribute source uniformly with D2Q9 weights (rest weight 4/9)
    let weights = [
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
    for i in 0..9 {
        f_out[i] += weights[i] * source_strength;
    }
    f_out
}

// ============================================================================
// 4. Acoustic Wave Propagation (D2Q9 grid)
// ============================================================================

/// A 2-D acoustic LBM grid using D2Q9.
///
/// Stores distribution functions at each grid point.
#[allow(dead_code)]
pub struct AcousticGrid2D {
    /// Number of grid cells in x.
    pub nx: usize,
    /// Number of grid cells in y.
    pub ny: usize,
    /// Distribution functions f\[i\]\[j\]\[q\].
    pub f: Vec<Vec<[f64; 9]>>,
    /// Post-collision distributions.
    pub f_star: Vec<Vec<[f64; 9]>>,
    /// BGK relaxation frequency omega = 1/tau.
    pub omega: f64,
    /// Mean background flow velocity in x.
    pub u0: f64,
    /// Mean background flow velocity in y.
    pub v0: f64,
}

impl AcousticGrid2D {
    /// Create a new acoustic 2-D grid of size nx × ny.
    ///
    /// `omega` is the relaxation parameter (1/tau), `u0` and `v0` are the
    /// mean flow velocities (for convected acoustics).
    pub fn new(nx: usize, ny: usize, omega: f64, u0: f64, v0: f64) -> Self {
        Self {
            nx,
            ny,
            f: vec![vec![[0.0; 9]; ny]; nx],
            f_star: vec![vec![[0.0; 9]; ny]; nx],
            omega,
            u0,
            v0,
        }
    }

    /// Initialize the grid with a Gaussian pressure pulse centered at (xc, yc).
    ///
    /// `amp` is the pulse amplitude, `sigma` the Gaussian width.
    pub fn init_gaussian_pulse(&mut self, xc: f64, yc: f64, amp: f64, sigma: f64) {
        for i in 0..self.nx {
            for j in 0..self.ny {
                let r2 = (i as f64 - xc).powi(2) + (j as f64 - yc).powi(2);
                let rho_prime = amp * (-r2 / (2.0 * sigma * sigma)).exp();
                let state = AcousticState::new(rho_prime, 0.0, 0.0);
                self.f[i][j] = acoustic_equilibrium(&state, self.u0, self.v0);
            }
        }
    }

    /// Perform collision step on the entire grid.
    pub fn collide(&mut self) {
        for i in 0..self.nx {
            for j in 0..self.ny {
                let (rho_p, ux_p, uy_p) = acoustic_moments(&self.f[i][j]);
                let state = AcousticState::new(rho_p, ux_p, uy_p);
                let feq = acoustic_equilibrium(&state, self.u0, self.v0);
                self.f_star[i][j] = acoustic_bgk_collision(&self.f[i][j], &feq, self.omega);
            }
        }
    }

    /// Perform streaming step with periodic boundary conditions.
    pub fn stream_periodic(&mut self) {
        let ex: [i64; 9] = [0, 1, 0, -1, 0, 1, -1, -1, 1];
        let ey: [i64; 9] = [0, 0, 1, 0, -1, 1, 1, -1, -1];
        let nx = self.nx as i64;
        let ny = self.ny as i64;
        let mut f_new = vec![vec![[0.0f64; 9]; self.ny]; self.nx];
        for i in 0..self.nx {
            for j in 0..self.ny {
                for q in 0..9 {
                    let src_i = ((i as i64 - ex[q]).rem_euclid(nx)) as usize;
                    let src_j = ((j as i64 - ey[q]).rem_euclid(ny)) as usize;
                    f_new[i][j][q] = self.f_star[src_i][src_j][q];
                }
            }
        }
        self.f = f_new;
    }

    /// Extract the pressure perturbation field as a flat Vec (row-major).
    pub fn pressure_field(&self) -> Vec<f64> {
        let mut p = Vec::with_capacity(self.nx * self.ny);
        for i in 0..self.nx {
            for j in 0..self.ny {
                let (rho_p, _ux, _uy) = acoustic_moments(&self.f[i][j]);
                p.push(CS2 * rho_p);
            }
        }
        p
    }

    /// Extract the density perturbation at grid cell (i, j).
    pub fn rho_at(&self, i: usize, j: usize) -> f64 {
        acoustic_moments(&self.f[i][j]).0
    }
}

// ============================================================================
// 5. Perfectly Matched Layer (PML) Absorbing Boundary
// ============================================================================

/// Parameters for a PML absorbing layer.
///
/// The PML attenuates outgoing waves without spurious reflections.
#[allow(dead_code)]
#[derive(Debug, Clone, Copy)]
pub struct PmlParams {
    /// Thickness of the PML layer in grid cells.
    pub thickness: usize,
    /// Maximum absorption coefficient σ_max.
    pub sigma_max: f64,
    /// Polynomial grading order (typically 2 or 3).
    pub grading_order: f64,
}

impl PmlParams {
    /// Create a new PML parameter set.
    pub fn new(thickness: usize, sigma_max: f64, grading_order: f64) -> Self {
        Self {
            thickness,
            sigma_max,
            grading_order,
        }
    }

    /// Compute the absorption coefficient σ at distance `d` from PML start.
    ///
    /// Uses polynomial grading: σ(d) = σ_max * (d/L)^m.
    pub fn sigma(&self, d: f64) -> f64 {
        let l = self.thickness as f64;
        if d <= 0.0 {
            0.0
        } else if d >= l {
            self.sigma_max
        } else {
            self.sigma_max * (d / l).powf(self.grading_order)
        }
    }

    /// Compute the PML damping factor exp(-σ * dt) for a time step dt.
    pub fn damping_factor(&self, d: f64, dt: f64) -> f64 {
        (-self.sigma(d) * dt).exp()
    }
}

/// Apply PML damping to the acoustic pressure field in-place.
///
/// Damps cells within `pml_thickness` cells of each boundary.
#[allow(dead_code)]
pub fn apply_pml_2d(f: &mut Vec<Vec<[f64; 9]>>, nx: usize, ny: usize, pml: &PmlParams, dt: f64) {
    for i in 0..nx {
        for j in 0..ny {
            let d_left = i as f64;
            let d_right = (nx - 1 - i) as f64;
            let d_bottom = j as f64;
            let d_top = (ny - 1 - j) as f64;
            let d_min = d_left.min(d_right).min(d_bottom).min(d_top);

            // PML region: within `thickness` cells of any boundary
            let d_from_pml_start = (pml.thickness as f64 - d_min).max(0.0);
            let damp = pml.damping_factor(d_from_pml_start, dt);
            if damp < 1.0 {
                for q in 0..9 {
                    f[i][j][q] *= damp;
                }
            }
        }
    }
}

/// Compute the optimal σ_max for a PML of given thickness such that the
/// theoretical reflection coefficient R < `target_r`.
///
/// Uses: R = exp(-2 * σ_max * L / (c0 * (m+1))) => σ_max = -ln(R)*(m+1)*c0/(2L).
#[allow(dead_code)]
pub fn optimal_sigma_max(thickness: f64, c0: f64, grading_order: f64, target_r: f64) -> f64 {
    -(target_r.ln()) * (grading_order + 1.0) * c0 / (2.0 * thickness)
}

// ============================================================================
// 6. Acoustic Scattering
// ============================================================================

/// Compute the scattering cross-section for a rigid cylinder of radius `a`
/// at wavenumber `k` (2-D, leading-order Born approximation).
///
/// σ_scat ≈ π * k * a² (small ka limit).
#[allow(dead_code)]
pub fn scattering_cross_section_cylinder_2d(k: f64, a: f64) -> f64 {
    PI * k * a * a
}

/// Compute the acoustic scattered pressure amplitude from a rigid sphere of
/// radius `a` at distance `r` and wavenumber `k` (3-D far field, leading order).
///
/// p_scat ≈ (k² a³) / (3 r) * p_inc at backscatter.
#[allow(dead_code)]
pub fn sphere_scattered_pressure(k: f64, a: f64, r: f64, p_inc: f64) -> f64 {
    if r < 1.0e-10 {
        return 0.0;
    }
    (k * k * a * a * a) / (3.0 * r) * p_inc
}

/// T-matrix element for a single D2Q9 scatterer node (bounce-back obstacle).
///
/// Returns the reflected distribution index for D2Q9 bounce-back.
/// The bounce-back rule maps direction q -> opposite direction.
#[allow(dead_code)]
pub fn d2q9_bounce_back_index(q: usize) -> usize {
    // D2Q9 opposite indices: 0<->0, 1<->3, 2<->4, 5<->7, 6<->8
    const OPPOSITE: [usize; 9] = [0, 3, 4, 1, 2, 7, 8, 5, 6];
    OPPOSITE[q]
}

/// Apply rigid-wall (bounce-back) boundary condition to a distribution.
///
/// Swaps post-streaming populations with their opposite directions.
#[allow(dead_code)]
pub fn apply_bounce_back(f: &mut [f64; 9]) {
    let f_copy = *f;
    for q in 0..9 {
        f[q] = f_copy[d2q9_bounce_back_index(q)];
    }
}

/// Compute the diffraction attenuation using Maekawa's formula.
///
/// `n_delta` is the Fresnel number N = 2δ/λ (path length difference over wavelength).
/// Returns attenuation in dB.
#[allow(dead_code)]
pub fn maekawa_diffraction_db(n_delta: f64) -> f64 {
    if n_delta <= -0.2 {
        0.0
    } else {
        5.0 + 20.0 * (3_f64.sqrt() * n_delta).log10().max(0.0)
    }
}

// ============================================================================
// 7. Aero-Acoustics: Lighthill's Equation (Discretized)
// ============================================================================

/// Lighthill stress tensor T_ij = ρ u_i u_j + (p - c0² ρ) δ_ij.
///
/// Returns the 2×2 tensor for 2-D flow.
#[allow(dead_code)]
pub fn lighthill_tensor_2d(rho: f64, u: [f64; 2], p: f64, c0_sq: f64) -> [[f64; 2]; 2] {
    let excess = p - c0_sq * rho;
    [
        [rho * u[0] * u[0] + excess, rho * u[0] * u[1]],
        [rho * u[1] * u[0], rho * u[1] * u[1] + excess],
    ]
}

/// Compute the divergence of the Lighthill tensor on a uniform 2-D grid.
///
/// Returns the source term vector \[∂T_xi/∂x_j\] for row `i=0` (x-component).
/// Uses second-order central differences with grid spacing `dx`.
#[allow(dead_code)]
pub fn lighthill_divergence_x(
    t: &[Vec<[[f64; 2]; 2]>],
    i: usize,
    j: usize,
    nx: usize,
    ny: usize,
    dx: f64,
) -> f64 {
    let im = if i == 0 { nx - 1 } else { i - 1 };
    let ip = if i == nx - 1 { 0 } else { i + 1 };
    let jm = if j == 0 { ny - 1 } else { j - 1 };
    let jp = if j == ny - 1 { 0 } else { j + 1 };

    let dt_xx = (t[ip][j][0][0] - t[im][j][0][0]) / (2.0 * dx);
    let dt_xy = (t[i][jp][0][1] - t[i][jm][0][1]) / (2.0 * dx);
    dt_xx + dt_xy
}

/// Compute the Lighthill source power spectral density estimate (simplified).
///
/// Integrates T_ij² over the source volume. Returns total acoustic power ∝ ρ c0^-5 U^8.
#[allow(dead_code)]
pub fn lighthill_acoustic_power(
    t_field: &[Vec<[[f64; 2]; 2]>],
    rho0: f64,
    c0: f64,
    dx: f64,
) -> f64 {
    let mut sum = 0.0;
    for row in t_field {
        for t in row {
            for i in 0..2 {
                for j in 0..2 {
                    sum += t[i][j] * t[i][j];
                }
            }
        }
    }
    sum * dx * dx * rho0 / (c0.powi(5))
}

// ============================================================================
// 8. Room Acoustics (Reflection / Diffraction — Simplified)
// ============================================================================

/// Room geometry for simplified image-source method.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct RoomGeometry {
    /// Room width in x \[m\].
    pub lx: f64,
    /// Room width in y \[m\].
    pub ly: f64,
    /// Wall absorption coefficient (0 = perfectly reflective, 1 = fully absorbing).
    pub absorption: f64,
}

impl RoomGeometry {
    /// Create a new rectangular room.
    pub fn new(lx: f64, ly: f64, absorption: f64) -> Self {
        Self { lx, ly, absorption }
    }

    /// Compute reverberation time T60 via Sabine's formula.
    ///
    /// T60 = 0.161 * V / (α * S) where V = volume (lx*ly*1), S = surface area.
    pub fn t60_sabine(&self) -> f64 {
        let volume = self.lx * self.ly; // 2-D: area as volume proxy
        let surface = 2.0 * (self.lx + self.ly);
        0.161 * volume / (self.absorption * surface)
    }

    /// Generate image source positions for first-order reflections.
    ///
    /// Returns up to 4 image positions (one per wall) as `Vec<(f64,f64)>`.
    pub fn image_sources_first_order(&self, sx: f64, sy: f64) -> Vec<(f64, f64)> {
        vec![
            (-sx, sy),                // reflect through x=0 wall
            (2.0 * self.lx - sx, sy), // reflect through x=Lx wall
            (sx, -sy),                // reflect through y=0 wall
            (sx, 2.0 * self.ly - sy), // reflect through y=Ly wall
        ]
    }

    /// Compute total SPL at receiver (rx, ry) from source (sx, sy) including
    /// first-order image sources.
    ///
    /// Uses 1/r amplitude attenuation and applies absorption coefficient.
    pub fn room_spl(&self, sx: f64, sy: f64, rx: f64, ry: f64, p0: f64, freq: f64, c0: f64) -> f64 {
        let k = 2.0 * PI * freq / c0;
        // Direct path
        let r0 = ((rx - sx).powi(2) + (ry - sy).powi(2)).sqrt().max(1.0e-6);
        let mut p_total_sq = (p0 / r0).powi(2);

        // First-order images
        let images = self.image_sources_first_order(sx, sy);
        for (ix, iy) in &images {
            let r = ((rx - ix).powi(2) + (ry - iy).powi(2)).sqrt().max(1.0e-6);
            let p_img = p0 * (1.0 - self.absorption) / r * (k * r).cos();
            p_total_sq += p_img * p_img;
        }
        sound_pressure_level(p_total_sq.sqrt())
    }
}

/// Compute edge diffraction coefficient using the Uniform Theory of Diffraction
/// (UTD) half-plane formula (simplified, 2-D).
///
/// `theta_i` is incident angle, `theta_d` diffraction angle, both from edge normal.
#[allow(dead_code)]
pub fn utd_diffraction_coefficient(theta_i: f64, theta_d: f64, k: f64, l: f64) -> f64 {
    // Simplified UTD: D = -1 / (2 * sqrt(2*pi*k) * cos((theta_d - theta_i)/2))
    let denom = 2.0 * (2.0 * PI * k * l).sqrt() * ((theta_d - theta_i) / 2.0).cos();
    if denom.abs() < 1.0e-10 {
        0.0
    } else {
        -1.0 / denom
    }
}

// ============================================================================
// 9. Acoustic Streaming
// ============================================================================

/// Compute the acoustic streaming body force (Eckart streaming).
///
/// F_stream = (2 * α * I) / c0  where α is the absorption coefficient,
/// I is the acoustic intensity.  Returns force vector \[Fx, Fy\].
#[allow(dead_code)]
pub fn eckart_streaming_force(intensity: [f64; 2], alpha_abs: f64, c0: f64) -> [f64; 2] {
    [
        2.0 * alpha_abs * intensity[0] / c0,
        2.0 * alpha_abs * intensity[1] / c0,
    ]
}

/// Compute the Rayleigh streaming velocity estimate between two parallel walls.
///
/// U_stream ~ (3/8) * (v_0² / c0) * k * sin(2kx) * sinh(2αy) / sinh(2αH)
/// Simplified to peak amplitude: U_max = (3 * v0² * k) / (8 * c0).
#[allow(dead_code)]
pub fn rayleigh_streaming_peak_velocity(v0: f64, k: f64, c0: f64) -> f64 {
    3.0 * v0 * v0 * k / (8.0 * c0)
}

/// Add acoustic streaming force to LBM distributions using Guo's forcing scheme.
///
/// Guo forcing: f_i += w_i * (e_i - u) / cs² * F · e_i * (1 - omega/2) * dt
#[allow(dead_code)]
#[allow(clippy::too_many_arguments)]
pub fn guo_forcing_streaming(
    f: &mut [f64; 9],
    force: [f64; 2],
    ux: f64,
    uy: f64,
    omega: f64,
    dt: f64,
) {
    let weights = [
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
    let ex = [0.0f64, 1.0, 0.0, -1.0, 0.0, 1.0, -1.0, -1.0, 1.0];
    let ey = [0.0f64, 0.0, 1.0, 0.0, -1.0, 1.0, 1.0, -1.0, -1.0];
    let coeff = (1.0 - 0.5 * omega) * dt;
    for i in 0..9 {
        let eu = ex[i] * ux + ey[i] * uy;
        let ef = ex[i] * force[0] + ey[i] * force[1];
        let ue_term = ef / CS2 + (ef * eu - (force[0] * ux + force[1] * uy)) / (CS2 * CS2);
        // Simplified: f_i += w_i * (e_i·F / cs²) * (1 - ω/2) * dt
        f[i] += weights[i] * coeff * ef / CS2;
        let _ = ue_term; // higher-order term stored for reference
    }
}

// ============================================================================
// 10. Acoustic Radiation Pressure
// ============================================================================

/// Compute the Langevin radiation pressure on a perfect reflector.
///
/// P_rad = 2 * `E` = 2 * p_rms² / (ρ0 * c0²)  \[Pa\].
#[allow(dead_code)]
pub fn langevin_radiation_pressure(p_rms: f64, rho0: f64, c0: f64) -> f64 {
    2.0 * p_rms * p_rms / (rho0 * c0 * c0)
}

/// Compute the Gorkov potential U for a small compressible sphere in a
/// standing wave field.
///
/// U = 2π a³ \[ f1 * `p²` / (3 ρ0 c0²) - f2 * ρ0 * `v²` / 2 \]
///
/// where f1 = 1 - κ_p/κ_0,  f2 = 2(ρ_p - ρ0)/(2ρ_p + ρ0),
/// a = sphere radius, κ = compressibility.
#[allow(dead_code)]
pub fn gorkov_potential(
    a: f64,
    p_sq_mean: f64,
    v_sq_mean: f64,
    rho0: f64,
    c0: f64,
    kappa_ratio: f64,
    density_ratio: f64,
) -> f64 {
    let f1 = 1.0 - kappa_ratio;
    let f2 = 2.0 * (density_ratio - 1.0) / (2.0 * density_ratio + 1.0);
    let c0_sq = c0 * c0;
    2.0 * PI * a * a * a * (f1 * p_sq_mean / (3.0 * rho0 * c0_sq) - f2 * rho0 * v_sq_mean / 2.0)
}

/// Compute the acoustic radiation force on a particle as -∇U (gradient of Gorkov potential).
///
/// For a 1-D standing wave p = P0 cos(kx), the force is:
/// F = -dU/dx ∝ sin(2kx).
#[allow(dead_code)]
pub fn gorkov_force_standing_wave_1d(
    a: f64,
    p0: f64,
    x: f64,
    k: f64,
    rho0: f64,
    c0: f64,
    kappa_ratio: f64,
    density_ratio: f64,
) -> f64 {
    let f1 = 1.0 - kappa_ratio;
    let f2 = 2.0 * (density_ratio - 1.0) / (2.0 * density_ratio + 1.0);
    let c0_sq = c0 * c0;
    let p_sq_grad = -p0 * p0 * k * (2.0 * k * x).sin();
    let v_sq_grad = p0 * p0 * k * (2.0 * k * x).sin() / (rho0 * rho0 * c0_sq);
    -2.0 * PI * a * a * a * (f1 * p_sq_grad / (3.0 * rho0 * c0_sq) - f2 * rho0 * v_sq_grad / 2.0)
}

// ============================================================================
// 11. Additional Utilities
// ============================================================================

/// Compute the acoustic absorption coefficient α for a plane wave in air.
///
/// Uses Stokes-Kirchhoff formula:
/// α ≈ ω² / (2 ρ0 c0³) * (4/3 η + κ(1/Cv - 1/Cp))
/// Here a simplified form is used: α ≈ (2 * η * ω²) / (3 * ρ0 * c0³).
#[allow(dead_code)]
pub fn classical_absorption(omega_rad: f64, eta: f64, rho0: f64, c0: f64) -> f64 {
    2.0 * eta * omega_rad * omega_rad / (3.0 * rho0 * c0 * c0 * c0)
}

/// Compute the acoustic impedance of a medium: Z = ρ0 * c0.
#[allow(dead_code)]
pub fn acoustic_impedance(rho0: f64, c0: f64) -> f64 {
    rho0 * c0
}

/// Compute the transmission coefficient for a plane wave at a flat interface.
///
/// T = 2 Z2 / (Z1 + Z2), where Z = ρ c0 is the acoustic impedance.
#[allow(dead_code)]
pub fn transmission_coefficient(z1: f64, z2: f64) -> f64 {
    2.0 * z2 / (z1 + z2)
}

/// Compute the reflection coefficient for a plane wave at a flat interface.
///
/// R = (Z2 - Z1) / (Z2 + Z1).
#[allow(dead_code)]
pub fn reflection_coefficient(z1: f64, z2: f64) -> f64 {
    (z2 - z1) / (z2 + z1)
}

/// Compute the LBM relaxation time tau from the physical kinematic viscosity ν.
///
/// ν = cs² * (tau - 0.5) * dt/dx², so tau = ν/(cs² * dt/dx²) + 0.5.
#[allow(dead_code)]
pub fn tau_from_viscosity(nu: f64, dt: f64, dx: f64) -> f64 {
    nu / (CS2 * dt / (dx * dx)) + 0.5
}

/// Compute the Mach number for a given flow velocity and speed of sound.
#[allow(dead_code)]
pub fn mach_number(u: f64, c0: f64) -> f64 {
    u / c0
}

/// Estimate the acoustic power radiated by a monopole source of volume velocity Q.
///
/// W = ρ0 * c0 * k² * Q² / (4π) for a 3-D monopole.
#[allow(dead_code)]
pub fn monopole_radiated_power(q: f64, freq: f64, rho0: f64, c0: f64) -> f64 {
    let k = 2.0 * PI * freq / c0;
    rho0 * c0 * k * k * q * q / (4.0 * PI)
}

/// Compute the near-field correction factor for a monopole source.
///
/// Factor = 1 + 1/(k*r)².
#[allow(dead_code)]
pub fn near_field_correction(k: f64, r: f64) -> f64 {
    1.0 + 1.0 / (k * r * k * r)
}

/// Frequency from wavenumber and speed of sound: f = k * c0 / (2π).
#[allow(dead_code)]
pub fn freq_from_wavenumber(k: f64, c0: f64) -> f64 {
    k * c0 / (2.0 * PI)
}

/// Wavenumber from frequency: k = 2π f / c0.
#[allow(dead_code)]
pub fn wavenumber_from_freq(freq: f64, c0: f64) -> f64 {
    2.0 * PI * freq / c0
}

/// Compute the acoustic energy density: E = p²/(ρ0 c0²).
#[allow(dead_code)]
pub fn acoustic_energy_density(p: f64, rho0: f64, c0: f64) -> f64 {
    p * p / (rho0 * c0 * c0)
}

/// Compute the acoustic power flux from intensity magnitude.
///
/// W = |I| * A where A is the area element.
#[allow(dead_code)]
pub fn acoustic_power_flux(intensity_mag: f64, area: f64) -> f64 {
    intensity_mag * area
}

/// LBM D2Q9 streaming index bounce-back list (for reference in tests).
#[allow(dead_code)]
pub const D2Q9_OPPOSITE: [usize; 9] = [0, 3, 4, 1, 2, 7, 8, 5, 6];

// ============================================================================
// Unit Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1.0e-10;

    #[test]
    fn test_acoustic_state_pressure() {
        let s = AcousticState::new(3.0, 0.1, 0.2);
        let expected = CS2 * 3.0;
        assert!((s.pressure_perturbation() - expected).abs() < EPS);
    }

    #[test]
    fn test_acoustic_state_intensity() {
        let s = AcousticState::new(3.0, 1.0, 2.0);
        let i = s.intensity();
        let p = CS2 * 3.0;
        assert!((i[0] - p * 1.0).abs() < EPS);
        assert!((i[1] - p * 2.0).abs() < EPS);
    }

    #[test]
    fn test_acoustic_equilibrium_zero_state() {
        let state = AcousticState::new(0.0, 0.0, 0.0);
        let feq = acoustic_equilibrium(&state, 0.0, 0.0);
        for &v in feq.iter() {
            assert!(v.abs() < EPS, "expected zero equilibrium for zero state");
        }
    }

    #[test]
    fn test_acoustic_moments_roundtrip() {
        // Initialize a state, compute equilibrium, extract moments
        let state = AcousticState::new(1.0, 0.05, -0.03);
        let feq = acoustic_equilibrium(&state, 0.0, 0.0);
        let (rho_p, ux_p, uy_p) = acoustic_moments(&feq);
        assert!((rho_p - state.rho_prime).abs() < 1.0e-12);
        assert!((ux_p - state.ux_prime).abs() < 1.0e-12);
        assert!((uy_p - state.uy_prime).abs() < 1.0e-12);
    }

    #[test]
    fn test_bgk_collision_relaxes_to_equilibrium() {
        let state = AcousticState::new(0.5, 0.01, -0.01);
        let feq = acoustic_equilibrium(&state, 0.0, 0.0);
        // Perturb f slightly
        let mut f = feq;
        f[0] += 0.01;
        let omega = 1.0; // tau = 1 => one-step full relaxation
        let f_out = acoustic_bgk_collision(&f, &feq, omega);
        // After one step with omega=1: f_out = feq
        for i in 0..9 {
            assert!((f_out[i] - feq[i]).abs() < EPS);
        }
    }

    #[test]
    fn test_sound_pressure_level_reference() {
        // SPL at reference pressure should be 0 dB
        let spl = sound_pressure_level(P_REF);
        assert!(spl.abs() < 1.0e-9);
    }

    #[test]
    fn test_spl_double_pressure() {
        // Doubling pressure gives +6.02 dB
        let spl1 = sound_pressure_level(1.0e-3);
        let spl2 = sound_pressure_level(2.0e-3);
        assert!((spl2 - spl1 - 20.0 * 2.0_f64.log10()).abs() < 1.0e-9);
    }

    #[test]
    fn test_rms_pressure_empty() {
        assert_eq!(rms_pressure(&[]), 0.0);
    }

    #[test]
    fn test_rms_pressure_constant() {
        let p = vec![2.0; 100];
        let rms = rms_pressure(&p);
        assert!((rms - 2.0).abs() < EPS);
    }

    #[test]
    fn test_a_weighting_1khz() {
        // A-weighting at 1 kHz should be ≈ 0 dB (by definition)
        let aw = a_weighting_db(1000.0);
        assert!(aw.abs() < 1.0, "A-weighting at 1 kHz: {aw}");
    }

    #[test]
    fn test_monopole_source_term_zero_distance() {
        let s = monopole_source_term(1.0, 0.0, 0.0, 0.0, 0.0, C0_AIR, RHO0_AIR);
        assert_eq!(s, 0.0);
    }

    #[test]
    fn test_dipole_source_zero_force() {
        let p = dipole_source_pressure([0.0, 0.0], 0.0, 0.0, 1.0, 1.0, 1000.0, C0_AIR, RHO0_AIR);
        assert_eq!(p, 0.0);
    }

    #[test]
    fn test_quadrupole_source_zero_tensor() {
        let t = [[0.0; 2]; 2];
        let p = quadrupole_source_pressure(t, 0.0, 0.0, 1.0, 1.0, 1000.0, C0_AIR, RHO0_AIR);
        assert_eq!(p, 0.0);
    }

    #[test]
    fn test_inject_monopole_lbm_mass_conservation_relative() {
        // Injecting source_strength should increase total density by source_strength
        let f = [0.0; 9];
        let s = 0.1;
        let f_out = inject_monopole_lbm(&f, s);
        let total: f64 = f_out.iter().sum();
        assert!((total - s).abs() < EPS);
    }

    #[test]
    fn test_pml_sigma_grading() {
        let pml = PmlParams::new(10, 1.0, 2.0);
        assert_eq!(pml.sigma(0.0), 0.0);
        assert_eq!(pml.sigma(10.0), 1.0);
        let mid = pml.sigma(5.0);
        assert!((mid - 0.25).abs() < EPS); // (0.5)^2 = 0.25
    }

    #[test]
    fn test_pml_damping_factor_no_absorption() {
        let pml = PmlParams::new(10, 1.0, 2.0);
        // At d=0, sigma=0 => damping = exp(0) = 1
        assert!((pml.damping_factor(0.0, 0.01) - 1.0).abs() < EPS);
    }

    #[test]
    fn test_optimal_sigma_max() {
        let sigma = optimal_sigma_max(10.0, C0_AIR, 2.0, 0.001);
        assert!(sigma > 0.0);
        // Verify: exp(-2*sigma*10 / (c0*3)) = 0.001
        let r = (-2.0 * sigma * 10.0 / (C0_AIR * 3.0)).exp();
        assert!((r - 0.001).abs() < 1.0e-6);
    }

    #[test]
    fn test_scattering_cross_section_positive() {
        let sigma = scattering_cross_section_cylinder_2d(10.0, 0.1);
        assert!(sigma > 0.0);
    }

    #[test]
    fn test_d2q9_bounce_back_involution() {
        // Applying bounce-back twice should return to original
        for q in 0..9 {
            assert_eq!(d2q9_bounce_back_index(d2q9_bounce_back_index(q)), q);
        }
    }

    #[test]
    fn test_apply_bounce_back_rest_direction() {
        let mut f = [0.0; 9];
        f[0] = 1.0;
        apply_bounce_back(&mut f);
        // Rest direction maps to itself
        assert!((f[0] - 1.0).abs() < EPS);
    }

    #[test]
    fn test_lighthill_tensor_symmetry() {
        let t = lighthill_tensor_2d(1.2, [0.3, 0.1], 101325.0, C0_AIR * C0_AIR);
        // T should be symmetric
        assert!((t[0][1] - t[1][0]).abs() < EPS);
    }

    #[test]
    fn test_room_t60_sabine_positive() {
        let room = RoomGeometry::new(10.0, 8.0, 0.2);
        assert!(room.t60_sabine() > 0.0);
    }

    #[test]
    fn test_room_image_sources_count() {
        let room = RoomGeometry::new(10.0, 8.0, 0.1);
        let imgs = room.image_sources_first_order(3.0, 4.0);
        assert_eq!(imgs.len(), 4);
    }

    #[test]
    fn test_eckart_streaming_force_zero_absorption() {
        let f = eckart_streaming_force([1.0, 0.0], 0.0, C0_AIR);
        assert_eq!(f[0], 0.0);
        assert_eq!(f[1], 0.0);
    }

    #[test]
    fn test_rayleigh_streaming_velocity_positive() {
        let u = rayleigh_streaming_peak_velocity(1.0, 10.0, C0_AIR);
        assert!(u > 0.0);
    }

    #[test]
    fn test_langevin_radiation_pressure_positive() {
        let p_rad = langevin_radiation_pressure(1.0, RHO0_AIR, C0_AIR);
        assert!(p_rad > 0.0);
    }

    #[test]
    fn test_gorkov_potential_zero_amplitude() {
        let u = gorkov_potential(1e-4, 0.0, 0.0, RHO0_AIR, C0_AIR, 0.5, 2.0);
        assert_eq!(u, 0.0);
    }

    #[test]
    fn test_acoustic_impedance() {
        let z = acoustic_impedance(RHO0_AIR, C0_AIR);
        // Expected ~413 Pa·s/m for air
        assert!((z - 413.0).abs() < 2.0);
    }

    #[test]
    fn test_reflection_transmission_sum() {
        let z1 = acoustic_impedance(RHO0_AIR, C0_AIR);
        let z2 = 1.5e6; // water-like
        let r = reflection_coefficient(z1, z2);
        let t = transmission_coefficient(z1, z2);
        // Energy conservation: R + T·Z1/Z2 ≈ 1 ... here just check bounds
        assert!(r < 1.0 && r > -1.0);
        assert!(t > 0.0);
    }

    #[test]
    fn test_tau_from_viscosity() {
        // nu = 1.5e-5 (air), dt=1, dx=1 => tau = 1.5e-5/CS2 + 0.5
        let tau = tau_from_viscosity(1.5e-5, 1.0, 1.0);
        let expected = 1.5e-5 / CS2 + 0.5;
        assert!((tau - expected).abs() < EPS);
    }

    #[test]
    fn test_wavenumber_roundtrip() {
        let freq = 440.0;
        let k = wavenumber_from_freq(freq, C0_AIR);
        let f2 = freq_from_wavenumber(k, C0_AIR);
        assert!((f2 - freq).abs() < 1.0e-10);
    }

    #[test]
    fn test_acoustic_energy_density_positive() {
        let e = acoustic_energy_density(1.0, RHO0_AIR, C0_AIR);
        assert!(e > 0.0);
    }

    #[test]
    fn test_monopole_radiated_power_positive() {
        let w = monopole_radiated_power(1.0e-3, 1000.0, RHO0_AIR, C0_AIR);
        assert!(w > 0.0);
    }

    #[test]
    fn test_acoustic_grid_init_gaussian() {
        let mut grid = AcousticGrid2D::new(32, 32, 1.0, 0.0, 0.0);
        grid.init_gaussian_pulse(16.0, 16.0, 1.0e-3, 3.0);
        let p = grid.pressure_field();
        // Center should have highest pressure
        let center_idx = 16 * 32 + 16;
        let edge_idx = 0;
        assert!(p[center_idx] > p[edge_idx]);
    }

    #[test]
    fn test_acoustic_grid_collide_stream() {
        let mut grid = AcousticGrid2D::new(16, 16, 1.5, 0.0, 0.0);
        grid.init_gaussian_pulse(8.0, 8.0, 1.0e-4, 2.0);
        // Should not panic
        grid.collide();
        grid.stream_periodic();
    }

    #[test]
    fn test_sphere_scattered_pressure_zero_distance() {
        let p = sphere_scattered_pressure(10.0, 0.05, 0.0, 1.0);
        assert_eq!(p, 0.0);
    }

    #[test]
    fn test_near_field_correction_large_r() {
        // For large kr, correction -> 1
        let c = near_field_correction(100.0, 1.0); // k*r = 100
        assert!((c - 1.0).abs() < 1.0e-3);
    }

    #[test]
    fn test_classical_absorption_zero_frequency() {
        let alpha = classical_absorption(0.0, 1.81e-5, RHO0_AIR, C0_AIR);
        assert_eq!(alpha, 0.0);
    }

    #[test]
    fn test_mach_number() {
        let ma = mach_number(C0_AIR / 2.0, C0_AIR);
        assert!((ma - 0.5).abs() < EPS);
    }
}
