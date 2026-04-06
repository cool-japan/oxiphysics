#![allow(clippy::needless_range_loop)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Acoustic Lattice Boltzmann Method (acoustic LBM).
//!
//! Provides linearized LBM for low-amplitude acoustics using a D2Q9 lattice
//! with linearized collision, perfectly matched layer (PML) absorbing
//! boundaries, point/line acoustic sources and receivers, sound speed
//! derivation, acoustic energy density, frequency response, acoustic
//! impedance, far-field approximation, and broadband noise prediction.

#![allow(dead_code)]
#![allow(clippy::too_many_arguments)]

use std::f64::consts::PI;

// ---------------------------------------------------------------------------
// D2Q9 velocity set constants
// ---------------------------------------------------------------------------

/// D2Q9 velocity directions (ex, ey).
const D2Q9_EX: [f64; 9] = [0.0, 1.0, 0.0, -1.0, 0.0, 1.0, -1.0, -1.0, 1.0];
const D2Q9_EY: [f64; 9] = [0.0, 0.0, 1.0, 0.0, -1.0, 1.0, 1.0, -1.0, -1.0];
/// D2Q9 equilibrium weights.
const D2Q9_W: [f64; 9] = [
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
/// D2Q9 opposite direction indices.
const D2Q9_OPP: [usize; 9] = [0, 3, 4, 1, 2, 7, 8, 5, 6];

// ---------------------------------------------------------------------------
// Helper functions
// ---------------------------------------------------------------------------

/// Compute D2Q9 equilibrium distribution for linearized acoustics.
///
/// Uses linearized equilibrium: `f_eq = w_i * (rho + 3*(ex*ux + ey*uy)/cs2)`.
/// For low-amplitude acoustics, `rho = rho0 + rho'`, `|rho'| << rho0`.
#[inline]
pub fn linearized_equilibrium(rho: f64, ux: f64, uy: f64, w: f64, ex: f64, ey: f64) -> f64 {
    let cs2 = 1.0 / 3.0;
    w * (rho + (ex * ux + ey * uy) / cs2)
}

// ---------------------------------------------------------------------------
// 1. AcousticLBM
// ---------------------------------------------------------------------------

/// Linearized LBM solver for low-amplitude acoustics on a D2Q9 lattice.
///
/// This implements the linearized LBM collision operator for acoustic wave
/// propagation in the regime where pressure fluctuations are small compared
/// to the background pressure: p' << p_0.
#[derive(Debug, Clone)]
pub struct AcousticLBM {
    /// Grid width (number of cells in x).
    pub nx: usize,
    /// Grid height (number of cells in y).
    pub ny: usize,
    /// Relaxation parameter omega (related to viscosity).
    pub omega: f64,
    /// Background density (LBM units).
    pub rho0: f64,
    /// Distribution functions f\[y * nx * 9 + x * 9 + q\].
    pub f: Vec<f64>,
    /// Temporary buffer for streaming step.
    pub f_tmp: Vec<f64>,
    /// Grid spacing in physical units (m).
    pub dx: f64,
    /// Time step in physical units (s).
    pub dt: f64,
    /// Current simulation time step count.
    pub step: usize,
}

impl AcousticLBM {
    /// Create a new acoustic LBM grid of size `nx` × `ny`.
    ///
    /// `omega` is the BGK relaxation frequency (0 < omega < 2).
    /// `rho0` is the background density.
    /// `dx` and `dt` are the lattice spacing and time step.
    pub fn new(nx: usize, ny: usize, omega: f64, rho0: f64, dx: f64, dt: f64) -> Self {
        let n = nx * ny * 9;
        let mut f = vec![0.0_f64; n];
        // Initialize to equilibrium at rest
        for y in 0..ny {
            for x in 0..nx {
                for q in 0..9 {
                    f[y * nx * 9 + x * 9 + q] = D2Q9_W[q] * rho0;
                }
            }
        }
        Self {
            nx,
            ny,
            omega,
            rho0,
            f: f.clone(),
            f_tmp: f,
            dx,
            dt,
            step: 0,
        }
    }

    /// Get index into distribution array.
    #[inline]
    fn idx(&self, x: usize, y: usize, q: usize) -> usize {
        y * self.nx * 9 + x * 9 + q
    }

    /// Compute density (pressure fluctuation) at cell (x, y).
    pub fn density(&self, x: usize, y: usize) -> f64 {
        (0..9).map(|q| self.f[self.idx(x, y, q)]).sum()
    }

    /// Compute x-velocity at cell (x, y).
    pub fn velocity_x(&self, x: usize, y: usize) -> f64 {
        let rho = self.density(x, y);
        if rho < 1e-30 {
            return 0.0;
        }
        let mom: f64 = (0..9).map(|q| self.f[self.idx(x, y, q)] * D2Q9_EX[q]).sum();
        mom / rho
    }

    /// Compute y-velocity at cell (x, y).
    pub fn velocity_y(&self, x: usize, y: usize) -> f64 {
        let rho = self.density(x, y);
        if rho < 1e-30 {
            return 0.0;
        }
        let mom: f64 = (0..9).map(|q| self.f[self.idx(x, y, q)] * D2Q9_EY[q]).sum();
        mom / rho
    }

    /// Perform one collision step using linearized BGK operator.
    ///
    /// The collision is: `f_q = f_q - omega*(f_q - f_eq_q)` where `f_eq`
    /// is the linearized equilibrium distribution.
    pub fn collide(&mut self) {
        for y in 0..self.ny {
            for x in 0..self.nx {
                let rho = self.density(x, y);
                let ux = self.velocity_x(x, y);
                let uy = self.velocity_y(x, y);
                for q in 0..9 {
                    let feq =
                        linearized_equilibrium(rho, ux, uy, D2Q9_W[q], D2Q9_EX[q], D2Q9_EY[q]);
                    let idx = self.idx(x, y, q);
                    self.f[idx] -= self.omega * (self.f[idx] - feq);
                }
            }
        }
    }

    /// Perform one streaming step with periodic boundaries.
    pub fn stream(&mut self) {
        let nx = self.nx;
        let ny = self.ny;
        self.f_tmp.copy_from_slice(&self.f);
        for y in 0..ny {
            for x in 0..nx {
                for q in 0..9 {
                    let ex = D2Q9_EX[q] as isize;
                    let ey = D2Q9_EY[q] as isize;
                    let xd = ((x as isize + ex).rem_euclid(nx as isize)) as usize;
                    let yd = ((y as isize + ey).rem_euclid(ny as isize)) as usize;
                    let dst = self.idx(xd, yd, q);
                    let src = self.idx(x, y, q);
                    self.f_tmp[dst] = self.f[src];
                }
            }
        }
        std::mem::swap(&mut self.f, &mut self.f_tmp);
    }

    /// Perform one full LBM step: collide then stream.
    pub fn step(&mut self) {
        self.collide();
        self.stream();
        self.step += 1;
    }

    /// Add a pressure perturbation `dp` at cell (x, y) by distributing
    /// it according to D2Q9 weights (monopole source injection).
    pub fn add_pressure_source(&mut self, x: usize, y: usize, dp: f64) {
        for q in 0..9 {
            let idx = self.idx(x, y, q);
            self.f[idx] += D2Q9_W[q] * dp;
        }
    }

    /// Get pressure fluctuation at (x, y): `p' = cs^2 * (rho - rho0)`.
    pub fn pressure_fluctuation(&self, x: usize, y: usize) -> f64 {
        let cs2 = 1.0 / 3.0;
        cs2 * (self.density(x, y) - self.rho0)
    }

    /// Return a 2D pressure field as a flat Vec of length nx*ny.
    pub fn pressure_field(&self) -> Vec<f64> {
        (0..self.ny)
            .flat_map(|y| (0..self.nx).map(move |x| (x, y)))
            .map(|(x, y)| self.pressure_fluctuation(x, y))
            .collect()
    }
}

// ---------------------------------------------------------------------------
// 2. PmlLayer
// ---------------------------------------------------------------------------

/// Perfectly Matched Layer (PML) absorbing boundary condition.
///
/// Applies exponential damping in the PML region to absorb outgoing waves
/// without spurious reflections. The damping profile is quadratic.
#[derive(Debug, Clone)]
pub struct PmlLayer {
    /// Thickness of PML in number of cells.
    pub thickness: usize,
    /// Maximum damping coefficient at the outer boundary (s^-1).
    pub sigma_max: f64,
    /// Damping in x-direction per cell, indexed by cell x-coordinate.
    pub d_x: Vec<f64>,
    /// Damping in y-direction per cell, indexed by cell y-coordinate.
    pub d_y: Vec<f64>,
}

impl PmlLayer {
    /// Create a PML layer for a grid of size `nx` × `ny`.
    ///
    /// The PML of given `thickness` cells is applied to all four boundaries.
    /// `sigma_max` controls absorption strength.
    pub fn new(nx: usize, ny: usize, thickness: usize, sigma_max: f64) -> Self {
        let mut d_x = vec![0.0_f64; nx];
        let mut d_y = vec![0.0_f64; ny];
        let d = thickness as f64;
        // Left boundary
        for i in 0..thickness.min(nx) {
            let xi = (thickness - i) as f64 / d;
            d_x[i] = sigma_max * xi * xi;
        }
        // Right boundary
        for i in (nx.saturating_sub(thickness))..nx {
            let xi = (i + thickness + 1 - nx) as f64 / d;
            d_x[i] = sigma_max * xi * xi;
        }
        // Bottom boundary
        for j in 0..thickness.min(ny) {
            let xi = (thickness - j) as f64 / d;
            d_y[j] = sigma_max * xi * xi;
        }
        // Top boundary
        for j in (ny.saturating_sub(thickness))..ny {
            let xi = (j + thickness + 1 - ny) as f64 / d;
            d_y[j] = sigma_max * xi * xi;
        }
        Self {
            thickness,
            sigma_max,
            d_x,
            d_y,
        }
    }

    /// Apply PML damping to distribution functions `f` for a cell at (x, y).
    ///
    /// The damping factor is: `f *= exp(-sigma * dt)`.
    pub fn apply(&self, f: &mut [f64], x: usize, y: usize, dt: f64) {
        let sigma = (self.d_x[x] + self.d_y[y]) * 0.5;
        if sigma > 0.0 {
            let decay = (-sigma * dt).exp();
            for q in f.iter_mut() {
                *q *= decay;
            }
        }
    }

    /// Return the total damping coefficient at position (x, y).
    pub fn sigma_at(&self, x: usize, y: usize) -> f64 {
        (self.d_x[x] + self.d_y[y]) * 0.5
    }

    /// Check if cell (x, y) is inside the PML region.
    pub fn is_pml(&self, x: usize, y: usize, nx: usize, ny: usize) -> bool {
        x < self.thickness
            || x >= nx.saturating_sub(self.thickness)
            || y < self.thickness
            || y >= ny.saturating_sub(self.thickness)
    }
}

// ---------------------------------------------------------------------------
// 3. AcousticSource
// ---------------------------------------------------------------------------

/// Type of acoustic source.
#[derive(Debug, Clone, PartialEq)]
pub enum SourceType {
    /// Point source (monopole) at a single grid cell.
    Point,
    /// Line source along x-direction at fixed y.
    LineX,
    /// Line source along y-direction at fixed x.
    LineY,
}

/// Point or line acoustic source emitting a sinusoidal pressure pulse.
///
/// Injects a harmonic pressure perturbation at the specified location.
#[derive(Debug, Clone)]
pub struct AcousticSource {
    /// Grid x-position of source.
    pub x: usize,
    /// Grid y-position of source.
    pub y: usize,
    /// Source type (point, line-x, line-y).
    pub source_type: SourceType,
    /// Angular frequency of the source (rad/s).
    pub omega: f64,
    /// Amplitude of pressure fluctuation (LBM units).
    pub amplitude: f64,
    /// Phase offset (radians).
    pub phase: f64,
}

impl AcousticSource {
    /// Create a sinusoidal point source at (x, y).
    ///
    /// `freq` is the physical frequency in Hz, converted using `dt`.
    pub fn new_point(x: usize, y: usize, freq: f64, amplitude: f64, dt: f64) -> Self {
        Self {
            x,
            y,
            source_type: SourceType::Point,
            omega: 2.0 * PI * freq * dt,
            amplitude,
            phase: 0.0,
        }
    }

    /// Create a line source along x at y-position `y`.
    pub fn new_line_x(y: usize, freq: f64, amplitude: f64, dt: f64) -> Self {
        Self {
            x: 0,
            y,
            source_type: SourceType::LineX,
            omega: 2.0 * PI * freq * dt,
            amplitude,
            phase: 0.0,
        }
    }

    /// Compute the pressure injection value at time step `t`.
    pub fn pressure_at(&self, t: usize) -> f64 {
        self.amplitude * (self.omega * t as f64 + self.phase).sin()
    }

    /// Apply source injection to a LBM solver at step `t`.
    pub fn inject(&self, lbm: &mut AcousticLBM, t: usize) {
        let dp = self.pressure_at(t);
        match self.source_type {
            SourceType::Point => lbm.add_pressure_source(self.x, self.y, dp),
            SourceType::LineX => {
                for x in 0..lbm.nx {
                    lbm.add_pressure_source(x, self.y, dp);
                }
            }
            SourceType::LineY => {
                for y in 0..lbm.ny {
                    lbm.add_pressure_source(self.x, y, dp);
                }
            }
        }
    }
}

// ---------------------------------------------------------------------------
// 4. AcousticReceiver
// ---------------------------------------------------------------------------

/// Records the pressure time history at a single grid point.
///
/// Use multiple receivers to sample the acoustic field at different locations.
#[derive(Debug, Clone)]
pub struct AcousticReceiver {
    /// Grid x-position of receiver.
    pub x: usize,
    /// Grid y-position of receiver.
    pub y: usize,
    /// Recorded pressure fluctuation time series.
    pub history: Vec<f64>,
}

impl AcousticReceiver {
    /// Create a new receiver at (x, y).
    pub fn new(x: usize, y: usize) -> Self {
        Self {
            x,
            y,
            history: Vec::new(),
        }
    }

    /// Sample the current pressure from the LBM solver and record it.
    pub fn record(&mut self, lbm: &AcousticLBM) {
        self.history.push(lbm.pressure_fluctuation(self.x, self.y));
    }

    /// Clear the recorded history.
    pub fn reset(&mut self) {
        self.history.clear();
    }

    /// Return the peak pressure recorded.
    pub fn peak_pressure(&self) -> f64 {
        self.history.iter().cloned().fold(0.0_f64, f64::max)
    }

    /// Return the RMS pressure.
    pub fn rms_pressure(&self) -> f64 {
        if self.history.is_empty() {
            return 0.0;
        }
        let sum_sq: f64 = self.history.iter().map(|p| p * p).sum();
        (sum_sq / self.history.len() as f64).sqrt()
    }
}

// ---------------------------------------------------------------------------
// 5. SoundSpeed
// ---------------------------------------------------------------------------

/// Sound speed derived from LBM lattice parameters.
///
/// In D2Q9/D3Q19 LBM, the lattice sound speed is: `c_s = dx / (sqrt(3) * dt)`.
#[derive(Debug, Clone, Copy)]
pub struct SoundSpeed {
    /// Lattice spacing (m).
    pub dx: f64,
    /// Time step (s).
    pub dt: f64,
}

impl SoundSpeed {
    /// Create a new SoundSpeed from lattice parameters.
    pub fn new(dx: f64, dt: f64) -> Self {
        Self { dx, dt }
    }

    /// Return the LBM lattice sound speed in m/s.
    ///
    /// Formula: `c_s = dx / (sqrt(3) * dt)`.
    pub fn lattice_sound_speed(&self) -> f64 {
        self.dx / (3.0_f64.sqrt() * self.dt)
    }

    /// Return the lattice Mach number for a given physical velocity.
    pub fn mach_number(&self, u_phys: f64) -> f64 {
        u_phys / self.lattice_sound_speed()
    }

    /// Return the LBM speed of sound squared: `cs2 = 1/3` in lattice units.
    pub fn cs2_lattice() -> f64 {
        1.0 / 3.0
    }

    /// Convert physical frequency to LBM frequency (per timestep).
    pub fn physical_to_lbm_freq(&self, freq_hz: f64) -> f64 {
        freq_hz * self.dt
    }

    /// Return wavelength in LBM cells for a given physical frequency.
    pub fn wavelength_cells(&self, freq_hz: f64) -> f64 {
        self.lattice_sound_speed() / freq_hz / self.dx
    }
}

// ---------------------------------------------------------------------------
// 6. AcousticEnergy
// ---------------------------------------------------------------------------

/// Acoustic energy density calculator.
///
/// Computes: `E = p'^2 / (2*rho0*c_s^2) + rho0*(ux^2 + uy^2)/2`.
#[derive(Debug, Clone, Copy)]
pub struct AcousticEnergy {
    /// Background density (kg/m³ in physical, or LBM units).
    pub rho0: f64,
    /// Sound speed (m/s or LBM units).
    pub c_s: f64,
}

impl AcousticEnergy {
    /// Create with given background density and sound speed.
    pub fn new(rho0: f64, c_s: f64) -> Self {
        Self { rho0, c_s }
    }

    /// Compute acoustic potential energy density.
    ///
    /// `E_pot = p'^2 / (2 * rho0 * c_s^2)`.
    pub fn potential_energy(&self, p_prime: f64) -> f64 {
        p_prime * p_prime / (2.0 * self.rho0 * self.c_s * self.c_s)
    }

    /// Compute acoustic kinetic energy density.
    ///
    /// `E_kin = rho0 * (ux^2 + uy^2) / 2`.
    pub fn kinetic_energy(&self, ux: f64, uy: f64) -> f64 {
        0.5 * self.rho0 * (ux * ux + uy * uy)
    }

    /// Compute total acoustic energy density at a point.
    pub fn total_energy(&self, p_prime: f64, ux: f64, uy: f64) -> f64 {
        self.potential_energy(p_prime) + self.kinetic_energy(ux, uy)
    }

    /// Integrate total acoustic energy over the entire LBM grid.
    pub fn integrate_field(&self, lbm: &AcousticLBM) -> f64 {
        let mut total = 0.0;
        for y in 0..lbm.ny {
            for x in 0..lbm.nx {
                let p = lbm.pressure_fluctuation(x, y);
                let ux = lbm.velocity_x(x, y);
                let uy = lbm.velocity_y(x, y);
                total += self.total_energy(p, ux, uy);
            }
        }
        total * lbm.dx * lbm.dx
    }
}

// ---------------------------------------------------------------------------
// 7. FrequencyResponse
// ---------------------------------------------------------------------------

/// Frequency response computed via DFT of recorded pressure time history.
///
/// Provides magnitude and phase spectrum of acoustic signals recorded
/// at receiver positions.
#[derive(Debug, Clone)]
pub struct FrequencyResponse {
    /// Sample rate (1/dt in Hz or per step).
    pub sample_rate: f64,
    /// Number of frequency bins.
    pub n_fft: usize,
    /// Magnitude spectrum.
    pub magnitude: Vec<f64>,
    /// Phase spectrum (radians).
    pub phase: Vec<f64>,
    /// Frequency axis values.
    pub frequencies: Vec<f64>,
}

impl FrequencyResponse {
    /// Compute frequency response from a pressure time series using DFT.
    ///
    /// `signal` is the pressure time history.
    /// `dt` is the time step in seconds.
    pub fn from_signal(signal: &[f64], dt: f64) -> Self {
        let n = signal.len();
        let n_fft = n / 2 + 1;
        let sample_rate = 1.0 / dt;
        let mut magnitude = vec![0.0_f64; n_fft];
        let mut phase = vec![0.0_f64; n_fft];
        let mut frequencies = vec![0.0_f64; n_fft];

        for k in 0..n_fft {
            let mut re = 0.0_f64;
            let mut im = 0.0_f64;
            for (t, &s) in signal.iter().enumerate() {
                let angle = -2.0 * PI * k as f64 * t as f64 / n as f64;
                re += s * angle.cos();
                im += s * angle.sin();
            }
            magnitude[k] = (re * re + im * im).sqrt() / n as f64;
            phase[k] = im.atan2(re);
            frequencies[k] = k as f64 * sample_rate / n as f64;
        }

        Self {
            sample_rate,
            n_fft,
            magnitude,
            phase,
            frequencies,
        }
    }

    /// Return the dominant frequency (peak in magnitude spectrum).
    pub fn dominant_frequency(&self) -> f64 {
        let peak_idx = self
            .magnitude
            .iter()
            .enumerate()
            .max_by(|a, b| a.1.partial_cmp(b.1).unwrap_or(std::cmp::Ordering::Equal))
            .map(|(i, _)| i)
            .unwrap_or(0);
        self.frequencies[peak_idx]
    }

    /// Return the -3 dB bandwidth around the dominant frequency.
    pub fn bandwidth_3db(&self) -> f64 {
        let peak = self.magnitude.iter().cloned().fold(0.0_f64, f64::max);
        let threshold = peak / 2.0_f64.sqrt();
        let mut f_low = self.frequencies[0];
        let mut f_high = self.frequencies[self.n_fft - 1];
        for i in 1..self.n_fft {
            if self.magnitude[i] >= threshold {
                f_low = self.frequencies[i - 1];
                break;
            }
        }
        for i in (0..self.n_fft - 1).rev() {
            if self.magnitude[i] >= threshold {
                f_high = self.frequencies[i + 1].min(self.frequencies[self.n_fft - 1]);
                break;
            }
        }
        f_high - f_low
    }

    /// Compute sound pressure level (SPL) in dB at a given frequency index.
    ///
    /// `p_ref` is the reference pressure (e.g. 20 μPa in air).
    pub fn spl_db(&self, idx: usize, p_ref: f64) -> f64 {
        if p_ref < 1e-30 || self.magnitude[idx] < 1e-30 {
            return -f64::INFINITY;
        }
        20.0 * (self.magnitude[idx] / p_ref).log10()
    }
}

// ---------------------------------------------------------------------------
// 8. AcousticImpedance
// ---------------------------------------------------------------------------

/// Characteristic acoustic impedance and reflection/transmission coefficients.
///
/// The characteristic impedance of a medium is `Z = rho * c`.
/// The pressure reflection coefficient at an interface is:
/// `R = (Z2 - Z1) / (Z2 + Z1)`.
#[derive(Debug, Clone, Copy)]
pub struct AcousticImpedance {
    /// Density of medium 1 (kg/m³).
    pub rho1: f64,
    /// Sound speed of medium 1 (m/s).
    pub c1: f64,
    /// Density of medium 2 (kg/m³).
    pub rho2: f64,
    /// Sound speed of medium 2 (m/s).
    pub c2: f64,
}

impl AcousticImpedance {
    /// Create with properties of two media.
    pub fn new(rho1: f64, c1: f64, rho2: f64, c2: f64) -> Self {
        Self { rho1, c1, rho2, c2 }
    }

    /// Characteristic impedance of medium 1: `Z1 = rho1 * c1`.
    pub fn z1(&self) -> f64 {
        self.rho1 * self.c1
    }

    /// Characteristic impedance of medium 2: `Z2 = rho2 * c2`.
    pub fn z2(&self) -> f64 {
        self.rho2 * self.c2
    }

    /// Pressure reflection coefficient: `R = (Z2 - Z1) / (Z2 + Z1)`.
    pub fn reflection_coefficient(&self) -> f64 {
        let z1 = self.z1();
        let z2 = self.z2();
        (z2 - z1) / (z2 + z1)
    }

    /// Pressure transmission coefficient: `T = 2*Z2 / (Z2 + Z1)`.
    pub fn transmission_coefficient(&self) -> f64 {
        let z1 = self.z1();
        let z2 = self.z2();
        2.0 * z2 / (z2 + z1)
    }

    /// Intensity reflection coefficient: `R_I = R^2`.
    pub fn intensity_reflection(&self) -> f64 {
        self.reflection_coefficient().powi(2)
    }

    /// Intensity transmission coefficient: `T_I = 1 - R_I`.
    pub fn intensity_transmission(&self) -> f64 {
        1.0 - self.intensity_reflection()
    }
}

// ---------------------------------------------------------------------------
// 9. FarFieldApproximation
// ---------------------------------------------------------------------------

/// Far-field pressure approximation using the Ffowcs Williams-Hawkings analogy.
///
/// Computes the acoustic pressure at a far-field observer location from
/// compact-source surface pressure and velocity data.
#[derive(Debug, Clone)]
pub struct FarFieldApproximation {
    /// Observer position (m).
    pub observer: [f64; 2],
    /// Reference sound speed (m/s).
    pub c0: f64,
    /// Reference density (kg/m³).
    pub rho0: f64,
}

impl FarFieldApproximation {
    /// Create a far-field approximation for an observer at position `obs`.
    pub fn new(observer: [f64; 2], c0: f64, rho0: f64) -> Self {
        Self { observer, c0, rho0 }
    }

    /// Compute the direction cosine vector from source to observer.
    ///
    /// Returns normalized vector `r_hat` from `source` to `observer`.
    pub fn direction_to_observer(&self, source: [f64; 2]) -> [f64; 2] {
        let dx = self.observer[0] - source[0];
        let dy = self.observer[1] - source[1];
        let r = (dx * dx + dy * dy).sqrt();
        if r < 1e-30 {
            [1.0, 0.0]
        } else {
            [dx / r, dy / r]
        }
    }

    /// Distance from source to observer (m).
    pub fn distance(&self, source: [f64; 2]) -> f64 {
        let dx = self.observer[0] - source[0];
        let dy = self.observer[1] - source[1];
        (dx * dx + dy * dy).sqrt()
    }

    /// Far-field pressure from a monopole source with volume velocity `Q`.
    ///
    /// `Q` is the acoustic volume velocity (m²/s in 2D).
    /// `r` is the distance to observer (m).
    /// Returns `p' ≈ rho0 * dQ/dt / (2*pi*r)` in 2D cylindrical.
    pub fn monopole_pressure(&self, dqdt: f64, r: f64) -> f64 {
        self.rho0 * dqdt / (2.0 * PI * r)
    }

    /// Far-field pressure from a dipole source with force `F` (N/m).
    ///
    /// `theta` is the angle between force direction and observer direction.
    pub fn dipole_pressure(&self, f_mag: f64, r: f64, theta: f64) -> f64 {
        -f_mag * theta.cos() / (2.0 * PI * r * self.c0)
    }

    /// FWH surface integral contribution from pressure on a surface element.
    ///
    /// `p_s` is surface pressure, `n_dot_r` is surface normal dotted with
    /// direction to observer, `ds` is the surface element length.
    pub fn fwh_surface_contribution(&self, p_s: f64, n_dot_r: f64, r: f64, ds: f64) -> f64 {
        p_s * n_dot_r * ds / (2.0 * PI * r)
    }

    /// Accumulate far-field pressure from multiple surface points.
    pub fn integrate_surface(&self, pressures: &[f64], sources: &[[f64; 2]], ds: f64) -> f64 {
        pressures
            .iter()
            .zip(sources.iter())
            .map(|(&p_s, &src)| {
                let r = self.distance(src);
                let rhat = self.direction_to_observer(src);
                // Assume outward normal points from source to observer
                let n_dot_r = rhat[0] * rhat[0] + rhat[1] * rhat[1]; // = 1 for outward
                self.fwh_surface_contribution(p_s, n_dot_r, r.max(1e-10), ds)
            })
            .sum()
    }
}

// ---------------------------------------------------------------------------
// 10. NoisePrediction
// ---------------------------------------------------------------------------

/// Broadband noise prediction using Proudman's formula for volume sources.
///
/// Proudman's formula gives the acoustic power per unit volume from turbulence:
/// `P = alpha_e * rho * u^5 / (c_s * L)` where `alpha_e ≈ 0.1` (Proudman constant).
#[derive(Debug, Clone, Copy)]
pub struct NoisePrediction {
    /// Reference density (kg/m³).
    pub rho0: f64,
    /// Reference sound speed (m/s).
    pub c0: f64,
    /// Proudman's constant (typically 0.1).
    pub alpha_e: f64,
    /// Turbulence length scale (m).
    pub length_scale: f64,
}

impl NoisePrediction {
    /// Create a NoisePrediction model.
    ///
    /// Uses Proudman's constant `alpha_e = 0.1` by default.
    pub fn new(rho0: f64, c0: f64, length_scale: f64) -> Self {
        Self {
            rho0,
            c0,
            alpha_e: 0.1,
            length_scale,
        }
    }

    /// Acoustic power per unit volume via Proudman's formula.
    ///
    /// `u_rms` is the RMS turbulent velocity fluctuation.
    /// Returns power in W/m³ (or W/m² in 2D).
    pub fn proudman_power(&self, u_rms: f64) -> f64 {
        self.alpha_e * self.rho0 * u_rms.powi(5) / (self.c0 * self.length_scale)
    }

    /// Acoustic efficiency (ratio of acoustic to turbulent kinetic energy flux).
    ///
    /// `M_t` is the turbulent Mach number `u_rms / c0`.
    pub fn acoustic_efficiency(&self, u_rms: f64) -> f64 {
        let m_t = u_rms / self.c0;
        self.alpha_e * m_t.powi(5) / (u_rms.powi(3) / self.length_scale)
    }

    /// Estimate the peak acoustic frequency from turbulence.
    ///
    /// Using Strouhal number St ≈ 0.2: `f_peak = St * u_rms / L`.
    pub fn peak_frequency(&self, u_rms: f64) -> f64 {
        0.2 * u_rms / self.length_scale
    }

    /// Compute the overall sound pressure level (OASPL) in dB.
    ///
    /// `p_ref` is the reference pressure (20 μPa in air).
    pub fn oaspl(&self, u_rms: f64, p_ref: f64) -> f64 {
        let power = self.proudman_power(u_rms);
        // Convert power to RMS pressure via p = sqrt(rho0 * c0 * power * L)
        let p_rms = (self.rho0 * self.c0 * power * self.length_scale).sqrt();
        20.0 * (p_rms / p_ref).log10()
    }

    /// Compute the A-weighted SPL correction at a given frequency (Hz).
    ///
    /// Uses simplified A-weighting approximation.
    pub fn a_weighting_db(freq_hz: f64) -> f64 {
        // Simplified A-weighting formula
        let f2 = freq_hz * freq_hz;
        let f4 = f2 * f2;
        let num = 12194.0_f64.powi(2) * f4;
        let d1 = f2 + 20.6_f64.powi(2);
        let d2 = (f2 + 107.7_f64.powi(2)).sqrt() * (f2 + 737.9_f64.powi(2)).sqrt();
        let d3 = f2 + 12194.0_f64.powi(2);
        20.0 * (num / (d1 * d2 * d3)).log10() + 2.0
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // --- AcousticLBM tests ---

    #[test]
    fn test_acoustic_lbm_init_density() {
        let lbm = AcousticLBM::new(8, 8, 1.0, 1.0, 1.0, 1.0);
        let rho = lbm.density(4, 4);
        assert!((rho - 1.0).abs() < 1e-12, "Initial density should be rho0");
    }

    #[test]
    fn test_acoustic_lbm_zero_velocity_at_rest() {
        let lbm = AcousticLBM::new(8, 8, 1.0, 1.0, 1.0, 1.0);
        assert!(lbm.velocity_x(4, 4).abs() < 1e-12);
        assert!(lbm.velocity_y(4, 4).abs() < 1e-12);
    }

    #[test]
    fn test_acoustic_lbm_pressure_fluctuation_zero() {
        let lbm = AcousticLBM::new(8, 8, 1.0, 1.0, 1.0, 1.0);
        let p = lbm.pressure_fluctuation(3, 3);
        assert!(
            p.abs() < 1e-12,
            "Initial pressure fluctuation should be zero"
        );
    }

    #[test]
    fn test_acoustic_lbm_step_conserves_mass() {
        let mut lbm = AcousticLBM::new(10, 10, 1.0, 1.0, 1.0, 1.0);
        let mass_before: f64 = (0..10)
            .flat_map(|y| (0..10).map(move |x| (x, y)))
            .map(|(x, y)| lbm.density(x, y))
            .sum();
        lbm.step();
        let mass_after: f64 = (0..10)
            .flat_map(|y| (0..10).map(move |x| (x, y)))
            .map(|(x, y)| lbm.density(x, y))
            .sum();
        assert!((mass_before - mass_after).abs() < 1e-10);
    }

    #[test]
    fn test_acoustic_lbm_add_pressure_source() {
        let mut lbm = AcousticLBM::new(8, 8, 1.0, 1.0, 1.0, 1.0);
        lbm.add_pressure_source(4, 4, 0.01);
        let rho = lbm.density(4, 4);
        assert!((rho - 1.01).abs() < 1e-12);
    }

    #[test]
    fn test_acoustic_lbm_step_increments_counter() {
        let mut lbm = AcousticLBM::new(8, 8, 1.0, 1.0, 1.0, 1.0);
        lbm.step();
        lbm.step();
        assert_eq!(lbm.step, 2);
    }

    #[test]
    fn test_acoustic_lbm_pressure_field_length() {
        let lbm = AcousticLBM::new(6, 7, 1.0, 1.0, 1.0, 1.0);
        assert_eq!(lbm.pressure_field().len(), 6 * 7);
    }

    #[test]
    fn test_acoustic_lbm_linearized_eq_sum_to_rho() {
        let rho = 1.05;
        let ux = 0.01;
        let uy = -0.005;
        let sum: f64 = (0..9)
            .map(|q| linearized_equilibrium(rho, ux, uy, D2Q9_W[q], D2Q9_EX[q], D2Q9_EY[q]))
            .sum();
        assert!(
            (sum - rho).abs() < 1e-12,
            "Linearized equilibrium should sum to rho"
        );
    }

    // --- PmlLayer tests ---

    #[test]
    fn test_pml_zero_interior() {
        let pml = PmlLayer::new(20, 20, 4, 10.0);
        // Interior cells should have zero damping
        assert_eq!(pml.sigma_at(10, 10), 0.0);
    }

    #[test]
    fn test_pml_nonzero_boundary() {
        let pml = PmlLayer::new(20, 20, 4, 10.0);
        // Corner cell should have high damping
        assert!(pml.sigma_at(0, 0) > 0.0);
    }

    #[test]
    fn test_pml_is_pml_corner() {
        let pml = PmlLayer::new(20, 20, 4, 10.0);
        assert!(pml.is_pml(0, 0, 20, 20));
        assert!(pml.is_pml(19, 19, 20, 20));
    }

    #[test]
    fn test_pml_is_not_pml_center() {
        let pml = PmlLayer::new(20, 20, 4, 10.0);
        assert!(!pml.is_pml(10, 10, 20, 20));
    }

    #[test]
    fn test_pml_apply_reduces_amplitude() {
        let pml = PmlLayer::new(20, 20, 4, 10.0);
        let mut f = vec![1.0_f64; 9];
        pml.apply(&mut f, 0, 0, 0.1);
        assert!(f[0] < 1.0, "PML should damp the distribution");
    }

    #[test]
    fn test_pml_no_damping_interior() {
        let pml = PmlLayer::new(20, 20, 4, 10.0);
        let mut f = vec![1.0_f64; 9];
        pml.apply(&mut f, 10, 10, 0.1);
        assert!((f[0] - 1.0).abs() < 1e-12, "No damping in interior");
    }

    // --- AcousticSource tests ---

    #[test]
    fn test_acoustic_source_pressure_sinusoidal() {
        let src = AcousticSource::new_point(5, 5, 100.0, 0.01, 1e-4);
        let p = src.pressure_at(0);
        assert_eq!(p, 0.0, "Sine wave starts at 0");
    }

    #[test]
    fn test_acoustic_source_pressure_nonzero_t() {
        let src = AcousticSource::new_point(5, 5, 100.0, 0.01, 1e-4);
        let p = src.pressure_at(25); // quarter period
        assert!(p.abs() > 0.0);
    }

    #[test]
    fn test_acoustic_source_inject_point() {
        let mut lbm = AcousticLBM::new(10, 10, 1.0, 1.0, 1.0, 1.0);
        let src = AcousticSource::new_point(5, 5, 100.0, 0.01, 1.0);
        // At t=25, sin(omega*25) != 0
        src.inject(&mut lbm, 25);
        // Pressure at source should change (unless sin=0 exactly)
        // Just check no panic
        let _ = lbm.density(5, 5);
    }

    #[test]
    fn test_acoustic_source_line_x_injects_row() {
        let mut lbm = AcousticLBM::new(10, 10, 1.0, 1.0, 1.0, 1.0);
        let src = AcousticSource::new_line_x(5, 100.0, 0.1, 1.0);
        src.inject(&mut lbm, 25);
        // All cells in row y=5 should have same density change
        let rho0 = lbm.density(0, 5);
        let rho1 = lbm.density(9, 5);
        assert!((rho0 - rho1).abs() < 1e-12);
    }

    // --- AcousticReceiver tests ---

    #[test]
    fn test_acoustic_receiver_record() {
        let lbm = AcousticLBM::new(8, 8, 1.0, 1.0, 1.0, 1.0);
        let mut recv = AcousticReceiver::new(4, 4);
        recv.record(&lbm);
        assert_eq!(recv.history.len(), 1);
    }

    #[test]
    fn test_acoustic_receiver_reset() {
        let lbm = AcousticLBM::new(8, 8, 1.0, 1.0, 1.0, 1.0);
        let mut recv = AcousticReceiver::new(4, 4);
        recv.record(&lbm);
        recv.reset();
        assert!(recv.history.is_empty());
    }

    #[test]
    fn test_acoustic_receiver_rms_zero_at_rest() {
        let lbm = AcousticLBM::new(8, 8, 1.0, 1.0, 1.0, 1.0);
        let mut recv = AcousticReceiver::new(4, 4);
        for _ in 0..10 {
            recv.record(&lbm);
        }
        assert!(recv.rms_pressure().abs() < 1e-12);
    }

    #[test]
    fn test_acoustic_receiver_empty_rms() {
        let recv = AcousticReceiver::new(0, 0);
        assert_eq!(recv.rms_pressure(), 0.0);
    }

    // --- SoundSpeed tests ---

    #[test]
    fn test_sound_speed_lattice() {
        let ss = SoundSpeed::new(1.0, 1.0);
        let c = ss.lattice_sound_speed();
        assert!((c - 1.0 / 3.0_f64.sqrt()).abs() < 1e-12);
    }

    #[test]
    fn test_sound_speed_mach_number() {
        let ss = SoundSpeed::new(1.0, 1.0);
        let c = ss.lattice_sound_speed();
        let mach = ss.mach_number(c);
        assert!((mach - 1.0).abs() < 1e-12, "Mach=1 at sound speed");
    }

    #[test]
    fn test_sound_speed_cs2_lattice() {
        assert!((SoundSpeed::cs2_lattice() - 1.0 / 3.0).abs() < 1e-14);
    }

    #[test]
    fn test_sound_speed_wavelength() {
        let ss = SoundSpeed::new(1.0, 1.0);
        // wavelength = c / f / dx; c = 1/sqrt(3), f=1/sqrt(3), dx=1 → wavelength = 1 cell
        let wl = ss.wavelength_cells(ss.lattice_sound_speed());
        assert!((wl - 1.0).abs() < 1e-12);
    }

    // --- AcousticEnergy tests ---

    #[test]
    fn test_acoustic_energy_potential_zero_at_rest() {
        let ae = AcousticEnergy::new(1.0, 1.0 / 3.0_f64.sqrt());
        assert_eq!(ae.potential_energy(0.0), 0.0);
    }

    #[test]
    fn test_acoustic_energy_kinetic_zero_at_rest() {
        let ae = AcousticEnergy::new(1.0, 1.0 / 3.0_f64.sqrt());
        assert_eq!(ae.kinetic_energy(0.0, 0.0), 0.0);
    }

    #[test]
    fn test_acoustic_energy_total_positive() {
        let ae = AcousticEnergy::new(1.0, 1.0 / 3.0_f64.sqrt());
        let e = ae.total_energy(0.01, 0.001, 0.0);
        assert!(e > 0.0);
    }

    #[test]
    fn test_acoustic_energy_integrate_zero_field() {
        let lbm = AcousticLBM::new(8, 8, 1.0, 1.0, 1.0, 1.0);
        let ae = AcousticEnergy::new(1.0, 1.0 / 3.0_f64.sqrt());
        let total = ae.integrate_field(&lbm);
        assert!(total.abs() < 1e-20);
    }

    // --- FrequencyResponse tests ---

    #[test]
    fn test_frequency_response_dc_signal() {
        let signal = vec![1.0_f64; 64];
        let fr = FrequencyResponse::from_signal(&signal, 1.0);
        // DC component should be the mean
        assert!((fr.magnitude[0] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_frequency_response_zero_signal() {
        let signal = vec![0.0_f64; 64];
        let fr = FrequencyResponse::from_signal(&signal, 1.0);
        assert!(fr.magnitude.iter().all(|&m| m.abs() < 1e-15));
    }

    #[test]
    fn test_frequency_response_length() {
        let signal = vec![0.0_f64; 64];
        let fr = FrequencyResponse::from_signal(&signal, 1.0);
        assert_eq!(fr.n_fft, 33);
        assert_eq!(fr.magnitude.len(), 33);
    }

    #[test]
    fn test_frequency_response_dominant_frequency() {
        // Pure DC: dominant freq = 0
        let signal = vec![1.0_f64; 64];
        let fr = FrequencyResponse::from_signal(&signal, 1.0);
        assert!((fr.dominant_frequency()).abs() < 1e-10);
    }

    #[test]
    fn test_frequency_response_spl_db() {
        let signal = vec![1.0_f64; 64];
        let fr = FrequencyResponse::from_signal(&signal, 1.0);
        let spl = fr.spl_db(0, 1e-6);
        assert!(spl.is_finite());
        assert!(spl > 0.0);
    }

    // --- AcousticImpedance tests ---

    #[test]
    fn test_impedance_z_values() {
        let imp = AcousticImpedance::new(1.2, 343.0, 1000.0, 1500.0);
        assert!((imp.z1() - 1.2 * 343.0).abs() < 1e-10);
        assert!((imp.z2() - 1000.0 * 1500.0).abs() < 1e-10);
    }

    #[test]
    fn test_impedance_reflection_same_media() {
        let imp = AcousticImpedance::new(1.2, 343.0, 1.2, 343.0);
        assert!((imp.reflection_coefficient()).abs() < 1e-12);
    }

    #[test]
    fn test_impedance_transmission_same_media() {
        let imp = AcousticImpedance::new(1.2, 343.0, 1.2, 343.0);
        assert!((imp.transmission_coefficient() - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_impedance_energy_conservation() {
        let imp = AcousticImpedance::new(1.2, 343.0, 1000.0, 1500.0);
        let r_i = imp.intensity_reflection();
        let t_i = imp.intensity_transmission();
        // Note: intensity transmission is 1 - R^2 (energy conservation)
        assert!((r_i + t_i - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_impedance_reflection_hard_wall() {
        // Hard wall: Z2 >> Z1
        let imp = AcousticImpedance::new(1.2, 343.0, 1e12, 1e12);
        let r = imp.reflection_coefficient();
        assert!((r - 1.0).abs() < 1e-3, "Hard wall: R ≈ 1");
    }

    // --- FarFieldApproximation tests ---

    #[test]
    fn test_far_field_direction_unit_vector() {
        let ff = FarFieldApproximation::new([100.0, 0.0], 343.0, 1.2);
        let dir = ff.direction_to_observer([0.0, 0.0]);
        let mag = (dir[0] * dir[0] + dir[1] * dir[1]).sqrt();
        assert!((mag - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_far_field_distance() {
        let ff = FarFieldApproximation::new([3.0, 4.0], 343.0, 1.2);
        assert!((ff.distance([0.0, 0.0]) - 5.0).abs() < 1e-12);
    }

    #[test]
    fn test_far_field_monopole_pressure_positive() {
        let ff = FarFieldApproximation::new([100.0, 0.0], 343.0, 1.2);
        let p = ff.monopole_pressure(1.0, 100.0);
        assert!(p > 0.0);
    }

    #[test]
    fn test_far_field_integrate_surface_empty() {
        let ff = FarFieldApproximation::new([100.0, 0.0], 343.0, 1.2);
        let p_total = ff.integrate_surface(&[], &[], 1.0);
        assert_eq!(p_total, 0.0);
    }

    #[test]
    fn test_far_field_fwh_contribution_positive() {
        let ff = FarFieldApproximation::new([100.0, 0.0], 343.0, 1.2);
        let contrib = ff.fwh_surface_contribution(1.0, 1.0, 100.0, 1.0);
        assert!(contrib > 0.0);
    }

    // --- NoisePrediction tests ---

    #[test]
    fn test_noise_proudman_power_positive() {
        let np = NoisePrediction::new(1.2, 343.0, 0.1);
        assert!(np.proudman_power(1.0) > 0.0);
    }

    #[test]
    fn test_noise_proudman_scales_u5() {
        let np = NoisePrediction::new(1.2, 343.0, 0.1);
        let p1 = np.proudman_power(1.0);
        let p2 = np.proudman_power(2.0);
        assert!(
            (p2 / p1 - 32.0).abs() < 1e-10,
            "Proudman power scales as u^5"
        );
    }

    #[test]
    fn test_noise_peak_frequency() {
        let np = NoisePrediction::new(1.2, 343.0, 0.1);
        let f = np.peak_frequency(10.0);
        assert!((f - 0.2 * 10.0 / 0.1).abs() < 1e-10);
    }

    #[test]
    fn test_noise_oaspl_finite() {
        let np = NoisePrediction::new(1.2, 343.0, 0.1);
        let spl = np.oaspl(10.0, 2e-5);
        assert!(spl.is_finite());
    }

    #[test]
    fn test_noise_a_weighting_1khz() {
        // A-weighting at 1 kHz should be approximately -3 to +3 dB
        let aw = NoisePrediction::a_weighting_db(1000.0);
        assert!(aw.abs() < 10.0, "A-weighting at 1kHz should be small: {aw}");
    }

    #[test]
    fn test_noise_a_weighting_low_freq_negative() {
        // A-weighting at very low frequencies should be strongly negative
        let aw = NoisePrediction::a_weighting_db(10.0);
        assert!(
            aw < -20.0,
            "A-weighting at 10 Hz should be very negative: {aw}"
        );
    }
}
