// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Conserved phase-field crystal (PFC) model on a periodic 2-D grid.
//!
//! The model of Elder & Grant (Phys. Rev. E 70, 051605, 2004) describes a crystal
//! by a periodic order-parameter field `φ(r)` that minimises the Swift-Hohenberg
//! free energy
//!
//! ```text
//! F[φ] = ∫ dr [ ½ φ (r + (1 + ∇²)²) φ + ¼ φ⁴ ]
//! ```
//!
//! and relaxes by the *conserved* dynamics
//!
//! ```text
//! ∂φ/∂t = M ∇² [ (r + (1 + ∇²)²) φ + φ³ ] + ∇·ζ,    ⟨ζ_i ζ_j⟩ = 2 M T δ_ij δ(r − r') δ(t − t').
//! ```
//!
//! The equation is integrated pseudo-spectrally with a semi-implicit scheme: the
//! linear part is advanced implicitly and the cubic term explicitly, in Fourier
//! space,
//!
//! ```text
//! φ̂' = ( φ̂ − Δt M k² FFT(φ³) + Δt FFT(∇·ζ) ) / ( 1 + Δt M k² (r + (1 − k²)²) ),
//! ```
//!
//! which leaves the `k = 0` mode, and with it the mean of `φ`, exactly unchanged.

use oxifft::{Complex, irfft2d, rfft2d};
use oxiphysics_core::random::{Xoshiro256, box_muller};
use std::f64::consts::PI;

/// Default undercooling parameter `r` of the Swift-Hohenberg operator.
pub const DEFAULT_PFC_R: f64 = -0.25;

/// Default seed of the noise generator.
pub const DEFAULT_PFC_SEED: u64 = 0x5EED_C0FF_EE15_600D;

/// Signed frequency of bin `index` of an `n`-point transform.
fn signed_frequency(index: usize, n: usize) -> f64 {
    if index <= n / 2 {
        index as f64
    } else {
        index as f64 - n as f64
    }
}

/// Phase-field crystal (PFC) model on a 2D periodic regular grid.
///
/// The grid is row-major with `nx` rows (x) and `ny` columns (y) of spacing `dx`;
/// lengths are in units of the PFC lattice scale `1/k₀` (the preferred wave
/// number of the Swift-Hohenberg operator is `k₀ = 1`), so a triangular crystal has
/// the lattice constant `4π/√3`.
#[derive(Debug, Clone)]
pub struct PhaseFieldCrystal {
    /// Order parameter field (row-major, size `nx * ny`).
    pub phi: Vec<f64>,
    /// Number of grid points in x.
    pub nx: usize,
    /// Number of grid points in y.
    pub ny: usize,
    /// Grid spacing.
    pub dx: f64,
    /// Undercooling parameter `r` of the operator `r + (1 + ∇²)²`
    /// (default [`DEFAULT_PFC_R`]); more negative values favour the crystal.
    pub r: f64,
    /// Seed of the noise generator.  Every noisy [`PhaseFieldCrystal::step`] draws
    /// its noise from this seed and replaces it by the next one of the chain, so
    /// each step sees a fresh noise realisation and a run is reproducible from the
    /// initial seed.
    pub seed: u64,
}

impl PhaseFieldCrystal {
    /// Create a new `PhaseFieldCrystal` initialised to `phi_0` everywhere, with the
    /// default undercooling [`DEFAULT_PFC_R`] and noise seed [`DEFAULT_PFC_SEED`].
    pub fn new(nx: usize, ny: usize, dx: f64, phi_0: f64) -> Self {
        Self {
            phi: vec![phi_0; nx * ny],
            nx,
            ny,
            dx,
            r: DEFAULT_PFC_R,
            seed: DEFAULT_PFC_SEED,
        }
    }

    /// Replace the undercooling parameter `r`.
    pub fn with_r(mut self, r: f64) -> Self {
        self.r = r;
        self
    }

    /// Replace the noise seed.
    pub fn with_seed(mut self, seed: u64) -> Self {
        self.seed = seed;
        self
    }

    /// Linear index from 2D grid coordinates.
    fn idx(&self, i: usize, j: usize) -> usize {
        i * self.ny + j
    }

    /// Periodic neighbour in x.
    fn wrap_x(&self, i: i64) -> usize {
        i.rem_euclid(self.nx as i64) as usize
    }

    /// Periodic neighbour in y.
    fn wrap_y(&self, j: i64) -> usize {
        j.rem_euclid(self.ny as i64) as usize
    }

    /// Compute the 5-point discrete Laplacian at grid point `(i, j)`.
    pub fn laplacian(&self, i: usize, j: usize) -> f64 {
        let dx2 = self.dx * self.dx;
        let center = self.phi[self.idx(i, j)];
        let ip = self.phi[self.idx(self.wrap_x(i as i64 + 1), j)];
        let im = self.phi[self.idx(self.wrap_x(i as i64 - 1), j)];
        let jp = self.phi[self.idx(i, self.wrap_y(j as i64 + 1))];
        let jm = self.phi[self.idx(i, self.wrap_y(j as i64 - 1))];
        (ip + im + jp + jm - 4.0 * center) / dx2
    }

    /// Whether the field matches the grid dimensions and the grid is non-empty.
    fn is_consistent(&self) -> bool {
        self.nx > 0 && self.ny > 0 && self.phi.len() == self.nx * self.ny
    }

    /// Wave number `2π f / (n dx)` of the signed frequency `f`.
    fn wave_number(&self, frequency: f64, n: usize) -> f64 {
        2.0 * PI * frequency / (n as f64 * self.dx)
    }

    /// Squared wave number `k²` of the half-spectrum bin `(i, jh)` of an
    /// `rfft2d` of the `nx × ny` field.
    fn k_squared(&self, i: usize, jh: usize) -> f64 {
        let kx = self.wave_number(signed_frequency(i, self.nx), self.nx);
        let ky = self.wave_number(jh as f64, self.ny);
        kx * kx + ky * ky
    }

    /// Wave number used for a first derivative along x or y: the usual `2π f/(n dx)`,
    /// but zero for the Nyquist bin of an even-sized axis, whose derivative of a real
    /// field is not defined.
    fn derivative_wave_number(&self, frequency: f64, n: usize) -> f64 {
        if n.is_multiple_of(2) && frequency.abs() == (n / 2) as f64 {
            0.0
        } else {
            self.wave_number(frequency, n)
        }
    }

    /// Symbol `r + (1 − k²)²` of the Swift-Hohenberg operator at wave number `k`.
    fn linear_symbol(&self, k_squared: f64) -> f64 {
        self.r + (1.0 - k_squared) * (1.0 - k_squared)
    }

    /// Compute the total free energy of the current state.
    ///
    /// This is the Swift-Hohenberg functional
    /// `F = ∫ [ ½ φ (r + (1 + ∇²)²) φ + ¼ φ⁴ ] dA`, evaluated with the same spectral
    /// operator `r + (1 − k²)²` that drives [`PhaseFieldCrystal::step`] (Parseval's
    /// theorem for the quadratic part, a grid sum for the local quartic part), so
    /// that the dynamics is the gradient flow of exactly this energy.
    pub fn energy(&self) -> f64 {
        if !self.is_consistent() {
            return 0.0;
        }
        let n = self.nx * self.ny;
        let ny_half = self.ny / 2 + 1;
        let phi_hat = rfft2d(&self.phi, self.nx, self.ny);
        let mut quadratic = 0.0_f64;
        for i in 0..self.nx {
            for jh in 0..ny_half {
                // The half spectrum stores each conjugate pair once.
                let multiplicity = if jh == 0 || (self.ny.is_multiple_of(2) && jh == self.ny / 2) {
                    1.0
                } else {
                    2.0
                };
                let bin = phi_hat[i * ny_half + jh];
                let symbol = self.linear_symbol(self.k_squared(i, jh));
                quadratic += multiplicity * symbol * bin.norm_sqr();
            }
        }
        let quadratic = 0.5 * quadratic / n as f64;
        let quartic = 0.25 * self.phi.iter().map(|p| p.powi(4)).sum::<f64>();
        (quadratic + quartic) * self.dx * self.dx
    }

    /// Spectrum (rfft2d half layout) of the conserved thermal noise `∇·ζ`.
    ///
    /// `ζ_x` and `ζ_y` are independent Gaussian white-noise fields with the
    /// fluctuation-dissipation variance `2 M T / (dx² Δt)` per cell and step.  Both
    /// are drawn from the generator seeded with `self.seed`, which is advanced to the
    /// next seed of the chain afterwards.  The divergence is taken in Fourier space,
    /// so it has no `k = 0` component: the mean of `φ` stays exact.
    fn conserved_noise_spectrum(
        &mut self,
        dt: f64,
        mobility: f64,
        noise_temperature: f64,
    ) -> Vec<Complex<f64>> {
        let n = self.nx * self.ny;
        let sigma = (2.0 * mobility * noise_temperature / (self.dx * self.dx * dt)).sqrt();
        let mut rng = Xoshiro256::new(self.seed);
        let gaussian_field = |rng: &mut Xoshiro256| -> Vec<f64> {
            let mut field = Vec::with_capacity(n + 1);
            while field.len() < n {
                let (a, b) = box_muller(rng);
                field.push(sigma * a);
                field.push(sigma * b);
            }
            field.truncate(n);
            field
        };
        let zeta_x = gaussian_field(&mut rng);
        let zeta_y = gaussian_field(&mut rng);
        self.seed = rng.next_u64();

        let hat_x = rfft2d(&zeta_x, self.nx, self.ny);
        let hat_y = rfft2d(&zeta_y, self.nx, self.ny);
        let ny_half = self.ny / 2 + 1;
        let mut spectrum = vec![Complex::zero(); self.nx * ny_half];
        for i in 0..self.nx {
            let kx = self.derivative_wave_number(signed_frequency(i, self.nx), self.nx);
            for jh in 0..ny_half {
                let ky = self.derivative_wave_number(jh as f64, self.ny);
                let flat = i * ny_half + jh;
                // i k·ζ̂ = i (kx ζ̂x + ky ζ̂y)
                spectrum[flat] = Complex::new(
                    -(kx * hat_x[flat].im + ky * hat_y[flat].im),
                    kx * hat_x[flat].re + ky * hat_y[flat].re,
                );
            }
        }
        spectrum
    }

    /// Advance the PFC order parameter by one time step of the conserved
    /// Elder-Grant PFC equation with a semi-implicit pseudo-spectral scheme.
    ///
    /// `∂φ/∂t = M ∇² [ (r + (1 + ∇²)²) φ + φ³ ] + ∇·ζ` is advanced as
    /// `φ̂' = (φ̂ − Δt M k² FFT(φ³) + Δt FFT(∇·ζ)) / (1 + Δt M k² (r + (1 − k²)²))`.
    /// The mean of `φ` is conserved to rounding, with and without noise.  The
    /// implicit linear part keeps the scheme stable for stiff `k⁶` modes; the
    /// explicit cubic term and the unstable band around `k = 1` require
    /// `Δt M |r| ≲ 1` and a `Δt` of order one at most.
    ///
    /// # Arguments
    /// * `dt`         - time step (must be positive)
    /// * `mobility`   - kinetic mobility M (must be non-negative)
    /// * `temp_noise` - noise temperature `T` of the conserved thermal noise
    ///   `⟨ζ_i ζ_j⟩ = 2 M T δ_ij δ δ` (set to 0 for deterministic evolution)
    ///
    /// A non-positive or non-finite `dt`, a negative or non-finite `mobility`, a zero
    /// `mobility` and an empty or inconsistent grid leave the field unchanged.
    pub fn step(&mut self, dt: f64, mobility: f64, temp_noise: f64) {
        let valid_step = dt.is_finite() && dt > 0.0 && mobility.is_finite() && mobility >= 0.0;
        // Zero mobility freezes the dynamics (and scales the noise to zero as well).
        if !self.is_consistent() || !valid_step || mobility == 0.0 {
            return;
        }
        let ny_half = self.ny / 2 + 1;
        let cubic: Vec<f64> = self.phi.iter().map(|p| p * p * p).collect();
        let phi_hat = rfft2d(&self.phi, self.nx, self.ny);
        let cubic_hat = rfft2d(&cubic, self.nx, self.ny);
        let noise_hat = if temp_noise.is_finite() && temp_noise > 0.0 && mobility > 0.0 {
            Some(self.conserved_noise_spectrum(dt, mobility, temp_noise))
        } else {
            None
        };

        let mut next = vec![Complex::zero(); self.nx * ny_half];
        for i in 0..self.nx {
            for jh in 0..ny_half {
                let flat = i * ny_half + jh;
                let k2 = self.k_squared(i, jh);
                let decay = dt * mobility * k2;
                let mut re = phi_hat[flat].re - decay * cubic_hat[flat].re;
                let mut im = phi_hat[flat].im - decay * cubic_hat[flat].im;
                if let Some(noise) = &noise_hat {
                    re += dt * noise[flat].re;
                    im += dt * noise[flat].im;
                }
                let denominator = 1.0 + decay * self.linear_symbol(k2);
                next[flat] = Complex::new(re / denominator, im / denominator);
            }
        }
        // `irfft2d` is already normalised by the number of grid points.
        self.phi = irfft2d(&next, self.nx, self.ny);
    }

    /// Return the mean order parameter over the grid (`0.0` for an empty grid).
    pub fn mean_phi(&self) -> f64 {
        if self.phi.is_empty() {
            return 0.0;
        }
        self.phi.iter().sum::<f64>() / self.phi.len() as f64
    }
}
