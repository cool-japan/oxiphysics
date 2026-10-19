// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Quantum-inspired LBM extensions for superfluids and quantum turbulence.
//!
//! # Overview
//!
//! This module provides tools for quantum fluid simulation via Lattice Boltzmann
//! methods including:
//!
//! - [`QuantumLbmParams`] — physical parameters for quantum LBM
//! - [`GrossPitaevskiiLbm`] — 1-D Gross-Pitaevskii equation for a superfluid BEC
//!   (Strang split-step Fourier integrator, exactly number conserving)
//! - [`QuantumVortex`] — quantized vortex with Biot-Savart velocity field
//! - [`VortexDynamics`] — multi-vortex Biot-Savart interactions
//! - [`SuperfluidLbm`] — two-fluid model: D1Q3 BGK lattice Boltzmann normal fluid,
//!   inviscid (Euler) superfluid and momentum-conserving mutual friction
//! - [`QuantumTurbulence`] — vortex tangle and Kolmogorov spectrum
//! - [`SchroedingerLbm`] — LBM solution of Schrödinger equation (Madelung)
//! - [`BoseEinsteinCondensate`] — BEC with Thomas-Fermi profile and solitons
//! - [`QuantumPressure`] — Bohm quantum potential Q = -ħ²∇²√ρ / (2m√ρ)
//!
//! # References
//! - Succi, S. (2015). Quantum Lattice Boltzmann. *Phil. Trans. R. Soc. A*, 373.
//! - Gross, E.P. (1961). *Nuovo Cimento*, 20, 454.
//! - Pitaevskii, L.P. (1961). *Sov. Phys. JETP*, 13, 451.

use crate::error::Error;
use oxifft::Complex;
use std::f64::consts::PI;

mod superfluid;

pub use superfluid::SuperfluidLbm;

// ---------------------------------------------------------------------------
// Physical constants
// ---------------------------------------------------------------------------

/// Reduced Planck constant ħ (J·s).
pub const HBAR: f64 = 1.054_571_817e-34;
/// Planck constant h (J·s).
pub const PLANCK_H: f64 = 6.626_070_15e-34;
/// Boltzmann constant k_B (J K⁻¹).
pub const K_B: f64 = 1.380_649e-23;
/// Mass of a ⁴He atom (kg).
pub const HELIUM4_MASS: f64 = 6.646_477_208_8e-27;
/// Quantum of circulation κ = h/m₄He (m² s⁻¹).
pub const QUANTUM_CIRCULATION_HE4: f64 = PLANCK_H / HELIUM4_MASS;

// ---------------------------------------------------------------------------
// Utility free functions
// ---------------------------------------------------------------------------

/// Compute the healing (coherence) length ξ = ħ / √(2m·g·n).
///
/// This is the characteristic length scale over which the superfluid order
/// parameter recovers from a perturbation.
///
/// # Arguments
/// * `hbar`    – reduced Planck constant (J·s)
/// * `mass`    – particle mass (kg)
/// * `g_n`     – coupling × density (g·n) (J m⁻³)
///
/// # Returns
/// Healing length in metres; returns 0 when `g_n ≤ 0`.
pub fn healing_length(hbar: f64, mass: f64, g_n: f64) -> f64 {
    if g_n <= 0.0 || mass <= 0.0 {
        return 0.0;
    }
    hbar / (2.0 * mass * g_n).sqrt()
}

/// Compute the coherence length ξ_c = ħ / √(2m·μ).
///
/// For a BEC with chemical potential μ this length sets the Thomas-Fermi radius.
///
/// # Arguments
/// * `hbar`  – reduced Planck constant (J·s)
/// * `mass`  – particle mass (kg)
/// * `mu`    – chemical potential (J)
pub fn coherence_length(hbar: f64, mass: f64, mu: f64) -> f64 {
    if mu <= 0.0 || mass <= 0.0 {
        return 0.0;
    }
    hbar / (2.0 * mass * mu).sqrt()
}

/// Compute the vortex core radius r_c ≈ ξ / √2.
///
/// The vortex core is the region of suppressed density around a quantized
/// vortex filament; its radius is approximately the healing length / √2.
///
/// # Arguments
/// * `xi` – healing length (m)
pub fn vortex_core_radius(xi: f64) -> f64 {
    xi / 2.0_f64.sqrt()
}

// ---------------------------------------------------------------------------
// QuantumLbmParams
// ---------------------------------------------------------------------------

/// Physical parameters shared by all quantum-LBM solvers.
///
/// Holds the fundamental length/energy scales that appear in the
/// Gross-Pitaevskii and Madelung equations.
#[derive(Debug, Clone)]
pub struct QuantumLbmParams {
    /// Reduced Planck constant ħ (J·s or dimensionless lattice units).
    pub hbar: f64,
    /// Particle mass (kg or lattice units).
    pub mass: f64,
    /// Healing length ξ (m or lattice units).
    pub healing_length: f64,
    /// Contact interaction coupling constant g (J m³ or lattice units).
    pub coupling: f64,
}

impl QuantumLbmParams {
    /// Construct a new parameter set.
    ///
    /// # Arguments
    /// * `hbar`           – reduced Planck constant
    /// * `mass`           – particle mass
    /// * `healing_length` – healing length ξ
    /// * `coupling`       – interaction coupling constant g
    pub fn new(hbar: f64, mass: f64, healing_length: f64, coupling: f64) -> Self {
        Self {
            hbar,
            mass,
            healing_length,
            coupling,
        }
    }

    /// Build from physical inputs: compute healing length from coupling and density.
    ///
    /// # Arguments
    /// * `hbar`    – reduced Planck constant
    /// * `mass`    – particle mass
    /// * `coupling`– interaction coupling g
    /// * `density` – background number density n₀
    pub fn from_physics(hbar: f64, mass: f64, coupling: f64, density: f64) -> Self {
        let g_n = coupling * density;
        let xi = healing_length(hbar, mass, g_n);
        Self {
            hbar,
            mass,
            healing_length: xi,
            coupling,
        }
    }

    /// Quantum of circulation κ = h/m = 2πħ/m.
    pub fn kappa(&self) -> f64 {
        2.0 * PI * self.hbar / self.mass
    }

    /// Sound speed of a weakly-interacting BEC: c_s = √(g·n / m).
    ///
    /// # Arguments
    /// * `density` – mean number density n₀
    pub fn sound_speed(&self, density: f64) -> f64 {
        if density <= 0.0 {
            return 0.0;
        }
        (self.coupling * density / self.mass).sqrt()
    }
}

// ---------------------------------------------------------------------------
// GrossPitaevskiiLbm
// ---------------------------------------------------------------------------

/// Wrap a phase difference into `(-π, π]`.
fn wrap_phase(delta: f64) -> f64 {
    delta - 2.0 * PI * (delta / (2.0 * PI)).round()
}

/// Solver for the one-dimensional Gross-Pitaevskii (GP) equation
///
/// ```text
/// i ħ ∂ψ/∂t = [ -ħ²/(2m) ∂²/∂x² + V(x) + g |ψ|² ] ψ
/// ```
///
/// on a periodic lattice of `nx` sites with spacing `dx`.
///
/// The condensate wavefunction ψ(x) = A(x)·exp(iθ(x)) is stored as the
/// amplitude `A` and the phase `θ` (wrapped to `(-π, π]`).  The density is
/// ρ = |ψ|² and the superfluid velocity is v_s = ħ∇θ/m.
///
/// # Time integration
///
/// One [`step`](Self::step) is a second-order Strang split-step Fourier step,
///
/// ```text
/// ψ ← exp(-i W dt/2ħ) · IFFT[ exp(-i ħ k² dt/2m) · FFT[ ψ ] ] · exp(-i W dt/2ħ),
/// W = V + g|ψ|²,
/// ```
///
/// where the potential/interaction half-steps are exact phase rotations (they
/// leave |ψ|² unchanged) and the kinetic step is exact in Fourier space.  Every
/// factor is unitary, so the particle number `N = Σ|ψ|² dx` is conserved to
/// rounding error, a uniform condensate evolves exactly as
/// `θ(t) = -g ρ t / ħ`, and the scheme is unconditionally stable; the error is
/// `O(dt²)` and the spatial discretisation is spectrally accurate for fields
/// resolved by the lattice.  The FFTs come from `oxifft`.
///
/// The type name is historical: an earlier version relaxed a D1Q3 lattice
/// Boltzmann distribution whose moments were inconsistent with the GP equation
/// and therefore did not solve it.  No lattice populations are involved any
/// more.
#[derive(Debug, Clone)]
pub struct GrossPitaevskiiLbm {
    /// Number of lattice sites.
    pub nx: usize,
    /// Amplitude |ψ(x)| at each site.
    pub amplitude: Vec<f64>,
    /// Phase θ(x) at each site (radians).
    pub phase: Vec<f64>,
    /// Physical parameters (ħ, particle mass m and contact coupling g).
    pub params: QuantumLbmParams,
    /// Lattice spacing Δx.
    pub dx: f64,
    /// Time step Δt.
    pub dt: f64,
    /// External potential V(x), one value per site.
    pub potential: Vec<f64>,
}

impl GrossPitaevskiiLbm {
    /// Create a new lattice with uniform unit density, zero phase and no
    /// external potential.
    ///
    /// # Arguments
    /// * `nx`     – number of lattice sites
    /// * `params` – quantum parameters (ħ, m, g)
    /// * `dx`     – lattice spacing
    /// * `dt`     – time step
    pub fn new(nx: usize, params: QuantumLbmParams, dx: f64, dt: f64) -> Self {
        Self {
            nx,
            amplitude: vec![1.0; nx],
            phase: vec![0.0; nx],
            params,
            dx,
            dt,
            potential: vec![0.0; nx],
        }
    }

    /// Set the wavefunction from amplitude and phase arrays.
    ///
    /// # Panics
    ///
    /// Panics if either array does not have `nx` entries.
    pub fn set_wavefunction(&mut self, amplitude: Vec<f64>, phase: Vec<f64>) {
        assert_eq!(amplitude.len(), self.nx);
        assert_eq!(phase.len(), self.nx);
        self.amplitude = amplitude;
        self.phase = phase;
    }

    /// Return density ρ(x) = |ψ(x)|².
    pub fn density(&self) -> Vec<f64> {
        self.amplitude.iter().map(|a| a * a).collect()
    }

    /// Return the superfluid velocity v_s = ħ∇θ/m.
    ///
    /// The phase gradient is the central difference of the phase with every
    /// difference wrapped into `(-π, π]`, so the 2π jumps of a wrapped phase
    /// do not create spurious velocities.  It is exact for phase differences
    /// below π between next-nearest sites, i.e. for fields resolved by the
    /// lattice.
    pub fn superfluid_velocity(&self) -> Vec<f64> {
        let n = self.nx;
        let prefactor = self.params.hbar / self.params.mass;
        (0..n)
            .map(|i| {
                let im = (i + n - 1) % n;
                let ip = (i + 1) % n;
                let dphi = wrap_phase(self.phase[ip] - self.phase[im]) / (2.0 * self.dx);
                prefactor * dphi
            })
            .collect()
    }

    /// Compute total particle number N = ∫|ψ|² dx ≈ Σ ρᵢ · Δx.
    pub fn total_number(&self) -> f64 {
        self.amplitude.iter().map(|a| a * a * self.dx).sum()
    }

    /// Complex wavefunction ψ = A·exp(iθ) at every site.
    fn wavefunction(&self) -> Vec<Complex<f64>> {
        self.amplitude
            .iter()
            .zip(&self.phase)
            .map(|(&a, &theta)| Complex::from_polar(a, theta))
            .collect()
    }

    /// Angular wavenumber k_j of Fourier mode `j` of an `nx`-point transform
    /// (negative for the upper half of the spectrum).
    fn wavenumber(&self, j: usize) -> f64 {
        let n = self.nx;
        let signed = if j <= n / 2 {
            j as f64
        } else {
            j as f64 - n as f64
        };
        2.0 * PI * signed / (n as f64 * self.dx)
    }

    /// Check that the lattice is well formed and the parameters are physical.
    fn validate(&self) -> Result<(), Error> {
        if self.nx == 0 {
            return Err(Error::InvalidGridDimensions { nx: 0, ny: 1 });
        }
        for (field, len) in [
            ("amplitude", self.amplitude.len()),
            ("phase", self.phase.len()),
            ("potential", self.potential.len()),
        ] {
            if len != self.nx {
                return Err(Error::ComponentCountMismatch {
                    expected: self.nx,
                    found: len,
                    field: field.to_string(),
                });
            }
        }
        for (name, value) in [
            ("hbar", self.params.hbar),
            ("mass", self.params.mass),
            ("dx", self.dx),
        ] {
            if !value.is_finite() || value <= 0.0 {
                return Err(Error::ParameterOutOfBounds {
                    name: name.to_string(),
                    value,
                    min: f64::MIN_POSITIVE,
                    max: f64::MAX,
                });
            }
        }
        for (name, value) in [("dt", self.dt), ("coupling", self.params.coupling)] {
            if !value.is_finite() {
                return Err(Error::ParameterOutOfBounds {
                    name: name.to_string(),
                    value,
                    min: f64::MIN,
                    max: f64::MAX,
                });
            }
        }
        Ok(())
    }

    /// Multiply ψ by `exp(-i (V + g|ψ|²) dt / (2ħ))`, half of the
    /// potential/interaction phase rotation of a Strang step.
    fn apply_potential_half_step(&self, psi: &mut [Complex<f64>]) {
        let scale = 0.5 * self.dt / self.params.hbar;
        let g = self.params.coupling;
        for (z, &v) in psi.iter_mut().zip(&self.potential) {
            let w = v + g * z.norm_sqr();
            *z *= Complex::cis(-w * scale);
        }
    }

    /// Advance the wavefunction by one time step `dt` (see the type docs).
    ///
    /// # Errors
    ///
    /// * [`Error::InvalidGridDimensions`] for an empty lattice,
    /// * [`Error::ComponentCountMismatch`] if `amplitude`, `phase` or
    ///   `potential` do not have `nx` entries,
    /// * [`Error::ParameterOutOfBounds`] if ħ, m or Δx are not positive and
    ///   finite, or Δt / g are not finite.
    ///
    /// The state is not modified when an error is returned.
    pub fn try_step(&mut self) -> Result<(), Error> {
        self.validate()?;
        let mut psi = self.wavefunction();

        // Half step of the potential and the mean-field interaction.
        self.apply_potential_half_step(&mut psi);

        // Kinetic step: exact free propagation in Fourier space.
        let mut spectrum = oxifft::fft(&psi);
        let kinetic_scale = 0.5 * self.params.hbar * self.dt / self.params.mass;
        for (j, c) in spectrum.iter_mut().enumerate() {
            let k = self.wavenumber(j);
            *c *= Complex::cis(-kinetic_scale * k * k);
        }
        psi = oxifft::ifft(&spectrum);

        // Second half step, with the density after the kinetic step.
        self.apply_potential_half_step(&mut psi);

        for ((a, theta), z) in self
            .amplitude
            .iter_mut()
            .zip(self.phase.iter_mut())
            .zip(&psi)
        {
            *a = z.norm();
            *theta = z.arg();
        }
        Ok(())
    }

    /// Advance the wavefunction by one time step.
    ///
    /// Equivalent to [`try_step`](Self::try_step); an invalid configuration
    /// leaves the state untouched.
    pub fn step(&mut self) {
        // The configuration error, if any, is reported by `try_step`.
        let _ = self.try_step();
    }

    /// Total energy
    /// `E = ∫ [ ħ²/(2m)|∂ₓψ|² + V|ψ|² + (g/2)|ψ|⁴ ] dx`.
    ///
    /// The gradient term is evaluated spectrally.  `E` is conserved by the
    /// exact dynamics; the split-step integrator keeps it constant up to
    /// `O(dt²)` oscillations.  Returns `0.0` for an empty lattice.
    pub fn total_energy(&self) -> f64 {
        let n = self.nx;
        if n == 0 || self.potential.len() != n {
            return 0.0;
        }
        let psi = self.wavefunction();
        let spectrum = oxifft::fft(&psi);
        // Parseval: Σ_j |∂ₓψ_j|² Δx = (Δx/N) Σ_k k² |ψ̂_k|².
        let gradient: f64 = spectrum
            .iter()
            .enumerate()
            .map(|(j, c)| {
                let k = self.wavenumber(j);
                k * k * c.norm_sqr()
            })
            .sum::<f64>()
            * self.dx
            / n as f64;
        let kinetic = 0.5 * self.params.hbar * self.params.hbar / self.params.mass * gradient;
        let local: f64 = psi
            .iter()
            .zip(&self.potential)
            .map(|(z, &v)| {
                let rho = z.norm_sqr();
                (v + 0.5 * self.params.coupling * rho) * rho
            })
            .sum::<f64>()
            * self.dx;
        kinetic + local
    }
}

// ---------------------------------------------------------------------------
// QuantumVortex
// ---------------------------------------------------------------------------

/// A single quantized vortex with position and circulation quantum.
///
/// Carries circulation Γ = n·κ where κ = h/m and n is the winding number.
#[derive(Debug, Clone)]
pub struct QuantumVortex {
    /// x-position of the vortex core (m).
    pub x: f64,
    /// y-position of the vortex core (m).
    pub y: f64,
    /// Winding number n (signed integer, ±1 for singly-quantized vortex).
    pub winding: i32,
    /// Mass of the superfluid particle (kg).
    pub mass: f64,
}

impl QuantumVortex {
    /// Create a new quantum vortex.
    ///
    /// # Arguments
    /// * `x`, `y`  – 2-D core position (m)
    /// * `winding` – topological winding number (typically ±1)
    /// * `mass`    – particle mass (kg)
    pub fn new(x: f64, y: f64, winding: i32, mass: f64) -> Self {
        Self {
            x,
            y,
            winding,
            mass,
        }
    }

    /// Quantum of circulation κ = h/m for this vortex type.
    pub fn kappa(&self) -> f64 {
        PLANCK_H / self.mass
    }

    /// Total circulation Γ = n·κ.
    pub fn circulation(&self) -> f64 {
        self.winding as f64 * self.kappa()
    }

    /// Induced velocity (vx, vy) at position (px, py) via Biot-Savart (2-D).
    ///
    /// v = Γ/(2π r²) × r̂_perp  where r = (px-x, py-y)
    ///
    /// Returns `[0.0, 0.0]` if the evaluation point is at the vortex core.
    pub fn velocity_at(&self, px: f64, py: f64) -> [f64; 2] {
        let dx = px - self.x;
        let dy = py - self.y;
        let r2 = dx * dx + dy * dy;
        if r2 < 1e-30 {
            return [0.0, 0.0];
        }
        let gamma = self.circulation();
        // v = Γ/(2π) · (-dy, dx) / r²
        let factor = gamma / (2.0 * PI * r2);
        [-factor * dy, factor * dx]
    }
}

// ---------------------------------------------------------------------------
// VortexDynamics
// ---------------------------------------------------------------------------

/// Multi-vortex Biot-Savart dynamics for a 2-D array of quantum vortices.
///
/// Implements point-vortex evolution: each vortex moves with the velocity
/// induced by all other vortices (mutual induction approximation).
#[derive(Debug, Clone)]
pub struct VortexDynamics {
    /// List of quantum vortices.
    pub vortices: Vec<QuantumVortex>,
    /// Time step for Euler integration.
    pub dt: f64,
    /// Cutoff radius to avoid divergence (m).
    pub cutoff: f64,
}

impl VortexDynamics {
    /// Create a new vortex dynamics solver.
    ///
    /// # Arguments
    /// * `vortices` – initial vortex configuration
    /// * `dt`       – integration time step (s)
    /// * `cutoff`   – minimum distance for Biot-Savart (m)
    pub fn new(vortices: Vec<QuantumVortex>, dt: f64, cutoff: f64) -> Self {
        Self {
            vortices,
            dt,
            cutoff,
        }
    }

    /// Compute the velocity at (px, py) due to all vortices.
    pub fn velocity_at(&self, px: f64, py: f64) -> [f64; 2] {
        let mut vx = 0.0f64;
        let mut vy = 0.0f64;
        for v in &self.vortices {
            let dx = px - v.x;
            let dy = py - v.y;
            let r2 = (dx * dx + dy * dy).max(self.cutoff * self.cutoff);
            let gamma = v.circulation();
            let factor = gamma / (2.0 * PI * r2);
            vx += -factor * dy;
            vy += factor * dx;
        }
        [vx, vy]
    }

    /// Advance all vortices one time step (forward Euler).
    pub fn step(&mut self) {
        let n = self.vortices.len();
        let mut vels = vec![[0.0f64; 2]; n];
        for (i, vel_i) in vels.iter_mut().enumerate() {
            let px = self.vortices[i].x;
            let py = self.vortices[i].y;
            for j in 0..n {
                if i == j {
                    continue;
                }
                let v = &self.vortices[j];
                let dx = px - v.x;
                let dy = py - v.y;
                let r2 = (dx * dx + dy * dy).max(self.cutoff * self.cutoff);
                let gamma = v.circulation();
                let factor = gamma / (2.0 * PI * r2);
                vel_i[0] += -factor * dy;
                vel_i[1] += factor * dx;
            }
        }
        for (i, v) in self.vortices.iter_mut().enumerate() {
            v.x += vels[i][0] * self.dt;
            v.y += vels[i][1] * self.dt;
        }
    }

    /// Total enstrophy Σ Γᵢ².
    pub fn enstrophy(&self) -> f64 {
        self.vortices.iter().map(|v| v.circulation().powi(2)).sum()
    }

    /// Number of vortex-antivortex pairs (count vortices with +1 and -1).
    pub fn pair_count(&self) -> usize {
        let pos: usize = self.vortices.iter().filter(|v| v.winding > 0).count();
        let neg: usize = self.vortices.iter().filter(|v| v.winding < 0).count();
        pos.min(neg)
    }
}

// ---------------------------------------------------------------------------
// SuperfluidLbm
// ---------------------------------------------------------------------------

// The two-fluid model lives in `superfluid.rs`; it is re-exported at the top of
// this file so that `quantum_lbm::SuperfluidLbm` remains the public path.

/// Compute the superfluid fraction using the two-fluid model approximation.
///
/// Uses the empirical formula f_s = 1 - (T/T_c)^4 for T < T_c.
///
/// # Arguments
/// * `temperature` – current temperature (K)
/// * `t_critical`  – lambda-point temperature (K)
pub fn superfluid_fraction(temperature: f64, t_critical: f64) -> f64 {
    if temperature <= 0.0 || t_critical <= 0.0 {
        return 1.0;
    }
    if temperature >= t_critical {
        return 0.0;
    }
    let reduced_t = temperature / t_critical;
    1.0 - reduced_t.powi(4)
}

// ---------------------------------------------------------------------------
// QuantumTurbulence
// ---------------------------------------------------------------------------

/// Quantum turbulence model: vortex tangle with Kolmogorov-type spectrum.
///
/// Maintains a vortex tangle characterized by its line density L (m⁻²)
/// and simulates energy cascade and reconnection statistics.
#[derive(Debug, Clone)]
pub struct QuantumTurbulence {
    /// Vortex line density L (m⁻²).
    pub vortex_line_density: f64,
    /// Quantum of circulation κ (m² s⁻¹).
    pub kappa: f64,
    /// Temperature (K) for mutual friction.
    pub temperature: f64,
    /// Critical temperature T_c (K).
    pub t_critical: f64,
    /// Mutual friction parameter α.
    pub alpha: f64,
    /// Reconnection event count.
    pub reconnection_count: u64,
    /// Energy dissipation rate ε (m² s⁻³).
    pub dissipation_rate: f64,
}

impl QuantumTurbulence {
    /// Create a new quantum turbulence state.
    ///
    /// # Arguments
    /// * `vortex_line_density` – initial L (m⁻²)
    /// * `kappa`               – quantum of circulation (m² s⁻¹)
    /// * `temperature`         – temperature (K)
    /// * `t_critical`          – lambda point (K)
    /// * `alpha`               – mutual friction coefficient
    pub fn new(
        vortex_line_density: f64,
        kappa: f64,
        temperature: f64,
        t_critical: f64,
        alpha: f64,
    ) -> Self {
        Self {
            vortex_line_density,
            kappa,
            temperature,
            t_critical,
            alpha,
            reconnection_count: 0,
            dissipation_rate: 0.0,
        }
    }

    /// Mean inter-vortex spacing δ ≈ 1/√L.
    pub fn inter_vortex_spacing(&self) -> f64 {
        if self.vortex_line_density <= 0.0 {
            return f64::INFINITY;
        }
        1.0 / self.vortex_line_density.sqrt()
    }

    /// Kolmogorov wavenumber k_K = (ε/ν³)^(1/4) (quantum analogue).
    ///
    /// Uses ε = α·κ³·L² and effective viscosity ν_eff = κ/4π.
    pub fn kolmogorov_wavenumber(&self) -> f64 {
        let eps = self.dissipation_rate.max(1e-30);
        let nu = self.kappa / (4.0 * PI);
        (eps / nu.powi(3)).powf(0.25)
    }

    /// Vinen equation: time evolution of vortex line density.
    ///
    /// dL/dt = α·κ·L² - χ₂·κ·L^(3/2)
    ///
    /// # Arguments
    /// * `chi2` – dimensionless Vinen coefficient
    /// * `dt`   – time step (s)
    pub fn vinen_step(&mut self, chi2: f64, dt: f64) {
        let l = self.vortex_line_density;
        let alpha = self.alpha;
        let kappa = self.kappa;
        let dl = alpha * kappa * l * l - chi2 * kappa * l.powf(1.5);
        self.vortex_line_density = (l + dl * dt).max(0.0);
        // update dissipation
        self.dissipation_rate = alpha * kappa.powi(3) * self.vortex_line_density.powi(2);
    }

    /// Register a reconnection event.
    pub fn register_reconnection(&mut self) {
        self.reconnection_count += 1;
    }

    /// Energy spectrum E(k) ~ ε^(2/3) k^(-5/3) at wavenumber k (Kolmogorov).
    pub fn kolmogorov_spectrum(&self, k: f64) -> f64 {
        if k <= 0.0 {
            return 0.0;
        }
        let c_k = 1.5; // Kolmogorov constant
        let eps = self.dissipation_rate.max(1e-30);
        c_k * eps.powf(2.0 / 3.0) * k.powf(-5.0 / 3.0)
    }
}

// ---------------------------------------------------------------------------
// SchroedingerLbm
// ---------------------------------------------------------------------------

/// LBM solution of the Schrödinger equation via the Madelung transform.
///
/// The Madelung transform ψ = √ρ · exp(iS/ħ) maps the Schrödinger equation
/// to hydrodynamic equations: continuity + momentum with quantum pressure.
/// This struct evolves amplitude √ρ and phase S on a 1-D periodic lattice.
#[derive(Debug, Clone)]
pub struct SchroedingerLbm {
    /// Number of lattice sites.
    pub nx: usize,
    /// Amplitude √ρ(x).
    pub amplitude: Vec<f64>,
    /// Reduced action S(x)/ħ = phase θ(x).
    pub phase: Vec<f64>,
    /// External potential V(x).
    pub potential: Vec<f64>,
    /// Particle mass m.
    pub mass: f64,
    /// Reduced Planck constant ħ.
    pub hbar: f64,
    /// Lattice spacing Δx.
    pub dx: f64,
    /// Time step Δt.
    pub dt: f64,
}

impl SchroedingerLbm {
    /// Create a Schrödinger LBM lattice.
    ///
    /// # Arguments
    /// * `nx`   – number of lattice sites
    /// * `mass` – particle mass
    /// * `hbar` – reduced Planck constant
    /// * `dx`   – lattice spacing
    /// * `dt`   – time step
    pub fn new(nx: usize, mass: f64, hbar: f64, dx: f64, dt: f64) -> Self {
        Self {
            nx,
            amplitude: vec![1.0; nx],
            phase: vec![0.0; nx],
            potential: vec![0.0; nx],
            mass,
            hbar,
            dx,
            dt,
        }
    }

    /// Set wavefunction from real and imaginary parts: ψ = re + i·im.
    pub fn set_from_ri(&mut self, re: &[f64], im: &[f64]) {
        for i in 0..self.nx {
            let a = (re[i] * re[i] + im[i] * im[i]).sqrt();
            self.amplitude[i] = a;
            self.phase[i] = im[i].atan2(re[i]);
        }
    }

    /// Probability density ρ(x) = |ψ|².
    pub fn probability_density(&self) -> Vec<f64> {
        self.amplitude.iter().map(|a| a * a).collect()
    }

    /// Total probability ∫|ψ|² dx (should equal 1 for normalized state).
    pub fn norm(&self) -> f64 {
        self.amplitude.iter().map(|a| a * a * self.dx).sum()
    }

    /// Advance one step via split-operator Madelung LBM.
    pub fn step(&mut self) {
        let nx = self.nx;
        let hbar = self.hbar;
        let mass = self.mass;
        let dx = self.dx;
        let dt = self.dt;

        // quantum pressure and phase evolution
        for i in 0..nx {
            let im_idx = if i == 0 { nx - 1 } else { i - 1 };
            let ip_idx = if i == nx - 1 { 0 } else { i + 1 };
            let a = self.amplitude[i].max(1e-30);
            let a_m = self.amplitude[im_idx].max(1e-30);
            let a_p = self.amplitude[ip_idx].max(1e-30);
            // Laplacian of amplitude
            let lap_a = (a_p - 2.0 * a + a_m) / (dx * dx);
            // quantum potential Q = -ħ²/(2m) · ∇²A/A
            let q = -hbar * hbar / (2.0 * mass) * lap_a / a;
            // phase evolution: ∂S/∂t = -V - Q - (∇S)²/(2m)
            let dphi = (self.phase[ip_idx] - self.phase[im_idx]) / (2.0 * dx);
            let kinetic = hbar * hbar * dphi * dphi / (2.0 * mass);
            let dphase = -(self.potential[i] + q + kinetic) / hbar * dt;
            self.phase[i] += dphase;
        }

        // amplitude evolution via continuity
        let mut new_amp = self.amplitude.clone();
        for (i, new_amp_i) in new_amp.iter_mut().enumerate() {
            let im_idx = if i == 0 { nx - 1 } else { i - 1 };
            let ip_idx = if i == nx - 1 { 0 } else { i + 1 };
            let dphi_p = (self.phase[ip_idx] - self.phase[i]) / dx;
            let dphi_m = (self.phase[i] - self.phase[im_idx]) / dx;
            let a_p = self.amplitude[ip_idx];
            let a_m = self.amplitude[im_idx];
            let a = self.amplitude[i];
            // ∂A/∂t = -ħ/(2m) · (A · ∇²S + ∇S · ∇A)
            let lap_phi = (dphi_p - dphi_m) / dx;
            let da_dx = (a_p - a_m) / (2.0 * dx);
            let dphi_c = (self.phase[ip_idx] - self.phase[im_idx]) / (2.0 * dx);
            let da = -hbar / (2.0 * mass) * (a * lap_phi + dphi_c * da_dx) * dt;
            *new_amp_i = (a + da).max(0.0);
        }
        self.amplitude = new_amp;
    }

    /// Expectation value ⟨V⟩ = ∫|ψ|²·V dx.
    pub fn expectation_potential(&self) -> f64 {
        self.amplitude
            .iter()
            .zip(self.potential.iter())
            .map(|(a, v)| a * a * v * self.dx)
            .sum()
    }
}

// ---------------------------------------------------------------------------
// BoseEinsteinCondensate
// ---------------------------------------------------------------------------

/// BEC simulator with Thomas-Fermi ground state and dark soliton solutions.
///
/// The Thomas-Fermi (TF) approximation neglects kinetic energy, yielding the
/// ground-state density: ρ_TF = max(0, (μ - V) / g).
#[derive(Debug, Clone)]
pub struct BoseEinsteinCondensate {
    /// Number of lattice sites.
    pub nx: usize,
    /// Density profile ρ(x).
    pub density: Vec<f64>,
    /// Phase field θ(x).
    pub phase: Vec<f64>,
    /// Chemical potential μ.
    pub mu: f64,
    /// Interaction coupling g.
    pub coupling: f64,
    /// External trapping potential V(x).
    pub trap: Vec<f64>,
    /// Temperature.
    pub temperature: f64,
    /// Lattice spacing.
    pub dx: f64,
    /// Particle mass.
    pub mass: f64,
    /// Reduced Planck constant ħ.
    pub hbar: f64,
}

impl BoseEinsteinCondensate {
    /// Create a BEC with Thomas-Fermi ground state.
    ///
    /// # Arguments
    /// * `nx`         – number of sites
    /// * `mu`         – chemical potential
    /// * `coupling`   – contact interaction g
    /// * `trap`       – trapping potential V(x)
    /// * `temperature`– temperature
    /// * `dx`         – lattice spacing
    /// * `mass`       – particle mass
    /// * `hbar`       – reduced Planck constant
    pub fn new(
        nx: usize,
        mu: f64,
        coupling: f64,
        trap: Vec<f64>,
        temperature: f64,
        dx: f64,
        mass: f64,
        hbar: f64,
    ) -> Self {
        let density: Vec<f64> = trap
            .iter()
            .map(|&v| ((mu - v) / coupling).max(0.0))
            .collect();
        Self {
            nx,
            density,
            phase: vec![0.0; nx],
            mu,
            coupling,
            trap,
            temperature,
            dx,
            mass,
            hbar,
        }
    }

    /// Compute Thomas-Fermi density profile: ρ_TF(x) = max(0, (μ - V(x))/g).
    pub fn thomas_fermi_density(&self) -> Vec<f64> {
        self.trap
            .iter()
            .map(|&v| ((self.mu - v) / self.coupling).max(0.0))
            .collect()
    }

    /// Total number of particles N = ∫ρ dx.
    pub fn total_number(&self) -> f64 {
        self.density.iter().map(|&r| r * self.dx).sum()
    }

    /// Condensate fraction (simplified: approaches 1 at T=0).
    pub fn condensate_fraction(&self) -> f64 {
        superfluid_fraction(self.temperature, 1.0) // normalized to T_c=1
    }

    /// Place a dark soliton at position x_0 (sets a notch in the density).
    ///
    /// A dark soliton has ρ(x) ≈ ρ₀ · tanh²((x - x₀)/ξ).
    ///
    /// # Arguments
    /// * `x0`  – soliton center position (in lattice units)
    /// * `xi`  – healing length (lattice units)
    pub fn place_dark_soliton(&mut self, x0: f64, xi: f64) {
        let xi_safe = xi.max(1e-10);
        for i in 0..self.nx {
            let x = i as f64 * self.dx;
            let tanh_val = ((x - x0) / xi_safe).tanh();
            let tf = ((self.mu - self.trap[i]) / self.coupling).max(0.0);
            self.density[i] = tf * tanh_val * tanh_val;
            // dark soliton has a π phase jump at x₀
            self.phase[i] = if x < x0 { 0.0 } else { PI };
        }
    }

    /// Mean-field interaction energy E_int = (g/2) ∫ρ² dx.
    pub fn interaction_energy(&self) -> f64 {
        self.density
            .iter()
            .map(|&r| 0.5 * self.coupling * r * r * self.dx)
            .sum()
    }

    /// Healing length ξ = ħ / √(2m·g·n_max).
    pub fn healing_length(&self) -> f64 {
        let n_max = self.density.iter().cloned().fold(0.0_f64, f64::max);
        healing_length(self.hbar, self.mass, self.coupling * n_max)
    }
}

// ---------------------------------------------------------------------------
// QuantumPressure
// ---------------------------------------------------------------------------

/// Bohm quantum pressure (quantum potential) tensor.
///
/// The Bohm potential is Q = -ħ²∇²√ρ / (2m√ρ).
/// Its gradient gives the quantum force on the fluid.
#[derive(Debug, Clone)]
pub struct QuantumPressure {
    /// Density profile ρ(x).
    pub density: Vec<f64>,
    /// Reduced Planck constant ħ.
    pub hbar: f64,
    /// Particle mass m.
    pub mass: f64,
    /// Lattice spacing Δx.
    pub dx: f64,
}

impl QuantumPressure {
    /// Create a QuantumPressure calculator for a given density profile.
    ///
    /// # Arguments
    /// * `density` – density ρ(x) at each site
    /// * `hbar`    – reduced Planck constant
    /// * `mass`    – particle mass
    /// * `dx`      – lattice spacing
    pub fn new(density: Vec<f64>, hbar: f64, mass: f64, dx: f64) -> Self {
        Self {
            density,
            hbar,
            mass,
            dx,
        }
    }

    /// Compute Q(x) = -ħ²∇²√ρ / (2m√ρ) at each site.
    pub fn quantum_potential(&self) -> Vec<f64> {
        let n = self.density.len();
        let dx = self.dx;
        let hbar = self.hbar;
        let mass = self.mass;
        (0..n)
            .map(|i| {
                let im = if i == 0 { n - 1 } else { i - 1 };
                let ip = if i == n - 1 { 0 } else { i + 1 };
                let sqrt_rho = self.density[i].max(0.0).sqrt();
                let sqrt_rho_m = self.density[im].max(0.0).sqrt();
                let sqrt_rho_p = self.density[ip].max(0.0).sqrt();
                if sqrt_rho < 1e-30 {
                    return 0.0;
                }
                let lap = (sqrt_rho_p - 2.0 * sqrt_rho + sqrt_rho_m) / (dx * dx);
                -hbar * hbar / (2.0 * mass) * lap / sqrt_rho
            })
            .collect()
    }

    /// Quantum force F_Q(x) = -∇Q at each site.
    pub fn quantum_force(&self) -> Vec<f64> {
        let q = self.quantum_potential();
        let n = q.len();
        let dx = self.dx;
        (0..n)
            .map(|i| {
                let im = if i == 0 { n - 1 } else { i - 1 };
                let ip = if i == n - 1 { 0 } else { i + 1 };
                -(q[ip] - q[im]) / (2.0 * dx)
            })
            .collect()
    }

    /// Quantum pressure vanishes identically for uniform density.
    ///
    /// Returns `true` if max|Q| < tolerance.
    pub fn is_uniform(&self, tol: f64) -> bool {
        self.quantum_potential().iter().all(|&q| q.abs() < tol)
    }

    /// Fisher information I_F = ∫(∇√ρ)² dx.
    pub fn fisher_information(&self) -> f64 {
        let n = self.density.len();
        let dx = self.dx;
        (0..n)
            .map(|i| {
                let ip = if i == n - 1 { 0 } else { i + 1 };
                let dsqrt = (self.density[ip].sqrt() - self.density[i].sqrt()) / dx;
                dsqrt * dsqrt * dx
            })
            .sum()
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // ---- QuantumLbmParams --------------------------------------------------

    #[test]
    fn test_params_new() {
        let p = QuantumLbmParams::new(1.0, 1.0, 0.5, 2.0);
        assert_eq!(p.hbar, 1.0);
        assert_eq!(p.coupling, 2.0);
    }

    #[test]
    fn test_params_kappa() {
        let p = QuantumLbmParams::new(1.0, 1.0, 0.5, 1.0);
        assert!((p.kappa() - 2.0 * PI).abs() < 1e-12);
    }

    #[test]
    fn test_params_sound_speed() {
        let p = QuantumLbmParams::new(1.0, 4.0, 0.5, 2.0);
        let cs = p.sound_speed(3.0);
        let expected = (2.0 * 3.0 / 4.0_f64).sqrt();
        assert!((cs - expected).abs() < 1e-12);
    }

    #[test]
    fn test_params_from_physics() {
        let p = QuantumLbmParams::from_physics(1.0, 1.0, 1.0, 1.0);
        assert!(p.healing_length > 0.0);
    }

    #[test]
    fn test_healing_length_positive() {
        let xi = healing_length(HBAR, HELIUM4_MASS, 1e-40);
        assert!(xi > 0.0);
    }

    #[test]
    fn test_healing_length_zero_coupling() {
        let xi = healing_length(1.0, 1.0, 0.0);
        assert_eq!(xi, 0.0);
    }

    #[test]
    fn test_coherence_length_positive() {
        let xi = coherence_length(1.0, 1.0, 0.5);
        assert!(xi > 0.0);
    }

    #[test]
    fn test_coherence_length_zero_mu() {
        let xi = coherence_length(1.0, 1.0, 0.0);
        assert_eq!(xi, 0.0);
    }

    #[test]
    fn test_vortex_core_radius() {
        let xi = 1.0;
        let rc = vortex_core_radius(xi);
        assert!((rc - xi / 2.0_f64.sqrt()).abs() < 1e-14);
    }

    // ---- GrossPitaevskiiLbm -----------------------------------------------
    // (the physics regression tests live in `gp_tests.rs`)

    #[test]
    fn test_gp_lbm_new() {
        let p = QuantumLbmParams::new(1.0, 1.0, 1.0, 0.1);
        let gp = GrossPitaevskiiLbm::new(16, p, 1.0, 0.1);
        assert_eq!(gp.nx, 16);
        assert_eq!(gp.amplitude.len(), 16);
        assert_eq!(gp.phase.len(), 16);
        assert_eq!(gp.potential.len(), 16);
    }

    #[test]
    fn test_gp_density_from_amplitude() {
        let p = QuantumLbmParams::new(1.0, 1.0, 1.0, 0.1);
        let mut gp = GrossPitaevskiiLbm::new(8, p, 1.0, 0.1);
        gp.amplitude = vec![2.0; 8];
        let rho = gp.density();
        assert!(rho.iter().all(|&r| (r - 4.0).abs() < 1e-12));
    }

    #[test]
    fn test_gp_total_number_conserved() {
        // The split-step scheme is unitary: N is conserved to rounding error.
        let p = QuantumLbmParams::new(1.0, 1.0, 1.0, 0.001);
        let mut gp = GrossPitaevskiiLbm::new(32, p, 1.0, 0.01);
        let n0 = gp.total_number();
        gp.step();
        let n1 = gp.total_number();
        assert!((n1 - n0).abs() / n0 < 1e-13, "{n0} -> {n1}");
    }

    #[test]
    fn test_gp_step_runs() {
        let p = QuantumLbmParams::new(1.0, 1.0, 1.0, 0.1);
        let mut gp = GrossPitaevskiiLbm::new(16, p, 1.0, 0.01);
        gp.step();
        assert!(gp.amplitude.iter().all(|a| a.is_finite()));
        assert!(gp.phase.iter().all(|a| a.is_finite()));
    }

    #[test]
    fn test_gp_superfluid_velocity_zero_for_uniform_phase() {
        let p = QuantumLbmParams::new(1.0, 1.0, 1.0, 0.1);
        let gp = GrossPitaevskiiLbm::new(16, p, 1.0, 0.01);
        let v = gp.superfluid_velocity();
        assert!(v.iter().all(|&vi| vi.abs() < 1e-12));
    }

    // ---- QuantumVortex ----------------------------------------------------

    #[test]
    fn test_vortex_kappa() {
        let v = QuantumVortex::new(0.0, 0.0, 1, HELIUM4_MASS);
        assert!((v.kappa() - QUANTUM_CIRCULATION_HE4).abs() < 1e-40);
    }

    #[test]
    fn test_vortex_circulation() {
        let v = QuantumVortex::new(0.0, 0.0, 2, HELIUM4_MASS);
        assert!((v.circulation() - 2.0 * QUANTUM_CIRCULATION_HE4).abs() < 1e-40);
    }

    #[test]
    fn test_vortex_velocity_at_origin_is_zero() {
        let v = QuantumVortex::new(0.0, 0.0, 1, HELIUM4_MASS);
        let vel = v.velocity_at(0.0, 0.0);
        assert_eq!(vel, [0.0, 0.0]);
    }

    #[test]
    fn test_vortex_velocity_1_over_r() {
        // |v| ∝ 1/r for a point vortex
        let v = QuantumVortex::new(0.0, 0.0, 1, 1.0);
        let r1 = 1.0f64;
        let r2 = 2.0f64;
        let vel1 = v.velocity_at(r1, 0.0);
        let vel2 = v.velocity_at(r2, 0.0);
        let spd1 = (vel1[0] * vel1[0] + vel1[1] * vel1[1]).sqrt();
        let spd2 = (vel2[0] * vel2[0] + vel2[1] * vel2[1]).sqrt();
        assert!((spd1 / spd2 - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_vortex_antivortex_cancel_on_axis() {
        // A vortex-antivortex pair on the y-axis: v at x=0 midpoint should be non-zero
        let v1 = QuantumVortex::new(0.0, 1.0, 1, 1.0);
        let v2 = QuantumVortex::new(0.0, -1.0, -1, 1.0);
        let vel1 = v1.velocity_at(0.0, 0.0);
        let vel2 = v2.velocity_at(0.0, 0.0);
        let total_vx = vel1[0] + vel2[0];
        // By symmetry vx contributions should add
        assert!(total_vx.abs() > 0.0 || vel1[1].abs() + vel2[1].abs() < 1e-10);
    }

    // ---- VortexDynamics ---------------------------------------------------

    #[test]
    fn test_vortex_dynamics_new() {
        let v = QuantumVortex::new(0.0, 0.0, 1, 1.0);
        let dyn_ = VortexDynamics::new(vec![v], 0.01, 0.01);
        assert_eq!(dyn_.vortices.len(), 1);
    }

    #[test]
    fn test_vortex_dynamics_enstrophy() {
        let v1 = QuantumVortex::new(0.0, 0.0, 1, 1.0);
        let v2 = QuantumVortex::new(1.0, 0.0, -1, 1.0);
        let dyn_ = VortexDynamics::new(vec![v1, v2], 0.01, 0.01);
        let kappa = PLANCK_H;
        let expected = 2.0 * kappa * kappa;
        assert!((dyn_.enstrophy() - expected).abs() < 1e-50);
    }

    #[test]
    fn test_vortex_dynamics_pair_count() {
        let v1 = QuantumVortex::new(0.0, 0.0, 1, 1.0);
        let v2 = QuantumVortex::new(1.0, 0.0, -1, 1.0);
        let dyn_ = VortexDynamics::new(vec![v1, v2], 0.01, 0.01);
        assert_eq!(dyn_.pair_count(), 1);
    }

    #[test]
    fn test_vortex_dynamics_step() {
        let v1 = QuantumVortex::new(-1.0, 0.0, 1, 1.0);
        let v2 = QuantumVortex::new(1.0, 0.0, -1, 1.0);
        let mut dyn_ = VortexDynamics::new(vec![v1, v2], 0.001, 0.01);
        dyn_.step();
        // Vortices should have moved
        assert!(dyn_.vortices[0].x.is_finite());
        assert!(dyn_.vortices[1].x.is_finite());
    }

    // ---- SuperfluidLbm ----------------------------------------------------

    #[test]
    fn test_superfluid_fraction_at_zero_temp() {
        let f = superfluid_fraction(0.0, 2.17);
        assert!((f - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_superfluid_fraction_at_tc() {
        let f = superfluid_fraction(2.17, 2.17);
        assert_eq!(f, 0.0);
    }

    #[test]
    fn test_superfluid_fraction_above_tc() {
        let f = superfluid_fraction(3.0, 2.17);
        assert_eq!(f, 0.0);
    }

    #[test]
    fn test_superfluid_lbm_new() {
        let sf = SuperfluidLbm::new(32, 1.0, 2.17, 0.1, 0.6, 1.0, 0.01);
        assert_eq!(sf.nx, 32);
        assert!(sf.superfluid_fraction() > 0.0);
    }

    #[test]
    fn test_superfluid_lbm_fractions_sum_to_one() {
        let sf = SuperfluidLbm::new(16, 1.0, 2.17, 0.1, 0.6, 1.0, 0.01);
        for i in 0..sf.nx {
            let total = sf.rho_s[i] + sf.rho_n[i];
            assert!((total - 1.0).abs() < 1e-12);
        }
    }

    #[test]
    fn test_superfluid_lbm_step() {
        let mut sf = SuperfluidLbm::new(16, 1.0, 2.17, 0.1, 0.6, 1.0, 0.01);
        sf.vel_s.fill(0.05);
        sf.vel_n.fill(-0.02);
        let slip_before = sf.vel_s[0] - sf.vel_n[0];
        sf.step();
        assert!(sf.vel_s.iter().chain(&sf.vel_n).all(|v| v.is_finite()));
        // Mutual friction drags the two components towards each other: the
        // counter-flow decays (it used to grow) but does not reverse.
        let slip_after = sf.vel_s[0] - sf.vel_n[0];
        assert!(
            slip_after > 0.0 && slip_after < slip_before,
            "slip {slip_before} -> {slip_after}"
        );
    }

    // ---- QuantumTurbulence ------------------------------------------------

    #[test]
    fn test_qt_new() {
        let qt = QuantumTurbulence::new(1e6, QUANTUM_CIRCULATION_HE4, 1.5, 2.17, 0.1);
        assert!(qt.vortex_line_density > 0.0);
    }

    #[test]
    fn test_qt_inter_vortex_spacing() {
        let qt = QuantumTurbulence::new(1e4, 1.0, 1.0, 2.17, 0.1);
        let delta = qt.inter_vortex_spacing();
        assert!((delta - 0.01).abs() < 1e-10); // 1/sqrt(1e4) = 0.01
    }

    #[test]
    fn test_qt_vinen_step() {
        let mut qt = QuantumTurbulence::new(1e6, 1e-7, 1.5, 2.17, 0.1);
        qt.vinen_step(0.2, 1e-3);
        assert!(qt.vortex_line_density.is_finite());
        assert!(qt.vortex_line_density >= 0.0);
    }

    #[test]
    fn test_qt_kolmogorov_spectrum_decreasing() {
        let mut qt = QuantumTurbulence::new(1e8, QUANTUM_CIRCULATION_HE4, 1.5, 2.17, 0.1);
        qt.dissipation_rate = 1e-6;
        let e1 = qt.kolmogorov_spectrum(1.0);
        let e2 = qt.kolmogorov_spectrum(10.0);
        assert!(e1 > e2);
    }

    #[test]
    fn test_qt_reconnection_count() {
        let mut qt = QuantumTurbulence::new(1e6, 1.0, 1.5, 2.17, 0.1);
        qt.register_reconnection();
        qt.register_reconnection();
        assert_eq!(qt.reconnection_count, 2);
    }

    // ---- SchroedingerLbm --------------------------------------------------

    #[test]
    fn test_schrodinger_lbm_new() {
        let slbm = SchroedingerLbm::new(32, 1.0, 1.0, 1.0, 0.01);
        assert_eq!(slbm.nx, 32);
    }

    #[test]
    fn test_schrodinger_lbm_norm_positive() {
        let slbm = SchroedingerLbm::new(32, 1.0, 1.0, 1.0, 0.01);
        assert!(slbm.norm() > 0.0);
    }

    #[test]
    fn test_schrodinger_lbm_step_runs() {
        let mut slbm = SchroedingerLbm::new(32, 1.0, 1.0, 1.0, 0.001);
        slbm.step();
        assert!(slbm.amplitude.iter().all(|a| a.is_finite()));
    }

    #[test]
    fn test_schrodinger_lbm_probability_density_nonneg() {
        let slbm = SchroedingerLbm::new(16, 1.0, 1.0, 1.0, 0.01);
        assert!(slbm.probability_density().iter().all(|&r| r >= 0.0));
    }

    #[test]
    fn test_schrodinger_lbm_expectation_flat_potential() {
        let mut slbm = SchroedingerLbm::new(16, 1.0, 1.0, 1.0, 0.01);
        // Flat potential = 5.0, expectation = 5 * norm
        slbm.potential = vec![5.0; 16];
        let ev = slbm.expectation_potential();
        let expected = 5.0 * slbm.norm();
        assert!((ev - expected).abs() < 1e-10);
    }

    // ---- BoseEinsteinCondensate -------------------------------------------

    #[test]
    fn test_bec_new() {
        let trap = vec![0.0; 32];
        let bec = BoseEinsteinCondensate::new(32, 1.0, 1.0, trap, 0.0, 1.0, 1.0, 1.0);
        assert_eq!(bec.nx, 32);
    }

    #[test]
    fn test_bec_thomas_fermi_density_shape() {
        let nx = 64;
        let mu = 2.0;
        let g = 1.0;
        let trap: Vec<f64> = (0..nx)
            .map(|i| {
                let x = (i as f64 - 32.0) * 0.1;
                0.5 * x * x
            })
            .collect();
        let bec = BoseEinsteinCondensate::new(nx, mu, g, trap, 0.0, 0.1, 1.0, 1.0);
        let tf = bec.thomas_fermi_density();
        // Density at center should be highest
        assert!(tf[32] >= tf[0]);
        assert!(tf[32] >= tf[63]);
    }

    #[test]
    fn test_bec_total_number_positive() {
        let trap = vec![0.0; 32];
        let bec = BoseEinsteinCondensate::new(32, 1.0, 0.5, trap, 0.0, 1.0, 1.0, 1.0);
        assert!(bec.total_number() > 0.0);
    }

    #[test]
    fn test_bec_healing_length_positive() {
        let trap = vec![0.0; 32];
        let bec = BoseEinsteinCondensate::new(32, 1.0, 0.5, trap, 0.0, 1.0, 1.0, 1.0);
        let xi = bec.healing_length();
        assert!(xi > 0.0);
    }

    #[test]
    fn test_bec_dark_soliton_notch() {
        let trap = vec![0.0; 64];
        let mut bec = BoseEinsteinCondensate::new(64, 1.0, 0.5, trap, 0.0, 1.0, 1.0, 1.0);
        bec.place_dark_soliton(32.0, 2.0);
        // Density at soliton center should be near zero
        let center_rho = bec.density[32];
        let far_rho = bec.density[0].max(bec.density[63]);
        assert!(center_rho < far_rho);
    }

    #[test]
    fn test_bec_interaction_energy_positive() {
        let trap = vec![0.0; 32];
        let bec = BoseEinsteinCondensate::new(32, 1.0, 0.5, trap, 0.0, 1.0, 1.0, 1.0);
        assert!(bec.interaction_energy() > 0.0);
    }

    // ---- QuantumPressure --------------------------------------------------

    #[test]
    fn test_quantum_pressure_new() {
        let rho = vec![1.0; 16];
        let qp = QuantumPressure::new(rho, 1.0, 1.0, 0.1);
        assert_eq!(qp.density.len(), 16);
    }

    #[test]
    fn test_quantum_pressure_uniform_is_zero() {
        let rho = vec![2.0; 32];
        let qp = QuantumPressure::new(rho, 1.0, 1.0, 0.1);
        assert!(qp.is_uniform(1e-10));
    }

    #[test]
    fn test_quantum_pressure_nonuniform_nonzero() {
        let mut rho = vec![1.0; 32];
        rho[16] = 4.0;
        let qp = QuantumPressure::new(rho, 1.0, 1.0, 0.1);
        let q = qp.quantum_potential();
        assert!(q.iter().any(|&qi| qi.abs() > 1e-10));
    }

    #[test]
    fn test_quantum_force_finite() {
        let rho: Vec<f64> = (0..32).map(|i| 1.0 + 0.1 * (i as f64).sin()).collect();
        let qp = QuantumPressure::new(rho, 1.0, 1.0, 0.1);
        let f = qp.quantum_force();
        assert!(f.iter().all(|fi| fi.is_finite()));
    }

    #[test]
    fn test_fisher_information_positive() {
        let rho: Vec<f64> = (0..32)
            .map(|i| 1.0 + 0.1 * (i as f64 * 0.2).sin())
            .collect();
        let qp = QuantumPressure::new(rho, 1.0, 1.0, 1.0);
        assert!(qp.fisher_information() > 0.0);
    }

    #[test]
    fn test_fisher_information_zero_for_uniform() {
        let rho = vec![1.0; 32];
        let qp = QuantumPressure::new(rho, 1.0, 1.0, 1.0);
        assert!(qp.fisher_information() < 1e-14);
    }
}

#[cfg(test)]
mod gp_tests;
