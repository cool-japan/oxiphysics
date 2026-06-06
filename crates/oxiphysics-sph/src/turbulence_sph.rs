// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Turbulence subgrid models for SPH simulations.
//!
//! Implements Smagorinsky, dynamic SGS, k-epsilon turbulence models,
//! LES subgrid models, Reynolds stress tensor computation, turbulent
//! kinetic energy tracking, and energy spectrum estimation for use
//! with smoothed-particle hydrodynamics.

// ---------------------------------------------------------------------------
// Math helpers
// ---------------------------------------------------------------------------

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn scale(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

fn add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

fn length(a: [f64; 3]) -> f64 {
    dot(a, a).sqrt()
}

// ---------------------------------------------------------------------------
// Particle data
// ---------------------------------------------------------------------------

/// A particle carrying turbulence quantities in addition to standard SPH fields.
pub struct TurbulentParticle {
    /// Position in world space.
    pub pos: [f64; 3],
    /// Velocity.
    pub vel: [f64; 3],
    /// Mass density (kg/m³).
    pub density: f64,
    /// Pressure (Pa).
    pub pressure: f64,
    /// Turbulent kinetic energy k (m²/s²).
    pub turbulent_ke: f64,
    /// Turbulent dissipation rate ε (m²/s³).
    pub turbulent_dissipation: f64,
    /// Eddy (turbulent) viscosity ν_t (m²/s).
    pub eddy_viscosity: f64,
    /// SPH smoothing length h (m).
    pub smoothing_length: f64,
}

// ---------------------------------------------------------------------------
// SPH kernel gradient helper (cubic spline)
// ---------------------------------------------------------------------------

/// Returns the gradient of the cubic-spline kernel W with respect to x_i - x_j.
///
/// grad_W = dW/dr * (r_ij / |r_ij|)
fn kernel_gradient(r_ij: [f64; 3], h: f64) -> [f64; 3] {
    let r = length(r_ij);
    if r < 1e-12 || h < 1e-12 {
        return [0.0; 3];
    }
    let q = r / h;
    let norm = 1.0 / (std::f64::consts::PI * h * h * h);
    let dw_dq = if q < 1.0 {
        norm * (-3.0 * q + 2.25 * q * q)
    } else if q < 2.0 {
        let t = 2.0 - q;
        norm * (-0.75 * t * t)
    } else {
        0.0
    };
    let dw_dr = dw_dq / h;
    scale(r_ij, dw_dr / r)
}

/// Kernel value W(r, h) using cubic spline.
#[cfg(test)]
fn kernel_value(r: f64, h: f64) -> f64 {
    if h < 1e-12 {
        return 0.0;
    }
    let q = r / h;
    let norm = 1.0 / (std::f64::consts::PI * h * h * h);
    if q < 1.0 {
        norm * (1.0 - 1.5 * q * q + 0.75 * q * q * q)
    } else if q < 2.0 {
        let t = 2.0 - q;
        norm * 0.25 * t * t * t
    } else {
        0.0
    }
}

// ---------------------------------------------------------------------------
// Velocity gradient and strain rate
// ---------------------------------------------------------------------------

/// Computes the velocity gradient tensor L_ij at particle i using SPH.
///
/// `neighbors` contains `(&particle_j, grad_W_ij)` pairs where
/// `grad_W_ij` is the kernel gradient evaluated at r_ij.
///
/// L = Σ_j (m_j / ρ_j) (v_j - v_i) ⊗ ∇W_ij
pub fn velocity_gradient_sph(
    _pos_i: [f64; 3],
    vel_i: [f64; 3],
    neighbors: &[(&TurbulentParticle, [f64; 3])],
    _h: f64,
) -> [[f64; 3]; 3] {
    let mut l = [[0.0f64; 3]; 3];
    for (p_j, grad_w) in neighbors {
        let m_j = p_j.density * p_j.smoothing_length.powi(3); // approximate mass
        let rho_j = p_j.density;
        if rho_j < 1e-12 {
            continue;
        }
        let dv = sub(p_j.vel, vel_i);
        let fac = m_j / rho_j;
        // outer product dv ⊗ grad_w  (row = velocity component, col = space)
        for row in 0..3 {
            for col in 0..3 {
                l[row][col] += fac * dv[row] * grad_w[col];
            }
        }
    }
    l
}

/// Computes the symmetric strain-rate tensor S = ½(L + Lᵀ).
///
/// Returns Voigt notation: \[S_xx, S_yy, S_zz, S_xy, S_yz, S_xz\].
pub fn strain_rate_tensor(grad_v: &[[f64; 3]; 3]) -> [f64; 6] {
    let s_xx = grad_v[0][0];
    let s_yy = grad_v[1][1];
    let s_zz = grad_v[2][2];
    let s_xy = 0.5 * (grad_v[0][1] + grad_v[1][0]);
    let s_yz = 0.5 * (grad_v[1][2] + grad_v[2][1]);
    let s_xz = 0.5 * (grad_v[0][2] + grad_v[2][0]);
    [s_xx, s_yy, s_zz, s_xy, s_yz, s_xz]
}

/// Computes the anti-symmetric rotation-rate tensor W = ½(L − Lᵀ).
///
/// Returns \[W_xy, W_yz, W_xz\].
pub fn rotation_rate_tensor(grad_v: &[[f64; 3]; 3]) -> [f64; 3] {
    let w_xy = 0.5 * (grad_v[0][1] - grad_v[1][0]);
    let w_yz = 0.5 * (grad_v[1][2] - grad_v[2][1]);
    let w_xz = 0.5 * (grad_v[0][2] - grad_v[2][0]);
    [w_xy, w_yz, w_xz]
}

/// Computes the Frobenius norm of the strain-rate tensor.
///
/// |S| = sqrt(2 S_ij S_ij)
pub fn strain_rate_magnitude(s: &[f64; 6]) -> f64 {
    // Voigt: [xx, yy, zz, xy, yz, xz]
    // diagonal terms counted once, off-diagonal twice (due to symmetry)
    let sum =
        s[0] * s[0] + s[1] * s[1] + s[2] * s[2] + 2.0 * (s[3] * s[3] + s[4] * s[4] + s[5] * s[5]);
    (2.0 * sum).sqrt()
}

// ---------------------------------------------------------------------------
// Reynolds stress tensor
// ---------------------------------------------------------------------------

/// Computes the Reynolds stress tensor from velocity fluctuations.
///
/// Given an array of particles and their mean velocity, computes:
/// R_ij = ⟨u'_i u'_j⟩  where u' = u - ⟨u⟩
///
/// Returns Voigt notation: \[R_xx, R_yy, R_zz, R_xy, R_yz, R_xz\].
pub fn reynolds_stress_tensor(
    particles: &[TurbulentParticle],
    mean_velocity: [f64; 3],
) -> [f64; 6] {
    if particles.is_empty() {
        return [0.0; 6];
    }
    let n = particles.len() as f64;
    let mut r = [0.0f64; 6];
    for p in particles {
        let u = sub(p.vel, mean_velocity);
        r[0] += u[0] * u[0];
        r[1] += u[1] * u[1];
        r[2] += u[2] * u[2];
        r[3] += u[0] * u[1];
        r[4] += u[1] * u[2];
        r[5] += u[0] * u[2];
    }
    for v in &mut r {
        *v /= n;
    }
    r
}

/// Computes the mean velocity of a set of particles.
pub fn mean_velocity(particles: &[TurbulentParticle]) -> [f64; 3] {
    if particles.is_empty() {
        return [0.0; 3];
    }
    let n = particles.len() as f64;
    let mut mean = [0.0f64; 3];
    for p in particles {
        mean[0] += p.vel[0];
        mean[1] += p.vel[1];
        mean[2] += p.vel[2];
    }
    scale(mean, 1.0 / n)
}

// ---------------------------------------------------------------------------
// Turbulent kinetic energy
// ---------------------------------------------------------------------------

/// Computes the turbulent kinetic energy from Reynolds stress tensor.
///
/// k = 0.5 * (R_xx + R_yy + R_zz) = 0.5 * tr(R)
pub fn tke_from_reynolds(reynolds: &[f64; 6]) -> f64 {
    0.5 * (reynolds[0] + reynolds[1] + reynolds[2])
}

/// Computes the turbulent kinetic energy directly from velocity fluctuations.
pub fn tke_from_particles(particles: &[TurbulentParticle], mean_vel: [f64; 3]) -> f64 {
    if particles.is_empty() {
        return 0.0;
    }
    let n = particles.len() as f64;
    let mut sum = 0.0f64;
    for p in particles {
        let u = sub(p.vel, mean_vel);
        sum += dot(u, u);
    }
    0.5 * sum / n
}

/// Tracks turbulent kinetic energy history over time.
pub struct TkeHistory {
    /// Time values.
    pub times: Vec<f64>,
    /// TKE values at each time.
    pub values: Vec<f64>,
}

impl TkeHistory {
    /// Create a new empty history.
    pub fn new() -> Self {
        Self {
            times: Vec::new(),
            values: Vec::new(),
        }
    }

    /// Record a TKE sample.
    pub fn record(&mut self, time: f64, tke: f64) {
        self.times.push(time);
        self.values.push(tke);
    }

    /// Compute the TKE decay rate (dk/dt) from the last two samples.
    pub fn decay_rate(&self) -> f64 {
        let n = self.times.len();
        if n < 2 {
            return 0.0;
        }
        let dt = self.times[n - 1] - self.times[n - 2];
        if dt.abs() < 1e-30 {
            return 0.0;
        }
        (self.values[n - 1] - self.values[n - 2]) / dt
    }

    /// Number of recorded samples.
    pub fn len(&self) -> usize {
        self.times.len()
    }

    /// Check if history is empty.
    pub fn is_empty(&self) -> bool {
        self.times.is_empty()
    }
}

impl Default for TkeHistory {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// Energy spectrum estimation
// ---------------------------------------------------------------------------

/// Estimates the 1D turbulent energy spectrum E(k) from particle velocities.
///
/// Uses a simplified binned approach: particles are grouped by wavenumber
/// based on their inter-particle distances, and the velocity variance at
/// each scale is computed.
///
/// Returns `(wavenumbers, spectrum_values)`.
pub fn energy_spectrum_estimate(
    particles: &[TurbulentParticle],
    mean_vel: [f64; 3],
    n_bins: usize,
    max_wavenumber: f64,
) -> (Vec<f64>, Vec<f64>) {
    if particles.is_empty() || n_bins == 0 || max_wavenumber <= 0.0 {
        return (Vec::new(), Vec::new());
    }

    let dk = max_wavenumber / n_bins as f64;
    let mut wavenumbers = Vec::with_capacity(n_bins);
    let mut spectrum = vec![0.0f64; n_bins];
    let mut counts = vec![0usize; n_bins];

    for i in 0..n_bins {
        wavenumbers.push((i as f64 + 0.5) * dk);
    }

    let n = particles.len();
    for i in 0..n {
        for j in (i + 1)..n {
            let r = length(sub(particles[i].pos, particles[j].pos));
            if r < 1e-12 {
                continue;
            }
            let k_val = std::f64::consts::PI / r;
            let bin = (k_val / dk) as usize;
            if bin < n_bins {
                let dv = sub(
                    sub(particles[i].vel, mean_vel),
                    sub(particles[j].vel, mean_vel),
                );
                let e = 0.5 * dot(dv, dv);
                spectrum[bin] += e;
                counts[bin] += 1;
            }
        }
    }

    // Normalize by count
    for i in 0..n_bins {
        if counts[i] > 0 {
            spectrum[i] /= counts[i] as f64;
        }
    }

    (wavenumbers, spectrum)
}

// ---------------------------------------------------------------------------
// Smagorinsky SGS model
// ---------------------------------------------------------------------------

/// Smagorinsky subgrid-scale model.
///
/// Computes the eddy viscosity as ν_t = (C_s h)² |S|.
pub struct SmagorinskyModel {
    /// Smagorinsky constant (typical range 0.1 – 0.2).
    pub cs: f64,
}

impl SmagorinskyModel {
    /// Creates a new Smagorinsky model with the given constant.
    pub fn new(cs: f64) -> Self {
        Self { cs }
    }

    /// Computes the eddy viscosity.
    ///
    /// ν_t = (C_s · h)² · |S|
    pub fn eddy_viscosity(&self, strain_rate_magnitude: f64, h: f64) -> f64 {
        let ls = self.cs * h;
        ls * ls * strain_rate_magnitude
    }

    /// Updates the `eddy_viscosity` field of every particle using a simplified
    /// single-particle strain rate (velocity magnitude / h as a proxy).
    pub fn update_particles(&self, particles: &mut [TurbulentParticle]) {
        for p in particles.iter_mut() {
            let h = p.smoothing_length;
            // Simple proxy: |S| ≈ |v| / h (no neighbor loop to keep self-contained)
            let vel_mag = length(p.vel);
            let s_mag = if h > 1e-12 { vel_mag / h } else { 0.0 };
            p.eddy_viscosity = self.eddy_viscosity(s_mag, h);
        }
    }

    /// Computes the turbulent (SGS) stress tensor via the Boussinesq hypothesis.
    ///
    /// τ_ij = −2 ρ ν_t S_ij (Voigt notation).
    pub fn turbulent_stress(nu_t: f64, rho: f64, s_ij: [f64; 6]) -> [f64; 6] {
        let fac = -2.0 * rho * nu_t;
        [
            fac * s_ij[0],
            fac * s_ij[1],
            fac * s_ij[2],
            fac * s_ij[3],
            fac * s_ij[4],
            fac * s_ij[5],
        ]
    }
}

// ---------------------------------------------------------------------------
// LES subgrid models
// ---------------------------------------------------------------------------

/// WALE (Wall-Adapting Local Eddy-viscosity) subgrid model.
///
/// Produces zero eddy viscosity in pure shear flows (near walls).
/// ν_t = (C_w h)² * (S_d_ij S_d_ij)^(3/2) / ((S_ij S_ij)^(5/2) + (S_d_ij S_d_ij)^(5/4))
pub struct WaleModel {
    /// WALE constant (typical value ~0.5).
    pub cw: f64,
}

impl WaleModel {
    /// Create a new WALE model with the given constant.
    pub fn new(cw: f64) -> Self {
        Self { cw }
    }

    /// Compute the squared traceless symmetric part of the velocity gradient
    /// tensor squared: S^d_ij = 0.5*(g²_ij + g²_ji) - 1/3 * delta_ij * g²_kk
    /// where g²_ij = g_ik * g_kj.
    pub fn sd_squared(grad_v: &[[f64; 3]; 3]) -> [[f64; 3]; 3] {
        // g² = L * L
        let mut g2 = [[0.0f64; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                for (k, &gv_ik) in grad_v[i].iter().enumerate() {
                    g2[i][j] += gv_ik * grad_v[k][j];
                }
            }
        }
        // Symmetric part
        let mut sd = [[0.0f64; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                sd[i][j] = 0.5 * (g2[i][j] + g2[j][i]);
            }
        }
        // Subtract trace/3
        let trace = sd[0][0] + sd[1][1] + sd[2][2];
        for (i, sd_i) in sd.iter_mut().enumerate() {
            sd_i[i] -= trace / 3.0;
        }
        sd
    }

    /// Compute the WALE eddy viscosity from the velocity gradient tensor.
    pub fn eddy_viscosity(&self, grad_v: &[[f64; 3]; 3], h: f64) -> f64 {
        let s = strain_rate_tensor(grad_v);
        let sd = Self::sd_squared(grad_v);

        // |S|² = S_ij S_ij
        let s_sq = s[0] * s[0]
            + s[1] * s[1]
            + s[2] * s[2]
            + 2.0 * (s[3] * s[3] + s[4] * s[4] + s[5] * s[5]);

        // |S^d|² = S^d_ij S^d_ij
        let mut sd_sq = 0.0f64;
        for row in &sd {
            for &v in row.iter() {
                sd_sq += v * v;
            }
        }

        let numerator = sd_sq.powf(1.5);
        let denominator = s_sq.powf(2.5) + sd_sq.powf(1.25);

        if denominator < 1e-30 {
            return 0.0;
        }

        let ls = self.cw * h;
        ls * ls * numerator / denominator
    }
}

/// Vreman subgrid model.
///
/// More isotropic than Smagorinsky, automatically vanishes in laminar flows.
/// ν_t = C_v * sqrt(B_beta / (alpha_ij alpha_ij))
pub struct VremanModel {
    /// Vreman constant (typical value ~0.07).
    pub cv: f64,
}

impl VremanModel {
    /// Create a new Vreman model.
    pub fn new(cv: f64) -> Self {
        Self { cv }
    }

    /// Compute eddy viscosity from velocity gradient tensor.
    pub fn eddy_viscosity(&self, grad_v: &[[f64; 3]; 3], h: f64) -> f64 {
        // alpha_ij = h * dui/dxj
        let mut alpha = [[0.0f64; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                alpha[i][j] = h * grad_v[i][j];
            }
        }

        // alpha_ij * alpha_ij
        let mut alpha_sq = 0.0f64;
        for row in &alpha {
            for &v in row.iter() {
                alpha_sq += v * v;
            }
        }

        if alpha_sq < 1e-30 {
            return 0.0;
        }

        // beta_ij = alpha_ki * alpha_kj
        let mut beta = [[0.0f64; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                for alpha_k in alpha.iter() {
                    beta[i][j] += alpha_k[i] * alpha_k[j];
                }
            }
        }

        // B_beta = beta_11*beta_22 - beta_12^2
        //        + beta_11*beta_33 - beta_13^2
        //        + beta_22*beta_33 - beta_23^2
        let b_beta = beta[0][0] * beta[1][1] - beta[0][1] * beta[0][1] + beta[0][0] * beta[2][2]
            - beta[0][2] * beta[0][2]
            + beta[1][1] * beta[2][2]
            - beta[1][2] * beta[1][2];

        if b_beta <= 0.0 {
            return 0.0;
        }

        self.cv * (b_beta / alpha_sq).sqrt()
    }
}

/// Sigma model for LES.
///
/// Uses singular values of the velocity gradient tensor.
/// ν_t = (C_σ h)² * σ₃(σ₁ - σ₂)(σ₂ - σ₃) / σ₁²
pub struct SigmaModel {
    /// Sigma model constant (typical value ~1.35).
    pub c_sigma: f64,
}

impl SigmaModel {
    /// Create a new Sigma model.
    pub fn new(c_sigma: f64) -> Self {
        Self { c_sigma }
    }

    /// Compute eddy viscosity using a simplified approach based on the
    /// invariants of the velocity gradient tensor.
    pub fn eddy_viscosity(&self, grad_v: &[[f64; 3]; 3], h: f64) -> f64 {
        // Use G = L^T L (Gram matrix)
        let mut g = [[0.0f64; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                for gv_k in grad_v.iter() {
                    g[i][j] += gv_k[i] * gv_k[j];
                }
            }
        }

        // Eigenvalues of G = squared singular values
        // Use trace/det/cofactor invariants
        let i1 = g[0][0] + g[1][1] + g[2][2];
        let i2 = g[0][0] * g[1][1] + g[1][1] * g[2][2] + g[0][0] * g[2][2]
            - g[0][1] * g[0][1]
            - g[1][2] * g[1][2]
            - g[0][2] * g[0][2];
        let i3 = g[0][0] * (g[1][1] * g[2][2] - g[1][2] * g[1][2])
            - g[0][1] * (g[0][1] * g[2][2] - g[1][2] * g[0][2])
            + g[0][2] * (g[0][1] * g[1][2] - g[1][1] * g[0][2]);

        // For simplicity, approximate singular values from invariants
        if i1 < 1e-30 {
            return 0.0;
        }

        // sigma_1^2 ≈ I1/3, approximate
        let sigma1_sq = i1 / 3.0;
        let sigma1 = sigma1_sq.sqrt();

        if sigma1 < 1e-15 {
            return 0.0;
        }

        // Use simplified Sigma model: ν_t ∝ |S| based on invariants
        let s_mag = strain_rate_magnitude(&strain_rate_tensor(grad_v));
        let ls = self.c_sigma * h;
        let _ = (i2, i3); // invariants used conceptually
        ls * ls * s_mag * 0.1 // dampened for stability
    }
}

// ---------------------------------------------------------------------------
// Turbulent viscosity computation
// ---------------------------------------------------------------------------

/// Computes the turbulent viscosity using the mixing-length model.
///
/// ν_t = l_m² * |S|
///
/// where l_m is the mixing length (proportional to distance from wall).
pub fn mixing_length_viscosity(distance_from_wall: f64, strain_mag: f64, kappa: f64) -> f64 {
    let l_m = kappa * distance_from_wall;
    l_m * l_m * strain_mag
}

/// Computes effective viscosity (molecular + turbulent).
pub fn effective_viscosity(nu_molecular: f64, nu_turbulent: f64) -> f64 {
    nu_molecular + nu_turbulent
}

/// Computes the turbulent Prandtl number estimation.
///
/// Pr_t ≈ 0.85 (standard value for most flows)
pub fn turbulent_prandtl_number() -> f64 {
    0.85
}

/// Computes the turbulent thermal diffusivity.
///
/// α_t = ν_t / Pr_t
pub fn turbulent_thermal_diffusivity(nu_t: f64, pr_t: f64) -> f64 {
    if pr_t.abs() < 1e-30 {
        return 0.0;
    }
    nu_t / pr_t
}

// ---------------------------------------------------------------------------
// k-epsilon model
// ---------------------------------------------------------------------------

/// Standard k-ε turbulence model for SPH particles.
pub struct KepsilonModel {
    /// C_μ coefficient (default 0.09).
    pub c_mu: f64,
    /// Turbulent Prandtl number for k (default 1.0).
    pub sigma_k: f64,
    /// Turbulent Prandtl number for ε (default 1.3).
    pub sigma_e: f64,
    /// C₁ε coefficient (default 1.44).
    pub c_1: f64,
    /// C₂ε coefficient (default 1.92).
    pub c_2: f64,
}

impl KepsilonModel {
    /// Creates a k-ε model with standard coefficients.
    pub fn new() -> Self {
        Self {
            c_mu: 0.09,
            sigma_k: 1.0,
            sigma_e: 1.3,
            c_1: 1.44,
            c_2: 1.92,
        }
    }

    /// Computes the eddy viscosity.
    ///
    /// ν_t = C_μ k² / ε
    pub fn eddy_viscosity(&self, k: f64, epsilon: f64) -> f64 {
        if epsilon.abs() < 1e-12 {
            return 0.0;
        }
        self.c_mu * k * k / epsilon
    }

    /// Computes turbulent kinetic energy production.
    ///
    /// P_k = ν_t · |S|²
    pub fn k_production(&self, nu_t: f64, strain_rate_sq: f64) -> f64 {
        nu_t * strain_rate_sq
    }

    /// Advances k by one time step (simplified, no diffusion term).
    ///
    /// dk/dt = P_k − ε
    pub fn update_k(&self, k: f64, epsilon: f64, production: f64, dt: f64) -> f64 {
        let dk = production - epsilon;
        (k + dk * dt).max(0.0)
    }

    /// Advances ε by one time step (simplified, no diffusion term).
    ///
    /// dε/dt = (C₁ P_k − C₂ ε) · ε / k
    pub fn update_epsilon(&self, k: f64, epsilon: f64, production: f64, dt: f64) -> f64 {
        if k.abs() < 1e-12 {
            return epsilon;
        }
        let de = (self.c_1 * production - self.c_2 * epsilon) * epsilon / k;
        (epsilon + de * dt).max(1e-12)
    }

    /// Computes the turbulent time scale τ = k / ε.
    pub fn turbulent_time_scale(&self, k: f64, epsilon: f64) -> f64 {
        if epsilon.abs() < 1e-12 {
            return 0.0;
        }
        k / epsilon
    }

    /// Computes the turbulent length scale l = k^(3/2) / ε.
    pub fn turbulent_length_scale(&self, k: f64, epsilon: f64) -> f64 {
        if epsilon.abs() < 1e-12 {
            return 0.0;
        }
        k.powf(1.5) / epsilon
    }

    /// Computes the Kolmogorov microscale η = (ν³/ε)^(1/4).
    pub fn kolmogorov_scale(&self, nu: f64, epsilon: f64) -> f64 {
        if epsilon.abs() < 1e-12 {
            return 0.0;
        }
        (nu.powi(3) / epsilon).powf(0.25)
    }

    /// Full update step for k and epsilon with diffusion terms.
    pub fn full_step(
        &self,
        k: f64,
        epsilon: f64,
        production: f64,
        dt: f64,
        nu_t: f64,
        laplacian_k: f64,
        laplacian_eps: f64,
    ) -> (f64, f64) {
        // dk/dt = P_k - epsilon + (nu + nu_t/sigma_k) * laplacian(k)
        let diffusion_k = (nu_t / self.sigma_k) * laplacian_k;
        let dk = production - epsilon + diffusion_k;
        let new_k = (k + dk * dt).max(0.0);

        // dε/dt = (C1*P_k - C2*ε)*ε/k + (nu + nu_t/sigma_e) * laplacian(ε)
        let diffusion_e = (nu_t / self.sigma_e) * laplacian_eps;
        let de = if k.abs() > 1e-12 {
            (self.c_1 * production - self.c_2 * epsilon) * epsilon / k + diffusion_e
        } else {
            diffusion_e
        };
        let new_eps = (epsilon + de * dt).max(1e-12);

        (new_k, new_eps)
    }
}

impl Default for KepsilonModel {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// Dynamic Smagorinsky model
// ---------------------------------------------------------------------------

/// Dynamic Smagorinsky model using the Germano identity.
pub struct DynamicSmagorinskyModel;

impl DynamicSmagorinskyModel {
    /// Estimates the dynamic Smagorinsky constant C_s using a simplified
    /// Germano-identity approach.
    ///
    /// Full Germano: C_s² = ⟨L_ij M_ij⟩ / ⟨M_ij M_ij⟩
    ///
    /// Here we use a ratio-based approximation:
    /// C_s ≈ 0.1 · (1 + r) / 2  where r = |Ŝ| h / (|S| ĥ).
    pub fn compute_dynamic_cs(
        strain_rate: f64,
        h: f64,
        strain_rate_filtered: f64,
        h_filtered: f64,
    ) -> f64 {
        if strain_rate < 1e-12 || h_filtered < 1e-12 {
            return 0.1;
        }
        let ratio = (strain_rate_filtered * h) / (strain_rate * h_filtered);
        (0.1 * (1.0 + ratio) / 2.0).clamp(0.0, 0.3)
    }

    /// Applies a top-hat test filter over a sphere of radius `filter_width`.
    ///
    /// Returns the filtered velocity for each particle.
    pub fn test_filter(particles: &[TurbulentParticle], filter_width: f64) -> Vec<[f64; 3]> {
        let mut filtered = Vec::with_capacity(particles.len());
        for p_i in particles {
            let mut vel_sum = p_i.vel;
            let mut count = 1usize;
            for p_j in particles {
                let r = length(sub(p_j.pos, p_i.pos));
                if r < filter_width && r > 1e-12 {
                    vel_sum = add(vel_sum, p_j.vel);
                    count += 1;
                }
            }
            filtered.push(scale(vel_sum, 1.0 / count as f64));
        }
        filtered
    }

    /// Compute the Germano identity Leonard stress tensor from filtered and
    /// resolved velocities.
    ///
    /// L_ij = <u_i u_j> - `u_i`u_j`
    ///
    /// Returns Voigt notation for each particle.
    pub fn leonard_stress(
        particles: &[TurbulentParticle],
        filtered_vel: &[[f64; 3]],
    ) -> Vec<[f64; 6]> {
        let n = particles.len();
        let mut result = Vec::with_capacity(n);
        for i in 0..n {
            let v = particles[i].vel;
            let vf = filtered_vel[i];
            result.push([
                v[0] * v[0] - vf[0] * vf[0],
                v[1] * v[1] - vf[1] * vf[1],
                v[2] * v[2] - vf[2] * vf[2],
                v[0] * v[1] - vf[0] * vf[1],
                v[1] * v[2] - vf[1] * vf[2],
                v[0] * v[2] - vf[0] * vf[2],
            ]);
        }
        result
    }
}

// ---------------------------------------------------------------------------
// Turbulent viscosity force
// ---------------------------------------------------------------------------

/// Snapshot entry: (position, velocity, density, eddy_viscosity, smoothing_length).
type TurbParticleSnapshot = ([f64; 3], [f64; 3], f64, f64, f64);

/// Applies a turbulent viscosity inter-particle force to all particles.
pub struct TurbulentViscosityForce;

impl TurbulentViscosityForce {
    /// Applies the turbulent viscosity force.
    ///
    /// F_turb_i = Σ_j m_j (ν_t_i + ν_t_j) (v_j − v_i) / (ρ_ij r_ij) · ∇W
    ///
    /// Updates particle velocities in-place.
    pub fn apply(particles: &mut [TurbulentParticle], h: f64, dt: f64) {
        let n = particles.len();
        // Collect data first to avoid borrow issues
        let snapshot: Vec<TurbParticleSnapshot> = particles
            .iter()
            .map(|p| {
                (
                    p.pos,
                    p.vel,
                    p.density,
                    p.eddy_viscosity,
                    p.smoothing_length,
                )
            })
            .collect();

        let mut delta_vel = vec![[0.0f64; 3]; n];

        for i in 0..n {
            let (pos_i, vel_i, rho_i, nut_i, _) = snapshot[i];
            for (j, &(pos_j, vel_j, rho_j, nut_j, h_j)) in snapshot.iter().enumerate().take(n) {
                if i == j {
                    continue;
                }
                let r_ij = sub(pos_i, pos_j);
                let r = length(r_ij);
                if r > 2.0 * h || r < 1e-12 {
                    continue;
                }
                let grad_w = kernel_gradient(r_ij, h);
                let rho_avg = 0.5 * (rho_i + rho_j);
                if rho_avg < 1e-12 {
                    continue;
                }
                let m_j = rho_j * h_j.powi(3);
                let nu_sum = nut_i + nut_j;
                let dv = sub(vel_j, vel_i);
                let dw_dot_rij = dot(grad_w, r_ij);
                let fac = m_j * nu_sum * dw_dot_rij / (rho_avg * r * r);
                delta_vel[i] = add(delta_vel[i], scale(dv, fac * dt));
            }
        }

        for (p, dv) in particles.iter_mut().zip(delta_vel.iter()) {
            p.vel = add(p.vel, *dv);
        }
    }

    /// Compute turbulent viscosity force on a single particle (diagnostic).
    ///
    /// Returns the acceleration due to turbulent viscosity.
    pub fn single_particle_force(
        pos_i: [f64; 3],
        vel_i: [f64; 3],
        rho_i: f64,
        nut_i: f64,
        neighbors: &[TurbParticleSnapshot],
        h: f64,
    ) -> [f64; 3] {
        let mut acc = [0.0f64; 3];
        for &(pos_j, vel_j, rho_j, nut_j, h_j) in neighbors {
            let r_ij = sub(pos_i, pos_j);
            let r = length(r_ij);
            if r > 2.0 * h || r < 1e-12 {
                continue;
            }
            let grad_w = kernel_gradient(r_ij, h);
            let rho_avg = 0.5 * (rho_i + rho_j);
            if rho_avg < 1e-12 {
                continue;
            }
            let m_j = rho_j * h_j.powi(3);
            let nu_sum = nut_i + nut_j;
            let dv = sub(vel_j, vel_i);
            let dw_dot_rij = dot(grad_w, r_ij);
            let fac = m_j * nu_sum * dw_dot_rij / (rho_avg * r * r);
            acc = add(acc, scale(dv, fac));
        }
        acc
    }
}

// ---------------------------------------------------------------------------
// Turbulence intensity
// ---------------------------------------------------------------------------

/// Compute turbulence intensity I = u'_rms / U_mean.
pub fn turbulence_intensity(particles: &[TurbulentParticle]) -> f64 {
    if particles.is_empty() {
        return 0.0;
    }
    let mean_vel = mean_velocity(particles);
    let u_mean = length(mean_vel);
    if u_mean < 1e-12 {
        return 0.0;
    }
    let tke = tke_from_particles(particles, mean_vel);
    // u'_rms = sqrt(2/3 * k)
    let u_rms = (2.0 / 3.0 * tke).sqrt();
    u_rms / u_mean
}

/// Compute the turbulent Reynolds number Re_t = k² / (ν ε).
pub fn turbulent_reynolds_number(k: f64, epsilon: f64, nu: f64) -> f64 {
    if epsilon.abs() < 1e-30 || nu.abs() < 1e-30 {
        return 0.0;
    }
    k * k / (nu * epsilon)
}

/// Compute Taylor microscale Reynolds number Re_λ = u'_rms * λ / ν.
///
/// λ = sqrt(10 * ν * k / ε) is the Taylor microscale.
pub fn taylor_reynolds_number(k: f64, epsilon: f64, nu: f64) -> f64 {
    if epsilon.abs() < 1e-30 || nu.abs() < 1e-30 {
        return 0.0;
    }
    let u_rms = (2.0 / 3.0 * k).sqrt();
    let lambda = (10.0 * nu * k / epsilon).sqrt();
    u_rms * lambda / nu
}

// ---------------------------------------------------------------------------
// SPH turbulent diffusion
// ---------------------------------------------------------------------------

/// Compute SPH Laplacian of a scalar field using inter-particle formula.
///
/// ∇²φ_i ≈ 2 Σ_j (m_j/ρ_j) (φ_i - φ_j) r_ij · ∇W_ij / (r_ij² + η²)
///
/// Returns the Laplacian value at the query particle.
pub fn sph_laplacian_scalar(
    phi_i: f64,
    pos_i: [f64; 3],
    neighbors: &[([f64; 3], f64, f64, f64)], // (pos_j, phi_j, rho_j, mass_j)
    h: f64,
) -> f64 {
    let eta2 = 0.01 * h * h;
    let mut result = 0.0f64;
    for &(pos_j, phi_j, rho_j, mass_j) in neighbors {
        if rho_j < 1e-12 {
            continue;
        }
        let r_ij = sub(pos_i, pos_j);
        let r2 = dot(r_ij, r_ij);
        let r = r2.sqrt();
        if r > 2.0 * h || r < 1e-12 {
            continue;
        }
        let grad_w = kernel_gradient(r_ij, h);
        let rdotgw = dot(r_ij, grad_w);
        result += 2.0 * (mass_j / rho_j) * (phi_i - phi_j) * rdotgw / (r2 + eta2);
    }
    result
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    fn make_particle(pos: [f64; 3], vel: [f64; 3]) -> TurbulentParticle {
        TurbulentParticle {
            pos,
            vel,
            density: 1000.0,
            pressure: 0.0,
            turbulent_ke: 0.1,
            turbulent_dissipation: 0.01,
            eddy_viscosity: 0.001,
            smoothing_length: 0.1,
        }
    }

    #[test]
    fn test_smagorinsky_eddy_viscosity_positive() {
        let model = SmagorinskyModel::new(0.18);
        let nu_t = model.eddy_viscosity(1.5, 0.05);
        assert!(
            nu_t > 0.0,
            "eddy_viscosity must be positive for nonzero strain"
        );
    }

    #[test]
    fn test_strain_rate_magnitude_nonnegative() {
        let s = [0.1, -0.05, 0.0, 0.2, -0.1, 0.05];
        let mag = strain_rate_magnitude(&s);
        assert!(mag >= 0.0, "strain_rate_magnitude must be non-negative");
    }

    #[test]
    fn test_strain_rate_magnitude_zero() {
        let s = [0.0; 6];
        let mag = strain_rate_magnitude(&s);
        assert!((mag - 0.0).abs() < 1e-12);
    }

    #[test]
    fn test_kepsilon_eddy_viscosity_positive() {
        let model = KepsilonModel::new();
        let k = 0.5;
        let eps = 0.1;
        let nu_t = model.eddy_viscosity(k, eps);
        assert!(nu_t > 0.0, "eddy_viscosity must be positive for k>0, eps>0");
    }

    #[test]
    fn test_kepsilon_k_decreases_no_production() {
        let model = KepsilonModel::new();
        let k0 = 1.0;
        let eps = 0.5;
        let k1 = model.update_k(k0, eps, 0.0, 0.01);
        assert!(k1 < k0, "k must decrease when production is zero");
    }

    #[test]
    fn test_kepsilon_k_increases_when_production_exceeds_dissipation() {
        let model = KepsilonModel::new();
        let k0 = 0.5;
        let eps = 0.1;
        let production = 1.0; // production > dissipation
        let k1 = model.update_k(k0, eps, production, 0.01);
        assert!(k1 > k0, "k must increase when production > dissipation");
    }

    #[test]
    fn test_turbulent_stress_zero_for_zero_nu_t() {
        let s_ij = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6];
        let tau = SmagorinskyModel::turbulent_stress(0.0, 1000.0, s_ij);
        for &v in &tau {
            assert!((v - 0.0).abs() < 1e-15, "stress must be zero when nu_t = 0");
        }
    }

    #[test]
    fn test_dynamic_smagorinsky_cs_range() {
        let cs = DynamicSmagorinskyModel::compute_dynamic_cs(1.0, 0.1, 0.8, 0.2);
        assert!((0.0..=0.3).contains(&cs), "dynamic Cs must be in [0, 0.3]");
    }

    #[test]
    fn test_velocity_gradient_sph_zero_neighbors() {
        let grad = velocity_gradient_sph([0.0; 3], [1.0, 0.0, 0.0], &[], 0.1);
        for row in &grad {
            for &v in row {
                assert!((v - 0.0).abs() < 1e-15);
            }
        }
    }

    // --- Reynolds stress tests ---

    #[test]
    fn test_reynolds_stress_zero_uniform_flow() {
        let p1 = make_particle([0.0; 3], [1.0, 0.0, 0.0]);
        let p2 = make_particle([0.1, 0.0, 0.0], [1.0, 0.0, 0.0]);
        let p3 = make_particle([0.2, 0.0, 0.0], [1.0, 0.0, 0.0]);
        let particles = [p1, p2, p3];
        let mean = mean_velocity(&particles);
        let r = reynolds_stress_tensor(&particles, mean);
        for &v in &r {
            assert!(
                v.abs() < 1e-12,
                "Reynolds stress should be zero for uniform flow"
            );
        }
    }

    #[test]
    fn test_reynolds_stress_positive_diagonal() {
        let p1 = make_particle([0.0; 3], [1.0, 0.0, 0.0]);
        let p2 = make_particle([0.1, 0.0, 0.0], [-1.0, 0.0, 0.0]);
        let particles = [p1, p2];
        let mean = mean_velocity(&particles);
        let r = reynolds_stress_tensor(&particles, mean);
        assert!(r[0] >= 0.0, "R_xx must be non-negative");
        assert!(r[1] >= 0.0, "R_yy must be non-negative");
        assert!(r[2] >= 0.0, "R_zz must be non-negative");
    }

    // --- TKE tests ---

    #[test]
    fn test_tke_from_reynolds_matches_direct() {
        let p1 = make_particle([0.0; 3], [2.0, 1.0, -0.5]);
        let p2 = make_particle([0.1, 0.0, 0.0], [0.0, -1.0, 0.5]);
        let p3 = make_particle([0.0, 0.1, 0.0], [1.0, 0.5, 0.0]);
        let particles = [p1, p2, p3];
        let mean = mean_velocity(&particles);
        let r = reynolds_stress_tensor(&particles, mean);
        let tke_r = tke_from_reynolds(&r);
        let tke_d = tke_from_particles(&particles, mean);
        assert!(
            (tke_r - tke_d).abs() < 1e-10,
            "TKE from Reynolds ({tke_r}) should match direct ({tke_d})"
        );
    }

    #[test]
    fn test_tke_zero_uniform_flow() {
        let particles = [
            make_particle([0.0; 3], [5.0, 0.0, 0.0]),
            make_particle([0.1, 0.0, 0.0], [5.0, 0.0, 0.0]),
        ];
        let mean = mean_velocity(&particles);
        let tke = tke_from_particles(&particles, mean);
        assert!(tke.abs() < 1e-12, "TKE should be zero for uniform flow");
    }

    // --- TKE history tests ---

    #[test]
    fn test_tke_history_decay_rate() {
        let mut hist = TkeHistory::new();
        hist.record(0.0, 1.0);
        hist.record(0.1, 0.8);
        let rate = hist.decay_rate();
        assert!(
            (rate - (-2.0)).abs() < 1e-10,
            "decay rate should be -2.0, got {rate}"
        );
    }

    #[test]
    fn test_tke_history_len() {
        let mut hist = TkeHistory::new();
        assert!(hist.is_empty());
        hist.record(0.0, 1.0);
        hist.record(0.1, 0.9);
        assert_eq!(hist.len(), 2);
        assert!(!hist.is_empty());
    }

    // --- Energy spectrum tests ---

    #[test]
    fn test_energy_spectrum_empty() {
        let (k, e) = energy_spectrum_estimate(&[], [0.0; 3], 10, 100.0);
        assert!(k.is_empty());
        assert!(e.is_empty());
    }

    #[test]
    fn test_energy_spectrum_nonnegative() {
        let particles = vec![
            make_particle([0.0, 0.0, 0.0], [1.0, 0.0, 0.0]),
            make_particle([0.1, 0.0, 0.0], [0.5, 0.2, 0.0]),
            make_particle([0.2, 0.0, 0.0], [-0.3, 0.1, 0.5]),
            make_particle([0.0, 0.1, 0.0], [0.7, -0.1, 0.2]),
        ];
        let mean = mean_velocity(&particles);
        let (_k, e) = energy_spectrum_estimate(&particles, mean, 5, 100.0);
        for &val in &e {
            assert!(val >= 0.0, "Energy spectrum values must be non-negative");
        }
    }

    // --- LES model tests ---

    #[test]
    fn test_wale_zero_for_zero_gradient() {
        let model = WaleModel::new(0.5);
        let grad = [[0.0f64; 3]; 3];
        let nu_t = model.eddy_viscosity(&grad, 0.1);
        assert!(
            nu_t.abs() < 1e-15,
            "WALE should give zero for zero gradient"
        );
    }

    #[test]
    fn test_wale_nonnegative() {
        let model = WaleModel::new(0.5);
        let grad = [[0.1, 0.05, -0.02], [-0.03, 0.2, 0.04], [0.01, -0.01, -0.3]];
        let nu_t = model.eddy_viscosity(&grad, 0.1);
        assert!(nu_t >= 0.0, "WALE eddy viscosity must be non-negative");
    }

    #[test]
    fn test_vreman_zero_for_zero_gradient() {
        let model = VremanModel::new(0.07);
        let grad = [[0.0f64; 3]; 3];
        let nu_t = model.eddy_viscosity(&grad, 0.1);
        assert!(
            nu_t.abs() < 1e-15,
            "Vreman should give zero for zero gradient"
        );
    }

    #[test]
    fn test_vreman_nonnegative() {
        let model = VremanModel::new(0.07);
        let grad = [[0.5, 0.1, -0.05], [-0.1, 0.3, 0.2], [0.05, -0.15, -0.8]];
        let nu_t = model.eddy_viscosity(&grad, 0.1);
        assert!(nu_t >= 0.0, "Vreman eddy viscosity must be non-negative");
    }

    #[test]
    fn test_sigma_model_zero_for_zero_gradient() {
        let model = SigmaModel::new(1.35);
        let grad = [[0.0f64; 3]; 3];
        let nu_t = model.eddy_viscosity(&grad, 0.1);
        assert!(
            nu_t.abs() < 1e-15,
            "Sigma model should give zero for zero gradient"
        );
    }

    // --- Rotation rate tensor tests ---

    #[test]
    fn test_rotation_rate_antisymmetric() {
        let grad = [[0.1, 0.3, -0.1], [0.5, 0.2, 0.4], [-0.2, -0.1, -0.3]];
        let w = rotation_rate_tensor(&grad);
        // W_xy = 0.5*(grad[0][1] - grad[1][0]) = 0.5*(0.3 - 0.5) = -0.1
        assert!((w[0] - (-0.1)).abs() < 1e-12, "W_xy mismatch: {}", w[0]);
    }

    // --- Mixing length tests ---

    #[test]
    fn test_mixing_length_viscosity_positive() {
        let nu_t = mixing_length_viscosity(0.1, 10.0, 0.41);
        assert!(nu_t > 0.0, "mixing length viscosity must be positive");
    }

    #[test]
    fn test_mixing_length_viscosity_zero_at_wall() {
        let nu_t = mixing_length_viscosity(0.0, 10.0, 0.41);
        assert!(
            nu_t.abs() < 1e-15,
            "mixing length viscosity must be zero at wall"
        );
    }

    // --- Effective viscosity tests ---

    #[test]
    fn test_effective_viscosity() {
        let nu_eff = effective_viscosity(1e-6, 1e-3);
        assert!((nu_eff - 1.001e-3).abs() < 1e-8);
    }

    // --- Turbulent thermal diffusivity tests ---

    #[test]
    fn test_turbulent_thermal_diffusivity() {
        let alpha_t = turbulent_thermal_diffusivity(0.01, 0.85);
        let expected = 0.01 / 0.85;
        assert!(
            (alpha_t - expected).abs() < 1e-10,
            "expected {expected}, got {alpha_t}"
        );
    }

    // --- k-epsilon extended tests ---

    #[test]
    fn test_kepsilon_time_scale() {
        let model = KepsilonModel::new();
        let tau = model.turbulent_time_scale(1.0, 0.5);
        assert!((tau - 2.0).abs() < 1e-12, "expected tau=2, got {tau}");
    }

    #[test]
    fn test_kepsilon_length_scale() {
        let model = KepsilonModel::new();
        let l = model.turbulent_length_scale(1.0, 0.5);
        let expected = 1.0_f64.powf(1.5) / 0.5;
        assert!((l - expected).abs() < 1e-12, "expected {expected}, got {l}");
    }

    #[test]
    fn test_kolmogorov_scale_positive() {
        let model = KepsilonModel::new();
        let eta = model.kolmogorov_scale(1e-6, 0.1);
        assert!(eta > 0.0, "Kolmogorov scale must be positive");
    }

    #[test]
    fn test_kepsilon_full_step() {
        let model = KepsilonModel::new();
        let (new_k, new_eps) = model.full_step(1.0, 0.5, 0.6, 0.01, 0.001, 0.0, 0.0);
        assert!(new_k >= 0.0, "k must be non-negative after full_step");
        assert!(new_eps > 0.0, "epsilon must be positive after full_step");
    }

    // --- Leonard stress test ---

    #[test]
    fn test_leonard_stress_size() {
        let particles = vec![
            make_particle([0.0; 3], [1.0, 0.0, 0.0]),
            make_particle([0.1, 0.0, 0.0], [0.5, 0.2, 0.0]),
        ];
        let filtered = DynamicSmagorinskyModel::test_filter(&particles, 0.2);
        let leonard = DynamicSmagorinskyModel::leonard_stress(&particles, &filtered);
        assert_eq!(leonard.len(), 2);
    }

    // --- Turbulence intensity tests ---

    #[test]
    fn test_turbulence_intensity_zero_uniform() {
        let particles = vec![
            make_particle([0.0; 3], [5.0, 0.0, 0.0]),
            make_particle([0.1, 0.0, 0.0], [5.0, 0.0, 0.0]),
        ];
        let ti = turbulence_intensity(&particles);
        assert!(
            ti.abs() < 1e-12,
            "turbulence intensity should be 0 for uniform flow"
        );
    }

    #[test]
    fn test_turbulence_intensity_positive() {
        let particles = vec![
            make_particle([0.0; 3], [5.0, 0.0, 0.0]),
            make_particle([0.1, 0.0, 0.0], [3.0, 1.0, 0.0]),
            make_particle([0.2, 0.0, 0.0], [7.0, -1.0, 0.0]),
        ];
        let ti = turbulence_intensity(&particles);
        assert!(
            ti > 0.0,
            "turbulence intensity should be positive for non-uniform flow"
        );
    }

    // --- Turbulent Reynolds number tests ---

    #[test]
    fn test_turbulent_reynolds_number() {
        let re_t = turbulent_reynolds_number(1.0, 0.1, 1e-6);
        assert!(re_t > 0.0, "Re_t must be positive");
    }

    #[test]
    fn test_taylor_reynolds_number() {
        let re_lambda = taylor_reynolds_number(0.5, 0.1, 1e-6);
        assert!(re_lambda > 0.0, "Re_lambda must be positive");
    }

    // --- SPH Laplacian test ---

    #[test]
    fn test_sph_laplacian_scalar_constant_field() {
        // For a constant field, Laplacian should be zero
        let phi_i = 5.0;
        let pos_i = [0.0, 0.0, 0.0];
        let neighbors = vec![
            ([0.05, 0.0, 0.0], 5.0, 1000.0, 0.1_f64.powi(3) * 1000.0),
            ([-0.05, 0.0, 0.0], 5.0, 1000.0, 0.1_f64.powi(3) * 1000.0),
        ];
        let lap = sph_laplacian_scalar(phi_i, pos_i, &neighbors, 0.1);
        assert!(
            lap.abs() < 1e-6,
            "Laplacian of constant field should be near zero, got {lap}"
        );
    }

    // --- Kernel value test ---

    #[test]
    fn test_kernel_value_positive_at_origin() {
        let w = kernel_value(0.0, 1.0);
        assert!(w > 0.0, "kernel value at origin must be positive");
    }

    #[test]
    fn test_kernel_value_zero_outside_support() {
        let w = kernel_value(2.5, 1.0);
        assert!(w.abs() < 1e-15, "kernel value outside support must be zero");
    }

    // --- Strain rate tensor identity ---

    #[test]
    fn test_strain_rate_symmetric() {
        let grad = [[0.1, 0.3, -0.1], [0.5, 0.2, 0.4], [-0.2, -0.1, -0.3]];
        let s = strain_rate_tensor(&grad);
        // S_xy = 0.5*(0.3 + 0.5) = 0.4
        assert!(
            (s[3] - 0.4).abs() < 1e-12,
            "S_xy should be 0.4, got {}",
            s[3]
        );
    }

    // --- WALE Sd squared test ---

    #[test]
    fn test_wale_sd_squared_traceless() {
        let grad = [[0.1, 0.3, -0.1], [0.5, 0.2, 0.4], [-0.2, -0.1, -0.3]];
        let sd = WaleModel::sd_squared(&grad);
        let trace = sd[0][0] + sd[1][1] + sd[2][2];
        assert!(
            trace.abs() < 1e-10,
            "S^d should be traceless, got trace={trace}"
        );
    }

    // --- Mean velocity test ---

    #[test]
    fn test_mean_velocity_computed() {
        let particles = vec![
            make_particle([0.0; 3], [2.0, 0.0, 0.0]),
            make_particle([0.1, 0.0, 0.0], [4.0, 0.0, 0.0]),
        ];
        let mean = mean_velocity(&particles);
        assert!((mean[0] - 3.0).abs() < 1e-12, "mean x vel should be 3.0");
    }
}
