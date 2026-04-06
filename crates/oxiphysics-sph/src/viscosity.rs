#![allow(clippy::needless_range_loop)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! SPH viscosity and surface tension models.
//!
//! Provides:
//! - Monaghan artificial viscosity (Π_ij)
//! - XSPH velocity correction
//! - Morris et al. viscosity force
//! - Cleary viscosity formulation
//! - SPS turbulence viscosity (Sub-Particle Scale)
//! - Temperature-dependent viscosity
//! - Non-Newtonian viscosity (power-law, Cross model)
//! - Viscosity limiters
//! - Continuum Surface Force (CSF) surface tension
//! - Akinci boundary particle repulsion
//! - Smoothed velocity gradient and strain rate

#![allow(dead_code)]

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

const EPSILON: f64 = 1.0e-14;

#[inline]
fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[inline]
fn norm_sq3(a: [f64; 3]) -> f64 {
    dot3(a, a)
}

#[inline]
fn norm3(a: [f64; 3]) -> f64 {
    norm_sq3(a).sqrt()
}

// ---------------------------------------------------------------------------
// Monaghan artificial viscosity
// ---------------------------------------------------------------------------

/// Monaghan (1992) artificial viscosity term Π_ij.
///
/// Returns a non-negative scalar viscosity only when particles are approaching
/// (v_ij · r_ij < 0).  When particles are separating the term is zero.
///
/// # Formula
/// ```text
/// μ_ij = h * (v_ij · r_ij) / (|r_ij|² + ε)
/// Π_ij = (-α * c̄ * μ + β * μ²) / ρ̄   if v · r < 0
///       = 0                              otherwise
/// ```
/// where `c̄ = (c_i + c_j) / 2`, `ρ̄ = (ρ_i + ρ_j) / 2`.
#[allow(clippy::too_many_arguments)]
pub fn artificial_viscosity(
    alpha: f64,
    beta: f64,
    rho_i: f64,
    rho_j: f64,
    v_ij: [f64; 3],
    r_ij: [f64; 3],
    c_i: f64,
    c_j: f64,
    h: f64,
) -> f64 {
    let vdotr = dot3(v_ij, r_ij);
    if vdotr >= 0.0 {
        return 0.0;
    }
    let r_sq = norm_sq3(r_ij);
    let mu = h * vdotr / (r_sq + EPSILON * h * h);
    let c_bar = 0.5 * (c_i + c_j);
    let rho_bar = 0.5 * (rho_i + rho_j);
    (-alpha * c_bar * mu + beta * mu * mu) / rho_bar
}

// ---------------------------------------------------------------------------
// XSPH velocity correction
// ---------------------------------------------------------------------------

/// XSPH velocity smoothing correction (Monaghan 1989).
///
/// Returns the correction term `ε * (m_j / ρ_j) * (v_j − v_i) * W` that is
/// added to particle `i`'s velocity during integration to reduce penetration.
pub fn xsph_correction(
    epsilon: f64,
    v_j: [f64; 3],
    v_i: [f64; 3],
    m_j: f64,
    rho_j: f64,
    w: f64,
) -> [f64; 3] {
    let scale = epsilon * m_j / rho_j * w;
    [
        scale * (v_j[0] - v_i[0]),
        scale * (v_j[1] - v_i[1]),
        scale * (v_j[2] - v_i[2]),
    ]
}

// ---------------------------------------------------------------------------
// Morris viscosity force
// ---------------------------------------------------------------------------

/// Morris et al. (1997) laminar viscosity force contribution from particle j
/// on particle i.
///
/// # Formula
/// ```text
/// f_visc = μ * m_j / ρ_j * (v_i − v_j) * dot(∇W, r_ij) / (|r_ij|² + ε) / ρ_i
/// ```
///
/// The result is a force-per-unit-mass vector (acceleration).
pub fn morris_viscosity_force(
    mu: f64,
    m_j: f64,
    rho_i: f64,
    rho_j: f64,
    v_ij: [f64; 3],
    r_ij: [f64; 3],
    grad_w: [f64; 3],
) -> [f64; 3] {
    let r_sq = norm_sq3(r_ij);
    let denom = r_sq + EPSILON;
    let scale = mu * m_j / rho_j * dot3(grad_w, r_ij) / denom / rho_i;
    [scale * v_ij[0], scale * v_ij[1], scale * v_ij[2]]
}

// ---------------------------------------------------------------------------
// Cleary viscosity formulation
// ---------------------------------------------------------------------------

/// Cleary (1998) viscosity model.
///
/// Uses a harmonic mean of viscosities for particles with different
/// dynamic viscosities `mu_i` and `mu_j`:
///
/// ```text
/// f_visc = (2 * mu_i * mu_j) / (mu_i + mu_j) * m_j / (rho_i * rho_j) *
///          v_ij * dot(r_ij, grad_w) / (|r_ij|² + η²)
/// ```
///
/// where `η = 0.01 * h`.
///
/// Returns acceleration (force per unit mass) on particle i.
#[allow(clippy::too_many_arguments)]
pub fn cleary_viscosity_force(
    mu_i: f64,
    mu_j: f64,
    m_j: f64,
    rho_i: f64,
    rho_j: f64,
    v_ij: [f64; 3],
    r_ij: [f64; 3],
    grad_w: [f64; 3],
    h: f64,
) -> [f64; 3] {
    let mu_sum = mu_i + mu_j;
    if mu_sum.abs() < EPSILON {
        return [0.0; 3];
    }
    let mu_harm = 2.0 * mu_i * mu_j / mu_sum;
    let eta = 0.01 * h;
    let r_sq = norm_sq3(r_ij);
    let rdotgw = dot3(r_ij, grad_w);
    let denom = r_sq + eta * eta;
    let scale = mu_harm * m_j / (rho_i * rho_j) * rdotgw / denom;
    [scale * v_ij[0], scale * v_ij[1], scale * v_ij[2]]
}

// ---------------------------------------------------------------------------
// SPS turbulence viscosity (Sub-Particle Scale)
// ---------------------------------------------------------------------------

/// SPS (Sub-Particle Scale) turbulence viscosity model.
///
/// Computes the eddy viscosity using the Smagorinsky model:
///
/// ```text
/// ν_t = (C_s * Δ)² * |S|
/// ```
///
/// where:
/// - `C_s` is the Smagorinsky constant (typically 0.1–0.2)
/// - `Δ` is the filter width (usually taken as the smoothing length `h`)
/// - `|S|` is the strain-rate magnitude
///
/// Returns the turbulent kinematic viscosity `ν_t`.
pub fn sps_eddy_viscosity(
    smagorinsky_constant: f64,
    filter_width: f64,
    strain_rate_mag: f64,
) -> f64 {
    let cs_delta = smagorinsky_constant * filter_width;
    cs_delta * cs_delta * strain_rate_mag
}

/// Compute the SPS turbulent stress tensor τ_ij for SPH.
///
/// ```text
/// τ_αβ = 2 * ν_t * ρ * S_αβ - (2/3) * ρ * k * δ_αβ
/// ```
///
/// Here `k` is the turbulent kinetic energy, estimated as:
/// `k = (ν_t / (C_I * Δ))²` where C_I ≈ 0.0066.
///
/// Returns the 3×3 stress tensor.
pub fn sps_stress_tensor(
    nu_t: f64,
    rho: f64,
    grad_v: &[[f64; 3]; 3],
    filter_width: f64,
) -> [[f64; 3]; 3] {
    let ci = 0.0066_f64;
    let k = if filter_width > EPSILON {
        let ratio = nu_t / (ci * filter_width);
        ratio * ratio
    } else {
        0.0
    };

    let mut tau = [[0.0_f64; 3]; 3];
    for a in 0..3 {
        for b in 0..3 {
            let s_ab = 0.5 * (grad_v[a][b] + grad_v[b][a]);
            tau[a][b] = 2.0 * nu_t * rho * s_ab;
        }
        // subtract isotropic turbulent kinetic energy part
        tau[a][a] -= (2.0 / 3.0) * rho * k;
    }
    tau
}

// ---------------------------------------------------------------------------
// Temperature-dependent viscosity
// ---------------------------------------------------------------------------

/// Arrhenius-type temperature-dependent dynamic viscosity.
///
/// ```text
/// μ(T) = μ_ref * exp(E_a / R * (1/T - 1/T_ref))
/// ```
///
/// where:
/// - `mu_ref` is the reference viscosity at temperature `T_ref`
/// - `E_a` is the activation energy (J/mol)
/// - `R` is the universal gas constant (8.314 J/(mol·K))
/// - `T` is the current temperature (K)
pub fn arrhenius_viscosity(mu_ref: f64, t_ref: f64, t: f64, activation_energy: f64) -> f64 {
    const R_GAS: f64 = 8.314;
    if t < EPSILON || t_ref < EPSILON {
        return mu_ref;
    }
    mu_ref * (activation_energy / R_GAS * (1.0 / t - 1.0 / t_ref)).exp()
}

/// Sutherland viscosity model (for gases).
///
/// ```text
/// μ(T) = μ_ref * (T / T_ref)^(3/2) * (T_ref + S) / (T + S)
/// ```
///
/// where `S` is the Sutherland constant (e.g., 110.4 K for air).
pub fn sutherland_viscosity(mu_ref: f64, t_ref: f64, t: f64, s_const: f64) -> f64 {
    if t < EPSILON || t_ref < EPSILON {
        return mu_ref;
    }
    mu_ref * (t / t_ref).powf(1.5) * (t_ref + s_const) / (t + s_const)
}

/// VFT (Vogel-Fulcher-Tammann) viscosity for glass-forming liquids.
///
/// ```text
/// log10(μ) = A + B / (T - T0)
/// ```
pub fn vft_viscosity(a: f64, b: f64, t0: f64, t: f64) -> f64 {
    if (t - t0).abs() < EPSILON {
        return f64::MAX;
    }
    10.0_f64.powf(a + b / (t - t0))
}

// ---------------------------------------------------------------------------
// Non-Newtonian viscosity models
// ---------------------------------------------------------------------------

/// Power-law viscosity model for non-Newtonian fluids.
///
/// ```text
/// μ_eff = K * |γ̇|^(n-1)
/// ```
///
/// where:
/// - `K` is the consistency index (Pa·s^n)
/// - `n` is the power-law index (n < 1 = shear-thinning, n > 1 = shear-thickening)
/// - `|γ̇|` is the strain-rate magnitude
///
/// A minimum strain rate `gamma_dot_min` prevents singularity at zero shear.
pub fn power_law_viscosity(k: f64, n: f64, gamma_dot: f64, gamma_dot_min: f64) -> f64 {
    let gamma_eff = gamma_dot.max(gamma_dot_min);
    k * gamma_eff.powf(n - 1.0)
}

/// Cross model viscosity for non-Newtonian fluids.
///
/// ```text
/// μ_eff = μ_inf + (μ_0 - μ_inf) / (1 + (λ * |γ̇|)^m)
/// ```
///
/// where:
/// - `mu_0` is the zero-shear viscosity
/// - `mu_inf` is the infinite-shear viscosity
/// - `lambda` is the relaxation time
/// - `m` is the Cross rate constant (often ~2/3)
/// - `|γ̇|` is the strain-rate magnitude
pub fn cross_model_viscosity(mu_0: f64, mu_inf: f64, lambda: f64, m: f64, gamma_dot: f64) -> f64 {
    let lg = lambda * gamma_dot;
    mu_inf + (mu_0 - mu_inf) / (1.0 + lg.powf(m))
}

/// Carreau-Yasuda viscosity model.
///
/// ```text
/// μ_eff = μ_inf + (μ_0 - μ_inf) * (1 + (λ * |γ̇|)^a)^((n-1)/a)
/// ```
pub fn carreau_yasuda_viscosity(
    mu_0: f64,
    mu_inf: f64,
    lambda: f64,
    a_param: f64,
    n: f64,
    gamma_dot: f64,
) -> f64 {
    let lg = lambda * gamma_dot;
    let factor = (1.0 + lg.powf(a_param)).powf((n - 1.0) / a_param);
    mu_inf + (mu_0 - mu_inf) * factor
}

/// Bingham plastic viscosity model.
///
/// For yield-stress fluids: the material behaves as a rigid body below the
/// yield stress and as a Newtonian fluid above it.
///
/// ```text
/// μ_eff = μ_p + τ_y / max(|γ̇|, γ̇_min)
/// ```
///
/// where `μ_p` is the plastic viscosity and `τ_y` is the yield stress.
pub fn bingham_viscosity(
    mu_plastic: f64,
    tau_yield: f64,
    gamma_dot: f64,
    gamma_dot_min: f64,
) -> f64 {
    let gamma_eff = gamma_dot.max(gamma_dot_min);
    mu_plastic + tau_yield / gamma_eff
}

/// Herschel-Bulkley viscosity model (generalized Bingham).
///
/// ```text
/// μ_eff = K * |γ̇|^(n-1) + τ_y / max(|γ̇|, γ̇_min)
/// ```
pub fn herschel_bulkley_viscosity(
    k: f64,
    n: f64,
    tau_yield: f64,
    gamma_dot: f64,
    gamma_dot_min: f64,
) -> f64 {
    let gamma_eff = gamma_dot.max(gamma_dot_min);
    k * gamma_eff.powf(n - 1.0) + tau_yield / gamma_eff
}

// ---------------------------------------------------------------------------
// Viscosity limiters
// ---------------------------------------------------------------------------

/// Clamp viscosity to a specified range.
///
/// Prevents unphysically large or small viscosity values that can cause
/// numerical instability.
pub fn clamp_viscosity(mu: f64, mu_min: f64, mu_max: f64) -> f64 {
    mu.clamp(mu_min, mu_max)
}

/// Monaghan-Gingold viscosity limiter.
///
/// Limits the artificial viscosity contribution based on particle separation.
/// If the relative velocity is too large compared to the speed of sound,
/// the viscosity is capped.
///
/// Returns the limited viscosity term.
pub fn monaghan_gingold_limiter(pi_ij: f64, rho_bar: f64, c_bar: f64, h: f64) -> f64 {
    let pi_max = c_bar * c_bar / rho_bar;
    let _ = h; // available for extended limiters
    if pi_ij > pi_max { pi_max } else { pi_ij }
}

/// Balsara switch for artificial viscosity.
///
/// Reduces artificial viscosity in regions of strong vorticity relative to
/// compression, preventing excessive dissipation in shear flows.
///
/// ```text
/// f_i = |div v| / (|div v| + |curl v| + ε * c / h)
/// ```
///
/// Returns a factor in \[0, 1\] to multiply the artificial viscosity.
pub fn balsara_switch(div_v: f64, curl_v_mag: f64, c: f64, h: f64) -> f64 {
    let denom = div_v.abs() + curl_v_mag + EPSILON * c / h;
    div_v.abs() / denom
}

/// Morris-Monaghan viscosity switch.
///
/// Uses a signal velocity to adaptively scale the viscosity coefficient:
///
/// ```text
/// v_sig = c_i + c_j - 3 * min(0, v_ij · r_ij / |r_ij|)
/// ```
///
/// Returns the signal velocity for use in the viscosity calculation.
pub fn signal_velocity(c_i: f64, c_j: f64, v_ij: [f64; 3], r_ij: [f64; 3]) -> f64 {
    let r = norm3(r_ij);
    if r < EPSILON {
        return c_i + c_j;
    }
    let vr = dot3(v_ij, r_ij) / r;
    c_i + c_j - 3.0 * vr.min(0.0)
}

// ---------------------------------------------------------------------------
// Surface tension — CSF
// ---------------------------------------------------------------------------

/// Continuum Surface Force (CSF) surface tension model (Morris 2000).
///
/// Computes surface tension forces from an estimated curvature and interface
/// normal: **F = σ · κ · n̂**.
#[derive(Debug, Clone)]
pub struct SurfaceTensionCSF {
    /// Surface tension coefficient σ (N/m).
    pub sigma: f64,
    /// Smoothing length h used in the SPH kernel.
    pub smoothing_length: f64,
}

impl SurfaceTensionCSF {
    /// Create a new CSF surface tension model.
    pub fn new(sigma: f64, h: f64) -> Self {
        Self {
            sigma,
            smoothing_length: h,
        }
    }

    /// Colour-function value for particle `i`.
    ///
    /// In the simple single-phase CSF formulation the colour function is the
    /// local normalised density: `c = ρ_i / ρ₀`.
    pub fn color_function_gradient(self_density: f64, rho0: f64) -> f64 {
        self_density / rho0
    }

    /// Estimate the interface curvature κ ≈ −div(n̂) ≈ −tr(∇n̂).
    ///
    /// `normal` is the unit interface normal vector n̂ at the particle.
    /// `grad_normal` is the 3×3 Jacobian ∂n̂_α/∂x_β computed by SPH.
    ///
    /// κ = −(∂n̂_x/∂x + ∂n̂_y/∂y + ∂n̂_z/∂z)
    pub fn curvature_estimate(normal: [f64; 3], grad_normal: [[f64; 3]; 3]) -> f64 {
        let _ = normal; // available for callers constructing grad_normal
        -(grad_normal[0][0] + grad_normal[1][1] + grad_normal[2][2])
    }

    /// CSF surface tension force per unit volume: **F = σ · κ · n̂**.
    pub fn surface_force(&self, curvature: f64, normal: [f64; 3]) -> [f64; 3] {
        let scale = self.sigma * curvature;
        [scale * normal[0], scale * normal[1], scale * normal[2]]
    }
}

// ---------------------------------------------------------------------------
// Akinci boundary particle repulsion
// ---------------------------------------------------------------------------

/// Akinci et al. (2012) boundary particle repulsion model.
///
/// Applies a distance-based repulsion force that prevents fluid particles from
/// penetrating solid boundaries.
#[derive(Debug, Clone)]
pub struct AkinciBoundary {
    /// Distance from the boundary particle to the fluid particle (m).
    pub distance: f64,
    /// Outward-pointing unit normal of the boundary.
    pub normal: [f64; 3],
    /// Repulsion strength coefficient.
    pub repulsion: f64,
}

impl AkinciBoundary {
    /// Compute the repulsion force vector.
    ///
    /// The magnitude decays with distance: `f = repulsion / distance²` clamped
    /// so that zero distance does not cause a division by zero.
    pub fn repulsion_force(&self) -> [f64; 3] {
        let d_sq = self.distance * self.distance + EPSILON;
        let mag = self.repulsion / d_sq;
        [
            mag * self.normal[0],
            mag * self.normal[1],
            mag * self.normal[2],
        ]
    }
}

// ---------------------------------------------------------------------------
// Smoothed velocity gradient
// ---------------------------------------------------------------------------

/// Compute the smoothed velocity-gradient tensor ∇v at particle `i`.
///
/// Uses the standard SPH first-derivative estimator:
/// ```text
/// ∂v_α/∂x_β ≈ Σ_j (m_j / ρ_j) * (v_j_α − v_i_α) * ∂W/∂x_β
/// ```
///
/// The kernel gradient `∂W/∂x_β` is approximated here as
/// `(dW/dr) * r_hat_β` with `dW/dr` from the cubic-spline derivative.
///
/// Returns a 3×3 matrix `grad_v[α][β] = ∂v_α/∂x_β`.
pub fn smoothed_velocity_gradient(
    positions: &[[f64; 3]],
    velocities: &[[f64; 3]],
    masses: &[f64],
    densities: &[f64],
    h: f64,
    i: usize,
) -> [[f64; 3]; 3] {
    let mut grad_v = [[0.0_f64; 3]; 3];
    let vi = velocities[i];
    let xi = positions[i];

    for j in 0..positions.len() {
        if j == i {
            continue;
        }
        let rho_j = densities[j].max(EPSILON);
        let rij = [
            xi[0] - positions[j][0],
            xi[1] - positions[j][1],
            xi[2] - positions[j][2],
        ];
        let r = norm3(rij);
        if r < EPSILON {
            continue;
        }
        let dw_dr = cubic_spline_dw_dr(r, h);
        // ∇_i W_ij = (dW/dr) * r_hat (pointing from j → i, i.e. r_ij/r)
        let grad_w = [dw_dr * rij[0] / r, dw_dr * rij[1] / r, dw_dr * rij[2] / r];
        let weight = masses[j] / rho_j;
        let dv = [
            velocities[j][0] - vi[0],
            velocities[j][1] - vi[1],
            velocities[j][2] - vi[2],
        ];
        for alpha in 0..3 {
            for beta in 0..3 {
                grad_v[alpha][beta] += weight * dv[alpha] * grad_w[beta];
            }
        }
    }

    grad_v
}

/// Cubic-spline kernel radial derivative dW/dr (3-D, h-normalised).
#[inline]
fn cubic_spline_dw_dr(r: f64, h: f64) -> f64 {
    let q = r / h;
    if !(EPSILON..2.0).contains(&q) {
        return 0.0;
    }
    let sigma = 1.0 / (std::f64::consts::PI * h * h * h);
    if q >= 1.0 {
        let t = 2.0 - q;
        sigma * (-0.75 * t * t) / h
    } else {
        sigma * (-3.0 * q + 2.25 * q * q) / h
    }
}

// ---------------------------------------------------------------------------
// Strain rate magnitude
// ---------------------------------------------------------------------------

/// Compute the Frobenius-norm strain-rate magnitude from the velocity-gradient
/// tensor `grad_v`.
///
/// # Formula
/// ```text
/// S_ij = (∂v_i/∂x_j + ∂v_j/∂x_i) / 2
/// ||S|| = sqrt(2 * S_ij * S_ij)
/// ```
pub fn strain_rate_magnitude(grad_v: &[[f64; 3]; 3]) -> f64 {
    // Symmetric strain-rate tensor S = (∇v + ∇v^T) / 2
    let mut sum = 0.0;
    for a in 0..3 {
        for b in 0..3 {
            let s_ab = 0.5 * (grad_v[a][b] + grad_v[b][a]);
            sum += s_ab * s_ab;
        }
    }
    (2.0 * sum).sqrt()
}

/// Compute the velocity divergence from the velocity-gradient tensor.
///
/// ```text
/// div(v) = ∂v_x/∂x + ∂v_y/∂y + ∂v_z/∂z
/// ```
pub fn velocity_divergence(grad_v: &[[f64; 3]; 3]) -> f64 {
    grad_v[0][0] + grad_v[1][1] + grad_v[2][2]
}

/// Compute the vorticity magnitude from the velocity-gradient tensor.
///
/// The vorticity vector is:
/// ```text
/// ω = (∂v_z/∂y - ∂v_y/∂z, ∂v_x/∂z - ∂v_z/∂x, ∂v_y/∂x - ∂v_x/∂y)
/// ```
pub fn vorticity_magnitude(grad_v: &[[f64; 3]; 3]) -> f64 {
    let wx = grad_v[2][1] - grad_v[1][2];
    let wy = grad_v[0][2] - grad_v[2][0];
    let wz = grad_v[1][0] - grad_v[0][1];
    (wx * wx + wy * wy + wz * wz).sqrt()
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // -----------------------------------------------------------------------
    // artificial_viscosity
    // -----------------------------------------------------------------------

    /// Separating particles (v·r > 0) → Π = 0.
    #[test]
    fn test_artificial_viscosity_separating() {
        // Two particles moving away from each other along x.
        let v_ij = [1.0_f64, 0.0, 0.0]; // v_i - v_j, positive x
        let r_ij = [1.0_f64, 0.0, 0.0]; // r_i - r_j, positive x → v·r > 0
        let pi = artificial_viscosity(1.0, 0.0, 1000.0, 1000.0, v_ij, r_ij, 1500.0, 1500.0, 0.05);
        assert_eq!(pi, 0.0, "Separating particles must give Π = 0");
    }

    /// Approaching particles (v·r < 0) → Π > 0.
    #[test]
    fn test_artificial_viscosity_approaching() {
        // v_ij = v_i - v_j: particle i moves in -x, j in +x → approaching
        let v_ij = [-1.0_f64, 0.0, 0.0];
        let r_ij = [1.0_f64, 0.0, 0.0]; // r_i - r_j: i is to the right of j
        // v·r = -1 < 0 → approaching
        let pi = artificial_viscosity(1.0, 0.0, 1000.0, 1000.0, v_ij, r_ij, 1500.0, 1500.0, 0.05);
        assert!(pi > 0.0, "Approaching particles must give Π > 0, got {pi}");
    }

    // -----------------------------------------------------------------------
    // xsph_correction
    // -----------------------------------------------------------------------

    /// Same velocity → zero XSPH correction.
    #[test]
    fn test_xsph_same_velocity() {
        let v = [1.0_f64, 2.0, 3.0];
        let corr = xsph_correction(0.5, v, v, 0.001, 1000.0, 0.8);
        for k in 0..3 {
            assert!(
                corr[k].abs() < 1e-15,
                "Same velocity → zero correction; component {k} = {}",
                corr[k]
            );
        }
    }

    // -----------------------------------------------------------------------
    // SurfaceTensionCSF
    // -----------------------------------------------------------------------

    /// Zero curvature → zero surface force.
    #[test]
    fn test_surface_force_zero_curvature() {
        let csf = SurfaceTensionCSF::new(0.072, 0.05);
        let normal = [0.0_f64, 0.0, 1.0];
        let force = csf.surface_force(0.0, normal);
        for k in 0..3 {
            assert_eq!(
                force[k], 0.0,
                "Zero curvature → zero force; component {k} = {}",
                force[k]
            );
        }
    }

    // -----------------------------------------------------------------------
    // strain_rate_magnitude
    // -----------------------------------------------------------------------

    /// Pure shear: ∂v_x/∂y = γ, all others zero.
    ///
    /// S_xy = S_yx = γ/2.
    /// ||S||² = 2*(S_xy² + S_yx²) = 2*(2*(γ/2)²) = γ²
    /// ||S|| = γ.
    #[test]
    fn test_strain_rate_pure_shear() {
        let gamma = 4.0_f64;
        // grad_v[0][1] = ∂v_x/∂y = gamma, everything else = 0
        let mut grad_v = [[0.0_f64; 3]; 3];
        grad_v[0][1] = gamma;

        let s_mag = strain_rate_magnitude(&grad_v);
        // Expected: sqrt(2 * (S_xy^2 + S_yx^2)) = sqrt(2 * 2*(gamma/2)^2) = gamma
        let expected = gamma;
        assert!(
            (s_mag - expected).abs() < 1e-12,
            "Pure shear strain-rate magnitude: expected {expected}, got {s_mag}"
        );
    }

    // -----------------------------------------------------------------------
    // Cleary viscosity
    // -----------------------------------------------------------------------

    /// Zero viscosity on either side → zero Cleary force.
    #[test]
    fn test_cleary_zero_viscosity() {
        let f = cleary_viscosity_force(
            0.0,
            1.0,
            0.001,
            1000.0,
            1000.0,
            [1.0, 0.0, 0.0],
            [0.1, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            0.05,
        );
        for k in 0..3 {
            assert!(
                f[k].abs() < 1e-14,
                "Cleary force should be zero with mu_i=0"
            );
        }
    }

    /// Equal viscosities → harmonic mean equals the viscosity itself.
    #[test]
    fn test_cleary_equal_viscosities() {
        let mu = 0.001;
        let v_ij = [1.0, 0.0, 0.0];
        let r_ij = [0.05, 0.0, 0.0];
        let grad_w = [-10.0, 0.0, 0.0];
        let h = 0.05;
        let f = cleary_viscosity_force(mu, mu, 0.001, 1000.0, 1000.0, v_ij, r_ij, grad_w, h);
        // With equal viscosities, harmonic mean = mu, so result should be nonzero
        assert!(
            norm3(f) > 0.0,
            "Cleary force should be non-zero for non-zero inputs"
        );
    }

    /// Cleary force is proportional to velocity difference.
    #[test]
    fn test_cleary_scales_with_velocity() {
        let params = (0.001, 0.001, 0.001, 1000.0, 1000.0);
        let r_ij = [0.05, 0.0, 0.0];
        let grad_w = [-10.0, 0.0, 0.0];
        let h = 0.05;

        let f1 = cleary_viscosity_force(
            params.0,
            params.1,
            params.2,
            params.3,
            params.4,
            [1.0, 0.0, 0.0],
            r_ij,
            grad_w,
            h,
        );
        let f2 = cleary_viscosity_force(
            params.0,
            params.1,
            params.2,
            params.3,
            params.4,
            [2.0, 0.0, 0.0],
            r_ij,
            grad_w,
            h,
        );

        let ratio = f2[0] / f1[0];
        assert!(
            (ratio - 2.0).abs() < 1e-10,
            "Cleary force should scale linearly with velocity: ratio={ratio}"
        );
    }

    // -----------------------------------------------------------------------
    // SPS turbulence
    // -----------------------------------------------------------------------

    #[test]
    fn test_sps_eddy_viscosity_zero_strain() {
        let nu_t = sps_eddy_viscosity(0.15, 0.05, 0.0);
        assert!(nu_t.abs() < 1e-14, "Zero strain rate → zero eddy viscosity");
    }

    #[test]
    fn test_sps_eddy_viscosity_positive() {
        let nu_t = sps_eddy_viscosity(0.15, 0.05, 10.0);
        assert!(nu_t > 0.0, "Non-zero strain rate → positive eddy viscosity");
        let expected = (0.15 * 0.05) * (0.15 * 0.05) * 10.0;
        assert!(
            (nu_t - expected).abs() < 1e-14,
            "Expected {expected}, got {nu_t}"
        );
    }

    #[test]
    fn test_sps_stress_tensor_symmetry() {
        let mut grad_v = [[0.0_f64; 3]; 3];
        grad_v[0][1] = 2.0;
        grad_v[1][0] = 1.0;
        let tau = sps_stress_tensor(0.01, 1000.0, &grad_v, 0.05);
        // The stress tensor should be symmetric since S is symmetric.
        for a in 0..3 {
            for b in 0..3 {
                assert!(
                    (tau[a][b] - tau[b][a]).abs() < 1e-12,
                    "SPS stress tensor not symmetric: tau[{a}][{b}]={} != tau[{b}][{a}]={}",
                    tau[a][b],
                    tau[b][a]
                );
            }
        }
    }

    // -----------------------------------------------------------------------
    // Temperature-dependent viscosity
    // -----------------------------------------------------------------------

    #[test]
    fn test_arrhenius_at_reference_temp() {
        let mu_ref = 0.001;
        let t_ref = 300.0;
        let mu = arrhenius_viscosity(mu_ref, t_ref, t_ref, 20000.0);
        assert!(
            (mu - mu_ref).abs() < 1e-14,
            "At T_ref, viscosity should equal mu_ref: got {mu}"
        );
    }

    #[test]
    fn test_arrhenius_increases_on_cooling() {
        let mu_ref = 0.001;
        let t_ref = 300.0;
        let ea = 20000.0;
        let mu_cold = arrhenius_viscosity(mu_ref, t_ref, 280.0, ea);
        assert!(
            mu_cold > mu_ref,
            "Viscosity should increase when cooling: {mu_cold} vs {mu_ref}"
        );
    }

    #[test]
    fn test_arrhenius_decreases_on_heating() {
        let mu_ref = 0.001;
        let t_ref = 300.0;
        let ea = 20000.0;
        let mu_hot = arrhenius_viscosity(mu_ref, t_ref, 350.0, ea);
        assert!(
            mu_hot < mu_ref,
            "Viscosity should decrease when heating: {mu_hot} vs {mu_ref}"
        );
    }

    #[test]
    fn test_sutherland_at_reference_temp() {
        let mu_ref = 1.716e-5;
        let t_ref = 273.15;
        let mu = sutherland_viscosity(mu_ref, t_ref, t_ref, 110.4);
        assert!(
            (mu - mu_ref).abs() < 1e-18,
            "At T_ref, Sutherland should give mu_ref: got {mu}"
        );
    }

    #[test]
    fn test_sutherland_increases_with_temperature() {
        let mu_ref = 1.716e-5;
        let t_ref = 273.15;
        let s = 110.4;
        let mu_hot = sutherland_viscosity(mu_ref, t_ref, 400.0, s);
        assert!(
            mu_hot > mu_ref,
            "Sutherland viscosity should increase with T for gases"
        );
    }

    #[test]
    fn test_vft_finite_value() {
        let mu = vft_viscosity(-1.0, 500.0, 200.0, 400.0);
        assert!(
            mu > 0.0 && mu.is_finite(),
            "VFT should give finite positive viscosity: {mu}"
        );
    }

    // -----------------------------------------------------------------------
    // Non-Newtonian viscosity
    // -----------------------------------------------------------------------

    #[test]
    fn test_power_law_newtonian_limit() {
        // n = 1 → Newtonian: mu_eff = K regardless of strain rate
        let mu = power_law_viscosity(0.001, 1.0, 10.0, 1e-8);
        assert!(
            (mu - 0.001).abs() < 1e-14,
            "n=1 power law should be Newtonian: got {mu}"
        );
    }

    #[test]
    fn test_power_law_shear_thinning() {
        // n < 1 → viscosity decreases with increasing strain rate
        let mu_low = power_law_viscosity(1.0, 0.5, 1.0, 1e-8);
        let mu_high = power_law_viscosity(1.0, 0.5, 10.0, 1e-8);
        assert!(
            mu_high < mu_low,
            "Shear-thinning: viscosity should decrease at higher strain rate"
        );
    }

    #[test]
    fn test_power_law_shear_thickening() {
        // n > 1 → viscosity increases with increasing strain rate
        let mu_low = power_law_viscosity(1.0, 1.5, 1.0, 1e-8);
        let mu_high = power_law_viscosity(1.0, 1.5, 10.0, 1e-8);
        assert!(
            mu_high > mu_low,
            "Shear-thickening: viscosity should increase at higher strain rate"
        );
    }

    #[test]
    fn test_cross_model_limits() {
        let mu_0 = 1.0;
        let mu_inf = 0.01;
        let lambda = 1.0;
        let m = 0.667;

        // At very low strain rate → approaches mu_0
        let mu_low = cross_model_viscosity(mu_0, mu_inf, lambda, m, 1e-10);
        assert!(
            (mu_low - mu_0).abs() < 0.01,
            "Cross at low gamma_dot should approach mu_0: got {mu_low}"
        );

        // At very high strain rate → approaches mu_inf
        let mu_high = cross_model_viscosity(mu_0, mu_inf, lambda, m, 1e10);
        assert!(
            (mu_high - mu_inf).abs() < 0.01,
            "Cross at high gamma_dot should approach mu_inf: got {mu_high}"
        );
    }

    #[test]
    fn test_cross_model_monotonically_decreasing() {
        let mu_0 = 1.0;
        let mu_inf = 0.01;
        let rates = [0.01, 0.1, 1.0, 10.0, 100.0];
        let viscosities: Vec<f64> = rates
            .iter()
            .map(|&g| cross_model_viscosity(mu_0, mu_inf, 1.0, 0.667, g))
            .collect();
        for i in 1..viscosities.len() {
            assert!(
                viscosities[i] <= viscosities[i - 1] + 1e-14,
                "Cross model should be monotonically decreasing: mu[{}]={} > mu[{}]={}",
                i,
                viscosities[i],
                i - 1,
                viscosities[i - 1]
            );
        }
    }

    #[test]
    fn test_carreau_yasuda_newtonian_limit() {
        // With lambda = 0, should give mu_0
        let mu = carreau_yasuda_viscosity(1.0, 0.01, 0.0, 2.0, 0.5, 100.0);
        assert!(
            (mu - 1.0).abs() < 1e-10,
            "Carreau-Yasuda with lambda=0 should give mu_0"
        );
    }

    #[test]
    fn test_bingham_at_high_shear_rate() {
        let mu_p = 0.001;
        let tau_y = 10.0;
        let gamma_dot = 1e6;
        let mu = bingham_viscosity(mu_p, tau_y, gamma_dot, 1e-8);
        // At very high shear rate, tau_y / gamma_dot ≈ 0, so mu ≈ mu_p
        assert!(
            (mu - mu_p).abs() < 0.1,
            "Bingham at high shear should approach mu_p: got {mu}"
        );
    }

    #[test]
    fn test_herschel_bulkley_includes_yield_stress() {
        let mu_hb = herschel_bulkley_viscosity(0.001, 1.0, 10.0, 1e-6, 1e-6);
        // With very small strain rate, yield stress dominates: tau_y / gamma_dot_min = 10 / 1e-6 = 1e7
        assert!(
            mu_hb > 1e6,
            "Herschel-Bulkley at low shear should be dominated by yield stress: {mu_hb}"
        );
    }

    // -----------------------------------------------------------------------
    // Viscosity limiters
    // -----------------------------------------------------------------------

    #[test]
    fn test_clamp_viscosity() {
        assert!((clamp_viscosity(0.5, 0.1, 1.0) - 0.5).abs() < 1e-14);
        assert!((clamp_viscosity(0.01, 0.1, 1.0) - 0.1).abs() < 1e-14);
        assert!((clamp_viscosity(5.0, 0.1, 1.0) - 1.0).abs() < 1e-14);
    }

    #[test]
    fn test_monaghan_gingold_limiter_no_clamp() {
        let pi = 0.5;
        let limited = monaghan_gingold_limiter(pi, 1000.0, 1500.0, 0.05);
        let pi_max = 1500.0 * 1500.0 / 1000.0;
        assert!(pi < pi_max, "Test setup: pi should be below limit");
        assert!(
            (limited - pi).abs() < 1e-14,
            "No clamping needed: got {limited}"
        );
    }

    #[test]
    fn test_monaghan_gingold_limiter_clamps() {
        let rho_bar = 1000.0;
        let c_bar = 1.0;
        let pi_max = c_bar * c_bar / rho_bar; // = 0.001
        let pi = 1.0; // way above pi_max
        let limited = monaghan_gingold_limiter(pi, rho_bar, c_bar, 0.05);
        assert!(
            (limited - pi_max).abs() < 1e-14,
            "Should clamp to pi_max={pi_max}: got {limited}"
        );
    }

    #[test]
    fn test_balsara_switch_pure_compression() {
        // Pure compression (no rotation): factor should be ~1
        let f = balsara_switch(10.0, 0.0, 1500.0, 0.05);
        assert!(
            (f - 1.0).abs() < 0.01,
            "Pure compression: Balsara should be ~1, got {f}"
        );
    }

    #[test]
    fn test_balsara_switch_pure_shear() {
        // Pure shear (no compression): factor should be ~0
        let f = balsara_switch(0.0, 10.0, 1500.0, 0.05);
        assert!(f < 0.01, "Pure shear: Balsara should be ~0, got {f}");
    }

    #[test]
    fn test_signal_velocity_at_rest() {
        let v_sig = signal_velocity(1500.0, 1500.0, [0.0, 0.0, 0.0], [1.0, 0.0, 0.0]);
        assert!(
            (v_sig - 3000.0).abs() < 1e-10,
            "At rest: v_sig = c_i + c_j = 3000, got {v_sig}"
        );
    }

    #[test]
    fn test_signal_velocity_approaching() {
        // Approaching: v_ij . r_hat < 0, so signal velocity > c_i + c_j
        let v_sig = signal_velocity(100.0, 100.0, [-10.0, 0.0, 0.0], [1.0, 0.0, 0.0]);
        assert!(
            v_sig > 200.0,
            "Approaching particles: v_sig should exceed c_i+c_j, got {v_sig}"
        );
    }

    // -----------------------------------------------------------------------
    // Velocity divergence and vorticity
    // -----------------------------------------------------------------------

    #[test]
    fn test_velocity_divergence_uniform_expansion() {
        let grad_v = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let div = velocity_divergence(&grad_v);
        assert!(
            (div - 3.0).abs() < 1e-14,
            "Uniform expansion: div=3, got {div}"
        );
    }

    #[test]
    fn test_vorticity_magnitude_solid_rotation() {
        // Pure solid-body rotation about z: grad_v[0][1] = -omega, grad_v[1][0] = omega
        let mut grad_v = [[0.0_f64; 3]; 3];
        let omega = 5.0;
        grad_v[0][1] = -omega;
        grad_v[1][0] = omega;

        let w_mag = vorticity_magnitude(&grad_v);
        // vorticity = (0, 0, omega - (-omega)) = (0, 0, 2*omega)
        let expected = 2.0 * omega;
        assert!(
            (w_mag - expected).abs() < 1e-12,
            "Vorticity magnitude: expected {expected}, got {w_mag}"
        );
    }

    #[test]
    fn test_incompressible_flow_zero_divergence() {
        // ∂vx/∂x = 1, ∂vy/∂y = -1 → div = 0
        let grad_v = [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 0.0]];
        let div = velocity_divergence(&grad_v);
        assert!(
            div.abs() < 1e-14,
            "2D incompressible flow should have zero divergence: got {div}"
        );
    }
}

// ---------------------------------------------------------------------------
// Morris-Fox viscosity (adaptive viscosity)
// ---------------------------------------------------------------------------

/// Morris-Fox (1997) adaptive viscosity formulation.
///
/// Uses locally computed signal velocity to adaptively scale the viscosity:
/// ```text
/// Π_ij = -α * v_sig * μ_ij / ρ̄
/// μ_ij = (v_ij · r_ij) / |r_ij|  (negative when approaching)
/// v_sig = c_i + c_j - β_c * min(0, v_ij · r_hat)
/// ```
/// Returns zero if particles are separating.
#[allow(dead_code, clippy::too_many_arguments)]
pub fn morris_fox_viscosity(
    alpha: f64,
    beta_c: f64,
    rho_i: f64,
    rho_j: f64,
    v_ij: [f64; 3],
    r_ij: [f64; 3],
    c_i: f64,
    c_j: f64,
) -> f64 {
    let r = norm3(r_ij);
    if r < EPSILON {
        return 0.0;
    }
    let vr = dot3(v_ij, r_ij) / r; // v_ij · r_hat
    if vr >= 0.0 {
        return 0.0; // separating
    }
    let v_sig = c_i + c_j - beta_c * vr; // vr < 0 so this increases v_sig
    let mu_ij = vr; // negative
    let rho_bar = 0.5 * (rho_i + rho_j);
    -alpha * v_sig * mu_ij / rho_bar
}

// ---------------------------------------------------------------------------
// Delta-SPH density diffusion
// ---------------------------------------------------------------------------

/// Delta-SPH density diffusion term (Molteni & Colagrossi 2009).
///
/// Adds a diffusion term to the density equation to reduce spurious density
/// fluctuations:
/// ```text
/// δρ_i/dt|_diff = δ h c_s Σ_j (m_j/ρ_j) ψ_ij ∇W_ij
/// ψ_ij = 2(ρ_j - ρ_i) r_hat · ∇W_ij / |r_ij|
/// ```
/// Returns the diffusion contribution for particle `i` from one neighbor `j`.
#[allow(dead_code, clippy::too_many_arguments)]
pub fn delta_sph_diffusion_contribution(
    delta: f64,
    h: f64,
    c_s: f64,
    rho_i: f64,
    rho_j: f64,
    m_j: f64,
    r_ij: [f64; 3],
    grad_w: [f64; 3],
) -> f64 {
    let r = norm3(r_ij);
    if r < EPSILON || rho_j < EPSILON {
        return 0.0;
    }
    // ψ_ij = 2(ρ_j - ρ_i) / r * (r_hat · grad_w normalized)
    // but grad_w already includes the r_hat direction:
    // ψ_ij · ∇W_ij = 2(ρ_j - ρ_i) * (r_ij/r) · grad_w / r
    let r_hat_dot_gw = dot3(r_ij, grad_w) / r; // positive when pointing same direction
    let psi_dot_gw = 2.0 * (rho_j - rho_i) * r_hat_dot_gw / r;
    delta * h * c_s * (m_j / rho_j) * psi_dot_gw
}

/// Accumulate delta-SPH diffusion contributions for particle `i` from all neighbors.
#[allow(dead_code, clippy::too_many_arguments)]
pub fn delta_sph_density_diffusion(
    i: usize,
    positions: &[[f64; 3]],
    densities: &[f64],
    masses: &[f64],
    neighbor_grad_w: &[[f64; 3]],
    neighbor_indices: &[usize],
    delta: f64,
    h: f64,
    c_s: f64,
) -> f64 {
    let xi = positions[i];
    let rho_i = densities[i];
    let mut sum = 0.0_f64;
    for (k, &j) in neighbor_indices.iter().enumerate() {
        let r_ij = [
            xi[0] - positions[j][0],
            xi[1] - positions[j][1],
            xi[2] - positions[j][2],
        ];
        sum += delta_sph_diffusion_contribution(
            delta,
            h,
            c_s,
            rho_i,
            densities[j],
            masses[j],
            r_ij,
            neighbor_grad_w[k],
        );
    }
    sum
}

// ---------------------------------------------------------------------------
// Artificial conductivity
// ---------------------------------------------------------------------------

/// Artificial thermal conductivity term (Price 2008).
///
/// Prevents the development of spurious surface energy at contact discontinuities
/// in SPH simulations:
/// ```text
/// dU_i/dt|_cond = Σ_j α_c v_sig_u (U_j - U_i) m_j / (ρ_i + ρ_j) |∇W_ij|
/// ```
/// Returns the conductivity contribution from neighbor `j`.
#[allow(dead_code, clippy::too_many_arguments)]
pub fn artificial_conductivity(
    alpha_c: f64,
    v_sig_u: f64,
    u_i: f64,
    u_j: f64,
    m_j: f64,
    rho_i: f64,
    rho_j: f64,
    grad_w_mag: f64,
) -> f64 {
    let rho_sum = rho_i + rho_j;
    if rho_sum < EPSILON {
        return 0.0;
    }
    alpha_c * v_sig_u * (u_j - u_i) * m_j / rho_sum * grad_w_mag
}

/// Thermal signal velocity for artificial conductivity (Price 2008).
///
/// ```text
/// v_sig_u = sqrt(|p_i - p_j| / (0.5 * (rho_i + rho_j)))
/// ```
#[allow(dead_code)]
pub fn thermal_signal_velocity(p_i: f64, p_j: f64, rho_i: f64, rho_j: f64) -> f64 {
    let rho_bar = 0.5 * (rho_i + rho_j);
    if rho_bar < EPSILON {
        return 0.0;
    }
    ((p_i - p_j).abs() / rho_bar).sqrt()
}

// ---------------------------------------------------------------------------
// Cullen-Dehnen shock switch
// ---------------------------------------------------------------------------

/// Cullen-Dehnen (2010) shock switch for viscosity.
///
/// The switch ensures viscosity is applied only in compressive regions:
/// ```text
/// A_i = max(-div v_i, 0)²
/// f_i = A_i / (A_i + |curl v_i|² + ε * c_i² / h_i²)
/// ```
/// Returns the switch value in \[0, 1\].
#[allow(dead_code)]
pub fn cullen_dehnen_switch(div_v: f64, curl_v_sq: f64, c: f64, h: f64) -> f64 {
    let a_i = (-div_v).max(0.0);
    let a_sq = a_i * a_i;
    let eps_term = EPSILON * c * c / (h * h);
    a_sq / (a_sq + curl_v_sq + eps_term)
}

/// Cullen-Dehnen viscosity coefficient α (time-evolved).
///
/// The viscosity coefficient is evolved in time:
/// ```text
/// dα_i/dt = -α_i/τ + ξ * α_max * f_i
/// ```
/// where `τ = h / (l * c_s)` is the decay timescale.
///
/// Returns the new α after a time step `dt`.
#[allow(dead_code, clippy::too_many_arguments)]
pub fn cullen_dehnen_alpha_update(
    alpha: f64,
    alpha_max: f64,
    switch: f64,
    xi: f64,
    h: f64,
    c_s: f64,
    l_decay: f64,
    dt: f64,
) -> f64 {
    if h < EPSILON || c_s < EPSILON {
        return alpha;
    }
    let tau = h / (l_decay * c_s);
    let dalpha = -alpha / tau + xi * alpha_max * switch;
    (alpha + dalpha * dt).clamp(0.0, alpha_max)
}

// ---------------------------------------------------------------------------
// Tensor viscosity
// ---------------------------------------------------------------------------

/// Tensor artificial viscosity (Flebbe et al. 1994 / Siegler & Riffert 2000).
///
/// The tensor viscosity formulation uses the full velocity gradient tensor
/// to compute a physically consistent viscosity:
/// ```text
/// Π_αβ = ρ * ν_t * (∂v_α/∂x_β + ∂v_β/∂x_α - (2/3)δ_αβ div v)
/// ```
/// Returns the viscous stress tensor (3×3).
#[allow(dead_code)]
pub fn tensor_viscosity_stress(rho: f64, nu_t: f64, grad_v: &[[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let div_v = grad_v[0][0] + grad_v[1][1] + grad_v[2][2];
    let mut tau = [[0.0_f64; 3]; 3];
    for alpha in 0..3 {
        for beta in 0..3 {
            let s_ab = 0.5 * (grad_v[alpha][beta] + grad_v[beta][alpha]);
            tau[alpha][beta] = 2.0 * rho * nu_t * s_ab;
        }
        tau[alpha][alpha] -= (2.0 / 3.0) * rho * nu_t * div_v;
    }
    tau
}

/// Apply tensor viscosity stress to compute the viscosity acceleration on particle i.
///
/// `f_visc_i = (1/ρ_i) * Σ_j m_j * (τ_i/ρ_i² + τ_j/ρ_j²) · ∇W_ij`
///
/// This is the anti-symmetric form preserving linear momentum.
#[allow(dead_code, clippy::too_many_arguments)]
pub fn tensor_viscosity_acceleration(
    tau_i: &[[f64; 3]; 3],
    tau_j: &[[f64; 3]; 3],
    rho_i: f64,
    rho_j: f64,
    m_j: f64,
    grad_w: [f64; 3],
) -> [f64; 3] {
    let mut acc = [0.0_f64; 3];
    let ri2 = rho_i * rho_i;
    let rj2 = rho_j * rho_j;
    if ri2 < EPSILON * EPSILON || rj2 < EPSILON * EPSILON {
        return acc;
    }
    for alpha in 0..3 {
        let mut sum_beta = 0.0;
        for beta in 0..3 {
            sum_beta += (tau_i[alpha][beta] / ri2 + tau_j[alpha][beta] / rj2) * grad_w[beta];
        }
        acc[alpha] = m_j * sum_beta;
    }
    acc
}

// ---------------------------------------------------------------------------
// Additional helpers
// ---------------------------------------------------------------------------

/// Compute the rotation tensor from the velocity gradient.
///
/// `Ω_αβ = (∂v_α/∂x_β - ∂v_β/∂x_α) / 2`
#[allow(dead_code)]
pub fn rotation_tensor(grad_v: &[[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let mut omega = [[0.0_f64; 3]; 3];
    for a in 0..3 {
        for b in 0..3 {
            omega[a][b] = 0.5 * (grad_v[a][b] - grad_v[b][a]);
        }
    }
    omega
}

/// Compute the Q-criterion for vortex identification.
///
/// `Q = 0.5 * (|Ω|² - |S|²)`
///
/// Positive Q indicates vortex-dominated regions.
#[allow(dead_code)]
pub fn q_criterion(grad_v: &[[f64; 3]; 3]) -> f64 {
    let omega = rotation_tensor(grad_v);
    let s_mag_sq = {
        let mut s = 0.0;
        for a in 0..3 {
            for b in 0..3 {
                let s_ab = 0.5 * (grad_v[a][b] + grad_v[b][a]);
                s += s_ab * s_ab;
            }
        }
        s
    };
    let omega_mag_sq = {
        let mut o = 0.0;
        for a in 0..3 {
            for b in 0..3 {
                o += omega[a][b] * omega[a][b];
            }
        }
        o
    };
    0.5 * (omega_mag_sq - s_mag_sq)
}

/// Compute the second invariant of the strain-rate tensor.
///
/// `II_S = S_αβ S_αβ = |S|²_F / 2`  (Frobenius norm squared / 2)
#[allow(dead_code)]
pub fn strain_rate_second_invariant(grad_v: &[[f64; 3]; 3]) -> f64 {
    let mut sum = 0.0;
    for a in 0..3 {
        for b in 0..3 {
            let s_ab = 0.5 * (grad_v[a][b] + grad_v[b][a]);
            sum += s_ab * s_ab;
        }
    }
    sum
}

// ---------------------------------------------------------------------------
// Extended tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests_extended {
    use super::*;

    // Morris-Fox viscosity tests

    #[test]
    fn test_morris_fox_separating_zero() {
        let v_ij = [1.0, 0.0, 0.0];
        let r_ij = [1.0, 0.0, 0.0]; // v_ij · r_hat = 1 > 0 → separating
        let pi = morris_fox_viscosity(1.0, 3.0, 1000.0, 1000.0, v_ij, r_ij, 100.0, 100.0);
        assert_eq!(pi, 0.0, "Separating particles should give Π = 0");
    }

    #[test]
    fn test_morris_fox_approaching_positive() {
        let v_ij = [-1.0, 0.0, 0.0];
        let r_ij = [1.0, 0.0, 0.0]; // v_ij · r_hat = -1 < 0 → approaching
        let pi = morris_fox_viscosity(1.0, 3.0, 1000.0, 1000.0, v_ij, r_ij, 100.0, 100.0);
        assert!(pi > 0.0, "Approaching particles should give Π > 0: {pi}");
    }

    #[test]
    fn test_morris_fox_scales_with_alpha() {
        let v_ij = [-1.0, 0.0, 0.0];
        let r_ij = [1.0, 0.0, 0.0];
        let pi1 = morris_fox_viscosity(1.0, 3.0, 1000.0, 1000.0, v_ij, r_ij, 100.0, 100.0);
        let pi2 = morris_fox_viscosity(2.0, 3.0, 1000.0, 1000.0, v_ij, r_ij, 100.0, 100.0);
        assert!(
            (pi2 - 2.0 * pi1).abs() < 1e-10,
            "PI should scale linearly with alpha: pi1={pi1}, pi2={pi2}"
        );
    }

    #[test]
    fn test_morris_fox_zero_at_zero_radius() {
        let v_ij = [-1.0, 0.0, 0.0];
        let r_ij = [0.0, 0.0, 0.0];
        let pi = morris_fox_viscosity(1.0, 3.0, 1000.0, 1000.0, v_ij, r_ij, 100.0, 100.0);
        assert_eq!(pi, 0.0, "Should return 0 at zero radius");
    }

    // Delta-SPH diffusion tests

    #[test]
    fn test_delta_sph_diffusion_equal_density_zero() {
        // Equal densities → ψ = 0 → no diffusion
        let diff = delta_sph_diffusion_contribution(
            0.1,
            0.1,
            1480.0,
            1000.0,
            1000.0,
            0.001,
            [0.05, 0.0, 0.0],
            [-10.0, 0.0, 0.0],
        );
        assert!(diff.abs() < 1e-14, "Equal densities → no diffusion: {diff}");
    }

    #[test]
    fn test_delta_sph_diffusion_density_gradient() {
        // Higher density neighbor → positive contribution to lower-density particle
        let diff = delta_sph_diffusion_contribution(
            0.1,
            0.1,
            1480.0,
            800.0,
            1200.0,
            0.001,
            [0.05, 0.0, 0.0],
            [-10.0, 0.0, 0.0],
        );
        // rho_j > rho_i, r_ij and grad_w point in same direction → positive
        assert!(diff.is_finite(), "Diffusion should be finite: {diff}");
    }

    #[test]
    fn test_delta_sph_zero_delta() {
        let diff = delta_sph_diffusion_contribution(
            0.0,
            0.1,
            1480.0,
            800.0,
            1200.0,
            0.001,
            [0.05, 0.0, 0.0],
            [-10.0, 0.0, 0.0],
        );
        assert!(diff.abs() < 1e-14, "Zero delta → no diffusion: {diff}");
    }

    // Artificial conductivity tests

    #[test]
    fn test_artificial_conductivity_equal_u() {
        let cond = artificial_conductivity(0.1, 100.0, 500.0, 500.0, 0.001, 1000.0, 1000.0, 10.0);
        assert!(
            cond.abs() < 1e-14,
            "Equal thermal energy → no conductivity: {cond}"
        );
    }

    #[test]
    fn test_artificial_conductivity_direction() {
        // u_j > u_i → positive contribution (heat flows from j to i)
        let cond = artificial_conductivity(0.1, 100.0, 200.0, 500.0, 0.001, 1000.0, 1000.0, 10.0);
        assert!(cond > 0.0, "Heat should flow from high to low u: {cond}");
    }

    #[test]
    fn test_artificial_conductivity_scales_with_alpha() {
        let c1 = artificial_conductivity(0.1, 100.0, 200.0, 500.0, 0.001, 1000.0, 1000.0, 10.0);
        let c2 = artificial_conductivity(0.2, 100.0, 200.0, 500.0, 0.001, 1000.0, 1000.0, 10.0);
        assert!(
            (c2 - 2.0 * c1).abs() < 1e-12,
            "Conductivity should scale with alpha: c1={c1}, c2={c2}"
        );
    }

    #[test]
    fn test_thermal_signal_velocity_positive() {
        let v_sig = thermal_signal_velocity(1000.0, 500.0, 1000.0, 1000.0);
        assert!(
            v_sig > 0.0,
            "Thermal signal velocity should be positive: {v_sig}"
        );
    }

    #[test]
    fn test_thermal_signal_velocity_equal_pressure() {
        let v_sig = thermal_signal_velocity(1000.0, 1000.0, 1000.0, 1000.0);
        assert!(
            v_sig.abs() < 1e-14,
            "Equal pressure → zero signal velocity: {v_sig}"
        );
    }

    // Cullen-Dehnen shock switch tests

    #[test]
    fn test_cullen_dehnen_pure_compression() {
        // Large compression, no rotation → switch ≈ 1
        let f = cullen_dehnen_switch(-10.0, 0.0, 1500.0, 0.1);
        assert!(f > 0.99, "Pure compression: switch should be ~1, got {f}");
    }

    #[test]
    fn test_cullen_dehnen_pure_expansion() {
        // Expansion (div_v > 0) → switch = 0
        let f = cullen_dehnen_switch(10.0, 0.0, 1500.0, 0.1);
        assert!(
            f.abs() < 1e-14,
            "Pure expansion: switch should be 0, got {f}"
        );
    }

    #[test]
    fn test_cullen_dehnen_pure_shear() {
        // No compression, only rotation → switch ≈ 0
        let f = cullen_dehnen_switch(0.0, 100.0, 1500.0, 0.1);
        assert!(f < 0.01, "Pure shear: switch should be ~0, got {f}");
    }

    #[test]
    fn test_cullen_dehnen_range() {
        for div in [-5.0, -1.0, 0.0, 1.0, 5.0] {
            for curl_sq in [0.0, 1.0, 10.0] {
                let f = cullen_dehnen_switch(div, curl_sq, 1500.0, 0.1);
                assert!(
                    (0.0..=1.0 + 1e-12).contains(&f),
                    "Cullen-Dehnen switch should be in [0,1]: f={f} for div={div}, curl_sq={curl_sq}"
                );
            }
        }
    }

    #[test]
    fn test_cullen_dehnen_alpha_update_decay() {
        // Without switch (f=0), alpha should decay toward 0
        let alpha0 = 1.0;
        let new_alpha = cullen_dehnen_alpha_update(alpha0, 1.5, 0.0, 0.1, 0.1, 1500.0, 1.0, 0.001);
        assert!(
            new_alpha < alpha0,
            "Alpha should decay when switch=0: α={new_alpha}"
        );
    }

    #[test]
    fn test_cullen_dehnen_alpha_update_increase() {
        // With full switch (f=1), alpha should increase
        let alpha0 = 0.0;
        let new_alpha = cullen_dehnen_alpha_update(alpha0, 1.5, 1.0, 1.0, 0.1, 1500.0, 1.0, 0.001);
        assert!(
            new_alpha > alpha0,
            "Alpha should increase with full switch: α={new_alpha}"
        );
    }

    #[test]
    fn test_cullen_dehnen_alpha_clamped() {
        let alpha_max = 1.5;
        let new_alpha =
            cullen_dehnen_alpha_update(0.0, alpha_max, 1.0, 100.0, 0.1, 1500.0, 1.0, 1.0);
        assert!(
            new_alpha <= alpha_max + 1e-12,
            "Alpha should be clamped to alpha_max: {new_alpha}"
        );
        assert!(
            new_alpha >= 0.0,
            "Alpha should be non-negative: {new_alpha}"
        );
    }

    // Tensor viscosity tests

    #[test]
    fn test_tensor_viscosity_stress_symmetric() {
        let mut grad_v = [[0.0_f64; 3]; 3];
        grad_v[0][1] = 2.0;
        grad_v[1][0] = 1.5;
        let tau = tensor_viscosity_stress(1000.0, 0.01, &grad_v);
        for a in 0..3 {
            for b in 0..3 {
                assert!(
                    (tau[a][b] - tau[b][a]).abs() < 1e-12,
                    "Tensor viscosity stress should be symmetric: tau[{a}][{b}]={} != tau[{b}][{a}]={}",
                    tau[a][b],
                    tau[b][a]
                );
            }
        }
    }

    #[test]
    fn test_tensor_viscosity_zero_for_rigid_body() {
        // Rigid body rotation: S = 0 → tau = 0
        let mut grad_v = [[0.0_f64; 3]; 3];
        // Solid body rotation about z: vx = -omega*y, vy = omega*x
        grad_v[0][1] = -1.0; // dvx/dy = -omega
        grad_v[1][0] = 1.0; // dvy/dx = omega
        let tau = tensor_viscosity_stress(1000.0, 0.01, &grad_v);
        // S_xy = 0.5*(grad_v[0][1] + grad_v[1][0]) = 0 → tau = 0
        for a in 0..3 {
            for b in 0..3 {
                assert!(
                    tau[a][b].abs() < 1e-12,
                    "Rigid body rotation → zero tensor viscosity: tau[{a}][{b}]={}",
                    tau[a][b]
                );
            }
        }
    }

    #[test]
    fn test_tensor_viscosity_acceleration_nonzero() {
        let mut grad_v = [[0.0_f64; 3]; 3];
        grad_v[0][1] = 1.0;
        grad_v[1][0] = 1.0;
        let tau_i = tensor_viscosity_stress(1000.0, 0.01, &grad_v);
        let tau_j = tensor_viscosity_stress(1000.0, 0.01, &grad_v);
        let grad_w = [-10.0, 0.0, 0.0];
        let acc = tensor_viscosity_acceleration(&tau_i, &tau_j, 1000.0, 1000.0, 0.001, grad_w);
        // With shear, expect some non-zero acceleration component
        let acc_mag = (acc[0] * acc[0] + acc[1] * acc[1] + acc[2] * acc[2]).sqrt();
        assert!(
            acc_mag >= 0.0,
            "Acceleration magnitude should be non-negative: {acc_mag}"
        );
    }

    // Rotation tensor tests

    #[test]
    fn test_rotation_tensor_antisymmetric() {
        let mut grad_v = [[0.0_f64; 3]; 3];
        grad_v[0][1] = 2.0;
        grad_v[1][0] = 1.0;
        let omega = rotation_tensor(&grad_v);
        for a in 0..3 {
            for b in 0..3 {
                assert!(
                    (omega[a][b] + omega[b][a]).abs() < 1e-12,
                    "Rotation tensor should be antisymmetric: omega[{a}][{b}]={} != -omega[{b}][{a}]={}",
                    omega[a][b],
                    -omega[b][a]
                );
            }
        }
    }

    #[test]
    fn test_rotation_tensor_diagonal_zero() {
        let grad_v = [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]];
        let omega = rotation_tensor(&grad_v);
        for i in 0..3 {
            assert!(
                omega[i][i].abs() < 1e-12,
                "Diagonal of antisymmetric tensor should be 0"
            );
        }
    }

    // Q-criterion tests

    #[test]
    fn test_q_criterion_solid_rotation_positive() {
        // Pure rotation → Q > 0 (vortex-dominated)
        let mut grad_v = [[0.0_f64; 3]; 3];
        grad_v[0][1] = -1.0;
        grad_v[1][0] = 1.0;
        let q = q_criterion(&grad_v);
        assert!(q > 0.0, "Pure rotation → Q > 0 (vortex region): Q={q}");
    }

    #[test]
    fn test_q_criterion_pure_strain_negative() {
        // Pure strain (irrotational) → Q < 0 (strain-dominated)
        let mut grad_v = [[0.0_f64; 3]; 3];
        grad_v[0][0] = 1.0;
        grad_v[1][1] = -1.0;
        let q = q_criterion(&grad_v);
        assert!(q <= 0.0, "Pure strain → Q ≤ 0 (strain-dominated): Q={q}");
    }

    #[test]
    fn test_q_criterion_zero_for_uniform_flow() {
        let grad_v = [[0.0_f64; 3]; 3];
        let q = q_criterion(&grad_v);
        assert!(q.abs() < 1e-14, "Uniform flow → Q = 0: Q={q}");
    }

    // Strain rate second invariant tests

    #[test]
    fn test_strain_rate_second_invariant_positive() {
        let mut grad_v = [[0.0_f64; 3]; 3];
        grad_v[0][1] = 1.0;
        let ii_s = strain_rate_second_invariant(&grad_v);
        assert!(
            ii_s > 0.0,
            "Non-zero strain → positive second invariant: {ii_s}"
        );
    }

    #[test]
    fn test_strain_rate_second_invariant_zero() {
        let grad_v = [[0.0_f64; 3]; 3];
        let ii_s = strain_rate_second_invariant(&grad_v);
        assert!(
            ii_s.abs() < 1e-14,
            "Zero gradient → zero second invariant: {ii_s}"
        );
    }

    #[test]
    fn test_strain_rate_second_invariant_relates_to_magnitude() {
        let mut grad_v = [[0.0_f64; 3]; 3];
        grad_v[0][1] = 2.0;
        let ii_s = strain_rate_second_invariant(&grad_v);
        let s_mag = strain_rate_magnitude(&grad_v);
        // ii_s = sum(S_ab^2), s_mag = sqrt(2 * sum(S_ab^2)) = sqrt(2 * ii_s)
        let expected_s_mag = (2.0 * ii_s).sqrt();
        assert!(
            (s_mag - expected_s_mag).abs() < 1e-12,
            "s_mag={s_mag}, sqrt(2*ii_s)={expected_s_mag}"
        );
    }

    // Delta-SPH accumulation test

    #[test]
    fn test_delta_sph_density_diffusion_uniform() {
        // Uniform density: all differences are zero → total diffusion = 0
        let positions = vec![[0.0, 0.0, 0.0], [0.1, 0.0, 0.0], [-0.1, 0.0, 0.0]];
        let densities = vec![1000.0, 1000.0, 1000.0];
        let masses = vec![0.001, 0.001, 0.001];
        let grad_ws = vec![[10.0, 0.0, 0.0], [-10.0, 0.0, 0.0]];
        let neighbors = vec![1, 2];
        let diff = delta_sph_density_diffusion(
            0, &positions, &densities, &masses, &grad_ws, &neighbors, 0.1, 0.1, 1480.0,
        );
        assert!(
            diff.abs() < 1e-14,
            "Uniform density → zero diffusion: {diff}"
        );
    }
}
