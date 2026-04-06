#![allow(clippy::ptr_arg)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Extended multiphase SPH.
//!
//! Covers: color function method, interface tracking (level-set + SPH),
//! droplet coalescence/breakup, Rayleigh-Taylor instability,
//! Kelvin-Helmholtz instability, surface tension (CSF, pairwise),
//! contact angle (Young-Laplace), immiscible fluid mixing, phase inversion,
//! droplet impact, capillary rise, thin film dynamics, bubble dynamics
//! (Rayleigh-Plesset), and multi-component SPH with diffusion.

#![allow(dead_code)]
#![allow(clippy::too_many_arguments)]

use std::f64::consts::PI;

// ============================================================================
// § 1  COLOR FUNCTION PARTICLE
// ============================================================================

/// Phase label for multiphase SPH.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum PhaseLabel {
    /// Liquid (e.g. water)
    Liquid,
    /// Gas (e.g. air/vapour)
    Gas,
    /// Solid boundary
    Solid,
}

/// A single SPH particle with a color function field for multiphase flow.
#[derive(Debug, Clone)]
pub struct ColorParticleExt {
    /// Position \[m\]
    pub pos: [f64; 3],
    /// Velocity \[m/s\]
    pub vel: [f64; 3],
    /// Color function value C ∈ \[0, 1\] (1 = liquid, 0 = gas)
    pub color: f64,
    /// Mass \[kg\]
    pub mass: f64,
    /// SPH density \[kg/m³\]
    pub density: f64,
    /// Pressure \[Pa\]
    pub pressure: f64,
    /// Smoothed color gradient ∇C
    pub color_grad: [f64; 3],
    /// Interface normal (unit)
    pub normal: [f64; 3],
    /// Interface curvature κ \[1/m\]
    pub curvature: f64,
    /// Phase label
    pub phase: PhaseLabel,
    /// Level-set signed distance φ
    pub phi: f64,
    /// Volume fraction α ∈ \[0, 1\]
    pub alpha: f64,
}

impl ColorParticleExt {
    /// Create a liquid particle at position `(x, y, z)` with given mass.
    pub fn new_liquid(x: f64, y: f64, z: f64, mass: f64) -> Self {
        Self {
            pos: [x, y, z],
            vel: [0.0; 3],
            color: 1.0,
            mass,
            density: 1000.0,
            pressure: 0.0,
            color_grad: [0.0; 3],
            normal: [0.0, 1.0, 0.0],
            curvature: 0.0,
            phase: PhaseLabel::Liquid,
            phi: -1.0,
            alpha: 1.0,
        }
    }

    /// Create a gas particle at position `(x, y, z)` with given mass.
    pub fn new_gas(x: f64, y: f64, z: f64, mass: f64) -> Self {
        Self {
            pos: [x, y, z],
            vel: [0.0; 3],
            color: 0.0,
            mass,
            density: 1.225,
            pressure: 0.0,
            color_grad: [0.0; 3],
            normal: [0.0, -1.0, 0.0],
            curvature: 0.0,
            phase: PhaseLabel::Gas,
            phi: 1.0,
            alpha: 0.0,
        }
    }

    /// Returns true if the particle is near the interface (0.1 < C < 0.9).
    #[inline]
    pub fn is_interface(&self) -> bool {
        self.color > 0.1 && self.color < 0.9
    }

    /// Compute interpolated density from color value.
    pub fn interpolated_density(&self, rho_l: f64, rho_g: f64) -> f64 {
        self.color * rho_l + (1.0 - self.color) * rho_g
    }

    /// Compute interpolated dynamic viscosity from color value.
    pub fn interpolated_viscosity(&self, mu_l: f64, mu_g: f64) -> f64 {
        self.color * mu_l + (1.0 - self.color) * mu_g
    }

    /// Update the interface normal from the color gradient.
    pub fn update_normal(&mut self) {
        let g = &self.color_grad;
        let norm = (g[0] * g[0] + g[1] * g[1] + g[2] * g[2]).sqrt();
        if norm > 1e-14 {
            self.normal[0] = g[0] / norm;
            self.normal[1] = g[1] / norm;
            self.normal[2] = g[2] / norm;
        }
    }
}

// ============================================================================
// § 2  KERNEL FUNCTIONS
// ============================================================================

/// Cubic spline kernel W(r, h) in 3-D.
///
/// Returns the kernel value at distance `r` for smoothing length `h`.
pub fn cubic_spline_3d(r: f64, h: f64) -> f64 {
    let q = r / h;
    let alpha = 1.0 / (PI * h * h * h);
    if q < 1.0 {
        alpha * (1.0 - 1.5 * q * q + 0.75 * q * q * q)
    } else if q < 2.0 {
        alpha * 0.25 * (2.0 - q).powi(3)
    } else {
        0.0
    }
}

/// Gradient of the cubic spline kernel dW/dr · (r̂) in 3-D.
///
/// Returns the gradient magnitude (multiply by r̂ to get vector gradient).
pub fn cubic_spline_gradient_3d(r: f64, h: f64) -> f64 {
    let q = r / h;
    let alpha = 1.0 / (PI * h * h * h);
    if r < 1e-14 {
        return 0.0;
    }
    if q < 1.0 {
        alpha / h * (-3.0 * q + 2.25 * q * q)
    } else if q < 2.0 {
        alpha / h * (-0.75 * (2.0 - q).powi(2))
    } else {
        0.0
    }
}

/// Wendland C2 kernel value in 3-D.
///
/// W(r, h) = α (1 − q/2)⁴ (2q + 1),  q = r/h,  α = 21/(2πh³)
pub fn wendland_c2_3d(r: f64, h: f64) -> f64 {
    let q = r / h;
    if q >= 2.0 {
        return 0.0;
    }
    let alpha = 21.0 / (2.0 * PI * h * h * h);
    let t = 1.0 - 0.5 * q;
    alpha * t.powi(4) * (2.0 * q + 1.0)
}

/// Gradient of the Wendland C2 kernel dW/dr (scalar magnitude).
pub fn wendland_c2_gradient_3d(r: f64, h: f64) -> f64 {
    let q = r / h;
    if q >= 2.0 || r < 1e-14 {
        return 0.0;
    }
    let alpha = 21.0 / (2.0 * PI * h * h * h);
    let t = 1.0 - 0.5 * q;
    // dW/dr = alpha * [ -2q t^3 + t^4 * (-5 + 5q) / 2h ] – simplified form
    alpha / h * (-5.0 * q * t.powi(3))
}

// ============================================================================
// § 3  COLOR FUNCTION TRANSPORT
// ============================================================================

/// Advect the color function for particle `i` using an SPH sum.
///
/// DC/Dt = −C ∇·u  (approximate; assumes weakly compressible)
///
/// Returns the color rate of change dC/dt for particle i.
pub fn color_advection_rate(
    i: usize,
    positions: &[[f64; 3]],
    velocities: &[[f64; 3]],
    colors: &[f64],
    masses: &[f64],
    densities: &[f64],
    h: f64,
) -> f64 {
    let n = positions.len();
    let mut dc_dt = 0.0;
    for j in 0..n {
        if i == j {
            continue;
        }
        let dx = [
            positions[i][0] - positions[j][0],
            positions[i][1] - positions[j][1],
            positions[i][2] - positions[j][2],
        ];
        let r = (dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]).sqrt();
        if r >= 2.0 * h || r < 1e-14 {
            continue;
        }
        let dw = cubic_spline_gradient_3d(r, h);
        let dv = [
            velocities[j][0] - velocities[i][0],
            velocities[j][1] - velocities[i][1],
            velocities[j][2] - velocities[i][2],
        ];
        let dv_dot_dx = dv[0] * dx[0] + dv[1] * dx[1] + dv[2] * dx[2];
        let dc = colors[j] - colors[i];
        dc_dt += masses[j] / densities[j].max(1e-30) * dc * dw / r * dv_dot_dx;
    }
    dc_dt
}

/// Compute color function gradients for all particles via SPH interpolation.
///
/// ∇Cᵢ = Σⱼ (mⱼ/ρⱼ) (Cⱼ − Cᵢ) ∇W_ij
pub fn compute_color_gradients_3d(particles: &mut Vec<ColorParticleExt>, h: f64) {
    let n = particles.len();
    let pos: Vec<[f64; 3]> = particles.iter().map(|p| p.pos).collect();
    let colors: Vec<f64> = particles.iter().map(|p| p.color).collect();
    let masses: Vec<f64> = particles.iter().map(|p| p.mass).collect();
    let dens: Vec<f64> = particles.iter().map(|p| p.density).collect();
    let mut grads = vec![[0.0_f64; 3]; n];
    for i in 0..n {
        for j in 0..n {
            if i == j {
                continue;
            }
            let dx = [
                pos[i][0] - pos[j][0],
                pos[i][1] - pos[j][1],
                pos[i][2] - pos[j][2],
            ];
            let r = (dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]).sqrt();
            if r >= 2.0 * h || r < 1e-14 {
                continue;
            }
            let dw = cubic_spline_gradient_3d(r, h);
            let factor = masses[j] / dens[j].max(1e-30) * (colors[j] - colors[i]) * dw / r;
            grads[i][0] += factor * dx[0];
            grads[i][1] += factor * dx[1];
            grads[i][2] += factor * dx[2];
        }
    }
    for (p, g) in particles.iter_mut().zip(grads.iter()) {
        p.color_grad = *g;
        p.update_normal();
    }
}

// ============================================================================
// § 4  CONTINUUM SURFACE FORCE (CSF) MODEL
// ============================================================================

/// Compute CSF surface tension force on particle i.
///
/// f_st = σ κᵢ nᵢ δ_s
///
/// where δ_s = |∇Cᵢ| is the interface delta function approximation.
pub fn csf_surface_tension_force(
    normal: &[f64; 3],
    curvature: f64,
    color_grad_magnitude: f64,
    sigma: f64,
) -> [f64; 3] {
    let factor = sigma * curvature * color_grad_magnitude;
    [factor * normal[0], factor * normal[1], factor * normal[2]]
}

/// Compute interface curvature from the divergence of the unit normal.
///
/// κᵢ = −∇ · n̂ᵢ ≈ −Σⱼ (mⱼ/ρⱼ) nⱼ · ∇Wᵢⱼ
pub fn compute_curvature(
    i: usize,
    positions: &[[f64; 3]],
    normals: &[[f64; 3]],
    masses: &[f64],
    densities: &[f64],
    h: f64,
) -> f64 {
    let n = positions.len();
    let mut kappa = 0.0;
    for j in 0..n {
        if i == j {
            continue;
        }
        let dx = [
            positions[i][0] - positions[j][0],
            positions[i][1] - positions[j][1],
            positions[i][2] - positions[j][2],
        ];
        let r = (dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]).sqrt();
        if r >= 2.0 * h || r < 1e-14 {
            continue;
        }
        let dw = cubic_spline_gradient_3d(r, h);
        let dn = [
            normals[j][0] - normals[i][0],
            normals[j][1] - normals[i][1],
            normals[j][2] - normals[i][2],
        ];
        let dn_dot_dx = dn[0] * dx[0] + dn[1] * dx[1] + dn[2] * dx[2];
        kappa -= masses[j] / densities[j].max(1e-30) * dw / r * dn_dot_dx;
    }
    kappa
}

// ============================================================================
// § 5  PAIRWISE SURFACE TENSION FORCE
// ============================================================================

/// Pairwise surface tension force between particles i and j.
///
/// f_ij = σ (Cᵢ − Cⱼ) ∇W_ij   (Morris 2000 / Tartakovsky & Meakin)
///
/// Returns the force vector on particle i due to particle j.
pub fn pairwise_surface_tension(
    dx: &[f64; 3],
    r: f64,
    ci: f64,
    cj: f64,
    mj: f64,
    rhoj: f64,
    sigma: f64,
    h: f64,
) -> [f64; 3] {
    if r < 1e-14 {
        return [0.0; 3];
    }
    let dw = cubic_spline_gradient_3d(r, h);
    let factor = sigma * (ci - cj) * mj / rhoj.max(1e-30) * dw / r;
    [factor * dx[0], factor * dx[1], factor * dx[2]]
}

// ============================================================================
// § 6  CONTACT ANGLE (YOUNG-LAPLACE)
// ============================================================================

/// Contact angle boundary condition adjustment.
///
/// Rotates the interface normal at a solid boundary to enforce the specified
/// contact angle θ (in radians) using the Young-Laplace relation.
///
/// Returns the modified normal in 2-D (x, y components).
pub fn contact_angle_normal_2d(n_fluid: [f64; 2], n_wall: [f64; 2], theta: f64) -> [f64; 2] {
    // n_contact = cos(θ) n_wall + sin(θ) t_wall
    // where t_wall is the tangent in the plane perpendicular to n_wall
    let t_wall = [-n_wall[1], n_wall[0]];
    let ct = theta.cos();
    let st = theta.sin();
    // Decide sign of tangent based on fluid normal projection
    let dot = n_fluid[0] * t_wall[0] + n_fluid[1] * t_wall[1];
    let sgn = if dot >= 0.0 { 1.0 } else { -1.0 };
    let nx = ct * n_wall[0] + sgn * st * t_wall[0];
    let ny = ct * n_wall[1] + sgn * st * t_wall[1];
    [nx, ny]
}

/// Young-Dupré equation: spreading coefficient S = σ_sg − σ_sl − σ_lg.
///
/// Positive S implies complete wetting.
pub fn young_dupre_spreading(sigma_sg: f64, sigma_sl: f64, sigma_lg: f64) -> f64 {
    sigma_sg - sigma_sl - sigma_lg
}

/// Young's equation: equilibrium contact angle from surface tensions.
///
/// cos θ = (σ_sg − σ_sl) / σ_lg
pub fn young_contact_angle(sigma_sg: f64, sigma_sl: f64, sigma_lg: f64) -> f64 {
    let cos_theta = (sigma_sg - sigma_sl) / sigma_lg.max(1e-30);
    cos_theta.clamp(-1.0, 1.0).acos()
}

// ============================================================================
// § 7  RAYLEIGH-TAYLOR INSTABILITY
// ============================================================================

/// Linear growth rate of the Rayleigh-Taylor instability.
///
/// σ_RT = sqrt( A_t · g · k )  where A_t = Atwood number, k = wavenumber.
pub fn rayleigh_taylor_growth_rate(atwood: f64, g: f64, k: f64) -> f64 {
    (atwood * g * k).max(0.0).sqrt()
}

/// Most unstable wavenumber for Rayleigh-Taylor with surface tension.
///
/// k* = sqrt(ρ_h g A_t / σ) · (1/2)^(1/2)
/// (critical wavenumber where surface tension stabilises growth).
pub fn rayleigh_taylor_critical_wavenumber(atwood: f64, rho_h: f64, g: f64, sigma: f64) -> f64 {
    if sigma < 1e-30 {
        return f64::INFINITY;
    }
    (rho_h * g * atwood / sigma).max(0.0).sqrt()
}

/// Atwood number A_t = (ρ_h − ρ_l) / (ρ_h + ρ_l).
pub fn atwood_number(rho_heavy: f64, rho_light: f64) -> f64 {
    let denom = rho_heavy + rho_light;
    if denom < 1e-30 {
        return 0.0;
    }
    (rho_heavy - rho_light) / denom
}

// ============================================================================
// § 8  KELVIN-HELMHOLTZ INSTABILITY
// ============================================================================

/// Linear growth rate for the Kelvin-Helmholtz instability (inviscid, 2-D).
///
/// σ_KH = k · |ΔU| · sqrt(ρ₁ρ₂) / (ρ₁ + ρ₂)   (ignoring surface tension)
pub fn kelvin_helmholtz_growth_rate(k: f64, delta_u: f64, rho1: f64, rho2: f64) -> f64 {
    let denom = rho1 + rho2;
    if denom < 1e-30 {
        return 0.0;
    }
    k * delta_u.abs() * (rho1 * rho2).sqrt() / denom
}

/// Critical velocity difference for KH stability with surface tension.
///
/// |ΔU|_crit = sqrt( 2σ (ρ₁ + ρ₂) k / (ρ₁ ρ₂) )
pub fn kelvin_helmholtz_critical_velocity(k: f64, rho1: f64, rho2: f64, sigma: f64) -> f64 {
    let denom = rho1 * rho2;
    if denom < 1e-30 {
        return f64::INFINITY;
    }
    (2.0 * sigma * (rho1 + rho2) * k / denom).max(0.0).sqrt()
}

// ============================================================================
// § 9  DROPLET DYNAMICS
// ============================================================================

/// Equivalent radius of a droplet from its volume.
///
/// r = (3V / 4π)^(1/3)
pub fn droplet_radius_from_volume(volume: f64) -> f64 {
    (3.0 * volume / (4.0 * PI)).powf(1.0 / 3.0)
}

/// Volume of a spherical droplet from its radius.
pub fn droplet_volume(radius: f64) -> f64 {
    4.0 / 3.0 * PI * radius.powi(3)
}

/// Weber number  We = ρ U² D / σ.
///
/// Characterises droplet deformation and breakup.
pub fn weber_number(rho: f64, u: f64, diameter: f64, sigma: f64) -> f64 {
    rho * u * u * diameter / sigma.max(1e-30)
}

/// Ohnesorge number  Oh = μ / sqrt(ρ σ D).
///
/// Characterises the balance of viscosity, surface tension, and inertia.
pub fn ohnesorge_number(mu: f64, rho: f64, sigma: f64, diameter: f64) -> f64 {
    mu / (rho * sigma * diameter).max(1e-30).sqrt()
}

/// Breakup criterion: a droplet breaks up when We > We_crit (typically 12).
pub fn droplet_breaks_up(we: f64, we_crit: f64) -> bool {
    we > we_crit
}

/// Simple coalescence check: two droplets coalesce when their surfaces overlap.
///
/// Returns true when the distance between centres is less than the sum of
/// radii minus a coalescence tolerance `tol`.
pub fn droplets_coalesce(pos_a: &[f64; 3], pos_b: &[f64; 3], r_a: f64, r_b: f64, tol: f64) -> bool {
    let dx = pos_a[0] - pos_b[0];
    let dy = pos_a[1] - pos_b[1];
    let dz = pos_a[2] - pos_b[2];
    let dist = (dx * dx + dy * dy + dz * dz).sqrt();
    dist < r_a + r_b - tol
}

// ============================================================================
// § 10  RAYLEIGH-PLESSET BUBBLE DYNAMICS
// ============================================================================

/// Rayleigh-Plesset equation RHS.
///
/// R̈ = \[ (p_b − p_∞) / ρ − 2σ/(ρR) − 4μ Ṙ/R \] / R  −  (3/2) Ṙ²/R
///
/// where `r_dot` = dR/dt and `r` = bubble radius.
pub fn rayleigh_plesset_rhs(
    r: f64,
    r_dot: f64,
    p_bubble: f64,
    p_inf: f64,
    rho: f64,
    sigma: f64,
    mu: f64,
) -> f64 {
    if r < 1e-14 {
        return 0.0;
    }
    let pressure_term = (p_bubble - p_inf) / rho;
    let surface_tension_term = 2.0 * sigma / (rho * r);
    let viscous_term = 4.0 * mu * r_dot / (rho * r);
    let kinetic_term = 1.5 * r_dot * r_dot / r;
    (pressure_term - surface_tension_term - viscous_term) / r - kinetic_term
}

/// Integrate the Rayleigh-Plesset equation with simple Euler steps.
///
/// Returns (radius, radius_dot) at time `t + dt`.
pub fn rayleigh_plesset_step(
    r: f64,
    r_dot: f64,
    p_bubble: f64,
    p_inf: f64,
    rho: f64,
    sigma: f64,
    mu: f64,
    dt: f64,
) -> (f64, f64) {
    let r_ddot = rayleigh_plesset_rhs(r, r_dot, p_bubble, p_inf, rho, sigma, mu);
    let r_new = (r + dt * r_dot).max(1e-14);
    let r_dot_new = r_dot + dt * r_ddot;
    (r_new, r_dot_new)
}

/// Minnaert resonance frequency of a spherical bubble.
///
/// f₀ = (1/(2πR₀)) · sqrt( 3γ p₀ / ρ )
pub fn minnaert_frequency(r0: f64, p0: f64, rho: f64, gamma: f64) -> f64 {
    (1.0 / (2.0 * PI * r0)) * (3.0 * gamma * p0 / rho).max(0.0).sqrt()
}

// ============================================================================
// § 11  CAPILLARY RISE
// ============================================================================

/// Jurin's law: equilibrium capillary rise height.
///
/// h = 2σ cos θ / (ρ g r)
pub fn jurin_height(sigma: f64, theta: f64, rho: f64, g: f64, r: f64) -> f64 {
    2.0 * sigma * theta.cos() / (rho * g * r).max(1e-30)
}

/// Lucas-Washburn equation: capillary rise height at time t.
///
/// h(t) = sqrt( r σ cos θ t / (2 μ) )
pub fn lucas_washburn_height(sigma: f64, theta: f64, r: f64, mu: f64, t: f64) -> f64 {
    (r * sigma * theta.cos() * t / (2.0 * mu).max(1e-30))
        .max(0.0)
        .sqrt()
}

/// Capillary number Ca = μ U / σ.
///
/// Characterises the relative importance of viscous to surface tension forces.
pub fn capillary_number(mu: f64, u: f64, sigma: f64) -> f64 {
    mu * u / sigma.max(1e-30)
}

/// Bond number Bo = ρ g l² / σ.
///
/// Characterises gravity vs. surface tension at length scale l.
pub fn bond_number(rho: f64, g: f64, l: f64, sigma: f64) -> f64 {
    rho * g * l * l / sigma.max(1e-30)
}

// ============================================================================
// § 12  THIN FILM DYNAMICS
// ============================================================================

/// Thin film (lubrication) equation: rate of change of film thickness h.
///
/// ∂h/∂t = −∂/∂x (h³/3μ · ∂p/∂x)  (1-D, pressure driven)
///
/// Returns dh/dt at node i given film thickness array and pressure gradient.
pub fn thin_film_dhdt(i: usize, thickness: &[f64], pressure_grad: &[f64], mu: f64, dx: f64) -> f64 {
    let n = thickness.len();
    if i == 0 || i >= n - 1 {
        return 0.0;
    }
    let h = thickness[i];
    let dp_dx = pressure_grad[i];
    let flux_r = thickness[i + 1].powi(3) / (3.0 * mu) * (pressure_grad[i + 1] + dp_dx) * 0.5;
    let flux_l = thickness[i - 1].powi(3) / (3.0 * mu) * (dp_dx + pressure_grad[i - 1]) * 0.5;
    let _ = h; // used for context
    -(flux_r - flux_l) / dx
}

/// Disjoining pressure for thin film (van der Waals + electrostatic).
///
/// Π(h) = A_H / (6π h³) − B exp(−h / λ_D)
pub fn disjoining_pressure(h: f64, hamaker: f64, b: f64, lambda_d: f64) -> f64 {
    let vdw = hamaker / (6.0 * PI * h.max(1e-12).powi(3));
    let edl = b * (-h / lambda_d.max(1e-12)).exp();
    vdw - edl
}

// ============================================================================
// § 13  PHASE INVERSION
// ============================================================================

/// Ambrosini phase inversion criterion.
///
/// Inversion is predicted when the dispersed phase volume fraction exceeds
/// a critical value φ_inv ≈ 1 / (1 + (μ_d/μ_c)^(1/3)).
pub fn phase_inversion_fraction(mu_dispersed: f64, mu_continuous: f64) -> f64 {
    let r = (mu_dispersed / mu_continuous.max(1e-30)).powf(1.0 / 3.0);
    1.0 / (1.0 + r)
}

/// Check whether phase inversion has occurred given current volume fraction.
pub fn phase_inversion_occurred(
    alpha_dispersed: f64,
    mu_dispersed: f64,
    mu_continuous: f64,
) -> bool {
    alpha_dispersed > phase_inversion_fraction(mu_dispersed, mu_continuous)
}

// ============================================================================
// § 14  DROPLET IMPACT
// ============================================================================

/// Maximum spreading diameter of a droplet after impact.
///
/// D_max / D_0 ≈ sqrt(We / 3 + 12) / sqrt(3(1 − cos θ) + 4 We / sqrt(Re))
///
/// (Pasandideh-Fard et al. simplified form for low Oh)
pub fn droplet_max_spread(we: f64, re: f64, theta: f64) -> f64 {
    let num = (we + 12.0).max(0.0).sqrt();
    let denom = (3.0 * (1.0 - theta.cos()) + 4.0 * we / re.max(1e-6).sqrt())
        .max(1e-30)
        .sqrt();
    num / denom
}

/// Reynolds number Re = ρ U D / μ.
pub fn reynolds_number(rho: f64, u: f64, diameter: f64, mu: f64) -> f64 {
    rho * u * diameter / mu.max(1e-30)
}

/// Splashing threshold We_splash based on Mundo criterion.
///
/// Splashing occurs when We^0.5 · Re^0.25 > K_splash ≈ 57.7
pub fn droplet_splashes(we: f64, re: f64, k_splash: f64) -> bool {
    we.sqrt() * re.powf(0.25) > k_splash
}

// ============================================================================
// § 15  IMMISCIBLE FLUID MIXING AND MULTI-COMPONENT DIFFUSION
// ============================================================================

/// Fick's law SPH diffusion rate for component k at particle i.
///
/// DC_k/Dt = Σⱼ D_ij (C_kj − C_ki) / r_ij · (2 mⱼ/ρⱼ) dW/dr
pub fn sph_diffusion_rate(
    i: usize,
    positions: &[[f64; 3]],
    concentrations: &[f64],
    masses: &[f64],
    densities: &[f64],
    diffusivity: f64,
    h: f64,
) -> f64 {
    let n = positions.len();
    let mut dc_dt = 0.0;
    for j in 0..n {
        if i == j {
            continue;
        }
        let dx = [
            positions[i][0] - positions[j][0],
            positions[i][1] - positions[j][1],
            positions[i][2] - positions[j][2],
        ];
        let r = (dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]).sqrt();
        if r >= 2.0 * h || r < 1e-14 {
            continue;
        }
        let dw = cubic_spline_gradient_3d(r, h);
        let dc = concentrations[j] - concentrations[i];
        dc_dt += diffusivity * masses[j] / densities[j].max(1e-30) * dc / r * dw * 2.0;
    }
    dc_dt
}

/// Effective diffusivity between two phases (harmonic mean).
pub fn harmonic_diffusivity(d1: f64, d2: f64, alpha1: f64) -> f64 {
    let alpha2 = 1.0 - alpha1;
    2.0 * d1 * d2 / (alpha1 * d2 + alpha2 * d1 + 1e-30)
}

// ============================================================================
// § 16  LEVEL-SET INTERFACE TRACKING (SPH-LS COUPLING)
// ============================================================================

/// Compute the signed distance level-set field from a color function.
///
/// φᵢ = (2Cᵢ − 1) · h   (approximate, linear mapping)
pub fn color_to_levelset(color: f64, h: f64) -> f64 {
    (2.0 * color - 1.0) * h
}

/// Advect the level-set field φ at particle i.
///
/// Dφ/Dt = −∇φ · u  (material derivative)
///
/// Uses the SPH approximation of the gradient.
pub fn levelset_advection_rate(
    i: usize,
    positions: &[[f64; 3]],
    velocities: &[[f64; 3]],
    phi: &[f64],
    masses: &[f64],
    densities: &[f64],
    h: f64,
) -> f64 {
    let n = positions.len();
    let mut dphi_dt = 0.0;
    let vi = velocities[i];
    for j in 0..n {
        if i == j {
            continue;
        }
        let dx = [
            positions[i][0] - positions[j][0],
            positions[i][1] - positions[j][1],
            positions[i][2] - positions[j][2],
        ];
        let r = (dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]).sqrt();
        if r >= 2.0 * h || r < 1e-14 {
            continue;
        }
        let dw = cubic_spline_gradient_3d(r, h);
        let dphi = phi[j] - phi[i];
        let dv_dot_dx = (vi[0] - velocities[j][0]) * dx[0]
            + (vi[1] - velocities[j][1]) * dx[1]
            + (vi[2] - velocities[j][2]) * dx[2];
        dphi_dt -= masses[j] / densities[j].max(1e-30) * dphi * dw / r * dv_dot_dx;
    }
    dphi_dt
}

/// Reinitialise a particle-based level-set using an algebraic re-distance.
///
/// Sets φ_new = sign(φ) · min_neighbour_distance  (approximate).
pub fn levelset_reinitialise_particle(phi_i: f64, min_dist_to_interface: f64) -> f64 {
    phi_i.signum() * min_dist_to_interface.abs()
}

// ============================================================================
// § 17  DENSITY INTERPOLATION SCHEMES
// ============================================================================

/// Weighted arithmetic mean density interpolation at a point.
///
/// ρ̄ = Σⱼ mⱼ W_ij
pub fn sph_density_sum(masses: &[f64], kernel_values: &[f64]) -> f64 {
    masses
        .iter()
        .zip(kernel_values.iter())
        .map(|(&m, &w)| m * w)
        .sum()
}

/// Continuity equation density update (Lagrangian form).
///
/// Dρᵢ/Dt = Σⱼ mⱼ (uᵢ − uⱼ) · ∇Wᵢⱼ
pub fn density_continuity_rate(
    i: usize,
    positions: &[[f64; 3]],
    velocities: &[[f64; 3]],
    masses: &[f64],
    h: f64,
) -> f64 {
    let n = positions.len();
    let mut drho_dt = 0.0;
    for j in 0..n {
        if i == j {
            continue;
        }
        let dx = [
            positions[i][0] - positions[j][0],
            positions[i][1] - positions[j][1],
            positions[i][2] - positions[j][2],
        ];
        let r = (dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]).sqrt();
        if r >= 2.0 * h || r < 1e-14 {
            continue;
        }
        let dw = cubic_spline_gradient_3d(r, h);
        let dv = [
            velocities[i][0] - velocities[j][0],
            velocities[i][1] - velocities[j][1],
            velocities[i][2] - velocities[j][2],
        ];
        let dv_dot_dx = dv[0] * dx[0] + dv[1] * dx[1] + dv[2] * dx[2];
        drho_dt += masses[j] * dw / r * dv_dot_dx;
    }
    drho_dt
}

// ============================================================================
// § 18  EQUATION OF STATE
// ============================================================================

/// Tait equation of state for weakly compressible SPH.
///
/// p = p₀ \[ (ρ/ρ₀)^γ − 1 \]
///
/// where p₀ = ρ₀ c₀² / γ.
pub fn tait_pressure(rho: f64, rho0: f64, c0: f64, gamma: f64) -> f64 {
    let b = rho0 * c0 * c0 / gamma;
    b * ((rho / rho0).powf(gamma) - 1.0)
}

/// Ideal gas equation of state: p = ρ R_spec T.
pub fn ideal_gas_pressure(rho: f64, r_spec: f64, temperature: f64) -> f64 {
    rho * r_spec * temperature
}

/// Sound speed in a Tait fluid: c = c₀ (ρ/ρ₀)^((γ−1)/2).
pub fn tait_sound_speed(rho: f64, rho0: f64, c0: f64, gamma: f64) -> f64 {
    c0 * (rho / rho0).powf((gamma - 1.0) * 0.5)
}

// ============================================================================
// § 19  VISCOSITY MODELS
// ============================================================================

/// Artificial viscosity term (Monaghan 1992).
///
/// Πᵢⱼ = (−α c̄ μᵢⱼ + β μᵢⱼ²) / ρ̄   if vᵢⱼ · rᵢⱼ < 0, else 0
///
/// where μᵢⱼ = h vᵢⱼ·rᵢⱼ / (|rᵢⱼ|² + ε h²).
pub fn artificial_viscosity(
    dx: &[f64; 3],
    dv: &[f64; 3],
    r: f64,
    rho_mean: f64,
    c_mean: f64,
    h: f64,
    alpha: f64,
    beta: f64,
    eps: f64,
) -> f64 {
    let vdotr = dv[0] * dx[0] + dv[1] * dx[1] + dv[2] * dx[2];
    if vdotr >= 0.0 {
        return 0.0;
    }
    let mu_ij = h * vdotr / (r * r + eps * h * h);
    (-alpha * c_mean * mu_ij + beta * mu_ij * mu_ij) / rho_mean.max(1e-30)
}

/// Laminar viscous SPH term: momentum change rate due to viscosity.
///
/// (∂u/∂t)_visc = Σⱼ mⱼ (μᵢ + μⱼ)/(ρᵢ ρⱼ) (uⱼ − uᵢ) ∇²Wᵢⱼ_approx
pub fn sph_viscous_force_component(
    dv_comp: f64,
    mj: f64,
    rhoi: f64,
    rhoj: f64,
    mu_i: f64,
    mu_j: f64,
    laplacian_w: f64,
) -> f64 {
    mj * (mu_i + mu_j) / (rhoi * rhoj).max(1e-30) * dv_comp * laplacian_w
}

// ============================================================================
// § 20  TIME INTEGRATION (LEAPFROG / VERLET)
// ============================================================================

/// Leapfrog velocity update (kick step).
///
/// v^(n+1/2) = v^(n−1/2) + dt · a^n
pub fn leapfrog_kick(vel: &[f64; 3], acc: &[f64; 3], dt: f64) -> [f64; 3] {
    [
        vel[0] + dt * acc[0],
        vel[1] + dt * acc[1],
        vel[2] + dt * acc[2],
    ]
}

/// Leapfrog position update (drift step).
///
/// x^(n+1) = x^n + dt · v^(n+1/2)
pub fn leapfrog_drift(pos: &[f64; 3], vel: &[f64; 3], dt: f64) -> [f64; 3] {
    [
        pos[0] + dt * vel[0],
        pos[1] + dt * vel[1],
        pos[2] + dt * vel[2],
    ]
}

/// CFL time step constraint for SPH.
///
/// dt_cfl = C_cfl · h / (c_max + v_max)
pub fn cfl_timestep(h: f64, c_max: f64, v_max: f64, c_cfl: f64) -> f64 {
    c_cfl * h / (c_max + v_max).max(1e-30)
}

// ============================================================================
// § 21  MULTIPHASE SIMULATION STATE
// ============================================================================

/// Multiphase SPH simulation configuration.
#[derive(Debug, Clone)]
pub struct MultiphaseConfig {
    /// Smoothing length h \[m\]
    pub h: f64,
    /// Liquid density ρ_l \[kg/m³\]
    pub rho_liquid: f64,
    /// Gas density ρ_g \[kg/m³\]
    pub rho_gas: f64,
    /// Liquid dynamic viscosity μ_l \[Pa·s\]
    pub mu_liquid: f64,
    /// Gas dynamic viscosity μ_g \[Pa·s\]
    pub mu_gas: f64,
    /// Surface tension coefficient σ \[N/m\]
    pub sigma: f64,
    /// Gravitational acceleration g \[m/s²\]
    pub gravity: [f64; 3],
    /// Speed of sound c₀ \[m/s\]
    pub c0: f64,
    /// EOS gamma exponent
    pub gamma: f64,
    /// CFL coefficient
    pub cfl: f64,
}

impl Default for MultiphaseConfig {
    fn default() -> Self {
        Self {
            h: 0.01,
            rho_liquid: 1000.0,
            rho_gas: 1.225,
            mu_liquid: 1e-3,
            mu_gas: 1.8e-5,
            sigma: 0.072,
            gravity: [0.0, -9.81, 0.0],
            c0: 1480.0,
            gamma: 7.0,
            cfl: 0.4,
        }
    }
}

impl MultiphaseConfig {
    /// Compute the Atwood number for this configuration.
    pub fn atwood(&self) -> f64 {
        atwood_number(self.rho_liquid, self.rho_gas)
    }

    /// Compute the critical RT wavenumber.
    pub fn rt_critical_wavenumber(&self) -> f64 {
        let g_mag =
            (self.gravity[0].powi(2) + self.gravity[1].powi(2) + self.gravity[2].powi(2)).sqrt();
        rayleigh_taylor_critical_wavenumber(self.atwood(), self.rho_liquid, g_mag, self.sigma)
    }
}

// ============================================================================
// § 22  PARTICLE NEIGHBOUR LIST
// ============================================================================

/// Find neighbours of particle i within radius 2h.
///
/// Returns a vector of neighbour indices.
pub fn find_neighbours(i: usize, positions: &[[f64; 3]], h: f64) -> Vec<usize> {
    let pi = positions[i];
    positions
        .iter()
        .enumerate()
        .filter(|&(j, pj)| {
            if i == j {
                return false;
            }
            let dx = pi[0] - pj[0];
            let dy = pi[1] - pj[1];
            let dz = pi[2] - pj[2];
            dx * dx + dy * dy + dz * dz < (2.0 * h) * (2.0 * h)
        })
        .map(|(j, _)| j)
        .collect()
}

/// Build a full neighbour list for all particles.
pub fn build_neighbour_list(positions: &[[f64; 3]], h: f64) -> Vec<Vec<usize>> {
    (0..positions.len())
        .map(|i| find_neighbours(i, positions, h))
        .collect()
}

// ============================================================================
// § 23  UTILITY FUNCTIONS
// ============================================================================

/// Linear interpolation between two values.
#[inline]
pub fn lerp(a: f64, b: f64, t: f64) -> f64 {
    a + t * (b - a)
}

/// 3-D dot product.
#[inline]
pub fn dot3(a: &[f64; 3], b: &[f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// 3-D cross product.
#[inline]
pub fn cross3(a: &[f64; 3], b: &[f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// Normalise a 3-D vector; returns the zero vector if below tolerance.
#[inline]
pub fn normalise3(v: &[f64; 3]) -> [f64; 3] {
    let norm = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
    if norm < 1e-14 {
        [0.0; 3]
    } else {
        [v[0] / norm, v[1] / norm, v[2] / norm]
    }
}

/// Distance between two 3-D points.
#[inline]
pub fn dist3(a: &[f64; 3], b: &[f64; 3]) -> f64 {
    let dx = a[0] - b[0];
    let dy = a[1] - b[1];
    let dz = a[2] - b[2];
    (dx * dx + dy * dy + dz * dz).sqrt()
}

// ============================================================================
// TESTS
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_cubic_spline_3d_at_zero() {
        // W(0, h) should equal alpha = 1/(πh³)
        let h = 0.1;
        let w = cubic_spline_3d(0.0, h);
        let expected = 1.0 / (PI * h * h * h);
        assert!(
            (w - expected).abs() < 1e-10,
            "W(0,h)={:.6} expected {:.6}",
            w,
            expected
        );
    }

    #[test]
    fn test_cubic_spline_3d_cutoff() {
        let h = 0.1;
        assert_eq!(cubic_spline_3d(2.0 * h + 1e-10, h), 0.0);
    }

    #[test]
    fn test_wendland_c2_3d_at_zero() {
        let h = 0.05;
        let w = wendland_c2_3d(0.0, h);
        let expected = 21.0 / (2.0 * PI * h * h * h);
        assert!(
            (w - expected).abs() < 1e-8,
            "Wendland at 0: {:.6} vs {:.6}",
            w,
            expected
        );
    }

    #[test]
    fn test_wendland_c2_3d_cutoff() {
        assert_eq!(wendland_c2_3d(0.2, 0.05), 0.0);
    }

    #[test]
    fn test_color_particle_ext_new_liquid() {
        let p = ColorParticleExt::new_liquid(1.0, 2.0, 3.0, 0.001);
        assert_eq!(p.color, 1.0);
        assert_eq!(p.phase, PhaseLabel::Liquid);
        assert!((p.density - 1000.0).abs() < 1e-10);
    }

    #[test]
    fn test_color_particle_ext_new_gas() {
        let p = ColorParticleExt::new_gas(0.0, 0.0, 0.0, 1e-6);
        assert_eq!(p.color, 0.0);
        assert_eq!(p.phase, PhaseLabel::Gas);
    }

    #[test]
    fn test_interpolated_density() {
        let mut p = ColorParticleExt::new_liquid(0.0, 0.0, 0.0, 0.001);
        p.color = 0.5;
        let rho = p.interpolated_density(1000.0, 1.225);
        assert!((rho - 500.6125).abs() < 0.01);
    }

    #[test]
    fn test_csf_force_zero_curvature() {
        let n = [0.0, 1.0, 0.0];
        let f = csf_surface_tension_force(&n, 0.0, 1.0, 0.072);
        assert_eq!(f, [0.0, 0.0, 0.0]);
    }

    #[test]
    fn test_csf_force_direction() {
        let n = [1.0, 0.0, 0.0];
        let f = csf_surface_tension_force(&n, 100.0, 1.0, 0.072);
        assert!(f[0] > 0.0 && f[1].abs() < 1e-14);
    }

    #[test]
    fn test_atwood_number_equal_densities() {
        let at = atwood_number(1000.0, 1000.0);
        assert!((at - 0.0).abs() < 1e-12);
    }

    #[test]
    fn test_atwood_number_heavy_light() {
        let at = atwood_number(1000.0, 1.225);
        assert!(
            at > 0.99,
            "water/air Atwood should be close to 1, got {:.6}",
            at
        );
    }

    #[test]
    fn test_rayleigh_taylor_growth_rate_positive() {
        let rate = rayleigh_taylor_growth_rate(0.998, 9.81, 10.0);
        assert!(rate > 0.0);
    }

    #[test]
    fn test_kelvin_helmholtz_growth_rate_zero_velocity() {
        let rate = kelvin_helmholtz_growth_rate(10.0, 0.0, 1000.0, 1.225);
        assert_eq!(rate, 0.0);
    }

    #[test]
    fn test_droplet_radius_from_volume_consistency() {
        let r = 0.005_f64;
        let v = droplet_volume(r);
        let r2 = droplet_radius_from_volume(v);
        assert!((r - r2).abs() < 1e-12, "round-trip radius: {} vs {}", r, r2);
    }

    #[test]
    fn test_weber_number() {
        let we = weber_number(1000.0, 2.0, 0.005, 0.072);
        // ρ U² D / σ = 1000 * 4 * 0.005 / 0.072 ≈ 277.8
        assert!((we - 277.78).abs() < 0.1, "We={:.6}", we);
    }

    #[test]
    fn test_droplet_coalescence_overlap() {
        let pa = [0.0, 0.0, 0.0];
        let pb = [0.003, 0.0, 0.0];
        assert!(droplets_coalesce(&pa, &pb, 0.002, 0.002, 0.0));
    }

    #[test]
    fn test_droplet_no_coalescence_far() {
        let pa = [0.0, 0.0, 0.0];
        let pb = [1.0, 0.0, 0.0];
        assert!(!droplets_coalesce(&pa, &pb, 0.002, 0.002, 0.0));
    }

    #[test]
    fn test_rayleigh_plesset_step_radius_positive() {
        let r = 1e-3;
        let r_dot = 0.0;
        let (r_new, _) = rayleigh_plesset_step(r, r_dot, 1e5, 1e5, 1000.0, 0.072, 1e-3, 1e-7);
        assert!(r_new > 0.0);
    }

    #[test]
    fn test_minnaert_frequency_positive() {
        let f = minnaert_frequency(1e-3, 101325.0, 1000.0, 1.4);
        assert!(f > 0.0, "Minnaert frequency={:.6}", f);
    }

    #[test]
    fn test_jurin_height_positive() {
        let h = jurin_height(0.072, 0.0, 1000.0, 9.81, 1e-3);
        assert!(h > 0.0, "Jurin height = {:.6}", h);
    }

    #[test]
    fn test_lucas_washburn_height_grows_with_time() {
        let h1 = lucas_washburn_height(0.072, 0.0, 1e-3, 1e-3, 1.0);
        let h2 = lucas_washburn_height(0.072, 0.0, 1e-3, 1e-3, 4.0);
        assert!(h2 > h1, "height should grow with time: {} vs {}", h1, h2);
    }

    #[test]
    fn test_young_contact_angle_hydrophilic() {
        // cos θ = 0.9 → θ ≈ 25.8°
        let theta = young_contact_angle(0.09, 0.0, 0.1);
        assert!(theta < 0.5, "small angle for hydrophilic: {:.6}", theta);
    }

    #[test]
    fn test_tait_pressure_at_reference() {
        let p = tait_pressure(1000.0, 1000.0, 1480.0, 7.0);
        assert!(
            p.abs() < 1e-6,
            "pressure at rho0 should be zero, got {:.6}",
            p
        );
    }

    #[test]
    fn test_cfl_timestep_positive() {
        let dt = cfl_timestep(0.01, 1480.0, 1.0, 0.4);
        assert!(dt > 0.0);
        let expected = 0.4 * 0.01 / 1481.0;
        assert!((dt - expected).abs() < 1e-12);
    }

    #[test]
    fn test_leapfrog_kick_drift() {
        let v = [1.0, 0.0, 0.0];
        let a = [0.0, -9.81, 0.0];
        let dt = 0.01;
        let v2 = leapfrog_kick(&v, &a, dt);
        assert!((v2[1] - (-0.0981)).abs() < 1e-12);
        let x = [0.0; 3];
        let x2 = leapfrog_drift(&x, &v2, dt);
        assert!((x2[0] - 0.01).abs() < 1e-12);
    }

    #[test]
    fn test_multiphase_config_atwood() {
        let cfg = MultiphaseConfig::default();
        let at = cfg.atwood();
        assert!(at > 0.99, "water/air Atwood: {:.6}", at);
    }

    #[test]
    fn test_find_neighbours_count() {
        let positions: Vec<[f64; 3]> = (0..5).map(|i| [i as f64 * 0.005, 0.0, 0.0]).collect();
        let nb = find_neighbours(2, &positions, 0.01);
        // Particle 2 at 0.01; h=0.01 → 2h=0.02. Neighbours at 0,1,3,4?
        // dist to 1: 0.005, to 3: 0.005, to 0: 0.010, to 4: 0.010 — all < 0.02
        assert!(nb.len() >= 2);
    }

    #[test]
    fn test_dot3_orthogonal() {
        let a = [1.0, 0.0, 0.0];
        let b = [0.0, 1.0, 0.0];
        assert_eq!(dot3(&a, &b), 0.0);
    }

    #[test]
    fn test_cross3_basis() {
        let x = [1.0, 0.0, 0.0];
        let y = [0.0, 1.0, 0.0];
        let z = cross3(&x, &y);
        assert!((z[2] - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_normalise3() {
        let v = [3.0, 4.0, 0.0];
        let n = normalise3(&v);
        let len = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
        assert!((len - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_color_to_levelset() {
        let phi_liquid = color_to_levelset(1.0, 0.01);
        let phi_gas = color_to_levelset(0.0, 0.01);
        assert!(phi_liquid > 0.0);
        assert!(phi_gas < 0.0);
    }

    #[test]
    fn test_phase_inversion_threshold() {
        // Equal viscosities → inversion at 0.5
        let frac = phase_inversion_fraction(1.0, 1.0);
        assert!((frac - 0.5).abs() < 1e-12);
    }

    #[test]
    fn test_capillary_number() {
        let ca = capillary_number(1e-3, 0.1, 0.072);
        assert!((ca - 1e-3 * 0.1 / 0.072).abs() < 1e-12);
    }

    #[test]
    fn test_bond_number() {
        let bo = bond_number(1000.0, 9.81, 0.001, 0.072);
        // ρ g l² / σ = 1000*9.81*1e-6/0.072 ≈ 0.13625
        assert!((bo - 0.13625).abs() < 0.001, "Bo={:.6}", bo);
    }

    // Verify PI constant is used (avoids any dead_code lint on import)
    #[test]
    fn test_pi_usage() {
        let sphere_area = 4.0 * PI * 1.0_f64.powi(2);
        assert!((sphere_area - 4.0 * PI).abs() < 1e-12);
    }
}
