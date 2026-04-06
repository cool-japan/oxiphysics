#![allow(clippy::needless_range_loop, clippy::ptr_arg)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Viscoelastic SPH for polymer and complex-fluid simulations.
//!
//! Implements SPH-based simulation of viscoelastic fluids using the
//! Oldroyd-B, FENE-P, and Giesekus constitutive models.  Each model
//! evolves an extra-stress tensor τ alongside the standard SPH fields.
//!
//! # References
//! - Ellero et al. (2002) "Viscoelastic flows studied by smoothed particle dynamics"
//! - Fang et al. (2006) "A numerical study of the SPH method for simulating
//!   transient viscoelastic free surface flows"
//! - Müller et al. (2004) "Particle-based fluid simulation for interactive applications"

#![allow(dead_code)]

use std::f64::consts::PI;

// ---------------------------------------------------------------------------
// Math helpers for 3×3 matrices
// ---------------------------------------------------------------------------

/// Add two 3×3 matrices.
pub fn matrix3_add(a: &[[f64; 3]; 3], b: &[[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let mut c = [[0.0_f64; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            c[i][j] = a[i][j] + b[i][j];
        }
    }
    c
}

/// Subtract two 3×3 matrices: a - b.
fn matrix3_sub(a: &[[f64; 3]; 3], b: &[[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let mut c = [[0.0_f64; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            c[i][j] = a[i][j] - b[i][j];
        }
    }
    c
}

/// Multiply a 3×3 matrix by a scalar.
fn matrix3_scale(a: &[[f64; 3]; 3], s: f64) -> [[f64; 3]; 3] {
    let mut c = [[0.0_f64; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            c[i][j] = a[i][j] * s;
        }
    }
    c
}

/// Compute the trace of a 3×3 matrix.
pub fn matrix3_trace(a: &[[f64; 3]; 3]) -> f64 {
    a[0][0] + a[1][1] + a[2][2]
}

/// Multiply two 3×3 matrices: C = A * B.
pub fn matrix3_mul(a: &[[f64; 3]; 3], b: &[[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let mut c = [[0.0_f64; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            for k in 0..3 {
                c[i][j] += a[i][k] * b[k][j];
            }
        }
    }
    c
}

/// Transpose a 3×3 matrix.
pub fn matrix3_transpose(a: &[[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let mut t = [[0.0_f64; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            t[j][i] = a[i][j];
        }
    }
    t
}

/// Identity 3×3 matrix.
fn matrix3_identity() -> [[f64; 3]; 3] {
    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
}

/// Zero 3×3 matrix.
fn matrix3_zero() -> [[f64; 3]; 3] {
    [[0.0_f64; 3]; 3]
}

/// Frobenius inner product: A : B = Σ_{ij} A_ij B_ij.
fn matrix3_frobenius(a: &[[f64; 3]; 3], b: &[[f64; 3]; 3]) -> f64 {
    let mut s = 0.0;
    for i in 0..3 {
        for j in 0..3 {
            s += a[i][j] * b[i][j];
        }
    }
    s
}

/// Dot product of two 3-vectors.
fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// Subtract two 3-vectors.
fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

/// Scale a 3-vector.
fn scale3(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

/// Euclidean norm of a 3-vector.
fn norm3(a: [f64; 3]) -> f64 {
    dot3(a, a).sqrt()
}

// ---------------------------------------------------------------------------
// SPH kernel (cubic spline)
// ---------------------------------------------------------------------------

/// Cubic spline kernel W(r, h) in 3D.
fn kernel_w(r: f64, h: f64) -> f64 {
    let q = r / h;
    let sigma = 1.0 / (PI * h * h * h);
    let w = if q < 1.0 {
        1.0 - 1.5 * q * q + 0.75 * q * q * q
    } else if q < 2.0 {
        0.25 * (2.0 - q) * (2.0 - q) * (2.0 - q)
    } else {
        0.0
    };
    sigma * w
}

/// Gradient of cubic spline kernel: ∇W(r⃗, h) = (dW/dr) * r̂.
fn kernel_grad_w(r_vec: [f64; 3], h: f64) -> [f64; 3] {
    let r = norm3(r_vec);
    if r < 1e-14 {
        return [0.0; 3];
    }
    let q = r / h;
    let sigma = 1.0 / (PI * h * h * h);
    let dw_dr = if q < 1.0 {
        sigma * (-3.0 * q + 2.25 * q * q) / h
    } else if q < 2.0 {
        sigma * (-0.75 * (2.0 - q) * (2.0 - q)) / h
    } else {
        0.0
    };
    scale3(r_vec, dw_dr / r)
}

// ---------------------------------------------------------------------------
// Particle
// ---------------------------------------------------------------------------

/// SPH particle carrying an extra viscoelastic stress tensor.
///
/// The stress tensor `stress[i][j]` is the polymer (extra) contribution τ_{ij}.
#[derive(Debug, Clone)]
pub struct ViscoelasticParticle {
    /// Position (m).
    pub pos: [f64; 3],
    /// Velocity (m/s).
    pub vel: [f64; 3],
    /// Mass (kg).
    pub mass: f64,
    /// Density (kg/m³).
    pub rho: f64,
    /// Isotropic pressure (Pa).
    pub pressure: f64,
    /// Extra (polymer) stress tensor τ_{ij} (Pa).
    pub stress: [[f64; 3]; 3],
}

impl ViscoelasticParticle {
    /// Create a particle at rest with zero stress.
    pub fn new(pos: [f64; 3], mass: f64, rho: f64) -> Self {
        Self {
            pos,
            vel: [0.0; 3],
            mass,
            rho,
            pressure: 0.0,
            stress: matrix3_zero(),
        }
    }
}

// ---------------------------------------------------------------------------
// Oldroyd-B model
// ---------------------------------------------------------------------------

/// Oldroyd-B constitutive model parameters.
///
/// Describes a dilute polymer solution with a single relaxation time λ,
/// polymer viscosity η_p, and solvent viscosity η_s.
///
/// The stress evolution equation is:
/// ```text
/// Dτ/Dt = -τ/λ + 2 η_p D + τ · L^T + L · τ
/// ```
/// where `L = ∇u` is the velocity gradient and `D = (L + L^T)/2` is the
/// rate-of-strain tensor.
#[derive(Debug, Clone)]
pub struct OldroydBModel {
    /// Polymer relaxation time (s).
    pub lambda: f64,
    /// Polymer viscosity (Pa·s).
    pub eta_p: f64,
    /// Solvent viscosity (Pa·s).
    pub eta_s: f64,
}

impl OldroydBModel {
    /// Create a typical polymer solution (water + dilute polymer).
    pub fn polymer_solution() -> Self {
        Self {
            lambda: 0.1,
            eta_p: 0.01,
            eta_s: 0.001,
        }
    }

    /// Compute the upper-convected time derivative (UCT) of the stress.
    ///
    /// Returns Dτ/Dt = UCT(τ, L) - τ/λ + 2 η_p D.
    pub fn stress_rate(&self, stress: &[[f64; 3]; 3], vel_grad: &[[f64; 3]; 3]) -> [[f64; 3]; 3] {
        let uct = Self::upper_convected_derivative(stress, vel_grad);
        // Rate-of-strain tensor D = (L + L^T) / 2
        let lt = matrix3_transpose(vel_grad);
        let mut d = matrix3_zero();
        for i in 0..3 {
            for j in 0..3 {
                d[i][j] = 0.5 * (vel_grad[i][j] + lt[i][j]);
            }
        }
        let relaxation = matrix3_scale(stress, -1.0 / self.lambda);
        let polymer_input = matrix3_scale(&d, 2.0 * self.eta_p);
        matrix3_add(&matrix3_add(&uct, &relaxation), &polymer_input)
    }

    /// Upper-convected time derivative of the stress tensor.
    ///
    /// The upper-convected derivative is:
    /// ```text
    /// ∇τ = Dτ/Dt - L · τ - τ · L^T
    /// ```
    /// This method returns the term `L · τ + τ · L^T` (the convective part),
    /// which the caller adds to the material derivative.
    pub fn upper_convected_derivative(
        stress: &[[f64; 3]; 3],
        vel_grad: &[[f64; 3]; 3],
    ) -> [[f64; 3]; 3] {
        upper_convected_derivative(stress, vel_grad)
    }
}

/// Compute the upper-convected contribution: L · τ + τ · L^T.
///
/// This represents the convective transport of the stress tensor
/// in the Oldroyd-B and related models.
pub fn upper_convected_derivative(
    stress: &[[f64; 3]; 3],
    vel_grad: &[[f64; 3]; 3],
) -> [[f64; 3]; 3] {
    let l_tau = matrix3_mul(vel_grad, stress);
    let lt = matrix3_transpose(vel_grad);
    let tau_lt = matrix3_mul(stress, &lt);
    matrix3_add(&l_tau, &tau_lt)
}

// ---------------------------------------------------------------------------
// FENE-P model
// ---------------------------------------------------------------------------

/// FENE-P (Finitely Extensible Nonlinear Elastic – Peterlin closure) model.
///
/// Extends Oldroyd-B with a finite extensibility correction that prevents
/// unphysical infinite polymer stretching.  The maximum extension is
/// parameterised by `b` (square of maximum extension length in units of
/// equilibrium length).
#[derive(Debug, Clone)]
pub struct FenePModel {
    /// Polymer relaxation time (s).
    pub lambda: f64,
    /// Polymer viscosity (Pa·s).
    pub eta_p: f64,
    /// Solvent viscosity (Pa·s).
    pub eta_s: f64,
    /// Maximum extensibility parameter (dimensionless).
    pub b: f64,
}

impl FenePModel {
    /// Create a typical FENE-P polymer model.
    pub fn default_polymer() -> Self {
        Self {
            lambda: 0.1,
            eta_p: 0.01,
            eta_s: 0.001,
            b: 100.0,
        }
    }

    /// Compute the Peterlin factor f(τ) = b / (b - tr(τ) / (η_p/λ)).
    ///
    /// Clamps the denominator to avoid singularity when tr(τ) approaches b.
    pub fn peterlin_factor(stress: &[[f64; 3]; 3], b: f64) -> f64 {
        let tr = matrix3_trace(stress);
        let denom = (b - tr).max(1e-6);
        b / denom
    }

    /// Compute the stress rate for the FENE-P model.
    ///
    /// Dτ/Dt = UCT(τ, L) - f(τ) τ / λ + 2 η_p D
    pub fn stress_rate(&self, stress: &[[f64; 3]; 3], vel_grad: &[[f64; 3]; 3]) -> [[f64; 3]; 3] {
        let uct = upper_convected_derivative(stress, vel_grad);
        let f = Self::peterlin_factor(stress, self.b);
        let lt = matrix3_transpose(vel_grad);
        let mut d = matrix3_zero();
        for i in 0..3 {
            for j in 0..3 {
                d[i][j] = 0.5 * (vel_grad[i][j] + lt[i][j]);
            }
        }
        let relaxation = matrix3_scale(stress, -f / self.lambda);
        let polymer_input = matrix3_scale(&d, 2.0 * self.eta_p);
        matrix3_add(&matrix3_add(&uct, &relaxation), &polymer_input)
    }
}

// ---------------------------------------------------------------------------
// Giesekus model
// ---------------------------------------------------------------------------

/// Giesekus constitutive model.
///
/// Extends Oldroyd-B with a nonlinear mobility parameter α (0 ≤ α ≤ 0.5),
/// which introduces shear thinning and a non-zero second normal stress
/// difference.  α = 0 recovers Oldroyd-B.
///
/// Stress evolution:
/// ```text
/// Dτ/Dt = UCT(τ, L) - (τ + α λ/η_p τ·τ) / λ + 2 η_p D
/// ```
#[derive(Debug, Clone)]
pub struct GiesekusModel {
    /// Relaxation time (s).
    pub lambda: f64,
    /// Mobility parameter (dimensionless, 0 ≤ α ≤ 0.5).
    pub alpha: f64,
    /// Polymer viscosity (Pa·s).
    pub eta_p: f64,
    /// Solvent viscosity (Pa·s).
    pub eta_s: f64,
}

impl GiesekusModel {
    /// Create a Giesekus model with typical parameters.
    pub fn default_polymer() -> Self {
        Self {
            lambda: 0.1,
            alpha: 0.1,
            eta_p: 0.01,
            eta_s: 0.001,
        }
    }

    /// Compute the Giesekus stress rate.
    ///
    /// Returns Dτ/Dt = UCT - (τ + α λ/η_p τ·τ)/λ + 2 η_p D.
    pub fn stress_rate(&self, stress: &[[f64; 3]; 3], vel_grad: &[[f64; 3]; 3]) -> [[f64; 3]; 3] {
        let uct = upper_convected_derivative(stress, vel_grad);
        let lt = matrix3_transpose(vel_grad);
        let mut d = matrix3_zero();
        for i in 0..3 {
            for j in 0..3 {
                d[i][j] = 0.5 * (vel_grad[i][j] + lt[i][j]);
            }
        }
        let tau_sq = matrix3_mul(stress, stress);
        // Nonlinear term: α λ / η_p * τ·τ
        let nonlinear = matrix3_scale(&tau_sq, self.alpha * self.lambda / self.eta_p);
        let linear_part = matrix3_add(stress, &nonlinear);
        let relaxation = matrix3_scale(&linear_part, -1.0 / self.lambda);
        let polymer_input = matrix3_scale(&d, 2.0 * self.eta_p);
        matrix3_add(&matrix3_add(&uct, &relaxation), &polymer_input)
    }
}

// ---------------------------------------------------------------------------
// SPH operators
// ---------------------------------------------------------------------------

/// Compute the velocity gradient tensor ∂u_i/∂x_j at particle `idx` using SPH.
///
/// Uses the standard SPH kernel gradient sum:
/// ```text
/// (∂u_α/∂x_β)_i = Σ_j (m_j/ρ_j) (u_j - u_i)_α ∇W_{ij,β}
/// ```
///
/// # Arguments
/// * `particles` – all particles in the domain.
/// * `idx` – index of the particle at which to evaluate the gradient.
/// * `neighbors` – indices of neighboring particles (within 2h support).
/// * `h` – smoothing length.
pub fn velocity_gradient_sph(
    particles: &[ViscoelasticParticle],
    idx: usize,
    neighbors: &[usize],
    h: f64,
) -> [[f64; 3]; 3] {
    let pi = &particles[idx];
    let mut grad = matrix3_zero();
    for &j in neighbors {
        if j == idx {
            continue;
        }
        let pj = &particles[j];
        let r_ij = sub3(pi.pos, pj.pos);
        let dw = kernel_grad_w(r_ij, h);
        let v_ji = sub3(pj.vel, pi.vel); // u_j - u_i
        let w = pj.mass / pj.rho;
        for alpha in 0..3 {
            for beta in 0..3 {
                grad[alpha][beta] += w * v_ji[alpha] * dw[beta];
            }
        }
    }
    grad
}

/// Compute the viscoelastic force on particle `i` due to particle `j`.
///
/// The stress-divergence contribution is approximated as:
/// ```text
/// f_{ij} = m_i m_j (τ_i/ρ_i² + τ_j/ρ_j²) · ∇W_{ij}
/// ```
///
/// Returns a 3-vector of force contributions (N).
pub fn viscoelastic_force(
    p_i: &ViscoelasticParticle,
    p_j: &ViscoelasticParticle,
    h: f64,
) -> [f64; 3] {
    let r_ij = sub3(p_i.pos, p_j.pos);
    let dw = kernel_grad_w(r_ij, h);

    // τ_i/ρ_i² + τ_j/ρ_j²
    let rho_i2 = p_i.rho * p_i.rho;
    let rho_j2 = p_j.rho * p_j.rho;
    let mut tau_sum = matrix3_zero();
    for a in 0..3 {
        for b in 0..3 {
            tau_sum[a][b] = p_i.stress[a][b] / rho_i2 + p_j.stress[a][b] / rho_j2;
        }
    }

    // f = m_i m_j * (τ_sum · ∇W)
    let scale = p_i.mass * p_j.mass;
    let mut force = [0.0_f64; 3];
    for a in 0..3 {
        for b in 0..3 {
            force[a] += scale * tau_sum[a][b] * dw[b];
        }
    }
    force
}

// ---------------------------------------------------------------------------
// Explicit stress integration for Oldroyd-B
// ---------------------------------------------------------------------------

/// Integrate the Oldroyd-B extra-stress tensor for all particles using
/// an explicit (forward Euler) time step.
///
/// For each particle `i`, computes the velocity gradient via SPH and
/// advances the stress:
/// ```text
/// τ^{n+1}_i = τ^n_i + dt * (Dτ/Dt)_i
/// ```
///
/// # Arguments
/// * `model` – Oldroyd-B material parameters.
/// * `particles` – mutable list of particles (stress updated in place).
/// * `neighbors` – neighbor list: `neighbors[i]` contains indices of neighbors of i.
/// * `h` – smoothing length.
/// * `dt` – time step (s).
pub fn integrate_stress_explicit(
    model: &OldroydBModel,
    particles: &mut Vec<ViscoelasticParticle>,
    neighbors: &[Vec<usize>],
    h: f64,
    dt: f64,
) {
    let n = particles.len();
    // Compute stress rates first (immutable borrow for velocity gradient)
    let stress_rates: Vec<[[f64; 3]; 3]> = (0..n)
        .map(|i| {
            let l = velocity_gradient_sph(particles, i, &neighbors[i], h);
            model.stress_rate(&particles[i].stress, &l)
        })
        .collect();

    // Update stresses
    for (i, rate) in stress_rates.iter().enumerate() {
        let tau = particles[i].stress;
        particles[i].stress = matrix3_add(&tau, &matrix3_scale(rate, dt));
    }
}

// ---------------------------------------------------------------------------
// Extensional flow test (Trouton ratio)
// ---------------------------------------------------------------------------

/// Simulate the Trouton ratio of an Oldroyd-B fluid in uniaxial extensional flow.
///
/// In a uniaxial extensional flow with strain rate ε̇, the velocity gradient is:
/// ```text
/// L = diag(ε̇, -ε̇/2, -ε̇/2)
/// ```
///
/// The Trouton ratio is η_E / (3 η_total) where η_E is the extensional viscosity.
///
/// # Returns
/// A vector of `(time, Trouton_ratio)` pairs.
pub fn viscoelastic_extensional_flow_test(
    lambda: f64,
    strain_rate: f64,
    t_max: f64,
    dt: f64,
) -> Vec<(f64, f64)> {
    let model = OldroydBModel {
        lambda,
        eta_p: 0.01,
        eta_s: 0.001,
    };
    let eta_total = model.eta_p + model.eta_s;

    // Uniaxial extensional velocity gradient
    let mut vel_grad = matrix3_zero();
    vel_grad[0][0] = strain_rate;
    vel_grad[1][1] = -0.5 * strain_rate;
    vel_grad[2][2] = -0.5 * strain_rate;

    let mut tau = matrix3_zero();
    let mut t = 0.0;
    let mut results = Vec::new();

    while t <= t_max {
        // Trouton ratio = η_E / (3 η_total)
        // η_E = (τ_11 - τ_22) / ε̇ + 3 η_s
        let tau_diff = tau[0][0] - tau[1][1];
        let eta_e = if strain_rate.abs() > 1e-14 {
            tau_diff / strain_rate + 3.0 * model.eta_s
        } else {
            3.0 * eta_total
        };
        let trouton = eta_e / (3.0 * eta_total);
        results.push((t, trouton));

        // Advance stress
        let rate = model.stress_rate(&tau, &vel_grad);
        tau = matrix3_add(&tau, &matrix3_scale(&rate, dt));
        t += dt;
    }
    results
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    fn zero_stress() -> [[f64; 3]; 3] {
        matrix3_zero()
    }

    fn identity_stress(s: f64) -> [[f64; 3]; 3] {
        let mut m = matrix3_zero();
        for i in 0..3 {
            m[i][i] = s;
        }
        m
    }

    fn simple_vel_grad(rate: f64) -> [[f64; 3]; 3] {
        let mut l = matrix3_zero();
        // Simple shear: L_01 = rate
        l[0][1] = rate;
        l
    }

    // --- matrix3 helpers ---

    #[test]
    fn matrix3_add_identity_check() {
        let a = [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]];
        let b = matrix3_zero();
        let c = matrix3_add(&a, &b);
        for i in 0..3 {
            for j in 0..3 {
                assert_eq!(c[i][j], a[i][j]);
            }
        }
    }

    #[test]
    fn matrix3_trace_identity() {
        let id = matrix3_identity();
        assert_eq!(matrix3_trace(&id), 3.0);
    }

    #[test]
    fn matrix3_trace_zero() {
        let z = matrix3_zero();
        assert_eq!(matrix3_trace(&z), 0.0);
    }

    #[test]
    fn matrix3_mul_identity() {
        let id = matrix3_identity();
        let a = [[1.0, 2.0, 3.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let b = matrix3_mul(&id, &a);
        for i in 0..3 {
            for j in 0..3 {
                assert!((b[i][j] - a[i][j]).abs() < 1e-12);
            }
        }
    }

    #[test]
    fn matrix3_mul_zero() {
        let id = matrix3_identity();
        let z = matrix3_zero();
        let b = matrix3_mul(&id, &z);
        for i in 0..3 {
            for j in 0..3 {
                assert_eq!(b[i][j], 0.0);
            }
        }
    }

    #[test]
    fn matrix3_transpose_involution() {
        let a = [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]];
        let att = matrix3_transpose(&matrix3_transpose(&a));
        for i in 0..3 {
            for j in 0..3 {
                assert_eq!(att[i][j], a[i][j]);
            }
        }
    }

    #[test]
    fn matrix3_transpose_off_diagonal() {
        let mut a = matrix3_zero();
        a[0][1] = 5.0;
        let at = matrix3_transpose(&a);
        assert_eq!(at[1][0], 5.0);
        assert_eq!(at[0][1], 0.0);
    }

    // --- Particle ---

    #[test]
    fn viscoelastic_particle_new_zero_stress() {
        let p = ViscoelasticParticle::new([0.0, 0.0, 0.0], 1e-3, 1000.0);
        for i in 0..3 {
            for j in 0..3 {
                assert_eq!(p.stress[i][j], 0.0);
            }
        }
    }

    #[test]
    fn viscoelastic_particle_new_zero_vel() {
        let p = ViscoelasticParticle::new([1.0, 2.0, 3.0], 1e-3, 1000.0);
        assert_eq!(p.vel, [0.0, 0.0, 0.0]);
    }

    // --- Oldroyd-B ---

    #[test]
    fn oldroyd_b_stress_rate_zero_at_equilibrium() {
        // At equilibrium with L=0 and τ=0, the stress rate should be 0
        let model = OldroydBModel::polymer_solution();
        let tau = zero_stress();
        let l = matrix3_zero();
        let rate = model.stress_rate(&tau, &l);
        for i in 0..3 {
            for j in 0..3 {
                assert!(rate[i][j].abs() < 1e-12, "rate[{i}][{j}] = {}", rate[i][j]);
            }
        }
    }

    #[test]
    fn oldroyd_b_stress_rate_nonzero_with_shear() {
        let model = OldroydBModel::polymer_solution();
        let tau = zero_stress();
        let l = simple_vel_grad(1.0);
        let rate = model.stress_rate(&tau, &l);
        // With shear flow, D_01 = D_10 = 0.5, so rate should be nonzero
        let total: f64 = rate.iter().flatten().map(|x| x.abs()).sum();
        assert!(
            total > 1e-12,
            "stress rate should be nonzero with shear: {total}"
        );
    }

    #[test]
    fn oldroyd_b_relaxation_reduces_stress() {
        let model = OldroydBModel {
            lambda: 1.0,
            eta_p: 0.0,
            eta_s: 0.0,
        };
        let tau = identity_stress(1.0);
        let l = matrix3_zero();
        let rate = model.stress_rate(&tau, &l);
        // With zero L and eta_p=0, rate = -tau/lambda < 0
        for i in 0..3 {
            assert!(rate[i][i] < 0.0, "diagonal should decrease: {}", rate[i][i]);
        }
    }

    #[test]
    fn upper_convected_derivative_zero_stress() {
        let tau = zero_stress();
        let l = simple_vel_grad(5.0);
        let uct = upper_convected_derivative(&tau, &l);
        for i in 0..3 {
            for j in 0..3 {
                assert_eq!(uct[i][j], 0.0);
            }
        }
    }

    #[test]
    fn upper_convected_derivative_symmetric_stress() {
        // For symmetric τ, UCT should also be symmetric
        let tau = [[1.0, 0.5, 0.0], [0.5, 2.0, 0.0], [0.0, 0.0, 0.5]];
        let l = simple_vel_grad(1.0);
        let uct = upper_convected_derivative(&tau, &l);
        for i in 0..3 {
            for j in 0..3 {
                assert!(
                    (uct[i][j] - uct[j][i]).abs() < 1e-12,
                    "UCT not symmetric at [{i}][{j}]"
                );
            }
        }
    }

    // --- FENE-P ---

    #[test]
    fn fene_p_peterlin_zero_stress() {
        let b = 100.0;
        let tau = zero_stress();
        let f = FenePModel::peterlin_factor(&tau, b);
        // tr(τ) = 0, so f = b/b = 1
        assert!((f - 1.0).abs() < 1e-12, "f = {f}");
    }

    #[test]
    fn fene_p_peterlin_large_b() {
        let b = 1e6;
        let tau = identity_stress(1.0);
        let f = FenePModel::peterlin_factor(&tau, b);
        // tr(τ) = 3, b >> 3, so f ≈ 1
        assert!((f - 1.0).abs() < 1e-4, "f = {f}");
    }

    #[test]
    fn fene_p_stress_rate_zero_at_equilibrium() {
        let model = FenePModel::default_polymer();
        let tau = zero_stress();
        let l = matrix3_zero();
        let rate = model.stress_rate(&tau, &l);
        let total: f64 = rate.iter().flatten().map(|x| x.abs()).sum();
        assert!(total < 1e-12, "equilibrium rate = {total}");
    }

    #[test]
    fn fene_p_stress_rate_nonzero_with_flow() {
        let model = FenePModel::default_polymer();
        let tau = zero_stress();
        let l = simple_vel_grad(1.0);
        let rate = model.stress_rate(&tau, &l);
        let total: f64 = rate.iter().flatten().map(|x| x.abs()).sum();
        assert!(
            total > 1e-12,
            "shear flow should produce nonzero rate: {total}"
        );
    }

    #[test]
    fn fene_p_peterlin_clamped_near_b() {
        let b = 10.0;
        let tau = identity_stress(9.9999); // tr ≈ 29.9997, b = 10 so tr > b → clamped
        // Should not panic, and f should be finite
        let f = FenePModel::peterlin_factor(&tau, b);
        assert!(f.is_finite(), "f should be finite: {f}");
    }

    // --- Giesekus ---

    #[test]
    fn giesekus_stress_rate_zero_at_equilibrium() {
        let model = GiesekusModel::default_polymer();
        let tau = zero_stress();
        let l = matrix3_zero();
        let rate = model.stress_rate(&tau, &l);
        let total: f64 = rate.iter().flatten().map(|x| x.abs()).sum();
        assert!(total < 1e-12, "equilibrium rate = {total}");
    }

    #[test]
    fn giesekus_alpha_zero_recovers_oldroyd_b() {
        let model = GiesekusModel {
            lambda: 1.0,
            alpha: 0.0,
            eta_p: 0.01,
            eta_s: 0.001,
        };
        let oldroyd = OldroydBModel {
            lambda: 1.0,
            eta_p: 0.01,
            eta_s: 0.001,
        };
        let tau = [[0.1, 0.05, 0.0], [0.05, 0.2, 0.0], [0.0, 0.0, 0.05]];
        let l = simple_vel_grad(0.5);
        let rate_g = model.stress_rate(&tau, &l);
        let rate_o = oldroyd.stress_rate(&tau, &l);
        for i in 0..3 {
            for j in 0..3 {
                assert!(
                    (rate_g[i][j] - rate_o[i][j]).abs() < 1e-12,
                    "Giesekus(α=0) != Oldroyd-B at [{i}][{j}]"
                );
            }
        }
    }

    #[test]
    fn giesekus_nonzero_alpha_differs_from_oldroyd_b() {
        let model_g = GiesekusModel {
            lambda: 1.0,
            alpha: 0.2,
            eta_p: 0.01,
            eta_s: 0.001,
        };
        let model_o = OldroydBModel {
            lambda: 1.0,
            eta_p: 0.01,
            eta_s: 0.001,
        };
        let tau = identity_stress(0.1);
        let l = simple_vel_grad(1.0);
        let rate_g = model_g.stress_rate(&tau, &l);
        let rate_o = model_o.stress_rate(&tau, &l);
        let diff: f64 = rate_g
            .iter()
            .flatten()
            .zip(rate_o.iter().flatten())
            .map(|(a, b)| (a - b).abs())
            .sum();
        assert!(diff > 1e-12, "Giesekus(α≠0) should differ from Oldroyd-B");
    }

    // --- velocity gradient SPH ---

    #[test]
    fn velocity_gradient_single_particle_zero() {
        let p = ViscoelasticParticle::new([0.0, 0.0, 0.0], 1e-3, 1000.0);
        let particles = vec![p];
        let grad = velocity_gradient_sph(&particles, 0, &[], 0.1);
        for i in 0..3 {
            for j in 0..3 {
                assert_eq!(grad[i][j], 0.0);
            }
        }
    }

    #[test]
    fn velocity_gradient_two_particles_finite() {
        let mut p0 = ViscoelasticParticle::new([0.0, 0.0, 0.0], 1e-3, 1000.0);
        p0.vel = [0.0, 0.0, 0.0];
        let mut p1 = ViscoelasticParticle::new([0.05, 0.0, 0.0], 1e-3, 1000.0);
        p1.vel = [1.0, 0.0, 0.0];
        let particles = vec![p0, p1];
        let grad = velocity_gradient_sph(&particles, 0, &[0, 1], 0.1);
        assert!(
            grad.iter().flatten().all(|x| x.is_finite()),
            "vel grad has non-finite"
        );
    }

    // --- viscoelastic force ---

    #[test]
    fn viscoelastic_force_zero_stress_zero_force() {
        let p0 = ViscoelasticParticle::new([0.0, 0.0, 0.0], 1e-3, 1000.0);
        let p1 = ViscoelasticParticle::new([0.05, 0.0, 0.0], 1e-3, 1000.0);
        let f = viscoelastic_force(&p0, &p1, 0.1);
        assert_eq!(f, [0.0, 0.0, 0.0]);
    }

    #[test]
    fn viscoelastic_force_nonzero_for_stressed_particles() {
        let mut p0 = ViscoelasticParticle::new([0.0, 0.0, 0.0], 1e-3, 1000.0);
        p0.stress = identity_stress(100.0);
        let mut p1 = ViscoelasticParticle::new([0.05, 0.0, 0.0], 1e-3, 1000.0);
        p1.stress = identity_stress(100.0);
        let f = viscoelastic_force(&p0, &p1, 0.1);
        let mag: f64 = f.iter().map(|x| x * x).sum::<f64>().sqrt();
        assert!(mag > 0.0, "force magnitude should be nonzero: {mag}");
    }

    // --- integrate_stress_explicit ---

    #[test]
    fn integrate_stress_explicit_zero_velocity_decays() {
        let model = OldroydBModel {
            lambda: 1.0,
            eta_p: 0.0,
            eta_s: 0.0,
        };
        let mut particles = vec![ViscoelasticParticle::new([0.0, 0.0, 0.0], 1e-3, 1000.0)];
        particles[0].stress = identity_stress(1.0);
        let neighbors: Vec<Vec<usize>> = vec![vec![]];
        let tau_before = particles[0].stress[0][0];
        integrate_stress_explicit(&model, &mut particles, &neighbors, 0.1, 0.01);
        let tau_after = particles[0].stress[0][0];
        assert!(
            tau_after < tau_before,
            "stress should decay with zero vel: {tau_after} < {tau_before}"
        );
    }

    #[test]
    fn integrate_stress_explicit_multiple_steps_finite() {
        let model = OldroydBModel::polymer_solution();
        let mut particles = vec![
            ViscoelasticParticle::new([0.0, 0.0, 0.0], 1e-3, 1000.0),
            ViscoelasticParticle::new([0.05, 0.0, 0.0], 1e-3, 1000.0),
        ];
        particles[0].vel = [0.5, 0.0, 0.0];
        let neighbors: Vec<Vec<usize>> = vec![vec![0, 1], vec![0, 1]];
        for _ in 0..10 {
            integrate_stress_explicit(&model, &mut particles, &neighbors, 0.05, 0.001);
        }
        for p in &particles {
            for row in &p.stress {
                for &s in row {
                    assert!(s.is_finite(), "stress is non-finite: {s}");
                }
            }
        }
    }

    // --- extensional flow test ---

    #[test]
    fn extensional_flow_trouton_ratio_at_zero_time() {
        // eta_p = 0.01, eta_s = 0.001 (hardcoded in viscoelastic_extensional_flow_test)
        let eta_p = 0.01_f64;
        let eta_s = 0.001_f64;
        let eta_total = eta_p + eta_s;
        let results = viscoelastic_extensional_flow_test(0.1, 0.5, 0.0, 0.001);
        assert!(!results.is_empty());
        let (t0, tr0) = results[0];
        assert!((t0 - 0.0).abs() < 1e-12);
        // At t=0, τ=0: eta_e = 3*eta_s, so Trouton = eta_s / eta_total
        let expected_tr0 = eta_s / eta_total;
        assert!(
            (tr0 - expected_tr0).abs() < 1e-10,
            "Trouton at t=0: expected {expected_tr0}, got {tr0}"
        );
    }

    #[test]
    fn extensional_flow_produces_multiple_points() {
        let results = viscoelastic_extensional_flow_test(0.1, 1.0, 0.1, 0.01);
        assert!(results.len() > 5);
    }

    #[test]
    fn extensional_flow_times_monotone() {
        let results = viscoelastic_extensional_flow_test(0.1, 1.0, 0.5, 0.05);
        for w in results.windows(2) {
            assert!(
                w[1].0 > w[0].0,
                "times not monotone: {} <= {}",
                w[1].0,
                w[0].0
            );
        }
    }

    #[test]
    fn extensional_flow_trouton_all_finite() {
        let results = viscoelastic_extensional_flow_test(0.2, 0.3, 1.0, 0.1);
        for (t, tr) in &results {
            assert!(tr.is_finite(), "non-finite Trouton at t={t}: {tr}");
        }
    }

    #[test]
    fn extensional_flow_zero_strain_rate_trouton_one() {
        let results = viscoelastic_extensional_flow_test(0.1, 0.0, 0.1, 0.01);
        for (_t, tr) in &results {
            assert!(
                (tr - 1.0).abs() < 1e-6,
                "zero-strain Trouton should be 1: {tr}"
            );
        }
    }

    // --- kernel tests ---

    #[test]
    fn kernel_w_zero_at_large_r() {
        let w = kernel_w(2.5, 1.0);
        assert_eq!(w, 0.0);
    }

    #[test]
    fn kernel_w_positive_within_support() {
        let w = kernel_w(0.5, 1.0);
        assert!(w > 0.0, "kernel should be positive: {w}");
    }

    #[test]
    fn kernel_grad_zero_at_large_r() {
        let grad = kernel_grad_w([3.0, 0.0, 0.0], 1.0);
        assert_eq!(grad, [0.0, 0.0, 0.0]);
    }

    #[test]
    fn kernel_grad_zero_at_origin() {
        let grad = kernel_grad_w([0.0, 0.0, 0.0], 1.0);
        assert_eq!(grad, [0.0, 0.0, 0.0]);
    }

    #[test]
    fn matrix3_frobenius_identity_trace() {
        let id = matrix3_identity();
        let fr = matrix3_frobenius(&id, &id);
        assert!((fr - 3.0).abs() < 1e-12, "Frobenius of I:I = 3, got {fr}");
    }

    #[test]
    fn oldroyd_b_polymer_solution_params() {
        let m = OldroydBModel::polymer_solution();
        assert!(m.lambda > 0.0);
        assert!(m.eta_p > 0.0);
        assert!(m.eta_s > 0.0);
    }

    #[test]
    fn fene_p_default_polymer_params() {
        let m = FenePModel::default_polymer();
        assert!(m.b > 0.0);
        assert!(m.lambda > 0.0);
    }

    #[test]
    fn giesekus_alpha_in_valid_range() {
        let m = GiesekusModel::default_polymer();
        assert!(m.alpha >= 0.0 && m.alpha <= 0.5);
    }
}
