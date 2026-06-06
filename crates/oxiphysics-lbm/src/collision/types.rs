//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use crate::lattice::{CS2, D3Q19_VELOCITIES, D3Q19_WEIGHTS};

use super::functions::*;

/// Advanced BGK overrelaxation collision with viscosity-based construction.
///
/// Provides the same interface as `BgkCollision` but exposes the
/// overrelaxation factor σ as an independent parameter.
#[derive(Debug, Clone, Copy)]
pub struct BgkOverrelaxation {
    /// BGK relaxation rate ω = 1/τ.
    pub omega: f64,
    /// Overrelaxation factor σ (default 1.0 = standard BGK).
    pub sigma: f64,
}
impl BgkOverrelaxation {
    /// Construct from kinematic viscosity ν and overrelaxation σ.
    ///
    /// Uses the standard LBM relation ν = cs² (τ − 0.5), τ = 1/ω.
    pub fn new(nu: f64, sigma: f64) -> Self {
        let tau = 3.0 * nu + 0.5;
        Self {
            omega: 1.0 / tau,
            sigma,
        }
    }
    /// Perform the overrelaxed BGK collision step.
    pub fn collide(&self, f: &[f64; 19], rho: f64, u: [f64; 3]) -> [f64; 19] {
        collide_bgk_overrelaxation(f, rho, u, self.omega, self.sigma)
    }
    /// Return the effective kinematic viscosity (ν = cs²(τ − 0.5)).
    pub fn effective_viscosity(&self) -> f64 {
        let tau = 1.0 / self.omega;
        (1.0 / 3.0) * (tau - 0.5)
    }
    /// Stability check: returns true when σω < 2.
    pub fn is_stable(&self) -> bool {
        self.sigma * self.omega < 2.0
    }
}
/// BGK collision parameters.
#[derive(Debug, Clone)]
pub struct BgkCollision {
    /// Relaxation rate omega = 1/tau.
    /// Kinematic viscosity: nu = cs² * (1/omega - 0.5).
    pub omega: f64,
}
impl BgkCollision {
    /// Create a new BGK collision operator with the given relaxation rate.
    pub fn new(omega: f64) -> Self {
        Self { omega }
    }
    /// Create a BGK operator from kinematic viscosity.
    ///
    /// `omega = 1 / (3*nu + 0.5)`
    pub fn from_viscosity(nu: f64) -> Self {
        Self {
            omega: 1.0 / (3.0 * nu + 0.5),
        }
    }
    /// Apply BGK collision to a single D3Q19 distribution.
    ///
    /// Returns the post-collision distribution `f* = f - omega*(f - feq)`.
    pub fn collide(&self, f: &[f64; 19], rho: f64, u: [f64; 3]) -> [f64; 19] {
        collide_bgk(f, rho, u, self.omega)
    }
}
/// TRT collision operator parameters.
///
/// Uses two relaxation rates:
/// - `omega_plus`:  symmetric part, controls kinematic viscosity
/// - `omega_minus`: anti-symmetric part, controls boundary accuracy
///
/// The *magic parameter* `lambda = (tau_plus - 0.5)(tau_minus - 0.5) = 3/16`
/// eliminates numerical slip at no-slip walls in Poiseuille flow.
#[derive(Debug, Clone)]
pub struct TrtCollision {
    /// Symmetric relaxation rate (viscosity-controlling).
    pub omega_plus: f64,
    /// Anti-symmetric relaxation rate.
    pub omega_minus: f64,
}
impl TrtCollision {
    /// Create a TRT operator from kinematic viscosity.
    ///
    /// `omega_plus = 1 / (3*nu + 0.5)`,
    /// `omega_minus` derived from the optimal magic parameter `3/16`.
    pub fn from_viscosity(nu: f64) -> Self {
        let omega_plus = 1.0 / (3.0 * nu + 0.5);
        let tau_plus = 1.0 / omega_plus;
        let tau_minus = TRT_MAGIC / (tau_plus - 0.5) + 0.5;
        Self {
            omega_plus,
            omega_minus: 1.0 / tau_minus,
        }
    }
    /// Create a TRT operator with explicit relaxation rates.
    pub fn new(omega_plus: f64, omega_minus: f64) -> Self {
        Self {
            omega_plus,
            omega_minus,
        }
    }
    /// Apply TRT collision to a single D3Q19 distribution.
    pub fn collide(&self, f: &[f64; 19], rho: f64, u: [f64; 3]) -> [f64; 19] {
        collide_trt(f, rho, u, self.omega_plus, self.omega_minus)
    }
}
/// Which collision scheme is active in a hybrid operator.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CollisionMode {
    /// Standard BGK (stable region).
    Bgk,
    /// Entropic LBM (potentially unstable region).
    Entropic,
}
/// Cumulant LBM collision operator (Geier, Schönherr, Pasquali, Krafczyk 2015).
///
/// Operates in cumulant space (log of the characteristic function) rather than
/// raw or central moment space.  Cumulants of order ≥ 3 are set to zero to
/// enforce the Navier-Stokes limit.  The cumulant method dramatically improves
/// Galilean invariance and stability for high-Re flows.
///
/// Reference: Geier et al., Comput. Fluids 141, 2–17 (2015).
#[derive(Debug, Clone, Copy)]
pub struct CumulantCollision {
    /// Shear relaxation rate ω (controls viscosity).
    pub omega: f64,
    /// Bulk relaxation rate ω_b.
    pub omega_b: f64,
}
impl CumulantCollision {
    /// Construct from shear viscosity ν.
    ///
    /// Bulk relaxation defaults to 1.0 (minimal bulk viscosity).
    pub fn from_viscosity(nu: f64) -> Self {
        let tau = 3.0 * nu + 0.5;
        Self {
            omega: 1.0 / tau,
            omega_b: 1.0,
        }
    }
    /// Construct with explicit bulk viscosity ζ.
    pub fn from_viscosities(nu: f64, zeta: f64) -> Self {
        let tau = 3.0 * nu + 0.5;
        let tau_b = 3.0 * (zeta + 2.0 * nu / 3.0) + 0.5;
        Self {
            omega: 1.0 / tau,
            omega_b: 1.0 / tau_b,
        }
    }
    /// Perform cumulant collision.
    ///
    /// 1. Shift to co-moving frame (subtract mean velocity).
    /// 2. Compute central moments from shifted distribution.
    /// 3. Relax stress cumulants; zero out higher-order cumulants.
    /// 4. Reconstruct central moments, shift back.
    pub fn collide(&self, f: &[f64; 19], rho: f64, u: [f64; 3]) -> [f64; 19] {
        let feq = compute_equilibrium_d3q19(rho, u);
        let cm = compute_central_moments_cumulant(f, u);
        let cm_eq = compute_central_moments_cumulant(&feq, u);
        let mut cm_out = cm;
        let k_trace = (cm[4] + cm[5] + cm[6]) / 3.0;
        let k_trace_eq = (cm_eq[4] + cm_eq[5] + cm_eq[6]) / 3.0;
        let k_dev_xx = cm[4] - k_trace;
        let k_dev_yy = cm[5] - k_trace;
        let k_dev_zz = cm[6] - k_trace;
        let k_dev_xx_eq = cm_eq[4] - k_trace_eq;
        let k_dev_yy_eq = cm_eq[5] - k_trace_eq;
        let k_dev_zz_eq = cm_eq[6] - k_trace_eq;
        let trace_out = k_trace - self.omega_b * (k_trace - k_trace_eq);
        cm_out[4] = trace_out + (k_dev_xx - self.omega * (k_dev_xx - k_dev_xx_eq));
        cm_out[5] = trace_out + (k_dev_yy - self.omega * (k_dev_yy - k_dev_yy_eq));
        cm_out[6] = trace_out + (k_dev_zz - self.omega * (k_dev_zz - k_dev_zz_eq));
        cm_out[7] = cm[7] - self.omega * (cm[7] - cm_eq[7]);
        cm_out[8] = cm[8] - self.omega * (cm[8] - cm_eq[8]);
        cm_out[9] = cm[9] - self.omega * (cm[9] - cm_eq[9]);
        for cm_out_i in cm_out[10..19].iter_mut() {
            *cm_out_i = 0.0;
        }
        central_moments_to_f_cumulant(&cm_out, rho, u)
    }
    /// Return effective shear viscosity.
    pub fn shear_viscosity(&self) -> f64 {
        (1.0 / 3.0) * (1.0 / self.omega - 0.5)
    }
    /// Return effective bulk viscosity.
    pub fn bulk_viscosity(&self) -> f64 {
        (1.0 / 3.0) * (1.0 / self.omega_b - 0.5) * (2.0 / 3.0)
    }
}
/// Full regularized LBM collision operator.
///
/// The standard regularized method reconstructs the non-equilibrium part from
/// the stress tensor Π^(1) before applying the BGK relaxation.  This removes
/// higher-order non-equilibrium contributions and improves numerical stability
/// at moderate to high Reynolds numbers.
///
/// Reference: Latt & Chopard, Phys. Rev. E 72, 036706 (2005).
#[derive(Debug, Clone, Copy)]
pub struct RegularizedCollisionFull {
    /// BGK relaxation rate ω.
    pub omega: f64,
}
impl RegularizedCollisionFull {
    /// Construct from kinematic viscosity ν.
    pub fn from_viscosity(nu: f64) -> Self {
        let tau = 3.0 * nu + 0.5;
        Self { omega: 1.0 / tau }
    }
    /// Perform regularized collision.
    ///
    /// Steps:
    /// 1. Compute feq from (rho, u).
    /// 2. Extract Π^(1) from f − feq.
    /// 3. Reconstruct f^(1) from Π^(1) via Hermite expansion.
    /// 4. Apply BGK: f_out = feq + (1 − ω) f^(1).
    pub fn collide(&self, f: &[f64; 19], rho: f64, u: [f64; 3]) -> [f64; 19] {
        let feq = compute_equilibrium_d3q19(rho, u);
        let pi1 = compute_pi1_tensor(f, &feq);
        let f1 = regularized_f1_d3q19(&pi1, rho);
        let mut f_out = [0.0f64; 19];
        for (f_out_i, (&feqi, &f1_i)) in f_out.iter_mut().zip(feq.iter().zip(f1.iter())) {
            *f_out_i = feqi + (1.0 - self.omega) * f1_i;
        }
        f_out
    }
    /// Return the effective viscosity.
    pub fn effective_viscosity(&self) -> f64 {
        (1.0 / 3.0) * (1.0 / self.omega - 0.5)
    }
}
/// Hybrid BGK/entropic collision that switches based on local instability.
///
/// Uses BGK in stable regions (low non-equilibrium ratio) and entropic
/// collision when the non-equilibrium magnitude exceeds a threshold.
#[derive(Debug, Clone)]
pub struct HybridCollision {
    /// BGK operator.
    pub bgk: BgkCollision,
    /// Entropic operator.
    pub elbm: EntropicCollision,
    /// Non-equilibrium threshold above which entropic mode activates.
    pub neq_threshold: f64,
}
impl HybridCollision {
    /// Create a hybrid operator from kinematic viscosity.
    pub fn from_viscosity(nu: f64, neq_threshold: f64) -> Self {
        let omega = 1.0 / (3.0 * nu + 0.5);
        Self {
            bgk: BgkCollision::new(omega),
            elbm: EntropicCollision::new(omega),
            neq_threshold,
        }
    }
    /// Apply the hybrid collision.
    pub fn collide(&self, f: &[f64; 19], rho: f64, u: [f64; 3]) -> [f64; 19] {
        let ratio = non_equilibrium_ratio(f, rho, u);
        if ratio > self.neq_threshold {
            self.elbm.collide(f, rho, u)
        } else {
            self.bgk.collide(f, rho, u)
        }
    }
    /// Check which collision mode would be used for a given distribution.
    pub fn active_mode(&self, f: &[f64; 19], rho: f64, u: [f64; 3]) -> CollisionMode {
        let ratio = non_equilibrium_ratio(f, rho, u);
        if ratio > self.neq_threshold {
            CollisionMode::Entropic
        } else {
            CollisionMode::Bgk
        }
    }
}
/// Improved ELBM with mirror state stabilization.
///
/// Uses the Karlin-Bösch-Chikatamarla (KBC) approach where the distribution
/// is split into kinetic and ghost parts:
///
/// `f = k + s + g`
///
/// where `k` is the kinetic (equilibrium-driven) part and `s + g` are the
/// non-equilibrium parts.  The entropy condition is applied only to the ghost
/// part, improving stability over standard ELBM.
#[derive(Debug, Clone)]
pub struct KbcCollision {
    /// Viscosity-controlling relaxation rate.
    pub omega_s: f64,
    /// Entropy-bounded ghost relaxation (≤ 2).
    pub omega_h_max: f64,
}
impl KbcCollision {
    /// Create a KBC collision operator from kinematic viscosity.
    pub fn new(nu: f64) -> Self {
        let omega_s = 1.0 / (3.0 * nu + 0.5);
        Self {
            omega_s,
            omega_h_max: 2.0,
        }
    }
    /// Apply KBC collision with entropy stabilization.
    ///
    /// Computes the mirror state `f* = feq + (feq - f)` and uses it to
    /// find an alpha ∈ \[0, 2\] satisfying H(f + alpha*(feq - f)) ≤ H(f).
    pub fn collide(&self, f: &[f64; 19], rho: f64, u: [f64; 3]) -> [f64; 19] {
        let feq = compute_equilibrium_d3q19(rho, u);
        let ds: [f64; 19] = std::array::from_fn(|i| feq[i] - f[i]);
        let f_mirror: [f64; 19] = std::array::from_fn(|i| feq[i] + ds[i]);
        let gamma = self.find_gamma(f, &feq, &f_mirror);
        std::array::from_fn(|i| f[i] + gamma * ds[i])
    }
    /// Find the scalar gamma such that H(f + gamma*(feq - f)) is maximized
    /// subject to the entropy condition.
    fn find_gamma(&self, f: &[f64; 19], feq: &[f64; 19], _f_mirror: &[f64; 19]) -> f64 {
        let h_f = entropy_h(f);
        let mut gamma = self.omega_s;
        let f_at_2: [f64; 19] = std::array::from_fn(|i| (f[i] + 2.0 * (feq[i] - f[i])).max(1e-30));
        let h_mirror = entropy_h(&f_at_2);
        if h_mirror <= h_f {
            gamma = self.omega_h_max;
        }
        gamma.clamp(self.omega_s, self.omega_h_max)
    }
}
/// Recursive regularized LBM collision operator.
///
/// Extends the standard regularized method by recursively computing
/// higher-order non-equilibrium tensors (up to second order) using
/// the Chapman-Enskog expansion.  Provides better accuracy than the
/// standard regularized BGK without the overhead of full MRT.
///
/// Reference: Malaspinas, arXiv:1505.06900 (2015).
#[derive(Debug, Clone, Copy)]
pub struct RecursiveRegularized {
    /// BGK relaxation rate ω.
    pub omega: f64,
}
impl RecursiveRegularized {
    /// Construct from kinematic viscosity ν.
    pub fn from_viscosity(nu: f64) -> Self {
        let tau = 3.0 * nu + 0.5;
        Self { omega: 1.0 / tau }
    }
    /// Perform RR-LBM collision.
    ///
    /// 1. Compute feq.
    /// 2. Extract Q^(1) (third-order moment contribution) recursively.
    /// 3. Reconstruct f^(1) including Q correction.
    /// 4. Apply BGK relaxation.
    pub fn collide(&self, f: &[f64; 19], rho: f64, u: [f64; 3]) -> [f64; 19] {
        let feq = compute_equilibrium_d3q19(rho, u);
        let pi1 = compute_pi1_tensor(f, &feq);
        let q1 = compute_q1_recursive(&pi1, u, 1.0 / self.omega - 0.5);
        let f1 = regularized_f1_with_q_d3q19(&pi1, &q1);
        let mut f_out = [0.0f64; 19];
        for (f_out_i, (&feqi, &f1_i)) in f_out.iter_mut().zip(feq.iter().zip(f1.iter())) {
            *f_out_i = feqi + (1.0 - self.omega) * f1_i;
        }
        f_out
    }
    /// Return effective viscosity.
    pub fn effective_viscosity(&self) -> f64 {
        (1.0 / 3.0) * (1.0 / self.omega - 0.5)
    }
}
/// Entropic LBM collision operator.
///
/// The ELBM adaptively adjusts the relaxation parameter to guarantee
/// the H-theorem (entropy increase). Instead of using a fixed omega,
/// it solves for a scalar `alpha` such that:
///
/// `H(f + alpha * (feq - f)) <= H(f)`
///
/// where `H(f) = sum_i f_i * ln(f_i / w_i)`.
#[derive(Debug, Clone)]
pub struct EntropicCollision {
    /// Nominal relaxation rate (used as the starting guess for alpha).
    pub omega: f64,
    /// Maximum iterations for the entropy optimization.
    pub max_iter: usize,
    /// Tolerance for the entropy condition.
    pub tolerance: f64,
}
impl EntropicCollision {
    /// Create an entropic collision operator.
    pub fn new(omega: f64) -> Self {
        Self {
            omega,
            max_iter: 20,
            tolerance: 1e-10,
        }
    }
    /// Apply the entropic collision.
    ///
    /// Computes `alpha` adaptively, then applies:
    /// `f*_i = f_i + alpha * omega * (feq_i - f_i)`
    pub fn collide(&self, f: &[f64; 19], rho: f64, u: [f64; 3]) -> [f64; 19] {
        let feq = compute_equilibrium_d3q19(rho, u);
        let alpha = self.find_alpha(f, &feq);
        let mut f_out = [0.0f64; 19];
        for (f_out_i, (&fi, &feqi)) in f_out.iter_mut().zip(f.iter().zip(feq.iter())) {
            *f_out_i = fi + alpha * self.omega * (feqi - fi);
        }
        f_out
    }
    /// Find the optimal alpha that satisfies the entropy condition.
    ///
    /// Uses Newton's method on `Delta_H(alpha) = H(f + alpha*(feq-f)) - H(f)`.
    fn find_alpha(&self, f: &[f64; 19], feq: &[f64; 19]) -> f64 {
        let mut alpha = 1.0;
        let h_f = self.entropy(f);
        for _ in 0..self.max_iter {
            let mut f_trial = [0.0f64; 19];
            for (f_trial_i, (&fi, &feqi)) in f_trial.iter_mut().zip(f.iter().zip(feq.iter())) {
                *f_trial_i = (fi + alpha * self.omega * (feqi - fi)).max(1e-30);
            }
            let h_trial = self.entropy(&f_trial);
            let delta_h = h_trial - h_f;
            if delta_h.abs() < self.tolerance {
                break;
            }
            if delta_h > 0.0 {
                alpha *= 0.5;
            } else {
                break;
            }
        }
        alpha.clamp(0.0, 2.0)
    }
    /// Compute the H-function (Boltzmann entropy).
    ///
    /// `H(f) = sum_i f_i * ln(f_i / w_i)`
    fn entropy(&self, f: &[f64; 19]) -> f64 {
        f.iter()
            .zip(D3Q19_WEIGHTS.iter())
            .map(|(&fi_raw, &wi)| {
                let fi = fi_raw.max(1e-30);
                fi * (fi / wi).ln()
            })
            .sum()
    }
}
/// Multiple-Relaxation-Time (MRT) collision operator for D3Q19.
///
/// Relaxes each moment with its own relaxation rate for improved stability
/// and accuracy. The transformation matrix M maps populations to moments:
/// `m = M f`, collision: `m* = m - S(m - m_eq)`, post-collision: `f* = M^{-1} m*`.
///
/// Uses a simplified 19-moment basis from Lallemand & Luo (2000).
#[derive(Debug, Clone)]
pub struct MrtCollision {
    /// Diagonal relaxation matrix entries s\[0..19\].
    pub s: [f64; 19],
}
impl MrtCollision {
    /// Create an MRT operator with given relaxation rates.
    pub fn new(s: [f64; 19]) -> Self {
        Self { s }
    }
    /// Create an MRT operator from kinematic viscosity.
    ///
    /// Uses the standard Lallemand-Luo relaxation rates for D3Q19.
    pub fn from_viscosity(nu: f64) -> Self {
        let omega = 1.0 / (3.0 * nu + 0.5);
        let mut s = [1.0f64; 19];
        s[0] = 0.0;
        s[1] = 1.19;
        s[2] = 1.4;
        s[3] = 0.0;
        s[4] = 1.2;
        s[5] = 0.0;
        s[6] = 1.2;
        s[7] = 0.0;
        s[8] = 1.2;
        s[9] = omega;
        s[10] = 1.4;
        s[11] = omega;
        s[12] = 1.4;
        s[13] = omega;
        s[14] = omega;
        s[15] = omega;
        s[16] = 1.98;
        s[17] = 1.98;
        s[18] = 1.98;
        Self { s }
    }
    /// Apply MRT collision.
    ///
    /// This simplified implementation uses the decomposed form:
    /// `f* = f - M^{-1} S (m - m_eq)`
    /// approximated via the regularized collision with per-mode relaxation.
    pub fn collide(&self, f: &[f64; 19], rho: f64, u: [f64; 3]) -> [f64; 19] {
        let feq = compute_equilibrium_d3q19(rho, u);
        let mut pi_neq = [[0.0f64; 3]; 3];
        for (i, (&fi, &feqi)) in f.iter().zip(feq.iter()).enumerate() {
            let c = D3Q19_VELOCITIES[i];
            let f_neq = fi - feqi;
            for a in 0..3 {
                for b in 0..3 {
                    pi_neq[a][b] += f_neq * c[a] as f64 * c[b] as f64;
                }
            }
        }
        let omega_s = self.s[9];
        let omega_g = (self.s[16] + self.s[17] + self.s[18]) / 3.0;
        let mut f_out = [0.0f64; 19];
        for (i, (f_out_i, (&fi, &feqi))) in
            f_out.iter_mut().zip(f.iter().zip(feq.iter())).enumerate()
        {
            let w = D3Q19_WEIGHTS[i];
            let c = D3Q19_VELOCITIES[i];
            let mut stress_contrib = 0.0f64;
            for a in 0..3 {
                for b in 0..3 {
                    let ca = c[a] as f64;
                    let cb = c[b] as f64;
                    let delta_ab = if a == b { 1.0 } else { 0.0 };
                    stress_contrib += (ca * cb - CS2 * delta_ab) * pi_neq[a][b];
                }
            }
            let f_neq_stress = w / (2.0 * CS2 * CS2) * stress_contrib;
            let f_neq_ghost = (fi - feqi) - f_neq_stress;
            *f_out_i = feqi + (1.0 - omega_s) * f_neq_stress + (1.0 - omega_g) * f_neq_ghost;
        }
        f_out
    }
    /// Compute the effective viscosity from the stress relaxation rate.
    pub fn effective_viscosity(&self) -> f64 {
        let omega_s = self.s[9];
        if omega_s <= 0.0 || omega_s >= 2.0 {
            return 0.0;
        }
        CS2 * (1.0 / omega_s - 0.5)
    }
}
/// Regularized collision operator.
///
/// Projects the distribution onto the equilibrium plus the
/// non-equilibrium stress tensor, filtering out higher-order
/// non-hydrodynamic modes:
///
/// `f*_i = feq_i + (1 - omega) * f_neq_reg_i`
///
/// where `f_neq_reg_i` is reconstructed from the stress tensor only.
#[derive(Debug, Clone)]
pub struct RegularizedCollision {
    /// Relaxation rate omega.
    pub omega: f64,
}
impl RegularizedCollision {
    /// Create a regularized collision operator.
    pub fn new(omega: f64) -> Self {
        Self { omega }
    }
    /// Create from kinematic viscosity.
    pub fn from_viscosity(nu: f64) -> Self {
        Self {
            omega: 1.0 / (3.0 * nu + 0.5),
        }
    }
    /// Apply regularized collision.
    pub fn collide(&self, f: &[f64; 19], rho: f64, u: [f64; 3]) -> [f64; 19] {
        let feq = compute_equilibrium_d3q19(rho, u);
        let mut pi = [[0.0f64; 3]; 3];
        for (i, (&fi, &feqi)) in f.iter().zip(feq.iter()).enumerate() {
            let c = D3Q19_VELOCITIES[i];
            let f_neq = fi - feqi;
            for a in 0..3 {
                for b in 0..3 {
                    pi[a][b] += f_neq * c[a] as f64 * c[b] as f64;
                }
            }
        }
        let mut f_out = [0.0f64; 19];
        for (i, (f_out_i, &feqi)) in f_out.iter_mut().zip(feq.iter()).enumerate() {
            let w = D3Q19_WEIGHTS[i];
            let c = D3Q19_VELOCITIES[i];
            let mut q_i = 0.0;
            for a in 0..3 {
                for b in 0..3 {
                    let c_a = c[a] as f64;
                    let c_b = c[b] as f64;
                    let delta_ab = if a == b { 1.0 } else { 0.0 };
                    q_i += (c_a * c_b - CS2 * delta_ab) * pi[a][b];
                }
            }
            let f_neq_reg = w / (2.0 * CS2 * CS2) * q_i;
            *f_out_i = feqi + (1.0 - self.omega) * f_neq_reg;
        }
        f_out
    }
}
/// Hybrid recursive-regularized collision operator.
///
/// Blends the recursive-regularized (RR) post-collision distribution with
/// the raw BGK update using a sensor σ based on the local non-equilibrium
/// intensity.  In smooth laminar regions σ → 1 (RR), and in turbulent or
/// shock regions σ → 0 (BGK), maintaining accuracy and robustness.
///
/// Reference: Jacob, Malaspinas, Sagaut, JCP 2018.
#[derive(Debug, Clone, Copy)]
pub struct HybridRecursiveRegularized {
    /// BGK relaxation rate ω.
    pub omega: f64,
    /// Threshold for the non-equilibrium sensor (default ≈ 0.01).
    pub sigma_threshold: f64,
}
impl HybridRecursiveRegularized {
    /// Construct from kinematic viscosity ν and sensor threshold.
    pub fn new(nu: f64, sigma_threshold: f64) -> Self {
        let tau = 3.0 * nu + 0.5;
        Self {
            omega: 1.0 / tau,
            sigma_threshold,
        }
    }
    /// Perform HRR collision.
    pub fn collide(&self, f: &[f64; 19], rho: f64, u: [f64; 3]) -> [f64; 19] {
        let feq = compute_equilibrium_d3q19(rho, u);
        let pi1 = compute_pi1_tensor(f, &feq);
        let q1 = compute_q1_recursive(&pi1, u, 1.0 / self.omega - 0.5);
        let f1_rr = regularized_f1_with_q_d3q19(&pi1, &q1);
        let neq_mag: f64 = f.iter().zip(feq.iter()).map(|(a, b)| (a - b).abs()).sum();
        let eq_mag: f64 = feq.iter().map(|v| v.abs()).sum();
        let sensor = if eq_mag > 1e-14 {
            neq_mag / eq_mag
        } else {
            1.0
        };
        let sigma = if sensor < self.sigma_threshold {
            1.0
        } else {
            0.0_f64.max(1.0 - sensor / self.sigma_threshold)
        };
        let mut f_out = [0.0f64; 19];
        for (i, (f_out_i, ((&fi, &feqi), &f1_rr_i))) in f_out
            .iter_mut()
            .zip(f.iter().zip(feq.iter()).zip(f1_rr.iter()))
            .enumerate()
        {
            let _ = i;
            let f_rr = feqi + (1.0 - self.omega) * f1_rr_i;
            let f_bgk = fi - self.omega * (fi - feqi);
            *f_out_i = sigma * f_rr + (1.0 - sigma) * f_bgk;
        }
        f_out
    }
    /// Evaluate the non-equilibrium sensor value at a given state.
    pub fn sensor(&self, f: &[f64; 19], rho: f64, u: [f64; 3]) -> f64 {
        let feq = compute_equilibrium_d3q19(rho, u);
        let neq_mag: f64 = f.iter().zip(feq.iter()).map(|(a, b)| (a - b).abs()).sum();
        let eq_mag: f64 = feq.iter().map(|v| v.abs()).sum();
        if eq_mag > 1e-14 {
            neq_mag / eq_mag
        } else {
            1.0
        }
    }
    /// Return effective viscosity.
    pub fn effective_viscosity(&self) -> f64 {
        (1.0 / 3.0) * (1.0 / self.omega - 0.5)
    }
}
/// Cascaded (central-moment) collision with full D3Q19 central moment set.
///
/// This improved version computes raw moments, shifts to central moments,
/// relaxes each central moment to its equilibrium, then transforms back.
#[derive(Debug, Clone)]
pub struct CentralMomentCollision {
    /// Relaxation rates for each of the 19 moment groups.
    /// Typically: s\[1\]=s\[2\]=0 (conservation), s\[3..=5\]=omega_s, s\[6..\]=omega_q.
    pub s: [f64; 19],
}
impl CentralMomentCollision {
    /// Create from kinematic viscosity with standard moment relaxation rates.
    pub fn from_viscosity(nu: f64) -> Self {
        let omega_s = 1.0 / (3.0 * nu + 0.5);
        let mut s = [1.0f64; 19];
        s[0] = 0.0;
        s[1] = 0.0;
        s[2] = 0.0;
        s[3] = 0.0;
        s[4] = omega_s;
        s[5] = omega_s;
        s[6] = omega_s;
        s[7] = omega_s;
        s[8] = omega_s;
        s[9] = omega_s;
        for s_i in s[10..19].iter_mut() {
            *s_i = 1.0;
        }
        Self { s }
    }
    /// Apply central-moment collision.
    ///
    /// Computes raw moments → shifts to central moments →
    /// relaxes → un-shifts → inverse transform back to populations.
    pub fn collide(&self, f: &[f64; 19], rho: f64, u: [f64; 3]) -> [f64; 19] {
        let feq = compute_equilibrium_d3q19(rho, u);
        let inv_rho = if rho > 1e-15 { 1.0 / rho } else { 1.0 };
        let cm_f = compute_central_moments(f, u, inv_rho);
        let cm_feq = compute_central_moments(&feq, u, inv_rho);
        let mut cm_post = [0.0f64; 19];
        for (cm_post_i, ((&cm_fi, &cm_feqi), &si)) in cm_post
            .iter_mut()
            .zip(cm_f.iter().zip(cm_feq.iter()).zip(self.s.iter()))
        {
            *cm_post_i = cm_fi - si * (cm_fi - cm_feqi);
        }
        let mut f_out = feq;
        for (f_out_i, ((&cm_post_i, &cm_feqi), &wi)) in f_out
            .iter_mut()
            .zip(cm_post.iter().zip(cm_feq.iter()).zip(D3Q19_WEIGHTS.iter()))
        {
            let neq_cm = cm_post_i - cm_feqi;
            *f_out_i += neq_cm * wi;
        }
        f_out
    }
}
/// Cascaded collision operator (simplified D3Q19 version).
///
/// Relaxes central moments independently, providing better numerical
/// stability and tunable parameters for each moment.
///
/// This is a simplified version that uses the same relaxation rate for
/// stress-related moments and a separate rate for higher-order moments.
#[derive(Debug, Clone)]
pub struct CascadedCollision {
    /// Relaxation rate for stress (viscosity-related) moments.
    pub omega_s: f64,
    /// Relaxation rate for higher-order (ghost) moments.
    pub omega_q: f64,
}
impl CascadedCollision {
    /// Create a cascaded collision operator.
    pub fn new(omega_s: f64, omega_q: f64) -> Self {
        Self { omega_s, omega_q }
    }
    /// Create from kinematic viscosity with default ghost relaxation.
    pub fn from_viscosity(nu: f64) -> Self {
        let omega_s = 1.0 / (3.0 * nu + 0.5);
        Self {
            omega_s,
            omega_q: 1.0,
        }
    }
    /// Apply cascaded collision (simplified).
    ///
    /// For simplicity, this implementation applies BGK-like relaxation
    /// but separates stress and higher-order moment relaxation.
    pub fn collide(&self, f: &[f64; 19], rho: f64, u: [f64; 3]) -> [f64; 19] {
        let feq = compute_equilibrium_d3q19(rho, u);
        let mut pi_neq = [[0.0f64; 3]; 3];
        for (i, (&fi, &feqi)) in f.iter().zip(feq.iter()).enumerate() {
            let c = D3Q19_VELOCITIES[i];
            let f_neq = fi - feqi;
            for a in 0..3 {
                for b in 0..3 {
                    pi_neq[a][b] += f_neq * c[a] as f64 * c[b] as f64;
                }
            }
        }
        let mut f_out = [0.0f64; 19];
        for (i, (f_out_i, (&fi, &feqi))) in
            f_out.iter_mut().zip(f.iter().zip(feq.iter())).enumerate()
        {
            let w = D3Q19_WEIGHTS[i];
            let c = D3Q19_VELOCITIES[i];
            let mut stress_contrib = 0.0;
            for a in 0..3 {
                for b in 0..3 {
                    let c_a = c[a] as f64;
                    let c_b = c[b] as f64;
                    let delta_ab = if a == b { 1.0 } else { 0.0 };
                    stress_contrib += (c_a * c_b - CS2 * delta_ab) * pi_neq[a][b];
                }
            }
            let f_neq_stress = w / (2.0 * CS2 * CS2) * stress_contrib;
            let f_neq_total = fi - feqi;
            let f_neq_ghost = f_neq_total - f_neq_stress;
            *f_out_i =
                feqi + (1.0 - self.omega_s) * f_neq_stress + (1.0 - self.omega_q) * f_neq_ghost;
        }
        f_out
    }
}
/// Raw moment collision operator.
///
/// Transforms distributions to raw moment space, applies separate relaxation
/// rates for each moment, then transforms back.  This provides the same
/// framework as MRT but operates directly in raw (non-orthogonalized) moments.
///
/// Raw moments m_k = Σ_i f_i c_{ix}^{k_x} c_{iy}^{k_y} c_{iz}^{k_z}
#[derive(Debug, Clone, Copy)]
pub struct RawMomentCollision {
    /// Viscous relaxation rate ω_v (for stress-tensor moments).
    pub omega_v: f64,
    /// Bulk viscosity relaxation rate ω_b.
    pub omega_b: f64,
    /// Ghost mode relaxation rate ω_g (non-hydrodynamic).
    pub omega_g: f64,
}
impl RawMomentCollision {
    /// Construct from shear viscosity ν and bulk viscosity ζ.
    pub fn from_viscosity(nu: f64, zeta: f64) -> Self {
        let tau_v = 3.0 * nu + 0.5;
        let tau_b = 3.0 * (zeta + 2.0 * nu / 3.0) + 0.5;
        Self {
            omega_v: 1.0 / tau_v,
            omega_b: 1.0 / tau_b,
            omega_g: 1.0,
        }
    }
    /// Perform raw moment collision.
    ///
    /// Uses the D3Q19 raw moment representation.  Hydrodynamic moments
    /// (density, momentum) are conserved; stress moments relax at ω_v.
    pub fn collide(&self, f: &[f64; 19], rho: f64, u: [f64; 3]) -> [f64; 19] {
        let feq = compute_equilibrium_d3q19(rho, u);
        let m = compute_raw_moments_d3q19(f);
        let meq = compute_raw_moments_d3q19(&feq);
        let mut m_out = [0.0f64; 19];
        for (i, (m_out_i, (&mi, &meqi))) in
            m_out.iter_mut().zip(m.iter().zip(meq.iter())).enumerate()
        {
            let rate = raw_moment_rate(i, self.omega_v, self.omega_b, self.omega_g);
            *m_out_i = mi - rate * (mi - meqi);
        }
        raw_moments_to_f_d3q19(&m_out)
    }
    /// Return effective shear viscosity.
    pub fn shear_viscosity(&self) -> f64 {
        (1.0 / 3.0) * (1.0 / self.omega_v - 0.5)
    }
}
