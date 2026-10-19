//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use crate::lattice::{CS2, D3Q19_VELOCITIES, D3Q19_WEIGHTS};

use super::central_moments::*;
use super::functions::*;
use super::kbc::{kbc_decompose_d3q19, kbc_stabilizer};

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
/// Cumulant LBM collision operator on the D3Q19 lattice
/// (after Geier, Schönherr, Pasquali, Krafczyk 2015).
///
/// The populations are mapped to the 19 central moments about the local velocity (the central
/// moments that exist on D3Q19, see [`central_moments`](super::central_moments)) and from there
/// to cumulants, i.e. to the moments of the logarithm of the characteristic function. The
/// cumulants relax towards their equilibrium values:
///
/// - the second-order cumulants relax with `omega` (shear stress and normal-stress differences)
///   and with `omega_b` (trace, i.e. bulk viscosity),
/// - the third-order cumulants and the three fourth-order cumulants `C_220, C_202, C_022`
///   relax with rate 1, i.e. they are reset to their equilibrium values,
/// - density and momentum are conserved.
///
/// The fourth-order cumulants are the central moments minus products of the second-order
/// moments, `C_220 = kappa_220 - (kappa_200 kappa_020 + 2 kappa_110^2) / rho`. Resetting the
/// cumulant therefore makes the post-collision fourth-order central moment follow the
/// *relaxed* second-order moments, which is what gives the cumulant operator its improved
/// Galilean invariance compared with relaxing the central moments independently.
///
/// The 27-velocity operator of Geier et al. also relaxes the cumulants that involve all three
/// axes. On D3Q19 (no corner velocities) those are not independent quantities, so this operator
/// is the restriction of the cumulant method to the 19 representable moments. The equilibrium
/// is an exact fixed point, and with `omega == omega_b == 1` the output is the equilibrium.
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
    ///
    /// The relaxation times are `τ = 3ν + 1/2` and `τ_b = 9ζ/2 + 1/2`, so that
    /// [`shear_viscosity`](Self::shear_viscosity) returns `ν` and
    /// [`bulk_viscosity`](Self::bulk_viscosity) returns `ζ`.
    pub fn from_viscosities(nu: f64, zeta: f64) -> Self {
        let tau = 3.0 * nu + 0.5;
        let tau_b = 4.5 * zeta + 0.5;
        Self {
            omega: 1.0 / tau,
            omega_b: 1.0 / tau_b,
        }
    }
    /// Perform cumulant collision.
    ///
    /// `rho` and `u` must be the density and velocity of `f`.
    ///
    /// 1. Compute the central moments of `f` and of its equilibrium in the frame moving with `u`.
    /// 2. Convert them to cumulants.
    /// 3. Relax the second-order cumulants (trace with `omega_b`, trace-free part and shear with
    ///    `omega`) and reset the third- and fourth-order cumulants to their equilibrium values.
    /// 4. Convert back to central moments and reconstruct the populations with the exact
    ///    inverse transform, so mass and momentum are conserved to round-off.
    pub fn collide(&self, f: &[f64; 19], rho: f64, u: [f64; 3]) -> [f64; 19] {
        let feq = compute_equilibrium_d3q19(rho, u);
        let cumulants = central_moments_to_cumulants_d3q19(&compute_central_moments_cumulant(f, u));
        let cumulants_eq =
            central_moments_to_cumulants_d3q19(&compute_central_moments_cumulant(&feq, u));
        let mut relaxed = cumulants;

        let trace = (cumulants[4] + cumulants[5] + cumulants[6]) / 3.0;
        let trace_eq = (cumulants_eq[4] + cumulants_eq[5] + cumulants_eq[6]) / 3.0;
        let trace_out = trace - self.omega_b * (trace - trace_eq);
        for k in 4..7 {
            let deviator = cumulants[k] - trace;
            let deviator_eq = cumulants_eq[k] - trace_eq;
            relaxed[k] = trace_out + deviator - self.omega * (deviator - deviator_eq);
        }
        for k in 7..10 {
            relaxed[k] = cumulants[k] - self.omega * (cumulants[k] - cumulants_eq[k]);
        }
        relaxed[10..19].copy_from_slice(&cumulants_eq[10..19]);

        central_moments_to_populations_d3q19(&cumulants_to_central_moments_d3q19(&relaxed), u)
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
/// Karlin-Bösch-Chikatamarla (KBC) entropic collision on D3Q19.
///
/// The non-equilibrium part `Δf = f − feq` is split into a *shear* part `Δs` and a
/// *higher-order* part `Δh = Δf − Δs`:
///
/// - `Δs_i = w_i / (2 cs⁴) (c_i c_i) : Π_dev` carries the traceless part of the
///   non-equilibrium stress `Π = Σ_i Δf_i c_i c_i`,
/// - `Δh` carries the trace of `Π` (the bulk mode) and all moments of higher order.
///
/// Both parts have no mass and no momentum. The collision relaxes them with different rates,
///
/// `f* = f − β (2 Δs + γ Δh)`, `β = ω_s / 2`,
///
/// so the shear stress always relaxes with `ω_s = 2β` (the kinematic viscosity
/// `ν = cs² (1/ω_s − 1/2)`), whereas the stabiliser `γ` is chosen from the entropy estimate
/// `γ = 1/β − (2 − 1/β) ⟨Δs|Δh⟩ / ⟨Δh|Δh⟩` with the entropic scalar product
/// `⟨x|y⟩ = Σ_i x_i y_i / feq_i` (it extremises the quadratic entropy estimate
/// `⟨f* − feq|f* − feq⟩`). For `γ = 2` the operator is BGK with rate `ω_s`; it is also BGK for
/// every state when `ω_s = 1`, and whenever `Δh = 0`. Because the trace of `Π` belongs to `Δh`,
/// the bulk mode relaxes with the entropic rate `β γ` instead of a prescribed bulk viscosity.
///
/// Reference: Karlin, Bösch, Chikatamarla, Phys. Rev. E 90, 031302(R) (2014); Bösch, Chikatamarla,
/// Karlin, Phys. Rev. E 92, 043309 (2015).
#[derive(Debug, Clone)]
pub struct KbcCollision {
    /// Shear relaxation rate `ω_s = 2β`, which sets the kinematic viscosity.
    pub omega_s: f64,
}
impl KbcCollision {
    /// Create a KBC collision operator from kinematic viscosity.
    ///
    /// `ω_s = 1 / (3 ν + 1/2)`.
    pub fn new(nu: f64) -> Self {
        let omega_s = 1.0 / (3.0 * nu + 0.5);
        Self { omega_s }
    }
    /// Kinematic viscosity `ν = cs² (1/ω_s − 1/2)` of the shear relaxation.
    pub fn effective_viscosity(&self) -> f64 {
        (1.0 / 3.0) * (1.0 / self.omega_s - 0.5)
    }
    /// Entropic stabiliser `γ` for the state `f` with density `rho` and velocity `u`.
    ///
    /// `γ = 2` reproduces BGK; see the type documentation for the definition.
    pub fn stabilizer(&self, f: &[f64; 19], rho: f64, u: [f64; 3]) -> f64 {
        let feq = compute_equilibrium_d3q19(rho, u);
        let (delta_s, delta_h) = kbc_decompose_d3q19(f, &feq);
        kbc_stabilizer(0.5 * self.omega_s, &delta_s, &delta_h, &feq)
    }
    /// Apply the KBC collision `f* = f − β (2 Δs + γ Δh)`.
    ///
    /// `rho` and `u` must be the density and velocity of `f`; mass and momentum are then
    /// conserved to round-off for every `γ`.
    pub fn collide(&self, f: &[f64; 19], rho: f64, u: [f64; 3]) -> [f64; 19] {
        let feq = compute_equilibrium_d3q19(rho, u);
        let (delta_s, delta_h) = kbc_decompose_d3q19(f, &feq);
        let beta = 0.5 * self.omega_s;
        let gamma = kbc_stabilizer(beta, &delta_s, &delta_h, &feq);
        std::array::from_fn(|i| f[i] - beta * (2.0 * delta_s[i] + gamma * delta_h[i]))
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
/// Cascaded (central-moment) collision with the full D3Q19 central moment set.
///
/// The populations are mapped to the 19 central moments about the local velocity, each moment
/// relaxes towards its equilibrium value with its own rate, and the populations are recovered
/// with the exact inverse of the central-moment transform (see
/// [`central_moments`](super::central_moments)). Mass and momentum are conserved to round-off
/// and a uniform rate reproduces BGK collision.
#[derive(Debug, Clone)]
pub struct CentralMomentCollision {
    /// Relaxation rates of the 19 central moments, in the order of
    /// [`MOMENT_EXPONENTS`]: `s\[0\]` density, `s\[1..=3\]` momentum
    /// (both stay at their equilibrium for consistent input), `s\[4..=6\]` the normal second-order
    /// moments `kappa_200, kappa_020, kappa_002`, `s\[7..=9\]` the shear moments
    /// `kappa_110, kappa_101, kappa_011`, `s\[10..=15\]` the third-order moments and
    /// `s\[16..=18\]` the fourth-order moments `kappa_220, kappa_202, kappa_022`.
    pub s: [f64; 19],
}
impl CentralMomentCollision {
    /// Create from kinematic viscosity with standard moment relaxation rates.
    ///
    /// The second-order moments relax with `omega = 1 / (3 nu + 1/2)`; all higher-order moments
    /// relax with rate 1 (they are reset to their equilibrium values).
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
    /// Computes the central moments of `f` and of its equilibrium about `u`, relaxes
    /// `kappa* = kappa - s (kappa - kappa_eq)` and transforms back to populations with the exact
    /// inverse (binomial shift to raw moments, then the closed-form raw-moment inverse).
    /// `rho` and `u` must be the density and velocity of `f`.
    pub fn collide(&self, f: &[f64; 19], rho: f64, u: [f64; 3]) -> [f64; 19] {
        let feq = compute_equilibrium_d3q19(rho, u);
        let cm_f = compute_central_moments_cumulant(f, u);
        let cm_feq = compute_central_moments_cumulant(&feq, u);
        let mut cm_post = [0.0f64; 19];
        for (cm_post_i, ((&cm_fi, &cm_feqi), &si)) in cm_post
            .iter_mut()
            .zip(cm_f.iter().zip(cm_feq.iter()).zip(self.s.iter()))
        {
            *cm_post_i = cm_fi - si * (cm_fi - cm_feqi);
        }
        central_moments_to_populations_d3q19(&cm_post, u)
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
/// rates for each moment group, then transforms back with the exact inverse of the 19×19 raw
/// moment matrix (see [`raw_moments_to_f_d3q19`]). This provides the same
/// framework as MRT but operates directly in raw (non-orthogonalized) moments.
///
/// Raw moments m_k = Σ_i f_i c_{ix}^{k_x} c_{iy}^{k_y} c_{iz}^{k_z}
///
/// The moment groups are updated as follows:
///
/// - density and momentum (`m_000`, `m_100`, `m_010`, `m_001`) are conserved (rate 1 towards
///   the equilibrium of the supplied `rho`, `u`),
/// - the trace `m_200 + m_020 + m_002` relaxes with the bulk rate `ω_b`,
/// - the trace-free normal stresses `m_aa − trace / 3` and the shear stresses
///   `m_110, m_101, m_011` relax with the shear rate `ω_v`,
/// - the third- and fourth-order moments (ghosts, `m_210 … m_022`) relax with `ω_g`.
///
/// With `ω_v = ω_b = ω_g = ω` the operator is exactly BGK with rate `ω`.
#[derive(Debug, Clone, Copy)]
pub struct RawMomentCollision {
    /// Viscous relaxation rate ω_v (for the shear and trace-free normal stresses).
    pub omega_v: f64,
    /// Bulk viscosity relaxation rate ω_b (for the trace of the stress tensor).
    pub omega_b: f64,
    /// Ghost mode relaxation rate ω_g (non-hydrodynamic).
    pub omega_g: f64,
}
impl RawMomentCollision {
    /// Construct from shear viscosity ν and bulk viscosity ζ.
    ///
    /// `ω_v = 1 / (3 ν + 1/2)` and `ω_b = 1 / (9 ζ / 2 + 1/2)`, so that
    /// [`shear_viscosity`](Self::shear_viscosity) returns `ν` and
    /// [`bulk_viscosity`](Self::bulk_viscosity) returns `ζ`. The ghost moments are reset to their
    /// equilibrium values (`ω_g = 1`). `ζ = 0` gives `ω_b = 2`, the edge of linear stability.
    pub fn from_viscosity(nu: f64, zeta: f64) -> Self {
        let tau_v = 3.0 * nu + 0.5;
        let tau_b = 4.5 * zeta + 0.5;
        Self {
            omega_v: 1.0 / tau_v,
            omega_b: 1.0 / tau_b,
            omega_g: 1.0,
        }
    }
    /// Perform raw moment collision.
    ///
    /// Uses the D3Q19 raw moment representation. The hydrodynamic moments (density, momentum)
    /// are conserved, the trace of the stress relaxes at `ω_b`, the trace-free stress at `ω_v`
    /// and the ghost moments at `ω_g`; `rho` and `u` must be the density and velocity of `f`.
    pub fn collide(&self, f: &[f64; 19], rho: f64, u: [f64; 3]) -> [f64; 19] {
        let feq = compute_equilibrium_d3q19(rho, u);
        let m = compute_raw_moments_d3q19(f);
        let meq = compute_raw_moments_d3q19(&feq);
        let mut m_out = [0.0f64; 19];
        // Conserved moments: density and momentum take the equilibrium values.
        m_out[..4].copy_from_slice(&meq[..4]);
        // Second order: the trace relaxes with the bulk rate, the deviators with the shear rate.
        let trace = m[4] + m[5] + m[6];
        let trace_eq = meq[4] + meq[5] + meq[6];
        let trace_out = trace - self.omega_b * (trace - trace_eq);
        for k in 4..7 {
            let deviator = m[k] - trace / 3.0;
            let deviator_eq = meq[k] - trace_eq / 3.0;
            m_out[k] = trace_out / 3.0 + deviator - self.omega_v * (deviator - deviator_eq);
        }
        for k in 7..10 {
            m_out[k] = m[k] - self.omega_v * (m[k] - meq[k]);
        }
        // Ghost moments: third and fourth order.
        for k in 10..19 {
            m_out[k] = m[k] - self.omega_g * (m[k] - meq[k]);
        }
        raw_moments_to_f_d3q19(&m_out)
    }
    /// Return effective shear viscosity.
    pub fn shear_viscosity(&self) -> f64 {
        (1.0 / 3.0) * (1.0 / self.omega_v - 0.5)
    }
    /// Return the effective bulk viscosity `ζ = (2/9) (1/ω_b − 1/2)` of the D3Q19 trace mode.
    pub fn bulk_viscosity(&self) -> f64 {
        (2.0 / 9.0) * (1.0 / self.omega_b - 0.5)
    }
}
