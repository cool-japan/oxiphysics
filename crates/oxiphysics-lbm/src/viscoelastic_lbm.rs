// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Viscoelastic fluid simulation with the Lattice Boltzmann Method.
//!
//! This module provides models for viscoelastic fluids including:
//! - **Maxwell fluid**: solvent + polymer viscosity decomposition
//! - **Oldroyd-B model**: upper-convected Maxwell with viscosity ratio β
//! - **FENE-P model**: finitely-extensible nonlinear elastic (Peterlin closure)
//! - **Viscoelastic LBM solver**: distribution functions + polymer stress coupling
//! - **Extensional flow**: Trouton ratio characterization
//!
//! Key dimensionless numbers:
//! - Weissenberg number: Wi = λ · γ̇  (elastic vs viscous forces)
//! - Deborah number:     De = λ / t_p (relaxation time vs process time)

// ─── Weissenberg / Deborah numbers ────────────────────────────────────────────

/// Compute the Weissenberg number Wi = λ · γ̇.
///
/// Wi compares the polymer relaxation time to the inverse shear rate.
/// Wi > 1 means elastic effects dominate.
///
/// # Arguments
/// * `lambda`     – polymer relaxation time (s)
/// * `shear_rate` – characteristic shear rate (1/s)
pub fn weissenberg_number(lambda: f64, shear_rate: f64) -> f64 {
    lambda * shear_rate
}

/// Compute the Deborah number De = λ / t_process.
///
/// De < 1 → fluid-like; De > 1 → solid-like (elastic) response.
///
/// # Arguments
/// * `lambda`    – polymer relaxation time (s)
/// * `t_process` – characteristic process time (s)
pub fn deborah_number(lambda: f64, t_process: f64) -> f64 {
    if t_process == 0.0 {
        return f64::INFINITY;
    }
    lambda / t_process
}

/// Return the Trouton ratio for a purely Newtonian fluid (= 3).
///
/// For a Newtonian fluid, extensional viscosity = 3 × shear viscosity.
pub fn trouton_ratio_newtonian() -> f64 {
    3.0
}

// ─── MaxwellFluid ─────────────────────────────────────────────────────────────

/// Maxwell fluid with solvent and polymer viscosity contributions.
///
/// The total zero-shear viscosity is η = η_s + η_p.
/// The Weissenberg number requires an external shear rate.
pub struct MaxwellFluid {
    /// Solvent (Newtonian) dynamic viscosity (Pa·s).
    pub viscosity_s: f64,
    /// Polymer (elastic) dynamic viscosity contribution (Pa·s).
    pub viscosity_p: f64,
    /// Polymer stress relaxation time λ (s).
    pub relaxation_time: f64,
}

impl MaxwellFluid {
    /// Create a new Maxwell fluid.
    ///
    /// # Arguments
    /// * `viscosity_s`    – solvent viscosity (Pa·s)
    /// * `viscosity_p`    – polymer viscosity (Pa·s)
    /// * `relaxation_time` – relaxation time λ (s)
    pub fn new(viscosity_s: f64, viscosity_p: f64, relaxation_time: f64) -> Self {
        Self {
            viscosity_s,
            viscosity_p,
            relaxation_time,
        }
    }

    /// Compute the Weissenberg number for a given shear rate.
    ///
    /// Wi = λ · γ̇
    pub fn weissenberg_number(&self, shear_rate: f64) -> f64 {
        self.relaxation_time * shear_rate
    }

    /// Total zero-shear viscosity η = η_s + η_p.
    pub fn total_viscosity(&self) -> f64 {
        self.viscosity_s + self.viscosity_p
    }

    /// Viscosity ratio β = η_s / (η_s + η_p).
    pub fn viscosity_ratio(&self) -> f64 {
        let total = self.total_viscosity();
        if total == 0.0 {
            0.0
        } else {
            self.viscosity_s / total
        }
    }
}

// ─── OldroydBModel ────────────────────────────────────────────────────────────

/// Oldroyd-B constitutive model for dilute polymer solutions.
///
/// Combines an upper-convected Maxwell (UCM) polymer contribution with a
/// Newtonian solvent via the viscosity ratio β = η_s / η_total.
/// The polymer stress tensor τ_p is stored in Voigt-like row-major order
/// (3×3 flattened to 9 components for 3D, with indices 0..8).
pub struct OldroydBModel {
    /// Underlying Maxwell fluid (gives viscosities + relaxation time).
    pub upper_convected_maxwell: MaxwellFluid,
    /// Viscosity ratio β = η_s / η_total ∈ (0, 1).
    pub beta: f64,
    /// Current polymer stress tensor (3×3 row-major, Pa).
    pub stress_tensor: [f64; 9],
}

impl OldroydBModel {
    /// Create an Oldroyd-B model.
    ///
    /// # Arguments
    /// * `ucm`  – underlying Maxwell fluid
    /// * `beta` – solvent viscosity fraction β
    pub fn new(ucm: MaxwellFluid, beta: f64) -> Self {
        Self {
            upper_convected_maxwell: ucm,
            beta,
            stress_tensor: [0.0; 9],
        }
    }

    /// Evolve the polymer stress tensor with a first-order explicit scheme.
    ///
    /// Upper-convected Maxwell evolution (simplified, affine motion):
    ///
    /// ```text
    /// dτ/dt = (η_p/λ)(κ + κᵀ) − τ/λ
    /// ```
    ///
    /// where κ is the velocity gradient tensor.
    ///
    /// # Arguments
    /// * `velocity_gradient` – 3×3 velocity gradient κ (row-major, 1/s)
    /// * `dt`                – time step (s)
    pub fn evolve_stress(&mut self, velocity_gradient: &[f64; 9], dt: f64) {
        let lambda = self.upper_convected_maxwell.relaxation_time;
        let eta_p = self.upper_convected_maxwell.viscosity_p;
        if lambda == 0.0 {
            return;
        }
        let prefactor = eta_p / lambda;
        for i in 0..3 {
            for j in 0..3 {
                let kappa_ij = velocity_gradient[i * 3 + j];
                let kappa_ji = velocity_gradient[j * 3 + i];
                // Symmetric rate-of-strain contribution 2 * D_ij = κ_ij + κ_ji
                let d_ij = kappa_ij + kappa_ji;
                let tau_ij = self.stress_tensor[i * 3 + j];
                let dtau = prefactor * d_ij - tau_ij / lambda;
                self.stress_tensor[i * 3 + j] += dt * dtau;
            }
        }
    }

    /// Return the first normal stress difference N1 = τ_xx − τ_yy.
    pub fn normal_stress_difference_1(&self) -> f64 {
        self.stress_tensor[0] - self.stress_tensor[4]
    }

    /// Return the second normal stress difference N2 = τ_yy − τ_zz.
    pub fn normal_stress_difference_2(&self) -> f64 {
        self.stress_tensor[4] - self.stress_tensor[8]
    }
}

// ─── FenePModel ───────────────────────────────────────────────────────────────

/// FENE-P (Finitely Extensible Nonlinear Elastic – Peterlin) dumbbell model.
///
/// Accounts for finite extensibility of polymer chains.  The Peterlin
/// approximation replaces the exact nonlinear spring with a mean-field closure.
///
/// Spring force: **F** = H · f(Q²) · **Q**
/// with Peterlin factor f(Q²) = 1 / (1 − Q²/L²).
pub struct FenePModel {
    /// Maximum extensibility parameter L² (dimensionless).
    pub extensibility_l: f64,
    /// Spring constant H (N/m).
    pub spring_constant: f64,
}

impl FenePModel {
    /// Create a FENE-P model.
    ///
    /// # Arguments
    /// * `extensibility_l`  – maximum extensibility L (chain fully stretched)
    /// * `spring_constant`  – spring stiffness H (N/m)
    pub fn new(extensibility_l: f64, spring_constant: f64) -> Self {
        Self {
            extensibility_l,
            spring_constant,
        }
    }

    /// Compute the FENE spring force vector for a given end-to-end vector Q.
    ///
    /// **F** = H · f(|Q|²) · **Q**,   f(r) = 1/(1 − r/L²)
    ///
    /// Returns zero vector if Q is fully extended (|Q|² ≥ L²).
    pub fn compute_spring_force(&self, q: [f64; 3]) -> [f64; 3] {
        let l2 = self.extensibility_l * self.extensibility_l;
        let q2 = q[0] * q[0] + q[1] * q[1] + q[2] * q[2];
        if q2 >= l2 {
            return [0.0; 3];
        }
        let f_peterlin = 1.0 / (1.0 - q2 / l2);
        let h = self.spring_constant;
        [
            h * f_peterlin * q[0],
            h * f_peterlin * q[1],
            h * f_peterlin * q[2],
        ]
    }

    /// Compute the polymer contribution to the stress tensor from
    /// a single dumbbell configuration Q at number concentration `conc`.
    ///
    /// Kramers-Kirkwood expression (outer product, Peterlin-weighted):
    ///   τ_p = conc · H · f(Q²) · Q ⊗ Q
    ///
    /// Returns a 3×3 row-major stress tensor (Pa).
    pub fn stress_from_distribution(&self, q: [f64; 3], conc: f64) -> [f64; 9] {
        let l2 = self.extensibility_l * self.extensibility_l;
        let q2 = q[0] * q[0] + q[1] * q[1] + q[2] * q[2];
        if q2 >= l2 {
            return [0.0; 9];
        }
        let f_peterlin = 1.0 / (1.0 - q2 / l2);
        let prefactor = conc * self.spring_constant * f_peterlin;
        let mut tau = [0.0f64; 9];
        for i in 0..3 {
            for j in 0..3 {
                tau[i * 3 + j] = prefactor * q[i] * q[j];
            }
        }
        tau
    }
}

// ─── ViscoelasticLbm ──────────────────────────────────────────────────────────

/// Viscoelastic Lattice Boltzmann solver (D2Q9 layout, periodic box).
///
/// Couples a BGK lattice Boltzmann fluid, which carries the Newtonian solvent,
/// with an **Oldroyd-B** polymer stress field.  Everything is in lattice units
/// (grid spacing 1, one time step = 1, density ≈ 1):
///
/// * The solvent kinematic viscosity is `ν_s = (1/ω − 1/2)/3`.
/// * The polymer stress `τ` obeys the upper-convected Maxwell equation
///   `∂τ/∂t + u·∇τ = κτ + τκᵀ + (η_p/λ)(κ + κᵀ) − τ/λ`, `κ_ij = ∂u_i/∂x_j`.
///   `κ` is the central-difference velocity gradient, `u·∇τ` is advected with a
///   MC-limited second-order upwind scheme (sub-cycled to CFL ≤ 0.4), and
///   the local right-hand side is integrated with the second-order
///   exponential integrator ETDRK2 ([`oldroyd_b_stress_update`]), so the
///   relaxation term is unconditionally stable and `λ < 1` is allowed.
///   `λ ≤ 0` is the Newtonian limit `τ = η_p (κ + κᵀ)`.
/// * The stress acts on the fluid as the force density `F = ∇·τ` (central
///   differences) through the Guo forcing scheme, with the half-force velocity
///   `u = (Σ c_i f_i + F/2)/ρ`.  The scheme adds exactly `F` of momentum per
///   step and no mass, and the discrete divergence sums to zero over the
///   periodic box, so total mass and momentum are conserved.
///
/// In simple shear the stress relaxes to `τ_xy = η_p γ̇`,
/// `N₁ = τ_xx − τ_yy = 2 η_p λ γ̇²`, and the total shear viscosity is
/// `η_s + η_p`.
///
/// # Stability
///
/// The polymer force enters the collision explicitly.  When the stress follows
/// the strain rate almost instantly (`λ ≲ 1`) the coupling therefore behaves
/// like an explicit viscous term and is only conditionally stable: a
/// small-amplitude vortex on a 32 × 32 box stayed bounded for `η_p ≲ 0.4`
/// with `0.8 ≤ ω ≤ 1.8`.  For `λ ≳ 5` much larger `η_p` are admissible (at
/// least 3.4 at `λ = 5`, and at least 8 at `λ = 50`).  These limits are
/// empirical: they were found by bisection on that test case, not derived.
///
/// Distribution functions are stored as one `Vec<f64>` of 9 populations per
/// lattice site (node index `iy * nx + ix`).
pub struct ViscoelasticLbm {
    /// Number of lattice nodes in x-direction.
    pub nx: usize,
    /// Number of lattice nodes in y-direction.
    pub ny: usize,
    /// Distribution functions f_i at each node.  Layout: `[node][velocity]`.
    pub f_distributions: Vec<Vec<f64>>,
    /// Polymer stress tensor at each node (3×3 row-major).
    pub stress_field: Vec<[f64; 9]>,
    /// BGK relaxation parameter ω = 1/τ.
    pub omega: f64,
    /// Polymer relaxation time λ (lattice units).
    pub lambda: f64,
    /// Polymer viscosity contribution η_p (lattice units).
    pub eta_p: f64,
}

/// D2Q9 equilibrium weights.
const W9: [f64; 9] = [
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

/// D2Q9 x-velocity components.
const CX9: [f64; 9] = [0.0, 1.0, 0.0, -1.0, 0.0, 1.0, -1.0, -1.0, 1.0];

/// D2Q9 y-velocity components.
const CY9: [f64; 9] = [0.0, 0.0, 1.0, 0.0, -1.0, 1.0, 1.0, -1.0, -1.0];

/// Largest sum of the Courant numbers `(|u_x| + |u_y|) Δt/Δx` of one
/// stress-advection sub-step (the MC-limited reconstruction with forward Euler
/// is total-variation diminishing up to 1/2 and unstable beyond it).
const ADVECTION_CFL: f64 = 0.4;

/// Upper bound on the number of stress-advection sub-steps per time step.
///
/// A valid lattice Boltzmann state has `|u| ≪ 1`, i.e. a Courant number far
/// below one, so a healthy simulation never gets near this bound.  For a
/// diverging one (velocities of order 10² and beyond) it keeps the cost of a
/// step finite instead of sub-cycling once per [`ADVECTION_CFL`] of the
/// blown-up Courant number.
const MAX_ADVECTION_SUBSTEPS: usize = 64;

/// Number of sub-steps that keep the Courant number of every sub-step at or
/// below [`ADVECTION_CFL`] when advecting over `dt` with the largest Courant
/// rate `max_courant` (at least one, at most [`MAX_ADVECTION_SUBSTEPS`]).
fn advection_substeps(max_courant: f64, dt: f64) -> usize {
    (((max_courant * dt) / ADVECTION_CFL).ceil().max(1.0) as usize).min(MAX_ADVECTION_SUBSTEPS)
}

/// `φ₁(z) = (eᶻ − 1)/z`, accurate for all `z` including `z → 0`.
fn phi1(z: f64) -> f64 {
    if z == 0.0 { 1.0 } else { z.exp_m1() / z }
}

/// `φ₂(z) = (eᶻ − 1 − z)/z²`; the Taylor series replaces the direct formula
/// where the latter loses digits to cancellation.
fn phi2(z: f64) -> f64 {
    if z.abs() < 0.05 {
        let mut term = 0.5;
        let mut sum = term;
        for k in 1..10 {
            term *= z / (k + 2) as f64;
            sum += term;
        }
        sum
    } else {
        (z.exp_m1() - z) / (z * z)
    }
}

/// The part of the Oldroyd-B right-hand side that does not relax:
/// `κτ + τκᵀ + (η_p/λ)(κ + κᵀ)` with `prefactor = η_p/λ`.
///
/// `(κτ)_ij = Σ_k κ_ik τ_kj` and `(τκᵀ)_ij = Σ_k τ_ik κ_jk`, with
/// `κ_ij = ∂u_i/∂x_j`.
fn upper_convected_source(tau: &[f64; 9], kappa: &[f64; 9], prefactor: f64) -> [f64; 9] {
    std::array::from_fn(|index| {
        let (i, j) = (index / 3, index % 3);
        let mut rate = prefactor * (kappa[i * 3 + j] + kappa[j * 3 + i]);
        for k in 0..3 {
            rate += kappa[i * 3 + k] * tau[k * 3 + j] + tau[i * 3 + k] * kappa[j * 3 + k];
        }
        rate
    })
}

/// Time-step dependent weights of the Oldroyd-B stress update, computed once
/// per step so that the per-node work contains no transcendental functions.
#[derive(Debug, Clone, Copy)]
enum StressIntegrator {
    /// `λ ≤ 0`: instantaneous relaxation, `τ = η_p (κ + κᵀ)`.
    Newtonian {
        /// Polymer viscosity η_p.
        eta_p: f64,
    },
    /// ETDRK2 for the Maxwell relaxation `−τ/λ`.
    Maxwell {
        /// η_p / λ.
        prefactor: f64,
        /// `exp(−h/λ)`.
        decay: f64,
        /// `h φ₁(−h/λ)`.
        w1: f64,
        /// `h φ₂(−h/λ)`.
        w2: f64,
    },
}

impl StressIntegrator {
    fn new(lambda: f64, eta_p: f64, h: f64) -> Self {
        if lambda <= 0.0 {
            return Self::Newtonian { eta_p };
        }
        let z = -h / lambda;
        Self::Maxwell {
            prefactor: eta_p / lambda,
            decay: z.exp(),
            w1: h * phi1(z),
            w2: h * phi2(z),
        }
    }

    fn advance(&self, tau: &[f64; 9], kappa: &[f64; 9]) -> [f64; 9] {
        match *self {
            // With τ = 0 only the source term η_p (κ + κᵀ) is left.
            Self::Newtonian { eta_p } => upper_convected_source(&[0.0; 9], kappa, eta_p),
            Self::Maxwell {
                prefactor,
                decay,
                w1,
                w2,
            } => {
                let n0 = upper_convected_source(tau, kappa, prefactor);
                let predictor: [f64; 9] = std::array::from_fn(|k| decay * tau[k] + w1 * n0[k]);
                let n1 = upper_convected_source(&predictor, kappa, prefactor);
                std::array::from_fn(|k| predictor[k] + w2 * (n1[k] - n0[k]))
            }
        }
    }
}

/// Advance an Oldroyd-B polymer stress over a time `h` in a velocity gradient
/// that is constant during the step.
///
/// Integrates `dτ/dt = κτ + τκᵀ + (η_p/λ)(κ + κᵀ) − τ/λ`
/// (`κ_ij = ∂u_i/∂x_j`, row-major 3×3) with the second-order exponential
/// Runge–Kutta scheme ETDRK2: the stiff relaxation `−τ/λ` is integrated
/// exactly, the stretching and source terms to second order.  The scheme is
/// stable for any `h/λ`, and a stress that already satisfies the steady
/// equation is reproduced exactly.
///
/// For `λ ≤ 0` the relaxation is instantaneous and the result is the
/// Newtonian limit `τ = η_p (κ + κᵀ)`.
///
/// In steady simple shear `κ_xy = γ̇` this converges to `τ_xy = η_p γ̇`,
/// `τ_xx = 2 η_p λ γ̇²`, `τ_yy = 0`.
pub fn oldroyd_b_stress_update(
    tau: &[f64; 9],
    kappa: &[f64; 9],
    lambda: f64,
    eta_p: f64,
    h: f64,
) -> [f64; 9] {
    StressIntegrator::new(lambda, eta_p, h).advance(tau, kappa)
}

/// Periodic wrap of `value` into `0..n` (`n > 0`).
#[inline]
fn wrap_index(value: isize, n: usize) -> usize {
    value.rem_euclid(n as isize) as usize
}

/// Monotonized-central slope limiter: the central slope, clipped so that the
/// reconstruction stays within twice the smaller one-sided difference (zero at
/// local extrema).
#[inline]
fn mc_slope(backward: f64, forward: f64) -> f64 {
    if backward * forward <= 0.0 {
        return 0.0;
    }
    let central = 0.5 * (backward + forward);
    let bound = 2.0 * backward.abs().min(forward.abs());
    central.signum() * central.abs().min(bound)
}

/// Second-order upwind difference `τ_{i+1/2} − τ_{i−1/2}` with MC-limited
/// slopes, from the stencil `[τ_{i−2}, τ_{i−1}, τ_i, τ_{i+1}, τ_{i+2}]` and the
/// local advection velocity (the upwind side follows its sign).
#[inline]
fn limited_upwind_difference(v: [f64; 5], velocity: f64) -> f64 {
    if velocity >= 0.0 {
        let right = v[2] + 0.5 * mc_slope(v[2] - v[1], v[3] - v[2]);
        let left = v[1] + 0.5 * mc_slope(v[1] - v[0], v[2] - v[1]);
        right - left
    } else {
        let right = v[3] - 0.5 * mc_slope(v[4] - v[3], v[3] - v[2]);
        let left = v[2] - 0.5 * mc_slope(v[3] - v[2], v[2] - v[1]);
        right - left
    }
}

/// Advect the stress field with the velocity field over a time `dt`
/// (`∂τ/∂t + u·∇τ = 0`, periodic), sub-cycling so that every sub-step has a
/// Courant number of at most [`ADVECTION_CFL`] (and the number of sub-steps is
/// bounded by [`MAX_ADVECTION_SUBSTEPS`]).  Components that vanish everywhere
/// are skipped.
fn advect_stress(stress: &mut [[f64; 9]], ux: &[f64], uy: &[f64], nx: usize, ny: usize, dt: f64) {
    // The two directions are updated together, so the Courant numbers add.
    let max_courant = ux
        .iter()
        .zip(uy.iter())
        .fold(0.0_f64, |m, (u, v)| m.max(u.abs() + v.abs()));
    if max_courant == 0.0 || nx == 0 || ny == 0 {
        return;
    }
    let substeps = advection_substeps(max_courant, dt);
    let h = dt / substeps as f64;
    let mut active = [false; 9];
    for node_stress in stress.iter() {
        for (flag, value) in active.iter_mut().zip(node_stress.iter()) {
            *flag |= *value != 0.0;
        }
    }
    let mut next = stress.to_vec();
    for _ in 0..substeps {
        for y in 0..ny {
            for x in 0..nx {
                let node = y * nx + x;
                let along_x: [usize; 5] =
                    std::array::from_fn(|k| y * nx + wrap_index(x as isize + k as isize - 2, nx));
                let along_y: [usize; 5] =
                    std::array::from_fn(|k| wrap_index(y as isize + k as isize - 2, ny) * nx + x);
                for c in (0..9).filter(|&c| active[c]) {
                    let dx = limited_upwind_difference(
                        std::array::from_fn(|k| stress[along_x[k]][c]),
                        ux[node],
                    );
                    let dy = limited_upwind_difference(
                        std::array::from_fn(|k| stress[along_y[k]][c]),
                        uy[node],
                    );
                    next[node][c] = stress[node][c] - h * (ux[node] * dx + uy[node] * dy);
                }
            }
        }
        stress.copy_from_slice(&next);
    }
}

/// Central-difference velocity gradient `κ_ij = ∂u_i/∂x_j` (periodic, unit
/// spacing) as row-major 3×3 tensors; the out-of-plane entries are zero.
fn velocity_gradient(ux: &[f64], uy: &[f64], nx: usize, ny: usize) -> Vec<[f64; 9]> {
    let mut kappa = vec![[0.0_f64; 9]; nx * ny];
    for y in 0..ny {
        for x in 0..nx {
            let node = y * nx + x;
            let east = y * nx + wrap_index(x as isize + 1, nx);
            let west = y * nx + wrap_index(x as isize - 1, nx);
            let north = wrap_index(y as isize + 1, ny) * nx + x;
            let south = wrap_index(y as isize - 1, ny) * nx + x;
            kappa[node][0] = 0.5 * (ux[east] - ux[west]);
            kappa[node][1] = 0.5 * (ux[north] - ux[south]);
            kappa[node][3] = 0.5 * (uy[east] - uy[west]);
            kappa[node][4] = 0.5 * (uy[north] - uy[south]);
        }
    }
    kappa
}

/// Velocity `u = (m + F/2)/ρ` of every node from its density `ρ`, momentum
/// `m = Σ c_i f_i` and force density `F` (zero where `ρ ≤ 0`), returned as
/// separate x and y components.
fn half_force_velocity(
    rho: &[f64],
    momentum: &[[f64; 2]],
    force: &[[f64; 2]],
) -> (Vec<f64>, Vec<f64>) {
    let mut ux = vec![0.0_f64; rho.len()];
    let mut uy = vec![0.0_f64; rho.len()];
    for (node, &density) in rho.iter().enumerate() {
        if density > 0.0 {
            ux[node] = (momentum[node][0] + 0.5 * force[node][0]) / density;
            uy[node] = (momentum[node][1] + 0.5 * force[node][1]) / density;
        }
    }
    (ux, uy)
}

impl ViscoelasticLbm {
    /// Create a new viscoelastic LBM domain initialised to rest.
    ///
    /// # Arguments
    /// * `nx`, `ny` – domain dimensions
    /// * `omega`    – BGK relaxation frequency (1/τ); the solvent viscosity is
    ///   `(1/ω − 1/2)/3`
    /// * `lambda`   – polymer relaxation time (lattice units); `λ ≤ 0` is the
    ///   Newtonian limit
    /// * `eta_p`    – polymer viscosity (lattice units)
    pub fn new(nx: usize, ny: usize, omega: f64, lambda: f64, eta_p: f64) -> Self {
        let n = nx * ny;
        // Equilibrium at rest: f_i = w_i * rho_0 with rho_0 = 1
        let f0: Vec<f64> = W9.to_vec();
        Self {
            nx,
            ny,
            f_distributions: vec![f0; n],
            stress_field: vec![[0.0; 9]; n],
            omega,
            lambda,
            eta_p,
        }
    }

    /// Polymer force density `F = ∇·τ` at node `(ix, iy)`.
    ///
    /// `F_i = ∂_x τ_xi + ∂_y τ_yi` with periodic central differences.  It sums
    /// to zero over the box, so it never changes the total momentum.
    pub fn stress_force(&self, ix: usize, iy: usize) -> [f64; 2] {
        let (nx, ny) = (self.nx, self.ny);
        if nx == 0 || ny == 0 {
            return [0.0; 2];
        }
        let east = &self.stress_field[iy * nx + wrap_index(ix as isize + 1, nx)];
        let west = &self.stress_field[iy * nx + wrap_index(ix as isize - 1, nx)];
        let north = &self.stress_field[wrap_index(iy as isize + 1, ny) * nx + ix];
        let south = &self.stress_field[wrap_index(iy as isize - 1, ny) * nx + ix];
        [
            0.5 * ((east[0] - west[0]) + (north[3] - south[3])),
            0.5 * ((east[1] - west[1]) + (north[4] - south[4])),
        ]
    }

    /// Force density `F = ∇·τ` of the current stress field at every node
    /// (node order `iy * nx + ix`).
    fn stress_force_field(&self) -> Vec<[f64; 2]> {
        (0..self.nx * self.ny)
            .map(|node| self.stress_force(node % self.nx, node / self.nx))
            .collect()
    }

    /// Perform a single LBM time step (one lattice time unit).
    ///
    /// 1. Evaluate the macroscopic density and the half-force-corrected
    ///    velocity `u = (Σ c_i f_i + F/2)/ρ` for the force `F = ∇·τ` of the
    ///    current stress.
    /// 2. Differentiate the velocity (`κ_ij = ∂u_i/∂x_j`) and advance the
    ///    polymer stress: advection by `u`, then the Oldroyd-B stretching,
    ///    source and relaxation terms ([`oldroyd_b_stress_update`]).
    /// 3. Take the force `F = ∇·τ` of the *updated* stress, and the velocity
    ///    that belongs to it.  Using the new stress instead of the old one
    ///    removes a one-step delay from the feedback loop between velocity and
    ///    stress; the range of `η_p` for which the coupling stays stable grows
    ///    by a factor of 1.2 to 6, depending on `ω` (see the stability notes of
    ///    [`ViscoelasticLbm`]).
    /// 4. BGK collision with Guo forcing `F`, with every post-collision
    ///    population streamed straight to its neighbour (periodic boundaries)
    ///    in the same sweep.
    pub fn step(&mut self) {
        let nx = self.nx;
        let ny = self.ny;
        let n = nx * ny;
        if n == 0 {
            return;
        }

        // ── 1. Macroscopic fields with the force of the current stress ────────
        let mut rho = vec![0.0_f64; n];
        let mut momentum = vec![[0.0_f64; 2]; n];
        for node in 0..n {
            let f = &self.f_distributions[node];
            rho[node] = f.iter().sum();
            momentum[node] = [
                f.iter().zip(CX9.iter()).map(|(fi, ci)| fi * ci).sum(),
                f.iter().zip(CY9.iter()).map(|(fi, ci)| fi * ci).sum(),
            ];
        }
        let (mut ux, mut uy) = half_force_velocity(&rho, &momentum, &self.stress_force_field());

        // ── 2. Polymer stress ─────────────────────────────────────────────────
        let kappa = velocity_gradient(&ux, &uy, nx, ny);
        advect_stress(&mut self.stress_field, &ux, &uy, nx, ny, 1.0);
        let integrator = StressIntegrator::new(self.lambda, self.eta_p, 1.0);
        for (tau, kappa_node) in self.stress_field.iter_mut().zip(kappa.iter()) {
            *tau = integrator.advance(tau, kappa_node);
        }

        // ── 3. Force of the updated stress ────────────────────────────────────
        let force = self.stress_force_field();
        (ux, uy) = half_force_velocity(&rho, &momentum, &force);

        // ── 4. Collision with Guo forcing, streamed to the neighbours ─────────
        // Each post-collision population is written straight to its periodic
        // destination, so a step needs one scratch lattice instead of two.
        let guo_factor = 1.0 - 0.5 * self.omega;
        let mut f_stream = vec![[0.0_f64; 9]; n];
        for y in 0..ny {
            for x in 0..nx {
                let node = y * nx + x;
                let f = &self.f_distributions[node];
                let live = rho[node] > 0.0;
                let (u_x, u_y) = (ux[node], uy[node]);
                let [fx, fy] = force[node];
                let u2 = u_x * u_x + u_y * u_y;
                for i in 0..9 {
                    let post = if live {
                        let cu = CX9[i] * u_x + CY9[i] * u_y;
                        let feq = W9[i] * rho[node] * (1.0 + 3.0 * cu + 4.5 * cu * cu - 1.5 * u2);
                        // Guo: w_i [ (c_i − u)/c_s² + (c_i·u) c_i / c_s⁴ ] · F
                        let source = W9[i]
                            * (3.0 * ((CX9[i] - u_x) * fx + (CY9[i] - u_y) * fy)
                                + 9.0 * cu * (CX9[i] * fx + CY9[i] * fy));
                        f[i] - self.omega * (f[i] - feq) + guo_factor * source
                    } else {
                        f[i]
                    };
                    let xd = wrap_index(x as isize + CX9[i] as isize, nx);
                    let yd = wrap_index(y as isize + CY9[i] as isize, ny);
                    f_stream[yd * nx + xd][i] = post;
                }
            }
        }
        // Reuse the per-node buffers of the lattice instead of replacing them.
        for (dst, post) in self.f_distributions.iter_mut().zip(f_stream.iter()) {
            dst.clear();
            dst.extend_from_slice(post);
        }
    }

    /// Return macroscopic density at node (ix, iy).
    pub fn density(&self, ix: usize, iy: usize) -> f64 {
        self.f_distributions[iy * self.nx + ix].iter().sum()
    }

    /// Return macroscopic velocity (ux, uy) at node (ix, iy).
    ///
    /// This is the physical (half-force-corrected) velocity
    /// `u = (Σ c_i f_i + F/2)/ρ` with the polymer force `F = ∇·τ` of the
    /// current stress field, i.e. the velocity the next [`ViscoelasticLbm::step`]
    /// differentiates.
    pub fn velocity(&self, ix: usize, iy: usize) -> (f64, f64) {
        let f = &self.f_distributions[iy * self.nx + ix];
        let rho: f64 = f.iter().sum();
        if rho == 0.0 {
            return (0.0, 0.0);
        }
        let force = self.stress_force(ix, iy);
        let mx: f64 = f.iter().zip(CX9.iter()).map(|(fi, ci)| fi * ci).sum();
        let my: f64 = f.iter().zip(CY9.iter()).map(|(fi, ci)| fi * ci).sum();
        ((mx + 0.5 * force[0]) / rho, (my + 0.5 * force[1]) / rho)
    }
}

// ─── ExtensionalFlow ─────────────────────────────────────────────────────────

/// Extensional flow characterization for viscoelastic fluids.
///
/// Extensional flows (e.g. uniaxial, planar elongation) are dominated by
/// elastic effects.  The Trouton ratio η_E / η compares extensional to shear
/// viscosity.
pub struct ExtensionalFlow {
    /// Extension rate ε̇ (1/s).
    pub extension_rate: f64,
}

impl ExtensionalFlow {
    /// Create an extensional flow with the given extension rate ε̇.
    pub fn new(extension_rate: f64) -> Self {
        Self { extension_rate }
    }

    /// Compute the Trouton ratio Tr = η_E / η.
    ///
    /// # Arguments
    /// * `viscosity`            – shear viscosity η (Pa·s)
    /// * `extensional_viscosity` – extensional (Trouton) viscosity η_E (Pa·s)
    pub fn compute_trouton_ratio(&self, viscosity: f64, extensional_viscosity: f64) -> f64 {
        if viscosity == 0.0 {
            return 0.0;
        }
        extensional_viscosity / viscosity
    }

    /// Predict the Hencky strain ε = ε̇ · t accumulated after time t.
    pub fn hencky_strain(&self, t: f64) -> f64 {
        self.extension_rate * t
    }
}

// ─── Tests ────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // ── Dimensionless numbers ─────────────────────────────────────────────────

    #[test]
    fn test_weissenberg_basic() {
        let wi = weissenberg_number(1.0, 2.0);
        assert!((wi - 2.0).abs() < 1e-12);
    }

    #[test]
    fn test_weissenberg_zero_shear() {
        assert_eq!(weissenberg_number(1.0, 0.0), 0.0);
    }

    #[test]
    fn test_weissenberg_zero_lambda() {
        assert_eq!(weissenberg_number(0.0, 5.0), 0.0);
    }

    #[test]
    fn test_deborah_basic() {
        let de = deborah_number(0.1, 1.0);
        assert!((de - 0.1).abs() < 1e-12);
    }

    #[test]
    fn test_deborah_fast_process() {
        // De < 1: fluid-like
        let de = deborah_number(0.01, 1.0);
        assert!(de < 1.0);
    }

    #[test]
    fn test_deborah_slow_process() {
        // De > 1: solid-like
        let de = deborah_number(10.0, 1.0);
        assert!(de > 1.0);
    }

    #[test]
    fn test_deborah_zero_process_time() {
        let de = deborah_number(1.0, 0.0);
        assert!(de.is_infinite());
    }

    #[test]
    fn test_trouton_newtonian() {
        assert!((trouton_ratio_newtonian() - 3.0).abs() < 1e-12);
    }

    // ── MaxwellFluid ──────────────────────────────────────────────────────────

    #[test]
    fn test_maxwell_total_viscosity() {
        let f = MaxwellFluid::new(0.001, 0.009, 0.1);
        assert!((f.total_viscosity() - 0.01).abs() < 1e-12);
    }

    #[test]
    fn test_maxwell_viscosity_ratio() {
        let f = MaxwellFluid::new(0.001, 0.009, 0.1);
        assert!((f.viscosity_ratio() - 0.1).abs() < 1e-12);
    }

    #[test]
    fn test_maxwell_weissenberg() {
        let f = MaxwellFluid::new(0.001, 0.009, 2.0);
        assert!((f.weissenberg_number(5.0) - 10.0).abs() < 1e-12);
    }

    #[test]
    fn test_maxwell_zero_total_viscosity() {
        let f = MaxwellFluid::new(0.0, 0.0, 1.0);
        assert_eq!(f.viscosity_ratio(), 0.0);
    }

    #[test]
    fn test_maxwell_pure_solvent() {
        let f = MaxwellFluid::new(0.01, 0.0, 0.5);
        assert!((f.viscosity_ratio() - 1.0).abs() < 1e-12);
    }

    // ── OldroydBModel ─────────────────────────────────────────────────────────

    #[test]
    fn test_oldroyd_b_initial_stress_zero() {
        let ucm = MaxwellFluid::new(0.001, 0.009, 0.1);
        let model = OldroydBModel::new(ucm, 0.1);
        assert_eq!(model.stress_tensor, [0.0; 9]);
    }

    #[test]
    fn test_oldroyd_b_stress_evolution_positive() {
        let ucm = MaxwellFluid::new(0.001, 0.009, 1.0);
        let mut model = OldroydBModel::new(ucm, 0.1);
        let kappa = [0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0f64];
        model.evolve_stress(&kappa, 0.01);
        // τ_xy should become non-zero after shear
        assert!(model.stress_tensor[1].abs() > 0.0);
    }

    #[test]
    fn test_oldroyd_b_stress_relaxes() {
        let ucm = MaxwellFluid::new(0.001, 0.009, 0.1);
        let mut model = OldroydBModel::new(ucm, 0.1);
        // Inject stress manually
        model.stress_tensor[0] = 1.0;
        let kappa = [0.0f64; 9];
        model.evolve_stress(&kappa, 0.01);
        // Stress should decay
        assert!(model.stress_tensor[0] < 1.0);
    }

    #[test]
    fn test_oldroyd_b_n1_initial() {
        let ucm = MaxwellFluid::new(0.001, 0.009, 0.1);
        let model = OldroydBModel::new(ucm, 0.1);
        assert_eq!(model.normal_stress_difference_1(), 0.0);
    }

    #[test]
    fn test_oldroyd_b_n2_initial() {
        let ucm = MaxwellFluid::new(0.001, 0.009, 0.1);
        let model = OldroydBModel::new(ucm, 0.1);
        assert_eq!(model.normal_stress_difference_2(), 0.0);
    }

    #[test]
    fn test_oldroyd_b_zero_lambda_no_change() {
        let ucm = MaxwellFluid::new(0.001, 0.009, 0.0);
        let mut model = OldroydBModel::new(ucm, 0.1);
        let kappa = [0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0f64];
        model.evolve_stress(&kappa, 0.01);
        assert_eq!(model.stress_tensor, [0.0; 9]);
    }

    // ── FenePModel ────────────────────────────────────────────────────────────

    #[test]
    fn test_fene_zero_extension() {
        let m = FenePModel::new(10.0, 1.0);
        let f = m.compute_spring_force([0.0, 0.0, 0.0]);
        assert_eq!(f, [0.0, 0.0, 0.0]);
    }

    #[test]
    fn test_fene_small_extension() {
        let m = FenePModel::new(10.0, 1.0);
        // Q = [1, 0, 0], Q² = 1, L² = 100
        let f = m.compute_spring_force([1.0, 0.0, 0.0]);
        // f = H * 1/(1 - 1/100) * Q = 1 * (100/99) ≈ 1.0101
        assert!((f[0] - 100.0 / 99.0).abs() < 1e-6);
        assert_eq!(f[1], 0.0);
        assert_eq!(f[2], 0.0);
    }

    #[test]
    fn test_fene_at_full_extension() {
        let m = FenePModel::new(5.0, 1.0);
        // Q² = L² → fully extended → zero force
        let f = m.compute_spring_force([5.0, 0.0, 0.0]);
        assert_eq!(f, [0.0, 0.0, 0.0]);
    }

    #[test]
    fn test_fene_beyond_extension() {
        let m = FenePModel::new(3.0, 1.0);
        let f = m.compute_spring_force([4.0, 0.0, 0.0]);
        assert_eq!(f, [0.0, 0.0, 0.0]);
    }

    #[test]
    fn test_fene_stress_from_distribution_zero() {
        let m = FenePModel::new(10.0, 1.0);
        let tau = m.stress_from_distribution([0.0, 0.0, 0.0], 1.0);
        assert_eq!(tau, [0.0; 9]);
    }

    #[test]
    fn test_fene_stress_from_distribution_symmetry() {
        let m = FenePModel::new(10.0, 1.0);
        let tau = m.stress_from_distribution([1.0, 0.5, 0.0], 1.0);
        // Should be symmetric: tau[1] == tau[3] (xy == yx)
        assert!((tau[1] - tau[3]).abs() < 1e-12);
    }

    #[test]
    fn test_fene_stress_off_diagonal_zero_if_component_zero() {
        let m = FenePModel::new(10.0, 1.0);
        // Q along x only
        let tau = m.stress_from_distribution([1.0, 0.0, 0.0], 1.0);
        assert!((tau[1]).abs() < 1e-12); // xy component
        assert!((tau[2]).abs() < 1e-12); // xz component
    }

    // ── ViscoelasticLbm ───────────────────────────────────────────────────────

    #[test]
    fn test_lbm_initialisation_density() {
        let lbm = ViscoelasticLbm::new(8, 8, 1.0, 0.1, 0.01);
        // Density at rest = sum of D2Q9 weights = 1
        let rho = lbm.density(0, 0);
        assert!((rho - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_lbm_initialisation_velocity_zero() {
        let lbm = ViscoelasticLbm::new(4, 4, 1.0, 0.1, 0.01);
        let (ux, uy) = lbm.velocity(2, 2);
        assert!(ux.abs() < 1e-12);
        assert!(uy.abs() < 1e-12);
    }

    #[test]
    fn test_lbm_step_does_not_crash() {
        let mut lbm = ViscoelasticLbm::new(8, 8, 1.0, 0.1, 0.01);
        lbm.step();
    }

    #[test]
    fn test_lbm_mass_conservation() {
        let mut lbm = ViscoelasticLbm::new(4, 4, 1.0, 0.0, 0.0);
        let total_before: f64 = (0..4 * 4)
            .map(|n| lbm.f_distributions[n].iter().sum::<f64>())
            .sum();
        lbm.step();
        let total_after: f64 = (0..4 * 4)
            .map(|n| lbm.f_distributions[n].iter().sum::<f64>())
            .sum();
        assert!((total_before - total_after).abs() < 1e-8);
    }

    #[test]
    fn test_lbm_stress_field_size() {
        let lbm = ViscoelasticLbm::new(6, 4, 1.0, 0.1, 0.01);
        assert_eq!(lbm.stress_field.len(), 24);
    }

    // ── ExtensionalFlow ───────────────────────────────────────────────────────

    #[test]
    fn test_extensional_trouton_newtonian_value() {
        let ef = ExtensionalFlow::new(1.0);
        let tr = ef.compute_trouton_ratio(1.0, 3.0);
        assert!((tr - 3.0).abs() < 1e-12);
    }

    #[test]
    fn test_extensional_trouton_zero_viscosity() {
        let ef = ExtensionalFlow::new(1.0);
        let tr = ef.compute_trouton_ratio(0.0, 3.0);
        assert_eq!(tr, 0.0);
    }

    #[test]
    fn test_extensional_trouton_viscoelastic() {
        // Viscoelastic fluid has Tr >> 3
        let ef = ExtensionalFlow::new(10.0);
        let tr = ef.compute_trouton_ratio(0.1, 5.0);
        assert!(tr > 3.0);
    }

    #[test]
    fn test_extensional_hencky_strain() {
        let ef = ExtensionalFlow::new(2.0);
        assert!((ef.hencky_strain(3.0) - 6.0).abs() < 1e-12);
    }

    #[test]
    fn test_extensional_hencky_zero_time() {
        let ef = ExtensionalFlow::new(5.0);
        assert_eq!(ef.hencky_strain(0.0), 0.0);
    }

    // ── Combined / integration ────────────────────────────────────────────────

    #[test]
    fn test_wi_de_ratio() {
        // Wi/De = t_process * gamma_dot
        let lambda = 0.5;
        let gamma = 2.0;
        let t_proc = 1.0;
        let wi = weissenberg_number(lambda, gamma);
        let de = deborah_number(lambda, t_proc);
        // wi/de = gamma * t_process
        assert!((wi / de - gamma * t_proc).abs() < 1e-10);
    }

    #[test]
    fn test_maxwell_and_oldroyd_b_consistency() {
        let ucm = MaxwellFluid::new(0.001, 0.009, 1.0);
        let model = OldroydBModel::new(ucm, 0.1);
        // beta should match
        assert!((model.beta - 0.1).abs() < 1e-12);
    }

    #[test]
    fn test_fene_force_direction() {
        // Force should be along Q
        let m = FenePModel::new(10.0, 1.0);
        let q = [2.0, 1.0, 0.0];
        let f = m.compute_spring_force(q);
        // f[0]/q[0] == f[1]/q[1]
        assert!((f[0] / q[0] - f[1] / q[1]).abs() < 1e-10);
    }

    #[test]
    fn test_lbm_polymer_stress_evolves() {
        let mut lbm = ViscoelasticLbm::new(4, 4, 1.0, 1.0, 0.01);
        // Set a non-zero distribution to create non-zero velocity
        lbm.f_distributions[5][1] = 0.5;
        lbm.step();
        // At least one stress component should be non-zero
        let any_nonzero = lbm.stress_field.iter().any(|s| s.iter().any(|&v| v != 0.0));
        assert!(any_nonzero);
    }

    #[test]
    fn test_maxwell_relaxation_time() {
        let f = MaxwellFluid::new(0.001, 0.009, 0.5);
        assert!((f.relaxation_time - 0.5).abs() < 1e-12);
    }

    #[test]
    fn test_oldroyd_b_n1_after_shear() {
        let ucm = MaxwellFluid::new(0.001, 0.009, 1.0);
        let mut model = OldroydBModel::new(ucm, 0.1);
        // Apply simple shear repeatedly
        let kappa = [0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0f64];
        for _ in 0..20 {
            model.evolve_stress(&kappa, 0.05);
        }
        // N1 = tau_xx - tau_yy; under shear, tau_xx grows
        let _n1 = model.normal_stress_difference_1();
        // Just ensure no NaN
        assert!(!model.stress_tensor[0].is_nan());
    }

    #[test]
    fn test_deborah_large() {
        assert!(deborah_number(100.0, 1.0) > 1.0);
    }

    #[test]
    fn test_deborah_small() {
        assert!(deborah_number(0.001, 1.0) < 1.0);
    }
}

// ─── Oldroyd-B coupling regression tests ──────────────────────────────────────

#[cfg(test)]
mod tests_oldroyd_lbm {
    use super::*;

    /// Deterministic pseudo-random numbers in [0, 1) (64-bit LCG).
    struct Lcg(u64);
    impl Lcg {
        fn next(&mut self) -> f64 {
            self.0 = self
                .0
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            (self.0 >> 11) as f64 / (1u64 << 53) as f64
        }
    }

    fn equilibrium(rho: f64, ux: f64, uy: f64, i: usize) -> f64 {
        let cu = CX9[i] * ux + CY9[i] * uy;
        W9[i] * rho * (1.0 + 3.0 * cu + 4.5 * cu * cu - 1.5 * (ux * ux + uy * uy))
    }

    /// Fill every node with the equilibrium of the given velocity field.
    fn set_velocity_field(lbm: &mut ViscoelasticLbm, field: impl Fn(usize, usize) -> (f64, f64)) {
        for y in 0..lbm.ny {
            for x in 0..lbm.nx {
                let (ux, uy) = field(x, y);
                for i in 0..9 {
                    lbm.f_distributions[y * lbm.nx + x][i] = equilibrium(1.0, ux, uy, i);
                }
            }
        }
    }

    /// Total mass and raw first moment of all populations.
    fn mass_and_momentum(lbm: &ViscoelasticLbm) -> (f64, f64, f64) {
        let mut mass = 0.0;
        let mut px = 0.0;
        let mut py = 0.0;
        for f in &lbm.f_distributions {
            for i in 0..9 {
                mass += f[i];
                px += CX9[i] * f[i];
                py += CY9[i] * f[i];
            }
        }
        (mass, px, py)
    }

    fn max_abs_stress(lbm: &ViscoelasticLbm) -> f64 {
        lbm.stress_field
            .iter()
            .flat_map(|s| s.iter())
            .fold(0.0_f64, |m, v| m.max(v.abs()))
    }

    // ── Constitutive update ───────────────────────────────────────────────────

    #[test]
    fn steady_simple_shear_reaches_the_oldroyd_b_values() {
        let (lambda, eta_p, gamma_dot) = (20.0, 0.3, 0.01);
        // u_x(y): kappa_xy = du_x/dy.
        let mut kappa = [0.0; 9];
        kappa[1] = gamma_dot;
        let mut tau = [0.0; 9];
        for _ in 0..(40.0 * lambda) as usize {
            tau = oldroyd_b_stress_update(&tau, &kappa, lambda, eta_p, 1.0);
        }
        let shear = eta_p * gamma_dot;
        let n1 = 2.0 * eta_p * lambda * gamma_dot * gamma_dot;
        assert!((tau[1] - shear).abs() < 1e-6 * shear, "tau_xy {}", tau[1]);
        assert!((tau[3] - shear).abs() < 1e-6 * shear, "tau_yx {}", tau[3]);
        assert!((tau[0] - n1).abs() < 1e-6 * n1, "tau_xx {} vs {n1}", tau[0]);
        assert!(tau[4].abs() < 1e-15, "tau_yy {}", tau[4]);
        assert!(tau[8].abs() < 1e-15, "tau_zz {}", tau[8]);
        let n1_measured = tau[0] - tau[4];
        assert!((n1_measured - n1).abs() < 1e-6 * n1);
        let n2_measured = tau[4] - tau[8];
        assert!(n2_measured.abs() < 1e-15);
    }

    #[test]
    fn stretching_follows_the_velocity_gradient_convention() {
        // kappa_ij = du_i/dx_j: shearing u_y along x (kappa_yx) must stretch
        // along y; the transposed gradient would stretch along x instead.
        let (lambda, eta_p, gamma_dot) = (20.0, 0.3, 0.01);
        let mut kappa = [0.0; 9];
        kappa[3] = gamma_dot;
        let mut tau = [0.0; 9];
        for _ in 0..800 {
            tau = oldroyd_b_stress_update(&tau, &kappa, lambda, eta_p, 1.0);
        }
        let n1 = 2.0 * eta_p * lambda * gamma_dot * gamma_dot;
        assert!((tau[4] - n1).abs() < 1e-6 * n1, "tau_yy {} vs {n1}", tau[4]);
        assert!(tau[0].abs() < 1e-15, "tau_xx {}", tau[0]);
    }

    #[test]
    fn shear_start_up_matches_the_analytic_transient() {
        let (lambda, eta_p, gamma_dot) = (25.0, 0.2, 0.02);
        let mut kappa = [0.0; 9];
        kappa[1] = gamma_dot;
        let mut tau = [0.0; 9];
        let mut time = 0.0;
        for target in [10.0, 25.0, 60.0, 150.0] {
            while time < target {
                tau = oldroyd_b_stress_update(&tau, &kappa, lambda, eta_p, 1.0);
                time += 1.0;
            }
            let decay = (-time / lambda).exp();
            let shear = eta_p * gamma_dot * (1.0 - decay);
            let n1 = 2.0
                * eta_p
                * lambda
                * gamma_dot
                * gamma_dot
                * (1.0 - decay * (1.0 + time / lambda));
            assert!(
                (tau[1] - shear).abs() < 5e-4 * eta_p * gamma_dot,
                "t = {time}: tau_xy"
            );
            assert!(
                (tau[0] - n1).abs() < 5e-4 * 2.0 * eta_p * lambda * gamma_dot * gamma_dot,
                "t = {time}: N1 {} vs {n1}",
                tau[0]
            );
        }
    }

    #[test]
    fn zero_relaxation_time_is_the_newtonian_limit() {
        let mut rng = Lcg(5);
        let kappa: [f64; 9] = std::array::from_fn(|_| rng.next() - 0.5);
        let tau = [9.0; 9];
        let result = oldroyd_b_stress_update(&tau, &kappa, 0.0, 0.7, 1.0);
        for i in 0..3 {
            for j in 0..3 {
                let expected = 0.7 * (kappa[i * 3 + j] + kappa[j * 3 + i]);
                assert!((result[i * 3 + j] - expected).abs() < 1e-15);
            }
        }
        // A vanishing relaxation time approaches the same limit smoothly.
        let near = oldroyd_b_stress_update(&tau, &kappa, 1e-9, 0.7, 1.0);
        for k in 0..9 {
            assert!((near[k] - result[k]).abs() < 1e-6, "component {k}");
        }
    }

    #[test]
    fn relaxation_is_stable_for_steps_much_longer_than_lambda() {
        // lambda = 0.1 with unit steps (h/lambda = 10) used to be a stiff explicit step.
        let kappa = [0.0; 9];
        let mut tau: [f64; 9] = std::array::from_fn(|k| 1.0e3 * (k as f64 + 1.0));
        for _ in 0..6 {
            tau = oldroyd_b_stress_update(&tau, &kappa, 0.1, 0.5, 1.0);
        }
        assert!(tau.iter().all(|v| v.abs() < 1e-10), "{tau:?}");
        // The exact decay over one step is exp(-h/lambda).
        let once = oldroyd_b_stress_update(&[2.0; 9], &kappa, 4.0, 0.5, 1.0);
        for v in once {
            assert!((v - 2.0 * (-0.25_f64).exp()).abs() < 1e-14);
        }
    }

    #[test]
    fn phi_functions_match_their_taylor_series_on_both_branches() {
        // phi_n(z) = sum_k z^k / (k + n)!, evaluated with enough terms for |z| <= 3.
        let series = |z: f64, n: usize| {
            let mut term = (1..=n).fold(1.0, |acc, m| acc / m as f64);
            let mut sum = term;
            for k in 1..40 {
                term *= z / (k + n) as f64;
                sum += term;
            }
            sum
        };
        for z in [
            -3.0_f64, -0.2, -0.0500001, -0.0499999, -1e-6, 1e-6, 0.0499999, 0.0500001, 0.3,
        ] {
            assert!((phi1(z) - series(z, 1)).abs() < 1e-14, "phi1({z})");
            assert!((phi2(z) - series(z, 2)).abs() < 1e-14, "phi2({z})");
        }
        assert_eq!(phi1(0.0), 1.0);
        assert!((phi2(0.0) - 0.5).abs() < 1e-16);
    }

    // ── Advection of the stress ───────────────────────────────────────────────

    #[test]
    fn stress_pulse_advects_with_the_flow_without_new_extrema() {
        let (nx, ny) = (64, 3);
        let mut stress = vec![[0.0_f64; 9]; nx * ny];
        for y in 0..ny {
            for x in 0..nx {
                let d = x as f64 - 16.0;
                stress[y * nx + x][0] = (-d * d / (2.0 * 3.0 * 3.0)).exp();
            }
        }
        let ux = vec![0.1; nx * ny];
        let uy = vec![0.0; nx * ny];
        let total0: f64 = stress.iter().map(|s| s[0]).sum();
        let centroid = |s: &[[f64; 9]]| {
            let mut m0 = 0.0;
            let mut m1 = 0.0;
            for y in 0..ny {
                for x in 0..nx {
                    let value = s[y * nx + x][0];
                    m0 += value;
                    m1 += value * x as f64;
                }
            }
            m1 / m0
        };
        let c0 = centroid(&stress);
        for _ in 0..200 {
            advect_stress(&mut stress, &ux, &uy, nx, ny, 1.0);
        }
        let c1 = centroid(&stress);
        let total1: f64 = stress.iter().map(|s| s[0]).sum();
        assert!(
            (c1 - c0 - 20.0).abs() < 0.05,
            "pulse moved {} cells",
            c1 - c0
        );
        assert!(
            (total1 - total0).abs() < 1e-12 * total0,
            "total {total0} -> {total1}"
        );
        let peak = stress
            .iter()
            .map(|s| s[0])
            .fold(f64::NEG_INFINITY, f64::max);
        let trough = stress.iter().map(|s| s[0]).fold(f64::INFINITY, f64::min);
        assert!(peak <= 1.0 + 1e-12, "overshoot: {peak}");
        assert!(trough >= -1e-12, "undershoot: {trough}");
        assert!(
            peak > 0.9,
            "the second-order scheme must not flatten the pulse: {peak}"
        );
        // Components that were zero stay zero.
        assert!(stress.iter().all(|s| s[1] == 0.0 && s[4] == 0.0));
    }

    #[test]
    fn stress_advection_subcycles_when_the_courant_number_exceeds_one() {
        let (nx, ny) = (32, 4);
        let mut stress = vec![[0.0_f64; 9]; nx * ny];
        for y in 0..ny {
            for x in 8..16 {
                stress[y * nx + x][1] = 1.0;
            }
        }
        let ux = vec![0.7; nx * ny];
        let uy = vec![0.35; nx * ny];
        let total0: f64 = stress.iter().map(|s| s[1]).sum();
        for _ in 0..10 {
            advect_stress(&mut stress, &ux, &uy, nx, ny, 1.9);
        }
        let total1: f64 = stress.iter().map(|s| s[1]).sum();
        assert!((total1 - total0).abs() < 1e-10 * total0);
        let peak = stress
            .iter()
            .map(|s| s[1])
            .fold(f64::NEG_INFINITY, f64::max);
        let trough = stress.iter().map(|s| s[1]).fold(f64::INFINITY, f64::min);
        assert!(
            peak <= 1.0 + 1e-12 && trough >= -1e-12,
            "{trough} .. {peak}"
        );
    }

    #[test]
    fn advection_substeps_follow_the_courant_number_and_stay_bounded() {
        // Courant sum 0.4 is the limit of a single sub-step.
        assert_eq!(advection_substeps(0.1, 1.0), 1);
        assert_eq!(advection_substeps(0.4, 1.0), 1);
        assert_eq!(advection_substeps(0.41, 1.0), 2);
        assert_eq!(advection_substeps(1.0, 1.0), 3);
        assert_eq!(advection_substeps(1.05, 1.9), 5);
        assert_eq!(advection_substeps(0.2, 3.9), 2);
        // A blown-up velocity must not turn one step into billions of sub-steps.
        assert_eq!(advection_substeps(1.0e6, 1.0), MAX_ADVECTION_SUBSTEPS);
        assert_eq!(advection_substeps(f64::MAX, 1.0), MAX_ADVECTION_SUBSTEPS);
        assert_eq!(
            advection_substeps(f64::INFINITY, 1.0),
            MAX_ADVECTION_SUBSTEPS
        );
        // The bound is reached exactly at the Courant number it stands for.
        let at_bound = ADVECTION_CFL * MAX_ADVECTION_SUBSTEPS as f64;
        assert_eq!(advection_substeps(at_bound, 1.0), MAX_ADVECTION_SUBSTEPS);
        assert_eq!(
            advection_substeps(0.5 * at_bound, 1.0),
            MAX_ADVECTION_SUBSTEPS / 2
        );
    }

    #[test]
    fn stress_advection_of_a_diverged_velocity_field_terminates() {
        let (nx, ny) = (8, 8);
        let mut stress = vec![[0.0_f64; 9]; nx * ny];
        stress[3 * nx + 3][1] = 1.0;
        let ux = vec![1.0e30; nx * ny];
        let uy = vec![-2.0e30; nx * ny];
        // Without the bound this would sub-cycle `usize::MAX` times.
        advect_stress(&mut stress, &ux, &uy, nx, ny, 1.0);
        // Components that were zero everywhere are never touched.
        assert!(stress.iter().all(|s| s[0] == 0.0 && s[4] == 0.0));
    }

    // ── Coupled solver ────────────────────────────────────────────────────────

    #[test]
    fn uniform_translation_creates_no_polymer_stress() {
        // The old "shear rate" |u_x| + |u_y| turned any translation into stress.
        let mut lbm = ViscoelasticLbm::new(8, 6, 1.2, 5.0, 0.1);
        set_velocity_field(&mut lbm, |_, _| (0.05, -0.02));
        for _ in 0..60 {
            lbm.step();
        }
        assert!(
            max_abs_stress(&lbm) < 1e-16,
            "stress {}",
            max_abs_stress(&lbm)
        );
        for y in 0..lbm.ny {
            for x in 0..lbm.nx {
                let (ux, uy) = lbm.velocity(x, y);
                assert!((ux - 0.05).abs() < 1e-12 && (uy + 0.02).abs() < 1e-12);
                assert!((lbm.density(x, y) - 1.0).abs() < 1e-12);
            }
        }
    }

    #[test]
    fn solver_advects_the_stress_with_the_fluid() {
        // τ_zz does not enter the in-plane force and has no source term in a
        // uniform flow (κ = 0), so the solver only has to carry it along with the
        // fluid while it relaxes as exp(−t/λ).
        let (nx, ny) = (64, 3);
        let (omega, lambda, eta_p, speed) = (1.2, 500.0, 0.1, 0.1);
        let mut lbm = ViscoelasticLbm::new(nx, ny, omega, lambda, eta_p);
        set_velocity_field(&mut lbm, |_, _| (speed, 0.0));
        for y in 0..ny {
            for x in 0..nx {
                let d = x as f64 - 16.0;
                lbm.stress_field[y * nx + x][8] = (-d * d / (2.0 * 3.0 * 3.0)).exp();
            }
        }
        let moments = |lbm: &ViscoelasticLbm| {
            let mut total = 0.0;
            let mut first = 0.0;
            for y in 0..ny {
                for x in 0..nx {
                    let value = lbm.stress_field[y * nx + x][8];
                    total += value;
                    first += value * x as f64;
                }
            }
            (total, first / total)
        };
        let (total_start, centroid_start) = moments(&lbm);
        let steps = 100;
        for _ in 0..steps {
            lbm.step();
        }
        let (total_end, centroid_end) = moments(&lbm);
        let moved = centroid_end - centroid_start;
        assert!(
            (moved - speed * steps as f64).abs() < 0.05,
            "the pulse moved {moved} cells, expected u t = {}",
            speed * steps as f64
        );
        let expected_decay = (-(steps as f64) / lambda).exp();
        assert!(
            (total_end / total_start - expected_decay).abs() < 1e-9,
            "relaxation {} vs exp(-t/lambda) = {expected_decay}",
            total_end / total_start
        );
        // Nothing else was created, and the flow stayed uniform.
        assert!(
            lbm.stress_field
                .iter()
                .all(|s| s[..8].iter().all(|&v| v == 0.0))
        );
        for y in 0..ny {
            for x in 0..nx {
                let (ux, uy) = lbm.velocity(x, y);
                assert!((ux - speed).abs() < 1e-12 && uy.abs() < 1e-12);
            }
        }
    }

    #[test]
    fn rest_state_is_a_fixed_point() {
        let mut lbm = ViscoelasticLbm::new(5, 7, 1.7, 3.0, 0.4);
        for _ in 0..20 {
            lbm.step();
        }
        assert_eq!(max_abs_stress(&lbm), 0.0);
        for f in &lbm.f_distributions {
            for i in 0..9 {
                assert!((f[i] - W9[i]).abs() < 1e-14);
            }
        }
    }

    #[test]
    fn mass_and_momentum_are_conserved_for_an_arbitrary_stress_field() {
        // The old stress term injected (tau_xx + tau_yy) of mass per node and step.
        let (nx, ny) = (9, 7);
        let mut rng = Lcg(2024);
        let mut lbm = ViscoelasticLbm::new(nx, ny, 1.1, 6.0, 0.2);
        for f in lbm.f_distributions.iter_mut() {
            for (i, fi) in f.iter_mut().enumerate() {
                *fi = W9[i] * (0.9 + 0.2 * rng.next());
            }
        }
        for s in lbm.stress_field.iter_mut() {
            let (xx, xy, yy) = (
                0.01 * (rng.next() - 0.5),
                0.01 * (rng.next() - 0.5),
                0.01 * (rng.next() - 0.5),
            );
            *s = [xx, xy, 0.0, xy, yy, 0.0, 0.0, 0.0, 0.0];
        }
        let (m0, px0, py0) = mass_and_momentum(&lbm);
        for _ in 0..40 {
            lbm.step();
        }
        let (m1, px1, py1) = mass_and_momentum(&lbm);
        assert!((m1 - m0).abs() < 1e-12 * m0, "mass {m0} -> {m1}");
        assert!((px1 - px0).abs() < 1e-12 * m0, "x momentum {px0} -> {px1}");
        assert!((py1 - py0).abs() < 1e-12 * m0, "y momentum {py0} -> {py1}");
    }

    #[test]
    fn stress_force_is_the_central_difference_divergence() {
        let (nx, ny) = (6, 8);
        let mut lbm = ViscoelasticLbm::new(nx, ny, 1.0, 1.0, 0.1);
        let mut rng = Lcg(77);
        for s in lbm.stress_field.iter_mut() {
            let (xx, xy, yy) = (rng.next(), rng.next(), rng.next());
            *s = [xx, xy, 0.0, xy, yy, 0.0, 0.0, 0.0, 0.0];
        }
        let at = |x: isize, y: isize, c: usize| {
            lbm.stress_field[wrap_index(y, ny) * nx + wrap_index(x, nx)][c]
        };
        let mut sum = [0.0, 0.0];
        for y in 0..ny {
            for x in 0..nx {
                let (xi, yi) = (x as isize, y as isize);
                let fx = 0.5 * (at(xi + 1, yi, 0) - at(xi - 1, yi, 0))
                    + 0.5 * (at(xi, yi + 1, 3) - at(xi, yi - 1, 3));
                let fy = 0.5 * (at(xi + 1, yi, 1) - at(xi - 1, yi, 1))
                    + 0.5 * (at(xi, yi + 1, 4) - at(xi, yi - 1, 4));
                let force = lbm.stress_force(x, y);
                assert!((force[0] - fx).abs() < 1e-15 && (force[1] - fy).abs() < 1e-15);
                sum[0] += force[0];
                sum[1] += force[1];
            }
        }
        assert!(
            sum[0].abs() < 1e-13 && sum[1].abs() < 1e-13,
            "net force {sum:?}"
        );
    }

    #[test]
    fn velocity_accessor_includes_the_half_force() {
        let (nx, ny) = (4, 8);
        let mut lbm = ViscoelasticLbm::new(nx, ny, 1.0, 10.0, 0.1);
        for y in 0..ny {
            for x in 0..nx {
                lbm.stress_field[y * nx + x][1] = 0.01 * (y as f64).sin();
                lbm.stress_field[y * nx + x][3] = 0.01 * (y as f64).sin();
            }
        }
        for y in 0..ny {
            let force = lbm.stress_force(1, y);
            assert!(force[0].abs() > 1e-4, "set-up must produce a force");
            let (ux, uy) = lbm.velocity(1, y);
            assert!((ux - 0.5 * force[0]).abs() < 1e-15);
            assert!((uy - 0.5 * force[1]).abs() < 1e-15);
        }
    }

    #[test]
    fn newtonian_limit_stress_follows_the_strain_rate() {
        // lambda = 0: tau = eta_p (kappa + kappa^T) from the current velocity.
        let (nx, ny) = (4, 32);
        let k = 2.0 * std::f64::consts::PI / ny as f64;
        let mut lbm = ViscoelasticLbm::new(nx, ny, 1.0, 0.0, 0.25);
        set_velocity_field(&mut lbm, |_, y| (0.01 * (k * y as f64).sin(), 0.0));
        lbm.step();
        for y in 0..ny {
            let up = 0.01 * (k * (y as f64 + 1.0)).sin();
            let down = 0.01 * (k * (y as f64 - 1.0)).sin();
            let gamma = 0.5 * (up - down);
            let tau = lbm.stress_field[y * nx];
            assert!(
                (tau[1] - 0.25 * gamma).abs() < 1e-9,
                "row {y}: tau_xy {}",
                tau[1]
            );
            assert!((tau[3] - tau[1]).abs() < 1e-15);
            assert!(tau[0].abs() < 1e-12 && tau[4].abs() < 1e-12);
        }
    }

    /// Fourth-order Runge-Kutta for the linearised shear-wave pair (U, T):
    /// dU/dt = -nu_s k^2 U - k T,  dT/dt = (eta_p k U - T) / lambda.
    fn shear_wave_reference(
        nu_s: f64,
        k: f64,
        eta_p: f64,
        lambda: f64,
        u0: f64,
        t_end: f64,
    ) -> [f64; 2] {
        let rhs = |s: [f64; 2]| {
            [
                -nu_s * k * k * s[0] - k * s[1],
                (eta_p * k * s[0] - s[1]) / lambda,
            ]
        };
        let steps = (t_end / 0.05).round() as usize;
        let h = t_end / steps as f64;
        let mut s = [u0, 0.0];
        for _ in 0..steps {
            let k1 = rhs(s);
            let k2 = rhs([s[0] + 0.5 * h * k1[0], s[1] + 0.5 * h * k1[1]]);
            let k3 = rhs([s[0] + 0.5 * h * k2[0], s[1] + 0.5 * h * k2[1]]);
            let k4 = rhs([s[0] + h * k3[0], s[1] + h * k3[1]]);
            for c in 0..2 {
                s[c] += h / 6.0 * (k1[c] + 2.0 * k2[c] + 2.0 * k3[c] + k4[c]);
            }
        }
        s
    }

    #[test]
    fn shear_wave_follows_the_linear_oldroyd_b_solution() {
        // u_x = U sin(k y), tau_xy = T cos(k y).  The elastic restoring force
        // turns the decay into a damped oscillation that changes sign; neither
        // a Newtonian fluid nor a mass-source stress term reproduces it.
        let (nx, ny) = (1, 64);
        let (omega, lambda, eta_p, amp) = (1.8, 300.0, 0.5, 0.02);
        let nu_s = (1.0 / omega - 0.5) / 3.0;
        let k = 2.0 * std::f64::consts::PI / ny as f64;
        let mut lbm = ViscoelasticLbm::new(nx, ny, omega, lambda, eta_p);
        set_velocity_field(&mut lbm, |_, y| (amp * (k * y as f64).sin(), 0.0));
        let project = |lbm: &ViscoelasticLbm| {
            let mut u_amp = 0.0;
            let mut t_amp = 0.0;
            for y in 0..ny {
                u_amp += lbm.velocity(0, y).0 * (k * y as f64).sin();
                t_amp += lbm.stress_field[y * nx][1] * (k * y as f64).cos();
            }
            [2.0 * u_amp / ny as f64, 2.0 * t_amp / ny as f64]
        };
        let start = project(&lbm);
        assert!((start[0] - amp).abs() < 1e-12 && start[1].abs() < 1e-15);
        let mut time = 0;
        let mut most_negative = 0.0_f64;
        for target in [150, 300, 450, 600, 750, 900, 1050, 1200, 1350] {
            while time < target {
                lbm.step();
                time += 1;
            }
            let reference = shear_wave_reference(nu_s, k, eta_p, lambda, amp, time as f64);
            let measured = project(&lbm);
            most_negative = most_negative.min(reference[0]);
            assert!(
                (measured[0] - reference[0]).abs() < 0.005 * amp,
                "t = {time}: U {} vs {}",
                measured[0],
                reference[0]
            );
            let t_scale = eta_p * k * amp;
            assert!(
                (measured[1] - reference[1]).abs() < 0.01 * t_scale,
                "t = {time}: T {} vs {}",
                measured[1],
                reference[1]
            );
        }
        // The reference oscillates through zero; a Newtonian decay cannot.
        assert!(
            most_negative < -0.05 * amp,
            "reference minimum {most_negative}"
        );
        let newtonian = amp * (-nu_s * k * k * 1050.0).exp();
        assert!(newtonian > 0.5 * amp);
    }

    #[test]
    fn solver_applies_the_constitutive_update_with_its_own_velocity_gradient() {
        // Replay the stress ODE of one row with the velocity the solver sees:
        // the nonlinear stretching (N1) must come out of the solver, not only
        // out of the helper.
        let (nx, ny) = (1, 48);
        let (omega, lambda, eta_p, amp) = (1.6, 40.0, 0.4, 0.06);
        let k = 2.0 * std::f64::consts::PI / ny as f64;
        let mut lbm = ViscoelasticLbm::new(nx, ny, omega, lambda, eta_p);
        set_velocity_field(&mut lbm, |_, y| (amp * (k * y as f64).sin(), 0.0));
        let mut shadow = vec![[0.0_f64; 9]; ny];
        for _ in 0..150 {
            let u: Vec<f64> = (0..ny).map(|y| lbm.velocity(0, y).0).collect();
            for (y, tau) in shadow.iter_mut().enumerate() {
                let mut kappa = [0.0; 9];
                kappa[1] = 0.5 * (u[(y + 1) % ny] - u[(y + ny - 1) % ny]);
                *tau = oldroyd_b_stress_update(tau, &kappa, lambda, eta_p, 1.0);
            }
            lbm.step();
        }
        let n1_scale = 2.0 * eta_p * lambda * (amp * k).powi(2);
        let mut peak_n1 = 0.0_f64;
        for (y, expected) in shadow.iter().enumerate() {
            let tau = &lbm.stress_field[y * nx];
            for c in 0..9 {
                assert!(
                    (tau[c] - expected[c]).abs() < 1e-9 * (1.0 + n1_scale),
                    "row {y}, component {c}: {} vs {}",
                    tau[c],
                    expected[c]
                );
            }
            peak_n1 = peak_n1.max(tau[0] - tau[4]);
        }
        assert!(
            peak_n1 > 0.05 * n1_scale,
            "N1 must develop: {peak_n1} vs {n1_scale}"
        );
        assert!(lbm.stress_field.iter().all(|s| s[4].abs() < 1e-14));
    }

    /// Kinetic energy `Σ ρ|u|²/2` (with the half-force velocity) of the solver.
    fn kinetic_energy(lbm: &ViscoelasticLbm) -> f64 {
        let mut energy = 0.0;
        for y in 0..lbm.ny {
            for x in 0..lbm.nx {
                let (ux, uy) = lbm.velocity(x, y);
                energy += 0.5 * lbm.density(x, y) * (ux * ux + uy * uy);
            }
        }
        energy
    }

    #[test]
    fn short_relaxation_time_coupling_stays_stable() {
        // ω = 1.4 (ν_s = 0.071), λ = 0.3, η_p = 0.2: the polymer viscosity is
        // 2.8 times the solvent viscosity and the stress follows the strain rate
        // almost instantly.  Taking the force from the stress of the previous
        // step made this vortex diverge (the kinetic energy grew by orders of
        // magnitude and the run ended in NaN); the force of the updated stress
        // keeps it decaying.
        let n = 24;
        let k = 2.0 * std::f64::consts::PI / n as f64;
        let amp = 0.05;
        let mut lbm = ViscoelasticLbm::new(n, n, 1.4, 0.3, 0.2);
        set_velocity_field(&mut lbm, |x, y| {
            let (x, y) = (x as f64, y as f64);
            (
                -amp * (k * x).cos() * (k * y).sin(),
                amp * (k * x).sin() * (k * y).cos(),
            )
        });
        let energy0 = kinetic_energy(&lbm);
        assert!(energy0 > 0.0);
        for step in 0..600 {
            lbm.step();
            let energy = kinetic_energy(&lbm);
            assert!(
                energy.is_finite() && energy <= energy0,
                "step {step}: kinetic energy {energy} (initial {energy0})"
            );
        }
        assert!(kinetic_energy(&lbm) < 1e-3 * energy0);
        for y in 0..n {
            for x in 0..n {
                assert!((lbm.density(x, y) - 1.0).abs() < 1e-3);
            }
        }
    }

    #[test]
    fn newtonian_limit_adds_the_polymer_viscosity_to_the_shear_decay() {
        // λ = 0 gives τ = η_p (κ + κᵀ), so a shear wave decays with the total
        // viscosity ν_s + η_p: ν_s k² from the lattice Boltzmann fluid and
        // η_p sin²k from the central-difference stress force.
        let (nx, ny) = (1, 64);
        let (omega, eta_p, amp) = (1.8, 0.1, 0.01);
        let nu_s = (1.0 / omega - 0.5) / 3.0;
        let k = 2.0 * std::f64::consts::PI / ny as f64;
        let mut lbm = ViscoelasticLbm::new(nx, ny, omega, 0.0, eta_p);
        set_velocity_field(&mut lbm, |_, y| (amp * (k * y as f64).sin(), 0.0));
        let amplitude = |lbm: &ViscoelasticLbm| {
            let projection: f64 = (0..ny)
                .map(|y| lbm.velocity(0, y).0 * (k * y as f64).sin())
                .sum();
            2.0 * projection / ny as f64
        };
        for _ in 0..100 {
            lbm.step();
        }
        let early = amplitude(&lbm);
        let steps = 400;
        for _ in 0..steps {
            lbm.step();
        }
        let late = amplitude(&lbm);
        let measured = (early / late).ln() / steps as f64;
        let expected = nu_s * k * k + eta_p * k.sin().powi(2);
        assert!(
            (measured - expected).abs() < 0.01 * expected,
            "decay rate {measured} vs (nu_s + eta_p) k^2 = {expected}"
        );
        // The solvent alone would decay several times more slowly.
        assert!(measured > 4.0 * nu_s * k * k);
    }
}
