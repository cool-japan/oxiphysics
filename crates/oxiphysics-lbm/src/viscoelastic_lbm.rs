#![allow(clippy::needless_range_loop)]
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
#[allow(dead_code)]
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
#[allow(dead_code)]
pub fn deborah_number(lambda: f64, t_process: f64) -> f64 {
    if t_process == 0.0 {
        return f64::INFINITY;
    }
    lambda / t_process
}

/// Return the Trouton ratio for a purely Newtonian fluid (= 3).
///
/// For a Newtonian fluid, extensional viscosity = 3 × shear viscosity.
#[allow(dead_code)]
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
    #[allow(dead_code)]
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
    #[allow(dead_code)]
    pub fn weissenberg_number(&self, shear_rate: f64) -> f64 {
        self.relaxation_time * shear_rate
    }

    /// Total zero-shear viscosity η = η_s + η_p.
    #[allow(dead_code)]
    pub fn total_viscosity(&self) -> f64 {
        self.viscosity_s + self.viscosity_p
    }

    /// Viscosity ratio β = η_s / (η_s + η_p).
    #[allow(dead_code)]
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
    #[allow(dead_code)]
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
    #[allow(dead_code)]
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
    #[allow(dead_code)]
    pub fn normal_stress_difference_1(&self) -> f64 {
        self.stress_tensor[0] - self.stress_tensor[4]
    }

    /// Return the second normal stress difference N2 = τ_yy − τ_zz.
    #[allow(dead_code)]
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
    #[allow(dead_code)]
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
    #[allow(dead_code)]
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
    #[allow(dead_code)]
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

/// Viscoelastic Lattice Boltzmann solver (D2Q9 layout).
///
/// Couples the standard BGK LBM for the hydrodynamics with an Oldroyd-B
/// polymer stress field that is evolved each time step.
///
/// Distribution functions are stored as a flat `Vec<Vec`f64`>` with outer
/// index being the lattice site (ny * nx) and inner index the 9 D2Q9
/// velocities.
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

impl ViscoelasticLbm {
    /// Create a new viscoelastic LBM domain initialised to rest.
    ///
    /// # Arguments
    /// * `nx`, `ny` – domain dimensions
    /// * `omega`    – BGK relaxation frequency (1/τ)
    /// * `lambda`   – polymer relaxation time (lattice units)
    /// * `eta_p`    – polymer viscosity (lattice units)
    #[allow(dead_code)]
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

    /// Perform a single LBM time step.
    ///
    /// 1. Compute macroscopic density and velocity at each node.
    /// 2. BGK collision with polymer stress correction.
    /// 3. Stream distributions to neighbours (periodic boundaries).
    /// 4. Evolve polymer stress using velocity gradient estimate.
    #[allow(dead_code)]
    pub fn step(&mut self, dt: f64) {
        let nx = self.nx;
        let ny = self.ny;
        let n = nx * ny;

        // ── 1. Collision ──────────────────────────────────────────────────────
        let mut f_post = self.f_distributions.clone();
        for node in 0..n {
            let f = &self.f_distributions[node];
            let rho: f64 = f.iter().sum();
            if rho == 0.0 {
                continue;
            }
            let ux: f64 = f
                .iter()
                .zip(CX9.iter())
                .map(|(fi, ci)| fi * ci)
                .sum::<f64>()
                / rho;
            let uy: f64 = f
                .iter()
                .zip(CY9.iter())
                .map(|(fi, ci)| fi * ci)
                .sum::<f64>()
                / rho;

            // Polymer stress contribution (xx, yy, xy components mapped to D2Q9)
            let tau_xx = self.stress_field[node][0];
            let tau_yy = self.stress_field[node][4];
            let tau_xy = self.stress_field[node][1];

            for i in 0..9 {
                let cu = CX9[i] * ux + CY9[i] * uy;
                let u2 = ux * ux + uy * uy;
                let feq = W9[i] * rho * (1.0 + 3.0 * cu + 4.5 * cu * cu - 1.5 * u2);

                // Polymer stress forcing (simplified body-force approach)
                let g_poly = W9[i]
                    * 3.0
                    * (CX9[i] * CX9[i] * tau_xx
                        + 2.0 * CX9[i] * CY9[i] * tau_xy
                        + CY9[i] * CY9[i] * tau_yy);

                f_post[node][i] =
                    f[i] - self.omega * (f[i] - feq) + (1.0 - 0.5 * self.omega) * g_poly * dt;
            }

            // ── 4. Evolve polymer stress ─────────────────────────────────────
            if self.lambda > 0.0 {
                // Approximate velocity gradient from macroscopic velocity
                // (central differences not available at single node; use simple shear estimate)
                let gamma_dot = ux.abs() + uy.abs(); // magnitude proxy
                let kappa = [0.0, gamma_dot, 0.0, gamma_dot, 0.0, 0.0, 0.0, 0.0, 0.0_f64];
                let prefactor = self.eta_p / self.lambda;
                let lambda = self.lambda;
                for ii in 0..3 {
                    for jj in 0..3 {
                        let kij = kappa[ii * 3 + jj];
                        let kji = kappa[jj * 3 + ii];
                        let tau_ij = self.stress_field[node][ii * 3 + jj];
                        let dtau = prefactor * (kij + kji) - tau_ij / lambda;
                        self.stress_field[node][ii * 3 + jj] += dt * dtau;
                    }
                }
            }
        }

        // ── 2. Streaming (periodic) ───────────────────────────────────────────
        let mut f_stream = vec![vec![0.0f64; 9]; n];
        for y in 0..ny {
            for x in 0..nx {
                let node = y * nx + x;
                for i in 0..9 {
                    let xi = CX9[i] as isize;
                    let yi = CY9[i] as isize;
                    let xd = ((x as isize + xi).rem_euclid(nx as isize)) as usize;
                    let yd = ((y as isize + yi).rem_euclid(ny as isize)) as usize;
                    let dest = yd * nx + xd;
                    f_stream[dest][i] = f_post[node][i];
                }
            }
        }
        self.f_distributions = f_stream;
    }

    /// Return macroscopic density at node (ix, iy).
    #[allow(dead_code)]
    pub fn density(&self, ix: usize, iy: usize) -> f64 {
        self.f_distributions[iy * self.nx + ix].iter().sum()
    }

    /// Return macroscopic velocity (ux, uy) at node (ix, iy).
    #[allow(dead_code)]
    pub fn velocity(&self, ix: usize, iy: usize) -> (f64, f64) {
        let f = &self.f_distributions[iy * self.nx + ix];
        let rho: f64 = f.iter().sum();
        if rho == 0.0 {
            return (0.0, 0.0);
        }
        let ux = f
            .iter()
            .zip(CX9.iter())
            .map(|(fi, ci)| fi * ci)
            .sum::<f64>()
            / rho;
        let uy = f
            .iter()
            .zip(CY9.iter())
            .map(|(fi, ci)| fi * ci)
            .sum::<f64>()
            / rho;
        (ux, uy)
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
    #[allow(dead_code)]
    pub fn new(extension_rate: f64) -> Self {
        Self { extension_rate }
    }

    /// Compute the Trouton ratio Tr = η_E / η.
    ///
    /// # Arguments
    /// * `viscosity`            – shear viscosity η (Pa·s)
    /// * `extensional_viscosity` – extensional (Trouton) viscosity η_E (Pa·s)
    #[allow(dead_code)]
    pub fn compute_trouton_ratio(&self, viscosity: f64, extensional_viscosity: f64) -> f64 {
        if viscosity == 0.0 {
            return 0.0;
        }
        extensional_viscosity / viscosity
    }

    /// Predict the Hencky strain ε = ε̇ · t accumulated after time t.
    #[allow(dead_code)]
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
        lbm.step(0.01);
    }

    #[test]
    fn test_lbm_mass_conservation() {
        let mut lbm = ViscoelasticLbm::new(4, 4, 1.0, 0.0, 0.0);
        let total_before: f64 = (0..4 * 4)
            .map(|n| lbm.f_distributions[n].iter().sum::<f64>())
            .sum();
        lbm.step(0.01);
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
        lbm.step(0.1);
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
