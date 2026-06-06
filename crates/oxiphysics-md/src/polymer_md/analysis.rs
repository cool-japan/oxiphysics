// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Analysis and constraint models: KremerGrestBond, ShakeConstraints,
//! PolymerPressureTensor, EndToEndCorrelation, RadiusOfGyrationTracker,
//! BondBreakingKinetics, RadiusOfGyration, PolymerEntanglement, DendrimerModel.

use super::{dot3_arr, sub3_arr};
use std::f64::consts::PI;

/// Kremer-Grest bead-spring model: FENE + WCA combined bond potential.
///
/// The standard model for coarse-grained polymer simulations.
///
/// References:
/// - Kremer, K. & Grest, G.S. (1990). *J. Chem. Phys.* 92, 5057.
#[derive(Debug, Clone)]
pub struct KremerGrestBond {
    /// FENE spring constant K (in kT/σ²).
    pub k_fene: f64,
    /// Maximum extension R₀ (in σ).
    pub r0: f64,
    /// LJ energy parameter ε.
    pub epsilon_lj: f64,
    /// LJ size parameter σ.
    pub sigma_lj: f64,
}

impl KremerGrestBond {
    /// Create a Kremer-Grest bond with standard parameters.
    ///
    /// Default: K=30, R₀=1.5σ, ε=1.0, σ=1.0.
    pub fn new(k_fene: f64, r0: f64, epsilon_lj: f64, sigma_lj: f64) -> Self {
        Self {
            k_fene,
            r0,
            epsilon_lj,
            sigma_lj,
        }
    }

    /// FENE potential energy at bond extension r.
    ///
    /// U_FENE(r) = -½ K R₀² ln(1 - (r/R₀)²)
    pub fn fene_potential(&self, r: f64) -> f64 {
        let x = r / self.r0;
        if x >= 1.0 {
            return f64::INFINITY;
        }
        -0.5 * self.k_fene * self.r0 * self.r0 * (1.0 - x * x).ln()
    }

    /// WCA (truncated + shifted LJ) potential.
    ///
    /// U_WCA = 4ε\[(σ/r)¹² - (σ/r)⁶\] + ε for r < r_cut = 2^(1/6)σ, else 0.
    pub fn wca_potential(&self, r: f64) -> f64 {
        let r_cut = 2.0_f64.powf(1.0 / 6.0) * self.sigma_lj;
        if r >= r_cut {
            return 0.0;
        }
        let s6 = (self.sigma_lj / r.max(1e-30)).powi(6);
        4.0 * self.epsilon_lj * (s6 * s6 - s6) + self.epsilon_lj
    }

    /// Total Kremer-Grest potential U = U_FENE + U_WCA.
    pub fn total_potential(&self, r: f64) -> f64 {
        self.fene_potential(r) + self.wca_potential(r)
    }

    /// FENE force magnitude (restoring): F_FENE = K r / (1 - (r/R₀)²).
    pub fn fene_force_magnitude(&self, r: f64) -> f64 {
        let x = r / self.r0;
        if x >= 1.0 {
            return f64::INFINITY;
        }
        self.k_fene * r / (1.0 - x * x)
    }

    /// WCA force magnitude.
    pub fn wca_force_magnitude(&self, r: f64) -> f64 {
        let r_cut = 2.0_f64.powf(1.0 / 6.0) * self.sigma_lj;
        if r >= r_cut {
            return 0.0;
        }
        let r_safe = r.max(1e-30);
        let s6 = (self.sigma_lj / r_safe).powi(6);
        24.0 * self.epsilon_lj / r_safe * (2.0 * s6 * s6 - s6)
    }

    /// Total bond force magnitude.
    pub fn total_force_magnitude(&self, r: f64) -> f64 {
        self.fene_force_magnitude(r) + self.wca_force_magnitude(r)
    }

    /// Equilibrium separation r_eq (minimum of total potential).
    ///
    /// For pure WCA r_eq = 2^(1/6) σ. With FENE the minimum shifts slightly.
    pub fn equilibrium_approx(&self) -> f64 {
        2.0_f64.powf(1.0 / 6.0) * self.sigma_lj
    }

    /// Check if bond is over-extended (r > 0.9 R₀).
    pub fn is_over_extended(&self, r: f64) -> bool {
        r > 0.9 * self.r0
    }
}

/// Bond constraint solver using SHAKE algorithm.
///
/// Enforces rigid bond length constraints using iterative Lagrange multipliers.
/// Implements the SHAKE algorithm (Ryckaert et al., 1977).
///
/// References:
/// - Ryckaert, J.-P. et al. (1977). *J. Comput. Phys.* 23, 327.
/// - Andersen, H.C. (1983). RATTLE. *J. Comput. Phys.* 52, 24.
#[derive(Debug, Clone)]
pub struct ShakeConstraints {
    /// Number of beads.
    pub n_beads: usize,
    /// Target bond length d₀.
    pub d0: f64,
    /// Convergence tolerance.
    pub tol: f64,
    /// Maximum iterations.
    pub max_iter: usize,
    /// Bead masses.
    pub masses: Vec<f64>,
}

impl ShakeConstraints {
    /// Create a SHAKE constraint solver.
    pub fn new(n_beads: usize, d0: f64, tol: f64, max_iter: usize) -> Self {
        let masses = vec![1.0f64; n_beads];
        Self {
            n_beads,
            d0,
            tol,
            max_iter,
            masses,
        }
    }

    /// Apply SHAKE constraints to positions given unconstrained new positions.
    ///
    /// Iterates until all bond lengths satisfy |r_{i+1} - r_i| = d₀ ± tol.
    pub fn apply(&self, positions: &mut [[f64; 3]], old_positions: &[[f64; 3]]) -> usize {
        let n = self.n_beads;
        let mut iter = 0usize;
        for _k in 0..self.max_iter {
            iter = _k + 1;
            let mut converged = true;
            for i in 0..n.saturating_sub(1) {
                let r12 = sub3_arr(positions[i + 1], positions[i]);
                let d_sq = dot3_arr(r12, r12);
                let d0_sq = self.d0 * self.d0;
                let err = (d_sq - d0_sq).abs();
                if err > self.tol * self.tol {
                    converged = false;
                    // SHAKE correction using old bond vector
                    let r12_old = sub3_arr(old_positions[i + 1], old_positions[i]);
                    let dot = dot3_arr(r12, r12_old);
                    if dot.abs() < 1e-30 {
                        continue;
                    }
                    let m_i = self.masses[i];
                    let m_j = self.masses[i + 1];
                    let lambda = (d0_sq - d_sq) / (2.0 * dot * (1.0 / m_i + 1.0 / m_j));
                    for d in 0..3 {
                        positions[i][d] -= lambda / m_i * r12_old[d];
                        positions[i + 1][d] += lambda / m_j * r12_old[d];
                    }
                }
            }
            if converged {
                break;
            }
        }
        iter
    }

    /// Compute bond length violations.
    pub fn bond_violations(&self, positions: &[[f64; 3]]) -> Vec<f64> {
        (0..self.n_beads.saturating_sub(1))
            .map(|i| {
                let r = sub3_arr(positions[i + 1], positions[i]);
                let d = dot3_arr(r, r).sqrt();
                (d - self.d0).abs()
            })
            .collect()
    }

    /// Maximum bond violation.
    pub fn max_violation(&self, positions: &[[f64; 3]]) -> f64 {
        self.bond_violations(positions)
            .iter()
            .cloned()
            .fold(0.0f64, f64::max)
    }
}

/// Polymer pressure tensor via the Kramers-Kirkwood / Rouse virial.
///
/// Computes the stress contribution from polymer chain bonds using
/// the Kramers-Kirkwood formula:
/// P_αβ = -ρ kT δ_αβ + (1/V) Σ_bonds r_α F_β
#[derive(Debug, Clone)]
pub struct PolymerPressureTensor {
    /// Number of beads N.
    pub n_beads: usize,
    /// System volume V.
    pub volume: f64,
    /// Thermal energy kT.
    pub kt: f64,
    /// Density ρ = N / V.
    pub density: f64,
}

impl PolymerPressureTensor {
    /// Create a polymer pressure tensor calculator.
    pub fn new(n_beads: usize, volume: f64, kt: f64) -> Self {
        let density = n_beads as f64 / volume;
        Self {
            n_beads,
            volume,
            kt,
            density,
        }
    }

    /// Kinetic (ideal gas) contribution to pressure tensor: P_kin = ρ kT δ_αβ.
    pub fn kinetic_pressure(&self) -> f64 {
        self.density * self.kt
    }

    /// Virial contribution from bond i→i+1: σ_αβ += r_α * F_β.
    ///
    /// Returns \[P_xx, P_xy, P_xz, P_yy, P_yz, P_zz\] bond contribution.
    pub fn bond_virial(&self, r_ij: [f64; 3], force: [f64; 3]) -> [f64; 6] {
        let scale = 1.0 / self.volume;
        [
            scale * r_ij[0] * force[0], // xx
            scale * r_ij[0] * force[1], // xy
            scale * r_ij[0] * force[2], // xz
            scale * r_ij[1] * force[1], // yy
            scale * r_ij[1] * force[2], // yz
            scale * r_ij[2] * force[2], // zz
        ]
    }

    /// Total Kramers-Kirkwood virial stress for a chain.
    ///
    /// P_αβ = (1/V) Σ_{bonds} r_α^(bond) F_β^(FENE)
    pub fn kramers_kirkwood_virial(
        &self,
        positions: &[[f64; 3]],
        k_fene: f64,
        r0_fene: f64,
    ) -> [f64; 6] {
        let mut sigma = [0.0f64; 6];
        let n = positions.len();
        for i in 0..n.saturating_sub(1) {
            let r_ij = sub3_arr(positions[i + 1], positions[i]);
            let rmag = dot3_arr(r_ij, r_ij).sqrt().max(1e-30);
            let x = rmag / r0_fene;
            if x >= 1.0 {
                continue;
            }
            let fmag = k_fene * rmag / (1.0 - x * x);
            let rhat = [r_ij[0] / rmag, r_ij[1] / rmag, r_ij[2] / rmag];
            let force = [fmag * rhat[0], fmag * rhat[1], fmag * rhat[2]];
            let contrib = self.bond_virial(r_ij, force);
            for d in 0..6 {
                sigma[d] += contrib[d];
            }
        }
        sigma
    }

    /// Isotropic pressure P = (P_xx + P_yy + P_zz) / 3 + P_kin.
    pub fn isotropic_pressure(&self, virial: &[f64; 6]) -> f64 {
        let p_vir = (virial[0] + virial[3] + virial[5]) / 3.0;
        self.kinetic_pressure() + p_vir
    }

    /// Shear stress: σ_xy (off-diagonal element).
    pub fn shear_stress(&self, virial: &[f64; 6]) -> f64 {
        virial[1]
    }
}

/// End-to-end distance correlation function tracker.
///
/// Computes C(t) = <R(t)·R(0)> / <R²(0)> to characterize chain dynamics.
#[derive(Debug, Clone)]
pub struct EndToEndCorrelation {
    /// Chain length N.
    pub n_beads: usize,
    /// Segment length b.
    pub b: f64,
    /// Thermal energy kT.
    pub kt: f64,
    /// Stored end-to-end vector history.
    pub r_history: Vec<[f64; 3]>,
    /// Maximum history length.
    pub max_history: usize,
}

impl EndToEndCorrelation {
    /// Create a new end-to-end correlation tracker.
    pub fn new(n_beads: usize, b: f64, kt: f64, max_history: usize) -> Self {
        Self {
            n_beads,
            b,
            kt,
            r_history: Vec::new(),
            max_history,
        }
    }

    /// Record the current end-to-end vector.
    pub fn record(&mut self, r_end: [f64; 3]) {
        if self.r_history.len() >= self.max_history {
            self.r_history.remove(0);
        }
        self.r_history.push(r_end);
    }

    /// Compute the correlation at lag τ (in stored steps).
    ///
    /// C(τ) = <R(t+τ)·R(t)> / <R(0)²>
    pub fn correlation(&self, lag: usize) -> f64 {
        let n = self.r_history.len();
        if n == 0 || lag >= n {
            return 0.0;
        }
        let r0 = self.r_history[0];
        let r0_sq = dot3_arr(r0, r0).max(1e-30);
        let count = n - lag;
        let sum: f64 = (0..count)
            .map(|t| dot3_arr(self.r_history[t], self.r_history[t + lag]))
            .sum();
        sum / (count as f64 * r0_sq)
    }

    /// Theoretical Rouse relaxation: C(t) = exp(-t / τ₁).
    pub fn rouse_theory(&self, t: f64, tau1: f64) -> f64 {
        (-t / tau1.max(1e-30)).exp()
    }

    /// Mean-square end-to-end distance from history.
    pub fn mean_r2(&self) -> f64 {
        if self.r_history.is_empty() {
            return 0.0;
        }
        let sum: f64 = self.r_history.iter().map(|&r| dot3_arr(r, r)).sum();
        sum / self.r_history.len() as f64
    }

    /// Theoretical `R²` = N b² (ideal chain).
    pub fn ideal_r2_theory(&self) -> f64 {
        self.n_beads as f64 * self.b * self.b
    }
}

/// Radius of gyration tracker with running statistics.
///
/// Tracks Rg over multiple configurations and provides statistics.
#[derive(Debug, Clone)]
pub struct RadiusOfGyrationTracker {
    /// Recorded Rg values.
    pub rg_values: Vec<f64>,
    /// Theoretical Rg scaling exponent ν.
    pub nu: f64,
    /// Bond length b.
    pub b: f64,
    /// Number of beads N.
    pub n_beads: usize,
}

impl RadiusOfGyrationTracker {
    /// Create a new Rg tracker.
    pub fn new(n_beads: usize, b: f64, nu: f64) -> Self {
        Self {
            rg_values: Vec::new(),
            nu,
            b,
            n_beads,
        }
    }

    /// Compute Rg from a configuration.
    ///
    /// Rg² = (1/N) Σ_i |r_i - r_cm|²
    pub fn compute_rg(positions: &[[f64; 3]]) -> f64 {
        let n = positions.len();
        if n == 0 {
            return 0.0;
        }
        let nf = n as f64;
        let cm: [f64; 3] = positions.iter().fold([0.0; 3], |acc, &p| {
            [acc[0] + p[0] / nf, acc[1] + p[1] / nf, acc[2] + p[2] / nf]
        });
        let rg2: f64 = positions
            .iter()
            .map(|&p| {
                let d = sub3_arr(p, cm);
                dot3_arr(d, d)
            })
            .sum::<f64>()
            / nf;
        rg2.sqrt()
    }

    /// Record a Rg measurement.
    pub fn record(&mut self, positions: &[[f64; 3]]) {
        self.rg_values.push(Self::compute_rg(positions));
    }

    /// Mean Rg over all recorded configurations.
    pub fn mean_rg(&self) -> f64 {
        if self.rg_values.is_empty() {
            return 0.0;
        }
        self.rg_values.iter().sum::<f64>() / self.rg_values.len() as f64
    }

    /// Variance of Rg.
    pub fn variance_rg(&self) -> f64 {
        let mean = self.mean_rg();
        if self.rg_values.len() < 2 {
            return 0.0;
        }
        self.rg_values
            .iter()
            .map(|&x| (x - mean).powi(2))
            .sum::<f64>()
            / (self.rg_values.len() - 1) as f64
    }

    /// Theoretical Rg = b N^ν / sqrt(6) for Gaussian chain (ν=0.5).
    pub fn theoretical_rg(&self) -> f64 {
        self.b * (self.n_beads as f64).powf(self.nu) / 6.0_f64.sqrt()
    }

    /// Ratio of measured to theoretical Rg.
    pub fn rg_ratio(&self) -> f64 {
        let th = self.theoretical_rg();
        if th < 1e-30 {
            return 0.0;
        }
        self.mean_rg() / th
    }
}

/// Bond breaking kinetics model.
///
/// Models mechanochemical bond breaking using Bell model and
/// Kramers escape theory for a polymer chain under tension.
///
/// References:
/// - Bell, G.I. (1978). *Science* 200, 618.
/// - Wiita, A.P. et al. (2006). *Nature* 440, 598.
#[derive(Debug, Clone)]
pub struct BondBreakingKinetics {
    /// Natural (zero-force) off-rate k_off,0 (s⁻¹).
    pub k_off_0: f64,
    /// Characteristic force for bond breaking F_β = kT / Δx (pN).
    pub f_beta: f64,
    /// Transition state distance Δx (nm).
    pub delta_x: f64,
    /// Thermal energy kT.
    pub kt: f64,
    /// Cumulative number of broken bonds.
    pub n_broken: usize,
}

impl BondBreakingKinetics {
    /// Create a bond breaking kinetics model.
    ///
    /// # Arguments
    /// * `k_off_0`  - zero-force off-rate (1/time)
    /// * `delta_x`  - transition state distance
    /// * `kt`       - thermal energy
    pub fn new(k_off_0: f64, delta_x: f64, kt: f64) -> Self {
        let f_beta = kt / delta_x.max(1e-30);
        Self {
            k_off_0,
            f_beta,
            delta_x,
            kt,
            n_broken: 0,
        }
    }

    /// Bell model: force-dependent off-rate.
    ///
    /// k_off(F) = k_off,0 * exp(F * Δx / kT) = k_off,0 * exp(F / F_β)
    pub fn k_off_bell(&self, force: f64) -> f64 {
        self.k_off_0 * (force / self.f_beta).exp()
    }

    /// Dudko-Hummer-Szabo model for bond breaking.
    ///
    /// k(F) = k₀ * (1 - ν F Δx / ΔG)^(1/ν - 1) * exp(ΔG/kT * \[1 - (1 - ν F Δx/ΔG)^(1/ν)\])
    /// with ν = 1/2 (cusp PES).
    pub fn k_off_dhs(&self, force: f64, delta_g: f64) -> f64 {
        let nu = 0.5f64;
        let ratio = nu * force * self.delta_x / delta_g.max(1e-30);
        if ratio >= 1.0 {
            return f64::INFINITY;
        }
        let base = 1.0 - ratio;
        let power_1 = base.powf(1.0 / nu - 1.0);
        let exponent = delta_g / self.kt * (1.0 - base.powf(1.0 / nu));
        self.k_off_0 * power_1 * exponent.exp()
    }

    /// Bond lifetime at constant force: τ(F) = 1 / k_off(F).
    pub fn lifetime(&self, force: f64) -> f64 {
        1.0 / self.k_off_bell(force).max(1e-30)
    }

    /// Probability that bond survives until time t at force F.
    ///
    /// P(t, F) = exp(-k_off(F) * t)
    pub fn survival_probability(&self, force: f64, t: f64) -> f64 {
        (-self.k_off_bell(force) * t).exp()
    }

    /// Mean rupture force in dynamic force spectroscopy.
    ///
    /// F* = (kT / Δx) * ln(r * Δx / (k_off,0 * kT))
    /// where r is the force loading rate (pN/s).
    pub fn mean_rupture_force(&self, loading_rate: f64) -> f64 {
        let arg = loading_rate * self.delta_x / (self.k_off_0 * self.kt).max(1e-30);
        if arg <= 0.0 {
            return 0.0;
        }
        (self.kt / self.delta_x) * arg.ln()
    }

    /// Estimate number of bonds broken in time dt under force.
    pub fn bonds_broken(&self, n_intact: usize, force: f64, dt: f64) -> usize {
        let rate = self.k_off_bell(force);
        let p_break = (rate * dt).min(1.0);
        (n_intact as f64 * p_break) as usize
    }

    /// Simulate bond lifetime at constant force using deterministic rate.
    pub fn expected_broken_fraction(&self, force: f64, t: f64) -> f64 {
        1.0 - self.survival_probability(force, t)
    }
}

/// Radius of gyration tensor analysis for a polymer chain.
///
/// Computes the full gyration tensor, its principal radii (eigenvalues),
/// and shape anisotropy descriptors: asphericity and prolateness.
#[derive(Debug, Clone)]
pub struct RadiusOfGyration {
    /// Number of beads in the chain.
    pub n_beads: usize,
    /// Bead positions \[x, y, z\].
    pub positions: Vec<[f64; 3]>,
    /// Bead masses (equal-mass by default).
    pub masses: Vec<f64>,
}

impl RadiusOfGyration {
    /// Construct from positions with equal masses.
    pub fn new(positions: Vec<[f64; 3]>) -> Self {
        let n = positions.len();
        Self {
            n_beads: n,
            masses: vec![1.0; n],
            positions,
        }
    }

    /// Construct from positions with individual masses.
    pub fn with_masses(positions: Vec<[f64; 3]>, masses: Vec<f64>) -> Self {
        let n = positions.len();
        Self {
            n_beads: n,
            positions,
            masses,
        }
    }

    /// Centre of mass \[x, y, z\].
    pub fn center_of_mass(&self) -> [f64; 3] {
        let total_mass: f64 = self.masses.iter().sum();
        let mut com = [0.0f64; 3];
        for (pos, &m) in self.positions.iter().zip(&self.masses) {
            com[0] += m * pos[0];
            com[1] += m * pos[1];
            com[2] += m * pos[2];
        }
        let inv = 1.0 / total_mass.max(1e-30);
        [com[0] * inv, com[1] * inv, com[2] * inv]
    }

    /// Gyration tensor S (3x3 symmetric, stored as \[Sxx, Sxy, Sxz, Syy, Syz, Szz\]).
    ///
    /// `S_ab = (1/M) sum m_i (r_i_a - r_cm_a)(r_i_b - r_cm_b)`
    pub fn gyration_tensor(&self) -> [f64; 6] {
        let com = self.center_of_mass();
        let total_mass: f64 = self.masses.iter().sum();
        let mut s = [0.0f64; 6]; // xx, xy, xz, yy, yz, zz
        for (pos, &m) in self.positions.iter().zip(&self.masses) {
            let dx = pos[0] - com[0];
            let dy = pos[1] - com[1];
            let dz = pos[2] - com[2];
            s[0] += m * dx * dx; // Sxx
            s[1] += m * dx * dy; // Sxy
            s[2] += m * dx * dz; // Sxz
            s[3] += m * dy * dy; // Syy
            s[4] += m * dy * dz; // Syz
            s[5] += m * dz * dz; // Szz
        }
        let inv = 1.0 / total_mass.max(1e-30);
        for x in &mut s {
            *x *= inv;
        }
        s
    }

    /// Squared radius of gyration Rg^2 = Tr(S) = l1 + l2 + l3.
    pub fn rg_squared(&self) -> f64 {
        let s = self.gyration_tensor();
        s[0] + s[3] + s[5] // Sxx + Syy + Szz
    }

    /// Radius of gyration Rg = sqrt(Rg^2).
    pub fn rg(&self) -> f64 {
        self.rg_squared().max(0.0).sqrt()
    }

    /// Eigenvalues of the gyration tensor (principal radii squared): \[l1, l2, l3\] sorted ascending.
    ///
    /// Uses Jacobi iteration for the 3x3 symmetric eigenvalue problem.
    pub fn principal_radii_squared(&self) -> [f64; 3] {
        let s = self.gyration_tensor();
        // Pack into 3x3: [Sxx,Sxy,Sxz; Sxy,Syy,Syz; Sxz,Syz,Szz]
        let mut mat = [[s[0], s[1], s[2]], [s[1], s[3], s[4]], [s[2], s[4], s[5]]];
        // Jacobi eigenvalue algorithm (symmetric 3x3)
        let mut vals = [mat[0][0], mat[1][1], mat[2][2]];
        for _ in 0..50 {
            // Find off-diagonal element with largest absolute value
            let (mut p, mut q) = (0, 1);
            let mut max_off = mat[0][1].abs();
            if mat[0][2].abs() > max_off {
                max_off = mat[0][2].abs();
                p = 0;
                q = 2;
            }
            if mat[1][2].abs() > max_off {
                p = 1;
                q = 2;
            }
            let a_pq = mat[p][q];
            if a_pq.abs() < 1e-15 {
                break;
            }
            let diff = mat[q][q] - mat[p][p];
            let theta = if diff.abs() < 1e-14 {
                std::f64::consts::FRAC_PI_4
            } else {
                0.5 * (2.0 * a_pq / diff).atan()
            };
            let c = theta.cos();
            let sn = theta.sin();
            // Apply Jacobi rotation
            let a_pp = mat[p][p];
            let a_qq = mat[q][q];
            mat[p][p] = c * c * a_pp - 2.0 * sn * c * a_pq + sn * sn * a_qq;
            mat[q][q] = sn * sn * a_pp + 2.0 * sn * c * a_pq + c * c * a_qq;
            mat[p][q] = 0.0;
            mat[q][p] = 0.0;
            // Update off-diagonal elements for row/col r != p, q
            let rows: Vec<usize> = (0..3).filter(|&r| r != p && r != q).collect();
            for &r in &rows {
                let a_rp = mat[r][p];
                let a_rq = mat[r][q];
                mat[r][p] = c * a_rp - sn * a_rq;
                mat[p][r] = mat[r][p];
                mat[r][q] = sn * a_rp + c * a_rq;
                mat[q][r] = mat[r][q];
            }
            vals = [mat[0][0], mat[1][1], mat[2][2]];
        }
        // Clamp negative values (numerical noise)
        for v in &mut vals {
            if *v < 0.0 {
                *v = 0.0;
            }
        }
        vals.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        vals
    }

    /// Principal radii \[l1^(1/2), l2^(1/2), l3^(1/2)\] sorted ascending.
    pub fn principal_radii(&self) -> [f64; 3] {
        let lambdas = self.principal_radii_squared();
        [lambdas[0].sqrt(), lambdas[1].sqrt(), lambdas[2].sqrt()]
    }

    /// Asphericity parameter b (Aronovitz-Nelson).
    ///
    /// `b = l3 - (l1 + l2)/2`
    ///
    /// b = 0 for sphere, b > 0 for prolate, b < 0 for oblate.
    pub fn asphericity(&self) -> f64 {
        let lam = self.principal_radii_squared();
        lam[2] - 0.5 * (lam[0] + lam[1])
    }

    /// Acylindricity parameter c.
    ///
    /// `c = l2 - l1`
    ///
    /// c = 0 for cylindrically symmetric shapes.
    pub fn acylindricity(&self) -> f64 {
        let lam = self.principal_radii_squared();
        lam[1] - lam[0]
    }

    /// Relative shape anisotropy k^2 in \[0, 1\].
    ///
    /// `k^2 = (b^2 + (3/4) c^2) / Rg^4`
    ///
    /// k^2 = 0 for sphere, k^2 = 1 for rod.
    pub fn relative_shape_anisotropy(&self) -> f64 {
        let rg2 = self.rg_squared();
        if rg2 < 1e-30 {
            return 0.0;
        }
        let b = self.asphericity();
        let c = self.acylindricity();
        (b * b + 0.75 * c * c) / (rg2 * rg2)
    }

    /// Prolateness S (third invariant of traceless tensor).
    ///
    /// S > 0: prolate (cigar-like), S < 0: oblate (disk-like).
    pub fn prolateness(&self) -> f64 {
        let rg2 = self.rg_squared();
        if rg2 < 1e-30 {
            return 0.0;
        }
        let mean = rg2 / 3.0;
        let lam = self.principal_radii_squared();
        let tl: Vec<f64> = lam.iter().map(|&l| l - mean).collect();
        let num = 27.0 * tl[0] * tl[1] * tl[2];
        let rg6 = rg2 * rg2 * rg2;
        num / rg6.max(1e-60)
    }
}

/// Polymer entanglement analysis using the Z1 primitive path algorithm.
///
/// The Z1 algorithm finds the minimum contour length of the primitive path
/// (the shortest path connecting chain ends while preserving topology).
/// From this, tube diameter and entanglement length are derived.
#[derive(Debug, Clone)]
pub struct PolymerEntanglement {
    /// Number of beads per chain.
    pub n_beads: usize,
    /// Number of chains.
    pub n_chains: usize,
    /// Bond length (Kuhn segment length) b \[m or reduced units\].
    pub bond_length: f64,
    /// Estimated entanglement length Ne (number of monomers between entanglements).
    pub ne: f64,
    /// Tube diameter d_T = b * sqrt(Ne) \[same units as b\].
    pub tube_diameter: f64,
    /// Plateau modulus G_N^0 = rho k_T / Ne \[Pa if consistent units\].
    pub plateau_modulus: f64,
    /// Primitive path positions (simplified: endpoints only).
    pub primitive_path_length: f64,
}

impl PolymerEntanglement {
    /// Construct from chain parameters and estimate Ne via the Z1 algorithm approximation.
    pub fn new(
        n_beads: usize,
        n_chains: usize,
        bond_length: f64,
        density_reduced: f64,
        kt: f64,
    ) -> Self {
        // Mean-field estimate: Ne ~ (b^6/rho)^(1/(3nu-1)) for nu=0.588 -> simplify
        // Use Kavassalis-Noolandi: Ne ~ 13.5 / (rho b^3) in reduced units
        let ne = (13.5 / (density_reduced * bond_length * bond_length * bond_length)).max(2.0);
        let tube_diameter = bond_length * ne.sqrt();
        // Primitive path contour: L_pp ~ N b / sqrt(Ne)
        let primitive_path_length = n_beads as f64 * bond_length / ne.sqrt().max(1e-15);
        // Plateau modulus G_N^0 = rho kT / Ne
        let plateau_modulus = density_reduced * kt / ne.max(1e-15);
        Self {
            n_beads,
            n_chains,
            bond_length,
            ne,
            tube_diameter,
            plateau_modulus,
            primitive_path_length,
        }
    }

    /// Construct with explicit Ne (e.g., from experiments or direct Z1 output).
    pub fn with_ne(
        n_beads: usize,
        n_chains: usize,
        bond_length: f64,
        ne: f64,
        kt: f64,
        density: f64,
    ) -> Self {
        let tube_diameter = bond_length * ne.sqrt();
        let primitive_path_length = n_beads as f64 * bond_length / ne.sqrt().max(1e-15);
        let plateau_modulus = density * kt / ne.max(1e-15);
        Self {
            n_beads,
            n_chains,
            bond_length,
            ne,
            tube_diameter,
            plateau_modulus,
            primitive_path_length,
        }
    }

    /// Number of entanglements per chain: Z = N / Ne.
    pub fn z_entanglements(&self) -> f64 {
        self.n_beads as f64 / self.ne.max(1e-15)
    }

    /// Tube diameter d_T = b * sqrt(Ne).
    pub fn tube_diameter(&self) -> f64 {
        self.tube_diameter
    }

    /// Entanglement length Ne.
    pub fn entanglement_length(&self) -> f64 {
        self.ne
    }

    /// Primitive path contour length estimate L_pp.
    pub fn primitive_path_length(&self) -> f64 {
        self.primitive_path_length
    }

    /// Z1 algorithm simplified: minimise primitive path for a single linear chain.
    pub fn z1_min_path(positions: &[[f64; 3]]) -> f64 {
        if positions.len() < 2 {
            return 0.0;
        }
        let n = positions.len();
        // Min path = straight line from first to last bead
        let dx = positions[n - 1][0] - positions[0][0];
        let dy = positions[n - 1][1] - positions[0][1];
        let dz = positions[n - 1][2] - positions[0][2];
        (dx * dx + dy * dy + dz * dz).sqrt()
    }

    /// Rouse time of entanglement strand tau_e = zeta Ne^2 b^2 / (3 pi^2 kT).
    pub fn rouse_time_entanglement_strand(&self, friction: f64, kt: f64) -> f64 {
        friction * self.ne * self.ne * self.bond_length * self.bond_length
            / (3.0 * PI * PI * kt.max(1e-30))
    }

    /// Reptation (disentanglement) time tau_d ~ 3 Z^3 tau_e (bare theory).
    pub fn reptation_time(&self, friction: f64, kt: f64) -> f64 {
        let z = self.z_entanglements();
        3.0 * z * z * z * self.rouse_time_entanglement_strand(friction, kt)
    }

    /// Viscosity estimate: eta ~ G_N^0 * tau_d.
    pub fn zero_shear_viscosity(&self, friction: f64, kt: f64) -> f64 {
        self.plateau_modulus * self.reptation_time(friction, kt)
    }

    /// Ratio of primitive path length to full chain contour (entanglement reduction factor).
    pub fn path_reduction_ratio(&self) -> f64 {
        let full_contour = (self.n_beads - 1) as f64 * self.bond_length;
        if full_contour < 1e-30 {
            return 1.0;
        }
        self.primitive_path_length / full_contour
    }
}

/// Branched generation-n dendrimer with fractal dimension and hydrodynamic radius.
///
/// Models ideal dendrimer trees: a core with `f` arms (core functionality),
/// each generation adding `b` new branches per terminal group.
#[derive(Debug, Clone)]
pub struct DendrimerModel {
    /// Generation number g (0 = core only).
    pub generation: usize,
    /// Core functionality (number of arms from the core).
    pub core_functionality: usize,
    /// Branch functionality (new branches per terminal group per generation).
    pub branch_functionality: usize,
    /// Spacer length (number of monomers between branch points).
    pub spacer_length: usize,
    /// Monomer size / statistical segment length \[nm\].
    pub segment_length: f64,
    /// Temperature kT \[J or reduced\].
    pub kt: f64,
}

impl DendrimerModel {
    /// Construct a new `DendrimerModel`.
    pub fn new(
        generation: usize,
        core_functionality: usize,
        branch_functionality: usize,
        spacer_length: usize,
        segment_length: f64,
        kt: f64,
    ) -> Self {
        Self {
            generation,
            core_functionality,
            branch_functionality,
            spacer_length,
            segment_length,
            kt,
        }
    }

    /// Construct a PAMAM-like dendrimer (f=3, b=2).
    pub fn pamam(generation: usize, segment_length: f64, kt: f64) -> Self {
        Self::new(generation, 3, 2, 3, segment_length, kt)
    }

    /// Total number of terminal groups (end groups).
    ///
    /// `N_end = f * b^g`
    pub fn n_terminal_groups(&self) -> usize {
        self.core_functionality * self.branch_functionality.pow(self.generation as u32)
    }

    /// Total number of monomers.
    pub fn n_total_monomers(&self) -> usize {
        let mut total = 0usize;
        let ns = self.spacer_length;
        for g in 0..=self.generation {
            let n_spacers = self.core_functionality * self.branch_functionality.pow(g as u32);
            total += n_spacers * ns;
        }
        total
    }

    /// Dendrimer mass (in monomer units).
    pub fn mass(&self) -> f64 {
        self.n_total_monomers() as f64
    }

    /// Radius of gyration estimate (de Gennes/Lescanec).
    pub fn radius_of_gyration(&self) -> f64 {
        let g = self.generation as f64;
        let ns = self.spacer_length as f64;
        self.segment_length * (g * ns + 1.0).sqrt()
    }

    /// Hydrodynamic radius estimate (Stokes-Einstein approximation).
    pub fn hydrodynamic_radius(&self) -> f64 {
        let rg = self.radius_of_gyration();
        let g = self.generation as f64;
        let correction = 1.0 + 1.0 / (g + 1.0);
        rg * correction
    }

    /// Fractal dimension Df of the dendrimer.
    pub fn fractal_dimension(&self) -> f64 {
        let n = self.n_total_monomers() as f64;
        let rg = self.radius_of_gyration();
        let b = self.segment_length;
        if rg < b || n <= 1.0 {
            return 1.0;
        }
        n.ln() / (rg / b).ln().max(1e-10)
    }

    /// Intrinsic viscosity estimate.
    pub fn intrinsic_viscosity(&self) -> f64 {
        let rh = self.hydrodynamic_radius();
        let mass = self.mass().max(1.0);
        10.0 * PI * rh * rh * rh / (3.0 * mass)
    }

    /// Critical generation at which dendrimer becomes space-filling (de Gennes limit).
    pub fn critical_generation(&self) -> usize {
        for g_try in 1..=20usize {
            let n_monomers = self.core_functionality
                * self.branch_functionality.pow(g_try as u32)
                * self.spacer_length;
            let rg = self.segment_length * (g_try as f64 * self.spacer_length as f64).sqrt();
            let n_sphere = (rg / self.segment_length).powi(3) as usize;
            if n_monomers >= n_sphere {
                return g_try;
            }
        }
        20
    }

    /// End-to-end distance of a single arm (spacer chain).
    pub fn arm_end_to_end_distance(&self) -> f64 {
        self.segment_length * (self.spacer_length as f64).sqrt()
    }

    /// Diffusion coefficient via Stokes-Einstein: `D = kT / (6 pi eta R_h)`.
    pub fn diffusion_coefficient(&self, viscosity_solvent: f64) -> f64 {
        let rh = self.hydrodynamic_radius();
        self.kt / (6.0 * PI * viscosity_solvent * rh).max(1e-30)
    }

    /// Check if the dendrimer is in the dense-packing regime (g >= g_crit).
    pub fn is_dense_packed(&self) -> bool {
        self.generation >= self.critical_generation()
    }
}
