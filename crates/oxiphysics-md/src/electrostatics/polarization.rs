// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Polarizable force field models: Drude oscillator, charge equilibration (QEq),
//! and induced dipole solver.

use super::coulomb::COULOMB_K;
use thiserror::Error;

// ---------------------------------------------------------------------------
// Drude oscillator polarizable model
// ---------------------------------------------------------------------------

/// Drude oscillator model for polarizable force fields.
///
/// Each heavy atom carries a Drude particle (charge -q_D) connected via a
/// harmonic spring with force constant k_D. The polarizability is:
///   alpha = q_D^2 / k_D
///
/// References: MacKerell & Roux, J. Comput. Chem. 2003.
#[derive(Debug, Clone)]
pub struct DrudeOscillator {
    /// Spring constant k_D (kJ mol^-1 Å^-2).
    pub k_drude: f64,
    /// Drude charge magnitude q_D (e).
    pub q_drude: f64,
}

impl DrudeOscillator {
    /// Create a new Drude oscillator.
    pub fn new(k_drude: f64, q_drude: f64) -> Self {
        Self { k_drude, q_drude }
    }

    /// Polarizability (Å^3): alpha = q_D^2 / k_D.
    pub fn polarizability(&self) -> f64 {
        if self.k_drude < 1e-30 {
            return 0.0;
        }
        self.q_drude * self.q_drude / self.k_drude
    }

    /// Harmonic spring energy for a displacement `d` (Å) between core and Drude.
    pub fn spring_energy(&self, d: f64) -> f64 {
        0.5 * self.k_drude * d * d
    }

    /// Induced dipole moment (e·Å) given a displacement vector.
    pub fn induced_dipole(&self, displacement: [f64; 3]) -> [f64; 3] {
        [
            self.q_drude * displacement[0],
            self.q_drude * displacement[1],
            self.q_drude * displacement[2],
        ]
    }

    /// Optimal Drude displacement given external electric field E_ext (kJ mol^-1 e^-1 Å^-1).
    ///
    /// Minimising 0.5*k*d^2 - q*E*d gives d = q*E/k.
    pub fn optimal_displacement(&self, e_field: [f64; 3]) -> [f64; 3] {
        if self.k_drude < 1e-30 {
            return [0.0; 3];
        }
        let fac = self.q_drude / self.k_drude;
        [fac * e_field[0], fac * e_field[1], fac * e_field[2]]
    }

    /// Interaction energy between core (charge q_core) and Drude charge
    /// at displacement `d_vec` (Å) from core.
    pub fn core_drude_interaction(&self, q_core: f64, d_vec: [f64; 3]) -> f64 {
        let r2 = d_vec[0] * d_vec[0] + d_vec[1] * d_vec[1] + d_vec[2] * d_vec[2];
        if r2 < 1e-20 {
            return 0.0;
        }
        let r = r2.sqrt();
        COULOMB_K * q_core * (-self.q_drude) / r
    }
}

// ---------------------------------------------------------------------------
// Charge Equilibration (QEq) algorithm
// ---------------------------------------------------------------------------

/// Default atomic shielding parameter γ (Å⁻¹) of the QEq Coulomb kernel.
///
/// The screening length 1/γ = 1.25 Å is of the magnitude of the ReaxFF/EEM
/// shielding parameters.
pub const DEFAULT_QEQ_SHIELDING: f64 = 0.8;

/// Errors from the charge-equilibration solver.
#[derive(Debug, Clone, PartialEq, Error)]
pub enum QeqError {
    /// The per-atom arrays and the coordinate list differ in length.
    #[error(
        "QEq length mismatch: {chi} electronegativities, {eta} hardnesses, \
         {gamma} shielding parameters, {positions} positions"
    )]
    LengthMismatch {
        /// Number of electronegativities `chi`.
        chi: usize,
        /// Number of hardnesses `eta`.
        eta: usize,
        /// Number of shielding parameters `gamma`.
        gamma: usize,
        /// Number of atomic positions.
        positions: usize,
    },
    /// A parameter is not finite, or a shielding parameter is not positive.
    #[error("invalid QEq parameter: {0}")]
    InvalidParameter(String),
    /// The equalisation system has no unique solution (a pivot of the
    /// elimination vanished).
    #[error("QEq equalisation system is singular")]
    Singular,
}

/// Result of a charge equilibration.
#[derive(Debug, Clone, PartialEq)]
pub struct QeqSolution {
    /// Equilibrated atomic charges q_i (e); they sum to the requested total charge.
    pub charges: Vec<f64>,
    /// The common electronegativity μ = χ_i + Σ_j J_ij q_j (kJ mol^-1 e^-1) that
    /// every atom has after equilibration.
    pub electronegativity: f64,
}

/// Charge equilibration (QEq) algorithm for computing atomic partial charges.
///
/// Minimises the electrostatic energy
/// E(q) = Σ_i χ_i q_i + ½ Σ_ij J_ij q_i q_j subject to Σ_i q_i = Q_total.
/// Stationarity gives electronegativity equalisation, χ_i + Σ_j J_ij q_j = μ for
/// every atom, so the charges follow from the (n+1)×(n+1) linear system
///
/// ```text
/// [ J   -1 ] [ q ]   [ -chi ]
/// [ 1ᵀ   0 ] [ μ ] = [ Q    ]
/// ```
///
/// which is solved directly by LU factorisation with partial pivoting.  The
/// solution is a stationary point of E; for non-physical parameters whose
/// hardness matrix is not positive definite on the constraint subspace it is a
/// saddle rather than a minimum.
///
/// Here J\[i\]\[i\] = 2*eta_i (idempotential) and J\[i\]\[j\] is the shielded
/// Coulomb kernel K / (r_ij³ + γ_ij⁻³)^(1/3) with γ_ij = √(γ_i γ_j), which
/// stays finite for r_ij → 0 and tends to K / r_ij at large distance.
///
/// Reference: Rappé & Goddard, J. Phys. Chem. 95, 3358 (1991).
#[derive(Debug, Clone)]
pub struct ChargeEquilibration {
    /// Electronegativity chi_i for each atom (kJ mol^-1 e^-1).
    pub chi: Vec<f64>,
    /// Idempotential eta_i for each atom (kJ mol^-1 e^-2).
    pub eta: Vec<f64>,
    /// Total charge constraint (e).
    pub total_charge: f64,
    /// Shielding parameter γ_i of each atom (Å^-1).
    pub gamma: Vec<f64>,
}

impl ChargeEquilibration {
    /// Create a new QEq solver with the default shielding [`DEFAULT_QEQ_SHIELDING`]
    /// on every atom.
    ///
    /// The array lengths are validated when the system is built or solved.
    pub fn new(chi: Vec<f64>, eta: Vec<f64>, total_charge: f64) -> Self {
        let gamma = vec![DEFAULT_QEQ_SHIELDING; chi.len()];
        Self {
            chi,
            eta,
            total_charge,
            gamma,
        }
    }

    /// Replace the per-atom shielding parameters γ_i (Å^-1).
    pub fn with_shielding(mut self, gamma: Vec<f64>) -> Self {
        self.gamma = gamma;
        self
    }

    /// Check array lengths against `n_positions` and the parameters for validity.
    fn validate(&self, n_positions: usize) -> Result<(), QeqError> {
        let n = self.chi.len();
        if self.eta.len() != n || self.gamma.len() != n || n_positions != n {
            return Err(QeqError::LengthMismatch {
                chi: n,
                eta: self.eta.len(),
                gamma: self.gamma.len(),
                positions: n_positions,
            });
        }
        if !self.total_charge.is_finite() {
            return Err(QeqError::InvalidParameter(format!(
                "total charge must be finite, got {}",
                self.total_charge
            )));
        }
        if self
            .chi
            .iter()
            .chain(self.eta.iter())
            .any(|v| !v.is_finite())
        {
            return Err(QeqError::InvalidParameter(
                "electronegativities and hardnesses must be finite".to_string(),
            ));
        }
        if let Some(bad) = self.gamma.iter().find(|g| !g.is_finite() || **g <= 0.0) {
            return Err(QeqError::InvalidParameter(format!(
                "shielding parameters must be finite and positive, got {bad}"
            )));
        }
        Ok(())
    }

    /// Build the electronegativity equalization matrix J.
    ///
    /// J\[i\]\[i\] = 2*eta_i
    /// J\[i\]\[j\] = COULOMB_K / (r_ij³ + γ_ij⁻³)^(1/3)  for i ≠ j, with
    /// γ_ij = √(γ_i γ_j) (shielded Coulomb kernel, finite at r_ij = 0).
    ///
    /// # Errors
    /// Returns [`QeqError`] if the array lengths differ from `positions.len()` or a
    /// parameter or coordinate is not finite.
    pub fn build_j_matrix(&self, positions: &[[f64; 3]]) -> Result<Vec<Vec<f64>>, QeqError> {
        self.validate(positions.len())?;
        if positions.iter().flatten().any(|x| !x.is_finite()) {
            return Err(QeqError::InvalidParameter(
                "atom positions must be finite".to_string(),
            ));
        }
        let n = positions.len();
        let mut j = vec![vec![0.0_f64; n]; n];
        for i in 0..n {
            j[i][i] = 2.0 * self.eta[i];
            for k in (i + 1)..n {
                let dx = positions[k][0] - positions[i][0];
                let dy = positions[k][1] - positions[i][1];
                let dz = positions[k][2] - positions[i][2];
                let r = (dx * dx + dy * dy + dz * dz).sqrt();
                let gamma_ij = (self.gamma[i] * self.gamma[k]).sqrt();
                let j_ij = COULOMB_K / (r * r * r + gamma_ij.powi(-3)).cbrt();
                j[i][k] = j_ij;
                j[k][i] = j_ij;
            }
        }
        Ok(j)
    }

    /// Solve the equalisation system for the charges and the common
    /// electronegativity.
    ///
    /// # Errors
    /// Returns [`QeqError::LengthMismatch`] / [`QeqError::InvalidParameter`] for
    /// inconsistent or non-finite input and [`QeqError::Singular`] if the system has
    /// no unique solution (for example two atoms whose hardness matrix satisfies
    /// J₁₁ + J₂₂ = 2 J₁₂, or a non-zero total charge on an empty system).
    pub fn solve(&self, positions: &[[f64; 3]]) -> Result<QeqSolution, QeqError> {
        let j = self.build_j_matrix(positions)?;
        let n = positions.len();
        if n == 0 {
            return if self.total_charge == 0.0 {
                Ok(QeqSolution {
                    charges: Vec::new(),
                    electronegativity: 0.0,
                })
            } else {
                Err(QeqError::Singular)
            };
        }
        // Bordered system [[J, -1], [1ᵀ, 0]] [q; μ] = [-χ; Q].
        let mut matrix = vec![vec![0.0_f64; n + 1]; n + 1];
        let mut rhs = vec![0.0_f64; n + 1];
        for (i, row) in j.iter().enumerate() {
            matrix[i][..n].copy_from_slice(row);
            matrix[i][n] = -1.0;
            rhs[i] = -self.chi[i];
        }
        for entry in matrix[n][..n].iter_mut() {
            *entry = 1.0;
        }
        rhs[n] = self.total_charge;

        let mut solution = solve_with_partial_pivoting(matrix, rhs)?;
        let electronegativity = solution.pop().ok_or(QeqError::Singular)?;
        Ok(QeqSolution {
            charges: solution,
            electronegativity,
        })
    }

    /// Compute the QEq charges (e) by solving the equalisation system exactly.
    ///
    /// The charges sum to `total_charge` and every atom ends up with the same
    /// electronegativity χ_i + Σ_j J_ij q_j.
    ///
    /// # Errors
    /// See [`ChargeEquilibration::solve`].
    pub fn compute_charges(&self, positions: &[[f64; 3]]) -> Result<Vec<f64>, QeqError> {
        Ok(self.solve(positions)?.charges)
    }
}

/// Solve the dense system `a x = b` by Gaussian elimination (LU factorisation)
/// with partial pivoting.
///
/// A pivot whose magnitude does not exceed `n · ε · max|a|` marks the matrix as
/// singular; unlike a silent zeroing of such pivots this reports
/// [`QeqError::Singular`].
fn solve_with_partial_pivoting(
    mut a: Vec<Vec<f64>>,
    mut b: Vec<f64>,
) -> Result<Vec<f64>, QeqError> {
    let n = b.len();
    let scale = a
        .iter()
        .flatten()
        .fold(0.0_f64, |largest, value| largest.max(value.abs()));
    let tolerance = n as f64 * f64::EPSILON * scale;

    for col in 0..n {
        let mut pivot_row = col;
        let mut pivot_abs = a[col][col].abs();
        for (row, a_row) in a.iter().enumerate().skip(col + 1) {
            let candidate = a_row[col].abs();
            if candidate > pivot_abs {
                pivot_abs = candidate;
                pivot_row = row;
            }
        }
        if !pivot_abs.is_finite() || pivot_abs <= tolerance {
            return Err(QeqError::Singular);
        }
        a.swap(col, pivot_row);
        b.swap(col, pivot_row);

        let pivot_tail: Vec<f64> = a[col][col..].to_vec();
        let pivot_rhs = b[col];
        for row in (col + 1)..n {
            let factor = a[row][col] / pivot_tail[0];
            if factor == 0.0 {
                continue;
            }
            for (entry, pivot_entry) in a[row][col..].iter_mut().zip(pivot_tail.iter()) {
                *entry -= factor * pivot_entry;
            }
            b[row] -= factor * pivot_rhs;
        }
    }

    // Back substitution on the upper-triangular factor.
    let mut x = vec![0.0_f64; n];
    for row in (0..n).rev() {
        let tail: f64 = a[row][(row + 1)..]
            .iter()
            .zip(x[(row + 1)..].iter())
            .map(|(coefficient, value)| coefficient * value)
            .sum();
        x[row] = (b[row] - tail) / a[row][row];
    }
    if x.iter().any(|v| !v.is_finite()) {
        return Err(QeqError::Singular);
    }
    Ok(x)
}

// ---------------------------------------------------------------------------
// Induced dipole iterative solver
// ---------------------------------------------------------------------------

/// Iterative solver for induced atomic dipoles under a polarizable force field.
///
/// Solves self-consistently: mu_i = alpha_i * (E_ext_i + sum_{j≠i} T_ij * mu_j)
/// where T_ij is the dipole-dipole interaction tensor.
#[derive(Debug, Clone)]
pub struct InducedDipoleSolver {
    /// Atomic polarizabilities alpha_i (Å^3).
    pub polarizabilities: Vec<f64>,
    /// Maximum iterations.
    pub max_iter: usize,
    /// Convergence tolerance (e·Å).
    pub tol: f64,
}

impl InducedDipoleSolver {
    /// Create a new solver.
    pub fn new(polarizabilities: Vec<f64>, max_iter: usize, tol: f64) -> Self {
        Self {
            polarizabilities,
            max_iter,
            tol,
        }
    }

    /// Dipole-dipole interaction tensor T_ij (Å^-3).
    ///
    /// T_ij\[a\]\[b\] = (3 r_a r_b / r^5 - delta_ab / r^3) * COULOMB_K
    fn dipole_tensor(r_vec: [f64; 3]) -> [[f64; 3]; 3] {
        let r2 = r_vec[0] * r_vec[0] + r_vec[1] * r_vec[1] + r_vec[2] * r_vec[2];
        if r2 < 1e-20 {
            return [[0.0; 3]; 3];
        }
        let r = r2.sqrt();
        let r3 = r2 * r;
        let r5 = r3 * r2;
        let mut t = [[0.0_f64; 3]; 3];
        for a in 0..3 {
            for b in 0..3 {
                let delta = if a == b { 1.0 } else { 0.0 };
                t[a][b] = COULOMB_K * (3.0 * r_vec[a] * r_vec[b] / r5 - delta / r3);
            }
        }
        t
    }

    /// Solve for induced dipoles given external electric fields.
    ///
    /// `e_ext[i]` is the external electric field at atom i (kJ mol^-1 e^-1 Å^-1).
    /// Returns induced dipoles mu\[i\] (e·Å).
    pub fn solve(&self, positions: &[[f64; 3]], e_ext: &[[f64; 3]]) -> Vec<[f64; 3]> {
        let n = self.polarizabilities.len();
        assert_eq!(positions.len(), n, "positions length mismatch");
        assert_eq!(e_ext.len(), n, "e_ext length mismatch");

        let mut mu = vec![[0.0_f64; 3]; n];

        for _iter in 0..self.max_iter {
            let mut max_change = 0.0_f64;

            for i in 0..n {
                let alpha_i = self.polarizabilities[i];
                let mut e_total = e_ext[i];

                // Add field from all other induced dipoles
                for j in 0..n {
                    if i == j {
                        continue;
                    }
                    let r_vec = [
                        positions[i][0] - positions[j][0],
                        positions[i][1] - positions[j][1],
                        positions[i][2] - positions[j][2],
                    ];
                    let r2 = r_vec[0] * r_vec[0] + r_vec[1] * r_vec[1] + r_vec[2] * r_vec[2];
                    if r2 < 1.0 {
                        continue;
                    } // avoid near-singularity
                    let t = Self::dipole_tensor(r_vec);
                    for a in 0..3 {
                        for b in 0..3 {
                            e_total[a] += t[a][b] * mu[j][b];
                        }
                    }
                }

                let mu_new = [
                    alpha_i * e_total[0],
                    alpha_i * e_total[1],
                    alpha_i * e_total[2],
                ];
                let change = ((mu_new[0] - mu[i][0]).powi(2)
                    + (mu_new[1] - mu[i][1]).powi(2)
                    + (mu_new[2] - mu[i][2]).powi(2))
                .sqrt();
                if change > max_change {
                    max_change = change;
                }
                mu[i] = mu_new;
            }

            if max_change < self.tol {
                break;
            }
        }
        mu
    }

    /// Polarization energy: E_pol = -0.5 * sum_i mu_i . E_ext_i.
    pub fn polarization_energy(&self, mu: &[[f64; 3]], e_ext: &[[f64; 3]]) -> f64 {
        mu.iter()
            .zip(e_ext.iter())
            .map(|(m, e)| -(m[0] * e[0] + m[1] * e[1] + m[2] * e[2]))
            .sum::<f64>()
            * 0.5
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_drude_polarizability_formula() {
        let d = DrudeOscillator::new(1000.0, 0.1);
        let alpha = d.polarizability();
        let expected = 0.01 / 1000.0;
        assert!(
            (alpha - expected).abs() < 1e-15,
            "alpha = {alpha}, expected {expected}"
        );
    }

    #[test]
    fn test_drude_spring_energy_zero_at_zero() {
        let d = DrudeOscillator::new(500.0, 0.08);
        assert_eq!(d.spring_energy(0.0), 0.0);
    }

    #[test]
    fn test_drude_spring_energy_positive() {
        let d = DrudeOscillator::new(500.0, 0.08);
        let e = d.spring_energy(0.1);
        assert!(e > 0.0, "spring energy positive: {e}");
    }

    #[test]
    fn test_drude_induced_dipole_direction() {
        let d = DrudeOscillator::new(500.0, 0.1);
        let mu = d.induced_dipole([1.0, 0.0, 0.0]);
        assert!(mu[0] > 0.0, "dipole along +x: {}", mu[0]);
        assert_eq!(mu[1], 0.0);
        assert_eq!(mu[2], 0.0);
    }

    #[test]
    fn test_drude_optimal_displacement_direction() {
        let d = DrudeOscillator::new(1000.0, 0.1);
        let disp = d.optimal_displacement([1.0, 0.0, 0.0]);
        assert!(disp[0] > 0.0, "displacement along +x: {}", disp[0]);
    }

    #[test]
    fn test_drude_zero_spring_constant_zero_displacement() {
        let d = DrudeOscillator::new(0.0, 0.1);
        let disp = d.optimal_displacement([5.0, 0.0, 0.0]);
        assert_eq!(disp, [0.0; 3]);
    }

    /// Realistic QEq inputs: χ and η of C/N/O-like atoms in kJ/mol (≈ 5-8 eV).
    fn molecule_qeq() -> (ChargeEquilibration, Vec<[f64; 3]>) {
        let chi = vec![480.0, 640.0, 760.0];
        let eta = vec![470.0, 520.0, 560.0];
        let positions = vec![[0.0, 0.0, 0.0], [1.4, 0.0, 0.0], [1.9, 1.2, 0.3]];
        (ChargeEquilibration::new(chi, eta, 0.0), positions)
    }

    /// μ_i = χ_i + Σ_j J_ij q_j for every atom.
    fn atomic_electronegativities(
        qeq: &ChargeEquilibration,
        positions: &[[f64; 3]],
        charges: &[f64],
    ) -> Vec<f64> {
        let j = qeq.build_j_matrix(positions).expect("valid system");
        (0..charges.len())
            .map(|i| {
                qeq.chi[i]
                    + (0..charges.len())
                        .map(|k| j[i][k] * charges[k])
                        .sum::<f64>()
            })
            .collect()
    }

    #[test]
    fn test_qeq_output_length() {
        let (qeq, positions) = molecule_qeq();
        let charges = qeq.compute_charges(&positions).expect("solvable");
        assert_eq!(charges.len(), 3);
    }

    #[test]
    fn test_qeq_charges_sum_to_total() {
        for total_charge in [0.0, 1.0, -0.75] {
            let (mut qeq, positions) = molecule_qeq();
            qeq.total_charge = total_charge;
            let charges = qeq.compute_charges(&positions).expect("solvable");
            let sum: f64 = charges.iter().sum();
            assert!(
                (sum - total_charge).abs() < 1e-12,
                "charges sum to {sum}, expected {total_charge}"
            );
        }
    }

    #[test]
    fn test_qeq_two_atoms_match_closed_form() {
        // q1 = (chi2 - chi1) / (J11 + J22 - 2 J12), q2 = -q1 for a neutral pair.
        let qeq = ChargeEquilibration::new(vec![450.0, 520.0], vec![480.0, 520.0], 0.0);
        let positions = [[0.0; 3], [2.0, 0.0, 0.0]];
        let j = qeq.build_j_matrix(&positions).expect("valid system");
        let expected_q1 = (520.0 - 450.0) / (j[0][0] + j[1][1] - 2.0 * j[0][1]);
        let charges = qeq.compute_charges(&positions).expect("solvable");
        assert!(
            (charges[0] - expected_q1).abs() < 1e-12,
            "q1 = {} expected {expected_q1}",
            charges[0]
        );
        assert!((charges[1] + expected_q1).abs() < 1e-12);
        // The less electronegative atom gives up electron density.
        assert!(charges[0] > 0.0 && charges[1] < 0.0);
        assert!((expected_q1 - 70.0 / (960.0 + 1040.0 - 2.0 * j[0][1])).abs() < 1e-15);
    }

    #[test]
    fn test_qeq_equalizes_electronegativity_of_all_atoms() {
        let (mut qeq, positions) = molecule_qeq();
        qeq.total_charge = 0.4;
        let solution = qeq.solve(&positions).expect("solvable");
        let mu = atomic_electronegativities(&qeq, &positions, &solution.charges);
        for (i, m) in mu.iter().enumerate() {
            assert!(
                (m - solution.electronegativity).abs() < 1e-8,
                "atom {i}: mu = {m}, common mu = {}",
                solution.electronegativity
            );
        }
        assert!((mu[0] - mu[1]).abs() < 1e-8 && (mu[1] - mu[2]).abs() < 1e-8);
    }

    #[test]
    fn test_qeq_nonphysical_hardness_is_solved_exactly_instead_of_diverging() {
        // eta = 10 kJ/mol with r = 3 A made the old damped iteration blow up (x3.5 per
        // sweep); the direct solve returns the exact stationary point.
        let qeq = ChargeEquilibration::new(vec![3.5, 4.5, 3.0], vec![10.0, 11.0, 9.0], 0.0);
        let positions = [[0.0; 3], [3.0, 0.0, 0.0], [0.0, 3.0, 0.0]];
        let solution = qeq.solve(&positions).expect("solvable");
        assert!(
            solution
                .charges
                .iter()
                .all(|q| q.is_finite() && q.abs() < 1e3)
        );
        let sum: f64 = solution.charges.iter().sum();
        assert!(sum.abs() < 1e-9, "sum {sum}");
        let mu = atomic_electronegativities(&qeq, &positions, &solution.charges);
        let scale = mu.iter().fold(1.0_f64, |m, v| m.max(v.abs()));
        for m in &mu {
            assert!((m - mu[0]).abs() < 1e-9 * scale, "mu = {mu:?}");
        }
    }

    #[test]
    fn test_qeq_identical_symmetric_atoms_stay_neutral() {
        let qeq = ChargeEquilibration::new(vec![500.0; 3], vec![480.0; 3], 0.0);
        // Equilateral triangle: all three atoms are equivalent.
        let h = 3.0_f64.sqrt() / 2.0 * 2.0;
        let positions = [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [1.0, h, 0.0]];
        let charges = qeq.compute_charges(&positions).expect("solvable");
        for q in &charges {
            assert!(q.abs() < 1e-12, "charges {charges:?}");
        }
    }

    #[test]
    fn test_qeq_single_atom_carries_the_total_charge() {
        let qeq = ChargeEquilibration::new(vec![500.0], vec![480.0], -1.0);
        let solution = qeq.solve(&[[0.0; 3]]).expect("solvable");
        assert!((solution.charges[0] + 1.0).abs() < 1e-14);
        // mu = chi + J_ii q = 500 + 960 * (-1)
        assert!((solution.electronegativity - (500.0 - 960.0)).abs() < 1e-10);
    }

    #[test]
    fn test_qeq_empty_system() {
        let neutral = ChargeEquilibration::new(Vec::new(), Vec::new(), 0.0);
        assert_eq!(
            neutral.compute_charges(&[]).expect("empty"),
            Vec::<f64>::new()
        );
        let charged = ChargeEquilibration::new(Vec::new(), Vec::new(), 1.0);
        assert_eq!(charged.compute_charges(&[]), Err(QeqError::Singular));
    }

    #[test]
    fn test_qeq_reports_a_singular_system() {
        // For two atoms the bordered matrix is singular exactly when
        // J11 + J22 = 2 J12; choose the hardness so that J11 = J22 = J12.
        let positions = [[0.0; 3], [2.5, 0.0, 0.0]];
        let probe = ChargeEquilibration::new(vec![0.0; 2], vec![0.0; 2], 0.0);
        let j12 = probe.build_j_matrix(&positions).expect("valid")[0][1];
        let qeq = ChargeEquilibration::new(vec![450.0, 520.0], vec![0.5 * j12; 2], 0.0);
        assert_eq!(qeq.compute_charges(&positions), Err(QeqError::Singular));
    }

    #[test]
    fn test_qeq_rejects_inconsistent_or_invalid_input() {
        let positions = [[0.0; 3], [2.0, 0.0, 0.0]];
        let mismatch = ChargeEquilibration::new(vec![1.0, 2.0], vec![1.0], 0.0);
        assert!(matches!(
            mismatch.compute_charges(&positions),
            Err(QeqError::LengthMismatch { chi: 2, eta: 1, .. })
        ));
        let (qeq, three) = molecule_qeq();
        assert!(matches!(
            qeq.compute_charges(&three[..2]),
            Err(QeqError::LengthMismatch { positions: 2, .. })
        ));
        let nan_chi = ChargeEquilibration::new(vec![f64::NAN, 1.0], vec![1.0, 1.0], 0.0);
        assert!(matches!(
            nan_chi.compute_charges(&positions),
            Err(QeqError::InvalidParameter(_))
        ));
        let bad_gamma = ChargeEquilibration::new(vec![1.0, 2.0], vec![1.0, 1.0], 0.0)
            .with_shielding(vec![0.8, 0.0]);
        assert!(matches!(
            bad_gamma.build_j_matrix(&positions),
            Err(QeqError::InvalidParameter(_))
        ));
        let nan_position = [[f64::NAN, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let valid = ChargeEquilibration::new(vec![1.0, 2.0], vec![1.0, 1.0], 0.0);
        assert!(matches!(
            valid.build_j_matrix(&nan_position),
            Err(QeqError::InvalidParameter(_))
        ));
    }

    #[test]
    fn test_qeq_j_matrix_symmetric() {
        let (qeq, positions) = molecule_qeq();
        let j = qeq.build_j_matrix(&positions).expect("valid system");
        for (a, row) in j.iter().enumerate() {
            for (b, &value) in row.iter().enumerate() {
                assert!((value - j[b][a]).abs() < 1e-12, "J should be symmetric");
            }
        }
    }

    #[test]
    fn test_qeq_diagonal_twice_eta() {
        let chi = vec![3.5, 4.5];
        let eta = vec![10.0, 11.0];
        let qeq = ChargeEquilibration::new(chi, eta.clone(), 0.0);
        let positions = [[0.0; 3], [3.0, 0.0, 0.0]];
        let j = qeq.build_j_matrix(&positions).expect("valid system");
        assert!((j[0][0] - 2.0 * eta[0]).abs() < 1e-12, "J[0][0] = 2*eta[0]");
        assert!((j[1][1] - 2.0 * eta[1]).abs() < 1e-12, "J[1][1] = 2*eta[1]");
    }

    #[test]
    fn test_qeq_shielded_kernel_limits() {
        let qeq = ChargeEquilibration::new(vec![0.0; 2], vec![0.0; 2], 0.0)
            .with_shielding(vec![0.5, 2.0]);
        // gamma_ij = sqrt(0.5 * 2.0) = 1.0 / A.
        let at = |r: f64| {
            qeq.build_j_matrix(&[[0.0; 3], [r, 0.0, 0.0]])
                .expect("valid")[0][1]
        };
        // Finite at r = 0: K * gamma_ij.
        assert!(
            (at(0.0) - COULOMB_K).abs() < 1e-9 * COULOMB_K,
            "{}",
            at(0.0)
        );
        // Closed form at r = 2 A: K / (8 + 1)^(1/3).
        assert!((at(2.0) - COULOMB_K / 9.0_f64.cbrt()).abs() < 1e-9);
        // Bare Coulomb K / r at large distance.
        assert!((at(200.0) * 200.0 / COULOMB_K - 1.0).abs() < 1e-6);
        // Monotonically decreasing and below the bare kernel.
        let mut previous = at(0.0);
        for step in 1..=40 {
            let r = 0.25 * f64::from(step);
            let value = at(r);
            assert!(value < previous && value < COULOMB_K / r, "r = {r}");
            previous = value;
        }
    }

    #[test]
    fn test_qeq_larger_system_satisfies_the_equalisation_conditions() {
        // 64 atoms on a jittered cubic lattice (spacing 2 A): checks the pivoting
        // on a dense, ill-scaled bordered system.
        let n_side = 4;
        let mut positions = Vec::new();
        let mut chi = Vec::new();
        let mut eta = Vec::new();
        for ix in 0..n_side {
            for iy in 0..n_side {
                for iz in 0..n_side {
                    let idx = (ix * n_side + iy) * n_side + iz;
                    let jitter = 0.07 * ((idx * 37 % 11) as f64 - 5.0);
                    positions.push([
                        2.0 * ix as f64 + jitter,
                        2.0 * iy as f64 - 0.5 * jitter,
                        2.0 * iz as f64 + 0.3 * jitter,
                    ]);
                    chi.push(450.0 + 15.0 * ((idx * 13 % 9) as f64));
                    eta.push(430.0 + 10.0 * ((idx * 7 % 5) as f64));
                }
            }
        }
        let qeq = ChargeEquilibration::new(chi, eta, 1.5);
        let solution = qeq.solve(&positions).expect("solvable");
        assert_eq!(solution.charges.len(), n_side * n_side * n_side);
        let sum: f64 = solution.charges.iter().sum();
        assert!((sum - 1.5).abs() < 1e-9, "sum {sum}");
        let mu = atomic_electronegativities(&qeq, &positions, &solution.charges);
        for m in &mu {
            assert!(
                (m - solution.electronegativity).abs() < 1e-7,
                "mu = {m}, common = {}",
                solution.electronegativity
            );
        }
    }

    #[test]
    fn test_induced_dipole_length() {
        let alphas = vec![1.5, 1.5];
        let solver = InducedDipoleSolver::new(alphas, 100, 1e-8);
        let positions = [[0.0; 3], [5.0, 0.0, 0.0]];
        let e_ext = [[1.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
        let mu = solver.solve(&positions, &e_ext);
        assert_eq!(mu.len(), 2);
    }

    #[test]
    fn test_induced_dipole_zero_field_zero_mu() {
        let alphas = vec![1.5];
        let solver = InducedDipoleSolver::new(alphas, 50, 1e-8);
        let positions = [[0.0; 3]];
        let e_ext = [[0.0; 3]];
        let mu = solver.solve(&positions, &e_ext);
        for &v in &mu[0] {
            assert!(v.abs() < 1e-10, "zero field → zero mu: {v}");
        }
    }

    #[test]
    fn test_polarization_energy_negative_for_positive_field() {
        let alphas = vec![1.5];
        let solver = InducedDipoleSolver::new(alphas, 50, 1e-8);
        let positions = [[0.0; 3]];
        let e_ext = [[2.0, 0.0, 0.0]];
        let mu = solver.solve(&positions, &e_ext);
        let e_pol = solver.polarization_energy(&mu, &e_ext);
        assert!(e_pol < 0.0, "polarization energy is negative: {e_pol}");
    }
}
