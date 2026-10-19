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
/// harmonic spring with force constant k_D. In a field E the Drude particle
/// settles at the displacement d = q_D E / k_D, so the induced dipole is
/// mu = q_D d = q_D^2 E / k_D. Expressed as a polarizability volume alpha (Å^3,
/// with mu = alpha E / K and K = [`COULOMB_K`]) this is
///   alpha = K q_D^2 / k_D
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

    /// Polarizability volume (Å^3): alpha = K q_D^2 / k_D with K = [`COULOMB_K`].
    ///
    /// This is the `alpha` of [`InducedDipoleSolver`]: an isolated oscillator in a
    /// field `E` carries the dipole `alpha E / K`, which equals
    /// `q_D * optimal_displacement(E)`.
    pub fn polarizability(&self) -> f64 {
        if self.k_drude < 1e-30 {
            return 0.0;
        }
        COULOMB_K * self.q_drude * self.q_drude / self.k_drude
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

/// Default Thole damping parameter `a` of the exponentially screened
/// dipole-dipole interaction (the AMOEBA value).
pub const DEFAULT_THOLE_DAMPING: f64 = 0.39;

/// Squared separation (Å²) below which two sites count as coincident.
const COINCIDENT_SITES_R2: f64 = 1.0e-24;

/// Errors reported by [`InducedDipoleSolver`].
#[derive(Debug, Clone, PartialEq, Error)]
pub enum PolarizationError {
    /// The position or field arrays differ in length from the polarizability list.
    #[error(
        "polarization length mismatch: {polarizabilities} polarizabilities, \
         {positions} positions, {fields} external fields"
    )]
    LengthMismatch {
        /// Number of atomic polarizabilities.
        polarizabilities: usize,
        /// Number of atomic positions.
        positions: usize,
        /// Number of external field vectors.
        fields: usize,
    },
    /// A polarizability, coordinate, field, tolerance or Thole parameter is
    /// non-finite or out of range, or two undamped sites coincide.
    #[error("invalid polarization input: {0}")]
    InvalidParameter(String),
    /// The self-consistent iteration reached `max_iter` sweeps without meeting
    /// the tolerance, or its residual became non-finite (divergence).
    #[error(
        "induced-dipole iteration did not converge after {iters} sweeps \
         (residual {residual:e} e·Å)"
    )]
    NotConverged {
        /// Number of sweeps performed.
        iters: usize,
        /// Largest change of an induced dipole (e·Å) in the last sweep.
        residual: f64,
    },
}

/// Thole exponential screening factors `(λ3, λ5)` for `x = a·u³`:
/// `λ3 = 1 − e^{−x}` and `λ5 = 1 − (1 + x)·e^{−x}`.
///
/// Both are evaluated without cancellation for small `x`.
fn thole_screening(x: f64) -> (f64, f64) {
    let lambda3 = -(-x).exp_m1();
    let lambda5 = if x < 1.0e-2 {
        // 1 − (1 + x)e^{−x} = x²/2 − x³/3 + x⁴/8 − x⁵/30 + …
        x * x * (0.5 - x / 3.0 + x * x / 8.0 - x * x * x / 30.0)
    } else {
        lambda3 - x * (-x).exp()
    };
    (lambda3, lambda5)
}

/// Iterative solver for induced atomic dipoles under a polarizable force field.
///
/// Solves the self-consistent induction equations
///
/// ```text
/// μ_i = (α_i / K) · ( E_ext,i + Σ_{j≠i} T_ij μ_j ),     K = COULOMB_K
/// ```
///
/// by Gauss-Seidel sweeps.  The atomic polarizabilities `α_i` are *polarizability
/// volumes* in Å³ (so that an isolated atom acquires `μ = α E / K`), the fields
/// `E` are in kJ mol⁻¹ e⁻¹ Å⁻¹ and the dipoles `μ` in e·Å.
///
/// `T_ij` is the dipole field tensor of site `j` at site `i`,
///
/// ```text
/// T_ij[a][b] = K · ( 3 λ5 r_a r_b / r⁵ − λ3 δ_ab / r³ ),
/// λ3 = 1 − e^{−a u³},  λ5 = 1 − (1 + a u³) e^{−a u³},  u = r / (α_i α_j)^{1/6},
/// ```
///
/// with Thole's exponential screening (parameter [`InducedDipoleSolver::thole_a`],
/// default [`DEFAULT_THOLE_DAMPING`]).  The screening keeps the tensor finite as two
/// sites approach each other and prevents the polarization catastrophe at short range;
/// `thole_a = 0` switches it off and leaves the bare point-dipole tensor.
///
/// Reference: Thole, Chem. Phys. 59, 341 (1981).
#[derive(Debug, Clone)]
pub struct InducedDipoleSolver {
    /// Atomic polarizabilities alpha_i (Å^3).
    pub polarizabilities: Vec<f64>,
    /// Maximum number of Gauss-Seidel sweeps.
    pub max_iter: usize,
    /// Convergence tolerance on the largest dipole change per sweep (e·Å).
    pub tol: f64,
    /// Thole damping parameter `a` (dimensionless); `0.0` disables the screening.
    pub thole_a: f64,
}

impl InducedDipoleSolver {
    /// Create a new solver with the default Thole damping
    /// [`DEFAULT_THOLE_DAMPING`].
    pub fn new(polarizabilities: Vec<f64>, max_iter: usize, tol: f64) -> Self {
        Self {
            polarizabilities,
            max_iter,
            tol,
            thole_a: DEFAULT_THOLE_DAMPING,
        }
    }

    /// Replace the Thole damping parameter (`0.0` disables the screening).
    pub fn with_thole_damping(mut self, thole_a: f64) -> Self {
        self.thole_a = thole_a;
        self
    }

    /// Dipole field tensor `T_ij` (kJ mol⁻¹ e⁻¹ Å⁻¹ per e·Å, i.e. K·Å⁻³) of a
    /// dipole at `r_j` acting at `r_i = r_j + r_vec`.
    ///
    /// Returns `None` for coincident sites without Thole damping, where the
    /// point-dipole tensor is singular.
    fn dipole_tensor(
        r_vec: [f64; 3],
        alpha_i: f64,
        alpha_j: f64,
        thole_a: f64,
    ) -> Option<[[f64; 3]; 3]> {
        let r2 = r_vec[0] * r_vec[0] + r_vec[1] * r_vec[1] + r_vec[2] * r_vec[2];
        let damped = thole_a > 0.0 && alpha_i > 0.0 && alpha_j > 0.0;
        if r2 < COINCIDENT_SITES_R2 {
            if !damped {
                return None;
            }
            // r → 0 limit of the screened tensor: the anisotropic part vanishes and
            // λ3 / r³ → a / √(α_i α_j).
            let isotropic = -COULOMB_K * thole_a / (alpha_i * alpha_j).sqrt();
            let mut t = [[0.0_f64; 3]; 3];
            for (a, row) in t.iter_mut().enumerate() {
                row[a] = isotropic;
            }
            return Some(t);
        }
        let r = r2.sqrt();
        let r3 = r2 * r;
        let r5 = r3 * r2;
        let (lambda3, lambda5) = if damped {
            thole_screening(thole_a * r3 / (alpha_i * alpha_j).sqrt())
        } else {
            (1.0, 1.0)
        };
        let mut t = [[0.0_f64; 3]; 3];
        for (a, row) in t.iter_mut().enumerate() {
            for (b, entry) in row.iter_mut().enumerate() {
                let delta = if a == b { 1.0 } else { 0.0 };
                *entry =
                    COULOMB_K * (3.0 * lambda5 * r_vec[a] * r_vec[b] / r5 - lambda3 * delta / r3);
            }
        }
        Some(t)
    }

    /// Check the array lengths and that every input is finite and in range.
    fn validate(
        &self,
        positions: &[[f64; 3]],
        e_ext: &[[f64; 3]],
    ) -> Result<(), PolarizationError> {
        let n = self.polarizabilities.len();
        if positions.len() != n || e_ext.len() != n {
            return Err(PolarizationError::LengthMismatch {
                polarizabilities: n,
                positions: positions.len(),
                fields: e_ext.len(),
            });
        }
        if let Some(bad) = self
            .polarizabilities
            .iter()
            .find(|a| !a.is_finite() || **a < 0.0)
        {
            return Err(PolarizationError::InvalidParameter(format!(
                "polarizabilities must be finite and non-negative, got {bad}"
            )));
        }
        if !self.tol.is_finite() || self.tol <= 0.0 {
            return Err(PolarizationError::InvalidParameter(format!(
                "tolerance must be finite and positive, got {}",
                self.tol
            )));
        }
        if !self.thole_a.is_finite() || self.thole_a < 0.0 {
            return Err(PolarizationError::InvalidParameter(format!(
                "Thole damping parameter must be finite and non-negative, got {}",
                self.thole_a
            )));
        }
        if positions
            .iter()
            .chain(e_ext.iter())
            .flatten()
            .any(|x| !x.is_finite())
        {
            return Err(PolarizationError::InvalidParameter(
                "positions and external fields must be finite".to_string(),
            ));
        }
        Ok(())
    }

    /// Solve for the induced dipoles given the external electric fields.
    ///
    /// `e_ext[i]` is the external electric field at atom `i`
    /// (kJ mol^-1 e^-1 Å^-1).  Returns the induced dipoles `mu[i]` (e·Å); a site
    /// with `α_i = 0` carries no dipole and radiates none.
    ///
    /// # Errors
    /// * [`PolarizationError::LengthMismatch`] / [`PolarizationError::InvalidParameter`]
    ///   for inconsistent or non-finite input (including coincident sites when the
    ///   Thole damping is off).
    /// * [`PolarizationError::NotConverged`] if the largest dipole change of a sweep
    ///   does not fall below `tol` within `max_iter` sweeps, or becomes non-finite.
    pub fn solve(
        &self,
        positions: &[[f64; 3]],
        e_ext: &[[f64; 3]],
    ) -> Result<Vec<[f64; 3]>, PolarizationError> {
        self.validate(positions, e_ext)?;
        let n = self.polarizabilities.len();
        let mut mu = vec![[0.0_f64; 3]; n];
        if n == 0 {
            return Ok(mu);
        }

        let mut residual = f64::INFINITY;
        for sweep in 1..=self.max_iter {
            let mut max_change = 0.0_f64;

            for i in 0..n {
                let alpha_i = self.polarizabilities[i];
                if alpha_i == 0.0 {
                    continue;
                }
                let mut e_total = e_ext[i];

                // Field at site i from the current induced dipoles of all other sites.
                for j in 0..n {
                    let alpha_j = self.polarizabilities[j];
                    if i == j || alpha_j == 0.0 {
                        continue;
                    }
                    let r_vec = [
                        positions[i][0] - positions[j][0],
                        positions[i][1] - positions[j][1],
                        positions[i][2] - positions[j][2],
                    ];
                    let t = Self::dipole_tensor(r_vec, alpha_i, alpha_j, self.thole_a).ok_or_else(
                        || {
                            PolarizationError::InvalidParameter(format!(
                                "polarizable sites {i} and {j} coincide and Thole damping is off"
                            ))
                        },
                    )?;
                    for (a, e_a) in e_total.iter_mut().enumerate() {
                        for (b, mu_b) in mu[j].iter().enumerate() {
                            *e_a += t[a][b] * mu_b;
                        }
                    }
                }

                // mu = alpha E / K  (alpha in Å^3, E in kJ mol^-1 e^-1 Å^-1).
                let scale = alpha_i / COULOMB_K;
                let mu_new = [scale * e_total[0], scale * e_total[1], scale * e_total[2]];
                // hypot: no spurious overflow of the squares for a diverging iteration.
                let change = (mu_new[0] - mu[i][0])
                    .hypot(mu_new[1] - mu[i][1])
                    .hypot(mu_new[2] - mu[i][2]);
                if change.is_nan() || change > max_change {
                    max_change = change;
                }
                mu[i] = mu_new;
            }

            residual = max_change;
            if !residual.is_finite() {
                return Err(PolarizationError::NotConverged {
                    iters: sweep,
                    residual,
                });
            }
            if residual < self.tol {
                return Ok(mu);
            }
        }
        Err(PolarizationError::NotConverged {
            iters: self.max_iter,
            residual,
        })
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
        // alpha = K q^2 / k_D in Å^3 (the bare q^2 / k_D is in e^2 Å^2 mol kJ^-1).
        let expected = COULOMB_K * 0.01 / 1000.0;
        assert!(
            (alpha - expected).abs() < 1e-12,
            "alpha = {alpha}, expected {expected}"
        );
    }

    #[test]
    fn test_drude_polarizability_is_a_realistic_volume() {
        // CHARMM/SWM4-NDP Drude oxygen: k_D = 1000 kcal mol^-1 Å^-2, q_D = -1.71636 e,
        // alpha_O = 0.978 Å^3.
        let d = DrudeOscillator::new(4184.0, -1.71636);
        let alpha = d.polarizability();
        assert!((alpha - 0.978).abs() < 0.005, "alpha = {alpha} Å^3");
    }

    #[test]
    fn test_drude_polarizability_matches_the_induced_dipole() {
        // The dipole of the relaxed oscillator, mu = q_D d, equals alpha E / K.
        let d = DrudeOscillator::new(4184.0, 1.1);
        let field = [0.3, -0.2, 0.5];
        let mu = d.induced_dipole(d.optimal_displacement(field));
        let alpha = d.polarizability();
        for axis in 0..3 {
            let expected = alpha * field[axis] / COULOMB_K;
            assert!(
                (mu[axis] - expected).abs() < 1e-15,
                "axis {axis}: mu = {}, alpha E / K = {expected}",
                mu[axis]
            );
        }
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

    /// Reference Thole-screened dipole tensor written directly from the formulas
    /// (independent of the solver's implementation).
    fn thole_tensor_reference(
        r_vec: [f64; 3],
        alpha_i: f64,
        alpha_j: f64,
        thole_a: f64,
    ) -> [[f64; 3]; 3] {
        let r = (r_vec[0].powi(2) + r_vec[1].powi(2) + r_vec[2].powi(2)).sqrt();
        let u = r / (alpha_i * alpha_j).powf(1.0 / 6.0);
        let x = thole_a * u.powi(3);
        let lambda3 = 1.0 - (-x).exp();
        let lambda5 = 1.0 - (1.0 + x) * (-x).exp();
        let mut t = [[0.0; 3]; 3];
        for a in 0..3 {
            for b in 0..3 {
                let delta = if a == b { 1.0 } else { 0.0 };
                t[a][b] = COULOMB_K
                    * (3.0 * lambda5 * r_vec[a] * r_vec[b] / r.powi(5)
                        - lambda3 * delta / r.powi(3));
            }
        }
        t
    }

    #[test]
    fn test_induced_dipole_two_sites_along_the_field_match_the_closed_form() {
        // Two alpha = 1.5 Å^3 sites 5 Å apart along E = 1: each dipole feels the field
        // of the other, mu = alpha E / (K (1 - 2 alpha / r^3)) (the old solver diverged
        // to ~1e104 on exactly this input because it lacked the 1/K).
        let (alpha, r) = (1.5, 5.0);
        let positions = [[0.0; 3], [r, 0.0, 0.0]];
        let e_ext = [[1.0, 0.0, 0.0]; 2];
        let expected = alpha / (COULOMB_K * (1.0 - 2.0 * alpha / (r * r * r)));
        assert!(
            (expected - 1.1062e-3).abs() < 1e-6,
            "closed form {expected}"
        );
        for thole in [0.0, DEFAULT_THOLE_DAMPING] {
            let solver =
                InducedDipoleSolver::new(vec![alpha; 2], 200, 1e-14).with_thole_damping(thole);
            let mu = solver.solve(&positions, &e_ext).expect("converges");
            assert_eq!(mu.len(), 2);
            for dipole in &mu {
                assert!(
                    (dipole[0] - expected).abs() <= 1e-12 * expected,
                    "a = {thole}: mu_x = {}, expected {expected}",
                    dipole[0]
                );
                assert!(dipole[1].abs() < 1e-18 && dipole[2].abs() < 1e-18);
            }
        }
    }

    #[test]
    fn test_induced_dipole_two_sites_across_the_field_match_the_closed_form() {
        // Field perpendicular to the axis: T_yy = -K / r^3, so
        // mu = alpha E / (K (1 + alpha / r^3)).
        let (alpha, r) = (1.5, 5.0);
        let positions = [[0.0; 3], [r, 0.0, 0.0]];
        let e_ext = [[0.0, 1.0, 0.0]; 2];
        let expected = alpha / (COULOMB_K * (1.0 + alpha / (r * r * r)));
        let solver = InducedDipoleSolver::new(vec![alpha; 2], 200, 1e-14).with_thole_damping(0.0);
        let mu = solver.solve(&positions, &e_ext).expect("converges");
        for dipole in &mu {
            assert!(
                (dipole[1] - expected).abs() <= 1e-12 * expected,
                "mu_y = {}, expected {expected}",
                dipole[1]
            );
            assert!(dipole[0].abs() < 1e-18 && dipole[2].abs() < 1e-18);
        }
    }

    #[test]
    fn test_induced_dipole_isolated_site_is_alpha_e_over_k() {
        let solver = InducedDipoleSolver::new(vec![2.3], 10, 1e-12);
        let mu = solver
            .solve(&[[1.0, 2.0, 3.0]], &[[0.4, -0.2, 0.9]])
            .expect("converges");
        for (axis, field) in [0.4, -0.2, 0.9].into_iter().enumerate() {
            assert!((mu[0][axis] - 2.3 * field / COULOMB_K).abs() < 1e-15);
        }
    }

    #[test]
    fn test_induced_dipole_solution_satisfies_the_induction_equations() {
        // Four sites with different polarizabilities in a general position: the solution
        // must be a fixed point mu_i = alpha_i / K (E_i + sum_j T_ij mu_j) of the
        // reference (independently coded) Thole tensor.
        let alphas = vec![1.4, 0.9, 2.1, 1.7];
        let positions = [
            [0.0, 0.0, 0.0],
            [2.6, 0.4, -0.3],
            [0.7, 3.1, 0.8],
            [-1.9, 1.2, 2.2],
        ];
        let e_ext = [
            [0.8, -0.1, 0.3],
            [0.2, 0.5, -0.4],
            [-0.6, 0.3, 0.1],
            [0.1, 0.1, 0.7],
        ];
        let solver = InducedDipoleSolver::new(alphas.clone(), 500, 1e-13);
        let mu = solver.solve(&positions, &e_ext).expect("converges");
        for i in 0..4 {
            let mut field = e_ext[i];
            for j in 0..4 {
                if i == j {
                    continue;
                }
                let r_vec = [
                    positions[i][0] - positions[j][0],
                    positions[i][1] - positions[j][1],
                    positions[i][2] - positions[j][2],
                ];
                let t = thole_tensor_reference(r_vec, alphas[i], alphas[j], DEFAULT_THOLE_DAMPING);
                for a in 0..3 {
                    for b in 0..3 {
                        field[a] += t[a][b] * mu[j][b];
                    }
                }
            }
            for axis in 0..3 {
                let rhs = alphas[i] * field[axis] / COULOMB_K;
                assert!(
                    (mu[i][axis] - rhs).abs() < 1e-11 * rhs.abs().max(1e-4),
                    "site {i} axis {axis}: mu = {}, alpha (E + T mu) / K = {rhs}",
                    mu[i][axis]
                );
            }
        }
    }

    #[test]
    fn test_induced_dipole_close_pairs_are_screened_not_skipped() {
        // Two alpha = 1.5 sites only 1 Å apart (the old solver silently dropped every pair
        // with r < 1 Å). The bare tensor would give a catastrophic coupling there, the Thole
        // screening keeps it finite and the pair clearly feels each other.
        let (alpha, r) = (1.5, 0.9);
        let positions = [[0.0; 3], [r, 0.0, 0.0]];
        let e_ext = [[1.0, 0.0, 0.0]; 2];
        let solver = InducedDipoleSolver::new(vec![alpha; 2], 500, 1e-14);
        let mu = solver.solve(&positions, &e_ext).expect("converges");
        let t = thole_tensor_reference([r, 0.0, 0.0], alpha, alpha, DEFAULT_THOLE_DAMPING);
        let expected = alpha / COULOMB_K / (1.0 - alpha * t[0][0] / COULOMB_K);
        for dipole in &mu {
            assert!(
                (dipole[0] - expected).abs() < 1e-10 * expected,
                "mu_x = {}, expected {expected}",
                dipole[0]
            );
        }
        let isolated = alpha / COULOMB_K;
        assert!(
            (mu[0][0] - isolated).abs() > 0.05 * isolated,
            "a pair at 0.9 Å must be coupled: {} vs isolated {isolated}",
            mu[0][0]
        );
        // The bare point-dipole tensor is singular for coincident sites.
        let bare = InducedDipoleSolver::new(vec![alpha; 2], 50, 1e-12).with_thole_damping(0.0);
        let coincident = [[0.5, 0.5, 0.5]; 2];
        assert!(matches!(
            bare.solve(&coincident, &e_ext),
            Err(PolarizationError::InvalidParameter(_))
        ));
        // With Thole damping coincident sites are finite and well defined.
        let damped = InducedDipoleSolver::new(vec![alpha; 2], 500, 1e-14);
        let mu = damped.solve(&coincident, &e_ext).expect("finite limit");
        let limit = -COULOMB_K * DEFAULT_THOLE_DAMPING / alpha;
        let expected_coincident = alpha / COULOMB_K / (1.0 - alpha * limit / COULOMB_K);
        // The sweep stops once the dipole changes by less than 1e-14 e·Å; the contraction
        // factor per sweep is the Thole parameter a = 0.39.
        assert!((mu[0][0] - expected_coincident).abs() < 1e-9 * expected_coincident);
    }

    #[test]
    fn test_induced_dipole_reports_the_polarization_catastrophe() {
        // alpha = 100 Å^3 at r = 3 Å without screening: 2 alpha / r^3 = 7.4 > 1, the
        // iteration diverges and must say so instead of returning garbage.
        let positions = [[0.0; 3], [3.0, 0.0, 0.0]];
        let e_ext = [[1.0, 0.0, 0.0]; 2];
        let solver = InducedDipoleSolver::new(vec![100.0; 2], 100, 1e-10).with_thole_damping(0.0);
        match solver.solve(&positions, &e_ext) {
            Err(PolarizationError::NotConverged { iters, residual }) => {
                assert_eq!(iters, 100);
                assert!(residual > 1e100, "residual {residual}");
            }
            other => panic!("expected NotConverged, got {other:?}"),
        }
        // A long run overflows to infinity and stops as soon as the residual is non-finite.
        let overflowing =
            InducedDipoleSolver::new(vec![100.0; 2], 5000, 1e-10).with_thole_damping(0.0);
        match overflowing.solve(&positions, &e_ext) {
            Err(PolarizationError::NotConverged { iters, residual }) => {
                assert!(iters < 5000, "stopped after {iters} sweeps");
                assert!(!residual.is_finite());
            }
            other => panic!("expected NotConverged, got {other:?}"),
        }
    }

    #[test]
    fn test_induced_dipole_reports_an_exhausted_iteration_budget() {
        let positions = [[0.0; 3], [5.0, 0.0, 0.0]];
        let e_ext = [[1.0, 0.0, 0.0]; 2];
        let solver = InducedDipoleSolver::new(vec![1.5; 2], 1, 1e-14);
        match solver.solve(&positions, &e_ext) {
            Err(PolarizationError::NotConverged { iters, residual }) => {
                assert_eq!(iters, 1);
                assert!(residual > 1e-14 && residual.is_finite());
            }
            other => panic!("expected NotConverged, got {other:?}"),
        }
        let no_sweeps = InducedDipoleSolver::new(vec![1.5; 2], 0, 1e-14);
        assert!(matches!(
            no_sweeps.solve(&positions, &e_ext),
            Err(PolarizationError::NotConverged { iters: 0, .. })
        ));
    }

    #[test]
    fn test_induced_dipole_rejects_inconsistent_or_invalid_input() {
        let solver = InducedDipoleSolver::new(vec![1.5, 1.5], 50, 1e-10);
        let positions = [[0.0; 3], [5.0, 0.0, 0.0]];
        let e_ext = [[1.0, 0.0, 0.0]; 2];
        assert_eq!(
            solver.solve(&positions[..1], &e_ext),
            Err(PolarizationError::LengthMismatch {
                polarizabilities: 2,
                positions: 1,
                fields: 2
            })
        );
        assert!(matches!(
            solver.solve(&positions, &e_ext[..1]),
            Err(PolarizationError::LengthMismatch { fields: 1, .. })
        ));
        let negative = InducedDipoleSolver::new(vec![1.5, -1.0], 50, 1e-10);
        assert!(matches!(
            negative.solve(&positions, &e_ext),
            Err(PolarizationError::InvalidParameter(_))
        ));
        let nan_field = [[f64::NAN, 0.0, 0.0], [1.0, 0.0, 0.0]];
        assert!(matches!(
            solver.solve(&positions, &nan_field),
            Err(PolarizationError::InvalidParameter(_))
        ));
        let zero_tol = InducedDipoleSolver::new(vec![1.5, 1.5], 50, 0.0);
        assert!(matches!(
            zero_tol.solve(&positions, &e_ext),
            Err(PolarizationError::InvalidParameter(_))
        ));
        let bad_thole =
            InducedDipoleSolver::new(vec![1.5, 1.5], 50, 1e-10).with_thole_damping(-0.1);
        assert!(matches!(
            bad_thole.solve(&positions, &e_ext),
            Err(PolarizationError::InvalidParameter(_))
        ));
    }

    #[test]
    fn test_induced_dipole_non_polarizable_site_neither_responds_nor_radiates() {
        let positions = [[0.0; 3], [2.0, 1.0, 0.0], [5.0, 0.0, 0.0]];
        let e_ext = [[1.0, 0.0, 0.0]; 3];
        let with_dummy = InducedDipoleSolver::new(vec![1.5, 0.0, 1.5], 200, 1e-14)
            .solve(&positions, &e_ext)
            .expect("converges");
        assert_eq!(with_dummy[1], [0.0; 3]);
        let pair = InducedDipoleSolver::new(vec![1.5, 1.5], 200, 1e-14)
            .solve(&[positions[0], positions[2]], &[e_ext[0], e_ext[2]])
            .expect("converges");
        for axis in 0..3 {
            assert!((with_dummy[0][axis] - pair[0][axis]).abs() < 1e-16);
            assert!((with_dummy[2][axis] - pair[1][axis]).abs() < 1e-16);
        }
    }

    #[test]
    fn test_induced_dipole_empty_system_and_zero_field() {
        let empty = InducedDipoleSolver::new(Vec::new(), 10, 1e-8);
        assert_eq!(empty.solve(&[], &[]), Ok(Vec::new()));
        let solver = InducedDipoleSolver::new(vec![1.5], 50, 1e-8);
        let mu = solver.solve(&[[0.0; 3]], &[[0.0; 3]]).expect("converges");
        for &v in &mu[0] {
            assert!(v.abs() < 1e-10, "zero field → zero mu: {v}");
        }
    }

    #[test]
    fn test_polarization_energy_is_half_the_field_dipole_product() {
        let (alpha, r) = (1.5, 5.0);
        let positions = [[0.0; 3], [r, 0.0, 0.0]];
        let e_ext = [[1.0, 0.0, 0.0]; 2];
        let solver = InducedDipoleSolver::new(vec![alpha; 2], 200, 1e-14).with_thole_damping(0.0);
        let mu = solver.solve(&positions, &e_ext).expect("converges");
        let e_pol = solver.polarization_energy(&mu, &e_ext);
        let expected = -alpha / (COULOMB_K * (1.0 - 2.0 * alpha / r.powi(3)));
        assert!(e_pol < 0.0, "polarization energy is negative: {e_pol}");
        assert!((e_pol - expected).abs() < 1e-12 * expected.abs());
    }

    #[test]
    fn test_thole_screened_tensor_matches_the_reference_and_has_a_finite_limit() {
        let (alpha_i, alpha_j, a) = (1.2, 2.4, DEFAULT_THOLE_DAMPING);
        for r_vec in [[2.0, 0.5, -0.3], [0.6, 0.2, 0.1], [4.0, -3.0, 2.0]] {
            let t = InducedDipoleSolver::dipole_tensor(r_vec, alpha_i, alpha_j, a)
                .expect("damped tensor exists");
            let reference = thole_tensor_reference(r_vec, alpha_i, alpha_j, a);
            for row in 0..3 {
                for col in 0..3 {
                    assert!(
                        (t[row][col] - reference[row][col]).abs()
                            <= 1e-9 * reference[row][col].abs().max(1.0),
                        "T[{row}][{col}] = {}, reference {}",
                        t[row][col],
                        reference[row][col]
                    );
                }
            }
        }
        // r -> 0 limit: -K a / sqrt(alpha_i alpha_j) on the diagonal.
        let limit = -COULOMB_K * a / (alpha_i * alpha_j).sqrt();
        let at_zero = InducedDipoleSolver::dipole_tensor([0.0; 3], alpha_i, alpha_j, a)
            .expect("limit exists");
        let near_zero = InducedDipoleSolver::dipole_tensor([1.0e-5, 0.0, 0.0], alpha_i, alpha_j, a)
            .expect("finite");
        for axis in 0..3 {
            assert!((at_zero[axis][axis] - limit).abs() < 1e-12 * limit.abs());
            assert!(
                (near_zero[axis][axis] - limit).abs() < 1e-6 * limit.abs(),
                "T[{axis}][{axis}] at 1e-5 Å: {} vs limit {limit}",
                near_zero[axis][axis]
            );
        }
        // Without screening the point-dipole tensor is singular at r = 0.
        assert!(InducedDipoleSolver::dipole_tensor([0.0; 3], alpha_i, alpha_j, 0.0).is_none());
    }

    #[test]
    fn test_thole_screening_series_matches_the_direct_formula() {
        // For moderate and large x the textbook expressions are accurate in f64.
        for x in [9.9e-3, 1.0e-2, 0.2, 3.0, 40.0] {
            let (lambda3, lambda5) = thole_screening(x);
            let direct3 = 1.0 - (-x).exp();
            let direct5 = 1.0 - (1.0 + x) * (-x).exp();
            assert!((lambda3 - direct3).abs() <= 1e-9 * direct3, "lambda3({x})");
            assert!(
                (lambda5 - direct5).abs() <= 1e-9 * direct5,
                "lambda5({x}) = {lambda5}, direct {direct5}"
            );
            assert!(lambda5 < lambda3 && lambda3 < 1.0 + 1e-15);
        }
        // For tiny x the direct difference loses all digits; compare with the Taylor
        // series x^2/2 - x^3/3 + x^4/8 - x^5/30 + x^6/144 instead.
        for x in [1.0e-8, 1.0e-6, 1.0e-4, 1.0e-3] {
            let (lambda3, lambda5) = thole_screening(x);
            let taylor3 = x - x * x / 2.0 + x.powi(3) / 6.0 - x.powi(4) / 24.0;
            let taylor5 = x * x / 2.0 - x.powi(3) / 3.0 + x.powi(4) / 8.0 - x.powi(5) / 30.0
                + x.powi(6) / 144.0;
            assert!((lambda3 - taylor3).abs() <= 1e-12 * taylor3, "lambda3({x})");
            assert!((lambda5 - taylor5).abs() <= 1e-9 * taylor5, "lambda5({x})");
        }
        // Continuous across the switch from the series to the direct form.
        let (_, below) = thole_screening(1.0e-2 - 1e-12);
        let (_, above) = thole_screening(1.0e-2 + 1e-12);
        assert!((above - below).abs() < 1e-12);
    }
}
