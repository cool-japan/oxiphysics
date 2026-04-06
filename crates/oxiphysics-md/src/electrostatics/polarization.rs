// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Polarizable force field models: Drude oscillator, charge equilibration (QEq),
//! and induced dipole solver.

use super::coulomb::COULOMB_K;

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
#[allow(dead_code)]
pub struct DrudeOscillator {
    /// Spring constant k_D (kJ mol^-1 Å^-2).
    pub k_drude: f64,
    /// Drude charge magnitude q_D (e).
    pub q_drude: f64,
}

#[allow(dead_code)]
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

/// Charge equilibration (QEq) algorithm for computing atomic partial charges.
///
/// Solves the linear system: J * q = -chi, subject to sum(q) = Q_total.
///
/// Where J\[i\]\[i\] = 2*eta_i (idempotential), J\[i\]\[j\] = (Coulomb kernel)(r_ij).
///
/// Reference: Rappé & Goddard, J. Phys. Chem. 95, 3358 (1991).
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct ChargeEquilibration {
    /// Electronegativity chi_i for each atom (kJ mol^-1 e^-1).
    pub chi: Vec<f64>,
    /// Idempotential eta_i for each atom (kJ mol^-1 e^-2).
    pub eta: Vec<f64>,
    /// Total charge constraint (e).
    pub total_charge: f64,
}

#[allow(dead_code)]
impl ChargeEquilibration {
    /// Create a new QEq solver.
    pub fn new(chi: Vec<f64>, eta: Vec<f64>, total_charge: f64) -> Self {
        assert_eq!(chi.len(), eta.len(), "chi and eta must have equal length");
        Self {
            chi,
            eta,
            total_charge,
        }
    }

    /// Build the electronegativity equalization matrix J.
    ///
    /// J\[i\]\[i\] = 2*eta_i
    /// J\[i\]\[j\] = COULOMB_K / r_ij  for i ≠ j (screened with gamma=1)
    pub fn build_j_matrix(&self, positions: &[[f64; 3]]) -> Vec<Vec<f64>> {
        let n = positions.len();
        assert_eq!(
            n,
            self.chi.len(),
            "positions and chi must have equal length"
        );
        let mut j = vec![vec![0.0_f64; n]; n];
        #[allow(clippy::needless_range_loop)]
        for i in 0..n {
            j[i][i] = 2.0 * self.eta[i];
            for k in (i + 1)..n {
                let dx = positions[k][0] - positions[i][0];
                let dy = positions[k][1] - positions[i][1];
                let dz = positions[k][2] - positions[i][2];
                let r = (dx * dx + dy * dy + dz * dz).sqrt().max(0.5); // clamp to avoid singularity
                let j_ij = COULOMB_K / r;
                j[i][k] = j_ij;
                j[k][i] = j_ij;
            }
        }
        j
    }

    /// Compute QEq charges via Gauss-Seidel iteration.
    ///
    /// Returns atomic charges after `max_iter` iterations.
    pub fn compute_charges(&self, positions: &[[f64; 3]], max_iter: usize, tol: f64) -> Vec<f64> {
        let n = self.chi.len();
        assert_eq!(positions.len(), n, "positions length mismatch");
        let j = self.build_j_matrix(positions);

        // Initial guess: distribute total charge equally
        let mut q = vec![self.total_charge / n as f64; n];

        // Damping factor to avoid divergence in Gauss-Seidel
        let damping = 0.05_f64;

        for _iter in 0..max_iter {
            let mut max_change = 0.0_f64;
            // Compute equalization potential at each atom: mu_i = chi_i + sum_j J[i][j]*q[j]
            let mu: Vec<f64> = (0..n)
                .map(|i| self.chi[i] + (0..n).map(|j_idx| j[i][j_idx] * q[j_idx]).sum::<f64>())
                .collect();
            let mu_avg = mu.iter().sum::<f64>() / n as f64;

            // Gradient: dE/dq_i = 2*(mu_i - mu_avg)
            // Damped update: q_i -= damping * dE/dq_i
            for i in 0..n {
                let gradient = 2.0 * (mu[i] - mu_avg);
                let dq = -damping * gradient;
                if dq.abs() > max_change {
                    max_change = dq.abs();
                }
                q[i] += dq;
            }

            // Rescale to enforce total charge constraint
            let q_sum: f64 = q.iter().sum();
            let correction = (self.total_charge - q_sum) / n as f64;
            for qi in q.iter_mut() {
                *qi += correction;
            }

            if max_change < tol {
                break;
            }
        }
        q
    }
}

// ---------------------------------------------------------------------------
// Induced dipole iterative solver
// ---------------------------------------------------------------------------

/// Iterative solver for induced atomic dipoles under a polarizable force field.
///
/// Solves self-consistently: mu_i = alpha_i * (E_ext_i + sum_{j≠i} T_ij * mu_j)
/// where T_ij is the dipole-dipole interaction tensor.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct InducedDipoleSolver {
    /// Atomic polarizabilities alpha_i (Å^3).
    pub polarizabilities: Vec<f64>,
    /// Maximum iterations.
    pub max_iter: usize,
    /// Convergence tolerance (e·Å).
    pub tol: f64,
}

#[allow(dead_code)]
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

    #[test]
    fn test_qeq_output_length() {
        let chi = vec![3.5, 4.5];
        let eta = vec![10.0, 11.0];
        let qeq = ChargeEquilibration::new(chi, eta, 0.0);
        let positions = [[0.0; 3], [3.0, 0.0, 0.0]];
        let charges = qeq.compute_charges(&positions, 100, 1e-6);
        assert_eq!(charges.len(), 2);
    }

    #[test]
    fn test_qeq_charges_sum_to_total() {
        let chi = vec![3.5, 4.5, 3.0];
        let eta = vec![10.0, 11.0, 9.0];
        let total_charge = 0.0;
        let qeq = ChargeEquilibration::new(chi, eta, total_charge);
        let positions = [[0.0; 3], [3.0, 0.0, 0.0], [0.0, 3.0, 0.0]];
        let charges = qeq.compute_charges(&positions, 200, 1e-6);
        let sum: f64 = charges.iter().sum();
        // The charge-rescaling step enforces sum == total_charge; use loose tolerance
        // because the Gauss-Seidel update may not fully converge in all cases.
        assert!(
            (sum - total_charge).abs() < 1e-6,
            "charges sum to {total_charge}: {sum}"
        );
    }

    #[test]
    fn test_qeq_j_matrix_symmetric() {
        let chi = vec![3.5, 4.5];
        let eta = vec![10.0, 11.0];
        let qeq = ChargeEquilibration::new(chi, eta, 0.0);
        let positions = [[0.0; 3], [3.0, 0.0, 0.0]];
        let j = qeq.build_j_matrix(&positions);
        assert!((j[0][1] - j[1][0]).abs() < 1e-12, "J should be symmetric");
    }

    #[test]
    fn test_qeq_diagonal_twice_eta() {
        let chi = vec![3.5, 4.5];
        let eta = vec![10.0, 11.0];
        let qeq = ChargeEquilibration::new(chi, eta.clone(), 0.0);
        let positions = [[0.0; 3], [3.0, 0.0, 0.0]];
        let j = qeq.build_j_matrix(&positions);
        assert!((j[0][0] - 2.0 * eta[0]).abs() < 1e-12, "J[0][0] = 2*eta[0]");
        assert!((j[1][1] - 2.0 * eta[1]).abs() < 1e-12, "J[1][1] = 2*eta[1]");
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
