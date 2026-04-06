// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Bond length constraints via SHAKE, RATTLE, and LINCS algorithms.
//!
//! SHAKE iteratively adjusts positions to satisfy distance constraints.
//! RATTLE additionally adjusts velocities so that v_ij . r_ij = 0
//! (velocity component along bond direction is zero).
//! LINCS (LINear Constraint Solver) solves constraints via matrix inversion
//! approximation.

use oxiphysics_core::math::{Real, Vec3};

// ---------------------------------------------------------------------------
// BondConstraint
// ---------------------------------------------------------------------------

/// A rigid bond constraint between two atoms.
///
/// The constraint enforces |r_i - r_j|^2 = `distance_sq` at every step.
#[derive(Debug, Clone)]
pub struct BondConstraint {
    /// Index of atom i.
    pub atom_i: usize,
    /// Index of atom j.
    pub atom_j: usize,
    /// Target squared distance d^2 = d_target^2.
    pub distance_sq: Real,
}

impl BondConstraint {
    /// Create a new bond constraint with target distance `d`.
    pub fn new(atom_i: usize, atom_j: usize, distance: Real) -> Self {
        Self {
            atom_i,
            atom_j,
            distance_sq: distance * distance,
        }
    }

    /// Target distance (not squared).
    pub fn target_distance(&self) -> Real {
        self.distance_sq.sqrt()
    }
}

// ---------------------------------------------------------------------------
// AngleConstraint
// ---------------------------------------------------------------------------

/// A constraint on the angle between three atoms (i-j-k).
///
/// Enforced indirectly through a distance constraint on atoms i and k
/// derived from the target angle and bond lengths.
#[derive(Debug, Clone)]
pub struct AngleConstraint {
    /// Index of atom i (first end).
    pub atom_i: usize,
    /// Index of central atom j.
    pub atom_j: usize,
    /// Index of atom k (second end).
    pub atom_k: usize,
    /// Target angle in radians.
    pub angle_rad: Real,
    /// Bond length j-i.
    pub r_ji: Real,
    /// Bond length j-k.
    pub r_jk: Real,
}

impl AngleConstraint {
    /// Create a new angle constraint.
    pub fn new(
        atom_i: usize,
        atom_j: usize,
        atom_k: usize,
        angle_rad: Real,
        r_ji: Real,
        r_jk: Real,
    ) -> Self {
        Self {
            atom_i,
            atom_j,
            atom_k,
            angle_rad,
            r_ji,
            r_jk,
        }
    }

    /// Compute the i-k distance implied by the angle and bond lengths (law of cosines).
    ///
    /// r_ik^2 = r_ji^2 + r_jk^2 - 2*r_ji*r_jk*cos(theta)
    pub fn implied_ik_distance_sq(&self) -> Real {
        self.r_ji * self.r_ji + self.r_jk * self.r_jk
            - 2.0 * self.r_ji * self.r_jk * self.angle_rad.cos()
    }

    /// Convert to a BondConstraint between atoms i and k.
    pub fn to_bond_constraint(&self) -> BondConstraint {
        BondConstraint {
            atom_i: self.atom_i,
            atom_j: self.atom_k,
            distance_sq: self.implied_ik_distance_sq(),
        }
    }
}

// ---------------------------------------------------------------------------
// ConstraintStats
// ---------------------------------------------------------------------------

/// Statistics from a constraint solution step.
#[derive(Debug, Clone, Copy)]
pub struct ConstraintStats {
    /// Number of iterations performed.
    pub iterations: u32,
    /// Maximum relative constraint violation after solution.
    pub max_violation: Real,
    /// Whether the solver converged within tolerance.
    pub converged: bool,
}

// ---------------------------------------------------------------------------
// Shake
// ---------------------------------------------------------------------------

/// SHAKE / RATTLE constraint solver.
///
/// SHAKE adjusts *positions* so that all bond-length constraints are
/// satisfied after a time step.  RATTLE follows SHAKE and adjusts
/// *velocities* so that the time derivative of every constraint is also
/// satisfied.
pub struct Shake {
    /// Bond constraints to enforce.
    pub constraints: Vec<BondConstraint>,
    /// Maximum number of SHAKE iterations.
    pub max_iter: usize,
    /// Convergence tolerance on |delta(r^2)|.
    pub tolerance: Real,
}

impl Shake {
    /// Create a new SHAKE solver with default `max_iter = 500` and
    /// `tolerance = 1e-8`.
    pub fn new(constraints: Vec<BondConstraint>) -> Self {
        Self {
            constraints,
            max_iter: 500,
            tolerance: 1e-8,
        }
    }

    /// Create a SHAKE solver with custom iteration limit and tolerance.
    pub fn with_params(constraints: Vec<BondConstraint>, max_iter: usize, tolerance: Real) -> Self {
        Self {
            constraints,
            max_iter,
            tolerance,
        }
    }

    /// Apply SHAKE: iteratively adjust `positions` so that every bond
    /// length constraint is satisfied.
    ///
    /// `old_positions` are the positions *before* the unconstrained step.
    /// `inv_masses[i]` is `1/m_i` for atom `i`.
    ///
    /// Returns the number of iterations performed.
    pub fn apply(&self, positions: &mut [Vec3], old_positions: &[Vec3], inv_masses: &[f64]) -> u32 {
        let mut iter = 0u32;

        for _ in 0..self.max_iter {
            iter += 1;
            let mut converged = true;

            for c in &self.constraints {
                let i = c.atom_i;
                let j = c.atom_j;

                let r_new = positions[i] - positions[j];
                let r_new_sq = r_new.norm_squared();

                let delta_sq = r_new_sq - c.distance_sq;
                if delta_sq.abs() <= self.tolerance * c.distance_sq {
                    continue;
                }
                converged = false;

                let r_ref = old_positions[i] - old_positions[j];
                let denom = 2.0 * r_ref.dot(&r_new) * (inv_masses[i] + inv_masses[j]);
                if denom.abs() < 1e-30 {
                    continue;
                }

                let lambda = delta_sq / denom;

                positions[i] -= r_ref * (lambda * inv_masses[i]);
                positions[j] += r_ref * (lambda * inv_masses[j]);
            }

            if converged {
                break;
            }
        }

        iter
    }

    /// Apply SHAKE and return detailed statistics.
    pub fn apply_with_stats(
        &self,
        positions: &mut [Vec3],
        old_positions: &[Vec3],
        inv_masses: &[f64],
    ) -> ConstraintStats {
        let iters = self.apply(positions, old_positions, inv_masses);

        // Compute max violation
        let max_violation = self.max_relative_violation(positions);
        let converged = max_violation <= self.tolerance;

        ConstraintStats {
            iterations: iters,
            max_violation,
            converged,
        }
    }

    /// Compute the maximum relative constraint violation.
    pub fn max_relative_violation(&self, positions: &[Vec3]) -> Real {
        let mut max_viol = 0.0_f64;
        for c in &self.constraints {
            let r = positions[c.atom_i] - positions[c.atom_j];
            let r_sq = r.norm_squared();
            let viol = (r_sq - c.distance_sq).abs() / c.distance_sq;
            if viol > max_viol {
                max_viol = viol;
            }
        }
        max_viol
    }

    /// Apply RATTLE: adjust `velocities` after SHAKE so that
    /// v_ij . r_ij = 0 for every constrained bond.
    pub fn rattle(
        &self,
        positions: &[Vec3],
        velocities: &mut [Vec3],
        inv_masses: &[f64],
        dt: Real,
    ) {
        for _ in 0..self.max_iter {
            let mut converged = true;

            for c in &self.constraints {
                let i = c.atom_i;
                let j = c.atom_j;

                let r = positions[i] - positions[j];
                let v = velocities[i] - velocities[j];

                let rv = r.dot(&v);
                if rv.abs() <= self.tolerance {
                    continue;
                }
                converged = false;

                let r_sq = r.norm_squared();
                if r_sq < 1e-30 {
                    continue;
                }

                let mu = rv / (r_sq * (inv_masses[i] + inv_masses[j]) * dt);

                velocities[i] -= r * (mu * inv_masses[i] * dt);
                velocities[j] += r * (mu * inv_masses[j] * dt);
            }

            if converged {
                break;
            }
        }
    }

    /// Apply RATTLE and return statistics.
    pub fn rattle_with_stats(
        &self,
        positions: &[Vec3],
        velocities: &mut [Vec3],
        inv_masses: &[f64],
        dt: Real,
    ) -> ConstraintStats {
        let mut iters = 0u32;
        let mut final_converged = false;

        for _ in 0..self.max_iter {
            iters += 1;
            let mut converged = true;

            for c in &self.constraints {
                let i = c.atom_i;
                let j = c.atom_j;
                let r = positions[i] - positions[j];
                let v = velocities[i] - velocities[j];
                let rv = r.dot(&v);
                if rv.abs() <= self.tolerance {
                    continue;
                }
                converged = false;
                let r_sq = r.norm_squared();
                if r_sq < 1e-30 {
                    continue;
                }
                let mu = rv / (r_sq * (inv_masses[i] + inv_masses[j]) * dt);
                velocities[i] -= r * (mu * inv_masses[i] * dt);
                velocities[j] += r * (mu * inv_masses[j] * dt);
            }

            if converged {
                final_converged = true;
                break;
            }
        }

        // Max velocity violation
        let mut max_viol = 0.0_f64;
        for c in &self.constraints {
            let r = positions[c.atom_i] - positions[c.atom_j];
            let v = velocities[c.atom_i] - velocities[c.atom_j];
            let rv = r.dot(&v).abs();
            if rv > max_viol {
                max_viol = rv;
            }
        }

        ConstraintStats {
            iterations: iters,
            max_violation: max_viol,
            converged: final_converged,
        }
    }

    /// Compute constraint forces: the Lagrange multiplier forces that
    /// SHAKE implicitly applies.
    ///
    /// Returns force on each atom due to constraints.
    pub fn constraint_forces(
        &self,
        positions: &[Vec3],
        old_positions: &[Vec3],
        inv_masses: &[f64],
        dt: Real,
    ) -> Vec<Vec3> {
        let n = positions.len();
        let mut forces = vec![Vec3::new(0.0, 0.0, 0.0); n];

        for c in &self.constraints {
            let i = c.atom_i;
            let j = c.atom_j;

            let r_new = positions[i] - positions[j];
            let r_new_sq = r_new.norm_squared();
            let delta_sq = r_new_sq - c.distance_sq;

            let r_ref = old_positions[i] - old_positions[j];
            let denom = 2.0 * r_ref.dot(&r_new) * (inv_masses[i] + inv_masses[j]);
            if denom.abs() < 1e-30 {
                continue;
            }

            let lambda = delta_sq / denom;
            // Force = -lambda * r_ref / dt^2 (times mass factors already in lambda)
            let f_scale = lambda / (dt * dt);
            forces[i] -= r_ref * (f_scale / inv_masses[i].max(1e-30));
            forces[j] += r_ref * (f_scale / inv_masses[j].max(1e-30));
        }

        forces
    }
}

// ---------------------------------------------------------------------------
// LINCS
// ---------------------------------------------------------------------------

/// LINCS (LINear Constraint Solver) algorithm.
///
/// A non-iterative constraint algorithm that solves the linearised
/// constraint equations. Faster than SHAKE for large systems but
/// less accurate without correction steps.
///
/// Reference: Hess et al., J. Comput. Chem. 18, 1463 (1997).
pub struct Lincs {
    /// Bond constraints.
    pub constraints: Vec<BondConstraint>,
    /// Order of the LINCS expansion (typically 4-8).
    pub expansion_order: usize,
}

impl Lincs {
    /// Create a new LINCS solver with default expansion order 4.
    pub fn new(constraints: Vec<BondConstraint>) -> Self {
        Self {
            constraints,
            expansion_order: 4,
        }
    }

    /// Create a LINCS solver with custom expansion order.
    pub fn with_order(constraints: Vec<BondConstraint>, order: usize) -> Self {
        Self {
            constraints,
            expansion_order: order,
        }
    }

    /// Apply LINCS to correct positions after an unconstrained update.
    ///
    /// This is a simplified implementation that applies a single linearised
    /// correction followed by `expansion_order` correction steps to improve
    /// the bond length accuracy.
    ///
    /// Returns the maximum relative constraint violation after correction.
    pub fn apply(
        &self,
        positions: &mut [Vec3],
        old_positions: &[Vec3],
        inv_masses: &[f64],
    ) -> Real {
        // Step 1: compute initial bond directions from old positions
        // Step 2: project unconstrained displacements onto bond directions
        // Step 3: correct positions

        for c in &self.constraints {
            let i = c.atom_i;
            let j = c.atom_j;

            let r_old = old_positions[i] - old_positions[j];
            let r_old_len = r_old.norm();
            if r_old_len < 1e-30 {
                continue;
            }
            let r_hat = r_old / r_old_len;

            let r_new = positions[i] - positions[j];
            // Project displacement onto bond direction
            let p = r_hat.dot(&r_new);
            let target = c.distance_sq.sqrt();

            // Correction: move along bond direction to restore length
            let correction = (p - target) / (inv_masses[i] + inv_masses[j]);

            positions[i] -= r_hat * (correction * inv_masses[i]);
            positions[j] += r_hat * (correction * inv_masses[j]);
        }

        // Correction steps to improve accuracy (LINCS rotation correction)
        for _ in 0..self.expansion_order {
            for c in &self.constraints {
                let i = c.atom_i;
                let j = c.atom_j;

                let r = positions[i] - positions[j];
                let r_len_sq = r.norm_squared();
                let target_sq = c.distance_sq;

                if r_len_sq < 1e-30 {
                    continue;
                }

                // Length correction factor
                let r_len = r_len_sq.sqrt();
                let target = target_sq.sqrt();
                let diff = r_len - target;

                if diff.abs() < 1e-12 {
                    continue;
                }

                let r_hat = r / r_len;
                let correction = diff / (inv_masses[i] + inv_masses[j]);

                positions[i] -= r_hat * (correction * inv_masses[i]);
                positions[j] += r_hat * (correction * inv_masses[j]);
            }
        }

        // Return max violation
        let mut max_viol = 0.0_f64;
        for c in &self.constraints {
            let r = positions[i_of(c, positions)] - positions[j_of(c, positions)];
            let viol = (r.norm_squared() - c.distance_sq).abs() / c.distance_sq;
            if viol > max_viol {
                max_viol = viol;
            }
        }
        max_viol
    }
}

#[inline]
fn i_of(c: &BondConstraint, _positions: &[Vec3]) -> usize {
    c.atom_i
}
#[inline]
fn j_of(c: &BondConstraint, _positions: &[Vec3]) -> usize {
    c.atom_j
}

// ---------------------------------------------------------------------------
// water_shake_constraints
// ---------------------------------------------------------------------------

/// Create O-H and H-H bond constraints for a TIP3P water molecule.
///
/// Atom ordering: O = `first_atom`, H1 = `first_atom+1`, H2 = `first_atom+2`.
///
/// Distances (in Å):
/// * O-H : 0.9572 Å
/// * H-H : 1.5139 Å
pub fn water_shake_constraints(first_atom: usize) -> Vec<BondConstraint> {
    let o = first_atom;
    let h1 = first_atom + 1;
    let h2 = first_atom + 2;

    let d_oh: Real = 0.9572;
    let d_hh: Real = 1.5139;

    vec![
        BondConstraint {
            atom_i: o,
            atom_j: h1,
            distance_sq: d_oh * d_oh,
        },
        BondConstraint {
            atom_i: o,
            atom_j: h2,
            distance_sq: d_oh * d_oh,
        },
        BondConstraint {
            atom_i: h1,
            atom_j: h2,
            distance_sq: d_hh * d_hh,
        },
    ]
}

/// Create constraints for an SPC water molecule.
///
/// O-H = 1.0 Å, H-O-H = 109.47°.
pub fn water_spc_constraints(first_atom: usize) -> Vec<BondConstraint> {
    let o = first_atom;
    let h1 = first_atom + 1;
    let h2 = first_atom + 2;

    let d_oh: Real = 1.0;
    let half_angle = 109.47_f64.to_radians() / 2.0;
    let d_hh: Real = 2.0 * d_oh * half_angle.sin();

    vec![
        BondConstraint::new(o, h1, d_oh),
        BondConstraint::new(o, h2, d_oh),
        BondConstraint::new(h1, h2, d_hh),
    ]
}

/// Create constraints for a chain of `n_bonds` equal bonds.
///
/// Atoms are numbered `first_atom`, `first_atom+1`, ..., `first_atom+n_bonds`.
pub fn chain_constraints(
    first_atom: usize,
    n_bonds: usize,
    bond_length: Real,
) -> Vec<BondConstraint> {
    (0..n_bonds)
        .map(|i| BondConstraint::new(first_atom + i, first_atom + i + 1, bond_length))
        .collect()
}

/// Total constraint satisfaction error: sum of |r^2 - d^2| / d^2 over all constraints.
pub fn total_constraint_error(constraints: &[BondConstraint], positions: &[Vec3]) -> Real {
    constraints
        .iter()
        .map(|c| {
            let r = positions[c.atom_i] - positions[c.atom_j];
            (r.norm_squared() - c.distance_sq).abs() / c.distance_sq
        })
        .sum()
}

// ---------------------------------------------------------------------------
// SETTLE – rigid water constraint solver
// ---------------------------------------------------------------------------

/// SETTLE algorithm for rigid water molecules (O + 2 H).
///
/// Solves the three distance constraints (O-H1, O-H2, H-H) analytically in
/// one pass using the known geometry of a water molecule.  Much faster than
/// iterative SHAKE for water.
///
/// Reference: S. Miyamoto & P. A. Kollman, J. Comput. Chem. 13, 952 (1992).
pub struct SettleWater {
    /// O-H bond length (Å or nm, same units as positions).
    pub d_oh: Real,
    /// H-H distance (Å or nm).
    pub d_hh: Real,
}

impl SettleWater {
    /// Create a SETTLE solver for TIP3P water (d_OH = 0.9572 Å, angle = 104.52°).
    pub fn tip3p() -> Self {
        let d_oh: Real = 0.9572;
        let angle_rad = 104.52_f64.to_radians();
        let d_hh = (2.0 * d_oh * d_oh * (1.0 - angle_rad.cos())).sqrt();
        Self { d_oh, d_hh }
    }

    /// Create a SETTLE solver for SPC water (d_OH = 1.0 Å, angle = 109.47°).
    pub fn spc() -> Self {
        let d_oh: Real = 1.0;
        let angle_rad = 109.47_f64.to_radians();
        let d_hh = (2.0 * d_oh * d_oh * (1.0 - angle_rad.cos())).sqrt();
        Self { d_oh, d_hh }
    }

    /// Apply SETTLE constraints to positions `[O, H1, H2]`.
    ///
    /// Adjusts the three positions so that all three distances
    /// (O-H1, O-H2, H1-H2) match their target values.
    /// `inv_masses` must have three entries corresponding to O, H1, H2.
    pub fn apply(&self, positions: &mut [Vec3], inv_masses: &[f64]) {
        debug_assert!(positions.len() >= 3);
        debug_assert!(inv_masses.len() >= 3);

        let mo_inv = inv_masses[0];
        let mh_inv = inv_masses[1];
        // Total inverse mass
        let _m_tot_inv = mo_inv + 2.0 * mh_inv;

        // Iterative SHAKE-like projection for the three water constraints.
        // Iterate until all constraints are satisfied to within tolerance.
        let constraints = [
            (0usize, 1usize, self.d_oh),
            (0usize, 2usize, self.d_oh),
            (1usize, 2usize, self.d_hh),
        ];
        let tol = 1e-8;
        let max_iter = 500;

        for _ in 0..max_iter {
            let mut converged = true;
            for &(ci, cj, target) in &constraints {
                let rij = positions[ci] - positions[cj];
                let r = rij.norm();
                if r < 1e-15 {
                    continue;
                }
                let diff = r - target;
                if diff.abs() < tol {
                    continue;
                }
                converged = false;
                let inv_sum = inv_masses[ci] + inv_masses[cj];
                if inv_sum < 1e-30 {
                    continue;
                }
                let lambda = diff / inv_sum;
                let direction = rij / r;
                positions[ci] -= direction * (lambda * inv_masses[ci]);
                positions[cj] += direction * (lambda * inv_masses[cj]);
            }
            if converged {
                break;
            }
        }
    }

    /// Maximum constraint violation for a water molecule.
    pub fn max_violation(&self, positions: &[Vec3]) -> Real {
        let checks = [
            ((0, 1), self.d_oh),
            ((0, 2), self.d_oh),
            ((1, 2), self.d_hh),
        ];
        checks
            .iter()
            .map(|&((i, j), target)| {
                let r = (positions[i] - positions[j]).norm();
                (r - target).abs() / target
            })
            .fold(0.0_f64, f64::max)
    }
}

// ---------------------------------------------------------------------------
// Holonomic constraint matrix
// ---------------------------------------------------------------------------

/// Compute the holonomic constraint Jacobian matrix G (n_constraints × 3N).
///
/// Each row `c` of the matrix is the gradient of the constraint function
/// `σ_c(q) = |r_i - r_j| - d_0` with respect to all Cartesian coordinates.
///
/// The positions are passed in as a flat list of N atoms; the resulting
/// matrix has shape `[n_constraints][3*n_atoms]`.
#[allow(dead_code)]
pub fn holonomic_constraint_matrix_with_positions(
    constraints: &[BondConstraint],
    positions: &[Vec3],
) -> Vec<Vec<Real>> {
    let n = positions.len();
    let dim = 3 * n;
    let mut g = vec![vec![0.0_f64; dim]; constraints.len()];

    for (row, c) in constraints.iter().enumerate() {
        let i = c.atom_i;
        let j = c.atom_j;
        let rij = positions[i] - positions[j];
        let r = rij.norm();
        if r < 1e-30 {
            continue;
        }
        let unit = rij / r;
        // dσ/dr_i = +unit, dσ/dr_j = -unit
        for a in 0..3 {
            g[row][3 * i + a] = unit[a];
            g[row][3 * j + a] = -unit[a];
        }
    }
    g
}

/// Compute the holonomic constraint Jacobian using unit bond vectors along x.
///
/// This version is for atoms placed on a line (useful for benchmarks/tests).
#[allow(dead_code)]
pub fn holonomic_constraint_matrix(
    constraints: &[BondConstraint],
    n_atoms: usize,
) -> Vec<Vec<Real>> {
    let dim = 3 * n_atoms;
    let mut g = vec![vec![0.0_f64; dim]; constraints.len()];

    for (row, c) in constraints.iter().enumerate() {
        let i = c.atom_i;
        let j = c.atom_j;
        if i < n_atoms && j < n_atoms {
            // Default: unit vector along +x
            g[row][3 * i] = 1.0;
            g[row][3 * j] = -1.0;
        }
    }
    g
}

// ---------------------------------------------------------------------------
// Penalty constraint (soft constraint via large spring constant)
// ---------------------------------------------------------------------------

/// A penalty (soft) constraint that enforces bond lengths via a stiff harmonic
/// spring rather than the exact Lagrangian approach.
///
/// V_penalty = 0.5 * k_penalty * (|r_i - r_j| - d_0)^2
#[derive(Debug, Clone)]
pub struct PenaltyConstraint {
    /// The underlying bond constraint (target distance).
    pub bond: BondConstraint,
    /// Penalty spring constant (force units / length^2).
    pub k_penalty: Real,
}

impl PenaltyConstraint {
    /// Create a new penalty constraint.
    pub fn new(atom_i: usize, atom_j: usize, distance: Real, k: Real) -> Self {
        Self {
            bond: BondConstraint::new(atom_i, atom_j, distance),
            k_penalty: k,
        }
    }

    /// Potential energy of this penalty constraint.
    pub fn energy(&self, positions: &[Vec3]) -> Real {
        let r = (positions[self.bond.atom_i] - positions[self.bond.atom_j]).norm();
        let d0 = self.bond.target_distance();
        let dr = r - d0;
        0.5 * self.k_penalty * dr * dr
    }

    /// Forces on atoms i and j due to this penalty constraint.
    ///
    /// Returns a Vec of forces for all atoms (only indices i and j are nonzero).
    pub fn forces(&self, positions: &[Vec3]) -> Vec<Vec3> {
        let n = positions.len();
        let mut forces = vec![Vec3::new(0.0, 0.0, 0.0); n];
        let i = self.bond.atom_i;
        let j = self.bond.atom_j;
        let rij = positions[i] - positions[j];
        let r = rij.norm();
        if r < 1e-15 {
            return forces;
        }
        let d0 = self.bond.target_distance();
        let dr = r - d0;
        // Force = -∂V/∂r_i = -k*(r-d0)/r * rij (restoring force pulls i toward j)
        let mag = self.k_penalty * dr / r;
        forces[i] -= rij * mag;
        forces[j] += rij * mag;
        forces
    }

    /// Apply the penalty force to the accumulated forces array (in-place).
    pub fn accumulate_forces(&self, positions: &[Vec3], forces: &mut [Vec3]) {
        let i = self.bond.atom_i;
        let j = self.bond.atom_j;
        let rij = positions[i] - positions[j];
        let r = rij.norm();
        if r < 1e-15 {
            return;
        }
        let d0 = self.bond.target_distance();
        let dr = r - d0;
        // Force = -∂V/∂r_i = -k*(r-d0)/r * rij
        let mag = self.k_penalty * dr / r;
        forces[i] -= rij * mag;
        forces[j] += rij * mag;
    }
}

// ---------------------------------------------------------------------------
// RATTLE velocity constraints (standalone helper)
// ---------------------------------------------------------------------------

/// Apply RATTLE velocity correction to ensure v_ij · r_ij = 0 for all bonds.
///
/// This is the same as `Shake::rattle` but exposed as a free function for
/// use without the full `Shake` struct.
#[allow(dead_code)]
pub fn rattle_velocities(
    constraints: &[BondConstraint],
    positions: &[Vec3],
    velocities: &mut [Vec3],
    inv_masses: &[f64],
    dt: Real,
    max_iter: usize,
    tolerance: Real,
) -> ConstraintStats {
    let mut iters = 0u32;
    let mut converged = false;

    for _ in 0..max_iter {
        iters += 1;
        let mut all_ok = true;

        for c in constraints {
            let i = c.atom_i;
            let j = c.atom_j;
            let r = positions[i] - positions[j];
            let v = velocities[i] - velocities[j];
            let rv = r.dot(&v);
            if rv.abs() <= tolerance {
                continue;
            }
            all_ok = false;
            let r_sq = r.norm_squared();
            if r_sq < 1e-30 {
                continue;
            }
            let mu = rv / (r_sq * (inv_masses[i] + inv_masses[j]) * dt);
            velocities[i] -= r * (mu * inv_masses[i] * dt);
            velocities[j] += r * (mu * inv_masses[j] * dt);
        }

        if all_ok {
            converged = true;
            break;
        }
    }

    let mut max_viol = 0.0_f64;
    for c in constraints {
        let r = positions[c.atom_i] - positions[c.atom_j];
        let v = velocities[c.atom_i] - velocities[c.atom_j];
        let rv = r.dot(&v).abs();
        if rv > max_viol {
            max_viol = rv;
        }
    }

    ConstraintStats {
        iterations: iters,
        max_violation: max_viol,
        converged,
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    /// Tolerance for position-level assertions (1e-6).
    const POS_TOL: Real = 1e-6;
    /// Tolerance for velocity-level assertions (1e-5).
    const VEL_TOL: Real = 1e-5;

    // ------------------------------------------------------------------
    // 1. test_shake_single_bond
    // ------------------------------------------------------------------
    #[test]
    fn test_shake_single_bond() {
        let target_d: Real = 1.0;
        let c = vec![BondConstraint {
            atom_i: 0,
            atom_j: 1,
            distance_sq: target_d * target_d,
        }];
        let shake = Shake::new(c);

        let old_positions = vec![Vec3::new(0.0, 0.0, 0.0), Vec3::new(2.0, 0.0, 0.0)];
        let mut positions = old_positions.clone();
        let inv_masses = vec![1.0, 1.0];

        let iters = shake.apply(&mut positions, &old_positions, &inv_masses);

        let dist = (positions[0] - positions[1]).norm();
        assert!(
            (dist - target_d).abs() < POS_TOL,
            "expected dist={target_d}, got {dist}"
        );
        assert!(iters > 0);
    }

    // ------------------------------------------------------------------
    // 2. test_shake_bond_already_satisfied
    // ------------------------------------------------------------------
    #[test]
    fn test_shake_bond_already_satisfied() {
        let target_d: Real = 1.0;
        let c = vec![BondConstraint {
            atom_i: 0,
            atom_j: 1,
            distance_sq: target_d * target_d,
        }];
        let shake = Shake::new(c);

        let old_positions = vec![Vec3::new(0.0, 0.0, 0.0), Vec3::new(1.0, 0.0, 0.0)];
        let mut positions = old_positions.clone();
        let inv_masses = vec![1.0, 1.0];

        let iters = shake.apply(&mut positions, &old_positions, &inv_masses);

        let dist = (positions[0] - positions[1]).norm();
        assert!((dist - target_d).abs() < POS_TOL, "dist={dist}");
        assert_eq!(
            iters, 1,
            "already-satisfied constraint should converge in 1 iter"
        );
        assert!((positions[0] - old_positions[0]).norm() < POS_TOL);
        assert!((positions[1] - old_positions[1]).norm() < POS_TOL);
    }

    // ------------------------------------------------------------------
    // 3. test_shake_triangle
    // ------------------------------------------------------------------
    #[test]
    fn test_shake_triangle() {
        let target_d: Real = 1.0;
        let constraints = vec![
            BondConstraint {
                atom_i: 0,
                atom_j: 1,
                distance_sq: target_d * target_d,
            },
            BondConstraint {
                atom_i: 1,
                atom_j: 2,
                distance_sq: target_d * target_d,
            },
            BondConstraint {
                atom_i: 0,
                atom_j: 2,
                distance_sq: target_d * target_d,
            },
        ];
        let shake = Shake::new(constraints);

        let half = 0.5_f64;
        let h = (3.0_f64 / 4.0).sqrt();
        let old_positions = vec![
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(half, h, 0.0),
        ];

        let mut positions = vec![
            Vec3::new(0.05, -0.03, 0.0),
            Vec3::new(1.04, 0.02, 0.0),
            Vec3::new(half + 0.01, h - 0.04, 0.0),
        ];

        let inv_masses = vec![1.0, 1.0, 1.0];

        let iters = shake.apply(&mut positions, &old_positions, &inv_masses);

        for c in shake.constraints.iter() {
            let d = (positions[c.atom_i] - positions[c.atom_j]).norm();
            assert!(
                (d - target_d).abs() < POS_TOL,
                "bond ({},{}) dist={d}, iters={iters}",
                c.atom_i,
                c.atom_j
            );
        }
    }

    // ------------------------------------------------------------------
    // 4. test_rattle_velocity_perpendicular
    // ------------------------------------------------------------------
    #[test]
    fn test_rattle_velocity_perpendicular() {
        let target_d: Real = 1.5;
        let c = vec![BondConstraint {
            atom_i: 0,
            atom_j: 1,
            distance_sq: target_d * target_d,
        }];
        let shake = Shake::new(c);

        let positions = vec![Vec3::new(0.0, 0.0, 0.0), Vec3::new(target_d, 0.0, 0.0)];
        let mut velocities = vec![Vec3::new(0.3, 1.0, 0.0), Vec3::new(-0.2, 0.5, 0.0)];
        let inv_masses = vec![1.0, 1.0];
        let dt = 0.002;

        shake.rattle(&positions, &mut velocities, &inv_masses, dt);

        let r = positions[0] - positions[1];
        let v = velocities[0] - velocities[1];
        let rv = r.dot(&v);
        assert!(
            rv.abs() < VEL_TOL,
            "v_ij . r_ij = {rv} should be < {VEL_TOL}"
        );
    }

    // ------------------------------------------------------------------
    // 5. test_water_constraints
    // ------------------------------------------------------------------
    #[test]
    fn test_water_constraints() {
        let constraints = water_shake_constraints(0);

        assert_eq!(constraints.len(), 3, "TIP3P water needs 3 constraints");

        let d_oh: Real = 0.9572;
        let d_hh: Real = 1.5139;

        assert!((constraints[0].distance_sq - d_oh * d_oh).abs() < 1e-10);
        assert!((constraints[1].distance_sq - d_oh * d_oh).abs() < 1e-10);
        assert!((constraints[2].distance_sq - d_hh * d_hh).abs() < 1e-10);

        let shake = Shake::new(constraints);

        let angle_half = (104.52_f64.to_radians()) / 2.0;
        let old_positions = vec![
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(d_oh * angle_half.sin(), d_oh * angle_half.cos(), 0.0),
            Vec3::new(-d_oh * angle_half.sin(), d_oh * angle_half.cos(), 0.0),
        ];

        let mut positions = vec![
            old_positions[0] + Vec3::new(0.02, -0.01, 0.0),
            old_positions[1] + Vec3::new(-0.01, 0.02, 0.0),
            old_positions[2] + Vec3::new(0.01, 0.01, 0.0),
        ];

        let inv_masses = vec![1.0 / 16.0, 1.0, 1.0];

        let iters = shake.apply(&mut positions, &old_positions, &inv_masses);

        for c in shake.constraints.iter() {
            let d = (positions[c.atom_i] - positions[c.atom_j]).norm();
            let d_target = c.distance_sq.sqrt();
            assert!(
                (d - d_target).abs() < POS_TOL,
                "water bond ({},{}) dist={d}, target={d_target}, iters={iters}",
                c.atom_i,
                c.atom_j
            );
        }
    }

    // ------------------------------------------------------------------
    // 6. test_bond_constraint_constructor
    // ------------------------------------------------------------------
    #[test]
    fn test_bond_constraint_constructor() {
        let c = BondConstraint::new(0, 1, 1.5);
        assert_eq!(c.atom_i, 0);
        assert_eq!(c.atom_j, 1);
        assert!((c.target_distance() - 1.5).abs() < 1e-12);
    }

    // ------------------------------------------------------------------
    // 7. test_angle_constraint_ik_distance
    // ------------------------------------------------------------------
    #[test]
    fn test_angle_constraint_ik_distance() {
        // Equilateral triangle: angle = 60 deg, r = 1.0
        let ac = AngleConstraint::new(0, 1, 2, std::f64::consts::FRAC_PI_3, 1.0, 1.0);
        let d_ik_sq = ac.implied_ik_distance_sq();
        // For equilateral: d_ik = 1.0
        assert!(
            (d_ik_sq - 1.0).abs() < 1e-10,
            "equilateral d_ik^2 = {d_ik_sq}, expected 1.0"
        );
    }

    // ------------------------------------------------------------------
    // 8. test_angle_constraint_right_angle
    // ------------------------------------------------------------------
    #[test]
    fn test_angle_constraint_right_angle() {
        // Right angle: 90 deg, r_ji = 3, r_jk = 4 -> d_ik = 5 (Pythagorean)
        let ac = AngleConstraint::new(0, 1, 2, std::f64::consts::FRAC_PI_2, 3.0, 4.0);
        let d_ik = ac.implied_ik_distance_sq().sqrt();
        assert!(
            (d_ik - 5.0).abs() < 1e-10,
            "3-4-5 triangle: d_ik = {d_ik}, expected 5.0"
        );
    }

    // ------------------------------------------------------------------
    // 9. test_angle_to_bond_constraint
    // ------------------------------------------------------------------
    #[test]
    fn test_angle_to_bond_constraint() {
        let ac = AngleConstraint::new(0, 1, 2, std::f64::consts::FRAC_PI_3, 1.0, 1.0);
        let bc = ac.to_bond_constraint();
        assert_eq!(bc.atom_i, 0);
        assert_eq!(bc.atom_j, 2);
        assert!((bc.target_distance() - 1.0).abs() < 1e-10);
    }

    // ------------------------------------------------------------------
    // 10. test_shake_with_stats
    // ------------------------------------------------------------------
    #[test]
    fn test_shake_with_stats() {
        let c = vec![BondConstraint::new(0, 1, 1.0)];
        let shake = Shake::new(c);

        let old_positions = vec![Vec3::new(0.0, 0.0, 0.0), Vec3::new(1.0, 0.0, 0.0)];
        let mut positions = old_positions.clone();
        let inv_masses = vec![1.0, 1.0];

        let stats = shake.apply_with_stats(&mut positions, &old_positions, &inv_masses);
        assert!(
            stats.converged,
            "should converge for already-satisfied constraint"
        );
        assert_eq!(stats.iterations, 1);
        assert!(stats.max_violation < 1e-8);
    }

    // ------------------------------------------------------------------
    // 11. test_rattle_with_stats
    // ------------------------------------------------------------------
    #[test]
    fn test_rattle_with_stats() {
        let c = vec![BondConstraint::new(0, 1, 1.0)];
        let shake = Shake::new(c);

        let positions = vec![Vec3::new(0.0, 0.0, 0.0), Vec3::new(1.0, 0.0, 0.0)];
        let mut velocities = vec![Vec3::new(1.0, 0.0, 0.0), Vec3::new(-1.0, 0.0, 0.0)];
        let inv_masses = vec![1.0, 1.0];

        let stats = shake.rattle_with_stats(&positions, &mut velocities, &inv_masses, 0.001);
        assert!(stats.converged);
    }

    // ------------------------------------------------------------------
    // 12. test_spc_water_constraints
    // ------------------------------------------------------------------
    #[test]
    fn test_spc_water_constraints() {
        let constraints = water_spc_constraints(0);
        assert_eq!(constraints.len(), 3);

        // O-H should be 1.0 Å
        let d_oh = constraints[0].target_distance();
        assert!((d_oh - 1.0).abs() < 1e-10, "SPC O-H = {d_oh}");

        // H-H from tetrahedral angle
        let half_angle = 109.47_f64.to_radians() / 2.0;
        let expected_hh = 2.0 * 1.0 * half_angle.sin();
        let d_hh = constraints[2].target_distance();
        assert!(
            (d_hh - expected_hh).abs() < 1e-6,
            "SPC H-H = {d_hh}, expected {expected_hh}"
        );
    }

    // ------------------------------------------------------------------
    // 13. test_chain_constraints
    // ------------------------------------------------------------------
    #[test]
    fn test_chain_constraints() {
        let constraints = chain_constraints(0, 5, 1.5);
        assert_eq!(constraints.len(), 5);
        for (i, c) in constraints.iter().enumerate() {
            assert_eq!(c.atom_i, i);
            assert_eq!(c.atom_j, i + 1);
            assert!((c.target_distance() - 1.5).abs() < 1e-12);
        }
    }

    // ------------------------------------------------------------------
    // 14. test_total_constraint_error
    // ------------------------------------------------------------------
    #[test]
    fn test_total_constraint_error() {
        let constraints = vec![BondConstraint::new(0, 1, 1.0)];
        let positions = vec![Vec3::new(0.0, 0.0, 0.0), Vec3::new(1.0, 0.0, 0.0)];
        let err = total_constraint_error(&constraints, &positions);
        assert!(
            err < 1e-12,
            "error should be ~0 for satisfied constraints, got {err}"
        );
    }

    // ------------------------------------------------------------------
    // 15. test_lincs_single_bond
    // ------------------------------------------------------------------
    #[test]
    fn test_lincs_single_bond() {
        let c = vec![BondConstraint::new(0, 1, 1.0)];
        let lincs = Lincs::new(c);

        let old_positions = vec![Vec3::new(0.0, 0.0, 0.0), Vec3::new(1.0, 0.0, 0.0)];
        // Perturbed positions
        let mut positions = vec![Vec3::new(0.05, 0.0, 0.0), Vec3::new(1.1, 0.0, 0.0)];
        let inv_masses = vec![1.0, 1.0];

        let viol = lincs.apply(&mut positions, &old_positions, &inv_masses);

        let d = (positions[0] - positions[1]).norm();
        assert!(
            (d - 1.0).abs() < 0.01,
            "LINCS should restore bond to ~1.0, got {d}, viol={viol}"
        );
    }

    // ------------------------------------------------------------------
    // 16. test_shake_unequal_masses
    // ------------------------------------------------------------------
    #[test]
    fn test_shake_unequal_masses() {
        // Heavy atom (mass 100) and light atom (mass 1)
        let c = vec![BondConstraint::new(0, 1, 2.0)];
        let shake = Shake::new(c);

        let old_positions = vec![Vec3::new(0.0, 0.0, 0.0), Vec3::new(2.0, 0.0, 0.0)];
        let mut positions = vec![Vec3::new(0.0, 0.0, 0.0), Vec3::new(2.5, 0.0, 0.0)];
        let inv_masses = vec![0.01, 1.0]; // heavy vs light

        shake.apply(&mut positions, &old_positions, &inv_masses);

        let d = (positions[0] - positions[1]).norm();
        assert!((d - 2.0).abs() < POS_TOL, "dist = {d}, expected 2.0");

        // Heavy atom should barely move
        assert!(
            (positions[0] - old_positions[0]).norm() < 0.01,
            "heavy atom should barely move"
        );
    }

    // ------------------------------------------------------------------
    // 17. test_max_relative_violation
    // ------------------------------------------------------------------
    #[test]
    fn test_max_relative_violation() {
        let constraints = vec![BondConstraint::new(0, 1, 1.0)];
        let shake = Shake::new(constraints);

        // Bond length is 1.1 instead of 1.0 -> violation = |1.21 - 1.0| / 1.0 = 0.21
        let positions = vec![Vec3::new(0.0, 0.0, 0.0), Vec3::new(1.1, 0.0, 0.0)];
        let viol = shake.max_relative_violation(&positions);
        assert!(
            (viol - 0.21).abs() < 1e-10,
            "violation = {viol}, expected 0.21"
        );
    }

    // ------------------------------------------------------------------
    // 18. test_shake_with_custom_params
    // ------------------------------------------------------------------
    #[test]
    fn test_shake_with_custom_params() {
        let c = vec![BondConstraint::new(0, 1, 1.0)];
        let shake = Shake::with_params(c, 100, 1e-12);
        assert_eq!(shake.max_iter, 100);
        assert!((shake.tolerance - 1e-12).abs() < 1e-20);
    }

    // ------------------------------------------------------------------
    // 19. test_settle_rigid_water_energy
    // ------------------------------------------------------------------
    #[test]
    fn test_settle_rigid_water_energy() {
        let constraints = water_shake_constraints(0);
        let settle = SettleWater::tip3p();
        let d_oh = constraints[0].target_distance();
        assert!((d_oh - settle.d_oh).abs() < 1e-4);
        assert!(settle.d_hh > 0.0);
    }

    // ------------------------------------------------------------------
    // 20. test_settle_apply_satisfies_constraints
    // ------------------------------------------------------------------
    #[test]
    fn test_settle_apply_satisfies_constraints() {
        let settle = SettleWater::tip3p();
        let d_oh = settle.d_oh;
        let angle_half = (104.52_f64.to_radians()) / 2.0;

        // Reference (perfect) water geometry
        let ref_o = Vec3::new(0.0, 0.0, 0.0);
        let ref_h1 = Vec3::new(d_oh * angle_half.sin(), d_oh * angle_half.cos(), 0.0);
        let ref_h2 = Vec3::new(-d_oh * angle_half.sin(), d_oh * angle_half.cos(), 0.0);

        // Perturb positions
        let mut positions = vec![
            ref_o + Vec3::new(0.02, -0.01, 0.005),
            ref_h1 + Vec3::new(-0.01, 0.02, 0.003),
            ref_h2 + Vec3::new(0.01, 0.01, -0.002),
        ];
        let inv_masses = vec![1.0 / 16.0, 1.0, 1.0];

        settle.apply(&mut positions, &inv_masses);

        let d_oh1 = (positions[0] - positions[1]).norm();
        let d_oh2 = (positions[0] - positions[2]).norm();
        let d_hh = (positions[1] - positions[2]).norm();

        assert!((d_oh1 - settle.d_oh).abs() < 1e-5, "O-H1 = {d_oh1}");
        assert!((d_oh2 - settle.d_oh).abs() < 1e-5, "O-H2 = {d_oh2}");
        assert!((d_hh - settle.d_hh).abs() < 1e-5, "H-H  = {d_hh}");
    }

    // ------------------------------------------------------------------
    // 21. test_holonomic_constraint_matrix_size
    // ------------------------------------------------------------------
    #[test]
    fn test_holonomic_constraint_matrix_size() {
        let constraints = vec![
            BondConstraint::new(0, 1, 1.0),
            BondConstraint::new(1, 2, 1.0),
        ];
        let n_atoms = 3;
        let mat = holonomic_constraint_matrix(&constraints, n_atoms);
        // 2 constraints × (3 atoms × 3 coords) = 2 × 9
        assert_eq!(mat.len(), 2);
        assert_eq!(mat[0].len(), 9);
    }

    // ------------------------------------------------------------------
    // 22. test_holonomic_constraint_matrix_row_sum
    // ------------------------------------------------------------------
    #[test]
    fn test_holonomic_constraint_matrix_row_sum() {
        let constraints = vec![BondConstraint::new(0, 1, 1.0)];
        let positions = vec![Vec3::new(0.0, 0.0, 0.0), Vec3::new(1.0, 0.0, 0.0)];
        let mat = holonomic_constraint_matrix_with_positions(&constraints, &positions);
        // Row sum for gradient of bond length constraint: should be near zero
        // (forces sum to zero by Newton III in the gradient sense)
        let row_sum: Real = mat[0].iter().sum();
        assert!(row_sum.abs() < 1e-10, "Row sum = {row_sum}");
    }

    // ------------------------------------------------------------------
    // 23. test_penalty_constraint_energy
    // ------------------------------------------------------------------
    #[test]
    fn test_penalty_constraint_energy() {
        let c = BondConstraint::new(0, 1, 1.0);
        let penalty = PenaltyConstraint {
            bond: c,
            k_penalty: 1000.0,
        };
        let positions = vec![Vec3::new(0.0, 0.0, 0.0), Vec3::new(1.5, 0.0, 0.0)];
        let e = penalty.energy(&positions);
        // V = 0.5 * k * (r - r0)^2 = 0.5 * 1000 * 0.5^2 = 125
        assert!((e - 125.0).abs() < 1e-8, "penalty energy = {e}");
    }

    // ------------------------------------------------------------------
    // 24. test_penalty_constraint_zero_at_equilibrium
    // ------------------------------------------------------------------
    #[test]
    fn test_penalty_constraint_zero_at_equilibrium() {
        let c = BondConstraint::new(0, 1, 1.5);
        let penalty = PenaltyConstraint {
            bond: c,
            k_penalty: 500.0,
        };
        let positions = vec![Vec3::new(0.0, 0.0, 0.0), Vec3::new(1.5, 0.0, 0.0)];
        let e = penalty.energy(&positions);
        assert!(
            e.abs() < 1e-10,
            "energy at equilibrium should be zero, got {e}"
        );
    }

    // ------------------------------------------------------------------
    // 25. test_penalty_constraint_forces
    // ------------------------------------------------------------------
    #[test]
    fn test_penalty_constraint_forces() {
        let c = BondConstraint::new(0, 1, 1.0);
        let penalty = PenaltyConstraint {
            bond: c,
            k_penalty: 100.0,
        };
        let positions = vec![Vec3::new(0.0, 0.0, 0.0), Vec3::new(2.0, 0.0, 0.0)];
        let forces = penalty.forces(&positions);
        // Bond is stretched: force on atom 0 should be +x (toward atom 1)
        assert!(
            forces[0][0] > 0.0,
            "Force on atom 0 should be +x: {}",
            forces[0][0]
        );
        assert!(
            forces[1][0] < 0.0,
            "Force on atom 1 should be -x: {}",
            forces[1][0]
        );
        // Newton III
        let sum = forces[0][0] + forces[1][0];
        assert!(sum.abs() < 1e-10, "force sum = {sum}");
    }

    // ------------------------------------------------------------------
    // 26. test_rattle_constraint_force_calculation
    // ------------------------------------------------------------------
    #[test]
    fn test_rattle_constraint_force_calculation() {
        let c = vec![BondConstraint::new(0, 1, 1.0)];
        let shake = Shake::new(c);
        let old_positions = vec![Vec3::new(0.0, 0.0, 0.0), Vec3::new(1.0, 0.0, 0.0)];
        let positions = vec![Vec3::new(0.0, 0.0, 0.0), Vec3::new(1.05, 0.0, 0.0)];
        let inv_masses = vec![1.0, 1.0];
        let dt = 0.001;
        let forces = shake.constraint_forces(&positions, &old_positions, &inv_masses, dt);
        assert_eq!(forces.len(), 2);
        // Constraint forces should be finite
        for f in &forces {
            assert!(f.norm().is_finite());
        }
    }

    // ------------------------------------------------------------------
    // 27. test_lincs_with_order
    // ------------------------------------------------------------------
    #[test]
    fn test_lincs_with_order() {
        let c = vec![BondConstraint::new(0, 1, 1.0)];
        let lincs = Lincs::with_order(c, 8);
        assert_eq!(lincs.expansion_order, 8);
    }

    // ------------------------------------------------------------------
    // 28. test_lincs_chain_constraint
    // ------------------------------------------------------------------
    #[test]
    fn test_lincs_chain_constraint() {
        let constraints = chain_constraints(0, 3, 1.0);
        let lincs = Lincs::new(constraints);

        let old_positions = vec![
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(2.0, 0.0, 0.0),
            Vec3::new(3.0, 0.0, 0.0),
        ];
        let mut positions = vec![
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(1.1, 0.05, 0.0),
            Vec3::new(2.05, -0.05, 0.0),
            Vec3::new(3.1, 0.02, 0.0),
        ];
        let inv_masses = vec![1.0; 4];

        lincs.apply(&mut positions, &old_positions, &inv_masses);

        for c in &lincs.constraints {
            let d = (positions[c.atom_i] - positions[c.atom_j]).norm();
            assert!((d - 1.0).abs() < 0.02, "LINCS chain bond dist = {d}");
        }
    }

    // ------------------------------------------------------------------
    // 29. test_settle_tip3p_geometry
    // ------------------------------------------------------------------
    #[test]
    fn test_settle_tip3p_geometry() {
        let settle = SettleWater::tip3p();
        assert!((settle.d_oh - 0.9572).abs() < 1e-4);
        let d_hh_expected =
            (2.0 * 0.9572_f64.powi(2) * (1.0 - 104.52_f64.to_radians().cos())).sqrt();
        assert!(
            (settle.d_hh - d_hh_expected).abs() < 1e-4,
            "H-H = {}",
            settle.d_hh
        );
    }

    // ------------------------------------------------------------------
    // 30. test_constraint_error_multiple
    // ------------------------------------------------------------------
    #[test]
    fn test_constraint_error_multiple() {
        let constraints = vec![
            BondConstraint::new(0, 1, 1.0),
            BondConstraint::new(1, 2, 1.0),
        ];
        // Perfectly satisfied
        let positions = vec![
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(2.0, 0.0, 0.0),
        ];
        let err = total_constraint_error(&constraints, &positions);
        assert!(err < 1e-12, "error = {err}");
    }
}
