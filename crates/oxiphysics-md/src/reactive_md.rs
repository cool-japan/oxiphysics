// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Reactive MD extensions: ReaxFF bond order, AIREBO, FIRE minimizer, NEB, KMC.
//!
//! This module provides infrastructure for reactive molecular dynamics where
//! bonds can form and break dynamically during simulation.

use crate::MdError;

/// Boltzmann constant in eV/K.
pub const KB_EV: f64 = 8.617333262e-5;

/// Coulomb constant `k_e = e² / (4π ε₀)` in eV·Å/e² (CODATA 2018).
///
/// This is the prefactor that turns the charge product of two atoms (in units
/// of the elementary charge) divided by their distance in Å into an energy in
/// eV: `E = k_e q_i q_j / r_ij`.  Every charge-equilibration routine of the
/// reactive modules shares this constant.
pub const KE_EV_ANG: f64 = 14.399_645_478_4;

// ---------------------------------------------------------------------------
// ChemicalBond
// ---------------------------------------------------------------------------

/// A dynamic chemical bond that can form and break during simulation.
#[derive(Clone, Debug, PartialEq)]
pub struct ChemicalBond {
    /// Index of atom i.
    pub atom_i: usize,
    /// Index of atom j.
    pub atom_j: usize,
    /// Current bond order (0 = broken, 1 = single, 2 = double, 3 = triple).
    pub bond_order: f64,
    /// Bond type string (e.g., "C-C", "C=O").
    pub bond_type: String,
    /// Distance threshold for bond formation.
    pub form_cutoff: f64,
    /// Distance threshold for bond breaking.
    pub break_cutoff: f64,
}

impl ChemicalBond {
    /// Create a new chemical bond.
    pub fn new(
        atom_i: usize,
        atom_j: usize,
        bond_type: &str,
        form_cutoff: f64,
        break_cutoff: f64,
    ) -> Self {
        Self {
            atom_i,
            atom_j,
            bond_order: 0.0,
            bond_type: bond_type.to_string(),
            form_cutoff,
            break_cutoff,
        }
    }

    /// Update bond order based on current interatomic distance.
    pub fn update(&mut self, distance: f64) {
        if distance < self.form_cutoff {
            // Smooth bond formation using Fermi-Dirac-like function
            let r_mid = (self.form_cutoff + self.break_cutoff) / 2.0;
            let sigma = (self.break_cutoff - self.form_cutoff) / 6.0;
            self.bond_order = 1.0 / (1.0 + ((distance - r_mid) / sigma).exp());
        } else if distance > self.break_cutoff {
            self.bond_order = 0.0;
        }
    }

    /// Returns true if the bond is considered formed (bond_order > 0.5).
    pub fn is_formed(&self) -> bool {
        self.bond_order > 0.5
    }
}

// ---------------------------------------------------------------------------
// BondOrderTable
// ---------------------------------------------------------------------------

/// Lookup table for bond order parameters for common atom pairs.
#[derive(Clone, Debug)]
pub struct BondOrderTable {
    /// Atom pair labels (e.g., \["C", "C"\]).
    pub pair_labels: Vec<[String; 2]>,
    /// Equilibrium bond lengths for each order \[single, double, triple\] (Angstrom).
    pub eq_lengths: Vec<[f64; 3]>,
    /// Formation cutoffs for each pair (Angstrom).
    pub form_cutoffs: Vec<f64>,
    /// Breaking cutoffs for each pair (Angstrom).
    pub break_cutoffs: Vec<f64>,
    /// Bond energy for each order \[single, double, triple\] (eV).
    pub bond_energies: Vec<[f64; 3]>,
}

impl BondOrderTable {
    /// Create a new empty bond order table.
    pub fn new() -> Self {
        Self {
            pair_labels: Vec::new(),
            eq_lengths: Vec::new(),
            form_cutoffs: Vec::new(),
            break_cutoffs: Vec::new(),
            bond_energies: Vec::new(),
        }
    }

    /// Add a common C-C bond entry.
    pub fn add_cc_entry(&mut self) {
        self.pair_labels.push(["C".to_string(), "C".to_string()]);
        self.eq_lengths.push([1.54, 1.34, 1.20]);
        self.form_cutoffs.push(1.80);
        self.break_cutoffs.push(2.20);
        self.bond_energies.push([3.6, 6.3, 8.7]);
    }

    /// Add a C-H bond entry.
    pub fn add_ch_entry(&mut self) {
        self.pair_labels.push(["C".to_string(), "H".to_string()]);
        self.eq_lengths.push([1.09, 1.09, 1.09]);
        self.form_cutoffs.push(1.40);
        self.break_cutoffs.push(1.80);
        self.bond_energies.push([4.3, 4.3, 4.3]);
    }

    /// Add a C-O bond entry.
    pub fn add_co_entry(&mut self) {
        self.pair_labels.push(["C".to_string(), "O".to_string()]);
        self.eq_lengths.push([1.43, 1.23, 1.16]);
        self.form_cutoffs.push(1.70);
        self.break_cutoffs.push(2.10);
        self.bond_energies.push([3.6, 7.3, 10.6]);
    }

    /// Lookup the formation cutoff for a pair, returns None if not found.
    pub fn get_form_cutoff(&self, elem_i: &str, elem_j: &str) -> Option<f64> {
        for (idx, pair) in self.pair_labels.iter().enumerate() {
            if (pair[0] == elem_i && pair[1] == elem_j) || (pair[0] == elem_j && pair[1] == elem_i)
            {
                return Some(self.form_cutoffs[idx]);
            }
        }
        None
    }

    /// Lookup the breaking cutoff for a pair.
    pub fn get_break_cutoff(&self, elem_i: &str, elem_j: &str) -> Option<f64> {
        for (idx, pair) in self.pair_labels.iter().enumerate() {
            if (pair[0] == elem_i && pair[1] == elem_j) || (pair[0] == elem_j && pair[1] == elem_i)
            {
                return Some(self.break_cutoffs[idx]);
            }
        }
        None
    }
}

impl Default for BondOrderTable {
    fn default() -> Self {
        let mut t = Self::new();
        t.add_cc_entry();
        t.add_ch_entry();
        t.add_co_entry();
        t
    }
}

// ---------------------------------------------------------------------------
// ReaxFFBondOrder
// ---------------------------------------------------------------------------

/// ReaxFF bond order computation: σ/π/ππ contributions with
/// over/under-coordination correction.
#[derive(Clone, Debug)]
pub struct ReaxFFBondOrder {
    /// ReaxFF p_bo1 parameter.
    pub p_bo1: f64,
    /// ReaxFF p_bo2 parameter.
    pub p_bo2: f64,
    /// ReaxFF p_bo3 parameter.
    pub p_bo3: f64,
    /// ReaxFF p_bo4 parameter.
    pub p_bo4: f64,
    /// ReaxFF p_bo5 parameter.
    pub p_bo5: f64,
    /// ReaxFF p_bo6 parameter.
    pub p_bo6: f64,
    /// Equilibrium sigma bond length r_o_sigma (Angstrom).
    pub r_o_sigma: f64,
    /// Equilibrium pi bond length r_o_pi (Angstrom).
    pub r_o_pi: f64,
    /// Equilibrium double-pi bond length r_o_pipi (Angstrom).
    pub r_o_pipi: f64,
}

impl ReaxFFBondOrder {
    /// Create ReaxFF parameters for C-C bonds (typical values).
    pub fn new_cc() -> Self {
        Self {
            p_bo1: -0.097,
            p_bo2: 6.375,
            p_bo3: -0.198,
            p_bo4: 9.000,
            p_bo5: -0.288,
            p_bo6: 20.00,
            r_o_sigma: 1.3850,
            r_o_pi: 1.1000,
            r_o_pipi: 1.0300,
        }
    }

    /// Compute sigma bond order contribution.
    pub fn sigma_bond_order(&self, r: f64) -> f64 {
        (self.p_bo1 * (r / self.r_o_sigma).powf(self.p_bo2)).exp()
    }

    /// Compute pi bond order contribution.
    pub fn pi_bond_order(&self, r: f64) -> f64 {
        (self.p_bo3 * (r / self.r_o_pi).powf(self.p_bo4)).exp()
    }

    /// Compute double-pi (pipi) bond order contribution.
    pub fn pipi_bond_order(&self, r: f64) -> f64 {
        (self.p_bo5 * (r / self.r_o_pipi).powf(self.p_bo6)).exp()
    }

    /// Total bond order BO' = BO_sigma + BO_pi + BO_pipi (uncorrected).
    pub fn total_bond_order_uncorrected(&self, r: f64) -> f64 {
        self.sigma_bond_order(r) + self.pi_bond_order(r) + self.pipi_bond_order(r)
    }

    /// Over-coordination correction factor (simplified).
    ///
    /// `delta_i` = valence deviation for atom i (sum of bond orders - valence).
    pub fn over_coord_correction(&self, bo_prime: f64, delta_i: f64) -> f64 {
        let f1 = 0.5;
        let corrected = bo_prime - f1 * (bo_prime * delta_i.abs()).sqrt();
        corrected.max(0.0)
    }

    /// Under-coordination correction (simplified linear model).
    pub fn under_coord_correction(&self, bo_prime: f64, delta_i: f64) -> f64 {
        if delta_i < 0.0 {
            bo_prime * (1.0 + 0.1 * delta_i.abs())
        } else {
            bo_prime
        }
    }
}

// ---------------------------------------------------------------------------
// AireboPotential
// ---------------------------------------------------------------------------

/// AIREBO (Adaptive Intermolecular Reactive Empirical Bond Order) potential
/// for carbon and hydrogen systems.
#[derive(Clone, Debug)]
pub struct AireboPotential {
    /// REBO cutoff inner radius (Angstrom).
    pub r_inner: f64,
    /// REBO cutoff outer radius (Angstrom).
    pub r_outer: f64,
    /// LJ epsilon for C-C (eV).
    pub lj_epsilon_cc: f64,
    /// LJ sigma for C-C (Angstrom).
    pub lj_sigma_cc: f64,
    /// LJ epsilon for C-H (eV).
    pub lj_epsilon_ch: f64,
    /// LJ sigma for C-H (Angstrom).
    pub lj_sigma_ch: f64,
    /// Torsional potential strength (eV).
    pub torsion_k: f64,
}

impl AireboPotential {
    /// Create AIREBO with standard C/H parameters.
    pub fn new() -> Self {
        Self {
            r_inner: 1.7,
            r_outer: 2.0,
            lj_epsilon_cc: 0.00284,
            lj_sigma_cc: 3.40,
            lj_epsilon_ch: 0.001497,
            lj_sigma_ch: 2.94,
            torsion_k: 0.003079,
        }
    }

    /// Switching function S(r): 0 at r_outer, 1 at r_inner.
    pub fn switching_function(&self, r: f64) -> f64 {
        if r <= self.r_inner {
            1.0
        } else if r >= self.r_outer {
            0.0
        } else {
            let t = (r - self.r_inner) / (self.r_outer - self.r_inner);
            1.0 - t * t * (3.0 - 2.0 * t)
        }
    }

    /// 12-6 Lennard-Jones energy between two C atoms.
    pub fn lj_energy_cc(&self, r: f64) -> f64 {
        let sr = self.lj_sigma_cc / r;
        4.0 * self.lj_epsilon_cc * (sr.powi(12) - sr.powi(6))
    }

    /// 12-6 Lennard-Jones energy for C-H pair.
    pub fn lj_energy_ch(&self, r: f64) -> f64 {
        let sr = self.lj_sigma_ch / r;
        4.0 * self.lj_epsilon_ch * (sr.powi(12) - sr.powi(6))
    }

    /// Torsional potential energy for dihedral angle phi (radians).
    pub fn torsion_energy(&self, phi: f64) -> f64 {
        self.torsion_k
            * (1.0 - (256.0 / 405.0) * (phi / 2.0).cos().powi(10)
                + (1.0 / 10.0) * (phi / 2.0).cos().powi(2))
    }

    /// REBO-style Tersoff bond-order function (simplified Brenner).
    ///
    /// Returns the bond-order factor for atom pair i-j given coordination `n_coord`.
    pub fn rebo_bond_order(&self, n_coord: f64) -> f64 {
        let delta = 0.5;
        let n = 3.7;
        (1.0 + (n_coord / n).powf(2.0 * delta)).powf(-1.0 / (2.0 * delta))
    }
}

impl Default for AireboPotential {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// ReactiveMdConfig
// ---------------------------------------------------------------------------

/// Configuration for reactive MD simulations.
#[derive(Clone, Debug)]
pub struct ReactiveMdConfig {
    /// Global cutoff for bonded interactions (Angstrom).
    pub bond_cutoff: f64,
    /// Global cutoff for non-bonded interactions (Angstrom).
    pub nonbond_cutoff: f64,
    /// Number of charge equilibration iterations per MD step.
    pub charge_eq_steps: usize,
    /// Charge equilibration convergence tolerance.
    pub charge_eq_tol: f64,
    /// Charge damping factor for EEM.
    pub charge_damping: f64,
    /// Maximum bond order included in energy.
    pub max_bond_order: f64,
}

impl ReactiveMdConfig {
    /// Create default reactive MD configuration.
    pub fn new() -> Self {
        Self {
            bond_cutoff: 2.5,
            nonbond_cutoff: 10.0,
            charge_eq_steps: 300,
            charge_eq_tol: 1e-6,
            charge_damping: 0.7,
            max_bond_order: 3.0,
        }
    }
}

impl Default for ReactiveMdConfig {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// EemCharges
// ---------------------------------------------------------------------------

/// Electronegativity Equalization Method (EEM) for charge assignment.
///
/// The charges minimise the electrostatic energy
///
/// ```text
/// E(q) = Σ_i (χ_i q_i + η_i q_i²) + Σ_{i<j} J_ij q_i q_j      subject to Σ_i q_i = Q,
/// ```
///
/// i.e. they solve the bordered linear system, whose Lagrange multiplier `μ` is
/// minus the common, equalised electronegativity (`χ_i + Σ_j H_ij q_j = −μ` for
/// every atom)
///
/// ```text
/// [ 2 η_i δ_ij + J_ij (i ≠ j)    1 ] [ q ]   [ −χ ]
/// [ 1ᵀ                           0 ] [ μ ] = [  Q ].
/// ```
///
/// Electronegativities `χ` are in eV, hardnesses `η` in eV/e², distances in Å and
/// charges in units of the elementary charge, so the off-diagonal Coulomb kernel
/// carries the constant [`KE_EV_ANG`] (eV·Å/e²):
///
/// * bare Coulomb (default): `J_ij = k_e / r_ij`;
/// * ReaxFF shielded Coulomb (after [`EemCharges::with_shielding`]):
///   `J_ij = k_e / (r_ij³ + γ_ij⁻³)^{1/3}` with `γ_ij = √(γ_i γ_j)`, which stays
///   finite (`k_e γ_ij`) when two atoms approach each other.
///
/// The diagonal is `2 η_i` (`η` as in the energy expression above, half the
/// self-Coulomb `J_ii`), whereas [`crate::reactive::eem_charges`] writes the
/// energy as `Σ(χ q + ½ η q²)` and takes `η_i` itself on the diagonal: pass
/// `2 η_i` there to obtain the same charges.
#[derive(Clone, Debug)]
pub struct EemCharges {
    /// Electronegativities for each atom (eV).
    pub electronegativities: Vec<f64>,
    /// Chemical hardnesses for each atom (eV/e^2).
    pub hardnesses: Vec<f64>,
    /// Total charge constraint.
    pub total_charge: f64,
    /// Optional ReaxFF shielding parameters `γ_i` per atom (1/Å); `None` selects
    /// the bare Coulomb kernel `k_e / r`.
    pub shielding: Option<Vec<f64>>,
}

impl EemCharges {
    /// Create a new EEM charge solver with the bare Coulomb kernel.
    pub fn new(electronegativities: Vec<f64>, hardnesses: Vec<f64>, total_charge: f64) -> Self {
        Self {
            electronegativities,
            hardnesses,
            total_charge,
            shielding: None,
        }
    }

    /// Switch to the ReaxFF shielded Coulomb kernel with per-atom parameters
    /// `gammas` (1/Å, one per atom, non-negative and finite).
    pub fn with_shielding(mut self, gammas: Vec<f64>) -> Self {
        self.shielding = Some(gammas);
        self
    }

    /// Number of atoms.
    pub fn n_atoms(&self) -> usize {
        self.electronegativities.len()
    }

    /// Check that all per-atom arrays and the positions describe the same atoms.
    fn validate(&self, positions: &[[f64; 3]]) -> Result<(), MdError> {
        let n = self.n_atoms();
        if self.hardnesses.len() != n || positions.len() != n {
            return Err(MdError::General(format!(
                "EEM needs one entry per atom: {n} electronegativities, {} hardnesses, \
                 {} positions",
                self.hardnesses.len(),
                positions.len()
            )));
        }
        if let Some(gammas) = &self.shielding {
            if gammas.len() != n {
                return Err(MdError::General(format!(
                    "EEM shielding needs one parameter per atom: {n} atoms, {} parameters",
                    gammas.len()
                )));
            }
            if gammas.iter().any(|g| !(g.is_finite() && *g >= 0.0)) {
                return Err(MdError::General(
                    "EEM shielding parameters must be finite and non-negative".into(),
                ));
            }
        }
        Ok(())
    }

    /// Coulomb kernel `J_ij` (eV/e²) for atoms `i` and `j` at distance `distance` (Å).
    fn coupling(&self, i: usize, j: usize, distance: f64) -> f64 {
        let r = distance.max(1e-6);
        let gammas = self
            .shielding
            .as_ref()
            .and_then(|gammas| gammas.get(i).zip(gammas.get(j)));
        match gammas {
            Some((&gamma_i, &gamma_j)) => {
                let gamma_ij = (gamma_i * gamma_j).sqrt();
                KE_EV_ANG / (r.powi(3) + gamma_ij.powi(-3)).cbrt()
            }
            None => KE_EV_ANG / r,
        }
    }

    /// Build the EEM matrix (n+1) x (n+1) including charge constraint row.
    ///
    /// The diagonal holds `2 η_i`, the off-diagonal the Coulomb kernel
    /// `J_ij` in eV/e² (see the type-level documentation), and the last
    /// row/column the Lagrange constraint `Σ q = Q`.
    ///
    /// # Errors
    /// [`MdError::General`] if the per-atom arrays, the shielding parameters and
    /// `positions` do not all describe the same number of atoms, or a shielding
    /// parameter is negative or not finite.
    pub fn build_matrix(&self, positions: &[[f64; 3]]) -> Result<Vec<Vec<f64>>, MdError> {
        self.validate(positions)?;
        let n = self.n_atoms();
        let mut mat = vec![vec![0.0f64; n + 1]; n + 1];

        // Diagonal: 2 * hardness
        for (i, row) in mat.iter_mut().enumerate().take(n) {
            row[i] = 2.0 * self.hardnesses[i];
        }

        // Off-diagonal: Coulomb kernel in eV/e^2
        for i in 0..n {
            for j in (i + 1)..n {
                let dx = positions[i][0] - positions[j][0];
                let dy = positions[i][1] - positions[j][1];
                let dz = positions[i][2] - positions[j][2];
                let r = (dx * dx + dy * dy + dz * dz).sqrt();
                let coupling = self.coupling(i, j, r);
                mat[i][j] = coupling;
                mat[j][i] = coupling;
            }
        }

        // Lagrange constraint for total charge
        for row in mat.iter_mut().take(n) {
            row[n] = 1.0;
        }
        for v in mat[n].iter_mut().take(n) {
            *v = 1.0;
        }
        mat[n][n] = 0.0;

        Ok(mat)
    }

    /// Solve for the charges by LU decomposition with partial pivoting.
    ///
    /// Returns per-atom charges satisfying `Σ q = total_charge`.
    ///
    /// # Errors
    /// [`MdError::General`] for inconsistent input (see
    /// [`EemCharges::build_matrix`]), a singular EEM matrix, or a non-finite
    /// solution.
    pub fn solve(&self, positions: &[[f64; 3]]) -> Result<Vec<f64>, MdError> {
        self.validate(positions)?;
        let n = self.n_atoms();
        if n == 0 {
            return Ok(Vec::new());
        }
        let matrix = self.build_matrix(positions)?;
        let mut rhs: Vec<f64> = self.electronegativities.iter().map(|chi| -chi).collect();
        rhs.push(self.total_charge);
        let mut solution =
            crate::reactive::functions::solve_bordered_system(&matrix, &rhs, 1e-10, 4)?;
        solution.truncate(n);
        Ok(solution)
    }
}

// ---------------------------------------------------------------------------
// FireMinimizer
// ---------------------------------------------------------------------------

/// FIRE (Fast Inertial Relaxation Engine) energy minimizer.
///
/// Reference: Bitzek et al., Phys. Rev. Lett. 97, 170201 (2006).
#[derive(Clone, Debug)]
pub struct FireMinimizer {
    /// Time step (fs).
    pub dt: f64,
    /// Maximum allowed time step (fs).
    pub dt_max: f64,
    /// Time step increase factor.
    pub f_inc: f64,
    /// Time step decrease factor.
    pub f_dec: f64,
    /// Velocity mixing parameter alpha_start.
    pub alpha_start: f64,
    /// Alpha decrease factor.
    pub f_alpha: f64,
    /// Number of steps before allowing dt increase.
    pub n_min: usize,
    /// Force convergence criterion (eV/Angstrom).
    pub force_tol: f64,
    /// Energy convergence criterion (eV).
    pub energy_tol: f64,
}

impl FireMinimizer {
    /// Create FIRE minimizer with default parameters.
    pub fn new() -> Self {
        Self {
            dt: 1.0,
            dt_max: 10.0,
            f_inc: 1.1,
            f_dec: 0.5,
            alpha_start: 0.1,
            f_alpha: 0.99,
            n_min: 5,
            force_tol: 1e-4,
            energy_tol: 1e-7,
        }
    }

    /// FIRE velocity update: mix velocity toward force direction.
    ///
    /// v_new = (1 - alpha) * v + alpha * F/|F| * |v|
    pub fn mix_velocity(velocities: &mut [[f64; 3]], forces: &[[f64; 3]], alpha: f64) {
        let n = velocities.len();
        // Compute |v| and |F|
        let mut v_norm2 = 0.0f64;
        let mut f_norm2 = 0.0f64;
        for i in 0..n {
            for d in 0..3 {
                v_norm2 += velocities[i][d] * velocities[i][d];
                f_norm2 += forces[i][d] * forces[i][d];
            }
        }
        let v_norm = v_norm2.sqrt();
        let f_norm = f_norm2.sqrt().max(1e-30);

        // v = (1-alpha)*v + alpha*(F/|F|)*|v|
        for i in 0..n {
            for d in 0..3 {
                velocities[i][d] =
                    (1.0 - alpha) * velocities[i][d] + alpha * (forces[i][d] / f_norm) * v_norm;
            }
        }
    }

    /// Check if converged: max force component < force_tol.
    pub fn is_converged(&self, forces: &[[f64; 3]]) -> bool {
        for f in forces {
            for &c in f {
                if c.abs() > self.force_tol {
                    return false;
                }
            }
        }
        true
    }

    /// Compute power P = F · v.
    pub fn compute_power(forces: &[[f64; 3]], velocities: &[[f64; 3]]) -> f64 {
        let mut p = 0.0;
        for (f, v) in forces.iter().zip(velocities.iter()) {
            for d in 0..3 {
                p += f[d] * v[d];
            }
        }
        p
    }

    /// Perform one FIRE step given forces, returns updated (dt, alpha, n_steps_positive).
    pub fn step(
        &self,
        positions: &mut [[f64; 3]],
        velocities: &mut [[f64; 3]],
        forces: &[[f64; 3]],
        masses: &[f64],
        dt: f64,
        alpha: f64,
        n_pos: usize,
    ) -> (f64, f64, usize) {
        let n = positions.len();
        let p = Self::compute_power(forces, velocities);

        // Velocity Verlet step
        for i in 0..n {
            let m = masses[i];
            for d in 0..3 {
                velocities[i][d] += 0.5 * forces[i][d] / m * dt;
                positions[i][d] += velocities[i][d] * dt;
            }
        }

        let (new_dt, new_alpha, new_n_pos) = if p > 0.0 {
            let np = n_pos + 1;
            let mut ndt = dt;
            let mut na = alpha;
            if np > self.n_min {
                ndt = (dt * self.f_inc).min(self.dt_max);
                na = alpha * self.f_alpha;
            }
            (ndt, na, np)
        } else {
            // Reset velocities to zero
            for vel in velocities.iter_mut().take(n) {
                *vel = [0.0; 3];
            }
            (dt * self.f_dec, self.alpha_start, 0)
        };

        // Mix velocity
        let mut vel_copy = velocities.to_vec();
        Self::mix_velocity(&mut vel_copy, forces, new_alpha);
        velocities.copy_from_slice(&vel_copy);

        (new_dt, new_alpha, new_n_pos)
    }
}

impl Default for FireMinimizer {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// NudgedElasticBand
// ---------------------------------------------------------------------------

/// Nudged Elastic Band (NEB) method for minimum energy path computation.
///
/// Reference: Henkelman & Jonsson, J. Chem. Phys. 113, 9901 (2000).
#[derive(Clone, Debug)]
pub struct NudgedElasticBand {
    /// Number of images (including endpoints).
    pub n_images: usize,
    /// Spring constant k (eV/Angstrom^2).
    pub spring_k: f64,
    /// Whether to use climbing image NEB.
    pub climbing_image: bool,
    /// Index of image at highest energy (for climbing image).
    pub climbing_idx: Option<usize>,
}

impl NudgedElasticBand {
    /// Create NEB with given images and spring constant.
    pub fn new(n_images: usize, spring_k: f64, climbing_image: bool) -> Self {
        Self {
            n_images,
            spring_k,
            climbing_image,
            climbing_idx: None,
        }
    }

    /// Compute the tangent vector for image `i` given neighboring image positions and energies.
    ///
    /// Uses the improved tangent estimate of Henkelman et al.
    pub fn tangent(
        &self,
        pos_prev: &[[f64; 3]],
        pos_curr: &[[f64; 3]],
        pos_next: &[[f64; 3]],
        e_prev: f64,
        e_curr: f64,
        e_next: f64,
    ) -> Vec<[f64; 3]> {
        let n = pos_curr.len();
        let mut tau = vec![[0.0f64; 3]; n];

        let tau_plus: Vec<[f64; 3]> = (0..n)
            .map(|k| {
                [
                    pos_next[k][0] - pos_curr[k][0],
                    pos_next[k][1] - pos_curr[k][1],
                    pos_next[k][2] - pos_curr[k][2],
                ]
            })
            .collect();
        let tau_minus: Vec<[f64; 3]> = (0..n)
            .map(|k| {
                [
                    pos_curr[k][0] - pos_prev[k][0],
                    pos_curr[k][1] - pos_prev[k][1],
                    pos_curr[k][2] - pos_prev[k][2],
                ]
            })
            .collect();

        if e_next > e_curr && e_curr > e_prev {
            tau = tau_plus;
        } else if e_next < e_curr && e_curr < e_prev {
            tau = tau_minus;
        } else {
            let de_max = (e_next - e_curr).abs().max((e_curr - e_prev).abs());
            let de_min = (e_next - e_curr).abs().min((e_curr - e_prev).abs());
            if e_next > e_prev {
                for k in 0..n {
                    for d in 0..3 {
                        tau[k][d] = de_max * tau_plus[k][d] + de_min * tau_minus[k][d];
                    }
                }
            } else {
                for k in 0..n {
                    for d in 0..3 {
                        tau[k][d] = de_min * tau_plus[k][d] + de_max * tau_minus[k][d];
                    }
                }
            }
        }

        // Normalize
        let norm: f64 = tau
            .iter()
            .map(|t| t[0] * t[0] + t[1] * t[1] + t[2] * t[2])
            .sum::<f64>()
            .sqrt()
            .max(1e-30);
        for t in &mut tau {
            for td in t.iter_mut() {
                *td /= norm;
            }
        }
        tau
    }

    /// Compute NEB spring force for image `i`.
    pub fn spring_force(
        &self,
        pos_prev: &[[f64; 3]],
        pos_curr: &[[f64; 3]],
        pos_next: &[[f64; 3]],
        tau: &[[f64; 3]],
    ) -> Vec<[f64; 3]> {
        let n = pos_curr.len();

        // |R_{i+1} - R_i| - |R_i - R_{i-1}|
        let d_plus: f64 = (0..n)
            .map(|k| {
                let dx = pos_next[k][0] - pos_curr[k][0];
                let dy = pos_next[k][1] - pos_curr[k][1];
                let dz = pos_next[k][2] - pos_curr[k][2];
                dx * dx + dy * dy + dz * dz
            })
            .sum::<f64>()
            .sqrt();
        let d_minus: f64 = (0..n)
            .map(|k| {
                let dx = pos_curr[k][0] - pos_prev[k][0];
                let dy = pos_curr[k][1] - pos_prev[k][1];
                let dz = pos_curr[k][2] - pos_prev[k][2];
                dx * dx + dy * dy + dz * dz
            })
            .sum::<f64>()
            .sqrt();

        let scalar = self.spring_k * (d_plus - d_minus);
        (0..n)
            .map(|k| [scalar * tau[k][0], scalar * tau[k][1], scalar * tau[k][2]])
            .collect()
    }

    /// Set the climbing image to the highest-energy image (excluding endpoints).
    pub fn set_climbing_image(&mut self, energies: &[f64]) {
        if !self.climbing_image || energies.len() < 3 {
            return;
        }
        let max_idx = energies[1..energies.len() - 1]
            .iter()
            .enumerate()
            .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal))
            .map(|(i, _)| i + 1)
            .unwrap_or(1);
        self.climbing_idx = Some(max_idx);
    }

    /// Compute maximum force on all images (convergence criterion).
    pub fn max_force(forces_all: &[Vec<[f64; 3]>]) -> f64 {
        forces_all
            .iter()
            .flat_map(|img| img.iter())
            .flat_map(|f| f.iter())
            .map(|c| c.abs())
            .fold(0.0f64, f64::max)
    }
}

// ---------------------------------------------------------------------------
// TransitionStateSearch (Dimer Method)
// ---------------------------------------------------------------------------

mod dimer;
pub use dimer::{DimerRotation, DimerSearchResult, TransitionStateSearch};

// ---------------------------------------------------------------------------
// KineticMonteCarlo
// ---------------------------------------------------------------------------

/// Kinetic Monte Carlo (KMC) for rare event dynamics.
///
/// Uses the residence time algorithm (Bortz-Kalos-Lebowitz).
#[derive(Clone, Debug)]
pub struct KineticMonteCarlo {
    /// List of event rates (s^-1).
    pub rates: Vec<f64>,
    /// Names of events (for logging).
    pub event_names: Vec<String>,
    /// Current simulation time (s).
    pub time: f64,
    /// Total number of KMC steps performed.
    pub n_steps: u64,
    /// Temperature (K) for Arrhenius rate computation.
    pub temperature: f64,
}

impl KineticMonteCarlo {
    /// Create a new KMC engine.
    pub fn new(temperature: f64) -> Self {
        Self {
            rates: Vec::new(),
            event_names: Vec::new(),
            time: 0.0,
            n_steps: 0,
            temperature,
        }
    }

    /// Add an event with given Arrhenius parameters.
    ///
    /// rate = nu_0 * exp(-E_a / (k_B * T))
    pub fn add_event(&mut self, name: &str, prefactor: f64, activation_energy_ev: f64) {
        let rate = prefactor * (-activation_energy_ev / (KB_EV * self.temperature)).exp();
        self.rates.push(rate);
        self.event_names.push(name.to_string());
    }

    /// Add an event with a directly specified rate.
    pub fn add_event_rate(&mut self, name: &str, rate: f64) {
        self.rates.push(rate);
        self.event_names.push(name.to_string());
    }

    /// Total rate (sum of all event rates).
    pub fn total_rate(&self) -> f64 {
        self.rates.iter().sum()
    }

    /// Select the next event using a uniform random number `u1` in \[0,1).
    ///
    /// Returns the event index.
    pub fn select_event(&self, u1: f64) -> usize {
        let r_total = self.total_rate();
        let threshold = u1 * r_total;
        let mut cumulative = 0.0;
        for (idx, &rate) in self.rates.iter().enumerate() {
            cumulative += rate;
            if cumulative >= threshold {
                return idx;
            }
        }
        self.rates.len().saturating_sub(1)
    }

    /// Advance time using the residence time algorithm.
    ///
    /// dt = -ln(u2) / R_total
    pub fn advance_time(&mut self, u2: f64) {
        let r_total = self.total_rate();
        if r_total > 0.0 {
            self.time += -(u2.max(1e-300)).ln() / r_total;
        }
        self.n_steps += 1;
    }

    /// Execute one KMC step: select event, advance time.
    ///
    /// Returns the selected event index.
    pub fn step(&mut self, u1: f64, u2: f64) -> usize {
        let event = self.select_event(u1);
        self.advance_time(u2);
        event
    }

    /// Compute Arrhenius rate for given E_a (eV).
    pub fn arrhenius_rate(&self, prefactor: f64, e_activation: f64) -> f64 {
        prefactor * (-e_activation / (KB_EV * self.temperature)).exp()
    }

    /// Update a rate for an existing event.
    pub fn update_rate(&mut self, event_idx: usize, new_rate: f64) {
        if event_idx < self.rates.len() {
            self.rates[event_idx] = new_rate;
        }
    }
}

// ---------------------------------------------------------------------------
// Helper functions
// ---------------------------------------------------------------------------

/// Compute Lennard-Jones energy between two atoms at distance r.
pub fn lj_energy(epsilon: f64, sigma: f64, r: f64) -> f64 {
    let sr = sigma / r.max(1e-10);
    4.0 * epsilon * (sr.powi(12) - sr.powi(6))
}

/// Compute Lennard-Jones force magnitude (positive = repulsive).
pub fn lj_force(epsilon: f64, sigma: f64, r: f64) -> f64 {
    let sr = sigma / r.max(1e-10);
    24.0 * epsilon / r * (2.0 * sr.powi(12) - sr.powi(6))
}

/// Compute Morse potential energy.
///
/// V(r) = D_e * \[1 - exp(-a*(r-r_e))\]^2 - D_e
pub fn morse_energy(d_e: f64, a: f64, r_e: f64, r: f64) -> f64 {
    let x = 1.0 - (-a * (r - r_e)).exp();
    d_e * x * x - d_e
}

/// Compute Morse force magnitude.
pub fn morse_force(d_e: f64, a: f64, r_e: f64, r: f64) -> f64 {
    let e = (-a * (r - r_e)).exp();
    -2.0 * d_e * a * e * (1.0 - e)
}

/// Compute harmonic spring energy: 0.5 * k * (r - r0)^2.
pub fn harmonic_energy(k: f64, r0: f64, r: f64) -> f64 {
    0.5 * k * (r - r0) * (r - r0)
}

/// Compute harmonic spring force: -k * (r - r0).
pub fn harmonic_force(k: f64, r0: f64, r: f64) -> f64 {
    -k * (r - r0)
}

/// Euclidean distance between two 3D points.
pub fn distance_3d(a: &[f64; 3], b: &[f64; 3]) -> f64 {
    let dx = a[0] - b[0];
    let dy = a[1] - b[1];
    let dz = a[2] - b[2];
    (dx * dx + dy * dy + dz * dz).sqrt()
}

/// Normalize a 3D vector.
pub fn normalize_3d(v: &[f64; 3]) -> [f64; 3] {
    let norm = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt().max(1e-30);
    [v[0] / norm, v[1] / norm, v[2] / norm]
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_chemical_bond_new() {
        let b = ChemicalBond::new(0, 1, "C-C", 1.8, 2.2);
        assert_eq!(b.atom_i, 0);
        assert_eq!(b.atom_j, 1);
        assert!(!b.is_formed());
    }

    #[test]
    fn test_chemical_bond_update_formed() {
        let mut b = ChemicalBond::new(0, 1, "C-C", 1.8, 2.2);
        b.update(1.5);
        assert!(b.bond_order > 0.5);
        assert!(b.is_formed());
    }

    #[test]
    fn test_chemical_bond_update_broken() {
        let mut b = ChemicalBond::new(0, 1, "C-C", 1.8, 2.2);
        b.update(3.0);
        assert_eq!(b.bond_order, 0.0);
        assert!(!b.is_formed());
    }

    #[test]
    fn test_bond_order_table_default() {
        let t = BondOrderTable::default();
        assert_eq!(t.pair_labels.len(), 3);
    }

    #[test]
    fn test_bond_order_table_lookup_cc() {
        let t = BondOrderTable::default();
        let fc = t.get_form_cutoff("C", "C");
        assert!(fc.is_some());
        assert!((fc.unwrap() - 1.80).abs() < 1e-10);
    }

    #[test]
    fn test_bond_order_table_lookup_hc() {
        let t = BondOrderTable::default();
        let fc = t.get_form_cutoff("H", "C");
        assert!(fc.is_some());
    }

    #[test]
    fn test_bond_order_table_lookup_missing() {
        let t = BondOrderTable::default();
        let fc = t.get_form_cutoff("N", "N");
        assert!(fc.is_none());
    }

    #[test]
    fn test_reaxff_sigma_bond_order() {
        let reax = ReaxFFBondOrder::new_cc();
        let bo_sigma = reax.sigma_bond_order(1.54);
        assert!(bo_sigma > 0.0 && bo_sigma < 2.0);
    }

    #[test]
    fn test_reaxff_pi_bond_order() {
        let reax = ReaxFFBondOrder::new_cc();
        let bo_pi = reax.pi_bond_order(1.34);
        assert!(bo_pi >= 0.0);
    }

    #[test]
    fn test_reaxff_total_bond_order() {
        let reax = ReaxFFBondOrder::new_cc();
        let bo = reax.total_bond_order_uncorrected(1.54);
        assert!(bo > 0.0);
        // All contributions positive
        assert!(bo >= reax.sigma_bond_order(1.54));
    }

    #[test]
    fn test_reaxff_over_coord_correction() {
        let reax = ReaxFFBondOrder::new_cc();
        let corrected = reax.over_coord_correction(1.0, 0.5);
        assert!((0.0..=1.0).contains(&corrected));
    }

    #[test]
    fn test_airebo_switching_function() {
        let airebo = AireboPotential::new();
        assert_eq!(airebo.switching_function(1.0), 1.0);
        assert_eq!(airebo.switching_function(3.0), 0.0);
        let mid = airebo.switching_function((airebo.r_inner + airebo.r_outer) / 2.0);
        assert!(mid > 0.0 && mid < 1.0);
    }

    #[test]
    fn test_airebo_lj_energy_cc() {
        let airebo = AireboPotential::new();
        let e_min = airebo.lj_energy_cc(airebo.lj_sigma_cc * 2.0f64.powf(1.0 / 6.0));
        assert!((e_min + airebo.lj_epsilon_cc).abs() < 1e-10);
    }

    #[test]
    fn test_airebo_rebo_bond_order() {
        let airebo = AireboPotential::new();
        let bo0 = airebo.rebo_bond_order(0.0);
        let bo3 = airebo.rebo_bond_order(3.0);
        assert!((bo0 - 1.0).abs() < 1e-10);
        assert!(bo3 < bo0);
    }

    #[test]
    fn test_reactive_md_config_default() {
        let cfg = ReactiveMdConfig::default();
        assert!(cfg.bond_cutoff > 0.0);
        assert!(cfg.charge_eq_steps > 0);
        assert!(cfg.charge_damping > 0.0);
    }

    #[test]
    fn test_eem_charges_build_matrix_size() {
        let chi = vec![6.96, 8.74];
        let eta = vec![12.1, 13.4];
        let eem = EemCharges::new(chi, eta, 0.0);
        let pos = vec![[0.0, 0.0, 0.0], [1.4, 0.0, 0.0]];
        let mat = eem.build_matrix(&pos).expect("two atoms, two positions");
        assert_eq!(mat.len(), 3);
        assert_eq!(mat[0].len(), 3);
    }

    #[test]
    fn test_eem_charges_solve_neutral() {
        let chi = vec![6.96, 8.74, 7.5];
        let eta = vec![12.1, 13.4, 11.0];
        let eem = EemCharges::new(chi, eta, 0.0);
        let pos = vec![[0.0, 0.0, 0.0], [1.4, 0.0, 0.0], [2.8, 0.0, 0.0]];
        let charges = eem.solve(&pos).expect("EEM matrix is regular");
        assert_eq!(charges.len(), 3);
        let total: f64 = charges.iter().sum();
        assert!(total.abs() < 1e-10, "total charge {total}");
        // The more electronegative second atom is negative, the first positive.
        assert!(charges[1] < 0.0 && charges[0] > 0.0, "{charges:?}");
    }

    #[test]
    fn test_eem_matrix_carries_the_coulomb_constant() {
        // Pre-fix the off-diagonal was 1/r (Å⁻¹) next to eV/e² hardnesses: charges were
        // wrong by about 14x.  The kernel is k_e / r in eV/e².
        let eem = EemCharges::new(vec![5.0, 7.0], vec![9.0, 11.0], 0.0);
        let r = 2.5;
        let mat = eem
            .build_matrix(&[[0.0, 0.0, 0.0], [r, 0.0, 0.0]])
            .expect("two atoms, two positions");
        assert!(
            (mat[0][1] - KE_EV_ANG / r).abs() < 1e-12,
            "J = {}",
            mat[0][1]
        );
        assert_eq!(mat[0][1], mat[1][0]);
        assert!((KE_EV_ANG - 14.399_645).abs() < 1e-6);
        // Diagonal 2 η, Lagrange border of ones with a zero corner.
        assert!((mat[0][0] - 18.0).abs() < 1e-14 && (mat[1][1] - 22.0).abs() < 1e-14);
        assert_eq!(
            (mat[0][2], mat[1][2], mat[2][0], mat[2][1], mat[2][2]),
            (1.0, 1.0, 1.0, 1.0, 0.0)
        );
    }

    #[test]
    fn test_eem_two_atom_charge_matches_the_closed_form() {
        // q_1 = (χ_2 − χ_1) / (2η_1 + 2η_2 − 2 k_e / r),  q_2 = −q_1  (neutral pair).
        let (chi, eta) = (vec![5.0, 7.0], vec![9.0, 11.0]);
        let eem = EemCharges::new(chi.clone(), eta.clone(), 0.0);
        for &r in &[3.0, 5.0, 8.0] {
            let charges = eem
                .solve(&[[0.0, 0.0, 0.0], [0.0, r, 0.0]])
                .expect("EEM matrix is regular");
            let expected = (chi[1] - chi[0]) / (2.0 * eta[0] + 2.0 * eta[1] - 2.0 * KE_EV_ANG / r);
            assert!(
                (charges[0] - expected).abs() < 1e-10,
                "r = {r}: {} vs {expected}",
                charges[0]
            );
            assert!((charges[1] + expected).abs() < 1e-10);
        }
    }

    #[test]
    fn test_eem_charges_equalise_the_electronegativity_and_sum_to_the_total() {
        let chi = vec![4.53, 6.93, 4.53, 5.2];
        let eta = vec![5.0, 8.0, 5.0, 6.5];
        let pos = vec![
            [-3.0, 0.0, 0.0],
            [0.0, 0.2, 0.0],
            [3.0, -0.1, 0.4],
            [0.5, 3.5, 1.0],
        ];
        for &total in &[0.0, 0.7, -2.0] {
            let eem = EemCharges::new(chi.clone(), eta.clone(), total);
            let charges = eem.solve(&pos).expect("EEM matrix is regular");
            let mat = eem.build_matrix(&pos).expect("consistent input");
            let effective: Vec<f64> = (0..4)
                .map(|i| chi[i] + (0..4).map(|j| mat[i][j] * charges[j]).sum::<f64>())
                .collect();
            let spread = effective.iter().cloned().fold(f64::NEG_INFINITY, f64::max)
                - effective.iter().cloned().fold(f64::INFINITY, f64::min);
            assert!(
                spread < 1e-9,
                "Q = {total}: electronegativity spread {spread}"
            );
            assert!((charges.iter().sum::<f64>() - total).abs() < 1e-10);
            assert!(
                charges.iter().any(|q| q.abs() > 1e-3),
                "charges collapsed: {charges:?}"
            );
        }
    }

    #[test]
    fn test_eem_matches_the_free_function_with_doubled_hardness() {
        // EemCharges puts 2 η on the diagonal (E = Σ χ q + η q²), crate::reactive::eem_charges
        // puts η (E = Σ χ q + ½ η q²): the two agree for η_free = 2 η.
        let chi = vec![6.96, 8.74, 7.5];
        let eta = vec![12.1, 13.4, 11.0];
        let pos = vec![[0.0, 0.0, 0.0], [1.4, 0.0, 0.0], [2.8, 0.0, 0.0]];
        let class = EemCharges::new(chi.clone(), eta.clone(), 1.0)
            .solve(&pos)
            .expect("EEM matrix is regular");
        let doubled: Vec<f64> = eta.iter().map(|e| 2.0 * e).collect();
        let free = crate::reactive::eem_charges(&pos, &chi, &doubled, 1.0, 1e-12, 10)
            .expect("EEM matrix is regular");
        for (a, b) in class.iter().zip(free.iter()) {
            assert!((a - b).abs() < 1e-9, "{a} vs {b}");
        }
        assert!((class.iter().sum::<f64>() - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_eem_shielded_kernel_follows_the_reaxff_form() {
        // J(r) = k_e / (r³ + γ_ij⁻³)^(1/3),  γ_ij = sqrt(γ_i γ_j).
        let eem =
            EemCharges::new(vec![5.0, 7.0], vec![9.0, 11.0], 0.0).with_shielding(vec![0.8, 1.2]);
        let gamma_ij = (0.8_f64 * 1.2).sqrt();
        for &r in &[0.5, 1.0, 2.0, 5.0] {
            let mat = eem
                .build_matrix(&[[0.0, 0.0, 0.0], [0.0, 0.0, r]])
                .expect("consistent input");
            let expected = KE_EV_ANG / (r.powi(3) + gamma_ij.powi(-3)).cbrt();
            assert!((mat[0][1] - expected).abs() < 1e-12, "r = {r}");
            assert!(
                mat[0][1] < KE_EV_ANG / r,
                "shielding must weaken the bare kernel"
            );
        }
        // Finite when the atoms coincide (bare Coulomb diverges): J(0) = k_e γ_ij.
        let contact = eem
            .build_matrix(&[[1.0, 1.0, 1.0], [1.0, 1.0, 1.0]])
            .expect("consistent input");
        assert!((contact[0][1] - KE_EV_ANG * gamma_ij).abs() < 1e-6 * KE_EV_ANG);
        // Far apart the shielded kernel is the bare k_e / r.
        let far = eem
            .build_matrix(&[[0.0, 0.0, 0.0], [60.0, 0.0, 0.0]])
            .expect("consistent input");
        assert!((far[0][1] * 60.0 / KE_EV_ANG - 1.0).abs() < 1e-5);
        // γ = 0 shields completely.
        let off =
            EemCharges::new(vec![5.0, 7.0], vec![9.0, 11.0], 0.0).with_shielding(vec![0.0, 1.0]);
        let mat = off
            .build_matrix(&[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
            .expect("consistent input");
        assert_eq!(mat[0][1], 0.0);
    }

    #[test]
    fn test_eem_shielding_regularises_close_atoms() {
        // At 0.4 Å the bare coupling k_e / r = 36 eV exceeds 2η, making the matrix
        // indefinite and the charge transfer unphysically large; the shielded kernel keeps
        // the neutral pair stable and the physical sign (less electronegative atom positive).
        let (chi, eta) = (vec![5.0, 7.0], vec![9.0, 11.0]);
        let pos = [[0.0, 0.0, 0.0], [0.4, 0.0, 0.0]];
        let bare = EemCharges::new(chi.clone(), eta.clone(), 0.0)
            .solve(&pos)
            .expect("EEM matrix is regular");
        let shielded = EemCharges::new(chi, eta, 0.0)
            .with_shielding(vec![0.9, 0.9])
            .solve(&pos)
            .expect("EEM matrix is regular");
        assert!(
            bare[0] < 0.0,
            "bare kernel gives the unphysical sign: {bare:?}"
        );
        assert!(shielded[0] > 0.0 && shielded[1] < 0.0, "{shielded:?}");
        let gamma_ij = 0.9_f64;
        let j = KE_EV_ANG / (0.4_f64.powi(3) + gamma_ij.powi(-3)).cbrt();
        let expected = 2.0 / (18.0 + 22.0 - 2.0 * j);
        assert!((shielded[0] - expected).abs() < 1e-10);
    }

    #[test]
    fn test_eem_handles_single_atom_and_empty_systems() {
        let single = EemCharges::new(vec![5.0], vec![9.0], 0.6)
            .solve(&[[0.0, 0.0, 0.0]])
            .expect("single atom");
        assert_eq!(single.len(), 1);
        assert!((single[0] - 0.6).abs() < 1e-12);
        let empty = EemCharges::new(Vec::new(), Vec::new(), 0.0)
            .solve(&[])
            .expect("empty system");
        assert!(empty.is_empty());
    }

    #[test]
    fn test_eem_rejects_inconsistent_input() {
        let eem = EemCharges::new(vec![5.0, 7.0], vec![9.0, 11.0], 0.0);
        assert!(eem.solve(&[[0.0; 3]]).is_err(), "too few positions");
        assert!(
            eem.build_matrix(&[[0.0; 3], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]])
                .is_err()
        );
        let short_eta = EemCharges::new(vec![5.0, 7.0], vec![9.0], 0.0);
        assert!(short_eta.solve(&[[0.0; 3], [1.0, 0.0, 0.0]]).is_err());
        let two_atoms = [[0.0; 3], [1.0, 0.0, 0.0]];
        assert!(
            eem.clone()
                .with_shielding(vec![1.0])
                .solve(&two_atoms)
                .is_err()
        );
        assert!(
            eem.clone()
                .with_shielding(vec![1.0, -1.0])
                .solve(&two_atoms)
                .is_err()
        );
        assert!(
            eem.with_shielding(vec![1.0, f64::NAN])
                .solve(&two_atoms)
                .is_err()
        );
        // Zero hardness with the coupling suppressed is a singular system, not silent zeros.
        let singular =
            EemCharges::new(vec![5.0, 7.0], vec![0.0, 0.0], 0.0).with_shielding(vec![0.0, 0.0]);
        assert!(singular.solve(&two_atoms).is_err());
    }

    #[test]
    fn test_fire_minimizer_default() {
        let fire = FireMinimizer::default();
        assert!(fire.force_tol > 0.0);
        assert!(fire.dt < fire.dt_max);
    }

    #[test]
    fn test_fire_mix_velocity_zero_force() {
        let mut vel = vec![[1.0, 0.0, 0.0]];
        let force = vec![[0.0, 0.0, 0.0]];
        FireMinimizer::mix_velocity(&mut vel, &force, 0.1);
        // With zero force, vel stays direction-mixed but might become 0
        // just check it doesn't panic
    }

    #[test]
    fn test_fire_power_positive() {
        let force = vec![[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let vel = vec![[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let p = FireMinimizer::compute_power(&force, &vel);
        assert!((p - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_fire_is_converged_true() {
        let fire = FireMinimizer::new();
        let forces = vec![[0.0, 0.0, 0.0]; 3];
        assert!(fire.is_converged(&forces));
    }

    #[test]
    fn test_fire_is_converged_false() {
        let fire = FireMinimizer::new();
        let forces = vec![[1.0, 0.0, 0.0]];
        assert!(!fire.is_converged(&forces));
    }

    #[test]
    fn test_neb_new() {
        let neb = NudgedElasticBand::new(5, 1.0, false);
        assert_eq!(neb.n_images, 5);
        assert!(neb.climbing_idx.is_none());
    }

    #[test]
    fn test_neb_set_climbing_image() {
        let mut neb = NudgedElasticBand::new(5, 1.0, true);
        let energies = vec![0.0, 0.5, 1.0, 0.8, 0.0];
        neb.set_climbing_image(&energies);
        assert_eq!(neb.climbing_idx, Some(2));
    }

    #[test]
    fn test_neb_tangent_uphill() {
        let neb = NudgedElasticBand::new(3, 1.0, false);
        let prev = vec![[0.0, 0.0, 0.0]];
        let curr = vec![[1.0, 0.0, 0.0]];
        let next = vec![[2.0, 0.0, 0.0]];
        let tau = neb.tangent(&prev, &curr, &next, 0.0, 0.5, 1.0);
        assert!((tau[0][0].abs() - 1.0).abs() < 0.01);
    }

    #[test]
    fn test_neb_spring_force() {
        let neb = NudgedElasticBand::new(3, 1.0, false);
        let prev = vec![[0.0, 0.0, 0.0]];
        let curr = vec![[1.0, 0.0, 0.0]];
        let next = vec![[2.0, 0.0, 0.0]];
        let tau = vec![[1.0, 0.0, 0.0]];
        let sf = neb.spring_force(&prev, &curr, &next, &tau);
        // Evenly spaced => zero spring force
        assert!(sf[0][0].abs() < 1e-10);
    }

    #[test]
    fn test_kmc_new() {
        let kmc = KineticMonteCarlo::new(300.0);
        assert_eq!(kmc.time, 0.0);
        assert_eq!(kmc.n_steps, 0);
    }

    #[test]
    fn test_kmc_add_event() {
        let mut kmc = KineticMonteCarlo::new(300.0);
        kmc.add_event("diffusion", 1e13, 0.5);
        assert_eq!(kmc.rates.len(), 1);
        assert!(kmc.rates[0] > 0.0);
    }

    #[test]
    fn test_kmc_total_rate() {
        let mut kmc = KineticMonteCarlo::new(300.0);
        kmc.add_event_rate("A", 1.0);
        kmc.add_event_rate("B", 2.0);
        assert!((kmc.total_rate() - 3.0).abs() < 1e-10);
    }

    #[test]
    fn test_kmc_select_event() {
        let mut kmc = KineticMonteCarlo::new(300.0);
        kmc.add_event_rate("A", 1.0);
        kmc.add_event_rate("B", 1.0);
        let e = kmc.select_event(0.1);
        assert_eq!(e, 0);
        let e2 = kmc.select_event(0.9);
        assert_eq!(e2, 1);
    }

    #[test]
    fn test_kmc_advance_time() {
        let mut kmc = KineticMonteCarlo::new(300.0);
        kmc.add_event_rate("A", 1e10);
        let t0 = kmc.time;
        kmc.advance_time(0.5);
        assert!(kmc.time > t0);
        assert_eq!(kmc.n_steps, 1);
    }

    #[test]
    fn test_kmc_step() {
        let mut kmc = KineticMonteCarlo::new(300.0);
        kmc.add_event_rate("A", 1.0);
        kmc.add_event_rate("B", 3.0);
        let e = kmc.step(0.9, 0.5);
        assert!(e < 2);
        assert_eq!(kmc.n_steps, 1);
    }

    #[test]
    fn test_kmc_arrhenius_rate() {
        let kmc = KineticMonteCarlo::new(300.0);
        let rate = kmc.arrhenius_rate(1e13, 0.0);
        assert!((rate - 1e13).abs() < 1.0);
    }

    #[test]
    fn test_lj_energy_minimum() {
        let e_min = lj_energy(1.0, 1.0, 2.0f64.powf(1.0 / 6.0));
        assert!((e_min + 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_lj_force_sign() {
        // At r < sigma, repulsive (positive force)
        let f = lj_force(1.0, 1.0, 0.9);
        assert!(f > 0.0);
        // At r > sigma, attractive
        let f2 = lj_force(1.0, 1.0, 1.2);
        assert!(f2 < 0.0);
    }

    #[test]
    fn test_morse_energy_minimum() {
        let e = morse_energy(1.0, 2.0, 1.5, 1.5);
        assert!((e - (-1.0)).abs() < 1e-10);
    }

    #[test]
    fn test_harmonic_energy_minimum() {
        let e = harmonic_energy(10.0, 1.5, 1.5);
        assert_eq!(e, 0.0);
    }

    #[test]
    fn test_harmonic_force_direction() {
        let f = harmonic_force(10.0, 1.5, 2.0);
        assert!(f < 0.0); // stretched: restoring toward equilibrium
    }

    #[test]
    fn test_distance_3d() {
        let a = [0.0, 0.0, 0.0];
        let b = [3.0, 4.0, 0.0];
        assert!((distance_3d(&a, &b) - 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_normalize_3d() {
        let v = [3.0, 4.0, 0.0];
        let n = normalize_3d(&v);
        let norm = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
        assert!((norm - 1.0).abs() < 1e-10);
    }
}
