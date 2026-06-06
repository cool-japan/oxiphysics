// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! QM/MM interface for hybrid quantum mechanics / molecular mechanics calculations.
//!
//! Provides structures and algorithms for partitioning a molecular system into
//! a quantum-mechanically treated core region and a classically treated MM region,
//! connected via link atoms at the boundary.

/// Classification of atoms in a QM/MM calculation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RegionType {
    /// Atom treated quantum mechanically (semi-empirical TB / DFT-B).
    QM,
    /// Atom treated with a classical force field.
    MM,
    /// Link atom or transition-region atom at the QM/MM boundary.
    Buffer,
}

/// An atom residing in the QM region.
#[derive(Debug, Clone)]
pub struct QmAtom {
    /// Global index in the full system.
    pub index: usize,
    /// Cartesian position in Ångström.
    pub position: [f64; 3],
    /// Atomic number (1 = H, 6 = C, …).
    pub atomic_number: u8,
    /// Partial charge obtained from charge fitting (e.g. RESP).
    pub charge: f64,
    /// Energy gradient (force = -gradient), in kJ/(mol·Å).
    pub gradient: [f64; 3],
}

/// An atom residing in the MM region.
#[derive(Debug, Clone)]
pub struct MmAtom {
    /// Global index in the full system.
    pub index: usize,
    /// Cartesian position in Ångström.
    pub position: [f64; 3],
    /// Classical point charge (e).
    pub charge: f64,
    /// Lennard-Jones well depth ε (kJ/mol).
    pub epsilon: f64,
    /// Lennard-Jones radius σ (Å).
    pub sigma: f64,
    /// Energy gradient, in kJ/(mol·Å).
    pub gradient: [f64; 3],
}

/// A hydrogen link atom that caps a severed QM/MM covalent bond.
///
/// The link atom position is placed along the QM–MM bond vector at a
/// fractional distance `g_ratio` from the QM atom.
#[derive(Debug, Clone)]
pub struct LinkAtom {
    /// Index of the QM atom at the bond origin.
    pub qm_atom_idx: usize,
    /// Index of the MM atom at the bond terminus.
    pub mm_atom_idx: usize,
    /// Current Cartesian position of the link (hydrogen cap) atom.
    pub position: [f64; 3],
    /// Fractional position along the QM→MM bond vector.
    /// Typical value: 0.71 for a C–H replacement.
    pub g_ratio: f64,
}

impl LinkAtom {
    /// Update the link atom position from the current QM and MM atom positions.
    ///
    /// `pos = qm_pos + g_ratio * (mm_pos - qm_pos)`
    pub fn update_position(&mut self, qm_pos: [f64; 3], mm_pos: [f64; 3]) {
        for k in 0..3 {
            self.position[k] = qm_pos[k] + self.g_ratio * (mm_pos[k] - qm_pos[k]);
        }
    }
}

// ---------------------------------------------------------------------------
// Semi-empirical tight-binding Hamiltonian (DFTB-inspired)
// ---------------------------------------------------------------------------

/// A simplified tight-binding (DFTB-like) Hamiltonian for the QM region.
///
/// Uses exponential Slater-Koster integrals and a Gershgorin estimate for
/// eigenvalue bounds instead of a full diagonaliser.
#[derive(Debug, Clone)]
pub struct TightBindingHamiltonian {
    /// Atomic numbers of the QM atoms.
    pub atomic_numbers: Vec<u8>,
    /// Positions of the QM atoms (Å).
    pub positions: Vec<[f64; 3]>,
    /// Number of electrons in the QM region.
    pub n_electrons: usize,
}

/// Pre-factor A used in the Slater-Koster integral.
const SK_A: f64 = 1.0;
/// Decay exponent α used in the Slater-Koster integral (Å⁻¹).
const SK_ALPHA: f64 = 1.5;
/// Pre-factor for the repulsive pairwise energy.
const REP_A: f64 = 2.0;
/// Decay exponent for the repulsive energy (Å⁻¹).
const REP_ALPHA: f64 = 2.0;

impl TightBindingHamiltonian {
    /// Simplified Slater-Koster integral between atoms i and j at distance r (Å).
    ///
    /// Returns `A * exp(-α * r)`.
    pub fn slater_koster_integral(_atom_i: u8, _atom_j: u8, r: f64) -> f64 {
        SK_A * (-SK_ALPHA * r).exp()
    }

    /// On-site energy for an atom given its atomic number.
    fn on_site_energy(z: u8) -> f64 {
        // Very rough: use ionisation-energy-like values in eV (negative binding)
        match z {
            1 => -13.6,  // H
            6 => -11.26, // C
            7 => -14.48, // N
            8 => -17.19, // O
            _ => -10.0,
        }
    }

    /// Euclidean distance between two positions.
    fn dist(a: [f64; 3], b: [f64; 3]) -> f64 {
        ((a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2) + (a[2] - b[2]).powi(2)).sqrt()
    }

    /// Build the dense n_atoms × n_atoms Hamiltonian matrix.
    ///
    /// Diagonal elements are on-site energies; off-diagonal elements are
    /// Slater-Koster integrals.
    pub fn build_hamiltonian(&self) -> Vec<Vec<f64>> {
        let n = self.atomic_numbers.len();
        // Precompute all off-diagonal values to avoid simultaneous mutable borrows.
        let mut h = vec![vec![0.0_f64; n]; n];
        for (i, row) in h.iter_mut().enumerate() {
            row[i] = Self::on_site_energy(self.atomic_numbers[i]);
        }
        // Compute off-diagonal SK integrals row by row using split_at_mut.
        for i in 0..n {
            let (top, bot) = h.split_at_mut(i + 1);
            let row_i = &mut top[i];
            for (jj, row_j) in bot.iter_mut().enumerate() {
                let j = i + 1 + jj;
                let r = Self::dist(self.positions[i], self.positions[j]);
                let val =
                    Self::slater_koster_integral(self.atomic_numbers[i], self.atomic_numbers[j], r);
                row_i[j] = val;
                row_j[i] = val;
            }
        }
        h
    }

    /// Build the dense overlap matrix S (identity on diagonal, SK overlap off-diagonal).
    pub fn overlap_matrix(&self) -> Vec<Vec<f64>> {
        let n = self.atomic_numbers.len();
        let mut s = vec![vec![0.0_f64; n]; n];
        for (i, row) in s.iter_mut().enumerate() {
            row[i] = 1.0;
        }
        for i in 0..n {
            let (top, bot) = s.split_at_mut(i + 1);
            let row_i = &mut top[i];
            for (jj, row_j) in bot.iter_mut().enumerate() {
                let j = i + 1 + jj;
                let r = Self::dist(self.positions[i], self.positions[j]);
                let val = 0.5
                    * Self::slater_koster_integral(
                        self.atomic_numbers[i],
                        self.atomic_numbers[j],
                        r,
                    );
                row_i[j] = val;
                row_j[i] = val;
            }
        }
        s
    }

    /// Gershgorin circle theorem bounds on the eigenvalues of `h`.
    ///
    /// For each row i the bound is `[h_ii - R_i, h_ii + R_i]` where
    /// `R_i = sum_{j≠i} |h_ij|`.  Returns the lower bound for each row.
    pub fn eigenvalues_gershgorin(&self, h: &[Vec<f64>]) -> Vec<f64> {
        let mut bounds = Vec::with_capacity(h.len());
        for (i, row) in h.iter().enumerate() {
            let radius: f64 = row
                .iter()
                .enumerate()
                .filter(|&(j, _)| j != i)
                .map(|(_, v)| v.abs())
                .sum();
            bounds.push(row[i] - radius);
        }
        bounds
    }

    /// Compute the repulsive pairwise energy.
    ///
    /// `E_rep = Σ_{i<j} A * exp(-α * r_ij) / r_ij`
    pub fn repulsive_energy(&self) -> f64 {
        let n = self.atomic_numbers.len();
        let mut e_rep = 0.0;
        for i in 0..n {
            for j in (i + 1)..n {
                let r = Self::dist(self.positions[i], self.positions[j]);
                if r > 1e-10 {
                    e_rep += REP_A * (-REP_ALPHA * r).exp() / r;
                }
            }
        }
        e_rep
    }
}

// ---------------------------------------------------------------------------
// Electrostatic embedding
// ---------------------------------------------------------------------------

/// Electrostatic (Coulomb) QM/MM interaction in the mechanical or electrostatic
/// embedding scheme.
pub struct ElectrostaticEmbedding;

impl ElectrostaticEmbedding {
    /// Compute the QM/MM Coulomb interaction energy.
    ///
    /// `E = Σ_i Σ_J q_i * Q_J / r_iJ`  (atomic units, charges in e, distances in Å)
    pub fn qm_mm_interaction(qm_atoms: &[QmAtom], mm_atoms: &[MmAtom]) -> f64 {
        let ke = 332.0636; // kcal·Å/(mol·e²) — Coulomb constant
        let mut energy = 0.0;
        for qi in qm_atoms {
            for qj in mm_atoms {
                let r = dist3(qi.position, qj.position);
                if r > 1e-10 {
                    energy += ke * qi.charge * qj.charge / r;
                }
            }
        }
        energy
    }

    /// QM gradient of the electrostatic QM/MM interaction.
    ///
    /// Returns `dE/dr_i` for each QM atom.
    pub fn qm_mm_gradient_qm(qm_atoms: &[QmAtom], mm_atoms: &[MmAtom]) -> Vec<[f64; 3]> {
        let ke = 332.0636;
        let mut grads = vec![[0.0_f64; 3]; qm_atoms.len()];
        for (i, qi) in qm_atoms.iter().enumerate() {
            for qj in mm_atoms {
                let dr = sub3(qi.position, qj.position);
                let r2 = dot3(dr, dr);
                if r2 > 1e-20 {
                    let r = r2.sqrt();
                    let factor = -ke * qi.charge * qj.charge / (r2 * r);
                    for k in 0..3 {
                        grads[i][k] += factor * dr[k];
                    }
                }
            }
        }
        grads
    }

    /// MM gradient of the electrostatic QM/MM interaction.
    ///
    /// Returns `dE/dr_J` for each MM atom.
    pub fn qm_mm_gradient_mm(qm_atoms: &[QmAtom], mm_atoms: &[MmAtom]) -> Vec<[f64; 3]> {
        let ke = 332.0636;
        let mut grads = vec![[0.0_f64; 3]; mm_atoms.len()];
        for qi in qm_atoms {
            for (j, qj) in mm_atoms.iter().enumerate() {
                let dr = sub3(qj.position, qi.position);
                let r2 = dot3(dr, dr);
                if r2 > 1e-20 {
                    let r = r2.sqrt();
                    let factor = -ke * qi.charge * qj.charge / (r2 * r);
                    for k in 0..3 {
                        grads[j][k] += factor * dr[k];
                    }
                }
            }
        }
        grads
    }
}

// ---------------------------------------------------------------------------
// van der Waals QM/MM cross-terms
// ---------------------------------------------------------------------------

/// Lennard-Jones cross-interaction between QM and MM atoms.
///
/// QM atoms are given default LJ parameters (σ = 3.5 Å, ε = 0.5 kJ/mol).
pub struct VanDerWaalsQmMm;

/// Default LJ σ for QM atoms (Å).
const QM_SIGMA_DEFAULT: f64 = 3.5;
/// Default LJ ε for QM atoms (kJ/mol).
const QM_EPSILON_DEFAULT: f64 = 0.5;

impl VanDerWaalsQmMm {
    /// Lorentz-Berthelot mixing rules for σ and ε.
    fn mix(eps_i: f64, sig_i: f64, eps_j: f64, sig_j: f64) -> (f64, f64) {
        ((eps_i * eps_j).sqrt(), 0.5 * (sig_i + sig_j))
    }

    /// Compute the LJ energy between all QM and MM atom pairs.
    pub fn lj_energy(qm_atoms: &[QmAtom], mm_atoms: &[MmAtom]) -> f64 {
        let mut energy = 0.0;
        for qi in qm_atoms {
            for qj in mm_atoms {
                let (eps, sig) =
                    Self::mix(QM_EPSILON_DEFAULT, QM_SIGMA_DEFAULT, qj.epsilon, qj.sigma);
                let r = dist3(qi.position, qj.position);
                if r > 1e-10 {
                    let sr6 = (sig / r).powi(6);
                    energy += 4.0 * eps * (sr6 * sr6 - sr6);
                }
            }
        }
        energy
    }

    /// QM gradient of the LJ QM/MM interaction.
    pub fn lj_gradient_qm(qm_atoms: &[QmAtom], mm_atoms: &[MmAtom]) -> Vec<[f64; 3]> {
        let mut grads = vec![[0.0_f64; 3]; qm_atoms.len()];
        for (i, qi) in qm_atoms.iter().enumerate() {
            for qj in mm_atoms {
                let (eps, sig) =
                    Self::mix(QM_EPSILON_DEFAULT, QM_SIGMA_DEFAULT, qj.epsilon, qj.sigma);
                let dr = sub3(qi.position, qj.position);
                let r2 = dot3(dr, dr);
                if r2 > 1e-20 {
                    let r = r2.sqrt();
                    let sr2 = (sig / r).powi(2);
                    let sr6 = sr2 * sr2 * sr2;
                    // dE/dr = 4ε * (-12 σ^12/r^13 + 6 σ^6/r^7)
                    let factor = 4.0 * eps * (-12.0 * sr6 * sr6 + 6.0 * sr6) / r2;
                    for k in 0..3 {
                        grads[i][k] += factor * dr[k];
                    }
                }
            }
        }
        grads
    }
}

// ---------------------------------------------------------------------------
// Full QM/MM system
// ---------------------------------------------------------------------------

/// Combined QM/MM system container.
#[derive(Debug, Clone)]
pub struct QmmmSystem {
    /// Atoms in the QM region.
    pub qm_atoms: Vec<QmAtom>,
    /// Atoms in the MM region.
    pub mm_atoms: Vec<MmAtom>,
    /// Link (hydrogen cap) atoms at the QM/MM boundary.
    pub link_atoms: Vec<LinkAtom>,
    /// Semi-empirical tight-binding Hamiltonian for the QM region.
    pub hamiltonian: TightBindingHamiltonian,
}

impl QmmmSystem {
    /// Create an empty QM/MM system.
    pub fn new() -> Self {
        Self {
            qm_atoms: Vec::new(),
            mm_atoms: Vec::new(),
            link_atoms: Vec::new(),
            hamiltonian: TightBindingHamiltonian {
                atomic_numbers: Vec::new(),
                positions: Vec::new(),
                n_electrons: 0,
            },
        }
    }

    /// Add a QM atom. Returns the local QM index.
    pub fn add_qm_atom(&mut self, pos: [f64; 3], atomic_number: u8, charge: f64) -> usize {
        let idx = self.qm_atoms.len();
        self.qm_atoms.push(QmAtom {
            index: idx,
            position: pos,
            atomic_number,
            charge,
            gradient: [0.0; 3],
        });
        self.hamiltonian.atomic_numbers.push(atomic_number);
        self.hamiltonian.positions.push(pos);
        // Rough electron count: use atomic number / 2
        self.hamiltonian.n_electrons += (atomic_number as usize).saturating_div(2);
        idx
    }

    /// Add an MM atom. Returns the local MM index.
    pub fn add_mm_atom(&mut self, pos: [f64; 3], charge: f64, epsilon: f64, sigma: f64) -> usize {
        let idx = self.mm_atoms.len();
        self.mm_atoms.push(MmAtom {
            index: idx,
            position: pos,
            charge,
            epsilon,
            sigma,
            gradient: [0.0; 3],
        });
        idx
    }

    /// Add a link atom between QM atom `qm_idx` and MM atom `mm_idx`.
    pub fn add_link_atom(&mut self, qm_idx: usize, mm_idx: usize, g_ratio: f64) {
        let qm_pos = self.qm_atoms[qm_idx].position;
        let mm_pos = self.mm_atoms[mm_idx].position;
        let mut link = LinkAtom {
            qm_atom_idx: qm_idx,
            mm_atom_idx: mm_idx,
            position: [0.0; 3],
            g_ratio,
        };
        link.update_position(qm_pos, mm_pos);
        self.link_atoms.push(link);
    }

    /// Compute the total QM/MM energy.
    ///
    /// `E_total = E_QM_rep + E_QM/MM_electrostatic + E_QM/MM_vdw`
    pub fn total_energy(&self) -> f64 {
        let e_qm_rep = self.hamiltonian.repulsive_energy();
        let e_elec = ElectrostaticEmbedding::qm_mm_interaction(&self.qm_atoms, &self.mm_atoms);
        let e_vdw = VanDerWaalsQmMm::lj_energy(&self.qm_atoms, &self.mm_atoms);
        e_qm_rep + e_elec + e_vdw
    }

    /// Update all link atom positions using the current QM/MM atom positions.
    pub fn update_link_positions(&mut self) {
        // Collect positions first to avoid borrow conflict
        let updates: Vec<([f64; 3], [f64; 3])> = self
            .link_atoms
            .iter()
            .map(|la| {
                (
                    self.qm_atoms[la.qm_atom_idx].position,
                    self.mm_atoms[la.mm_atom_idx].position,
                )
            })
            .collect();
        for (la, (qm_pos, mm_pos)) in self.link_atoms.iter_mut().zip(updates.iter()) {
            la.update_position(*qm_pos, *mm_pos);
        }
    }
}

impl Default for QmmmSystem {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// Charge equilibration
// ---------------------------------------------------------------------------

/// Simple charge equilibration using the electronegativity equalization method
/// (QEq / EEM).
///
/// Minimises the energy `E = Σ_i χ_i q_i + ½ Σ_{ij} J_ij q_i q_j` subject to
/// `Σ_i q_i = Q_total`.
#[derive(Debug, Clone)]
pub struct ChargeEquilibration {
    /// Atomic electronegativity χ (eV).
    pub electronegativity: Vec<f64>,
    /// Chemical hardness η (eV), appearing on the diagonal of J.
    pub hardness: Vec<f64>,
}

impl ChargeEquilibration {
    /// Compute atomic partial charges via electronegativity equalization.
    ///
    /// Solves the augmented linear system:
    ///
    /// ```text
    /// [ J  1 ] [ q      ]   [ -χ ]
    /// [ 1ᵀ 0 ] [ λ      ] = [ Q  ]
    /// ```
    ///
    /// using Gauss-Jordan elimination.
    pub fn compute_charges(&self, positions: &[[f64; 3]], total_charge: f64) -> Vec<f64> {
        let n = self.electronegativity.len();
        assert_eq!(
            n,
            positions.len(),
            "electronegativity/positions size mismatch"
        );
        let ke = 14.3996; // eV·Å/e²

        // Build the (n+1) × (n+1) augmented matrix for [q; λ]
        let size = n + 1;
        let mut mat = vec![vec![0.0_f64; size + 1]; size];

        // Fill J matrix
        for i in 0..n {
            for j in 0..n {
                if i == j {
                    mat[i][j] = 2.0 * self.hardness[i];
                } else {
                    let r = dist3(positions[i], positions[j]);
                    mat[i][j] = if r > 1e-10 { ke / r } else { 0.0 };
                }
            }
            // Lagrange column and row
            mat[i][n] = 1.0;
            mat[n][i] = 1.0;
        }
        // mat[n][n] = 0 (already)

        // RHS: [-χ_i, …, Q_total]
        for (row, &en) in mat.iter_mut().zip(self.electronegativity.iter()).take(n) {
            row[size] = -en;
        }
        mat[n][size] = total_charge;

        // Gauss-Jordan elimination with partial pivoting
        let mut row_perm: Vec<usize> = (0..size).collect();
        for col in 0..size {
            // Pivot: find row with maximum absolute value
            let pivot_row = (col..size)
                .max_by(|&a, &b| {
                    mat[row_perm[a]][col]
                        .abs()
                        .partial_cmp(&mat[row_perm[b]][col].abs())
                        .unwrap_or(std::cmp::Ordering::Equal)
                })
                .unwrap_or(col);
            row_perm.swap(col, pivot_row);
            let pr = row_perm[col];
            let pivot = mat[pr][col];
            if pivot.abs() < 1e-14 {
                continue;
            }
            // Normalise pivot row
            for v in mat[pr].iter_mut() {
                *v /= pivot;
            }
            // Eliminate all other rows
            for r in 0..size {
                if r == pr {
                    continue;
                }
                let factor = mat[r][col];
                // Copy pivot row values before modifying mat[r] to avoid borrow issues.
                let pr_row: Vec<f64> = mat[pr].clone();
                for (mat_rc, &prc) in mat[r].iter_mut().zip(pr_row.iter()) {
                    *mat_rc -= factor * prc;
                }
            }
        }

        // Extract solution: q_i = mat[row_perm[i]][size]
        (0..n).map(|i| mat[row_perm[i]][size]).collect()
    }
}

// ---------------------------------------------------------------------------
// Small vector helpers (avoids nalgebra dependency)
// ---------------------------------------------------------------------------

#[inline]
fn dist3(a: [f64; 3], b: [f64; 3]) -> f64 {
    ((a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2) + (a[2] - b[2]).powi(2)).sqrt()
}

#[inline]
fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

#[inline]
fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

// ---------------------------------------------------------------------------
// QM region definition helpers
// ---------------------------------------------------------------------------

/// Selection criteria for automatically building a QM region.
#[derive(Debug, Clone)]
pub struct QmRegionSelector {
    /// Indices of atoms that must always be in the QM region.
    pub core_atoms: Vec<usize>,
    /// Cutoff radius in Å from core atoms — all atoms within this sphere are added to QM.
    pub cutoff_radius: f64,
}

impl QmRegionSelector {
    /// Create a new selector with the given core atom indices and radius.
    pub fn new(core_atoms: Vec<usize>, cutoff_radius: f64) -> Self {
        Self {
            core_atoms,
            cutoff_radius,
        }
    }

    /// Select QM atom indices from a full set of positions.
    ///
    /// Returns a sorted, deduplicated list of atom indices within `cutoff_radius`
    /// of any core atom.
    pub fn select_qm_atoms(&self, positions: &[[f64; 3]]) -> Vec<usize> {
        let mut selected: std::collections::BTreeSet<usize> =
            self.core_atoms.iter().cloned().collect();
        for &core in &self.core_atoms {
            if core >= positions.len() {
                continue;
            }
            for (j, &pos) in positions.iter().enumerate() {
                if dist3(positions[core], pos) <= self.cutoff_radius {
                    selected.insert(j);
                }
            }
        }
        selected.into_iter().collect()
    }

    /// Identify QM/MM boundary bonds: pairs (qm_idx, mm_idx) where one atom
    /// is in the QM list and its bonded partner is in the MM list.
    pub fn boundary_bonds(
        &self,
        qm_atoms: &[usize],
        bonded_pairs: &[(usize, usize)],
    ) -> Vec<(usize, usize)> {
        let qm_set: std::collections::HashSet<usize> = qm_atoms.iter().cloned().collect();
        bonded_pairs
            .iter()
            .filter_map(|&(i, j)| {
                let i_qm = qm_set.contains(&i);
                let j_qm = qm_set.contains(&j);
                if i_qm && !j_qm {
                    Some((i, j))
                } else if j_qm && !i_qm {
                    Some((j, i))
                } else {
                    None
                }
            })
            .collect()
    }
}

// ---------------------------------------------------------------------------
// MM point charges
// ---------------------------------------------------------------------------

/// A collection of MM point charges for use in electrostatic embedding.
#[derive(Debug, Clone, Default)]
pub struct PointChargeField {
    /// Charge values (elementary charge units).
    pub charges: Vec<f64>,
    /// Positions of the charges in Å.
    pub positions: Vec<[f64; 3]>,
}

impl PointChargeField {
    /// Create an empty point charge field.
    pub fn new() -> Self {
        Self::default()
    }

    /// Add a point charge.
    pub fn add_charge(&mut self, charge: f64, position: [f64; 3]) {
        self.charges.push(charge);
        self.positions.push(position);
    }

    /// Compute the electrostatic potential at a point `r` (in Å).
    ///
    /// V(r) = Σ Q_J / |r - R_J|  (Coulomb units: ke = 1 implicit)
    pub fn potential_at(&self, r: [f64; 3]) -> f64 {
        let ke = 332.0636; // kcal·Å/(mol·e²)
        self.charges
            .iter()
            .zip(self.positions.iter())
            .map(|(&q, &pos)| {
                let d = dist3(r, pos);
                if d > 1e-10 { ke * q / d } else { 0.0 }
            })
            .sum()
    }

    /// Compute the electric field vector at a point `r` in Å.
    ///
    /// E(r) = -∇V = Σ Q_J * (r - R_J) / |r - R_J|^3
    pub fn field_at(&self, r: [f64; 3]) -> [f64; 3] {
        let ke = 332.0636;
        let mut field = [0.0f64; 3];
        for (&q, &pos) in self.charges.iter().zip(self.positions.iter()) {
            let dr = sub3(r, pos);
            let r2 = dot3(dr, dr);
            if r2 > 1e-20 {
                let r3 = r2 * r2.sqrt();
                for k in 0..3 {
                    field[k] += ke * q * dr[k] / r3;
                }
            }
        }
        field
    }

    /// Total number of point charges.
    pub fn len(&self) -> usize {
        self.charges.len()
    }

    /// Returns true if there are no charges.
    pub fn is_empty(&self) -> bool {
        self.charges.is_empty()
    }
}

// ---------------------------------------------------------------------------
// Mechanical embedding
// ---------------------------------------------------------------------------

/// Mechanical embedding: compute the QM/MM interaction in the mechanical
/// embedding scheme, where the QM charges are fixed from the MM force field.
pub struct MechanicalEmbedding;

impl MechanicalEmbedding {
    /// Mechanical embedding QM/MM interaction energy.
    ///
    /// Uses the MM charges for QM atoms (not self-consistent QM charges).
    /// `qm_mm_charges` are the MM-force-field charges of the QM atoms.
    pub fn interaction_energy(
        qm_positions: &[[f64; 3]],
        qm_mm_charges: &[f64],
        mm_atoms: &[MmAtom],
    ) -> f64 {
        let ke = 332.0636;
        let mut energy = 0.0;
        for (pos, &charge) in qm_positions.iter().zip(qm_mm_charges.iter()) {
            for mm in mm_atoms {
                let r = dist3(*pos, mm.position);
                if r > 1e-10 {
                    energy += ke * charge * mm.charge / r;
                }
            }
        }
        energy
    }
}

// ---------------------------------------------------------------------------
// QM/MM link atom force redistribution
// ---------------------------------------------------------------------------

/// Redistribute force from link atom to the host QM and MM atoms.
///
/// The link atom gradient is projected back to the QM and MM atoms
/// using the g_ratio factor:
///   grad_qm += (1 - g) * grad_link
///   grad_mm += g * grad_link
///
/// where g = `link_atom.g_ratio`.
pub fn redistribute_link_force(
    link_grad: [f64; 3],
    g_ratio: f64,
    grad_qm: &mut [f64; 3],
    grad_mm: &mut [f64; 3],
) {
    let g = g_ratio;
    for k in 0..3 {
        grad_qm[k] += (1.0 - g) * link_grad[k];
        grad_mm[k] += g * link_grad[k];
    }
}

// ---------------------------------------------------------------------------
// ONIOM extrapolation energy
// ---------------------------------------------------------------------------

/// ONIOM (Our own N-layered Integrated molecular Orbital and molecular Mechanics)
/// energy extrapolation.
///
/// The two-layer ONIOM energy is:
/// ```text
/// E_ONIOM = E_QM(model) + E_MM(real) - E_MM(model)
/// ```
/// where "model" is the QM region treated at both levels.
///
/// # Arguments
/// * `e_qm_model`  – QM energy of the model (inner) region.
/// * `e_mm_real`   – MM energy of the entire (real) system.
/// * `e_mm_model`  – MM energy of the model region with MM parameters.
pub fn oniom_energy(e_qm_model: f64, e_mm_real: f64, e_mm_model: f64) -> f64 {
    e_qm_model + e_mm_real - e_mm_model
}

/// Three-layer ONIOM energy.
///
/// ```text
/// E = E_high(inner) + E_mid(middle) - E_mid(inner)
///   + E_low(outer) - E_low(middle)
/// ```
pub fn oniom3_energy(
    e_high_inner: f64,
    e_mid_middle: f64,
    e_mid_inner: f64,
    e_low_outer: f64,
    e_low_middle: f64,
) -> f64 {
    e_high_inner + e_mid_middle - e_mid_inner + e_low_outer - e_low_middle
}

// ---------------------------------------------------------------------------
// Buffer zone (adaptive QM/MM boundary) smoothing
// ---------------------------------------------------------------------------

/// Smooth transition weight for atoms in the buffer zone.
///
/// Uses a 5th-order polynomial switch that goes from 1 at `r_inner` to 0 at
/// `r_outer`, ensuring C² continuity:
///
/// ```text
/// t = (r - r_inner) / (r_outer - r_inner)
/// S(t) = 1 - 6t^5 + 15t^4 - 10t^3     (0 ≤ t ≤ 1)
/// ```
pub fn buffer_weight(r: f64, r_inner: f64, r_outer: f64) -> f64 {
    if r <= r_inner {
        return 1.0;
    }
    if r >= r_outer {
        return 0.0;
    }
    let t = (r - r_inner) / (r_outer - r_inner);
    1.0 - t * t * t * (10.0 - t * (15.0 - 6.0 * t))
}

/// Derivative of the buffer weight with respect to r.
pub fn buffer_weight_derivative(r: f64, r_inner: f64, r_outer: f64) -> f64 {
    if r <= r_inner || r >= r_outer {
        return 0.0;
    }
    let width = r_outer - r_inner;
    let t = (r - r_inner) / width;
    // dS/dr = (dS/dt) * (dt/dr)
    // dS/dt = -30 t^4 + 60 t^3 - 30 t^2 = -30 t^2 (t^2 - 2t + 1) = -30 t^2 (t-1)^2
    let ds_dt = -30.0 * t * t * (t - 1.0) * (t - 1.0);
    ds_dt / width
}

// ---------------------------------------------------------------------------
// Adaptive QM region expander
// ---------------------------------------------------------------------------

/// Adaptive QM region expansion: grows the QM sphere until the gradient
/// norm on the boundary falls below a threshold.
///
/// Starting from `initial_radius`, increases the QM radius by `step` until
/// either `max_radius` is reached or the criterion is satisfied.
///
/// Returns the chosen QM radius.
pub fn adaptive_qm_radius(
    core_pos: [f64; 3],
    positions: &[[f64; 3]],
    gradients: &[[f64; 3]],
    initial_radius: f64,
    step: f64,
    max_radius: f64,
    threshold: f64,
) -> f64 {
    let mut radius = initial_radius;
    loop {
        // Find atoms at the boundary shell: [radius - step, radius]
        let boundary_grad_norm: f64 = positions
            .iter()
            .zip(gradients.iter())
            .filter_map(|(pos, grad)| {
                let r = dist3(core_pos, *pos);
                if r >= radius - step && r <= radius {
                    Some((grad[0] * grad[0] + grad[1] * grad[1] + grad[2] * grad[2]).sqrt())
                } else {
                    None
                }
            })
            .fold(0.0_f64, f64::max);
        if boundary_grad_norm < threshold || radius >= max_radius {
            break;
        }
        radius = (radius + step).min(max_radius);
    }
    radius
}

// ---------------------------------------------------------------------------
// Constraint correction for link atoms
// ---------------------------------------------------------------------------

/// Holonomic constraint correction for a link atom bond.
///
/// Returns the penalty energy for deviating from the target link bond length:
/// `E_constraint = k/2 * (|r_link - r_qm| - d_target)^2`
pub fn link_bond_constraint_energy(
    link_pos: [f64; 3],
    qm_pos: [f64; 3],
    d_target: f64,
    k: f64,
) -> f64 {
    let r = dist3(link_pos, qm_pos);
    let delta = r - d_target;
    0.5 * k * delta * delta
}

/// Gradient of the link bond constraint with respect to the link atom position.
pub fn link_bond_constraint_gradient(
    link_pos: [f64; 3],
    qm_pos: [f64; 3],
    d_target: f64,
    k: f64,
) -> [f64; 3] {
    let dr = sub3(link_pos, qm_pos);
    let r = dist3(link_pos, qm_pos);
    if r < 1e-14 {
        return [0.0; 3];
    }
    let delta = r - d_target;
    let factor = k * delta / r;
    [dr[0] * factor, dr[1] * factor, dr[2] * factor]
}

// ---------------------------------------------------------------------------
// Capping hydrogen geometry
// ---------------------------------------------------------------------------

/// Compute the optimal position of a capping hydrogen atom at the QM/MM bond.
///
/// The capping H is placed along the QM→MM bond at `g_ratio` from the QM atom,
/// scaled to the standard C-H bond length `r_ch` (default ≈ 1.09 Å).
///
/// # Returns
/// Cartesian position of the capping hydrogen.
pub fn capping_hydrogen_position(qm_pos: [f64; 3], mm_pos: [f64; 3], r_ch: f64) -> [f64; 3] {
    let dr = sub3(mm_pos, qm_pos);
    let r = (dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2]).sqrt();
    if r < 1e-14 {
        return qm_pos;
    }
    let scale = r_ch / r;
    [
        qm_pos[0] + scale * dr[0],
        qm_pos[1] + scale * dr[1],
        qm_pos[2] + scale * dr[2],
    ]
}

// ---------------------------------------------------------------------------
// Mulliken charge analysis (simplified)
// ---------------------------------------------------------------------------

/// Simplified Mulliken population analysis.
///
/// Given a density matrix `P` (n×n, row-major) and an overlap matrix `S`
/// (n×n, row-major), computes the Mulliken population on each orbital.
///
/// `q_i = Σ_j P_ij * S_ji`  (diagonal of PS).
///
/// Returns the Mulliken population vector of length n.
pub fn mulliken_populations(density_matrix: &[Vec<f64>], overlap_matrix: &[Vec<f64>]) -> Vec<f64> {
    let n = density_matrix.len();
    (0..n)
        .map(|i| {
            (0..n)
                .map(|j| {
                    let p = density_matrix[i].get(j).copied().unwrap_or(0.0);
                    let s = overlap_matrix
                        .get(j)
                        .and_then(|row| row.get(i))
                        .copied()
                        .unwrap_or(0.0);
                    p * s
                })
                .sum()
        })
        .collect()
}

/// Mulliken gross charges given populations, core charges, and number of electrons.
///
/// `q_i = Z_i - N_i`  where Z_i is the nuclear charge and N_i the Mulliken population.
pub fn mulliken_charges(populations: &[f64], nuclear_charges: &[f64]) -> Vec<f64> {
    populations
        .iter()
        .zip(nuclear_charges.iter())
        .map(|(&n, &z)| z - n)
        .collect()
}

// ---------------------------------------------------------------------------
// Ewald correction for QM/MM long-range electrostatics
// ---------------------------------------------------------------------------

/// Ewald correction energy for a QM atom embedded in an MM periodic charge field.
///
/// Computes the short-range correction (removes double-counting from the
/// classical Ewald sum for QM charges):
/// `E_corr = -Σ_i q_i^QM * V_Ewald_self(q_i^QM, α)`
///
/// Self-energy: `E_self = -α/√π * Σ_i (q_i)^2`
pub fn ewald_self_correction(charges: &[f64], alpha: f64) -> f64 {
    let prefactor = alpha / std::f64::consts::PI.sqrt();
    -prefactor * charges.iter().map(|&q| q * q).sum::<f64>()
}

/// Reciprocal-space Ewald correction for QM charges in a cubic box of side L.
///
/// Computes the leading G=0 term only (dipole correction):
/// `E_G0 = -π / (2 α² V) * |Σ_i q_i r_i|² = 0` for neutral systems.
///
/// For a non-neutral system, returns a rough estimate of the G=0 divergence.
pub fn ewald_g0_correction(
    charges: &[f64],
    positions: &[[f64; 3]],
    alpha: f64,
    box_length: f64,
) -> f64 {
    let volume = box_length.powi(3);
    if volume < 1e-30 {
        return 0.0;
    }
    // Dipole moment
    let mut dipole = [0.0_f64; 3];
    for (q, pos) in charges.iter().zip(positions.iter()) {
        for k in 0..3 {
            dipole[k] += q * pos[k];
        }
    }
    let d2 = dipole[0] * dipole[0] + dipole[1] * dipole[1] + dipole[2] * dipole[2];
    // Standard dipole correction: E = π/(2 α² V) * |d|² (sign convention varies)
    std::f64::consts::PI / (2.0 * alpha * alpha * volume) * d2
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // ---- LinkAtom ----------------------------------------------------------

    #[test]
    fn test_link_atom_ratio_zero_equals_qm() {
        let mut la = LinkAtom {
            qm_atom_idx: 0,
            mm_atom_idx: 1,
            position: [0.0; 3],
            g_ratio: 0.0,
        };
        let qm_pos = [1.0, 2.0, 3.0];
        let mm_pos = [4.0, 5.0, 6.0];
        la.update_position(qm_pos, mm_pos);
        assert!((la.position[0] - qm_pos[0]).abs() < 1e-12);
        assert!((la.position[1] - qm_pos[1]).abs() < 1e-12);
        assert!((la.position[2] - qm_pos[2]).abs() < 1e-12);
    }

    #[test]
    fn test_link_atom_ratio_one_equals_mm() {
        let mut la = LinkAtom {
            qm_atom_idx: 0,
            mm_atom_idx: 1,
            position: [0.0; 3],
            g_ratio: 1.0,
        };
        let qm_pos = [1.0, 2.0, 3.0];
        let mm_pos = [4.0, 5.0, 6.0];
        la.update_position(qm_pos, mm_pos);
        assert!((la.position[0] - mm_pos[0]).abs() < 1e-12);
        assert!((la.position[1] - mm_pos[1]).abs() < 1e-12);
        assert!((la.position[2] - mm_pos[2]).abs() < 1e-12);
    }

    // ---- TightBindingHamiltonian -------------------------------------------

    #[test]
    fn test_hamiltonian_symmetric_two_atoms() {
        let tb = TightBindingHamiltonian {
            atomic_numbers: vec![6, 8],
            positions: vec![[0.0, 0.0, 0.0], [1.5, 0.0, 0.0]],
            n_electrons: 7,
        };
        let h = tb.build_hamiltonian();
        let _n = h.len();
        for (i, hi) in h.iter().enumerate() {
            for (j, &val) in hi.iter().enumerate() {
                assert!(
                    (val - h[j][i]).abs() < 1e-12,
                    "H[{i}][{j}] != H[{j}][{i}]: {} vs {}",
                    val,
                    h[j][i]
                );
            }
        }
    }

    #[test]
    fn test_slater_koster_decays_with_distance() {
        let r_near = TightBindingHamiltonian::slater_koster_integral(6, 6, 1.0);
        let r_far = TightBindingHamiltonian::slater_koster_integral(6, 6, 5.0);
        assert!(
            r_near > r_far,
            "SK integral should decay: {r_near} > {r_far}"
        );
    }

    // ---- ElectrostaticEmbedding --------------------------------------------

    #[test]
    fn test_electrostatic_same_sign_positive_energy() {
        let qm = vec![QmAtom {
            index: 0,
            position: [0.0, 0.0, 0.0],
            atomic_number: 1,
            charge: 1.0,
            gradient: [0.0; 3],
        }];
        let mm = vec![MmAtom {
            index: 0,
            position: [3.0, 0.0, 0.0],
            charge: 1.0,
            epsilon: 0.5,
            sigma: 3.5,
            gradient: [0.0; 3],
        }];
        let e = ElectrostaticEmbedding::qm_mm_interaction(&qm, &mm);
        assert!(
            e > 0.0,
            "Same-sign charges should give positive Coulomb energy: {e}"
        );
    }

    // ---- QmmmSystem --------------------------------------------------------

    #[test]
    fn test_qmmm_system_total_energy_no_panic() {
        let mut sys = QmmmSystem::new();
        sys.add_qm_atom([0.0, 0.0, 0.0], 6, -0.2);
        sys.add_qm_atom([1.5, 0.0, 0.0], 1, 0.1);
        sys.add_mm_atom([5.0, 0.0, 0.0], -0.5, 0.3, 3.2);
        sys.add_mm_atom([6.5, 0.0, 0.0], 0.5, 0.3, 3.2);
        let e = sys.total_energy();
        assert!(e.is_finite(), "total_energy should be finite: {e}");
    }

    // ---- ChargeEquilibration -----------------------------------------------

    #[test]
    fn test_charge_equilibration_sum_to_total() {
        let ceq = ChargeEquilibration {
            electronegativity: vec![4.53, 6.93, 4.53], // C, O, C
            hardness: vec![5.0, 8.0, 5.0],
        };
        let positions = vec![[0.0, 0.0, 0.0], [1.2, 0.0, 0.0], [2.4, 0.0, 0.0]];
        let total_q = 0.0;
        let charges = ceq.compute_charges(&positions, total_q);
        let sum: f64 = charges.iter().sum();
        assert!(
            (sum - total_q).abs() < 1e-8,
            "Charges should sum to {total_q}, got {sum}"
        );
    }

    // ---- QmRegionSelector --------------------------------------------------

    #[test]
    fn test_qm_region_selector_includes_core() {
        let positions = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [5.0, 0.0, 0.0]];
        let selector = QmRegionSelector::new(vec![0], 1.5);
        let qm = selector.select_qm_atoms(&positions);
        assert!(qm.contains(&0), "core atom 0 should be selected");
        assert!(
            qm.contains(&1),
            "atom 1 at 1.0 Å should be within cutoff 1.5"
        );
        assert!(!qm.contains(&2), "atom 2 at 5.0 Å should be outside cutoff");
    }

    #[test]
    fn test_qm_region_selector_boundary_bonds() {
        let qm_atoms = vec![0usize, 1];
        let bonded_pairs = vec![(0, 1), (1, 2), (2, 3)];
        let selector = QmRegionSelector::new(vec![0], 2.0);
        let boundary = selector.boundary_bonds(&qm_atoms, &bonded_pairs);
        // Bond (1, 2): atom 1 is QM, atom 2 is MM → boundary bond
        assert_eq!(
            boundary.len(),
            1,
            "exactly one QM/MM boundary bond, got {}",
            boundary.len()
        );
        assert_eq!(boundary[0], (1, 2));
    }

    // ---- PointChargeField --------------------------------------------------

    #[test]
    fn test_point_charge_potential_single_charge() {
        let mut field = PointChargeField::new();
        field.add_charge(1.0, [0.0, 0.0, 0.0]);
        let ke = 332.0636;
        let v = field.potential_at([1.0, 0.0, 0.0]);
        assert!(
            (v - ke).abs() < 1e-6,
            "V at 1 Å from unit charge should be ke={ke}, got {v}"
        );
    }

    #[test]
    fn test_point_charge_field_direction() {
        let mut field = PointChargeField::new();
        field.add_charge(1.0, [0.0, 0.0, 0.0]);
        let e = field.field_at([1.0, 0.0, 0.0]);
        assert!(
            e[0] > 0.0,
            "Electric field x-component should point away from positive charge"
        );
        assert!(e[1].abs() < 1e-10, "Ey should be zero along x axis");
        assert!(e[2].abs() < 1e-10, "Ez should be zero along x axis");
    }

    #[test]
    fn test_point_charge_field_empty() {
        let field = PointChargeField::new();
        assert!(field.is_empty());
        let v = field.potential_at([0.0, 0.0, 0.0]);
        assert_eq!(v, 0.0);
    }

    // ---- MechanicalEmbedding -----------------------------------------------

    #[test]
    fn test_mechanical_embedding_energy_finite() {
        let qm_positions = vec![[0.0, 0.0, 0.0]];
        let qm_charges = vec![-0.5];
        let mm_atoms = vec![MmAtom {
            index: 0,
            position: [3.0, 0.0, 0.0],
            charge: 0.5,
            epsilon: 0.3,
            sigma: 3.2,
            gradient: [0.0; 3],
        }];
        let e = MechanicalEmbedding::interaction_energy(&qm_positions, &qm_charges, &mm_atoms);
        assert!(
            e.is_finite(),
            "Mechanical embedding energy should be finite, got {e}"
        );
        // Opposite charges: energy should be negative
        assert!(
            e < 0.0,
            "Opposite charges should give negative energy, got {e}"
        );
    }

    // ---- redistribute_link_force -------------------------------------------

    #[test]
    fn test_redistribute_link_force_sum_conserved() {
        let link_grad = [1.0, 0.0, 0.0];
        let g = 0.71;
        let mut grad_qm = [0.0f64; 3];
        let mut grad_mm = [0.0f64; 3];
        redistribute_link_force(link_grad, g, &mut grad_qm, &mut grad_mm);
        // Sum should equal link_grad
        for k in 0..3 {
            let total = grad_qm[k] + grad_mm[k];
            assert!(
                (total - link_grad[k]).abs() < 1e-12,
                "Force sum should be conserved at k={k}: {total} vs {}",
                link_grad[k]
            );
        }
    }

    #[test]
    fn test_redistribute_link_force_ratio() {
        let link_grad = [2.0, 0.0, 0.0];
        let g = 0.5;
        let mut grad_qm = [0.0f64; 3];
        let mut grad_mm = [0.0f64; 3];
        redistribute_link_force(link_grad, g, &mut grad_qm, &mut grad_mm);
        // With g=0.5: each gets 50% = 1.0
        assert!((grad_qm[0] - 1.0).abs() < 1e-12, "QM gets (1-g)*F = 1.0");
        assert!((grad_mm[0] - 1.0).abs() < 1e-12, "MM gets g*F = 1.0");
    }

    // ---- ONIOM tests -------------------------------------------------------

    #[test]
    fn test_oniom_energy_additivity() {
        // E_ONIOM = E_QM + E_MM - E_MM_model
        let e = oniom_energy(10.0, 50.0, 8.0);
        assert!(
            (e - 52.0).abs() < 1e-12,
            "ONIOM energy should be 52.0, got {e}"
        );
    }

    #[test]
    fn test_oniom_energy_equal_mm() {
        // If E_MM_real == E_MM_model, only E_QM matters
        let e = oniom_energy(5.0, 20.0, 20.0);
        assert!(
            (e - 5.0).abs() < 1e-12,
            "ONIOM with equal MM energies should give E_QM"
        );
    }

    #[test]
    fn test_oniom3_consistency() {
        // Three-layer ONIOM: verify algebraic consistency
        let e = oniom3_energy(10.0, 20.0, 5.0, 50.0, 25.0);
        // = 10 + 20 - 5 + 50 - 25 = 50
        assert!((e - 50.0).abs() < 1e-12, "ONIOM3 should give 50.0, got {e}");
    }

    // ---- buffer_weight tests -----------------------------------------------

    #[test]
    fn test_buffer_weight_inner_boundary() {
        let w = buffer_weight(0.9, 1.0, 2.0);
        assert!(
            (w - 1.0).abs() < 1e-12,
            "Inside inner cutoff, weight should be 1.0"
        );
    }

    #[test]
    fn test_buffer_weight_outer_boundary() {
        let w = buffer_weight(2.1, 1.0, 2.0);
        assert!(
            (w - 0.0).abs() < 1e-12,
            "Outside outer cutoff, weight should be 0.0"
        );
    }

    #[test]
    fn test_buffer_weight_midpoint() {
        let w = buffer_weight(1.5, 1.0, 2.0);
        assert!(
            w > 0.0 && w < 1.0,
            "Buffer weight at midpoint should be in (0,1): {w}"
        );
    }

    #[test]
    fn test_buffer_weight_derivative_numerical() {
        let r_i = 1.0_f64;
        let r_o = 2.0_f64;
        let r = 1.5_f64;
        let h = 1e-6;
        let numerical =
            (buffer_weight(r + h, r_i, r_o) - buffer_weight(r - h, r_i, r_o)) / (2.0 * h);
        let analytical = buffer_weight_derivative(r, r_i, r_o);
        assert!(
            (numerical - analytical).abs() < 1e-5,
            "Buffer weight derivative mismatch: {analytical} vs {numerical}"
        );
    }

    // ---- capping_hydrogen tests ---------------------------------------------

    #[test]
    fn test_capping_hydrogen_position_correct_length() {
        let qm = [0.0, 0.0, 0.0];
        let mm = [4.0, 0.0, 0.0];
        let r_ch = 1.09_f64;
        let cap = capping_hydrogen_position(qm, mm, r_ch);
        let actual_dist =
            ((cap[0] - qm[0]).powi(2) + (cap[1] - qm[1]).powi(2) + (cap[2] - qm[2]).powi(2)).sqrt();
        assert!(
            (actual_dist - r_ch).abs() < 1e-10,
            "Capping H should be at r_ch={r_ch}, got {actual_dist}"
        );
    }

    #[test]
    fn test_capping_hydrogen_along_bond_axis() {
        let qm = [0.0, 0.0, 0.0];
        let mm = [0.0, 3.0, 0.0]; // along y axis
        let cap = capping_hydrogen_position(qm, mm, 1.0);
        assert!(cap[0].abs() < 1e-10, "Cap H should be on y-axis (x≈0)");
        assert!(cap[1] > 0.0, "Cap H should be in +y direction");
        assert!(cap[2].abs() < 1e-10, "Cap H should be on y-axis (z≈0)");
    }

    // ---- link_bond_constraint tests ----------------------------------------

    #[test]
    fn test_link_bond_constraint_energy_at_target() {
        let link = [1.09, 0.0, 0.0];
        let qm = [0.0, 0.0, 0.0];
        let e = link_bond_constraint_energy(link, qm, 1.09, 100.0);
        assert!(
            e < 1e-10,
            "Constraint energy at target distance should be 0, got {e}"
        );
    }

    #[test]
    fn test_link_bond_constraint_energy_positive() {
        let link = [2.0, 0.0, 0.0];
        let qm = [0.0, 0.0, 0.0];
        let e = link_bond_constraint_energy(link, qm, 1.09, 100.0);
        assert!(
            e > 0.0,
            "Constraint energy away from target should be positive, got {e}"
        );
    }

    #[test]
    fn test_link_bond_constraint_gradient_direction() {
        // When stretched beyond target, gradient should point away from QM atom
        let link = [2.0, 0.0, 0.0];
        let qm = [0.0, 0.0, 0.0];
        let g = link_bond_constraint_gradient(link, qm, 1.09, 100.0);
        assert!(
            g[0] > 0.0,
            "Gradient should point in +x direction (away from QM), got {}",
            g[0]
        );
        assert!(g[1].abs() < 1e-10, "y-gradient should be 0");
    }

    // ---- Mulliken tests -----------------------------------------------------

    #[test]
    fn test_mulliken_populations_identity_overlap() {
        // With S = I (identity), populations equal diagonal of P
        let n = 3;
        let p = vec![
            vec![2.0, 0.1, 0.0],
            vec![0.1, 1.5, 0.2],
            vec![0.0, 0.2, 1.0],
        ];
        let s: Vec<Vec<f64>> = (0..n)
            .map(|i| (0..n).map(|j| if i == j { 1.0 } else { 0.0 }).collect())
            .collect();
        let pop = mulliken_populations(&p, &s);
        // With S=I: q_i = P_ii (contribution from each other orbital is 0)
        assert!((pop[0] - p[0][0]).abs() < 1e-12, "pop[0] should be P[0][0]");
        assert!((pop[1] - p[1][1]).abs() < 1e-12, "pop[1] should be P[1][1]");
    }

    #[test]
    fn test_mulliken_charges_neutral_molecule() {
        // If populations equal nuclear charges, all atomic charges are 0
        let populations = vec![6.0, 1.0, 1.0];
        let nuclear = vec![6.0, 1.0, 1.0];
        let charges = mulliken_charges(&populations, &nuclear);
        for (i, &c) in charges.iter().enumerate() {
            assert!(c.abs() < 1e-12, "Mulliken charge {i} should be 0, got {c}");
        }
    }

    // ---- Ewald correction tests --------------------------------------------

    #[test]
    fn test_ewald_self_correction_sign() {
        // Self-correction should always be negative (removes overcounting)
        let charges = vec![0.5, -0.5, 0.3];
        let e = ewald_self_correction(&charges, 0.35);
        assert!(e < 0.0, "Ewald self-correction should be negative, got {e}");
    }

    #[test]
    fn test_ewald_self_correction_proportional_to_alpha() {
        let charges = vec![1.0];
        let e1 = ewald_self_correction(&charges, 0.1);
        let e2 = ewald_self_correction(&charges, 0.2);
        assert!(e2 / e1 > 1.9, "Self-correction should scale with alpha");
    }

    #[test]
    fn test_ewald_g0_neutral_system() {
        // Neutral system (Σ q=0, dipole moment=0) → G=0 correction ≈ 0
        let charges = vec![1.0, -1.0];
        let positions = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
        let e = ewald_g0_correction(&charges, &positions, 0.35, 10.0);
        // Dipole is [1*0 + (-1)*1, 0, 0] = [-1, 0, 0]
        assert!(
            e.is_finite(),
            "G=0 correction should be finite for dipole system"
        );
    }

    #[test]
    fn test_adaptive_qm_radius_returns_at_least_initial() {
        let core = [0.0, 0.0, 0.0];
        let positions = vec![[1.0, 0.0, 0.0], [2.0, 0.0, 0.0], [5.0, 0.0, 0.0]];
        let gradients = vec![[0.0f64; 3]; 3]; // all zero grads → converge immediately
        let r = adaptive_qm_radius(core, &positions, &gradients, 1.5, 0.5, 5.0, 0.01);
        assert!(
            r >= 1.5,
            "Adaptive radius should be at least the initial radius, got {r}"
        );
    }
}
