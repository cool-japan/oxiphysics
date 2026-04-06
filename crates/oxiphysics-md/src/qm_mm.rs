// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! QM/MM hybrid quantum mechanics / molecular mechanics methods.
//!
//! Provides a framework for QM/MM simulations where a small quantum region
//! is embedded in a larger classical MM region.  The QM method is a
//! placeholder (PM3/SCC-DFTB interface), while the MM region uses standard
//! force-field point charges.
//!
//! ## Energy decomposition
//!
//! E_total = E_QM + E_MM + E_QM/MM
//!
//! where E_QM/MM can be computed via:
//! - **Mechanical embedding**: simple force mixing at the boundary
//! - **Electrostatic embedding**: QM wavefunction polarised by MM point charges
//!
//! ## Boundary treatment
//!
//! Covalent bonds that cross the QM/MM boundary are capped with hydrogen
//! link atoms to saturate the QM region valence.
//!
//! References:
//! - Field, M. J., Bash, P. A., & Karplus, M. (1990). *J. Comput. Chem.* 11, 700.
//! - Bakowies, D., & Thiel, W. (1996). *J. Phys. Chem.* 100, 10580.

// ---------------------------------------------------------------------------
// QmMethod — placeholder enum for QM level of theory
// ---------------------------------------------------------------------------

/// QM level of theory (placeholder implementations).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum QmMethod {
    /// PM3 semi-empirical method.
    Pm3,
    /// Self-consistent-charge density-functional tight-binding (SCC-DFTB).
    SccDftb,
    /// Hartree-Fock (minimal basis, placeholder).
    Hf,
    /// DFT/B3LYP (placeholder).
    DftB3lyp,
}

impl QmMethod {
    /// Human-readable name for this method.
    pub fn name(&self) -> &'static str {
        match self {
            Self::Pm3 => "PM3",
            Self::SccDftb => "SCC-DFTB",
            Self::Hf => "HF",
            Self::DftB3lyp => "DFT/B3LYP",
        }
    }
}

// ---------------------------------------------------------------------------
// QmRegion
// ---------------------------------------------------------------------------

/// The quantum-mechanical region of a QM/MM system.
///
/// Contains the atom indices that are treated quantum mechanically, the
/// chosen QM method, and the total charge of the QM region.
#[allow(dead_code)]
pub struct QmRegion {
    /// Indices of atoms in the QM region (into the parent system's atom list).
    pub atom_indices: Vec<usize>,
    /// QM method to use.
    pub method: QmMethod,
    /// Total formal charge of the QM region.
    pub charge: i32,
    /// Spin multiplicity (2S+1).
    pub multiplicity: u32,
    /// Positions of QM atoms (x, y, z) in Å.
    pub positions: Vec<[f64; 3]>,
    /// Atomic numbers of QM atoms.
    pub atomic_numbers: Vec<u32>,
    /// QM energy (in kcal/mol, filled after a QM calculation).
    pub energy: f64,
    /// QM forces on QM atoms \[\[fx, fy, fz\\]; n_qm].
    pub forces: Vec<[f64; 3]>,
}

impl QmRegion {
    /// Create a new QM region.
    ///
    /// # Arguments
    /// * `atom_indices` — indices into the full system
    /// * `method` — QM level of theory
    /// * `charge` — total charge
    /// * `multiplicity` — spin multiplicity (1 = singlet)
    pub fn new(atom_indices: Vec<usize>, method: QmMethod, charge: i32, multiplicity: u32) -> Self {
        let n = atom_indices.len();
        Self {
            atom_indices,
            method,
            charge,
            multiplicity,
            positions: vec![[0.0; 3]; n],
            atomic_numbers: vec![1; n],
            energy: 0.0,
            forces: vec![[0.0; 3]; n],
        }
    }

    /// Number of QM atoms.
    pub fn n_atoms(&self) -> usize {
        self.atom_indices.len()
    }

    /// Set positions of QM atoms (in Å).
    pub fn set_positions(&mut self, positions: Vec<[f64; 3]>) {
        self.positions = positions;
    }

    /// Placeholder: run the QM calculation and store energy/forces.
    ///
    /// Uses a simple harmonic model as a stand-in for a real QM solver.
    pub fn run_placeholder(&mut self) {
        // Placeholder: E = sum of |r_i|^2 (harmonic well)
        let mut e = 0.0;
        for (i, pos) in self.positions.iter().enumerate() {
            let r2 = pos[0] * pos[0] + pos[1] * pos[1] + pos[2] * pos[2];
            e += r2;
            self.forces[i] = [-2.0 * pos[0], -2.0 * pos[1], -2.0 * pos[2]];
        }
        self.energy = e;
    }
}

// ---------------------------------------------------------------------------
// MmRegion
// ---------------------------------------------------------------------------

/// The classical molecular mechanics region.
#[allow(dead_code)]
pub struct MmRegion {
    /// Indices of atoms in the MM region.
    pub atom_indices: Vec<usize>,
    /// Positions of MM atoms (x, y, z) in Å.
    pub positions: Vec<[f64; 3]>,
    /// Partial charges of MM atoms in units of elementary charge.
    pub charges: Vec<f64>,
    /// Epsilon (LJ well depth) for each MM atom in kcal/mol.
    pub epsilon: Vec<f64>,
    /// Sigma (LJ diameter) for each MM atom in Å.
    pub sigma: Vec<f64>,
    /// MM energy (kcal/mol).
    pub energy: f64,
    /// MM forces on MM atoms.
    pub forces: Vec<[f64; 3]>,
}

impl MmRegion {
    /// Create a new MM region.
    pub fn new(atom_indices: Vec<usize>) -> Self {
        let n = atom_indices.len();
        Self {
            atom_indices,
            positions: vec![[0.0; 3]; n],
            charges: vec![0.0; n],
            epsilon: vec![0.1; n],
            sigma: vec![3.5; n],
            energy: 0.0,
            forces: vec![[0.0; 3]; n],
        }
    }

    /// Number of MM atoms.
    pub fn n_atoms(&self) -> usize {
        self.atom_indices.len()
    }

    /// Placeholder MM energy: pairwise LJ between MM atoms.
    pub fn run_placeholder(&mut self) {
        let n = self.positions.len();
        let mut e = 0.0;
        for i in 0..n {
            for j in (i + 1)..n {
                let dx = self.positions[i][0] - self.positions[j][0];
                let dy = self.positions[i][1] - self.positions[j][1];
                let dz = self.positions[i][2] - self.positions[j][2];
                let r2 = dx * dx + dy * dy + dz * dz;
                if r2 < 1e-10 {
                    continue;
                }
                let sigma_ij = 0.5 * (self.sigma[i] + self.sigma[j]);
                let eps_ij = (self.epsilon[i] * self.epsilon[j]).sqrt();
                let sr2 = (sigma_ij * sigma_ij) / r2;
                let sr6 = sr2 * sr2 * sr2;
                let sr12 = sr6 * sr6;
                e += 4.0 * eps_ij * (sr12 - sr6);
            }
        }
        self.energy = e;
    }
}

// ---------------------------------------------------------------------------
// QmMmSystem
// ---------------------------------------------------------------------------

/// Combined QM/MM system.
#[allow(dead_code)]
pub struct QmMmSystem {
    /// The QM region.
    pub qm: QmRegion,
    /// The MM region.
    pub mm: MmRegion,
    /// QM/MM interaction energy (kcal/mol).
    pub energy_qm_mm: f64,
    /// Link atoms inserted at QM/MM boundaries.
    pub link_atoms: Vec<LinkAtom>,
}

impl QmMmSystem {
    /// Create a new QM/MM system from QM and MM regions.
    pub fn new(qm: QmRegion, mm: MmRegion) -> Self {
        Self {
            qm,
            mm,
            energy_qm_mm: 0.0,
            link_atoms: Vec::new(),
        }
    }

    /// Total energy: E_QM + E_MM + E_QM/MM.
    pub fn total_energy(&self) -> f64 {
        qm_mm_energy(self.qm.energy, self.mm.energy, self.energy_qm_mm)
    }

    /// Add a link atom at the QM/MM boundary.
    pub fn add_link_atom(&mut self, la: LinkAtom) {
        self.link_atoms.push(la);
    }
}

// ---------------------------------------------------------------------------
// LinkAtom
// ---------------------------------------------------------------------------

/// A hydrogen capping link atom at a QM/MM boundary bond.
///
/// The link atom (always hydrogen) is placed along the bond between a
/// QM atom and a MM atom, at a fixed distance g × r_{QM-MM}.
#[allow(dead_code)]
pub struct LinkAtom {
    /// QM atom index (in QM region) adjacent to the boundary.
    pub qm_atom: usize,
    /// MM atom index (in MM region) that the QM atom is bonded to.
    pub mm_atom: usize,
    /// Scale factor g ∈ (0, 1): link atom placed at g × r_{QM-MM} from QM atom.
    pub scale: f64,
    /// Position of the link atom (Å).
    pub position: [f64; 3],
}

impl LinkAtom {
    /// Create a new link atom.
    pub fn new(qm_atom: usize, mm_atom: usize, scale: f64) -> Self {
        Self {
            qm_atom,
            mm_atom,
            scale,
            position: [0.0; 3],
        }
    }
}

// ---------------------------------------------------------------------------
// Standalone functions
// ---------------------------------------------------------------------------

/// Insert/update a hydrogen capping link atom at a QM/MM boundary bond.
///
/// The link atom H is placed along the Q1–M1 bond at:
/// r_H = r_Q1 + g × (r_M1 − r_Q1)
///
/// where g ≈ r_{Q-H} / r_{Q-M} ≈ 0.709 for C–C bonds.
///
/// # Arguments
/// * `pos_qm` — position of the QM boundary atom Q1 (Å)
/// * `pos_mm` — position of the MM boundary atom M1 (Å)
/// * `scale` — scale factor g
///
/// # Returns
/// Position of the link hydrogen atom.
pub fn link_atom_correction(pos_qm: [f64; 3], pos_mm: [f64; 3], scale: f64) -> [f64; 3] {
    [
        pos_qm[0] + scale * (pos_mm[0] - pos_qm[0]),
        pos_qm[1] + scale * (pos_mm[1] - pos_qm[1]),
        pos_qm[2] + scale * (pos_mm[2] - pos_qm[2]),
    ]
}

/// Compute electrostatic embedding energy: Coulomb interaction between
/// QM atom charges and MM point charges.
///
/// E_elec = Σ_{i ∈ QM} Σ_{j ∈ MM} q_i q_j / r_ij
///
/// (in units consistent with input; no dielectric constant applied here)
///
/// # Arguments
/// * `qm_positions` — positions of QM atoms
/// * `qm_charges` — partial charges of QM atoms
/// * `mm_positions` — positions of MM atoms
/// * `mm_charges` — partial charges of MM atoms
pub fn electrostatic_embedding(
    qm_positions: &[[f64; 3]],
    qm_charges: &[f64],
    mm_positions: &[[f64; 3]],
    mm_charges: &[f64],
) -> f64 {
    let n_qm = qm_positions.len().min(qm_charges.len());
    let n_mm = mm_positions.len().min(mm_charges.len());
    let mut energy = 0.0;
    for i in 0..n_qm {
        for j in 0..n_mm {
            let dx = qm_positions[i][0] - mm_positions[j][0];
            let dy = qm_positions[i][1] - mm_positions[j][1];
            let dz = qm_positions[i][2] - mm_positions[j][2];
            let r = (dx * dx + dy * dy + dz * dz).sqrt();
            if r < 1e-10 {
                continue;
            }
            energy += qm_charges[i] * mm_charges[j] / r;
        }
    }
    energy
}

/// Mechanical embedding: combine QM and MM forces at the boundary by
/// simple linear interpolation.
///
/// The QM/MM interaction energy is estimated as the Lennard-Jones + Coulomb
/// interaction between boundary QM atoms and all MM atoms using MM parameters.
///
/// # Arguments
/// * `qm_positions` — QM atom positions
/// * `qm_charges` — QM charges (for Coulomb)
/// * `mm_positions` — MM atom positions
/// * `mm_charges` — MM charges
/// * `epsilon_mm` — LJ epsilon for MM atoms
/// * `sigma_mm` — LJ sigma for MM atoms
pub fn mechanical_embedding(
    qm_positions: &[[f64; 3]],
    qm_charges: &[f64],
    mm_positions: &[[f64; 3]],
    mm_charges: &[f64],
    epsilon_mm: &[f64],
    sigma_mm: &[f64],
) -> f64 {
    let n_qm = qm_positions.len().min(qm_charges.len());
    let n_mm = mm_positions
        .len()
        .min(mm_charges.len())
        .min(epsilon_mm.len())
        .min(sigma_mm.len());
    let mut energy = 0.0;

    for i in 0..n_qm {
        for j in 0..n_mm {
            let dx = qm_positions[i][0] - mm_positions[j][0];
            let dy = qm_positions[i][1] - mm_positions[j][1];
            let dz = qm_positions[i][2] - mm_positions[j][2];
            let r2 = dx * dx + dy * dy + dz * dz;
            if r2 < 1e-10 {
                continue;
            }
            let r = r2.sqrt();
            // Coulomb
            energy += qm_charges[i] * mm_charges[j] / r;
            // LJ with QM atom treated as Lennard-Jones particle (sigma=1, eps=0.1)
            let eps = (0.1_f64 * epsilon_mm[j]).sqrt();
            let sig = 0.5 * (1.0 + sigma_mm[j]);
            let sr2 = (sig * sig) / r2;
            let sr6 = sr2 * sr2 * sr2;
            let sr12 = sr6 * sr6;
            energy += 4.0 * eps * (sr12 - sr6);
        }
    }
    energy
}

/// Total QM/MM energy: E_total = E_QM + E_MM + E_QM/MM.
///
/// # Arguments
/// * `e_qm` — QM energy
/// * `e_mm` — MM energy
/// * `e_qm_mm` — QM/MM coupling energy
pub fn qm_mm_energy(e_qm: f64, e_mm: f64, e_qm_mm: f64) -> f64 {
    e_qm + e_mm + e_qm_mm
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-12;

    // ── QmMethod ──────────────────────────────────────────────────────────

    #[test]
    fn test_qm_method_name() {
        assert_eq!(QmMethod::Pm3.name(), "PM3");
        assert_eq!(QmMethod::SccDftb.name(), "SCC-DFTB");
        assert_eq!(QmMethod::Hf.name(), "HF");
        assert_eq!(QmMethod::DftB3lyp.name(), "DFT/B3LYP");
    }

    #[test]
    fn test_qm_method_eq() {
        assert_eq!(QmMethod::Pm3, QmMethod::Pm3);
        assert_ne!(QmMethod::Pm3, QmMethod::Hf);
    }

    // ── QmRegion ──────────────────────────────────────────────────────────

    #[test]
    fn test_qm_region_new() {
        let qm = QmRegion::new(vec![0, 1, 2], QmMethod::Pm3, 0, 1);
        assert_eq!(qm.n_atoms(), 3);
        assert_eq!(qm.charge, 0);
        assert_eq!(qm.multiplicity, 1);
    }

    #[test]
    fn test_qm_region_placeholder_energy_zero_at_origin() {
        let mut qm = QmRegion::new(vec![0], QmMethod::Pm3, 0, 1);
        qm.set_positions(vec![[0.0, 0.0, 0.0]]);
        qm.run_placeholder();
        assert!(qm.energy.abs() < EPS);
    }

    #[test]
    fn test_qm_region_placeholder_energy_positive() {
        let mut qm = QmRegion::new(vec![0], QmMethod::Pm3, 0, 1);
        qm.set_positions(vec![[1.0, 0.0, 0.0]]);
        qm.run_placeholder();
        assert!((qm.energy - 1.0).abs() < EPS);
    }

    #[test]
    fn test_qm_region_forces_direction() {
        // Forces should point toward origin (restoring force)
        let mut qm = QmRegion::new(vec![0], QmMethod::Pm3, 0, 1);
        qm.set_positions(vec![[1.0, 0.0, 0.0]]);
        qm.run_placeholder();
        assert!(qm.forces[0][0] < 0.0);
    }

    // ── MmRegion ──────────────────────────────────────────────────────────

    #[test]
    fn test_mm_region_new() {
        let mm = MmRegion::new(vec![3, 4, 5]);
        assert_eq!(mm.n_atoms(), 3);
    }

    #[test]
    fn test_mm_region_placeholder_two_atoms() {
        let mut mm = MmRegion::new(vec![0, 1]);
        mm.positions = vec![[0.0, 0.0, 0.0], [4.0, 0.0, 0.0]];
        mm.run_placeholder();
        // LJ at r=4Å with sigma=3.5, eps=0.1: should be a finite value
        assert!(mm.energy.is_finite());
    }

    #[test]
    fn test_mm_region_placeholder_single_atom_zero() {
        let mut mm = MmRegion::new(vec![0]);
        mm.positions = vec![[0.0, 0.0, 0.0]];
        mm.run_placeholder();
        assert!(mm.energy.abs() < EPS);
    }

    // ── Link atom correction ──────────────────────────────────────────────

    #[test]
    fn test_link_atom_scale_zero() {
        // scale=0 → link atom at QM atom position
        let pos_qm = [1.0, 2.0, 3.0];
        let pos_mm = [5.0, 6.0, 7.0];
        let la = link_atom_correction(pos_qm, pos_mm, 0.0);
        for k in 0..3 {
            assert!((la[k] - pos_qm[k]).abs() < EPS);
        }
    }

    #[test]
    fn test_link_atom_scale_one() {
        // scale=1 → link atom at MM atom position
        let pos_qm = [1.0, 0.0, 0.0];
        let pos_mm = [3.0, 0.0, 0.0];
        let la = link_atom_correction(pos_qm, pos_mm, 1.0);
        assert!((la[0] - 3.0).abs() < EPS);
    }

    #[test]
    fn test_link_atom_scale_half() {
        // scale=0.5 → midpoint
        let pos_qm = [0.0, 0.0, 0.0];
        let pos_mm = [2.0, 0.0, 0.0];
        let la = link_atom_correction(pos_qm, pos_mm, 0.5);
        assert!((la[0] - 1.0).abs() < EPS);
    }

    #[test]
    fn test_link_atom_typical_scale() {
        // Typical C-H / C-C scale ≈ 0.709
        let pos_qm = [0.0, 0.0, 0.0];
        let pos_mm = [1.0, 0.0, 0.0];
        let la = link_atom_correction(pos_qm, pos_mm, 0.709);
        assert!((la[0] - 0.709).abs() < EPS);
    }

    // ── Electrostatic embedding ───────────────────────────────────────────

    #[test]
    fn test_electrostatic_embedding_single_pair() {
        // Two unit charges at distance 1.0 → E = 1.0
        let qm_pos = vec![[0.0, 0.0, 0.0]];
        let qm_chg = vec![1.0];
        let mm_pos = vec![[1.0, 0.0, 0.0]];
        let mm_chg = vec![1.0];
        let e = electrostatic_embedding(&qm_pos, &qm_chg, &mm_pos, &mm_chg);
        assert!((e - 1.0).abs() < EPS);
    }

    #[test]
    fn test_electrostatic_embedding_opposite_charges() {
        let qm_pos = vec![[0.0, 0.0, 0.0]];
        let qm_chg = vec![1.0];
        let mm_pos = vec![[1.0, 0.0, 0.0]];
        let mm_chg = vec![-1.0];
        let e = electrostatic_embedding(&qm_pos, &qm_chg, &mm_pos, &mm_chg);
        assert!((e + 1.0).abs() < EPS);
    }

    #[test]
    fn test_electrostatic_embedding_zero_charges() {
        let qm_pos = vec![[0.0, 0.0, 0.0]];
        let qm_chg = vec![0.0];
        let mm_pos = vec![[1.0, 0.0, 0.0]];
        let mm_chg = vec![1.0];
        let e = electrostatic_embedding(&qm_pos, &qm_chg, &mm_pos, &mm_chg);
        assert!(e.abs() < EPS);
    }

    #[test]
    fn test_electrostatic_embedding_inverse_distance() {
        let qm_pos = vec![[0.0, 0.0, 0.0]];
        let qm_chg = vec![1.0];
        let mm_pos1 = vec![[1.0, 0.0, 0.0]];
        let mm_pos2 = vec![[2.0, 0.0, 0.0]];
        let mm_chg = vec![1.0];
        let e1 = electrostatic_embedding(&qm_pos, &qm_chg, &mm_pos1, &mm_chg);
        let e2 = electrostatic_embedding(&qm_pos, &qm_chg, &mm_pos2, &mm_chg);
        assert!((e2 * 2.0 - e1).abs() < EPS);
    }

    // ── Mechanical embedding ──────────────────────────────────────────────

    #[test]
    fn test_mechanical_embedding_single_pair() {
        let qm_pos = vec![[0.0, 0.0, 0.0]];
        let qm_chg = vec![0.0];
        let mm_pos = vec![[10.0, 0.0, 0.0]]; // far away, LJ → 0
        let mm_chg = vec![0.0];
        let eps = vec![0.1];
        let sig = vec![3.5];
        let e = mechanical_embedding(&qm_pos, &qm_chg, &mm_pos, &mm_chg, &eps, &sig);
        // At r=10Å, LJ is tiny and Coulomb is zero
        assert!(e.abs() < 0.01);
    }

    #[test]
    fn test_mechanical_embedding_finite() {
        let qm_pos = vec![[0.0, 0.0, 0.0]];
        let qm_chg = vec![0.5];
        let mm_pos = vec![[5.0, 0.0, 0.0]];
        let mm_chg = vec![-0.5];
        let eps = vec![0.1];
        let sig = vec![3.5];
        let e = mechanical_embedding(&qm_pos, &qm_chg, &mm_pos, &mm_chg, &eps, &sig);
        assert!(e.is_finite());
    }

    // ── qm_mm_energy ──────────────────────────────────────────────────────

    #[test]
    fn test_qm_mm_energy_additive() {
        let e = qm_mm_energy(10.0, 5.0, 2.0);
        assert!((e - 17.0).abs() < EPS);
    }

    #[test]
    fn test_qm_mm_energy_zero() {
        assert!(qm_mm_energy(0.0, 0.0, 0.0).abs() < EPS);
    }

    #[test]
    fn test_qm_mm_energy_negative_interaction() {
        let e = qm_mm_energy(10.0, 5.0, -3.0);
        assert!((e - 12.0).abs() < EPS);
    }

    // ── QmMmSystem ────────────────────────────────────────────────────────

    #[test]
    fn test_qm_mm_system_total_energy() {
        let mut qm = QmRegion::new(vec![0], QmMethod::Pm3, 0, 1);
        qm.energy = 5.0;
        let mut mm = MmRegion::new(vec![1]);
        mm.energy = 3.0;
        let mut sys = QmMmSystem::new(qm, mm);
        sys.energy_qm_mm = 1.0;
        assert!((sys.total_energy() - 9.0).abs() < EPS);
    }

    #[test]
    fn test_qm_mm_system_add_link_atom() {
        let qm = QmRegion::new(vec![0], QmMethod::Pm3, 0, 1);
        let mm = MmRegion::new(vec![1]);
        let mut sys = QmMmSystem::new(qm, mm);
        sys.add_link_atom(LinkAtom::new(0, 0, 0.709));
        assert_eq!(sys.link_atoms.len(), 1);
    }

    #[test]
    fn test_link_atom_new() {
        let la = LinkAtom::new(2, 5, 0.72);
        assert_eq!(la.qm_atom, 2);
        assert_eq!(la.mm_atom, 5);
        assert!((la.scale - 0.72).abs() < EPS);
    }

    #[test]
    fn test_qm_region_charged() {
        let qm = QmRegion::new(vec![0, 1], QmMethod::SccDftb, -1, 2);
        assert_eq!(qm.charge, -1);
        assert_eq!(qm.multiplicity, 2);
    }
}
