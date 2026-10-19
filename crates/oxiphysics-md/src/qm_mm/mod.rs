// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! QM/MM hybrid quantum mechanics / molecular mechanics methods.
//!
//! Provides QM/MM hybrid calculations where the QM region uses PM3 (semi-empirical),
//! SCC-DFTB (density-functional tight-binding), HF/STO-3G (Hartree-Fock), or
//! B3LYP/STO-3G (restricted Kohn-Sham DFT).  The MM region uses standard force-field point charges.
//!
//! ## Energy decomposition
//!
//! E_total = E_QM + E_MM + E_QM/MM
//!
//! where E_QM/MM can be computed via:
//! - **Mechanical embedding**: classical Coulomb + Lennard-Jones interaction between
//!   fixed QM and MM charges ([`mechanical_embedding`]).
//! - **Electrostatic embedding**: the MM point charges enter the QM Hamiltonian
//!   ([`QmRegion::external_charges`]) and polarise the QM wavefunction of every
//!   engine (PM3, SCC-DFTB, HF, KS-DFT).  [`QmMmSystem::compute`] fills the external
//!   charges from the MM region, runs the QM engine and evaluates E_QM/MM together
//!   with the QM/MM forces on both regions.
//!
//! ## Boundary treatment
//!
//! Covalent bonds that cross the QM/MM boundary are capped with hydrogen
//! link atoms to saturate the QM region valence.  The charge of the MM boundary
//! atom bonded to a link atom is deleted from the embedding (charge-deletion
//! scheme) and the link-atom force is redistributed onto the QM and MM boundary
//! atoms by the chain rule.
//!
//! References:
//! - Field, M. J., Bash, P. A., & Karplus, M. (1990). *J. Comput. Chem.* 11, 700.
//! - Bakowies, D., & Thiel, W. (1996). *J. Phys. Chem.* 100, 10580.

pub mod dft;
pub(crate) mod embedding;
pub mod error;
pub mod hf;
pub mod pm3;
pub mod scc_dftb;
pub(crate) mod sp_overlap;

pub use error::QmError;

/// Coulomb constant in kcal·Å/(mol·e²): E = 332.0637 · q_i q_j / r.
pub const COULOMB_KCAL_PER_MOL_ANG: f64 = 332.063_713;

// ---------------------------------------------------------------------------
// QmMethod — enum for QM level of theory
// ---------------------------------------------------------------------------

/// QM level of theory.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum QmMethod {
    /// PM3 semi-empirical NDDO method (Stewart 1989).
    Pm3,
    /// Self-consistent-charge DFTB (Elstner 1998, mio-0-1 parameters).
    SccDftb,
    /// Restricted Hartree-Fock with STO-3G minimal basis.
    Hf,
    /// Restricted Kohn-Sham DFT with the B3LYP hybrid functional (STO-3G basis).
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
///
/// ## Energy decomposition with external charges
///
/// When [`external_charges`](Self::external_charges) is non-empty (electrostatic
/// embedding) the QM energy is split as
///
/// E = `energy` + `embedding_energy`
///
/// where `embedding_energy` is the electrostatic interaction of the polarised QM
/// electrons and cores with the external charges and `energy` is the remaining
/// (polarised) QM energy.  The reported `forces` and `external_charge_forces`
/// are the negative gradients of the sum.  Without external charges
/// `embedding_energy` is zero and `energy` is the plain QM energy.
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
    /// External point charges `(position in Å, charge in e)` that polarise the
    /// QM wavefunction (electrostatic embedding).  Empty for a gas-phase run.
    pub external_charges: Vec<([f64; 3], f64)>,
    /// QM energy (in kcal/mol, filled after a QM calculation), excluding the
    /// electrostatic interaction with `external_charges` (see `embedding_energy`).
    pub energy: f64,
    /// Electrostatic interaction energy (kcal/mol) between the polarised QM region
    /// (electrons and cores) and `external_charges`.
    pub embedding_energy: f64,
    /// QM forces on QM atoms \[\[fx, fy, fz\\]; n_qm] in kcal/(mol·Å):
    /// −∇(energy + embedding_energy).
    pub forces: Vec<[f64; 3]>,
    /// Mulliken net atomic charges (e) of the QM atoms after the last run.
    pub partial_charges: Vec<f64>,
    /// Forces (kcal/(mol·Å)) exerted by the QM region on each external charge.
    pub external_charge_forces: Vec<[f64; 3]>,
    /// Lennard-Jones well depth (kcal/mol) of each QM atom for QM/MM van der Waals
    /// interactions in [`QmMmSystem::compute`]; 0 disables the term for that atom.
    pub lj_epsilon: Vec<f64>,
    /// Lennard-Jones diameter (Å) of each QM atom for QM/MM van der Waals interactions.
    pub lj_sigma: Vec<f64>,
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
            external_charges: Vec::new(),
            energy: 0.0,
            embedding_energy: 0.0,
            forces: vec![[0.0; 3]; n],
            partial_charges: vec![0.0; n],
            external_charge_forces: Vec::new(),
            lj_epsilon: vec![0.0; n],
            lj_sigma: vec![0.0; n],
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

    /// Total potential energy of the QM region in the field of the external
    /// charges (kcal/mol): `energy + embedding_energy`.
    pub fn total_energy_in_field(&self) -> f64 {
        self.energy + self.embedding_energy
    }

    /// Dipole moment Σ_A q_A R_A (e·Å) of the Mulliken partial charges from the
    /// last run.  Origin independent only for a neutral region.
    pub fn dipole_moment(&self) -> [f64; 3] {
        let mut mu = [0.0_f64; 3];
        for (q, r) in self.partial_charges.iter().zip(self.positions.iter()) {
            for (m, x) in mu.iter_mut().zip(r.iter()) {
                *m += q * x;
            }
        }
        mu
    }

    /// Run the QM calculation using the configured method.
    ///
    /// Populates `energy`, `embedding_energy` (kcal/mol), `forces`,
    /// `external_charge_forces` (kcal/(mol·Å)) and `partial_charges` (e).
    /// All outputs are cleared before the calculation starts, so a failed run
    /// never leaves the results of a previous call behind.
    ///
    /// # Errors
    /// Returns a [`QmError`] if the region is empty or inconsistent, contains an
    /// element the method does not support, is open-shell, has an invalid
    /// electron count or external charge, or if the SCF fails to converge.
    pub fn run(&mut self) -> Result<(), QmError> {
        match self.method {
            QmMethod::Pm3 => pm3::run_pm3(self),
            QmMethod::SccDftb => scc_dftb::run_scc_dftb(self),
            QmMethod::Hf => hf::run_hf(self),
            QmMethod::DftB3lyp => dft::run_dft_b3lyp(self),
        }
    }

    /// Former placeholder QM calculation; now simply forwards to [`run`](QmRegion::run).
    ///
    /// # Deprecated
    /// Use [`run`](QmRegion::run) instead, which dispatches to the real QM engine.
    ///
    /// # Errors
    /// Same as [`run`](QmRegion::run).
    #[deprecated(since = "0.1.2", note = "Use `run()` which calls the real QM solver")]
    pub fn run_placeholder(&mut self) -> Result<(), QmError> {
        self.run()
    }
}

// ---------------------------------------------------------------------------
// MmRegion
// ---------------------------------------------------------------------------

/// The classical molecular mechanics region.
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
///
/// [`compute`](Self::compute) evaluates the QM region with electrostatic
/// embedding in the MM point charges and fills the QM/MM coupling terms:
///
/// E_total = E_QM + E_MM + E_QM/MM,  E_QM/MM = E_elec + E_vdW
///
/// with E_QM = `qm.energy` (polarised QM energy), E_elec = `qm.embedding_energy`
/// and E_vdW the QM–MM Lennard-Jones energy.  E_MM (`mm.energy`) and the
/// MM–MM forces are the responsibility of the MM force field.
pub struct QmMmSystem {
    /// The QM region.
    pub qm: QmRegion,
    /// The MM region.
    pub mm: MmRegion,
    /// QM/MM interaction energy (kcal/mol): electrostatic embedding + van der Waals.
    pub energy_qm_mm: f64,
    /// QM–MM Lennard-Jones (van der Waals) part of `energy_qm_mm` (kcal/mol).
    pub energy_qm_mm_vdw: f64,
    /// Link atoms inserted at QM/MM boundaries.
    pub link_atoms: Vec<LinkAtom>,
    /// Total forces on the QM atoms (kcal/(mol·Å)): QM (including embedding and
    /// redistributed link-atom forces) plus QM–MM van der Waals.
    pub qm_total_forces: Vec<[f64; 3]>,
    /// QM/MM coupling forces on the MM atoms (kcal/(mol·Å)): electrostatic
    /// embedding, QM–MM van der Waals and redistributed link-atom forces.  Add
    /// these to the MM force-field forces.
    pub mm_coupling_forces: Vec<[f64; 3]>,
}

impl QmMmSystem {
    /// Create a new QM/MM system from QM and MM regions.
    pub fn new(qm: QmRegion, mm: MmRegion) -> Self {
        let n_qm = qm.atomic_numbers.len();
        let n_mm = mm.positions.len();
        Self {
            qm,
            mm,
            energy_qm_mm: 0.0,
            energy_qm_mm_vdw: 0.0,
            link_atoms: Vec::new(),
            qm_total_forces: vec![[0.0; 3]; n_qm],
            mm_coupling_forces: vec![[0.0; 3]; n_mm],
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

    /// Evaluate the QM region in the field of the MM point charges and the QM/MM
    /// coupling energy and forces.
    ///
    /// 1. Every link atom is placed at r_Q + g (r_M − r_Q) and added to the QM
    ///    calculation as a hydrogen atom.
    /// 2. The MM charges (except those of MM boundary atoms bonded to a link atom)
    ///    become `qm.external_charges`.
    /// 3. The QM engine is run; `qm.energy`, `qm.embedding_energy`, `qm.forces`,
    ///    `qm.partial_charges` and `qm.external_charge_forces` are filled.  The
    ///    link-atom force F_L is redistributed as (1 − g) F_L onto the QM and
    ///    g F_L onto the MM boundary atom.
    /// 4. The QM–MM Lennard-Jones energy is evaluated with Lorentz–Berthelot
    ///    combination of `qm.lj_epsilon/lj_sigma` and `mm.epsilon/sigma`, excluding
    ///    the bonded boundary pairs.
    /// 5. `energy_qm_mm = qm.embedding_energy + energy_qm_mm_vdw`.
    ///
    /// # Errors
    /// Returns a [`QmError`] if the input arrays are inconsistent, a link atom
    /// references a non-existent atom, or the QM calculation fails.  All outputs
    /// are cleared first, so a failure never leaves stale values behind.
    pub fn compute(&mut self) -> Result<(), QmError> {
        let n_qm = self.qm.atomic_numbers.len();
        let n_mm = self.mm.positions.len();
        self.energy_qm_mm = 0.0;
        self.energy_qm_mm_vdw = 0.0;
        self.qm_total_forces = vec![[0.0; 3]; n_qm];
        self.mm_coupling_forces = vec![[0.0; 3]; n_mm];
        self.validate_inputs(n_qm, n_mm)?;

        // 1. Link atoms.
        let mut boundary_mm = vec![false; n_mm];
        for la in &mut self.link_atoms {
            la.position = link_atom_correction(
                self.qm.positions[la.qm_atom],
                self.mm.positions[la.mm_atom],
                la.scale,
            );
            boundary_mm[la.mm_atom] = true;
        }

        // 2. External charges from the MM region (charge deletion at the boundary).
        let embedded_mm: Vec<usize> = (0..n_mm).filter(|&j| !boundary_mm[j]).collect();
        let external: Vec<([f64; 3], f64)> = embedded_mm
            .iter()
            .map(|&j| (self.mm.positions[j], self.mm.charges[j]))
            .collect();

        // 3. QM calculation on the link-atom-capped region.
        let mut work_indices = self.qm.atom_indices.clone();
        let mut work_numbers = self.qm.atomic_numbers.clone();
        let mut work_positions = self.qm.positions.clone();
        for la in &self.link_atoms {
            work_indices.push(
                self.mm
                    .atom_indices
                    .get(la.mm_atom)
                    .copied()
                    .unwrap_or(la.mm_atom),
            );
            work_numbers.push(1);
            work_positions.push(la.position);
        }
        let mut work = QmRegion::new(
            work_indices,
            self.qm.method,
            self.qm.charge,
            self.qm.multiplicity,
        );
        work.atomic_numbers = work_numbers;
        work.positions = work_positions;
        work.external_charges = external.clone();
        let run_result = work.run();

        self.qm.external_charges = external;
        self.qm.energy = work.energy;
        self.qm.embedding_energy = work.embedding_energy;
        self.qm.forces = work.forces.iter().take(n_qm).copied().collect();
        self.qm.forces.resize(n_qm, [0.0; 3]);
        self.qm.partial_charges = work.partial_charges.iter().take(n_qm).copied().collect();
        self.qm.partial_charges.resize(n_qm, 0.0);
        self.qm.external_charge_forces = work.external_charge_forces.clone();
        run_result?;

        for (k, &j) in embedded_mm.iter().enumerate() {
            if let Some(f) = work.external_charge_forces.get(k) {
                add3(&mut self.mm_coupling_forces[j], f, 1.0);
            }
        }
        for (l, la) in self.link_atoms.iter().enumerate() {
            if let Some(f_link) = work.forces.get(n_qm + l) {
                add3(&mut self.qm.forces[la.qm_atom], f_link, 1.0 - la.scale);
                add3(&mut self.mm_coupling_forces[la.mm_atom], f_link, la.scale);
            }
        }
        self.qm_total_forces = self.qm.forces.clone();

        // 4. QM–MM van der Waals.
        let mut e_vdw = 0.0_f64;
        for a in 0..n_qm {
            for j in 0..n_mm {
                if self
                    .link_atoms
                    .iter()
                    .any(|la| la.qm_atom == a && la.mm_atom == j)
                {
                    continue;
                }
                let eps = (self.qm.lj_epsilon[a] * self.mm.epsilon[j]).sqrt();
                if eps <= 0.0 {
                    continue;
                }
                let sigma = 0.5 * (self.qm.lj_sigma[a] + self.mm.sigma[j]);
                let d = [
                    self.qm.positions[a][0] - self.mm.positions[j][0],
                    self.qm.positions[a][1] - self.mm.positions[j][1],
                    self.qm.positions[a][2] - self.mm.positions[j][2],
                ];
                let r2 = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
                if r2 < 1e-20 {
                    return Err(QmError::InvalidExternalCharge { index: j });
                }
                let sr2 = sigma * sigma / r2;
                let sr6 = sr2 * sr2 * sr2;
                let sr12 = sr6 * sr6;
                e_vdw += 4.0 * eps * (sr12 - sr6);
                // F_a = −dE/dr_a = 24 ε (2 sr¹² − sr⁶) / r² · (r_a − r_j)
                let f_scale = 24.0 * eps * (2.0 * sr12 - sr6) / r2;
                add3(&mut self.qm_total_forces[a], &d, f_scale);
                add3(&mut self.mm_coupling_forces[j], &d, -f_scale);
            }
        }
        self.energy_qm_mm_vdw = e_vdw;

        // 5. Interaction energy.
        self.energy_qm_mm = self.qm.embedding_energy + e_vdw;
        Ok(())
    }

    /// Check the consistency of all per-atom arrays and link-atom indices.
    fn validate_inputs(&self, n_qm: usize, n_mm: usize) -> Result<(), QmError> {
        let checks: [(&'static str, usize, usize); 6] = [
            ("QM positions", n_qm, self.qm.positions.len()),
            ("QM lj_epsilon", n_qm, self.qm.lj_epsilon.len()),
            ("QM lj_sigma", n_qm, self.qm.lj_sigma.len()),
            ("MM charges", n_mm, self.mm.charges.len()),
            ("MM epsilon", n_mm, self.mm.epsilon.len()),
            ("MM sigma", n_mm, self.mm.sigma.len()),
        ];
        for (what, expected, found) in checks {
            if expected != found {
                return Err(QmError::DimensionMismatch {
                    what,
                    expected,
                    found,
                });
            }
        }
        for la in &self.link_atoms {
            if la.qm_atom >= n_qm {
                return Err(QmError::InvalidLinkAtom {
                    region: "QM",
                    index: la.qm_atom,
                    len: n_qm,
                });
            }
            if la.mm_atom >= n_mm {
                return Err(QmError::InvalidLinkAtom {
                    region: "MM",
                    index: la.mm_atom,
                    len: n_mm,
                });
            }
        }
        Ok(())
    }
}

/// `target += scale · v`.
fn add3(target: &mut [f64; 3], v: &[f64; 3], scale: f64) {
    for (t, x) in target.iter_mut().zip(v.iter()) {
        *t += scale * x;
    }
}

// ---------------------------------------------------------------------------
// LinkAtom
// ---------------------------------------------------------------------------

/// A hydrogen capping link atom at a QM/MM boundary bond.
///
/// The link atom (always hydrogen) is placed along the bond between a
/// QM atom and a MM atom, at a fixed distance g × r_{QM-MM}.
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

/// Classical point-charge Coulomb interaction between fixed QM atom charges and
/// MM point charges.
///
/// E_elec = Σ_{i ∈ QM} Σ_{j ∈ MM} q_i q_j / r_ij
///
/// The result is in e²/length of the input (no Coulomb constant or dielectric
/// applied); multiply by [`COULOMB_KCAL_PER_MOL_ANG`] for kcal/mol with Å and e.
/// This does **not** polarise the QM wavefunction — for true electrostatic
/// embedding set [`QmRegion::external_charges`] or use [`QmMmSystem::compute`].
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
    fn test_qm_region_run_pm3_h2() {
        // PM3: H2 at 0.74 Å is a bound molecule.
        let mut qm = QmRegion::new(vec![0, 1], QmMethod::Pm3, 0, 1);
        qm.atomic_numbers = vec![1, 1];
        qm.set_positions(vec![[0.0, 0.0, 0.0], [0.74, 0.0, 0.0]]);
        qm.run().expect("PM3 H2");
        assert!(
            qm.energy.is_finite() && qm.energy < 0.0,
            "PM3 H2: {}",
            qm.energy
        );
    }

    #[test]
    fn test_qm_region_run_hf_single_h_is_open_shell_error() {
        // A single H atom has one electron: the closed-shell HF engine must
        // report it instead of returning a meaningless energy.
        let mut qm = QmRegion::new(vec![0], QmMethod::Hf, 0, 1);
        qm.atomic_numbers = vec![1];
        qm.set_positions(vec![[0.0, 0.0, 0.0]]);
        assert_eq!(qm.run(), Err(QmError::OpenShellNotSupported));
        assert_eq!(qm.energy, 0.0);
    }

    #[test]
    fn test_qm_region_run_dftb_h2() {
        let mut qm = QmRegion::new(vec![0, 1], QmMethod::SccDftb, 0, 1);
        qm.atomic_numbers = vec![1, 1];
        qm.set_positions(vec![[0.0, 0.0, 0.0], [0.74, 0.0, 0.0]]);
        qm.run().expect("DFTB H2");
        assert!(
            qm.energy.is_finite() && qm.energy < 0.0,
            "DFTB H2: {}",
            qm.energy
        );
    }

    #[test]
    fn test_qm_region_run_unsupported_element_is_error_for_every_method() {
        for method in [
            QmMethod::Pm3,
            QmMethod::SccDftb,
            QmMethod::Hf,
            QmMethod::DftB3lyp,
        ] {
            let mut qm = QmRegion::new(vec![0, 1], method, 0, 1);
            qm.atomic_numbers = vec![92, 92];
            qm.set_positions(vec![[0.0, 0.0, 0.0], [2.5, 0.0, 0.0]]);
            qm.energy = 42.0;
            qm.forces = vec![[1.0; 3]; 2];
            assert_eq!(qm.run(), Err(QmError::UnsupportedElement(92)), "{method:?}");
            assert_eq!(qm.energy, 0.0, "{method:?} left a stale energy");
            assert!(
                qm.forces.iter().all(|f| *f == [0.0; 3]),
                "{method:?} stale forces"
            );
        }
    }

    #[test]
    fn test_qm_region_empty_and_inconsistent_inputs() {
        let mut empty = QmRegion::new(vec![], QmMethod::Pm3, 0, 1);
        assert_eq!(empty.run(), Err(QmError::EmptyRegion));
        let mut bad = QmRegion::new(vec![0, 1], QmMethod::SccDftb, 0, 1);
        bad.positions = vec![[0.0; 3]];
        assert!(matches!(bad.run(), Err(QmError::DimensionMismatch { .. })));
        let mut on_top = QmRegion::new(vec![0, 1], QmMethod::Pm3, 0, 1);
        on_top.set_positions(vec![[0.0; 3], [0.74, 0.0, 0.0]]);
        on_top.external_charges = vec![([0.74, 0.0, 0.0], 0.5)];
        assert_eq!(
            on_top.run(),
            Err(QmError::InvalidExternalCharge { index: 0 })
        );
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

    /// QM water (SCC-DFTB) next to an MM TIP3P-like water.
    fn water_dimer_system() -> QmMmSystem {
        let mut qm = QmRegion::new(vec![0, 1, 2], QmMethod::SccDftb, 0, 1);
        qm.atomic_numbers = vec![8, 1, 1];
        qm.set_positions(vec![
            [0.0, 0.0, 0.0],
            [0.0, 0.757, 0.586],
            [0.0, -0.757, 0.586],
        ]);
        qm.lj_epsilon = vec![0.1521, 0.0, 0.0];
        qm.lj_sigma = vec![3.1507, 0.0, 0.0];
        let mut mm = MmRegion::new(vec![3, 4, 5]);
        mm.positions = vec![
            [0.2, 0.1, -2.9],
            [0.2, 0.857, -3.486],
            [0.2, -0.657, -3.486],
        ];
        mm.charges = vec![-0.834, 0.417, 0.417];
        mm.epsilon = vec![0.1521, 0.0, 0.0];
        mm.sigma = vec![3.1507, 0.0, 0.0];
        mm.energy = -1.25;
        QmMmSystem::new(qm, mm)
    }

    /// E_QM + E_QM/MM (the MM–MM energy is held fixed by these tests).
    fn coupled_energy(sys: &QmMmSystem) -> f64 {
        sys.qm.energy + sys.energy_qm_mm
    }

    #[test]
    fn test_qm_mm_compute_energy_decomposition() {
        let mut sys = water_dimer_system();
        sys.compute().expect("QM/MM water dimer");
        // MM charges were embedded and polarise the QM region.
        assert_eq!(sys.qm.external_charges.len(), 3);
        assert!(sys.qm.embedding_energy.abs() > 1e-3);
        assert!(sys.energy_qm_mm_vdw.abs() > 1e-6);
        assert!(
            (sys.energy_qm_mm - (sys.qm.embedding_energy + sys.energy_qm_mm_vdw)).abs() < 1e-12
        );
        // No double counting: the embedding energy is not part of qm.energy.
        assert!(
            (sys.total_energy() - (sys.qm.energy + sys.mm.energy + sys.energy_qm_mm)).abs() < 1e-12
        );
        let mut vacuum = water_dimer_system().qm;
        vacuum.run().expect("vacuum");
        // SCC-DFTB is variational, so polarising the QM water costs internal
        // energy, which the embedding interaction more than pays back.
        assert!(sys.qm.energy > vacuum.energy);
        // The dipole of the QM water changes in the field of the MM water.
        let mu_vac = vacuum.dipole_moment();
        let mu_emb = sys.qm.dipole_moment();
        let diff: f64 = (0..3).map(|k| (mu_vac[k] - mu_emb[k]).powi(2)).sum::<f64>();
        assert!(diff.sqrt() > 1e-3);
    }

    #[test]
    fn test_qm_mm_compute_forces_match_finite_differences() {
        let mut sys = water_dimer_system();
        sys.compute().expect("QM/MM water dimer");
        let h = 1e-3;
        // MM oxygen (embedding + van der Waals) and an MM hydrogen (embedding only).
        for (j, k) in [(0usize, 2usize), (0, 0), (1, 1)] {
            let mut e = [0.0; 2];
            for (slot, sign) in e.iter_mut().zip([1.0, -1.0]) {
                let mut d = water_dimer_system();
                d.mm.positions[j][k] += sign * h;
                d.compute().expect("displaced MM atom");
                *slot = coupled_energy(&d);
            }
            let fd = -(e[0] - e[1]) / (2.0 * h);
            let f = sys.mm_coupling_forces[j][k];
            assert!(
                (fd - f).abs() < 5e-3 + 1e-2 * fd.abs(),
                "MM atom {j} component {k}: {f} vs FD {fd}"
            );
        }
        // QM oxygen.
        for k in 0..3 {
            let mut e = [0.0; 2];
            for (slot, sign) in e.iter_mut().zip([1.0, -1.0]) {
                let mut d = water_dimer_system();
                d.qm.positions[0][k] += sign * h;
                d.compute().expect("displaced QM atom");
                *slot = coupled_energy(&d);
            }
            let fd = -(e[0] - e[1]) / (2.0 * h);
            let f = sys.qm_total_forces[0][k];
            assert!(
                (fd - f).abs() < 2e-2 + 1e-2 * fd.abs(),
                "QM O component {k}: {f} vs FD {fd}"
            );
        }
        // Newton's third law over both regions.
        for k in 0..3 {
            let total: f64 = sys.qm_total_forces.iter().map(|f| f[k]).sum::<f64>()
                + sys.mm_coupling_forces.iter().map(|f| f[k]).sum::<f64>();
            assert!(total.abs() < 2e-2, "net force {k}: {total}");
        }
    }

    /// QM methyl group (PM3) bonded to an MM methyl group through a link atom.
    fn ethane_link_system() -> QmMmSystem {
        let b = 1.09;
        let tet = [
            [1.0 / 3.0, 0.942_809, 0.0],
            [1.0 / 3.0, -0.471_405, 0.816_497],
            [1.0 / 3.0, -0.471_405, -0.816_497],
        ];
        let c_mm = [1.53, 0.0, 0.0];
        let mut qm = QmRegion::new(vec![0, 1, 2, 3], QmMethod::Pm3, 0, 1);
        qm.atomic_numbers = vec![6, 1, 1, 1];
        let mut qm_pos = vec![[0.0; 3]];
        for t in &tet {
            qm_pos.push([-b * t[0], b * t[1], b * t[2]]);
        }
        qm.set_positions(qm_pos);
        let mut mm = MmRegion::new(vec![4, 5, 6, 7]);
        let mut mm_pos = vec![c_mm];
        for t in &tet {
            mm_pos.push([c_mm[0] + b * t[0], b * t[1] + 0.05, b * t[2]]);
        }
        mm.positions = mm_pos;
        mm.charges = vec![-0.27, 0.09, 0.09, 0.09];
        mm.epsilon = vec![0.0; 4];
        mm.sigma = vec![0.0; 4];
        let mut sys = QmMmSystem::new(qm, mm);
        sys.add_link_atom(LinkAtom::new(0, 0, 0.709));
        sys
    }

    #[test]
    fn test_qm_mm_compute_with_link_atom() {
        let mut sys = ethane_link_system();
        sys.compute().expect("QM/MM ethane");
        // Link hydrogen placed along the C–C bond; boundary MM charge deleted.
        assert!((sys.link_atoms[0].position[0] - 0.709 * 1.53).abs() < 1e-12);
        assert_eq!(sys.qm.external_charges.len(), 3);
        assert_eq!(sys.qm.forces.len(), 4);
        // The MM boundary carbon only feels the redistributed link-atom force.
        let h = 1e-3;
        for k in 0..3 {
            let mut e = [0.0; 2];
            for (slot, sign) in e.iter_mut().zip([1.0, -1.0]) {
                let mut d = ethane_link_system();
                d.mm.positions[0][k] += sign * h;
                d.compute().expect("displaced boundary atom");
                *slot = coupled_energy(&d);
            }
            let fd = -(e[0] - e[1]) / (2.0 * h);
            let f = sys.mm_coupling_forces[0][k];
            assert!(
                (fd - f).abs() < 1e-2 + 1e-2 * fd.abs(),
                "boundary MM C component {k}: {f} vs FD {fd}"
            );
        }
        // QM boundary carbon: QM force + (1 − g) F_link.
        let mut e = [0.0; 2];
        for (slot, sign) in e.iter_mut().zip([1.0, -1.0]) {
            let mut d = ethane_link_system();
            d.qm.positions[0][0] += sign * h;
            d.compute().expect("displaced QM boundary atom");
            *slot = coupled_energy(&d);
        }
        let fd = -(e[0] - e[1]) / (2.0 * h);
        assert!(
            (fd - sys.qm_total_forces[0][0]).abs() < 2e-2 + 1e-2 * fd.abs(),
            "QM C: {} vs FD {fd}",
            sys.qm_total_forces[0][0]
        );
    }

    #[test]
    fn test_qm_mm_compute_rejects_bad_input_and_clears_outputs() {
        let mut sys = water_dimer_system();
        sys.compute().expect("valid system");
        assert!(sys.energy_qm_mm != 0.0);
        sys.add_link_atom(LinkAtom::new(0, 7, 0.7));
        assert_eq!(
            sys.compute(),
            Err(QmError::InvalidLinkAtom {
                region: "MM",
                index: 7,
                len: 3
            })
        );
        assert_eq!(sys.energy_qm_mm, 0.0);
        assert!(sys.mm_coupling_forces.iter().all(|f| *f == [0.0; 3]));

        let mut sys = water_dimer_system();
        sys.mm.charges.pop();
        assert!(matches!(
            sys.compute(),
            Err(QmError::DimensionMismatch { .. })
        ));

        let mut sys = water_dimer_system();
        sys.qm.atomic_numbers[0] = 54;
        assert_eq!(sys.compute(), Err(QmError::UnsupportedElement(54)));
        assert_eq!(sys.qm.energy, 0.0);
        assert_eq!(sys.energy_qm_mm, 0.0);
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
