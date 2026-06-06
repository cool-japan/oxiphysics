//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

/// A link atom used to cap the QM–MM boundary.
#[derive(Debug, Clone)]
pub struct LinkAtomQmmm {
    /// Position of the link atom \[Bohr\].
    pub position: [f64; 3],
    /// QM atom index it is bonded to.
    pub qm_atom: usize,
    /// MM atom index it replaces.
    pub mm_atom: usize,
    /// Scale factor for position interpolation (g = R_QM–MM / R_QM–LA).
    pub g_factor: f64,
}
impl LinkAtomQmmm {
    /// Create a link atom by interpolating between QM and MM atom positions.
    pub fn from_bond(
        pos_qm: [f64; 3],
        pos_mm: [f64; 3],
        qm_atom: usize,
        mm_atom: usize,
        g: f64,
    ) -> Self {
        let position = [
            pos_qm[0] + g * (pos_mm[0] - pos_qm[0]),
            pos_qm[1] + g * (pos_mm[1] - pos_qm[1]),
            pos_qm[2] + g * (pos_mm[2] - pos_qm[2]),
        ];
        Self {
            position,
            qm_atom,
            mm_atom,
            g_factor: g,
        }
    }
    /// Force correction on QM atom due to link atom force.
    pub fn qm_force_correction(&self, f_link: [f64; 3]) -> [f64; 3] {
        let g = self.g_factor;
        [
            (1.0 - g) * f_link[0],
            (1.0 - g) * f_link[1],
            (1.0 - g) * f_link[2],
        ]
    }
    /// Force correction on MM atom.
    pub fn mm_force_correction(&self, f_link: [f64; 3]) -> [f64; 3] {
        let g = self.g_factor;
        [g * f_link[0], g * f_link[1], g * f_link[2]]
    }
}
/// Convergence status of the SCF iterations.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ScfStatus {
    /// SCF converged within the threshold.
    Converged,
    /// SCF did not converge within the maximum number of iterations.
    NotConverged,
}
/// Kohn-Sham DFT result.
#[derive(Debug, Clone)]
pub struct DftResult {
    /// Total DFT energy \[Hartree\].
    pub energy: f64,
    /// Exchange-correlation energy \[Hartree\].
    pub xc_energy: f64,
    /// Kinetic energy \[Hartree\].
    pub kinetic_energy: f64,
    /// Kohn-Sham orbital energies \[Hartree\].
    pub orbital_energies: Vec<f64>,
    /// Converged density matrix.
    pub density_matrix: Vec<Vec<f64>>,
    /// SCF status.
    pub status: ScfStatus,
}
/// A molecular orbital described in the LCAO basis.
#[derive(Debug, Clone)]
pub struct MolecularOrbital {
    /// Orbital index (0-based).
    pub index: usize,
    /// Orbital energy \[Hartree\].
    pub energy: f64,
    /// Occupation number (0, 1, or 2).
    pub occupation: f64,
    /// LCAO coefficients C_{mu,i} for this MO.
    pub lcao_coefficients: Vec<f64>,
    /// Symmetry label (e.g., "a1", "b2").
    pub symmetry: String,
}
impl MolecularOrbital {
    /// Create a molecular orbital.
    pub fn new(index: usize, energy: f64, occupation: f64, lcao_coefficients: Vec<f64>) -> Self {
        Self {
            index,
            energy,
            occupation,
            lcao_coefficients,
            symmetry: String::from("a"),
        }
    }
    /// Bond order between atoms a and b (Mayer bond order approximation).
    /// `basis_a` and `basis_b` are the indices of AOs on atoms a and b.
    pub fn mayer_bond_order(
        density: &[Vec<f64>],
        overlap: &[Vec<f64>],
        basis_a: &[usize],
        basis_b: &[usize],
    ) -> f64 {
        let mut bo = 0.0;
        for &mu in basis_a {
            for &nu in basis_b {
                if mu < density.len()
                    && nu < density.len()
                    && mu < overlap.len()
                    && nu < overlap.len()
                {
                    let ps_mn = density[mu]
                        .iter()
                        .zip(overlap[mu].iter())
                        .map(|(p, s)| p * s)
                        .sum::<f64>();
                    let ps_nm = density[nu]
                        .iter()
                        .zip(overlap[nu].iter())
                        .map(|(p, s)| p * s)
                        .sum::<f64>();
                    bo += density[mu][nu] * overlap[mu][nu] * ps_mn * ps_nm;
                }
            }
        }
        bo
    }
    /// Wiberg bond index between atoms a and b.
    pub fn wiberg_bond_index(
        density: &[Vec<f64>],
        overlap: &[Vec<f64>],
        basis_a: &[usize],
        basis_b: &[usize],
    ) -> f64 {
        let mut wbi = 0.0;
        for &mu in basis_a {
            for &nu in basis_b {
                if mu < density.len() && nu < density[mu].len() {
                    let p_mn = density[mu][nu];
                    let s_mn = if mu < overlap.len() && nu < overlap[mu].len() {
                        overlap[mu][nu]
                    } else {
                        0.0
                    };
                    wbi += (p_mn * s_mn).powi(2);
                }
            }
        }
        wbi
    }
    /// Orbital participation ratio (delocalization measure).
    pub fn participation_ratio(&self) -> f64 {
        let sum_sq: f64 = self.lcao_coefficients.iter().map(|c| c * c).sum();
        let sum_4: f64 = self.lcao_coefficients.iter().map(|c| c * c * c * c).sum();
        if sum_4 < 1e-20 {
            return 0.0;
        }
        sum_sq * sum_sq / sum_4
    }
    /// Check if this is an occupied orbital.
    pub fn is_occupied(&self) -> bool {
        self.occupation > 0.5
    }
    /// Compute dipole moment contribution along axis from this MO.
    pub fn dipole_contribution(&self, basis_centers: &[[f64; 3]], axis: usize) -> f64 {
        let mut d = 0.0;
        for (mu, c) in self.lcao_coefficients.iter().enumerate() {
            if mu < basis_centers.len() {
                d += c * c * self.occupation * basis_centers[mu][axis];
            }
        }
        d
    }
}
/// Exchange-correlation functional type.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum XcFunctional {
    /// Local density approximation (Slater exchange + VWN correlation).
    Lda,
    /// Generalized gradient approximation (PBE).
    Pbe,
    /// Hybrid B3LYP functional.
    B3lyp,
    /// Meta-GGA TPSS.
    Tpss,
}
/// Configuration Interaction Singles (CIS) excited state.
#[derive(Debug, Clone)]
pub struct CisState {
    /// Excitation energy \[Hartree\].
    pub excitation_energy: f64,
    /// Oscillator strength (dimensionless).
    pub oscillator_strength: f64,
    /// Dominant excitation amplitude (i → a).
    pub dominant_occ: usize,
    /// Virtual orbital index in dominant excitation.
    pub dominant_virt: usize,
    /// CIS amplitudes vector.
    pub amplitudes: Vec<f64>,
}
/// Result of a Hartree-Fock calculation.
#[derive(Debug, Clone)]
pub struct HartreeFockResult {
    /// Total electronic energy \[Hartree\].
    pub energy: f64,
    /// Nuclear repulsion energy \[Hartree\].
    pub nuclear_repulsion: f64,
    /// Orbital energies (eigenvalues) \[Hartree\].
    pub orbital_energies: Vec<f64>,
    /// MO coefficient matrix (columns = MOs, rows = AOs).
    pub mo_coefficients: Vec<Vec<f64>>,
    /// Density matrix elements.
    pub density_matrix: Vec<Vec<f64>>,
    /// Number of SCF iterations performed.
    pub iterations: usize,
    /// SCF convergence status.
    pub status: ScfStatus,
}
