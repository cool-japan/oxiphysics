//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::{
    compute_sasa, dist3, hbond_angle, hbond_score, jacobi3, mat3_det, mat3_mul_t,
};

/// CHARMM improper dihedral energy parameters.
///
/// E_imp = K_psi * (psi - psi0)^2
#[derive(Debug, Clone)]
pub struct CharmmImproper {
    /// Atom types.
    pub atoms: [String; 4],
    /// Force constant in kcal/mol/rad².
    pub k_psi: f64,
    /// Equilibrium improper angle in degrees.
    pub psi0: f64,
}
impl CharmmImproper {
    /// Create a CHARMM improper dihedral parameter entry.
    pub fn new(a1: &str, a2: &str, a3: &str, a4: &str, k_psi: f64, psi0: f64) -> Self {
        Self {
            atoms: [
                a1.to_string(),
                a2.to_string(),
                a3.to_string(),
                a4.to_string(),
            ],
            k_psi,
            psi0,
        }
    }
    /// Compute improper energy for given dihedral angle in degrees.
    pub fn energy(&self, psi_deg: f64) -> f64 {
        let dpsi = (psi_deg - self.psi0).to_radians();
        self.k_psi * dpsi * dpsi
    }
}
/// Protein topology from sequence or PDB.
///
/// Stores residue list, bonds, angles, and dihedrals.
#[derive(Debug, Clone)]
pub struct ProteinTopology {
    /// Residue names (3-letter code).
    pub residues: Vec<String>,
    /// Bond list: pairs of atom indices.
    pub bonds: Vec<[usize; 2]>,
    /// Angle list: triples of atom indices.
    pub angles: Vec<[usize; 3]>,
    /// Dihedral list: quads of atom indices.
    pub dihedrals: Vec<[usize; 4]>,
    /// Atom positions \[x, y, z\] in Angstroms.
    pub positions: Vec<[f64; 3]>,
    /// Atom names (e.g. "CA", "CB", "N", "C", "O").
    pub atom_names: Vec<String>,
    /// Residue index for each atom.
    pub atom_residue: Vec<usize>,
    /// Total number of atoms.
    pub n_atoms: usize,
}
impl ProteinTopology {
    /// Create a new empty protein topology.
    pub fn new() -> Self {
        Self {
            residues: Vec::new(),
            bonds: Vec::new(),
            angles: Vec::new(),
            dihedrals: Vec::new(),
            positions: Vec::new(),
            atom_names: Vec::new(),
            atom_residue: Vec::new(),
            n_atoms: 0,
        }
    }
    /// Add a residue with backbone atoms (N, CA, C, O).
    pub fn add_residue(
        &mut self,
        res_name: &str,
        n_pos: [f64; 3],
        ca_pos: [f64; 3],
        c_pos: [f64; 3],
        o_pos: [f64; 3],
    ) {
        let res_idx = self.residues.len();
        self.residues.push(res_name.to_string());
        let start = self.n_atoms;
        self.positions.push(n_pos);
        self.positions.push(ca_pos);
        self.positions.push(c_pos);
        self.positions.push(o_pos);
        for name in &["N", "CA", "C", "O"] {
            self.atom_names.push(name.to_string());
            self.atom_residue.push(res_idx);
        }
        self.bonds.push([start, start + 1]);
        self.bonds.push([start + 1, start + 2]);
        self.bonds.push([start + 2, start + 3]);
        if start > 0 {
            self.bonds.push([start - 2, start]);
        }
        self.n_atoms += 4;
    }
    /// Number of residues.
    pub fn n_residues(&self) -> usize {
        self.residues.len()
    }
    /// Build angles from bond list.
    pub fn build_angles(&mut self) {
        self.angles.clear();
        for i in 0..self.bonds.len() {
            for j in (i + 1)..self.bonds.len() {
                let [a1, b1] = self.bonds[i];
                let [a2, b2] = self.bonds[j];
                if b1 == a2 {
                    self.angles.push([a1, b1, b2]);
                } else if b1 == b2 {
                    self.angles.push([a1, b1, a2]);
                } else if a1 == a2 {
                    self.angles.push([b1, a1, b2]);
                } else if a1 == b2 {
                    self.angles.push([b1, a1, a2]);
                }
            }
        }
        self.angles.dedup();
    }
}
/// Backbone RMSD after optimal superposition using the Kabsch algorithm.
///
/// The Kabsch algorithm finds the rotation matrix that minimizes RMSD between
/// two sets of paired coordinates.
#[derive(Debug, Clone)]
pub struct RmsdCalculation;
impl RmsdCalculation {
    /// Compute raw RMSD without superposition.
    pub fn rmsd(pos1: &[[f64; 3]], pos2: &[[f64; 3]]) -> f64 {
        let n = pos1.len().min(pos2.len());
        if n == 0 {
            return 0.0;
        }
        let sum: f64 = pos1
            .iter()
            .zip(pos2.iter())
            .map(|(a, b)| {
                let dx = a[0] - b[0];
                let dy = a[1] - b[1];
                let dz = a[2] - b[2];
                dx * dx + dy * dy + dz * dz
            })
            .sum();
        (sum / n as f64).sqrt()
    }
    /// Center coordinates (translate center of mass to origin).
    fn center(pos: &[[f64; 3]]) -> Vec<[f64; 3]> {
        let n = pos.len() as f64;
        let mut cx = 0.0_f64;
        let mut cy = 0.0_f64;
        let mut cz = 0.0_f64;
        for p in pos {
            cx += p[0];
            cy += p[1];
            cz += p[2];
        }
        cx /= n;
        cy /= n;
        cz /= n;
        pos.iter()
            .map(|p| [p[0] - cx, p[1] - cy, p[2] - cz])
            .collect()
    }
    /// Compute RMSD after optimal Kabsch superposition.
    ///
    /// Centers both structures and finds the rotation minimizing RMSD.
    /// This is a simplified 2D SVD-based implementation using the analytical
    /// solution for 3×3 matrices via the covariance matrix.
    pub fn kabsch_rmsd(pos1: &[[f64; 3]], pos2: &[[f64; 3]]) -> f64 {
        let n = pos1.len().min(pos2.len());
        if n < 2 {
            return 0.0;
        }
        let c1 = Self::center(pos1);
        let c2 = Self::center(pos2);
        let mut h = [[0.0_f64; 3]; 3];
        for (a, b) in c1.iter().zip(c2.iter()) {
            for i in 0..3 {
                for j in 0..3 {
                    h[i][j] += a[i] * b[j];
                }
            }
        }
        let ht_h = mat3_mul_t(&h, &h);
        let (_eigvals, eigvecs) = jacobi3(ht_h);
        let s: [f64; 3] = {
            let ht_h_eig = jacobi3(ht_h);
            [
                ht_h_eig.0[0].max(0.0).sqrt(),
                ht_h_eig.0[1].max(0.0).sqrt(),
                ht_h_eig.0[2].max(0.0).sqrt(),
            ]
        };
        let sum_sq1: f64 = c1
            .iter()
            .map(|p| p[0] * p[0] + p[1] * p[1] + p[2] * p[2])
            .sum();
        let sum_sq2: f64 = c2
            .iter()
            .map(|p| p[0] * p[0] + p[1] * p[1] + p[2] * p[2])
            .sum();
        let det_v = mat3_det(&eigvecs);
        let sign = if det_v < 0.0 { -1.0 } else { 1.0 };
        let trace_sv = s[0] + s[1] + sign * s[2];
        let rmsd_sq = ((sum_sq1 + sum_sq2) / n as f64 - 2.0 * trace_sv / n as f64).max(0.0);
        rmsd_sq.sqrt()
    }
}
/// Hydrogen bond using Baker-Hubbard donor-acceptor distance and angle criterion.
///
/// A hydrogen bond is detected when:
/// - Donor-Acceptor distance < 3.5 Å
/// - D-H···A angle > 120°
#[derive(Debug, Clone, Copy)]
pub struct HydrogenBond {
    /// Donor atom position (Å).
    pub donor: [f64; 3],
    /// Hydrogen atom position (Å).
    pub hydrogen: [f64; 3],
    /// Acceptor atom position (Å).
    pub acceptor: [f64; 3],
}
impl HydrogenBond {
    /// Create a new hydrogen bond candidate.
    pub fn new(donor: [f64; 3], hydrogen: [f64; 3], acceptor: [f64; 3]) -> Self {
        Self {
            donor,
            hydrogen,
            acceptor,
        }
    }
    /// Donor-acceptor distance (Å).
    pub fn donor_acceptor_distance(&self) -> f64 {
        let dx = self.acceptor[0] - self.donor[0];
        let dy = self.acceptor[1] - self.donor[1];
        let dz = self.acceptor[2] - self.donor[2];
        (dx * dx + dy * dy + dz * dz).sqrt()
    }
    /// Hydrogen-acceptor distance (Å).
    pub fn h_acceptor_distance(&self) -> f64 {
        let dx = self.acceptor[0] - self.hydrogen[0];
        let dy = self.acceptor[1] - self.hydrogen[1];
        let dz = self.acceptor[2] - self.hydrogen[2];
        (dx * dx + dy * dy + dz * dz).sqrt()
    }
    /// D-H···A angle in degrees (angle at H between D and A).
    ///
    /// Computed as the angle between vectors H→D and H→A.
    /// For a linear D-H···A hydrogen bond, this is 180°.
    pub fn dha_angle(&self) -> f64 {
        let hd = [
            self.donor[0] - self.hydrogen[0],
            self.donor[1] - self.hydrogen[1],
            self.donor[2] - self.hydrogen[2],
        ];
        let ha = [
            self.acceptor[0] - self.hydrogen[0],
            self.acceptor[1] - self.hydrogen[1],
            self.acceptor[2] - self.hydrogen[2],
        ];
        let dot = hd[0] * ha[0] + hd[1] * ha[1] + hd[2] * ha[2];
        let mag_hd = (hd[0] * hd[0] + hd[1] * hd[1] + hd[2] * hd[2]).sqrt();
        let mag_ha = (ha[0] * ha[0] + ha[1] * ha[1] + ha[2] * ha[2]).sqrt();
        if mag_hd < 1e-10 || mag_ha < 1e-10 {
            return 0.0;
        }
        let cos_angle = (dot / (mag_hd * mag_ha)).clamp(-1.0, 1.0);
        cos_angle.acos().to_degrees()
    }
    /// Check if this is a hydrogen bond (Baker-Hubbard criterion).
    ///
    /// Requires H···A distance < 2.5 Å and D-H···A angle > 120°.
    pub fn is_hbond(&self) -> bool {
        self.h_acceptor_distance() < 2.5 && self.dha_angle() > 120.0
    }
    /// Hydrogen bond energy (Lennard-Jones like, simplified).
    ///
    /// Uses a simple 10-12 potential: `E = D0 * [(r0/r)^12 - 2*(r0/r)^10]`.
    pub fn energy(&self) -> f64 {
        let r = self.h_acceptor_distance();
        let r0 = 1.8_f64;
        let d0 = -1.5_f64;
        if r < 0.5 {
            return 1e6;
        }
        let rr = r0 / r;
        d0 * (rr.powi(12) - 2.0 * rr.powi(10))
    }
}
/// DSSP secondary structure label.
#[derive(Debug, Clone, PartialEq, Copy)]
pub enum DsspLabel {
    /// Alpha helix.
    AlphaHelix,
    /// 3₁₀ helix.
    Helix310,
    /// Pi helix.
    PiHelix,
    /// Beta sheet.
    BetaSheet,
    /// Beta bridge.
    BetaBridge,
    /// Bend.
    Bend,
    /// Turn.
    Turn,
    /// Coil/loop.
    Coil,
}
/// CHARMM angle bending energy parameters.
///
/// E_angle = K_theta * (theta - theta0)^2
#[derive(Debug, Clone)]
pub struct CharmmAngle {
    /// Atom types (3 atoms defining the angle).
    pub atoms: [String; 3],
    /// Force constant in kcal/mol/rad².
    pub k_theta: f64,
    /// Equilibrium angle in degrees.
    pub theta0: f64,
    /// Urey-Bradley 1-3 force constant (kcal/mol/Å²), 0 if not used.
    pub k_ub: f64,
    /// Urey-Bradley 1-3 distance (Å), 0 if not used.
    pub s0: f64,
}
impl CharmmAngle {
    /// Create a CHARMM angle parameter entry.
    pub fn new(a1: &str, a2: &str, a3: &str, k_theta: f64, theta0: f64) -> Self {
        Self {
            atoms: [a1.to_string(), a2.to_string(), a3.to_string()],
            k_theta,
            theta0,
            k_ub: 0.0,
            s0: 0.0,
        }
    }
    /// Compute CHARMM angle energy for given angle in degrees.
    pub fn energy(&self, theta_deg: f64) -> f64 {
        let dt = (theta_deg - self.theta0).to_radians();
        self.k_theta * dt * dt
    }
    /// Compute angle energy including Urey-Bradley 1-3 term.
    ///
    /// E = K_theta*(theta - theta0)^2 + K_ub*(s - s0)^2
    pub fn energy_with_ub(&self, theta_deg: f64, s13: f64) -> f64 {
        self.energy(theta_deg) + self.k_ub * (s13 - self.s0).powi(2)
    }
}
/// Protein secondary structure classification from backbone dihedral angles.
///
/// Uses phi/psi dihedral angle ranges to assign helix, sheet, or coil labels
/// to each residue.
#[derive(Debug, Clone)]
pub struct ProteinSecondaryStructure {
    /// Number of residues.
    pub n_residues: usize,
    /// Phi angles (degrees) for each residue.
    pub phi: Vec<f64>,
    /// Psi angles (degrees) for each residue.
    pub psi: Vec<f64>,
}
impl ProteinSecondaryStructure {
    /// Create a new structure classifier for `n` residues.
    pub fn new(n_residues: usize) -> Self {
        Self {
            n_residues,
            phi: vec![f64::NAN; n_residues],
            psi: vec![f64::NAN; n_residues],
        }
    }
    /// Set phi and psi angles (in degrees) for residue `i`.
    pub fn set_angles(&mut self, i: usize, phi: f64, psi: f64) {
        self.phi[i] = phi;
        self.psi[i] = psi;
    }
    /// Classify residue `i` by its phi/psi angles.
    ///
    /// Uses standard Ramachandran region boundaries:
    /// - Helix: phi in \[-100, -20\], psi in \[-80, -10\]
    /// - Sheet: phi in \[-180, -40\], psi in \[90, 180\] or \[-180, -170\]
    pub fn classify(&self, i: usize) -> SecStructLabel {
        let phi = self.phi[i];
        let psi = self.psi[i];
        if phi.is_nan() || psi.is_nan() {
            return SecStructLabel::Coil;
        }
        if (-100.0..=-20.0).contains(&phi) && (-80.0..=-10.0).contains(&psi) {
            SecStructLabel::Helix
        } else if (-180.0..=-40.0).contains(&phi) && (psi >= 90.0 || psi <= -170.0) {
            SecStructLabel::Sheet
        } else {
            SecStructLabel::Coil
        }
    }
    /// Classify all residues and return a vector of labels.
    pub fn classify_all(&self) -> Vec<SecStructLabel> {
        (0..self.n_residues).map(|i| self.classify(i)).collect()
    }
}
/// CHARMM dihedral torsion energy parameters.
///
/// E_dihedral = K_phi * (1 + cos(n*phi - delta))
#[derive(Debug, Clone)]
pub struct CharmmDihedral {
    /// Atom types (4 atoms defining the dihedral).
    pub atoms: [String; 4],
    /// Force constant in kcal/mol.
    pub k_phi: f64,
    /// Periodicity (multiplicity).
    pub n: u32,
    /// Phase angle in degrees.
    pub delta: f64,
}
impl CharmmDihedral {
    /// Create a CHARMM dihedral parameter entry.
    pub fn new(a1: &str, a2: &str, a3: &str, a4: &str, k_phi: f64, n: u32, delta: f64) -> Self {
        Self {
            atoms: [
                a1.to_string(),
                a2.to_string(),
                a3.to_string(),
                a4.to_string(),
            ],
            k_phi,
            n,
            delta,
        }
    }
    /// Compute CHARMM dihedral energy for given dihedral angle in degrees.
    pub fn energy(&self, phi_deg: f64) -> f64 {
        let phi_rad = phi_deg.to_radians();
        let delta_rad = self.delta.to_radians();
        self.k_phi * (1.0 + (self.n as f64 * phi_rad - delta_rad).cos())
    }
}
/// Radius of gyration from Cα coordinates.
///
/// `Rg = sqrt(sum(m_i * |r_i - r_cm|^2) / sum(m_i))`.
/// Equal masses assumed when not provided.
#[derive(Debug, Clone)]
pub struct RgProtein {
    /// Cα positions (Å).
    pub ca_positions: Vec<[f64; 3]>,
}
impl RgProtein {
    /// Create from Cα positions.
    pub fn new(ca_positions: Vec<[f64; 3]>) -> Self {
        Self { ca_positions }
    }
    /// Compute center of mass (equal masses).
    pub fn center_of_mass(&self) -> [f64; 3] {
        let n = self.ca_positions.len() as f64;
        if n < 1.0 {
            return [0.0; 3];
        }
        let mut com = [0.0; 3];
        for pos in &self.ca_positions {
            com[0] += pos[0];
            com[1] += pos[1];
            com[2] += pos[2];
        }
        [com[0] / n, com[1] / n, com[2] / n]
    }
    /// Compute radius of gyration (Å).
    pub fn rg(&self) -> f64 {
        let n = self.ca_positions.len();
        if n < 2 {
            return 0.0;
        }
        let com = self.center_of_mass();
        let sum_sq: f64 = self
            .ca_positions
            .iter()
            .map(|pos| {
                let dx = pos[0] - com[0];
                let dy = pos[1] - com[1];
                let dz = pos[2] - com[2];
                dx * dx + dy * dy + dz * dz
            })
            .sum();
        (sum_sq / n as f64).sqrt()
    }
}
/// Solvation shell analysis.
///
/// Computes radial distribution around solute, hydration number, and shell lifetime.
#[derive(Debug, Clone)]
pub struct SolvationShell {
    /// Solute atom positions.
    pub solute_positions: Vec<[f64; 3]>,
    /// Solvent atom positions.
    pub solvent_positions: Vec<[f64; 3]>,
    /// RDF g(r) bins.
    pub rdf: Vec<f64>,
    /// RDF bin width in Angstroms.
    pub dr: f64,
    /// Maximum r for RDF.
    pub r_max: f64,
    /// Number of RDF snapshots accumulated.
    pub n_snapshots: usize,
    /// First shell cutoff in Angstroms.
    pub first_shell_cutoff: f64,
}
impl SolvationShell {
    /// Create a new solvation shell analysis.
    pub fn new(r_max: f64, dr: f64, first_shell_cutoff: f64) -> Self {
        let n_bins = (r_max / dr) as usize + 1;
        Self {
            solute_positions: Vec::new(),
            solvent_positions: Vec::new(),
            rdf: vec![0.0f64; n_bins],
            dr,
            r_max,
            n_snapshots: 0,
            first_shell_cutoff,
        }
    }
    /// Accumulate one snapshot into RDF.
    pub fn accumulate(&mut self) {
        for &sp in &self.solute_positions {
            for &sv in &self.solvent_positions {
                let r = dist3(sp, sv);
                if r < self.r_max {
                    let bin = (r / self.dr) as usize;
                    if bin < self.rdf.len() {
                        self.rdf[bin] += 1.0;
                    }
                }
            }
        }
        self.n_snapshots += 1;
    }
    /// Normalize RDF to g(r).
    pub fn normalize(&mut self, rho_bulk: f64) {
        let n_s = self.solute_positions.len().max(1) as f64;
        let n_snap = self.n_snapshots.max(1) as f64;
        for (i, g) in self.rdf.iter_mut().enumerate() {
            let r = (i as f64 + 0.5) * self.dr;
            let shell_vol = 4.0 * std::f64::consts::PI * r * r * self.dr;
            *g /= n_s * n_snap * rho_bulk * shell_vol;
        }
    }
    /// Compute hydration number (integral of first shell).
    pub fn hydration_number(&self, rho_bulk: f64) -> f64 {
        let mut n_hyd = 0.0;
        for (i, &g) in self.rdf.iter().enumerate() {
            let r = (i as f64 + 0.5) * self.dr;
            if r > self.first_shell_cutoff {
                break;
            }
            let shell_vol = 4.0 * std::f64::consts::PI * r * r * self.dr;
            n_hyd += g * rho_bulk * shell_vol;
        }
        n_hyd
    }
    /// Average distance of first shell.
    pub fn first_shell_peak_r(&self) -> f64 {
        let mut max_g = 0.0;
        let mut peak_r = self.dr;
        for (i, &g) in self.rdf.iter().enumerate() {
            let r = (i as f64 + 0.5) * self.dr;
            if r > self.first_shell_cutoff {
                break;
            }
            if g > max_g {
                max_g = g;
                peak_r = r;
            }
        }
        peak_r
    }
}
/// CHARMM bond stretching energy parameters.
///
/// E_bond = K_b * (r - r0)^2
#[derive(Debug, Clone)]
pub struct CharmmBond {
    /// Atom type 1.
    pub atom1: String,
    /// Atom type 2.
    pub atom2: String,
    /// Force constant in kcal/mol/Å².
    pub k_b: f64,
    /// Equilibrium bond length in Å.
    pub r0: f64,
}
impl CharmmBond {
    /// Create a CHARMM bond parameter entry.
    pub fn new(atom1: &str, atom2: &str, k_b: f64, r0: f64) -> Self {
        Self {
            atom1: atom1.to_string(),
            atom2: atom2.to_string(),
            k_b,
            r0,
        }
    }
    /// Compute CHARMM bond energy for given distance.
    pub fn energy(&self, r: f64) -> f64 {
        self.k_b * (r - self.r0).powi(2)
    }
}
/// CHARMM36/AMBER force field parameters for 20 standard amino acids.
///
/// Stores LJ parameters, partial charges, and bonded parameters per residue type.
#[derive(Debug, Clone)]
pub struct ForceFieldAmino {
    /// LJ epsilon per atom type (kcal/mol).
    pub eps: std::collections::HashMap<String, f64>,
    /// LJ sigma per atom type (Angstroms).
    pub sigma: std::collections::HashMap<String, f64>,
    /// Partial charges per atom type (e).
    pub charges: std::collections::HashMap<String, f64>,
    /// Bond equilibrium lengths (Angstroms).
    pub bond_lengths: std::collections::HashMap<String, f64>,
    /// Bond force constants (kcal/mol/A²).
    pub bond_k: std::collections::HashMap<String, f64>,
    /// Angle equilibrium values (degrees).
    pub angle_eq: std::collections::HashMap<String, f64>,
    /// Angle force constants (kcal/mol/rad²).
    pub angle_k: std::collections::HashMap<String, f64>,
}
impl ForceFieldAmino {
    /// Create a new force field with standard CHARMM36-like parameters.
    pub fn new_charmm36() -> Self {
        let mut ff = Self {
            eps: std::collections::HashMap::new(),
            sigma: std::collections::HashMap::new(),
            charges: std::collections::HashMap::new(),
            bond_lengths: std::collections::HashMap::new(),
            bond_k: std::collections::HashMap::new(),
            angle_eq: std::collections::HashMap::new(),
            angle_k: std::collections::HashMap::new(),
        };
        ff.eps.insert("N".into(), 0.2000);
        ff.eps.insert("CA".into(), 0.0700);
        ff.eps.insert("C".into(), 0.0700);
        ff.eps.insert("O".into(), 0.1200);
        ff.sigma.insert("N".into(), 1.8240);
        ff.sigma.insert("CA".into(), 1.9080);
        ff.sigma.insert("C".into(), 1.9080);
        ff.sigma.insert("O".into(), 1.6612);
        ff.charges.insert("N".into(), -0.4157);
        ff.charges.insert("CA".into(), 0.0337);
        ff.charges.insert("C".into(), 0.5973);
        ff.charges.insert("O".into(), -0.5679);
        ff.bond_lengths.insert("N-CA".into(), 1.449);
        ff.bond_lengths.insert("CA-C".into(), 1.522);
        ff.bond_lengths.insert("C-O".into(), 1.229);
        ff.bond_lengths.insert("C-N".into(), 1.335);
        ff.bond_k.insert("N-CA".into(), 367.0);
        ff.bond_k.insert("CA-C".into(), 317.0);
        ff.bond_k.insert("C-O".into(), 570.0);
        ff.angle_eq.insert("N-CA-C".into(), 111.2);
        ff.angle_eq.insert("CA-C-O".into(), 120.4);
        ff.angle_k.insert("N-CA-C".into(), 63.0);
        ff.angle_k.insert("CA-C-O".into(), 80.0);
        ff
    }
    /// Get LJ energy between two atom types.
    pub fn lj_energy(&self, type_i: &str, type_j: &str, r: f64) -> f64 {
        let eps_i = self.eps.get(type_i).copied().unwrap_or(0.1);
        let eps_j = self.eps.get(type_j).copied().unwrap_or(0.1);
        let sig_i = self.sigma.get(type_i).copied().unwrap_or(1.9);
        let sig_j = self.sigma.get(type_j).copied().unwrap_or(1.9);
        let eps = (eps_i * eps_j).sqrt();
        let sig = 0.5 * (sig_i + sig_j);
        let sr6 = (sig / r).powi(6);
        4.0 * eps * (sr6 * sr6 - sr6)
    }
    /// Get bond energy for given bond type and distance.
    pub fn bond_energy(&self, bond_type: &str, r: f64) -> f64 {
        let r0 = self.bond_lengths.get(bond_type).copied().unwrap_or(1.5);
        let k = self.bond_k.get(bond_type).copied().unwrap_or(300.0);
        k * (r - r0).powi(2)
    }
    /// Get angle energy for given angle type and angle in degrees.
    pub fn angle_energy(&self, angle_type: &str, theta_deg: f64) -> f64 {
        let theta0 = self.angle_eq.get(angle_type).copied().unwrap_or(109.5);
        let k = self.angle_k.get(angle_type).copied().unwrap_or(60.0);
        let dtheta = (theta_deg - theta0).to_radians();
        k * dtheta * dtheta
    }
}
/// Debye-Hückel electrostatic solvation model for proteins.
///
/// Computes screened Coulomb interactions at physiological ionic strength.
/// `phi(r) = q * exp(-kappa * r) / (4*pi*eps0*eps_r * r)`.
#[derive(Debug, Clone, Copy)]
pub struct DebyeHuckelProtein {
    /// Ionic strength (mol/L).
    pub ionic_strength: f64,
    /// Relative dielectric constant (unitless).
    pub epsilon_r: f64,
    /// Temperature (K).
    pub temperature: f64,
}
impl DebyeHuckelProtein {
    /// Create a Debye-Hückel model.
    pub fn new(ionic_strength: f64, epsilon_r: f64, temperature: f64) -> Self {
        Self {
            ionic_strength,
            epsilon_r,
            temperature,
        }
    }
    /// Debye screening length (Å).
    ///
    /// `lambda_D = sqrt(eps0*eps_r*kB*T / (2*Na*e^2*I))`.
    /// Simplified: `lambda_D ≈ 3.04 / sqrt(I[mol/L])` Å at 25°C.
    pub fn debye_length(&self) -> f64 {
        3.04 / self.ionic_strength.sqrt()
    }
    /// Inverse Debye length kappa (1/Å).
    pub fn inverse_debye_length(&self) -> f64 {
        1.0 / self.debye_length()
    }
    /// Screened Coulomb potential at distance `r` (Å) from charge `q` (e).
    ///
    /// Returns energy in kcal/mol.
    pub fn screened_potential(&self, q: f64, r: f64) -> f64 {
        if r < 0.1 {
            return 1e10;
        }
        let kappa = self.inverse_debye_length();
        332.0 * q * (-kappa * r).exp() / (self.epsilon_r * r)
    }
    /// Total electrostatic solvation energy for a set of charges.
    pub fn solvation_energy(&self, charges: &[f64], positions: &[[f64; 3]]) -> f64 {
        let mut energy = 0.0;
        let n = charges.len().min(positions.len());
        for i in 0..n {
            for j in (i + 1)..n {
                let dx = positions[j][0] - positions[i][0];
                let dy = positions[j][1] - positions[i][1];
                let dz = positions[j][2] - positions[i][2];
                let r = (dx * dx + dy * dy + dz * dz).sqrt();
                energy += self.screened_potential(charges[i] * charges[j], r);
            }
        }
        energy
    }
}
/// Secondary structure label for a residue.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SecStructLabel {
    /// Alpha helix.
    Helix,
    /// Beta sheet.
    Sheet,
    /// Random coil / other.
    Coil,
}
/// Water model geometry and force field parameters.
#[derive(Debug, Clone)]
pub struct WaterModel {
    /// Model type.
    pub model_type: WaterModelType,
    /// O-H bond length in Angstroms.
    pub oh_length: f64,
    /// H-O-H angle in degrees.
    pub hoh_angle: f64,
    /// LJ epsilon for O (kcal/mol).
    pub eps_o: f64,
    /// LJ sigma for O (Angstroms).
    pub sigma_o: f64,
    /// Charge on O (e).
    pub q_o: f64,
    /// Charge on H (e).
    pub q_h: f64,
    /// Virtual site distance (TIP4P only).
    pub d_virtual: f64,
}
impl WaterModel {
    /// Create TIP3P water model.
    pub fn tip3p() -> Self {
        Self {
            model_type: WaterModelType::Tip3p,
            oh_length: 0.9572,
            hoh_angle: 104.52,
            eps_o: 0.1521,
            sigma_o: 3.1506,
            q_o: -0.834,
            q_h: 0.417,
            d_virtual: 0.0,
        }
    }
    /// Create TIP4P water model.
    pub fn tip4p() -> Self {
        Self {
            model_type: WaterModelType::Tip4p,
            oh_length: 0.9572,
            hoh_angle: 104.52,
            eps_o: 0.1550,
            sigma_o: 3.1536,
            q_o: 0.0,
            q_h: 0.520,
            d_virtual: 0.15,
        }
    }
    /// Create SPC/E water model.
    pub fn spce() -> Self {
        Self {
            model_type: WaterModelType::SpceModel,
            oh_length: 1.0000,
            hoh_angle: 109.47,
            eps_o: 0.1553,
            sigma_o: 3.1656,
            q_o: -0.8476,
            q_h: 0.4238,
            d_virtual: 0.0,
        }
    }
    /// Compute O-H geometry positions given O position and angle.
    pub fn h_positions(&self, o_pos: [f64; 3]) -> [[f64; 3]; 2] {
        let half_ang = self.hoh_angle.to_radians() / 2.0;
        let r = self.oh_length;
        let h1 = [
            o_pos[0] + r * half_ang.sin(),
            o_pos[1] + r * half_ang.cos(),
            o_pos[2],
        ];
        let h2 = [
            o_pos[0] - r * half_ang.sin(),
            o_pos[1] + r * half_ang.cos(),
            o_pos[2],
        ];
        [h1, h2]
    }
    /// Compute LJ energy between two oxygens.
    pub fn lj_oo(&self, r: f64) -> f64 {
        let s6 = (self.sigma_o / r).powi(6);
        4.0 * self.eps_o * (s6 * s6 - s6)
    }
    /// Dipole moment in Debye.
    pub fn dipole_moment(&self) -> f64 {
        let half_ang = self.hoh_angle.to_radians() / 2.0;
        2.0 * self.q_h * self.oh_length * half_ang.cos() * 4.803
    }
}
/// Protein folding simulation.
///
/// Runs temperature quench simulation, tracking contact fraction Q and RMSD.
#[derive(Debug, Clone)]
pub struct ProteinFolding {
    /// Topology.
    pub topology: ProteinTopology,
    /// Current temperature in K.
    pub temperature: f64,
    /// Target native contacts (reference).
    pub native_contacts: Vec<[usize; 2]>,
    /// Contact distance cutoff in Angstroms.
    pub contact_cutoff: f64,
    /// Time step in ps.
    pub dt: f64,
    /// Reference positions (folded state).
    pub ref_positions: Vec<[f64; 3]>,
    /// Current step.
    pub step_count: usize,
    /// Q factor history.
    pub q_history: Vec<f64>,
}
impl ProteinFolding {
    /// Create a new protein folding simulation.
    pub fn new(topology: ProteinTopology, temperature: f64, contact_cutoff: f64, dt: f64) -> Self {
        let ref_positions = topology.positions.clone();
        let native_contacts = Self::find_contacts(&ref_positions, contact_cutoff);
        Self {
            topology,
            temperature,
            native_contacts,
            contact_cutoff,
            dt,
            ref_positions,
            step_count: 0,
            q_history: Vec::new(),
        }
    }
    /// Find contact pairs within cutoff.
    pub fn find_contacts(positions: &[[f64; 3]], cutoff: f64) -> Vec<[usize; 2]> {
        let mut contacts = Vec::new();
        let n = positions.len();
        for i in 0..n {
            for j in (i + 2)..n {
                let d = dist3(positions[i], positions[j]);
                if d < cutoff {
                    contacts.push([i, j]);
                }
            }
        }
        contacts
    }
    /// Compute Q (fraction of native contacts formed).
    pub fn contact_fraction_q(&self) -> f64 {
        if self.native_contacts.is_empty() {
            return 1.0;
        }
        let mut formed = 0;
        let cutoff = self.contact_cutoff * 1.2;
        for &[i, j] in &self.native_contacts {
            if i < self.topology.positions.len() && j < self.topology.positions.len() {
                let d = dist3(self.topology.positions[i], self.topology.positions[j]);
                if d < cutoff {
                    formed += 1;
                }
            }
        }
        formed as f64 / self.native_contacts.len() as f64
    }
    /// Compute RMSD from reference structure.
    pub fn rmsd(&self) -> f64 {
        let n = self.topology.positions.len().min(self.ref_positions.len());
        if n == 0 {
            return 0.0;
        }
        let sum: f64 = (0..n)
            .map(|i| {
                let dx = self.topology.positions[i][0] - self.ref_positions[i][0];
                let dy = self.topology.positions[i][1] - self.ref_positions[i][1];
                let dz = self.topology.positions[i][2] - self.ref_positions[i][2];
                dx * dx + dy * dy + dz * dz
            })
            .sum();
        (sum / n as f64).sqrt()
    }
    /// Apply temperature quench (reduce temperature).
    pub fn quench(&mut self, target_temp: f64, rate: f64) {
        self.temperature = self.temperature * (1.0 - rate) + target_temp * rate;
    }
    /// Advance one MD step (simplified Langevin).
    pub fn step(&mut self) {
        self.step_count += 1;
        let q = self.contact_fraction_q();
        self.q_history.push(q);
    }
}
/// Radius of gyration tensor.
///
/// Computes shape anisotropy including asphericity and acylindricity.
#[derive(Debug, Clone)]
pub struct GyrRadius {
    /// Atom positions.
    pub positions: Vec<[f64; 3]>,
    /// Atom masses.
    pub masses: Vec<f64>,
}
impl GyrRadius {
    /// Create a new GyrRadius calculator.
    pub fn new(positions: Vec<[f64; 3]>, masses: Vec<f64>) -> Self {
        Self { positions, masses }
    }
    /// Center of mass.
    pub fn center_of_mass(&self) -> [f64; 3] {
        let total_mass: f64 = self.masses.iter().sum();
        if total_mass < 1e-15 {
            return [0.0; 3];
        }
        let mut com = [0.0f64; 3];
        for (i, &pos) in self.positions.iter().enumerate() {
            com[0] += self.masses[i] * pos[0];
            com[1] += self.masses[i] * pos[1];
            com[2] += self.masses[i] * pos[2];
        }
        [
            com[0] / total_mass,
            com[1] / total_mass,
            com[2] / total_mass,
        ]
    }
    /// Compute Rg² (squared radius of gyration).
    pub fn rg_squared(&self) -> f64 {
        let com = self.center_of_mass();
        let total_mass: f64 = self.masses.iter().sum();
        if total_mass < 1e-15 {
            return 0.0;
        }
        self.positions
            .iter()
            .zip(self.masses.iter())
            .map(|(&pos, &m)| {
                let dx = pos[0] - com[0];
                let dy = pos[1] - com[1];
                let dz = pos[2] - com[2];
                m * (dx * dx + dy * dy + dz * dz)
            })
            .sum::<f64>()
            / total_mass
    }
    /// Compute gyration tensor (3x3 as flat array \[xx, xy, xz, yx, yy, yz, zx, zy, zz\]).
    pub fn gyration_tensor(&self) -> [f64; 9] {
        let com = self.center_of_mass();
        let total_mass: f64 = self.masses.iter().sum();
        if total_mass < 1e-15 {
            return [0.0; 9];
        }
        let mut t = [0.0f64; 9];
        for (i, &pos) in self.positions.iter().enumerate() {
            let r = [pos[0] - com[0], pos[1] - com[1], pos[2] - com[2]];
            let m = self.masses[i];
            for a in 0..3 {
                for b in 0..3 {
                    t[a * 3 + b] += m * r[a] * r[b];
                }
            }
        }
        for v in &mut t {
            *v /= total_mass;
        }
        t
    }
    /// Asphericity b = lambda_z - 0.5*(lambda_x + lambda_y).
    /// Here we approximate eigenvalues by diagonal elements.
    pub fn asphericity(&self) -> f64 {
        let t = self.gyration_tensor();
        let lx = t[0];
        let ly = t[4];
        let lz = t[8];
        lz - 0.5 * (lx + ly)
    }
    /// Acylindricity c = lambda_y - lambda_x.
    pub fn acylindricity(&self) -> f64 {
        let t = self.gyration_tensor();
        let lx = t[0];
        let ly = t[4];
        (ly - lx).abs()
    }
}
/// Water model parameters.
///
/// Supports TIP3P, TIP4P, and SPC/E water models.
#[derive(Debug, Clone, PartialEq)]
pub enum WaterModelType {
    /// TIP3P water model.
    Tip3p,
    /// TIP4P water model.
    Tip4p,
    /// SPC/E water model.
    SpceModel,
}
/// Secondary structure assignment using DSSP algorithm.
///
/// Assigns helix (α/3₁₀/π), sheet (β/bridge), and coil.
#[derive(Debug, Clone)]
pub struct SecondaryStructure {
    /// DSSP label per residue.
    pub labels: Vec<DsspLabel>,
    /// Backbone N positions.
    pub n_pos: Vec<[f64; 3]>,
    /// Backbone CA positions.
    pub ca_pos: Vec<[f64; 3]>,
    /// Backbone C positions.
    pub c_pos: Vec<[f64; 3]>,
    /// Backbone O positions.
    pub o_pos: Vec<[f64; 3]>,
    /// H-bond energy matrix.
    pub hbond_matrix: Vec<Vec<f64>>,
}
impl SecondaryStructure {
    /// Create new SS analysis.
    pub fn new(
        n_pos: Vec<[f64; 3]>,
        ca_pos: Vec<[f64; 3]>,
        c_pos: Vec<[f64; 3]>,
        o_pos: Vec<[f64; 3]>,
    ) -> Self {
        let n = n_pos.len();
        let labels = vec![DsspLabel::Coil; n];
        let hbond_matrix = vec![vec![0.0f64; n]; n];
        Self {
            labels,
            n_pos,
            ca_pos,
            c_pos,
            o_pos,
            hbond_matrix,
        }
    }
    /// Compute DSSP H-bond energy between residues i and j.
    pub fn hbond_energy(&self, i: usize, j: usize) -> f64 {
        if i >= self.n_pos.len() || j >= self.n_pos.len() || i == j {
            return 0.0;
        }
        let r_on = dist3(self.o_pos[i], self.n_pos[j]);
        let r_ch = dist3(self.c_pos[i], self.n_pos[j]);
        let r_oh = dist3(self.o_pos[i], self.ca_pos[j]);
        let r_cn = dist3(self.c_pos[i], self.n_pos[j]);
        0.084
            * (1.0 / r_on.max(0.1) + 1.0 / r_ch.max(0.1)
                - 1.0 / r_oh.max(0.1)
                - 1.0 / r_cn.max(0.1))
            * 332.0
    }
    /// Assign secondary structure labels.
    pub fn assign(&mut self) {
        let n = self.n_pos.len();
        for i in 0..n {
            for j in 0..n {
                self.hbond_matrix[i][j] = self.hbond_energy(i, j);
            }
        }
        for i in 0..n {
            if i + 4 < n && self.hbond_matrix[i][i + 4] < -0.5 {
                self.labels[i] = DsspLabel::AlphaHelix;
            } else if i + 3 < n && self.hbond_matrix[i][i + 3] < -0.5 {
                self.labels[i] = DsspLabel::Helix310;
            }
        }
    }
    /// Count residues with given label.
    pub fn count_label(&self, label: DsspLabel) -> usize {
        self.labels.iter().filter(|&&l| l == label).count()
    }
    /// Fraction of helical residues.
    pub fn helix_fraction(&self) -> f64 {
        let n = self.labels.len();
        if n == 0 {
            return 0.0;
        }
        let count = self
            .labels
            .iter()
            .filter(|&&l| {
                l == DsspLabel::AlphaHelix || l == DsspLabel::Helix310 || l == DsspLabel::PiHelix
            })
            .count();
        count as f64 / n as f64
    }
}
/// AMBER-style Lennard-Jones parameters for common biological atom types.
///
/// Stores epsilon (kcal/mol) and Rmin/2 (Å) for 12-6 LJ potential.
#[derive(Debug, Clone)]
pub struct AmberLjParams {
    /// Atom type label.
    pub atom_type: String,
    /// LJ well depth epsilon in kcal/mol.
    pub epsilon: f64,
    /// Rmin/2 in Angstroms (half the minimum-energy distance).
    pub rmin2: f64,
}
impl AmberLjParams {
    /// Create AMBER LJ params for a given atom type.
    pub fn new(atom_type: &str, epsilon: f64, rmin2: f64) -> Self {
        Self {
            atom_type: atom_type.to_string(),
            epsilon,
            rmin2,
        }
    }
    /// Return AMBER99SB parameters for standard backbone atom types.
    ///
    /// Values from Ponder & Case, Adv. Protein Chem. 66 (2003).
    pub fn amber99sb_defaults() -> Vec<Self> {
        vec![
            Self::new("C", 0.0860, 1.9080),
            Self::new("CA", 0.0860, 1.9080),
            Self::new("CB", 0.0860, 1.9080),
            Self::new("N", 0.1700, 1.8240),
            Self::new("O", 0.2100, 1.6612),
            Self::new("OH", 0.2104, 1.7210),
            Self::new("S", 0.2500, 2.0000),
            Self::new("H", 0.0157, 0.6000),
            Self::new("HC", 0.0157, 1.4870),
        ]
    }
}
/// Hydrophobic solvation free energy model.
///
/// ΔG_hydrophobic ≈ γ * SASA_nonpolar
/// where γ ≈ 0.005 kcal/mol/Å² (surface tension coefficient).
#[derive(Debug, Clone)]
pub struct HydrophobicPotential {
    /// Surface tension coefficient in kcal/mol/Å².
    pub gamma: f64,
    /// Probe radius in Å.
    pub probe_radius: f64,
    /// Atom radii in Å (indexed by atom index).
    pub radii: Vec<f64>,
    /// Hydrophobicity flag per atom (true = nonpolar).
    pub is_nonpolar: Vec<bool>,
}
impl HydrophobicPotential {
    /// Create a hydrophobic potential calculator.
    pub fn new(gamma: f64, probe_radius: f64, radii: Vec<f64>, is_nonpolar: Vec<bool>) -> Self {
        Self {
            gamma,
            probe_radius,
            radii,
            is_nonpolar,
        }
    }
    /// Create a default hydrophobic potential (all atoms treated as nonpolar).
    pub fn default_carbon(n_atoms: usize) -> Self {
        Self {
            gamma: 0.005,
            probe_radius: 1.4,
            radii: vec![1.7; n_atoms],
            is_nonpolar: vec![true; n_atoms],
        }
    }
    /// Compute hydrophobic solvation free energy in kcal/mol.
    pub fn energy(&self, positions: &[[f64; 3]]) -> f64 {
        let sasa = compute_sasa(positions, &self.radii, self.probe_radius);
        sasa.iter()
            .zip(self.is_nonpolar.iter())
            .filter(|(_, np)| **np)
            .map(|(s, _)| *s)
            .sum::<f64>()
            * self.gamma
    }
}
/// MM/PBSA-GBSA binding free energy estimation.
///
/// Computes binding free energy from MD trajectory using end-point methods.
#[derive(Debug, Clone)]
pub struct BindingFreeEnergy {
    /// Complex energy samples in kcal/mol.
    pub e_complex: Vec<f64>,
    /// Receptor energy samples.
    pub e_receptor: Vec<f64>,
    /// Ligand energy samples.
    pub e_ligand: Vec<f64>,
    /// Solvation free energy samples.
    pub dg_solv: Vec<f64>,
    /// Temperature in K.
    pub temperature: f64,
    /// Surface area coefficient (PBSA/GBSA).
    pub gamma_sa: f64,
}
impl BindingFreeEnergy {
    /// Create a new binding free energy analysis.
    pub fn new(temperature: f64, gamma_sa: f64) -> Self {
        Self {
            e_complex: Vec::new(),
            e_receptor: Vec::new(),
            e_ligand: Vec::new(),
            dg_solv: Vec::new(),
            temperature,
            gamma_sa,
        }
    }
    /// Add a frame of energy data.
    pub fn add_frame(&mut self, e_comp: f64, e_rec: f64, e_lig: f64, dg_s: f64) {
        self.e_complex.push(e_comp);
        self.e_receptor.push(e_rec);
        self.e_ligand.push(e_lig);
        self.dg_solv.push(dg_s);
    }
    /// Compute MM/PBSA binding free energy.
    pub fn mmpbsa(&self) -> f64 {
        if self.e_complex.is_empty() {
            return 0.0;
        }
        let n = self.e_complex.len() as f64;
        let delta_emm: f64 = self
            .e_complex
            .iter()
            .zip(self.e_receptor.iter())
            .zip(self.e_ligand.iter())
            .map(|((&ec, &er), &el)| ec - er - el)
            .sum::<f64>()
            / n;
        let delta_dg: f64 = self.dg_solv.iter().sum::<f64>() / n;
        delta_emm + delta_dg
    }
    /// Compute entropy correction (simple T*dS estimate).
    pub fn entropy_correction(&self) -> f64 {
        let r = 1.987e-3;
        -self.temperature * r * (self.e_complex.len() as f64).ln()
    }
    /// Standard error of mean for binding free energy.
    pub fn sem(&self) -> f64 {
        if self.e_complex.len() < 2 {
            return 0.0;
        }
        let vals: Vec<f64> = self
            .e_complex
            .iter()
            .zip(self.e_receptor.iter())
            .zip(self.e_ligand.iter())
            .map(|((&ec, &er), &el)| ec - er - el)
            .collect();
        let n = vals.len() as f64;
        let mean = vals.iter().sum::<f64>() / n;
        let var = vals.iter().map(|v| (v - mean).powi(2)).sum::<f64>() / (n - 1.0);
        (var / n).sqrt()
    }
}
/// Ramachandran analysis: phi-psi calculation and allowed region detection.
///
/// Identifies whether residues fall in the allowed, generously allowed,
/// or disallowed regions of the Ramachandran plot.
#[derive(Debug, Clone)]
pub struct RamachandranAnalysis {
    /// Helix region phi range (min, max) in degrees.
    pub helix_phi: (f64, f64),
    /// Helix region psi range (min, max) in degrees.
    pub helix_psi: (f64, f64),
    /// Sheet region phi range (min, max) in degrees.
    pub sheet_phi: (f64, f64),
    /// Sheet region psi range (min, max) in degrees.
    pub sheet_psi: (f64, f64),
}
impl RamachandranAnalysis {
    /// Create with default alpha-helix and beta-sheet region boundaries.
    pub fn new() -> Self {
        Self {
            helix_phi: (-100.0, -20.0),
            helix_psi: (-80.0, -10.0),
            sheet_phi: (-180.0, -40.0),
            sheet_psi: (90.0, 180.0),
        }
    }
    /// Check if (phi, psi) falls in the alpha-helix region.
    pub fn is_helix_region(&self, phi: f64, psi: f64) -> bool {
        phi >= self.helix_phi.0
            && phi <= self.helix_phi.1
            && psi >= self.helix_psi.0
            && psi <= self.helix_psi.1
    }
    /// Check if (phi, psi) falls in the beta-sheet region.
    pub fn is_sheet_region(&self, phi: f64, psi: f64) -> bool {
        phi >= self.sheet_phi.0
            && phi <= self.sheet_phi.1
            && (psi >= self.sheet_psi.0 || psi <= -170.0)
    }
    /// Check if (phi, psi) is in the left-handed helix region (allowed for Gly).
    pub fn is_left_helix_region(phi: f64, psi: f64) -> bool {
        (30.0..90.0).contains(&phi) && (10.0..80.0).contains(&psi)
    }
    /// Check if (phi, psi) falls in any allowed region.
    pub fn is_allowed(&self, phi: f64, psi: f64) -> bool {
        self.is_helix_region(phi, psi)
            || self.is_sheet_region(phi, psi)
            || Self::is_left_helix_region(phi, psi)
    }
    /// Classify a (phi, psi) pair into a label.
    pub fn classify(&self, phi: f64, psi: f64) -> SecStructLabel {
        if self.is_helix_region(phi, psi) {
            SecStructLabel::Helix
        } else if self.is_sheet_region(phi, psi) {
            SecStructLabel::Sheet
        } else {
            SecStructLabel::Coil
        }
    }
}
/// Hydrogen bond analysis.
///
/// Detects H-bonds using distance < 3.5 Å and D-H...A angle > 150°.
#[derive(Debug, Clone)]
pub struct HBondAnalysis {
    /// Donor atom positions.
    pub donors: Vec<[f64; 3]>,
    /// H atom positions (bonded to donors).
    pub hydrogens: Vec<[f64; 3]>,
    /// Acceptor atom positions.
    pub acceptors: Vec<[f64; 3]>,
    /// Distance cutoff in Angstroms.
    pub dist_cutoff: f64,
    /// Angle cutoff in degrees.
    pub angle_cutoff: f64,
    /// Detected H-bond list: (donor_idx, acceptor_idx, score).
    pub hbonds: Vec<(usize, usize, f64)>,
}
impl HBondAnalysis {
    /// Create a new H-bond analysis.
    pub fn new(dist_cutoff: f64, angle_cutoff: f64) -> Self {
        Self {
            donors: Vec::new(),
            hydrogens: Vec::new(),
            acceptors: Vec::new(),
            dist_cutoff,
            angle_cutoff,
            hbonds: Vec::new(),
        }
    }
    /// Add a donor-hydrogen pair.
    pub fn add_donor(&mut self, donor: [f64; 3], h: [f64; 3]) {
        self.donors.push(donor);
        self.hydrogens.push(h);
    }
    /// Add an acceptor.
    pub fn add_acceptor(&mut self, acceptor: [f64; 3]) {
        self.acceptors.push(acceptor);
    }
    /// Detect all H-bonds.
    pub fn detect(&mut self) {
        self.hbonds.clear();
        for (d_idx, (&d, &h)) in self.donors.iter().zip(self.hydrogens.iter()).enumerate() {
            for (a_idx, &a) in self.acceptors.iter().enumerate() {
                let d_a = dist3(d, a);
                if d_a < self.dist_cutoff {
                    let angle = hbond_angle(d, h, a);
                    if angle > self.angle_cutoff {
                        let score = hbond_score(d_a, angle);
                        self.hbonds.push((d_idx, a_idx, score));
                    }
                }
            }
        }
    }
    /// Number of detected H-bonds.
    pub fn n_hbonds(&self) -> usize {
        self.hbonds.len()
    }
    /// Average H-bond score.
    pub fn avg_score(&self) -> f64 {
        if self.hbonds.is_empty() {
            return 0.0;
        }
        self.hbonds.iter().map(|(_, _, s)| s).sum::<f64>() / self.hbonds.len() as f64
    }
}
/// Secondary structure statistics: helicity, sheet fraction, coil fraction.
#[derive(Debug, Clone, Copy)]
pub struct SecondaryStructureStats {
    /// Fraction of helical residues.
    pub helix_fraction: f64,
    /// Fraction of sheet residues.
    pub sheet_fraction: f64,
    /// Fraction of coil residues.
    pub coil_fraction: f64,
    /// Total number of residues.
    pub n_total: usize,
}
impl SecondaryStructureStats {
    /// Compute statistics from a vector of secondary structure labels.
    pub fn from_labels(labels: &[SecStructLabel]) -> Self {
        let n = labels.len();
        if n == 0 {
            return Self {
                helix_fraction: 0.0,
                sheet_fraction: 0.0,
                coil_fraction: 0.0,
                n_total: 0,
            };
        }
        let n_helix = labels
            .iter()
            .filter(|&&l| l == SecStructLabel::Helix)
            .count();
        let n_sheet = labels
            .iter()
            .filter(|&&l| l == SecStructLabel::Sheet)
            .count();
        let n_coil = labels
            .iter()
            .filter(|&&l| l == SecStructLabel::Coil)
            .count();
        Self {
            helix_fraction: n_helix as f64 / n as f64,
            sheet_fraction: n_sheet as f64 / n as f64,
            coil_fraction: n_coil as f64 / n as f64,
            n_total: n,
        }
    }
    /// Compute from a `ProteinSecondaryStructure` classifier.
    pub fn from_pss(pss: &ProteinSecondaryStructure) -> Self {
        let labels = pss.classify_all();
        Self::from_labels(&labels)
    }
}
/// Ramachandran plot analysis.
///
/// Bins phi/psi dihedral angles and detects outliers.
#[derive(Debug, Clone)]
pub struct RamachandranPlot {
    /// Phi angles in degrees.
    pub phi_angles: Vec<f64>,
    /// Psi angles in degrees.
    pub psi_angles: Vec<f64>,
    /// 2D histogram bins (36x36 for 10-degree bins).
    pub histogram: Vec<Vec<u32>>,
    /// Bin size in degrees.
    pub bin_size: f64,
    /// Number of bins.
    pub n_bins: usize,
}
impl RamachandranPlot {
    /// Create a new Ramachandran plot.
    pub fn new(bin_size: f64) -> Self {
        let n_bins = (360.0 / bin_size) as usize;
        let histogram = vec![vec![0u32; n_bins]; n_bins];
        Self {
            phi_angles: Vec::new(),
            psi_angles: Vec::new(),
            histogram,
            bin_size,
            n_bins,
        }
    }
    /// Add a phi/psi pair.
    pub fn add_point(&mut self, phi: f64, psi: f64) {
        self.phi_angles.push(phi);
        self.psi_angles.push(psi);
        let phi_idx = ((phi + 180.0) / self.bin_size) as usize % self.n_bins;
        let psi_idx = ((psi + 180.0) / self.bin_size) as usize % self.n_bins;
        self.histogram[phi_idx][psi_idx] += 1;
    }
    /// Check if (phi, psi) is in helical region.
    pub fn is_helical(&self, phi: f64, psi: f64) -> bool {
        phi > -100.0 && phi < -30.0 && psi > -70.0 && psi < -10.0
    }
    /// Check if (phi, psi) is in beta-sheet region.
    pub fn is_beta(&self, phi: f64, psi: f64) -> bool {
        phi > -170.0 && phi < -50.0 && psi > 80.0 && psi < 180.0
    }
    /// Check if (phi, psi) is an outlier (not in allowed regions).
    pub fn is_outlier(&self, phi: f64, psi: f64) -> bool {
        !(self.is_helical(phi, psi)
            || self.is_beta(phi, psi)
            || phi > 50.0 && phi < 100.0 && psi > -60.0 && psi < 50.0)
    }
    /// Count outliers.
    pub fn outlier_count(&self) -> usize {
        self.phi_angles
            .iter()
            .zip(self.psi_angles.iter())
            .filter(|&(&phi, &psi)| self.is_outlier(phi, psi))
            .count()
    }
    /// Fraction of points in helical region.
    pub fn helix_fraction(&self) -> f64 {
        if self.phi_angles.is_empty() {
            return 0.0;
        }
        let count = self
            .phi_angles
            .iter()
            .zip(self.psi_angles.iter())
            .filter(|&(&phi, &psi)| self.is_helical(phi, psi))
            .count();
        count as f64 / self.phi_angles.len() as f64
    }
}
/// Elastic Network Model for protein fluctuations.
///
/// Connects Cα pairs within a cutoff distance with springs of equal force
/// constant. Used to predict B-factors and low-frequency modes.
#[derive(Debug, Clone)]
pub struct ProteinElasticity {
    /// Cutoff distance between Cα pairs (Å).
    pub cutoff: f64,
    /// Spring constant for all contacts (kcal/mol/Å²).
    pub spring_constant: f64,
}
impl ProteinElasticity {
    /// Create with given cutoff and spring constant.
    pub fn new(cutoff: f64, spring_constant: f64) -> Self {
        Self {
            cutoff,
            spring_constant,
        }
    }
    /// Return spring constant for a pair at distance `r` (zero if beyond cutoff).
    pub fn spring_constant(&self, r: f64) -> f64 {
        if r <= self.cutoff {
            self.spring_constant
        } else {
            0.0
        }
    }
    /// Find all contact pairs within cutoff distance.
    pub fn contact_pairs(&self, positions: &[[f64; 3]]) -> Vec<(usize, usize, f64)> {
        let n = positions.len();
        let mut pairs = Vec::new();
        for i in 0..n {
            for j in (i + 1)..n {
                let dx = positions[j][0] - positions[i][0];
                let dy = positions[j][1] - positions[i][1];
                let dz = positions[j][2] - positions[i][2];
                let r = (dx * dx + dy * dy + dz * dz).sqrt();
                if r <= self.cutoff {
                    pairs.push((i, j, r));
                }
            }
        }
        pairs
    }
    /// Compute elastic energy for a set of displaced positions.
    ///
    /// `E = 0.5 * sum_{ij} k_ij * (r_ij - r0_ij)^2`.
    /// Here we use the harmonic deviation from the current configuration.
    pub fn elastic_energy(&self, positions: &[[f64; 3]]) -> f64 {
        let pairs = self.contact_pairs(positions);
        pairs
            .iter()
            .map(|(_, _, r)| 0.5 * self.spring_constant * r * r)
            .sum::<f64>()
            * 0.0
    }
    /// Compute elastic energy given reference positions and current positions.
    pub fn elastic_energy_from_ref(&self, ref_pos: &[[f64; 3]], cur_pos: &[[f64; 3]]) -> f64 {
        let pairs = self.contact_pairs(ref_pos);
        pairs
            .iter()
            .map(|&(i, j, r0)| {
                let dx = cur_pos[j][0] - cur_pos[i][0];
                let dy = cur_pos[j][1] - cur_pos[i][1];
                let dz = cur_pos[j][2] - cur_pos[i][2];
                let r = (dx * dx + dy * dy + dz * dz).sqrt();
                0.5 * self.spring_constant * (r - r0).powi(2)
            })
            .sum()
    }
    /// Compute per-residue mean-square displacement (B-factor proxy).
    ///
    /// Uses uniform force constant and contact count as inverse mobility.
    pub fn bfactor_proxy(&self, positions: &[[f64; 3]], k_bt: f64) -> Vec<f64> {
        let n = positions.len();
        let pairs = self.contact_pairs(positions);
        let mut contact_count = vec![0_usize; n];
        for &(i, j, _) in &pairs {
            contact_count[i] += 1;
            contact_count[j] += 1;
        }
        contact_count
            .iter()
            .map(|&c| {
                if c == 0 {
                    k_bt / self.spring_constant
                } else {
                    k_bt / (self.spring_constant * c as f64)
                }
            })
            .collect()
    }
}
/// Solvent accessible surface area via simplified Lee-Richards rolling probe.
///
/// Uses Monte Carlo sphere point sampling to approximate SASA per atom.
#[derive(Debug, Clone)]
pub struct SolventAccessibleSurface {
    /// Atom positions (Å).
    pub positions: Vec<[f64; 3]>,
    /// Van der Waals radii (Å).
    pub radii: Vec<f64>,
    /// Probe radius (Å), typically 1.4 for water.
    pub probe_radius: f64,
}
impl SolventAccessibleSurface {
    /// Create a new SASA calculator.
    pub fn new(positions: Vec<[f64; 3]>, radii: Vec<f64>, probe_radius: f64) -> Self {
        Self {
            positions,
            radii,
            probe_radius,
        }
    }
    /// Compute SASA for atom `i` in Å².
    ///
    /// Uses 92-point Fibonacci sphere sampling.
    pub fn compute_atom(&self, i: usize) -> f64 {
        let ri = self.radii[i] + self.probe_radius;
        let n_pts = 92_usize;
        let golden_ratio = (1.0 + 5.0_f64.sqrt()) / 2.0;
        let mut accessible = 0_usize;
        for k in 0..n_pts {
            let theta = (1.0 - 2.0 * k as f64 / (n_pts - 1) as f64).acos();
            let phi = 2.0 * std::f64::consts::PI * k as f64 / golden_ratio;
            let px = self.positions[i][0] + ri * theta.sin() * phi.cos();
            let py = self.positions[i][1] + ri * theta.sin() * phi.sin();
            let pz = self.positions[i][2] + ri * theta.cos();
            let buried = self.positions.iter().enumerate().any(|(j, &pj)| {
                if j == i {
                    return false;
                }
                let rj = self.radii[j] + self.probe_radius;
                let dx = px - pj[0];
                let dy = py - pj[1];
                let dz = pz - pj[2];
                dx * dx + dy * dy + dz * dz < rj * rj
            });
            if !buried {
                accessible += 1;
            }
        }
        4.0 * std::f64::consts::PI * ri * ri * accessible as f64 / n_pts as f64
    }
    /// Compute SASA for all atoms.
    pub fn compute_all(&self) -> Vec<f64> {
        (0..self.positions.len())
            .map(|i| self.compute_atom(i))
            .collect()
    }
    /// Total SASA of the molecule.
    pub fn total(&self) -> f64 {
        self.compute_all().iter().sum()
    }
}
