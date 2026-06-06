//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use crate::potential::Potential;

use super::functions::PotentialClone;
use super::functions::*;

/// A 1-4 pair entry.
#[derive(Debug, Clone, Copy)]
pub struct OneFourPair {
    /// First atom index.
    pub i: usize,
    /// Fourth atom index.
    pub l: usize,
}
/// CGenFF force field container.
#[derive(Debug, Clone)]
pub struct CgenffForceField {
    /// Atom type definitions.
    pub atom_types: Vec<CgenffAtomType>,
    /// Bonds (CHARMM-style: k, r0).
    pub bonds: Vec<Bond>,
    /// Angles (CHARMM-style: k_theta, theta0).
    pub angles: Vec<AngleTerm>,
    /// Urey-Bradley terms.
    pub urey_bradley: Vec<UreyBradleyTerm>,
    /// Torsion (dihedral) terms.
    pub torsions: Vec<CgenffTorsion>,
    /// Improper dihedrals.
    pub impropers: Vec<ImproperDihedral>,
    /// Maximum charge penalty acceptable for reliable parameters.
    pub charge_penalty_threshold: f64,
    /// Maximum LJ penalty acceptable.
    pub lj_penalty_threshold: f64,
}
impl CgenffForceField {
    /// Create a new empty CGenFF force field.
    pub fn new() -> Self {
        Self {
            atom_types: Vec::new(),
            bonds: Vec::new(),
            angles: Vec::new(),
            urey_bradley: Vec::new(),
            torsions: Vec::new(),
            impropers: Vec::new(),
            charge_penalty_threshold: 10.0,
            lj_penalty_threshold: 50.0,
        }
    }
    /// Add a CGenFF atom type.
    pub fn add_atom_type(
        &mut self,
        name: &str,
        element: &str,
        charge: f64,
        rmin_half: f64,
        epsilon: f64,
        mass: f64,
    ) {
        self.atom_types.push(CgenffAtomType {
            name: name.to_string(),
            element: element.to_string(),
            charge,
            rmin_half,
            epsilon,
            mass,
            charge_penalty: 0.0,
            lj_penalty: 0.0,
        });
    }
    /// Find atom type by name.
    pub fn find_atom_type(&self, name: &str) -> Option<&CgenffAtomType> {
        self.atom_types.iter().find(|t| t.name == name)
    }
    /// Check whether any atom type has a high penalty (requires parameterization review).
    pub fn has_high_penalty_types(&self) -> bool {
        self.atom_types.iter().any(|t| {
            t.charge_penalty > self.charge_penalty_threshold
                || t.lj_penalty > self.lj_penalty_threshold
        })
    }
    /// Compute the CGenFF dihedral energy for a set of Fourier terms.
    ///
    /// V = Σ V_n * (1 + cos(n*phi - delta))
    pub fn cgenff_torsion_energy(phi: f64, terms: &[(f64, u32, f64)]) -> f64 {
        terms
            .iter()
            .map(|(v_n, n, delta)| v_n * (1.0 + (((*n) as f64) * phi - delta).cos()))
            .sum()
    }
    /// Number of defined atom types.
    pub fn n_atom_types(&self) -> usize {
        self.atom_types.len()
    }
}
/// OPLS dihedral term using the Fourier form:
/// V = V1/2*(1+cos(phi)) + V2/2*(1-cos(2*phi)) + V3/2*(1+cos(3*phi)) + V4/2*(1-cos(4*phi))
#[derive(Debug, Clone, Copy)]
pub struct OplsDihedral {
    /// First atom.
    pub i: usize,
    /// Second atom.
    pub j: usize,
    /// Third atom.
    pub k: usize,
    /// Fourth atom.
    pub l: usize,
    /// V1 coefficient (kJ/mol).
    pub v1: f64,
    /// V2 coefficient.
    pub v2: f64,
    /// V3 coefficient.
    pub v3: f64,
    /// V4 coefficient.
    pub v4: f64,
}
/// A bond between two atoms: (atom_i, atom_j, spring_constant, equilibrium_distance).
#[derive(Debug, Clone, Copy)]
pub struct Bond {
    /// First atom index.
    pub i: usize,
    /// Second atom index.
    pub j: usize,
    /// Spring constant.
    pub k: f64,
    /// Equilibrium distance.
    pub r0: f64,
}
/// 1-4 interaction scaling parameters.
///
/// In many force fields, non-bonded interactions between atoms separated
/// by exactly 3 bonds (1-4 pairs) are scaled by constant factors.
#[derive(Debug, Clone, Copy)]
pub struct OneFourScaling {
    /// Scaling factor for LJ 1-4 interactions (typically 0.5 for AMBER).
    pub lj_scale: f64,
    /// Scaling factor for electrostatic 1-4 interactions (typically 1/1.2 for AMBER).
    pub coulomb_scale: f64,
}
impl OneFourScaling {
    /// AMBER-style 1-4 scaling: LJ = 0.5, Coulomb = 1/1.2.
    pub fn amber() -> Self {
        Self {
            lj_scale: 0.5,
            coulomb_scale: 1.0 / 1.2,
        }
    }
    /// OPLS-style 1-4 scaling: LJ = 0.5, Coulomb = 0.5.
    pub fn opls() -> Self {
        Self {
            lj_scale: 0.5,
            coulomb_scale: 0.5,
        }
    }
    /// CHARMM-style: full LJ and Coulomb with special 1-4 LJ parameters.
    pub fn charmm() -> Self {
        Self {
            lj_scale: 1.0,
            coulomb_scale: 1.0,
        }
    }
    /// No scaling (full interactions).
    pub fn full() -> Self {
        Self {
            lj_scale: 1.0,
            coulomb_scale: 1.0,
        }
    }
    /// Scale LJ energy by the 1-4 factor.
    pub fn scale_lj(&self, energy: f64) -> f64 {
        energy * self.lj_scale
    }
    /// Scale Coulomb energy by the 1-4 factor.
    pub fn scale_coulomb(&self, energy: f64) -> f64 {
        energy * self.coulomb_scale
    }
}
/// Lennard-Jones mixing rule types.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum MixingRule {
    /// Lorentz-Berthelot: sigma = (s_i + s_j)/2, epsilon = sqrt(e_i * e_j)
    LorentzBerthelot,
    /// Geometric: sigma = sqrt(s_i * s_j), epsilon = sqrt(e_i * e_j)
    Geometric,
    /// Arithmetic: sigma = (s_i + s_j)/2, epsilon = (e_i + e_j)/2
    Arithmetic,
    /// Waldman-Hagler: sigma = ((s_i^6 + s_j^6)/2)^(1/6), eps = 2*sqrt(e_i*e_j)*s_i^3*s_j^3/(s_i^6+s_j^6)
    WaldmanHagler,
}
/// A harmonic angle between three atoms (i-j-k), with j at the vertex.
#[derive(Debug, Clone, Copy)]
pub struct AngleTerm {
    /// First atom.
    pub i: usize,
    /// Vertex atom.
    pub j: usize,
    /// Third atom.
    pub k: usize,
    /// Force constant.
    pub k_theta: f64,
    /// Equilibrium angle in radians.
    pub theta0: f64,
}
/// CHARMM force field container.
#[derive(Debug, Clone)]
pub struct CharmmForceField {
    /// Atom type definitions.
    pub atom_types: Vec<CharmmAtomType>,
    /// Bonds.
    pub bonds: Vec<Bond>,
    /// Angles.
    pub angles: Vec<AngleTerm>,
    /// Urey-Bradley terms.
    pub urey_bradley: Vec<UreyBradleyTerm>,
    /// Proper dihedrals.
    pub dihedrals: Vec<DihedralTerm>,
    /// Improper dihedrals.
    pub impropers: Vec<ImproperDihedral>,
    /// 1-4 scaling.
    pub scaling_14: OneFourScaling,
}
impl CharmmForceField {
    /// Create a new empty CHARMM force field.
    pub fn new() -> Self {
        Self {
            atom_types: Vec::new(),
            bonds: Vec::new(),
            angles: Vec::new(),
            urey_bradley: Vec::new(),
            dihedrals: Vec::new(),
            impropers: Vec::new(),
            scaling_14: OneFourScaling::charmm(),
        }
    }
    /// Add an atom type.
    pub fn add_atom_type(
        &mut self,
        name: &str,
        element: &str,
        charge: f64,
        rmin_half: f64,
        epsilon: f64,
        mass: f64,
    ) {
        self.atom_types.push(CharmmAtomType {
            name: name.to_string(),
            element: element.to_string(),
            charge,
            rmin_half,
            epsilon,
            mass,
            rmin_half_14: None,
            epsilon_14: None,
        });
    }
    /// Convert CHARMM Rmin/2 to sigma: sigma = 2*Rmin/2 / 2^(1/6).
    pub fn rmin_half_to_sigma(rmin_half: f64) -> f64 {
        2.0 * rmin_half / 2.0_f64.powf(1.0 / 6.0)
    }
    /// Compute Urey-Bradley energy for a term given the 1-3 distance.
    pub fn urey_bradley_energy(r: f64, k_ub: f64, r0: f64) -> f64 {
        let delta = r - r0;
        k_ub * delta * delta
    }
    /// Compute improper dihedral energy.
    pub fn improper_energy(phi: f64, k_imp: f64, phi0: f64) -> f64 {
        let delta = phi - phi0;
        k_imp * delta * delta
    }
    /// Find an atom type by name.
    pub fn find_atom_type(&self, name: &str) -> Option<&CharmmAtomType> {
        self.atom_types.iter().find(|t| t.name == name)
    }
    /// Number of defined atom types.
    pub fn n_atom_types(&self) -> usize {
        self.atom_types.len()
    }
}
/// Severity level for validation issues.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ValidationSeverity {
    /// Non-critical issue.
    Warning,
    /// Critical issue that prevents simulation.
    Error,
}
/// GROMOS atom type parameters (united-atom model).
#[derive(Debug, Clone)]
pub struct GromosAtomType {
    /// Atom type name (e.g. "CH2").
    pub name: String,
    /// Atomic mass (unified atom group mass in Da).
    pub mass: f64,
    /// Partial charge (e).
    pub charge: f64,
    /// C6 parameter (attractive LJ coefficient, kJ·mol^{-1}·nm^6).
    pub c6: f64,
    /// C12 parameter (repulsive LJ coefficient, kJ·mol^{-1}·nm^{12}).
    pub c12: f64,
}
/// CGenFF torsion term (supports multiple multiplicities for one dihedral).
#[derive(Debug, Clone)]
pub struct CgenffTorsion {
    /// Atom type names (4 types).
    pub types: [String; 4],
    /// Fourier terms: (V_n, n, delta).
    pub terms: Vec<(f64, u32, f64)>,
    /// Penalty score for this torsion.
    pub penalty: f64,
}
/// Decomposed intramolecular energy contributions.
#[derive(Debug, Clone, Copy, Default)]
pub struct IntraEnergyDecomposition {
    /// Harmonic bond stretching energy (kJ/mol).
    pub bond_energy: f64,
    /// Harmonic angle bending energy (kJ/mol).
    pub angle_energy: f64,
    /// Proper dihedral (torsion) energy (kJ/mol).
    pub dihedral_energy: f64,
    /// Improper dihedral energy (kJ/mol).
    pub improper_energy: f64,
    /// Urey-Bradley energy (kJ/mol).
    pub urey_bradley_energy: f64,
    /// 1-4 Lennard-Jones energy (kJ/mol).
    pub lj14_energy: f64,
    /// 1-4 Coulomb energy (kJ/mol).
    pub coulomb14_energy: f64,
}
impl IntraEnergyDecomposition {
    /// Total intramolecular energy.
    pub fn total(&self) -> f64 {
        self.bond_energy
            + self.angle_energy
            + self.dihedral_energy
            + self.improper_energy
            + self.urey_bradley_energy
            + self.lj14_energy
            + self.coulomb14_energy
    }
    /// Bonded energy only (bond + angle + dihedral + improper + UB).
    pub fn bonded_total(&self) -> f64 {
        self.bond_energy
            + self.angle_energy
            + self.dihedral_energy
            + self.improper_energy
            + self.urey_bradley_energy
    }
    /// Non-bonded 1-4 energy (LJ14 + Coulomb14).
    pub fn nonbonded14_total(&self) -> f64 {
        self.lj14_energy + self.coulomb14_energy
    }
    /// Fraction of bond energy relative to total (0..1).
    pub fn bond_fraction(&self) -> f64 {
        let t = self.total().abs();
        if t < 1e-30 {
            return 0.0;
        }
        self.bond_energy.abs() / t
    }
}
/// AMBER-like force field combining pair, bond, angle, and dihedral terms.
#[derive(Debug, Clone)]
pub struct AmberForceField {
    /// Pair (non-bonded) force field.
    pub pair: PairForceField,
    /// Bonded force field.
    pub bonded: BondedForceField,
    /// Angle force field.
    pub angle: AngleForceField,
    /// Dihedral force field.
    pub dihedral: DihedralForceField,
}
impl AmberForceField {
    /// Create a new empty AMBER force field.
    pub fn new() -> Self {
        Self {
            pair: PairForceField::new(),
            bonded: BondedForceField::new(),
            angle: AngleForceField::new(),
            dihedral: DihedralForceField::new(),
        }
    }
}
/// General-purpose force field wrapper providing improper dihedral energy,
/// Urey-Bradley 1-3 stretching, and CMAP backbone correction routines.
///
/// These supplement the existing CHARMM/GROMOS/CGenFF containers and can be
/// combined freely with any of those force fields.
#[derive(Debug, Clone, Default)]
pub struct Forcefield {
    /// Name label (e.g. "CHARMM36", "AMBER99SB").
    pub name: String,
}
impl Forcefield {
    /// Create a new empty force field wrapper.
    pub fn new(name: &str) -> Self {
        Self {
            name: name.to_string(),
        }
    }
    /// Compute the improper dihedral energy for a planar-restraint term.
    ///
    /// The CHARMM-style improper dihedral uses a harmonic potential in the
    /// improper angle φ (the dihedral formed by atoms I–J–K–L):
    ///
    /// ```text
    /// E_imp = k_imp * (φ − φ₀)²
    /// ```
    ///
    /// # Arguments
    /// * `ri`, `rj`, `rk`, `rl` — positions of the four atoms (any consistent unit)
    /// * `k_imp` — force constant (energy / rad²)
    /// * `phi0`  — equilibrium improper angle (radians)
    ///
    /// Returns the energy in the same energy unit as `k_imp`.
    pub fn compute_improper_dihedral(
        ri: [f64; 3],
        rj: [f64; 3],
        rk: [f64; 3],
        rl: [f64; 3],
        k_imp: f64,
        phi0: f64,
    ) -> f64 {
        let phi = dihedral_angle(ri, rj, rk, rl);
        let delta = phi - phi0;
        k_imp * delta * delta
    }
    /// Compute Urey-Bradley 1-3 stretching energy.
    ///
    /// The Urey-Bradley term penalises deviation of the 1-3 distance (the
    /// distance between the two end atoms of an angle) from its equilibrium
    /// value:
    ///
    /// ```text
    /// E_UB = k_ub * (r₁₃ − r₁₃⁰)²
    /// ```
    ///
    /// # Arguments
    /// * `r1`, `r3` — positions of atoms 1 and 3 (any consistent unit)
    /// * `k_ub`     — force constant
    /// * `r0`       — equilibrium 1-3 distance
    ///
    /// Returns energy in the same units as `k_ub * distance²`.
    pub fn compute_urey_bradley(r1: [f64; 3], r3: [f64; 3], k_ub: f64, r0: f64) -> f64 {
        let dx = r3[0] - r1[0];
        let dy = r3[1] - r1[1];
        let dz = r3[2] - r1[2];
        let r13 = (dx * dx + dy * dy + dz * dz).sqrt();
        let delta = r13 - r0;
        k_ub * delta * delta
    }
    /// Compute a CMAP backbone correction energy.
    ///
    /// CMAP (Cross-term MAP) provides a 2-D grid-based correction to the
    /// backbone dihedral energy surface.  The correction is tabulated on an
    /// (N×N) grid of (φ, ψ) angles and interpolated bilinearly.
    ///
    /// # Arguments
    /// * `phi`   — backbone φ dihedral (radians, −π … +π)
    /// * `psi`   — backbone ψ dihedral (radians, −π … +π)
    /// * `grid`  — flattened N×N energy grid, row-major.  `grid[i*n + j]`
    ///   is the energy at the i-th φ bin and j-th ψ bin.
    /// * `n`     — grid dimension (must satisfy `grid.len() == n * n` and `n ≥ 2`)
    ///
    /// Returns the bilinearly interpolated CMAP energy at (φ, ψ).
    pub fn compute_cmap_correction(phi: f64, psi: f64, grid: &[f64], n: usize) -> f64 {
        if n < 2 || grid.len() != n * n {
            return 0.0;
        }
        use std::f64::consts::PI;
        let to_frac = |angle: f64| -> f64 {
            let a = ((angle + PI) / (2.0 * PI)) * n as f64;
            a.rem_euclid(n as f64)
        };
        let fi = to_frac(phi);
        let fj = to_frac(psi);
        let i0 = fi.floor() as usize % n;
        let j0 = fj.floor() as usize % n;
        let i1 = (i0 + 1) % n;
        let j1 = (j0 + 1) % n;
        let di = fi - fi.floor();
        let dj = fj - fj.floor();
        let e00 = grid[i0 * n + j0];
        let e01 = grid[i0 * n + j1];
        let e10 = grid[i1 * n + j0];
        let e11 = grid[i1 * n + j1];
        (1.0 - di) * (1.0 - dj) * e00
            + (1.0 - di) * dj * e01
            + di * (1.0 - dj) * e10
            + di * dj * e11
    }
}
/// Validation error for force field parameters.
#[derive(Debug, Clone)]
pub struct ValidationError {
    /// Description of the error.
    pub message: String,
    /// Severity (warning or error).
    pub severity: ValidationSeverity,
}
/// CHARMM Urey-Bradley term: harmonic potential between atoms i and k (1-3 pair).
///
/// V = k_ub * (r_ik - r_ub0)^2
#[derive(Debug, Clone, Copy)]
pub struct UreyBradleyTerm {
    /// First atom.
    pub i: usize,
    /// Third atom (across the angle vertex).
    pub k: usize,
    /// Force constant.
    pub k_ub: f64,
    /// Equilibrium distance.
    pub r0: f64,
}
/// CHARMM atom type parameters.
#[derive(Debug, Clone)]
pub struct CharmmAtomType {
    /// Atom type name (e.g. "CT1").
    pub name: String,
    /// Element symbol.
    pub element: String,
    /// Partial charge in units of e.
    pub charge: f64,
    /// LJ Rmin/2 in Angstrom (half of the LJ minimum distance).
    pub rmin_half: f64,
    /// LJ epsilon in kcal/mol.
    pub epsilon: f64,
    /// Mass in Da.
    pub mass: f64,
    /// Special 1-4 Rmin/2 (CHARMM uses separate 1-4 LJ params).
    pub rmin_half_14: Option<f64>,
    /// Special 1-4 epsilon.
    pub epsilon_14: Option<f64>,
}
/// Pair force field that applies pair potentials using cell list neighbors.
#[derive(Debug, Clone)]
pub struct PairForceField {
    /// List of pair interactions.
    pub interactions: Vec<PairInteraction>,
}
impl PairForceField {
    /// Create a new pair force field.
    pub fn new() -> Self {
        Self {
            interactions: Vec::new(),
        }
    }
    /// Add a pair interaction.
    pub fn add_interaction<P: Potential + Clone + 'static>(
        &mut self,
        type_i: u32,
        type_j: u32,
        potential: P,
    ) {
        self.interactions.push(PairInteraction {
            type_i,
            type_j,
            potential: Box::new(potential),
        });
    }
    /// Find the potential for a pair of atom types.
    pub(super) fn find_potential(&self, ti: u32, tj: u32) -> Option<&dyn PotentialClone> {
        for inter in &self.interactions {
            if (inter.type_i == ti && inter.type_j == tj)
                || (inter.type_i == tj && inter.type_j == ti)
            {
                return Some(inter.potential.as_ref());
            }
        }
        None
    }
    /// Maximum cutoff across all pair interactions.
    pub(super) fn max_cutoff(&self) -> f64 {
        self.interactions
            .iter()
            .map(|i| i.potential.cutoff())
            .fold(0.0_f64, f64::max)
    }
}
/// A pair interaction entry: maps a (type_i, type_j) pair to a potential.
#[derive(Clone)]
pub struct PairInteraction {
    /// Atom type i.
    pub type_i: u32,
    /// Atom type j.
    pub type_j: u32,
    /// The pair potential.
    pub potential: Box<dyn PotentialClone>,
}
/// GROMOS force field container.
#[derive(Debug, Clone)]
pub struct GromosForceField {
    /// Atom type definitions.
    pub atom_types: Vec<GromosAtomType>,
    /// Bonds.
    pub bonds: Vec<Bond>,
    /// Angles.
    pub angles: Vec<AngleTerm>,
    /// Dihedrals.
    pub dihedrals: Vec<DihedralTerm>,
}
impl GromosForceField {
    /// Create a new empty GROMOS force field.
    pub fn new() -> Self {
        Self {
            atom_types: Vec::new(),
            bonds: Vec::new(),
            angles: Vec::new(),
            dihedrals: Vec::new(),
        }
    }
    /// Add a GROMOS atom type.
    pub fn add_atom_type(&mut self, name: &str, mass: f64, charge: f64, c6: f64, c12: f64) {
        self.atom_types.push(GromosAtomType {
            name: name.to_string(),
            mass,
            charge,
            c6,
            c12,
        });
    }
    /// Compute GROMOS LJ energy: V = C12/r^12 - C6/r^6.
    pub fn gromos_lj_energy(r: f64, c6: f64, c12: f64) -> f64 {
        if r < 1e-15 {
            return 0.0;
        }
        let r6 = r.powi(6);
        let r12 = r6 * r6;
        c12 / r12 - c6 / r6
    }
    /// Compute GROMOS LJ force magnitude: dV/dr (negative for attractive).
    pub fn gromos_lj_force(r: f64, c6: f64, c12: f64) -> f64 {
        if r < 1e-15 {
            return 0.0;
        }
        let r7 = r.powi(7);
        let r13 = r.powi(13);
        -12.0 * c12 / r13 + 6.0 * c6 / r7
    }
    /// Convert GROMOS C6/C12 parameters to sigma/epsilon.
    ///
    /// sigma = (C12/C6)^(1/6), epsilon = C6^2 / (4 * C12).
    pub fn c6c12_to_sigma_epsilon(c6: f64, c12: f64) -> (f64, f64) {
        if c6 < 1e-30 || c12 < 1e-30 {
            return (0.0, 0.0);
        }
        let sigma = (c12 / c6).powf(1.0 / 6.0);
        let epsilon = c6 * c6 / (4.0 * c12);
        (sigma, epsilon)
    }
    /// Convert sigma/epsilon to GROMOS C6/C12.
    pub fn sigma_epsilon_to_c6c12(sigma: f64, epsilon: f64) -> (f64, f64) {
        let s6 = sigma.powi(6);
        let c6 = 4.0 * epsilon * s6;
        let c12 = 4.0 * epsilon * s6 * s6;
        (c6, c12)
    }
    /// Find an atom type by name.
    pub fn find_atom_type(&self, name: &str) -> Option<&GromosAtomType> {
        self.atom_types.iter().find(|t| t.name == name)
    }
    /// Number of defined atom types.
    pub fn n_atom_types(&self) -> usize {
        self.atom_types.len()
    }
}
/// Dihedral (torsion) force field.
#[derive(Debug, Clone)]
pub struct DihedralForceField {
    /// List of dihedral terms.
    pub dihedrals: Vec<DihedralTerm>,
}
impl DihedralForceField {
    /// Create a new dihedral force field.
    pub fn new() -> Self {
        Self {
            dihedrals: Vec::new(),
        }
    }
    /// Add a dihedral term.
    pub fn add_dihedral(
        &mut self,
        i: usize,
        j: usize,
        k: usize,
        l: usize,
        v_n: f64,
        gamma: f64,
        n: u32,
    ) {
        self.dihedrals.push(DihedralTerm {
            i,
            j,
            k,
            l,
            v_n,
            gamma,
            n,
        });
    }
}
/// Harmonic bond force field.
#[derive(Debug, Clone)]
pub struct BondedForceField {
    /// List of bonds.
    pub bonds: Vec<Bond>,
}
impl BondedForceField {
    /// Create a new bonded force field.
    pub fn new() -> Self {
        Self { bonds: Vec::new() }
    }
    /// Add a bond.
    pub fn add_bond(&mut self, i: usize, j: usize, k: f64, r0: f64) {
        self.bonds.push(Bond { i, j, k, r0 });
    }
}
/// CHARMM improper dihedral: V = k_imp * (phi - phi0)^2.
#[derive(Debug, Clone, Copy)]
pub struct ImproperDihedral {
    /// Central atom.
    pub i: usize,
    /// Second atom.
    pub j: usize,
    /// Third atom.
    pub k: usize,
    /// Fourth atom.
    pub l: usize,
    /// Force constant.
    pub k_imp: f64,
    /// Equilibrium angle (radians).
    pub phi0: f64,
}
/// CGenFF atom type with CHARMM-style parameters plus CGenFF penalty scores.
#[derive(Debug, Clone)]
pub struct CgenffAtomType {
    /// Atom type name (e.g. "CG2R61").
    pub name: String,
    /// Element symbol.
    pub element: String,
    /// Partial charge (e).
    pub charge: f64,
    /// LJ Rmin/2 (Å).
    pub rmin_half: f64,
    /// LJ epsilon (kcal/mol).
    pub epsilon: f64,
    /// Atomic mass (Da).
    pub mass: f64,
    /// CGenFF penalty score for charge (0 = from analogy, higher = less reliable).
    pub charge_penalty: f64,
    /// CGenFF penalty score for LJ parameters.
    pub lj_penalty: f64,
}
/// A proper dihedral angle term (i-j-k-l).
///
/// V = V_n * (1 + cos(n*phi - gamma))
#[derive(Debug, Clone, Copy)]
pub struct DihedralTerm {
    /// First atom.
    pub i: usize,
    /// Second atom.
    pub j: usize,
    /// Third atom.
    pub k: usize,
    /// Fourth atom.
    pub l: usize,
    /// Barrier height.
    pub v_n: f64,
    /// Phase offset in radians.
    pub gamma: f64,
    /// Multiplicity.
    pub n: u32,
}
/// Harmonic angle force field.
#[derive(Debug, Clone)]
pub struct AngleForceField {
    /// List of angle terms.
    pub angles: Vec<AngleTerm>,
}
impl AngleForceField {
    /// Create a new angle force field.
    pub fn new() -> Self {
        Self { angles: Vec::new() }
    }
    /// Add an angle term.
    pub fn add_angle(&mut self, i: usize, j: usize, k: usize, k_theta: f64, theta0: f64) {
        self.angles.push(AngleTerm {
            i,
            j,
            k,
            k_theta,
            theta0,
        });
    }
}
/// OPLS-AA atom type parameters.
#[derive(Debug, Clone)]
pub struct OplsAtomType {
    /// Atom type name (e.g. "opls_135").
    pub name: String,
    /// Element symbol.
    pub element: String,
    /// Partial charge in units of e.
    pub charge: f64,
    /// LJ sigma in nm.
    pub sigma: f64,
    /// LJ epsilon in kJ/mol.
    pub epsilon: f64,
    /// Atomic mass in Da.
    pub mass: f64,
}
/// LJ atom type parameters.
#[derive(Debug, Clone, Copy)]
pub struct LjAtomType {
    /// Atom type identifier.
    pub type_id: u32,
    /// LJ sigma parameter.
    pub sigma: f64,
    /// LJ epsilon parameter.
    pub epsilon: f64,
}
/// OPLS force field container.
#[derive(Debug, Clone)]
pub struct OplsForceField {
    /// Atom type definitions.
    pub atom_types: Vec<OplsAtomType>,
    /// Bonds (harmonic).
    pub bonds: Vec<Bond>,
    /// Angles (harmonic).
    pub angles: Vec<AngleTerm>,
    /// Dihedrals (Fourier series with OPLS functional form).
    pub dihedrals: Vec<OplsDihedral>,
    /// 1-4 scaling factors.
    pub scaling_14: OneFourScaling,
}
impl OplsForceField {
    /// Create a new empty OPLS force field.
    pub fn new() -> Self {
        Self {
            atom_types: Vec::new(),
            bonds: Vec::new(),
            angles: Vec::new(),
            dihedrals: Vec::new(),
            scaling_14: OneFourScaling::opls(),
        }
    }
    /// Add an atom type.
    pub fn add_atom_type(
        &mut self,
        name: &str,
        element: &str,
        charge: f64,
        sigma: f64,
        epsilon: f64,
        mass: f64,
    ) {
        self.atom_types.push(OplsAtomType {
            name: name.to_string(),
            element: element.to_string(),
            charge,
            sigma,
            epsilon,
            mass,
        });
    }
    /// Compute the OPLS dihedral energy for one term.
    pub fn opls_dihedral_energy(phi: f64, v1: f64, v2: f64, v3: f64, v4: f64) -> f64 {
        v1 / 2.0 * (1.0 + phi.cos())
            + v2 / 2.0 * (1.0 - (2.0 * phi).cos())
            + v3 / 2.0 * (1.0 + (3.0 * phi).cos())
            + v4 / 2.0 * (1.0 - (4.0 * phi).cos())
    }
    /// Find an atom type by name.
    pub fn find_atom_type(&self, name: &str) -> Option<&OplsAtomType> {
        self.atom_types.iter().find(|t| t.name == name)
    }
    /// Number of defined atom types.
    pub fn n_atom_types(&self) -> usize {
        self.atom_types.len()
    }
}
