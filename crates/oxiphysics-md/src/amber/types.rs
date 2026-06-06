//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::*;

/// Coulomb electrostatic interaction between partial charges.
pub struct CoulombInteraction {
    /// Coulomb constant (kJ mol^-1 nm e^-2).
    /// In AMBER units: 138.935485 kJ mol^-1 nm e^-2.
    pub coulomb_k: f64,
}
impl CoulombInteraction {
    /// Create with AMBER default Coulomb constant.
    pub fn new() -> Self {
        Self {
            coulomb_k: 138.935485,
        }
    }
    /// Coulomb energy between charges q_i and q_j at distance r.
    pub fn energy(&self, q_i: f64, q_j: f64, r: f64) -> f64 {
        self.coulomb_k * q_i * q_j / r
    }
    /// Coulomb force magnitude (positive = repulsive for like charges).
    pub fn force(&self, q_i: f64, q_j: f64, r: f64) -> f64 {
        -self.coulomb_k * q_i * q_j / (r * r)
    }
    /// Force vectors on two charged particles.
    pub fn force_vectors(
        &self,
        q_i: f64,
        q_j: f64,
        ri: [f64; 3],
        rj: [f64; 3],
    ) -> ([f64; 3], [f64; 3]) {
        let dr = sub(rj, ri);
        let r = norm(dr);
        if r < 1e-15 {
            return ([0.0; 3], [0.0; 3]);
        }
        let f_mag = self.force(q_i, q_j, r);
        let fi = scale(dr, f_mag / r);
        let fj = neg(fi);
        (fi, fj)
    }
}
/// AMBER non-bonded exclusion and 1-4 scaling bookkeeper.
///
/// * 1-2 pairs (directly bonded) → fully excluded.
/// * 1-3 pairs (angle neighbours) → fully excluded.
/// * 1-4 pairs (dihedral neighbours) → scaled by `scale_14_elec` and
///   `scale_14_vdw`.
pub struct NonBondedExclusions {
    /// 1-2 bonded pairs (excluded from non-bonded interactions).
    pub pairs_12: Vec<(usize, usize)>,
    /// 1-3 angle pairs (excluded from non-bonded interactions).
    pub pairs_13: Vec<(usize, usize)>,
    /// 1-4 dihedral pairs (scaled non-bonded interactions).
    pub pairs_14: Vec<(usize, usize)>,
    /// Electrostatic 1-4 scaling factor (AMBER default: 1/1.2).
    pub scale_14_elec: f64,
    /// van der Waals 1-4 scaling factor (AMBER default: 0.5).
    pub scale_14_vdw: f64,
}
impl NonBondedExclusions {
    /// Create a new exclusion list with AMBER default 1-4 scaling.
    pub fn new() -> Self {
        Self {
            pairs_12: Vec::new(),
            pairs_13: Vec::new(),
            pairs_14: Vec::new(),
            scale_14_elec: 1.0 / 1.2,
            scale_14_vdw: 0.5,
        }
    }
    /// Add a 1-2 (directly bonded) pair.
    pub fn add_bond(&mut self, i: usize, j: usize) {
        let key = Self::ordered(i, j);
        if !self.pairs_12.contains(&key) {
            self.pairs_12.push(key);
        }
    }
    /// Add a 1-3 (angle) pair (i, k) for the angle i–j–k.
    pub fn add_angle(&mut self, i: usize, _j: usize, k: usize) {
        let key = Self::ordered(i, k);
        if !self.pairs_13.contains(&key) {
            self.pairs_13.push(key);
        }
    }
    /// Add a 1-4 (dihedral) pair (i, l) for the dihedral i–j–k–l.
    pub fn add_dihedral(&mut self, i: usize, _j: usize, _k: usize, l: usize) {
        let key = Self::ordered(i, l);
        if !self.pairs_14.contains(&key) {
            self.pairs_14.push(key);
        }
    }
    /// Returns `true` if the pair (i, j) is fully excluded (1-2 or 1-3).
    pub fn is_excluded(&self, i: usize, j: usize) -> bool {
        let key = Self::ordered(i, j);
        self.pairs_12.contains(&key) || self.pairs_13.contains(&key)
    }
    /// Returns `true` if the pair (i, j) is a 1-4 scaled pair (and not also
    /// a 1-2/1-3 excluded pair).
    pub fn is_14(&self, i: usize, j: usize) -> bool {
        let key = Self::ordered(i, j);
        self.pairs_14.contains(&key) && !self.is_excluded(i, j)
    }
    #[inline]
    fn ordered(i: usize, j: usize) -> (usize, usize) {
        if i <= j { (i, j) } else { (j, i) }
    }
}
/// Representative bond parameter for AMBER99SB.
#[derive(Clone, Debug)]
pub struct Amber99sbBondParam {
    /// First atom type.
    pub type_i: &'static str,
    /// Second atom type.
    pub type_j: &'static str,
    /// Force constant (kJ mol⁻¹ nm⁻²).
    pub k: f64,
    /// Equilibrium bond length (nm).
    pub r0: f64,
}
/// A single Fourier term in the AMBER dihedral potential:
/// Vₙ/2 · (1 + cos(n·φ − γ)).
///
/// `vn` is the barrier height (kJ mol⁻¹), `n` the periodicity, and
/// `gamma` the phase offset (radians).
#[derive(Clone)]
pub struct DihedralTerm {
    /// Barrier height (kJ mol⁻¹).
    pub vn: f64,
    /// Periodicity.
    pub n: u32,
    /// Phase offset (radians).
    pub gamma: f64,
}
/// Parameters for a single atom in a Generalized Born (GB) calculation.
#[derive(Debug, Clone)]
pub struct GbAtom {
    /// Partial charge (elementary charge units, e).
    pub charge: f64,
    /// Effective Born radius (nm).
    pub born_radius: f64,
    /// Van der Waals radius used for SA calculation (nm).
    pub vdw_radius: f64,
}
/// GAFF atom type categories.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum GaffType {
    /// sp3 carbon (gaff: c3)
    C3,
    /// sp2 carbon in a ring (gaff: cc)
    Cc,
    /// sp2 carbon in non-ring (gaff: c2)
    C2,
    /// sp carbon (gaff: c1)
    C1,
    /// carbonyl carbon (gaff: c)
    CarbonylC,
    /// aromatic carbon (gaff: ca)
    Ca,
    /// sp3 nitrogen (gaff: n3)
    N3,
    /// sp2 nitrogen (gaff: n2)
    N2,
    /// aromatic nitrogen (gaff: na)
    Na,
    /// sp3 oxygen (gaff: os or oh)
    O,
    /// carbonyl oxygen (gaff: o)
    CarbonylO,
    /// hydroxyl oxygen (gaff: oh)
    Oh,
    /// sp3 sulfur (gaff: s3)
    S3,
    /// hydrogen bonded to sp3 carbon (gaff: hc)
    Hc,
    /// hydrogen bonded to nitrogen (gaff: hn)
    Hn,
    /// hydrogen bonded to oxygen (gaff: ho)
    Ho,
    /// Other / unknown
    Unknown,
}
/// AMBER force field with Generalized Born and surface area solvation routines.
///
/// This wrapper extends the standard [`AmberForceField`] with implicit solvent
/// energy terms used in the GB/SA model.
#[derive(Default)]
pub struct AmberFF {
    /// The underlying explicit-topology AMBER force field.
    pub ff: AmberForceField,
    /// GB/SA parameters.
    pub gb_params: GbParams,
}
impl AmberFF {
    /// Create a new AmberFF with default GB parameters.
    pub fn new() -> Self {
        Self::default()
    }
    /// Compute the Generalised Born (GB) implicit solvation energy.
    ///
    /// Uses the *Still* GB model:
    ///
    /// ```text
    /// E_GB = −½ · k_e · (1/ε_s − 1/ε_p) · Σ_{i,j} q_i q_j / f_{GB}(r_{ij}, R_i, R_j)
    /// ```
    ///
    /// where the GB function is:
    ///
    /// ```text
    /// f_{GB} = √(r_{ij}² + R_i·R_j·exp(−r_{ij}²/(4·R_i·R_j)))
    /// ```
    ///
    /// # Arguments
    /// * `atoms`     — slice of GB atom parameters (charge, Born radius, …)
    /// * `positions` — Cartesian positions in nm, one per atom
    ///
    /// Returns the GB electrostatic solvation energy in kJ mol⁻¹.
    pub fn compute_generalized_born(&self, atoms: &[GbAtom], positions: &[[f64; 3]]) -> f64 {
        let n = atoms.len().min(positions.len());
        if n == 0 {
            return 0.0;
        }
        let p = &self.gb_params;
        let prefactor = -0.5 * p.coulomb_k * (1.0 / p.eps_solvent - 1.0 / p.eps_solute);
        let mut e_gb = 0.0_f64;
        for i in 0..n {
            for j in 0..n {
                let r2 = if i == j {
                    0.0
                } else {
                    let dx = positions[i][0] - positions[j][0];
                    let dy = positions[i][1] - positions[j][1];
                    let dz = positions[i][2] - positions[j][2];
                    dx * dx + dy * dy + dz * dz
                };
                let ri = atoms[i].born_radius;
                let rj = atoms[j].born_radius;
                let ri_rj = ri * rj;
                let f_gb = (r2 + ri_rj * (-r2 / (4.0 * ri_rj)).exp()).sqrt();
                if f_gb < 1e-15 {
                    continue;
                }
                e_gb += atoms[i].charge * atoms[j].charge / f_gb;
            }
        }
        prefactor * e_gb
    }
    /// Compute the surface area (SA) contribution to GB/SA solvation energy.
    ///
    /// Uses the solvent-accessible surface area (SASA) approximation:
    ///
    /// ```text
    /// E_SA = Σ_i γ · 4π · (R_i + R_probe)²
    /// ```
    ///
    /// where R_probe = 0.14 nm (water probe radius) and γ is the surface tension.
    ///
    /// # Arguments
    /// * `atoms` — slice of GB atoms (only `vdw_radius` is used here)
    ///
    /// Returns the SA energy contribution in kJ mol⁻¹.
    pub fn compute_sa_term(&self, atoms: &[GbAtom]) -> f64 {
        const R_PROBE: f64 = 0.14;
        use std::f64::consts::PI;
        let gamma = self.gb_params.surface_tension;
        let offset = self.gb_params.sa_offset;
        let mut e_sa = offset;
        for atom in atoms {
            let r = atom.vdw_radius + R_PROBE;
            e_sa += gamma * 4.0 * PI * r * r;
        }
        e_sa
    }
    /// Compute the neck GB correction for buried atom pairs.
    ///
    /// The neck correction accounts for the additional exclusion volume
    /// between pairs of overlapping van der Waals spheres.  A simplified
    /// analytical approximation is used:
    ///
    /// ```text
    /// E_neck = −½ · k_e · (1/ε_s − 1/ε_p) · Σ_{i<j} q_i q_j · neck(r_{ij}, R_i, R_j)
    /// ```
    ///
    /// where the neck function is:
    ///
    /// ```text
    /// neck(r, Ri, Rj) = m_neck / (1 + n_neck · (r − d_neck)²)
    /// ```
    ///
    /// with empirical constants m_neck = 0.3516, n_neck = 2.909, and
    /// d_neck = |Ri − Rj| (the offset distance).
    ///
    /// This correction is applied only when r < Ri + Rj (overlapping spheres).
    ///
    /// # Arguments
    /// * `atoms`     — GB atom parameters
    /// * `positions` — atom positions (nm)
    ///
    /// Returns the neck correction energy in kJ mol⁻¹.
    pub fn compute_neck_correction(&self, atoms: &[GbAtom], positions: &[[f64; 3]]) -> f64 {
        let n = atoms.len().min(positions.len());
        if n < 2 {
            return 0.0;
        }
        const M_NECK: f64 = 0.3516;
        const N_NECK: f64 = 2.909;
        let p = &self.gb_params;
        let prefactor = -0.5 * p.coulomb_k * (1.0 / p.eps_solvent - 1.0 / p.eps_solute);
        let mut e_neck = 0.0_f64;
        for i in 0..n {
            for j in i + 1..n {
                let ri = atoms[i].born_radius;
                let rj = atoms[j].born_radius;
                let dx = positions[i][0] - positions[j][0];
                let dy = positions[i][1] - positions[j][1];
                let dz = positions[i][2] - positions[j][2];
                let r = (dx * dx + dy * dy + dz * dz).sqrt();
                if r >= ri + rj {
                    continue;
                }
                let d_neck = (ri - rj).abs();
                let delta = r - d_neck;
                let neck = M_NECK / (1.0 + N_NECK * delta * delta);
                e_neck += atoms[i].charge * atoms[j].charge * neck;
            }
        }
        prefactor * e_neck
    }
}
/// Partial charge record for a single atom in a residue.
#[derive(Clone, Debug)]
pub struct ResidueChargeRecord {
    /// Residue name (e.g., "ALA", "GLY").
    pub residue: &'static str,
    /// Atom name within the residue (e.g., "N", "CA", "C", "O").
    pub atom: &'static str,
    /// Partial charge in elementary-charge units.
    pub charge: f64,
    /// AMBER atom type.
    pub amber_type: &'static str,
}
/// Lorentz-Berthelot combination rules for LJ cross-interactions.
pub struct LorentzBerthelot;
impl LorentzBerthelot {
    /// Combined sigma: sigma_ij = (sigma_i + sigma_j) / 2.
    pub fn sigma(sigma_i: f64, sigma_j: f64) -> f64 {
        (sigma_i + sigma_j) * 0.5
    }
    /// Combined epsilon: epsilon_ij = sqrt(epsilon_i * epsilon_j).
    pub fn epsilon(eps_i: f64, eps_j: f64) -> f64 {
        (eps_i * eps_j).sqrt()
    }
    /// LJ energy between two different atom types at distance r.
    pub fn lj_energy(type_i: &AmberAtomType, type_j: &AmberAtomType, r: f64) -> f64 {
        let sig = Self::sigma(type_i.sigma, type_j.sigma);
        let eps = Self::epsilon(type_i.epsilon, type_j.epsilon);
        let sr6 = (sig / r).powi(6);
        4.0 * eps * (sr6 * sr6 - sr6)
    }
    /// LJ force between two different atom types at distance r.
    pub fn lj_force(type_i: &AmberAtomType, type_j: &AmberAtomType, r: f64) -> f64 {
        let sig = Self::sigma(type_i.sigma, type_j.sigma);
        let eps = Self::epsilon(type_i.epsilon, type_j.epsilon);
        let sr6 = (sig / r).powi(6);
        24.0 * eps / r * (2.0 * sr6 * sr6 - sr6)
    }
}
/// Generalised Born (Still) parameters shared for the whole system.
#[derive(Debug, Clone)]
pub struct GbParams {
    /// Solvent dielectric constant (typically 78.5 for water).
    pub eps_solvent: f64,
    /// Solute dielectric constant (typically 1.0 or 4.0).
    pub eps_solute: f64,
    /// Coulomb constant k_e = 1/(4πε₀) in kJ mol⁻¹ nm e⁻² (138.935 in AMBER units).
    pub coulomb_k: f64,
    /// SASA-based surface tension coefficient γ (kJ mol⁻¹ nm⁻²).
    pub surface_tension: f64,
    /// Offset constant for SA term (kJ mol⁻¹).
    pub sa_offset: f64,
}
/// Harmonic bond potential: V = k·(r − r₀)².
///
/// Force constant `k` in kJ mol⁻¹ nm⁻², equilibrium length `r0` in nm.
/// Note: AMBER convention stores the *full* k (not k/2), so the potential
/// is V = k·(r−r₀)², and the force magnitude is −dV/dr = −2k·(r−r₀).
pub struct HarmonicBond {
    /// Force constant (kJ mol⁻¹ nm⁻²).
    pub k: f64,
    /// Equilibrium bond length (nm).
    pub r0: f64,
}
impl HarmonicBond {
    /// Create a new harmonic bond.
    pub fn new(k: f64, r0: f64) -> Self {
        Self { k, r0 }
    }
    /// Potential energy for distance `r`.
    pub fn energy(&self, r: f64) -> f64 {
        let d = r - self.r0;
        self.k * d * d
    }
    /// Force magnitude (negative gradient): −dV/dr = −2k·(r−r₀).
    pub fn force(&self, r: f64) -> f64 {
        -2.0 * self.k * (r - self.r0)
    }
    /// Force vectors on atoms i and j.
    ///
    /// Returns `(f_i, f_j)` where `f_i = -f_j` (Newton's third law).
    /// When the bond is stretched (`r > r0`) atom i is pulled toward j and
    /// vice-versa.
    pub fn force_vectors(&self, ri: [f64; 3], rj: [f64; 3]) -> ([f64; 3], [f64; 3]) {
        let dr = sub(rj, ri);
        let r = norm(dr);
        if r < 1e-15 {
            return ([0.0; 3], [0.0; 3]);
        }
        let f_mag = -self.force(r);
        let fi = scale(dr, f_mag / r);
        let fj = neg(fi);
        (fi, fj)
    }
}
/// AMBER99SB atom-type record: LJ sigma (nm), epsilon (kJ/mol), typical mass (amu).
///
/// Sources: parm99.dat / ff99SB AMBER parameter file (Cornell et al. 1995,
/// Hornak et al. 2006).
#[derive(Clone, Debug)]
pub struct Amber99sbType {
    /// AMBER atom-type symbol.
    pub symbol: &'static str,
    /// Lennard-Jones sigma in nm (converted from Å: sigma_nm = rmin/2 * 2 / 10).
    pub sigma: f64,
    /// Lennard-Jones well depth epsilon in kJ mol⁻¹.
    pub epsilon: f64,
    /// Typical mass (amu).
    pub mass: f64,
}
/// AMBER proper dihedral: V = Σ Vₙ/2 · (1 + cos(n·φ − γ)).
pub struct ProperDihedral {
    /// Sum of Fourier terms.
    pub terms: Vec<DihedralTerm>,
}
impl ProperDihedral {
    /// Create a new proper dihedral from a list of Fourier terms.
    pub fn new(terms: Vec<DihedralTerm>) -> Self {
        Self { terms }
    }
    /// Compute the dihedral angle φ for atoms i–j–k–l (IUPAC convention).
    ///
    /// Uses atan2 to return an angle in (−π, π].
    pub fn dihedral_angle(ri: [f64; 3], rj: [f64; 3], rk: [f64; 3], rl: [f64; 3]) -> f64 {
        let b1 = sub(rj, ri);
        let b2 = sub(rk, rj);
        let b3 = sub(rl, rk);
        let n1 = cross(b1, b2);
        let n2 = cross(b2, b3);
        let n1_norm = norm(n1);
        let n2_norm = norm(n2);
        if n1_norm < 1e-15 || n2_norm < 1e-15 {
            return 0.0;
        }
        let n1_hat = scale(n1, 1.0 / n1_norm);
        let n2_hat = scale(n2, 1.0 / n2_norm);
        let b2_hat = scale(b2, 1.0 / norm(b2).max(1e-15));
        let m1 = cross(n1_hat, n2_hat);
        let cos_phi = dot(n1_hat, n2_hat).clamp(-1.0, 1.0);
        let sin_phi = dot(m1, b2_hat);
        sin_phi.atan2(cos_phi)
    }
    /// Potential energy for atoms i–j–k–l.
    pub fn energy(&self, ri: [f64; 3], rj: [f64; 3], rk: [f64; 3], rl: [f64; 3]) -> f64 {
        let phi = Self::dihedral_angle(ri, rj, rk, rl);
        self.terms
            .iter()
            .map(|t| 0.5 * t.vn * (1.0 + (t.n as f64 * phi - t.gamma).cos()))
            .sum()
    }
    /// Force vectors on the four atoms i, j, k, l.
    ///
    /// Forces are derived analytically using the Blondel–Karplus method.
    /// Returns `(f_i, f_j, f_k, f_l)`.
    pub fn force_vectors(
        &self,
        ri: [f64; 3],
        rj: [f64; 3],
        rk: [f64; 3],
        rl: [f64; 3],
    ) -> ([f64; 3], [f64; 3], [f64; 3], [f64; 3]) {
        let b1 = sub(rj, ri);
        let b2 = sub(rk, rj);
        let b3 = sub(rl, rk);
        let n1 = cross(b1, b2);
        let n2 = cross(b2, b3);
        let n1_sq = dot(n1, n1);
        let n2_sq = dot(n2, n2);
        if n1_sq < 1e-30 || n2_sq < 1e-30 {
            return ([0.0; 3], [0.0; 3], [0.0; 3], [0.0; 3]);
        }
        let b2_len = norm(b2).max(1e-15);
        let b2_hat = scale(b2, 1.0 / b2_len);
        let phi = Self::dihedral_angle(ri, rj, rk, rl);
        let dv_dphi: f64 = self
            .terms
            .iter()
            .map(|t| {
                let nf = t.n as f64;
                -0.5 * t.vn * nf * (nf * phi - t.gamma).sin()
            })
            .sum();
        let dphi_dri = scale(n1, -b2_len / n1_sq);
        let dphi_drl = scale(n2, b2_len / n2_sq);
        let b1_dot_b2 = dot(b1, b2) / (b2_len * b2_len);
        let b3_dot_b2 = dot(b3, b2) / (b2_len * b2_len);
        let dphi_drj = sub(scale(dphi_dri, b1_dot_b2 - 1.0), scale(dphi_drl, b3_dot_b2));
        let dphi_drk = neg(add(add(dphi_dri, dphi_drj), dphi_drl));
        let _ = b2_hat;
        let fi = scale(dphi_dri, -dv_dphi);
        let fj = scale(dphi_drj, -dv_dphi);
        let fk = scale(dphi_drk, -dv_dphi);
        let fl = scale(dphi_drl, -dv_dphi);
        (fi, fj, fk, fl)
    }
}
/// Representative angle parameter for AMBER99SB.
#[derive(Clone, Debug)]
pub struct Amber99sbAngleParam {
    /// First atom type i.
    pub type_i: &'static str,
    /// Central atom type j.
    pub type_j: &'static str,
    /// Third atom type k.
    pub type_k: &'static str,
    /// Force constant (kJ mol⁻¹ rad⁻²).
    pub k: f64,
    /// Equilibrium angle (degrees).
    pub theta0_deg: f64,
}
/// Improper torsion (out-of-plane bend) using the same Fourier form as
/// proper dihedrals.
///
/// V = Vₙ/2 · (1 + cos(n·φ − γ))
///
/// Typically n = 2, γ = π (phase = 180°) in AMBER, which penalises
/// deviation from planarity.
pub struct ImproperTorsion {
    /// Barrier height (kJ mol⁻¹).
    pub vn: f64,
    /// Periodicity.
    pub n: u32,
    /// Phase offset (radians).
    pub gamma: f64,
}
impl ImproperTorsion {
    /// Create a new improper torsion.
    pub fn new(vn: f64, n: u32, gamma_deg: f64) -> Self {
        Self {
            vn,
            n,
            gamma: gamma_deg.to_radians(),
        }
    }
    /// Energy for the i–j–k–l improper: V = Vₙ/2 · (1 + cos(n·φ − γ)).
    pub fn energy(&self, ri: [f64; 3], rj: [f64; 3], rk: [f64; 3], rl: [f64; 3]) -> f64 {
        let phi = ProperDihedral::dihedral_angle(ri, rj, rk, rl);
        0.5 * self.vn * (1.0 + (self.n as f64 * phi - self.gamma).cos())
    }
    /// dV/dφ for the improper torsion.
    pub fn dv_dphi(&self, phi: f64) -> f64 {
        -0.5 * self.vn * (self.n as f64) * (self.n as f64 * phi - self.gamma).sin()
    }
}
/// Dihedral scan results: energy as a function of dihedral angle.
pub struct DihedralScan {
    /// Scanned angles (radians).
    pub angles: Vec<f64>,
    /// Corresponding energies.
    pub energies: Vec<f64>,
}
impl DihedralScan {
    /// Scan the dihedral energy profile from -pi to +pi with `n_points` steps.
    ///
    /// Uses a simple geometry: atoms placed in a plane with the dihedral
    /// angle varied by rotating atom l around the j-k bond axis.
    pub fn scan(dihedral: &ProperDihedral, n_points: usize) -> Self {
        let mut angles = Vec::with_capacity(n_points);
        let mut energies = Vec::with_capacity(n_points);
        let ri = [0.0, 1.0, 0.0];
        let rj = [0.0, 0.0, 0.0];
        let rk = [1.0, 0.0, 0.0];
        for i in 0..n_points {
            let phi =
                -std::f64::consts::PI + 2.0 * std::f64::consts::PI * (i as f64) / (n_points as f64);
            let rl = [1.0, phi.sin(), -phi.cos()];
            let e = dihedral.energy(ri, rj, rk, rl);
            angles.push(phi);
            energies.push(e);
        }
        Self { angles, energies }
    }
    /// Find the minimum energy and corresponding angle.
    pub fn minimum(&self) -> (f64, f64) {
        let mut min_e = f64::INFINITY;
        let mut min_phi = 0.0;
        for (i, &e) in self.energies.iter().enumerate() {
            if e < min_e {
                min_e = e;
                min_phi = self.angles[i];
            }
        }
        (min_phi, min_e)
    }
    /// Find the maximum energy and corresponding angle.
    pub fn maximum(&self) -> (f64, f64) {
        let mut max_e = f64::NEG_INFINITY;
        let mut max_phi = 0.0;
        for (i, &e) in self.energies.iter().enumerate() {
            if e > max_e {
                max_e = e;
                max_phi = self.angles[i];
            }
        }
        (max_phi, max_e)
    }
    /// Barrier height (max energy - min energy).
    pub fn barrier_height(&self) -> f64 {
        let (_, e_min) = self.minimum();
        let (_, e_max) = self.maximum();
        e_max - e_min
    }
    /// Number of scan points.
    pub fn n_points(&self) -> usize {
        self.angles.len()
    }
}
/// Simple AMBER energy evaluator for a set of bonded interactions.
pub struct AmberEnergyEvaluator {
    /// Bond list: (atom_i, atom_j, HarmonicBond).
    pub bond_list: Vec<(usize, usize, HarmonicBond)>,
    /// Angle list: (atom_i, atom_j, atom_k, HarmonicAngle).
    pub angle_list: Vec<(usize, usize, usize, HarmonicAngle)>,
    /// Dihedral list: (atom_i, atom_j, atom_k, atom_l, ProperDihedral).
    pub dihedral_list: Vec<(usize, usize, usize, usize, ProperDihedral)>,
}
impl AmberEnergyEvaluator {
    /// Create an empty evaluator.
    pub fn new() -> Self {
        Self {
            bond_list: Vec::new(),
            angle_list: Vec::new(),
            dihedral_list: Vec::new(),
        }
    }
    /// Total bond energy.
    pub fn bond_energy(&self, positions: &[[f64; 3]]) -> f64 {
        self.bond_list
            .iter()
            .map(|(i, j, bond)| {
                let r = norm(sub(positions[*j], positions[*i]));
                bond.energy(r)
            })
            .sum()
    }
    /// Total angle energy.
    pub fn angle_energy(&self, positions: &[[f64; 3]]) -> f64 {
        self.angle_list
            .iter()
            .map(|(i, j, k, angle)| angle.energy(positions[*i], positions[*j], positions[*k]))
            .sum()
    }
    /// Total dihedral energy.
    pub fn dihedral_energy(&self, positions: &[[f64; 3]]) -> f64 {
        self.dihedral_list
            .iter()
            .map(|(i, j, k, l, dih)| {
                dih.energy(positions[*i], positions[*j], positions[*k], positions[*l])
            })
            .sum()
    }
    /// Total bonded energy (bonds + angles + dihedrals).
    pub fn total_energy(&self, positions: &[[f64; 3]]) -> f64 {
        self.bond_energy(positions) + self.angle_energy(positions) + self.dihedral_energy(positions)
    }
    /// Compute bond forces and accumulate into force array.
    pub fn accumulate_bond_forces(&self, positions: &[[f64; 3]], forces: &mut [[f64; 3]]) {
        for (i, j, bond) in &self.bond_list {
            let (fi, fj) = bond.force_vectors(positions[*i], positions[*j]);
            forces[*i] = add(forces[*i], fi);
            forces[*j] = add(forces[*j], fj);
        }
    }
    /// Compute angle forces and accumulate into force array.
    pub fn accumulate_angle_forces(&self, positions: &[[f64; 3]], forces: &mut [[f64; 3]]) {
        for (i, j, k, angle) in &self.angle_list {
            let (fi, fj, fk) = angle.force_vectors(positions[*i], positions[*j], positions[*k]);
            forces[*i] = add(forces[*i], fi);
            forces[*j] = add(forces[*j], fj);
            forces[*k] = add(forces[*k], fk);
        }
    }
    /// Number of interactions.
    pub fn num_bonds(&self) -> usize {
        self.bond_list.len()
    }
    /// Number of angle interactions.
    pub fn num_angles(&self) -> usize {
        self.angle_list.len()
    }
    /// Number of dihedral interactions.
    pub fn num_dihedrals(&self) -> usize {
        self.dihedral_list.len()
    }
}
/// Harmonic angle potential: V = k·(θ − θ₀)².
///
/// `k` in kJ mol⁻¹ rad⁻², `theta0` stored internally in radians.
pub struct HarmonicAngle {
    /// Force constant (kJ mol⁻¹ rad⁻²).
    pub k: f64,
    /// Equilibrium angle (radians).
    pub theta0: f64,
}
impl HarmonicAngle {
    /// Create a new harmonic angle.  `theta0_deg` is supplied in degrees and
    /// converted to radians internally.
    pub fn new(k: f64, theta0_deg: f64) -> Self {
        Self {
            k,
            theta0: theta0_deg.to_radians(),
        }
    }
    /// Compute the angle θ for atoms i–j–k where j is the central atom.
    ///
    /// θ = acos( r̂ⱼᵢ · r̂ⱼₖ )
    pub fn angle(ri: [f64; 3], rj: [f64; 3], rk: [f64; 3]) -> f64 {
        let rji = sub(ri, rj);
        let rjk = sub(rk, rj);
        let dji = norm(rji);
        let djk = norm(rjk);
        if dji < 1e-15 || djk < 1e-15 {
            return 0.0;
        }
        let cos_theta = (dot(rji, rjk) / (dji * djk)).clamp(-1.0, 1.0);
        cos_theta.acos()
    }
    /// Potential energy for the i–j–k configuration.
    pub fn energy(&self, ri: [f64; 3], rj: [f64; 3], rk: [f64; 3]) -> f64 {
        let theta = Self::angle(ri, rj, rk);
        let d = theta - self.theta0;
        self.k * d * d
    }
    /// Force vectors on atoms i, j, k.
    ///
    /// Forces are computed via the chain rule: ∂V/∂r = (∂V/∂θ)·(∂θ/∂r).
    /// Returns `(f_i, f_j, f_k)`.
    pub fn force_vectors(
        &self,
        ri: [f64; 3],
        rj: [f64; 3],
        rk: [f64; 3],
    ) -> ([f64; 3], [f64; 3], [f64; 3]) {
        let rji = sub(ri, rj);
        let rjk = sub(rk, rj);
        let dji = norm(rji);
        let djk = norm(rjk);
        if dji < 1e-15 || djk < 1e-15 {
            return ([0.0; 3], [0.0; 3], [0.0; 3]);
        }
        let rji_n = scale(rji, 1.0 / dji);
        let rjk_n = scale(rjk, 1.0 / djk);
        let cos_theta = dot(rji_n, rjk_n).clamp(-1.0, 1.0);
        let theta = cos_theta.acos();
        let sin_theta = (1.0 - cos_theta * cos_theta).sqrt().max(1e-15);
        let dv_dtheta = 2.0 * self.k * (theta - self.theta0);
        let fi_dir = scale(sub(rjk_n, scale(rji_n, cos_theta)), 1.0 / (dji * sin_theta));
        let fi = scale(fi_dir, dv_dtheta);
        let fk_dir = scale(sub(rji_n, scale(rjk_n, cos_theta)), 1.0 / (djk * sin_theta));
        let fk = scale(fk_dir, dv_dtheta);
        let fj = neg(add(fi, fk));
        (fi, fj, fk)
    }
}
/// An AMBER atom type with Lennard-Jones parameters.
#[derive(Clone, Debug)]
pub struct AmberAtomType {
    /// Atom type name (e.g., "CT", "N", "O").
    pub name: String,
    /// Atomic mass (amu).
    pub mass: f64,
    /// Lennard-Jones well depth epsilon (kJ/mol).
    pub epsilon: f64,
    /// Lennard-Jones radius sigma (nm).
    pub sigma: f64,
    /// Partial charge (elementary charge units).
    pub charge: f64,
}
impl AmberAtomType {
    /// Create a new atom type.
    pub fn new(name: &str, mass: f64, epsilon: f64, sigma: f64, charge: f64) -> Self {
        Self {
            name: name.to_string(),
            mass,
            epsilon,
            sigma,
            charge,
        }
    }
    /// Lennard-Jones interaction energy at distance r.
    ///
    /// V_LJ = 4 * epsilon * \[(sigma/r)^12 - (sigma/r)^6\]
    pub fn lj_energy(&self, r: f64) -> f64 {
        let sr6 = (self.sigma / r).powi(6);
        4.0 * self.epsilon * (sr6 * sr6 - sr6)
    }
    /// Lennard-Jones force magnitude (positive = repulsive).
    ///
    /// F_LJ = 24 * epsilon / r * \[2*(sigma/r)^12 - (sigma/r)^6\]
    pub fn lj_force(&self, r: f64) -> f64 {
        let sr6 = (self.sigma / r).powi(6);
        24.0 * self.epsilon / r * (2.0 * sr6 * sr6 - sr6)
    }
    /// Distance at which LJ potential is zero (r = sigma).
    pub fn lj_zero_crossing(&self) -> f64 {
        self.sigma
    }
    /// Minimum energy distance (r_min = 2^(1/6) * sigma).
    pub fn lj_r_min(&self) -> f64 {
        self.sigma * 2.0_f64.powf(1.0 / 6.0)
    }
}
/// A collection of AMBER force-field parameters.
pub struct AmberForceField {
    /// Atom types.
    pub atom_types: Vec<AmberAtomType>,
    /// Bond parameters: (type_i, type_j) -> HarmonicBond.
    pub bonds: Vec<(String, String, HarmonicBond)>,
    /// Angle parameters: (type_i, type_j, type_k) -> HarmonicAngle.
    pub angles: Vec<(String, String, String, HarmonicAngle)>,
    /// Dihedral parameters: (type_i, type_j, type_k, type_l) -> ProperDihedral.
    pub dihedrals: Vec<(String, String, String, String, ProperDihedral)>,
}
impl AmberForceField {
    /// Create an empty force field.
    pub fn new() -> Self {
        Self {
            atom_types: Vec::new(),
            bonds: Vec::new(),
            angles: Vec::new(),
            dihedrals: Vec::new(),
        }
    }
    /// Add an atom type.
    pub fn add_atom_type(&mut self, at: AmberAtomType) {
        self.atom_types.push(at);
    }
    /// Find an atom type by name.
    pub fn find_atom_type(&self, name: &str) -> Option<&AmberAtomType> {
        self.atom_types.iter().find(|at| at.name == name)
    }
    /// Add bond parameters.
    pub fn add_bond_params(&mut self, type_i: &str, type_j: &str, k: f64, r0: f64) {
        self.bonds.push((
            type_i.to_string(),
            type_j.to_string(),
            HarmonicBond::new(k, r0),
        ));
    }
    /// Find bond parameters for a given pair of atom types.
    pub fn find_bond(&self, type_i: &str, type_j: &str) -> Option<&HarmonicBond> {
        self.bonds
            .iter()
            .find(|(ti, tj, _)| (ti == type_i && tj == type_j) || (ti == type_j && tj == type_i))
            .map(|(_, _, b)| b)
    }
    /// Add angle parameters.
    pub fn add_angle_params(
        &mut self,
        type_i: &str,
        type_j: &str,
        type_k: &str,
        k: f64,
        theta0_deg: f64,
    ) {
        self.angles.push((
            type_i.to_string(),
            type_j.to_string(),
            type_k.to_string(),
            HarmonicAngle::new(k, theta0_deg),
        ));
    }
    /// Number of atom types defined.
    pub fn num_atom_types(&self) -> usize {
        self.atom_types.len()
    }
    /// Number of bond types defined.
    pub fn num_bond_types(&self) -> usize {
        self.bonds.len()
    }
}
