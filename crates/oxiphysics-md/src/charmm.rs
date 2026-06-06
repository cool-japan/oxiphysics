// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! CHARMM force field implementation.
//!
//! Provides harmonic bonds, Urey-Bradley angles, n-fold dihedral terms,
//! improper dihedrals, and van der Waals (vdW) parameters in the CHARMM
//! parameterisation convention.
//!
//! # Unit conventions
//! - Energies: kcal/mol
//! - Lengths: Å (Angstroms)
//! - Angles: radians (internally); parameters stored in radians
//! - Force constants: kcal/(mol·Å²) for bonds, kcal/(mol·rad²) for angles

use std::collections::HashMap;

// ── Constants ─────────────────────────────────────────────────────────────────

const PI: f64 = std::f64::consts::PI;

// ── CharmmBondParams ──────────────────────────────────────────────────────────

/// Harmonic bond parameters for CHARMM.
///
/// V_bond = k * (r - r0)²
#[derive(Debug, Clone, PartialEq)]
pub struct CharmmBondParams {
    /// Force constant (kcal/mol/Å²).
    pub k: f64,
    /// Equilibrium bond length (Å).
    pub r0: f64,
}

impl CharmmBondParams {
    /// Create new bond parameters.
    pub fn new(k: f64, r0: f64) -> Self {
        Self { k, r0 }
    }
}

// ── CharmmAngleParams ─────────────────────────────────────────────────────────

/// Angle parameters for CHARMM, including optional Urey-Bradley term.
///
/// V_angle = k * (θ - θ₀)²  +  k_UB * (r_13 - r_UB₀)²
#[derive(Debug, Clone, PartialEq)]
pub struct CharmmAngleParams {
    /// Harmonic angle force constant (kcal/mol/rad²).
    pub k: f64,
    /// Equilibrium angle (radians).
    pub theta0: f64,
    /// Urey-Bradley force constant (kcal/mol/Å²); 0 if absent.
    pub urey_bradley_k: f64,
    /// Urey-Bradley equilibrium 1,3-distance (Å); 0 if absent.
    pub urey_bradley_r0: f64,
}

impl CharmmAngleParams {
    /// Create new angle parameters with optional Urey-Bradley correction.
    pub fn new(k: f64, theta0: f64, urey_bradley_k: f64, urey_bradley_r0: f64) -> Self {
        Self {
            k,
            theta0,
            urey_bradley_k,
            urey_bradley_r0,
        }
    }

    /// Create angle parameters without a Urey-Bradley term.
    pub fn simple(k: f64, theta0: f64) -> Self {
        Self::new(k, theta0, 0.0, 0.0)
    }
}

// ── CharmmDihedralParams ──────────────────────────────────────────────────────

/// N-fold (cosine series) dihedral parameters for CHARMM.
///
/// V_dihedral = k * (1 + cos(n·φ − δ))
#[derive(Debug, Clone, PartialEq)]
pub struct CharmmDihedralParams {
    /// Dihedral force constant (kcal/mol).
    pub k: f64,
    /// Periodicity (integer multiplicity).
    pub n: u32,
    /// Phase shift δ (radians).
    pub delta: f64,
}

impl CharmmDihedralParams {
    /// Create new dihedral parameters.
    pub fn new(k: f64, n: u32, delta: f64) -> Self {
        Self { k, n, delta }
    }
}

// ── CharmmImproperParams ──────────────────────────────────────────────────────

/// Improper dihedral parameters for CHARMM.
///
/// V_improper = k * (ψ − ψ₀)²
#[derive(Debug, Clone, PartialEq)]
pub struct CharmmImproperParams {
    /// Improper force constant (kcal/mol/rad²).
    pub k: f64,
    /// Equilibrium improper angle ψ₀ (radians).
    pub psi0: f64,
}

impl CharmmImproperParams {
    /// Create new improper parameters.
    pub fn new(k: f64, psi0: f64) -> Self {
        Self { k, psi0 }
    }
}

// ── CharmmVdwParams ───────────────────────────────────────────────────────────

/// Lennard-Jones vdW parameters in CHARMM notation (ε, Rmin/2).
///
/// The full pair potential uses the combination rules:
///   ε_ij = sqrt(ε_i · ε_j),  Rmin_ij = Rmin/2_i + Rmin/2_j
///
/// V_vdW(r) = ε_ij * \[(Rmin_ij/r)^12 − 2·(Rmin_ij/r)^6\]
#[derive(Debug, Clone, PartialEq)]
pub struct CharmmVdwParams {
    /// Well depth ε (kcal/mol, positive).
    pub epsilon: f64,
    /// Half the minimum-energy radius Rmin/2 (Å).
    pub rmin_half: f64,
}

impl CharmmVdwParams {
    /// Create new vdW parameters.
    pub fn new(epsilon: f64, rmin_half: f64) -> Self {
        Self { epsilon, rmin_half }
    }
}

// ── CharmmForceField ──────────────────────────────────────────────────────────

/// CHARMM force field parameter database.
///
/// Stores maps from atom-type string pairs (or quadruplets) to the
/// corresponding CHARMM parameters.  Atom types follow the CHARMM naming
/// convention (e.g., `"CT1"`, `"CT2"`, `"NH1"`, …).
#[derive(Debug, Default, Clone)]
pub struct CharmmForceField {
    /// Bond parameters keyed by sorted atom-type pair `(type_a, type_b)`.
    pub bonds: HashMap<(String, String), CharmmBondParams>,
    /// Angle parameters keyed by atom-type triple `(type_a, type_b, type_c)`.
    pub angles: HashMap<(String, String, String), CharmmAngleParams>,
    /// Dihedral parameters keyed by atom-type quadruple.
    /// Multiple terms (different n) are supported via `Vec`.
    pub dihedrals: HashMap<(String, String, String, String), Vec<CharmmDihedralParams>>,
    /// Improper parameters keyed by atom-type quadruple.
    pub impropers: HashMap<(String, String, String, String), CharmmImproperParams>,
    /// Atom vdW parameters keyed by atom type.
    pub vdw: HashMap<String, CharmmVdwParams>,
}

impl CharmmForceField {
    /// Create an empty force field ready to receive parameters.
    pub fn new() -> Self {
        Self::default()
    }

    /// Add or overwrite bond parameters for the pair (type_a, type_b).
    ///
    /// The pair is stored in sorted order so lookup is order-independent.
    pub fn add_bond(&mut self, type_a: &str, type_b: &str, params: CharmmBondParams) {
        let key = sorted_pair(type_a, type_b);
        self.bonds.insert(key, params);
    }

    /// Look up bond parameters (order-independent).
    pub fn get_bond(&self, type_a: &str, type_b: &str) -> Option<&CharmmBondParams> {
        let key = sorted_pair(type_a, type_b);
        self.bonds.get(&key)
    }

    /// Add angle parameters.
    pub fn add_angle(&mut self, ta: &str, tb: &str, tc: &str, params: CharmmAngleParams) {
        let key = (ta.to_string(), tb.to_string(), tc.to_string());
        self.angles.insert(key, params);
    }

    /// Look up angle parameters (tries original and reversed order).
    pub fn get_angle(&self, ta: &str, tb: &str, tc: &str) -> Option<&CharmmAngleParams> {
        let key = (ta.to_string(), tb.to_string(), tc.to_string());
        if let Some(p) = self.angles.get(&key) {
            return Some(p);
        }
        let rev_key = (tc.to_string(), tb.to_string(), ta.to_string());
        self.angles.get(&rev_key)
    }

    /// Add one dihedral term (appends to existing terms for that quadruple).
    pub fn add_dihedral(
        &mut self,
        ta: &str,
        tb: &str,
        tc: &str,
        td: &str,
        params: CharmmDihedralParams,
    ) {
        let key = (
            ta.to_string(),
            tb.to_string(),
            tc.to_string(),
            td.to_string(),
        );
        self.dihedrals.entry(key).or_default().push(params);
    }

    /// Look up all dihedral terms for a quadruple (tries forward and reversed).
    pub fn get_dihedral(
        &self,
        ta: &str,
        tb: &str,
        tc: &str,
        td: &str,
    ) -> Option<&Vec<CharmmDihedralParams>> {
        let key = (
            ta.to_string(),
            tb.to_string(),
            tc.to_string(),
            td.to_string(),
        );
        if let Some(p) = self.dihedrals.get(&key) {
            return Some(p);
        }
        let rev = (
            td.to_string(),
            tc.to_string(),
            tb.to_string(),
            ta.to_string(),
        );
        self.dihedrals.get(&rev)
    }

    /// Add improper dihedral parameters.
    pub fn add_improper(
        &mut self,
        ta: &str,
        tb: &str,
        tc: &str,
        td: &str,
        params: CharmmImproperParams,
    ) {
        let key = (
            ta.to_string(),
            tb.to_string(),
            tc.to_string(),
            td.to_string(),
        );
        self.impropers.insert(key, params);
    }

    /// Look up improper parameters.
    pub fn get_improper(
        &self,
        ta: &str,
        tb: &str,
        tc: &str,
        td: &str,
    ) -> Option<&CharmmImproperParams> {
        let key = (
            ta.to_string(),
            tb.to_string(),
            tc.to_string(),
            td.to_string(),
        );
        self.impropers.get(&key)
    }

    /// Add vdW parameters for an atom type.
    pub fn add_vdw(&mut self, atom_type: &str, params: CharmmVdwParams) {
        self.vdw.insert(atom_type.to_string(), params);
    }

    /// Look up vdW parameters for an atom type.
    pub fn get_vdw(&self, atom_type: &str) -> Option<&CharmmVdwParams> {
        self.vdw.get(atom_type)
    }

    /// Compute combined vdW parameters for an i-j pair using CHARMM combining rules.
    ///
    /// ε_ij = √(ε_i · ε_j),  Rmin_ij = Rmin/2_i + Rmin/2_j
    ///
    /// Returns `None` if either atom type is not found.
    pub fn combined_vdw(&self, ta: &str, tb: &str) -> Option<(f64, f64)> {
        let pi = self.vdw.get(ta)?;
        let pj = self.vdw.get(tb)?;
        let epsilon_ij = (pi.epsilon * pj.epsilon).sqrt();
        let rmin_ij = pi.rmin_half + pj.rmin_half;
        Some((epsilon_ij, rmin_ij))
    }

    /// Build a minimal CHARMM-like parameter set for a backbone peptide unit.
    ///
    /// Atom types used: CT1, CT2, CT3, NH1, C, O, HN
    /// (representative CHARMM22 values, simplified).
    pub fn charmm22_backbone() -> Self {
        let mut ff = Self::new();

        // Bonds (k in kcal/mol/Å², r0 in Å)
        ff.add_bond("CT1", "C", CharmmBondParams::new(250.0, 1.522));
        ff.add_bond("CT2", "C", CharmmBondParams::new(250.0, 1.522));
        ff.add_bond("CT3", "C", CharmmBondParams::new(250.0, 1.507));
        ff.add_bond("CT1", "NH1", CharmmBondParams::new(320.0, 1.430));
        ff.add_bond("CT2", "NH1", CharmmBondParams::new(320.0, 1.430));
        ff.add_bond("C", "NH1", CharmmBondParams::new(370.0, 1.345));
        ff.add_bond("C", "O", CharmmBondParams::new(620.0, 1.230));
        ff.add_bond("NH1", "HN", CharmmBondParams::new(480.0, 1.000));
        ff.add_bond("CT1", "CT2", CharmmBondParams::new(222.5, 1.538));
        ff.add_bond("CT2", "CT3", CharmmBondParams::new(222.5, 1.528));

        // Angles (k in kcal/mol/rad², θ₀ in radians)
        let deg = |d: f64| d * PI / 180.0;
        ff.add_angle(
            "NH1",
            "CT1",
            "C",
            CharmmAngleParams::simple(50.0, deg(111.5)),
        );
        ff.add_angle(
            "NH1",
            "CT2",
            "C",
            CharmmAngleParams::simple(50.0, deg(111.5)),
        );
        ff.add_angle(
            "CT1",
            "C",
            "NH1",
            CharmmAngleParams::simple(80.0, deg(116.5)),
        );
        ff.add_angle(
            "CT2",
            "C",
            "NH1",
            CharmmAngleParams::simple(80.0, deg(116.5)),
        );
        ff.add_angle("CT1", "C", "O", CharmmAngleParams::simple(80.0, deg(121.0)));
        ff.add_angle("CT2", "C", "O", CharmmAngleParams::simple(80.0, deg(121.0)));
        ff.add_angle(
            "C",
            "NH1",
            "HN",
            CharmmAngleParams::simple(30.0, deg(119.8)),
        );
        ff.add_angle(
            "C",
            "NH1",
            "CT1",
            CharmmAngleParams::simple(50.0, deg(121.0)),
        );
        ff.add_angle(
            "C",
            "NH1",
            "CT2",
            CharmmAngleParams::simple(50.0, deg(121.0)),
        );

        // Dihedrals  (phi: C-NH1-CT1-C,  psi: NH1-CT1-C-NH1)
        ff.add_dihedral(
            "C",
            "NH1",
            "CT1",
            "C",
            CharmmDihedralParams::new(0.2, 1, deg(0.0)),
        );
        ff.add_dihedral(
            "C",
            "NH1",
            "CT1",
            "C",
            CharmmDihedralParams::new(0.4, 3, deg(0.0)),
        );
        ff.add_dihedral(
            "NH1",
            "CT1",
            "C",
            "NH1",
            CharmmDihedralParams::new(0.2, 1, deg(0.0)),
        );
        ff.add_dihedral(
            "NH1",
            "CT1",
            "C",
            "O",
            CharmmDihedralParams::new(0.4, 2, deg(180.0)),
        );

        // Improper (planarity of amide group)
        ff.add_improper(
            "C",
            "CT2",
            "NH1",
            "O",
            CharmmImproperParams::new(10.5, deg(0.0)),
        );

        // vdW parameters (ε in kcal/mol, Rmin/2 in Å)
        ff.add_vdw("C", CharmmVdwParams::new(0.0700, 2.0000));
        ff.add_vdw("CT1", CharmmVdwParams::new(0.0200, 2.2750));
        ff.add_vdw("CT2", CharmmVdwParams::new(0.0550, 2.1750));
        ff.add_vdw("CT3", CharmmVdwParams::new(0.0800, 2.0600));
        ff.add_vdw("NH1", CharmmVdwParams::new(0.2000, 1.8500));
        ff.add_vdw("O", CharmmVdwParams::new(0.1200, 1.7000));
        ff.add_vdw("HN", CharmmVdwParams::new(0.0046, 0.9000));

        ff
    }
}

// ── Force and energy computation functions ────────────────────────────────────

/// Compute bond energy and force magnitude from distance and parameters.
///
/// V = k·(r − r₀)²
/// F = −dV/dr = −2·k·(r − r₀)   (positive = repulsive, negative = attractive)
///
/// Returns `(energy, force_magnitude)`.  A positive force magnitude means
/// the force along the bond vector points from atom j toward atom i (restoring).
pub fn compute_bond_force(r: f64, params: &CharmmBondParams) -> (f64, f64) {
    let d = r - params.r0;
    let energy = params.k * d * d;
    let force_mag = -2.0 * params.k * d;
    (energy, force_mag)
}

/// Compute angle energy and force (torque) from angle θ and parameters.
///
/// V_harmonic = k·(θ − θ₀)²
/// F = −dV/dθ = −2·k·(θ − θ₀)   (in kcal/mol/rad)
///
/// The Urey-Bradley contribution is NOT included here (it requires the 1-3
/// distance separately); use [`compute_urey_bradley`] for that term.
///
/// Returns `(energy, force_angular)`.
pub fn compute_angle_force(theta: f64, params: &CharmmAngleParams) -> (f64, f64) {
    let d = theta - params.theta0;
    let energy = params.k * d * d;
    let force_angular = -2.0 * params.k * d;
    (energy, force_angular)
}

/// Compute the Urey-Bradley energy and force magnitude.
///
/// V_UB = k_UB·(r₁₃ − r_UB₀)²
///
/// Returns `(energy, force_magnitude)` analogous to [`compute_bond_force`].
pub fn compute_urey_bradley(r13: f64, params: &CharmmAngleParams) -> (f64, f64) {
    if params.urey_bradley_k == 0.0 {
        return (0.0, 0.0);
    }
    let d = r13 - params.urey_bradley_r0;
    let energy = params.urey_bradley_k * d * d;
    let force_mag = -2.0 * params.urey_bradley_k * d;
    (energy, force_mag)
}

/// Compute the combined angle + Urey-Bradley energy and angular force.
///
/// Requires both the angle `theta` and the 1-3 distance `r13`.
/// Returns `(total_energy, angular_force_in_kcal_per_mol_per_rad)`.
pub fn compute_angle_and_ub_force(theta: f64, r13: f64, params: &CharmmAngleParams) -> (f64, f64) {
    let (e_ang, f_ang) = compute_angle_force(theta, params);
    let (e_ub, _) = compute_urey_bradley(r13, params);
    (e_ang + e_ub, f_ang)
}

/// Compute dihedral energy and force from angle φ and CHARMM parameters.
///
/// V = k·(1 + cos(n·φ − δ))
/// F = −dV/dφ = k·n·sin(n·φ − δ)   (in kcal/mol/rad)
///
/// Returns `(energy, force_angular)`.
pub fn compute_dihedral_force(phi: f64, params: &CharmmDihedralParams) -> (f64, f64) {
    let arg = params.n as f64 * phi - params.delta;
    let energy = params.k * (1.0 + arg.cos());
    let force_angular = params.k * params.n as f64 * arg.sin();
    (energy, force_angular)
}

/// Compute total dihedral energy and force for multiple CHARMM terms.
///
/// Returns `(total_energy, total_angular_force)`.
pub fn compute_dihedral_force_multi(phi: f64, terms: &[CharmmDihedralParams]) -> (f64, f64) {
    let (mut e, mut f) = (0.0, 0.0);
    for term in terms {
        let (de, df) = compute_dihedral_force(phi, term);
        e += de;
        f += df;
    }
    (e, f)
}

/// Compute improper dihedral energy and force.
///
/// V = k·(ψ − ψ₀)²
/// F = −dV/dψ = −2·k·(ψ − ψ₀)   (in kcal/mol/rad)
///
/// Returns `(energy, force_angular)`.
pub fn compute_improper_force(psi: f64, params: &CharmmImproperParams) -> (f64, f64) {
    let d = psi - params.psi0;
    let energy = params.k * d * d;
    let force_angular = -2.0 * params.k * d;
    (energy, force_angular)
}

/// Compute Lennard-Jones vdW energy between a pair of atoms given combined ε
/// and Rmin.
///
/// V_LJ = ε * \[(Rmin/r)^12 − 2·(Rmin/r)^6\]
///
/// This form has a minimum of −ε at r = Rmin.
///
/// Returns `(energy, force_magnitude)` where a positive `force_magnitude`
/// means the force along the interatomic vector is repulsive.
pub fn compute_vdw_force(r: f64, epsilon: f64, rmin: f64) -> (f64, f64) {
    if r < 1e-10 {
        return (f64::MAX, f64::MAX);
    }
    let ratio = rmin / r;
    let r6 = ratio.powi(6);
    let r12 = r6 * r6;
    let energy = epsilon * (r12 - 2.0 * r6);
    // dV/dr = ε * (-12·Rmin^12/r^13 + 12·Rmin^6/r^7)
    //       = 12·ε/r * (-(Rmin/r)^12 + (Rmin/r)^6)
    let force_mag = 12.0 * epsilon / r * (r6 - r12); // positive → repulsive
    (energy, force_mag)
}

/// Evaluate the full bonded contribution for a list of bonds.
///
/// `bond_list` contains tuples `(i, j, params)`.
/// `positions` is a flat array of atom positions.
///
/// Returns `(total_bond_energy, forces)` where `forces[i]` is the force
/// vector on atom `i`.
pub fn compute_bonded_forces(
    positions: &[[f64; 3]],
    bond_list: &[(usize, usize, CharmmBondParams)],
) -> (f64, Vec<[f64; 3]>) {
    let n = positions.len();
    let mut forces = vec![[0.0f64; 3]; n];
    let mut total_energy = 0.0;

    for (i, j, params) in bond_list {
        let ri = positions[*i];
        let rj = positions[*j];
        let dx = rj[0] - ri[0];
        let dy = rj[1] - ri[1];
        let dz = rj[2] - ri[2];
        let r = (dx * dx + dy * dy + dz * dz).sqrt();
        if r < 1e-15 {
            continue;
        }
        let (energy, force_mag) = compute_bond_force(r, params);
        total_energy += energy;

        // Force on i points toward j when bond is stretched (force_mag < 0 → attractive)
        let fx = -force_mag * dx / r;
        let fy = -force_mag * dy / r;
        let fz = -force_mag * dz / r;
        forces[*i][0] += fx;
        forces[*i][1] += fy;
        forces[*i][2] += fz;
        forces[*j][0] -= fx;
        forces[*j][1] -= fy;
        forces[*j][2] -= fz;
    }

    (total_energy, forces)
}

/// Compute the valence angle θ between atoms i-j-k (j is the vertex).
pub fn compute_angle(positions: &[[f64; 3]], i: usize, j: usize, k: usize) -> f64 {
    let ri = positions[i];
    let rj = positions[j];
    let rk = positions[k];

    let u = [ri[0] - rj[0], ri[1] - rj[1], ri[2] - rj[2]];
    let v = [rk[0] - rj[0], rk[1] - rj[1], rk[2] - rj[2]];

    let u_norm = (u[0] * u[0] + u[1] * u[1] + u[2] * u[2]).sqrt();
    let v_norm = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();

    if u_norm < 1e-15 || v_norm < 1e-15 {
        return 0.0;
    }

    let cos_theta = (u[0] * v[0] + u[1] * v[1] + u[2] * v[2]) / (u_norm * v_norm);
    cos_theta.clamp(-1.0, 1.0).acos()
}

/// Evaluate angle bending forces for a list of angle terms.
///
/// `angle_list` contains tuples `(i, j, k, params)`.
///
/// Returns `(total_angle_energy, forces)`.
pub fn compute_angle_forces(
    positions: &[[f64; 3]],
    angle_list: &[(usize, usize, usize, CharmmAngleParams)],
) -> (f64, Vec<[f64; 3]>) {
    let n = positions.len();
    let mut forces = vec![[0.0f64; 3]; n];
    let mut total_energy = 0.0;

    for (ai, aj, ak, params) in angle_list {
        let theta = compute_angle(positions, *ai, *aj, *ak);
        let (energy, f_angular) = compute_angle_force(theta, params);
        total_energy += energy;

        // Urey-Bradley
        let ri = positions[*ai];
        let rk = positions[*ak];
        let dr13 = [rk[0] - ri[0], rk[1] - ri[1], rk[2] - ri[2]];
        let r13 = (dr13[0] * dr13[0] + dr13[1] * dr13[1] + dr13[2] * dr13[2]).sqrt();
        let (e_ub, f_ub) = compute_urey_bradley(r13, params);
        total_energy += e_ub;

        // Analytic gradient of angle with respect to atom positions
        let rj = positions[*aj];
        let u = [ri[0] - rj[0], ri[1] - rj[1], ri[2] - rj[2]];
        let v = [rk[0] - rj[0], rk[1] - rj[1], rk[2] - rj[2]];
        let u_norm = (u[0] * u[0] + u[1] * u[1] + u[2] * u[2]).sqrt();
        let v_norm = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();

        if u_norm > 1e-15 && v_norm > 1e-15 {
            let cos_t = (u[0] * v[0] + u[1] * v[1] + u[2] * v[2]) / (u_norm * v_norm);
            let cos_t = cos_t.clamp(-1.0, 1.0);
            let sin_t = (1.0 - cos_t * cos_t).sqrt().max(1e-15);

            let u_hat = [u[0] / u_norm, u[1] / u_norm, u[2] / u_norm];
            let v_hat = [v[0] / v_norm, v[1] / v_norm, v[2] / v_norm];

            // dV/dri = (dV/dθ) * (dθ/dri)
            // dθ/dri = -(v_hat - cos_t * u_hat) / (u_norm * sin_t)
            let scale_i = f_angular / (u_norm * sin_t);
            let scale_k = f_angular / (v_norm * sin_t);

            let fi = [
                scale_i * (v_hat[0] - cos_t * u_hat[0]),
                scale_i * (v_hat[1] - cos_t * u_hat[1]),
                scale_i * (v_hat[2] - cos_t * u_hat[2]),
            ];
            let fk = [
                scale_k * (u_hat[0] - cos_t * v_hat[0]),
                scale_k * (u_hat[1] - cos_t * v_hat[1]),
                scale_k * (u_hat[2] - cos_t * v_hat[2]),
            ];
            let fj = [-fi[0] - fk[0], -fi[1] - fk[1], -fi[2] - fk[2]];

            for d in 0..3 {
                forces[*ai][d] += fi[d];
                forces[*aj][d] += fj[d];
                forces[*ak][d] += fk[d];
            }
        }

        // Urey-Bradley force contribution
        if params.urey_bradley_k > 0.0 && r13 > 1e-15 {
            for d in 0..3 {
                let f_ub_d = -f_ub * dr13[d] / r13;
                forces[*ai][d] += f_ub_d;
                forces[*ak][d] -= f_ub_d;
            }
        }
    }

    (total_energy, forces)
}

/// Compute the dihedral angle for atoms i-j-k-l.
pub fn compute_dihedral(positions: &[[f64; 3]], i: usize, j: usize, k: usize, l: usize) -> f64 {
    let ri = positions[i];
    let rj = positions[j];
    let rk = positions[k];
    let rl = positions[l];

    let b1 = [rj[0] - ri[0], rj[1] - ri[1], rj[2] - ri[2]];
    let b2 = [rk[0] - rj[0], rk[1] - rj[1], rk[2] - rj[2]];
    let b3 = [rl[0] - rk[0], rl[1] - rk[1], rl[2] - rk[2]];

    let n1 = cross3(b1, b2);
    let n2 = cross3(b2, b3);
    let m1 = cross3(n1, b2);

    let n1_norm = norm3(n1);
    let n2_norm = norm3(n2);
    let m1_norm = norm3(m1);

    if n1_norm < 1e-15 || n2_norm < 1e-15 || m1_norm < 1e-15 {
        return 0.0;
    }

    let x = dot3(n1, n2) / (n1_norm * n2_norm);
    let y = dot3(m1, n2) / (m1_norm * n2_norm);
    x.clamp(-1.0, 1.0).acos().copysign(y)
}

// ── Internal helpers ──────────────────────────────────────────────────────────

fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn norm3(a: [f64; 3]) -> f64 {
    (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt()
}

/// Return a sorted pair of owned strings for use as a map key.
fn sorted_pair(a: &str, b: &str) -> (String, String) {
    if a <= b {
        (a.to_string(), b.to_string())
    } else {
        (b.to_string(), a.to_string())
    }
}

// ── CMAP correction ───────────────────────────────────────────────────────────

/// CMAP (Coupled torsion correction MAP) backbone correction.
///
/// The CMAP correction adds a 2D grid-based correction to backbone
/// phi/psi dihedral combinations.  The grid is bilinearly interpolated.
///
/// Reference: MacKerell et al., J. Comput. Chem. 25, 1400 (2004).
#[derive(Debug, Clone)]
pub struct CharmmCmap {
    /// Grid resolution (number of points per axis, typically 24).
    pub n_grid: usize,
    /// Correction energy values on the \[phi\]\[psi\] grid (kcal/mol).
    /// Stored row-major with phi varying along rows.
    pub grid: Vec<f64>,
    /// Grid spacing (radians): 2*pi / n_grid.
    pub d_phi: f64,
}

impl CharmmCmap {
    /// Create a new CMAP correction from a flat grid of values.
    ///
    /// `n_grid` must be > 1.  The grid covers \[-pi, pi\] x \[-pi, pi\].
    pub fn new(n_grid: usize, grid: Vec<f64>) -> Self {
        assert_eq!(grid.len(), n_grid * n_grid, "CMAP grid size mismatch");
        let d_phi = 2.0 * PI / (n_grid as f64);
        Self {
            n_grid,
            grid,
            d_phi,
        }
    }

    /// Create a zero CMAP (no correction).
    pub fn zero(n_grid: usize) -> Self {
        Self::new(n_grid, vec![0.0; n_grid * n_grid])
    }

    /// Bilinear interpolation of the CMAP energy at (phi, psi).
    ///
    /// Both angles are wrapped to \[-pi, pi\] before lookup.
    pub fn energy(&self, phi: f64, psi: f64) -> f64 {
        let _n = self.n_grid as f64;
        let dp = self.d_phi;

        // Map angles to [0, n) indices
        let phi_w = phi.rem_euclid(2.0 * PI);
        let psi_w = psi.rem_euclid(2.0 * PI);

        let i0 = (phi_w / dp) as usize % self.n_grid;
        let j0 = (psi_w / dp) as usize % self.n_grid;
        let i1 = (i0 + 1) % self.n_grid;
        let j1 = (j0 + 1) % self.n_grid;

        let t = (phi_w / dp) - (i0 as f64);
        let u = (psi_w / dp) - (j0 as f64);

        let g = |i: usize, j: usize| self.grid[i * self.n_grid + j];

        (1.0 - t) * (1.0 - u) * g(i0, j0)
            + t * (1.0 - u) * g(i1, j0)
            + (1.0 - t) * u * g(i0, j1)
            + t * u * g(i1, j1)
    }

    /// Gradient (d E/d phi, d E/d psi) via bilinear differentiation.
    pub fn gradient(&self, phi: f64, psi: f64) -> (f64, f64) {
        let dp = self.d_phi;
        let n = self.n_grid;

        let phi_w = phi.rem_euclid(2.0 * PI);
        let psi_w = psi.rem_euclid(2.0 * PI);

        let i0 = (phi_w / dp) as usize % n;
        let j0 = (psi_w / dp) as usize % n;
        let i1 = (i0 + 1) % n;
        let j1 = (j0 + 1) % n;

        let t = (phi_w / dp) - (i0 as f64);
        let u = (psi_w / dp) - (j0 as f64);

        let g = |i: usize, j: usize| self.grid[i * n + j];

        let de_dt = (1.0 - u) * (g(i1, j0) - g(i0, j0)) + u * (g(i1, j1) - g(i0, j1));
        let de_du = (1.0 - t) * (g(i0, j1) - g(i0, j0)) + t * (g(i1, j1) - g(i1, j0));

        (de_dt / dp, de_du / dp)
    }
}

// ── Nonbonded 1-4 scaling ─────────────────────────────────────────────────────

/// CHARMM 1-4 nonbonded scaling factors.
///
/// In CHARMM, atoms separated by exactly 3 bonds (1-4 pair) use reduced
/// vdW and electrostatic interactions to avoid double-counting with
/// the torsion term.
#[derive(Debug, Clone, Copy)]
pub struct Nonbonded14Scaling {
    /// vdW scaling factor for 1-4 pairs (CHARMM: 1.0).
    pub vdw_scale: f64,
    /// Electrostatic scaling factor for 1-4 pairs (CHARMM: 1.0).
    pub elec_scale: f64,
}

impl Nonbonded14Scaling {
    /// CHARMM default: both vdW and electrostatic unscaled (factor = 1.0).
    pub fn charmm_default() -> Self {
        Self {
            vdw_scale: 1.0,
            elec_scale: 1.0,
        }
    }

    /// AMBER convention: vdW scaled by 0.5, electrostatic by 0.8333.
    pub fn amber() -> Self {
        Self {
            vdw_scale: 0.5,
            elec_scale: 1.0 / 1.2,
        }
    }

    /// Apply scaling to a raw vdW energy for a 1-4 pair.
    pub fn scale_vdw(&self, raw_energy: f64) -> f64 {
        self.vdw_scale * raw_energy
    }

    /// Apply scaling to a raw electrostatic energy for a 1-4 pair.
    pub fn scale_elec(&self, raw_energy: f64) -> f64 {
        self.elec_scale * raw_energy
    }
}

/// Compute Coulomb electrostatic energy between two point charges.
///
/// V_elec = k_e * q1 * q2 / r
/// where k_e = 332.0 kcal·Å/(mol·e²) is the CHARMM electrostatic conversion.
pub fn compute_coulomb_energy(q1: f64, q2: f64, r: f64) -> f64 {
    const KE_CHARMM: f64 = 332.0; // kcal·Å/(mol·e²)
    if r < 1e-10 {
        return f64::MAX;
    }
    KE_CHARMM * q1 * q2 / r
}

// ── Water models ──────────────────────────────────────────────────────────────

/// Rigid water model geometry parameters.
///
/// Defines O-H bond length and H-O-H angle for rigid water models.
#[derive(Debug, Clone, Copy)]
pub struct WaterModel {
    /// O-H bond length (Å).
    pub oh_bond: f64,
    /// H-O-H angle (radians).
    pub hoh_angle: f64,
    /// Oxygen partial charge (e).
    pub charge_o: f64,
    /// Hydrogen partial charge (e).
    pub charge_h: f64,
    /// vdW parameters for oxygen.
    pub vdw_o: (f64, f64), // (epsilon, rmin_half)
    /// Optional 4th site offset for TIP4P (Å along bisector, 0 = TIP3P).
    pub m_site_offset: f64,
}

impl WaterModel {
    /// TIP3P water model (Jorgensen et al. 1983).
    pub fn tip3p() -> Self {
        Self {
            oh_bond: 0.9572,
            hoh_angle: 104.52_f64.to_radians(),
            charge_o: -0.834,
            charge_h: 0.417,
            vdw_o: (0.1521, 1.7682),
            m_site_offset: 0.0,
        }
    }

    /// TIP4P water model (Jorgensen et al. 1983) — 4-site.
    pub fn tip4p() -> Self {
        Self {
            oh_bond: 0.9572,
            hoh_angle: 104.52_f64.to_radians(),
            charge_o: 0.0,
            charge_h: 0.52,
            vdw_o: (0.1550, 1.7700),
            m_site_offset: 0.15, // M site offset along H-O-H bisector
        }
    }

    /// TIP4P/2005 water model (Abascal & Vega 2005).
    pub fn tip4p_2005() -> Self {
        Self {
            oh_bond: 0.9572,
            hoh_angle: 104.52_f64.to_radians(),
            charge_o: 0.0,
            charge_h: 0.5564,
            vdw_o: (0.1852, 1.7769),
            m_site_offset: 0.1546,
        }
    }

    /// Water molecule dipole moment (Debye) for this model.
    ///
    /// mu = 2 * q_H * r_OH * sin(theta/2) * (1/0.2082)  → convert to Debye
    pub fn dipole_moment_debye(&self) -> f64 {
        // effective charge is charge_h if tip3p, else use total charge balance
        let q = self.charge_h;
        let r = self.oh_bond; // Å
        let half_angle = self.hoh_angle / 2.0;
        // H-H separation contribution
        let mu_au = 2.0 * q * r * half_angle.sin();
        mu_au / 0.2082 // convert from e·Å to Debye
    }

    /// Check if this is a 4-site model.
    pub fn is_four_site(&self) -> bool {
        self.m_site_offset.abs() > 1e-10
    }
}

// ── CHARMM Energy Decomposition ───────────────────────────────────────────────

/// Full energy decomposition for a CHARMM simulation.
#[derive(Debug, Clone, Default)]
pub struct CharmmEnergyDecomposition {
    /// Bond stretching energy (kcal/mol).
    pub bonds: f64,
    /// Angle bending energy (kcal/mol).
    pub angles: f64,
    /// Urey-Bradley energy (kcal/mol).
    pub urey_bradley: f64,
    /// Proper dihedral energy (kcal/mol).
    pub dihedrals: f64,
    /// Improper dihedral energy (kcal/mol).
    pub impropers: f64,
    /// CMAP correction energy (kcal/mol).
    pub cmap: f64,
    /// Van der Waals energy (kcal/mol).
    pub vdw: f64,
    /// Electrostatic energy (kcal/mol).
    pub electrostatics: f64,
}

impl CharmmEnergyDecomposition {
    /// Total potential energy (sum of all terms).
    pub fn total(&self) -> f64 {
        self.bonds
            + self.angles
            + self.urey_bradley
            + self.dihedrals
            + self.impropers
            + self.cmap
            + self.vdw
            + self.electrostatics
    }

    /// Bonded energy only.
    pub fn bonded(&self) -> f64 {
        self.bonds + self.angles + self.urey_bradley + self.dihedrals + self.impropers + self.cmap
    }

    /// Nonbonded energy only.
    pub fn nonbonded(&self) -> f64 {
        self.vdw + self.electrostatics
    }
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // ── Bond force ────────────────────────────────────────────────────────────

    #[test]
    fn test_bond_energy_at_equilibrium_is_zero() {
        let params = CharmmBondParams::new(250.0, 1.52);
        let (energy, force) = compute_bond_force(1.52, &params);
        assert!(
            energy.abs() < 1e-12,
            "Bond energy at equilibrium must be zero, got {energy}"
        );
        assert!(
            force.abs() < 1e-12,
            "Bond force at equilibrium must be zero, got {force}"
        );
    }

    #[test]
    fn test_bond_energy_increases_with_stretch() {
        let params = CharmmBondParams::new(250.0, 1.52);
        let (e0, _) = compute_bond_force(1.52, &params);
        let (e1, _) = compute_bond_force(1.60, &params);
        assert!(e1 > e0, "Energy must increase when bond is stretched");
    }

    #[test]
    fn test_bond_force_sign_at_stretch() {
        let params = CharmmBondParams::new(250.0, 1.52);
        let (_, force) = compute_bond_force(1.60, &params);
        // stretched: force should be negative (attractive, toward r0)
        assert!(
            force < 0.0,
            "Force on stretched bond should be attractive (negative), got {force}"
        );
    }

    #[test]
    fn test_bond_force_sign_at_compression() {
        let params = CharmmBondParams::new(250.0, 1.52);
        let (_, force) = compute_bond_force(1.40, &params);
        // compressed: force should be positive (repulsive, away from r0)
        assert!(
            force > 0.0,
            "Force on compressed bond should be repulsive (positive), got {force}"
        );
    }

    #[test]
    fn test_bond_energy_symmetry() {
        let params = CharmmBondParams::new(250.0, 1.52);
        let (e_stretch, _) = compute_bond_force(1.62, &params);
        let (e_compress, _) = compute_bond_force(1.42, &params);
        assert!(
            (e_stretch - e_compress).abs() < 1e-10,
            "Bond energy should be symmetric around equilibrium"
        );
    }

    // ── Angle force ───────────────────────────────────────────────────────────

    #[test]
    fn test_angle_energy_at_equilibrium_is_zero() {
        let theta0 = PI * 109.5 / 180.0;
        let params = CharmmAngleParams::simple(50.0, theta0);
        let (energy, force) = compute_angle_force(theta0, &params);
        assert!(
            energy.abs() < 1e-12,
            "Angle energy at equilibrium must be zero, got {energy}"
        );
        assert!(
            force.abs() < 1e-12,
            "Angle force at equilibrium must be zero, got {force}"
        );
    }

    #[test]
    fn test_angle_energy_increases_with_deformation() {
        let theta0 = PI * 109.5 / 180.0;
        let params = CharmmAngleParams::simple(50.0, theta0);
        let (e0, _) = compute_angle_force(theta0, &params);
        let (e1, _) = compute_angle_force(theta0 + 0.1, &params);
        assert!(e1 > e0, "Angle energy must increase with deformation");
    }

    #[test]
    fn test_angle_with_urey_bradley_at_equilibrium() {
        let theta0 = PI * 109.5 / 180.0;
        let r_ub0 = 2.55;
        let params = CharmmAngleParams::new(50.0, theta0, 10.0, r_ub0);
        let (e_ang, _) = compute_angle_force(theta0, &params);
        let (e_ub, f_ub) = compute_urey_bradley(r_ub0, &params);
        assert!(e_ang.abs() < 1e-12, "Angle energy at equilibrium: {e_ang}");
        assert!(e_ub.abs() < 1e-12, "UB energy at equilibrium: {e_ub}");
        assert!(f_ub.abs() < 1e-12, "UB force at equilibrium: {f_ub}");
    }

    // ── Dihedral force ────────────────────────────────────────────────────────

    #[test]
    fn test_dihedral_energy_at_minimum() {
        // For n=1, δ=0: V = k*(1 + cos(φ)) → minimum at φ=π where cos(π)=-1 → V=0
        let params = CharmmDihedralParams::new(1.0, 1, 0.0);
        let (energy, _) = compute_dihedral_force(PI, &params);
        assert!(
            energy.abs() < 1e-12,
            "Dihedral energy at minimum (φ=π, n=1, δ=0) must be zero, got {energy}"
        );
    }

    #[test]
    fn test_dihedral_energy_at_maximum() {
        // For n=1, δ=0: V = k*(1 + cos(φ)) → maximum at φ=0 where cos(0)=1 → V=2k
        let k = 2.5;
        let params = CharmmDihedralParams::new(k, 1, 0.0);
        let (energy, _) = compute_dihedral_force(0.0, &params);
        assert!(
            (energy - 2.0 * k).abs() < 1e-10,
            "Dihedral energy at maximum must be 2k={}, got {energy}",
            2.0 * k
        );
    }

    #[test]
    fn test_dihedral_energy_bounded() {
        let params = CharmmDihedralParams::new(1.0, 3, PI / 3.0);
        for i in 0..360 {
            let phi = i as f64 * PI / 180.0;
            let (e, _) = compute_dihedral_force(phi, &params);
            assert!(
                e >= 0.0,
                "Dihedral energy must be >= 0, got {e} at phi={phi}"
            );
            assert!(
                e <= 2.0 * 1.0 + 1e-10,
                "Dihedral energy must be <= 2k, got {e}"
            );
        }
    }

    #[test]
    fn test_dihedral_force_zero_at_minimum() {
        // For n=2, δ=π: V = k*(1 + cos(2φ - π)) → minimum at φ=π/2
        let params = CharmmDihedralParams::new(1.0, 2, PI);
        let (_, force) = compute_dihedral_force(PI / 2.0, &params);
        assert!(
            force.abs() < 1e-10,
            "Dihedral force at minimum must be zero, got {force}"
        );
    }

    #[test]
    fn test_dihedral_multi_terms_sum() {
        let terms = vec![
            CharmmDihedralParams::new(0.2, 1, 0.0),
            CharmmDihedralParams::new(0.4, 3, 0.0),
        ];
        let phi = 1.0;
        let (e1, f1) = compute_dihedral_force(phi, &terms[0]);
        let (e2, f2) = compute_dihedral_force(phi, &terms[1]);
        let (e_total, f_total) = compute_dihedral_force_multi(phi, &terms);
        assert!(
            (e_total - e1 - e2).abs() < 1e-12,
            "Multi-term energy mismatch: {e_total} vs {}",
            e1 + e2
        );
        assert!(
            (f_total - f1 - f2).abs() < 1e-12,
            "Multi-term force mismatch: {f_total} vs {}",
            f1 + f2
        );
    }

    // ── Improper force ────────────────────────────────────────────────────────

    #[test]
    fn test_improper_energy_at_equilibrium_is_zero() {
        let params = CharmmImproperParams::new(10.5, 0.0);
        let (energy, force) = compute_improper_force(0.0, &params);
        assert!(
            energy.abs() < 1e-12,
            "Improper energy at equilibrium must be zero, got {energy}"
        );
        assert!(
            force.abs() < 1e-12,
            "Improper force at equilibrium must be zero, got {force}"
        );
    }

    #[test]
    fn test_improper_energy_positive_away_from_eq() {
        let params = CharmmImproperParams::new(10.5, 0.0);
        let (e1, _) = compute_improper_force(0.2, &params);
        let (e2, _) = compute_improper_force(-0.2, &params);
        assert!(e1 > 0.0);
        assert!(
            (e1 - e2).abs() < 1e-10,
            "Improper energy should be symmetric"
        );
    }

    // ── vdW force ─────────────────────────────────────────────────────────────

    #[test]
    fn test_vdw_energy_at_rmin_is_minus_epsilon() {
        let epsilon = 0.07;
        let rmin = 2.0 * 2.0; // rmin = 2 * rmin_half
        let (energy, _) = compute_vdw_force(rmin, epsilon, rmin);
        assert!(
            (energy - (-epsilon)).abs() < 1e-10,
            "vdW energy at Rmin must equal -epsilon={epsilon}, got {energy}"
        );
    }

    #[test]
    fn test_vdw_energy_positive_at_short_range() {
        let (energy, _) = compute_vdw_force(1.0, 0.1, 3.0);
        assert!(
            energy > 0.0,
            "vdW energy must be repulsive at short range, got {energy}"
        );
    }

    #[test]
    fn test_vdw_energy_approaches_zero_at_large_r() {
        let (energy, _) = compute_vdw_force(100.0, 0.1, 2.0);
        assert!(
            energy.abs() < 1e-6,
            "vdW energy must approach zero at large r, got {energy}"
        );
    }

    // ── CharmmForceField ──────────────────────────────────────────────────────

    #[test]
    fn test_charmm_ff_add_get_bond() {
        let mut ff = CharmmForceField::new();
        ff.add_bond("CT1", "C", CharmmBondParams::new(250.0, 1.52));
        let p = ff
            .get_bond("C", "CT1")
            .expect("Bond should be found in reversed order");
        assert!((p.k - 250.0).abs() < 1e-10);
        assert!((p.r0 - 1.52).abs() < 1e-10);
    }

    #[test]
    fn test_charmm_ff_add_get_angle() {
        let mut ff = CharmmForceField::new();
        ff.add_angle("NH1", "CT1", "C", CharmmAngleParams::simple(50.0, 1.94));
        let p = ff
            .get_angle("NH1", "CT1", "C")
            .expect("Angle should be found");
        assert!((p.k - 50.0).abs() < 1e-10);
    }

    #[test]
    fn test_charmm_ff_add_get_dihedral() {
        let mut ff = CharmmForceField::new();
        ff.add_dihedral(
            "C",
            "NH1",
            "CT1",
            "C",
            CharmmDihedralParams::new(0.2, 1, 0.0),
        );
        ff.add_dihedral(
            "C",
            "NH1",
            "CT1",
            "C",
            CharmmDihedralParams::new(0.4, 3, 0.0),
        );
        let terms = ff
            .get_dihedral("C", "NH1", "CT1", "C")
            .expect("Dihedral should be found");
        assert_eq!(terms.len(), 2, "Should have 2 dihedral terms");
    }

    #[test]
    fn test_charmm_ff_combined_vdw() {
        let mut ff = CharmmForceField::new();
        ff.add_vdw("CT1", CharmmVdwParams::new(0.02, 2.275));
        ff.add_vdw("NH1", CharmmVdwParams::new(0.20, 1.85));
        let (eps, rmin) = ff.combined_vdw("CT1", "NH1").expect("Combined vdW params");
        let expected_eps = (0.02_f64 * 0.20).sqrt();
        let expected_rmin = 2.275 + 1.85;
        assert!(
            (eps - expected_eps).abs() < 1e-10,
            "Combined epsilon mismatch"
        );
        assert!(
            (rmin - expected_rmin).abs() < 1e-10,
            "Combined Rmin mismatch"
        );
    }

    #[test]
    fn test_charmm22_backbone_has_expected_params() {
        let ff = CharmmForceField::charmm22_backbone();
        // Should contain backbone bond
        assert!(ff.get_bond("CT1", "C").is_some(), "CT1-C bond should exist");
        assert!(
            ff.get_angle("NH1", "CT1", "C").is_some(),
            "NH1-CT1-C angle should exist"
        );
        assert!(ff.get_vdw("C").is_some(), "C vdW params should exist");
    }

    #[test]
    fn test_compute_bonded_forces_collinear_bond() {
        // Two atoms along x-axis, stretched from r0=1.0 to r=2.0
        let positions = vec![[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let bonds = vec![(0usize, 1usize, CharmmBondParams::new(100.0, 1.0))];
        let (energy, forces) = compute_bonded_forces(&positions, &bonds);

        // energy = k*(2-1)^2 = 100
        assert!(
            (energy - 100.0).abs() < 1e-8,
            "Bond energy should be 100, got {energy}"
        );

        // Force on atom 0 should be positive x (toward atom 1, attractive)
        assert!(
            forces[0][0] > 0.0,
            "Force on atom 0 should pull it toward atom 1"
        );
        // Newton's 3rd law
        assert!(
            (forces[0][0] + forces[1][0]).abs() < 1e-10,
            "Net force must be zero"
        );
    }

    #[test]
    fn test_compute_angle_forces_linear_geometry() {
        // Atoms at [-1,0,0], [0,0,0], [1,0,0] → angle = 180°
        let positions = vec![[-1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
        let theta = compute_angle(&positions, 0, 1, 2);
        assert!(
            (theta - PI).abs() < 1e-10,
            "Angle should be PI, got {theta}"
        );

        let params = CharmmAngleParams::simple(50.0, PI);
        let angle_list = vec![(0usize, 1usize, 2usize, params)];
        let (energy, _forces) = compute_angle_forces(&positions, &angle_list);
        assert!(
            energy.abs() < 1e-10,
            "Energy at equilibrium should be zero, got {energy}"
        );
    }

    // ── CMAP tests ─────────────────────────────────────────────────────────────

    #[test]
    fn test_cmap_zero_energy() {
        let cmap = CharmmCmap::zero(6);
        let e = cmap.energy(0.5, 1.2);
        assert!(
            e.abs() < 1e-12,
            "zero CMAP should give zero energy, got {e}"
        );
    }

    #[test]
    fn test_cmap_energy_bilinear_corners() {
        // A 2x2 grid with value 1.0 everywhere → energy should be 1.0 everywhere
        let cmap = CharmmCmap::new(4, vec![1.0; 16]);
        let e = cmap.energy(0.1, 0.2);
        assert!(
            (e - 1.0).abs() < 1e-10,
            "uniform grid should give 1.0, got {e}"
        );
    }

    #[test]
    fn test_cmap_energy_finite() {
        let n = 8;
        let mut grid = vec![0.0; n * n];
        // Set a non-trivial grid
        for i in 0..n {
            for j in 0..n {
                grid[i * n + j] = (i as f64) * 0.1 + (j as f64) * 0.05;
            }
        }
        let cmap = CharmmCmap::new(n, grid);
        let e = cmap.energy(1.0, 2.0);
        assert!(e.is_finite(), "CMAP energy should be finite, got {e}");
    }

    #[test]
    fn test_cmap_gradient_zero_for_uniform_grid() {
        let cmap = CharmmCmap::new(4, vec![1.0; 16]);
        let (de_dphi, de_dpsi) = cmap.gradient(0.5, 0.5);
        assert!(
            de_dphi.abs() < 1e-10,
            "uniform grid gradient wrt phi = 0, got {de_dphi}"
        );
        assert!(
            de_dpsi.abs() < 1e-10,
            "uniform grid gradient wrt psi = 0, got {de_dpsi}"
        );
    }

    // ── Nonbonded 1-4 scaling tests ───────────────────────────────────────────

    #[test]
    fn test_nonbonded14_charmm_default_unity() {
        let s = Nonbonded14Scaling::charmm_default();
        let e = 10.0;
        assert!((s.scale_vdw(e) - e).abs() < 1e-12);
        assert!((s.scale_elec(e) - e).abs() < 1e-12);
    }

    #[test]
    fn test_nonbonded14_amber_vdw() {
        let s = Nonbonded14Scaling::amber();
        let e = 10.0;
        assert!(
            (s.scale_vdw(e) - 5.0).abs() < 1e-10,
            "AMBER vdW scale = 0.5"
        );
    }

    #[test]
    fn test_coulomb_energy_positive_same_sign() {
        let e = compute_coulomb_energy(1.0, 1.0, 3.0);
        assert!(e > 0.0, "same-sign charges should repel, got {e}");
    }

    #[test]
    fn test_coulomb_energy_negative_opposite_sign() {
        let e = compute_coulomb_energy(1.0, -1.0, 3.0);
        assert!(e < 0.0, "opposite charges should attract, got {e}");
    }

    #[test]
    fn test_coulomb_energy_inverse_r() {
        let e1 = compute_coulomb_energy(1.0, 1.0, 2.0);
        let e2 = compute_coulomb_energy(1.0, 1.0, 4.0);
        assert!((e1 / e2 - 2.0).abs() < 1e-10, "Coulomb scales as 1/r");
    }

    // ── Water model tests ─────────────────────────────────────────────────────

    #[test]
    fn test_tip3p_bond_length() {
        let w = WaterModel::tip3p();
        assert!((w.oh_bond - 0.9572).abs() < 1e-6);
    }

    #[test]
    fn test_tip3p_charge_neutral() {
        let w = WaterModel::tip3p();
        let total_charge = w.charge_o + 2.0 * w.charge_h;
        assert!(
            total_charge.abs() < 1e-10,
            "water should be neutral, got {total_charge}"
        );
    }

    #[test]
    fn test_tip4p_four_site() {
        let w = WaterModel::tip4p();
        assert!(w.is_four_site(), "TIP4P should be a 4-site model");
    }

    #[test]
    fn test_tip3p_three_site() {
        let w = WaterModel::tip3p();
        assert!(!w.is_four_site(), "TIP3P should NOT be a 4-site model");
    }

    #[test]
    fn test_water_dipole_moment_reasonable() {
        let w = WaterModel::tip3p();
        let mu = w.dipole_moment_debye();
        // TIP3P geometric dipole from q_H * r_OH * sin(θ/2) ≈ 2.3–3.1 D
        assert!(
            mu > 1.5 && mu < 4.0,
            "TIP3P dipole should be ~2-3 D, got {mu}"
        );
    }

    // ── Energy decomposition tests ────────────────────────────────────────────

    #[test]
    fn test_energy_decomposition_total() {
        let e = CharmmEnergyDecomposition {
            bonds: 1.0,
            angles: 2.0,
            urey_bradley: 0.5,
            dihedrals: 3.0,
            impropers: 0.1,
            cmap: 0.2,
            vdw: 5.0,
            electrostatics: -4.0,
        };
        let total = e.total();
        assert!((total - 7.8).abs() < 1e-10, "total = {total}");
    }

    #[test]
    fn test_energy_decomposition_bonded_nonbonded() {
        let e = CharmmEnergyDecomposition {
            bonds: 1.0,
            vdw: 2.0,
            electrostatics: 3.0,
            ..Default::default()
        };
        assert!((e.bonded() - 1.0).abs() < 1e-12);
        assert!((e.nonbonded() - 5.0).abs() < 1e-12);
    }

    #[test]
    fn test_energy_decomposition_zero() {
        let e = CharmmEnergyDecomposition::default();
        assert!(
            e.total().abs() < 1e-12,
            "all-zero decomposition should sum to 0"
        );
    }
}
