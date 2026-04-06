// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Force field construction and validation utilities.
//!
//! Provides parameter containers for atoms, bonds, angles, and dihedrals,
//! together with Lorentz-Berthelot combining rules, force-field comparison,
//! and validation routines.

#![allow(dead_code)]

use std::collections::HashMap;

// ---------------------------------------------------------------------------
// Free functions — combining rules
// ---------------------------------------------------------------------------

/// Geometric combining rule: √(x₁ · x₂).
pub fn combining_rule_geometric(x1: f64, x2: f64) -> f64 {
    (x1 * x2).sqrt()
}

/// Arithmetic combining rule: (x₁ + x₂) / 2.
pub fn combining_rule_arithmetic(x1: f64, x2: f64) -> f64 {
    (x1 + x2) / 2.0
}

/// Reduced mass μ = m₁ m₂ / (m₁ + m₂).
pub fn reduced_mass(m1: f64, m2: f64) -> f64 {
    m1 * m2 / (m1 + m2)
}

// ---------------------------------------------------------------------------
// AtomType
// ---------------------------------------------------------------------------

/// Atom type definition used in a force field.
#[derive(Debug, Clone)]
pub struct AtomType {
    /// IUPAC-style name, e.g. `"CT"`, `"OH"`.
    pub name: String,
    /// Atomic mass (g/mol or u).
    pub mass: f64,
    /// Lennard-Jones σ parameter (Å or nm).
    pub sigma: f64,
    /// Lennard-Jones ε parameter (kcal/mol or kJ/mol).
    pub epsilon: f64,
    /// Partial charge (elementary charges).
    pub charge: f64,
    /// Atomic number (proton count).
    pub atom_number: usize,
}

impl AtomType {
    /// Create a new atom type.
    pub fn new(
        name: &str,
        mass: f64,
        sigma: f64,
        epsilon: f64,
        charge: f64,
        atom_number: usize,
    ) -> Self {
        Self {
            name: name.to_string(),
            mass,
            sigma,
            epsilon,
            charge,
            atom_number,
        }
    }
}

// ---------------------------------------------------------------------------
// BondType
// ---------------------------------------------------------------------------

/// Harmonic (or optionally Morse) bond type.
#[derive(Debug, Clone)]
pub struct BondType {
    /// Name of first atom type.
    pub atom1_type: String,
    /// Name of second atom type.
    pub atom2_type: String,
    /// Equilibrium bond length r₀ (Å or nm).
    pub r0: f64,
    /// Harmonic force constant k_bond (kcal mol⁻¹ Å⁻² or kJ mol⁻¹ nm⁻²).
    pub k_bond: f64,
    /// Optional Morse dissociation energy D_e.
    pub morse_d: Option<f64>,
}

impl BondType {
    /// Create a harmonic bond type (no Morse term).
    pub fn new_harmonic(atom1_type: &str, atom2_type: &str, r0: f64, k_bond: f64) -> Self {
        Self {
            atom1_type: atom1_type.to_string(),
            atom2_type: atom2_type.to_string(),
            r0,
            k_bond,
            morse_d: None,
        }
    }

    /// Create a Morse bond type.
    pub fn new_morse(
        atom1_type: &str,
        atom2_type: &str,
        r0: f64,
        k_bond: f64,
        morse_d: f64,
    ) -> Self {
        Self {
            atom1_type: atom1_type.to_string(),
            atom2_type: atom2_type.to_string(),
            r0,
            k_bond,
            morse_d: Some(morse_d),
        }
    }

    /// Bond key in canonical order `"A-B"` where A ≤ B lexicographically.
    pub fn canonical_key(&self) -> String {
        if self.atom1_type <= self.atom2_type {
            format!("{}-{}", self.atom1_type, self.atom2_type)
        } else {
            format!("{}-{}", self.atom2_type, self.atom1_type)
        }
    }
}

// ---------------------------------------------------------------------------
// AngleType
// ---------------------------------------------------------------------------

/// Harmonic angle type for a three-atom sequence.
#[derive(Debug, Clone)]
pub struct AngleType {
    /// Atom type names for the angle: \[i, j, k\].
    pub atom_types: [String; 3],
    /// Equilibrium angle θ₀ (radians).
    pub theta0: f64,
    /// Harmonic force constant k_angle.
    pub k_angle: f64,
}

impl AngleType {
    /// Create a new harmonic angle type.
    pub fn new(atom_types: [&str; 3], theta0: f64, k_angle: f64) -> Self {
        Self {
            atom_types: atom_types.map(|s| s.to_string()),
            theta0,
            k_angle,
        }
    }

    /// Angle energy U = k_angle / 2 · (θ − θ₀)².
    pub fn energy(&self, theta: f64) -> f64 {
        let dt = theta - self.theta0;
        0.5 * self.k_angle * dt * dt
    }
}

// ---------------------------------------------------------------------------
// DihedralType
// ---------------------------------------------------------------------------

/// OPLS/CHARMM-style dihedral (torsion) with four Fourier coefficients.
#[derive(Debug, Clone)]
pub struct DihedralType {
    /// Atom type names for the dihedral: \[i, j, k, l\].
    pub atom_types: [String; 4],
    /// First Fourier coefficient V1 (kcal/mol or kJ/mol).
    pub v1: f64,
    /// Second Fourier coefficient V2.
    pub v2: f64,
    /// Third Fourier coefficient V3.
    pub v3: f64,
    /// Fourth Fourier coefficient V4.
    pub v4: f64,
}

impl DihedralType {
    /// Create a new dihedral type.
    pub fn new(atom_types: [&str; 4], v1: f64, v2: f64, v3: f64, v4: f64) -> Self {
        Self {
            atom_types: atom_types.map(|s| s.to_string()),
            v1,
            v2,
            v3,
            v4,
        }
    }

    /// OPLS Fourier dihedral energy:
    /// U(φ) = V1/2·(1+cosφ) + V2/2·(1−cos2φ) + V3/2·(1+cos3φ) + V4/2·(1−cos4φ)
    pub fn energy(&self, phi: f64) -> f64 {
        0.5 * self.v1 * (1.0 + phi.cos())
            + 0.5 * self.v2 * (1.0 - (2.0 * phi).cos())
            + 0.5 * self.v3 * (1.0 + (3.0 * phi).cos())
            + 0.5 * self.v4 * (1.0 - (4.0 * phi).cos())
    }
}

// ---------------------------------------------------------------------------
// ForceFieldSpec
// ---------------------------------------------------------------------------

/// Complete force field specification: atom types plus bonded parameters.
#[derive(Debug, Clone, Default)]
pub struct ForceFieldSpec {
    /// Atom types indexed by name.
    pub atom_types: HashMap<String, AtomType>,
    /// Bond types indexed by canonical key.
    pub bond_types: HashMap<String, BondType>,
    /// Angle types stored as a list.
    pub angle_types: Vec<AngleType>,
    /// Dihedral types stored as a list.
    pub dihedral_types: Vec<DihedralType>,
}

impl ForceFieldSpec {
    /// Create an empty force field specification.
    pub fn new() -> Self {
        Self::default()
    }

    /// Add an atom type.  Returns `false` if the name was already present.
    pub fn add_atom(&mut self, atom: AtomType) -> bool {
        if self.atom_types.contains_key(&atom.name) {
            return false;
        }
        self.atom_types.insert(atom.name.clone(), atom);
        true
    }

    /// Add a bond type by canonical key.  Returns `false` if already present.
    pub fn add_bond(&mut self, bond: BondType) -> bool {
        let key = bond.canonical_key();
        if self.bond_types.contains_key(&key) {
            return false;
        }
        self.bond_types.insert(key, bond);
        true
    }

    /// Add an angle type.
    pub fn add_angle(&mut self, angle: AngleType) {
        self.angle_types.push(angle);
    }

    /// Add a dihedral type.
    pub fn add_dihedral(&mut self, dihedral: DihedralType) {
        self.dihedral_types.push(dihedral);
    }

    /// Validate the force field and return a list of error messages.
    ///
    /// Checks performed:
    /// 1. All atom types have positive mass.
    /// 2. All atom types have non-negative σ and ε.
    /// 3. Bond r₀ and k_bond are positive.
    /// 4. Angle θ₀ is in (0, π) and k_angle ≥ 0.
    /// 5. Bond atom types exist in the atom type table.
    pub fn validate(&self) -> Vec<String> {
        let mut errors = Vec::new();

        for (name, at) in &self.atom_types {
            if at.mass <= 0.0 {
                errors.push(format!(
                    "AtomType '{name}': mass must be positive (got {})",
                    at.mass
                ));
            }
            if at.sigma < 0.0 {
                errors.push(format!(
                    "AtomType '{name}': sigma must be ≥ 0 (got {})",
                    at.sigma
                ));
            }
            if at.epsilon < 0.0 {
                errors.push(format!(
                    "AtomType '{name}': epsilon must be ≥ 0 (got {})",
                    at.epsilon
                ));
            }
        }

        for (key, bt) in &self.bond_types {
            if bt.r0 <= 0.0 {
                errors.push(format!(
                    "BondType '{key}': r0 must be positive (got {})",
                    bt.r0
                ));
            }
            if bt.k_bond <= 0.0 {
                errors.push(format!(
                    "BondType '{key}': k_bond must be positive (got {})",
                    bt.k_bond
                ));
            }
            if !self.atom_types.contains_key(&bt.atom1_type) {
                errors.push(format!(
                    "BondType '{key}': atom type '{}' not found",
                    bt.atom1_type
                ));
            }
            if !self.atom_types.contains_key(&bt.atom2_type) {
                errors.push(format!(
                    "BondType '{key}': atom type '{}' not found",
                    bt.atom2_type
                ));
            }
        }

        for (i, ang) in self.angle_types.iter().enumerate() {
            if ang.theta0 <= 0.0 || ang.theta0 >= std::f64::consts::PI {
                errors.push(format!(
                    "AngleType[{i}]: theta0 must be in (0, π) (got {})",
                    ang.theta0
                ));
            }
            if ang.k_angle < 0.0 {
                errors.push(format!(
                    "AngleType[{i}]: k_angle must be ≥ 0 (got {})",
                    ang.k_angle
                ));
            }
        }

        errors
    }
}

// ---------------------------------------------------------------------------
// LorentzBerthelotMixer
// ---------------------------------------------------------------------------

/// Lorentz-Berthelot combining rules for unlike LJ pairs.
#[derive(Debug, Clone, Default)]
pub struct LorentzBerthelotMixer;

impl LorentzBerthelotMixer {
    /// Create a new mixer instance.
    pub fn new() -> Self {
        Self
    }

    /// Cross σ_ij = (σ_i + σ_j) / 2  (arithmetic mean).
    pub fn compute_cross_sigma(&self, s1: f64, s2: f64) -> f64 {
        combining_rule_arithmetic(s1, s2)
    }

    /// Cross ε_ij = √(ε_i · ε_j)  (geometric mean).
    pub fn compute_cross_epsilon(&self, e1: f64, e2: f64) -> f64 {
        combining_rule_geometric(e1, e2)
    }

    /// Generate combined atom type parameters for an unlike pair.
    ///
    /// Returns `(sigma_cross, epsilon_cross)`.
    pub fn mix(&self, at1: &AtomType, at2: &AtomType) -> (f64, f64) {
        (
            self.compute_cross_sigma(at1.sigma, at2.sigma),
            self.compute_cross_epsilon(at1.epsilon, at2.epsilon),
        )
    }
}

// ---------------------------------------------------------------------------
// ForceFieldComparison
// ---------------------------------------------------------------------------

/// Compare the LJ energy surfaces of two force fields over a set of pair distances.
#[derive(Debug, Clone, Default)]
pub struct ForceFieldComparison;

impl ForceFieldComparison {
    /// Create a new comparison instance.
    pub fn new() -> Self {
        Self
    }

    /// Compute the mean absolute difference in LJ pair energies between two
    /// force fields evaluated at the positions given in `geometry`.
    ///
    /// Each row in `geometry` is `[x, y, z]` (unused here — only pair
    /// distances between consecutive atoms are evaluated).  Returns the mean
    /// absolute energy difference over all evaluated pairs.
    pub fn compare_energies(
        &self,
        ff1: &ForceFieldSpec,
        ff2: &ForceFieldSpec,
        geometry: &[[f64; 3]],
    ) -> f64 {
        if geometry.len() < 2 {
            return 0.0;
        }

        let mixer1 = LorentzBerthelotMixer::new();
        let mixer2 = LorentzBerthelotMixer::new();

        let atoms1: Vec<&AtomType> = ff1.atom_types.values().collect();
        let atoms2: Vec<&AtomType> = ff2.atom_types.values().collect();

        if atoms1.is_empty() || atoms2.is_empty() {
            return 0.0;
        }

        let n = geometry.len();
        let mut total = 0.0;
        let mut count = 0usize;

        for i in 0..n - 1 {
            let dx = geometry[i + 1][0] - geometry[i][0];
            let dy = geometry[i + 1][1] - geometry[i][1];
            let dz = geometry[i + 1][2] - geometry[i][2];
            let r = (dx * dx + dy * dy + dz * dz).sqrt();
            if r < 1e-10 {
                continue;
            }

            let at1a = atoms1[i % atoms1.len()];
            let at1b = atoms1[(i + 1) % atoms1.len()];
            let at2a = atoms2[i % atoms2.len()];
            let at2b = atoms2[(i + 1) % atoms2.len()];

            let (s1, e1) = mixer1.mix(at1a, at1b);
            let (s2, e2) = mixer2.mix(at2a, at2b);

            let e_lj1 = lj_energy(r, e1, s1);
            let e_lj2 = lj_energy(r, e2, s2);
            total += (e_lj1 - e_lj2).abs();
            count += 1;
        }

        if count == 0 {
            0.0
        } else {
            total / count as f64
        }
    }
}

// ---------------------------------------------------------------------------
// Internal helpers
// ---------------------------------------------------------------------------

/// LJ 12-6 energy: 4 ε \[(σ/r)¹² − (σ/r)⁶\].
fn lj_energy(r: f64, epsilon: f64, sigma: f64) -> f64 {
    let sr = sigma / r;
    let sr6 = sr.powi(6);
    4.0 * epsilon * (sr6 * sr6 - sr6)
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    // ── combining rules ─────────────────────────────────────────────────────

    #[test]
    fn test_geometric_equal() {
        assert!((combining_rule_geometric(4.0, 4.0) - 4.0).abs() < 1e-14);
    }

    #[test]
    fn test_geometric_distinct() {
        let result = combining_rule_geometric(9.0, 4.0);
        assert!((result - 6.0).abs() < 1e-12);
    }

    #[test]
    fn test_geometric_zero() {
        assert_eq!(combining_rule_geometric(0.0, 5.0), 0.0);
    }

    #[test]
    fn test_arithmetic_equal() {
        assert!((combining_rule_arithmetic(3.0, 3.0) - 3.0).abs() < 1e-14);
    }

    #[test]
    fn test_arithmetic_distinct() {
        assert!((combining_rule_arithmetic(1.0, 3.0) - 2.0).abs() < 1e-14);
    }

    #[test]
    fn test_arithmetic_zero() {
        assert_eq!(combining_rule_arithmetic(0.0, 0.0), 0.0);
    }

    #[test]
    fn test_reduced_mass_equal() {
        // μ(m, m) = m/2
        assert!((reduced_mass(12.0, 12.0) - 6.0).abs() < 1e-12);
    }

    #[test]
    fn test_reduced_mass_hydrogen_carbon() {
        // C-H: m1=12, m2=1  → μ = 12/13
        let mu = reduced_mass(12.0, 1.0);
        let expected = 12.0 / 13.0;
        assert!((mu - expected).abs() < 1e-12);
    }

    #[test]
    fn test_reduced_mass_symmetry() {
        assert!((reduced_mass(14.0, 16.0) - reduced_mass(16.0, 14.0)).abs() < 1e-14);
    }

    // ── AtomType ─────────────────────────────────────────────────────────────

    #[test]
    fn test_atom_type_new() {
        let at = AtomType::new("CT", 12.011, 3.5, 0.066, -0.18, 6);
        assert_eq!(at.name, "CT");
        assert_eq!(at.atom_number, 6);
    }

    #[test]
    fn test_atom_type_fields() {
        let at = AtomType::new("OH", 15.999, 3.12, 0.21, -0.64, 8);
        assert!((at.charge - (-0.64)).abs() < 1e-14);
        assert_eq!(at.atom_number, 8);
    }

    // ── BondType ─────────────────────────────────────────────────────────────

    #[test]
    fn test_bond_harmonic_new() {
        let bt = BondType::new_harmonic("CT", "OH", 1.41, 450.0);
        assert!(bt.morse_d.is_none());
    }

    #[test]
    fn test_bond_morse_new() {
        let bt = BondType::new_morse("CT", "OH", 1.41, 450.0, 100.0);
        assert_eq!(bt.morse_d, Some(100.0));
    }

    #[test]
    fn test_bond_canonical_key_order() {
        let bt = BondType::new_harmonic("OH", "CT", 1.41, 450.0);
        assert_eq!(bt.canonical_key(), "CT-OH");
    }

    #[test]
    fn test_bond_canonical_key_same() {
        let bt = BondType::new_harmonic("CT", "CT", 1.52, 400.0);
        assert_eq!(bt.canonical_key(), "CT-CT");
    }

    // ── AngleType ─────────────────────────────────────────────────────────────

    #[test]
    fn test_angle_energy_at_equilibrium() {
        let at = AngleType::new(["CT", "CT", "OH"], PI / 2.0, 100.0);
        assert!((at.energy(PI / 2.0)).abs() < 1e-14);
    }

    #[test]
    fn test_angle_energy_displaced() {
        let at = AngleType::new(["CT", "CT", "OH"], PI / 2.0, 100.0);
        let dtheta = 0.1_f64;
        let expected = 0.5 * 100.0 * dtheta * dtheta;
        assert!((at.energy(PI / 2.0 + dtheta) - expected).abs() < 1e-12);
    }

    // ── DihedralType ─────────────────────────────────────────────────────────

    #[test]
    fn test_dihedral_energy_zero_phi() {
        // V1=1, V2=V3=V4=0, phi=0 → U = V1/2*(1+1) = 1
        let dt = DihedralType::new(["C", "C", "C", "C"], 1.0, 0.0, 0.0, 0.0);
        assert!((dt.energy(0.0) - 1.0).abs() < 1e-14);
    }

    #[test]
    fn test_dihedral_energy_pi() {
        // V1=1, phi=π → U = V1/2*(1+cos(π)) = 0
        let dt = DihedralType::new(["C", "C", "C", "C"], 1.0, 0.0, 0.0, 0.0);
        assert!(dt.energy(PI).abs() < 1e-14);
    }

    #[test]
    fn test_dihedral_v2_barrier() {
        // V2 term: 0.5*V2*(1 - cos(2φ)), at φ = π/2: cos(π) = -1 → 0.5*2*(1-(-1)) = 2
        let dt = DihedralType::new(["C", "C", "C", "C"], 0.0, 2.0, 0.0, 0.0);
        let u = dt.energy(PI / 2.0);
        assert!((u - 2.0).abs() < 1e-12, "u = {u}");
    }

    // ── ForceFieldSpec ────────────────────────────────────────────────────────

    #[test]
    fn test_force_field_add_atom() {
        let mut ff = ForceFieldSpec::new();
        let at = AtomType::new("CT", 12.011, 3.5, 0.066, -0.18, 6);
        assert!(ff.add_atom(at));
    }

    #[test]
    fn test_force_field_add_atom_duplicate() {
        let mut ff = ForceFieldSpec::new();
        ff.add_atom(AtomType::new("CT", 12.011, 3.5, 0.066, -0.18, 6));
        let result = ff.add_atom(AtomType::new("CT", 12.011, 3.5, 0.066, -0.18, 6));
        assert!(!result);
    }

    #[test]
    fn test_force_field_add_bond() {
        let mut ff = ForceFieldSpec::new();
        ff.add_atom(AtomType::new("CT", 12.011, 3.5, 0.066, -0.18, 6));
        ff.add_atom(AtomType::new("OH", 15.999, 3.12, 0.21, -0.64, 8));
        let bt = BondType::new_harmonic("CT", "OH", 1.41, 450.0);
        assert!(ff.add_bond(bt));
    }

    #[test]
    fn test_force_field_add_bond_duplicate() {
        let mut ff = ForceFieldSpec::new();
        ff.add_bond(BondType::new_harmonic("CT", "OH", 1.41, 450.0));
        let result = ff.add_bond(BondType::new_harmonic("OH", "CT", 1.41, 450.0));
        assert!(!result);
    }

    #[test]
    fn test_validate_empty_ff() {
        let ff = ForceFieldSpec::new();
        assert!(ff.validate().is_empty());
    }

    #[test]
    fn test_validate_negative_mass() {
        let mut ff = ForceFieldSpec::new();
        ff.add_atom(AtomType::new("X", -1.0, 3.5, 0.1, 0.0, 1));
        let errs = ff.validate();
        assert!(!errs.is_empty());
        assert!(errs[0].contains("mass"));
    }

    #[test]
    fn test_validate_missing_atom_in_bond() {
        let mut ff = ForceFieldSpec::new();
        ff.add_atom(AtomType::new("CT", 12.011, 3.5, 0.066, 0.0, 6));
        ff.bond_types.insert(
            "CT-OH".to_string(),
            BondType::new_harmonic("CT", "OH", 1.41, 450.0),
        );
        let errs = ff.validate();
        assert!(errs.iter().any(|e| e.contains("OH")));
    }

    #[test]
    fn test_validate_bad_r0() {
        let mut ff = ForceFieldSpec::new();
        ff.add_atom(AtomType::new("CT", 12.011, 3.5, 0.066, 0.0, 6));
        ff.bond_types.insert(
            "CT-CT".to_string(),
            BondType {
                atom1_type: "CT".to_string(),
                atom2_type: "CT".to_string(),
                r0: -1.0,
                k_bond: 400.0,
                morse_d: None,
            },
        );
        let errs = ff.validate();
        assert!(errs.iter().any(|e| e.contains("r0")));
    }

    #[test]
    fn test_validate_angle_out_of_range() {
        let mut ff = ForceFieldSpec::new();
        ff.add_angle(AngleType::new(["CT", "CT", "OH"], 0.0, 100.0));
        let errs = ff.validate();
        assert!(errs.iter().any(|e| e.contains("theta0")));
    }

    // ── LorentzBerthelotMixer ─────────────────────────────────────────────────

    #[test]
    fn test_lb_mixer_sigma() {
        let mixer = LorentzBerthelotMixer::new();
        assert!((mixer.compute_cross_sigma(3.0, 5.0) - 4.0).abs() < 1e-14);
    }

    #[test]
    fn test_lb_mixer_epsilon() {
        let mixer = LorentzBerthelotMixer::new();
        assert!((mixer.compute_cross_epsilon(4.0, 9.0) - 6.0).abs() < 1e-12);
    }

    #[test]
    fn test_lb_mixer_mix() {
        let mixer = LorentzBerthelotMixer::new();
        let at1 = AtomType::new("CT", 12.0, 3.0, 4.0, 0.0, 6);
        let at2 = AtomType::new("OH", 16.0, 5.0, 9.0, 0.0, 8);
        let (sigma, epsilon) = mixer.mix(&at1, &at2);
        assert!((sigma - 4.0).abs() < 1e-14);
        assert!((epsilon - 6.0).abs() < 1e-12);
    }

    #[test]
    fn test_lb_mixer_same_type() {
        let mixer = LorentzBerthelotMixer::new();
        let at = AtomType::new("CT", 12.0, 3.5, 0.066, 0.0, 6);
        let (sigma, epsilon) = mixer.mix(&at, &at.clone());
        assert!((sigma - 3.5).abs() < 1e-14);
        assert!((epsilon - 0.066).abs() < 1e-14);
    }

    // ── ForceFieldComparison ──────────────────────────────────────────────────

    #[test]
    fn test_compare_energies_identical_ffs() {
        let mut ff = ForceFieldSpec::new();
        ff.add_atom(AtomType::new("CT", 12.0, 3.5, 0.066, 0.0, 6));
        let geometry = [[0.0, 0.0, 0.0], [3.5, 0.0, 0.0], [7.0, 0.0, 0.0]];
        let cmp = ForceFieldComparison::new();
        let diff = cmp.compare_energies(&ff, &ff.clone(), &geometry);
        assert!((diff).abs() < 1e-10);
    }

    #[test]
    fn test_compare_energies_empty_geometry() {
        let ff = ForceFieldSpec::new();
        let cmp = ForceFieldComparison::new();
        assert_eq!(cmp.compare_energies(&ff, &ff, &[]), 0.0);
    }

    #[test]
    fn test_compare_energies_single_point() {
        let ff = ForceFieldSpec::new();
        let cmp = ForceFieldComparison::new();
        // Only 1 point → no pairs → 0
        assert_eq!(cmp.compare_energies(&ff, &ff, &[[0.0, 0.0, 0.0]]), 0.0);
    }

    #[test]
    fn test_compare_energies_different_epsilon() {
        let mut ff1 = ForceFieldSpec::new();
        ff1.add_atom(AtomType::new("CT", 12.0, 3.5, 0.066, 0.0, 6));
        let mut ff2 = ForceFieldSpec::new();
        ff2.add_atom(AtomType::new("CT", 12.0, 3.5, 1.0, 0.0, 6));
        let geometry = [[0.0, 0.0, 0.0], [4.0, 0.0, 0.0]];
        let cmp = ForceFieldComparison::new();
        let diff = cmp.compare_energies(&ff1, &ff2, &geometry);
        // Different ε → non-zero difference
        assert!(diff > 0.0);
    }

    // ── Integration tests ─────────────────────────────────────────────────────

    #[test]
    fn test_full_ff_round_trip() {
        let mut ff = ForceFieldSpec::new();
        ff.add_atom(AtomType::new("CT", 12.011, 3.5, 0.066, -0.18, 6));
        ff.add_atom(AtomType::new("HC", 1.008, 2.5, 0.030, 0.06, 1));
        ff.add_bond(BondType::new_harmonic("CT", "HC", 1.09, 340.0));
        ff.add_angle(AngleType::new(["HC", "CT", "HC"], 1.9111, 33.0));
        ff.add_dihedral(DihedralType::new(
            ["HC", "CT", "CT", "HC"],
            0.0,
            0.0,
            0.3,
            0.0,
        ));
        let errors = ff.validate();
        assert!(errors.is_empty(), "Errors: {errors:?}");
        assert_eq!(ff.atom_types.len(), 2);
        assert_eq!(ff.bond_types.len(), 1);
    }

    #[test]
    fn test_reduced_mass_proton() {
        // Proton (1 Da) homodiatomic: μ = 0.5 Da
        assert!((reduced_mass(1.0, 1.0) - 0.5).abs() < 1e-14);
    }

    #[test]
    fn test_dihedral_symmetry() {
        let dt = DihedralType::new(["C", "C", "C", "C"], 1.0, 1.0, 1.0, 1.0);
        // OPLS is symmetric around φ = 0: U(φ) = U(-φ)
        assert!((dt.energy(0.5) - dt.energy(-0.5)).abs() < 1e-12);
    }

    #[test]
    fn test_angle_type_new_fields() {
        let at = AngleType::new(["CA", "CB", "CC"], 2.094, 80.0);
        assert_eq!(at.atom_types[0], "CA");
        assert_eq!(at.atom_types[1], "CB");
        assert!((at.theta0 - 2.094).abs() < 1e-14);
        assert!((at.k_angle - 80.0).abs() < 1e-14);
    }
}
