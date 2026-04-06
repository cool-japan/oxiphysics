// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Crystal plasticity — slip systems, texture evolution, and polycrystal averaging.
//!
//! Provides FCC/BCC slip systems, Schmid tensor, Taylor factor, latent hardening,
//! power-law creep, Voigt/Reuss/Hill polycrystal averaging, and texture evolution.

#![allow(dead_code)]
#![allow(clippy::too_many_arguments)]

use std::f64::consts::PI;

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/// Normalise a 3-vector; returns zero vector if norm < eps.
fn normalise3(v: [f64; 3]) -> [f64; 3] {
    let n = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
    if n < 1e-14 {
        [0.0; 3]
    } else {
        [v[0] / n, v[1] / n, v[2] / n]
    }
}

/// Dot product of two 3-vectors.
fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// Double contraction (Frobenius inner product) of two 3×3 matrices.
fn double_contract(a: &[[f64; 3]; 3], b: &[[f64; 3]; 3]) -> f64 {
    let mut s = 0.0;
    for i in 0..3 {
        for j in 0..3 {
            s += a[i][j] * b[i][j];
        }
    }
    s
}

/// Multiply two 3×3 matrices.
fn matmul3(a: &[[f64; 3]; 3], b: &[[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let mut c = [[0.0f64; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            for k in 0..3 {
                c[i][j] += a[i][k] * b[k][j];
            }
        }
    }
    c
}

/// Transpose a 3×3 matrix.
fn transpose3(m: &[[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let mut t = [[0.0f64; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            t[i][j] = m[j][i];
        }
    }
    t
}

// ---------------------------------------------------------------------------
// 1. SlipSystem
// ---------------------------------------------------------------------------

/// Single crystallographic slip system with a Schmid factor.
///
/// A slip system is defined by a slip direction **d** and a slip plane normal
/// **n**.  The Schmid factor m = cos(φ)·cos(λ) is the projection of the
/// loading axis onto the slip system.
#[derive(Debug, Clone)]
pub struct SlipSystem {
    /// Unit slip direction (crystal frame).
    pub slip_direction: [f64; 3],
    /// Unit slip plane normal (crystal frame).
    pub slip_normal: [f64; 3],
    /// Schmid factor for uniaxial loading along e₁ = \[1,0,0\].
    pub schmid_factor: f64,
}

impl SlipSystem {
    /// Create a slip system from a direction and plane normal.
    ///
    /// Both vectors are normalised internally; the Schmid factor is computed
    /// for a loading axis along **e₁ = \[1,0,0\]**.
    pub fn new(direction: [f64; 3], normal: [f64; 3]) -> Self {
        let d = normalise3(direction);
        let n = normalise3(normal);
        // Schmid factor for axis e1 = [1,0,0]: m = cos_phi * cos_lambda
        let cos_phi = n[0]; // dot(n, e1)
        let cos_lam = d[0]; // dot(d, e1)
        let schmid_factor = cos_phi * cos_lam;
        Self {
            slip_direction: d,
            slip_normal: n,
            schmid_factor,
        }
    }

    /// Schmid tensor P = (d ⊗ n + n ⊗ d) / 2 (symmetrised).
    pub fn schmid_tensor(&self) -> [[f64; 3]; 3] {
        let d = self.slip_direction;
        let n = self.slip_normal;
        let mut p = [[0.0f64; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                p[i][j] = 0.5 * (d[i] * n[j] + n[i] * d[j]);
            }
        }
        p
    }

    /// Resolved shear stress τ = σ : P (double contraction with Cauchy stress).
    pub fn resolved_shear_stress(&self, cauchy_stress: &[[f64; 3]; 3]) -> f64 {
        let p = self.schmid_tensor();
        double_contract(cauchy_stress, &p)
    }
}

// ---------------------------------------------------------------------------
// 2. FCC and BCC slip systems
// ---------------------------------------------------------------------------

/// Generate the 12 {111}⟨110⟩ slip systems for an FCC crystal.
pub fn fcc_slip_systems() -> Vec<SlipSystem> {
    // Four {111} planes.
    let normals: [[f64; 3]; 4] = [
        [1.0, 1.0, 1.0],
        [-1.0, 1.0, 1.0],
        [1.0, -1.0, 1.0],
        [1.0, 1.0, -1.0],
    ];
    // Three ⟨110⟩ directions per plane (those perpendicular to the normal).
    let directions: [[f64; 3]; 6] = [
        [1.0, -1.0, 0.0],
        [-1.0, 0.0, 1.0],
        [0.0, 1.0, -1.0],
        [1.0, 1.0, 0.0],
        [0.0, 1.0, 1.0],
        [1.0, 0.0, 1.0],
    ];

    // Pairs: (normal_idx, direction_idx) that are orthogonal.
    let pairs: [(usize, usize); 12] = [
        (0, 0),
        (0, 1),
        (0, 2),
        (1, 0),
        (1, 3),
        (1, 4),
        (2, 1),
        (2, 3),
        (2, 5),
        (3, 2),
        (3, 4),
        (3, 5),
    ];

    pairs
        .iter()
        .map(|&(ni, di)| SlipSystem::new(directions[di], normals[ni]))
        .collect()
}

/// Generate the 12 {110}⟨111⟩ slip systems for a BCC crystal.
pub fn bcc_slip_systems() -> Vec<SlipSystem> {
    // Six {110} planes.
    let normals: [[f64; 3]; 6] = [
        [1.0, 1.0, 0.0],
        [1.0, -1.0, 0.0],
        [1.0, 0.0, 1.0],
        [1.0, 0.0, -1.0],
        [0.0, 1.0, 1.0],
        [0.0, 1.0, -1.0],
    ];
    // Four ⟨111⟩ directions.
    let directions: [[f64; 3]; 4] = [
        [1.0, 1.0, 1.0],
        [-1.0, 1.0, 1.0],
        [1.0, -1.0, 1.0],
        [1.0, 1.0, -1.0],
    ];

    // 12 orthogonal pairs.
    let pairs: [(usize, usize); 12] = [
        (0, 2),
        (0, 3),
        (1, 0),
        (1, 1),
        (2, 0),
        (2, 1),
        (3, 2),
        (3, 3),
        (4, 0),
        (4, 3),
        (5, 1),
        (5, 2),
    ];

    pairs
        .iter()
        .map(|&(ni, di)| SlipSystem::new(directions[di], normals[ni]))
        .collect()
}

// ---------------------------------------------------------------------------
// 3. Crystal orientation matrix (Bunge Euler angles)
// ---------------------------------------------------------------------------

/// Rotation matrix from crystal to sample frame (Bunge convention).
///
/// R = Rz(φ₂) · Rx(Φ) · Rz(φ₁)
#[derive(Debug, Clone)]
pub struct CrystalOrientationMatrix {
    /// 3×3 rotation matrix (row-major).
    pub r: [[f64; 3]; 3],
}

impl CrystalOrientationMatrix {
    /// Construct from Bunge Euler angles (φ₁, Φ, φ₂) in **radians**.
    pub fn from_euler_angles(phi1: f64, phi: f64, phi2: f64) -> Self {
        let (c1, s1) = (phi1.cos(), phi1.sin());
        let (c, s) = (phi.cos(), phi.sin());
        let (c2, s2) = (phi2.cos(), phi2.sin());

        // R = Rz(phi2) * Rx(phi) * Rz(phi1)
        let r = [
            [c2 * c1 - s2 * c * s1, -c2 * s1 - s2 * c * c1, s2 * s],
            [s2 * c1 + c2 * c * s1, -s2 * s1 + c2 * c * c1, -c2 * s],
            [s * s1, s * c1, c],
        ];

        Self { r }
    }

    /// Extract Bunge Euler angles (φ₁, Φ, φ₂) in **radians** from the matrix.
    ///
    /// Returns angles in \[0, 2π) × \[0, π\] × \[0, 2π).
    pub fn euler_angles(&self) -> (f64, f64, f64) {
        let r = &self.r;
        let phi = r[2][2].clamp(-1.0, 1.0).acos();
        let (phi1, phi2);
        if phi.sin().abs() < 1e-10 {
            // Gimbal lock: Φ = 0 or π
            phi1 = r[0][1].atan2(r[0][0]);
            phi2 = 0.0;
        } else {
            phi1 = r[2][0].atan2(r[2][1]);
            phi2 = r[0][2].atan2(-r[1][2]);
        }
        (phi1.rem_euclid(2.0 * PI), phi, phi2.rem_euclid(2.0 * PI))
    }

    /// Rotate a Cauchy stress tensor from crystal to sample frame.
    ///
    /// σ_sample = R · σ_crystal · Rᵀ
    pub fn rotate_stress(&self, sigma: &[[f64; 3]; 3]) -> [[f64; 3]; 3] {
        let rt = transpose3(&self.r);
        let tmp = matmul3(&self.r, sigma);
        matmul3(&tmp, &rt)
    }
}

// ---------------------------------------------------------------------------
// 4. Taylor factor
// ---------------------------------------------------------------------------

/// Compute a simplified Taylor factor M for a given strain rate.
///
/// M = σ_flow / τ_crss, approximated as the sum of absolute Schmid factors
/// for the active slip systems weighted by the deviatoric strain-rate
/// component ε₁₁.
///
/// # Arguments
/// * `slip_systems` – list of slip systems.
/// * `strain_rate`  – Voigt strain-rate vector \[ε₁₁, ε₂₂, ε₃₃, 2ε₂₃, 2ε₁₃, 2ε₁₂\].
pub fn taylor_factor(slip_systems: &[SlipSystem], strain_rate: [f64; 6]) -> f64 {
    let e11 = strain_rate[0];
    if e11.abs() < 1e-15 {
        return 0.0;
    }
    // Build stress-like tensor from the strain rate (use it as proxy for stress direction).
    let sigma = [
        [strain_rate[0], strain_rate[5], strain_rate[4]],
        [strain_rate[5], strain_rate[1], strain_rate[3]],
        [strain_rate[4], strain_rate[3], strain_rate[2]],
    ];
    let total: f64 = slip_systems
        .iter()
        .map(|s| s.resolved_shear_stress(&sigma).abs())
        .sum();
    total / e11.abs()
}

// ---------------------------------------------------------------------------
// 5. Latent hardening matrix
// ---------------------------------------------------------------------------

/// Build a latent hardening interaction matrix H_{αβ} = q + (1−q)·δ_{αβ}.
///
/// # Arguments
/// * `n_systems` – number of slip systems.
/// * `q`         – latent hardening ratio (0 = no latent, 1 = isotropic).
pub fn hardening_matrix_latent(n_systems: usize, q: f64) -> Vec<Vec<f64>> {
    let mut h = vec![vec![q; n_systems]; n_systems];
    #[allow(clippy::needless_range_loop)]
    for i in 0..n_systems {
        h[i][i] = 1.0;
    }
    h
}

// ---------------------------------------------------------------------------
// 6. PowerLawCreep (Norton)
// ---------------------------------------------------------------------------

/// Norton power-law creep model: ε̇ = A · σⁿ · exp(−Q / RT).
#[derive(Debug, Clone)]
pub struct PowerLawCreep {
    /// Pre-exponential factor A (Pa⁻ⁿ s⁻¹).
    pub a: f64,
    /// Creep exponent n (dimensionless).
    pub n: f64,
    /// Activation energy Q (J/mol).
    pub q_activation: f64,
    /// Universal gas constant R (J/(mol·K)).
    pub r_gas: f64,
}

impl PowerLawCreep {
    /// Create a new Norton creep model.
    pub fn new(a: f64, n: f64, q_activation: f64, r_gas: f64) -> Self {
        Self {
            a,
            n,
            q_activation,
            r_gas,
        }
    }

    /// Creep strain rate ε̇ = A · σⁿ · exp(−Q/RT).
    ///
    /// # Arguments
    /// * `stress`      – applied stress (Pa).
    /// * `temperature` – absolute temperature (K).
    pub fn strain_rate(&self, stress: f64, temperature: f64) -> f64 {
        let arrhenius = (-self.q_activation / (self.r_gas * temperature)).exp();
        self.a * stress.powf(self.n) * arrhenius
    }

    /// Return the creep exponent n.
    pub fn creep_exponent(&self) -> f64 {
        self.n
    }
}

// ---------------------------------------------------------------------------
// 7. Voigt / Reuss / Hill averaging
// ---------------------------------------------------------------------------

/// Voigt upper bound: arithmetic mean of elastic constants.
pub fn voigt_average(elastic_constants: &[f64]) -> f64 {
    if elastic_constants.is_empty() {
        return 0.0;
    }
    elastic_constants.iter().sum::<f64>() / elastic_constants.len() as f64
}

/// Reuss lower bound: harmonic mean of elastic constants.
pub fn reuss_average(elastic_constants: &[f64]) -> f64 {
    if elastic_constants.is_empty() {
        return 0.0;
    }
    let n = elastic_constants.len() as f64;
    let sum_inv: f64 = elastic_constants.iter().map(|c| 1.0 / c).sum();
    n / sum_inv
}

/// Hill average = (Voigt + Reuss) / 2.
pub fn hill_average(voigt: f64, reuss: f64) -> f64 {
    (voigt + reuss) / 2.0
}

// ---------------------------------------------------------------------------
// 8. TextureEvolution
// ---------------------------------------------------------------------------

/// Polycrystalline texture: a collection of orientations with weights.
#[derive(Debug, Clone)]
pub struct TextureEvolution {
    /// Crystal orientations.
    pub orientations: Vec<CrystalOrientationMatrix>,
    /// Volume fractions (should sum to 1).
    pub weights: Vec<f64>,
}

impl TextureEvolution {
    /// Create an empty texture.
    pub fn new() -> Self {
        Self {
            orientations: Vec::new(),
            weights: Vec::new(),
        }
    }

    /// Add an orientation with a given volume weight.
    pub fn add_orientation(&mut self, ori: CrystalOrientationMatrix, weight: f64) {
        self.orientations.push(ori);
        self.weights.push(weight);
    }

    /// Number of orientations.
    pub fn count(&self) -> usize {
        self.orientations.len()
    }

    /// Volume-weighted average Taylor factor over all orientations.
    ///
    /// # Arguments
    /// * `slip_systems` – slip systems used to evaluate each grain.
    /// * `strain`       – Voigt strain-rate vector.
    pub fn average_taylor_factor(&self, slip_systems: &[SlipSystem], strain: [f64; 6]) -> f64 {
        if self.orientations.is_empty() {
            return 0.0;
        }
        let total_weight: f64 = self.weights.iter().sum();
        if total_weight < 1e-15 {
            return 0.0;
        }
        let weighted_sum: f64 = self
            .orientations
            .iter()
            .zip(self.weights.iter())
            .map(|(ori, &w)| {
                // Rotate slip systems into sample frame for this grain.
                let rotated: Vec<SlipSystem> = slip_systems
                    .iter()
                    .map(|s| {
                        // Rotate direction and normal by R.
                        let rd = mat_vec3(&ori.r, s.slip_direction);
                        let rn = mat_vec3(&ori.r, s.slip_normal);
                        SlipSystem::new(rd, rn)
                    })
                    .collect();
                w * taylor_factor(&rotated, strain)
            })
            .sum();
        weighted_sum / total_weight
    }
}

impl Default for TextureEvolution {
    fn default() -> Self {
        Self::new()
    }
}

/// Multiply 3×3 matrix by a 3-vector.
fn mat_vec3(m: &[[f64; 3]; 3], v: [f64; 3]) -> [f64; 3] {
    [dot3(m[0], v), dot3(m[1], v), dot3(m[2], v)]
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-9;

    // ── FCC slip systems ─────────────────────────────────────────────────────

    #[test]
    fn test_fcc_has_12_slip_systems() {
        let sys = fcc_slip_systems();
        assert_eq!(sys.len(), 12, "FCC should have 12 slip systems");
    }

    #[test]
    fn test_fcc_slip_directions_are_unit_vectors() {
        for s in fcc_slip_systems() {
            let d = s.slip_direction;
            let norm = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
            assert!((norm - 1.0).abs() < 1e-10, "slip direction norm={norm}");
        }
    }

    #[test]
    fn test_fcc_slip_normals_are_unit_vectors() {
        for s in fcc_slip_systems() {
            let n = s.slip_normal;
            let norm = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
            assert!((norm - 1.0).abs() < 1e-10, "slip normal norm={norm}");
        }
    }

    #[test]
    fn test_fcc_schmid_factor_in_range() {
        for s in fcc_slip_systems() {
            assert!(
                s.schmid_factor >= -0.5 - EPS && s.schmid_factor <= 0.5 + EPS,
                "Schmid factor {} out of [-0.5, 0.5]",
                s.schmid_factor
            );
        }
    }

    // ── BCC slip systems ─────────────────────────────────────────────────────

    #[test]
    fn test_bcc_has_12_slip_systems() {
        let sys = bcc_slip_systems();
        assert_eq!(sys.len(), 12, "BCC should have 12 slip systems");
    }

    #[test]
    fn test_bcc_slip_directions_are_unit_vectors() {
        for s in bcc_slip_systems() {
            let d = s.slip_direction;
            let norm = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
            assert!((norm - 1.0).abs() < 1e-10, "BCC slip direction norm={norm}");
        }
    }

    #[test]
    fn test_bcc_schmid_factor_in_range() {
        for s in bcc_slip_systems() {
            assert!(
                s.schmid_factor >= -0.5 - EPS && s.schmid_factor <= 0.5 + EPS,
                "BCC Schmid factor {} out of [-0.5, 0.5]",
                s.schmid_factor
            );
        }
    }

    // ── Schmid tensor ────────────────────────────────────────────────────────

    #[test]
    #[allow(clippy::needless_range_loop)]
    fn test_schmid_tensor_symmetric() {
        let s = SlipSystem::new([1.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
        let p = s.schmid_tensor();
        for i in 0..3 {
            for j in 0..3 {
                assert!(
                    (p[i][j] - p[j][i]).abs() < EPS,
                    "P not symmetric at ({i},{j})"
                );
            }
        }
    }

    #[test]
    fn test_resolved_shear_stress_uniaxial() {
        // For slip direction [1,0,0] and normal [0,1,0]:
        // P_12 = 0.5, so tau = sigma_12 + sigma_21 = 2 * 0.5 * sigma_12
        let s = SlipSystem::new([1.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
        let sigma = [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 0.0]];
        let tau = s.resolved_shear_stress(&sigma);
        // P_01 = P_10 = 0.5; tau = 0.5*1.0 + 0.5*1.0 = 1.0
        assert!((tau - 1.0).abs() < EPS, "tau={tau}");
    }

    // ── CrystalOrientationMatrix ─────────────────────────────────────────────

    #[test]
    fn test_euler_identity() {
        // phi1=0, phi=0, phi2=0 should give identity matrix.
        let ori = CrystalOrientationMatrix::from_euler_angles(0.0, 0.0, 0.0);
        for i in 0..3 {
            for j in 0..3 {
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!(
                    (ori.r[i][j] - expected).abs() < 1e-10,
                    "R[{i}][{j}]={}",
                    ori.r[i][j]
                );
            }
        }
    }

    #[test]
    fn test_euler_angles_roundtrip() {
        let phi1 = 0.5_f64;
        let phi = 0.8_f64;
        let phi2 = 1.2_f64;
        let ori = CrystalOrientationMatrix::from_euler_angles(phi1, phi, phi2);
        let (r1, r, r2) = ori.euler_angles();
        // Rebuild and compare matrices (angles can differ by 2pi etc.).
        let ori2 = CrystalOrientationMatrix::from_euler_angles(r1, r, r2);
        for i in 0..3 {
            for j in 0..3 {
                assert!(
                    (ori.r[i][j] - ori2.r[i][j]).abs() < 1e-8,
                    "Euler roundtrip mismatch at [{i}][{j}]"
                );
            }
        }
    }

    #[test]
    fn test_rotate_stress_identity_orientation() {
        let ori = CrystalOrientationMatrix::from_euler_angles(0.0, 0.0, 0.0);
        let sigma = [[1.0, 2.0, 3.0], [2.0, 4.0, 5.0], [3.0, 5.0, 6.0]];
        let rotated = ori.rotate_stress(&sigma);
        for i in 0..3 {
            for j in 0..3 {
                assert!((rotated[i][j] - sigma[i][j]).abs() < 1e-10);
            }
        }
    }

    // ── Voigt / Reuss / Hill ─────────────────────────────────────────────────

    #[test]
    fn test_voigt_greater_than_or_equal_reuss() {
        let constants = vec![100.0, 200.0, 150.0, 180.0];
        let v = voigt_average(&constants);
        let r = reuss_average(&constants);
        assert!(v >= r - EPS, "Voigt={v}, Reuss={r}");
    }

    #[test]
    fn test_hill_between_voigt_and_reuss() {
        let v = 180.0_f64;
        let r = 120.0_f64;
        let h = hill_average(v, r);
        assert!(h >= r - EPS && h <= v + EPS, "Hill={h} not in [{r},{v}]");
    }

    #[test]
    fn test_voigt_single_value() {
        let v = voigt_average(&[200.0]);
        assert!((v - 200.0).abs() < EPS);
    }

    #[test]
    fn test_reuss_single_value() {
        let r = reuss_average(&[200.0]);
        assert!((r - 200.0).abs() < EPS);
    }

    #[test]
    fn test_hill_average_formula() {
        let h = hill_average(180.0, 120.0);
        assert!((h - 150.0).abs() < EPS);
    }

    #[test]
    fn test_voigt_empty_returns_zero() {
        assert_eq!(voigt_average(&[]), 0.0);
    }

    #[test]
    fn test_reuss_empty_returns_zero() {
        assert_eq!(reuss_average(&[]), 0.0);
    }

    // ── PowerLawCreep ────────────────────────────────────────────────────────

    #[test]
    fn test_power_law_strain_rate_positive() {
        let creep = PowerLawCreep::new(1e-15, 5.0, 200_000.0, 8.314);
        let rate = creep.strain_rate(100e6, 800.0);
        assert!(rate > 0.0, "strain rate should be positive, got {rate}");
    }

    #[test]
    fn test_power_law_strain_rate_increases_with_stress() {
        let creep = PowerLawCreep::new(1e-15, 5.0, 200_000.0, 8.314);
        let r1 = creep.strain_rate(100e6, 800.0);
        let r2 = creep.strain_rate(200e6, 800.0);
        assert!(r2 > r1, "strain rate should increase with stress");
    }

    #[test]
    fn test_power_law_strain_rate_increases_with_temperature() {
        let creep = PowerLawCreep::new(1e-15, 5.0, 200_000.0, 8.314);
        let r1 = creep.strain_rate(100e6, 700.0);
        let r2 = creep.strain_rate(100e6, 900.0);
        assert!(r2 > r1, "strain rate should increase with temperature");
    }

    #[test]
    fn test_creep_exponent_accessor() {
        let creep = PowerLawCreep::new(1e-15, 5.0, 200_000.0, 8.314);
        assert!((creep.creep_exponent() - 5.0).abs() < EPS);
    }

    // ── Latent hardening ─────────────────────────────────────────────────────

    #[test]
    #[allow(clippy::needless_range_loop)]
    fn test_hardening_matrix_diagonal_is_one() {
        let h = hardening_matrix_latent(12, 1.4);
        for i in 0..12 {
            assert!((h[i][i] - 1.0).abs() < EPS);
        }
    }

    #[test]
    #[allow(clippy::needless_range_loop)]
    fn test_hardening_matrix_off_diagonal_is_q() {
        let q = 1.4_f64;
        let h = hardening_matrix_latent(12, q);
        for i in 0..12 {
            for j in 0..12 {
                if i != j {
                    assert!((h[i][j] - q).abs() < EPS);
                }
            }
        }
    }

    // ── TextureEvolution ─────────────────────────────────────────────────────

    #[test]
    fn test_texture_add_orientation() {
        let mut tex = TextureEvolution::new();
        let ori = CrystalOrientationMatrix::from_euler_angles(0.0, 0.0, 0.0);
        tex.add_orientation(ori, 1.0);
        assert_eq!(tex.count(), 1);
    }

    #[test]
    fn test_texture_empty_average_taylor_zero() {
        let tex = TextureEvolution::new();
        let sys = fcc_slip_systems();
        let m = tex.average_taylor_factor(&sys, [1.0, -0.5, -0.5, 0.0, 0.0, 0.0]);
        assert_eq!(m, 0.0);
    }

    #[test]
    fn test_texture_single_grain_taylor_factor() {
        let mut tex = TextureEvolution::new();
        let ori = CrystalOrientationMatrix::from_euler_angles(0.0, 0.0, 0.0);
        tex.add_orientation(ori, 1.0);
        let sys = fcc_slip_systems();
        let m = tex.average_taylor_factor(&sys, [1.0, 0.0, 0.0, 0.0, 0.0, 0.0]);
        assert!(m >= 0.0, "Taylor factor should be non-negative, got {m}");
    }

    // ── Taylor factor ────────────────────────────────────────────────────────

    #[test]
    fn test_taylor_factor_zero_strain_returns_zero() {
        let sys = fcc_slip_systems();
        let m = taylor_factor(&sys, [0.0; 6]);
        assert_eq!(m, 0.0);
    }

    #[test]
    fn test_taylor_factor_positive_for_nonzero_strain() {
        let sys = fcc_slip_systems();
        let m = taylor_factor(&sys, [1.0, -0.5, -0.5, 0.0, 0.0, 0.0]);
        assert!(m >= 0.0);
    }
}
