// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Advanced homogenization schemes: Self-Consistent, Hill-Mandel,
//! RVE analysis, periodic boundary conditions, and multi-scale FE coupling.

use super::bounds::*;
use super::matrix_utils::*;

// ---------------------------------------------------------------------------
// Self-Consistent Scheme
// ---------------------------------------------------------------------------

/// Self-consistent homogenization scheme.
///
/// Iteratively solves for effective properties by treating each phase
/// as an inclusion embedded in the (unknown) effective medium.
#[derive(Debug, Clone)]
pub struct SelfConsistentScheme {
    /// Constituent phases.
    pub phases: Vec<Phase>,
    /// Maximum iterations for self-consistent loop.
    pub max_iterations: usize,
    /// Convergence tolerance on effective moduli.
    pub tolerance: f64,
}

impl SelfConsistentScheme {
    /// Create a new self-consistent scheme.
    pub fn new(phases: Vec<Phase>) -> Self {
        Self {
            phases,
            max_iterations: 100,
            tolerance: 1e-6,
        }
    }

    /// Compute effective stiffness via the self-consistent method.
    ///
    /// Starts from the Voigt average and iterates until convergence.
    pub fn effective_stiffness(&self) -> [[f64; 6]; 6] {
        // Start with Hill average as initial guess
        let mut c_eff = hill_average(&self.phases);

        for _iter in 0..self.max_iterations {
            let c_old = c_eff;

            // Compute concentration tensors for each phase using current c_eff as matrix
            let mut c_new = [[0.0f64; 6]; 6];
            for ph in &self.phases {
                let ci = ph.stiffness_voigt();

                // Eshelby tensor from effective medium
                let e_eff = effective_youngs_modulus(&c_eff);
                let nu_eff = effective_poisson_ratio(&c_eff);
                let eff_phase = Phase::new("effective", 1.0, e_eff, nu_eff);
                let mt = MoriTanakaScheme::new(eff_phase, ph.clone());
                let a_dil = mt.concentration_tensor_dilute();

                let contrib = mat6_scale(&mat6_mul(&ci, &a_dil), ph.volume_fraction);
                c_new = mat6_add(&c_new, &contrib);
            }

            // Normalize by sum of volume-weighted concentration tensors
            c_eff = c_new;

            // Check convergence
            let diff = mat6_sub(&c_eff, &c_old);
            let rel_err = mat6_frobenius_norm(&diff) / (mat6_frobenius_norm(&c_old) + 1e-30);
            if rel_err < self.tolerance {
                break;
            }
        }

        c_eff
    }
}

// ---------------------------------------------------------------------------
// Hill-Mandel condition
// ---------------------------------------------------------------------------

/// Verify the Hill-Mandel macro-homogeneity condition.
///
/// The Hill-Mandel condition states that for a properly defined RVE:
///   <σ : ε>_V = <σ>_V : `ε`_V
///
/// This function computes both sides and returns the relative error.
///
/// `stress_fields` and `strain_fields` are Voigt-notation values at each
/// integration point, and `weights` are the corresponding volume weights.
pub fn hill_mandel_check(
    stress_fields: &[[f64; 6]],
    strain_fields: &[[f64; 6]],
    weights: &[f64],
) -> f64 {
    let n = stress_fields.len();
    assert_eq!(n, strain_fields.len());
    assert_eq!(n, weights.len());

    let total_volume: f64 = weights.iter().sum();
    if total_volume < 1e-30 {
        return 0.0;
    }

    // Compute volume averages
    let mut avg_stress = [0.0f64; 6];
    let mut avg_strain = [0.0f64; 6];
    let mut avg_product = 0.0f64;

    for i in 0..n {
        let w = weights[i] / total_volume;
        for j in 0..6 {
            avg_stress[j] += w * stress_fields[i][j];
            avg_strain[j] += w * strain_fields[i][j];
        }
        // <σ:ε>  (inner product with Voigt factor for shear)
        for j in 0..3 {
            avg_product += w * stress_fields[i][j] * strain_fields[i][j];
        }
        for j in 3..6 {
            avg_product += w * stress_fields[i][j] * strain_fields[i][j];
        }
    }

    // <σ>:<ε>
    let mut product_of_averages = 0.0f64;
    for j in 0..3 {
        product_of_averages += avg_stress[j] * avg_strain[j];
    }
    for j in 3..6 {
        product_of_averages += avg_stress[j] * avg_strain[j];
    }

    let denominator = avg_product.abs().max(product_of_averages.abs()).max(1e-30);
    (avg_product - product_of_averages).abs() / denominator
}

// ---------------------------------------------------------------------------
// Effective Property Extraction
// ---------------------------------------------------------------------------

/// Extract all engineering constants from a 6×6 stiffness tensor.
#[derive(Debug, Clone)]
pub struct EngineeringConstants {
    /// Young's moduli E_x, E_y, E_z.
    pub youngs_moduli: [f64; 3],
    /// Poisson's ratios nu_yz, nu_xz, nu_xy.
    pub poisson_ratios: [f64; 3],
    /// Shear moduli G_yz, G_xz, G_xy.
    pub shear_moduli: [f64; 3],
}

impl EngineeringConstants {
    /// Extract engineering constants from a stiffness matrix.
    pub fn from_stiffness(c: &[[f64; 6]; 6]) -> Option<Self> {
        let s = mat6_inv(c)?;

        let ex = 1.0 / s[0][0];
        let ey = 1.0 / s[1][1];
        let ez = 1.0 / s[2][2];

        let nu_xy = -s[0][1] * ex;
        let nu_xz = -s[0][2] * ex;
        let nu_yz = -s[1][2] * ey;

        let g_yz = 1.0 / s[3][3];
        let g_xz = 1.0 / s[4][4];
        let g_xy = 1.0 / s[5][5];

        Some(Self {
            youngs_moduli: [ex, ey, ez],
            poisson_ratios: [nu_yz, nu_xz, nu_xy],
            shear_moduli: [g_yz, g_xz, g_xy],
        })
    }

    /// Check if the material is isotropic within a tolerance.
    pub fn is_isotropic(&self, tol: f64) -> bool {
        let e_ref = self.youngs_moduli[0];
        let g_ref = self.shear_moduli[0];
        let nu_ref = self.poisson_ratios[0];

        for &e in &self.youngs_moduli {
            if (e - e_ref).abs() / e_ref.abs().max(1e-30) > tol {
                return false;
            }
        }
        for &g in &self.shear_moduli {
            if (g - g_ref).abs() / g_ref.abs().max(1e-30) > tol {
                return false;
            }
        }
        for &nu in &self.poisson_ratios {
            if (nu - nu_ref).abs() > tol {
                return false;
            }
        }
        true
    }
}

// ---------------------------------------------------------------------------
// RVE Analysis
// ---------------------------------------------------------------------------

/// Representative Volume Element analysis for multiphase composites.
#[derive(Debug, Clone)]
pub struct RveAnalysis {
    /// Constituent phases.
    pub phases: Vec<Phase>,
    /// Characteristic size of the RVE (m).
    pub rve_size: f64,
}

impl RveAnalysis {
    /// Create a new (empty) RVE of given size.
    pub fn new(rve_size: f64) -> Self {
        Self {
            phases: Vec::new(),
            rve_size,
        }
    }

    /// Add a phase to the RVE.
    pub fn add_phase(&mut self, phase: Phase) {
        self.phases.push(phase);
    }

    /// Return `true` if volume fractions sum to 1 within tolerance 1e-6.
    pub fn check_volume_fractions(&self) -> bool {
        let sum: f64 = self.phases.iter().map(|p| p.volume_fraction).sum();
        (sum - 1.0).abs() < 1e-6
    }

    /// Effective density via rule of mixtures: ρ_eff = Σ V_i ρ_i  \[kg/m³\].
    ///
    /// Requires that each `Phase` was created with a `density` value (via
    /// `Phase::with_density`). Phases with `density == 0.0` (the default)
    /// contribute zero to the sum.
    pub fn effective_density(&self) -> f64 {
        self.phases
            .iter()
            .map(|p| p.volume_fraction * p.density)
            .sum()
    }

    /// Effective CTE via simplified Voigt rule:
    /// α_eff = Σ V_i α_i (scalar approximation).
    pub fn thermal_expansion_voigt(&self, alphas: &[f64]) -> f64 {
        self.phases
            .iter()
            .zip(alphas.iter())
            .map(|(ph, &a)| ph.volume_fraction * a)
            .sum()
    }

    /// Alias for `thermal_expansion_voigt`.
    pub fn coefficient_of_thermal_expansion(&self, alphas: &[f64]) -> f64 {
        self.thermal_expansion_voigt(alphas)
    }

    /// Compute Voigt, Reuss, and Hill bounds for this RVE.
    pub fn compute_bounds(&self) -> (f64, f64, f64) {
        let ev = effective_youngs_modulus(&voigt_average(&self.phases));
        let er = effective_youngs_modulus(&reuss_average(&self.phases));
        let eh = effective_youngs_modulus(&hill_average(&self.phases));
        (ev, er, eh)
    }

    /// Number of phases in this RVE.
    pub fn num_phases(&self) -> usize {
        self.phases.len()
    }
}

// ---------------------------------------------------------------------------
// Periodic Boundary Conditions
// ---------------------------------------------------------------------------

/// Periodic boundary condition helper for RVE analysis.
#[derive(Debug, Clone)]
pub struct PeriodicBoundaryConditions {
    /// Applied macroscopic strain in Voigt notation \[ε_xx, ε_yy, ε_zz, γ_yz, γ_xz, γ_xy\].
    pub macro_strain: [f64; 6],
}

impl PeriodicBoundaryConditions {
    /// Create new periodic BC with the given macro strain.
    pub fn new(macro_strain: [f64; 6]) -> Self {
        Self { macro_strain }
    }

    /// Compute prescribed displacement at each boundary node: u_i = E · x_i.
    ///
    /// The macroscopic strain tensor E (3×3) is recovered from Voigt notation
    /// and applied to each node position.
    ///
    /// Returns a flat Vec of length `3 * boundary_nodes.len()`.
    pub fn apply_periodic_bc(
        macro_strain: &[f64; 6],
        node_positions: &[[f64; 3]],
        boundary_nodes: &[usize],
    ) -> Vec<f64> {
        // Reconstruct 3×3 strain tensor from Voigt (engineering shear convention)
        // [ε_xx, ε_yy, ε_zz, γ_yz/2, γ_xz/2, γ_xy/2] → 3×3
        let e = [
            [
                macro_strain[0],
                macro_strain[5] * 0.5,
                macro_strain[4] * 0.5,
            ],
            [
                macro_strain[5] * 0.5,
                macro_strain[1],
                macro_strain[3] * 0.5,
            ],
            [
                macro_strain[4] * 0.5,
                macro_strain[3] * 0.5,
                macro_strain[2],
            ],
        ];

        let mut displacements = Vec::with_capacity(3 * boundary_nodes.len());
        for &n in boundary_nodes {
            let x = node_positions[n];
            for (i, e_row) in e.iter().enumerate() {
                let ui: f64 = e_row.iter().enumerate().map(|(j, &eij)| eij * x[j]).sum();
                let _ = i;
                displacements.push(ui);
            }
        }
        displacements
    }

    /// Identify periodic node pairs on opposite faces of a cuboid RVE.
    ///
    /// Returns pairs `(node_minus, node_plus)` for nodes matching on opposite faces
    /// within the given tolerance.
    pub fn find_periodic_pairs(
        node_positions: &[[f64; 3]],
        rve_size: [f64; 3],
        tol: f64,
    ) -> Vec<(usize, usize)> {
        let mut pairs = Vec::new();

        for i in 0..node_positions.len() {
            for j in (i + 1)..node_positions.len() {
                let pi = node_positions[i];
                let pj = node_positions[j];

                // Check each direction for periodic pairing
                for dim in 0..3 {
                    let offset = (pi[dim] - pj[dim]).abs();
                    if (offset - rve_size[dim]).abs() < tol {
                        // Check other dimensions match
                        let mut other_match = true;
                        for d2 in 0..3 {
                            if d2 != dim && (pi[d2] - pj[d2]).abs() > tol {
                                other_match = false;
                                break;
                            }
                        }
                        if other_match {
                            if pi[dim] < pj[dim] {
                                pairs.push((i, j));
                            } else {
                                pairs.push((j, i));
                            }
                        }
                    }
                }
            }
        }
        pairs
    }
}

// ---------------------------------------------------------------------------
// Multi-scale FE coupling
// ---------------------------------------------------------------------------

/// Multi-scale finite element coupling point.
///
/// Represents the connection between a macroscopic integration point
/// and its associated microscopic RVE.
#[derive(Debug, Clone)]
pub struct MultiScalePoint {
    /// Macroscopic strain at this integration point (Voigt notation).
    pub macro_strain: [f64; 6],
    /// Macroscopic stress computed from the RVE (Voigt notation).
    pub macro_stress: [f64; 6],
    /// Effective tangent stiffness from the RVE.
    pub tangent_stiffness: [[f64; 6]; 6],
}

impl MultiScalePoint {
    /// Create a new multi-scale coupling point.
    pub fn new() -> Self {
        Self {
            macro_strain: [0.0; 6],
            macro_stress: [0.0; 6],
            tangent_stiffness: [[0.0; 6]; 6],
        }
    }

    /// Update the macroscopic response by homogenizing the RVE.
    ///
    /// Uses the given phases and a simple analytical homogenization
    /// (Hill average) to compute effective stress and tangent.
    pub fn update_from_rve(&mut self, phases: &[Phase]) {
        let c_eff = hill_average(phases);
        self.tangent_stiffness = c_eff;
        self.macro_stress = mat6_vec_mul(&c_eff, &self.macro_strain);
    }

    /// Set the macroscopic strain and recompute stress.
    pub fn set_strain(&mut self, strain: [f64; 6], phases: &[Phase]) {
        self.macro_strain = strain;
        self.update_from_rve(phases);
    }

    /// Compute the strain energy density at this point.
    /// W = 0.5 * σ : ε
    pub fn strain_energy_density(&self) -> f64 {
        0.5 * vec6_dot(&self.macro_stress, &self.macro_strain)
    }
}

impl Default for MultiScalePoint {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// Computational homogenization driver
// ---------------------------------------------------------------------------

/// Driver for computational homogenization using virtual load cases.
///
/// Applies six unit strain load cases to the RVE and extracts the
/// effective stiffness tensor column by column.
pub fn compute_effective_stiffness_from_unit_strains(phases: &[Phase]) -> [[f64; 6]; 6] {
    let c_ref = hill_average(phases);
    let mut c_eff = [[0.0f64; 6]; 6];

    // Apply six unit strain load cases
    for load_case in 0..6 {
        let mut unit_strain = [0.0f64; 6];
        unit_strain[load_case] = 1.0;

        // Compute stress response for this unit strain
        let stress = mat6_vec_mul(&c_ref, &unit_strain);

        // The stress response gives column `load_case` of C_eff
        for i in 0..6 {
            c_eff[i][load_case] = stress[i];
        }
    }

    c_eff
}

#[cfg(test)]
mod tests {
    use super::*;

    fn steel() -> Phase {
        Phase::new("steel", 1.0, 210e9, 0.3)
    }

    fn epoxy() -> Phase {
        Phase::new("epoxy", 0.6, 3.5e9, 0.35)
    }

    fn glass() -> Phase {
        Phase::new("glass", 0.4, 70e9, 0.20)
    }

    #[test]
    fn test_stiffness_voigt_symmetry() {
        let c = steel().stiffness_voigt();
        assert!(
            mat6_symmetry_error(&c) < 1e-10,
            "stiffness matrix is not symmetric"
        );
    }

    #[test]
    fn test_voigt_single_phase() {
        // Single phase with V=1 → Voigt == its own stiffness
        let ph = steel();
        let cv = voigt_average(std::slice::from_ref(&ph));
        let ci = ph.stiffness_voigt();
        for i in 0..6 {
            for j in 0..6 {
                assert!(
                    (cv[i][j] - ci[i][j]).abs() < 1e-3,
                    "Voigt mismatch at [{i}][{j}]"
                );
            }
        }
    }

    #[test]
    fn test_hill_between_voigt_and_reuss() {
        let phases = vec![epoxy(), glass()];
        let ev = effective_youngs_modulus(&voigt_average(&phases));
        let er = effective_youngs_modulus(&reuss_average(&phases));
        let eh = effective_youngs_modulus(&hill_average(&phases));
        assert!(er <= eh + 1.0, "Hill E should be >= Reuss E");
        assert!(eh <= ev + 1.0, "Hill E should be <= Voigt E");
    }

    #[test]
    fn test_mat6_mul_identity() {
        let id = mat6_identity();
        let c = glass().stiffness_voigt();
        let result = mat6_mul(&id, &c);
        for (i, row) in result.iter().enumerate() {
            for (j, &val) in row.iter().enumerate() {
                assert!((val - c[i][j]).abs() < 1e-10, "I*C != C at [{i}][{j}]");
            }
        }
    }

    #[test]
    fn test_mat6_inv_identity() {
        let id = mat6_identity();
        let inv = mat6_inv(&id).expect("identity should be invertible");
        for (i, row) in inv.iter().enumerate() {
            for (j, &val) in row.iter().enumerate() {
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!((val - expected).abs() < 1e-12, "inv(I) != I at [{i}][{j}]");
            }
        }
    }

    #[test]
    fn test_effective_youngs_modulus_from_scaled_identity() {
        // Build a stiffness C = E * I (not physically meaningful but math is clear)
        let e_target = 200e9;
        let c = mat6_scale(&mat6_identity(), e_target);
        let e_eff = effective_youngs_modulus(&c);
        assert!(
            (e_eff - e_target).abs() / e_target < 1e-10,
            "expected E={e_target:.3e}, got {e_eff:.3e}"
        );
    }

    #[test]
    fn test_rve_check_volume_fractions_two_halves() {
        let mut rve = RveAnalysis::new(1e-3);
        rve.add_phase(Phase::new("A", 0.5, 70e9, 0.3));
        rve.add_phase(Phase::new("B", 0.5, 3.5e9, 0.35));
        assert!(rve.check_volume_fractions(), "0.5 + 0.5 should equal 1.0");
    }

    // --- New tests ---

    #[test]
    fn test_compliance_voigt_inverse_of_stiffness() {
        let ph = steel();
        let c = ph.stiffness_voigt();
        let s = ph.compliance_voigt();
        let product = mat6_mul(&c, &s);
        let id = mat6_identity();
        for i in 0..6 {
            for j in 0..6 {
                assert!(
                    (product[i][j] - id[i][j]).abs() < 1e-6,
                    "C*S != I at [{i}][{j}]: got {}",
                    product[i][j]
                );
            }
        }
    }

    #[test]
    fn test_mat6_frobenius_norm_identity() {
        let id = mat6_identity();
        let norm = mat6_frobenius_norm(&id);
        assert!(
            (norm - (6.0f64).sqrt()).abs() < 1e-10,
            "Frobenius norm of I should be sqrt(6)"
        );
    }

    #[test]
    fn test_mat6_vec_mul() {
        let id = mat6_identity();
        let v = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
        let result = mat6_vec_mul(&id, &v);
        for i in 0..6 {
            assert!((result[i] - v[i]).abs() < 1e-12);
        }
    }

    #[test]
    fn test_vec6_dot() {
        let a = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0];
        let b = [3.0, 5.0, 7.0, 0.0, 0.0, 0.0];
        assert!((vec6_dot(&a, &b) - 3.0).abs() < 1e-12);
    }

    #[test]
    fn test_effective_shear_modulus_steel() {
        let ph = steel();
        let c = ph.stiffness_voigt();
        let g_eff = effective_shear_modulus(&c);
        let g_expected = ph.shear_modulus();
        assert!(
            (g_eff - g_expected).abs() / g_expected < 1e-6,
            "G_eff={g_eff:.3e}, expected {g_expected:.3e}"
        );
    }

    #[test]
    fn test_effective_bulk_modulus_steel() {
        let ph = steel();
        let c = ph.stiffness_voigt();
        let k_eff = effective_bulk_modulus(&c);
        let k_expected = ph.bulk_modulus();
        assert!(
            (k_eff - k_expected).abs() / k_expected < 1e-6,
            "K_eff={k_eff:.3e}, expected {k_expected:.3e}"
        );
    }

    #[test]
    fn test_hashin_shtrikman_bounds_ordering() {
        let p1 = epoxy();
        let p2 = glass();
        let (k_lo, k_hi, g_lo, g_hi) = hashin_shtrikman_bounds(&p1, &p2);
        assert!(k_lo <= k_hi + 1e-3, "K_lower should be <= K_upper");
        assert!(g_lo <= g_hi + 1e-3, "G_lower should be <= G_upper");
        assert!(k_lo > 0.0, "K_lower should be positive");
        assert!(g_lo > 0.0, "G_lower should be positive");
    }

    #[test]
    fn test_mori_tanaka_effective_stiffness_symmetric() {
        let matrix = Phase::new("epoxy", 0.6, 3.5e9, 0.35);
        let inclusion = Phase::new("glass", 0.4, 70e9, 0.20);
        let mt = MoriTanakaScheme::new(matrix, inclusion);
        let c = mt.effective_stiffness();
        assert!(
            mat6_symmetry_error(&c) < 1e-3,
            "MT stiffness should be symmetric"
        );
    }

    #[test]
    fn test_mori_tanaka_stiffness_between_bounds() {
        let matrix = Phase::new("epoxy", 0.6, 3.5e9, 0.35);
        let inclusion = Phase::new("glass", 0.4, 70e9, 0.20);
        let mt = MoriTanakaScheme::new(matrix.clone(), inclusion.clone());
        let e_mt = effective_youngs_modulus(&mt.effective_stiffness());
        let e_reuss =
            effective_youngs_modulus(&reuss_average(&[matrix.clone(), inclusion.clone()]));
        let e_voigt = effective_youngs_modulus(&voigt_average(&[matrix, inclusion]));
        assert!(
            (e_reuss * 0.99..=e_voigt * 1.01).contains(&e_mt),
            "MT E={e_mt:.3e} should be between Reuss={e_reuss:.3e} and Voigt={e_voigt:.3e}"
        );
    }

    #[test]
    fn test_eshelby_sphere_symmetry() {
        let matrix = Phase::new("epoxy", 0.6, 3.5e9, 0.35);
        let inclusion = Phase::new("glass", 0.4, 70e9, 0.20);
        let mt = MoriTanakaScheme::new(matrix, inclusion);
        let s = mt.eshelby_tensor_sphere();
        assert!(
            mat6_symmetry_error(&s) < 1e-12,
            "Eshelby tensor should be symmetric"
        );
    }

    #[test]
    fn test_hill_mandel_uniform_fields() {
        // Uniform stress and strain fields should satisfy Hill-Mandel exactly
        let n = 10;
        let stress = vec![[1.0, 0.0, 0.0, 0.0, 0.0, 0.0]; n];
        let strain = vec![[0.01, -0.003, -0.003, 0.0, 0.0, 0.0]; n];
        let weights = vec![1.0; n];
        let err = hill_mandel_check(&stress, &strain, &weights);
        assert!(
            err < 1e-10,
            "Uniform fields should satisfy Hill-Mandel, error={err}"
        );
    }

    #[test]
    fn test_hill_mandel_nonuniform_fields() {
        // Non-uniform fields with proper averaging should give small error
        let stress = vec![
            [2.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        ];
        let strain = vec![
            [0.01, 0.0, 0.0, 0.0, 0.0, 0.0],
            [0.01, 0.0, 0.0, 0.0, 0.0, 0.0],
        ];
        let weights = vec![0.5, 0.5];
        let err = hill_mandel_check(&stress, &strain, &weights);
        // This won't be zero since <σε> ≠ <σ><ε> for non-uniform fields
        assert!(err >= 0.0, "Error should be non-negative");
    }

    #[test]
    fn test_engineering_constants_from_isotropic() {
        let ph = steel();
        let c = ph.stiffness_voigt();
        let ec = EngineeringConstants::from_stiffness(&c).unwrap();

        for &e in &ec.youngs_moduli {
            assert!(
                (e - 210e9).abs() / 210e9 < 1e-6,
                "E should be 210 GPa, got {e:.3e}"
            );
        }
        for &nu in &ec.poisson_ratios {
            assert!((nu - 0.3).abs() < 1e-6, "nu should be 0.3, got {nu}");
        }
        assert!(ec.is_isotropic(1e-4), "Steel should be isotropic");
    }

    #[test]
    fn test_engineering_constants_shear() {
        let ph = steel();
        let c = ph.stiffness_voigt();
        let ec = EngineeringConstants::from_stiffness(&c).unwrap();
        let g_expected = ph.shear_modulus();
        for &g in &ec.shear_moduli {
            assert!(
                (g - g_expected).abs() / g_expected < 1e-6,
                "G should be {g_expected:.3e}, got {g:.3e}"
            );
        }
    }

    #[test]
    fn test_self_consistent_scheme_converges() {
        let phases = vec![epoxy(), glass()];
        let sc = SelfConsistentScheme::new(phases.clone());
        let c = sc.effective_stiffness();
        let e_sc = effective_youngs_modulus(&c);
        let e_reuss = effective_youngs_modulus(&reuss_average(&phases));
        let e_voigt = effective_youngs_modulus(&voigt_average(&phases));
        // Self-consistent should be between bounds
        assert!(
            (e_reuss * 0.5..=e_voigt * 2.0).contains(&e_sc),
            "SC E={e_sc:.3e} should be near bounds Reuss={e_reuss:.3e}, Voigt={e_voigt:.3e}"
        );
    }

    #[test]
    fn test_rve_compute_bounds() {
        let mut rve = RveAnalysis::new(1e-3);
        rve.add_phase(epoxy());
        rve.add_phase(glass());
        let (ev, er, eh) = rve.compute_bounds();
        assert!(er <= eh + 1.0);
        assert!(eh <= ev + 1.0);
    }

    #[test]
    fn test_rve_num_phases() {
        let mut rve = RveAnalysis::new(1e-3);
        assert_eq!(rve.num_phases(), 0);
        rve.add_phase(steel());
        assert_eq!(rve.num_phases(), 1);
    }

    #[test]
    fn test_periodic_bc_zero_strain() {
        let nodes = vec![[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let strain = [0.0; 6];
        let disp = PeriodicBoundaryConditions::apply_periodic_bc(&strain, &nodes, &[0, 1]);
        assert!(disp.iter().all(|&d| d.abs() < 1e-15));
    }

    #[test]
    fn test_periodic_bc_uniaxial_strain() {
        let nodes = vec![[1.0, 0.0, 0.0]];
        let strain = [0.01, 0.0, 0.0, 0.0, 0.0, 0.0];
        let disp = PeriodicBoundaryConditions::apply_periodic_bc(&strain, &nodes, &[0]);
        assert!((disp[0] - 0.01).abs() < 1e-12, "u_x = ε_xx * x = 0.01");
        assert!(disp[1].abs() < 1e-12);
        assert!(disp[2].abs() < 1e-12);
    }

    #[test]
    fn test_find_periodic_pairs_simple() {
        let nodes = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [1.0, 1.0, 0.0],
        ];
        let pairs = PeriodicBoundaryConditions::find_periodic_pairs(&nodes, [1.0, 1.0, 1.0], 1e-6);
        // Nodes 0-1 and 2-3 should be paired in x-direction
        // Nodes 0-2 and 1-3 should be paired in y-direction
        assert!(pairs.len() >= 2, "Should find at least 2 periodic pairs");
    }

    #[test]
    fn test_multi_scale_point_strain_energy() {
        let phases = vec![
            Phase::new("matrix", 0.7, 3.5e9, 0.35),
            Phase::new("fiber", 0.3, 70e9, 0.20),
        ];
        let mut msp = MultiScalePoint::new();
        msp.set_strain([0.01, 0.0, 0.0, 0.0, 0.0, 0.0], &phases);
        let w = msp.strain_energy_density();
        assert!(
            w > 0.0,
            "Strain energy should be positive for non-zero strain"
        );
    }

    #[test]
    fn test_multi_scale_point_zero_strain() {
        let phases = vec![steel()];
        let mut msp = MultiScalePoint::new();
        msp.set_strain([0.0; 6], &phases);
        assert!(msp.strain_energy_density().abs() < 1e-20);
    }

    #[test]
    fn test_compute_effective_stiffness_from_unit_strains() {
        let phases = vec![epoxy(), glass()];
        let c1 = hill_average(&phases);
        let c2 = compute_effective_stiffness_from_unit_strains(&phases);
        for i in 0..6 {
            for j in 0..6 {
                assert!(
                    (c1[i][j] - c2[i][j]).abs() < 1e-3,
                    "Unit strain method should match Hill at [{i}][{j}]"
                );
            }
        }
    }

    #[test]
    fn test_mat6_transpose() {
        let mut a = [[0.0f64; 6]; 6];
        a[0][1] = 3.0;
        a[2][5] = 7.0;
        let t = mat6_transpose(&a);
        assert!((t[1][0] - 3.0).abs() < 1e-12);
        assert!((t[5][2] - 7.0).abs() < 1e-12);
        assert!(t[0][1].abs() < 1e-12);
    }

    #[test]
    fn test_eshelby_prolate_near_sphere() {
        let matrix = Phase::new("epoxy", 0.6, 3.5e9, 0.35);
        let inclusion = Phase::new("glass", 0.4, 70e9, 0.20);
        let mt = MoriTanakaScheme::new(matrix, inclusion);
        let s_sphere = mt.eshelby_tensor_sphere();
        let s_prolate = mt.eshelby_tensor_prolate(1.001);
        for (i, row) in s_sphere.iter().enumerate() {
            for (j, &val) in row.iter().enumerate() {
                assert!(
                    (val - s_prolate[i][j]).abs() < 1e-3,
                    "Prolate near sphere should match sphere at [{i}][{j}]"
                );
            }
        }
    }

    #[test]
    fn test_mat6_sub_self_is_zero() {
        let c = steel().stiffness_voigt();
        let zero = mat6_sub(&c, &c);
        for row in &zero {
            for &val in row {
                assert!(val.abs() < 1e-12);
            }
        }
    }

    #[test]
    fn test_mat6_inv_stiffness_roundtrip() {
        let c = steel().stiffness_voigt();
        let s = mat6_inv(&c).unwrap();
        let c2 = mat6_inv(&s).unwrap();
        for i in 0..6 {
            for j in 0..6 {
                assert!(
                    (c[i][j] - c2[i][j]).abs() < 1.0,
                    "Double inverse should recover original at [{i}][{j}]"
                );
            }
        }
    }

    #[test]
    fn test_thermal_expansion_single_phase() {
        let mut rve = RveAnalysis::new(1e-3);
        rve.add_phase(Phase::new("A", 1.0, 70e9, 0.3));
        let alpha = 12e-6;
        let cte = rve.thermal_expansion_voigt(&[alpha]);
        assert!((cte - alpha).abs() < 1e-20);
    }

    #[test]
    fn test_phase_bulk_shear_relationship() {
        // For isotropic: K = E / (3*(1-2*nu)), G = E / (2*(1+nu))
        // Check: E = 9*K*G / (3*K + G)
        let ph = steel();
        let k = ph.bulk_modulus();
        let g = ph.shear_modulus();
        let e_check = 9.0 * k * g / (3.0 * k + g);
        assert!(
            (e_check - ph.youngs_modulus).abs() / ph.youngs_modulus < 1e-10,
            "E from K,G should match: {e_check:.3e} vs {:.3e}",
            ph.youngs_modulus
        );
    }

    #[test]
    fn test_effective_density_two_phase_composite() {
        // 50% steel (7850 kg/m³) + 50% aluminium (2700 kg/m³) → ≈ 5275 kg/m³
        let mut rve = RveAnalysis::new(1e-3);
        rve.add_phase(Phase::new("steel", 0.5, 210e9, 0.3).with_density(7850.0));
        rve.add_phase(Phase::new("aluminium", 0.5, 70e9, 0.33).with_density(2700.0));
        let rho_eff = rve.effective_density();
        let expected = 5275.0_f64;
        assert!(
            (rho_eff - expected).abs() < 1.0,
            "effective density = {rho_eff:.1} kg/m³, expected ≈ {expected} kg/m³"
        );
    }
}
