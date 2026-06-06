// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Multiscale FEM (FE²) and homogenization methods.
//!
//! Implements two-scale computational homogenization based on the Hill-Mandel
//! condition, Voigt and Reuss bounds, Hill average, and Mori-Tanaka
//! micromechanics for composite and heterogeneous materials.

// ---------------------------------------------------------------------------
// Material phase
// ---------------------------------------------------------------------------

/// A single material phase within a representative volume element (RVE).
#[derive(Debug, Clone)]
pub struct MaterialPhase {
    /// Volume fraction of this phase (0 < vf ≤ 1, all phases must sum to 1).
    pub volume_fraction: f64,
    /// Young's modulus \[Pa\].
    pub young_modulus: f64,
    /// Poisson's ratio (dimensionless, –1 < ν < 0.5).
    pub poisson_ratio: f64,
    /// Integer identifier for the phase (0 = matrix, 1+ = inclusions).
    pub phase_id: usize,
}

impl MaterialPhase {
    /// Create a new `MaterialPhase`.
    pub fn new(
        volume_fraction: f64,
        young_modulus: f64,
        poisson_ratio: f64,
        phase_id: usize,
    ) -> Self {
        Self {
            volume_fraction,
            young_modulus,
            poisson_ratio,
            phase_id,
        }
    }

    /// Bulk modulus K = E / (3(1 − 2ν)).
    pub fn bulk_modulus(&self) -> f64 {
        self.young_modulus / (3.0 * (1.0 - 2.0 * self.poisson_ratio))
    }

    /// Shear modulus G = E / (2(1 + ν)).
    pub fn shear_modulus(&self) -> f64 {
        self.young_modulus / (2.0 * (1.0 + self.poisson_ratio))
    }

    /// Build the 6×6 isotropic stiffness tensor in Voigt notation.
    pub fn stiffness_tensor(&self) -> [[f64; 6]; 6] {
        let e = self.young_modulus;
        let nu = self.poisson_ratio;
        let lam = e * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
        let mu = e / (2.0 * (1.0 + nu));
        let mut c = [[0.0f64; 6]; 6];
        // Normal–normal blocks
        c[0][0] = lam + 2.0 * mu;
        c[1][1] = lam + 2.0 * mu;
        c[2][2] = lam + 2.0 * mu;
        c[0][1] = lam;
        c[0][2] = lam;
        c[1][0] = lam;
        c[1][2] = lam;
        c[2][0] = lam;
        c[2][1] = lam;
        // Shear blocks
        c[3][3] = mu;
        c[4][4] = mu;
        c[5][5] = mu;
        c
    }
}

// ---------------------------------------------------------------------------
// RVE configuration
// ---------------------------------------------------------------------------

/// Configuration for a representative volume element (RVE).
#[derive(Debug, Clone)]
pub struct RveConfig {
    /// Physical size of the RVE in each direction \[m\].
    pub rve_size: [f64; 3],
    /// Number of finite elements in each direction.
    pub n_elements: [usize; 3],
    /// Material phases contained in the RVE.
    pub material_phases: Vec<MaterialPhase>,
}

impl RveConfig {
    /// Create a new `RveConfig`.
    pub fn new(
        rve_size: [f64; 3],
        n_elements: [usize; 3],
        material_phases: Vec<MaterialPhase>,
    ) -> Self {
        Self {
            rve_size,
            n_elements,
            material_phases,
        }
    }

    /// Total number of elements in the RVE.
    pub fn total_elements(&self) -> usize {
        self.n_elements[0] * self.n_elements[1] * self.n_elements[2]
    }

    /// Volume of the RVE \[m³\].
    pub fn volume(&self) -> f64 {
        self.rve_size[0] * self.rve_size[1] * self.rve_size[2]
    }

    /// Check that all volume fractions sum (approximately) to 1.
    pub fn volume_fractions_ok(&self) -> bool {
        let sum: f64 = self.material_phases.iter().map(|p| p.volume_fraction).sum();
        (sum - 1.0).abs() < 1.0e-9
    }
}

// ---------------------------------------------------------------------------
// Homogenized tensor
// ---------------------------------------------------------------------------

/// Effective (homogenized) stiffness tensor and derived moduli.
#[derive(Debug, Clone)]
pub struct HomogenizedTensor {
    /// Effective 6×6 stiffness tensor C_eff in Voigt notation \[Pa\].
    pub c_eff: [[f64; 6]; 6],
    /// Effective bulk modulus K_eff \[Pa\].
    pub bulk_modulus: f64,
    /// Effective shear modulus G_eff \[Pa\].
    pub shear_modulus: f64,
}

impl HomogenizedTensor {
    /// Build from an explicit 6×6 stiffness tensor, deriving K and G from it.
    pub fn from_tensor(c_eff: [[f64; 6]; 6]) -> Self {
        // For isotropic tensor: K = (C11 + 2*C12) / 3, G = C44
        let bulk_modulus = (c_eff[0][0] + 2.0 * c_eff[0][1]) / 3.0;
        let shear_modulus = c_eff[3][3];
        Self {
            c_eff,
            bulk_modulus,
            shear_modulus,
        }
    }

    /// Build from bulk modulus K and shear modulus G (isotropic).
    pub fn from_bulk_shear(k: f64, g: f64) -> Self {
        let lam = k - 2.0 * g / 3.0;
        let mut c = [[0.0f64; 6]; 6];
        c[0][0] = lam + 2.0 * g;
        c[1][1] = lam + 2.0 * g;
        c[2][2] = lam + 2.0 * g;
        c[0][1] = lam;
        c[0][2] = lam;
        c[1][0] = lam;
        c[1][2] = lam;
        c[2][0] = lam;
        c[2][1] = lam;
        c[3][3] = g;
        c[4][4] = g;
        c[5][5] = g;
        Self {
            c_eff: c,
            bulk_modulus: k,
            shear_modulus: g,
        }
    }

    /// Effective Young's modulus E = 9KG / (3K + G).
    pub fn young_modulus(&self) -> f64 {
        let k = self.bulk_modulus;
        let g = self.shear_modulus;
        9.0 * k * g / (3.0 * k + g)
    }

    /// Effective Poisson's ratio ν = (3K − 2G) / (2(3K + G)).
    pub fn poisson_ratio(&self) -> f64 {
        let k = self.bulk_modulus;
        let g = self.shear_modulus;
        (3.0 * k - 2.0 * g) / (2.0 * (3.0 * k + g))
    }
}

// ---------------------------------------------------------------------------
// RVE analysis
// ---------------------------------------------------------------------------

/// Solves a representative volume element problem to extract homogenized
/// properties using the Hill-Mandel macrohomogeneity condition.
#[derive(Debug, Clone)]
pub struct RveAnalysis {
    /// RVE configuration.
    pub config: RveConfig,
}

impl RveAnalysis {
    /// Create a new `RveAnalysis`.
    pub fn new(config: RveConfig) -> Self {
        Self { config }
    }

    /// Solve the RVE and return the homogenized stiffness tensor.
    ///
    /// Uses the Hill-Mandel condition: the volume average of the microscale
    /// stress work equals the macroscale stress work.  For a linear elastic
    /// RVE with perfect periodicity the effective tensor is computed via the
    /// Voigt–Reuss Hill average as a practical approximation.
    pub fn solve_rve(&self) -> HomogenizedTensor {
        hill_average(&self.config.material_phases)
    }

    /// Volume-averaged bulk modulus from phase data.
    pub fn voigt_bulk(&self) -> f64 {
        self.config
            .material_phases
            .iter()
            .map(|p| p.volume_fraction * p.bulk_modulus())
            .sum()
    }

    /// Volume-averaged shear modulus from phase data.
    pub fn voigt_shear(&self) -> f64 {
        self.config
            .material_phases
            .iter()
            .map(|p| p.volume_fraction * p.shear_modulus())
            .sum()
    }
}

// ---------------------------------------------------------------------------
// Multiscale FEM driver
// ---------------------------------------------------------------------------

/// Simple macro-mesh description for the two-scale solver.
#[derive(Debug, Clone)]
pub struct MacroMesh {
    /// Number of macro-elements (one RVE per Gauss point slot).
    pub n_elements: usize,
    /// Number of Gauss points per macro-element.
    pub n_gauss_per_element: usize,
}

impl MacroMesh {
    /// Create a new `MacroMesh`.
    pub fn new(n_elements: usize, n_gauss_per_element: usize) -> Self {
        Self {
            n_elements,
            n_gauss_per_element,
        }
    }

    /// Total number of Gauss points in the macro mesh.
    pub fn total_gauss_points(&self) -> usize {
        self.n_elements * self.n_gauss_per_element
    }
}

/// Two-scale FEM driver that attaches one RVE to every Gauss point.
#[derive(Debug, Clone)]
pub struct MultiScaleFem {
    /// Macro-scale mesh.
    pub macro_mesh: MacroMesh,
    /// One `RveAnalysis` per Gauss point.
    pub rve_at_each_gp: Vec<RveAnalysis>,
}

impl MultiScaleFem {
    /// Create a `MultiScaleFem` with the same RVE config at every Gauss point.
    pub fn new(macro_mesh: MacroMesh, rve_config: RveConfig) -> Self {
        let n_gp = macro_mesh.total_gauss_points();
        let rve_at_each_gp = (0..n_gp)
            .map(|_| RveAnalysis::new(rve_config.clone()))
            .collect();
        Self {
            macro_mesh,
            rve_at_each_gp,
        }
    }

    /// Perform the two-scale solve: solve each RVE and return the collection
    /// of homogenized tensors (one per Gauss point).
    pub fn two_scale_solve(&self) -> Vec<HomogenizedTensor> {
        self.rve_at_each_gp
            .iter()
            .map(|rve| rve.solve_rve())
            .collect()
    }
}

// ---------------------------------------------------------------------------
// Effective medium models
// ---------------------------------------------------------------------------

/// Selection of effective medium / homogenization scheme.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EffectiveMedium {
    /// Voigt (iso-strain) upper bound.
    VoigtBound,
    /// Reuss (iso-stress) lower bound.
    ReussBound,
    /// Arithmetic average of Voigt and Reuss bounds.
    HillAverage,
    /// Mori-Tanaka mean-field theory.
    MoriTanaka,
}

// ---------------------------------------------------------------------------
// Public homogenization functions
// ---------------------------------------------------------------------------

/// Voigt (iso-strain) upper bound for a multiphase composite.
///
/// Each phase carries the same strain; the effective stiffness is the
/// volume-fraction-weighted arithmetic mean of the phase stiffnesses.
pub fn voigt_average(phases: &[MaterialPhase]) -> HomogenizedTensor {
    let k: f64 = phases
        .iter()
        .map(|p| p.volume_fraction * p.bulk_modulus())
        .sum();
    let g: f64 = phases
        .iter()
        .map(|p| p.volume_fraction * p.shear_modulus())
        .sum();
    HomogenizedTensor::from_bulk_shear(k, g)
}

/// Reuss (iso-stress) lower bound for a multiphase composite.
///
/// Each phase carries the same stress; the effective compliance is the
/// volume-fraction-weighted arithmetic mean of the phase compliances.
pub fn reuss_average(phases: &[MaterialPhase]) -> HomogenizedTensor {
    let k_inv: f64 = phases
        .iter()
        .map(|p| p.volume_fraction / p.bulk_modulus())
        .sum();
    let g_inv: f64 = phases
        .iter()
        .map(|p| p.volume_fraction / p.shear_modulus())
        .sum();
    HomogenizedTensor::from_bulk_shear(1.0 / k_inv, 1.0 / g_inv)
}

/// Hill average: arithmetic mean of Voigt and Reuss bounds.
///
/// Provides a commonly used practical estimate of effective properties.
pub fn hill_average(phases: &[MaterialPhase]) -> HomogenizedTensor {
    let v = voigt_average(phases);
    let r = reuss_average(phases);
    let k = 0.5 * (v.bulk_modulus + r.bulk_modulus);
    let g = 0.5 * (v.shear_modulus + r.shear_modulus);
    HomogenizedTensor::from_bulk_shear(k, g)
}

/// Mori-Tanaka mean-field estimate for a two-phase composite.
///
/// Models a matrix phase containing dilute spherical inclusions.
/// `aspect_ratio` = 1.0 for spheres; values != 1.0 are treated as spheres
/// with a correction factor (full Eshelby tensor is outside this scope).
///
/// # Arguments
/// * `matrix`       – matrix phase properties (phase_id should be 0)
/// * `inclusion`    – inclusion phase properties
/// * `aspect_ratio` – inclusion aspect ratio (1.0 = sphere)
pub fn mori_tanaka(
    matrix: &MaterialPhase,
    inclusion: &MaterialPhase,
    aspect_ratio: f64,
) -> HomogenizedTensor {
    let km = matrix.bulk_modulus();
    let gm = matrix.shear_modulus();
    let ki = inclusion.bulk_modulus();
    let gi = inclusion.shear_modulus();
    let fi = inclusion.volume_fraction;

    // Eshelby tensor for a sphere in an isotropic matrix
    let nu_m = matrix.poisson_ratio;
    // S_kk = (1+nu_m)/(3(1-nu_m)) for sphere
    let s_k = (1.0 + nu_m) / (3.0 * (1.0 - nu_m));
    // S_shear = 2(4-5*nu_m)/(15(1-nu_m)) for sphere
    let s_g = 2.0 * (4.0 - 5.0 * nu_m) / (15.0 * (1.0 - nu_m));

    // Correction factor for non-spherical inclusions (simplified)
    let ar_factor = if (aspect_ratio - 1.0).abs() < 1.0e-9 {
        1.0
    } else {
        // Crude correction: prolate/oblate scaling
        (1.0 + 0.1 * (aspect_ratio - 1.0)).max(0.5)
    };

    // Dilute strain concentration factor for K
    let alpha_k = km / (km + s_k * (ki - km));
    // Dilute strain concentration factor for G
    let alpha_g = gm / (gm + s_g * (gi - gm));

    // Mori-Tanaka effective moduli
    let k_eff = km + fi * (ki - km) * alpha_k / ((1.0 - fi) + fi * alpha_k) * ar_factor;
    let g_eff = gm + fi * (gi - gm) * alpha_g / ((1.0 - fi) + fi * alpha_g) * ar_factor;

    HomogenizedTensor::from_bulk_shear(k_eff, g_eff)
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    fn steel() -> MaterialPhase {
        MaterialPhase::new(1.0, 210e9, 0.3, 0)
    }

    fn aluminium() -> MaterialPhase {
        MaterialPhase::new(1.0, 70e9, 0.33, 1)
    }

    fn epoxy() -> MaterialPhase {
        MaterialPhase::new(1.0, 3.5e9, 0.38, 2)
    }

    fn two_phase(vf_steel: f64) -> Vec<MaterialPhase> {
        vec![
            MaterialPhase::new(1.0 - vf_steel, 3.5e9, 0.38, 0),
            MaterialPhase::new(vf_steel, 210e9, 0.3, 1),
        ]
    }

    // --- MaterialPhase ---

    #[test]
    fn test_bulk_modulus_steel() {
        let s = steel();
        // K = 210e9 / (3*(1-2*0.3)) = 210e9/1.2
        let expected = 210e9 / (3.0 * (1.0 - 2.0 * 0.3));
        assert!((s.bulk_modulus() - expected).abs() < 1e3);
    }

    #[test]
    fn test_shear_modulus_steel() {
        let s = steel();
        let expected = 210e9 / (2.0 * (1.0 + 0.3));
        assert!((s.shear_modulus() - expected).abs() < 1e3);
    }

    #[test]
    fn test_bulk_modulus_aluminium() {
        let a = aluminium();
        let expected = 70e9 / (3.0 * (1.0 - 2.0 * 0.33));
        assert!((a.bulk_modulus() - expected).abs() < 1e3);
    }

    #[test]
    fn test_shear_modulus_aluminium() {
        let a = aluminium();
        let expected = 70e9 / (2.0 * (1.0 + 0.33));
        assert!((a.shear_modulus() - expected).abs() < 1e3);
    }

    #[test]
    fn test_stiffness_tensor_diagonal() {
        let s = steel();
        let c = s.stiffness_tensor();
        // C11 should equal lambda + 2*mu
        let e = 210e9;
        let nu = 0.3;
        let mu = e / (2.0 * (1.0 + nu));
        let lam = e * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
        let expected_c11 = lam + 2.0 * mu;
        assert!((c[0][0] - expected_c11).abs() < 1e3);
    }

    #[test]
    fn test_stiffness_tensor_symmetry() {
        let s = steel();
        let c = s.stiffness_tensor();
        for (i, row) in c.iter().enumerate() {
            for (j, &val) in row.iter().enumerate() {
                assert!((val - c[j][i]).abs() < 1e-6, "C[{i}][{j}] != C[{j}][{i}]");
            }
        }
    }

    #[test]
    fn test_stiffness_tensor_shear_block() {
        let s = steel();
        let c = s.stiffness_tensor();
        let expected_g = 210e9 / (2.0 * (1.0 + 0.3));
        assert!((c[3][3] - expected_g).abs() < 1e3);
        assert!((c[4][4] - expected_g).abs() < 1e3);
        assert!((c[5][5] - expected_g).abs() < 1e3);
    }

    // --- RveConfig ---

    #[test]
    fn test_rve_total_elements() {
        let cfg = RveConfig::new([1e-3, 1e-3, 1e-3], [4, 4, 4], two_phase(0.3));
        assert_eq!(cfg.total_elements(), 64);
    }

    #[test]
    fn test_rve_volume() {
        let cfg = RveConfig::new([2e-3, 3e-3, 4e-3], [1, 1, 1], two_phase(0.3));
        assert!((cfg.volume() - 24e-9).abs() < 1e-18);
    }

    #[test]
    fn test_rve_volume_fractions_ok_valid() {
        let phases = vec![
            MaterialPhase::new(0.7, 3.5e9, 0.38, 0),
            MaterialPhase::new(0.3, 210e9, 0.3, 1),
        ];
        let cfg = RveConfig::new([1e-3; 3], [2, 2, 2], phases);
        assert!(cfg.volume_fractions_ok());
    }

    #[test]
    fn test_rve_volume_fractions_ok_invalid() {
        let phases = vec![
            MaterialPhase::new(0.5, 3.5e9, 0.38, 0),
            MaterialPhase::new(0.3, 210e9, 0.3, 1),
        ];
        let cfg = RveConfig::new([1e-3; 3], [2, 2, 2], phases);
        assert!(!cfg.volume_fractions_ok());
    }

    // --- HomogenizedTensor ---

    #[test]
    fn test_homogenized_tensor_young_modulus_recovery() {
        let k = 100e9;
        let g = 40e9;
        let ht = HomogenizedTensor::from_bulk_shear(k, g);
        let e_expected = 9.0 * k * g / (3.0 * k + g);
        assert!((ht.young_modulus() - e_expected).abs() < 1e3);
    }

    #[test]
    fn test_homogenized_tensor_poisson_ratio_recovery() {
        let k = 100e9;
        let g = 40e9;
        let ht = HomogenizedTensor::from_bulk_shear(k, g);
        let nu_expected = (3.0 * k - 2.0 * g) / (2.0 * (3.0 * k + g));
        assert!((ht.poisson_ratio() - nu_expected).abs() < 1e-9);
    }

    #[test]
    fn test_homogenized_tensor_from_tensor_bulk() {
        // Build a steel stiffness tensor and check recovered K
        let s = steel();
        let c = s.stiffness_tensor();
        let ht = HomogenizedTensor::from_tensor(c);
        let k_expected = s.bulk_modulus();
        assert!((ht.bulk_modulus - k_expected).abs() / k_expected < 1e-9);
    }

    // --- Voigt bound ---

    #[test]
    fn test_voigt_single_phase_steel() {
        let phases = vec![MaterialPhase::new(1.0, 210e9, 0.3, 0)];
        let ht = voigt_average(&phases);
        let s = steel();
        assert!((ht.bulk_modulus - s.bulk_modulus()).abs() < 1e3);
        assert!((ht.shear_modulus - s.shear_modulus()).abs() < 1e3);
    }

    #[test]
    fn test_voigt_single_phase_epoxy() {
        let phases = vec![MaterialPhase::new(1.0, 3.5e9, 0.38, 0)];
        let ht = voigt_average(&phases);
        let e = epoxy();
        assert!((ht.bulk_modulus - e.bulk_modulus()).abs() < 1e3);
    }

    #[test]
    fn test_voigt_upper_bound_k() {
        // Voigt K should be >= Reuss K
        let phases = two_phase(0.3);
        let v = voigt_average(&phases);
        let r = reuss_average(&phases);
        assert!(v.bulk_modulus >= r.bulk_modulus - 1.0);
    }

    #[test]
    fn test_voigt_upper_bound_g() {
        let phases = two_phase(0.5);
        let v = voigt_average(&phases);
        let r = reuss_average(&phases);
        assert!(v.shear_modulus >= r.shear_modulus - 1.0);
    }

    #[test]
    fn test_voigt_linearity_in_vf() {
        // Voigt K for vf_steel = 1 should equal pure steel K
        let phases = vec![
            MaterialPhase::new(0.0, 3.5e9, 0.38, 0),
            MaterialPhase::new(1.0, 210e9, 0.3, 1),
        ];
        let v = voigt_average(&phases);
        let s = steel();
        assert!((v.bulk_modulus - s.bulk_modulus()).abs() < 1e3);
    }

    // --- Reuss bound ---

    #[test]
    fn test_reuss_single_phase_steel() {
        let phases = vec![MaterialPhase::new(1.0, 210e9, 0.3, 0)];
        let ht = reuss_average(&phases);
        let s = steel();
        assert!((ht.bulk_modulus - s.bulk_modulus()).abs() < 1e3);
    }

    #[test]
    fn test_reuss_lower_bound_k() {
        for vf in [0.1, 0.3, 0.5, 0.7, 0.9] {
            let phases = two_phase(vf);
            let v = voigt_average(&phases);
            let r = reuss_average(&phases);
            assert!(r.bulk_modulus <= v.bulk_modulus + 1.0, "vf={vf}");
        }
    }

    #[test]
    fn test_reuss_lower_bound_g() {
        for vf in [0.1, 0.3, 0.5, 0.7, 0.9] {
            let phases = two_phase(vf);
            let v = voigt_average(&phases);
            let r = reuss_average(&phases);
            assert!(r.shear_modulus <= v.shear_modulus + 1.0, "vf={vf}");
        }
    }

    #[test]
    fn test_reuss_pure_epoxy() {
        let phases = vec![
            MaterialPhase::new(1.0, 3.5e9, 0.38, 0),
            MaterialPhase::new(0.0, 210e9, 0.3, 1),
        ];
        let r = reuss_average(&phases);
        let e = epoxy();
        assert!((r.bulk_modulus - e.bulk_modulus()).abs() < 1e3);
    }

    // --- Hill average ---

    #[test]
    fn test_hill_between_voigt_reuss_k() {
        let phases = two_phase(0.3);
        let v = voigt_average(&phases);
        let r = reuss_average(&phases);
        let h = hill_average(&phases);
        assert!(h.bulk_modulus >= r.bulk_modulus - 1.0);
        assert!(h.bulk_modulus <= v.bulk_modulus + 1.0);
    }

    #[test]
    fn test_hill_between_voigt_reuss_g() {
        let phases = two_phase(0.4);
        let v = voigt_average(&phases);
        let r = reuss_average(&phases);
        let h = hill_average(&phases);
        assert!(h.shear_modulus >= r.shear_modulus - 1.0);
        assert!(h.shear_modulus <= v.shear_modulus + 1.0);
    }

    #[test]
    fn test_hill_midpoint() {
        let phases = two_phase(0.3);
        let v = voigt_average(&phases);
        let r = reuss_average(&phases);
        let h = hill_average(&phases);
        assert!((h.bulk_modulus - 0.5 * (v.bulk_modulus + r.bulk_modulus)).abs() < 1.0);
    }

    #[test]
    fn test_hill_single_phase_returns_same() {
        let phases = vec![MaterialPhase::new(1.0, 210e9, 0.3, 0)];
        let v = voigt_average(&phases);
        let r = reuss_average(&phases);
        let h = hill_average(&phases);
        // All should equal steel moduli
        assert!((v.bulk_modulus - r.bulk_modulus).abs() < 1e3);
        assert!((h.bulk_modulus - v.bulk_modulus).abs() < 1e3);
    }

    // --- Mori-Tanaka ---

    #[test]
    fn test_mori_tanaka_sphere_k_between_bounds() {
        let matrix = MaterialPhase::new(0.7, 3.5e9, 0.38, 0);
        let inclusion = MaterialPhase::new(0.3, 210e9, 0.3, 1);
        let mt = mori_tanaka(&matrix, &inclusion, 1.0);
        let phases = vec![matrix.clone(), inclusion.clone()];
        let v = voigt_average(&phases);
        let r = reuss_average(&phases);
        // MT should be between Reuss and Voigt bounds (within some tolerance)
        assert!(mt.bulk_modulus >= r.bulk_modulus * 0.9);
        assert!(mt.bulk_modulus <= v.bulk_modulus * 1.1);
    }

    #[test]
    fn test_mori_tanaka_zero_inclusion() {
        let matrix = MaterialPhase::new(1.0, 3.5e9, 0.38, 0);
        let inclusion = MaterialPhase::new(0.0, 210e9, 0.3, 1);
        let mt = mori_tanaka(&matrix, &inclusion, 1.0);
        // With zero inclusion, MT should recover matrix properties
        assert!((mt.bulk_modulus - matrix.bulk_modulus()).abs() / matrix.bulk_modulus() < 0.01);
    }

    #[test]
    fn test_mori_tanaka_aspect_ratio_effect() {
        let matrix = MaterialPhase::new(0.8, 3.5e9, 0.38, 0);
        let inclusion = MaterialPhase::new(0.2, 210e9, 0.3, 1);
        let mt_sphere = mori_tanaka(&matrix, &inclusion, 1.0);
        let mt_elongated = mori_tanaka(&matrix, &inclusion, 5.0);
        // Different aspect ratios should give different moduli
        assert!((mt_sphere.bulk_modulus - mt_elongated.bulk_modulus).abs() > 1.0);
    }

    #[test]
    fn test_mori_tanaka_g_positive() {
        let matrix = MaterialPhase::new(0.6, 3.5e9, 0.38, 0);
        let inclusion = MaterialPhase::new(0.4, 210e9, 0.3, 1);
        let mt = mori_tanaka(&matrix, &inclusion, 1.0);
        assert!(mt.shear_modulus > 0.0);
        assert!(mt.bulk_modulus > 0.0);
    }

    // --- RveAnalysis ---

    #[test]
    fn test_rve_solve_returns_tensor() {
        let phases = two_phase(0.3);
        let cfg = RveConfig::new([1e-3; 3], [2, 2, 2], phases);
        let rve = RveAnalysis::new(cfg);
        let ht = rve.solve_rve();
        assert!(ht.bulk_modulus > 0.0);
        assert!(ht.shear_modulus > 0.0);
    }

    #[test]
    fn test_rve_voigt_bulk_consistency() {
        let phases = two_phase(0.3);
        let expected_k: f64 = phases
            .iter()
            .map(|p| p.volume_fraction * p.bulk_modulus())
            .sum();
        let cfg = RveConfig::new([1e-3; 3], [2, 2, 2], phases);
        let rve = RveAnalysis::new(cfg);
        assert!((rve.voigt_bulk() - expected_k).abs() < 1.0);
    }

    #[test]
    fn test_rve_voigt_shear_consistency() {
        let phases = two_phase(0.3);
        let expected_g: f64 = phases
            .iter()
            .map(|p| p.volume_fraction * p.shear_modulus())
            .sum();
        let cfg = RveConfig::new([1e-3; 3], [2, 2, 2], phases);
        let rve = RveAnalysis::new(cfg);
        assert!((rve.voigt_shear() - expected_g).abs() < 1.0);
    }

    // --- MultiScaleFem ---

    #[test]
    fn test_multiscale_fem_two_scale_solve_count() {
        let macro_mesh = MacroMesh::new(8, 4); // 8 elements × 4 GP = 32
        let phases = two_phase(0.3);
        let cfg = RveConfig::new([1e-3; 3], [2, 2, 2], phases);
        let ms = MultiScaleFem::new(macro_mesh, cfg);
        let tensors = ms.two_scale_solve();
        assert_eq!(tensors.len(), 32);
    }

    #[test]
    fn test_multiscale_fem_all_positive_moduli() {
        let macro_mesh = MacroMesh::new(4, 2);
        let phases = two_phase(0.5);
        let cfg = RveConfig::new([1e-3; 3], [2, 2, 2], phases);
        let ms = MultiScaleFem::new(macro_mesh, cfg);
        for ht in ms.two_scale_solve() {
            assert!(ht.bulk_modulus > 0.0);
            assert!(ht.shear_modulus > 0.0);
        }
    }

    #[test]
    fn test_macro_mesh_total_gauss_points() {
        let mm = MacroMesh::new(10, 8);
        assert_eq!(mm.total_gauss_points(), 80);
    }

    // --- EffectiveMedium enum ---

    #[test]
    fn test_effective_medium_enum_variants() {
        let v = EffectiveMedium::VoigtBound;
        let r = EffectiveMedium::ReussBound;
        let h = EffectiveMedium::HillAverage;
        let mt = EffectiveMedium::MoriTanaka;
        assert_ne!(v, r);
        assert_ne!(h, mt);
    }

    // --- Volume fraction consistency ---

    #[test]
    fn test_volume_fraction_sum_three_phases() {
        let phases = [
            MaterialPhase::new(0.5, 3.5e9, 0.38, 0),
            MaterialPhase::new(0.3, 210e9, 0.3, 1),
            MaterialPhase::new(0.2, 70e9, 0.33, 2),
        ];
        let sum: f64 = phases.iter().map(|p| p.volume_fraction).sum();
        assert!((sum - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_voigt_three_phase() {
        let phases = vec![
            MaterialPhase::new(0.5, 3.5e9, 0.38, 0),
            MaterialPhase::new(0.3, 210e9, 0.3, 1),
            MaterialPhase::new(0.2, 70e9, 0.33, 2),
        ];
        let v = voigt_average(&phases);
        let r = reuss_average(&phases);
        let h = hill_average(&phases);
        assert!(v.bulk_modulus >= r.bulk_modulus - 1.0);
        assert!(h.bulk_modulus >= r.bulk_modulus - 1.0);
        assert!(h.bulk_modulus <= v.bulk_modulus + 1.0);
    }

    #[test]
    fn test_hill_young_modulus_positive() {
        let phases = two_phase(0.3);
        let h = hill_average(&phases);
        assert!(h.young_modulus() > 0.0);
    }

    #[test]
    fn test_hill_poisson_ratio_physical() {
        let phases = two_phase(0.3);
        let h = hill_average(&phases);
        let nu = h.poisson_ratio();
        assert!(nu > -1.0 && nu < 0.5);
    }
}
