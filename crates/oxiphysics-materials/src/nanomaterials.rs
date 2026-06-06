// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Nanomaterial mechanical, thermal, and quantum properties.
//!
//! Covers:
//! - Carbon nanotube (CNT) mechanics for armchair, zigzag, and chiral types
//! - Graphene 2D elastic stiffness tensor
//! - Quantum confinement effects on band gaps and Young's modulus
//! - Nanocomposite rule-of-mixtures (Voigt, Reuss, Halpin-Tsai)
//! - Surface/interface energy effects at the nanoscale
//! - Size-dependent plasticity and Hall-Petch breakdown
//! - Nano-indentation Oliver-Pharr analysis
//! - Molecular mechanics of polymer chains (freely-jointed and worm-like chain)
//! - Nanoparticle surface-area-to-volume ratio
//! - Thermal conductivity phonon mean-free-path model

use rand::RngExt;

// ---------------------------------------------------------------------------
// Constants
// ---------------------------------------------------------------------------

/// Boltzmann constant (J/K).
pub const K_B: f64 = 1.380_649e-23;

/// Planck constant (J·s).
pub const H_PLANCK: f64 = 6.626_070_15e-34;

/// Reduced Planck constant (J·s).
pub const H_BAR: f64 = H_PLANCK / (2.0 * std::f64::consts::PI);

/// Elementary charge (C).
pub const E_CHARGE: f64 = 1.602_176_634e-19;

/// Carbon–carbon bond length in graphene/CNT (nm).
pub const CC_BOND_NM: f64 = 0.142;

/// Graphene lattice parameter (nm).
pub const A_GRAPHENE: f64 = 0.246;

// ---------------------------------------------------------------------------
// CNT chirality
// ---------------------------------------------------------------------------

/// Chirality type of a single-walled carbon nanotube (SWCNT).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CntChirality {
    /// Armchair nanotube: indices (n, n). Always metallic.
    Armchair,
    /// Zigzag nanotube: indices (n, 0). Metallic when n mod 3 == 0.
    Zigzag,
    /// Chiral nanotube: general (n, m) with n ≠ m and m ≠ 0.
    Chiral,
}

impl CntChirality {
    /// Returns whether this chirality type is always metallic.
    pub fn is_metallic(&self) -> bool {
        matches!(self, CntChirality::Armchair)
    }
}

// ---------------------------------------------------------------------------
// CNT geometry
// ---------------------------------------------------------------------------

/// Diameter of an (n, m) SWCNT in nanometres.
pub fn cnt_diameter_nm(n: u32, m: u32) -> f64 {
    let n = n as f64;
    let m = m as f64;
    let a = A_GRAPHENE;
    a * (n * n + n * m + m * m).sqrt() / std::f64::consts::PI
}

/// Chiral angle of an (n, m) SWCNT in radians.
pub fn cnt_chiral_angle(n: u32, m: u32) -> f64 {
    let n = n as f64;
    let m = m as f64;
    (3.0_f64.sqrt() * m / (2.0 * n + m)).atan()
}

/// Returns the chirality classification for (n, m) indices.
pub fn classify_chirality(n: u32, m: u32) -> CntChirality {
    if n == m {
        CntChirality::Armchair
    } else if m == 0 {
        CntChirality::Zigzag
    } else {
        CntChirality::Chiral
    }
}

/// Determines whether an (n, m) CNT is metallic.
///
/// Rule: metallic if `(n - m) mod 3 == 0`.
pub fn cnt_is_metallic(n: u32, m: u32) -> bool {
    let diff = (n as i32 - m as i32).unsigned_abs();
    diff.is_multiple_of(3)
}

// ---------------------------------------------------------------------------
// CNT mechanical properties
// ---------------------------------------------------------------------------

/// Mechanical and thermal properties of a single-walled CNT.
#[derive(Debug, Clone)]
pub struct CntProperties {
    /// Chiral index n.
    pub n: u32,
    /// Chiral index m.
    pub m: u32,
    /// Chirality type.
    pub chirality: CntChirality,
    /// Diameter (nm).
    pub diameter_nm: f64,
    /// Chiral angle (radians).
    pub chiral_angle: f64,
    /// Young's modulus (TPa). Typical value ~1 TPa.
    pub youngs_modulus_tpa: f64,
    /// Poisson's ratio (dimensionless).
    pub poisson_ratio: f64,
    /// Tensile strength (GPa).
    pub tensile_strength_gpa: f64,
    /// Thermal conductivity along tube axis (W/m/K).
    pub thermal_conductivity_axial: f64,
    /// Electronic band gap (eV); 0 for metallic.
    pub band_gap_ev: f64,
}

impl CntProperties {
    /// Constructs CNT properties for a given (n, m) pair using empirical models.
    ///
    /// Young's modulus is modelled as slightly diameter-dependent based on
    /// atomistic simulation fits: `E ≈ 1.06 - 0.02/d²` TPa.
    pub fn new(n: u32, m: u32) -> Self {
        let chirality = classify_chirality(n, m);
        let diameter_nm = cnt_diameter_nm(n, m);
        let chiral_angle = cnt_chiral_angle(n, m);

        // Empirical Young's modulus (TPa) — slight curvature correction
        let d_sq = diameter_nm * diameter_nm;
        let youngs_modulus_tpa = (1.06 - 0.02 / (d_sq + 1e-6)).clamp(0.9, 1.10);

        // Poisson ratio: armchair ~0.16, zigzag ~0.19, chiral interpolated
        let poisson_ratio = match chirality {
            CntChirality::Armchair => 0.16,
            CntChirality::Zigzag => 0.19,
            CntChirality::Chiral => {
                let theta = chiral_angle;
                let frac = theta / (std::f64::consts::PI / 6.0);
                0.16 + 0.03 * frac
            }
        };

        // Tensile strength: ~100 GPa for ideal SWCNTs
        let tensile_strength_gpa = 100.0 - 5.0 / (diameter_nm + 0.5);

        // Thermal conductivity: ballistic at small diameter, diffusive at large
        let thermal_conductivity_axial = 3500.0 * (1.0 - (-diameter_nm).exp());

        // Band gap: metallic => 0, semiconducting => ~0.9 eV/d(nm)
        let band_gap_ev = if cnt_is_metallic(n, m) {
            0.0
        } else {
            0.9 / diameter_nm
        };

        Self {
            n,
            m,
            chirality,
            diameter_nm,
            chiral_angle,
            youngs_modulus_tpa,
            poisson_ratio,
            tensile_strength_gpa,
            thermal_conductivity_axial,
            band_gap_ev,
        }
    }

    /// Returns the wall thickness used in equivalent continuum models (nm).
    ///
    /// Standard value: 0.34 nm (graphene interlayer spacing).
    pub fn wall_thickness_nm(&self) -> f64 {
        0.34
    }

    /// Axial stiffness per unit length `EA` (N) of the equivalent shell.
    pub fn axial_stiffness_n(&self) -> f64 {
        let e_pa = self.youngs_modulus_tpa * 1e12; // Pa
        let h = self.wall_thickness_nm() * 1e-9; // m
        let d = self.diameter_nm * 1e-9; // m
        e_pa * std::f64::consts::PI * d * h
    }

    /// Bending stiffness `EI` (N·m²) of the equivalent hollow cylinder.
    pub fn bending_stiffness_nm2(&self) -> f64 {
        let e_pa = self.youngs_modulus_tpa * 1e12;
        let r = self.diameter_nm * 0.5e-9;
        let h = self.wall_thickness_nm() * 1e-9;
        // I ≈ π r³ h for thin-walled tube
        let i = std::f64::consts::PI * r * r * r * h;
        e_pa * i
    }

    /// Critical Euler buckling load (N) for a tube of length `length_nm`.
    pub fn euler_buckling_load_n(&self, length_nm: f64) -> f64 {
        let ei = self.bending_stiffness_nm2();
        let l = length_nm * 1e-9;
        std::f64::consts::PI * std::f64::consts::PI * ei / (l * l)
    }
}

// ---------------------------------------------------------------------------
// Graphene 2D elastic tensor
// ---------------------------------------------------------------------------

/// In-plane elastic stiffness tensor for monolayer graphene.
///
/// Graphene is a 2D hexagonal crystal with two independent elastic constants:
/// `C₁₁` and `C₁₂`. The tensor is in Voigt notation `[C₁₁, C₁₂, C₆₆]`.
#[derive(Debug, Clone, Copy)]
pub struct GrapheneElasticTensor {
    /// C₁₁ = C₂₂ component (N/m, 2D stiffness).
    pub c11: f64,
    /// C₁₂ component (N/m).
    pub c12: f64,
    /// C₆₆ = (C₁₁ − C₁₂)/2 (N/m), shear component.
    pub c66: f64,
}

impl GrapheneElasticTensor {
    /// Constructs the graphene elastic tensor from ab initio values.
    ///
    /// Reference values: C₁₁ ≈ 357 N/m, C₁₂ ≈ 60 N/m (DFT).
    pub fn default_graphene() -> Self {
        let c11 = 357.0;
        let c12 = 60.0;
        let c66 = (c11 - c12) / 2.0;
        Self { c11, c12, c66 }
    }

    /// Returns the 2D Young's modulus (N/m) in the zigzag direction.
    pub fn youngs_modulus_2d(&self) -> f64 {
        (self.c11 * self.c11 - self.c12 * self.c12) / self.c11
    }

    /// Returns the 2D Young's modulus as bulk value (Pa) given thickness `h` (m).
    pub fn youngs_modulus_3d(&self, thickness_m: f64) -> f64 {
        self.youngs_modulus_2d() / thickness_m
    }

    /// Returns the Poisson ratio `ν = C₁₂/C₁₁`.
    pub fn poisson_ratio(&self) -> f64 {
        self.c12 / self.c11
    }

    /// Returns the biaxial modulus `B = C₁₁ + C₁₂` (N/m).
    pub fn biaxial_modulus(&self) -> f64 {
        self.c11 + self.c12
    }

    /// Computes in-plane stress vector `[σ₁, σ₂, σ₆]` (N/m) from strain `[ε₁, ε₂, γ₆]`.
    pub fn stress_from_strain(&self, strain: [f64; 3]) -> [f64; 3] {
        let s1 = self.c11 * strain[0] + self.c12 * strain[1];
        let s2 = self.c12 * strain[0] + self.c11 * strain[1];
        let s6 = self.c66 * strain[2];
        [s1, s2, s6]
    }

    /// Computes strain from stress using the compliance tensor.
    pub fn strain_from_stress(&self, stress: [f64; 3]) -> [f64; 3] {
        let det = self.c11 * self.c11 - self.c12 * self.c12;
        let e1 = (self.c11 * stress[0] - self.c12 * stress[1]) / det;
        let e2 = (self.c11 * stress[1] - self.c12 * stress[0]) / det;
        let g6 = stress[2] / self.c66;
        [e1, e2, g6]
    }
}

// ---------------------------------------------------------------------------
// Quantum confinement
// ---------------------------------------------------------------------------

/// Quantum confinement model for a spherical quantum dot.
///
/// Uses the particle-in-a-sphere model for the ground-state energy shift.
#[derive(Debug, Clone, Copy)]
pub struct QuantumDot {
    /// Dot radius (nm).
    pub radius_nm: f64,
    /// Effective electron mass (relative to free electron mass).
    pub electron_mass_ratio: f64,
    /// Bulk band gap of the semiconductor (eV).
    pub bulk_bandgap_ev: f64,
}

impl QuantumDot {
    /// Creates a new quantum dot.
    pub fn new(radius_nm: f64, electron_mass_ratio: f64, bulk_bandgap_ev: f64) -> Self {
        Self {
            radius_nm,
            electron_mass_ratio,
            bulk_bandgap_ev,
        }
    }

    /// Returns the quantum confinement energy shift ΔE (eV).
    ///
    /// Uses the Brus equation: `ΔE = ħ²π²/(2m*r²)`.
    pub fn confinement_energy_ev(&self) -> f64 {
        let m_star = self.electron_mass_ratio * 9.109_383_7e-31; // kg
        let r = self.radius_nm * 1e-9; // m
        let pi = std::f64::consts::PI;
        H_BAR * H_BAR * pi * pi / (2.0 * m_star * r * r * E_CHARGE) // eV
    }

    /// Returns the effective band gap of the quantum dot (eV).
    pub fn effective_bandgap_ev(&self) -> f64 {
        self.bulk_bandgap_ev + self.confinement_energy_ev()
    }

    /// Returns the emission wavelength (nm) corresponding to the quantum dot band gap.
    pub fn emission_wavelength_nm(&self) -> f64 {
        let e_j = self.effective_bandgap_ev() * E_CHARGE; // J
        H_PLANCK * 3e8 / e_j * 1e9 // nm
    }

    /// Estimates the size-dependent Young's modulus relative to bulk.
    ///
    /// Uses a surface-stress model: `E(r) = E_bulk * (1 + α/r)` where `α` ~ 0.1 nm.
    pub fn size_dependent_youngs_ratio(&self, alpha_nm: f64) -> f64 {
        1.0 + alpha_nm / (self.radius_nm + 1e-15)
    }
}

// ---------------------------------------------------------------------------
// Nanocomposite rule-of-mixtures
// ---------------------------------------------------------------------------

/// Nanocomposite mechanical model using rule of mixtures and Halpin-Tsai.
#[derive(Debug, Clone)]
pub struct Nanocomposite {
    /// Matrix Young's modulus (GPa).
    pub e_matrix_gpa: f64,
    /// Filler Young's modulus (GPa).
    pub e_filler_gpa: f64,
    /// Matrix Poisson ratio.
    pub nu_matrix: f64,
    /// Filler Poisson ratio.
    pub nu_filler: f64,
    /// Volume fraction of filler (0–1).
    pub volume_fraction: f64,
    /// Filler aspect ratio (length/diameter).
    pub aspect_ratio: f64,
}

impl Nanocomposite {
    /// Creates a new nanocomposite.
    pub fn new(
        e_matrix_gpa: f64,
        e_filler_gpa: f64,
        nu_matrix: f64,
        nu_filler: f64,
        volume_fraction: f64,
        aspect_ratio: f64,
    ) -> Self {
        Self {
            e_matrix_gpa,
            e_filler_gpa,
            nu_matrix,
            nu_filler,
            volume_fraction,
            aspect_ratio,
        }
    }

    /// Voigt upper bound Young's modulus (GPa).
    ///
    /// `E_Voigt = Vf * Ef + Vm * Em`
    pub fn voigt_modulus_gpa(&self) -> f64 {
        let vf = self.volume_fraction;
        let vm = 1.0 - vf;
        vf * self.e_filler_gpa + vm * self.e_matrix_gpa
    }

    /// Reuss lower bound Young's modulus (GPa).
    ///
    /// `1/E_Reuss = Vf/Ef + Vm/Em`
    pub fn reuss_modulus_gpa(&self) -> f64 {
        let vf = self.volume_fraction;
        let vm = 1.0 - vf;
        1.0 / (vf / (self.e_filler_gpa + 1e-15) + vm / (self.e_matrix_gpa + 1e-15))
    }

    /// Halpin-Tsai longitudinal modulus (GPa).
    ///
    /// Uses parameter ξ = 2 × aspect_ratio for aligned fibers.
    pub fn halpin_tsai_longitudinal_gpa(&self) -> f64 {
        let xi = 2.0 * self.aspect_ratio;
        let ef = self.e_filler_gpa;
        let em = self.e_matrix_gpa;
        let eta = (ef / em - 1.0) / (ef / em + xi);
        let vf = self.volume_fraction;
        em * (1.0 + xi * eta * vf) / (1.0 - eta * vf)
    }

    /// Halpin-Tsai transverse modulus (GPa).
    ///
    /// Uses ξ = 2 for transverse stiffness of discontinuous fibers.
    pub fn halpin_tsai_transverse_gpa(&self) -> f64 {
        let xi = 2.0;
        let ef = self.e_filler_gpa;
        let em = self.e_matrix_gpa;
        let eta = (ef / em - 1.0) / (ef / em + xi);
        let vf = self.volume_fraction;
        em * (1.0 + xi * eta * vf) / (1.0 - eta * vf)
    }

    /// Effective Poisson ratio via rule of mixtures.
    pub fn effective_poisson_ratio(&self) -> f64 {
        let vf = self.volume_fraction;
        vf * self.nu_filler + (1.0 - vf) * self.nu_matrix
    }

    /// Bounds check: returns true if Voigt ≥ Reuss (as required by theory).
    pub fn bounds_valid(&self) -> bool {
        self.voigt_modulus_gpa() >= self.reuss_modulus_gpa() - 1e-9
    }
}

// ---------------------------------------------------------------------------
// Surface/interface energy at nanoscale
// ---------------------------------------------------------------------------

/// Surface energy model for nanoparticles.
///
/// Accounts for the increasing importance of surface energy
/// relative to bulk energy as particle radius decreases.
#[derive(Debug, Clone, Copy)]
pub struct NanoparticleSurface {
    /// Particle radius (nm).
    pub radius_nm: f64,
    /// Bulk surface energy density (J/m²).
    pub surface_energy_j_m2: f64,
    /// Bulk cohesive energy per unit volume (J/m³).
    pub cohesive_energy_j_m3: f64,
}

impl NanoparticleSurface {
    /// Creates a new nanoparticle surface model.
    pub fn new(radius_nm: f64, surface_energy_j_m2: f64, cohesive_energy_j_m3: f64) -> Self {
        Self {
            radius_nm,
            surface_energy_j_m2,
            cohesive_energy_j_m3,
        }
    }

    /// Surface area (m²) of a spherical nanoparticle.
    pub fn surface_area_m2(&self) -> f64 {
        let r = self.radius_nm * 1e-9;
        4.0 * std::f64::consts::PI * r * r
    }

    /// Volume (m³) of a spherical nanoparticle.
    pub fn volume_m3(&self) -> f64 {
        let r = self.radius_nm * 1e-9;
        4.0 / 3.0 * std::f64::consts::PI * r * r * r
    }

    /// Surface-area-to-volume ratio (1/m).
    pub fn sa_to_volume_ratio(&self) -> f64 {
        self.surface_area_m2() / self.volume_m3()
    }

    /// Total surface energy (J).
    pub fn total_surface_energy_j(&self) -> f64 {
        self.surface_energy_j_m2 * self.surface_area_m2()
    }

    /// Ratio of surface energy to bulk cohesive energy (dimensionless).
    ///
    /// This ratio → 1 for very small particles, indicating surface dominance.
    pub fn surface_energy_ratio(&self) -> f64 {
        let e_surf = self.total_surface_energy_j();
        let e_bulk = self.cohesive_energy_j_m3 * self.volume_m3();
        e_surf / (e_bulk + 1e-30)
    }

    /// Laplace pressure inside the nanoparticle (Pa).
    ///
    /// `ΔP = 2γ/r` (Young-Laplace).
    pub fn laplace_pressure_pa(&self) -> f64 {
        let r = self.radius_nm * 1e-9;
        2.0 * self.surface_energy_j_m2 / (r + 1e-30)
    }

    /// Size-dependent melting point suppression factor (Gibbs-Thomson effect).
    ///
    /// `ΔTm/Tm_bulk ≈ −4γ_sl/(ρ L d)` — returns the ratio `ΔTm/Tm_bulk`.
    pub fn melting_point_suppression_ratio(
        &self,
        density_kg_m3: f64,
        latent_heat_j_kg: f64,
    ) -> f64 {
        let d = 2.0 * self.radius_nm * 1e-9;
        -4.0 * self.surface_energy_j_m2 / (density_kg_m3 * latent_heat_j_kg * d + 1e-30)
    }
}

// ---------------------------------------------------------------------------
// Hall-Petch and its breakdown
// ---------------------------------------------------------------------------

/// Hall-Petch and inverse Hall-Petch (breakdown at nanoscale).
///
/// Describes grain-size-dependent yield strength.
#[derive(Debug, Clone, Copy)]
pub struct HallPetchModel {
    /// Friction stress (yield strength of a single crystal) σ₀ (MPa).
    pub sigma_0_mpa: f64,
    /// Hall-Petch slope k (MPa·√nm).
    pub k_hp_mpa_sqrt_nm: f64,
    /// Critical grain size below which inverse HP (softening) occurs (nm).
    pub d_critical_nm: f64,
}

impl HallPetchModel {
    /// Creates a Hall-Petch model.
    pub fn new(sigma_0_mpa: f64, k_hp_mpa_sqrt_nm: f64, d_critical_nm: f64) -> Self {
        Self {
            sigma_0_mpa,
            k_hp_mpa_sqrt_nm,
            d_critical_nm,
        }
    }

    /// Yield strength (MPa) using the classical Hall-Petch equation.
    ///
    /// `σy = σ₀ + k / √d`
    pub fn yield_strength_classical_mpa(&self, grain_size_nm: f64) -> f64 {
        self.sigma_0_mpa + self.k_hp_mpa_sqrt_nm / grain_size_nm.sqrt()
    }

    /// Yield strength (MPa) with inverse Hall-Petch softening for d < d_critical.
    ///
    /// Uses a bilinear model: classical HP for d ≥ d_c, then linear softening.
    pub fn yield_strength_mpa(&self, grain_size_nm: f64) -> f64 {
        if grain_size_nm >= self.d_critical_nm {
            self.yield_strength_classical_mpa(grain_size_nm)
        } else {
            // Inverse HP: softening proportional to (d_c - d)
            let sigma_at_dc = self.yield_strength_classical_mpa(self.d_critical_nm);
            let slope = sigma_at_dc / (self.d_critical_nm * 2.0);
            (sigma_at_dc - slope * (self.d_critical_nm - grain_size_nm)).max(0.0)
        }
    }

    /// Grain size (nm) at which yield strength is maximised.
    pub fn optimal_grain_size_nm(&self) -> f64 {
        self.d_critical_nm
    }
}

// ---------------------------------------------------------------------------
// Oliver-Pharr nano-indentation
// ---------------------------------------------------------------------------

/// Oliver-Pharr method for nano-indentation analysis.
///
/// Extracts hardness and reduced modulus from load-displacement data.
#[derive(Debug, Clone)]
pub struct OliverPharr {
    /// Peak load (mN).
    pub p_max_mn: f64,
    /// Contact depth at peak load (nm).
    pub h_max_nm: f64,
    /// Contact stiffness S = dP/dh at unloading (mN/nm).
    pub stiffness_mn_per_nm: f64,
    /// Indenter geometry constant β (1.034 for Berkovich).
    pub beta: f64,
    /// Indenter area function coefficient C₀ (nm²/nm²).
    ///
    /// The projected contact area is `A = C₀ * h_c²`.
    pub c0: f64,
}

impl OliverPharr {
    /// Creates an Oliver-Pharr analysis instance.
    pub fn new(p_max_mn: f64, h_max_nm: f64, stiffness_mn_per_nm: f64, beta: f64, c0: f64) -> Self {
        Self {
            p_max_mn,
            h_max_nm,
            stiffness_mn_per_nm,
            beta,
            c0,
        }
    }

    /// Returns the default Berkovich indenter instance placeholder.
    ///
    /// The C₀ for a perfect Berkovich tip is 24.5 (θ = 65.27°).
    pub fn berkovich() -> Self {
        Self {
            p_max_mn: 1.0,
            h_max_nm: 100.0,
            stiffness_mn_per_nm: 0.05,
            beta: 1.034,
            c0: 24.5,
        }
    }

    /// Contact depth `h_c` (nm).
    ///
    /// `h_c = h_max - ε * P_max / S` where `ε = 0.75` for paraboloid.
    pub fn contact_depth_nm(&self) -> f64 {
        let eps = 0.75;
        self.h_max_nm - eps * self.p_max_mn / (self.stiffness_mn_per_nm + 1e-30)
    }

    /// Projected contact area `A` (nm²).
    pub fn contact_area_nm2(&self) -> f64 {
        let hc = self.contact_depth_nm();
        self.c0 * hc * hc
    }

    /// Hardness H (GPa).
    ///
    /// `H = P_max / A`
    pub fn hardness_gpa(&self) -> f64 {
        let a_m2 = self.contact_area_nm2() * 1e-18; // nm² → m²
        let p_n = self.p_max_mn * 1e-3; // mN → N
        p_n / (a_m2 + 1e-30) / 1e9 // Pa → GPa
    }

    /// Reduced modulus E_r (GPa).
    ///
    /// `E_r = (√π / (2β)) * S / √A`
    pub fn reduced_modulus_gpa(&self) -> f64 {
        let a_nm2 = self.contact_area_nm2();
        let sqrt_a_nm = a_nm2.sqrt() * 1e-9; // m
        let s_n_per_m = self.stiffness_mn_per_nm * 1e-3 / 1e-9; // N/m
        let pi = std::f64::consts::PI;
        let er_pa = (pi.sqrt() / (2.0 * self.beta)) * s_n_per_m / (sqrt_a_nm + 1e-30);
        er_pa / 1e9
    }

    /// Extracts the sample Young's modulus (GPa) given indenter properties.
    ///
    /// `1/E_r = (1-ν²)/E + (1-νi²)/Ei`
    pub fn sample_youngs_modulus_gpa(
        &self,
        poisson_sample: f64,
        indenter_modulus_gpa: f64,
        indenter_poisson: f64,
    ) -> f64 {
        let er = self.reduced_modulus_gpa();
        let inv_er = 1.0 / (er + 1e-15);
        let inv_ei = (1.0 - indenter_poisson * indenter_poisson) / (indenter_modulus_gpa + 1e-15);
        let inv_e = inv_er - inv_ei;
        if inv_e < 1e-15 {
            return 0.0;
        }
        (1.0 - poisson_sample * poisson_sample) / inv_e
    }
}

// ---------------------------------------------------------------------------
// Polymer chain molecular mechanics
// ---------------------------------------------------------------------------

/// Freely-jointed chain (FJC) model for a linear polymer.
#[derive(Debug, Clone, Copy)]
pub struct FreelyJointedChain {
    /// Number of statistical segments.
    pub n_segments: u64,
    /// Kuhn segment length (nm).
    pub b_nm: f64,
}

impl FreelyJointedChain {
    /// Creates a new FJC.
    pub fn new(n_segments: u64, b_nm: f64) -> Self {
        Self { n_segments, b_nm }
    }

    /// End-to-end distance RMS `<R²>^(1/2)` (nm).
    pub fn rms_end_to_end_nm(&self) -> f64 {
        (self.n_segments as f64).sqrt() * self.b_nm
    }

    /// Contour length `L` (nm).
    pub fn contour_length_nm(&self) -> f64 {
        self.n_segments as f64 * self.b_nm
    }

    /// Radius of gyration `Rg` (nm).
    pub fn radius_of_gyration_nm(&self) -> f64 {
        self.rms_end_to_end_nm() / 6.0_f64.sqrt()
    }

    /// Entropic restoring force (pN) at extension `r_nm`.
    ///
    /// Uses the Langevin approximation: `f = (k_B T / b) L^{-1}(r/L)`
    /// where `L^{-1}` is the inverse Langevin, approximated by `3x/(1-x²)`.
    pub fn restoring_force_pn(&self, r_nm: f64, temp_k: f64) -> f64 {
        let l = self.contour_length_nm() * 1e-9;
        let r = (r_nm * 1e-9).min(l * 0.9999);
        let x = r / l;
        // Cohen's Padé approximation to inverse Langevin
        let inv_l = x * (3.0 - x * x) / (1.0 - x * x);
        let b = self.b_nm * 1e-9;
        let f_n = K_B * temp_k / b * inv_l;
        f_n * 1e12 // pN
    }
}

/// Worm-like chain (WLC) model for semi-flexible polymers (e.g. DNA).
#[derive(Debug, Clone, Copy)]
pub struct WormLikeChain {
    /// Persistence length `Lp` (nm).
    pub persistence_length_nm: f64,
    /// Contour length `L` (nm).
    pub contour_length_nm: f64,
}

impl WormLikeChain {
    /// Creates a new WLC.
    pub fn new(persistence_length_nm: f64, contour_length_nm: f64) -> Self {
        Self {
            persistence_length_nm,
            contour_length_nm,
        }
    }

    /// End-to-end distance RMS (nm) using the WLC result.
    pub fn rms_end_to_end_nm(&self) -> f64 {
        let lp = self.persistence_length_nm;
        let l = self.contour_length_nm;
        (2.0 * lp * l * (1.0 - lp / l * (1.0 - (-l / lp).exp()))).sqrt()
    }

    /// Restoring force (pN) at extension `r_nm` (Marko-Siggia interpolation).
    ///
    /// `f ≈ (k_B T / Lp) [1/(4(1-x)²) − 1/4 + x]` where `x = r/L`.
    pub fn restoring_force_pn(&self, r_nm: f64, temp_k: f64) -> f64 {
        let l = self.contour_length_nm * 1e-9;
        let lp = self.persistence_length_nm * 1e-9;
        let r = (r_nm * 1e-9).min(l * 0.9999);
        let x = r / l;
        let factor = K_B * temp_k / lp;
        let f_n = factor * (0.25 / ((1.0 - x) * (1.0 - x)) - 0.25 + x);
        f_n * 1e12 // pN
    }

    /// Bending stiffness `κ = k_B T * Lp` (J·m).
    pub fn bending_stiffness_jm(&self, temp_k: f64) -> f64 {
        K_B * temp_k * self.persistence_length_nm * 1e-9
    }
}

// ---------------------------------------------------------------------------
// Nanoparticle surface-area-to-volume
// ---------------------------------------------------------------------------

/// Shape of a nanoparticle for SA/V computation.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum NanoparticleShape {
    /// Spherical nanoparticle.
    Sphere,
    /// Cubic nanoparticle.
    Cube,
    /// Cylindrical nanoparticle (rod).
    Cylinder,
    /// Disk-shaped nanoparticle.
    Disk,
}

/// Surface-area-to-volume ratio (nm⁻¹) for a nanoparticle.
///
/// - Sphere: `6/d` (d = diameter)
/// - Cube: `6/a` (a = edge)
/// - Cylinder: `2(r+l)/(rl)` with aspect ratio `l/d`
/// - Disk: `2(r+h)/(rh)`
pub fn sa_to_volume_ratio_nm(
    shape: NanoparticleShape,
    characteristic_dim_nm: f64,
    aspect_ratio: f64,
) -> f64 {
    match shape {
        NanoparticleShape::Sphere => {
            let d = characteristic_dim_nm;
            6.0 / d
        }
        NanoparticleShape::Cube => {
            let a = characteristic_dim_nm;
            6.0 / a
        }
        NanoparticleShape::Cylinder => {
            let r = characteristic_dim_nm / 2.0;
            let l = r * 2.0 * aspect_ratio;
            (2.0 * std::f64::consts::PI * r * (r + l)) / (std::f64::consts::PI * r * r * l)
        }
        NanoparticleShape::Disk => {
            let r = characteristic_dim_nm / 2.0;
            let h = characteristic_dim_nm / aspect_ratio;
            2.0 * (r + h) / (r * h)
        }
    }
}

/// Fraction of atoms on the surface of a spherical nanoparticle.
///
/// Approximation: `f_surf ≈ 6 * d_atom / d_particle`.
pub fn surface_atom_fraction(particle_diameter_nm: f64, atom_diameter_nm: f64) -> f64 {
    (6.0 * atom_diameter_nm / particle_diameter_nm).min(1.0)
}

// ---------------------------------------------------------------------------
// Nanoscale thermal conductivity — phonon mean free path
// ---------------------------------------------------------------------------

/// Phonon-limited thermal conductivity model at the nanoscale.
///
/// Based on the kinetic theory expression `κ = (1/3) C_v v Λ`
/// with a Matthiessen-rule combination of bulk and size-limited MFP.
#[derive(Debug, Clone, Copy)]
pub struct PhononThermalModel {
    /// Bulk thermal conductivity (W/m/K).
    pub kappa_bulk_w_mk: f64,
    /// Phonon mean free path in bulk (nm).
    pub mfp_bulk_nm: f64,
    /// Speed of sound (group velocity) (m/s).
    pub v_sound_m_s: f64,
    /// Volumetric heat capacity (J/m³/K).
    pub cv_j_m3_k: f64,
}

impl PhononThermalModel {
    /// Creates a phonon thermal model.
    pub fn new(kappa_bulk_w_mk: f64, mfp_bulk_nm: f64, v_sound_m_s: f64, cv_j_m3_k: f64) -> Self {
        Self {
            kappa_bulk_w_mk,
            mfp_bulk_nm,
            v_sound_m_s,
            cv_j_m3_k,
        }
    }

    /// Silicon-like defaults at 300 K.
    pub fn silicon_300k() -> Self {
        Self {
            kappa_bulk_w_mk: 150.0,
            mfp_bulk_nm: 300.0,
            v_sound_m_s: 8_430.0,
            cv_j_m3_k: 1.63e6,
        }
    }

    /// Effective phonon MFP (nm) in a nanostructure of characteristic size `d_nm`.
    ///
    /// Matthiessen's rule: `1/Λ = 1/Λ_bulk + 1/d`.
    pub fn effective_mfp_nm(&self, d_nm: f64) -> f64 {
        1.0 / (1.0 / self.mfp_bulk_nm + 1.0 / d_nm)
    }

    /// Effective thermal conductivity (W/m/K) for a nanostructure of size `d_nm`.
    ///
    /// Scales as `κ_eff = κ_bulk * Λ_eff / Λ_bulk`.
    pub fn effective_conductivity_w_mk(&self, d_nm: f64) -> f64 {
        let lam_eff = self.effective_mfp_nm(d_nm);
        self.kappa_bulk_w_mk * lam_eff / self.mfp_bulk_nm
    }

    /// Knudsen number `Kn = Λ/d` (dimensionless).
    ///
    /// `Kn >> 1`: ballistic regime; `Kn << 1`: diffusive.
    pub fn knudsen_number(&self, d_nm: f64) -> f64 {
        self.mfp_bulk_nm / d_nm
    }

    /// Thermal boundary (Kapitza) resistance (m²·K/W) at a grain boundary,
    /// using the acoustic mismatch model approximation.
    ///
    /// `R_K ≈ 4 / (C_v v)`.
    pub fn kapitza_resistance_m2kw(&self) -> f64 {
        4.0 / (self.cv_j_m3_k * self.v_sound_m_s + 1e-30)
    }
}

// ---------------------------------------------------------------------------
// Size-dependent elastic modulus (combined model)
// ---------------------------------------------------------------------------

/// Computes a size-dependent Young's modulus (GPa) using the core-shell model.
///
/// `E_nano = E_bulk * (1 - f_surf) + E_surf * f_surf`
///
/// where `f_surf` is the surface atom fraction.
pub fn size_dependent_modulus_gpa(
    e_bulk_gpa: f64,
    e_surface_gpa: f64,
    particle_diameter_nm: f64,
    atom_diameter_nm: f64,
) -> f64 {
    let f = surface_atom_fraction(particle_diameter_nm, atom_diameter_nm);
    e_bulk_gpa * (1.0 - f) + e_surface_gpa * f
}

// ---------------------------------------------------------------------------
// Multi-walled CNT (MWCNT) model
// ---------------------------------------------------------------------------

/// Properties of a multi-walled CNT approximated as concentric SWCNTs.
#[derive(Debug, Clone)]
pub struct MwcntProperties {
    /// Number of walls.
    pub n_walls: u32,
    /// Outermost diameter (nm).
    pub outer_diameter_nm: f64,
    /// Interlayer spacing (nm), typically 0.34 nm.
    pub interlayer_spacing_nm: f64,
    /// Effective Young's modulus (TPa) (lower than SWCNT due to interlayer shear).
    pub youngs_modulus_tpa: f64,
}

impl MwcntProperties {
    /// Creates a MWCNT model from outer diameter and wall count.
    pub fn new(outer_diameter_nm: f64, n_walls: u32) -> Self {
        let interlayer_spacing_nm = 0.34;
        // Effective E decreases due to interlayer coupling
        let youngs_modulus_tpa = 0.9 - 0.05 * (n_walls as f64 - 1.0).max(0.0) / 10.0;
        Self {
            n_walls,
            outer_diameter_nm,
            interlayer_spacing_nm,
            youngs_modulus_tpa: youngs_modulus_tpa.max(0.5),
        }
    }

    /// Inner diameter (nm).
    pub fn inner_diameter_nm(&self) -> f64 {
        (self.outer_diameter_nm - 2.0 * self.n_walls as f64 * self.interlayer_spacing_nm).max(0.5)
    }

    /// Cross-sectional area of tube wall (nm²).
    pub fn wall_area_nm2(&self) -> f64 {
        let pi = std::f64::consts::PI;
        let ro = self.outer_diameter_nm / 2.0;
        let ri = self.inner_diameter_nm() / 2.0;
        pi * (ro * ro - ri * ri)
    }

    /// Axial stiffness EA (nN) = E * A.
    pub fn axial_stiffness_nn(&self) -> f64 {
        let e_pa = self.youngs_modulus_tpa * 1e12;
        let a_m2 = self.wall_area_nm2() * 1e-18;
        e_pa * a_m2 * 1e9 // nN
    }
}

// ---------------------------------------------------------------------------
// Monte-Carlo sampling of CNT bundles
// ---------------------------------------------------------------------------

/// Result of Monte-Carlo sampling of nanotube orientation effects.
#[derive(Debug, Clone)]
pub struct CntBundleSample {
    /// Number of CNTs sampled.
    pub n_cnt: usize,
    /// Mean axial Young's modulus (TPa).
    pub mean_modulus_tpa: f64,
    /// Standard deviation of modulus (TPa).
    pub std_modulus_tpa: f64,
    /// Mean diameter (nm).
    pub mean_diameter_nm: f64,
}

/// Samples a bundle of CNTs with random (n, m) indices and returns aggregate properties.
pub fn sample_cnt_bundle(n_samples: usize, n_max: u32) -> CntBundleSample {
    let mut rng = rand::rng();
    let mut moduli = Vec::with_capacity(n_samples);
    let mut diams = Vec::with_capacity(n_samples);
    for _ in 0..n_samples {
        let n = rng.random_range(5..=n_max);
        let m = rng.random_range(0..=n);
        let cnt = CntProperties::new(n, m);
        moduli.push(cnt.youngs_modulus_tpa);
        diams.push(cnt.diameter_nm);
    }
    let n = moduli.len() as f64;
    let mean_mod = moduli.iter().sum::<f64>() / n;
    let std_mod = {
        let var = moduli
            .iter()
            .map(|&x| (x - mean_mod) * (x - mean_mod))
            .sum::<f64>()
            / n;
        var.sqrt()
    };
    let mean_diam = diams.iter().sum::<f64>() / n;
    CntBundleSample {
        n_cnt: n_samples,
        mean_modulus_tpa: mean_mod,
        std_modulus_tpa: std_mod,
        mean_diameter_nm: mean_diam,
    }
}

// ---------------------------------------------------------------------------
// Graphene grain boundary model
// ---------------------------------------------------------------------------

/// Poly-crystalline graphene model with grain boundary effects.
#[derive(Debug, Clone, Copy)]
pub struct PolycrystallineGraphene {
    /// Average grain size (nm).
    pub grain_size_nm: f64,
    /// Grain boundary energy density (J/m²).
    pub grain_boundary_energy_j_m2: f64,
    /// Reduction factor of Young's modulus due to grain boundaries.
    pub modulus_reduction_factor: f64,
}

impl PolycrystallineGraphene {
    /// Creates a polycrystalline graphene model.
    pub fn new(grain_size_nm: f64, grain_boundary_energy_j_m2: f64) -> Self {
        // Larger grains → less reduction
        let modulus_reduction_factor = 1.0 - 0.05 * (10.0 / grain_size_nm).min(1.0);
        Self {
            grain_size_nm,
            grain_boundary_energy_j_m2,
            modulus_reduction_factor,
        }
    }

    /// Effective Young's modulus (2D, N/m) including grain boundary reduction.
    pub fn effective_modulus_2d_nm(&self) -> f64 {
        let pristine = GrapheneElasticTensor::default_graphene();
        pristine.youngs_modulus_2d() * self.modulus_reduction_factor
    }

    /// Grain boundary area fraction (per unit area of graphene).
    ///
    /// Approximated as `4 * δ / d_g` where δ = 0.5 nm boundary width.
    pub fn grain_boundary_area_fraction(&self) -> f64 {
        let delta_nm = 0.5;
        (4.0 * delta_nm / self.grain_size_nm).min(1.0)
    }
}

// ---------------------------------------------------------------------------
// Summary / Report
// ---------------------------------------------------------------------------

/// A compact summary of nanomaterial properties for reporting.
#[derive(Debug, Clone)]
pub struct NanomaterialSummary {
    /// Name of the material.
    pub name: String,
    /// Young's modulus (GPa).
    pub youngs_modulus_gpa: f64,
    /// Tensile strength (GPa).
    pub tensile_strength_gpa: f64,
    /// Thermal conductivity (W/m/K).
    pub thermal_conductivity_w_mk: f64,
    /// Characteristic dimension (nm).
    pub characteristic_dim_nm: f64,
}

impl NanomaterialSummary {
    /// Creates a summary from a SWCNT.
    pub fn from_cnt(cnt: &CntProperties) -> Self {
        Self {
            name: format!("SWCNT({},{})", cnt.n, cnt.m),
            youngs_modulus_gpa: cnt.youngs_modulus_tpa * 1000.0,
            tensile_strength_gpa: cnt.tensile_strength_gpa,
            thermal_conductivity_w_mk: cnt.thermal_conductivity_axial,
            characteristic_dim_nm: cnt.diameter_nm,
        }
    }
}

// ---------------------------------------------------------------------------
// Unit tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    const PI: f64 = std::f64::consts::PI;

    // --- CNT geometry ---

    #[test]
    fn test_cnt_diameter_armchair_10_10() {
        let d = cnt_diameter_nm(10, 10);
        // Expected: a * sqrt(300) / pi ≈ 0.246 * 17.32 / 3.14159 ≈ 1.356 nm
        assert!((d - 1.356).abs() < 0.01, "Diameter: {d}");
    }

    #[test]
    fn test_cnt_chirality_armchair() {
        assert_eq!(classify_chirality(5, 5), CntChirality::Armchair);
        assert_eq!(classify_chirality(10, 10), CntChirality::Armchair);
    }

    #[test]
    fn test_cnt_chirality_zigzag() {
        assert_eq!(classify_chirality(8, 0), CntChirality::Zigzag);
    }

    #[test]
    fn test_cnt_chirality_chiral() {
        assert_eq!(classify_chirality(6, 4), CntChirality::Chiral);
    }

    #[test]
    fn test_cnt_is_metallic_armchair() {
        // (n,n): n-n=0 → metallic
        assert!(cnt_is_metallic(5, 5));
        assert!(cnt_is_metallic(10, 10));
    }

    #[test]
    fn test_cnt_is_metallic_rule() {
        // (6,3): 6-3=3 → metallic
        assert!(cnt_is_metallic(6, 3));
        // (7,3): 7-3=4 → not metallic
        assert!(!cnt_is_metallic(7, 3));
    }

    #[test]
    fn test_cnt_chiral_angle_armchair() {
        let theta = cnt_chiral_angle(5, 5);
        // Armchair angle = 30° = pi/6
        assert!((theta - PI / 6.0).abs() < 1e-10);
    }

    #[test]
    fn test_cnt_chiral_angle_zigzag() {
        let theta = cnt_chiral_angle(8, 0);
        assert!(theta.abs() < 1e-12);
    }

    // --- CntProperties ---

    #[test]
    fn test_cnt_properties_modulus_range() {
        let cnt = CntProperties::new(10, 10);
        assert!(cnt.youngs_modulus_tpa > 0.8 && cnt.youngs_modulus_tpa < 1.2);
    }

    #[test]
    fn test_cnt_properties_band_gap_metallic() {
        let cnt = CntProperties::new(5, 5);
        assert!(cnt.band_gap_ev.abs() < 1e-12);
    }

    #[test]
    fn test_cnt_properties_band_gap_semiconducting() {
        let cnt = CntProperties::new(7, 3);
        assert!(cnt.band_gap_ev > 0.0);
    }

    #[test]
    fn test_cnt_axial_stiffness_positive() {
        let cnt = CntProperties::new(10, 10);
        assert!(cnt.axial_stiffness_n() > 0.0);
    }

    #[test]
    fn test_cnt_bending_stiffness_positive() {
        let cnt = CntProperties::new(10, 10);
        assert!(cnt.bending_stiffness_nm2() > 0.0);
    }

    #[test]
    fn test_cnt_euler_buckling_decreases_with_length() {
        let cnt = CntProperties::new(10, 10);
        let f1 = cnt.euler_buckling_load_n(10.0);
        let f2 = cnt.euler_buckling_load_n(20.0);
        assert!(f1 > f2);
    }

    // --- Graphene ---

    #[test]
    fn test_graphene_modulus_positive() {
        let g = GrapheneElasticTensor::default_graphene();
        assert!(g.youngs_modulus_2d() > 0.0);
    }

    #[test]
    fn test_graphene_poisson_ratio() {
        let g = GrapheneElasticTensor::default_graphene();
        let nu = g.poisson_ratio();
        assert!(nu > 0.0 && nu < 0.5);
    }

    #[test]
    fn test_graphene_stress_strain_roundtrip() {
        let g = GrapheneElasticTensor::default_graphene();
        let strain = [1e-3, 0.5e-3, 0.2e-3];
        let stress = g.stress_from_strain(strain);
        let recovered = g.strain_from_stress(stress);
        for i in 0..3 {
            assert!((recovered[i] - strain[i]).abs() < 1e-15, "Mismatch at {i}");
        }
    }

    #[test]
    fn test_graphene_biaxial_modulus() {
        let g = GrapheneElasticTensor::default_graphene();
        assert!(g.biaxial_modulus() > g.c11);
    }

    // --- Quantum dot ---

    #[test]
    fn test_quantum_dot_confinement_positive() {
        let qd = QuantumDot::new(3.0, 0.067, 1.42);
        assert!(qd.confinement_energy_ev() > 0.0);
    }

    #[test]
    fn test_quantum_dot_bandgap_exceeds_bulk() {
        let qd = QuantumDot::new(3.0, 0.067, 1.42);
        assert!(qd.effective_bandgap_ev() > qd.bulk_bandgap_ev);
    }

    #[test]
    fn test_quantum_dot_wavelength_positive() {
        let qd = QuantumDot::new(3.0, 0.067, 1.42);
        assert!(qd.emission_wavelength_nm() > 0.0);
    }

    #[test]
    fn test_quantum_dot_size_modulus_larger() {
        let qd = QuantumDot::new(2.0, 0.067, 1.42);
        let ratio = qd.size_dependent_youngs_ratio(0.1);
        assert!(ratio > 1.0);
    }

    // --- Nanocomposite ---

    #[test]
    fn test_nanocomposite_voigt_bounds() {
        let nc = Nanocomposite::new(3.0, 1000.0, 0.35, 0.2, 0.05, 100.0);
        assert!(nc.voigt_modulus_gpa() > nc.reuss_modulus_gpa());
    }

    #[test]
    fn test_nanocomposite_bounds_valid() {
        let nc = Nanocomposite::new(3.0, 1000.0, 0.35, 0.2, 0.05, 100.0);
        assert!(nc.bounds_valid());
    }

    #[test]
    fn test_nanocomposite_halpin_tsai_between_bounds() {
        let nc = Nanocomposite::new(3.0, 1000.0, 0.35, 0.2, 0.05, 100.0);
        let ht = nc.halpin_tsai_longitudinal_gpa();
        assert!(ht > nc.reuss_modulus_gpa());
        assert!(ht < nc.voigt_modulus_gpa() + 1.0);
    }

    #[test]
    fn test_nanocomposite_zero_vf() {
        let nc = Nanocomposite::new(3.0, 1000.0, 0.35, 0.2, 0.0, 100.0);
        assert!((nc.voigt_modulus_gpa() - 3.0).abs() < 1e-10);
    }

    // --- Surface energy ---

    #[test]
    fn test_surface_sa_v_ratio_increases_small_r() {
        let p1 = NanoparticleSurface::new(10.0, 1.0, 1e9);
        let p2 = NanoparticleSurface::new(1.0, 1.0, 1e9);
        assert!(p2.sa_to_volume_ratio() > p1.sa_to_volume_ratio());
    }

    #[test]
    fn test_laplace_pressure_positive() {
        let p = NanoparticleSurface::new(5.0, 1.5, 1e9);
        assert!(p.laplace_pressure_pa() > 0.0);
    }

    #[test]
    fn test_surface_energy_ratio_small_particle() {
        let p = NanoparticleSurface::new(0.5, 1.0, 1e9);
        // Should be large fraction
        assert!(p.surface_energy_ratio() > 0.01);
    }

    // --- Hall-Petch ---

    #[test]
    fn test_hall_petch_classical_increases_fine_grain() {
        let hp = HallPetchModel::new(100.0, 500.0, 10.0);
        let sy_100 = hp.yield_strength_classical_mpa(100.0);
        let sy_10 = hp.yield_strength_classical_mpa(10.0);
        assert!(sy_10 > sy_100);
    }

    #[test]
    fn test_hall_petch_breakdown_softens() {
        let hp = HallPetchModel::new(100.0, 500.0, 20.0);
        let sy_dc = hp.yield_strength_mpa(20.0);
        let sy_below = hp.yield_strength_mpa(5.0);
        assert!(sy_below < sy_dc);
    }

    #[test]
    fn test_hall_petch_optimal_at_critical() {
        let hp = HallPetchModel::new(100.0, 500.0, 20.0);
        assert!((hp.optimal_grain_size_nm() - 20.0).abs() < 1e-12);
    }

    // --- Oliver-Pharr ---

    #[test]
    fn test_oliver_pharr_contact_depth_less_than_hmax() {
        let op = OliverPharr::new(2.0, 150.0, 0.08, 1.034, 24.5);
        assert!(op.contact_depth_nm() < op.h_max_nm);
    }

    #[test]
    fn test_oliver_pharr_hardness_positive() {
        let op = OliverPharr::new(2.0, 150.0, 0.08, 1.034, 24.5);
        assert!(op.hardness_gpa() > 0.0);
    }

    #[test]
    fn test_oliver_pharr_reduced_modulus_positive() {
        let op = OliverPharr::new(2.0, 150.0, 0.08, 1.034, 24.5);
        assert!(op.reduced_modulus_gpa() > 0.0);
    }

    #[test]
    fn test_oliver_pharr_sample_modulus_positive() {
        let op = OliverPharr::new(2.0, 150.0, 0.08, 1.034, 24.5);
        let e = op.sample_youngs_modulus_gpa(0.25, 1140.0, 0.07);
        assert!(e > 0.0, "Sample modulus: {e}");
    }

    // --- Polymer chains ---

    #[test]
    fn test_fjc_rms_end_to_end() {
        let fjc = FreelyJointedChain::new(100, 1.0);
        let r = fjc.rms_end_to_end_nm();
        assert!((r - 10.0).abs() < 1e-10); // sqrt(100) * 1.0
    }

    #[test]
    fn test_fjc_contour_length() {
        let fjc = FreelyJointedChain::new(50, 2.0);
        assert!((fjc.contour_length_nm() - 100.0).abs() < 1e-10);
    }

    #[test]
    fn test_fjc_restoring_force_positive() {
        let fjc = FreelyJointedChain::new(1000, 0.5);
        let f = fjc.restoring_force_pn(100.0, 300.0);
        assert!(f > 0.0);
    }

    #[test]
    fn test_wlc_rms_end_to_end_positive() {
        let wlc = WormLikeChain::new(50.0, 1000.0);
        assert!(wlc.rms_end_to_end_nm() > 0.0);
    }

    #[test]
    fn test_wlc_restoring_force_positive() {
        let wlc = WormLikeChain::new(50.0, 1000.0);
        let f = wlc.restoring_force_pn(200.0, 300.0);
        assert!(f > 0.0);
    }

    #[test]
    fn test_wlc_bending_stiffness_positive() {
        let wlc = WormLikeChain::new(50.0, 1000.0);
        assert!(wlc.bending_stiffness_jm(300.0) > 0.0);
    }

    // --- SA/V ratio ---

    #[test]
    fn test_sa_to_v_sphere() {
        // 6/d for d=10 nm → 0.6 nm^{-1}
        let r = sa_to_volume_ratio_nm(NanoparticleShape::Sphere, 10.0, 1.0);
        assert!((r - 0.6).abs() < 1e-12);
    }

    #[test]
    fn test_sa_to_v_cube() {
        let r = sa_to_volume_ratio_nm(NanoparticleShape::Cube, 10.0, 1.0);
        assert!((r - 0.6).abs() < 1e-12);
    }

    #[test]
    fn test_surface_atom_fraction_clamped() {
        let f = surface_atom_fraction(0.1, 0.3); // tiny particle
        assert!(f <= 1.0);
    }

    // --- Phonon thermal model ---

    #[test]
    fn test_phonon_effective_conductivity_less_than_bulk() {
        let m = PhononThermalModel::silicon_300k();
        let k_eff = m.effective_conductivity_w_mk(30.0);
        assert!(k_eff < m.kappa_bulk_w_mk);
    }

    #[test]
    fn test_phonon_knudsen_large_for_small_structure() {
        let m = PhononThermalModel::silicon_300k();
        let kn = m.knudsen_number(1.0); // 1 nm structure
        assert!(kn > 1.0);
    }

    #[test]
    fn test_kapitza_resistance_positive() {
        let m = PhononThermalModel::silicon_300k();
        assert!(m.kapitza_resistance_m2kw() > 0.0);
    }

    // --- Size-dependent modulus ---

    #[test]
    fn test_size_dependent_modulus_bulk_limit() {
        // When particle is huge, f_surf → 0, so E_nano → E_bulk
        let e = size_dependent_modulus_gpa(200.0, 50.0, 1000.0, 0.25);
        assert!((e - 200.0).abs() < 1.0);
    }

    // --- MWCNT ---

    #[test]
    fn test_mwcnt_inner_diameter_positive() {
        let m = MwcntProperties::new(10.0, 3);
        assert!(m.inner_diameter_nm() > 0.0);
    }

    #[test]
    fn test_mwcnt_axial_stiffness_positive() {
        let m = MwcntProperties::new(10.0, 3);
        assert!(m.axial_stiffness_nn() > 0.0);
    }

    // --- Bundle sampling ---

    #[test]
    fn test_cnt_bundle_sample_count() {
        let s = sample_cnt_bundle(50, 15);
        assert_eq!(s.n_cnt, 50);
    }

    #[test]
    fn test_cnt_bundle_sample_modulus_range() {
        let s = sample_cnt_bundle(200, 20);
        assert!(s.mean_modulus_tpa > 0.8 && s.mean_modulus_tpa < 1.2);
    }

    // --- Polycrystalline graphene ---

    #[test]
    fn test_polycrystalline_graphene_modulus_less_than_pristine() {
        let pg = PolycrystallineGraphene::new(10.0, 0.5);
        let pristine = GrapheneElasticTensor::default_graphene().youngs_modulus_2d();
        assert!(pg.effective_modulus_2d_nm() < pristine + 1.0);
    }

    #[test]
    fn test_gb_area_fraction_increases_small_grain() {
        let pg1 = PolycrystallineGraphene::new(100.0, 0.5);
        let pg2 = PolycrystallineGraphene::new(5.0, 0.5);
        assert!(pg2.grain_boundary_area_fraction() > pg1.grain_boundary_area_fraction());
    }

    // --- Summary ---

    #[test]
    fn test_nanomaterial_summary_from_cnt() {
        let cnt = CntProperties::new(10, 10);
        let summary = NanomaterialSummary::from_cnt(&cnt);
        assert!(summary.youngs_modulus_gpa > 0.0);
        assert!(!summary.name.is_empty());
    }
}
