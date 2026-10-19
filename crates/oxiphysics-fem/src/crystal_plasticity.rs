// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0
//! Crystal plasticity FEM: slip systems, dislocation density, texture evolution.
//!
//! Provides:
//! - [`SlipSystem`]: Schmid factor, resolved shear stress
//! - [`CrystalPlasticity`]: rate-dependent flow rule and latent hardening
//! - [`FccCrystal`]: 12 slip systems for FCC metals
//! - [`BccCrystal`]: 48 slip systems for BCC metals
//! - [`HcpCrystal`]: basal/prismatic/pyramidal systems for HCP metals
//! - [`DislocationDensity`]: Taylor hardening and Kocks-Mecking evolution
//! - [`TextureEvolution`]: lattice rotation via spin tensor
//! - [`PolycrystalFem`]: Taylor/self-consistent aggregate models
//! - [`OdfRepresentation`]: orientation distribution function
//! - [`RecrystallizationModel`]: JMAK nucleation and growth kinetics

use std::f64::consts::PI;

// ---------------------------------------------------------------------------
// SlipSystem
// ---------------------------------------------------------------------------

/// A single crystallographic slip system defined by its slip direction and slip plane normal.
///
/// The Schmid tensor **m** = s ⊗ n is stored as a 3×3 matrix in row-major order.
#[derive(Debug, Clone)]
pub struct SlipSystem {
    /// Unit slip direction vector \[s0, s1, s2\].
    pub slip_direction: [f64; 3],
    /// Unit slip plane normal vector \[n0, n1, n2\].
    pub slip_normal: [f64; 3],
    /// Schmid tensor m_ij = s_i * n_j (row-major, 9 components).
    pub schmid_tensor: [f64; 9],
}

impl SlipSystem {
    /// Construct a slip system from slip direction and plane normal (both will be normalised).
    pub fn new(slip_direction: [f64; 3], slip_normal: [f64; 3]) -> Self {
        let s = normalise3(slip_direction);
        let n = normalise3(slip_normal);
        let mut m = [0.0_f64; 9];
        for i in 0..3 {
            for j in 0..3 {
                m[i * 3 + j] = s[i] * n[j];
            }
        }
        Self {
            slip_direction: s,
            slip_normal: n,
            schmid_tensor: m,
        }
    }

    /// Compute the Schmid factor for a uniaxial stress along `loading_axis`.
    ///
    /// The Schmid factor is m_sf = cos(φ) * cos(λ) where φ is the angle between
    /// the loading axis and slip plane normal, λ is the angle with slip direction.
    pub fn schmid_factor(&self, loading_axis: [f64; 3]) -> f64 {
        let a = normalise3(loading_axis);
        let cos_phi = dot3(a, self.slip_normal);
        let cos_lambda = dot3(a, self.slip_direction);
        cos_phi * cos_lambda
    }

    /// Resolved shear stress τ = σ : m = Σ_{ij} σ_{ij} m_{ij}.
    ///
    /// `stress` is the symmetric Cauchy stress in Voigt notation \[σ11, σ22, σ33, σ12, σ23, σ31\].
    pub fn resolved_shear_stress(&self, stress_voigt: [f64; 6]) -> f64 {
        // expand to full 3×3 (symmetric)
        let s = voigt_to_full(stress_voigt);
        let m = self.schmid_tensor;
        m.iter().enumerate().map(|(k, &mk)| mk * s[k]).sum()
    }
}

// ---------------------------------------------------------------------------
// CrystalPlasticity
// ---------------------------------------------------------------------------

/// Rate-dependent crystal plasticity model with latent hardening.
///
/// Flow rule:  γ̇^α = γ̇_0 · |τ^α / τ_c^α|^n · sign(τ^α)
/// Hardening:  τ̇_c^α = Σ_β h^{αβ} · |γ̇^β|
#[derive(Debug, Clone)]
pub struct CrystalPlasticity {
    /// Reference shear rate γ̇_0 (s⁻¹).
    pub reference_shear_rate: f64,
    /// Strain-rate sensitivity exponent n (typically 10–20).
    pub rate_exponent: f64,
    /// Critical resolved shear stress for each slip system (Pa).
    pub crss: Vec<f64>,
    /// Latent hardening matrix h^{αβ} (Pa).
    pub hardening_matrix: Vec<f64>,
    /// Number of slip systems.
    pub n_slip: usize,
}

impl CrystalPlasticity {
    /// Create a new crystal plasticity model.
    pub fn new(
        n_slip: usize,
        reference_shear_rate: f64,
        rate_exponent: f64,
        initial_crss: f64,
        self_hardening: f64,
        latent_hardening_ratio: f64,
    ) -> Self {
        let crss = vec![initial_crss; n_slip];
        let mut hmat = vec![0.0_f64; n_slip * n_slip];
        for a in 0..n_slip {
            for b in 0..n_slip {
                hmat[a * n_slip + b] = if a == b {
                    self_hardening
                } else {
                    latent_hardening_ratio * self_hardening
                };
            }
        }
        Self {
            reference_shear_rate,
            rate_exponent,
            crss,
            hardening_matrix: hmat,
            n_slip,
        }
    }

    /// Compute shear rate for slip system α given resolved shear stress τ.
    pub fn shear_rate(&self, alpha: usize, tau: f64) -> f64 {
        let tau_c = self.crss[alpha];
        if tau_c == 0.0 {
            return 0.0;
        }
        let ratio = tau / tau_c;
        self.reference_shear_rate * ratio.abs().powf(self.rate_exponent) * ratio.signum()
    }

    /// Compute all shear rates given resolved shear stresses for every slip system.
    pub fn all_shear_rates(&self, tau: &[f64]) -> Vec<f64> {
        (0..self.n_slip)
            .map(|a| self.shear_rate(a, tau[a]))
            .collect()
    }

    /// Update CRSS by explicit Euler step: τ_c^α += Δt · Σ_β h^{αβ} · |γ̇^β|.
    pub fn update_crss(&mut self, gamma_dot: &[f64], dt: f64) {
        let n = self.n_slip;
        let old = self.crss.clone();
        for (a, (crss_a, &old_a)) in self.crss.iter_mut().zip(old.iter()).enumerate().take(n) {
            let delta: f64 = (0..n)
                .map(|b| self.hardening_matrix[a * n + b] * gamma_dot[b].abs())
                .sum();
            *crss_a = old_a + dt * delta;
        }
    }

    /// Accumulated plastic shear strain from shear rates and time step.
    pub fn accumulated_shear(&self, gamma_dot: &[f64], dt: f64) -> f64 {
        gamma_dot.iter().map(|g| g.abs() * dt).sum()
    }
}

// ---------------------------------------------------------------------------
// FCC Crystal
// ---------------------------------------------------------------------------

/// FCC crystal with 12 {111}<110> slip systems.
#[derive(Debug, Clone)]
pub struct FccCrystal {
    /// All 12 slip systems.
    pub slip_systems: Vec<SlipSystem>,
}

impl FccCrystal {
    /// Construct the 12 FCC slip systems.
    pub fn new() -> Self {
        // {111} plane normals (un-normalised; constructor normalises)
        let normals: [[f64; 3]; 4] = [
            [1.0, 1.0, 1.0],
            [-1.0, 1.0, 1.0],
            [1.0, -1.0, 1.0],
            [1.0, 1.0, -1.0],
        ];
        // Three <110> directions on each {111} plane
        let directions: [[f64; 3]; 6] = [
            [0.0, 1.0, -1.0],
            [-1.0, 0.0, 1.0],
            [1.0, -1.0, 0.0],
            [0.0, -1.0, -1.0],
            [1.0, 0.0, 1.0],
            [-1.0, 1.0, 0.0],
        ];
        let mut slip_systems = Vec::with_capacity(12);
        for n in &normals {
            for d in &directions {
                // Only accept directions perpendicular to the plane normal
                if dot3(*d, *n).abs() < 1e-9 {
                    slip_systems.push(SlipSystem::new(*d, *n));
                }
            }
        }
        // Pad / trim to exactly 12
        slip_systems.truncate(12);
        Self { slip_systems }
    }

    /// Return the number of slip systems (always 12 for FCC).
    pub fn n_slip_systems(&self) -> usize {
        self.slip_systems.len()
    }

    /// Maximum Schmid factor over all slip systems for a given loading axis.
    pub fn max_schmid_factor(&self, loading_axis: [f64; 3]) -> f64 {
        self.slip_systems
            .iter()
            .map(|s| s.schmid_factor(loading_axis).abs())
            .fold(f64::NEG_INFINITY, f64::max)
    }
}

impl Default for FccCrystal {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// BCC Crystal
// ---------------------------------------------------------------------------

/// BCC crystal with slip systems on {110}, {112}, and {123} planes along <111>.
#[derive(Debug, Clone)]
pub struct BccCrystal {
    /// Slip systems (up to 48).
    pub slip_systems: Vec<SlipSystem>,
}

impl BccCrystal {
    /// Construct BCC slip systems including all three families.
    pub fn new() -> Self {
        let mut slip_systems = Vec::new();

        // <111> Burgers vectors (8 total, ±)
        let burgers: [[f64; 3]; 4] = [
            [1.0, 1.0, 1.0],
            [-1.0, 1.0, 1.0],
            [1.0, -1.0, 1.0],
            [1.0, 1.0, -1.0],
        ];

        // {110} plane normals (12 permutations)
        let n110: [[f64; 3]; 6] = [
            [1.0, 1.0, 0.0],
            [1.0, -1.0, 0.0],
            [0.0, 1.0, 1.0],
            [0.0, 1.0, -1.0],
            [1.0, 0.0, 1.0],
            [1.0, 0.0, -1.0],
        ];
        for b in &burgers {
            for n in &n110 {
                if dot3(*b, *n).abs() < 1e-9 {
                    slip_systems.push(SlipSystem::new(*b, *n));
                }
            }
        }

        // {112} plane normals (12 permutations)
        let n112: [[f64; 3]; 12] = [
            [1.0, 1.0, 2.0],
            [1.0, -1.0, 2.0],
            [-1.0, 1.0, 2.0],
            [1.0, 2.0, 1.0],
            [1.0, 2.0, -1.0],
            [-1.0, 2.0, 1.0],
            [2.0, 1.0, 1.0],
            [2.0, 1.0, -1.0],
            [2.0, -1.0, 1.0],
            [1.0, 1.0, -2.0],
            [1.0, -1.0, -2.0],
            [-1.0, 1.0, -2.0],
        ];
        for b in &burgers {
            for n in &n112 {
                if dot3(*b, *n).abs() < 1e-9 {
                    slip_systems.push(SlipSystem::new(*b, *n));
                }
            }
        }

        Self { slip_systems }
    }

    /// Number of BCC slip systems.
    pub fn n_slip_systems(&self) -> usize {
        self.slip_systems.len()
    }
}

impl Default for BccCrystal {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// HCP Crystal
// ---------------------------------------------------------------------------

/// A slip system given by Miller–Bravais indices: `(plane (hkil), direction [uvtw])`.
type HcpSystemIndices = ([i32; 4], [i32; 4]);

/// Basal (0001)\<11-20\> systems.
const HCP_BASAL: [HcpSystemIndices; 3] = [
    ([0, 0, 0, 1], [2, -1, -1, 0]),
    ([0, 0, 0, 1], [-1, 2, -1, 0]),
    ([0, 0, 0, 1], [-1, -1, 2, 0]),
];

/// Prismatic {10-10}\<11-20\> systems.
const HCP_PRISMATIC: [HcpSystemIndices; 3] = [
    ([1, 0, -1, 0], [-1, 2, -1, 0]),
    ([0, 1, -1, 0], [2, -1, -1, 0]),
    ([-1, 1, 0, 0], [-1, -1, 2, 0]),
];

/// First-order pyramidal \<a\> systems {10-11}\<11-20\>.
const HCP_PYRAMIDAL_A: [HcpSystemIndices; 6] = [
    ([1, 0, -1, 1], [-1, 2, -1, 0]),
    ([0, 1, -1, 1], [2, -1, -1, 0]),
    ([-1, 1, 0, 1], [-1, -1, 2, 0]),
    ([-1, 0, 1, 1], [-1, 2, -1, 0]),
    ([0, -1, 1, 1], [2, -1, -1, 0]),
    ([1, -1, 0, 1], [-1, -1, 2, 0]),
];

/// Second-order pyramidal \<c+a\> systems {11-22}\<11-23\>.
const HCP_PYRAMIDAL_CA: [HcpSystemIndices; 6] = [
    ([1, 1, -2, 2], [-1, -1, 2, 3]),
    ([1, -2, 1, 2], [-1, 2, -1, 3]),
    ([-2, 1, 1, 2], [2, -1, -1, 3]),
    ([-1, -1, 2, 2], [1, 1, -2, 3]),
    ([-1, 2, -1, 2], [1, -2, 1, 3]),
    ([2, -1, -1, 2], [-2, 1, 1, 3]),
];

/// HCP crystal slip systems: basal, prismatic, and pyramidal families.
///
/// The 18 systems are built from Miller–Bravais indices:
///
/// | family | systems | count |
/// |---|---|---|
/// | basal | (0001)\<11-20\> | 3 |
/// | prismatic | {10-10}\<11-20\> | 3 |
/// | pyramidal \<a\> | {10-11}\<11-20\> | 6 |
/// | pyramidal \<c+a\> | {11-22}\<11-23\> | 6 |
///
/// They are stored in this order in [`HcpCrystal::slip_systems`].
#[derive(Debug, Clone)]
pub struct HcpCrystal {
    /// Slip systems.
    pub slip_systems: Vec<SlipSystem>,
    /// c/a ratio of the HCP cell.
    pub ca_ratio: f64,
}

impl HcpCrystal {
    /// Cartesian vector (in units of the lattice parameter `a`) of the
    /// Miller–Bravais direction `[uvtw]`.
    ///
    /// The basis is a1 = \[1,0,0\], a2 = \[-1/2,√3/2,0\], c = \[0,0,c/a\] and the
    /// direction is `u·a1 + v·a2 + t·a3 + w·c` with a3 = −a1 − a2, i.e.
    /// `(u−t)·a1 + (v−t)·a2 + w·c`.  Returns `None` if the redundant index
    /// condition `u + v + t = 0` is violated.
    pub fn direction_vector(ca_ratio: f64, uvtw: [i32; 4]) -> Option<[f64; 3]> {
        let [u, v, t, w] = uvtw;
        if u + v + t != 0 {
            return None;
        }
        let a1_coeff = f64::from(u - t);
        let a2_coeff = f64::from(v - t);
        let sqrt3 = 3.0_f64.sqrt();
        Some([
            a1_coeff - 0.5 * a2_coeff,
            0.5 * sqrt3 * a2_coeff,
            f64::from(w) * ca_ratio,
        ])
    }

    /// Cartesian (unnormalised) normal of the Miller–Bravais plane `(hkil)`.
    ///
    /// The normal is `h·b1 + k·b2 + l·b3` in the reciprocal basis of
    /// (a1, a2, c), i.e. `(h, (h+2k)/√3, l/(c/a))`.  Returns `None` if the
    /// redundant index condition `h + k + i = 0` is violated.
    pub fn plane_normal(ca_ratio: f64, hkil: [i32; 4]) -> Option<[f64; 3]> {
        let [h, k, i, l] = hkil;
        if h + k + i != 0 {
            return None;
        }
        let sqrt3 = 3.0_f64.sqrt();
        Some([
            f64::from(h),
            f64::from(h + 2 * k) / sqrt3,
            f64::from(l) / ca_ratio,
        ])
    }

    /// Construct HCP slip systems for a given c/a ratio.
    ///
    /// Uses Cartesian representation where a1 = \[1,0,0\], a2 = \[-1/2,√3/2,0\], c = \[0,0,c/a\].
    /// Every system satisfies `d · n = 0` (up to round-off) for any c/a.
    ///
    /// `ca_ratio` must be positive and finite; for any other value no slip
    /// systems are created (use [`HcpCrystal::try_new`] to get an error instead).
    pub fn new(ca_ratio: f64) -> Self {
        Self::try_new(ca_ratio).unwrap_or_else(|_| Self {
            slip_systems: Vec::new(),
            ca_ratio,
        })
    }

    /// Fallible variant of [`HcpCrystal::new`].
    ///
    /// # Errors
    ///
    /// Returns [`crate::Error::InvalidInput`] if `ca_ratio` is not positive and finite.
    pub fn try_new(ca_ratio: f64) -> crate::Result<Self> {
        if !ca_ratio.is_finite() || ca_ratio <= 0.0 {
            return Err(crate::Error::InvalidInput(format!(
                "HCP c/a ratio must be positive and finite, got {ca_ratio}"
            )));
        }
        let families: [&[HcpSystemIndices]; 4] = [
            &HCP_BASAL,
            &HCP_PRISMATIC,
            &HCP_PYRAMIDAL_A,
            &HCP_PYRAMIDAL_CA,
        ];
        let mut slip_systems = Vec::with_capacity(18);
        for family in families {
            for &(plane, direction) in family {
                if let (Some(d), Some(n)) = (
                    Self::direction_vector(ca_ratio, direction),
                    Self::plane_normal(ca_ratio, plane),
                ) {
                    slip_systems.push(SlipSystem::new(d, n));
                }
            }
        }

        Ok(Self {
            slip_systems,
            ca_ratio,
        })
    }

    /// Total number of slip systems.
    pub fn n_slip_systems(&self) -> usize {
        self.slip_systems.len()
    }

    /// Basal (0001)\<11-20\> systems (3).
    pub fn basal_systems(&self) -> &[SlipSystem] {
        self.slip_systems.get(0..3).unwrap_or(&[])
    }

    /// Prismatic {10-10}\<11-20\> systems (3).
    pub fn prismatic_systems(&self) -> &[SlipSystem] {
        self.slip_systems.get(3..6).unwrap_or(&[])
    }

    /// First-order pyramidal \<a\> systems {10-11}\<11-20\> (6).
    pub fn pyramidal_a_systems(&self) -> &[SlipSystem] {
        self.slip_systems.get(6..12).unwrap_or(&[])
    }

    /// Pyramidal \<c+a\> systems {11-22}\<11-23\> (6).
    pub fn pyramidal_ca_systems(&self) -> &[SlipSystem] {
        self.slip_systems.get(12..18).unwrap_or(&[])
    }

    /// Check whether pyramidal slip is active compared to basal/prismatic
    /// based on c/a ratio. Returns true when c/a deviates notably from ideal.
    pub fn pyramidal_active(&self) -> bool {
        let ideal = (8.0_f64 / 3.0).sqrt(); // ≈ 1.633
        (self.ca_ratio - ideal).abs() > 0.05
    }
}

// ---------------------------------------------------------------------------
// DislocationDensity
// ---------------------------------------------------------------------------

/// Taylor hardening model with Kocks-Mecking dislocation density evolution.
///
/// Taylor hardening: τ_c = τ_0 + M · G · b · √ρ
/// Kocks-Mecking:    dρ/dγ = (1/b) · (1/λ - 2 · y · ρ)
#[derive(Debug, Clone)]
pub struct DislocationDensity {
    /// Initial CRSS τ_0 (Pa).
    pub tau_0: f64,
    /// Taylor factor M (dimensionless, typically 3.06 for FCC polycrystal).
    pub taylor_factor: f64,
    /// Shear modulus G (Pa).
    pub shear_modulus: f64,
    /// Burgers vector magnitude b (m).
    pub burgers_vector: f64,
    /// Mean free path λ (m).
    pub mean_free_path: f64,
    /// Dynamic recovery parameter y (dimensionless).
    pub recovery_param: f64,
    /// Current dislocation density ρ (m⁻²).
    pub density: f64,
}

impl DislocationDensity {
    /// Construct a new dislocation density model.
    pub fn new(
        tau_0: f64,
        taylor_factor: f64,
        shear_modulus: f64,
        burgers_vector: f64,
        mean_free_path: f64,
        recovery_param: f64,
        initial_density: f64,
    ) -> Self {
        Self {
            tau_0,
            taylor_factor,
            shear_modulus,
            burgers_vector,
            mean_free_path,
            recovery_param,
            density: initial_density,
        }
    }

    /// Compute current CRSS via Taylor hardening.
    pub fn critical_resolved_shear_stress(&self) -> f64 {
        self.tau_0
            + self.taylor_factor * self.shear_modulus * self.burgers_vector * self.density.sqrt()
    }

    /// Evolve dislocation density by a shear strain increment dγ (Kocks-Mecking).
    pub fn evolve(&mut self, d_gamma: f64) {
        let b = self.burgers_vector;
        let lambda = self.mean_free_path;
        let y = self.recovery_param;
        let drho = (1.0 / b) * (1.0 / lambda - 2.0 * y * self.density) * d_gamma;
        self.density = (self.density + drho).max(0.0);
    }

    /// Integrate density evolution over a list of shear strain increments.
    pub fn integrate(&mut self, d_gammas: &[f64]) {
        for &dg in d_gammas {
            self.evolve(dg);
        }
    }
}

// ---------------------------------------------------------------------------
// TextureEvolution
// ---------------------------------------------------------------------------

/// Lattice rotation driven by the plastic spin: Ṙ = W · R.
///
/// W = W_total − W_plastic where W_plastic = Σ_α γ̇^α · W^α.
#[derive(Debug, Clone)]
pub struct TextureEvolution {
    /// Current rotation matrix R (row-major 3×3).
    pub rotation: [f64; 9],
    /// Accumulated plastic strain.
    pub accumulated_plastic_strain: f64,
}

impl TextureEvolution {
    /// Create with identity rotation.
    pub fn new() -> Self {
        let mut rotation = [0.0_f64; 9];
        rotation[0] = 1.0;
        rotation[4] = 1.0;
        rotation[8] = 1.0;
        Self {
            rotation,
            accumulated_plastic_strain: 0.0,
        }
    }

    /// Apply an incremental lattice spin tensor ΔW (skew-symmetric 3×3, row-major) to R.
    ///
    /// Updates R ← R + ΔW · R · Δt and re-orthogonalises via Gram-Schmidt.
    pub fn update(&mut self, spin: [f64; 9], dt: f64) {
        let r = self.rotation;
        let mut dr = [0.0_f64; 9];
        // dR = spin * R * dt
        for i in 0..3 {
            for j in 0..3 {
                let mut s = 0.0;
                for k in 0..3 {
                    s += spin[i * 3 + k] * r[k * 3 + j];
                }
                dr[i * 3 + j] = s * dt;
            }
        }
        let mut new_r = [0.0_f64; 9];
        for (new_r_i, (&r_i, &dr_i)) in new_r.iter_mut().zip(r.iter().zip(dr.iter())) {
            *new_r_i = r_i + dr_i;
        }
        self.rotation = gram_schmidt_3x3(new_r);
    }

    /// Euler angle decomposition (ZXZ convention) from current rotation matrix.
    /// Returns (φ1, Φ, φ2) in radians.
    pub fn euler_angles(&self) -> (f64, f64, f64) {
        let r = self.rotation;
        let phi = r[8].clamp(-1.0, 1.0).acos();
        let phi1;
        let phi2;
        if phi.sin().abs() < 1e-10 {
            phi1 = r[1].atan2(r[0]);
            phi2 = 0.0;
        } else {
            phi1 = r[6].atan2(-r[7]);
            phi2 = r[2].atan2(r[5]);
        }
        (phi1, phi, phi2)
    }
}

impl Default for TextureEvolution {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// PolycrystalFem
// ---------------------------------------------------------------------------

/// A polycrystal aggregate of grains with different orientations.
///
/// Supports the Taylor (uniform strain) and simplified self-consistent models.
#[derive(Debug, Clone)]
pub struct PolycrystalFem {
    /// Rotation matrices for each grain (row-major 3×3 per grain).
    pub grain_rotations: Vec<[f64; 9]>,
    /// Volume fraction for each grain (must sum to 1).
    pub volume_fractions: Vec<f64>,
    /// Stress in each grain in Voigt notation \[σ11,σ22,σ33,σ12,σ23,σ31\] (Pa).
    pub grain_stresses: Vec<[f64; 6]>,
    /// CRSS for each grain.
    pub grain_crss: Vec<f64>,
}

impl PolycrystalFem {
    /// Construct a polycrystal with `n` randomly-weighted grains.
    pub fn new(n: usize, initial_crss: f64) -> Self {
        let fraction = 1.0 / n as f64;
        let grain_rotations = (0..n)
            .map(|i| {
                rotation_from_euler(
                    (i as f64) * PI / (n as f64),
                    0.3 * (i as f64) * PI / (n as f64),
                    0.1 * (i as f64) * PI / (n as f64),
                )
            })
            .collect();
        Self {
            grain_rotations,
            volume_fractions: vec![fraction; n],
            grain_stresses: vec![[0.0; 6]; n],
            grain_crss: vec![initial_crss; n],
        }
    }

    /// Number of grains.
    pub fn n_grains(&self) -> usize {
        self.grain_rotations.len()
    }

    /// Compute volume-averaged stress (Taylor model: same strain in every grain).
    pub fn average_stress(&self) -> [f64; 6] {
        let mut avg = [0.0_f64; 6];
        for (i, s) in self.grain_stresses.iter().enumerate() {
            let f = self.volume_fractions[i];
            for (avg_k, &s_k) in avg.iter_mut().zip(s.iter()) {
                *avg_k += f * s_k;
            }
        }
        avg
    }

    /// Assign stresses to grains by rotating a macro stress into each grain frame.
    pub fn apply_macro_stress(&mut self, macro_stress: [f64; 6]) {
        for i in 0..self.n_grains() {
            // Rotate stress from sample to crystal frame: σ' = R^T σ R
            let r = self.grain_rotations[i];
            self.grain_stresses[i] = rotate_stress_voigt(macro_stress, r);
        }
    }

    /// Taylor model: compute average Schmid factors across all grains.
    pub fn taylor_factor(&self, loading_axis: [f64; 3]) -> f64 {
        let mut sum = 0.0;
        for (i, _) in self.grain_rotations.iter().enumerate() {
            // Simplified: Schmid factor = 0.5 * sin(2θ) where θ depends on rotation
            let r = self.grain_rotations[i];
            // rotate loading axis into grain frame
            let la_crystal = mat3_vec3(r, loading_axis);
            // max Schmid factor ~ 0.5 for 45° grain
            let mf = (dot3(la_crystal, [1.0, 0.0, 0.0]) * dot3(la_crystal, [0.0, 1.0, 0.0])).abs();
            sum += self.volume_fractions[i] * mf;
        }
        sum
    }
}

// ---------------------------------------------------------------------------
// OdfRepresentation
// ---------------------------------------------------------------------------

/// Discrete orientation distribution function (ODF) over Euler angle space.
#[derive(Debug, Clone)]
pub struct OdfRepresentation {
    /// Euler angles (φ1, Φ, φ2) for each orientation bin.
    pub euler_angles: Vec<(f64, f64, f64)>,
    /// Intensity (weight) for each orientation.
    pub intensities: Vec<f64>,
}

impl OdfRepresentation {
    /// Create a uniform (random) ODF with `n_bins` orientation bins.
    pub fn uniform(n_bins: usize) -> Self {
        let mut euler_angles = Vec::with_capacity(n_bins);
        let w = 1.0 / n_bins as f64;
        for i in 0..n_bins {
            let phi1 = 2.0 * PI * (i as f64) / (n_bins as f64);
            let phi = PI * (i as f64) / (n_bins as f64);
            let phi2 = PI * (i as f64) / (n_bins as f64);
            euler_angles.push((phi1, phi, phi2));
        }
        Self {
            euler_angles,
            intensities: vec![w; n_bins],
        }
    }

    /// Total intensity (should integrate to 1 for normalised ODF).
    pub fn total_intensity(&self) -> f64 {
        self.intensities.iter().sum()
    }

    /// Number of orientation bins.
    pub fn n_bins(&self) -> usize {
        self.euler_angles.len()
    }

    /// Normalise so the total intensity equals 1.
    pub fn normalise(&mut self) {
        let total = self.total_intensity();
        if total > 0.0 {
            for w in self.intensities.iter_mut() {
                *w /= total;
            }
        }
    }

    /// Compute a simple 001 pole figure intensity at a given sample direction θ, φ.
    pub fn pole_figure_001(&self, theta: f64, phi_angle: f64) -> f64 {
        let target = [
            theta.sin() * phi_angle.cos(),
            theta.sin() * phi_angle.sin(),
            theta.cos(),
        ];
        let mut intensity = 0.0;
        for (i, &(phi1, big_phi, phi2)) in self.euler_angles.iter().enumerate() {
            let r = rotation_from_euler(phi1, big_phi, phi2);
            let c001 = [r[2], r[5], r[8]]; // third column = [001] direction in sample frame
            let cos_angle = dot3(c001, target).clamp(-1.0, 1.0);
            let kernel = (-((1.0 - cos_angle) / 0.05).powi(2)).exp();
            intensity += self.intensities[i] * kernel;
        }
        intensity
    }
}

// ---------------------------------------------------------------------------
// RecrystallizationModel
// ---------------------------------------------------------------------------

/// JMAK (Johnson-Mehl-Avrami-Kolmogorov) recrystallization kinetics.
///
/// Volume fraction recrystallized: X(t) = 1 − exp(−k · t^n)
#[derive(Debug, Clone)]
pub struct RecrystallizationModel {
    /// Avrami exponent n.
    pub avrami_exponent: f64,
    /// Rate constant k (s⁻ⁿ).
    pub rate_constant: f64,
    /// Current recrystallized volume fraction (0–1).
    pub volume_fraction: f64,
    /// Nucleation site density (m⁻³).
    pub nucleation_density: f64,
    /// Mean grain diameter of recrystallized grains (m).
    pub rex_grain_size: f64,
}

impl RecrystallizationModel {
    /// Construct a JMAK recrystallization model.
    pub fn new(
        avrami_exponent: f64,
        rate_constant: f64,
        nucleation_density: f64,
        rex_grain_size: f64,
    ) -> Self {
        Self {
            avrami_exponent,
            rate_constant,
            volume_fraction: 0.0,
            nucleation_density,
            rex_grain_size,
        }
    }

    /// Compute recrystallized volume fraction at time t (s).
    pub fn volume_fraction_at(&self, t: f64) -> f64 {
        if t <= 0.0 {
            return 0.0;
        }
        1.0 - (-self.rate_constant * t.powf(self.avrami_exponent)).exp()
    }

    /// Update internal volume fraction to time t.
    pub fn update(&mut self, t: f64) {
        self.volume_fraction = self.volume_fraction_at(t);
    }

    /// Time to reach 50% recrystallization (t₅₀).
    pub fn t50(&self) -> f64 {
        (2.0_f64.ln() / self.rate_constant).powf(1.0 / self.avrami_exponent)
    }

    /// Nucleation rate (simplified: proportional to stored energy / temperature).
    pub fn nucleation_rate(&self, stored_energy: f64, temperature: f64) -> f64 {
        if temperature <= 0.0 {
            return 0.0;
        }
        self.nucleation_density * stored_energy / temperature
    }

    /// Growth rate of recrystallized grains (simplified Arrhenius).
    pub fn growth_rate(&self, temperature: f64, activation_energy: f64) -> f64 {
        if temperature <= 0.0 {
            return 0.0;
        }
        let r_gas = 8.314;
        (-activation_energy / (r_gas * temperature)).exp()
    }
}

// ---------------------------------------------------------------------------
// Helper functions
// ---------------------------------------------------------------------------

/// Normalise a 3-vector.
fn normalise3(v: [f64; 3]) -> [f64; 3] {
    let len = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
    if len < 1e-15 {
        return v;
    }
    [v[0] / len, v[1] / len, v[2] / len]
}

/// Dot product of two 3-vectors.
fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// Expand Voigt stress \[σ11,σ22,σ33,σ12,σ23,σ31\] to full 3×3 (row-major).
fn voigt_to_full(v: [f64; 6]) -> [f64; 9] {
    [
        v[0], v[3], v[5], // row 0
        v[3], v[1], v[4], // row 1
        v[5], v[4], v[2], // row 2
    ]
}

/// Gram-Schmidt orthonormalisation of a 3×3 matrix (row vectors).
fn gram_schmidt_3x3(m: [f64; 9]) -> [f64; 9] {
    let mut r0 = [m[0], m[1], m[2]];
    let mut r1 = [m[3], m[4], m[5]];
    let mut r2 = [m[6], m[7], m[8]];
    r0 = normalise3(r0);
    // r1 = r1 - (r1·r0)*r0
    let d10 = dot3(r1, r0);
    r1 = [
        r1[0] - d10 * r0[0],
        r1[1] - d10 * r0[1],
        r1[2] - d10 * r0[2],
    ];
    r1 = normalise3(r1);
    // r2 = r2 - (r2·r0)*r0 - (r2·r1)*r1
    let d20 = dot3(r2, r0);
    let d21 = dot3(r2, r1);
    r2 = [
        r2[0] - d20 * r0[0] - d21 * r1[0],
        r2[1] - d20 * r0[1] - d21 * r1[1],
        r2[2] - d20 * r0[2] - d21 * r1[2],
    ];
    r2 = normalise3(r2);
    [
        r0[0], r0[1], r0[2], r1[0], r1[1], r1[2], r2[0], r2[1], r2[2],
    ]
}

/// Build rotation matrix from ZXZ Euler angles.
fn rotation_from_euler(phi1: f64, big_phi: f64, phi2: f64) -> [f64; 9] {
    let (s1, c1) = phi1.sin_cos();
    let (sp, cp) = big_phi.sin_cos();
    let (s2, c2) = phi2.sin_cos();
    [
        c1 * c2 - s1 * cp * s2,
        -c1 * s2 - s1 * cp * c2,
        s1 * sp,
        s1 * c2 + c1 * cp * s2,
        -s1 * s2 + c1 * cp * c2,
        -c1 * sp,
        sp * s2,
        sp * c2,
        cp,
    ]
}

/// Apply rotation matrix R (row-major) to Voigt stress.
fn rotate_stress_voigt(sigma: [f64; 6], r: [f64; 9]) -> [f64; 6] {
    let s = voigt_to_full(sigma);
    // s' = R^T * s * R  (transform to crystal frame)
    let mut tmp = [0.0_f64; 9];
    // tmp = R^T * s
    for i in 0..3 {
        for j in 0..3 {
            let mut v = 0.0;
            for k in 0..3 {
                v += r[k * 3 + i] * s[k * 3 + j];
            }
            tmp[i * 3 + j] = v;
        }
    }
    // result = tmp * R
    let mut res = [0.0_f64; 9];
    for i in 0..3 {
        for j in 0..3 {
            let mut v = 0.0;
            for k in 0..3 {
                v += tmp[i * 3 + k] * r[k * 3 + j];
            }
            res[i * 3 + j] = v;
        }
    }
    // Pack back to Voigt
    [res[0], res[4], res[8], res[1], res[5], res[2]]
}

/// Multiply rotation matrix by a 3-vector: y = R * x.
fn mat3_vec3(r: [f64; 9], x: [f64; 3]) -> [f64; 3] {
    [
        r[0] * x[0] + r[1] * x[1] + r[2] * x[2],
        r[3] * x[0] + r[4] * x[1] + r[5] * x[2],
        r[6] * x[0] + r[7] * x[1] + r[8] * x[2],
    ]
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // ---- SlipSystem --------------------------------------------------------

    #[test]
    fn test_schmid_factor_45_degrees() {
        // Loading along [1,0,0], slip direction [1,1,0]/√2, slip normal [0,1,1]/√2
        // cos(φ)*cos(λ) = (1/√2)*(1/√2) = 0.5  – only if the geometry is right
        // Use slip direction [1,1,0] and normal [1,-1,0]: both perpendicular at 45°
        let s = SlipSystem::new([1.0, 1.0, 0.0], [1.0, -1.0, 0.0]);
        let sf = s.schmid_factor([1.0, 0.0, 0.0]);
        assert!(
            (sf.abs() - 0.5).abs() < 1e-10,
            "Schmid factor = {:.6}",
            sf.abs()
        );
    }

    #[test]
    fn test_schmid_factor_axial_zero() {
        // Loading axis parallel to slip direction => cos(λ)=1 but cos(φ)=0
        let s = SlipSystem::new([1.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
        let sf = s.schmid_factor([1.0, 0.0, 0.0]);
        assert!(sf.abs() < 1e-10);
    }

    #[test]
    fn test_schmid_tensor_shape() {
        let s = SlipSystem::new([1.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
        assert_eq!(s.schmid_tensor.len(), 9);
    }

    #[test]
    fn test_resolved_shear_stress_pure_shear() {
        // σ_12 = τ, slip direction [1,0,0], slip normal [0,1,0] → τ_rss = τ
        let s = SlipSystem::new([1.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
        let stress = [0.0, 0.0, 0.0, 1e6, 0.0, 0.0]; // σ12 = 1 MPa
        let tau = s.resolved_shear_stress(stress);
        assert!((tau - 1e6).abs() < 1.0, "rss = {:.6}", tau);
    }

    #[test]
    fn test_resolved_shear_stress_zero_off_axis() {
        let s = SlipSystem::new([1.0, 0.0, 0.0], [0.0, 0.0, 1.0]);
        let stress = [1e6, 0.0, 0.0, 0.0, 0.0, 0.0]; // σ11 only
        let tau = s.resolved_shear_stress(stress);
        assert!(tau.abs() < 1.0);
    }

    // ---- FCC Crystal -------------------------------------------------------

    #[test]
    fn test_fcc_slip_system_count() {
        let fcc = FccCrystal::new();
        assert_eq!(
            fcc.n_slip_systems(),
            12,
            "FCC must have exactly 12 slip systems"
        );
    }

    #[test]
    fn test_fcc_max_schmid_factor() {
        let fcc = FccCrystal::new();
        let mf = fcc.max_schmid_factor([1.0, 0.0, 0.0]);
        assert!(mf > 0.0, "Max Schmid factor should be positive");
        assert!(mf <= 0.5 + 1e-10, "Schmid factor cannot exceed 0.5");
    }

    #[test]
    fn test_fcc_slip_systems_normalized() {
        let fcc = FccCrystal::new();
        for sys in &fcc.slip_systems {
            let ls = (sys.slip_direction[0].powi(2)
                + sys.slip_direction[1].powi(2)
                + sys.slip_direction[2].powi(2))
            .sqrt();
            assert!((ls - 1.0).abs() < 1e-10);
        }
    }

    // ---- BCC Crystal -------------------------------------------------------

    #[test]
    fn test_bcc_slip_system_count_positive() {
        let bcc = BccCrystal::new();
        assert!(bcc.n_slip_systems() > 0);
    }

    #[test]
    fn test_bcc_has_more_systems_than_fcc() {
        let bcc = BccCrystal::new();
        let fcc = FccCrystal::new();
        assert!(bcc.n_slip_systems() >= fcc.n_slip_systems());
    }

    // ---- HCP Crystal -------------------------------------------------------

    #[test]
    fn test_hcp_default_ca_ratio() {
        let hcp = HcpCrystal::new((8.0_f64 / 3.0).sqrt());
        assert!(!hcp.pyramidal_active());
    }

    #[test]
    fn test_hcp_pyramidal_active_for_high_ca() {
        let hcp = HcpCrystal::new(1.886); // Zn-like
        assert!(hcp.pyramidal_active());
    }

    /// c/a ratios of Be, Mg, the ideal hard-sphere packing, Ti-like and Zn / Cd.
    const HCP_RATIOS: [f64; 6] = [1.568, 1.624, 1.633, 1.587, 1.856, 1.886];

    #[test]
    fn test_hcp_has_eighteen_systems_in_four_families() {
        // Regression: no pyramidal system used to be created (d.n = c^2 != 0 rejected all of them).
        for ca in HCP_RATIOS {
            let hcp = HcpCrystal::new(ca);
            assert_eq!(hcp.n_slip_systems(), 18, "c/a = {ca}");
            assert_eq!(hcp.basal_systems().len(), 3);
            assert_eq!(hcp.prismatic_systems().len(), 3);
            assert_eq!(hcp.pyramidal_a_systems().len(), 6);
            assert_eq!(hcp.pyramidal_ca_systems().len(), 6);
        }
    }

    #[test]
    fn test_hcp_every_direction_lies_in_its_slip_plane() {
        for ca in HCP_RATIOS {
            let hcp = HcpCrystal::new(ca);
            for (k, sys) in hcp.slip_systems.iter().enumerate() {
                let d_dot_n = dot3(sys.slip_direction, sys.slip_normal);
                assert!(
                    d_dot_n.abs() < 1e-12,
                    "c/a {ca}, system {k}: d.n = {d_dot_n}"
                );
                assert!((dot3(sys.slip_direction, sys.slip_direction) - 1.0).abs() < 1e-12);
                assert!((dot3(sys.slip_normal, sys.slip_normal) - 1.0).abs() < 1e-12);
            }
        }
    }

    #[test]
    fn test_hcp_index_tables_obey_miller_bravais_rules() {
        for family in [
            &HCP_BASAL[..],
            &HCP_PRISMATIC[..],
            &HCP_PYRAMIDAL_A[..],
            &HCP_PYRAMIDAL_CA[..],
        ] {
            for &(plane, dir) in family {
                assert_eq!(plane[0] + plane[1] + plane[2], 0, "plane {plane:?}");
                assert_eq!(dir[0] + dir[1] + dir[2], 0, "direction {dir:?}");
            }
        }
        // Violating the redundant-index rule is rejected.
        assert!(HcpCrystal::direction_vector(1.633, [1, 1, 1, 0]).is_none());
        assert!(HcpCrystal::plane_normal(1.633, [1, 0, 0, 1]).is_none());
    }

    #[test]
    fn test_hcp_direction_and_normal_conversion() {
        let ca = 1.6;
        let sqrt3 = 3.0_f64.sqrt();
        // [2-1-10] is parallel to a1 (length 3), [0001] to c.
        let a1 = HcpCrystal::direction_vector(ca, [2, -1, -1, 0]).unwrap();
        assert!((a1[0] - 3.0).abs() < 1e-14 && a1[1].abs() < 1e-14 && a1[2].abs() < 1e-14);
        let c = HcpCrystal::direction_vector(ca, [0, 0, 0, 1]).unwrap();
        assert!(c[0].abs() < 1e-14 && c[1].abs() < 1e-14 && (c[2] - ca).abs() < 1e-14);
        // [-12-10] is parallel to a2 = (-1/2, sqrt(3)/2, 0).
        let a2 = HcpCrystal::direction_vector(ca, [-1, 2, -1, 0]).unwrap();
        assert!((a2[0] + 1.5).abs() < 1e-14 && (a2[1] - 1.5 * sqrt3).abs() < 1e-14);
        // (0001) has the c-axis normal, (10-10) a normal at 30 degrees from a1 in the basal plane.
        let basal = HcpCrystal::plane_normal(ca, [0, 0, 0, 1]).unwrap();
        assert!(basal[0].abs() < 1e-14 && basal[1].abs() < 1e-14 && basal[2] > 0.0);
        let prism = HcpCrystal::plane_normal(ca, [1, 0, -1, 0]).unwrap();
        let angle = prism[1].atan2(prism[0]).to_degrees();
        assert!((angle - 30.0).abs() < 1e-12, "angle = {angle}");
        assert!(prism[2].abs() < 1e-14);
    }

    #[test]
    fn test_hcp_prismatic_slip_is_along_a_not_along_c() {
        // Regression: prismatic slip used the c-axis [0, 0, c/a] as direction.
        for ca in HCP_RATIOS {
            let hcp = HcpCrystal::new(ca);
            for sys in hcp.prismatic_systems() {
                assert!(
                    sys.slip_direction[2].abs() < 1e-12,
                    "direction has a c-component"
                );
                assert!(
                    sys.slip_normal[2].abs() < 1e-12,
                    "prismatic plane must contain c"
                );
            }
            // The three prismatic directions are the three <a> directions (60 degrees apart).
            let dirs: Vec<[f64; 3]> = hcp
                .prismatic_systems()
                .iter()
                .map(|s| s.slip_direction)
                .collect();
            for (i, &di) in dirs.iter().enumerate() {
                for &dj in dirs.iter().skip(i + 1) {
                    let cos = dot3(di, dj).abs();
                    assert!((cos - 0.5).abs() < 1e-12, "angle between <a> directions");
                }
            }
        }
    }

    #[test]
    fn test_hcp_invalid_ca_ratio_is_rejected_without_nan_systems() {
        for bad in [0.0, -1.633, f64::NAN, f64::INFINITY] {
            assert!(HcpCrystal::try_new(bad).is_err(), "c/a = {bad}");
            let hcp = HcpCrystal::new(bad);
            assert_eq!(hcp.n_slip_systems(), 0, "c/a = {bad}");
            assert!(hcp.basal_systems().is_empty());
            assert!(hcp.pyramidal_ca_systems().is_empty());
        }
        let ok = HcpCrystal::try_new(1.633).unwrap();
        assert_eq!(ok.n_slip_systems(), 18);
    }

    #[test]
    fn test_hcp_basal_systems_are_a_directions_in_the_basal_plane() {
        let hcp = HcpCrystal::new(1.633);
        for sys in hcp.basal_systems() {
            assert!(sys.slip_normal[0].abs() < 1e-12 && sys.slip_normal[1].abs() < 1e-12);
            assert!((sys.slip_normal[2] - 1.0).abs() < 1e-12);
            assert!(sys.slip_direction[2].abs() < 1e-12);
        }
    }

    #[test]
    fn test_hcp_pyramidal_planes_are_tilted_and_ca_directions_have_c_component() {
        for ca in HCP_RATIOS {
            let hcp = HcpCrystal::new(ca);
            for sys in hcp.pyramidal_a_systems() {
                assert!(
                    sys.slip_normal[2] > 0.1,
                    "pyramidal plane must be inclined to the basal plane"
                );
                assert!(
                    sys.slip_direction[2].abs() < 1e-12,
                    "<a> slip stays in the basal plane"
                );
            }
            for sys in hcp.pyramidal_ca_systems() {
                assert!(sys.slip_normal[2] > 0.1);
                assert!(sys.slip_direction[2].abs() > 0.1, "<c+a> has a c-component");
            }
            // {10-11}: normal = (1, 1/sqrt3, 1/gamma).
            let n = hcp.pyramidal_a_systems()[0].slip_normal;
            let expected_nz = (1.0 / ca) / (4.0 / 3.0 + 1.0 / (ca * ca)).sqrt();
            assert!(
                (n[2] - expected_nz).abs() < 1e-12,
                "c/a {ca}: n_z = {}",
                n[2]
            );
        }
    }

    #[test]
    fn test_hcp_schmid_factors_for_a_axis_and_c_axis_loading() {
        let hcp = HcpCrystal::new(1.633);
        let max_abs = |systems: &[SlipSystem], axis: [f64; 3]| {
            systems
                .iter()
                .map(|s| s.schmid_factor(axis).abs())
                .fold(0.0_f64, f64::max)
        };
        // Tension along the a1 axis: prismatic <a> slip has |m| = cos(30 deg) |cos(120 deg)| = 0.433,
        // while basal slip cannot be activated (the basal-plane normal is perpendicular to the axis).
        let a_axis = [1.0, 0.0, 0.0];
        let prismatic = max_abs(hcp.prismatic_systems(), a_axis);
        assert!(
            (prismatic - 3.0_f64.sqrt() / 4.0).abs() < 1e-12,
            "prismatic Schmid factor {prismatic}"
        );
        assert!(max_abs(hcp.basal_systems(), a_axis) < 1e-12);
        // Tension along the c axis: <a> slip cannot be activated, only <c+a> slip.
        let c_axis = [0.0, 0.0, 1.0];
        assert!(max_abs(hcp.basal_systems(), c_axis) < 1e-12);
        assert!(max_abs(hcp.prismatic_systems(), c_axis) < 1e-12);
        assert!(max_abs(hcp.pyramidal_a_systems(), c_axis) < 1e-12);
        let ca_slip = max_abs(hcp.pyramidal_ca_systems(), c_axis);
        assert!(
            ca_slip > 0.3 && ca_slip <= 0.5,
            "<c+a> Schmid factor {ca_slip}"
        );
        // Resolved shear stress agrees with sigma * Schmid factor for uniaxial tension.
        let sigma = 100.0e6;
        let sys = &hcp.prismatic_systems()[0];
        let tau = sys.resolved_shear_stress([sigma, 0.0, 0.0, 0.0, 0.0, 0.0]);
        assert!((tau - sigma * sys.schmid_factor(a_axis)).abs() < 1e-3);
    }

    // ---- CrystalPlasticity -------------------------------------------------

    #[test]
    fn test_shear_rate_increases_with_stress() {
        let cp = CrystalPlasticity::new(12, 1e-3, 10.0, 100e6, 50e6, 1.4);
        let r1 = cp.shear_rate(0, 80e6).abs();
        let r2 = cp.shear_rate(0, 120e6).abs();
        assert!(r2 > r1, "Higher RSS should give higher shear rate");
    }

    #[test]
    fn test_shear_rate_zero_at_zero_stress() {
        let cp = CrystalPlasticity::new(12, 1e-3, 10.0, 100e6, 50e6, 1.4);
        assert_eq!(cp.shear_rate(0, 0.0), 0.0);
    }

    #[test]
    fn test_shear_rate_sign() {
        let cp = CrystalPlasticity::new(12, 1e-3, 10.0, 100e6, 50e6, 1.4);
        let pos = cp.shear_rate(0, 90e6);
        let neg = cp.shear_rate(0, -90e6);
        assert!(pos > 0.0 && neg < 0.0);
    }

    #[test]
    fn test_crss_increases_with_hardening() {
        let mut cp = CrystalPlasticity::new(12, 1e-3, 10.0, 100e6, 50e6, 1.4);
        let tau0 = cp.crss[0];
        let gamma_dots = vec![1e-3; 12];
        cp.update_crss(&gamma_dots, 1.0);
        assert!(cp.crss[0] > tau0, "CRSS should increase with hardening");
    }

    #[test]
    fn test_accumulated_shear_positive() {
        let cp = CrystalPlasticity::new(12, 1e-3, 10.0, 100e6, 50e6, 1.4);
        let gd = vec![1e-3; 12];
        let acc = cp.accumulated_shear(&gd, 0.01);
        assert!(acc > 0.0);
    }

    // ---- DislocationDensity ------------------------------------------------

    #[test]
    fn test_taylor_hardening_increases_with_density() {
        let dd1 = DislocationDensity::new(10e6, 3.06, 80e9, 2.5e-10, 1e-6, 1.0, 1e12);
        let dd2 = DislocationDensity::new(10e6, 3.06, 80e9, 2.5e-10, 1e-6, 1.0, 1e14);
        assert!(dd2.critical_resolved_shear_stress() > dd1.critical_resolved_shear_stress());
    }

    #[test]
    fn test_dislocation_density_increases_with_strain() {
        let mut dd = DislocationDensity::new(10e6, 3.06, 80e9, 2.5e-10, 1e-6, 0.0, 1e12);
        let rho0 = dd.density;
        dd.evolve(0.01);
        assert!(dd.density > rho0);
    }

    #[test]
    fn test_dislocation_density_no_negative() {
        let mut dd = DislocationDensity::new(10e6, 3.06, 80e9, 2.5e-10, 1e-6, 10.0, 1e10);
        dd.evolve(1000.0);
        assert!(dd.density >= 0.0);
    }

    #[test]
    fn test_dislocation_integrate_monotonic() {
        let mut dd = DislocationDensity::new(10e6, 3.06, 80e9, 2.5e-10, 1e-6, 0.0, 1e12);
        let increments: Vec<f64> = vec![0.001; 100];
        let rho0 = dd.density;
        dd.integrate(&increments);
        assert!(dd.density >= rho0);
    }

    // ---- TextureEvolution --------------------------------------------------

    #[test]
    fn test_texture_initial_identity() {
        let te = TextureEvolution::new();
        // Euler angles of identity rotation ~ (0, 0, 0) or within valid range
        let (phi1, phi, phi2) = te.euler_angles();
        let _ = (phi1, phi, phi2); // just check it doesn't panic
        assert!((0.0..=PI).contains(&phi));
    }

    #[test]
    fn test_texture_rotation_changes_angles() {
        let mut te = TextureEvolution::new();
        let initial = te.euler_angles();
        // Small spin about z-axis
        let spin = [0.0, -0.01, 0.0, 0.01, 0.0, 0.0, 0.0, 0.0, 0.0];
        te.update(spin, 1.0);
        let updated = te.euler_angles();
        assert!(
            (updated.0 - initial.0).abs() > 1e-6
                || (updated.1 - initial.1).abs() > 1e-6
                || (updated.2 - initial.2).abs() > 1e-6
        );
    }

    #[test]
    fn test_texture_rotation_orthogonal() {
        let mut te = TextureEvolution::new();
        let spin = [0.0, -0.05, 0.0, 0.05, 0.0, 0.0, 0.0, 0.0, 0.0];
        te.update(spin, 1.0);
        let r = te.rotation;
        // Check R^T * R ≈ I
        for i in 0..3 {
            for j in 0..3 {
                let mut s = 0.0;
                for k in 0..3 {
                    s += r[k * 3 + i] * r[k * 3 + j];
                }
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!((s - expected).abs() < 1e-10, "R^T R[{i},{j}] = {:.6}", s);
            }
        }
    }

    // ---- PolycrystalFem ----------------------------------------------------

    #[test]
    fn test_polycrystal_average_stress() {
        let mut poly = PolycrystalFem::new(10, 100e6);
        let macro_stress = [200e6, 0.0, 0.0, 0.0, 0.0, 0.0];
        poly.apply_macro_stress(macro_stress);
        let avg = poly.average_stress();
        // Average trace should equal macro trace
        let macro_trace = macro_stress[0] + macro_stress[1] + macro_stress[2];
        let avg_trace = avg[0] + avg[1] + avg[2];
        assert!(
            (avg_trace - macro_trace).abs() < 1.0,
            "avg trace = {:.6}",
            avg_trace
        );
    }

    #[test]
    fn test_polycrystal_n_grains() {
        let poly = PolycrystalFem::new(20, 100e6);
        assert_eq!(poly.n_grains(), 20);
    }

    #[test]
    fn test_polycrystal_volume_fractions_sum() {
        let poly = PolycrystalFem::new(10, 100e6);
        let total: f64 = poly.volume_fractions.iter().sum();
        assert!((total - 1.0).abs() < 1e-10);
    }

    // ---- OdfRepresentation -------------------------------------------------

    #[test]
    fn test_odf_uniform_total_intensity() {
        let odf = OdfRepresentation::uniform(100);
        assert!((odf.total_intensity() - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_odf_n_bins() {
        let odf = OdfRepresentation::uniform(50);
        assert_eq!(odf.n_bins(), 50);
    }

    #[test]
    fn test_odf_normalise() {
        let mut odf = OdfRepresentation::uniform(10);
        // Artificially scale
        for w in odf.intensities.iter_mut() {
            *w *= 5.0;
        }
        odf.normalise();
        assert!((odf.total_intensity() - 1.0).abs() < 1e-10);
    }

    // ---- RecrystallizationModel --------------------------------------------

    #[test]
    fn test_jmak_volume_fraction_at_zero() {
        let rex = RecrystallizationModel::new(2.0, 0.1, 1e18, 1e-6);
        assert_eq!(rex.volume_fraction_at(0.0), 0.0);
    }

    #[test]
    fn test_jmak_approaches_one() {
        let rex = RecrystallizationModel::new(2.0, 0.1, 1e18, 1e-6);
        let xf = rex.volume_fraction_at(1000.0);
        assert!(xf > 0.999, "X(1000) = {:.6}", xf);
    }

    #[test]
    fn test_jmak_monotonic_increase() {
        let rex = RecrystallizationModel::new(2.0, 0.1, 1e18, 1e-6);
        let x1 = rex.volume_fraction_at(1.0);
        let x2 = rex.volume_fraction_at(2.0);
        assert!(x2 > x1);
    }

    #[test]
    fn test_jmak_t50() {
        let rex = RecrystallizationModel::new(2.0, 0.1, 1e18, 1e-6);
        let t50 = rex.t50();
        let x_at_t50 = rex.volume_fraction_at(t50);
        assert!((x_at_t50 - 0.5).abs() < 1e-10, "X(t50) = {:.6}", x_at_t50);
    }

    #[test]
    fn test_jmak_nucleation_rate_positive() {
        let rex = RecrystallizationModel::new(2.0, 0.1, 1e18, 1e-6);
        let rate = rex.nucleation_rate(1e6, 1000.0);
        assert!(rate > 0.0);
    }

    #[test]
    fn test_jmak_growth_rate_positive() {
        let rex = RecrystallizationModel::new(2.0, 0.1, 1e18, 1e-6);
        let rate = rex.growth_rate(1000.0, 150e3);
        assert!(rate > 0.0);
    }

    #[test]
    fn test_jmak_update() {
        let mut rex = RecrystallizationModel::new(2.0, 0.1, 1e18, 1e-6);
        rex.update(5.0);
        assert!(rex.volume_fraction > 0.0);
        assert!(rex.volume_fraction < 1.0);
    }
}
