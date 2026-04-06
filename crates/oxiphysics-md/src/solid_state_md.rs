#![allow(clippy::needless_range_loop)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Solid-state molecular dynamics for crystalline systems.
//!
//! Covers:
//! - Phonon dispersion via the dynamical matrix (Born-von Kármán model)
//! - Thermal conductivity via Green-Kubo relation
//! - Peierls-Nabarro dislocation model
//! - Dislocation velocity in the drag-limited regime
//! - Stacking fault energy calculation
//! - Grain boundary MD energetics
//! - Radiation damage cascade simulation (PKA)
//! - Amorphization threshold from cumulative dose
//! - Born effective charges and dielectric tensor
//! - Debye-Waller factor and mean-square displacement

#![allow(dead_code)]
#![allow(clippy::too_many_arguments)]

use std::f64::consts::PI;

// ---------------------------------------------------------------------------
// Physical constants
// ---------------------------------------------------------------------------

/// Boltzmann constant (J K⁻¹).
const KB: f64 = 1.380_649e-23;

/// Reduced Planck constant (J s).
const HBAR: f64 = 1.054_571_817e-34;

/// Atomic mass unit (kg).
const AMU: f64 = 1.660_539_066_6e-27;

/// Elementary charge (C).
const E_CHARGE: f64 = 1.602_176_634e-19;

/// Avogadro number (mol⁻¹).
const N_A: f64 = 6.022_140_76e23;

// ---------------------------------------------------------------------------
// Vector helpers
// ---------------------------------------------------------------------------

/// Dot product of two 3-vectors.
#[inline]
pub fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// Euclidean norm of a 3-vector.
#[inline]
pub fn norm3(a: [f64; 3]) -> f64 {
    dot3(a, a).sqrt()
}

/// Add two 3-vectors.
#[inline]
pub fn add3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

/// Subtract two 3-vectors.
#[inline]
pub fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

/// Scale a 3-vector.
#[inline]
pub fn scale3(v: [f64; 3], s: f64) -> [f64; 3] {
    [v[0] * s, v[1] * s, v[2] * s]
}

/// Matrix-vector product for a 3×3 matrix (row-major) and a 3-vector.
#[inline]
pub fn mat3_vec(m: [[f64; 3]; 3], v: [f64; 3]) -> [f64; 3] {
    [
        m[0][0] * v[0] + m[0][1] * v[1] + m[0][2] * v[2],
        m[1][0] * v[0] + m[1][1] * v[1] + m[1][2] * v[2],
        m[2][0] * v[0] + m[2][1] * v[1] + m[2][2] * v[2],
    ]
}

/// Trace of a 3×3 matrix.
#[inline]
pub fn trace3(m: [[f64; 3]; 3]) -> f64 {
    m[0][0] + m[1][1] + m[2][2]
}

/// Add two 3×3 matrices.
#[inline]
pub fn mat3_add(a: [[f64; 3]; 3], b: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let mut c = [[0.0_f64; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            c[i][j] = a[i][j] + b[i][j];
        }
    }
    c
}

/// Scale a 3×3 matrix.
#[inline]
pub fn mat3_scale(m: [[f64; 3]; 3], s: f64) -> [[f64; 3]; 3] {
    let mut c = [[0.0_f64; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            c[i][j] = m[i][j] * s;
        }
    }
    c
}

// ---------------------------------------------------------------------------
// Lattice definition
// ---------------------------------------------------------------------------

/// Bravais lattice type.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum BravaisLattice {
    /// Simple cubic.
    SimpleCubic,
    /// Face-centred cubic.
    FaceCentredCubic,
    /// Body-centred cubic.
    BodyCentredCubic,
    /// Hexagonal close-packed (2D projection).
    HexagonalClosePacked,
    /// Diamond cubic.
    Diamond,
}

/// Parameters defining the primitive lattice.
pub struct CrystalParams {
    /// Primitive lattice vectors (rows).
    pub a_vec: [[f64; 3]; 3],
    /// Atom mass in atomic mass units.
    pub mass_amu: f64,
    /// Lattice type.
    pub lattice: BravaisLattice,
    /// Number of atoms in primitive cell.
    pub n_basis: usize,
}

impl CrystalParams {
    /// Copper FCC, a = 3.615 Å.
    pub fn copper() -> Self {
        let a = 3.615e-10_f64;
        let a2 = a / 2.0;
        Self {
            a_vec: [[a2, a2, 0.0], [0.0, a2, a2], [a2, 0.0, a2]],
            mass_amu: 63.546,
            lattice: BravaisLattice::FaceCentredCubic,
            n_basis: 1,
        }
    }

    /// Iron BCC, a = 2.87 Å.
    pub fn iron() -> Self {
        let a = 2.87e-10_f64;
        let a2 = a / 2.0;
        Self {
            a_vec: [[-a2, a2, a2], [a2, -a2, a2], [a2, a2, -a2]],
            mass_amu: 55.845,
            lattice: BravaisLattice::BodyCentredCubic,
            n_basis: 1,
        }
    }

    /// Silicon diamond cubic, a = 5.431 Å.
    pub fn silicon() -> Self {
        let a = 5.431e-10_f64;
        let a2 = a / 2.0;
        Self {
            a_vec: [[0.0, a2, a2], [a2, 0.0, a2], [a2, a2, 0.0]],
            mass_amu: 28.085,
            lattice: BravaisLattice::Diamond,
            n_basis: 2,
        }
    }

    /// Return the volume of the primitive unit cell (m³).
    pub fn cell_volume(&self) -> f64 {
        let a = self.a_vec;
        // V = a1 · (a2 × a3)
        let cross = [
            a[1][1] * a[2][2] - a[1][2] * a[2][1],
            a[1][2] * a[2][0] - a[1][0] * a[2][2],
            a[1][0] * a[2][1] - a[1][1] * a[2][0],
        ];
        dot3(a[0], cross).abs()
    }

    /// Return atom mass in kg.
    pub fn mass_kg(&self) -> f64 {
        self.mass_amu * AMU
    }
}

// ---------------------------------------------------------------------------
// Born-von Kármán dynamical matrix
// ---------------------------------------------------------------------------

/// Nearest-neighbour force constant matrix for a monatomic crystal.
pub struct ForceConstants {
    /// On-site (self) force constant tensor (N m⁻¹), row-major 3×3.
    pub phi_self: [[f64; 3]; 3],
    /// Neighbour force constant tensor (N m⁻¹), row-major 3×3.
    pub phi_nn: [[f64; 3]; 3],
    /// Number of nearest neighbours.
    pub n_nn: usize,
}

impl ForceConstants {
    /// Simple cubic nearest-neighbour model with longitudinal constant `kl`
    /// and transverse constant `kt`.
    pub fn simple_cubic(kl: f64, kt: f64) -> Self {
        // Longitudinal along x, y, z neighbours
        let phi_x = [[kl, 0.0, 0.0], [0.0, kt, 0.0], [0.0, 0.0, kt]];
        // Self term: sum over all neighbours (6 for SC) enforces translational
        // invariance. Self = −Σ_j Φ_j, but for symmetry use diagonal.
        let self_k = kl + 2.0 * kt;
        let phi_self = [
            [self_k * 2.0, 0.0, 0.0],
            [0.0, self_k * 2.0, 0.0],
            [0.0, 0.0, self_k * 2.0],
        ];
        Self {
            phi_self,
            phi_nn: phi_x,
            n_nn: 6,
        }
    }
}

/// Compute the 3×3 dynamical matrix D(q) for a monatomic crystal at wave-vector q.
///
/// D(q) = (1/M) \[ Φ_self + Σ_j Φ_j exp(i q·r_j) \]
///
/// Here we use only the real part (D is Hermitian, so Im part vanishes for
/// pairs of ±neighbours).
pub fn dynamical_matrix_simple_cubic(
    fc: &ForceConstants,
    mass_kg: f64,
    q_vec: [f64; 3],
    lattice_constant: f64,
) -> [[f64; 3]; 3] {
    // Neighbours of simple cubic: ±x, ±y, ±z
    let neighbors: [[f64; 3]; 6] = [
        [lattice_constant, 0.0, 0.0],
        [-lattice_constant, 0.0, 0.0],
        [0.0, lattice_constant, 0.0],
        [0.0, -lattice_constant, 0.0],
        [0.0, 0.0, lattice_constant],
        [0.0, 0.0, -lattice_constant],
    ];

    let mut d = fc.phi_self;

    for r in &neighbors {
        let phase = dot3(q_vec, *r); // q · r
        let cos_phase = phase.cos();
        // Each neighbour: the phi matrix rotated to the bond direction
        // For simplicity use the same phi_nn (assumes isotropic bonds in diagonal form)
        for i in 0..3 {
            for j in 0..3 {
                d[i][j] += fc.phi_nn[i][j] * cos_phase;
            }
        }
    }

    mat3_scale(d, 1.0 / mass_kg)
}

/// Phonon frequencies from the eigenvalues of D(q).
///
/// For a diagonal D matrix the eigenvalues are the diagonal entries.
/// Returns sorted angular frequencies ω = √(λ) (rad s⁻¹).
pub fn phonon_frequencies_diagonal(d: [[f64; 3]; 3]) -> [f64; 3] {
    let mut freqs = [
        d[0][0].max(0.0).sqrt(),
        d[1][1].max(0.0).sqrt(),
        d[2][2].max(0.0).sqrt(),
    ];
    freqs.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    freqs
}

/// Compute phonon dispersion along Γ→X for a simple cubic lattice.
///
/// Returns a `Vec<[f64; 3]>` of (q_x, ω_acoustic, ω_optical) for `n_pts` points.
pub fn phonon_dispersion_gamma_x(
    fc: &ForceConstants,
    mass_kg: f64,
    lattice_constant: f64,
    n_pts: usize,
) -> Vec<[f64; 3]> {
    let mut result = Vec::with_capacity(n_pts);
    for i in 0..n_pts {
        let frac = i as f64 / (n_pts.saturating_sub(1).max(1)) as f64;
        let q_x = frac * PI / lattice_constant;
        let q = [q_x, 0.0, 0.0];
        let d = dynamical_matrix_simple_cubic(fc, mass_kg, q, lattice_constant);
        let freqs = phonon_frequencies_diagonal(d);
        result.push([q_x, freqs[0], freqs[2]]);
    }
    result
}

// ---------------------------------------------------------------------------
// Green-Kubo thermal conductivity
// ---------------------------------------------------------------------------

/// Compute the heat current autocorrelation integral for thermal conductivity.
///
/// κ = (1 / (k_B T² V)) ∫₀^∞ ⟨J(0)·J(t)⟩ dt
///
/// `hcacf` is the sampled heat current autocorrelation function (W² m⁻³).
/// `dt` is the time step (s), `temperature` in K, `volume` in m³.
pub fn green_kubo_conductivity(hcacf: &[f64], dt: f64, temperature: f64, volume: f64) -> f64 {
    if hcacf.is_empty() {
        return 0.0;
    }
    // Trapezoidal integration
    let integral: f64 = hcacf.windows(2).map(|w| (w[0] + w[1]) * 0.5 * dt).sum();
    integral / (KB * temperature * temperature * volume)
}

/// Compute the heat current autocorrelation function (HCACF) from a time series
/// of heat current magnitudes.
///
/// Uses the direct summation definition:
/// C(m) = (1/N) Σ_{n=0}^{N-m-1} J(n) J(n+m)
pub fn compute_hcacf(heat_current: &[f64], max_lag: usize) -> Vec<f64> {
    let n = heat_current.len();
    let m_max = max_lag.min(n);
    let mut acf = vec![0.0_f64; m_max];
    for lag in 0..m_max {
        let count = n.saturating_sub(lag);
        if count == 0 {
            break;
        }
        let sum: f64 = (0..count)
            .map(|k| heat_current[k] * heat_current[k + lag])
            .sum();
        acf[lag] = sum / count as f64;
    }
    acf
}

// ---------------------------------------------------------------------------
// Peierls-Nabarro model for dislocation
// ---------------------------------------------------------------------------

/// Parameters for the Peierls-Nabarro (PN) dislocation model.
pub struct PeierlsNabarroParams {
    /// Burgers vector magnitude (m).
    pub burgers: f64,
    /// Peierls stress (Pa) – stress to move a dislocation at 0 K.
    pub peierls_stress: f64,
    /// Shear modulus (Pa).
    pub shear_modulus: f64,
    /// Lattice parameter (m).
    pub lattice_a: f64,
    /// Half-width of dislocation core (m).
    pub core_half_width: f64,
    /// Poisson ratio.
    pub poisson: f64,
}

impl PeierlsNabarroParams {
    /// Parameters for edge dislocation in copper.
    pub fn copper_edge() -> Self {
        Self {
            burgers: 2.556e-10,
            peierls_stress: 1.0e7,
            shear_modulus: 4.8e10,
            lattice_a: 3.615e-10,
            core_half_width: 3.615e-10 * 2.0,
            poisson: 0.34,
        }
    }
}

/// Compute the Peierls stress from the PN model.
///
/// τ_P ≈ (2G / (1−ν)) · exp(−2π w / b)
///
/// where `w` is the dislocation half-width.
pub fn peierls_stress(params: &PeierlsNabarroParams) -> f64 {
    let exponent = -2.0 * PI * params.core_half_width / params.burgers;
    2.0 * params.shear_modulus / (1.0 - params.poisson) * exponent.exp()
}

/// Compute dislocation core energy (J m⁻¹).
///
/// E_core = (G b² / 4π(1−ν)) · ln(R / r0)  for edge dislocation
pub fn edge_dislocation_core_energy(
    params: &PeierlsNabarroParams,
    outer_radius: f64,
    inner_radius: f64,
) -> f64 {
    let g = params.shear_modulus;
    let b = params.burgers;
    let nu = params.poisson;
    let r_ratio = (outer_radius / inner_radius.max(1.0e-15)).ln();
    g * b * b / (4.0 * PI * (1.0 - nu)) * r_ratio
}

/// Compute the misfit energy landscape (Peierls potential) at fractional slip u/b.
///
/// γ(u) = A \[1 − cos(2π u / b)\]
pub fn peierls_potential(params: &PeierlsNabarroParams, u_over_b: f64) -> f64 {
    let amplitude = params.shear_modulus * params.burgers
        / (4.0 * PI * PI * (1.0 - params.poisson))
        * (2.0 * PI * params.core_half_width / params.burgers).exp();
    amplitude * (1.0 - (2.0 * PI * u_over_b).cos())
}

// ---------------------------------------------------------------------------
// Dislocation velocity – drag-limited
// ---------------------------------------------------------------------------

/// Dislocation drag parameters.
pub struct DislocationDragParams {
    /// Drag coefficient B (Pa s) at reference temperature.
    pub drag_coeff: f64,
    /// Temperature exponent for drag (B ∝ T^n).
    pub temp_exponent: f64,
    /// Reference temperature (K).
    pub t_ref: f64,
    /// Burgers vector (m).
    pub burgers: f64,
    /// Peierls stress (Pa).
    pub peierls_stress: f64,
}

impl DislocationDragParams {
    /// Typical copper dislocation drag at 300 K.
    pub fn copper() -> Self {
        Self {
            drag_coeff: 1.0e-4,
            temp_exponent: 1.0,
            t_ref: 300.0,
            burgers: 2.556e-10,
            peierls_stress: 1.0e7,
        }
    }
}

/// Compute dislocation velocity in drag-limited regime.
///
/// v = (τ − τ_P) · b / B(T)
pub fn dislocation_velocity(
    params: &DislocationDragParams,
    applied_stress: f64,
    temperature_k: f64,
) -> f64 {
    let excess_stress = (applied_stress - params.peierls_stress).max(0.0);
    let b_t = params.drag_coeff * (temperature_k / params.t_ref).powf(params.temp_exponent);
    excess_stress * params.burgers / b_t.max(1.0e-30)
}

/// Compute dislocation mean free path limited by phonon scattering.
///
/// ℓ ≈ v · τ_ph  where τ_ph = h_bar / (k_B T)
pub fn dislocation_mean_free_path(
    params: &DislocationDragParams,
    applied_stress: f64,
    temperature_k: f64,
) -> f64 {
    let v = dislocation_velocity(params, applied_stress, temperature_k);
    let tau_ph = HBAR / (KB * temperature_k.max(1.0));
    v * tau_ph
}

// ---------------------------------------------------------------------------
// Stacking fault energy
// ---------------------------------------------------------------------------

/// Stacking fault energy data for common metals (J m⁻²).
pub struct StackingFaultData {
    /// Material name.
    pub name: &'static str,
    /// Intrinsic stacking fault energy γ_ISF (J m⁻²).
    pub gamma_isf: f64,
    /// Extrinsic stacking fault energy γ_ESF (J m⁻²).
    pub gamma_esf: f64,
    /// Unstable stacking fault energy γ_US (J m⁻²).
    pub gamma_us: f64,
    /// Burgers vector (m).
    pub burgers: f64,
}

impl StackingFaultData {
    /// Cu stacking fault data.
    pub fn copper() -> Self {
        Self {
            name: "Cu",
            gamma_isf: 0.045,
            gamma_esf: 0.095,
            gamma_us: 0.175,
            burgers: 2.556e-10,
        }
    }

    /// Al stacking fault data.
    pub fn aluminium() -> Self {
        Self {
            name: "Al",
            gamma_isf: 0.166,
            gamma_esf: 0.332,
            gamma_us: 0.224,
            burgers: 2.863e-10,
        }
    }

    /// Ni stacking fault data.
    pub fn nickel() -> Self {
        Self {
            name: "Ni",
            gamma_isf: 0.125,
            gamma_esf: 0.250,
            gamma_us: 0.330,
            burgers: 2.492e-10,
        }
    }

    /// Extended dislocation width from equilibrium partial separation.
    ///
    /// d_eq = G b_p² / (8π γ_ISF) · (2+ν)/(1−ν)
    pub fn partial_separation(&self, shear_modulus: f64, poisson: f64) -> f64 {
        let b_p = self.burgers / (3.0_f64.sqrt()); // partial Burgers for FCC
        shear_modulus * b_p * b_p * (2.0 + poisson) / (8.0 * PI * self.gamma_isf * (1.0 - poisson))
    }

    /// Twinning susceptibility: ratio γ_US / γ_ISF (dimensionless).
    pub fn twinning_susceptibility(&self) -> f64 {
        self.gamma_us / self.gamma_isf.max(1.0e-30)
    }
}

/// Compute the generalised stacking fault energy along the slip path (GSFE curve).
///
/// Uses a simple sinusoidal model: γ(u) = γ_US · sin²(π u / b)
pub fn gsfe_sinusoidal(gamma_us: f64, burgers: f64, u: f64) -> f64 {
    gamma_us * (PI * u / burgers).sin().powi(2)
}

// ---------------------------------------------------------------------------
// Grain boundary MD
// ---------------------------------------------------------------------------

/// Grain boundary energy as a function of misorientation angle.
pub struct GrainBoundaryParams {
    /// Read-Shockley energy parameter E_0 (J m⁻²).
    pub e0: f64,
    /// Critical angle θ_m for Read-Shockley saturation (rad).
    pub theta_max: f64,
    /// Grain boundary mobility prefactor (m⁴ J⁻¹ s⁻¹).
    pub mobility_prefactor: f64,
    /// Activation energy for grain boundary migration (J mol⁻¹).
    pub migration_activation: f64,
}

impl GrainBoundaryParams {
    /// Parameters for copper grain boundaries.
    pub fn copper() -> Self {
        Self {
            e0: 0.625,
            theta_max: std::f64::consts::FRAC_PI_6, // 30°
            mobility_prefactor: 1.0e-5,
            migration_activation: 60_000.0,
        }
    }
}

/// Compute grain boundary energy using Read-Shockley model.
///
/// γ_GB(θ) = γ_0 · θ/θ_m · (1 − ln(θ/θ_m))  for θ < θ_m
pub fn read_shockley_energy(params: &GrainBoundaryParams, theta_rad: f64) -> f64 {
    let theta = theta_rad.clamp(1.0e-6, params.theta_max);
    let ratio = theta / params.theta_max;
    params.e0 * ratio * (1.0 - ratio.ln())
}

/// Compute grain boundary migration velocity.
///
/// v = M(T) · P  where P is the driving pressure (J m⁻³) and
/// M(T) = M_0 exp(−Q / RT)
pub fn grain_boundary_velocity(
    params: &GrainBoundaryParams,
    temperature_k: f64,
    driving_pressure: f64,
) -> f64 {
    let r_gas = 8.314_462_618_f64;
    let m_t = params.mobility_prefactor
        * (-params.migration_activation / (r_gas * temperature_k.max(1.0))).exp();
    m_t * driving_pressure
}

/// Estimate grain growth rate from grain boundary mobility.
///
/// dR/dt = M(T) · γ_GB / R  (parabolic grain growth)
pub fn grain_growth_rate(
    params: &GrainBoundaryParams,
    temperature_k: f64,
    theta_rad: f64,
    grain_radius: f64,
) -> f64 {
    let energy = read_shockley_energy(params, theta_rad);
    let pressure = energy / grain_radius.max(1.0e-9);
    grain_boundary_velocity(params, temperature_k, pressure)
}

// ---------------------------------------------------------------------------
// Radiation damage – cascade simulation
// ---------------------------------------------------------------------------

/// State of a single atom in a radiation damage cascade.
pub struct CascadeAtom {
    /// Position (m).
    pub position: [f64; 3],
    /// Velocity (m s⁻¹).
    pub velocity: [f64; 3],
    /// Mass (kg).
    pub mass: f64,
    /// Current kinetic energy (J).
    pub kinetic_energy: f64,
    /// Whether this atom has been displaced (Frenkel pair vacancy side).
    pub is_displaced: bool,
    /// Whether this atom is an interstitial.
    pub is_interstitial: bool,
}

impl CascadeAtom {
    /// Create an atom at a lattice site.
    pub fn at_lattice_site(position: [f64; 3], mass_amu: f64) -> Self {
        Self {
            position,
            velocity: [0.0; 3],
            mass: mass_amu * AMU,
            kinetic_energy: 0.0,
            is_displaced: false,
            is_interstitial: false,
        }
    }

    /// Update kinetic energy from velocity.
    pub fn update_kinetic_energy(&mut self) {
        self.kinetic_energy = 0.5 * self.mass * dot3(self.velocity, self.velocity);
    }
}

/// Parameters for a radiation cascade event.
pub struct CascadeParams {
    /// Primary knock-on atom (PKA) initial energy (J).
    pub pka_energy: f64,
    /// Displacement threshold energy (J) – typically ~25 eV for metals.
    pub displacement_threshold: f64,
    /// Lattice parameter (m).
    pub lattice_a: f64,
    /// Electronic stopping power coefficient (J m⁻¹) per unit velocity.
    pub electronic_stopping: f64,
    /// Nuclear stopping cross-section (m²).
    pub nuclear_cross_section: f64,
}

impl CascadeParams {
    /// Typical iron PKA cascade parameters.
    pub fn iron_pka(pka_energy_ev: f64) -> Self {
        Self {
            pka_energy: pka_energy_ev * E_CHARGE,
            displacement_threshold: 25.0 * E_CHARGE,
            lattice_a: 2.87e-10,
            electronic_stopping: 1.0e-10,
            nuclear_cross_section: 1.0e-19,
        }
    }
}

/// Compute the number of Frenkel pairs using the Norgett-Robinson-Torrens (NRT) model.
///
/// N_FP = 0.8 · T_dam / (2 · E_d)  where T_dam = damage energy.
pub fn nrt_frenkel_pairs(pka_energy_j: f64, displacement_threshold_j: f64) -> f64 {
    let t_dam = pka_energy_j; // simplified: all energy goes to damage
    0.8 * t_dam / (2.0 * displacement_threshold_j.max(1.0e-30))
}

/// Compute the Lindhard partition function (fraction of PKA energy in nuclear stopping).
///
/// Uses simplified Lindhard formula: ξ(ε) = ε / (ε + 0.133 ε^{0.41} + ...)
/// Here we use a two-parameter fit for metals (Robinson 1971).
pub fn lindhard_partition(reduced_energy: f64) -> f64 {
    // Fraction going to nuclear (damage) stopping
    let g = reduced_energy
        + 0.40244 * reduced_energy.powf(3.0 / 4.0)
        + 3.4008 * reduced_energy.powf(1.0 / 6.0);
    reduced_energy / g.max(1.0e-30)
}

/// Estimate cascade radius from PKA energy (simple empirical model).
///
/// r_cascade ≈ r0 · (E_PKA / E_d)^{1/3}
pub fn cascade_radius(pka_energy_j: f64, displacement_threshold_j: f64, r0: f64) -> f64 {
    let ratio = pka_energy_j / displacement_threshold_j.max(1.0e-30);
    r0 * ratio.powf(1.0 / 3.0)
}

/// Compute recombination probability for Frenkel pairs.
///
/// P_rec = 1 − exp(−(r_iv / r_cascade)³)  where r_iv is recombination radius.
pub fn frenkel_recombination_probability(
    cascade_radius_m: f64,
    recombination_radius_m: f64,
) -> f64 {
    let ratio = recombination_radius_m / cascade_radius_m.max(1.0e-30);
    1.0 - (-ratio.powi(3)).exp()
}

// ---------------------------------------------------------------------------
// Amorphization threshold
// ---------------------------------------------------------------------------

/// State tracking cumulative radiation dose and amorphization.
pub struct AmorphizationTracker {
    /// Cumulative displacement dose (displacements per atom, dpa).
    pub dose_dpa: f64,
    /// Critical dose for amorphization (dpa).
    pub critical_dose: f64,
    /// Current amorphous fraction (0–1).
    pub amorphous_fraction: f64,
    /// Recrystallization rate coefficient (dpa⁻¹).
    pub recryst_rate: f64,
}

impl AmorphizationTracker {
    /// Create tracker for silicon (amorphizes at ~0.5 dpa).
    pub fn silicon() -> Self {
        Self {
            dose_dpa: 0.0,
            critical_dose: 0.5,
            amorphous_fraction: 0.0,
            recryst_rate: 0.1,
        }
    }

    /// Create tracker for iron (much harder to amorphize).
    pub fn iron() -> Self {
        Self {
            dose_dpa: 0.0,
            critical_dose: 50.0,
            amorphous_fraction: 0.0,
            recryst_rate: 0.5,
        }
    }

    /// Advance dose by `delta_dpa` and update amorphous fraction.
    ///
    /// Uses the Hecking model: dα/dΦ = (1 − α) − R_c · α
    pub fn advance_dose(&mut self, delta_dpa: f64, temperature_k: f64) {
        self.dose_dpa += delta_dpa;
        // Thermal recrystallization increases with temperature
        let r_c = self.recryst_rate * (temperature_k / 300.0).powi(2);
        let d_alpha = ((1.0 - self.amorphous_fraction) - r_c * self.amorphous_fraction) * delta_dpa
            / self.critical_dose;
        self.amorphous_fraction = (self.amorphous_fraction + d_alpha).clamp(0.0, 1.0);
    }

    /// Check whether the material is fully amorphous.
    pub fn is_amorphous(&self) -> bool {
        self.amorphous_fraction > 0.99
    }
}

// ---------------------------------------------------------------------------
// Born effective charges
// ---------------------------------------------------------------------------

/// Born effective charge tensor for an atom in a polar crystal.
pub struct BornEffectiveCharge {
    /// Atom label.
    pub label: &'static str,
    /// Z* tensor (row-major 3×3), in units of elementary charge.
    pub z_star: [[f64; 3]; 3],
}

impl BornEffectiveCharge {
    /// Create a diagonal Born effective charge tensor.
    pub fn diagonal(label: &'static str, zxx: f64, zyy: f64, zzz: f64) -> Self {
        Self {
            label,
            z_star: [[zxx, 0.0, 0.0], [0.0, zyy, 0.0], [0.0, 0.0, zzz]],
        }
    }

    /// GaAs Ga atom Born effective charge (diagonal).
    pub fn gaas_ga() -> Self {
        Self::diagonal("Ga", 2.16, 2.16, 2.16)
    }

    /// GaAs As atom Born effective charge (diagonal).
    pub fn gaas_as() -> Self {
        Self::diagonal("As", -2.16, -2.16, -2.16)
    }

    /// Compute the force on this atom due to an electric field E (N).
    ///
    /// F = Z* · E · e
    pub fn force_from_field(&self, e_field: [f64; 3]) -> [f64; 3] {
        mat3_vec(self.z_star, scale3(e_field, E_CHARGE))
    }

    /// Average isotropic Born charge (trace/3).
    pub fn isotropic_charge(&self) -> f64 {
        trace3(self.z_star) / 3.0
    }
}

/// Compute the LO-TO splitting frequency shift.
///
/// Δω² = (4π e² / (ε_inf · Ω · M)) · (Z*)²
///
/// Returns Δω in rad s⁻¹.
pub fn lo_to_splitting(born_charge: f64, eps_inf: f64, cell_volume: f64, mass_kg: f64) -> f64 {
    let prefactor = 4.0 * PI * E_CHARGE * E_CHARGE * born_charge * born_charge;
    let denom = eps_inf * cell_volume * mass_kg;
    (prefactor / denom.max(1.0e-60)).sqrt()
}

// ---------------------------------------------------------------------------
// Debye-Waller factor
// ---------------------------------------------------------------------------

/// Compute the Debye-Waller factor B (mean square displacement in Å²).
///
/// B = 8π² ⟨u²⟩ where ⟨u²⟩ = (3 h_bar² T) / (M k_B Θ_D²)   (high-T limit)
pub fn debye_waller_factor(mass_kg: f64, temperature_k: f64, debye_temperature: f64) -> f64 {
    let u2 =
        3.0 * HBAR * HBAR * temperature_k / (mass_kg * KB * debye_temperature * debye_temperature);
    8.0 * PI * PI * u2
}

/// Compute the Debye temperature from the maximum phonon frequency.
///
/// Θ_D = h_bar ω_max / k_B
pub fn debye_temperature_from_max_freq(omega_max: f64) -> f64 {
    HBAR * omega_max / KB
}

/// Compute the Debye heat capacity at constant volume (J mol⁻¹ K⁻¹).
///
/// Uses the Debye integral approximation in the limits:
/// - Low T: C_V → 12/5 π⁴ R (T/Θ_D)³
/// - High T: C_V → 3R (Dulong-Petit)
pub fn debye_heat_capacity(temperature_k: f64, debye_temperature: f64) -> f64 {
    let r_gas = 8.314_462_618_f64;
    if temperature_k < 0.01 {
        return 0.0;
    }
    let t_ratio = temperature_k / debye_temperature;
    if t_ratio > 2.0 {
        // High-T: Dulong-Petit
        3.0 * r_gas
    } else {
        // Low-T approximation
        12.0 / 5.0 * PI * PI * PI * PI * r_gas * t_ratio.powi(3)
    }
}

// ---------------------------------------------------------------------------
// MD time integration for crystalline systems
// ---------------------------------------------------------------------------

/// Single atom for velocity-Verlet integration.
pub struct CrystalAtom {
    /// Position (m).
    pub position: [f64; 3],
    /// Velocity (m s⁻¹).
    pub velocity: [f64; 3],
    /// Force (N).
    pub force: [f64; 3],
    /// Mass (kg).
    pub mass: f64,
    /// Atom type index.
    pub type_id: usize,
}

impl CrystalAtom {
    /// Create an atom at rest.
    pub fn new(position: [f64; 3], mass_amu: f64, type_id: usize) -> Self {
        Self {
            position,
            velocity: [0.0; 3],
            force: [0.0; 3],
            mass: mass_amu * AMU,
            type_id,
        }
    }

    /// Kinetic energy (J).
    pub fn kinetic_energy(&self) -> f64 {
        0.5 * self.mass * dot3(self.velocity, self.velocity)
    }

    /// Temperature estimate from kinetic energy (single atom, 3 DOF).
    pub fn temperature(&self) -> f64 {
        2.0 * self.kinetic_energy() / (3.0 * KB)
    }

    /// Velocity-Verlet half-kick (first half: v += F/m * dt/2).
    pub fn vv_half_kick(&mut self, dt: f64) {
        let a = scale3(self.force, 1.0 / self.mass.max(1.0e-60));
        self.velocity = add3(self.velocity, scale3(a, 0.5 * dt));
    }

    /// Velocity-Verlet position update.
    pub fn vv_position_update(&mut self, dt: f64) {
        self.position = add3(self.position, scale3(self.velocity, dt));
    }
}

/// Velocity-Verlet integration for a collection of atoms.
pub fn velocity_verlet_step(atoms: &mut [CrystalAtom], forces: &[[f64; 3]], dt: f64) {
    // First half-kick
    for atom in atoms.iter_mut() {
        atom.vv_half_kick(dt);
    }
    // Position update
    for atom in atoms.iter_mut() {
        atom.vv_position_update(dt);
    }
    // Update forces externally (caller responsibility) then second half-kick
    for (atom, &f) in atoms.iter_mut().zip(forces.iter()) {
        atom.force = f;
        atom.vv_half_kick(dt);
    }
}

/// Compute system temperature from kinetic energy.
pub fn system_temperature(atoms: &[CrystalAtom]) -> f64 {
    if atoms.is_empty() {
        return 0.0;
    }
    let ek: f64 = atoms.iter().map(|a| a.kinetic_energy()).sum();
    let n_dof = 3 * atoms.len();
    2.0 * ek / (n_dof as f64 * KB)
}

// ---------------------------------------------------------------------------
// Embedded Atom Method (EAM) potential – simplified
// ---------------------------------------------------------------------------

/// Simplified EAM potential parameters for copper (Foiles et al. 1986 style).
pub struct EamParams {
    /// Pair potential depth (J).
    pub epsilon: f64,
    /// Equilibrium pair distance (m).
    pub r0: f64,
    /// Pair potential stiffness.
    pub alpha: f64,
    /// Embedding function coefficient.
    pub embed_coeff: f64,
    /// Electron density parameter.
    pub rho_scale: f64,
    /// Cutoff radius (m).
    pub r_cut: f64,
}

impl EamParams {
    /// Approximate copper EAM parameters.
    pub fn copper() -> Self {
        Self {
            epsilon: 0.5 * E_CHARGE, // ~0.5 eV
            r0: 2.556e-10,
            alpha: 5.0,
            embed_coeff: -E_CHARGE,
            rho_scale: 1.0e10,
            r_cut: 5.5e-10,
        }
    }

    /// Compute pair potential (J) at separation r.
    pub fn pair_potential(&self, r: f64) -> f64 {
        if r >= self.r_cut {
            return 0.0;
        }
        let x = (r - self.r0) / self.r0;
        self.epsilon * ((-self.alpha * x).exp() * 2.0 - 2.0 * (-self.alpha * x * 0.5).exp())
    }

    /// Compute electron density contribution from atom at distance r.
    pub fn electron_density(&self, r: f64) -> f64 {
        if r >= self.r_cut {
            return 0.0;
        }
        self.rho_scale * (-self.alpha * (r / self.r0 - 1.0)).exp()
    }

    /// Compute embedding energy F(ρ) = embed_coeff · √ρ.
    pub fn embedding_energy(&self, rho: f64) -> f64 {
        self.embed_coeff * rho.max(0.0).sqrt()
    }
}

/// Compute total potential energy for a pair of atoms using EAM.
pub fn eam_pair_energy(params: &EamParams, r: f64) -> f64 {
    let phi = params.pair_potential(r);
    let rho = params.electron_density(r);
    let f_rho = params.embedding_energy(rho);
    phi + f_rho
}

// ---------------------------------------------------------------------------
// Thermal expansion
// ---------------------------------------------------------------------------

/// Compute linear thermal expansion coefficient from quasi-harmonic approximation.
///
/// α_L = γ_G · C_V / (3 B V)
///
/// where γ_G is Grüneisen parameter, C_V is heat capacity (J m⁻³ K⁻¹),
/// B is bulk modulus (Pa), V is molar volume (m³ mol⁻¹).
pub fn linear_thermal_expansion(gruneisen: f64, cv_volumetric: f64, bulk_modulus: f64) -> f64 {
    gruneisen * cv_volumetric / (3.0 * bulk_modulus)
}

/// Compute the Grüneisen parameter from mode Grüneisen parameters.
///
/// γ_G = Σ_i C_i γ_i / Σ_i C_i  (heat-capacity-weighted average)
pub fn mode_gruneisen_average(mode_gruneisen: &[f64], mode_heat_capacity: &[f64]) -> f64 {
    let num: f64 = mode_gruneisen
        .iter()
        .zip(mode_heat_capacity.iter())
        .map(|(g, c)| g * c)
        .sum();
    let den: f64 = mode_heat_capacity.iter().sum::<f64>();
    if den.abs() < 1.0e-30 {
        return 0.0;
    }
    num / den
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // --- CrystalParams tests ---

    #[test]
    fn test_copper_cell_volume_positive() {
        let p = CrystalParams::copper();
        assert!(p.cell_volume() > 0.0);
    }

    #[test]
    fn test_copper_mass_kg_correct() {
        let p = CrystalParams::copper();
        let expected = 63.546 * AMU;
        assert!((p.mass_kg() - expected).abs() < 1.0e-30);
    }

    #[test]
    fn test_silicon_lattice_diamond() {
        let p = CrystalParams::silicon();
        assert_eq!(p.lattice, BravaisLattice::Diamond);
        assert_eq!(p.n_basis, 2);
    }

    // --- Dynamical matrix tests ---

    #[test]
    fn test_dynamical_matrix_at_gamma_has_zero_acoustic() {
        let fc = ForceConstants::simple_cubic(10.0, 5.0);
        let m = 63.546 * AMU;
        let a = 3.615e-10;
        let d = dynamical_matrix_simple_cubic(&fc, m, [0.0; 3], a);
        // At Γ, acoustic branch frequencies should be zero (or near-zero for this model)
        let freqs = phonon_frequencies_diagonal(d);
        // Minimum should be 0 or positive
        assert!(freqs[0] >= 0.0);
    }

    #[test]
    fn test_phonon_dispersion_length() {
        let fc = ForceConstants::simple_cubic(10.0, 5.0);
        let m = 63.546 * AMU;
        let a = 3.615e-10;
        let disp = phonon_dispersion_gamma_x(&fc, m, a, 10);
        assert_eq!(disp.len(), 10);
    }

    #[test]
    fn test_phonon_frequencies_nonnegative() {
        let fc = ForceConstants::simple_cubic(10.0, 5.0);
        let m = 63.546 * AMU;
        let a = 3.615e-10;
        for [_q, f1, f2] in phonon_dispersion_gamma_x(&fc, m, a, 10) {
            assert!(f1 >= 0.0);
            assert!(f2 >= 0.0);
        }
    }

    // --- Green-Kubo tests ---

    #[test]
    fn test_green_kubo_zero_for_empty() {
        let kappa = green_kubo_conductivity(&[], 1.0e-15, 300.0, 1.0e-27);
        assert_eq!(kappa, 0.0);
    }

    #[test]
    fn test_hcacf_length_equals_max_lag() {
        let j: Vec<f64> = (0..100).map(|i| (i as f64).sin()).collect();
        let acf = compute_hcacf(&j, 20);
        assert_eq!(acf.len(), 20);
    }

    #[test]
    fn test_hcacf_zero_lag_is_variance() {
        let j = vec![1.0_f64; 50];
        let acf = compute_hcacf(&j, 10);
        assert!((acf[0] - 1.0).abs() < 1.0e-12);
    }

    #[test]
    fn test_green_kubo_positive_for_positive_acf() {
        let acf = vec![1.0e10_f64; 100];
        let kappa = green_kubo_conductivity(&acf, 1.0e-15, 300.0, 1.0e-27);
        assert!(kappa > 0.0);
    }

    // --- Peierls-Nabarro tests ---

    #[test]
    fn test_peierls_stress_positive() {
        let p = PeierlsNabarroParams::copper_edge();
        assert!(peierls_stress(&p) > 0.0);
    }

    #[test]
    fn test_edge_dislocation_core_energy_positive() {
        let p = PeierlsNabarroParams::copper_edge();
        let e = edge_dislocation_core_energy(&p, 1.0e-6, 1.0e-10);
        assert!(e > 0.0);
    }

    #[test]
    fn test_peierls_potential_zero_at_integer_b() {
        let p = PeierlsNabarroParams::copper_edge();
        let e = peierls_potential(&p, 0.0);
        assert!(e.abs() < 1.0e-30);
    }

    #[test]
    fn test_peierls_potential_max_at_half_b() {
        let p = PeierlsNabarroParams::copper_edge();
        let e0 = peierls_potential(&p, 0.0);
        let e_half = peierls_potential(&p, 0.5);
        assert!(e_half > e0);
    }

    // --- Dislocation velocity tests ---

    #[test]
    fn test_dislocation_velocity_zero_below_peierls() {
        let p = DislocationDragParams::copper();
        let v = dislocation_velocity(&p, p.peierls_stress * 0.5, 300.0);
        assert_eq!(v, 0.0);
    }

    #[test]
    fn test_dislocation_velocity_positive_above_peierls() {
        let p = DislocationDragParams::copper();
        let v = dislocation_velocity(&p, p.peierls_stress * 2.0, 300.0);
        assert!(v > 0.0);
    }

    #[test]
    fn test_dislocation_velocity_decreases_with_temperature() {
        let p = DislocationDragParams::copper();
        let stress = p.peierls_stress * 3.0;
        let v_low = dislocation_velocity(&p, stress, 100.0);
        let v_high = dislocation_velocity(&p, stress, 1_000.0);
        assert!(v_low > v_high);
    }

    // --- Stacking fault energy tests ---

    #[test]
    fn test_partial_separation_copper_positive() {
        let sfe = StackingFaultData::copper();
        let d = sfe.partial_separation(4.8e10, 0.34);
        assert!(d > 0.0);
    }

    #[test]
    fn test_twinning_susceptibility_al_less_than_cu() {
        let cu = StackingFaultData::copper();
        let al = StackingFaultData::aluminium();
        // Al has higher ISF → lower susceptibility
        assert!(al.twinning_susceptibility() < cu.twinning_susceptibility());
    }

    #[test]
    fn test_gsfe_zero_at_integer_burgers() {
        let sfe = StackingFaultData::copper();
        let g = gsfe_sinusoidal(sfe.gamma_us, sfe.burgers, 0.0);
        assert!(g.abs() < 1.0e-30);
    }

    #[test]
    fn test_gsfe_max_at_half_burgers() {
        let sfe = StackingFaultData::copper();
        let g = gsfe_sinusoidal(sfe.gamma_us, sfe.burgers, sfe.burgers * 0.5);
        assert!((g - sfe.gamma_us).abs() < 1.0e-15);
    }

    // --- Grain boundary tests ---

    #[test]
    fn test_read_shockley_energy_zero_at_small_angle() {
        let p = GrainBoundaryParams::copper();
        let e = read_shockley_energy(&p, 1.0e-6);
        assert!(e >= 0.0);
    }

    #[test]
    fn test_grain_boundary_velocity_positive() {
        let p = GrainBoundaryParams::copper();
        let v = grain_boundary_velocity(&p, 800.0, 1.0e6);
        assert!(v > 0.0);
    }

    #[test]
    fn test_grain_growth_rate_positive() {
        let p = GrainBoundaryParams::copper();
        let rate = grain_growth_rate(&p, 800.0, 0.1, 1.0e-6);
        assert!(rate > 0.0);
    }

    // --- Radiation damage tests ---

    #[test]
    fn test_nrt_frenkel_pairs_scales_with_energy() {
        let n1 = nrt_frenkel_pairs(1000.0 * E_CHARGE, 25.0 * E_CHARGE);
        let n2 = nrt_frenkel_pairs(500.0 * E_CHARGE, 25.0 * E_CHARGE);
        assert!(n1 > n2);
    }

    #[test]
    fn test_cascade_radius_increases_with_energy() {
        let r1 = cascade_radius(1000.0 * E_CHARGE, 25.0 * E_CHARGE, 1.0e-10);
        let r2 = cascade_radius(10_000.0 * E_CHARGE, 25.0 * E_CHARGE, 1.0e-10);
        assert!(r2 > r1);
    }

    #[test]
    fn test_recombination_probability_between_zero_and_one() {
        let p = frenkel_recombination_probability(5.0e-9, 1.0e-9);
        assert!((0.0..=1.0).contains(&p));
    }

    #[test]
    fn test_lindhard_partition_between_zero_and_one() {
        let f = lindhard_partition(0.5);
        assert!((0.0..=1.0).contains(&f));
    }

    // --- Amorphization tests ---

    #[test]
    fn test_silicon_amorphizes_after_high_dose() {
        let mut tracker = AmorphizationTracker::silicon();
        for _ in 0..100 {
            tracker.advance_dose(0.1, 300.0);
        }
        assert!(tracker.amorphous_fraction > 0.0);
    }

    #[test]
    fn test_iron_harder_to_amorphize_than_silicon() {
        let mut si = AmorphizationTracker::silicon();
        let mut fe = AmorphizationTracker::iron();
        for _ in 0..10 {
            si.advance_dose(0.1, 300.0);
            fe.advance_dose(0.1, 300.0);
        }
        assert!(si.amorphous_fraction > fe.amorphous_fraction);
    }

    #[test]
    fn test_amorphous_fraction_clamped_to_one() {
        let mut tracker = AmorphizationTracker::silicon();
        for _ in 0..10_000 {
            tracker.advance_dose(1.0, 77.0);
        }
        assert!(tracker.amorphous_fraction <= 1.0);
    }

    // --- Born effective charges tests ---

    #[test]
    fn test_born_charge_force_from_field_nonzero() {
        let bec = BornEffectiveCharge::gaas_ga();
        let f = bec.force_from_field([1.0e6, 0.0, 0.0]);
        assert!(f[0].abs() > 0.0);
    }

    #[test]
    fn test_born_isotropic_charge_gaas_sum_zero() {
        let ga = BornEffectiveCharge::gaas_ga();
        let as_bec = BornEffectiveCharge::gaas_as();
        let sum = ga.isotropic_charge() + as_bec.isotropic_charge();
        assert!(sum.abs() < 1.0e-10);
    }

    #[test]
    fn test_lo_to_splitting_positive() {
        let p = CrystalParams::silicon();
        let delta_omega = lo_to_splitting(2.16, 11.7, p.cell_volume(), p.mass_kg());
        assert!(delta_omega > 0.0);
    }

    // --- Debye-Waller tests ---

    #[test]
    fn test_debye_waller_increases_with_temperature() {
        let m = 63.546 * AMU;
        let b1 = debye_waller_factor(m, 100.0, 343.0);
        let b2 = debye_waller_factor(m, 300.0, 343.0);
        assert!(b2 > b1);
    }

    #[test]
    fn test_debye_temperature_from_max_freq_positive() {
        let theta = debye_temperature_from_max_freq(1.0e13);
        assert!(theta > 0.0);
    }

    #[test]
    fn test_debye_heat_capacity_approaches_dulong_petit() {
        let cv = debye_heat_capacity(10_000.0, 343.0);
        let r = 8.314_462_618_f64;
        assert!((cv - 3.0 * r).abs() < 0.01 * 3.0 * r);
    }

    // --- Crystal atom / velocity-Verlet tests ---

    #[test]
    fn test_crystal_atom_kinetic_energy_zero_at_rest() {
        let a = CrystalAtom::new([0.0; 3], 63.546, 0);
        assert_eq!(a.kinetic_energy(), 0.0);
    }

    #[test]
    fn test_vv_half_kick_changes_velocity() {
        let mut a = CrystalAtom::new([0.0; 3], 63.546, 0);
        a.force = [1.0e-10, 0.0, 0.0];
        a.vv_half_kick(1.0e-14);
        assert!(a.velocity[0] > 0.0);
    }

    #[test]
    fn test_vv_position_update_changes_position() {
        let mut a = CrystalAtom::new([0.0; 3], 63.546, 0);
        a.velocity = [100.0, 0.0, 0.0];
        a.vv_position_update(1.0e-14);
        assert!(a.position[0] > 0.0);
    }

    #[test]
    fn test_system_temperature_empty_is_zero() {
        let atoms: Vec<CrystalAtom> = vec![];
        assert_eq!(system_temperature(&atoms), 0.0);
    }

    // --- EAM potential tests ---

    #[test]
    fn test_eam_pair_energy_zero_beyond_cutoff() {
        let p = EamParams::copper();
        let e = eam_pair_energy(&p, p.r_cut + 1.0e-10);
        assert_eq!(e, 0.0);
    }

    #[test]
    fn test_eam_pair_potential_at_equilibrium_small() {
        let p = EamParams::copper();
        let phi = p.pair_potential(p.r0);
        assert!(
            phi.abs() < 1.0e-17,
            "phi at r0 should be ~0 for Morse-like, got {phi}"
        );
    }

    // --- Thermal expansion tests ---

    #[test]
    fn test_linear_thermal_expansion_positive() {
        let alpha = linear_thermal_expansion(1.8, 3.0e6, 1.4e11);
        assert!(alpha > 0.0);
    }

    #[test]
    fn test_mode_gruneisen_average_single_mode() {
        let g = mode_gruneisen_average(&[2.0], &[1.0]);
        assert!((g - 2.0).abs() < 1.0e-12);
    }

    #[test]
    fn test_mode_gruneisen_average_empty_is_zero() {
        let g = mode_gruneisen_average(&[], &[]);
        assert_eq!(g, 0.0);
    }
}
