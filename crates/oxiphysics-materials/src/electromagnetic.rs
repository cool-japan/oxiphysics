// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Electromagnetic material properties.
//!
//! Covers permittivity / permeability / conductivity tensors,
//! metamaterial properties (negative-index, epsilon-near-zero),
//! the Drude model for metals, the Lorentz oscillator model,
//! the Sellmeier equation for optical glass, magneto-optic effects
//! (Faraday rotation), the Meissner effect in superconductors,
//! ferroelectric and magnetic hysteresis (P-E loop, B-H loop),
//! the Jiles-Atherton hysteresis model, and electromagnetic shielding
//! effectiveness.

use std::f64::consts::PI;

// ---------------------------------------------------------------------------
// Physical constants
// ---------------------------------------------------------------------------

/// Speed of light in vacuum \[m/s\].
pub const C0: f64 = 2.997_924_58e8;

/// Permittivity of free space ε₀ \[F/m\].
pub const EPS0: f64 = 8.854_187_817e-12;

/// Permeability of free space μ₀ \[H/m\].
pub const MU0: f64 = 1.256_637_061_4e-6;

/// Impedance of free space Z₀ \[Ω\].
pub const Z0: f64 = 376.730_313_668;

/// Electron charge \[C\].
pub const E_CHARGE: f64 = 1.602_176_634e-19;

/// Electron mass \[kg\].
pub const E_MASS: f64 = 9.109_383_701_5e-31;

/// Boltzmann constant \[J/K\].
pub const KB: f64 = 1.380_649e-23;

/// Planck constant \[J·s\].
pub const HBAR: f64 = 1.054_571_817e-34;

// ---------------------------------------------------------------------------
// 3×3 tensor type
// ---------------------------------------------------------------------------

/// A 3×3 real-valued tensor stored in row-major order.
///
/// Used for permittivity, permeability, and conductivity of anisotropic
/// media.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Tensor3x3 {
    /// Row-major elements: `data[row][col]`.
    pub data: [[f64; 3]; 3],
}

impl Tensor3x3 {
    /// Create a diagonal tensor with equal diagonal entries `d`.
    pub fn diagonal(d: f64) -> Self {
        let mut t = Self::zero();
        t.data[0][0] = d;
        t.data[1][1] = d;
        t.data[2][2] = d;
        t
    }

    /// Create a zero tensor.
    pub fn zero() -> Self {
        Self {
            data: [[0.0; 3]; 3],
        }
    }

    /// Create a tensor from three diagonal entries.
    pub fn from_diag(d0: f64, d1: f64, d2: f64) -> Self {
        let mut t = Self::zero();
        t.data[0][0] = d0;
        t.data[1][1] = d1;
        t.data[2][2] = d2;
        t
    }

    /// Multiply tensor by a 3-vector: y = T·x.
    pub fn mul_vec(&self, x: [f64; 3]) -> [f64; 3] {
        let mut y = [0.0_f64; 3];
        for (y_i, row) in y.iter_mut().zip(self.data.iter()) {
            *y_i = row
                .iter()
                .zip(x.iter())
                .map(|(&t_ij, &x_j)| t_ij * x_j)
                .sum();
        }
        y
    }

    /// Return the trace of the tensor.
    pub fn trace(&self) -> f64 {
        self.data[0][0] + self.data[1][1] + self.data[2][2]
    }

    /// Return the transpose of the tensor.
    pub fn transpose(&self) -> Self {
        let mut t = Self::zero();
        for i in 0..3 {
            for j in 0..3 {
                t.data[i][j] = self.data[j][i];
            }
        }
        t
    }

    /// Check whether the tensor is symmetric within tolerance `eps`.
    pub fn is_symmetric(&self, eps: f64) -> bool {
        for i in 0..3 {
            for j in 0..3 {
                if (self.data[i][j] - self.data[j][i]).abs() > eps {
                    return false;
                }
            }
        }
        true
    }
}

// ---------------------------------------------------------------------------
// 1. GENERAL ELECTROMAGNETIC MATERIAL
// ---------------------------------------------------------------------------

/// Full set of electromagnetic material properties for a given medium.
#[derive(Debug, Clone)]
pub struct ElectromagneticMaterial {
    /// Human-readable label.
    pub name: String,
    /// Relative permittivity tensor ε_r (dimensionless).
    pub permittivity: Tensor3x3,
    /// Relative permeability tensor μ_r (dimensionless).
    pub permeability: Tensor3x3,
    /// DC conductivity tensor σ \[S/m\].
    pub conductivity: Tensor3x3,
    /// Loss tangent tan(δ) at the reference frequency.
    pub loss_tangent: f64,
    /// Reference frequency for loss tangent \[Hz\].
    pub ref_freq: f64,
}

impl ElectromagneticMaterial {
    /// Create an isotropic non-magnetic dielectric.
    pub fn dielectric(name: &str, eps_r: f64, loss_tangent: f64, ref_freq: f64) -> Self {
        Self {
            name: name.to_owned(),
            permittivity: Tensor3x3::diagonal(eps_r),
            permeability: Tensor3x3::diagonal(1.0),
            conductivity: Tensor3x3::zero(),
            loss_tangent,
            ref_freq,
        }
    }

    /// Create an isotropic conductor (σ >> ωε).
    pub fn conductor(name: &str, sigma: f64) -> Self {
        Self {
            name: name.to_owned(),
            permittivity: Tensor3x3::diagonal(1.0),
            permeability: Tensor3x3::diagonal(1.0),
            conductivity: Tensor3x3::diagonal(sigma),
            loss_tangent: 0.0,
            ref_freq: 1.0e9,
        }
    }

    /// Scalar relative permittivity (trace / 3) for isotropic media.
    pub fn eps_r_scalar(&self) -> f64 {
        self.permittivity.trace() / 3.0
    }

    /// Scalar relative permeability (trace / 3) for isotropic media.
    pub fn mu_r_scalar(&self) -> f64 {
        self.permeability.trace() / 3.0
    }

    /// Intrinsic impedance Z = Z₀ √(μ_r / ε_r) \[Ω\] (isotropic approximation).
    pub fn intrinsic_impedance(&self) -> f64 {
        let eps = self.eps_r_scalar().max(1e-30);
        let mu = self.mu_r_scalar().max(1e-30);
        Z0 * (mu / eps).sqrt()
    }

    /// Phase velocity v_p = c / √(ε_r μ_r) \[m/s\] (isotropic).
    pub fn phase_velocity(&self) -> f64 {
        let eps = self.eps_r_scalar().max(1e-30);
        let mu = self.mu_r_scalar().max(1e-30);
        C0 / (eps * mu).sqrt()
    }

    /// Skin depth δ = 1 / √(π f μ σ) \[m\] at frequency `f` \[Hz\].
    pub fn skin_depth(&self, f: f64) -> f64 {
        let sigma = self.conductivity.trace() / 3.0;
        let mu = self.mu_r_scalar();
        if sigma < 1e-30 || f < 1e-30 {
            return f64::INFINITY;
        }
        1.0 / (PI * f * mu * MU0 * sigma).sqrt()
    }
}

// ---------------------------------------------------------------------------
// 2. METAMATERIALS
// ---------------------------------------------------------------------------

/// Classification of metamaterial index regime.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MetamaterialClass {
    /// ε > 0 and μ > 0: ordinary double-positive (DPS) medium.
    DoublePosive,
    /// ε < 0 and μ < 0: double-negative (DNG) / left-handed medium.
    DoubleNegative,
    /// ε < 0 and μ > 0: epsilon-negative (ENG) medium.
    EpsilonNegative,
    /// ε > 0 and μ < 0: mu-negative (MNG) medium.
    MuNegative,
    /// |ε| < threshold: epsilon-near-zero (ENZ) medium.
    EpsilonNearZero,
    /// |μ| < threshold: mu-near-zero (MNZ) medium.
    MuNearZero,
}

/// Frequency-dispersive metamaterial described by scalar Drude-Lorentz
/// parameters at a single operating frequency.
#[derive(Debug, Clone)]
pub struct Metamaterial {
    /// Label.
    pub name: String,
    /// Real part of relative permittivity ε'_r.
    pub eps_r_real: f64,
    /// Imaginary part of relative permittivity ε''_r (loss).
    pub eps_r_imag: f64,
    /// Real part of relative permeability μ'_r.
    pub mu_r_real: f64,
    /// Imaginary part of relative permeability μ''_r (loss).
    pub mu_r_imag: f64,
    /// Operating frequency \[Hz\].
    pub freq: f64,
}

impl Metamaterial {
    /// Classify the metamaterial based on ε and μ signs.
    pub fn classify(&self) -> MetamaterialClass {
        let enz_thresh = 0.01;
        let mnz_thresh = 0.01;
        if self.eps_r_real.abs() < enz_thresh {
            return MetamaterialClass::EpsilonNearZero;
        }
        if self.mu_r_real.abs() < mnz_thresh {
            return MetamaterialClass::MuNearZero;
        }
        match (self.eps_r_real < 0.0, self.mu_r_real < 0.0) {
            (true, true) => MetamaterialClass::DoubleNegative,
            (true, false) => MetamaterialClass::EpsilonNegative,
            (false, true) => MetamaterialClass::MuNegative,
            (false, false) => MetamaterialClass::DoublePosive,
        }
    }

    /// Refractive index n = √(ε μ).  Negative for DNG media.
    ///
    /// Uses the convention that the imaginary part is negative for lossy DNG.
    pub fn refractive_index(&self) -> f64 {
        let n2 = self.eps_r_real * self.mu_r_real;
        if n2 >= 0.0 {
            // DNG (both negative): product is positive but n is defined negative.
            if self.eps_r_real < 0.0 && self.mu_r_real < 0.0 {
                -(n2.sqrt())
            } else {
                n2.sqrt()
            }
        } else {
            // Single-negative: n is negative real (no imaginary part in this simplified form).
            -((-n2).sqrt())
        }
    }

    /// Group velocity v_g = c / (n + ω dn/dω).  Here we return -C0/n for
    /// ideal DNG as a qualitative indicator.
    pub fn group_velocity_approx(&self) -> f64 {
        let n = self.refractive_index();
        if n.abs() < 1e-30 { C0 } else { C0 / n }
    }
}

// ---------------------------------------------------------------------------
// 3. DRUDE MODEL
// ---------------------------------------------------------------------------

/// Drude model for the frequency-dependent permittivity of a free-electron
/// metal.
///
/// ε(ω) = ε∞ − ωp² / (ω² + iγω)
#[derive(Debug, Clone, Copy)]
pub struct DrudeModel {
    /// High-frequency (background) permittivity ε∞.
    pub eps_inf: f64,
    /// Plasma frequency ωp \[rad/s\].
    pub omega_p: f64,
    /// Collision (damping) rate γ \[rad/s\].
    pub gamma: f64,
}

impl DrudeModel {
    /// Construct a Drude model from material parameters.
    pub fn new(eps_inf: f64, omega_p: f64, gamma: f64) -> Self {
        Self {
            eps_inf,
            omega_p,
            gamma,
        }
    }

    /// Typical Drude parameters for gold at optical frequencies.
    pub fn gold() -> Self {
        Self {
            eps_inf: 9.5,
            omega_p: 1.37e16, // rad/s
            gamma: 1.05e14,   // rad/s
        }
    }

    /// Typical Drude parameters for silver.
    pub fn silver() -> Self {
        Self {
            eps_inf: 5.0,
            omega_p: 1.39e16,
            gamma: 2.73e13,
        }
    }

    /// Real part of the complex permittivity ε'(ω).
    pub fn eps_real(&self, omega: f64) -> f64 {
        self.eps_inf - self.omega_p * self.omega_p / (omega * omega + self.gamma * self.gamma)
    }

    /// Imaginary part of the complex permittivity ε''(ω) (sign convention: -iωt).
    pub fn eps_imag(&self, omega: f64) -> f64 {
        self.omega_p * self.omega_p * self.gamma
            / (omega * (omega * omega + self.gamma * self.gamma))
    }

    /// DC conductivity σ₀ = ε₀ ωp² / γ \[S/m\].
    pub fn dc_conductivity(&self) -> f64 {
        EPS0 * self.omega_p * self.omega_p / self.gamma
    }

    /// Skin depth at angular frequency ω \[rad/s\].
    pub fn skin_depth(&self, omega: f64) -> f64 {
        let eps_im = self.eps_imag(omega);
        if eps_im < 1e-30 {
            return f64::INFINITY;
        }
        // δ ≈ c / (ω √(ε''/2)) for good conductors.
        C0 / (omega * (eps_im / 2.0).sqrt())
    }

    /// Plasma frequency in Hz.
    pub fn plasma_freq_hz(&self) -> f64 {
        self.omega_p / (2.0 * PI)
    }
}

// ---------------------------------------------------------------------------
// 4. LORENTZ OSCILLATOR MODEL
// ---------------------------------------------------------------------------

/// Single-resonance Lorentz oscillator contribution to permittivity.
///
/// ε(ω) = ε∞ + Δε ω₀² / (ω₀² − ω² − iΓω)
#[derive(Debug, Clone, Copy)]
pub struct LorentzOscillator {
    /// Background permittivity.
    pub eps_inf: f64,
    /// Oscillator strength Δε (static − high-freq permittivity).
    pub delta_eps: f64,
    /// Resonance angular frequency ω₀ \[rad/s\].
    pub omega0: f64,
    /// Damping rate Γ \[rad/s\].
    pub gamma: f64,
}

impl LorentzOscillator {
    /// Construct a Lorentz oscillator.
    pub fn new(eps_inf: f64, delta_eps: f64, omega0: f64, gamma: f64) -> Self {
        Self {
            eps_inf,
            delta_eps,
            omega0,
            gamma,
        }
    }

    /// Real part of ε(ω).
    pub fn eps_real(&self, omega: f64) -> f64 {
        let d = (self.omega0 * self.omega0 - omega * omega).powi(2) + (self.gamma * omega).powi(2);
        self.eps_inf
            + self.delta_eps
                * self.omega0
                * self.omega0
                * (self.omega0 * self.omega0 - omega * omega)
                / d
    }

    /// Imaginary part of ε(ω) (absorptive, > 0 for lossy material).
    pub fn eps_imag(&self, omega: f64) -> f64 {
        let d = (self.omega0 * self.omega0 - omega * omega).powi(2) + (self.gamma * omega).powi(2);
        self.delta_eps * self.omega0 * self.omega0 * self.gamma * omega / d
    }

    /// Static permittivity ε(0) = ε∞ + Δε.
    pub fn eps_static(&self) -> f64 {
        self.eps_inf + self.delta_eps
    }

    /// Resonant absorption peak ε''(ω₀).
    pub fn eps_imag_at_resonance(&self) -> f64 {
        self.eps_imag(self.omega0)
    }
}

// ---------------------------------------------------------------------------
// 5. SELLMEIER EQUATION
// ---------------------------------------------------------------------------

/// Sellmeier dispersion model for the refractive index of optical glass.
///
/// n²(λ) = 1 + Σ_i B_i λ² / (λ² − C_i)
///
/// where λ is in micrometres and C_i are resonance wavelengths squared \[μm²\].
#[derive(Debug, Clone)]
pub struct SellmeierModel {
    /// Sellmeier B coefficients (dimensionless oscillator strengths).
    pub b: Vec<f64>,
    /// Sellmeier C coefficients \[μm²\].
    pub c: Vec<f64>,
}

impl SellmeierModel {
    /// Construct a Sellmeier model.
    pub fn new(b: Vec<f64>, c: Vec<f64>) -> Self {
        assert_eq!(b.len(), c.len(), "Sellmeier B and C must have equal length");
        Self { b, c }
    }

    /// Standard Sellmeier coefficients for BK7 borosilicate glass.
    pub fn bk7() -> Self {
        Self::new(
            vec![1.039_612_12, 0.231_792_344, 1.010_469_45],
            vec![0.006_000_699_5, 0.020_017_914, 103.560_653],
        )
    }

    /// Standard Sellmeier coefficients for fused silica (SiO₂).
    pub fn fused_silica() -> Self {
        Self::new(
            vec![0.696_166_3, 0.407_942_6, 0.897_479_4],
            vec![0.004_679_148_6, 0.013_512_063, 97.934_002_5],
        )
    }

    /// Refractive index n(λ) for wavelength `lambda_um` in micrometres.
    pub fn refractive_index(&self, lambda_um: f64) -> f64 {
        let l2 = lambda_um * lambda_um;
        let n2 = 1.0
            + self
                .b
                .iter()
                .zip(self.c.iter())
                .map(|(&bi, &ci)| bi * l2 / (l2 - ci))
                .sum::<f64>();
        n2.max(0.0).sqrt()
    }

    /// Group refractive index n_g = n − λ dn/dλ, estimated by finite difference.
    pub fn group_index(&self, lambda_um: f64) -> f64 {
        let dl = 1e-5_f64;
        let n1 = self.refractive_index(lambda_um - dl);
        let n2 = self.refractive_index(lambda_um + dl);
        let dn_dl = (n2 - n1) / (2.0 * dl);
        self.refractive_index(lambda_um) - lambda_um * dn_dl
    }

    /// Group velocity dispersion GVD = λ³/(2πc²) · d²n/dλ² \[s²/m\] at `lambda_um`.
    pub fn gvd(&self, lambda_um: f64) -> f64 {
        let dl = 1e-4_f64;
        let n0 = self.refractive_index(lambda_um);
        let np = self.refractive_index(lambda_um + dl);
        let nm = self.refractive_index(lambda_um - dl);
        let d2n = (np - 2.0 * n0 + nm) / (dl * dl);
        let lambda_m = lambda_um * 1e-6;
        lambda_m.powi(3) / (2.0 * PI * C0 * C0) * d2n
    }
}

// ---------------------------------------------------------------------------
// 6. MAGNETO-OPTIC EFFECTS  –  Faraday rotation
// ---------------------------------------------------------------------------

/// Faraday rotation parameters for a magneto-optic material.
#[derive(Debug, Clone, Copy)]
pub struct FaradayRotation {
    /// Verdet constant V \[rad/(T·m)\].
    pub verdet_constant: f64,
    /// Applied magnetic field B \[T\].
    pub b_field: f64,
    /// Interaction length L \[m\].
    pub length: f64,
}

impl FaradayRotation {
    /// Construct Faraday rotation parameters.
    pub fn new(verdet_constant: f64, b_field: f64, length: f64) -> Self {
        Self {
            verdet_constant,
            b_field,
            length,
        }
    }

    /// Rotation angle θ = V · B · L \[rad\].
    pub fn rotation_angle(&self) -> f64 {
        self.verdet_constant * self.b_field * self.length
    }

    /// Rotation angle in degrees.
    pub fn rotation_angle_deg(&self) -> f64 {
        self.rotation_angle().to_degrees()
    }

    /// Verdet constant for terbium gallium garnet (TGG) at 1064 nm \[rad/(T·m)\].
    pub fn tgg_verdet_1064nm() -> f64 {
        -40.0 // approx value
    }

    /// Off-diagonal permittivity element ε_xy for a gyromagnetic medium.
    ///
    /// ε_xy = i n Δn where Δn is the circular birefringence.
    /// Here we return a qualitative estimate.
    pub fn eps_xy_estimate(&self, n: f64) -> f64 {
        // θ ≈ π Δn L / λ  →  Δn ≈ θ λ / (π L)
        // ε_xy ≈ 2 n Δn (rough).
        let theta = self.rotation_angle();
        let lambda_approx = 1.064e-6_f64;
        let delta_n = theta * lambda_approx / (PI * self.length.max(1e-30));
        2.0 * n * delta_n
    }
}

// ---------------------------------------------------------------------------
// 7. SUPERCONDUCTOR  –  Meissner effect
// ---------------------------------------------------------------------------

/// Superconductor electromagnetic parameters.
#[derive(Debug, Clone, Copy)]
pub struct Superconductor {
    /// Critical temperature Tc \[K\].
    pub tc: f64,
    /// London penetration depth λ_L at T = 0 \[m\].
    pub lambda_london_0: f64,
    /// BCS coherence length ξ₀ \[m\].
    pub xi0: f64,
}

impl Superconductor {
    /// Construct superconductor parameters.
    pub fn new(tc: f64, lambda_london_0: f64, xi0: f64) -> Self {
        Self {
            tc,
            lambda_london_0,
            xi0,
        }
    }

    /// Typical parameters for YBCO (YBa₂Cu₃O₇).
    pub fn ybco() -> Self {
        Self {
            tc: 92.0,
            lambda_london_0: 140e-9,
            xi0: 1.2e-9,
        }
    }

    /// Typical parameters for niobium.
    pub fn niobium() -> Self {
        Self {
            tc: 9.26,
            lambda_london_0: 39e-9,
            xi0: 38e-9,
        }
    }

    /// Temperature-dependent London penetration depth \[m\] using the
    /// two-fluid model: λ(T) = λ_L(0) / √(1 − (T/Tc)⁴).
    pub fn penetration_depth(&self, temp_k: f64) -> f64 {
        if temp_k >= self.tc {
            return f64::INFINITY; // Normal state.
        }
        let t_ratio = temp_k / self.tc;
        self.lambda_london_0 / (1.0 - t_ratio.powi(4)).sqrt()
    }

    /// Ginzburg-Landau parameter κ = λ / ξ.
    ///
    /// κ > 1/√2 → type-II superconductor.
    pub fn gl_parameter(&self, temp_k: f64) -> f64 {
        self.penetration_depth(temp_k) / self.xi0
    }

    /// Surface resistance Rs(f) \[Ω\] in the Mattis-Bardeen limit
    /// (low T, f << 2Δ/h): R_s ≈ μ₀² σ_n ω² λ³_L / 2.
    pub fn surface_resistance(&self, freq: f64, sigma_normal: f64, temp_k: f64) -> f64 {
        let omega = 2.0 * PI * freq;
        let lam = self.penetration_depth(temp_k);
        if lam.is_infinite() {
            return f64::INFINITY;
        }
        MU0 * MU0 * sigma_normal * omega * omega * lam.powi(3) / 2.0
    }

    /// Whether the material is in the superconducting state at `temp_k`.
    pub fn is_superconducting(&self, temp_k: f64) -> bool {
        temp_k < self.tc
    }

    /// Lower critical field Hc1 estimate \[A/m\] (from London theory).
    pub fn hc1_estimate(&self, temp_k: f64) -> f64 {
        let lam = self.penetration_depth(temp_k);
        if lam.is_infinite() {
            return 0.0;
        }
        // Hc1 ≈ Φ₀ ln(κ) / (4π μ₀ λ²)
        let phi0 = 2.067_833_848e-15_f64; // Flux quantum [Wb]
        let kappa = self.gl_parameter(temp_k);
        if kappa <= 1.0 {
            return 0.0;
        }
        phi0 * kappa.ln() / (4.0 * PI * MU0 * lam * lam)
    }
}

// ---------------------------------------------------------------------------
// 8. FERROELECTRIC HYSTERESIS  –  P-E loop
// ---------------------------------------------------------------------------

/// Simple Preisach-inspired P-E loop model for ferroelectrics.
///
/// Uses the phenomenological tanh saturation model:
/// P(E) = P_sat · tanh((E − E_c · sign(E)) / E₀)
#[derive(Debug, Clone, Copy)]
pub struct FerroelectricPE {
    /// Saturation polarisation P_sat \[C/m²\].
    pub p_sat: f64,
    /// Coercive field E_c \[V/m\].
    pub e_coercive: f64,
    /// Slope parameter E₀ \[V/m\] (steepness of the loop wall).
    pub e0: f64,
}

impl FerroelectricPE {
    /// Construct a P-E hysteresis model.
    pub fn new(p_sat: f64, e_coercive: f64, e0: f64) -> Self {
        Self {
            p_sat,
            e_coercive,
            e0,
        }
    }

    /// Typical barium titanate parameters.
    pub fn batio3() -> Self {
        Self {
            p_sat: 0.26,
            e_coercive: 1.0e5,
            e0: 5.0e4,
        }
    }

    /// Polarisation on the ascending branch P⁺(E).
    pub fn polarisation_ascending(&self, e_field: f64) -> f64 {
        self.p_sat * ((e_field - self.e_coercive) / self.e0).tanh()
    }

    /// Polarisation on the descending branch P⁻(E).
    pub fn polarisation_descending(&self, e_field: f64) -> f64 {
        self.p_sat * ((e_field + self.e_coercive) / self.e0).tanh()
    }

    /// Remnant polarisation P_r at E = 0 on the ascending branch \[C/m²\].
    pub fn remnant_polarisation(&self) -> f64 {
        self.polarisation_ascending(0.0).abs()
    }

    /// Hysteresis loop area (energy loss per cycle per unit volume) \[J/m³\].
    ///
    /// Estimated by integrating E dP over one cycle.
    pub fn hysteresis_energy_loss(&self, e_max: f64, n_steps: usize) -> f64 {
        let de = 2.0 * e_max / n_steps as f64;
        let mut loss = 0.0_f64;
        // Ascending: E from -e_max to +e_max.
        let mut e = -e_max;
        let mut p_prev = self.polarisation_ascending(e);
        for _ in 0..n_steps {
            e += de;
            let p = self.polarisation_ascending(e);
            loss += e * (p - p_prev);
            p_prev = p;
        }
        // Descending: E from +e_max to -e_max.
        e = e_max;
        p_prev = self.polarisation_descending(e);
        for _ in 0..n_steps {
            e -= de;
            let p = self.polarisation_descending(e);
            loss -= e * (p - p_prev); // sign for descending direction.
            p_prev = p;
        }
        loss.abs()
    }
}

// ---------------------------------------------------------------------------
// 9. MAGNETIC HYSTERESIS  –  B-H loop & Jiles-Atherton model
// ---------------------------------------------------------------------------

/// Nominal RK4 steps per characteristic field scale `min(a, k)` of the JA integrator.
const JA_STEPS_PER_SCALE: f64 = 20.0;

/// Local error tolerance of the step-doubling RK4, relative to `M_s`.
const JA_RELATIVE_TOLERANCE: f64 = 1.0e-10;

/// Smallest step, as a fraction of the nominal step, before the integration is
/// abandoned.
const JA_MIN_STEP_FRACTION: f64 = 1.0e-9;

/// Smallest admissible value of the irreversible denominator
/// `(1 − c)·k − α·e` as a fraction of the pinning term `(1 − c)·k`.
const JA_POLE_FLOOR: f64 = 1.0e-6;

/// Upper bound on the integration steps (accepted and rejected) of one sweep.
const JA_MAX_STEPS: usize = 4_000_000;

/// Jiles-Atherton (JA) magnetic hysteresis model parameters.
///
/// The model splits the magnetisation into a reversible and an irreversible
/// (domain-wall pinning) part,
/// `M = (1 − c)·M_irr + c·M_an`, with
///
/// - effective field `H_e = H + α·M`,
/// - anhysteretic magnetisation `M_an = M_s·L(H_e / a)`, `L(x) = coth x − 1/x`,
/// - irreversible susceptibility
///   `dM_irr/dH = (M_an − M_irr) / (δ·k − α·(M_an − M_irr))` for
///   `δ·(M_an − M_irr) > 0` and `0` otherwise, where `δ = sign(dH)`.
///
/// Eliminating `M_irr = (M − c·M_an)/(1 − c)` and the implicit `H_e` gives the
/// explicit ODE in `(H, M)` integrated by [`JilesAthertonParams::compute_bh_loop`].
#[derive(Debug, Clone, Copy)]
pub struct JilesAthertonParams {
    /// Saturation magnetisation M_s \[A/m\].
    pub m_sat: f64,
    /// Shape parameter a \[A/m\] of the Langevin function.
    pub a: f64,
    /// Interdomain coupling coefficient α (dimensionless).
    pub alpha: f64,
    /// Pinning loss coefficient k \[A/m\].
    pub k: f64,
    /// Reversibility coefficient c (0..1).
    pub c: f64,
}

impl JilesAthertonParams {
    /// Construct JA parameters.
    pub fn new(m_sat: f64, a: f64, alpha: f64, k: f64, c: f64) -> Self {
        Self {
            m_sat,
            a,
            alpha,
            k,
            c,
        }
    }

    /// Typical soft iron parameters.
    pub fn soft_iron() -> Self {
        Self {
            m_sat: 1.6e6,
            a: 470.0,
            alpha: 1.0e-4,
            k: 400.0,
            c: 0.1,
        }
    }

    /// Checks that the parameters describe a physical, well-posed JA model.
    ///
    /// Requires finite `m_sat > 0`, `a > 0`, `k > 0`, `alpha ≥ 0`, `0 ≤ c ≤ 1`
    /// and the stability condition `c·α·m_sat / (3a) < 1` (otherwise the
    /// reversible feedback `1 − c·α·dM_an/dH_e` can vanish).  The pole of the
    /// irreversible susceptibility, `α·(M_an − M) = (1 − c)·k`, depends on the
    /// state and is detected while integrating (see
    /// [`JilesAthertonParams::try_compute_bh_loop`]).
    ///
    /// # Errors
    /// [`crate::Error::ParameterOutOfRange`] naming the offending parameter,
    /// or [`crate::Error::General`] for a non-finite value.
    pub fn validate(&self) -> Result<(), crate::Error> {
        crate::require_finite("m_sat", self.m_sat)?;
        crate::require_finite("a", self.a)?;
        crate::require_finite("alpha", self.alpha)?;
        crate::require_finite("k", self.k)?;
        crate::require_finite("c", self.c)?;
        crate::require_positive("m_sat", self.m_sat)?;
        crate::require_positive("a", self.a)?;
        crate::require_positive("k", self.k)?;
        crate::require_non_negative("alpha", self.alpha)?;
        crate::require_in_range("c", self.c, 0.0, 1.0)?;
        let alpha_limit = if self.c > 0.0 {
            3.0 * self.a / (self.c * self.m_sat)
        } else {
            f64::INFINITY
        };
        if self.alpha >= alpha_limit {
            return Err(crate::Error::ParameterOutOfRange {
                param: "alpha",
                value: self.alpha,
                min: 0.0,
                max: alpha_limit,
            });
        }
        Ok(())
    }

    /// Langevin function L(x) = coth(x) − 1/x (Taylor series near the origin,
    /// where the closed form suffers from cancellation).
    fn langevin(x: f64) -> f64 {
        if x.abs() < 0.1 {
            let x2 = x * x;
            x * (1.0 / 3.0
                - x2 * (1.0 / 45.0
                    - x2 * (2.0 / 945.0 - x2 * (1.0 / 4725.0 - x2 * 2.0 / 93_555.0))))
        } else {
            1.0 / x.tanh() - 1.0 / x
        }
    }

    /// Derivative of the Langevin function, L′(x) = 1/x² − 1/sinh²(x).
    fn langevin_derivative(x: f64) -> f64 {
        if x.abs() < 0.05 {
            let x2 = x * x;
            1.0 / 3.0
                - x2 * (1.0 / 15.0 - x2 * (2.0 / 189.0 - x2 * (1.0 / 675.0 - x2 * 2.0 / 10_395.0)))
        } else {
            let s = x.sinh();
            1.0 / (x * x) - 1.0 / (s * s)
        }
    }

    /// Effective field H_eff = H + α·M.
    pub fn h_effective(&self, h: f64, m: f64) -> f64 {
        h + self.alpha * m
    }

    /// Anhysteretic magnetisation M_an(H, M).
    pub fn anhysteretic(&self, h: f64, m: f64) -> f64 {
        let h_eff = self.h_effective(h, m);
        self.m_sat * Self::langevin(h_eff / self.a)
    }

    /// Slope dM/dH of the JA model at `(h, m)` on a branch with sweep
    /// direction `delta = ±1`, or `None` where the irreversible susceptibility
    /// has its pole.
    ///
    /// With `e = δ·(M_an − M)` the irreversible contribution is
    /// `(1 − c)·e / ((1 − c)·k − α·e)` for `e > 0` and zero otherwise (no
    /// negative susceptibility); the reversible contribution is
    /// `c·dM_an/dH_e`, and the effective-field feedback divides the sum by
    /// `1 − α·c·dM_an/dH_e`.  The denominator `(1 − c)·k − α·e` vanishes when
    /// the coupling α is too strong for the pinning (a Barkhausen-type jump
    /// that an integration in H cannot follow); such states are reported as
    /// `None`.
    fn slope(&self, h: f64, m: f64, delta: f64) -> Option<f64> {
        let x = self.h_effective(h, m) / self.a;
        let m_an = self.m_sat * Self::langevin(x);
        let dm_an_dh_eff = self.m_sat * Self::langevin_derivative(x) / self.a;
        let reversible = self.c * dm_an_dh_eff;
        let feedback = 1.0 - self.alpha * reversible;
        let drive = delta * (m_an - m);
        let irreversible = if drive > 0.0 && self.c < 1.0 {
            let pinning = (1.0 - self.c) * self.k;
            let denominator = pinning - self.alpha * drive;
            if denominator <= JA_POLE_FLOOR * pinning {
                return None;
            }
            (1.0 - self.c) * drive / denominator
        } else {
            0.0
        };
        Some((irreversible + reversible) / feedback)
    }

    /// One classical 4th-order Runge-Kutta step of signed width `dh` from
    /// `(h, m)`; M stays within ±M_s.  `None` if a stage reaches the pole.
    fn rk4_step(&self, h: f64, m: f64, dh: f64, delta: f64) -> Option<f64> {
        let k1 = self.slope(h, m, delta)?;
        let k2 = self.slope(h + 0.5 * dh, m + 0.5 * dh * k1, delta)?;
        let k3 = self.slope(h + 0.5 * dh, m + 0.5 * dh * k2, delta)?;
        let k4 = self.slope(h + dh, m + dh * k3, delta)?;
        Some((m + dh / 6.0 * (k1 + 2.0 * k2 + 2.0 * k3 + k4)).clamp(-self.m_sat, self.m_sat))
    }

    /// Integrates M(H) from `(h0, m0)` to the field `h1` (one monotone sweep)
    /// with the classical RK4 scheme and step doubling: every step is also
    /// taken as two half steps, the difference controls the step size and the
    /// result is the Richardson-extrapolated value.  The nominal step is
    /// `min(a, k) / 20`; it shrinks near kinks and near the susceptibility
    /// pole.
    ///
    /// # Errors
    /// [`crate::Error::Singularity`] when no admissible step exists (the
    /// pole of the irreversible susceptibility is reached) and
    /// [`crate::Error::NumericalDivergence`] when the step budget is spent.
    fn integrate(&self, h0: f64, m0: f64, h1: f64) -> Result<f64, crate::Error> {
        let span = h1 - h0;
        if span == 0.0 {
            return Ok(m0);
        }
        let delta = span.signum();
        let nominal = self.a.min(self.k) / JA_STEPS_PER_SCALE;
        let minimum = nominal * JA_MIN_STEP_FRACTION;
        let tolerance = JA_RELATIVE_TOLERANCE * self.m_sat;
        let mut h = h0;
        let mut m = m0;
        let mut trial = nominal.min(span.abs());
        for _ in 0..JA_MAX_STEPS {
            let remaining = (h1 - h) * delta;
            if remaining <= 0.0 {
                return Ok(m);
            }
            let width = trial.min(remaining);
            let step = delta * width;
            let attempt = self.rk4_step(h, m, step, delta).and_then(|full| {
                let half = self.rk4_step(h, m, 0.5 * step, delta)?;
                let two = self.rk4_step(h + 0.5 * step, half, 0.5 * step, delta)?;
                Some((full, two))
            });
            let Some((full, two)) = attempt else {
                trial = 0.5 * width;
                if trial < minimum {
                    return Err(self.pole_error(h, m));
                }
                continue;
            };
            let error = (two - full).abs() / 15.0;
            let factor = if error > 0.0 {
                (0.9 * (tolerance / error).powf(0.2)).clamp(0.2, 5.0)
            } else {
                5.0
            };
            if error <= tolerance {
                m = (two + (two - full) / 15.0).clamp(-self.m_sat, self.m_sat);
                h = if width >= remaining { h1 } else { h + step };
                trial = (width * factor).min(4.0 * nominal);
            } else {
                trial = width * factor;
                if trial < minimum {
                    return Err(self.pole_error(h, m));
                }
            }
        }
        Err(crate::Error::diverged(
            "jiles_atherton",
            JA_MAX_STEPS,
            ((h1 - h) * delta).abs(),
        ))
    }

    /// Error for a state at which the integration cannot proceed.
    fn pole_error(&self, h: f64, m: f64) -> crate::Error {
        crate::Error::singular(
            "jiles_atherton",
            format!(
                "no admissible integration step at H = {h:.6e} A/m, M = {m:.6e} A/m: the \
                 irreversible susceptibility has a pole (alpha·(M_an − M) reaches (1 − c)·k), \
                 i.e. the coupling alpha is too strong for the pinning k"
            ),
        )
    }

    /// Computes the major M(H) loop between ±`h_max`.
    ///
    /// Starting from the demagnetised state (H = 0, M = 0) the virgin curve is
    /// followed to +`h_max`, the field is swept down to −`h_max` and one full
    /// conditioning cycle is run so that the recorded loop is the periodic
    /// steady state.  The loop is then recorded as `2·n_steps` samples:
    /// the ascending branch (H = −`h_max` … `h_max` − ΔH) followed by the
    /// descending branch (H = `h_max` … −`h_max` + ΔH) with
    /// ΔH = 2·`h_max` / `n_steps`.
    ///
    /// # Errors
    /// Fails if [`JilesAthertonParams::validate`] rejects the parameters or
    /// `h_max` is not finite, with [`crate::Error::Singularity`] when the
    /// coupling α is so strong relative to the pinning `k` that the
    /// irreversible susceptibility has a pole inside the swept loop (the loop
    /// then contains a Barkhausen jump that cannot be integrated in H), and
    /// with [`crate::Error::NumericalDivergence`] if the step budget is spent.
    pub fn try_compute_bh_loop(
        &self,
        h_max: f64,
        n_steps: usize,
    ) -> Result<(Vec<f64>, Vec<f64>), crate::Error> {
        self.validate()?;
        crate::require_finite("h_max", h_max)?;
        if n_steps == 0 {
            return Ok((Vec::new(), Vec::new()));
        }
        let h_max = h_max.abs();
        let dh = 2.0 * h_max / n_steps as f64;

        // Virgin curve to the upper tip, descent to the lower tip, then one
        // conditioning cycle back to the lower tip.
        let mut m = self.integrate(0.0, 0.0, h_max)?;
        m = self.integrate(h_max, m, -h_max)?;
        m = self.integrate(-h_max, m, h_max)?;
        m = self.integrate(h_max, m, -h_max)?;

        let mut h_vals = Vec::with_capacity(2 * n_steps);
        let mut m_vals = Vec::with_capacity(2 * n_steps);

        // Ascending branch: -h_max -> +h_max.
        let mut h = -h_max;
        for i in 0..n_steps {
            h_vals.push(h);
            m_vals.push(m);
            let h_next = if i + 1 == n_steps {
                h_max
            } else {
                -h_max + (i + 1) as f64 * dh
            };
            m = self.integrate(h, m, h_next)?;
            h = h_next;
        }
        // Descending branch: +h_max -> -h_max.
        for i in 0..n_steps {
            h_vals.push(h);
            m_vals.push(m);
            let h_next = if i + 1 == n_steps {
                -h_max
            } else {
                h_max - (i + 1) as f64 * dh
            };
            m = self.integrate(h, m, h_next)?;
            h = h_next;
        }
        Ok((h_vals, m_vals))
    }

    /// Compute the M(H) loop for a triangular H waveform (see
    /// [`JilesAthertonParams::try_compute_bh_loop`]).
    ///
    /// Returns `(h_vals, m_vals)` over the cycle; both vectors are empty when
    /// the parameters are not physical, `h_max` is not finite or the model has
    /// no integrable loop (use [`JilesAthertonParams::try_compute_bh_loop`] to
    /// obtain the reason).
    pub fn compute_bh_loop(&self, h_max: f64, n_steps: usize) -> (Vec<f64>, Vec<f64>) {
        self.try_compute_bh_loop(h_max, n_steps).unwrap_or_default()
    }

    /// Estimate coercive field Hc from the BH loop \[A/m\].
    ///
    /// Hc is |H| where the descending branch crosses M = 0 (linear
    /// interpolation).  Returns 0 when the loop does not reach M = 0 within
    /// ±`h_max` (Hc > `h_max`) or the parameters are not physical.
    pub fn coercive_field(&self, h_max: f64, n_steps: usize) -> f64 {
        let (h_vals, m_vals) = self.compute_bh_loop(h_max, n_steps);
        let n = h_vals.len() / 2;
        if n == 0 {
            return 0.0;
        }
        // Descending branch followed by the lower tip that closes the loop.
        let h_desc = h_vals[n..].iter().chain(h_vals[..1].iter());
        let m_desc = m_vals[n..].iter().chain(m_vals[..1].iter());
        let pairs: Vec<(f64, f64)> = h_desc.copied().zip(m_desc.copied()).collect();
        for w in pairs.windows(2) {
            let (h0, m0) = w[0];
            let (h1, m1) = w[1];
            if m0 > 0.0 && m1 <= 0.0 {
                let t = m0 / (m0 - m1);
                return (h0 + t * (h1 - h0)).abs();
            }
        }
        0.0
    }

    /// Remanent magnetisation M_r \[A/m\]: M at H = 0 on the descending
    /// branch of the major loop (linear interpolation).
    ///
    /// Returns 0 if `h_max` is 0 or the parameters are not physical.
    pub fn remanent_magnetisation(&self, h_max: f64, n_steps: usize) -> f64 {
        let (h_vals, m_vals) = self.compute_bh_loop(h_max, n_steps);
        let n = h_vals.len() / 2;
        if n == 0 {
            return 0.0;
        }
        let h_desc = &h_vals[n..];
        let m_desc = &m_vals[n..];
        for (hw, mw) in h_desc.windows(2).zip(m_desc.windows(2)) {
            if hw[0] > 0.0 && hw[1] <= 0.0 {
                let t = hw[0] / (hw[0] - hw[1]);
                return mw[0] + t * (mw[1] - mw[0]);
            }
        }
        0.0
    }
}

/// Compute B = μ₀(H + M) for a JA model at field H.
pub fn b_from_jiles_atherton(_ja: &JilesAthertonParams, h: f64, m: f64) -> f64 {
    MU0 * (h + m)
}

// ---------------------------------------------------------------------------
// 10. ELECTROMAGNETIC SHIELDING EFFECTIVENESS
// ---------------------------------------------------------------------------

/// Electromagnetic shielding effectiveness of a conductive enclosure.
///
/// SE = R + A + M  (dB)
/// where R = reflection loss, A = absorption loss, M = re-reflection correction.
#[derive(Debug, Clone, Copy)]
pub struct ShieldingEffectiveness {
    /// Shield material relative permeability μ_r.
    pub mu_r: f64,
    /// Shield material conductivity σ \[S/m\].
    pub sigma: f64,
    /// Shield thickness d \[m\].
    pub thickness: f64,
}

impl ShieldingEffectiveness {
    /// Construct shielding effectiveness parameters.
    pub fn new(mu_r: f64, sigma: f64, thickness: f64) -> Self {
        Self {
            mu_r,
            sigma,
            thickness,
        }
    }

    /// Copper shield.
    pub fn copper(thickness: f64) -> Self {
        Self {
            mu_r: 1.0,
            sigma: 5.8e7,
            thickness,
        }
    }

    /// Mild steel shield.
    pub fn mild_steel(thickness: f64) -> Self {
        Self {
            mu_r: 200.0,
            sigma: 1.0e7,
            thickness,
        }
    }

    /// Skin depth δ at frequency f \[Hz\].
    pub fn skin_depth(&self, f: f64) -> f64 {
        1.0 / (PI * f * self.mu_r * MU0 * self.sigma).sqrt()
    }

    /// Absorption loss A \[dB\] = 20 log₁₀(e^(d/δ)).
    pub fn absorption_loss_db(&self, f: f64) -> f64 {
        let delta = self.skin_depth(f);
        if delta < 1e-30 {
            return f64::INFINITY;
        }
        20.0 * (self.thickness / delta) * std::f64::consts::LOG10_E
    }

    /// Reflection loss R \[dB\] for plane-wave incidence in free space.
    ///
    /// R ≈ 168 + 10 log₁₀(σ / (μ_r f)) \[dB\].
    pub fn reflection_loss_db(&self, f: f64) -> f64 {
        168.0 + 10.0 * (self.sigma / (self.mu_r * f)).log10()
    }

    /// Re-reflection correction M \[dB\] (only significant when A < 10 dB).
    pub fn re_reflection_correction_db(&self, f: f64) -> f64 {
        let a = self.absorption_loss_db(f);
        if a > 10.0 {
            return 0.0;
        }
        // M ≈ 20 log₁₀|1 − 10^(−A/10)| (simplified).
        let ratio = 10.0_f64.powf(-a / 10.0);
        20.0 * (1.0 - ratio).abs().log10()
    }

    /// Total shielding effectiveness SE \[dB\] = R + A + M.
    pub fn total_se_db(&self, f: f64) -> f64 {
        self.reflection_loss_db(f)
            + self.absorption_loss_db(f)
            + self.re_reflection_correction_db(f)
    }
}

// ---------------------------------------------------------------------------
// Convenience constructors for common materials
// ---------------------------------------------------------------------------

/// Return electromagnetic properties for FR4 PCB substrate.
pub fn fr4() -> ElectromagneticMaterial {
    ElectromagneticMaterial::dielectric("FR4", 4.5, 0.02, 1.0e9)
}

/// Return electromagnetic properties for PTFE (Teflon).
pub fn ptfe() -> ElectromagneticMaterial {
    ElectromagneticMaterial::dielectric("PTFE", 2.1, 0.0002, 1.0e9)
}

/// Return electromagnetic properties for distilled water at 2.45 GHz.
pub fn water_2_45ghz() -> ElectromagneticMaterial {
    ElectromagneticMaterial::dielectric("Water (2.45 GHz)", 80.0, 0.157, 2.45e9)
}

/// Return electromagnetic properties for copper.
pub fn copper_em() -> ElectromagneticMaterial {
    ElectromagneticMaterial::conductor("Copper", 5.8e7)
}

// ---------------------------------------------------------------------------
// Unit tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-9;

    // ---- Tensor3x3 ---------------------------------------------------------

    #[test]
    fn test_tensor_diagonal_trace() {
        let t = Tensor3x3::diagonal(3.0);
        assert!((t.trace() - 9.0).abs() < EPS);
    }

    #[test]
    fn test_tensor_from_diag_mul_vec() {
        let t = Tensor3x3::from_diag(2.0, 3.0, 4.0);
        let y = t.mul_vec([1.0, 1.0, 1.0]);
        assert!((y[0] - 2.0).abs() < EPS);
        assert!((y[1] - 3.0).abs() < EPS);
        assert!((y[2] - 4.0).abs() < EPS);
    }

    #[test]
    fn test_tensor_transpose_diagonal_unchanged() {
        let t = Tensor3x3::diagonal(5.0);
        let tr = t.transpose();
        for i in 0..3 {
            for j in 0..3 {
                assert!((t.data[i][j] - tr.data[i][j]).abs() < EPS);
            }
        }
    }

    #[test]
    fn test_tensor_is_symmetric() {
        let t = Tensor3x3::diagonal(1.0);
        assert!(t.is_symmetric(1e-12));
    }

    #[test]
    fn test_tensor_zero_trace() {
        let t = Tensor3x3::zero();
        assert_eq!(t.trace(), 0.0);
    }

    // ---- ElectromagneticMaterial -------------------------------------------

    #[test]
    fn test_em_material_eps_r_scalar() {
        let m = ElectromagneticMaterial::dielectric("test", 4.0, 0.01, 1e9);
        assert!((m.eps_r_scalar() - 4.0).abs() < EPS);
    }

    #[test]
    fn test_em_material_intrinsic_impedance_free_space() {
        // Vacuum: eps_r = 1, mu_r = 1 → Z = Z0.
        let m = ElectromagneticMaterial::dielectric("vacuum", 1.0, 0.0, 1e9);
        assert!((m.intrinsic_impedance() - Z0).abs() < 1e-6);
    }

    #[test]
    fn test_em_material_phase_velocity_free_space() {
        let m = ElectromagneticMaterial::dielectric("vacuum", 1.0, 0.0, 1e9);
        assert!((m.phase_velocity() - C0).abs() < 1.0); // within 1 m/s
    }

    #[test]
    fn test_em_material_skin_depth_conductor() {
        let m = copper_em();
        let d = m.skin_depth(1.0e9); // 1 GHz
        assert!(d > 0.0 && d < 1e-3, "skin depth = {d}");
    }

    #[test]
    fn test_em_material_insulator_infinite_skin_depth() {
        let m = ptfe();
        let d = m.skin_depth(1.0e9);
        assert!(d.is_infinite());
    }

    // ---- Metamaterial ------------------------------------------------------

    #[test]
    fn test_metamaterial_double_negative() {
        let m = Metamaterial {
            name: "DNG".to_owned(),
            eps_r_real: -1.0,
            eps_r_imag: 0.01,
            mu_r_real: -1.0,
            mu_r_imag: 0.01,
            freq: 1.0e10,
        };
        assert_eq!(m.classify(), MetamaterialClass::DoubleNegative);
        assert!(m.refractive_index() < 0.0);
    }

    #[test]
    fn test_metamaterial_enz() {
        let m = Metamaterial {
            name: "ENZ".to_owned(),
            eps_r_real: 0.005,
            eps_r_imag: 0.001,
            mu_r_real: 1.0,
            mu_r_imag: 0.0,
            freq: 1.0e12,
        };
        assert_eq!(m.classify(), MetamaterialClass::EpsilonNearZero);
    }

    #[test]
    fn test_metamaterial_eps_negative() {
        let m = Metamaterial {
            name: "ENG".to_owned(),
            eps_r_real: -2.0,
            eps_r_imag: 0.1,
            mu_r_real: 1.0,
            mu_r_imag: 0.0,
            freq: 1.0e10,
        };
        assert_eq!(m.classify(), MetamaterialClass::EpsilonNegative);
    }

    // ---- Drude model -------------------------------------------------------

    #[test]
    fn test_drude_gold_eps_negative_visible() {
        let d = DrudeModel::gold();
        // Gold at ω = 2π × 500 THz (green light).
        let omega = 2.0 * PI * 500.0e12;
        assert!(
            d.eps_real(omega) < 0.0,
            "gold eps_real should be negative in visible"
        );
    }

    #[test]
    fn test_drude_dc_conductivity_positive() {
        let d = DrudeModel::silver();
        assert!(d.dc_conductivity() > 0.0);
    }

    #[test]
    fn test_drude_plasma_freq_hz() {
        let d = DrudeModel::gold();
        let fp = d.plasma_freq_hz();
        assert!(fp > 1.0e14, "plasma freq should be >100 THz for gold: {fp}");
    }

    #[test]
    fn test_drude_skin_depth_positive() {
        let d = DrudeModel::silver();
        let omega = 2.0 * PI * 1.0e15;
        let delta = d.skin_depth(omega);
        assert!(delta > 0.0);
    }

    // ---- Lorentz oscillator ------------------------------------------------

    #[test]
    fn test_lorentz_static_eps() {
        let l = LorentzOscillator::new(2.0, 3.0, 1.0e15, 1.0e13);
        assert!((l.eps_static() - 5.0).abs() < EPS);
    }

    #[test]
    fn test_lorentz_eps_imag_positive() {
        let l = LorentzOscillator::new(2.0, 3.0, 1.0e15, 1.0e13);
        assert!(l.eps_imag(1.0e15) > 0.0);
    }

    #[test]
    fn test_lorentz_eps_real_below_resonance() {
        // Far below resonance eps_real → eps_inf + delta_eps = eps_static.
        let l = LorentzOscillator::new(2.0, 3.0, 1.0e15, 1.0e13);
        let er_low = l.eps_real(1.0e10); // much below resonance
        assert!((er_low - l.eps_static()).abs() < 0.1);
    }

    #[test]
    fn test_lorentz_resonance_peak() {
        let l = LorentzOscillator::new(2.0, 3.0, 1.0e15, 1.0e13);
        let peak = l.eps_imag_at_resonance();
        assert!(peak > 0.0);
    }

    // ---- Sellmeier model ---------------------------------------------------

    #[test]
    fn test_sellmeier_bk7_visible() {
        let s = SellmeierModel::bk7();
        let n = s.refractive_index(0.587); // 587 nm
        assert!(n > 1.5 && n < 1.6, "BK7 n at 587 nm = {n}");
    }

    #[test]
    fn test_sellmeier_fused_silica() {
        let s = SellmeierModel::fused_silica();
        let n = s.refractive_index(1.0); // 1000 nm
        assert!(n > 1.44 && n < 1.46, "SiO2 n at 1 μm = {n}");
    }

    #[test]
    fn test_sellmeier_group_index_greater_than_phase_index() {
        let s = SellmeierModel::bk7();
        let n = s.refractive_index(0.8);
        let ng = s.group_index(0.8);
        // In normal dispersion, n_g > n.
        assert!(ng > n, "n_g={ng}, n={n}");
    }

    #[test]
    fn test_sellmeier_gvd_finite() {
        let s = SellmeierModel::fused_silica();
        let gvd = s.gvd(1.3); // near zero-dispersion wavelength
        assert!(gvd.is_finite());
    }

    // ---- Faraday rotation --------------------------------------------------

    #[test]
    fn test_faraday_rotation_angle() {
        let f = FaradayRotation::new(100.0, 1.0, 0.01); // V=100 rad/(T·m), B=1 T, L=1 cm
        let theta = f.rotation_angle_deg();
        assert!((theta - 100.0_f64.to_degrees() * 0.01).abs() < 1e-6);
    }

    #[test]
    fn test_faraday_rotation_zero_field() {
        let f = FaradayRotation::new(100.0, 0.0, 0.01);
        assert!((f.rotation_angle()).abs() < EPS);
    }

    #[test]
    fn test_faraday_rotation_eps_xy_finite() {
        let f = FaradayRotation::new(40.0, 0.5, 0.05);
        let eps_xy = f.eps_xy_estimate(2.3);
        assert!(eps_xy.is_finite());
    }

    // ---- Superconductor ----------------------------------------------------

    #[test]
    fn test_superconductor_below_tc() {
        let s = Superconductor::ybco();
        assert!(s.is_superconducting(77.0));
        assert!(!s.is_superconducting(100.0));
    }

    #[test]
    fn test_superconductor_penetration_depth_increases_with_temp() {
        let s = Superconductor::niobium();
        let d_low = s.penetration_depth(1.0);
        let d_high = s.penetration_depth(8.0);
        assert!(d_high > d_low, "λ should increase with temperature");
    }

    #[test]
    fn test_superconductor_gl_parameter_type_ii() {
        let s = Superconductor::ybco();
        let kappa = s.gl_parameter(77.0);
        assert!(
            kappa > std::f64::consts::FRAC_1_SQRT_2,
            "YBCO should be type-II (κ={kappa})"
        );
    }

    #[test]
    fn test_superconductor_hc1_positive() {
        let s = Superconductor::ybco();
        let hc1 = s.hc1_estimate(77.0);
        assert!(hc1 > 0.0);
    }

    // ---- Ferroelectric hysteresis ------------------------------------------

    #[test]
    fn test_ferroelectric_remnant_polarisation_positive() {
        let fe = FerroelectricPE::batio3();
        assert!(fe.remnant_polarisation() >= 0.0);
    }

    #[test]
    fn test_ferroelectric_saturation_at_large_field() {
        let fe = FerroelectricPE::batio3();
        let p_big = fe.polarisation_ascending(1.0e7);
        assert!(
            (p_big - fe.p_sat).abs() < 1e-4 * fe.p_sat,
            "should be near saturation at large E"
        );
    }

    #[test]
    fn test_ferroelectric_hysteresis_energy_positive() {
        let fe = FerroelectricPE::batio3();
        let e_loss = fe.hysteresis_energy_loss(5.0e5, 100);
        assert!(e_loss >= 0.0);
    }

    #[test]
    fn test_ferroelectric_ascending_descending_differ() {
        let fe = FerroelectricPE::batio3();
        let pa = fe.polarisation_ascending(0.0);
        let pd = fe.polarisation_descending(0.0);
        assert!(
            (pa - pd).abs() > 1e-10,
            "ascending and descending should differ at E=0"
        );
    }

    // ---- Jiles-Atherton model ----------------------------------------------

    #[test]
    fn test_ja_langevin_small_x() {
        let ja = JilesAthertonParams::soft_iron();
        // Anhysteretic at H=0, M=0 should be 0.
        let man = ja.anhysteretic(0.0, 0.0);
        assert!(man.abs() < 1e-10, "M_an(0) should be ~0, got {man}");
    }

    #[test]
    fn test_ja_langevin_is_continuous_across_series_switch() {
        for x0 in [0.1, 0.05] {
            let below = JilesAthertonParams::langevin(x0 - 1.0e-9);
            let above = JilesAthertonParams::langevin(x0 + 1.0e-9);
            assert!((below - above).abs() < 1.0e-8, "L jumps at {x0}");
            let d_below = JilesAthertonParams::langevin_derivative(x0 - 1.0e-9);
            let d_above = JilesAthertonParams::langevin_derivative(x0 + 1.0e-9);
            assert!((d_below - d_above).abs() < 1.0e-8, "L' jumps at {x0}");
        }
        // L(x) -> 1 - 1/x for large x and L is odd.
        let big = JilesAthertonParams::langevin(200.0);
        assert!((big - (1.0 - 1.0 / 200.0)).abs() < 1.0e-12);
        let l = JilesAthertonParams::langevin(0.37);
        assert!((l + JilesAthertonParams::langevin(-0.37)).abs() < 1.0e-15);
    }

    #[test]
    fn test_ja_langevin_derivative_matches_finite_difference() {
        for x in [0.01_f64, 0.049, 0.051, 0.3, 1.0, 2.5, 8.0, 25.0] {
            let h = 1.0e-6 * x.max(1.0);
            let numeric = (JilesAthertonParams::langevin(x + h)
                - JilesAthertonParams::langevin(x - h))
                / (2.0 * h);
            let analytic = JilesAthertonParams::langevin_derivative(x);
            assert!(
                (numeric - analytic).abs() <= 1.0e-6 * analytic.abs().max(1.0e-3),
                "L'({x}): analytic {analytic} vs numeric {numeric}"
            );
        }
        assert!((JilesAthertonParams::langevin_derivative(0.0) - 1.0 / 3.0).abs() < 1.0e-15);
    }

    #[test]
    fn test_ja_compute_bh_loop_length() {
        let ja = JilesAthertonParams::soft_iron();
        let (h, m) = ja.compute_bh_loop(2000.0, 50);
        assert_eq!(h.len(), 100);
        assert_eq!(m.len(), 100);
    }

    #[test]
    fn test_ja_loop_sampling_layout() {
        // Ascending branch -h_max .. h_max - dh, then descending h_max .. -h_max + dh.
        let ja = JilesAthertonParams::soft_iron();
        let (n, h_max) = (50_usize, 2000.0);
        let dh = 2.0 * h_max / n as f64;
        let (h, _m) = ja.compute_bh_loop(h_max, n);
        assert!((h[0] + h_max).abs() < 1.0e-9);
        assert!((h[n - 1] - (h_max - dh)).abs() < 1.0e-9);
        assert!((h[n] - h_max).abs() < 1.0e-9);
        assert!((h[2 * n - 1] - (-h_max + dh)).abs() < 1.0e-9);
        for i in 1..n {
            assert!(h[i] > h[i - 1] && h[n + i] < h[n + i - 1]);
        }
    }

    #[test]
    fn test_ja_magnetisation_bounded() {
        let ja = JilesAthertonParams::soft_iron();
        let (_h, m) = ja.compute_bh_loop(5000.0, 100);
        for mi in &m {
            assert!(mi.abs() <= ja.m_sat * 1.01, "M exceeds M_sat: {mi}");
        }
    }

    #[test]
    fn test_ja_loop_is_odd_symmetric() {
        // The periodic loop satisfies M_asc(H) = -M_desc(-H); in particular
        // M(+H_max) = -M(-H_max).  (The old integrator ran M up to +M_sat on
        // the descending branch.)
        let ja = JilesAthertonParams::soft_iron();
        let n = 100;
        let (h, m) = ja.compute_bh_loop(5000.0, n);
        assert!(
            (m[n] + m[0]).abs() < 1.0e-6 * ja.m_sat,
            "tips: {} {}",
            m[n],
            m[0]
        );
        for i in 0..n {
            assert!((h[i] + h[n + i]).abs() < 1.0e-9);
            assert!(
                (m[i] + m[n + i]).abs() < 1.0e-6 * ja.m_sat,
                "M asc({}) = {} vs -M desc({}) = {}",
                h[i],
                m[i],
                h[n + i],
                -m[n + i]
            );
        }
        // Tip magnetisation lags just below the anhysteretic value.
        let m_tip = m[n];
        let m_an = ja.anhysteretic(5000.0, m_tip);
        assert!(
            m_tip > 0.95 * m_an && m_tip < m_an,
            "{m_tip} vs M_an = {m_an}"
        );
        assert!(
            m[0] < -0.8 * ja.m_sat,
            "lower tip must be negative: {}",
            m[0]
        );
    }

    #[test]
    fn test_ja_branches_are_monotone_and_ordered() {
        let ja = JilesAthertonParams::soft_iron();
        let n = 100;
        let (_h, m) = ja.compute_bh_loop(3000.0, n);
        for i in 1..n {
            assert!(m[i] >= m[i - 1], "ascending M must not decrease at {i}");
            assert!(
                m[n + i] <= m[n + i - 1],
                "descending M must not increase at {i}"
            );
        }
        // At equal H the descending branch lies above the ascending one.
        for i in 1..n {
            let ascending = m[i];
            let descending = m[2 * n - i];
            assert!(
                descending >= ascending - 1.0e-9 * ja.m_sat,
                "H = {}: descending {descending} below ascending {ascending}",
                -3000.0 + i as f64 * 60.0
            );
        }
    }

    #[test]
    fn test_ja_coercive_field_positive() {
        let ja = JilesAthertonParams::soft_iron();
        let hc = ja.coercive_field(2000.0, 100);
        assert!(hc > 0.0, "coercive field must be positive");
        // Hc < (1 - c) k (the thin-loop value) and of the order of k.
        assert!(hc > 300.0 && hc < 0.9 * ja.k, "Hc = {hc}");
    }

    #[test]
    fn test_ja_remanence_positive_and_bounded() {
        let ja = JilesAthertonParams::soft_iron();
        let br = ja.remanent_magnetisation(3000.0, 200);
        assert!(br > 0.2 * ja.m_sat && br < 0.3 * ja.m_sat, "Br = {br}");
        // Remanence on the ascending branch is -Br (odd symmetry).
        let n = 200;
        let (h, m) = ja.compute_bh_loop(3000.0, n);
        let i_zero = h[..n]
            .iter()
            .position(|&hv| hv.abs() < 1.0e-9)
            .expect("H = 0 is a sample for even n");
        assert!(
            (m[i_zero] + br).abs() < 1.0e-3 * br,
            "{} vs {br}",
            m[i_zero]
        );
    }

    #[test]
    fn test_ja_thin_loop_coercivity_matches_analytic_limit() {
        // Without interdomain coupling and for k << a the descending branch
        // lags the anhysteretic curve by (1 - c) k, so Hc = (1 - c) k.
        for (c, k) in [(0.1, 5.0), (0.3, 5.0)] {
            let ja = JilesAthertonParams::new(1.6e6, 470.0, 0.0, k, c);
            let hc = ja.coercive_field(3000.0, 200);
            let expected = (1.0 - c) * k;
            assert!(
                (hc - expected).abs() < 0.01 * expected,
                "c = {c}: Hc = {hc}, expected {expected}"
            );
        }
    }

    fn ja_loop_area(h: &[f64], m: &[f64]) -> f64 {
        // -oint M dH over the closed loop (ascending, then descending branch).
        let n = h.len();
        let mut area = 0.0;
        for i in 0..n {
            let j = (i + 1) % n;
            area -= 0.5 * (m[i] + m[j]) * (h[j] - h[i]);
        }
        area
    }

    #[test]
    fn test_ja_loss_and_coercivity_increase_with_pinning() {
        let mut previous_hc = 0.0;
        let mut previous_area = 0.0;
        for k in [100.0, 200.0, 400.0, 800.0] {
            let ja = JilesAthertonParams::new(1.6e6, 470.0, 1.0e-4, k, 0.1);
            let hc = ja.coercive_field(3000.0, 200);
            let (h, m) = ja.compute_bh_loop(3000.0, 200);
            let area = ja_loop_area(&h, &m);
            assert!(hc > previous_hc, "Hc must grow with k: k = {k}, Hc = {hc}");
            assert!(area > previous_area, "loop area must grow with k: k = {k}");
            previous_hc = hc;
            previous_area = area;
        }
        assert!(previous_area > 0.0);
    }

    #[test]
    fn test_ja_fully_reversible_material_has_no_hysteresis() {
        // c = 1: M follows the anhysteretic curve M = M_an(H + alpha M).
        let ja = JilesAthertonParams::new(1.6e6, 470.0, 1.0e-4, 400.0, 1.0);
        let n = 200;
        let (h, m) = ja.compute_bh_loop(3000.0, n);
        for i in 1..n {
            assert!(
                (m[i] - m[2 * n - i]).abs() < 1.0e-8 * ja.m_sat,
                "branches differ at H = {}",
                h[i]
            );
        }
        for (hv, mv) in h.iter().zip(m.iter()) {
            let m_an = ja.anhysteretic(*hv, *mv);
            assert!(
                (mv - m_an).abs() < 1.0e-7 * ja.m_sat,
                "H = {hv}: M = {mv}, M_an = {m_an}"
            );
        }
        assert!(ja.coercive_field(3000.0, n) < 1.0e-3);
    }

    #[test]
    fn test_ja_loop_is_independent_of_sampling_resolution() {
        let ja = JilesAthertonParams::soft_iron();
        let coarse = ja.compute_bh_loop(3000.0, 50);
        let fine = ja.compute_bh_loop(3000.0, 400);
        assert!((coarse.1[50] - fine.1[400]).abs() < 1.0e-6 * ja.m_sat);
        let hc_coarse = ja.coercive_field(3000.0, 50);
        let hc_fine = ja.coercive_field(3000.0, 400);
        assert!(
            (hc_coarse - hc_fine).abs() < 0.005 * hc_fine,
            "{hc_coarse} vs {hc_fine}"
        );
    }

    #[test]
    fn test_ja_invalid_parameters_are_rejected() {
        let good = JilesAthertonParams::soft_iron();
        assert!(good.validate().is_ok());
        let bad_cases = [
            JilesAthertonParams::new(0.0, 470.0, 1.0e-4, 400.0, 0.1),
            JilesAthertonParams::new(1.6e6, 0.0, 1.0e-4, 400.0, 0.1),
            JilesAthertonParams::new(1.6e6, 470.0, -1.0e-4, 400.0, 0.1),
            JilesAthertonParams::new(1.6e6, 470.0, 1.0e-4, 0.0, 0.1),
            JilesAthertonParams::new(1.6e6, 470.0, 1.0e-4, 400.0, 1.5),
            JilesAthertonParams::new(1.6e6, 470.0, 1.0e-4, 400.0, -0.1),
            JilesAthertonParams::new(1.6e6, f64::NAN, 1.0e-4, 400.0, 0.1),
            // c alpha M_s / (3a) = 11.3 > 1: unstable reversible feedback.
            JilesAthertonParams::new(1.6e6, 470.0, 0.1, 400.0, 0.1),
        ];
        for ja in &bad_cases {
            assert!(ja.validate().is_err(), "{ja:?} must be rejected");
            assert!(ja.try_compute_bh_loop(1000.0, 10).is_err());
            let (h, m) = ja.compute_bh_loop(1000.0, 10);
            assert!(h.is_empty() && m.is_empty());
            assert_eq!(ja.coercive_field(1000.0, 10), 0.0);
        }
        assert!(good.try_compute_bh_loop(f64::NAN, 10).is_err());
        let (h, m) = good
            .try_compute_bh_loop(1000.0, 0)
            .expect("n_steps = 0 is valid");
        assert!(h.is_empty() && m.is_empty());
    }

    #[test]
    fn test_ja_strong_but_admissible_coupling_is_integrated_accurately() {
        // alpha M_s is several times k here and the irreversible denominator
        // (1 - c) k - alpha e comes within 6 % (alpha = 9e-4) or 1 % (9.5e-4) of
        // its pole, where a fixed-step RK4 overshoots into the M_s clamp.
        for (alpha, m_tip_expected) in [(9.0e-4, 1_404_627.98), (9.5e-4, 1_408_373.40)] {
            let ja = JilesAthertonParams::new(1.6e6, 470.0, alpha, 400.0, 0.1);
            let n = 100;
            let (h, m) = ja
                .try_compute_bh_loop(3000.0, n)
                .expect("loop must be integrable");
            assert!(
                (m[n] - m_tip_expected).abs() < 1.0e-4 * m_tip_expected,
                "alpha = {alpha}: M_tip = {} vs {m_tip_expected}",
                m[n]
            );
            // The tip lies below the anhysteretic curve (M lags M_an).
            assert!(m[n] < ja.anhysteretic(3000.0, m[n]), "alpha = {alpha}");
            for i in 0..n {
                assert!((h[i] + h[n + i]).abs() < 1.0e-9);
                assert!(
                    (m[i] + m[n + i]).abs() < 1.0e-6 * ja.m_sat,
                    "alpha = {alpha}: odd symmetry broken at H = {}",
                    h[i]
                );
            }
            for i in 1..n {
                assert!(m[i] >= m[i - 1] && m[n + i] <= m[n + i - 1]);
                assert!(m[2 * n - i] >= m[i], "alpha = {alpha}: branch order at {i}");
            }
        }
    }

    #[test]
    fn test_ja_coupling_beyond_the_susceptibility_pole_is_reported() {
        // alpha = 1e-3 makes alpha (M_an - M) exceed (1 - c) k on the major loop:
        // the model then has a Barkhausen jump that cannot be integrated in H.
        let strong = JilesAthertonParams::new(1.6e6, 470.0, 1.0e-3, 400.0, 0.1);
        assert!(
            strong.validate().is_ok(),
            "the parameters are individually valid"
        );
        assert!(matches!(
            strong.try_compute_bh_loop(1000.0, 100),
            Err(crate::Error::Singularity { .. })
        ));
        let (h, m) = strong.compute_bh_loop(1000.0, 100);
        assert!(h.is_empty() && m.is_empty());
        assert_eq!(strong.coercive_field(1000.0, 100), 0.0);
        assert_eq!(strong.remanent_magnetisation(1000.0, 100), 0.0);
        // A small minor loop never reaches the pole and is still computed.
        let (h, m) = strong
            .try_compute_bh_loop(100.0, 100)
            .expect("minor loop is integrable");
        assert_eq!(h.len(), 200);
        assert!(m.iter().all(|v| v.is_finite() && v.abs() < strong.m_sat));
        // Strongly coupled, partly reversible material.
        let steel = JilesAthertonParams::new(1.6e6, 470.0, 1.6e-3, 400.0, 0.2);
        assert!(matches!(
            steel.try_compute_bh_loop(500.0, 50),
            Err(crate::Error::Singularity { .. })
        ));
    }

    #[test]
    fn test_b_from_ja_positive_at_positive_h() {
        let ja = JilesAthertonParams::soft_iron();
        let b = b_from_jiles_atherton(&ja, 1000.0, 0.5e6);
        assert!(b > 0.0);
    }

    // ---- Shielding effectiveness -------------------------------------------

    #[test]
    fn test_shielding_skin_depth_copper() {
        let s = ShieldingEffectiveness::copper(1e-3);
        let d = s.skin_depth(1.0e6); // 1 MHz
        assert!(d > 0.0 && d < 1e-3, "copper skin depth at 1 MHz: {d}");
    }

    #[test]
    fn test_shielding_absorption_loss_increases_with_freq() {
        let s = ShieldingEffectiveness::copper(1e-3);
        let a1 = s.absorption_loss_db(1.0e6);
        let a2 = s.absorption_loss_db(1.0e9);
        assert!(a2 > a1, "absorption loss should increase with frequency");
    }

    #[test]
    fn test_shielding_total_se_positive() {
        let s = ShieldingEffectiveness::mild_steel(2e-3);
        let se = s.total_se_db(1.0e6);
        assert!(se > 0.0, "SE should be positive: {se}");
    }

    #[test]
    fn test_shielding_reflection_loss_copper_1mhz() {
        let s = ShieldingEffectiveness::copper(1e-3);
        let r = s.reflection_loss_db(1.0e6);
        // Copper at 1 MHz has high reflection loss (> 100 dB).
        assert!(r > 100.0, "R = {r} dB");
    }

    #[test]
    fn test_fr4_eps_r() {
        let m = fr4();
        assert!((m.eps_r_scalar() - 4.5).abs() < EPS);
    }

    #[test]
    fn test_water_high_permittivity() {
        let m = water_2_45ghz();
        assert!(m.eps_r_scalar() > 70.0);
    }
}
