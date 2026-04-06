#![allow(clippy::needless_range_loop)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Elastic wave propagation in soft bodies.
//!
//! Provides wave equation modelling including P-waves, S-waves, Rayleigh waves,
//! Love waves, dispersion relations, Zoeppritz reflection/transmission coefficients,
//! standing waves, 1D FEM wave solver, attenuation models, and phononic crystal
//! band-gap analysis.

// ---------------------------------------------------------------------------
// WaveMode
// ---------------------------------------------------------------------------

/// Elastic wave mode classification.
///
/// Distinguishes the four principal elastic wave types encountered in
/// soft-body and geomechanical simulations.
#[allow(dead_code)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum WaveMode {
    /// Compressional (P) wave — particle motion parallel to propagation.
    Longitudinal,
    /// Shear (S) wave — particle motion perpendicular to propagation.
    Transverse,
    /// Rayleigh surface wave — retrograde elliptical particle motion.
    Rayleigh,
    /// Love surface wave — horizontally polarised shear.
    Love,
}

// ---------------------------------------------------------------------------
// ElasticWaveParams
// ---------------------------------------------------------------------------

/// Material parameters for elastic wave propagation.
///
/// Encapsulates the Lamé constants and density required to compute
/// wave speeds and acoustic impedances.
#[allow(dead_code)]
#[derive(Debug, Clone, Copy)]
pub struct ElasticWaveParams {
    /// Material density (kg/m³).
    pub density: f64,
    /// First Lamé parameter λ (Pa).
    pub lame_lambda: f64,
    /// Shear modulus μ (Pa).
    pub lame_mu: f64,
}

impl ElasticWaveParams {
    /// Create new elastic wave parameters.
    ///
    /// # Arguments
    /// * `density` - Material density in kg/m³.
    /// * `lame_lambda` - First Lamé constant in Pa.
    /// * `lame_mu` - Shear modulus in Pa.
    pub fn new(density: f64, lame_lambda: f64, lame_mu: f64) -> Self {
        Self {
            density,
            lame_lambda,
            lame_mu,
        }
    }

    /// Compute the P-wave (compressional) speed.
    ///
    /// c_P = sqrt((λ + 2μ) / ρ)
    pub fn p_wave_speed(&self) -> f64 {
        ((self.lame_lambda + 2.0 * self.lame_mu) / self.density).sqrt()
    }

    /// Compute the S-wave (shear) speed.
    ///
    /// c_S = sqrt(μ / ρ)
    pub fn s_wave_speed(&self) -> f64 {
        (self.lame_mu / self.density).sqrt()
    }

    /// Estimate the Rayleigh wave speed (approximately 0.92–0.95 × c_S).
    ///
    /// Uses the approximation c_R ≈ c_S * (0.862 + 1.14ν) / (1 + ν)
    /// where ν is Poisson's ratio.
    pub fn rayleigh_wave_speed(&self) -> f64 {
        let cs = self.s_wave_speed();
        let nu = self.poissons_ratio();
        cs * (0.862 + 1.14 * nu) / (1.0 + nu)
    }

    /// Compute the acoustic impedance for P-waves.
    ///
    /// Z_P = ρ * c_P
    pub fn impedance_p(&self) -> f64 {
        self.density * self.p_wave_speed()
    }

    /// Compute the acoustic impedance for S-waves.
    ///
    /// Z_S = ρ * c_S
    pub fn impedance_s(&self) -> f64 {
        self.density * self.s_wave_speed()
    }

    /// Compute Young's modulus E from Lamé parameters.
    ///
    /// E = μ(3λ + 2μ) / (λ + μ)
    pub fn youngs_modulus(&self) -> f64 {
        let mu = self.lame_mu;
        let lambda = self.lame_lambda;
        mu * (3.0 * lambda + 2.0 * mu) / (lambda + mu)
    }

    /// Compute Poisson's ratio ν from Lamé parameters.
    ///
    /// ν = λ / (2(λ + μ))
    pub fn poissons_ratio(&self) -> f64 {
        let lambda = self.lame_lambda;
        let mu = self.lame_mu;
        lambda / (2.0 * (lambda + mu))
    }

    /// Compute the bulk modulus K.
    ///
    /// K = λ + 2μ/3
    pub fn bulk_modulus(&self) -> f64 {
        self.lame_lambda + 2.0 * self.lame_mu / 3.0
    }
}

// ---------------------------------------------------------------------------
// WavePacket
// ---------------------------------------------------------------------------

/// A localised wave packet with amplitude, frequency, wavenumber, and velocities.
///
/// Models a Gaussian-modulated sinusoidal wave packet for wave propagation
/// studies.
#[allow(dead_code)]
pub struct WavePacket {
    /// Wave amplitude (metres or Pa depending on context).
    pub amplitude: f64,
    /// Angular frequency ω (rad/s).
    pub frequency: f64,
    /// Wavenumber k (rad/m).
    pub wavenumber: f64,
    /// Phase velocity v_ph = ω/k (m/s).
    pub phase_velocity: f64,
    /// Group velocity v_g = dω/dk (m/s).
    pub group_velocity: f64,
    /// Packet envelope width (metres).
    pub envelope_width: f64,
    /// Initial position of packet centre (metres).
    pub center_position: f64,
}

impl WavePacket {
    /// Create a new wave packet.
    ///
    /// # Arguments
    /// * `amplitude` - Wave amplitude.
    /// * `frequency` - Angular frequency (rad/s).
    /// * `wavenumber` - Wavenumber (rad/m).
    /// * `group_velocity` - Group velocity (m/s).
    /// * `envelope_width` - Gaussian envelope width (metres).
    /// * `center_position` - Initial envelope centre (metres).
    pub fn new(
        amplitude: f64,
        frequency: f64,
        wavenumber: f64,
        group_velocity: f64,
        envelope_width: f64,
        center_position: f64,
    ) -> Self {
        let phase_velocity = if wavenumber.abs() > 1e-30 {
            frequency / wavenumber
        } else {
            0.0
        };
        Self {
            amplitude,
            frequency,
            wavenumber,
            phase_velocity,
            group_velocity,
            envelope_width,
            center_position,
        }
    }

    /// Evaluate the wave packet displacement at position x and time t.
    ///
    /// Uses a Gaussian envelope: u(x, t) = A * exp(-((x - x_c - v_g*t)²) / (2σ²)) * cos(kx - ωt)
    ///
    /// # Arguments
    /// * `x` - Spatial coordinate (metres).
    /// * `t` - Time (seconds).
    pub fn evaluate_at(&self, x: f64, t: f64) -> f64 {
        let envelope_center = self.center_position + self.group_velocity * t;
        let arg = (x - envelope_center) / self.envelope_width;
        let envelope = (-0.5 * arg * arg).exp();
        let phase = self.wavenumber * x - self.frequency * t;
        self.amplitude * envelope * phase.cos()
    }

    /// Compute the energy density of the wave packet.
    ///
    /// Approximated as E = A² * ω² / 2 (kinetic + potential equipartition).
    pub fn energy_density(&self) -> f64 {
        0.5 * self.amplitude * self.amplitude * self.frequency * self.frequency
    }

    /// Compute the instantaneous frequency (identical to carrier ω for linear waves).
    pub fn instantaneous_frequency(&self) -> f64 {
        self.frequency
    }

    /// Compute the wavelength λ = 2π/k.
    pub fn wavelength(&self) -> f64 {
        if self.wavenumber.abs() > 1e-30 {
            std::f64::consts::TAU / self.wavenumber
        } else {
            f64::INFINITY
        }
    }
}

// ---------------------------------------------------------------------------
// DispersionRelation
// ---------------------------------------------------------------------------

/// Dispersion relations for elastic waves in various geometries.
///
/// Provides functions to compute wavenumber k for a given angular frequency ω
/// for bulk, plate, and Love waves.
#[allow(dead_code)]
pub struct DispersionRelation;

impl DispersionRelation {
    /// Compute the wavenumber k for a bulk P-wave at angular frequency ω.
    ///
    /// k = ω / c_P
    ///
    /// # Arguments
    /// * `omega` - Angular frequency (rad/s).
    /// * `params` - Elastic material parameters.
    pub fn bulk_p_wave(omega: f64, params: &ElasticWaveParams) -> f64 {
        omega / params.p_wave_speed()
    }

    /// Compute the wavenumber k for a bulk S-wave at angular frequency ω.
    ///
    /// k = ω / c_S
    ///
    /// # Arguments
    /// * `omega` - Angular frequency (rad/s).
    /// * `params` - Elastic material parameters.
    pub fn bulk_s_wave(omega: f64, params: &ElasticWaveParams) -> f64 {
        omega / params.s_wave_speed()
    }

    /// Compute the wavenumber k for flexural waves in a Kirchhoff plate.
    ///
    /// k = (ρh / D)^(1/4) * √ω  where D = Eh³ / (12(1-ν²)) is the plate stiffness.
    ///
    /// # Arguments
    /// * `omega` - Angular frequency (rad/s).
    /// * `h` - Plate thickness (metres).
    /// * `params` - Elastic material parameters.
    pub fn plate_flexural_wave(omega: f64, h: f64, params: &ElasticWaveParams) -> f64 {
        let nu = params.poissons_ratio();
        let e = params.youngs_modulus();
        let d = e * h * h * h / (12.0 * (1.0 - nu * nu));
        let rho_h = params.density * h;
        (rho_h / d).sqrt().sqrt() * omega.sqrt()
    }

    /// Compute the effective wavenumber for a Love wave in a layer over half-space.
    ///
    /// Love wave dispersion requires the horizontal phase velocity c to satisfy:
    /// tan(k_1 * h * sqrt(c²/cs1² - 1)) = μ2*k_2 / (μ1 * k_1)
    /// where k_i = ω * sqrt(1/csi² - 1/c²).
    ///
    /// This implementation returns an approximate k using the long-wave limit.
    ///
    /// # Arguments
    /// * `omega` - Angular frequency (rad/s).
    /// * `h` - Layer thickness (metres).
    /// * `params1` - Layer material parameters.
    /// * `params2` - Half-space material parameters.
    pub fn love_wave_dispersion(
        omega: f64,
        h: f64,
        params1: &ElasticWaveParams,
        params2: &ElasticWaveParams,
    ) -> f64 {
        let cs1 = params1.s_wave_speed();
        let cs2 = params2.s_wave_speed();
        // Love waves exist only when cs1 < cs2
        if cs1 >= cs2 {
            // No guided Love wave; return S-wave wavenumber in layer
            return omega / cs1;
        }
        // Approximate: c ≈ cs1 + (cs2 - cs1) * (ω * h / π)^(-2) for high frequencies
        // For low frequencies, interpolate between cs1 and cs2
        let cutoff = std::f64::consts::PI * cs1 / h;
        let c = if omega > cutoff {
            cs1 + (cs2 - cs1) / (1.0 + (omega * h / std::f64::consts::PI).powi(2))
        } else {
            cs2
        };
        omega / c
    }

    /// Compute the group velocity dω/dk for a bulk P-wave (non-dispersive).
    ///
    /// # Arguments
    /// * `_omega` - Angular frequency (not used; P-waves are non-dispersive in bulk).
    /// * `params` - Elastic material parameters.
    pub fn bulk_p_group_velocity(_omega: f64, params: &ElasticWaveParams) -> f64 {
        params.p_wave_speed()
    }

    /// Compute the group velocity for a Kirchhoff plate flexural wave.
    ///
    /// v_g = 2 * v_ph (flexural waves have v_g = 2 v_ph in Kirchhoff theory).
    ///
    /// # Arguments
    /// * `omega` - Angular frequency (rad/s).
    /// * `h` - Plate thickness (metres).
    /// * `params` - Elastic material parameters.
    pub fn plate_flexural_group_velocity(omega: f64, h: f64, params: &ElasticWaveParams) -> f64 {
        let k = Self::plate_flexural_wave(omega, h, params);
        if k.abs() > 1e-30 {
            2.0 * omega / k
        } else {
            0.0
        }
    }
}

// ---------------------------------------------------------------------------
// WaveReflection
// ---------------------------------------------------------------------------

/// Zoeppritz-based reflection and transmission coefficients for elastic interfaces.
///
/// Computes P-P, P-S, S-P, and S-S coefficients at a planar interface between
/// two elastic half-spaces.
#[allow(dead_code)]
pub struct WaveReflection;

impl WaveReflection {
    /// Compute the P-P reflection coefficient using simplified Zoeppritz equations.
    ///
    /// Uses the approximation valid for small impedance contrasts:
    /// R_PP ≈ (Z2 - Z1) / (Z2 + Z1)
    ///
    /// # Arguments
    /// * `incident_angle` - Angle of incidence in radians (measured from normal).
    /// * `p1` - Material parameters of the incident half-space.
    /// * `p2` - Material parameters of the transmitted half-space.
    pub fn reflection_coefficient_pp(
        incident_angle: f64,
        p1: &ElasticWaveParams,
        p2: &ElasticWaveParams,
    ) -> f64 {
        let z1 = p1.impedance_p() * incident_angle.cos();
        let snell_sin = p1.p_wave_speed() * incident_angle.sin() / p2.p_wave_speed();
        if snell_sin.abs() > 1.0 {
            // Total internal reflection
            return 1.0;
        }
        let transmitted_angle = snell_sin.asin();
        let z2 = p2.impedance_p() * transmitted_angle.cos();
        (z2 - z1) / (z2 + z1)
    }

    /// Compute the P-P transmission coefficient at an elastic interface.
    ///
    /// T_PP = 1 + R_PP  (energy normalised by impedance ratio)
    ///
    /// # Arguments
    /// * `incident_angle` - Angle of incidence in radians.
    /// * `p1` - Material parameters of the incident half-space.
    /// * `p2` - Material parameters of the transmitted half-space.
    pub fn transmission_coefficient_pp(
        incident_angle: f64,
        p1: &ElasticWaveParams,
        p2: &ElasticWaveParams,
    ) -> f64 {
        let rpp = Self::reflection_coefficient_pp(incident_angle, p1, p2);
        // Amplitude transmission coefficient
        1.0 + rpp
    }

    /// Compute the S-S reflection coefficient at an elastic interface.
    ///
    /// # Arguments
    /// * `incident_angle` - Angle of incidence of S-wave (from normal).
    /// * `p1` - Incident medium parameters.
    /// * `p2` - Transmitted medium parameters.
    pub fn reflection_coefficient_ss(
        incident_angle: f64,
        p1: &ElasticWaveParams,
        p2: &ElasticWaveParams,
    ) -> f64 {
        let z1 = p1.impedance_s() * incident_angle.cos();
        let snell_sin = p1.s_wave_speed() * incident_angle.sin() / p2.s_wave_speed();
        if snell_sin.abs() > 1.0 {
            return 1.0;
        }
        let transmitted_angle = snell_sin.asin();
        let z2 = p2.impedance_s() * transmitted_angle.cos();
        (z2 - z1) / (z2 + z1)
    }

    /// Compute the S-S transmission coefficient.
    ///
    /// # Arguments
    /// * `incident_angle` - Angle of incidence of S-wave.
    /// * `p1` - Incident medium parameters.
    /// * `p2` - Transmitted medium parameters.
    pub fn transmission_coefficient_ss(
        incident_angle: f64,
        p1: &ElasticWaveParams,
        p2: &ElasticWaveParams,
    ) -> f64 {
        1.0 + Self::reflection_coefficient_ss(incident_angle, p1, p2)
    }

    /// Estimate the P-S mode conversion reflection coefficient.
    ///
    /// Approximated using Aki-Richards linearised Zoeppritz for small contrasts.
    ///
    /// # Arguments
    /// * `incident_angle` - Angle of P-wave incidence.
    /// * `p1` - Incident medium.
    /// * `p2` - Transmitted medium.
    pub fn reflection_coefficient_ps(
        incident_angle: f64,
        p1: &ElasticWaveParams,
        p2: &ElasticWaveParams,
    ) -> f64 {
        // Simplified: scale by sin(2θ) * impedance contrast
        let delta_vs = p2.s_wave_speed() - p1.s_wave_speed();
        let mean_vs = 0.5 * (p1.s_wave_speed() + p2.s_wave_speed());

        -2.0 * (p1.s_wave_speed() / p1.p_wave_speed()).powi(2)
            * (delta_vs / mean_vs)
            * (2.0 * incident_angle).sin()
    }

    /// Check for total internal reflection of P-waves.
    ///
    /// # Arguments
    /// * `incident_angle` - Angle of incidence in radians.
    /// * `p1` - Incident medium.
    /// * `p2` - Transmitted medium.
    pub fn is_total_internal_reflection_p(
        incident_angle: f64,
        p1: &ElasticWaveParams,
        p2: &ElasticWaveParams,
    ) -> bool {
        if p1.p_wave_speed() <= p2.p_wave_speed() {
            return false;
        }
        let critical_angle = (p2.p_wave_speed() / p1.p_wave_speed()).asin();
        incident_angle >= critical_angle
    }
}

// ---------------------------------------------------------------------------
// StandingWave
// ---------------------------------------------------------------------------

/// Standing wave solutions for 1D elastic bars and strings.
///
/// Computes mode shapes and natural frequencies for pinned-pinned,
/// clamped-free, and clamped-clamped boundary conditions.
#[allow(dead_code)]
pub struct StandingWave;

impl StandingWave {
    /// Compute the mode shape for a pinned-pinned bar.
    ///
    /// u(x, n) = sin(n π x / L)
    ///
    /// # Arguments
    /// * `x` - Position along the bar (metres).
    /// * `n` - Mode number (1 = fundamental).
    /// * `length` - Bar length (metres).
    pub fn mode_shape(x: f64, n: usize, length: f64) -> f64 {
        ((n as f64 * std::f64::consts::PI * x) / length).sin()
    }

    /// Compute the natural frequency for a pinned-pinned bar.
    ///
    /// f_n = n / (2L) * c
    ///
    /// # Arguments
    /// * `n` - Mode number.
    /// * `length` - Bar length (metres).
    /// * `c` - Wave speed (m/s).
    pub fn natural_frequency(n: usize, length: f64, c: f64) -> f64 {
        n as f64 * c / (2.0 * length)
    }

    /// Compute the mode shape for a clamped-free (cantilever) bar.
    ///
    /// u(x, n) = cos((2n-1)πx / (2L)) for longitudinal modes.
    ///
    /// # Arguments
    /// * `x` - Position (metres).
    /// * `n` - Mode number.
    /// * `length` - Bar length (metres).
    pub fn cantilever_mode_shape(x: f64, n: usize, length: f64) -> f64 {
        let k = (2 * n - 1) as f64 * std::f64::consts::PI / (2.0 * length);
        (k * x).cos()
    }

    /// Compute the natural frequency of a cantilever bar.
    ///
    /// f_n = (2n-1) * c / (4L)
    ///
    /// # Arguments
    /// * `n` - Mode number.
    /// * `length` - Bar length (metres).
    /// * `c` - Wave speed (m/s).
    pub fn cantilever_frequency(n: usize, length: f64, c: f64) -> f64 {
        (2 * n - 1) as f64 * c / (4.0 * length)
    }

    /// Compute time-dependent standing wave displacement.
    ///
    /// u(x, t, n) = sin(nπx/L) * cos(ω_n * t)
    ///
    /// # Arguments
    /// * `x` - Position (metres).
    /// * `t` - Time (seconds).
    /// * `n` - Mode number.
    /// * `length` - Bar length (metres).
    /// * `c` - Wave speed (m/s).
    pub fn evaluate_at_time(x: f64, t: f64, n: usize, length: f64, c: f64) -> f64 {
        let omega = 2.0 * std::f64::consts::PI * Self::natural_frequency(n, length, c);
        Self::mode_shape(x, n, length) * (omega * t).cos()
    }
}

// ---------------------------------------------------------------------------
// FiniteElementWave1D
// ---------------------------------------------------------------------------

/// 1D finite element model for the wave equation on an elastic bar.
///
/// Solves the equation of motion: M ü + K u = f using the central difference
/// explicit time integration scheme.
#[allow(dead_code)]
pub struct FiniteElementWave1D {
    /// Number of elements.
    pub n_elements: usize,
    /// Total bar length (metres).
    pub length: f64,
    /// Elastic material parameters.
    pub params: ElasticWaveParams,
}

impl FiniteElementWave1D {
    /// Create a new 1D FEM wave model.
    ///
    /// # Arguments
    /// * `n_elements` - Number of linear bar elements.
    /// * `length` - Total bar length (metres).
    /// * `params` - Elastic material parameters.
    pub fn new(n_elements: usize, length: f64, params: ElasticWaveParams) -> Self {
        Self {
            n_elements,
            length,
            params,
        }
    }

    /// Number of DOF (nodes = n_elements + 1).
    pub fn n_dof(&self) -> usize {
        self.n_elements + 1
    }

    /// Compute the consistent mass matrix (n_dof × n_dof).
    ///
    /// Uses the consistent mass matrix for linear bar elements:
    /// M_e = (ρ A L_e / 6) * \[2 1; 1 2\]
    ///
    /// The element cross-sectional area A is assumed to be 1 m².
    pub fn mass_matrix(&self) -> Vec<Vec<f64>> {
        let n = self.n_dof();
        let mut m = vec![vec![0.0; n]; n];
        let le = self.length / self.n_elements as f64;
        let rho = self.params.density;
        // Element consistent mass: each element contributes to 2 DOFs
        for e in 0..self.n_elements {
            let scale = rho * le / 6.0;
            m[e][e] += 2.0 * scale;
            m[e][e + 1] += scale;
            m[e + 1][e] += scale;
            m[e + 1][e + 1] += 2.0 * scale;
        }
        m
    }

    /// Compute the stiffness matrix (n_dof × n_dof).
    ///
    /// Uses the bar element stiffness:
    /// K_e = (E A / L_e) * \[1 -1; -1 1\]
    ///
    /// The element cross-sectional area A is assumed to be 1 m².
    pub fn stiffness_matrix(&self) -> Vec<Vec<f64>> {
        let n = self.n_dof();
        let mut k = vec![vec![0.0; n]; n];
        let le = self.length / self.n_elements as f64;
        let e_mod = self.params.youngs_modulus();
        let ke = e_mod / le;
        for e in 0..self.n_elements {
            k[e][e] += ke;
            k[e][e + 1] -= ke;
            k[e + 1][e] -= ke;
            k[e + 1][e + 1] += ke;
        }
        k
    }

    /// Advance the wave solution one step using central-difference integration.
    ///
    /// u_new = (2 u - u_prev + dt² M⁻¹ (f - K u))
    ///
    /// Uses a lumped (diagonal) mass matrix for efficiency.
    ///
    /// # Arguments
    /// * `u` - Current displacement vector (modified in place).
    /// * `u_prev` - Previous displacement vector (modified to hold `u` before update).
    /// * `dt` - Time step (seconds).
    /// * `force` - External force vector.
    pub fn step_central_diff(
        &self,
        u: &mut Vec<f64>,
        u_prev: &mut Vec<f64>,
        dt: f64,
        force: &[f64],
    ) {
        let n = self.n_dof();
        let k = self.stiffness_matrix();
        // Lumped mass: row sum of consistent mass
        let le = self.length / self.n_elements as f64;
        let rho = self.params.density;
        let m_lump: Vec<f64> = (0..n)
            .map(|i| {
                if i == 0 || i == n - 1 {
                    rho * le * 0.5
                } else {
                    rho * le
                }
            })
            .collect();
        // Residual = f - K u
        let mut res = vec![0.0; n];
        for i in 0..n {
            for j in 0..n {
                res[i] -= k[i][j] * u[j];
            }
            if i < force.len() {
                res[i] += force[i];
            }
        }
        let mut u_new = vec![0.0; n];
        for i in 0..n {
            u_new[i] = 2.0 * u[i] - u_prev[i] + dt * dt * res[i] / m_lump[i];
        }
        // Apply fixed BCs: pin both ends
        u_new[0] = 0.0;
        u_new[n - 1] = 0.0;
        *u_prev = u.clone();
        *u = u_new;
    }

    /// Compute the critical time step for stability (CFL condition).
    ///
    /// dt_crit = h / c_P
    pub fn critical_time_step(&self) -> f64 {
        let h = self.length / self.n_elements as f64;
        h / self.params.p_wave_speed()
    }

    /// Compute the natural frequencies of the FEM bar model.
    ///
    /// Uses the exact bar theory: f_n = n * c_P / (2 * L).
    pub fn natural_frequencies(&self, n_modes: usize) -> Vec<f64> {
        let cp = self.params.p_wave_speed();
        (1..=n_modes)
            .map(|n| n as f64 * cp / (2.0 * self.length))
            .collect()
    }
}

// ---------------------------------------------------------------------------
// AttenuationModel
// ---------------------------------------------------------------------------

/// Frequency-dependent wave attenuation model based on quality factor Q.
///
/// Models visco-elastic energy dissipation in soft biological tissues
/// and geological media.
#[allow(dead_code)]
pub struct AttenuationModel {
    /// Quality factor Q (dimensionless). Higher Q → less attenuation.
    pub quality_factor: f64,
    /// Frequency power-law exponent α; attenuation ~ f^α.
    pub frequency_dependence: f64,
    /// Reference wave speed at reference frequency (m/s).
    pub reference_speed: f64,
    /// Reference frequency for the power-law (Hz).
    pub reference_frequency: f64,
}

impl AttenuationModel {
    /// Create a new attenuation model.
    ///
    /// # Arguments
    /// * `quality_factor` - Q factor.
    /// * `frequency_dependence` - Power-law exponent for frequency dependence.
    /// * `reference_speed` - Reference wave speed (m/s).
    /// * `reference_frequency` - Reference frequency (Hz).
    pub fn new(
        quality_factor: f64,
        frequency_dependence: f64,
        reference_speed: f64,
        reference_frequency: f64,
    ) -> Self {
        Self {
            quality_factor,
            frequency_dependence,
            reference_speed,
            reference_frequency,
        }
    }

    /// Compute amplitude decay factor after propagating a given distance.
    ///
    /// A(x) = A0 * exp(-α_att * x)  where  α_att = π f / (Q * c)
    ///
    /// # Arguments
    /// * `distance` - Propagation distance (metres).
    /// * `frequency` - Wave frequency (Hz).
    pub fn amplitude_decay(&self, distance: f64, frequency: f64) -> f64 {
        let alpha_att =
            std::f64::consts::PI * frequency / (self.quality_factor * self.reference_speed);
        (-alpha_att * distance).exp()
    }

    /// Compute frequency-dependent damping ratio ζ(f).
    ///
    /// ζ(f) = 1 / (2Q) * (f / f_ref)^α
    ///
    /// # Arguments
    /// * `f` - Frequency (Hz).
    pub fn frequency_dependent_damping(&self, f: f64) -> f64 {
        let freq_ratio = f / self.reference_frequency;
        (1.0 / (2.0 * self.quality_factor)) * freq_ratio.powf(self.frequency_dependence)
    }

    /// Compute the phase velocity at frequency f accounting for dispersion due to Q.
    ///
    /// Uses Kramers-Kronig approximation: c(f) ≈ c_ref * (1 + ln(f/f_ref) / (π Q))
    ///
    /// # Arguments
    /// * `f` - Frequency (Hz).
    pub fn dispersive_phase_velocity(&self, f: f64) -> f64 {
        let ratio = (f / self.reference_frequency).max(1e-30);
        let correction = ratio.ln() / (std::f64::consts::PI * self.quality_factor);
        self.reference_speed * (1.0 + correction)
    }

    /// Compute the spatial attenuation coefficient α_att at frequency f.
    ///
    /// α_att = π f / (Q c)  (Np/m)
    ///
    /// # Arguments
    /// * `f` - Frequency (Hz).
    pub fn spatial_attenuation_coefficient(&self, f: f64) -> f64 {
        std::f64::consts::PI * f / (self.quality_factor * self.reference_speed)
    }
}

// ---------------------------------------------------------------------------
// PhononicCrystal
// ---------------------------------------------------------------------------

/// Phononic crystal analysis for periodic elastic composites.
///
/// Provides band-gap prediction and effective medium parameter estimation
/// for 1D phononic crystals (alternating material layers).
#[allow(dead_code)]
pub struct PhononicCrystal;

impl PhononicCrystal {
    /// Check whether a given angular frequency falls within a phononic band gap.
    ///
    /// Uses the transfer matrix method for a 1D bilayer unit cell.
    /// A band gap exists when |trace(M)| > 2 for the unit-cell transfer matrix M.
    ///
    /// # Arguments
    /// * `omega` - Angular frequency (rad/s).
    /// * `a` - Lattice constant (unit-cell period in metres).
    /// * `params1` - Material 1 parameters (filling fraction = 0.5).
    /// * `params2` - Material 2 parameters.
    pub fn band_gap_check(
        omega: f64,
        a: f64,
        params1: &ElasticWaveParams,
        params2: &ElasticWaveParams,
    ) -> bool {
        let h1 = a * 0.5;
        let h2 = a * 0.5;
        let c1 = params1.p_wave_speed();
        let c2 = params2.p_wave_speed();
        let z1 = params1.impedance_p();
        let z2 = params2.impedance_p();
        let phi1 = omega * h1 / c1;
        let phi2 = omega * h2 / c2;
        // Transfer matrix trace for bilayer unit cell
        let trace = 2.0 * phi1.cos() * phi2.cos() - (z1 / z2 + z2 / z1) * phi1.sin() * phi2.sin();
        trace.abs() > 2.0
    }

    /// Compute effective medium parameters for a two-phase composite.
    ///
    /// Uses the Voigt–Reuss Hill average for elastic constants.
    ///
    /// # Arguments
    /// * `vf1` - Volume fraction of material 1 (0 ≤ vf1 ≤ 1).
    /// * `p1` - Material 1 parameters.
    /// * `p2` - Material 2 parameters.
    pub fn effective_medium_params(
        vf1: f64,
        p1: &ElasticWaveParams,
        p2: &ElasticWaveParams,
    ) -> ElasticWaveParams {
        let vf2 = 1.0 - vf1;
        // Voigt bound (upper)
        let k1 = p1.bulk_modulus();
        let k2 = p2.bulk_modulus();
        let mu1 = p1.lame_mu;
        let mu2 = p2.lame_mu;
        let k_voigt = vf1 * k1 + vf2 * k2;
        let mu_voigt = vf1 * mu1 + vf2 * mu2;
        // Reuss bound (lower)
        let k_reuss = 1.0 / (vf1 / k1.max(1e-30) + vf2 / k2.max(1e-30));
        let mu_reuss = 1.0 / (vf1 / mu1.max(1e-30) + vf2 / mu2.max(1e-30));
        // Hill average
        let k_eff = 0.5 * (k_voigt + k_reuss);
        let mu_eff = 0.5 * (mu_voigt + mu_reuss);
        let rho_eff = vf1 * p1.density + vf2 * p2.density;
        // Convert K, μ to Lamé λ, μ
        let lambda_eff = k_eff - 2.0 * mu_eff / 3.0;
        ElasticWaveParams::new(rho_eff, lambda_eff, mu_eff)
    }

    /// Compute the Brillouin zone boundary frequency for a periodic bilayer.
    ///
    /// f_BZ = c_eff / (2a) where c_eff is the effective wave speed.
    ///
    /// # Arguments
    /// * `a` - Lattice constant (metres).
    /// * `params1` - Material 1.
    /// * `params2` - Material 2.
    pub fn brillouin_zone_boundary_frequency(
        a: f64,
        params1: &ElasticWaveParams,
        params2: &ElasticWaveParams,
    ) -> f64 {
        let eff = Self::effective_medium_params(0.5, params1, params2);
        eff.p_wave_speed() / (2.0 * a)
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    fn soft_tissue() -> ElasticWaveParams {
        // Approximate soft tissue: ρ=1000, λ=2.2 GPa (bulk), μ=1 kPa
        ElasticWaveParams::new(1000.0, 2.2e9, 1e3)
    }

    fn steel() -> ElasticWaveParams {
        // Steel: ρ=7800, E=200 GPa, ν=0.3 → λ=115 GPa, μ=77 GPa
        ElasticWaveParams::new(7800.0, 1.154e11, 7.69e10)
    }

    fn soft_layer() -> ElasticWaveParams {
        ElasticWaveParams::new(1200.0, 1e6, 5e5)
    }

    fn stiff_halfspace() -> ElasticWaveParams {
        ElasticWaveParams::new(2000.0, 1e9, 5e8)
    }

    // -----------------------------------------------------------------------
    // WaveMode tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_wave_mode_variants() {
        let modes = [
            WaveMode::Longitudinal,
            WaveMode::Transverse,
            WaveMode::Rayleigh,
            WaveMode::Love,
        ];
        assert_eq!(modes.len(), 4);
        assert_eq!(modes[0], WaveMode::Longitudinal);
    }

    #[test]
    fn test_wave_mode_debug() {
        let m = WaveMode::Rayleigh;
        let s = format!("{m:?}");
        assert!(s.contains("Rayleigh"));
    }

    // -----------------------------------------------------------------------
    // ElasticWaveParams tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_p_wave_speed_positive() {
        let p = soft_tissue();
        assert!(p.p_wave_speed() > 0.0);
    }

    #[test]
    fn test_s_wave_speed_positive() {
        let p = steel();
        assert!(p.s_wave_speed() > 0.0);
    }

    #[test]
    fn test_p_wave_faster_than_s() {
        let p = steel();
        assert!(p.p_wave_speed() > p.s_wave_speed());
    }

    #[test]
    fn test_rayleigh_wave_speed_less_than_s() {
        let p = steel();
        let cr = p.rayleigh_wave_speed();
        let cs = p.s_wave_speed();
        assert!(
            cr < cs,
            "Rayleigh wave {cr} should be slower than S-wave {cs}"
        );
    }

    #[test]
    fn test_impedance_p_equals_rho_times_cp() {
        let p = steel();
        let expected = p.density * p.p_wave_speed();
        let got = p.impedance_p();
        assert!((got / expected - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_impedance_s_equals_rho_times_cs() {
        let p = steel();
        let expected = p.density * p.s_wave_speed();
        let got = p.impedance_s();
        assert!((got / expected - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_youngs_modulus_consistency() {
        // For steel: E ≈ 200 GPa
        let p = steel();
        let e = p.youngs_modulus();
        assert!(
            e > 1.5e11 && e < 2.5e11,
            "Steel E should be ~200 GPa, got {e}"
        );
    }

    #[test]
    fn test_poissons_ratio_range() {
        let p = steel();
        let nu = p.poissons_ratio();
        assert!(
            nu > 0.0 && nu < 0.5,
            "Poisson's ratio should be in (0, 0.5): {nu}"
        );
    }

    #[test]
    fn test_bulk_modulus_positive() {
        let p = steel();
        assert!(p.bulk_modulus() > 0.0);
    }

    #[test]
    fn test_steel_p_wave_speed_range() {
        let p = steel();
        let cp = p.p_wave_speed();
        // Steel P-wave ~ 5000-6000 m/s
        assert!(
            cp > 4000.0 && cp < 7000.0,
            "Steel P-wave speed should be ~5000-6000 m/s, got {cp}"
        );
    }

    // -----------------------------------------------------------------------
    // WavePacket tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_wave_packet_creation() {
        let wp = WavePacket::new(1.0, 100.0, 10.0, 10.0, 0.1, 0.0);
        assert!((wp.amplitude - 1.0).abs() < 1e-12);
        assert!((wp.phase_velocity - 10.0).abs() < 1e-10);
    }

    #[test]
    fn test_wave_packet_evaluate_at_center() {
        // At t=0, x=center, envelope=1 → u = A * cos(k * x_0)
        let wp = WavePacket::new(2.0, 10.0, 1.0, 5.0, 0.5, 0.0);
        let u = wp.evaluate_at(0.0, 0.0);
        // x=0, t=0: envelope = exp(0) = 1, cos(0) = 1 → u = 2.0
        assert!(
            (u - 2.0).abs() < 1e-10,
            "u at center should be 2.0, got {u}"
        );
    }

    #[test]
    fn test_wave_packet_energy_positive() {
        let wp = WavePacket::new(1.0, 50.0, 5.0, 10.0, 0.1, 0.0);
        assert!(wp.energy_density() > 0.0);
    }

    #[test]
    fn test_wave_packet_wavelength() {
        let k = 2.0 * PI;
        let wp = WavePacket::new(1.0, 10.0, k, 5.0, 0.5, 0.0);
        let lambda = wp.wavelength();
        assert!(
            (lambda - 1.0).abs() < 1e-10,
            "Wavelength should be 1.0, got {lambda}"
        );
    }

    #[test]
    fn test_wave_packet_zero_wavenumber() {
        let wp = WavePacket::new(1.0, 10.0, 0.0, 5.0, 0.5, 0.0);
        // Phase velocity should be 0 when k=0
        assert_eq!(wp.phase_velocity, 0.0);
        assert!(wp.wavelength().is_infinite());
    }

    // -----------------------------------------------------------------------
    // DispersionRelation tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_dispersion_bulk_p_wave() {
        let p = steel();
        let omega = 1000.0;
        let k = DispersionRelation::bulk_p_wave(omega, &p);
        let expected = omega / p.p_wave_speed();
        assert!((k / expected - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_dispersion_bulk_s_wave() {
        let p = steel();
        let omega = 1000.0;
        let k = DispersionRelation::bulk_s_wave(omega, &p);
        let expected = omega / p.s_wave_speed();
        assert!((k / expected - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_dispersion_bulk_p_faster_than_s() {
        let p = steel();
        let omega = 1000.0;
        let kp = DispersionRelation::bulk_p_wave(omega, &p);
        let ks = DispersionRelation::bulk_s_wave(omega, &p);
        // Smaller k → faster wave (v = ω/k)
        assert!(
            kp < ks,
            "P-wave k ({kp}) should be smaller than S-wave k ({ks})"
        );
    }

    #[test]
    fn test_dispersion_plate_flexural_positive() {
        let p = steel();
        let k = DispersionRelation::plate_flexural_wave(1000.0, 0.01, &p);
        assert!(k > 0.0, "Flexural wavenumber should be positive: {k}");
    }

    #[test]
    fn test_love_wave_dispersion_exists_when_cs1_lt_cs2() {
        let layer = soft_layer();
        let halfspace = stiff_halfspace();
        assert!(layer.s_wave_speed() < halfspace.s_wave_speed());
        let k = DispersionRelation::love_wave_dispersion(1000.0, 0.1, &layer, &halfspace);
        assert!(k > 0.0, "Love wave k should be positive: {k}");
    }

    #[test]
    fn test_love_wave_no_guiding_inverted() {
        // cs1 > cs2 → no Love waves; returns S-wave k in layer
        let layer = stiff_halfspace();
        let halfspace = soft_layer();
        let k = DispersionRelation::love_wave_dispersion(1000.0, 0.1, &layer, &halfspace);
        let expected = 1000.0 / layer.s_wave_speed();
        assert!((k / expected - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_bulk_p_group_velocity() {
        let p = steel();
        let vg = DispersionRelation::bulk_p_group_velocity(1000.0, &p);
        assert!((vg - p.p_wave_speed()).abs() < 1e-10);
    }

    // -----------------------------------------------------------------------
    // WaveReflection tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_reflection_pp_normal_incidence() {
        let p1 = ElasticWaveParams::new(1000.0, 2.0e9, 1e3);
        let p2 = ElasticWaveParams::new(2000.0, 4.0e9, 2e3);
        let r = WaveReflection::reflection_coefficient_pp(0.0, &p1, &p2);
        // Should be positive (harder medium) and in (-1, 1)
        assert!(r > -1.0 && r < 1.0, "R_PP should be in (-1,1): {r}");
    }

    #[test]
    fn test_transmission_plus_reflection_pp() {
        // For normal incidence: T = 1 + R
        let p1 = ElasticWaveParams::new(1000.0, 2.0e9, 1e3);
        let p2 = ElasticWaveParams::new(2000.0, 4.0e9, 2e3);
        let r = WaveReflection::reflection_coefficient_pp(0.0, &p1, &p2);
        let t = WaveReflection::transmission_coefficient_pp(0.0, &p1, &p2);
        assert!(
            (t - (1.0 + r)).abs() < 1e-12,
            "T should equal 1 + R: T={t}, R={r}"
        );
    }

    #[test]
    fn test_reflection_identical_media() {
        let p = steel();
        let r = WaveReflection::reflection_coefficient_pp(0.0, &p, &p);
        assert!(
            r.abs() < 1e-10,
            "Identical media: R_PP should be 0, got {r}"
        );
    }

    #[test]
    fn test_reflection_ss_identical_media() {
        let p = steel();
        let r = WaveReflection::reflection_coefficient_ss(0.0, &p, &p);
        assert!(
            r.abs() < 1e-10,
            "Identical media: R_SS should be 0, got {r}"
        );
    }

    #[test]
    fn test_total_internal_reflection_p() {
        // Put a fast layer over a slow one
        let fast = ElasticWaveParams::new(1000.0, 1e10, 5e9);
        let slow = ElasticWaveParams::new(1000.0, 1e8, 5e7);
        let critical = (slow.p_wave_speed() / fast.p_wave_speed()).asin();
        assert!(
            WaveReflection::is_total_internal_reflection_p(critical + 0.01, &fast, &slow),
            "Should be total internal reflection above critical angle"
        );
        assert!(
            !WaveReflection::is_total_internal_reflection_p(critical - 0.01, &fast, &slow),
            "Should not be total internal reflection below critical angle"
        );
    }

    #[test]
    fn test_reflection_ps_finite() {
        let p1 = soft_tissue();
        let p2 = steel();
        let r = WaveReflection::reflection_coefficient_ps(PI / 6.0, &p1, &p2);
        assert!(r.is_finite(), "R_PS should be finite: {r}");
    }

    // -----------------------------------------------------------------------
    // StandingWave tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_standing_wave_mode_shape_zeros_at_ends() {
        // For a pinned-pinned bar, mode shape = 0 at x=0 and x=L
        let l = 1.0;
        assert!(StandingWave::mode_shape(0.0, 1, l).abs() < 1e-12);
        assert!(StandingWave::mode_shape(l, 1, l).abs() < 1e-12);
    }

    #[test]
    fn test_standing_wave_mode1_peak_at_midspan() {
        let l = 1.0;
        let mid = StandingWave::mode_shape(l / 2.0, 1, l);
        assert!(
            (mid - 1.0).abs() < 1e-10,
            "Mode 1 peak at midspan = 1: {mid}"
        );
    }

    #[test]
    fn test_standing_wave_natural_frequency_scales() {
        let c = 5000.0;
        let l = 1.0;
        let f1 = StandingWave::natural_frequency(1, l, c);
        let f2 = StandingWave::natural_frequency(2, l, c);
        assert!(
            (f2 / f1 - 2.0).abs() < 1e-10,
            "Mode 2 twice mode 1: f1={f1}, f2={f2}"
        );
    }

    #[test]
    fn test_standing_wave_evaluate_at_time() {
        let u = StandingWave::evaluate_at_time(0.5, 0.0, 1, 1.0, 5000.0);
        // t=0: cos(0)=1, so u = mode_shape(0.5, 1, 1.0)
        let expected = StandingWave::mode_shape(0.5, 1, 1.0);
        assert!((u - expected).abs() < 1e-10);
    }

    #[test]
    fn test_cantilever_frequency_increases_with_mode() {
        let c = 5000.0;
        let l = 1.0;
        let f1 = StandingWave::cantilever_frequency(1, l, c);
        let f2 = StandingWave::cantilever_frequency(2, l, c);
        assert!(
            f2 > f1,
            "Higher modes must have higher frequencies: f1={f1}, f2={f2}"
        );
    }

    // -----------------------------------------------------------------------
    // FiniteElementWave1D tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_fem_wave_mass_matrix_size() {
        let fem = FiniteElementWave1D::new(4, 1.0, steel());
        let m = fem.mass_matrix();
        assert_eq!(m.len(), 5);
        assert_eq!(m[0].len(), 5);
    }

    #[test]
    fn test_fem_wave_stiffness_matrix_size() {
        let fem = FiniteElementWave1D::new(4, 1.0, steel());
        let k = fem.stiffness_matrix();
        assert_eq!(k.len(), 5);
        assert_eq!(k[0].len(), 5);
    }

    #[test]
    fn test_fem_wave_mass_matrix_positive_diagonal() {
        let fem = FiniteElementWave1D::new(4, 1.0, steel());
        let m = fem.mass_matrix();
        for i in 0..fem.n_dof() {
            assert!(
                m[i][i] > 0.0,
                "Diagonal of mass matrix should be positive: m[{i}][{i}]={}",
                m[i][i]
            );
        }
    }

    #[test]
    fn test_fem_wave_stiffness_matrix_symmetric() {
        let fem = FiniteElementWave1D::new(4, 1.0, steel());
        let k = fem.stiffness_matrix();
        let n = fem.n_dof();
        for i in 0..n {
            for j in 0..n {
                assert!(
                    (k[i][j] - k[j][i]).abs() < 1e-10,
                    "K should be symmetric at [{i}][{j}]"
                );
            }
        }
    }

    #[test]
    fn test_fem_wave_step_no_crash() {
        let p = ElasticWaveParams::new(1000.0, 1e6, 5e5);
        let fem = FiniteElementWave1D::new(10, 1.0, p);
        let n = fem.n_dof();
        let mut u = vec![0.0; n];
        let mut u_prev = vec![0.0; n];
        let mut force = vec![0.0; n];
        force[n / 2] = 1000.0; // apply force at midpoint
        let dt = fem.critical_time_step() * 0.5;
        for _ in 0..50 {
            fem.step_central_diff(&mut u, &mut u_prev, dt, &force);
        }
    }

    #[test]
    fn test_fem_wave_natural_frequencies_increasing() {
        let fem = FiniteElementWave1D::new(10, 1.0, steel());
        let freqs = fem.natural_frequencies(5);
        for i in 1..freqs.len() {
            assert!(
                freqs[i] > freqs[i - 1],
                "Natural frequencies should increase"
            );
        }
    }

    #[test]
    fn test_fem_wave_critical_dt_positive() {
        let fem = FiniteElementWave1D::new(10, 1.0, steel());
        assert!(fem.critical_time_step() > 0.0);
    }

    // -----------------------------------------------------------------------
    // AttenuationModel tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_attenuation_decay_at_zero_distance() {
        let att = AttenuationModel::new(100.0, 1.0, 1500.0, 1000.0);
        let decay = att.amplitude_decay(0.0, 1000.0);
        assert!(
            (decay - 1.0).abs() < 1e-12,
            "Zero distance → no decay, got {decay}"
        );
    }

    #[test]
    fn test_attenuation_decay_decreases_with_distance() {
        let att = AttenuationModel::new(50.0, 1.0, 1500.0, 1000.0);
        let d1 = att.amplitude_decay(0.1, 1000.0);
        let d2 = att.amplitude_decay(1.0, 1000.0);
        assert!(
            d2 < d1,
            "Amplitude should decay with distance: d1={d1}, d2={d2}"
        );
    }

    #[test]
    fn test_attenuation_high_q_less_damping() {
        let att_low_q = AttenuationModel::new(10.0, 1.0, 1500.0, 1000.0);
        let att_high_q = AttenuationModel::new(1000.0, 1.0, 1500.0, 1000.0);
        let d_low = att_low_q.amplitude_decay(1.0, 1000.0);
        let d_high = att_high_q.amplitude_decay(1.0, 1000.0);
        assert!(
            d_high > d_low,
            "Higher Q → less attenuation: d_low={d_low}, d_high={d_high}"
        );
    }

    #[test]
    fn test_frequency_dependent_damping_at_ref_freq() {
        // At reference frequency, ratio = 1 → damping = 1/(2Q)
        let att = AttenuationModel::new(50.0, 1.0, 1500.0, 1000.0);
        let zeta = att.frequency_dependent_damping(1000.0);
        let expected = 1.0 / (2.0 * 50.0);
        assert!(
            (zeta - expected).abs() < 1e-12,
            "Damping at ref freq: {zeta} vs {expected}"
        );
    }

    #[test]
    fn test_spatial_attenuation_coefficient_positive() {
        let att = AttenuationModel::new(50.0, 1.0, 1500.0, 1000.0);
        let alpha = att.spatial_attenuation_coefficient(1000.0);
        assert!(
            alpha > 0.0,
            "Attenuation coefficient should be positive: {alpha}"
        );
    }

    #[test]
    fn test_dispersive_phase_velocity_at_ref() {
        let att = AttenuationModel::new(50.0, 1.0, 1500.0, 1000.0);
        let c = att.dispersive_phase_velocity(1000.0);
        // At reference: ln(1) = 0, so c = reference_speed
        assert!(
            (c - 1500.0).abs() < 1e-10,
            "At reference freq c should be 1500, got {c}"
        );
    }

    // -----------------------------------------------------------------------
    // PhononicCrystal tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_phononic_band_gap_check_returns_bool() {
        let p1 = ElasticWaveParams::new(1000.0, 2e9, 1e9);
        let p2 = ElasticWaveParams::new(8000.0, 1e11, 5e10);
        // Just check it runs without panic
        let _in_gap = PhononicCrystal::band_gap_check(1e6, 0.001, &p1, &p2);
    }

    #[test]
    fn test_phononic_identical_media_no_gap() {
        let p = steel();
        // Identical media → transfer matrix trace = 2*cos(phi1)*cos(phi2) - 2*sin(phi1)*sin(phi2)
        //                                          = 2*cos(phi1+phi2) ≤ 2
        let in_gap = PhononicCrystal::band_gap_check(1e6, 0.001, &p, &p);
        assert!(!in_gap, "Identical media should not have a band gap");
    }

    #[test]
    fn test_phononic_effective_medium_vf0() {
        // vf1 = 0 → all material 2
        let p1 = steel();
        let p2 = soft_tissue();
        let eff = PhononicCrystal::effective_medium_params(0.0, &p1, &p2);
        // Effective density should equal p2 density
        assert!(
            (eff.density - p2.density).abs() < 1e-6,
            "vf1=0 → density should be p2: {}",
            eff.density
        );
    }

    #[test]
    fn test_phononic_effective_medium_vf1() {
        // vf1 = 1 → all material 1
        let p1 = steel();
        let p2 = soft_tissue();
        let eff = PhononicCrystal::effective_medium_params(1.0, &p1, &p2);
        assert!(
            (eff.density - p1.density).abs() < 1e-6,
            "vf1=1 → density should be p1: {}",
            eff.density
        );
    }

    #[test]
    fn test_phononic_effective_medium_density_average() {
        let p1 = ElasticWaveParams::new(1000.0, 1e9, 5e8);
        let p2 = ElasticWaveParams::new(2000.0, 2e9, 1e9);
        let eff = PhononicCrystal::effective_medium_params(0.5, &p1, &p2);
        let expected_rho = 1500.0;
        assert!(
            (eff.density - expected_rho).abs() < 1e-6,
            "Average density: {}",
            eff.density
        );
    }

    #[test]
    fn test_phononic_brillouin_zone_frequency_positive() {
        let p1 = steel();
        let p2 = soft_tissue();
        let f = PhononicCrystal::brillouin_zone_boundary_frequency(0.001, &p1, &p2);
        assert!(f > 0.0, "BZ boundary frequency should be positive: {f}");
    }

    // -----------------------------------------------------------------------
    // Integration tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_wave_round_trip_consistency() {
        // R + T amplitude should be consistent with impedance
        let p1 = ElasticWaveParams::new(1200.0, 1e9, 5e8);
        let p2 = ElasticWaveParams::new(2400.0, 2e9, 1e9);
        let r = WaveReflection::reflection_coefficient_pp(0.0, &p1, &p2);
        let t = WaveReflection::transmission_coefficient_pp(0.0, &p1, &p2);
        assert!(r.is_finite() && t.is_finite());
        // Power reflection + transmission ~ 1
        let z1 = p1.impedance_p();
        let z2 = p2.impedance_p();
        let pr = r * r;
        let pt = t * t * z1 / z2;
        assert!(
            (pr + pt - 1.0).abs() < 1e-10,
            "Energy conservation: R²+T²Z1/Z2={}",
            pr + pt
        );
    }

    #[test]
    fn test_fem_wave_boundary_conditions_enforced() {
        let p = ElasticWaveParams::new(1000.0, 1e6, 5e5);
        let fem = FiniteElementWave1D::new(5, 1.0, p);
        let n = fem.n_dof();
        let mut u = vec![0.0; n];
        let mut u_prev = vec![0.0; n];
        let mut force = vec![0.0; n];
        force[n / 2] = 1.0;
        let dt = fem.critical_time_step() * 0.4;
        for _ in 0..20 {
            fem.step_central_diff(&mut u, &mut u_prev, dt, &force);
        }
        // Both ends must be pinned (zero displacement)
        assert!(u[0].abs() < 1e-12, "Left end must be zero: {}", u[0]);
        assert!(
            u[n - 1].abs() < 1e-12,
            "Right end must be zero: {}",
            u[n - 1]
        );
    }
}
