//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use std::f64::consts::PI;

/// Frequency weighting filters for noise measurement.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum WeightingFilter {
    /// A-weighting (most common, approximates human hearing sensitivity)
    A,
    /// B-weighting (intermediate, less used)
    B,
    /// C-weighting (flat response above ~30 Hz, used for peak measurements)
    C,
}
/// Porous absorber model: rigid-frame porous material (Delany-Bazley model).
///
/// Empirical complex impedance and propagation constant for fibrous materials
/// as a function of resistivity R_f (Pa·s/m²) and frequency f (Hz).
///
/// Reference: Delany & Bazley, Appl. Acoust. 3, 105 (1970).
#[derive(Debug, Clone)]
pub struct PorousAbsorber {
    /// Flow resistivity (Pa·s/m²), typically 5000–50000 for common absorbers.
    pub flow_resistivity: f64,
    /// Thickness of the absorber layer (m).
    pub thickness: f64,
}
impl PorousAbsorber {
    /// Create a new porous absorber.
    pub fn new(flow_resistivity: f64, thickness: f64) -> Self {
        Self {
            flow_resistivity,
            thickness,
        }
    }
    /// Dimensionless frequency parameter X = ρ_air * f / R_f.
    fn x_param(&self, freq: f64) -> f64 {
        1.204 * freq / self.flow_resistivity
    }
    /// Real part of normalised characteristic impedance (Delany-Bazley).
    pub fn impedance_real(&self, freq: f64) -> f64 {
        let x = self.x_param(freq);
        1.0 + 9.08 * x.powf(-0.75)
    }
    /// Imaginary part of normalised characteristic impedance (Delany-Bazley).
    pub fn impedance_imag(&self, freq: f64) -> f64 {
        let x = self.x_param(freq);
        -11.9 * x.powf(-0.73)
    }
    /// Real part of propagation constant k_r = ω/c * α_r.
    pub fn propagation_real(&self, freq: f64) -> f64 {
        let x = self.x_param(freq);
        let omega = 2.0 * PI * freq;
        let c_air = 343.0;
        (omega / c_air) * (1.0 + 10.8 * x.powf(-0.70))
    }
    /// Imaginary (attenuation) part of propagation constant.
    pub fn propagation_imag(&self, freq: f64) -> f64 {
        let x = self.x_param(freq);
        let omega = 2.0 * PI * freq;
        let c_air = 343.0;
        (omega / c_air) * (10.3 * x.powf(-0.59))
    }
    /// Normal incidence absorption coefficient (approximate).
    ///
    /// Uses the impedance at the face of the absorber layer (rigid backing).
    /// α = 1 − |R|²
    pub fn absorption_coefficient(&self, freq: f64) -> f64 {
        let z_r = self.impedance_real(freq);
        let z_i = self.impedance_imag(freq);
        let ki = self.propagation_imag(freq);
        let kd = ki * self.thickness;
        let coth_kd = if kd > 50.0 { 1.0 } else { 1.0 / kd.tanh() };
        let zs_r = z_r * coth_kd;
        let zs_i = z_i * coth_kd;
        let num_re = zs_r - 1.0;
        let denom_re = zs_r + 1.0;
        let r_sq = (num_re * num_re + zs_i * zs_i) / (denom_re * denom_re + zs_i * zs_i);
        (1.0 - r_sq).clamp(0.0, 1.0)
    }
}
/// Locally resonant acoustic metamaterial (LRAM) unit cell model.
///
/// A heavy mass m coated with a soft layer in a matrix medium.
/// Near the local resonance frequency f_0 = (1/2π) * sqrt(k/m),
/// the effective mass density becomes negative.
///
/// Reference: Liu et al., Science 289, 1734 (2000).
#[derive(Debug, Clone)]
pub struct LocallyResonantUnit {
    /// Mass of the heavy inclusion (kg).
    pub mass: f64,
    /// Stiffness of the soft coating (N/m).
    pub stiffness: f64,
    /// Volume of the unit cell (m³).
    pub cell_volume: f64,
    /// Background medium density (kg/m³).
    pub rho_matrix: f64,
}
impl LocallyResonantUnit {
    /// Create a new locally resonant unit cell.
    pub fn new(mass: f64, stiffness: f64, cell_volume: f64, rho_matrix: f64) -> Self {
        Self {
            mass,
            stiffness,
            cell_volume,
            rho_matrix,
        }
    }
    /// Local resonance frequency (Hz).
    pub fn resonance_frequency(&self) -> f64 {
        (1.0 / (2.0 * PI)) * (self.stiffness / self.mass).sqrt()
    }
    /// Effective mass density at frequency f (Hz).
    ///
    /// rho_eff = rho_matrix + (mass / V_cell) * (omega_0² / (omega_0² - omega²))
    pub fn effective_density(&self, freq: f64) -> f64 {
        let omega_0 = 2.0 * PI * self.resonance_frequency();
        let omega = 2.0 * PI * freq;
        let omega_0_sq = omega_0 * omega_0;
        let omega_sq = omega * omega;
        let denom = omega_0_sq - omega_sq;
        if denom.abs() < 1e-30 {
            return f64::INFINITY;
        }
        let rho_incl = self.mass / self.cell_volume;
        self.rho_matrix + rho_incl * omega_0_sq / denom
    }
    /// Whether the effective density is negative at frequency f.
    pub fn is_negative_mass(&self, freq: f64) -> bool {
        let rho_eff = self.effective_density(freq);
        rho_eff < 0.0
    }
}
/// Frequency-dependent absorption model using multiple relaxation mechanisms.
///
/// Models absorption in water including classical viscous losses and
/// a magnesium sulfate (MgSO₄) relaxation term.
///
/// Reference: Francois & Garrison (1982), simplified.
#[derive(Debug, Clone)]
pub struct WaterAbsorption {
    /// Temperature \[°C\]
    pub temperature_c: f64,
    /// Salinity \[ppt (parts per thousand)\]
    pub salinity_ppt: f64,
    /// Depth \[m\]
    pub depth_m: f64,
}
impl WaterAbsorption {
    /// Create a new WaterAbsorption model.
    pub fn new(temperature_c: f64, salinity_ppt: f64, depth_m: f64) -> Self {
        Self {
            temperature_c,
            salinity_ppt,
            depth_m,
        }
    }
    /// Freshwater (zero salinity) at 20°C, sea surface.
    pub fn freshwater() -> Self {
        Self::new(20.0, 0.0, 0.0)
    }
    /// Seawater at typical ocean conditions (15°C, 35 ppt, 0 m).
    pub fn seawater() -> Self {
        Self::new(15.0, 35.0, 0.0)
    }
    /// Absorption coefficient \[dB/km\] at frequency f \[kHz\].
    ///
    /// Uses the Francois-Garrison formula (simplified).
    pub fn absorption_db_per_km(&self, freq_khz: f64) -> f64 {
        let t = self.temperature_c;
        let s = self.salinity_ppt;
        let d = self.depth_m / 1000.0;
        let f = freq_khz;
        let f2 = f * f;
        let f1 = 0.78 * (s / 35.0).sqrt() * (t / 26.0).exp();
        let a1 = 0.106 * f1 * f2 / (f1 * f1 + f2) * (-d / 6.0).exp();
        let f2_rel = 42.0 * (t / 17.0).exp();
        let a2 = 0.52 * (s / 35.0) * (1.0 + t / 40.0) * f2_rel * f2 / (f2_rel * f2_rel + f2)
            * (-d / 6.0).exp();
        let a3 = 4.9e-4 * f2 * (-t / 27.0 + d / 17.0).exp();
        a1 + a2 + a3
    }
}
/// Material properties relevant to acoustic wave propagation.
#[derive(Debug, Clone)]
pub struct AcousticMaterial {
    /// Mass density in kg/m³
    pub density: f64,
    /// Bulk modulus in Pa
    pub bulk_modulus: f64,
    /// Shear modulus in Pa (0 for fluids)
    pub shear_modulus: f64,
    /// Structural loss factor (dimensionless)
    pub loss_factor: f64,
}
impl AcousticMaterial {
    /// Create a new acoustic material.
    pub fn new(density: f64, bulk_modulus: f64, shear_modulus: f64, loss_factor: f64) -> Self {
        Self {
            density,
            bulk_modulus,
            shear_modulus,
            loss_factor,
        }
    }
    /// Water at 20°C.
    pub fn water() -> Self {
        Self {
            density: 998.2,
            bulk_modulus: 2.18e9,
            shear_modulus: 0.0,
            loss_factor: 0.0,
        }
    }
    /// Steel (structural).
    pub fn steel() -> Self {
        Self {
            density: 7850.0,
            bulk_modulus: 166.67e9,
            shear_modulus: 80.0e9,
            loss_factor: 0.001,
        }
    }
    /// Air at 20°C, 1 atm.
    pub fn air() -> Self {
        Self {
            density: 1.204,
            bulk_modulus: 141.9e3,
            shear_modulus: 0.0,
            loss_factor: 0.0,
        }
    }
    /// Concrete (normal weight).
    pub fn concrete() -> Self {
        Self {
            density: 2300.0,
            bulk_modulus: 13.33e9,
            shear_modulus: 10.0e9,
            loss_factor: 0.02,
        }
    }
    /// Longitudinal (compressional) wave velocity: c_L = sqrt((K + 4G/3) / rho).
    pub fn longitudinal_velocity(&self) -> f64 {
        let m_modulus = self.bulk_modulus + 4.0 * self.shear_modulus / 3.0;
        (m_modulus / self.density).sqrt()
    }
    /// Shear wave velocity: c_S = sqrt(G / rho). Returns 0 for fluids.
    pub fn shear_velocity(&self) -> f64 {
        if self.shear_modulus == 0.0 {
            0.0
        } else {
            (self.shear_modulus / self.density).sqrt()
        }
    }
    /// Characteristic acoustic impedance: Z = rho * c_L (Pa·s/m).
    pub fn acoustic_impedance(&self) -> f64 {
        self.density * self.longitudinal_velocity()
    }
    /// Compute the insertion loss (IL) of a barrier made of this material.
    ///
    /// The insertion loss is the reduction in sound power level (SPL) achieved
    /// by inserting a barrier between a source and receiver. Using the
    /// simplified mass-law / Maekawa formula for a thin barrier:
    ///
    /// IL = 20 · log10(1 + π · N · f · (m_s / (ρ_air · c_air)) )
    ///
    /// where N = Fresnel number characterising diffraction over the barrier,
    /// and m_s = ρ · h is the surface mass density of the barrier panel.
    ///
    /// For thin panels the dominant mechanism is mass-law transmission loss,
    /// so we combine barrier TL with Maekawa diffraction:
    ///
    /// IL ≈ TL_mass + ΔL_diffraction
    ///
    /// Here we use the practical formula:
    ///
    /// IL = 10 · log10( 1 + 20 · N )   \[dB\]  (Maekawa 1968, simplified)
    ///
    /// # Arguments
    /// * `frequency`     - Frequency f \[Hz\]
    /// * `thickness`     - Barrier panel thickness h \[m\]
    /// * `fresnel_number` - Fresnel number N = 2·δ/λ where δ is path-length
    ///   difference (dimensionless, N ≥ 0)
    ///
    /// # Returns
    /// Insertion loss IL \[dB\]
    pub fn compute_insertion_loss(
        &self,
        frequency: f64,
        thickness: f64,
        fresnel_number: f64,
    ) -> f64 {
        let m_s = self.density * thickness;
        let rho_air = 1.204_f64;
        let c_air = 343.0_f64;
        let tl_mass = 20.0 * (1.0 + PI * frequency * m_s / (rho_air * c_air)).log10();
        let il_diffraction = if fresnel_number > 0.0 {
            10.0 * (1.0 + 20.0 * fresnel_number).log10()
        } else {
            0.0
        };
        tl_mass + il_diffraction
    }
    /// Compute the Noise Reduction Coefficient (NRC) for this material.
    ///
    /// The NRC is a single-number rating of sound absorption, defined as the
    /// arithmetic mean of the Sabine absorption coefficients at 250, 500,
    /// 1000, and 2000 Hz, rounded to the nearest 0.05.
    ///
    /// For a material characterised by its loss factor η, the statistical
    /// absorption coefficient at each frequency band is estimated from the
    /// relationship between loss factor and random-incidence absorption
    /// coefficient:
    ///
    /// α_stat(f) ≈ η(f) · m_s · ω / (ρ_air · c_air)
    ///
    /// where m_s = ρ · h is the surface mass density and ω = 2πf. When η
    /// is independent of frequency (as stored in `self.loss_factor`), the
    /// coefficients scale with f, so we compute and average them.
    ///
    /// The result is clamped to \[0, 1\].
    ///
    /// # Arguments
    /// * `thickness` - Panel thickness h \[m\]
    ///
    /// # Returns
    /// NRC \[dimensionless, 0–1\]
    pub fn compute_noise_reduction_coefficient(&self, thickness: f64) -> f64 {
        let bands = [250.0_f64, 500.0, 1000.0, 2000.0];
        let rho_air = 1.204_f64;
        let c_air = 343.0_f64;
        let m_s = self.density * thickness;
        let sum: f64 = bands
            .iter()
            .map(|&f| {
                let omega = 2.0 * PI * f;
                let alpha = self.loss_factor * m_s * omega / (rho_air * c_air);
                alpha.clamp(0.0, 1.0)
            })
            .sum();
        let nrc_raw = sum / bands.len() as f64;
        let nrc_rounded = (nrc_raw / 0.05).round() * 0.05;
        nrc_rounded.clamp(0.0, 1.0)
    }
    /// Compute the transmission loss (TL) by the mass law.
    ///
    /// The mass law is the classical result for a limp, homogeneous panel:
    ///
    /// TL = 20 · log10(m_s · f) − 47   \[dB\]  (for normal incidence, SI units)
    ///
    /// or equivalently for random-incidence (field incidence):
    ///
    /// TL = 20 · log10(m_s · f) − 47.3  \[dB\]
    ///
    /// where m_s = ρ · h \[kg/m²\] is the surface mass density.
    /// The constant 47 arises from 20·log10(2π / (ρ_air · c_air)).
    ///
    /// At double the frequency or double the mass, TL increases by ~6 dB
    /// ("mass law doubling rule").
    ///
    /// # Arguments
    /// * `frequency` - Frequency f \[Hz\]
    /// * `thickness` - Panel thickness h \[m\]
    ///
    /// # Returns
    /// Transmission loss TL \[dB\]
    pub fn compute_transmission_loss_mass_law(&self, frequency: f64, thickness: f64) -> f64 {
        let m_s = self.density * thickness;
        let rho_air = 1.204_f64;
        let c_air = 343.0_f64;
        let arg = PI * m_s * frequency / (rho_air * c_air);
        20.0 * arg.log10()
    }
}
/// Rubber (natural rubber, acoustic grade).
impl AcousticMaterial {
    /// Natural rubber (acoustic grade).
    pub fn rubber() -> Self {
        Self {
            density: 1200.0,
            bulk_modulus: 2.0e9,
            shear_modulus: 0.5e6,
            loss_factor: 0.1,
        }
    }
    /// Aluminium alloy (2024-T3).
    pub fn aluminium() -> Self {
        Self {
            density: 2780.0,
            bulk_modulus: 73.1e9,
            shear_modulus: 27.0e9,
            loss_factor: 0.0002,
        }
    }
    /// Glass (borosilicate).
    pub fn glass() -> Self {
        Self {
            density: 2230.0,
            bulk_modulus: 46.0e9,
            shear_modulus: 26.0e9,
            loss_factor: 0.0005,
        }
    }
    /// Sound speed through the material \[m/s\] (alias for longitudinal_velocity).
    pub fn sound_speed(&self) -> f64 {
        self.longitudinal_velocity()
    }
    /// Specific acoustic impedance Z = ρ·c normalised to air impedance.
    ///
    /// Z_air = 413 Pa·s/m (20°C, 1 atm).
    pub fn relative_impedance(&self) -> f64 {
        self.acoustic_impedance() / 413.0
    }
}
/// Ultrasonic non-destructive evaluation (NDE) signal model.
///
/// Simulates pulse-echo time-of-flight and amplitude for a pulse reflected
/// from a planar reflector.
#[derive(Debug, Clone)]
pub struct UltrasonicNDE {
    /// Wave velocity in the test material \[m/s\]
    pub velocity: f64,
    /// Ultrasonic frequency \[Hz\]
    pub frequency: f64,
    /// Initial pulse amplitude \[Pa\] (peak pressure)
    pub pulse_amplitude: f64,
    /// Attenuation coefficient \[dB/m\]
    pub attenuation_db_m: f64,
}
impl UltrasonicNDE {
    /// Create a new `UltrasonicNDE` model.
    pub fn new(velocity: f64, frequency: f64, pulse_amplitude: f64, attenuation_db_m: f64) -> Self {
        Self {
            velocity,
            frequency,
            pulse_amplitude,
            attenuation_db_m,
        }
    }
    /// Round-trip time-of-flight \[s\] for a reflector at depth `depth_m` \[m\].
    pub fn time_of_flight(&self, depth_m: f64) -> f64 {
        2.0 * depth_m / self.velocity
    }
    /// Echo amplitude \[Pa\] after round-trip path of `2·depth` \[m\].
    ///
    /// Accounts for beam spreading (1/r law) and material attenuation.
    pub fn echo_amplitude(&self, depth_m: f64, reflection_coeff: f64) -> f64 {
        if depth_m <= 0.0 {
            return self.pulse_amplitude * reflection_coeff.abs();
        }
        let d = 2.0 * depth_m;
        let geo = 1.0 / (1.0 + d);
        let alpha_np = self.attenuation_db_m / 8.686;
        let att = (-alpha_np * d).exp();
        self.pulse_amplitude * reflection_coeff.abs() * geo * att
    }
    /// Wavelength in the material \[m\].
    pub fn wavelength(&self) -> f64 {
        self.velocity / self.frequency
    }
    /// Near-field (Fresnel) distance \[m\] for a circular transducer of radius `r` \[m\].
    ///
    /// `N = r² / λ`
    pub fn near_field_distance(&self, transducer_radius: f64) -> f64 {
        transducer_radius * transducer_radius / self.wavelength()
    }
    /// Lateral resolution at depth `d` \[m\] (Rayleigh criterion for focused beam):
    ///
    /// `x_r ≈ 1.22 · λ · d / (2·r)`
    pub fn lateral_resolution(&self, depth_m: f64, transducer_radius: f64) -> f64 {
        1.22 * self.wavelength() * depth_m / (2.0 * transducer_radius)
    }
}
/// Simplified 1D phononic crystal band gap estimator.
///
/// For a 1D bi-material phononic crystal with alternating layers A and B:
///
/// Band gap frequency: f_gap = c_eff / (2 * a)  where a is the lattice constant.
/// The gap width depends on the impedance mismatch.
///
/// This uses a simplified formula from Sigalas & Economou (1992):
/// Δω/ω_0 ≈ (2/π) * |arcsin((Z_b - Z_a) / (Z_b + Z_a))|
#[derive(Debug, Clone)]
pub struct PhononicCrystal1D {
    /// Material A impedance (Pa·s/m).
    pub z_a: f64,
    /// Material B impedance (Pa·s/m).
    pub z_b: f64,
    /// Wave speed in material A (m/s).
    pub c_a: f64,
    /// Wave speed in material B (m/s).
    pub c_b: f64,
    /// Thickness of material A layer (m).
    pub d_a: f64,
    /// Thickness of material B layer (m).
    pub d_b: f64,
}
impl PhononicCrystal1D {
    /// Create a new 1D phononic crystal.
    pub fn new(z_a: f64, z_b: f64, c_a: f64, c_b: f64, d_a: f64, d_b: f64) -> Self {
        Self {
            z_a,
            z_b,
            c_a,
            c_b,
            d_a,
            d_b,
        }
    }
    /// Lattice constant a = d_a + d_b.
    pub fn lattice_constant(&self) -> f64 {
        self.d_a + self.d_b
    }
    /// First band gap centre frequency (Hz).
    ///
    /// Based on the Bragg condition: 2a = λ → f_Bragg = c_eff / (2a).
    pub fn bragg_frequency(&self) -> f64 {
        let c_eff = (self.c_a * self.d_a + self.c_b * self.d_b) / self.lattice_constant();
        c_eff / (2.0 * self.lattice_constant())
    }
    /// Estimated relative band gap width Δf/f_0.
    pub fn relative_gap_width(&self) -> f64 {
        let r = (self.z_b - self.z_a).abs() / (self.z_a + self.z_b);
        (2.0 / PI) * r.asin()
    }
    /// Absolute band gap width (Hz).
    pub fn gap_width_hz(&self) -> f64 {
        self.relative_gap_width() * self.bragg_frequency()
    }
    /// Whether a given frequency lies in the first band gap.
    pub fn is_in_gap(&self, freq: f64) -> bool {
        let f0 = self.bragg_frequency();
        let half_gap = 0.5 * self.gap_width_hz();
        (freq - f0).abs() < half_gap
    }
}
