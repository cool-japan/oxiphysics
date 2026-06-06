//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use std::f64::consts::PI;

use super::types_impl::SmaPhase;

/// Shape memory alloy model using the Tanaka-Liang-Rogers cosine model.
/// Tracks martensite fraction ξ ∈ \[0, 1\], where 1 = full martensite.
#[derive(Debug, Clone)]
pub struct ShapeMemoryAlloy {
    /// Martensite start temperature \[K\].
    pub ms: f64,
    /// Martensite finish temperature \[K\].
    pub mf: f64,
    /// Austenite start temperature \[K\].
    pub a_s: f64,
    /// Austenite finish temperature \[K\].
    pub af: f64,
    /// Current martensite fraction ξ ∈ \[0, 1\].
    pub xi: f64,
    /// Austenite elastic modulus \[Pa\].
    pub e_a: f64,
    /// Martensite elastic modulus \[Pa\].
    pub e_m: f64,
    /// Maximum recoverable strain ε_L.
    pub max_strain: f64,
    /// Stress influence coefficient for martensite transformation \[Pa/K\].
    pub cm: f64,
    /// Stress influence coefficient for austenite transformation \[Pa/K\].
    pub ca: f64,
}
impl ShapeMemoryAlloy {
    /// Create an SMA model with transformation temperatures and elastic moduli.
    pub fn new(
        ms: f64,
        mf: f64,
        a_s: f64,
        af: f64,
        e_a: f64,
        e_m: f64,
        max_strain: f64,
        cm: f64,
        ca: f64,
    ) -> Self {
        Self {
            ms,
            mf,
            a_s,
            af,
            xi: 1.0,
            e_a,
            e_m,
            max_strain,
            cm,
            ca,
        }
    }
    /// Create a standard NiTi Nitinol model with typical parameters.
    pub fn nitinol() -> Self {
        Self::new(291.0, 273.0, 307.0, 325.0, 75e9, 28e9, 0.08, 8e6, 13e6)
    }
    /// Elastic modulus at current martensite fraction.
    pub fn elastic_modulus(&self) -> f64 {
        self.e_a + self.xi * (self.e_m - self.e_a)
    }
    /// Update martensite fraction for given temperature T \[K\] and stress σ \[Pa\].
    /// Returns updated ξ.
    pub fn update_phase(&mut self, temperature: f64, stress: f64) -> f64 {
        let ms_stress = self.ms + stress / self.cm;
        let mf_stress = self.mf + stress / self.cm;
        if temperature <= ms_stress && temperature >= mf_stress {
            let xi_m =
                (1.0 + (PI * (temperature - ms_stress) / (mf_stress - ms_stress)).cos()) / 2.0;
            self.xi = xi_m.clamp(self.xi, 1.0);
        }
        let as_stress = self.a_s + stress / self.ca;
        let af_stress = self.af + stress / self.ca;
        if temperature >= as_stress && temperature <= af_stress {
            let xi_a = self.xi
                * (1.0 + (PI * (temperature - as_stress) / (af_stress - as_stress)).cos())
                / 2.0;
            self.xi = xi_a.clamp(0.0, self.xi);
        }
        if temperature < mf_stress {
            self.xi = 1.0;
        }
        if temperature > af_stress {
            self.xi = 0.0;
        }
        self.xi
    }
    /// Current phase classification.
    pub fn phase(&self) -> SmaPhase {
        if self.xi > 0.99 {
            SmaPhase::Martensite
        } else if self.xi < 0.01 {
            SmaPhase::Austenite
        } else {
            SmaPhase::Mixed
        }
    }
    /// Recovery stress for given strain \[Pa\].
    pub fn recovery_stress(&self, strain: f64) -> f64 {
        self.elastic_modulus() * (strain - self.xi * self.max_strain)
    }
}
/// Type of electroactive polymer.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum EapType {
    /// Ionic EAP: actuation by ion migration (IPMC, conductive polymer).
    Ionic,
    /// Electronic EAP: actuation by electrostatic forces (dielectric elastomer, PVDF).
    Electronic,
}
/// Magnetic shape memory effect in Heusler alloys (Ni-Mn-Ga).
///
/// Models variant reorientation driven by magnetic field via energy minimization.
/// Strain output up to ~12 % for Ni2MnGa.
#[derive(Debug, Clone, Copy)]
pub struct MagneticShape {
    /// Maximum transformation strain ε_max.
    pub max_strain: f64,
    /// Magnetic anisotropy energy density K_u \[J/m³\].
    pub k_u: f64,
    /// Saturation magnetization M_s \[A/m\].
    pub m_sat: f64,
    /// Critical field for reorientation H_cr \[A/m\].
    pub h_critical: f64,
    /// Current variant fraction η ∈ \[0, 1\].
    pub eta: f64,
    /// Elastic modulus \[Pa\].
    pub elastic_modulus: f64,
}
impl MagneticShape {
    /// Create a magnetic shape memory material.
    pub fn new(
        max_strain: f64,
        k_u: f64,
        m_sat: f64,
        h_critical: f64,
        elastic_modulus: f64,
    ) -> Self {
        Self {
            max_strain,
            k_u,
            m_sat,
            h_critical,
            eta: 0.0,
            elastic_modulus,
        }
    }
    /// Standard Ni2MnGa model.
    pub fn ni2mnga() -> Self {
        Self::new(0.06, 1.65e5, 600e3, 300e3, 2.0e9)
    }
    /// Zeeman energy density at field H (field along easy axis) \[J/m³\].
    pub fn zeeman_energy(&self, h_field: f64, eta: f64) -> f64 {
        const MU0: f64 = 4.0 * PI * 1e-7;
        -MU0 * self.m_sat * h_field * eta
    }
    /// Anisotropy energy density \[J/m³\].
    pub fn anisotropy_energy(&self, eta: f64) -> f64 {
        self.k_u * eta * (1.0 - eta)
    }
    /// Update variant fraction η from applied field H \[A/m\].
    pub fn update(&mut self, h_field: f64) {
        if h_field.abs() > self.h_critical {
            self.eta = if h_field > 0.0 { 1.0 } else { 0.0 };
        } else {
            let fraction = (h_field / self.h_critical).clamp(-1.0, 1.0);
            self.eta = (0.5 + 0.5 * fraction).clamp(0.0, 1.0);
        }
    }
    /// Macroscopic strain at current η.
    pub fn strain(&self) -> f64 {
        self.eta * self.max_strain
    }
    /// Blocking stress (stress required to prevent strain) \[Pa\].
    pub fn blocking_stress(&self) -> f64 {
        self.elastic_modulus * self.strain()
    }
}
