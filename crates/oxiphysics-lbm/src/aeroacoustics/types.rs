//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::fwh_acoustic_pressure;
#[allow(unused_imports)]
use super::functions::*;
use super::functions_2::spl_from_pressure_fluctuation;
use std::f64::consts::PI;

/// A-weighting filter coefficients structure.
///
/// Stores the numerator and denominator coefficients of the A-weighting
/// transfer function in continuous-time bilinear form.
#[allow(dead_code)]
#[derive(Clone, Debug)]
pub struct AWeightingFilter {
    /// Poles of the A-weighting filter (rad/s)
    pub poles: [f64; 4],
    /// Gain normalization constant
    pub gain: f64,
}
impl AWeightingFilter {
    /// Create the standard A-weighting filter.
    ///
    /// Based on IEC 61672-1 / ANSI S1.42 specification.
    #[allow(dead_code)]
    pub fn new() -> Self {
        Self {
            poles: [
                2.0 * PI * 20.6,
                2.0 * PI * 107.7,
                2.0 * PI * 737.9,
                2.0 * PI * 12194.0,
            ],
            gain: 7.39705e9,
        }
    }
    /// Compute A-weighting magnitude (linear, not dB) at frequency f (Hz).
    #[allow(dead_code)]
    pub fn magnitude(&self, f: f64) -> f64 {
        if f < 1e-6 {
            return 0.0;
        }
        let w = 2.0 * PI * f;
        let w2 = w * w;
        let p1 = self.poles[0];
        let p2 = self.poles[1];
        let p3 = self.poles[2];
        let p4 = self.poles[3];
        let num = self.gain * w2 * w2;
        let den = (w2 + p1 * p1) * ((w2 + p2 * p2) * (w2 + p3 * p3)).sqrt() * (w2 + p4 * p4);
        if den < 1e-30 {
            return 0.0;
        }
        num / den
    }
    /// Compute A-weighting in dB at frequency f (Hz).
    #[allow(dead_code)]
    pub fn magnitude_db(&self, f: f64) -> f64 {
        let mag = self.magnitude(f);
        if mag < 1e-30 {
            return -200.0;
        }
        20.0 * mag.log10()
    }
    /// Apply A-weighting correction to a spectrum.
    ///
    /// `freqs` and `spls` must have the same length.
    /// Returns A-weighted SPL values.
    #[allow(dead_code)]
    pub fn apply_to_spectrum(&self, freqs: &[f64], spls: &[f64]) -> Vec<f64> {
        freqs
            .iter()
            .zip(spls.iter())
            .map(|(&f, &spl)| spl + self.magnitude_db(f))
            .collect()
    }
}
/// FWH surface integral data: surface normal, velocity, and pressure at one panel.
#[allow(dead_code)]
pub struct FwhSurfacePanel {
    /// Panel area (m^2).
    pub area: f64,
    /// Outward unit normal \[nx, ny, nz\].
    pub normal: [f64; 3],
    /// Surface velocity \[ux, uy, uz\] (m/s).
    pub velocity: [f64; 3],
    /// Surface pressure fluctuation (Pa).
    pub delta_p: f64,
    /// Observer-to-source unit vector \[rx, ry, rz\].
    pub r_hat: [f64; 3],
    /// Observer distance (m).
    pub r: f64,
}
/// Source type identification based on Mach number scaling.
#[allow(dead_code)]
#[derive(Debug, Clone, PartialEq)]
pub enum AcousticSourceType {
    /// W ∝ Ma² (compact volume fluctuation)
    Monopole,
    /// W ∝ Ma⁴ (unsteady force)
    Dipole,
    /// W ∝ Ma⁸ (turbulent flow, Lighthill 8th power)
    Quadrupole,
}
/// A single FW-H (Ffowcs Williams-Hawkings) surface panel.
#[allow(dead_code)]
#[derive(Clone)]
pub struct FwhPanel {
    /// Panel center position \[x, y, z\].
    pub center: [f64; 3],
    /// Outward unit normal \[nx, ny, nz\].
    pub normal: [f64; 3],
    /// Panel area.
    pub area: f64,
    /// Surface pressure (gauge).
    pub pressure: f64,
    /// Surface velocity (normal component of fluid velocity at surface).
    pub un: f64,
}
/// LBM near-field to far-field pressure extrapolation.
///
/// Uses FWH surface integral on a near-field control surface surrounding
/// the source region. The surface encloses the near-field LBM data.
#[allow(dead_code)]
pub struct NearToFarField {
    /// Control surface panels (reuse FwhPanel)
    pub panels: Vec<FwhPanel>,
    /// Observer positions for far-field calculation
    pub observers: Vec<[f64; 3]>,
}
impl NearToFarField {
    /// Create a new near-to-far-field transformer.
    #[allow(dead_code)]
    pub fn new() -> Self {
        Self {
            panels: Vec::new(),
            observers: Vec::new(),
        }
    }
    /// Add a control surface panel.
    #[allow(dead_code)]
    pub fn add_panel(&mut self, panel: FwhPanel) {
        self.panels.push(panel);
    }
    /// Add an observer position.
    #[allow(dead_code)]
    pub fn add_observer(&mut self, pos: [f64; 3]) {
        self.observers.push(pos);
    }
    /// Compute far-field pressure at all observer positions.
    #[allow(dead_code)]
    pub fn compute_far_field(&self) -> Vec<f64> {
        self.observers
            .iter()
            .map(|&obs| fwh_acoustic_pressure(&self.panels, obs))
            .collect()
    }
    /// Compute SPL at all observer positions.
    #[allow(dead_code)]
    pub fn compute_spl(&self) -> Vec<f64> {
        self.compute_far_field()
            .iter()
            .map(|&p| spl_from_pressure_fluctuation(p.abs()))
            .collect()
    }
}
/// Helmholtz resonator model with volume, neck geometry and speed of sound.
#[allow(dead_code)]
pub struct HelmholtzResonator {
    /// Cavity volume (m³ or lattice units³)
    pub volume: f64,
    /// Neck length (m or lattice units)
    pub neck_length: f64,
    /// Neck cross-sectional area (m² or lattice units²)
    pub neck_area: f64,
    /// Speed of sound
    pub c0: f64,
}
impl HelmholtzResonator {
    /// Create a new Helmholtz resonator.
    #[allow(dead_code)]
    pub fn new(v: f64, l: f64, a: f64, c: f64) -> Self {
        Self {
            volume: v,
            neck_length: l,
            neck_area: a,
            c0: c,
        }
    }
    /// Resonance frequency with end correction:
    /// f = c/(2π) * sqrt(A / (V*(L + 0.85*sqrt(A/π))))
    #[allow(dead_code)]
    pub fn resonance_frequency(&self) -> f64 {
        let l_eff = self.neck_length + 0.85 * (self.neck_area / PI).sqrt();
        self.c0 / (2.0 * PI) * (self.neck_area / (self.volume * l_eff)).sqrt()
    }
    /// Quality factor: Q = sqrt(V*L) / (A^0.25 * sqrt(2*nu/omega))
    ///
    /// `viscosity` is the kinematic viscosity ν.
    #[allow(dead_code)]
    pub fn quality_factor(&self, viscosity: f64) -> f64 {
        let omega = 2.0 * PI * self.resonance_frequency();
        let denom = self.neck_area.powf(0.25) * (2.0 * viscosity / omega).sqrt();
        (self.volume * self.neck_length).sqrt() / denom
    }
    /// Peak transmission loss: TL_max ≈ 10*log10(1 + V/(4*L*A))
    #[allow(dead_code)]
    pub fn transmission_loss_peak(&self) -> f64 {
        10.0 * (1.0 + self.volume / (4.0 * self.neck_length * self.neck_area)).log10()
    }
}
/// Observer position for FWH calculation.
#[allow(dead_code)]
#[derive(Clone, Debug)]
pub struct FwhObserver {
    /// Observer position \[x, y, z\]
    pub position: [f64; 3],
    /// Speed of sound at observer location
    pub c0: f64,
    /// Mean fluid density at observer location
    pub rho0: f64,
}
impl FwhObserver {
    /// Create a new FWH observer.
    #[allow(dead_code)]
    pub fn new(x: f64, y: f64, z: f64, c0: f64, rho0: f64) -> Self {
        Self {
            position: [x, y, z],
            c0,
            rho0,
        }
    }
    /// Compute distance from this observer to a source point.
    #[allow(dead_code)]
    pub fn distance_to(&self, source: [f64; 3]) -> f64 {
        let dx = self.position[0] - source[0];
        let dy = self.position[1] - source[1];
        let dz = self.position[2] - source[2];
        (dx * dx + dy * dy + dz * dz).sqrt()
    }
    /// Unit vector from source to observer.
    #[allow(dead_code)]
    pub fn unit_vector_from(&self, source: [f64; 3]) -> [f64; 3] {
        let dx = self.position[0] - source[0];
        let dy = self.position[1] - source[1];
        let dz = self.position[2] - source[2];
        let r = (dx * dx + dy * dy + dz * dz).sqrt();
        if r < 1e-15 {
            return [0.0, 0.0, 0.0];
        }
        [dx / r, dy / r, dz / r]
    }
}
