//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::needless_range_loop)]
use std::f64::consts::PI;

use super::functions::vec3_norm;
#[allow(unused_imports)]
use super::functions::*;

/// Regular and irregular wave load models.
///
/// Supports Airy (linear) waves and JONSWAP spectral models.
pub struct WaveLoading {
    /// Significant wave height Hs (m).
    pub hs: f64,
    /// Peak spectral period Tp (s).
    pub tp: f64,
    /// Water depth d (m).
    pub water_depth: f64,
    /// Fluid density ρ (kg/m³).
    pub fluid_density: f64,
    /// Gravitational acceleration g (m/s²).
    pub gravity: f64,
}
impl WaveLoading {
    /// Create a wave loading model.
    pub fn new(hs: f64, tp: f64, water_depth: f64) -> Self {
        WaveLoading {
            hs,
            tp,
            water_depth,
            fluid_density: 1025.0,
            gravity: 9.81,
        }
    }
    /// Peak angular frequency ωp = 2π / Tp (rad/s).
    pub fn peak_frequency(&self) -> f64 {
        2.0 * PI / self.tp
    }
    /// Airy wave horizontal velocity at depth z (m below surface) and time t (s).
    /// Uses linear wave theory: u = (π*H/T) * cosh(k*(d+z)) / sinh(k*d) * cos(kx - ωt).
    pub fn airy_velocity(&self, z: f64, t: f64, x: f64) -> f64 {
        let h = self.hs;
        let omega = self.peak_frequency();
        let k = self.wave_number();
        let d = self.water_depth;
        let num = PI * h / self.tp;
        let depth_factor = ((k * (d + z)).cosh()) / ((k * d).sinh().max(1e-10));
        num * depth_factor * (k * x - omega * t).cos()
    }
    /// Dispersion relation: ω² = g * k * tanh(k * d). Solved iteratively.
    pub fn wave_number(&self) -> f64 {
        let omega = self.peak_frequency();
        let g = self.gravity;
        let d = self.water_depth;
        let mut k = omega * omega / g;
        for _ in 0..50 {
            let f = g * k * (k * d).tanh() - omega * omega;
            let df = g * (k * d).tanh() + g * k * d / (k * d).cosh().powi(2);
            k -= f / df.max(1e-10);
            if f.abs() < 1e-12 {
                break;
            }
        }
        k
    }
    /// JONSWAP spectral density S(ω) at angular frequency ω.
    pub fn jonswap_spectrum(&self, omega: f64) -> f64 {
        let wp = self.peak_frequency();
        let hs = self.hs;
        let g = self.gravity;
        let gamma: f64 = 3.3;
        let alpha: f64 = 0.0081;
        let sigma = if omega <= wp { 0.07 } else { 0.09 };
        let exp_term = (-5.0 / 4.0 * (wp / omega).powi(4)).exp();
        let s_pm = alpha * g * g / omega.powi(5) * exp_term;
        let r = (-(omega - wp).powi(2) / (2.0 * sigma * sigma * wp * wp)).exp();
        let _hs_scale = hs * hs;
        s_pm * gamma.powf(r)
    }
}
/// Lagrangian marker point for the immersed boundary method.
#[derive(Clone, Debug)]
pub struct IbMarker {
    /// Marker position (m).
    pub position: [f64; 3],
    /// Target (desired) position (m) — for elastic IB forces.
    pub target_position: [f64; 3],
    /// Stiffness coefficient κ for the IB spring force (N/m).
    pub stiffness: f64,
    /// Force applied from IB to fluid (N).
    pub ib_force: [f64; 3],
}
impl IbMarker {
    /// Create an IB marker at a given position.
    pub fn new(position: [f64; 3], stiffness: f64) -> Self {
        Self {
            position,
            target_position: position,
            stiffness,
            ib_force: [0.0; 3],
        }
    }
    /// Compute the IB spring force: F = κ * (X_target - X).
    pub fn compute_force(&mut self) {
        for i in 0..3 {
            self.ib_force[i] = self.stiffness * (self.target_position[i] - self.position[i]);
        }
    }
    /// Displacement of marker from target (m).
    pub fn displacement(&self) -> f64 {
        let d: [f64; 3] = [
            self.target_position[0] - self.position[0],
            self.target_position[1] - self.position[1],
            self.target_position[2] - self.position[2],
        ];
        vec3_norm(d)
    }
}
/// Hydrodynamic drag on a body moving through a fluid.
pub struct HydrodynamicDrag {
    /// Fluid density ρ (kg/m³).
    pub fluid_density: f64,
    /// Reference area A_ref (m²).
    pub reference_area: f64,
    /// Drag coefficient Cd (dimensionless, depends on shape and Re).
    pub cd: f64,
    /// Viscosity μ (Pa·s) for Reynolds number calculation.
    pub viscosity: f64,
    /// Characteristic length L (m).
    pub char_length: f64,
}
impl HydrodynamicDrag {
    /// Create a new hydrodynamic drag model.
    pub fn new(fluid_density: f64, reference_area: f64, cd: f64) -> Self {
        HydrodynamicDrag {
            fluid_density,
            reference_area,
            cd,
            viscosity: 1e-3,
            char_length: 1.0,
        }
    }
    /// Reynolds number: Re = ρ * v * L / μ.
    pub fn reynolds_number(&self, velocity: f64) -> f64 {
        self.fluid_density * velocity * self.char_length / self.viscosity
    }
    /// Drag force magnitude F_d = 0.5 * ρ * v² * Cd * A (N).
    pub fn drag_force(&self, velocity: f64) -> f64 {
        0.5 * self.fluid_density * velocity * velocity * self.cd * self.reference_area
    }
    /// Drag vector opposing motion.
    pub fn drag_vector(&self, velocity_vec: [f64; 3]) -> [f64; 3] {
        let v_mag = vec3_norm(velocity_vec);
        if v_mag < 1e-15 {
            return [0.0; 3];
        }
        let fd = self.drag_force(v_mag);
        [
            -fd * velocity_vec[0] / v_mag,
            -fd * velocity_vec[1] / v_mag,
            -fd * velocity_vec[2] / v_mag,
        ]
    }
    /// Cd vs Reynolds number (simplified empirical for a sphere).
    pub fn cd_from_reynolds(re: f64) -> f64 {
        if re < 0.5 {
            24.0 / re.max(1e-10)
        } else if re < 1000.0 {
            24.0 / re + 6.0 / (1.0 + re.sqrt()) + 0.4
        } else {
            0.44
        }
    }
}
/// Wind pressure on structures (ASCE 7 log-law model).
pub struct WindPressure {
    /// Basic wind speed V (m/s) at 10 m elevation.
    pub v_ref: f64,
    /// Air density ρ (kg/m³).
    pub air_density: f64,
    /// Surface roughness length z0 (m).
    pub roughness_length: f64,
    /// Exposure category (A, B, C, D).
    pub exposure_category: char,
    /// Shape/pressure coefficient Cp for the structure face.
    pub cp: f64,
    /// Gust factor G.
    pub gust_factor: f64,
}
impl WindPressure {
    /// Create a wind pressure model for an urban environment (Exposure B).
    pub fn new(v_ref: f64) -> Self {
        WindPressure {
            v_ref,
            air_density: 1.225,
            roughness_length: 0.3,
            exposure_category: 'B',
            cp: 0.8,
            gust_factor: 0.85,
        }
    }
    /// Wind speed at height z (m) using log-law: V(z) = V_ref * ln(z/z0) / ln(10/z0).
    pub fn wind_speed_at(&self, z: f64) -> f64 {
        if z <= 0.0 {
            return 0.0;
        }
        let z_safe = z.max(self.roughness_length);
        let z0 = self.roughness_length;
        self.v_ref * (z_safe / z0).ln() / (10.0 / z0).ln()
    }
    /// Dynamic pressure q(z) = 0.5 * ρ * V(z)² (Pa).
    pub fn dynamic_pressure(&self, z: f64) -> f64 {
        let v = self.wind_speed_at(z);
        0.5 * self.air_density * v * v
    }
    /// Design wind pressure p(z) = G * Cp * q(z) (Pa).
    pub fn design_pressure(&self, z: f64) -> f64 {
        self.gust_factor * self.cp * self.dynamic_pressure(z)
    }
    /// Wind force on a surface of area A at height z (N).
    pub fn wind_force(&self, z: f64, area: f64) -> f64 {
        self.design_pressure(z) * area
    }
}
/// Immersed boundary (IB) body represented by a collection of Lagrangian markers.
pub struct ImmersedBoundaryBody {
    /// Lagrangian markers on the body surface.
    pub markers: Vec<IbMarker>,
    /// Fluid grid spacing Δx (m) — used for spreading/interpolation.
    pub grid_spacing: f64,
    /// Fluid density ρ (kg/m³).
    pub fluid_density: f64,
}
impl ImmersedBoundaryBody {
    /// Create an IB body.
    pub fn new(grid_spacing: f64, fluid_density: f64) -> Self {
        Self {
            markers: Vec::new(),
            grid_spacing,
            fluid_density,
        }
    }
    /// Add a marker to the body.
    pub fn add_marker(&mut self, position: [f64; 3], stiffness: f64) {
        self.markers.push(IbMarker::new(position, stiffness));
    }
    /// Compute all IB forces.
    pub fn compute_forces(&mut self) {
        for m in &mut self.markers {
            m.compute_force();
        }
    }
    /// Total IB force on the fluid (sum over all markers).
    pub fn total_force(&self) -> [f64; 3] {
        let mut f = [0.0f64; 3];
        for m in &self.markers {
            for i in 0..3 {
                f[i] += m.ib_force[i];
            }
        }
        f
    }
    /// Regularised delta function (Roma-Peskin 3-pt), spread from marker to grid.
    ///
    /// `r` = distance from marker to grid node (normalised by Δx).
    pub fn delta_function(r: f64) -> f64 {
        let ra = r.abs();
        if ra <= 0.5 {
            (1.0 + (4.0 * ra * ra).sqrt()) / 3.0
        } else if ra <= 1.5 {
            (5.0 - 3.0 * ra - (1.0 - 3.0 * (1.0 - ra) * (1.0 - ra)).sqrt()) / 6.0
        } else {
            0.0
        }
    }
    /// Number of markers.
    pub fn n_markers(&self) -> usize {
        self.markers.len()
    }
}
/// Random sea state described by a JONSWAP spectrum, for computing
/// statistical wave loads on offshore structures.
pub struct RandomSeaState {
    /// Significant wave height Hs (m).
    pub hs: f64,
    /// Peak period Tp (s).
    pub tp: f64,
    /// JONSWAP peak enhancement factor γ (typically 3.3).
    pub gamma: f64,
    /// Water depth d (m).
    pub water_depth: f64,
    /// Fluid density ρ (kg/m³).
    pub fluid_density: f64,
    /// Gravity g (m/s²).
    pub gravity: f64,
    /// Number of frequency components for irregular wave synthesis.
    pub n_components: usize,
}
impl RandomSeaState {
    /// Create a random sea state with JONSWAP parameters.
    pub fn new(hs: f64, tp: f64, water_depth: f64) -> Self {
        Self {
            hs,
            tp,
            gamma: 3.3,
            water_depth,
            fluid_density: 1025.0,
            gravity: 9.81,
            n_components: 100,
        }
    }
    /// Peak angular frequency ωp (rad/s).
    pub fn omega_p(&self) -> f64 {
        2.0 * PI / self.tp
    }
    /// JONSWAP spectral density S(ω) (m²·s/rad).
    pub fn jonswap(&self, omega: f64) -> f64 {
        let wp = self.omega_p();
        let g = self.gravity;
        let alpha = 5.0 * wp.powi(4) * self.hs.powi(2) / (16.0 * g * g);
        let sigma: f64 = if omega <= wp { 0.07 } else { 0.09 };
        let exp1 = (-5.0 / 4.0 * (wp / omega.max(1e-10)).powi(4)).exp();
        let r = (-(omega - wp).powi(2) / (2.0 * sigma.powi(2) * wp.powi(2))).exp();
        alpha * g * g / omega.max(1e-10).powi(5) * exp1 * self.gamma.powf(r)
    }
    /// Approximate zero-crossing period Tz from JONSWAP parameters.
    ///
    /// `Tz ≈ Tp / 1.408` for γ = 3.3.
    pub fn zero_crossing_period(&self) -> f64 {
        self.tp / 1.408
    }
    /// Most probable maximum wave height in `n_waves` waves.
    ///
    /// `H_max ≈ Hs * sqrt(ln(n_waves) / 2)` (Rayleigh distribution).
    pub fn most_probable_max_wave(&self, n_waves: f64) -> f64 {
        self.hs * (n_waves.ln() / 2.0).sqrt()
    }
    /// Spectral moment m_n (∫ ωⁿ S(ω) dω) estimated over frequency range.
    pub fn spectral_moment(&self, n: i32, omega_lo: f64, omega_hi: f64, n_pts: usize) -> f64 {
        let dw = (omega_hi - omega_lo) / n_pts as f64;
        let mut sum = 0.0;
        for i in 0..n_pts {
            let w = omega_lo + (i as f64 + 0.5) * dw;
            sum += w.powi(n) * self.jonswap(w) * dw;
        }
        sum
    }
    /// Significant wave height from spectral moment m0.
    ///
    /// `Hs = 4 * sqrt(m0)`.
    pub fn hs_from_spectrum(&self, omega_lo: f64, omega_hi: f64, n_pts: usize) -> f64 {
        let m0 = self.spectral_moment(0, omega_lo, omega_hi, n_pts);
        4.0 * m0.sqrt()
    }
    /// Peak wave force on a vertical cylinder using Morison equation with spectral approach.
    ///
    /// Returns the RMS force per unit length (N/m).
    pub fn rms_morison_force(&self, diameter: f64, cm: f64, cd: f64) -> f64 {
        let wp = self.omega_p();
        let s_wp = self.jonswap(wp);
        let k = wp * wp / self.gravity;
        let eta_rms = s_wp.sqrt();
        let u_rms = wp * eta_rms;
        let inertia = self.fluid_density * cm * PI / 4.0 * diameter.powi(2) * wp * u_rms;
        let drag = 0.5 * self.fluid_density * cd * diameter * u_rms.powi(2);
        let _ = k;
        (inertia.powi(2) + drag.powi(2)).sqrt()
    }
    /// Synthesise a wave surface elevation time series using linear superposition.
    ///
    /// Returns a vector of (time, elevation) pairs.
    pub fn synthesise(
        &self,
        t_start: f64,
        t_end: f64,
        dt: f64,
        omega_lo: f64,
        omega_hi: f64,
    ) -> Vec<[f64; 2]> {
        let dw = (omega_hi - omega_lo) / self.n_components as f64;
        let n_t = ((t_end - t_start) / dt).ceil() as usize;
        let mut result = Vec::with_capacity(n_t);
        let components: Vec<(f64, f64, f64)> = (0..self.n_components)
            .map(|i| {
                let w = omega_lo + (i as f64 + 0.5) * dw;
                let s = self.jonswap(w);
                let amp = (2.0 * s * dw).sqrt();
                let phase = (i as f64 * 0.31416 * 7.0) % (2.0 * PI);
                (w, amp, phase)
            })
            .collect();
        for it in 0..n_t {
            let t = t_start + it as f64 * dt;
            let eta: f64 = components
                .iter()
                .map(|(w, amp, phi)| amp * (w * t - phi).cos())
                .sum();
            result.push([t, eta]);
        }
        result
    }
}
/// Galloping instability (Den Hartog criterion) for bluff bodies.
pub struct GallopingAnalysis {
    /// Air density ρ (kg/m³).
    pub air_density: f64,
    /// Mean wind speed U (m/s).
    pub wind_speed: f64,
    /// Cross-section dimension D (m) — across-wind.
    pub dimension: f64,
    /// Mass per unit length m (kg/m).
    pub mass_per_length: f64,
    /// Structural damping ratio ζ.
    pub damping_ratio: f64,
    /// Natural frequency ωn (rad/s).
    pub natural_frequency: f64,
    /// Slope of lift + drag coefficient dCl/dα + Cd at α=0 (Den Hartog coefficient A1).
    pub den_hartog_coeff: f64,
}
impl GallopingAnalysis {
    /// Create a galloping analysis model.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        air_density: f64,
        wind_speed: f64,
        dimension: f64,
        mass_per_length: f64,
        damping_ratio: f64,
        natural_frequency: f64,
        den_hartog_coeff: f64,
    ) -> Self {
        Self {
            air_density,
            wind_speed,
            dimension,
            mass_per_length,
            damping_ratio,
            natural_frequency,
            den_hartog_coeff,
        }
    }
    /// Reduced velocity Ur = U / (fn * D).
    pub fn reduced_velocity(&self) -> f64 {
        let fn_hz = self.natural_frequency / (2.0 * PI);
        self.wind_speed / (fn_hz * self.dimension).max(1e-15)
    }
    /// Den Hartog galloping criterion: unstable if A1 < 0 and |A1| > 4*ζ*m/(ρ*D²).
    ///
    /// Returns `true` if galloping is predicted.
    pub fn is_galloping_unstable(&self) -> bool {
        let rho = self.air_density;
        let d = self.dimension;
        let m = self.mass_per_length;
        let zeta = self.damping_ratio;
        if self.den_hartog_coeff >= 0.0 {
            return false;
        }
        let threshold = 4.0 * zeta * m / (rho * d).max(1e-15);
        (-self.den_hartog_coeff) > threshold
    }
    /// Galloping amplitude at steady state (Parkinson-Smith model, leading term, m).
    ///
    /// Approximate: `Y_0 = D * sqrt(-A1 / A3)` where A3 > 0.
    pub fn galloping_amplitude(&self, a3: f64) -> f64 {
        if self.den_hartog_coeff >= 0.0 || a3 <= 0.0 {
            return 0.0;
        }
        self.dimension * (-self.den_hartog_coeff / a3).sqrt()
    }
    /// Critical wind speed for galloping onset (m/s).
    ///
    /// `U_crit = 4 * ζ * m * ωn / (ρ * D * |A1|)`.
    pub fn critical_wind_speed(&self) -> f64 {
        let rho = self.air_density;
        let d = self.dimension;
        let m = self.mass_per_length;
        let zeta = self.damping_ratio;
        let omega_n = self.natural_frequency;
        let a1 = self.den_hartog_coeff.abs().max(1e-15);
        4.0 * zeta * m * omega_n / (rho * d * a1)
    }
}
/// Flutter analysis using the p-k method (simplified 2-DOF aeroelastic).
pub struct FlutterAnalysis {
    /// Bending stiffness EI (N·m²).
    pub ei: f64,
    /// Torsional stiffness GJ (N·m²/rad).
    pub gj: f64,
    /// Mass per span m (kg/m).
    pub mass_per_span: f64,
    /// Moment of inertia per span Iα (kg·m).
    pub i_alpha_per_span: f64,
    /// Semi-chord b (m).
    pub semi_chord: f64,
    /// Span L (m).
    pub span: f64,
    /// Air density ρ (kg/m³).
    pub air_density: f64,
    /// Static unbalance Sα = m * x_α (kg).
    pub static_unbalance: f64,
}
impl FlutterAnalysis {
    /// Create a flutter analysis model.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        ei: f64,
        gj: f64,
        mass_per_span: f64,
        i_alpha_per_span: f64,
        semi_chord: f64,
        span: f64,
        air_density: f64,
        static_unbalance: f64,
    ) -> Self {
        Self {
            ei,
            gj,
            mass_per_span,
            i_alpha_per_span,
            semi_chord,
            span,
            air_density,
            static_unbalance,
        }
    }
    /// Uncoupled bending frequency ωh (rad/s).
    pub fn bending_frequency(&self) -> f64 {
        let m = self.mass_per_span * self.span;
        (3.0 * self.ei / (m * self.span.powi(3))).sqrt()
    }
    /// Uncoupled torsional frequency ωα (rad/s).
    pub fn torsional_frequency(&self) -> f64 {
        let i_total = self.i_alpha_per_span * self.span;
        (self.gj / (i_total * self.span)).sqrt()
    }
    /// Mass ratio μ = m / (π * ρ * b²).
    pub fn mass_ratio(&self) -> f64 {
        self.mass_per_span / (PI * self.air_density * self.semi_chord.powi(2))
    }
    /// Frequency ratio r_f = ωh / ωα.
    pub fn frequency_ratio(&self) -> f64 {
        self.bending_frequency() / self.torsional_frequency().max(1e-15)
    }
    /// Simplified flutter speed estimate (Frazer-Duncan-Collar method, m/s).
    ///
    /// `U_flutter ≈ ωα * b * sqrt(μ)`.
    pub fn flutter_speed(&self) -> f64 {
        let wa = self.torsional_frequency();
        let mu = self.mass_ratio();
        wa * self.semi_chord * mu.sqrt()
    }
    /// Reduced flutter speed Vf* = U_flutter / (ωα * b).
    pub fn reduced_flutter_speed(&self) -> f64 {
        self.flutter_speed() / (self.torsional_frequency() * self.semi_chord).max(1e-15)
    }
    /// Check if static unbalance destabilises flutter.
    pub fn is_destabilised_by_unbalance(&self) -> bool {
        self.static_unbalance > 0.0
    }
}
/// Along-wind, across-wind and torsional response of tall buildings.
pub struct TallBuildingWindResponse {
    /// Building height H (m).
    pub height: f64,
    /// Building width B (m) — across-wind.
    pub width: f64,
    /// Building depth D (m) — along-wind.
    pub depth: f64,
    /// Fundamental natural frequency f1 (Hz).
    pub natural_frequency: f64,
    /// Structural damping ratio ζ.
    pub damping_ratio: f64,
    /// Design wind speed V_H at roof height (m/s).
    pub wind_speed_at_roof: f64,
    /// Air density ρ (kg/m³).
    pub air_density: f64,
    /// Total building mass M (kg).
    pub total_mass: f64,
    /// Mean drag coefficient Cd.
    pub cd: f64,
}
impl TallBuildingWindResponse {
    /// Create a tall building wind response model.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        height: f64,
        width: f64,
        depth: f64,
        natural_frequency: f64,
        damping_ratio: f64,
        wind_speed_at_roof: f64,
        total_mass: f64,
    ) -> Self {
        Self {
            height,
            width,
            depth,
            natural_frequency,
            damping_ratio,
            wind_speed_at_roof,
            air_density: 1.225,
            total_mass,
            cd: 1.3,
        }
    }
    /// Mean along-wind base shear (N).
    ///
    /// `F_mean = 0.5 * ρ * V_H² * Cd * B * H`.
    pub fn mean_along_wind_force(&self) -> f64 {
        0.5 * self.air_density
            * self.wind_speed_at_roof.powi(2)
            * self.cd
            * self.width
            * self.height
    }
    /// Mean along-wind base overturning moment (N·m).
    pub fn mean_overturning_moment(&self) -> f64 {
        self.mean_along_wind_force() * 2.0 / 3.0 * self.height
    }
    /// Gust factor G for along-wind response (ASCE 7 simplified).
    ///
    /// `G ≈ 1 + 2 * g * σ_u / V_H`.
    pub fn gust_factor(&self, turbulence_intensity: f64) -> f64 {
        let g = 3.4;
        1.0 + 2.0 * g * turbulence_intensity
    }
    /// Design along-wind base shear including gust effects (N).
    pub fn design_along_wind_force(&self, turbulence_intensity: f64) -> f64 {
        self.mean_along_wind_force() * self.gust_factor(turbulence_intensity)
    }
    /// Across-wind RMS base shear using empirical formula (N).
    ///
    /// `F_rms = 0.5 * ρ * V_H² * B * H * CL_rms`.
    pub fn across_wind_rms_force(&self, cl_rms: f64) -> f64 {
        0.5 * self.air_density * self.wind_speed_at_roof.powi(2) * self.width * self.height * cl_rms
    }
    /// Peak across-wind displacement at roof (m) using random vibration theory.
    ///
    /// `y_peak = g * σ_y` where `σ_y = F_rms / (M * ω_n² * 2 * ζ)`.
    pub fn peak_across_wind_displacement(&self, cl_rms: f64) -> f64 {
        let g_peak = 3.5;
        let omega_n = 2.0 * PI * self.natural_frequency;
        let f_rms = self.across_wind_rms_force(cl_rms);
        let sigma_y =
            f_rms / (self.total_mass * omega_n.powi(2) * 2.0 * self.damping_ratio).max(1e-30);
        g_peak * sigma_y
    }
    /// Strouhal shedding frequency for rectangular cross-section (Hz).
    ///
    /// Uses St ≈ 0.1 for tall buildings.
    pub fn shedding_frequency(&self) -> f64 {
        let st = 0.1;
        st * self.wind_speed_at_roof / self.width
    }
    /// Check if across-wind resonance is likely (lock-in condition).
    pub fn is_across_wind_resonant(&self) -> bool {
        let fs = self.shedding_frequency();
        (self.natural_frequency / fs - 1.0).abs() < 0.2
    }
    /// Torsional base moment from asymmetric wind loading (N·m).
    ///
    /// Simplified: `M_T = e * F_along` where `e ≈ 0.15 * B`.
    pub fn torsional_moment(&self, eccentricity_ratio: f64) -> f64 {
        let e = eccentricity_ratio * self.width;
        self.mean_along_wind_force() * e
    }
    /// Reduced natural frequency f1 * H / V_H (dimensionless).
    pub fn reduced_natural_frequency(&self) -> f64 {
        self.natural_frequency * self.height / self.wind_speed_at_roof.max(1e-15)
    }
    /// Peak accelerations at roof level (m/s²).
    ///
    /// `a_peak = g * σ_a = g * ω_n² * σ_y`.
    pub fn peak_acceleration(&self, cl_rms: f64) -> f64 {
        let omega_n = 2.0 * PI * self.natural_frequency;
        let disp = self.peak_across_wind_displacement(cl_rms);
        omega_n.powi(2) * disp
    }
}
/// ALE domain with a collection of mesh nodes.
pub struct AleDomain {
    /// Mesh nodes.
    pub nodes: Vec<AleMeshNode>,
    /// Domain identifier.
    pub domain_id: u64,
}
impl AleDomain {
    /// Create an ALE domain.
    pub fn new(domain_id: u64) -> Self {
        Self {
            nodes: Vec::new(),
            domain_id,
        }
    }
    /// Add a node.
    pub fn add_node(&mut self, position: [f64; 3]) {
        self.nodes.push(AleMeshNode::new(position));
    }
    /// Advance all nodes in time.
    pub fn advance(&mut self, dt: f64) {
        for n in &mut self.nodes {
            n.advance_mesh(dt);
        }
    }
    /// Maximum mesh velocity magnitude.
    pub fn max_mesh_velocity(&self) -> f64 {
        self.nodes
            .iter()
            .map(|n| vec3_norm(n.mesh_velocity))
            .fold(0.0_f64, f64::max)
    }
    /// Number of nodes.
    pub fn n_nodes(&self) -> usize {
        self.nodes.len()
    }
}
/// Beam in cross-flow: hydroelastic response model.
///
/// Models a flexible beam subjected to hydrodynamic loading using
/// Euler-Bernoulli beam theory with added mass and damping.
pub struct HydroelasticBeam {
    /// Beam length L (m).
    pub length: f64,
    /// Bending stiffness EI (N·m²).
    pub ei: f64,
    /// Mass per unit length m (kg/m).
    pub mass_per_length: f64,
    /// Added mass per unit length m_a (kg/m).
    pub added_mass_per_length: f64,
    /// Structural damping ratio ζ.
    pub damping_ratio: f64,
    /// Fluid velocity U (m/s).
    pub fluid_velocity: f64,
    /// Fluid density ρ (kg/m³).
    pub fluid_density: f64,
    /// Cross-section diameter D (m).
    pub diameter: f64,
    /// Drag coefficient Cd.
    pub cd: f64,
}
impl HydroelasticBeam {
    /// Create a hydroelastic beam model.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        length: f64,
        ei: f64,
        mass_per_length: f64,
        fluid_density: f64,
        diameter: f64,
    ) -> Self {
        let added_mass_per_length = fluid_density * PI / 4.0 * diameter.powi(2);
        Self {
            length,
            ei,
            mass_per_length,
            added_mass_per_length,
            damping_ratio: 0.02,
            fluid_velocity: 0.0,
            fluid_density,
            diameter,
            cd: 1.0,
        }
    }
    /// Total mass per unit length including added mass (kg/m).
    pub fn total_mass_per_length(&self) -> f64 {
        self.mass_per_length + self.added_mass_per_length
    }
    /// Natural frequency of nth mode for a cantilevered beam (rad/s).
    ///
    /// Uses exact cantilever eigenvalues: β_n * L = \[1.875, 4.694, 7.855, ...\].
    pub fn natural_frequency_mode(&self, mode: u32) -> f64 {
        let beta_l = match mode {
            1 => 1.8751,
            2 => 4.6941,
            3 => 7.8548,
            4 => 10.9955,
            _ => PI * (mode as f64 - 0.5),
        };
        let beta = beta_l / self.length;
        let m_eff = self.total_mass_per_length();
        (beta.powi(4) * self.ei / m_eff).sqrt()
    }
    /// Added mass coefficient Ca = m_a / (ρ * A).
    pub fn added_mass_coefficient(&self) -> f64 {
        let rho_a = PI / 4.0 * self.diameter.powi(2) * self.fluid_density;
        self.added_mass_per_length / rho_a.max(1e-30)
    }
    /// Static tip deflection under uniform drag load (m).
    ///
    /// `δ = q * L⁴ / (8 * EI)` where `q = 0.5 * ρ * U² * Cd * D`.
    pub fn static_tip_deflection(&self) -> f64 {
        let q = 0.5 * self.fluid_density * self.fluid_velocity.powi(2) * self.cd * self.diameter;
        q * self.length.powi(4) / (8.0 * self.ei).max(1e-30)
    }
    /// Dynamic magnification factor at frequency ω.
    pub fn dynamic_magnification(&self, omega: f64) -> f64 {
        let omega_n = self.natural_frequency_mode(1);
        let r = omega / omega_n.max(1e-15);
        let zeta = self.damping_ratio;
        1.0 / ((1.0 - r.powi(2)).powi(2) + (2.0 * zeta * r).powi(2)).sqrt()
    }
    /// Reduced velocity Ur = U / (f1 * D).
    pub fn reduced_velocity(&self) -> f64 {
        let f1 = self.natural_frequency_mode(1) / (2.0 * PI);
        self.fluid_velocity / (f1 * self.diameter).max(1e-15)
    }
}
/// VIV lock-in analysis with mass-damping parameter (Skop-Griffin model).
pub struct VivLockIn {
    /// Strouhal number St.
    pub strouhal: f64,
    /// Cylinder diameter D (m).
    pub diameter: f64,
    /// Mass ratio m* = m / (ρ D² L) (dimensionless).
    pub mass_ratio: f64,
    /// Structural damping ratio ζ.
    pub damping_ratio: f64,
    /// Fluid velocity U (m/s).
    pub fluid_velocity: f64,
    /// Natural frequency fn (Hz).
    pub natural_frequency_hz: f64,
    /// Maximum lift coefficient CL in lock-in.
    pub cl_max: f64,
}
impl VivLockIn {
    /// Create a VIV lock-in model.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        diameter: f64,
        mass_ratio: f64,
        damping_ratio: f64,
        fluid_velocity: f64,
        natural_frequency_hz: f64,
    ) -> Self {
        Self {
            strouhal: 0.2,
            diameter,
            mass_ratio,
            damping_ratio,
            fluid_velocity,
            natural_frequency_hz,
            cl_max: 0.3,
        }
    }
    /// Vortex shedding frequency fs (Hz).
    pub fn shedding_frequency(&self) -> f64 {
        self.strouhal * self.fluid_velocity / self.diameter
    }
    /// Reduced velocity Ur = U / (fn * D).
    pub fn reduced_velocity(&self) -> f64 {
        self.fluid_velocity / (self.natural_frequency_hz * self.diameter).max(1e-15)
    }
    /// Mass-damping parameter (Scruton number) Sc = m* * ζ.
    pub fn scruton_number(&self) -> f64 {
        self.mass_ratio * self.damping_ratio
    }
    /// Skop-Griffin parameter SG = 2 * π³ * St² * Sc.
    pub fn skop_griffin_parameter(&self) -> f64 {
        2.0 * PI.powi(3) * self.strouhal.powi(2) * self.scruton_number()
    }
    /// Griffin plot: maximum VIV amplitude A/D from Skop-Griffin parameter.
    ///
    /// `A/D = 1.29 / (1 + 0.43 * SG)^(0.54)`.
    pub fn max_amplitude_ratio(&self) -> f64 {
        let sg = self.skop_griffin_parameter();
        1.29 / (1.0 + 0.43 * sg).powf(0.54)
    }
    /// Check if lock-in condition is met.
    ///
    /// Lock-in when `|Ur - 1/St| / (1/St) < 0.3`.
    pub fn is_locked_in(&self) -> bool {
        let ur = self.reduced_velocity();
        let ur_lock = 1.0 / self.strouhal;
        (ur / ur_lock - 1.0).abs() < 0.3
    }
    /// Peak transverse amplitude at lock-in (m).
    pub fn peak_amplitude(&self) -> f64 {
        self.diameter * self.max_amplitude_ratio()
    }
    /// RMS lift force per unit length (N/m) at lock-in.
    pub fn rms_lift_force(&self, fluid_density: f64) -> f64 {
        let fd = 0.5 * fluid_density * self.fluid_velocity.powi(2) * self.diameter;
        fd * self.cl_max / (2.0_f64).sqrt()
    }
    /// Bandwidth of lock-in region in terms of reduced velocity.
    pub fn lock_in_bandwidth(&self) -> f64 {
        0.6 / self.strouhal
    }
}
/// Water entry (slamming) force using von Karman model.
pub struct SlammingForce {
    /// Fluid density ρ (kg/m³).
    pub fluid_density: f64,
    /// Entry velocity V (m/s).
    pub entry_velocity: f64,
    /// Deadrise angle β (degrees).
    pub deadrise_deg: f64,
    /// Maximum beam width B (m).
    pub beam_width: f64,
    /// Section length L (m).
    pub length: f64,
}
impl SlammingForce {
    /// Create a slamming force model for a wedge-shaped section.
    pub fn new(
        fluid_density: f64,
        entry_velocity: f64,
        deadrise_deg: f64,
        beam_width: f64,
        length: f64,
    ) -> Self {
        SlammingForce {
            fluid_density,
            entry_velocity,
            deadrise_deg,
            beam_width,
            length,
        }
    }
    /// Von Karman slamming pressure p_s (Pa): p_s = 0.5 * ρ * Cs * V².
    pub fn slamming_pressure(&self) -> f64 {
        let beta = self.deadrise_deg.to_radians();
        let cs = (PI / (2.0 * beta.tan())).powi(2);
        0.5 * self.fluid_density * cs * self.entry_velocity.powi(2)
    }
    /// Peak slamming force F_s (N) = p_s * B * L.
    pub fn peak_slamming_force(&self) -> f64 {
        self.slamming_pressure() * self.beam_width * self.length
    }
    /// Impact duration estimate Δt = B / (2 * V * tan(β)).
    pub fn impact_duration(&self) -> f64 {
        let beta = self.deadrise_deg.to_radians();
        self.beam_width / (2.0 * self.entry_velocity * beta.tan().max(1e-10))
    }
}
/// Rectangular tank sloshing model using linear wave theory.
///
/// The fundamental sloshing mode is modelled as an equivalent pendulum.
pub struct TankSloshing {
    /// Tank length L (m) in the direction of sloshing.
    pub tank_length: f64,
    /// Liquid depth h (m).
    pub liquid_depth: f64,
    /// Liquid density ρ (kg/m³).
    pub liquid_density: f64,
    /// Gravitational acceleration g (m/s²).
    pub gravity: f64,
    /// Tank width W (m) perpendicular to sloshing.
    pub tank_width: f64,
    /// Structural damping from liquid viscosity / baffles ζ_s.
    pub sloshing_damping: f64,
}
impl TankSloshing {
    /// Create a tank sloshing model.
    pub fn new(tank_length: f64, liquid_depth: f64, tank_width: f64, liquid_density: f64) -> Self {
        Self {
            tank_length,
            liquid_depth,
            liquid_density,
            gravity: 9.81,
            tank_width,
            sloshing_damping: 0.02,
        }
    }
    /// Fundamental sloshing frequency (rad/s): ω₁ = √(π*g/L * tanh(π*h/L)).
    pub fn natural_frequency(&self) -> f64 {
        let l = self.tank_length;
        let h = self.liquid_depth;
        let g = self.gravity;
        (PI * g / l * (PI * h / l).tanh()).sqrt()
    }
    /// Equivalent pendulum length for TLD analogy (m).
    ///
    /// `l_eq = g / ω₁²`.
    pub fn equivalent_pendulum_length(&self) -> f64 {
        let omega1 = self.natural_frequency();
        self.gravity / omega1.powi(2).max(1e-15)
    }
    /// Effective liquid mass participating in first mode (kg).
    ///
    /// `m_eff = (8/π²) * ρ * L * W * h * tanh(π*h/L)`.
    pub fn effective_mass(&self) -> f64 {
        let l = self.tank_length;
        let h = self.liquid_depth;
        let term = (PI * h / l).tanh();
        8.0 / PI.powi(2) * self.liquid_density * l * self.tank_width * h * term
    }
    /// Total liquid mass in tank (kg).
    pub fn total_liquid_mass(&self) -> f64 {
        self.liquid_density * self.tank_length * self.tank_width * self.liquid_depth
    }
    /// Mass ratio: effective/total.
    pub fn mass_ratio(&self) -> f64 {
        self.effective_mass() / self.total_liquid_mass().max(1e-30)
    }
    /// Optimal TLD frequency ratio for maximum damping of structure.
    ///
    /// `f_opt = 1 / (1 + m_eff/m_struct)` — Den Hartog absorber tuning.
    pub fn optimal_frequency_ratio(&self, structural_natural_freq: f64) -> f64 {
        let omega1 = self.natural_frequency();
        omega1 / structural_natural_freq.max(1e-15)
    }
    /// Peak free-surface displacement amplitude for harmonic base excitation (m).
    ///
    /// `η_max = X_0 * ω² / (ω₁² - ω²)` (resonance excluded).
    pub fn free_surface_amplitude(&self, base_displacement: f64, excitation_freq: f64) -> f64 {
        let omega1 = self.natural_frequency();
        let omega = excitation_freq;
        let denom = (omega1.powi(2) - omega.powi(2)).abs().max(1e-10);
        base_displacement * omega.powi(2) / denom
    }
}
/// Vortex-induced vibration (VIV) transverse force model.
pub struct VortexSheddingForce {
    /// Strouhal number St (typically 0.2 for a cylinder).
    pub strouhal: f64,
    /// Fluid velocity U (m/s).
    pub fluid_velocity: f64,
    /// Cylinder diameter D (m).
    pub diameter: f64,
    /// Fluid density ρ (kg/m³).
    pub fluid_density: f64,
    /// Lift coefficient CL for transverse force.
    pub cl: f64,
    /// Reference area A = D * L (m²).
    pub reference_area: f64,
}
impl VortexSheddingForce {
    /// Create a vortex shedding model for a cylinder.
    pub fn new(fluid_velocity: f64, diameter: f64, length: f64, fluid_density: f64) -> Self {
        VortexSheddingForce {
            strouhal: 0.2,
            fluid_velocity,
            diameter,
            fluid_density,
            cl: 0.3,
            reference_area: diameter * length,
        }
    }
    /// Shedding frequency f_s = St * U / D (Hz).
    pub fn shedding_frequency(&self) -> f64 {
        self.strouhal * self.fluid_velocity / self.diameter
    }
    /// Shedding angular frequency ω_s (rad/s).
    pub fn shedding_angular_frequency(&self) -> f64 {
        2.0 * PI * self.shedding_frequency()
    }
    /// Peak transverse (lift) force amplitude F_L = 0.5 * ρ * U² * CL * A (N).
    pub fn peak_transverse_force(&self) -> f64 {
        0.5 * self.fluid_density * self.fluid_velocity.powi(2) * self.cl * self.reference_area
    }
    /// Instantaneous transverse force at time t (sinusoidal model).
    pub fn transverse_force_at(&self, t: f64) -> f64 {
        self.peak_transverse_force() * (self.shedding_angular_frequency() * t).sin()
    }
    /// Lock-in condition: check if natural frequency fn is near shedding frequency.
    pub fn is_locked_in(&self, natural_frequency: f64) -> bool {
        let fs = self.shedding_frequency();
        (natural_frequency / fs - 1.0).abs() < 0.1
    }
}
/// Shape type for added mass calculation.
#[derive(Debug, Clone, PartialEq)]
pub enum BodyShape {
    /// Sphere of radius r.
    Sphere,
    /// Infinite cylinder of radius r (per unit length).
    Cylinder,
    /// Prolate spheroid with semi-axes a > b = c.
    ProlateSpheroid,
    /// General ellipsoid with semi-axes a, b, c.
    Ellipsoid,
}
/// Morison equation: hydrodynamic force on slender cylindrical structures.
///
/// F = ρ * Cm * V * du/dt + 0.5 * ρ * Cd * A * u * |u|
pub struct MorisonForce {
    /// Fluid density ρ (kg/m³).
    pub fluid_density: f64,
    /// Inertia coefficient Cm (typically 2.0 for a cylinder).
    pub cm: f64,
    /// Drag coefficient Cd (typically 1.0–1.2).
    pub cd: f64,
    /// Cylinder diameter D (m).
    pub diameter: f64,
    /// Cylinder length L (m).
    pub length: f64,
}
impl MorisonForce {
    /// Create a Morison force model for a cylinder.
    pub fn new(fluid_density: f64, diameter: f64, length: f64) -> Self {
        MorisonForce {
            fluid_density,
            cm: 2.0,
            cd: 1.0,
            diameter,
            length,
        }
    }
    /// Volume of cylinder V = π/4 * D² * L.
    pub fn volume(&self) -> f64 {
        PI / 4.0 * self.diameter.powi(2) * self.length
    }
    /// Projected area A = D * L.
    pub fn projected_area(&self) -> f64 {
        self.diameter * self.length
    }
    /// Inertia force component F_i = ρ * Cm * V * du_dt (N).
    pub fn inertia_force(&self, fluid_acceleration: f64) -> f64 {
        self.fluid_density * self.cm * self.volume() * fluid_acceleration
    }
    /// Drag force component F_d = 0.5 * ρ * Cd * A * u * |u| (N).
    pub fn drag_force_component(&self, fluid_velocity: f64) -> f64 {
        0.5 * self.fluid_density
            * self.cd
            * self.projected_area()
            * fluid_velocity
            * fluid_velocity.abs()
    }
    /// Total Morison force (N) given fluid velocity and acceleration.
    pub fn total_force(&self, fluid_velocity: f64, fluid_acceleration: f64) -> f64 {
        self.inertia_force(fluid_acceleration) + self.drag_force_component(fluid_velocity)
    }
}
/// Full 6×6 added mass matrix for an arbitrary body in potential flow.
///
/// The matrix couples translational (surge, sway, heave) and rotational
/// (roll, pitch, yaw) degrees of freedom.
pub struct AddedMassMatrix {
    /// 6×6 added mass matrix \[M_ij\] (kg or kg·m²).
    pub matrix: [[f64; 6]; 6],
    /// Fluid density (kg/m³).
    pub fluid_density: f64,
    /// Reference volume (m³).
    pub reference_volume: f64,
}
impl AddedMassMatrix {
    /// Create a zero added mass matrix.
    pub fn new(fluid_density: f64, reference_volume: f64) -> Self {
        Self {
            matrix: [[0.0; 6]; 6],
            fluid_density,
            reference_volume,
        }
    }
    /// Create an added mass matrix for a sphere (diagonal, isotropic).
    pub fn sphere(fluid_density: f64, radius: f64) -> Self {
        let vol = 4.0 / 3.0 * PI * radius.powi(3);
        let ma = 0.5 * fluid_density * vol;
        let mut m = Self::new(fluid_density, vol);
        for i in 0..3 {
            m.matrix[i][i] = ma;
        }
        m
    }
    /// Create added mass matrix for an ellipsoid with semi-axes a, b, c.
    pub fn ellipsoid(fluid_density: f64, a: f64, b: f64, c: f64) -> Self {
        let vol = 4.0 / 3.0 * PI * a * b * c;
        let mut m = Self::new(fluid_density, vol);
        let ma_x = 0.5 * fluid_density * vol;
        let ma_y = 0.5 * fluid_density * vol * (b / a).powi(2);
        let ma_z = 0.5 * fluid_density * vol * (c / a).powi(2);
        m.matrix[0][0] = ma_x;
        m.matrix[1][1] = ma_y;
        m.matrix[2][2] = ma_z;
        m
    }
    /// Multiply added mass matrix by acceleration vector (6-DOF).
    ///
    /// `a_vec` = \[ax, ay, az, αx, αy, αz\].  Returns force/moment vector.
    pub fn force_from_acceleration(&self, a_vec: [f64; 6]) -> [f64; 6] {
        let mut result = [0.0f64; 6];
        for i in 0..6 {
            for j in 0..6 {
                result[i] += self.matrix[i][j] * a_vec[j];
            }
        }
        result
    }
    /// Diagonal elements (translational added masses).
    pub fn translational_diagonal(&self) -> [f64; 3] {
        [self.matrix[0][0], self.matrix[1][1], self.matrix[2][2]]
    }
    /// Trace of the full matrix.
    pub fn trace(&self) -> f64 {
        (0..6).map(|i| self.matrix[i][i]).sum()
    }
}
/// Potential flow added mass tensor for common body shapes.
///
/// The added mass represents the inertia of the surrounding fluid accelerated
/// with the body (potential flow approximation).
pub struct AddedMass {
    /// Fluid density ρ (kg/m³).
    pub fluid_density: f64,
    /// Body shape type.
    pub shape: BodyShape,
    /// Primary dimension (radius for sphere/cylinder, semi-axis a for ellipsoid).
    pub dim_a: f64,
    /// Secondary dimension (length for cylinder, semi-axis b for ellipsoid).
    pub dim_b: f64,
    /// Tertiary dimension (semi-axis c for ellipsoid).
    pub dim_c: f64,
}
impl AddedMass {
    /// Create added mass for a sphere.
    pub fn sphere(fluid_density: f64, radius: f64) -> Self {
        AddedMass {
            fluid_density,
            shape: BodyShape::Sphere,
            dim_a: radius,
            dim_b: radius,
            dim_c: radius,
        }
    }
    /// Create added mass for a cylinder (per unit length).
    pub fn cylinder(fluid_density: f64, radius: f64, length: f64) -> Self {
        AddedMass {
            fluid_density,
            shape: BodyShape::Cylinder,
            dim_a: radius,
            dim_b: length,
            dim_c: radius,
        }
    }
    /// Translational added mass (kg). For a sphere: m_a = (2/3) * π * ρ * r³.
    pub fn translational_added_mass(&self) -> f64 {
        match self.shape {
            BodyShape::Sphere => (2.0 / 3.0) * PI * self.fluid_density * self.dim_a.powi(3),
            BodyShape::Cylinder => PI * self.fluid_density * self.dim_a.powi(2) * self.dim_b,
            BodyShape::ProlateSpheroid | BodyShape::Ellipsoid => {
                let vol = (4.0 / 3.0) * PI * self.dim_a * self.dim_b * self.dim_c;
                0.5 * self.fluid_density * vol
            }
        }
    }
    /// Added mass coefficient Ca = m_added / m_displaced.
    pub fn added_mass_coefficient(&self) -> f64 {
        match self.shape {
            BodyShape::Sphere => 0.5,
            BodyShape::Cylinder => 1.0,
            BodyShape::ProlateSpheroid | BodyShape::Ellipsoid => 0.5,
        }
    }
}
/// ALE mesh node with position and velocity fields.
#[derive(Clone, Debug)]
pub struct AleMeshNode {
    /// Reference (material) position (m).
    pub reference_position: [f64; 3],
    /// Current (deformed) position (m).
    pub current_position: [f64; 3],
    /// Mesh velocity ṽ (m/s) — velocity of the mesh frame.
    pub mesh_velocity: [f64; 3],
    /// Material velocity u (m/s) — velocity of the material.
    pub material_velocity: [f64; 3],
}
impl AleMeshNode {
    /// Create an ALE node at the given reference position.
    pub fn new(reference_position: [f64; 3]) -> Self {
        Self {
            reference_position,
            current_position: reference_position,
            mesh_velocity: [0.0; 3],
            material_velocity: [0.0; 3],
        }
    }
    /// Convective velocity: u - ṽ (relative velocity of material w.r.t. mesh).
    pub fn convective_velocity(&self) -> [f64; 3] {
        [
            self.material_velocity[0] - self.mesh_velocity[0],
            self.material_velocity[1] - self.mesh_velocity[1],
            self.material_velocity[2] - self.mesh_velocity[2],
        ]
    }
    /// Advance mesh position by `dt` seconds using mesh velocity.
    pub fn advance_mesh(&mut self, dt: f64) {
        for i in 0..3 {
            self.current_position[i] += self.mesh_velocity[i] * dt;
        }
    }
    /// Displacement from reference position (m).
    pub fn displacement(&self) -> [f64; 3] {
        [
            self.current_position[0] - self.reference_position[0],
            self.current_position[1] - self.reference_position[1],
            self.current_position[2] - self.reference_position[2],
        ]
    }
}
/// Thermally-driven buoyancy and plume rise model.
pub struct ThermalFluidForce {
    /// Ambient fluid density ρ_0 (kg/m³).
    pub ambient_density: f64,
    /// Ambient temperature T_0 (K).
    pub ambient_temp: f64,
    /// Heat source temperature T_s (K).
    pub source_temp: f64,
    /// Gravitational acceleration g (m/s²).
    pub gravity: f64,
    /// Thermal expansion coefficient β = 1/T_0 for ideal gas (1/K).
    pub beta: f64,
}
impl ThermalFluidForce {
    /// Create a thermal fluid force model.
    pub fn new(ambient_density: f64, ambient_temp: f64, source_temp: f64) -> Self {
        ThermalFluidForce {
            ambient_density,
            ambient_temp,
            source_temp,
            gravity: 9.81,
            beta: 1.0 / ambient_temp,
        }
    }
    /// Buoyancy force per unit volume: f_b = ρ_0 * g * β * ΔT (N/m³).
    pub fn buoyancy_per_volume(&self) -> f64 {
        let dt = self.source_temp - self.ambient_temp;
        self.ambient_density * self.gravity * self.beta * dt
    }
    /// Hot fluid density: ρ_hot = ρ_0 * T_0 / T_s (ideal gas).
    pub fn hot_density(&self) -> f64 {
        self.ambient_density * self.ambient_temp / self.source_temp
    }
    /// Archimedes buoyancy per unit volume: f = (ρ_cold - ρ_hot) * g.
    pub fn archimedes_buoyancy(&self) -> f64 {
        let rho_hot = self.hot_density();
        (self.ambient_density - rho_hot) * self.gravity
    }
    /// Morton-Taylor plume rise Δh (m): Δh = 1.6 * F^(1/4) / u_bar * x^(3/4).
    /// `heat_flux_w` = buoyancy flux F (m⁴/s³), `mean_wind` = u_bar (m/s), `x` = downwind distance (m).
    pub fn plume_rise(&self, heat_flux_w: f64, mean_wind: f64, x: f64) -> f64 {
        if mean_wind < 1e-10 {
            return f64::INFINITY;
        }
        let cp = 1005.0;
        let f = self.gravity * heat_flux_w / (self.ambient_density * cp * self.ambient_temp);
        1.6 * f.powf(0.25) / mean_wind * x.powf(0.75)
    }
    /// Grashof number for natural convection: Gr = g * β * ΔT * L³ / ν².
    pub fn grashof(&self, length: f64, kinematic_viscosity: f64) -> f64 {
        let dt = (self.source_temp - self.ambient_temp).abs();
        self.gravity * self.beta * dt * length.powi(3) / kinematic_viscosity.powi(2)
    }
}
/// Aeroelastic flutter and divergence analysis (Theodorsen-based).
pub struct AeroelasticCoupling {
    /// Air density ρ (kg/m³).
    pub air_density: f64,
    /// Wing span b (m).
    pub span: f64,
    /// Chord c (m).
    pub chord: f64,
    /// Bending stiffness EI (N·m²).
    pub ei: f64,
    /// Torsional stiffness GJ (N·m²/rad).
    pub gj: f64,
    /// Mass per unit span m (kg/m).
    pub mass_per_span: f64,
    /// Moment of inertia per unit span Iα (kg·m²/m).
    pub i_alpha: f64,
}
impl AeroelasticCoupling {
    /// Create an aeroelastic analysis model.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        air_density: f64,
        span: f64,
        chord: f64,
        ei: f64,
        gj: f64,
        mass_per_span: f64,
        i_alpha: f64,
    ) -> Self {
        AeroelasticCoupling {
            air_density,
            span,
            chord,
            ei,
            gj,
            mass_per_span,
            i_alpha,
        }
    }
    /// Natural bending frequency ωh (rad/s) — simplified cantilever.
    pub fn bending_frequency(&self) -> f64 {
        let m_total = self.mass_per_span * self.span;
        (3.0 * self.ei / (m_total * self.span.powi(3))).sqrt()
    }
    /// Natural torsional frequency ωα (rad/s) — simplified.
    pub fn torsional_frequency(&self) -> f64 {
        let i_total = self.i_alpha * self.span;
        (self.gj / (i_total * self.span)).sqrt()
    }
    /// Theodorsen flutter speed estimate (simplified Frazer-Duncan):
    /// U_flutter ≈ ωα * b * sqrt(m / (π * ρ * b²)).
    pub fn flutter_speed(&self) -> f64 {
        let b = self.chord / 2.0;
        let wa = self.torsional_frequency();
        let mu = self.mass_per_span / (PI * self.air_density * b * b);
        wa * b * mu.sqrt()
    }
    /// Divergence speed: U_div = sqrt(2 * GJ / (ρ * c² * span * dCl_dalpha)).
    pub fn divergence_speed(&self, dcl_dalpha: f64) -> f64 {
        let num = 2.0 * self.gj;
        let den = self.air_density * self.chord.powi(2) * self.span * dcl_dalpha;
        (num / den.max(1e-10)).sqrt()
    }
    /// Simple gust response peak load factor n_g = 1 + (ρ * U * a * de) / (2 * w_s).
    /// `u` = airspeed, `a` = lift curve slope, `de` = gust intensity, `w_s` = wing loading.
    pub fn gust_response(
        &self,
        airspeed: f64,
        a: f64,
        gust_intensity: f64,
        wing_loading: f64,
    ) -> f64 {
        1.0 + self.air_density * airspeed * a * gust_intensity / (2.0 * wing_loading)
    }
}
/// Tuned Liquid Column Damper (TLCD) model.
///
/// Liquid oscillates in a U-tube and dissipates energy via head losses.
pub struct TuneableLiquidColumnDamper {
    /// Total liquid column length L (m).
    pub column_length: f64,
    /// Horizontal section length B (m).
    pub horizontal_length: f64,
    /// Tube cross-section area A (m²).
    pub cross_section_area: f64,
    /// Liquid density ρ (kg/m³).
    pub liquid_density: f64,
    /// Head loss coefficient ξ (from orifice).
    pub head_loss_coeff: f64,
    /// Gravitational acceleration g (m/s²).
    pub gravity: f64,
}
impl TuneableLiquidColumnDamper {
    /// Create a TLCD.
    pub fn new(
        column_length: f64,
        horizontal_length: f64,
        cross_section_area: f64,
        liquid_density: f64,
        head_loss_coeff: f64,
    ) -> Self {
        Self {
            column_length,
            horizontal_length,
            cross_section_area,
            liquid_density,
            head_loss_coeff,
            gravity: 9.81,
        }
    }
    /// Natural frequency of liquid column (rad/s): ω = sqrt(2*g/L).
    pub fn natural_frequency(&self) -> f64 {
        (2.0 * self.gravity / self.column_length).sqrt()
    }
    /// Liquid mass in damper (kg).
    pub fn liquid_mass(&self) -> f64 {
        self.liquid_density * self.cross_section_area * self.column_length
    }
    /// Effective damping ratio from head loss coefficient.
    ///
    /// `ζ_tlcd = ξ * A * x_dot_rms / (2 * L * ω)`.
    pub fn damping_ratio_approx(&self, x_dot_rms: f64) -> f64 {
        let omega = self.natural_frequency();
        self.head_loss_coeff * x_dot_rms / (2.0 * self.column_length * omega).max(1e-15)
    }
    /// Frequency tuning ratio β = ω_tlcd / ω_struct.
    pub fn tuning_ratio(&self, structural_natural_freq: f64) -> f64 {
        self.natural_frequency() / structural_natural_freq.max(1e-15)
    }
}
/// Buoyancy force calculation for a partially or fully submerged body.
pub struct BuoyancyForce {
    /// Fluid density ρ_f (kg/m³).
    pub fluid_density: f64,
    /// Gravitational acceleration g (m/s²).
    pub gravity: f64,
    /// Submerged volume V_sub (m³).
    pub submerged_volume: f64,
    /// Total volume V (m³).
    pub total_volume: f64,
    /// Center of buoyancy position \[x, y, z\] (m).
    pub center_of_buoyancy: [f64; 3],
}
impl BuoyancyForce {
    /// Create a new buoyancy calculation.
    pub fn new(fluid_density: f64, submerged_volume: f64, center_of_buoyancy: [f64; 3]) -> Self {
        BuoyancyForce {
            fluid_density,
            gravity: 9.81,
            submerged_volume,
            total_volume: submerged_volume,
            center_of_buoyancy,
        }
    }
    /// Displaced fluid mass (kg): m_d = ρ_f * V_sub.
    pub fn displaced_mass(&self) -> f64 {
        self.fluid_density * self.submerged_volume
    }
    /// Buoyancy force magnitude (N): F_b = ρ_f * g * V_sub.
    pub fn buoyancy_magnitude(&self) -> f64 {
        self.fluid_density * self.gravity * self.submerged_volume
    }
    /// Buoyancy force vector (acts upward in global Y): \[0, F_b, 0\].
    pub fn buoyancy_vector(&self) -> [f64; 3] {
        [0.0, self.buoyancy_magnitude(), 0.0]
    }
    /// Fraction of body submerged (0–1).
    pub fn submersion_fraction(&self) -> f64 {
        if self.total_volume < 1e-30 {
            return 0.0;
        }
        (self.submerged_volume / self.total_volume).min(1.0)
    }
    /// Net vertical force on body: F_b - W. Positive = upward net.
    pub fn net_vertical_force(&self, body_mass: f64) -> f64 {
        self.buoyancy_magnitude() - body_mass * self.gravity
    }
}
/// Monolithic FSI coupling — assembles a combined fluid-structure system.
///
/// Uses a simplified 1-DOF model where the structure and fluid are solved
/// together with a single effective stiffness matrix.
pub struct MonolithicFsiCoupling {
    /// Added mass from the fluid (kg).
    pub added_mass: f64,
    /// Fluid damping coefficient (N·s/m).
    pub fluid_damping: f64,
    /// Fluid stiffness coefficient (N/m).
    pub fluid_stiffness: f64,
    /// Structural mass (kg).
    pub structural_mass: f64,
    /// Structural stiffness (N/m).
    pub structural_stiffness: f64,
    /// Structural damping (N·s/m).
    pub structural_damping: f64,
    /// External forcing amplitude (N).
    pub forcing_amplitude: f64,
    /// Forcing frequency (rad/s).
    pub forcing_frequency: f64,
}
impl MonolithicFsiCoupling {
    /// Create a monolithic FSI coupler.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        structural_mass: f64,
        structural_stiffness: f64,
        structural_damping: f64,
        added_mass: f64,
        fluid_damping: f64,
        fluid_stiffness: f64,
    ) -> Self {
        Self {
            added_mass,
            fluid_damping,
            fluid_stiffness,
            structural_mass,
            structural_stiffness,
            structural_damping,
            forcing_amplitude: 0.0,
            forcing_frequency: 0.0,
        }
    }
    /// Effective total mass = structural + added.
    pub fn total_mass(&self) -> f64 {
        self.structural_mass + self.added_mass
    }
    /// Effective total stiffness.
    pub fn total_stiffness(&self) -> f64 {
        self.structural_stiffness + self.fluid_stiffness
    }
    /// Effective total damping.
    pub fn total_damping(&self) -> f64 {
        self.structural_damping + self.fluid_damping
    }
    /// Combined natural frequency (rad/s).
    pub fn natural_frequency(&self) -> f64 {
        (self.total_stiffness() / self.total_mass().max(1e-30)).sqrt()
    }
    /// Combined damping ratio ζ.
    pub fn damping_ratio(&self) -> f64 {
        let omega_n = self.natural_frequency();
        self.total_damping() / (2.0 * self.total_mass() * omega_n).max(1e-30)
    }
    /// Steady-state displacement amplitude under harmonic forcing (m).
    pub fn steady_state_amplitude(&self) -> f64 {
        let omega = self.forcing_frequency;
        let k = self.total_stiffness();
        let m = self.total_mass();
        let c = self.total_damping();
        let f0 = self.forcing_amplitude;
        let denom = ((k - m * omega * omega).powi(2) + (c * omega).powi(2)).sqrt();
        f0 / denom.max(1e-30)
    }
}
/// Radiation damping matrix for a floating body.
///
/// Arises from wave generation due to body oscillations (frequency-dependent).
pub struct RadiationDampingMatrix {
    /// 6×6 damping matrix \[B_ij\] at a given frequency (N·s/m or N·m·s/rad).
    pub matrix: [[f64; 6]; 6],
    /// Angular frequency ω at which the matrix is defined (rad/s).
    pub omega: f64,
}
impl RadiationDampingMatrix {
    /// Create a zero radiation damping matrix at frequency `omega`.
    pub fn new(omega: f64) -> Self {
        Self {
            matrix: [[0.0; 6]; 6],
            omega,
        }
    }
    /// Set diagonal radiation damping (simplified diagonal model).
    pub fn set_diagonal(&mut self, b: [f64; 6]) {
        for i in 0..6 {
            self.matrix[i][i] = b[i];
        }
    }
    /// Damping force from velocity vector.
    pub fn force_from_velocity(&self, vel: [f64; 6]) -> [f64; 6] {
        let mut result = [0.0f64; 6];
        for i in 0..6 {
            for j in 0..6 {
                result[i] += self.matrix[i][j] * vel[j];
            }
        }
        result
    }
    /// Power dissipated: P = vᵀ B v.
    pub fn dissipated_power(&self, vel: [f64; 6]) -> f64 {
        let fv = self.force_from_velocity(vel);
        vel.iter().zip(fv.iter()).map(|(v, f)| v * f).sum()
    }
}
/// Partitioned FSI coupling scheme (weakly or strongly coupled).
///
/// In a partitioned approach the fluid and structural solvers are called
/// separately and exchange boundary-condition information at the interface.
pub struct PartitionedFsiCoupling {
    /// Maximum number of sub-iterations per time step (strongly coupled).
    pub max_iterations: u32,
    /// Convergence tolerance on interface displacement (metres).
    pub displacement_tol: f64,
    /// Under-relaxation factor ω ∈ (0, 1].
    pub relaxation: f64,
    /// Current interface displacement (metres).
    pub interface_displacement: f64,
    /// Current interface velocity (m/s).
    pub interface_velocity: f64,
    /// Fluid force on interface (N).
    pub fluid_force: f64,
    /// Structural stiffness at interface (N/m).
    pub structural_stiffness: f64,
    /// Structural damping at interface (N·s/m).
    pub structural_damping: f64,
    /// Structural mass at interface (kg).
    pub structural_mass: f64,
}
impl PartitionedFsiCoupling {
    /// Create a partitioned FSI coupler.
    #[allow(clippy::too_many_arguments)]
    pub fn new(structural_mass: f64, structural_stiffness: f64, structural_damping: f64) -> Self {
        Self {
            max_iterations: 20,
            displacement_tol: 1e-6,
            relaxation: 0.5,
            interface_displacement: 0.0,
            interface_velocity: 0.0,
            fluid_force: 0.0,
            structural_stiffness,
            structural_damping,
            structural_mass,
        }
    }
    /// Perform one sub-iteration using Aitken acceleration.
    ///
    /// `fluid_force` is the new fluid force predicted for this sub-iteration.
    /// Returns the displacement residual.
    pub fn sub_iterate(&mut self, fluid_force: f64, dt: f64) -> f64 {
        self.fluid_force = fluid_force;
        let beta = 0.25_f64;
        let gamma = 0.5_f64;
        let k = self.structural_stiffness;
        let c = self.structural_damping;
        let m = self.structural_mass;
        let k_eff = k + gamma / (beta * dt) * c + 1.0 / (beta * dt * dt) * m;
        let d_new = (fluid_force) / k_eff.max(1e-15);
        let residual = (d_new - self.interface_displacement).abs();
        self.interface_displacement =
            self.relaxation * d_new + (1.0 - self.relaxation) * self.interface_displacement;
        self.interface_velocity = (d_new) / dt;
        residual
    }
    /// Check convergence.
    pub fn is_converged(&self, residual: f64) -> bool {
        residual < self.displacement_tol
    }
    /// Natural frequency of the structural subsystem (rad/s).
    pub fn natural_frequency(&self) -> f64 {
        (self.structural_stiffness / self.structural_mass.max(1e-30)).sqrt()
    }
    /// Critical damping ratio ζ.
    pub fn damping_ratio(&self) -> f64 {
        let omega_n = self.natural_frequency();
        self.structural_damping / (2.0 * self.structural_mass * omega_n).max(1e-30)
    }
}
