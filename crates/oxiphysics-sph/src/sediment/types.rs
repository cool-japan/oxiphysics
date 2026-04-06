//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::needless_range_loop)]
use super::functions::*;
use super::types_advanced::ShieldsParam;

/// Exner equation bed evolution with bed-load divergence.
///
/// ∂η/∂t = -1/(1-n) ∂q_b/∂x
pub struct ExnerEquation {
    /// Bed elevation grid η (m).
    pub eta: Vec<f64>,
    /// Grid spacing Δx (m).
    pub dx: f64,
    /// Bed porosity n (−).
    pub porosity: f64,
}
impl ExnerEquation {
    /// Construct an Exner equation model.
    pub fn new(eta: Vec<f64>, dx: f64, porosity: f64) -> Self {
        ExnerEquation { eta, dx, porosity }
    }
    /// Update bed elevation one time step from bed-load flux array q_b (m²/s).
    ///
    /// Uses upwind finite difference for ∂q_b/∂x.
    pub fn update(&mut self, q_b: &[f64], dt: f64) {
        let n = self.eta.len();
        if n < 2 || q_b.len() != n {
            return;
        }
        let scale = 1.0 / (1.0 - self.porosity);
        for i in 1..n {
            let dq = (q_b[i] - q_b[i - 1]) / self.dx;
            self.eta[i] -= dt * scale * dq;
        }
    }
    /// Volume of sediment eroded (positive) or deposited (negative) per unit length.
    pub fn volume_change(&self, eta_0: &[f64]) -> f64 {
        self.eta
            .iter()
            .zip(eta_0.iter())
            .map(|(e, e0)| (e - e0) * self.dx)
            .sum()
    }
    /// Maximum scour depth.
    pub fn max_scour(&self, eta_0: &[f64]) -> f64 {
        self.eta
            .iter()
            .zip(eta_0.iter())
            .map(|(e, e0)| e0 - e)
            .fold(0.0f64, f64::max)
    }
}
/// Meyer–Peter–Müller bed-load transport model.
pub struct BedLoadTransport {
    /// Shields parameter object.
    pub shields: ShieldsParam,
    /// Calibration coefficient (default 8.0 in MPM).
    pub mpm_coeff: f64,
}
impl BedLoadTransport {
    /// Construct with default MPM coefficient 8.0.
    pub fn new(grain_diameter: f64, rho_s: f64, rho_f: f64, theta_cr: f64) -> Self {
        BedLoadTransport {
            shields: ShieldsParam::new(grain_diameter, rho_s, rho_f, theta_cr),
            mpm_coeff: 8.0,
        }
    }
    /// Bed-load volumetric transport rate per unit width q_b (m²/s).
    pub fn transport_rate(&self, tau_b: f64) -> f64 {
        let theta = self.shields.theta(tau_b);
        if theta <= self.shields.theta_cr {
            return 0.0;
        }
        let excess = (theta - self.shields.theta_cr).powf(1.5);
        let d = self.shields.grain_diameter;
        let s = self.shields.rho_s / self.shields.rho_f;
        self.mpm_coeff * excess * ((s - 1.0) * G * d * d * d).sqrt()
    }
    /// Transport vector aligned with the bed shear stress direction.
    pub fn transport_vector(&self, tau_b_vec: [f64; 3]) -> [f64; 3] {
        let mag = len3(tau_b_vec);
        if mag < 1e-300 {
            return [0.0; 3];
        }
        let q = self.transport_rate(mag);
        scale3(
            [tau_b_vec[0] / mag, tau_b_vec[1] / mag, tau_b_vec[2] / mag],
            q,
        )
    }
}
/// Partheniades-type erosion model.
pub struct ErosionModel {
    /// Critical bed shear stress τ_cr (Pa).
    pub tau_cr: f64,
    /// Erosion rate coefficient M (kg/m²/s).
    pub erosion_coeff: f64,
    /// Exponent in excess-shear formula.
    pub exponent: f64,
}
impl ErosionModel {
    /// Create a new erosion model.
    pub fn new(tau_cr: f64, erosion_coeff: f64, exponent: f64) -> Self {
        ErosionModel {
            tau_cr,
            erosion_coeff,
            exponent,
        }
    }
    /// Erosion rate E (kg/m²/s) for bed shear stress τ_b.
    pub fn erosion_rate(&self, tau_b: f64) -> f64 {
        if tau_b <= self.tau_cr {
            return 0.0;
        }
        self.erosion_coeff * ((tau_b / self.tau_cr - 1.0).max(0.0)).powf(self.exponent)
    }
    /// Sediment flux vector (kg/m/s) aligned with τ_b_vec.
    pub fn sediment_flux(&self, tau_b_vec: [f64; 3]) -> [f64; 3] {
        let mag = len3(tau_b_vec);
        let rate = self.erosion_rate(mag);
        if mag < 1e-300 {
            return [0.0; 3];
        }
        scale3(
            [tau_b_vec[0] / mag, tau_b_vec[1] / mag, tau_b_vec[2] / mag],
            rate,
        )
    }
}
/// Scour depth estimator for bridge piers and abutments.
///
/// Implements HEC-18 / CSU pier scour equation and
/// the Froehlich abutment scour equation.
pub struct ScourModel {
    /// Pier width b (m).
    pub pier_width: f64,
    /// Water depth upstream y_0 (m).
    pub upstream_depth: f64,
    /// Pier shape factor K1 (1.1 for round, 1.3 for square).
    pub shape_factor: f64,
    /// Attack angle factor K2.
    pub angle_factor: f64,
    /// Bed condition factor K3.
    pub bed_factor: f64,
    /// Armouring factor K4.
    pub armouring_factor: f64,
}
impl ScourModel {
    /// Construct a pier scour model (CSU / HEC-18).
    pub fn new_pier(pier_width: f64, upstream_depth: f64) -> Self {
        ScourModel {
            pier_width,
            upstream_depth,
            shape_factor: 1.1,
            angle_factor: 1.0,
            bed_factor: 1.1,
            armouring_factor: 1.0,
        }
    }
    /// CSU pier scour depth y_s (m).
    ///
    /// y_s = 2.0 K1 K2 K3 K4 (y_0/b)^0.35 Fr^0.43 b
    pub fn csu_pier_scour(&self, froude_number: f64) -> f64 {
        let ratio = (self.upstream_depth / self.pier_width).powf(0.35);
        2.0 * self.shape_factor
            * self.angle_factor
            * self.bed_factor
            * self.armouring_factor
            * ratio
            * froude_number.powf(0.43)
            * self.pier_width
    }
    /// Melville-Coleman pier scour formula (m).
    ///
    /// d_s = K_I K_d K_y K_s K_a y_0 where K factors depend on flow and geometry.
    pub fn melville_pier_scour(&self, froude_number: f64, d50: f64) -> f64 {
        let ky = if self.upstream_depth / self.pier_width < 0.7 {
            (self.upstream_depth * self.pier_width).sqrt()
        } else {
            2.4 * self.pier_width
        };
        let u_star_c = G.sqrt() * (self.upstream_depth * self.bed_factor * d50 * 0.001).sqrt();
        let u_c = 5.75 * u_star_c * (self.upstream_depth / d50).log10();
        let fr_local = froude_number * u_c.max(1e-300);
        let _fr_local_unused = fr_local;
        self.shape_factor * ky
    }
    /// Froehlich abutment scour depth y_s (m).
    ///
    /// y_s = 2.27 K1 K2 (a'/y_0)^0.43 Fr^0.61 y_0
    pub fn froehlich_abutment_scour(&self, flow_length: f64, froude_number: f64) -> f64 {
        let ratio = (flow_length / self.upstream_depth.max(1e-300)).powf(0.43);
        2.27 * self.shape_factor
            * self.angle_factor
            * ratio
            * froude_number.powf(0.61)
            * self.upstream_depth
    }
    /// Time scale of scour development (hours, empirical).
    pub fn scour_time_scale(&self) -> f64 {
        self.upstream_depth / 1e-5 / 3600.0
    }
}
/// Van Rijn (1984) bed-load transport formula.
///
/// q_b = 0.053 d T*^{2.1} D*^{-0.3} sqrt((s-1)g d³)
pub struct VanRijnBedLoad {
    /// Grain diameter d₅₀ (m).
    pub d50: f64,
    /// Sediment specific gravity s = ρ_s / ρ_f.
    pub specific_gravity: f64,
    /// Kinematic viscosity ν (m²/s).
    pub nu: f64,
}
impl VanRijnBedLoad {
    /// Construct a van Rijn bed-load model.
    pub fn new(d50: f64, rho_s: f64, rho_f: f64, nu: f64) -> Self {
        VanRijnBedLoad {
            d50,
            specific_gravity: rho_s / rho_f,
            nu,
        }
    }
    /// Dimensionless grain size D* = d \[(s-1)g/ν²\]^{1/3}.
    pub fn dimensionless_grain_size(&self) -> f64 {
        let s = self.specific_gravity;
        self.d50 * ((s - 1.0) * G / (self.nu * self.nu)).powf(1.0 / 3.0)
    }
    /// Transport stage T* = (u*² - u*_cr²) / u*_cr² (excess shear).
    pub fn transport_stage(&self, u_star: f64, u_star_cr: f64) -> f64 {
        if u_star_cr < 1e-300 {
            return 0.0;
        }
        ((u_star * u_star - u_star_cr * u_star_cr) / (u_star_cr * u_star_cr)).max(0.0)
    }
    /// Van Rijn bed-load transport rate q_b (m²/s).
    pub fn transport_rate(&self, u_star: f64, u_star_cr: f64) -> f64 {
        let t_star = self.transport_stage(u_star, u_star_cr);
        if t_star < 1e-300 {
            return 0.0;
        }
        let d_star = self.dimensionless_grain_size();
        let s = self.specific_gravity;
        let d = self.d50;
        if d_star < 1e-300 || d < 1e-300 {
            return 0.0;
        }
        0.053 * d * t_star.powf(2.1) / d_star.powf(0.3) * ((s - 1.0) * G * d * d * d).sqrt()
    }
}
/// Rouse sediment concentration profile.
pub struct RouseProfile {
    /// Total water depth H (m).
    pub depth: f64,
    /// Reference height a (m) above bed.
    pub ref_height: f64,
    /// Reference concentration C_a (−).
    pub ref_concentration: f64,
    /// Rouse number Ro = w_s / (κ u*).
    pub rouse_number: f64,
    /// von Kármán constant κ (≈ 0.41).
    pub kappa: f64,
    /// Shear velocity u* (m/s).
    pub shear_velocity: f64,
}
impl RouseProfile {
    /// Construct a Rouse profile.
    pub fn new(
        depth: f64,
        ref_height: f64,
        ref_concentration: f64,
        settling_velocity: f64,
        shear_velocity: f64,
    ) -> Self {
        let kappa = 0.41;
        let rouse_number = if shear_velocity.abs() > 1e-300 {
            settling_velocity / (kappa * shear_velocity)
        } else {
            0.0
        };
        RouseProfile {
            depth,
            ref_height,
            ref_concentration,
            rouse_number,
            kappa,
            shear_velocity,
        }
    }
    /// Evaluate concentration at elevation `z` above bed.
    pub fn concentration_at(&self, z: f64) -> f64 {
        rouse_concentration(
            z,
            self.depth,
            self.ref_height,
            self.ref_concentration,
            self.rouse_number,
        )
    }
    /// Depth-averaged suspended sediment concentration.
    ///
    /// Uses simple trapezoidal rule with `n_points` quadrature points.
    pub fn depth_averaged_concentration(&self, n_points: usize) -> f64 {
        let a = self.ref_height;
        let b = self.depth * 0.99;
        if a >= b || n_points < 2 {
            return 0.0;
        }
        let dz = (b - a) / (n_points - 1) as f64;
        let mut sum = 0.0;
        for i in 0..n_points {
            let z = a + i as f64 * dz;
            let c = self.concentration_at(z);
            let weight = if i == 0 || i == n_points - 1 {
                0.5
            } else {
                1.0
            };
            sum += weight * c * dz;
        }
        sum / (self.depth - self.ref_height)
    }
}
/// Engelund-Hansen total sediment transport formula.
///
/// q_t = 0.05 u^5 / (g^0.5 d^{3/2} (s-1)^2 C^3)
///
/// where u is depth-averaged velocity, C is Chézy coefficient.
pub struct EngelundHansen {
    /// Grain diameter d (m).
    pub grain_diameter: f64,
    /// Sediment specific gravity s = ρ_s / ρ_f.
    pub specific_gravity: f64,
    /// Chézy friction coefficient C (m^{1/2}/s).
    pub chezy_coeff: f64,
}
impl EngelundHansen {
    /// Construct an Engelund-Hansen transport model.
    pub fn new(grain_diameter: f64, rho_s: f64, rho_f: f64, chezy_coeff: f64) -> Self {
        EngelundHansen {
            grain_diameter,
            specific_gravity: rho_s / rho_f,
            chezy_coeff,
        }
    }
    /// Total sediment transport rate q_t (m²/s) for depth-averaged velocity u (m/s).
    pub fn total_transport(&self, u: f64) -> f64 {
        let d = self.grain_diameter;
        let s = self.specific_gravity;
        let c = self.chezy_coeff;
        let g = G;
        if d < 1e-300 || (s - 1.0).abs() < 1e-300 || c < 1e-300 {
            return 0.0;
        }
        0.05 * u.powi(5) / (g.sqrt() * d.powf(1.5) * (s - 1.0).powi(2) * c.powi(3))
    }
    /// Shields parameter from Chézy velocity: θ = u² / (C² (s-1) d).
    pub fn shields_from_velocity(&self, u: f64) -> f64 {
        let c = self.chezy_coeff;
        let s = self.specific_gravity;
        let d = self.grain_diameter;
        if c < 1e-300 || (s - 1.0).abs() < 1e-300 || d < 1e-300 {
            return 0.0;
        }
        u * u / (c * c * (s - 1.0) * d)
    }
    /// Bed shear stress from Chézy velocity (Pa) with fluid density.
    pub fn bed_shear_stress(&self, u: f64, rho_f: f64, depth: f64) -> f64 {
        let c = self.chezy_coeff;
        if c < 1e-300 {
            return 0.0;
        }
        rho_f * G * depth * u * u / (c * c)
    }
}
/// Sediment transport mode classification by Rouse number.
#[derive(Clone, Debug, PartialEq)]
pub enum SedimentTransportMode {
    /// Ro > 2.5: bed-load (rolling/sliding).
    BedLoad,
    /// 1.2 < Ro ≤ 2.5: saltation load (bouncing).
    SaltationLoad,
    /// 0.8 < Ro ≤ 1.2: suspended load.
    SuspendedLoad,
    /// Ro ≤ 0.8: wash load (fully suspended).
    WashLoad,
}
/// Multi-formula settling velocity calculator.
pub struct SettlingVelocity {
    /// Grain diameter d (m).
    pub diameter: f64,
    /// Sediment density ρ_s (kg/m³).
    pub rho_s: f64,
    /// Fluid density ρ_f (kg/m³).
    pub rho_f: f64,
    /// Kinematic viscosity ν (m²/s).
    pub nu: f64,
}
impl SettlingVelocity {
    /// Construct a settling velocity calculator.
    pub fn new(diameter: f64, rho_s: f64, rho_f: f64, nu: f64) -> Self {
        SettlingVelocity {
            diameter,
            rho_s,
            rho_f,
            nu,
        }
    }
    /// Stokes settling velocity (valid for small Reynolds numbers).
    pub fn stokes(&self) -> f64 {
        let mu = self.rho_f * self.nu;
        settling_velocity_stokes(self.diameter, self.rho_s, self.rho_f, mu)
    }
    /// Rubey (1933) settling velocity formula.
    ///
    /// w_s = √((2/3 + 36 ν²/(g d³ (s-1)))^{1/2} − √(36 ν²/(g d³ (s-1)))) * √((s-1)*g*d)
    pub fn rubey(&self) -> f64 {
        let s = self.rho_s / self.rho_f;
        let d = self.diameter;
        let nu = self.nu;
        if (s - 1.0).abs() < 1e-300 || d < 1e-300 {
            return 0.0;
        }
        let base = (s - 1.0) * G * d;
        let xi = 36.0 * nu * nu / (base * d);
        let inner = (2.0 / 3.0 + xi).sqrt() - xi.sqrt();
        if inner < 0.0 {
            return 0.0;
        }
        inner * base.sqrt()
    }
    /// Drag-corrected settling velocity for non-spherical particles.
    ///
    /// Uses Corey shape factor `csf` ∈ (0, 1] to correct Stokes velocity.
    pub fn non_spherical(&self, csf: f64) -> f64 {
        let vs_stokes = self.stokes();
        vs_stokes * csf.clamp(0.0, 1.0).sqrt()
    }
    /// Particle Reynolds number Re_p = w_s * d / ν for Stokes settling.
    pub fn reynolds_number(&self) -> f64 {
        let vs = self.stokes();
        vs * self.diameter / self.nu.max(1e-300)
    }
}
/// Multi-regime settling velocity: Stokes (Re_p < 0.5), transitional, Newton (Re_p > 500).
///
/// Uses the Dietrich (1982) / Ferguson & Church (2004) smooth transition formula.
pub struct MultiRegimeSettling {
    /// Grain diameter d (m).
    pub diameter: f64,
    /// Sediment density ρ_s (kg/m³).
    pub rho_s: f64,
    /// Fluid density ρ_f (kg/m³).
    pub rho_f: f64,
    /// Kinematic viscosity ν (m²/s).
    pub nu: f64,
}
impl MultiRegimeSettling {
    /// Construct a multi-regime settling velocity calculator.
    pub fn new(diameter: f64, rho_s: f64, rho_f: f64, nu: f64) -> Self {
        MultiRegimeSettling {
            diameter,
            rho_s,
            rho_f,
            nu,
        }
    }
    /// Ferguson & Church (2004) settling velocity formula.
    ///
    /// w_s = R g d² / (C_1 ν + (0.75 C_2 R g d³)^{0.5})
    ///
    /// C_1 = 18, C_2 = 1.0 (for natural grains).
    pub fn ferguson_church(&self) -> f64 {
        let r = (self.rho_s - self.rho_f) / self.rho_f.max(1e-300);
        let d = self.diameter;
        let nu = self.nu;
        if d < 1e-300 || nu < 1e-300 {
            return 0.0;
        }
        let c1 = 18.0_f64;
        let c2 = 1.0_f64;
        let numer = r * G * d * d;
        let denom = c1 * nu + (0.75 * c2 * r * G * d * d * d).sqrt();
        if denom < 1e-300 {
            return 0.0;
        }
        numer / denom
    }
    /// Stokes settling velocity (valid Re_p < 0.5).
    pub fn stokes(&self) -> f64 {
        let mu = self.rho_f * self.nu;
        settling_velocity_stokes(self.diameter, self.rho_s, self.rho_f, mu)
    }
    /// Newton-regime terminal velocity: w_s = sqrt(4/3 * R g d / C_D).
    ///
    /// C_D ≈ 0.44 in the Newton regime.
    pub fn newton_regime(&self) -> f64 {
        let r = (self.rho_s - self.rho_f) / self.rho_f.max(1e-300);
        let d = self.diameter;
        let cd = 0.44_f64;
        if d < 1e-300 {
            return 0.0;
        }
        (4.0 / 3.0 * r * G * d / cd).sqrt()
    }
    /// Particle Reynolds number Re_p = w_s d / ν.
    pub fn particle_reynolds_number(&self) -> f64 {
        let ws = self.ferguson_church();
        ws * self.diameter / self.nu.max(1e-300)
    }
    /// Select the appropriate regime and return the settling velocity.
    pub fn settling_velocity(&self) -> f64 {
        let re = self.particle_reynolds_number();
        if re < 0.5 {
            self.stokes()
        } else {
            self.ferguson_church()
        }
    }
}
/// Krone-Partheniades cohesive sediment model.
///
/// Handles deposition (Krone, 1962) and erosion (Partheniades, 1965) of
/// fine-grained cohesive sediment (clay, fine silt).
pub struct CohesiveSediment {
    /// Critical shear stress for deposition τ_d (Pa).
    pub tau_deposition: f64,
    /// Critical shear stress for erosion τ_e (Pa).
    pub tau_erosion: f64,
    /// Erosion rate coefficient M (kg/m²/s).
    pub erosion_coeff: f64,
    /// Floc settling velocity w_s (m/s).
    pub settling_velocity: f64,
    /// Bed surface concentration (kg/m²).
    pub bed_mass: f64,
}
impl CohesiveSediment {
    /// Construct a cohesive sediment model.
    pub fn new(
        tau_deposition: f64,
        tau_erosion: f64,
        erosion_coeff: f64,
        settling_velocity: f64,
    ) -> Self {
        CohesiveSediment {
            tau_deposition,
            tau_erosion,
            erosion_coeff,
            settling_velocity,
            bed_mass: 0.0,
        }
    }
    /// Krone deposition flux D (kg/m²/s).
    ///
    /// D = w_s C (1 − τ_b/τ_d) for τ_b < τ_d, else 0.
    pub fn deposition_flux(&self, tau_b: f64, c_near_bed: f64) -> f64 {
        if tau_b >= self.tau_deposition {
            return 0.0;
        }
        let prob = 1.0 - tau_b / self.tau_deposition.max(1e-300);
        self.settling_velocity * c_near_bed * prob.max(0.0)
    }
    /// Partheniades erosion flux E (kg/m²/s).
    ///
    /// E = M (τ_b/τ_e − 1) for τ_b > τ_e, else 0.
    pub fn erosion_flux(&self, tau_b: f64) -> f64 {
        if tau_b <= self.tau_erosion {
            return 0.0;
        }
        let excess = tau_b / self.tau_erosion.max(1e-300) - 1.0;
        self.erosion_coeff * excess.max(0.0)
    }
    /// Net flux (positive = net deposition, negative = net erosion).
    pub fn net_flux(&self, tau_b: f64, c_near_bed: f64) -> f64 {
        self.deposition_flux(tau_b, c_near_bed) - self.erosion_flux(tau_b)
    }
    /// Update bed mass over time step `dt` (s).
    pub fn update_bed(&mut self, tau_b: f64, c_near_bed: f64, dt: f64) {
        let net = self.net_flux(tau_b, c_near_bed);
        self.bed_mass = (self.bed_mass + net * dt).max(0.0);
    }
    /// Sediment concentration change dC/dt in water column.
    ///
    /// `depth` is flow depth (m). Positive net deposition removes from water column.
    pub fn dc_dt(&self, tau_b: f64, c_near_bed: f64, depth: f64) -> f64 {
        let net = self.net_flux(tau_b, c_near_bed);
        if depth < 1e-300 {
            return 0.0;
        }
        -net / depth
    }
}
/// 1-D turbulent diffusion of suspended sediment with Fickian model.
///
/// ∂C/∂t = ∂/∂z (ε_s ∂C/∂z) − w_s ∂C/∂z
///
/// Represents the vertical concentration profile evolution.
pub struct VerticalDiffusion {
    /// Concentration profile C(z_i) (volumetric −).
    pub profile: Vec<f64>,
    /// Vertical grid spacing Δz (m).
    pub dz: f64,
    /// Turbulent diffusivity profile ε_s(z_i) (m²/s).
    pub diffusivity: Vec<f64>,
    /// Settling velocity w_s (m/s, positive = downward).
    pub settling_velocity: f64,
}
impl VerticalDiffusion {
    /// Construct a vertical diffusion model with uniform diffusivity.
    pub fn new(n_cells: usize, dz: f64, eps: f64, ws: f64) -> Self {
        VerticalDiffusion {
            profile: vec![0.0; n_cells],
            dz,
            diffusivity: vec![eps; n_cells],
            settling_velocity: ws,
        }
    }
    /// Advance one time step using explicit finite differences.
    pub fn advance(&mut self, dt: f64) {
        let n = self.profile.len();
        if n < 3 {
            return;
        }
        let mut new_c = self.profile.clone();
        let dz = self.dz;
        let ws = self.settling_velocity;
        for i in 1..n - 1 {
            let eps_p = (self.diffusivity[i] + self.diffusivity[i + 1]) * 0.5;
            let eps_m = (self.diffusivity[i] + self.diffusivity[i - 1]) * 0.5;
            let diff = (eps_p * (self.profile[i + 1] - self.profile[i])
                - eps_m * (self.profile[i] - self.profile[i - 1]))
                / (dz * dz);
            let settle = ws * (self.profile[i + 1] - self.profile[i - 1]) / (2.0 * dz);
            new_c[i] = (self.profile[i] + dt * diff - dt * settle).max(0.0);
        }
        self.profile = new_c;
    }
    /// Depth-averaged concentration.
    pub fn mean_concentration(&self) -> f64 {
        if self.profile.is_empty() {
            return 0.0;
        }
        self.profile.iter().sum::<f64>() / self.profile.len() as f64
    }
    /// Stability criterion: Δt < Δz² / (2 ε_s) and Δt < Δz / w_s.
    pub fn max_dt(&self) -> f64 {
        let eps_max = self
            .diffusivity
            .iter()
            .cloned()
            .fold(0.0_f64, f64::max)
            .max(1e-300);
        let diff_dt = self.dz * self.dz / (2.0 * eps_max);
        let settle_dt = self.dz / self.settling_velocity.max(1e-300);
        diff_dt.min(settle_dt)
    }
}
/// 1D bed profile for erosion/deposition modelling.
pub struct SedimentBed {
    /// Horizontal positions of bed nodes (m).
    pub x: Vec<f64>,
    /// Bed elevation at each node (m).
    pub z_bed: Vec<f64>,
    /// Horizontal spacing (m).
    pub dx: f64,
}
impl SedimentBed {
    /// Create a flat bed with `n` nodes and spacing `dx`.
    pub fn new(n: usize, dx: f64) -> Self {
        Self {
            x: (0..n).map(|i| i as f64 * dx).collect(),
            z_bed: vec![0.0; n],
            dx,
        }
    }
    /// Bed elevation at node index `i`.
    pub fn depth_at(&self, i: usize) -> f64 {
        if i < self.z_bed.len() {
            self.z_bed[i]
        } else {
            0.0
        }
    }
    /// Erode the bed: z_bed\[i\] -= rates\[i\] * dt.
    pub fn erode(&mut self, rates: &[f64], dt: f64) {
        let n = self.z_bed.len().min(rates.len());
        for i in 0..n {
            self.z_bed[i] -= rates[i] * dt;
        }
    }
    /// Deposit on the bed: z_bed\[i\] += rates\[i\] * dt.
    pub fn deposit(&mut self, rates: &[f64], dt: f64) {
        let n = self.z_bed.len().min(rates.len());
        for i in 0..n {
            self.z_bed[i] += rates[i] * dt;
        }
    }
    /// Total bed volume: sum of z_bed\[i\] * dx.
    pub fn total_volume(&self) -> f64 {
        self.z_bed.iter().sum::<f64>() * self.dx
    }
}
/// Morphodynamic acceleration factor for SPH bed evolution.
///
/// Allows morphological time scales to be accelerated relative to
/// hydrodynamic time scales (commonly used in coastal engineering).
pub struct MorphodynamicFeedback {
    /// Morphological acceleration factor f_morph (−).
    pub f_morph: f64,
    /// Bed roughness height k_s (m).
    pub roughness: f64,
    /// Active layer thickness d_a (m).
    pub active_layer: f64,
}
impl MorphodynamicFeedback {
    /// Construct a morphodynamic feedback model.
    pub fn new(f_morph: f64, roughness: f64, active_layer: f64) -> Self {
        MorphodynamicFeedback {
            f_morph,
            roughness,
            active_layer,
        }
    }
    /// Effective bed-load flux for morphodynamic update: q_eff = f_morph q_b.
    pub fn effective_flux(&self, q_b: f64) -> f64 {
        self.f_morph * q_b
    }
    /// Bed resistance: Manning's n from roughness height k_s.
    ///
    /// n = k_s^{1/6} / 26 (empirical, Strickler relation).
    pub fn manning_n(&self) -> f64 {
        self.roughness.powf(1.0 / 6.0) / 26.0
    }
    /// Bed shear stress from Manning's equation (Pa).
    ///
    /// τ_b = ρ g n² u² / h^{1/3}
    pub fn manning_shear_stress(&self, velocity: f64, depth: f64, rho_f: f64) -> f64 {
        let n = self.manning_n();
        if depth < 1e-300 {
            return 0.0;
        }
        rho_f * G * n * n * velocity * velocity / depth.powf(1.0 / 3.0)
    }
}
/// SPH advection-diffusion operator for suspended sediment concentration.
pub struct SphSedimentAdvection {
    /// SPH smoothing length (m).
    pub h: f64,
    /// Settling velocity w_s (m/s).
    pub settling_velocity: f64,
    /// Turbulent diffusivity ε_s (m²/s).
    pub diffusivity: f64,
}
impl SphSedimentAdvection {
    /// Construct an SPH sediment advection-diffusion model.
    pub fn new(h: f64, settling_velocity: f64, diffusivity: f64) -> Self {
        SphSedimentAdvection {
            h,
            settling_velocity,
            diffusivity,
        }
    }
    /// SPH gradient of concentration: ∇C_i = Σ_j m_j/ρ_j (C_j - C_i) ∇W_ij.
    pub fn concentration_gradient(
        &self,
        c_i: f64,
        neighbors: &[([f64; 3], f64, f64, f64)],
    ) -> [f64; 3] {
        let mut grad = [0.0f64; 3];
        for &(r_ij, c_j, m_j, rho_j) in neighbors {
            if rho_j < 1e-300 {
                continue;
            }
            let gw = cubic_kernel_grad_sed(r_ij, self.h);
            let dc = c_j - c_i;
            grad[0] += m_j / rho_j * dc * gw[0];
            grad[1] += m_j / rho_j * dc * gw[1];
            grad[2] += m_j / rho_j * dc * gw[2];
        }
        grad
    }
    /// dC/dt from diffusion: ε_s ∇²C (Brookshaw form).
    pub fn diffusion_rate(&self, c_i: f64, neighbors: &[([f64; 3], f64, f64, f64)]) -> f64 {
        let mut sum = 0.0;
        for &(r_ij, c_j, m_j, rho_j) in neighbors {
            let r2 = dot3_sed(r_ij, r_ij);
            if r2 < 1e-300 || rho_j < 1e-300 {
                continue;
            }
            let gw = cubic_kernel_grad_sed(r_ij, self.h);
            let dot_rg = dot3_sed(r_ij, gw);
            sum += m_j / rho_j * 2.0 * (c_i - c_j) / r2 * dot_rg;
        }
        -self.diffusivity * sum
    }
    /// Settling contribution to dC/dt: -w_s ∂C/∂z ≈ -w_s (∇C)_z.
    pub fn settling_rate(&self, c_i: f64, neighbors: &[([f64; 3], f64, f64, f64)]) -> f64 {
        let grad = self.concentration_gradient(c_i, neighbors);
        -self.settling_velocity * grad[2]
    }
}
/// 1-D depth-averaged sediment concentration transport equation.
///
/// ∂C/∂t + u ∂C/∂x = ε_s ∂²C/∂x² + (E − D) / h
///
/// Discretised with upwind advection and central diffusion on a uniform grid.
pub struct SedimentConcentration1D {
    /// Concentration at each grid cell C_i (kg/m³ or volumetric −).
    pub concentration: Vec<f64>,
    /// Grid spacing Δx (m).
    pub dx: f64,
    /// Depth-averaged flow velocity u (m/s, positive = rightward).
    pub velocity: f64,
    /// Turbulent diffusivity ε_s (m²/s).
    pub diffusivity: f64,
    /// Settling velocity w_s (m/s).
    pub settling_velocity: f64,
    /// Water depth h (m).
    pub depth: f64,
}
impl SedimentConcentration1D {
    /// Construct a sediment concentration model.
    pub fn new(
        n_cells: usize,
        dx: f64,
        velocity: f64,
        diffusivity: f64,
        settling_velocity: f64,
        depth: f64,
    ) -> Self {
        SedimentConcentration1D {
            concentration: vec![0.0; n_cells],
            dx,
            velocity,
            diffusivity,
            settling_velocity,
            depth,
        }
    }
    /// Advance one time step with explicit upwind advection + central diffusion.
    ///
    /// `erosion_rate` and `deposition_rate` are (E − D)/h source terms per cell.
    pub fn advance(&mut self, dt: f64, erosion_rates: &[f64], deposition_rates: &[f64]) {
        let n = self.concentration.len();
        if n < 2 {
            return;
        }
        let mut new_c = self.concentration.clone();
        let u = self.velocity;
        let eps = self.diffusivity;
        let dx = self.dx;
        for i in 1..n - 1 {
            let adv = if u >= 0.0 {
                u * (self.concentration[i] - self.concentration[i - 1]) / dx
            } else {
                u * (self.concentration[i + 1] - self.concentration[i]) / dx
            };
            let diff = eps
                * (self.concentration[i + 1] - 2.0 * self.concentration[i]
                    + self.concentration[i - 1])
                / (dx * dx);
            let src = erosion_rates.get(i).copied().unwrap_or(0.0)
                - deposition_rates.get(i).copied().unwrap_or(0.0);
            new_c[i] = (self.concentration[i] - dt * adv + dt * diff + dt * src).max(0.0);
        }
        self.concentration = new_c;
    }
    /// CFL stability criterion for advection.
    pub fn cfl_condition(&self) -> f64 {
        self.dx / self.velocity.abs().max(1e-300)
    }
    /// Diffusion stability criterion.
    pub fn diffusion_stability(&self) -> f64 {
        self.dx * self.dx / (2.0 * self.diffusivity.max(1e-300))
    }
    /// Maximum stable time step.
    pub fn max_dt(&self) -> f64 {
        self.cfl_condition().min(self.diffusion_stability())
    }
    /// Total mass per unit width (kg/m or −).
    pub fn total_mass(&self) -> f64 {
        self.concentration.iter().sum::<f64>() * self.dx * self.depth
    }
}
/// Suspended sediment load: advection-diffusion with deposition/erosion.
pub struct SuspendedLoad {
    /// Turbulent Schmidt number (typically ~0.7 for sediment).
    pub schmidt_number: f64,
    /// Turbulent diffusivity ε_s (m²/s).
    pub diffusivity: f64,
    /// Settling velocity w_s (m/s).
    pub settling_velocity: f64,
    /// Erosion rate coefficient E_0 (kg/m²/s).
    pub erosion_coeff: f64,
}
impl SuspendedLoad {
    /// Construct a suspended load model.
    pub fn new(
        schmidt_number: f64,
        turbulent_viscosity: f64,
        settling_velocity: f64,
        erosion_coeff: f64,
    ) -> Self {
        let diffusivity = turbulent_viscosity / schmidt_number.max(1e-10);
        SuspendedLoad {
            schmidt_number,
            diffusivity,
            settling_velocity,
            erosion_coeff,
        }
    }
    /// Erosion flux E (kg/m²/s) from Partheniades-type formula:
    ///
    /// E = M * (τ_b / τ_cr − 1)^n for τ_b > τ_cr, else 0.
    pub fn erosion_flux(&self, tau_b: f64, tau_cr: f64, n_exp: f64) -> f64 {
        if tau_b <= tau_cr {
            return 0.0;
        }
        self.erosion_coeff * ((tau_b / tau_cr - 1.0).max(0.0)).powf(n_exp)
    }
    /// Deposition flux D (kg/m²/s) = w_s * C_near_bed.
    pub fn deposition_flux(&self, c_near_bed: f64) -> f64 {
        self.settling_velocity * c_near_bed.max(0.0)
    }
    /// Net sediment flux per unit area: D − E (positive = net deposition).
    pub fn net_flux(&self, c_near_bed: f64, tau_b: f64, tau_cr: f64) -> f64 {
        self.deposition_flux(c_near_bed) - self.erosion_flux(tau_b, tau_cr, 1.0)
    }
}
/// Bed level evolution from net sediment flux.
pub struct BedEvolution {
    /// Bed elevation η (m) at each grid cell.
    pub bed_elevation: Vec<f64>,
    /// Horizontal cell spacings Δx (m).
    pub dx: Vec<f64>,
    /// Porosity of bed material (0–1).
    pub porosity: f64,
    /// Avalanching angle of repose (degrees).
    pub angle_of_repose: f64,
}
impl BedEvolution {
    /// Construct a bed evolution model.
    pub fn new(
        bed_elevation: Vec<f64>,
        dx: Vec<f64>,
        porosity: f64,
        angle_of_repose_deg: f64,
    ) -> Self {
        BedEvolution {
            bed_elevation,
            dx,
            porosity,
            angle_of_repose: angle_of_repose_deg.to_radians(),
        }
    }
    /// Update bed elevation using Exner equation: (1-n) ∂η/∂t = -∂q_b/∂x.
    ///
    /// `q_b` is the bed-load flux at each cell face (m²/s), `dt` time step (s).
    pub fn update(&mut self, q_b: &[f64], dt: f64) {
        let n = self.bed_elevation.len();
        if n == 0 || q_b.len() != n {
            return;
        }
        for i in 0..n {
            let dq = if i == 0 { q_b[0] } else { q_b[i] - q_b[i - 1] };
            let dx = self.dx.get(i).copied().unwrap_or(1.0);
            self.bed_elevation[i] -= dt * dq / (dx * (1.0 - self.porosity));
        }
    }
    /// Apply avalanching: limit bed slope to angle of repose.
    pub fn avalanche(&mut self) {
        let n = self.bed_elevation.len();
        if n < 2 {
            return;
        }
        for _ in 0..10 {
            let mut any_change = false;
            for i in 0..n - 1 {
                let dx = self.dx.get(i).copied().unwrap_or(1.0);
                let dz = self.bed_elevation[i + 1] - self.bed_elevation[i];
                let slope = dz / dx;
                let max_slope = self.angle_of_repose.tan();
                if slope.abs() > max_slope {
                    let correction = (slope.abs() - max_slope) * dx * 0.5;
                    if slope > 0.0 {
                        self.bed_elevation[i] += correction;
                        self.bed_elevation[i + 1] -= correction;
                    } else {
                        self.bed_elevation[i] -= correction;
                        self.bed_elevation[i + 1] += correction;
                    }
                    any_change = true;
                }
            }
            if !any_change {
                break;
            }
        }
    }
    /// Net volume change (m³/m) from initial elevation.
    pub fn volume_change(&self, initial: &[f64]) -> f64 {
        self.bed_elevation
            .iter()
            .zip(initial.iter())
            .zip(self.dx.iter())
            .map(|((z, z0), dx)| (z - z0) * dx)
            .sum()
    }
}
