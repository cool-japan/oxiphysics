// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Tire force models for vehicle dynamics.
//!
//! Provides the [`TireModel`] trait and implementations:
//! - [`PacejkaTire`]: Pacejka "Magic Formula" (lateral, longitudinal, combined slip)
//! - [`FialaTire`]: Fiala brush model with linear region + saturation
//! - [`LinearTire`]: Simple proportional model for testing
//! - [`TireRelaxation`]: First-order relaxation length model
//! - Load sensitivity, camber thrust, self-aligning torque improvements

use oxiphysics_core::math::Real;

/// Trait for computing tire lateral and longitudinal forces.
pub trait TireModel {
    /// Compute lateral and longitudinal forces from slip conditions.
    ///
    /// # Returns
    /// `(lateral_force, longitudinal_force)` in Newtons.
    fn compute_forces(
        &self,
        slip_angle: Real,
        slip_ratio: Real,
        normal_force: Real,
        friction: Real,
    ) -> (Real, Real);
}

// ---------------------------------------------------------------------------
// Pacejka "Magic Formula" tire model
// ---------------------------------------------------------------------------

/// Coefficients for one axis of the Pacejka Magic Formula.
///
/// `F = D * sin(C * atan(B*x - E*(B*x - atan(B*x))))`
#[derive(Debug, Clone)]
pub struct PacejkaCoeffs {
    /// Stiffness factor.
    pub b: Real,
    /// Shape factor.
    pub c: Real,
    /// Peak factor (scales with normal force * friction).
    pub d: Real,
    /// Curvature factor.
    pub e: Real,
}

impl Default for PacejkaCoeffs {
    fn default() -> Self {
        Self {
            b: 10.0,
            c: 1.9,
            d: 1.0,
            e: 0.97,
        }
    }
}

/// Pacejka "Magic Formula" tire model.
///
/// Models both lateral (cornering) and longitudinal (traction/braking) forces
/// with optional combined-slip coupling via friction ellipse.
#[derive(Debug, Clone)]
pub struct PacejkaTire {
    /// Lateral (cornering) coefficients.
    pub lateral: PacejkaCoeffs,
    /// Longitudinal (traction/braking) coefficients.
    pub longitudinal: PacejkaCoeffs,
    /// Enable combined slip weight (0 = no interaction, 1 = full friction ellipse).
    pub combined_slip_weight: Real,
}

impl Default for PacejkaTire {
    fn default() -> Self {
        Self {
            lateral: PacejkaCoeffs {
                b: 10.0,
                c: 1.9,
                d: 1.0,
                e: 0.97,
            },
            longitudinal: PacejkaCoeffs {
                b: 12.0,
                c: 2.3,
                d: 1.0,
                e: 0.97,
            },
            combined_slip_weight: 0.5,
        }
    }
}

impl PacejkaTire {
    /// Create a new Pacejka tire model with default coefficients.
    pub fn new() -> Self {
        Self::default()
    }

    /// Create with custom lateral and longitudinal coefficients.
    pub fn with_coefficients(
        lateral: PacejkaCoeffs,
        longitudinal: PacejkaCoeffs,
        combined_slip_weight: Real,
    ) -> Self {
        Self {
            lateral,
            longitudinal,
            combined_slip_weight,
        }
    }

    /// Evaluate the Magic Formula for a single axis.
    fn magic_formula(coeffs: &PacejkaCoeffs, x: Real, normal_force: Real, friction: Real) -> Real {
        let peak = coeffs.d * normal_force * friction;
        let bx = coeffs.b * x;
        let inner = bx - coeffs.e * (bx - bx.atan());
        peak * (coeffs.c * inner.atan()).sin()
    }

    /// Combined slip via friction ellipse: returns (Fy, Fx).
    ///
    /// The friction ellipse scales the pure-slip forces so that the resultant
    /// magnitude does not exceed mu*Fz.
    #[allow(dead_code)]
    pub fn combined_slip_friction_ellipse(
        &self,
        slip_angle: Real,
        slip_ratio: Real,
        normal_force: Real,
        friction: Real,
    ) -> (Real, Real) {
        if normal_force <= 0.0 {
            return (0.0, 0.0);
        }
        let fy_pure = Self::magic_formula(&self.lateral, slip_angle.tan(), normal_force, friction);
        let fx_pure = Self::magic_formula(&self.longitudinal, slip_ratio, normal_force, friction);

        let max_force = friction * normal_force;
        let resultant = (fy_pure * fy_pure + fx_pure * fx_pure).sqrt();
        if resultant > max_force && resultant > 1e-10 {
            let scale = max_force / resultant;
            (fy_pure * scale, fx_pure * scale)
        } else {
            (fy_pure, fx_pure)
        }
    }
}

impl TireModel for PacejkaTire {
    fn compute_forces(
        &self,
        slip_angle: Real,
        slip_ratio: Real,
        normal_force: Real,
        friction: Real,
    ) -> (Real, Real) {
        if normal_force <= 0.0 {
            return (0.0, 0.0);
        }

        let fy_pure = Self::magic_formula(&self.lateral, slip_angle.tan(), normal_force, friction);
        let fx_pure = Self::magic_formula(&self.longitudinal, slip_ratio, normal_force, friction);

        if self.combined_slip_weight > 0.0 {
            let sigma_lat = slip_angle.tan().abs();
            let sigma_lon = slip_ratio.abs();
            let sigma_total = (sigma_lat * sigma_lat + sigma_lon * sigma_lon).sqrt();

            if sigma_total > 1e-10 {
                let lat_weight = sigma_lat / sigma_total;
                let lon_weight = sigma_lon / sigma_total;

                let w = self.combined_slip_weight;
                let fy = fy_pure * (1.0 - w + w * lat_weight);
                let fx = fx_pure * (1.0 - w + w * lon_weight);
                (fy, fx)
            } else {
                (fy_pure, fx_pure)
            }
        } else {
            (fy_pure, fx_pure)
        }
    }
}

// ---------------------------------------------------------------------------
// Fiala (brush) tire model
// ---------------------------------------------------------------------------

/// Fiala brush tire model.
#[derive(Debug, Clone)]
pub struct FialaTire {
    /// Lateral cornering stiffness (N/rad).
    pub cornering_stiffness: Real,
    /// Longitudinal stiffness (N/unit-slip-ratio).
    pub longitudinal_stiffness: Real,
}

impl Default for FialaTire {
    fn default() -> Self {
        Self {
            cornering_stiffness: 50000.0,
            longitudinal_stiffness: 80000.0,
        }
    }
}

impl FialaTire {
    /// Create a new Fiala tire model with custom stiffness values.
    pub fn new(cornering_stiffness: Real, longitudinal_stiffness: Real) -> Self {
        Self {
            cornering_stiffness,
            longitudinal_stiffness,
        }
    }

    /// Compute force for one axis using Fiala brush model.
    fn fiala_force(stiffness: Real, slip: Real, max_force: Real) -> Real {
        let slip_abs = slip.abs();
        if max_force <= 0.0 {
            return 0.0;
        }
        let s_crit = if stiffness > 1e-10 {
            3.0 * max_force / stiffness
        } else {
            return 0.0;
        };

        let force_magnitude = if slip_abs < s_crit {
            stiffness * slip_abs - stiffness * stiffness * slip_abs * slip_abs / (3.0 * max_force)
                + stiffness * stiffness * stiffness * slip_abs * slip_abs * slip_abs
                    / (27.0 * max_force * max_force)
        } else {
            max_force
        };

        force_magnitude.min(max_force) * slip.signum()
    }
}

impl TireModel for FialaTire {
    fn compute_forces(
        &self,
        slip_angle: Real,
        slip_ratio: Real,
        normal_force: Real,
        friction: Real,
    ) -> (Real, Real) {
        if normal_force <= 0.0 {
            return (0.0, 0.0);
        }
        let max_force = friction * normal_force;

        let fy = Self::fiala_force(self.cornering_stiffness, slip_angle, max_force);
        let fx = Self::fiala_force(self.longitudinal_stiffness, slip_ratio, max_force);

        (fy, fx)
    }
}

// ---------------------------------------------------------------------------
// Linear tire model
// ---------------------------------------------------------------------------

/// Simple linear (proportional) tire model for testing.
#[derive(Debug, Clone)]
pub struct LinearTire {
    /// Lateral cornering stiffness (N/rad).
    pub cornering_stiffness: Real,
    /// Longitudinal stiffness (N/unit-slip-ratio).
    pub longitudinal_stiffness: Real,
}

impl Default for LinearTire {
    fn default() -> Self {
        Self {
            cornering_stiffness: 40000.0,
            longitudinal_stiffness: 60000.0,
        }
    }
}

impl LinearTire {
    /// Create a new linear tire model.
    pub fn new(cornering_stiffness: Real, longitudinal_stiffness: Real) -> Self {
        Self {
            cornering_stiffness,
            longitudinal_stiffness,
        }
    }
}

impl TireModel for LinearTire {
    fn compute_forces(
        &self,
        slip_angle: Real,
        slip_ratio: Real,
        normal_force: Real,
        friction: Real,
    ) -> (Real, Real) {
        if normal_force <= 0.0 {
            return (0.0, 0.0);
        }
        let max_force = friction * normal_force;

        let fy_raw = self.cornering_stiffness * slip_angle;
        let fx_raw = self.longitudinal_stiffness * slip_ratio;

        let fy = fy_raw.clamp(-max_force, max_force);
        let fx = fx_raw.clamp(-max_force, max_force);

        (fy, fx)
    }
}

// ---------------------------------------------------------------------------
// Tire relaxation length model
// ---------------------------------------------------------------------------

/// First-order tire relaxation length model.
///
/// The tire does not develop force instantaneously; the relaxation length
/// controls the lag between slip input and force output.
///
/// dF/dt = (v / sigma) * (F_steady - F)
///
/// where sigma is the relaxation length and v is the forward speed.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct TireRelaxation {
    /// Lateral relaxation length (m).
    pub sigma_lat: f64,
    /// Longitudinal relaxation length (m).
    pub sigma_lon: f64,
    /// Current lateral force (N) — state variable.
    pub fy: f64,
    /// Current longitudinal force (N) — state variable.
    pub fx: f64,
}

#[allow(dead_code)]
impl TireRelaxation {
    /// Create a new relaxation model with typical relaxation lengths.
    pub fn new(sigma_lat: f64, sigma_lon: f64) -> Self {
        Self {
            sigma_lat,
            sigma_lon,
            fy: 0.0,
            fx: 0.0,
        }
    }

    /// Default relaxation lengths for a passenger car tire.
    pub fn default_passenger() -> Self {
        Self::new(0.3, 0.15)
    }

    /// Update the relaxation model given the steady-state forces and dt.
    ///
    /// `speed` is the forward vehicle speed (m/s).
    pub fn step(&mut self, fy_steady: f64, fx_steady: f64, speed: f64, dt: f64) {
        let v = speed.abs().max(0.1); // prevent division by near-zero

        let tau_lat = self.sigma_lat / v;
        let alpha_lat = if tau_lat > 1e-12 {
            (dt / tau_lat).min(1.0)
        } else {
            1.0
        };

        let tau_lon = self.sigma_lon / v;
        let alpha_lon = if tau_lon > 1e-12 {
            (dt / tau_lon).min(1.0)
        } else {
            1.0
        };

        self.fy += alpha_lat * (fy_steady - self.fy);
        self.fx += alpha_lon * (fx_steady - self.fx);
    }

    /// Current relaxed forces: (Fy, Fx).
    pub fn forces(&self) -> (f64, f64) {
        (self.fy, self.fx)
    }

    /// Reset forces to zero.
    pub fn reset(&mut self) {
        self.fy = 0.0;
        self.fx = 0.0;
    }
}

// ---------------------------------------------------------------------------
// Load sensitivity
// ---------------------------------------------------------------------------

/// Load sensitivity model: the friction coefficient decreases with increasing
/// normal load.
///
/// mu(Fz) = mu_ref * (Fz_ref / Fz)^load_sensitivity_exponent
///
/// This captures the well-known effect that tire grip per unit load degrades
/// as the tire is loaded more heavily.
#[allow(dead_code)]
pub struct LoadSensitivity {
    /// Reference friction coefficient at `fz_ref`.
    pub mu_ref: f64,
    /// Reference normal load (N).
    pub fz_ref: f64,
    /// Load sensitivity exponent (typically 0.05..0.15).
    pub exponent: f64,
}

#[allow(dead_code)]
impl LoadSensitivity {
    /// Create a new load sensitivity model.
    pub fn new(mu_ref: f64, fz_ref: f64, exponent: f64) -> Self {
        Self {
            mu_ref,
            fz_ref,
            exponent,
        }
    }

    /// Default load sensitivity for a sport tire.
    pub fn default_sport() -> Self {
        Self {
            mu_ref: 1.4,
            fz_ref: 4000.0,
            exponent: 0.10,
        }
    }

    /// Compute friction coefficient at the given normal load.
    pub fn friction_at_load(&self, fz: f64) -> f64 {
        if fz <= 0.0 || self.fz_ref <= 0.0 {
            return self.mu_ref;
        }
        self.mu_ref * (self.fz_ref / fz).powf(self.exponent)
    }

    /// Compute the peak lateral force (mu * Fz) at the given load.
    pub fn peak_force_at_load(&self, fz: f64) -> f64 {
        self.friction_at_load(fz) * fz
    }
}

// ---------------------------------------------------------------------------
// Camber thrust
// ---------------------------------------------------------------------------

/// Camber thrust model: lateral force contribution from tire camber angle.
///
/// Fy_camber = C_gamma * gamma * Fz / Fz_ref
///
/// where gamma is the camber angle and C_gamma is the camber stiffness.
#[allow(dead_code)]
pub struct CamberThrust {
    /// Camber stiffness (N/rad at reference load).
    pub c_gamma: f64,
    /// Reference normal load (N).
    pub fz_ref: f64,
}

#[allow(dead_code)]
impl CamberThrust {
    /// Create a new camber thrust model.
    pub fn new(c_gamma: f64, fz_ref: f64) -> Self {
        Self { c_gamma, fz_ref }
    }

    /// Default camber thrust for a typical road tire.
    pub fn default_road() -> Self {
        Self {
            c_gamma: 800.0,
            fz_ref: 4000.0,
        }
    }

    /// Compute camber thrust lateral force.
    ///
    /// `gamma` is the camber angle in radians (positive = top tilted outward).
    /// `fz` is the normal load (N).
    pub fn lateral_force(&self, gamma: f64, fz: f64) -> f64 {
        if self.fz_ref <= 0.0 || fz <= 0.0 {
            return 0.0;
        }
        self.c_gamma * gamma * fz / self.fz_ref
    }

    /// Combined lateral force with camber: Fy_total = Fy_slip + Fy_camber.
    pub fn combined_lateral(&self, fy_slip: f64, gamma: f64, fz: f64) -> f64 {
        fy_slip + self.lateral_force(gamma, fz)
    }
}

// ---------------------------------------------------------------------------
// Self-aligning torque
// ---------------------------------------------------------------------------

/// Self-aligning torque model.
///
/// The pneumatic trail decreases with slip angle, producing a torque that
/// tries to reduce the slip angle (hence "self-aligning").
///
/// Mz = -t_p * Fy
///
/// where t_p is the pneumatic trail (m).
#[allow(dead_code)]
pub struct SelfAligningTorque {
    /// Maximum pneumatic trail at zero slip (m).
    pub trail_max: f64,
    /// Slip angle (radians) at which the trail has decayed to half.
    pub alpha_half: f64,
    /// Residual mechanical trail (m).
    pub mechanical_trail: f64,
}

#[allow(dead_code)]
impl SelfAligningTorque {
    /// Create a new self-aligning torque model.
    pub fn new(trail_max: f64, alpha_half: f64, mechanical_trail: f64) -> Self {
        Self {
            trail_max,
            alpha_half,
            mechanical_trail,
        }
    }

    /// Default parameters for a passenger car tire.
    pub fn default_passenger() -> Self {
        Self {
            trail_max: 0.03,
            alpha_half: 0.08, // ~4.6 degrees
            mechanical_trail: 0.005,
        }
    }

    /// Pneumatic trail as a function of slip angle.
    ///
    /// Decays as: t_p = t_max * cos(atan(alpha / alpha_half))
    pub fn pneumatic_trail(&self, slip_angle: f64) -> f64 {
        let ratio = slip_angle.abs() / self.alpha_half;
        self.trail_max * (ratio.atan()).cos()
    }

    /// Total trail (pneumatic + mechanical).
    pub fn total_trail(&self, slip_angle: f64) -> f64 {
        self.pneumatic_trail(slip_angle) + self.mechanical_trail
    }

    /// Compute self-aligning torque (N*m).
    ///
    /// Negative sign: torque acts to reduce slip angle.
    pub fn torque(&self, slip_angle: f64, lateral_force: f64) -> f64 {
        let trail = self.total_trail(slip_angle);
        -trail * lateral_force
    }

    /// Torque with load dependence: trail scales with sqrt(Fz/Fz_ref).
    pub fn torque_with_load(
        &self,
        slip_angle: f64,
        lateral_force: f64,
        fz: f64,
        fz_ref: f64,
    ) -> f64 {
        let load_factor = if fz_ref > 0.0 {
            (fz / fz_ref).sqrt()
        } else {
            1.0
        };
        let trail = self.total_trail(slip_angle) * load_factor;
        -trail * lateral_force
    }
}

// ---------------------------------------------------------------------------
// Combined slip Pacejka (Fx + Fy coupling) — standalone function
// ---------------------------------------------------------------------------

/// Compute combined-slip Pacejka forces using the Similarity Method.
///
/// The combined lateral slip is: sigma_alpha' = tan(alpha) / (1 + kappa)
/// The combined longitudinal slip is: kappa' = kappa / (1 + kappa)
/// The resultant force is limited by the friction ellipse.
#[allow(dead_code)]
pub fn pacejka_combined_slip(
    tire: &PacejkaTire,
    slip_angle: Real,
    slip_ratio: Real,
    normal_force: Real,
    friction: Real,
) -> (Real, Real) {
    if normal_force <= 0.0 {
        return (0.0, 0.0);
    }

    // Similarity method: equivalent slips
    let denom = (1.0 + slip_ratio.abs()).max(1e-10);
    let alpha_eq = slip_angle.tan() / denom;
    let kappa_eq = slip_ratio / denom;

    let fy = PacejkaTire::magic_formula(&tire.lateral, alpha_eq, normal_force, friction);
    let fx = PacejkaTire::magic_formula(&tire.longitudinal, kappa_eq, normal_force, friction);

    // Friction ellipse limit
    let max_force = friction * normal_force;
    let resultant = (fy * fy + fx * fx).sqrt();
    if resultant > max_force && resultant > 1e-10 {
        let s = max_force / resultant;
        (fy * s, fx * s)
    } else {
        (fy, fx)
    }
}

// ---------------------------------------------------------------------------
// Tire temperature influence on stiffness
// ---------------------------------------------------------------------------

/// Models how tire temperature affects the friction coefficient and stiffness.
///
/// Tires have an optimal operating temperature window.  Below optimum the
/// compound is cold and hard (low grip); above it the compound degrades.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct TireTemperatureModel {
    /// Optimal (peak-grip) temperature in °C.
    pub optimal_temp: f64,
    /// Temperature window half-width below which grip is reduced (°C).
    pub cold_window: f64,
    /// Temperature window half-width above which overheating occurs (°C).
    pub hot_window: f64,
    /// Minimum friction factor (at extreme cold).
    pub cold_friction_min: f64,
    /// Minimum friction factor (at overheating).
    pub hot_friction_min: f64,
    /// Thermal capacity (J / °C) — governs how fast the tire heats up.
    pub thermal_capacity: f64,
    /// Cooling coefficient (W / °C) — heat loss to environment.
    pub cooling_coefficient: f64,
}

#[allow(dead_code)]
impl TireTemperatureModel {
    /// Typical sport slick tire parameters.
    pub fn default_slick() -> Self {
        Self {
            optimal_temp: 90.0,
            cold_window: 40.0,
            hot_window: 30.0,
            cold_friction_min: 0.6,
            hot_friction_min: 0.75,
            thermal_capacity: 800.0,
            cooling_coefficient: 50.0,
        }
    }

    /// Friction factor at a given temperature (1.0 = maximum).
    ///
    /// Uses a piecewise quadratic model centred on `optimal_temp`.
    pub fn friction_factor(&self, temperature: f64) -> f64 {
        let delta = temperature - self.optimal_temp;
        if delta < 0.0 {
            // Cold side
            let t = (-delta / self.cold_window).min(1.0);
            1.0 - (1.0 - self.cold_friction_min) * t * t
        } else {
            // Hot side
            let t = (delta / self.hot_window).min(1.0);
            1.0 - (1.0 - self.hot_friction_min) * t * t
        }
    }

    /// Effective stiffness scaling factor.
    ///
    /// Stiffness is highest at the optimal temperature and reduced both
    /// below (cold compound) and above (degraded compound).
    pub fn stiffness_factor(&self, temperature: f64) -> f64 {
        // Stiffness degrades more slowly than friction — use friction^0.5 approx.
        self.friction_factor(temperature).sqrt()
    }

    /// Update tire temperature using a first-order thermal model.
    ///
    /// `dT/dt = (Q_gen - Q_cool) / C_th`
    ///
    /// where `Q_gen = slip_force · slip_speed` (heat generation) and
    /// `Q_cool = k_cool · (T - T_ambient)`.
    ///
    /// `lateral_force` is in N, `slip_speed` is the sliding speed in m/s.
    pub fn update_temperature(&self, current_temp: f64, heat_generation: f64, dt: f64) -> f64 {
        let ambient = 25.0;
        let q_cool = self.cooling_coefficient * (current_temp - ambient).max(0.0);
        let delta_t = (heat_generation - q_cool) / self.thermal_capacity;
        current_temp + delta_t * dt
    }
}

// ---------------------------------------------------------------------------
// Tire wear model
// ---------------------------------------------------------------------------

/// Simple tire wear model that accumulates rubber abrasion.
///
/// Wear rate increases with slip speed and normal force, similar to
/// Archard's law of abrasive wear.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct TireWearModel {
    /// Current wear level (0 = new, 1 = fully worn).
    pub wear_level: f64,
    /// Wear coefficient (dimensionless); higher = faster wear.
    pub wear_coefficient: f64,
    /// Tire hardness factor: softer compounds wear faster.
    pub hardness: f64,
}

#[allow(dead_code)]
impl TireWearModel {
    /// Typical soft-compound tire (fast wear, high grip new).
    pub fn new_soft_compound() -> Self {
        Self {
            wear_level: 0.0,
            wear_coefficient: 1e-9,
            hardness: 0.6,
        }
    }

    /// Typical hard-compound tire (slow wear, consistent grip).
    pub fn new_hard_compound() -> Self {
        Self {
            wear_level: 0.0,
            wear_coefficient: 3e-10,
            hardness: 0.9,
        }
    }

    /// Accumulate wear over a time step.
    ///
    /// `wear_rate = k_wear * Fz * slip_speed / hardness`
    ///
    /// `slip_speed` is the magnitude of the sliding velocity at the
    /// contact patch (m/s).
    pub fn accumulate_wear(&mut self, normal_force: f64, slip_speed: f64, dt: f64) {
        let rate = self.wear_coefficient * normal_force * slip_speed / self.hardness;
        self.wear_level = (self.wear_level + rate * dt).min(1.0);
    }

    /// Remaining tire life as a fraction.
    pub fn remaining_life(&self) -> f64 {
        (1.0 - self.wear_level).max(0.0)
    }

    /// Friction multiplier based on current wear.
    ///
    /// Fresh tire has full grip; worn tire has reduced grip due to loss of
    /// microscopic texture.  Modelled as a quadratic decay.
    pub fn friction_multiplier(&self) -> f64 {
        let w = self.wear_level.min(1.0);
        1.0 - 0.3 * w * w
    }

    /// Effective stiffness reduction from wear (thinning of the carcass).
    pub fn stiffness_multiplier(&self) -> f64 {
        let w = self.wear_level.min(1.0);
        1.0 - 0.15 * w
    }
}

// ---------------------------------------------------------------------------
// Transient brush model
// ---------------------------------------------------------------------------

/// Transient (dynamic) brush tire model.
///
/// The contact patch deforms elastically before the rubber slides.  The
/// transient slip state `σ` evolves as:
///
/// `dσ/ds = (σ_input − σ) / σ_length`
///
/// where `s` is arc-length travelled, `σ_input` is the kinematic slip, and
/// `σ_length` is the contact patch half-length.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct TransientBrushModel {
    /// Lateral cornering stiffness (N/rad).
    pub cornering_stiffness: f64,
    /// Longitudinal stiffness (N/unit-slip-ratio).
    pub longitudinal_stiffness: f64,
    /// Lateral relaxation length / contact half-length (m).
    pub sigma_lat: f64,
    /// Longitudinal relaxation length (m).
    pub sigma_lon: f64,
}

/// State of the transient brush model.
#[allow(dead_code)]
#[derive(Debug, Clone, Default)]
pub struct TransientBrushState {
    /// Current transient lateral slip.
    pub slip_lat: f64,
    /// Current transient longitudinal slip.
    pub slip_lon: f64,
}

#[allow(dead_code)]
impl TransientBrushModel {
    /// Create a transient brush model with typical passenger-car parameters.
    pub fn new(
        cornering_stiffness: f64,
        longitudinal_stiffness: f64,
        sigma_lat: f64,
        sigma_lon: f64,
    ) -> Self {
        Self {
            cornering_stiffness,
            longitudinal_stiffness,
            sigma_lat,
            sigma_lon,
        }
    }
}

#[allow(dead_code)]
impl TransientBrushState {
    /// Update the transient slip state.
    ///
    /// `slip_angle_input` and `slip_ratio_input` are the kinematic slip inputs
    /// (from vehicle kinematics).  `speed` is the forward speed (m/s) and
    /// `dt` is the time step (s).
    pub fn update(
        &mut self,
        slip_angle_input: f64,
        slip_ratio_input: f64,
        speed: f64,
        dt: f64,
        model: &TransientBrushModel,
    ) {
        let v = speed.abs().max(0.1);
        let alpha_lat = (v * dt / model.sigma_lat).min(1.0);
        let alpha_lon = (v * dt / model.sigma_lon).min(1.0);
        self.slip_lat += alpha_lat * (slip_angle_input - self.slip_lat);
        self.slip_lon += alpha_lon * (slip_ratio_input - self.slip_lon);
    }

    /// Compute lateral and longitudinal forces from the current transient slip.
    pub fn forces(
        &self,
        model: &TransientBrushModel,
        normal_force: f64,
        friction: f64,
    ) -> (f64, f64) {
        let max_force = friction * normal_force;
        let fy_raw = model.cornering_stiffness * self.slip_lat;
        let fx_raw = model.longitudinal_stiffness * self.slip_lon;
        let fy = fy_raw.clamp(-max_force, max_force);
        let fx = fx_raw.clamp(-max_force, max_force);
        (fy, fx)
    }
}

// ---------------------------------------------------------------------------
// Loaded radius computation
// ---------------------------------------------------------------------------

/// Computes the loaded (dynamic) rolling radius of a tire as a function of
/// normal force.
///
/// The loaded radius decreases under load due to carcass deflection.
///
/// `r_load = r_free - Fz / (2 · k_radial)`
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct LoadedRadiusModel {
    /// Unloaded (free) rolling radius in m.
    pub free_radius: f64,
    /// Radial stiffness in N/m.
    pub radial_stiffness: f64,
    /// Minimum allowed loaded radius (structural limit) in m.
    pub min_radius: f64,
}

#[allow(dead_code)]
impl LoadedRadiusModel {
    /// Create a model for a typical 245/40 R18 passenger tire.
    pub fn new_passenger_tire() -> Self {
        Self {
            free_radius: 0.330,
            radial_stiffness: 200_000.0,
            min_radius: 0.270,
        }
    }

    /// Loaded radius at normal force `fz` (N).
    pub fn loaded_radius(&self, fz: f64) -> f64 {
        let deflection = fz / (2.0 * self.radial_stiffness);
        (self.free_radius - deflection).max(self.min_radius)
    }

    /// Contact patch half-length approximation.
    ///
    /// `a ≈ sqrt(2 · r_free · deflection)` — simple Hertzian geometry.
    pub fn contact_half_length(&self, fz: f64) -> f64 {
        let deflection = (self.free_radius - self.loaded_radius(fz)).max(0.0);
        (2.0 * self.free_radius * deflection).sqrt()
    }

    /// Effective rolling circumference.
    pub fn rolling_circumference(&self, fz: f64) -> f64 {
        2.0 * std::f64::consts::PI * self.loaded_radius(fz)
    }
}

// ---------------------------------------------------------------------------
// Contact patch pressure distribution
// ---------------------------------------------------------------------------

/// Contact patch pressure distribution model.
///
/// Approximates the 2-D pressure distribution in the contact patch as an
/// elliptical (Hertzian) profile.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct ContactPatchPressure {
    /// Free rolling radius in m.
    pub free_radius: f64,
    /// Tire width in m.
    pub tire_width: f64,
    /// Radial stiffness N/m.
    pub radial_stiffness: f64,
    /// Inflation pressure in Pa.
    pub inflation_pressure: f64,
}

#[allow(dead_code)]
impl ContactPatchPressure {
    /// Typical passenger tire (245/40 R18).
    pub fn new_passenger_tire() -> Self {
        Self {
            free_radius: 0.330,
            tire_width: 0.245,
            radial_stiffness: 200_000.0,
            inflation_pressure: 220_000.0, // 2.2 bar
        }
    }

    /// Approximate contact patch area (m²).
    ///
    /// `A_contact ≈ Fz / p_inflation`
    pub fn contact_area(&self, fz: f64) -> f64 {
        if self.inflation_pressure < 1e-10 {
            return 0.0;
        }
        fz / self.inflation_pressure
    }

    /// Contact patch dimensions `(length, width)` in m.
    ///
    /// Assumes an elliptical footprint with aspect ratio ~1.5:1 (length:width).
    pub fn patch_dimensions(&self, fz: f64) -> (f64, f64) {
        let area = self.contact_area(fz);
        let aspect = 1.5_f64;
        // A = π * a * b with b = a / aspect → A = π * a² / aspect
        let a = (area * aspect / std::f64::consts::PI).sqrt(); // half-length
        let b = a / aspect; // half-width
        (2.0 * a, 2.0 * b)
    }

    /// Peak and average pressure in the contact patch (Pa).
    ///
    /// For an elliptical Hertzian distribution, peak = 1.5 × average.
    pub fn peak_and_average_pressure(&self, fz: f64) -> (f64, f64) {
        let area = self.contact_area(fz).max(1e-10);
        let avg = fz / area;
        let peak = 1.5 * avg;
        (peak, avg)
    }

    /// Lateral shear stress at the patch edge.
    ///
    /// Approximates the traction coefficient at the edge of the contact ellipse.
    pub fn edge_shear_stress(&self, fz: f64, friction: f64) -> f64 {
        let (_, avg) = self.peak_and_average_pressure(fz);
        friction * avg
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_pacejka_peak_at_reasonable_slip_angle() {
        let tire = PacejkaTire::new();
        let normal = 5000.0;
        let friction = 1.0;

        let mut max_force: Real = 0.0;
        let mut peak_angle: Real = 0.0;
        for i in 1..200 {
            let angle = i as Real * 0.001;
            let (fy, _) = tire.compute_forces(angle, 0.0, normal, friction);
            if fy.abs() > max_force {
                max_force = fy.abs();
                peak_angle = angle;
            }
        }

        assert!(
            peak_angle > 0.02 && peak_angle < 0.3,
            "peak at {peak_angle} rad"
        );
        assert!(max_force > 0.5 * normal, "peak force {max_force} too low");
    }

    #[test]
    fn test_pacejka_zero_normal_force() {
        let tire = PacejkaTire::new();
        let (fy, fx) = tire.compute_forces(0.1, 0.1, 0.0, 1.0);
        assert!((fy).abs() < 1e-10);
        assert!((fx).abs() < 1e-10);
    }

    #[test]
    fn test_fiala_linear_region() {
        let tire = FialaTire::new(50000.0, 80000.0);
        let normal = 5000.0;
        let friction = 1.0;

        let small_angle = 0.001;
        let (fy, _) = tire.compute_forces(small_angle, 0.0, normal, friction);
        let expected = 50000.0 * small_angle;
        assert!(
            (fy - expected).abs() / expected < 0.05,
            "fy={fy}, expected~{expected}"
        );
    }

    #[test]
    fn test_fiala_saturation() {
        let tire = FialaTire::new(50000.0, 80000.0);
        let normal = 5000.0;
        let friction = 1.0;
        let max_force = friction * normal;

        let (fy, _) = tire.compute_forces(1.0, 0.0, normal, friction);
        assert!(
            (fy.abs() - max_force).abs() < 1.0,
            "fy={fy}, max={max_force}"
        );
    }

    #[test]
    fn test_linear_tire_proportional() {
        let tire = LinearTire::new(40000.0, 60000.0);
        let normal = 5000.0;
        let friction = 1.0;

        let (fy, fx) = tire.compute_forces(0.01, 0.05, normal, friction);
        assert!((fy - 400.0).abs() < 1e-10);
        assert!((fx - 3000.0).abs() < 1e-10);
    }

    #[test]
    fn test_linear_tire_clamped_to_friction() {
        let tire = LinearTire::new(40000.0, 60000.0);
        let normal = 1000.0;
        let friction = 1.0;
        let max_force = friction * normal;

        let (fy, _) = tire.compute_forces(1.0, 0.0, normal, friction);
        assert!(
            (fy.abs() - max_force).abs() < 1e-10,
            "fy={fy}, max={max_force}"
        );
    }

    #[test]
    fn test_pacejka_peak_slip() {
        let tire = PacejkaTire::new();
        let normal = 5000.0;
        let friction = 1.0;

        let mut max_force: Real = 0.0;
        let mut peak_slip: Real = 0.0;
        for i in 1..=50 {
            let slip = i as Real * 0.01;
            let (_, fx) = tire.compute_forces(0.0, slip, normal, friction);
            if fx.abs() > max_force {
                max_force = fx.abs();
                peak_slip = slip;
            }
        }

        assert!(
            (0.05..=0.30).contains(&peak_slip),
            "peak_slip={peak_slip}, expected 0.05..0.30"
        );
        assert!(
            max_force > 0.5 * normal,
            "peak longitudinal force {max_force} too low (mu*N={normal})"
        );
    }

    #[test]
    fn test_pacejka_zero_slip() {
        let tire = PacejkaTire::new();
        let normal = 5000.0;
        let friction = 1.0;

        let (fy, fx) = tire.compute_forces(0.0, 0.0, normal, friction);
        assert!(
            fx.abs() < 1e-6,
            "longitudinal force at zero slip should be 0, got {fx}"
        );
        assert!(
            fy.abs() < 1e-6,
            "lateral force at zero slip angle should be 0, got {fy}"
        );
    }

    // --- Combined slip friction ellipse tests ---

    #[test]
    fn test_combined_slip_friction_ellipse() {
        let tire = PacejkaTire::new();
        let normal = 5000.0;
        let friction = 1.0;
        let (fy, fx) = tire.combined_slip_friction_ellipse(0.1, 0.1, normal, friction);
        let resultant = (fy * fy + fx * fx).sqrt();
        let max_force = friction * normal;
        assert!(
            resultant <= max_force + 1e-6,
            "resultant {resultant} should not exceed mu*Fz={max_force}"
        );
    }

    #[test]
    fn test_combined_slip_friction_ellipse_zero_normal() {
        let tire = PacejkaTire::new();
        let (fy, fx) = tire.combined_slip_friction_ellipse(0.1, 0.1, 0.0, 1.0);
        assert!(fy.abs() < 1e-10);
        assert!(fx.abs() < 1e-10);
    }

    // --- Pacejka combined slip (similarity method) tests ---

    #[test]
    fn test_pacejka_combined_slip_within_ellipse() {
        let tire = PacejkaTire::new();
        let normal = 5000.0;
        let friction = 1.0;
        let (fy, fx) = pacejka_combined_slip(&tire, 0.1, 0.15, normal, friction);
        let resultant = (fy * fy + fx * fx).sqrt();
        let max_force = friction * normal;
        assert!(
            resultant <= max_force + 1e-6,
            "resultant {resultant} exceeds friction limit {max_force}"
        );
    }

    #[test]
    fn test_pacejka_combined_slip_pure_lateral() {
        let tire = PacejkaTire::new();
        let normal = 5000.0;
        let friction = 1.0;
        // Zero slip ratio → should behave like pure lateral
        let (fy, fx) = pacejka_combined_slip(&tire, 0.05, 0.0, normal, friction);
        assert!(fy.abs() > 0.0, "should have lateral force");
        assert!(fx.abs() < 1e-6, "no longitudinal force at zero slip ratio");
    }

    // --- Tire relaxation tests ---

    #[test]
    fn test_relaxation_converges() {
        let mut relax = TireRelaxation::new(0.3, 0.15);
        let fy_steady = 3000.0;
        let fx_steady = 1500.0;
        let speed = 20.0;

        // Run for many steps
        for _ in 0..1000 {
            relax.step(fy_steady, fx_steady, speed, 0.001);
        }
        let (fy, fx) = relax.forces();
        assert!(
            (fy - fy_steady).abs() < 1.0,
            "should converge to steady state: fy={fy}"
        );
        assert!(
            (fx - fx_steady).abs() < 1.0,
            "should converge to steady state: fx={fx}"
        );
    }

    #[test]
    fn test_relaxation_lag() {
        let mut relax = TireRelaxation::new(0.3, 0.15);
        relax.step(3000.0, 1500.0, 20.0, 0.001);
        let (fy, fx) = relax.forces();
        // After one tiny step, force should be small (still lagging)
        assert!(
            fy < 3000.0,
            "after one step, fy should lag behind steady state"
        );
        assert!(
            fx < 1500.0,
            "after one step, fx should lag behind steady state"
        );
    }

    #[test]
    fn test_relaxation_reset() {
        let mut relax = TireRelaxation::new(0.3, 0.15);
        relax.step(3000.0, 1500.0, 20.0, 0.01);
        relax.reset();
        let (fy, fx) = relax.forces();
        assert!(fy.abs() < 1e-10);
        assert!(fx.abs() < 1e-10);
    }

    #[test]
    fn test_relaxation_high_speed_faster() {
        let mut relax_slow = TireRelaxation::new(0.3, 0.15);
        let mut relax_fast = TireRelaxation::new(0.3, 0.15);

        relax_slow.step(3000.0, 1500.0, 5.0, 0.01);
        relax_fast.step(3000.0, 1500.0, 50.0, 0.01);

        // Higher speed → faster convergence
        assert!(
            relax_fast.fy > relax_slow.fy,
            "higher speed should converge faster"
        );
    }

    // --- Load sensitivity tests ---

    #[test]
    fn test_load_sensitivity_at_reference() {
        let ls = LoadSensitivity::new(1.4, 4000.0, 0.10);
        let mu = ls.friction_at_load(4000.0);
        assert!((mu - 1.4).abs() < 1e-10, "at ref load, mu should be mu_ref");
    }

    #[test]
    fn test_load_sensitivity_degrades() {
        let ls = LoadSensitivity::new(1.4, 4000.0, 0.10);
        let mu_heavy = ls.friction_at_load(8000.0);
        let mu_ref = ls.friction_at_load(4000.0);
        assert!(
            mu_heavy < mu_ref,
            "heavier load should reduce friction: {mu_heavy} vs {mu_ref}"
        );
    }

    #[test]
    fn test_load_sensitivity_light_load() {
        let ls = LoadSensitivity::new(1.4, 4000.0, 0.10);
        let mu_light = ls.friction_at_load(2000.0);
        let mu_ref = ls.friction_at_load(4000.0);
        assert!(
            mu_light > mu_ref,
            "lighter load should increase friction: {mu_light} vs {mu_ref}"
        );
    }

    #[test]
    fn test_load_sensitivity_peak_force() {
        let ls = LoadSensitivity::new(1.4, 4000.0, 0.10);
        let f1 = ls.peak_force_at_load(4000.0);
        let f2 = ls.peak_force_at_load(8000.0);
        // Even though friction drops, total force should still increase with load
        assert!(
            f2 > f1,
            "peak force should increase with load: {f2} vs {f1}"
        );
        // But less than linear: f2 < 2 * f1
        assert!(f2 < 2.0 * f1, "load sensitivity means sub-linear scaling");
    }

    // --- Camber thrust tests ---

    #[test]
    fn test_camber_thrust_zero_camber() {
        let ct = CamberThrust::new(800.0, 4000.0);
        let fy = ct.lateral_force(0.0, 4000.0);
        assert!(fy.abs() < 1e-10);
    }

    #[test]
    fn test_camber_thrust_positive_camber() {
        let ct = CamberThrust::new(800.0, 4000.0);
        let fy = ct.lateral_force(0.05, 4000.0); // ~2.9 degrees
        // Fy = 800 * 0.05 * 4000/4000 = 40 N
        assert!((fy - 40.0).abs() < 1e-10);
    }

    #[test]
    fn test_camber_thrust_load_scaling() {
        let ct = CamberThrust::new(800.0, 4000.0);
        let fy_ref = ct.lateral_force(0.05, 4000.0);
        let fy_heavy = ct.lateral_force(0.05, 8000.0);
        assert!(
            (fy_heavy - 2.0 * fy_ref).abs() < 1e-10,
            "force should scale linearly with load"
        );
    }

    #[test]
    fn test_camber_thrust_combined() {
        let ct = CamberThrust::new(800.0, 4000.0);
        let fy_slip = 2000.0;
        let fy_total = ct.combined_lateral(fy_slip, 0.05, 4000.0);
        assert!((fy_total - 2040.0).abs() < 1e-10);
    }

    // --- Self-aligning torque tests ---

    #[test]
    fn test_sat_max_trail_at_zero_slip() {
        let sat = SelfAligningTorque::default_passenger();
        let trail = sat.pneumatic_trail(0.0);
        // At zero slip, trail should be close to max
        assert!(
            (trail - sat.trail_max).abs() < 1e-10,
            "trail at zero slip: got {trail}"
        );
    }

    #[test]
    fn test_sat_trail_decreases_with_slip() {
        let sat = SelfAligningTorque::default_passenger();
        let t0 = sat.pneumatic_trail(0.0);
        let t1 = sat.pneumatic_trail(0.1);
        assert!(
            t1 < t0,
            "trail should decrease with slip angle: {t1} vs {t0}"
        );
    }

    #[test]
    fn test_sat_torque_sign() {
        let sat = SelfAligningTorque::default_passenger();
        // Positive slip angle with positive lateral force → negative torque
        let mz = sat.torque(0.05, 2000.0);
        assert!(
            mz < 0.0,
            "self-aligning torque should oppose slip: got {mz}"
        );
    }

    #[test]
    fn test_sat_torque_with_load() {
        let sat = SelfAligningTorque::default_passenger();
        let mz_ref = sat.torque_with_load(0.05, 2000.0, 4000.0, 4000.0);
        let mz_heavy = sat.torque_with_load(0.05, 2000.0, 16000.0, 4000.0);
        // Heavier load → larger trail → larger torque magnitude
        assert!(
            mz_heavy.abs() > mz_ref.abs(),
            "heavier load should give larger SAT"
        );
        // Load factor = sqrt(16000/4000) = 2, so |mz_heavy| ≈ 2 * |mz_ref|
        assert!(
            (mz_heavy.abs() - 2.0 * mz_ref.abs()).abs() < 1e-6,
            "mz_heavy={mz_heavy}, mz_ref={mz_ref}"
        );
    }

    #[test]
    fn test_sat_total_trail() {
        let sat = SelfAligningTorque::new(0.03, 0.08, 0.005);
        let total = sat.total_trail(0.0);
        assert!(
            (total - 0.035).abs() < 1e-10,
            "total = pneumatic + mechanical"
        );
    }

    // ── TireTemperatureModel ────────────────────────────────────────────────

    #[test]
    fn test_temperature_friction_at_optimal() {
        let model = TireTemperatureModel::default_slick();
        let mu = model.friction_factor(model.optimal_temp);
        assert!(
            (mu - 1.0).abs() < 1e-9,
            "at optimal temp friction factor should be 1.0"
        );
    }

    #[test]
    fn test_temperature_cold_reduces_friction() {
        let model = TireTemperatureModel::default_slick();
        let mu_cold = model.friction_factor(20.0);
        let mu_opt = model.friction_factor(model.optimal_temp);
        assert!(
            mu_cold < mu_opt,
            "cold tire should have less grip: {mu_cold} vs {mu_opt}"
        );
    }

    #[test]
    fn test_temperature_heating_converges() {
        let model = TireTemperatureModel::default_slick();
        let mut temp = 20.0_f64;
        // Use a higher heat generation to overcome cooling and reach a substantial temperature
        for _ in 0..5000 {
            temp = model.update_temperature(temp, 5000.0, 0.01);
        }
        assert!(
            temp > 50.0,
            "temperature should rise significantly under heavy load: {temp}"
        );
    }

    // ── TireWearModel ───────────────────────────────────────────────────────

    #[test]
    fn test_wear_starts_at_zero() {
        let wear = TireWearModel::new_soft_compound();
        assert!((wear.wear_level).abs() < 1e-10);
    }

    #[test]
    fn test_wear_increases_with_slip() {
        let mut wear = TireWearModel::new_soft_compound();
        wear.accumulate_wear(3000.0, 0.15, 1.0);
        assert!(wear.wear_level > 0.0, "wear should accumulate under slip");
    }

    #[test]
    fn test_wear_capped_at_one() {
        let mut wear = TireWearModel::new_soft_compound();
        for _ in 0..10000 {
            wear.accumulate_wear(5000.0, 0.5, 0.1);
        }
        assert!(wear.wear_level <= 1.0, "wear level must not exceed 1.0");
    }

    #[test]
    fn test_worn_tire_reduced_friction() {
        let wear_new = TireWearModel {
            wear_level: 0.0,
            ..TireWearModel::new_soft_compound()
        };
        let wear_worn = TireWearModel {
            wear_level: 0.8,
            ..TireWearModel::new_soft_compound()
        };
        assert!(
            wear_worn.friction_multiplier() < wear_new.friction_multiplier(),
            "worn tire should have less friction"
        );
    }

    // ── TransientBrushModel ─────────────────────────────────────────────────

    #[test]
    fn test_transient_brush_converges_to_steady() {
        let brush = TransientBrushModel::new(50000.0, 80000.0, 0.3, 0.15);
        let mut state = TransientBrushState::default();
        let target_lat = 0.1_f64;
        for _ in 0..500 {
            state.update(target_lat, 0.0, 20.0, 0.001, &brush);
        }
        // After many steps the transient slip should be close to target
        assert!(
            (state.slip_lat - target_lat).abs() < 0.01,
            "transient slip should converge: got {}",
            state.slip_lat
        );
    }

    #[test]
    fn test_transient_brush_lags_step_input() {
        let brush = TransientBrushModel::new(50000.0, 80000.0, 0.3, 0.15);
        let mut state = TransientBrushState::default();
        state.update(0.2, 0.0, 20.0, 0.001, &brush);
        assert!(
            state.slip_lat < 0.2,
            "after one step the transient slip should lag behind target"
        );
    }

    // ── LoadedRadiusModel ───────────────────────────────────────────────────

    #[test]
    fn test_loaded_radius_decreases_under_load() {
        let model = LoadedRadiusModel::new_passenger_tire();
        let r_unloaded = model.loaded_radius(0.0);
        let r_loaded = model.loaded_radius(5000.0);
        assert!(
            r_loaded < r_unloaded,
            "loaded radius should shrink under normal force"
        );
    }

    #[test]
    fn test_loaded_radius_positive() {
        let model = LoadedRadiusModel::new_passenger_tire();
        assert!(model.loaded_radius(0.0) > 0.0);
        assert!(model.loaded_radius(10000.0) > 0.0);
    }

    // ── ContactPatchPressure ─────────────────────────────────────────────────

    #[test]
    fn test_contact_patch_area_increases_with_load() {
        let model = ContactPatchPressure::new_passenger_tire();
        let a1 = model.contact_area(3000.0);
        let a2 = model.contact_area(6000.0);
        assert!(a2 > a1, "heavier load → larger contact area");
    }

    #[test]
    fn test_contact_patch_peak_pressure() {
        let model = ContactPatchPressure::new_passenger_tire();
        let (peak, avg) = model.peak_and_average_pressure(4000.0);
        assert!(peak > avg, "peak pressure should exceed average");
        assert!(avg > 0.0, "average pressure must be positive");
    }

    #[test]
    fn test_contact_patch_dimensions_positive() {
        let model = ContactPatchPressure::new_passenger_tire();
        let (l, w) = model.patch_dimensions(4000.0);
        assert!(
            l > 0.0 && w > 0.0,
            "contact patch must have positive dimensions"
        );
    }
}
