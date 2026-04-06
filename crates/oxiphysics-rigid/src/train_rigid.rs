// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Railway/train rigid-body dynamics.
//!
//! This module provides physically-based models for railway vehicle dynamics
//! from a rigid-body perspective, including:
//!
//! - **Hertzian wheel–rail contact** – elliptical contact patch and normal force
//! - **Creepage forces** – Kalker linear theory and FASTSIM tangential forces
//! - **Bogie dynamics** – two-axle bogie with primary/secondary suspension
//! - **Track geometry** – cant (superelevation), curvature, transition curves
//! - **Hunting oscillation** – lateral stability and critical speed
//! - **Ride comfort indices** – ISO 2631 / Sperling Wz
//! - **Braking distance** – service and emergency braking on grade
//! - **Traction curves** – adhesion-limited tractive effort vs speed
//! - **Coupler forces** – draft gear / buffer spring model
//! - **Derailment criteria** – Nadal L/V ratio and wheel unloading

#![allow(dead_code)]

use std::f64::consts::PI;

// ═══════════════════════════════════════════════════════════════════════════════
// Constants
// ═══════════════════════════════════════════════════════════════════════════════

/// Standard gravitational acceleration (m/s²).
const G: f64 = 9.81;

/// Standard rail gauge (m) — UIC / 1435 mm.
const STANDARD_GAUGE: f64 = 1.435;

/// Typical wheel radius for passenger rolling stock (m).
const DEFAULT_WHEEL_RADIUS: f64 = 0.46;

/// Typical wheel conicity (tread angle tangent) for new S1002 profile.
const DEFAULT_CONICITY: f64 = 0.05;

// ═══════════════════════════════════════════════════════════════════════════════
// Hertzian wheel–rail contact
// ═══════════════════════════════════════════════════════════════════════════════

/// Parameters for Hertzian wheel–rail contact calculation.
#[derive(Debug, Clone)]
pub struct HertzianContactParams {
    /// Wheel radius (m).
    pub wheel_radius: f64,
    /// Rail crown radius (m) — transverse curvature of rail head.
    pub rail_crown_radius: f64,
    /// Combined elastic modulus E* = E / (2(1 - nu²)) for steel-on-steel (Pa).
    pub combined_modulus: f64,
    /// Poisson ratio of wheel/rail material (dimensionless).
    pub poisson_ratio: f64,
}

impl Default for HertzianContactParams {
    fn default() -> Self {
        let e = 210.0e9; // Steel Young's modulus (Pa)
        let nu = 0.3;
        let e_star = e / (2.0 * (1.0 - nu * nu));
        Self {
            wheel_radius: DEFAULT_WHEEL_RADIUS,
            rail_crown_radius: 0.30,
            combined_modulus: e_star,
            poisson_ratio: nu,
        }
    }
}

/// Result of a Hertzian contact computation.
#[derive(Debug, Clone)]
pub struct HertzianContactResult {
    /// Semi-axis of the contact ellipse along the rolling direction (m).
    pub a: f64,
    /// Semi-axis of the contact ellipse in the lateral direction (m).
    pub b: f64,
    /// Maximum contact pressure (Pa).
    pub p0: f64,
    /// Contact area (m²).
    pub area: f64,
    /// Elastic approach / penetration depth (m).
    pub penetration: f64,
}

/// Compute Hertzian contact between a cylindrical wheel and a crowned rail.
///
/// Uses the simplified circular contact model (a ≈ b) which is acceptable for
/// many wheel–rail geometries. For a given normal force `normal_force` (N),
/// returns the contact patch geometry and peak pressure.
pub fn hertzian_contact(
    params: &HertzianContactParams,
    normal_force: f64,
) -> HertzianContactResult {
    // Equivalent radius: 1/R = 1/R_wheel + 1/R_rail
    let r_eq = (params.wheel_radius * params.rail_crown_radius)
        / (params.wheel_radius + params.rail_crown_radius);

    // Circular contact radius: a = (3 F R / (4 E*))^(1/3)
    let a = (3.0 * normal_force * r_eq / (4.0 * params.combined_modulus))
        .max(0.0)
        .cbrt();

    // For wheel–rail the ellipse is slightly elongated along rolling direction.
    // Use an aspect ratio factor based on typical geometries.
    let aspect = 1.1;
    let b = a / aspect;

    let area = PI * a * b;

    // Peak pressure: p0 = 3F / (2 pi a b) = 1.5 F / (pi a b)
    let p0 = if area > 1e-20 {
        1.5 * normal_force / area
    } else {
        0.0
    };

    // Elastic approach: delta = a² / R_eq
    let penetration = if r_eq > 1e-10 { a * a / r_eq } else { 0.0 };

    HertzianContactResult {
        a,
        b,
        p0,
        area,
        penetration,
    }
}

// ═══════════════════════════════════════════════════════════════════════════════
// Creepage and tangential forces
// ═══════════════════════════════════════════════════════════════════════════════

/// Creepage values at the wheel–rail contact patch.
#[derive(Debug, Clone, Copy, Default)]
pub struct Creepage {
    /// Longitudinal creepage (dimensionless): (v_wheel - v_rail) / v.
    pub xi: f64,
    /// Lateral creepage (dimensionless).
    pub eta: f64,
    /// Spin creepage (1/m).
    pub phi: f64,
}

/// Result of creep force calculation.
#[derive(Debug, Clone, Copy)]
pub struct CreepForceResult {
    /// Longitudinal tangential force (N).
    pub fx: f64,
    /// Lateral tangential force (N).
    pub fy: f64,
    /// Spin moment (N·m).
    pub mz: f64,
}

/// Kalker linear creep theory coefficients for circular contact.
///
/// Returns `(f11, f22, f23)` — the Kalker coefficients scaled by contact patch
/// size and shear modulus, giving forces directly when multiplied by creepages.
///
/// * `a`, `b` — contact ellipse semi-axes (m).
/// * `shear_modulus` — G (Pa), typically ~80 GPa for steel.
pub fn kalker_linear_coefficients(a: f64, b: f64, shear_modulus: f64) -> (f64, f64, f64) {
    // Kalker's tabulated coefficients for a/b ≈ 1 (circular contact):
    // C11 ≈ 4.12, C22 ≈ 3.67, C23 ≈ 1.47
    let c11 = 4.12;
    let c22 = 3.67;
    let c23 = 1.47;

    let ab = a * b;
    let f11 = shear_modulus * ab * c11;
    let f22 = shear_modulus * ab * c22;
    let f23 = shear_modulus * ab.sqrt() * ab * c23;

    (f11, f22, f23)
}

/// Compute creep forces using Kalker's linear theory.
///
/// Valid for small creepages. For larger creepages, use [`fastsim_creep_forces`].
pub fn kalker_creep_forces(
    creepage: &Creepage,
    a: f64,
    b: f64,
    shear_modulus: f64,
) -> CreepForceResult {
    let (f11, f22, f23) = kalker_linear_coefficients(a, b, shear_modulus);

    let fx = -f11 * creepage.xi;
    let fy = -f22 * creepage.eta - f23 * creepage.phi;
    let mz = f23 * creepage.eta - f22 * creepage.phi * (a * b);

    CreepForceResult { fx, fy, mz }
}

/// Simplified FASTSIM tangential force computation.
///
/// Uses a saturation law to limit the linear Kalker forces to the friction
/// cone `mu * N`. This provides a nonlinear creep force model suitable for
/// traction and braking studies.
///
/// * `creepage` — the creepage triple.
/// * `a`, `b` — contact ellipse semi-axes (m).
/// * `shear_modulus` — G (Pa).
/// * `normal_force` — contact normal force N (N).
/// * `mu` — friction coefficient.
pub fn fastsim_creep_forces(
    creepage: &Creepage,
    a: f64,
    b: f64,
    shear_modulus: f64,
    normal_force: f64,
    mu: f64,
) -> CreepForceResult {
    let linear = kalker_creep_forces(creepage, a, b, shear_modulus);

    let f_linear = (linear.fx * linear.fx + linear.fy * linear.fy).sqrt();
    let f_limit = mu * normal_force;

    if f_linear < 1e-20 {
        return CreepForceResult {
            fx: 0.0,
            fy: 0.0,
            mz: 0.0,
        };
    }

    // Johnson–Vermeulen saturation curve
    let ratio = f_linear / f_limit;
    let factor = if ratio <= 3.0 {
        let r3 = ratio / 3.0;
        let r3_sq = r3 * r3;
        ratio * (1.0 - r3 + r3_sq / 3.0) / ratio
    } else {
        f_limit / f_linear
    };

    let saturated_f = f_linear * factor;
    let scale = if f_linear > 1e-20 {
        saturated_f / f_linear
    } else {
        0.0
    };

    CreepForceResult {
        fx: linear.fx * scale,
        fy: linear.fy * scale,
        mz: linear.mz * scale,
    }
}

// ═══════════════════════════════════════════════════════════════════════════════
// Track geometry
// ═══════════════════════════════════════════════════════════════════════════════

/// Description of track geometry at a given chainage.
#[derive(Debug, Clone, Copy)]
pub struct TrackGeometry {
    /// Horizontal curvature (1/m). Positive = curves left.
    pub curvature: f64,
    /// Cant (superelevation) angle (radians).
    pub cant: f64,
    /// Gradient (rise/run, positive uphill).
    pub gradient: f64,
    /// Track gauge (m).
    pub gauge: f64,
}

impl Default for TrackGeometry {
    fn default() -> Self {
        Self {
            curvature: 0.0,
            cant: 0.0,
            gradient: 0.0,
            gauge: STANDARD_GAUGE,
        }
    }
}

/// Compute cant (superelevation) for equilibrium speed on a curve.
///
/// * `speed` — vehicle speed (m/s).
/// * `radius` — curve radius (m).
/// * `gauge` — track gauge (m).
///
/// Returns the cant in metres (height difference between outer and inner rail).
pub fn equilibrium_cant(speed: f64, radius: f64, gauge: f64) -> f64 {
    if radius.abs() < 1e-6 {
        return 0.0;
    }
    gauge * speed * speed / (G * radius.abs())
}

/// Compute the equilibrium cant angle (radians).
pub fn equilibrium_cant_angle(speed: f64, radius: f64) -> f64 {
    if radius.abs() < 1e-6 {
        return 0.0;
    }
    (speed * speed / (G * radius.abs())).atan()
}

/// Compute the cant deficiency (m) — the difference between required and
/// actual cant.
///
/// * `speed` — vehicle speed (m/s).
/// * `radius` — curve radius (m).
/// * `actual_cant` — installed cant (m).
/// * `gauge` — track gauge (m).
pub fn cant_deficiency(speed: f64, radius: f64, actual_cant: f64, gauge: f64) -> f64 {
    let eq_cant = equilibrium_cant(speed, radius, gauge);
    eq_cant - actual_cant
}

/// Compute lateral acceleration experienced by passengers due to cant
/// deficiency (m/s²).
///
/// * `speed` — vehicle speed (m/s).
/// * `radius` — curve radius (m).
/// * `cant_angle` — actual cant angle (radians).
pub fn uncompensated_lateral_acceleration(speed: f64, radius: f64, cant_angle: f64) -> f64 {
    if radius.abs() < 1e-6 {
        return 0.0;
    }
    speed * speed / radius.abs() * cant_angle.cos() - G * cant_angle.sin()
}

/// Clothoid (Euler spiral) transition curve: curvature varies linearly with
/// arc length.
///
/// * `s` — arc length from start of transition (m).
/// * `transition_length` — total length of transition (m).
/// * `target_curvature` — curvature at end of transition (1/m).
///
/// Returns curvature at position `s`.
pub fn clothoid_curvature(s: f64, transition_length: f64, target_curvature: f64) -> f64 {
    if transition_length < 1e-10 {
        return target_curvature;
    }
    let t = (s / transition_length).clamp(0.0, 1.0);
    t * target_curvature
}

/// Compute the x,y offset from a clothoid transition at arc length `s`.
///
/// Uses the Fresnel integral approximation (first few terms of the series).
///
/// * `s` — arc length (m).
/// * `a_param` — clothoid parameter A = sqrt(R * L) (m).
///
/// Returns `[x, y]` offset from the transition origin.
pub fn clothoid_offset(s: f64, a_param: f64) -> [f64; 2] {
    if a_param < 1e-10 {
        return [s, 0.0];
    }
    let tau = s * s / (2.0 * a_param * a_param);
    // Fresnel-like series for x and y
    let x = s * (1.0 - tau * tau / 10.0 + tau.powi(4) / 216.0);
    let y = s * (tau / 3.0 - tau.powi(3) / 42.0 + tau.powi(5) / 1320.0);
    [x, y]
}

// ═══════════════════════════════════════════════════════════════════════════════
// Bogie dynamics
// ═══════════════════════════════════════════════════════════════════════════════

/// Parameters for a two-axle bogie.
#[derive(Debug, Clone)]
pub struct BogieParams {
    /// Bogie frame mass (kg).
    pub frame_mass: f64,
    /// Wheelset mass (kg), per wheelset.
    pub wheelset_mass: f64,
    /// Semi-spacing between wheelsets (m). Full spacing = 2 * semi_spacing.
    pub semi_spacing: f64,
    /// Wheel radius (m).
    pub wheel_radius: f64,
    /// Wheel conicity (tread angle tangent).
    pub conicity: f64,
    /// Primary suspension lateral stiffness per side (N/m).
    pub primary_ky: f64,
    /// Primary suspension vertical stiffness per side (N/m).
    pub primary_kz: f64,
    /// Primary suspension lateral damping per side (N·s/m).
    pub primary_cy: f64,
    /// Secondary suspension lateral stiffness per side (N/m).
    pub secondary_ky: f64,
    /// Secondary suspension vertical stiffness per side (N/m).
    pub secondary_kz: f64,
    /// Secondary suspension lateral damping per side (N·s/m).
    pub secondary_cy: f64,
    /// Track gauge (m).
    pub gauge: f64,
}

impl Default for BogieParams {
    fn default() -> Self {
        Self {
            frame_mass: 2500.0,
            wheelset_mass: 1200.0,
            semi_spacing: 1.25,
            wheel_radius: DEFAULT_WHEEL_RADIUS,
            conicity: DEFAULT_CONICITY,
            primary_ky: 5.0e6,
            primary_kz: 1.0e6,
            primary_cy: 20_000.0,
            secondary_ky: 0.3e6,
            secondary_kz: 0.5e6,
            secondary_cy: 50_000.0,
            gauge: STANDARD_GAUGE,
        }
    }
}

/// State vector for lateral dynamics of a single wheelset.
#[derive(Debug, Clone, Copy, Default)]
pub struct WheelsetState {
    /// Lateral displacement (m).
    pub y: f64,
    /// Yaw angle (rad).
    pub psi: f64,
    /// Lateral velocity (m/s).
    pub dy: f64,
    /// Yaw rate (rad/s).
    pub dpsi: f64,
}

/// State of a two-axle bogie (leading and trailing wheelsets + frame).
#[derive(Debug, Clone, Default)]
pub struct BogieState {
    /// Leading wheelset.
    pub leading: WheelsetState,
    /// Trailing wheelset.
    pub trailing: WheelsetState,
    /// Bogie frame lateral displacement (m).
    pub frame_y: f64,
    /// Bogie frame yaw angle (rad).
    pub frame_psi: f64,
    /// Frame lateral velocity (m/s).
    pub frame_dy: f64,
    /// Frame yaw rate (rad/s).
    pub frame_dpsi: f64,
}

/// Compute lateral gravitational stiffness of a wheelset on the rail.
///
/// Due to the conical tread profile, lateral displacement causes a restoring
/// force proportional to conicity, weight, and the inverse of the wheel radius.
///
/// * `conicity` — effective tread conicity (lambda).
/// * `weight` — axle load (N).
/// * `wheel_radius` — wheel radius (m).
pub fn gravitational_stiffness(conicity: f64, weight: f64, wheel_radius: f64) -> f64 {
    if wheel_radius < 1e-10 {
        return 0.0;
    }
    2.0 * conicity * weight / wheel_radius
}

/// Compute the equivalent lateral creep stiffness for Kalker linear theory.
///
/// * `f22` — Kalker lateral creep coefficient (N).
/// * `speed` — vehicle speed (m/s).
///
/// Returns the effective stiffness contribution (N/m) for lateral dynamics.
pub fn lateral_creep_stiffness(f22: f64, speed: f64) -> f64 {
    if speed.abs() < 1e-6 {
        return 0.0;
    }
    f22 / speed
}

// ═══════════════════════════════════════════════════════════════════════════════
// Hunting oscillation
// ═══════════════════════════════════════════════════════════════════════════════

/// Parameters for hunting oscillation analysis.
#[derive(Debug, Clone)]
pub struct HuntingParams {
    /// Equivalent conicity (dimensionless).
    pub conicity: f64,
    /// Wheel radius (m).
    pub wheel_radius: f64,
    /// Track gauge (half-gauge = distance from center to rail, m).
    pub half_gauge: f64,
    /// Wheelset mass (kg).
    pub wheelset_mass: f64,
    /// Yaw moment of inertia of wheelset (kg·m²).
    pub yaw_inertia: f64,
    /// Primary yaw stiffness (N·m/rad).
    pub yaw_stiffness: f64,
    /// Lateral creep coefficient f22 (N).
    pub f22: f64,
}

impl Default for HuntingParams {
    fn default() -> Self {
        Self {
            conicity: DEFAULT_CONICITY,
            wheel_radius: DEFAULT_WHEEL_RADIUS,
            half_gauge: STANDARD_GAUGE / 2.0,
            wheelset_mass: 1200.0,
            yaw_inertia: 700.0,
            yaw_stiffness: 5.0e6,
            f22: 10.0e6,
        }
    }
}

/// Klingel frequency (rad/s) — the kinematic hunting frequency of a conical
/// wheelset without creep forces.
///
/// omega = v * sqrt(2 * lambda / (r0 * gauge_half))
pub fn klingel_frequency(speed: f64, conicity: f64, wheel_radius: f64, half_gauge: f64) -> f64 {
    if wheel_radius < 1e-10 || half_gauge < 1e-10 {
        return 0.0;
    }
    let arg = 2.0 * conicity / (wheel_radius * half_gauge);
    if arg < 0.0 {
        return 0.0;
    }
    speed * arg.sqrt()
}

/// Klingel wavelength (m) — the spatial wavelength of kinematic hunting.
///
/// L = 2 * pi * sqrt(r0 * gauge_half / (2 * lambda))
pub fn klingel_wavelength(conicity: f64, wheel_radius: f64, half_gauge: f64) -> f64 {
    if conicity < 1e-10 {
        return f64::INFINITY;
    }
    2.0 * PI * (wheel_radius * half_gauge / (2.0 * conicity)).sqrt()
}

/// Estimate the critical speed for hunting instability (m/s).
///
/// Uses a simplified formula from the eigenvalue analysis of a single
/// wheelset with creep forces and primary yaw restraint:
///
/// V_crit ≈ sqrt(k_yaw * r0 * gauge / (2 * f22 * lambda))
pub fn critical_hunting_speed(params: &HuntingParams) -> f64 {
    let num = params.yaw_stiffness * params.wheel_radius * 2.0 * params.half_gauge;
    let den = 2.0 * params.f22 * params.conicity;
    if den < 1e-10 {
        return f64::INFINITY;
    }
    (num / den).sqrt()
}

/// Check if a vehicle speed exceeds the critical hunting speed.
pub fn is_hunting_unstable(speed: f64, params: &HuntingParams) -> bool {
    speed > critical_hunting_speed(params)
}

// ═══════════════════════════════════════════════════════════════════════════════
// Ride comfort
// ═══════════════════════════════════════════════════════════════════════════════

/// Compute the Sperling ride comfort index (Wz) from a set of acceleration
/// samples.
///
/// * `accelerations` — lateral or vertical acceleration samples (m/s²).
/// * `sample_rate` — sampling frequency (Hz).
///
/// Returns the Wz index. Values < 2.5 are considered comfortable; > 3.25 is
/// very uncomfortable.
pub fn sperling_wz(accelerations: &[f64], sample_rate: f64) -> f64 {
    if accelerations.is_empty() || sample_rate < 1e-6 {
        return 0.0;
    }

    // RMS acceleration
    let n = accelerations.len() as f64;
    let rms = (accelerations.iter().map(|a| a * a).sum::<f64>() / n).sqrt();

    // Sperling Wz = 0.896 * (a_rms / f^0.277)^(1/3)
    // Use dominant frequency estimate from zero-crossing rate or default 1 Hz
    let _period = 2.0 * n / sample_rate;
    // Simplified: assume dominant frequency ~ sample_rate / (2*N) * zero_crossings
    let zero_crossings = accelerations
        .windows(2)
        .filter(|w| w[0] * w[1] < 0.0)
        .count() as f64;
    let dominant_freq = if n > 1.0 {
        (zero_crossings * sample_rate / (2.0 * n)).max(0.5)
    } else {
        1.0
    };

    // Frequency weighting factor (B filter)
    let b_weight = frequency_weight_b(dominant_freq);
    let weighted_rms = rms * b_weight;

    // Wz = (a^3 * 0.588)^0.1 — simplified power law
    // More standard: Wz = 4.42 * (a_rms_weighted)^(1/3)
    4.42 * weighted_rms.cbrt()
}

/// B-weighting filter magnitude for ride comfort (simplified).
///
/// Approximation of the Sperling weighting function for vertical vibrations.
fn frequency_weight_b(freq: f64) -> f64 {
    // The B-filter peaks around 5–6 Hz and rolls off at low and high frequencies.
    let f2 = freq * freq;
    let num = f2;
    let den = f2 + 0.4 * 0.4;
    let high_roll = 1.0 / (1.0 + (freq / 20.0).powi(4));
    (num / den).sqrt() * high_roll
}

/// ISO 2631 weighted RMS acceleration for ride quality.
///
/// * `accelerations` — acceleration samples (m/s²).
///
/// Returns the RMS acceleration in m/s². The ISO 2631 comfort scale:
/// - < 0.315: not uncomfortable
/// - 0.315–0.63: a little uncomfortable
/// - 0.5–1.0: fairly uncomfortable
/// - > 2.0: extremely uncomfortable
pub fn iso2631_rms(accelerations: &[f64]) -> f64 {
    if accelerations.is_empty() {
        return 0.0;
    }
    let n = accelerations.len() as f64;
    (accelerations.iter().map(|a| a * a).sum::<f64>() / n).sqrt()
}

// ═══════════════════════════════════════════════════════════════════════════════
// Braking distance
// ═══════════════════════════════════════════════════════════════════════════════

/// Compute braking distance from initial speed to final speed.
///
/// Uses the energy balance: 0.5 * m * (v0² - vf²) = F_brake * d - m*g*grade*d
///
/// * `mass` — train mass (kg).
/// * `initial_speed` — speed at brake application (m/s).
/// * `final_speed` — target speed (m/s), usually 0.
/// * `brake_force` — average braking force (N).
/// * `gradient` — track gradient (positive uphill). Downhill grade increases distance.
/// * `resistance` — additional running resistance force (N), always retarding.
///
/// Returns the braking distance in metres. Returns `f64::INFINITY` if braking
/// is insufficient to decelerate on the given grade.
#[allow(clippy::too_many_arguments)]
pub fn braking_distance(
    mass: f64,
    initial_speed: f64,
    final_speed: f64,
    brake_force: f64,
    gradient: f64,
    resistance: f64,
) -> f64 {
    let v0_sq = initial_speed * initial_speed;
    let vf_sq = final_speed * final_speed;
    let delta_ke = 0.5 * mass * (v0_sq - vf_sq);
    let net_retarding = brake_force + resistance + mass * G * gradient;

    if net_retarding < 1e-6 {
        return f64::INFINITY;
    }

    (delta_ke / net_retarding).max(0.0)
}

/// Emergency braking distance with reaction time.
///
/// * `mass` — train mass (kg).
/// * `speed` — speed at driver awareness (m/s).
/// * `brake_force` — emergency braking force (N).
/// * `reaction_time` — time before brakes are fully applied (s).
/// * `gradient` — track gradient.
pub fn emergency_braking_distance(
    mass: f64,
    speed: f64,
    brake_force: f64,
    reaction_time: f64,
    gradient: f64,
) -> f64 {
    let reaction_distance = speed * reaction_time;
    let stopping_distance = braking_distance(mass, speed, 0.0, brake_force, gradient, 0.0);
    reaction_distance + stopping_distance
}

/// Compute time to stop from given speed.
///
/// * `mass` — train mass (kg).
/// * `speed` — initial speed (m/s).
/// * `brake_force` — constant braking force (N).
/// * `gradient` — track gradient.
pub fn braking_time(mass: f64, speed: f64, brake_force: f64, gradient: f64) -> f64 {
    let net_decel = (brake_force + mass * G * gradient) / mass;
    if net_decel < 1e-10 {
        return f64::INFINITY;
    }
    speed / net_decel
}

// ═══════════════════════════════════════════════════════════════════════════════
// Traction curves
// ═══════════════════════════════════════════════════════════════════════════════

/// Traction motor characteristics.
#[derive(Debug, Clone)]
pub struct TractionMotor {
    /// Continuous rated power (W).
    pub rated_power: f64,
    /// Maximum tractive effort at standstill (N).
    pub max_effort: f64,
    /// Speed at which power limit takes over from effort limit (m/s).
    pub transition_speed: f64,
}

impl TractionMotor {
    /// Create a new traction motor.
    pub fn new(rated_power: f64, max_effort: f64) -> Self {
        let transition_speed = if max_effort > 1e-6 {
            rated_power / max_effort
        } else {
            0.0
        };
        Self {
            rated_power,
            max_effort,
            transition_speed,
        }
    }

    /// Compute the available tractive effort (N) at a given speed (m/s).
    ///
    /// Below the transition speed, effort is constant at `max_effort`.
    /// Above it, effort = power / speed (hyperbolic).
    pub fn tractive_effort(&self, speed: f64) -> f64 {
        let v = speed.abs().max(0.01);
        if v <= self.transition_speed {
            self.max_effort
        } else {
            self.rated_power / v
        }
    }
}

/// Adhesion-limited tractive effort.
///
/// * `normal_force` — vertical force per wheel (N).
/// * `mu` — wheel–rail friction coefficient.
/// * `num_driven_axles` — number of powered axles.
pub fn adhesion_limit(normal_force: f64, mu: f64, num_driven_axles: u32) -> f64 {
    mu * normal_force * num_driven_axles as f64
}

/// Curtius–Kniffler formula for available adhesion coefficient.
///
/// mu = 0.161 + 7.5 / (44 + 3.6 * V_km_h)
///
/// * `speed_kmh` — speed in km/h.
pub fn curtius_kniffler_adhesion(speed_kmh: f64) -> f64 {
    let v = speed_kmh.abs();
    0.161 + 7.5 / (44.0 + 3.6 * v)
}

/// Generate a traction curve: tractive effort vs speed.
///
/// Returns a vector of `(speed_m_s, effort_N)` tuples.
pub fn traction_curve(
    motor: &TractionMotor,
    axle_load: f64,
    num_driven_axles: u32,
    max_speed: f64,
    num_points: usize,
) -> Vec<(f64, f64)> {
    let mut points = Vec::with_capacity(num_points);
    let n = num_points.max(2);
    for i in 0..n {
        let speed = (i as f64 / (n - 1) as f64) * max_speed;
        let motor_effort = motor.tractive_effort(speed);
        let speed_kmh = speed * 3.6;
        let mu = curtius_kniffler_adhesion(speed_kmh);
        let adhesion = adhesion_limit(axle_load, mu, num_driven_axles);
        let effort = motor_effort.min(adhesion);
        points.push((speed, effort));
    }
    points
}

// ═══════════════════════════════════════════════════════════════════════════════
// Coupler forces
// ═══════════════════════════════════════════════════════════════════════════════

/// Draft gear / coupler model parameters.
#[derive(Debug, Clone)]
pub struct CouplerParams {
    /// Stiffness in compression (N/m).
    pub stiffness_compression: f64,
    /// Stiffness in tension (N/m).
    pub stiffness_tension: f64,
    /// Damping coefficient (N·s/m).
    pub damping: f64,
    /// Free slack (m) — gap before the coupler engages.
    pub slack: f64,
    /// Maximum compressive force before failure (N).
    pub max_compression: f64,
    /// Maximum tensile force before failure (N).
    pub max_tension: f64,
}

impl Default for CouplerParams {
    fn default() -> Self {
        Self {
            stiffness_compression: 10.0e6,
            stiffness_tension: 8.0e6,
            damping: 100_000.0,
            slack: 0.01,
            max_compression: 2_000_000.0,
            max_tension: 1_500_000.0,
        }
    }
}

/// Compute coupler force between two cars.
///
/// * `displacement` — relative displacement (positive = extension, negative = compression).
/// * `velocity` — relative velocity along the coupler axis (m/s).
/// * `params` — coupler parameters.
///
/// Returns the force in Newtons (positive = tension, negative = compression).
pub fn coupler_force(displacement: f64, velocity: f64, params: &CouplerParams) -> f64 {
    let eff_disp = if displacement.abs() < params.slack {
        0.0
    } else if displacement > 0.0 {
        displacement - params.slack
    } else {
        displacement + params.slack
    };

    let spring_force = if eff_disp >= 0.0 {
        eff_disp * params.stiffness_tension
    } else {
        eff_disp * params.stiffness_compression
    };

    let damping_force = velocity * params.damping;

    let total = spring_force + damping_force;

    // Clamp to failure limits
    total.clamp(-params.max_compression, params.max_tension)
}

/// Compute coupler forces for an entire train consist.
///
/// * `positions` — car center positions along the track (m).
/// * `velocities` — car velocities (m/s).
/// * `car_lengths` — lengths of each car (m).
/// * `params` — coupler parameters (same for all couplers).
///
/// Returns a vector of coupler forces, length = `n_cars - 1`.
pub fn train_coupler_forces(
    positions: &[f64],
    velocities: &[f64],
    car_lengths: &[f64],
    params: &CouplerParams,
) -> Vec<f64> {
    let n = positions.len();
    if n < 2 {
        return vec![];
    }

    let mut forces = Vec::with_capacity(n - 1);
    for i in 0..n - 1 {
        let nominal_gap = (car_lengths[i] + car_lengths[i + 1]) / 2.0;
        let actual_gap = positions[i + 1] - positions[i];
        let displacement = actual_gap - nominal_gap;
        let rel_velocity = velocities[i + 1] - velocities[i];
        forces.push(coupler_force(displacement, rel_velocity, params));
    }
    forces
}

// ═══════════════════════════════════════════════════════════════════════════════
// Derailment criteria
// ═══════════════════════════════════════════════════════════════════════════════

/// Nadal derailment criterion: L/V ratio.
///
/// Computes the ratio of lateral (L) to vertical (V) wheel force.
/// If L/V exceeds the critical value, derailment is likely.
///
/// * `lateral_force` — lateral force on the wheel (N).
/// * `vertical_force` — vertical force on the wheel (N).
///
/// Returns the L/V ratio.
pub fn nadal_lv_ratio(lateral_force: f64, vertical_force: f64) -> f64 {
    if vertical_force.abs() < 1e-6 {
        return f64::INFINITY;
    }
    lateral_force.abs() / vertical_force.abs()
}

/// Critical L/V ratio from Nadal's formula.
///
/// L/V_crit = (tan(alpha) - mu) / (1 + mu * tan(alpha))
///
/// * `flange_angle` — wheel flange angle (radians), typically 60–70°.
/// * `mu` — friction coefficient.
pub fn nadal_critical_lv(flange_angle: f64, mu: f64) -> f64 {
    let tan_a = flange_angle.tan();
    (tan_a - mu) / (1.0 + mu * tan_a)
}

/// Check if the Nadal derailment criterion is exceeded.
pub fn check_nadal_derailment(
    lateral_force: f64,
    vertical_force: f64,
    flange_angle: f64,
    mu: f64,
) -> bool {
    let lv = nadal_lv_ratio(lateral_force, vertical_force);
    let critical = nadal_critical_lv(flange_angle, mu);
    lv > critical
}

/// Weinstock wheel unloading criterion.
///
/// Derailment risk increases when the vertical load on one wheel drops
/// significantly below the static load. The ratio Q/Q0 < 0.1 is dangerous.
///
/// * `dynamic_load` — current vertical wheel load (N).
/// * `static_load` — nominal static wheel load (N).
///
/// Returns the load ratio Q/Q0.
pub fn wheel_unloading_ratio(dynamic_load: f64, static_load: f64) -> f64 {
    if static_load.abs() < 1e-6 {
        return 0.0;
    }
    dynamic_load / static_load
}

/// Check if wheel unloading exceeds the danger threshold.
///
/// * `dynamic_load` — current vertical wheel load (N).
/// * `static_load` — nominal static wheel load (N).
/// * `threshold` — minimum acceptable ratio (typically 0.1).
pub fn check_wheel_unloading(dynamic_load: f64, static_load: f64, threshold: f64) -> bool {
    wheel_unloading_ratio(dynamic_load, static_load) < threshold
}

/// Prud'homme limit for lateral track force (N).
///
/// The lateral force on the track must not exceed:
/// F_lat <= alpha + beta * Q_static   (where alpha = 10 kN, beta = 1/3)
///
/// * `static_axle_load` — static axle load (N).
///
/// Returns the maximum allowable lateral track force (N).
pub fn prudhomme_limit(static_axle_load: f64) -> f64 {
    10_000.0 + static_axle_load / 3.0
}

// ═══════════════════════════════════════════════════════════════════════════════
// Slope / grade helpers
// ═══════════════════════════════════════════════════════════════════════════════

/// Convert grade from per-mille (‰) to ratio (rise/run).
pub fn permille_to_ratio(permille: f64) -> f64 {
    permille / 1000.0
}

/// Convert grade from ratio (rise/run) to per-mille (‰).
pub fn ratio_to_permille(ratio: f64) -> f64 {
    ratio * 1000.0
}

/// Compute the gravitational component of resistance force on a gradient.
///
/// * `mass` — vehicle mass (kg).
/// * `gradient` — grade as ratio (rise/run, positive uphill).
pub fn grade_resistance(mass: f64, gradient: f64) -> f64 {
    mass * G * gradient
}

// ═══════════════════════════════════════════════════════════════════════════════
// Davis equation resistance
// ═══════════════════════════════════════════════════════════════════════════════

/// Davis equation running resistance (N).
///
/// R = A + B*v + C*v²
///
/// * `mass` — mass (kg), used to derive default A if coefficients are per-tonne.
/// * `speed` — speed (m/s).
/// * `a`, `b`, `c` — Davis coefficients (N, N·s/m, N·s²/m²).
#[allow(clippy::too_many_arguments)]
pub fn davis_resistance_force(mass: f64, speed: f64, a: f64, b: f64, c: f64) -> f64 {
    let _mass = mass; // included for API symmetry / future per-tonne scaling
    let v = speed.abs();
    a + b * v + c * v * v
}

// ═══════════════════════════════════════════════════════════════════════════════
// Train simulation step
// ═══════════════════════════════════════════════════════════════════════════════

/// State of a single train car for longitudinal dynamics.
#[derive(Debug, Clone)]
pub struct TrainCarState {
    /// Mass (kg).
    pub mass: f64,
    /// Position along track (m).
    pub position: f64,
    /// Velocity (m/s).
    pub velocity: f64,
    /// Length (m).
    pub length: f64,
}

/// Simulate one timestep of longitudinal train dynamics.
///
/// Updates positions and velocities using semi-implicit Euler.
///
/// * `cars` — mutable slice of car states.
/// * `traction_forces` — external traction/braking force on each car (N).
/// * `gradient` — track gradient (rise/run).
/// * `coupler_params` — coupler model.
/// * `dt` — timestep (s).
pub fn train_longitudinal_step(
    cars: &mut [TrainCarState],
    traction_forces: &[f64],
    gradient: f64,
    coupler_params: &CouplerParams,
    dt: f64,
) {
    let n = cars.len();
    if n == 0 {
        return;
    }

    // Compute coupler forces
    let positions: Vec<f64> = cars.iter().map(|c| c.position).collect();
    let velocities: Vec<f64> = cars.iter().map(|c| c.velocity).collect();
    let lengths: Vec<f64> = cars.iter().map(|c| c.length).collect();
    let cforces = train_coupler_forces(&positions, &velocities, &lengths, coupler_params);

    // Update each car
    for i in 0..n {
        let mut net_force = traction_forces.get(i).copied().unwrap_or(0.0);

        // Grade resistance
        net_force -= grade_resistance(cars[i].mass, gradient);

        // Davis resistance (simplified)
        let resistance = davis_resistance_force(cars[i].mass, cars[i].velocity, 500.0, 5.0, 0.6);
        net_force -= resistance * cars[i].velocity.signum();

        // Coupler forces: coupler[i-1] pulls forward, coupler[i] pulls backward
        if i > 0 {
            net_force += cforces[i - 1];
        }
        if i < n - 1 {
            net_force -= cforces[i];
        }

        let accel = net_force / cars[i].mass;
        cars[i].velocity += accel * dt;
        cars[i].position += cars[i].velocity * dt;
    }
}

// ═══════════════════════════════════════════════════════════════════════════════
// Curve negotiation
// ═══════════════════════════════════════════════════════════════════════════════

/// Maximum permissible speed through a curve.
///
/// Based on uncompensated lateral acceleration limit.
///
/// * `radius` — curve radius (m).
/// * `cant` — actual cant (m).
/// * `gauge` — track gauge (m).
/// * `max_cant_deficiency` — maximum allowed cant deficiency (m).
///
/// Returns maximum speed in m/s.
pub fn max_curve_speed(radius: f64, cant: f64, gauge: f64, max_cant_deficiency: f64) -> f64 {
    if radius.abs() < 1e-6 || gauge < 1e-6 {
        return 0.0;
    }
    let total_cant = cant + max_cant_deficiency;
    (G * radius.abs() * total_cant / gauge).sqrt()
}

/// Compute the centrifugal force on a vehicle in a curve.
///
/// * `mass` — vehicle mass (kg).
/// * `speed` — speed (m/s).
/// * `radius` — curve radius (m).
pub fn centrifugal_force(mass: f64, speed: f64, radius: f64) -> f64 {
    if radius.abs() < 1e-6 {
        return 0.0;
    }
    mass * speed * speed / radius.abs()
}

/// Compute lateral load transfer in a curve (N) — the difference in vertical
/// load between outer and inner wheels.
///
/// * `mass` — vehicle mass (kg).
/// * `speed` — speed (m/s).
/// * `radius` — curve radius (m).
/// * `cant_angle` — actual cant angle (rad).
/// * `cg_height` — centre of gravity height (m).
/// * `gauge` — track gauge (m).
pub fn lateral_load_transfer(
    mass: f64,
    speed: f64,
    radius: f64,
    cant_angle: f64,
    cg_height: f64,
    gauge: f64,
) -> f64 {
    if gauge < 1e-6 || radius.abs() < 1e-6 {
        return 0.0;
    }
    let centripetal = mass * speed * speed / radius.abs();
    let uncompensated = centripetal * cant_angle.cos() - mass * G * cant_angle.sin();
    uncompensated * cg_height / gauge
}

// ═══════════════════════════════════════════════════════════════════════════════
// Wheel–rail profile
// ═══════════════════════════════════════════════════════════════════════════════

/// Compute the rolling radius difference (Delta_r) for a conical wheelset at
/// lateral displacement `y`.
///
/// Delta_r = 2 * lambda * y
///
/// * `conicity` — effective conicity (lambda).
/// * `y` — lateral displacement from centre (m).
pub fn rolling_radius_difference(conicity: f64, y: f64) -> f64 {
    2.0 * conicity * y
}

/// Compute the rolling radii of left and right wheels at lateral displacement `y`.
///
/// * `r0` — nominal wheel radius (m).
/// * `conicity` — effective conicity (lambda).
/// * `y` — lateral displacement from centre (m, positive = towards left rail).
///
/// Returns `(r_left, r_right)`.
pub fn wheel_radii(r0: f64, conicity: f64, y: f64) -> (f64, f64) {
    let dr = conicity * y;
    (r0 + dr, r0 - dr)
}

/// Compute effective conicity from measured rolling radius difference data.
///
/// * `delta_r` — measured rolling radius difference (m).
/// * `y` — lateral displacement at which it was measured (m).
pub fn effective_conicity(delta_r: f64, y: f64) -> f64 {
    if y.abs() < 1e-10 {
        return 0.0;
    }
    delta_r / (2.0 * y)
}

// ═══════════════════════════════════════════════════════════════════════════════
// Wheelset yaw torque
// ═══════════════════════════════════════════════════════════════════════════════

/// Compute the longitudinal creep-induced yaw torque on a wheelset due to
/// rolling radius difference.
///
/// When the wheelset is displaced laterally, the different rolling radii on
/// left and right wheels create different circumferential speeds, generating
/// a yaw moment.
///
/// * `f11` — Kalker longitudinal creep coefficient (N).
/// * `delta_r` — rolling radius difference (m).
/// * `speed` — forward speed (m/s).
/// * `wheel_radius` — nominal wheel radius (m).
/// * `half_gauge` — half of the track gauge (m).
pub fn creep_yaw_torque(
    f11: f64,
    delta_r: f64,
    speed: f64,
    wheel_radius: f64,
    half_gauge: f64,
) -> f64 {
    if speed.abs() < 1e-10 || wheel_radius < 1e-10 {
        return 0.0;
    }
    // Longitudinal creepage from rolling radius difference
    let xi = delta_r / (2.0 * wheel_radius);
    // Force on one wheel
    let f_long = f11 * xi;
    // Yaw torque = force * lever arm (half-gauge)
    f_long * half_gauge
}

// ═══════════════════════════════════════════════════════════════════════════════
// Speed conversion helpers
// ═══════════════════════════════════════════════════════════════════════════════

/// Convert speed from km/h to m/s.
pub fn kmh_to_ms(kmh: f64) -> f64 {
    kmh / 3.6
}

/// Convert speed from m/s to km/h.
pub fn ms_to_kmh(ms: f64) -> f64 {
    ms * 3.6
}

// ═══════════════════════════════════════════════════════════════════════════════
// Tests
// ═══════════════════════════════════════════════════════════════════════════════

#[cfg(test)]
mod tests {
    use super::*;

    const TOL: f64 = 1e-6;

    // ── Hertzian contact ──────────────────────────────────────────────────

    #[test]
    fn test_hertzian_contact_default_params() {
        let params = HertzianContactParams::default();
        let result = hertzian_contact(&params, 50_000.0);
        assert!(result.a > 0.0, "contact semi-axis a must be positive");
        assert!(result.b > 0.0, "contact semi-axis b must be positive");
        assert!(result.p0 > 0.0, "peak pressure must be positive");
        assert!(result.area > 0.0, "contact area must be positive");
    }

    #[test]
    fn test_hertzian_contact_zero_force() {
        let params = HertzianContactParams::default();
        let result = hertzian_contact(&params, 0.0);
        assert!(result.a.abs() < TOL, "zero force => zero contact radius");
        assert!(result.p0.abs() < TOL, "zero force => zero pressure");
    }

    #[test]
    fn test_hertzian_contact_pressure_scales_with_force() {
        let params = HertzianContactParams::default();
        let r1 = hertzian_contact(&params, 50_000.0);
        let r2 = hertzian_contact(&params, 100_000.0);
        assert!(
            r2.p0 > r1.p0,
            "doubling force should increase peak pressure"
        );
        assert!(r2.a > r1.a, "doubling force should increase contact size");
    }

    // ── Creepage / Kalker ─────────────────────────────────────────────────

    #[test]
    fn test_kalker_zero_creepage_zero_force() {
        let creepage = Creepage::default();
        let result = kalker_creep_forces(&creepage, 0.005, 0.004, 80.0e9);
        assert!(result.fx.abs() < TOL);
        assert!(result.fy.abs() < TOL);
    }

    #[test]
    fn test_kalker_longitudinal_creepage_produces_fx() {
        let creepage = Creepage {
            xi: 0.001,
            eta: 0.0,
            phi: 0.0,
        };
        let result = kalker_creep_forces(&creepage, 0.005, 0.004, 80.0e9);
        assert!(
            result.fx.abs() > 0.0,
            "longitudinal creepage should produce fx"
        );
        assert!(
            result.fx < 0.0,
            "positive creepage should produce negative (retarding) force"
        );
    }

    #[test]
    fn test_fastsim_saturates_at_friction_limit() {
        let creepage = Creepage {
            xi: 0.1,
            eta: 0.0,
            phi: 0.0,
        };
        let normal_force = 50_000.0;
        let mu = 0.3;
        let result = fastsim_creep_forces(&creepage, 0.005, 0.004, 80.0e9, normal_force, mu);
        let total_f = (result.fx * result.fx + result.fy * result.fy).sqrt();
        let limit = mu * normal_force;
        assert!(
            total_f <= limit + 1.0,
            "saturated force {} should not exceed mu*N = {}",
            total_f,
            limit
        );
    }

    // ── Track geometry ────────────────────────────────────────────────────

    #[test]
    fn test_equilibrium_cant_straight_track() {
        let cant = equilibrium_cant(30.0, f64::INFINITY, STANDARD_GAUGE);
        // Infinite radius => zero cant needed
        // Actually with f64::INFINITY the division gives ~0
        assert!(cant.abs() < 1e-6, "straight track needs no cant");
    }

    #[test]
    fn test_equilibrium_cant_curve() {
        let cant = equilibrium_cant(30.0, 500.0, STANDARD_GAUGE);
        assert!(cant > 0.0, "cant should be positive for a curve");
        // Expected: 1.435 * 900 / (9.81 * 500) ≈ 0.263 m
        let expected = STANDARD_GAUGE * 30.0 * 30.0 / (G * 500.0);
        assert!(
            (cant - expected).abs() < TOL,
            "cant={cant}, expected={expected}"
        );
    }

    #[test]
    fn test_cant_deficiency() {
        let cd = cant_deficiency(40.0, 300.0, 0.10, STANDARD_GAUGE);
        assert!(cd > 0.0, "at high speed there should be cant deficiency");
    }

    #[test]
    fn test_clothoid_curvature_linear() {
        let target = 1.0 / 500.0; // curvature for R=500m
        let mid = clothoid_curvature(50.0, 100.0, target);
        assert!(
            (mid - 0.5 * target).abs() < TOL,
            "curvature at midpoint should be half of target"
        );
    }

    #[test]
    fn test_clothoid_offset_straight() {
        let [x, y] = clothoid_offset(100.0, 1e15);
        // Very large A => essentially straight
        assert!((x - 100.0).abs() < 1.0, "x should be ~100 for straight");
        assert!(y.abs() < 1.0, "y should be ~0 for straight");
    }

    // ── Hunting oscillation ───────────────────────────────────────────────

    #[test]
    fn test_klingel_frequency_increases_with_speed() {
        let f1 = klingel_frequency(10.0, 0.05, 0.46, STANDARD_GAUGE / 2.0);
        let f2 = klingel_frequency(20.0, 0.05, 0.46, STANDARD_GAUGE / 2.0);
        assert!(f2 > f1, "frequency should increase with speed");
        assert!(
            (f2 / f1 - 2.0).abs() < TOL,
            "frequency should be proportional to speed"
        );
    }

    #[test]
    fn test_klingel_wavelength_positive() {
        let wl = klingel_wavelength(0.05, 0.46, STANDARD_GAUGE / 2.0);
        assert!(wl > 0.0, "wavelength must be positive");
        assert!(wl < 100.0, "wavelength should be realistic (< 100m)");
    }

    #[test]
    fn test_critical_hunting_speed_positive() {
        let params = HuntingParams::default();
        let v_crit = critical_hunting_speed(&params);
        assert!(v_crit > 0.0, "critical speed must be positive");
        // Should be in a realistic range for passenger trains
        assert!(
            v_crit > 0.5 && v_crit < 500.0,
            "v_crit={} should be positive and finite",
            v_crit
        );
    }

    #[test]
    fn test_hunting_stability_below_critical() {
        let params = HuntingParams::default();
        let v_crit = critical_hunting_speed(&params);
        assert!(!is_hunting_unstable(v_crit * 0.5, &params));
        assert!(is_hunting_unstable(v_crit * 1.5, &params));
    }

    // ── Ride comfort ──────────────────────────────────────────────────────

    #[test]
    fn test_sperling_wz_zero_for_no_vibration() {
        let accel = vec![0.0; 100];
        let wz = sperling_wz(&accel, 100.0);
        assert!(wz.abs() < TOL, "no vibration => Wz=0");
    }

    #[test]
    fn test_sperling_wz_increases_with_amplitude() {
        let mild: Vec<f64> = (0..200)
            .map(|i| 0.1 * (2.0 * PI * 5.0 * i as f64 / 200.0).sin())
            .collect();
        let rough: Vec<f64> = (0..200)
            .map(|i| 1.0 * (2.0 * PI * 5.0 * i as f64 / 200.0).sin())
            .collect();
        let wz_mild = sperling_wz(&mild, 200.0);
        let wz_rough = sperling_wz(&rough, 200.0);
        assert!(
            wz_rough > wz_mild,
            "rougher ride should have higher Wz: {} vs {}",
            wz_rough,
            wz_mild
        );
    }

    #[test]
    fn test_iso2631_rms() {
        let accel = vec![1.0, -1.0, 1.0, -1.0];
        let rms = iso2631_rms(&accel);
        assert!((rms - 1.0).abs() < TOL, "RMS of ±1 should be 1.0");
    }

    // ── Braking ───────────────────────────────────────────────────────────

    #[test]
    fn test_braking_distance_flat() {
        let d = braking_distance(100_000.0, 30.0, 0.0, 50_000.0, 0.0, 0.0);
        // Energy = 0.5 * 100000 * 900 = 45e6, d = 45e6 / 50000 = 900 m
        assert!((d - 900.0).abs() < 1.0, "braking distance on flat = {d}");
    }

    #[test]
    fn test_braking_distance_downhill_longer() {
        let d_flat = braking_distance(100_000.0, 30.0, 0.0, 50_000.0, 0.0, 0.0);
        let d_down = braking_distance(100_000.0, 30.0, 0.0, 50_000.0, -0.01, 0.0);
        assert!(
            d_down > d_flat,
            "downhill braking should be longer: {} vs {}",
            d_down,
            d_flat
        );
    }

    #[test]
    fn test_emergency_braking_includes_reaction() {
        let d_no_react = emergency_braking_distance(100_000.0, 30.0, 50_000.0, 0.0, 0.0);
        let d_react = emergency_braking_distance(100_000.0, 30.0, 50_000.0, 3.0, 0.0);
        assert!(
            d_react > d_no_react,
            "reaction time adds to braking distance"
        );
        assert!(
            (d_react - d_no_react - 90.0).abs() < 1.0,
            "reaction adds 30*3=90m"
        );
    }

    #[test]
    fn test_braking_time() {
        let t = braking_time(100_000.0, 30.0, 50_000.0, 0.0);
        // decel = 50000/100000 = 0.5 m/s², time = 30/0.5 = 60 s
        assert!((t - 60.0).abs() < 0.1, "braking time = {t}");
    }

    // ── Traction ──────────────────────────────────────────────────────────

    #[test]
    fn test_traction_motor_below_transition() {
        let motor = TractionMotor::new(1_000_000.0, 200_000.0);
        let effort = motor.tractive_effort(1.0);
        assert!(
            (effort - 200_000.0).abs() < TOL,
            "below transition speed effort should be max"
        );
    }

    #[test]
    fn test_traction_motor_above_transition() {
        let motor = TractionMotor::new(1_000_000.0, 200_000.0);
        // transition speed = 1e6 / 2e5 = 5 m/s
        let effort = motor.tractive_effort(10.0);
        // power / speed = 1e6 / 10 = 100000
        assert!(
            (effort - 100_000.0).abs() < 1.0,
            "above transition effort should be P/v"
        );
    }

    #[test]
    fn test_curtius_kniffler_decreases_with_speed() {
        let mu_low = curtius_kniffler_adhesion(0.0);
        let mu_high = curtius_kniffler_adhesion(200.0);
        assert!(
            mu_low > mu_high,
            "adhesion should decrease with speed: {} vs {}",
            mu_low,
            mu_high
        );
    }

    #[test]
    fn test_traction_curve_length() {
        let motor = TractionMotor::new(1_000_000.0, 200_000.0);
        let curve = traction_curve(&motor, 100_000.0, 4, 50.0, 20);
        assert_eq!(curve.len(), 20, "should have 20 points");
        assert!((curve[0].0).abs() < TOL, "first point should be at speed 0");
        assert!(
            (curve[19].0 - 50.0).abs() < TOL,
            "last point should be at max speed"
        );
    }

    // ── Coupler ───────────────────────────────────────────────────────────

    #[test]
    fn test_coupler_force_within_slack() {
        let params = CouplerParams::default();
        let f = coupler_force(0.005, 0.0, &params); // within 10mm slack
        assert!(f.abs() < TOL, "within slack => zero spring force");
    }

    #[test]
    fn test_coupler_force_tension() {
        let params = CouplerParams::default();
        let f = coupler_force(0.05, 0.0, &params); // 50mm extension
        assert!(f > 0.0, "extension should produce tension: {f}");
    }

    #[test]
    fn test_coupler_force_compression() {
        let params = CouplerParams::default();
        let f = coupler_force(-0.05, 0.0, &params); // 50mm compression
        assert!(f < 0.0, "compression should produce negative force: {f}");
    }

    #[test]
    fn test_train_coupler_forces_count() {
        let positions = vec![0.0, 25.0, 50.0];
        let velocities = vec![10.0, 10.0, 10.0];
        let lengths = vec![25.0, 25.0, 25.0];
        let params = CouplerParams::default();
        let forces = train_coupler_forces(&positions, &velocities, &lengths, &params);
        assert_eq!(forces.len(), 2, "3 cars => 2 couplers");
    }

    // ── Derailment criteria ───────────────────────────────────────────────

    #[test]
    fn test_nadal_lv_ratio() {
        let lv = nadal_lv_ratio(30_000.0, 50_000.0);
        assert!((lv - 0.6).abs() < TOL, "L/V = {lv}");
    }

    #[test]
    fn test_nadal_critical_lv_typical() {
        // 70 degree flange angle, mu=0.3
        let critical = nadal_critical_lv(70.0_f64.to_radians(), 0.3);
        assert!(
            critical > 0.5 && critical < 2.0,
            "critical L/V should be realistic: {}",
            critical
        );
    }

    #[test]
    fn test_derailment_check() {
        let flange = 70.0_f64.to_radians();
        let mu = 0.3;
        let critical = nadal_critical_lv(flange, mu);
        // Safe case
        assert!(!check_nadal_derailment(10_000.0, 50_000.0, flange, mu));
        // Dangerous case: L/V = 2.0 > critical
        let high_lateral = critical * 50_000.0 * 1.5;
        assert!(check_nadal_derailment(high_lateral, 50_000.0, flange, mu));
    }

    #[test]
    fn test_wheel_unloading_ratio() {
        let ratio = wheel_unloading_ratio(40_000.0, 50_000.0);
        assert!((ratio - 0.8).abs() < TOL);
    }

    #[test]
    fn test_wheel_unloading_danger() {
        assert!(check_wheel_unloading(2_000.0, 50_000.0, 0.1));
        assert!(!check_wheel_unloading(40_000.0, 50_000.0, 0.1));
    }

    #[test]
    fn test_prudhomme_limit() {
        let limit = prudhomme_limit(150_000.0);
        // 10000 + 150000/3 = 60000 N
        assert!((limit - 60_000.0).abs() < TOL);
    }

    // ── Speed conversion ──────────────────────────────────────────────────

    #[test]
    fn test_speed_conversion_roundtrip() {
        let v = 100.0;
        assert!((ms_to_kmh(kmh_to_ms(v)) - v).abs() < TOL);
    }

    // ── Grade helpers ─────────────────────────────────────────────────────

    #[test]
    fn test_permille_conversion() {
        assert!((permille_to_ratio(25.0) - 0.025).abs() < TOL);
        assert!((ratio_to_permille(0.025) - 25.0).abs() < TOL);
    }

    #[test]
    fn test_grade_resistance() {
        let r = grade_resistance(100_000.0, 0.01);
        assert!((r - 100_000.0 * G * 0.01).abs() < 1.0);
    }

    // ── Curve negotiation ─────────────────────────────────────────────────

    #[test]
    fn test_max_curve_speed() {
        let v = max_curve_speed(500.0, 0.10, STANDARD_GAUGE, 0.15);
        assert!(v > 20.0 && v < 100.0, "max curve speed = {v} m/s");
    }

    #[test]
    fn test_centrifugal_force() {
        let f = centrifugal_force(50_000.0, 30.0, 500.0);
        // m*v²/R = 50000*900/500 = 90000
        assert!((f - 90_000.0).abs() < 1.0);
    }

    // ── Gravitational stiffness ───────────────────────────────────────────

    #[test]
    fn test_gravitational_stiffness() {
        let k = gravitational_stiffness(0.05, 100_000.0, 0.46);
        assert!(k > 0.0, "stiffness must be positive");
        let expected = 2.0 * 0.05 * 100_000.0 / 0.46;
        assert!((k - expected).abs() < 1.0);
    }

    // ── Wheel radii ───────────────────────────────────────────────────────

    #[test]
    fn test_wheel_radii_centered() {
        let (rl, rr) = wheel_radii(0.46, 0.05, 0.0);
        assert!((rl - 0.46).abs() < TOL);
        assert!((rr - 0.46).abs() < TOL);
    }

    #[test]
    fn test_wheel_radii_displaced() {
        let (rl, rr) = wheel_radii(0.46, 0.05, 0.005);
        assert!(
            rl > rr,
            "left wheel should be larger for positive displacement"
        );
    }

    #[test]
    fn test_effective_conicity() {
        let lambda = effective_conicity(0.002, 0.02);
        assert!((lambda - 0.05).abs() < TOL);
    }

    // ── Train simulation step ─────────────────────────────────────────────

    #[test]
    fn test_train_longitudinal_step_single_car() {
        let mut cars = vec![TrainCarState {
            mass: 50_000.0,
            position: 0.0,
            velocity: 0.0,
            length: 25.0,
        }];
        let traction = vec![100_000.0];
        let params = CouplerParams::default();
        train_longitudinal_step(&mut cars, &traction, 0.0, &params, 0.1);
        assert!(cars[0].velocity > 0.0, "car should accelerate");
        assert!(cars[0].position > 0.0, "car should move forward");
    }

    // ── Davis resistance ──────────────────────────────────────────────────

    #[test]
    fn test_davis_resistance_increases_with_speed() {
        let r1 = davis_resistance_force(100_000.0, 10.0, 500.0, 5.0, 0.6);
        let r2 = davis_resistance_force(100_000.0, 30.0, 500.0, 5.0, 0.6);
        assert!(r2 > r1, "resistance should increase with speed");
    }

    // ── Creep yaw torque ──────────────────────────────────────────────────

    #[test]
    fn test_creep_yaw_torque_zero_displacement() {
        let torque = creep_yaw_torque(10.0e6, 0.0, 30.0, 0.46, STANDARD_GAUGE / 2.0);
        assert!(torque.abs() < TOL, "no displacement => no yaw torque");
    }

    #[test]
    fn test_creep_yaw_torque_nonzero() {
        let delta_r = rolling_radius_difference(0.05, 0.005);
        let torque = creep_yaw_torque(10.0e6, delta_r, 30.0, 0.46, STANDARD_GAUGE / 2.0);
        assert!(
            torque.abs() > 0.0,
            "displaced wheelset should have yaw torque"
        );
    }
}
