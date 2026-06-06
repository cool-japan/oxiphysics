// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Advanced suspension analysis: kinematic geometry, damper models, comfort
//! analysis, and NVH natural-frequency calculation.

use std::f64::consts::PI;

// ---------------------------------------------------------------------------
// Helper vector math (no nalgebra)
// ---------------------------------------------------------------------------

#[cfg(test)]
fn vec3_add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

#[cfg(test)]
fn vec3_sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

#[cfg(test)]
fn vec3_scale(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

#[cfg(test)]
fn vec3_dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[cfg(test)]
fn vec3_norm(a: [f64; 3]) -> f64 {
    vec3_dot(a, a).sqrt()
}

#[cfg(test)]
fn vec3_cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

// ---------------------------------------------------------------------------
// SuspensionKinematic
// ---------------------------------------------------------------------------

/// Basic suspension kinematic parameters shared across all geometry types.
#[derive(Debug, Clone)]
pub struct SuspensionKinematic {
    /// Track width (m) — distance between left and right wheel centres.
    pub track_width: f64,
    /// Wheelbase (m) — distance between front and rear axle centrelines.
    pub wheelbase: f64,
    /// Roll centre height above ground (m).
    pub roll_centre_height: f64,
    /// Scrub radius (m) — lateral offset of kingpin axis at ground plane.
    pub scrub_radius: f64,
    /// Caster angle (rad) — longitudinal tilt of the kingpin axis.
    pub caster_angle: f64,
    /// Kingpin inclination angle (rad) — lateral tilt of the kingpin axis.
    pub kingpin_inclination: f64,
    /// Static camber angle (rad) — positive = top leaning outward.
    pub static_camber: f64,
    /// Static toe angle (rad) — positive = toe-in.
    pub static_toe: f64,
}

impl SuspensionKinematic {
    /// Create a default front-axle kinematic configuration.
    pub fn default_front() -> Self {
        Self {
            track_width: 1.60,
            wheelbase: 2.70,
            roll_centre_height: 0.07,
            scrub_radius: 0.025,
            caster_angle: 4.0_f64.to_radians(),
            kingpin_inclination: 12.0_f64.to_radians(),
            static_camber: (-1.5_f64).to_radians(),
            static_toe: 0.5_f64.to_radians(),
        }
    }

    /// Create a default rear-axle kinematic configuration.
    pub fn default_rear() -> Self {
        Self {
            track_width: 1.55,
            wheelbase: 2.70,
            roll_centre_height: 0.12,
            scrub_radius: 0.0,
            caster_angle: 0.0,
            kingpin_inclination: 0.0,
            static_camber: (-2.0_f64).to_radians(),
            static_toe: (-0.2_f64).to_radians(),
        }
    }

    /// Compute mechanical trail from caster angle and kingpin offset height.
    ///
    /// `h_kingpin` — height of the kingpin pivot above ground (m).
    pub fn mechanical_trail(&self, h_kingpin: f64) -> f64 {
        h_kingpin * self.caster_angle.tan()
    }

    /// Compute the pneumatic trail contribution to steering returnability.
    ///
    /// `contact_patch_length` — longitudinal half-length of the contact patch (m).
    pub fn pneumatic_trail(contact_patch_length: f64) -> f64 {
        contact_patch_length / 3.0
    }

    /// Total steering offset at ground = scrub radius + kingpin offset lateral.
    pub fn total_steering_offset(&self) -> f64 {
        self.scrub_radius
    }
}

// ---------------------------------------------------------------------------
// MacPhersonStrut
// ---------------------------------------------------------------------------

/// MacPherson strut suspension geometry.
#[derive(Debug, Clone)]
pub struct MacPhersonStrut {
    /// Strut length at design ride height (m).
    pub strut_length: f64,
    /// Height of the spring perch above the wheel centre (m).
    pub spring_perch_height: f64,
    /// Steer axis inclination angle (rad).
    pub steer_axis_inclination: f64,
    /// Static camber at design ride height (rad).
    pub static_camber: f64,
    /// Camber change rate (rad/m of bump travel).
    pub camber_rate: f64,
}

impl MacPhersonStrut {
    /// Default MacPherson strut representative of a compact front suspension.
    pub fn default_compact() -> Self {
        Self {
            strut_length: 0.42,
            spring_perch_height: 0.18,
            steer_axis_inclination: 13.0_f64.to_radians(),
            static_camber: (-1.0_f64).to_radians(),
            camber_rate: 0.8_f64.to_radians(),
        }
    }

    /// Compute camber at a given bump travel (m). Positive bump = compression.
    pub fn camber_at_bump(&self, bump: f64) -> f64 {
        self.static_camber + self.camber_rate * bump
    }

    /// Estimate the instant-centre height for this strut geometry.
    ///
    /// Uses the half-track width `half_track` (m) and the lower balljoint
    /// inboard pivot height `pivot_height` (m).
    pub fn instant_centre_height(&self, half_track: f64, pivot_height: f64) -> f64 {
        // Simplified: IC at intersection of strut axis and imaginary lower arm
        let strut_angle = self.static_camber.abs() + self.steer_axis_inclination;
        pivot_height + half_track * strut_angle.tan()
    }
}

// ---------------------------------------------------------------------------
// DoubleWishbone
// ---------------------------------------------------------------------------

/// Double wishbone suspension geometry parameters.
#[derive(Debug, Clone)]
pub struct DoubleWishbone {
    /// Upper wishbone length (m).
    pub upper_arm_length: f64,
    /// Lower wishbone length (m).
    pub lower_arm_length: f64,
    /// Upper arm angle relative to horizontal (rad).
    pub upper_arm_angle: f64,
    /// Lower arm angle relative to horizontal (rad).
    pub lower_arm_angle: f64,
    /// Vertical distance between upper and lower ball joints (m).
    pub ball_joint_span: f64,
    /// Static camber (rad).
    pub static_camber: f64,
}

impl DoubleWishbone {
    /// Default double-wishbone representative of a sports car front axle.
    pub fn default_sports() -> Self {
        Self {
            upper_arm_length: 0.32,
            lower_arm_length: 0.40,
            upper_arm_angle: 5.0_f64.to_radians(),
            lower_arm_angle: 3.0_f64.to_radians(),
            ball_joint_span: 0.30,
            static_camber: (-2.0_f64).to_radians(),
        }
    }

    /// Approximate camber gain per metre of bump travel (rad/m).
    ///
    /// Uses the difference in arm lengths: shorter upper arm → negative camber gain.
    pub fn camber_gain(&self) -> f64 {
        // Simplified linear approximation
        let arm_ratio = self.upper_arm_length / self.lower_arm_length;
        -(1.0 - arm_ratio) / self.ball_joint_span
    }

    /// Compute roll-centre height using the Force-Line method (simplified).
    ///
    /// `half_track` — half the track width (m).
    pub fn roll_centre_height(&self, half_track: f64) -> f64 {
        // IC height at intersection of upper/lower arm extension lines
        let ic_height = half_track * (self.lower_arm_angle.tan() - self.upper_arm_angle.tan())
            / (1.0 + self.lower_arm_length / self.upper_arm_length);
        ic_height.max(0.0)
    }

    /// Compute camber at a given bump travel (m).
    pub fn camber_at_bump(&self, bump: f64) -> f64 {
        self.static_camber + self.camber_gain() * bump
    }

    /// Instant-centre distance from the wheel (m) in the front view.
    pub fn instant_centre_distance(&self) -> f64 {
        let arm_diff = (self.lower_arm_length - self.upper_arm_length).abs();
        self.ball_joint_span / (arm_diff + 1e-6) * self.upper_arm_length
    }
}

// ---------------------------------------------------------------------------
// MultiLink
// ---------------------------------------------------------------------------

/// 5-link rear suspension geometry.
#[derive(Debug, Clone)]
pub struct MultiLink {
    /// Static camber (rad).
    pub static_camber: f64,
    /// Static toe (rad) — positive = toe-in.
    pub static_toe: f64,
    /// Camber gradient with bump travel (rad/m).
    pub camber_gradient: f64,
    /// Toe gradient with bump travel (rad/m).
    pub toe_gradient: f64,
    /// Anti-squat geometry factor (fraction, 0–1+).
    pub anti_squat: f64,
    /// Anti-lift geometry factor (fraction, 0–1+).
    pub anti_lift: f64,
    /// Roll-centre height (m).
    pub roll_centre_height: f64,
}

impl MultiLink {
    /// Default 5-link rear geometry for a sports saloon.
    pub fn default_sports_rear() -> Self {
        Self {
            static_camber: (-2.5_f64).to_radians(),
            static_toe: (-0.3_f64).to_radians(),
            camber_gradient: 0.6_f64.to_radians(),
            toe_gradient: 0.15_f64.to_radians(),
            anti_squat: 0.65,
            anti_lift: 0.40,
            roll_centre_height: 0.13,
        }
    }

    /// Camber angle at a given bump travel (m).
    pub fn camber_at_bump(&self, bump: f64) -> f64 {
        self.static_camber + self.camber_gradient * bump
    }

    /// Toe angle at a given bump travel (m).
    pub fn toe_at_bump(&self, bump: f64) -> f64 {
        self.static_toe + self.toe_gradient * bump
    }

    /// Anti-squat force transfer ratio under acceleration.
    ///
    /// `cog_height` — centre of gravity height (m).
    /// `wheelbase` — vehicle wheelbase (m).
    pub fn anti_squat_ratio(&self, cog_height: f64, wheelbase: f64) -> f64 {
        self.anti_squat * wheelbase / cog_height.max(1e-3)
    }
}

// ---------------------------------------------------------------------------
// AntiRollBar
// ---------------------------------------------------------------------------

/// Anti-roll bar (stabiliser bar) model.
#[derive(Debug, Clone)]
pub struct AntiRollBar {
    /// Torsional stiffness of the bar (N·m/rad).
    pub stiffness: f64,
    /// Motion ratio (suspension travel to bar twist, dimensionless).
    pub motion_ratio: f64,
    /// Effective roll stiffness contribution (N·m/rad of body roll).
    pub effective_roll_stiffness: f64,
}

impl AntiRollBar {
    /// Construct from bar stiffness and motion ratio.
    pub fn new(stiffness: f64, motion_ratio: f64) -> Self {
        let effective = stiffness * motion_ratio * motion_ratio;
        Self {
            stiffness,
            motion_ratio,
            effective_roll_stiffness: effective,
        }
    }

    /// Front/rear roll stiffness balance factor.
    ///
    /// Returns the fraction of total roll stiffness at this axle
    /// given the other axle's effective roll stiffness.
    pub fn roll_balance(&self, other_effective: f64) -> f64 {
        let total = self.effective_roll_stiffness + other_effective;
        if total < 1e-6 {
            0.5
        } else {
            self.effective_roll_stiffness / total
        }
    }

    /// Roll gradient (deg/g) from spring stiffnesses and geometry.
    ///
    /// `spring_stiffness` — combined axle spring stiffness (N/m).
    /// `half_track` — half track width (m).
    /// `mass` — sprung mass (kg).
    pub fn roll_gradient(&self, spring_stiffness: f64, half_track: f64, mass: f64) -> f64 {
        let g = 9.81;
        let spring_moment = 2.0 * spring_stiffness * half_track * half_track;
        let roll_stiffness = spring_moment + self.effective_roll_stiffness;
        if roll_stiffness < 1e-6 {
            return f64::INFINITY;
        }
        (mass * g * half_track / roll_stiffness).to_degrees()
    }
}

// ---------------------------------------------------------------------------
// BumpStop
// ---------------------------------------------------------------------------

/// Progressive bump-stop model with rebound buffer.
#[derive(Debug, Clone)]
pub struct BumpStop {
    /// Travel at which the bump stop begins to engage (m).
    pub engagement_travel: f64,
    /// Progressive stiffness coefficient (N/m²) — F = k * x².
    pub progressive_coeff: f64,
    /// Rebound buffer travel limit (m) — negative from design ride height.
    pub rebound_limit: f64,
    /// Linear bump travel limit (m).
    pub bump_limit: f64,
}

impl BumpStop {
    /// Default bump-stop parameters for a road car.
    pub fn default_road() -> Self {
        Self {
            engagement_travel: 0.040,
            progressive_coeff: 80_000.0,
            rebound_limit: -0.060,
            bump_limit: 0.075,
        }
    }

    /// Compute bump-stop force at `travel` (m) beyond design ride height.
    ///
    /// Returns zero below the engagement point, progressive beyond it.
    pub fn force(&self, travel: f64) -> f64 {
        if travel <= self.engagement_travel {
            0.0
        } else {
            let excess = travel - self.engagement_travel;
            self.progressive_coeff * excess * excess
        }
    }

    /// Returns `true` if travel is outside the legal range.
    pub fn is_limit_hit(&self, travel: f64) -> bool {
        travel >= self.bump_limit || travel <= self.rebound_limit
    }
}

// ---------------------------------------------------------------------------
// SuspensionDamper
// ---------------------------------------------------------------------------

/// Valved twin-tube damper model with velocity-dependent bump/rebound asymmetry.
#[derive(Debug, Clone)]
pub struct SuspensionDamper {
    /// Low-speed bump (compression) damping coefficient (N·s/m).
    pub bump_low_speed: f64,
    /// High-speed bump damping coefficient (N·s/m).
    pub bump_high_speed: f64,
    /// Low-speed rebound damping coefficient (N·s/m).
    pub rebound_low_speed: f64,
    /// High-speed rebound damping coefficient (N·s/m).
    pub rebound_high_speed: f64,
    /// Velocity threshold separating low/high speed (m/s).
    pub speed_threshold: f64,
    /// Dry-friction / stiction force (N).
    pub stiction: f64,
}

impl SuspensionDamper {
    /// Default damper parameters for a sports road car.
    pub fn default_sports() -> Self {
        Self {
            bump_low_speed: 3_500.0,
            bump_high_speed: 1_200.0,
            rebound_low_speed: 5_000.0,
            rebound_high_speed: 1_800.0,
            speed_threshold: 0.10,
            stiction: 25.0,
        }
    }

    /// Compute the damping force (N) for a given piston velocity (m/s).
    ///
    /// Positive velocity = bump (compression), negative = rebound (extension).
    pub fn force(&self, velocity: f64) -> f64 {
        let (c_low, c_high) = if velocity >= 0.0 {
            (self.bump_low_speed, self.bump_high_speed)
        } else {
            (self.rebound_low_speed, self.rebound_high_speed)
        };
        let speed = velocity.abs();
        let coeff = if speed <= self.speed_threshold {
            c_low
        } else {
            // Blend at threshold, then high-speed slope
            let excess = speed - self.speed_threshold;
            (c_low * self.speed_threshold + c_high * excess) / speed
        };
        let raw_force = coeff * speed;
        let sign = if velocity >= 0.0 { 1.0 } else { -1.0 };
        sign * (raw_force + self.stiction)
    }

    /// Bump-to-rebound damping ratio at low speed.
    pub fn asymmetry_ratio(&self) -> f64 {
        self.rebound_low_speed / self.bump_low_speed.max(1.0)
    }
}

// ---------------------------------------------------------------------------
// KinematicSolver
// ---------------------------------------------------------------------------

/// Result from the suspension kinematic solver at a given bump travel.
#[derive(Debug, Clone)]
pub struct KinematicResult {
    /// Camber angle (rad).
    pub camber: f64,
    /// Toe angle (rad).
    pub toe: f64,
    /// Effective caster angle (rad).
    pub caster: f64,
    /// Kingpin offset at the wheel centre plane (m).
    pub kingpin_offset: f64,
    /// Scrub radius (m).
    pub scrub: f64,
}

/// Iterative kinematic solver over bump travel for a given suspension geometry.
#[derive(Debug, Clone)]
pub struct KinematicSolver {
    /// Base kinematic parameters.
    pub kinematics: SuspensionKinematic,
    /// Optional double-wishbone geometry.
    pub double_wishbone: Option<DoubleWishbone>,
    /// Optional multi-link geometry.
    pub multi_link: Option<MultiLink>,
}

impl KinematicSolver {
    /// Construct from base kinematic parameters only.
    pub fn new(kinematics: SuspensionKinematic) -> Self {
        Self {
            kinematics,
            double_wishbone: None,
            multi_link: None,
        }
    }

    /// Attach double-wishbone geometry for more detailed camber/toe computation.
    pub fn with_double_wishbone(mut self, dw: DoubleWishbone) -> Self {
        self.double_wishbone = Some(dw);
        self
    }

    /// Attach multi-link geometry.
    pub fn with_multi_link(mut self, ml: MultiLink) -> Self {
        self.multi_link = Some(ml);
        self
    }

    /// Iterate geometry at `travel` (m) — positive = bump (compression).
    pub fn iterate_geometry(&self, travel: f64) -> KinematicResult {
        let base_camber = self.kinematics.static_camber;
        let base_toe = self.kinematics.static_toe;
        let base_caster = self.kinematics.caster_angle;

        let (camber, toe) = if let Some(ref dw) = self.double_wishbone {
            (dw.camber_at_bump(travel), base_toe)
        } else if let Some(ref ml) = self.multi_link {
            (ml.camber_at_bump(travel), ml.toe_at_bump(travel))
        } else {
            (base_camber, base_toe)
        };

        KinematicResult {
            camber,
            toe,
            caster: base_caster,
            kingpin_offset: self.kinematics.kingpin_inclination * 0.10, // simplified
            scrub: self.kinematics.scrub_radius,
        }
    }

    /// Generate a vector of results over `n` steps of bump from `min` to `max` (m).
    pub fn sweep_travel(
        &self,
        min_travel: f64,
        max_travel: f64,
        n: usize,
    ) -> Vec<(f64, KinematicResult)> {
        if n < 2 {
            return vec![(min_travel, self.iterate_geometry(min_travel))];
        }
        let step = (max_travel - min_travel) / (n - 1) as f64;
        (0..n)
            .map(|i| {
                let t = min_travel + i as f64 * step;
                (t, self.iterate_geometry(t))
            })
            .collect()
    }
}

// ---------------------------------------------------------------------------
// ComfortAnalysis (Quarter-car model)
// ---------------------------------------------------------------------------

/// ISO 8608 road roughness class.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RoadClass {
    /// Class A — very smooth motorway.
    A,
    /// Class B — good quality road.
    B,
    /// Class C — average road.
    C,
    /// Class D — poor road.
    D,
    /// Class E — very poor / off road.
    E,
}

impl RoadClass {
    /// Roughness coefficient Gd (m² / (rad/m)) at n0 = 0.1 cycles/m.
    pub fn roughness_coeff(self) -> f64 {
        match self {
            RoadClass::A => 16e-6,
            RoadClass::B => 64e-6,
            RoadClass::C => 256e-6,
            RoadClass::D => 1024e-6,
            RoadClass::E => 4096e-6,
        }
    }
}

/// Quarter-car comfort analysis model.
///
/// Masses: sprung (body) and unsprung (wheel + hub).
#[derive(Debug, Clone)]
pub struct ComfortAnalysis {
    /// Sprung mass per corner (kg).
    pub sprung_mass: f64,
    /// Unsprung mass per corner (kg).
    pub unsprung_mass: f64,
    /// Suspension spring stiffness (N/m).
    pub spring_stiffness: f64,
    /// Suspension damping coefficient (N·s/m).
    pub damping: f64,
    /// Tyre vertical stiffness (N/m).
    pub tyre_stiffness: f64,
    /// Road class for excitation spectrum.
    pub road_class: RoadClass,
}

impl ComfortAnalysis {
    /// Default quarter-car comfort parameters for a road saloon.
    pub fn default_saloon() -> Self {
        Self {
            sprung_mass: 375.0,
            unsprung_mass: 45.0,
            spring_stiffness: 22_000.0,
            damping: 2_200.0,
            tyre_stiffness: 180_000.0,
            road_class: RoadClass::B,
        }
    }

    /// Natural frequency of the sprung mass (Hz).
    pub fn body_natural_freq(&self) -> f64 {
        let omega = (self.spring_stiffness / self.sprung_mass).sqrt();
        omega / (2.0 * PI)
    }

    /// Natural frequency of the unsprung mass (Hz).
    pub fn wheel_natural_freq(&self) -> f64 {
        let k_total = self.spring_stiffness + self.tyre_stiffness;
        let omega = (k_total / self.unsprung_mass).sqrt();
        omega / (2.0 * PI)
    }

    /// Damping ratio of the suspension.
    pub fn damping_ratio(&self) -> f64 {
        let cc = 2.0 * (self.spring_stiffness * self.sprung_mass).sqrt();
        self.damping / cc
    }

    /// Approximate RMS body acceleration (m/s²) at vehicle speed `v` (m/s)
    /// using ISO 8608 power spectral density.
    ///
    /// A simplified single-DOF approximation is used.
    pub fn rms_body_acceleration(&self, vehicle_speed: f64) -> f64 {
        let gd = self.road_class.roughness_coeff();
        let fn_hz = self.body_natural_freq();
        let zeta = self.damping_ratio();
        // PSD approximation: S_z(n) = Gd / n^2, integrated through resonance
        let n0 = fn_hz / vehicle_speed.max(1.0);
        let psd_at_resonance = gd / (n0 * n0).max(1e-12);
        let bandwidth = fn_hz * zeta;
        (psd_at_resonance * bandwidth * vehicle_speed).sqrt()
    }

    /// Wertman discomfort index (weighted RMS acceleration normalised to 0.315 m/s²).
    pub fn discomfort_index(&self, vehicle_speed: f64) -> f64 {
        let rms = self.rms_body_acceleration(vehicle_speed);
        rms / 0.315
    }
}

// ---------------------------------------------------------------------------
// NVH — Natural frequencies, mode shapes, transfer functions
// ---------------------------------------------------------------------------

/// Vehicle NVH (Noise, Vibration, Harshness) natural frequency set.
#[derive(Debug, Clone)]
pub struct Nvh {
    /// Sprung mass (kg).
    pub sprung_mass: f64,
    /// Pitch moment of inertia (kg·m²).
    pub pitch_inertia: f64,
    /// Roll moment of inertia (kg·m²).
    pub roll_inertia: f64,
    /// Front axle spring stiffness (N/m).
    pub front_stiffness: f64,
    /// Rear axle spring stiffness (N/m).
    pub rear_stiffness: f64,
    /// Half wheelbase (m).
    pub half_wheelbase: f64,
    /// Half track width (m).
    pub half_track: f64,
    /// Front damping coefficient (N·s/m).
    pub front_damping: f64,
    /// Rear damping coefficient (N·s/m).
    pub rear_damping: f64,
}

impl Nvh {
    /// Default NVH parameters for a road saloon.
    pub fn default_saloon() -> Self {
        Self {
            sprung_mass: 1400.0,
            pitch_inertia: 2_600.0,
            roll_inertia: 650.0,
            front_stiffness: 22_000.0,
            rear_stiffness: 24_000.0,
            half_wheelbase: 1.35,
            half_track: 0.80,
            front_damping: 2_200.0,
            rear_damping: 2_400.0,
        }
    }

    /// Bounce (heave) natural frequency (Hz).
    pub fn bounce_frequency(&self) -> f64 {
        let k_total = 2.0 * self.front_stiffness + 2.0 * self.rear_stiffness;
        (k_total / self.sprung_mass).sqrt() / (2.0 * PI)
    }

    /// Pitch natural frequency (Hz).
    pub fn pitch_frequency(&self) -> f64 {
        let k_pitch = 2.0 * self.front_stiffness * self.half_wheelbase.powi(2)
            + 2.0 * self.rear_stiffness * self.half_wheelbase.powi(2);
        (k_pitch / self.pitch_inertia).sqrt() / (2.0 * PI)
    }

    /// Roll natural frequency (Hz).
    pub fn roll_frequency(&self) -> f64 {
        let k_roll = 2.0 * (self.front_stiffness + self.rear_stiffness) * self.half_track.powi(2);
        (k_roll / self.roll_inertia).sqrt() / (2.0 * PI)
    }

    /// Bounce damping ratio.
    pub fn bounce_damping_ratio(&self) -> f64 {
        let k_total = 2.0 * self.front_stiffness + 2.0 * self.rear_stiffness;
        let c_total = 2.0 * self.front_damping + 2.0 * self.rear_damping;
        c_total / (2.0 * (k_total * self.sprung_mass).sqrt())
    }

    /// Simple transmissibility at frequency `f` (Hz) for the bounce mode.
    ///
    /// Returns the ratio of output to input acceleration amplitude.
    pub fn transmissibility(&self, freq: f64) -> f64 {
        let fn_ = self.bounce_frequency();
        let zeta = self.bounce_damping_ratio();
        let r = freq / fn_.max(1e-6);
        let num = (1.0 + (2.0 * zeta * r).powi(2)).sqrt();
        let den = ((1.0 - r * r).powi(2) + (2.0 * zeta * r).powi(2)).sqrt();
        num / den.max(1e-12)
    }

    /// Resonance peak transmissibility (approximate).
    pub fn peak_transmissibility(&self) -> f64 {
        let zeta = self.bounce_damping_ratio().max(1e-6);
        1.0 / (2.0 * zeta * (1.0 - zeta * zeta).sqrt().max(1e-6))
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // --- SuspensionKinematic ---

    #[test]
    fn mechanical_trail_positive() {
        let k = SuspensionKinematic::default_front();
        let trail = k.mechanical_trail(0.30);
        assert!(trail > 0.0, "trail should be positive: {trail}");
    }

    #[test]
    fn pneumatic_trail_one_third_patch() {
        let trail = SuspensionKinematic::pneumatic_trail(0.12);
        assert!((trail - 0.04).abs() < 1e-10);
    }

    #[test]
    fn scrub_radius_accessible() {
        let k = SuspensionKinematic::default_front();
        assert!(k.scrub_radius > 0.0);
    }

    #[test]
    fn default_rear_roll_centre_higher_than_front() {
        let front = SuspensionKinematic::default_front();
        let rear = SuspensionKinematic::default_rear();
        assert!(rear.roll_centre_height > front.roll_centre_height);
    }

    // --- MacPhersonStrut ---

    #[test]
    fn macpherson_camber_increases_with_bump() {
        let s = MacPhersonStrut::default_compact();
        let c0 = s.camber_at_bump(0.0);
        let c1 = s.camber_at_bump(0.05);
        assert!(
            c1 > c0,
            "camber should increase (become less negative) with bump"
        );
    }

    #[test]
    fn macpherson_instant_centre_height_positive() {
        let s = MacPhersonStrut::default_compact();
        let h = s.instant_centre_height(0.80, 0.25);
        assert!(h > 0.0, "IC height should be positive: {h}");
    }

    // --- DoubleWishbone ---

    #[test]
    fn double_wishbone_camber_gain_negative() {
        let dw = DoubleWishbone::default_sports();
        assert!(
            dw.camber_gain() < 0.0,
            "camber gain should be negative (more camber in bump)"
        );
    }

    #[test]
    fn double_wishbone_roll_centre_non_negative() {
        let dw = DoubleWishbone::default_sports();
        let rc = dw.roll_centre_height(0.80);
        assert!(rc >= 0.0);
    }

    #[test]
    fn double_wishbone_camber_at_bump() {
        let dw = DoubleWishbone::default_sports();
        let c = dw.camber_at_bump(0.03);
        // Should differ from static
        assert!((c - dw.static_camber).abs() > 1e-10);
    }

    #[test]
    fn double_wishbone_ic_distance_positive() {
        let dw = DoubleWishbone::default_sports();
        assert!(dw.instant_centre_distance() > 0.0);
    }

    // --- MultiLink ---

    #[test]
    fn multilink_camber_at_bump() {
        let ml = MultiLink::default_sports_rear();
        let c0 = ml.camber_at_bump(0.0);
        let c1 = ml.camber_at_bump(0.02);
        assert!((c1 - c0).abs() > 1e-10);
    }

    #[test]
    fn multilink_toe_at_bump() {
        let ml = MultiLink::default_sports_rear();
        let t0 = ml.toe_at_bump(0.0);
        let t1 = ml.toe_at_bump(0.02);
        assert!((t1 - t0).abs() > 1e-10);
    }

    #[test]
    fn multilink_anti_squat_ratio_reasonable() {
        let ml = MultiLink::default_sports_rear();
        let r = ml.anti_squat_ratio(0.50, 2.70);
        assert!(r > 0.0 && r < 10.0, "anti-squat ratio out of range: {r}");
    }

    // --- AntiRollBar ---

    #[test]
    fn anti_roll_bar_effective_stiffness() {
        let arb = AntiRollBar::new(10_000.0, 0.5);
        assert!((arb.effective_roll_stiffness - 2_500.0).abs() < 1e-6);
    }

    #[test]
    fn anti_roll_bar_balance_sums_to_one() {
        let front = AntiRollBar::new(10_000.0, 0.5);
        let rear_eff = 3_000.0;
        let balance = front.roll_balance(rear_eff);
        assert!(
            (balance + (rear_eff / (front.effective_roll_stiffness + rear_eff)) - 1.0).abs()
                < 1e-10
        );
    }

    #[test]
    fn anti_roll_bar_roll_gradient_positive() {
        let arb = AntiRollBar::new(10_000.0, 0.5);
        let rg = arb.roll_gradient(22_000.0, 0.80, 400.0);
        assert!(rg > 0.0, "roll gradient should be positive: {rg}");
    }

    // --- BumpStop ---

    #[test]
    fn bump_stop_zero_below_engagement() {
        let bs = BumpStop::default_road();
        assert_eq!(bs.force(0.010), 0.0);
        assert_eq!(bs.force(bs.engagement_travel), 0.0);
    }

    #[test]
    fn bump_stop_progressive_above_engagement() {
        let bs = BumpStop::default_road();
        let f1 = bs.force(bs.engagement_travel + 0.01);
        let f2 = bs.force(bs.engagement_travel + 0.02);
        assert!(f2 > f1, "force should increase progressively: {f1} < {f2}");
    }

    #[test]
    fn bump_stop_limit_detection() {
        let bs = BumpStop::default_road();
        assert!(bs.is_limit_hit(bs.bump_limit));
        assert!(bs.is_limit_hit(bs.rebound_limit));
        assert!(!bs.is_limit_hit(0.0));
    }

    // --- SuspensionDamper ---

    #[test]
    fn damper_zero_velocity_returns_stiction() {
        let d = SuspensionDamper::default_sports();
        let f = d.force(0.0);
        assert!(
            (f.abs() - d.stiction).abs() < 1e-6,
            "f={f}, stiction={}",
            d.stiction
        );
    }

    #[test]
    fn damper_rebound_greater_than_bump_at_low_speed() {
        let d = SuspensionDamper::default_sports();
        let f_bump = d.force(0.05).abs();
        let f_rebound = d.force(-0.05).abs();
        assert!(
            f_rebound > f_bump,
            "rebound should be stiffer: {f_rebound} vs {f_bump}"
        );
    }

    #[test]
    fn damper_asymmetry_ratio_greater_than_one() {
        let d = SuspensionDamper::default_sports();
        assert!(d.asymmetry_ratio() > 1.0);
    }

    #[test]
    fn damper_high_speed_force_increases() {
        let d = SuspensionDamper::default_sports();
        let f_low = d.force(d.speed_threshold * 0.5).abs();
        let f_high = d.force(d.speed_threshold * 5.0).abs();
        assert!(f_high > f_low);
    }

    // --- KinematicSolver ---

    #[test]
    fn kinematic_solver_returns_result() {
        let k = SuspensionKinematic::default_front();
        let solver = KinematicSolver::new(k);
        let res = solver.iterate_geometry(0.02);
        assert!(res.scrub >= 0.0);
    }

    #[test]
    fn kinematic_solver_with_double_wishbone() {
        let k = SuspensionKinematic::default_front();
        let dw = DoubleWishbone::default_sports();
        let solver = KinematicSolver::new(k).with_double_wishbone(dw.clone());
        let res = solver.iterate_geometry(0.03);
        let expected_camber = dw.camber_at_bump(0.03);
        assert!((res.camber - expected_camber).abs() < 1e-12);
    }

    #[test]
    fn kinematic_solver_sweep_travel_length() {
        let k = SuspensionKinematic::default_front();
        let solver = KinematicSolver::new(k);
        let results = solver.sweep_travel(-0.05, 0.05, 11);
        assert_eq!(results.len(), 11);
    }

    #[test]
    fn kinematic_solver_sweep_first_last_travel() {
        let k = SuspensionKinematic::default_front();
        let solver = KinematicSolver::new(k);
        let results = solver.sweep_travel(-0.05, 0.05, 11);
        assert!((results[0].0 - (-0.05)).abs() < 1e-10);
        assert!((results[10].0 - 0.05).abs() < 1e-10);
    }

    // --- ComfortAnalysis ---

    #[test]
    fn comfort_body_natural_freq_range() {
        let c = ComfortAnalysis::default_saloon();
        let f = c.body_natural_freq();
        // Typical road car: 1–2 Hz
        assert!(f > 0.5 && f < 3.0, "body natural freq out of range: {f}");
    }

    #[test]
    fn comfort_wheel_freq_higher_than_body() {
        let c = ComfortAnalysis::default_saloon();
        assert!(c.wheel_natural_freq() > c.body_natural_freq());
    }

    #[test]
    fn comfort_damping_ratio_reasonable() {
        let c = ComfortAnalysis::default_saloon();
        let zeta = c.damping_ratio();
        assert!(
            zeta > 0.1 && zeta < 1.0,
            "damping ratio out of range: {zeta}"
        );
    }

    #[test]
    fn comfort_rms_acceleration_positive() {
        let c = ComfortAnalysis::default_saloon();
        let rms = c.rms_body_acceleration(30.0);
        assert!(rms > 0.0, "RMS acceleration should be positive: {rms}");
    }

    #[test]
    fn comfort_discomfort_index_increases_with_speed() {
        let c = ComfortAnalysis::default_saloon();
        let d1 = c.discomfort_index(20.0);
        let d2 = c.discomfort_index(40.0);
        assert!(
            d2 > d1,
            "discomfort should increase with speed: {d1} vs {d2}"
        );
    }

    #[test]
    fn road_class_roughness_increases() {
        assert!(RoadClass::B.roughness_coeff() > RoadClass::A.roughness_coeff());
        assert!(RoadClass::C.roughness_coeff() > RoadClass::B.roughness_coeff());
        assert!(RoadClass::E.roughness_coeff() > RoadClass::D.roughness_coeff());
    }

    // --- NVH ---

    #[test]
    fn nvh_bounce_frequency_range() {
        let nvh = Nvh::default_saloon();
        let f = nvh.bounce_frequency();
        assert!(f > 0.5 && f < 3.0, "bounce frequency out of range: {f}");
    }

    #[test]
    fn nvh_pitch_frequency_reasonable() {
        let nvh = Nvh::default_saloon();
        let f = nvh.pitch_frequency();
        assert!(f > 0.5 && f < 5.0, "pitch frequency out of range: {f}");
    }

    #[test]
    fn nvh_roll_frequency_reasonable() {
        let nvh = Nvh::default_saloon();
        let f = nvh.roll_frequency();
        assert!(f > 0.0 && f < 10.0, "roll frequency out of range: {f}");
    }

    #[test]
    fn nvh_bounce_damping_ratio_sub_critical() {
        let nvh = Nvh::default_saloon();
        let zeta = nvh.bounce_damping_ratio();
        assert!(zeta < 1.0, "bounce should be under-damped: {zeta}");
    }

    #[test]
    fn nvh_transmissibility_at_resonance_above_one() {
        let nvh = Nvh::default_saloon();
        let fn_ = nvh.bounce_frequency();
        let t = nvh.transmissibility(fn_);
        assert!(
            t > 1.0,
            "transmissibility at resonance should exceed 1: {t}"
        );
    }

    #[test]
    fn nvh_transmissibility_drops_at_high_freq() {
        let nvh = Nvh::default_saloon();
        let t_low = nvh.transmissibility(1.0);
        let t_high = nvh.transmissibility(20.0);
        assert!(
            t_high < t_low,
            "transmissibility should drop at high frequency: {t_low} vs {t_high}"
        );
    }

    #[test]
    fn nvh_peak_transmissibility_positive() {
        let nvh = Nvh::default_saloon();
        let pt = nvh.peak_transmissibility();
        assert!(pt > 1.0, "peak transmissibility should exceed 1: {pt}");
    }

    // --- vec3 helpers ---

    #[test]
    fn vec3_helpers_basic() {
        let a = [1.0, 2.0, 3.0];
        let b = [4.0, 5.0, 6.0];
        let sum = vec3_add(a, b);
        assert_eq!(sum, [5.0, 7.0, 9.0]);
        let diff = vec3_sub(b, a);
        assert_eq!(diff, [3.0, 3.0, 3.0]);
        assert!((vec3_dot(a, b) - 32.0).abs() < 1e-10);
        let cross = vec3_cross(a, b);
        // a × b = [-3, 6, -3]
        assert!((cross[0] - (-3.0)).abs() < 1e-10);
        assert!((cross[1] - 6.0).abs() < 1e-10);
        assert!((cross[2] - (-3.0)).abs() < 1e-10);
    }

    #[test]
    fn vec3_norm_basic() {
        let a = [3.0, 4.0, 0.0];
        assert!((vec3_norm(a) - 5.0).abs() < 1e-10);
    }

    #[test]
    fn vec3_scale_basic() {
        let a = [1.0, 2.0, 3.0];
        let scaled = vec3_scale(a, 2.0);
        assert_eq!(scaled, [2.0, 4.0, 6.0]);
    }
}
