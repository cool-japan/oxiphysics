//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::*;
use oxiphysics_core::math::Real;

/// Anti-roll bar stiffness model — computes the roll-restoring torque for a
/// given roll angle.
///
/// This extends the basic `AntiRollBar` with a nonlinear stiffness curve.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct AntiRollBarNonlinear {
    /// Linear stiffness (N·m/rad).
    pub linear_rate: f64,
    /// Cubic stiffness hardening coefficient (N·m/rad³).
    pub cubic_rate: f64,
    /// Maximum effective torque (N·m) — clipped at this value.
    pub max_torque: f64,
}
#[allow(dead_code)]
impl AntiRollBarNonlinear {
    /// Construct from rates and torque limit.
    pub fn new(linear_rate: f64, cubic_rate: f64, max_torque: f64) -> Self {
        Self {
            linear_rate,
            cubic_rate,
            max_torque,
        }
    }
    /// Compute the anti-roll torque (N·m) for a given roll angle (rad).
    pub fn torque(&self, roll_rad: f64) -> f64 {
        let t = self.linear_rate * roll_rad + self.cubic_rate * roll_rad.powi(3);
        t.clamp(-self.max_torque, self.max_torque)
    }
    /// Equivalent linear stiffness at a given roll angle.
    pub fn effective_stiffness(&self, roll_rad: f64) -> f64 {
        self.linear_rate + 3.0 * self.cubic_rate * roll_rad * roll_rad
    }
    /// Maximum available torque gradient (tangent stiffness at max torque point).
    pub fn peak_stiffness(&self) -> f64 {
        self.linear_rate
    }
}
/// Asymmetric damper with different compression and rebound coefficients.
///
/// Optionally supports a piecewise-linear high-speed knee point.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct AsymmetricDamper {
    /// Low-speed compression damping coefficient (N*s/m).
    pub compression_coeff: f64,
    /// Low-speed rebound damping coefficient (N*s/m).
    pub rebound_coeff: f64,
    /// Velocity knee point for high-speed compression (m/s).
    /// Above this speed, a different coefficient is used.
    pub compression_knee: f64,
    /// High-speed compression coefficient (N*s/m).
    pub compression_high: f64,
    /// Velocity knee point for high-speed rebound (m/s).
    pub rebound_knee: f64,
    /// High-speed rebound coefficient (N*s/m).
    pub rebound_high: f64,
}
#[allow(dead_code)]
impl AsymmetricDamper {
    /// Create a simple asymmetric damper with no high-speed transition.
    pub fn new(compression: f64, rebound: f64) -> Self {
        Self {
            compression_coeff: compression,
            rebound_coeff: rebound,
            compression_knee: f64::MAX,
            compression_high: compression,
            rebound_knee: f64::MAX,
            rebound_high: rebound,
        }
    }
    /// Create a damper with high-speed knee points.
    #[allow(clippy::too_many_arguments)]
    pub fn with_knee(
        compression: f64,
        rebound: f64,
        comp_knee: f64,
        comp_high: f64,
        reb_knee: f64,
        reb_high: f64,
    ) -> Self {
        Self {
            compression_coeff: compression,
            rebound_coeff: rebound,
            compression_knee: comp_knee,
            compression_high: comp_high,
            rebound_knee: reb_knee,
            rebound_high: reb_high,
        }
    }
    /// Compute the damping force.
    ///
    /// Positive velocity = extension (rebound); negative = compression.
    /// Force opposes the velocity direction.
    pub fn force(&self, velocity: f64) -> f64 {
        if velocity >= 0.0 {
            if velocity <= self.rebound_knee {
                -self.rebound_coeff * velocity
            } else {
                let base = -self.rebound_coeff * self.rebound_knee;
                base - self.rebound_high * (velocity - self.rebound_knee)
            }
        } else {
            let speed = velocity.abs();
            if speed <= self.compression_knee {
                -self.compression_coeff * velocity
            } else {
                let base = self.compression_coeff * self.compression_knee;
                base + self.compression_high * (speed - self.compression_knee)
            }
        }
    }
    /// Ratio of rebound to compression damping (low-speed).
    pub fn rebound_ratio(&self) -> f64 {
        if self.compression_coeff > 1e-10 {
            self.rebound_coeff / self.compression_coeff
        } else {
            1.0
        }
    }
}
/// Parameters for the classic quarter-car 2-DOF suspension model.
///
/// Two masses:
/// - Sprung mass `ms` (body, kg)
/// - Unsprung mass `mu` (wheel+hub, kg)
///
/// Connected by:
/// - Suspension spring `ks` (N/m) + damper `cs` (N·s/m)
/// - Tyre stiffness `kt` (N/m, no tyre damping)
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct QuarterCarModel {
    /// Sprung mass (kg).
    pub ms: f64,
    /// Unsprung mass (kg).
    pub mu: f64,
    /// Suspension spring stiffness (N/m).
    pub ks: f64,
    /// Tyre stiffness (N/m).
    pub kt: f64,
    /// Suspension damping coefficient (N·s/m).
    pub cs: f64,
}
impl QuarterCarModel {
    /// Create a new quarter-car model.
    #[allow(dead_code)]
    pub fn new(ms: f64, mu: f64, ks: f64, kt: f64, cs: f64) -> Self {
        Self { ms, mu, ks, kt, cs }
    }
    /// Static deflections due to gravity `g`.
    ///
    /// Returns `(xs_static, xu_static)` — the static compressions of the
    /// suspension spring and tyre spring respectively.
    #[allow(dead_code)]
    pub fn static_equilibrium(&self, g: f64) -> (f64, f64) {
        let xs = (self.ms + self.mu) * g / self.ks;
        let xu = (self.ms + self.mu) * g / self.kt;
        (xs, xu)
    }
    /// Natural frequency of the ride mode (sprung mass on suspension), Hz.
    ///
    /// `f_ride = (1/2π) * √(ks / ms)`
    #[allow(dead_code)]
    pub fn ride_frequency(&self) -> f64 {
        if self.ms <= 0.0 || self.ks <= 0.0 {
            return 0.0;
        }
        (self.ks / self.ms).sqrt() / (2.0 * std::f64::consts::PI)
    }
    /// Natural frequency of the wheel-hop mode (unsprung mass on tyre), Hz.
    ///
    /// `f_hop = (1/2π) * √((ks + kt) / mu)`
    #[allow(dead_code)]
    pub fn wheel_hop_frequency(&self) -> f64 {
        if self.mu <= 0.0 {
            return 0.0;
        }
        let k_eff = self.ks + self.kt;
        (k_eff / self.mu).sqrt() / (2.0 * std::f64::consts::PI)
    }
    /// Damping ratio of the sprung-mass mode.
    ///
    /// `ζ = cs / (2 * √(ks * ms))`
    #[allow(dead_code)]
    pub fn damping_ratio_sprung(&self) -> f64 {
        let c_crit = 2.0 * (self.ks * self.ms).sqrt();
        if c_crit < 1e-15 {
            return 0.0;
        }
        self.cs / c_crit
    }
    /// Sprung-mass transmissibility |H(ω)| at excitation frequency ω (rad/s).
    ///
    /// Uses the simplified 2-DOF transfer function magnitude (road displacement
    /// to sprung mass displacement), ignoring tyre damping.
    ///
    /// This is computed numerically via the equations of motion in the
    /// frequency domain.  The result is the ratio |Xs / Xr|.
    #[allow(dead_code)]
    pub fn sprung_mass_transmissibility(&self, omega: f64) -> f64 {
        let ms = self.ms;
        let mu = self.mu;
        let ks = self.ks;
        let kt = self.kt;
        let cs = self.cs;
        let w2 = omega * omega;
        let a11_r = -ms * w2 + ks;
        let a11_i = cs * omega;
        let a12_r = -ks;
        let a12_i = -cs * omega;
        let a21_r = -ks;
        let a21_i = -cs * omega;
        let a22_r = -mu * w2 + ks + kt;
        let a22_i = cs * omega;
        let det_r = a11_r * a22_r - a11_i * a22_i - (a12_r * a21_r - a12_i * a21_i);
        let det_i = a11_r * a22_i + a11_i * a22_r - (a12_r * a21_i + a12_i * a21_r);
        let det_mag2 = det_r * det_r + det_i * det_i;
        if det_mag2 < 1e-30 {
            return 0.0;
        }
        let num_r = -(kt * a12_r);
        let num_i = -(kt * a12_i);
        let xs_r = (num_r * det_r + num_i * det_i) / det_mag2;
        let xs_i = (num_i * det_r - num_r * det_i) / det_mag2;
        (xs_r * xs_r + xs_i * xs_i).sqrt()
    }
    /// Perform one explicit Euler integration step of the 2-DOF quarter-car.
    ///
    /// # Arguments
    /// * `state`        – mutable state (velocities and positions)
    /// * `dt`           – time step (s)
    /// * `road_velocity` – velocity of the road input at the tyre contact (m/s)
    #[allow(dead_code)]
    pub fn step(&self, state: &mut QuarterCarState, dt: f64, road_velocity: f64) {
        let f_susp = self.ks * (state.xu - state.xs) + self.cs * (state.vu - state.vs);
        let f_tyre = self.kt * state.xu;
        let as_ = f_susp / self.ms;
        let au = (-f_tyre + f_susp + self.kt * road_velocity * dt / dt) / self.mu;
        state.vs += as_ * dt;
        state.vu += au * dt;
        state.xs += state.vs * dt;
        state.xu += state.vu * dt;
    }
}
/// Progressive (nonlinear) suspension model.
///
/// The stiffness increases as the spring compresses further.
///
/// `F = k * d * (1 + progressive_factor * d / rest_length) - c * v`
#[derive(Debug, Clone)]
pub struct ProgressiveSuspension {
    /// Factor controlling how much stiffer the spring becomes at full compression.
    pub progressive_factor: Real,
}
impl ProgressiveSuspension {
    /// Create a new progressive suspension with the given factor.
    pub fn new(progressive_factor: Real) -> Self {
        Self { progressive_factor }
    }
}
/// Skyhook active suspension controller.
///
/// Applies a damping force proportional to the absolute velocity of the sprung
/// mass (as if damping against the sky, not the road).  When the absolute and
/// relative velocities have the same sign the sky-hook force is used; otherwise
/// the passive damper is used to avoid energy injection.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct SkyhookDamper {
    /// Sky-hook damping coefficient (N·s/m).
    pub c_skyhook: f64,
    /// Passive (fallback) damping coefficient (N·s/m).
    pub c_passive: f64,
}
impl SkyhookDamper {
    /// Create a skyhook damper.
    #[allow(dead_code)]
    pub fn new(c_skyhook: f64, c_passive: f64) -> Self {
        Self {
            c_skyhook,
            c_passive,
        }
    }
    /// Compute the actuator force.
    ///
    /// # Arguments
    /// * `abs_velocity` – absolute velocity of the sprung mass (m/s)
    /// * `rel_velocity` – relative velocity between sprung and unsprung (m/s)
    ///
    /// # Returns
    /// The controlled damping force (N).  Positive = pushing sprung mass up.
    #[allow(dead_code)]
    pub fn force(&self, abs_velocity: f64, rel_velocity: f64) -> f64 {
        if abs_velocity.abs() < 1e-12 {
            return 0.0;
        }
        if abs_velocity * rel_velocity > 0.0 {
            -self.c_skyhook * abs_velocity
        } else {
            -self.c_passive * rel_velocity
        }
    }
}
/// Progressive bump stop that engages softly and hardens with further
/// compression, with an optional velocity-softening factor.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct ProgressiveBumpStop {
    /// Travel (m) at which the bump stop engages (measured from design position).
    pub engage_travel: f64,
    /// Initial (soft) rate (N/m).
    pub soft_rate: f64,
    /// Hard rate (N/m) at full engagement.
    pub hard_rate: f64,
    /// Travel over which the rate transitions from soft to hard (m).
    pub transition_length: f64,
}
#[allow(dead_code)]
impl ProgressiveBumpStop {
    /// Construct from engagement parameters.
    pub fn new(engage_travel: f64, soft_rate: f64, hard_rate: f64, transition_length: f64) -> Self {
        Self {
            engage_travel,
            soft_rate,
            hard_rate,
            transition_length,
        }
    }
    /// Compute the bump-stop force (N) for a given suspension travel (m).
    ///
    /// Positive force resists further compression.
    pub fn force(&self, travel: f64) -> f64 {
        if travel <= self.engage_travel {
            return 0.0;
        }
        let delta = travel - self.engage_travel;
        let t = (delta / self.transition_length.max(1e-6)).min(1.0);
        let effective_rate = self.soft_rate + t * (self.hard_rate - self.soft_rate);
        effective_rate * delta
    }
    /// Energy stored in the bump stop (J).
    pub fn stored_energy(&self, travel: f64) -> f64 {
        if travel <= self.engage_travel {
            return 0.0;
        }
        let f = self.force(travel);
        let delta = travel - self.engage_travel;
        0.5 * f * delta
    }
}
/// Bump stop model with progressive stiffness.
///
/// The bump stop engages when jounce exceeds a gap distance and provides
/// a sharply rising force.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct BumpStop {
    /// Distance from rest before the bump stop engages (m).
    pub gap: f64,
    /// Initial stiffness of the bump stop (N/m).
    pub stiffness: f64,
    /// Progressive exponent: F = k * (penetration)^exponent.
    pub exponent: f64,
}
#[allow(dead_code)]
impl BumpStop {
    /// Create a bump stop with the given gap, stiffness, and exponent.
    pub fn new(gap: f64, stiffness: f64, exponent: f64) -> Self {
        Self {
            gap,
            stiffness,
            exponent,
        }
    }
    /// Default bump stop: 0.10 m gap, 100 kN/m, exponent 2.0.
    pub fn default_bump() -> Self {
        Self {
            gap: 0.10,
            stiffness: 100_000.0,
            exponent: 2.0,
        }
    }
    /// Compute the bump stop force for a given jounce (compression).
    ///
    /// Returns zero if jounce < gap; otherwise returns k * (jounce - gap)^exp.
    pub fn force(&self, jounce: f64) -> f64 {
        let penetration = jounce - self.gap;
        if penetration <= 0.0 {
            0.0
        } else {
            self.stiffness * penetration.powf(self.exponent)
        }
    }
    /// Combined spring + bump stop force.
    pub fn combined_force(&self, jounce: f64, spring_force_val: f64) -> f64 {
        spring_force_val + self.force(jounce)
    }
}
/// Double-wishbone suspension parameterised by scalar arm lengths.
///
/// A simpler, scalar representation suitable for calculating camber gain
/// without full 3-D geometry.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct DoubleWishboneSimple {
    /// Length of the upper A-arm (m).
    pub upper_arm_length: f64,
    /// Length of the lower A-arm (m).
    pub lower_arm_length: f64,
    /// Spring rate (N/m).
    pub spring_rate: f64,
    /// Damping coefficient (N·s/m).
    pub damping: f64,
    /// Static camber angle (radians; negative = top tilted inward).
    pub camber_angle: f64,
}
#[allow(dead_code)]
impl DoubleWishboneSimple {
    /// Create a new scalar double-wishbone model.
    pub fn new(
        upper_arm_length: f64,
        lower_arm_length: f64,
        spring_rate: f64,
        damping: f64,
        camber_angle: f64,
    ) -> Self {
        Self {
            upper_arm_length,
            lower_arm_length,
            spring_rate,
            damping,
            camber_angle,
        }
    }
    /// Camber gain (radians) as a function of wheel travel (m).
    ///
    /// With unequal arm lengths the roll-centre rises on bump, producing
    /// camber gain.  The simplified formula is:
    ///
    /// Δcamber ≈ static_camber − (1/lower_arm − 1/upper_arm) * travel
    pub fn wheel_travel_to_camber(&self, travel: f64) -> f64 {
        let lower = self.lower_arm_length.max(1e-9);
        let upper = self.upper_arm_length.max(1e-9);
        let gain = 1.0 / upper - 1.0 / lower;
        self.camber_angle - gain * travel
    }
}
/// Double-wishbone suspension geometry.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct DoubleWishbone {
    /// Pivot point of the upper A-arm (body frame, metres).
    pub upper_a_arm: [f64; 3],
    /// Pivot point of the lower A-arm (body frame, metres).
    pub lower_a_arm: [f64; 3],
    /// Wheel centre at design (zero-jounce) position (body frame, metres).
    pub wheel_center: [f64; 3],
    /// Static camber angle (radians; negative = top tilted inward).
    pub camber_angle: f64,
}
impl DoubleWishbone {
    /// Create a double-wishbone with zero static camber.
    #[allow(dead_code)]
    pub fn new(upper: [f64; 3], lower: [f64; 3], wheel: [f64; 3]) -> Self {
        Self {
            upper_a_arm: upper,
            lower_a_arm: lower,
            wheel_center: wheel,
            camber_angle: 0.0,
        }
    }
    /// Simplified camber gain: camber change (radians) per metre of jounce.
    #[allow(dead_code)]
    pub fn compute_camber_gain(&self, jounce: f64) -> f64 {
        let arm_spread = (self.upper_a_arm[1] - self.lower_a_arm[1]).abs().max(1e-6);
        let gain_per_metre = -1.0 / arm_spread;
        self.camber_angle + gain_per_metre * jounce
    }
    /// Return the wheel-centre position in body frame at a given jounce.
    #[allow(dead_code)]
    pub fn wheel_center_at_jounce(&self, jounce: f64) -> [f64; 3] {
        let camber = self.compute_camber_gain(jounce);
        [
            self.wheel_center[0],
            self.wheel_center[1] + camber * jounce * 0.1,
            self.wheel_center[2] + jounce,
        ]
    }
    /// Instantaneous roll centre height (simplified).
    ///
    /// For a symmetric double-wishbone, the roll centre is at the intersection
    /// of the lines from the tyre contact patch through the instant centre.
    /// This simplified model returns the average Z of the two arm pivots.
    #[allow(dead_code)]
    pub fn roll_center_height(&self) -> f64 {
        (self.upper_a_arm[2] + self.lower_a_arm[2]) * 0.5
    }
    /// Effective arm length ratio (upper / lower).
    #[allow(dead_code)]
    pub fn arm_length_ratio(&self) -> f64 {
        let upper_len = ((self.upper_a_arm[0] - self.wheel_center[0]).powi(2)
            + (self.upper_a_arm[1] - self.wheel_center[1]).powi(2)
            + (self.upper_a_arm[2] - self.wheel_center[2]).powi(2))
        .sqrt();
        let lower_len = ((self.lower_a_arm[0] - self.wheel_center[0]).powi(2)
            + (self.lower_a_arm[1] - self.wheel_center[1]).powi(2)
            + (self.lower_a_arm[2] - self.wheel_center[2]).powi(2))
        .sqrt();
        if lower_len > 1e-10 {
            upper_len / lower_len
        } else {
            1.0
        }
    }
}
/// Simplified vehicle roll model combining front and rear suspension anti-roll
/// bars and lateral load transfer.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct VehicleRollModel {
    /// Front anti-roll bar stiffness (N·m/rad).
    pub front_arb_rate: f64,
    /// Rear anti-roll bar stiffness (N·m/rad).
    pub rear_arb_rate: f64,
    /// Front suspension roll stiffness (N·m/rad, from springs + geometry).
    pub front_roll_stiffness: f64,
    /// Rear suspension roll stiffness (N·m/rad).
    pub rear_roll_stiffness: f64,
    /// Vehicle total mass (kg).
    pub mass: f64,
    /// CG height (m).
    pub cg_height: f64,
    /// Front track width (m).
    pub track: f64,
}
#[allow(dead_code)]
impl VehicleRollModel {
    /// Construct from vehicle and suspension parameters.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        front_arb_rate: f64,
        rear_arb_rate: f64,
        front_roll_stiffness: f64,
        rear_roll_stiffness: f64,
        mass: f64,
        cg_height: f64,
        track: f64,
    ) -> Self {
        Self {
            front_arb_rate,
            rear_arb_rate,
            front_roll_stiffness,
            rear_roll_stiffness,
            mass,
            cg_height,
            track,
        }
    }
    /// Total roll stiffness (N·m/rad) — sum of all sources.
    pub fn total_roll_stiffness(&self) -> f64 {
        self.front_arb_rate
            + self.rear_arb_rate
            + self.front_roll_stiffness
            + self.rear_roll_stiffness
    }
    /// Steady-state roll angle (rad) for a given lateral acceleration (m/s²).
    pub fn roll_angle(&self, lat_accel: f64) -> f64 {
        let k_total = self.total_roll_stiffness();
        if k_total < 1e-6 {
            return 0.0;
        }
        let rolling_moment = self.mass * lat_accel * self.cg_height;
        rolling_moment / k_total
    }
    /// Lateral load transfer at the front axle (N) due to roll.
    pub fn front_load_transfer(&self, lat_accel: f64) -> f64 {
        if self.track < 1e-6 {
            return 0.0;
        }
        let phi = self.roll_angle(lat_accel);
        let front_moment = (self.front_arb_rate + self.front_roll_stiffness) * phi;
        front_moment / self.track
    }
    /// Lateral load transfer at the rear axle (N).
    pub fn rear_load_transfer(&self, lat_accel: f64) -> f64 {
        if self.track < 1e-6 {
            return 0.0;
        }
        let phi = self.roll_angle(lat_accel);
        let rear_moment = (self.rear_arb_rate + self.rear_roll_stiffness) * phi;
        rear_moment / self.track
    }
    /// Front/rear load transfer split ratio.
    ///
    /// Returns `(front_fraction, rear_fraction)` summing to 1.
    pub fn load_transfer_split(&self) -> (f64, f64) {
        let k_f = self.front_arb_rate + self.front_roll_stiffness;
        let k_r = self.rear_arb_rate + self.rear_roll_stiffness;
        let k_tot = k_f + k_r;
        if k_tot < 1e-9 {
            return (0.5, 0.5);
        }
        (k_f / k_tot, k_r / k_tot)
    }
}
/// Parameters for a spring-damper suspension unit.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct SuspensionParams {
    /// Spring stiffness (N/m).
    pub spring_stiffness: f64,
    /// Damper coefficient (N*s/m).
    pub damper_coeff: f64,
    /// Natural (unloaded) length of the spring (m).
    pub rest_length: f64,
    /// Maximum jounce travel from rest (m).
    pub max_travel: f64,
    /// Maximum rebound travel from rest (negative, m).
    pub min_travel: f64,
    /// Anti-roll bar stiffness contribution (N/m).
    pub anti_roll_stiffness: f64,
}
impl SuspensionParams {
    /// Create suspension parameters with sensible defaults.
    #[allow(dead_code)]
    pub fn new(stiffness: f64, damping: f64, rest_length: f64) -> Self {
        Self {
            spring_stiffness: stiffness,
            damper_coeff: damping,
            rest_length,
            max_travel: 0.15,
            min_travel: -0.15,
            anti_roll_stiffness: 0.0,
        }
    }
}
/// Describes the state of a wheel's contact with the ground.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct WheelContact {
    /// Contact point in world space (metres).
    pub position: [f64; 3],
    /// Surface normal at the contact point (unit vector).
    pub normal: [f64; 3],
    /// Longitudinal slip ratio (dimensionless).
    pub slip_ratio: f64,
    /// Lateral slip angle (radians).
    pub slip_angle: f64,
    /// Normal (vertical) load on the tyre (N).
    pub vertical_load: f64,
    /// Whether the wheel is in contact with the ground.
    pub in_contact: bool,
}
/// McPherson strut suspension geometry.
///
/// Models the strut as a single pivot at the lower control arm and a top mount.
/// The wheel centre moves along the strut axis; camber change comes from the
/// angle between the strut axis and the vertical.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct McPhersonStrut {
    /// Lower ball-joint position (body frame, metres).
    pub lower_pivot: [f64; 3],
    /// Top mount (strut tower) position (body frame, metres).
    pub top_mount: [f64; 3],
    /// Wheel centre at design ride height (body frame, metres).
    pub wheel_center: [f64; 3],
    /// Static camber angle (radians; negative = top tilted inward).
    pub static_camber: f64,
}
#[allow(dead_code)]
impl McPhersonStrut {
    /// Create a McPherson strut with zero static camber.
    pub fn new(lower: [f64; 3], top: [f64; 3], wheel: [f64; 3]) -> Self {
        Self {
            lower_pivot: lower,
            top_mount: top,
            wheel_center: wheel,
            static_camber: 0.0,
        }
    }
    /// Strut axis unit direction vector (from lower to top).
    pub fn strut_axis(&self) -> [f64; 3] {
        let dx = self.top_mount[0] - self.lower_pivot[0];
        let dy = self.top_mount[1] - self.lower_pivot[1];
        let dz = self.top_mount[2] - self.lower_pivot[2];
        let len = (dx * dx + dy * dy + dz * dz).sqrt().max(1e-10);
        [dx / len, dy / len, dz / len]
    }
    /// Strut length (design position).
    pub fn strut_length(&self) -> f64 {
        let dx = self.top_mount[0] - self.lower_pivot[0];
        let dy = self.top_mount[1] - self.lower_pivot[1];
        let dz = self.top_mount[2] - self.lower_pivot[2];
        (dx * dx + dy * dy + dz * dz).sqrt()
    }
    /// Camber gain per metre of jounce.
    ///
    /// For a McPherson strut the camber gain is approximately the inverse of the
    /// strut length times the sine of the strut inclination angle.
    pub fn camber_gain(&self, jounce: f64) -> f64 {
        let axis = self.strut_axis();
        let cos_incl = axis[2].abs().min(1.0);
        let sin_incl = (1.0 - cos_incl * cos_incl).sqrt();
        let strut_len = self.strut_length().max(1e-6);
        self.static_camber - sin_incl * jounce / strut_len
    }
    /// Wheel centre position at a given jounce (positive = compression).
    pub fn wheel_center_at_jounce(&self, jounce: f64) -> [f64; 3] {
        let _axis = self.strut_axis();
        let camber = self.camber_gain(jounce);
        [
            self.wheel_center[0] + camber * jounce * 0.05,
            self.wheel_center[1] + camber * jounce * 0.1,
            self.wheel_center[2] + jounce,
        ]
    }
    /// Kingpin inclination angle (radians).
    ///
    /// Approximated from the strut axis projected onto the YZ plane.
    pub fn kingpin_inclination(&self) -> f64 {
        let axis = self.strut_axis();
        axis[1].atan2(axis[2]).abs()
    }
}
/// Represents an anti-roll (stabiliser) bar connecting left and right wheels.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct AntiRollBar {
    /// Torsional stiffness (N/m).
    pub stiffness: f64,
}
impl AntiRollBar {
    /// Compute the anti-roll forces acting on each side.
    ///
    /// Returns `(left_force, right_force)` which are equal and opposite.
    #[allow(dead_code)]
    pub fn compute_force(&self, left_travel: f64, right_travel: f64) -> (f64, f64) {
        let f = self.stiffness * (left_travel - right_travel);
        (f, -f)
    }
    /// Compute the anti-roll moment (N*m) about the roll axis.
    ///
    /// `M = k * track_width * (left_travel - right_travel)`
    #[allow(dead_code)]
    pub fn roll_moment(&self, left_travel: f64, right_travel: f64, track_width: f64) -> f64 {
        self.stiffness * track_width * (left_travel - right_travel)
    }
}
impl AntiRollBar {
    /// Anti-roll moment (N·m) for given left and right suspension travels.
    ///
    /// `M = k * (left_travel − right_travel)`
    #[allow(dead_code)]
    pub fn moment(&self, left_travel: f64, right_travel: f64) -> f64 {
        self.stiffness * (left_travel - right_travel)
    }
}
/// Dynamic state for the quarter-car 2-DOF model.
///
/// Positions are measured from the static equilibrium position.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct QuarterCarState {
    /// Sprung mass displacement from equilibrium (m, positive = up).
    pub xs: f64,
    /// Unsprung mass displacement from equilibrium (m).
    pub xu: f64,
    /// Sprung mass velocity (m/s).
    pub vs: f64,
    /// Unsprung mass velocity (m/s).
    pub vu: f64,
}
impl QuarterCarState {
    /// Create state at rest at the static equilibrium position.
    #[allow(dead_code)]
    pub fn at_equilibrium(model: &QuarterCarModel, _g: f64) -> Self {
        let _ = model;
        Self {
            xs: 0.0,
            xu: 0.0,
            vs: 0.0,
            vu: 0.0,
        }
    }
}
/// McPherson strut spring-damper unit (spring/damper focus).
///
/// Captures the spring rate, damping, rest length, caster angle, and
/// kingpin inclination as scalar parameters.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct McPhersonSuspension {
    /// Spring stiffness (N/m).
    pub spring_rate: f64,
    /// Damping coefficient (N·s/m).
    pub damping: f64,
    /// Unloaded spring length (m).
    pub rest_length: f64,
    /// Caster angle (radians).
    pub caster_angle: f64,
    /// Kingpin inclination angle (radians).
    pub kingpin_inclination: f64,
}
#[allow(dead_code)]
impl McPhersonSuspension {
    /// Create a new McPherson suspension unit.
    pub fn new(
        spring_rate: f64,
        damping: f64,
        rest_length: f64,
        caster_angle: f64,
        kingpin_inclination: f64,
    ) -> Self {
        Self {
            spring_rate,
            damping,
            rest_length,
            caster_angle,
            kingpin_inclination,
        }
    }
    /// Linear spring force for a given compression (m).
    ///
    /// F_spring = k * compression
    pub fn spring_force(&self, compression: f64) -> f64 {
        self.spring_rate * compression
    }
    /// Linear damper force for a given velocity (m/s).
    ///
    /// F_damper = -c * velocity  (opposes motion)
    pub fn damper_force(&self, velocity: f64) -> f64 {
        -self.damping * velocity
    }
}
/// Full McPherson strut kinematic model.
///
/// Tracks the strut mount point, lower ball joint, wheel centre, and the
/// spring/damper aligned along the strut axis.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct McPhersonKinematic {
    /// Strut mount point on the body (unloaded position), m.
    pub mount_body: [f64; 3],
    /// Lower ball-joint position on the lower arm (unloaded), m.
    pub lower_ball_joint: [f64; 3],
    /// Wheel centre (unloaded), m.
    pub wheel_centre: [f64; 3],
    /// Spring free length (m).
    pub spring_free_length: f64,
    /// Spring stiffness (N/m).
    pub spring_k: f64,
    /// Damper rate (N·s/m).
    pub damper_c: f64,
}
#[allow(dead_code)]
impl McPhersonKinematic {
    /// Construct from geometric and compliance parameters.
    pub fn new(
        mount_body: [f64; 3],
        lower_ball_joint: [f64; 3],
        wheel_centre: [f64; 3],
        spring_free_length: f64,
        spring_k: f64,
        damper_c: f64,
    ) -> Self {
        Self {
            mount_body,
            lower_ball_joint,
            wheel_centre,
            spring_free_length,
            spring_k,
            damper_c,
        }
    }
    /// Compute current strut length from mount and lower ball-joint positions
    /// given a wheel travel `delta_z` (positive = bump/compression, m).
    pub fn strut_length(&self, delta_z: f64) -> f64 {
        let bj = [
            self.lower_ball_joint[0],
            self.lower_ball_joint[1],
            self.lower_ball_joint[2] + delta_z,
        ];
        let dx = bj[0] - self.mount_body[0];
        let dy = bj[1] - self.mount_body[1];
        let dz = bj[2] - self.mount_body[2];
        (dx * dx + dy * dy + dz * dz).sqrt()
    }
    /// Spring deflection (m) for a given wheel travel delta_z.
    ///
    /// Positive = compression.
    pub fn spring_deflection(&self, delta_z: f64) -> f64 {
        let current = self.strut_length(delta_z);
        self.spring_free_length - current
    }
    /// Spring force (N) for a given wheel travel.
    ///
    /// Positive = pushing apart (rebound force).
    pub fn spring_force(&self, delta_z: f64) -> f64 {
        self.spring_k * self.spring_deflection(delta_z)
    }
    /// Damper force (N) for given strut velocity (m/s, positive = extension).
    pub fn damper_force(&self, strut_velocity: f64) -> f64 {
        -self.damper_c * strut_velocity
    }
    /// Total strut force (spring + damper).
    pub fn total_force(&self, delta_z: f64, strut_velocity: f64) -> f64 {
        self.spring_force(delta_z) + self.damper_force(strut_velocity)
    }
    /// Caster angle (radians) — angle of the strut axis projected into the XZ plane.
    ///
    /// Positive caster = mount leans rearward of lower ball joint.
    pub fn caster_angle(&self) -> f64 {
        let dx = self.lower_ball_joint[0] - self.mount_body[0];
        let dz = self.lower_ball_joint[2] - self.mount_body[2];
        dx.atan2(dz.abs().max(1e-9))
    }
    /// Camber angle (radians) — angle of the strut axis projected into the YZ plane.
    pub fn camber_angle(&self) -> f64 {
        let dy = self.lower_ball_joint[1] - self.mount_body[1];
        let dz = self.lower_ball_joint[2] - self.mount_body[2];
        dy.atan2(dz.abs().max(1e-9))
    }
    /// Wheel-centre height change due to bump travel.
    ///
    /// Uses a linear approximation based on the lower arm geometry.
    pub fn wheel_centre_height(&self, delta_z: f64) -> f64 {
        self.wheel_centre[2] + delta_z
    }
}
/// Linear spring-damper suspension model.
///
/// Computes force as `F = -k * (x - x0) - c * v`.
#[derive(Debug, Clone, Default)]
pub struct LinearSuspension;
/// Geometric kinematics parameters derived from suspension geometry.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct SuspensionKinematics {
    /// Roll-centre height above the ground plane (m).
    pub roll_center_height: f64,
    /// Scrub radius at the contact patch (m).
    pub scrub_radius: f64,
    /// Caster trail (mechanical trail) (m).
    pub caster_trail: f64,
}
/// Runtime state of a single suspension unit.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct SuspensionState {
    /// Current length of the spring (m).
    pub current_length: f64,
    /// Rate of change of spring length (m/s); positive = extending.
    pub velocity: f64,
    /// Computed spring-damper force (N).
    pub force: f64,
}
/// Computes pitch and roll moment distribution for a 4-wheeled vehicle.
///
/// Uses a simplified quasi-static model based on the centre of mass height,
/// track width, and front/rear axle distances.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct PitchRollDistributor {
    /// Total vehicle mass (kg).
    pub mass: f64,
    /// Centre-of-mass height above the ground (m).
    pub h_cg: f64,
    /// Front track width (m).
    pub track_f: f64,
    /// Front axle distance from CoM (m).
    pub l_f: f64,
    /// Rear axle distance from CoM (m).
    pub l_r: f64,
    /// Gravitational acceleration (m/s²).
    pub g: f64,
}
impl PitchRollDistributor {
    /// Create a new distributor.
    #[allow(dead_code)]
    pub fn new(mass: f64, h_cg: f64, track_f: f64, l_f: f64, l_r: f64, g: f64) -> Self {
        Self {
            mass,
            h_cg,
            track_f,
            l_f,
            l_r,
            g,
        }
    }
    /// Pitch moment (N·m) due to longitudinal acceleration.
    ///
    /// `M_pitch = m * a_long * h_cg`
    #[allow(dead_code)]
    pub fn pitch_moment(&self, long_accel: f64) -> f64 {
        self.mass * long_accel * self.h_cg
    }
    /// Roll moment (N·m) due to lateral acceleration.
    ///
    /// `M_roll = m * a_lat * h_cg`
    #[allow(dead_code)]
    pub fn roll_moment(&self, lat_accel: f64) -> f64 {
        self.mass * lat_accel.abs() * self.h_cg
    }
    /// Front/rear load transfer due to pitch (longitudinal acceleration).
    ///
    /// Returns `(delta_front_N, delta_rear_N)`.
    ///
    /// - Front axle gains load during braking (negative `long_accel`).
    /// - Rear axle gains load during acceleration.
    #[allow(dead_code)]
    pub fn front_rear_load_transfer_pitch(&self, long_accel: f64) -> (f64, f64) {
        let wb = self.l_f + self.l_r;
        if wb < 1e-12 {
            return (0.0, 0.0);
        }
        let m_pitch = self.pitch_moment(long_accel);
        let delta_f = m_pitch * self.l_r / (wb * wb);
        let delta_r = -m_pitch * self.l_f / (wb * wb);
        (delta_f, delta_r)
    }
    /// Left/right load transfer due to roll (lateral acceleration).
    ///
    /// Returns `(delta_left_N, delta_right_N)` where positive = additional load.
    #[allow(dead_code)]
    pub fn left_right_load_transfer_roll(&self, lat_accel: f64) -> (f64, f64) {
        if self.track_f <= 0.0 {
            return (0.0, 0.0);
        }
        let m_roll = self.mass * lat_accel * self.h_cg;
        let delta = m_roll / self.track_f;
        (-delta, delta)
    }
}
/// Full double-wishbone suspension kinematic model.
///
/// Tracks upper and lower control arms, their inboard pivot points, and the
/// wheel spindle location.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct DoubleWishboneKinematic {
    /// Upper inboard pivot (body side), m.
    pub upper_inboard: [f64; 3],
    /// Upper outboard point (at wheel hub), m.
    pub upper_outboard: [f64; 3],
    /// Lower inboard pivot, m.
    pub lower_inboard: [f64; 3],
    /// Lower outboard point (lower ball joint), m.
    pub lower_outboard: [f64; 3],
    /// Wheel-spindle reference position (unloaded), m.
    pub spindle: [f64; 3],
    /// Spring free length (m).
    pub spring_free_length: f64,
    /// Spring stiffness (N/m).
    pub spring_k: f64,
    /// Damper rate (N·s/m).
    pub damper_c: f64,
}
#[allow(dead_code)]
impl DoubleWishboneKinematic {
    /// Construct from geometry and compliance data.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        upper_inboard: [f64; 3],
        upper_outboard: [f64; 3],
        lower_inboard: [f64; 3],
        lower_outboard: [f64; 3],
        spindle: [f64; 3],
        spring_free_length: f64,
        spring_k: f64,
        damper_c: f64,
    ) -> Self {
        Self {
            upper_inboard,
            upper_outboard,
            lower_inboard,
            lower_outboard,
            spindle,
            spring_free_length,
            spring_k,
            damper_c,
        }
    }
    /// Upper arm length (m).
    pub fn upper_arm_length(&self) -> f64 {
        let d = sub3(self.upper_outboard, self.upper_inboard);
        norm3(d)
    }
    /// Lower arm length (m).
    pub fn lower_arm_length(&self) -> f64 {
        let d = sub3(self.lower_outboard, self.lower_inboard);
        norm3(d)
    }
    /// Instant centre height (m) from ground using the 4-bar analogy.
    ///
    /// Approximated as the intersection of lines extending through the upper
    /// and lower arm pivot axes in the YZ plane.
    pub fn instant_centre_height(&self) -> f64 {
        let uy = self.upper_outboard[1] - self.upper_inboard[1];
        let uz = self.upper_outboard[2] - self.upper_inboard[2];
        let ly = self.lower_outboard[1] - self.lower_inboard[1];
        let lz = self.lower_outboard[2] - self.lower_inboard[2];
        let det = uy * lz - uz * ly;
        if det.abs() < 1e-9 {
            return 0.0;
        }
        let dy = self.lower_inboard[1] - self.upper_inboard[1];
        let dz = self.lower_inboard[2] - self.upper_inboard[2];
        let t = (dy * lz - dz * ly) / det;
        self.upper_inboard[2] + t * uz
    }
    /// Roll centre height (m) — estimated as the instant centre height.
    pub fn roll_centre_height(&self) -> f64 {
        self.instant_centre_height()
    }
    /// Camber change per unit of wheel travel (deg/mm).
    ///
    /// Uses a finite-difference approximation with Δz = 1 mm.
    pub fn camber_gain(&self) -> f64 {
        let dz = 0.001;
        let camber_base = self.camber_at_travel(0.0);
        let camber_bump = self.camber_at_travel(dz);
        (camber_bump - camber_base) / (dz * 1000.0)
    }
    /// Camber angle (degrees) at a given wheel travel (m, positive = bump).
    pub fn camber_at_travel(&self, delta_z: f64) -> f64 {
        let upper_z = self.upper_outboard[2] + delta_z;
        let lower_z = self.lower_outboard[2] + delta_z;
        let dz = upper_z - lower_z;
        let dy = self.upper_outboard[1] - self.lower_outboard[1];
        dz.atan2(dy.abs().max(1e-9)).to_degrees()
    }
    /// Spring length at a given wheel travel (m).
    pub fn spring_length(&self, delta_z: f64) -> f64 {
        self.spring_free_length - delta_z
    }
    /// Spring force (N) at a given wheel travel.
    pub fn spring_force(&self, delta_z: f64) -> f64 {
        let defl = delta_z;
        self.spring_k * defl
    }
    /// Total suspension force (N) including damper.
    pub fn total_force(&self, delta_z: f64, velocity: f64) -> f64 {
        self.spring_force(delta_z) - self.damper_c * velocity
    }
}
/// Bump stop with preload, linear rate, and maximum force clamp.
///
/// The force engages as soon as compression exceeds zero (preload already
/// present), rising linearly at `rate` until it reaches `max_force`.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct BumpStopLinear {
    /// Pre-load force present at zero compression (N).
    pub preload: f64,
    /// Linear stiffness rate (N/m).
    pub rate: f64,
    /// Maximum force clamp (N).
    pub max_force: f64,
}
#[allow(dead_code)]
impl BumpStopLinear {
    /// Create a new linear bump stop.
    pub fn new(preload: f64, rate: f64, max_force: f64) -> Self {
        Self {
            preload,
            rate,
            max_force,
        }
    }
    /// Force at a given compression (m).
    ///
    /// Returns `clamp(preload + rate * compression, 0, max_force)`.
    pub fn force(&self, compression: f64) -> f64 {
        if compression <= 0.0 {
            return 0.0;
        }
        (self.preload + self.rate * compression)
            .min(self.max_force)
            .max(0.0)
    }
}
