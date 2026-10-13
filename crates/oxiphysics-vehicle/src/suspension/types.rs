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
pub struct AntiRollBarNonlinear {
    /// Linear stiffness (N·m/rad).
    pub linear_rate: f64,
    /// Cubic stiffness hardening coefficient (N·m/rad³).
    pub cubic_rate: f64,
    /// Maximum effective torque (N·m) — clipped at this value.
    pub max_torque: f64,
}
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
    pub fn new(ms: f64, mu: f64, ks: f64, kt: f64, cs: f64) -> Self {
        Self { ms, mu, ks, kt, cs }
    }
    /// Static deflections due to gravity `g`.
    ///
    /// Returns `(xs_static, xu_static)` — the static compressions of the
    /// suspension spring and tyre spring respectively.  The suspension spring
    /// carries only the sprung mass (`ms·g / ks`) while the tyre carries the
    /// sprung and the unsprung mass (`(ms + mu)·g / kt`).
    pub fn static_equilibrium(&self, g: f64) -> (f64, f64) {
        let xs = self.ms * g / self.ks;
        let xu = (self.ms + self.mu) * g / self.kt;
        (xs, xu)
    }
    /// Natural frequency of the ride mode (sprung mass on suspension), Hz.
    ///
    /// `f_ride = (1/2π) * √(ks / ms)`
    pub fn ride_frequency(&self) -> f64 {
        if self.ms <= 0.0 || self.ks <= 0.0 {
            return 0.0;
        }
        (self.ks / self.ms).sqrt() / (2.0 * std::f64::consts::PI)
    }
    /// Natural frequency of the wheel-hop mode (unsprung mass on tyre), Hz.
    ///
    /// `f_hop = (1/2π) * √((ks + kt) / mu)`
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
    /// Advance the 2-DOF quarter-car by one semi-implicit (symplectic) Euler step.
    ///
    /// All positions are measured from the static-equilibrium configuration, so
    /// gravity is already balanced and does not appear in the equations of
    /// motion.  With `f_susp = ks·(xu − xs) + cs·(vu − vs)` (force the suspension
    /// exerts on the sprung mass, positive upwards) and
    /// `f_tyre = kt·(xu − zr)` (compression of the tyre spring against the road
    /// height `zr`):
    ///
    /// ```text
    /// ms · a_s = +f_susp
    /// mu · a_u = −f_susp − f_tyre
    /// ```
    ///
    /// The suspension force acts with opposite signs on the two masses, so the
    /// damper dissipates energy; the tyre is the only element tied to the road.
    /// The road height is integrated from the supplied road velocity before the
    /// forces are evaluated (`zr ← zr + road_velocity·dt`).
    ///
    /// The damper force is treated explicitly, so `dt` must stay below
    /// `2·m_red / cs` with `m_red = ms·mu / (ms + mu)` and well below the
    /// wheel-hop period for the integration to remain stable.
    ///
    /// # Arguments
    /// * `state`         – mutable state (positions, velocities and road height)
    /// * `dt`            – time step (s)
    /// * `road_velocity` – vertical velocity of the road surface under the tyre
    ///   (m/s, positive = road rising)
    pub fn step(&self, state: &mut QuarterCarState, dt: f64, road_velocity: f64) {
        state.zr += road_velocity * dt;
        let f_susp = self.ks * (state.xu - state.xs) + self.cs * (state.vu - state.vs);
        let f_tyre = self.kt * (state.xu - state.zr);
        let as_ = if self.ms > 0.0 { f_susp / self.ms } else { 0.0 };
        let au = if self.mu > 0.0 {
            (-f_susp - f_tyre) / self.mu
        } else {
            0.0
        };
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
#[derive(Debug, Clone)]
pub struct SkyhookDamper {
    /// Sky-hook damping coefficient (N·s/m).
    pub c_skyhook: f64,
    /// Passive (fallback) damping coefficient (N·s/m).
    pub c_passive: f64,
}
impl SkyhookDamper {
    /// Create a skyhook damper.
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
#[derive(Debug, Clone)]
pub struct BumpStop {
    /// Distance from rest before the bump stop engages (m).
    pub gap: f64,
    /// Initial stiffness of the bump stop (N/m).
    pub stiffness: f64,
    /// Progressive exponent: F = k * (penetration)^exponent.
    pub exponent: f64,
}
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
    pub fn new(upper: [f64; 3], lower: [f64; 3], wheel: [f64; 3]) -> Self {
        Self {
            upper_a_arm: upper,
            lower_a_arm: lower,
            wheel_center: wheel,
            camber_angle: 0.0,
        }
    }
    /// Simplified camber gain: camber change (radians) per metre of jounce.
    pub fn compute_camber_gain(&self, jounce: f64) -> f64 {
        let arm_spread = (self.upper_a_arm[1] - self.lower_a_arm[1]).abs().max(1e-6);
        let gain_per_metre = -1.0 / arm_spread;
        self.camber_angle + gain_per_metre * jounce
    }
    /// Return the wheel-centre position in body frame at a given jounce.
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
    pub fn roll_center_height(&self) -> f64 {
        (self.upper_a_arm[2] + self.lower_a_arm[2]) * 0.5
    }
    /// Effective arm length ratio (upper / lower).
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
impl VehicleRollModel {
    /// Construct from vehicle and suspension parameters.
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
#[derive(Debug, Clone)]
pub struct AntiRollBar {
    /// Torsional stiffness (N/m).
    pub stiffness: f64,
}
impl AntiRollBar {
    /// Compute the anti-roll forces acting on each side.
    ///
    /// Returns `(left_force, right_force)` which are equal and opposite.
    pub fn compute_force(&self, left_travel: f64, right_travel: f64) -> (f64, f64) {
        let f = self.stiffness * (left_travel - right_travel);
        (f, -f)
    }
    /// Compute the anti-roll moment (N*m) about the roll axis.
    ///
    /// `M = k * track_width * (left_travel - right_travel)`
    pub fn roll_moment(&self, left_travel: f64, right_travel: f64, track_width: f64) -> f64 {
        self.stiffness * track_width * (left_travel - right_travel)
    }
}
impl AntiRollBar {
    /// Anti-roll moment (N·m) for given left and right suspension travels.
    ///
    /// `M = k * (left_travel − right_travel)`
    pub fn moment(&self, left_travel: f64, right_travel: f64) -> f64 {
        self.stiffness * (left_travel - right_travel)
    }
}
/// Dynamic state for the quarter-car 2-DOF model.
///
/// Positions are measured from the static equilibrium position.
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
    /// Road surface height under the tyre, relative to the equilibrium road
    /// level (m, positive = road raised).  It is integrated from the road
    /// velocity passed to [`QuarterCarModel::step`].
    pub zr: f64,
}
impl QuarterCarState {
    /// Create state at rest at the static equilibrium position on a flat road
    /// (`zr = 0`).
    pub fn at_equilibrium(model: &QuarterCarModel, _g: f64) -> Self {
        let _ = model;
        Self {
            xs: 0.0,
            xu: 0.0,
            vs: 0.0,
            vu: 0.0,
            zr: 0.0,
        }
    }
}
/// McPherson strut spring-damper unit (spring/damper focus).
///
/// Captures the spring rate, damping, rest length, caster angle, and
/// kingpin inclination as scalar parameters.
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
    pub fn pitch_moment(&self, long_accel: f64) -> f64 {
        self.mass * long_accel * self.h_cg
    }
    /// Roll moment (N·m) due to lateral acceleration.
    ///
    /// `M_roll = m * a_lat * h_cg`
    pub fn roll_moment(&self, lat_accel: f64) -> f64 {
        self.mass * lat_accel.abs() * self.h_cg
    }
    /// Front/rear load transfer due to pitch (longitudinal acceleration).
    ///
    /// Returns `(delta_front_N, delta_rear_N)`.
    ///
    /// - Front axle gains load during braking (negative `long_accel`).
    /// - Rear axle gains load during acceleration.
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
    pub fn left_right_load_transfer_roll(&self, lat_accel: f64) -> (f64, f64) {
        if self.track_f <= 0.0 {
            return (0.0, 0.0);
        }
        let m_roll = self.mass * lat_accel * self.h_cg;
        let delta = m_roll / self.track_f;
        (-delta, delta)
    }
}
/// Relative tolerance below which two directions are treated as parallel and
/// below which two points are treated as coincident in the front-view solver.
const PARALLEL_TOLERANCE: f64 = 1e-12;
/// Relative tolerance of the circle-intersection tests of the front-view
/// solver; it absorbs rounding at the boundary of the reachable range.
const CIRCLE_TOLERANCE: f64 = 1e-12;
/// Half-length of the chord at distance `offset` from the centre of a circle
/// of radius `radius`: `sqrt(radius² - offset²)`.
///
/// Returns `None` when the line lies outside the circle.  A line that misses by
/// no more than rounding noise is treated as tangent.
fn circle_half_chord(radius: f64, offset: f64) -> Option<f64> {
    let squared = radius * radius - offset * offset;
    if squared >= 0.0 {
        Some(squared.sqrt())
    } else if -squared <= CIRCLE_TOLERANCE * radius * radius {
        Some(0.0)
    } else {
        None
    }
}
/// Front-view configuration of a [`DoubleWishboneKinematic`] linkage at one
/// wheel travel, returned by [`DoubleWishboneKinematic::pose_at_travel`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DoubleWishbonePose {
    /// Upper ball joint `[x, y, z]` (m).
    pub upper_outboard: [f64; 3],
    /// Lower ball joint `[x, y, z]` (m).
    pub lower_outboard: [f64; 3],
    /// Upright inclination from the vertical (degrees); positive when the top
    /// of the upright leans outboard.
    pub camber_deg: f64,
}
/// Full double-wishbone suspension kinematic model.
///
/// Tracks upper and lower control arms, their inboard pivot points, and the
/// wheel spindle location.
#[derive(Debug, Clone)]
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
impl DoubleWishboneKinematic {
    /// Construct from geometry and compliance data.
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
    /// Sign of the outboard direction along the lateral `y` axis.
    ///
    /// `+1.0` when the control arms reach towards `+y` (the usual right-hand
    /// side of the car), `-1.0` for a mirrored (left-hand) layout.  It makes the
    /// camber sign convention independent of the side of the car.
    fn outboard_sign(&self) -> f64 {
        let reach = (self.upper_outboard[1] - self.upper_inboard[1])
            + (self.lower_outboard[1] - self.lower_inboard[1]);
        if reach < 0.0 { -1.0 } else { 1.0 }
    }
    /// Lengths `(upper arm, lower arm, upright)` of the three moving links,
    /// measured in the front-view (`YZ`) plane.
    ///
    /// The mechanism moves in the front-view plane with every `x` coordinate
    /// held fixed, so the effective link lengths are the `YZ` projections of
    /// the 3D members; for a member with a fore-aft offset this is the radius
    /// of the circle its outboard end sweeps in the plane `x = const`.
    fn front_view_link_lengths(&self) -> (f64, f64, f64) {
        let projected = |p: [f64; 3], q: [f64; 3]| (q[1] - p[1]).hypot(q[2] - p[2]);
        (
            projected(self.upper_inboard, self.upper_outboard),
            projected(self.lower_inboard, self.lower_outboard),
            projected(self.lower_outboard, self.upper_outboard),
        )
    }
    /// Front-view instant centre `[y, z]` (m): the intersection of the lines
    /// through the upper and the lower control arm.
    ///
    /// Returns `None` when the arms are parallel (the instant centre is then at
    /// infinity) or when an arm has zero length.
    pub fn instant_centre(&self) -> Option<[f64; 2]> {
        let uy = self.upper_outboard[1] - self.upper_inboard[1];
        let uz = self.upper_outboard[2] - self.upper_inboard[2];
        let ly = self.lower_outboard[1] - self.lower_inboard[1];
        let lz = self.lower_outboard[2] - self.lower_inboard[2];
        let det = uy * lz - uz * ly;
        let scale = uy.hypot(uz) * ly.hypot(lz);
        if scale == 0.0 || det.abs() <= PARALLEL_TOLERANCE * scale {
            return None;
        }
        let dy = self.lower_inboard[1] - self.upper_inboard[1];
        let dz = self.lower_inboard[2] - self.upper_inboard[2];
        let t = (dy * lz - dz * ly) / det;
        Some([
            self.upper_inboard[1] + t * uy,
            self.upper_inboard[2] + t * uz,
        ])
    }
    /// Instant centre height (m) in the suspension coordinate frame, whose
    /// ground plane is `z = 0`.
    ///
    /// It is the height of [`instant_centre`](Self::instant_centre).  For
    /// parallel arms the instant centre is at infinity and `0.0` is returned;
    /// use [`roll_centre_height`](Self::roll_centre_height), which handles that
    /// case exactly, when the instant centre is only a means to an end.
    pub fn instant_centre_height(&self) -> f64 {
        self.instant_centre().map_or(0.0, |centre| centre[1])
    }
    /// Roll centre height (m) at the design position, by the instant-centre
    /// method.
    ///
    /// The ground plane is `z = 0` and the contact patch is taken directly
    /// below the wheel centre, at `(spindle y, 0)`.  The roll centre is where
    /// the line from the contact patch through the instant centre crosses the
    /// vehicle centre plane `y = 0`.  With parallel control arms the instant
    /// centre is at infinity along the arm direction, so the line runs parallel
    /// to the arms: horizontal parallel arms put the roll centre on the ground,
    /// arms rising outboard put it below the ground.
    ///
    /// Returns `NaN` when the roll centre is undefined: the instant centre
    /// lies vertically above the contact patch, the parallel arms are vertical,
    /// or both arms have zero length.
    pub fn roll_centre_height(&self) -> f64 {
        let contact_y = self.spindle[1];
        match self.instant_centre() {
            Some([centre_y, centre_z]) => {
                let run = contact_y - centre_y;
                let reference = contact_y.abs().max(centre_y.abs());
                if run.abs() <= PARALLEL_TOLERANCE * reference {
                    return f64::NAN;
                }
                // Adding 0.0 turns a negative zero into +0.0.
                contact_y * centre_z / run + 0.0
            }
            None => {
                let upper = (
                    self.upper_outboard[1] - self.upper_inboard[1],
                    self.upper_outboard[2] - self.upper_inboard[2],
                );
                let lower = (
                    self.lower_outboard[1] - self.lower_inboard[1],
                    self.lower_outboard[2] - self.lower_inboard[2],
                );
                let (run, rise) = if upper.0.hypot(upper.1) > 0.0 {
                    upper
                } else {
                    lower
                };
                if run.abs() <= PARALLEL_TOLERANCE * run.hypot(rise) {
                    return f64::NAN;
                }
                -contact_y * rise / run + 0.0
            }
        }
    }
    /// Camber change per unit of wheel travel (deg/mm) at the design position.
    ///
    /// Central finite difference of [`camber_at_travel`](Self::camber_at_travel)
    /// over ±1 mm of lower-ball-joint travel; a one-sided difference is used
    /// when the linkage cannot move 1 mm in one direction.  `NaN` if the
    /// linkage cannot move in either direction.
    ///
    /// To first order the upright turns about the instant centre, so the gain
    /// is `-s / (y_lower - y_ic)` radians per metre (`s` = `+1` for a
    /// right-hand layout, `-1` for a mirrored one): negative, i.e. the wheel
    /// gains negative camber in bump, when the instant centre lies inboard of
    /// the lower ball joint, positive when it lies outboard, and zero when the
    /// arms are parallel (instant centre at infinity).  The length ratio of the
    /// arms only matters through the instant centre it produces: horizontal
    /// parallel arms of unequal length still have zero gain at the design
    /// position and only curve at second order.
    pub fn camber_gain(&self) -> f64 {
        let step = 0.001;
        let per_mm = |delta_deg: f64, travel_m: f64| delta_deg / (travel_m * 1000.0);
        let bump = self.camber_at_travel(step);
        let droop = self.camber_at_travel(-step);
        match (bump.is_finite(), droop.is_finite()) {
            (true, true) => per_mm(bump - droop, 2.0 * step),
            (true, false) => per_mm(bump - self.camber_at_travel(0.0), step),
            (false, true) => per_mm(self.camber_at_travel(0.0) - droop, step),
            (false, false) => f64::NAN,
        }
    }
    /// Camber angle (degrees) at a given wheel travel.
    ///
    /// `delta_z` is the vertical displacement of the lower ball joint from its
    /// design position (m, positive = bump).  The camber is the inclination of
    /// the upright (the line from the lower to the upper ball joint) from the
    /// vertical in the front view, positive when the top leans outboard.  At
    /// `delta_z = 0` it equals the static inclination of the given geometry;
    /// the change from that value is the camber change of the wheel.
    ///
    /// Returns `NaN` if the linkage cannot reach the requested travel; see
    /// [`pose_at_travel`](Self::pose_at_travel).
    pub fn camber_at_travel(&self, delta_z: f64) -> f64 {
        self.pose_at_travel(delta_z)
            .map_or(f64::NAN, |pose| pose.camber_deg)
    }
    /// Solve the front-view four-bar linkage at a given wheel travel.
    ///
    /// The body is the fixed link between the two inboard pivots, the control
    /// arms are rigid links of constant (front-view) length and the upright is
    /// the coupler between the two ball joints.  The lower ball joint is moved
    /// to `z + delta_z` on the circle swept by the lower arm about its inboard
    /// pivot; the upper ball joint is the intersection of the circle swept by
    /// the upper arm and the circle of the upright's length around the lower
    /// ball joint.  Of the two intersections the one that keeps the assembly
    /// mode of the design position (the same side of the line from the upper
    /// inboard pivot to the lower ball joint) is taken, which is the continuous
    /// continuation of the motion.  If the design position itself is a dead
    /// centre (upper arm and upright collinear) both modes coincide there and
    /// the root nearest to the design position of the upper ball joint is used.
    ///
    /// `delta_z = 0` returns the design position exactly.  Returns `None` for a
    /// non-finite `delta_z`, when the lower arm cannot reach the height, or when
    /// the upper arm and the upright cannot close the loop (no real
    /// intersection).
    pub fn pose_at_travel(&self, delta_z: f64) -> Option<DoubleWishbonePose> {
        if !delta_z.is_finite() {
            return None;
        }
        let sign = self.outboard_sign();
        let camber_of = |upper: [f64; 3], lower: [f64; 3]| {
            (sign * (upper[1] - lower[1]))
                .atan2(upper[2] - lower[2])
                .to_degrees()
        };
        if delta_z == 0.0 {
            return Some(DoubleWishbonePose {
                upper_outboard: self.upper_outboard,
                lower_outboard: self.lower_outboard,
                camber_deg: camber_of(self.upper_outboard, self.lower_outboard),
            });
        }
        let (upper_len, lower_len, upright_len) = self.front_view_link_lengths();
        let upper_pivot = [self.upper_inboard[1], self.upper_inboard[2]];
        let lower_pivot = [self.lower_inboard[1], self.lower_inboard[2]];
        let upper_ref = [self.upper_outboard[1], self.upper_outboard[2]];
        let lower_ref = [self.lower_outboard[1], self.lower_outboard[2]];

        // Lower ball joint on its circle at the requested height; it stays on the
        // side of the inboard pivot where the design position has it.
        let lower_z = lower_ref[1] + delta_z;
        let rise = lower_z - lower_pivot[1];
        let lower_reach = circle_half_chord(lower_len, rise)?;
        let lower_side = if lower_ref[0] < lower_pivot[0] {
            -1.0
        } else {
            1.0
        };
        let lower = [lower_pivot[0] + lower_side * lower_reach, lower_z];

        // Upper ball joint: intersection of the circle (upper pivot, upper arm)
        // and the circle (lower ball joint, upright).
        let chord = [lower[0] - upper_pivot[0], lower[1] - upper_pivot[1]];
        let distance = chord[0].hypot(chord[1]);
        let scale = upper_len.max(upright_len).max(distance);
        if distance <= CIRCLE_TOLERANCE * scale
            || distance > upper_len + upright_len + CIRCLE_TOLERANCE * scale
            || distance < (upper_len - upright_len).abs() - CIRCLE_TOLERANCE * scale
        {
            return None;
        }
        let along = (upper_len * upper_len - upright_len * upright_len + distance * distance)
            / (2.0 * distance);
        let offset = circle_half_chord(upper_len, along)?;
        let unit = [chord[0] / distance, chord[1] / distance];
        let normal = [-unit[1], unit[0]];
        let foot = [
            upper_pivot[0] + along * unit[0],
            upper_pivot[1] + along * unit[1],
        ];
        let plus = [foot[0] + offset * normal[0], foot[1] + offset * normal[1]];
        let minus = [foot[0] - offset * normal[0], foot[1] - offset * normal[1]];

        // Assembly mode of the design position: the sign of the cross product
        // (lower ball joint - upper pivot) x (upper ball joint - upper pivot).
        let ref_chord = [lower_ref[0] - upper_pivot[0], lower_ref[1] - upper_pivot[1]];
        let ref_arm = [upper_ref[0] - upper_pivot[0], upper_ref[1] - upper_pivot[1]];
        let ref_cross = ref_chord[0] * ref_arm[1] - ref_chord[1] * ref_arm[0];
        let ref_scale = ref_chord[0].hypot(ref_chord[1]) * upper_len;
        let upper = if ref_cross > CIRCLE_TOLERANCE * ref_scale {
            plus
        } else if ref_cross < -CIRCLE_TOLERANCE * ref_scale {
            minus
        } else {
            // The design position is a dead centre; fall back to the root nearest
            // to the design position of the upper ball joint.
            let gap = |p: [f64; 2]| (p[0] - upper_ref[0]).hypot(p[1] - upper_ref[1]);
            if gap(plus) <= gap(minus) { plus } else { minus }
        };

        let upper_outboard = [self.upper_outboard[0], upper[0], upper[1]];
        let lower_outboard = [self.lower_outboard[0], lower[0], lower[1]];
        Some(DoubleWishbonePose {
            upper_outboard,
            lower_outboard,
            camber_deg: camber_of(upper_outboard, lower_outboard),
        })
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
#[derive(Debug, Clone)]
pub struct BumpStopLinear {
    /// Pre-load force present at zero compression (N).
    pub preload: f64,
    /// Linear stiffness rate (N/m).
    pub rate: f64,
    /// Maximum force clamp (N).
    pub max_force: f64,
}
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
