//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)
use oxiphysics_core::math::Real;

/// Extended steer-by-wire simulator with column torque estimation and
/// feedback force rendering.
#[derive(Debug, Clone)]
pub struct SteerByWireSimulator {
    /// Inner SBW servo.
    pub servo: SteerByWire,
    /// Steering feel model.
    pub feel: SteeringFeel,
    /// Column compliance.
    pub column: SteeringColumnCompliance,
    /// Current road wheel angle (rad).
    pub road_wheel_angle: Real,
}
impl SteerByWireSimulator {
    /// Create a new steer-by-wire simulator.
    pub fn new(servo: SteerByWire, feel: SteeringFeel, column: SteeringColumnCompliance) -> Self {
        Self {
            road_wheel_angle: 0.0,
            servo,
            feel,
            column,
        }
    }
    /// Update the simulator one step.
    ///
    /// * `driver_input` — target steering angle from driver (rad).
    /// * `lateral_force` — tire lateral force (N).
    /// * `slip_angle` — tire slip angle (rad).
    /// * `dt` — time step (s).
    ///
    /// Returns `(road_wheel_angle, feedback_torque)`.
    pub fn step(
        &mut self,
        driver_input: Real,
        lateral_force: Real,
        slip_angle: Real,
        dt: Real,
    ) -> (Real, Real) {
        self.servo.target_angle = driver_input;
        self.road_wheel_angle = self.servo.update(self.road_wheel_angle, dt);
        let aligning = self.feel.aligning_torque(lateral_force, slip_angle);
        let rack_torque = aligning;
        let pinion = self
            .column
            .step(0.0, rack_torque, self.road_wheel_angle, dt);
        self.road_wheel_angle = pinion;
        let feedback = aligning.clamp(-self.servo.torque_limit, self.servo.torque_limit);
        (self.road_wheel_angle, feedback)
    }
    /// Column torque estimation: the torque felt by the driver through the wheel.
    pub fn column_torque(&self, lateral_force: Real, slip_angle: Real) -> Real {
        let aligning = self.feel.aligning_torque(lateral_force, slip_angle);
        (aligning / self.column.gear_ratio.max(1e-9))
            .clamp(-self.servo.torque_limit, self.servo.torque_limit)
    }
}
/// Steering column torsional compliance model.
///
/// Models the twist in the steering column shaft between the steering wheel
/// and the rack input.  The compliance introduces a phase lag in steering
/// response and affects the feedback feel.
#[derive(Debug, Clone)]
pub struct SteeringColumnCompliance {
    /// Torsional stiffness (N·m/rad).
    pub torsional_stiffness: Real,
    /// Torsional damping coefficient (N·m·s/rad).
    pub torsional_damping: Real,
    /// Gear ratio: steering wheel turns per pinion turn.
    pub gear_ratio: Real,
    /// Current column twist angle (radians, state variable).
    pub(super) twist_angle: Real,
    /// Current twist rate (rad/s, state variable).
    pub(super) twist_rate: Real,
}
impl SteeringColumnCompliance {
    /// Create a new compliance model.
    pub fn new(torsional_stiffness: Real, torsional_damping: Real, gear_ratio: Real) -> Self {
        Self {
            torsional_stiffness,
            torsional_damping,
            gear_ratio,
            twist_angle: 0.0,
            twist_rate: 0.0,
        }
    }
    /// Advance the compliance model one time step.
    ///
    /// * `driver_torque` — torque applied by driver at steering wheel (N·m).
    /// * `rack_torque` — reaction torque from the rack (N·m, referred to pinion).
    /// * `dt` — time step (s).
    ///
    /// Returns the updated pinion angle (radians).
    pub fn step(
        &mut self,
        driver_torque: Real,
        rack_torque: Real,
        pinion_angle_prev: Real,
        dt: Real,
    ) -> Real {
        let spring_torque = self.torsional_stiffness * self.twist_angle;
        let damping_torque = self.torsional_damping * self.twist_rate;
        let rack_referred = rack_torque / self.gear_ratio.max(1e-9);
        let net = driver_torque - rack_referred;
        let target_twist = net / self.torsional_stiffness.max(1e-9);
        let tau = self.torsional_damping / self.torsional_stiffness.max(1e-9);
        let alpha = dt / (tau + dt).max(dt);
        self.twist_angle += alpha * (target_twist - self.twist_angle);
        let _ = spring_torque;
        let _ = damping_torque;
        self.twist_rate = (target_twist - self.twist_angle) / dt.max(1e-9);
        pinion_angle_prev + self.twist_angle / self.gear_ratio.max(1e-9)
    }
    /// Return current twist angle (radians).
    pub fn twist_angle(&self) -> Real {
        self.twist_angle
    }
    /// Reset column state.
    pub fn reset(&mut self) {
        self.twist_angle = 0.0;
        self.twist_rate = 0.0;
    }
    /// Torsional natural frequency (Hz).
    pub fn natural_frequency(&self) -> Real {
        let omega_n = (self.torsional_stiffness / 1.0_f64).sqrt();
        omega_n / (2.0 * std::f64::consts::PI)
    }
}
/// Rack-and-pinion steering mechanism model.
///
/// Converts rotational pinion motion to linear rack travel, and then to
/// steering angle via the tie-rod geometry.
///
/// The fundamental relationship is:
///   rack_travel = pinion_angle * pinion_radius
///   steer_angle ≈ atan(rack_travel / steer_arm_length)
#[derive(Debug, Clone)]
pub struct RackAndPinion {
    /// Pinion pitch radius (meters).
    pub pinion_radius: Real,
    /// Effective steering arm length from rack attachment to wheel pivot (meters).
    pub steer_arm_length: Real,
    /// Maximum rack travel (meters, half total travel).
    pub max_rack_travel: Real,
    /// Rack compliance: lateral deflection per unit force (m / N).
    pub rack_compliance: Real,
}
impl RackAndPinion {
    /// Create a new rack-and-pinion model.
    pub fn new(
        pinion_radius: Real,
        steer_arm_length: Real,
        max_rack_travel: Real,
        rack_compliance: Real,
    ) -> Self {
        Self {
            pinion_radius,
            steer_arm_length,
            max_rack_travel,
            rack_compliance,
        }
    }
    /// Rack travel (m) for a given pinion rotation angle (radians).
    pub fn rack_travel(&self, pinion_angle: Real) -> Real {
        let travel = pinion_angle * self.pinion_radius;
        travel.clamp(-self.max_rack_travel, self.max_rack_travel)
    }
    /// Steering angle (radians) from rack travel (meters).
    pub fn steer_angle_from_travel(&self, rack_travel: Real) -> Real {
        (rack_travel / self.steer_arm_length.max(1e-9)).atan()
    }
    /// Overall steer angle (radians) from pinion angle (radians).
    pub fn steer_angle(&self, pinion_angle: Real) -> Real {
        let travel = self.rack_travel(pinion_angle);
        self.steer_angle_from_travel(travel)
    }
    /// Rack-to-wheel force ratio (mechanical advantage of the steering arm).
    ///
    /// Returns the ratio: lateral force at wheel / rack force.
    pub fn mechanical_advantage(&self, rack_travel: Real) -> Real {
        let ratio = rack_travel / self.steer_arm_length.max(1e-9);
        1.0 / (self.steer_arm_length * (1.0 + ratio * ratio)).max(1e-12)
    }
    /// Lateral force on the rack for a given wheel aligning torque (N·m)
    /// and current rack travel (m).
    ///
    /// Returns the rack force in Newtons.
    pub fn rack_force(&self, aligning_torque: Real, rack_travel: Real) -> Real {
        let ma = self.mechanical_advantage(rack_travel);
        if ma.abs() < 1e-12 {
            return 0.0;
        }
        aligning_torque / (ma * self.steer_arm_length).max(1e-12)
    }
    /// Pinion torque required to overcome rack force `rack_force` (N).
    pub fn required_pinion_torque(&self, rack_force: Real) -> Real {
        rack_force * self.pinion_radius
    }
    /// Deflection of rack under lateral force `force` (N).
    pub fn rack_deflection(&self, force: Real) -> Real {
        force * self.rack_compliance
    }
    /// Total rack travel (meters) for a given steering wheel angle (radians)
    /// with `steering_ratio` = steering_wheel_angle / pinion_angle.
    pub fn rack_travel_from_wheel(&self, steering_wheel_angle: Real, steering_ratio: Real) -> Real {
        let pinion_angle = steering_wheel_angle / steering_ratio.max(1e-9);
        self.rack_travel(pinion_angle)
    }
}
/// Real-time oversteer/understeer detector.
///
/// Compares the measured yaw rate to the expected yaw rate from the bicycle
/// model.  The difference indicates understeer (positive) or oversteer (negative).
#[derive(Debug, Clone)]
pub struct OversteerUndersteerDetector {
    /// Wheelbase (m).
    pub wheelbase: Real,
    /// Dead-band for neutral classification (rad/s).
    pub dead_band: Real,
}
impl OversteerUndersteerDetector {
    /// Create a new detector.
    pub fn new(wheelbase: Real, dead_band: Real) -> Self {
        Self {
            wheelbase,
            dead_band,
        }
    }
    /// Expected yaw rate from bicycle model: ψ̇_ref = v * tan(δ) / L.
    pub fn expected_yaw_rate(&self, speed: Real, steer_angle: Real) -> Real {
        speed * steer_angle.tan() / self.wheelbase.max(1e-6)
    }
    /// Yaw rate deviation: measured - expected.
    ///
    /// Positive = vehicle yaws less than expected (understeer).
    /// Negative = vehicle yaws more than expected (oversteer).
    pub fn yaw_rate_deviation(&self, measured: Real, speed: Real, steer_angle: Real) -> Real {
        self.expected_yaw_rate(speed, steer_angle) - measured
    }
    /// Classify the current handling state.
    ///
    /// Returns `"understeer"`, `"oversteer"`, or `"neutral"`.
    pub fn classify(&self, measured: Real, speed: Real, steer_angle: Real) -> &'static str {
        let dev = self.yaw_rate_deviation(measured, speed, steer_angle);
        if dev > self.dead_band {
            "understeer"
        } else if dev < -self.dead_band {
            "oversteer"
        } else {
            "neutral"
        }
    }
    /// Normalised oversteer index: 0 = neutral, +1 = full understeer, -1 = full oversteer.
    pub fn oversteer_index(&self, measured: Real, speed: Real, steer_angle: Real) -> Real {
        let ref_rate = self.expected_yaw_rate(speed, steer_angle);
        if ref_rate.abs() < 1e-6 {
            return 0.0;
        }
        (ref_rate - measured) / ref_rate.abs()
    }
}
/// Pure Ackermann geometry helper (no blending factor).
///
/// Given wheelbase and track_width, computes ideal inner/outer angles
/// and turning radius.
#[derive(Debug, Clone)]
pub struct AckermannGeometry {
    /// Distance between front and rear axles (meters).
    pub wheelbase: Real,
    /// Distance between left and right wheels on the same axle (meters).
    pub track_width: Real,
}
impl AckermannGeometry {
    /// Create a new geometry helper.
    pub fn new(wheelbase: Real, track_width: Real) -> Self {
        Self {
            wheelbase,
            track_width,
        }
    }
    /// Compute inner and outer wheel steering angles from a mid-axle angle.
    ///
    /// Returns `(inner_angle, outer_angle)` in radians.
    /// The inner wheel has the larger angle.
    pub fn steering_angles(&self, steering_input: Real) -> (Real, Real) {
        if steering_input.abs() < 1e-10 {
            return (0.0, 0.0);
        }
        let sign = steering_input.signum();
        let base = steering_input.abs();
        let r = self.wheelbase / base.tan();
        let inner = (self.wheelbase / (r - self.track_width * 0.5)).atan();
        let outer = (self.wheelbase / (r + self.track_width * 0.5)).atan();
        (sign * inner, sign * outer)
    }
    /// Compute turning radius for a given (average) steer angle.
    ///
    /// Returns `f64::INFINITY` when angle is near zero.
    pub fn turning_radius(&self, steer_angle: Real) -> Real {
        if steer_angle.abs() < 1e-10 {
            return f64::INFINITY;
        }
        (self.wheelbase / steer_angle.tan()).abs()
    }
    /// Maximum steering angle derivable from geometry (inner wheel limit
    /// when the inner turning circle radius equals track_width / 2).
    ///
    /// Returns the angle in radians.
    pub fn max_steering_angle(&self) -> Real {
        let r_min = self.track_width * 0.5 + 0.01;
        (self.wheelbase / r_min).atan()
    }
}
/// Torque steer compensation model.
///
/// Unequal drive shaft lengths in FWD vehicles cause steering pull under
/// power. This model computes a corrective steering offset.
#[derive(Debug, Clone)]
pub struct TorqueSteerCompensation {
    /// Sensitivity: angular offset per unit drive torque (rad / N*m).
    pub sensitivity: Real,
    /// Maximum correction angle (radians).
    pub max_correction: Real,
    /// Asymmetry factor: positive means pull to the right under power.
    pub asymmetry: Real,
}
impl TorqueSteerCompensation {
    /// Create a new torque steer compensation model.
    pub fn new(sensitivity: Real, max_correction: Real, asymmetry: Real) -> Self {
        Self {
            sensitivity,
            max_correction,
            asymmetry,
        }
    }
    /// Compute the corrective steering offset to counteract torque steer.
    ///
    /// `drive_torque` is the total drive torque delivered to the front wheels (N*m).
    pub fn correction(&self, drive_torque: Real) -> Real {
        let raw = -drive_torque * self.sensitivity * self.asymmetry;
        raw.clamp(-self.max_correction, self.max_correction)
    }
    /// Apply the correction to a base steering angle.
    pub fn apply(&self, base_angle: Real, drive_torque: Real) -> Real {
        base_angle + self.correction(drive_torque)
    }
}
/// Rear-wheel-only steering model.
///
/// In some forklifts and specialized vehicles, only the rear wheels steer.
/// This model computes the rear angle directly from the driver input.
#[derive(Debug, Clone)]
pub struct RearWheelSteering {
    /// Maximum rear steering angle (radians).
    pub max_angle: Real,
    /// Steering gain (ratio of rear angle to input).
    pub gain: Real,
}
impl RearWheelSteering {
    /// Create a new rear-wheel steering model.
    pub fn new(max_angle: Real, gain: Real) -> Self {
        Self { max_angle, gain }
    }
    /// Compute the rear steering angle from a normalized input \[-1, 1\].
    pub fn compute_rear_angle(&self, input: Real) -> Real {
        let raw = input.clamp(-1.0, 1.0) * self.gain * self.max_angle;
        raw.clamp(-self.max_angle, self.max_angle)
    }
    /// Compute the turning radius for rear-wheel steering.
    ///
    /// `wheelbase` is the axle-to-axle distance (m).
    pub fn turning_radius(&self, input: Real, wheelbase: Real) -> Option<Real> {
        let angle = self.compute_rear_angle(input);
        if angle.abs() < 1e-10 {
            None
        } else {
            Some((wheelbase / angle.tan()).abs())
        }
    }
}
/// Front/rear steering correlation model for 4WS vehicles.
///
/// Tracks the correlation between front and rear steering angles and
/// detects anomalies (e.g. sensor faults).
#[derive(Debug, Clone)]
pub struct FrontRearSteeringCorrelation {
    /// Expected ratio of rear angle to front angle at low speed.
    pub low_speed_ratio: Real,
    /// Expected ratio of rear angle to front angle at high speed.
    pub high_speed_ratio: Real,
    /// Speed threshold (m/s).
    pub speed_threshold: Real,
    /// Tolerance for fault detection.
    pub fault_tolerance: Real,
}
impl FrontRearSteeringCorrelation {
    /// Create a new correlation model.
    pub fn new(
        low_speed_ratio: Real,
        high_speed_ratio: Real,
        speed_threshold: Real,
        fault_tolerance: Real,
    ) -> Self {
        Self {
            low_speed_ratio,
            high_speed_ratio,
            speed_threshold,
            fault_tolerance,
        }
    }
    /// Expected rear angle for a given front angle and speed.
    pub fn expected_rear(&self, front_angle: Real, speed: Real) -> Real {
        let t = (speed / self.speed_threshold.max(1e-6)).clamp(0.0, 1.0);
        let ratio = self.low_speed_ratio + t * (self.high_speed_ratio - self.low_speed_ratio);
        front_angle * ratio
    }
    /// Check if the actual rear angle is within tolerance of the expected value.
    pub fn is_correlated(&self, front_angle: Real, rear_angle: Real, speed: Real) -> bool {
        let expected = self.expected_rear(front_angle, speed);
        (rear_angle - expected).abs() <= self.fault_tolerance
    }
    /// Correlation error: actual - expected.
    pub fn correlation_error(&self, front_angle: Real, rear_angle: Real, speed: Real) -> Real {
        rear_angle - self.expected_rear(front_angle, speed)
    }
}
/// Speed-dependent steering ratio that varies the mechanical advantage
/// between the steering wheel and road wheels.
///
/// At low speed, the ratio is `low_speed_ratio` (quick steering).
/// At high speed, the ratio is `high_speed_ratio` (slower, more stable).
#[derive(Debug, Clone)]
pub struct SpeedDependentRatio {
    /// Steering ratio at zero speed.
    pub low_speed_ratio: Real,
    /// Steering ratio at `transition_speed`.
    pub high_speed_ratio: Real,
    /// Speed (m/s) at which the ratio reaches the high-speed value.
    pub transition_speed: Real,
}
impl SpeedDependentRatio {
    /// Create a new speed-dependent ratio.
    pub fn new(low_ratio: Real, high_ratio: Real, transition_speed: Real) -> Self {
        Self {
            low_speed_ratio: low_ratio,
            high_speed_ratio: high_ratio,
            transition_speed,
        }
    }
    /// Compute the steering ratio at a given speed.
    pub fn ratio_at_speed(&self, speed: Real) -> Real {
        let t = (speed / self.transition_speed.max(1e-6)).clamp(0.0, 1.0);
        self.low_speed_ratio + t * (self.high_speed_ratio - self.low_speed_ratio)
    }
    /// Apply the speed-dependent ratio to a steering input.
    ///
    /// Returns the effective road wheel angle for a given steering-wheel angle
    /// and vehicle speed.
    pub fn apply(&self, steering_wheel_angle: Real, speed: Real) -> Real {
        let ratio = self.ratio_at_speed(speed);
        if ratio.abs() < 1e-10 {
            return 0.0;
        }
        steering_wheel_angle / ratio
    }
}
/// Steer-by-wire actuator model with servo rate limiting and torque feedback.
#[derive(Debug, Clone)]
pub struct SteerByWire {
    /// Desired target angle in radians.
    pub target_angle: Real,
    /// Maximum rate the servo can move (rad/s).
    pub servo_rate: Real,
    /// Maximum torque the servo can apply (N·m).
    pub torque_limit: Real,
}
impl SteerByWire {
    /// Create a new steer-by-wire actuator.
    pub fn new(target_angle: Real, servo_rate: Real, torque_limit: Real) -> Self {
        Self {
            target_angle,
            servo_rate,
            torque_limit,
        }
    }
    /// Step the servo toward `target_angle` with rate limiting.
    ///
    /// Returns the new actuated angle after `dt` seconds.
    pub fn update(&self, current_angle: Real, dt: Real) -> Real {
        let delta = self.target_angle - current_angle;
        let max_delta = self.servo_rate * dt;
        current_angle + delta.clamp(-max_delta, max_delta)
    }
    /// Compute the driver-feedback torque given road self-aligning torque.
    ///
    /// Returns a signed torque value clamped to `torque_limit`.
    pub fn apply_torque_feedback(&self, road_torque: Real, driver_torque: Real) -> Real {
        let combined = road_torque + driver_torque;
        combined.clamp(-self.torque_limit, self.torque_limit)
    }
}
/// A complete 4-wheel steering system combining Ackermann front geometry
/// with speed-dependent rear steering.
#[derive(Debug, Clone)]
pub struct FullFourWheelSteering {
    /// Front Ackermann geometry.
    pub front: AckermannSteering,
    /// Four-wheel steering model.
    pub rear: FourWheelSteering,
    /// Speed-dependent ratio.
    pub ratio: SpeedDependentRatio,
}
impl FullFourWheelSteering {
    /// Create a full 4WS system.
    pub fn new(
        front: AckermannSteering,
        rear: FourWheelSteering,
        ratio: SpeedDependentRatio,
    ) -> Self {
        Self { front, rear, ratio }
    }
    /// Compute all four wheel angles.
    ///
    /// Returns `(front_left, front_right, rear_left, rear_right)`.
    pub fn compute_all_angles(
        &self,
        steering_wheel_angle: Real,
        speed: Real,
    ) -> (Real, Real, Real, Real) {
        let road_angle = self.ratio.apply(steering_wheel_angle, speed);
        let input = (road_angle / self.front.max_steer_angle).clamp(-1.0, 1.0);
        let (fl, fr) = self.front.compute_wheel_angles(input);
        let front_avg = (fl + fr) * 0.5;
        let rear_angle = self.rear.rear_angle(front_avg, speed);
        (fl, fr, rear_angle, rear_angle)
    }
}
/// Enhanced steering feel model with on-centre weighting and speed scaling.
///
/// Adds a velocity-dependent damping gradient and on-centre notch.
#[derive(Debug, Clone)]
pub struct SteeringFeelEnhanced {
    /// Base steering feel model.
    pub base: SteeringFeel,
    /// On-centre notch width (radians).  Feel is lightened near centre.
    pub on_centre_width: Real,
    /// On-centre gain reduction factor (0 = no feel at centre, 1 = unchanged).
    pub on_centre_gain: Real,
    /// Speed scaling factor — feel increases with speed up to `speed_ref` (m/s).
    pub speed_ref: Real,
}
impl SteeringFeelEnhanced {
    /// Create a new enhanced steering feel model.
    pub fn new(
        base: SteeringFeel,
        on_centre_width: Real,
        on_centre_gain: Real,
        speed_ref: Real,
    ) -> Self {
        Self {
            base,
            on_centre_width,
            on_centre_gain: on_centre_gain.clamp(0.0, 1.0),
            speed_ref,
        }
    }
    /// Compute the enhanced aligning torque (N·m).
    pub fn aligning_torque(
        &self,
        lateral_force: Real,
        slip_angle: Real,
        steer_angle: Real,
        speed: Real,
    ) -> Real {
        let base_torque = self.base.aligning_torque(lateral_force, slip_angle);
        let oc_t = (steer_angle.abs() / self.on_centre_width.max(1e-9)).clamp(0.0, 1.0);
        let oc_gain = self.on_centre_gain + (1.0 - self.on_centre_gain) * oc_t;
        let speed_scale = (speed / self.speed_ref.max(1e-9)).clamp(0.0, 2.0);
        base_torque * oc_gain * speed_scale
    }
}
/// Steering feel model — computes self-aligning torque from caster, trail,
/// kingpin offset and lateral tire force.
#[derive(Debug, Clone)]
pub struct SteeringFeel {
    /// Caster angle in radians.
    pub caster: Real,
    /// Pneumatic trail in meters (distance behind contact patch centre).
    pub trail: Real,
    /// Kingpin offset (scrub radius) in meters.
    pub kingpin_offset: Real,
}
impl SteeringFeel {
    /// Create a new steering feel model.
    pub fn new(caster: Real, trail: Real, kingpin_offset: Real) -> Self {
        Self {
            caster,
            trail,
            kingpin_offset,
        }
    }
    /// Compute the self-aligning torque (N·m).
    ///
    /// Positive torque acts to return the wheel to centre.
    ///
    /// `lateral_force` — lateral tire force in Newtons.
    /// `slip_angle`    — tire slip angle in radians.
    pub fn aligning_torque(&self, lateral_force: Real, slip_angle: Real) -> Real {
        let effective_arm = self.trail + self.caster * slip_angle.sin().abs();
        -lateral_force * effective_arm - lateral_force * self.kingpin_offset * 0.1
    }
}
/// Dynamic (speed and angle dependent) steering ratio.
///
/// Extends `SpeedDependentRatio` by adding an angle-dependent component:
/// at large steering angles the effective ratio decreases (quicker response).
#[derive(Debug, Clone)]
pub struct DynamicSteeringRatio {
    /// Base speed-dependent ratio.
    pub base: SpeedDependentRatio,
    /// Maximum additional ratio reduction at full steer lock (dimensionless).
    pub angle_reduction: Real,
    /// Steering angle at which reduction is fully applied (radians).
    pub reduction_angle: Real,
}
impl DynamicSteeringRatio {
    /// Create a new dynamic steering ratio.
    pub fn new(base: SpeedDependentRatio, angle_reduction: Real, reduction_angle: Real) -> Self {
        Self {
            base,
            angle_reduction,
            reduction_angle,
        }
    }
    /// Compute the effective ratio at a given speed and steering wheel angle.
    pub fn ratio_at(&self, speed: Real, steering_wheel_angle: Real) -> Real {
        let base_ratio = self.base.ratio_at_speed(speed);
        let angle_t = (steering_wheel_angle.abs() / self.reduction_angle.max(1e-9)).clamp(0.0, 1.0);
        let reduction = angle_t * self.angle_reduction;
        (base_ratio - reduction).max(1.0)
    }
    /// Apply the dynamic ratio to convert steering wheel angle to road wheel angle.
    pub fn apply(&self, steering_wheel_angle: Real, speed: Real) -> Real {
        let ratio = self.ratio_at(speed, steering_wheel_angle);
        if ratio.abs() < 1e-9 {
            return 0.0;
        }
        steering_wheel_angle / ratio
    }
}
/// Ackermann steering geometry.
///
/// Computes inner and outer wheel steering angles so that all wheels
/// trace concentric arcs during a turn, reducing tire scrub.
///
/// The geometry is based on the bicycle model extended to two front wheels:
/// - Inner wheel turns more than the steering input
/// - Outer wheel turns less
#[derive(Debug, Clone)]
pub struct AckermannSteering {
    /// Distance between front and rear axles (meters).
    pub wheelbase: Real,
    /// Distance between left and right wheels on the same axle (meters).
    pub track_width: Real,
    /// Maximum steering angle in radians.
    pub max_steer_angle: Real,
    /// Ackermann factor: 0.0 = parallel, 1.0 = full Ackermann.
    pub ackermann_factor: Real,
}
impl AckermannSteering {
    /// Create a new Ackermann steering geometry.
    pub fn new(wheelbase: Real, track_width: Real, max_steer_angle: Real) -> Self {
        Self {
            wheelbase,
            track_width,
            max_steer_angle,
            ackermann_factor: 1.0,
        }
    }
    /// Create with a custom Ackermann factor.
    pub fn with_ackermann_factor(mut self, factor: Real) -> Self {
        self.ackermann_factor = factor.clamp(0.0, 1.0);
        self
    }
    /// Compute left and right wheel angles from a steering input.
    ///
    /// # Arguments
    /// * `steer_input` - Normalized steering input (-1..1), positive = turn right
    ///
    /// # Returns
    /// `(left_angle, right_angle)` in radians.
    /// Positive angle = wheel points right.
    pub fn compute_wheel_angles(&self, steer_input: Real) -> (Real, Real) {
        let steer_input = steer_input.clamp(-1.0, 1.0);
        let base_angle = steer_input * self.max_steer_angle;
        if base_angle.abs() < 1e-10 {
            return (0.0, 0.0);
        }
        let turn_radius = self.wheelbase / base_angle.tan();
        let inner_angle = (self.wheelbase / (turn_radius.abs() - self.track_width * 0.5)).atan();
        let outer_angle = (self.wheelbase / (turn_radius.abs() + self.track_width * 0.5)).atan();
        let parallel_inner = base_angle.abs();
        let parallel_outer = base_angle.abs();
        let inner_blended =
            parallel_inner * (1.0 - self.ackermann_factor) + inner_angle * self.ackermann_factor;
        let outer_blended =
            parallel_outer * (1.0 - self.ackermann_factor) + outer_angle * self.ackermann_factor;
        if base_angle > 0.0 {
            (outer_blended, inner_blended)
        } else {
            (-inner_blended, -outer_blended)
        }
    }
    /// Compute the turning radius for a given steering input.
    ///
    /// Returns `None` if steering is near zero (straight line).
    pub fn turning_radius(&self, steer_input: Real) -> Option<Real> {
        let base_angle = steer_input.clamp(-1.0, 1.0) * self.max_steer_angle;
        if base_angle.abs() < 1e-10 {
            None
        } else {
            Some((self.wheelbase / base_angle.tan()).abs())
        }
    }
}
/// Steering column model with dead zone, hysteresis and return-to-centre spring.
#[derive(Debug, Clone)]
pub struct SteeringColumn {
    /// Half-width of the input dead zone (dimensionless, 0..1).
    pub dead_zone: Real,
    /// Hysteresis band half-width (dimensionless).
    pub hysteresis: Real,
    /// Return-to-centre spring strength (applied as damping of offset, 0..1 per step).
    pub spring_strength: Real,
    /// Current internal angle state (tracks hysteresis history).
    pub(super) current_state: Real,
}
impl SteeringColumn {
    /// Create a new steering column.
    pub fn new(dead_zone: Real, hysteresis: Real, spring_strength: Real) -> Self {
        Self {
            dead_zone,
            hysteresis,
            spring_strength: spring_strength.clamp(0.0, 1.0),
            current_state: 0.0,
        }
    }
    /// Process a raw steering input through the column model.
    ///
    /// `raw_input` is in -1..1; `speed` modulates spring return strength.
    /// Returns the filtered output in -1..1.
    pub fn process_input(&mut self, raw_input: Real, speed: Real) -> Real {
        let after_dz = if raw_input.abs() < self.dead_zone {
            0.0
        } else {
            let sign = raw_input.signum();
            sign * (raw_input.abs() - self.dead_zone) / (1.0 - self.dead_zone).max(1e-6)
        };
        if (after_dz - self.current_state).abs() > self.hysteresis {
            self.current_state = after_dz - self.hysteresis * after_dz.signum();
        }
        let speed_factor = (speed / 30.0).clamp(0.0, 1.0);
        let spring = self.spring_strength * speed_factor;
        self.current_state *= 1.0 - spring;
        self.current_state
    }
}
/// Steering angle limiter with speed-dependent tightening.
///
/// At higher speeds the maximum allowed steering angle is reduced to
/// prevent excessive yaw rates.
#[derive(Debug, Clone)]
pub struct SteerLimiter {
    /// Maximum steering angle at rest (radians).
    pub max_angle: Real,
    /// Minimum steering angle at `speed_limit` (radians).
    pub min_angle: Real,
    /// Speed at which `min_angle` is reached (m/s).
    pub speed_limit: Real,
}
impl SteerLimiter {
    /// Create a new steering limiter.
    pub fn new(max_angle: Real, min_angle: Real, speed_limit: Real) -> Self {
        Self {
            max_angle,
            min_angle,
            speed_limit,
        }
    }
    /// Return the speed-limited maximum angle.
    pub fn effective_max(&self, speed: Real) -> Real {
        let t = (speed / self.speed_limit.max(1e-6)).clamp(0.0, 1.0);
        self.max_angle + t * (self.min_angle - self.max_angle)
    }
    /// Clamp `angle` to the speed-dependent limit.
    pub fn limit(&self, angle: Real, speed: Real) -> Real {
        let bound = self.effective_max(speed);
        angle.clamp(-bound, bound)
    }
}
/// Steering feedback force model that computes the force the driver feels
/// through the steering wheel.
#[derive(Debug, Clone)]
pub struct SteeringFeedback {
    /// Self-aligning torque gain.
    pub aligning_gain: Real,
    /// Damping coefficient for steering rate.
    pub damping: Real,
    /// Return-to-centre spring rate (N*m/rad).
    pub centering_spring: Real,
    /// Maximum feedback force (N).
    pub max_force: Real,
}
impl SteeringFeedback {
    /// Create a new steering feedback model.
    pub fn new(
        aligning_gain: Real,
        damping: Real,
        centering_spring: Real,
        max_force: Real,
    ) -> Self {
        Self {
            aligning_gain,
            damping,
            centering_spring,
            max_force,
        }
    }
    /// Compute the feedback force felt by the driver.
    ///
    /// * `lateral_force` — tire lateral force (N).
    /// * `steer_angle` — current road wheel angle (rad).
    /// * `steer_rate` — rate of change of steering angle (rad/s).
    pub fn compute_force(&self, lateral_force: Real, steer_angle: Real, steer_rate: Real) -> Real {
        let aligning = -lateral_force * self.aligning_gain;
        let centering = -steer_angle * self.centering_spring;
        let damping_force = -steer_rate * self.damping;
        let total = aligning + centering + damping_force;
        total.clamp(-self.max_force, self.max_force)
    }
}
/// Four-wheel steering model with speed-dependent rear steering.
///
/// At low speed the rear wheels steer out-of-phase (counter-steer) for
/// tighter turning; at high speed they steer in-phase for stability.
#[derive(Debug, Clone)]
pub struct FourWheelSteering {
    /// Scaling factor applied to the front angle for the rear (can be negative).
    pub front_rate: Real,
    /// Speed threshold above which rear steers in-phase (m/s).
    pub rear_rate: Real,
}
impl FourWheelSteering {
    /// Create a new four-wheel steering model.
    pub fn new(front_rate: Real, rear_rate: Real) -> Self {
        Self {
            front_rate,
            rear_rate,
        }
    }
    /// Compute the rear wheel steering angle.
    ///
    /// At low speed (`speed` < `rear_rate`), rear steers counter to front.
    /// Above `rear_rate`, rear steers in phase.
    ///
    /// Returns the rear steering angle in radians.
    pub fn rear_angle(&self, front_angle: Real, speed: Real) -> Real {
        let t = (speed / self.rear_rate.max(1e-6)).clamp(0.0, 2.0);
        let ratio = self.front_rate * (t - 1.0);
        front_angle * ratio
    }
}
/// Understeer / oversteer gradient analyser.
///
/// The understeer gradient K (rad/g) describes how the required steering
/// angle changes with lateral acceleration:
///
///   δ = L/R + K · Ay
///
/// where δ is front steer angle, L is wheelbase, R is turn radius, Ay is
/// lateral acceleration.  K > 0 → understeer; K < 0 → oversteer.
#[derive(Debug, Clone)]
pub struct UndersteerAnalyzer {
    /// Wheelbase (meters).
    pub wheelbase: Real,
    /// Front axle cornering stiffness (N/rad).
    pub cf: Real,
    /// Rear axle cornering stiffness (N/rad).
    pub cr: Real,
    /// Vehicle mass (kg).
    pub mass: Real,
    /// Front axle load fraction (0..1).
    pub front_load_fraction: Real,
}
impl UndersteerAnalyzer {
    /// Create a new analyzer.
    pub fn new(wheelbase: Real, cf: Real, cr: Real, mass: Real, front_load_fraction: Real) -> Self {
        Self {
            wheelbase,
            cf,
            cr,
            mass,
            front_load_fraction: front_load_fraction.clamp(0.0, 1.0),
        }
    }
    /// Compute the understeer gradient K (rad/g = rad·s²/m).
    ///
    /// Positive → understeer, negative → oversteer, zero → neutral steer.
    ///
    /// Formula (linear bicycle model):
    ///   K = (m / L²) * (a/Cf - b/Cr)
    /// where a = front axle to CoM, b = rear axle to CoM.
    pub fn understeer_gradient(&self) -> Real {
        let g = 9.81_f64;
        let a = self.wheelbase * self.front_load_fraction;
        let b = self.wheelbase * (1.0 - self.front_load_fraction);
        (self.mass / (self.wheelbase * self.wheelbase))
            * (a / self.cf.max(1.0) - b / self.cr.max(1.0))
            * g
    }
    /// Predicted steering angle for a given speed and turn radius.
    ///
    /// δ = L/R + K · (v²/R) / g
    pub fn predicted_steer_angle(&self, speed: Real, turn_radius: Real) -> Real {
        let g = 9.81_f64;
        if turn_radius.abs() < 1e-6 {
            return 0.0;
        }
        let ackermann_term = self.wheelbase / turn_radius;
        let ay = speed * speed / turn_radius;
        let k = self.understeer_gradient();
        ackermann_term + k * ay / g
    }
    /// Critical speed above which an oversteering vehicle becomes unstable.
    ///
    /// Returns `None` for understeer or neutral steer vehicles (K >= 0).
    pub fn critical_speed(&self) -> Option<Real> {
        let k = self.understeer_gradient();
        if k >= 0.0 {
            return None;
        }
        let g = 9.81_f64;
        let v2 = -g * self.wheelbase / k;
        if v2 > 0.0 { Some(v2.sqrt()) } else { None }
    }
    /// Characteristic speed for an understeering vehicle.
    ///
    /// The speed at which the steering response is 50 % of its low-speed value.
    /// Returns `None` for oversteering vehicles (K < 0).
    pub fn characteristic_speed(&self) -> Option<Real> {
        let k = self.understeer_gradient();
        if k <= 0.0 {
            return None;
        }
        let g = 9.81_f64;
        let v2 = g * self.wheelbase / k;
        if v2 > 0.0 { Some(v2.sqrt()) } else { None }
    }
    /// Classify the vehicle handling: "understeer", "oversteer", or "neutral".
    pub fn classify(&self) -> &'static str {
        let k = self.understeer_gradient();
        if k > 1e-6 {
            "understeer"
        } else if k < -1e-6 {
            "oversteer"
        } else {
            "neutral"
        }
    }
}
/// Simplified model-predictive steering controller for path tracking.
///
/// Uses a finite horizon bicycle model to compute the steering angle that
/// minimises future cross-track error and heading error.
#[derive(Debug, Clone)]
pub struct ModelPredictiveSteering {
    /// Prediction horizon length (seconds).
    pub horizon: Real,
    /// Number of prediction steps.
    pub steps: usize,
    /// Wheelbase (m).
    pub wheelbase: Real,
    /// Maximum steering angle (rad).
    pub max_steer: Real,
    /// Cross-track error weight.
    pub w_cte: Real,
    /// Heading error weight.
    pub w_heading: Real,
    /// Steering effort weight.
    pub w_steer: Real,
}
impl ModelPredictiveSteering {
    /// Create a new MPC steering controller.
    pub fn new(
        horizon: Real,
        steps: usize,
        wheelbase: Real,
        max_steer: Real,
        w_cte: Real,
        w_heading: Real,
        w_steer: Real,
    ) -> Self {
        Self {
            horizon,
            steps,
            wheelbase,
            max_steer,
            w_cte,
            w_heading,
            w_steer,
        }
    }
    /// Simulate one step of the bicycle model.
    ///
    /// * `x, y, psi` — current state (m, m, rad).
    /// * `v`         — speed (m/s).
    /// * `steer`     — steering angle (rad).
    /// * `dt`        — time step (s).
    ///
    /// Returns `(x_new, y_new, psi_new)`.
    pub fn bicycle_step(
        &self,
        x: Real,
        y: Real,
        psi: Real,
        v: Real,
        steer: Real,
        dt: Real,
    ) -> (Real, Real, Real) {
        let x_new = x + v * psi.cos() * dt;
        let y_new = y + v * psi.sin() * dt;
        let psi_new = psi + v * steer.tan() / self.wheelbase.max(1e-9) * dt;
        (x_new, y_new, psi_new)
    }
    /// Compute cost for a given steering angle over the prediction horizon.
    ///
    /// `target_y` is the desired lateral position (0 = on track).
    pub fn compute_cost(
        &self,
        x0: Real,
        y0: Real,
        psi0: Real,
        v: Real,
        steer: Real,
        target_y: Real,
    ) -> Real {
        let dt = self.horizon / self.steps.max(1) as Real;
        let mut x = x0;
        let mut y = y0;
        let mut psi = psi0;
        let mut cost = 0.0;
        for _ in 0..self.steps {
            let (nx, ny, np) = self.bicycle_step(x, y, psi, v, steer, dt);
            let cte = ny - target_y;
            cost +=
                self.w_cte * cte * cte + self.w_heading * psi * psi + self.w_steer * steer * steer;
            x = nx;
            y = ny;
            psi = np;
        }
        cost
    }
    /// Grid search for the optimal steering angle (simplified MPC).
    ///
    /// Samples `n_samples` candidate angles and returns the one with minimum cost.
    pub fn optimal_steer(
        &self,
        x0: Real,
        y0: Real,
        psi0: Real,
        v: Real,
        target_y: Real,
        n_samples: usize,
    ) -> Real {
        let n = n_samples.max(3);
        let mut best_steer = 0.0;
        let mut best_cost = f64::INFINITY;
        for i in 0..n {
            let t = i as Real / (n - 1) as Real;
            let steer = -self.max_steer + t * 2.0 * self.max_steer;
            let cost = self.compute_cost(x0, y0, psi0, v, steer, target_y);
            if cost < best_cost {
                best_cost = cost;
                best_steer = steer;
            }
        }
        best_steer.clamp(-self.max_steer, self.max_steer)
    }
}
