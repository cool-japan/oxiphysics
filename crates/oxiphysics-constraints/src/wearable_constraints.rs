// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0
//! Wearable robotics constraints: exoskeleton, prosthetics, assistive devices.
//!
//! Implements joint-level and system-level models for powered exoskeletons,
//! prosthetic limbs, rehabilitation robots, and fall-prevention controllers.

// ── Series Elastic Actuator / Exoskeleton Joint ───────────────────────────────

/// Powered exoskeleton joint with a compliant series elastic actuator (SEA).
///
/// The SEA torque model: `τ = k_s * (θ_motor - θ_output)`.
/// An impedance controller tracks a reference trajectory:
/// `τ = k * (θ_ref - θ) + b * (dθ_ref - dθ)`.
#[derive(Debug, Clone)]
pub struct ExoskeletonJoint {
    /// Spring stiffness of the series elastic element \[N·m/rad\].
    pub spring_stiffness: f64,
    /// Motor angle \[rad\].
    pub motor_angle: f64,
    /// Output (joint) angle \[rad\].
    pub output_angle: f64,
    /// Reference angle for impedance control \[rad\].
    pub ref_angle: f64,
    /// Reference angular velocity for impedance control \[rad/s\].
    pub ref_velocity: f64,
    /// Impedance stiffness \[N·m/rad\].
    pub impedance_k: f64,
    /// Impedance damping \[N·m·s/rad\].
    pub impedance_b: f64,
    /// Minimum output angle limit \[rad\].
    pub angle_min: f64,
    /// Maximum output angle limit \[rad\].
    pub angle_max: f64,
    /// Maximum actuator torque \[N·m\].
    pub torque_limit: f64,
}

impl ExoskeletonJoint {
    /// Create a new `ExoskeletonJoint` with given SEA stiffness and impedance parameters.
    pub fn new(
        spring_stiffness: f64,
        impedance_k: f64,
        impedance_b: f64,
        angle_min: f64,
        angle_max: f64,
        torque_limit: f64,
    ) -> Self {
        Self {
            spring_stiffness,
            motor_angle: 0.0,
            output_angle: 0.0,
            ref_angle: 0.0,
            ref_velocity: 0.0,
            impedance_k,
            impedance_b,
            angle_min,
            angle_max,
            torque_limit,
        }
    }

    /// Compute the SEA torque from the spring deflection.
    ///
    /// `τ_sea = k_s * (θ_motor - θ_output)`
    pub fn sea_torque(&self) -> f64 {
        self.spring_stiffness * (self.motor_angle - self.output_angle)
    }

    /// Compute the impedance control torque.
    ///
    /// `τ_imp = k * (θ_ref - θ) + b * (dθ_ref - dθ)`
    pub fn impedance_torque(&self, current_velocity: f64) -> f64 {
        let pos_error = self.ref_angle - self.output_angle;
        let vel_error = self.ref_velocity - current_velocity;
        let tau = self.impedance_k * pos_error + self.impedance_b * vel_error;
        tau.max(-self.torque_limit).min(self.torque_limit)
    }

    /// Check whether the joint is within its range of motion.
    pub fn within_rom(&self) -> bool {
        self.output_angle >= self.angle_min && self.output_angle <= self.angle_max
    }
}

// ── Assistance Torque / Gait-Phase Profile ────────────────────────────────────

/// Gait phase \[0, 1\] and stance/swing detection from foot pressure.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum GaitPhaseState {
    /// Foot is in contact with the ground (stance phase).
    Stance,
    /// Foot is off the ground (swing phase).
    Swing,
}

/// Gait-phase-dependent assistive torque profile using a cubic spline approximation.
///
/// Torque is defined by four control points over the normalized gait cycle `[0, 1]`.
/// Between points the torque is interpolated with a cubic Hermite curve.
#[derive(Debug, Clone)]
pub struct AssistanceTorque {
    /// Foot pressure threshold for stance detection \[N\].
    pub pressure_threshold: f64,
    /// Control-point phases (4 values in \[0, 1\]).
    pub phase_points: [f64; 4],
    /// Torque values at the control points \[N·m\].
    pub torque_points: [f64; 4],
    /// Current gait phase \[0, 1\].
    pub current_phase: f64,
}

impl AssistanceTorque {
    /// Create a new `AssistanceTorque` profile.
    pub fn new(pressure_threshold: f64, phase_points: [f64; 4], torque_points: [f64; 4]) -> Self {
        Self {
            pressure_threshold,
            phase_points,
            torque_points,
            current_phase: 0.0,
        }
    }

    /// Detect stance or swing from foot pressure.
    pub fn detect_phase(&self, foot_pressure: f64) -> GaitPhaseState {
        if foot_pressure >= self.pressure_threshold {
            GaitPhaseState::Stance
        } else {
            GaitPhaseState::Swing
        }
    }

    /// Evaluate the torque profile at a given gait phase via piecewise cubic Hermite.
    ///
    /// Returns 0 for phases outside the defined range.
    pub fn torque_at_phase(&self, phase: f64) -> f64 {
        let p = &self.phase_points;
        let t = &self.torque_points;
        if phase <= p[0] {
            return t[0];
        }
        if phase >= p[3] {
            return t[3];
        }
        // Find segment
        for i in 0..3 {
            if phase <= p[i + 1] {
                let u = (phase - p[i]) / (p[i + 1] - p[i]);
                return cubic_hermite(t[i], t[i + 1], u);
            }
        }
        t[3]
    }
}

/// Cubic Hermite interpolation between `a` and `b` at parameter `u` in \[0, 1\].
fn cubic_hermite(a: f64, b: f64, u: f64) -> f64 {
    let u2 = u * u;
    let u3 = u2 * u;
    a * (2.0 * u3 - 3.0 * u2 + 1.0) + b * (-2.0 * u3 + 3.0 * u2)
}

// ── Gait Analysis ─────────────────────────────────────────────────────────────

/// Gait analysis: step detection, phase estimation, cadence, walking speed.
#[derive(Debug, Clone)]
pub struct GaitAnalysis {
    /// Stride time (time for a full cycle) \[s\].
    pub stride_time: f64,
    /// Step length \[m\].
    pub step_length: f64,
    /// Time elapsed in current stride \[s\].
    pub elapsed: f64,
    /// Cadence \[steps/min\].
    pub cadence: f64,
    /// Estimated walking speed \[m/s\].
    pub walking_speed: f64,
    /// Number of steps counted.
    pub step_count: u32,
}

impl GaitAnalysis {
    /// Create a new `GaitAnalysis` with given stride time and step length.
    pub fn new(stride_time: f64, step_length: f64) -> Self {
        let cadence = if stride_time > 1e-9 {
            2.0 * 60.0 / stride_time
        } else {
            0.0
        };
        let walking_speed = if stride_time > 1e-9 {
            2.0 * step_length / stride_time
        } else {
            0.0
        };
        Self {
            stride_time,
            step_length,
            elapsed: 0.0,
            cadence,
            walking_speed,
            step_count: 0,
        }
    }

    /// Advance time by `dt` and return the normalized gait phase \[0, 1\].
    pub fn advance(&mut self, dt: f64) -> f64 {
        if self.stride_time < 1e-9 {
            return 0.0;
        }
        self.elapsed += dt;
        if self.elapsed >= self.stride_time {
            self.elapsed -= self.stride_time;
            self.step_count += 2; // two steps per stride
        }
        self.elapsed / self.stride_time
    }

    /// Current gait phase without advancing time.
    pub fn phase(&self) -> f64 {
        if self.stride_time < 1e-9 {
            return 0.0;
        }
        (self.elapsed / self.stride_time).min(1.0)
    }
}

// ── Human Motion Model / Musculoskeletal ─────────────────────────────────────

/// Simplified musculoskeletal joint constraint: ROM limits, torque limits, metabolic cost.
#[derive(Debug, Clone)]
pub struct HumanMotionModel {
    /// Joint name/label.
    pub joint_name: String,
    /// Minimum angle (ROM lower bound) \[rad\].
    pub rom_min: f64,
    /// Maximum angle (ROM upper bound) \[rad\].
    pub rom_max: f64,
    /// Maximum voluntary torque \[N·m\].
    pub max_torque: f64,
    /// Metabolic cost coefficient \[W·s²/m²\] for speed-squared model.
    pub metabolic_coefficient: f64,
}

impl HumanMotionModel {
    /// Create a new `HumanMotionModel`.
    pub fn new(
        joint_name: &str,
        rom_min: f64,
        rom_max: f64,
        max_torque: f64,
        metabolic_coefficient: f64,
    ) -> Self {
        Self {
            joint_name: joint_name.to_string(),
            rom_min,
            rom_max,
            max_torque,
            metabolic_coefficient,
        }
    }

    /// Clamp an angle to the ROM.
    pub fn clamp_angle(&self, angle: f64) -> f64 {
        angle.max(self.rom_min).min(self.rom_max)
    }

    /// Clamp a torque to the maximum voluntary torque.
    pub fn clamp_torque(&self, torque: f64) -> f64 {
        torque.max(-self.max_torque).min(self.max_torque)
    }

    /// Estimate metabolic cost \[W\] given walking speed \[m/s\].
    ///
    /// Simple quadratic model: `P = c * v²`.
    pub fn metabolic_cost(&self, speed: f64) -> f64 {
        self.metabolic_coefficient * speed * speed
    }
}

// ── Cable Actuator / Bowden Cable ─────────────────────────────────────────────

/// Bowden cable transmission with capstan friction model.
///
/// The capstan equation: `T_out = T_in * exp(μ * θ_wrap)`.
/// Backlash introduces a dead-zone before force transmission.
#[derive(Debug, Clone)]
pub struct CableActuator {
    /// Cable stiffness \[N/m\].
    pub stiffness: f64,
    /// Coefficient of friction (Coulomb, dimensionless).
    pub friction_coeff: f64,
    /// Wrap angle of cable around pulley/sheath \[rad\].
    pub wrap_angle: f64,
    /// Backlash (dead-zone in cable elongation) \[m\].
    pub backlash: f64,
    /// Current cable elongation \[m\].
    pub elongation: f64,
}

impl CableActuator {
    /// Create a new `CableActuator`.
    pub fn new(stiffness: f64, friction_coeff: f64, wrap_angle: f64, backlash: f64) -> Self {
        Self {
            stiffness,
            friction_coeff,
            wrap_angle,
            backlash,
            elongation: 0.0,
        }
    }

    /// Compute output tension from input tension using the capstan equation.
    ///
    /// `T_out = T_in * e^{μ * θ}`
    pub fn capstan_tension(&self, t_in: f64) -> f64 {
        t_in * (self.friction_coeff * self.wrap_angle).exp()
    }

    /// Compute cable force accounting for backlash and stiffness.
    ///
    /// Force is zero within the backlash dead-zone.
    pub fn cable_force(&self) -> f64 {
        let effective = (self.elongation.abs() - self.backlash).max(0.0);
        self.stiffness * effective * self.elongation.signum()
    }
}

// ── Prosthetic Hand ──────────────────────────────────────────────────────────

/// Under-actuated prosthetic hand: 1 motor drives N compliant finger joints.
///
/// Adaptive grasping through joint compliance allows passive shape-matching.
#[derive(Debug, Clone)]
pub struct ProstheticHand {
    /// Number of finger joints.
    pub num_joints: usize,
    /// Finger joint angles \[rad\].
    pub joint_angles: Vec<f64>,
    /// Compliance of each joint spring \[N·m/rad\].
    pub joint_stiffness: Vec<f64>,
    /// Maximum grip force \[N\].
    pub max_grip_force: f64,
    /// Current motor position \[rad\].
    pub motor_position: f64,
}

impl ProstheticHand {
    /// Create a new `ProstheticHand` with `num_joints` joints and uniform stiffness.
    pub fn new(num_joints: usize, joint_stiffness: f64, max_grip_force: f64) -> Self {
        Self {
            num_joints,
            joint_angles: vec![0.0; num_joints],
            joint_stiffness: vec![joint_stiffness; num_joints],
            max_grip_force,
            motor_position: 0.0,
        }
    }

    /// Compute grip force for joint `i` from its spring deflection.
    ///
    /// `F_grip = k_i * (θ_motor - θ_joint_i)`.
    pub fn grip_force(&self, joint_idx: usize) -> f64 {
        if joint_idx >= self.num_joints {
            return 0.0;
        }
        let deflection = self.motor_position - self.joint_angles[joint_idx];
        let force = self.joint_stiffness[joint_idx] * deflection;
        force.max(-self.max_grip_force).min(self.max_grip_force)
    }

    /// Estimate total grip force as sum across all joints.
    pub fn total_grip_force(&self) -> f64 {
        (0..self.num_joints)
            .map(|i| self.grip_force(i).abs())
            .sum::<f64>()
            .min(self.max_grip_force)
    }

    /// Move the motor to a new position (open/close hand).
    pub fn set_motor_position(&mut self, pos: f64) {
        self.motor_position = pos;
    }
}

// ── Fall Detection ────────────────────────────────────────────────────────────

/// Fall detection via inertial sensor thresholding.
///
/// Triggers if angular velocity magnitude exceeds `omega_threshold` or
/// linear acceleration magnitude exceeds `accel_threshold`.
#[derive(Debug, Clone)]
pub struct FallDetection {
    /// Angular velocity threshold \[rad/s\].
    pub omega_threshold: f64,
    /// Acceleration threshold \[m/s²\].
    pub accel_threshold: f64,
    /// Number of consecutive frames above threshold needed to trigger.
    pub trigger_frames: u32,
    /// Current count of consecutive frames above threshold.
    pub frame_count: u32,
    /// Whether a fall has been detected.
    pub fall_detected: bool,
    /// Stiffening gain applied to exoskeleton joints when fall detected.
    pub stiffening_gain: f64,
}

impl FallDetection {
    /// Create a new `FallDetection` detector.
    pub fn new(
        omega_threshold: f64,
        accel_threshold: f64,
        trigger_frames: u32,
        stiffening_gain: f64,
    ) -> Self {
        Self {
            omega_threshold,
            accel_threshold,
            trigger_frames,
            frame_count: 0,
            fall_detected: false,
            stiffening_gain,
        }
    }

    /// Update detector with new inertial measurements.
    ///
    /// `omega` is angular velocity magnitude \[rad/s\], `accel` is linear
    /// acceleration magnitude \[m/s²\]. Returns `true` if a fall is detected.
    pub fn update(&mut self, omega: f64, accel: f64) -> bool {
        if omega > self.omega_threshold || accel > self.accel_threshold {
            self.frame_count += 1;
        } else {
            self.frame_count = 0;
        }
        if self.frame_count >= self.trigger_frames {
            self.fall_detected = true;
        }
        self.fall_detected
    }

    /// Reset the fall detection state.
    pub fn reset(&mut self) {
        self.frame_count = 0;
        self.fall_detected = false;
    }

    /// Compute effective joint stiffness during fall prevention.
    ///
    /// Multiplies the nominal stiffness by the stiffening gain when a fall is detected.
    pub fn effective_stiffness(&self, nominal_stiffness: f64) -> f64 {
        if self.fall_detected {
            nominal_stiffness * self.stiffening_gain
        } else {
            nominal_stiffness
        }
    }
}

// ── Rehabilitation Controller ─────────────────────────────────────────────────

/// Assist-as-needed rehabilitation controller with error deadband.
///
/// The controller applies assistive torque only when the patient's tracking
/// error exceeds a deadband (the "patient effort zone"). This encourages
/// active participation.
#[derive(Debug, Clone)]
pub struct RehabilitationController {
    /// Error deadband width (half-width, symmetric) \[rad\].
    pub deadband: f64,
    /// Feedback gain \[N·m/rad\].
    pub gain: f64,
    /// Maximum assistive torque \[N·m\].
    pub max_torque: f64,
    /// Integral gain \[N·m/(rad·s)\].
    pub integral_gain: f64,
    /// Current integral accumulator \[N·m·s\].
    pub integral: f64,
    /// Integration anti-windup limit \[N·m·s\].
    pub integral_limit: f64,
}

impl RehabilitationController {
    /// Create a new `RehabilitationController`.
    pub fn new(
        deadband: f64,
        gain: f64,
        max_torque: f64,
        integral_gain: f64,
        integral_limit: f64,
    ) -> Self {
        Self {
            deadband,
            gain,
            max_torque,
            integral_gain,
            integral: 0.0,
            integral_limit,
        }
    }

    /// Compute the assist torque given tracking error and time step.
    ///
    /// No torque is applied within the deadband. Outside the deadband the
    /// torque is `gain * (|error| - deadband) * sign(error)` plus an
    /// integral term.
    pub fn compute_torque(&mut self, error: f64, dt: f64) -> f64 {
        let excess = error.abs() - self.deadband;
        if excess <= 0.0 {
            // Within deadband — let the patient work
            self.integral = 0.0;
            return 0.0;
        }
        let signed_excess = excess * error.signum();
        self.integral = (self.integral + signed_excess * dt)
            .max(-self.integral_limit)
            .min(self.integral_limit);
        let tau = self.gain * signed_excess + self.integral_gain * self.integral;
        tau.max(-self.max_torque).min(self.max_torque)
    }

    /// Reset the integral state.
    pub fn reset_integral(&mut self) {
        self.integral = 0.0;
    }
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // --- ExoskeletonJoint / SEA ---

    #[test]
    fn test_sea_torque_proportional_to_deflection() {
        let mut j = ExoskeletonJoint::new(100.0, 50.0, 5.0, -1.5, 1.5, 200.0);
        j.motor_angle = 0.1;
        j.output_angle = 0.0;
        let tau = j.sea_torque();
        assert!((tau - 10.0).abs() < 1e-10, "τ={tau}");
    }

    #[test]
    fn test_sea_torque_zero_deflection() {
        let j = ExoskeletonJoint::new(100.0, 50.0, 5.0, -1.5, 1.5, 200.0);
        assert!(j.sea_torque().abs() < 1e-10);
    }

    #[test]
    fn test_sea_torque_negative_deflection() {
        let mut j = ExoskeletonJoint::new(200.0, 50.0, 5.0, -1.5, 1.5, 300.0);
        j.motor_angle = -0.05;
        let tau = j.sea_torque();
        assert!((tau + 10.0).abs() < 1e-10, "τ={tau}");
    }

    #[test]
    fn test_impedance_torque_pure_position_error() {
        let mut j = ExoskeletonJoint::new(100.0, 50.0, 5.0, -1.5, 1.5, 200.0);
        j.ref_angle = 0.2;
        j.output_angle = 0.0;
        // no velocity error
        let tau = j.impedance_torque(0.0);
        // k * 0.2 = 50 * 0.2 = 10 N·m
        assert!((tau - 10.0).abs() < 1e-10, "τ={tau}");
    }

    #[test]
    fn test_impedance_torque_damping() {
        let mut j = ExoskeletonJoint::new(100.0, 0.0, 10.0, -1.5, 1.5, 200.0);
        j.ref_velocity = 1.0;
        // current velocity 0 → velocity error = 1
        let tau = j.impedance_torque(0.0);
        assert!((tau - 10.0).abs() < 1e-10, "τ={tau}");
    }

    #[test]
    fn test_impedance_torque_clamped() {
        let mut j = ExoskeletonJoint::new(100.0, 1000.0, 0.0, -1.5, 1.5, 5.0);
        j.ref_angle = 1.0;
        let tau = j.impedance_torque(0.0);
        assert!((tau - 5.0).abs() < 1e-10, "τ clamped={tau}");
    }

    #[test]
    fn test_joint_within_rom() {
        let mut j = ExoskeletonJoint::new(100.0, 50.0, 5.0, -1.5, 1.5, 200.0);
        j.output_angle = 1.0;
        assert!(j.within_rom());
        j.output_angle = 2.0;
        assert!(!j.within_rom());
    }

    // --- AssistanceTorque ---

    #[test]
    fn test_detect_stance() {
        let at = AssistanceTorque::new(20.0, [0.0, 0.3, 0.6, 1.0], [0.0, 10.0, 5.0, 0.0]);
        assert_eq!(at.detect_phase(50.0), GaitPhaseState::Stance);
        assert_eq!(at.detect_phase(5.0), GaitPhaseState::Swing);
    }

    #[test]
    fn test_torque_at_endpoints() {
        let at = AssistanceTorque::new(20.0, [0.0, 0.3, 0.6, 1.0], [1.0, 10.0, 5.0, 2.0]);
        let t0 = at.torque_at_phase(0.0);
        let t1 = at.torque_at_phase(1.0);
        assert!((t0 - 1.0).abs() < 1e-9, "t0={t0}");
        assert!((t1 - 2.0).abs() < 1e-9, "t1={t1}");
    }

    #[test]
    fn test_torque_monotone_between_rising_points() {
        let at = AssistanceTorque::new(20.0, [0.0, 0.5, 0.75, 1.0], [0.0, 20.0, 20.0, 0.0]);
        // Between 0 and 0.5 torque should be non-decreasing
        let mut prev = at.torque_at_phase(0.0);
        for i in 1..=10 {
            let phase = i as f64 * 0.05;
            let curr = at.torque_at_phase(phase);
            assert!(
                curr >= prev - 1e-9,
                "phase={phase:.2} curr={curr} prev={prev}"
            );
            prev = curr;
        }
    }

    // --- GaitAnalysis ---

    #[test]
    fn test_gait_phase_monotonically_increasing() {
        let mut ga = GaitAnalysis::new(1.0, 0.75);
        let mut prev = ga.advance(0.0);
        for _ in 0..9 {
            let phase = ga.advance(0.05);
            assert!(
                phase >= prev - 1e-9,
                "non-monotone: phase={phase} prev={prev}"
            );
            prev = phase;
        }
    }

    #[test]
    fn test_gait_phase_wraps_to_zero() {
        let mut ga = GaitAnalysis::new(1.0, 0.75);
        // advance past one full stride
        for _ in 0..25 {
            ga.advance(0.05);
        }
        // after 1.25 s of a 1.0 s stride, elapsed should be 0.25
        let phase = ga.phase();
        assert!(phase < 0.5, "phase after wrap={phase}");
    }

    #[test]
    fn test_cadence_calculation() {
        let ga = GaitAnalysis::new(1.0, 0.75);
        // 2 steps per stride, 60 s/min → 120 steps/min
        assert!((ga.cadence - 120.0).abs() < 1e-9, "cadence={}", ga.cadence);
    }

    #[test]
    fn test_walking_speed_calculation() {
        let ga = GaitAnalysis::new(1.0, 0.75);
        // speed = 2 * 0.75 / 1.0 = 1.5 m/s
        assert!((ga.walking_speed - 1.5).abs() < 1e-9);
    }

    // --- HumanMotionModel ---

    #[test]
    fn test_metabolic_cost_increases_with_speed() {
        let model = HumanMotionModel::new("hip", -0.5, 1.5, 80.0, 3.0);
        let p1 = model.metabolic_cost(1.0);
        let p2 = model.metabolic_cost(2.0);
        assert!(p2 > p1, "p2={p2} p1={p1}");
    }

    #[test]
    fn test_metabolic_cost_zero_speed() {
        let model = HumanMotionModel::new("knee", -0.2, 1.8, 60.0, 3.0);
        assert!(model.metabolic_cost(0.0).abs() < 1e-10);
    }

    #[test]
    fn test_clamp_angle_within_rom() {
        let model = HumanMotionModel::new("ankle", -0.3, 0.5, 40.0, 2.5);
        assert!((model.clamp_angle(0.2) - 0.2).abs() < 1e-10);
        assert!((model.clamp_angle(1.0) - 0.5).abs() < 1e-10);
        assert!((model.clamp_angle(-1.0) + 0.3).abs() < 1e-10);
    }

    #[test]
    fn test_clamp_torque_within_limit() {
        let model = HumanMotionModel::new("elbow", -1.5, 1.5, 30.0, 2.0);
        assert!((model.clamp_torque(20.0) - 20.0).abs() < 1e-10);
        assert!((model.clamp_torque(50.0) - 30.0).abs() < 1e-10);
        assert!((model.clamp_torque(-50.0) + 30.0).abs() < 1e-10);
    }

    // --- CableActuator (Capstan) ---

    #[test]
    fn test_capstan_friction_amplifies_force() {
        let cable = CableActuator::new(500.0, 0.2, std::f64::consts::PI, 0.005);
        let t_in = 10.0;
        let t_out = cable.capstan_tension(t_in);
        assert!(t_out > t_in, "t_out={t_out} should > t_in={t_in}");
    }

    #[test]
    fn test_capstan_zero_friction() {
        let cable = CableActuator::new(500.0, 0.0, std::f64::consts::PI, 0.0);
        let t_out = cable.capstan_tension(10.0);
        assert!((t_out - 10.0).abs() < 1e-9, "t_out={t_out}");
    }

    #[test]
    fn test_cable_force_zero_elongation() {
        let cable = CableActuator::new(500.0, 0.2, 1.0, 0.002);
        assert!(cable.cable_force().abs() < 1e-10);
    }

    #[test]
    fn test_cable_force_within_backlash() {
        let mut cable = CableActuator::new(500.0, 0.2, 1.0, 0.01);
        cable.elongation = 0.005; // within backlash
        assert!(cable.cable_force().abs() < 1e-10);
    }

    #[test]
    fn test_cable_force_beyond_backlash() {
        let mut cable = CableActuator::new(1000.0, 0.2, 1.0, 0.01);
        cable.elongation = 0.02; // 0.01 beyond backlash
        let f = cable.cable_force();
        assert!((f - 10.0).abs() < 1e-9, "f={f}");
    }

    // --- ProstheticHand ---

    #[test]
    fn test_grip_force_within_limits() {
        let mut hand = ProstheticHand::new(3, 50.0, 100.0);
        hand.set_motor_position(0.5);
        for i in 0..3 {
            let f = hand.grip_force(i);
            assert!(f.abs() <= 100.0 + 1e-9, "grip_force[{i}]={f}");
        }
    }

    #[test]
    fn test_grip_force_proportional_to_deflection() {
        let mut hand = ProstheticHand::new(1, 50.0, 200.0);
        hand.set_motor_position(0.4);
        hand.joint_angles[0] = 0.2;
        let f = hand.grip_force(0);
        // deflection = 0.4 - 0.2 = 0.2, force = 50 * 0.2 = 10
        assert!((f - 10.0).abs() < 1e-9, "f={f}");
    }

    #[test]
    fn test_grip_force_invalid_index() {
        let hand = ProstheticHand::new(2, 50.0, 100.0);
        assert!(hand.grip_force(10).abs() < 1e-10);
    }

    #[test]
    fn test_total_grip_force_bounded() {
        let mut hand = ProstheticHand::new(5, 100.0, 50.0);
        hand.set_motor_position(1.0);
        let total = hand.total_grip_force();
        assert!(total <= 50.0 + 1e-9, "total={total}");
    }

    // --- FallDetection ---

    #[test]
    fn test_fall_detection_trigger_above_threshold() {
        let mut fd = FallDetection::new(2.0, 15.0, 2, 3.0);
        fd.update(3.0, 5.0);
        let detected = fd.update(3.0, 5.0);
        assert!(detected, "should detect fall");
    }

    #[test]
    fn test_fall_detection_no_trigger_below_threshold() {
        let mut fd = FallDetection::new(2.0, 15.0, 3, 3.0);
        let detected = fd.update(1.0, 5.0);
        assert!(!detected, "should not trigger");
    }

    #[test]
    fn test_fall_detection_reset() {
        let mut fd = FallDetection::new(1.0, 5.0, 1, 3.0);
        fd.update(2.0, 10.0);
        fd.reset();
        assert!(!fd.fall_detected);
        assert_eq!(fd.frame_count, 0);
    }

    #[test]
    fn test_stiffening_when_fall_detected() {
        let mut fd = FallDetection::new(1.0, 5.0, 1, 5.0);
        fd.update(2.0, 10.0);
        let stiff = fd.effective_stiffness(100.0);
        assert!((stiff - 500.0).abs() < 1e-9, "stiff={stiff}");
    }

    #[test]
    fn test_no_stiffening_without_fall() {
        let fd = FallDetection::new(1.0, 5.0, 1, 5.0);
        let stiff = fd.effective_stiffness(100.0);
        assert!((stiff - 100.0).abs() < 1e-9);
    }

    // --- RehabilitationController ---

    #[test]
    fn test_rehab_no_torque_within_deadband() {
        let mut ctrl = RehabilitationController::new(0.05, 100.0, 50.0, 10.0, 20.0);
        let tau = ctrl.compute_torque(0.02, 0.01);
        assert!(tau.abs() < 1e-10, "tau={tau}");
    }

    #[test]
    fn test_rehab_torque_outside_deadband() {
        let mut ctrl = RehabilitationController::new(0.0, 100.0, 200.0, 0.0, 0.0);
        let tau = ctrl.compute_torque(0.1, 0.01);
        assert!((tau - 10.0).abs() < 1e-9, "tau={tau}");
    }

    #[test]
    fn test_rehab_torque_clamped() {
        let mut ctrl = RehabilitationController::new(0.0, 10000.0, 50.0, 0.0, 0.0);
        let tau = ctrl.compute_torque(1.0, 0.01);
        assert!((tau - 50.0).abs() < 1e-9, "tau={tau}");
    }

    #[test]
    fn test_rehab_integral_resets_on_enter_deadband() {
        let mut ctrl = RehabilitationController::new(0.05, 100.0, 200.0, 50.0, 100.0);
        ctrl.compute_torque(0.2, 0.1); // outside deadband
        ctrl.compute_torque(0.01, 0.1); // inside deadband → resets integral
        assert!(ctrl.integral.abs() < 1e-10, "integral={}", ctrl.integral);
    }

    #[test]
    fn test_rehab_reset_integral() {
        let mut ctrl = RehabilitationController::new(0.0, 100.0, 200.0, 50.0, 100.0);
        ctrl.compute_torque(0.5, 0.1);
        ctrl.reset_integral();
        assert!(ctrl.integral.abs() < 1e-10);
    }

    // --- Stiffness / Damping trade-off (impedance) ---

    #[test]
    fn test_impedance_stiffness_dominates_at_low_velocity() {
        let mut j = ExoskeletonJoint::new(100.0, 200.0, 1.0, -2.0, 2.0, 1000.0);
        j.ref_angle = 0.5;
        j.output_angle = 0.0;
        j.ref_velocity = 0.0;
        // velocity error is small; position error large
        let tau_high_k = j.impedance_torque(0.0);
        assert!((tau_high_k - 100.0).abs() < 1e-9, "tau={tau_high_k}");
    }

    #[test]
    fn test_impedance_damping_dominates_at_high_velocity_error() {
        let mut j = ExoskeletonJoint::new(100.0, 0.0, 50.0, -2.0, 2.0, 1000.0);
        j.ref_velocity = 2.0;
        let tau = j.impedance_torque(0.0);
        assert!((tau - 100.0).abs() < 1e-9, "tau={tau}");
    }
}
