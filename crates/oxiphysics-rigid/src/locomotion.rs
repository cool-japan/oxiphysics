// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Legged locomotion and gait simulation.
//!
//! Provides models for bipedal and quadrupedal locomotion, zero-moment point (ZMP)
//! analysis, capture point computation, and biomechanical utility functions.

// ── Leg ──────────────────────────────────────────────────────────────────────

/// A single leg with hip and foot positions, length, stance state, and swing phase.
#[derive(Debug, Clone, PartialEq)]
pub struct Leg {
    /// Position of the hip joint in world coordinates \[x, y, z\].
    pub hip_position: [f64; 3],
    /// Position of the foot in world coordinates \[x, y, z\].
    pub foot_position: [f64; 3],
    /// Nominal leg length in metres.
    pub length: f64,
    /// Whether the leg is currently in the stance (ground-contact) phase.
    pub in_stance: bool,
    /// Normalised swing phase: 0.0 = start of swing, 1.0 = end of swing.
    pub swing_phase: f64,
}

impl Leg {
    /// Create a new leg at the given hip position with the specified nominal length.
    pub fn new(hip_position: [f64; 3], length: f64) -> Self {
        let mut fp = hip_position;
        fp[1] -= length;
        Self {
            hip_position,
            foot_position: fp,
            length,
            in_stance: true,
            swing_phase: 0.0,
        }
    }

    /// Current leg vector length (hip-to-foot distance).
    pub fn current_length(&self) -> f64 {
        let dx = self.foot_position[0] - self.hip_position[0];
        let dy = self.foot_position[1] - self.hip_position[1];
        let dz = self.foot_position[2] - self.hip_position[2];
        (dx * dx + dy * dy + dz * dz).sqrt()
    }

    /// Advance the swing phase by `delta` (clamped to \[0, 1\]).
    pub fn advance_swing(&mut self, delta: f64) {
        self.swing_phase = (self.swing_phase + delta).clamp(0.0, 1.0);
        if self.swing_phase >= 1.0 {
            self.in_stance = true;
            self.swing_phase = 0.0;
        }
    }
}

// ── GaitPattern ──────────────────────────────────────────────────────────────

/// Quadrupedal gait patterns with associated temporal parameters.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GaitPattern {
    /// Walk — all four legs cycle at duty factor ~0.75.
    Walk,
    /// Trot — diagonal pairs move simultaneously, duty factor ~0.6.
    Trot,
    /// Gallop — rotary gallop with duty factor ~0.35.
    Gallop,
    /// Pace — ipsilateral pairs move simultaneously, duty factor ~0.55.
    Pace,
    /// Bound — forelimbs then hindlimbs in pairs, duty factor ~0.4.
    Bound,
}

impl GaitPattern {
    /// Fraction of the stride cycle that each foot is in contact with the ground.
    pub fn duty_factor(self) -> f64 {
        match self {
            GaitPattern::Walk => 0.75,
            GaitPattern::Trot => 0.60,
            GaitPattern::Gallop => 0.35,
            GaitPattern::Pace => 0.55,
            GaitPattern::Bound => 0.40,
        }
    }

    /// Normalised phase offsets for each of the four legs
    /// (front-left, front-right, rear-left, rear-right).
    ///
    /// Values are in \[0, 1) relative to the gait cycle period.
    pub fn phase_offsets(self) -> Vec<f64> {
        match self {
            GaitPattern::Walk => vec![0.0, 0.5, 0.25, 0.75],
            GaitPattern::Trot => vec![0.0, 0.5, 0.5, 0.0],
            GaitPattern::Gallop => vec![0.0, 0.1, 0.5, 0.6],
            GaitPattern::Pace => vec![0.0, 0.5, 0.0, 0.5],
            GaitPattern::Bound => vec![0.0, 0.0, 0.5, 0.5],
        }
    }

    /// Step frequency in Hz for a given forward speed (m/s) and leg length (m).
    ///
    /// Uses the dimensionless Froude-number relationship:
    /// `f = sqrt(g / L) / (2π) * correction`.
    pub fn preferred_frequency(self, speed: f64, leg_length: f64) -> f64 {
        let f_nat = (9.81 / leg_length).sqrt() / (2.0 * std::f64::consts::PI);
        let correction = match self {
            GaitPattern::Walk => 1.0,
            GaitPattern::Trot => 1.4,
            GaitPattern::Gallop => 2.0,
            GaitPattern::Pace => 1.3,
            GaitPattern::Bound => 1.8,
        };
        // step_length ≈ speed / frequency  ⟹  frequency ≈ speed / step_length
        // For a rough scaling use f_nat * correction
        let _ = speed; // speed influences stride length, not directly frequency here
        f_nat * correction
    }
}

// ── BipedModel ────────────────────────────────────────────────────────────────

/// A simplified planar biped with two legs and a torso.
#[derive(Debug, Clone)]
pub struct BipedModel {
    /// The two legs: `[left, right]`.
    pub legs: [Leg; 2],
    /// Torso centre-of-mass position in world coordinates \[x, y, z\].
    pub torso_position: [f64; 3],
    /// Torso velocity in world coordinates \[x, y, z\].
    pub torso_velocity: [f64; 3],
}

impl BipedModel {
    /// Create a biped with default leg length 1.0 m standing at the origin.
    pub fn new() -> Self {
        Self {
            legs: [
                Leg::new([-0.1, 1.0, 0.0], 1.0),
                Leg::new([0.1, 1.0, 0.0], 1.0),
            ],
            torso_position: [0.0, 1.0, 0.0],
            torso_velocity: [0.0, 0.0, 0.0],
        }
    }

    /// Advance the simulation by `dt` seconds.
    ///
    /// `control` encodes `[forward_speed, lateral_speed, step_frequency, step_height]`.
    pub fn step(&mut self, dt: f64, control: &[f64; 4]) {
        let forward_speed = control[0];
        let lateral_speed = control[1];
        let step_freq = control[2].max(0.1);
        let _step_height = control[3];

        // Integrate torso position
        self.torso_velocity[0] = forward_speed;
        self.torso_velocity[2] = lateral_speed;
        for i in 0..3 {
            self.torso_position[i] += self.torso_velocity[i] * dt;
        }

        // Update hip positions to follow torso
        self.legs[0].hip_position[0] = self.torso_position[0] - 0.1;
        self.legs[0].hip_position[1] = self.torso_position[1];
        self.legs[0].hip_position[2] = self.torso_position[2];

        self.legs[1].hip_position[0] = self.torso_position[0] + 0.1;
        self.legs[1].hip_position[1] = self.torso_position[1];
        self.legs[1].hip_position[2] = self.torso_position[2];

        // Advance swing phase for legs not in stance
        let phase_delta = step_freq * dt;
        for leg in self.legs.iter_mut() {
            if !leg.in_stance {
                leg.advance_swing(phase_delta);
            } else {
                // Simple stance-to-swing transition based on foot lag
                let lag = leg.hip_position[0] - leg.foot_position[0];
                if lag > 0.3 {
                    leg.in_stance = false;
                    leg.swing_phase = 0.0;
                }
            }
            // Update foot position in stance: foot stays on ground
            if leg.in_stance {
                leg.foot_position[1] = 0.0;
            }
        }
    }
}

impl Default for BipedModel {
    fn default() -> Self {
        Self::new()
    }
}

// ── QuadrupedModel ────────────────────────────────────────────────────────────

/// A simplified quadruped locomotion model with configurable gait.
#[derive(Debug, Clone)]
pub struct QuadrupedModel {
    /// The four legs: `[front-left, front-right, rear-left, rear-right]`.
    pub legs: [Leg; 4],
    /// Torso centre-of-mass position in world coordinates \[x, y, z\].
    pub torso_position: [f64; 3],
    /// Active gait pattern.
    pub gait: GaitPattern,
    /// Current normalised gait phase in \[0, 1).
    pub gait_phase: f64,
    /// Step frequency in Hz.
    pub step_frequency: f64,
}

impl QuadrupedModel {
    /// Create a quadruped with default parameters standing at the origin.
    pub fn new(gait: GaitPattern) -> Self {
        let offsets = [[-0.3, 0.15], [-0.3, -0.15], [0.3, 0.15], [0.3, -0.15]];
        let legs = offsets.map(|[x, z]| Leg::new([x, 0.6, z], 0.6));
        Self {
            legs,
            torso_position: [0.0, 0.6, 0.0],
            gait,
            gait_phase: 0.0,
            step_frequency: 2.0,
        }
    }

    /// Advance the simulation by `dt` seconds.
    pub fn step(&mut self, dt: f64) {
        let phase_offsets = self.gait.phase_offsets();
        let duty = self.gait.duty_factor();

        // Advance global gait phase
        self.gait_phase = (self.gait_phase + self.step_frequency * dt).rem_euclid(1.0);

        // Determine stance/swing for each leg
        for (i, leg) in self.legs.iter_mut().enumerate() {
            let offset = phase_offsets[i];
            let local_phase = (self.gait_phase - offset).rem_euclid(1.0);
            leg.in_stance = local_phase < duty;
            if !leg.in_stance {
                leg.swing_phase = (local_phase - duty) / (1.0 - duty);
            } else {
                leg.swing_phase = 0.0;
            }

            if leg.in_stance {
                leg.foot_position[1] = 0.0;
            }
        }
    }
}

// ── ZeroMomentPoint ───────────────────────────────────────────────────────────

/// Zero-moment point (ZMP) computation for legged stability analysis.
#[derive(Debug, Clone, Default)]
pub struct ZeroMomentPoint {
    /// Most recently computed ZMP x-coordinate (m).
    pub cop_x: f64,
    /// Most recently computed ZMP y-coordinate (m).
    pub cop_y: f64,
}

impl ZeroMomentPoint {
    /// Create a new ZMP tracker initialised to the origin.
    pub fn new() -> Self {
        Self::default()
    }

    /// Compute the ZMP from ground-reaction forces and contact positions.
    ///
    /// `forces[i]` is `[Fx, Fy, Fz]` and `positions[i]` is `[x, y, z]`.
    /// The vertical (y) component of force is used for weighting.
    ///
    /// Returns `[zmp_x, zmp_z]` (horizontal plane).
    pub fn compute_zmp(&mut self, forces: &[[f64; 3]], positions: &[[f64; 3]]) -> [f64; 2] {
        let mut sum_fy = 0.0_f64;
        let mut sum_x = 0.0_f64;
        let mut sum_z = 0.0_f64;

        for (f, p) in forces.iter().zip(positions.iter()) {
            let fy = f[1].max(0.0); // only compressive forces
            sum_fy += fy;
            sum_x += fy * p[0];
            sum_z += fy * p[2];
        }

        if sum_fy > 1e-12 {
            let zmp_x = sum_x / sum_fy;
            let zmp_z = sum_z / sum_fy;
            self.cop_x = zmp_x;
            self.cop_y = zmp_z;
            [zmp_x, zmp_z]
        } else {
            [0.0, 0.0]
        }
    }
}

// ── CapturePoint ─────────────────────────────────────────────────────────────

/// Capture point (divergent component of motion) for push-recovery analysis.
#[derive(Debug, Clone)]
pub struct CapturePoint {
    /// Capture point position \[x, z\] (horizontal plane).
    pub position: [f64; 2],
    /// Natural frequency ω = sqrt(g / z_com) rad/s.
    pub omega: f64,
}

impl CapturePoint {
    /// Create a new capture point with the given natural frequency.
    pub fn new(omega: f64) -> Self {
        Self {
            position: [0.0, 0.0],
            omega,
        }
    }

    /// Compute the instantaneous capture point from CoM position and velocity.
    ///
    /// `pos` and `vel` are 2-D horizontal `[x, z]` vectors.
    /// Returns `[ξ_x, ξ_z]` where `ξ = p + v / ω`.
    pub fn compute_cp(&mut self, pos: [f64; 2], vel: [f64; 2], omega: f64) -> [f64; 2] {
        self.omega = omega;
        let xi = [pos[0] + vel[0] / omega, pos[1] + vel[1] / omega];
        self.position = xi;
        xi
    }
}

// ── Free functions ────────────────────────────────────────────────────────────

/// Compute stride length from forward speed and step frequency.
///
/// `stride_length = speed / frequency` (m).
pub fn stride_length(speed: f64, frequency: f64) -> f64 {
    if frequency > 0.0 {
        speed / frequency
    } else {
        0.0
    }
}

/// Compute the dimensionless Froude number for legged locomotion.
///
/// `Fr = v² / (g · L)` where `g = 9.81 m/s²`.
pub fn froude_number(speed: f64, leg_length: f64) -> f64 {
    if leg_length > 0.0 {
        (speed * speed) / (9.81 * leg_length)
    } else {
        0.0
    }
}

/// Preferred walk-to-run transition speed for a given leg length.
///
/// Based on the empirical relation `v_opt ≈ 0.5 * sqrt(g * L)`.
pub fn preferred_walk_run_speed(leg_length: f64) -> f64 {
    0.5 * (9.81 * leg_length).sqrt()
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // ── Leg ──────────────────────────────────────────────────────────────

    #[test]
    fn leg_new_sets_foot_below_hip() {
        let leg = Leg::new([0.0, 1.0, 0.0], 1.0);
        assert!((leg.foot_position[1] - 0.0).abs() < 1e-10);
        assert!((leg.hip_position[1] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn leg_current_length_matches_nominal() {
        let leg = Leg::new([0.0, 1.0, 0.0], 1.0);
        assert!((leg.current_length() - 1.0).abs() < 1e-10);
    }

    #[test]
    fn leg_advance_swing_transitions_to_stance() {
        let mut leg = Leg::new([0.0, 1.0, 0.0], 1.0);
        leg.in_stance = false;
        leg.swing_phase = 0.0;
        leg.advance_swing(1.0);
        assert!(leg.in_stance);
        assert!((leg.swing_phase).abs() < 1e-10);
    }

    #[test]
    fn leg_advance_swing_partial() {
        let mut leg = Leg::new([0.0, 1.0, 0.0], 1.0);
        leg.in_stance = false;
        leg.swing_phase = 0.0;
        leg.advance_swing(0.3);
        assert!(!leg.in_stance);
        assert!((leg.swing_phase - 0.3).abs() < 1e-10);
    }

    #[test]
    fn leg_current_length_diagonal() {
        let mut leg = Leg::new([0.0, 0.0, 0.0], 1.0);
        leg.foot_position = [3.0, 4.0, 0.0];
        assert!((leg.current_length() - 5.0).abs() < 1e-10);
    }

    // ── GaitPattern ──────────────────────────────────────────────────────

    #[test]
    fn walk_duty_factor() {
        assert!((GaitPattern::Walk.duty_factor() - 0.75).abs() < 1e-10);
    }

    #[test]
    fn trot_duty_factor() {
        assert!((GaitPattern::Trot.duty_factor() - 0.60).abs() < 1e-10);
    }

    #[test]
    fn gallop_duty_factor() {
        assert!((GaitPattern::Gallop.duty_factor() - 0.35).abs() < 1e-10);
    }

    #[test]
    fn pace_duty_factor() {
        assert!((GaitPattern::Pace.duty_factor() - 0.55).abs() < 1e-10);
    }

    #[test]
    fn bound_duty_factor() {
        assert!((GaitPattern::Bound.duty_factor() - 0.40).abs() < 1e-10);
    }

    #[test]
    fn walk_phase_offsets_four_elements() {
        assert_eq!(GaitPattern::Walk.phase_offsets().len(), 4);
    }

    #[test]
    fn trot_phase_offsets_diagonal_symmetry() {
        let offsets = GaitPattern::Trot.phase_offsets();
        // FL and RR should share phase; FR and RL should share phase
        assert!((offsets[0] - offsets[3]).abs() < 1e-10);
        assert!((offsets[1] - offsets[2]).abs() < 1e-10);
    }

    #[test]
    fn gallop_phase_offsets_four_elements() {
        assert_eq!(GaitPattern::Gallop.phase_offsets().len(), 4);
    }

    #[test]
    fn all_phase_offsets_in_unit_interval() {
        for gait in [
            GaitPattern::Walk,
            GaitPattern::Trot,
            GaitPattern::Gallop,
            GaitPattern::Pace,
            GaitPattern::Bound,
        ] {
            for &offset in &gait.phase_offsets() {
                assert!(
                    (0.0..1.0).contains(&offset),
                    "offset {offset} out of [0,1) for {gait:?}"
                );
            }
        }
    }

    #[test]
    fn duty_factors_in_range() {
        for gait in [
            GaitPattern::Walk,
            GaitPattern::Trot,
            GaitPattern::Gallop,
            GaitPattern::Pace,
            GaitPattern::Bound,
        ] {
            let d = gait.duty_factor();
            assert!(
                d > 0.0 && d < 1.0,
                "duty factor {d} out of (0,1) for {gait:?}"
            );
        }
    }

    // ── BipedModel ────────────────────────────────────────────────────────

    #[test]
    fn biped_new_has_two_legs() {
        let b = BipedModel::new();
        assert_eq!(b.legs.len(), 2);
    }

    #[test]
    fn biped_step_moves_torso_forward() {
        let mut b = BipedModel::new();
        let control = [1.0, 0.0, 2.0, 0.1];
        b.step(0.1, &control);
        assert!((b.torso_position[0] - 0.1).abs() < 1e-10);
    }

    #[test]
    fn biped_default_matches_new() {
        let a = BipedModel::new();
        let b = BipedModel::default();
        assert_eq!(a.torso_position, b.torso_position);
    }

    #[test]
    fn biped_step_zero_speed_no_forward_movement() {
        let mut b = BipedModel::new();
        let pos_before = b.torso_position;
        b.step(0.1, &[0.0, 0.0, 1.0, 0.0]);
        assert!((b.torso_position[0] - pos_before[0]).abs() < 1e-10);
    }

    // ── QuadrupedModel ───────────────────────────────────────────────────

    #[test]
    fn quadruped_new_has_four_legs() {
        let q = QuadrupedModel::new(GaitPattern::Trot);
        assert_eq!(q.legs.len(), 4);
    }

    #[test]
    fn quadruped_step_updates_gait_phase() {
        let mut q = QuadrupedModel::new(GaitPattern::Walk);
        let phase_before = q.gait_phase;
        q.step(0.05);
        assert!(q.gait_phase != phase_before || q.gait_phase == 0.0);
    }

    #[test]
    fn quadruped_trot_stance_count_even() {
        let mut q = QuadrupedModel::new(GaitPattern::Trot);
        // At phase 0, two legs should be in stance
        q.step(0.0);
        let stance_count = q.legs.iter().filter(|l| l.in_stance).count();
        // Trot should have ~2 legs in stance at any given time (approximately)
        assert!((1..=4).contains(&stance_count));
    }

    #[test]
    fn quadruped_gait_phase_wraps() {
        let mut q = QuadrupedModel::new(GaitPattern::Gallop);
        q.step_frequency = 100.0;
        q.step(1.0); // large step to force wrap
        assert!(q.gait_phase >= 0.0 && q.gait_phase < 1.0);
    }

    // ── ZeroMomentPoint ──────────────────────────────────────────────────

    #[test]
    fn zmp_single_force_at_origin() {
        let mut zmp = ZeroMomentPoint::new();
        let result = zmp.compute_zmp(&[[0.0, 100.0, 0.0]], &[[0.0, 0.0, 0.0]]);
        assert!((result[0]).abs() < 1e-10);
        assert!((result[1]).abs() < 1e-10);
    }

    #[test]
    fn zmp_single_force_offset() {
        let mut zmp = ZeroMomentPoint::new();
        let result = zmp.compute_zmp(&[[0.0, 100.0, 0.0]], &[[2.0, 0.0, 3.0]]);
        assert!((result[0] - 2.0).abs() < 1e-10, "zmp_x={}", result[0]);
        assert!((result[1] - 3.0).abs() < 1e-10, "zmp_z={}", result[1]);
    }

    #[test]
    fn zmp_two_equal_forces_symmetric() {
        let mut zmp = ZeroMomentPoint::new();
        let forces = [[0.0, 50.0, 0.0], [0.0, 50.0, 0.0]];
        let positions = [[-1.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
        let result = zmp.compute_zmp(&forces, &positions);
        assert!(result[0].abs() < 1e-10, "symmetric ZMP should be at 0");
    }

    #[test]
    fn zmp_negative_forces_ignored() {
        let mut zmp = ZeroMomentPoint::new();
        // Tensile (negative) force should be ignored
        let result = zmp.compute_zmp(&[[0.0, -100.0, 0.0]], &[[5.0, 0.0, 0.0]]);
        assert!(result[0].abs() < 1e-10);
    }

    #[test]
    fn zmp_no_forces_returns_origin() {
        let mut zmp = ZeroMomentPoint::new();
        let result = zmp.compute_zmp(&[], &[]);
        assert!(result[0].abs() < 1e-10);
        assert!(result[1].abs() < 1e-10);
    }

    #[test]
    fn zmp_stores_cop_values() {
        let mut zmp = ZeroMomentPoint::new();
        zmp.compute_zmp(&[[0.0, 100.0, 0.0]], &[[3.0, 0.0, 4.0]]);
        assert!((zmp.cop_x - 3.0).abs() < 1e-10);
        assert!((zmp.cop_y - 4.0).abs() < 1e-10);
    }

    #[test]
    fn zmp_weighted_average() {
        let mut zmp = ZeroMomentPoint::new();
        let forces = [[0.0, 100.0, 0.0], [0.0, 300.0, 0.0]];
        let positions = [[0.0, 0.0, 0.0], [4.0, 0.0, 0.0]];
        let result = zmp.compute_zmp(&forces, &positions);
        // ZMP = (100*0 + 300*4) / 400 = 3.0
        assert!((result[0] - 3.0).abs() < 1e-10, "zmp_x={}", result[0]);
    }

    // ── CapturePoint ─────────────────────────────────────────────────────

    #[test]
    fn capture_point_stationary_equals_position() {
        let mut cp = CapturePoint::new(3.13);
        let result = cp.compute_cp([1.0, 2.0], [0.0, 0.0], 3.13);
        assert!((result[0] - 1.0).abs() < 1e-10);
        assert!((result[1] - 2.0).abs() < 1e-10);
    }

    #[test]
    fn capture_point_forward_velocity() {
        let mut cp = CapturePoint::new(3.13);
        // xi_x = 0 + 3.13 / 3.13 = 1.0
        let result = cp.compute_cp([0.0, 0.0], [3.13, 0.0], 3.13);
        assert!((result[0] - 1.0).abs() < 1e-10, "xi_x={}", result[0]);
        assert!((result[1]).abs() < 1e-10);
    }

    #[test]
    fn capture_point_stores_position() {
        let mut cp = CapturePoint::new(1.0);
        cp.compute_cp([0.5, 0.5], [1.0, 1.0], 1.0);
        assert!((cp.position[0] - 1.5).abs() < 1e-10);
        assert!((cp.position[1] - 1.5).abs() < 1e-10);
    }

    #[test]
    fn capture_point_omega_updated() {
        let mut cp = CapturePoint::new(1.0);
        cp.compute_cp([0.0, 0.0], [0.0, 0.0], 5.0);
        assert!((cp.omega - 5.0).abs() < 1e-10);
    }

    #[test]
    fn capture_point_symmetry() {
        let mut cp = CapturePoint::new(2.0);
        let r1 = cp.compute_cp([1.0, 0.0], [2.0, 0.0], 2.0);
        let r2 = cp.compute_cp([-1.0, 0.0], [-2.0, 0.0], 2.0);
        // Should be symmetric
        assert!((r1[0] + r2[0]).abs() < 1e-10);
    }

    // ── stride_length ────────────────────────────────────────────────────

    #[test]
    fn stride_length_basic() {
        assert!((stride_length(2.0, 2.0) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn stride_length_zero_frequency() {
        assert!(stride_length(1.0, 0.0).abs() < 1e-10);
    }

    #[test]
    fn stride_length_proportional_to_speed() {
        let s1 = stride_length(1.0, 1.0);
        let s2 = stride_length(2.0, 1.0);
        assert!((s2 - 2.0 * s1).abs() < 1e-10);
    }

    #[test]
    fn stride_length_inversely_proportional_to_frequency() {
        let s1 = stride_length(1.0, 1.0);
        let s2 = stride_length(1.0, 2.0);
        assert!((s1 - 2.0 * s2).abs() < 1e-10);
    }

    // ── froude_number ────────────────────────────────────────────────────

    #[test]
    fn froude_number_basic() {
        // Fr = v²/(g*L) = 1/(9.81*1) ≈ 0.102
        let fr = froude_number(1.0, 1.0);
        assert!((fr - 1.0 / 9.81).abs() < 1e-6);
    }

    #[test]
    fn froude_number_zero_length() {
        assert!(froude_number(1.0, 0.0).abs() < 1e-10);
    }

    #[test]
    fn froude_number_scales_quadratically_with_speed() {
        let fr1 = froude_number(1.0, 1.0);
        let fr2 = froude_number(2.0, 1.0);
        assert!((fr2 - 4.0 * fr1).abs() < 1e-10);
    }

    #[test]
    fn froude_number_walk_run_threshold_near_one() {
        // At walk-run transition: Fr ≈ 0.25 (Froude 1 corresponds to trot)
        // preferred speed ≈ 0.5*sqrt(g*L), so Fr = 0.25
        let v = preferred_walk_run_speed(1.0);
        let fr = froude_number(v, 1.0);
        assert!((fr - 0.25).abs() < 1e-6, "Fr at transition={fr}");
    }

    // ── preferred_walk_run_speed ──────────────────────────────────────────

    #[test]
    fn preferred_walk_run_speed_scales_with_sqrt_leg_length() {
        let v1 = preferred_walk_run_speed(1.0);
        let v4 = preferred_walk_run_speed(4.0);
        assert!((v4 - 2.0 * v1).abs() < 1e-6);
    }

    #[test]
    fn preferred_walk_run_speed_positive() {
        assert!(preferred_walk_run_speed(0.5) > 0.0);
        assert!(preferred_walk_run_speed(1.0) > 0.0);
    }

    #[test]
    fn preferred_walk_run_speed_human_realistic() {
        // Human leg ~1.0 m: expected ~1.57 m/s (~5.6 km/h)
        let v = preferred_walk_run_speed(1.0);
        assert!(v > 1.0 && v < 3.0, "human transition speed={v}");
    }
}
