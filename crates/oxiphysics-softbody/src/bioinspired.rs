// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Bio-inspired soft robotics and locomotion.
//!
//! Models pneumatic actuators, soft tentacles, peristaltic locomotion,
//! cephalopod arm mechanics, crawling robots, and a central pattern generator
//! (CPG) oscillator for rhythmic gait control.

use std::f64::consts::PI;

// ---------------------------------------------------------------------------
// Math helpers
// ---------------------------------------------------------------------------

#[inline]
fn clamp(x: f64, lo: f64, hi: f64) -> f64 {
    x.max(lo).min(hi)
}

#[cfg(test)]
#[inline]
fn norm3(v: [f64; 3]) -> f64 {
    (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
}

#[inline]
fn add3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

#[inline]
fn scale3(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

#[inline]
fn rotate_y(v: [f64; 3], angle: f64) -> [f64; 3] {
    let c = angle.cos();
    let s = angle.sin();
    [c * v[0] + s * v[2], v[1], -s * v[0] + c * v[2]]
}

// ---------------------------------------------------------------------------
// PneumaticActuator
// ---------------------------------------------------------------------------

/// McKibben-style pneumatic actuator (soft muscle).
///
/// Geometry follows the braided McKibben model:
/// * rest length `L₀` = `length`
/// * tube radius  `R₀` = `radius`
/// * wall thickness `t`
/// * elastic modulus `E` (Pa)
pub struct PneumaticActuator {
    /// Rest (zero-pressure) length of the actuator in metres.
    pub length: f64,
    /// Outer radius of the actuator in metres.
    pub radius: f64,
    /// Wall thickness of the elastomeric tube in metres.
    pub wall_thickness: f64,
    /// Young's modulus of the elastomeric material (Pa).
    pub modulus: f64,
}

impl PneumaticActuator {
    /// Create a new pneumatic actuator.
    pub fn new(length: f64, radius: f64, wall_thickness: f64, modulus: f64) -> Self {
        Self {
            length,
            radius,
            wall_thickness,
            modulus,
        }
    }

    /// Predicted elongation (m) at gauge pressure `p` (Pa).
    ///
    /// Uses a simplified pressure-volume energy balance:
    /// `δL = p * π R² * L / (E * 2πR * t)` (thin-wall cylinder).
    pub fn pressure_extension(&self, pressure: f64) -> f64 {
        if self.modulus < 1e-9 || self.wall_thickness < 1e-9 {
            return 0.0;
        }
        let area = PI * self.radius * self.radius;
        let wall_stiffness = self.modulus * 2.0 * PI * self.radius * self.wall_thickness;
        pressure * area * self.length / wall_stiffness
    }

    /// Blocking force (N) at gauge pressure `p` (Pa).
    ///
    /// `F_block = p * π R²` — force when extension is fully constrained.
    pub fn blocking_force(&self, pressure: f64) -> f64 {
        pressure * PI * self.radius * self.radius
    }

    /// Free stroke (m) at maximum operating pressure `p_max` (Pa).
    pub fn free_stroke(&self, pressure_max: f64) -> f64 {
        self.pressure_extension(pressure_max)
    }

    /// Curvature (1/m) of a bending actuator at pressure `p`.
    ///
    /// Modelled as `κ = δL / (R * L)` where `δL` is the differential
    /// elongation of the outer fibre relative to the inner fibre.
    pub fn curvature_actuation(&self, pressure: f64) -> f64 {
        let dl = self.pressure_extension(pressure);
        if self.radius < 1e-15 || self.length < 1e-15 {
            return 0.0;
        }
        dl / (self.radius * self.length)
    }
}

// ---------------------------------------------------------------------------
// TentacleSegment
// ---------------------------------------------------------------------------

/// A single rigid segment of a soft tentacle described by constant-curvature
/// kinematics.
pub struct TentacleSegment {
    /// Rest length of this segment (m).
    pub length: f64,
    /// Signed curvature of this segment (1/m).  Positive bends in the +Y plane.
    pub curvature: f64,
    /// Twist about the segment centreline (rad).
    pub twist: f64,
}

impl TentacleSegment {
    /// Create a new tentacle segment.
    pub fn new(length: f64, curvature: f64, twist: f64) -> Self {
        Self {
            length,
            curvature,
            twist,
        }
    }
}

// ---------------------------------------------------------------------------
// SoftTentacle
// ---------------------------------------------------------------------------

/// Multi-segment constant-curvature soft tentacle.
pub struct SoftTentacle {
    /// Ordered segments from base to tip.
    pub segments: Vec<TentacleSegment>,
    /// Number of segments.
    pub n_seg: usize,
}

impl SoftTentacle {
    /// Create a tentacle with `n` segments each of rest length `seg_length`.
    pub fn new(n: usize, seg_length: f64) -> Self {
        let segments = (0..n)
            .map(|_| TentacleSegment::new(seg_length, 0.0, 0.0))
            .collect();
        Self { segments, n_seg: n }
    }

    /// Compute the 3-D position of the tip using constant-curvature forward
    /// kinematics.  Returns `[x, y, z]` in the base frame.
    pub fn tip_position(&self) -> [f64; 3] {
        let pts = self.forward_kinematics_internal();
        *pts.last().unwrap_or(&[0.0; 3])
    }

    /// Approximate workspace radius (m) — maximum reachable tip distance.
    pub fn workspace_radius(&self) -> f64 {
        self.segments.iter().map(|s| s.length).sum()
    }

    /// Forward kinematics for given per-segment bending angles (rad).
    ///
    /// Returns a list of joint positions `[base, joint1, ..., tip]`.
    pub fn forward_kinematics(&self, angles: &[f64]) -> Vec<[f64; 3]> {
        let mut pos = [0.0f64; 3];
        let mut dir = [0.0f64, 0.0, 1.0]; // initial direction: +Z
        let mut pts = vec![pos];

        for (i, seg) in self.segments.iter().enumerate() {
            let angle = angles.get(i).copied().unwrap_or(0.0);
            // Rotate direction by angle about the Y axis.
            dir = rotate_y(dir, angle);
            // Advance along the (possibly twisted) direction.
            let step = scale3(dir, seg.length);
            pos = add3(pos, step);
            pts.push(pos);
        }
        pts
    }

    /// Internal FK using segment curvatures and lengths.
    fn forward_kinematics_internal(&self) -> Vec<[f64; 3]> {
        let mut pos = [0.0f64; 3];
        let mut dir = [0.0f64, 0.0, 1.0];
        let mut pts = vec![pos];

        for seg in &self.segments {
            let angle = seg.curvature * seg.length;
            dir = rotate_y(dir, angle);
            let step = scale3(dir, seg.length);
            pos = add3(pos, step);
            pts.push(pos);
        }
        pts
    }

    /// Apply gravity: increase curvature of each segment proportional to
    /// the weight of all downstream segments.
    ///
    /// Uses a linear stiffness model `κ = W * L / EI` where `EI` is set to
    /// `1.0` (normalised).  Call with `gravity = 9.81` m/s².
    pub fn apply_gravity(&mut self) {
        let n = self.n_seg;
        for i in 0..n {
            // Weight of all downstream segments (assumed unit mass per metre).
            let downstream_length: f64 = self.segments[i..].iter().map(|s| s.length).sum();
            // Simple bending due to gravity.
            self.segments[i].curvature += 0.01 * downstream_length;
        }
    }

    /// Advance the tentacle by one timestep `dt` under the given actuation
    /// signals (one per segment, normalised `[-1, 1]`).
    pub fn step(&mut self, actuation: &[f64], dt: f64) {
        for (i, seg) in self.segments.iter_mut().enumerate() {
            let act = actuation.get(i).copied().unwrap_or(0.0);
            // First-order actuation dynamics.
            let target_curvature = act * 10.0; // max curvature 10 m⁻¹
            seg.curvature += (target_curvature - seg.curvature) * dt * 5.0;
        }
    }
}

// ---------------------------------------------------------------------------
// PeristalticLocomotion
// ---------------------------------------------------------------------------

/// Peristaltic locomotion model for a worm-like soft robot.
///
/// The body is divided into `n` segments, each with a rest length stored in
/// `segments`.  A travelling wave of contraction/extension propagates from
/// head to tail.
pub struct PeristalticLocomotion {
    /// Rest lengths of each body segment (m).
    pub segments: Vec<f64>,
    /// Wave propagation speed (m/s).
    pub wave_speed: f64,
    /// Isotropic friction coefficient.
    pub friction_coeff: f64,
}

impl PeristalticLocomotion {
    /// Create a peristaltic locomotion model.
    ///
    /// # Arguments
    /// * `n_segments`     – number of body segments.
    /// * `segment_length` – rest length of each segment (m).
    /// * `wave_speed`     – wave speed (m/s).
    /// * `friction_coeff` – Coulomb friction coefficient.
    pub fn new(
        n_segments: usize,
        segment_length: f64,
        wave_speed: f64,
        friction_coeff: f64,
    ) -> Self {
        Self {
            segments: vec![segment_length; n_segments],
            wave_speed,
            friction_coeff,
        }
    }

    /// Net forward body displacement (m) at time `t` (s).
    ///
    /// Uses a simplified analytical estimate based on the wave phase:
    /// `x(t) = A * sin(2π f t)` where `A = L_body / (2π)` and
    /// `f = wave_speed / L_body`.
    pub fn body_displacement(&self, t: f64) -> f64 {
        let l_body: f64 = self.segments.iter().sum();
        if l_body < 1e-15 {
            return 0.0;
        }
        let freq = self.wave_speed / l_body;
        let amplitude = l_body / (2.0 * PI);
        amplitude * (2.0 * PI * freq * t).sin()
    }

    /// Average locomotion speed (m/s).
    ///
    /// `v = f * λ * η_friction` where `λ = L_body / n_segments`.
    pub fn locomotion_speed(&self) -> f64 {
        let n = self.segments.len() as f64;
        let l_body: f64 = self.segments.iter().sum();
        if l_body < 1e-15 || n < 1.0 {
            return 0.0;
        }
        let lambda = l_body / n;
        let freq = self.wave_speed / l_body;
        freq * lambda * (1.0 - self.friction_coeff * 0.1)
    }

    /// Locomotion efficiency (dimensionless, 0–1).
    ///
    /// Defined as `η = v_net / v_wave` — ratio of net body speed to wave speed.
    pub fn efficiency(&self) -> f64 {
        if self.wave_speed < 1e-15 {
            return 0.0;
        }
        clamp(self.locomotion_speed() / self.wave_speed, 0.0, 1.0)
    }
}

// ---------------------------------------------------------------------------
// CephalopodArm
// ---------------------------------------------------------------------------

/// Simplified muscle-like model for a cephalopod arm (octopus / squid).
///
/// The arm behaves like a muscular hydrostat: constant volume with transverse
/// muscles and longitudinal muscles coupled by incompressibility.
pub struct CephalopodArm {
    /// Rest (unloaded) length of the arm (m).
    pub rest_length: f64,
    /// Maximum active muscle force (N).
    pub muscle_force: f64,
    /// Passive elastic stiffness (N/m).
    pub passive_stiffness: f64,
}

impl CephalopodArm {
    /// Create a new cephalopod arm model.
    pub fn new(rest_length: f64, muscle_force: f64, passive_stiffness: f64) -> Self {
        Self {
            rest_length,
            muscle_force,
            passive_stiffness,
        }
    }

    /// Predicted elongation (m) under applied tensile force `force` (N).
    ///
    /// `δ = F / k_passive`  (linear elastic model).
    pub fn elongation(&self, force: f64) -> f64 {
        if self.passive_stiffness < 1e-15 {
            return 0.0;
        }
        force / self.passive_stiffness
    }

    /// Muscle contraction ratio — fraction of rest length the arm can shorten.
    ///
    /// Estimated from the Hill-type relation: `ε_max = F_muscle / k_passive / L₀`.
    pub fn contraction_ratio(&self) -> f64 {
        if self.rest_length < 1e-15 || self.passive_stiffness < 1e-15 {
            return 0.0;
        }
        clamp(
            self.muscle_force / (self.passive_stiffness * self.rest_length),
            0.0,
            1.0,
        )
    }
}

// ---------------------------------------------------------------------------
// CrawlingRobot
// ---------------------------------------------------------------------------

/// Simplified 2-D model of a crawling soft robot driven by peristaltic waves
/// or inchworm-like gaits.
pub struct CrawlingRobot {
    /// Bending stiffness of the robot body (N·m).
    pub body_stiffness: f64,
    /// Friction anisotropy ratio: `μ_backward / μ_forward`.
    /// Values > 1 mean the robot resists backward slip.
    pub friction_anisotropy: f64,
}

impl CrawlingRobot {
    /// Create a new crawling robot model.
    pub fn new(body_stiffness: f64, friction_anisotropy: f64) -> Self {
        Self {
            body_stiffness,
            friction_anisotropy,
        }
    }

    /// Anisotropic Coulomb friction force `[Fx, Fy]` (N) for body velocity
    /// `vel = [vx, vy]` and normal force `normal` (N).
    ///
    /// Forward direction is +X.  Backward motion (`vx < 0`) experiences
    /// friction scaled by `friction_anisotropy`.
    pub fn friction_force(&self, vel: [f64; 2], normal: f64) -> [f64; 2] {
        let base_mu = 0.3; // baseline friction coefficient
        let mu_x = if vel[0] < 0.0 {
            base_mu * self.friction_anisotropy
        } else {
            base_mu
        };
        let mu_y = base_mu;

        let fx = if vel[0].abs() > 1e-12 {
            -mu_x * normal * vel[0].signum()
        } else {
            0.0
        };
        let fy = if vel[1].abs() > 1e-12 {
            -mu_y * normal * vel[1].signum()
        } else {
            0.0
        };
        [fx, fy]
    }

    /// Estimate steady-state crawling speed (m/s) for a given actuation
    /// frequency `actuation_freq` (Hz).
    ///
    /// Uses a simple empirical model:
    /// `v = f * L_step * η`  where `L_step = 0.01 m` and
    /// `η = (anisotropy - 1) / anisotropy` (ratchet efficiency).
    pub fn crawl_speed(&self, actuation_freq: f64) -> f64 {
        let l_step = 0.01; // 10 mm per actuation cycle (typical for cm-scale robots)
        let eta = if self.friction_anisotropy > 1.0 {
            (self.friction_anisotropy - 1.0) / self.friction_anisotropy
        } else {
            0.0
        };
        actuation_freq * l_step * eta
    }
}

// ---------------------------------------------------------------------------
// Central Pattern Generator
// ---------------------------------------------------------------------------

/// Simulate one timestep of a Hopf-oscillator-based central pattern generator
/// (CPG).
///
/// The CPG drives rhythmic locomotion gaits.  The oscillator equation is:
///
/// ```text
/// dφ/dt = ω
/// dx/dt = α (μ - r²) x - ω y
/// dy/dt = α (μ - r²) y + ω x
/// ```
///
/// For this simplified version we use an explicit Euler step on the
/// oscillator.
///
/// # Arguments
/// * `phase`     – current oscillator phase (rad).
/// * `freq`      – desired oscillation frequency (Hz).
/// * `amplitude` – desired steady-state amplitude.
/// * `dt`        – timestep (s).
///
/// # Returns
/// `(new_phase, output_signal)` where `output_signal = amplitude * sin(phase)`.
pub fn cpg_oscillator(phase: f64, freq: f64, amplitude: f64, dt: f64) -> (f64, f64) {
    let omega = 2.0 * PI * freq;
    let new_phase = phase + omega * dt;
    let output = amplitude * new_phase.sin();
    (new_phase, output)
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // -----------------------------------------------------------------------
    // PneumaticActuator tests
    // -----------------------------------------------------------------------

    #[test]
    fn pneumatic_extension_zero_at_zero_pressure() {
        let a = PneumaticActuator::new(0.1, 0.01, 0.001, 1e6);
        assert!((a.pressure_extension(0.0)).abs() < 1e-15);
    }

    #[test]
    fn pneumatic_extension_positive_at_positive_pressure() {
        let a = PneumaticActuator::new(0.1, 0.01, 0.001, 1e6);
        assert!(a.pressure_extension(100e3) > 0.0);
    }

    #[test]
    fn pneumatic_extension_scales_linearly_with_pressure() {
        let a = PneumaticActuator::new(0.1, 0.01, 0.001, 1e6);
        let e1 = a.pressure_extension(100e3);
        let e2 = a.pressure_extension(200e3);
        assert!((e2 - 2.0 * e1).abs() < 1e-12, "e1={e1} e2={e2}");
    }

    #[test]
    fn pneumatic_blocking_force_scales_with_pressure() {
        let a = PneumaticActuator::new(0.1, 0.01, 0.001, 1e6);
        let f1 = a.blocking_force(100e3);
        let f2 = a.blocking_force(200e3);
        assert!((f2 - 2.0 * f1).abs() < 1e-8, "f1={f1} f2={f2}");
    }

    #[test]
    fn pneumatic_blocking_force_equals_pressure_times_area() {
        let r = 0.01;
        let a = PneumaticActuator::new(0.1, r, 0.001, 1e6);
        let p = 50e3;
        let expected = p * PI * r * r;
        let got = a.blocking_force(p);
        assert!(
            (got - expected).abs() < 1e-8,
            "got={got} expected={expected}"
        );
    }

    #[test]
    fn pneumatic_free_stroke_equals_extension_at_max_pressure() {
        let a = PneumaticActuator::new(0.1, 0.01, 0.001, 1e6);
        let p_max = 300e3;
        let fs = a.free_stroke(p_max);
        let ext = a.pressure_extension(p_max);
        assert!((fs - ext).abs() < 1e-15);
    }

    #[test]
    fn pneumatic_curvature_zero_at_zero_pressure() {
        let a = PneumaticActuator::new(0.1, 0.01, 0.001, 1e6);
        assert!((a.curvature_actuation(0.0)).abs() < 1e-15);
    }

    #[test]
    fn pneumatic_curvature_positive_at_positive_pressure() {
        let a = PneumaticActuator::new(0.1, 0.01, 0.001, 1e6);
        assert!(a.curvature_actuation(100e3) > 0.0);
    }

    // -----------------------------------------------------------------------
    // SoftTentacle / TentacleSegment tests
    // -----------------------------------------------------------------------

    #[test]
    fn tentacle_tip_at_origin_for_zero_segments() {
        let t = SoftTentacle::new(0, 0.05);
        let tip = t.tip_position();
        assert!(norm3(tip) < 1e-12, "tip = {tip:?}");
    }

    #[test]
    fn tentacle_workspace_radius_correct() {
        let t = SoftTentacle::new(5, 0.1);
        let wr = t.workspace_radius();
        assert!((wr - 0.5).abs() < 1e-10);
    }

    #[test]
    fn tentacle_forward_kinematics_length_correct() {
        let t = SoftTentacle::new(3, 0.1);
        let angles = [0.0; 3];
        let pts = t.forward_kinematics(&angles);
        assert_eq!(pts.len(), 4, "n+1 points");
        let tip = pts[3];
        // Straight: tip should be at z = 0.3
        assert!((tip[2] - 0.3).abs() < 1e-10, "tip z = {}", tip[2]);
    }

    #[test]
    fn tentacle_tip_changes_after_step() {
        let mut t = SoftTentacle::new(4, 0.1);
        let tip_before = t.tip_position();
        let actuation = [1.0, 1.0, 1.0, 1.0];
        t.step(&actuation, 0.1);
        let tip_after = t.tip_position();
        // At least one coordinate should have changed.
        let delta = norm3([
            tip_after[0] - tip_before[0],
            tip_after[1] - tip_before[1],
            tip_after[2] - tip_before[2],
        ]);
        assert!(delta > 1e-8, "tip should move after step: delta={delta}");
    }

    #[test]
    fn tentacle_apply_gravity_increases_curvature() {
        let mut t = SoftTentacle::new(3, 0.1);
        let kappa_before = t.segments[0].curvature;
        t.apply_gravity();
        let kappa_after = t.segments[0].curvature;
        assert!(
            kappa_after > kappa_before,
            "gravity should increase curvature"
        );
    }

    #[test]
    fn tentacle_step_moves_toward_target_curvature() {
        let mut t = SoftTentacle::new(2, 0.1);
        // Set initial curvature to zero, drive with full positive actuation.
        let kappa_before = t.segments[0].curvature;
        let actuation = [1.0, 1.0];
        for _ in 0..20 {
            t.step(&actuation, 0.05);
        }
        let kappa_after = t.segments[0].curvature;
        assert!(
            kappa_after > kappa_before,
            "curvature should increase with positive actuation"
        );
    }

    // -----------------------------------------------------------------------
    // PeristalticLocomotion tests
    // -----------------------------------------------------------------------

    #[test]
    fn peristaltic_locomotion_speed_positive() {
        let p = PeristalticLocomotion::new(10, 0.05, 0.1, 0.3);
        assert!(p.locomotion_speed() > 0.0);
    }

    #[test]
    fn peristaltic_efficiency_between_zero_and_one() {
        let p = PeristalticLocomotion::new(10, 0.05, 0.1, 0.3);
        let e = p.efficiency();
        assert!((0.0..=1.0).contains(&e), "efficiency = {e}");
    }

    #[test]
    fn peristaltic_body_displacement_zero_at_t0() {
        let p = PeristalticLocomotion::new(8, 0.05, 0.1, 0.2);
        assert!((p.body_displacement(0.0)).abs() < 1e-15);
    }

    #[test]
    fn peristaltic_displacement_periodic() {
        let p = PeristalticLocomotion::new(8, 0.05, 0.2, 0.2);
        let l_body: f64 = p.segments.iter().sum();
        let freq = p.wave_speed / l_body;
        let period = 1.0 / freq;
        let d0 = p.body_displacement(0.0);
        let d1 = p.body_displacement(period);
        // After one full period the displacement should return to the initial value.
        assert!((d1 - d0).abs() < 1e-10, "d0={d0} d1={d1}");
    }

    #[test]
    fn peristaltic_speed_increases_with_wave_speed() {
        let p1 = PeristalticLocomotion::new(5, 0.05, 0.1, 0.2);
        let p2 = PeristalticLocomotion::new(5, 0.05, 0.2, 0.2);
        assert!(p2.locomotion_speed() > p1.locomotion_speed());
    }

    // -----------------------------------------------------------------------
    // CephalopodArm tests
    // -----------------------------------------------------------------------

    #[test]
    fn cephalopod_elongation_zero_at_zero_force() {
        let arm = CephalopodArm::new(0.3, 5.0, 100.0);
        assert!((arm.elongation(0.0)).abs() < 1e-15);
    }

    #[test]
    fn cephalopod_elongation_linear_with_force() {
        let arm = CephalopodArm::new(0.3, 5.0, 100.0);
        let e1 = arm.elongation(10.0);
        let e2 = arm.elongation(20.0);
        assert!((e2 - 2.0 * e1).abs() < 1e-12);
    }

    #[test]
    fn cephalopod_contraction_ratio_in_range() {
        let arm = CephalopodArm::new(0.3, 5.0, 100.0);
        let cr = arm.contraction_ratio();
        assert!((0.0..=1.0).contains(&cr), "cr = {cr}");
    }

    #[test]
    fn cephalopod_contraction_ratio_increases_with_muscle_force() {
        let arm1 = CephalopodArm::new(0.3, 5.0, 100.0);
        let arm2 = CephalopodArm::new(0.3, 20.0, 100.0);
        assert!(arm2.contraction_ratio() > arm1.contraction_ratio());
    }

    // -----------------------------------------------------------------------
    // CrawlingRobot tests
    // -----------------------------------------------------------------------

    #[test]
    fn crawling_speed_zero_at_zero_frequency() {
        let r = CrawlingRobot::new(1.0, 2.0);
        assert!((r.crawl_speed(0.0)).abs() < 1e-15);
    }

    #[test]
    fn crawling_speed_positive_with_anisotropy() {
        let r = CrawlingRobot::new(1.0, 3.0);
        assert!(r.crawl_speed(1.0) > 0.0);
    }

    #[test]
    fn crawling_speed_scales_with_frequency() {
        let r = CrawlingRobot::new(1.0, 2.0);
        let v1 = r.crawl_speed(1.0);
        let v2 = r.crawl_speed(2.0);
        assert!((v2 - 2.0 * v1).abs() < 1e-12);
    }

    #[test]
    fn crawling_speed_zero_with_isotropic_friction() {
        let r = CrawlingRobot::new(1.0, 1.0); // anisotropy = 1
        assert!((r.crawl_speed(5.0)).abs() < 1e-15);
    }

    #[test]
    fn crawling_friction_opposes_velocity() {
        let r = CrawlingRobot::new(1.0, 2.0);
        let f = r.friction_force([1.0, 0.0], 10.0);
        assert!(f[0] < 0.0, "friction should oppose +X motion");
    }

    #[test]
    fn crawling_friction_backward_larger_than_forward() {
        let r = CrawlingRobot::new(1.0, 3.0);
        let f_fwd = r.friction_force([1.0, 0.0], 10.0);
        let f_bwd = r.friction_force([-1.0, 0.0], 10.0);
        // |f_bwd| > |f_fwd| due to anisotropy.
        assert!(
            f_bwd[0].abs() > f_fwd[0].abs(),
            "backward friction ({}) should exceed forward ({})",
            f_bwd[0].abs(),
            f_fwd[0].abs()
        );
    }

    #[test]
    fn crawling_friction_zero_for_zero_velocity() {
        let r = CrawlingRobot::new(1.0, 2.0);
        let f = r.friction_force([0.0, 0.0], 10.0);
        assert!(f[0].abs() < 1e-15 && f[1].abs() < 1e-15);
    }

    // -----------------------------------------------------------------------
    // CPG oscillator tests
    // -----------------------------------------------------------------------

    #[test]
    fn cpg_output_bounded_by_amplitude() {
        let amplitude = 2.5;
        let mut phase = 0.0;
        for _ in 0..1000 {
            let (new_phase, out) = cpg_oscillator(phase, 1.0, amplitude, 0.001);
            phase = new_phase;
            assert!(
                out.abs() <= amplitude + 1e-10,
                "CPG output {out} exceeds amplitude {amplitude}"
            );
        }
    }

    #[test]
    fn cpg_phase_advances_at_correct_rate() {
        let freq = 2.0;
        let dt = 0.01;
        let (new_phase, _) = cpg_oscillator(0.0, freq, 1.0, dt);
        let expected = 2.0 * PI * freq * dt;
        assert!(
            (new_phase - expected).abs() < 1e-12,
            "new_phase = {new_phase}"
        );
    }

    #[test]
    fn cpg_zero_amplitude_zero_output() {
        let (_, out) = cpg_oscillator(0.5, 1.0, 0.0, 0.01);
        assert!(out.abs() < 1e-15, "output = {out}");
    }

    #[test]
    fn cpg_zero_frequency_constant_phase() {
        let phase = 0.7;
        let (new_phase, _) = cpg_oscillator(phase, 0.0, 1.0, 0.1);
        assert!((new_phase - phase).abs() < 1e-12);
    }

    #[test]
    fn cpg_output_oscillates_over_time() {
        let mut phase = 0.0;
        let freq = 1.0;
        let amplitude = 1.0;
        let dt = 0.01;
        let mut has_positive = false;
        let mut has_negative = false;
        for _ in 0..300 {
            let (new_phase, out) = cpg_oscillator(phase, freq, amplitude, dt);
            phase = new_phase;
            if out > 0.5 {
                has_positive = true;
            }
            if out < -0.5 {
                has_negative = true;
            }
        }
        assert!(has_positive, "CPG should produce positive outputs");
        assert!(has_negative, "CPG should produce negative outputs");
    }

    // -----------------------------------------------------------------------
    // Misc helper tests
    // -----------------------------------------------------------------------

    #[test]
    fn clamp_within_bounds() {
        assert!((clamp(5.0, 0.0, 10.0) - 5.0).abs() < 1e-15);
        assert!((clamp(-1.0, 0.0, 10.0) - 0.0).abs() < 1e-15);
        assert!((clamp(15.0, 0.0, 10.0) - 10.0).abs() < 1e-15);
    }

    #[test]
    fn rotate_y_quarter_turn() {
        let v = [1.0, 0.0, 0.0];
        let rotated = rotate_y(v, PI / 2.0);
        // [1,0,0] rotated +90° about Y → [0, 0, -1]
        assert!(rotated[0].abs() < 1e-10, "x = {}", rotated[0]);
        assert!(rotated[1].abs() < 1e-10, "y = {}", rotated[1]);
        assert!((rotated[2] + 1.0).abs() < 1e-10, "z = {}", rotated[2]);
    }
}
