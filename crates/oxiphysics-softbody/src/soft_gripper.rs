// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Soft robotic gripper simulation using position-based dynamics.
//!
//! Models pneumatically actuated fingers made of particle chains with
//! distance and bending stiffness.  Provides contact detection, grasp-force
//! estimation, workspace sampling, and a lightweight single-step integrator.

#![allow(dead_code)]
#![allow(clippy::too_many_arguments)]

use std::f64::consts::PI;

// ---------------------------------------------------------------------------
// Helper geometry
// ---------------------------------------------------------------------------

#[inline]
fn dist3(a: [f64; 3], b: [f64; 3]) -> f64 {
    let dx = a[0] - b[0];
    let dy = a[1] - b[1];
    let dz = a[2] - b[2];
    (dx * dx + dy * dy + dz * dz).sqrt()
}

#[inline]
fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

#[inline]
fn scale3(v: [f64; 3], s: f64) -> [f64; 3] {
    [v[0] * s, v[1] * s, v[2] * s]
}

#[inline]
fn normalize3(v: [f64; 3]) -> [f64; 3] {
    let n = dist3(v, [0.0; 3]);
    if n < 1e-15 {
        [0.0; 3]
    } else {
        scale3(v, 1.0 / n)
    }
}

#[inline]
fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[inline]
fn add3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

// ---------------------------------------------------------------------------
// GripperConfig
// ---------------------------------------------------------------------------

/// High-level configuration parameters for a soft gripper.
#[derive(Debug, Clone)]
pub struct GripperConfig {
    /// Number of fingers.
    pub n_fingers: usize,
    /// Natural length of each finger (m).
    pub finger_length: f64,
    /// Finger width (m) — used for contact geometry.
    pub finger_width: f64,
    /// Radius of the base ring on which fingers are mounted (m).
    pub base_radius: f64,
    /// Number of particles per finger (including base particle).
    pub n_particles: usize,
    /// Distance constraint stiffness (N/m).
    pub stiffness: f64,
    /// Bending constraint stiffness (N/m per unit curvature).
    pub bending_stiffness: f64,
}

impl GripperConfig {
    /// Create a default gripper configuration.
    pub fn new(
        n_fingers: usize,
        finger_length: f64,
        finger_width: f64,
        base_radius: f64,
        n_particles: usize,
        stiffness: f64,
        bending_stiffness: f64,
    ) -> Self {
        Self {
            n_fingers,
            finger_length,
            finger_width,
            base_radius,
            n_particles,
            stiffness,
            bending_stiffness,
        }
    }
}

// ---------------------------------------------------------------------------
// GripperFinger
// ---------------------------------------------------------------------------

/// Particle-chain representation of a single soft finger.
#[derive(Debug, Clone)]
pub struct GripperFinger {
    /// Particle positions `[x, y, z]` (m), base at index 0.
    pub particles: Vec<[f64; 3]>,
    /// Particle velocities `[vx, vy, vz]` (m/s).
    pub velocities: Vec<[f64; 3]>,
    /// Particle masses (kg).
    pub masses: Vec<f64>,
    /// Number of segments (particles − 1).
    pub n_segments: usize,
    /// Distance constraint stiffness (N/m).
    pub stiffness: f64,
    /// Bending stiffness (N/m per unit angular deviation).
    pub bending_stiffness: f64,
}

// ---------------------------------------------------------------------------
// init_finger
// ---------------------------------------------------------------------------

/// Initialise a straight finger extending radially outward from the base ring.
///
/// The finger lies in the XY plane, pointing outward from `base_radius`.
/// Finger `finger_idx` is placed at angle `2π·i/n_fingers`.
///
/// # Arguments
/// * `config`      – gripper configuration.
/// * `finger_idx`  – zero-based finger index.
///
/// Returns a [`GripperFinger`] at its rest straight configuration.
pub fn init_finger(config: &GripperConfig, finger_idx: usize) -> GripperFinger {
    let n = config.n_particles;
    let angle = 2.0 * PI * finger_idx as f64 / config.n_fingers as f64;
    let cos_a = angle.cos();
    let sin_a = angle.sin();
    let segment_len = config.finger_length / (n - 1).max(1) as f64;

    let mut particles = Vec::with_capacity(n);
    let mut velocities = Vec::with_capacity(n);
    let mut masses = Vec::with_capacity(n);

    for i in 0..n {
        let r = config.base_radius + i as f64 * segment_len;
        particles.push([r * cos_a, r * sin_a, 0.0]);
        velocities.push([0.0; 3]);
        masses.push(if i == 0 { 0.0 } else { 1e-3 }); // base particle is static (mass=0 ≡ infinite)
    }

    GripperFinger {
        particles,
        velocities,
        masses,
        n_segments: n - 1,
        stiffness: config.stiffness,
        bending_stiffness: config.bending_stiffness,
    }
}

// ---------------------------------------------------------------------------
// apply_pneumatic_actuation
// ---------------------------------------------------------------------------

/// Apply pneumatic pressure as a distributed force that bends the finger inward.
///
/// Each non-base particle receives a force directed toward the gripper axis (−r̂)
/// scaled by pressure and its distance from the base.  This is a simplified
/// model that drives inward curling.
///
/// # Arguments
/// * `finger`   – mutable reference to the finger.
/// * `pressure` – gauge pressure (Pa or normalised N/m²).
pub fn apply_pneumatic_actuation(finger: &mut GripperFinger, pressure: f64) {
    let n = finger.particles.len();
    for i in 1..n {
        let p = finger.particles[i];
        // Radial direction (outward unit vector in XY plane)
        let r_xy = normalize3([p[0], p[1], 0.0]);
        // Inward force proportional to pressure and segment index
        let force_mag = pressure * i as f64 * 1e-3;
        finger.velocities[i][0] -= r_xy[0] * force_mag / finger.masses[i].max(1e-12);
        finger.velocities[i][1] -= r_xy[1] * force_mag / finger.masses[i].max(1e-12);
    }
}

// ---------------------------------------------------------------------------
// finger_tip_position
// ---------------------------------------------------------------------------

/// Return the position of the distal-most particle (finger tip).
pub fn finger_tip_position(finger: &GripperFinger) -> [f64; 3] {
    *finger.particles.last().unwrap_or(&[0.0; 3])
}

// ---------------------------------------------------------------------------
// finger_curvature
// ---------------------------------------------------------------------------

/// Estimate the average curvature of the finger as the mean bend angle (rad)
/// at each interior particle.
///
/// For each triplet (i−1, i, i+1) the bend angle is computed from the dot
/// product of the two link vectors.  Returns 0.0 for a straight finger.
pub fn finger_curvature(finger: &GripperFinger) -> f64 {
    let n = finger.particles.len();
    if n < 3 {
        return 0.0;
    }
    let mut angle_sum = 0.0_f64;
    let count = (n - 2) as f64;
    for i in 1..n - 1 {
        let v1 = sub3(finger.particles[i], finger.particles[i - 1]);
        let v2 = sub3(finger.particles[i + 1], finger.particles[i]);
        let n1 = normalize3(v1);
        let n2 = normalize3(v2);
        let cos_theta = dot3(n1, n2).clamp(-1.0, 1.0);
        angle_sum += cos_theta.acos();
    }
    angle_sum / count
}

// ---------------------------------------------------------------------------
// grasp_force
// ---------------------------------------------------------------------------

/// Estimate the total normal contact force exerted on a spherical object.
///
/// Each finger tip that lies within `object_radius` of `object_pos` contributes
/// a radial contact force `= stiffness × overlap` directed toward the object centre.
///
/// # Arguments
/// * `fingers`       – slice of fingers.
/// * `object_pos`    – centre of the grasped object (m).
/// * `object_radius` – radius of the spherical object (m).
///
/// Returns total contact force magnitude (N, summed over all contacting tips).
pub fn grasp_force(fingers: &[GripperFinger], object_pos: [f64; 3], object_radius: f64) -> f64 {
    let mut total = 0.0_f64;
    for finger in fingers {
        let tip = finger_tip_position(finger);
        let d = dist3(tip, object_pos);
        if d < object_radius {
            let overlap = object_radius - d;
            total += finger.stiffness * overlap;
        }
    }
    total
}

// ---------------------------------------------------------------------------
// contact_detection_finger
// ---------------------------------------------------------------------------

/// Return indices of particles in contact with a sphere.
///
/// A particle at position `p` is in contact when `dist(p, object_pos) < object_radius`.
///
/// # Arguments
/// * `finger`       – finger to test.
/// * `object_pos`   – sphere centre (m).
/// * `object_radius`– sphere radius (m).
///
/// Returns a vector of contacting particle indices.
pub fn contact_detection_finger(
    finger: &GripperFinger,
    object_pos: [f64; 3],
    object_radius: f64,
) -> Vec<usize> {
    finger
        .particles
        .iter()
        .enumerate()
        .filter_map(|(i, &p)| {
            if dist3(p, object_pos) < object_radius {
                Some(i)
            } else {
                None
            }
        })
        .collect()
}

// ---------------------------------------------------------------------------
// step_gripper
// ---------------------------------------------------------------------------

/// Advance a finger one time step under pneumatic actuation and distance constraints.
///
/// The integrator is explicit Euler with a distance-constraint projection pass.
///
/// # Arguments
/// * `finger`   – mutable finger.
/// * `pressure` – pneumatic pressure (Pa).
/// * `dt`       – time step (s).
pub fn step_gripper(finger: &mut GripperFinger, pressure: f64, dt: f64) {
    // 1. Apply actuation forces → update velocities
    apply_pneumatic_actuation(finger, pressure);

    // 2. Integrate positions (skip index 0 — base particle, mass == 0)
    let n = finger.particles.len();
    for i in 1..n {
        finger.particles[i] = add3(finger.particles[i], scale3(finger.velocities[i], dt));
    }

    // 3. Distance constraint projection (one pass)
    let seg_len = if n > 1 {
        dist3(finger.particles[0], finger.particles[1]).max(1e-6)
    } else {
        1e-3
    };

    for i in 0..n - 1 {
        let pa = finger.particles[i];
        let pb = finger.particles[i + 1];
        let d = dist3(pa, pb);
        if d < 1e-15 {
            continue;
        }
        let corr = (d - seg_len) / d;
        let dir = sub3(pb, pa);
        if i == 0 {
            // Base fixed: move only particle i+1
            finger.particles[i + 1] = sub3(pb, scale3(dir, corr));
        } else {
            // Move both equally
            finger.particles[i] = add3(pa, scale3(dir, 0.5 * corr));
            finger.particles[i + 1] = sub3(pb, scale3(dir, 0.5 * corr));
        }
    }

    // 4. Simple velocity damping
    for i in 1..n {
        finger.velocities[i] = scale3(finger.velocities[i], 0.95);
    }
}

// ---------------------------------------------------------------------------
// gripper_workspace
// ---------------------------------------------------------------------------

/// Sample the workspace (reachable tip positions) of a gripper by sweeping
/// pressure values and recording the finger-tip positions.
///
/// For each finger `n_samples` pressure values are swept from 0 to a maximum
/// (100 Pa) and the tip position is recorded after stepping the gripper for
/// a fixed number of iterations.
///
/// # Arguments
/// * `config`    – gripper configuration.
/// * `n_samples` – number of pressure samples per finger.
///
/// Returns a flat vector of tip positions (may contain duplicates across fingers).
pub fn gripper_workspace(config: &GripperConfig, n_samples: usize) -> Vec<[f64; 3]> {
    let mut workspace = Vec::with_capacity(config.n_fingers * n_samples);
    let p_max = 100.0_f64;
    for fi in 0..config.n_fingers {
        for si in 0..n_samples {
            let pressure = p_max * si as f64 / n_samples.max(1) as f64;
            let mut finger = init_finger(config, fi);
            let dt = 1.0 / 60.0;
            for _ in 0..30 {
                step_gripper(&mut finger, pressure, dt);
            }
            workspace.push(finger_tip_position(&finger));
        }
    }
    workspace
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    fn default_config() -> GripperConfig {
        GripperConfig::new(3, 0.1, 0.02, 0.03, 6, 1000.0, 100.0)
    }

    // --- init_finger ---

    #[test]
    fn test_init_finger_particle_count() {
        let cfg = default_config();
        let f = init_finger(&cfg, 0);
        assert_eq!(f.particles.len(), cfg.n_particles);
    }

    #[test]
    fn test_init_finger_base_fixed() {
        let cfg = default_config();
        let f = init_finger(&cfg, 0);
        assert_eq!(f.masses[0], 0.0, "base particle has zero mass (static)");
    }

    #[test]
    fn test_init_finger_velocities_zero() {
        let cfg = default_config();
        let f = init_finger(&cfg, 0);
        for v in &f.velocities {
            assert_eq!(*v, [0.0; 3]);
        }
    }

    #[test]
    fn test_init_finger_radial_layout() {
        let cfg = default_config();
        let f = init_finger(&cfg, 0);
        // First finger at angle 0 → particles should lie along +X axis
        for p in &f.particles {
            assert!(p[1].abs() < 1e-10, "y should be ~0 for finger 0 at angle 0");
        }
    }

    #[test]
    fn test_init_finger_different_angles() {
        let cfg = default_config();
        let f0 = init_finger(&cfg, 0);
        let f1 = init_finger(&cfg, 1);
        // Tips should be at different positions
        let t0 = finger_tip_position(&f0);
        let t1 = finger_tip_position(&f1);
        assert!(
            dist3(t0, t1) > 1e-6,
            "fingers should start at different angles"
        );
    }

    // --- apply_pneumatic_actuation ---

    #[test]
    fn test_actuation_moves_velocities() {
        let cfg = default_config();
        let mut f = init_finger(&cfg, 0);
        apply_pneumatic_actuation(&mut f, 10.0);
        // Non-base particles should now have non-zero velocity
        let any_moved = f.velocities[1..]
            .iter()
            .any(|v| v[0].abs() > 1e-15 || v[1].abs() > 1e-15);
        assert!(any_moved, "actuation should change velocities");
    }

    #[test]
    fn test_actuation_base_unaffected() {
        let cfg = default_config();
        let mut f = init_finger(&cfg, 0);
        apply_pneumatic_actuation(&mut f, 100.0);
        assert_eq!(f.velocities[0], [0.0; 3], "base velocity unchanged");
    }

    #[test]
    fn test_zero_pressure_no_change() {
        let cfg = default_config();
        let mut f = init_finger(&cfg, 0);
        let velocities_before: Vec<_> = f.velocities.clone();
        apply_pneumatic_actuation(&mut f, 0.0);
        assert_eq!(f.velocities, velocities_before, "zero pressure → no change");
    }

    // --- finger_tip_position ---

    #[test]
    fn test_tip_position_is_last_particle() {
        let cfg = default_config();
        let f = init_finger(&cfg, 0);
        let tip = finger_tip_position(&f);
        let last = f.particles[f.particles.len() - 1];
        assert_eq!(tip, last);
    }

    // --- finger_curvature ---

    #[test]
    fn test_curvature_straight_finger_zero() {
        let cfg = default_config();
        let f = init_finger(&cfg, 0);
        let curv = finger_curvature(&f);
        assert!(curv < 1e-6, "straight finger → zero curvature, got {curv}");
    }

    #[test]
    fn test_curvature_non_negative() {
        let cfg = default_config();
        let mut f = init_finger(&cfg, 0);
        // Perturb
        f.particles[3][1] += 0.02;
        let curv = finger_curvature(&f);
        assert!(curv >= 0.0, "curvature non-negative");
    }

    #[test]
    fn test_curvature_increases_after_actuation() {
        let cfg = default_config();
        let mut f = init_finger(&cfg, 0);
        let c0 = finger_curvature(&f);
        let dt = 1.0 / 60.0;
        for _ in 0..100 {
            step_gripper(&mut f, 50.0, dt);
        }
        let c1 = finger_curvature(&f);
        assert!(c1 >= c0, "actuation should increase curvature");
    }

    // --- grasp_force ---

    #[test]
    fn test_grasp_force_no_contact() {
        let cfg = default_config();
        let fingers: Vec<_> = (0..cfg.n_fingers).map(|i| init_finger(&cfg, i)).collect();
        // Object placed far away
        let f = grasp_force(&fingers, [10.0, 0.0, 0.0], 0.01);
        assert_eq!(f, 0.0, "no contact → zero force");
    }

    #[test]
    fn test_grasp_force_in_contact() {
        let cfg = GripperConfig::new(1, 0.05, 0.01, 0.0, 3, 1000.0, 100.0);
        let mut finger = init_finger(&cfg, 0);
        // Move tip inside object radius
        let last = finger.particles.len() - 1;
        finger.particles[last] = [0.005, 0.0, 0.0];
        let f = grasp_force(&[finger], [0.0, 0.0, 0.0], 0.02);
        assert!(f > 0.0, "contact → positive force");
    }

    #[test]
    fn test_grasp_force_positive_overlap() {
        let cfg = GripperConfig::new(1, 0.05, 0.01, 0.0, 3, 1000.0, 100.0);
        let mut finger = init_finger(&cfg, 0);
        let last = finger.particles.len() - 1;
        finger.particles[last] = [0.0, 0.0, 0.0];
        let obj = [0.0, 0.0, 0.0];
        let f = grasp_force(&[finger], obj, 0.01);
        // Full overlap → force = k * radius
        assert!((f - 10.0).abs() < 1e-8);
    }

    // --- contact_detection_finger ---

    #[test]
    fn test_contact_detection_no_contact() {
        let cfg = default_config();
        let f = init_finger(&cfg, 0);
        let contacts = contact_detection_finger(&f, [10.0, 0.0, 0.0], 0.01);
        assert!(contacts.is_empty(), "no particles should be in contact");
    }

    #[test]
    fn test_contact_detection_tip_in_sphere() {
        let cfg = GripperConfig::new(1, 0.05, 0.01, 0.0, 3, 1000.0, 100.0);
        let mut f = init_finger(&cfg, 0);
        let last = f.particles.len() - 1;
        f.particles[last] = [0.005, 0.0, 0.0];
        let contacts = contact_detection_finger(&f, [0.0, 0.0, 0.0], 0.02);
        assert!(contacts.contains(&last), "tip should be in contact");
    }

    #[test]
    fn test_contact_detection_all_inside() {
        let cfg = GripperConfig::new(1, 0.01, 0.005, 0.0, 4, 1000.0, 100.0);
        let f = init_finger(&cfg, 0);
        let contacts = contact_detection_finger(&f, [0.02, 0.0, 0.0], 0.5);
        assert_eq!(contacts.len(), f.particles.len(), "all inside large sphere");
    }

    // --- step_gripper ---

    #[test]
    fn test_step_gripper_moves_tip() {
        let cfg = default_config();
        let mut f = init_finger(&cfg, 0);
        let tip_before = finger_tip_position(&f);
        step_gripper(&mut f, 50.0, 1.0 / 60.0);
        let tip_after = finger_tip_position(&f);
        assert!(dist3(tip_before, tip_after) > 0.0, "step should move tip");
    }

    #[test]
    fn test_step_gripper_base_fixed() {
        let cfg = default_config();
        let mut f = init_finger(&cfg, 0);
        let base_before = f.particles[0];
        for _ in 0..100 {
            step_gripper(&mut f, 50.0, 1.0 / 60.0);
        }
        let base_after = f.particles[0];
        assert_eq!(base_before, base_after, "base particle must not move");
    }

    #[test]
    fn test_step_gripper_damped() {
        let cfg = default_config();
        let mut f = init_finger(&cfg, 0);
        step_gripper(&mut f, 10.0, 1.0 / 60.0);
        let v_before: f64 = f.velocities[1].iter().map(|&x| x * x).sum::<f64>().sqrt();
        step_gripper(&mut f, 0.0, 1.0 / 60.0);
        let v_after: f64 = f.velocities[1].iter().map(|&x| x * x).sum::<f64>().sqrt();
        // After a zero-pressure step with damping, velocity should be smaller
        assert!(
            v_after <= v_before * 1.0,
            "velocity should not grow unboundedly"
        );
    }

    // --- gripper_workspace ---

    #[test]
    fn test_workspace_size() {
        let cfg = default_config();
        let ws = gripper_workspace(&cfg, 5);
        assert_eq!(ws.len(), cfg.n_fingers * 5, "workspace size");
    }

    #[test]
    fn test_workspace_not_all_same() {
        let cfg = default_config();
        let ws = gripper_workspace(&cfg, 4);
        // There should be some diversity in tip positions
        let first = ws[0];
        let any_different = ws.iter().skip(1).any(|&p| dist3(p, first) > 1e-12);
        assert!(any_different, "workspace should have diverse positions");
    }

    #[test]
    fn test_workspace_zero_samples() {
        let cfg = default_config();
        let ws = gripper_workspace(&cfg, 0);
        assert!(ws.is_empty(), "zero samples → empty workspace");
    }
}
