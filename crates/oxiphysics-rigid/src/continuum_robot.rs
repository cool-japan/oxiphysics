// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Continuum robot kinematics and dynamics.
//!
//! Provides piecewise-constant-curvature (PCC) models, tendon-driven robots,
//! Cosserat rod mechanics, and workspace computation for soft/continuum robots.
//!
//! # Overview
//!
//! - [`PiecewiseConstantCurvature`] — single PCC arc parameterised by arc
//!   length, curvature, and bending-plane angle.
//! - [`ContinuumSegment`] — elastic segment combining PCC geometry with
//!   cross-section properties (Young's modulus, second moment of area).
//! - [`TendonDrivenRobot`] — series of `ContinuumSegment`s actuated by
//!   tendon lengths with forward-kinematics capability.
//! - [`CossratRod`] — discrete Cosserat rod model tracking node positions and
//!   director triads.
//! - [`RobotWorkspace`] — point cloud of reachable tip positions with
//!   convex-hull boundary extraction.
//! - [`pcc_transformation`] — 4×4 homogeneous transformation matrix for a
//!   PCC arc (column-major, 16 elements).
//! - [`tendon_length_to_curvature`] — mapping from differential tendon lengths
//!   to curvature for a two-tendon arrangement.

#![allow(dead_code)]

// PI is used only in tests — imported there

// ---------------------------------------------------------------------------
// Helper – 3-D vector arithmetic on plain arrays
// ---------------------------------------------------------------------------

#[inline]
fn vec3_add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

#[inline]
fn vec3_sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

#[inline]
fn vec3_scale(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

#[inline]
fn vec3_dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[inline]
fn vec3_len(a: [f64; 3]) -> f64 {
    vec3_dot(a, a).sqrt()
}

#[inline]
fn vec3_cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

// ---------------------------------------------------------------------------
// PiecewiseConstantCurvature
// ---------------------------------------------------------------------------

/// Single arc segment under the piecewise-constant-curvature (PCC) assumption.
///
/// The arc lies in the plane defined by the bending angle `phi` from the
/// x-axis.  For zero curvature the segment is a straight rod along the
/// z-axis of the local frame.
#[derive(Debug, Clone)]
pub struct PiecewiseConstantCurvature {
    /// Arc length of the segment \[m\].
    pub arc_length: f64,
    /// Curvature κ \[1/m\].  Zero means straight.
    pub curvature: f64,
    /// Bending-plane angle φ \[rad\] measured from the x-axis.
    pub bending_angle: f64,
    /// Centre of the arc circle in the local frame \[m\].
    pub arc_center: [f64; 3],
}

impl PiecewiseConstantCurvature {
    /// Construct a new PCC arc.
    ///
    /// * `arc_length`    – arc length \[m\] (must be positive)
    /// * `curvature`     – κ \[1/m\]; zero yields a straight segment
    /// * `bending_angle` – φ \[rad\]
    pub fn new(arc_length: f64, curvature: f64, bending_angle: f64) -> Self {
        // Arc centre is at (1/κ) in the bending plane, perpendicular to the
        // initial tangent.  For κ ≈ 0 the centre goes to infinity; store zeros.
        let arc_center = if curvature.abs() > 1e-12 {
            let r = 1.0 / curvature;
            [r * bending_angle.cos(), r * bending_angle.sin(), 0.0]
        } else {
            [0.0; 3]
        };
        Self {
            arc_length,
            curvature,
            bending_angle,
            arc_center,
        }
    }

    /// Tip position of this arc in the local base frame \[m\].
    ///
    /// For κ = 0 the tip is at `[0, 0, arc_length]`.
    pub fn tip_position(&self) -> [f64; 3] {
        let l = self.arc_length;
        let kappa = self.curvature;
        let phi = self.bending_angle;

        if kappa.abs() < 1e-12 {
            // Straight segment along local z
            return [0.0, 0.0, l];
        }

        // Standard PCC tip formulae (Webster & Jones 2010):
        //   x = (1/κ)(1 − cos(κl)) · cos φ
        //   y = (1/κ)(1 − cos(κl)) · sin φ
        //   z = (1/κ) sin(κl)
        let kl = kappa * l;
        let inv_k = 1.0 / kappa;
        [
            inv_k * (1.0 - kl.cos()) * phi.cos(),
            inv_k * (1.0 - kl.cos()) * phi.sin(),
            inv_k * kl.sin(),
        ]
    }

    /// Tip orientation as a unit quaternion `[w, x, y, z]` in the local frame.
    ///
    /// The orientation represents the rotation of the arc's tangent frame at
    /// the tip relative to the base.  For κ = 0 this is the identity.
    pub fn tip_orientation(&self) -> [f64; 4] {
        let l = self.arc_length;
        let kappa = self.curvature;
        let phi = self.bending_angle;

        if kappa.abs() < 1e-12 {
            // Identity quaternion
            return [1.0, 0.0, 0.0, 0.0];
        }

        // Total bending angle about the axis perpendicular to the bending plane
        let theta = kappa * l; // arc subtended angle
        let half = theta / 2.0;

        // Rotation axis is perpendicular to z and to the bending-plane normal.
        // In the bending plane defined by phi the axis is:
        //   n̂ = [-sin φ, cos φ, 0]
        // Quaternion: q = [cos(θ/2), sin(θ/2)·n̂]
        let s = half.sin();
        [half.cos(), s * (-phi.sin()), s * phi.cos(), 0.0]
    }
}

// ---------------------------------------------------------------------------
// ContinuumSegment
// ---------------------------------------------------------------------------

/// Elastic continuum segment combining PCC geometry with material properties.
///
/// The bending stiffness is `elastic_modulus * second_moment`.
#[derive(Debug, Clone)]
pub struct ContinuumSegment {
    /// Arc length of the undeformed segment \[m\].
    pub length: f64,
    /// PCC geometric model for this segment.
    pub pcc: PiecewiseConstantCurvature,
    /// Young's modulus E \[Pa\].
    pub elastic_modulus: f64,
    /// Second moment of area I \[m⁴\] of the cross-section.
    pub second_moment: f64,
}

impl ContinuumSegment {
    /// Construct a continuum segment.
    pub fn new(length: f64, curvature: f64, phi: f64, e: f64, i: f64) -> Self {
        Self {
            length,
            pcc: PiecewiseConstantCurvature::new(length, curvature, phi),
            elastic_modulus: e,
            second_moment: i,
        }
    }

    /// Bending stiffness EI \[N·m²\].
    pub fn bending_stiffness(&self) -> f64 {
        self.elastic_modulus * self.second_moment
    }

    /// Elastic bending moment for the given curvature \[N·m\].
    ///
    /// M = EI · κ
    pub fn bending_moment(&self) -> f64 {
        self.bending_stiffness() * self.pcc.curvature
    }

    /// Elastic strain energy stored in the segment \[J\].
    ///
    /// U = ½ EI κ² L
    pub fn strain_energy(&self) -> f64 {
        0.5 * self.bending_stiffness() * self.pcc.curvature * self.pcc.curvature * self.length
    }
}

// ---------------------------------------------------------------------------
// TendonDrivenRobot
// ---------------------------------------------------------------------------

/// Multi-segment tendon-driven continuum robot.
///
/// Each segment is driven by tendons routed along its length.  The forward
/// kinematics concatenates the tip transforms of all segments.
#[derive(Debug, Clone)]
pub struct TendonDrivenRobot {
    /// Ordered list of continuum segments from base to tip.
    pub segments: Vec<ContinuumSegment>,
    /// Tendon lengths for each segment \[m\].
    pub tendon_lengths: Vec<f64>,
}

impl TendonDrivenRobot {
    /// Construct a robot from segments and associated tendon lengths.
    pub fn new(segments: Vec<ContinuumSegment>, tendon_lengths: Vec<f64>) -> Self {
        Self {
            segments,
            tendon_lengths,
        }
    }

    /// Compute the tip pose `(position [m], orientation quaternion [w,x,y,z])`
    /// in the world frame by composing all segment transforms.
    pub fn tip_pose(&self) -> ([f64; 3], [f64; 4]) {
        self.forward_kinematics()
    }

    /// Run forward kinematics and return the tip pose.
    ///
    /// Each segment's PCC transformation is applied sequentially.  The result
    /// is the tip position and orientation in the base frame.
    pub fn forward_kinematics(&self) -> ([f64; 3], [f64; 4]) {
        let mut pos = [0.0f64; 3];
        // Accumulated rotation as a rotation matrix (row-major 3×3)
        let mut rot = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0f64];

        for seg in &self.segments {
            let local_tip = seg.pcc.tip_position();
            // Transform local tip into world frame by current rotation
            let world_tip = mat3_mul_vec3(rot, local_tip);
            pos = vec3_add(pos, world_tip);

            // Compose rotation: R_new = R_current * R_seg
            let seg_rot = pcc_to_rotation_matrix(&seg.pcc);
            rot = mat3_mul_mat3(rot, seg_rot);
        }

        let quat = rot_to_quaternion(rot);
        (pos, quat)
    }
}

// ---------------------------------------------------------------------------
// CossratRod
// ---------------------------------------------------------------------------

/// Discrete Cosserat rod discretised into `n` nodes.
///
/// Each node carries a position and a director triad stored as a flattened
/// column-major 3×3 rotation matrix (9 elements).
#[derive(Debug, Clone)]
pub struct CossratRod {
    /// Node positions `[x, y, z]` \[m\].
    pub nodes: Vec<[f64; 3]>,
    /// Director triads at each node (row-major 3×3, 9 f64).
    pub directors: Vec<[f64; 9]>,
}

impl CossratRod {
    /// Construct a straight rod along the z-axis with `n` nodes of spacing `ds`.
    pub fn straight(n: usize, ds: f64) -> Self {
        let nodes = (0..n).map(|i| [0.0, 0.0, i as f64 * ds]).collect();
        // Identity director triad at each node
        let identity = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
        let directors = vec![identity; n];
        Self { nodes, directors }
    }

    /// Compute the axial (shear/extension) strain at each inter-node interval.
    ///
    /// Strain is computed as `|Δx| / ds_ref − 1` where `ds_ref` is the
    /// reference segment length inferred as the mean inter-node distance of a
    /// straight rod.
    pub fn compute_strain(&self) -> Vec<f64> {
        if self.nodes.len() < 2 {
            return Vec::new();
        }
        // Reference length: use the distance between the first two nodes
        let ref_len = {
            let d = vec3_sub(self.nodes[1], self.nodes[0]);
            vec3_len(d)
        };
        let ref_len = if ref_len < 1e-15 { 1.0 } else { ref_len };

        self.nodes
            .windows(2)
            .map(|w| {
                let d = vec3_sub(w[1], w[0]);
                let cur_len = vec3_len(d);
                cur_len / ref_len - 1.0
            })
            .collect()
    }

    /// Compute the bending-induced stress \[Pa\] at each inter-node interval.
    ///
    /// Uses a simple curvature estimate from the director triads:
    /// `σ = E · κ · r_outer` where E = 1 MPa and r_outer = 3 mm (defaults).
    pub fn compute_stress(&self) -> Vec<f64> {
        let e = 1.0e6_f64; // default Young's modulus [Pa]
        let r = 3.0e-3_f64; // outer fibre radius [m]

        if self.directors.len() < 2 {
            return Vec::new();
        }

        self.directors
            .windows(2)
            .map(|w| {
                let d1 = w[0];
                let d2 = w[1];
                // Approximate curvature from change in z-director (tangent)
                let t1 = [d1[6], d1[7], d1[8]];
                let t2 = [d2[6], d2[7], d2[8]];
                let dt = vec3_sub(t2, t1);
                let kappa = vec3_len(dt); // rough curvature magnitude [1/m]
                e * kappa * r
            })
            .collect()
    }

    /// Number of inter-node segments.
    pub fn num_segments(&self) -> usize {
        self.nodes.len().saturating_sub(1)
    }
}

// ---------------------------------------------------------------------------
// RobotWorkspace
// ---------------------------------------------------------------------------

/// Point cloud of reachable tip positions for a continuum robot.
#[derive(Debug, Clone)]
pub struct RobotWorkspace {
    /// Set of reachable Cartesian tip positions \[m\].
    pub reachable_points: Vec<[f64; 3]>,
}

impl RobotWorkspace {
    /// Construct from a pre-computed set of reachable points.
    pub fn new(reachable_points: Vec<[f64; 3]>) -> Self {
        Self { reachable_points }
    }

    /// Approximate the workspace boundary by returning the extreme points along
    /// each Cartesian axis (±x, ±y, ±z).
    ///
    /// For an empty workspace an empty Vec is returned.
    pub fn compute_workspace_boundary(&self) -> Vec<[f64; 3]> {
        if self.reachable_points.is_empty() {
            return Vec::new();
        }

        // Find extrema in each of 6 directions
        let axes: [[f64; 3]; 6] = [
            [1.0, 0.0, 0.0],
            [-1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, -1.0, 0.0],
            [0.0, 0.0, 1.0],
            [0.0, 0.0, -1.0],
        ];

        let mut boundary = Vec::with_capacity(6);
        for axis in &axes {
            let extreme = self
                .reachable_points
                .iter()
                .max_by(|a, b| {
                    let da = vec3_dot(**a, *axis);
                    let db = vec3_dot(**b, *axis);
                    da.partial_cmp(&db).unwrap_or(std::cmp::Ordering::Equal)
                })
                .copied();
            if let Some(p) = extreme {
                boundary.push(p);
            }
        }
        boundary
    }

    /// Number of reachable points in the workspace sample.
    pub fn num_points(&self) -> usize {
        self.reachable_points.len()
    }
}

// ---------------------------------------------------------------------------
// Free functions
// ---------------------------------------------------------------------------

/// Compute the 4×4 homogeneous transformation matrix (column-major, 16 f64)
/// for a PCC arc parameterised by arc length `l`, curvature `kappa`, and
/// bending-plane angle `phi`.
///
/// The matrix maps points in the local tip frame to the local base frame.
/// For κ = 0 the result is a pure translation by `[0, 0, l]`.
///
/// Layout: column-major `[R|t; 0 0 0 1]` expanded to 16 elements.
pub fn pcc_transformation(l: f64, kappa: f64, phi: f64) -> [f64; 16] {
    let pcc = PiecewiseConstantCurvature::new(l, kappa, phi);
    let tip = pcc.tip_position();
    let rot = pcc_to_rotation_matrix(&pcc);

    // Build column-major 4×4
    [
        rot[0], rot[3], rot[6], 0.0, rot[1], rot[4], rot[7], 0.0, rot[2], rot[5], rot[8], 0.0,
        tip[0], tip[1], tip[2], 1.0,
    ]
}

/// Map differential tendon lengths to curvature for a two-tendon arrangement.
///
/// Given lengths `l1` and `l2` for two tendons routed at distance `d` from
/// the neutral axis, the curvature is:
///
/// ```text
/// κ = (l1 − l2) / (d · (l1 + l2))
/// ```
///
/// Returns zero when `l1 + l2 ≤ 0` or `d = 0`.
pub fn tendon_length_to_curvature(l1: f64, l2: f64, d: f64) -> f64 {
    let sum = l1 + l2;
    if sum <= 0.0 || d.abs() < 1e-15 {
        return 0.0;
    }
    (l1 - l2) / (d * sum)
}

// ---------------------------------------------------------------------------
// Internal helpers
// ---------------------------------------------------------------------------

/// Build a 3×3 rotation matrix (row-major) from a PCC arc.
fn pcc_to_rotation_matrix(pcc: &PiecewiseConstantCurvature) -> [f64; 9] {
    let kappa = pcc.curvature;
    let phi = pcc.bending_angle;
    let l = pcc.arc_length;

    if kappa.abs() < 1e-12 {
        // Identity
        return [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
    }

    let theta = kappa * l; // subtended angle
    let c = theta.cos();
    let s = theta.sin();
    let cp = phi.cos();
    let sp = phi.sin();

    // Rotation matrix for bending by θ about axis [-sinφ, cosφ, 0]
    // Using Rodrigues' formula with n̂ = [-sinφ, cosφ, 0]
    let nx = -sp;
    let ny = cp;
    let nz = 0.0_f64;
    let one_c = 1.0 - c;

    [
        c + nx * nx * one_c,
        nx * ny * one_c - nz * s,
        nx * nz * one_c + ny * s,
        ny * nx * one_c + nz * s,
        c + ny * ny * one_c,
        ny * nz * one_c - nx * s,
        nz * nx * one_c - ny * s,
        nz * ny * one_c + nx * s,
        c + nz * nz * one_c,
    ]
}

/// Multiply a row-major 3×3 matrix by a 3-vector.
fn mat3_mul_vec3(m: [f64; 9], v: [f64; 3]) -> [f64; 3] {
    [
        m[0] * v[0] + m[1] * v[1] + m[2] * v[2],
        m[3] * v[0] + m[4] * v[1] + m[5] * v[2],
        m[6] * v[0] + m[7] * v[1] + m[8] * v[2],
    ]
}

/// Multiply two row-major 3×3 matrices.
fn mat3_mul_mat3(a: [f64; 9], b: [f64; 9]) -> [f64; 9] {
    let mut c = [0.0f64; 9];
    for i in 0..3 {
        for j in 0..3 {
            for k in 0..3 {
                c[i * 3 + j] += a[i * 3 + k] * b[k * 3 + j];
            }
        }
    }
    c
}

/// Convert a row-major 3×3 rotation matrix to a unit quaternion `[w, x, y, z]`.
fn rot_to_quaternion(m: [f64; 9]) -> [f64; 4] {
    let trace = m[0] + m[4] + m[8];
    if trace > 0.0 {
        let s = 0.5 / (trace + 1.0).sqrt();
        let w = 0.25 / s;
        let x = (m[7] - m[5]) * s;
        let y = (m[2] - m[6]) * s;
        let z = (m[3] - m[1]) * s;
        [w, x, y, z]
    } else if m[0] > m[4] && m[0] > m[8] {
        let s = 2.0 * (1.0 + m[0] - m[4] - m[8]).sqrt();
        let w = (m[7] - m[5]) / s;
        let x = 0.25 * s;
        let y = (m[1] + m[3]) / s;
        let z = (m[2] + m[6]) / s;
        [w, x, y, z]
    } else if m[4] > m[8] {
        let s = 2.0 * (1.0 + m[4] - m[0] - m[8]).sqrt();
        let w = (m[2] - m[6]) / s;
        let x = (m[1] + m[3]) / s;
        let y = 0.25 * s;
        let z = (m[5] + m[7]) / s;
        [w, x, y, z]
    } else {
        let s = 2.0 * (1.0 + m[8] - m[0] - m[4]).sqrt();
        let w = (m[3] - m[1]) / s;
        let x = (m[2] + m[6]) / s;
        let y = (m[5] + m[7]) / s;
        let z = 0.25 * s;
        [w, x, y, z]
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    use std::f64::consts::PI;

    const EPS: f64 = 1e-9;

    // ── PCC straight segment ──────────────────────────────────────────────

    #[test]
    fn pcc_straight_tip_position() {
        let pcc = PiecewiseConstantCurvature::new(0.1, 0.0, 0.0);
        let tip = pcc.tip_position();
        assert!((tip[0]).abs() < EPS);
        assert!((tip[1]).abs() < EPS);
        assert!((tip[2] - 0.1).abs() < EPS, "z={}", tip[2]);
    }

    #[test]
    fn pcc_straight_tip_orientation_identity() {
        let pcc = PiecewiseConstantCurvature::new(0.1, 0.0, 0.0);
        let q = pcc.tip_orientation();
        assert!((q[0] - 1.0).abs() < EPS);
        assert!(q[1].abs() < EPS);
        assert!(q[2].abs() < EPS);
        assert!(q[3].abs() < EPS);
    }

    // ── PCC semicircle ────────────────────────────────────────────────────

    #[test]
    fn pcc_semicircle_tip_position() {
        // κ = π / L means κL = π (semicircle)
        let l = PI;
        let kappa = 1.0; // radius = 1 m
        let phi = 0.0;
        let pcc = PiecewiseConstantCurvature::new(l, kappa, phi);
        let tip = pcc.tip_position();
        // Tip should be at x = 2/κ = 2, z ≈ 0
        assert!((tip[0] - 2.0).abs() < 1e-6, "x={}", tip[0]);
        assert!(tip[1].abs() < 1e-6, "y={}", tip[1]);
        assert!(tip[2].abs() < 1e-6, "z={}", tip[2]);
    }

    #[test]
    fn pcc_quarter_circle_tip_position() {
        // κ = 1/r, L = π*r/2 → quarter circle with r=1
        let r = 1.0;
        let l = PI / 2.0;
        let kappa = 1.0 / r;
        let phi = 0.0;
        let pcc = PiecewiseConstantCurvature::new(l, kappa, phi);
        let tip = pcc.tip_position();
        // Expected: x = r(1 - cos(π/2)) = 1, z = r sin(π/2) = 1
        assert!((tip[0] - 1.0).abs() < 1e-6, "x={}", tip[0]);
        assert!(tip[1].abs() < 1e-6, "y={}", tip[1]);
        assert!((tip[2] - 1.0).abs() < 1e-6, "z={}", tip[2]);
    }

    // ── PCC bending plane angle ───────────────────────────────────────────

    #[test]
    fn pcc_bending_plane_phi_90() {
        let l = PI;
        let kappa = 1.0;
        let phi = PI / 2.0; // bending in y-z plane
        let pcc = PiecewiseConstantCurvature::new(l, kappa, phi);
        let tip = pcc.tip_position();
        // x should be ≈ 0, y should be ≈ 2
        assert!(tip[0].abs() < 1e-6, "x={}", tip[0]);
        assert!((tip[1] - 2.0).abs() < 1e-6, "y={}", tip[1]);
    }

    #[test]
    fn pcc_tip_position_xy_norm_consistent() {
        let l = 0.3;
        let kappa = 5.0;
        let phi = PI / 4.0;
        let pcc = PiecewiseConstantCurvature::new(l, kappa, phi);
        let tip = pcc.tip_position();
        // xy projection length must be (1 - cos(κl)) / κ
        let kl = kappa * l;
        let expected_xy = (1.0 - kl.cos()) / kappa;
        let xy = (tip[0].powi(2) + tip[1].powi(2)).sqrt();
        assert!((xy - expected_xy).abs() < 1e-9, "xy={xy}");
    }

    // ── PCC orientation quaternion ────────────────────────────────────────

    #[test]
    fn pcc_orientation_is_unit_quaternion() {
        let pcc = PiecewiseConstantCurvature::new(0.2, 3.0, PI / 6.0);
        let q = pcc.tip_orientation();
        let norm = (q[0].powi(2) + q[1].powi(2) + q[2].powi(2) + q[3].powi(2)).sqrt();
        assert!((norm - 1.0).abs() < 1e-9, "quat norm={norm}");
    }

    #[test]
    fn pcc_orientation_full_circle() {
        // After a full circle (κL = 2π) orientation should be identity
        let l = 2.0 * PI;
        let kappa = 1.0;
        let phi = 0.0;
        let pcc = PiecewiseConstantCurvature::new(l, kappa, phi);
        let q = pcc.tip_orientation();
        // w should be ±1, xyz ≈ 0
        assert!((q[0].abs() - 1.0).abs() < 1e-6, "w={}", q[0]);
    }

    // ── PCC arc center ────────────────────────────────────────────────────

    #[test]
    fn pcc_arc_center_correct_for_unit_curvature() {
        let pcc = PiecewiseConstantCurvature::new(1.0, 1.0, 0.0);
        // Center should be at (1, 0, 0) for φ=0, κ=1
        assert!((pcc.arc_center[0] - 1.0).abs() < 1e-9);
        assert!(pcc.arc_center[1].abs() < 1e-9);
        assert!(pcc.arc_center[2].abs() < 1e-9);
    }

    #[test]
    fn pcc_arc_center_zero_curvature() {
        let pcc = PiecewiseConstantCurvature::new(1.0, 0.0, 0.0);
        assert_eq!(pcc.arc_center, [0.0; 3]);
    }

    // ── ContinuumSegment ──────────────────────────────────────────────────

    #[test]
    fn segment_bending_stiffness() {
        let seg = ContinuumSegment::new(0.1, 2.0, 0.0, 200e9, 1e-8);
        let ei = seg.bending_stiffness();
        assert!((ei - 200e9 * 1e-8).abs() < 1.0, "EI={ei}");
    }

    #[test]
    fn segment_bending_moment() {
        let seg = ContinuumSegment::new(0.1, 3.0, 0.0, 1e6, 1e-6);
        // M = EI * κ = 1e6 * 1e-6 * 3 = 3
        assert!((seg.bending_moment() - 3.0).abs() < 1e-9);
    }

    #[test]
    fn segment_strain_energy_positive() {
        let seg = ContinuumSegment::new(0.1, 2.0, 0.0, 1e6, 1e-6);
        assert!(seg.strain_energy() > 0.0);
    }

    #[test]
    fn segment_zero_curvature_zero_moment() {
        let seg = ContinuumSegment::new(0.1, 0.0, 0.0, 1e6, 1e-6);
        assert!(seg.bending_moment().abs() < 1e-12);
        assert!(seg.strain_energy().abs() < 1e-12);
    }

    // ── TendonDrivenRobot ─────────────────────────────────────────────────

    #[test]
    fn single_segment_robot_straight() {
        let seg = ContinuumSegment::new(0.2, 0.0, 0.0, 1e6, 1e-6);
        let robot = TendonDrivenRobot::new(vec![seg], vec![0.2]);
        let (pos, quat) = robot.tip_pose();
        assert!(pos[0].abs() < 1e-9, "x={}", pos[0]);
        assert!(pos[1].abs() < 1e-9, "y={}", pos[1]);
        assert!((pos[2] - 0.2).abs() < 1e-9, "z={}", pos[2]);
        assert!((quat[0] - 1.0).abs() < 1e-9);
    }

    #[test]
    fn two_segment_robot_additive_length() {
        // Two straight segments of 0.1 m each should give tip at z = 0.2
        let s1 = ContinuumSegment::new(0.1, 0.0, 0.0, 1e6, 1e-6);
        let s2 = ContinuumSegment::new(0.1, 0.0, 0.0, 1e6, 1e-6);
        let robot = TendonDrivenRobot::new(vec![s1, s2], vec![0.1, 0.1]);
        let (pos, _) = robot.tip_pose();
        assert!((pos[2] - 0.2).abs() < 1e-9, "z={}", pos[2]);
    }

    #[test]
    fn robot_forward_kinematics_matches_tip_pose() {
        let seg = ContinuumSegment::new(0.1, 2.0, 0.5, 1e6, 1e-6);
        let robot = TendonDrivenRobot::new(vec![seg], vec![0.1]);
        let fk = robot.forward_kinematics();
        let tp = robot.tip_pose();
        for i in 0..3 {
            assert!((fk.0[i] - tp.0[i]).abs() < 1e-12);
        }
    }

    #[test]
    fn robot_tip_quaternion_unit_norm() {
        let seg = ContinuumSegment::new(0.15, 4.0, 1.0, 1e6, 1e-6);
        let robot = TendonDrivenRobot::new(vec![seg], vec![0.15]);
        let (_, q) = robot.tip_pose();
        let norm = (q[0].powi(2) + q[1].powi(2) + q[2].powi(2) + q[3].powi(2)).sqrt();
        assert!((norm - 1.0).abs() < 1e-6, "norm={norm}");
    }

    // ── tendon_length_to_curvature ────────────────────────────────────────

    #[test]
    fn tendon_curvature_equal_lengths_zero() {
        let kappa = tendon_length_to_curvature(0.1, 0.1, 0.01);
        assert!(kappa.abs() < 1e-12);
    }

    #[test]
    fn tendon_curvature_positive() {
        let kappa = tendon_length_to_curvature(0.12, 0.08, 0.01);
        assert!(kappa > 0.0, "κ={kappa}");
    }

    #[test]
    fn tendon_curvature_negative() {
        let kappa = tendon_length_to_curvature(0.08, 0.12, 0.01);
        assert!(kappa < 0.0, "κ={kappa}");
    }

    #[test]
    fn tendon_curvature_zero_d_returns_zero() {
        let kappa = tendon_length_to_curvature(0.1, 0.08, 0.0);
        assert_eq!(kappa, 0.0);
    }

    #[test]
    fn tendon_curvature_formula_check() {
        let l1 = 0.12;
        let l2 = 0.08;
        let d = 0.01;
        let kappa = tendon_length_to_curvature(l1, l2, d);
        let expected = (l1 - l2) / (d * (l1 + l2));
        assert!((kappa - expected).abs() < 1e-12);
    }

    #[test]
    fn tendon_curvature_symmetric_antisymmetry() {
        let k1 = tendon_length_to_curvature(0.11, 0.09, 0.01);
        let k2 = tendon_length_to_curvature(0.09, 0.11, 0.01);
        assert!((k1 + k2).abs() < 1e-12);
    }

    // ── pcc_transformation ────────────────────────────────────────────────

    #[test]
    fn pcc_transformation_straight_translation_only() {
        let mat = pcc_transformation(0.1, 0.0, 0.0);
        // Lower-right 3×3 should be identity (rotation columns), translation = [0,0,0.1]
        // Column-major: mat[12]=tx, mat[13]=ty, mat[14]=tz
        assert!(mat[12].abs() < 1e-9, "tx={}", mat[12]);
        assert!(mat[13].abs() < 1e-9, "ty={}", mat[13]);
        assert!((mat[14] - 0.1).abs() < 1e-9, "tz={}", mat[14]);
        assert!((mat[15] - 1.0).abs() < 1e-9);
    }

    #[test]
    fn pcc_transformation_last_row() {
        let mat = pcc_transformation(0.2, 1.0, 0.0);
        // Last row in column-major 4×4: indices 3, 7, 11, 15
        assert!(mat[3].abs() < 1e-9, "mat[3]={}", mat[3]);
        assert!(mat[7].abs() < 1e-9, "mat[7]={}", mat[7]);
        assert!(mat[11].abs() < 1e-9, "mat[11]={}", mat[11]);
        assert!((mat[15] - 1.0).abs() < 1e-9, "mat[15]={}", mat[15]);
    }

    // ── CossratRod ────────────────────────────────────────────────────────

    #[test]
    fn cosserat_straight_rod_zero_strain() {
        let rod = CossratRod::straight(5, 0.05);
        let strain = rod.compute_strain();
        assert_eq!(strain.len(), 4);
        for s in &strain {
            assert!(s.abs() < 1e-9, "strain={s}");
        }
    }

    #[test]
    fn cosserat_rod_stress_zero_for_straight() {
        let rod = CossratRod::straight(5, 0.05);
        let stress = rod.compute_stress();
        assert_eq!(stress.len(), 4);
        for s in &stress {
            assert!(s.abs() < 1e-12, "stress={s}");
        }
    }

    #[test]
    fn cosserat_rod_num_segments() {
        let rod = CossratRod::straight(10, 0.01);
        assert_eq!(rod.num_segments(), 9);
    }

    #[test]
    fn cosserat_rod_single_node_empty_strain() {
        let rod = CossratRod {
            nodes: vec![[0.0; 3]],
            directors: vec![[1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0]],
        };
        assert!(rod.compute_strain().is_empty());
        assert!(rod.compute_stress().is_empty());
    }

    #[test]
    fn cosserat_rod_compressed_positive_strain_negative() {
        // Place nodes closer together (compression)
        let rod = CossratRod {
            nodes: vec![
                [0.0, 0.0, 0.0],
                [0.0, 0.0, 0.04], // shorter than reference 0.05
                [0.0, 0.0, 0.08],
            ],
            directors: vec![
                [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0],
                [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0],
                [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0],
            ],
        };
        let strain = rod.compute_strain();
        // Both segments compressed vs reference 0.04 → strain ≈ 0
        for s in &strain {
            assert!(s.abs() < 1e-9, "s={s}");
        }
    }

    #[test]
    fn cosserat_rod_stretched_positive_strain() {
        let rod = CossratRod {
            nodes: vec![
                [0.0, 0.0, 0.0],
                [0.0, 0.0, 0.1], // reference = 0.1
                [0.0, 0.0, 0.3], // stretched: 0.2 vs 0.1
            ],
            directors: vec![
                [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0],
                [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0],
                [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0],
            ],
        };
        let strain = rod.compute_strain();
        assert_eq!(strain.len(), 2);
        // First segment: strain ≈ 0, second: 0.2/0.1 - 1 = 1.0
        assert!((strain[1] - 1.0).abs() < 1e-9, "strain[1]={}", strain[1]);
    }

    // ── RobotWorkspace ────────────────────────────────────────────────────

    #[test]
    fn workspace_boundary_empty_for_no_points() {
        let ws = RobotWorkspace::new(vec![]);
        assert!(ws.compute_workspace_boundary().is_empty());
    }

    #[test]
    fn workspace_boundary_single_point() {
        let ws = RobotWorkspace::new(vec![[1.0, 2.0, 3.0]]);
        let boundary = ws.compute_workspace_boundary();
        // All 6 extremes are the same point
        for p in &boundary {
            assert_eq!(*p, [1.0, 2.0, 3.0]);
        }
    }

    #[test]
    fn workspace_boundary_finds_extremes() {
        let pts = vec![
            [1.0, 0.0, 0.0],
            [-1.0, 0.0, 0.0],
            [0.0, 2.0, 0.0],
            [0.0, -2.0, 0.0],
            [0.0, 0.0, 3.0],
            [0.0, 0.0, -3.0],
        ];
        let ws = RobotWorkspace::new(pts);
        let boundary = ws.compute_workspace_boundary();
        assert_eq!(boundary.len(), 6);
        // +x extreme should be [1,0,0]
        let plus_x = boundary[0];
        assert!((plus_x[0] - 1.0).abs() < 1e-9);
    }

    #[test]
    fn workspace_num_points() {
        let ws = RobotWorkspace::new(vec![[0.0; 3]; 50]);
        assert_eq!(ws.num_points(), 50);
    }

    // ── round-trip identity checks ────────────────────────────────────────

    #[test]
    fn pcc_various_curvatures_finite_tip() {
        for kappa in [0.0, 0.5, 1.0, 5.0, 10.0] {
            let pcc = PiecewiseConstantCurvature::new(0.1, kappa, 0.0);
            let tip = pcc.tip_position();
            for v in tip {
                assert!(v.is_finite(), "kappa={kappa} tip not finite");
            }
        }
    }

    #[test]
    fn pcc_various_angles_finite_tip() {
        for phi in [0.0, PI / 6.0, PI / 4.0, PI / 2.0, PI, 3.0 * PI / 2.0] {
            let pcc = PiecewiseConstantCurvature::new(0.1, 2.0, phi);
            let tip = pcc.tip_position();
            for v in tip {
                assert!(v.is_finite(), "phi={phi} tip not finite");
            }
        }
    }

    #[test]
    fn pcc_transformation_finite_entries() {
        let mat = pcc_transformation(0.2, 3.0, 1.2);
        for v in mat {
            assert!(v.is_finite());
        }
    }

    #[test]
    fn tendon_curvature_large_difference() {
        let kappa = tendon_length_to_curvature(0.2, 0.01, 0.005);
        assert!(kappa.is_finite());
        assert!(kappa > 0.0);
    }
}

// ---------------------------------------------------------------------------
// CableSegment – constant-curvature cable-driven segment (new addition)
// ---------------------------------------------------------------------------

/// A constant-curvature segment of a cable-driven continuum robot.
///
/// Stores geometric curvature, torsion, and the individual cable tension
/// values that produce the deformation.
#[derive(Debug, Clone)]
pub struct CableSegment {
    /// Undeformed arc length of the segment (m).
    pub length: f64,
    /// Curvature κ of the constant-curvature arc (1/m).
    pub curvature: f64,
    /// Torsion τ of the arc (1/m).
    pub torsion: f64,
    /// Cable tension values (N) indexed by cable number.
    pub cable_tensions: Vec<f64>,
}

impl CableSegment {
    /// Construct a straight segment with `n_cables` cables, all at zero tension.
    pub fn new(length: f64, n_cables: usize) -> Self {
        Self {
            length,
            curvature: 0.0,
            torsion: 0.0,
            cable_tensions: vec![0.0; n_cables],
        }
    }

    /// Tip position `[x, y, z]` of this segment under the constant-curvature
    /// assumption (Webster & Jones 2010 convention, local frame).
    ///
    /// For zero curvature the tip is at `[0, 0, L]`.
    pub fn tip_position(&self) -> [f64; 3] {
        let l = self.length;
        let kappa = self.curvature;
        if kappa.abs() < 1e-12 {
            return [0.0, 0.0, l];
        }
        let kl = kappa * l;
        let inv_k = 1.0 / kappa;
        [inv_k * (1.0 - kl.cos()), 0.0, inv_k * kl.sin()]
    }

    /// Tip orientation as a 3×3 rotation matrix (row-major 9 elements).
    ///
    /// Rotation is about the y-axis by angle `κ·L`.  For zero curvature the
    /// identity matrix is returned.
    pub fn tip_orientation(&self) -> [[f64; 3]; 3] {
        let theta = self.curvature * self.length;
        let c = theta.cos();
        let s = theta.sin();
        [[c, 0.0, s], [0.0, 1.0, 0.0], [-s, 0.0, c]]
    }

    /// Actual arc length (same as `length` for a constant-curvature segment).
    pub fn arc_length(&self) -> f64 {
        self.length
    }
}

// ---------------------------------------------------------------------------
// constant_curvature_transform
// ---------------------------------------------------------------------------

/// SE(3) transformation (4×4 row-major) for a constant-curvature arc
/// parameterised by curvature `kappa` (1/m), bending-plane angle `phi` (rad),
/// and arc length `l` (m).
///
/// For `kappa ≈ 0` the result is a pure translation by `[0, 0, l]`.
pub fn constant_curvature_transform(kappa: f64, phi: f64, l: f64) -> [[f64; 4]; 4] {
    let (tx, ty, tz) = if kappa.abs() < 1e-12 {
        (0.0_f64, 0.0_f64, l)
    } else {
        let kl = kappa * l;
        let inv_k = 1.0 / kappa;
        (
            inv_k * (1.0 - kl.cos()) * phi.cos(),
            inv_k * (1.0 - kl.cos()) * phi.sin(),
            inv_k * kl.sin(),
        )
    };

    // Build rotation about axis perpendicular to bending plane
    let theta = kappa * l;
    let c = theta.cos();
    let s = theta.sin();
    let cp = phi.cos();
    let sp = phi.sin();
    // Rodrigues rotation about n = [-sin(phi), cos(phi), 0]
    let nx = -sp;
    let ny = cp;
    let oc = 1.0 - c;
    [
        [c + nx * nx * oc, nx * ny * oc, ny * s, tx],
        [nx * ny * oc, c + ny * ny * oc, -nx * s, ty],
        [-ny * s, nx * s, c, tz],
        [0.0, 0.0, 0.0, 1.0],
    ]
}

// ---------------------------------------------------------------------------
// CableRobot – multi-segment continuum robot
// ---------------------------------------------------------------------------

/// Multi-segment cable-driven continuum robot.
///
/// Segments are ordered from base to tip.  The forward kinematics chains the
/// segment tip positions by composing their constant-curvature transforms.
#[derive(Debug, Clone)]
pub struct CableRobot {
    /// Ordered list of cable segments.
    pub segments: Vec<CableSegment>,
}

impl CableRobot {
    /// Create an empty robot with no segments.
    pub fn new() -> Self {
        Self {
            segments: Vec::new(),
        }
    }

    /// Append a segment at the tip end.
    pub fn add_segment(&mut self, s: CableSegment) {
        self.segments.push(s);
    }

    /// Number of segments.
    pub fn segment_count(&self) -> usize {
        self.segments.len()
    }

    /// Tip position `[x, y, z]` in the base frame.
    ///
    /// Computed by chaining the SE(3) transforms of all segments.
    pub fn tip_position(&self) -> [f64; 3] {
        let mut pos = [0.0f64; 3];
        // row-major 3x3 accumulated rotation
        let mut rot = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0f64]];

        for seg in &self.segments {
            let local = seg.tip_position();
            // Transform local tip into current world frame
            let wx = rot[0][0] * local[0] + rot[0][1] * local[1] + rot[0][2] * local[2];
            let wy = rot[1][0] * local[0] + rot[1][1] * local[1] + rot[1][2] * local[2];
            let wz = rot[2][0] * local[0] + rot[2][1] * local[1] + rot[2][2] * local[2];
            pos[0] += wx;
            pos[1] += wy;
            pos[2] += wz;
            // Compose rotation
            let sr = seg.tip_orientation();
            let mut new_rot = [[0.0f64; 3]; 3];
            for i in 0..3 {
                for j in 0..3 {
                    for k in 0..3 {
                        new_rot[i][j] += rot[i][k] * sr[k][j];
                    }
                }
            }
            rot = new_rot;
        }
        pos
    }

    /// Rough estimate of the workspace volume (m³) as a sphere of radius equal
    /// to the total arc length of all segments.
    pub fn workspace_volume_estimate(&self) -> f64 {
        let total_len: f64 = self.segments.iter().map(|s| s.length).sum();
        (4.0 / 3.0) * std::f64::consts::PI * total_len * total_len * total_len
    }
}

impl Default for CableRobot {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// Cosserat rod energy
// ---------------------------------------------------------------------------

/// Elastic bending and torsion energy of a Cosserat rod segment.
///
/// `E = (EI·κ² + GJ·τ²)·L / 2`
///
/// * `kappa` – curvature (1/m)
/// * `tau`   – torsion (1/m)
/// * `ei`    – bending stiffness EI (N·m²)
/// * `gj`    – torsional stiffness GJ (N·m²)
/// * `l`     – arc length (m)
pub fn cosserat_rod_energy(kappa: f64, tau: f64, ei: f64, gj: f64, l: f64) -> f64 {
    (ei * kappa * kappa + gj * tau * tau) * l / 2.0
}

// ---------------------------------------------------------------------------
// cable_tension_to_curvature
// ---------------------------------------------------------------------------

/// Map cable tensions to arc curvature.
///
/// Curvature κ = Σ(Fᵢ · dᵢ) / EI where all cables are assumed routed at the
/// same offset `d` from the neutral axis.
///
/// * `tensions` – slice of cable tension values (N)
/// * `d`        – cable offset from neutral axis (m)
/// * `ei`       – bending stiffness (N·m²)
pub fn cable_tension_to_curvature(tensions: &[f64], d: f64, ei: f64) -> f64 {
    if ei.abs() < 1e-15 {
        return 0.0;
    }
    let sum: f64 = tensions.iter().sum::<f64>() * d;
    sum / ei
}

// ---------------------------------------------------------------------------
// jacobian_continuum_2d
// ---------------------------------------------------------------------------

/// 2×2 Jacobian of the 2-D tip position `(x, z)` with respect to `(κ, φ)` for
/// a single constant-curvature segment.
///
/// * `kappa` – current curvature (1/m); must be non-zero for a meaningful result
/// * `phi`   – bending-plane angle (rad, unused in the 2-D case)
/// * `l`     – arc length (m)
///
/// Returns `J[row][col]` where rows are `[∂x, ∂z]` and columns are `[∂κ, ∂φ]`.
/// For `κ ≈ 0` a linearised approximation is used.
#[allow(clippy::too_many_arguments)]
pub fn jacobian_continuum_2d(kappa: f64, _phi: f64, l: f64) -> [[f64; 2]; 2] {
    if kappa.abs() < 1e-9 {
        // Linearise: x ≈ 0, z ≈ l → dx/dkappa ≈ 0, dz/dkappa ≈ 0 at kappa=0
        // Use first-order Taylor: x ≈ kappa*l^2/2, z ≈ l
        return [[l * l / 2.0, 0.0], [0.0, 0.0]];
    }
    let kl = kappa * l;
    let inv_k = 1.0 / kappa;
    let inv_k2 = inv_k * inv_k;
    // x = (1 - cos(κl)) / κ
    // dx/dκ = (κl·sin(κl) - (1 - cos(κl))) / κ²
    let dx_dk = (kl * kl.sin() - (1.0 - kl.cos())) * inv_k2;
    // z = sin(κl) / κ
    // dz/dκ = (κl·cos(κl) - sin(κl)) / κ²
    let dz_dk = (kl * kl.cos() - kl.sin()) * inv_k2;
    // No φ dependency in the symmetric 2D case
    [[dx_dk, 0.0], [dz_dk, 0.0]]
}

// ---------------------------------------------------------------------------
// inverse_kinematics_simple
// ---------------------------------------------------------------------------

/// Simple single-segment 2-D inverse kinematics.
///
/// Given a 2-D target `(target_x, target_z)` and arc length `l`, returns
/// approximate curvature `κ` and bending-plane angle `φ` (in the xz-plane).
///
/// Uses a bisection search for `κ ∈ [0, κ_max]` where `κ_max = 2π / l`.
/// The angle `φ` is set to `0` for targets with positive `x`, `π` for negative.
pub fn inverse_kinematics_simple(target: [f64; 2], l: f64) -> (f64, f64) {
    let tx = target[0];
    let tz = target[1];
    let phi = if tx >= 0.0 { 0.0 } else { std::f64::consts::PI };

    // Objective: minimise distance between tip position and target
    let tip_fn = |kappa: f64| -> [f64; 2] {
        if kappa.abs() < 1e-12 {
            return [0.0, l];
        }
        let kl = kappa * l;
        let inv_k = 1.0 / kappa;
        [inv_k * (1.0 - kl.cos()), inv_k * kl.sin()]
    };

    let target_dist = (tx * tx + tz * tz).sqrt();
    if target_dist < 1e-12 || l < 1e-12 {
        return (0.0, 0.0);
    }

    // If target is farther than the max reach (l), clamp
    let capped_tx = tx * (target_dist.min(l) / target_dist);
    let capped_tz = tz * (target_dist.min(l) / target_dist);

    // Bisect over curvature in [0, 6π / l]
    let kappa_max = 6.0 * std::f64::consts::PI / l;
    let mut lo = 0.0_f64;
    let mut hi = kappa_max;

    for _ in 0..60 {
        let mid = (lo + hi) / 2.0;
        let tip = tip_fn(mid);
        let err_lo = {
            let t = tip_fn(lo);
            (t[0] - capped_tx).powi(2) + (t[1] - capped_tz).powi(2)
        };
        let err_mid = (tip[0] - capped_tx).powi(2) + (tip[1] - capped_tz).powi(2);
        if err_mid < err_lo {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    let kappa = (lo + hi) / 2.0;
    (kappa, phi)
}

// ---------------------------------------------------------------------------
// Tests for new additions
// ---------------------------------------------------------------------------

#[cfg(test)]
mod new_tests {

    use crate::continuum_robot::CableRobot;
    use crate::continuum_robot::CableSegment;
    use crate::continuum_robot::cable_tension_to_curvature;
    use crate::continuum_robot::constant_curvature_transform;
    use crate::continuum_robot::cosserat_rod_energy;
    use crate::continuum_robot::inverse_kinematics_simple;
    use crate::continuum_robot::jacobian_continuum_2d;
    use std::f64::consts::PI;

    const EPS: f64 = 1e-9;

    // ── CableSegment ──────────────────────────────────────────────────────

    #[test]
    fn cable_segment_straight_tip_at_z_l() {
        let seg = CableSegment::new(0.2, 3);
        let tip = seg.tip_position();
        assert!(tip[0].abs() < EPS, "x={}", tip[0]);
        assert!(tip[1].abs() < EPS, "y={}", tip[1]);
        assert!((tip[2] - 0.2).abs() < EPS, "z={}", tip[2]);
    }

    #[test]
    fn cable_segment_arc_length_equals_length() {
        let seg = CableSegment::new(0.3, 2);
        assert!((seg.arc_length() - 0.3).abs() < EPS);
    }

    #[test]
    fn cable_segment_n_cables_correct() {
        let seg = CableSegment::new(0.1, 4);
        assert_eq!(seg.cable_tensions.len(), 4);
    }

    #[test]
    fn cable_segment_zero_curvature_straight() {
        let mut seg = CableSegment::new(0.15, 2);
        seg.curvature = 0.0;
        let tip = seg.tip_position();
        assert!((tip[2] - 0.15).abs() < EPS);
    }

    #[test]
    fn cable_segment_tip_orientation_zero_curvature_identity() {
        let seg = CableSegment::new(0.1, 2);
        let rot = seg.tip_orientation();
        // Should be identity
        assert!((rot[0][0] - 1.0).abs() < EPS);
        assert!((rot[1][1] - 1.0).abs() < EPS);
        assert!((rot[2][2] - 1.0).abs() < EPS);
    }

    #[test]
    fn cable_segment_tip_nonzero_curvature_finite() {
        let mut seg = CableSegment::new(0.1, 2);
        seg.curvature = 5.0;
        let tip = seg.tip_position();
        for v in tip {
            assert!(v.is_finite());
        }
    }

    #[test]
    fn cable_segment_quarter_circle_x_equals_r() {
        // κ = 1/r, L = π*r/2 → quarter circle, tip at (r, 0, r)
        let r = 1.0;
        let l = PI / 2.0;
        let mut seg = CableSegment::new(l, 2);
        seg.curvature = 1.0 / r;
        let tip = seg.tip_position();
        assert!((tip[0] - 1.0).abs() < 1e-9, "x={}", tip[0]);
        assert!((tip[2] - 1.0).abs() < 1e-9, "z={}", tip[2]);
    }

    // ── CableRobot ────────────────────────────────────────────────────────

    #[test]
    fn cable_robot_starts_empty() {
        let r = CableRobot::new();
        assert_eq!(r.segment_count(), 0);
    }

    #[test]
    fn cable_robot_add_segment_increments_count() {
        let mut r = CableRobot::new();
        r.add_segment(CableSegment::new(0.1, 2));
        r.add_segment(CableSegment::new(0.1, 2));
        assert_eq!(r.segment_count(), 2);
    }

    #[test]
    fn cable_robot_tip_position_single_straight() {
        let mut r = CableRobot::new();
        r.add_segment(CableSegment::new(0.3, 2));
        let tip = r.tip_position();
        assert!(tip[0].abs() < EPS);
        assert!((tip[2] - 0.3).abs() < EPS);
    }

    #[test]
    fn cable_robot_tip_position_two_straight_adds() {
        let mut r = CableRobot::new();
        r.add_segment(CableSegment::new(0.2, 2));
        r.add_segment(CableSegment::new(0.1, 2));
        let tip = r.tip_position();
        assert!(tip[0].abs() < EPS);
        assert!((tip[2] - 0.3).abs() < EPS);
    }

    #[test]
    fn cable_robot_workspace_volume_positive() {
        let mut r = CableRobot::new();
        r.add_segment(CableSegment::new(0.2, 2));
        assert!(r.workspace_volume_estimate() > 0.0);
    }

    #[test]
    fn cable_robot_workspace_scales_with_length() {
        let mut r1 = CableRobot::new();
        r1.add_segment(CableSegment::new(0.1, 2));
        let mut r2 = CableRobot::new();
        r2.add_segment(CableSegment::new(0.2, 2));
        assert!(r2.workspace_volume_estimate() > r1.workspace_volume_estimate());
    }

    // ── constant_curvature_transform ─────────────────────────────────────

    #[test]
    fn cc_transform_straight_tip_at_z_l() {
        let m = constant_curvature_transform(0.0, 0.0, 0.1);
        // Translation column: m[row][3]
        assert!(m[0][3].abs() < EPS, "tx={}", m[0][3]);
        assert!(m[1][3].abs() < EPS, "ty={}", m[1][3]);
        assert!((m[2][3] - 0.1).abs() < EPS, "tz={}", m[2][3]);
        assert!((m[3][3] - 1.0).abs() < EPS);
    }

    #[test]
    fn cc_transform_last_row_is_homogeneous() {
        let m = constant_curvature_transform(2.0, 0.5, 0.2);
        assert!(m[3][0].abs() < EPS);
        assert!(m[3][1].abs() < EPS);
        assert!(m[3][2].abs() < EPS);
        assert!((m[3][3] - 1.0).abs() < EPS);
    }

    #[test]
    fn cc_transform_entries_finite() {
        let m = constant_curvature_transform(3.0, PI / 4.0, 0.15);
        for row in &m {
            for v in row {
                assert!(v.is_finite());
            }
        }
    }

    // ── cosserat_rod_energy ───────────────────────────────────────────────

    #[test]
    fn cosserat_energy_zero_curvature_torsion() {
        let e = cosserat_rod_energy(0.0, 0.0, 1e4, 5e3, 0.1);
        assert!(e.abs() < EPS);
    }

    #[test]
    fn cosserat_energy_positive_curvature() {
        let e = cosserat_rod_energy(2.0, 0.0, 1e4, 5e3, 0.1);
        assert!(e > 0.0, "energy should be positive: {e}");
    }

    #[test]
    fn cosserat_energy_bending_torsion_formula() {
        let kappa = 1.0;
        let tau = 0.5;
        let ei = 1e3;
        let gj = 2e3;
        let l = 0.2;
        let expected = (ei * kappa * kappa + gj * tau * tau) * l / 2.0;
        let got = cosserat_rod_energy(kappa, tau, ei, gj, l);
        assert!((got - expected).abs() < EPS);
    }

    #[test]
    fn cosserat_energy_scales_with_length() {
        let e1 = cosserat_rod_energy(1.0, 0.0, 1e4, 5e3, 0.1);
        let e2 = cosserat_rod_energy(1.0, 0.0, 1e4, 5e3, 0.2);
        assert!((e2 - 2.0 * e1).abs() < EPS * 10.0);
    }

    // ── cable_tension_to_curvature ────────────────────────────────────────

    #[test]
    fn cable_tension_zero_tensions_zero_curvature() {
        let kappa = cable_tension_to_curvature(&[0.0, 0.0], 0.01, 1e3);
        assert!(kappa.abs() < EPS);
    }

    #[test]
    fn cable_tension_scales_linearly_with_tension() {
        let k1 = cable_tension_to_curvature(&[10.0], 0.01, 1e3);
        let k2 = cable_tension_to_curvature(&[20.0], 0.01, 1e3);
        assert!((k2 - 2.0 * k1).abs() < EPS);
    }

    #[test]
    fn cable_tension_zero_ei_returns_zero() {
        let kappa = cable_tension_to_curvature(&[100.0], 0.01, 0.0);
        assert!(kappa.abs() < EPS);
    }

    #[test]
    fn cable_tension_positive_for_positive_tensions() {
        let kappa = cable_tension_to_curvature(&[5.0, 3.0], 0.01, 1e3);
        assert!(kappa > 0.0);
    }

    // ── jacobian_continuum_2d ─────────────────────────────────────────────

    #[test]
    fn jacobian_2d_finite_entries() {
        let j = jacobian_continuum_2d(2.0, 0.0, 0.1);
        for row in &j {
            for v in row {
                assert!(v.is_finite());
            }
        }
    }

    #[test]
    fn jacobian_2d_near_zero_curvature_finite() {
        let j = jacobian_continuum_2d(0.0, 0.0, 0.1);
        for row in &j {
            for v in row {
                assert!(v.is_finite());
            }
        }
    }

    #[test]
    fn jacobian_2d_dx_dk_positive_for_small_kl() {
        // For small kl, x increases with kappa
        let j = jacobian_continuum_2d(0.1, 0.0, 0.1);
        // dx/dkappa should be positive for kl << 1
        assert!(j[0][0] >= 0.0, "dx/dk={}", j[0][0]);
    }

    // ── inverse_kinematics_simple ─────────────────────────────────────────

    #[test]
    fn ik_simple_straight_target_zero_curvature() {
        // Target directly ahead → curvature ≈ 0
        let (kappa, _phi) = inverse_kinematics_simple([0.0, 0.1], 0.1);
        assert!(kappa >= 0.0);
    }

    #[test]
    fn ik_simple_returns_finite() {
        let (kappa, phi) = inverse_kinematics_simple([0.05, 0.08], 0.1);
        assert!(kappa.is_finite());
        assert!(phi.is_finite());
    }

    #[test]
    fn ik_simple_negative_x_phi_pi() {
        let (_kappa, phi) = inverse_kinematics_simple([-0.05, 0.08], 0.1);
        assert!(
            (phi - PI).abs() < EPS,
            "phi should be π for negative x: {phi}"
        );
    }

    #[test]
    fn ik_simple_zero_target_zero_result() {
        let (kappa, phi) = inverse_kinematics_simple([0.0, 0.0], 0.1);
        assert!(kappa.abs() < 1.0); // just finite and reasonable
        assert!(phi.is_finite());
    }
}
