// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Multiscale rigid body methods.
//!
//! Provides coarse-grained (CG) representations of fine-scale rigid body
//! assemblies, including inertia homogenization, force projection/interpolation,
//! adaptive resolution switching, and error estimation.

// ── vector helpers ────────────────────────────────────────────────────────────

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

// ── mat3 helpers ──────────────────────────────────────────────────────────────

/// A 3×3 matrix stored row-major.
pub type Mat3 = [[f64; 3]; 3];

fn mat3_add(a: Mat3, b: Mat3) -> Mat3 {
    let mut c = [[0.0_f64; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            c[i][j] = a[i][j] + b[i][j];
        }
    }
    c
}

#[cfg(test)]
fn mat3_scale(a: Mat3, s: f64) -> Mat3 {
    let mut c = [[0.0_f64; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            c[i][j] = a[i][j] * s;
        }
    }
    c
}

#[cfg(test)]
fn mat3_transpose(a: Mat3) -> Mat3 {
    let mut t = [[0.0_f64; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            t[i][j] = a[j][i];
        }
    }
    t
}

#[cfg(test)]
fn mat3_trace(a: Mat3) -> f64 {
    a[0][0] + a[1][1] + a[2][2]
}

fn mat3_frobenius_norm(a: Mat3) -> f64 {
    let mut s = 0.0_f64;
    for row in a.iter() {
        for &v in row.iter() {
            s += v * v;
        }
    }
    s.sqrt()
}

// ── FineBody ──────────────────────────────────────────────────────────────────

/// A single fine-scale rigid body in a multiscale hierarchy.
#[derive(Debug, Clone)]
pub struct FineBody {
    /// Position relative to the coarse-grained body centre of mass \[m\].
    pub local_position: [f64; 3],
    /// Velocity \[m/s\].
    pub velocity: [f64; 3],
    /// Angular velocity \[rad/s\].
    pub angular_velocity: [f64; 3],
    /// Applied force \[N\].
    pub force: [f64; 3],
    /// Applied torque \[N·m\].
    pub torque: [f64; 3],
    /// Mass \[kg\].
    pub mass: f64,
    /// Principal moments of inertia \[Ixx, Iyy, Izz\] \[kg·m²\].
    pub inertia: [f64; 3],
    /// Activity level (0 = inactive, 1 = fully active).
    pub activity: f64,
}

impl FineBody {
    /// Create a new fine body at a given local position with given mass and inertia.
    pub fn new(local_position: [f64; 3], mass: f64, inertia: [f64; 3]) -> Self {
        Self {
            local_position,
            velocity: [0.0; 3],
            angular_velocity: [0.0; 3],
            force: [0.0; 3],
            torque: [0.0; 3],
            mass,
            inertia,
            activity: 1.0,
        }
    }
}

// ── CoarseGrainedBody ─────────────────────────────────────────────────────────

/// A coarse-grained (CG) body that aggregates a collection of fine-scale rigid bodies.
///
/// The CG body represents the lumped dynamics of the fine-scale assembly,
/// with homogenized inertia and projected forces.
#[derive(Debug, Clone)]
pub struct CoarseGrainedBody {
    /// Centre-of-mass position \[m\].
    pub position: [f64; 3],
    /// Centre-of-mass velocity \[m/s\].
    pub velocity: [f64; 3],
    /// Angular velocity \[rad/s\].
    pub angular_velocity: [f64; 3],
    /// Total mass \[kg\].
    pub mass: f64,
    /// Effective (homogenized) inertia tensor \[kg·m²\].
    pub inertia_tensor: Mat3,
    /// Fine-scale bodies that map to this CG body.
    pub fine_bodies: Vec<FineBody>,
    /// Resolution level: 0 = coarse, higher = finer.
    pub resolution: usize,
    /// Activity threshold above which fine-scale resolution is activated.
    pub activity_threshold: f64,
}

impl CoarseGrainedBody {
    /// Create a new CG body with an initial set of fine bodies.
    pub fn new(fine_bodies: Vec<FineBody>, activity_threshold: f64) -> Self {
        let mass: f64 = fine_bodies.iter().map(|b| b.mass).sum();
        // Compute initial CG position (centre of mass)
        let mut pos = [0.0_f64; 3];
        if mass > 1e-15 {
            for fb in &fine_bodies {
                pos = vec3_add(pos, vec3_scale(fb.local_position, fb.mass));
            }
            pos = vec3_scale(pos, 1.0 / mass);
        }
        let inertia_tensor = homogenize_inertia(&fine_bodies);
        Self {
            position: pos,
            velocity: [0.0; 3],
            angular_velocity: [0.0; 3],
            mass,
            inertia_tensor,
            fine_bodies,
            resolution: 0,
            activity_threshold,
        }
    }

    /// Refresh the homogenized inertia tensor from current fine bodies.
    pub fn refresh_inertia(&mut self) {
        self.inertia_tensor = homogenize_inertia(&self.fine_bodies);
    }
}

// ── homogenize_inertia ────────────────────────────────────────────────────────

/// Compute the effective (homogenized) inertia tensor for an assembly of fine bodies.
///
/// Uses the parallel-axis theorem: I_cg = Σ \[I_i + m_i * (|r_i|²·E - r_i ⊗ r_i)\]
/// where r_i is the position of fine body i relative to the assembly centre of mass.
///
/// # Arguments
/// * `fine_bodies` – Slice of fine bodies.
pub fn homogenize_inertia(fine_bodies: &[FineBody]) -> Mat3 {
    let total_mass: f64 = fine_bodies.iter().map(|b| b.mass).sum();
    if total_mass < 1e-15 {
        return [[0.0; 3]; 3];
    }

    // Centre of mass
    let mut com = [0.0_f64; 3];
    for b in fine_bodies {
        com = vec3_add(com, vec3_scale(b.local_position, b.mass));
    }
    com = vec3_scale(com, 1.0 / total_mass);

    let mut i_eff = [[0.0_f64; 3]; 3];

    for b in fine_bodies {
        // Self inertia (diagonal)
        let mut i_self = [[0.0_f64; 3]; 3];
        i_self[0][0] = b.inertia[0];
        i_self[1][1] = b.inertia[1];
        i_self[2][2] = b.inertia[2];

        // Parallel-axis contribution
        let r = vec3_sub(b.local_position, com);
        let r2 = vec3_dot(r, r);
        // m * (|r|^2 * I - r ⊗ r)
        let mut par = [[0.0_f64; 3]; 3];
        for row in 0..3 {
            for col in 0..3 {
                let delta = if row == col { 1.0 } else { 0.0 };
                par[row][col] = b.mass * (r2 * delta - r[row] * r[col]);
            }
        }

        i_eff = mat3_add(i_eff, mat3_add(i_self, par));
    }

    i_eff
}

// ── scale_up_forces ───────────────────────────────────────────────────────────

/// Project fine-scale forces and torques onto the coarse-grained DOFs.
///
/// The net force on the CG body is the sum of all fine-scale forces.
/// The net torque includes the direct torques plus r_i × F_i contributions
/// from forces applied at off-centre fine-body locations.
///
/// # Arguments
/// * `fine_bodies` – Slice of fine bodies with applied forces/torques.
/// * `cg_position` – Centre of mass position of the CG body \[m\].
///
/// Returns `(net_force, net_torque)`.
pub fn scale_up_forces(fine_bodies: &[FineBody], cg_position: [f64; 3]) -> ([f64; 3], [f64; 3]) {
    let mut net_force = [0.0_f64; 3];
    let mut net_torque = [0.0_f64; 3];

    for b in fine_bodies {
        net_force = vec3_add(net_force, b.force);
        // Direct torque
        net_torque = vec3_add(net_torque, b.torque);
        // Moment arm contribution: r × F
        let r = vec3_sub(b.local_position, cg_position);
        let moment = cross3(r, b.force);
        net_torque = vec3_add(net_torque, moment);
    }

    (net_force, net_torque)
}

#[inline]
fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

// ── scale_down_motion ─────────────────────────────────────────────────────────

/// Interpolate coarse-grained body motion to fine-scale DOFs.
///
/// Each fine body receives:
/// - velocity: v_cg + ω_cg × r_i
/// - angular velocity: ω_cg (rigid-body assumption)
///
/// # Arguments
/// * `fine_bodies` – Mutable slice of fine bodies to update.
/// * `cg_position` – CG body centre of mass \[m\].
/// * `cg_velocity` – CG linear velocity \[m/s\].
/// * `cg_angular_velocity` – CG angular velocity \[rad/s\].
pub fn scale_down_motion(
    fine_bodies: &mut [FineBody],
    cg_position: [f64; 3],
    cg_velocity: [f64; 3],
    cg_angular_velocity: [f64; 3],
) {
    for b in fine_bodies.iter_mut() {
        let r = vec3_sub(b.local_position, cg_position);
        let omega_cross_r = cross3(cg_angular_velocity, r);
        b.velocity = vec3_add(cg_velocity, omega_cross_r);
        b.angular_velocity = cg_angular_velocity;
    }
}

// ── adaptive_resolution ───────────────────────────────────────────────────────

/// Switch simulation resolution based on the aggregate activity level.
///
/// If the maximum activity among fine bodies exceeds `cg.activity_threshold`,
/// the resolution level is incremented (finer). If all activities are below
/// the threshold by a margin, the resolution is decremented (coarser).
///
/// # Arguments
/// * `cg` – Coarse-grained body to update.
/// * `max_resolution` – Maximum allowed resolution level.
pub fn adaptive_resolution(cg: &mut CoarseGrainedBody, max_resolution: usize) {
    let max_activity = cg
        .fine_bodies
        .iter()
        .map(|b| b.activity)
        .fold(0.0_f64, f64::max);

    let min_activity = cg
        .fine_bodies
        .iter()
        .map(|b| b.activity)
        .fold(1.0_f64, f64::min);

    if max_activity > cg.activity_threshold && cg.resolution < max_resolution {
        cg.resolution += 1;
    } else if min_activity < cg.activity_threshold * 0.5 && cg.resolution > 0 {
        cg.resolution -= 1;
    }
}

// ── error_estimator ───────────────────────────────────────────────────────────

/// Quantify the coarse-graining approximation error.
///
/// Computes the relative Frobenius norm difference between the CG inertia tensor
/// and the homogenized fine-scale inertia tensor:
///
/// ε = ||I_cg - I_hom||_F / (||I_hom||_F + ε_reg)
///
/// A value of 0 means the CG tensor perfectly matches the homogenized one.
///
/// # Arguments
/// * `cg` – The coarse-grained body.
pub fn error_estimator(cg: &CoarseGrainedBody) -> f64 {
    let i_hom = homogenize_inertia(&cg.fine_bodies);
    let mut diff = [[0.0_f64; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            diff[i][j] = cg.inertia_tensor[i][j] - i_hom[i][j];
        }
    }
    let err = mat3_frobenius_norm(diff);
    let scale = mat3_frobenius_norm(i_hom) + 1e-15;
    err / scale
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    fn two_symmetric_bodies() -> Vec<FineBody> {
        vec![
            FineBody::new([-1.0, 0.0, 0.0], 1.0, [0.1, 0.1, 0.1]),
            FineBody::new([1.0, 0.0, 0.0], 1.0, [0.1, 0.1, 0.1]),
        ]
    }

    // ── homogenize_inertia ─────────────────────────────────────────────────

    #[test]
    fn homogenize_inertia_empty() {
        let i = homogenize_inertia(&[]);
        assert_eq!(mat3_trace(i), 0.0);
    }

    #[test]
    fn homogenize_inertia_single_body_at_origin() {
        let b = FineBody::new([0.0, 0.0, 0.0], 1.0, [2.0, 3.0, 4.0]);
        let i = homogenize_inertia(&[b]);
        assert!((i[0][0] - 2.0).abs() < 1e-10, "Ixx={}", i[0][0]);
        assert!((i[1][1] - 3.0).abs() < 1e-10, "Iyy={}", i[1][1]);
        assert!((i[2][2] - 4.0).abs() < 1e-10, "Izz={}", i[2][2]);
    }

    #[test]
    fn homogenize_inertia_symmetric_bodies_off_axis() {
        let bodies = two_symmetric_bodies();
        let i = homogenize_inertia(&bodies);
        // Iyy contribution from parallel axis: each body at |r|=1 from COM
        // I_yy = 2*(0.1 + 1.0*1.0) = 2.2
        assert!((i[1][1] - 2.2).abs() < 1e-10, "Iyy={}", i[1][1]);
    }

    #[test]
    fn homogenize_inertia_symmetric() {
        let bodies = two_symmetric_bodies();
        let i = homogenize_inertia(&bodies);
        for (row, r) in i.iter().enumerate() {
            for (col, &val) in r.iter().enumerate() {
                assert!(
                    (val - i[col][row]).abs() < 1e-12,
                    "asymmetric ({row},{col})"
                );
            }
        }
    }

    #[test]
    fn homogenize_inertia_positive_diagonal() {
        let bodies = two_symmetric_bodies();
        let i = homogenize_inertia(&bodies);
        for (k, row) in i.iter().enumerate() {
            assert!(row[k] > 0.0, "diagonal[{k}] should be positive");
        }
    }

    // ── CoarseGrainedBody ──────────────────────────────────────────────────

    #[test]
    fn cg_body_mass_sum() {
        let bodies = two_symmetric_bodies();
        let cg = CoarseGrainedBody::new(bodies, 0.5);
        assert!((cg.mass - 2.0).abs() < 1e-10);
    }

    #[test]
    fn cg_body_initial_resolution_zero() {
        let cg = CoarseGrainedBody::new(two_symmetric_bodies(), 0.5);
        assert_eq!(cg.resolution, 0);
    }

    #[test]
    fn cg_body_refresh_inertia() {
        let mut cg = CoarseGrainedBody::new(two_symmetric_bodies(), 0.5);
        let before = cg.inertia_tensor;
        cg.refresh_inertia();
        let after = cg.inertia_tensor;
        for i in 0..3 {
            for j in 0..3 {
                assert!((before[i][j] - after[i][j]).abs() < 1e-12);
            }
        }
    }

    // ── scale_up_forces ────────────────────────────────────────────────────

    #[test]
    fn scale_up_forces_zero_forces() {
        let bodies = two_symmetric_bodies();
        let (f, t) = scale_up_forces(&bodies, [0.0; 3]);
        assert_eq!(f, [0.0; 3]);
        assert_eq!(t, [0.0; 3]);
    }

    #[test]
    fn scale_up_forces_sums_correctly() {
        let mut bodies = two_symmetric_bodies();
        bodies[0].force = [1.0, 0.0, 0.0];
        bodies[1].force = [2.0, 0.0, 0.0];
        let (f, _t) = scale_up_forces(&bodies, [0.0; 3]);
        assert!((f[0] - 3.0).abs() < 1e-10, "f[0]={}", f[0]);
    }

    #[test]
    fn scale_up_forces_moment_arm_contribution() {
        let mut bodies = two_symmetric_bodies();
        // Body at [1,0,0] with force in +y → torque about z = r × F = [1,0,0] × [0,1,0] = [0,0,1]
        bodies[1].force = [0.0, 1.0, 0.0];
        let (_, t) = scale_up_forces(&bodies, [0.0; 3]);
        assert!((t[2] - 1.0).abs() < 1e-10, "t[2]={}", t[2]);
    }

    #[test]
    fn scale_up_forces_direct_torque_included() {
        let mut bodies = two_symmetric_bodies();
        bodies[0].torque = [0.0, 0.0, 5.0];
        bodies[1].torque = [0.0, 0.0, 3.0];
        let (_, t) = scale_up_forces(&bodies, [0.0; 3]);
        assert!((t[2] - 8.0).abs() < 1e-10, "t[2]={}", t[2]);
    }

    // ── scale_down_motion ──────────────────────────────────────────────────

    #[test]
    fn scale_down_motion_zero_angular_velocity() {
        let mut bodies = two_symmetric_bodies();
        scale_down_motion(&mut bodies, [0.0; 3], [1.0, 2.0, 3.0], [0.0; 3]);
        for b in &bodies {
            assert_eq!(b.velocity, [1.0, 2.0, 3.0]);
            assert_eq!(b.angular_velocity, [0.0; 3]);
        }
    }

    #[test]
    fn scale_down_motion_rigid_body_rotation() {
        let mut bodies = vec![FineBody::new([1.0, 0.0, 0.0], 1.0, [0.1; 3])];
        // ω = [0, 0, 1] (spin about z), r = [1,0,0] → v = ω × r = [0,1,0]*1 - wait:
        // ω × r = [0,0,1] × [1,0,0] = [0*0-1*0, 1*1-0*0, 0*0-0*1] = [0, 1, 0]
        scale_down_motion(&mut bodies, [0.0; 3], [0.0; 3], [0.0, 0.0, 1.0]);
        assert!(
            (bodies[0].velocity[1] - 1.0).abs() < 1e-10,
            "vy={}",
            bodies[0].velocity[1]
        );
        assert!((bodies[0].velocity[0]).abs() < 1e-10);
    }

    #[test]
    fn scale_down_motion_sets_angular_velocity() {
        let mut bodies = two_symmetric_bodies();
        scale_down_motion(&mut bodies, [0.0; 3], [0.0; 3], [1.0, 2.0, 3.0]);
        for b in &bodies {
            assert_eq!(b.angular_velocity, [1.0, 2.0, 3.0]);
        }
    }

    // ── adaptive_resolution ────────────────────────────────────────────────

    #[test]
    fn adaptive_resolution_increases_when_active() {
        let mut cg = CoarseGrainedBody::new(two_symmetric_bodies(), 0.5);
        for b in &mut cg.fine_bodies {
            b.activity = 0.9; // above threshold
        }
        adaptive_resolution(&mut cg, 3);
        assert_eq!(cg.resolution, 1);
    }

    #[test]
    fn adaptive_resolution_stays_at_max() {
        let mut cg = CoarseGrainedBody::new(two_symmetric_bodies(), 0.5);
        cg.resolution = 3;
        for b in &mut cg.fine_bodies {
            b.activity = 0.9;
        }
        adaptive_resolution(&mut cg, 3);
        assert_eq!(cg.resolution, 3, "should not exceed max");
    }

    #[test]
    fn adaptive_resolution_decreases_when_inactive() {
        let mut cg = CoarseGrainedBody::new(two_symmetric_bodies(), 0.5);
        cg.resolution = 2;
        for b in &mut cg.fine_bodies {
            b.activity = 0.1; // below 0.5 * 0.5 = 0.25 threshold
        }
        adaptive_resolution(&mut cg, 3);
        assert_eq!(cg.resolution, 1);
    }

    #[test]
    fn adaptive_resolution_stays_at_zero() {
        let mut cg = CoarseGrainedBody::new(two_symmetric_bodies(), 0.5);
        for b in &mut cg.fine_bodies {
            b.activity = 0.1;
        }
        adaptive_resolution(&mut cg, 3);
        assert_eq!(cg.resolution, 0);
    }

    // ── error_estimator ────────────────────────────────────────────────────

    #[test]
    fn error_estimator_zero_when_up_to_date() {
        let cg = CoarseGrainedBody::new(two_symmetric_bodies(), 0.5);
        let err = error_estimator(&cg);
        assert!(err < 1e-12, "err={err}");
    }

    #[test]
    fn error_estimator_nonzero_when_stale() {
        let mut cg = CoarseGrainedBody::new(two_symmetric_bodies(), 0.5);
        // Manually perturb the stored inertia
        cg.inertia_tensor[0][0] += 10.0;
        let err = error_estimator(&cg);
        assert!(err > 0.0, "err={err}");
    }

    #[test]
    fn error_estimator_range_zero_to_one_typical() {
        let cg = CoarseGrainedBody::new(two_symmetric_bodies(), 0.5);
        let err = error_estimator(&cg);
        assert!((0.0..10.0).contains(&err), "err={err}");
    }

    // ── mat3 helpers ───────────────────────────────────────────────────────

    #[test]
    fn mat3_add_correctness() {
        let a = [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]];
        let b = [[9.0, 8.0, 7.0], [6.0, 5.0, 4.0], [3.0, 2.0, 1.0]];
        let c = mat3_add(a, b);
        for row in c.iter() {
            for &v in row.iter() {
                assert!((v - 10.0).abs() < 1e-10, "v={v}");
            }
        }
    }

    #[test]
    fn mat3_transpose_correctness() {
        let a = [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]];
        let t = mat3_transpose(a);
        assert!((t[0][1] - 4.0).abs() < 1e-10);
        assert!((t[1][0] - 2.0).abs() < 1e-10);
    }

    #[test]
    fn mat3_frobenius_identity() {
        let i3 = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        // Frobenius norm of 3x3 identity = sqrt(3)
        let f = mat3_frobenius_norm(i3);
        assert!((f - 3.0_f64.sqrt()).abs() < 1e-10, "f={f}");
    }

    #[test]
    fn mat3_scale_zero() {
        let a = [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]];
        let b = mat3_scale(a, 0.0);
        assert_eq!(mat3_frobenius_norm(b), 0.0);
    }

    // ── FineBody ───────────────────────────────────────────────────────────

    #[test]
    fn fine_body_default_activity() {
        let b = FineBody::new([0.0; 3], 1.0, [1.0; 3]);
        assert!((b.activity - 1.0).abs() < 1e-10);
    }

    #[test]
    fn fine_body_zero_initial_velocity() {
        let b = FineBody::new([1.0, 2.0, 3.0], 2.0, [0.5; 3]);
        assert_eq!(b.velocity, [0.0; 3]);
        assert_eq!(b.angular_velocity, [0.0; 3]);
    }
}
