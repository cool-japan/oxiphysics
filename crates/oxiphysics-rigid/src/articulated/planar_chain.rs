// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Exact rigid-body dynamics of a planar serial chain of revolute joints.
//!
//! [`PlanarChain`] is the single model behind the planar-chain functions of
//! `articulated::functions` (`rnea_3d_revolute_z`, `rnea_revolute_chain`,
//! `aba_revolute_chain`, `joint_space_inertia_matrix_planar`,
//! `coriolis_bias_forces`, `gravity_compensation_torques`,
//! `forward_dynamics_full`, `kinetic_energy_planar`), so that all of them
//! describe the same mechanism and agree with each other.
//!
//! # Model
//!
//! * `n` rigid links move in the xy-plane. Joint `i` is a revolute joint about
//!   the z axis at the proximal end of link `i`; joint 0 attaches link 0 to the
//!   fixed base, which sits at the origin.
//! * `q[i]` is the *relative* angle between link `i` and its parent, so link `i`
//!   points along the absolute angle `theta_i = q_0 + ... + q_i` (measured
//!   counter-clockwise from +x) and its distal end is the joint origin of link
//!   `i + 1`.
//! * A link is a rod of length `l_i` whose centre of mass lies at its midpoint.
//!   It has mass `m_i` and a moment of inertia `I_i` about the axis through its
//!   centre of mass parallel to z (`m l^2 / 12` for a uniform rod).
//! * Gravity is the acceleration vector `[gx, gy, gz]`; only the in-plane
//!   components produce joint torques.
//! * Missing entries of the parameter slices default to 1 (mass, length,
//!   inertia) and missing entries of the state/torque slices to 0, so slices
//!   shorter than `n` never panic.
//!
//! # Equations of motion
//!
//! With the generalised coordinates `q` the dynamics are
//! `H(q) q'' + C(q, q') + G(q) = tau`, where
//!
//! * `H_kl = sum_{j >= max(k,l)} ( m_j (c_j - p_k) . (c_j - p_l) + I_j )` is the
//!   joint-space inertia matrix (`p_k` is the origin of joint `k`, `c_j` the
//!   centre of mass of link `j`), which depends on the configuration through
//!   the relative positions of the joints and centres of mass;
//! * `C` collects the Coriolis and centrifugal terms and `G` the gravity
//!   torques.
//!
//! The recursive Newton-Euler algorithm (inverse dynamics), the composite rigid
//! body algorithm (`H`) and the articulated body algorithm (forward dynamics)
//! are all implemented, and cross-checked against each other in the tests.

use nalgebra::{DMatrix, DVector, Matrix3, Vector3};

use crate::error::{Error, Result};

/// Joints whose articulated-body inertia is not above this value are treated as
/// locked by the articulated body algorithm (their acceleration is reported as
/// zero) instead of dividing by (nearly) zero.
const MIN_ARTICULATED_INERTIA: f64 = 1.0e-30;

/// Entry `i` of a state/torque slice, or 0 when the slice is too short.
#[inline]
fn coordinate(values: &[f64], i: usize) -> f64 {
    values.get(i).copied().unwrap_or(0.0)
}

/// Difference of two planar vectors.
#[inline]
fn sub2(a: [f64; 2], b: [f64; 2]) -> [f64; 2] {
    [a[0] - b[0], a[1] - b[1]]
}

/// Dot product of two planar vectors.
#[inline]
fn dot2(a: [f64; 2], b: [f64; 2]) -> f64 {
    a[0] * b[0] + a[1] * b[1]
}

/// z component of the cross product of two planar vectors.
#[inline]
fn cross2(a: [f64; 2], b: [f64; 2]) -> f64 {
    a[0] * b[1] - a[1] * b[0]
}

/// Acceleration of a point at offset `r` from a reference point of a body that
/// rotates with angular velocity `omega` and angular acceleration `alpha`:
/// `a_ref + alpha z x r - omega^2 r`.
#[inline]
fn point_acceleration(a_ref: [f64; 2], alpha: f64, omega: f64, r: [f64; 2]) -> [f64; 2] {
    [
        a_ref[0] - alpha * r[1] - omega * omega * r[0],
        a_ref[1] + alpha * r[0] - omega * omega * r[1],
    ]
}

/// Physical description of a planar serial revolute chain.
#[derive(Debug, Clone)]
pub(super) struct PlanarChain {
    /// Mass of every link \[kg\].
    mass: Vec<f64>,
    /// Length of every link \[m\].
    length: Vec<f64>,
    /// Moment of inertia of every link about its centre of mass \[kg m^2\].
    inertia_com: Vec<f64>,
}

/// Planar kinematics of a chain in one configuration.
#[derive(Debug, Clone)]
struct Pose {
    /// Joint origin `p_i` of every link; entry `n` is the tip of the last link.
    origin: Vec<[f64; 2]>,
    /// Centre of mass `c_i` of every link.
    com: Vec<[f64; 2]>,
}

impl Pose {
    /// Vector from the joint origin of link `i` to its centre of mass.
    #[inline]
    fn com_offset(&self, i: usize) -> [f64; 2] {
        sub2(self.com[i], self.origin[i])
    }

    /// Vector from the joint origin of link `i` to its distal end.
    #[inline]
    fn link_vector(&self, i: usize) -> [f64; 2] {
        sub2(self.origin[i + 1], self.origin[i])
    }
}

/// Per-link working data of the articulated body algorithm. All spatial vectors
/// are `[angular z, linear x, linear y]` quantities expressed about the joint
/// origin of the link with world-aligned axes.
struct AbaLink {
    /// Motion transform from the parent joint origin to this joint origin (a
    /// translation by the parent link vector).
    x: Matrix3<f64>,
    /// Velocity-product (centripetal) acceleration of this joint origin as a
    /// point of the parent link.
    c: Vector3<f64>,
    /// Articulated-body inertia about this joint origin.
    ia: Matrix3<f64>,
    /// Articulated-body bias force about this joint origin.
    pa: Vector3<f64>,
    /// `ia * S`, the coupling of the joint axis with the spatial motion.
    u_vec: Vector3<f64>,
    /// `S^T ia S`, the inertia seen by the joint.
    d: f64,
    /// Net generalised force left for the joint, `tau - S^T pa`.
    u: f64,
}

/// Spatial inertia matrix about a joint origin for a link of mass `m`, inertia
/// `i_com` about its centre of mass, whose centre of mass is at offset `r`.
///
/// Maps the spatial acceleration `[alpha, a_x, a_y]` of the joint origin to the
/// moment about the origin and the force of the Newton-Euler equations (without
/// the centripetal term).
fn spatial_inertia(m: f64, i_com: f64, r: [f64; 2]) -> Matrix3<f64> {
    Matrix3::new(
        i_com + m * dot2(r, r),
        -m * r[1],
        m * r[0],
        -m * r[1],
        m,
        0.0,
        m * r[0],
        0.0,
        m,
    )
}

/// Motion transform that moves the reference point of a spatial acceleration
/// by the planar vector `l` (world-aligned axes).
fn translation_transform(l: [f64; 2]) -> Matrix3<f64> {
    Matrix3::new(1.0, 0.0, 0.0, -l[1], 1.0, 0.0, l[0], 0.0, 1.0)
}

impl PlanarChain {
    /// Chain of `n` links; masses, lengths and centre-of-mass inertias missing
    /// from the slices default to 1.
    pub(super) fn new(n: usize, masses: &[f64], lengths: &[f64], inertias_com: &[f64]) -> Self {
        let take = |values: &[f64]| -> Vec<f64> {
            (0..n)
                .map(|i| values.get(i).copied().unwrap_or(1.0))
                .collect()
        };
        Self {
            mass: take(masses),
            length: take(lengths),
            inertia_com: take(inertias_com),
        }
    }

    /// Chain of `n` uniform slender rods, `I_com = m l^2 / 12`.
    pub(super) fn uniform_rods(n: usize, masses: &[f64], lengths: &[f64]) -> Self {
        let mut chain = Self::new(n, masses, lengths, &[]);
        chain.inertia_com = chain
            .mass
            .iter()
            .zip(&chain.length)
            .map(|(m, l)| m * l * l / 12.0)
            .collect();
        chain
    }

    /// Number of links.
    #[inline]
    fn len(&self) -> usize {
        self.mass.len()
    }

    /// Joint origins and centres of mass for the relative joint angles `q`.
    fn pose(&self, q: &[f64]) -> Pose {
        let n = self.len();
        let mut origin = Vec::with_capacity(n + 1);
        let mut com = Vec::with_capacity(n);
        let mut p = [0.0_f64; 2];
        let mut theta = 0.0_f64;
        origin.push(p);
        for (i, &l) in self.length.iter().enumerate() {
            theta += coordinate(q, i);
            let (sin_t, cos_t) = theta.sin_cos();
            com.push([p[0] + 0.5 * l * cos_t, p[1] + 0.5 * l * sin_t]);
            p = [p[0] + l * cos_t, p[1] + l * sin_t];
            origin.push(p);
        }
        Pose { origin, com }
    }

    /// Inverse dynamics (recursive Newton-Euler): the joint torques that
    /// produce the accelerations `qdd` in state `(q, qd)` under `gravity`.
    ///
    /// Forward pass: absolute angular velocity/acceleration of every link and the
    /// acceleration of every joint origin and centre of mass, with gravity
    /// folded in as a fictitious upward base acceleration `-g`. Backward pass:
    /// the force `F_i = m_i a_ci + F_{i+1}` and the moment about the joint
    /// origin `N_i = I_i alpha_i + r_ci x m_i a_ci + N_{i+1} + L_i x F_{i+1}`
    /// (`L_i` is the link vector); the joint torque is `tau_i = N_i`.
    pub(super) fn inverse_dynamics(
        &self,
        q: &[f64],
        qd: &[f64],
        qdd: &[f64],
        gravity: [f64; 3],
    ) -> Vec<f64> {
        let n = self.len();
        if n == 0 {
            return Vec::new();
        }
        let pose = self.pose(q);

        let mut alpha = vec![0.0_f64; n];
        let mut a_com = vec![[0.0_f64; 2]; n];
        let mut omega_abs = 0.0_f64;
        let mut alpha_abs = 0.0_f64;
        let mut a_origin = [-gravity[0], -gravity[1]];
        for i in 0..n {
            omega_abs += coordinate(qd, i);
            alpha_abs += coordinate(qdd, i);
            alpha[i] = alpha_abs;
            a_com[i] = point_acceleration(a_origin, alpha_abs, omega_abs, pose.com_offset(i));
            a_origin = point_acceleration(a_origin, alpha_abs, omega_abs, pose.link_vector(i));
        }

        let mut torques = vec![0.0_f64; n];
        let mut child_force = [0.0_f64; 2];
        let mut child_moment = 0.0_f64;
        for i in (0..n).rev() {
            let inertial_force = [self.mass[i] * a_com[i][0], self.mass[i] * a_com[i][1]];
            let moment = self.inertia_com[i] * alpha[i]
                + cross2(pose.com_offset(i), inertial_force)
                + child_moment
                + cross2(pose.link_vector(i), child_force);
            child_force = [
                inertial_force[0] + child_force[0],
                inertial_force[1] + child_force[1],
            ];
            child_moment = moment;
            torques[i] = moment;
        }
        torques
    }

    /// Composite terms of the sub-chains `l..n` about their joint origins:
    /// `J_l = sum_{j >= l} (m_j |c_j - p_l|^2 + I_j)` and the first moment
    /// `F_l = sum_{j >= l} m_j (c_j - p_l)`.
    fn composite_terms(&self, pose: &Pose) -> (Vec<f64>, Vec<[f64; 2]>) {
        let n = self.len();
        let mut inertia = vec![0.0_f64; n];
        let mut moment = vec![[0.0_f64; 2]; n];
        for l in 0..n {
            for j in l..n {
                let d = sub2(pose.com[j], pose.origin[l]);
                inertia[l] += self.mass[j] * dot2(d, d) + self.inertia_com[j];
                moment[l][0] += self.mass[j] * d[0];
                moment[l][1] += self.mass[j] * d[1];
            }
        }
        (inertia, moment)
    }

    /// Joint-space inertia matrix `H(q)` (row-major `n x n`), by the composite
    /// rigid body method: for `k <= l`,
    /// `H_kl = J_l + F_l . (p_l - p_k)` with the composite terms of
    /// [`PlanarChain::composite_terms`].
    pub(super) fn mass_matrix(&self, q: &[f64]) -> Vec<f64> {
        let n = self.len();
        let pose = self.pose(q);
        let (inertia, moment) = self.composite_terms(&pose);
        let mut h = vec![0.0_f64; n * n];
        for l in 0..n {
            for k in 0..=l {
                let offset = sub2(pose.origin[l], pose.origin[k]);
                let value = inertia[l] + dot2(moment[l], offset);
                h[k * n + l] = value;
                h[l * n + k] = value;
            }
        }
        h
    }

    /// Gravity torques `G(q)`: the joint torques that hold the chain still
    /// against `gravity`, `G_k = F_k x (-g)` with the first moments of
    /// [`PlanarChain::composite_terms`].
    pub(super) fn gravity_torques(&self, q: &[f64], gravity: [f64; 3]) -> Vec<f64> {
        let pose = self.pose(q);
        let (_, moment) = self.composite_terms(&pose);
        let weight_direction = [-gravity[0], -gravity[1]];
        moment
            .iter()
            .map(|&first_moment| cross2(first_moment, weight_direction))
            .collect()
    }

    /// Coriolis and centrifugal joint forces `C(q, q')`: the inverse dynamics
    /// with zero joint acceleration and no gravity. They do not depend on the
    /// rotational inertias of the links.
    pub(super) fn bias_forces(&self, q: &[f64], qd: &[f64]) -> Vec<f64> {
        self.inverse_dynamics(q, qd, &[], [0.0; 3])
    }

    /// Kinetic energy `0.5 q'^T H(q) q'`, evaluated from the centre-of-mass
    /// velocities `v_cj = sum_{k <= j} q'_k z x (c_j - p_k)` and the absolute
    /// angular velocities.
    pub(super) fn kinetic_energy(&self, q: &[f64], qd: &[f64]) -> f64 {
        let pose = self.pose(q);
        let mut energy = 0.0_f64;
        let mut omega_abs = 0.0_f64;
        for j in 0..self.len() {
            omega_abs += coordinate(qd, j);
            let mut velocity = [0.0_f64; 2];
            for k in 0..=j {
                let d = sub2(pose.com[j], pose.origin[k]);
                let rate = coordinate(qd, k);
                velocity[0] -= rate * d[1];
                velocity[1] += rate * d[0];
            }
            energy += 0.5 * self.mass[j] * dot2(velocity, velocity)
                + 0.5 * self.inertia_com[j] * omega_abs * omega_abs;
        }
        energy
    }

    /// Forward dynamics by the articulated body algorithm (Featherstone), in
    /// `O(n)`: the joint accelerations produced by the joint torques `tau` in
    /// state `(q, qd)` under `gravity`.
    ///
    /// The algorithm is run in planar spatial algebra: motion vectors are
    /// `[alpha, a_x, a_y]` and force vectors `[moment, f_x, f_y]` about a joint
    /// origin, with axes parallel to the world axes, so that moving to the next
    /// joint is a pure translation by the link vector.
    ///
    /// 1. Outward pass: the angular rates, the per-link inertia about the joint
    ///    origin, the centripetal bias force `-m omega^2 r_c` and the
    ///    velocity-product acceleration of the next joint origin.
    /// 2. Inward pass: the articulated-body inertia `I^A` and bias force `p^A`;
    ///    with `U = I^A S`, `D = S^T U`, `u = tau - S^T p^A` the joint is
    ///    eliminated (`I^a = I^A - U U^T / D`, `p^a = p^A + I^a c + U u / D`) and
    ///    the result is added to the parent.
    /// 3. Outward pass: `a' = X a_parent + c`, `q'' = (u - U^T a') / D` and
    ///    `a = a' + S q''`, starting from the base acceleration `-g`.
    ///
    /// A joint whose articulated inertia is not above `1e-30` cannot be
    /// accelerated and is reported with zero acceleration.
    pub(super) fn forward_dynamics(
        &self,
        q: &[f64],
        qd: &[f64],
        tau: &[f64],
        gravity: [f64; 3],
    ) -> Vec<f64> {
        let n = self.len();
        if n == 0 {
            return Vec::new();
        }
        let pose = self.pose(q);

        let mut omega = Vec::with_capacity(n);
        let mut omega_abs = 0.0_f64;
        for i in 0..n {
            omega_abs += coordinate(qd, i);
            omega.push(omega_abs);
        }

        let mut links: Vec<AbaLink> = (0..n)
            .map(|i| {
                let r = pose.com_offset(i);
                let m = self.mass[i];
                let (x, c) = if i == 0 {
                    (Matrix3::identity(), Vector3::zeros())
                } else {
                    let l = pose.link_vector(i - 1);
                    let centripetal = -omega[i - 1] * omega[i - 1];
                    (
                        translation_transform(l),
                        Vector3::new(0.0, centripetal * l[0], centripetal * l[1]),
                    )
                };
                AbaLink {
                    x,
                    c,
                    ia: spatial_inertia(m, self.inertia_com[i], r),
                    pa: Vector3::new(
                        0.0,
                        -m * omega[i] * omega[i] * r[0],
                        -m * omega[i] * omega[i] * r[1],
                    ),
                    u_vec: Vector3::zeros(),
                    d: 0.0,
                    u: 0.0,
                }
            })
            .collect();

        for i in (0..n).rev() {
            let ia = links[i].ia;
            let pa = links[i].pa;
            let u_vec = Vector3::new(ia[(0, 0)], ia[(1, 0)], ia[(2, 0)]);
            let d = ia[(0, 0)];
            let u = coordinate(tau, i) - pa[0];
            links[i].u_vec = u_vec;
            links[i].d = d;
            links[i].u = u;
            if i == 0 {
                continue;
            }
            let (ia_reduced, pa_reduced) = if d <= MIN_ARTICULATED_INERTIA {
                (ia, pa + ia * links[i].c)
            } else {
                let ia_reduced = ia - u_vec * u_vec.transpose() / d;
                let pa_reduced = pa + ia_reduced * links[i].c + u_vec * (u / d);
                (ia_reduced, pa_reduced)
            };
            let x = links[i].x;
            links[i - 1].ia += x.transpose() * ia_reduced * x;
            links[i - 1].pa += x.transpose() * pa_reduced;
        }

        let mut result = Vec::with_capacity(n);
        let mut a_parent = Vector3::new(0.0, -gravity[0], -gravity[1]);
        for link in &links {
            let a_prime = link.x * a_parent + link.c;
            let qdd = if link.d <= MIN_ARTICULATED_INERTIA {
                0.0
            } else {
                (link.u - link.u_vec.dot(&a_prime)) / link.d
            };
            a_parent = a_prime + Vector3::new(qdd, 0.0, 0.0);
            result.push(qdd);
        }
        result
    }

    /// Forward dynamics by solving `H(q) q'' = tau - C(q, q') - G(q)` with a
    /// Cholesky factorisation of the joint-space inertia matrix.
    ///
    /// # Errors
    ///
    /// Returns [`Error::ForwardDynamicsFailure`] when `H` is not symmetric
    /// positive definite (for example zero-mass links or non-finite state) or
    /// when the resulting accelerations are not finite.
    pub(super) fn forward_dynamics_crba(
        &self,
        q: &[f64],
        qd: &[f64],
        tau: &[f64],
        gravity: [f64; 3],
    ) -> Result<Vec<f64>> {
        let n = self.len();
        if n == 0 {
            return Ok(Vec::new());
        }
        let h = self.mass_matrix(q);
        let bias = self.bias_forces(q, qd);
        let gravity_torques = self.gravity_torques(q, gravity);
        let rhs = DVector::from_iterator(
            n,
            (0..n).map(|i| coordinate(tau, i) - bias[i] - gravity_torques[i]),
        );
        let factorisation = DMatrix::from_row_slice(n, n, &h)
            .cholesky()
            .ok_or_else(|| {
                let body_index = (0..n)
                    .find(|&i| !(h[i * n + i].is_finite() && h[i * n + i] > 0.0))
                    .unwrap_or(0);
                Error::forward_dynamics_failure(
                    body_index,
                    "joint-space inertia matrix is not positive definite",
                )
            })?;
        let qdd = factorisation.solve(&rhs);
        match qdd.iter().position(|value| !value.is_finite()) {
            Some(body_index) => Err(Error::forward_dynamics_failure(
                body_index,
                "joint acceleration is not finite",
            )),
            None => Ok(qdd.iter().copied().collect()),
        }
    }
}

/// Joint-space inertia matrix `H_ij = sum_{k >= max(i,j)} I_k` of a coaxial
/// rotor chain, row-major `n x n`. Missing inertias count as zero.
pub(super) fn coaxial_inertia_matrix(n: usize, inertias: &[f64]) -> Vec<f64> {
    let mut h = vec![0.0_f64; n * n];
    for i in 0..n {
        for j in 0..n {
            h[i * n + j] = (i.max(j)..n).map(|k| coordinate(inertias, k)).sum();
        }
    }
    h
}

/// Kinetic energy `0.5 q'^T H q'` of a coaxial rotor chain (see
/// [`coaxial_inertia_matrix`]).
pub(super) fn coaxial_kinetic_energy(n: usize, inertias: &[f64], qd: &[f64]) -> f64 {
    let h = coaxial_inertia_matrix(n, inertias);
    let mut energy = 0.0_f64;
    for i in 0..n {
        for j in 0..n {
            energy += h[i * n + j] * coordinate(qd, i) * coordinate(qd, j);
        }
    }
    0.5 * energy
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::articulated::{
        aba_revolute_chain, coriolis_bias_forces, end_effector_jacobian_2d, forward_dynamics_full,
        gravity_compensation_torques, joint_space_inertia_matrix_planar, kinetic_energy_planar,
        potential_energy_chain, rnea_3d_revolute_z, rnea_revolute_chain,
    };
    use std::f64::consts::{FRAC_PI_2, PI};

    const G: f64 = 9.81;

    /// Small deterministic 64-bit LCG for reproducible test data.
    struct Lcg(u64);

    impl Lcg {
        fn new(seed: u64) -> Self {
            Self(seed)
        }

        fn uniform(&mut self) -> f64 {
            self.0 = self
                .0
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            (self.0 >> 11) as f64 / (1u64 << 53) as f64
        }

        fn range(&mut self, lo: f64, hi: f64) -> f64 {
            lo + (hi - lo) * self.uniform()
        }

        fn vec(&mut self, n: usize, lo: f64, hi: f64) -> Vec<f64> {
            (0..n).map(|_| self.range(lo, hi)).collect()
        }
    }

    /// A random chain together with a random state, accelerations and torques.
    struct Case {
        n: usize,
        masses: Vec<f64>,
        lengths: Vec<f64>,
        inertias: Vec<f64>,
        q: Vec<f64>,
        qd: Vec<f64>,
        qdd: Vec<f64>,
        tau: Vec<f64>,
        gravity: [f64; 3],
    }

    impl Case {
        fn new(n: usize, seed: u64) -> Self {
            let mut rng = Lcg::new(seed ^ ((n as u64) << 32));
            let masses = rng.vec(n, 0.4, 3.0);
            let lengths = rng.vec(n, 0.3, 1.5);
            let inertias: Vec<f64> = masses
                .iter()
                .zip(&lengths)
                .map(|(m, l)| m * l * l / 12.0 * rng.range(0.5, 1.5))
                .collect();
            Self {
                n,
                masses,
                lengths,
                inertias,
                q: rng.vec(n, -PI, PI),
                qd: rng.vec(n, -2.0, 2.0),
                qdd: rng.vec(n, -3.0, 3.0),
                tau: rng.vec(n, -5.0, 5.0),
                gravity: [
                    rng.range(-2.0, 2.0),
                    -rng.range(5.0, 12.0),
                    rng.range(-1.0, 1.0),
                ],
            }
        }

        fn chain(&self) -> PlanarChain {
            PlanarChain::new(self.n, &self.masses, &self.lengths, &self.inertias)
        }
    }

    fn assert_close(actual: &[f64], expected: &[f64], tol: f64, what: &str) {
        assert_eq!(actual.len(), expected.len(), "{what}: length mismatch");
        for (k, (a, e)) in actual.iter().zip(expected).enumerate() {
            assert!(
                (a - e).abs() <= tol * (1.0 + e.abs()),
                "{what}[{k}]: {a} vs {e}"
            );
        }
    }

    fn mat_vec(h: &[f64], x: &[f64]) -> Vec<f64> {
        let n = x.len();
        (0..n)
            .map(|i| (0..n).map(|j| h[i * n + j] * x[j]).sum())
            .collect()
    }

    /// Textbook planar two-link arm (relative joint angles, centres of mass at
    /// half the link lengths, y up, gravity `(0, -g)`): an independent oracle.
    struct TwoLinkArm {
        m1: f64,
        m2: f64,
        l1: f64,
        l2: f64,
        i1: f64,
        i2: f64,
    }

    impl TwoLinkArm {
        fn inertia(&self, q2: f64) -> [[f64; 2]; 2] {
            let (lc1, lc2) = (0.5 * self.l1, 0.5 * self.l2);
            let h11 = self.m1 * lc1 * lc1
                + self.i1
                + self.m2 * (self.l1 * self.l1 + lc2 * lc2 + 2.0 * self.l1 * lc2 * q2.cos())
                + self.i2;
            let h12 = self.m2 * (lc2 * lc2 + self.l1 * lc2 * q2.cos()) + self.i2;
            let h22 = self.m2 * lc2 * lc2 + self.i2;
            [[h11, h12], [h12, h22]]
        }

        fn bias(&self, q2: f64, qd: [f64; 2]) -> [f64; 2] {
            let h = -self.m2 * self.l1 * 0.5 * self.l2 * q2.sin();
            [h * qd[1] * (2.0 * qd[0] + qd[1]), -h * qd[0] * qd[0]]
        }

        fn gravity_torque(&self, q: [f64; 2], g: f64) -> [f64; 2] {
            let (lc1, lc2) = (0.5 * self.l1, 0.5 * self.l2);
            let c12 = (q[0] + q[1]).cos();
            [
                (self.m1 * lc1 + self.m2 * self.l1) * g * q[0].cos() + self.m2 * lc2 * g * c12,
                self.m2 * lc2 * g * c12,
            ]
        }

        fn forward(&self, q: [f64; 2], qd: [f64; 2], tau: [f64; 2], g: f64) -> [f64; 2] {
            let h = self.inertia(q[1]);
            let c = self.bias(q[1], qd);
            let gv = self.gravity_torque(q, g);
            let rhs = [tau[0] - c[0] - gv[0], tau[1] - c[1] - gv[1]];
            let det = h[0][0] * h[1][1] - h[0][1] * h[1][0];
            [
                (h[1][1] * rhs[0] - h[0][1] * rhs[1]) / det,
                (-h[1][0] * rhs[0] + h[0][0] * rhs[1]) / det,
            ]
        }
    }

    fn two_link_case() -> (TwoLinkArm, Vec<f64>, Vec<f64>, Vec<f64>) {
        let arm = TwoLinkArm {
            m1: 2.0,
            m2: 1.5,
            l1: 1.2,
            l2: 0.8,
            i1: 2.0 * 1.2 * 1.2 / 12.0,
            i2: 1.5 * 0.8 * 0.8 / 12.0,
        };
        let masses = vec![arm.m1, arm.m2];
        let lengths = vec![arm.l1, arm.l2];
        let inertias = vec![arm.i1, arm.i2];
        (arm, masses, lengths, inertias)
    }

    /// Gravitational potential energy `-sum m_j g . c_j`.
    fn potential_energy(chain: &PlanarChain, q: &[f64], gravity: [f64; 3]) -> f64 {
        let pose = chain.pose(q);
        chain
            .mass
            .iter()
            .zip(&pose.com)
            .map(|(m, c)| -m * (gravity[0] * c[0] + gravity[1] * c[1]))
            .sum()
    }

    // ── Statics ────────────────────────────────────────────────────────────

    #[test]
    fn static_horizontal_two_link_torques_include_the_distal_lever_arm() {
        let (m, l) = ([2.0, 1.5], [1.2, 0.8]);
        let gravity = [0.0, -G, 0.0];
        // tau_0 holds both links (including m_1 g l_0, the weight of the distal
        // link acting at the full length of link 0); tau_1 holds link 1 alone.
        let expected = [
            G * (m[0] * l[0] / 2.0 + m[1] * (l[0] + l[1] / 2.0)),
            G * m[1] * l[1] / 2.0,
        ];
        let inertias = [m[0] * l[0] * l[0] / 12.0, m[1] * l[1] * l[1] / 12.0];
        let q = [0.0, 0.0];
        let zeros = [0.0, 0.0];
        let from_3d = rnea_3d_revolute_z(2, &m, &inertias, &l, &q, &zeros, &zeros, gravity);
        let from_rods = rnea_revolute_chain(2, &m, &l, &q, &zeros, &zeros, gravity);
        let from_gravity = gravity_compensation_torques(2, &m, &l, &q, gravity);
        assert_close(&from_3d, &expected, 1e-12, "rnea_3d_revolute_z");
        assert_close(&from_rods, &expected, 1e-12, "rnea_revolute_chain");
        assert_close(
            &from_gravity,
            &expected,
            1e-12,
            "gravity_compensation_torques",
        );
    }

    #[test]
    fn inverse_dynamics_at_rest_equals_the_gravity_torques_of_bent_chains() {
        for n in 1..=6 {
            for seed in 0..4 {
                let c = Case::new(n, seed);
                let at_rest = rnea_3d_revolute_z(
                    n,
                    &c.masses,
                    &c.inertias,
                    &c.lengths,
                    &c.q,
                    &[],
                    &[],
                    c.gravity,
                );
                let held = gravity_compensation_torques(n, &c.masses, &c.lengths, &c.q, c.gravity);
                assert_close(&at_rest, &held, 1e-10, "rnea at rest vs gravity torques");
            }
        }
    }

    #[test]
    fn gravity_torques_are_the_gradient_of_the_potential_energy() {
        // Independent check of the lever arms: tau_k = dV/dq_k, with V built
        // directly from the centre-of-mass positions.
        for n in 1..=5 {
            let c = Case::new(n, 7);
            let chain = c.chain();
            let eps = 1e-6;
            let analytic = chain.gravity_torques(&c.q, c.gravity);
            for k in 0..n {
                let mut plus = c.q.clone();
                let mut minus = c.q.clone();
                plus[k] += eps;
                minus[k] -= eps;
                let numeric = (potential_energy(&chain, &plus, c.gravity)
                    - potential_energy(&chain, &minus, c.gravity))
                    / (2.0 * eps);
                assert!(
                    (numeric - analytic[k]).abs() < 1e-6 * (1.0 + analytic[k].abs()),
                    "n={n} k={k}: dV/dq = {numeric}, torque = {}",
                    analytic[k]
                );
            }
        }
    }

    #[test]
    fn bent_chain_gravity_torque_uses_the_true_lever_arm() {
        // Link 0 horizontal, link 1 vertical (up): the weight of link 1 acts at
        // horizontal distance l_0 from joint 0 and at zero distance from joint 1.
        let (m, l) = ([1.0, 1.0], [1.0, 1.0]);
        let gravity = [0.0, -G, 0.0];
        let tau = gravity_compensation_torques(2, &m, &l, &[0.0, FRAC_PI_2], gravity);
        assert!((tau[0] - 1.5 * G).abs() < 1e-12, "tau_0 = {}", tau[0]);
        assert!(tau[1].abs() < 1e-12, "tau_1 = {}", tau[1]);
    }

    // ── Joint-space inertia matrix ────────────────────────────────────────

    #[test]
    fn mass_matrix_of_a_single_link_is_the_inertia_about_the_joint() {
        let h = joint_space_inertia_matrix_planar(1, &[2.0], &[1.5], &[0.3], &[0.9]);
        assert_eq!(h.len(), 1);
        assert!((h[0] - (0.3 + 2.0 * 1.5 * 1.5 / 4.0)).abs() < 1e-12);
    }

    #[test]
    fn mass_matrix_of_a_two_link_arm_matches_the_closed_form() {
        let (arm, m, l, i) = two_link_case();
        for q2 in [-2.5, -0.7, 0.0, 0.4, 1.9, 3.0] {
            let h = joint_space_inertia_matrix_planar(2, &m, &l, &i, &[0.6, q2]);
            let expected = arm.inertia(q2);
            assert_close(
                &h,
                &[
                    expected[0][0],
                    expected[0][1],
                    expected[1][0],
                    expected[1][1],
                ],
                1e-12,
                "H",
            );
        }
    }

    #[test]
    fn mass_matrix_depends_on_the_configuration() {
        // The coupling H_01 contains m_1 l_0 l_c1 cos q_1; the old suffix-sum
        // matrix was constant.
        let (_, m, l, i) = two_link_case();
        let straight = joint_space_inertia_matrix_planar(2, &m, &l, &i, &[0.0, 0.0]);
        let bent = joint_space_inertia_matrix_planar(2, &m, &l, &i, &[0.0, FRAC_PI_2]);
        let expected_difference = m[1] * l[0] * 0.5 * l[1];
        assert!(
            ((straight[1] - bent[1]) - expected_difference).abs() < 1e-12,
            "H_01 straight {} bent {}",
            straight[1],
            bent[1]
        );
        assert!(straight[0] > bent[0], "reach changes H_00 as well");
        assert!((straight[3] - bent[3]).abs() < 1e-14, "H_11 is constant");
    }

    #[test]
    fn mass_matrix_is_symmetric_positive_definite() {
        for n in 1..=8 {
            let c = Case::new(n, 3);
            let h = joint_space_inertia_matrix_planar(n, &c.masses, &c.lengths, &c.inertias, &c.q);
            for i in 0..n {
                for j in 0..n {
                    assert_eq!(h[i * n + j], h[j * n + i], "H not symmetric at ({i},{j})");
                }
            }
            assert!(
                DMatrix::from_row_slice(n, n, &h).cholesky().is_some(),
                "H must be positive definite for n = {n}"
            );
        }
    }

    #[test]
    fn mass_matrix_columns_are_the_inverse_dynamics_of_unit_accelerations() {
        for n in 1..=5 {
            let c = Case::new(n, 11);
            let chain = c.chain();
            let h = chain.mass_matrix(&c.q);
            for k in 0..n {
                let mut unit = vec![0.0; n];
                unit[k] = 1.0;
                let column = chain.inverse_dynamics(&c.q, &[], &unit, [0.0; 3]);
                let expected: Vec<f64> = (0..n).map(|i| h[i * n + k]).collect();
                assert_close(&column, &expected, 1e-10, "H column");
            }
        }
    }

    #[test]
    fn mass_matrix_is_invariant_under_rotation_of_the_whole_chain() {
        let c = Case::new(4, 5);
        let chain = c.chain();
        let mut rotated = c.q.clone();
        rotated[0] += 0.83;
        assert_close(
            &chain.mass_matrix(&rotated),
            &chain.mass_matrix(&c.q),
            1e-12,
            "H",
        );
    }

    // ── Kinetic energy ────────────────────────────────────────────────────

    #[test]
    fn kinetic_energy_is_half_qd_h_qd() {
        for n in 1..=6 {
            let c = Case::new(n, 13);
            let h = joint_space_inertia_matrix_planar(n, &c.masses, &c.lengths, &c.inertias, &c.q);
            let hq = mat_vec(&h, &c.qd);
            let expected = 0.5 * c.qd.iter().zip(&hq).map(|(a, b)| a * b).sum::<f64>();
            let energy = kinetic_energy_planar(n, &c.masses, &c.lengths, &c.inertias, &c.q, &c.qd);
            assert!(
                (energy - expected).abs() < 1e-10 * (1.0 + expected),
                "n={n}: T = {energy}, 0.5 qd.H.qd = {expected}"
            );
        }
    }

    #[test]
    fn kinetic_energy_polarisation_reproduces_the_mass_matrix() {
        // T is quadratic in qd, so H_kl = T(e_k + e_l) - T(e_k) - T(e_l): this
        // checks H against the velocity-based energy, which shares no code with it.
        let c = Case::new(4, 17);
        let h = joint_space_inertia_matrix_planar(4, &c.masses, &c.lengths, &c.inertias, &c.q);
        let energy =
            |qd: &[f64]| kinetic_energy_planar(4, &c.masses, &c.lengths, &c.inertias, &c.q, qd);
        for k in 0..4 {
            let mut ek = [0.0; 4];
            ek[k] = 1.0;
            assert!((2.0 * energy(&ek) - h[k * 4 + k]).abs() < 1e-12);
            for l in (k + 1)..4 {
                let mut el = [0.0; 4];
                el[l] = 1.0;
                let mut both = ek;
                both[l] = 1.0;
                let coupling = energy(&both) - energy(&ek) - energy(&el);
                assert!(
                    (coupling - h[k * 4 + l]).abs() < 1e-12,
                    "H[{k}][{l}] = {} vs polarisation {coupling}",
                    h[k * 4 + l]
                );
            }
        }
    }

    #[test]
    fn kinetic_energy_depends_on_the_configuration() {
        // Both links spin at the same relative rate, once straight and once bent.
        let (_, m, l, i) = two_link_case();
        let qd = [1.0, 1.5];
        let straight = kinetic_energy_planar(2, &m, &l, &i, &[0.0, 0.0], &qd);
        let bent = kinetic_energy_planar(2, &m, &l, &i, &[0.0, 2.0], &qd);
        assert!((straight - bent).abs() > 1e-3, "{straight} vs {bent}");
    }

    // ── Coriolis and centrifugal forces ───────────────────────────────────

    #[test]
    fn single_spinning_link_has_no_coriolis_force() {
        for (q, qd) in [(0.0, 3.0), (0.7, -2.0), (2.2, 5.0)] {
            let bias = coriolis_bias_forces(1, &[2.0], &[1.3], &[q], &[qd]);
            assert!(bias[0].abs() < 1e-12, "C = {} at q = {q}", bias[0]);
        }
    }

    #[test]
    fn coriolis_of_a_two_link_arm_matches_the_closed_form() {
        let (arm, m, l, _) = two_link_case();
        for (q2, qd) in [(0.5, [1.0, -0.5]), (-1.3, [2.0, 0.7]), (2.4, [-1.5, 1.1])] {
            let bias = coriolis_bias_forces(2, &m, &l, &[0.4, q2], &qd);
            let expected = arm.bias(q2, qd);
            assert_close(&bias, &expected, 1e-12, "C");
        }
    }

    #[test]
    fn coriolis_does_not_depend_on_the_rotational_inertias() {
        let c = Case::new(4, 19);
        let light = PlanarChain::new(4, &c.masses, &c.lengths, &[0.0; 4]);
        let heavy = PlanarChain::new(4, &c.masses, &c.lengths, &[5.0, 6.0, 7.0, 8.0]);
        assert_close(
            &light.bias_forces(&c.q, &c.qd),
            &heavy.bias_forces(&c.q, &c.qd),
            1e-12,
            "C",
        );
    }

    #[test]
    fn coriolis_satisfies_the_power_balance_with_the_mass_matrix_rate() {
        // For an unforced, gravity-free system dT/dt = 0, i.e.
        // qd . C(q, qd) = 0.5 qd^T (dH/dt) qd with dH/dt = dH/dq . qd.
        for n in 2..=5 {
            let c = Case::new(n, 23);
            let chain = c.chain();
            let eps = 1e-6;
            let shift = |sign: f64| -> Vec<f64> {
                c.q.iter()
                    .zip(&c.qd)
                    .map(|(q, v)| q + sign * eps * v)
                    .collect()
            };
            let h_plus = chain.mass_matrix(&shift(1.0));
            let h_minus = chain.mass_matrix(&shift(-1.0));
            let h_dot: Vec<f64> = h_plus
                .iter()
                .zip(&h_minus)
                .map(|(a, b)| (a - b) / (2.0 * eps))
                .collect();
            let h_dot_qd = mat_vec(&h_dot, &c.qd);
            let half_qd_hdot_qd = 0.5 * c.qd.iter().zip(&h_dot_qd).map(|(a, b)| a * b).sum::<f64>();
            let bias = chain.bias_forces(&c.q, &c.qd);
            let power: f64 = c.qd.iter().zip(&bias).map(|(a, b)| a * b).sum();
            assert!(
                (power - half_qd_hdot_qd).abs() < 1e-6 * (1.0 + power.abs()),
                "n={n}: qd.C = {power}, 0.5 qd.Hdot.qd = {half_qd_hdot_qd}"
            );
        }
    }

    // ── Inverse dynamics ──────────────────────────────────────────────────

    #[test]
    fn inverse_dynamics_matches_h_qdd_plus_c_plus_g() {
        for n in 1..=6 {
            for seed in 0..3 {
                let c = Case::new(n, seed + 100);
                let chain = c.chain();
                let h = chain.mass_matrix(&c.q);
                let mut expected = mat_vec(&h, &c.qdd);
                let bias = chain.bias_forces(&c.q, &c.qd);
                let gravity = chain.gravity_torques(&c.q, c.gravity);
                for k in 0..n {
                    expected[k] += bias[k] + gravity[k];
                }
                let torques = rnea_3d_revolute_z(
                    n,
                    &c.masses,
                    &c.inertias,
                    &c.lengths,
                    &c.q,
                    &c.qd,
                    &c.qdd,
                    c.gravity,
                );
                assert_close(&torques, &expected, 1e-9, "tau");
            }
        }
    }

    #[test]
    fn inverse_dynamics_is_invariant_under_rotation_of_chain_and_gravity() {
        // Rotating the whole mechanism about the base, together with gravity,
        // leaves the joint torques unchanged.
        let c = Case::new(4, 29);
        let chain = c.chain();
        let phi = 0.9_f64;
        let (s, k) = phi.sin_cos();
        let rotated_gravity = [
            k * c.gravity[0] - s * c.gravity[1],
            s * c.gravity[0] + k * c.gravity[1],
            c.gravity[2],
        ];
        let mut rotated_q = c.q.clone();
        rotated_q[0] += phi;
        assert_close(
            &chain.inverse_dynamics(&rotated_q, &c.qd, &c.qdd, rotated_gravity),
            &chain.inverse_dynamics(&c.q, &c.qd, &c.qdd, c.gravity),
            1e-10,
            "tau",
        );
        assert_close(
            &chain.forward_dynamics(&rotated_q, &c.qd, &c.tau, rotated_gravity),
            &chain.forward_dynamics(&c.q, &c.qd, &c.tau, c.gravity),
            1e-10,
            "q_dd",
        );
    }

    #[test]
    fn rnea_revolute_chain_models_uniform_rods() {
        let c = Case::new(3, 31);
        let rods: Vec<f64> = c
            .masses
            .iter()
            .zip(&c.lengths)
            .map(|(m, l)| m * l * l / 12.0)
            .collect();
        let expected = rnea_3d_revolute_z(
            3, &c.masses, &rods, &c.lengths, &c.q, &c.qd, &c.qdd, c.gravity,
        );
        let torques = rnea_revolute_chain(3, &c.masses, &c.lengths, &c.q, &c.qd, &c.qdd, c.gravity);
        assert_close(&torques, &expected, 1e-12, "tau");
        // The end-pivot inertia m l^2 / 3 must not be used as the COM inertia
        // (the old code double-counted the COM offset).
        let end_pivot: Vec<f64> = c
            .masses
            .iter()
            .zip(&c.lengths)
            .map(|(m, l)| m * l * l / 3.0)
            .collect();
        let wrong = rnea_3d_revolute_z(
            3, &c.masses, &end_pivot, &c.lengths, &c.q, &c.qd, &c.qdd, c.gravity,
        );
        assert!(
            torques
                .iter()
                .zip(&wrong)
                .any(|(a, b)| (a - b).abs() > 1e-3),
            "the rod model must differ from the double-counted inertia"
        );
    }

    #[test]
    fn rnea_revolute_chain_single_link_known_values() {
        // tau = I_joint q_dd + m g (l/2) cos(q) with I_joint = m l^2 / 3.
        let (m, l, q, qdd) = (2.0, 1.5, 0.4_f64, 1.7);
        let tau = rnea_revolute_chain(1, &[m], &[l], &[q], &[0.0], &[qdd], [0.0, -G, 0.0]);
        let expected = m * l * l / 3.0 * qdd + m * G * 0.5 * l * q.cos();
        assert!(
            (tau[0] - expected).abs() < 1e-12,
            "{} vs {expected}",
            tau[0]
        );
    }

    #[test]
    fn inverse_dynamics_of_an_empty_chain_is_empty() {
        assert!(rnea_revolute_chain(0, &[], &[], &[], &[], &[], [0.0, -G, 0.0]).is_empty());
        assert!(rnea_3d_revolute_z(0, &[], &[], &[], &[], &[], &[], [0.0; 3]).is_empty());
        assert!(aba_revolute_chain(0, &[], &[], &[], &[], &[], &[], [0.0; 3]).is_empty());
        assert!(coriolis_bias_forces(0, &[], &[], &[], &[]).is_empty());
        assert!(gravity_compensation_torques(0, &[], &[], &[], [0.0; 3]).is_empty());
        assert!(joint_space_inertia_matrix_planar(0, &[], &[], &[], &[]).is_empty());
        assert_eq!(kinetic_energy_planar(0, &[], &[], &[], &[], &[]), 0.0);
        assert!(
            forward_dynamics_full(0, &[], &[], &[], &[], &[], &[], [0.0; 3])
                .map(|v| v.is_empty())
                .unwrap_or(false)
        );
    }

    /// Evaluates `f(masses, lengths, q, qd)` on slices that are shorter than the
    /// chain (`n = 3`) and on the same data padded with the documented defaults,
    /// and compares the two results.
    fn assert_short_slices_match_padded<F>(name: &str, f: F, tolerance: f64)
    where
        F: Fn(&[f64], &[f64], &[f64], &[f64]) -> Vec<f64>,
    {
        let short = f(&[2.0], &[0.5, 2.0], &[0.3], &[0.5]);
        let padded = f(
            &[2.0, 1.0, 1.0],
            &[0.5, 2.0, 1.0],
            &[0.3, 0.0, 0.0],
            &[0.5, 0.0, 0.0],
        );
        assert_close(&short, &padded, tolerance, name);
    }

    #[test]
    fn short_slices_default_instead_of_panicking() {
        // n = 3 but every slice is shorter: masses/lengths/inertias default to 1,
        // q/qd/qdd/torques to 0, exactly like explicitly padded slices.
        let gravity = [0.0, -G, 0.0];
        assert_short_slices_match_padded(
            "rnea_3d_revolute_z",
            |m, l, q, qd| rnea_3d_revolute_z(3, m, m, l, q, qd, qd, gravity),
            1e-14,
        );
        assert_short_slices_match_padded(
            "rnea_revolute_chain",
            |m, l, q, qd| rnea_revolute_chain(3, m, l, q, qd, qd, gravity),
            1e-14,
        );
        assert_short_slices_match_padded(
            "aba_revolute_chain",
            |m, l, q, qd| aba_revolute_chain(3, m, m, l, q, qd, qd, gravity),
            1e-14,
        );
        assert_short_slices_match_padded(
            "gravity_compensation_torques",
            |m, l, q, _qd| gravity_compensation_torques(3, m, l, q, gravity),
            1e-14,
        );
        assert_short_slices_match_padded(
            "coriolis_bias_forces",
            |m, l, q, qd| coriolis_bias_forces(3, m, l, q, qd),
            1e-14,
        );
        assert_short_slices_match_padded(
            "joint_space_inertia_matrix_planar",
            |m, l, q, _qd| joint_space_inertia_matrix_planar(3, m, l, m, q),
            1e-14,
        );
        assert_short_slices_match_padded(
            "forward_dynamics_full",
            |m, l, q, qd| {
                let qdd = forward_dynamics_full(3, m, m, l, q, qd, qd, gravity)
                    .expect("positive definite inertia matrix");
                assert_eq!(qdd.len(), 3);
                qdd
            },
            1e-12,
        );
        let energy_short = kinetic_energy_planar(3, &[2.0], &[0.5, 2.0], &[2.0], &[0.3], &[0.5]);
        let energy_padded = kinetic_energy_planar(
            3,
            &[2.0, 1.0, 1.0],
            &[0.5, 2.0, 1.0],
            &[2.0, 1.0, 1.0],
            &[0.3, 0.0, 0.0],
            &[0.5, 0.0, 0.0],
        );
        assert!((energy_short - energy_padded).abs() < 1e-14);
    }

    // ── Forward dynamics: articulated body algorithm ──────────────────────

    #[test]
    fn aba_single_link_known_values() {
        // I_joint = I_com + m l^2 / 4 = m l^2 / 3 for a uniform rod.
        let (m, l) = (2.0, 1.5);
        let i_com = m * l * l / 12.0;
        let i_joint = m * l * l / 3.0;
        let no_gravity =
            aba_revolute_chain(1, &[i_com], &[m], &[l], &[0.0], &[0.0], &[3.0], [0.0; 3]);
        assert!(
            (no_gravity[0] - 3.0 / i_joint).abs() < 1e-12,
            "{}",
            no_gravity[0]
        );
        let q = 0.4_f64;
        let with_gravity = aba_revolute_chain(
            1,
            &[i_com],
            &[m],
            &[l],
            &[q],
            &[0.0],
            &[3.0],
            [0.0, -G, 0.0],
        );
        let expected = (3.0 - m * G * 0.5 * l * q.cos()) / i_joint;
        assert!(
            (with_gravity[0] - expected).abs() < 1e-12,
            "{}",
            with_gravity[0]
        );
    }

    #[test]
    fn aba_with_gravity_compensation_torques_holds_the_chain_still() {
        for n in 1..=6 {
            for seed in 0..4 {
                let c = Case::new(n, seed + 200);
                let hold = gravity_compensation_torques(n, &c.masses, &c.lengths, &c.q, c.gravity);
                let qdd = aba_revolute_chain(
                    n,
                    &c.inertias,
                    &c.masses,
                    &c.lengths,
                    &c.q,
                    &[],
                    &hold,
                    c.gravity,
                );
                for (k, a) in qdd.iter().enumerate() {
                    assert!(a.abs() < 1e-9, "n={n} seed={seed}: q_dd[{k}] = {a}");
                }
            }
        }
    }

    #[test]
    fn aba_matches_the_cholesky_forward_dynamics() {
        for n in 1..=7 {
            for seed in 0..4 {
                let c = Case::new(n, seed + 300);
                let by_aba = aba_revolute_chain(
                    n,
                    &c.inertias,
                    &c.masses,
                    &c.lengths,
                    &c.q,
                    &c.qd,
                    &c.tau,
                    c.gravity,
                );
                let by_crba = forward_dynamics_full(
                    n,
                    &c.inertias,
                    &c.masses,
                    &c.lengths,
                    &c.q,
                    &c.qd,
                    &c.tau,
                    c.gravity,
                )
                .expect("random chains have a positive definite inertia matrix");
                assert_close(&by_aba, &by_crba, 1e-9, "q_dd");
            }
        }
    }

    #[test]
    fn aba_inverts_the_inverse_dynamics() {
        for n in 1..=6 {
            let c = Case::new(n, 400);
            let torques = rnea_3d_revolute_z(
                n,
                &c.masses,
                &c.inertias,
                &c.lengths,
                &c.q,
                &c.qd,
                &c.qdd,
                c.gravity,
            );
            let recovered = aba_revolute_chain(
                n,
                &c.inertias,
                &c.masses,
                &c.lengths,
                &c.q,
                &c.qd,
                &torques,
                c.gravity,
            );
            assert_close(&recovered, &c.qdd, 1e-9, "q_dd");
        }
    }

    #[test]
    fn aba_two_link_matches_the_explicit_lagrangian() {
        let (arm, m, l, i) = two_link_case();
        let cases = [
            ([0.3, 0.8], [0.0, 0.0], [0.0, 0.0]),
            ([0.3, 0.8], [1.0, -2.0], [0.0, 0.0]),
            ([-1.1, 2.0], [2.5, 1.5], [3.0, -1.0]),
            ([2.4, -0.6], [-1.0, 3.0], [-2.0, 0.5]),
        ];
        for (q, qd, tau) in cases {
            let expected = arm.forward(q, qd, tau, G);
            let by_aba = aba_revolute_chain(2, &i, &m, &l, &q, &qd, &tau, [0.0, -G, 0.0]);
            assert_close(&by_aba, &expected, 1e-9, "ABA vs Lagrangian");
            let by_crba = forward_dynamics_full(2, &i, &m, &l, &q, &qd, &tau, [0.0, -G, 0.0])
                .expect("positive definite");
            assert_close(&by_crba, &expected, 1e-9, "Cholesky vs Lagrangian");
        }
    }

    #[test]
    fn aba_includes_the_velocity_dependent_terms() {
        // No torque and no gravity: only Coriolis and centrifugal effects move
        // the bent arm, and the base joint pushes the elbow outwards.
        let (arm, m, l, i) = two_link_case();
        let q = [0.0, 1.2];
        let qd = [2.0, 0.0];
        let qdd = aba_revolute_chain(2, &i, &m, &l, &q, &qd, &[0.0, 0.0], [0.0; 3]);
        assert!(qdd.iter().any(|a| a.abs() > 1e-3), "{qdd:?}");
        assert_close(&qdd, &arm.forward(q, qd, [0.0; 2], 0.0), 1e-9, "q_dd");
        // Unlike the old code the answer must change with the velocity.
        let still = aba_revolute_chain(2, &i, &m, &l, &q, &[0.0, 0.0], &[0.0, 0.0], [0.0; 3]);
        assert!(still.iter().all(|a| a.abs() < 1e-12), "{still:?}");
    }

    #[test]
    fn aba_locks_joints_without_inertia() {
        // A massless chain has no articulated inertia: the result must stay
        // finite (zero) instead of dividing by zero.
        let qdd = aba_revolute_chain(
            2,
            &[0.0, 0.0],
            &[0.0, 0.0],
            &[1.0, 1.0],
            &[0.2, 0.4],
            &[1.0, 1.0],
            &[1.0, 1.0],
            [0.0, -G, 0.0],
        );
        assert_eq!(qdd, vec![0.0, 0.0]);
    }

    /// Total energy `T + V` of the state `(q, qd)` under gravity `(0, -g)`.
    fn total_energy(chain: &PlanarChain, q: &[f64], qd: &[f64]) -> f64 {
        chain.kinetic_energy(q, qd) + potential_energy(chain, q, [0.0, -G, 0.0])
    }

    /// Integrates the unforced chain of `case` with the classical Runge-Kutta
    /// method, `steps` steps of size `dt`, using the articulated body algorithm
    /// for the accelerations. Returns the initial total energy and the largest
    /// deviation of the total energy from it along the trajectory.
    fn rk4_energy_drift(case: &Case, steps: usize, dt: f64) -> (f64, f64) {
        fn rhs(chain: &PlanarChain, q: &[f64], qd: &[f64]) -> (Vec<f64>, Vec<f64>) {
            let n = q.len();
            (
                qd.to_vec(),
                chain.forward_dynamics(q, qd, &vec![0.0; n], [0.0, -G, 0.0]),
            )
        }
        fn axpy(x: &[f64], k: &[f64], h: f64) -> Vec<f64> {
            x.iter().zip(k).map(|(a, b)| a + h * b).collect()
        }
        let chain = case.chain();
        let (mut q, mut qd) = (case.q.clone(), case.qd.clone());
        let e0 = total_energy(&chain, &q, &qd);
        let mut max_drift = 0.0_f64;
        for _ in 0..steps {
            let (k1q, k1v) = rhs(&chain, &q, &qd);
            let (k2q, k2v) = rhs(
                &chain,
                &axpy(&q, &k1q, 0.5 * dt),
                &axpy(&qd, &k1v, 0.5 * dt),
            );
            let (k3q, k3v) = rhs(
                &chain,
                &axpy(&q, &k2q, 0.5 * dt),
                &axpy(&qd, &k2v, 0.5 * dt),
            );
            let (k4q, k4v) = rhs(&chain, &axpy(&q, &k3q, dt), &axpy(&qd, &k3v, dt));
            for k in 0..case.n {
                q[k] += dt / 6.0 * (k1q[k] + 2.0 * k2q[k] + 2.0 * k3q[k] + k4q[k]);
                qd[k] += dt / 6.0 * (k1v[k] + 2.0 * k2v[k] + 2.0 * k3v[k] + k4v[k]);
            }
            max_drift = max_drift.max((total_energy(&chain, &q, &qd) - e0).abs());
        }
        (e0, max_drift)
    }

    #[test]
    fn unforced_two_link_chain_conserves_energy_under_the_aba() {
        // One second at dt = 1e-4: kinetic + potential energy must stay constant,
        // which holds only if H, C and G (and the ABA built on them) are
        // consistent.
        let (e0, drift) = rk4_energy_drift(&Case::new(2, 500), 10_000, 1e-4);
        assert!(
            drift < 1e-8 * (1.0 + e0.abs()),
            "energy drifted by {drift} from {e0}"
        );
    }

    #[test]
    fn energy_drift_of_unforced_chains_shrinks_with_the_fourth_power_of_the_step() {
        // For a consistent model the only energy error is the RK4 truncation
        // error, which shrinks about 16-fold (measured 14 to 16) each time the
        // step is halved. An inconsistency between H, C and G would instead leave
        // a drift that does not vanish with the step size. One second of motion
        // of random 2- to 4-link chains with steps of 4e-3 down to 5e-4.
        for n in 2..=4 {
            let case = Case::new(n, 500);
            let mut e0 = 0.0_f64;
            let mut drifts = Vec::new();
            for halvings in 0..4_u32 {
                let steps = 250_usize << halvings;
                let dt = 4e-3 / f64::from(1_u32 << halvings);
                let (initial_energy, drift) = rk4_energy_drift(&case, steps, dt);
                e0 = initial_energy;
                drifts.push(drift);
            }
            for pair in drifts.windows(2) {
                assert!(
                    pair[1] < pair[0] / 8.0,
                    "n={n}: halving the step must cut the drift at least 8-fold, got {drifts:?}"
                );
            }
            assert!(
                drifts[3] < 1e-5 * (1.0 + e0.abs()),
                "n={n}: energy drifted by {} from {e0} at dt = 5e-4",
                drifts[3]
            );
        }
    }

    // ── Forward dynamics: Cholesky solve of H q_dd = tau - C - G ───────────

    #[test]
    fn forward_dynamics_full_two_link_matches_the_explicit_inverse() {
        // No velocity and no gravity: q_dd = H^-1 tau with the 2x2 inverse
        // written out, to 1e-12.
        let (arm, m, l, i) = two_link_case();
        let q = [0.7, 1.1];
        let tau = [2.0, -1.0];
        let h = arm.inertia(q[1]);
        let det = h[0][0] * h[1][1] - h[0][1] * h[1][0];
        let expected = [
            (h[1][1] * tau[0] - h[0][1] * tau[1]) / det,
            (-h[1][0] * tau[0] + h[0][0] * tau[1]) / det,
        ];
        let qdd = forward_dynamics_full(2, &i, &m, &l, &q, &[0.0, 0.0], &tau, [0.0; 3])
            .expect("positive definite");
        assert_close(&qdd, &expected, 1e-12, "q_dd");
    }

    #[test]
    fn forward_dynamics_full_couples_the_joints() {
        // A torque on the elbow alone accelerates the shoulder as well; a solve
        // that divides by diag(H) would leave q_dd[0] at zero.
        let (_, m, l, i) = two_link_case();
        let qdd = forward_dynamics_full(
            2,
            &i,
            &m,
            &l,
            &[0.0, 0.5],
            &[0.0, 0.0],
            &[0.0, 1.0],
            [0.0; 3],
        )
        .expect("positive definite");
        assert!(
            qdd[0] < -1e-3,
            "shoulder must react to the elbow torque: {qdd:?}"
        );
        assert!(qdd[1] > 0.0);
    }

    #[test]
    fn forward_dynamics_full_single_spinning_link_does_not_accelerate() {
        let qdd =
            forward_dynamics_full(1, &[0.2], &[1.0], &[1.0], &[0.3], &[1.0], &[0.0], [0.0; 3])
                .expect("positive definite");
        assert!(qdd[0].abs() < 1e-12, "q_dd = {}", qdd[0]);
    }

    #[test]
    fn forward_dynamics_full_reports_a_singular_inertia_matrix() {
        let err = forward_dynamics_full(
            2,
            &[0.0, 0.0],
            &[0.0, 0.0],
            &[1.0, 1.0],
            &[0.0, 0.0],
            &[0.0, 0.0],
            &[1.0, 1.0],
            [0.0, -G, 0.0],
        )
        .expect_err("a massless chain has a singular inertia matrix");
        assert!(
            matches!(err, Error::ForwardDynamicsFailure { body_index: 0, .. }),
            "{err:?}"
        );
        // A non-finite configuration is reported on the offending joint.
        let err = forward_dynamics_full(
            2,
            &[0.1, 0.1],
            &[1.0, 1.0],
            &[1.0, 1.0],
            &[0.0, f64::NAN],
            &[0.0, 0.0],
            &[0.0, 0.0],
            [0.0, -G, 0.0],
        )
        .expect_err("NaN joint angle");
        assert!(
            matches!(err, Error::ForwardDynamicsFailure { .. }),
            "{err:?}"
        );
    }

    // ── Kinematics shared with the other planar helpers ───────────────────

    #[test]
    fn pose_agrees_with_the_planar_kinematics_helpers() {
        let c = Case::new(4, 41);
        let chain = c.chain();
        let pose = chain.pose(&c.q);
        // Potential energy with unit gravity is the mass-weighted COM height.
        let expected: f64 = c.masses.iter().zip(&pose.com).map(|(m, p)| m * p[1]).sum();
        let helper = potential_energy_chain(4, &c.masses, &c.lengths, &c.q, 1.0, [0.0, 0.0]);
        assert!((helper - expected).abs() < 1e-12, "{helper} vs {expected}");
        // The tip of the last link is where the end-effector Jacobian puts it:
        // its columns are the tip velocity per unit joint rate, so the tip moves
        // by the Jacobian times a small joint step.
        let jac = end_effector_jacobian_2d(4, &c.lengths, &c.q, [0.0, 0.0]);
        let step = [1e-6, -2e-6, 1.5e-6, 0.5e-6];
        let moved: Vec<f64> = c.q.iter().zip(&step).map(|(q, d)| q + d).collect();
        let tip_before = chain.pose(&c.q).origin[4];
        let tip_after = chain.pose(&moved).origin[4];
        for axis in 0..2 {
            let predicted: f64 = jac.iter().zip(&step).map(|(col, d)| col[axis] * d).sum();
            assert!(
                (tip_after[axis] - tip_before[axis] - predicted).abs() < 1e-9,
                "axis {axis}"
            );
        }
    }

    // ── Coaxial rotor chain (deprecated public helpers) ───────────────────

    #[test]
    fn coaxial_inertia_matrix_is_the_suffix_sum_of_the_inertias() {
        let h = coaxial_inertia_matrix(1, &[5.0]);
        assert_eq!(h, vec![5.0]);
        let h = coaxial_inertia_matrix(2, &[3.0, 2.0]);
        assert_eq!(h, vec![5.0, 2.0, 2.0, 2.0]);
        let h = coaxial_inertia_matrix(3, &[1.0, 2.0, 3.0]);
        for i in 0..3 {
            for j in 0..3 {
                assert_eq!(h[i * 3 + j], h[j * 3 + i]);
            }
        }
        assert_eq!(coaxial_inertia_matrix(2, &[3.0]), vec![3.0, 0.0, 0.0, 0.0]);
    }

    #[test]
    fn coaxial_kinetic_energy_is_half_the_sum_of_i_omega_squared() {
        // Absolute rates 1 and 3: T = 0.5 (3 * 1^2 + 2 * 3^2) = 10.5.
        let energy = coaxial_kinetic_energy(2, &[3.0, 2.0], &[1.0, 2.0]);
        assert!((energy - 10.5).abs() < 1e-12, "{energy}");
        assert_eq!(coaxial_kinetic_energy(2, &[1.0, 0.5], &[0.0, 0.0]), 0.0);
    }
}
