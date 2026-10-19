// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Stretch / shear coupling and elastic loads of the [`DiscreteCosseratRod`].
//!
//! The nodes of the discrete rod carry a position `x_i` and a material frame
//! `R_i`.  A segment between two nodes has the Simo–Reissner strain measures
//!
//! ```text
//! φ = log(R_iᵀ R_{i+1}),   κ = φ / ds                 (bending κ₁, κ₂ and twist τ)
//! R_m = R_i · exp(φ / 2)                              (midpoint frame)
//! γ = R_mᵀ (x_{i+1} − x_i) / ds − e₃                  (shear γ₁, γ₂ and stretch γ₃)
//! ```
//!
//! and the stored energy
//!
//! ```text
//! W = Σ ds · [ ½ (GA₁γ₁² + GA₂γ₂² + EAγ₃²) + ½ (EI₁Δκ₁² + EI₂Δκ₂² + GJΔτ²) ].
//! ```
//!
//! The nodal loads are its exact negative gradient (the rotations enter through the
//! left Jacobian of SO(3), see [`left_jacobian`]):
//!
//! - the force `n = R_m · (GA₁γ₁, GA₂γ₂, EAγ₃)` acts as `+n` on node `i` and as
//!   `−n` on node `i + 1`: a stretched segment draws its end nodes together;
//! - the force pair exerts the torque `−Δx × n`, which the frames take up so that
//!   angular momentum is conserved.  The couple `T = Δx × n` turns the midpoint
//!   frame towards the tangent; it is shared between the two end frames by the
//!   exact derivative `W` of the midpoint frame, `T_i = (I − Wᵀ) T`,
//!   `T_{i+1} = Wᵀ T` with `W ≈ ½ I` for small relative rotations.  This is the
//!   `d₃ × n` term of the Cosserat balance law, and it is what couples bending to
//!   the positions of the nodes: a bending moment rotates the frames, the frames
//!   leave the tangent, and the shear force moves the nodes;
//! - the joint moment `M = EI (κ − κ₀)` of the bending / twist strain acts as
//!   `+R_i J⁻ᵀ M` on node `i` and as `−R_i J⁻ᵀ M` on node `i + 1`.
//!
//! Forces and moments together satisfy `Σ f = 0` and `Σ (x × f + m) = 0` exactly,
//! so a free rod conserves its linear and angular momentum.

use super::{
    DirectorTriad, DiscreteCosseratRod, add3, cross3, internal_force, internal_moment,
    kirchhoff_energy_density, mat3_identity, mat3_mul, mat3_mul_vec, mat3_transpose, norm3,
    rodrigues, scale3, simo_reissner_energy_density, stretch_shear_strain, sub3,
};

/// Row-major 3×3 matrix.
type Mat3 = [f64; 9];

/// Below this rotation angle the SO(3) Jacobians use their Taylor series.
const SERIES_ANGLE: f64 = 1.0e-2;

fn mat3_add(a: Mat3, b: Mat3) -> Mat3 {
    let mut out = [0.0_f64; 9];
    for (o, (x, y)) in out.iter_mut().zip(a.iter().zip(b.iter())) {
        *o = x + y;
    }
    out
}

fn mat3_scale(a: Mat3, s: f64) -> Mat3 {
    a.map(|x| x * s)
}

/// Skew-symmetric matrix `[v]×` with `[v]× w = v × w`.
fn skew(v: [f64; 3]) -> Mat3 {
    [0.0, -v[2], v[1], v[2], 0.0, -v[0], -v[1], v[0], 0.0]
}

/// Unit quaternion `[w, x, y, z]` with `w >= 0` of the rotation matrix `m`
/// (row-major), by Shepperd's method (accurate for every rotation angle).
fn quaternion_from_matrix(m: Mat3) -> [f64; 4] {
    let trace = m[0] + m[4] + m[8];
    let mut q = if trace > 0.0 {
        let s = 2.0 * (1.0 + trace).sqrt();
        [
            0.25 * s,
            (m[7] - m[5]) / s,
            (m[2] - m[6]) / s,
            (m[3] - m[1]) / s,
        ]
    } else if m[0] > m[4] && m[0] > m[8] {
        let s = 2.0 * (1.0 + m[0] - m[4] - m[8]).sqrt();
        [
            (m[7] - m[5]) / s,
            0.25 * s,
            (m[1] + m[3]) / s,
            (m[2] + m[6]) / s,
        ]
    } else if m[4] > m[8] {
        let s = 2.0 * (1.0 + m[4] - m[0] - m[8]).sqrt();
        [
            (m[2] - m[6]) / s,
            (m[1] + m[3]) / s,
            0.25 * s,
            (m[5] + m[7]) / s,
        ]
    } else {
        let s = 2.0 * (1.0 + m[8] - m[0] - m[4]).sqrt();
        [
            (m[3] - m[1]) / s,
            (m[2] + m[6]) / s,
            (m[5] + m[7]) / s,
            0.25 * s,
        ]
    };
    // q and -q describe the same rotation: take the shorter one.
    if q[0] < 0.0 {
        for c in &mut q {
            *c = -*c;
        }
    }
    let norm = q.iter().map(|c| c * c).sum::<f64>().sqrt();
    if norm > 0.0 {
        for c in &mut q {
            *c /= norm;
        }
    }
    q
}

/// Rotation vector `φ = θ u` (`θ ∈ [0, π]`) of the rotation matrix `m`: the
/// logarithm map of SO(3), accurate for every angle (it goes through the unit
/// quaternion rather than `acos((tr R − 1) / 2)`).
pub(super) fn rotation_vector(m: Mat3) -> [f64; 3] {
    let [w, x, y, z] = quaternion_from_matrix(m);
    let sin_half = (x * x + y * y + z * z).sqrt();
    if sin_half < 1.0e-9 {
        // log(q) = 2 atan2(|v|, w) v / |v| -> 2 v / w for tiny |v|.
        let factor = 2.0 / w.max(f64::MIN_POSITIVE);
        return [factor * x, factor * y, factor * z];
    }
    let factor = 2.0 * sin_half.atan2(w) / sin_half;
    [factor * x, factor * y, factor * z]
}

/// Rotation vector of the relative rotation `R_aᵀ R_b` of two frames, expressed in
/// the frame of `a`.
pub(super) fn relative_rotation_vector(ra: Mat3, rb: Mat3) -> [f64; 3] {
    rotation_vector(mat3_mul(mat3_transpose(ra), rb))
}

/// Rotation matrix `exp([φ]×)` of the rotation vector `φ`.
pub(super) fn exp_rotation(phi: [f64; 3]) -> Mat3 {
    let angle = norm3(phi);
    if angle < 1.0e-300 {
        mat3_identity()
    } else {
        rodrigues(phi, angle)
    }
}

/// Left Jacobian of SO(3),
/// `J(φ) = I + (1 − cos θ) / θ² [φ]× + (θ − sin θ) / θ³ [φ]×²`, which relates a
/// perturbation of the rotation vector to a left perturbation of the rotation,
/// `exp(φ + δ) = exp(J(φ) δ) exp(φ)`.
pub(super) fn left_jacobian(phi: [f64; 3]) -> Mat3 {
    let theta = norm3(phi);
    let (a, b) = if theta < SERIES_ANGLE {
        let t2 = theta * theta;
        (
            0.5 - t2 / 24.0 + t2 * t2 / 720.0,
            1.0 / 6.0 - t2 / 120.0 + t2 * t2 / 5040.0,
        )
    } else {
        (
            (1.0 - theta.cos()) / (theta * theta),
            (theta - theta.sin()) / (theta * theta * theta),
        )
    };
    let k = skew(phi);
    mat3_add(
        mat3_add(mat3_identity(), mat3_scale(k, a)),
        mat3_scale(mat3_mul(k, k), b),
    )
}

/// Inverse of the left Jacobian,
/// `J⁻¹(φ) = I − ½ [φ]× + (1 / θ² − cot(θ / 2) / (2θ)) [φ]×²`, which maps a left
/// perturbation of the rotation to the change of its rotation vector,
/// `log(exp(δ) exp(φ)) = φ + J⁻¹(φ) δ`.  It is finite up to `θ = π`.
pub(super) fn left_jacobian_inverse(phi: [f64; 3]) -> Mat3 {
    let theta = norm3(phi);
    let c = if theta < SERIES_ANGLE {
        let t2 = theta * theta;
        1.0 / 12.0 + t2 / 720.0 + t2 * t2 / 30240.0
    } else {
        1.0 / (theta * theta) - 0.5 / (theta * (0.5 * theta).tan())
    };
    let k = skew(phi);
    mat3_add(
        mat3_add(mat3_identity(), mat3_scale(k, -0.5)),
        mat3_scale(mat3_mul(k, k), c),
    )
}

/// Geodesic midpoint `R_m = R_a · exp(½ log(R_aᵀ R_b))` of two rotation matrices:
/// the rotation halfway along the shortest rotation from `R_a` to `R_b`.
pub fn midpoint_rotation(ra: [f64; 9], rb: [f64; 9]) -> [f64; 9] {
    let half = scale3(relative_rotation_vector(ra, rb), 0.5);
    mat3_mul(ra, exp_rotation(half))
}

/// Geodesic midpoint of two material frames (see [`midpoint_rotation`]).
pub fn midpoint_frame(a: &DirectorTriad, b: &DirectorTriad) -> DirectorTriad {
    DirectorTriad::from_rotation_matrix(midpoint_rotation(
        a.to_rotation_matrix(),
        b.to_rotation_matrix(),
    ))
}

/// Everything one segment contributes to the nodal loads.
struct SegmentLoads {
    /// Internal force `n` in the global frame (`+n` on node `i`, `−n` on node `i + 1`).
    force: [f64; 3],
    /// Joint moment `R_i J⁻ᵀ M` (`+` on node `i`, `−` on node `i + 1`).
    joint_moment: [f64; 3],
    /// Couple of the force pair taken up by the frame of node `i`.
    couple_first: [f64; 3],
    /// Couple of the force pair taken up by the frame of node `i + 1`.
    couple_second: [f64; 3],
}

impl DiscreteCosseratRod {
    /// Material frame of segment `i`: the geodesic midpoint of the frames of its
    /// two end nodes.
    pub fn segment_frame(&self, i: usize) -> DirectorTriad {
        midpoint_frame(&self.nodes[i].frame, &self.nodes[i + 1].frame)
    }

    /// Stretch-shear strain `γ = R_mᵀ (x_{i+1} − x_i) / ds − e₃` of segment `i` in
    /// the midpoint frame: `γ₁`, `γ₂` shear along `d₁`, `d₂` and `γ₃` stretch along
    /// `d₃`.
    pub fn stretch_shear(&self, i: usize) -> [f64; 3] {
        let ds = self.rest_lengths[i];
        let dr = sub3(self.nodes[i + 1].position, self.nodes[i].position);
        stretch_shear_strain(&self.segment_frame(i), scale3(dr, 1.0 / ds))
    }

    /// Internal force `n = R_m · (GA₁γ₁, GA₂γ₂, EAγ₃)` of segment `i` in the global
    /// frame: the pull the segment exerts on its first node (`+n`; the second node
    /// feels `−n`).
    pub fn segment_force(&self, i: usize) -> [f64; 3] {
        self.segment_loads(i).force
    }

    fn segment_loads(&self, i: usize) -> SegmentLoads {
        let ds = self.rest_lengths[i];
        let ra = self.nodes[i].frame.to_rotation_matrix();
        let rb = self.nodes[i + 1].frame.to_rotation_matrix();
        let phi = relative_rotation_vector(ra, rb);
        let rm = mat3_mul(ra, exp_rotation(scale3(phi, 0.5)));

        // Stretch / shear force in the global frame.
        let dr = sub3(self.nodes[i + 1].position, self.nodes[i].position);
        let gamma = stretch_shear_strain(
            &DirectorTriad::from_rotation_matrix(rm),
            scale3(dr, 1.0 / ds),
        );
        let force = mat3_mul_vec(rm, internal_force(&self.props, gamma));

        // Bending / twist joint moment: dW = M . dφ and dφ = J⁻¹ R_aᵀ (δθ_b − δθ_a).
        let kappa = scale3(phi, 1.0 / ds);
        let moment = internal_moment(&self.props, kappa, self.intrinsic_kappa[i]);
        let j_inv = left_jacobian_inverse(phi);
        let joint_moment = mat3_mul_vec(ra, mat3_mul_vec(mat3_transpose(j_inv), moment));

        // The stretch / shear energy depends on the rotations through R_m only,
        // dW = -(Δx × n) . δθ_m with δθ_m = (I - W) δθ_a + W δθ_b and
        // W = ½ R_a J(φ/2) J⁻¹(φ) R_aᵀ  (W = ½ I for a straight rod).
        let w = mat3_scale(
            mat3_mul(
                mat3_mul(ra, mat3_mul(left_jacobian(scale3(phi, 0.5)), j_inv)),
                mat3_transpose(ra),
            ),
            0.5,
        );
        let couple = cross3(dr, force);
        let couple_second = mat3_mul_vec(mat3_transpose(w), couple);
        let couple_first = sub3(couple, couple_second);

        SegmentLoads {
            force,
            joint_moment,
            couple_first,
            couple_second,
        }
    }

    /// Elastic stretch and shear forces on every node (Simo–Reissner force law):
    /// `+n` on the first and `−n` on the second node of each segment.  A stretched
    /// segment draws its end nodes together, and a segment whose tangent leaves the
    /// directors pushes its nodes back in line with them.
    pub fn elastic_forces(&self) -> Vec<[f64; 3]> {
        let mut forces = vec![[0.0_f64; 3]; self.nodes.len()];
        for i in 0..self.num_segments() {
            let n = self.segment_loads(i).force;
            forces[i] = add3(forces[i], n);
            forces[i + 1] = sub3(forces[i + 1], n);
        }
        forces
    }

    /// Elastic restoring moments the joints exert on every node, in the global
    /// frame.
    ///
    /// The Darboux vector `κ = log(R_iᵀ R_{i+1}) / ds` of segment `i` is the
    /// rotation that takes frame `i` to frame `i + 1`, expressed in frame `i`.  The
    /// joint moment `R_i J⁻ᵀ EI (κ − κ₀)` therefore turns node `i` *towards* node
    /// `i + 1` (`+M`) and node `i + 1` back towards node `i` (`−M`): action and
    /// reaction, so the moments of all nodes sum to zero and a curvature mismatch
    /// relaxes instead of producing a net torque.  (`J⁻ᵀ ≈ I + ½ [φ]×` is the
    /// correction of the exact energy gradient for a finite relative rotation `φ`.)
    pub fn elastic_moments(&self) -> Vec<[f64; 3]> {
        let mut moments = vec![[0.0_f64; 3]; self.nodes.len()];
        for i in 0..self.num_segments() {
            let moment = self.segment_loads(i).joint_moment;
            moments[i] = add3(moments[i], moment);
            moments[i + 1] = sub3(moments[i + 1], moment);
        }
        moments
    }

    /// Couples of the stretch / shear force pairs on the node frames.
    ///
    /// The pair `+n` at `x_i`, `−n` at `x_{i+1}` exerts the torque `−Δx × n` about
    /// the origin, so the frames take up `+Δx × n` (half on each node for a small
    /// relative rotation) to conserve angular momentum; this is the `d₃ × n` term
    /// of the Cosserat balance law and the only way the shear force can turn the
    /// frames back towards the tangent.
    pub fn shear_couple_moments(&self) -> Vec<[f64; 3]> {
        let mut moments = vec![[0.0_f64; 3]; self.nodes.len()];
        for i in 0..self.num_segments() {
            let loads = self.segment_loads(i);
            moments[i] = add3(moments[i], loads.couple_first);
            moments[i + 1] = add3(moments[i + 1], loads.couple_second);
        }
        moments
    }

    /// All elastic loads on the nodes in the global frame, as `(forces, moments)`:
    /// the stretch / shear forces ([`elastic_forces`](Self::elastic_forces)) and the
    /// sum of the bending / twist joint moments
    /// ([`elastic_moments`](Self::elastic_moments)) and the couples of the force
    /// pairs ([`shear_couple_moments`](Self::shear_couple_moments)).
    ///
    /// They are the exact negative gradient of [`elastic_energy`](Self::elastic_energy)
    /// with respect to the node positions and to (global) rotations of the node
    /// frames.  The loads of a free rod add up to zero force and, about any point,
    /// zero torque: `Σ f = 0` and `Σ (x × f + m) = 0`.
    pub fn elastic_loads(&self) -> (Vec<[f64; 3]>, Vec<[f64; 3]>) {
        let mut forces = vec![[0.0_f64; 3]; self.nodes.len()];
        let mut moments = vec![[0.0_f64; 3]; self.nodes.len()];
        for i in 0..self.num_segments() {
            let loads = self.segment_loads(i);
            forces[i] = add3(forces[i], loads.force);
            forces[i + 1] = sub3(forces[i + 1], loads.force);
            moments[i] = add3(add3(moments[i], loads.joint_moment), loads.couple_first);
            moments[i + 1] = add3(
                sub3(moments[i + 1], loads.joint_moment),
                loads.couple_second,
            );
        }
        (forces, moments)
    }

    /// Stretch / shear energy `Σ ½ ds (GA₁γ₁² + GA₂γ₂² + EAγ₃²)` \[J\].
    pub fn stretch_shear_energy(&self) -> f64 {
        (0..self.num_segments())
            .map(|i| {
                simo_reissner_energy_density(&self.props, self.stretch_shear(i), [0.0; 3])
                    * self.rest_lengths[i]
            })
            .sum()
    }

    /// Bending / twist (Kirchhoff) energy `Σ ½ ds (EI₁Δκ₁² + EI₂Δκ₂² + GJΔτ²)` \[J\].
    pub fn bending_twist_energy(&self) -> f64 {
        (0..self.num_segments())
            .map(|i| {
                kirchhoff_energy_density(&self.props, self.darboux(i), self.intrinsic_kappa[i])
                    * self.rest_lengths[i]
            })
            .sum()
    }
}
