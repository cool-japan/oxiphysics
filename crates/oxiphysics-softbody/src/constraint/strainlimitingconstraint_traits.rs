//! # StrainLimitingConstraint - Trait Implementations
//!
//! This module contains trait implementations for `StrainLimitingConstraint`.
//!
//! ## Implemented Traits
//!
//! - `SoftConstraint`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use crate::particle::SoftParticle;
use oxiphysics_core::math::{Real, Vec3};

use super::functions::SoftConstraint;
use super::types::StrainLimitingConstraint;

/// A 2x2 matrix, row-major: `m[row][col]`.
type Mat2 = [[Real; 2]; 2];

fn mat2_mul(a: Mat2, b: Mat2) -> Mat2 {
    [
        [
            a[0][0] * b[0][0] + a[0][1] * b[1][0],
            a[0][0] * b[0][1] + a[0][1] * b[1][1],
        ],
        [
            a[1][0] * b[0][0] + a[1][1] * b[1][0],
            a[1][0] * b[0][1] + a[1][1] * b[1][1],
        ],
    ]
}

/// Rotation matrix `[[cos a, -sin a], [sin a, cos a]]`.
fn rotation2(angle: Real) -> Mat2 {
    let (s, c) = angle.sin_cos();
    [[c, -s], [s, c]]
}

/// Singular value decomposition `F = U diag(s) V^T` of a 2x2 matrix with `U` and
/// `V^T` rotations and `s[0] >= |s[1]|`; `s[1]` is negative only if `det F < 0`.
struct Svd2 {
    u: Mat2,
    s: [Real; 2],
    vt: Mat2,
}

/// Closed-form 2x2 SVD (Blinn): split `F` into a rotation-scaling part
/// `[[e, -k], [k, e]]` and a reflection-scaling part `[[g, h], [h, -g]]`; then
/// `F = Rot(b) diag(q + r, q - r) Rot(c)` with `q = |(e, k)|`, `r = |(g, h)|`,
/// `b = (atan2(k, e) + atan2(h, g)) / 2` and `c = (atan2(k, e) - atan2(h, g)) / 2`.
fn svd2(f: Mat2) -> Svd2 {
    let e = 0.5 * (f[0][0] + f[1][1]);
    let g = 0.5 * (f[0][0] - f[1][1]);
    let h = 0.5 * (f[0][1] + f[1][0]);
    let k = 0.5 * (f[1][0] - f[0][1]);
    let q = e.hypot(k);
    let r = g.hypot(h);
    let a_conformal = k.atan2(e);
    let a_reflecting = h.atan2(g);
    Svd2 {
        u: rotation2(0.5 * (a_conformal + a_reflecting)),
        s: [q + r, q - r],
        vt: rotation2(0.5 * (a_conformal - a_reflecting)),
    }
}

/// Inverse of a 2x2 matrix, or `None` if it is (numerically) singular.
fn mat2_inverse(m: Mat2) -> Option<Mat2> {
    let det = m[0][0] * m[1][1] - m[0][1] * m[1][0];
    if det.abs() < 1e-14 {
        return None;
    }
    Some([
        [m[1][1] / det, -m[0][1] / det],
        [-m[1][0] / det, m[0][0] / det],
    ])
}

/// Inverse-mass weighted correction that moves the two edges of a triangle by
/// `delta1 = e1* - e1` and `delta2 = e2* - e2`.
///
/// With `w` the inverse masses, the vertex displacements `d0, d1, d2` minimise
/// `sum |d_i|^2 / w_i` subject to `d1 - d0 = delta1` and `d2 - d0 = delta2` (so the
/// mass-weighted centroid stays fixed).  Writing the Lagrange multipliers as
/// `mu1, mu2` gives `d0 = -w0 (mu1 + mu2)`, `d1 = w1 mu1`, `d2 = w2 mu2` and the
/// 2x2 system `(S + alpha I) mu = delta` with
/// `S = [[w0 + w1, w0], [w0, w0 + w2]]`; `alpha` is the XPBD compliance term
/// `compliance / dt^2` (infinite compliance leaves the triangle alone).
///
/// `S` is inverted through its eigen decomposition, dropping zero eigenvalues
/// (pseudo-inverse): if two vertices are pinned the single free vertex takes the
/// least-squares solution on the edges it can still influence, and if all are
/// pinned nothing moves.
fn edge_corrections(w: [Real; 3], alpha: Real, delta1: Vec3, delta2: Vec3) -> [Vec3; 3] {
    let s00 = w[0] + w[1];
    let s01 = w[0];
    let s11 = w[0] + w[2];
    let mean = 0.5 * (s00 + s11);
    let radius = (0.5 * (s00 - s11)).hypot(s01);
    // Principal axis of the symmetric matrix [[s00, s01], [s01, s11]].
    let angle = 0.5 * (2.0 * s01).atan2(s00 - s11);
    let (sin_a, cos_a) = angle.sin_cos();
    let eigen = [
        (mean + radius, [cos_a, sin_a]),
        (mean - radius, [-sin_a, cos_a]),
    ];
    let cutoff = 1e-12 * (mean + radius).max(0.0);
    let mut mu1 = Vec3::zeros();
    let mut mu2 = Vec3::zeros();
    for (lambda, q) in eigen {
        if alpha == 0.0 && lambda <= cutoff {
            continue;
        }
        let denom = lambda + alpha;
        if !denom.is_finite() || denom <= 0.0 {
            continue;
        }
        let projected = (delta1 * q[0] + delta2 * q[1]) / denom;
        mu1 += projected * q[0];
        mu2 += projected * q[1];
    }
    [(mu1 + mu2) * -w[0], mu1 * w[1], mu2 * w[2]]
}

impl SoftConstraint for StrainLimitingConstraint {
    /// Clamp the principal stretches of the triangle into
    /// `[min_compression, max_stretch]` and move the vertices to the clamped shape.
    ///
    /// With `alpha = compliance / dt_sub^2` the displacements are those of the XPBD
    /// system `(S + alpha I) mu = delta` (see `edge_corrections`): exact for
    /// `compliance == 0`, and ever softer as the compliance grows.
    fn project(&mut self, particles: &mut [SoftParticle], dt_sub: Real) {
        let [i0, i1, i2] = self.indices;
        let p0 = particles[i0].position;
        let p1 = particles[i1].position;
        let p2 = particles[i2].position;
        let d1 = p1 - p0;
        let d2 = p2 - p0;
        let e1_len = d1.norm();
        if e1_len < 1e-14 {
            return;
        }
        let u = d1 / e1_len;
        let normal = d1.cross(&d2);
        if normal.norm() < 1e-14 {
            return;
        }
        let v = normal.cross(&d1).normalize();
        // Current edge matrix Ds in the (u, v) frame of the triangle (columns d1, d2)
        // and the deformation gradient F = Ds Dm^-1.
        let ds: Mat2 = [[e1_len, d2.dot(&u)], [0.0, d2.dot(&v)]];
        let inv_rest: Mat2 = [
            [self.inv_rest[0], self.inv_rest[1]],
            [self.inv_rest[2], self.inv_rest[3]],
        ];
        let Some(dm) = mat2_inverse(inv_rest) else {
            return;
        };
        let svd = svd2(mat2_mul(ds, inv_rest));

        // Clamp the singular values into the (ordered) limits.
        let lower = self.min_compression.min(self.max_stretch);
        let upper = self.min_compression.max(self.max_stretch);
        let clamped = [
            svd.s[0].max(lower).min(upper),
            svd.s[1].max(lower).min(upper),
        ];
        let violation = (svd.s[0] - clamped[0]).abs() + (svd.s[1] - clamped[1]).abs();
        if violation < 1e-14 {
            return;
        }

        // Target edges E* = F* Dm with F* = U diag(clamped) V^T, in the (u, v) frame.
        let sigma_vt = mat2_mul([[clamped[0], 0.0], [0.0, clamped[1]]], svd.vt);
        let target = mat2_mul(mat2_mul(svd.u, sigma_vt), dm);
        let delta1 = u * (target[0][0] - ds[0][0]) + v * (target[1][0] - ds[1][0]);
        let delta2 = u * (target[0][1] - ds[0][1]) + v * (target[1][1] - ds[1][1]);

        let dt_sq = dt_sub * dt_sub;
        let alpha_tilde = if self.compliance <= 0.0 {
            0.0
        } else if dt_sq > 0.0 {
            self.compliance / dt_sq
        } else {
            Real::INFINITY
        };
        let weights = [
            particles[i0].inverse_mass,
            particles[i1].inverse_mass,
            particles[i2].inverse_mass,
        ];
        let corrections = edge_corrections(weights, alpha_tilde, delta1, delta2);
        for (&index, correction) in self.indices.iter().zip(corrections.iter()) {
            particles[index].position += correction;
        }
    }

    fn reset_lambda(&mut self) {
        Self::reset_lambda(self);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const DT: Real = 1.0e-2;

    fn particle(x: Real, y: Real, z: Real, mass: Real) -> SoftParticle {
        SoftParticle::new(Vec3::new(x, y, z), mass)
    }

    /// Right triangle with unit legs.
    fn rest_triangle() -> Vec<SoftParticle> {
        vec![
            particle(0.0, 0.0, 0.0, 1.0),
            particle(1.0, 0.0, 0.0, 1.0),
            particle(0.0, 1.0, 0.0, 1.0),
        ]
    }

    fn scaled(rest: &[SoftParticle], sx: Real, sy: Real) -> Vec<SoftParticle> {
        rest.iter()
            .map(|p| {
                SoftParticle::new(
                    Vec3::new(p.position.x * sx, p.position.y * sy, p.position.z),
                    1.0 / p.inverse_mass.max(1e-300),
                )
            })
            .collect()
    }

    fn constraint(
        rest: &[SoftParticle],
        max: Real,
        min: Real,
        compliance: Real,
    ) -> StrainLimitingConstraint {
        StrainLimitingConstraint::from_particles([0, 1, 2], rest, max, min, compliance)
    }

    /// Edge matrix (columns are the edges) in an orthonormal frame of the triangle.
    fn edge_matrix(t: &[SoftParticle]) -> Mat2 {
        let e1 = t[1].position - t[0].position;
        let e2 = t[2].position - t[0].position;
        let u = e1.normalize();
        let v = e1.cross(&e2).cross(&e1).normalize();
        [[e1.dot(&u), e2.dot(&u)], [e1.dot(&v), e2.dot(&v)]]
    }

    /// Singular values of the deformation gradient of `current` against `rest`,
    /// from the eigenvalues of F^T F (independent of `svd2`).
    fn stretches(rest: &[SoftParticle], current: &[SoftParticle]) -> [Real; 2] {
        let dm = edge_matrix(rest);
        let ds = edge_matrix(current);
        let det = dm[0][0] * dm[1][1] - dm[0][1] * dm[1][0];
        let dm_inv = [
            [dm[1][1] / det, -dm[0][1] / det],
            [-dm[1][0] / det, dm[0][0] / det],
        ];
        let f = mat2_mul(ds, dm_inv);
        let a = f[0][0] * f[0][0] + f[1][0] * f[1][0];
        let b = f[0][0] * f[0][1] + f[1][0] * f[1][1];
        let d = f[0][1] * f[0][1] + f[1][1] * f[1][1];
        let half_trace = 0.5 * (a + d);
        let disc = (half_trace * half_trace - (a * d - b * b)).max(0.0).sqrt();
        [
            (half_trace + disc).sqrt(),
            (half_trace - disc).max(0.0).sqrt(),
        ]
    }

    fn mass_weighted_centroid(t: &[SoftParticle]) -> Vec3 {
        let mut sum = Vec3::zeros();
        let mut total = 0.0;
        for p in t {
            let m = 1.0 / p.inverse_mass.max(1e-300);
            sum += p.position * m;
            total += m;
        }
        sum / total
    }

    fn normal(t: &[SoftParticle]) -> Vec3 {
        (t[1].position - t[0].position).cross(&(t[2].position - t[0].position))
    }

    fn assert_stretches(s: [Real; 2], expected: [Real; 2], tolerance: Real) {
        assert!(
            (s[0] - expected[0]).abs() < tolerance && (s[1] - expected[1]).abs() < tolerance,
            "singular values {s:?}, expected {expected:?}"
        );
    }

    #[test]
    fn svd_reconstructs_general_matrices() {
        for f in [
            [[2.0, 0.0], [0.0, 0.5]],
            [[0.0, -1.0], [1.0, 0.0]],
            [[1.0, 1.0], [0.0, 1.0]],
            [[0.3, -1.7], [2.2, 0.9]],
            [[1.0, 0.0], [0.0, -1.0]],
        ] {
            let svd = svd2(f);
            let sigma = [[svd.s[0], 0.0], [0.0, svd.s[1]]];
            let rebuilt = mat2_mul(mat2_mul(svd.u, sigma), svd.vt);
            for i in 0..2 {
                for j in 0..2 {
                    assert!(
                        (rebuilt[i][j] - f[i][j]).abs() < 1e-12,
                        "F = {f:?}, U S V^T = {rebuilt:?}"
                    );
                }
            }
            assert!(svd.s[0] >= svd.s[1].abs() - 1e-15);
        }
    }

    #[test]
    fn compressed_triangle_recovers_the_minimum_compression() {
        let rest = rest_triangle();
        let mut c = constraint(&rest, 1.1, 0.9, 0.0);
        let mut current = scaled(&rest, 0.5, 0.5);
        let area_before = normal(&current).norm();
        let centroid_before = mass_weighted_centroid(&current);
        c.project(&mut current, DT);
        let s = stretches(&rest, &current);
        assert!(
            s[0] >= 0.9 - 1e-12 && s[1] >= 0.9 - 1e-12,
            "compression must be limited to 0.9: {s:?}"
        );
        assert_stretches(s, [0.9, 0.9], 1e-12);
        assert!(
            normal(&current).norm() > area_before,
            "an over-compressed triangle must grow, not shrink further"
        );
        assert!((mass_weighted_centroid(&current) - centroid_before).norm() < 1e-12);
    }

    #[test]
    fn stretched_triangle_shrinks_exactly_to_the_maximum_stretch() {
        let rest = rest_triangle();
        let mut c = constraint(&rest, 1.1, 0.9, 0.0);
        let mut current = scaled(&rest, 2.0, 2.0);
        let centroid_before = mass_weighted_centroid(&current);
        c.project(&mut current, DT);
        assert_stretches(stretches(&rest, &current), [1.1, 1.1], 1e-12);
        assert!((mass_weighted_centroid(&current) - centroid_before).norm() < 1e-12);
    }

    #[test]
    fn large_violations_do_not_overshoot_or_invert_the_triangle() {
        let rest = rest_triangle();
        let mut c = constraint(&rest, 1.1, 0.9, 0.0);
        let mut current = scaled(&rest, 6.0, 6.0);
        let normal_before = normal(&current);
        c.project(&mut current, DT);
        assert_stretches(stretches(&rest, &current), [1.1, 1.1], 1e-12);
        assert!(
            normal(&current).dot(&normal_before) > 0.0,
            "the triangle must keep its orientation"
        );
    }

    #[test]
    fn only_the_violating_principal_direction_is_corrected() {
        let rest = rest_triangle();
        let mut c = constraint(&rest, 1.1, 0.9, 0.0);
        let mut current = scaled(&rest, 2.0, 1.0);
        c.project(&mut current, DT);
        assert_stretches(stretches(&rest, &current), [1.1, 1.0], 1e-12);
        // Equal masses: the x edge shrinks to 1.1, the centroid stays at (2/3, 1/3).
        let expected = [[0.3, 0.0], [1.4, 0.0], [0.3, 1.0]];
        for (p, want) in current.iter().zip(expected.iter()) {
            assert!(
                (p.position.x - want[0]).abs() < 1e-12 && (p.position.y - want[1]).abs() < 1e-12,
                "position {:?}, expected {want:?}",
                p.position
            );
            assert!(p.position.z.abs() < 1e-14, "the plane must be preserved");
        }
    }

    #[test]
    fn triangle_within_limits_is_untouched_under_rigid_motion() {
        // Skewed, non-axis-aligned rest triangle moved rigidly in space.
        let rest = vec![
            particle(0.0, 0.0, 0.0, 1.0),
            particle(1.2, 0.1, 0.3, 1.0),
            particle(0.2, 0.9, -0.4, 1.0),
        ];
        let mut c = constraint(&rest, 1.1, 0.9, 0.0);
        let axis = Vec3::new(1.0, 2.0, 2.0) / 3.0;
        let (sin_a, cos_a) = 0.7_f64.sin_cos();
        let mut moved: Vec<SoftParticle> = rest
            .iter()
            .map(|p| {
                let x = p.position;
                let r = x * cos_a + axis.cross(&x) * sin_a + axis * (axis.dot(&x) * (1.0 - cos_a));
                SoftParticle::new(r + Vec3::new(3.0, -2.0, 1.0), 1.0)
            })
            .collect();
        let before: Vec<Vec3> = moved.iter().map(|p| p.position).collect();
        c.project(&mut moved, DT);
        for (p, b) in moved.iter().zip(before.iter()) {
            assert!(
                (p.position - b).norm() < 1e-12,
                "rigid motion changed the triangle"
            );
        }
    }

    #[test]
    fn projection_commutes_with_rigid_motions() {
        let rest = vec![
            particle(0.0, 0.0, 0.0, 1.0),
            particle(1.2, 0.1, 0.3, 1.0),
            particle(0.2, 0.9, -0.4, 1.0),
        ];
        let mut c = constraint(&rest, 1.15, 0.85, 0.0);
        // Anisotropic stretch along the first edge, then a rigid motion.
        let stretch = |x: Vec3| -> Vec3 {
            x + (rest[1].position - rest[0].position) * 0.6 * (x.dot(&Vec3::new(1.0, 0.0, 0.0)))
        };
        let base: Vec<SoftParticle> = rest
            .iter()
            .map(|p| SoftParticle::new(stretch(p.position), 1.0))
            .collect();
        let axis = Vec3::new(-1.0, 0.5, 2.0).normalize();
        let (sin_a, cos_a) = 1.1_f64.sin_cos();
        let rotate =
            |x: Vec3| x * cos_a + axis.cross(&x) * sin_a + axis * (axis.dot(&x) * (1.0 - cos_a));
        let shift = Vec3::new(0.4, -1.5, 2.5);

        let mut plain = base.clone();
        c.clone().project(&mut plain, DT);
        let mut moved: Vec<SoftParticle> = base
            .iter()
            .map(|p| SoftParticle::new(rotate(p.position) + shift, 1.0))
            .collect();
        c.project(&mut moved, DT);
        for (a, b) in plain.iter().zip(moved.iter()) {
            let expected = rotate(a.position) + shift;
            assert!(
                (b.position - expected).norm() < 1e-11,
                "projection is not frame invariant"
            );
        }
        assert!(
            stretches(&rest, &plain)[0] <= 1.15 + 1e-10,
            "the stretch must end up inside its limit"
        );
    }

    #[test]
    fn heavy_vertices_move_less_and_the_weighted_centroid_is_conserved() {
        let mut rest = rest_triangle();
        rest[2] = particle(0.0, 1.0, 0.0, 100.0);
        let mut c = constraint(&rest, 1.1, 0.9, 0.0);
        let mut current = scaled(&rest, 1.5, 1.5);
        let before: Vec<Vec3> = current.iter().map(|p| p.position).collect();
        let centroid_before = mass_weighted_centroid(&current);
        c.project(&mut current, DT);
        assert_stretches(stretches(&rest, &current), [1.1, 1.1], 1e-12);
        assert!((mass_weighted_centroid(&current) - centroid_before).norm() < 1e-12);
        let light = (current[0].position - before[0]).norm();
        let heavy = (current[2].position - before[2]).norm();
        assert!(
            heavy < 0.05 * light,
            "the heavy vertex must barely move: light {light}, heavy {heavy}"
        );
    }

    #[test]
    fn pinned_vertex_stays_and_the_limit_is_met_exactly() {
        let mut rest = rest_triangle();
        rest[0] = SoftParticle::new_static(Vec3::zeros());
        let mut c = constraint(&rest, 1.1, 0.9, 0.0);
        let mut current = scaled(&rest, 1.5, 1.5);
        current[0] = SoftParticle::new_static(Vec3::zeros());
        c.project(&mut current, DT);
        assert_eq!(
            current[0].position,
            Vec3::zeros(),
            "the pinned vertex moved"
        );
        assert_stretches(stretches(&rest, &current), [1.1, 1.1], 1e-12);
    }

    #[test]
    fn single_free_vertex_between_two_pins_takes_the_least_squares_step() {
        let mut rest = rest_triangle();
        rest[0] = SoftParticle::new_static(Vec3::zeros());
        rest[1] = SoftParticle::new_static(Vec3::new(1.0, 0.0, 0.0));
        let mut c = constraint(&rest, 1.1, 0.9, 0.0);
        let mut current = scaled(&rest, 1.5, 1.5);
        current[0] = SoftParticle::new_static(Vec3::zeros());
        current[1] = SoftParticle::new_static(Vec3::new(1.5, 0.0, 0.0));
        let before = stretches(&rest, &current);
        c.project(&mut current, DT);
        assert_eq!(current[0].position, Vec3::zeros());
        assert_eq!(current[1].position, Vec3::new(1.5, 0.0, 0.0));
        // The free vertex jumps to p0 + e2*, the target of the edge it controls.
        assert!(
            (current[2].position - Vec3::new(0.0, 1.1, 0.0)).norm() < 1e-12,
            "free vertex at {:?}",
            current[2].position
        );
        let after = stretches(&rest, &current);
        assert!(
            after[0] + after[1] < before[0] + before[1],
            "the violation must shrink: {before:?} -> {after:?}"
        );
        assert!(current[2].position.norm().is_finite());
    }

    #[test]
    fn fully_pinned_triangle_does_not_move() {
        let rest: Vec<SoftParticle> = [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]]
            .iter()
            .map(|p| SoftParticle::new_static(Vec3::new(p[0], p[1], 0.0)))
            .collect();
        let mut c = constraint(&rest, 1.1, 0.9, 0.0);
        let mut current: Vec<SoftParticle> = [[0.0, 0.0], [2.0, 0.0], [0.0, 2.0]]
            .iter()
            .map(|p| SoftParticle::new_static(Vec3::new(p[0], p[1], 0.0)))
            .collect();
        c.project(&mut current, DT);
        assert_eq!(current[1].position, Vec3::new(2.0, 0.0, 0.0));
        assert_eq!(current[2].position, Vec3::new(0.0, 2.0, 0.0));
    }

    #[test]
    fn compliance_softens_the_projection() {
        let rest = rest_triangle();
        let remaining = |compliance: Real, dt: Real| {
            let mut c = constraint(&rest, 1.1, 0.9, compliance);
            let mut current = scaled(&rest, 1.5, 1.5);
            c.project(&mut current, dt);
            stretches(&rest, &current)[0]
        };
        let rigid = remaining(0.0, DT);
        let soft = remaining(1.0e-4, DT);
        let softer = remaining(1.0e-2, DT);
        let softest = remaining(1.0, DT);
        assert!((rigid - 1.1).abs() < 1e-12);
        assert!(
            rigid < soft && soft < softer && softer < softest && softest < 1.5,
            "{rigid} {soft} {softer} {softest}"
        );
        assert!(
            softest > 1.499,
            "huge compliance leaves the triangle alone: {softest}"
        );
        // alpha = compliance / dt^2: a smaller sub-step is softer.
        assert!(remaining(1.0e-4, 5.0e-3) > soft);
        // Zero sub-step: infinite alpha with compliance, no NaN without it.
        assert!((remaining(1.0e-4, 0.0) - 1.5).abs() < 1e-12);
        assert!((remaining(0.0, 0.0) - 1.1).abs() < 1e-12);
    }

    #[test]
    fn degenerate_and_misordered_inputs_do_not_panic() {
        // Collinear current triangle: nothing to project.
        let rest = rest_triangle();
        let mut c = constraint(&rest, 1.1, 0.9, 0.0);
        let mut line = vec![
            particle(0.0, 0.0, 0.0, 1.0),
            particle(1.0, 0.0, 0.0, 1.0),
            particle(2.0, 0.0, 0.0, 1.0),
        ];
        c.project(&mut line, DT);
        assert_eq!(line[2].position, Vec3::new(2.0, 0.0, 0.0));
        // Degenerate rest triangle.
        let mut degenerate = constraint(&line, 1.1, 0.9, 0.0);
        let mut current = scaled(&rest, 1.5, 1.5);
        degenerate.project(&mut current, DT);
        assert!(current.iter().all(|p| p.position.norm().is_finite()));
        // Swapped limits behave as the ordered interval.
        let mut swapped = constraint(&rest, 0.9, 1.1, 0.0);
        let mut stretched = scaled(&rest, 1.5, 1.5);
        swapped.project(&mut stretched, DT);
        assert_stretches(stretches(&rest, &stretched), [1.1, 1.1], 1e-12);
        swapped.reset_lambda();
    }

    #[test]
    fn a_pair_of_triangles_converges_inside_the_limits() {
        // Unit square split along its diagonal, stretched by 50 %.
        let rest = vec![
            particle(0.0, 0.0, 0.0, 1.0),
            particle(1.0, 0.0, 0.0, 1.0),
            particle(1.0, 1.0, 0.0, 1.0),
            particle(0.0, 1.0, 0.0, 1.0),
        ];
        let triangles = [[0usize, 1, 2], [0, 2, 3]];
        let mut constraints: Vec<StrainLimitingConstraint> = triangles
            .iter()
            .map(|&t| StrainLimitingConstraint::from_particles(t, &rest, 1.1, 0.9, 0.0))
            .collect();
        let mut current: Vec<SoftParticle> = rest
            .iter()
            .map(|p| SoftParticle::new(p.position * 1.5, 1.0))
            .collect();
        let centroid_before = mass_weighted_centroid(&current);
        for _ in 0..100 {
            for c in constraints.iter_mut() {
                c.project(&mut current, DT);
            }
        }
        for t in triangles {
            let tri_rest: Vec<SoftParticle> = t.iter().map(|&i| rest[i].clone()).collect();
            let tri_now: Vec<SoftParticle> = t.iter().map(|&i| current[i].clone()).collect();
            let s = stretches(&tri_rest, &tri_now);
            assert!(
                s[0] < 1.1 + 1e-4 && s[1] > 0.9 - 1e-4,
                "triangle {t:?}: {s:?}"
            );
        }
        assert!((mass_weighted_centroid(&current) - centroid_before).norm() < 1e-9);
    }
}
