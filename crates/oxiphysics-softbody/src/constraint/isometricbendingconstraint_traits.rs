//! # IsometricBendingConstraint - Trait Implementations
//!
//! This module contains trait implementations for `IsometricBendingConstraint`.
//!
//! ## Implemented Traits
//!
//! - `SoftConstraint`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use crate::particle::SoftParticle;
use oxiphysics_core::math::{Real, Vec3};

use super::functions::SoftConstraint;
use super::types::IsometricBendingConstraint;

impl SoftConstraint for IsometricBendingConstraint {
    fn project(&mut self, particles: &mut [SoftParticle], dt_sub: Real) {
        let [i0, i1, i2, i3] = self.indices;
        let positions = [
            particles[i0].position,
            particles[i1].position,
            particles[i2].position,
            particles[i3].position,
        ];
        let indices = [i0, i1, i2, i3];
        let c_val = self.constraint_value(particles);
        if c_val.abs() < 1e-14 {
            return;
        }
        let mut grads = [Vec3::zeros(); 4];
        for (a, grad) in grads.iter_mut().enumerate() {
            for (b, pos) in positions.iter().enumerate() {
                *grad += self.q_matrix[a * 4 + b] * pos;
            }
            *grad *= 2.0 * self.stiffness;
        }
        let mut w_sum = 0.0_f64;
        for (a, grad) in grads.iter().enumerate() {
            w_sum += particles[indices[a]].inverse_mass * grad.norm_squared();
        }
        if w_sum < 1e-14 {
            return;
        }
        let alpha_tilde = self.compliance / (dt_sub * dt_sub);
        let delta_lambda = (-c_val - alpha_tilde * self.lambda) / (w_sum + alpha_tilde);
        self.lambda += delta_lambda;
        for a in 0..4 {
            particles[indices[a]].position +=
                grads[a] * (delta_lambda * particles[indices[a]].inverse_mass);
        }
    }

    fn reset_lambda(&mut self) {
        Self::reset_lambda(self);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::constraint::IsometricBendingConstraintV2;

    fn particle(x: Real, y: Real, z: Real) -> SoftParticle {
        SoftParticle::new(Vec3::new(x, y, z), 1.0)
    }

    /// Symmetric kite: p0 = (0,0), p1 = (1,0), wings at (0.5, +-1).
    fn kite() -> Vec<SoftParticle> {
        vec![
            particle(0.0, 0.0, 0.0),
            particle(1.0, 0.0, 0.0),
            particle(0.5, 1.0, 0.0),
            particle(0.5, -1.0, 0.0),
        ]
    }

    /// Planar hinge with scalene triangles.
    fn skewed_flat() -> Vec<SoftParticle> {
        vec![
            particle(0.0, 0.0, 0.0),
            particle(1.3, 0.2, 0.0),
            particle(0.4, 1.1, 0.0),
            particle(1.0, -0.8, 0.0),
        ]
    }

    /// Rotates every position about the axis (1, 2, 2)/3 by `angle` and shifts it.
    fn rigidly_moved(particles: &[SoftParticle], angle: Real, shift: Vec3) -> Vec<SoftParticle> {
        let axis = Vec3::new(1.0, 2.0, 2.0) / 3.0;
        let (s, c) = angle.sin_cos();
        particles
            .iter()
            .map(|p| {
                let x = p.position;
                let rotated = x * c + axis.cross(&x) * s + axis * (axis.dot(&x) * (1.0 - c));
                SoftParticle::new(rotated + shift, 1.0)
            })
            .collect()
    }

    fn centroid(particles: &[SoftParticle]) -> Vec3 {
        particles
            .iter()
            .fold(Vec3::zeros(), |acc, p| acc + p.position)
            / particles.len() as Real
    }

    fn lifted(mut particles: Vec<SoftParticle>, z: Real) -> Vec<SoftParticle> {
        particles[3].position.z += z;
        particles
    }

    fn weights(particles: &[SoftParticle]) -> ([Real; 4], Real) {
        IsometricBendingConstraint::hinge_weights(
            &particles[0].position,
            &particles[1].position,
            &particles[2].position,
            &particles[3].position,
        )
        .expect("non-degenerate hinge")
    }

    #[test]
    fn hinge_weights_of_symmetric_kite() {
        // Both triangles have angles of atan(2) at p0 and p1, so every cotangent is
        // 0.5: K = [1, 1, -1, -1] (the wing angles of 0.75 must not enter), and the
        // total area is 1.
        let (k, scale) = weights(&kite());
        let expected = [1.0, 1.0, -1.0, -1.0];
        for (got, want) in k.iter().zip(expected.iter()) {
            assert!((got - want).abs() < 1e-12, "K = {k:?}");
        }
        assert!((scale - 3.0).abs() < 1e-12, "3 / (A0 + A1) = {scale}");
    }

    #[test]
    fn hinge_weights_sum_to_zero_and_annihilate_planar_positions() {
        for hinge in [kite(), skewed_flat()] {
            let (k, _) = weights(&hinge);
            let sum: Real = k.iter().sum();
            assert!(sum.abs() < 1e-12, "sum(K) = {sum}");
            let mut weighted = Vec3::zeros();
            for (k_a, p) in k.iter().zip(hinge.iter()) {
                weighted += p.position * *k_a;
            }
            assert!(weighted.norm() < 1e-12, "K . p = {weighted:?}");
        }
        // The same holds for a planar hinge in an arbitrary orientation.
        let moved = rigidly_moved(&skewed_flat(), 1.1, Vec3::new(3.0, -2.0, 5.0));
        let (k, _) = weights(&moved);
        let mut weighted = Vec3::zeros();
        for (k_a, p) in k.iter().zip(moved.iter()) {
            weighted += p.position * *k_a;
        }
        assert!(weighted.norm() < 1e-11, "K . p = {weighted:?}");
    }

    #[test]
    fn flat_rest_hinge_has_zero_constraint_value() {
        for hinge in [kite(), skewed_flat()] {
            let c = IsometricBendingConstraint::from_particles([0, 1, 2, 3], &hinge, 2.5, 0.0);
            assert!(c.constraint_value(&hinge).abs() < 1e-12);
        }
    }

    #[test]
    fn flat_rest_hinge_is_not_corrected() {
        let mut hinge = skewed_flat();
        let before: Vec<Vec3> = hinge.iter().map(|p| p.position).collect();
        let mut c = IsometricBendingConstraint::from_particles([0, 1, 2, 3], &hinge, 1.0, 0.0);
        c.project(&mut hinge, 0.01);
        for (p, b) in hinge.iter().zip(before.iter()) {
            assert!(
                (p.position - b).norm() < 1e-12,
                "a flat hinge must stay put"
            );
        }
    }

    #[test]
    fn constraint_value_is_invariant_under_rigid_motion() {
        let rest = kite();
        let c = IsometricBendingConstraint::from_particles([0, 1, 2, 3], &rest, 1.0, 0.0);
        let folded = lifted(kite(), 0.7);
        let reference = c.constraint_value(&folded);
        assert!(
            reference > 1e-3,
            "a folded hinge must violate the constraint"
        );
        let shifted = rigidly_moved(&folded, 0.0, Vec3::new(3.0, -2.0, 5.0));
        let moved = rigidly_moved(&folded, 1.234, Vec3::new(-4.0, 1.5, 0.25));
        assert!(
            (c.constraint_value(&shifted) - reference).abs() < 1e-9 * reference,
            "translation changed C: {} vs {reference}",
            c.constraint_value(&shifted)
        );
        assert!(
            (c.constraint_value(&moved) - reference).abs() < 1e-9 * reference,
            "rigid motion changed C: {} vs {reference}",
            c.constraint_value(&moved)
        );
    }

    #[test]
    fn constraint_value_grows_with_the_fold() {
        let c = IsometricBendingConstraint::from_particles([0, 1, 2, 3], &kite(), 1.0, 0.0);
        let small = c.constraint_value(&lifted(kite(), 0.2));
        let large = c.constraint_value(&lifted(kite(), 0.4));
        assert!(0.0 < small && small < large, "C = {small}, {large}");
    }

    #[test]
    fn projection_flattens_a_folded_hinge_and_conserves_the_centroid() {
        let mut hinge = lifted(kite(), 0.5);
        let mut c = IsometricBendingConstraint::from_particles([0, 1, 2, 3], &kite(), 1.0, 0.0);
        let initial = c.constraint_value(&hinge);
        let centroid_before = centroid(&hinge);
        for _ in 0..60 {
            c.reset_lambda();
            c.project(&mut hinge, 0.01);
        }
        let remaining = c.constraint_value(&hinge);
        assert!(
            remaining < 1e-6 * initial,
            "the hinge must flatten: {initial} -> {remaining}"
        );
        let centroid_after = centroid(&hinge);
        assert!(
            (centroid_after - centroid_before).norm() < 1e-12,
            "bending must conserve linear momentum"
        );
    }

    #[test]
    fn v2_shares_the_corrected_quadratic_form() {
        let hinge = skewed_flat();
        let v1 = IsometricBendingConstraint::from_particles([0, 1, 2, 3], &hinge, 1.0, 0.0);
        let v2 = IsometricBendingConstraintV2::from_particles([0, 1, 2, 3], &hinge, 1.0, 0.0);
        for (a, b) in v1.q_matrix.iter().zip(v2.q_matrix.iter()) {
            assert!((a - b).abs() < 1e-12);
        }
        assert!(
            v2.bending_energy(&hinge).abs() < 1e-12,
            "a flat sheet has no bending energy"
        );
        let energy_folded = v2.bending_energy(&lifted(skewed_flat(), 0.5));
        assert!(energy_folded > 1e-3, "folded hinge energy {energy_folded}");
    }

    #[test]
    fn v2_projection_step_does_not_depend_on_the_stiffness_scale() {
        // C = stiffness * c: the Newton step on C is independent of the scale.
        let rest = kite();
        let mut soft = lifted(kite(), 0.5);
        let mut stiff = soft.clone();
        let mut c_soft =
            IsometricBendingConstraintV2::from_particles([0, 1, 2, 3], &rest, 1.0, 0.0);
        let mut c_stiff =
            IsometricBendingConstraintV2::from_particles([0, 1, 2, 3], &rest, 7.0, 0.0);
        c_soft.project(&mut soft, 0.01);
        c_stiff.project(&mut stiff, 0.01);
        for (a, b) in soft.iter().zip(stiff.iter()) {
            assert!((a.position - b.position).norm() < 1e-12);
        }
    }
}
