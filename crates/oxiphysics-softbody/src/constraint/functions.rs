//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use crate::particle::SoftParticle;
use oxiphysics_core::math::{Real, Vec3};

/// A constraint that can be projected onto a set of particles.
pub trait SoftConstraint: Send + Sync {
    /// Project (correct) particle positions to satisfy this constraint.
    ///
    /// `dt_sub` is the sub-step time step used by XPBD for compliance scaling.
    fn project(&mut self, particles: &mut [SoftParticle], dt_sub: Real);
}
/// Compute the unit normal of a triangle (p0, p1, p2).
///
/// Returns zero vector if the triangle is degenerate.
#[allow(dead_code)]
pub(super) fn triangle_normal(p0: &Vec3, p1: &Vec3, p2: &Vec3) -> Vec3 {
    let n = (p1 - p0).cross(&(p2 - p0));
    let len = n.norm();
    if len < 1e-14 { Vec3::zeros() } else { n / len }
}
/// Compute the area of a triangle (p0, p1, p2).
#[allow(dead_code)]
pub(super) fn triangle_area(p0: &Vec3, p1: &Vec3, p2: &Vec3) -> Real {
    0.5 * (p1 - p0).cross(&(p2 - p0)).norm()
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::BendingConstraint;
    use crate::CollisionConstraint;
    use crate::DistanceConstraint;
    use crate::SoftParticle;
    use crate::VolumeConstraint;
    use crate::constraint::AreaConstraint;

    use crate::constraint::IsometricBendingConstraint;

    use crate::constraint::ShapeMatchingConstraint;
    use crate::constraint::StrainLimitingConstraint;

    use crate::constraint::functions::Vec3;

    fn make_particle(pos: Vec3, mass: Real) -> SoftParticle {
        SoftParticle::new(pos, mass)
    }
    fn make_static(pos: Vec3) -> SoftParticle {
        SoftParticle::new_static(pos)
    }
    #[test]
    fn distance_constraint_maintains_rest_length() {
        let mut particles = vec![
            make_particle(Vec3::new(0.0, 0.0, 0.0), 1.0),
            make_particle(Vec3::new(2.0, 0.0, 0.0), 1.0),
        ];
        let rest = 1.0;
        let mut c = DistanceConstraint::new(0, 1, rest, 0.0);
        for _ in 0..20 {
            c.reset_lambda();
            c.project(&mut particles, 0.01);
        }
        let dist = (particles[0].position - particles[1].position).norm();
        assert!(
            (dist - rest).abs() < 1e-3,
            "Expected dist ≈ {rest}, got {dist}"
        );
    }
    #[test]
    fn distance_constraint_static_particle() {
        let mut particles = vec![
            make_static(Vec3::new(0.0, 0.0, 0.0)),
            make_particle(Vec3::new(2.0, 0.0, 0.0), 1.0),
        ];
        let mut c = DistanceConstraint::new(0, 1, 1.0, 0.0);
        let p0_before = particles[0].position;
        for _ in 0..20 {
            c.reset_lambda();
            c.project(&mut particles, 0.01);
        }
        assert!(
            (particles[0].position - p0_before).norm() < 1e-14,
            "Static particle should not move"
        );
        let dist = (particles[0].position - particles[1].position).norm();
        assert!((dist - 1.0).abs() < 1e-3, "Dist should be ~1: {dist}");
    }
    #[test]
    fn distance_from_particles_sets_rest_length() {
        let particles = vec![
            make_particle(Vec3::new(0.0, 0.0, 0.0), 1.0),
            make_particle(Vec3::new(3.0, 4.0, 0.0), 1.0),
        ];
        let c = DistanceConstraint::from_particles(0, 1, &particles, 0.0);
        assert!((c.rest_length - 5.0).abs() < 1e-10);
    }
    #[test]
    fn bending_constraint_flat_rest_angle() {
        let particles = vec![
            make_particle(Vec3::new(0.0, 0.0, 0.0), 1.0),
            make_particle(Vec3::new(1.0, 0.0, 0.0), 1.0),
            make_particle(Vec3::new(0.5, 1.0, 0.0), 1.0),
            make_particle(Vec3::new(0.5, -1.0, 0.0), 1.0),
        ];
        let c = BendingConstraint::from_particles([0, 1, 2, 3], &particles, 0.0);
        assert!(
            c.rest_angle.abs() < 1e-6 || (c.rest_angle - std::f64::consts::PI).abs() < 1e-6,
            "Coplanar rest angle should be 0 or π: got {}",
            c.rest_angle
        );
    }
    #[test]
    fn bending_constraint_projects_toward_rest() {
        let mut particles = vec![
            make_particle(Vec3::new(0.0, 0.0, 0.0), 1.0),
            make_particle(Vec3::new(1.0, 0.0, 0.0), 1.0),
            make_particle(Vec3::new(0.5, 1.0, 0.0), 1.0),
            make_particle(Vec3::new(0.5, -1.0, 0.5), 1.0),
        ];
        let rest_angle = BendingConstraint::compute_dihedral(
            &Vec3::new(0.0, 0.0, 0.0),
            &Vec3::new(1.0, 0.0, 0.0),
            &Vec3::new(0.5, 1.0, 0.0),
            &Vec3::new(0.5, -1.0, 0.0),
        );
        let mut c = BendingConstraint::new([0, 1, 2, 3], rest_angle, 0.0);
        let initial_angle = BendingConstraint::compute_dihedral(
            &particles[0].position,
            &particles[1].position,
            &particles[2].position,
            &particles[3].position,
        );
        let initial_err = (initial_angle - rest_angle).abs();
        for _ in 0..50 {
            c.reset_lambda();
            c.project(&mut particles, 0.01);
        }
        let final_angle = BendingConstraint::compute_dihedral(
            &particles[0].position,
            &particles[1].position,
            &particles[2].position,
            &particles[3].position,
        );
        let final_err = (final_angle - rest_angle).abs();
        assert!(
            final_err < initial_err,
            "Bending should reduce error: initial={initial_err}, final={final_err}"
        );
    }
    #[test]
    fn isometric_bending_no_correction_at_rest() {
        let particles = vec![
            make_particle(Vec3::new(-1.0, 0.0, 0.0), 1.0),
            make_particle(Vec3::new(1.0, 0.0, 0.0), 1.0),
            make_particle(Vec3::new(0.0, 1.0, 0.0), 1.0),
            make_particle(Vec3::new(0.0, -1.0, 0.0), 1.0),
        ];
        let mut c = IsometricBendingConstraint::from_particles([0, 1, 2, 3], &particles, 1.0, 0.0);
        let mut test_particles = particles.clone();
        let pos_before: Vec<Vec3> = test_particles.iter().map(|p| p.position).collect();
        c.project(&mut test_particles, 0.01);
        let total_move: Real = test_particles
            .iter()
            .zip(pos_before.iter())
            .map(|(p, &b)| (p.position - b).norm())
            .sum();
        assert!(
            total_move < 1e-6,
            "No correction at rest: total move = {total_move}"
        );
    }
    #[test]
    fn isometric_bending_corrects_deformation() {
        let particles = vec![
            make_particle(Vec3::new(0.0, 0.0, 0.0), 1.0),
            make_particle(Vec3::new(1.0, 0.0, 0.0), 1.0),
            make_particle(Vec3::new(0.5, 1.0, 0.0), 1.0),
            make_particle(Vec3::new(0.5, -1.0, 0.0), 1.0),
        ];
        let mut c = IsometricBendingConstraint::from_particles([0, 1, 2, 3], &particles, 1.0, 0.0);
        let mut deformed = particles;
        deformed[3].position = Vec3::new(0.5, -1.0, 1.0);
        let p3_before = deformed[3].position;
        for _ in 0..10 {
            c.reset_lambda();
            c.project(&mut deformed, 0.01);
        }
        let p3_after = deformed[3].position;
        let moved = (p3_after - p3_before).norm();
        assert!(
            moved > 1e-6,
            "Isometric bending should correct deformation: moved={moved}"
        );
    }
    #[test]
    fn area_constraint_preserves_area() {
        let mut particles = vec![
            make_particle(Vec3::new(0.0, 0.0, 0.0), 1.0),
            make_particle(Vec3::new(1.0, 0.0, 0.0), 1.0),
            make_particle(Vec3::new(0.0, 1.0, 0.0), 1.0),
        ];
        let rest_area = AreaConstraint::compute_triangle_area(
            &particles[0].position,
            &particles[1].position,
            &particles[2].position,
        );
        assert!((rest_area - 0.5).abs() < 1e-10, "Rest area = {rest_area}");
        particles[1].position = Vec3::new(2.0, 0.0, 0.0);
        particles[2].position = Vec3::new(0.0, 2.0, 0.0);
        let stretched_area = AreaConstraint::compute_triangle_area(
            &particles[0].position,
            &particles[1].position,
            &particles[2].position,
        );
        assert!(
            stretched_area > rest_area * 1.5,
            "Should be stretched: {stretched_area}"
        );
        let mut c = AreaConstraint::new([0, 1, 2], rest_area, 0.0);
        for _ in 0..50 {
            c.reset_lambda();
            c.project(&mut particles, 0.01);
        }
        let final_area = AreaConstraint::compute_triangle_area(
            &particles[0].position,
            &particles[1].position,
            &particles[2].position,
        );
        let error = (final_area - rest_area).abs();
        assert!(
            error < 0.1,
            "Area should be closer to rest: final={final_area}, rest={rest_area}, error={error}"
        );
    }
    #[test]
    fn area_constraint_from_particles() {
        let particles = vec![
            make_particle(Vec3::new(0.0, 0.0, 0.0), 1.0),
            make_particle(Vec3::new(2.0, 0.0, 0.0), 1.0),
            make_particle(Vec3::new(0.0, 3.0, 0.0), 1.0),
        ];
        let c = AreaConstraint::from_particles([0, 1, 2], &particles, 0.0);
        assert!(
            (c.rest_area - 3.0).abs() < 1e-10,
            "Rest area = {}",
            c.rest_area
        );
    }
    #[test]
    fn shape_matching_rigid_body() {
        let mut particles = vec![
            make_particle(Vec3::new(0.0, 0.0, 0.0), 1.0),
            make_particle(Vec3::new(1.0, 0.0, 0.0), 1.0),
            make_particle(Vec3::new(0.0, 1.0, 0.0), 1.0),
            make_particle(Vec3::new(0.0, 0.0, 1.0), 1.0),
        ];
        let indices = vec![0, 1, 2, 3];
        let mut c = ShapeMatchingConstraint::from_particles(indices, &particles, 1.0);
        particles[1].position = Vec3::new(1.5, 0.3, -0.1);
        let p1_before = particles[1].position;
        c.project(&mut particles, 0.01);
        let p1_after = particles[1].position;
        let moved = (p1_after - p1_before).norm();
        assert!(
            moved > 1e-4,
            "Shape matching should correct deformation: moved={moved}"
        );
    }
    #[test]
    fn shape_matching_preserves_com() {
        let particles = vec![
            make_particle(Vec3::new(0.0, 0.0, 0.0), 1.0),
            make_particle(Vec3::new(1.0, 0.0, 0.0), 1.0),
            make_particle(Vec3::new(0.0, 1.0, 0.0), 1.0),
        ];
        let indices = vec![0, 1, 2];
        let mut c = ShapeMatchingConstraint::from_particles(indices, &particles, 0.5);
        let mut test_particles = particles;
        test_particles[0].position = Vec3::new(-0.3, 0.2, 0.1);
        let com_before = c.current_com(&test_particles);
        c.project(&mut test_particles, 0.01);
        let com_after = c.current_com(&test_particles);
        let com_diff = (com_after - com_before).norm();
        assert!(
            com_diff < 0.5,
            "CoM should be approximately preserved: diff={com_diff}"
        );
    }
    #[test]
    fn shape_matching_static_particles_unchanged() {
        let mut particles = vec![
            make_static(Vec3::new(0.0, 0.0, 0.0)),
            make_particle(Vec3::new(1.0, 0.0, 0.0), 1.0),
            make_particle(Vec3::new(0.0, 1.0, 0.0), 1.0),
        ];
        let indices = vec![0, 1, 2];
        let mut c = ShapeMatchingConstraint::from_particles(indices, &particles, 1.0);
        let p0_before = particles[0].position;
        particles[1].position = Vec3::new(1.5, 0.3, 0.0);
        c.project(&mut particles, 0.01);
        assert!(
            (particles[0].position - p0_before).norm() < 1e-14,
            "Static particle should not move"
        );
    }
    #[test]
    fn strain_limiting_no_correction_within_limits() {
        let particles = vec![
            make_particle(Vec3::new(0.0, 0.0, 0.0), 1.0),
            make_particle(Vec3::new(1.0, 0.0, 0.0), 1.0),
            make_particle(Vec3::new(0.0, 1.0, 0.0), 1.0),
        ];
        let mut c = StrainLimitingConstraint::from_particles([0, 1, 2], &particles, 1.5, 0.5, 0.0);
        let mut test_particles = particles;
        let pos_before: Vec<Vec3> = test_particles.iter().map(|p| p.position).collect();
        c.project(&mut test_particles, 0.01);
        let total_move: Real = test_particles
            .iter()
            .zip(pos_before.iter())
            .map(|(p, &b)| (p.position - b).norm())
            .sum();
        assert!(
            total_move < 1e-10,
            "No correction when strain within limits: total_move={total_move}"
        );
    }
    #[test]
    fn strain_limiting_corrects_excessive_stretch() {
        let particles = vec![
            make_particle(Vec3::new(0.0, 0.0, 0.0), 1.0),
            make_particle(Vec3::new(1.0, 0.0, 0.0), 1.0),
            make_particle(Vec3::new(0.0, 1.0, 0.0), 1.0),
        ];
        let mut c = StrainLimitingConstraint::from_particles([0, 1, 2], &particles, 1.1, 0.9, 0.0);
        let mut stretched = particles;
        stretched[1].position = Vec3::new(3.0, 0.0, 0.0);
        stretched[2].position = Vec3::new(0.0, 3.0, 0.0);
        let p1_before = stretched[1].position;
        for _ in 0..10 {
            c.reset_lambda();
            c.project(&mut stretched, 0.01);
        }
        let p1_after = stretched[1].position;
        let moved = (p1_after - p1_before).norm();
        assert!(
            moved > 1e-4,
            "Should correct excessive stretch: moved={moved}"
        );
    }
    #[test]
    fn volume_constraint_preserves_volume() {
        let mut particles = vec![
            make_particle(Vec3::new(0.0, 0.0, 0.0), 1.0),
            make_particle(Vec3::new(1.0, 0.0, 0.0), 1.0),
            make_particle(Vec3::new(0.0, 1.0, 0.0), 1.0),
            make_particle(Vec3::new(0.0, 0.0, 1.0), 1.0),
        ];
        let rest_vol = VolumeConstraint::compute_tet_volume(
            &particles[0].position,
            &particles[1].position,
            &particles[2].position,
            &particles[3].position,
        );
        assert!(
            (rest_vol - 1.0 / 6.0).abs() < 1e-10,
            "Rest vol = {rest_vol}"
        );
        particles[1].position *= 1.5;
        particles[2].position *= 1.5;
        particles[3].position *= 1.5;
        let mut c = VolumeConstraint::new([0, 1, 2, 3], rest_vol, 0.0);
        for _ in 0..50 {
            c.reset_lambda();
            c.project(&mut particles, 0.01);
        }
        let final_vol = VolumeConstraint::compute_tet_volume(
            &particles[0].position,
            &particles[1].position,
            &particles[2].position,
            &particles[3].position,
        );
        assert!(
            (final_vol - rest_vol).abs() < 0.05,
            "Volume should be near rest: final={final_vol}, rest={rest_vol}"
        );
    }
    #[test]
    fn volume_constraint_from_particles() {
        let particles = vec![
            make_particle(Vec3::new(0.0, 0.0, 0.0), 1.0),
            make_particle(Vec3::new(2.0, 0.0, 0.0), 1.0),
            make_particle(Vec3::new(0.0, 3.0, 0.0), 1.0),
            make_particle(Vec3::new(0.0, 0.0, 4.0), 1.0),
        ];
        let c = VolumeConstraint::from_particles([0, 1, 2, 3], &particles, 0.0);
        assert!(
            (c.rest_volume - 4.0).abs() < 1e-10,
            "Rest volume = {}",
            c.rest_volume
        );
    }
    #[test]
    fn collision_prevents_penetration() {
        let mut particles = vec![make_particle(Vec3::new(0.0, -0.5, 0.0), 1.0)];
        let mut c = CollisionConstraint::new(0, Vec3::new(0.0, 1.0, 0.0), Vec3::new(0.0, 0.0, 0.0));
        c.project(&mut particles, 0.01);
        assert!(
            particles[0].position.y >= -1e-14,
            "Should be above plane: y={}",
            particles[0].position.y
        );
    }
    #[test]
    fn collision_static_particle_unchanged() {
        let mut particles = vec![make_static(Vec3::new(0.0, -0.5, 0.0))];
        let mut c = CollisionConstraint::new(0, Vec3::new(0.0, 1.0, 0.0), Vec3::new(0.0, 0.0, 0.0));
        let before = particles[0].position;
        c.project(&mut particles, 0.01);
        assert!(
            (particles[0].position - before).norm() < 1e-14,
            "Static particle should not move"
        );
    }
    #[test]
    fn collision_above_plane_no_correction() {
        let mut particles = vec![make_particle(Vec3::new(0.0, 1.0, 0.0), 1.0)];
        let mut c = CollisionConstraint::new(0, Vec3::new(0.0, 1.0, 0.0), Vec3::new(0.0, 0.0, 0.0));
        let before = particles[0].position;
        c.project(&mut particles, 0.01);
        assert!(
            (particles[0].position - before).norm() < 1e-14,
            "Above plane: no correction needed"
        );
    }
    #[test]
    fn test_triangle_normal() {
        let n = triangle_normal(
            &Vec3::new(0.0, 0.0, 0.0),
            &Vec3::new(1.0, 0.0, 0.0),
            &Vec3::new(0.0, 1.0, 0.0),
        );
        assert!(
            (n.z.abs() - 1.0).abs() < 1e-10,
            "Normal should be ±z: {:?}",
            n
        );
    }
    #[test]
    fn test_triangle_area() {
        let a = triangle_area(
            &Vec3::new(0.0, 0.0, 0.0),
            &Vec3::new(2.0, 0.0, 0.0),
            &Vec3::new(0.0, 3.0, 0.0),
        );
        assert!((a - 3.0).abs() < 1e-10, "Area should be 3: {a}");
    }
    #[test]
    fn test_triangle_degenerate() {
        let n = triangle_normal(
            &Vec3::new(0.0, 0.0, 0.0),
            &Vec3::new(1.0, 0.0, 0.0),
            &Vec3::new(2.0, 0.0, 0.0),
        );
        assert!(
            n.norm() < 1e-10,
            "Degenerate triangle normal should be zero"
        );
        let a = triangle_area(
            &Vec3::new(0.0, 0.0, 0.0),
            &Vec3::new(1.0, 0.0, 0.0),
            &Vec3::new(2.0, 0.0, 0.0),
        );
        assert!(a.abs() < 1e-10, "Degenerate triangle area should be zero");
    }
    #[test]
    fn multiple_constraints_converge() {
        let mut particles = vec![
            make_particle(Vec3::new(0.0, 0.0, 0.0), 1.0),
            make_particle(Vec3::new(1.0, 0.0, 0.0), 1.0),
            make_particle(Vec3::new(0.5, 1.0, 0.0), 1.0),
        ];
        let mut d01 = DistanceConstraint::from_particles(0, 1, &particles, 0.0);
        let mut d12 = DistanceConstraint::from_particles(1, 2, &particles, 0.0);
        let mut d02 = DistanceConstraint::from_particles(0, 2, &particles, 0.0);
        particles[2].position = Vec3::new(0.5, 1.5, 0.3);
        for _ in 0..100 {
            d01.reset_lambda();
            d12.reset_lambda();
            d02.reset_lambda();
            d01.project(&mut particles, 0.01);
            d12.project(&mut particles, 0.01);
            d02.project(&mut particles, 0.01);
        }
        let dist01 = (particles[0].position - particles[1].position).norm();
        let dist12 = (particles[1].position - particles[2].position).norm();
        let dist02 = (particles[0].position - particles[2].position).norm();
        assert!(
            (dist01 - d01.rest_length).abs() < 0.05,
            "d01: {} vs {}",
            dist01,
            d01.rest_length
        );
        assert!(
            (dist12 - d12.rest_length).abs() < 0.05,
            "d12: {} vs {}",
            dist12,
            d12.rest_length
        );
        assert!(
            (dist02 - d02.rest_length).abs() < 0.05,
            "d02: {} vs {}",
            dist02,
            d02.rest_length
        );
    }
}
/// 3×3 matrix determinant (column-major flat array).
#[allow(dead_code)]
pub(super) fn mat3_det(m: [Real; 9]) -> Real {
    m[0] * (m[4] * m[8] - m[5] * m[7]) - m[3] * (m[1] * m[8] - m[2] * m[7])
        + m[6] * (m[1] * m[5] - m[2] * m[4])
}
/// 3×3 matrix inverse (column-major flat array). Returns None if singular.
#[allow(dead_code)]
pub(super) fn mat3_inverse(m: [Real; 9]) -> Option<[Real; 9]> {
    let det = mat3_det(m);
    if det.abs() < 1e-30 {
        return None;
    }
    let inv_det = 1.0 / det;
    Some([
        (m[4] * m[8] - m[5] * m[7]) * inv_det,
        (m[2] * m[7] - m[1] * m[8]) * inv_det,
        (m[1] * m[5] - m[2] * m[4]) * inv_det,
        (m[5] * m[6] - m[3] * m[8]) * inv_det,
        (m[0] * m[8] - m[2] * m[6]) * inv_det,
        (m[2] * m[3] - m[0] * m[5]) * inv_det,
        (m[3] * m[7] - m[4] * m[6]) * inv_det,
        (m[1] * m[6] - m[0] * m[7]) * inv_det,
        (m[0] * m[4] - m[1] * m[3]) * inv_det,
    ])
}
/// 3×3 matrix multiply (column-major flat arrays): C = A * B.
#[allow(dead_code)]
pub(super) fn mat3_mul(a: [Real; 9], b: [Real; 9]) -> [Real; 9] {
    let mut c = [0.0f64; 9];
    for i in 0..3 {
        for j in 0..3 {
            for k in 0..3 {
                c[i + j * 3] += a[i + k * 3] * b[k + j * 3];
            }
        }
    }
    c
}
#[cfg(test)]
mod new_constraint_tests {
    use super::*;

    use crate::SoftParticle;

    use crate::constraint::BalloonPressureConstraint;

    use crate::constraint::IsometricBendingConstraintV2;
    use crate::constraint::LraConstraint;
    use crate::constraint::MotorConstraint;
    use crate::constraint::NeoHookeanConstraint;
    use crate::constraint::RigidBindingConstraint;

    use crate::constraint::SurfaceTensionConstraint;
    use crate::constraint::functions::Vec3;
    use crate::constraint::mat3_det;
    use crate::constraint::mat3_inverse;
    use crate::constraint::mat3_mul;
    fn p(x: f64, y: f64, z: f64, mass: f64) -> SoftParticle {
        SoftParticle::new(Vec3::new(x, y, z), mass)
    }
    fn static_p(x: f64, y: f64, z: f64) -> SoftParticle {
        SoftParticle::new_static(Vec3::new(x, y, z))
    }
    #[test]
    fn lra_no_correction_within_range() {
        let mut particles = vec![p(0.5, 0.0, 0.0, 1.0)];
        let mut c = LraConstraint::new(0, Vec3::zeros(), 1.0, 0.0);
        let before = particles[0].position;
        c.project(&mut particles, 0.01);
        assert!(
            (particles[0].position - before).norm() < 1e-14,
            "within range: no move"
        );
    }
    #[test]
    fn lra_corrects_when_too_far() {
        let mut particles = vec![p(2.0, 0.0, 0.0, 1.0)];
        let mut c = LraConstraint::new(0, Vec3::zeros(), 1.0, 0.0);
        c.project(&mut particles, 0.01);
        let dist = particles[0].position.norm();
        assert!(dist <= 1.0 + 1e-6, "should be pulled back: dist={dist}");
    }
    #[test]
    fn lra_static_particle_unchanged() {
        let mut particles = vec![static_p(3.0, 0.0, 0.0)];
        let mut c = LraConstraint::new(0, Vec3::zeros(), 1.0, 0.0);
        let before = particles[0].position;
        c.project(&mut particles, 0.01);
        assert!((particles[0].position - before).norm() < 1e-14);
    }
    #[test]
    fn lra_constraint_value_positive_when_too_far() {
        let particles = vec![p(3.0, 0.0, 0.0, 1.0)];
        let c = LraConstraint::new(0, Vec3::zeros(), 1.0, 0.0);
        assert!(c.constraint_value(&particles) > 0.0);
    }
    #[test]
    fn lra_constraint_value_zero_when_within() {
        let particles = vec![p(0.5, 0.0, 0.0, 1.0)];
        let c = LraConstraint::new(0, Vec3::zeros(), 1.0, 0.0);
        assert_eq!(c.constraint_value(&particles), 0.0);
    }
    #[test]
    fn lra_from_particle_sets_correct_distance() {
        let particles = vec![p(2.5, 0.0, 0.0, 1.0)];
        let c = LraConstraint::from_particle(0, &particles, Vec3::zeros(), 0.0);
        assert!((c.max_distance - 2.5).abs() < 1e-10);
    }
    #[test]
    fn lra_reset_lambda() {
        let mut c = LraConstraint::new(0, Vec3::zeros(), 1.0, 0.0);
        c.lambda = 3.125;
        c.reset_lambda();
        assert_eq!(c.lambda, 0.0);
    }
    #[test]
    fn neo_hookean_deformation_gradient_rest_is_identity_ish() {
        let particles = vec![
            p(0.0, 0.0, 0.0, 1.0),
            p(1.0, 0.0, 0.0, 1.0),
            p(0.0, 1.0, 0.0, 1.0),
            p(0.0, 0.0, 1.0, 1.0),
        ];
        let c = NeoHookeanConstraint::from_particles([0, 1, 2, 3], &particles, 1.0, 1.0, 0.0);
        let f = c.deformation_gradient(&particles);
        let j = NeoHookeanConstraint::compute_j(f);
        assert!((j - 1.0).abs() < 1e-8, "rest J should be 1: {j}");
    }
    #[test]
    fn neo_hookean_compute_i1_identity() {
        let f = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0f64];
        let i1 = NeoHookeanConstraint::compute_i1(f);
        assert!((i1 - 3.0).abs() < 1e-10, "I1 of identity is 3: {i1}");
    }
    #[test]
    fn neo_hookean_compute_j_identity() {
        let f = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0f64];
        let j = NeoHookeanConstraint::compute_j(f);
        assert!((j - 1.0).abs() < 1e-10, "det(I) = 1: {j}");
    }
    #[test]
    fn neo_hookean_project_doesnt_panic() {
        let mut particles = vec![
            p(0.0, 0.0, 0.0, 1.0),
            p(1.0, 0.0, 0.0, 1.0),
            p(0.0, 1.0, 0.0, 1.0),
            p(0.0, 0.0, 1.0, 1.0),
        ];
        let mut c = NeoHookeanConstraint::from_particles([0, 1, 2, 3], &particles, 1.0, 1.0, 0.0);
        particles[1].position = Vec3::new(1.5, 0.0, 0.0);
        c.project(&mut particles, 0.01);
    }
    #[test]
    fn neo_hookean_rest_volume_correct() {
        let particles = vec![
            p(0.0, 0.0, 0.0, 1.0),
            p(1.0, 0.0, 0.0, 1.0),
            p(0.0, 1.0, 0.0, 1.0),
            p(0.0, 0.0, 1.0, 1.0),
        ];
        let c = NeoHookeanConstraint::from_particles([0, 1, 2, 3], &particles, 1.0, 1.0, 0.0);
        assert!(
            (c.rest_volume - 1.0 / 6.0).abs() < 1e-8,
            "rest_volume={}",
            c.rest_volume
        );
    }
    #[test]
    fn balloon_triangle_volume_contribution_sign() {
        let p0 = Vec3::new(1.0, 0.0, 0.0);
        let p1 = Vec3::new(0.0, 1.0, 0.0);
        let p2 = Vec3::new(0.0, 0.0, 1.0);
        let v = BalloonPressureConstraint::triangle_volume_contribution(&p0, &p1, &p2);
        assert!(v.abs() > 1e-10, "non-zero volume contribution: {v}");
    }
    #[test]
    fn balloon_project_moves_particles() {
        let mut particles = vec![
            p(1.0, 0.0, 0.0, 1.0),
            p(0.0, 1.0, 0.0, 1.0),
            p(0.0, 0.0, 1.0, 1.0),
        ];
        let before: Vec<Vec3> = particles.iter().map(|p| p.position).collect();
        let mut c = BalloonPressureConstraint::new([0, 1, 2], 2.0, 0.5, 0.0);
        c.project(&mut particles, 0.01);
        let moved: f64 = particles
            .iter()
            .zip(before.iter())
            .map(|(a, &b)| (a.position - b).norm())
            .sum();
        assert!(moved.is_finite(), "positions should remain finite");
    }
    #[test]
    fn balloon_reset_lambda() {
        let mut c = BalloonPressureConstraint::new([0, 1, 2], 1.0, 0.1, 0.0);
        c.lambda = 5.0;
        c.reset_lambda();
        assert_eq!(c.lambda, 0.0);
    }
    #[test]
    fn surface_tension_energy_zero_at_rest() {
        let particles = vec![
            p(0.0, 0.0, 0.0, 1.0),
            p(1.0, 0.0, 0.0, 1.0),
            p(0.0, 1.0, 0.0, 1.0),
        ];
        let c = SurfaceTensionConstraint::from_particles([0, 1, 2], &particles, 0.07, 0.0);
        let energy = c.surface_energy(&particles);
        assert!(
            energy.abs() < 1e-10,
            "energy at rest should be zero: {energy}"
        );
    }
    #[test]
    fn surface_tension_energy_positive_when_stretched() {
        let particles = vec![
            p(0.0, 0.0, 0.0, 1.0),
            p(1.0, 0.0, 0.0, 1.0),
            p(0.0, 1.0, 0.0, 1.0),
        ];
        let c = SurfaceTensionConstraint::from_particles([0, 1, 2], &particles, 0.07, 0.0);
        let mut stretched = particles;
        stretched[1].position = Vec3::new(2.0, 0.0, 0.0);
        stretched[2].position = Vec3::new(0.0, 2.0, 0.0);
        let energy = c.surface_energy(&stretched);
        assert!(
            energy > 0.0,
            "stretched area should have positive energy: {energy}"
        );
    }
    #[test]
    fn surface_tension_rest_area_correct() {
        let particles = vec![
            p(0.0, 0.0, 0.0, 1.0),
            p(2.0, 0.0, 0.0, 1.0),
            p(0.0, 3.0, 0.0, 1.0),
        ];
        let c = SurfaceTensionConstraint::from_particles([0, 1, 2], &particles, 0.07, 0.0);
        assert!(
            (c.rest_area - 3.0).abs() < 1e-10,
            "rest_area={}",
            c.rest_area
        );
    }
    #[test]
    fn surface_tension_current_area_unit_triangle() {
        let p0 = Vec3::new(0.0, 0.0, 0.0);
        let p1 = Vec3::new(1.0, 0.0, 0.0);
        let p2 = Vec3::new(0.0, 1.0, 0.0);
        let area = SurfaceTensionConstraint::current_area(&p0, &p1, &p2);
        assert!(
            (area - 0.5).abs() < 1e-10,
            "unit triangle area = 0.5: {area}"
        );
    }
    #[test]
    fn surface_tension_project_finite() {
        let mut particles = vec![
            p(0.0, 0.0, 0.0, 1.0),
            p(2.0, 0.0, 0.0, 1.0),
            p(0.0, 2.0, 0.0, 1.0),
        ];
        let mut c = SurfaceTensionConstraint::from_particles([0, 1, 2], &particles, 0.07, 0.001);
        particles[1].position = Vec3::new(3.0, 0.0, 0.0);
        c.project(&mut particles, 0.01);
        for pp in &particles {
            assert!(pp.position.x.is_finite());
            assert!(pp.position.y.is_finite());
        }
    }
    #[test]
    fn motor_contraction_ratio_at_zero_actuation() {
        let mut m = MotorConstraint::new(0, 1, 1.0, 0.5, 0.0);
        m.set_actuation(0.0);
        assert!((m.contraction_ratio() - 1.0).abs() < 1e-10);
    }
    #[test]
    fn motor_contraction_ratio_at_full_actuation() {
        let mut m = MotorConstraint::new(0, 1, 1.0, 0.5, 0.0);
        m.set_actuation(1.0);
        assert!(
            (m.contraction_ratio() - 0.5).abs() < 1e-10,
            "ratio={}",
            m.contraction_ratio()
        );
    }
    #[test]
    fn motor_target_length_between_rest_and_min() {
        let mut m = MotorConstraint::new(0, 1, 1.0, 0.6, 0.0);
        m.set_actuation(0.5);
        assert!(m.target_length >= m.rest_length * m.min_ratio - 1e-10);
        assert!(m.target_length <= m.rest_length + 1e-10);
    }
    #[test]
    fn motor_project_drives_towards_target() {
        let mut particles = vec![p(0.0, 0.0, 0.0, 1.0), p(2.0, 0.0, 0.0, 1.0)];
        let mut m = MotorConstraint::from_particles(0, 1, &particles, 0.5, 0.0);
        m.set_actuation(1.0);
        for _ in 0..100 {
            m.reset_lambda();
            m.project(&mut particles, 0.01);
        }
        let dist = (particles[0].position - particles[1].position).norm();
        assert!(dist < 1.5, "motor should reduce gap: dist={dist}");
    }
    #[test]
    fn motor_from_particles_rest_length_correct() {
        let particles = vec![p(0.0, 0.0, 0.0, 1.0), p(3.0, 0.0, 0.0, 1.0)];
        let m = MotorConstraint::from_particles(0, 1, &particles, 0.5, 0.0);
        assert!((m.rest_length - 3.0).abs() < 1e-10);
    }
    #[test]
    fn motor_actuation_clamps_to_unit() {
        let mut m = MotorConstraint::new(0, 1, 1.0, 0.5, 0.0);
        m.set_actuation(5.0);
        assert!((m.actuation - 1.0).abs() < 1e-10);
        m.set_actuation(-1.0);
        assert!(m.actuation.abs() < 1e-10);
    }
    #[test]
    fn rigid_binding_pulls_particle_to_target() {
        let mut particles = vec![p(2.0, 0.0, 0.0, 1.0)];
        let mut c = RigidBindingConstraint::new(0, Vec3::zeros(), 0.0);
        for _ in 0..100 {
            c.reset_lambda();
            c.project(&mut particles, 0.01);
        }
        assert!(
            particles[0].position.norm() < 0.1,
            "should converge to target"
        );
    }
    #[test]
    fn rigid_binding_residual_at_target_is_zero() {
        let particles = vec![p(0.0, 0.0, 0.0, 1.0)];
        let c = RigidBindingConstraint::new(0, Vec3::zeros(), 0.0);
        assert!(c.residual(&particles) < 1e-10);
    }
    #[test]
    fn rigid_binding_residual_when_far() {
        let particles = vec![p(3.0, 4.0, 0.0, 1.0)];
        let c = RigidBindingConstraint::new(0, Vec3::zeros(), 0.0);
        assert!((c.residual(&particles) - 5.0).abs() < 1e-10);
    }
    #[test]
    fn rigid_binding_set_target() {
        let mut c = RigidBindingConstraint::new(0, Vec3::zeros(), 0.0);
        c.set_target(Vec3::new(1.0, 2.0, 3.0));
        assert!((c.target.x - 1.0).abs() < 1e-10);
    }
    #[test]
    fn rigid_binding_static_particle_unchanged() {
        let mut particles = vec![static_p(3.0, 0.0, 0.0)];
        let mut c = RigidBindingConstraint::new(0, Vec3::zeros(), 0.0);
        let before = particles[0].position;
        c.project(&mut particles, 0.01);
        assert!((particles[0].position - before).norm() < 1e-14);
    }
    #[test]
    fn iso_bending_v2_energy_near_zero_at_rest() {
        let particles = vec![
            p(-1.0, 0.0, 0.0, 1.0),
            p(1.0, 0.0, 0.0, 1.0),
            p(0.0, 1.0, 0.0, 1.0),
            p(0.0, -1.0, 0.0, 1.0),
        ];
        let c = IsometricBendingConstraintV2::from_particles([0, 1, 2, 3], &particles, 1.0, 0.001);
        let energy = c.bending_energy(&particles);
        assert!(energy.is_finite(), "energy should be finite: {energy}");
    }
    #[test]
    fn iso_bending_v2_project_finite() {
        let mut particles = vec![
            p(-1.0, 0.0, 0.0, 1.0),
            p(1.0, 0.0, 0.0, 1.0),
            p(0.0, 1.0, 0.5, 1.0),
            p(0.0, -1.0, -0.5, 1.0),
        ];
        let mut c =
            IsometricBendingConstraintV2::from_particles([0, 1, 2, 3], &particles, 1.0, 0.001);
        for _ in 0..10 {
            c.reset_lambda();
            c.project(&mut particles, 0.01);
        }
        for pp in &particles {
            assert!(pp.position.norm().is_finite());
        }
    }
    #[test]
    fn iso_bending_v2_reset_lambda() {
        let particles = vec![
            p(-1.0, 0.0, 0.0, 1.0),
            p(1.0, 0.0, 0.0, 1.0),
            p(0.0, 1.0, 0.0, 1.0),
            p(0.0, -1.0, 0.0, 1.0),
        ];
        let mut c =
            IsometricBendingConstraintV2::from_particles([0, 1, 2, 3], &particles, 1.0, 0.0);
        c.lambda = 99.0;
        c.reset_lambda();
        assert_eq!(c.lambda, 0.0);
    }
    #[test]
    fn mat3_det_identity() {
        let m = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0f64];
        assert!((mat3_det(m) - 1.0).abs() < 1e-10);
    }
    #[test]
    fn mat3_det_zero_for_singular() {
        let m = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0f64];
        assert!(mat3_det(m).abs() < 1e-10);
    }
    #[test]
    fn mat3_inverse_identity() {
        let m = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0f64];
        let inv = mat3_inverse(m).unwrap();
        for i in 0..9 {
            assert!((inv[i] - m[i]).abs() < 1e-10, "inv[{i}]={}", inv[i]);
        }
    }
    #[test]
    fn mat3_inverse_singular_returns_none() {
        let m = [1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0f64];
        assert!(mat3_inverse(m).is_none());
    }
    #[test]
    fn mat3_mul_identity() {
        let id = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0f64];
        let a = [2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0f64];
        let c = mat3_mul(a, id);
        for i in 0..9 {
            assert!((c[i] - a[i]).abs() < 1e-10, "c[{i}]={}", c[i]);
        }
    }
    #[test]
    fn mat3_mul_returns_correct_shape() {
        let a = [1.0f64; 9];
        let b = [1.0f64; 9];
        let c = mat3_mul(a, b);
        assert_eq!(c.len(), 9);
    }
    #[test]
    fn multiple_constraint_types_convergence() {
        let mut particles = vec![
            p(0.0, 0.0, 0.0, 1.0),
            p(1.0, 0.0, 0.0, 1.0),
            p(0.5, 1.0, 0.0, 1.0),
        ];
        let anchor = Vec3::zeros();
        let mut lra = LraConstraint::new(0, anchor, 0.5, 0.0);
        let mut motor = MotorConstraint::from_particles(0, 1, &particles, 0.8, 0.01);
        motor.set_actuation(0.5);
        for _ in 0..50 {
            lra.reset_lambda();
            motor.reset_lambda();
            lra.project(&mut particles, 0.01);
            motor.project(&mut particles, 0.01);
        }
        for pp in &particles {
            assert!(
                pp.position.norm().is_finite(),
                "positions must remain finite"
            );
        }
    }
    #[test]
    fn rigid_binding_convergence_with_compliance() {
        let mut particles = vec![p(5.0, 0.0, 0.0, 1.0)];
        let target = Vec3::new(1.0, 0.0, 0.0);
        let mut c = RigidBindingConstraint::new(0, target, 0.001);
        for _ in 0..200 {
            c.reset_lambda();
            c.project(&mut particles, 0.01);
        }
        let dist = (particles[0].position - target).norm();
        assert!(
            dist < 4.5,
            "should move closer with compliance: dist={dist}"
        );
    }
    #[test]
    fn surface_tension_area_conservation_iterations() {
        let particles_orig = vec![
            p(0.0, 0.0, 0.0, 1.0),
            p(1.0, 0.0, 0.0, 1.0),
            p(0.0, 1.0, 0.0, 1.0),
        ];
        let c_ref = SurfaceTensionConstraint::from_particles([0, 1, 2], &particles_orig, 0.07, 0.0);
        let mut particles = particles_orig;
        particles[2].position = Vec3::new(0.0, 2.0, 0.0);
        let mut c = SurfaceTensionConstraint {
            rest_area: c_ref.rest_area,
            ..c_ref
        };
        for _ in 0..100 {
            c.reset_lambda();
            c.project(&mut particles, 0.01);
        }
        let final_area = SurfaceTensionConstraint::current_area(
            &particles[0].position,
            &particles[1].position,
            &particles[2].position,
        );
        assert!(final_area.is_finite(), "area={final_area}");
    }
}
