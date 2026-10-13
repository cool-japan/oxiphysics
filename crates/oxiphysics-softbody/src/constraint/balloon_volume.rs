// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Whole-balloon enclosed-volume constraint (XPBD).
//!
//! The pressure state of a closed membrane is a *global* quantity: the volume
//! enclosed by the whole surface.  [`BalloonVolumeConstraint`] therefore owns
//! every surface triangle and enforces
//!
//! ```text
//! C(x) = V(x) - pressure * rest_volume,      V(x) = 1/6 * sum_t p0 . (p1 x p2)
//! ```
//!
//! with one XPBD multiplier for the entire balloon.  The vertex gradients are
//! `dV/dp0 = (p1 x p2) / 6`, `dV/dp1 = (p2 x p0) / 6` and
//! `dV/dp2 = (p0 x p1) / 6`, summed over all triangles incident to a vertex.

use super::functions::SoftConstraint;
use crate::particle::SoftParticle;
use crate::{Error, Result};
use oxiphysics_core::math::{Real, Vec3};

/// Weighted squared gradient norm (plus compliance) below which no correction
/// is applied because the constraint cannot move any vertex.
const MIN_GRADIENT_ENERGY: Real = 1e-14;

/// `true` when `value` is a number strictly greater than `floor` (`NaN` is not).
fn exceeds(value: Real, floor: Real) -> bool {
    value > floor
}

/// Enclosed-volume (pressure) constraint for a closed triangle mesh.
///
/// `triangles` must form a closed surface whose triangles are wound
/// counter-clockwise when seen from outside, so that the signed volume
/// `V = 1/6 * sum p0 . (p1 x p2)` is positive and the vertex gradients point
/// outward.  For a closed surface `V` and all its gradients are independent of
/// the origin, so the constraint is translation invariant.
///
/// The target volume is `pressure * rest_volume`: `pressure = 1` keeps the rest
/// volume, `pressure = 1.2` inflates the balloon to 120 % of it.  `rest_volume`
/// is a *signed* volume: use the value returned by
/// [`BalloonVolumeConstraint::enclosed_volume`] on the rest shape, which is
/// negative for an inward-wound mesh and then still inflates correctly.
///
/// This replaces the per-triangle
/// [`BalloonPressureConstraint`](super::types::BalloonPressureConstraint), which
/// compared the volume of a single triangle with the volume of the whole
/// balloon.
#[derive(Debug, Clone)]
pub struct BalloonVolumeConstraint {
    /// Surface triangles (indices into the particle array).
    triangles: Vec<[usize; 3]>,
    /// Sorted, unique particle indices referenced by `triangles`.
    vertices: Vec<usize>,
    /// For every triangle corner, the position of its particle in `vertices`.
    corner_slots: Vec<[usize; 3]>,
    /// Scratch: volume gradient of every entry of `vertices`.
    gradients: Vec<Vec3>,
    /// Signed enclosed volume of the rest shape.
    pub rest_volume: Real,
    /// Ratio of the target volume to `rest_volume`.
    pub pressure: Real,
    /// XPBD compliance of the volume constraint in `m^3 / Pa` (inverse of the
    /// stiffness of the enclosed gas); `0` enforces the target volume rigidly.
    pub compliance: Real,
    /// Accumulated XPBD Lagrange multiplier of the whole balloon.
    lambda: Real,
}

impl BalloonVolumeConstraint {
    /// Create a constraint with an explicit signed rest volume.
    pub fn new(
        triangles: Vec<[usize; 3]>,
        rest_volume: Real,
        pressure: Real,
        compliance: Real,
    ) -> Self {
        let mut vertices: Vec<usize> = triangles.iter().flatten().copied().collect();
        vertices.sort_unstable();
        vertices.dedup();
        let slot_of = |index: usize| vertices.binary_search(&index).unwrap_or_else(|slot| slot);
        let corner_slots = triangles
            .iter()
            .map(|tri| [slot_of(tri[0]), slot_of(tri[1]), slot_of(tri[2])])
            .collect();
        let gradients = vec![Vec3::zeros(); vertices.len()];
        Self {
            triangles,
            vertices,
            corner_slots,
            gradients,
            rest_volume,
            pressure,
            compliance,
            lambda: 0.0,
        }
    }

    /// Create a constraint whose rest volume is the volume currently enclosed
    /// by `triangles`.
    ///
    /// # Errors
    ///
    /// Fails when a triangle references a particle outside `particles`.
    pub fn from_particles(
        triangles: Vec<[usize; 3]>,
        particles: &[SoftParticle],
        pressure: Real,
        compliance: Real,
    ) -> Result<Self> {
        let mut constraint = Self::new(triangles, 0.0, pressure, compliance);
        constraint.rest_volume = constraint.enclosed_volume(particles).ok_or_else(|| {
            Error::General(
                "balloon triangle references a particle outside the particle array".to_string(),
            )
        })?;
        Ok(constraint)
    }

    /// Surface triangles of the balloon.
    pub fn triangles(&self) -> &[[usize; 3]] {
        &self.triangles
    }

    /// Accumulated Lagrange multiplier of the current sub-step.
    pub fn lambda(&self) -> Real {
        self.lambda
    }

    /// Reset the Lagrange multiplier.
    pub fn reset_lambda(&mut self) {
        self.lambda = 0.0;
    }

    /// Volume the balloon is driven towards: `pressure * rest_volume`.
    pub fn target_volume(&self) -> Real {
        self.pressure * self.rest_volume
    }

    /// Whether every referenced particle exists in a particle array of length `len`.
    fn indices_in_range(&self, len: usize) -> bool {
        self.vertices.last().is_none_or(|&largest| largest < len)
    }

    /// Signed volume enclosed by the surface, or `None` when a triangle
    /// references a particle outside `particles`.
    ///
    /// The divergence theorem is evaluated relative to the centroid of the
    /// referenced particles.  For a closed surface that equals the textbook
    /// `1/6 * sum p0 . (p1 x p2)`, but it keeps the sum free of cancellation
    /// when the balloon is far from the origin.
    pub fn enclosed_volume(&self, particles: &[SoftParticle]) -> Option<Real> {
        if !self.indices_in_range(particles.len()) {
            return None;
        }
        let centroid = self.centroid(particles)?;
        let volume: Real = self
            .triangles
            .iter()
            .map(|tri| {
                let a = particles[tri[0]].position - centroid;
                let b = particles[tri[1]].position - centroid;
                let c = particles[tri[2]].position - centroid;
                a.dot(&b.cross(&c))
            })
            .sum();
        Some(volume / 6.0)
    }

    /// Centroid of the referenced particles, `None` for an empty constraint.
    fn centroid(&self, particles: &[SoftParticle]) -> Option<Vec3> {
        if self.vertices.is_empty() {
            return None;
        }
        let sum = self
            .vertices
            .iter()
            .fold(Vec3::zeros(), |acc, &i| acc + particles[i].position);
        Some(sum / self.vertices.len() as Real)
    }

    /// Fill `self.gradients` with `dV/dp` for every referenced particle and
    /// return the signed volume.  The particle indices must be in range.
    fn accumulate_gradients(&mut self, particles: &[SoftParticle]) -> Real {
        let Some(centroid) = self.centroid(particles) else {
            return 0.0;
        };
        self.gradients.clear();
        self.gradients.resize(self.vertices.len(), Vec3::zeros());
        let mut volume = 0.0;
        for (tri, slots) in self.triangles.iter().zip(&self.corner_slots) {
            let a = particles[tri[0]].position - centroid;
            let b = particles[tri[1]].position - centroid;
            let c = particles[tri[2]].position - centroid;
            let bc = b.cross(&c);
            volume += a.dot(&bc);
            self.gradients[slots[0]] += bc;
            self.gradients[slots[1]] += c.cross(&a);
            self.gradients[slots[2]] += a.cross(&b);
        }
        for gradient in &mut self.gradients {
            *gradient /= 6.0;
        }
        volume / 6.0
    }
}

impl SoftConstraint for BalloonVolumeConstraint {
    /// Apply one XPBD update that moves every balloon vertex along its volume
    /// gradient.  The constraint does nothing for a non-positive `dt_sub`, when
    /// a triangle references a missing particle, or when no vertex can move.
    fn project(&mut self, particles: &mut [SoftParticle], dt_sub: Real) {
        let usable = dt_sub.is_finite() && dt_sub > 0.0 && self.indices_in_range(particles.len());
        if !usable {
            return;
        }
        let volume = self.accumulate_gradients(particles);
        let c = volume - self.target_volume();
        let weighted_gradient_norm: Real = self
            .vertices
            .iter()
            .zip(&self.gradients)
            .map(|(&i, gradient)| particles[i].inverse_mass * gradient.norm_squared())
            .sum();
        let alpha_tilde = self.compliance.max(0.0) / (dt_sub * dt_sub);
        let denominator = weighted_gradient_norm + alpha_tilde;
        if !exceeds(denominator, MIN_GRADIENT_ENERGY) {
            return;
        }
        let delta_lambda = (-c - alpha_tilde * self.lambda) / denominator;
        if !delta_lambda.is_finite() {
            return;
        }
        self.lambda += delta_lambda;
        for (&i, gradient) in self.vertices.iter().zip(&self.gradients) {
            let weight = particles[i].inverse_mass;
            particles[i].position += gradient * (weight * delta_lambda);
        }
    }

    fn reset_lambda(&mut self) {
        Self::reset_lambda(self);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::BTreeMap;

    /// Unit icosphere subdivided `levels` times; all triangles wound outward.
    fn icosphere(levels: usize) -> (Vec<Vec3>, Vec<[usize; 3]>) {
        let t = (1.0 + 5.0_f64.sqrt()) / 2.0;
        let mut points: Vec<Vec3> = [
            [-1.0, t, 0.0],
            [1.0, t, 0.0],
            [-1.0, -t, 0.0],
            [1.0, -t, 0.0],
            [0.0, -1.0, t],
            [0.0, 1.0, t],
            [0.0, -1.0, -t],
            [0.0, 1.0, -t],
            [t, 0.0, -1.0],
            [t, 0.0, 1.0],
            [-t, 0.0, -1.0],
            [-t, 0.0, 1.0],
        ]
        .iter()
        .map(|p| Vec3::new(p[0], p[1], p[2]).normalize())
        .collect();
        let mut faces: Vec<[usize; 3]> = vec![
            [0, 11, 5],
            [0, 5, 1],
            [0, 1, 7],
            [0, 7, 10],
            [0, 10, 11],
            [1, 5, 9],
            [5, 11, 4],
            [11, 10, 2],
            [10, 7, 6],
            [7, 1, 8],
            [3, 9, 4],
            [3, 4, 2],
            [3, 2, 6],
            [3, 6, 8],
            [3, 8, 9],
            [4, 9, 5],
            [2, 4, 11],
            [6, 2, 10],
            [8, 6, 7],
            [9, 8, 1],
        ];
        for _ in 0..levels {
            let mut midpoints: BTreeMap<(usize, usize), usize> = BTreeMap::new();
            let mut next_faces = Vec::with_capacity(faces.len() * 4);
            for &[a, b, c] in &faces {
                let mut mid = |u: usize, v: usize| -> usize {
                    *midpoints.entry((u.min(v), u.max(v))).or_insert_with(|| {
                        points.push(((points[u] + points[v]) * 0.5).normalize());
                        points.len() - 1
                    })
                };
                let ab = mid(a, b);
                let bc = mid(b, c);
                let ca = mid(c, a);
                next_faces.extend([[a, ab, ca], [b, bc, ab], [c, ca, bc], [ab, bc, ca]]);
            }
            faces = next_faces;
        }
        for face in &mut faces {
            let [a, b, c] = *face;
            let outward = (points[b] - points[a]).cross(&(points[c] - points[a]));
            if outward.dot(&(points[a] + points[b] + points[c])) < 0.0 {
                face.swap(1, 2);
            }
        }
        (points, faces)
    }

    fn particles_from(points: &[Vec3], offset: Vec3) -> Vec<SoftParticle> {
        points
            .iter()
            .map(|&p| SoftParticle::new(p + offset, 1.0))
            .collect()
    }

    fn cube() -> (Vec<Vec3>, Vec<[usize; 3]>) {
        let points = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [1.0, 0.0, 1.0],
            [1.0, 1.0, 1.0],
            [0.0, 1.0, 1.0],
        ]
        .iter()
        .map(|p| Vec3::new(p[0], p[1], p[2]))
        .collect();
        let faces = vec![
            [0, 3, 2],
            [0, 2, 1],
            [4, 5, 6],
            [4, 6, 7],
            [0, 1, 5],
            [0, 5, 4],
            [2, 3, 7],
            [2, 7, 6],
            [1, 2, 6],
            [1, 6, 5],
            [3, 0, 4],
            [3, 4, 7],
        ];
        (points, faces)
    }

    fn run(
        constraint: &mut BalloonVolumeConstraint,
        particles: &mut [SoftParticle],
        iterations: usize,
    ) {
        for _ in 0..iterations {
            constraint.project(particles, 1.0 / 60.0);
        }
    }

    #[test]
    fn enclosed_volume_of_unit_cube_is_one_and_translation_invariant() {
        let (points, faces) = cube();
        for offset in [Vec3::zeros(), Vec3::new(1.0e3, -2.0e3, 5.0e2)] {
            let particles = particles_from(&points, offset);
            let constraint = BalloonVolumeConstraint::new(faces.clone(), 1.0, 1.0, 0.0);
            let volume = constraint.enclosed_volume(&particles).expect("valid mesh");
            assert!(
                (volume - 1.0).abs() < 1e-9,
                "cube volume {volume} at {offset:?}"
            );
        }
    }

    #[test]
    fn icosphere_volume_approaches_sphere_volume() {
        let (points, faces) = icosphere(3);
        let particles = particles_from(&points, Vec3::zeros());
        let constraint = BalloonVolumeConstraint::new(faces, 1.0, 1.0, 0.0);
        let volume = constraint.enclosed_volume(&particles).expect("valid mesh");
        let sphere = 4.0 / 3.0 * std::f64::consts::PI;
        assert!(
            volume > 0.0 && (sphere - volume) / sphere < 0.01,
            "volume {volume}"
        );
    }

    #[test]
    fn gradients_match_finite_differences_and_sum_to_zero() {
        let (points, faces) = icosphere(1);
        let mut particles = particles_from(&points, Vec3::new(0.4, -0.3, 0.9));
        // Deform the sphere into a generic closed shape.
        for (k, p) in particles.iter_mut().enumerate() {
            let phase = k as f64;
            p.position += Vec3::new(
                0.05 * phase.sin(),
                0.04 * (2.0 * phase).cos(),
                0.03 * (3.0 * phase).sin(),
            );
        }
        let mut constraint = BalloonVolumeConstraint::new(faces, 1.0, 1.0, 0.0);
        let volume = constraint.accumulate_gradients(&particles);
        let analytic = constraint.gradients.clone();
        assert!((volume - constraint.enclosed_volume(&particles).expect("valid")).abs() < 1e-12);
        let total = analytic.iter().fold(Vec3::zeros(), |acc, g| acc + g);
        assert!(
            total.norm() < 1e-12,
            "translation invariance violated: {total:?}"
        );
        let h = 1e-6;
        for (slot, &index) in constraint.vertices.clone().iter().enumerate() {
            for (axis, &expected) in analytic[slot].iter().enumerate() {
                let original = particles[index].position[axis];
                particles[index].position[axis] = original + h;
                let plus = constraint.enclosed_volume(&particles).expect("valid");
                particles[index].position[axis] = original - h;
                let minus = constraint.enclosed_volume(&particles).expect("valid");
                particles[index].position[axis] = original;
                let numeric = (plus - minus) / (2.0 * h);
                assert!(
                    (numeric - expected).abs() < 1e-7,
                    "vertex {index} axis {axis}: numeric {numeric} vs analytic {expected}"
                );
            }
        }
    }

    #[test]
    fn closed_icosphere_inflates_to_requested_pressure_within_one_percent() {
        let (points, faces) = icosphere(2);
        let rest = particles_from(&points, Vec3::zeros());
        let rest_volume = BalloonVolumeConstraint::new(faces.clone(), 1.0, 1.0, 0.0)
            .enclosed_volume(&rest)
            .expect("valid mesh");
        let mut results = Vec::new();
        for offset in [
            Vec3::zeros(),
            Vec3::new(7.3, -4.1, 12.9),
            Vec3::new(-250.0, 80.0, 31.0),
        ] {
            let mut particles = particles_from(&points, offset);
            let mut constraint = BalloonVolumeConstraint::new(faces.clone(), rest_volume, 1.2, 0.0);
            run(&mut constraint, &mut particles, 40);
            let volume = constraint.enclosed_volume(&particles).expect("valid mesh");
            assert!(
                (volume / rest_volume - 1.2).abs() < 0.01 * 1.2,
                "V / V0 = {} should be 1.2 (offset {offset:?})",
                volume / rest_volume
            );
            results.push(
                particles
                    .iter()
                    .map(|p| p.position - offset)
                    .collect::<Vec<_>>(),
            );
        }
        // The result does not depend on where the balloon is.
        for other in &results[1..] {
            for (a, b) in results[0].iter().zip(other) {
                assert!(
                    (a - b).norm() < 1e-9,
                    "inflation depends on the position: {a:?} vs {b:?}"
                );
            }
        }
        // A uniform inflation moves every vertex outward and keeps it spherical.
        let scale = 1.2_f64.cbrt();
        for p in &results[0] {
            assert!(
                (p.norm() - scale).abs() < 0.02,
                "radius {} vs {scale}",
                p.norm()
            );
        }
    }

    #[test]
    fn deflation_and_rest_pressure() {
        let (points, faces) = icosphere(2);
        let rest = particles_from(&points, Vec3::zeros());
        let rest_volume = BalloonVolumeConstraint::new(faces.clone(), 1.0, 1.0, 0.0)
            .enclosed_volume(&rest)
            .expect("valid mesh");

        let mut particles = particles_from(&points, Vec3::zeros());
        let mut deflate = BalloonVolumeConstraint::new(faces.clone(), rest_volume, 0.6, 0.0);
        run(&mut deflate, &mut particles, 40);
        let volume = deflate.enclosed_volume(&particles).expect("valid mesh");
        assert!(
            (volume / rest_volume - 0.6).abs() < 0.006,
            "V / V0 = {}",
            volume / rest_volume
        );

        let mut particles = particles_from(&points, Vec3::zeros());
        let mut hold = BalloonVolumeConstraint::new(faces, rest_volume, 1.0, 0.0);
        run(&mut hold, &mut particles, 10);
        for (p, q) in particles.iter().zip(&rest) {
            assert!(
                (p.position - q.position).norm() < 1e-12,
                "the rest volume is a fixed point"
            );
        }
    }

    #[test]
    fn from_particles_uses_the_enclosed_volume_and_handles_inward_winding() {
        let (points, mut faces) = icosphere(2);
        let particles = particles_from(&points, Vec3::new(2.0, 1.0, -3.0));
        let outward = BalloonVolumeConstraint::from_particles(faces.clone(), &particles, 1.3, 0.0)
            .expect("valid mesh");
        assert!(outward.rest_volume > 0.0);
        assert!((outward.target_volume() - 1.3 * outward.rest_volume).abs() < 1e-12);

        for face in &mut faces {
            face.swap(1, 2);
        }
        let mut inward = BalloonVolumeConstraint::from_particles(faces, &particles, 1.3, 0.0)
            .expect("valid mesh");
        assert!(
            (inward.rest_volume + outward.rest_volume).abs() < 1e-9,
            "inward winding gives a negative signed volume"
        );
        let mut inflated = particles.clone();
        run(&mut inward, &mut inflated, 40);
        let ratio = inward.enclosed_volume(&inflated).expect("valid mesh") / inward.rest_volume;
        assert!((ratio - 1.3).abs() < 0.013, "|V| / |V0| = {ratio}");
        let radius: f64 = inflated
            .iter()
            .map(|p| (p.position - Vec3::new(2.0, 1.0, -3.0)).norm())
            .sum::<f64>()
            / inflated.len() as f64;
        assert!(
            radius > 1.05,
            "the inward-wound balloon must still inflate, radius {radius}"
        );
    }

    #[test]
    fn centroid_is_preserved_for_equal_masses_and_momentum_for_unequal_masses() {
        let (points, faces) = icosphere(2);
        let rest = particles_from(&points, Vec3::zeros());
        let rest_volume = BalloonVolumeConstraint::new(faces.clone(), 1.0, 1.0, 0.0)
            .enclosed_volume(&rest)
            .expect("valid mesh");
        let mut particles = particles_from(&points, Vec3::new(3.0, 3.0, 3.0));
        let before = particles
            .iter()
            .fold(Vec3::zeros(), |acc, p| acc + p.position)
            / particles.len() as f64;
        let mut constraint = BalloonVolumeConstraint::new(faces.clone(), rest_volume, 1.5, 0.0);
        run(&mut constraint, &mut particles, 20);
        let after = particles
            .iter()
            .fold(Vec3::zeros(), |acc, p| acc + p.position)
            / particles.len() as f64;
        assert!(
            (after - before).norm() < 1e-10,
            "centroid moved by {:?}",
            after - before
        );

        // Unequal masses: the mass-weighted displacement sums to zero.
        let mut heavy = particles_from(&points, Vec3::zeros());
        for (k, p) in heavy.iter_mut().enumerate() {
            let mass = 1.0 + (k % 5) as f64;
            *p = SoftParticle::new(p.position, mass);
        }
        let start: Vec<Vec3> = heavy.iter().map(|p| p.position).collect();
        let mut constraint = BalloonVolumeConstraint::new(faces, rest_volume, 1.2, 0.0);
        run(&mut constraint, &mut heavy, 1);
        let momentum = heavy.iter().zip(&start).fold(Vec3::zeros(), |acc, (p, s)| {
            acc + (p.position - s) / p.inverse_mass
        });
        assert!(
            momentum.norm() < 1e-10,
            "mass-weighted displacement {momentum:?}"
        );
    }

    #[test]
    fn pinned_vertices_stay_and_the_rest_compensates() {
        let (points, faces) = icosphere(2);
        let rest = particles_from(&points, Vec3::zeros());
        let rest_volume = BalloonVolumeConstraint::new(faces.clone(), 1.0, 1.0, 0.0)
            .enclosed_volume(&rest)
            .expect("valid mesh");
        let mut particles = particles_from(&points, Vec3::zeros());
        let pinned: Vec<usize> = (0..particles.len()).filter(|k| k % 7 == 0).collect();
        for &k in &pinned {
            let position = particles[k].position;
            particles[k] = SoftParticle::new_static(position);
        }
        let mut constraint = BalloonVolumeConstraint::new(faces, rest_volume, 1.1, 0.0);
        run(&mut constraint, &mut particles, 60);
        for &k in &pinned {
            assert_eq!(
                particles[k].position, rest[k].position,
                "pinned vertex {k} moved"
            );
        }
        let ratio = constraint.enclosed_volume(&particles).expect("valid mesh") / rest_volume;
        assert!((ratio - 1.1).abs() < 0.011, "V / V0 = {ratio}");
    }

    #[test]
    fn compliance_softens_the_correction_and_lambda_resets() {
        let (points, faces) = icosphere(2);
        let rest = particles_from(&points, Vec3::zeros());
        let rest_volume = BalloonVolumeConstraint::new(faces.clone(), 1.0, 1.0, 0.0)
            .enclosed_volume(&rest)
            .expect("valid mesh");
        let mut stiff_particles = particles_from(&points, Vec3::zeros());
        let mut soft_particles = particles_from(&points, Vec3::zeros());
        let mut stiff = BalloonVolumeConstraint::new(faces.clone(), rest_volume, 1.3, 0.0);
        let mut soft = BalloonVolumeConstraint::new(faces, rest_volume, 1.3, 1.0e-2);
        stiff.project(&mut stiff_particles, 1.0 / 60.0);
        soft.project(&mut soft_particles, 1.0 / 60.0);
        let stiff_volume = stiff.enclosed_volume(&stiff_particles).expect("valid mesh");
        let soft_volume = soft.enclosed_volume(&soft_particles).expect("valid mesh");
        assert!(soft_volume > rest_volume, "the soft balloon still inflates");
        assert!(
            soft_volume < stiff_volume,
            "compliance must reduce the first correction: {soft_volume} vs {stiff_volume}"
        );
        assert!(soft.lambda() > 0.0 && stiff.lambda() > 0.0);
        soft.reset_lambda();
        assert_eq!(soft.lambda(), 0.0);
        let boxed: &mut dyn SoftConstraint = &mut stiff;
        boxed.reset_lambda();
        assert_eq!(stiff.lambda(), 0.0);
    }

    #[test]
    fn invalid_input_is_ignored_without_panicking() {
        let (points, faces) = icosphere(1);
        let mut particles = particles_from(&points, Vec3::zeros());
        let snapshot: Vec<Vec3> = particles.iter().map(|p| p.position).collect();

        let mut constraint = BalloonVolumeConstraint::new(faces.clone(), 1.0, 2.0, 0.0);
        for bad_dt in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            constraint.project(&mut particles, bad_dt);
        }
        for (p, s) in particles.iter().zip(&snapshot) {
            assert_eq!(p.position, *s);
        }

        let mut out_of_range = vec![[0, 1, 2], [0, 2, 999]];
        let mut broken = BalloonVolumeConstraint::new(out_of_range.clone(), 1.0, 2.0, 0.0);
        broken.project(&mut particles, 1.0 / 60.0);
        for (p, s) in particles.iter().zip(&snapshot) {
            assert_eq!(p.position, *s);
        }
        assert!(broken.enclosed_volume(&particles).is_none());
        assert!(
            BalloonVolumeConstraint::from_particles(out_of_range.clone(), &particles, 1.0, 0.0)
                .is_err()
        );
        out_of_range.clear();
        let mut empty = BalloonVolumeConstraint::new(out_of_range, 1.0, 2.0, 0.0);
        empty.project(&mut particles, 1.0 / 60.0);
        assert_eq!(empty.triangles().len(), 0);

        // Every vertex pinned: nothing can move.
        let mut frozen: Vec<SoftParticle> = particles
            .iter()
            .map(|p| SoftParticle::new_static(p.position))
            .collect();
        let mut constraint = BalloonVolumeConstraint::new(faces, 1.0, 2.0, 0.0);
        constraint.project(&mut frozen, 1.0 / 60.0);
        for (p, s) in frozen.iter().zip(&snapshot) {
            assert_eq!(p.position, *s);
        }
        assert_eq!(constraint.lambda(), 0.0);
    }

    #[test]
    fn elastic_membrane_balloon_inflates_through_the_xpbd_solver() {
        use crate::constraint::DistanceConstraint;
        use crate::particle::SoftBody;
        use crate::solver::XpbdSolver;
        use std::collections::BTreeSet;

        // A stretchable membrane (distance constraints along every edge) is
        // inflated to 120 % of its volume by the balloon constraint, driven
        // through the dyn-dispatching XPBD solver that restarts the multiplier
        // in every sub-step.
        let (points, faces) = icosphere(2);
        let centre = Vec3::new(4.0, -1.5, 2.0);
        let particles = particles_from(&points, centre);
        let volume_probe = BalloonVolumeConstraint::new(faces.clone(), 1.0, 1.0, 0.0);
        let rest_volume = volume_probe
            .enclosed_volume(&particles)
            .expect("valid mesh");
        let edges: BTreeSet<(usize, usize)> = faces
            .iter()
            .flat_map(|&[a, b, c]| [(a, b), (b, c), (c, a)])
            .map(|(u, v)| (u.min(v), u.max(v)))
            .collect();
        let mut constraints: Vec<Box<dyn SoftConstraint>> = edges
            .iter()
            .map(|&(u, v)| {
                Box::new(DistanceConstraint::from_particles(u, v, &particles, 1.0e-3))
                    as Box<dyn SoftConstraint>
            })
            .collect();
        constraints.push(Box::new(BalloonVolumeConstraint::new(
            faces,
            rest_volume,
            1.2,
            0.0,
        )));
        let mut body = SoftBody::from_particles(particles);
        body.damping = 0.02;
        let mut solver = XpbdSolver::new(10);
        // Keep the solver awake: the test checks the settled state, not sleeping.
        solver.sleep_threshold = 0.0;
        for _ in 0..120 {
            solver.solve(&mut body, &mut constraints, 1.0 / 60.0);
        }

        let volume = volume_probe
            .enclosed_volume(&body.particles)
            .expect("valid mesh");
        assert!(
            (volume / rest_volume - 1.2).abs() < 0.006,
            "V / V0 = {} should settle at 1.2",
            volume / rest_volume
        );
        let radii: Vec<f64> = body
            .particles
            .iter()
            .map(|p| (p.position - centre).norm())
            .collect();
        let mean_radius = radii.iter().sum::<f64>() / radii.len() as f64;
        let expected_radius = 1.2_f64.cbrt();
        assert!(
            (mean_radius - expected_radius).abs() < 0.01 * expected_radius,
            "mean radius {mean_radius} vs {expected_radius}"
        );
        let spread = radii
            .iter()
            .map(|r| (r - mean_radius).abs())
            .fold(0.0_f64, f64::max);
        assert!(
            spread < 0.02,
            "an inflated spherical balloon must stay spherical, radius spread {spread}"
        );
        let centroid = body
            .particles
            .iter()
            .fold(Vec3::zeros(), |acc, p| acc + p.position)
            / body.particles.len() as f64;
        assert!(
            (centroid - centre).norm() < 1e-9,
            "internal pressure must not move the centroid, moved by {}",
            (centroid - centre).norm()
        );
        assert!(
            body.particles
                .iter()
                .all(|p| p.position.iter().all(|c| c.is_finite())),
            "positions must stay finite"
        );
    }

    #[test]
    fn extreme_pressure_ratios_converge() {
        // The enclosed volume is cubic in the vertex positions, so a Newton
        // step overshoots for large ratios; the iteration must still settle.
        let (points, faces) = icosphere(2);
        let rest = particles_from(&points, Vec3::zeros());
        let rest_volume = BalloonVolumeConstraint::new(faces.clone(), 1.0, 1.0, 0.0)
            .enclosed_volume(&rest)
            .expect("valid mesh");
        for pressure in [0.1, 0.2, 3.0, 8.0] {
            let mut particles = particles_from(&points, Vec3::new(1.0, 2.0, 3.0));
            let mut constraint =
                BalloonVolumeConstraint::new(faces.clone(), rest_volume, pressure, 0.0);
            run(&mut constraint, &mut particles, 20);
            let ratio = constraint.enclosed_volume(&particles).expect("valid mesh") / rest_volume;
            assert!(
                (ratio / pressure - 1.0).abs() < 0.01,
                "V / V0 = {ratio} should settle at {pressure}"
            );
            assert!(
                particles
                    .iter()
                    .all(|p| p.position.iter().all(|c| c.is_finite())),
                "pressure {pressure}: positions must stay finite"
            );
        }
    }
}
