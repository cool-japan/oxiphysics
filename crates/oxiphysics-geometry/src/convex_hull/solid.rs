// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Polyhedral solid of a [`ConvexHull`](super::types::ConvexHull).
//!
//! A [`ConvexHull`](super::types::ConvexHull) is defined by a vertex cloud.
//! [`HullSolid`] turns that cloud into the closed triangulated polyhedron it
//! spans (outward counter-clockwise facets and their supporting planes) and
//! derives from it, by exact polyhedral integration, the volume, centre of
//! mass and inertia tensor of the uniform solid, as well as ray casts.

use oxiphysics_core::math::{Mat3, Real, Vec3};

use super::types::ConvexHull3D;
use crate::shape::RayHit;

/// Rays whose direction has a component along a face normal smaller than
/// this (relative to the direction length) are treated as parallel to it.
const PARALLEL_TOLERANCE: Real = 1e-12;

/// A supporting plane `normal · x = offset` with outward unit `normal`.
#[derive(Debug, Clone, Copy)]
struct Plane {
    normal: Vec3,
    offset: Real,
}

/// Closed triangulated polyhedron spanned by a vertex cloud, with the mass
/// properties of the uniform solid.
#[derive(Debug, Clone)]
pub(super) struct HullSolid {
    /// Supporting planes of the (non-degenerate) facets.
    planes: Vec<Plane>,
    /// Volume of the solid.
    volume: Real,
    /// Centre of mass of the uniform solid.
    center_of_mass: Vec3,
    /// Inertia tensor about the centre of mass for unit mass.
    unit_inertia: Mat3,
}

impl HullSolid {
    /// Build the solid spanned by `vertices`.
    ///
    /// Returns `None` if a coordinate is not finite or if the points do not
    /// span a volume (fewer than four points, or all points coplanar or
    /// collinear).
    pub(super) fn build(vertices: &[Vec3]) -> Option<Self> {
        let points: Vec<[f64; 3]> = vertices.iter().map(|v| [v.x, v.y, v.z]).collect();
        let hull = ConvexHull3D::build(&points)?;
        let hull_vertices: Vec<Vec3> = hull
            .vertices
            .iter()
            .map(|v| Vec3::new(v[0], v[1], v[2]))
            .collect();
        // The mean of the hull vertices is strictly interior for a hull with
        // volume; it is the reference apex of the facet tetrahedra.
        let reference = hull_vertices.iter().sum::<Vec3>() / hull_vertices.len() as Real;

        // Outward facets: the hull builder winds facets counter-clockwise
        // seen from outside; flip any facet that still faces the interior.
        let mut facets: Vec<[Vec3; 3]> = Vec::with_capacity(hull.faces.len());
        for face in &hull.faces {
            let (a, mut b, mut c) = (
                hull_vertices[face[0]],
                hull_vertices[face[1]],
                hull_vertices[face[2]],
            );
            let normal = (b - a).cross(&(c - a));
            if normal.dot(&(reference - a)) > 0.0 {
                std::mem::swap(&mut b, &mut c);
            }
            facets.push([a, b, c]);
        }

        // Volume and centre of mass from signed tetrahedra (reference, a, b, c).
        let mut six_volume = 0.0;
        let mut weighted_centroid = Vec3::zeros();
        for [a, b, c] in &facets {
            let (pa, pb, pc) = (a - reference, b - reference, c - reference);
            let det = pa.dot(&pb.cross(&pc));
            six_volume += det;
            weighted_centroid += (pa + pb + pc) * (det / 4.0);
        }
        let volume = six_volume / 6.0;
        if !(volume.is_finite() && volume > 0.0) {
            return None;
        }
        let center_of_mass = reference + weighted_centroid / six_volume;

        // Second moments C = ∫ x xᵀ dV about the centre of mass.  For the
        // tetrahedron with apex at the origin and base (a, b, c):
        //   ∫ x xᵀ dV = V/20 (a aᵀ + b bᵀ + c cᵀ + s sᵀ),  s = a + b + c.
        let mut second_moment = Mat3::zeros();
        let mut planes: Vec<Plane> = Vec::with_capacity(facets.len());
        for [a, b, c] in &facets {
            let (pa, pb, pc) = (a - center_of_mass, b - center_of_mass, c - center_of_mass);
            let tet_volume = pa.dot(&pb.cross(&pc)) / 6.0;
            let sum = pa + pb + pc;
            second_moment += (pa * pa.transpose()
                + pb * pb.transpose()
                + pc * pc.transpose()
                + sum * sum.transpose())
                * (tet_volume / 20.0);

            let normal = (b - a).cross(&(c - a));
            let length = normal.norm();
            if length > 0.0 && length.is_finite() {
                let unit = normal / length;
                planes.push(Plane {
                    normal: unit,
                    offset: unit.dot(a),
                });
            }
        }
        // I = ρ (tr(C) 1 − C) with ρ = m / V; symmetrise away rounding noise.
        let inertia_unit_density = Mat3::identity() * second_moment.trace() - second_moment;
        let unit_inertia =
            (inertia_unit_density + inertia_unit_density.transpose()) * (0.5 / volume);

        Some(Self {
            planes,
            volume,
            center_of_mass,
            unit_inertia,
        })
    }

    /// Volume of the solid.
    pub(super) fn volume(&self) -> Real {
        self.volume
    }

    /// Centre of mass of the uniform solid.
    pub(super) fn center_of_mass(&self) -> Vec3 {
        self.center_of_mass
    }

    /// Inertia tensor about the centre of mass for a body of the given mass.
    pub(super) fn inertia_tensor(&self, mass: Real) -> Mat3 {
        self.unit_inertia * mass
    }

    /// Cast a ray against the solid (Cyrus–Beck clipping against the facet
    /// planes).
    ///
    /// `toi` is measured in units of `direction` (not normalised).  A ray
    /// that starts inside the solid reports the point where it leaves, with
    /// the outward normal there.
    pub(super) fn ray_cast(
        &self,
        origin: &Vec3,
        direction: &Vec3,
        max_toi: Real,
    ) -> Option<RayHit> {
        let direction_length = direction.norm();
        if !(direction_length > 0.0 && direction_length.is_finite()) {
            return None;
        }
        let mut t_enter = Real::NEG_INFINITY;
        let mut t_exit = Real::INFINITY;
        let mut enter_normal: Option<Vec3> = None;
        let mut exit_normal: Option<Vec3> = None;
        for plane in &self.planes {
            let slope = plane.normal.dot(direction);
            // Distance of the origin below the plane (>= 0 on the inner side).
            let gap = plane.offset - plane.normal.dot(origin);
            if slope.abs() <= PARALLEL_TOLERANCE * direction_length {
                if gap < 0.0 {
                    return None;
                }
                continue;
            }
            let t = gap / slope;
            if slope < 0.0 {
                if t > t_enter {
                    t_enter = t;
                    enter_normal = Some(plane.normal);
                }
            } else if t < t_exit {
                t_exit = t;
                exit_normal = Some(plane.normal);
            }
            if t_enter > t_exit {
                return None;
            }
        }
        let (toi, normal) = if t_enter >= 0.0 {
            (t_enter, enter_normal?)
        } else if t_exit >= 0.0 {
            (t_exit, exit_normal?)
        } else {
            return None;
        };
        if toi > max_toi {
            return None;
        }
        Some(RayHit {
            point: origin + direction * toi,
            normal,
            toi,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::FRAC_PI_4;

    fn cube_vertices(lo: f64, hi: f64) -> Vec<Vec3> {
        let mut v = Vec::new();
        for &x in &[lo, hi] {
            for &y in &[lo, hi] {
                for &z in &[lo, hi] {
                    v.push(Vec3::new(x, y, z));
                }
            }
        }
        v
    }

    fn tetra_vertices() -> Vec<Vec3> {
        vec![
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
            Vec3::new(0.0, 0.0, 1.0),
        ]
    }

    fn assert_mat_close(got: &Mat3, want: &Mat3, tol: f64) {
        for i in 0..3 {
            for j in 0..3 {
                assert!(
                    (got[(i, j)] - want[(i, j)]).abs() <= tol,
                    "[{i}][{j}] = {} expected {}",
                    got[(i, j)],
                    want[(i, j)]
                );
            }
        }
    }

    /// Deterministic permutation generator for vertex-order tests.
    fn shuffled<T: Clone>(items: &[T], seed: u64) -> Vec<T> {
        let mut state = seed | 1;
        let mut out = items.to_vec();
        for i in (1..out.len()).rev() {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            out.swap(i, (state % (i as u64 + 1)) as usize);
        }
        out
    }

    #[test]
    fn test_unit_cube_mass_properties_for_every_vertex_order() {
        for seed in 0..40 {
            let solid = HullSolid::build(&shuffled(&cube_vertices(0.0, 1.0), seed))
                .expect("a cube spans a volume");
            assert!((solid.volume() - 1.0).abs() < 1e-12, "seed {seed}");
            let com = solid.center_of_mass();
            assert!(
                (com - Vec3::new(0.5, 0.5, 0.5)).norm() < 1e-12,
                "seed {seed}"
            );
            // Unit-mass cube: I = m a² / 6 on the diagonal, zero products.
            let want = Mat3::identity() * (1.0 / 6.0);
            assert_mat_close(&solid.inertia_tensor(1.0), &want, 1e-12);
            assert_mat_close(&solid.inertia_tensor(6.0), &Mat3::identity(), 1e-11);
        }
    }

    #[test]
    fn test_unit_tetrahedron_mass_properties() {
        let solid = HullSolid::build(&tetra_vertices()).expect("tetrahedron");
        assert!((solid.volume() - 1.0 / 6.0).abs() < 1e-14);
        assert!((solid.center_of_mass() - Vec3::new(0.25, 0.25, 0.25)).norm() < 1e-14);
        // About the centre of mass: I_xx = 3m/40 and I_xy = m/80.
        let mass = 2.0;
        let inertia = solid.inertia_tensor(mass);
        for i in 0..3 {
            assert!((inertia[(i, i)] - 3.0 * mass / 40.0).abs() < 1e-14);
            for j in 0..3 {
                if i != j {
                    assert!((inertia[(i, j)] - mass / 80.0).abs() < 1e-14);
                }
            }
        }
    }

    #[test]
    fn test_box_hull_matches_closed_form_box_inertia() {
        // Half-extents (1, 2, 3) centred at an offset: inertia about the
        // centre of mass is the textbook box tensor.
        let (hx, hy, hz) = (1.0, 2.0, 3.0);
        let offset = Vec3::new(5.0, -3.0, 2.0);
        let mut vertices = Vec::new();
        for &sx in &[-1.0, 1.0] {
            for &sy in &[-1.0, 1.0] {
                for &sz in &[-1.0, 1.0] {
                    vertices.push(offset + Vec3::new(sx * hx, sy * hy, sz * hz));
                }
            }
        }
        let solid = HullSolid::build(&vertices).expect("box");
        assert!((solid.volume() - 8.0 * hx * hy * hz).abs() < 1e-10);
        assert!((solid.center_of_mass() - offset).norm() < 1e-12);
        let mass = 12.0;
        let k = mass / 12.0;
        let want = Mat3::new(
            k * (4.0 * hy * hy + 4.0 * hz * hz),
            0.0,
            0.0,
            0.0,
            k * (4.0 * hx * hx + 4.0 * hz * hz),
            0.0,
            0.0,
            0.0,
            k * (4.0 * hx * hx + 4.0 * hy * hy),
        );
        assert_mat_close(&solid.inertia_tensor(mass), &want, 1e-9);
    }

    #[test]
    fn test_rotated_box_hull_has_rotated_inertia_with_cross_terms() {
        let (hx, hy, hz) = (1.0, 2.0, 3.0);
        let rotation =
            oxiphysics_core::math::quat_from_axis_angle(&Vec3::new(0.0, 0.0, 1.0), FRAC_PI_4)
                .to_rotation_matrix()
                .into_inner();
        let mut vertices = Vec::new();
        for &sx in &[-1.0, 1.0] {
            for &sy in &[-1.0, 1.0] {
                for &sz in &[-1.0, 1.0] {
                    vertices.push(rotation * Vec3::new(sx * hx, sy * hy, sz * hz));
                }
            }
        }
        let solid = HullSolid::build(&vertices).expect("rotated box");
        let mass = 7.0;
        let k = mass / 12.0;
        let local = Mat3::new(
            k * (4.0 * hy * hy + 4.0 * hz * hz),
            0.0,
            0.0,
            0.0,
            k * (4.0 * hx * hx + 4.0 * hz * hz),
            0.0,
            0.0,
            0.0,
            k * (4.0 * hx * hx + 4.0 * hy * hy),
        );
        let want = rotation * local * rotation.transpose();
        let got = solid.inertia_tensor(mass);
        assert_mat_close(&got, &want, 1e-9);
        // The 45° rotation mixes x and y: the product of inertia is
        // (I_xx - I_yy) / 2 = k (4 hy² - 4 hx²) / 2 ... in magnitude.
        let expected_xy = 0.5 * (local[(0, 0)] - local[(1, 1)]).abs();
        assert!(expected_xy > 1.0, "test box must be anisotropic");
        assert!(
            (got[(0, 1)].abs() - expected_xy).abs() < 1e-9,
            "{}",
            got[(0, 1)]
        );
    }

    #[test]
    fn test_interior_and_duplicate_points_do_not_change_the_solid() {
        let mut vertices = cube_vertices(0.0, 2.0);
        vertices.push(Vec3::new(1.0, 1.0, 1.0));
        vertices.push(Vec3::new(0.3, 1.7, 0.9));
        vertices.push(Vec3::new(2.0, 2.0, 2.0));
        vertices.push(Vec3::new(1.0, 1.0, 0.0)); // centre of a face
        let solid = HullSolid::build(&vertices).expect("cube with extras");
        assert!((solid.volume() - 8.0).abs() < 1e-10);
        assert!((solid.center_of_mass() - Vec3::new(1.0, 1.0, 1.0)).norm() < 1e-10);
    }

    #[test]
    fn test_degenerate_or_non_finite_clouds_have_no_solid() {
        assert!(HullSolid::build(&[]).is_none());
        let flat: Vec<Vec3> = cube_vertices(0.0, 1.0)
            .into_iter()
            .map(|v| Vec3::new(v.x, v.y, 0.25))
            .collect();
        assert!(HullSolid::build(&flat).is_none());
        let line: Vec<Vec3> = (0..6)
            .map(|i| Vec3::new(i as f64, 2.0 * i as f64, 0.0))
            .collect();
        assert!(HullSolid::build(&line).is_none());
        let mut bad = cube_vertices(0.0, 1.0);
        bad[3].y = f64::NAN;
        assert!(HullSolid::build(&bad).is_none());
        bad[3].y = f64::INFINITY;
        assert!(HullSolid::build(&bad).is_none());
    }

    #[test]
    fn test_small_and_large_scales() {
        for &scale in &[1e-4, 1.0, 1e5] {
            let vertices: Vec<Vec3> = cube_vertices(0.0, 1.0)
                .into_iter()
                .map(|v| v * scale)
                .collect();
            let solid = HullSolid::build(&vertices).expect("scaled cube");
            let want = scale * scale * scale;
            assert!(
                (solid.volume() - want).abs() <= 1e-10 * want,
                "scale {scale}"
            );
        }
    }

    #[test]
    fn test_ray_cast_cube_axis_ray_hits_front_face_at_minus_half() {
        let solid = HullSolid::build(&cube_vertices(-0.5, 0.5)).expect("cube");
        let hit = solid
            .ray_cast(&Vec3::new(-3.0, 0.0, 0.0), &Vec3::new(1.0, 0.0, 0.0), 100.0)
            .expect("hit");
        // The ray starts at x = -3, so the front face x = -0.5 is 2.5 away.
        assert!((hit.toi - 2.5).abs() < 1e-12, "toi = {}", hit.toi);
        assert!((hit.point - Vec3::new(-0.5, 0.0, 0.0)).norm() < 1e-12);
        assert!((hit.normal - Vec3::new(-1.0, 0.0, 0.0)).norm() < 1e-12);
        // The direction is not normalised: doubling it halves the toi.
        let hit = solid
            .ray_cast(&Vec3::new(-3.0, 0.0, 0.0), &Vec3::new(2.0, 0.0, 0.0), 100.0)
            .expect("hit");
        assert!((hit.toi - 1.25).abs() < 1e-12);
        assert!((hit.point - Vec3::new(-0.5, 0.0, 0.0)).norm() < 1e-12);
    }

    #[test]
    fn test_ray_cast_misses_and_limits() {
        let solid = HullSolid::build(&cube_vertices(-0.5, 0.5)).expect("cube");
        let dir = Vec3::new(1.0, 0.0, 0.0);
        // Passes above the cube.
        assert!(
            solid
                .ray_cast(&Vec3::new(-3.0, 0.6, 0.0), &dir, 100.0)
                .is_none()
        );
        // Pointing away.
        assert!(
            solid
                .ray_cast(&Vec3::new(3.0, 0.0, 0.0), &dir, 100.0)
                .is_none()
        );
        // Hit farther than max_toi.
        assert!(
            solid
                .ray_cast(&Vec3::new(-3.0, 0.0, 0.0), &dir, 2.0)
                .is_none()
        );
        assert!(
            solid
                .ray_cast(&Vec3::new(-3.0, 0.0, 0.0), &dir, 2.5)
                .is_some()
        );
        // Zero or non-finite direction.
        assert!(
            solid
                .ray_cast(&Vec3::new(-3.0, 0.0, 0.0), &Vec3::zeros(), 10.0)
                .is_none()
        );
        assert!(
            solid
                .ray_cast(
                    &Vec3::new(-3.0, 0.0, 0.0),
                    &Vec3::new(f64::NAN, 0.0, 0.0),
                    10.0
                )
                .is_none()
        );
        // Parallel to a face pair and outside its slab.
        assert!(
            solid
                .ray_cast(&Vec3::new(-3.0, 0.0, 0.7), &dir, 100.0)
                .is_none()
        );
    }

    #[test]
    fn test_ray_cast_from_inside_reports_exit_with_outward_normal() {
        let solid = HullSolid::build(&cube_vertices(-0.5, 0.5)).expect("cube");
        let hit = solid
            .ray_cast(&Vec3::new(0.1, 0.0, 0.0), &Vec3::new(0.0, 1.0, 0.0), 10.0)
            .expect("exit hit");
        assert!((hit.toi - 0.5).abs() < 1e-12);
        assert!((hit.normal - Vec3::new(0.0, 1.0, 0.0)).norm() < 1e-12);
    }

    #[test]
    fn test_ray_cast_oblique_ray_on_tetrahedron() {
        let solid = HullSolid::build(&tetra_vertices()).expect("tetrahedron");
        // Through the slanted face x + y + z = 1 at (1/3, 1/3, 1/3).
        let hit = solid
            .ray_cast(
                &Vec3::new(2.0, 2.0, 2.0),
                &Vec3::new(-1.0, -1.0, -1.0),
                10.0,
            )
            .expect("hit");
        assert!((hit.point - Vec3::new(1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0)).norm() < 1e-12);
        assert!((hit.toi - 5.0 / 3.0).abs() < 1e-12);
        let n = 1.0 / 3.0_f64.sqrt();
        assert!((hit.normal - Vec3::new(n, n, n)).norm() < 1e-12);
    }
}
