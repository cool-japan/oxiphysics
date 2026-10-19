// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Mass distribution of a triangle mesh: centre of mass and inertia tensor.
//!
//! For a closed, consistently oriented mesh the enclosed solid is integrated
//! exactly with the divergence theorem.  Every triangle `(a, b, c)` forms a
//! signed tetrahedron with a reference point `r`; with `a`, `b`, `c` taken
//! relative to `r` the tetrahedron has
//!
//! * volume `det / 6` with `det = a . (b x c)`,
//! * first moment `det / 24 * (a + b + c)`,
//! * second moment `det / 120 * (a a^T + b b^T + c c^T + s s^T)` with
//!   `s = a + b + c`.
//!
//! Signed tetrahedra outside the solid cancel, so for a closed surface the sums
//! are exact for any reference point.  The centroid is found first and the
//! second moment is summed relative to it, which keeps the result accurate for
//! meshes far from the origin.
//!
//! The sums depend on the reference point when the surface is open: moving it
//! by `d` changes the signed volume by `-d . S / 3`, where
//! `S = 1/2 sum (b - a) x (c - a)` is the vector area of the triangles, which
//! vanishes exactly when the surface is closed.  A mesh is taken to enclose a
//! volume when its signed volume is well above rounding noise and moving the
//! reference point across the whole extent of the mesh would change it by at
//! most [`MAX_RELATIVE_VOLUME_UNCERTAINTY`].  A closed mesh qualifies, with or
//! without shared vertices, and so does a closed mesh with a few missing
//! triangles.
//!
//! A mesh that encloses no volume (an open surface such as a flat plate or a
//! box without a lid, or a surface of zero volume) has no solid to integrate.
//! Its mass is then spread uniformly over the surface (a thin shell),
//! integrated exactly over the triangles with the second moment
//! `A / 12 * (a a^T + b b^T + c c^T + s s^T)`.  The inertia tensor and the
//! centre of mass always describe the same distribution, so the tensor is the
//! one about the reported centre of mass.

use super::TriangleMesh;
use oxiphysics_core::math::{Mat3, Real, Vec3};

/// A mesh whose absolute signed volume is below this fraction of the cubed
/// bounding-box diagonal encloses no volume.
const RELATIVE_VOLUME_TOLERANCE: f64 = 1e-12;

/// Largest change of the signed volume, relative to the volume itself, that
/// moving the reference point across the extent of the mesh may cause for the
/// mesh to still count as the boundary of a solid.
const MAX_RELATIVE_VOLUME_UNCERTAINTY: f64 = 1e-2;

/// How the mass of a triangle mesh is distributed.
#[derive(Debug, Clone)]
pub(super) enum MassDistribution {
    /// The mass fills the volume enclosed by the mesh.
    Solid {
        /// Signed volume (negative when the triangles face inwards).
        signed_volume: Real,
        /// Centre of mass.
        centroid: Vec3,
        /// `integral (x - centroid)(x - centroid)^T dV`, scaled by the sign of
        /// the orientation like `signed_volume`.
        second_moment: Mat3,
    },
    /// The mesh encloses no volume; the mass is spread over its surface.
    Shell {
        /// Total surface area.
        area: Real,
        /// Area centroid.
        centroid: Vec3,
        /// `integral (x - centroid)(x - centroid)^T dA`.
        second_moment: Mat3,
    },
    /// No usable triangle.
    Empty,
}

impl MassDistribution {
    /// Centre of mass (the origin for [`MassDistribution::Empty`]).
    pub(super) fn centroid(&self) -> Vec3 {
        match self {
            Self::Solid { centroid, .. } | Self::Shell { centroid, .. } => *centroid,
            Self::Empty => Vec3::zeros(),
        }
    }

    /// Inertia tensor about the centre of mass for the given total `mass`.
    ///
    /// The zero matrix for a non-positive or non-finite mass and for an empty
    /// distribution.
    pub(super) fn inertia_tensor(&self, mass: Real) -> Mat3 {
        if !(mass.is_finite() && mass > 0.0) {
            return Mat3::zeros();
        }
        let (measure, second_moment) = match self {
            Self::Solid {
                signed_volume,
                second_moment,
                ..
            } => (*signed_volume, second_moment),
            Self::Shell {
                area,
                second_moment,
                ..
            } => (*area, second_moment),
            Self::Empty => return Mat3::zeros(),
        };
        // `measure` carries the same orientation sign as `second_moment`, so
        // the density `mass / measure` makes the tensor positive either way.
        let density = mass / measure;
        let tensor = (Mat3::identity() * second_moment.trace() - second_moment) * density;
        // Remove rounding asymmetry.
        (tensor + tensor.transpose()) * 0.5
    }
}

/// `a a^T + b b^T + c c^T + s s^T` with `s = a + b + c`.
fn corner_second_moment(a: &Vec3, b: &Vec3, c: &Vec3) -> Mat3 {
    let s = a + b + c;
    a * a.transpose() + b * b.transpose() + c * c.transpose() + s * s.transpose()
}

impl TriangleMesh {
    /// The triangle corners of every triangle whose indices are in range.
    fn corner_triples(&self) -> impl Iterator<Item = [Vec3; 3]> + '_ {
        self.indices.iter().filter_map(|tri| {
            Some([
                *self.vertices.get(tri[0])?,
                *self.vertices.get(tri[1])?,
                *self.vertices.get(tri[2])?,
            ])
        })
    }

    /// Integrate the mass distribution of the mesh; see the module docs.
    ///
    /// Triangles that index past the vertex list are ignored.
    pub(super) fn mass_distribution(&self) -> MassDistribution {
        let mut lower = Vec3::repeat(f64::INFINITY);
        let mut upper = Vec3::repeat(f64::NEG_INFINITY);
        let mut triangle_count = 0usize;
        for corners in self.corner_triples() {
            triangle_count += 1;
            for corner in &corners {
                lower = lower.inf(corner);
                upper = upper.sup(corner);
            }
        }
        if triangle_count == 0 {
            return MassDistribution::Empty;
        }
        let reference = (lower + upper) * 0.5;
        let diagonal = (upper - lower).norm();

        // Pass 1: signed volume, volume first moment, area, area first moment
        // and vector area.
        let mut six_volume = 0.0;
        let mut volume_moment = Vec3::zeros();
        let mut area = 0.0;
        let mut area_moment = Vec3::zeros();
        let mut vector_area = Vec3::zeros();
        for corners in self.corner_triples() {
            let a = corners[0] - reference;
            let b = corners[1] - reference;
            let c = corners[2] - reference;
            let det = a.dot(&b.cross(&c));
            six_volume += det;
            volume_moment += (a + b + c) * (det / 24.0);
            let normal = (b - a).cross(&(c - a));
            vector_area += normal * 0.5;
            let triangle_area = 0.5 * normal.norm();
            area += triangle_area;
            area_moment += (a + b + c) * (triangle_area / 3.0);
        }
        let signed_volume = six_volume / 6.0;

        // The change of the signed volume when the reference point moves across
        // the mesh; zero for a closed surface.
        let volume_uncertainty = vector_area.norm() * diagonal / 3.0;
        let encloses_volume = signed_volume.abs() > RELATIVE_VOLUME_TOLERANCE * diagonal.powi(3)
            && volume_uncertainty <= MAX_RELATIVE_VOLUME_UNCERTAINTY * signed_volume.abs();
        if encloses_volume {
            let centroid = reference + volume_moment / signed_volume;
            let mut second_moment = Mat3::zeros();
            for corners in self.corner_triples() {
                let a = corners[0] - centroid;
                let b = corners[1] - centroid;
                let c = corners[2] - centroid;
                let det = a.dot(&b.cross(&c));
                second_moment += corner_second_moment(&a, &b, &c) * (det / 120.0);
            }
            return MassDistribution::Solid {
                signed_volume,
                centroid,
                second_moment,
            };
        }

        if area > 0.0 && area.is_finite() {
            let centroid = reference + area_moment / area;
            let mut second_moment = Mat3::zeros();
            for corners in self.corner_triples() {
                let a = corners[0] - centroid;
                let b = corners[1] - centroid;
                let c = corners[2] - centroid;
                let triangle_area = 0.5 * (b - a).cross(&(c - a)).norm();
                second_moment += corner_second_moment(&a, &b, &c) * (triangle_area / 12.0);
            }
            return MassDistribution::Shell {
                area,
                centroid,
                second_moment,
            };
        }
        MassDistribution::Empty
    }
}

#[cfg(test)]
mod tests {
    use super::super::{tetrahedron_mesh, unit_cube_mesh};
    use super::*;
    use crate::shape::Shape;
    use oxiphysics_core::math::{Mat3, Vec3};

    /// Axis-aligned box `[0, sx] x [0, sy] x [0, sz]` with outward triangles.
    fn box_mesh(sx: f64, sy: f64, sz: f64) -> TriangleMesh {
        let mut mesh = unit_cube_mesh();
        for v in &mut mesh.vertices {
            *v = Vec3::new(v.x * sx, v.y * sy, v.z * sz);
        }
        mesh
    }

    fn transformed(mesh: &TriangleMesh, rotation: &Mat3, shift: Vec3) -> TriangleMesh {
        TriangleMesh::new(
            mesh.vertices.iter().map(|v| rotation * v + shift).collect(),
            mesh.indices.clone(),
        )
    }

    fn rotation_z(angle: f64) -> Mat3 {
        let (s, c) = angle.sin_cos();
        Mat3::new(c, -s, 0.0, s, c, 0.0, 0.0, 0.0, 1.0)
    }

    fn assert_matrix_close(actual: &Mat3, expected: &Mat3, tolerance: f64) {
        for i in 0..3 {
            for j in 0..3 {
                assert!(
                    (actual[(i, j)] - expected[(i, j)]).abs() <= tolerance,
                    "entry ({i},{j}): got {}, expected {}\nactual: {actual}\nexpected: {expected}",
                    actual[(i, j)],
                    expected[(i, j)]
                );
            }
        }
    }

    /// Inertia tensor of a solid box with side lengths `(sx, sy, sz)`.
    fn box_inertia(mass: f64, sx: f64, sy: f64, sz: f64) -> Mat3 {
        let k = mass / 12.0;
        Mat3::new(
            k * (sy * sy + sz * sz),
            0.0,
            0.0,
            0.0,
            k * (sx * sx + sz * sz),
            0.0,
            0.0,
            0.0,
            k * (sx * sx + sy * sy),
        )
    }

    #[test]
    fn unit_cube_inertia_is_mass_over_six_on_the_diagonal() {
        let mesh = unit_cube_mesh();
        let mass = 3.0;
        let inertia = mesh.inertia_tensor(mass);
        assert_matrix_close(&inertia, &(Mat3::identity() * (mass / 6.0)), 1e-14);
        assert!((mesh.center_of_mass() - Vec3::new(0.5, 0.5, 0.5)).norm() < 1e-14);
    }

    #[test]
    fn offset_cube_inertia_is_still_about_its_centre_of_mass() {
        let mass = 2.0;
        let shift = Vec3::new(10.0, -7.0, 3.0);
        let mesh = transformed(&unit_cube_mesh(), &Mat3::identity(), shift);
        let inertia = mesh.inertia_tensor(mass);
        assert_matrix_close(&inertia, &(Mat3::identity() * (mass / 6.0)), 1e-12);
        let expected_com = shift + Vec3::new(0.5, 0.5, 0.5);
        assert!((mesh.center_of_mass() - expected_com).norm() < 1e-12);
    }

    #[test]
    fn box_inertia_is_anisotropic_and_not_the_bounding_box_of_a_different_shape() {
        let mesh = box_mesh(2.0, 1.0, 0.5);
        let mass = 5.0;
        assert_matrix_close(
            &mesh.inertia_tensor(mass),
            &box_inertia(mass, 2.0, 1.0, 0.5),
            1e-13,
        );
    }

    #[test]
    fn rotated_box_inertia_has_products_of_inertia() {
        let mass = 5.0;
        let rotation = rotation_z(30.0_f64.to_radians());
        let mesh = transformed(
            &box_mesh(2.0, 1.0, 0.5),
            &rotation,
            Vec3::new(1.0, 2.0, 3.0),
        );
        let expected = rotation * box_inertia(mass, 2.0, 1.0, 0.5) * rotation.transpose();
        let inertia = mesh.inertia_tensor(mass);
        assert_matrix_close(&inertia, &expected, 1e-12);
        // Ixx - Iyy = (m/12)(sy^2 - sx^2) = -(5/12) * 3 along the box axes, so
        // the rotated tensor has a non-zero product of inertia.
        assert!(inertia[(0, 1)].abs() > 0.3, "got Ixy = {}", inertia[(0, 1)]);
        assert!((inertia - inertia.transpose()).norm() < 1e-14);
        assert!(inertia[(0, 2)].abs() < 1e-12 && inertia[(1, 2)].abs() < 1e-12);
    }

    #[test]
    fn inward_oriented_mesh_gives_the_same_positive_inertia() {
        let mass = 4.0;
        let outward = box_mesh(1.0, 2.0, 3.0);
        let inward = TriangleMesh::new(
            outward.vertices.clone(),
            outward.indices.iter().map(|t| [t[0], t[2], t[1]]).collect(),
        );
        let expected = box_inertia(mass, 1.0, 2.0, 3.0);
        assert_matrix_close(&inward.inertia_tensor(mass), &expected, 1e-13);
        assert!((inward.center_of_mass() - outward.center_of_mass()).norm() < 1e-14);
    }

    #[test]
    fn regular_tetrahedron_inertia_is_mass_edge_squared_over_twenty() {
        // Edge length 2 * sqrt(2); I = m a^2 / 20 = 0.4 m about any axis
        // through the centroid.
        let mesh = tetrahedron_mesh();
        let mass = 1.5;
        let expected = Mat3::identity() * (0.4 * mass);
        assert_matrix_close(&mesh.inertia_tensor(mass), &expected, 1e-13);
        assert!(mesh.center_of_mass().norm() < 1e-14);
    }

    #[test]
    fn mass_properties_use_one_consistent_centre_of_mass() {
        let mesh = transformed(
            &box_mesh(2.0, 2.0, 2.0),
            &Mat3::identity(),
            Vec3::new(1.0, 1.0, 1.0),
        );
        let props = mesh.mass_properties(3.0);
        // density 3 * volume 8 = 24; I = m (a^2 + a^2) / 12 = 24 * 8 / 12 = 16
        // about the COM (2, 2, 2).
        assert!((props.mass - 24.0).abs() < 1e-12);
        assert!((props.center_of_mass - Vec3::new(2.0, 2.0, 2.0)).norm() < 1e-12);
        assert_matrix_close(&props.local_inertia, &(Mat3::identity() * 16.0), 1e-10);
    }

    #[test]
    fn tiny_mesh_keeps_its_centre_of_mass_and_inertia() {
        // A cube of side 1e-5 has volume 1e-15, far below an absolute 1e-12
        // threshold, yet is a perfectly well-defined solid.
        let side = 1.0e-5;
        let mesh = box_mesh(side, side, side);
        let mass = 1.0e-3;
        let expected = Mat3::identity() * (mass * side * side / 6.0);
        assert_matrix_close(
            &mesh.inertia_tensor(mass),
            &expected,
            1e-9 * expected[(0, 0)],
        );
        let half = 0.5 * side;
        assert!((mesh.center_of_mass() - Vec3::new(half, half, half)).norm() < 1e-18);
    }

    #[test]
    fn flat_plate_uses_the_thin_shell_inertia_about_its_area_centroid() {
        // Unit square in the plane y = 0.5 (two triangles); no enclosed volume.
        let mesh = TriangleMesh::new(
            vec![
                Vec3::new(0.0, 0.5, 0.0),
                Vec3::new(1.0, 0.5, 0.0),
                Vec3::new(1.0, 0.5, 1.0),
                Vec3::new(0.0, 0.5, 1.0),
            ],
            vec![[0, 1, 2], [0, 2, 3]],
        );
        let mass = 2.0;
        let com = mesh.center_of_mass();
        assert!((com - Vec3::new(0.5, 0.5, 0.5)).norm() < 1e-14, "COM {com}");
        // Lamina: Ixx = Izz = m / 12, Iyy = m / 6 (perpendicular axis theorem).
        let expected = Mat3::new(
            mass / 12.0,
            0.0,
            0.0,
            0.0,
            mass / 6.0,
            0.0,
            0.0,
            0.0,
            mass / 12.0,
        );
        assert_matrix_close(&mesh.inertia_tensor(mass), &expected, 1e-14);
    }

    /// Triangle soup of the unit cube `[0, 1]^3` whose faces are `n x n` grids
    /// of quads (two outward triangles each) that share no vertices.
    fn subdivided_cube_soup(n: usize) -> TriangleMesh {
        // (origin, u, v) of every face; `u x v` points out of the cube.
        let faces = [
            ([0.0, 0.0, 0.0], [0.0, 1.0, 0.0], [1.0, 0.0, 0.0]),
            ([0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]),
            ([0.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, 1.0, 0.0]),
            ([1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]),
            ([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]),
            ([0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]),
        ];
        let step = 1.0 / n as f64;
        let mut vertices = Vec::new();
        let mut indices = Vec::new();
        for (origin, u, v) in faces {
            let (origin, u, v) = (Vec3::from(origin), Vec3::from(u), Vec3::from(v));
            for i in 0..n {
                for j in 0..n {
                    let corner = |di: usize, dj: usize| {
                        origin + u * ((i + di) as f64 * step) + v * ((j + dj) as f64 * step)
                    };
                    let base = vertices.len();
                    vertices.extend([corner(0, 0), corner(1, 0), corner(1, 1), corner(0, 1)]);
                    indices.push([base, base + 1, base + 2]);
                    indices.push([base, base + 2, base + 3]);
                }
            }
        }
        TriangleMesh::new(vertices, indices)
    }

    #[test]
    fn closed_triangle_soup_without_shared_vertices_is_a_solid() {
        let mesh = subdivided_cube_soup(4);
        assert_eq!(mesh.indices.len(), 192);
        let mass = 2.4;
        assert_matrix_close(
            &mesh.inertia_tensor(mass),
            &(Mat3::identity() * (mass / 6.0)),
            1e-13,
        );
        assert!((mesh.center_of_mass() - Vec3::new(0.5, 0.5, 0.5)).norm() < 1e-14);
    }

    #[test]
    fn closed_mesh_with_a_missing_triangle_is_still_a_solid() {
        // One hole of area 0.005 in a mesh of 1200 triangles: moving the
        // reference point across the mesh changes the signed volume by far
        // less than 1 %, so the mesh is integrated as a solid and its inertia
        // stays close to that of the closed cube.
        let mut mesh = subdivided_cube_soup(10);
        assert_eq!(mesh.indices.len(), 1200);
        mesh.indices.remove(0);
        let mass = 1.0;
        assert_matrix_close(
            &mesh.inertia_tensor(mass),
            &(Mat3::identity() * (mass / 6.0)),
            0.01 * mass / 6.0,
        );
        let com = mesh.center_of_mass();
        assert!((com - Vec3::new(0.5, 0.5, 0.5)).norm() < 0.01, "COM {com}");
    }

    #[test]
    fn open_box_without_lid_is_a_thin_shell() {
        // Drop the +Z face (triangles 2 and 3) of the unit cube.  The five
        // remaining faces enclose no volume: the mass sits on the surface and
        // the tensor is the sum of five square laminae moved to the common
        // centre of mass.
        let mut mesh = unit_cube_mesh();
        mesh.indices.drain(2..4);
        assert_eq!(mesh.indices.len(), 10);
        let mass = 5.0;
        let com = mesh.center_of_mass();
        assert!((com - Vec3::new(0.5, 0.5, 0.4)).norm() < 1e-14, "COM {com}");
        // (face centre, axis of the face normal)
        let faces = [
            (Vec3::new(0.5, 0.5, 0.0), 2),
            (Vec3::new(0.0, 0.5, 0.5), 0),
            (Vec3::new(1.0, 0.5, 0.5), 0),
            (Vec3::new(0.5, 0.0, 0.5), 1),
            (Vec3::new(0.5, 1.0, 0.5), 1),
        ];
        let face_mass = mass / 5.0;
        let mut expected = Mat3::zeros();
        for (centre, normal_axis) in faces {
            // Unit square lamina: m / 12 about the in-plane axes, m / 6 about
            // the normal.
            let mut lamina = Mat3::identity() * (face_mass / 12.0);
            lamina[(normal_axis, normal_axis)] = face_mass / 6.0;
            let d = centre - com;
            let transfer = (Mat3::identity() * d.norm_squared() - d * d.transpose()) * face_mass;
            expected += lamina + transfer;
        }
        let inertia = mesh.inertia_tensor(mass);
        assert_matrix_close(&inertia, &expected, 1e-13);
        // The solid cube would have mass / 6 on the diagonal.
        assert!(inertia[(2, 2)] > 1.5 * mass / 6.0);
        // The result does not depend on where the mesh is.
        let moved = transformed(&mesh, &Mat3::identity(), Vec3::new(100.0, -50.0, 25.0));
        assert_matrix_close(&moved.inertia_tensor(mass), &inertia, 1e-10);
        let shifted_com = com + Vec3::new(100.0, -50.0, 25.0);
        assert!((moved.center_of_mass() - shifted_com).norm() < 1e-10);
    }

    #[test]
    fn tilted_flat_plate_is_a_thin_shell_about_its_area_centroid() {
        // The centre of the bounding box of this triangle is off its plane, so a
        // signed volume measured from the box centre would be spurious.
        let a = Vec3::new(0.0, 0.0, 0.0);
        let b = Vec3::new(1.0, 0.0, 0.5);
        let c = Vec3::new(0.0, 1.0, 0.2);
        let mesh = TriangleMesh::new(vec![a, b, c], vec![[0, 1, 2]]);
        let mass = 3.0;
        let centroid = (a + b + c) / 3.0;
        assert!((mesh.center_of_mass() - centroid).norm() < 1e-14);
        // Triangular lamina: m / 12 (tr(Q) 1 - Q) with Q = sum (v - g)(v - g)^T.
        let mut q = Mat3::zeros();
        for v in [a, b, c] {
            let d = v - centroid;
            q += d * d.transpose();
        }
        let expected = (Mat3::identity() * q.trace() - q) * (mass / 12.0);
        let inertia = mesh.inertia_tensor(mass);
        assert_matrix_close(&inertia, &expected, 1e-14);
        // About the plate normal the moment is m (a^2 + b^2 + c^2) / 36 with the
        // squared side lengths, and the trace of a lamina is twice that.
        let normal = (b - a).cross(&(c - a)).normalize();
        let sides = (b - a).norm_squared() + (c - b).norm_squared() + (a - c).norm_squared();
        let about_normal = (normal.transpose() * inertia * normal)[(0, 0)];
        assert!((about_normal - mass * sides / 36.0).abs() < 1e-13);
        assert!((inertia.trace() - 2.0 * about_normal).abs() < 1e-13);
    }

    #[test]
    fn empty_mesh_and_invalid_mass_give_the_zero_tensor() {
        let empty = TriangleMesh::new(vec![], vec![]);
        assert_eq!(empty.inertia_tensor(1.0), Mat3::zeros());
        assert_eq!(empty.center_of_mass(), Vec3::zeros());
        let cube = unit_cube_mesh();
        assert_eq!(cube.inertia_tensor(0.0), Mat3::zeros());
        assert_eq!(cube.inertia_tensor(-1.0), Mat3::zeros());
        assert_eq!(cube.inertia_tensor(f64::NAN), Mat3::zeros());
    }

    #[test]
    fn out_of_range_triangle_is_ignored_without_panicking() {
        let mut mesh = unit_cube_mesh();
        mesh.indices.push([0, 1, 99]);
        let inertia = mesh.inertia_tensor(3.0);
        assert_matrix_close(&inertia, &(Mat3::identity() * 0.5), 1e-14);
    }

    #[test]
    fn explicit_centre_of_mass_matches_the_shape_trait() {
        let mesh = transformed(
            &box_mesh(1.0, 2.0, 3.0),
            &rotation_z(0.7),
            Vec3::new(-4.0, 5.0, 6.0),
        );
        let explicit = mesh.center_of_mass_explicit();
        let via_trait = mesh.center_of_mass();
        assert!((Vec3::new(explicit[0], explicit[1], explicit[2]) - via_trait).norm() < 1e-14);
    }
}
