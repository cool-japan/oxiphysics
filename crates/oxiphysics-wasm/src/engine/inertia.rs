// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Mass properties (volume and rotational inertia) of collider shapes.
//!
//! All formulas are for solid bodies of uniform density.  Capsules and
//! cylinders are symmetric about their local Y axis (see
//! [`shapes`](super::shapes)).

use std::f64::consts::PI;

use super::geom::{
    MAT3_IDENTITY, MAT3_ZERO, Mat3, Vec3, dot, length_sq, mat_add, mat_diag, mat_mul, mat_outer,
    mat_scale, mat_transpose,
};
use super::shapes::Shape;

/// Volume of a solid shape in cubic metres (zero for planes).
pub(super) fn volume(shape: &Shape) -> f64 {
    match *shape {
        Shape::Sphere { radius } => 4.0 / 3.0 * PI * radius.powi(3),
        Shape::Capsule {
            radius,
            half_height,
        } => PI * radius * radius * (2.0 * half_height) + 4.0 / 3.0 * PI * radius.powi(3),
        Shape::Box { half_extents } => 8.0 * half_extents[0] * half_extents[1] * half_extents[2],
        Shape::Cylinder {
            radius,
            half_height,
        } => PI * radius * radius * (2.0 * half_height),
        Shape::Plane { .. } => 0.0,
    }
}

/// Principal moments of inertia `[Ixx, Iyy, Izz]` of a solid shape of the
/// given `mass`, taken about the shape's centre along its local axes.
///
/// Returns `None` for planes, which have no volume and therefore no inertia.
pub(super) fn principal_moments(shape: &Shape, mass: f64) -> Option<Vec3> {
    match *shape {
        Shape::Sphere { radius } => {
            let i = 0.4 * mass * radius * radius;
            Some([i, i, i])
        }
        Shape::Box { half_extents: h } => Some([
            mass / 3.0 * (h[1] * h[1] + h[2] * h[2]),
            mass / 3.0 * (h[0] * h[0] + h[2] * h[2]),
            mass / 3.0 * (h[0] * h[0] + h[1] * h[1]),
        ]),
        Shape::Cylinder {
            radius,
            half_height,
        } => {
            let length = 2.0 * half_height;
            let axial = 0.5 * mass * radius * radius;
            let transverse = mass * (3.0 * radius * radius + length * length) / 12.0;
            Some([transverse, axial, transverse])
        }
        Shape::Capsule {
            radius,
            half_height,
        } => {
            let length = 2.0 * half_height;
            let cylinder_volume = PI * radius * radius * length;
            let sphere_volume = 4.0 / 3.0 * PI * radius.powi(3);
            let total = cylinder_volume + sphere_volume;
            let cylinder_mass = mass * cylinder_volume / total;
            let sphere_mass = mass * sphere_volume / total;
            // About the symmetry axis: solid cylinder plus a solid sphere
            // (the two hemispheres together).
            let axial = 0.5 * cylinder_mass * radius * radius + 0.4 * sphere_mass * radius * radius;
            // Transverse: cylinder about a diameter, plus the two hemispheres
            // shifted by the parallel-axis theorem (each hemisphere's centre
            // of mass sits 3r/8 from its flat face, which is L/2 from the
            // capsule centre).
            let transverse = cylinder_mass * (radius * radius / 4.0 + length * length / 12.0)
                + sphere_mass
                    * (0.4 * radius * radius + length * length / 4.0 + 3.0 * length * radius / 8.0);
            Some([transverse, axial, transverse])
        }
        Shape::Plane { .. } => None,
    }
}

/// A shape together with its pose in the body frame.
#[derive(Debug, Clone, Copy)]
pub(super) struct ShapePart {
    /// Local geometry.
    pub shape: Shape,
    /// Position of the shape origin in the body frame.
    pub position: Vec3,
    /// Rotation of the shape frame relative to the body frame.
    pub rotation: Mat3,
}

/// Body-frame inertia tensor of a body of total `mass` made of `parts`.
///
/// The mass is distributed over the parts in proportion to their volumes
/// (uniform density) and each part contributes its own rotated inertia plus
/// the parallel-axis term for its offset from the body origin.  The body
/// origin is treated as the centre of mass.
///
/// Returns `None` when no part encloses a volume (for instance a body that
/// only carries plane colliders) or the result is not finite.
pub(super) fn compound_inertia(mass: f64, parts: &[ShapePart]) -> Option<Mat3> {
    let total_volume: f64 = parts.iter().map(|p| volume(&p.shape)).sum();
    if !(total_volume.is_finite() && total_volume > 0.0 && mass.is_finite() && mass > 0.0) {
        return None;
    }
    let mut tensor = MAT3_ZERO;
    for part in parts {
        let part_volume = volume(&part.shape);
        if part_volume <= 0.0 {
            continue;
        }
        let part_mass = mass * part_volume / total_volume;
        let moments = principal_moments(&part.shape, part_mass)?;
        let own = mat_mul(
            &mat_mul(&part.rotation, &mat_diag(moments)),
            &mat_transpose(&part.rotation),
        );
        // Parallel-axis theorem: m * (|d|^2 * E - d d^T).
        let d = part.position;
        let shift = mat_scale(
            &mat_add(
                &mat_scale(&MAT3_IDENTITY, length_sq(d)),
                &mat_scale(&mat_outer(d, d), -1.0),
            ),
            part_mass,
        );
        tensor = mat_add(&tensor, &mat_add(&own, &shift));
    }
    if tensor.iter().flatten().all(|v| v.is_finite()) && dot(tensor[0], tensor[0]) > 0.0 {
        Some(tensor)
    } else {
        None
    }
}
