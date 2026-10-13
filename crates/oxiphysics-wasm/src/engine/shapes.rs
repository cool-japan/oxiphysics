// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Collider geometry: local shapes, their placement in the world and the
//! convex "cores" consumed by the generic GJK/EPA narrow phase.
//!
//! ## Conventions
//!
//! * Sphere - centred on the collider origin.
//! * Box - oriented box with the given half-extents along the collider's
//!   local X, Y and Z axes.
//! * Capsule / Cylinder - symmetric about the collider's local **Y** axis.
//!   `height` is the length of the cylindrical part (for a capsule the two
//!   hemispherical caps come on top of it), matching
//!   [`ColliderConfig::volume`](crate::types::ColliderConfig::volume).
//! * Plane - the half-space boundary `normal . x = offset` expressed in the
//!   collider's local frame; the normal points to the free (outside) side.

use super::geom::{
    self, Mat3, Vec3, add, any_perpendicular, dot, madd, mat_col, mat_mul, mat_vec, scale, sub,
    unit_perpendicular_part,
};
use crate::types::{ColliderConfig, ColliderShapeType};

/// Collider geometry in a single reference frame (local or world).
#[derive(Debug, Clone, Copy, PartialEq)]
pub(super) enum Shape {
    /// Sphere of the given radius.
    Sphere {
        /// Radius (metres).
        radius: f64,
    },
    /// Capsule around the local Y axis.
    Capsule {
        /// Radius of the cylinder and of both hemispherical caps.
        radius: f64,
        /// Half the length of the cylindrical part.
        half_height: f64,
    },
    /// Oriented box.
    Box {
        /// Half-extents along the local X, Y and Z axes.
        half_extents: Vec3,
    },
    /// Solid cylinder around the local Y axis.
    Cylinder {
        /// Radius of the circular cross-section.
        radius: f64,
        /// Half the cylinder length.
        half_height: f64,
    },
    /// Infinite plane `normal . x = offset` (normal is a unit vector).
    Plane {
        /// Unit normal pointing away from the solid half-space.
        normal: Vec3,
        /// Signed distance of the plane from the origin along `normal`.
        offset: f64,
    },
}

impl Shape {
    /// Build the shape described by a (validated) collider configuration.
    ///
    /// Returns `None` for [`ColliderShapeType::ConvexHull`], whose geometry is
    /// not carried by [`ColliderConfig`].
    pub(super) fn from_config(cfg: &ColliderConfig) -> Option<Self> {
        match cfg.shape_type {
            ColliderShapeType::Sphere => Some(Self::Sphere { radius: cfg.radius }),
            ColliderShapeType::Box => Some(Self::Box {
                half_extents: cfg.half_extents,
            }),
            ColliderShapeType::Capsule => Some(Self::Capsule {
                radius: cfg.radius,
                half_height: 0.5 * cfg.height,
            }),
            ColliderShapeType::Cylinder => Some(Self::Cylinder {
                radius: cfg.radius,
                half_height: 0.5 * cfg.height,
            }),
            ColliderShapeType::Plane => Some(Self::Plane {
                normal: cfg.plane_normal,
                offset: cfg.plane_offset,
            }),
            ColliderShapeType::ConvexHull => None,
        }
    }

    /// Ordering used by the narrow phase dispatcher: the lower-ranked shape of
    /// a pair always plays the role of body A in a reported contact.
    pub(super) fn rank(&self) -> u8 {
        match self {
            Self::Sphere { .. } => 0,
            Self::Capsule { .. } => 1,
            Self::Box { .. } => 2,
            Self::Cylinder { .. } => 3,
            Self::Plane { .. } => 4,
        }
    }

    /// Radius of a sphere centred on the shape origin that encloses the shape
    /// (`INFINITY` for planes).
    pub(super) fn bound_radius(&self) -> f64 {
        match *self {
            Self::Sphere { radius } => radius,
            Self::Capsule {
                radius,
                half_height,
            } => radius + half_height,
            Self::Box { half_extents } => geom::length(half_extents),
            Self::Cylinder {
                radius,
                half_height,
            } => radius.hypot(half_height),
            Self::Plane { .. } => f64::INFINITY,
        }
    }
}

/// A convex core that can be queried through its support function.
///
/// Every convex collider is modelled as a core inflated by a rounding margin
/// (see [`WorldShape::core`]): a sphere is a point with margin `radius`, a
/// capsule a segment with margin `radius`, boxes and cylinders have a zero
/// margin.  GJK and EPA operate on the cores only, which keeps them exact and
/// well conditioned for round shapes.
#[derive(Debug, Clone, Copy)]
pub(super) enum Core {
    /// A single point.
    Point(Vec3),
    /// A line segment between two points.
    Segment(Vec3, Vec3),
    /// An oriented box.
    Box {
        /// Box centre.
        center: Vec3,
        /// Unit direction of each box axis in world space.
        axes: [Vec3; 3],
        /// Half-extents along `axes`.
        half: Vec3,
    },
    /// A solid cylinder.
    Cylinder {
        /// Cylinder centre.
        center: Vec3,
        /// Unit symmetry axis in world space.
        axis: Vec3,
        /// Radius of the cross-section.
        radius: f64,
        /// Half the cylinder length.
        half_height: f64,
    },
}

impl Core {
    /// Point of the core that is farthest along `dir` (which need not be a
    /// unit vector).
    pub(super) fn support(&self, dir: Vec3) -> Vec3 {
        match *self {
            Self::Point(p) => p,
            Self::Segment(a, b) => {
                if dot(dir, a) >= dot(dir, b) {
                    a
                } else {
                    b
                }
            }
            Self::Box { center, axes, half } => {
                let mut p = center;
                for (axis, h) in axes.iter().zip(half.iter()) {
                    let s = if dot(dir, *axis) >= 0.0 { *h } else { -*h };
                    p = madd(p, *axis, s);
                }
                p
            }
            Self::Cylinder {
                center,
                axis,
                radius,
                half_height,
            } => {
                let along = dot(dir, axis);
                let end = if along >= 0.0 {
                    half_height
                } else {
                    -half_height
                };
                let mut p = madd(center, axis, end);
                // A direction along the axis has no radial part: the support
                // point is then the centre of the cap.
                if let Some(u) = unit_perpendicular_part(dir, axis) {
                    p = madd(p, u, radius);
                }
                p
            }
        }
    }

    /// A point strictly inside (or on) the core, used to seed GJK.
    pub(super) fn center(&self) -> Vec3 {
        match *self {
            Self::Point(p) => p,
            Self::Segment(a, b) => scale(add(a, b), 0.5),
            Self::Box { center, .. } | Self::Cylinder { center, .. } => center,
        }
    }
}

/// A collider placed in the world: its shape, pose and a bounding sphere.
///
/// For [`Shape::Plane`] the stored normal and offset are already expressed in
/// world coordinates.
#[derive(Debug, Clone, Copy)]
pub(super) struct WorldShape {
    /// World-space geometry (planes carry world-space normal and offset).
    pub shape: Shape,
    /// World-space position of the collider origin.
    pub position: Vec3,
    /// Collider-to-world rotation (columns are the local axes in world space).
    pub rotation: Mat3,
    /// Radius of a sphere around `position` enclosing the shape.
    pub bound_radius: f64,
}

impl WorldShape {
    /// Place a local `shape` into the world.
    ///
    /// The collider pose is the body pose composed with the collider's local
    /// offset: `world = body_position + body_rotation * local_position` and
    /// `rotation = body_rotation * local_rotation`.
    pub(super) fn place(
        shape: Shape,
        body_position: Vec3,
        body_rotation: &Mat3,
        local_position: Vec3,
        local_rotation: &Mat3,
    ) -> Self {
        let position = add(body_position, mat_vec(body_rotation, local_position));
        let rotation = mat_mul(body_rotation, local_rotation);
        let placed = match shape {
            Shape::Plane { normal, offset } => {
                let world_normal = mat_vec(&rotation, normal);
                Shape::Plane {
                    normal: world_normal,
                    offset: offset + dot(world_normal, position),
                }
            }
            other => other,
        };
        Self {
            shape: placed,
            position,
            rotation,
            bound_radius: placed.bound_radius(),
        }
    }

    /// World-space direction of the `i`-th local axis (`0 = X`, `1 = Y`, `2 = Z`).
    #[inline]
    pub(super) fn axis(&self, i: usize) -> Vec3 {
        mat_col(&self.rotation, i)
    }

    /// Transform a world-space point into the collider's local frame.
    #[inline]
    pub(super) fn to_local_point(self, p: Vec3) -> Vec3 {
        geom::mat_t_vec(&self.rotation, sub(p, self.position))
    }

    /// Transform a world-space direction into the collider's local frame.
    #[inline]
    pub(super) fn to_local_dir(self, d: Vec3) -> Vec3 {
        geom::mat_t_vec(&self.rotation, d)
    }

    /// Transform a local-frame direction into world space.
    #[inline]
    pub(super) fn to_world_dir(self, d: Vec3) -> Vec3 {
        mat_vec(&self.rotation, d)
    }

    /// Transform a local-frame point into world space.
    #[inline]
    pub(super) fn to_world_point(self, p: Vec3) -> Vec3 {
        add(self.position, mat_vec(&self.rotation, p))
    }

    /// Endpoints of the axis segment of a capsule (or the axis of a cylinder).
    ///
    /// Returns `None` for shapes without a Y-symmetry axis.
    pub(super) fn axis_endpoints(&self) -> Option<(Vec3, Vec3)> {
        let half_height = match self.shape {
            Shape::Capsule { half_height, .. } | Shape::Cylinder { half_height, .. } => half_height,
            _ => return None,
        };
        let axis = self.axis(1);
        Some((
            madd(self.position, axis, -half_height),
            madd(self.position, axis, half_height),
        ))
    }

    /// The convex core of the shape together with its rounding margin, or
    /// `None` for planes (which are handled analytically).
    pub(super) fn core(&self) -> Option<(Core, f64)> {
        match self.shape {
            Shape::Sphere { radius } => Some((Core::Point(self.position), radius)),
            Shape::Capsule {
                radius,
                half_height,
            } => {
                let axis = self.axis(1);
                Some((
                    Core::Segment(
                        madd(self.position, axis, -half_height),
                        madd(self.position, axis, half_height),
                    ),
                    radius,
                ))
            }
            Shape::Box { half_extents } => Some((
                Core::Box {
                    center: self.position,
                    axes: [self.axis(0), self.axis(1), self.axis(2)],
                    half: half_extents,
                },
                0.0,
            )),
            Shape::Cylinder {
                radius,
                half_height,
            } => Some((
                Core::Cylinder {
                    center: self.position,
                    axis: self.axis(1),
                    radius,
                    half_height,
                },
                0.0,
            )),
            Shape::Plane { .. } => None,
        }
    }

    /// A unit vector perpendicular to the symmetry axis (used for degenerate
    /// contact normals of round shapes).
    pub(super) fn any_perpendicular_to_axis(&self) -> Vec3 {
        any_perpendicular(self.axis(1))
    }
}
