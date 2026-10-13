// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Validation of user supplied body and collider configurations.
//!
//! Everything that crosses the WASM boundary is untrusted: a NaN position, a
//! negative radius or a zero plane normal would otherwise be stored silently
//! and surface much later as NaN poisoning or a body that never collides.
//! The functions here reject such input with a descriptive [`Error`] and hand
//! back a normalised copy (unit quaternions, unit plane normals) of the
//! accepted configuration.

use super::geom::{Vec3, length, normalized, quat_normalized};
use crate::error::{Error, Result};
use crate::types::{ColliderConfig, ColliderShapeType, RigidBodyConfig};

/// `true` when all three components are finite.
fn finite3(v: Vec3) -> bool {
    v.iter().all(|c| c.is_finite())
}

/// Require `value` to be finite and strictly positive.
fn positive_dimension(value: f64) -> Result<()> {
    if value.is_finite() && value > 0.0 {
        Ok(())
    } else {
        Err(Error::InvalidDimension(value))
    }
}

/// Require all three components of a vector argument to be finite.
///
/// Used by the body setters: a single `NaN` or infinity written into a
/// velocity or position would poison the body (and, through contacts, its
/// neighbours) for the rest of the simulation.
pub(super) fn require_finite3(name: &str, value: Vec3) -> Result<()> {
    if finite3(value) {
        Ok(())
    } else {
        Err(Error::invalid_param(name, "components must be finite"))
    }
}

/// Require a finite, non-negative scalar parameter.
pub(super) fn non_negative_param(name: &str, value: f64) -> Result<()> {
    if value.is_finite() && value >= 0.0 {
        Ok(())
    } else {
        Err(Error::invalid_param(
            name,
            format!("must be finite and >= 0, got {value}"),
        ))
    }
}

/// Validate a mass given for a *dynamic* body: finite and strictly positive.
pub(super) fn validate_dynamic_mass(mass: f64) -> Result<()> {
    if mass.is_finite() && mass > 0.0 {
        Ok(())
    } else {
        Err(Error::InvalidMass(mass))
    }
}

/// Validate a rigid-body configuration and return a normalised copy.
///
/// A body with `mass == 0` is accepted (static, or kinematic when flagged);
/// any negative or non-finite mass is rejected.
pub(super) fn validate_body_config(cfg: &RigidBodyConfig) -> Result<RigidBodyConfig> {
    if !(cfg.mass.is_finite() && cfg.mass >= 0.0) {
        return Err(Error::InvalidMass(cfg.mass));
    }
    if !finite3(cfg.position) {
        return Err(Error::invalid_param(
            "position",
            "components must be finite",
        ));
    }
    let rotation = quat_normalized(cfg.rotation)
        .ok_or_else(|| Error::invalid_param("rotation", "must be a finite, non-zero quaternion"))?;
    if !finite3(cfg.linear_velocity) {
        return Err(Error::invalid_param(
            "linear_velocity",
            "components must be finite",
        ));
    }
    if !finite3(cfg.angular_velocity) {
        return Err(Error::invalid_param(
            "angular_velocity",
            "components must be finite",
        ));
    }
    non_negative_param("linear_damping", cfg.linear_damping)?;
    non_negative_param("angular_damping", cfg.angular_damping)?;
    if !(cfg.restitution.is_finite() && (0.0..=1.0).contains(&cfg.restitution)) {
        return Err(Error::InvalidRestitution(cfg.restitution));
    }
    if !(cfg.friction.is_finite() && cfg.friction >= 0.0) {
        return Err(Error::InvalidFriction(cfg.friction));
    }
    if !cfg.gravity_scale.is_finite() {
        return Err(Error::invalid_param("gravity_scale", "must be finite"));
    }
    let mut checked = cfg.clone();
    checked.rotation = rotation;
    Ok(checked)
}

/// Validate a collider configuration and return a normalised copy.
///
/// * radii, half-extents and heights must be finite and strictly positive
///   ([`Error::InvalidDimension`]);
/// * a plane needs a finite, non-zero normal ([`Error::InvalidParameter`]
///   named `plane_normal`); the returned plane has a unit normal and the
///   offset is rescaled so that it describes the same plane;
/// * the local offset and rotation must be finite (the rotation is
///   normalised);
/// * convex hulls cannot be created because [`ColliderConfig`] does not carry
///   their vertices ([`Error::UnsupportedShape`]).
pub(super) fn validate_collider_config(cfg: &ColliderConfig) -> Result<ColliderConfig> {
    let mut checked = cfg.clone();
    match cfg.shape_type {
        ColliderShapeType::Sphere => positive_dimension(cfg.radius)?,
        ColliderShapeType::Box => {
            for h in cfg.half_extents {
                positive_dimension(h)?;
            }
        }
        ColliderShapeType::Capsule | ColliderShapeType::Cylinder => {
            positive_dimension(cfg.radius)?;
            positive_dimension(cfg.height)?;
        }
        ColliderShapeType::Plane => {
            let len = length(cfg.plane_normal);
            let unit = if finite3(cfg.plane_normal) {
                normalized(cfg.plane_normal)
            } else {
                None
            };
            let unit = unit.ok_or_else(|| {
                Error::invalid_param("plane_normal", "must be a finite, non-zero vector")
            })?;
            if !cfg.plane_offset.is_finite() {
                return Err(Error::invalid_param("plane_offset", "must be finite"));
            }
            checked.plane_normal = unit;
            checked.plane_offset = cfg.plane_offset / len;
        }
        ColliderShapeType::ConvexHull => {
            return Err(Error::UnsupportedShape(
                "ConvexHull colliders cannot be created: ColliderConfig carries no vertex data"
                    .to_string(),
            ));
        }
    }
    if !finite3(cfg.local_position) {
        return Err(Error::invalid_param(
            "local_position",
            "components must be finite",
        ));
    }
    checked.local_rotation = quat_normalized(cfg.local_rotation).ok_or_else(|| {
        Error::invalid_param("local_rotation", "must be a finite, non-zero quaternion")
    })?;
    if !(cfg.density.is_finite() && cfg.density > 0.0) {
        return Err(Error::invalid_param(
            "density",
            format!("must be finite and > 0, got {}", cfg.density),
        ));
    }
    Ok(checked)
}
