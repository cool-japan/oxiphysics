// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Rigid body types for articulated-body models.

use crate::spatial::{SpatialInertia, SpatialTransform};

/// A rigid body in an articulated chain.
///
/// Each body carries its inertia properties and the fixed (joint-independent)
/// transform from the parent body's coordinate frame to this body's joint frame.
pub struct RigidBody {
    /// Human-readable name (for debugging).
    pub name: String,
    /// Spatial inertia of this body.
    pub inertia: SpatialInertia,
    /// Index of the parent body, or `None` for the root body.
    pub parent_id: Option<usize>,
    /// Fixed transform from parent body frame to this body's joint frame (`X_T`).
    ///
    /// This is the part of the transform that does not depend on joint coordinates.
    /// The full parent-to-child transform is `X_T ∘ X_J(q)`.
    pub parent_transform: SpatialTransform,
}

impl RigidBody {
    /// Create a new rigid body.
    pub fn new(
        name: impl Into<String>,
        inertia: SpatialInertia,
        parent_id: Option<usize>,
        parent_transform: SpatialTransform,
    ) -> Self {
        Self {
            name: name.into(),
            inertia,
            parent_id,
            parent_transform,
        }
    }
}
