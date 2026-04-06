// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Prelude module for OxiPhysics.
//!
//! Import everything you need to get started with a single glob import:
//!
//! ```no_run
//! use oxiphysics::prelude::*;
//! ```
//!
//! The prelude re-exports the most commonly used types, traits, and functions
//! from all sub-crates so that users do not have to know the full module
//! hierarchy when prototyping or writing straightforward simulations.
//!
//! # Prelude Pattern
//!
//! The prelude pattern is widely used in the Rust ecosystem (e.g. `std::prelude`,
//! `nalgebra::prelude`).  It is intentionally *not* exhaustive — advanced or
//! niche items remain in their own sub-modules — but it covers the 80 % case.

// ---------------------------------------------------------------------------
// Core math re-exports
// ---------------------------------------------------------------------------

/// The default floating-point scalar type used throughout OxiPhysics.
pub use oxiphysics_core::math::Real;

/// 3-D vector type.
pub use oxiphysics_core::math::Vec3;

/// 4-D vector type.
pub use oxiphysics_core::math::Vec4;

/// 3×3 matrix type.
pub use oxiphysics_core::math::Mat3;

/// 4×4 matrix type.
pub use oxiphysics_core::math::Mat4;

/// Unit quaternion representing a 3-D rotation.
pub use oxiphysics_core::math::Quat;

/// Axis-aligned bounding box.
pub use oxiphysics_core::Aabb;

/// Rigid-body transform (position + orientation).
pub use oxiphysics_core::Transform;

/// Mass properties (mass, centre of mass, inertia tensor).
pub use oxiphysics_core::MassProperties;

/// Opaque handle to a rigid body inside a world.
pub use oxiphysics_core::BodyHandle;

/// Opaque handle to a collider inside a world.
pub use oxiphysics_core::ColliderHandle;

/// Top-level physics engine configuration.
pub use oxiphysics_core::PhysicsConfig;

/// Fixed-step time-step descriptor.
pub use oxiphysics_core::TimeStep;

// ---------------------------------------------------------------------------
// Core math utility re-exports
// ---------------------------------------------------------------------------

/// Construct a quaternion from an axis and angle (radians).
pub use oxiphysics_core::math::quat_from_axis_angle;

/// Convert a quaternion to Euler angles (roll, pitch, yaw).
pub use oxiphysics_core::math::quat_to_euler;

/// Construct a quaternion from Euler angles (roll, pitch, yaw).
pub use oxiphysics_core::math::quat_from_euler;

/// Spherical linear interpolation between two quaternions.
pub use oxiphysics_core::math::quat_slerp;

/// Build a perspective projection matrix.
pub use oxiphysics_core::math::perspective;

/// Build a look-at view matrix.
pub use oxiphysics_core::math::look_at;

// ---------------------------------------------------------------------------
// Rigid-body re-exports
// ---------------------------------------------------------------------------

/// Rigid body simulation world.
pub use oxiphysics_rigid::world::RigidWorld;

/// Configuration for a rigid body world.
pub use oxiphysics_rigid::world::WorldConfig;

/// World-level statistics (contact counts, sleep counts, …).
pub use oxiphysics_rigid::world::WorldStatistics;

/// The core rigid-body data structure.
pub use oxiphysics_rigid::RigidBody;

/// Distinguishes static, kinematic, and dynamic bodies.
pub use oxiphysics_rigid::BodyType;

/// Observable activation state of a rigid body.
pub use oxiphysics_rigid::BodyState;

/// Set of rigid bodies indexed by [`BodyHandle`].
pub use oxiphysics_rigid::RigidBodySet;

/// Set of colliders indexed by [`ColliderHandle`].
pub use oxiphysics_rigid::ColliderSet;

// ---------------------------------------------------------------------------
// Collision re-exports
// ---------------------------------------------------------------------------

/// A pair of collider indices produced by the broadphase.
pub use oxiphysics_collision::CollisionPair;

/// Contact manifold holding one or more contact points.
pub use oxiphysics_collision::ContactManifold;

/// Single contact point between two shapes.
pub use oxiphysics_collision::Contact;

/// Dispatch table for narrowphase contact generation.
pub use oxiphysics_collision::NarrowPhaseDispatcher;

/// Sweep-and-prune broadphase.
pub use oxiphysics_collision::SweepAndPrune;

/// Brute-force O(n²) broadphase (useful for debugging).
pub use oxiphysics_collision::BruteForceBroadPhase;

/// BVH-based broadphase.
pub use oxiphysics_collision::BvhBroadphase;

// ---------------------------------------------------------------------------
// Geometry re-exports
// ---------------------------------------------------------------------------

/// Abstract shape trait with ray-cast support.
pub use oxiphysics_geometry::Shape;

/// Sphere shape.
pub use oxiphysics_geometry::Sphere;

/// Axis-aligned box shape.
pub use oxiphysics_geometry::BoxShape;

/// Capsule shape (cylinder + hemispherical caps).
pub use oxiphysics_geometry::Capsule;

/// Cylinder shape.
pub use oxiphysics_geometry::Cylinder;

/// Cone shape.
pub use oxiphysics_geometry::Cone;

/// Triangle mesh / polyhedron shape.
pub use oxiphysics_geometry::TriangleMesh;

/// Height-field terrain shape.
pub use oxiphysics_geometry::HeightField;

/// Compound shape composed of multiple primitives.
pub use oxiphysics_geometry::Compound;

/// Torus shape.
pub use oxiphysics_geometry::Torus;

// ---------------------------------------------------------------------------
// Material re-exports
// ---------------------------------------------------------------------------

/// Isotropic linear elastic material.
pub use oxiphysics_materials::elastic::ElasticMaterial;

/// Isotropic linear-elastic constitutive model.
pub use oxiphysics_materials::elastic::IsotropicElastic;

/// Orthotropic elastic constitutive model.
pub use oxiphysics_materials::elastic::OrthotropicElastic;

/// Neo-Hookean hyperelastic constitutive model.
pub use oxiphysics_materials::elastic::NeoHookean;

// ---------------------------------------------------------------------------
// Constraint / joint re-exports
// ---------------------------------------------------------------------------

/// Abstract constraint trait.
pub use oxiphysics_constraints::Constraint;

/// Contact constraint between two rigid bodies.
pub use oxiphysics_constraints::ContactConstraint;

/// Projected Gauss-Seidel constraint solver.
pub use oxiphysics_constraints::PgsSolver;

/// Temporal Gauss-Seidel constraint solver.
pub use oxiphysics_constraints::TgsSolver;

/// Fixed (weld) joint — zero relative degrees of freedom.
pub use oxiphysics_constraints::FixedJoint;

/// Revolute (hinge) joint — one rotational degree of freedom.
pub use oxiphysics_constraints::RevoluteJoint;

/// Prismatic (slider) joint — one translational degree of freedom.
pub use oxiphysics_constraints::PrismaticJoint;

/// Ball-and-socket joint — three rotational degrees of freedom.
pub use oxiphysics_constraints::BallJoint;

/// Spring joint with configurable stiffness and damping.
pub use oxiphysics_constraints::SpringJoint;

/// Motorised joint that can track a target velocity or position.
pub use oxiphysics_constraints::MotorJoint;

/// Spherical joint identical to BallJoint but with limit cone support.
pub use oxiphysics_constraints::SphericalJoint;

/// Distance constraint keeping two bodies a fixed distance apart.
pub use oxiphysics_constraints::DistanceJoint;

/// Weld joint (rigid attachment).
pub use oxiphysics_constraints::WeldJoint;

// ---------------------------------------------------------------------------
// Convenience type aliases
// ---------------------------------------------------------------------------

/// Convenience alias for a physics 3-D vector stored as a plain array.
///
/// Useful when interfacing with external APIs that work with `[f64; 3]`.
pub type PhysicsVec3 = [f64; 3];

/// Convenience alias for the physics scalar type.
pub type PhysicsScalar = f64;

/// Convenience alias for a 3×3 matrix stored as a row-major flat array.
pub type PhysicsMat3 = [f64; 9];

/// Convenience alias for an index into a body collection.
pub type BodyIndex = usize;

/// Convenience alias for an index into a collider collection.
pub type ColliderIndex = usize;

// ---------------------------------------------------------------------------
// Utility helpers
// ---------------------------------------------------------------------------

/// Convert a `[f64; 3]` array to a [`Vec3`].
///
/// # Examples
/// ```no_run
/// use oxiphysics::prelude::{arr_to_vec3, PhysicsVec3};
/// let v = arr_to_vec3([1.0, 2.0, 3.0]);
/// assert!((v.x - 1.0).abs() < 1e-12);
/// ```
#[inline]
pub fn arr_to_vec3(a: PhysicsVec3) -> Vec3 {
    Vec3::new(a[0], a[1], a[2])
}

/// Convert a [`Vec3`] to a `[f64; 3]` array.
///
/// # Examples
/// ```no_run
/// use oxiphysics::prelude::{vec3_to_arr, Vec3};
/// let v = Vec3::new(4.0, 5.0, 6.0);
/// let a = vec3_to_arr(v);
/// assert_eq!(a[0], 4.0);
/// ```
#[inline]
pub fn vec3_to_arr(v: Vec3) -> PhysicsVec3 {
    [v.x, v.y, v.z]
}

/// Compute the Euclidean distance between two `[f64; 3]` points.
///
/// # Examples
/// ```no_run
/// use oxiphysics::prelude::arr_distance;
/// let d = arr_distance([0.0, 0.0, 0.0], [1.0, 0.0, 0.0]);
/// assert!((d - 1.0).abs() < 1e-12);
/// ```
#[inline]
pub fn arr_distance(a: PhysicsVec3, b: PhysicsVec3) -> PhysicsScalar {
    let dx = a[0] - b[0];
    let dy = a[1] - b[1];
    let dz = a[2] - b[2];
    (dx * dx + dy * dy + dz * dz).sqrt()
}

/// Linearly interpolate between two `[f64; 3]` arrays.
///
/// `t = 0.0` returns `a`, `t = 1.0` returns `b`.
///
/// # Examples
/// ```no_run
/// use oxiphysics::prelude::arr_lerp;
/// let m = arr_lerp([0.0, 0.0, 0.0], [2.0, 4.0, 6.0], 0.5);
/// assert!((m[0] - 1.0).abs() < 1e-12);
/// ```
#[inline]
pub fn arr_lerp(a: PhysicsVec3, b: PhysicsVec3, t: PhysicsScalar) -> PhysicsVec3 {
    [
        a[0] + t * (b[0] - a[0]),
        a[1] + t * (b[1] - a[1]),
        a[2] + t * (b[2] - a[2]),
    ]
}

/// Normalise a `[f64; 3]` vector. Returns `[0,0,0]` for zero-length input.
///
/// # Examples
/// ```no_run
/// use oxiphysics::prelude::arr_normalize;
/// let n = arr_normalize([3.0, 0.0, 0.0]);
/// assert!((n[0] - 1.0).abs() < 1e-12);
/// ```
#[inline]
pub fn arr_normalize(v: PhysicsVec3) -> PhysicsVec3 {
    let len = arr_distance([0.0, 0.0, 0.0], v);
    if len < 1e-15 {
        [0.0, 0.0, 0.0]
    } else {
        [v[0] / len, v[1] / len, v[2] / len]
    }
}

/// Dot product of two `[f64; 3]` arrays.
///
/// # Examples
/// ```no_run
/// use oxiphysics::prelude::arr_dot;
/// assert!((arr_dot([1.0, 0.0, 0.0], [1.0, 0.0, 0.0]) - 1.0).abs() < 1e-12);
/// ```
#[inline]
pub fn arr_dot(a: PhysicsVec3, b: PhysicsVec3) -> PhysicsScalar {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// Cross product of two `[f64; 3]` arrays.
///
/// # Examples
/// ```no_run
/// use oxiphysics::prelude::arr_cross;
/// let c = arr_cross([1.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
/// assert!((c[2] - 1.0).abs() < 1e-12);
/// ```
#[inline]
pub fn arr_cross(a: PhysicsVec3, b: PhysicsVec3) -> PhysicsVec3 {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // --- arr_to_vec3 / vec3_to_arr ---

    #[test]
    fn test_arr_to_vec3_basic() {
        let v = arr_to_vec3([1.0, 2.0, 3.0]);
        assert!((v.x - 1.0).abs() < 1e-12);
        assert!((v.y - 2.0).abs() < 1e-12);
        assert!((v.z - 3.0).abs() < 1e-12);
    }

    #[test]
    fn test_vec3_to_arr_basic() {
        let v = Vec3::new(7.0, 8.0, 9.0);
        let a = vec3_to_arr(v);
        assert_eq!(a, [7.0, 8.0, 9.0]);
    }

    #[test]
    fn test_round_trip_arr_vec3() {
        let orig = [1.5, -2.5, 3.14];
        let v = arr_to_vec3(orig);
        let back = vec3_to_arr(v);
        for (a, b) in orig.iter().zip(back.iter()) {
            assert!((a - b).abs() < 1e-12);
        }
    }

    #[test]
    fn test_arr_to_vec3_zero() {
        let v = arr_to_vec3([0.0, 0.0, 0.0]);
        assert_eq!(v.norm(), 0.0);
    }

    #[test]
    fn test_arr_to_vec3_negative() {
        let v = arr_to_vec3([-1.0, -2.0, -3.0]);
        assert!((v.x + 1.0).abs() < 1e-12);
    }

    // --- arr_distance ---

    #[test]
    fn test_arr_distance_unit() {
        let d = arr_distance([0.0, 0.0, 0.0], [1.0, 0.0, 0.0]);
        assert!((d - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_arr_distance_zero() {
        let d = arr_distance([1.0, 2.0, 3.0], [1.0, 2.0, 3.0]);
        assert!(d.abs() < 1e-12);
    }

    #[test]
    fn test_arr_distance_3d() {
        let d = arr_distance([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        assert!((d - 3f64.sqrt()).abs() < 1e-12);
    }

    #[test]
    fn test_arr_distance_negative_coords() {
        let d = arr_distance([-1.0, 0.0, 0.0], [1.0, 0.0, 0.0]);
        assert!((d - 2.0).abs() < 1e-12);
    }

    // --- arr_lerp ---

    #[test]
    fn test_arr_lerp_start() {
        let r = arr_lerp([0.0, 0.0, 0.0], [10.0, 20.0, 30.0], 0.0);
        assert_eq!(r, [0.0, 0.0, 0.0]);
    }

    #[test]
    fn test_arr_lerp_end() {
        let r = arr_lerp([0.0, 0.0, 0.0], [10.0, 20.0, 30.0], 1.0);
        assert!((r[0] - 10.0).abs() < 1e-12);
    }

    #[test]
    fn test_arr_lerp_mid() {
        let r = arr_lerp([0.0, 0.0, 0.0], [2.0, 4.0, 6.0], 0.5);
        assert!((r[0] - 1.0).abs() < 1e-12);
        assert!((r[1] - 2.0).abs() < 1e-12);
        assert!((r[2] - 3.0).abs() < 1e-12);
    }

    #[test]
    fn test_arr_lerp_quarter() {
        let r = arr_lerp([0.0, 0.0, 0.0], [4.0, 0.0, 0.0], 0.25);
        assert!((r[0] - 1.0).abs() < 1e-12);
    }

    // --- arr_normalize ---

    #[test]
    fn test_arr_normalize_x() {
        let n = arr_normalize([5.0, 0.0, 0.0]);
        assert!((n[0] - 1.0).abs() < 1e-12);
        assert!(n[1].abs() < 1e-12);
        assert!(n[2].abs() < 1e-12);
    }

    #[test]
    fn test_arr_normalize_diagonal() {
        let n = arr_normalize([1.0, 1.0, 1.0]);
        let expected = 1.0 / 3f64.sqrt();
        for c in n.iter() {
            assert!((c - expected).abs() < 1e-12);
        }
    }

    #[test]
    fn test_arr_normalize_zero_safe() {
        let n = arr_normalize([0.0, 0.0, 0.0]);
        assert_eq!(n, [0.0, 0.0, 0.0]);
    }

    #[test]
    fn test_arr_normalize_unit_length() {
        let n = arr_normalize([3.0, 4.0, 0.0]);
        let len = arr_distance([0.0, 0.0, 0.0], n);
        assert!((len - 1.0).abs() < 1e-12);
    }

    // --- arr_dot ---

    #[test]
    fn test_arr_dot_parallel() {
        assert!((arr_dot([1.0, 0.0, 0.0], [1.0, 0.0, 0.0]) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_arr_dot_perpendicular() {
        assert!(arr_dot([1.0, 0.0, 0.0], [0.0, 1.0, 0.0]).abs() < 1e-12);
    }

    #[test]
    fn test_arr_dot_anti_parallel() {
        assert!((arr_dot([1.0, 0.0, 0.0], [-1.0, 0.0, 0.0]) + 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_arr_dot_general() {
        let d = arr_dot([1.0, 2.0, 3.0], [4.0, 5.0, 6.0]);
        assert!((d - 32.0).abs() < 1e-12);
    }

    // --- arr_cross ---

    #[test]
    fn test_arr_cross_xy() {
        let c = arr_cross([1.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
        assert!((c[2] - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_arr_cross_anti_commutative() {
        let a = arr_cross([1.0, 2.0, 3.0], [4.0, 5.0, 6.0]);
        let b = arr_cross([4.0, 5.0, 6.0], [1.0, 2.0, 3.0]);
        for i in 0..3 {
            assert!((a[i] + b[i]).abs() < 1e-12);
        }
    }

    #[test]
    fn test_arr_cross_parallel_zero() {
        let c = arr_cross([1.0, 0.0, 0.0], [2.0, 0.0, 0.0]);
        for comp in c.iter() {
            assert!(comp.abs() < 1e-12);
        }
    }

    #[test]
    fn test_arr_cross_magnitude() {
        let c = arr_cross([1.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
        let mag = arr_distance([0.0, 0.0, 0.0], c);
        assert!((mag - 1.0).abs() < 1e-12);
    }

    // --- Type alias sanity checks ---

    #[test]
    fn test_physics_scalar_type() {
        let s: PhysicsScalar = 3.14;
        assert!((s - 3.14).abs() < 1e-12);
    }

    #[test]
    fn test_physics_vec3_type() {
        let v: PhysicsVec3 = [1.0, 2.0, 3.0];
        assert_eq!(v.len(), 3);
    }

    #[test]
    fn test_physics_mat3_type() {
        let m: PhysicsMat3 = [0.0; 9];
        assert_eq!(m.len(), 9);
    }

    #[test]
    fn test_body_index_type() {
        let idx: BodyIndex = 42;
        assert_eq!(idx, 42);
    }

    #[test]
    fn test_collider_index_type() {
        let idx: ColliderIndex = 7;
        assert_eq!(idx, 7);
    }

    // --- PhysicsConfig default construction ---

    #[test]
    fn test_physics_config_default() {
        let _cfg = PhysicsConfig::default();
    }

    // --- Vec3 round-trip ---

    #[test]
    fn test_vec3_negation() {
        let v = arr_to_vec3([1.0, 2.0, 3.0]);
        let neg = -v;
        assert!((neg.x + 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_vec3_addition() {
        let a = arr_to_vec3([1.0, 0.0, 0.0]);
        let b = arr_to_vec3([0.0, 1.0, 0.0]);
        let sum = a + b;
        assert!((sum.x - 1.0).abs() < 1e-12);
        assert!((sum.y - 1.0).abs() < 1e-12);
    }
}
