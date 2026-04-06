// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Collision detection algorithms including broad-phase and narrow-phase.
//!
//! Provides broadphase (brute force, sweep-and-prune, BVH), narrowphase
//! (GJK, EPA, specialized sphere/box/capsule tests), and contact manifold
//! generation.
#![allow(missing_docs)]
#![allow(ambiguous_glob_reexports)]
#![allow(dead_code)]

mod error;
pub use error::*;

pub mod types;
pub use types::{
    CollisionEvent, CollisionFilter, CollisionPair, Contact, ContactManifold, FeatureId,
    MAX_CONTACTS, PhysicsMaterial, RichContact, RichContactManifold,
};

pub mod broadphase;
pub use broadphase::{BroadPhase, BruteForceBroadPhase, BvhBroadphase, SweepAndPrune};

pub mod dbvt;
pub use dbvt::{BvhAabb, DynamicBvh};

pub mod narrowphase;
pub use narrowphase::{
    BatchNarrowPhase, ContactFilter, NarrowPhaseContact, ShapeKind, shape_shape_contact,
};
pub use narrowphase::{Epa, Gjk, NarrowPhaseDispatcher};

pub mod ccd;
pub use ccd::{
    CcdBodyEntry, CcdPair, CcdPipeline, CcdPipelineConfig, ConservativeAdvancement, RigidBodyState,
    ToiResult, advance_body_to, integrate_angular, integrate_linear, quat_mul, quat_normalize,
    quat_rotate, slerp_quat, time_of_impact,
};

pub mod sweep;
pub use sweep::{
    RayHit, SweepResult, aabb_aabb_overlap, point_aabb_dist_sq, ray_aabb, ray_sphere, ray_triangle,
    sphere_aabb_overlap, swept_aabb_aabb, swept_sphere_aabb, swept_sphere_sphere,
};

pub mod sap;

pub mod query;

pub mod gjk_epa;
pub use gjk_epa::{
    ConvexShape, EpaFace, EpaResult, GjkBox, GjkCapsule, GjkSphere, Simplex, add3, cross3,
    do_simplex, dot3, epa, gjk_distance, gjk_epa_contact, gjk_intersect, len3, minkowski_support,
    normalize3, scale3, sub3,
};

pub mod contact_graph;
pub use contact_graph::{
    ContactCache, ContactGraph, ContactKey, PersistedContact, SpeculativeContact,
    speculative_contact, speculative_impulse,
};

#[allow(missing_docs)]
pub mod manifold_cache;
pub use manifold_cache::*;

#[allow(missing_docs)]
pub mod sat_collision;
pub use sat_collision::*;

#[cfg(test)]
mod prop_tests {

    use oxiphysics_core::Transform;
    use oxiphysics_core::math::Vec3;
    use oxiphysics_geometry::Sphere;
    use proptest::prelude::*;

    use crate::broadphase::{BroadPhase, BruteForceBroadPhase};
    use crate::narrowphase::{Gjk, sphere_sphere};
    use oxiphysics_core::Aabb;

    fn positive_f64() -> impl Strategy<Value = f64> {
        0.01_f64..20.0_f64
    }

    fn coord_f64() -> impl Strategy<Value = f64> {
        -50.0_f64..50.0_f64
    }

    fn aabb_strategy() -> impl Strategy<Value = Aabb> {
        (
            coord_f64(),
            coord_f64(),
            coord_f64(),
            positive_f64(),
            positive_f64(),
            positive_f64(),
        )
            .prop_map(|(x, y, z, ex, ey, ez)| {
                Aabb::new(Vec3::new(x, y, z), Vec3::new(x + ex, y + ey, z + ez))
            })
    }

    proptest! {
        #[test]
        fn prop_gjk_separated_spheres_no_intersection(
            r1 in positive_f64(),
            r2 in positive_f64(),
            // separation > sum of radii
            extra_sep in 0.01_f64..10.0_f64,
        ) {
            let s1 = Sphere::new(r1);
            let s2 = Sphere::new(r2);
            let t1 = Transform::from_position(Vec3::zeros());
            let sep = r1 + r2 + extra_sep;
            let t2 = Transform::from_position(Vec3::new(sep, 0.0, 0.0));
            let intersects = Gjk::intersect(&s1, &t1, &s2, &t2);
            prop_assert!(!intersects, "separated spheres (r1={}, r2={}, sep={}) reported intersecting", r1, r2, sep);
        }

        #[test]
        fn prop_gjk_overlapping_spheres_intersect(
            r1 in positive_f64(),
            r2 in positive_f64(),
            // overlap: distance < sum of radii
            dist_frac in 0.0_f64..0.99_f64,
        ) {
            let s1 = Sphere::new(r1);
            let s2 = Sphere::new(r2);
            let t1 = Transform::from_position(Vec3::zeros());
            let dist = (r1 + r2) * dist_frac;
            let t2 = Transform::from_position(Vec3::new(dist, 0.0, 0.0));
            let intersects = Gjk::intersect(&s1, &t1, &s2, &t2);
            prop_assert!(intersects, "overlapping spheres (r1={}, r2={}, dist={}) not reported intersecting", r1, r2, dist);
        }

        #[test]
        fn prop_broadphase_pairs_have_overlapping_aabbs(
            aabbs in proptest::collection::vec(aabb_strategy(), 1..8),
        ) {
            let pairs = BruteForceBroadPhase.find_pairs(&aabbs);
            for pair in &pairs {
                let a = &aabbs[pair.a];
                let b = &aabbs[pair.b];
                prop_assert!(
                    a.intersects(b),
                    "pair ({},{}) reported but AABBs don't overlap", pair.a, pair.b
                );
            }
        }

        #[test]
        fn prop_contact_normal_unit_length(
            r1 in positive_f64(),
            r2 in positive_f64(),
            dist_frac in 0.01_f64..0.99_f64,
        ) {
            let s1 = Sphere::new(r1);
            let s2 = Sphere::new(r2);
            let t1 = Transform::from_position(Vec3::zeros());
            let dist = (r1 + r2) * dist_frac;
            let t2 = Transform::from_position(Vec3::new(dist, 0.0, 0.0));
            if let Some(contact) = sphere_sphere(&s1, &t1, &s2, &t2) {
                let nlen = contact.normal.norm();
                prop_assert!((nlen - 1.0).abs() < 1e-6, "contact normal length={}", nlen);
            }
        }

        #[test]
        fn prop_contact_depth_non_negative(
            r1 in positive_f64(),
            r2 in positive_f64(),
            dist_frac in 0.01_f64..0.99_f64,
        ) {
            let s1 = Sphere::new(r1);
            let s2 = Sphere::new(r2);
            let t1 = Transform::from_position(Vec3::zeros());
            let dist = (r1 + r2) * dist_frac;
            let t2 = Transform::from_position(Vec3::new(dist, 0.0, 0.0));
            if let Some(contact) = sphere_sphere(&s1, &t1, &s2, &t2) {
                prop_assert!(contact.depth >= 0.0, "contact depth={} is negative", contact.depth);
            }
        }
    }
}
pub mod compound_shapes;
pub mod contact_generation;
pub mod contact_manifold;
pub mod deformable_collision;
pub mod fluid_collision;
pub mod gjk_enhanced;
pub mod gjk_extended;
pub mod joint_collision;
pub mod kdtree_collision;
pub mod mesh_collision;
pub mod parallel_collision;
pub mod proximity;
pub mod proximity_query;
pub mod ray_casting;
pub mod shape_cast;
pub mod soft_body_collision;
pub mod softbody_collision;
pub mod spatial_queries;
pub mod terrain_collision;
pub mod voxel_collision;
