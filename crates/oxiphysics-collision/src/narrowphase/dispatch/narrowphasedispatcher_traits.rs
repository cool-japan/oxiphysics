//! # NarrowPhaseDispatcher - Trait Implementations
//!
//! This module contains trait implementations for `NarrowPhaseDispatcher`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use crate::types::ContactManifold;
use oxiphysics_geometry::{BoxShape, Capsule, Sphere};

use super::box_manifold::box_box_manifold;
use super::functions::{gjk_fallback_dispatch, sphere_capsule_dispatch, type_name};
use super::types::{NarrowPhaseDispatcher, NarrowPhaseResult, ShapeType};
use crate::narrowphase::specialized;

impl Default for NarrowPhaseDispatcher {
    /// Build a dispatcher pre-registered with sphere/box/capsule pair algorithms.
    fn default() -> Self {
        let mut d = NarrowPhaseDispatcher::empty();
        d.register_pair(
            ShapeType::Sphere,
            ShapeType::Sphere,
            |sa, ta, sb, tb, pair| {
                // Recover both `Sphere`s through safe `Any` downcasts. Registered
                // under the (Sphere, Sphere) key, so in canonical order both are
                // spheres; a mismatched direct call downcasts to `None` and falls
                // back to GJK/EPA rather than triggering UB.
                let (Some(s1), Some(s2)) = (
                    sa.as_any().downcast_ref::<Sphere>(),
                    sb.as_any().downcast_ref::<Sphere>(),
                ) else {
                    return gjk_fallback_dispatch(sa, ta, sb, tb, pair);
                };
                debug_assert!(
                    type_name(sa) == "Sphere",
                    "narrowphase dispatch: shape_a must be Sphere for this handler (registration/order invariant)"
                );
                debug_assert!(
                    type_name(sb) == "Sphere",
                    "narrowphase dispatch: shape_b must be Sphere for this handler (registration/order invariant)"
                );
                match specialized::sphere_sphere(s1, ta, s2, tb) {
                    Some(contact) => {
                        let mut m = ContactManifold::new(pair);
                        m.add_contact(contact);
                        NarrowPhaseResult::contact(m)
                    }
                    None => NarrowPhaseResult::separated(),
                }
            },
        );
        d.register_pair(ShapeType::Sphere, ShapeType::Box, |sa, ta, sb, tb, pair| {
            // Recover the concrete shapes through safe `Any` downcasts. Registered
            // under the (Sphere, Box) key, so in canonical order `sa` is a Sphere
            // and `sb` a BoxShape; a mismatched direct call downcasts to `None` and
            // falls back to GJK/EPA rather than triggering UB.
            let (Some(s), Some(b)) = (
                sa.as_any().downcast_ref::<Sphere>(),
                sb.as_any().downcast_ref::<BoxShape>(),
            ) else {
                return gjk_fallback_dispatch(sa, ta, sb, tb, pair);
            };
            debug_assert!(
                type_name(sa) == "Sphere",
                "narrowphase dispatch: shape_a must be Sphere for this handler (registration/order invariant)"
            );
            debug_assert!(
                type_name(sb) == "BoxShape",
                "narrowphase dispatch: shape_b must be BoxShape for this handler (registration/order invariant)"
            );
            match specialized::sphere_box(s, ta, b, tb) {
                Some(contact) => {
                    let mut m = ContactManifold::new(pair);
                    m.add_contact(contact);
                    NarrowPhaseResult::contact(m)
                }
                None => NarrowPhaseResult::separated(),
            }
        });
        d.register_pair(ShapeType::Box, ShapeType::Box, |sa, ta, sb, tb, pair| {
            // Recover both `BoxShape`s through safe `Any` downcasts. Registered
            // under the (Box, Box) key; a mismatched direct call downcasts to
            // `None` and falls back to GJK/EPA rather than triggering UB.
            let (Some(b1), Some(b2)) = (
                sa.as_any().downcast_ref::<BoxShape>(),
                sb.as_any().downcast_ref::<BoxShape>(),
            ) else {
                return gjk_fallback_dispatch(sa, ta, sb, tb, pair);
            };
            debug_assert!(
                type_name(sa) == "BoxShape",
                "narrowphase dispatch: shape_a must be BoxShape for this handler (registration/order invariant)"
            );
            debug_assert!(
                type_name(sb) == "BoxShape",
                "narrowphase dispatch: shape_b must be BoxShape for this handler (registration/order invariant)"
            );
            box_box_manifold(b1, ta, b2, tb, pair)
        });
        d.register_pair(
            ShapeType::Capsule,
            ShapeType::Capsule,
            |sa, ta, sb, tb, pair| {
                // Recover both `Capsule`s through safe `Any` downcasts. Registered
                // under the (Capsule, Capsule) key; a mismatched direct call
                // downcasts to `None` and falls back to GJK/EPA rather than UB.
                let (Some(c1), Some(c2)) = (
                    sa.as_any().downcast_ref::<Capsule>(),
                    sb.as_any().downcast_ref::<Capsule>(),
                ) else {
                    return gjk_fallback_dispatch(sa, ta, sb, tb, pair);
                };
                debug_assert!(
                    type_name(sa) == "Capsule",
                    "narrowphase dispatch: shape_a must be Capsule for this handler (registration/order invariant)"
                );
                debug_assert!(
                    type_name(sb) == "Capsule",
                    "narrowphase dispatch: shape_b must be Capsule for this handler (registration/order invariant)"
                );
                match specialized::capsule_capsule(c1, ta, c2, tb) {
                    Some(contact) => {
                        let mut m = ContactManifold::new(pair);
                        m.add_contact(contact);
                        NarrowPhaseResult::contact(m)
                    }
                    None => NarrowPhaseResult::separated(),
                }
            },
        );
        d.register_pair(
            ShapeType::Sphere,
            ShapeType::Capsule,
            sphere_capsule_dispatch,
        );
        d.register_pair(ShapeType::Box, ShapeType::Capsule, gjk_fallback_dispatch);
        d
    }
}
