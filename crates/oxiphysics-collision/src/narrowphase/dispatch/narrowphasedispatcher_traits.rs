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
use oxiphysics_geometry::{BoxShape, Capsule, Shape, Sphere};

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
                // SAFETY: this handler is registered in the dispatch table under the
                // (Sphere, Sphere) key, so the dispatcher only invokes it when both
                // `sa` and `sb` are concrete `Sphere` values (the caller pairs each
                // shape with its matching `ShapeType`). The reborrow keeps `sa`/`sb`'s
                // lifetime, and `Sphere` has no stricter alignment than `dyn Shape`.
                debug_assert!(
                    type_name(sa) == "Sphere",
                    "narrowphase dispatch: shape_a must be Sphere for this handler (registration/order invariant)"
                );
                debug_assert!(
                    type_name(sb) == "Sphere",
                    "narrowphase dispatch: shape_b must be Sphere for this handler (registration/order invariant)"
                );
                let s1 = unsafe { &*(sa as *const dyn Shape as *const Sphere) };
                let s2 = unsafe { &*(sb as *const dyn Shape as *const Sphere) };
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
            // SAFETY: registered under the (Sphere, Box) key, so the dispatcher only
            // invokes this handler when `sa` is a concrete `Sphere` and `sb` a
            // concrete `BoxShape` — the caller pairs each shape with its matching
            // `ShapeType`, and for this mixed pair the arguments arrive in the same
            // positional order they were registered (Sphere first, Box second). The
            // reborrows keep the input lifetimes; both targets are plain structs with
            // alignment no stricter than `dyn Shape`.
            debug_assert!(
                type_name(sa) == "Sphere",
                "narrowphase dispatch: shape_a must be Sphere for this handler (registration/order invariant)"
            );
            debug_assert!(
                type_name(sb) == "BoxShape",
                "narrowphase dispatch: shape_b must be BoxShape for this handler (registration/order invariant)"
            );
            let s = unsafe { &*(sa as *const dyn Shape as *const Sphere) };
            let b = unsafe { &*(sb as *const dyn Shape as *const BoxShape) };
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
            // SAFETY: registered under the (Box, Box) key, so the dispatcher only
            // invokes this handler when both `sa` and `sb` are concrete `BoxShape`
            // values (the caller pairs each shape with its matching `ShapeType`). The
            // reborrows keep the input lifetimes; `BoxShape` has no stricter alignment
            // than `dyn Shape`.
            debug_assert!(
                type_name(sa) == "BoxShape",
                "narrowphase dispatch: shape_a must be BoxShape for this handler (registration/order invariant)"
            );
            debug_assert!(
                type_name(sb) == "BoxShape",
                "narrowphase dispatch: shape_b must be BoxShape for this handler (registration/order invariant)"
            );
            let b1 = unsafe { &*(sa as *const dyn Shape as *const BoxShape) };
            let b2 = unsafe { &*(sb as *const dyn Shape as *const BoxShape) };
            match specialized::box_box_sat(b1, ta, b2, tb) {
                Some(contact) => {
                    let mut m = ContactManifold::new(pair);
                    m.add_contact(contact);
                    NarrowPhaseResult::contact(m)
                }
                None => NarrowPhaseResult::separated(),
            }
        });
        d.register_pair(
            ShapeType::Capsule,
            ShapeType::Capsule,
            |sa, ta, sb, tb, pair| {
                // SAFETY: registered under the (Capsule, Capsule) key, so the
                // dispatcher only invokes this handler when both `sa` and `sb` are
                // concrete `Capsule` values (the caller pairs each shape with its
                // matching `ShapeType`). The reborrows keep the input lifetimes;
                // `Capsule` has no stricter alignment than `dyn Shape`.
                debug_assert!(
                    type_name(sa) == "Capsule",
                    "narrowphase dispatch: shape_a must be Capsule for this handler (registration/order invariant)"
                );
                debug_assert!(
                    type_name(sb) == "Capsule",
                    "narrowphase dispatch: shape_b must be Capsule for this handler (registration/order invariant)"
                );
                let c1 = unsafe { &*(sa as *const dyn Shape as *const Capsule) };
                let c2 = unsafe { &*(sb as *const dyn Shape as *const Capsule) };
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
