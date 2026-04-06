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

#[allow(unused_imports)]
use super::functions::*;
use super::functions::{gjk_fallback_dispatch, sphere_capsule_dispatch};
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
