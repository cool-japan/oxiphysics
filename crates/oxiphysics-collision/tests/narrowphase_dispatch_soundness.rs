// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Soundness and warm-start tests for the narrow-phase dispatcher.
//!
//! Covers the v0.2.0 dispatch-soundness fix (safe `Any` downcasts + argument
//! canonicalisation, no raw-pointer transmutes) and the default-on GJK
//! warm-start cache wiring. All assertions go through the public API.

use oxiphysics_collision::narrowphase::{DispatchConfig, ShapeType};
use oxiphysics_collision::{CollisionPair, NarrowPhaseDispatcher};
use oxiphysics_core::Transform;
use oxiphysics_core::math::Vec3;
use oxiphysics_geometry::{BoxShape, Capsule, Shape, Sphere};

fn at(x: f64) -> Transform {
    Transform::from_position(Vec3::new(x, 0.0, 0.0))
}

/// Dispatch a pair in both argument orders and assert the results are mirror
/// images: opposite normals and (optionally) swapped witness points.
fn assert_order_independent(
    dispatcher: &NarrowPhaseDispatcher,
    shape_a: &dyn Shape,
    type_a: ShapeType,
    transform_a: &Transform,
    shape_b: &dyn Shape,
    type_b: ShapeType,
    transform_b: &Transform,
    strict_points: bool,
    label: &str,
) {
    let pair = CollisionPair::new(0, 1);
    let ab = dispatcher.dispatch(
        shape_a,
        type_a,
        transform_a,
        shape_b,
        type_b,
        transform_b,
        pair,
    );
    let ba = dispatcher.dispatch(
        shape_b,
        type_b,
        transform_b,
        shape_a,
        type_a,
        transform_a,
        pair,
    );

    assert_eq!(
        ab.has_contact(),
        ba.has_contact(),
        "{label}: contact presence must match in both argument orders"
    );

    if let (Some(mab), Some(mba)) = (ab.manifold, ba.manifold) {
        assert_eq!(
            mab.contacts.len(),
            mba.contacts.len(),
            "{label}: contact count must match in both orders"
        );
        let cab = &mab.contacts[0];
        let cba = &mba.contacts[0];

        let normal_sum = (cab.normal + cba.normal).norm();
        assert!(
            normal_sum < 1e-6,
            "{label}: normal sign must flip with argument order ({:?} vs {:?})",
            cab.normal,
            cba.normal
        );

        if strict_points {
            assert!(
                (cab.point_a - cba.point_b).norm() < 1e-6,
                "{label}: point_a(A,B) must equal point_b(B,A)"
            );
            assert!(
                (cab.point_b - cba.point_a).norm() < 1e-6,
                "{label}: point_b(A,B) must equal point_a(B,A)"
            );
        }
    }
}

#[test]
fn every_registered_pair_is_order_independent() {
    let d = NarrowPhaseDispatcher::default();
    let sphere = Sphere::new(1.0);
    let sphere2 = Sphere::new(1.0);
    let box_a = BoxShape::new(Vec3::new(1.0, 1.0, 1.0));
    let box_b = BoxShape::new(Vec3::new(1.0, 1.0, 1.0));
    let cap_a = Capsule::new(0.5, 1.0);
    let cap_b = Capsule::new(0.5, 1.0);

    // Equal-ordinal symmetric pairs: normal must flip; points need not strictly swap.
    assert_order_independent(
        &d,
        &sphere,
        ShapeType::Sphere,
        &at(0.0),
        &sphere2,
        ShapeType::Sphere,
        &at(1.5),
        true,
        "sphere-sphere",
    );
    assert_order_independent(
        &d,
        &box_a,
        ShapeType::Box,
        &at(0.0),
        &box_b,
        ShapeType::Box,
        &at(1.5),
        false,
        "box-box",
    );
    assert_order_independent(
        &d,
        &cap_a,
        ShapeType::Capsule,
        &at(0.0),
        &cap_b,
        ShapeType::Capsule,
        &at(0.5),
        false,
        "capsule-capsule",
    );

    // Different-ordinal pairs exercise the canonicalise-swap-and-flip path: the
    // reversed result must be an exact mirror (points swapped too).
    assert_order_independent(
        &d,
        &sphere,
        ShapeType::Sphere,
        &at(0.0),
        &box_b,
        ShapeType::Box,
        &at(1.5),
        true,
        "sphere-box",
    );
    assert_order_independent(
        &d,
        &sphere,
        ShapeType::Sphere,
        &at(0.0),
        &cap_b,
        ShapeType::Capsule,
        &at(1.0),
        true,
        "sphere-capsule",
    );
    assert_order_independent(
        &d,
        &box_a,
        ShapeType::Box,
        &at(0.0),
        &cap_b,
        ShapeType::Capsule,
        &at(1.0),
        true,
        "box-capsule (gjk fallback)",
    );
}

#[test]
fn reversed_order_does_not_panic_in_release_or_debug() {
    // Before the fix, a non-canonical pair invoked a registered handler with
    // mismatched concrete types, panicking (debug) or hitting UB (release).
    // Now it is canonicalised first, so both orders are safe.
    let d = NarrowPhaseDispatcher::default();
    let sphere = Sphere::new(1.0);
    let box_s = BoxShape::new(Vec3::new(1.0, 1.0, 1.0));
    let pair = CollisionPair::new(0, 1);

    // (Box, Sphere) — only (Sphere, Box) is registered.
    let r = d.dispatch(
        &box_s,
        ShapeType::Box,
        &at(0.0),
        &sphere,
        ShapeType::Sphere,
        &at(1.5),
        pair,
    );
    assert!(
        r.has_contact(),
        "overlapping box-sphere must produce contact"
    );
}

#[test]
fn wrong_concrete_type_falls_back_to_gjk_without_panic() {
    // Feed the (Sphere, Sphere) handler two boxes while claiming they are
    // spheres. The handler's `downcast_ref::<Sphere>()` returns `None`, so it
    // must fall back to the safe GJK/EPA path instead of transmuting/panicking.
    let d = NarrowPhaseDispatcher::default();
    let box_a = BoxShape::new(Vec3::new(1.0, 1.0, 1.0));
    let box_b = BoxShape::new(Vec3::new(1.0, 1.0, 1.0));
    let pair = CollisionPair::new(0, 1);

    let overlapping = d.dispatch(
        &box_a,
        ShapeType::Sphere, // deliberately wrong type tag
        &at(0.0),
        &box_b,
        ShapeType::Sphere,
        &at(1.0),
        pair,
    );
    assert!(
        overlapping.has_contact(),
        "GJK fallback must still detect the overlap of two boxes"
    );

    let separated = d.dispatch(
        &box_a,
        ShapeType::Sphere,
        &at(0.0),
        &box_b,
        ShapeType::Sphere,
        &at(50.0),
        pair,
    );
    assert!(
        !separated.has_contact(),
        "GJK fallback must report separation for far-apart boxes"
    );
}

#[test]
fn warm_start_enabled_by_default() {
    assert!(
        DispatchConfig::default().enable_warm_start,
        "warm start must be on by default"
    );
}

#[test]
fn warm_start_hit_rate_exceeds_half_on_coherent_motion() {
    // Empty dispatcher -> every pair takes the GJK fallback, which the cached
    // entry point warm-starts from the persistent per-pair simplex cache.
    let mut d = NarrowPhaseDispatcher::empty();
    let s1 = Sphere::new(1.0);
    let s2 = Sphere::new(1.0);
    let pair = CollisionPair::new(0, 1);

    let frames = 60;
    for frame in 0..frames {
        // Slow coherent drift while staying overlapping (centre distance < 2.0).
        let x = 1.0 + 0.0005 * frame as f64;
        let r = d.dispatch_cached(
            &s1,
            ShapeType::Sphere,
            &at(0.0),
            &s2,
            ShapeType::Sphere,
            &at(x),
            pair,
        );
        assert!(
            r.has_contact(),
            "overlapping spheres must contact at frame {frame}"
        );
    }

    let stats = d.gjk_cache_stats();
    assert!(
        stats.queries >= frames as u64,
        "expected at least {frames} cached queries, got {}",
        stats.queries
    );
    assert!(
        stats.hit_ratio() > 0.5,
        "warm-start hit ratio {} must exceed 0.5 on coherent motion",
        stats.hit_ratio()
    );
}

#[test]
fn cache_is_invalidated_on_body_removal() {
    let mut d = NarrowPhaseDispatcher::empty();
    let s1 = Sphere::new(1.0);
    let s2 = Sphere::new(1.0);
    let pair = CollisionPair::new(7, 9);

    for _ in 0..5 {
        let _ = d.dispatch_cached(
            &s1,
            ShapeType::Sphere,
            &at(0.0),
            &s2,
            ShapeType::Sphere,
            &at(1.0),
            pair,
        );
    }
    assert_eq!(
        d.gjk_cache_pair_count(),
        1,
        "the moving pair should have a cached simplex"
    );

    d.remove_body(7);
    assert_eq!(
        d.gjk_cache_pair_count(),
        0,
        "removing body 7 must drop every cached pair that references it"
    );
}
