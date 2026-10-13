// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Tests: generational body and collider handles.
//!
//! Covers the ABA problem of plain index handles (a stale handle silently
//! controlling whatever body later reuses its slot), recycling of collider
//! slots when a body is removed, and the slot-capacity limit.

#![cfg(test)]

use super::handles::{
    HANDLE_GENERATION_MASK, HANDLE_INDEX_BITS, HANDLE_INDEX_MASK, decode_handle, encode_handle,
    next_generation,
};
use super::*;

/// An engine without gravity, so that bodies stay where they are put.
fn quiet_engine() -> WasmPhysicsEngine {
    WasmPhysicsEngine::new(0.0, 0.0, 0.0)
}

// ---------------------------------------------------------------------------
// Handle encoding
// ---------------------------------------------------------------------------

#[test]
fn handle_layout_is_20_bit_index_and_12_bit_generation() {
    assert_eq!(HANDLE_INDEX_BITS, 20);
    assert_eq!(HANDLE_INDEX_MASK, 0x000F_FFFF);
    assert_eq!(HANDLE_GENERATION_MASK, 0x0000_0FFF);
    assert_eq!(MAX_SLOTS, 1_048_575);
}

#[test]
fn encode_decode_roundtrip() {
    for &(index, generation) in &[
        (0usize, 0u32),
        (1, 0),
        (17, 3),
        (1_000, 4_095),
        (MAX_SLOTS - 1, 1),
    ] {
        let handle = encode_handle(index, generation);
        assert_eq!(decode_handle(handle), (index, generation));
    }
}

#[test]
fn generation_zero_handles_equal_their_slot_index() {
    for index in [0usize, 1, 2, 999] {
        assert_eq!(encode_handle(index, 0), index as u32);
    }
}

#[test]
fn generation_wraps_at_twelve_bits() {
    assert_eq!(next_generation(0), 1);
    assert_eq!(next_generation(4_094), 4_095);
    assert_eq!(next_generation(4_095), 0);
    // A generation that does not fit is reduced, never leaking into the index.
    assert_eq!(decode_handle(encode_handle(7, 4_096)), (7, 0));
}

#[test]
fn no_allocatable_handle_equals_the_world_sentinel() {
    // `u32::MAX` is used throughout the crate as "no body"; the largest
    // allocatable slot index is MAX_SLOTS - 1, so even the largest generation
    // can never produce it.
    let largest = encode_handle(MAX_SLOTS - 1, HANDLE_GENERATION_MASK);
    assert_ne!(largest, u32::MAX);
    assert_eq!(decode_handle(u32::MAX).0, MAX_SLOTS);
}

#[test]
fn js_handle_helpers_decode_index_and_generation() {
    let handle = encode_handle(42, 5);
    assert_eq!(WasmPhysicsEngine::handle_index_js(handle), 42);
    assert_eq!(WasmPhysicsEngine::handle_generation_js(handle), 5);
}

// ---------------------------------------------------------------------------
// Stale body handles (F5653)
// ---------------------------------------------------------------------------

#[test]
fn fresh_engine_hands_out_sequential_handles() {
    let mut engine = quiet_engine();
    let a = engine.add_dynamic_body(1.0, 0.0, 0.0, 0.0).expect("body a");
    let b = engine.add_dynamic_body(1.0, 1.0, 0.0, 0.0).expect("body b");
    let c = engine.add_static_body(2.0, 0.0, 0.0).expect("body c");
    assert_eq!((a, b, c), (0, 1, 2));
    assert_eq!(engine.get_all_body_handles(), vec![0, 1, 2]);
}

#[test]
fn recycled_slot_gets_a_new_generation() {
    let mut engine = quiet_engine();
    let old = engine
        .add_dynamic_body(1.0, 0.0, 0.0, 0.0)
        .expect("old body");
    engine.remove_body(old).expect("remove old body");
    let new = engine
        .add_dynamic_body(1.0, 5.0, 5.0, 5.0)
        .expect("new body");

    assert_ne!(
        old, new,
        "a recycled slot must not hand out the same handle"
    );
    assert_eq!(decode_handle(old), (0, 0));
    assert_eq!(decode_handle(new), (0, 1));
    assert_eq!(engine.get_position(new), [5.0, 5.0, 5.0]);
}

#[test]
fn stale_handle_cannot_control_the_body_that_reuses_its_slot() {
    let mut engine = quiet_engine();
    let old = engine
        .add_dynamic_body(1.0, 0.0, 0.0, 0.0)
        .expect("old body");
    engine.remove_body(old).expect("remove old body");
    let new = engine
        .add_dynamic_body(2.0, 5.0, 0.0, 0.0)
        .expect("new body");

    assert_eq!(
        engine.set_velocity(old, 9.0, 0.0, 0.0),
        Err(Error::InvalidHandle(old))
    );
    assert_eq!(
        engine.set_angular_velocity(old, 1.0, 1.0, 1.0),
        Err(Error::InvalidHandle(old))
    );
    assert_eq!(
        engine.set_position(old, 1.0, 2.0, 3.0),
        Err(Error::InvalidHandle(old))
    );
    assert!(engine.apply_force(old, 1.0, 0.0, 0.0).is_err());
    assert!(engine.apply_torque(old, 1.0, 0.0, 0.0).is_err());
    assert!(engine.apply_impulse(old, 1.0, 0.0, 0.0).is_err());
    assert!(
        engine
            .apply_impulse_at_point(old, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0)
            .is_err()
    );
    assert!(engine.set_body_linear_damping(old, 0.5).is_err());
    assert!(engine.set_body_angular_damping(old, 0.5).is_err());

    // The new body is untouched by every one of the rejected calls.
    assert_eq!(engine.get_velocity(new), [0.0; 3]);
    assert_eq!(engine.get_angular_velocity(new), [0.0; 3]);
    assert_eq!(engine.get_position(new), [5.0, 0.0, 0.0]);
    let state = engine.get_body_state(new).expect("new body state");
    assert!(!state.is_sleeping);
    assert_eq!(state.kinetic_energy, 0.0);
}

#[test]
fn stale_handle_cannot_remove_the_body_that_reuses_its_slot() {
    let mut engine = quiet_engine();
    let old = engine
        .add_dynamic_body(1.0, 0.0, 0.0, 0.0)
        .expect("old body");
    engine.remove_body(old).expect("remove old body");
    let new = engine
        .add_dynamic_body(1.0, 0.0, 0.0, 0.0)
        .expect("new body");

    assert_eq!(engine.remove_body(old), Err(Error::InvalidHandle(old)));
    assert_eq!(engine.get_body_count(), 1);
    assert!(engine.contains_body(new));
    assert!(!engine.contains_body(old));
}

#[test]
fn stale_handle_queries_return_neutral_values() {
    let mut engine = quiet_engine();
    let old = engine
        .add_dynamic_body(1.0, 3.0, 4.0, 5.0)
        .expect("old body");
    engine.remove_body(old).expect("remove old body");
    let _new = engine
        .add_dynamic_body(1.0, 7.0, 8.0, 9.0)
        .expect("new body");

    assert_eq!(engine.get_position(old), [0.0; 3]);
    assert_eq!(engine.get_velocity(old), [0.0; 3]);
    assert_eq!(engine.get_rotation(old), [0.0, 0.0, 0.0, 1.0]);
    assert!(engine.get_body_state(old).is_none());
    assert_eq!(engine.body_inv_mass(old), 0.0);
    assert!(!engine.body_is_dynamic(old));
}

#[test]
fn removed_body_without_reuse_is_rejected() {
    let mut engine = quiet_engine();
    let handle = engine.add_dynamic_body(1.0, 0.0, 0.0, 0.0).expect("body");
    engine.remove_body(handle).expect("remove");
    assert_eq!(
        engine.remove_body(handle),
        Err(Error::InvalidHandle(handle))
    );
    assert_eq!(
        engine.set_velocity(handle, 1.0, 0.0, 0.0),
        Err(Error::InvalidHandle(handle))
    );
    assert_eq!(engine.get_body_count(), 0);
}

#[test]
fn handles_of_live_bodies_stay_valid_when_other_bodies_are_recycled() {
    let mut engine = quiet_engine();
    let keep_a = engine.add_dynamic_body(1.0, 0.0, 0.0, 0.0).expect("keep a");
    let churn = engine.add_dynamic_body(1.0, 1.0, 0.0, 0.0).expect("churn");
    let keep_b = engine.add_dynamic_body(1.0, 2.0, 0.0, 0.0).expect("keep b");
    for round in 0..10 {
        engine.remove_body(churn).ok();
        let again = engine
            .add_dynamic_body(1.0, 10.0 + round as f64, 0.0, 0.0)
            .expect("churn again");
        engine.remove_body(again).expect("remove churn again");
    }
    assert_eq!(engine.get_position(keep_a), [0.0, 0.0, 0.0]);
    assert_eq!(engine.get_position(keep_b), [2.0, 0.0, 0.0]);
    assert_eq!(engine.get_body_count(), 2);
}

#[test]
fn reset_invalidates_handles_obtained_before_it() {
    let mut engine = quiet_engine();
    let before = engine.add_dynamic_body(1.0, 0.0, 0.0, 0.0).expect("before");
    let collider = engine.add_sphere_collider(before, 0.5).expect("collider");
    engine.reset();
    assert_eq!(engine.get_body_count(), 0);
    assert_eq!(engine.get_active_collider_count(), 0);

    let after = engine.add_dynamic_body(1.0, 1.0, 1.0, 1.0).expect("after");
    assert_ne!(before, after);
    assert_eq!(
        engine.set_velocity(before, 5.0, 0.0, 0.0),
        Err(Error::InvalidHandle(before))
    );
    assert_eq!(engine.get_velocity(after), [0.0; 3]);
    assert_eq!(
        engine.remove_collider(collider),
        Err(Error::InvalidHandle(collider))
    );
}

#[test]
fn engine_keeps_working_after_thousands_of_recycles() {
    let mut engine = quiet_engine();
    let mut handles = Vec::new();
    for _ in 0..5_000 {
        let handle = engine.add_dynamic_body(1.0, 0.0, 0.0, 0.0).expect("body");
        handles.push(handle);
        engine.remove_body(handle).expect("remove");
    }
    // Only one slot was ever needed.
    assert_eq!(engine.bodies.len(), 1);
    // The generation wrapped around 4096 recycles; consecutive handles still
    // differ, and the live one is the last handle handed out.
    assert!(handles.windows(2).all(|pair| pair[0] != pair[1]));
    let live = engine.add_dynamic_body(1.0, 0.0, 0.0, 0.0).expect("live");
    assert!(engine.contains_body(live));
}

// ---------------------------------------------------------------------------
// Collider handles and free list (F5653)
// ---------------------------------------------------------------------------

#[test]
fn removing_a_body_recycles_its_collider_slots() {
    let mut engine = quiet_engine();
    let body = engine.add_dynamic_body(1.0, 0.0, 0.0, 0.0).expect("body");
    engine.add_sphere_collider(body, 0.5).expect("collider 0");
    engine
        .add_box_collider(body, 0.5, 0.5, 0.5)
        .expect("collider 1");
    engine
        .add_capsule_collider(body, 0.25, 1.0)
        .expect("collider 2");
    assert_eq!(engine.get_active_collider_count(), 3);
    assert_eq!(engine.colliders.len(), 3);

    engine.remove_body(body).expect("remove body");
    assert_eq!(engine.get_active_collider_count(), 0);
    assert_eq!(engine.collider_free_list.len(), 3, "slots must be recycled");

    // A new body with three colliders re-uses the slots instead of growing the
    // table.
    let next = engine
        .add_dynamic_body(1.0, 0.0, 0.0, 0.0)
        .expect("next body");
    let handles = [
        engine.add_sphere_collider(next, 0.5).expect("collider 0"),
        engine
            .add_box_collider(next, 0.5, 0.5, 0.5)
            .expect("collider 1"),
        engine
            .add_capsule_collider(next, 0.25, 1.0)
            .expect("collider 2"),
    ];
    assert_eq!(
        engine.colliders.len(),
        3,
        "the collider table must not grow"
    );
    assert!(engine.collider_free_list.is_empty());
    for handle in handles {
        assert_eq!(
            decode_handle(handle).1,
            1,
            "recycled collider slot, generation 1"
        );
    }
}

#[test]
fn collider_handles_are_generational() {
    let mut engine = quiet_engine();
    let body = engine.add_dynamic_body(1.0, 0.0, 0.0, 0.0).expect("body");
    let first = engine
        .add_sphere_collider(body, 0.5)
        .expect("first collider");
    assert_eq!(first, 0);
    engine.remove_collider(first).expect("remove collider");
    assert_eq!(engine.get_active_collider_count(), 0);

    let second = engine
        .add_sphere_collider(body, 0.75)
        .expect("second collider");
    assert_ne!(first, second);
    assert_eq!(decode_handle(second), (0, 1));
    assert_eq!(
        engine.remove_collider(first),
        Err(Error::InvalidHandle(first)),
        "a stale collider handle must not remove the collider that reuses its slot"
    );
    assert_eq!(engine.get_active_collider_count(), 1);
    engine
        .remove_collider(second)
        .expect("remove second collider");
    assert_eq!(engine.get_active_collider_count(), 0);
}

#[test]
fn remove_collider_with_unknown_handle_fails() {
    let mut engine = quiet_engine();
    assert_eq!(engine.remove_collider(0), Err(Error::InvalidHandle(0)));
    assert_eq!(engine.remove_collider(123), Err(Error::InvalidHandle(123)));
    assert_eq!(
        engine.remove_collider(u32::MAX),
        Err(Error::InvalidHandle(u32::MAX))
    );
}

#[test]
fn colliders_of_a_removed_body_do_not_follow_the_new_occupant_of_its_slot() {
    let mut engine = WasmPhysicsEngine::new(0.0, -10.0, 0.0);
    let ground = engine.add_static_body(0.0, 0.0, 0.0).expect("ground");
    engine
        .add_plane_collider(ground, 0.0, 1.0, 0.0, 0.0)
        .expect("ground plane");

    let old = engine
        .add_dynamic_body(1.0, 0.0, 5.0, 0.0)
        .expect("old body");
    engine.add_sphere_collider(old, 0.5).expect("old collider");
    engine.remove_body(old).expect("remove old body");

    // Same slot, but no collider: it must fall through the ground plane.
    let new = engine
        .add_dynamic_body(1.0, 0.0, 0.2, 0.0)
        .expect("new body");
    assert_eq!(decode_handle(new).0, decode_handle(old).0);
    assert_eq!(
        engine.get_active_collider_count(),
        1,
        "only the plane remains"
    );
    for _ in 0..30 {
        engine.step(1.0 / 60.0).expect("step");
    }
    assert!(
        engine.get_position(new)[1] < 0.0,
        "a body without colliders falls through the plane, got y = {}",
        engine.get_position(new)[1]
    );
    assert_eq!(engine.get_contact_count(), 0);
}

// ---------------------------------------------------------------------------
// Capacity (F5653: beyond the slot limit return CapacityExceeded)
// ---------------------------------------------------------------------------

#[test]
fn body_table_reports_capacity_exceeded() {
    let mut engine = quiet_engine();
    engine.slot_limit = 2;
    let a = engine.add_dynamic_body(1.0, 0.0, 0.0, 0.0).expect("body a");
    engine.add_dynamic_body(1.0, 1.0, 0.0, 0.0).expect("body b");
    assert_eq!(
        engine.add_dynamic_body(1.0, 2.0, 0.0, 0.0),
        Err(Error::CapacityExceeded { max: 2 })
    );
    assert_eq!(
        engine.add_static_body(2.0, 0.0, 0.0),
        Err(Error::CapacityExceeded { max: 2 })
    );
    assert_eq!(engine.get_body_count(), 2);

    // Freeing a slot makes room again.
    engine.remove_body(a).expect("remove body a");
    engine
        .add_dynamic_body(1.0, 2.0, 0.0, 0.0)
        .expect("recycled slot is available");
    assert_eq!(engine.get_body_count(), 2);
}

#[test]
fn collider_table_reports_capacity_exceeded() {
    let mut engine = quiet_engine();
    engine.slot_limit = 2;
    let body = engine.add_dynamic_body(1.0, 0.0, 0.0, 0.0).expect("body");
    let first = engine.add_sphere_collider(body, 0.5).expect("collider a");
    engine.add_sphere_collider(body, 0.25).expect("collider b");
    assert_eq!(
        engine.add_sphere_collider(body, 0.1),
        Err(Error::CapacityExceeded { max: 2 })
    );
    assert_eq!(engine.get_active_collider_count(), 2);

    engine.remove_collider(first).expect("remove collider a");
    engine
        .add_sphere_collider(body, 0.1)
        .expect("recycled collider slot is available");
}

#[test]
fn add_body_with_collider_rolls_back_the_body_when_the_collider_does_not_fit() {
    let mut engine = quiet_engine();
    engine.slot_limit = 2;
    let host = engine
        .add_dynamic_body(1.0, 0.0, 0.0, 0.0)
        .expect("host body");
    engine.add_sphere_collider(host, 0.5).expect("collider a");
    engine.add_sphere_collider(host, 0.25).expect("collider b");

    let body = RigidBodyConfig::new().with_position(5.0, 0.0, 0.0);
    let result = engine.add_body_with_collider(&body, &ColliderConfig::sphere(0.5));
    assert_eq!(result, Err(Error::CapacityExceeded { max: 2 }));
    assert_eq!(
        engine.get_body_count(),
        1,
        "the body must have been rolled back"
    );
    assert_eq!(engine.get_all_body_handles(), vec![host]);
    assert_eq!(engine.get_active_collider_count(), 2);
}

// ---------------------------------------------------------------------------
// JS wrappers that can run natively (success paths do not build a JsValue)
// ---------------------------------------------------------------------------

#[test]
fn js_wrappers_create_and_remove_bodies() {
    let mut engine = quiet_engine();
    let handle = engine
        .add_dynamic_body_js(1.0, 1.0, 2.0, 3.0)
        .expect("dynamic body");
    let ground = engine
        .add_static_body_js(0.0, 0.0, 0.0)
        .expect("static body");
    engine
        .add_sphere_collider_js(handle, 0.5)
        .expect("sphere collider");
    engine
        .add_box_collider_js(handle, 0.5, 0.5, 0.5)
        .expect("box collider");
    engine
        .add_capsule_collider_js(handle, 0.25, 1.0)
        .expect("capsule collider");
    engine
        .add_cylinder_collider_js(handle, 0.25, 1.0)
        .expect("cylinder collider");
    engine
        .add_plane_collider_js(ground, 0.0, 1.0, 0.0, 0.0)
        .expect("plane collider");
    assert_eq!(engine.get_body_count_js(), 2);
    assert!(engine.contains_body_js(handle));
    assert_eq!(engine.get_position_js(handle), vec![1.0, 2.0, 3.0]);

    assert!(engine.remove_body_js(handle));
    assert!(
        !engine.remove_body_js(handle),
        "second removal reports false"
    );
    assert!(!engine.contains_body_js(handle));
    assert!(!engine.set_velocity_js(handle, 1.0, 0.0, 0.0));
    assert!(!engine.apply_force_js(handle, 1.0, 0.0, 0.0));
    assert!(!engine.apply_torque_js(handle, 1.0, 0.0, 0.0));
    assert!(!engine.apply_impulse_js(handle, 1.0, 0.0, 0.0));
    assert!(!engine.apply_impulse_at_point_js(handle, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0));
    assert!(!engine.set_angular_velocity_js(handle, 1.0, 0.0, 0.0));
    assert!(!engine.set_position_js(handle, 1.0, 0.0, 0.0));
}
