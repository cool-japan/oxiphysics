// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0
//! Integration tests for the TGS-Soft (Catto / Box2D v3) contact response path
//! of [`PhysicsWorld`].
//!
//! The contact normal points from B toward A and every body in this simplified
//! world is a sphere of radius `mass.cbrt() * 0.5` (min `0.1`). A sphere-sphere
//! contact lies on the line of centres, so the angular Jacobian terms vanish and
//! these scenes exercise the linear soft normal response.
//!
//! Soft-off behaviour is byte-identical to the legacy paths; the soft path is
//! opt-in via `solver_config.use_soft_contacts`.

use oxiphysics_core::BodyHandle;
use oxiphysics_rigid::world::PhysicsWorld;
use oxiphysics_rigid::{BodyType, RigidBody};

/// Spawn an undamped dynamic sphere of `mass` (radius `mass.cbrt() * 0.5`) at
/// `(x, y, z)` and return its handle.
fn spawn_sphere(world: &mut PhysicsWorld, mass: f64, x: f64, y: f64, z: f64) -> BodyHandle {
    let mut body = RigidBody::new(mass);
    body.linear_damping = 0.0;
    body.angular_damping = 0.0;
    let handle = world.add_rigid_body(body);
    let b = world.get_body_mut(handle).expect("body just added");
    b.transform.position.x = x;
    b.transform.position.y = y;
    b.transform.position.z = z;
    handle
}

/// Spawn a static sphere of radius `mass.cbrt() * 0.5` centred at the origin.
/// `.mass` is set after `new_static()` because the broad/narrowphase radius
/// formula reads `body.mass.cbrt()`.
fn spawn_static_ground(world: &mut PhysicsWorld, mass: f64) -> BodyHandle {
    let mut ground = RigidBody::new_static();
    ground.mass = mass;
    world.add_rigid_body(ground)
}

/// Disable the sleeping system so settling dynamics stay observable (thresholds
/// of 0 make the `< threshold` sleep test always fail).
fn disable_sleep(world: &mut PhysicsWorld) {
    world.linear_sleep_threshold = 0.0;
    world.angular_sleep_threshold = 0.0;
}

/// Collect the full kinematic state (position + velocity bits) of `handles` for
/// exact (bit-for-bit) comparison.
fn state_bits(world: &PhysicsWorld, handles: &[BodyHandle]) -> Vec<u64> {
    let mut bits = Vec::with_capacity(handles.len() * 6);
    for &h in handles {
        let b = world.get_body(h).expect("body present");
        bits.push(b.transform.position.x.to_bits());
        bits.push(b.transform.position.y.to_bits());
        bits.push(b.transform.position.z.to_bits());
        bits.push(b.velocity.x.to_bits());
        bits.push(b.velocity.y.to_bits());
        bits.push(b.velocity.z.to_bits());
    }
    bits
}

/// Build a falling-sphere-on-ground scene; returns the dynamic-body handles.
fn build_scene(world: &mut PhysicsWorld) -> Vec<BodyHandle> {
    spawn_static_ground(world, 8000.0); // radius 10
    let mut handles = Vec::new();
    for &y in &[10.5_f64, 11.5, 12.5] {
        handles.push(spawn_sphere(world, 1.0, 0.0, y, 0.0));
    }
    handles
}

// ───────────────────────── Test 1: soft-off byte-identical ─────────────────────

/// With `use_soft_contacts == false` the solver must take the legacy path
/// regardless of the soft tuning fields, producing bit-for-bit identical state.
/// Verified for both `step` and `step_small_steps`.
#[test]
fn soft_off_is_byte_identical() {
    // ── step() path ──
    let mut reference = PhysicsWorld::new([0.0, -9.81, 0.0], 1.0 / 60.0);
    let ref_handles = build_scene(&mut reference);

    let mut probe = PhysicsWorld::new([0.0, -9.81, 0.0], 1.0 / 60.0);
    let probe_handles = build_scene(&mut probe);
    // Soft tuning set but DISABLED — must be ignored on the legacy path.
    probe.solver_config.use_soft_contacts = false;
    probe.solver_config.contact_hertz = 30.0;
    probe.solver_config.contact_damping_ratio = 0.5;

    for _ in 0..10 {
        reference.step();
        probe.step();
    }
    assert_eq!(
        state_bits(&reference, &ref_handles),
        state_bits(&probe, &probe_handles),
        "soft-off step() diverged from the legacy path"
    );
    // Sanity: the scene actually evolved (bodies fell).
    assert!(
        reference
            .get_body(ref_handles[0])
            .expect("present")
            .transform
            .position
            .y
            < 10.5,
        "scene did not evolve — test would be vacuous"
    );

    // ── step_small_steps() path ──
    let mut ref_ss = PhysicsWorld::new([0.0, -9.81, 0.0], 1.0 / 60.0);
    let ref_ss_handles = build_scene(&mut ref_ss);
    ref_ss.solver_config.substeps = 4;

    let mut probe_ss = PhysicsWorld::new([0.0, -9.81, 0.0], 1.0 / 60.0);
    let probe_ss_handles = build_scene(&mut probe_ss);
    probe_ss.solver_config.substeps = 4;
    probe_ss.solver_config.use_soft_contacts = false;
    probe_ss.solver_config.contact_hertz = 45.0;
    probe_ss.solver_config.contact_damping_ratio = 0.25;

    for _ in 0..10 {
        ref_ss.step_small_steps();
        probe_ss.step_small_steps();
    }
    assert_eq!(
        state_bits(&ref_ss, &ref_ss_handles),
        state_bits(&probe_ss, &probe_ss_handles),
        "soft-off step_small_steps() diverged from the legacy path"
    );
}

// ───────────────────────── Test 2: soft stack jitter-free ──────────────────────

/// A vertical column of 5 unit spheres on the ground dome settles under soft
/// contacts (30 Hz, ζ = 1.0) with bounded penetration, no NaN, and a
/// non-increasing kinetic-energy tail.
///
/// A vertical column (not a flat grid) is used because the ground is a convex
/// sphere: off-axis bodies on a flat grid roll off the dome (a geometry
/// artefact, see `small_steps.rs`). The column sits on the dome apex where the
/// normal is vertical.
#[test]
fn soft_stack_settles_without_jitter() {
    let mut world = PhysicsWorld::new([0.0, -9.81, 0.0], 1.0 / 60.0);
    disable_sleep(&mut world);
    world.solver_config.use_soft_contacts = true;
    world.solver_config.contact_hertz = 30.0;
    world.solver_config.contact_damping_ratio = 1.0;
    world.solver_config.restitution = 0.0;
    world.solver_config.velocity_iterations = 12;

    spawn_static_ground(&mut world, 8000.0); // radius 10
    let mut handles = Vec::new();
    for &y in &[10.5_f64, 11.5, 12.5, 13.5, 14.5] {
        handles.push(spawn_sphere(&mut world, 1.0, 0.0, y, 0.0));
    }

    let mut max_penetration = 0.0_f64;
    let mut ke_history = Vec::with_capacity(300);
    for _ in 0..300 {
        world.step();
        let step_pen = world
            .last_contacts
            .iter()
            .map(|c| c.depth)
            .fold(0.0_f64, f64::max);
        max_penetration = max_penetration.max(step_pen);
        ke_history.push(world.kinetic_energy());
    }

    // No NaN / tunnelling.
    for h in &handles {
        let b = world.get_body(*h).expect("dynamic body present");
        let p = &b.transform.position;
        assert!(
            p.x.is_finite() && p.y.is_finite() && p.z.is_finite(),
            "non-finite position"
        );
        assert!(p.y > 0.0, "body tunnelled below ground centre: y = {}", p.y);
    }

    // Penetration stays bounded (soft contacts at the per-step cap of 0.25/dt =
    // 15 Hz are springier than rigid; a stacked column loads the lower contacts,
    // so the bound reflects a settled — not zero — overlap).
    assert!(
        max_penetration < 0.05,
        "soft stack penetration too large: {}",
        max_penetration
    );

    // Energy: bounded peak and a non-increasing settled tail (gravity injects a
    // little KE between solves each step, hence the tolerance).
    let max_ke = ke_history.iter().cloned().fold(0.0_f64, f64::max);
    assert!(
        max_ke.is_finite() && max_ke < 500.0,
        "peak KE diverged: {}",
        max_ke
    );
    let n = ke_history.len();
    for i in (n - 40)..n {
        assert!(
            ke_history[i] <= ke_history[i - 1] + 0.05,
            "KE increased in settled tail at frame {}: {} -> {}",
            i,
            ke_history[i - 1],
            ke_history[i]
        );
    }
    // The tail must actually be small (settled), not merely flat.
    assert!(
        *ke_history.last().expect("ran frames") < 0.5,
        "stack never settled: final KE = {}",
        ke_history.last().expect("ran frames")
    );
}

// ───────────────────────── Test 3: dt-invariance ───────────────────────────────

/// Returns the settle time (seconds) of a single sphere dropped 0.3 m onto the
/// ground with soft contacts (hz = 10, ζ = 0.7) at the given `dt`. Settle time
/// is the last sim-time at which the sphere's speed exceeds 0.1 m/s.
fn settle_time_for_dt(dt: f64) -> f64 {
    let mut world = PhysicsWorld::new([0.0, -9.81, 0.0], dt);
    disable_sleep(&mut world);
    world.solver_config.use_soft_contacts = true;
    world.solver_config.contact_hertz = 10.0;
    world.solver_config.contact_damping_ratio = 0.7;
    world.solver_config.restitution = 0.0;
    world.solver_config.velocity_iterations = 12;

    spawn_static_ground(&mut world, 8000.0); // radius 10, rest contact at y = 10.5
    let ball = spawn_sphere(&mut world, 1.0, 0.0, 10.8, 0.0); // 0.3 m above rest

    let duration = 3.0_f64;
    let steps = (duration / dt).round() as usize;
    let mut last_moving_time = 0.0_f64;
    for k in 0..steps {
        world.step();
        let v = &world.get_body(ball).expect("ball present").velocity;
        let speed = (v.x * v.x + v.y * v.y + v.z * v.z).sqrt();
        assert!(speed.is_finite(), "speed diverged at dt = {}", dt);
        if speed > 0.1 {
            last_moving_time = (k + 1) as f64 * dt;
        }
    }
    last_moving_time
}

/// The soft contact parameterization makes settling behaviour (in seconds)
/// roughly independent of the time-step: settle times at dt ∈ {1/60, 1/120,
/// 1/240} agree within a factor of 2.
#[test]
fn soft_contact_is_dt_invariant() {
    let t60 = settle_time_for_dt(1.0 / 60.0);
    let t120 = settle_time_for_dt(1.0 / 120.0);
    let t240 = settle_time_for_dt(1.0 / 240.0);

    for (label, t) in [("1/60", t60), ("1/120", t120), ("1/240", t240)] {
        assert!(t > 0.0, "no settling transient observed at dt = {}", label);
        assert!(t < 3.0, "did not settle within 3 s at dt = {}", label);
    }
    let max_t = t60.max(t120).max(t240);
    let min_t = t60.min(t120).min(t240);
    assert!(
        max_t / min_t < 2.0,
        "settle time not dt-invariant: 1/60={:.3}s 1/120={:.3}s 1/240={:.3}s (ratio {:.2})",
        t60,
        t120,
        t240,
        max_t / min_t
    );
}

// ───────────────────────── Test 4: heavy-on-light mass ratio ───────────────────

/// A 10 kg sphere resting on a 0.1 kg sphere (100:1 ratio) on the ground stays
/// stable under soft contacts: no NaN and bounded penetration over 80 steps.
#[test]
fn soft_contact_heavy_on_light_mass_ratio() {
    let mut world = PhysicsWorld::new([0.0, -9.81, 0.0], 1.0 / 60.0);
    disable_sleep(&mut world);
    world.solver_config.use_soft_contacts = true;
    world.solver_config.contact_hertz = 30.0;
    world.solver_config.contact_damping_ratio = 1.0;
    world.solver_config.restitution = 0.0;
    world.solver_config.velocity_iterations = 16;

    spawn_static_ground(&mut world, 8000.0); // radius 10
    // Light sphere: 0.1 kg => radius 0.1.cbrt()*0.5 ≈ 0.232; rest on ground ≈ 10.232.
    let light = spawn_sphere(&mut world, 0.1, 0.0, 10.45, 0.0);
    // Heavy sphere: 10 kg => radius 10.cbrt()*0.5 ≈ 1.077; rest on light ≈ 11.54.
    let heavy = spawn_sphere(&mut world, 10.0, 0.0, 11.85, 0.0);

    let mut max_penetration = 0.0_f64;
    for _ in 0..80 {
        world.step();
        let step_pen = world
            .last_contacts
            .iter()
            .map(|c| c.depth)
            .fold(0.0_f64, f64::max);
        max_penetration = max_penetration.max(step_pen);
    }

    for (name, h) in [("light", light), ("heavy", heavy)] {
        let b = world.get_body(h).expect("body present");
        let p = &b.transform.position;
        assert!(
            p.x.is_finite() && p.y.is_finite() && p.z.is_finite(),
            "{name} sphere position non-finite"
        );
        assert!(p.y > 0.0, "{name} sphere tunnelled: y = {}", p.y);
        let v = &b.velocity;
        assert!(
            v.x.is_finite() && v.y.is_finite() && v.z.is_finite(),
            "{name} sphere velocity non-finite"
        );
    }
    assert!(
        max_penetration < 0.2,
        "100:1 mass-ratio penetration unbounded: {}",
        max_penetration
    );
}

// ─────────────────── Test 5: bouncing-ball restitution (soft OFF) ──────────────

/// With soft contacts OFF, the legacy small-steps restitution path is preserved:
/// a ball dropped 1.0 m with e = 0.7 rebounds to ≈ e² · H = 0.49 m (within ~15%).
///
/// `step_small_steps` is used (not `step`) because the legacy `step` velocity
/// solver hardcodes restitution = 0.3; the configurable-restitution behaviour
/// lives on the small-steps path, which this test confirms is unchanged.
///
/// `substeps = 1` makes the per-sub-step restitution exactly `e`
/// (`1 − (1 − e)^(1/1) = e`), so the single contact impulse reflects the full
/// coefficient and the rebound matches the analytic `e²` ratio; with more
/// sub-steps the per-sub-step coefficient is smaller and the realised
/// restitution depends on contact persistence.
#[test]
fn rigid_bouncing_ball_restitution_preserved() {
    let mut world = PhysicsWorld::new([0.0, -9.81, 0.0], 1.0 / 60.0);
    world.solver_config.substeps = 1;
    world.solver_config.restitution = 0.7;
    world.solver_config.use_soft_contacts = false; // legacy path

    spawn_static_ground(&mut world, 8000.0); // radius 10, rest at y = 10.5
    // Drop from H = 1.0 m above rest.
    let ball = spawn_sphere(&mut world, 1.0, 0.0, 11.5, 0.0);

    let mut engaged = false;
    let mut prev_vy = 0.0_f64;
    let mut tracking = false;
    let mut peak_y = f64::NEG_INFINITY;
    let mut h_prime = f64::NAN;

    for _ in 0..400 {
        world.step_small_steps();
        let b = world.get_body(ball).expect("ball present");
        let y = b.transform.position.y;
        let vy = b.velocity.y;

        if y < 10.7 {
            engaged = true;
        }
        if engaged && !tracking && prev_vy < 0.0 && vy >= 0.0 {
            tracking = true;
            peak_y = y;
        } else if tracking {
            if vy >= 0.0 {
                if y > peak_y {
                    peak_y = y;
                }
            } else if h_prime.is_nan() {
                h_prime = peak_y - 10.5;
            }
        }
        prev_vy = vy;
    }
    if h_prime.is_nan() && tracking {
        h_prime = peak_y - 10.5;
    }

    assert!(h_prime.is_finite(), "no rebound detected");
    // Analytic: h'/H = e² = 0.49. Tolerant window for sub-stepping / discrete
    // sampling / gravity-during-contact.
    let ratio = h_prime / 1.0;
    assert!(
        (0.40..0.58).contains(&ratio),
        "rebound ratio out of window: h'/H = {} (expected ≈ 0.49)",
        ratio
    );
}

// ─────────────────── Test 6: soft small-steps smoke test ───────────────────────

/// Exercises the soft path through `step_small_steps` (the cached-contact soft
/// solver): a 5-sphere column with substeps stays finite and bounded.
#[test]
fn soft_small_steps_stack_stable() {
    let mut world = PhysicsWorld::new([0.0, -9.81, 0.0], 1.0 / 60.0);
    world.solver_config.substeps = 4;
    world.solver_config.use_soft_contacts = true;
    world.solver_config.contact_hertz = 30.0;
    world.solver_config.contact_damping_ratio = 1.0;
    world.solver_config.restitution = 0.0;

    spawn_static_ground(&mut world, 8000.0); // radius 10
    let mut handles = Vec::new();
    for &y in &[10.5_f64, 11.5, 12.5, 13.5, 14.5] {
        handles.push(spawn_sphere(&mut world, 1.0, 0.0, y, 0.0));
    }

    for _ in 0..120 {
        world.step_small_steps();
    }

    let max_depth = world
        .last_contacts
        .iter()
        .map(|c| c.depth)
        .fold(0.0_f64, f64::max);
    for h in &handles {
        let b = world.get_body(*h).expect("dynamic body present");
        assert!(b.body_type == BodyType::Dynamic);
        let p = &b.transform.position;
        assert!(
            p.x.is_finite() && p.y.is_finite() && p.z.is_finite(),
            "non-finite position"
        );
        assert!(p.y > 0.0, "body tunnelled below ground centre: y = {}", p.y);
    }
    assert!(
        max_depth < 0.05,
        "soft small-steps penetration too large: {}",
        max_depth
    );
}
