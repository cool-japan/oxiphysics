// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0
//! Regression tests for `BuoyancyForce::cylinder_partial_buoyancy`.
//!
//! The wetted cross-section of a horizontal cylinder is the circular segment
//! of height `h`: `A(h) = r² acos((r − h)/r) − (r − h) sqrt(2 r h − h²)`.

use oxiphysics_rigid::BuoyancyForce;
use std::f64::consts::PI;

const RHO: f64 = 1000.0;
const G: [f64; 3] = [0.0, -9.81, 0.0];

/// Magnitude of the (upward) buoyancy force for a cylinder of radius `r` and
/// length `len` submerged to depth `h`.
fn buoyancy(r: f64, h: f64, len: f64) -> f64 {
    let f = BuoyancyForce::cylinder_partial_buoyancy(r, h, len, RHO, G);
    assert!(
        f[0].abs() < 1e-12 && f[2].abs() < 1e-12,
        "force is vertical"
    );
    f[1]
}

#[test]
fn dry_cylinder_has_no_buoyancy() {
    assert!(buoyancy(0.5, 0.0, 3.0).abs() < 1e-12);
}

#[test]
fn half_submerged_cylinder_displaces_half_the_volume() {
    let (r, len) = (0.5, 3.0);
    let expected = RHO * 9.81 * PI * r * r * len / 2.0;
    let got = buoyancy(r, r, len);
    assert!(
        (got - expected).abs() < 1e-9 * expected,
        "half-submerged: got {got}, expected {expected}"
    );
}

#[test]
fn fully_submerged_cylinder_displaces_the_full_volume() {
    let (r, len) = (0.5, 3.0);
    let expected = RHO * 9.81 * PI * r * r * len;
    let got = buoyancy(r, 2.0 * r, len);
    assert!(
        (got - expected).abs() < 1e-9 * expected,
        "fully submerged: got {got}, expected {expected}"
    );
}

#[test]
fn depth_beyond_the_diameter_is_clamped_to_full_submersion() {
    let (r, len) = (0.25, 2.0);
    let full = buoyancy(r, 2.0 * r, len);
    let deeper = buoyancy(r, 10.0 * r, len);
    assert!((full - deeper).abs() < 1e-9 * full);
    let above = buoyancy(r, -1.0, len);
    assert!(above.abs() < 1e-12, "negative depth means dry: {above}");
}

#[test]
fn buoyancy_increases_monotonically_with_depth() {
    let (r, len) = (0.4, 1.5);
    let mut previous = buoyancy(r, 0.0, len);
    for step in 1..=64 {
        let h = 2.0 * r * f64::from(step) / 64.0;
        let current = buoyancy(r, h, len);
        assert!(
            current > previous,
            "buoyancy must grow with depth: h={h}, {previous} -> {current}"
        );
        previous = current;
    }
}

#[test]
fn quarter_depth_matches_the_closed_form_segment_area() {
    // h = r/2: A = r² acos(1/2) − (r/2) sqrt(r² − r²/4) = r² (π/3 − sqrt(3)/4).
    let (r, len) = (2.0, 0.5);
    let area = r * r * (PI / 3.0 - 3.0_f64.sqrt() / 4.0);
    let expected = RHO * 9.81 * area * len;
    let got = buoyancy(r, 0.5 * r, len);
    assert!(
        (got - expected).abs() < 1e-9 * expected,
        "h = r/2: got {got}, expected {expected}"
    );
}

#[test]
fn buoyancy_scales_linearly_with_length_and_density() {
    let base = buoyancy(0.3, 0.2, 1.0);
    let doubled = buoyancy(0.3, 0.2, 2.0);
    assert!((doubled - 2.0 * base).abs() < 1e-9 * base);
    let heavy = BuoyancyForce::cylinder_partial_buoyancy(0.3, 0.2, 1.0, 2.0 * RHO, G);
    assert!((heavy[1] - 2.0 * base).abs() < 1e-9 * base);
}

#[test]
fn degenerate_radius_gives_zero_force_instead_of_panicking() {
    for radius in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        let f = BuoyancyForce::cylinder_partial_buoyancy(radius, 0.3, 1.0, RHO, G);
        assert!(
            f.iter().all(|c| c.abs() < 1e-12),
            "radius {radius} must yield zero force, got {f:?}"
        );
    }
}
