// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! SDF combinators: Boolean CSG operations and smooth blending.

use super::helpers::clamp;

// ─────────────────────────────────────────────────────────────────────────────
// Boolean CSG
// ─────────────────────────────────────────────────────────────────────────────

/// Boolean union of two SDFs: min(a, b).
#[inline]
pub fn sdf_union(a: f64, b: f64) -> f64 {
    a.min(b)
}

/// Boolean intersection of two SDFs: max(a, b).
#[inline]
pub fn sdf_intersection(a: f64, b: f64) -> f64 {
    a.max(b)
}

/// Boolean difference of two SDFs (a minus b): max(a, −b).
#[inline]
pub fn sdf_difference(a: f64, b: f64) -> f64 {
    a.max(-b)
}

// ─────────────────────────────────────────────────────────────────────────────
// Smooth blending combinators
// ─────────────────────────────────────────────────────────────────────────────

/// Smooth union of two SDFs with blend radius `k` (polynomial C¹).
///
/// Smoothly blends the two surfaces within distance `k` of their intersection.
pub fn sdf_smooth_union(a: f64, b: f64, k: f64) -> f64 {
    let h = clamp(0.5 + 0.5 * (b - a) / k, 0.0, 1.0);
    a * h + b * (1.0 - h) - k * h * (1.0 - h)
}

/// Smooth intersection of two SDFs with blend radius `k`.
pub fn sdf_smooth_intersection(a: f64, b: f64, k: f64) -> f64 {
    let h = clamp(0.5 - 0.5 * (b - a) / k, 0.0, 1.0);
    a * h + b * (1.0 - h) + k * h * (1.0 - h)
}

/// Smooth difference (a minus b) with blend radius `k`.
pub fn sdf_smooth_difference(a: f64, b: f64, k: f64) -> f64 {
    let h = clamp(0.5 - 0.5 * (b + a) / k, 0.0, 1.0);
    a * h + (-b) * (1.0 - h) + k * h * (1.0 - h)
}

/// Exponential smooth minimum.
///
/// Returns −(1/k) ln(e^(−k·a) + e^(−k·b)).
pub fn sdf_exp_smooth_union(a: f64, b: f64, k: f64) -> f64 {
    let ea = (-k * a).exp();
    let eb = (-k * b).exp();
    -(ea + eb).ln() / k
}
