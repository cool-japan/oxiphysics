// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Conversion utilities for WebAssembly / JavaScript interop.
//!
//! These helpers bridge between Rust physics types and the flat arrays /
//! `JsValue` objects that JavaScript expects when consuming wasm-bindgen APIs.

use serde::Serialize;
use wasm_bindgen::prelude::*;

/// Convert any `Display`-able error into a `JsValue` string for use as the
/// error type in `Result<_, JsValue>` wasm-bindgen function signatures.
pub fn err_to_jsvalue(e: impl std::fmt::Display) -> JsValue {
    JsValue::from_str(&e.to_string())
}

/// Flatten a slice of `[f64; 3]` vectors into a contiguous `Vec<f64>` suitable
/// for transferring to a JavaScript `Float64Array`.
pub fn flatten_vec3s(v: &[[f64; 3]]) -> Vec<f64> {
    let mut out = Vec::with_capacity(v.len() * 3);
    for arr in v {
        out.push(arr[0]);
        out.push(arr[1]);
        out.push(arr[2]);
    }
    out
}

/// Reconstruct a `Vec<[f64; 3]>` from a flat `f64` slice.
///
/// # Errors
///
/// Returns an error string if `flat.len()` is not divisible by 3.
pub fn unflatten_vec3s(flat: &[f64]) -> Result<Vec<[f64; 3]>, String> {
    if !flat.len().is_multiple_of(3) {
        return Err(format!(
            "flat slice length {} is not divisible by 3",
            flat.len()
        ));
    }
    Ok(flat.chunks_exact(3).map(|c| [c[0], c[1], c[2]]).collect())
}

/// Flatten a slice of `[f64; 4]` values (e.g. quaternions) into a contiguous
/// `Vec<f64>` for JavaScript consumption.
pub fn flatten_vec4s(v: &[[f64; 4]]) -> Vec<f64> {
    let mut out = Vec::with_capacity(v.len() * 4);
    for arr in v {
        out.push(arr[0]);
        out.push(arr[1]);
        out.push(arr[2]);
        out.push(arr[3]);
    }
    out
}

/// Serialize any `serde::Serialize` value to a `JsValue` using
/// `serde-wasm-bindgen`.
///
/// # Errors
///
/// Returns a `JsValue` error string if serialization fails.
pub fn to_js_value<T: Serialize>(v: &T) -> Result<JsValue, JsValue> {
    serde_wasm_bindgen::to_value(v).map_err(err_to_jsvalue)
}

/// Convert an `Option<[f64; 3]>` to a `Vec<f64>`.
///
/// Returns an empty vector when the option is `None`.
pub fn opt_vec3_to_js(opt: Option<[f64; 3]>) -> Vec<f64> {
    match opt {
        Some(arr) => arr.to_vec(),
        None => Vec::new(),
    }
}

/// Convert an `Option<[f64; 4]>` to a `Vec<f64>`.
///
/// Returns an empty vector when the option is `None`.
pub fn opt_vec4_to_js(opt: Option<[f64; 4]>) -> Vec<f64> {
    match opt {
        Some(arr) => arr.to_vec(),
        None => Vec::new(),
    }
}

// ---------------------------------------------------------------------------
// Monotonic wall-clock
// ---------------------------------------------------------------------------

/// Return a monotonic wall-clock timestamp in **milliseconds**, or `None` when
/// no real clock is available on the current target.
///
/// - On native targets this uses a process-lifetime [`std::time::Instant`]
///   baseline (`Instant` is monotonic but epoch-free, so timestamps are
///   relative to first call — only differences are meaningful).
/// - On `wasm32` with a browser/worker global this uses
///   [`web_sys::Performance::now`] (sub-millisecond resolution, already in ms).
/// - On `wasm32-unknown-unknown` without a `Performance` (e.g. a bare
///   non-browser host), it returns `None` so callers can report an honest
///   sentinel instead of a fabricated number. `std::time::Instant` is *not*
///   used on wasm because it panics there.
///
/// Use two calls and subtract to measure an elapsed interval.
#[cfg(not(target_arch = "wasm32"))]
pub fn now_ms() -> Option<f64> {
    use std::sync::OnceLock;
    use std::time::Instant;
    static BASE: OnceLock<Instant> = OnceLock::new();
    let base = BASE.get_or_init(Instant::now);
    Some(base.elapsed().as_secs_f64() * 1_000.0)
}

/// Monotonic wall-clock timestamp in milliseconds (wasm32 implementation).
///
/// See the native overload for semantics. Resolves `performance.now()` on the
/// global object via `js-sys` reflection so it works from both a `Window` and a
/// `WorkerGlobalScope` (off-thread physics) without requiring additional
/// `web-sys` features. On the main thread the typed
/// [`web_sys::Performance::now`] path is preferred.
#[cfg(target_arch = "wasm32")]
pub fn now_ms() -> Option<f64> {
    // Fast path: typed Window → Performance (Window + Performance features).
    if let Some(perf) = web_sys::window().and_then(|w| w.performance()) {
        return Some(perf.now());
    }
    // Worker / fallback path: globalThis.performance.now() via reflection.
    let global = js_sys::global();
    let performance = js_sys::Reflect::get(&global, &JsValue::from_str("performance")).ok()?;
    if performance.is_undefined() || performance.is_null() {
        return None;
    }
    let now_fn = js_sys::Reflect::get(&performance, &JsValue::from_str("now")).ok()?;
    let now_fn = now_fn.dyn_ref::<js_sys::Function>()?;
    now_fn.call0(&performance).ok()?.as_f64()
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn flatten_round_trip() {
        let pts = vec![[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]];
        let flat = flatten_vec3s(&pts);
        assert_eq!(flat, vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0]);
        let back = unflatten_vec3s(&flat).expect("round-trip");
        assert_eq!(back, pts);
    }

    #[test]
    fn unflatten_rejects_bad_length() {
        assert!(unflatten_vec3s(&[1.0, 2.0]).is_err());
    }

    #[test]
    fn opt_vec3_to_js_none_is_empty() {
        assert!(opt_vec3_to_js(None).is_empty());
    }

    #[test]
    fn opt_vec3_to_js_some() {
        assert_eq!(opt_vec3_to_js(Some([1.0, 2.0, 3.0])), vec![1.0, 2.0, 3.0]);
    }

    #[cfg(not(target_arch = "wasm32"))]
    #[test]
    fn now_ms_is_real_and_monotonic() {
        // On the host this must be a real Instant-based measurement, not a
        // sentinel — and it must be non-decreasing across calls.
        let t0 = now_ms().expect("native clock should be available");
        assert!(t0.is_finite() && t0 >= 0.0);
        let mut spins = 0u64;
        let mut t1 = now_ms().expect("native clock should be available");
        // Busy-spin until the monotonic clock visibly advances.
        while t1 <= t0 && spins < 100_000_000 {
            spins += 1;
            t1 = now_ms().expect("native clock should be available");
        }
        assert!(t1 >= t0, "clock went backwards: {t0} -> {t1}");
    }
}
