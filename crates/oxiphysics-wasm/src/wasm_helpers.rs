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
}
