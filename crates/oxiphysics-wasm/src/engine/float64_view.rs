// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! `Float64View` — zero-copy helper for flat f64 data at the WASM boundary.

use wasm_bindgen::prelude::*;

/// A lightweight wrapper around a `Vec<f64>` that represents a view into
/// a flat floating-point buffer suitable for JavaScript `Float64Array`.
///
/// Exposed to JavaScript with typed-array getters and index access.
#[wasm_bindgen]
#[derive(Debug, Clone)]
pub struct Float64View {
    data: Vec<f64>,
}

#[wasm_bindgen]
impl Float64View {
    /// Create a new view from a `Vec<f64>`.
    pub fn from_vec(data: Vec<f64>) -> Self {
        Self { data }
    }

    /// Create an empty view.
    pub fn empty() -> Self {
        Self { data: Vec::new() }
    }

    /// Number of elements (JS-compatible `u32`).
    pub fn len(&self) -> u32 {
        self.data.len() as u32
    }

    /// Whether the view contains no elements.
    pub fn is_empty(&self) -> bool {
        self.data.is_empty()
    }

    /// Consume this view and return the inner data as `Vec<f64>`.
    pub fn into_vec(self) -> Vec<f64> {
        self.data
    }

    /// Get element at index (returns `NaN` if out of bounds).
    pub fn get(&self, index: u32) -> f64 {
        self.data.get(index as usize).copied().unwrap_or(f64::NAN)
    }

    /// Return the data as a JavaScript `Float64Array` (copies the buffer).
    pub fn to_typed_array(&self) -> js_sys::Float64Array {
        js_sys::Float64Array::from(self.data.as_slice())
    }
}

impl Float64View {
    /// Return a slice of the underlying data (Rust-only).
    pub fn as_slice(&self) -> &[f64] {
        &self.data
    }

    /// Get element at index, or `None` if out of bounds (Rust-only).
    pub fn get_option(&self, index: usize) -> Option<f64> {
        self.data.get(index).copied()
    }
}

impl From<Vec<f64>> for Float64View {
    fn from(v: Vec<f64>) -> Self {
        Self::from_vec(v)
    }
}
