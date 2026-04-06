// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! `Float64View` — zero-copy helper for flat f64 data at the WASM boundary.

#![allow(missing_docs)]

/// A lightweight wrapper around a `Vec`f64` that represents a view into
/// a flat floating-point buffer suitable for JavaScript `Float64Array`.
///
/// In a real wasm-bindgen build this would be annotated with `#[wasm_bindgen]`
/// and return `Float64Array` from a typed-array view. Here it carries the data
/// as an owned `Vec`f64` so that host-side (non-WASM) tests can use it.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct Float64View {
    data: Vec<f64>,
}

impl Float64View {
    /// Create a new view from a `Vec`f64`.
    pub fn from_vec(data: Vec<f64>) -> Self {
        Self { data }
    }

    /// Create an empty view.
    pub fn empty() -> Self {
        Self { data: Vec::new() }
    }

    /// Return a slice of the underlying data.
    pub fn as_slice(&self) -> &[f64] {
        &self.data
    }

    /// Number of elements.
    pub fn len(&self) -> usize {
        self.data.len()
    }

    /// Whether the view contains no elements.
    pub fn is_empty(&self) -> bool {
        self.data.is_empty()
    }

    /// Consume this view and return the inner `Vec`f64`.
    pub fn into_vec(self) -> Vec<f64> {
        self.data
    }

    /// Get element at index, or `None` if out of bounds.
    pub fn get(&self, index: usize) -> Option<f64> {
        self.data.get(index).copied()
    }
}

impl From<Vec<f64>> for Float64View {
    fn from(v: Vec<f64>) -> Self {
        Self::from_vec(v)
    }
}
