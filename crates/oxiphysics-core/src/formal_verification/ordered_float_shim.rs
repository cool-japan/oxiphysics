// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

/// A wrapper that makes `f64` implement `Eq + Ord + Hash` by treating NaN
/// as larger than all non-NaN values (i.e., `NaN == NaN`, `NaN > finite`).
#[derive(Debug, Clone, Copy)]
pub struct NotNan(f64);

impl NotNan {
    /// Wraps `v`.
    pub fn new(v: f64) -> Self {
        Self(v)
    }
    /// Returns the inner value.
    pub fn into_inner(self) -> f64 {
        self.0
    }
}

impl PartialEq for NotNan {
    fn eq(&self, other: &Self) -> bool {
        self.0.to_bits() == other.0.to_bits()
    }
}
impl Eq for NotNan {}
impl PartialOrd for NotNan {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}
impl Ord for NotNan {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.0
            .partial_cmp(&other.0)
            .unwrap_or_else(|| self.0.is_nan().cmp(&other.0.is_nan()))
    }
}
impl std::hash::Hash for NotNan {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.0.to_bits().hash(state);
    }
}
