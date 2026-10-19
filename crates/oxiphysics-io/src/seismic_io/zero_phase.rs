// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Zero-phase (forward-backward) filtering with edge handling.

use oxiphysics_core::signal::ButterworthBandpass;

/// Extend `data` by `pad` samples at both ends with its odd reflection about
/// the end points: `x[-k] = 2 x[0] - x[k]` and `x[n-1+k] = 2 x[n-1] - x[n-1-k]`.
///
/// The extension continues a signal's value and slope across the boundary, so
/// the filter transient is spent inside the padding. `pad` is clamped to
/// `data.len() - 1`.
fn odd_extension(data: &[f64], pad: usize) -> Vec<f64> {
    let n = data.len();
    let (Some(&first), Some(&last)) = (data.first(), data.last()) else {
        return Vec::new();
    };
    let pad = pad.min(n - 1);
    let mut extended = Vec::with_capacity(n + 2 * pad);
    for k in (1..=pad).rev() {
        extended.push(2.0 * first - data[k]);
    }
    extended.extend_from_slice(data);
    for k in 1..=pad {
        extended.push(2.0 * last - data[n - 1 - k]);
    }
    extended
}

/// One causal pass of `filter` over `signal`, started from rest.
///
/// A band-pass has exactly zero gain at DC, so a constant offset in the input
/// only changes the start-up transient, never the steady-state output.
/// Subtracting the first sample therefore starts the filter in the steady
/// state of a constant input equal to that sample, which keeps a non-zero
/// first value from ringing through the beginning of the record.
fn causal_pass(filter: &mut ButterworthBandpass, signal: &[f64]) -> Vec<f64> {
    filter.reset();
    let offset = signal.first().copied().unwrap_or(0.0);
    let shifted: Vec<f64> = signal.iter().map(|&value| value - offset).collect();
    filter.apply(&shifted)
}

/// Zero-phase filtering: run `filter` forward over the odd-extended signal,
/// then backward over the result, and drop the padding.
///
/// The output has the length of `data`; its magnitude response is the square
/// of the filter's, and its phase response is zero.
pub(super) fn filtfilt(filter: &mut ButterworthBandpass, data: &[f64], pad: usize) -> Vec<f64> {
    let n = data.len();
    if n == 0 {
        return Vec::new();
    }
    let pad = pad.min(n - 1);
    let extended = odd_extension(data, pad);
    let mut reversed = causal_pass(filter, &extended);
    reversed.reverse();
    let mut result = causal_pass(filter, &reversed);
    result.reverse();
    result[pad..pad + n].to_vec()
}
