// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Small helpers shared by the binary readers of this crate.

/// Largest element count it is safe to pre-allocate for a collection whose
/// size was read from untrusted input.
///
/// `declared` is the count stored in the file header, `remaining_bytes` the
/// number of input bytes that are still unread and `min_record_bytes` the
/// smallest number of bytes one element can occupy in the input.  A valid file
/// cannot hold more than `remaining_bytes / min_record_bytes` elements, so the
/// returned value never exceeds that bound.  A corrupt header that declares
/// billions of elements therefore cannot make `Vec::with_capacity` abort the
/// process; the reader simply runs out of input and reports a truncation error
/// while the vector grows on demand.
///
/// A `min_record_bytes` of zero is treated as one byte per record.
pub(crate) fn bounded_capacity(
    declared: u64,
    remaining_bytes: usize,
    min_record_bytes: usize,
) -> usize {
    usize::try_from(declared)
        .unwrap_or(usize::MAX)
        .min(remaining_bytes / min_record_bytes.max(1))
}

#[cfg(test)]
mod tests {
    use super::bounded_capacity;

    #[test]
    fn honest_count_is_kept() {
        assert_eq!(bounded_capacity(10, 1000, 24), 10);
        assert_eq!(bounded_capacity(0, 1000, 24), 0);
    }

    #[test]
    fn count_is_limited_by_the_remaining_bytes() {
        assert_eq!(bounded_capacity(1_000_000, 240, 24), 10);
        assert_eq!(bounded_capacity(5, 23, 24), 0);
        assert_eq!(bounded_capacity(5, 24, 24), 1);
    }

    #[test]
    fn hostile_counts_never_exceed_the_input_size() {
        assert_eq!(bounded_capacity(u64::MAX, 4096, 8), 512);
        assert_eq!(bounded_capacity(u64::MAX, 0, 8), 0);
        assert_eq!(bounded_capacity(u64::MAX, usize::MAX, 1), usize::MAX);
    }

    #[test]
    fn zero_record_size_counts_as_one_byte() {
        assert_eq!(bounded_capacity(100, 7, 0), 7);
    }
}
