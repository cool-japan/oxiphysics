// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Generational handle encoding for bodies and colliders.
//!
//! A handle is a `u32` that packs a slot index into its low
//! [`HANDLE_INDEX_BITS`] bits and a generation counter into the remaining
//! high bits:
//!
//! ```text
//!  31            20 19                    0
//! +----------------+-----------------------+
//! |   generation   |      slot index       |
//! +----------------+-----------------------+
//! ```
//!
//! Every time a slot is recycled its generation is incremented, so a handle
//! that was returned for a previous occupant of the slot no longer resolves
//! (the classic ABA problem of plain index handles).  A freshly allocated slot
//! starts at generation `0`, hence on an engine that never recycles slots the
//! handle of the `n`-th body is simply `n`.
//!
//! The generation is 12 bits wide and wraps around, so a stale handle only
//! aliases a live one after a slot has been recycled 4096 times.

/// Number of low handle bits that hold the slot index.
pub(super) const HANDLE_INDEX_BITS: u32 = 20;

/// Mask selecting the slot-index bits of a handle.
pub(super) const HANDLE_INDEX_MASK: u32 = (1 << HANDLE_INDEX_BITS) - 1;

/// Mask selecting the generation after it has been shifted down.
pub(super) const HANDLE_GENERATION_MASK: u32 = u32::MAX >> HANDLE_INDEX_BITS;

/// Maximum number of body (or collider) slots an engine may allocate.
///
/// The all-ones index `0xF_FFFF` is never allocated, which guarantees that no
/// valid handle equals `u32::MAX` - the value used throughout the crate as the
/// "no body" / "world" sentinel.
pub const MAX_SLOTS: usize = HANDLE_INDEX_MASK as usize;

/// Pack a slot `index` and `generation` into a handle.
///
/// The generation is reduced modulo 2^12; the index must be below
/// [`MAX_SLOTS`] (callers enforce this when allocating slots).
#[inline]
pub(super) fn encode_handle(index: usize, generation: u32) -> u32 {
    ((generation & HANDLE_GENERATION_MASK) << HANDLE_INDEX_BITS)
        | ((index as u32) & HANDLE_INDEX_MASK)
}

/// Split a handle into its `(slot index, generation)` parts.
#[inline]
pub(super) fn decode_handle(handle: u32) -> (usize, u32) {
    (
        (handle & HANDLE_INDEX_MASK) as usize,
        handle >> HANDLE_INDEX_BITS,
    )
}

/// The generation that follows `generation` (wrapping at 12 bits).
#[inline]
pub(super) fn next_generation(generation: u32) -> u32 {
    generation.wrapping_add(1) & HANDLE_GENERATION_MASK
}
