// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Deterministic, platform-independent, non-cryptographic hashing.
//!
//! Hash values that are exchanged between peers (rollback desync detection,
//! replay verification, WASM <-> native parity checks) must be bit-identical on
//! every target.  This module is the single place in the workspace that
//! defines the 64-bit FNV-1a function used for that purpose, so that the native
//! engine, the Python bindings and the WASM bindings can never drift apart.
//!
//! FNV-1a is **not** a cryptographic hash: it detects accidental divergence of
//! simulation state, not tampering.

/// Offset basis of the 64-bit FNV-1a hash (`0xcbf29ce484222325`).
pub const FNV1A_64_OFFSET_BASIS: u64 = 0xcbf2_9ce4_8422_2325;

/// Prime of the 64-bit FNV-1a hash (`2^40 + 2^8 + 0xb3 = 0x100000001b3`).
pub const FNV1A_64_PRIME: u64 = 0x0000_0100_0000_01b3;

/// Compute the 64-bit FNV-1a hash of `bytes`.
///
/// The hash of the empty slice is [`FNV1A_64_OFFSET_BASIS`].  For every byte
/// the running hash is XOR-ed with the byte and then multiplied (wrapping) by
/// [`FNV1A_64_PRIME`].
///
/// # Examples
///
/// ```
/// use oxiphysics_core::fnv1a_64;
///
/// assert_eq!(fnv1a_64(b""), 0xcbf2_9ce4_8422_2325);
/// assert_eq!(fnv1a_64(b"a"), 0xaf63_dc4c_8601_ec8c);
/// assert_eq!(fnv1a_64(b"foobar"), 0x8594_4171_f739_67e8);
/// ```
#[must_use]
pub fn fnv1a_64(bytes: &[u8]) -> u64 {
    let mut hash = FNV1A_64_OFFSET_BASIS;
    for &byte in bytes {
        hash ^= u64::from(byte);
        hash = hash.wrapping_mul(FNV1A_64_PRIME);
    }
    hash
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn empty_input_hashes_to_the_offset_basis() {
        assert_eq!(fnv1a_64(b""), 0xcbf2_9ce4_8422_2325);
        assert_eq!(fnv1a_64(&[]), FNV1A_64_OFFSET_BASIS);
    }

    #[test]
    fn constants_match_the_fnv_specification() {
        assert_eq!(FNV1A_64_OFFSET_BASIS, 14_695_981_039_346_656_037);
        assert_eq!(FNV1A_64_PRIME, 1_099_511_628_211);
        assert_eq!(FNV1A_64_PRIME, (1_u64 << 40) + 0x1b3);
    }

    #[test]
    fn matches_published_fnv1a_64_vectors() {
        // Vectors from the reference FNV test suite.
        assert_eq!(fnv1a_64(b"a"), 0xaf63_dc4c_8601_ec8c);
        assert_eq!(fnv1a_64(b"b"), 0xaf63_df4c_8601_f1a5);
        assert_eq!(fnv1a_64(b"foo"), 0xdcb2_7518_fed9_d577);
        assert_eq!(fnv1a_64(b"foobar"), 0x8594_4171_f739_67e8);
    }

    #[test]
    fn matches_independent_reference_for_binary_snapshots() {
        // Reference values computed with an independent big-integer
        // implementation (Python) of the FNV-1a recurrence.
        assert_eq!(fnv1a_64(&[0xab, 0xcd]), 0x0997_9b07_b61b_0365);
        assert_eq!(fnv1a_64(&5_i64.to_le_bytes()), 0x0de2_1504_f16d_c720);
        assert_eq!(fnv1a_64(&42_i64.to_le_bytes()), 0xff3a_dd6b_3789_daef);
        assert_eq!(fnv1a_64(b"hello rollback"), 0x396b_7f0e_63df_8d71);
    }

    #[test]
    fn differs_from_the_previous_incorrect_prime() {
        // The pre-fix rollback hash multiplied by 0x100000003d811, which is not
        // the FNV prime; it hashed "a" to 0x396ec70ffdd1b784.
        assert_ne!(fnv1a_64(b"a"), 0x396e_c70f_fdd1_b784);
    }

    #[test]
    fn is_sensitive_to_byte_order_and_single_bit_flips() {
        assert_ne!(fnv1a_64(b"ab"), fnv1a_64(b"ba"));
        assert_ne!(fnv1a_64(&[0x00]), fnv1a_64(&[0x01]));
        assert_ne!(fnv1a_64(&[0x00, 0x00]), fnv1a_64(&[0x00]));
    }
}
