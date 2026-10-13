// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Sobol low-discrepancy sequence direction numbers.
//!
//! The direction numbers are the first 128 dimensions of the Joe-Kuo table
//! `new-joe-kuo-6.21201` (S. Joe and F. Y. Kuo, "Constructing Sobol sequences
//! with better two-dimensional projections", SIAM J. Sci. Comput. 30 (2008),
//! 2635-2654).  Every dimension is generated from a primitive polynomial over
//! GF(2) and a set of initial odd integers; the remaining direction numbers
//! follow from the Bratley-Fox recurrence.
//!
//! # Attribution
//!
//! The tabulated primitive polynomials and initial direction integers in
//! `JOE_KUO_ROWS` are taken from the `new-joe-kuo-6.21201` data file, which is
//! distributed under the BSD 3-clause licence:
//!
//! ```text
//! Copyright (c) 2008, Frances Y. Kuo and Stephen Joe
//! All rights reserved.
//!
//! Redistribution and use in source and binary forms, with or without
//! modification, are permitted provided that the following conditions are met:
//!   * Redistributions of source code must retain the above copyright notice,
//!     this list of conditions and the following disclaimer.
//!   * Redistributions in binary form must reproduce the above copyright
//!     notice, this list of conditions and the following disclaimer in the
//!     documentation and/or other materials provided with the distribution.
//!   * Neither the names of the copyright holders nor the names of the
//!     contributors may be used to endorse or promote products derived from
//!     this software without specific prior written permission.
//!
//! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
//! AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
//! IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
//! ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR CONTRIBUTORS BE
//! LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
//! CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
//! SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
//! INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
//! CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
//! ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
//! POSSIBILITY OF SUCH DAMAGE.
//! ```

use crate::error::{PhysicsError, PhysicsResult};

/// Number of bits of the integer state of a Sobol coordinate (resolution `2^-32`).
pub const SOBOL_BITS: usize = 32;

/// One row of the Joe-Kuo table (dimensions 2, 3, ... in the published numbering):
/// `(s, a, m)` with
///
/// * `s` the degree of the primitive polynomial `x^s + a_1 x^(s-1) + ... + a_(s-1) x + 1`,
/// * `a` the interior coefficients `a_1 ... a_(s-1)` packed as an integer (`a_1` in the highest bit),
/// * `m` the initial direction integers `m_1 ... m_s` (odd, `m_k < 2^k`).
type JoeKuoRow = (u32, u32, &'static [u32]);

/// Joe-Kuo rows for coordinates `1..SOBOL_MAX_DIM` (coordinate `0` is the van der Corput sequence).
const JOE_KUO_ROWS: [JoeKuoRow; 127] = [
    (1, 0, &[1]),
    (2, 1, &[1, 3]),
    (3, 1, &[1, 3, 1]),
    (3, 2, &[1, 1, 1]),
    (4, 1, &[1, 1, 3, 3]),
    (4, 4, &[1, 3, 5, 13]),
    (5, 2, &[1, 1, 5, 5, 17]),
    (5, 4, &[1, 1, 5, 5, 5]),
    (5, 7, &[1, 1, 7, 11, 19]),
    (5, 11, &[1, 1, 5, 1, 1]),
    (5, 13, &[1, 1, 1, 3, 11]),
    (5, 14, &[1, 3, 5, 5, 31]),
    (6, 1, &[1, 3, 3, 9, 7, 49]),
    (6, 13, &[1, 1, 1, 15, 21, 21]),
    (6, 16, &[1, 3, 1, 13, 27, 49]),
    (6, 19, &[1, 1, 1, 15, 7, 5]),
    (6, 22, &[1, 3, 1, 15, 13, 25]),
    (6, 25, &[1, 1, 5, 5, 19, 61]),
    (7, 1, &[1, 3, 7, 11, 23, 15, 103]),
    (7, 4, &[1, 3, 7, 13, 13, 15, 69]),
    (7, 7, &[1, 1, 3, 13, 7, 35, 63]),
    (7, 8, &[1, 3, 5, 9, 1, 25, 53]),
    (7, 14, &[1, 3, 1, 13, 9, 35, 107]),
    (7, 19, &[1, 3, 1, 5, 27, 61, 31]),
    (7, 21, &[1, 1, 5, 11, 19, 41, 61]),
    (7, 28, &[1, 3, 5, 3, 3, 13, 69]),
    (7, 31, &[1, 1, 7, 13, 1, 19, 1]),
    (7, 32, &[1, 3, 7, 5, 13, 19, 59]),
    (7, 37, &[1, 1, 3, 9, 25, 29, 41]),
    (7, 41, &[1, 3, 5, 13, 23, 1, 55]),
    (7, 42, &[1, 3, 7, 3, 13, 59, 17]),
    (7, 50, &[1, 3, 1, 3, 5, 53, 69]),
    (7, 55, &[1, 1, 5, 5, 23, 33, 13]),
    (7, 56, &[1, 1, 7, 7, 1, 61, 123]),
    (7, 59, &[1, 1, 7, 9, 13, 61, 49]),
    (7, 62, &[1, 3, 3, 5, 3, 55, 33]),
    (8, 14, &[1, 3, 1, 15, 31, 13, 49, 245]),
    (8, 21, &[1, 3, 5, 15, 31, 59, 63, 97]),
    (8, 22, &[1, 3, 1, 11, 11, 11, 77, 249]),
    (8, 38, &[1, 3, 1, 11, 27, 43, 71, 9]),
    (8, 47, &[1, 1, 7, 15, 21, 11, 81, 45]),
    (8, 49, &[1, 3, 7, 3, 25, 31, 65, 79]),
    (8, 50, &[1, 3, 1, 1, 19, 11, 3, 205]),
    (8, 52, &[1, 1, 5, 9, 19, 21, 29, 157]),
    (8, 56, &[1, 3, 7, 11, 1, 33, 89, 185]),
    (8, 67, &[1, 3, 3, 3, 15, 9, 79, 71]),
    (8, 70, &[1, 3, 7, 11, 15, 39, 119, 27]),
    (8, 84, &[1, 1, 3, 1, 11, 31, 97, 225]),
    (8, 97, &[1, 1, 1, 3, 23, 43, 57, 177]),
    (8, 103, &[1, 3, 7, 7, 17, 17, 37, 71]),
    (8, 115, &[1, 3, 1, 5, 27, 63, 123, 213]),
    (8, 122, &[1, 1, 3, 5, 11, 43, 53, 133]),
    (9, 8, &[1, 3, 5, 5, 29, 17, 47, 173, 479]),
    (9, 13, &[1, 3, 3, 11, 3, 1, 109, 9, 69]),
    (9, 16, &[1, 1, 1, 5, 17, 39, 23, 5, 343]),
    (9, 22, &[1, 3, 1, 5, 25, 15, 31, 103, 499]),
    (9, 25, &[1, 1, 1, 11, 11, 17, 63, 105, 183]),
    (9, 44, &[1, 1, 5, 11, 9, 29, 97, 231, 363]),
    (9, 47, &[1, 1, 5, 15, 19, 45, 41, 7, 383]),
    (9, 52, &[1, 3, 7, 7, 31, 19, 83, 137, 221]),
    (9, 55, &[1, 1, 1, 3, 23, 15, 111, 223, 83]),
    (9, 59, &[1, 1, 5, 13, 31, 15, 55, 25, 161]),
    (9, 62, &[1, 1, 3, 13, 25, 47, 39, 87, 257]),
    (9, 67, &[1, 1, 1, 11, 21, 53, 125, 249, 293]),
    (9, 74, &[1, 1, 7, 11, 11, 7, 57, 79, 323]),
    (9, 81, &[1, 1, 5, 5, 17, 13, 81, 3, 131]),
    (9, 82, &[1, 1, 7, 13, 23, 7, 65, 251, 475]),
    (9, 87, &[1, 3, 5, 1, 9, 43, 3, 149, 11]),
    (9, 91, &[1, 1, 3, 13, 31, 13, 13, 255, 487]),
    (9, 94, &[1, 3, 3, 1, 5, 63, 89, 91, 127]),
    (9, 103, &[1, 1, 3, 3, 1, 19, 123, 127, 237]),
    (9, 104, &[1, 1, 5, 7, 23, 31, 37, 243, 289]),
    (9, 109, &[1, 1, 5, 11, 17, 53, 117, 183, 491]),
    (9, 122, &[1, 1, 1, 5, 1, 13, 13, 209, 345]),
    (9, 124, &[1, 1, 3, 15, 1, 57, 115, 7, 33]),
    (9, 137, &[1, 3, 1, 11, 7, 43, 81, 207, 175]),
    (9, 138, &[1, 3, 1, 1, 15, 27, 63, 255, 49]),
    (9, 143, &[1, 3, 5, 3, 27, 61, 105, 171, 305]),
    (9, 145, &[1, 1, 5, 3, 1, 3, 57, 249, 149]),
    (9, 152, &[1, 1, 3, 5, 5, 57, 15, 13, 159]),
    (9, 157, &[1, 1, 1, 11, 7, 11, 105, 141, 225]),
    (9, 167, &[1, 3, 3, 5, 27, 59, 121, 101, 271]),
    (9, 173, &[1, 3, 5, 9, 11, 49, 51, 59, 115]),
    (9, 176, &[1, 1, 7, 1, 23, 45, 125, 71, 419]),
    (9, 181, &[1, 1, 3, 5, 23, 5, 105, 109, 75]),
    (9, 182, &[1, 1, 7, 15, 7, 11, 67, 121, 453]),
    (9, 185, &[1, 3, 7, 3, 9, 13, 31, 27, 449]),
    (9, 191, &[1, 3, 1, 15, 19, 39, 39, 89, 15]),
    (9, 194, &[1, 1, 1, 1, 1, 33, 73, 145, 379]),
    (9, 199, &[1, 3, 1, 15, 15, 43, 29, 13, 483]),
    (9, 218, &[1, 1, 7, 3, 19, 27, 85, 131, 431]),
    (9, 220, &[1, 3, 3, 3, 5, 35, 23, 195, 349]),
    (9, 227, &[1, 3, 3, 7, 9, 27, 39, 59, 297]),
    (9, 229, &[1, 1, 3, 9, 11, 17, 13, 241, 157]),
    (9, 230, &[1, 3, 7, 15, 25, 57, 33, 189, 213]),
    (9, 234, &[1, 1, 7, 1, 9, 55, 73, 83, 217]),
    (9, 236, &[1, 3, 3, 13, 19, 27, 23, 113, 249]),
    (9, 241, &[1, 3, 5, 3, 23, 43, 3, 253, 479]),
    (9, 244, &[1, 1, 5, 5, 11, 5, 45, 117, 217]),
    (9, 253, &[1, 3, 3, 7, 29, 37, 33, 123, 147]),
    (10, 4, &[1, 3, 1, 15, 5, 5, 37, 227, 223, 459]),
    (10, 13, &[1, 1, 7, 5, 5, 39, 63, 255, 135, 487]),
    (10, 19, &[1, 3, 1, 7, 9, 7, 87, 249, 217, 599]),
    (10, 22, &[1, 1, 3, 13, 9, 47, 7, 225, 363, 247]),
    (10, 50, &[1, 3, 7, 13, 19, 13, 9, 67, 9, 737]),
    (10, 55, &[1, 3, 5, 5, 19, 59, 7, 41, 319, 677]),
    (10, 64, &[1, 1, 5, 3, 31, 63, 15, 43, 207, 789]),
    (10, 69, &[1, 1, 7, 9, 13, 39, 3, 47, 497, 169]),
    (10, 98, &[1, 3, 1, 7, 21, 17, 97, 19, 415, 905]),
    (10, 107, &[1, 3, 7, 1, 3, 31, 71, 111, 165, 127]),
    (10, 115, &[1, 1, 5, 11, 1, 61, 83, 119, 203, 847]),
    (10, 121, &[1, 3, 3, 13, 9, 61, 19, 97, 47, 35]),
    (10, 127, &[1, 1, 7, 7, 15, 29, 63, 95, 417, 469]),
    (10, 134, &[1, 3, 1, 9, 25, 9, 71, 57, 213, 385]),
    (10, 140, &[1, 3, 5, 13, 31, 47, 101, 57, 39, 341]),
    (10, 145, &[1, 1, 3, 3, 31, 57, 125, 173, 365, 551]),
    (10, 152, &[1, 3, 7, 1, 13, 57, 67, 157, 451, 707]),
    (10, 158, &[1, 1, 1, 7, 21, 13, 105, 89, 429, 965]),
    (10, 161, &[1, 1, 5, 9, 17, 51, 45, 119, 157, 141]),
    (10, 171, &[1, 3, 7, 7, 13, 45, 91, 9, 129, 741]),
    (10, 181, &[1, 3, 7, 1, 23, 57, 67, 141, 151, 571]),
    (10, 194, &[1, 1, 3, 11, 17, 47, 93, 107, 375, 157]),
    (10, 199, &[1, 3, 3, 5, 11, 21, 43, 51, 169, 915]),
    (10, 203, &[1, 1, 5, 3, 15, 55, 101, 67, 455, 625]),
    (10, 208, &[1, 3, 5, 9, 1, 23, 29, 47, 345, 595]),
    (10, 227, &[1, 3, 7, 7, 5, 49, 29, 155, 323, 589]),
    (10, 242, &[1, 3, 3, 7, 5, 41, 127, 61, 261, 717]),
];

/// Largest supported number of dimensions.
pub const SOBOL_MAX_DIM: usize = JOE_KUO_ROWS.len() + 1;

/// Direction numbers `v_1 ... v_32` of the Sobol coordinate `dim` (0-based).
///
/// `v_k` is stored as a 32-bit integer scaled by `2^32`, i.e. `v_k = m_k * 2^(32-k)`.
/// Coordinate `0` is the van der Corput sequence in base 2 (`m_k = 1` for all `k`).
/// For the other coordinates `v_1 ... v_s` come from the tabulated initial
/// integers and the rest follow from
/// `v_k = a_1 v_(k-1) ^ a_2 v_(k-2) ^ ... ^ a_(s-1) v_(k-s+1) ^ v_(k-s) ^ (v_(k-s) >> s)`.
///
/// # Errors
///
/// Returns [`PhysicsError::OutOfBounds`] when `dim >= SOBOL_MAX_DIM`.
pub fn sobol_directions(dim: usize) -> PhysicsResult<[u32; SOBOL_BITS]> {
    let mut v = [0u32; SOBOL_BITS];
    if dim == 0 {
        for (k, slot) in v.iter_mut().enumerate() {
            *slot = 1u32 << (SOBOL_BITS - 1 - k);
        }
        return Ok(v);
    }
    let row = JOE_KUO_ROWS.get(dim - 1).ok_or(PhysicsError::OutOfBounds {
        index: dim,
        max: SOBOL_MAX_DIM,
    })?;
    let &(degree, coeff, initial) = row;
    let s = degree as usize;
    for (k, &m) in initial.iter().enumerate() {
        v[k] = m << (SOBOL_BITS - 1 - k);
    }
    for k in s..SOBOL_BITS {
        let mut x = v[k - s] ^ (v[k - s] >> s);
        for j in 1..s {
            if (coeff >> (s - 1 - j)) & 1 == 1 {
                x ^= v[k - j];
            }
        }
        v[k] = x;
    }
    Ok(v)
}

/// Integer state `x_n` of one Sobol coordinate at point index `n` (Gray-code construction).
///
/// `x_n` is the XOR of the direction numbers selected by the set bits of the
/// Gray code `n ^ (n >> 1)`; it is the value reached after `n` incremental
/// updates `x_i = x_(i-1) ^ v[ctz(i)]`, but is obtained in `O(32)` operations.
pub fn sobol_state_at(directions: &[u32; SOBOL_BITS], index: u32) -> u32 {
    let mut gray = index ^ (index >> 1);
    let mut state = 0u32;
    let mut bit = 0usize;
    while gray != 0 {
        if gray & 1 == 1 {
            state ^= directions[bit];
        }
        gray >>= 1;
        bit += 1;
    }
    state
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Integer states `x_1000` of coordinates `0..128`, cross-checked against SciPy's
    /// `qmc.Sobol(scramble=False, bits=32)` which uses the same Joe-Kuo data.
    const STATES_AT_1000: [u32; 128] = [
        0x38400000, 0x18c00000, 0x84c00000, 0xad400000, 0x47c00000, 0xe8400000, 0x0bc00000,
        0xe6400000, 0x80400000, 0x11c00000, 0x15c00000, 0x41400000, 0x29400000, 0x62400000,
        0x24c00000, 0x5ec00000, 0xb8400000, 0x58400000, 0xfdc00000, 0xb9c00000, 0x85c00000,
        0x8cc00000, 0xf3400000, 0x8a400000, 0x95400000, 0xe8400000, 0x0a400000, 0xfac00000,
        0x0f400000, 0x57400000, 0x25c00000, 0x25400000, 0x4bc00000, 0x97c00000, 0xcd400000,
        0xc5400000, 0xd9400000, 0xd4c00000, 0x4ec00000, 0x7ac00000, 0xe9c00000, 0x41400000,
        0xf5c00000, 0xb3c00000, 0xd7400000, 0x17c00000, 0xdf400000, 0xeb400000, 0xffc00000,
        0x7ac00000, 0x5a400000, 0x84400000, 0xc0c00000, 0x70400000, 0x20400000, 0xe4400000,
        0x74400000, 0x14400000, 0xe2400000, 0x92400000, 0x72400000, 0x9a400000, 0x42400000,
        0x72400000, 0x21400000, 0x7f400000, 0x8f400000, 0x0f400000, 0x31400000, 0xdb400000,
        0x27400000, 0x41400000, 0xe9400000, 0x3d400000, 0x2f400000, 0xe5400000, 0x87c00000,
        0x88c00000, 0xd8c00000, 0x1ec00000, 0x1bc00000, 0xbcc00000, 0xedc00000, 0xbfc00000,
        0x57c00000, 0xfbc00000, 0xf0c00000, 0x7ac00000, 0x95c00000, 0xb4c00000, 0xc4c00000,
        0x88c00000, 0x01c00000, 0x99c00000, 0x5dc00000, 0xddc00000, 0x3dc00000, 0xc3c00000,
        0xd0c00000, 0x2fc00000, 0xc4c00000, 0x8ac00000, 0xe1c00000, 0x8dc00000, 0xc5c00000,
        0x10400000, 0xc1400000, 0xad400000, 0x32400000, 0x1a400000, 0xf7c00000, 0xcbc00000,
        0xf0c00000, 0x9d400000, 0x18400000, 0xdd400000, 0x21c00000, 0x28c00000, 0x09400000,
        0x9b400000, 0x41400000, 0xc6c00000, 0x7f400000, 0x8cc00000, 0x74400000, 0xacc00000,
        0x2b400000, 0x8b400000,
    ];
    /// Integer states `x_0xDEADBEEF` of coordinates `0..128`
    /// (exercises direction numbers up to `v_32`).
    const STATES_AT_DEADBEEF: [u32; 128] = [
        0x1986df8d, 0xfbdd9a3b, 0xba6d6621, 0x56e9df9f, 0x6c52a639, 0x53fd6bd7, 0xf2aca34b,
        0x023e034d, 0x243f5d9d, 0x05bf8ad5, 0x84cbe049, 0x11ca7589, 0xe59f6ffb, 0x34d9d3df,
        0xf13272e1, 0xfb930a89, 0x00cc5afb, 0xb2c5ae51, 0x995849c7, 0xe214688f, 0x7aa46389,
        0x2a601381, 0xaaedb78d, 0xcab38989, 0xaf81fa81, 0x54f61587, 0x8b193587, 0x207a0281,
        0xa4426211, 0x366f76f5, 0x35747699, 0x877a3f37, 0xbff08177, 0x136b2579, 0xf2f2261b,
        0x2dbaccb5, 0x6392eff1, 0xc8eb2869, 0xa95dc1cd, 0x000674d5, 0xec52afa1, 0x2b2a1985,
        0x77b21dbf, 0x74d97e5d, 0x6e5340ed, 0x7d4ebdb9, 0xc99b929d, 0x22a14b69, 0xb82bb6e7,
        0xb2801df3, 0x1e6837dd, 0xfe7aca0b, 0xd7055b17, 0x82556bf1, 0x273518f7, 0xc171ce4d,
        0x1cfa5f95, 0x6614d437, 0x837f0f1d, 0xdf0ce887, 0x43dd80f3, 0x7c3e4a23, 0x56a9111b,
        0x6907e8c5, 0x96631b2d, 0x8f143d93, 0x6fe56439, 0x117ea2bf, 0x3769eb71, 0x95db8ed7,
        0xbd79ce1d, 0x48f7ae11, 0x51ee3167, 0xfd8f7fe1, 0xbfa5dfa1, 0x3ea667e1, 0xb136d155,
        0x5e3689c9, 0xae43c089, 0xf2e652e3, 0x529fb093, 0xc220d445, 0xadf61245, 0xdc800975,
        0xae41cbd1, 0x5f3e2689, 0x29f7ea1d, 0xd65ca1bb, 0xeb541e21, 0xdcd73013, 0x2f09d4c1,
        0x17fa1ef5, 0x04e0bb63, 0x18a9feaf, 0x66d84e61, 0xb648a0ef, 0x5b212453, 0xf0d89a81,
        0xb2c25ca9, 0xc3df88a9, 0xbceac54b, 0x622f34d7, 0x6b20c765, 0x7b6de5f7, 0x5b0ce265,
        0x8760e45f, 0x6d15e41f, 0x48eaa415, 0xa6aead65, 0x43802dd7, 0xceaf548f, 0xc9d1b5ed,
        0xf88c8b7f, 0x018d4b45, 0xcf8c74fb, 0x96ada32b, 0x3b3d4859, 0x0d51f9ab, 0x9d0cd259,
        0xfc257129, 0x8278a40b, 0xef5c841b, 0x4b330ab1, 0xceadc4b3, 0xb5ff1ff1, 0x9792d053,
        0x23d28903, 0xad4ab4e3,
    ];
    /// Integer states `x_0x7FFFFFFF` of coordinates `0..128`.
    const STATES_AT_7FFFFFFF: [u32; 128] = [
        0x00000002, 0xaaaaaaaa, 0x8e00eeee, 0x200200a2, 0x200fb67a, 0x80d6d6d6, 0x8020f866,
        0xcc95a082, 0x880a04a2, 0xb992e922, 0xdfa0041a, 0xce06a74a, 0xcdb6077a, 0xbee00042,
        0xd52d529a, 0xfd836266, 0xe186b2d6, 0xfb6db616, 0x180b600e, 0x2e00000e, 0xea02000e,
        0xa3a20006, 0x0788200a, 0x25db2002, 0xdeeb5e02, 0x31606a0a, 0xd1ab1e0a, 0x1c235a0e,
        0x2f61b32e, 0x341aa9e6, 0x37e26aea, 0x45ea452e, 0x08011262, 0xc959cd6a, 0x9603b36e,
        0x7ce10b6e, 0x976ad5a6, 0x56200202, 0x1d56574e, 0x52524a6a, 0x9758b816, 0x5d48b92a,
        0xe4f2037a, 0x6cf6880e, 0xf272c872, 0x9f22d8ea, 0x839a1d9a, 0x8fb24812, 0x9a000bce,
        0xe3669986, 0x34a09b66, 0x1a623a62, 0xa380038e, 0xd948174a, 0x91fa26d6, 0x81c61cea,
        0x4e6db8ea, 0xdef828b6, 0x4ca01346, 0x1895884e, 0x37563f3e, 0xba6b85f6, 0xcd88886a,
        0x30f9866a, 0x338e8d1e, 0xd58f8dde, 0xfb8816e2, 0x8855acf2, 0x5d3cb27a, 0x44dd8bb2,
        0xd41dbb3a, 0x338e209e, 0x4213a7d6, 0xc8bfa24e, 0x81c826f2, 0x49a58de6, 0xd9680e32,
        0xd4da8ee6, 0x589aa672, 0x9f75a0d2, 0xd2660c06, 0xb78a1992, 0x2c7233f6, 0x39262a36,
        0xba64be4e, 0x42158956, 0x4d603b52, 0xc8003d0a, 0xf335b68a, 0x9f5812ee, 0xfb883062,
        0xcb889c92, 0x76481612, 0x90e083aa, 0x299628a6, 0x01e19caa, 0x988a2de6, 0x25a586fe,
        0x7eca125a, 0xc2862dc6, 0x7642a30a, 0x48f28702, 0xd69f6b82, 0x66da0fc2, 0x0c8180c2,
        0xea000262, 0x81738722, 0xadeee532, 0x6c8ca172, 0xa3e0ef92, 0xb2800a12, 0x0d4a66d2,
        0x16620132, 0x1feb8652, 0xfdf8af4a, 0x30da8d1a, 0xd819288a, 0xd37a00ba, 0x93468e8a,
        0x80e528ea, 0xc698088a, 0x23d8e7fa, 0x0272862a, 0xcb3e079a, 0xc0f3eaaa, 0x15b884da,
        0x7e828b2a, 0x65e3853a,
    ];

    #[test]
    fn table_covers_at_least_21_dimensions() {
        assert_eq!(SOBOL_MAX_DIM, 128);
        const { assert!(SOBOL_MAX_DIM >= 21) };
    }

    #[test]
    fn dimension_zero_is_van_der_corput() {
        let v = sobol_directions(0).unwrap();
        for (k, &vk) in v.iter().enumerate() {
            assert_eq!(vk, 1u32 << (31 - k));
        }
    }

    #[test]
    fn dimension_one_uses_polynomial_x_plus_one() {
        // s = 1, a = 0, m_1 = 1: v_k = v_(k-1) ^ (v_(k-1) >> 1).
        let v = sobol_directions(1).unwrap();
        assert_eq!(v[0], 0x8000_0000);
        assert_eq!(v[1], 0xC000_0000);
        assert_eq!(v[2], 0xA000_0000);
        assert_eq!(v[3], 0xF000_0000);
        assert_eq!(v[4], 0x8800_0000);
    }

    #[test]
    fn dimension_two_uses_joe_kuo_row_three() {
        // d = 3 in the Joe-Kuo file: s = 2, a = 1, m = (1, 3).
        let v = sobol_directions(2).unwrap();
        assert_eq!(v[0], 0x8000_0000);
        assert_eq!(v[1], 0xC000_0000);
        assert_eq!(v[2], 0x6000_0000);
        assert_eq!(v[3], 0x9000_0000);
    }

    #[test]
    fn direction_numbers_are_distinct_between_dimensions() {
        let all: Vec<[u32; SOBOL_BITS]> = (0..SOBOL_MAX_DIM)
            .map(|d| sobol_directions(d).unwrap())
            .collect();
        for i in 0..all.len() {
            for j in (i + 1)..all.len() {
                assert_ne!(
                    all[i], all[j],
                    "dimensions {i} and {j} share direction numbers"
                );
            }
        }
    }

    #[test]
    fn direction_numbers_are_valid_sobol_integers() {
        // v_k = m_k * 2^(32-k) with m_k odd and below 2^k, i.e. the lowest set bit of v_k
        // is bit (32 - k).
        for dim in 0..SOBOL_MAX_DIM {
            let v = sobol_directions(dim).unwrap();
            for (k, &vk) in v.iter().enumerate() {
                assert_eq!(
                    vk.trailing_zeros() as usize,
                    SOBOL_BITS - 1 - k,
                    "dim {dim}, v_{} has the wrong lowest set bit",
                    k + 1
                );
            }
        }
    }

    #[test]
    fn dimension_beyond_table_is_an_error() {
        assert!(sobol_directions(SOBOL_MAX_DIM - 1).is_ok());
        assert_eq!(
            sobol_directions(SOBOL_MAX_DIM),
            Err(PhysicsError::OutOfBounds {
                index: SOBOL_MAX_DIM,
                max: SOBOL_MAX_DIM
            })
        );
    }

    #[test]
    fn states_match_reference_for_all_dimensions() {
        for (index, reference) in [
            (1000u32, &STATES_AT_1000),
            (0xDEAD_BEEF, &STATES_AT_DEADBEEF),
            (0x7FFF_FFFF, &STATES_AT_7FFFFFFF),
        ] {
            for (dim, &want) in reference.iter().enumerate() {
                let v = sobol_directions(dim).unwrap();
                assert_eq!(
                    sobol_state_at(&v, index),
                    want,
                    "dim {dim}, index {index:#x}"
                );
            }
        }
    }

    #[test]
    fn gray_code_state_equals_incremental_update() {
        for dim in [0usize, 1, 2, 5, 20, 63, 127] {
            let v = sobol_directions(dim).unwrap();
            let mut x = 0u32;
            for i in 1u32..4096 {
                x ^= v[i.trailing_zeros() as usize];
                assert_eq!(x, sobol_state_at(&v, i), "dim {dim}, index {i}");
            }
        }
    }
}
