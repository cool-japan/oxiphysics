// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Binary layout of the SAC (Seismic Analysis Code) header.
//!
//! The header is 632 bytes: 70 four-byte float words (bytes 0..280), 40
//! four-byte integer words (280..440) and 24 four-byte character words
//! (440..632). Offsets below are the word number times four. They follow the
//! SAC file-format description; for example `B` is float word 5 (byte 20),
//! `NVHDR` is integer word 76 (byte 304) and `NPTS` is integer word 79
//! (byte 316).

/// Size of the SAC header in bytes: 70 float words, 40 integer words and
/// 24 character words of four bytes each.
pub(super) const SAC_HEADER_BYTES: usize = 632;

/// SAC's "undefined" value for float header words.
pub(super) const SAC_UNDEFINED_FLOAT: f32 = -12345.0;

/// SAC's "undefined" value for integer header words.
pub(super) const SAC_UNDEFINED_INT: i32 = -12345;

/// SAC's "undefined" value for character header words (one eight-byte field).
pub(super) const SAC_UNDEFINED_CHARS: &[u8; 8] = b"-12345  ";

/// `IFTYPE` value of an evenly sampled time series (`ITIME`).
pub(super) const SAC_ITIME: i32 = 1;

/// Header versions (`NVHDR`) that share the 632-byte layout read here.
pub(super) const SAC_SUPPORTED_NVHDR: std::ops::RangeInclusive<i32> = 6..=7;

/// Header version written: the 632-byte version 6 layout without the
/// double-precision footer that version 7 appends to the data.
pub(super) const SAC_WRITTEN_NVHDR: i32 = 6;

/// Byte offsets of the SAC header words used here (word number x 4).
///
/// The header is 70 floats (bytes 0..280), 40 integers (280..440) and the
/// character words (440..632); see the SAC file-format description.
pub(super) const DELTA: usize = 0;
pub(super) const DEPMIN: usize = 4;
pub(super) const DEPMAX: usize = 8;
pub(super) const B: usize = 20;
pub(super) const E: usize = 24;
pub(super) const O: usize = 28;
pub(super) const A: usize = 32;
pub(super) const STLA: usize = 124;
pub(super) const STLO: usize = 128;
pub(super) const STEL: usize = 132;
pub(super) const EVLA: usize = 140;
pub(super) const EVLO: usize = 144;
pub(super) const EVDP: usize = 152;
pub(super) const MAG: usize = 156;
pub(super) const DIST: usize = 200;
pub(super) const AZ: usize = 204;
pub(super) const BAZ: usize = 208;
pub(super) const CMPAZ: usize = 228;
pub(super) const CMPINC: usize = 232;
/// First integer word (`NZYEAR`).
pub(super) const INT_START: usize = 280;
pub(super) const NZYEAR: usize = 280;
pub(super) const NZJDAY: usize = 284;
pub(super) const NZHOUR: usize = 288;
pub(super) const NZMIN: usize = 292;
pub(super) const NZSEC: usize = 296;
pub(super) const NZMSEC: usize = 300;
pub(super) const NVHDR: usize = 304;
pub(super) const NPTS: usize = 316;
pub(super) const IFTYPE: usize = 340;
pub(super) const LEVEN: usize = 420;
/// First character word (`KSTNM`).
pub(super) const CHAR_START: usize = 440;
pub(super) const KSTNM: usize = 440;
pub(super) const KCMPNM: usize = 600;
pub(super) const KNETWK: usize = 608;
/// Width of the character fields read and written here.
pub(super) const CHAR_FIELD_BYTES: usize = 8;

/// Bounds-checked access to the words of a SAC header in either byte order.
pub(super) struct SacWords<'a> {
    pub(super) bytes: &'a [u8],
    pub(super) big_endian: bool,
}

impl SacWords<'_> {
    fn word(&self, offset: usize) -> crate::Result<[u8; 4]> {
        self.bytes
            .get(offset..offset + 4)
            .and_then(|slice| <[u8; 4]>::try_from(slice).ok())
            .ok_or_else(|| {
                crate::Error::Parse(format!("SAC header word at byte {offset} is out of range"))
            })
    }

    pub(super) fn i32_at(&self, offset: usize) -> crate::Result<i32> {
        let word = self.word(offset)?;
        Ok(if self.big_endian {
            i32::from_be_bytes(word)
        } else {
            i32::from_le_bytes(word)
        })
    }

    pub(super) fn f32_at(&self, offset: usize) -> crate::Result<f32> {
        let word = self.word(offset)?;
        Ok(if self.big_endian {
            f32::from_be_bytes(word)
        } else {
            f32::from_le_bytes(word)
        })
    }

    /// An eight-byte character field with trailing blanks and NULs removed.
    pub(super) fn string_at(&self, offset: usize) -> crate::Result<String> {
        let field = self
            .bytes
            .get(offset..offset + CHAR_FIELD_BYTES)
            .ok_or_else(|| {
                crate::Error::Parse(format!("SAC header field at byte {offset} is out of range"))
            })?;
        let text = String::from_utf8_lossy(field);
        Ok(text.trim_end_matches([' ', '\0']).to_string())
    }
}

/// Write `value` into an eight-byte SAC character field, padded with blanks.
///
/// Longer strings are truncated and non-ASCII characters become `?`.
pub(super) fn write_sac_string(buf: &mut [u8], offset: usize, value: &str) {
    let field = &mut buf[offset..offset + CHAR_FIELD_BYTES];
    field.fill(b' ');
    for (slot, ch) in field.iter_mut().zip(value.chars()) {
        *slot = if ch.is_ascii() && ch != '\0' {
            ch as u8
        } else {
            b'?'
        };
    }
}
