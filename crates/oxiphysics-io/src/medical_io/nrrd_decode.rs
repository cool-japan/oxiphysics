// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Decoding of the data section of an NRRD file.
//!
//! [`NrrdReader::parse_header`] reads the header; [`NrrdReader::load_data`] turns the
//! bytes that follow it into `f64` values according to the header's `type`, `endian`,
//! `encoding` and `byte skip` fields.  The `raw`, `text` and `gzip` encodings are
//! decoded; `hex` and `bzip2` are rejected with an error.

use std::io::Read;

use oxiarc_deflate::InflateReader;

use super::types::{NrrdEncoding, NrrdReader};
use crate::{Error, Result};

/// Largest number of bytes a data section may decode to by default (1 GiB).
const DEFAULT_MAX_DECODED_BYTES: u64 = 1 << 30;

/// Element type of an NRRD data set.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum NrrdType {
    Int8,
    Uint8,
    Int16,
    Uint16,
    Int32,
    Uint32,
    Int64,
    Uint64,
    Float32,
    Float64,
}

impl NrrdType {
    /// Type named by an NRRD `type:` field, with all the spellings the format allows.
    fn from_name(name: &str) -> Option<Self> {
        let words: Vec<&str> = name.split_whitespace().collect();
        let normalized = words.join(" ").to_ascii_lowercase();
        Some(match normalized.as_str() {
            "signed char" | "int8" | "int8_t" => Self::Int8,
            "uchar" | "unsigned char" | "uint8" | "uint8_t" => Self::Uint8,
            "short" | "short int" | "signed short" | "signed short int" | "int16" | "int16_t" => {
                Self::Int16
            }
            "ushort" | "unsigned short" | "unsigned short int" | "uint16" | "uint16_t" => {
                Self::Uint16
            }
            "int" | "signed int" | "int32" | "int32_t" => Self::Int32,
            "uint" | "unsigned int" | "uint32" | "uint32_t" => Self::Uint32,
            "longlong"
            | "long long"
            | "long long int"
            | "signed long long"
            | "signed long long int"
            | "int64"
            | "int64_t" => Self::Int64,
            "ulonglong"
            | "unsigned long long"
            | "unsigned long long int"
            | "uint64"
            | "uint64_t" => Self::Uint64,
            "float" => Self::Float32,
            "double" => Self::Float64,
            _ => return None,
        })
    }

    /// Size of one element in bytes.
    fn size(self) -> usize {
        match self {
            Self::Int8 | Self::Uint8 => 1,
            Self::Int16 | Self::Uint16 => 2,
            Self::Int32 | Self::Uint32 | Self::Float32 => 4,
            Self::Int64 | Self::Uint64 | Self::Float64 => 8,
        }
    }

    /// Value of the element stored in `chunk`, which holds exactly [`Self::size`] bytes.
    fn decode(self, chunk: &[u8], big_endian: bool) -> f64 {
        macro_rules! number {
            ($ty:ty, $width:literal, $to_f64:expr) => {{
                let bytes = take::<$width>(chunk);
                let value = if big_endian {
                    <$ty>::from_be_bytes(bytes)
                } else {
                    <$ty>::from_le_bytes(bytes)
                };
                $to_f64(value)
            }};
        }
        match self {
            Self::Int8 => number!(i8, 1, f64::from),
            Self::Uint8 => number!(u8, 1, f64::from),
            Self::Int16 => number!(i16, 2, f64::from),
            Self::Uint16 => number!(u16, 2, f64::from),
            Self::Int32 => number!(i32, 4, f64::from),
            Self::Uint32 => number!(u32, 4, f64::from),
            Self::Int64 => number!(i64, 8, |v: i64| v as f64),
            Self::Uint64 => number!(u64, 8, |v: u64| v as f64),
            Self::Float32 => number!(f32, 4, f64::from),
            Self::Float64 => number!(f64, 8, |v: f64| v),
        }
    }
}

/// The first `N` bytes of `chunk` (zero padded when it is shorter).
fn take<const N: usize>(chunk: &[u8]) -> [u8; N] {
    let mut out = [0u8; N];
    for (dst, src) in out.iter_mut().zip(chunk) {
        *dst = *src;
    }
    out
}

/// How much of the decoded stream is skipped before the data (`byte skip`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum ByteSkip {
    /// Skip this many bytes from the start.
    Front(usize),
    /// `byte skip: -1`: the data is the last `n_elements * element_size` bytes.
    Tail,
}

/// Decompress a gzip stream that may not produce more than `max_output` bytes.
fn gunzip_bounded(data: &[u8], max_output: u64) -> Result<Vec<u8>> {
    let mut reader = InflateReader::gzip(data).with_max_output(max_output);
    let mut out = Vec::new();
    reader.read_to_end(&mut out).map_err(|e| {
        Error::Parse(format!(
            "NRRD gzip data cannot be decompressed within {max_output} bytes: {e}"
        ))
    })?;
    Ok(out)
}

fn too_large() -> Error {
    Error::Parse("NRRD sizes describe more data than can be addressed".to_string())
}

impl NrrdReader {
    /// Header value of `key` (keys are stored in lower case).
    fn header_value(&self, key: &str) -> Option<&str> {
        self.fields.get(key).map(String::as_str)
    }

    /// Element type of the data, from the `type` field.
    fn element_type(&self) -> Result<NrrdType> {
        let name = self
            .header_value("type")
            .ok_or_else(|| Error::Parse("NRRD header has no 'type' field".to_string()))?;
        NrrdType::from_name(name).ok_or_else(|| {
            Error::Parse(format!(
                "NRRD type '{name}' is not a numeric type this reader can decode"
            ))
        })
    }

    /// Whether multi-byte values are big endian, from the `endian` field.
    fn big_endian(&self, ty: NrrdType) -> Result<bool> {
        match self.header_value("endian").map(str::to_ascii_lowercase) {
            Some(value) if value == "little" => Ok(false),
            Some(value) if value == "big" => Ok(true),
            Some(value) => Err(Error::Parse(format!(
                "NRRD endian '{value}' is neither 'little' nor 'big'"
            ))),
            None if ty.size() == 1 => Ok(false),
            None => Err(Error::Parse(
                "NRRD header has no 'endian' field for a multi-byte type".to_string(),
            )),
        }
    }

    /// Number of elements the `sizes` field describes.
    fn element_count(&self) -> Result<usize> {
        if self.sizes.is_empty() {
            return Err(Error::Parse("NRRD header has no 'sizes' field".to_string()));
        }
        if let Some(dimension) = self.header_value("dimension") {
            let declared: usize = dimension.trim().parse().map_err(|_| {
                Error::Parse(format!("NRRD dimension '{dimension}' is not a number"))
            })?;
            if declared != self.sizes.len() {
                return Err(Error::Parse(format!(
                    "NRRD dimension is {declared} but sizes lists {} axes",
                    self.sizes.len()
                )));
            }
        }
        self.sizes
            .iter()
            .try_fold(1usize, |count, &size| count.checked_mul(size))
            .ok_or_else(too_large)
    }

    /// The `byte skip` field (`byteskip` is the same field); absent means 0.
    fn byte_skip(&self) -> Result<ByteSkip> {
        let Some(text) = self
            .header_value("byte skip")
            .or_else(|| self.header_value("byteskip"))
        else {
            return Ok(ByteSkip::Front(0));
        };
        let value: i64 = text
            .trim()
            .parse()
            .map_err(|_| Error::Parse(format!("NRRD byte skip '{text}' is not an integer")))?;
        match value {
            -1 => Ok(ByteSkip::Tail),
            n if n >= 0 => usize::try_from(n)
                .map(ByteSkip::Front)
                .map_err(|_| too_large()),
            n => Err(Error::Parse(format!(
                "NRRD byte skip {n} is invalid (only -1 or a non-negative count)"
            ))),
        }
    }

    /// Reject headers that ask for something this reader does not do.
    fn check_supported(&self) -> Result<()> {
        if let Some(name) = self.header_value("encoding")
            && NrrdEncoding::from_name(name).is_none()
        {
            return Err(Error::Parse(format!(
                "NRRD encoding '{name}' is not supported (raw, text and gzip are)"
            )));
        }
        let line_skip = self
            .header_value("line skip")
            .or_else(|| self.header_value("lineskip"));
        if let Some(text) = line_skip
            && !matches!(text.trim().parse::<i64>(), Ok(0))
        {
            return Err(Error::Parse(format!(
                "NRRD line skip '{text}' is not supported"
            )));
        }
        Ok(())
    }

    /// Decode the data section of an NRRD file into [`NrrdReader::data`].
    ///
    /// Call [`NrrdReader::parse_header`] first: the header's `type`, `sizes`,
    /// `endian`, `encoding` and `byte skip` fields describe `bytes`, which are the
    /// bytes after the blank line that ends the header (or the contents of the
    /// detached data file).  The decoded values are converted to `f64`.
    ///
    /// * `raw` and `gzip` data hold `sizes.product()` elements of the header's
    ///   `type` in the header's `endian` order; `gzip` data is inflated first.
    /// * `text` data holds that many whitespace-separated numbers.
    /// * `byte skip: N` drops the first `N` (decompressed) bytes, and `byte skip: -1`
    ///   takes the data from the end of the (decompressed) stream.
    ///
    /// A data section may decode to at most 1 GiB; see
    /// [`NrrdReader::load_data_limited`] for another bound.
    ///
    /// # Errors
    ///
    /// [`Error::Parse`] when the header lacks `type` / `sizes` (or `endian` for a
    /// multi-byte type), names an unsupported type or encoding (`block`, `hex`,
    /// `bzip2`), asks for a `line skip`, or when the data is corrupt, larger than the
    /// limit, or does not hold exactly the number of elements `sizes` describes.
    pub fn load_data(&mut self, bytes: &[u8]) -> Result<()> {
        self.load_data_limited(bytes, DEFAULT_MAX_DECODED_BYTES)
    }

    /// [`NrrdReader::load_data`] with an explicit bound on the number of bytes the
    /// data (before any `byte skip`) may decode to.
    ///
    /// The bound also stops a gzip stream that expands far beyond what the header
    /// declares, so untrusted files cannot exhaust memory.
    ///
    /// # Errors
    ///
    /// The errors of [`NrrdReader::load_data`], and [`Error::Parse`] when the
    /// declared data or the `byte skip` exceeds `max_decoded_bytes`.
    pub fn load_data_limited(&mut self, bytes: &[u8], max_decoded_bytes: u64) -> Result<()> {
        self.check_supported()?;
        let ty = self.element_type()?;
        let count = self.element_count()?;
        let payload = count.checked_mul(ty.size()).ok_or_else(too_large)?;
        if payload as u64 > max_decoded_bytes {
            return Err(Error::Parse(format!(
                "NRRD data of {payload} bytes exceeds the limit of {max_decoded_bytes} bytes"
            )));
        }
        let skip = self.byte_skip()?;
        if let ByteSkip::Front(n) = skip
            && n as u64 > max_decoded_bytes
        {
            return Err(Error::Parse(format!(
                "NRRD byte skip {n} exceeds the limit of {max_decoded_bytes} bytes"
            )));
        }
        let values = match self.encoding {
            NrrdEncoding::Text => decode_text(bytes, skip, count)?,
            NrrdEncoding::Raw => decode_binary(bytes, skip, ty, self.big_endian(ty)?, payload)?,
            NrrdEncoding::Gzip => {
                let big_endian = self.big_endian(ty)?;
                let limit = match skip {
                    ByteSkip::Front(n) => (n as u64).saturating_add(payload as u64),
                    ByteSkip::Tail => max_decoded_bytes,
                };
                let inflated = gunzip_bounded(bytes, limit)?;
                decode_binary(&inflated, skip, ty, big_endian, payload)?
            }
        };
        self.data = values;
        Ok(())
    }

    /// Load gzip-compressed raw little-endian `f32` data, the compressed counterpart
    /// of [`NrrdReader::load_raw_f32`].
    ///
    /// The stream may inflate to at most 1 GiB.
    ///
    /// # Errors
    ///
    /// [`Error::Parse`] when `bytes` is not a valid gzip stream or inflates beyond the
    /// limit.
    pub fn load_gzip_f32(&mut self, bytes: &[u8]) -> Result<()> {
        let inflated = gunzip_bounded(bytes, DEFAULT_MAX_DECODED_BYTES)?;
        self.load_raw_f32(&inflated);
        Ok(())
    }
}

/// Decode `payload` bytes of `ty` elements from a decompressed or raw stream.
fn decode_binary(
    stream: &[u8],
    skip: ByteSkip,
    ty: NrrdType,
    big_endian: bool,
    payload: usize,
) -> Result<Vec<f64>> {
    let body = match skip {
        ByteSkip::Front(n) => stream.get(n..),
        ByteSkip::Tail => stream
            .len()
            .checked_sub(payload)
            .and_then(|start| stream.get(start..)),
    }
    .ok_or_else(|| {
        Error::Parse(format!(
            "NRRD data section of {} bytes is too short for the byte skip and {payload} bytes of data",
            stream.len()
        ))
    })?;
    if body.len() != payload {
        return Err(Error::Parse(format!(
            "NRRD data section holds {} bytes but sizes and type describe {payload}",
            body.len()
        )));
    }
    Ok(body
        .chunks_exact(ty.size())
        .map(|chunk| ty.decode(chunk, big_endian))
        .collect())
}

/// Parse `count` whitespace-separated numbers from text data.
fn decode_text(bytes: &[u8], skip: ByteSkip, count: usize) -> Result<Vec<f64>> {
    let ByteSkip::Front(skipped) = skip else {
        return Err(Error::Parse(
            "NRRD byte skip -1 is only valid for raw and gzip data".to_string(),
        ));
    };
    let body = bytes.get(skipped..).ok_or_else(|| {
        Error::Parse(format!(
            "NRRD byte skip {skipped} is beyond the {} bytes of data",
            bytes.len()
        ))
    })?;
    let text = std::str::from_utf8(body)
        .map_err(|e| Error::Parse(format!("NRRD text data is not valid UTF-8: {e}")))?;
    let values = text
        .split_whitespace()
        .map(|token| {
            token
                .parse::<f64>()
                .map_err(|_| Error::Parse(format!("NRRD text value '{token}' is not a number")))
        })
        .collect::<Result<Vec<f64>>>()?;
    if values.len() != count {
        return Err(Error::Parse(format!(
            "NRRD text data holds {} values but sizes describe {count}",
            values.len()
        )));
    }
    Ok(values)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn header(extra: &str) -> NrrdReader {
        let mut reader = NrrdReader::new();
        reader.parse_header(&format!("NRRD0004\n{extra}\n"));
        reader
    }

    fn le_f32(values: &[f32]) -> Vec<u8> {
        values.iter().flat_map(|v| v.to_le_bytes()).collect()
    }

    fn gzip(bytes: &[u8]) -> Vec<u8> {
        oxiarc_deflate::gzip_compress(bytes, 6).unwrap()
    }

    fn message(result: Result<()>) -> String {
        match result {
            Err(Error::Parse(m)) => m,
            other => panic!("expected a parse error, got {other:?}"),
        }
    }

    #[test]
    fn gzip_float_volume_is_decoded() {
        let values: Vec<f32> = (0..24).map(|i| i as f32 * 0.5 - 3.0).collect();
        let mut reader =
            header("type: float\ndimension: 3\nsizes: 4 3 2\nendian: little\nencoding: gzip");
        assert_eq!(reader.encoding, NrrdEncoding::Gzip);
        reader.load_data(&gzip(&le_f32(&values))).unwrap();
        assert_eq!(reader.data.len(), 24);
        for (got, want) in reader.data.iter().zip(&values) {
            assert_eq!(*got, f64::from(*want));
        }
        assert_eq!(reader.n_voxels(), 24);
    }

    #[test]
    fn gzip_matches_raw_decoding_for_every_type() {
        let cases: [(&str, usize); 10] = [
            ("signed char", 1),
            ("uchar", 1),
            ("short", 2),
            ("ushort", 2),
            ("int", 4),
            ("uint", 4),
            ("longlong", 8),
            ("ulonglong", 8),
            ("float", 4),
            ("double", 8),
        ];
        for (name, width) in cases {
            let bytes: Vec<u8> = (0..(6 * width) as u8).map(|b| b.wrapping_mul(7)).collect();
            for endian in ["little", "big"] {
                let spec = format!("type: {name}\nsizes: 6\nendian: {endian}");
                let mut raw = header(&format!("{spec}\nencoding: raw"));
                raw.load_data(&bytes).unwrap();
                let mut gz = header(&format!("{spec}\nencoding: gzip"));
                gz.load_data(&gzip(&bytes)).unwrap();
                assert_eq!(raw.data.len(), 6, "{name} {endian}");
                for (a, b) in raw.data.iter().zip(&gz.data) {
                    assert!(a == b || (a.is_nan() && b.is_nan()), "{name} {endian}");
                }
            }
        }
    }

    #[test]
    fn integer_types_honour_signedness_and_endianness() {
        let mut reader = header("type: short\nsizes: 3\nendian: big\nencoding: raw");
        reader
            .load_data(&[0xff, 0xfe, 0x00, 0x01, 0x80, 0x00])
            .unwrap();
        assert_eq!(reader.data, vec![-2.0, 1.0, -32768.0]);

        let mut reader = header("type: uint16\nsizes: 3\nendian: little\nencoding: raw");
        reader
            .load_data(&[0xfe, 0xff, 0x01, 0x00, 0x00, 0x80])
            .unwrap();
        assert_eq!(reader.data, vec![65534.0, 1.0, 32768.0]);

        let mut reader = header("type: int8\nsizes: 2\nencoding: raw");
        reader.load_data(&[0xff, 0x7f]).unwrap();
        assert_eq!(reader.data, vec![-1.0, 127.0]);

        let mut reader = header("type: int64\nsizes: 1\nendian: big\nencoding: raw");
        reader
            .load_data(&(-5_000_000_000i64).to_be_bytes())
            .unwrap();
        assert_eq!(reader.data, vec![-5_000_000_000.0]);

        let mut reader = header("type: double\nsizes: 2\nendian: little\nencoding: raw");
        let mut bytes = 1.25f64.to_le_bytes().to_vec();
        bytes.extend_from_slice(&(-0.5f64).to_le_bytes());
        reader.load_data(&bytes).unwrap();
        assert_eq!(reader.data, vec![1.25, -0.5]);
    }

    #[test]
    fn text_encoding_accepts_txt_and_checks_the_count() {
        let mut reader = header("type: float\nsizes: 2 2\nencoding: txt");
        assert_eq!(reader.encoding, NrrdEncoding::Text);
        reader.load_data(b"1.5 -2\n3e2 4\n").unwrap();
        assert_eq!(reader.data, vec![1.5, -2.0, 300.0, 4.0]);

        let m = message(reader.load_data(b"1 2 3"));
        assert!(m.contains("3 values") && m.contains("4"), "{m}");
        let m = message(reader.load_data(b"1 2 three 4"));
        assert!(m.contains("three"), "{m}");
        let m = message(reader.load_data(&[0xff, 0xfe]));
        assert!(m.contains("UTF-8"), "{m}");
    }

    #[test]
    fn byte_skip_drops_leading_bytes_or_selects_the_tail() {
        let values = [1.0f32, 2.0, 3.0];
        let mut payload = vec![0xAA; 5];
        payload.extend(le_f32(&values));

        let mut reader =
            header("type: float\nsizes: 3\nendian: little\nencoding: raw\nbyte skip: 5");
        reader.load_data(&payload).unwrap();
        assert_eq!(reader.data, vec![1.0, 2.0, 3.0]);

        let mut reader =
            header("type: float\nsizes: 3\nendian: little\nencoding: raw\nbyte skip: -1");
        reader.load_data(&payload).unwrap();
        assert_eq!(reader.data, vec![1.0, 2.0, 3.0]);

        let mut reader =
            header("type: float\nsizes: 3\nendian: little\nencoding: gzip\nbyteskip: 5");
        reader.load_data(&gzip(&payload)).unwrap();
        assert_eq!(reader.data, vec![1.0, 2.0, 3.0]);

        let mut reader =
            header("type: float\nsizes: 3\nendian: little\nencoding: gzip\nbyte skip: -1");
        reader.load_data(&gzip(&payload)).unwrap();
        assert_eq!(reader.data, vec![1.0, 2.0, 3.0]);

        let mut reader = header("type: float\nsizes: 2\nencoding: text\nbyte skip: 3");
        reader.load_data(b"xx 4 5").unwrap();
        assert_eq!(reader.data, vec![4.0, 5.0]);

        let mut reader =
            header("type: float\nsizes: 3\nendian: little\nencoding: raw\nbyte skip: -2");
        assert!(message(reader.load_data(&payload)).contains("-2"));
        let mut reader = header("type: float\nsizes: 2\nencoding: text\nbyte skip: -1");
        assert!(message(reader.load_data(b"1 2")).contains("raw and gzip"));
        let mut reader =
            header("type: float\nsizes: 3\nendian: little\nencoding: raw\nbyte skip: 100");
        assert!(message(reader.load_data(&payload)).contains("too short"));
    }

    #[test]
    fn wrong_data_length_is_an_error() {
        let mut reader = header("type: float\nsizes: 4\nendian: little\nencoding: raw");
        let m = message(reader.load_data(&le_f32(&[1.0, 2.0, 3.0])));
        assert!(m.contains("12 bytes") && m.contains("16"), "{m}");
        let m = message(reader.load_data(&le_f32(&[1.0, 2.0, 3.0, 4.0, 5.0])));
        assert!(m.contains("20 bytes"), "{m}");
        assert!(reader.data.is_empty(), "failed loads must not change data");

        let mut gz = header("type: float\nsizes: 4\nendian: little\nencoding: gzip");
        let m = message(gz.load_data(&gzip(&le_f32(&[1.0, 2.0, 3.0]))));
        assert!(m.contains("12 bytes"), "{m}");
    }

    #[test]
    fn header_problems_are_reported() {
        let data = le_f32(&[1.0]);
        let cases = [
            ("sizes: 1\nendian: little\nencoding: raw", "'type'"),
            ("type: float\nendian: little\nencoding: raw", "'sizes'"),
            ("type: float\nsizes: 1\nencoding: raw", "'endian'"),
            (
                "type: float\nsizes: 1\nendian: middle\nencoding: raw",
                "middle",
            ),
            (
                "type: block\nsizes: 1\nendian: little\nencoding: raw",
                "block",
            ),
            (
                "type: float\nsizes: 1\nendian: little\nencoding: bzip2",
                "bzip2",
            ),
            (
                "type: float\nsizes: 1\nendian: little\nencoding: hex",
                "hex",
            ),
            (
                "type: float\nsizes: 1\nendian: little\nencoding: raw\nline skip: 2",
                "line skip",
            ),
            (
                "type: float\nsizes: 1\nendian: little\nencoding: raw\nbyte skip: x",
                "byte skip",
            ),
            (
                "type: float\ndimension: 2\nsizes: 1\nendian: little\nencoding: raw",
                "dimension",
            ),
            (
                "type: float\ndimension: two\nsizes: 1\nendian: little\nencoding: raw",
                "dimension",
            ),
        ];
        for (text, expected) in cases {
            let mut reader = header(text);
            let m = message(reader.load_data(&data));
            assert!(m.contains(expected), "{text}: {m}");
        }
    }

    #[test]
    fn sizes_that_overflow_or_exceed_the_limit_are_errors_not_allocations() {
        let mut reader =
            header("type: double\nsizes: 4294967296 4294967296 4\nendian: little\nencoding: raw");
        assert!(message(reader.load_data(&[])).contains("addressed"));

        let mut reader =
            header("type: double\nsizes: 100000 100000\nendian: little\nencoding: gzip");
        let m = message(reader.load_data(&gzip(&[0u8; 16])));
        assert!(m.contains("limit"), "{m}");

        let mut reader = header("type: float\nsizes: 6\nendian: little\nencoding: raw");
        let m = message(reader.load_data_limited(&le_f32(&[0.0; 6]), 16));
        assert!(m.contains("24 bytes") && m.contains("16"), "{m}");
        reader.load_data_limited(&le_f32(&[0.0; 6]), 24).unwrap();
        assert_eq!(reader.data.len(), 6);
    }

    #[test]
    fn a_gzip_bomb_is_stopped_at_the_declared_size() {
        // 8 MiB of zeros compresses to a few kilobytes but the header declares 4 floats.
        let bomb = gzip(&vec![0u8; 8 << 20]);
        assert!(bomb.len() < 64 * 1024);
        let mut reader = header("type: float\nsizes: 4\nendian: little\nencoding: gzip");
        let m = message(reader.load_data(&bomb));
        assert!(m.contains("cannot be decompressed within 16 bytes"), "{m}");
    }

    #[test]
    fn corrupt_gzip_data_is_an_error() {
        let mut reader = header("type: float\nsizes: 4\nendian: little\nencoding: gzip");
        let mut stream = gzip(&le_f32(&[1.0, 2.0, 3.0, 4.0]));
        let mid = stream.len() / 2;
        stream[mid] ^= 0xff;
        let last = stream.len() - 1;
        stream[last] ^= 0xff;
        assert!(reader.load_data(&stream).is_err());
        assert!(reader.load_data(b"not gzip at all").is_err());
        assert!(reader.load_data(&[]).is_err());
        assert!(reader.data.is_empty());
    }

    #[test]
    fn load_gzip_f32_decodes_headerless_compressed_floats() {
        let values = [0.25f32, -8.0, 1.0e6];
        let mut reader = NrrdReader::new();
        reader.load_gzip_f32(&gzip(&le_f32(&values))).unwrap();
        assert_eq!(reader.data, vec![0.25, -8.0, 1.0e6]);
        assert!(reader.load_gzip_f32(b"plain bytes").is_err());
    }

    #[test]
    fn encoding_names_map_to_encodings() {
        assert_eq!(NrrdEncoding::from_name("raw"), Some(NrrdEncoding::Raw));
        assert_eq!(NrrdEncoding::from_name("gz"), Some(NrrdEncoding::Gzip));
        assert_eq!(NrrdEncoding::from_name(" GZIP "), Some(NrrdEncoding::Gzip));
        for name in ["text", "txt", "ascii"] {
            assert_eq!(NrrdEncoding::from_name(name), Some(NrrdEncoding::Text));
        }
        for name in ["hex", "bzip2", "bz2", "", "gzip2"] {
            assert_eq!(NrrdEncoding::from_name(name), None, "{name}");
        }
    }

    #[test]
    fn spelled_out_type_names_are_understood() {
        for (name, ty) in [
            ("unsigned char", NrrdType::Uint8),
            ("uint8_t", NrrdType::Uint8),
            ("signed short int", NrrdType::Int16),
            ("unsigned   short", NrrdType::Uint16),
            ("int32_t", NrrdType::Int32),
            ("unsigned int", NrrdType::Uint32),
            ("long long int", NrrdType::Int64),
            ("unsigned long long", NrrdType::Uint64),
            ("FLOAT", NrrdType::Float32),
            ("double", NrrdType::Float64),
        ] {
            assert_eq!(NrrdType::from_name(name), Some(ty), "{name}");
        }
        assert_eq!(NrrdType::from_name("block"), None);
        assert_eq!(NrrdType::from_name("long"), None);
    }
}
