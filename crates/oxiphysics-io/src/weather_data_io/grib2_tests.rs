// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! GRIB2 reader tests.
//!
//! Every fixture is built octet by octet from the WMO tables (FM 92 GRIB
//! edition 2): section 1 identification (21 octets), template 3.0 grid (72
//! octets), 3.20 (65), 3.30 (81), template 4.0 product (34), template 5.0
//! data representation (21), 5.4 (12), 6 bit-map and 7 data. The builder
//! writes by 1-based octet position, so the fixtures do not share any offset
//! arithmetic with the decoder, which reads `payload[octet - 6]`.

use super::*;

// ── fixture builder ──────────────────────────────────────────────────────

/// A section built by 1-based WMO octet position.
struct SectionBuilder {
    bytes: Vec<u8>,
}

impl SectionBuilder {
    fn new(number: u8, length: usize) -> Self {
        let mut bytes = vec![0u8; length];
        bytes[0..4].copy_from_slice(&(length as u32).to_be_bytes());
        bytes[4] = number;
        Self { bytes }
    }

    fn put(mut self, octet: usize, value: &[u8]) -> Self {
        let start = octet - 1;
        self.bytes[start..start + value.len()].copy_from_slice(value);
        self
    }

    fn u8(self, octet: usize, value: u8) -> Self {
        self.put(octet, &[value])
    }

    fn u16(self, octet: usize, value: u16) -> Self {
        self.put(octet, &value.to_be_bytes())
    }

    fn u32(self, octet: usize, value: u32) -> Self {
        self.put(octet, &value.to_be_bytes())
    }

    /// 4-octet sign-magnitude integer.
    fn signed32(self, octet: usize, value: i64) -> Self {
        let magnitude = value.unsigned_abs() as u32;
        let raw = if value < 0 {
            magnitude | 0x8000_0000
        } else {
            magnitude
        };
        self.u32(octet, raw)
    }

    /// 2-octet sign-magnitude integer.
    fn signed16(self, octet: usize, value: i32) -> Self {
        let magnitude = value.unsigned_abs() as u16;
        let raw = if value < 0 {
            magnitude | 0x8000
        } else {
            magnitude
        };
        self.u16(octet, raw)
    }

    fn f32(self, octet: usize, value: f32) -> Self {
        self.put(octet, &value.to_be_bytes())
    }

    fn build(self) -> Vec<u8> {
        self.bytes
    }
}

/// Section 1: originating centre and reference time.
fn identification(centre: u16, date: (u16, u8, u8), time: (u8, u8, u8)) -> Vec<u8> {
    SectionBuilder::new(1, 21)
        .u16(6, centre)
        .u16(8, 0)
        .u8(10, 28)
        .u8(11, 0)
        .u8(12, 1)
        .u16(13, date.0)
        .u8(15, date.1)
        .u8(16, date.2)
        .u8(17, time.0)
        .u8(18, time.1)
        .u8(19, time.2)
        .u8(20, 0)
        .u8(21, 1)
        .build()
}

fn default_identification() -> Vec<u8> {
    identification(98, (2026, 10, 19), (12, 30, 15))
}

/// Template 3.0 parameters; angles are in the unit set by the basic angle
/// and its subdivisions (10^-6 degrees when both are zero).
#[derive(Clone)]
struct LatLonGrid {
    ni: u32,
    nj: u32,
    la1: i64,
    lo1: i64,
    la2: i64,
    lo2: i64,
    di: u32,
    dj: u32,
    flags: u8,
    basic_angle: u32,
    subdivisions: u32,
}

impl LatLonGrid {
    /// A 4 x 3 grid from 20.5 N / 340 E to 4.5 S / 347.5 E in 2.5 / 12.5
    /// degree steps. The negative latitude and the 12.5 degree step make a
    /// two's-complement or mis-offset decoder fail.
    fn standard() -> Self {
        Self {
            ni: 4,
            nj: 3,
            la1: 20_500_000,
            lo1: 340_000_000,
            la2: -4_500_000,
            lo2: 347_500_000,
            di: 2_500_000,
            dj: 12_500_000,
            flags: 0x30,
            basic_angle: 0,
            subdivisions: 0,
        }
    }

    fn section(&self) -> Vec<u8> {
        SectionBuilder::new(3, 72)
            .u8(6, 0)
            .u32(7, self.ni.wrapping_mul(self.nj))
            .u8(11, 0)
            .u8(12, 0)
            .u16(13, 0)
            .u8(15, 6)
            .u32(31, self.ni)
            .u32(35, self.nj)
            .u32(39, self.basic_angle)
            .u32(43, self.subdivisions)
            .signed32(47, self.la1)
            .signed32(51, self.lo1)
            .u8(55, self.flags)
            .signed32(56, self.la2)
            .signed32(60, self.lo2)
            .u32(64, self.di)
            .u32(68, self.dj)
            .u8(72, 0)
            .build()
    }
}

/// Template 3.20, polar stereographic (65 octets).
fn polar_stereographic_grid(nx: u32, ny: u32) -> Vec<u8> {
    SectionBuilder::new(3, 65)
        .u32(7, nx * ny)
        .u16(13, 20)
        .u8(15, 6)
        .u32(31, nx)
        .u32(35, ny)
        .signed32(39, 30_000_000)
        .signed32(43, 279_000_000)
        .u8(47, 0x30)
        .signed32(48, 60_000_000)
        .signed32(52, 255_000_000)
        .u32(56, 25_000_000)
        .u32(60, 25_000_000)
        .u8(64, 0)
        .u8(65, 0x40)
        .build()
}

/// Template 3.30, Lambert conformal; the HRRR 3 km grid (81 octets).
fn lambert_grid() -> Vec<u8> {
    SectionBuilder::new(3, 81)
        .u32(7, 1799 * 1059)
        .u16(13, 30)
        .u8(15, 6)
        .u32(31, 1799)
        .u32(35, 1059)
        .signed32(39, 21_138_123)
        .signed32(43, 237_280_472)
        .u8(47, 0x08)
        .signed32(48, 38_500_000)
        .signed32(52, 262_500_000)
        .u32(56, 3_000_000)
        .u32(60, 3_000_000)
        .u8(64, 0)
        .u8(65, 0x40)
        .signed32(66, 38_500_000)
        .signed32(70, 38_500_000)
        .signed32(74, -90_000_000)
        .signed32(78, 0)
        .build()
}

/// Section 4, template 4.0.
fn product(category: u8, number: u8, unit: u8, forecast: u32) -> Vec<u8> {
    SectionBuilder::new(4, 34)
        .u16(6, 0)
        .u16(8, 0)
        .u8(10, category)
        .u8(11, number)
        .u8(12, 2)
        .u8(13, 0)
        .u8(14, 96)
        .u16(15, 0)
        .u8(17, 0)
        .u8(18, unit)
        .u32(19, forecast)
        .u8(23, 103)
        .u8(24, 0)
        .u32(25, 2)
        .u8(29, 255)
        .u8(30, 255)
        .u32(31, u32::MAX)
        .build()
}

/// Section 4 with an arbitrary template number and otherwise template 4.0 content.
fn product_template(template: u16) -> Vec<u8> {
    SectionBuilder::new(4, 34)
        .u16(8, template)
        .u8(18, 1)
        .build()
}

/// Section 5, template 5.0: `Y = (R + X * 2^E) / 10^D`, `bits` per value.
fn simple_packing(n_values: u32, reference: f32, e: i32, d: i32, bits: u8) -> Vec<u8> {
    SectionBuilder::new(5, 21)
        .u32(6, n_values)
        .u16(10, 0)
        .f32(12, reference)
        .signed16(16, e)
        .signed16(18, d)
        .u8(20, bits)
        .u8(21, 0)
        .build()
}

/// Section 5, template 5.4 (IEEE floating point).
fn ieee_packing(n_values: u32, precision: u8) -> Vec<u8> {
    SectionBuilder::new(5, 12)
        .u32(6, n_values)
        .u16(10, 4)
        .u8(12, precision)
        .build()
}

/// Section 5 with a template this reader does not decode.
fn other_packing(n_values: u32, template: u16) -> Vec<u8> {
    SectionBuilder::new(5, 21)
        .u32(6, n_values)
        .u16(10, template)
        .build()
}

fn bitmap(indicator: u8, bits: &[u8]) -> Vec<u8> {
    SectionBuilder::new(6, 6 + bits.len())
        .u8(6, indicator)
        .put(7, bits)
        .build()
}

fn no_bitmap() -> Vec<u8> {
    bitmap(255, &[])
}

fn data(payload: &[u8]) -> Vec<u8> {
    SectionBuilder::new(7, 5 + payload.len())
        .put(6, payload)
        .build()
}

/// Sections 0 and 8 around `sections`.
fn message(discipline: u8, sections: &[Vec<u8>]) -> Vec<u8> {
    let body: usize = sections.iter().map(Vec::len).sum();
    let mut out = Vec::new();
    out.extend_from_slice(b"GRIB");
    out.extend_from_slice(&[0xFF, 0xFF]);
    out.push(discipline);
    out.push(2);
    out.extend_from_slice(&((16 + body + 4) as u64).to_be_bytes());
    for section in sections {
        out.extend_from_slice(section);
    }
    out.extend_from_slice(b"7777");
    out
}

/// `values` packed MSB first, `width` bits each, zero padded to a byte.
fn pack_msb_first(values: &[u64], width: usize) -> Vec<u8> {
    let mut bytes = vec![0u8; (values.len() * width).div_ceil(8)];
    let mut bit = 0usize;
    for &value in values {
        for shift in (0..width).rev() {
            if (value >> shift) & 1 == 1 {
                bytes[bit / 8] |= 0x80 >> (bit % 8);
            }
            bit += 1;
        }
    }
    bytes
}

/// The 4 x 3 standard grid with 12 simple-packed values of 12 bits.
fn standard_message(values: &[u64; 12], reference: f32, e: i32, d: i32) -> Vec<u8> {
    message(
        0,
        &[
            default_identification(),
            LatLonGrid::standard().section(),
            product(2, 3, 1, 6),
            simple_packing(12, reference, e, d, 12),
            no_bitmap(),
            data(&pack_msb_first(values, 12)),
        ],
    )
}

/// A one-field message on a `ni` x `nj` grid of 8-bit packed values.
fn small_message(ni: u32, nj: u32, values: &[u64]) -> Vec<u8> {
    let grid = LatLonGrid {
        ni,
        nj,
        ..LatLonGrid::standard()
    };
    message(
        0,
        &[
            default_identification(),
            grid.section(),
            product(0, 0, 1, 0),
            simple_packing(values.len() as u32, 0.0, 0, 0, 8),
            no_bitmap(),
            data(&pack_msb_first(values, 8)),
        ],
    )
}

fn read_one(bytes: &[u8]) -> Result<GribMessage, String> {
    match GribReader::new(bytes).next_message() {
        Ok(Some(message)) => Ok(message),
        Ok(None) => Err("no message".to_string()),
        Err(error) => Err(error),
    }
}

/// The error of a read that must fail (the message type has no `Debug`
/// requirement here, so `unwrap_err` is avoided).
fn read_error(bytes: &[u8]) -> String {
    match read_one(bytes) {
        Ok(message) => panic!(
            "expected an error, but decoded {} values",
            message.data.len()
        ),
        Err(error) => error,
    }
}

/// Drain `reader`, continuing past errors; returns (messages, errors).
/// Fails when the reader does not reach the end of the stream within
/// `max_calls` calls, i.e. when it stops making progress.
fn drain(reader: &mut GribReader, max_calls: usize) -> (usize, usize) {
    let (mut messages, mut errors) = (0, 0);
    let mut calls = 0;
    loop {
        assert!(
            calls < max_calls,
            "reader made no progress in {max_calls} calls"
        );
        calls += 1;
        match reader.next_message() {
            Ok(Some(_)) => messages += 1,
            Ok(None) => return (messages, errors),
            Err(_) => errors += 1,
        }
    }
}

fn close(a: f64, b: f64) -> bool {
    (a - b).abs() < 1e-9
}

// ── F2384: message length guard and forward progress ─────────────────────

/// "GRIB", reserved, discipline 0, edition 2 and a total length field.
fn bare_indicator(total_len: u64) -> Vec<u8> {
    let mut buf = b"GRIB".to_vec();
    buf.extend_from_slice(&[0, 0, 0, 2]);
    buf.extend_from_slice(&total_len.to_be_bytes());
    buf
}

#[test]
fn total_length_zero_terminates() {
    let buf = bare_indicator(0);
    let mut reader = GribReader::new(&buf);
    // The old loop spun forever on `while let Some(m) = r.next_message()`.
    let before = reader.pos;
    let error = reader.next_message().unwrap_err();
    assert!(error.contains("length 0"), "{error}");
    assert!(
        reader.pos > before,
        "position must advance past the bad header"
    );
    assert!(matches!(reader.next_message(), Ok(None)));
    let mut reader = GribReader::new(&buf);
    assert_eq!(drain(&mut reader, 10), (0, 1));
}

#[test]
fn total_length_u64_max_does_not_panic_or_wrap() {
    let buf = bare_indicator(u64::MAX);
    let mut reader = GribReader::new(&buf);
    let error = reader.next_message().unwrap_err();
    assert!(error.contains("declares"), "{error}");
    assert!(matches!(reader.next_message(), Ok(None)));
    assert_eq!(drain(&mut GribReader::new(&buf), 10), (0, 1));
}

#[test]
fn total_length_overflowing_the_position_is_an_error() {
    // 8 bytes of garbage before the message: `pos + total_len` overflows a
    // usize for a length of u64::MAX - 4 once pos is non-zero.
    let mut buf = vec![0u8; 8];
    buf.extend_from_slice(&bare_indicator(u64::MAX - 4));
    assert_eq!(drain(&mut GribReader::new(&buf), 10), (0, 1));
}

#[test]
fn total_length_below_the_minimum_is_an_error() {
    // 20 is the smallest possible message: 16-byte indicator plus "7777".
    for total_len in [0u64, 1, 15, 16, 19] {
        let mut buf = bare_indicator(total_len);
        buf.extend_from_slice(b"7777");
        let mut reader = GribReader::new(&buf);
        let error = reader.next_message().unwrap_err();
        assert!(
            error.contains(&format!("length {total_len}")),
            "length {total_len}: {error}"
        );
    }
}

#[test]
fn truncated_message_is_an_error_then_the_end() {
    let full = standard_message(&[0; 12], 0.0, 0, 0);
    let truncated = &full[..full.len() - 10];
    let mut reader = GribReader::new(truncated);
    let error = reader.next_message().unwrap_err();
    assert!(error.contains("truncated"), "{error}");
    assert!(matches!(reader.next_message(), Ok(None)));
    // Cut inside the 16-byte indicator.
    let mut reader = GribReader::new(&full[..10]);
    assert!(reader.next_message().is_err());
    assert!(matches!(reader.next_message(), Ok(None)));
}

#[test]
fn empty_and_magic_free_streams_have_no_messages() {
    assert!(matches!(GribReader::new(&[]).next_message(), Ok(None)));
    assert!(matches!(
        GribReader::new(b"no grib message in here..").next_message(),
        Ok(None)
    ));
    // Trailing padding after a message is not an error.
    let mut buf = standard_message(&[0; 12], 0.0, 0, 0);
    buf.extend_from_slice(&[0u8; 37]);
    assert_eq!(drain(&mut GribReader::new(&buf), 10), (1, 0));
}

#[test]
fn corrupt_message_is_skipped_and_the_next_one_is_read() {
    let good = standard_message(&[5; 12], 0.0, 0, 0);
    let mut stream = bare_indicator(0);
    stream.extend_from_slice(&good);
    let mut reader = GribReader::new(&stream);
    assert!(reader.next_message().is_err());
    let second = reader.next_message().unwrap().unwrap();
    assert_eq!(second.data.len(), 12);
    assert!(matches!(reader.next_message(), Ok(None)));
}

#[test]
fn unsupported_edition_is_an_error_and_scanning_resumes() {
    let good = standard_message(&[5; 12], 0.0, 0, 0);
    let mut edition1 = bare_indicator(100);
    edition1[7] = 1;
    let mut stream = edition1;
    stream.extend_from_slice(&good);
    let mut reader = GribReader::new(&stream);
    let error = reader.next_message().unwrap_err();
    assert!(error.contains("edition 1"), "{error}");
    assert!(reader.next_message().unwrap().is_some());
    assert!(matches!(reader.next_message(), Ok(None)));
}

#[test]
fn message_without_the_end_marker_is_skipped() {
    let good = standard_message(&[5; 12], 0.0, 0, 0);
    let mut broken = good.clone();
    let len = broken.len();
    broken[len - 4..].copy_from_slice(b"7776");
    let mut stream = broken.clone();
    stream.extend_from_slice(&good);
    let mut reader = GribReader::new(&stream);
    let error = reader.next_message().unwrap_err();
    assert!(error.contains("7777"), "{error}");
    // The reader jumped to the end of the bad message, not into its middle.
    assert_eq!(reader.pos, broken.len());
    assert!(reader.next_message().unwrap().is_some());
}

#[test]
fn undecodable_message_is_skipped_and_the_next_one_is_read() {
    let complex = message(
        0,
        &[
            default_identification(),
            LatLonGrid::standard().section(),
            product(0, 0, 1, 0),
            other_packing(12, 2),
            no_bitmap(),
            data(&[0; 8]),
        ],
    );
    let good = standard_message(&[1; 12], 0.0, 0, 0);
    let mut stream = complex.clone();
    stream.extend_from_slice(&good);
    let mut reader = GribReader::new(&stream);
    let error = reader.next_message().unwrap_err();
    assert!(error.contains("5.2"), "{error}");
    assert_eq!(reader.pos, complex.len());
    assert!(reader.next_message().unwrap().is_some());
    assert!(matches!(reader.next_message(), Ok(None)));
}

#[test]
fn reader_terminates_on_every_prefix_and_single_byte_corruption() {
    let good = standard_message(&[0, 1, 100, 2048, 4095, 77, 5, 6, 7, 8, 9, 10], -5.5, -1, 1);
    for cut in 0..=good.len() {
        let mut reader = GribReader::new(&good[..cut]);
        drain(&mut reader, good.len() + 4);
    }
    for index in 0..good.len() {
        for value in [0x00u8, 0x01, 0x7F, 0x80, 0xFF] {
            let mut corrupted = good.clone();
            corrupted[index] = value;
            let mut reader = GribReader::new(&corrupted);
            drain(&mut reader, corrupted.len() + 4);
        }
    }
}

#[test]
fn every_call_on_single_field_messages_advances_the_position() {
    // Single-field messages and bad headers: each call consumes input. (The
    // queued later fields of a multi-field message are covered by
    // `queued_fields_leave_the_position_behind_the_message`.)
    let good = standard_message(&[1; 12], 0.0, 0, 0);
    let mut stream = bare_indicator(0);
    stream.extend_from_slice(&good);
    stream.extend_from_slice(&bare_indicator(u64::MAX));
    stream.extend_from_slice(&good);
    let mut reader = GribReader::new(&stream);
    loop {
        let before = reader.pos;
        match reader.next_message() {
            Ok(None) => break,
            Ok(Some(_)) | Err(_) => assert!(reader.pos > before, "stuck at {before}"),
        }
    }
}

// ── section walk ─────────────────────────────────────────────────────────

#[test]
fn section_list_of_a_complete_message() {
    let bytes = standard_message(&[0; 12], 0.0, 0, 0);
    let mut reader = GribReader::new(&bytes);
    assert!(reader.next_message().unwrap().is_some());
    let types: Vec<GribSectionType> = reader.sections().iter().map(|s| s.section_type).collect();
    assert_eq!(
        types,
        vec![
            GribSectionType::Indicator,
            GribSectionType::Identification,
            GribSectionType::GridDefinitionSection,
            GribSectionType::ProductDefinition,
            GribSectionType::DataRepresentationSection,
            GribSectionType::Bitmap,
            GribSectionType::Data,
            GribSectionType::End,
        ]
    );
    assert_eq!(reader.sections()[2].length, 72);
    assert_eq!(reader.sections()[2].payload.len(), 67);
}

#[test]
fn bare_indicator_and_end_marker_hold_no_field() {
    let bytes = message(0, &[]);
    assert_eq!(bytes.len(), 20);
    let mut reader = GribReader::new(&bytes);
    let error = reader.next_message().unwrap_err();
    assert!(error.contains("no data section"), "{error}");
    let types: Vec<GribSectionType> = reader.sections().iter().map(|s| s.section_type).collect();
    assert_eq!(
        types,
        vec![GribSectionType::Indicator, GribSectionType::End]
    );
    assert!(matches!(reader.next_message(), Ok(None)));
}

#[test]
fn malformed_section_lengths_are_errors() {
    let good = standard_message(&[0; 12], 0.0, 0, 0);
    // The grid section starts after the 16-byte indicator and the 21-byte
    // identification section.
    let grid_len_offset = 16 + 21;
    for bad_len in [0u32, 4, 73, 5000, u32::MAX] {
        let mut broken = good.clone();
        broken[grid_len_offset..grid_len_offset + 4].copy_from_slice(&bad_len.to_be_bytes());
        let error = read_error(&broken);
        assert!(error.contains("section"), "length {bad_len}: {error}");
    }
    // An unknown section number is not silently treated as "local use".
    let mut broken = good.clone();
    broken[grid_len_offset + 4] = 9;
    let error = read_error(&broken);
    assert!(error.contains("section number 9"), "{error}");
}

#[test]
fn sections_missing_before_the_data_are_errors() {
    let identification = default_identification();
    let grid = LatLonGrid::standard().section();
    let payload = data(&pack_msb_first(&[0; 12], 12));
    let cases: [(&str, Vec<Vec<u8>>); 5] = [
        (
            "identification",
            vec![
                grid.clone(),
                product(0, 0, 1, 0),
                simple_packing(12, 0.0, 0, 0, 12),
                no_bitmap(),
                payload.clone(),
            ],
        ),
        (
            "grid definition",
            vec![
                identification.clone(),
                product(0, 0, 1, 0),
                simple_packing(12, 0.0, 0, 0, 12),
                no_bitmap(),
                payload.clone(),
            ],
        ),
        (
            "product definition",
            vec![
                identification.clone(),
                grid.clone(),
                simple_packing(12, 0.0, 0, 0, 12),
                no_bitmap(),
                payload.clone(),
            ],
        ),
        (
            "data representation",
            vec![
                identification.clone(),
                grid.clone(),
                product(0, 0, 1, 0),
                no_bitmap(),
                payload.clone(),
            ],
        ),
        (
            "bit-map",
            vec![
                identification,
                grid,
                product(0, 0, 1, 0),
                simple_packing(12, 0.0, 0, 0, 12),
                payload,
            ],
        ),
    ];
    for (missing, sections) in cases {
        let error = read_error(&message(0, &sections));
        assert!(
            error.contains(&format!("without a {missing}"))
                || error.contains(&format!("without an {missing}")),
            "{missing}: {error}"
        );
    }
}

// ── F2386: Section 1, 3 and 4 octet offsets ──────────────────────────────

#[test]
fn identification_fields_are_decoded() {
    let decoded = read_one(&small_message(2, 1, &[1, 2])).unwrap();
    assert_eq!(decoded.centre_id, 98);
    assert_eq!(decoded.discipline, 0);
    // 2026-10-19 12:30:15 UTC.
    assert_eq!(decoded.reference_time, 1_792_413_015);
}

#[test]
fn reference_time_handles_leap_years_and_century_rules() {
    // (date, time of day, seconds since the epoch)
    let cases = [
        ((1970, 1, 1), (0, 0, 0), 0),
        ((2000, 3, 1), (0, 0, 0), 951_868_800),
        ((2024, 2, 29), (0, 0, 0), 1_709_164_800),
        ((1999, 12, 31), (23, 59, 59), 946_684_799),
        ((2100, 3, 1), (0, 0, 0), 4_107_542_400),
    ];
    for (date, time, expected) in cases {
        let bytes = message(
            0,
            &[
                identification(7, date, time),
                LatLonGrid {
                    ni: 1,
                    nj: 1,
                    ..LatLonGrid::standard()
                }
                .section(),
                product(0, 0, 1, 0),
                simple_packing(1, 0.0, 0, 0, 8),
                no_bitmap(),
                data(&[0]),
            ],
        );
        let decoded = read_one(&bytes).unwrap();
        assert_eq!(decoded.reference_time, expected, "{date:?} {time:?}");
    }
}

#[test]
fn invalid_reference_times_are_errors() {
    for (date, time) in [
        ((2026u16, 13u8, 1u8), (0u8, 0u8, 0u8)),
        ((2026, 0, 1), (0, 0, 0)),
        ((2026, 2, 30), (0, 0, 0)),
        ((2023, 2, 29), (0, 0, 0)),
        ((2026, 4, 31), (0, 0, 0)),
        ((2026, 1, 1), (24, 0, 0)),
        ((2026, 1, 1), (0, 60, 0)),
    ] {
        let bytes = message(
            0,
            &[
                identification(7, date, time),
                LatLonGrid::standard().section(),
            ],
        );
        let error = read_error(&bytes);
        assert!(
            error.contains("invalid reference time"),
            "{date:?} {time:?}: {error}"
        );
    }
}

#[test]
fn lat_lon_grid_is_read_from_the_wmo_octets() {
    let decoded = read_one(&standard_message(&[0; 12], 0.0, 0, 0)).unwrap();
    match decoded.grid_definition {
        GridDefinition::LatLon {
            ni,
            nj,
            lat_first,
            lat_last,
            lon_first,
            lon_last,
            di,
            dj,
        } => {
            // Ni counts points along a parallel (octets 31-34), Nj along a
            // meridian (35-38); the old reader took them from octets 15-22.
            assert_eq!((ni, nj), (4, 3));
            assert!(close(lat_first, 20.5), "lat_first={lat_first}");
            assert!(close(lat_last, -4.5), "lat_last={lat_last}");
            assert!(close(lon_first, 340.0), "lon_first={lon_first}");
            assert!(close(lon_last, 347.5), "lon_last={lon_last}");
            assert!(close(di, 2.5), "di={di}");
            assert!(close(dj, 12.5), "dj={dj}");
        }
        other => panic!("expected a lat/lon grid, got {other:?}"),
    }
    assert_eq!(decoded.grid_point_count(), 12);
}

#[test]
fn product_parameter_is_read_from_the_wmo_octets() {
    let decoded = read_one(&standard_message(&[0; 12], 0.0, 0, 0)).unwrap();
    // Template 4.0, category 2 / number 3 (octets 10 and 11). The old reader
    // took the template number and decoded every message as category 0 / 0.
    assert_eq!(decoded.parameter_category, 2);
    assert_eq!(decoded.parameter_number, 3);
    assert_eq!(decoded.parameter_description(), "V-component of wind");
    // 6 hours (unit 1 at octet 18, value at octets 19-22).
    assert_eq!(decoded.forecast_time, 6 * 3600);
}

#[test]
fn forecast_time_units_follow_code_table_4_4() {
    let cases: [(u8, u32, i64); 8] = [
        (0, 90, 90 * 60),
        (1, 6, 6 * 3600),
        (2, 2, 2 * 86_400),
        (10, 2, 2 * 10_800),
        (11, 2, 2 * 21_600),
        (12, 1, 43_200),
        (13, 45, 45),
        // Missing unit: no forecast time.
        (255, 7, 0),
    ];
    for (unit, value, expected) in cases {
        let bytes = message(
            0,
            &[
                default_identification(),
                LatLonGrid {
                    ni: 1,
                    nj: 1,
                    ..LatLonGrid::standard()
                }
                .section(),
                product(0, 0, unit, value),
                simple_packing(1, 0.0, 0, 0, 8),
                no_bitmap(),
                data(&[0]),
            ],
        );
        let decoded = read_one(&bytes).unwrap();
        assert_eq!(decoded.forecast_time, expected, "unit {unit}");
    }
}

#[test]
fn calendar_forecast_units_are_errors() {
    for unit in [3u8, 4, 5, 6, 7, 8, 9, 14, 200] {
        let bytes = message(
            0,
            &[
                default_identification(),
                LatLonGrid::standard().section(),
                product(0, 0, unit, 1),
                simple_packing(12, 0.0, 0, 0, 12),
                no_bitmap(),
                data(&[0; 18]),
            ],
        );
        let error = read_error(&bytes);
        assert!(error.contains("forecast time unit"), "unit {unit}: {error}");
    }
}

#[test]
fn product_templates_outside_the_common_layout_are_errors() {
    for template in [16u16, 20, 40, 1000] {
        let bytes = message(
            0,
            &[
                default_identification(),
                LatLonGrid::standard().section(),
                product_template(template),
                simple_packing(12, 0.0, 0, 0, 12),
                no_bitmap(),
                data(&[0; 18]),
            ],
        );
        let error = read_error(&bytes);
        assert!(
            error.contains(&format!("product definition template 4.{template}")),
            "{error}"
        );
    }
    // Template 4.8 (accumulation over a time interval) shares the layout.
    let bytes = message(
        0,
        &[
            default_identification(),
            LatLonGrid::standard().section(),
            product_template(8),
            simple_packing(12, 0.0, 0, 0, 12),
            no_bitmap(),
            data(&[0; 18]),
        ],
    );
    assert!(read_one(&bytes).is_ok());
}

#[test]
fn angle_unit_comes_from_the_basic_angle_and_its_subdivisions() {
    // Basic angle 1 with 1000 subdivisions: angles are in millidegrees.
    let grid = LatLonGrid {
        la1: 20_500,
        lo1: 340_000,
        la2: -4_500,
        lo2: 347_500,
        di: 2_500,
        dj: 12_500,
        basic_angle: 1,
        subdivisions: 1000,
        ..LatLonGrid::standard()
    };
    let bytes = message(
        0,
        &[
            default_identification(),
            grid.section(),
            product(0, 0, 1, 0),
            simple_packing(12, 0.0, 0, 0, 8),
            no_bitmap(),
            data(&[0; 12]),
        ],
    );
    match read_one(&bytes).unwrap().grid_definition {
        GridDefinition::LatLon {
            lat_first,
            lat_last,
            lon_first,
            di,
            dj,
            ..
        } => {
            assert!(close(lat_first, 20.5));
            assert!(close(lat_last, -4.5));
            assert!(close(lon_first, 340.0));
            assert!(close(di, 2.5));
            assert!(close(dj, 12.5));
        }
        other => panic!("expected a lat/lon grid, got {other:?}"),
    }
}

#[test]
fn missing_direction_increments_are_nan_not_a_huge_negative_angle() {
    let decode = |flags: u8, di: u32, dj: u32| {
        let grid = LatLonGrid {
            flags,
            di,
            dj,
            ..LatLonGrid::standard()
        };
        let bytes = message(
            0,
            &[
                default_identification(),
                grid.section(),
                product(0, 0, 1, 0),
                simple_packing(12, 0.0, 0, 0, 8),
                no_bitmap(),
                data(&[0; 12]),
            ],
        );
        match read_one(&bytes).unwrap().grid_definition {
            GridDefinition::LatLon { di, dj, .. } => (di, dj),
            other => panic!("expected a lat/lon grid, got {other:?}"),
        }
    };
    // Flags clear: the increments are not given, whatever the octets hold.
    let (di, dj) = decode(0x00, 2_500_000, 12_500_000);
    assert!(di.is_nan() && dj.is_nan());
    // Flags set but all-ones (missing): -2147.483647 degrees as sign-magnitude.
    let (di, dj) = decode(0x30, u32::MAX, u32::MAX);
    assert!(di.is_nan() && dj.is_nan());
    // Only the i direction increment is given (bit 3 = 0x20).
    let (di, dj) = decode(0x20, 2_500_000, u32::MAX);
    assert!(close(di, 2.5));
    assert!(dj.is_nan());
}

#[test]
fn polar_stereographic_grid_is_read_from_the_wmo_octets() {
    let bytes = message(
        0,
        &[
            default_identification(),
            polar_stereographic_grid(3, 2),
            product(0, 0, 1, 0),
            simple_packing(6, 0.0, 0, 0, 8),
            no_bitmap(),
            data(&[0; 6]),
        ],
    );
    match read_one(&bytes).unwrap().grid_definition {
        GridDefinition::PolarStereographic {
            nx,
            ny,
            lat_first,
            lon_first,
            dx,
            dy,
        } => {
            assert_eq!((nx, ny), (3, 2));
            assert!(close(lat_first, 30.0), "lat_first={lat_first}");
            assert!(close(lon_first, 279.0), "lon_first={lon_first}");
            assert!(close(dx, 25_000.0), "dx={dx}");
            assert!(close(dy, 25_000.0), "dy={dy}");
        }
        other => panic!("expected a polar stereographic grid, got {other:?}"),
    }
}

#[test]
fn lambert_conformal_grid_is_read_from_the_wmo_octets() {
    // A constant field (0 bits per value) over the 1 905 141 HRRR points.
    let bytes = message(
        0,
        &[
            default_identification(),
            lambert_grid(),
            product(0, 0, 1, 0),
            simple_packing(1799 * 1059, 280.0, 0, 0, 0),
            no_bitmap(),
            data(&[]),
        ],
    );
    let decoded = read_one(&bytes).unwrap();
    match decoded.grid_definition {
        GridDefinition::LambertConformal {
            nx,
            ny,
            lat_ref,
            lon_ref,
            latin1,
            latin2,
            dx,
            dy,
        } => {
            assert_eq!((nx, ny), (1799, 1059));
            assert!(close(lat_ref, 38.5), "lat_ref={lat_ref}");
            assert!(close(lon_ref, 262.5), "lon_ref={lon_ref}");
            assert!(close(latin1, 38.5));
            assert!(close(latin2, 38.5));
            assert!(close(dx, 3000.0), "dx={dx}");
            assert!(close(dy, 3000.0), "dy={dy}");
        }
        other => panic!("expected a Lambert conformal grid, got {other:?}"),
    }
    assert_eq!(decoded.data.len(), 1_905_141);
    assert!(decoded.data.iter().all(|&v| v == 280.0));
}

#[test]
fn unknown_grid_templates_are_errors_not_a_fabricated_grid() {
    for template in [1u16, 10, 40, 65535] {
        let section = SectionBuilder::new(3, 72)
            .u32(7, 12)
            .u16(13, template)
            .build();
        let bytes = message(0, &[default_identification(), section]);
        let error = read_error(&bytes);
        assert!(
            error.contains(&format!("grid definition template 3.{template}")),
            "{error}"
        );
    }
}

#[test]
fn quasi_regular_and_inconsistent_grids_are_errors() {
    // A list of points per row (octet 11 != 0).
    let list = SectionBuilder::new(3, 72)
        .u32(7, 12)
        .u8(11, 2)
        .u16(13, 0)
        .build();
    let error = read_error(&message(0, &[default_identification(), list]));
    assert!(error.contains("quasi-regular"), "{error}");
    // Missing Ni / Nj.
    let missing = LatLonGrid {
        ni: u32::MAX,
        ..LatLonGrid::standard()
    }
    .section();
    let error = read_error(&message(0, &[default_identification(), missing]));
    assert!(error.contains("quasi-regular"), "{error}");
    // 12 declared points on a 4 x 3 grid is fine; 13 is not.
    let mut inconsistent = LatLonGrid::standard().section();
    inconsistent[6..10].copy_from_slice(&13u32.to_be_bytes());
    let error = read_error(&message(0, &[default_identification(), inconsistent]));
    assert!(error.contains("declares 13 data points"), "{error}");
    // A section too short for the template.
    let short = SectionBuilder::new(3, 40).u32(7, 12).u16(13, 0).build();
    let error = read_error(&message(0, &[default_identification(), short]));
    assert!(error.contains("too short"), "{error}");
}

// ── F2387: Section 5 / 6 / 7 ─────────────────────────────────────────────

#[test]
fn simple_packing_with_twelve_bits() {
    // R = -5.5, E = -1 and D = 1 are sign-magnitude: E = 0x8001, D = 0x0001.
    let values = [0, 1, 100, 2048, 4095, 77, 5, 6, 7, 8, 9, 10];
    let decoded = read_one(&standard_message(&values, -5.5, -1, 1)).unwrap();
    assert_eq!(
        decoded.data_representation,
        DataRepresentation::SimplePacking
    );
    assert!(decoded.bitmap.is_none());
    assert_eq!(decoded.data.len(), 12);
    // Y = (R + X * 2^E) / 10^D, first values worked out by hand.
    let by_hand = [-0.55, -0.5, 4.45, 101.85, 204.2, 3.3];
    for (got, expected) in decoded.data.iter().zip(by_hand) {
        assert!(close(*got, expected), "got {got} expected {expected}");
    }
    for (x, got) in values.iter().zip(&decoded.data) {
        let expected = (-5.5 + *x as f64 * 0.5) / 10.0;
        assert!(
            close(*got, expected),
            "x={x}: got {got} expected {expected}"
        );
    }
}

#[test]
fn simple_packing_with_a_negative_decimal_scale() {
    // D = -1 (sign-magnitude 0x8001): Y = (R + X * 2^E) * 10.
    let values = [0, 1, 100, 2048, 4095, 77, 5, 6, 7, 8, 9, 10];
    let decoded = read_one(&standard_message(&values, -5.5, -1, -1)).unwrap();
    let by_hand = [-55.0, -50.0, 445.0, 10185.0, 20420.0, 330.0];
    for (got, expected) in decoded.data.iter().zip(by_hand) {
        assert!(close(*got, expected), "got {got} expected {expected}");
    }
}

#[test]
fn simple_packing_with_a_positive_binary_scale() {
    // E = 2: Y = R + 4 X, D = 0.
    let decoded = read_one(&standard_message(
        &[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11],
        10.0,
        2,
        0,
    ))
    .unwrap();
    for (x, got) in decoded.data.iter().enumerate() {
        assert!(close(*got, 10.0 + 4.0 * x as f64), "x={x} got {got}");
    }
}

#[test]
fn simple_packing_at_every_byte_and_bit_width() {
    for width in [1usize, 3, 7, 8, 9, 13, 16, 17, 24, 31, 32, 33, 48, 53, 64] {
        let max = if width == 64 {
            u64::MAX
        } else {
            (1u64 << width) - 1
        };
        let values: Vec<u64> = [
            0,
            1,
            max,
            max / 2,
            max / 3 + 1,
            0,
            max,
            1,
            2,
            3,
            max - 1,
            max / 5,
        ]
        .iter()
        .map(|v| v & max)
        .collect();
        let bytes = message(
            0,
            &[
                default_identification(),
                LatLonGrid::standard().section(),
                product(0, 0, 1, 0),
                simple_packing(12, 0.0, 0, 0, width as u8),
                no_bitmap(),
                data(&pack_msb_first(&values, width)),
            ],
        );
        let decoded = read_one(&bytes).unwrap();
        for (x, got) in values.iter().zip(&decoded.data) {
            assert_eq!(*got, *x as f64, "width {width}");
        }
    }
}

#[test]
fn zero_bit_packing_is_a_constant_field() {
    // R = 273.15 as f32 and D = 2: every value is R / 100.
    let bytes = message(
        0,
        &[
            default_identification(),
            LatLonGrid::standard().section(),
            product(0, 0, 1, 0),
            simple_packing(12, 273.15, 0, 2, 0),
            no_bitmap(),
            data(&[]),
        ],
    );
    let decoded = read_one(&bytes).unwrap();
    assert_eq!(decoded.data.len(), 12);
    let expected = f64::from(273.15_f32) / 100.0;
    assert!(decoded.data.iter().all(|&v| close(v, expected)));
}

#[test]
fn bitmap_expands_the_data_with_nan() {
    // Six points, bit-map 1 0 1 1 0 1 + two padding bits set to 1: the
    // padding must not count as grid points.
    let bytes = message(
        0,
        &[
            default_identification(),
            LatLonGrid {
                ni: 3,
                nj: 2,
                ..LatLonGrid::standard()
            }
            .section(),
            product(0, 0, 1, 0),
            simple_packing(4, 0.0, 0, 0, 8),
            bitmap(0, &[0b1011_0111]),
            data(&[10, 20, 30, 40]),
        ],
    );
    let decoded = read_one(&bytes).unwrap();
    assert_eq!(decoded.data.len(), 6);
    assert_eq!(decoded.data[0], 10.0);
    assert!(decoded.data[1].is_nan());
    assert_eq!(decoded.data[2], 20.0);
    assert_eq!(decoded.data[3], 30.0);
    assert!(decoded.data[4].is_nan());
    assert_eq!(decoded.data[5], 40.0);
    assert_eq!(
        decoded.bitmap,
        Some(vec![true, false, true, true, false, true])
    );
    // data and bitmap line up, as `valid_data` relies on.
    assert_eq!(decoded.data.len(), decoded.grid_point_count());
    assert_eq!(decoded.valid_data(), vec![10.0, 20.0, 30.0, 40.0]);
}

#[test]
fn bitmap_must_match_the_grid_and_the_value_count() {
    let grid = LatLonGrid {
        ni: 3,
        nj: 2,
        ..LatLonGrid::standard()
    }
    .section();
    let build = |packed_values: u32, bits: &[u8], payload: &[u8]| {
        message(
            0,
            &[
                default_identification(),
                grid.clone(),
                product(0, 0, 1, 0),
                simple_packing(packed_values, 0.0, 0, 0, 8),
                bitmap(0, bits),
                data(payload),
            ],
        )
    };
    // The bit-map marks 4 points present but section 5 declares 5 values.
    let error = read_error(&build(5, &[0b1011_0100], &[1, 2, 3, 4, 5]));
    assert!(error.contains("marks 4 points present"), "{error}");
    // 6 grid points need at least one full bit-map byte; an empty one is short.
    let error = read_error(&build(0, &[], &[]));
    assert!(error.contains("bit-map holds 0 bits"), "{error}");
    // Section 5 declares 6 values without a bit-map for a 6-point grid: fine,
    // but 4 values with no bit-map is a mismatch.
    let no_map = message(
        0,
        &[
            default_identification(),
            grid.clone(),
            product(0, 0, 1, 0),
            simple_packing(4, 0.0, 0, 0, 8),
            no_bitmap(),
            data(&[1, 2, 3, 4]),
        ],
    );
    let error = read_error(&no_map);
    assert!(error.contains("no bit-map applies"), "{error}");
}

#[test]
fn bitmap_indicator_254_reuses_the_bitmap_of_the_message() {
    let grid = LatLonGrid {
        ni: 3,
        nj: 2,
        ..LatLonGrid::standard()
    }
    .section();
    let bytes = message(
        0,
        &[
            default_identification(),
            grid.clone(),
            product(0, 0, 1, 0),
            simple_packing(4, 0.0, 0, 0, 8),
            bitmap(0, &[0b1011_0100]),
            data(&[10, 20, 30, 40]),
            product(0, 1, 1, 0),
            simple_packing(4, 0.0, 0, 0, 8),
            bitmap(254, &[]),
            data(&[50, 60, 70, 80]),
        ],
    );
    let mut reader = GribReader::new(&bytes);
    let first = reader.next_message().unwrap().unwrap();
    let second = reader.next_message().unwrap().unwrap();
    assert!(matches!(reader.next_message(), Ok(None)));
    assert_eq!(first.bitmap, second.bitmap);
    assert_eq!(second.valid_data(), vec![50.0, 60.0, 70.0, 80.0]);
    assert!(second.data[1].is_nan() && second.data[4].is_nan());
    assert_eq!(second.parameter_number, 1);

    // Indicator 254 with no earlier bit-map is an error, not "no bit-map".
    let orphan = message(
        0,
        &[
            default_identification(),
            grid,
            product(0, 0, 1, 0),
            simple_packing(6, 0.0, 0, 0, 8),
            bitmap(254, &[]),
            data(&[1, 2, 3, 4, 5, 6]),
        ],
    );
    let error = read_error(&orphan);
    assert!(error.contains("254"), "{error}");
}

#[test]
fn predefined_bitmaps_are_errors() {
    let bytes = message(
        0,
        &[
            default_identification(),
            LatLonGrid::standard().section(),
            product(0, 0, 1, 0),
            simple_packing(12, 0.0, 0, 0, 8),
            bitmap(3, &[]),
            data(&[0; 12]),
        ],
    );
    let error = read_error(&bytes);
    assert!(error.contains("pre-defined bit-map 3"), "{error}");
}

#[test]
fn ieee_floating_point_data() {
    let values = [
        1.5f32, -2.25, 1.0e10, 0.0, 3.25, -0.5, 7.0, 8.0, 9.0, 10.5, 11.0, -12.0,
    ];
    let payload: Vec<u8> = values.iter().flat_map(|v| v.to_be_bytes()).collect();
    let bytes = message(
        0,
        &[
            default_identification(),
            LatLonGrid::standard().section(),
            product(0, 0, 1, 0),
            ieee_packing(12, 1),
            no_bitmap(),
            data(&payload),
        ],
    );
    let decoded = read_one(&bytes).unwrap();
    assert_eq!(decoded.data_representation, DataRepresentation::IeeeFloat);
    for (got, expected) in decoded.data.iter().zip(values) {
        assert_eq!(*got, f64::from(expected));
    }

    let doubles: Vec<f64> = (0..12).map(|i| f64::from(i) * 0.1 - 0.3).collect();
    let payload: Vec<u8> = doubles.iter().flat_map(|v| v.to_be_bytes()).collect();
    let bytes = message(
        0,
        &[
            default_identification(),
            LatLonGrid::standard().section(),
            product(0, 0, 1, 0),
            ieee_packing(12, 2),
            no_bitmap(),
            data(&payload),
        ],
    );
    let decoded = read_one(&bytes).unwrap();
    assert_eq!(decoded.data, doubles);
}

#[test]
fn ieee_errors() {
    let build = |precision: u8, payload: &[u8]| {
        message(
            0,
            &[
                default_identification(),
                LatLonGrid::standard().section(),
                product(0, 0, 1, 0),
                ieee_packing(12, precision),
                no_bitmap(),
                data(payload),
            ],
        )
    };
    assert!(read_error(&build(3, &[0; 192])).contains("128-bit"));
    assert!(read_error(&build(9, &[0; 48])).contains("unknown IEEE precision 9"));
    let error = read_error(&build(1, &[0; 47]));
    assert!(error.contains("need 48"), "{error}");
}

#[test]
fn unsupported_data_representation_templates_are_errors() {
    for (template, name) in [
        (2u16, "complex packing"),
        (3, "spatial differencing"),
        (40, "JPEG 2000"),
        (41, "PNG"),
        (200, "run length"),
    ] {
        let bytes = message(
            0,
            &[
                default_identification(),
                LatLonGrid::standard().section(),
                product(0, 0, 1, 0),
                other_packing(12, template),
                no_bitmap(),
                data(&[0; 24]),
            ],
        );
        let error = read_error(&bytes);
        assert!(error.contains(&format!("5.{template}")), "{error}");
        assert!(error.contains(name), "{error}");
    }
    let bytes = message(
        0,
        &[
            default_identification(),
            LatLonGrid::standard().section(),
            product(0, 0, 1, 0),
            other_packing(12, 99),
            no_bitmap(),
            data(&[0; 24]),
        ],
    );
    assert!(read_error(&bytes).contains("unknown data representation template 5.99"));
}

#[test]
fn data_section_shorter_than_the_packed_values_is_an_error() {
    let bytes = message(
        0,
        &[
            default_identification(),
            LatLonGrid::standard().section(),
            product(0, 0, 1, 0),
            simple_packing(12, 0.0, 0, 0, 12),
            no_bitmap(),
            data(&[0; 17]),
        ],
    );
    let error = read_error(&bytes);
    assert!(
        error.contains("12 values of 12 bits need 18 bytes"),
        "{error}"
    );
}

#[test]
fn more_than_64_bits_per_value_is_an_error() {
    let bytes = message(
        0,
        &[
            default_identification(),
            LatLonGrid::standard().section(),
            product(0, 0, 1, 0),
            simple_packing(12, 0.0, 0, 0, 65),
            no_bitmap(),
            data(&[0; 100]),
        ],
    );
    let error = read_error(&bytes);
    assert!(error.contains("65 bits per packed value"), "{error}");
}

// ── allocation bounds ────────────────────────────────────────────────────

#[test]
fn huge_declared_grids_are_rejected_before_allocating() {
    // 65535 x 65537 = 4 294 967 295 points with a constant field: there are
    // no data bytes to bound the allocation, so the point limit does.
    let grid = LatLonGrid {
        ni: 65_535,
        nj: 65_537,
        ..LatLonGrid::standard()
    }
    .section();
    let bytes = message(
        0,
        &[
            default_identification(),
            grid,
            product(0, 0, 1, 0),
            simple_packing(u32::MAX, 1.0, 0, 0, 0),
            no_bitmap(),
            data(&[]),
        ],
    );
    let error = read_error(&bytes);
    assert!(error.contains("exceeds the reader limit"), "{error}");
}

#[test]
fn the_point_limit_is_configurable() {
    let bytes = standard_message(&[0; 12], 0.0, 0, 0);
    let mut tight = GribReader::new(&bytes).with_max_points(11);
    let error = tight.next_message().unwrap_err();
    assert!(error.contains("limit of 11 points"), "{error}");
    let mut exact = GribReader::new(&bytes).with_max_points(12);
    assert_eq!(exact.next_message().unwrap().unwrap().data.len(), 12);
    assert_eq!(
        GribReader::DEFAULT_MAX_POINTS,
        1 << 27,
        "the default limit documents itself"
    );
}

#[test]
fn packed_values_beyond_the_data_section_do_not_allocate() {
    // A 1000 x 1000 grid of 16-bit values with a 10-byte data section.
    let grid = LatLonGrid {
        ni: 1000,
        nj: 1000,
        ..LatLonGrid::standard()
    }
    .section();
    let bytes = message(
        0,
        &[
            default_identification(),
            grid,
            product(0, 0, 1, 0),
            simple_packing(1_000_000, 0.0, 0, 0, 16),
            no_bitmap(),
            data(&[0; 10]),
        ],
    );
    let error = read_error(&bytes);
    assert!(error.contains("need 2000000 bytes"), "{error}");
}

// ── repeated sections ────────────────────────────────────────────────────

#[test]
fn a_message_with_repeated_sections_yields_one_message_per_field() {
    let grid = LatLonGrid {
        ni: 3,
        nj: 2,
        ..LatLonGrid::standard()
    }
    .section();
    let bytes = message(
        0,
        &[
            default_identification(),
            grid,
            product(0, 0, 1, 0),
            simple_packing(6, 0.0, 0, 0, 8),
            no_bitmap(),
            data(&[1, 2, 3, 4, 5, 6]),
            // Second field: sections 4-7 repeated, the grid is reused.
            product(0, 6, 1, 12),
            simple_packing(6, 100.0, 0, 0, 8),
            no_bitmap(),
            data(&[0, 0, 0, 7, 7, 7]),
        ],
    );
    let mut reader = GribReader::new(&bytes);
    let first = reader.next_message().unwrap().unwrap();
    let second = reader.next_message().unwrap().unwrap();
    assert!(matches!(reader.next_message(), Ok(None)));
    assert_eq!(first.data, vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0]);
    assert_eq!(first.parameter_number, 0);
    assert_eq!(first.forecast_time, 0);
    assert_eq!(second.data, vec![100.0, 100.0, 100.0, 107.0, 107.0, 107.0]);
    assert_eq!(second.parameter_number, 6);
    assert_eq!(second.parameter_description(), "Dew point temperature");
    assert_eq!(second.forecast_time, 12 * 3600);
    assert_eq!(first.grid_definition, second.grid_definition);
    assert_eq!(first.reference_time, second.reference_time);
    // All sections of the message stay visible.
    assert_eq!(reader.sections().len(), 1 + 1 + 1 + 4 + 4 + 1);
}

#[test]
fn queued_fields_leave_the_position_behind_the_message() {
    // The first call decodes every field of the message and moves the
    // position behind it. The later fields are handed out from the queue
    // without moving it, and a loop that stops at `Ok(None)` still ends.
    let grid = LatLonGrid {
        ni: 3,
        nj: 2,
        ..LatLonGrid::standard()
    }
    .section();
    let mut sections = vec![default_identification(), grid];
    for (parameter, values) in [(0u8, [1u64; 6]), (1, [2; 6]), (2, [3; 6])] {
        sections.push(product(0, parameter, 1, 0));
        sections.push(simple_packing(6, 0.0, 0, 0, 8));
        sections.push(no_bitmap());
        sections.push(data(&pack_msb_first(&values, 8)));
    }
    let bytes = message(0, &sections);
    let mut reader = GribReader::new(&bytes);
    let first = reader.next_message().unwrap().unwrap();
    assert_eq!(first.parameter_number, 0);
    assert_eq!(first.data, vec![1.0; 6]);
    assert_eq!(reader.pos, bytes.len());
    for (expected_parameter, expected_value) in [(1u8, 2.0), (2, 3.0)] {
        let queued = reader.next_message().unwrap().unwrap();
        assert_eq!(queued.parameter_number, expected_parameter);
        assert_eq!(queued.data, vec![expected_value; 6]);
        assert_eq!(
            reader.pos,
            bytes.len(),
            "a queued field must not move the position"
        );
    }
    assert!(matches!(reader.next_message(), Ok(None)));
    assert_eq!(drain(&mut GribReader::new(&bytes), 5), (3, 0));
}

#[test]
fn fields_decoded_before_a_bad_field_are_still_returned() {
    let grid = LatLonGrid {
        ni: 3,
        nj: 2,
        ..LatLonGrid::standard()
    }
    .section();
    let bytes = message(
        0,
        &[
            default_identification(),
            grid,
            product(0, 0, 1, 0),
            simple_packing(6, 0.0, 0, 0, 8),
            no_bitmap(),
            data(&[1, 2, 3, 4, 5, 6]),
            product(0, 1, 1, 0),
            other_packing(6, 40),
            no_bitmap(),
            data(&[0; 6]),
        ],
    );
    let mut reader = GribReader::new(&bytes);
    assert_eq!(reader.next_message().unwrap().unwrap().data.len(), 6);
    let error = reader.next_message().unwrap_err();
    assert!(error.contains("5.40"), "{error}");
    assert!(matches!(reader.next_message(), Ok(None)));
}

#[test]
fn two_messages_in_one_stream() {
    let mut stream = small_message(2, 2, &[1, 2, 3, 4]);
    stream.extend_from_slice(&small_message(3, 1, &[9, 8, 7]));
    let mut reader = GribReader::new(&stream);
    let first = reader.next_message().unwrap().unwrap();
    let second = reader.next_message().unwrap().unwrap();
    assert!(matches!(reader.next_message(), Ok(None)));
    assert_eq!(first.data, vec![1.0, 2.0, 3.0, 4.0]);
    assert_eq!(second.data, vec![9.0, 8.0, 7.0]);
}
