// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Tests of the DICOM reader: sequences and items of undefined length,
//! Implicit and Explicit VR, the long/short VR length forms, encapsulated and
//! deflated data, and the rejection of truncated or hostile input.

use super::dicom_parse::{MAX_SEQUENCE_DEPTH, PixelData, parse_dicom_limited};
use super::functions::build_test_dicom;
use super::types::{DicomReader, ValueRepresentation};

/// Small deterministic xorshift generator (the crate has no `rand` dependency).
struct XorShift(u64);

impl XorShift {
    fn next_u64(&mut self) -> u64 {
        let mut x = self.0;
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        self.0 = x;
        x
    }

    /// Uniform index in `0..n` (`n` must not be zero).
    fn below(&mut self, n: usize) -> usize {
        (self.next_u64() % n as u64) as usize
    }
}

const IMPLICIT_LE: &str = "1.2.840.10008.1.2";
const EXPLICIT_LE: &str = "1.2.840.10008.1.2.1";
const DEFLATED_LE: &str = "1.2.840.10008.1.2.1.99";
const ENCAPSULATED_UNCOMPRESSED_LE: &str = "1.2.840.10008.1.2.1.98";
const EXPLICIT_BE: &str = "1.2.840.10008.1.2.2";
const JPEG_BASELINE: &str = "1.2.840.10008.1.2.4.50";

/// VRs with the long header form (two reserved bytes and a 4-byte length).
const LONG_VRS: [&str; 13] = [
    "OB", "OD", "OF", "OL", "OV", "OW", "SQ", "SV", "UC", "UN", "UR", "UT", "UV",
];
/// VRs with the short header form (a 2-byte length).
const SHORT_VRS: [&str; 21] = [
    "AE", "AS", "AT", "CS", "DA", "DS", "DT", "FD", "FL", "IS", "LO", "LT", "PN", "SH", "SL", "SS",
    "ST", "TM", "UI", "UL", "US",
];

// ---------------------------------------------------------------------------
// builders (independent of the reader's own tables)
// ---------------------------------------------------------------------------

fn tag_bytes(group: u16, element: u16) -> Vec<u8> {
    let mut out = group.to_le_bytes().to_vec();
    out.extend_from_slice(&element.to_le_bytes());
    out
}

/// `value` padded to an even length (NUL for UIDs and binary data, space for text).
fn even(vr: &str, value: &[u8]) -> Vec<u8> {
    let mut out = value.to_vec();
    if out.len() % 2 == 1 {
        out.push(if matches!(vr, "UI" | "OB" | "UN") {
            0
        } else {
            b' '
        });
    }
    out
}

/// An Explicit VR element with a defined length.
fn explicit(group: u16, element: u16, vr: &str, value: &[u8]) -> Vec<u8> {
    let value = even(vr, value);
    let mut out = tag_bytes(group, element);
    out.extend_from_slice(vr.as_bytes());
    if LONG_VRS.contains(&vr) {
        out.extend_from_slice(&[0, 0]);
        out.extend_from_slice(&(value.len() as u32).to_le_bytes());
    } else {
        assert!(SHORT_VRS.contains(&vr), "unknown VR {vr}");
        out.extend_from_slice(&(value.len() as u16).to_le_bytes());
    }
    out.extend_from_slice(&value);
    out
}

/// An Implicit VR element with a defined length.
fn implicit(group: u16, element: u16, value: &[u8]) -> Vec<u8> {
    let value = even("LO", value);
    let mut out = tag_bytes(group, element);
    out.extend_from_slice(&(value.len() as u32).to_le_bytes());
    out.extend_from_slice(&value);
    out
}

/// An item of defined length.
fn item(content: &[u8]) -> Vec<u8> {
    let mut out = tag_bytes(0xFFFE, 0xE000);
    out.extend_from_slice(&(content.len() as u32).to_le_bytes());
    out.extend_from_slice(content);
    out
}

/// An item of undefined length, closed by an item delimitation item.
fn item_undefined(content: &[u8]) -> Vec<u8> {
    let mut out = tag_bytes(0xFFFE, 0xE000);
    out.extend_from_slice(&u32::MAX.to_le_bytes());
    out.extend_from_slice(content);
    out.extend_from_slice(&tag_bytes(0xFFFE, 0xE00D));
    out.extend_from_slice(&0u32.to_le_bytes());
    out
}

fn sequence_delimiter() -> Vec<u8> {
    let mut out = tag_bytes(0xFFFE, 0xE0DD);
    out.extend_from_slice(&0u32.to_le_bytes());
    out
}

/// An Explicit VR `SQ` of undefined length holding `items`.
fn explicit_sequence(group: u16, element: u16, items: &[Vec<u8>]) -> Vec<u8> {
    let mut out = tag_bytes(group, element);
    out.extend_from_slice(b"SQ");
    out.extend_from_slice(&[0, 0]);
    out.extend_from_slice(&u32::MAX.to_le_bytes());
    for it in items {
        out.extend_from_slice(it);
    }
    out.extend_from_slice(&sequence_delimiter());
    out
}

/// An Implicit VR sequence of undefined length holding `items`.
fn implicit_sequence(group: u16, element: u16, items: &[Vec<u8>]) -> Vec<u8> {
    let mut out = tag_bytes(group, element);
    out.extend_from_slice(&u32::MAX.to_le_bytes());
    for it in items {
        out.extend_from_slice(it);
    }
    out.extend_from_slice(&sequence_delimiter());
    out
}

/// Explicit VR encapsulated pixel data: an empty offset table and one item per fragment.
fn encapsulated_pixel_data(fragments: &[&[u8]]) -> Vec<u8> {
    encapsulated_pixel_data_with_table(&[], fragments)
}

/// Explicit VR encapsulated pixel data with the given Basic Offset Table bytes.
fn encapsulated_pixel_data_with_table(table: &[u8], fragments: &[&[u8]]) -> Vec<u8> {
    let mut out = tag_bytes(0x7FE0, 0x0010);
    out.extend_from_slice(b"OB");
    out.extend_from_slice(&[0, 0]);
    out.extend_from_slice(&u32::MAX.to_le_bytes());
    out.extend_from_slice(&item(table));
    for fragment in fragments {
        out.extend_from_slice(&item(fragment));
    }
    out.extend_from_slice(&sequence_delimiter());
    out
}

fn us(value: u16) -> [u8; 2] {
    value.to_le_bytes()
}

/// Preamble, `DICM`, the file meta information announcing `syntax` (if any),
/// and `data_set`.
fn file(syntax: Option<&str>, data_set: &[u8]) -> Vec<u8> {
    let mut out = vec![0u8; 128];
    out.extend_from_slice(b"DICM");
    if let Some(uid) = syntax {
        let uid_element = explicit(0x0002, 0x0010, "UI", uid.as_bytes());
        out.extend_from_slice(&explicit(
            0x0002,
            0x0000,
            "UL",
            &(uid_element.len() as u32).to_le_bytes(),
        ));
        out.extend_from_slice(&uid_element);
    }
    out.extend_from_slice(data_set);
    out
}

fn modality(value: &str) -> Vec<u8> {
    explicit(0x0008, 0x0060, "CS", value.as_bytes())
}

fn image_size(rows: u16, columns: u16) -> Vec<u8> {
    [
        explicit(0x0028, 0x0010, "US", &us(rows)),
        explicit(0x0028, 0x0011, "US", &us(columns)),
    ]
    .concat()
}

fn pixel_data(bytes: &[u8]) -> Vec<u8> {
    explicit(0x7FE0, 0x0010, "OW", bytes)
}

fn parse(bytes: &[u8]) -> Result<(DicomReader, super::types::DicomHeader), String> {
    let mut reader = DicomReader::new();
    let header = reader.parse(bytes)?;
    Ok((reader, header))
}

fn parse_ok(bytes: &[u8]) -> (DicomReader, super::types::DicomHeader) {
    parse(bytes).unwrap_or_else(|e| panic!("parse failed: {e}"))
}

fn parse_err(bytes: &[u8]) -> String {
    match parse(bytes) {
        Ok(_) => panic!("expected the DICOM parse to fail"),
        Err(message) => message,
    }
}

// ---------------------------------------------------------------------------
// sequences and items of undefined length
// ---------------------------------------------------------------------------

#[test]
fn sequence_of_undefined_length_before_pixel_data_still_yields_pixel_data() {
    let referenced = item(&explicit(
        0x0008,
        0x1150,
        "UI",
        b"1.2.840.10008.5.1.4.1.1.2",
    ));
    let sequence = explicit_sequence(0x0008, 0x1140, std::slice::from_ref(&referenced));
    let data_set = [
        modality("CT"),
        sequence,
        image_size(4, 2),
        explicit(0x0028, 0x0100, "US", &us(8)),
        pixel_data(&[1, 2, 3, 4, 5, 6, 7, 8]),
    ]
    .concat();
    let (reader, header) = parse_ok(&file(None, &data_set));
    assert_eq!(header.modality, "CT");
    assert_eq!(
        (header.rows, header.columns, header.bits_allocated),
        (4, 2, 8)
    );
    assert_eq!(reader.pixel_data, [1, 2, 3, 4, 5, 6, 7, 8]);
    // Modality, the sequence, Rows, Columns and Bits Allocated; not the pixel data.
    assert_eq!(reader.element_count(), 5);
    let sequence_element = reader.find_element(0x0008, 0x1140).unwrap();
    assert_eq!(sequence_element.tag.vr, ValueRepresentation::SQ);
    assert_eq!(
        sequence_element.value, referenced,
        "the raw items, without the delimiter"
    );
}

#[test]
fn undefined_length_items_and_nested_sequences_are_walked() {
    let innermost = explicit_sequence(
        0x0040,
        0xA730,
        &[item_undefined(&explicit(0x0040, 0xA043, "SH", b"AB"))],
    );
    let first = item_undefined(&[explicit(0x0008, 0x1150, "UI", b"1.2.3"), innermost].concat());
    let second = item(&explicit(0x0008, 0x1155, "UI", b"1.2.4"));
    let sequence = explicit_sequence(0x0008, 0x1115, &[first.clone(), second.clone()]);
    let data_set = [
        sequence,
        image_size(3, 5),
        pixel_data(&[9, 8, 7, 6]),
        modality("MR"),
    ]
    .concat();
    let (reader, header) = parse_ok(&file(None, &data_set));
    assert_eq!((header.rows, header.columns), (3, 5));
    assert_eq!(
        header.modality, "MR",
        "an element after the pixel data is read"
    );
    assert_eq!(reader.pixel_data, [9, 8, 7, 6]);
    let sequence_element = reader.find_element(0x0008, 0x1115).unwrap();
    assert_eq!(sequence_element.value, [first, second].concat());
    // Elements inside the sequence are not top-level elements.
    assert!(reader.find_element(0x0040, 0xA043).is_none());
    assert!(reader.find_element(0x0008, 0x1150).is_none());
}

#[test]
fn undefined_length_un_holds_an_implicit_vr_sequence() {
    // PS3.5 6.2.2: the content of an undefined-length UN is Implicit VR Little
    // Endian even when the data set is Explicit VR.
    let inner = implicit_sequence(0x0008, 0x1140, &[item(&implicit(0x0008, 0x1150, b"1.2"))]);
    let un_items = [item_undefined(
        &[implicit(0x0008, 0x0060, b"CT"), inner].concat(),
    )];
    let mut un = tag_bytes(0x0009, 0x1010);
    un.extend_from_slice(b"UN");
    un.extend_from_slice(&[0, 0]);
    un.extend_from_slice(&u32::MAX.to_le_bytes());
    un.extend_from_slice(&un_items.concat());
    un.extend_from_slice(&sequence_delimiter());
    let data_set = [un, image_size(7, 9), pixel_data(&[1, 1])].concat();
    let (reader, header) = parse_ok(&file(None, &data_set));
    assert_eq!((header.rows, header.columns), (7, 9));
    assert_eq!(reader.pixel_data, [1, 1]);
    assert_eq!(
        reader.find_element(0x0009, 0x1010).unwrap().tag.vr,
        ValueRepresentation::SQ
    );
}

#[test]
fn sequence_of_defined_length_keeps_its_raw_value() {
    let content = [
        item_undefined(&explicit(0x0008, 0x1150, "UI", b"1.2.3")),
        item(&explicit(0x0008, 0x1155, "UI", b"1.2.4")),
    ]
    .concat();
    let mut sequence = tag_bytes(0x0008, 0x1140);
    sequence.extend_from_slice(b"SQ");
    sequence.extend_from_slice(&[0, 0]);
    sequence.extend_from_slice(&(content.len() as u32).to_le_bytes());
    sequence.extend_from_slice(&content);
    let data_set = [sequence, image_size(2, 2), pixel_data(&[5, 5, 5, 5])].concat();
    let (reader, header) = parse_ok(&file(None, &data_set));
    assert_eq!(header.rows, 2);
    assert_eq!(reader.find_element(0x0008, 0x1140).unwrap().value, content);
    assert_eq!(reader.pixel_data, [5, 5, 5, 5]);
}

#[test]
fn malformed_sequences_are_errors() {
    // No sequence delimitation item.
    let mut open = explicit_sequence(
        0x0008,
        0x1140,
        &[item(&explicit(0x0008, 0x1150, "UI", b"1"))],
    );
    open.truncate(open.len() - 8);
    let message = parse_err(&file(None, &open));
    assert!(message.contains("sequence delimitation"), "{message}");
    // An undefined-length item without its delimitation item.
    let mut item_open = tag_bytes(0xFFFE, 0xE000);
    item_open.extend_from_slice(&u32::MAX.to_le_bytes());
    item_open.extend_from_slice(&explicit(0x0008, 0x1150, "UI", b"1.2"));
    let unclosed = explicit_sequence(0x0008, 0x1140, &[item_open]);
    parse_err(&file(None, &unclosed));
    // A tag other than an item inside a sequence.
    let mut stray = tag_bytes(0x0008, 0x1140);
    stray.extend_from_slice(b"SQ");
    stray.extend_from_slice(&[0, 0]);
    stray.extend_from_slice(&u32::MAX.to_le_bytes());
    stray.extend_from_slice(&modality("CT"));
    stray.extend_from_slice(&sequence_delimiter());
    let message = parse_err(&file(None, &stray));
    assert!(message.contains("expected an item tag"), "{message}");
    // A delimiter with a length.
    let mut sized = tag_bytes(0xFFFE, 0xE0DD);
    sized.extend_from_slice(&5u32.to_le_bytes());
    let mut bad_delimiter = tag_bytes(0x0008, 0x1140);
    bad_delimiter.extend_from_slice(b"SQ");
    bad_delimiter.extend_from_slice(&[0, 0]);
    bad_delimiter.extend_from_slice(&u32::MAX.to_le_bytes());
    bad_delimiter.extend_from_slice(&sized);
    let message = parse_err(&file(None, &bad_delimiter));
    assert!(message.contains("non-zero length"), "{message}");
    // A new item header inside an undefined-length item.
    let nested_item = item_undefined(&item(&[]));
    let message = parse_err(&file(
        None,
        &explicit_sequence(0x0008, 0x1140, &[nested_item]),
    ));
    assert!(message.contains("inside an item"), "{message}");
}

#[test]
fn undefined_length_is_only_allowed_for_sequences_and_pixel_data() {
    let mut other_data = tag_bytes(0x0009, 0x1010);
    other_data.extend_from_slice(b"OB");
    other_data.extend_from_slice(&[0, 0]);
    other_data.extend_from_slice(&u32::MAX.to_le_bytes());
    other_data.extend_from_slice(&sequence_delimiter());
    let message = parse_err(&file(None, &other_data));
    assert!(message.contains("undefined length"), "{message}");
}

#[test]
fn nesting_is_bounded_and_deep_input_does_not_overflow_the_stack() {
    let nested = |depth: usize| {
        let mut inner = modality("CT");
        for _ in 0..depth {
            inner = explicit_sequence(0x0008, 0x1140, &[item_undefined(&inner)]);
        }
        [inner, image_size(1, 1)].concat()
    };
    let (_, header) = parse_ok(&file(None, &nested(MAX_SEQUENCE_DEPTH)));
    assert_eq!(header.rows, 1);
    for depth in [MAX_SEQUENCE_DEPTH + 1, 2_000] {
        let message = parse_err(&file(None, &nested(depth)));
        assert!(message.contains("nested deeper"), "{depth}: {message}");
    }
}

// ---------------------------------------------------------------------------
// VR length forms
// ---------------------------------------------------------------------------

#[test]
fn lt_has_a_two_byte_length_so_the_elements_after_it_are_found() {
    let data_set = [
        explicit(0x0020, 0x4000, "LT", b"Image comments: lateral view"),
        image_size(4, 2),
    ]
    .concat();
    let (reader, header) = parse_ok(&file(None, &data_set));
    assert_eq!((header.rows, header.columns), (4, 2));
    let comments = reader.find_element(0x0020, 0x4000).unwrap();
    assert_eq!(comments.tag.vr, ValueRepresentation::LT);
    assert_eq!(
        comments.value_as_string().as_deref(),
        Some("Image comments: lateral view")
    );
}

#[test]
fn every_vr_is_read_with_its_own_header_form() {
    let mut data_set = Vec::new();
    for (i, vr) in LONG_VRS.iter().chain(&SHORT_VRS).enumerate() {
        // A private element per VR, with an even-length value, then a marker.
        let value: &[u8] = if *vr == "SQ" {
            &[]
        } else {
            &[1, 2, 3, 4, 5, 6, 7, 8]
        };
        data_set.extend_from_slice(&explicit(0x0011, 0x1000 + i as u16, vr, value));
    }
    data_set.extend_from_slice(&modality("PT"));
    let (reader, header) = parse_ok(&file(None, &data_set));
    assert_eq!(header.modality, "PT");
    assert_eq!(reader.element_count(), LONG_VRS.len() + SHORT_VRS.len() + 1);
    for (i, vr) in LONG_VRS.iter().chain(&SHORT_VRS).enumerate() {
        let element = reader.find_element(0x0011, 0x1000 + i as u16).unwrap();
        assert_eq!(element.tag.vr.as_str(), *vr);
        let expected: &[u8] = if *vr == "SQ" {
            &[]
        } else {
            &[1, 2, 3, 4, 5, 6, 7, 8]
        };
        assert_eq!(element.value, expected, "{vr}");
    }
}

#[test]
fn vr_table_matches_the_standard() {
    for vr in LONG_VRS.iter().chain(&SHORT_VRS) {
        let b = vr.as_bytes();
        let parsed = ValueRepresentation::from_bytes(b[0], b[1]).unwrap();
        assert_eq!(parsed.as_str(), *vr);
        assert_eq!(
            parsed.has_extended_length(),
            LONG_VRS.contains(vr),
            "long header form of {vr}"
        );
    }
    assert!(ValueRepresentation::from_bytes(b'Z', b'Z').is_none());
    assert!(ValueRepresentation::from_bytes(b'u', b'n').is_none());
}

#[test]
fn an_unknown_vr_of_two_letters_is_read_like_un() {
    let mut zz = tag_bytes(0x0011, 0x1000);
    zz.extend_from_slice(b"ZZ");
    zz.extend_from_slice(&[0, 0]);
    zz.extend_from_slice(&4u32.to_le_bytes());
    zz.extend_from_slice(&[9, 9, 9, 9]);
    let (reader, header) = parse_ok(&file(None, &[zz, modality("MG")].concat()));
    assert_eq!(header.modality, "MG");
    let element = reader.find_element(0x0011, 0x1000).unwrap();
    assert_eq!(element.tag.vr, ValueRepresentation::UN);
    assert_eq!(element.value, [9, 9, 9, 9]);
}

#[test]
fn invalid_vr_bytes_are_an_error_that_mentions_implicit_vr() {
    let mut bad = tag_bytes(0x0008, 0x0060);
    bad.extend_from_slice(&[0x01, 0x02, 0x02, 0x00]);
    bad.extend_from_slice(b"CT");
    let message = parse_err(&file(None, &bad));
    assert!(message.contains("invalid VR"), "{message}");
    // Implicit VR data without a transfer syntax element is not guessed at.
    let message = parse_err(&file(None, &implicit(0x0008, 0x0060, b"CT")));
    assert!(message.contains("Implicit VR"), "{message}");
}

// ---------------------------------------------------------------------------
// transfer syntaxes
// ---------------------------------------------------------------------------

#[test]
fn implicit_vr_little_endian_data_sets_are_read() {
    let sequence = implicit_sequence(0x0008, 0x1140, &[item(&implicit(0x0008, 0x1150, b"1.2.3"))]);
    let data_set = [
        implicit(0x0008, 0x0000, &20u32.to_le_bytes()),
        implicit(0x0008, 0x0060, b"MR"),
        sequence,
        implicit(0x0009, 0x1001, b"AB"),
        implicit(0x0010, 0x0010, b"Roe^Rick"),
        implicit(0x0020, 0x0013, b"7 "),
        implicit(0x0028, 0x0010, &us(6)),
        implicit(0x0028, 0x0011, &us(5)),
        implicit(0x0028, 0x0100, &us(16)),
        implicit(0x7FE0, 0x0010, &[0xAA; 6]),
    ]
    .concat();
    let (reader, header) = parse_ok(&file(Some(IMPLICIT_LE), &data_set));
    assert_eq!(header.modality, "MR");
    assert_eq!(header.patient_name, "Roe^Rick");
    assert_eq!(header.instance_number, 7);
    assert_eq!(
        (header.rows, header.columns, header.bits_allocated),
        (6, 5, 16)
    );
    assert_eq!(reader.pixel_data, [0xAA; 6]);
    let vr_of = |group, element| reader.find_element(group, element).unwrap().tag.vr;
    assert_eq!(vr_of(0x0008, 0x0000), ValueRepresentation::UL);
    assert_eq!(vr_of(0x0008, 0x0060), ValueRepresentation::CS);
    assert_eq!(vr_of(0x0008, 0x1140), ValueRepresentation::SQ);
    assert_eq!(vr_of(0x0009, 0x1001), ValueRepresentation::UN);
    assert_eq!(vr_of(0x0010, 0x0010), ValueRepresentation::PN);
    assert_eq!(vr_of(0x0028, 0x0010), ValueRepresentation::US);
    // The meta information is Explicit VR whatever the data set uses.
    assert_eq!(
        reader
            .find_element(0x0002, 0x0010)
            .unwrap()
            .value_as_string()
            .as_deref(),
        Some(IMPLICIT_LE)
    );
    assert_eq!(vr_of(0x0002, 0x0010), ValueRepresentation::UI);
}

#[test]
fn explicit_vr_little_endian_is_read_with_and_without_the_meta_information() {
    let data_set = [modality("CR"), image_size(10, 20), pixel_data(&[3, 4])].concat();
    for syntax in [None, Some(EXPLICIT_LE)] {
        let (reader, header) = parse_ok(&file(syntax, &data_set));
        assert_eq!(header.modality, "CR");
        assert_eq!((header.rows, header.columns), (10, 20));
        assert_eq!(reader.pixel_data, [3, 4]);
        assert_eq!(
            reader.find_element(0x0002, 0x0010).is_some(),
            syntax.is_some()
        );
    }
}

#[test]
fn the_older_test_builder_still_parses() {
    let bytes = build_test_dicom(&[
        (0x0008, 0x0060, "CS", b"CT"),
        (0x0028, 0x0010, "US", &512u16.to_le_bytes()),
    ]);
    let (_, header) = parse_ok(&bytes);
    assert_eq!((header.modality.as_str(), header.rows), ("CT", 512));
}

#[test]
fn big_endian_data_is_rejected() {
    let bytes = file(
        Some(EXPLICIT_BE),
        &[modality("CT"), image_size(2, 2)].concat(),
    );
    let message = parse_err(&bytes);
    assert!(message.contains("Big Endian"), "{message}");
    let mut reader = DicomReader::new();
    assert!(reader.parse_header(&bytes).is_err());
}

#[test]
fn compressed_pixel_data_is_rejected_by_parse_but_the_header_is_readable() {
    let data_set = [
        modality("CT"),
        image_size(8, 8),
        encapsulated_pixel_data(&[&[1, 2, 3, 4], &[5, 6]]),
        explicit(0xFFFC, 0xFFFC, "OB", &[0; 6]),
    ]
    .concat();
    let bytes = file(Some(JPEG_BASELINE), &data_set);
    let message = parse_err(&bytes);
    assert!(message.contains(JPEG_BASELINE), "{message}");
    assert!(message.contains("compressed"), "{message}");

    let mut reader = DicomReader::new();
    let header = reader.parse_header(&bytes).unwrap();
    assert_eq!((header.modality.as_str(), header.rows), ("CT", 8));
    assert!(
        reader.pixel_data.is_empty(),
        "compressed bytes are never raw samples"
    );
    assert_eq!(
        reader
            .find_element(0x0002, 0x0010)
            .unwrap()
            .value_as_string()
            .as_deref(),
        Some(JPEG_BASELINE)
    );
    // The walk over the fragments ended at the right place: the trailing
    // padding element after the pixel data was read.
    assert!(reader.find_element(0xFFFC, 0xFFFC).is_some());
}

#[test]
fn encapsulated_pixel_data_without_a_transfer_syntax_is_not_returned_as_samples() {
    let data_set = [
        modality("US"),
        image_size(2, 2),
        encapsulated_pixel_data(&[&[0xFF, 0xD8, 0xFF, 0xD9]]),
    ]
    .concat();
    let bytes = file(None, &data_set);
    let message = parse_err(&bytes);
    assert!(message.contains("encapsulated"), "{message}");
    let mut reader = DicomReader::new();
    let header = reader.parse_header(&bytes).unwrap();
    assert_eq!(header.modality, "US");
    assert!(reader.pixel_data.is_empty());
}

#[test]
fn encapsulated_uncompressed_fragments_are_joined_into_the_pixel_data() {
    // Two frames of 4 bytes, in three fragments; the offset table lists the
    // byte offsets of the frames (0 and 4) and is not pixel data.
    let table: Vec<u8> = [0u32, 4].iter().flat_map(|o| o.to_le_bytes()).collect();
    let data_set = [
        modality("OT"),
        image_size(2, 2),
        encapsulated_pixel_data_with_table(&table, &[&[1, 2, 3], &[4, 5, 6, 7], &[8]]),
        explicit(0xFFFC, 0xFFFC, "OB", &[0; 6]),
    ]
    .concat();
    let bytes = file(Some(ENCAPSULATED_UNCOMPRESSED_LE), &data_set);
    let (reader, header) = parse_ok(&bytes);
    assert_eq!((header.modality.as_str(), header.rows), ("OT", 2));
    assert_eq!(reader.pixel_data, [1, 2, 3, 4, 5, 6, 7, 8]);
    assert!(
        reader.find_element(0xFFFC, 0xFFFC).is_some(),
        "the walk over the fragments ended at the delimiter"
    );
    // The header-only read skips the samples as it skips every pixel data.
    let mut header_only = DicomReader::new();
    assert_eq!(header_only.parse_header(&bytes).unwrap().rows, 2);
    assert!(header_only.pixel_data.is_empty());
}

#[test]
fn encapsulated_uncompressed_without_fragments_has_no_pixel_data() {
    let data_set = [
        modality("OT"),
        encapsulated_pixel_data_with_table(&[], &[]),
        image_size(5, 6),
    ]
    .concat();
    let (reader, header) = parse_ok(&file(Some(ENCAPSULATED_UNCOMPRESSED_LE), &data_set));
    assert_eq!(header.columns, 6);
    assert!(reader.pixel_data.is_empty());
}

#[test]
fn the_same_fragments_under_a_compressed_syntax_are_not_samples() {
    let data_set = [
        modality("OT"),
        encapsulated_pixel_data(&[&[1, 2, 3, 4]]),
        image_size(2, 2),
    ]
    .concat();
    for syntax in [JPEG_BASELINE, EXPLICIT_LE] {
        let message = parse_err(&file(Some(syntax), &data_set));
        assert!(message.contains("encapsulated"), "{syntax}: {message}");
    }
}

#[test]
fn a_fragment_of_undefined_length_is_an_error() {
    let mut pixels = tag_bytes(0x7FE0, 0x0010);
    pixels.extend_from_slice(b"OB");
    pixels.extend_from_slice(&[0, 0]);
    pixels.extend_from_slice(&u32::MAX.to_le_bytes());
    pixels.extend_from_slice(&item(&[]));
    pixels.extend_from_slice(&item_undefined(&[1, 2, 3, 4]));
    pixels.extend_from_slice(&sequence_delimiter());
    let bytes = file(
        Some(ENCAPSULATED_UNCOMPRESSED_LE),
        &[modality("OT"), pixels].concat(),
    );
    let message = parse_err(&bytes);
    assert!(message.contains("undefined length"), "{message}");
    assert!(DicomReader::new().parse_header(&bytes).is_err());
}

#[test]
fn a_tag_other_than_an_item_in_encapsulated_pixel_data_is_an_error() {
    let mut pixels = tag_bytes(0x7FE0, 0x0010);
    pixels.extend_from_slice(b"OB");
    pixels.extend_from_slice(&[0, 0]);
    pixels.extend_from_slice(&u32::MAX.to_le_bytes());
    pixels.extend_from_slice(&item(&[]));
    // An item delimitation item has no place between the fragments.
    pixels.extend_from_slice(&tag_bytes(0xFFFE, 0xE00D));
    pixels.extend_from_slice(&0u32.to_le_bytes());
    pixels.extend_from_slice(&sequence_delimiter());
    let message = parse_err(&file(
        Some(ENCAPSULATED_UNCOMPRESSED_LE),
        &[modality("OT"), pixels].concat(),
    ));
    assert!(
        message.contains("inside encapsulated pixel data"),
        "{message}"
    );
    assert!(message.contains("E00D"), "{message}");
}

#[test]
fn encapsulated_pixel_data_cut_anywhere_is_an_error() {
    let head = file(
        Some(ENCAPSULATED_UNCOMPRESSED_LE),
        &[modality("OT"), image_size(2, 2)].concat(),
    );
    let mut bytes = head.clone();
    bytes.extend_from_slice(&encapsulated_pixel_data(&[&[1, 2], &[3, 4]]));
    assert_eq!(parse_ok(&bytes).0.pixel_data, [1, 2, 3, 4]);
    // A file that stops before the pixel data is fine; one that stops inside it is not.
    assert!(parse(&head).is_ok());
    for cut in head.len() + 1..bytes.len() {
        assert!(parse(&bytes[..cut]).is_err(), "cut at {cut}");
    }
}

#[test]
fn unknown_transfer_syntaxes_are_rejected_by_parse_only() {
    let bytes = file(
        Some("1.2.3.4.5"),
        &[modality("CT"), image_size(2, 2)].concat(),
    );
    let message = parse_err(&bytes);
    assert!(message.contains("1.2.3.4.5"), "{message}");
    let mut reader = DicomReader::new();
    assert_eq!(reader.parse_header(&bytes).unwrap().rows, 2);
}

#[test]
fn deflated_data_sets_are_inflated() {
    let sequence = explicit_sequence(
        0x0008,
        0x1140,
        &[item_undefined(&explicit(0x0008, 0x1150, "UI", b"1.2.3"))],
    );
    let data_set = [
        modality("XA"),
        sequence,
        image_size(16, 12),
        pixel_data(&(0..=255u8).collect::<Vec<u8>>()),
    ]
    .concat();
    let deflated = oxiarc_deflate::deflate(&data_set, 6).unwrap();
    assert!(deflated.len() < data_set.len(), "the test data compresses");
    let bytes = file(Some(DEFLATED_LE), &deflated);
    let (reader, header) = parse_ok(&bytes);
    assert_eq!(header.modality, "XA");
    assert_eq!((header.rows, header.columns), (16, 12));
    assert_eq!(reader.pixel_data, (0..=255u8).collect::<Vec<u8>>());
    assert!(reader.find_element(0x0008, 0x1140).is_some());

    // A cut stream and a stream larger than the limit are errors.
    let cut = file(Some(DEFLATED_LE), &deflated[..deflated.len() / 2]);
    assert!(parse_err(&cut).contains("inflate"));
    let message = parse_dicom_limited(&bytes, PixelData::Keep, 16).unwrap_err();
    assert!(message.contains("inflate"), "{message}");
    // Bytes that are not a deflate stream.
    let garbage = file(Some(DEFLATED_LE), &[0xFF; 64]);
    assert!(parse_err(&garbage).contains("inflate"));
}

// ---------------------------------------------------------------------------
// truncated and hostile input
// ---------------------------------------------------------------------------

#[test]
fn a_value_longer_than_the_data_is_an_error_not_a_truncated_value() {
    let mut cut = tag_bytes(0x0008, 0x0060);
    cut.extend_from_slice(b"CS");
    cut.extend_from_slice(&100u16.to_le_bytes());
    cut.extend_from_slice(&[b'C'; 10]);
    let message = parse_err(&file(None, &cut));
    assert!(message.contains("needs 100 bytes"), "{message}");
    assert!(message.contains("(0008,0060)"), "{message}");
    // The same for a value whose 32-bit length is far beyond the file.
    let mut huge = tag_bytes(0x7FE0, 0x0010);
    huge.extend_from_slice(b"OW");
    huge.extend_from_slice(&[0, 0]);
    huge.extend_from_slice(&0xFFFF_FFF0u32.to_le_bytes());
    huge.extend_from_slice(&[0; 16]);
    let message = parse_err(&file(None, &huge));
    assert!(message.contains("needs"), "{message}");
}

#[test]
fn trailing_bytes_are_padding_or_an_error() {
    let data_set = [modality("CT"), image_size(2, 2)].concat();
    let valid = file(None, &data_set);
    for padding in [3usize, 7, 8, 20] {
        let mut padded = valid.clone();
        padded.extend(std::iter::repeat_n(0u8, padding));
        let (_, header) = parse_ok(&padded);
        assert_eq!(header.rows, 2, "{padding} zero bytes");
    }
    let mut stray = valid.clone();
    stray.extend_from_slice(&[1, 2, 3]);
    let message = parse_err(&stray);
    assert!(message.contains("stray"), "{message}");
    let mut garbage = valid;
    garbage.extend_from_slice(&[0xFF; 12]);
    assert!(parse(&garbage).is_err());
}

#[test]
fn item_and_delimitation_tags_outside_a_sequence_are_errors() {
    for element in [0xE000u16, 0xE00D, 0xE0DD] {
        let mut stray = modality("CT");
        stray.extend_from_slice(&tag_bytes(0xFFFE, element));
        stray.extend_from_slice(&0u32.to_le_bytes());
        let message = parse_err(&file(None, &stray));
        assert!(
            message.contains("outside a sequence"),
            "{element:04X}: {message}"
        );
    }
}

#[test]
fn every_prefix_of_a_file_is_an_error_or_an_exact_prefix_of_the_elements() {
    let data_set = [
        modality("CT"),
        explicit_sequence(
            0x0008,
            0x1140,
            &[
                item_undefined(&explicit(0x0008, 0x1150, "UI", b"1.2.3")),
                item(&explicit(0x0008, 0x1155, "UI", b"1.2.4")),
            ],
        ),
        explicit(0x0010, 0x0010, "PN", b"Doe^Jane"),
        image_size(4, 4),
        pixel_data(&[7; 32]),
    ]
    .concat();
    let bytes = file(Some(EXPLICIT_LE), &data_set);
    let (full, _) = parse_ok(&bytes);
    for cut in 132..bytes.len() {
        let mut reader = DicomReader::new();
        if reader.parse(&bytes[..cut]).is_ok() {
            // A cut on an element boundary: the elements read are the first ones, unchanged.
            assert!(reader.elements.len() < full.elements.len() || reader.pixel_data.is_empty());
            for (got, want) in reader.elements.iter().zip(&full.elements) {
                assert_eq!(got.tag, want.tag, "cut at {cut}");
                assert_eq!(got.value, want.value, "cut at {cut}");
            }
        }
    }
    for cut in 0..132 {
        assert!(DicomReader::new().parse(&bytes[..cut]).is_err());
    }
}

#[test]
fn a_failed_parse_leaves_the_reader_unchanged() {
    let good = file(
        None,
        &[modality("CT"), image_size(3, 3), pixel_data(&[1, 2])].concat(),
    );
    let mut reader = DicomReader::new();
    let header = reader.parse(&good).unwrap();
    assert_eq!(header.rows, 3);
    let before: Vec<(u16, u16, Vec<u8>)> = reader
        .elements
        .iter()
        .map(|e| (e.tag.group, e.tag.element, e.value.clone()))
        .collect();
    assert!(
        reader
            .parse(&file(Some(EXPLICIT_BE), &[modality("MR")].concat()))
            .is_err()
    );
    assert!(reader.parse(&good[..good.len() - 3]).is_err());
    let after: Vec<(u16, u16, Vec<u8>)> = reader
        .elements
        .iter()
        .map(|e| (e.tag.group, e.tag.element, e.value.clone()))
        .collect();
    assert_eq!(before, after);
    assert_eq!(reader.pixel_data, [1, 2]);
}

#[test]
fn pixel_data_of_length_zero_and_repeated_parses_start_fresh() {
    let first = file(None, &[modality("CT"), pixel_data(&[1, 2, 3, 4])].concat());
    let second = file(None, &[modality("MR"), pixel_data(&[])].concat());
    let mut reader = DicomReader::new();
    reader.parse(&first).unwrap();
    assert_eq!(reader.pixel_data, [1, 2, 3, 4]);
    let header = reader.parse(&second).unwrap();
    assert_eq!(header.modality, "MR");
    assert!(reader.pixel_data.is_empty());
    assert_eq!(reader.element_count(), 1);
}

/// One valid file of every kind of encoding the reader walks through.
fn corruption_corpus() -> Vec<Vec<u8>> {
    let explicit_set = [
        modality("CT"),
        explicit_sequence(
            0x0008,
            0x1140,
            &[
                item_undefined(&explicit(0x0008, 0x1150, "UI", b"1.2.3")),
                item(&explicit(0x0008, 0x1155, "UI", b"1.2.4")),
            ],
        ),
        explicit(0x0010, 0x0010, "PN", b"Doe^Jane"),
        image_size(4, 4),
        pixel_data(&[7; 32]),
    ]
    .concat();
    let implicit_set = [
        implicit(0x0008, 0x0060, b"MR"),
        implicit_sequence(
            0x0008,
            0x1140,
            &[item_undefined(&implicit(0x0008, 0x1150, b"1.2"))],
        ),
        implicit(0x0028, 0x0010, &us(8)),
        implicit(0x7FE0, 0x0010, &[1, 2, 3, 4]),
    ]
    .concat();
    let encapsulated_set = [
        modality("OT"),
        image_size(2, 2),
        encapsulated_pixel_data(&[&[1, 2], &[3, 4]]),
    ]
    .concat();
    vec![
        file(None, &explicit_set),
        file(Some(EXPLICIT_LE), &explicit_set),
        file(Some(IMPLICIT_LE), &implicit_set),
        file(Some(ENCAPSULATED_UNCOMPRESSED_LE), &encapsulated_set),
        file(Some(JPEG_BASELINE), &encapsulated_set),
    ]
}

#[test]
fn corrupted_files_never_panic_and_never_return_more_than_they_hold() {
    let mut rng = XorShift(0x9E37_79B9_7F4A_7C15);
    let (mut accepted, mut rejected) = (0_usize, 0_usize);
    for base in corruption_corpus() {
        for _ in 0..3000 {
            let mut bytes = base.clone();
            for _ in 0..1 + rng.below(4) {
                let at = 128 + rng.below(bytes.len() - 128);
                bytes[at] = rng.next_u64() as u8;
            }
            if rng.below(4) == 0 {
                let keep = 128 + rng.below(bytes.len() - 127);
                bytes.truncate(keep);
            }
            let mut reader = DicomReader::new();
            if reader.parse(&bytes).is_ok() {
                accepted += 1;
                assert!(reader.pixel_data.len() <= bytes.len());
                // Every element has a tag and a length: at least 8 bytes of the file.
                assert!(reader.elements.len() * 8 <= bytes.len());
                for element in &reader.elements {
                    assert!(element.value.len() <= bytes.len());
                }
            }
            let mut header_only = DicomReader::new();
            if header_only.parse_header(&bytes).is_ok() {
                assert!(header_only.pixel_data.is_empty());
            } else {
                rejected += 1;
            }
        }
    }
    // The corruption reaches both outcomes, so the checks above are not vacuous.
    assert!(accepted > 500, "{accepted} files were accepted");
    assert!(rejected > 500, "{rejected} headers were rejected");
}
