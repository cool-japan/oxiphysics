// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Parser for DICOM Part 10 files.
//!
//! A file is a 128-byte preamble, the magic `DICM`, the *file meta information*
//! (group `0002`, always Explicit VR Little Endian) and the data set, whose
//! encoding is announced by the Transfer Syntax UID `(0002,0010)`:
//!
//! | Transfer syntax | Handling |
//! |---|---|
//! | none announced | Explicit VR Little Endian |
//! | `1.2.840.10008.1.2` Implicit VR Little Endian | tag + `u32` length; VRs of known tags come from a small dictionary |
//! | `1.2.840.10008.1.2.1` Explicit VR Little Endian | tag + VR + length |
//! | `1.2.840.10008.1.2.1.99` Deflated Explicit VR Little Endian | the data set is inflated first |
//! | `1.2.840.10008.1.2.1.98` Encapsulated Uncompressed Explicit VR Little Endian | the samples are cut into fragments; they are joined into the pixel data |
//! | `1.2.840.10008.1.2.2` Explicit VR Big Endian | rejected |
//! | anything else (JPEG, JPEG 2000, RLE, ...) | rejected when pixel data is requested, because the pixel data is *encapsulated* (compressed) and is not raw samples |
//!
//! Sequences (`SQ`) and encapsulated pixel data have an *undefined length*
//! (`0xFFFFFFFF`) in many files.  Their extent is found by walking the items
//! (`FFFE,E000`) up to the sequence delimitation item (`FFFE,E0DD`); an item of
//! undefined length is walked element by element up to its item delimitation
//! item (`FFFE,E00D`), recursively.  The walk is bounded by
//! [`MAX_SEQUENCE_DEPTH`] nesting levels, so hostile input cannot exhaust the
//! stack.
//!
//! Encapsulated pixel data is a list of items: the first is the Basic Offset
//! Table, each of the others is one fragment of a frame, and the sequence
//! delimitation item ends the list.  The fragments are raw samples only in the
//! Encapsulated Uncompressed transfer syntax; with any other syntax they are
//! compressed frames, which this module does not decode.
//!
//! Every length is checked against the bytes that remain: a value that runs past
//! the end of the data, a sequence without its delimiter, and bytes at the end
//! that cannot form an element are errors, never silently truncated data.

use super::types::{DicomElement, DicomTag, ValueRepresentation};
use oxiarc_deflate::InflateStream;

/// Implicit VR Little Endian.
const IMPLICIT_VR_LITTLE_ENDIAN: &str = "1.2.840.10008.1.2";
/// Explicit VR Little Endian.
const EXPLICIT_VR_LITTLE_ENDIAN: &str = "1.2.840.10008.1.2.1";
/// Deflated Explicit VR Little Endian.
const DEFLATED_EXPLICIT_VR_LITTLE_ENDIAN: &str = "1.2.840.10008.1.2.1.99";
/// Encapsulated Uncompressed Explicit VR Little Endian.
const ENCAPSULATED_UNCOMPRESSED_EXPLICIT_VR_LITTLE_ENDIAN: &str = "1.2.840.10008.1.2.1.98";
/// Explicit VR Big Endian (retired).
const EXPLICIT_VR_BIG_ENDIAN: &str = "1.2.840.10008.1.2.2";

/// Group of the file meta information elements.
const FILE_META_GROUP: u16 = 0x0002;
/// Transfer Syntax UID.
const TRANSFER_SYNTAX_UID: (u16, u16) = (0x0002, 0x0010);
/// Pixel Data.
const PIXEL_DATA: (u16, u16) = (0x7FE0, 0x0010);
/// Group of the item and delimitation tags.
const ITEM_GROUP: u16 = 0xFFFE;
/// Item (`FFFE,E000`).
const ITEM: u16 = 0xE000;
/// Item Delimitation Item (`FFFE,E00D`).
const ITEM_DELIMITATION: u16 = 0xE00D;
/// Sequence Delimitation Item (`FFFE,E0DD`).
const SEQUENCE_DELIMITATION: u16 = 0xE0DD;
/// The length value that marks an element of undefined length.
const UNDEFINED_LENGTH: u32 = u32::MAX;
/// Deepest nesting of undefined-length sequences that is walked.
pub(super) const MAX_SEQUENCE_DEPTH: usize = 32;
/// Largest data set the Deflated transfer syntax may inflate to.
const MAX_INFLATED_BYTES: u64 = 1 << 30;

/// What to do with the Pixel Data element.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum PixelData {
    /// Return native pixel samples; encapsulated (compressed) pixel data is an
    /// error.
    Keep,
    /// Skip the pixel data, whatever its encoding.
    Skip,
}

/// The elements of a parsed file.
#[derive(Debug, Default)]
pub(super) struct DicomDataSet {
    /// Every element in file order except Pixel Data: the file meta
    /// information first, then the data set.  Sequence elements hold their raw
    /// encoded items.
    pub(super) elements: Vec<DicomElement>,
    /// Native pixel samples; empty without pixel data or with
    /// [`PixelData::Skip`].
    pub(super) pixel_data: Vec<u8>,
}

/// How the elements of a data set are encoded (always little endian).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Encoding {
    /// Tag, `u32` length.
    Implicit,
    /// Tag, two VR bytes, then a `u16` length or (long VRs) two reserved bytes
    /// and a `u32` length.
    Explicit,
}

/// Parse a DICOM Part 10 file.
pub(super) fn parse_dicom(data: &[u8], pixels: PixelData) -> Result<DicomDataSet, String> {
    parse_dicom_limited(data, pixels, MAX_INFLATED_BYTES)
}

/// [`parse_dicom`] with a custom limit for the size of an inflated data set.
pub(super) fn parse_dicom_limited(
    data: &[u8],
    pixels: PixelData,
    max_inflated_bytes: u64,
) -> Result<DicomDataSet, String> {
    if data.len() < 132 {
        return Err("Data too short for DICOM preamble + magic".into());
    }
    if data.get(128..132) != Some(b"DICM".as_slice()) {
        return Err("Missing DICM magic number at offset 128".into());
    }
    let mut out = DicomDataSet::default();
    let mut cursor = Cursor::new(data, 132);
    // The file meta information is Explicit VR Little Endian whatever the data set uses.
    while cursor
        .peek_tag()
        .is_some_and(|(group, _)| group == FILE_META_GROUP)
    {
        read_element(
            &mut cursor,
            Encoding::Explicit,
            Fragments::Opaque,
            pixels,
            &mut out,
        )?;
    }
    let syntax = out
        .elements
        .iter()
        .find(|e| (e.tag.group, e.tag.element) == TRANSFER_SYNTAX_UID)
        .and_then(DicomElement::value_as_string)
        .filter(|uid| !uid.is_empty());
    match transfer_syntax_plan(syntax.as_deref(), pixels)? {
        Plan::Direct(encoding, fragments) => {
            parse_data_set(&mut cursor, encoding, fragments, pixels, &mut out)?;
        }
        Plan::Deflated => {
            let deflated = data.get(cursor.pos..).unwrap_or_default();
            let inflated = InflateStream::new()
                .with_max_output(max_inflated_bytes)
                .inflate_to_vec(deflated)
                .map_err(|e| format!("DICOM: cannot inflate the deflated data set: {e}"))?;
            let mut inner = Cursor::new(&inflated, 0);
            parse_data_set(
                &mut inner,
                Encoding::Explicit,
                Fragments::Opaque,
                pixels,
                &mut out,
            )?;
        }
    }
    Ok(out)
}

/// What the fragments of encapsulated (undefined-length) Pixel Data hold.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Fragments {
    /// Not raw samples: compressed frames, or a transfer syntax that does not
    /// encapsulate pixel data at all.
    Opaque,
    /// Raw samples cut into fragments (Encapsulated Uncompressed).
    Samples,
}

/// How the data set following the meta information is read.
enum Plan {
    /// Parse the bytes as they are, with this encoding.
    Direct(Encoding, Fragments),
    /// Inflate the bytes first; the result is Explicit VR Little Endian.
    Deflated,
}

fn transfer_syntax_plan(uid: Option<&str>, pixels: PixelData) -> Result<Plan, String> {
    match uid {
        None | Some(EXPLICIT_VR_LITTLE_ENDIAN) => {
            Ok(Plan::Direct(Encoding::Explicit, Fragments::Opaque))
        }
        Some(IMPLICIT_VR_LITTLE_ENDIAN) => Ok(Plan::Direct(Encoding::Implicit, Fragments::Opaque)),
        Some(DEFLATED_EXPLICIT_VR_LITTLE_ENDIAN) => Ok(Plan::Deflated),
        Some(ENCAPSULATED_UNCOMPRESSED_EXPLICIT_VR_LITTLE_ENDIAN) => {
            Ok(Plan::Direct(Encoding::Explicit, Fragments::Samples))
        }
        Some(EXPLICIT_VR_BIG_ENDIAN) => Err(format!(
            "DICOM: the Explicit VR Big Endian transfer syntax ({EXPLICIT_VR_BIG_ENDIAN}) is not \
             supported"
        )),
        // Every other syntax encodes the data set as Explicit VR Little Endian and
        // compresses the pixel data (encapsulated pixel data).
        Some(other) => match pixels {
            PixelData::Keep => Err(format!(
                "DICOM: unsupported transfer syntax {other}: its pixel data is compressed \
                 (encapsulated) and cannot be returned as raw samples; use parse_header to read \
                 the header of such a file"
            )),
            PixelData::Skip => Ok(Plan::Direct(Encoding::Explicit, Fragments::Opaque)),
        },
    }
}

/// Read elements until the data ends.
fn parse_data_set(
    cursor: &mut Cursor<'_>,
    encoding: Encoding,
    fragments: Fragments,
    pixels: PixelData,
    out: &mut DicomDataSet,
) -> Result<(), String> {
    while !at_end(cursor)? {
        read_element(cursor, encoding, fragments, pixels, out)?;
    }
    Ok(())
}

/// Whether only the end of the data remains.
///
/// Zero padding after the last element is accepted; other bytes too short to
/// hold an element are an error (the file is cut inside an element header).
fn at_end(cursor: &Cursor<'_>) -> Result<bool, String> {
    let rest = cursor.rest();
    if rest.is_empty() {
        return Ok(true);
    }
    if (rest.len() < 8 || rest.starts_with(&[0, 0, 0, 0])) && rest.iter().all(|&b| b == 0) {
        return Ok(true);
    }
    if rest.len() < 8 {
        return Err(format!(
            "DICOM: {} stray bytes at offset {} cannot hold a data element",
            rest.len(),
            cursor.pos
        ));
    }
    Ok(false)
}

/// Read one element (header and value) and record it.
///
/// `fragments` says what encapsulated Pixel Data holds; it matters only for a
/// Pixel Data element of undefined length.
fn read_element(
    cursor: &mut Cursor<'_>,
    encoding: Encoding,
    fragments: Fragments,
    pixels: PixelData,
    out: &mut DicomDataSet,
) -> Result<(), String> {
    let header = read_header(cursor, encoding)?;
    if header.group == ITEM_GROUP {
        return Err(format!(
            "DICOM: unexpected item or delimitation tag ({:04X},{:04X}) at offset {} outside a \
             sequence",
            header.group, header.element, header.offset
        ));
    }
    let tag = header.tag();
    if header.length == UNDEFINED_LENGTH {
        check_undefined_length(&header)?;
        if tag == PIXEL_DATA {
            let items = read_fragments(cursor)?;
            return match (pixels, fragments) {
                (PixelData::Skip, _) => Ok(()),
                (PixelData::Keep, Fragments::Samples) => {
                    // The first item is the Basic Offset Table, not pixel data.
                    out.pixel_data = items
                        .iter()
                        .skip(1)
                        .flat_map(|f| f.iter().copied())
                        .collect();
                    Ok(())
                }
                (PixelData::Keep, Fragments::Opaque) => Err(format!(
                    "DICOM: the pixel data at offset {} is encapsulated (compressed) and cannot \
                     be returned as raw samples; use parse_header to skip it",
                    header.offset
                )),
            };
        }
        let content = scan_sequence(cursor, header.nested_encoding(encoding), 1)?;
        out.elements.push(DicomElement::new(
            DicomTag::new(tag.0, tag.1, ValueRepresentation::SQ),
            content.to_vec(),
        ));
        return Ok(());
    }
    let value = cursor.take(header.length_bytes(), || {
        format!("the value of ({:04X},{:04X})", header.group, header.element)
    })?;
    if tag == PIXEL_DATA {
        if pixels == PixelData::Keep {
            out.pixel_data = value.to_vec();
        }
        return Ok(());
    }
    out.elements.push(DicomElement::new(
        DicomTag::new(tag.0, tag.1, header.value_representation()),
        value.to_vec(),
    ));
    Ok(())
}

/// Tag, VR and length of one element, as they appear in the stream.
#[derive(Clone, Copy, Debug)]
struct ElementHeader {
    group: u16,
    element: u16,
    /// The two VR bytes of an Explicit VR element; `None` in Implicit VR and for
    /// item and delimitation tags, which never carry a VR.
    vr: Option<[u8; 2]>,
    /// The raw length field.
    length: u32,
    /// Offset of the tag in the data, for messages.
    offset: usize,
}

impl ElementHeader {
    fn tag(&self) -> (u16, u16) {
        (self.group, self.element)
    }

    fn length_bytes(&self) -> usize {
        // A u32 length always fits `usize` on the 32/64-bit targets.
        self.length as usize
    }

    /// The VR from the stream, or from the dictionary in Implicit VR.  A VR code
    /// this reader does not know is `UN`.
    fn value_representation(&self) -> ValueRepresentation {
        match self.vr {
            Some(code) => {
                ValueRepresentation::from_bytes(code[0], code[1]).unwrap_or(ValueRepresentation::UN)
            }
            None => implicit_vr(self.group, self.element),
        }
    }

    /// The encoding of the items of this element when its length is undefined.
    ///
    /// The content of an undefined-length `UN` element is Implicit VR Little
    /// Endian whatever the surrounding encoding (PS3.5 6.2.2).
    fn nested_encoding(&self, outer: Encoding) -> Encoding {
        match self.vr {
            Some(code)
                if ValueRepresentation::from_bytes(code[0], code[1])
                    .is_none_or(|vr| vr == ValueRepresentation::UN) =>
            {
                Encoding::Implicit
            }
            _ => outer,
        }
    }
}

/// An undefined length is allowed for sequences (`SQ`, or `UN` holding one),
/// for encapsulated pixel data, and for every element of Implicit VR data.
fn check_undefined_length(header: &ElementHeader) -> Result<(), String> {
    let allowed = match header.vr {
        None => true,
        Some(code) => match ValueRepresentation::from_bytes(code[0], code[1]) {
            Some(ValueRepresentation::SQ | ValueRepresentation::UN) | None => true,
            Some(ValueRepresentation::OB | ValueRepresentation::OW) => header.tag() == PIXEL_DATA,
            Some(_) => false,
        },
    };
    if allowed {
        Ok(())
    } else {
        Err(format!(
            "DICOM: element ({:04X},{:04X}) at offset {} has an undefined length but is neither a \
             sequence nor encapsulated pixel data",
            header.group, header.element, header.offset
        ))
    }
}

/// Read the header of one element.
///
/// Item and delimitation tags (`FFFE,xxxx`) are a tag and a `u32` length in
/// every transfer syntax: they have no VR.
fn read_header(cursor: &mut Cursor<'_>, encoding: Encoding) -> Result<ElementHeader, String> {
    let offset = cursor.pos;
    let group = cursor.u16(|| "a tag".to_string())?;
    let element = cursor.u16(|| "a tag".to_string())?;
    if group == ITEM_GROUP || encoding == Encoding::Implicit {
        let length = cursor.u32(|| format!("the length of ({group:04X},{element:04X})"))?;
        return Ok(ElementHeader {
            group,
            element,
            vr: None,
            length,
            offset,
        });
    }
    let code: [u8; 2] = cursor.array(|| format!("the VR of ({group:04X},{element:04X})"))?;
    if !code.iter().all(u8::is_ascii_uppercase) {
        return Err(format!(
            "DICOM: invalid VR bytes {:#04x} {:#04x} for tag ({group:04X},{element:04X}) at \
             offset {offset}: the data set is corrupt, or it is Implicit VR without a transfer \
             syntax element",
            code[0], code[1]
        ));
    }
    // A VR code that is two letters but unknown is read like `UN`, which has the long form.
    let long_form =
        ValueRepresentation::from_bytes(code[0], code[1]).is_none_or(|vr| vr.has_extended_length());
    let length = if long_form {
        cursor.take(2, || "the reserved bytes of a long-form length".to_string())?;
        cursor.u32(|| format!("the length of ({group:04X},{element:04X})"))?
    } else {
        u32::from(cursor.u16(|| format!("the length of ({group:04X},{element:04X})"))?)
    };
    Ok(ElementHeader {
        group,
        element,
        vr: Some(code),
        length,
        offset,
    })
}

/// Walk the items of an undefined-length sequence (or of encapsulated pixel
/// data); the cursor is at the first item header.
///
/// Returns the bytes of the items, without the sequence delimitation item that
/// ends them, and leaves the cursor after that item.
fn scan_sequence<'a>(
    cursor: &mut Cursor<'a>,
    encoding: Encoding,
    depth: usize,
) -> Result<&'a [u8], String> {
    if depth > MAX_SEQUENCE_DEPTH {
        return Err(format!(
            "DICOM: sequences are nested deeper than {MAX_SEQUENCE_DEPTH} levels at offset {}",
            cursor.pos
        ));
    }
    let start = cursor.pos;
    loop {
        let item_offset = cursor.pos;
        let (element, length) = read_item_header(cursor)?;
        match element {
            SEQUENCE_DELIMITATION => {
                require_empty_delimiter(length, item_offset)?;
                return cursor
                    .slice(start, item_offset)
                    .ok_or_else(|| "DICOM: internal error: sequence range".to_string());
            }
            ITEM => {
                if length == UNDEFINED_LENGTH {
                    scan_item(cursor, encoding, depth)?;
                } else {
                    cursor.take(length as usize, || {
                        format!("the item at offset {item_offset}")
                    })?;
                }
            }
            other => {
                return Err(format!(
                    "DICOM: unexpected delimitation tag (FFFE,{other:04X}) at offset \
                     {item_offset} inside a sequence"
                ));
            }
        }
    }
}

/// Read the items of encapsulated pixel data: the Basic Offset Table, then one
/// item per fragment; the cursor is at the first item header.
///
/// Returns the value of every item in order (the offset table first) and leaves
/// the cursor after the sequence delimitation item.  A fragment always has a
/// defined length.
fn read_fragments<'a>(cursor: &mut Cursor<'a>) -> Result<Vec<&'a [u8]>, String> {
    let mut items = Vec::new();
    loop {
        let item_offset = cursor.pos;
        let (element, length) = read_item_header(cursor)?;
        match element {
            SEQUENCE_DELIMITATION => {
                require_empty_delimiter(length, item_offset)?;
                return Ok(items);
            }
            ITEM if length != UNDEFINED_LENGTH => {
                items.push(cursor.take(length as usize, || {
                    format!("the fragment at offset {item_offset}")
                })?);
            }
            ITEM => {
                return Err(format!(
                    "DICOM: the fragment at offset {item_offset} of the encapsulated pixel data \
                     has an undefined length"
                ));
            }
            other => {
                return Err(format!(
                    "DICOM: unexpected delimitation tag (FFFE,{other:04X}) at offset \
                     {item_offset} inside encapsulated pixel data"
                ));
            }
        }
    }
}

/// Read an item or delimitation tag and its length.
fn read_item_header(cursor: &mut Cursor<'_>) -> Result<(u16, u32), String> {
    let offset = cursor.pos;
    if cursor.rest().len() < 8 {
        return Err(format!(
            "DICOM: the sequence that ends at offset {offset} has no sequence delimitation item"
        ));
    }
    let group = cursor.u16(|| "an item tag".to_string())?;
    let element = cursor.u16(|| "an item tag".to_string())?;
    let length = cursor.u32(|| "an item length".to_string())?;
    if group != ITEM_GROUP {
        return Err(format!(
            "DICOM: expected an item tag (FFFE,E000) at offset {offset} but found \
             ({group:04X},{element:04X})"
        ));
    }
    Ok((element, length))
}

fn require_empty_delimiter(length: u32, offset: usize) -> Result<(), String> {
    if length == 0 {
        Ok(())
    } else {
        Err(format!(
            "DICOM: the delimitation item at offset {offset} has the non-zero length {length}"
        ))
    }
}

/// Walk the elements of an undefined-length item up to its item delimitation
/// item; the cursor is just after the item header.
fn scan_item(cursor: &mut Cursor<'_>, encoding: Encoding, depth: usize) -> Result<(), String> {
    loop {
        if cursor.rest().len() < 8 {
            return Err(format!(
                "DICOM: the item that ends at offset {} has no item delimitation item",
                cursor.pos
            ));
        }
        let header = read_header(cursor, encoding)?;
        if header.group == ITEM_GROUP {
            if header.element == ITEM_DELIMITATION {
                return require_empty_delimiter(header.length, header.offset);
            }
            return Err(format!(
                "DICOM: unexpected tag (FFFE,{:04X}) at offset {} inside an item",
                header.element, header.offset
            ));
        }
        if header.length == UNDEFINED_LENGTH {
            check_undefined_length(&header)?;
            scan_sequence(cursor, header.nested_encoding(encoding), depth + 1)?;
        } else {
            cursor.take(header.length_bytes(), || {
                format!("the value of ({:04X},{:04X})", header.group, header.element)
            })?;
        }
    }
}

/// VR of a tag in Implicit VR data, where the stream does not name it.
///
/// Only the tags this reader interprets and a few common ones are listed;
/// others are `UN`, which keeps their bytes untouched.  An element of
/// undefined length is a sequence and is recorded as `SQ` by the caller.
fn implicit_vr(group: u16, element: u16) -> ValueRepresentation {
    use ValueRepresentation as Vr;
    if element == 0x0000 {
        return Vr::UL;
    }
    match (group, element) {
        (0x0002, 0x0010) | (0x0008, 0x0016) | (0x0008, 0x0018) => Vr::UI,
        (0x0020, 0x000D) | (0x0020, 0x000E) => Vr::UI,
        (0x0008, 0x0020) => Vr::DA,
        (0x0008, 0x0030) => Vr::TM,
        (0x0008, 0x0060) | (0x0010, 0x0040) => Vr::CS,
        (0x0008, 0x0070) | (0x0008, 0x1030) | (0x0008, 0x103E) => Vr::LO,
        (0x0010, 0x0010) => Vr::PN,
        (0x0010, 0x0020) => Vr::LO,
        (0x0010, 0x0030) => Vr::DA,
        (0x0018, 0x0050) | (0x0018, 0x0088) => Vr::DS,
        (0x0020, 0x0010) => Vr::SH,
        (0x0020, 0x0011) | (0x0020, 0x0013) | (0x0028, 0x0008) => Vr::IS,
        (0x0020, 0x0032) | (0x0020, 0x0037) | (0x0020, 0x1041) => Vr::DS,
        (0x0028, 0x0002) | (0x0028, 0x0006) | (0x0028, 0x0010) | (0x0028, 0x0011) => Vr::US,
        (0x0028, 0x0100..=0x0103) => Vr::US,
        (0x0028, 0x0004) => Vr::CS,
        (0x0028, 0x0030) => Vr::DS,
        (0x0028, 0x1050..=0x1053) => Vr::DS,
        (0x7FE0, 0x0010) => Vr::OW,
        _ => Vr::UN,
    }
}

/// A read position in a byte slice; every read is bounds-checked.
struct Cursor<'a> {
    data: &'a [u8],
    pos: usize,
}

impl<'a> Cursor<'a> {
    fn new(data: &'a [u8], pos: usize) -> Self {
        Self { data, pos }
    }

    /// The bytes not read yet.
    fn rest(&self) -> &'a [u8] {
        self.data.get(self.pos..).unwrap_or_default()
    }

    /// `data[from..to]`, if it is in range.
    fn slice(&self, from: usize, to: usize) -> Option<&'a [u8]> {
        self.data.get(from..to)
    }

    /// The tag at the cursor, without consuming it.
    fn peek_tag(&self) -> Option<(u16, u16)> {
        let rest = self.rest();
        let bytes = rest.get(..4)?;
        Some((
            u16::from_le_bytes([bytes[0], bytes[1]]),
            u16::from_le_bytes([bytes[2], bytes[3]]),
        ))
    }

    /// Take `len` bytes; `what` names them in the error.
    fn take(&mut self, len: usize, what: impl FnOnce() -> String) -> Result<&'a [u8], String> {
        let end = self
            .pos
            .checked_add(len)
            .filter(|&end| end <= self.data.len());
        let Some(end) = end else {
            return Err(format!(
                "DICOM: {} at offset {} needs {len} bytes but only {} remain",
                what(),
                self.pos,
                self.rest().len()
            ));
        };
        let bytes = self.data.get(self.pos..end).unwrap_or_default();
        self.pos = end;
        Ok(bytes)
    }

    fn array<const N: usize>(&mut self, what: impl FnOnce() -> String) -> Result<[u8; N], String> {
        let bytes = self.take(N, what)?;
        <[u8; N]>::try_from(bytes).map_err(|_| "DICOM: internal error: array length".to_string())
    }

    fn u16(&mut self, what: impl FnOnce() -> String) -> Result<u16, String> {
        self.array(what).map(u16::from_le_bytes)
    }

    fn u32(&mut self, what: impl FnOnce() -> String) -> Result<u32, String> {
        self.array(what).map(u32::from_le_bytes)
    }
}
