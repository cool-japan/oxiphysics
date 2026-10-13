// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! LAS (ASPRS LiDAR exchange format) binary point-cloud reader and a minimal
//! writer.
//!
//! [`read_las_binary`] parses the public header block of every LAS revision
//! from 1.0 to 1.4 (227-, 235- and 375-byte headers) and decodes point data
//! record formats (PDRF) 0 to 10:
//!
//! | PDRF | min. record | extra fields beyond the base record |
//! |------|-------------|--------------------------------------|
//! | 0    | 20 bytes    | none                                 |
//! | 1    | 28 bytes    | GPS time                             |
//! | 2    | 26 bytes    | RGB                                  |
//! | 3    | 34 bytes    | GPS time, RGB                        |
//! | 4    | 57 bytes    | GPS time, waveform packet            |
//! | 5    | 63 bytes    | GPS time, RGB, waveform packet       |
//! | 6    | 30 bytes    | GPS time (LAS 1.4 base record)       |
//! | 7    | 36 bytes    | GPS time, RGB                        |
//! | 8    | 38 bytes    | GPS time, RGB, NIR                   |
//! | 9    | 59 bytes    | GPS time, waveform packet            |
//! | 10   | 67 bytes    | GPS time, RGB, NIR, waveform packet  |
//!
//! Waveform packet descriptors are skipped (the record length is honoured so
//! the following record is found correctly) but not decoded. LAZ-compressed
//! files are rejected.

// ── Public header layout (byte offsets from the LAS 1.4 specification) ──────

/// Size of the public header block of LAS 1.0 - 1.2 (the minimum for any file).
const LAS_MIN_HEADER_SIZE: usize = 227;
/// Size of the public header block of LAS 1.4.
const LAS_14_HEADER_SIZE: usize = 375;

const OFF_VERSION_MAJOR: usize = 24;
const OFF_VERSION_MINOR: usize = 25;
const OFF_CREATION_DOY: usize = 90;
const OFF_CREATION_YEAR: usize = 92;
const OFF_HEADER_SIZE: usize = 94;
const OFF_POINT_DATA_OFFSET: usize = 96;
const OFF_POINT_FORMAT: usize = 104;
const OFF_RECORD_LEN: usize = 105;
/// Legacy (u32) number of point records; the only count before LAS 1.4.
const OFF_LEGACY_COUNT: usize = 107;
const OFF_SCALE: usize = 131;
const OFF_OFFSET: usize = 155;
/// LAS 1.4 only: 64-bit number of point records.
const OFF_EXTENDED_COUNT: usize = 247;

/// Degrees per count of the 16-bit scan angle of PDRF 6 - 10.
const EXTENDED_SCAN_ANGLE_DEGREES: f64 = 0.006;

// ── Public types ────────────────────────────────────────────────────────────

/// Public header block parsed from a LAS 1.0 - 1.4 binary file.
#[derive(Debug, Clone)]
pub struct LasHeader {
    /// File creation year (GPS week or calendar year).
    pub file_creation_year: u16,
    /// File creation day of year.
    pub file_creation_doy: u16,
    /// LAS spec major version (always 1).
    pub version_major: u8,
    /// LAS spec minor version (0 - 4).
    pub version_minor: u8,
    /// Size of the public header block in bytes (227, 235 or 375).
    pub header_size: u16,
    /// Number of point records. For LAS 1.4 this is the 64-bit count (falling
    /// back to the legacy 32-bit count when that is zero); before LAS 1.4 it is
    /// the legacy 32-bit count.
    pub point_count: u64,
    /// X scale factor.
    pub x_scale: f64,
    /// Y scale factor.
    pub y_scale: f64,
    /// Z scale factor.
    pub z_scale: f64,
    /// X offset.
    pub x_offset: f64,
    /// Y offset.
    pub y_offset: f64,
    /// Z offset.
    pub z_offset: f64,
    /// Point data record format ID (0–10).
    pub point_format: u8,
    /// Number of bytes per point record.
    pub point_record_len: u16,
    /// Byte offset to the first point record.
    pub point_data_offset: u32,
}

/// One decoded LAS point.
///
/// The bit layout of [`return_flags`](Self::return_flags),
/// [`classification`](Self::classification) and [`flags`](Self::flags) depends
/// on the point data record format; the raw bytes are exposed unchanged.
#[derive(Debug, Clone, PartialEq)]
pub struct LasPoint {
    /// World-space X coordinate (metres or survey feet).
    pub x: f64,
    /// World-space Y coordinate.
    pub y: f64,
    /// World-space Z coordinate.
    pub z: f64,
    /// Return intensity (0–65535).
    pub intensity: u16,
    /// Raw return byte (record offset 14).
    ///
    /// PDRF 0–5: bits 0–2 return number, bits 3–5 number of returns, bit 6 scan
    /// direction flag, bit 7 edge of flight line. PDRF 6–10: bits 0–3 return
    /// number, bits 4–7 number of returns.
    pub return_flags: u8,
    /// Raw flag byte of PDRF 6–10 (record offset 15): bits 0–3 classification
    /// flags (synthetic, key-point, withheld, overlap), bits 4–5 scanner
    /// channel, bit 6 scan direction flag, bit 7 edge of flight line. Always `0`
    /// for PDRF 0–5, whose flags live in `return_flags` and in the top three
    /// bits of `classification`.
    pub flags: u8,
    /// Classification byte. PDRF 6–10: the classification code. PDRF 0–5: bits
    /// 0–4 are the code and bits 5–7 the synthetic / key-point / withheld
    /// flags.
    pub classification: u8,
    /// Scan angle in degrees. PDRF 0–5 store an integer scan angle rank;
    /// PDRF 6–10 store a signed 16-bit count of 0.006 degrees.
    pub scan_angle: f32,
    /// User data.
    pub user_data: u8,
    /// Point source ID.
    pub point_source_id: u16,
    /// GPS time (PDRF 1, 3, 4, 5 and 6–10), `None` for PDRF 0 and 2.
    pub gps_time: Option<f64>,
    /// 16-bit red / green / blue (PDRF 2, 3, 5, 7, 8 and 10).
    pub rgb: Option<[u16; 3]>,
    /// 16-bit near-infrared channel (PDRF 8 and 10).
    pub nir: Option<u16>,
}

/// A decoded LAS point cloud (header + points).
#[derive(Debug, Clone)]
pub struct LasPointCloud {
    /// Parsed LAS header.
    pub header: LasHeader,
    /// Decoded point records.
    pub points: Vec<LasPoint>,
}

// ── Point data record formats ───────────────────────────────────────────────

/// Static byte layout of one point data record format.
struct PointFormatLayout {
    /// Minimum record length in bytes (the record may carry extra bytes).
    min_len: u16,
    /// `true` for the LAS 1.4 base record of PDRF 6–10.
    extended: bool,
    /// Offset of the GPS time field.
    gps_time: Option<usize>,
    /// Offset of the first of three RGB channels.
    rgb: Option<usize>,
    /// Offset of the NIR channel.
    nir: Option<usize>,
}

/// Look up the layout of PDRF `format`, or `None` if it is not 0–10.
fn point_format_layout(format: u8) -> Option<PointFormatLayout> {
    let (min_len, extended, gps_time, rgb, nir) = match format {
        0 => (20, false, None, None, None),
        1 => (28, false, Some(20), None, None),
        2 => (26, false, None, Some(20), None),
        3 => (34, false, Some(20), Some(28), None),
        4 => (57, false, Some(20), None, None),
        5 => (63, false, Some(20), Some(28), None),
        6 => (30, true, Some(22), None, None),
        7 => (36, true, Some(22), Some(30), None),
        8 => (38, true, Some(22), Some(30), Some(36)),
        9 => (59, true, Some(22), None, None),
        10 => (67, true, Some(22), Some(30), Some(36)),
        _ => return None,
    };
    Some(PointFormatLayout {
        min_len,
        extended,
        gps_time,
        rgb,
        nir,
    })
}

// ── Internal little-endian readers (absolute offsets, bounds-checked) ───────

/// Borrow `N` bytes at `off`, or fail with an EOF error naming `what`.
fn le_bytes<const N: usize>(data: &[u8], off: usize, what: &str) -> Result<[u8; N], String> {
    let end = off
        .checked_add(N)
        .ok_or_else(|| format!("LAS: offset overflow reading {what}"))?;
    let slice = data
        .get(off..end)
        .ok_or_else(|| format!("LAS: unexpected EOF reading {what} at offset {off}"))?;
    <[u8; N]>::try_from(slice).map_err(|_| format!("LAS: cannot read {what} at offset {off}"))
}

fn le_u8(data: &[u8], off: usize, what: &str) -> Result<u8, String> {
    le_bytes::<1>(data, off, what).map(|b| b[0])
}

fn le_u16(data: &[u8], off: usize, what: &str) -> Result<u16, String> {
    le_bytes::<2>(data, off, what).map(u16::from_le_bytes)
}

fn le_i16(data: &[u8], off: usize, what: &str) -> Result<i16, String> {
    le_bytes::<2>(data, off, what).map(i16::from_le_bytes)
}

fn le_u32(data: &[u8], off: usize, what: &str) -> Result<u32, String> {
    le_bytes::<4>(data, off, what).map(u32::from_le_bytes)
}

fn le_i32(data: &[u8], off: usize, what: &str) -> Result<i32, String> {
    le_bytes::<4>(data, off, what).map(i32::from_le_bytes)
}

fn le_u64(data: &[u8], off: usize, what: &str) -> Result<u64, String> {
    le_bytes::<8>(data, off, what).map(u64::from_le_bytes)
}

fn le_f64(data: &[u8], off: usize, what: &str) -> Result<f64, String> {
    le_bytes::<8>(data, off, what).map(f64::from_le_bytes)
}

// ── Reader ──────────────────────────────────────────────────────────────────

/// Parse a LAS 1.0 – 1.4 binary point cloud from raw bytes.
///
/// Validates the `"LASF"` file signature and the version (`1.0` – `1.4`),
/// decodes the public header block whose size is given by its own *header
/// size* field (227 bytes for LAS 1.0 – 1.2, 235 for 1.3, 375 for 1.4), then
/// decodes every point record according to the file's point data record format
/// (0 – 10, see the [module documentation](self)). The 64-bit point count at
/// byte 247 is only read for LAS 1.4 files whose header is at least 375 bytes
/// long; older files use the 32-bit legacy count.
///
/// Records are read with the stride given by the header's *point data record
/// length*, so records carrying extra bytes are supported.
///
/// # Errors
///
/// Returns a descriptive `String` error if the data is shorter than the
/// header it declares, the signature or version is wrong, the point format is
/// unsupported (greater than 10, or LAZ-compressed), the record length is
/// smaller than the format requires, the point data offset lies inside the
/// header, or the file holds fewer point records than the header declares.
pub fn read_las_binary(data: &[u8]) -> Result<LasPointCloud, String> {
    if data.len() < LAS_MIN_HEADER_SIZE {
        return Err(format!(
            "LAS: data too short ({} bytes, need at least {LAS_MIN_HEADER_SIZE})",
            data.len()
        ));
    }

    // ── File signature "LASF" (bytes 0-3) ───────────────────────────────────
    let signature = le_bytes::<4>(data, 0, "file signature")?;
    if &signature != b"LASF" {
        return Err(format!(
            "LAS: bad signature {signature:?}, expected b\"LASF\""
        ));
    }

    let version_major = le_u8(data, OFF_VERSION_MAJOR, "version major")?;
    let version_minor = le_u8(data, OFF_VERSION_MINOR, "version minor")?;
    if version_major != 1 || version_minor > 4 {
        return Err(format!(
            "LAS: unsupported version {version_major}.{version_minor} (supported: 1.0 - 1.4)"
        ));
    }

    let header_size = le_u16(data, OFF_HEADER_SIZE, "header size")?;
    let header_len = usize::from(header_size);
    if header_len < LAS_MIN_HEADER_SIZE {
        return Err(format!(
            "LAS: header size {header_size} is below the {LAS_MIN_HEADER_SIZE}-byte minimum"
        ));
    }
    if header_len > data.len() {
        return Err(format!(
            "LAS: header declares {header_size} bytes but the data has only {}",
            data.len()
        ));
    }

    let file_creation_doy = le_u16(data, OFF_CREATION_DOY, "file creation day")?;
    let file_creation_year = le_u16(data, OFF_CREATION_YEAR, "file creation year")?;
    let point_data_offset = le_u32(data, OFF_POINT_DATA_OFFSET, "offset to point data")?;
    if (point_data_offset as usize) < header_len {
        return Err(format!(
            "LAS: offset to point data ({point_data_offset}) lies inside the {header_size}-byte header"
        ));
    }

    let point_format = le_u8(data, OFF_POINT_FORMAT, "point data record format")?;
    if point_format & 0x80 != 0 {
        return Err(format!(
            "LAS: LAZ-compressed point data (format byte 0x{point_format:02x}) is not supported"
        ));
    }
    let layout = point_format_layout(point_format).ok_or_else(|| {
        format!("LAS: unsupported point data record format {point_format} (expected 0-10)")
    })?;
    let point_record_len = le_u16(data, OFF_RECORD_LEN, "point data record length")?;
    if point_record_len < layout.min_len {
        return Err(format!(
            "LAS: record length {point_record_len} is below the {} bytes point format {point_format} requires",
            layout.min_len
        ));
    }

    // Before LAS 1.4 the 32-bit legacy count is the only count; the bytes at
    // offset 247 belong to the VLRs or the point data there.
    let legacy_count = le_u32(data, OFF_LEGACY_COUNT, "legacy point count")?;
    let extended_count = if version_minor >= 4 && header_len >= LAS_14_HEADER_SIZE {
        le_u64(data, OFF_EXTENDED_COUNT, "extended point count")?
    } else {
        0
    };
    let point_count = if extended_count > 0 {
        extended_count
    } else {
        u64::from(legacy_count)
    };

    let x_scale = le_f64(data, OFF_SCALE, "x scale")?;
    let y_scale = le_f64(data, OFF_SCALE + 8, "y scale")?;
    let z_scale = le_f64(data, OFF_SCALE + 16, "z scale")?;
    let x_offset = le_f64(data, OFF_OFFSET, "x offset")?;
    let y_offset = le_f64(data, OFF_OFFSET + 8, "y offset")?;
    let z_offset = le_f64(data, OFF_OFFSET + 16, "z offset")?;

    let header = LasHeader {
        file_creation_year,
        file_creation_doy,
        version_major,
        version_minor,
        header_size,
        point_count,
        x_scale,
        y_scale,
        z_scale,
        x_offset,
        y_offset,
        z_offset,
        point_format,
        point_record_len,
        point_data_offset,
    };

    // ── Point records ────────────────────────────────────────────────────────
    let rec_len = usize::from(point_record_len);
    let start = point_data_offset as usize;
    // Number of complete records the file can hold; bounds the declared count
    // before anything is allocated from it.
    let available = (data.len().saturating_sub(start) / rec_len) as u64;
    if point_count > available {
        return Err(format!(
            "LAS: truncated at point {available}/{point_count} (point data starts at byte {start}, {rec_len} bytes per record, file is {} bytes)",
            data.len()
        ));
    }

    let mut points: Vec<LasPoint> = Vec::with_capacity(point_count as usize);
    for index in 0..point_count as usize {
        let p_off = start + index * rec_len;
        let record = data
            .get(p_off..p_off + rec_len)
            .ok_or_else(|| format!("LAS: truncated at point {index}/{point_count}"))?;
        points.push(decode_point(record, &layout, &header)?);
    }

    Ok(LasPointCloud { header, points })
}

/// Decode one point record (`record.len()` is at least `layout.min_len`).
fn decode_point(
    record: &[u8],
    layout: &PointFormatLayout,
    header: &LasHeader,
) -> Result<LasPoint, String> {
    let raw_x = le_i32(record, 0, "point x")?;
    let raw_y = le_i32(record, 4, "point y")?;
    let raw_z = le_i32(record, 8, "point z")?;
    let intensity = le_u16(record, 12, "intensity")?;
    let return_flags = le_u8(record, 14, "return flags")?;

    let (flags, classification, scan_angle, user_data, point_source_id) = if layout.extended {
        let flags = le_u8(record, 15, "classification flags")?;
        let classification = le_u8(record, 16, "classification")?;
        let user_data = le_u8(record, 17, "user data")?;
        let angle_counts = le_i16(record, 18, "scan angle")?;
        let point_source_id = le_u16(record, 20, "point source id")?;
        (
            flags,
            classification,
            (f64::from(angle_counts) * EXTENDED_SCAN_ANGLE_DEGREES) as f32,
            user_data,
            point_source_id,
        )
    } else {
        let classification = le_u8(record, 15, "classification")?;
        let scan_angle_rank = le_u8(record, 16, "scan angle rank")? as i8;
        let user_data = le_u8(record, 17, "user data")?;
        let point_source_id = le_u16(record, 18, "point source id")?;
        (
            0,
            classification,
            f32::from(scan_angle_rank),
            user_data,
            point_source_id,
        )
    };

    let gps_time = layout
        .gps_time
        .map(|off| le_f64(record, off, "gps time"))
        .transpose()?;
    let rgb = layout
        .rgb
        .map(|off| -> Result<[u16; 3], String> {
            Ok([
                le_u16(record, off, "red")?,
                le_u16(record, off + 2, "green")?,
                le_u16(record, off + 4, "blue")?,
            ])
        })
        .transpose()?;
    let nir = layout
        .nir
        .map(|off| le_u16(record, off, "nir"))
        .transpose()?;

    Ok(LasPoint {
        x: f64::from(raw_x) * header.x_scale + header.x_offset,
        y: f64::from(raw_y) * header.y_scale + header.y_offset,
        z: f64::from(raw_z) * header.z_scale + header.z_offset,
        intensity,
        return_flags,
        flags,
        classification,
        scan_angle,
        user_data,
        point_source_id,
        gps_time,
        rgb,
        nir,
    })
}

// ── Writer ──────────────────────────────────────────────────────────────────

/// Write a minimal LAS 1.4 binary file from a slice of world-space points.
///
/// Uses Point Data Format 0 with a 375-byte header and a scale factor of
/// `scale` (default 0.001 m = 1 mm precision) and no offset.  Intended only
/// for round-trip testing of [`read_las_binary`].
pub fn write_las_binary_minimal(points: &[[f64; 3]], scale: f64) -> Vec<u8> {
    let n = points.len() as u64;
    let rec_len: u16 = 20; // Point Format 0 minimum
    let header_size: u16 = 375;
    let point_data_offset: u32 = header_size as u32;
    let mut buf = Vec::with_capacity(header_size as usize + n as usize * rec_len as usize);

    // File signature
    buf.extend_from_slice(b"LASF");
    // File source ID + global encoding
    buf.extend_from_slice(&0u16.to_le_bytes());
    buf.extend_from_slice(&0u16.to_le_bytes());
    // Project GUID (16 zero bytes)
    buf.extend_from_slice(&[0u8; 16]);
    // Version: 1.4
    buf.push(1u8);
    buf.push(4u8);
    // System identifier (32 bytes)
    buf.extend_from_slice(&[0u8; 32]);
    // Generating software (32 bytes) — fill with "oxiphysics"
    let mut sw = [b' '; 32];
    let tag = b"oxiphysics";
    sw[..tag.len()].copy_from_slice(tag);
    buf.extend_from_slice(&sw);
    // File creation DOY / year
    buf.extend_from_slice(&1u16.to_le_bytes()); // DOY
    buf.extend_from_slice(&2026u16.to_le_bytes()); // year
    // Header size
    buf.extend_from_slice(&header_size.to_le_bytes());
    // Offset to point data
    buf.extend_from_slice(&point_data_offset.to_le_bytes());
    // Number of VLRs
    buf.extend_from_slice(&0u32.to_le_bytes());
    // Point data format
    buf.push(0u8); // format 0
    // Point data record length
    buf.extend_from_slice(&rec_len.to_le_bytes());
    // Legacy point count
    let legacy_n = n.min(u32::MAX as u64) as u32;
    buf.extend_from_slice(&legacy_n.to_le_bytes());
    // Legacy point counts by return (5 × u32)
    buf.extend_from_slice(&[0u8; 20]);
    // Scale factors
    buf.extend_from_slice(&scale.to_le_bytes());
    buf.extend_from_slice(&scale.to_le_bytes());
    buf.extend_from_slice(&scale.to_le_bytes());
    // Offsets (0, 0, 0)
    buf.extend_from_slice(&0.0f64.to_le_bytes());
    buf.extend_from_slice(&0.0f64.to_le_bytes());
    buf.extend_from_slice(&0.0f64.to_le_bytes());
    // Max/min X, Y, Z — compute from points
    let (mut max_x, mut min_x) = (f64::NEG_INFINITY, f64::INFINITY);
    let (mut max_y, mut min_y) = (f64::NEG_INFINITY, f64::INFINITY);
    let (mut max_z, mut min_z) = (f64::NEG_INFINITY, f64::INFINITY);
    for &[x, y, z] in points {
        if x > max_x {
            max_x = x;
        }
        if x < min_x {
            min_x = x;
        }
        if y > max_y {
            max_y = y;
        }
        if y < min_y {
            min_y = y;
        }
        if z > max_z {
            max_z = z;
        }
        if z < min_z {
            min_z = z;
        }
    }
    if points.is_empty() {
        max_x = 0.0;
        min_x = 0.0;
        max_y = 0.0;
        min_y = 0.0;
        max_z = 0.0;
        min_z = 0.0;
    }
    buf.extend_from_slice(&max_x.to_le_bytes());
    buf.extend_from_slice(&min_x.to_le_bytes());
    buf.extend_from_slice(&max_y.to_le_bytes());
    buf.extend_from_slice(&min_y.to_le_bytes());
    buf.extend_from_slice(&max_z.to_le_bytes());
    buf.extend_from_slice(&min_z.to_le_bytes());
    // Start of Waveform Data Packet Record (u64)
    buf.extend_from_slice(&0u64.to_le_bytes());
    // Start of first EVLR (u64)
    buf.extend_from_slice(&0u64.to_le_bytes());
    // Number of EVLRs (u32)
    buf.extend_from_slice(&0u32.to_le_bytes());
    // LAS 1.4 extended point count (u64)
    buf.extend_from_slice(&n.to_le_bytes());
    // Extended point counts by return (15 × u64)
    buf.extend_from_slice(&[0u8; 120]);

    // Sanity-check: buf.len() must equal header_size at this point.
    assert_eq!(buf.len(), header_size as usize, "header length mismatch");

    // ── Point records (Point Format 0) ────────────────────────────────────────
    for &[x, y, z] in points {
        let raw_x = ((x) / scale).round() as i32;
        let raw_y = ((y) / scale).round() as i32;
        let raw_z = ((z) / scale).round() as i32;
        buf.extend_from_slice(&raw_x.to_le_bytes());
        buf.extend_from_slice(&raw_y.to_le_bytes());
        buf.extend_from_slice(&raw_z.to_le_bytes());
        buf.extend_from_slice(&0u16.to_le_bytes()); // intensity
        buf.push(0u8); // return flags
        buf.push(0u8); // classification
        buf.push(0u8); // scan angle (i8 → u8 cast fine for 0)
        buf.push(0u8); // user data
        buf.extend_from_slice(&0u16.to_le_bytes()); // point source ID
    }

    buf
}

// ── Tests ───────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // The fixtures below are laid out straight from the LAS specification
    // tables (hard-coded byte offsets), independently of the reader's constants.

    fn put_u16(buf: &mut [u8], off: usize, v: u16) {
        buf[off..off + 2].copy_from_slice(&v.to_le_bytes());
    }
    fn put_u32(buf: &mut [u8], off: usize, v: u32) {
        buf[off..off + 4].copy_from_slice(&v.to_le_bytes());
    }
    fn put_u64(buf: &mut [u8], off: usize, v: u64) {
        buf[off..off + 8].copy_from_slice(&v.to_le_bytes());
    }
    fn put_f64(buf: &mut [u8], off: usize, v: f64) {
        buf[off..off + 8].copy_from_slice(&v.to_le_bytes());
    }

    /// Scale 0.01 on every axis, offsets (100, 200, 300).
    const SCALE: f64 = 0.01;
    const OFFSETS: [f64; 3] = [100.0, 200.0, 300.0];

    /// Parameters of a synthetic public header block.
    #[derive(Clone, Copy)]
    struct Spec {
        version_minor: u8,
        header_size: usize,
        point_format: u8,
        record_len: u16,
        legacy_count: u32,
        extended_count: u64,
        point_data_offset: u32,
    }

    impl Spec {
        fn new(version_minor: u8, header_size: usize, point_format: u8, record_len: u16) -> Self {
            Self {
                version_minor,
                header_size,
                point_format,
                record_len,
                legacy_count: 0,
                extended_count: 0,
                point_data_offset: header_size as u32,
            }
        }
        fn count(mut self, legacy: u32, extended: u64) -> Self {
            self.legacy_count = legacy;
            self.extended_count = extended;
            self
        }
    }

    fn header_bytes(spec: Spec) -> Vec<u8> {
        let mut h = vec![0u8; spec.header_size];
        h[0..4].copy_from_slice(b"LASF");
        h[24] = 1; // version major
        h[25] = spec.version_minor;
        put_u16(&mut h, 90, 33); // day of year
        put_u16(&mut h, 92, 2025); // year
        put_u16(&mut h, 94, spec.header_size as u16);
        put_u32(&mut h, 96, spec.point_data_offset);
        h[104] = spec.point_format;
        put_u16(&mut h, 105, spec.record_len);
        put_u32(&mut h, 107, spec.legacy_count);
        for (axis, offset) in OFFSETS.iter().enumerate() {
            put_f64(&mut h, 131 + 8 * axis, SCALE);
            put_f64(&mut h, 155 + 8 * axis, *offset);
        }
        if spec.header_size >= 375 {
            put_u64(&mut h, 247, spec.extended_count);
        }
        h
    }

    /// PDRF 0 record (20 bytes).
    fn record0(
        xyz: [i32; 3],
        intensity: u16,
        return_byte: u8,
        class_byte: u8,
        angle_rank: i8,
        user: u8,
        source: u16,
    ) -> Vec<u8> {
        let mut r = Vec::new();
        for v in xyz {
            r.extend_from_slice(&v.to_le_bytes());
        }
        r.extend_from_slice(&intensity.to_le_bytes());
        r.push(return_byte);
        r.push(class_byte);
        r.push(angle_rank as u8);
        r.push(user);
        r.extend_from_slice(&source.to_le_bytes());
        assert_eq!(r.len(), 20);
        r
    }

    /// PDRF 6 base record (30 bytes).
    fn record6(
        xyz: [i32; 3],
        intensity: u16,
        return_byte: u8,
        flag_byte: u8,
        class: u8,
        user: u8,
        angle_counts: i16,
        source: u16,
        gps: f64,
    ) -> Vec<u8> {
        let mut r = Vec::new();
        for v in xyz {
            r.extend_from_slice(&v.to_le_bytes());
        }
        r.extend_from_slice(&intensity.to_le_bytes());
        r.push(return_byte);
        r.push(flag_byte);
        r.push(class);
        r.push(user);
        r.extend_from_slice(&angle_counts.to_le_bytes());
        r.extend_from_slice(&source.to_le_bytes());
        r.extend_from_slice(&gps.to_le_bytes());
        assert_eq!(r.len(), 30);
        r
    }

    fn world(raw: [i32; 3]) -> [f64; 3] {
        [
            f64::from(raw[0]) * SCALE + OFFSETS[0],
            f64::from(raw[1]) * SCALE + OFFSETS[1],
            f64::from(raw[2]) * SCALE + OFFSETS[2],
        ]
    }

    fn assert_xyz(p: &LasPoint, raw: [i32; 3]) {
        let w = world(raw);
        assert!((p.x - w[0]).abs() < 1e-9, "x {} vs {}", p.x, w[0]);
        assert!((p.y - w[1]).abs() < 1e-9, "y {} vs {}", p.y, w[1]);
        assert!((p.z - w[2]).abs() < 1e-9, "z {} vs {}", p.z, w[2]);
    }

    // -- original contract ----------------------------------------------------

    #[test]
    fn test_read_las_bad_signature() {
        let data = b"NOTL\x00\x00\x00\x00";
        assert!(read_las_binary(data).is_err());
    }

    #[test]
    fn test_read_las_bad_signature_full_header() {
        let mut data = header_bytes(Spec::new(2, 227, 0, 20));
        data[0..4].copy_from_slice(b"NOTL");
        let err = read_las_binary(&data).unwrap_err();
        assert!(err.contains("bad signature"), "{err}");
    }

    #[test]
    fn test_read_las_too_short() {
        let data = vec![0u8; 100];
        assert!(read_las_binary(&data).is_err());
    }

    #[test]
    fn test_las_roundtrip_empty() {
        let data = write_las_binary_minimal(&[], 0.001);
        let cloud = read_las_binary(&data).expect("should parse empty LAS");
        assert_eq!(cloud.points.len(), 0);
        assert_eq!(cloud.header.version_major, 1);
        assert_eq!(cloud.header.version_minor, 4);
        assert_eq!(cloud.header.header_size, 375);
    }

    #[test]
    fn test_las_roundtrip_single_point() {
        let pts = vec![[1.234, 5.678, -0.5]];
        let data = write_las_binary_minimal(&pts, 0.001);
        let cloud = read_las_binary(&data).expect("should parse single-point LAS");
        assert_eq!(cloud.points.len(), 1);
        assert!((cloud.points[0].x - 1.234).abs() < 0.002, "x mismatch");
        assert!((cloud.points[0].y - 5.678).abs() < 0.002, "y mismatch");
        assert!((cloud.points[0].z - (-0.5)).abs() < 0.002, "z mismatch");
    }

    #[test]
    fn test_las_roundtrip_multiple_points() {
        let pts: Vec<[f64; 3]> = (0..10)
            .map(|i| [i as f64 * 0.5, i as f64 * 0.25, i as f64 * 0.1])
            .collect();
        let data = write_las_binary_minimal(&pts, 0.001);
        let cloud = read_las_binary(&data).expect("should parse multi-point LAS");
        assert_eq!(cloud.points.len(), 10);
        for (i, pt) in cloud.points.iter().enumerate() {
            let expected_x = i as f64 * 0.5;
            assert!(
                (pt.x - expected_x).abs() < 0.002,
                "x mismatch at point {i}: expected {expected_x}, got {}",
                pt.x
            );
        }
    }

    #[test]
    fn test_las_writer_output_fields() {
        let data = write_las_binary_minimal(&[[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]], 0.001);
        let cloud = read_las_binary(&data).unwrap();
        let h = &cloud.header;
        assert_eq!(h.point_count, 2);
        assert_eq!(h.point_format, 0);
        assert_eq!(h.point_record_len, 20);
        assert_eq!(h.point_data_offset, 375);
        assert_eq!(h.file_creation_year, 2026);
        assert_eq!(h.file_creation_doy, 1);
        assert!((h.x_scale - 0.001).abs() < 1e-15);
        for p in &cloud.points {
            assert_eq!(p.gps_time, None);
            assert_eq!(p.rgb, None);
            assert_eq!(p.nir, None);
            assert_eq!(p.flags, 0);
        }
    }

    // -- pre-1.4 headers (F2435 / F2436) --------------------------------------

    #[test]
    fn las_v1_2_227_byte_header_is_read_with_the_legacy_count() {
        // 3 records straight after the 227-byte header. The record at index 1
        // starts at byte 247, where a LAS 1.4 header keeps its 64-bit point
        // count; reading it as a count gave (100000 | 200000 << 32), which
        // aborted in `Vec::with_capacity`.
        let spec = Spec::new(2, 227, 0, 20).count(3, 0);
        let raws = [[10, 20, 30], [100_000, 200_000, 300_000], [-5, -6, -7]];
        let mut data = header_bytes(spec);
        data.extend(record0(raws[0], 11, 0b0001_1001, 2, -12, 7, 1001));
        data.extend(record0(raws[1], 22, 0b0000_1010, 5, 90, 8, 1002));
        data.extend(record0(raws[2], 33, 0b0100_0011, 9, 0, 9, 1003));

        let cloud = read_las_binary(&data).expect("LAS 1.2 file must parse");
        assert_eq!(cloud.header.version_minor, 2);
        assert_eq!(cloud.header.header_size, 227);
        assert_eq!(cloud.header.point_count, 3);
        assert_eq!(cloud.header.point_data_offset, 227);
        assert_eq!(cloud.points.len(), 3);
        for (p, raw) in cloud.points.iter().zip(raws) {
            assert_xyz(p, raw);
        }
        let p0 = &cloud.points[0];
        assert_eq!(p0.intensity, 11);
        assert_eq!(p0.return_flags, 0b0001_1001);
        assert_eq!(p0.classification, 2);
        assert_eq!(p0.scan_angle, -12.0);
        assert_eq!(p0.user_data, 7);
        assert_eq!(p0.point_source_id, 1001);
        assert_eq!(cloud.points[1].scan_angle, 90.0);
        assert_eq!(cloud.points[2].point_source_id, 1003);
    }

    #[test]
    fn las_v1_2_with_a_huge_value_at_byte_247_still_uses_the_legacy_count() {
        // Same layout but with the bytes at 247..255 forced to u64::MAX.
        let spec = Spec::new(2, 227, 0, 20).count(2, 0);
        let mut data = header_bytes(spec);
        data.extend(record0([1, 2, 3], 0, 0, 0, 0, 0, 0));
        data.extend(record0([4, 5, 6], 0, 0, 0, 0, 0, 0));
        put_u64(&mut data, 247, u64::MAX);
        // Bytes 247..255 are the second record's x/y; the point count must not
        // come from there.
        let cloud = read_las_binary(&data).unwrap();
        assert_eq!(cloud.header.point_count, 2);
        assert_eq!(cloud.points.len(), 2);
    }

    #[test]
    fn las_v1_3_235_byte_header_with_gps_time_records() {
        let spec = Spec::new(3, 235, 1, 28).count(2, 0);
        let mut data = header_bytes(spec);
        for (i, raw) in [[1, 2, 3], [4, 5, 6]].into_iter().enumerate() {
            let mut rec = record0(raw, 5, 0, 1, 3, 4, 5);
            rec.extend_from_slice(&(1234.5 + i as f64).to_le_bytes());
            data.extend(rec);
        }
        let cloud = read_las_binary(&data).unwrap();
        assert_eq!(cloud.header.header_size, 235);
        assert_eq!(cloud.points.len(), 2);
        assert_eq!(cloud.points[0].gps_time, Some(1234.5));
        assert_eq!(cloud.points[1].gps_time, Some(1235.5));
        assert_eq!(cloud.points[1].rgb, None);
        assert_xyz(&cloud.points[1], [4, 5, 6]);
    }

    #[test]
    fn las_v1_2_pdrf2_and_pdrf3_carry_rgb() {
        let mut data = header_bytes(Spec::new(2, 227, 3, 34).count(1, 0));
        let mut rec = record0([7, 8, 9], 1, 2, 3, 4, 5, 6);
        rec.extend_from_slice(&99.25f64.to_le_bytes()); // gps time @20
        for c in [111u16, 222, 333] {
            rec.extend_from_slice(&c.to_le_bytes()); // rgb @28
        }
        data.extend(rec);
        let cloud = read_las_binary(&data).unwrap();
        let p = &cloud.points[0];
        assert_eq!(p.gps_time, Some(99.25));
        assert_eq!(p.rgb, Some([111, 222, 333]));

        // PDRF 2: RGB directly after the base record, no GPS time.
        let mut data = header_bytes(Spec::new(2, 227, 2, 26).count(1, 0));
        let mut rec = record0([7, 8, 9], 1, 2, 3, 4, 5, 6);
        for c in [10u16, 20, 30] {
            rec.extend_from_slice(&c.to_le_bytes());
        }
        data.extend(rec);
        let cloud = read_las_binary(&data).unwrap();
        assert_eq!(cloud.points[0].gps_time, None);
        assert_eq!(cloud.points[0].rgb, Some([10, 20, 30]));
    }

    // -- LAS 1.4 point formats 6-10 -------------------------------------------

    #[test]
    fn las_v1_4_pdrf6_decodes_the_extended_base_record() {
        // PDRF 6 has no legacy count (0) and a different field order than
        // PDRF 0: flags byte, classification byte, user data, a 16-bit scan
        // angle in 0.006 degree steps, then the source id and GPS time.
        let spec = Spec::new(4, 375, 6, 30).count(0, 2);
        let mut data = header_bytes(spec);
        data.extend(record6(
            [123, -456, 789],
            4321,
            0x32, // return 2 of 3
            0b1010_0101,
            6,
            77,
            -5000, // -30.0 degrees
            4242,
            86_400.125,
        ));
        data.extend(record6([1, 2, 3], 1, 0x11, 0, 2, 3, 15_000, 9, 0.5));

        let cloud = read_las_binary(&data).unwrap();
        assert_eq!(cloud.header.version_minor, 4);
        assert_eq!(cloud.header.point_format, 6);
        assert_eq!(cloud.header.point_count, 2);
        assert_eq!(cloud.points.len(), 2);

        let p = &cloud.points[0];
        assert_xyz(p, [123, -456, 789]);
        assert_eq!(p.intensity, 4321);
        assert_eq!(p.return_flags, 0x32);
        assert_eq!(p.flags, 0b1010_0101);
        assert_eq!(p.classification, 6);
        assert_eq!(p.user_data, 77);
        assert!((p.scan_angle + 30.0).abs() < 1e-4, "{}", p.scan_angle);
        assert_eq!(p.point_source_id, 4242);
        assert_eq!(p.gps_time, Some(86_400.125));
        assert_eq!(p.rgb, None);
        assert_eq!(p.nir, None);

        let q = &cloud.points[1];
        assert_xyz(q, [1, 2, 3]);
        assert!((q.scan_angle - 90.0).abs() < 1e-4, "{}", q.scan_angle);
        assert_eq!(q.gps_time, Some(0.5));
    }

    #[test]
    fn las_v1_4_pdrf7_and_pdrf8_carry_rgb_and_nir() {
        let spec7 = Spec::new(4, 375, 7, 36).count(0, 1);
        let mut data = header_bytes(spec7);
        let mut rec = record6([5, 6, 7], 8, 0x11, 0, 2, 0, 0, 1, 10.0);
        for c in [65535u16, 0, 321] {
            rec.extend_from_slice(&c.to_le_bytes());
        }
        data.extend(rec);
        let cloud = read_las_binary(&data).unwrap();
        assert_eq!(cloud.points[0].rgb, Some([65535, 0, 321]));
        assert_eq!(cloud.points[0].nir, None);

        let spec8 = Spec::new(4, 375, 8, 38).count(0, 1);
        let mut data = header_bytes(spec8);
        let mut rec = record6([5, 6, 7], 8, 0x11, 0, 2, 0, 0, 1, 10.0);
        for c in [1u16, 2, 3, 4444] {
            rec.extend_from_slice(&c.to_le_bytes());
        }
        data.extend(rec);
        let cloud = read_las_binary(&data).unwrap();
        assert_eq!(cloud.points[0].rgb, Some([1, 2, 3]));
        assert_eq!(cloud.points[0].nir, Some(4444));
    }

    #[test]
    fn las_v1_4_waveform_formats_advance_by_the_declared_record_length() {
        // PDRF 9 = PDRF 6 + 29 waveform bytes; the second record must be found
        // at +59, not +30.
        let mut data = header_bytes(Spec::new(4, 375, 9, 59).count(0, 2));
        for raw in [[11, 12, 13], [21, 22, 23]] {
            let mut rec = record6(raw, 0, 0x11, 0, 2, 0, 0, 0, 1.0);
            rec.extend_from_slice(&[0xAB; 29]);
            data.extend(rec);
        }
        let cloud = read_las_binary(&data).unwrap();
        assert_eq!(cloud.points.len(), 2);
        assert_xyz(&cloud.points[0], [11, 12, 13]);
        assert_xyz(&cloud.points[1], [21, 22, 23]);
    }

    #[test]
    fn las_v1_4_legacy_formats_still_use_the_64_bit_count() {
        // LAS 1.4 with PDRF 0: both counts set; the 64-bit one is used.
        let mut data = header_bytes(Spec::new(4, 375, 0, 20).count(2, 3));
        for raw in [[1, 1, 1], [2, 2, 2], [3, 3, 3]] {
            data.extend(record0(raw, 0, 0, 0, 0, 0, 0));
        }
        let cloud = read_las_binary(&data).unwrap();
        assert_eq!(cloud.header.point_count, 3);
        assert_eq!(cloud.points.len(), 3);
    }

    #[test]
    fn las_v1_4_zero_extended_count_falls_back_to_the_legacy_count() {
        let mut data = header_bytes(Spec::new(4, 375, 0, 20).count(2, 0));
        for raw in [[1, 1, 1], [2, 2, 2]] {
            data.extend(record0(raw, 0, 0, 0, 0, 0, 0));
        }
        let cloud = read_las_binary(&data).unwrap();
        assert_eq!(cloud.header.point_count, 2);
        assert_eq!(cloud.points.len(), 2);
    }

    #[test]
    fn las_scan_angle_resolution_of_extended_formats() {
        let mut data = header_bytes(Spec::new(4, 375, 6, 30).count(0, 3));
        for counts in [0i16, 15_000, -15_000] {
            data.extend(record6([0, 0, 0], 0, 0x11, 0, 0, 0, counts, 0, 0.0));
        }
        let cloud = read_las_binary(&data).unwrap();
        let angles: Vec<f32> = cloud.points.iter().map(|p| p.scan_angle).collect();
        assert_eq!(angles[0], 0.0);
        assert!((angles[1] - 90.0).abs() < 1e-4);
        assert!((angles[2] + 90.0).abs() < 1e-4);
    }

    // -- layout and framing ---------------------------------------------------

    #[test]
    fn las_point_format_table_matches_the_specification() {
        let min_len: Vec<u16> = (0..=10)
            .map(|f| point_format_layout(f).expect("formats 0-10 exist").min_len)
            .collect();
        assert_eq!(min_len, vec![20, 28, 26, 34, 57, 63, 30, 36, 38, 59, 67]);
        assert!(point_format_layout(11).is_none());
        assert!(point_format_layout(0x86).is_none());
    }

    #[test]
    fn las_records_with_extra_bytes_use_the_declared_stride() {
        // 4 user-defined extra bytes per record (stride 24 instead of 20).
        let mut data = header_bytes(Spec::new(2, 227, 0, 24).count(2, 0));
        for raw in [[1, 2, 3], [4, 5, 6]] {
            let mut rec = record0(raw, 9, 0, 0, 0, 0, 0);
            rec.extend_from_slice(&[0xEE; 4]);
            data.extend(rec);
        }
        let cloud = read_las_binary(&data).unwrap();
        assert_eq!(cloud.header.point_record_len, 24);
        assert_xyz(&cloud.points[1], [4, 5, 6]);
    }

    #[test]
    fn las_point_data_offset_after_vlrs_is_honoured() {
        // 50 bytes of variable length record data between header and points.
        let mut spec = Spec::new(2, 227, 0, 20).count(1, 0);
        spec.point_data_offset = 227 + 50;
        let mut data = header_bytes(spec);
        data.extend_from_slice(&[0xCD; 50]);
        data.extend(record0([3, 2, 1], 0, 0, 0, 0, 0, 0));
        let cloud = read_las_binary(&data).unwrap();
        assert_eq!(cloud.points.len(), 1);
        assert_xyz(&cloud.points[0], [3, 2, 1]);
    }

    // -- malformed input ------------------------------------------------------

    #[test]
    fn las_truncated_point_data_is_an_error() {
        let mut data = header_bytes(Spec::new(2, 227, 0, 20).count(5, 0));
        for raw in [[1, 1, 1], [2, 2, 2], [3, 3, 3]] {
            data.extend(record0(raw, 0, 0, 0, 0, 0, 0));
        }
        // One byte short of the 4th record as well.
        data.extend_from_slice(&[0u8; 19]);
        let err = read_las_binary(&data).unwrap_err();
        assert!(err.contains("truncated at point 3/5"), "{err}");
    }

    #[test]
    fn las_absurd_declared_count_is_an_error_not_an_allocation() {
        for count in [u64::MAX, u64::MAX / 20, 1 << 40] {
            let data = header_bytes(Spec::new(4, 375, 0, 20).count(0, count));
            let err = read_las_binary(&data).unwrap_err();
            assert!(err.contains("truncated"), "{err}");
        }
        // Legacy 32-bit count: 4 billion records in a header-only file.
        let data = header_bytes(Spec::new(2, 227, 0, 20).count(u32::MAX, 0));
        assert!(read_las_binary(&data).is_err());
    }

    #[test]
    fn las_header_longer_than_the_data_is_an_error() {
        let mut data = header_bytes(Spec::new(4, 375, 0, 20));
        data.truncate(300);
        let err = read_las_binary(&data).unwrap_err();
        assert!(err.contains("header declares 375 bytes"), "{err}");
    }

    #[test]
    fn las_header_size_below_the_minimum_is_an_error() {
        let mut data = header_bytes(Spec::new(2, 300, 0, 20));
        put_u16(&mut data, 94, 100);
        let err = read_las_binary(&data).unwrap_err();
        assert!(err.contains("below the 227-byte minimum"), "{err}");
    }

    #[test]
    fn las_unsupported_point_formats_are_errors() {
        let data = header_bytes(Spec::new(4, 375, 11, 67));
        let err = read_las_binary(&data).unwrap_err();
        assert!(
            err.contains("unsupported point data record format 11"),
            "{err}"
        );

        // LASzip marks compressed data with bit 7 of the format byte.
        let data = header_bytes(Spec::new(2, 227, 0x80 | 1, 28));
        let err = read_las_binary(&data).unwrap_err();
        assert!(err.contains("LAZ"), "{err}");
    }

    #[test]
    fn las_record_length_below_the_format_minimum_is_an_error() {
        let data = header_bytes(Spec::new(4, 375, 6, 20));
        let err = read_las_binary(&data).unwrap_err();
        assert!(err.contains("record length 20"), "{err}");
    }

    #[test]
    fn las_point_data_offset_inside_the_header_is_an_error() {
        let mut spec = Spec::new(2, 227, 0, 20);
        spec.point_data_offset = 100;
        let err = read_las_binary(&header_bytes(spec)).unwrap_err();
        assert!(err.contains("lies inside"), "{err}");
    }

    #[test]
    fn las_unsupported_versions_are_errors() {
        for (major, minor) in [(1u8, 5u8), (2, 0), (0, 4)] {
            let mut data = header_bytes(Spec::new(2, 227, 0, 20));
            data[24] = major;
            data[25] = minor;
            let err = read_las_binary(&data).unwrap_err();
            assert!(
                err.contains("unsupported version"),
                "{major}.{minor}: {err}"
            );
        }
    }

    #[test]
    fn las_single_byte_corruption_never_panics_or_over_allocates() {
        // Every header and record byte overwritten with extreme values: the
        // reader may accept or reject the file, but it must neither panic nor
        // abort on an allocation sized from a corrupted count.
        let mut v14 = header_bytes(Spec::new(4, 375, 6, 30).count(0, 2));
        v14.extend(record6([1, 2, 3], 4, 0x11, 0, 2, 0, 0, 0, 1.0));
        v14.extend(record6([4, 5, 6], 4, 0x11, 0, 2, 0, 0, 0, 2.0));
        let mut v12 = header_bytes(Spec::new(2, 227, 3, 34).count(2, 0));
        for raw in [[1, 2, 3], [4, 5, 6]] {
            let mut rec = record0(raw, 1, 2, 3, 4, 5, 6);
            rec.extend_from_slice(&1.5f64.to_le_bytes()); // gps time
            rec.extend_from_slice(&[1, 0, 2, 0, 3, 0]); // rgb
            v12.extend(rec);
        }
        for template in [&v14, &v12] {
            assert!(read_las_binary(template).is_ok());
            for i in 0..template.len() {
                for fill in [0x00u8, 0x01, 0x7F, 0x80, 0xFF] {
                    let mut data = template.clone();
                    data[i] = fill;
                    let _ = read_las_binary(&data);
                }
            }
        }
    }

    #[test]
    fn las_every_proper_prefix_of_a_valid_file_is_rejected() {
        let mut data = header_bytes(Spec::new(4, 375, 6, 30).count(0, 2));
        data.extend(record6([1, 2, 3], 4, 0x11, 0, 2, 0, 0, 0, 1.0));
        data.extend(record6([4, 5, 6], 4, 0x11, 0, 2, 0, 0, 0, 2.0));
        for len in 0..data.len() {
            assert!(
                read_las_binary(&data[..len]).is_err(),
                "prefix of {len} bytes must not parse"
            );
        }
        assert!(read_las_binary(&data).is_ok());
    }
}
