// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! GRIB2 section decoding for [`GribReader`](super::types::GribReader).
//!
//! Octet numbers follow the WMO Manual on Codes (FM 92 GRIB edition 2): they
//! are 1-based and counted from the start of the section, where octets 1-4
//! hold the section length and octet 5 the section number. A
//! [`GribSection`] payload omits those five octets, so octet `n` is
//! `payload[n - 6]`. [`SectionView`] is addressed by WMO octet numbers, which
//! lets every decoder below be checked line by line against the code tables.
//!
//! Supported content:
//!
//! * Section 1 (identification): reference time and originating centre.
//! * Section 3 (grid definition): templates 3.0 (regular latitude/longitude),
//!   3.20 (polar stereographic) and 3.30 (Lambert conformal). Any other
//!   template, a grid with a list of points per row, or a point count that
//!   disagrees with the grid dimensions is an error.
//! * Section 4 (product definition): templates 4.0 to 4.15, which share the
//!   layout of the parameter category/number and the forecast time.
//! * Section 5 (data representation): simple packing (5.0) and raw IEEE
//!   floating point (5.4). Complex packing (5.2, 5.3), JPEG 2000 (5.40), PNG
//!   (5.41) and the other templates are reported as errors instead of being
//!   decoded to garbage.
//! * Section 6 (bit-map): an explicit bit-map, "no bit-map" (255) and "use
//!   the bit-map defined earlier in the message" (254).
//!
//! A message with repeated sections 4-7 holds one field per repetition; every
//! field is decoded on its own.

use super::types::{DataRepresentation, GribMessage, GribSection, GribSectionType, GridDefinition};

/// First payload octet of a section: octets 1-4 are the length and octet 5
/// the section number.
const FIRST_PAYLOAD_OCTET: usize = 6;

/// Highest product definition template that shares the layout of template 4.0
/// up to and including the forecast time (octets 10 to 22).
const LAST_COMMON_PRODUCT_TEMPLATE: u16 = 15;

/// Largest packed-value width the bit reader accepts.
const MAX_PACKED_BITS: u8 = 64;

/// Degrees per unit of a template angle field when no basic angle is given.
const MICRODEGREE: f64 = 1e-6;

/// Metres per unit of a template grid length (millimetres).
const MILLIMETRE: f64 = 1e-3;

/// A section payload addressed by WMO octet numbers.
struct SectionView<'a> {
    number: u8,
    payload: &'a [u8],
}

impl<'a> SectionView<'a> {
    fn new(number: u8, section: &'a GribSection) -> Self {
        Self {
            number,
            payload: &section.payload,
        }
    }

    fn too_short(&self, first_octet: usize, last_octet: usize) -> String {
        format!(
            "GRIB2: section {} ({} payload bytes) is too short for octets {first_octet}-{last_octet}",
            self.number,
            self.payload.len()
        )
    }

    /// `N` bytes starting at the 1-based `octet`.
    fn array<const N: usize>(&self, octet: usize) -> Result<[u8; N], String> {
        let too_short = || self.too_short(octet, octet.saturating_add(N).saturating_sub(1));
        let start = octet
            .checked_sub(FIRST_PAYLOAD_OCTET)
            .ok_or_else(too_short)?;
        let end = start.checked_add(N).ok_or_else(too_short)?;
        let bytes = self.payload.get(start..end).ok_or_else(too_short)?;
        <[u8; N]>::try_from(bytes).map_err(|_| too_short())
    }

    /// Require that the section reaches at least `last_octet`.
    fn require_octets(&self, last_octet: usize) -> Result<(), String> {
        if self.payload.len() + FIRST_PAYLOAD_OCTET > last_octet {
            Ok(())
        } else {
            Err(self.too_short(FIRST_PAYLOAD_OCTET, last_octet))
        }
    }

    fn u8_at(&self, octet: usize) -> Result<u8, String> {
        self.array::<1>(octet).map(|[byte]| byte)
    }

    fn u16_at(&self, octet: usize) -> Result<u16, String> {
        self.array::<2>(octet).map(u16::from_be_bytes)
    }

    fn u32_at(&self, octet: usize) -> Result<u32, String> {
        self.array::<4>(octet).map(u32::from_be_bytes)
    }

    fn f32_at(&self, octet: usize) -> Result<f32, String> {
        self.array::<4>(octet).map(f32::from_be_bytes)
    }
}

/// GRIB2 signed integers are sign-magnitude: the most significant bit is the
/// sign and the remaining bits are the magnitude (not two's complement).
fn sign_magnitude_i32(raw: u32) -> i64 {
    let magnitude = i64::from(raw & 0x7FFF_FFFF);
    if raw & 0x8000_0000 != 0 {
        -magnitude
    } else {
        magnitude
    }
}

/// 16-bit sign-magnitude integer, see [`sign_magnitude_i32`].
fn sign_magnitude_i16(raw: u16) -> i32 {
    let magnitude = i32::from(raw & 0x7FFF);
    if raw & 0x8000 != 0 {
        -magnitude
    } else {
        magnitude
    }
}

/// Degrees of one unit of a template angle field.
///
/// The unit is the ratio of the basic angle and its subdivisions; a zero or
/// missing (all ones) basic angle means 1 and a zero or missing subdivision
/// count means 10^6, i.e. the ordinary unit of 10^-6 degrees.
fn angle_unit_degrees(basic_angle: u32, subdivisions: u32) -> f64 {
    let basic = if basic_angle == 0 || basic_angle == u32::MAX {
        1.0
    } else {
        f64::from(basic_angle)
    };
    let subdivisions = if subdivisions == 0 || subdivisions == u32::MAX {
        1_000_000.0
    } else {
        f64::from(subdivisions)
    };
    basic / subdivisions
}

/// A sign-magnitude angle field in degrees.
fn degrees(raw: u32, unit: f64) -> f64 {
    sign_magnitude_i32(raw) as f64 * unit
}

/// A direction increment: `NaN` when the message does not give it (the flag
/// is clear, or the field holds the all-ones missing value).
fn increment(raw: u32, given: bool, unit: f64) -> f64 {
    if given && raw != u32::MAX {
        f64::from(raw) * unit
    } else {
        f64::NAN
    }
}

// ── Section 1: identification ────────────────────────────────────────────

struct Identification {
    centre: u16,
    /// Seconds since 1970-01-01T00:00:00 UTC.
    reference_time: i64,
}

fn is_leap_year(year: i64) -> bool {
    (year % 4 == 0 && year % 100 != 0) || year % 400 == 0
}

fn days_in_month(year: i64, month: i64) -> i64 {
    match month {
        1 | 3 | 5 | 7 | 8 | 10 | 12 => 31,
        4 | 6 | 9 | 11 => 30,
        _ if is_leap_year(year) => 29,
        _ => 28,
    }
}

/// Days from 1970-01-01 to the given proleptic Gregorian date.
fn days_from_civil(year: i64, month: i64, day: i64) -> i64 {
    let year = if month <= 2 { year - 1 } else { year };
    let era = year.div_euclid(400);
    let year_of_era = year - era * 400;
    let month_from_march = (month + 9) % 12;
    let day_of_year = (153 * month_from_march + 2) / 5 + day - 1;
    let day_of_era = year_of_era * 365 + year_of_era / 4 - year_of_era / 100 + day_of_year;
    era * 146_097 + day_of_era - 719_468
}

fn decode_identification(sec: &SectionView<'_>) -> Result<Identification, String> {
    let centre = sec.u16_at(6)?;
    let year = i64::from(sec.u16_at(13)?);
    let month = i64::from(sec.u8_at(15)?);
    let day = i64::from(sec.u8_at(16)?);
    let hour = i64::from(sec.u8_at(17)?);
    let minute = i64::from(sec.u8_at(18)?);
    let second = i64::from(sec.u8_at(19)?);
    let valid = (1..=12).contains(&month)
        && (1..=days_in_month(year, month)).contains(&day)
        && hour <= 23
        && minute <= 59
        && second <= 60;
    if !valid {
        return Err(format!(
            "GRIB2: invalid reference time {year:04}-{month:02}-{day:02} {hour:02}:{minute:02}:{second:02}"
        ));
    }
    Ok(Identification {
        centre,
        reference_time: days_from_civil(year, month, day) * 86_400
            + hour * 3_600
            + minute * 60
            + second,
    })
}

// ── Section 3: grid definition ───────────────────────────────────────────

struct GridInfo {
    definition: GridDefinition,
    /// Number of data points of the grid (octets 7-10).
    n_points: usize,
}

/// Template 3.0, regular latitude/longitude grid (72 octets).
fn decode_lat_lon(sec: &SectionView<'_>) -> Result<GridDefinition, String> {
    sec.require_octets(72)?;
    let ni = sec.u32_at(31)?;
    let nj = sec.u32_at(35)?;
    if ni == u32::MAX || nj == u32::MAX {
        return Err(
            "GRIB2: lat/lon grid without Ni/Nj (quasi-regular grid) is not supported".to_string(),
        );
    }
    let unit = angle_unit_degrees(sec.u32_at(39)?, sec.u32_at(43)?);
    let lat_first = degrees(sec.u32_at(47)?, unit);
    let lon_first = degrees(sec.u32_at(51)?, unit);
    let flags = sec.u8_at(55)?;
    let lat_last = degrees(sec.u32_at(56)?, unit);
    let lon_last = degrees(sec.u32_at(60)?, unit);
    // Flag table 3.3: bit 3 (0x20) = i direction increments given, bit 4
    // (0x10) = j direction increments given.
    let di = increment(sec.u32_at(64)?, flags & 0x20 != 0, unit);
    let dj = increment(sec.u32_at(68)?, flags & 0x10 != 0, unit);
    Ok(GridDefinition::LatLon {
        ni,
        nj,
        lat_first,
        lat_last,
        lon_first,
        lon_last,
        di,
        dj,
    })
}

/// Template 3.20, polar stereographic projection (65 octets).
fn decode_polar_stereographic(sec: &SectionView<'_>) -> Result<GridDefinition, String> {
    sec.require_octets(65)?;
    Ok(GridDefinition::PolarStereographic {
        nx: sec.u32_at(31)?,
        ny: sec.u32_at(35)?,
        lat_first: degrees(sec.u32_at(39)?, MICRODEGREE),
        lon_first: degrees(sec.u32_at(43)?, MICRODEGREE),
        dx: f64::from(sec.u32_at(56)?) * MILLIMETRE,
        dy: f64::from(sec.u32_at(60)?) * MILLIMETRE,
    })
}

/// Template 3.30, Lambert conformal projection (81 octets).
fn decode_lambert_conformal(sec: &SectionView<'_>) -> Result<GridDefinition, String> {
    sec.require_octets(81)?;
    Ok(GridDefinition::LambertConformal {
        nx: sec.u32_at(31)?,
        ny: sec.u32_at(35)?,
        lat_ref: degrees(sec.u32_at(48)?, MICRODEGREE),
        lon_ref: degrees(sec.u32_at(52)?, MICRODEGREE),
        latin1: degrees(sec.u32_at(66)?, MICRODEGREE),
        latin2: degrees(sec.u32_at(70)?, MICRODEGREE),
        dx: f64::from(sec.u32_at(56)?) * MILLIMETRE,
        dy: f64::from(sec.u32_at(60)?) * MILLIMETRE,
    })
}

fn decode_grid(sec: &SectionView<'_>) -> Result<GridInfo, String> {
    let declared_points = sec.u32_at(7)?;
    let list_octets = sec.u8_at(11)?;
    let template = sec.u16_at(13)?;
    if list_octets != 0 {
        return Err(
            "GRIB2: grids with a list of points per row (quasi-regular grids) are not supported"
                .to_string(),
        );
    }
    let definition = match template {
        0 => decode_lat_lon(sec)?,
        20 => decode_polar_stereographic(sec)?,
        30 => decode_lambert_conformal(sec)?,
        other => {
            return Err(format!(
                "GRIB2: unsupported grid definition template 3.{other} (supported: 3.0, 3.20, 3.30)"
            ));
        }
    };
    let (first, second) = match &definition {
        GridDefinition::LatLon { ni, nj, .. } => (*ni, *nj),
        GridDefinition::LambertConformal { nx, ny, .. }
        | GridDefinition::PolarStereographic { nx, ny, .. } => (*nx, *ny),
    };
    if u64::from(first) * u64::from(second) != u64::from(declared_points) {
        return Err(format!(
            "GRIB2: section 3 declares {declared_points} data points but the grid is {first} x {second}"
        ));
    }
    let n_points = usize::try_from(declared_points).map_err(|_| {
        format!("GRIB2: grid of {declared_points} points does not fit this platform")
    })?;
    Ok(GridInfo {
        definition,
        n_points,
    })
}

// ── Section 4: product definition ────────────────────────────────────────

struct ProductInfo {
    category: u8,
    number: u8,
    /// Forecast time in seconds.
    forecast_seconds: i64,
}

/// Seconds per unit of Code Table 4.4. `None` is the missing value (255),
/// which carries no forecast time; units with no fixed length (month, year,
/// decade, normal, century) and reserved codes are errors.
fn forecast_unit_seconds(unit: u8) -> Result<Option<i64>, String> {
    match unit {
        0 => Ok(Some(60)),
        1 => Ok(Some(3_600)),
        2 => Ok(Some(86_400)),
        10 => Ok(Some(10_800)),
        11 => Ok(Some(21_600)),
        12 => Ok(Some(43_200)),
        13 => Ok(Some(1)),
        255 => Ok(None),
        other => Err(format!(
            "GRIB2: unsupported forecast time unit {other} (Code Table 4.4); calendar units have no fixed length in seconds"
        )),
    }
}

/// Templates 4.0-4.15: parameter category (octet 10), parameter number
/// (octet 11), time unit (octet 18) and forecast time (octets 19-22).
fn decode_product(sec: &SectionView<'_>) -> Result<ProductInfo, String> {
    let template = sec.u16_at(8)?;
    if template > LAST_COMMON_PRODUCT_TEMPLATE {
        return Err(format!(
            "GRIB2: unsupported product definition template 4.{template} (supported: 4.0 to 4.{LAST_COMMON_PRODUCT_TEMPLATE})"
        ));
    }
    sec.require_octets(22)?;
    let category = sec.u8_at(10)?;
    let number = sec.u8_at(11)?;
    let unit = sec.u8_at(18)?;
    let value = i64::from(sec.u32_at(19)?);
    let forecast_seconds = match forecast_unit_seconds(unit)? {
        Some(seconds_per_unit) => value * seconds_per_unit,
        None => 0,
    };
    Ok(ProductInfo {
        category,
        number,
        forecast_seconds,
    })
}

// ── Section 5: data representation ───────────────────────────────────────

/// How the values of a Section 7 are packed.
enum Packing {
    /// Template 5.0: `Y = (R + X * 2^E) / 10^D` with `nbits` bits per `X`.
    Simple {
        reference: f32,
        binary_scale: i32,
        decimal_scale: i32,
        bits: u8,
    },
    /// Template 5.4: raw IEEE floats (32-bit, or 64-bit when `double`).
    Ieee { double: bool },
}

struct DataRepresentationInfo {
    /// Number of values stored in Section 7: the grid points, or the points a
    /// bit-map marks as present.
    n_values: usize,
    representation: DataRepresentation,
    packing: Packing,
}

fn unsupported_template_name(template: u16) -> Option<&'static str> {
    match template {
        2 => Some("complex packing"),
        3 => Some("complex packing with spatial differencing"),
        40 => Some("JPEG 2000"),
        41 => Some("PNG"),
        42 => Some("CCSDS recommended lossless compression"),
        200 => Some("run length packing"),
        _ => None,
    }
}

fn decode_data_representation(sec: &SectionView<'_>) -> Result<DataRepresentationInfo, String> {
    let declared_values = sec.u32_at(6)?;
    let n_values = usize::try_from(declared_values)
        .map_err(|_| format!("GRIB2: {declared_values} data values do not fit this platform"))?;
    let template = sec.u16_at(10)?;
    match template {
        0 => {
            sec.require_octets(20)?;
            let bits = sec.u8_at(20)?;
            if bits > MAX_PACKED_BITS {
                return Err(format!(
                    "GRIB2: {bits} bits per packed value exceeds the supported {MAX_PACKED_BITS}"
                ));
            }
            Ok(DataRepresentationInfo {
                n_values,
                representation: DataRepresentation::SimplePacking,
                packing: Packing::Simple {
                    reference: sec.f32_at(12)?,
                    binary_scale: sign_magnitude_i16(sec.u16_at(16)?),
                    decimal_scale: sign_magnitude_i16(sec.u16_at(18)?),
                    bits,
                },
            })
        }
        4 => {
            let double = match sec.u8_at(12)? {
                1 => false,
                2 => true,
                3 => {
                    return Err(
                        "GRIB2: 128-bit IEEE floating point data (template 5.4) is not supported"
                            .to_string(),
                    );
                }
                other => {
                    return Err(format!(
                        "GRIB2: unknown IEEE precision {other} in data representation template 5.4"
                    ));
                }
            };
            Ok(DataRepresentationInfo {
                n_values,
                representation: DataRepresentation::IeeeFloat,
                packing: Packing::Ieee { double },
            })
        }
        other => match unsupported_template_name(other) {
            Some(name) => Err(format!(
                "GRIB2: data representation template 5.{other} ({name}) is not supported; only simple packing (5.0) and IEEE floating point (5.4) are decoded"
            )),
            None => Err(format!(
                "GRIB2: unknown data representation template 5.{other}"
            )),
        },
    }
}

// ── Section 7: data ──────────────────────────────────────────────────────

/// MSB-first reader of packed unsigned integers of up to 64 bits.
struct BitReader<'a> {
    data: &'a [u8],
    bit_pos: usize,
}

impl<'a> BitReader<'a> {
    fn new(data: &'a [u8]) -> Self {
        Self { data, bit_pos: 0 }
    }

    /// The next `width` bits (`width <= 64`), or `None` past the end.
    fn read(&mut self, width: usize) -> Option<u64> {
        let mut value = 0u64;
        let mut remaining = width;
        while remaining > 0 {
            let byte = u64::from(*self.data.get(self.bit_pos / 8)?);
            let available = 8 - self.bit_pos % 8;
            let take = remaining.min(available);
            let bits = (byte >> (available - take)) & ((1u64 << take) - 1);
            value = (value << take) | bits;
            self.bit_pos += take;
            remaining -= take;
        }
        Some(value)
    }
}

/// Bit `index` (0 = most significant bit of the first byte) of `bytes`.
fn bit_at(bytes: &[u8], index: usize) -> bool {
    bytes
        .get(index / 8)
        .is_some_and(|byte| (byte >> (7 - index % 8)) & 1 != 0)
}

/// Unpack `n` simple-packed values: `Y = (R + X * 2^E) / 10^D`.
fn unpack_simple(
    data: &[u8],
    n: usize,
    reference: f32,
    binary_scale: i32,
    decimal_scale: i32,
    bits: u8,
) -> Result<Vec<f64>, String> {
    let reference = f64::from(reference);
    let binary = 2f64.powi(binary_scale);
    let decimal = 10f64.powi(decimal_scale.abs());
    let decode = |x: f64| {
        let scaled = reference + x * binary;
        if decimal_scale >= 0 {
            scaled / decimal
        } else {
            scaled * decimal
        }
    };
    if bits == 0 {
        // A constant field: every value is the reference value.
        return Ok(vec![decode(0.0); n]);
    }
    let needed_bits = (n as u128) * u128::from(bits);
    if needed_bits > (data.len() as u128) * 8 {
        return Err(format!(
            "GRIB2: data section holds {} bytes but {n} values of {bits} bits need {} bytes",
            data.len(),
            needed_bits.div_ceil(8)
        ));
    }
    let mut reader = BitReader::new(data);
    let mut values = Vec::with_capacity(n);
    for _ in 0..n {
        let packed = reader
            .read(usize::from(bits))
            .ok_or_else(|| "GRIB2: data section ended inside a packed value".to_string())?;
        values.push(decode(packed as f64));
    }
    Ok(values)
}

/// Unpack `n` raw big-endian IEEE values (template 5.4).
fn unpack_ieee(data: &[u8], n: usize, double: bool) -> Result<Vec<f64>, String> {
    let width = if double { 8 } else { 4 };
    let needed = n
        .checked_mul(width)
        .ok_or_else(|| format!("GRIB2: {n} IEEE values overflow the data section size"))?;
    if data.len() < needed {
        return Err(format!(
            "GRIB2: data section holds {} bytes but {n} IEEE values of {width} bytes need {needed}",
            data.len()
        ));
    }
    let truncated = |_| "GRIB2: data section ended inside an IEEE value".to_string();
    let mut values = Vec::with_capacity(n);
    for chunk in data.chunks_exact(width).take(n) {
        let value = if double {
            f64::from_be_bytes(<[u8; 8]>::try_from(chunk).map_err(truncated)?)
        } else {
            f64::from(f32::from_be_bytes(
                <[u8; 4]>::try_from(chunk).map_err(truncated)?,
            ))
        };
        values.push(value);
    }
    Ok(values)
}

/// Decode a Section 7 payload into one value per grid point (`NaN` where the
/// bit-map marks a point missing) and the bit-map flags, when there is one.
fn decode_field_values(
    grid: &GridInfo,
    representation: &DataRepresentationInfo,
    bitmap: Option<&[u8]>,
    data: &[u8],
    max_points: usize,
) -> Result<(Vec<f64>, Option<Vec<bool>>), String> {
    if grid.n_points > max_points {
        return Err(format!(
            "GRIB2: grid of {} points exceeds the reader limit of {max_points} points (see GribReader::with_max_points)",
            grid.n_points
        ));
    }
    let mask = match bitmap {
        Some(bytes) => {
            let available = bytes.len().saturating_mul(8);
            if available < grid.n_points {
                return Err(format!(
                    "GRIB2: bit-map holds {available} bits but the grid has {} points",
                    grid.n_points
                ));
            }
            let mask: Vec<bool> = (0..grid.n_points).map(|i| bit_at(bytes, i)).collect();
            let present = mask.iter().filter(|&&flag| flag).count();
            if present != representation.n_values {
                return Err(format!(
                    "GRIB2: bit-map marks {present} points present but section 5 declares {} values",
                    representation.n_values
                ));
            }
            Some(mask)
        }
        None => {
            if representation.n_values != grid.n_points {
                return Err(format!(
                    "GRIB2: section 5 declares {} values but the grid has {} points and no bit-map applies",
                    representation.n_values, grid.n_points
                ));
            }
            None
        }
    };
    let values = match representation.packing {
        Packing::Simple {
            reference,
            binary_scale,
            decimal_scale,
            bits,
        } => unpack_simple(
            data,
            representation.n_values,
            reference,
            binary_scale,
            decimal_scale,
            bits,
        )?,
        Packing::Ieee { double } => unpack_ieee(data, representation.n_values, double)?,
    };
    let expanded = match &mask {
        Some(mask) => {
            let mut packed = values.into_iter();
            mask.iter()
                .map(|&present| {
                    if present {
                        packed.next().unwrap_or(f64::NAN)
                    } else {
                        f64::NAN
                    }
                })
                .collect()
        }
        None => values,
    };
    Ok((expanded, mask))
}

// ── Message assembly ─────────────────────────────────────────────────────

/// Whether the current field has seen its Section 6.
#[derive(Clone, Copy, Default)]
enum BitmapState {
    /// No Section 6 since the last data section.
    #[default]
    Undefined,
    /// Section 6 said that no bit-map applies (indicator 255).
    Absent,
    /// A bit-map applies (explicit, or reused with indicator 254).
    Present,
}

/// Sections collected while walking a message; every Section 7 completes one
/// field.
#[derive(Default)]
struct FieldState {
    identification: Option<Identification>,
    grid: Option<GridInfo>,
    product: Option<ProductInfo>,
    representation: Option<DataRepresentationInfo>,
    /// The last explicit bit-map of the message (valid for the current grid).
    bitmap: Option<Vec<u8>>,
    bitmap_state: BitmapState,
}

impl FieldState {
    fn apply(
        &mut self,
        section: &GribSection,
        discipline: u8,
        max_points: usize,
    ) -> Result<Option<GribMessage>, String> {
        match section.section_type {
            GribSectionType::Indicator | GribSectionType::LocalUse | GribSectionType::End => {}
            GribSectionType::Identification => {
                self.identification = Some(decode_identification(&SectionView::new(1, section))?);
            }
            GribSectionType::GridDefinitionSection => {
                self.grid = Some(decode_grid(&SectionView::new(3, section))?);
                // A bit-map belongs to the grid it was defined for.
                self.bitmap = None;
            }
            GribSectionType::ProductDefinition => {
                self.product = Some(decode_product(&SectionView::new(4, section))?);
            }
            GribSectionType::DataRepresentationSection => {
                self.representation =
                    Some(decode_data_representation(&SectionView::new(5, section))?);
            }
            GribSectionType::Bitmap => {
                let view = SectionView::new(6, section);
                match view.u8_at(6)? {
                    0 => {
                        self.bitmap = Some(section.payload.get(1..).unwrap_or_default().to_vec());
                        self.bitmap_state = BitmapState::Present;
                    }
                    254 => {
                        if self.bitmap.is_none() {
                            return Err("GRIB2: bit-map indicator 254 refers to a bit-map defined earlier in the message, but none was".to_string());
                        }
                        self.bitmap_state = BitmapState::Present;
                    }
                    255 => self.bitmap_state = BitmapState::Absent,
                    other => {
                        return Err(format!(
                            "GRIB2: pre-defined bit-map {other} is not supported"
                        ));
                    }
                }
            }
            GribSectionType::Data => return self.complete_field(section, discipline, max_points),
        }
        Ok(None)
    }

    /// Decode the field that the data section `section` completes.
    fn complete_field(
        &mut self,
        section: &GribSection,
        discipline: u8,
        max_points: usize,
    ) -> Result<Option<GribMessage>, String> {
        let missing = |name: &str| format!("GRIB2: data section without {name}");
        let identification = self
            .identification
            .as_ref()
            .ok_or_else(|| missing("an identification section (section 1)"))?;
        let grid = self
            .grid
            .as_ref()
            .ok_or_else(|| missing("a grid definition section (section 3)"))?;
        let product = self
            .product
            .as_ref()
            .ok_or_else(|| missing("a product definition section (section 4)"))?;
        let representation = self
            .representation
            .as_ref()
            .ok_or_else(|| missing("a data representation section (section 5)"))?;
        let bitmap = match self.bitmap_state {
            BitmapState::Undefined => {
                return Err(missing("a bit-map section (section 6)"));
            }
            BitmapState::Absent => None,
            BitmapState::Present => self.bitmap.as_deref(),
        };
        let (data, bitmap_flags) =
            decode_field_values(grid, representation, bitmap, &section.payload, max_points)?;
        let message = GribMessage {
            discipline,
            parameter_category: product.category,
            parameter_number: product.number,
            grid_definition: grid.definition.clone(),
            data_representation: representation.representation,
            reference_time: identification.reference_time,
            forecast_time: product.forecast_seconds,
            centre_id: identification.centre,
            data,
            bitmap: bitmap_flags,
        };
        // The next field needs its own sections 4-6; the grid is kept.
        self.product = None;
        self.representation = None;
        self.bitmap_state = BitmapState::Undefined;
        Ok(Some(message))
    }
}

/// Decode every field of one message from its parsed sections.
///
/// Fields that were decoded before a section fails are still returned, in
/// order, followed by that error; a message without any data section yields
/// a single error.
pub(super) fn assemble_fields(
    sections: &[GribSection],
    discipline: u8,
    max_points: usize,
) -> Vec<Result<GribMessage, String>> {
    let mut fields = Vec::new();
    let mut state = FieldState::default();
    for section in sections {
        match state.apply(section, discipline, max_points) {
            Ok(Some(message)) => fields.push(Ok(message)),
            Ok(None) => {}
            Err(error) => {
                fields.push(Err(error));
                return fields;
            }
        }
    }
    if fields.is_empty() {
        fields.push(Err(
            "GRIB2: message contains no data section (section 7)".to_string()
        ));
    }
    fields
}
