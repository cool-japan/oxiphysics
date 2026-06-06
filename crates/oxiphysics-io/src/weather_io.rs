// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Weather and atmospheric data I/O.
//!
//! This module provides parsers, readers, and models for common weather and
//! atmospheric data formats and computations:
//!
//! - **GRIB2 basic reader**: Parse GRIB2 binary format headers and extract metadata.
//! - **NetCDF atmospheric**: Read/write simple NetCDF-like atmospheric data.
//! - **METAR parser**: Parse METAR aviation weather reports.
//! - **Sounding (TEMP) data**: Parse upper-air radiosonde sounding data.
//! - **Wind profile**: Vertical wind profile models (power law, log law).
//! - **Temperature profile**: Vertical temperature profile models.
//! - **Pressure altitude**: Barometric altitude calculations.
//! - **ISA atmosphere model**: International Standard Atmosphere.
//! - **Humidity/dew point conversion**: Various humidity conversions.

use std::f64::consts::E;

// ════════════════════════════════════════════════════════════════════════════
// GRIB2 Basic Reader
// ════════════════════════════════════════════════════════════════════════════

/// GRIB2 section identifiers.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Grib2Section {
    /// Section 0: Indicator section.
    Indicator,
    /// Section 1: Identification section.
    Identification,
    /// Section 2: Local use section.
    LocalUse,
    /// Section 3: Grid definition section.
    GridDefinition,
    /// Section 4: Product definition section.
    ProductDefinition,
    /// Section 5: Data representation section.
    DataRepresentation,
    /// Section 6: Bit-map section.
    BitMap,
    /// Section 7: Data section.
    Data,
    /// Section 8: End section ("7777").
    End,
}

/// Metadata extracted from a GRIB2 file header.
#[derive(Debug, Clone)]
pub struct Grib2Header {
    /// GRIB edition number (should be 2).
    pub edition: u8,
    /// Discipline (0=Meteorological, 1=Hydrological, 2=Land surface, etc.).
    pub discipline: u8,
    /// Total message length in bytes.
    pub total_length: u64,
    /// Originating center ID.
    pub center_id: u16,
    /// Originating sub-center ID.
    pub subcenter_id: u16,
    /// Reference time (year, month, day, hour, minute, second).
    pub reference_time: [u16; 6],
    /// Production status (0=Operational, 1=Test, etc.).
    pub production_status: u8,
    /// Type of data (0=Analysis, 1=Forecast, etc.).
    pub data_type: u8,
}

/// Parse the indicator section (Section 0) of a GRIB2 message.
///
/// The indicator section is always 16 bytes:
/// - Bytes 0-3: "GRIB"
/// - Bytes 4-5: Reserved
/// - Byte 6: Discipline
/// - Byte 7: Edition (must be 2)
/// - Bytes 8-15: Total message length
///
/// Returns `None` if the data is too short or doesn't start with "GRIB".
pub fn parse_grib2_indicator(data: &[u8]) -> Option<(u8, u8, u64)> {
    if data.len() < 16 {
        return None;
    }
    // Check magic bytes "GRIB"
    if data[0] != b'G' || data[1] != b'R' || data[2] != b'I' || data[3] != b'B' {
        return None;
    }
    let discipline = data[6];
    let edition = data[7];
    if edition != 2 {
        return None;
    }
    let total_length = u64::from_be_bytes([
        data[8], data[9], data[10], data[11], data[12], data[13], data[14], data[15],
    ]);
    Some((discipline, edition, total_length))
}

/// Parse the identification section (Section 1) of a GRIB2 message.
///
/// Section 1 starts after the 16-byte indicator.
/// Returns a `Grib2Header` on success.
pub fn parse_grib2_header(data: &[u8]) -> Option<Grib2Header> {
    let (discipline, edition, total_length) = parse_grib2_indicator(data)?;

    if data.len() < 16 + 21 {
        return None;
    }

    let sec1 = &data[16..];
    // Section 1 length
    let _sec1_len = u32::from_be_bytes([sec1[0], sec1[1], sec1[2], sec1[3]]);
    // Section number (should be 1)
    if sec1[4] != 1 {
        return None;
    }

    let center_id = u16::from_be_bytes([sec1[5], sec1[6]]);
    let subcenter_id = u16::from_be_bytes([sec1[7], sec1[8]]);
    // Master tables version, local tables version...
    // Reference time
    let year = u16::from_be_bytes([sec1[12], sec1[13]]);
    let month = sec1[14] as u16;
    let day = sec1[15] as u16;
    let hour = sec1[16] as u16;
    let minute = sec1[17] as u16;
    let second = sec1[18] as u16;
    let production_status = sec1[19];
    let data_type = sec1[20];

    Some(Grib2Header {
        edition,
        discipline,
        total_length,
        center_id,
        subcenter_id,
        reference_time: [year, month, day, hour, minute, second],
        production_status,
        data_type,
    })
}

/// Check if raw bytes contain a valid GRIB2 end section marker ("7777").
pub fn has_grib2_end_marker(data: &[u8]) -> bool {
    if data.len() < 4 {
        return false;
    }
    let end = &data[data.len() - 4..];
    end == b"7777"
}

/// Get the GRIB2 discipline name.
pub fn grib2_discipline_name(discipline: u8) -> &'static str {
    match discipline {
        0 => "Meteorological",
        1 => "Hydrological",
        2 => "Land Surface",
        3 => "Satellite Remote Sensing",
        4 => "Space Weather",
        10 => "Oceanographic",
        _ => "Unknown",
    }
}

// ════════════════════════════════════════════════════════════════════════════
// NetCDF Atmospheric Data
// ════════════════════════════════════════════════════════════════════════════

/// A simple representation of a NetCDF atmospheric variable.
#[derive(Debug, Clone)]
pub struct NetcdfAtmVariable {
    /// Variable name (e.g., "temperature", "wind_speed").
    pub name: String,
    /// Units (e.g., "K", "m/s", "Pa").
    pub units: String,
    /// Dimension names.
    pub dimensions: Vec<String>,
    /// Data values (flattened).
    pub data: Vec<f64>,
    /// Shape of each dimension.
    pub shape: Vec<usize>,
    /// Fill/missing value.
    pub fill_value: f64,
    /// Scale factor for packed data.
    pub scale_factor: f64,
    /// Add offset for packed data.
    pub add_offset: f64,
}

impl NetcdfAtmVariable {
    /// Create a new atmospheric variable.
    pub fn new(name: &str, units: &str, dimensions: Vec<String>, shape: Vec<usize>) -> Self {
        let total_size: usize = shape.iter().product();
        Self {
            name: name.to_string(),
            units: units.to_string(),
            dimensions,
            data: vec![0.0; total_size],
            shape,
            fill_value: -9999.0,
            scale_factor: 1.0,
            add_offset: 0.0,
        }
    }

    /// Unpack a raw value using scale_factor and add_offset.
    pub fn unpack(&self, raw: f64) -> f64 {
        raw * self.scale_factor + self.add_offset
    }

    /// Pack a physical value for storage.
    pub fn pack(&self, value: f64) -> f64 {
        (value - self.add_offset) / self.scale_factor
    }

    /// Get a value at a flat index.
    pub fn get_value(&self, index: usize) -> Option<f64> {
        self.data.get(index).copied().map(|v| {
            if (v - self.fill_value).abs() < 1e-10 {
                f64::NAN
            } else {
                self.unpack(v)
            }
        })
    }

    /// Set a value at a flat index.
    pub fn set_value(&mut self, index: usize, value: f64) {
        if index < self.data.len() {
            self.data[index] = self.pack(value);
        }
    }

    /// Get the total number of elements.
    pub fn total_size(&self) -> usize {
        self.shape.iter().product()
    }

    /// Convert multi-dimensional indices to flat index.
    pub fn flat_index(&self, indices: &[usize]) -> Option<usize> {
        if indices.len() != self.shape.len() {
            return None;
        }
        let mut idx = 0;
        let mut stride = 1;
        for i in (0..indices.len()).rev() {
            if indices[i] >= self.shape[i] {
                return None;
            }
            idx += indices[i] * stride;
            stride *= self.shape[i];
        }
        Some(idx)
    }
}

/// A simple NetCDF-like atmospheric dataset.
#[derive(Debug, Clone)]
pub struct NetcdfAtmDataset {
    /// Global attributes (key-value pairs).
    pub global_attrs: Vec<(String, String)>,
    /// Variables in the dataset.
    pub variables: Vec<NetcdfAtmVariable>,
}

impl NetcdfAtmDataset {
    /// Create an empty dataset.
    pub fn new() -> Self {
        Self {
            global_attrs: Vec::new(),
            variables: Vec::new(),
        }
    }

    /// Add a global attribute.
    pub fn add_attribute(&mut self, key: &str, value: &str) {
        self.global_attrs.push((key.to_string(), value.to_string()));
    }

    /// Add a variable.
    pub fn add_variable(&mut self, var: NetcdfAtmVariable) {
        self.variables.push(var);
    }

    /// Find a variable by name.
    pub fn find_variable(&self, name: &str) -> Option<&NetcdfAtmVariable> {
        self.variables.iter().find(|v| v.name == name)
    }

    /// Find a mutable variable by name.
    pub fn find_variable_mut(&mut self, name: &str) -> Option<&mut NetcdfAtmVariable> {
        self.variables.iter_mut().find(|v| v.name == name)
    }

    /// Serialize the dataset to a simple text format for exchange.
    ///
    /// Format:
    /// ```text
    /// # NETCDF_ATM_SIMPLE v1
    /// @attr key=value
    /// !var name units dim1,dim2 shape1,shape2
    /// data_val1 data_val2 ...
    /// ```
    pub fn to_text(&self) -> String {
        let mut out = String::new();
        out.push_str("# NETCDF_ATM_SIMPLE v1\n");

        for (k, v) in &self.global_attrs {
            out.push_str(&format!("@attr {}={}\n", k, v));
        }

        for var in &self.variables {
            let dims = var.dimensions.join(",");
            let shape: Vec<String> = var.shape.iter().map(|s| s.to_string()).collect();
            let shape_str = shape.join(",");
            out.push_str(&format!(
                "!var {} {} {} {}\n",
                var.name, var.units, dims, shape_str
            ));

            let vals: Vec<String> = var.data.iter().map(|v| format!("{:.6}", v)).collect();
            // Write in chunks of 10 values per line
            for chunk in vals.chunks(10) {
                out.push_str(&chunk.join(" "));
                out.push('\n');
            }
        }
        out
    }

    /// Parse a dataset from the simple text format.
    pub fn from_text(text: &str) -> Option<Self> {
        let mut dataset = Self::new();
        let mut lines = text.lines().peekable();

        // Check header
        let header = lines.next()?;
        if !header.starts_with("# NETCDF_ATM_SIMPLE") {
            return None;
        }

        while let Some(line) = lines.next() {
            let line = line.trim();
            if line.is_empty() || line.starts_with('#') {
                continue;
            }

            if let Some(attr) = line.strip_prefix("@attr ") {
                if let Some((k, v)) = attr.split_once('=') {
                    dataset.add_attribute(k, v);
                }
            } else if let Some(var_def) = line.strip_prefix("!var ") {
                let parts: Vec<&str> = var_def.split_whitespace().collect();
                if parts.len() < 4 {
                    continue;
                }
                let name = parts[0];
                let units = parts[1];
                let dims: Vec<String> = parts[2].split(',').map(|s| s.to_string()).collect();
                let shape: Vec<usize> =
                    parts[3].split(',').filter_map(|s| s.parse().ok()).collect();

                let total: usize = shape.iter().product();
                let mut var = NetcdfAtmVariable::new(name, units, dims, shape);

                // Read data values
                let mut values = Vec::with_capacity(total);
                while values.len() < total {
                    if let Some(data_line) = lines.peek() {
                        let data_line = data_line.trim();
                        if data_line.starts_with('@')
                            || data_line.starts_with('!')
                            || data_line.starts_with('#')
                        {
                            break;
                        }
                        let parsed: Vec<f64> = data_line
                            .split_whitespace()
                            .filter_map(|s| s.parse().ok())
                            .collect();
                        values.extend(parsed);
                        lines.next();
                    } else {
                        break;
                    }
                }
                var.data = values;
                dataset.add_variable(var);
            }
        }
        Some(dataset)
    }
}

impl Default for NetcdfAtmDataset {
    fn default() -> Self {
        Self::new()
    }
}

// ════════════════════════════════════════════════════════════════════════════
// METAR Parser
// ════════════════════════════════════════════════════════════════════════════

/// Parsed METAR weather observation.
#[derive(Debug, Clone)]
pub struct MetarReport {
    /// Station identifier (ICAO code, e.g., "KJFK").
    pub station: String,
    /// Observation day of month.
    pub day: u8,
    /// Observation hour (UTC).
    pub hour: u8,
    /// Observation minute (UTC).
    pub minute: u8,
    /// Wind direction (degrees, 0-360). `None` if variable.
    pub wind_direction: Option<u16>,
    /// Wind speed (knots).
    pub wind_speed: u16,
    /// Wind gust speed (knots). `None` if no gusts.
    pub wind_gust: Option<u16>,
    /// Visibility in statute miles.
    pub visibility_sm: f64,
    /// Temperature (Celsius).
    pub temperature: Option<f64>,
    /// Dew point (Celsius).
    pub dew_point: Option<f64>,
    /// Altimeter setting (inches of mercury).
    pub altimeter: Option<f64>,
    /// Raw METAR string.
    pub raw: String,
    /// Cloud layers.
    pub clouds: Vec<CloudLayer>,
    /// Weather phenomena (e.g., "RA", "SN", "FG").
    pub weather: Vec<String>,
    /// Whether this is an automatic observation.
    pub auto: bool,
}

/// A cloud layer reported in METAR.
#[derive(Debug, Clone)]
pub struct CloudLayer {
    /// Coverage code: "FEW", "SCT", "BKN", "OVC", "CLR", "SKC".
    pub coverage: String,
    /// Height in hundreds of feet AGL. `None` for CLR/SKC.
    pub height_ft: Option<u32>,
    /// Cloud type modifier (e.g., "CB", "TCU"). Usually empty.
    pub cloud_type: String,
}

impl MetarReport {
    /// Create a new empty METAR report.
    pub fn new(station: &str) -> Self {
        Self {
            station: station.to_string(),
            day: 0,
            hour: 0,
            minute: 0,
            wind_direction: None,
            wind_speed: 0,
            wind_gust: None,
            visibility_sm: 10.0,
            temperature: None,
            dew_point: None,
            altimeter: None,
            raw: String::new(),
            clouds: Vec::new(),
            weather: Vec::new(),
            auto: false,
        }
    }
}

/// Parse a METAR string into a `MetarReport`.
///
/// Handles the most common METAR elements:
/// - Station ID, date/time group
/// - Wind (direction, speed, gusts)
/// - Visibility
/// - Temperature/dew point
/// - Altimeter setting
/// - Cloud layers
///
/// # Example
/// ```text
/// KJFK 301456Z 32015G25KT 10SM FEW050 SCT250 24/16 A3002
/// ```
pub fn parse_metar(raw: &str) -> Option<MetarReport> {
    let raw_trimmed = raw.trim();
    // Remove "METAR" or "SPECI" prefix if present
    let text = raw_trimmed
        .strip_prefix("METAR ")
        .or_else(|| raw_trimmed.strip_prefix("SPECI "))
        .unwrap_or(raw_trimmed);

    let tokens: Vec<&str> = text.split_whitespace().collect();
    if tokens.len() < 3 {
        return None;
    }

    let mut report = MetarReport::new(tokens[0]);
    report.raw = raw_trimmed.to_string();

    let mut idx = 1;

    // Parse date/time group (DDHHMMz)
    if idx < tokens.len() && tokens[idx].len() == 7 && tokens[idx].ends_with('Z') {
        let dtg = tokens[idx];
        report.day = dtg[0..2].parse().unwrap_or(0);
        report.hour = dtg[2..4].parse().unwrap_or(0);
        report.minute = dtg[4..6].parse().unwrap_or(0);
        idx += 1;
    }

    // Check for AUTO
    if idx < tokens.len() && tokens[idx] == "AUTO" {
        report.auto = true;
        idx += 1;
    }

    // Parse remaining tokens
    while idx < tokens.len() {
        let token = tokens[idx];
        idx += 1;

        // Wind: DDDSSKT or DDDSSGSKT
        if token.ends_with("KT") && token.len() >= 7 {
            parse_metar_wind(token, &mut report);
            continue;
        }

        // Visibility: XXXXSM
        if token.ends_with("SM") {
            if let Some(vis_str) = token.strip_suffix("SM") {
                if let Ok(vis) = vis_str.parse::<f64>() {
                    report.visibility_sm = vis;
                } else if vis_str.contains('/') {
                    // Handle fractions like "1/2SM" or "1 1/2SM"
                    let parts: Vec<&str> = vis_str.split('/').collect();
                    if parts.len() == 2 {
                        let num: f64 = parts[0].parse().unwrap_or(0.0);
                        let den: f64 = parts[1].parse().unwrap_or(1.0);
                        if den.abs() > 1e-10 {
                            report.visibility_sm = num / den;
                        }
                    }
                }
            }
            continue;
        }

        // Cloud layers: FEW/SCT/BKN/OVC + height, or CLR/SKC
        if token.starts_with("FEW")
            || token.starts_with("SCT")
            || token.starts_with("BKN")
            || token.starts_with("OVC")
        {
            let coverage = &token[..3];
            let height_str = &token[3..];
            let mut cloud_type = String::new();

            // Check for CB/TCU suffix
            let height_part = if let Some(stripped) = height_str.strip_suffix("CB") {
                cloud_type = "CB".to_string();
                stripped
            } else if let Some(stripped) = height_str.strip_suffix("TCU") {
                cloud_type = "TCU".to_string();
                stripped
            } else {
                height_str
            };

            let height_ft = height_part.parse::<u32>().ok().map(|h| h * 100);
            report.clouds.push(CloudLayer {
                coverage: coverage.to_string(),
                height_ft,
                cloud_type,
            });
            continue;
        }

        if token == "CLR" || token == "SKC" || token == "NCD" || token == "NSC" {
            report.clouds.push(CloudLayer {
                coverage: token.to_string(),
                height_ft: None,
                cloud_type: String::new(),
            });
            continue;
        }

        // Temperature/dew point: TT/TD (e.g., "24/16", "M01/M05")
        if token.contains('/') && !token.ends_with("SM") && !token.ends_with("KT") {
            let parts: Vec<&str> = token.split('/').collect();
            if parts.len() == 2 {
                if let Some(t) = parse_metar_temp(parts[0]) {
                    report.temperature = Some(t);
                }
                if let Some(d) = parse_metar_temp(parts[1]) {
                    report.dew_point = Some(d);
                }
            }
            continue;
        }

        // Altimeter: Annnn
        if token.starts_with('A') && token.len() == 5 {
            if let Ok(alt) = token.strip_prefix('A').unwrap_or("").parse::<f64>() {
                report.altimeter = Some(alt / 100.0);
            }
            continue;
        }

        // Weather phenomena (simplified)
        let weather_codes = [
            "RA", "SN", "DZ", "FG", "BR", "HZ", "TS", "SH", "GR", "GS", "FZ", "PL", "SG", "IC",
            "PE", "UP", "FU", "VA", "DU", "SA", "SS", "DS", "SQ", "FC", "PO",
        ];
        for code in &weather_codes {
            if token.contains(code) {
                report.weather.push(token.to_string());
                break;
            }
        }
    }

    Some(report)
}

/// Parse wind from a METAR token.
fn parse_metar_wind(token: &str, report: &mut MetarReport) {
    let wind_str = token.strip_suffix("KT").unwrap_or(token);

    if wind_str.contains('G') {
        // Gusts present
        let parts: Vec<&str> = wind_str.split('G').collect();
        if parts.len() == 2 && parts[0].len() >= 5 {
            let dir_str = &parts[0][..3];
            let spd_str = &parts[0][3..];
            if dir_str == "VRB" {
                report.wind_direction = None;
            } else {
                report.wind_direction = dir_str.parse().ok();
            }
            report.wind_speed = spd_str.parse().unwrap_or(0);
            report.wind_gust = parts[1].parse().ok();
        }
    } else if wind_str.len() >= 5 {
        let dir_str = &wind_str[..3];
        let spd_str = &wind_str[3..];
        if dir_str == "VRB" {
            report.wind_direction = None;
        } else {
            report.wind_direction = dir_str.parse().ok();
        }
        report.wind_speed = spd_str.parse().unwrap_or(0);
    }
}

/// Parse a METAR temperature string (e.g., "24", "M01" for -1).
fn parse_metar_temp(s: &str) -> Option<f64> {
    if let Some(rest) = s.strip_prefix('M') {
        rest.parse::<f64>().ok().map(|v| -v)
    } else {
        s.parse::<f64>().ok()
    }
}

// ════════════════════════════════════════════════════════════════════════════
// Sounding (TEMP) Data
// ════════════════════════════════════════════════════════════════════════════

/// A single level in a radiosonde sounding.
#[derive(Debug, Clone, Copy)]
pub struct SoundingLevel {
    /// Pressure (hPa / mb).
    pub pressure: f64,
    /// Height (m above sea level).
    pub height: f64,
    /// Temperature (Celsius).
    pub temperature: f64,
    /// Dew point (Celsius).
    pub dew_point: f64,
    /// Wind direction (degrees).
    pub wind_direction: f64,
    /// Wind speed (knots).
    pub wind_speed: f64,
}

/// A complete radiosonde sounding profile.
#[derive(Debug, Clone)]
pub struct SoundingProfile {
    /// Station identifier.
    pub station: String,
    /// Observation time (YYYYMMDDHH).
    pub observation_time: String,
    /// Station latitude (degrees).
    pub latitude: f64,
    /// Station longitude (degrees).
    pub longitude: f64,
    /// Station elevation (m).
    pub elevation: f64,
    /// Sounding levels (sorted by decreasing pressure / increasing height).
    pub levels: Vec<SoundingLevel>,
}

impl SoundingProfile {
    /// Create a new sounding profile.
    pub fn new(station: &str) -> Self {
        Self {
            station: station.to_string(),
            observation_time: String::new(),
            latitude: 0.0,
            longitude: 0.0,
            elevation: 0.0,
            levels: Vec::new(),
        }
    }

    /// Add a sounding level.
    pub fn add_level(&mut self, level: SoundingLevel) {
        self.levels.push(level);
    }

    /// Sort levels by pressure (descending = surface first).
    pub fn sort_by_pressure(&mut self) {
        self.levels.sort_by(|a, b| {
            b.pressure
                .partial_cmp(&a.pressure)
                .unwrap_or(std::cmp::Ordering::Equal)
        });
    }

    /// Sort levels by height (ascending = surface first).
    pub fn sort_by_height(&mut self) {
        self.levels.sort_by(|a, b| {
            a.height
                .partial_cmp(&b.height)
                .unwrap_or(std::cmp::Ordering::Equal)
        });
    }

    /// Interpolate temperature at a given pressure level using linear interpolation.
    pub fn interpolate_temperature(&self, pressure: f64) -> Option<f64> {
        if self.levels.len() < 2 {
            return None;
        }
        interpolate_sounding_value(&self.levels, pressure, |l| l.temperature)
    }

    /// Interpolate wind speed at a given pressure level.
    pub fn interpolate_wind_speed(&self, pressure: f64) -> Option<f64> {
        if self.levels.len() < 2 {
            return None;
        }
        interpolate_sounding_value(&self.levels, pressure, |l| l.wind_speed)
    }

    /// Interpolate height at a given pressure level.
    pub fn interpolate_height(&self, pressure: f64) -> Option<f64> {
        if self.levels.len() < 2 {
            return None;
        }
        interpolate_sounding_value(&self.levels, pressure, |l| l.height)
    }

    /// Compute the Lifting Condensation Level (LCL) pressure.
    ///
    /// Uses the surface level temperature and dew point.
    pub fn compute_lcl_pressure(&self) -> Option<f64> {
        if self.levels.is_empty() {
            return None;
        }
        let surface = &self.levels[0];
        let t = surface.temperature;
        let td = surface.dew_point;
        let p = surface.pressure;

        // Bolton (1980) approximation for LCL temperature
        let t_lcl = 1.0 / (1.0 / (td + 273.15 - 56.0) + (t - td).ln() / 800.0) + 56.0;
        // LCL pressure using Poisson's equation
        let p_lcl = p * ((t_lcl) / (t + 273.15)).powf(3.5);
        Some(p_lcl)
    }

    /// Compute the total precipitable water (mm).
    ///
    /// Integrates the mixing ratio over all levels.
    pub fn precipitable_water(&self) -> f64 {
        if self.levels.len() < 2 {
            return 0.0;
        }
        let g = 9.80665;
        let mut pw = 0.0;
        for i in 0..self.levels.len() - 1 {
            let dp = (self.levels[i].pressure - self.levels[i + 1].pressure).abs() * 100.0; // hPa -> Pa
            let w1 = mixing_ratio(self.levels[i].dew_point, self.levels[i].pressure);
            let w2 = mixing_ratio(self.levels[i + 1].dew_point, self.levels[i + 1].pressure);
            let w_avg = (w1 + w2) / 2.0;
            pw += w_avg * dp / g;
        }
        // Convert kg/m^2 to mm (1 kg/m^2 = 1 mm)
        pw
    }

    /// Serialize the sounding to a simple text format.
    pub fn to_text(&self) -> String {
        let mut out = String::new();
        out.push_str(&format!("# SOUNDING {}\n", self.station));
        out.push_str(&format!("# TIME {}\n", self.observation_time));
        out.push_str(&format!(
            "# LAT {:.4} LON {:.4} ELEV {:.1}\n",
            self.latitude, self.longitude, self.elevation
        ));
        out.push_str("# PRES(hPa) HGT(m) TEMP(C) DWPT(C) WDIR(deg) WSPD(kt)\n");
        for level in &self.levels {
            out.push_str(&format!(
                "{:.1} {:.1} {:.1} {:.1} {:.0} {:.0}\n",
                level.pressure,
                level.height,
                level.temperature,
                level.dew_point,
                level.wind_direction,
                level.wind_speed
            ));
        }
        out
    }

    /// Parse a sounding from text format.
    pub fn from_text(text: &str) -> Option<Self> {
        let mut profile = SoundingProfile::new("");
        for line in text.lines() {
            let line = line.trim();
            if let Some(rest) = line.strip_prefix("# SOUNDING ") {
                profile.station = rest.trim().to_string();
            } else if let Some(rest) = line.strip_prefix("# TIME ") {
                profile.observation_time = rest.trim().to_string();
            } else if line.starts_with("# LAT ") {
                let parts: Vec<&str> = line.split_whitespace().collect();
                if parts.len() >= 6 {
                    profile.latitude = parts[2].parse().unwrap_or(0.0);
                    profile.longitude = parts[4].parse().unwrap_or(0.0);
                    profile.elevation = parts[6].parse().unwrap_or(0.0);
                }
            } else if line.starts_with('#') || line.is_empty() {
                continue;
            } else {
                let vals: Vec<f64> = line
                    .split_whitespace()
                    .filter_map(|s| s.parse().ok())
                    .collect();
                if vals.len() >= 6 {
                    profile.add_level(SoundingLevel {
                        pressure: vals[0],
                        height: vals[1],
                        temperature: vals[2],
                        dew_point: vals[3],
                        wind_direction: vals[4],
                        wind_speed: vals[5],
                    });
                }
            }
        }
        if profile.station.is_empty() {
            return None;
        }
        Some(profile)
    }
}

/// Interpolate a sounding value at a given pressure.
fn interpolate_sounding_value(
    levels: &[SoundingLevel],
    pressure: f64,
    extract: fn(&SoundingLevel) -> f64,
) -> Option<f64> {
    // Find bracketing levels (levels sorted by decreasing pressure)
    for i in 0..levels.len() - 1 {
        let p1 = levels[i].pressure;
        let p2 = levels[i + 1].pressure;
        if (pressure <= p1 && pressure >= p2) || (pressure >= p1 && pressure <= p2) {
            let log_p = pressure.ln();
            let log_p1 = p1.ln();
            let log_p2 = p2.ln();
            let denom = log_p2 - log_p1;
            if denom.abs() < 1e-15 {
                return Some(extract(&levels[i]));
            }
            let t = (log_p - log_p1) / denom;
            let v1 = extract(&levels[i]);
            let v2 = extract(&levels[i + 1]);
            return Some(v1 + t * (v2 - v1));
        }
    }
    None
}

// ════════════════════════════════════════════════════════════════════════════
// Wind Profile Models
// ════════════════════════════════════════════════════════════════════════════

/// Compute wind speed at height using the power law profile.
///
/// V(z) = V_ref * (z / z_ref)^alpha
///
/// `v_ref` - Reference wind speed (m/s) at `z_ref`.
/// `z_ref` - Reference height (m).
/// `z` - Target height (m).
/// `alpha` - Power law exponent (typically 0.14 for open terrain, 0.25 for urban).
pub fn wind_power_law(v_ref: f64, z_ref: f64, z: f64, alpha: f64) -> f64 {
    if z_ref < 1e-10 || z < 0.0 {
        return 0.0;
    }
    v_ref * (z / z_ref).powf(alpha)
}

/// Compute wind speed at height using the logarithmic law profile.
///
/// V(z) = (u* / k) * ln(z / z0)
///
/// `u_star` - Friction velocity (m/s).
/// `z` - Target height (m).
/// `z0` - Roughness length (m).
/// `von_karman` - Von Karman constant (typically 0.41).
pub fn wind_log_law(u_star: f64, z: f64, z0: f64, von_karman: f64) -> f64 {
    if z <= z0 || z0 < 1e-15 || von_karman.abs() < 1e-15 {
        return 0.0;
    }
    (u_star / von_karman) * (z / z0).ln()
}

/// Compute the friction velocity from a reference wind speed.
///
/// u* = V_ref * k / ln(z_ref / z0)
pub fn friction_velocity(v_ref: f64, z_ref: f64, z0: f64, von_karman: f64) -> f64 {
    if z_ref <= z0 || z0 < 1e-15 {
        return 0.0;
    }
    v_ref * von_karman / (z_ref / z0).ln()
}

/// Compute the Ekman spiral wind direction change with height.
///
/// The wind direction rotates clockwise (Northern Hemisphere) with height
/// due to the decrease in friction.
///
/// `z` - Height (m).
/// `boundary_layer_height` - Atmospheric boundary layer height (m).
/// `surface_direction` - Surface wind direction (degrees).
/// `geostrophic_direction` - Geostrophic wind direction (degrees).
pub fn ekman_wind_direction(
    z: f64,
    boundary_layer_height: f64,
    surface_direction: f64,
    geostrophic_direction: f64,
) -> f64 {
    if boundary_layer_height < 1e-10 {
        return geostrophic_direction;
    }
    let ratio = (z / boundary_layer_height).min(1.0);
    let delta = geostrophic_direction - surface_direction;
    surface_direction + delta * ratio
}

// ════════════════════════════════════════════════════════════════════════════
// Temperature Profile Models
// ════════════════════════════════════════════════════════════════════════════

/// Compute temperature at altitude using a constant lapse rate.
///
/// T(h) = T_surface - lapse_rate * h
///
/// `t_surface` - Surface temperature (K or C).
/// `lapse_rate` - Temperature lapse rate (K/m or C/m, positive means decreasing with height).
/// `height` - Height above surface (m).
pub fn temperature_lapse(t_surface: f64, lapse_rate: f64, height: f64) -> f64 {
    t_surface - lapse_rate * height
}

/// Compute the potential temperature.
///
/// theta = T * (P0 / P)^(R/cp)
///
/// `temperature` - Temperature (K).
/// `pressure` - Pressure (hPa).
/// `p0` - Reference pressure (typically 1000 hPa).
pub fn potential_temperature(temperature: f64, pressure: f64, p0: f64) -> f64 {
    if pressure < 1e-10 {
        return temperature;
    }
    temperature * (p0 / pressure).powf(0.286)
}

/// Compute the virtual temperature.
///
/// Tv = T * (1 + 0.61 * w)
///
/// `temperature` - Temperature (K).
/// `mixing_ratio` - Water vapor mixing ratio (kg/kg).
pub fn virtual_temperature(temperature: f64, mixing_ratio_val: f64) -> f64 {
    temperature * (1.0 + 0.61 * mixing_ratio_val)
}

/// Compute the equivalent potential temperature (Bolton 1980).
///
/// `temperature` - Temperature (K).
/// `dew_point` - Dew point temperature (K).
/// `pressure` - Pressure (hPa).
pub fn equivalent_potential_temperature(temperature: f64, dew_point: f64, pressure: f64) -> f64 {
    let w = mixing_ratio(dew_point - 273.15, pressure); // Convert K to C for mixing_ratio
    let theta = potential_temperature(temperature, pressure, 1000.0);
    let t_lcl = 1.0 / (1.0 / (dew_point - 56.0) + (temperature / dew_point).ln() / 800.0) + 56.0;
    theta
        * (3.376 / t_lcl - 0.00254)
            .exp()
            .powf(w * 1000.0 * (1.0 + 0.81e-3 * w * 1000.0))
}

// ════════════════════════════════════════════════════════════════════════════
// Pressure-Altitude Computations
// ════════════════════════════════════════════════════════════════════════════

/// Compute pressure altitude from actual pressure and standard sea-level pressure.
///
/// Uses the hypsometric equation simplified for standard atmosphere.
///
/// `pressure` - Actual pressure (hPa).
/// `sea_level_pressure` - Sea-level pressure (hPa, standard = 1013.25).
///
/// Returns altitude in meters.
pub fn pressure_altitude(pressure: f64, sea_level_pressure: f64) -> f64 {
    if pressure < 1e-10 || sea_level_pressure < 1e-10 {
        return 0.0;
    }
    // International barometric formula for the troposphere
    let t0 = 288.15; // ISA sea level temperature (K)
    let lapse = 0.0065; // Standard lapse rate (K/m)
    let g = 9.80665;
    let r = 287.053;

    (t0 / lapse) * (1.0 - (pressure / sea_level_pressure).powf(r * lapse / g))
}

/// Compute density altitude.
///
/// `pressure_alt` - Pressure altitude (m).
/// `oat` - Outside air temperature (C).
///
/// Returns density altitude in meters.
pub fn density_altitude(pressure_alt: f64, oat: f64) -> f64 {
    let isa_temp = 15.0 - 0.0065 * pressure_alt; // ISA temp at this altitude
    let temp_dev = oat - isa_temp;
    pressure_alt + 120.0 * temp_dev
}

/// Convert between pressure and flight level.
///
/// `pressure` - Pressure in hPa.
///
/// Returns flight level (in hundreds of feet).
pub fn pressure_to_flight_level(pressure: f64) -> f64 {
    let alt_ft = pressure_altitude(pressure, 1013.25) * 3.28084; // meters to feet
    (alt_ft / 100.0).round()
}

/// Convert flight level to pressure (hPa).
///
/// `flight_level` - Flight level (hundreds of feet).
pub fn flight_level_to_pressure(flight_level: f64) -> f64 {
    let alt_m = flight_level * 100.0 / 3.28084; // feet to meters
    let t0 = 288.15;
    let lapse = 0.0065;
    let g = 9.80665;
    let r = 287.053;
    1013.25 * (1.0 - lapse * alt_m / t0).powf(g / (r * lapse))
}

// ════════════════════════════════════════════════════════════════════════════
// ISA Atmosphere Model
// ════════════════════════════════════════════════════════════════════════════

/// International Standard Atmosphere (ISA) properties at a given altitude.
#[derive(Debug, Clone, Copy)]
pub struct IsaProperties {
    /// Temperature (K).
    pub temperature: f64,
    /// Pressure (Pa).
    pub pressure: f64,
    /// Density (kg/m^3).
    pub density: f64,
    /// Speed of sound (m/s).
    pub speed_of_sound: f64,
    /// Dynamic viscosity (Pa*s).
    pub dynamic_viscosity: f64,
}

/// Compute ISA properties at a given geometric altitude.
///
/// Valid for altitudes from -2000m to 86000m. The atmosphere is divided
/// into layers with constant lapse rates.
///
/// `altitude` - Geometric altitude (m above sea level).
pub fn isa_properties(altitude: f64) -> IsaProperties {
    // ISA layer definitions: (base_altitude, base_temperature, lapse_rate)
    // Layer 0: Troposphere (0 to 11000m)
    // Layer 1: Tropopause (11000 to 20000m)
    // Layer 2: Stratosphere 1 (20000 to 32000m)
    // Layer 3: Stratosphere 2 (32000 to 47000m)
    // Layer 4: Stratopause (47000 to 51000m)
    // Layer 5: Mesosphere 1 (51000 to 71000m)
    // Layer 6: Mesosphere 2 (71000 to 86000m)

    let layers: [(f64, f64, f64, f64); 7] = [
        (0.0, 288.15, -0.0065, 101325.0),
        (11000.0, 216.65, 0.0, 22632.1),
        (20000.0, 216.65, 0.001, 5474.89),
        (32000.0, 228.65, 0.0028, 868.019),
        (47000.0, 270.65, 0.0, 110.906),
        (51000.0, 270.65, -0.0028, 66.9389),
        (71000.0, 214.65, -0.002, 3.95642),
    ];

    let r = 287.053; // Gas constant for dry air (J/(kg*K))
    let g = 9.80665;
    let gamma = 1.4; // Ratio of specific heats

    let alt = altitude.clamp(-2000.0, 86000.0);

    // Find the applicable layer
    let layer_idx = layers
        .iter()
        .enumerate()
        .rev()
        .find(|(_, layer)| alt >= layer.0)
        .map(|(i, _)| i)
        .unwrap_or(0);

    let (h_base, t_base, lapse, p_base) = layers[layer_idx];
    let dh = alt - h_base;

    let temperature;
    let pressure;

    if lapse.abs() < 1e-10 {
        // Isothermal layer
        temperature = t_base;
        pressure = p_base * E.powf(-g * dh / (r * t_base));
    } else {
        // Gradient layer
        temperature = t_base + lapse * dh;
        pressure = p_base * (temperature / t_base).powf(-g / (lapse * r));
    }

    let density = pressure / (r * temperature);
    let speed_of_sound = (gamma * r * temperature).sqrt();

    // Sutherland's law for dynamic viscosity
    let mu_ref = 1.716e-5;
    let t_ref = 273.15;
    let s = 110.4;
    let dynamic_viscosity =
        mu_ref * (temperature / t_ref).powf(1.5) * (t_ref + s) / (temperature + s);

    IsaProperties {
        temperature,
        pressure,
        density,
        speed_of_sound,
        dynamic_viscosity,
    }
}

/// Compute ISA temperature at altitude (simplified, troposphere only).
pub fn isa_temperature(altitude: f64) -> f64 {
    let alt = altitude.clamp(-2000.0, 86000.0);
    if alt <= 11000.0 {
        288.15 - 0.0065 * alt
    } else if alt <= 20000.0 {
        216.65
    } else if alt <= 32000.0 {
        216.65 + 0.001 * (alt - 20000.0)
    } else {
        228.65 + 0.0028 * (alt - 32000.0)
    }
}

/// Compute ISA pressure at altitude (simplified, troposphere only).
pub fn isa_pressure(altitude: f64) -> f64 {
    isa_properties(altitude).pressure
}

/// Compute ISA density at altitude.
pub fn isa_density(altitude: f64) -> f64 {
    isa_properties(altitude).density
}

// ════════════════════════════════════════════════════════════════════════════
// Humidity / Dew Point Conversions
// ════════════════════════════════════════════════════════════════════════════

/// Compute saturation vapor pressure using the Magnus formula.
///
/// es(T) = 6.112 * exp(17.67 * T / (T + 243.5))
///
/// `temperature` - Temperature in Celsius.
///
/// Returns saturation vapor pressure in hPa.
pub fn saturation_vapor_pressure(temperature: f64) -> f64 {
    6.112 * (17.67 * temperature / (temperature + 243.5)).exp()
}

/// Compute actual vapor pressure from dew point.
///
/// `dew_point` - Dew point temperature in Celsius.
pub fn vapor_pressure_from_dewpoint(dew_point: f64) -> f64 {
    saturation_vapor_pressure(dew_point)
}

/// Compute relative humidity from temperature and dew point.
///
/// RH = (es(Td) / es(T)) * 100
///
/// Returns relative humidity in percent (0-100).
pub fn relative_humidity(temperature: f64, dew_point: f64) -> f64 {
    let es_td = saturation_vapor_pressure(dew_point);
    let es_t = saturation_vapor_pressure(temperature);
    if es_t < 1e-10 {
        return 0.0;
    }
    (es_td / es_t * 100.0).clamp(0.0, 100.0)
}

/// Compute dew point from temperature and relative humidity.
///
/// Uses the inverse of the Magnus formula.
///
/// `temperature` - Temperature in Celsius.
/// `rh` - Relative humidity in percent (0-100).
///
/// Returns dew point in Celsius.
pub fn dew_point_from_rh(temperature: f64, rh: f64) -> f64 {
    let rh_frac = (rh / 100.0).clamp(0.001, 1.0);
    let a = 17.67;
    let b = 243.5;
    let gamma = a * temperature / (b + temperature) + rh_frac.ln();
    b * gamma / (a - gamma)
}

/// Compute mixing ratio from dew point and pressure.
///
/// w = 0.622 * e / (P - e)
///
/// `dew_point` - Dew point in Celsius.
/// `pressure` - Atmospheric pressure in hPa.
///
/// Returns mixing ratio in kg/kg.
pub fn mixing_ratio(dew_point: f64, pressure: f64) -> f64 {
    let e = saturation_vapor_pressure(dew_point);
    let denom = pressure - e;
    if denom < 1e-10 {
        return 0.0;
    }
    0.622 * e / denom
}

/// Compute specific humidity from mixing ratio.
///
/// q = w / (1 + w)
pub fn specific_humidity(mixing_ratio_val: f64) -> f64 {
    mixing_ratio_val / (1.0 + mixing_ratio_val)
}

/// Compute wet-bulb temperature (approximate, Stull 2011).
///
/// `temperature` - Temperature in Celsius.
/// `rh` - Relative humidity in percent.
///
/// Returns wet-bulb temperature in Celsius.
pub fn wet_bulb_temperature(temperature: f64, rh: f64) -> f64 {
    let t = temperature;
    let rh_val = rh.clamp(0.0, 100.0);
    t * (0.151977 * (rh_val + 8.313659_f64).sqrt()).atan() + (t + rh_val).atan()
        - (rh_val - 1.676331).atan()
        + 0.00391838 * rh_val.powf(1.5) * (0.023101 * rh_val).atan()
        - 4.686035
}

/// Compute heat index (apparent temperature, Rothfusz 1990).
///
/// `temperature` - Temperature in Fahrenheit.
/// `rh` - Relative humidity in percent.
///
/// Returns heat index in Fahrenheit.
pub fn heat_index_fahrenheit(temperature: f64, rh: f64) -> f64 {
    let t = temperature;
    let r = rh;

    if t < 80.0 {
        return t;
    }

    -42.379 + 2.04901523 * t + 10.14333127 * r
        - 0.22475541 * t * r
        - 6.83783e-3 * t * t
        - 5.481717e-2 * r * r
        + 1.22874e-3 * t * t * r
        + 8.5282e-4 * t * r * r
        - 1.99e-6 * t * t * r * r
}

/// Convert Celsius to Fahrenheit.
pub fn celsius_to_fahrenheit(c: f64) -> f64 {
    c * 9.0 / 5.0 + 32.0
}

/// Convert Fahrenheit to Celsius.
pub fn fahrenheit_to_celsius(f: f64) -> f64 {
    (f - 32.0) * 5.0 / 9.0
}

/// Convert Celsius to Kelvin.
pub fn celsius_to_kelvin(c: f64) -> f64 {
    c + 273.15
}

/// Convert Kelvin to Celsius.
pub fn kelvin_to_celsius(k: f64) -> f64 {
    k - 273.15
}

/// Compute wind chill (for temperatures below 10C and wind > 4.8 km/h).
///
/// Uses the North American wind chill index.
///
/// `temperature` - Temperature in Celsius.
/// `wind_speed_kmh` - Wind speed in km/h.
pub fn wind_chill(temperature: f64, wind_speed_kmh: f64) -> f64 {
    if temperature > 10.0 || wind_speed_kmh < 4.8 {
        return temperature;
    }
    let v = wind_speed_kmh;
    13.12 + 0.6215 * temperature - 11.37 * v.powf(0.16) + 0.3965 * temperature * v.powf(0.16)
}

/// Convert wind speed from knots to m/s.
pub fn knots_to_ms(knots: f64) -> f64 {
    knots * 0.514444
}

/// Convert wind speed from m/s to knots.
pub fn ms_to_knots(ms: f64) -> f64 {
    ms / 0.514444
}

/// Convert wind speed from km/h to m/s.
pub fn kmh_to_ms(kmh: f64) -> f64 {
    kmh / 3.6
}

/// Convert pressure from hPa to inHg.
pub fn hpa_to_inhg(hpa: f64) -> f64 {
    hpa * 0.02953
}

/// Convert pressure from inHg to hPa.
pub fn inhg_to_hpa(inhg: f64) -> f64 {
    inhg / 0.02953
}

// ════════════════════════════════════════════════════════════════════════════
// Stability Indices
// ════════════════════════════════════════════════════════════════════════════

/// Compute the Showalter Index.
///
/// SI = T_500 - T_parcel_500
///
/// Lifts a parcel from 850 hPa to 500 hPa.
///
/// `t850` - Temperature at 850 hPa (C).
/// `td850` - Dew point at 850 hPa (C).
/// `t500` - Temperature at 500 hPa (C).
///
/// Returns the Showalter Index (negative = unstable).
pub fn showalter_index(t850: f64, td850: f64, t500: f64) -> f64 {
    // Lift parcel dry-adiabatically from 850 to LCL, then moist-adiabatically to 500
    let _t_diff = t850 - td850;

    // Approximate LCL temperature
    let t_lcl = td850;
    // Approximate LCL pressure
    let p_lcl = 850.0 * ((t_lcl + 273.15) / (t850 + 273.15)).powf(3.5);

    // Dry adiabatic lapse from 850 to LCL
    let dry_lapse = 9.8 / 1000.0; // K/m
    let _h_lcl = (t850 - t_lcl) / dry_lapse;

    // Moist adiabatic lapse from LCL to 500 hPa (approximate 6 K/km)
    let moist_lapse = 6.0 / 1000.0;
    let _h_lcl_m = pressure_altitude(p_lcl, 1013.25);
    let _h_500 = pressure_altitude(500.0, 1013.25);

    // Simplified: lift moist adiabatically
    let t_parcel_500 = t_lcl
        - moist_lapse
            * (pressure_altitude(500.0, 1013.25) - pressure_altitude(p_lcl.max(500.0), 1013.25));

    t500 - t_parcel_500
}

/// Compute the K-Index (thunderstorm potential).
///
/// K = (T850 - T500) + Td850 - (T700 - Td700)
///
/// Higher K = more thunderstorm potential.
pub fn k_index(t850: f64, t500: f64, td850: f64, t700: f64, td700: f64) -> f64 {
    (t850 - t500) + td850 - (t700 - td700)
}

/// Compute the Total Totals index.
///
/// TT = (T850 - T500) + (Td850 - T500) = T850 + Td850 - 2*T500
pub fn total_totals(t850: f64, td850: f64, t500: f64) -> f64 {
    t850 + td850 - 2.0 * t500
}

/// Compute the SWEAT (Severe Weather Threat) index (simplified).
///
/// SWEAT = 12*Td850 + 20*(TT-49) + 2*f850 + f500 + 125*(S + 0.2)
///
/// (Simplified version - full version includes wind shear terms.)
pub fn sweat_index(td850: f64, total_totals_val: f64, wind850: f64, wind500: f64) -> f64 {
    let term1 = 12.0 * td850.max(0.0);
    let term2 = 20.0 * (total_totals_val - 49.0).max(0.0);
    let term3 = 2.0 * wind850;
    let term4 = wind500;
    // Simplified: omit wind shear direction term
    term1 + term2 + term3 + term4
}

// ════════════════════════════════════════════════════════════════════════════
// Beaufort Scale
// ════════════════════════════════════════════════════════════════════════════

/// Convert wind speed (m/s) to Beaufort scale number (0-12).
pub fn beaufort_number(wind_speed_ms: f64) -> u8 {
    let v = wind_speed_ms;
    if v < 0.3 {
        0
    } else if v < 1.6 {
        1
    } else if v < 3.4 {
        2
    } else if v < 5.5 {
        3
    } else if v < 8.0 {
        4
    } else if v < 10.8 {
        5
    } else if v < 13.9 {
        6
    } else if v < 17.2 {
        7
    } else if v < 20.8 {
        8
    } else if v < 24.5 {
        9
    } else if v < 28.5 {
        10
    } else if v < 32.7 {
        11
    } else {
        12
    }
}

/// Get the Beaufort scale description.
pub fn beaufort_description(number: u8) -> &'static str {
    match number {
        0 => "Calm",
        1 => "Light air",
        2 => "Light breeze",
        3 => "Gentle breeze",
        4 => "Moderate breeze",
        5 => "Fresh breeze",
        6 => "Strong breeze",
        7 => "High wind",
        8 => "Gale",
        9 => "Strong gale",
        10 => "Storm",
        11 => "Violent storm",
        12 => "Hurricane force",
        _ => "Unknown",
    }
}

// ════════════════════════════════════════════════════════════════════════════
// Cloud Base Estimation
// ════════════════════════════════════════════════════════════════════════════

/// Estimate cloud base height from surface temperature and dew point.
///
/// Uses the simple formula: cloud_base_ft = (T - Td) / 2.5 * 1000
///
/// `temperature` - Surface temperature (C).
/// `dew_point` - Surface dew point (C).
///
/// Returns estimated cloud base in feet AGL.
pub fn estimated_cloud_base_ft(temperature: f64, dew_point: f64) -> f64 {
    let spread = (temperature - dew_point).max(0.0);
    spread / 2.5 * 1000.0
}

/// Estimate cloud base height in meters AGL.
pub fn estimated_cloud_base_m(temperature: f64, dew_point: f64) -> f64 {
    estimated_cloud_base_ft(temperature, dew_point) * 0.3048
}

// ════════════════════════════════════════════════════════════════════════════
// Tests
// ════════════════════════════════════════════════════════════════════════════

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_grib2_indicator_valid() {
        let mut data = vec![0u8; 37];
        data[0] = b'G';
        data[1] = b'R';
        data[2] = b'I';
        data[3] = b'B';
        data[6] = 0; // discipline
        data[7] = 2; // edition
        // Total length = 37
        data[15] = 37;

        let result = parse_grib2_indicator(&data);
        assert!(result.is_some());
        let (disc, ed, len) = result.unwrap();
        assert_eq!(disc, 0);
        assert_eq!(ed, 2);
        assert_eq!(len, 37);
    }

    #[test]
    fn test_grib2_indicator_invalid_magic() {
        let data = vec![0u8; 16];
        assert!(parse_grib2_indicator(&data).is_none());
    }

    #[test]
    fn test_grib2_indicator_too_short() {
        let data = vec![b'G', b'R', b'I', b'B'];
        assert!(parse_grib2_indicator(&data).is_none());
    }

    #[test]
    fn test_grib2_end_marker() {
        let data = b"some_data7777";
        assert!(has_grib2_end_marker(data));
        let data2 = b"no_end_marker";
        assert!(!has_grib2_end_marker(data2));
    }

    #[test]
    fn test_grib2_discipline_name() {
        assert_eq!(grib2_discipline_name(0), "Meteorological");
        assert_eq!(grib2_discipline_name(10), "Oceanographic");
        assert_eq!(grib2_discipline_name(255), "Unknown");
    }

    #[test]
    fn test_metar_parse_basic() {
        let metar = "KJFK 301456Z 32015G25KT 10SM FEW050 SCT250 24/16 A3002";
        let report = parse_metar(metar).unwrap();
        assert_eq!(report.station, "KJFK");
        assert_eq!(report.day, 30);
        assert_eq!(report.hour, 14);
        assert_eq!(report.minute, 56);
        assert_eq!(report.wind_direction, Some(320));
        assert_eq!(report.wind_speed, 15);
        assert_eq!(report.wind_gust, Some(25));
        assert!((report.visibility_sm - 10.0).abs() < 0.1);
        assert!((report.temperature.unwrap() - 24.0).abs() < 0.1);
        assert!((report.dew_point.unwrap() - 16.0).abs() < 0.1);
        assert!((report.altimeter.unwrap() - 30.02).abs() < 0.01);
        assert_eq!(report.clouds.len(), 2);
    }

    #[test]
    fn test_metar_parse_negative_temp() {
        let metar = "CYUL 150000Z 36010KT 15SM OVC020 M05/M12 A2980";
        let report = parse_metar(metar).unwrap();
        assert!((report.temperature.unwrap() - (-5.0)).abs() < 0.1);
        assert!((report.dew_point.unwrap() - (-12.0)).abs() < 0.1);
    }

    #[test]
    fn test_metar_parse_auto() {
        let metar = "KORD 301556Z AUTO 28012KT 10SM CLR 22/10 A3010";
        let report = parse_metar(metar).unwrap();
        assert!(report.auto);
        assert_eq!(report.wind_direction, Some(280));
    }

    #[test]
    fn test_netcdf_variable() {
        let mut var = NetcdfAtmVariable::new(
            "temperature",
            "K",
            vec!["lat".to_string(), "lon".to_string()],
            vec![3, 4],
        );
        assert_eq!(var.total_size(), 12);

        var.set_value(0, 300.0);
        let v = var.get_value(0).unwrap();
        assert!((v - 300.0).abs() < 1e-6);

        let idx = var.flat_index(&[1, 2]).unwrap();
        assert_eq!(idx, 6); // 1*4 + 2
    }

    #[test]
    fn test_netcdf_dataset_roundtrip() {
        let mut ds = NetcdfAtmDataset::new();
        ds.add_attribute("source", "test");
        let mut var = NetcdfAtmVariable::new("temp", "K", vec!["level".to_string()], vec![3]);
        var.data = vec![300.0, 250.0, 200.0];
        ds.add_variable(var);

        let text = ds.to_text();
        let ds2 = NetcdfAtmDataset::from_text(&text).unwrap();
        assert_eq!(ds2.variables.len(), 1);
        assert_eq!(ds2.variables[0].name, "temp");
    }

    #[test]
    fn test_sounding_profile() {
        let mut profile = SoundingProfile::new("72451");
        profile.add_level(SoundingLevel {
            pressure: 1000.0,
            height: 100.0,
            temperature: 25.0,
            dew_point: 20.0,
            wind_direction: 180.0,
            wind_speed: 10.0,
        });
        profile.add_level(SoundingLevel {
            pressure: 850.0,
            height: 1500.0,
            temperature: 15.0,
            dew_point: 10.0,
            wind_direction: 200.0,
            wind_speed: 20.0,
        });
        profile.add_level(SoundingLevel {
            pressure: 500.0,
            height: 5500.0,
            temperature: -10.0,
            dew_point: -20.0,
            wind_direction: 270.0,
            wind_speed: 40.0,
        });

        let t = profile.interpolate_temperature(700.0);
        assert!(t.is_some());

        let pw = profile.precipitable_water();
        assert!(pw > 0.0);
    }

    #[test]
    fn test_sounding_text_roundtrip() {
        let mut profile = SoundingProfile::new("72451");
        profile.observation_time = "2026033012".to_string();
        profile.latitude = 40.78;
        profile.longitude = -73.97;
        profile.elevation = 10.0;
        profile.add_level(SoundingLevel {
            pressure: 1000.0,
            height: 10.0,
            temperature: 20.0,
            dew_point: 15.0,
            wind_direction: 180.0,
            wind_speed: 10.0,
        });
        let text = profile.to_text();
        let p2 = SoundingProfile::from_text(&text).unwrap();
        assert_eq!(p2.station, "72451");
        assert_eq!(p2.levels.len(), 1);
    }

    #[test]
    fn test_wind_power_law() {
        let v = wind_power_law(5.0, 10.0, 80.0, 0.14);
        // V = 5 * (80/10)^0.14 ≈ 5 * 1.337 ≈ 6.68
        assert!(v > 5.0 && v < 10.0, "v={v}");
    }

    #[test]
    fn test_wind_log_law() {
        let u_star = 0.5;
        let v = wind_log_law(u_star, 10.0, 0.03, 0.41);
        assert!(v > 0.0, "v={v}");
    }

    #[test]
    fn test_friction_velocity() {
        let u_star = friction_velocity(10.0, 10.0, 0.03, 0.41);
        assert!(u_star > 0.0);
        // Verify: wind_log_law with this u* at z=10 should give ~10
        let v = wind_log_law(u_star, 10.0, 0.03, 0.41);
        assert!((v - 10.0).abs() < 0.1, "v={v}");
    }

    #[test]
    fn test_temperature_lapse() {
        let t = temperature_lapse(288.15, 0.0065, 1000.0);
        assert!((t - 281.65).abs() < 1e-6);
    }

    #[test]
    fn test_potential_temperature() {
        let theta = potential_temperature(273.15, 1000.0, 1000.0);
        assert!((theta - 273.15).abs() < 0.01);
    }

    #[test]
    fn test_pressure_altitude_sea_level() {
        let alt = pressure_altitude(1013.25, 1013.25);
        assert!(alt.abs() < 1.0, "alt={alt}");
    }

    #[test]
    fn test_pressure_altitude_500hpa() {
        let alt = pressure_altitude(500.0, 1013.25);
        // ~5500m
        assert!(alt > 5000.0 && alt < 6000.0, "alt={alt}");
    }

    #[test]
    fn test_density_altitude() {
        let da = density_altitude(1000.0, 30.0);
        // ISA temp at 1000m = 15 - 6.5 = 8.5C, deviation = 30 - 8.5 = 21.5
        // DA = 1000 + 120*21.5 = 3580
        assert!((da - 3580.0).abs() < 1.0, "da={da}");
    }

    #[test]
    fn test_isa_sea_level() {
        let props = isa_properties(0.0);
        assert!((props.temperature - 288.15).abs() < 0.01);
        assert!((props.pressure - 101325.0).abs() < 1.0);
        assert!((props.density - 1.225).abs() < 0.01);
        assert!((props.speed_of_sound - 340.3).abs() < 1.0);
    }

    #[test]
    fn test_isa_tropopause() {
        let props = isa_properties(11000.0);
        assert!((props.temperature - 216.65).abs() < 0.1);
    }

    #[test]
    fn test_saturation_vapor_pressure() {
        let es = saturation_vapor_pressure(20.0);
        // Should be about 23.37 hPa
        assert!((es - 23.37).abs() < 0.5, "es={es}");
    }

    #[test]
    fn test_relative_humidity() {
        let rh = relative_humidity(20.0, 20.0);
        assert!((rh - 100.0).abs() < 0.1);
        let rh2 = relative_humidity(20.0, 10.0);
        assert!(rh2 < 100.0 && rh2 > 0.0);
    }

    #[test]
    fn test_dew_point_from_rh() {
        let td = dew_point_from_rh(20.0, 100.0);
        assert!((td - 20.0).abs() < 0.5, "td={td}");

        let td2 = dew_point_from_rh(20.0, 50.0);
        assert!(td2 < 20.0);
    }

    #[test]
    fn test_mixing_ratio() {
        let w = mixing_ratio(20.0, 1013.25);
        // Should be around 0.0147 kg/kg
        assert!(w > 0.01 && w < 0.02, "w={w}");
    }

    #[test]
    fn test_specific_humidity() {
        let q = specific_humidity(0.01);
        // q = 0.01 / 1.01 ≈ 0.0099
        assert!((q - 0.0099).abs() < 0.001, "q={q}");
    }

    #[test]
    fn test_celsius_fahrenheit_conversion() {
        assert!((celsius_to_fahrenheit(0.0) - 32.0).abs() < 1e-10);
        assert!((celsius_to_fahrenheit(100.0) - 212.0).abs() < 1e-10);
        assert!((fahrenheit_to_celsius(32.0) - 0.0).abs() < 1e-10);
    }

    #[test]
    fn test_celsius_kelvin_conversion() {
        assert!((celsius_to_kelvin(0.0) - 273.15).abs() < 1e-10);
        assert!((kelvin_to_celsius(273.15) - 0.0).abs() < 1e-10);
    }

    #[test]
    fn test_wind_chill() {
        let wc = wind_chill(-10.0, 30.0);
        assert!(wc < -10.0, "wc={wc}");
        // Above 10C, should return temperature
        let wc2 = wind_chill(15.0, 30.0);
        assert!((wc2 - 15.0).abs() < 1e-10);
    }

    #[test]
    fn test_knots_conversions() {
        let ms = knots_to_ms(1.0);
        assert!((ms - 0.514444).abs() < 1e-4);
        let kt = ms_to_knots(0.514444);
        assert!((kt - 1.0).abs() < 1e-4);
    }

    #[test]
    fn test_pressure_conversions() {
        let inhg = hpa_to_inhg(1013.25);
        assert!((inhg - 29.92).abs() < 0.1, "inhg={inhg}");
    }

    #[test]
    fn test_beaufort_scale() {
        assert_eq!(beaufort_number(0.0), 0);
        assert_eq!(beaufort_number(5.0), 3);
        assert_eq!(beaufort_number(35.0), 12);
        assert_eq!(beaufort_description(0), "Calm");
        assert_eq!(beaufort_description(12), "Hurricane force");
    }

    #[test]
    fn test_cloud_base_estimation() {
        let cb_ft = estimated_cloud_base_ft(25.0, 15.0);
        // spread = 10, base = 10/2.5*1000 = 4000 ft
        assert!((cb_ft - 4000.0).abs() < 1.0, "cb_ft={cb_ft}");
    }

    #[test]
    fn test_k_index() {
        let ki = k_index(20.0, -10.0, 15.0, 10.0, 5.0);
        // KI = (20-(-10)) + 15 - (10-5) = 30 + 15 - 5 = 40
        assert!((ki - 40.0).abs() < 1e-6, "ki={ki}");
    }

    #[test]
    fn test_total_totals() {
        let tt = total_totals(20.0, 15.0, -10.0);
        // TT = 20 + 15 - 2*(-10) = 55
        assert!((tt - 55.0).abs() < 1e-6);
    }

    #[test]
    fn test_flight_level() {
        let fl = pressure_to_flight_level(500.0);
        // 500 hPa ≈ FL180
        assert!(fl > 150.0 && fl < 200.0, "fl={fl}");
    }

    #[test]
    fn test_ekman_wind_direction() {
        let dir = ekman_wind_direction(500.0, 1000.0, 180.0, 270.0);
        // At z=500/1000=0.5, dir = 180 + 0.5*90 = 225
        assert!((dir - 225.0).abs() < 1e-6, "dir={dir}");
    }

    #[test]
    fn test_virtual_temperature() {
        let tv = virtual_temperature(300.0, 0.01);
        // Tv = 300 * (1 + 0.61*0.01) = 300 * 1.0061 = 301.83
        assert!((tv - 301.83).abs() < 0.01, "tv={tv}");
    }
}
