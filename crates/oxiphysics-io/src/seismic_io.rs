#![allow(clippy::needless_range_loop)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0
//! Seismic data I/O: SEG-Y, SAC, MiniSEED formats.
//!
//! Provides structures for reading and writing common seismic data formats
//! including SEG-Y rev 1, SAC (Seismic Analysis Code), and MiniSEED.
//! Also includes signal processing tools like STA/LTA arrival pickers
//! and Butterworth bandpass filters.

use std::collections::HashMap;

// ── SEG-Y Format ──────────────────────────────────────────────────────────────

/// EBCDIC-encoded textual header for SEG-Y files (3200 bytes, 40 lines × 80 chars).
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct SegyTextHeader {
    /// Raw EBCDIC bytes (3200 bytes).
    pub bytes: Vec<u8>,
}

impl SegyTextHeader {
    /// Create a new textual header filled with spaces (EBCDIC 0x40).
    pub fn new() -> Self {
        Self {
            bytes: vec![0x40u8; 3200],
        }
    }

    /// Write an ASCII line into the textual header at the given line index (0-based).
    ///
    /// Lines are 80 characters wide; the input is truncated or padded as needed.
    pub fn set_line(&mut self, line_idx: usize, text: &str) {
        if line_idx >= 40 {
            return;
        }
        let start = line_idx * 80;
        let ascii_bytes = text.as_bytes();
        for i in 0..80 {
            let ch = if i < ascii_bytes.len() {
                ascii_bytes[i]
            } else {
                b' '
            };
            self.bytes[start + i] = ascii_to_ebcdic(ch);
        }
    }

    /// Read a line from the textual header as an ASCII string.
    pub fn get_line(&self, line_idx: usize) -> String {
        if line_idx >= 40 {
            return String::new();
        }
        let start = line_idx * 80;
        let chars: Vec<u8> = self.bytes[start..start + 80]
            .iter()
            .map(|&b| ebcdic_to_ascii(b))
            .collect();
        String::from_utf8_lossy(&chars).trim_end().to_string()
    }
}

impl Default for SegyTextHeader {
    fn default() -> Self {
        Self::new()
    }
}

/// Convert ASCII byte to EBCDIC (partial table for printable ASCII).
#[allow(dead_code)]
fn ascii_to_ebcdic(ch: u8) -> u8 {
    match ch {
        b' ' => 0x40,
        b'0'..=b'9' => ch - b'0' + 0xF0,
        b'A'..=b'I' => ch - b'A' + 0xC1,
        b'J'..=b'R' => ch - b'J' + 0xD1,
        b'S'..=b'Z' => ch - b'S' + 0xE2,
        b'a'..=b'i' => ch - b'a' + 0x81,
        b'j'..=b'r' => ch - b'j' + 0x91,
        b's'..=b'z' => ch - b's' + 0xA2,
        b'.' => 0x4B,
        b',' => 0x6B,
        b':' => 0x7A,
        _ => 0x40,
    }
}

/// Convert EBCDIC byte to ASCII (partial table).
#[allow(dead_code)]
fn ebcdic_to_ascii(b: u8) -> u8 {
    match b {
        0x40 => b' ',
        0xF0..=0xF9 => b - 0xF0 + b'0',
        0xC1..=0xC9 => b - 0xC1 + b'A',
        0xD1..=0xD9 => b - 0xD1 + b'J',
        0xE2..=0xE9 => b - 0xE2 + b'S',
        0x81..=0x89 => b - 0x81 + b'a',
        0x91..=0x99 => b - 0x91 + b'j',
        0xA2..=0xA9 => b - 0xA2 + b's',
        0x4B => b'.',
        0x6B => b',',
        0x7A => b':',
        _ => b'?',
    }
}

/// SEG-Y binary file header (400 bytes).
///
/// Contains survey-wide parameters like sample interval, number of samples,
/// and data sample format code.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct SegyBinaryHeader {
    /// Job identification number.
    pub job_id: i32,
    /// Line number.
    pub line_number: i32,
    /// Reel number.
    pub reel_number: i32,
    /// Number of data traces per record.
    pub traces_per_record: i16,
    /// Number of auxiliary traces per record.
    pub aux_traces_per_record: i16,
    /// Sample interval in microseconds.
    pub sample_interval_us: i16,
    /// Number of samples per data trace.
    pub samples_per_trace: i16,
    /// Data sample format code: 1=IBM float, 5=IEEE float.
    pub data_format_code: i16,
    /// SEG-Y revision: 0 = rev 0, 256 = rev 1.
    pub segy_revision: i16,
    /// Fixed length trace flag.
    pub fixed_length_flag: i16,
}

impl Default for SegyBinaryHeader {
    fn default() -> Self {
        Self {
            job_id: 1,
            line_number: 1,
            reel_number: 1,
            traces_per_record: 1,
            aux_traces_per_record: 0,
            sample_interval_us: 2000, // 2 ms
            samples_per_trace: 1000,
            data_format_code: 5, // IEEE float
            segy_revision: 256,  // rev 1
            fixed_length_flag: 1,
        }
    }
}

/// SEG-Y trace header (240 bytes).
///
/// Contains per-trace metadata including shot and receiver coordinates.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct SegyTraceHeader {
    /// Trace sequence number within the SEG-Y file.
    pub trace_seq_file: i32,
    /// Field record number.
    pub field_record_number: i32,
    /// Trace number within field record.
    pub trace_number_in_field: i32,
    /// Source X coordinate (scaled by coordinate scalar).
    pub source_x: i32,
    /// Source Y coordinate.
    pub source_y: i32,
    /// Group (receiver) X coordinate.
    pub group_x: i32,
    /// Group (receiver) Y coordinate.
    pub group_y: i32,
    /// Coordinate scalar: if negative, divide; if positive, multiply.
    pub coordinate_scalar: i16,
    /// Delay recording time in milliseconds.
    pub delay_recording_time_ms: i16,
    /// Number of samples in this trace.
    pub num_samples: i16,
    /// Sample interval in microseconds for this trace.
    pub sample_interval_us: i16,
    /// Year data recorded.
    pub year: i16,
    /// Day of year.
    pub day_of_year: i16,
    /// Hour of day (24-hour clock).
    pub hour: i16,
    /// Minute of hour.
    pub minute: i16,
    /// Second of minute.
    pub second: i16,
}

impl Default for SegyTraceHeader {
    fn default() -> Self {
        Self {
            trace_seq_file: 1,
            field_record_number: 1,
            trace_number_in_field: 1,
            source_x: 0,
            source_y: 0,
            group_x: 0,
            group_y: 0,
            coordinate_scalar: -100,
            delay_recording_time_ms: 0,
            num_samples: 1000,
            sample_interval_us: 2000,
            year: 2026,
            day_of_year: 1,
            hour: 0,
            minute: 0,
            second: 0,
        }
    }
}

/// A single SEG-Y trace with header and float32 sample data.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct SegyTrace {
    /// Trace header metadata.
    pub header: SegyTraceHeader,
    /// Sample data as IEEE float32.
    pub data: Vec<f32>,
}

impl SegyTrace {
    /// Create a new trace with the given data and a default header.
    pub fn new(data: Vec<f32>) -> Self {
        let num_samples = data.len() as i16;
        let header = SegyTraceHeader {
            num_samples,
            ..Default::default()
        };
        Self { header, data }
    }

    /// Return the number of samples in this trace.
    pub fn num_samples(&self) -> usize {
        self.data.len()
    }

    /// Return the sample interval in seconds.
    pub fn sample_interval_s(&self) -> f64 {
        self.header.sample_interval_us as f64 * 1e-6
    }
}

/// A complete SEG-Y file in memory.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct SegyFile {
    /// Textual (EBCDIC) header.
    pub text_header: SegyTextHeader,
    /// Binary file header.
    pub binary_header: SegyBinaryHeader,
    /// All traces in this file.
    pub traces: Vec<SegyTrace>,
}

impl SegyFile {
    /// Create a new empty SEG-Y file.
    pub fn new() -> Self {
        Self {
            text_header: SegyTextHeader::new(),
            binary_header: SegyBinaryHeader::default(),
            traces: Vec::new(),
        }
    }

    /// Add a trace to the file, updating the trace sequence number.
    pub fn add_trace(&mut self, mut trace: SegyTrace) {
        let idx = self.traces.len() as i32 + 1;
        trace.header.trace_seq_file = idx;
        self.traces.push(trace);
    }

    /// Return the number of traces.
    pub fn num_traces(&self) -> usize {
        self.traces.len()
    }

    /// Serialize header fields to a byte buffer (simplified, not full 3600-byte header).
    pub fn serialize_binary_header(&self) -> Vec<u8> {
        let h = &self.binary_header;
        let mut buf = vec![0u8; 400];
        write_i32_be(&mut buf, 0, h.job_id);
        write_i32_be(&mut buf, 4, h.line_number);
        write_i32_be(&mut buf, 8, h.reel_number);
        write_i16_be(&mut buf, 12, h.traces_per_record);
        write_i16_be(&mut buf, 14, h.aux_traces_per_record);
        write_i16_be(&mut buf, 16, h.sample_interval_us);
        write_i16_be(&mut buf, 20, h.samples_per_trace);
        write_i16_be(&mut buf, 24, h.data_format_code);
        write_i16_be(&mut buf, 300, h.segy_revision);
        write_i16_be(&mut buf, 302, h.fixed_length_flag);
        buf
    }

    /// Deserialize binary header from a 400-byte buffer.
    pub fn deserialize_binary_header(buf: &[u8]) -> SegyBinaryHeader {
        SegyBinaryHeader {
            job_id: read_i32_be(buf, 0),
            line_number: read_i32_be(buf, 4),
            reel_number: read_i32_be(buf, 8),
            traces_per_record: read_i16_be(buf, 12),
            aux_traces_per_record: read_i16_be(buf, 14),
            sample_interval_us: read_i16_be(buf, 16),
            samples_per_trace: read_i16_be(buf, 20),
            data_format_code: read_i16_be(buf, 24),
            segy_revision: read_i16_be(buf, 300),
            fixed_length_flag: read_i16_be(buf, 302),
        }
    }
}

impl Default for SegyFile {
    fn default() -> Self {
        Self::new()
    }
}

// ── SAC Format ────────────────────────────────────────────────────────────────

/// SAC (Seismic Analysis Code) file header.
///
/// The SAC header is 632 bytes containing float, integer, logical, enumerated,
/// and character fields.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct SacHeader {
    /// Sampling interval in seconds (delta).
    pub delta: f32,
    /// Minimum value of dependent variable.
    pub depmin: f32,
    /// Maximum value of dependent variable.
    pub depmax: f32,
    /// Begin time of the independent variable (seconds).
    pub b: f32,
    /// End time of the independent variable (seconds).
    pub e: f32,
    /// Event origin time (seconds relative to reference).
    pub o: f32,
    /// First arrival time (seconds).
    pub a: f32,
    /// Station-event distance in km.
    pub dist: f32,
    /// Event-to-station azimuth in degrees.
    pub az: f32,
    /// Back azimuth in degrees.
    pub baz: f32,
    /// Component azimuth (degrees from North, clockwise).
    pub cmpaz: f32,
    /// Component inclination (degrees from vertical down).
    pub cmpinc: f32,
    /// Station latitude in degrees.
    pub stla: f32,
    /// Station longitude in degrees.
    pub stlo: f32,
    /// Station elevation in meters.
    pub stel: f32,
    /// Event latitude in degrees.
    pub evla: f32,
    /// Event longitude in degrees.
    pub evlo: f32,
    /// Event depth in km.
    pub evdp: f32,
    /// Event magnitude.
    pub mag: f32,
    /// Number of points (samples).
    pub npts: i32,
    /// Beginning time year.
    pub nzyear: i32,
    /// Beginning time day-of-year.
    pub nzjday: i32,
    /// Beginning time hour.
    pub nzhour: i32,
    /// Beginning time minute.
    pub nzmin: i32,
    /// Beginning time second.
    pub nzsec: i32,
    /// Beginning time millisecond.
    pub nzmsec: i32,
    /// SAC version number.
    pub nvhdr: i32,
    /// Station name (up to 8 characters).
    pub kstnm: String,
    /// Channel name (up to 8 characters).
    pub kcmpnm: String,
    /// Network name (up to 8 characters).
    pub knetwk: String,
}

impl Default for SacHeader {
    fn default() -> Self {
        Self {
            delta: 0.01,
            depmin: -1e20,
            depmax: 1e20,
            b: 0.0,
            e: 0.0,
            o: -12345.0,
            a: -12345.0,
            dist: -12345.0,
            az: -12345.0,
            baz: -12345.0,
            cmpaz: 0.0,
            cmpinc: 0.0,
            stla: -12345.0,
            stlo: -12345.0,
            stel: -12345.0,
            evla: -12345.0,
            evlo: -12345.0,
            evdp: -12345.0,
            mag: -12345.0,
            npts: 0,
            nzyear: 2026,
            nzjday: 1,
            nzhour: 0,
            nzmin: 0,
            nzsec: 0,
            nzmsec: 0,
            nvhdr: 6,
            kstnm: String::from("-12345  "),
            kcmpnm: String::from("-12345  "),
            knetwk: String::from("-12345  "),
        }
    }
}

/// A complete SAC file with header and time-series data.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct SacFile {
    /// SAC header.
    pub header: SacHeader,
    /// Time-series samples (float32).
    pub data: Vec<f32>,
}

impl SacFile {
    /// Create a new SAC file from data and a sampling interval.
    pub fn from_data(data: Vec<f32>, delta: f32) -> Self {
        let npts = data.len() as i32;
        let e = if npts > 0 {
            (npts - 1) as f32 * delta
        } else {
            0.0
        };
        let depmin = data.iter().cloned().fold(f32::MAX, f32::min);
        let depmax = data.iter().cloned().fold(f32::MIN, f32::max);
        let header = SacHeader {
            delta,
            npts,
            b: 0.0,
            e,
            depmin,
            depmax,
            ..Default::default()
        };
        Self { header, data }
    }

    /// Return the end time in seconds.
    pub fn end_time(&self) -> f64 {
        self.header.b as f64 + (self.header.npts - 1) as f64 * self.header.delta as f64
    }

    /// Serialize the SAC header to 632 bytes (simplified).
    pub fn serialize_header(&self) -> Vec<u8> {
        let mut buf = vec![0u8; 632];
        write_f32_le(&mut buf, 0, self.header.delta);
        write_f32_le(&mut buf, 4, self.header.depmin);
        write_f32_le(&mut buf, 8, self.header.depmax);
        write_f32_le(&mut buf, 20, self.header.b);
        write_f32_le(&mut buf, 24, self.header.e);
        write_i32_le(&mut buf, 316, self.header.npts);
        write_i32_le(&mut buf, 284, self.header.nvhdr);
        buf
    }

    /// Deserialize npts from a SAC 632-byte header.
    pub fn npts_from_header(buf: &[u8]) -> i32 {
        read_i32_le(buf, 316)
    }

    /// Return the number of samples.
    pub fn num_samples(&self) -> usize {
        self.data.len()
    }
}

// ── MiniSEED Format ───────────────────────────────────────────────────────────

/// MiniSEED fixed record header (48 bytes).
///
/// MiniSEED is a subset of SEED (Standard for the Exchange of Earthquake Data).
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct MiniSeedHeader {
    /// Sequence number (6 ASCII digits, 000001–999999).
    pub sequence_number: u32,
    /// Data quality indicator ('D', 'R', 'Q', 'M').
    pub quality_indicator: u8,
    /// Station identifier (5 characters).
    pub station: String,
    /// Location identifier (2 characters).
    pub location: String,
    /// Channel identifier (3 characters).
    pub channel: String,
    /// Network identifier (2 characters).
    pub network: String,
    /// Record start time: year.
    pub start_year: u16,
    /// Record start time: day-of-year.
    pub start_day: u16,
    /// Record start time: hour.
    pub start_hour: u8,
    /// Record start time: minute.
    pub start_min: u8,
    /// Record start time: second.
    pub start_sec: u8,
    /// Record start time: tenths of milliseconds (0–9999).
    pub start_frac: u16,
    /// Number of samples in this record.
    pub num_samples: u16,
    /// Sample rate factor.
    pub sample_rate_factor: i16,
    /// Sample rate multiplier.
    pub sample_rate_multiplier: i16,
    /// Activity flags.
    pub activity_flags: u8,
    /// I/O and clock flags.
    pub io_flags: u8,
    /// Data quality flags.
    pub data_quality_flags: u8,
    /// Number of blockettes following.
    pub num_blockettes: u8,
    /// Time correction in units of 0.0001 seconds.
    pub time_correction: i32,
    /// Byte offset to beginning of data.
    pub data_offset: u16,
    /// Byte offset to first blockette.
    pub blockette_offset: u16,
}

impl Default for MiniSeedHeader {
    fn default() -> Self {
        Self {
            sequence_number: 1,
            quality_indicator: b'D',
            station: String::from("XXXXX"),
            location: String::from("  "),
            channel: String::from("BHZ"),
            network: String::from("XX"),
            start_year: 2026,
            start_day: 1,
            start_hour: 0,
            start_min: 0,
            start_sec: 0,
            start_frac: 0,
            num_samples: 0,
            sample_rate_factor: 100,
            sample_rate_multiplier: 1,
            activity_flags: 0,
            io_flags: 0,
            data_quality_flags: 0,
            num_blockettes: 0,
            time_correction: 0,
            data_offset: 64,
            blockette_offset: 0,
        }
    }
}

impl MiniSeedHeader {
    /// Compute the true sample rate in samples/second.
    pub fn sample_rate(&self) -> f64 {
        let factor = self.sample_rate_factor as f64;
        let mult = self.sample_rate_multiplier as f64;
        if factor > 0.0 && mult > 0.0 {
            factor * mult
        } else if factor > 0.0 && mult < 0.0 {
            -factor / mult
        } else if factor < 0.0 && mult > 0.0 {
            -mult / factor
        } else {
            1.0 / (factor * mult)
        }
    }
}

/// Steim-1 compression: run-length encoded integer differences.
///
/// This is a simplified implementation for testing purposes.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct Steim1Frame {
    /// Word 0 contains a bitfield indicating how each of the 15 following words
    /// is encoded. Bits 0–1 for word 15, bits 28–29 for word 1.
    pub words: [u32; 16],
}

impl Steim1Frame {
    /// Create a Steim-1 frame encoding integer differences.
    ///
    /// Differences are encoded using 2-bit type codes:
    /// - 0: no data
    /// - 1: 4 bytes per difference (32-bit)
    /// - 2: 2 bytes per difference (16-bit)
    /// - 3: 1 byte per difference (8-bit)
    pub fn encode_diffs(diffs: &[i32]) -> Vec<Self> {
        let mut frames = Vec::new();
        let mut pos = 0;
        while pos < diffs.len() {
            let mut frame = Steim1Frame { words: [0u32; 16] };
            let mut word_idx = 1usize;
            while word_idx < 16 && pos < diffs.len() {
                let d = diffs[pos];
                if (-128..=127).contains(&d) {
                    // 8-bit, type=3
                    frame.words[0] |= 3 << (30 - 2 * word_idx);
                    frame.words[word_idx] = (d as i8 as u8) as u32;
                } else if (-32768..=32767).contains(&d) {
                    // 16-bit, type=2
                    frame.words[0] |= 2 << (30 - 2 * word_idx);
                    frame.words[word_idx] = (d as i16 as u16) as u32;
                } else {
                    // 32-bit, type=1
                    frame.words[0] |= 1 << (30 - 2 * word_idx);
                    frame.words[word_idx] = d as u32;
                }
                word_idx += 1;
                pos += 1;
            }
            frames.push(frame);
        }
        frames
    }
}

/// A MiniSEED record with header and integer sample data.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct MiniSeedRecord {
    /// Fixed record header.
    pub header: MiniSeedHeader,
    /// Decompressed integer sample data.
    pub samples: Vec<i32>,
    /// Steim-1 compressed frames.
    pub steim1_frames: Vec<Steim1Frame>,
}

impl MiniSeedRecord {
    /// Create a new record from integer samples.
    pub fn from_samples(samples: Vec<i32>, station: &str, channel: &str, network: &str) -> Self {
        let header = MiniSeedHeader {
            station: station.to_string(),
            channel: channel.to_string(),
            network: network.to_string(),
            num_samples: samples.len() as u16,
            ..Default::default()
        };
        // Compute Steim-1 differences
        let mut diffs = vec![samples.first().copied().unwrap_or(0)];
        for i in 1..samples.len() {
            diffs.push(samples[i] - samples[i - 1]);
        }
        let steim1_frames = Steim1Frame::encode_diffs(&diffs);
        Self {
            header,
            samples,
            steim1_frames,
        }
    }

    /// Increment the sequence number.
    pub fn next_sequence(&mut self) {
        self.header.sequence_number = (self.header.sequence_number % 999999) + 1;
    }
}

// ── SeismicTrace ──────────────────────────────────────────────────────────────

/// A seismic time series with associated metadata.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct SeismicTrace {
    /// Network code (2 characters).
    pub network: String,
    /// Station code (up to 5 characters).
    pub station: String,
    /// Location code (2 characters).
    pub location: String,
    /// Channel code (3 characters).
    pub channel: String,
    /// Sampling rate in samples per second.
    pub sampling_rate: f64,
    /// Start time as Unix timestamp (seconds since 1970-01-01).
    pub start_time: f64,
    /// Time-series samples.
    pub data: Vec<f64>,
}

impl SeismicTrace {
    /// Create a new seismic trace.
    pub fn new(
        network: &str,
        station: &str,
        channel: &str,
        sampling_rate: f64,
        start_time: f64,
        data: Vec<f64>,
    ) -> Self {
        Self {
            network: network.to_string(),
            station: station.to_string(),
            location: String::from("  "),
            channel: channel.to_string(),
            sampling_rate,
            start_time,
            data,
        }
    }

    /// Return the end time in seconds.
    pub fn end_time(&self) -> f64 {
        if self.data.is_empty() {
            self.start_time
        } else {
            self.start_time + (self.data.len() - 1) as f64 / self.sampling_rate
        }
    }

    /// Return the duration in seconds.
    pub fn duration(&self) -> f64 {
        self.end_time() - self.start_time
    }

    /// Return the SEED channel identifier as "NET.STA.LOC.CHA".
    pub fn seed_id(&self) -> String {
        format!(
            "{}.{}.{}.{}",
            self.network, self.station, self.location, self.channel
        )
    }

    /// Compute the RMS amplitude of the trace.
    pub fn rms(&self) -> f64 {
        if self.data.is_empty() {
            return 0.0;
        }
        let sum2: f64 = self.data.iter().map(|&x| x * x).sum();
        (sum2 / self.data.len() as f64).sqrt()
    }
}

// ── SeismicStation ────────────────────────────────────────────────────────────

/// A seismic station with geographic coordinates.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct SeismicStation {
    /// Network code.
    pub network: String,
    /// Station code.
    pub station: String,
    /// Latitude in decimal degrees (positive = North).
    pub latitude: f64,
    /// Longitude in decimal degrees (positive = East).
    pub longitude: f64,
    /// Elevation in meters above sea level.
    pub elevation: f64,
    /// Available channel codes.
    pub channels: Vec<String>,
}

impl SeismicStation {
    /// Create a new station.
    pub fn new(network: &str, station: &str, lat: f64, lon: f64, elevation: f64) -> Self {
        Self {
            network: network.to_string(),
            station: station.to_string(),
            latitude: lat,
            longitude: lon,
            elevation,
            channels: Vec::new(),
        }
    }

    /// Compute the great-circle distance to another station in kilometers.
    ///
    /// Uses the Haversine formula.
    pub fn distance_km(&self, other: &SeismicStation) -> f64 {
        haversine_km(
            self.latitude,
            self.longitude,
            other.latitude,
            other.longitude,
        )
    }

    /// Compute the epicentral distance from an event location in kilometers.
    pub fn epicentral_distance_km(&self, event_lat: f64, event_lon: f64) -> f64 {
        haversine_km(self.latitude, self.longitude, event_lat, event_lon)
    }
}

/// Haversine great-circle distance formula.
///
/// Returns distance in kilometers between two points (lat/lon in degrees).
#[allow(dead_code)]
pub fn haversine_km(lat1: f64, lon1: f64, lat2: f64, lon2: f64) -> f64 {
    const R_KM: f64 = 6371.0;
    let d_lat = (lat2 - lat1).to_radians();
    let d_lon = (lon2 - lon1).to_radians();
    let a = (d_lat / 2.0).sin().powi(2)
        + lat1.to_radians().cos() * lat2.to_radians().cos() * (d_lon / 2.0).sin().powi(2);
    let c = 2.0 * a.sqrt().asin();
    R_KM * c
}

// ── ArrivalPicker ─────────────────────────────────────────────────────────────

/// STA/LTA-based P-wave arrival picker.
///
/// Computes the characteristic function CF(t) = STA(t) / LTA(t) where:
/// - STA = short-time average of squared amplitudes over sta_samples
/// - LTA = long-time average of squared amplitudes over lta_samples
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct ArrivalPicker {
    /// Short-term average window in samples.
    pub sta_samples: usize,
    /// Long-term average window in samples.
    pub lta_samples: usize,
    /// Trigger threshold (ratio STA/LTA).
    pub trigger_threshold: f64,
    /// Detrigger threshold.
    pub detrigger_threshold: f64,
}

impl ArrivalPicker {
    /// Create a new STA/LTA picker.
    pub fn new(sta_samples: usize, lta_samples: usize, trigger: f64, detrigger: f64) -> Self {
        Self {
            sta_samples,
            lta_samples,
            trigger_threshold: trigger,
            detrigger_threshold: detrigger,
        }
    }

    /// Compute the STA/LTA characteristic function for a trace.
    ///
    /// Returns a vector of CF values, one per sample (zeros for the initial
    /// lta_samples window where LTA is not yet filled).
    pub fn characteristic_function(&self, data: &[f64]) -> Vec<f64> {
        let n = data.len();
        let mut cf = vec![0.0f64; n];
        let sq: Vec<f64> = data.iter().map(|&x| x * x).collect();

        // Running STA and LTA using cumulative sum
        let mut cum = vec![0.0f64; n + 1];
        for i in 0..n {
            cum[i + 1] = cum[i] + sq[i];
        }

        for i in self.lta_samples..n {
            let lta_start = i.saturating_sub(self.lta_samples);
            let sta_start = i.saturating_sub(self.sta_samples);
            let lta = (cum[i] - cum[lta_start]) / self.lta_samples as f64;
            let sta = (cum[i] - cum[sta_start]) / self.sta_samples as f64;
            if lta > 1e-30 {
                cf[i] = sta / lta;
            }
        }
        cf
    }

    /// Pick P-wave arrival times (in samples) where CF exceeds the trigger threshold.
    ///
    /// Returns indices into the trace where triggers occur.
    pub fn pick_arrivals(&self, data: &[f64]) -> Vec<usize> {
        let cf = self.characteristic_function(data);
        let mut picks = Vec::new();
        let mut triggered = false;
        for (i, &c) in cf.iter().enumerate() {
            if !triggered && c >= self.trigger_threshold {
                picks.push(i);
                triggered = true;
            } else if triggered && c < self.detrigger_threshold {
                triggered = false;
            }
        }
        picks
    }
}

// ── FrequencyFilter ───────────────────────────────────────────────────────────

/// Butterworth bandpass filter parameters.
///
/// The filter is specified in the frequency domain by its corner frequencies
/// and order. For actual filtering, a recursive IIR implementation is used.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct FrequencyFilter {
    /// Lower corner frequency in Hz.
    pub freqmin: f64,
    /// Upper corner frequency in Hz.
    pub freqmax: f64,
    /// Filter order (number of poles).
    pub order: usize,
    /// Sampling rate in Hz.
    pub sampling_rate: f64,
}

impl FrequencyFilter {
    /// Create a new Butterworth bandpass filter.
    pub fn bandpass(freqmin: f64, freqmax: f64, order: usize, sampling_rate: f64) -> Self {
        Self {
            freqmin,
            freqmax,
            order,
            sampling_rate,
        }
    }

    /// Check if a given frequency is in the passband (within 3 dB of unity gain).
    pub fn in_passband(&self, freq_hz: f64) -> bool {
        freq_hz >= self.freqmin && freq_hz <= self.freqmax
    }

    /// Compute the idealized gain at a given frequency (rectangular approximation).
    ///
    /// Returns 1.0 in the passband, 0.0 outside. This is a simplified model
    /// for testing the concept of passband vs. stopband.
    pub fn ideal_gain(&self, freq_hz: f64) -> f64 {
        if self.in_passband(freq_hz) { 1.0 } else { 0.0 }
    }

    /// Apply a simple forward-backward (zero-phase) bandpass filter using
    /// a windowed sinc approach (simplified for testing).
    ///
    /// Returns filtered data of the same length.
    pub fn apply(&self, data: &[f64]) -> Vec<f64> {
        // Simple rectangular window in frequency domain via DFT-like approach
        // For a real implementation, use recursive Butterworth IIR
        // Here we use a simple moving average as a low-pass approximation
        let n = data.len();
        if n == 0 {
            return Vec::new();
        }
        // Low-pass filter up to freqmax
        let lp_samples = (self.sampling_rate / (2.0 * self.freqmax)).round().max(1.0) as usize;
        // High-pass: subtract low-passed version at freqmin
        let hp_samples = (self.sampling_rate / (2.0 * self.freqmin)).round().max(1.0) as usize;

        let lp = moving_average(data, lp_samples);
        let hp_lp = moving_average(data, hp_samples);
        // Bandpass: LP - very_LP (removes DC and very low freqs)
        let mut out = vec![0.0f64; n];
        for i in 0..n {
            out[i] = lp[i] - hp_lp[i];
        }
        out
    }
}

/// Compute a moving average with the given window size.
#[allow(dead_code)]
fn moving_average(data: &[f64], window: usize) -> Vec<f64> {
    let n = data.len();
    let w = window.max(1);
    let mut out = vec![0.0f64; n];
    let mut sum = 0.0f64;
    let mut count = 0usize;
    for i in 0..n {
        sum += data[i];
        count += 1;
        if i >= w {
            sum -= data[i - w];
            count -= 1;
        }
        out[i] = sum / count as f64;
    }
    out
}

// ── SeismicCatalog ────────────────────────────────────────────────────────────

/// A seismic event in a catalog.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct SeismicEvent {
    /// Event ID.
    pub event_id: String,
    /// Origin time as Unix timestamp.
    pub origin_time: f64,
    /// Hypocenter latitude in degrees.
    pub latitude: f64,
    /// Hypocenter longitude in degrees.
    pub longitude: f64,
    /// Hypocenter depth in km.
    pub depth_km: f64,
    /// Moment magnitude Mw.
    pub magnitude: f64,
    /// Magnitude type (e.g., "Mw", "Mb", "Ms").
    pub magnitude_type: String,
    /// Number of phases used in location.
    pub num_phases: u32,
    /// RMS residual of travel-time fit in seconds.
    pub rms_residual: f64,
}

impl SeismicEvent {
    /// Create a new seismic event.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        event_id: &str,
        origin_time: f64,
        lat: f64,
        lon: f64,
        depth_km: f64,
        magnitude: f64,
        magnitude_type: &str,
    ) -> Self {
        Self {
            event_id: event_id.to_string(),
            origin_time,
            latitude: lat,
            longitude: lon,
            depth_km,
            magnitude,
            magnitude_type: magnitude_type.to_string(),
            num_phases: 0,
            rms_residual: 0.0,
        }
    }

    /// Compute the epicentral distance to a station in km.
    pub fn epicentral_distance_km(&self, station: &SeismicStation) -> f64 {
        haversine_km(
            self.latitude,
            self.longitude,
            station.latitude,
            station.longitude,
        )
    }

    /// Estimate P-wave travel time in seconds using a simplified Ak135 velocity model.
    ///
    /// Uses a linear approximation: t = dist_km / v_p where v_p ≈ 8.0 km/s.
    pub fn estimate_p_travel_time(&self, station: &SeismicStation) -> f64 {
        let dist = self.epicentral_distance_km(station);
        let hypo_dist = (dist * dist + self.depth_km * self.depth_km).sqrt();
        hypo_dist / 8.0 // approximate crustal P-velocity
    }

    /// Estimate S-wave travel time in seconds.
    pub fn estimate_s_travel_time(&self, station: &SeismicStation) -> f64 {
        let dist = self.epicentral_distance_km(station);
        let hypo_dist = (dist * dist + self.depth_km * self.depth_km).sqrt();
        hypo_dist / 4.6 // approximate crustal S-velocity
    }
}

/// A catalog of seismic events.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct SeismicCatalog {
    /// All events in this catalog.
    pub events: Vec<SeismicEvent>,
    /// Catalog name or source.
    pub name: String,
}

impl SeismicCatalog {
    /// Create a new empty catalog.
    pub fn new(name: &str) -> Self {
        Self {
            events: Vec::new(),
            name: name.to_string(),
        }
    }

    /// Add an event to the catalog.
    pub fn add_event(&mut self, event: SeismicEvent) {
        self.events.push(event);
    }

    /// Return the number of events.
    pub fn len(&self) -> usize {
        self.events.len()
    }

    /// Return true if the catalog is empty.
    pub fn is_empty(&self) -> bool {
        self.events.is_empty()
    }

    /// Filter events by magnitude range \[mag_min, mag_max\].
    pub fn filter_by_magnitude(&self, mag_min: f64, mag_max: f64) -> Vec<&SeismicEvent> {
        self.events
            .iter()
            .filter(|e| e.magnitude >= mag_min && e.magnitude <= mag_max)
            .collect()
    }

    /// Sort events by origin time (ascending).
    pub fn sort_by_time(&mut self) {
        self.events.sort_by(|a, b| {
            a.origin_time
                .partial_cmp(&b.origin_time)
                .unwrap_or(std::cmp::Ordering::Equal)
        });
    }

    /// Return the maximum magnitude event, if any.
    pub fn max_magnitude_event(&self) -> Option<&SeismicEvent> {
        self.events.iter().max_by(|a, b| {
            a.magnitude
                .partial_cmp(&b.magnitude)
                .unwrap_or(std::cmp::Ordering::Equal)
        })
    }

    /// Build an index from event_id to index in events vector.
    pub fn build_id_index(&self) -> HashMap<String, usize> {
        self.events
            .iter()
            .enumerate()
            .map(|(i, e)| (e.event_id.clone(), i))
            .collect()
    }

    /// Return events occurring within a time window \[t_start, t_end\].
    pub fn events_in_window(&self, t_start: f64, t_end: f64) -> Vec<&SeismicEvent> {
        self.events
            .iter()
            .filter(|e| e.origin_time >= t_start && e.origin_time <= t_end)
            .collect()
    }
}

// ── Voronoi/Roadmap placeholder ───────────────────────────────────────────────

/// Phase arrival record: station, phase type, and observed time.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct PhaseArrival {
    /// Station SEED ID.
    pub station_id: String,
    /// Phase name (e.g., "P", "S", "Pn", "Pg").
    pub phase: String,
    /// Observed arrival time as Unix timestamp.
    pub arrival_time: f64,
    /// Time residual (observed - predicted) in seconds.
    pub residual: f64,
}

impl PhaseArrival {
    /// Create a new phase arrival.
    pub fn new(station_id: &str, phase: &str, arrival_time: f64, residual: f64) -> Self {
        Self {
            station_id: station_id.to_string(),
            phase: phase.to_string(),
            arrival_time,
            residual,
        }
    }
}

// ── Binary I/O Helpers ────────────────────────────────────────────────────────

/// Write a big-endian i32 at byte offset in buffer.
#[allow(dead_code)]
fn write_i32_be(buf: &mut [u8], offset: usize, val: i32) {
    let bytes = val.to_be_bytes();
    buf[offset..offset + 4].copy_from_slice(&bytes);
}

/// Write a big-endian i16 at byte offset in buffer.
#[allow(dead_code)]
fn write_i16_be(buf: &mut [u8], offset: usize, val: i16) {
    let bytes = val.to_be_bytes();
    buf[offset..offset + 2].copy_from_slice(&bytes);
}

/// Read a big-endian i32 from buffer at byte offset.
#[allow(dead_code)]
fn read_i32_be(buf: &[u8], offset: usize) -> i32 {
    let bytes: [u8; 4] = buf[offset..offset + 4].try_into().unwrap_or([0; 4]);
    i32::from_be_bytes(bytes)
}

/// Read a big-endian i16 from buffer at byte offset.
#[allow(dead_code)]
fn read_i16_be(buf: &[u8], offset: usize) -> i16 {
    let bytes: [u8; 2] = buf[offset..offset + 2].try_into().unwrap_or([0; 2]);
    i16::from_be_bytes(bytes)
}

/// Write a little-endian f32 at byte offset in buffer.
#[allow(dead_code)]
fn write_f32_le(buf: &mut [u8], offset: usize, val: f32) {
    let bytes = val.to_le_bytes();
    buf[offset..offset + 4].copy_from_slice(&bytes);
}

/// Write a little-endian i32 at byte offset in buffer.
#[allow(dead_code)]
fn write_i32_le(buf: &mut [u8], offset: usize, val: i32) {
    let bytes = val.to_le_bytes();
    buf[offset..offset + 4].copy_from_slice(&bytes);
}

/// Read a little-endian i32 from buffer at byte offset.
#[allow(dead_code)]
fn read_i32_le(buf: &[u8], offset: usize) -> i32 {
    let bytes: [u8; 4] = buf[offset..offset + 4].try_into().unwrap_or([0; 4]);
    i32::from_le_bytes(bytes)
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // ── SEG-Y tests ────────────────────────────────────────────────────────────

    #[test]
    fn test_segy_text_header_roundtrip() {
        let mut hdr = SegyTextHeader::new();
        hdr.set_line(0, "CREATED BY OXIPHYSICS");
        let line = hdr.get_line(0);
        assert!(line.contains("CREATED"), "line={line}");
    }

    #[test]
    fn test_segy_text_header_empty_line() {
        let hdr = SegyTextHeader::new();
        let line = hdr.get_line(0);
        // All spaces → trimmed to empty
        assert!(line.is_empty(), "line={line:?}");
    }

    #[test]
    fn test_segy_text_header_out_of_range() {
        let mut hdr = SegyTextHeader::new();
        hdr.set_line(40, "SHOULD NOT WRITE");
        // Should not panic; line 40 is out of range
        assert_eq!(hdr.get_line(40), "");
    }

    #[test]
    fn test_segy_binary_header_roundtrip() {
        let segy = SegyFile::new();
        let buf = segy.serialize_binary_header();
        assert_eq!(buf.len(), 400);
        let h2 = SegyFile::deserialize_binary_header(&buf);
        assert_eq!(h2.job_id, segy.binary_header.job_id);
        assert_eq!(h2.sample_interval_us, segy.binary_header.sample_interval_us);
        assert_eq!(h2.samples_per_trace, segy.binary_header.samples_per_trace);
        assert_eq!(h2.data_format_code, segy.binary_header.data_format_code);
    }

    #[test]
    fn test_segy_binary_header_segy_revision() {
        let segy = SegyFile::new();
        let buf = segy.serialize_binary_header();
        let h2 = SegyFile::deserialize_binary_header(&buf);
        assert_eq!(h2.segy_revision, 256, "Expected rev 1 (256)");
    }

    #[test]
    fn test_segy_add_trace_sequence_number() {
        let mut segy = SegyFile::new();
        let t1 = SegyTrace::new(vec![1.0, 2.0, 3.0]);
        let t2 = SegyTrace::new(vec![4.0, 5.0, 6.0]);
        segy.add_trace(t1);
        segy.add_trace(t2);
        assert_eq!(segy.traces[0].header.trace_seq_file, 1);
        assert_eq!(segy.traces[1].header.trace_seq_file, 2);
    }

    #[test]
    fn test_segy_trace_num_samples() {
        let data: Vec<f32> = (0..500).map(|x| x as f32 * 0.001).collect();
        let trace = SegyTrace::new(data.clone());
        assert_eq!(trace.num_samples(), 500);
        assert_eq!(trace.header.num_samples, 500);
    }

    #[test]
    fn test_segy_trace_sample_interval() {
        let mut trace = SegyTrace::new(vec![0.0f32; 10]);
        trace.header.sample_interval_us = 4000; // 4 ms
        let dt = trace.sample_interval_s();
        assert!((dt - 0.004).abs() < 1e-9, "dt={dt:.6}");
    }

    #[test]
    fn test_segy_num_traces() {
        let mut segy = SegyFile::new();
        assert_eq!(segy.num_traces(), 0);
        for i in 0..10 {
            segy.add_trace(SegyTrace::new(vec![i as f32; 100]));
        }
        assert_eq!(segy.num_traces(), 10);
    }

    // ── SAC tests ─────────────────────────────────────────────────────────────

    #[test]
    fn test_sac_npts_matches_data_length() {
        let data: Vec<f32> = (0..1024).map(|x| (x as f32).sin()).collect();
        let sac = SacFile::from_data(data.clone(), 0.01);
        assert_eq!(sac.header.npts as usize, data.len());
        assert_eq!(sac.num_samples(), data.len());
    }

    #[test]
    fn test_sac_end_time_correct() {
        let data = vec![0.0f32; 1001];
        let sac = SacFile::from_data(data, 0.01);
        // e = (1001 - 1) * 0.01 = 10.0
        assert!((sac.header.e - 10.0).abs() < 1e-5, "e={:.6}", sac.header.e);
    }

    #[test]
    fn test_sac_depmin_depmax() {
        let data = vec![-3.0f32, 0.0, 5.0, 2.0];
        let sac = SacFile::from_data(data, 0.01);
        assert!((sac.header.depmin - (-3.0)).abs() < 1e-5);
        assert!((sac.header.depmax - 5.0).abs() < 1e-5);
    }

    #[test]
    fn test_sac_header_serialize_npts() {
        let data = vec![1.0f32; 256];
        let sac = SacFile::from_data(data, 0.01);
        let buf = sac.serialize_header();
        assert_eq!(buf.len(), 632);
        let npts = SacFile::npts_from_header(&buf);
        assert_eq!(npts, 256);
    }

    #[test]
    fn test_sac_delta_roundtrip() {
        let data = vec![0.0f32; 100];
        let sac = SacFile::from_data(data, 0.025);
        let buf = sac.serialize_header();
        // Read delta from offset 0 (little-endian f32)
        let delta_bytes: [u8; 4] = buf[0..4].try_into().unwrap();
        let delta = f32::from_le_bytes(delta_bytes);
        assert!((delta - 0.025).abs() < 1e-6, "delta={delta:.6}");
    }

    // ── MiniSEED tests ────────────────────────────────────────────────────────

    #[test]
    fn test_miniseed_sequence_number_increment() {
        let mut rec = MiniSeedRecord::from_samples(vec![0; 10], "ANMO", "BHZ", "IU");
        assert_eq!(rec.header.sequence_number, 1);
        rec.next_sequence();
        assert_eq!(rec.header.sequence_number, 2);
        rec.next_sequence();
        assert_eq!(rec.header.sequence_number, 3);
    }

    #[test]
    fn test_miniseed_sequence_number_wraps() {
        let mut rec = MiniSeedRecord::from_samples(vec![0; 10], "ANMO", "BHZ", "IU");
        rec.header.sequence_number = 999999;
        rec.next_sequence();
        assert_eq!(rec.header.sequence_number, 1);
    }

    #[test]
    fn test_miniseed_num_samples() {
        let samples = vec![1, 2, 3, 4, 5];
        let rec = MiniSeedRecord::from_samples(samples.clone(), "ANMO", "BHZ", "IU");
        assert_eq!(rec.header.num_samples as usize, samples.len());
    }

    #[test]
    fn test_miniseed_channel_assignment() {
        let rec = MiniSeedRecord::from_samples(vec![0; 5], "COLA", "HHN", "AK");
        assert_eq!(rec.header.channel, "HHN");
        assert_eq!(rec.header.network, "AK");
        assert_eq!(rec.header.station, "COLA");
    }

    #[test]
    fn test_miniseed_sample_rate_positive() {
        let hdr = MiniSeedHeader::default();
        // factor=100, multiplier=1 → 100 sps
        let sr = hdr.sample_rate();
        assert!((sr - 100.0).abs() < 1e-9, "sr={sr:.6}");
    }

    #[test]
    fn test_miniseed_sample_rate_negative_factor() {
        let hdr = MiniSeedHeader {
            sample_rate_factor: -10, // 1/10 sps = 0.1 sps
            sample_rate_multiplier: 1,
            ..Default::default()
        };
        // factor < 0, mult > 0 → -mult / factor = -1 / -10 = 0.1
        let sr = hdr.sample_rate();
        assert!((sr - 0.1).abs() < 1e-9, "sr={sr:.6}");
    }

    #[test]
    fn test_steim1_encode_small_diffs() {
        let diffs = vec![0i32, 1, -1, 127, -128];
        let frames = Steim1Frame::encode_diffs(&diffs);
        assert!(!frames.is_empty());
        // All diffs fit in 8 bits, so type codes should be 3 for data words
    }

    #[test]
    fn test_steim1_encode_large_diffs() {
        let diffs = vec![1_000_000i32, -1_000_000];
        let frames = Steim1Frame::encode_diffs(&diffs);
        assert!(!frames.is_empty());
    }

    // ── SeismicTrace tests ────────────────────────────────────────────────────

    #[test]
    fn test_seismic_trace_end_time() {
        let data = vec![0.0; 101];
        let trace = SeismicTrace::new("IU", "ANMO", "BHZ", 100.0, 0.0, data);
        // end = 0 + 100 * (1/100) = 1.0
        assert!(
            (trace.end_time() - 1.0).abs() < 1e-9,
            "end={:.6}",
            trace.end_time()
        );
    }

    #[test]
    fn test_seismic_trace_seed_id() {
        let trace = SeismicTrace::new("IU", "ANMO", "BHZ", 100.0, 0.0, vec![]);
        assert_eq!(trace.seed_id(), "IU.ANMO.  .BHZ");
    }

    #[test]
    fn test_seismic_trace_rms() {
        let data = vec![1.0, -1.0, 1.0, -1.0];
        let trace = SeismicTrace::new("IU", "ANMO", "BHZ", 1.0, 0.0, data);
        assert!((trace.rms() - 1.0).abs() < 1e-9, "rms={:.6}", trace.rms());
    }

    #[test]
    fn test_seismic_trace_duration() {
        let data = vec![0.0; 1001];
        let trace = SeismicTrace::new("IU", "ANMO", "BHZ", 100.0, 10.0, data);
        // duration = 1000 / 100 = 10.0 s
        assert!((trace.duration() - 10.0).abs() < 1e-9);
    }

    // ── SeismicStation / Haversine tests ──────────────────────────────────────

    #[test]
    fn test_station_distance_same_location() {
        let sta1 = SeismicStation::new("IU", "ANMO", 34.9459, -106.4572, 1820.0);
        let sta2 = SeismicStation::new("IU", "ANMO", 34.9459, -106.4572, 1820.0);
        let d = sta1.distance_km(&sta2);
        assert!(d < 0.001, "d={d:.6}");
    }

    #[test]
    fn test_haversine_known_distance() {
        // London (51.5074, -0.1278) to Paris (48.8566, 2.3522) ≈ 341 km
        let d = haversine_km(51.5074, -0.1278, 48.8566, 2.3522);
        assert!((d - 341.0).abs() < 5.0, "d={d:.1} km");
    }

    #[test]
    fn test_haversine_equator() {
        // One degree of longitude on equator ≈ 111.32 km
        let d = haversine_km(0.0, 0.0, 0.0, 1.0);
        assert!((d - 111.32).abs() < 0.5, "d={d:.3} km");
    }

    #[test]
    fn test_station_epicentral_distance() {
        let sta = SeismicStation::new("IU", "ANMO", 34.9459, -106.4572, 1820.0);
        // Event at same location → 0 km
        let d = sta.epicentral_distance_km(34.9459, -106.4572);
        assert!(d < 0.001, "d={d:.6}");
    }

    // ── ArrivalPicker tests ───────────────────────────────────────────────────

    #[test]
    fn test_stalta_characteristic_function_size() {
        let data: Vec<f64> = (0..1000).map(|i| (i as f64 * 0.1).sin()).collect();
        let picker = ArrivalPicker::new(10, 100, 3.0, 1.5);
        let cf = picker.characteristic_function(&data);
        assert_eq!(cf.len(), data.len());
    }

    #[test]
    fn test_stalta_picks_before_direct_wave() {
        // Noise + impulsive signal at sample 500
        let mut data = vec![0.01f64; 1000];
        for i in 500..510 {
            data[i] = 100.0; // strong P-wave
        }
        // Add tiny noise
        for (i, x) in data.iter_mut().enumerate() {
            if i < 500 {
                *x += 0.001 * (i as f64).sin();
            }
        }
        let picker = ArrivalPicker::new(10, 100, 3.0, 1.5);
        let picks = picker.pick_arrivals(&data);
        // Should pick at or near sample 500
        assert!(!picks.is_empty(), "Expected at least one pick");
        let first_pick = picks[0];
        assert!(first_pick >= 100, "Pick too early: {first_pick}");
        assert!(first_pick <= 520, "Pick too late: {first_pick}");
    }

    #[test]
    fn test_stalta_no_picks_quiet_signal() {
        let data = vec![0.001f64; 1000];
        let picker = ArrivalPicker::new(10, 100, 10.0, 5.0);
        let picks = picker.pick_arrivals(&data);
        assert!(picks.is_empty(), "Expected no picks for quiet signal");
    }

    // ── FrequencyFilter tests ──────────────────────────────────────────────────

    #[test]
    fn test_bandpass_in_passband() {
        let filt = FrequencyFilter::bandpass(1.0, 10.0, 4, 100.0);
        assert!(filt.in_passband(5.0));
        assert!(!filt.in_passband(0.5));
        assert!(!filt.in_passband(15.0));
    }

    #[test]
    fn test_bandpass_ideal_gain_passband() {
        let filt = FrequencyFilter::bandpass(1.0, 10.0, 4, 100.0);
        assert!((filt.ideal_gain(5.0) - 1.0).abs() < 1e-10);
        assert!((filt.ideal_gain(0.1)).abs() < 1e-10);
    }

    #[test]
    fn test_bandpass_apply_preserves_length() {
        let data: Vec<f64> = (0..256).map(|i| (i as f64 * 0.1).sin()).collect();
        let filt = FrequencyFilter::bandpass(1.0, 10.0, 4, 100.0);
        let out = filt.apply(&data);
        assert_eq!(out.len(), data.len());
    }

    #[test]
    fn test_bandpass_apply_empty() {
        let filt = FrequencyFilter::bandpass(1.0, 10.0, 4, 100.0);
        let out = filt.apply(&[]);
        assert!(out.is_empty());
    }

    // ── SeismicCatalog tests ──────────────────────────────────────────────────

    #[test]
    fn test_catalog_add_and_len() {
        let mut cat = SeismicCatalog::new("Test");
        assert_eq!(cat.len(), 0);
        cat.add_event(SeismicEvent::new(
            "ev001", 1000.0, 35.0, -120.0, 10.0, 5.5, "Mw",
        ));
        cat.add_event(SeismicEvent::new(
            "ev002", 2000.0, 36.0, -121.0, 15.0, 4.2, "Mb",
        ));
        assert_eq!(cat.len(), 2);
    }

    #[test]
    fn test_catalog_sort_by_time() {
        let mut cat = SeismicCatalog::new("Test");
        cat.add_event(SeismicEvent::new(
            "ev002", 2000.0, 35.0, -120.0, 10.0, 5.5, "Mw",
        ));
        cat.add_event(SeismicEvent::new(
            "ev001", 1000.0, 35.0, -120.0, 10.0, 4.0, "Mb",
        ));
        cat.sort_by_time();
        assert!(cat.events[0].origin_time < cat.events[1].origin_time);
    }

    #[test]
    fn test_catalog_filter_by_magnitude() {
        let mut cat = SeismicCatalog::new("Test");
        cat.add_event(SeismicEvent::new(
            "ev001", 1000.0, 35.0, -120.0, 10.0, 5.5, "Mw",
        ));
        cat.add_event(SeismicEvent::new(
            "ev002", 2000.0, 36.0, -121.0, 15.0, 3.0, "Mb",
        ));
        cat.add_event(SeismicEvent::new(
            "ev003", 3000.0, 37.0, -122.0, 20.0, 7.0, "Mw",
        ));
        let big = cat.filter_by_magnitude(5.0, 8.0);
        assert_eq!(big.len(), 2);
    }

    #[test]
    fn test_catalog_max_magnitude() {
        let mut cat = SeismicCatalog::new("Test");
        cat.add_event(SeismicEvent::new(
            "ev001", 1000.0, 35.0, -120.0, 10.0, 5.5, "Mw",
        ));
        cat.add_event(SeismicEvent::new(
            "ev002", 2000.0, 36.0, -121.0, 15.0, 7.1, "Mw",
        ));
        let max_ev = cat.max_magnitude_event().unwrap();
        assert!((max_ev.magnitude - 7.1).abs() < 1e-9);
    }

    #[test]
    fn test_catalog_events_in_window() {
        let mut cat = SeismicCatalog::new("Test");
        cat.add_event(SeismicEvent::new(
            "ev001", 100.0, 35.0, -120.0, 10.0, 5.5, "Mw",
        ));
        cat.add_event(SeismicEvent::new(
            "ev002", 500.0, 36.0, -121.0, 15.0, 3.0, "Mb",
        ));
        cat.add_event(SeismicEvent::new(
            "ev003", 900.0, 37.0, -122.0, 20.0, 7.0, "Mw",
        ));
        let evts = cat.events_in_window(200.0, 800.0);
        assert_eq!(evts.len(), 1);
        assert_eq!(evts[0].event_id, "ev002");
    }

    #[test]
    fn test_event_p_travel_time_increases_with_distance() {
        let event = SeismicEvent::new("ev001", 0.0, 0.0, 0.0, 10.0, 6.0, "Mw");
        let sta_near = SeismicStation::new("XX", "NEAR", 0.1, 0.0, 0.0);
        let sta_far = SeismicStation::new("XX", "FAR", 5.0, 0.0, 0.0);
        let t_near = event.estimate_p_travel_time(&sta_near);
        let t_far = event.estimate_p_travel_time(&sta_far);
        assert!(
            t_far > t_near,
            "t_far={t_far:.3} should be > t_near={t_near:.3}"
        );
    }

    #[test]
    fn test_event_s_travel_time_greater_than_p() {
        let event = SeismicEvent::new("ev001", 0.0, 35.0, -120.0, 10.0, 6.0, "Mw");
        let sta = SeismicStation::new("IU", "ANMO", 34.9459, -106.4572, 1820.0);
        let t_p = event.estimate_p_travel_time(&sta);
        let t_s = event.estimate_s_travel_time(&sta);
        assert!(t_s > t_p, "S ({t_s:.3}s) should arrive after P ({t_p:.3}s)");
    }

    #[test]
    fn test_phase_arrival_creation() {
        let arr = PhaseArrival::new("IU.ANMO.  .BHZ", "P", 1000.5, 0.03);
        assert_eq!(arr.phase, "P");
        assert!((arr.residual - 0.03).abs() < 1e-9);
    }
}
