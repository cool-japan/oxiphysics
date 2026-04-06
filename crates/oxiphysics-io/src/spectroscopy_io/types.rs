//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::needless_range_loop)]
use std::collections::HashMap;
use std::io::{BufRead, BufReader, Read, Write};

use super::functions::ANIF_MAGIC;
#[allow(unused_imports)]
use super::functions::*;

/// Result of a peak-finding operation.
#[derive(Debug, Clone)]
pub struct PeakResult {
    /// Index in the spectrum array.
    pub index: usize,
    /// X position of the peak.
    pub x: f64,
    /// Y (intensity) of the peak.
    pub y: f64,
    /// Estimated full-width at half-maximum (in X units).
    pub fwhm: f64,
}
/// Endianness of the ANIF file.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum AnifEndian {
    /// Little-endian byte order.
    Little,
    /// Big-endian byte order.
    Big,
}
/// In-memory spectral database with similarity search.
///
/// Spectra are resampled to a common X grid before comparison.
#[derive(Debug, Default)]
pub struct SpectralDatabase {
    pub(super) entries: Vec<DatabaseEntry>,
    pub(super) next_id: usize,
}
impl SpectralDatabase {
    /// Create an empty database.
    pub fn new() -> Self {
        Self::default()
    }
    /// Insert a spectrum and return its assigned ID.
    pub fn insert(&mut self, record: SpectrumRecord) -> usize {
        let id = self.next_id;
        self.entries.push(DatabaseEntry { id, record });
        self.next_id += 1;
        id
    }
    /// Retrieve a spectrum by ID.
    pub fn get(&self, id: usize) -> Option<&SpectrumRecord> {
        self.entries.iter().find(|e| e.id == id).map(|e| &e.record)
    }
    /// Remove a spectrum by ID. Returns `true` if it was present.
    pub fn remove(&mut self, id: usize) -> bool {
        let before = self.entries.len();
        self.entries.retain(|e| e.id != id);
        self.entries.len() < before
    }
    /// Total number of stored spectra.
    pub fn len(&self) -> usize {
        self.entries.len()
    }
    /// Return `true` if the database is empty.
    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }
    /// Search for the `top_k` most similar spectra to `query`.
    ///
    /// Spectra are compared over their raw Y arrays (resampled if needed via
    /// `resample_to`). If the arrays differ in length, the shorter is zero-padded.
    pub fn search(
        &self,
        query: &SpectrumRecord,
        metric: SimilarityMetric,
        top_k: usize,
    ) -> Vec<SearchResult> {
        let mut scores: Vec<(usize, &str, f64)> = self
            .entries
            .iter()
            .map(|e| {
                let score = Self::compute_similarity(query, &e.record, metric);
                (e.id, e.record.metadata.title.as_str(), score)
            })
            .collect();
        match metric {
            SimilarityMetric::Euclidean => {
                scores.sort_by(|a, b| a.2.partial_cmp(&b.2).unwrap_or(std::cmp::Ordering::Equal));
            }
            _ => {
                scores.sort_by(|a, b| b.2.partial_cmp(&a.2).unwrap_or(std::cmp::Ordering::Equal));
            }
        }
        scores
            .into_iter()
            .take(top_k)
            .map(|(id, title, score)| SearchResult {
                id,
                title: title.to_string(),
                score,
            })
            .collect()
    }
    /// Compute similarity between two spectra.
    fn compute_similarity(a: &SpectrumRecord, b: &SpectrumRecord, metric: SimilarityMetric) -> f64 {
        let len = a.y.len().min(b.y.len());
        if len == 0 {
            return 0.0;
        }
        let ya: &[f64] = &a.y[..len];
        let yb: &[f64] = &b.y[..len];
        match metric {
            SimilarityMetric::DotProduct => {
                let dot: f64 = ya.iter().zip(yb).map(|(x, y)| x * y).sum();
                let na: f64 = ya.iter().map(|v| v * v).sum::<f64>().sqrt();
                let nb: f64 = yb.iter().map(|v| v * v).sum::<f64>().sqrt();
                if na > 0.0 && nb > 0.0 {
                    dot / (na * nb)
                } else {
                    0.0
                }
            }
            SimilarityMetric::SpectralAngleMapper => {
                let dot: f64 = ya.iter().zip(yb).map(|(x, y)| x * y).sum();
                let na: f64 = ya.iter().map(|v| v * v).sum::<f64>().sqrt();
                let nb: f64 = yb.iter().map(|v| v * v).sum::<f64>().sqrt();
                if na > 0.0 && nb > 0.0 {
                    let cos = (dot / (na * nb)).clamp(-1.0, 1.0);
                    std::f64::consts::PI / 2.0 - cos.acos()
                } else {
                    0.0
                }
            }
            SimilarityMetric::Euclidean => ya
                .iter()
                .zip(yb)
                .map(|(x, y)| (x - y).powi(2))
                .sum::<f64>()
                .sqrt(),
            SimilarityMetric::Pearson => {
                let mean_a = ya.iter().sum::<f64>() / len as f64;
                let mean_b = yb.iter().sum::<f64>() / len as f64;
                let num: f64 = ya
                    .iter()
                    .zip(yb)
                    .map(|(x, y)| (x - mean_a) * (y - mean_b))
                    .sum();
                let sa: f64 = ya.iter().map(|x| (x - mean_a).powi(2)).sum::<f64>().sqrt();
                let sb: f64 = yb.iter().map(|y| (y - mean_b).powi(2)).sum::<f64>().sqrt();
                if sa > 0.0 && sb > 0.0 {
                    num / (sa * sb)
                } else {
                    0.0
                }
            }
        }
    }
    /// Iterate over all entries.
    pub fn iter(&self) -> impl Iterator<Item = &DatabaseEntry> {
        self.entries.iter()
    }
}
/// Metadata associated with a [`SpectrumRecord`].
#[derive(Debug, Clone, Default)]
pub struct SpectrumMetadata {
    /// Title or name of the compound / measurement.
    pub title: String,
    /// Source technique that produced the spectrum.
    pub source: SpectrumSource,
    /// Date string of the measurement (ISO-8601 recommended).
    pub date: String,
    /// Instrument identifier / model.
    pub instrument: String,
    /// X-axis label (e.g. "Wavenumber (cm-1)").
    pub x_label: String,
    /// Y-axis label (e.g. "Absorbance").
    pub y_label: String,
    /// CAS Registry Number of the compound, if known.
    pub cas_number: String,
    /// Molecular formula string.
    pub molecular_formula: String,
    /// Molecular weight in g/mol, or NaN if unknown.
    pub molecular_weight: f64,
    /// Arbitrary additional key/value pairs from the file header.
    pub extra: HashMap<String, String>,
}
/// Spectral processing and analysis utilities.
pub struct SpectralAnalysis;
impl SpectralAnalysis {
    /// Find peaks in the Y array using a simple local-maximum criterion.
    ///
    /// A point is a peak if it is strictly greater than its `window` neighbours
    /// on each side and above the given `threshold`.
    pub fn find_peaks(record: &SpectrumRecord, window: usize, threshold: f64) -> Vec<PeakResult> {
        let n = record.y.len();
        let half = window.max(1);
        let mut peaks = Vec::new();
        for i in half..n.saturating_sub(half) {
            let y_i = record.y[i];
            if y_i < threshold {
                continue;
            }
            let is_max = (1..=half).all(|k| y_i > record.y[i - k] && y_i > record.y[i + k]);
            if is_max {
                let fwhm = Self::estimate_fwhm(record, i);
                peaks.push(PeakResult {
                    index: i,
                    x: record.x[i],
                    y: y_i,
                    fwhm,
                });
            }
        }
        peaks
    }
    /// Estimate FWHM for a peak at index `peak_idx` by walking outward until
    /// the intensity drops below half-maximum.
    fn estimate_fwhm(record: &SpectrumRecord, peak_idx: usize) -> f64 {
        let half_max = record.y[peak_idx] * 0.5;
        let n = record.y.len();
        let mut left_x = record.x[peak_idx];
        for i in (0..peak_idx).rev() {
            if record.y[i] <= half_max {
                if i + 1 < n {
                    let y0 = record.y[i];
                    let y1 = record.y[i + 1];
                    let x0 = record.x[i];
                    let x1 = record.x[i + 1];
                    if (y1 - y0).abs() > 1e-15 {
                        left_x = x0 + (half_max - y0) * (x1 - x0) / (y1 - y0);
                    }
                }
                break;
            }
        }
        let mut right_x = record.x[peak_idx];
        for i in peak_idx + 1..n {
            if record.y[i] <= half_max {
                if i > 0 {
                    let y0 = record.y[i - 1];
                    let y1 = record.y[i];
                    let x0 = record.x[i - 1];
                    let x1 = record.x[i];
                    if (y1 - y0).abs() > 1e-15 {
                        right_x = x0 + (half_max - y0) * (x1 - x0) / (y1 - y0);
                    }
                }
                break;
            }
        }
        (right_x - left_x).abs()
    }
    /// Asymmetric Least Squares (ALS) baseline correction.
    ///
    /// Returns the corrected Y array (input Y minus estimated baseline).
    ///
    /// # Parameters
    /// - `y` — input intensity array.
    /// - `lam` — smoothing parameter (typically 1e4 – 1e7).
    /// - `p` — asymmetry parameter (typically 0.001 – 0.1).
    /// - `iterations` — number of ALS iterations (typically 10–30).
    pub fn als_baseline(y: &[f64], lam: f64, p: f64, iterations: usize) -> Vec<f64> {
        let n = y.len();
        if n < 3 {
            return y.to_vec();
        }
        let mut weights = vec![1.0_f64; n];
        let mut baseline = y.to_vec();
        for _iter in 0..iterations {
            let mut diag = vec![0.0_f64; n];
            let mut sub = vec![0.0_f64; n - 1];
            let mut sup = vec![0.0_f64; n - 1];
            diag[..n].copy_from_slice(&weights[..n]);
            for i in 0..n - 1 {
                diag[i] += lam;
                diag[i + 1] += lam;
                sub[i] -= lam;
                sup[i] -= lam;
            }
            let rhs: Vec<f64> = (0..n).map(|i| weights[i] * y[i]).collect();
            baseline = Self::tridiagonal_solve(&diag, &sub, &sup, &rhs);
            for i in 0..n {
                if y[i] > baseline[i] {
                    weights[i] = p;
                } else {
                    weights[i] = 1.0 - p;
                }
            }
        }
        y.iter()
            .zip(baseline.iter())
            .map(|(&yi, &bi)| yi - bi)
            .collect()
    }
    /// Solve a tridiagonal system Ax = b using the Thomas algorithm.
    /// `sub` is the sub-diagonal (length n-1), `sup` is the super-diagonal.
    fn tridiagonal_solve(diag: &[f64], sub: &[f64], sup: &[f64], rhs: &[f64]) -> Vec<f64> {
        let n = diag.len();
        if n == 0 {
            return Vec::new();
        }
        let mut c_prime = vec![0.0; n - 1];
        let mut d_prime = vec![0.0; n];
        let mut x = vec![0.0; n];
        c_prime[0] = if diag[0].abs() > 1e-15 {
            sup[0] / diag[0]
        } else {
            0.0
        };
        d_prime[0] = if diag[0].abs() > 1e-15 {
            rhs[0] / diag[0]
        } else {
            0.0
        };
        for i in 1..n {
            let denom = diag[i]
                - if i > 0 {
                    sub[i - 1] * c_prime[i - 1]
                } else {
                    0.0
                };
            if denom.abs() < 1e-15 {
                d_prime[i] = 0.0;
                if i < n - 1 {
                    c_prime[i] = 0.0;
                }
            } else {
                d_prime[i] = (rhs[i] - sub[i - 1] * d_prime[i - 1]) / denom;
                if i < n - 1 {
                    c_prime[i] = sup[i] / denom;
                }
            }
        }
        x[n - 1] = d_prime[n - 1];
        for i in (0..n - 1).rev() {
            x[i] = d_prime[i] - c_prime[i] * x[i + 1];
        }
        x
    }
    /// Savitzky-Golay smoothing filter applied to Y values.
    ///
    /// Uses a 2nd-order polynomial over the given window.
    /// Returns the smoothed Y array (same length as input).
    pub fn savitzky_golay(y: &[f64], window: SgWindow) -> Vec<f64> {
        let n = y.len();
        let half = window.half();
        if n <= 2 * half {
            return y.to_vec();
        }
        let coeffs = Self::sg_coefficients(half);
        let mut out = y.to_vec();
        for i in half..n - half {
            let mut val = 0.0;
            for (k, &c) in coeffs.iter().enumerate() {
                val += c * y[i + k - half];
            }
            out[i] = val;
        }
        out
    }
    /// Compute Savitzky-Golay convolution coefficients for a 2nd-order polynomial
    /// and given half-window size.
    fn sg_coefficients(half: usize) -> Vec<f64> {
        let m = 2 * half + 1;
        let h = half as f64;
        let norm = (2.0 * h + 1.0) * (2.0 * h + 3.0) * (2.0 * h - 1.0) / 3.0;
        let mut coeffs = Vec::with_capacity(m);
        for j in -(half as i64)..=(half as i64) {
            let c = (3.0 * h * (h + 1.0) - 1.0 - 5.0 * j as f64 * j as f64) / norm;
            coeffs.push(c);
        }
        coeffs
    }
    /// Simple moving-average smoothing.
    ///
    /// The output at index `i` is the average of `y[i-half..=i+half]`.
    /// Boundary points are left unchanged.
    pub fn moving_average(y: &[f64], half_window: usize) -> Vec<f64> {
        let n = y.len();
        let mut out = y.to_vec();
        for i in half_window..n - half_window {
            let sum: f64 = y[i - half_window..=i + half_window].iter().sum();
            out[i] = sum / (2 * half_window + 1) as f64;
        }
        out
    }
    /// Compute the first derivative of Y with respect to X using central differences.
    ///
    /// Boundary points use one-sided differences.
    pub fn derivative(record: &SpectrumRecord) -> Vec<f64> {
        let n = record.y.len();
        if n < 2 {
            return vec![0.0; n];
        }
        let mut deriv = vec![0.0; n];
        for i in 1..n - 1 {
            let dx = record.x[i + 1] - record.x[i - 1];
            if dx.abs() > 1e-15 {
                deriv[i] = (record.y[i + 1] - record.y[i - 1]) / dx;
            }
        }
        if n >= 2 {
            let dx = record.x[1] - record.x[0];
            if dx.abs() > 1e-15 {
                deriv[0] = (record.y[1] - record.y[0]) / dx;
            }
            let dx_end = record.x[n - 1] - record.x[n - 2];
            if dx_end.abs() > 1e-15 {
                deriv[n - 1] = (record.y[n - 1] - record.y[n - 2]) / dx_end;
            }
        }
        deriv
    }
    /// Integrate the spectrum using the trapezoidal rule between `x_start` and `x_end`.
    pub fn integrate(record: &SpectrumRecord, x_start: f64, x_end: f64) -> f64 {
        let mut integral = 0.0;
        let n = record.x.len();
        for i in 1..n {
            let x0 = record.x[i - 1];
            let x1 = record.x[i];
            if x0 >= x_end || x1 <= x_start {
                continue;
            }
            let xa = x0.max(x_start);
            let xb = x1.min(x_end);
            let t0 = (xa - x0) / (x1 - x0);
            let t1 = (xb - x0) / (x1 - x0);
            let ya = record.y[i - 1] + t0 * (record.y[i] - record.y[i - 1]);
            let yb = record.y[i - 1] + t1 * (record.y[i] - record.y[i - 1]);
            integral += 0.5 * (ya + yb) * (xb - xa);
        }
        integral
    }
}
/// A spectral record holding X (wavenumber / frequency / m/z) and Y (intensity)
/// arrays together with rich metadata.
#[derive(Debug, Clone, Default)]
pub struct SpectrumRecord {
    /// X-axis values (wavenumber in cm⁻¹, frequency in Hz, m/z, …).
    pub x: Vec<f64>,
    /// Intensity / absorbance / transmittance values.
    pub y: Vec<f64>,
    /// Associated metadata.
    pub metadata: SpectrumMetadata,
}
impl SpectrumRecord {
    /// Create an empty `SpectrumRecord`.
    pub fn new() -> Self {
        Self::default()
    }
    /// Construct from raw arrays and a metadata object.
    pub fn from_arrays(x: Vec<f64>, y: Vec<f64>, metadata: SpectrumMetadata) -> Self {
        Self { x, y, metadata }
    }
    /// Return the number of data points.
    pub fn len(&self) -> usize {
        self.x.len()
    }
    /// Return `true` if the record contains no data points.
    pub fn is_empty(&self) -> bool {
        self.x.is_empty()
    }
    /// Return the minimum X value, or `f64::NAN` if empty.
    pub fn x_min(&self) -> f64 {
        self.x.iter().cloned().fold(f64::NAN, f64::min)
    }
    /// Return the maximum X value, or `f64::NAN` if empty.
    pub fn x_max(&self) -> f64 {
        self.x.iter().cloned().fold(f64::NAN, f64::max)
    }
    /// Return the maximum intensity value, or `f64::NAN` if empty.
    pub fn y_max(&self) -> f64 {
        self.y.iter().cloned().fold(f64::NAN, f64::max)
    }
    /// Return the minimum intensity value, or `f64::NAN` if empty.
    pub fn y_min(&self) -> f64 {
        self.y.iter().cloned().fold(f64::NAN, f64::min)
    }
    /// Normalize intensity so the maximum equals 1.0.
    /// No-op if the record is empty or has zero maximum.
    pub fn normalize(&mut self) {
        let max_y = self.y_max();
        if max_y > 0.0 && max_y.is_finite() {
            for v in &mut self.y {
                *v /= max_y;
            }
        }
    }
    /// Interpolate the intensity at an arbitrary X value using linear interpolation.
    /// Returns `None` if the spectrum is empty or `x_val` is out of range.
    pub fn interpolate_at(&self, x_val: f64) -> Option<f64> {
        if self.x.is_empty() {
            return None;
        }
        let idx = self.x.partition_point(|&v| v < x_val);
        if idx == 0 {
            if (self.x[0] - x_val).abs() < f64::EPSILON {
                return Some(self.y[0]);
            }
            return None;
        }
        if idx >= self.x.len() {
            let last = self.x.len() - 1;
            if (self.x[last] - x_val).abs() < f64::EPSILON {
                return Some(self.y[last]);
            }
            return None;
        }
        let x0 = self.x[idx - 1];
        let x1 = self.x[idx];
        let y0 = self.y[idx - 1];
        let y1 = self.y[idx];
        let t = (x_val - x0) / (x1 - x0);
        Some(y0 + t * (y1 - y0))
    }
}
/// Reader for the ANIF binary spectral format.
///
/// The format layout (little-endian or big-endian based on header flag):
/// ```text
/// [4] Magic "ANIF"
/// [1] Major version
/// [1] Minor version
/// [1] Endian flag  (0 = LE, 1 = BE)
/// [1] Reserved
/// [4] u32  num_points
/// [8] f64  x_start
/// [8] f64  x_delta
/// [64] title (UTF-8, null-padded)
/// [32] instrument (UTF-8, null-padded)
/// [16] date (UTF-8, null-padded)
/// [num_points * 8] f64 Y values
/// ```
#[derive(Debug, Default)]
pub struct AnifReader {
    /// Parsed header.
    pub header: AnifHeader,
    /// Decoded spectrum record.
    pub record: SpectrumRecord,
}
impl AnifReader {
    /// Create a new, empty reader.
    pub fn new() -> Self {
        Self::default()
    }
    /// Parse an ANIF byte stream.
    pub fn parse<R: Read>(&mut self, mut reader: R) -> crate::Result<()> {
        let mut buf = Vec::new();
        reader.read_to_end(&mut buf).map_err(crate::Error::Io)?;
        self.parse_bytes(&buf)
    }
    /// Parse from a byte slice.
    pub fn parse_bytes(&mut self, buf: &[u8]) -> crate::Result<()> {
        if buf.len() < 4 || &buf[0..4] != ANIF_MAGIC {
            return Err(crate::Error::Parse("Not an ANIF file (bad magic)".into()));
        }
        if buf.len() < 132 {
            return Err(crate::Error::Parse("ANIF header too short".into()));
        }
        let major = buf[4];
        let minor = buf[5];
        let endian_flag = buf[6];
        let endian = if endian_flag == 0 {
            AnifEndian::Little
        } else {
            AnifEndian::Big
        };
        let read_u32 = |b: &[u8], off: usize| -> u32 {
            let arr: [u8; 4] = b[off..off + 4].try_into().unwrap_or([0u8; 4]);
            if endian == AnifEndian::Little {
                u32::from_le_bytes(arr)
            } else {
                u32::from_be_bytes(arr)
            }
        };
        let read_f64 = |b: &[u8], off: usize| -> f64 {
            let arr: [u8; 8] = b[off..off + 8].try_into().unwrap_or([0u8; 8]);
            if endian == AnifEndian::Little {
                f64::from_le_bytes(arr)
            } else {
                f64::from_be_bytes(arr)
            }
        };
        let read_str = |b: &[u8], off: usize, len: usize| -> String {
            let slice = &b[off..off + len];
            let end = slice.iter().position(|&c| c == 0).unwrap_or(len);
            String::from_utf8_lossy(&slice[..end]).into_owned()
        };
        let num_points = read_u32(buf, 8);
        let x_start = read_f64(buf, 12);
        let x_delta = read_f64(buf, 20);
        let title = read_str(buf, 28, 64);
        let instrument = read_str(buf, 92, 32);
        let date = read_str(buf, 124, 16);
        self.header = AnifHeader {
            version: (major, minor),
            num_points,
            x_start,
            x_delta,
            endian,
            title: title.clone(),
            instrument: instrument.clone(),
            date: date.clone(),
        };
        let data_offset = 140usize;
        let expected_len = data_offset + num_points as usize * 8;
        if buf.len() < expected_len {
            return Err(crate::Error::Parse(format!(
                "ANIF data block too short: expected {} bytes, got {}",
                expected_len,
                buf.len()
            )));
        }
        let mut x_vals = Vec::with_capacity(num_points as usize);
        let mut y_vals = Vec::with_capacity(num_points as usize);
        for i in 0..num_points as usize {
            let off = data_offset + i * 8;
            x_vals.push(x_start + x_delta * i as f64);
            y_vals.push(read_f64(buf, off));
        }
        let meta = SpectrumMetadata {
            title,
            instrument,
            date,
            ..Default::default()
        };
        self.record = SpectrumRecord::from_arrays(x_vals, y_vals, meta);
        Ok(())
    }
    /// Encode a `SpectrumRecord` to ANIF bytes (little-endian).
    pub fn encode(record: &SpectrumRecord) -> Vec<u8> {
        let n = record.len();
        let mut buf = Vec::with_capacity(140 + n * 8);
        buf.extend_from_slice(ANIF_MAGIC);
        buf.push(1);
        buf.push(0);
        buf.push(0);
        buf.push(0);
        buf.extend_from_slice(&(n as u32).to_le_bytes());
        let x_start = if record.x.is_empty() {
            0.0
        } else {
            record.x[0]
        };
        let x_delta = if record.x.len() >= 2 {
            record.x[1] - record.x[0]
        } else {
            1.0
        };
        buf.extend_from_slice(&x_start.to_le_bytes());
        buf.extend_from_slice(&x_delta.to_le_bytes());
        let title_bytes = record.metadata.title.as_bytes();
        let mut tmp = [0u8; 64];
        let copy_len = title_bytes.len().min(63);
        tmp[..copy_len].copy_from_slice(&title_bytes[..copy_len]);
        buf.extend_from_slice(&tmp);
        let instr_bytes = record.metadata.instrument.as_bytes();
        let mut tmp2 = [0u8; 32];
        let copy_len2 = instr_bytes.len().min(31);
        tmp2[..copy_len2].copy_from_slice(&instr_bytes[..copy_len2]);
        buf.extend_from_slice(&tmp2);
        let date_bytes = record.metadata.date.as_bytes();
        let mut tmp3 = [0u8; 16];
        let copy_len3 = date_bytes.len().min(15);
        tmp3[..copy_len3].copy_from_slice(&date_bytes[..copy_len3]);
        buf.extend_from_slice(&tmp3);
        for &yv in &record.y {
            buf.extend_from_slice(&yv.to_le_bytes());
        }
        buf
    }
}
/// Similarity metric for spectral search.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum SimilarityMetric {
    /// Normalized dot product between intensity vectors.
    DotProduct,
    /// Spectral Angle Mapper (SAM) — angle in radians between vectors.
    SpectralAngleMapper,
    /// Euclidean distance (smaller = more similar).
    Euclidean,
    /// Pearson correlation coefficient.
    Pearson,
}
/// Result of a similarity search.
#[derive(Debug, Clone)]
pub struct SearchResult {
    /// Database entry ID.
    pub id: usize,
    /// Title of the matching spectrum.
    pub title: String,
    /// Similarity score (interpretation depends on metric).
    pub score: f64,
}
/// Writer for various spectral file formats.
pub struct SpectralExport;
impl SpectralExport {
    /// Write a spectrum to a CSV file (`x,y` per line).
    pub fn write_csv<W: Write>(record: &SpectrumRecord, writer: &mut W) -> crate::Result<()> {
        writeln!(writer, "# {}", record.metadata.title).map_err(crate::Error::Io)?;
        writeln!(writer, "x,y").map_err(crate::Error::Io)?;
        for (&x, &y) in record.x.iter().zip(record.y.iter()) {
            writeln!(writer, "{},{}", x, y).map_err(crate::Error::Io)?;
        }
        Ok(())
    }
    /// Write a spectrum in JCAMP-DX format.
    pub fn write_jcamp_dx<W: Write>(record: &SpectrumRecord, writer: &mut W) -> crate::Result<()> {
        let title = &record.metadata.title;
        let n = record.len();
        let x_start = record.x_min();
        let x_end = record.x_max();
        let y_min = record.y_min();
        let y_max = record.y_max();
        writeln!(writer, "##TITLE={}", title).map_err(crate::Error::Io)?;
        writeln!(writer, "##JCAMP-DX=4.24").map_err(crate::Error::Io)?;
        writeln!(writer, "##DATA TYPE={}", record.metadata.source).map_err(crate::Error::Io)?;
        writeln!(writer, "##DATE={}", record.metadata.date).map_err(crate::Error::Io)?;
        writeln!(writer, "##INSTRUMENT={}", record.metadata.instrument)
            .map_err(crate::Error::Io)?;
        writeln!(writer, "##XUNITS={}", record.metadata.x_label).map_err(crate::Error::Io)?;
        writeln!(writer, "##YUNITS={}", record.metadata.y_label).map_err(crate::Error::Io)?;
        writeln!(writer, "##NPOINTS={}", n).map_err(crate::Error::Io)?;
        writeln!(writer, "##FIRSTX={}", x_start).map_err(crate::Error::Io)?;
        writeln!(writer, "##LASTX={}", x_end).map_err(crate::Error::Io)?;
        writeln!(writer, "##MINY={}", y_min).map_err(crate::Error::Io)?;
        writeln!(writer, "##MAXY={}", y_max).map_err(crate::Error::Io)?;
        writeln!(writer, "##XFACTOR=1.0").map_err(crate::Error::Io)?;
        writeln!(writer, "##YFACTOR=1.0").map_err(crate::Error::Io)?;
        writeln!(writer, "##XYDATA=(X++(Y..Y))").map_err(crate::Error::Io)?;
        let chunk_size = 8usize;
        let mut i = 0;
        while i < n {
            let x_prefix = record.x[i];
            write!(writer, "{:.6}", x_prefix).map_err(crate::Error::Io)?;
            let end = (i + chunk_size).min(n);
            for j in i..end {
                write!(writer, " {:.6}", record.y[j]).map_err(crate::Error::Io)?;
            }
            writeln!(writer).map_err(crate::Error::Io)?;
            i = end;
        }
        writeln!(writer, "##END=").map_err(crate::Error::Io)?;
        Ok(())
    }
    /// Write a spectrum as plain-text (space-separated x y columns).
    pub fn write_plain_text<W: Write>(
        record: &SpectrumRecord,
        writer: &mut W,
    ) -> crate::Result<()> {
        writeln!(writer, "# Spectrum: {}", record.metadata.title).map_err(crate::Error::Io)?;
        writeln!(writer, "# Instrument: {}", record.metadata.instrument)
            .map_err(crate::Error::Io)?;
        writeln!(writer, "# Date: {}", record.metadata.date).map_err(crate::Error::Io)?;
        writeln!(writer, "# X_label: {}", record.metadata.x_label).map_err(crate::Error::Io)?;
        writeln!(writer, "# Y_label: {}", record.metadata.y_label).map_err(crate::Error::Io)?;
        for (&x, &y) in record.x.iter().zip(record.y.iter()) {
            writeln!(writer, "{:.8e}  {:.8e}", x, y).map_err(crate::Error::Io)?;
        }
        Ok(())
    }
    /// Write multiple spectra to a single CSV with columns `x,y1,y2,...`.
    ///
    /// All records must share the same X grid. If they differ, the X from the
    /// first record is used.
    pub fn write_multi_csv<W: Write>(
        records: &[SpectrumRecord],
        writer: &mut W,
    ) -> crate::Result<()> {
        if records.is_empty() {
            return Ok(());
        }
        write!(writer, "x").map_err(crate::Error::Io)?;
        for r in records {
            write!(writer, ",{}", r.metadata.title).map_err(crate::Error::Io)?;
        }
        writeln!(writer).map_err(crate::Error::Io)?;
        let n = records[0].x.len();
        for i in 0..n {
            write!(writer, "{}", records[0].x[i]).map_err(crate::Error::Io)?;
            for r in records {
                if i < r.y.len() {
                    write!(writer, ",{}", r.y[i]).map_err(crate::Error::Io)?;
                } else {
                    write!(writer, ",").map_err(crate::Error::Io)?;
                }
            }
            writeln!(writer).map_err(crate::Error::Io)?;
        }
        Ok(())
    }
}
/// Entry stored in the spectral database.
#[derive(Debug, Clone)]
pub struct DatabaseEntry {
    /// Unique numeric key.
    pub id: usize,
    /// The stored spectrum.
    pub record: SpectrumRecord,
}
/// Parser for JCAMP-DX spectral files (`.jdx` / `.dx`).
///
/// Supports the following JCAMP-DX features:
/// - Single-block files with `##XYDATA=(X++(Y..Y))` data records.
/// - Multi-block compound files (multiple `##TITLE` blocks).
/// - `##NTUPLES` sections for multi-channel data.
/// - `##DELTAX` and `##FIRSTX` / `##LASTX` for compressed X grids.
/// - `##XFACTOR` / `##YFACTOR` scaling factors.
#[derive(Debug, Default)]
pub struct JcampDxReader {
    /// All spectrum blocks parsed from the file.
    pub records: Vec<SpectrumRecord>,
}
impl JcampDxReader {
    /// Create a new, empty reader.
    pub fn new() -> Self {
        Self::default()
    }
    /// Parse JCAMP-DX content from a `Read` source.
    pub fn parse<R: Read>(&mut self, reader: R) -> crate::Result<()> {
        let buf = BufReader::new(reader);
        let mut current = SpectrumRecord::new();
        let mut section = JdxSection::Header;
        let mut x_factor = 1.0_f64;
        let mut y_factor = 1.0_f64;
        let mut firstx = f64::NAN;
        let mut deltax = f64::NAN;
        let mut npoints: usize = 0;
        let mut raw_y: Vec<f64> = Vec::new();
        let mut in_block = false;
        let mut ntuple_x: Vec<f64> = Vec::new();
        let mut ntuple_y: Vec<f64> = Vec::new();
        let mut ntuple_var_names: Vec<String> = Vec::new();
        let mut ntuple_collecting = false;
        for line_res in buf.lines() {
            let line = line_res.map_err(crate::Error::Io)?;
            let trimmed = line.trim();
            if trimmed.is_empty() || trimmed.starts_with("$$") {
                continue;
            }
            if let Some(content) = trimmed.strip_prefix("##") {
                let eq_pos = content.find('=');
                let (label, value) = if let Some(pos) = eq_pos {
                    (
                        content[..pos].trim().to_uppercase(),
                        content[pos + 1..].trim().to_string(),
                    )
                } else {
                    (content.trim().to_uppercase(), String::new())
                };
                match label.as_str() {
                    "TITLE" => {
                        if in_block && !current.is_empty() {
                            self.records.push(current.clone());
                            current = SpectrumRecord::new();
                            x_factor = 1.0;
                            y_factor = 1.0;
                            firstx = f64::NAN;
                            deltax = f64::NAN;
                            npoints = 0;
                            raw_y.clear();
                        }
                        in_block = true;
                        current.metadata.title = value;
                        section = JdxSection::Header;
                    }
                    "END" => {
                        if !raw_y.is_empty() {
                            Self::build_xy_from_raw(
                                &mut current,
                                &raw_y,
                                firstx,
                                deltax,
                                x_factor,
                                y_factor,
                            );
                            raw_y.clear();
                        }
                        if !ntuple_x.is_empty() {
                            current.x = ntuple_x.clone();
                            current.y = ntuple_y.clone();
                            ntuple_x.clear();
                            ntuple_y.clear();
                        }
                        if !current.is_empty() {
                            self.records.push(current.clone());
                        }
                        current = SpectrumRecord::new();
                        in_block = false;
                        x_factor = 1.0;
                        y_factor = 1.0;
                        firstx = f64::NAN;
                        deltax = f64::NAN;
                        npoints = 0;
                        section = JdxSection::Header;
                    }
                    "DATA TYPE" | "DATATYPE" => {
                        let v = value.to_uppercase();
                        if v.contains("INFRARED") || v.contains("IR") {
                            current.metadata.source = SpectrumSource::Ftir;
                        } else if v.contains("RAMAN") {
                            current.metadata.source = SpectrumSource::Raman;
                        } else if v.contains("UV") || v.contains("VIS") {
                            current.metadata.source = SpectrumSource::UvVis;
                        } else if v.contains("NMR") {
                            current.metadata.source = SpectrumSource::Nmr;
                        } else if v.contains("MASS") {
                            current.metadata.source = SpectrumSource::MassSpec;
                        }
                    }
                    "DATE" => current.metadata.date = value,
                    "INSTRUMENT" | "SPECTROMETER/DATA SYSTEM" => {
                        current.metadata.instrument = value;
                    }
                    "XUNITS" => current.metadata.x_label = value,
                    "YUNITS" => current.metadata.y_label = value,
                    "CAS REGISTRY NO" | "CASREGNO" => current.metadata.cas_number = value,
                    "MOLFORM" | "MOLECULAR FORMULA" => {
                        current.metadata.molecular_formula = value;
                    }
                    "MW" | "MOLECULAR WEIGHT" => {
                        current.metadata.molecular_weight =
                            value.parse::<f64>().unwrap_or(f64::NAN);
                    }
                    "XFACTOR" => x_factor = value.parse::<f64>().unwrap_or(1.0),
                    "YFACTOR" => y_factor = value.parse::<f64>().unwrap_or(1.0),
                    "FIRSTX" => firstx = value.parse::<f64>().unwrap_or(f64::NAN),
                    "DELTAX" => deltax = value.parse::<f64>().unwrap_or(f64::NAN),
                    "NPOINTS" => npoints = value.parse::<usize>().unwrap_or(0),
                    "XYDATA" => {
                        section = JdxSection::XyData;
                        raw_y.clear();
                    }
                    "XYPOINTS" => {
                        section = JdxSection::XyData;
                        raw_y.clear();
                    }
                    "NTUPLES" => {
                        section = JdxSection::NtupleDecl;
                        ntuple_var_names.clear();
                        ntuple_x.clear();
                        ntuple_y.clear();
                        ntuple_collecting = false;
                    }
                    "VAR_NAME" => {
                        ntuple_var_names = value.split(',').map(|s| s.trim().to_string()).collect();
                    }
                    "PAGE" => {
                        ntuple_collecting = true;
                        section = JdxSection::NtupleData;
                    }
                    "END NTUPLES" => {
                        ntuple_collecting = false;
                        section = JdxSection::Header;
                    }
                    _ => {
                        current.metadata.extra.insert(label.to_string(), value);
                    }
                }
                continue;
            }
            match section {
                JdxSection::XyData => {
                    Self::parse_xydata_line(trimmed, &mut current, &mut raw_y, x_factor, y_factor);
                }
                JdxSection::NtupleData if ntuple_collecting => {
                    Self::parse_ntuple_line(trimmed, &mut ntuple_x, &mut ntuple_y);
                }
                _ => {}
            }
            let _ = npoints;
        }
        if in_block && !current.is_empty() {
            if !raw_y.is_empty() {
                Self::build_xy_from_raw(&mut current, &raw_y, firstx, deltax, x_factor, y_factor);
            }
            if !ntuple_x.is_empty() {
                current.x = ntuple_x;
                current.y = ntuple_y;
            }
            if !current.is_empty() {
                self.records.push(current);
            }
        }
        Ok(())
    }
    /// Parse a line in `XYDATA=(X++(Y..Y))` mode.
    fn parse_xydata_line(
        line: &str,
        record: &mut SpectrumRecord,
        raw_y: &mut Vec<f64>,
        x_factor: f64,
        y_factor: f64,
    ) {
        let tokens: Vec<&str> = line.split_whitespace().collect();
        if tokens.is_empty() {
            return;
        }
        if let Ok(x_val) = tokens[0].parse::<f64>() {
            record.x.push(x_val * x_factor);
            for tok in &tokens[1..] {
                if let Ok(y_val) = tok.parse::<f64>() {
                    raw_y.push(y_val * y_factor);
                }
            }
        }
    }
    /// Build the full X/Y arrays from compressed X++(Y..Y) raw data.
    fn build_xy_from_raw(
        record: &mut SpectrumRecord,
        raw_y: &[f64],
        firstx: f64,
        deltax: f64,
        x_factor: f64,
        y_factor: f64,
    ) {
        if !firstx.is_nan() && !deltax.is_nan() {
            record.x.clear();
            record.y.clear();
            for (i, &yv) in raw_y.iter().enumerate() {
                record.x.push((firstx + deltax * i as f64) * x_factor);
                record.y.push(yv * y_factor);
            }
        } else {
            record.y.clear();
            record.y.extend_from_slice(raw_y);
            let _ = (x_factor, y_factor);
        }
    }
    /// Parse NTUPLE data line (tab or space separated x,y pairs).
    fn parse_ntuple_line(line: &str, xs: &mut Vec<f64>, ys: &mut Vec<f64>) {
        let tokens: Vec<&str> = line.split_whitespace().collect();
        if tokens.len() >= 2
            && let (Ok(x), Ok(y)) = (tokens[0].parse::<f64>(), tokens[1].parse::<f64>())
        {
            xs.push(x);
            ys.push(y);
        }
    }
    /// Return the number of parsed blocks.
    pub fn count(&self) -> usize {
        self.records.len()
    }
}
/// Savitzky-Golay filter window size options.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum SgWindow {
    /// 5-point window.
    W5,
    /// 7-point window.
    W7,
    /// 11-point window.
    W11,
    /// 15-point window.
    W15,
    /// 25-point window.
    W25,
}
impl SgWindow {
    /// Half-window size (number of points on each side).
    pub fn half(self) -> usize {
        match self {
            SgWindow::W5 => 2,
            SgWindow::W7 => 3,
            SgWindow::W11 => 5,
            SgWindow::W15 => 7,
            SgWindow::W25 => 12,
        }
    }
}
/// Header parsed from an ANIF file.
#[derive(Debug, Clone)]
pub struct AnifHeader {
    /// Format version (major.minor).
    pub version: (u8, u8),
    /// Number of data points.
    pub num_points: u32,
    /// Starting X value.
    pub x_start: f64,
    /// X increment per point.
    pub x_delta: f64,
    /// Byte endianness.
    pub endian: AnifEndian,
    /// User-defined title string (up to 64 bytes).
    pub title: String,
    /// Instrument identifier (up to 32 bytes).
    pub instrument: String,
    /// Date string (up to 16 bytes).
    pub date: String,
}
/// Source type that generated the spectrum.
#[derive(Debug, Clone, PartialEq)]
pub enum SpectrumSource {
    /// Fourier-Transform Infrared spectroscopy.
    Ftir,
    /// Near-Infrared spectroscopy.
    Nir,
    /// Raman spectroscopy.
    Raman,
    /// UV-Visible spectroscopy.
    UvVis,
    /// Mass spectrometry.
    MassSpec,
    /// Nuclear Magnetic Resonance spectroscopy.
    Nmr,
    /// Unknown / unspecified source.
    Unknown(String),
}
/// State machine for JCAMP-DX section parsing.
#[derive(Debug, PartialEq)]
pub(super) enum JdxSection {
    Header,
    XyData,
    NtupleDecl,
    NtupleData,
}
