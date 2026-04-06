// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Experimental data I/O for the OxiPhysics engine.
//!
//! Provides readers and processors for common experimental data formats:
//! - Time series with metadata (sensor, units, sampling rate), outlier detection
//! - Digital Image Correlation (DIC) displacement/strain fields
//! - Particle Image Velocimetry (PIV) velocity fields
//! - Strain gauge (rosette) data and principal strains
//! - Dynamic Mechanical Analysis (DMA) frequency sweeps and master curves
//! - Simulation vs experiment comparison metrics

#![allow(dead_code)]
#![allow(clippy::too_many_arguments)]

use std::collections::HashMap;
use std::fmt;

// ---------------------------------------------------------------------------
// ExperimentalDataset — time series with metadata and outlier flags
// ---------------------------------------------------------------------------

/// Metadata associated with an experimental time series.
#[derive(Debug, Clone)]
pub struct SensorMetadata {
    /// Sensor identifier string.
    pub sensor_id: String,
    /// Physical quantity measured (e.g., "displacement", "force").
    pub quantity: String,
    /// Physical unit of the measured values.
    pub unit: String,
    /// Nominal sampling rate in Hz.
    pub sampling_rate_hz: f64,
    /// Optional location description.
    pub location: Option<String>,
    /// Arbitrary key-value pairs for additional metadata.
    pub extra: HashMap<String, String>,
}

impl SensorMetadata {
    /// Create new sensor metadata.
    pub fn new(
        sensor_id: impl Into<String>,
        quantity: impl Into<String>,
        unit: impl Into<String>,
        sampling_rate_hz: f64,
    ) -> Self {
        Self {
            sensor_id: sensor_id.into(),
            quantity: quantity.into(),
            unit: unit.into(),
            sampling_rate_hz,
            location: None,
            extra: HashMap::new(),
        }
    }

    /// Set the location description.
    pub fn with_location(mut self, location: impl Into<String>) -> Self {
        self.location = Some(location.into());
        self
    }

    /// Insert an extra key-value pair.
    pub fn with_extra(mut self, key: impl Into<String>, value: impl Into<String>) -> Self {
        self.extra.insert(key.into(), value.into());
        self
    }
}

/// A single time-series sample with an outlier flag.
#[derive(Debug, Clone, PartialEq)]
pub struct DataSample {
    /// Time stamp in seconds.
    pub time: f64,
    /// Measured value.
    pub value: f64,
    /// Whether this sample has been flagged as an outlier.
    pub is_outlier: bool,
}

impl DataSample {
    /// Create a new sample (not an outlier).
    pub fn new(time: f64, value: f64) -> Self {
        Self {
            time,
            value,
            is_outlier: false,
        }
    }
}

/// Experimental time series dataset.
#[derive(Debug, Clone)]
pub struct ExperimentalDataset {
    /// Sensor metadata.
    pub metadata: SensorMetadata,
    /// Ordered samples.
    pub samples: Vec<DataSample>,
}

impl ExperimentalDataset {
    /// Create an empty dataset with the given metadata.
    pub fn new(metadata: SensorMetadata) -> Self {
        Self {
            metadata,
            samples: Vec::new(),
        }
    }

    /// Append a sample.
    pub fn push(&mut self, time: f64, value: f64) {
        self.samples.push(DataSample::new(time, value));
    }

    /// Return the number of samples.
    pub fn len(&self) -> usize {
        self.samples.len()
    }

    /// Return `true` if there are no samples.
    pub fn is_empty(&self) -> bool {
        self.samples.is_empty()
    }

    /// Compute the mean of non-outlier values.
    pub fn mean(&self) -> f64 {
        let vals: Vec<f64> = self
            .samples
            .iter()
            .filter(|s| !s.is_outlier)
            .map(|s| s.value)
            .collect();
        if vals.is_empty() {
            return 0.0;
        }
        vals.iter().sum::<f64>() / vals.len() as f64
    }

    /// Compute the standard deviation of non-outlier values.
    pub fn std_dev(&self) -> f64 {
        let vals: Vec<f64> = self
            .samples
            .iter()
            .filter(|s| !s.is_outlier)
            .map(|s| s.value)
            .collect();
        if vals.len() < 2 {
            return 0.0;
        }
        let m = vals.iter().sum::<f64>() / vals.len() as f64;
        let var = vals.iter().map(|v| (v - m).powi(2)).sum::<f64>() / (vals.len() - 1) as f64;
        var.sqrt()
    }

    /// Flag outliers using the IQR method (fence factor `k`, typically 1.5).
    pub fn flag_outliers_iqr(&mut self, k: f64) {
        let mut vals: Vec<f64> = self.samples.iter().map(|s| s.value).collect();
        vals.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        let n = vals.len();
        if n < 4 {
            return;
        }
        let q1 = percentile(&vals, 25.0);
        let q3 = percentile(&vals, 75.0);
        let iqr = q3 - q1;
        let lo = q1 - k * iqr;
        let hi = q3 + k * iqr;
        for s in &mut self.samples {
            s.is_outlier = s.value < lo || s.value > hi;
        }
    }

    /// Flag outliers using the z-score method (threshold `z_thresh`, typically 3.0).
    pub fn flag_outliers_zscore(&mut self, z_thresh: f64) {
        let m = self.mean();
        let sd = self.std_dev();
        if sd == 0.0 {
            return;
        }
        for s in &mut self.samples {
            s.is_outlier = ((s.value - m) / sd).abs() > z_thresh;
        }
    }

    /// Return only the outlier samples.
    pub fn outliers(&self) -> Vec<&DataSample> {
        self.samples.iter().filter(|s| s.is_outlier).collect()
    }

    /// Parse from CSV text: `time,value` rows.
    pub fn from_csv(csv: &str, metadata: SensorMetadata) -> Result<Self, String> {
        let mut ds = Self::new(metadata);
        for (i, line) in csv.lines().enumerate() {
            let line = line.trim();
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            let parts: Vec<&str> = line.split(',').collect();
            if parts.len() < 2 {
                return Err(format!("line {}: expected 2 columns", i + 1));
            }
            let t = parts[0]
                .trim()
                .parse::<f64>()
                .map_err(|e| format!("line {}: time parse error: {}", i + 1, e))?;
            let v = parts[1]
                .trim()
                .parse::<f64>()
                .map_err(|e| format!("line {}: value parse error: {}", i + 1, e))?;
            ds.push(t, v);
        }
        Ok(ds)
    }

    /// Serialize to CSV text.
    pub fn to_csv(&self) -> String {
        let mut out = String::from("# time,value,is_outlier\n");
        for s in &self.samples {
            out.push_str(&format!("{},{},{}\n", s.time, s.value, s.is_outlier as u8));
        }
        out
    }

    /// Return min and max of non-outlier values, or `None` if empty.
    pub fn range(&self) -> Option<(f64, f64)> {
        let vals: Vec<f64> = self
            .samples
            .iter()
            .filter(|s| !s.is_outlier)
            .map(|s| s.value)
            .collect();
        if vals.is_empty() {
            return None;
        }
        let min = vals.iter().cloned().fold(f64::INFINITY, f64::min);
        let max = vals.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        Some((min, max))
    }

    /// Downsample by keeping every `n`-th sample.
    pub fn downsample(&self, n: usize) -> Self {
        let samples = self
            .samples
            .iter()
            .enumerate()
            .filter(|(i, _)| i % n.max(1) == 0)
            .map(|(_, s)| s.clone())
            .collect();
        Self {
            metadata: self.metadata.clone(),
            samples,
        }
    }
}

/// Compute percentile (0–100) of a sorted slice.
fn percentile(sorted: &[f64], p: f64) -> f64 {
    if sorted.is_empty() {
        return 0.0;
    }
    let idx = (p / 100.0 * (sorted.len() - 1) as f64).clamp(0.0, (sorted.len() - 1) as f64);
    let lo = idx.floor() as usize;
    let hi = idx.ceil() as usize;
    let frac = idx - lo as f64;
    sorted[lo] * (1.0 - frac) + sorted[hi] * frac
}

// ---------------------------------------------------------------------------
// DigitalImageCorrelation — DIC displacement/strain field
// ---------------------------------------------------------------------------

/// A 2-D point in a DIC field.
#[derive(Debug, Clone, PartialEq)]
pub struct DicPoint {
    /// Reference x coordinate (pixels or mm).
    pub x: f64,
    /// Reference y coordinate.
    pub y: f64,
    /// Displacement in x.
    pub u: f64,
    /// Displacement in y.
    pub v: f64,
    /// Correlation coefficient (quality metric, 0–1).
    pub correlation: f64,
}

impl DicPoint {
    /// Create a new DIC point.
    pub fn new(x: f64, y: f64, u: f64, v: f64, correlation: f64) -> Self {
        Self {
            x,
            y,
            u,
            v,
            correlation,
        }
    }

    /// Displacement magnitude.
    pub fn displacement_magnitude(&self) -> f64 {
        (self.u * self.u + self.v * self.v).sqrt()
    }
}

/// Computed strain components at a point.
#[derive(Debug, Clone, PartialEq)]
pub struct StrainTensor2D {
    /// Normal strain in x direction (εxx).
    pub exx: f64,
    /// Normal strain in y direction (εyy).
    pub eyy: f64,
    /// Shear strain (εxy = γxy / 2).
    pub exy: f64,
}

impl StrainTensor2D {
    /// Principal strains (ε1, ε2) and principal angle θ in radians.
    pub fn principal_strains(&self) -> (f64, f64, f64) {
        let avg = (self.exx + self.eyy) / 2.0;
        let r = ((self.exx - self.eyy).powi(2) / 4.0 + self.exy.powi(2)).sqrt();
        let e1 = avg + r;
        let e2 = avg - r;
        let theta = 0.5 * (2.0 * self.exy).atan2(self.exx - self.eyy);
        (e1, e2, theta)
    }

    /// Von Mises equivalent strain.
    pub fn von_mises(&self) -> f64 {
        let (e1, e2, _) = self.principal_strains();
        ((e1 - e2).powi(2) + e2.powi(2) + e1.powi(2)).sqrt() / std::f64::consts::SQRT_2
    }
}

/// Digital Image Correlation displacement field reader and processor.
#[derive(Debug, Clone)]
pub struct DigitalImageCorrelation {
    /// DIC data points.
    pub points: Vec<DicPoint>,
    /// Minimum accepted correlation coefficient.
    pub min_correlation: f64,
}

impl DigitalImageCorrelation {
    /// Create an empty DIC reader with a given minimum correlation threshold.
    pub fn new(min_correlation: f64) -> Self {
        Self {
            points: Vec::new(),
            min_correlation,
        }
    }

    /// Parse from CSV: `x,y,u,v,correlation` rows.
    pub fn from_csv(csv: &str, min_correlation: f64) -> Result<Self, String> {
        let mut dic = Self::new(min_correlation);
        for (i, line) in csv.lines().enumerate() {
            let line = line.trim();
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            let p: Vec<&str> = line.split(',').collect();
            if p.len() < 5 {
                return Err(format!("line {}: expected 5 columns", i + 1));
            }
            let parse = |s: &str| s.trim().parse::<f64>().map_err(|e| format!("parse: {}", e));
            let x = parse(p[0])?;
            let y = parse(p[1])?;
            let u = parse(p[2])?;
            let v = parse(p[3])?;
            let c = parse(p[4])?;
            dic.points.push(DicPoint::new(x, y, u, v, c));
        }
        Ok(dic)
    }

    /// Filter out points below the minimum correlation threshold.
    pub fn filter_by_correlation(&mut self) {
        self.points
            .retain(|p| p.correlation >= self.min_correlation);
    }

    /// Compute strain field using finite differences on a regular grid.
    ///
    /// Returns a vector of `(x, y, StrainTensor2D)` for interior points.
    pub fn compute_strain_field(&self, grid_spacing: f64) -> Vec<(f64, f64, StrainTensor2D)> {
        if self.points.len() < 4 || grid_spacing <= 0.0 {
            return Vec::new();
        }
        // Build lookup by grid indices (approximate)
        let mut map: HashMap<(i64, i64), &DicPoint> = HashMap::new();
        for p in &self.points {
            let ix = (p.x / grid_spacing).round() as i64;
            let iy = (p.y / grid_spacing).round() as i64;
            map.insert((ix, iy), p);
        }
        let mut result = Vec::new();
        for (&(ix, iy), &pt) in &map {
            let r = map.get(&(ix + 1, iy));
            let l = map.get(&(ix - 1, iy));
            let u = map.get(&(ix, iy + 1));
            let d = map.get(&(ix, iy - 1));
            let (du_dx, dv_dx) = match (r, l) {
                (Some(pr), Some(pl)) => (
                    (pr.u - pl.u) / (2.0 * grid_spacing),
                    (pr.v - pl.v) / (2.0 * grid_spacing),
                ),
                (Some(pr), None) => ((pr.u - pt.u) / grid_spacing, (pr.v - pt.v) / grid_spacing),
                (None, Some(pl)) => ((pt.u - pl.u) / grid_spacing, (pt.v - pl.v) / grid_spacing),
                _ => continue,
            };
            let (du_dy, dv_dy) = match (u, d) {
                (Some(pu), Some(pd)) => (
                    (pu.u - pd.u) / (2.0 * grid_spacing),
                    (pu.v - pd.v) / (2.0 * grid_spacing),
                ),
                (Some(pu), None) => ((pu.u - pt.u) / grid_spacing, (pu.v - pt.v) / grid_spacing),
                (None, Some(pd)) => ((pt.u - pd.u) / grid_spacing, (pt.v - pd.v) / grid_spacing),
                _ => continue,
            };
            let strain = StrainTensor2D {
                exx: du_dx,
                eyy: dv_dy,
                exy: 0.5 * (du_dy + dv_dx),
            };
            result.push((pt.x, pt.y, strain));
        }
        result
    }

    /// Average displacement magnitude across all (valid) points.
    pub fn mean_displacement(&self) -> f64 {
        if self.points.is_empty() {
            return 0.0;
        }
        let sum: f64 = self.points.iter().map(|p| p.displacement_magnitude()).sum();
        sum / self.points.len() as f64
    }

    /// Maximum displacement magnitude.
    pub fn max_displacement(&self) -> f64 {
        self.points
            .iter()
            .map(|p| p.displacement_magnitude())
            .fold(0.0_f64, f64::max)
    }
}

// ---------------------------------------------------------------------------
// PivDataReader — PIV velocity field import
// ---------------------------------------------------------------------------

/// A single velocity vector in a PIV field.
#[derive(Debug, Clone, PartialEq)]
pub struct PivVector {
    /// x position.
    pub x: f64,
    /// y position.
    pub y: f64,
    /// x velocity component.
    pub u: f64,
    /// v velocity component.
    pub v: f64,
    /// Signal-to-noise ratio.
    pub snr: f64,
}

impl PivVector {
    /// Create a new PIV vector.
    pub fn new(x: f64, y: f64, u: f64, v: f64, snr: f64) -> Self {
        Self { x, y, u, v, snr }
    }

    /// Velocity magnitude.
    pub fn speed(&self) -> f64 {
        (self.u * self.u + self.v * self.v).sqrt()
    }

    /// Vorticity (scalar, 2-D: ∂v/∂x − ∂u/∂y) requires neighbors.
    /// This method returns `(u, v)` as a convenience.
    pub fn components(&self) -> (f64, f64) {
        (self.u, self.v)
    }
}

/// Vector field statistics for PIV data.
#[derive(Debug, Clone)]
pub struct PivStatistics {
    /// Mean velocity magnitude.
    pub mean_speed: f64,
    /// RMS of the fluctuating velocity u'.
    pub u_rms: f64,
    /// RMS of the fluctuating velocity v'.
    pub v_rms: f64,
    /// Turbulence intensity (TI = sqrt((u_rms²+v_rms²)/2) / mean_speed).
    pub turbulence_intensity: f64,
    /// Mean u component.
    pub mean_u: f64,
    /// Mean v component.
    pub mean_v: f64,
}

/// PIV velocity field reader.
#[derive(Debug, Clone)]
pub struct PivDataReader {
    /// All velocity vectors.
    pub vectors: Vec<PivVector>,
    /// Minimum accepted SNR.
    pub min_snr: f64,
}

impl PivDataReader {
    /// Create an empty PIV reader.
    pub fn new(min_snr: f64) -> Self {
        Self {
            vectors: Vec::new(),
            min_snr,
        }
    }

    /// Parse from CSV: `x,y,u,v,snr` rows.
    pub fn from_csv(csv: &str, min_snr: f64) -> Result<Self, String> {
        let mut reader = Self::new(min_snr);
        for (i, line) in csv.lines().enumerate() {
            let line = line.trim();
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            let p: Vec<&str> = line.split(',').collect();
            if p.len() < 5 {
                return Err(format!("line {}: expected 5 columns", i + 1));
            }
            let parse = |s: &str| s.trim().parse::<f64>().map_err(|e| format!("{}", e));
            let x = parse(p[0])?;
            let y = parse(p[1])?;
            let u = parse(p[2])?;
            let v = parse(p[3])?;
            let snr = parse(p[4])?;
            reader.vectors.push(PivVector::new(x, y, u, v, snr));
        }
        Ok(reader)
    }

    /// Remove vectors with SNR below the threshold.
    pub fn filter_by_snr(&mut self) {
        self.vectors.retain(|v| v.snr >= self.min_snr);
    }

    /// Compute vector field statistics.
    pub fn statistics(&self) -> PivStatistics {
        let n = self.vectors.len();
        if n == 0 {
            return PivStatistics {
                mean_speed: 0.0,
                u_rms: 0.0,
                v_rms: 0.0,
                turbulence_intensity: 0.0,
                mean_u: 0.0,
                mean_v: 0.0,
            };
        }
        let nf = n as f64;
        let mean_u = self.vectors.iter().map(|v| v.u).sum::<f64>() / nf;
        let mean_v = self.vectors.iter().map(|v| v.v).sum::<f64>() / nf;
        let mean_speed = self.vectors.iter().map(|v| v.speed()).sum::<f64>() / nf;
        let u_rms = (self
            .vectors
            .iter()
            .map(|v| (v.u - mean_u).powi(2))
            .sum::<f64>()
            / nf)
            .sqrt();
        let v_rms = (self
            .vectors
            .iter()
            .map(|v| (v.v - mean_v).powi(2))
            .sum::<f64>()
            / nf)
            .sqrt();
        let ti = if mean_speed > 1e-12 {
            ((u_rms * u_rms + v_rms * v_rms) / 2.0).sqrt() / mean_speed
        } else {
            0.0
        };
        PivStatistics {
            mean_speed,
            u_rms,
            v_rms,
            turbulence_intensity: ti,
            mean_u,
            mean_v,
        }
    }

    /// Compute vorticity field using central differences on a regular grid.
    ///
    /// Returns a vector of `(x, y, vorticity)`.
    pub fn vorticity_field(&self, grid_spacing: f64) -> Vec<(f64, f64, f64)> {
        if self.vectors.len() < 4 || grid_spacing <= 0.0 {
            return Vec::new();
        }
        let mut map: HashMap<(i64, i64), &PivVector> = HashMap::new();
        for vec in &self.vectors {
            let ix = (vec.x / grid_spacing).round() as i64;
            let iy = (vec.y / grid_spacing).round() as i64;
            map.insert((ix, iy), vec);
        }
        let mut result = Vec::new();
        for (&(ix, iy), &pt) in &map {
            let r = map.get(&(ix + 1, iy));
            let l = map.get(&(ix - 1, iy));
            let u = map.get(&(ix, iy + 1));
            let d = map.get(&(ix, iy - 1));
            let dv_dx = match (r, l) {
                (Some(pr), Some(pl)) => (pr.v - pl.v) / (2.0 * grid_spacing),
                (Some(pr), None) => (pr.v - pt.v) / grid_spacing,
                (None, Some(pl)) => (pt.v - pl.v) / grid_spacing,
                _ => continue,
            };
            let du_dy = match (u, d) {
                (Some(pu), Some(pd)) => (pu.u - pd.u) / (2.0 * grid_spacing),
                (Some(pu), None) => (pu.u - pt.u) / grid_spacing,
                (None, Some(pd)) => (pt.u - pd.u) / grid_spacing,
                _ => continue,
            };
            result.push((pt.x, pt.y, dv_dx - du_dy));
        }
        result
    }
}

// ---------------------------------------------------------------------------
// StrainGaugeReader — rosette gauge data and principal strains
// ---------------------------------------------------------------------------

/// Rosette strain gauge reading: three gauges at 0°, 45°, 90°.
#[derive(Debug, Clone, PartialEq)]
pub struct RosetteReading {
    /// Strain at 0° gauge (micro-strain).
    pub e0: f64,
    /// Strain at 45° gauge (micro-strain).
    pub e45: f64,
    /// Strain at 90° gauge (micro-strain).
    pub e90: f64,
    /// Timestamp in seconds.
    pub time: f64,
}

impl RosetteReading {
    /// Create a new rosette reading.
    pub fn new(time: f64, e0: f64, e45: f64, e90: f64) -> Self {
        Self { e0, e45, e90, time }
    }

    /// Compute Cartesian strains (εxx, εyy, γxy) from rosette gauges.
    pub fn cartesian_strains(&self) -> (f64, f64, f64) {
        let exx = self.e0;
        let eyy = self.e90;
        let gamma_xy = 2.0 * self.e45 - self.e0 - self.e90;
        (exx, eyy, gamma_xy)
    }

    /// Compute principal strains (ε1, ε2) and principal angle θ (radians).
    pub fn principal_strains(&self) -> (f64, f64, f64) {
        let (exx, eyy, gamma_xy) = self.cartesian_strains();
        let avg = (exx + eyy) / 2.0;
        let r = ((exx - eyy).powi(2) / 4.0 + (gamma_xy / 2.0).powi(2)).sqrt();
        let e1 = avg + r;
        let e2 = avg - r;
        let theta = 0.5 * gamma_xy.atan2(exx - eyy);
        (e1, e2, theta)
    }

    /// Mohr circle: center and radius.
    pub fn mohr_circle(&self) -> (f64, f64) {
        let (exx, eyy, gamma_xy) = self.cartesian_strains();
        let center = (exx + eyy) / 2.0;
        let radius = ((exx - eyy).powi(2) / 4.0 + (gamma_xy / 2.0).powi(2)).sqrt();
        (center, radius)
    }

    /// Von Mises strain.
    pub fn von_mises(&self) -> f64 {
        let (e1, e2, _) = self.principal_strains();
        ((e1.powi(2) - e1 * e2 + e2.powi(2)) * 2.0 / 3.0).sqrt()
    }
}

/// Reader for strain gauge (rosette) data.
#[derive(Debug, Clone)]
pub struct StrainGaugeReader {
    /// All rosette readings in time order.
    pub readings: Vec<RosetteReading>,
    /// Gauge factor (default 2.0).
    pub gauge_factor: f64,
}

impl StrainGaugeReader {
    /// Create an empty reader.
    pub fn new(gauge_factor: f64) -> Self {
        Self {
            readings: Vec::new(),
            gauge_factor,
        }
    }

    /// Append a reading.
    pub fn push(&mut self, time: f64, e0: f64, e45: f64, e90: f64) {
        self.readings.push(RosetteReading::new(time, e0, e45, e90));
    }

    /// Parse from CSV: `time,e0,e45,e90` rows.
    pub fn from_csv(csv: &str, gauge_factor: f64) -> Result<Self, String> {
        let mut reader = Self::new(gauge_factor);
        for (i, line) in csv.lines().enumerate() {
            let line = line.trim();
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            let p: Vec<&str> = line.split(',').collect();
            if p.len() < 4 {
                return Err(format!("line {}: need 4 columns", i + 1));
            }
            let parse = |s: &str| s.trim().parse::<f64>().map_err(|e| format!("{}", e));
            let t = parse(p[0])?;
            let e0 = parse(p[1])?;
            let e45 = parse(p[2])?;
            let e90 = parse(p[3])?;
            reader.push(t, e0, e45, e90);
        }
        Ok(reader)
    }

    /// Compute peak principal strain across all readings.
    pub fn peak_principal_strain(&self) -> f64 {
        self.readings
            .iter()
            .map(|r| {
                let (e1, e2, _) = r.principal_strains();
                e1.abs().max(e2.abs())
            })
            .fold(0.0_f64, f64::max)
    }

    /// Mean von Mises strain.
    pub fn mean_von_mises(&self) -> f64 {
        if self.readings.is_empty() {
            return 0.0;
        }
        self.readings.iter().map(|r| r.von_mises()).sum::<f64>() / self.readings.len() as f64
    }

    /// Serialize to CSV.
    pub fn to_csv(&self) -> String {
        let mut out = String::from("# time,e0,e45,e90\n");
        for r in &self.readings {
            out.push_str(&format!("{},{},{},{}\n", r.time, r.e0, r.e45, r.e90));
        }
        out
    }
}

// ---------------------------------------------------------------------------
// DmaDataReader — DMA frequency sweep
// ---------------------------------------------------------------------------

/// A single DMA data point at a given frequency.
#[derive(Debug, Clone, PartialEq)]
pub struct DmaPoint {
    /// Frequency in Hz.
    pub frequency_hz: f64,
    /// Storage modulus E' (Pa).
    pub storage_modulus: f64,
    /// Loss modulus E'' (Pa).
    pub loss_modulus: f64,
    /// Temperature in °C.
    pub temperature_c: f64,
}

impl DmaPoint {
    /// Create a DMA point.
    pub fn new(
        frequency_hz: f64,
        storage_modulus: f64,
        loss_modulus: f64,
        temperature_c: f64,
    ) -> Self {
        Self {
            frequency_hz,
            storage_modulus,
            loss_modulus,
            temperature_c,
        }
    }

    /// Loss tangent (tan δ = E'' / E').
    pub fn tan_delta(&self) -> f64 {
        if self.storage_modulus.abs() < 1e-30 {
            return 0.0;
        }
        self.loss_modulus / self.storage_modulus
    }

    /// Complex modulus magnitude |E*| = √(E'² + E''²).
    pub fn complex_modulus_magnitude(&self) -> f64 {
        (self.storage_modulus.powi(2) + self.loss_modulus.powi(2)).sqrt()
    }
}

/// Master curve shift factor for time-temperature superposition.
#[derive(Debug, Clone)]
pub struct ShiftFactor {
    /// Temperature in °C.
    pub temperature_c: f64,
    /// Log₁₀ horizontal shift factor aT.
    pub log_at: f64,
}

/// DMA frequency sweep reader.
#[derive(Debug, Clone)]
pub struct DmaDataReader {
    /// DMA data points.
    pub points: Vec<DmaPoint>,
    /// Reference temperature for master curve (°C).
    pub reference_temperature_c: f64,
}

impl DmaDataReader {
    /// Create an empty DMA reader.
    pub fn new(reference_temperature_c: f64) -> Self {
        Self {
            points: Vec::new(),
            reference_temperature_c,
        }
    }

    /// Append a data point.
    pub fn push(
        &mut self,
        frequency_hz: f64,
        storage_modulus: f64,
        loss_modulus: f64,
        temperature_c: f64,
    ) {
        self.points.push(DmaPoint::new(
            frequency_hz,
            storage_modulus,
            loss_modulus,
            temperature_c,
        ));
    }

    /// Parse from CSV: `frequency_hz,storage_modulus,loss_modulus,temperature_c` rows.
    pub fn from_csv(csv: &str, reference_temperature_c: f64) -> Result<Self, String> {
        let mut reader = Self::new(reference_temperature_c);
        for (i, line) in csv.lines().enumerate() {
            let line = line.trim();
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            let p: Vec<&str> = line.split(',').collect();
            if p.len() < 4 {
                return Err(format!("line {}: need 4 columns", i + 1));
            }
            let parse = |s: &str| s.trim().parse::<f64>().map_err(|e| format!("{}", e));
            let f = parse(p[0])?;
            let ep = parse(p[1])?;
            let epp = parse(p[2])?;
            let t = parse(p[3])?;
            reader.push(f, ep, epp, t);
        }
        Ok(reader)
    }

    /// Mean tan δ across all points.
    pub fn mean_tan_delta(&self) -> f64 {
        if self.points.is_empty() {
            return 0.0;
        }
        self.points.iter().map(|p| p.tan_delta()).sum::<f64>() / self.points.len() as f64
    }

    /// Maximum storage modulus.
    pub fn max_storage_modulus(&self) -> f64 {
        self.points
            .iter()
            .map(|p| p.storage_modulus)
            .fold(0.0_f64, f64::max)
    }

    /// Compute WLF shift factors for master curve construction.
    ///
    /// Uses the WLF equation: log(aT) = -C1*(T - Tref) / (C2 + T - Tref).
    pub fn wlf_shift_factors(&self, c1: f64, c2: f64) -> Vec<ShiftFactor> {
        let tref = self.reference_temperature_c;
        let temps: Vec<f64> = {
            let mut ts: Vec<f64> = self.points.iter().map(|p| p.temperature_c).collect();
            ts.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
            ts.dedup_by(|a, b| (*a - *b).abs() < 0.01);
            ts
        };
        temps
            .iter()
            .map(|&t| {
                let dt = t - tref;
                let log_at = if (c2 + dt).abs() > 1e-12 {
                    -c1 * dt / (c2 + dt)
                } else {
                    0.0
                };
                ShiftFactor {
                    temperature_c: t,
                    log_at,
                }
            })
            .collect()
    }

    /// Build master curve: returns `(log10(reduced_freq), E')` pairs.
    pub fn master_curve_storage(&self, c1: f64, c2: f64) -> Vec<(f64, f64)> {
        let shifts = self.wlf_shift_factors(c1, c2);
        let shift_map: HashMap<i64, f64> = shifts
            .iter()
            .map(|sf| {
                let key = (sf.temperature_c * 100.0).round() as i64;
                (key, sf.log_at)
            })
            .collect();
        let mut curve: Vec<(f64, f64)> = self
            .points
            .iter()
            .map(|p| {
                let key = (p.temperature_c * 100.0).round() as i64;
                let log_at = shift_map.get(&key).copied().unwrap_or(0.0);
                let log_f_reduced = p.frequency_hz.log10() + log_at;
                (log_f_reduced, p.storage_modulus)
            })
            .collect();
        curve.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
        curve
    }

    /// Build master curve: returns `(log10(reduced_freq), tan_delta)` pairs.
    pub fn master_curve_tan_delta(&self, c1: f64, c2: f64) -> Vec<(f64, f64)> {
        let shifts = self.wlf_shift_factors(c1, c2);
        let shift_map: HashMap<i64, f64> = shifts
            .iter()
            .map(|sf| {
                let key = (sf.temperature_c * 100.0).round() as i64;
                (key, sf.log_at)
            })
            .collect();
        let mut curve: Vec<(f64, f64)> = self
            .points
            .iter()
            .map(|p| {
                let key = (p.temperature_c * 100.0).round() as i64;
                let log_at = shift_map.get(&key).copied().unwrap_or(0.0);
                let log_f_reduced = p.frequency_hz.log10() + log_at;
                (log_f_reduced, p.tan_delta())
            })
            .collect();
        curve.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
        curve
    }

    /// Filter points by temperature range \[t_min, t_max\].
    pub fn filter_by_temperature(&self, t_min: f64, t_max: f64) -> Vec<&DmaPoint> {
        self.points
            .iter()
            .filter(|p| p.temperature_c >= t_min && p.temperature_c <= t_max)
            .collect()
    }
}

// ---------------------------------------------------------------------------
// ExperimentalComparison — simulation vs experiment metrics
// ---------------------------------------------------------------------------

/// Comparison error type.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ComparisonError {
    /// Input slices have different lengths.
    LengthMismatch,
    /// Input slices are empty.
    EmptyData,
    /// Division by zero in normalization.
    ZeroDenominator,
}

impl fmt::Display for ComparisonError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::LengthMismatch => write!(f, "slice lengths do not match"),
            Self::EmptyData => write!(f, "empty data"),
            Self::ZeroDenominator => write!(f, "zero denominator in normalization"),
        }
    }
}

/// Collection of comparison metrics between a simulation and experiment.
#[derive(Debug, Clone)]
pub struct ComparisonMetrics {
    /// Root mean square error.
    pub rmse: f64,
    /// Normalized RMSE (RMSE / range of experimental data).
    pub nrmse: f64,
    /// Pearson correlation coefficient.
    pub correlation: f64,
    /// Mean absolute error.
    pub mae: f64,
    /// Maximum absolute error.
    pub max_error: f64,
    /// R² (coefficient of determination).
    pub r_squared: f64,
}

/// Simulation vs experiment comparison tool.
#[derive(Debug, Clone)]
pub struct ExperimentalComparison {
    /// Experimental (reference) values.
    pub experimental: Vec<f64>,
    /// Simulation (predicted) values.
    pub simulation: Vec<f64>,
}

impl ExperimentalComparison {
    /// Create a new comparison.
    pub fn new(experimental: Vec<f64>, simulation: Vec<f64>) -> Self {
        Self {
            experimental,
            simulation,
        }
    }

    /// Compute RMSE.
    pub fn rmse(&self) -> Result<f64, ComparisonError> {
        let (exp, sim) = self.check_lengths()?;
        let mse = exp
            .iter()
            .zip(sim)
            .map(|(e, s)| (e - s).powi(2))
            .sum::<f64>()
            / exp.len() as f64;
        Ok(mse.sqrt())
    }

    /// Compute normalized RMSE (divided by the experimental data range).
    pub fn nrmse(&self) -> Result<f64, ComparisonError> {
        let rmse = self.rmse()?;
        let min = self
            .experimental
            .iter()
            .cloned()
            .fold(f64::INFINITY, f64::min);
        let max = self
            .experimental
            .iter()
            .cloned()
            .fold(f64::NEG_INFINITY, f64::max);
        let range = max - min;
        if range.abs() < 1e-30 {
            return Err(ComparisonError::ZeroDenominator);
        }
        Ok(rmse / range)
    }

    /// Compute Pearson correlation coefficient.
    pub fn correlation(&self) -> Result<f64, ComparisonError> {
        let (exp, sim) = self.check_lengths()?;
        let n = exp.len() as f64;
        let me = exp.iter().sum::<f64>() / n;
        let ms = sim.iter().sum::<f64>() / n;
        let cov = exp
            .iter()
            .zip(sim)
            .map(|(e, s)| (e - me) * (s - ms))
            .sum::<f64>()
            / n;
        let std_e = (exp.iter().map(|e| (e - me).powi(2)).sum::<f64>() / n).sqrt();
        let std_s = (sim.iter().map(|s| (s - ms).powi(2)).sum::<f64>() / n).sqrt();
        if std_e < 1e-30 || std_s < 1e-30 {
            return Err(ComparisonError::ZeroDenominator);
        }
        Ok(cov / (std_e * std_s))
    }

    /// Compute mean absolute error.
    pub fn mae(&self) -> Result<f64, ComparisonError> {
        let (exp, sim) = self.check_lengths()?;
        Ok(exp.iter().zip(sim).map(|(e, s)| (e - s).abs()).sum::<f64>() / exp.len() as f64)
    }

    /// Compute maximum absolute error.
    pub fn max_error(&self) -> Result<f64, ComparisonError> {
        let (exp, sim) = self.check_lengths()?;
        Ok(exp
            .iter()
            .zip(sim)
            .map(|(e, s)| (e - s).abs())
            .fold(0.0_f64, f64::max))
    }

    /// Compute R² (coefficient of determination).
    pub fn r_squared(&self) -> Result<f64, ComparisonError> {
        let (exp, sim) = self.check_lengths()?;
        let n = exp.len() as f64;
        let me = exp.iter().sum::<f64>() / n;
        let ss_tot = exp.iter().map(|e| (e - me).powi(2)).sum::<f64>();
        let ss_res = exp
            .iter()
            .zip(sim)
            .map(|(e, s)| (e - s).powi(2))
            .sum::<f64>();
        if ss_tot.abs() < 1e-30 {
            return Err(ComparisonError::ZeroDenominator);
        }
        Ok(1.0 - ss_res / ss_tot)
    }

    /// Compute all metrics at once.
    pub fn all_metrics(&self) -> Result<ComparisonMetrics, ComparisonError> {
        Ok(ComparisonMetrics {
            rmse: self.rmse()?,
            nrmse: self.nrmse().unwrap_or(0.0),
            correlation: self.correlation().unwrap_or(0.0),
            mae: self.mae()?,
            max_error: self.max_error()?,
            r_squared: self.r_squared().unwrap_or(0.0),
        })
    }

    fn check_lengths(&self) -> Result<(&[f64], &[f64]), ComparisonError> {
        if self.experimental.is_empty() || self.simulation.is_empty() {
            return Err(ComparisonError::EmptyData);
        }
        if self.experimental.len() != self.simulation.len() {
            return Err(ComparisonError::LengthMismatch);
        }
        Ok((&self.experimental, &self.simulation))
    }

    /// Bias (mean of simulation - experiment).
    pub fn bias(&self) -> Result<f64, ComparisonError> {
        let (exp, sim) = self.check_lengths()?;
        let n = exp.len() as f64;
        Ok(sim.iter().zip(exp).map(|(s, e)| s - e).sum::<f64>() / n)
    }

    /// Relative error for each point: (sim - exp) / exp.
    pub fn relative_errors(&self) -> Result<Vec<f64>, ComparisonError> {
        let (exp, sim) = self.check_lengths()?;
        let v: Vec<f64> = exp
            .iter()
            .zip(sim)
            .map(|(e, s)| if e.abs() < 1e-30 { 0.0 } else { (s - e) / e })
            .collect();
        Ok(v)
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // --- ExperimentalDataset ---

    #[test]
    fn test_dataset_push_and_len() {
        let meta = SensorMetadata::new("S1", "force", "N", 100.0);
        let mut ds = ExperimentalDataset::new(meta);
        ds.push(0.0, 1.0);
        ds.push(0.01, 2.0);
        assert_eq!(ds.len(), 2);
        assert!(!ds.is_empty());
    }

    #[test]
    fn test_dataset_mean() {
        let meta = SensorMetadata::new("S1", "force", "N", 100.0);
        let mut ds = ExperimentalDataset::new(meta);
        ds.push(0.0, 1.0);
        ds.push(0.1, 3.0);
        assert!((ds.mean() - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_dataset_std_dev() {
        let meta = SensorMetadata::new("S1", "displacement", "mm", 50.0);
        let mut ds = ExperimentalDataset::new(meta);
        ds.push(0.0, 2.0);
        ds.push(0.1, 4.0);
        let sd = ds.std_dev();
        assert!((sd - std::f64::consts::SQRT_2).abs() < 1e-10);
    }

    #[test]
    fn test_dataset_flag_outliers_zscore() {
        let meta = SensorMetadata::new("S1", "temp", "C", 1.0);
        let mut ds = ExperimentalDataset::new(meta);
        for i in 0..20 {
            ds.push(i as f64, i as f64);
        }
        ds.push(20.0, 1000.0); // outlier
        ds.flag_outliers_zscore(2.5);
        let outliers = ds.outliers();
        assert!(!outliers.is_empty());
    }

    #[test]
    fn test_dataset_flag_outliers_iqr() {
        let meta = SensorMetadata::new("S1", "pressure", "Pa", 10.0);
        let mut ds = ExperimentalDataset::new(meta);
        for i in 0..20 {
            ds.push(i as f64, i as f64);
        }
        ds.push(20.0, 9999.0);
        ds.flag_outliers_iqr(1.5);
        assert!(!ds.outliers().is_empty());
    }

    #[test]
    fn test_dataset_from_csv() {
        let csv = "0.0,1.0\n0.1,2.0\n0.2,3.0\n";
        let meta = SensorMetadata::new("S1", "x", "m", 10.0);
        let ds = ExperimentalDataset::from_csv(csv, meta).unwrap();
        assert_eq!(ds.len(), 3);
        assert!((ds.mean() - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_dataset_to_csv_roundtrip() {
        let meta = SensorMetadata::new("S1", "x", "m", 10.0);
        let mut ds = ExperimentalDataset::new(meta.clone());
        ds.push(0.0, 1.5);
        ds.push(1.0, 2.5);
        let csv = ds.to_csv();
        let ds2 = ExperimentalDataset::from_csv(
            &csv.lines()
                .filter(|l| !l.starts_with('#'))
                .collect::<Vec<_>>()
                .join("\n"),
            meta,
        )
        .unwrap();
        assert_eq!(ds2.len(), 2);
    }

    #[test]
    fn test_dataset_range() {
        let meta = SensorMetadata::new("S1", "x", "m", 10.0);
        let mut ds = ExperimentalDataset::new(meta);
        ds.push(0.0, 5.0);
        ds.push(0.1, 10.0);
        ds.push(0.2, 3.0);
        let (min, max) = ds.range().unwrap();
        assert!((min - 3.0).abs() < 1e-10);
        assert!((max - 10.0).abs() < 1e-10);
    }

    #[test]
    fn test_dataset_downsample() {
        let meta = SensorMetadata::new("S1", "x", "m", 100.0);
        let mut ds = ExperimentalDataset::new(meta);
        for i in 0..10 {
            ds.push(i as f64, i as f64);
        }
        let ds2 = ds.downsample(2);
        assert_eq!(ds2.len(), 5);
    }

    #[test]
    fn test_sensor_metadata_builder() {
        let meta = SensorMetadata::new("A1", "strain", "με", 1000.0)
            .with_location("beam top")
            .with_extra("calibration", "2026-01");
        assert_eq!(meta.location.unwrap(), "beam top");
        assert_eq!(meta.extra["calibration"], "2026-01");
    }

    // --- DIC ---

    #[test]
    fn test_dic_point_displacement() {
        let p = DicPoint::new(0.0, 0.0, 3.0, 4.0, 0.95);
        assert!((p.displacement_magnitude() - 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_dic_from_csv() {
        let csv = "0,0,1,2,0.9\n1,0,3,4,0.85\n";
        let dic = DigitalImageCorrelation::from_csv(csv, 0.8).unwrap();
        assert_eq!(dic.points.len(), 2);
    }

    #[test]
    fn test_dic_filter_by_correlation() {
        let csv = "0,0,1,2,0.9\n1,0,3,4,0.3\n";
        let mut dic = DigitalImageCorrelation::from_csv(csv, 0.8).unwrap();
        dic.filter_by_correlation();
        assert_eq!(dic.points.len(), 1);
    }

    #[test]
    fn test_dic_mean_displacement() {
        let mut dic = DigitalImageCorrelation::new(0.5);
        dic.points.push(DicPoint::new(0.0, 0.0, 3.0, 4.0, 0.9)); // magnitude 5
        dic.points.push(DicPoint::new(1.0, 0.0, 0.0, 0.0, 0.9)); // magnitude 0
        assert!((dic.mean_displacement() - 2.5).abs() < 1e-10);
    }

    #[test]
    fn test_dic_max_displacement() {
        let mut dic = DigitalImageCorrelation::new(0.5);
        dic.points.push(DicPoint::new(0.0, 0.0, 3.0, 4.0, 0.9));
        dic.points.push(DicPoint::new(1.0, 0.0, 1.0, 0.0, 0.9));
        assert!((dic.max_displacement() - 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_dic_compute_strain_field() {
        let mut dic = DigitalImageCorrelation::new(0.5);
        let gs = 1.0;
        // Create a 3x3 grid so central-difference neighbors exist in both axes
        for iy in 0..3i64 {
            for ix in 0..3i64 {
                let x = ix as f64 * gs;
                let y = iy as f64 * gs;
                dic.points
                    .push(DicPoint::new(x, y, x * 0.01, y * 0.005, 0.9));
            }
        }
        let strains = dic.compute_strain_field(gs);
        assert!(!strains.is_empty());
    }

    #[test]
    fn test_strain_tensor_principal() {
        let st = StrainTensor2D {
            exx: 0.001,
            eyy: 0.0,
            exy: 0.0,
        };
        let (e1, e2, _) = st.principal_strains();
        assert!(e1 > e2);
    }

    #[test]
    fn test_strain_tensor_von_mises() {
        let st = StrainTensor2D {
            exx: 0.001,
            eyy: 0.0,
            exy: 0.0,
        };
        assert!(st.von_mises() > 0.0);
    }

    // --- PIV ---

    #[test]
    fn test_piv_from_csv() {
        let csv = "0,0,1.0,0.5,10.0\n1,0,1.5,0.3,12.0\n";
        let piv = PivDataReader::from_csv(csv, 5.0).unwrap();
        assert_eq!(piv.vectors.len(), 2);
    }

    #[test]
    fn test_piv_filter_snr() {
        let csv = "0,0,1.0,0.5,10.0\n1,0,1.5,0.3,2.0\n";
        let mut piv = PivDataReader::from_csv(csv, 5.0).unwrap();
        piv.filter_by_snr();
        assert_eq!(piv.vectors.len(), 1);
    }

    #[test]
    fn test_piv_statistics_basic() {
        let mut piv = PivDataReader::new(0.0);
        piv.vectors.push(PivVector::new(0.0, 0.0, 1.0, 0.0, 10.0));
        piv.vectors.push(PivVector::new(1.0, 0.0, 1.0, 0.0, 10.0));
        let stats = piv.statistics();
        assert!((stats.mean_u - 1.0).abs() < 1e-10);
        assert!((stats.u_rms).abs() < 1e-10);
    }

    #[test]
    fn test_piv_turbulence_intensity() {
        let mut piv = PivDataReader::new(0.0);
        piv.vectors.push(PivVector::new(0.0, 0.0, 1.0, 0.0, 10.0));
        piv.vectors.push(PivVector::new(1.0, 0.0, 2.0, 0.0, 10.0));
        let stats = piv.statistics();
        assert!(stats.turbulence_intensity >= 0.0);
    }

    #[test]
    fn test_piv_vorticity_field() {
        let mut piv = PivDataReader::new(0.0);
        let gs = 1.0;
        for ix in 0..3i64 {
            for iy in 0..3i64 {
                let x = ix as f64 * gs;
                let y = iy as f64 * gs;
                piv.vectors
                    .push(PivVector::new(x, y, -y * 0.1, x * 0.1, 10.0));
            }
        }
        let vort = piv.vorticity_field(gs);
        assert!(!vort.is_empty());
    }

    #[test]
    fn test_piv_speed() {
        let v = PivVector::new(0.0, 0.0, 3.0, 4.0, 10.0);
        assert!((v.speed() - 5.0).abs() < 1e-10);
    }

    // --- StrainGauge ---

    #[test]
    fn test_rosette_cartesian_strains() {
        // e0=100, e45=50, e90=0 → γxy = 2*50 - 100 - 0 = 0
        let r = RosetteReading::new(0.0, 100.0, 50.0, 0.0);
        let (exx, eyy, gxy) = r.cartesian_strains();
        assert!((exx - 100.0).abs() < 1e-10);
        assert!((eyy - 0.0).abs() < 1e-10);
        assert!((gxy - 0.0).abs() < 1e-10);
    }

    #[test]
    fn test_rosette_principal_strains() {
        let r = RosetteReading::new(0.0, 100.0, 0.0, 0.0);
        let (e1, e2, _) = r.principal_strains();
        assert!(e1 >= e2);
    }

    #[test]
    fn test_rosette_mohr_circle() {
        let r = RosetteReading::new(0.0, 100.0, 0.0, 0.0);
        let (center, radius) = r.mohr_circle();
        assert!(radius >= 0.0);
        assert!(center.is_finite());
    }

    #[test]
    fn test_strain_gauge_from_csv() {
        let csv = "0.0,100,50,10\n0.1,200,80,20\n";
        let sg = StrainGaugeReader::from_csv(csv, 2.0).unwrap();
        assert_eq!(sg.readings.len(), 2);
    }

    #[test]
    fn test_strain_gauge_peak() {
        let mut sg = StrainGaugeReader::new(2.0);
        sg.push(0.0, 100.0, 0.0, 0.0);
        sg.push(1.0, 500.0, 0.0, 0.0);
        assert!(sg.peak_principal_strain() > 0.0);
    }

    #[test]
    fn test_strain_gauge_mean_von_mises() {
        let mut sg = StrainGaugeReader::new(2.0);
        sg.push(0.0, 100.0, 50.0, 0.0);
        let vm = sg.mean_von_mises();
        assert!(vm >= 0.0);
    }

    #[test]
    fn test_strain_gauge_to_csv() {
        let mut sg = StrainGaugeReader::new(2.0);
        sg.push(0.0, 100.0, 50.0, 10.0);
        let csv = sg.to_csv();
        assert!(csv.contains("100"));
    }

    // --- DMA ---

    #[test]
    fn test_dma_tan_delta() {
        let p = DmaPoint::new(1.0, 1e6, 2e5, 25.0);
        assert!((p.tan_delta() - 0.2).abs() < 1e-10);
    }

    #[test]
    fn test_dma_complex_modulus() {
        let p = DmaPoint::new(1.0, 3.0, 4.0, 25.0);
        assert!((p.complex_modulus_magnitude() - 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_dma_from_csv() {
        let csv = "1.0,1e6,2e5,25.0\n10.0,1.1e6,2.1e5,25.0\n";
        let dma = DmaDataReader::from_csv(csv, 25.0).unwrap();
        assert_eq!(dma.points.len(), 2);
    }

    #[test]
    fn test_dma_mean_tan_delta() {
        let mut dma = DmaDataReader::new(25.0);
        dma.push(1.0, 1e6, 2e5, 25.0);
        dma.push(10.0, 1e6, 3e5, 25.0);
        let mtd = dma.mean_tan_delta();
        assert!(mtd > 0.0);
    }

    #[test]
    fn test_dma_max_storage_modulus() {
        let mut dma = DmaDataReader::new(25.0);
        dma.push(1.0, 1e6, 2e5, 25.0);
        dma.push(10.0, 2e6, 3e5, 25.0);
        assert!((dma.max_storage_modulus() - 2e6).abs() < 1.0);
    }

    #[test]
    fn test_dma_wlf_shift_factors() {
        let mut dma = DmaDataReader::new(25.0);
        dma.push(1.0, 1e6, 2e5, 0.0);
        dma.push(1.0, 1e6, 2e5, 25.0);
        dma.push(1.0, 1e6, 2e5, 50.0);
        let shifts = dma.wlf_shift_factors(8.86, 101.6);
        assert_eq!(shifts.len(), 3);
    }

    #[test]
    fn test_dma_master_curve_storage() {
        let mut dma = DmaDataReader::new(25.0);
        dma.push(1.0, 1e6, 2e5, 0.0);
        dma.push(10.0, 1.1e6, 2.1e5, 25.0);
        let curve = dma.master_curve_storage(8.86, 101.6);
        assert!(!curve.is_empty());
    }

    #[test]
    fn test_dma_master_curve_tan_delta() {
        let mut dma = DmaDataReader::new(25.0);
        dma.push(1.0, 1e6, 2e5, 25.0);
        dma.push(10.0, 1.2e6, 2.4e5, 40.0);
        let curve = dma.master_curve_tan_delta(8.86, 101.6);
        assert!(!curve.is_empty());
        for (_, td) in &curve {
            assert!(td.is_finite());
        }
    }

    #[test]
    fn test_dma_filter_by_temperature() {
        let mut dma = DmaDataReader::new(25.0);
        dma.push(1.0, 1e6, 2e5, 10.0);
        dma.push(1.0, 1e6, 2e5, 30.0);
        dma.push(1.0, 1e6, 2e5, 60.0);
        let filtered = dma.filter_by_temperature(0.0, 40.0);
        assert_eq!(filtered.len(), 2);
    }

    // --- ExperimentalComparison ---

    #[test]
    fn test_comparison_rmse() {
        let cmp = ExperimentalComparison::new(vec![1.0, 2.0, 3.0], vec![1.0, 2.0, 4.0]);
        let rmse = cmp.rmse().unwrap();
        // error at last = 1.0, rmse = sqrt(1/3)
        assert!((rmse - (1.0_f64 / 3.0).sqrt()).abs() < 1e-10);
    }

    #[test]
    fn test_comparison_mae() {
        let cmp = ExperimentalComparison::new(vec![0.0, 1.0, 2.0], vec![0.0, 2.0, 2.0]);
        let mae = cmp.mae().unwrap();
        assert!((mae - 1.0 / 3.0).abs() < 1e-10);
    }

    #[test]
    fn test_comparison_correlation_perfect() {
        let exp: Vec<f64> = (0..10).map(|i| i as f64).collect();
        let sim = exp.clone();
        let cmp = ExperimentalComparison::new(exp, sim);
        let r = cmp.correlation().unwrap();
        assert!((r - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_comparison_r_squared_perfect() {
        let exp: Vec<f64> = (0..10).map(|i| i as f64).collect();
        let sim = exp.clone();
        let cmp = ExperimentalComparison::new(exp, sim);
        let r2 = cmp.r_squared().unwrap();
        assert!((r2 - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_comparison_length_mismatch() {
        let cmp = ExperimentalComparison::new(vec![1.0, 2.0], vec![1.0]);
        assert_eq!(cmp.rmse(), Err(ComparisonError::LengthMismatch));
    }

    #[test]
    fn test_comparison_empty() {
        let cmp = ExperimentalComparison::new(vec![], vec![]);
        assert_eq!(cmp.rmse(), Err(ComparisonError::EmptyData));
    }

    #[test]
    fn test_comparison_max_error() {
        let cmp = ExperimentalComparison::new(vec![0.0, 1.0, 2.0], vec![0.0, 1.5, 2.0]);
        let me = cmp.max_error().unwrap();
        assert!((me - 0.5).abs() < 1e-10);
    }

    #[test]
    fn test_comparison_bias() {
        let cmp = ExperimentalComparison::new(vec![0.0, 1.0], vec![1.0, 2.0]);
        let b = cmp.bias().unwrap();
        assert!((b - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_comparison_relative_errors() {
        let cmp = ExperimentalComparison::new(vec![1.0, 2.0], vec![1.1, 2.2]);
        let re = cmp.relative_errors().unwrap();
        assert!((re[0] - 0.1).abs() < 1e-10);
        assert!((re[1] - 0.1).abs() < 1e-10);
    }

    #[test]
    fn test_comparison_nrmse() {
        let cmp = ExperimentalComparison::new(vec![0.0, 10.0], vec![1.0, 9.0]);
        let nrmse = cmp.nrmse().unwrap();
        assert!(nrmse.is_finite() && nrmse > 0.0);
    }

    #[test]
    fn test_comparison_all_metrics() {
        let exp: Vec<f64> = (0..5).map(|i| i as f64).collect();
        let sim: Vec<f64> = exp.iter().map(|&v| v + 0.1).collect();
        let cmp = ExperimentalComparison::new(exp, sim);
        let m = cmp.all_metrics().unwrap();
        assert!(m.rmse > 0.0);
        assert!(m.mae > 0.0);
        assert!(m.correlation > 0.0);
    }

    #[test]
    fn test_percentile() {
        let data: Vec<f64> = (1..=10).map(|i| i as f64).collect();
        let p50 = percentile(&data, 50.0);
        assert!((p50 - 5.5).abs() < 0.1);
    }

    #[test]
    fn test_piv_empty_statistics() {
        let piv = PivDataReader::new(0.0);
        let stats = piv.statistics();
        assert_eq!(stats.mean_speed, 0.0);
    }

    #[test]
    fn test_dic_empty_strain_field() {
        let dic = DigitalImageCorrelation::new(0.5);
        let strains = dic.compute_strain_field(1.0);
        assert!(strains.is_empty());
    }

    #[test]
    fn test_dataset_empty_mean() {
        let meta = SensorMetadata::new("S", "x", "m", 1.0);
        let ds = ExperimentalDataset::new(meta);
        assert_eq!(ds.mean(), 0.0);
    }

    #[test]
    fn test_dataset_empty_range() {
        let meta = SensorMetadata::new("S", "x", "m", 1.0);
        let ds = ExperimentalDataset::new(meta);
        assert!(ds.range().is_none());
    }

    #[test]
    fn test_rosette_von_mises_nonneg() {
        let r = RosetteReading::new(0.0, 200.0, 100.0, 50.0);
        assert!(r.von_mises() >= 0.0);
    }

    #[test]
    fn test_comparison_error_display() {
        let e = ComparisonError::LengthMismatch;
        assert!(!format!("{}", e).is_empty());
    }
}
