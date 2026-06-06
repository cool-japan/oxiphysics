// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Sensor data I/O for the OxiPhysics engine.
//!
//! Provides structured types and serialization/deserialization helpers for
//! common physical sensor data streams including:
//!
//! - **IMU** (accelerometer, gyroscope, magnetometer)
//! - **Pressure sensor** time series
//! - **Temperature sensor array**
//! - **Strain gauge** data
//! - **Force/torque sensor** (6-DOF)
//! - **Optical encoder** (position/velocity)
//! - **LVDT displacement** sensor
//! - **Thermocouple calibration** tables
//! - **Sensor fusion** data (Kalman-filtered state)
//! - **Calibration matrix** storage and retrieval

// ─────────────────────────────────────────────────────────────────────────────
// Common helpers
// ─────────────────────────────────────────────────────────────────────────────

/// Read a little-endian `f64` from a byte slice at the given offset.
///
/// Returns `None` when fewer than 8 bytes remain.
fn read_f64_le(data: &[u8], offset: usize) -> Option<f64> {
    if offset + 8 > data.len() {
        return None;
    }
    let arr: [u8; 8] = data[offset..offset + 8].try_into().ok()?;
    Some(f64::from_le_bytes(arr))
}

/// Append the little-endian encoding of `v` to `buf`.
fn push_f64_le(buf: &mut Vec<u8>, v: f64) {
    buf.extend_from_slice(&v.to_le_bytes());
}

// ─────────────────────────────────────────────────────────────────────────────
// ImuSample — inertial measurement unit
// ─────────────────────────────────────────────────────────────────────────────

/// A single sample from a 9-axis IMU.
///
/// All coordinate frames are body-fixed (right-hand rule, X forward, Y left,
/// Z up) unless otherwise stated.
#[derive(Debug, Clone, PartialEq)]
pub struct ImuSample {
    /// Timestamp in seconds.
    pub timestamp: f64,
    /// Accelerometer reading `[ax, ay, az]` in m/s².
    pub accel: [f64; 3],
    /// Gyroscope reading `[ωx, ωy, ωz]` in rad/s.
    pub gyro: [f64; 3],
    /// Magnetometer reading `[bx, by, bz]` in µT.
    pub mag: [f64; 3],
}

impl ImuSample {
    /// Construct a new sample from component arrays.
    pub fn new(timestamp: f64, accel: [f64; 3], gyro: [f64; 3], mag: [f64; 3]) -> Self {
        Self {
            timestamp,
            accel,
            gyro,
            mag,
        }
    }

    /// Serialize to 80 bytes of little-endian `f64`.
    ///
    /// Layout: `[t, ax, ay, az, ωx, ωy, ωz, bx, by, bz]` (10 × 8 bytes).
    pub fn to_bytes(&self) -> Vec<u8> {
        let mut buf = Vec::with_capacity(80);
        push_f64_le(&mut buf, self.timestamp);
        for &v in &self.accel {
            push_f64_le(&mut buf, v);
        }
        for &v in &self.gyro {
            push_f64_le(&mut buf, v);
        }
        for &v in &self.mag {
            push_f64_le(&mut buf, v);
        }
        buf
    }

    /// Deserialize from 80 bytes.  Returns `None` on truncated input.
    pub fn from_bytes(data: &[u8]) -> Option<Self> {
        if data.len() < 80 {
            return None;
        }
        let mut o = 0usize;
        let mut next = || {
            let v = read_f64_le(data, o);
            o += 8;
            v
        };
        Some(Self {
            timestamp: next()?,
            accel: [next()?, next()?, next()?],
            gyro: [next()?, next()?, next()?],
            mag: [next()?, next()?, next()?],
        })
    }

    /// Magnitude of the accelerometer vector in m/s².
    pub fn accel_magnitude(&self) -> f64 {
        let [ax, ay, az] = self.accel;
        (ax * ax + ay * ay + az * az).sqrt()
    }

    /// Magnitude of the angular-rate vector in rad/s.
    pub fn gyro_magnitude(&self) -> f64 {
        let [gx, gy, gz] = self.gyro;
        (gx * gx + gy * gy + gz * gz).sqrt()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// ImuStream — buffered IMU log
// ─────────────────────────────────────────────────────────────────────────────

/// A buffered sequence of [`ImuSample`] records with basic statistics.
#[derive(Debug, Clone, Default)]
pub struct ImuStream {
    /// Raw samples in chronological order.
    pub samples: Vec<ImuSample>,
}

impl ImuStream {
    /// Create an empty stream.
    pub fn new() -> Self {
        Self {
            samples: Vec::new(),
        }
    }

    /// Append a sample.
    pub fn push(&mut self, sample: ImuSample) {
        self.samples.push(sample);
    }

    /// Return the number of samples.
    pub fn len(&self) -> usize {
        self.samples.len()
    }

    /// Return `true` when the stream contains no samples.
    pub fn is_empty(&self) -> bool {
        self.samples.is_empty()
    }

    /// Compute the mean accelerometer vector over all samples.
    ///
    /// Returns `[0,0,0]` for an empty stream.
    pub fn mean_accel(&self) -> [f64; 3] {
        if self.samples.is_empty() {
            return [0.0; 3];
        }
        let n = self.samples.len() as f64;
        let mut sum = [0.0f64; 3];
        for s in &self.samples {
            for (acc, &a) in sum.iter_mut().zip(s.accel.iter()) {
                *acc += a;
            }
        }
        [sum[0] / n, sum[1] / n, sum[2] / n]
    }

    /// Compute the mean gyroscope vector over all samples.
    pub fn mean_gyro(&self) -> [f64; 3] {
        if self.samples.is_empty() {
            return [0.0; 3];
        }
        let n = self.samples.len() as f64;
        let mut sum = [0.0f64; 3];
        for s in &self.samples {
            for (acc, &g) in sum.iter_mut().zip(s.gyro.iter()) {
                *acc += g;
            }
        }
        [sum[0] / n, sum[1] / n, sum[2] / n]
    }

    /// Serialize the entire stream to a byte vector.
    ///
    /// Format: 8-byte little-endian `u64` count, followed by concatenated
    /// 80-byte sample records.
    pub fn to_bytes(&self) -> Vec<u8> {
        let n = self.samples.len() as u64;
        let mut buf = Vec::with_capacity(8 + self.samples.len() * 80);
        buf.extend_from_slice(&n.to_le_bytes());
        for s in &self.samples {
            buf.extend_from_slice(&s.to_bytes());
        }
        buf
    }

    /// Deserialize a stream from bytes produced by [`ImuStream::to_bytes`].
    pub fn from_bytes(data: &[u8]) -> Option<Self> {
        if data.len() < 8 {
            return None;
        }
        let n = u64::from_le_bytes(data[0..8].try_into().ok()?) as usize;
        if data.len() < 8 + n * 80 {
            return None;
        }
        let mut samples = Vec::with_capacity(n);
        for i in 0..n {
            samples.push(ImuSample::from_bytes(&data[8 + i * 80..])?);
        }
        Some(Self { samples })
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// PressureSample — barometric / differential / absolute pressure sensor
// ─────────────────────────────────────────────────────────────────────────────

/// A single pressure sensor reading with timestamp.
#[derive(Debug, Clone, PartialEq)]
pub struct PressureSample {
    /// Timestamp in seconds.
    pub timestamp: f64,
    /// Pressure reading in Pascal.
    pub pressure_pa: f64,
    /// Temperature of the pressure sensor die in °C.
    pub temperature_c: f64,
}

impl PressureSample {
    /// Create a new sample.
    pub fn new(timestamp: f64, pressure_pa: f64, temperature_c: f64) -> Self {
        Self {
            timestamp,
            pressure_pa,
            temperature_c,
        }
    }

    /// Altitude estimate from pressure using the international barometric formula.
    ///
    /// Reference pressure `p0` in Pa and reference altitude `h0` in metres.
    pub fn altitude_m(&self, p0: f64, h0: f64) -> f64 {
        const T0: f64 = 288.15; // sea-level temperature K
        const L: f64 = 0.0065; // lapse rate K/m
        const R: f64 = 8.3144598; // J/(mol·K)
        const M: f64 = 0.0289644; // kg/mol
        const G: f64 = 9.80665; // m/s²
        let exp = R * L / (G * M);
        h0 + (T0 / L) * (1.0 - (self.pressure_pa / p0).powf(exp))
    }

    /// Serialize to 24 bytes.
    pub fn to_bytes(&self) -> [u8; 24] {
        let mut buf = [0u8; 24];
        buf[0..8].copy_from_slice(&self.timestamp.to_le_bytes());
        buf[8..16].copy_from_slice(&self.pressure_pa.to_le_bytes());
        buf[16..24].copy_from_slice(&self.temperature_c.to_le_bytes());
        buf
    }

    /// Deserialize from 24 bytes.
    pub fn from_bytes(data: &[u8]) -> Option<Self> {
        if data.len() < 24 {
            return None;
        }
        Some(Self {
            timestamp: read_f64_le(data, 0)?,
            pressure_pa: read_f64_le(data, 8)?,
            temperature_c: read_f64_le(data, 16)?,
        })
    }
}

/// A time-ordered series of pressure samples.
#[derive(Debug, Clone, Default)]
pub struct PressureTimeSeries {
    /// The underlying samples.
    pub samples: Vec<PressureSample>,
}

impl PressureTimeSeries {
    /// Create an empty series.
    pub fn new() -> Self {
        Self {
            samples: Vec::new(),
        }
    }

    /// Append a sample.
    pub fn push(&mut self, s: PressureSample) {
        self.samples.push(s);
    }

    /// Minimum pressure over the series in Pa, or `f64::NAN` if empty.
    pub fn min_pressure(&self) -> f64 {
        self.samples
            .iter()
            .map(|s| s.pressure_pa)
            .fold(f64::INFINITY, f64::min)
    }

    /// Maximum pressure over the series in Pa, or `f64::NAN` if empty.
    pub fn max_pressure(&self) -> f64 {
        self.samples
            .iter()
            .map(|s| s.pressure_pa)
            .fold(f64::NEG_INFINITY, f64::max)
    }

    /// Mean pressure in Pa.
    pub fn mean_pressure(&self) -> f64 {
        if self.samples.is_empty() {
            return 0.0;
        }
        let sum: f64 = self.samples.iter().map(|s| s.pressure_pa).sum();
        sum / self.samples.len() as f64
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// TemperatureSensorArray — spatially distributed temperature probes
// ─────────────────────────────────────────────────────────────────────────────

/// A reading from one probe in a distributed temperature array.
#[derive(Debug, Clone, PartialEq)]
pub struct TempProbeReading {
    /// Probe identifier (0-based).
    pub probe_id: u32,
    /// Timestamp in seconds.
    pub timestamp: f64,
    /// Temperature in °C.
    pub temperature_c: f64,
    /// 3-D position of the probe `[x, y, z]` in metres.
    pub position: [f64; 3],
}

impl TempProbeReading {
    /// Create a new reading.
    pub fn new(probe_id: u32, timestamp: f64, temperature_c: f64, position: [f64; 3]) -> Self {
        Self {
            probe_id,
            timestamp,
            temperature_c,
            position,
        }
    }

    /// Serialize to 44 bytes (4 + 8 + 8 + 3×8).
    pub fn to_bytes(&self) -> Vec<u8> {
        let mut buf = Vec::with_capacity(44);
        buf.extend_from_slice(&self.probe_id.to_le_bytes());
        push_f64_le(&mut buf, self.timestamp);
        push_f64_le(&mut buf, self.temperature_c);
        for &v in &self.position {
            push_f64_le(&mut buf, v);
        }
        buf
    }

    /// Deserialize from bytes.
    pub fn from_bytes(data: &[u8]) -> Option<Self> {
        if data.len() < 44 {
            return None;
        }
        let probe_id = u32::from_le_bytes(data[0..4].try_into().ok()?);
        Some(Self {
            probe_id,
            timestamp: read_f64_le(data, 4)?,
            temperature_c: read_f64_le(data, 12)?,
            position: [
                read_f64_le(data, 20)?,
                read_f64_le(data, 28)?,
                read_f64_le(data, 36)?,
            ],
        })
    }
}

/// Snapshot of all probes in a temperature array at one instant.
#[derive(Debug, Clone, Default)]
pub struct TempArraySnapshot {
    /// Timestamp for this snapshot.
    pub timestamp: f64,
    /// Temperature readings, one per probe.
    pub readings: Vec<f64>,
}

impl TempArraySnapshot {
    /// Create a new snapshot.
    pub fn new(timestamp: f64, readings: Vec<f64>) -> Self {
        Self {
            timestamp,
            readings,
        }
    }

    /// Number of probes.
    pub fn num_probes(&self) -> usize {
        self.readings.len()
    }

    /// Maximum temperature in the array.
    pub fn max_temp(&self) -> f64 {
        self.readings
            .iter()
            .cloned()
            .fold(f64::NEG_INFINITY, f64::max)
    }

    /// Minimum temperature in the array.
    pub fn min_temp(&self) -> f64 {
        self.readings.iter().cloned().fold(f64::INFINITY, f64::min)
    }

    /// Mean temperature in the array.
    pub fn mean_temp(&self) -> f64 {
        if self.readings.is_empty() {
            return 0.0;
        }
        self.readings.iter().sum::<f64>() / self.readings.len() as f64
    }

    /// Root-mean-square deviation from the mean.
    pub fn rms_deviation(&self) -> f64 {
        if self.readings.is_empty() {
            return 0.0;
        }
        let mean = self.mean_temp();
        let var = self
            .readings
            .iter()
            .map(|&t| (t - mean).powi(2))
            .sum::<f64>()
            / self.readings.len() as f64;
        var.sqrt()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// StrainGaugeSample
// ─────────────────────────────────────────────────────────────────────────────

/// A reading from a resistive strain gauge.
#[derive(Debug, Clone, PartialEq)]
pub struct StrainGaugeSample {
    /// Gauge identifier (0-based).
    pub gauge_id: u32,
    /// Timestamp in seconds.
    pub timestamp: f64,
    /// Raw bridge voltage imbalance in Volts.
    pub bridge_voltage_v: f64,
    /// Calibrated strain in µm/m (microstrain).
    pub microstrain: f64,
    /// Gauge temperature in °C (for temperature compensation).
    pub temperature_c: f64,
}

impl StrainGaugeSample {
    /// Create a new sample.
    pub fn new(
        gauge_id: u32,
        timestamp: f64,
        bridge_voltage_v: f64,
        microstrain: f64,
        temperature_c: f64,
    ) -> Self {
        Self {
            gauge_id,
            timestamp,
            bridge_voltage_v,
            microstrain,
            temperature_c,
        }
    }

    /// Stress estimate in MPa from microstrain and Young's modulus `e_gpa`.
    ///
    /// Uses Hooke's law: σ = E · ε, with ε in µm/m → m/m.
    pub fn stress_mpa(&self, e_gpa: f64) -> f64 {
        self.microstrain * 1e-6 * e_gpa * 1e3
    }

    /// Serialize to 36 bytes (4-byte gauge_id + four f64 fields).
    pub fn to_bytes(&self) -> Vec<u8> {
        let mut buf = Vec::with_capacity(36);
        buf.extend_from_slice(&self.gauge_id.to_le_bytes());
        push_f64_le(&mut buf, self.timestamp);
        push_f64_le(&mut buf, self.bridge_voltage_v);
        push_f64_le(&mut buf, self.microstrain);
        push_f64_le(&mut buf, self.temperature_c);
        buf
    }

    /// Deserialize from bytes.
    pub fn from_bytes(data: &[u8]) -> Option<Self> {
        if data.len() < 36 {
            return None;
        }
        let gauge_id = u32::from_le_bytes(data[0..4].try_into().ok()?);
        Some(Self {
            gauge_id,
            timestamp: read_f64_le(data, 4)?,
            bridge_voltage_v: read_f64_le(data, 12)?,
            microstrain: read_f64_le(data, 20)?,
            temperature_c: read_f64_le(data, 28)?,
        })
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// ForceTorqueSample — 6-DOF force/torque sensor
// ─────────────────────────────────────────────────────────────────────────────

/// A single reading from a 6-DOF force/torque sensor.
///
/// Forces are in Newtons; torques are in Newton·metres.
#[derive(Debug, Clone, PartialEq)]
pub struct ForceTorqueSample {
    /// Timestamp in seconds.
    pub timestamp: f64,
    /// Force vector `[Fx, Fy, Fz]` in N.
    pub force_n: [f64; 3],
    /// Torque vector `[Tx, Ty, Tz]` in N·m.
    pub torque_nm: [f64; 3],
}

impl ForceTorqueSample {
    /// Create a new 6-DOF sample.
    pub fn new(timestamp: f64, force_n: [f64; 3], torque_nm: [f64; 3]) -> Self {
        Self {
            timestamp,
            force_n,
            torque_nm,
        }
    }

    /// Net force magnitude in N.
    pub fn force_magnitude(&self) -> f64 {
        let [fx, fy, fz] = self.force_n;
        (fx * fx + fy * fy + fz * fz).sqrt()
    }

    /// Net torque magnitude in N·m.
    pub fn torque_magnitude(&self) -> f64 {
        let [tx, ty, tz] = self.torque_nm;
        (tx * tx + ty * ty + tz * tz).sqrt()
    }

    /// Serialize to 56 bytes (8 + 3×8 + 3×8).
    pub fn to_bytes(&self) -> Vec<u8> {
        let mut buf = Vec::with_capacity(56);
        push_f64_le(&mut buf, self.timestamp);
        for &v in &self.force_n {
            push_f64_le(&mut buf, v);
        }
        for &v in &self.torque_nm {
            push_f64_le(&mut buf, v);
        }
        buf
    }

    /// Deserialize from 56 bytes.
    pub fn from_bytes(data: &[u8]) -> Option<Self> {
        if data.len() < 56 {
            return None;
        }
        let mut o = 0usize;
        let mut next = || {
            let v = read_f64_le(data, o);
            o += 8;
            v
        };
        Some(Self {
            timestamp: next()?,
            force_n: [next()?, next()?, next()?],
            torque_nm: [next()?, next()?, next()?],
        })
    }
}

/// A log of sequential [`ForceTorqueSample`] records.
#[derive(Debug, Clone, Default)]
pub struct ForceTorqueLog {
    /// All samples, newest appended last.
    pub samples: Vec<ForceTorqueSample>,
}

impl ForceTorqueLog {
    /// Create an empty log.
    pub fn new() -> Self {
        Self {
            samples: Vec::new(),
        }
    }

    /// Append a sample.
    pub fn push(&mut self, s: ForceTorqueSample) {
        self.samples.push(s);
    }

    /// Peak force magnitude seen across all samples.
    pub fn peak_force(&self) -> f64 {
        self.samples
            .iter()
            .map(|s| s.force_magnitude())
            .fold(0.0_f64, f64::max)
    }

    /// Peak torque magnitude seen across all samples.
    pub fn peak_torque(&self) -> f64 {
        self.samples
            .iter()
            .map(|s| s.torque_magnitude())
            .fold(0.0_f64, f64::max)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// OpticalEncoderSample — incremental / absolute rotary encoder
// ─────────────────────────────────────────────────────────────────────────────

/// A reading from an optical encoder.
#[derive(Debug, Clone, PartialEq)]
pub struct OpticalEncoderSample {
    /// Encoder axis ID.
    pub axis_id: u32,
    /// Timestamp in seconds.
    pub timestamp: f64,
    /// Absolute angle in radians.
    pub angle_rad: f64,
    /// Angular velocity in rad/s (computed or measured).
    pub velocity_rad_s: f64,
    /// Encoder count (raw ticks).
    pub ticks: i64,
}

impl OpticalEncoderSample {
    /// Create a new sample.
    pub fn new(
        axis_id: u32,
        timestamp: f64,
        angle_rad: f64,
        velocity_rad_s: f64,
        ticks: i64,
    ) -> Self {
        Self {
            axis_id,
            timestamp,
            angle_rad,
            velocity_rad_s,
            ticks,
        }
    }

    /// Angle in degrees.
    pub fn angle_deg(&self) -> f64 {
        self.angle_rad.to_degrees()
    }

    /// Angular velocity in RPM.
    pub fn velocity_rpm(&self) -> f64 {
        self.velocity_rad_s * 60.0 / (2.0 * std::f64::consts::PI)
    }

    /// Serialize to 36 bytes (4 + 8 + 8 + 8 + 8).
    pub fn to_bytes(&self) -> Vec<u8> {
        let mut buf = Vec::with_capacity(36);
        buf.extend_from_slice(&self.axis_id.to_le_bytes());
        push_f64_le(&mut buf, self.timestamp);
        push_f64_le(&mut buf, self.angle_rad);
        push_f64_le(&mut buf, self.velocity_rad_s);
        buf.extend_from_slice(&self.ticks.to_le_bytes());
        buf
    }

    /// Deserialize from bytes.
    pub fn from_bytes(data: &[u8]) -> Option<Self> {
        if data.len() < 36 {
            return None;
        }
        let axis_id = u32::from_le_bytes(data[0..4].try_into().ok()?);
        let ticks = i64::from_le_bytes(data[28..36].try_into().ok()?);
        Some(Self {
            axis_id,
            timestamp: read_f64_le(data, 4)?,
            angle_rad: read_f64_le(data, 12)?,
            velocity_rad_s: read_f64_le(data, 20)?,
            ticks,
        })
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// LvdtSample — Linear Variable Differential Transformer
// ─────────────────────────────────────────────────────────────────────────────

/// A reading from an LVDT displacement sensor.
#[derive(Debug, Clone, PartialEq)]
pub struct LvdtSample {
    /// Sensor ID.
    pub sensor_id: u32,
    /// Timestamp in seconds.
    pub timestamp: f64,
    /// Displacement in metres (positive = extension).
    pub displacement_m: f64,
    /// Raw output voltage in Volts.
    pub raw_voltage_v: f64,
}

impl LvdtSample {
    /// Create a new sample.
    pub fn new(sensor_id: u32, timestamp: f64, displacement_m: f64, raw_voltage_v: f64) -> Self {
        Self {
            sensor_id,
            timestamp,
            displacement_m,
            raw_voltage_v,
        }
    }

    /// Displacement in millimetres.
    pub fn displacement_mm(&self) -> f64 {
        self.displacement_m * 1e3
    }

    /// Serialize to 28 bytes.
    pub fn to_bytes(&self) -> Vec<u8> {
        let mut buf = Vec::with_capacity(28);
        buf.extend_from_slice(&self.sensor_id.to_le_bytes());
        push_f64_le(&mut buf, self.timestamp);
        push_f64_le(&mut buf, self.displacement_m);
        push_f64_le(&mut buf, self.raw_voltage_v);
        buf
    }

    /// Deserialize from bytes.
    pub fn from_bytes(data: &[u8]) -> Option<Self> {
        if data.len() < 28 {
            return None;
        }
        let sensor_id = u32::from_le_bytes(data[0..4].try_into().ok()?);
        Some(Self {
            sensor_id,
            timestamp: read_f64_le(data, 4)?,
            displacement_m: read_f64_le(data, 12)?,
            raw_voltage_v: read_f64_le(data, 20)?,
        })
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// ThermocoupleCalibration — Seebeck-coefficient polynomial table
// ─────────────────────────────────────────────────────────────────────────────

/// Thermocouple type classification.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ThermocoupleType {
    /// Type K (chromel–alumel), −200…+1350 °C.
    TypeK,
    /// Type J (iron–constantan), −210…+1200 °C.
    TypeJ,
    /// Type T (copper–constantan), −250…+400 °C.
    TypeT,
    /// Type E (chromel–constantan), −270…+1000 °C.
    TypeE,
    /// Type N (nicrosil–nisil), −270…+1300 °C.
    TypeN,
    /// Type R (platinum-13%rhodium / platinum), 0…+1768 °C.
    TypeR,
    /// Type S (platinum-10%rhodium / platinum), 0…+1768 °C.
    TypeS,
    /// Type B (platinum-30%rhodium / platinum-6%rhodium), 0…+1820 °C.
    TypeB,
}

/// Polynomial calibration record for a thermocouple.
///
/// The EMF (Volts) to temperature (°C) relationship is expressed as
/// `T = Σ cᵢ · V^i` evaluated at ambient-compensated EMF.
#[derive(Debug, Clone)]
pub struct ThermocoupleCalibration {
    /// Thermocouple type.
    pub tc_type: ThermocoupleType,
    /// Cold-junction (reference) temperature in °C.
    pub cold_junction_c: f64,
    /// Polynomial coefficients `[c0, c1, ..., cn]` for EMF (mV) → °C.
    pub poly_coeffs: Vec<f64>,
    /// Valid EMF range `[emf_min_mv, emf_max_mv]` in mV.
    pub valid_range_mv: [f64; 2],
}

impl ThermocoupleCalibration {
    /// Create a calibration record.
    pub fn new(
        tc_type: ThermocoupleType,
        cold_junction_c: f64,
        poly_coeffs: Vec<f64>,
        valid_range_mv: [f64; 2],
    ) -> Self {
        Self {
            tc_type,
            cold_junction_c,
            poly_coeffs,
            valid_range_mv,
        }
    }

    /// Evaluate the calibration polynomial at EMF `emf_mv` in mV.
    ///
    /// Returns `None` if `emf_mv` is outside the valid range.
    pub fn convert(&self, emf_mv: f64) -> Option<f64> {
        let [lo, hi] = self.valid_range_mv;
        if emf_mv < lo || emf_mv > hi {
            return None;
        }
        let mut result = 0.0f64;
        let mut power = 1.0f64;
        for &c in &self.poly_coeffs {
            result += c * power;
            power *= emf_mv;
        }
        Some(result + self.cold_junction_c)
    }

    /// Approximate type-K conversion for a raw EMF in mV (for quick checks).
    ///
    /// Uses a simplified linear sensitivity of 41 µV/°C.
    pub fn type_k_linear(emf_mv: f64) -> f64 {
        emf_mv / 0.041
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// KalmanState — 15-state inertial navigation Kalman filter output
// ─────────────────────────────────────────────────────────────────────────────

/// The fused navigation state from a Kalman filter.
///
/// State vector: position (3), velocity (3), orientation quaternion (4),
/// accelerometer bias (3), gyroscope bias (3) → total 16 scalars.
/// The covariance diagonal (15 terms) is stored separately.
#[derive(Debug, Clone)]
pub struct KalmanState {
    /// Timestamp in seconds.
    pub timestamp: f64,
    /// Position `[x, y, z]` in metres (ECEF or local NED depending on frame).
    pub position_m: [f64; 3],
    /// Velocity `[vx, vy, vz]` in m/s.
    pub velocity_ms: [f64; 3],
    /// Orientation quaternion `[qw, qx, qy, qz]` (unit quaternion).
    pub quaternion: [f64; 4],
    /// Accelerometer bias estimate `[bax, bay, baz]` in m/s².
    pub accel_bias: [f64; 3],
    /// Gyroscope bias estimate `[bgx, bgy, bgz]` in rad/s.
    pub gyro_bias: [f64; 3],
    /// Diagonal of the 15×15 error covariance matrix.
    pub cov_diagonal: [f64; 15],
}

impl KalmanState {
    /// Create a new Kalman state with identity quaternion and zero everything else.
    pub fn identity(timestamp: f64) -> Self {
        Self {
            timestamp,
            position_m: [0.0; 3],
            velocity_ms: [0.0; 3],
            quaternion: [1.0, 0.0, 0.0, 0.0],
            accel_bias: [0.0; 3],
            gyro_bias: [0.0; 3],
            cov_diagonal: [1.0; 15],
        }
    }

    /// Check that the quaternion is (approximately) normalised.
    pub fn quaternion_is_unit(&self) -> bool {
        let [qw, qx, qy, qz] = self.quaternion;
        let norm = (qw * qw + qx * qx + qy * qy + qz * qz).sqrt();
        (norm - 1.0).abs() < 1e-6
    }

    /// Roll angle in radians extracted from the quaternion.
    pub fn roll_rad(&self) -> f64 {
        let [qw, qx, qy, qz] = self.quaternion;
        let sinr_cosp = 2.0 * (qw * qx + qy * qz);
        let cosr_cosp = 1.0 - 2.0 * (qx * qx + qy * qy);
        sinr_cosp.atan2(cosr_cosp)
    }

    /// Pitch angle in radians extracted from the quaternion.
    pub fn pitch_rad(&self) -> f64 {
        let [qw, qx, qy, qz] = self.quaternion;
        let sinp = 2.0 * (qw * qy - qz * qx);
        sinp.clamp(-1.0, 1.0).asin()
    }

    /// Yaw angle in radians extracted from the quaternion.
    pub fn yaw_rad(&self) -> f64 {
        let [qw, qx, qy, qz] = self.quaternion;
        let siny_cosp = 2.0 * (qw * qz + qx * qy);
        let cosy_cosp = 1.0 - 2.0 * (qy * qy + qz * qz);
        siny_cosp.atan2(cosy_cosp)
    }

    /// Serialize to bytes (timestamp 8B + 3+3+4+3+3 f64s + 15 cov f64s = 8+95×8 = 768B).
    pub fn to_bytes(&self) -> Vec<u8> {
        let n_floats = 1 + 3 + 3 + 4 + 3 + 3 + 15; // = 32
        let mut buf = Vec::with_capacity(n_floats * 8);
        push_f64_le(&mut buf, self.timestamp);
        for &v in &self.position_m {
            push_f64_le(&mut buf, v);
        }
        for &v in &self.velocity_ms {
            push_f64_le(&mut buf, v);
        }
        for &v in &self.quaternion {
            push_f64_le(&mut buf, v);
        }
        for &v in &self.accel_bias {
            push_f64_le(&mut buf, v);
        }
        for &v in &self.gyro_bias {
            push_f64_le(&mut buf, v);
        }
        for &v in &self.cov_diagonal {
            push_f64_le(&mut buf, v);
        }
        buf
    }

    /// Deserialize from bytes.
    pub fn from_bytes(data: &[u8]) -> Option<Self> {
        let needed = 32 * 8;
        if data.len() < needed {
            return None;
        }
        let mut o = 0usize;
        let mut next = || {
            let v = read_f64_le(data, o);
            o += 8;
            v
        };
        let timestamp = next()?;
        let position_m = [next()?, next()?, next()?];
        let velocity_ms = [next()?, next()?, next()?];
        let quaternion = [next()?, next()?, next()?, next()?];
        let accel_bias = [next()?, next()?, next()?];
        let gyro_bias = [next()?, next()?, next()?];
        let mut cov = [0.0f64; 15];
        for c in &mut cov {
            *c = next()?;
        }
        Some(Self {
            timestamp,
            position_m,
            velocity_ms,
            quaternion,
            accel_bias,
            gyro_bias,
            cov_diagonal: cov,
        })
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// CalibrationMatrix — generic N×M matrix with metadata
// ─────────────────────────────────────────────────────────────────────────────

/// A labelled calibration matrix with row/column metadata.
///
/// Can represent misalignment corrections, scale-factor matrices, cross-axis
/// coupling matrices, or any linear mapping used for sensor compensation.
#[derive(Debug, Clone)]
pub struct CalibrationMatrix {
    /// Human-readable name (e.g., `"accel_misalignment"`).
    pub name: String,
    /// Number of rows.
    pub rows: usize,
    /// Number of columns.
    pub cols: usize,
    /// Elements in row-major order, length `rows × cols`.
    pub data: Vec<f64>,
    /// Unit string for the input quantities (e.g., `"m/s²"`).
    pub input_unit: String,
    /// Unit string for the output quantities (e.g., `"m/s²"`).
    pub output_unit: String,
    /// Calibration date as a Unix timestamp in seconds.
    pub calibration_timestamp: f64,
}

impl CalibrationMatrix {
    /// Create a new identity-like calibration matrix (diagonal 1, off-diagonal 0).
    ///
    /// The matrix must be square (`rows == cols`) for the identity to make sense;
    /// for non-square matrices the main diagonal up to `min(rows,cols)` is set to 1.
    pub fn identity(
        name: impl Into<String>,
        rows: usize,
        cols: usize,
        input_unit: impl Into<String>,
        output_unit: impl Into<String>,
    ) -> Self {
        let mut data = vec![0.0f64; rows * cols];
        for i in 0..rows.min(cols) {
            data[i * cols + i] = 1.0;
        }
        Self {
            name: name.into(),
            rows,
            cols,
            data,
            input_unit: input_unit.into(),
            output_unit: output_unit.into(),
            calibration_timestamp: 0.0,
        }
    }

    /// Get element `(row, col)`.
    pub fn get(&self, row: usize, col: usize) -> Option<f64> {
        if row >= self.rows || col >= self.cols {
            return None;
        }
        Some(self.data[row * self.cols + col])
    }

    /// Set element `(row, col)`.
    pub fn set(&mut self, row: usize, col: usize, value: f64) -> bool {
        if row >= self.rows || col >= self.cols {
            return false;
        }
        self.data[row * self.cols + col] = value;
        true
    }

    /// Apply the calibration matrix to input vector `v` (length == `cols`).
    ///
    /// Returns `None` when `v.len() != cols`.
    pub fn apply(&self, v: &[f64]) -> Option<Vec<f64>> {
        if v.len() != self.cols {
            return None;
        }
        let mut out = vec![0.0f64; self.rows];
        for (r, o) in out.iter_mut().enumerate() {
            for (c, &vi) in v.iter().enumerate() {
                *o += self.data[r * self.cols + c] * vi;
            }
        }
        Some(out)
    }

    /// Frobenius norm of the matrix.
    pub fn frobenius_norm(&self) -> f64 {
        self.data.iter().map(|&v| v * v).sum::<f64>().sqrt()
    }

    /// Serialize to bytes.
    ///
    /// Layout: `[rows u64][cols u64][calib_ts f64][data… f64s]`
    /// (name and units are not persisted in this binary format).
    pub fn to_bytes(&self) -> Vec<u8> {
        let n = self.rows * self.cols;
        let mut buf = Vec::with_capacity(24 + n * 8);
        buf.extend_from_slice(&(self.rows as u64).to_le_bytes());
        buf.extend_from_slice(&(self.cols as u64).to_le_bytes());
        push_f64_le(&mut buf, self.calibration_timestamp);
        for &v in &self.data {
            push_f64_le(&mut buf, v);
        }
        buf
    }

    /// Deserialize from bytes.
    pub fn from_bytes(data: &[u8]) -> Option<Self> {
        if data.len() < 24 {
            return None;
        }
        let rows = u64::from_le_bytes(data[0..8].try_into().ok()?) as usize;
        let cols = u64::from_le_bytes(data[8..16].try_into().ok()?) as usize;
        let calib_ts = read_f64_le(data, 16)?;
        let n = rows * cols;
        if data.len() < 24 + n * 8 {
            return None;
        }
        let mut elems = Vec::with_capacity(n);
        for i in 0..n {
            elems.push(read_f64_le(data, 24 + i * 8)?);
        }
        Some(Self {
            name: String::new(),
            rows,
            cols,
            data: elems,
            input_unit: String::new(),
            output_unit: String::new(),
            calibration_timestamp: calib_ts,
        })
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// SensorFusionRecord — fused multi-sensor output
// ─────────────────────────────────────────────────────────────────────────────

/// Fused output combining IMU, pressure, and encoder data after Kalman filtering.
#[derive(Debug, Clone)]
pub struct SensorFusionRecord {
    /// Timestamp in seconds.
    pub timestamp: f64,
    /// Estimated position `[x, y, z]` in metres.
    pub position_m: [f64; 3],
    /// Estimated velocity `[vx, vy, vz]` in m/s.
    pub velocity_ms: [f64; 3],
    /// Estimated Euler angles `[roll, pitch, yaw]` in radians.
    pub euler_rad: [f64; 3],
    /// Estimated altitude in metres (from barometric fusion).
    pub altitude_m: f64,
    /// Fused temperature in °C.
    pub temperature_c: f64,
    /// Innovation (residual) magnitude for IMU update.
    pub imu_innovation: f64,
}

impl SensorFusionRecord {
    /// Create a new fusion record with all-zero fields.
    pub fn zero(timestamp: f64) -> Self {
        Self {
            timestamp,
            position_m: [0.0; 3],
            velocity_ms: [0.0; 3],
            euler_rad: [0.0; 3],
            altitude_m: 0.0,
            temperature_c: 20.0,
            imu_innovation: 0.0,
        }
    }

    /// Horizontal speed in m/s (`sqrt(vx² + vy²)`).
    pub fn horizontal_speed(&self) -> f64 {
        let [vx, vy, _] = self.velocity_ms;
        (vx * vx + vy * vy).sqrt()
    }

    /// 3-D speed in m/s.
    pub fn speed_3d(&self) -> f64 {
        let [vx, vy, vz] = self.velocity_ms;
        (vx * vx + vy * vy + vz * vz).sqrt()
    }

    /// Serialize to 104 bytes.
    pub fn to_bytes(&self) -> Vec<u8> {
        let mut buf = Vec::with_capacity(104);
        push_f64_le(&mut buf, self.timestamp);
        for &v in &self.position_m {
            push_f64_le(&mut buf, v);
        }
        for &v in &self.velocity_ms {
            push_f64_le(&mut buf, v);
        }
        for &v in &self.euler_rad {
            push_f64_le(&mut buf, v);
        }
        push_f64_le(&mut buf, self.altitude_m);
        push_f64_le(&mut buf, self.temperature_c);
        push_f64_le(&mut buf, self.imu_innovation);
        buf
    }

    /// Deserialize from 104 bytes.
    pub fn from_bytes(data: &[u8]) -> Option<Self> {
        if data.len() < 104 {
            return None;
        }
        let mut o = 0usize;
        let mut next = || {
            let v = read_f64_le(data, o);
            o += 8;
            v
        };
        Some(Self {
            timestamp: next()?,
            position_m: [next()?, next()?, next()?],
            velocity_ms: [next()?, next()?, next()?],
            euler_rad: [next()?, next()?, next()?],
            altitude_m: next()?,
            temperature_c: next()?,
            imu_innovation: next()?,
        })
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// SensorDataHeader — binary file header
// ─────────────────────────────────────────────────────────────────────────────

/// Magic bytes used to identify OxiPhysics sensor data binary files.
pub const SENSOR_DATA_MAGIC: [u8; 8] = *b"OXISENS\0";

/// Version of the sensor data binary format.
pub const SENSOR_DATA_VERSION: u32 = 1;

/// Binary file header for sensor data archives.
#[derive(Debug, Clone)]
pub struct SensorDataHeader {
    /// Sensor node name.
    pub node_name: String,
    /// Sample rate in Hz.
    pub sample_rate_hz: f64,
    /// Number of records in the file.
    pub record_count: u64,
    /// File creation timestamp (Unix time, seconds).
    pub created_at: f64,
}

impl SensorDataHeader {
    /// Create a new header.
    pub fn new(
        node_name: impl Into<String>,
        sample_rate_hz: f64,
        record_count: u64,
        created_at: f64,
    ) -> Self {
        Self {
            node_name: node_name.into(),
            sample_rate_hz,
            record_count,
            created_at,
        }
    }

    /// Sample interval in seconds (reciprocal of sample rate).
    pub fn sample_interval_s(&self) -> f64 {
        if self.sample_rate_hz == 0.0 {
            f64::INFINITY
        } else {
            1.0 / self.sample_rate_hz
        }
    }

    /// Duration of the recording in seconds.
    pub fn duration_s(&self) -> f64 {
        self.record_count as f64 * self.sample_interval_s()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// SensorDataStats — running statistics helper
// ─────────────────────────────────────────────────────────────────────────────

/// Welford online-variance statistics accumulator for one scalar channel.
///
/// Computes mean and variance incrementally without storing all samples.
#[derive(Debug, Clone)]
pub struct ChannelStats {
    /// Number of samples observed.
    pub count: u64,
    /// Running mean.
    pub mean: f64,
    /// Running M2 (sum of squared deviations) for Welford's algorithm.
    pub m2: f64,
    /// Minimum value seen.
    pub min: f64,
    /// Maximum value seen.
    pub max: f64,
}

impl Default for ChannelStats {
    fn default() -> Self {
        Self {
            count: 0,
            mean: 0.0,
            m2: 0.0,
            min: f64::INFINITY,
            max: f64::NEG_INFINITY,
        }
    }
}

impl ChannelStats {
    /// Create an empty accumulator.
    pub fn new() -> Self {
        Self::default()
    }

    /// Update with a new sample value.
    pub fn update(&mut self, value: f64) {
        self.count += 1;
        if value < self.min {
            self.min = value;
        }
        if value > self.max {
            self.max = value;
        }
        let delta = value - self.mean;
        self.mean += delta / self.count as f64;
        let delta2 = value - self.mean;
        self.m2 += delta * delta2;
    }

    /// Sample variance (Bessel-corrected).
    pub fn variance(&self) -> f64 {
        if self.count < 2 {
            0.0
        } else {
            self.m2 / (self.count - 1) as f64
        }
    }

    /// Sample standard deviation.
    pub fn std_dev(&self) -> f64 {
        self.variance().sqrt()
    }

    /// Range `max - min`.
    pub fn range(&self) -> f64 {
        if self.count == 0 {
            0.0
        } else {
            self.max - self.min
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// MultiChannelSensorLog — generic timestamped multi-channel log
// ─────────────────────────────────────────────────────────────────────────────

/// A generic multi-channel sensor log, suitable for storing any mix of
/// floating-point channel data alongside a common timestamp.
#[derive(Debug, Clone)]
pub struct MultiChannelSensorLog {
    /// Channel names, one per column.
    pub channel_names: Vec<String>,
    /// Rows, each containing a timestamp followed by one value per channel.
    pub rows: Vec<Vec<f64>>,
}

impl MultiChannelSensorLog {
    /// Create an empty log with the given channel names.
    pub fn new(channel_names: Vec<String>) -> Self {
        Self {
            channel_names,
            rows: Vec::new(),
        }
    }

    /// Number of channels (excluding the timestamp column).
    pub fn num_channels(&self) -> usize {
        self.channel_names.len()
    }

    /// Append a row.  `timestamp` is prepended; `values` must have exactly
    /// `num_channels` elements.  Returns `false` when the length is wrong.
    pub fn push_row(&mut self, timestamp: f64, values: &[f64]) -> bool {
        if values.len() != self.num_channels() {
            return false;
        }
        let mut row = Vec::with_capacity(1 + values.len());
        row.push(timestamp);
        row.extend_from_slice(values);
        self.rows.push(row);
        true
    }

    /// Extract all values for channel `idx` (0-based, excludes timestamp).
    pub fn channel_values(&self, idx: usize) -> Vec<f64> {
        if idx >= self.num_channels() {
            return vec![];
        }
        self.rows.iter().map(|r| r[idx + 1]).collect()
    }

    /// Compute per-channel statistics.
    pub fn channel_stats(&self, idx: usize) -> ChannelStats {
        let mut s = ChannelStats::new();
        for v in self.channel_values(idx) {
            s.update(v);
        }
        s
    }

    /// Export to CSV lines (header + data rows).
    pub fn to_csv_lines(&self) -> Vec<String> {
        let header = {
            let mut h = String::from("timestamp");
            for name in &self.channel_names {
                h.push(',');
                h.push_str(name);
            }
            h
        };
        let mut lines = vec![header];
        for row in &self.rows {
            let parts: Vec<String> = row.iter().map(|v| format!("{:.6}", v)).collect();
            lines.push(parts.join(","));
        }
        lines
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Unit tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // ── ImuSample ────────────────────────────────────────────────────────────

    #[test]
    fn test_imu_sample_round_trip() {
        let s = ImuSample::new(
            1.23,
            [1.0, -2.0, 9.8],
            [0.01, -0.02, 0.003],
            [24.1, -5.0, 42.0],
        );
        let bytes = s.to_bytes();
        assert_eq!(bytes.len(), 80);
        let s2 = ImuSample::from_bytes(&bytes).unwrap();
        assert!((s2.timestamp - s.timestamp).abs() < 1e-12);
        for k in 0..3 {
            assert!((s2.accel[k] - s.accel[k]).abs() < 1e-12);
            assert!((s2.gyro[k] - s.gyro[k]).abs() < 1e-12);
            assert!((s2.mag[k] - s.mag[k]).abs() < 1e-12);
        }
    }

    #[test]
    fn test_imu_sample_accel_magnitude() {
        let s = ImuSample::new(0.0, [3.0, 4.0, 0.0], [0.0; 3], [0.0; 3]);
        assert!((s.accel_magnitude() - 5.0).abs() < 1e-12);
    }

    #[test]
    fn test_imu_sample_gyro_magnitude() {
        let s = ImuSample::new(0.0, [0.0; 3], [0.0, 0.0, 1.0], [0.0; 3]);
        assert!((s.gyro_magnitude() - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_imu_sample_from_bytes_truncated() {
        let short = [0u8; 10];
        assert!(ImuSample::from_bytes(&short).is_none());
    }

    // ── ImuStream ────────────────────────────────────────────────────────────

    #[test]
    fn test_imu_stream_round_trip() {
        let mut stream = ImuStream::new();
        for i in 0..5 {
            stream.push(ImuSample::new(
                i as f64 * 0.01,
                [0.0, 0.0, 9.81],
                [0.0; 3],
                [0.0; 3],
            ));
        }
        let bytes = stream.to_bytes();
        let s2 = ImuStream::from_bytes(&bytes).unwrap();
        assert_eq!(s2.len(), 5);
        assert!((s2.samples[0].timestamp).abs() < 1e-12);
    }

    #[test]
    fn test_imu_stream_mean_accel() {
        let mut stream = ImuStream::new();
        stream.push(ImuSample::new(0.0, [2.0, 0.0, 0.0], [0.0; 3], [0.0; 3]));
        stream.push(ImuSample::new(0.01, [4.0, 0.0, 0.0], [0.0; 3], [0.0; 3]));
        let m = stream.mean_accel();
        assert!((m[0] - 3.0).abs() < 1e-12);
    }

    #[test]
    fn test_imu_stream_empty_mean() {
        let stream = ImuStream::new();
        assert_eq!(stream.mean_accel(), [0.0; 3]);
        assert_eq!(stream.mean_gyro(), [0.0; 3]);
    }

    // ── PressureSample ───────────────────────────────────────────────────────

    #[test]
    fn test_pressure_sample_round_trip() {
        let s = PressureSample::new(5.0, 101325.0, 22.5);
        let bytes = s.to_bytes();
        let s2 = PressureSample::from_bytes(&bytes).unwrap();
        assert!((s2.pressure_pa - 101325.0).abs() < 1e-6);
        assert!((s2.temperature_c - 22.5).abs() < 1e-9);
    }

    #[test]
    fn test_pressure_altitude_sea_level() {
        let s = PressureSample::new(0.0, 101325.0, 15.0);
        let alt = s.altitude_m(101325.0, 0.0);
        assert!(
            alt.abs() < 1.0,
            "altitude at sea-level pressure should be ~0 m, got {}",
            alt
        );
    }

    #[test]
    fn test_pressure_altitude_decreases_with_pressure() {
        let s_low = PressureSample::new(0.0, 101325.0, 15.0);
        let s_high = PressureSample::new(0.0, 89875.0, 10.0);
        let alt_low = s_low.altitude_m(101325.0, 0.0);
        let alt_high = s_high.altitude_m(101325.0, 0.0);
        assert!(alt_high > alt_low, "lower pressure → higher altitude");
    }

    #[test]
    fn test_pressure_time_series_statistics() {
        let mut ts = PressureTimeSeries::new();
        for i in 0..10 {
            ts.push(PressureSample::new(
                i as f64,
                100000.0 + i as f64 * 100.0,
                20.0,
            ));
        }
        assert!((ts.min_pressure() - 100000.0).abs() < 1e-6);
        assert!((ts.max_pressure() - 100900.0).abs() < 1e-6);
        assert!((ts.mean_pressure() - 100450.0).abs() < 1e-6);
    }

    // ── TempProbeReading ─────────────────────────────────────────────────────

    #[test]
    fn test_temp_probe_round_trip() {
        let r = TempProbeReading::new(3, 1.5, 37.5, [1.0, 2.0, 3.0]);
        let bytes = r.to_bytes();
        assert_eq!(bytes.len(), 44);
        let r2 = TempProbeReading::from_bytes(&bytes).unwrap();
        assert_eq!(r2.probe_id, 3);
        assert!((r2.temperature_c - 37.5).abs() < 1e-9);
    }

    #[test]
    fn test_temp_array_snapshot_stats() {
        let snap = TempArraySnapshot::new(0.0, vec![10.0, 20.0, 30.0, 40.0]);
        assert!((snap.min_temp() - 10.0).abs() < 1e-9);
        assert!((snap.max_temp() - 40.0).abs() < 1e-9);
        assert!((snap.mean_temp() - 25.0).abs() < 1e-9);
    }

    #[test]
    fn test_temp_array_rms_deviation() {
        // All equal → rms_deviation should be 0
        let snap = TempArraySnapshot::new(0.0, vec![20.0, 20.0, 20.0]);
        assert!(snap.rms_deviation() < 1e-9);
    }

    // ── StrainGaugeSample ────────────────────────────────────────────────────

    #[test]
    fn test_strain_gauge_round_trip() {
        let s = StrainGaugeSample::new(0, 2.0, 0.001, 500.0, 23.5);
        let bytes = s.to_bytes();
        let s2 = StrainGaugeSample::from_bytes(&bytes).unwrap();
        assert!((s2.microstrain - 500.0).abs() < 1e-9);
    }

    #[test]
    fn test_strain_gauge_stress() {
        // 1000 µm/m at 200 GPa → 200 MPa
        let s = StrainGaugeSample::new(0, 0.0, 0.0, 1000.0, 20.0);
        let stress = s.stress_mpa(200.0);
        assert!((stress - 200.0).abs() < 1e-6);
    }

    // ── ForceTorqueSample ────────────────────────────────────────────────────

    #[test]
    fn test_force_torque_round_trip() {
        let s = ForceTorqueSample::new(3.125, [10.0, 0.0, -5.0], [0.0, 1.0, 0.0]);
        let bytes = s.to_bytes();
        assert_eq!(bytes.len(), 56);
        let s2 = ForceTorqueSample::from_bytes(&bytes).unwrap();
        assert!((s2.force_n[0] - 10.0).abs() < 1e-9);
        assert!((s2.torque_nm[1] - 1.0).abs() < 1e-9);
    }

    #[test]
    fn test_force_torque_magnitudes() {
        let s = ForceTorqueSample::new(0.0, [3.0, 4.0, 0.0], [0.0, 0.0, 5.0]);
        assert!((s.force_magnitude() - 5.0).abs() < 1e-12);
        assert!((s.torque_magnitude() - 5.0).abs() < 1e-12);
    }

    #[test]
    fn test_force_torque_log_peaks() {
        let mut log = ForceTorqueLog::new();
        log.push(ForceTorqueSample::new(
            0.0,
            [1.0, 0.0, 0.0],
            [0.0, 0.0, 0.5],
        ));
        log.push(ForceTorqueSample::new(
            1.0,
            [10.0, 0.0, 0.0],
            [0.0, 0.0, 2.0],
        ));
        assert!((log.peak_force() - 10.0).abs() < 1e-9);
        assert!((log.peak_torque() - 2.0).abs() < 1e-9);
    }

    // ── OpticalEncoderSample ─────────────────────────────────────────────────

    #[test]
    fn test_encoder_round_trip() {
        let s = OpticalEncoderSample::new(0, 0.5, std::f64::consts::PI, 1.0, 2048);
        let bytes = s.to_bytes();
        assert_eq!(bytes.len(), 36);
        let s2 = OpticalEncoderSample::from_bytes(&bytes).unwrap();
        assert!((s2.angle_rad - std::f64::consts::PI).abs() < 1e-9);
        assert_eq!(s2.ticks, 2048);
    }

    #[test]
    fn test_encoder_angle_conversion() {
        let s = OpticalEncoderSample::new(0, 0.0, std::f64::consts::PI, 0.0, 0);
        assert!((s.angle_deg() - 180.0).abs() < 1e-9);
    }

    #[test]
    fn test_encoder_rpm() {
        // 2π rad/s → 60 RPM
        let s = OpticalEncoderSample::new(0, 0.0, 0.0, 2.0 * std::f64::consts::PI, 0);
        assert!((s.velocity_rpm() - 60.0).abs() < 1e-9);
    }

    // ── LvdtSample ────────────────────────────────────────────────────────────

    #[test]
    fn test_lvdt_round_trip() {
        let s = LvdtSample::new(1, 0.1, 0.025, 2.5);
        let bytes = s.to_bytes();
        assert_eq!(bytes.len(), 28);
        let s2 = LvdtSample::from_bytes(&bytes).unwrap();
        assert!((s2.displacement_m - 0.025).abs() < 1e-12);
    }

    #[test]
    fn test_lvdt_mm_conversion() {
        let s = LvdtSample::new(0, 0.0, 0.005, 0.5);
        assert!((s.displacement_mm() - 5.0).abs() < 1e-9);
    }

    // ── ThermocoupleCalibration ───────────────────────────────────────────────

    #[test]
    fn test_thermocouple_linear_poly() {
        // Simple linear calibration: T = 0 + 25·V
        let cal = ThermocoupleCalibration::new(
            ThermocoupleType::TypeK,
            0.0,
            vec![0.0, 25.0],
            [0.0, 10.0],
        );
        let t = cal.convert(4.0).unwrap();
        assert!((t - 100.0).abs() < 1e-9);
    }

    #[test]
    fn test_thermocouple_out_of_range() {
        let cal = ThermocoupleCalibration::new(
            ThermocoupleType::TypeK,
            0.0,
            vec![0.0, 25.0],
            [0.0, 10.0],
        );
        assert!(cal.convert(-1.0).is_none());
        assert!(cal.convert(11.0).is_none());
    }

    #[test]
    fn test_thermocouple_type_k_linear() {
        // At 1 mV → ~24.4 °C (1/0.041)
        let t = ThermocoupleCalibration::type_k_linear(1.0);
        assert!((t - 1.0 / 0.041).abs() < 1e-6);
    }

    #[test]
    fn test_thermocouple_cold_junction_offset() {
        // constant polynomial (c0 = 5, c1 = 0) + cold_junction 20 → T = 25
        let cal =
            ThermocoupleCalibration::new(ThermocoupleType::TypeJ, 20.0, vec![5.0], [0.0, 100.0]);
        let t = cal.convert(50.0).unwrap();
        assert!((t - 25.0).abs() < 1e-9);
    }

    // ── KalmanState ───────────────────────────────────────────────────────────

    #[test]
    fn test_kalman_state_round_trip() {
        let mut k = KalmanState::identity(100.0);
        k.position_m = [1.0, 2.0, 3.0];
        k.velocity_ms = [0.5, -0.1, 0.0];
        let bytes = k.to_bytes();
        let k2 = KalmanState::from_bytes(&bytes).unwrap();
        assert!((k2.timestamp - 100.0).abs() < 1e-9);
        assert!((k2.position_m[0] - 1.0).abs() < 1e-9);
    }

    #[test]
    fn test_kalman_identity_quaternion() {
        let k = KalmanState::identity(0.0);
        assert!(k.quaternion_is_unit());
    }

    #[test]
    fn test_kalman_euler_zero_at_identity() {
        let k = KalmanState::identity(0.0);
        assert!(k.roll_rad().abs() < 1e-9);
        assert!(k.pitch_rad().abs() < 1e-9);
        assert!(k.yaw_rad().abs() < 1e-9);
    }

    // ── CalibrationMatrix ─────────────────────────────────────────────────────

    #[test]
    fn test_calibration_matrix_identity_apply() {
        let m = CalibrationMatrix::identity("test", 3, 3, "m/s²", "m/s²");
        let v = [1.0, 2.0, 3.0];
        let out = m.apply(&v).unwrap();
        for k in 0..3 {
            assert!((out[k] - v[k]).abs() < 1e-12);
        }
    }

    #[test]
    fn test_calibration_matrix_round_trip() {
        let mut m = CalibrationMatrix::identity("cm", 2, 2, "in", "out");
        m.set(0, 1, 0.5);
        let bytes = m.to_bytes();
        let m2 = CalibrationMatrix::from_bytes(&bytes).unwrap();
        assert_eq!(m2.rows, 2);
        assert_eq!(m2.cols, 2);
        assert!((m2.get(0, 1).unwrap() - 0.5).abs() < 1e-12);
    }

    #[test]
    fn test_calibration_matrix_frobenius_identity() {
        let m = CalibrationMatrix::identity("id", 3, 3, "a", "b");
        // Frobenius norm of I₃ = sqrt(3)
        assert!((m.frobenius_norm() - 3.0f64.sqrt()).abs() < 1e-9);
    }

    #[test]
    fn test_calibration_matrix_wrong_vector_size() {
        let m = CalibrationMatrix::identity("id", 3, 3, "a", "b");
        assert!(m.apply(&[1.0, 2.0]).is_none());
    }

    // ── SensorFusionRecord ────────────────────────────────────────────────────

    #[test]
    fn test_fusion_record_round_trip() {
        let mut r = SensorFusionRecord::zero(9.9);
        r.velocity_ms = [3.0, 4.0, 0.0];
        r.altitude_m = 150.0;
        let bytes = r.to_bytes();
        assert_eq!(bytes.len(), 104);
        let r2 = SensorFusionRecord::from_bytes(&bytes).unwrap();
        assert!((r2.altitude_m - 150.0).abs() < 1e-9);
        assert!((r2.horizontal_speed() - 5.0).abs() < 1e-9);
    }

    #[test]
    fn test_fusion_record_speed_3d() {
        let mut r = SensorFusionRecord::zero(0.0);
        r.velocity_ms = [1.0, 2.0, 2.0];
        assert!((r.speed_3d() - 3.0).abs() < 1e-9);
    }

    // ── ChannelStats ──────────────────────────────────────────────────────────

    #[test]
    fn test_channel_stats_basic() {
        // Three samples with sample std-dev exactly 2.0: {0, 2, 4} → mean=2, m2=8, var=8/2=4
        let mut s = ChannelStats::new();
        for v in [0.0_f64, 2.0, 4.0] {
            s.update(v);
        }
        assert!((s.mean - 2.0).abs() < 1e-9);
        assert!((s.std_dev() - 2.0).abs() < 1e-9);
    }

    #[test]
    fn test_channel_stats_single_sample() {
        let mut s = ChannelStats::new();
        s.update(42.0);
        assert!((s.mean - 42.0).abs() < 1e-12);
        assert!(s.variance() < 1e-12);
        assert!((s.range()).abs() < 1e-12);
    }

    #[test]
    fn test_channel_stats_empty_variance() {
        let s = ChannelStats::new();
        assert_eq!(s.variance(), 0.0);
        assert_eq!(s.range(), 0.0);
    }

    // ── MultiChannelSensorLog ─────────────────────────────────────────────────

    #[test]
    fn test_multi_channel_push_and_retrieve() {
        let mut log = MultiChannelSensorLog::new(vec!["ax".into(), "ay".into()]);
        assert!(log.push_row(0.0, &[1.0, 2.0]));
        assert!(log.push_row(0.01, &[3.0, 4.0]));
        let ax_vals = log.channel_values(0);
        assert_eq!(ax_vals, vec![1.0, 3.0]);
    }

    #[test]
    fn test_multi_channel_wrong_row_length() {
        let mut log = MultiChannelSensorLog::new(vec!["x".into(), "y".into()]);
        assert!(!log.push_row(0.0, &[1.0])); // too few values
    }

    #[test]
    fn test_multi_channel_csv_header() {
        let log = MultiChannelSensorLog::new(vec!["p".into(), "q".into()]);
        let lines = log.to_csv_lines();
        assert_eq!(lines[0], "timestamp,p,q");
    }

    #[test]
    fn test_multi_channel_stats() {
        let mut log = MultiChannelSensorLog::new(vec!["v".into()]);
        for i in 0..5 {
            log.push_row(i as f64, &[i as f64 * 2.0]);
        }
        let stats = log.channel_stats(0);
        assert!((stats.mean - 4.0).abs() < 1e-9);
        assert!((stats.min).abs() < 1e-9);
        assert!((stats.max - 8.0).abs() < 1e-9);
    }

    // ── SensorDataHeader ──────────────────────────────────────────────────────

    #[test]
    fn test_sensor_data_header_interval() {
        let hdr = SensorDataHeader::new("node1", 100.0, 1000, 0.0);
        assert!((hdr.sample_interval_s() - 0.01).abs() < 1e-12);
        assert!((hdr.duration_s() - 10.0).abs() < 1e-9);
    }

    #[test]
    fn test_sensor_data_header_zero_rate() {
        let hdr = SensorDataHeader::new("n", 0.0, 100, 0.0);
        assert!(hdr.sample_interval_s().is_infinite());
    }
}
