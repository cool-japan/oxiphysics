// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Sensor and telemetry I/O.
//!
//! Provides writers and readers for IMU, GPS, lidar, force gauge data,
//! MAVLink-inspired binary protocol, sensor fusion, and data logging.

#![allow(dead_code)]
#![allow(clippy::too_many_arguments)]

use std::collections::VecDeque;
use std::io::{BufWriter, Write};

/// IMU data record: timestamp, accelerometer, gyroscope, magnetometer.
#[derive(Debug, Clone, PartialEq)]
pub struct ImuDataRecord {
    /// Timestamp in seconds.
    pub timestamp: f64,
    /// Accelerometer \[ax, ay, az\] in m/s².
    pub accel: [f64; 3],
    /// Gyroscope \[gx, gy, gz\] in rad/s.
    pub gyro: [f64; 3],
    /// Magnetometer \[mx, my, mz\] in µT.
    pub mag: [f64; 3],
}

impl ImuDataRecord {
    /// Create a new IMU record.
    pub fn new(timestamp: f64, accel: [f64; 3], gyro: [f64; 3], mag: [f64; 3]) -> Self {
        Self {
            timestamp,
            accel,
            gyro,
            mag,
        }
    }

    /// Serialize to bytes (little-endian f64).
    pub fn to_bytes(&self) -> Vec<u8> {
        let mut buf = Vec::with_capacity(80);
        buf.extend_from_slice(&self.timestamp.to_le_bytes());
        for &v in &self.accel {
            buf.extend_from_slice(&v.to_le_bytes());
        }
        for &v in &self.gyro {
            buf.extend_from_slice(&v.to_le_bytes());
        }
        for &v in &self.mag {
            buf.extend_from_slice(&v.to_le_bytes());
        }
        buf
    }

    /// Deserialize from bytes.
    pub fn from_bytes(data: &[u8]) -> Option<Self> {
        if data.len() < 80 {
            return None;
        }
        let mut off = 0usize;
        let read_f64 = |d: &[u8], o: &mut usize| -> f64 {
            let v = f64::from_le_bytes(d[*o..*o + 8].try_into().unwrap_or([0u8; 8]));
            *o += 8;
            v
        };
        let timestamp = read_f64(data, &mut off);
        let accel = [
            read_f64(data, &mut off),
            read_f64(data, &mut off),
            read_f64(data, &mut off),
        ];
        let gyro = [
            read_f64(data, &mut off),
            read_f64(data, &mut off),
            read_f64(data, &mut off),
        ];
        let mag = [
            read_f64(data, &mut off),
            read_f64(data, &mut off),
            read_f64(data, &mut off),
        ];
        Some(Self {
            timestamp,
            accel,
            gyro,
            mag,
        })
    }
}

/// GPS data record: timestamp, position, velocity, HDOP.
#[derive(Debug, Clone, PartialEq)]
pub struct GpsDataRecord {
    /// Timestamp in seconds.
    pub timestamp: f64,
    /// Latitude in degrees.
    pub lat: f64,
    /// Longitude in degrees.
    pub lon: f64,
    /// Altitude in meters.
    pub alt: f64,
    /// Velocity \[vn, ve, vd\] in m/s (NED frame).
    pub velocity: [f64; 3],
    /// Horizontal dilution of precision.
    pub hdop: f64,
}

impl GpsDataRecord {
    /// Create a new GPS record.
    pub fn new(
        timestamp: f64,
        lat: f64,
        lon: f64,
        alt: f64,
        velocity: [f64; 3],
        hdop: f64,
    ) -> Self {
        Self {
            timestamp,
            lat,
            lon,
            alt,
            velocity,
            hdop,
        }
    }

    /// Serialize to bytes.
    pub fn to_bytes(&self) -> Vec<u8> {
        let mut buf = Vec::with_capacity(64);
        for v in [
            self.timestamp,
            self.lat,
            self.lon,
            self.alt,
            self.velocity[0],
            self.velocity[1],
            self.velocity[2],
            self.hdop,
        ] {
            buf.extend_from_slice(&v.to_le_bytes());
        }
        buf
    }

    /// Ground speed in m/s.
    pub fn ground_speed(&self) -> f64 {
        (self.velocity[0] * self.velocity[0] + self.velocity[1] * self.velocity[1]).sqrt()
    }
}

/// LiDAR frame: timestamp, point cloud, intensities.
#[derive(Debug, Clone)]
pub struct LidarFrame {
    /// Timestamp in seconds.
    pub timestamp: f64,
    /// Number of points.
    pub n_points: usize,
    /// 3D positions \[x, y, z\] in meters (f32 for compactness).
    pub positions: Vec<[f32; 3]>,
    /// Intensity values (0.0..1.0).
    pub intensities: Vec<f32>,
}

impl LidarFrame {
    /// Create a new empty LiDAR frame.
    pub fn new(timestamp: f64) -> Self {
        Self {
            timestamp,
            n_points: 0,
            positions: Vec::new(),
            intensities: Vec::new(),
        }
    }

    /// Add a point.
    pub fn add_point(&mut self, pos: [f32; 3], intensity: f32) {
        self.positions.push(pos);
        self.intensities.push(intensity);
        self.n_points += 1;
    }

    /// Serialize to bytes.
    pub fn to_bytes(&self) -> Vec<u8> {
        let mut buf = Vec::new();
        buf.extend_from_slice(&self.timestamp.to_le_bytes());
        buf.extend_from_slice(&(self.n_points as u32).to_le_bytes());
        for p in &self.positions {
            for &v in p {
                buf.extend_from_slice(&v.to_le_bytes());
            }
        }
        for &intensity in &self.intensities {
            buf.extend_from_slice(&intensity.to_le_bytes());
        }
        buf
    }

    /// Point cloud bounding box \[\[min_x, max_x\\], \[min_y, max_y\], \[min_z, max_z\]].
    pub fn bounding_box(&self) -> [[f32; 2]; 3] {
        let mut bb = [[f32::MAX, f32::MIN]; 3];
        for p in &self.positions {
            for axis in 0..3 {
                bb[axis][0] = bb[axis][0].min(p[axis]);
                bb[axis][1] = bb[axis][1].max(p[axis]);
            }
        }
        bb
    }
}

/// Force gauge record: timestamp, force, torque, temperature.
#[derive(Debug, Clone)]
pub struct ForceGaugeRecord {
    /// Timestamp in seconds.
    pub timestamp: f64,
    /// Force \[fx, fy, fz\] in N.
    pub force: [f64; 3],
    /// Torque \[tx, ty, tz\] in N·m.
    pub torque: [f64; 3],
    /// Temperature in Celsius.
    pub temperature: f64,
}

impl ForceGaugeRecord {
    /// Create a new force gauge record.
    pub fn new(timestamp: f64, force: [f64; 3], torque: [f64; 3], temperature: f64) -> Self {
        Self {
            timestamp,
            force,
            torque,
            temperature,
        }
    }

    /// Force magnitude in N.
    pub fn force_magnitude(&self) -> f64 {
        let f = &self.force;
        (f[0] * f[0] + f[1] * f[1] + f[2] * f[2]).sqrt()
    }

    /// Torque magnitude in N·m.
    pub fn torque_magnitude(&self) -> f64 {
        let t = &self.torque;
        (t[0] * t[0] + t[1] * t[1] + t[2] * t[2]).sqrt()
    }

    /// Serialize to bytes.
    pub fn to_bytes(&self) -> Vec<u8> {
        let mut buf = Vec::with_capacity(56);
        for v in [
            self.timestamp,
            self.force[0],
            self.force[1],
            self.force[2],
            self.torque[0],
            self.torque[1],
            self.torque[2],
            self.temperature,
        ] {
            buf.extend_from_slice(&v.to_le_bytes());
        }
        buf
    }
}

/// MAVLink-inspired telemetry protocol.
///
/// Binary protocol with magic byte, version, message_id, payload length, and CRC.
#[derive(Debug, Clone)]
pub struct TelemetryProtocol {
    /// Magic byte (0xFE for MAVLink v1-style).
    pub magic: u8,
    /// Protocol version.
    pub version: u8,
    /// System ID.
    pub sys_id: u8,
    /// Component ID.
    pub comp_id: u8,
}

/// Telemetry message.
#[derive(Debug, Clone)]
pub struct TelemetryMessage {
    /// Message ID.
    pub message_id: u16,
    /// Sequence number.
    pub seq: u8,
    /// Payload bytes.
    pub payload: Vec<u8>,
}

impl TelemetryProtocol {
    /// Create a new telemetry protocol instance.
    pub fn new(sys_id: u8, comp_id: u8) -> Self {
        Self {
            magic: 0xFE,
            version: 1,
            sys_id,
            comp_id,
        }
    }

    /// Encode a message to bytes.
    pub fn encode(&self, msg: &TelemetryMessage) -> Vec<u8> {
        let payload_len = msg.payload.len() as u8;
        let mut buf = vec![self.magic, payload_len, msg.seq, self.sys_id, self.comp_id];
        buf.extend_from_slice(&msg.message_id.to_le_bytes());
        buf.extend_from_slice(&msg.payload);
        let crc = self.crc16(&buf[1..]);
        buf.extend_from_slice(&crc.to_le_bytes());
        buf
    }

    /// Decode a message from bytes.
    pub fn decode(&self, data: &[u8]) -> Option<TelemetryMessage> {
        if data.len() < 8 {
            return None;
        }
        if data[0] != self.magic {
            return None;
        }
        let payload_len = data[1] as usize;
        if data.len() < 7 + payload_len + 2 {
            return None;
        }
        let seq = data[2];
        let _sys_id = data[3];
        let _comp_id = data[4];
        let message_id = u16::from_le_bytes([data[5], data[6]]);
        let payload = data[7..7 + payload_len].to_vec();
        Some(TelemetryMessage {
            message_id,
            seq,
            payload,
        })
    }

    /// Simple CRC-16.
    fn crc16(&self, data: &[u8]) -> u16 {
        let mut crc = 0xFFFFu16;
        for &byte in data {
            let tmp = byte as u16 ^ (crc & 0xFF);
            let tmp = tmp ^ (tmp << 4);
            crc = (crc >> 8) ^ (tmp << 8) ^ (tmp << 3) ^ (tmp >> 4);
        }
        crc
    }
}

/// Sensor fusion using complementary and Madgwick AHRS filters.
///
/// Estimates orientation from IMU data.
#[derive(Debug, Clone)]
pub struct SensorFusion {
    /// Quaternion \[w, x, y, z\].
    pub quaternion: [f64; 4],
    /// Complementary filter weight alpha (for gyro vs accel).
    pub alpha: f64,
    /// Madgwick beta parameter.
    pub beta: f64,
    /// Euler angles \[roll, pitch, yaw\] in degrees.
    pub euler: [f64; 3],
}

impl SensorFusion {
    /// Create a new sensor fusion instance.
    pub fn new(alpha: f64, beta: f64) -> Self {
        Self {
            quaternion: [1.0, 0.0, 0.0, 0.0],
            alpha,
            beta,
            euler: [0.0; 3],
        }
    }

    /// Update with complementary filter.
    pub fn update_complementary(&mut self, accel: [f64; 3], gyro: [f64; 3], dt: f64) {
        // Integrate gyro
        let [w, x, y, z] = self.quaternion;
        let [gx, gy, gz] = gyro;
        let dw = 0.5 * (-x * gx - y * gy - z * gz);
        let dx = 0.5 * (w * gx + y * gz - z * gy);
        let dy = 0.5 * (w * gy - x * gz + z * gx);
        let dz = 0.5 * (w * gz + x * gy - y * gx);
        let qw = w + dw * dt;
        let qx = x + dx * dt;
        let qy = y + dy * dt;
        let qz = z + dz * dt;
        let mag = (qw * qw + qx * qx + qy * qy + qz * qz).sqrt().max(1e-10);
        let q_gyro = [qw / mag, qx / mag, qy / mag, qz / mag];

        // Accelerometer-derived roll/pitch
        let ax = accel[0];
        let ay = accel[1];
        let az = accel[2];
        let norm_a = (ax * ax + ay * ay + az * az).sqrt().max(1e-10);
        let roll = (ay / norm_a).atan2((az / norm_a).max(1e-10));
        let pitch = -(ax / norm_a)
            .asin()
            .clamp(-std::f64::consts::FRAC_PI_2, std::f64::consts::FRAC_PI_2);
        let q_accel = euler_to_quat(roll, pitch, 0.0);

        // Complementary blend
        let a = self.alpha;
        let q_blend = quat_slerp(q_accel, q_gyro, a);
        self.quaternion = q_blend;
        self.euler = quat_to_euler(self.quaternion);
    }

    /// Update with Madgwick AHRS filter.
    pub fn update_madgwick(&mut self, accel: [f64; 3], gyro: [f64; 3], dt: f64) {
        let [mut qw, mut qx, mut qy, mut qz] = self.quaternion;
        let [gx, gy, gz] = gyro;
        let [ax, ay, az] = accel;

        let norm_a = (ax * ax + ay * ay + az * az).sqrt().max(1e-10);
        let (ax, ay, az) = (ax / norm_a, ay / norm_a, az / norm_a);

        // Gradient step
        let f1 = 2.0 * (qx * qz - qw * qy) - ax;
        let f2 = 2.0 * (qw * qx + qy * qz) - ay;
        let f3 = 1.0 - 2.0 * (qx * qx + qy * qy) - az;
        let j11 = 2.0 * qy;
        let j12 = 2.0 * qz;
        let j13 = 2.0 * qw;
        let j14 = 2.0 * qx;
        let j21 = 2.0 * qx;
        let j22 = 2.0 * qw;
        let j23 = 2.0 * qz;
        let j24 = 2.0 * qy;
        let j31 = 0.0;
        let j32 = 4.0 * qx;
        let j33 = 4.0 * qy;
        let j34 = 0.0;
        let step_w = j11 * f1 + j21 * f2 + j31 * f3;
        let step_x = j12 * f1 + j22 * f2 + j32 * f3;
        let step_y = j13 * f1 + j23 * f2 + j33 * f3;
        let step_z = j14 * f1 + j24 * f2 + j34 * f3;
        let norm_s = (step_w * step_w + step_x * step_x + step_y * step_y + step_z * step_z)
            .sqrt()
            .max(1e-10);
        let (sw, sx, sy, sz) = (
            step_w / norm_s,
            step_x / norm_s,
            step_y / norm_s,
            step_z / norm_s,
        );

        // Rate of change from gyro + correction
        let qdotw = 0.5 * (-qx * gx - qy * gy - qz * gz) - self.beta * sw;
        let qdotx = 0.5 * (qw * gx + qy * gz - qz * gy) - self.beta * sx;
        let qdoty = 0.5 * (qw * gy - qx * gz + qz * gx) - self.beta * sy;
        let qdotz = 0.5 * (qw * gz + qx * gy - qy * gx) - self.beta * sz;

        qw += qdotw * dt;
        qx += qdotx * dt;
        qy += qdoty * dt;
        qz += qdotz * dt;
        let norm_q = (qw * qw + qx * qx + qy * qy + qz * qz).sqrt().max(1e-10);
        self.quaternion = [qw / norm_q, qx / norm_q, qy / norm_q, qz / norm_q];
        self.euler = quat_to_euler(self.quaternion);
    }

    /// Get roll angle in degrees.
    pub fn roll(&self) -> f64 {
        self.euler[0]
    }
    /// Get pitch angle in degrees.
    pub fn pitch(&self) -> f64 {
        self.euler[1]
    }
    /// Get yaw angle in degrees.
    pub fn yaw(&self) -> f64 {
        self.euler[2]
    }
}

/// Ring-buffer data logger with configurable size.
///
/// Flushes to in-memory buffer when full.
#[derive(Debug, Clone)]
pub struct DataLogger<T: Clone> {
    /// Ring buffer.
    buffer: VecDeque<T>,
    /// Maximum capacity.
    capacity: usize,
    /// Flushed records.
    flushed: Vec<T>,
}

impl<T: Clone> DataLogger<T> {
    /// Create a new data logger with given capacity.
    pub fn new(capacity: usize) -> Self {
        Self {
            buffer: VecDeque::with_capacity(capacity),
            capacity,
            flushed: Vec::new(),
        }
    }

    /// Log a record.
    pub fn log(&mut self, record: T) {
        if self.buffer.len() >= self.capacity {
            self.flush();
        }
        self.buffer.push_back(record);
    }

    /// Flush buffer to flushed storage.
    pub fn flush(&mut self) {
        self.flushed.extend(self.buffer.drain(..));
    }

    /// Total number of logged records.
    pub fn total_count(&self) -> usize {
        self.buffer.len() + self.flushed.len()
    }

    /// Current buffer size.
    pub fn buffer_size(&self) -> usize {
        self.buffer.len()
    }

    /// Access all flushed records.
    pub fn flushed_records(&self) -> &[T] {
        &self.flushed
    }

    /// Pop oldest record.
    pub fn pop(&mut self) -> Option<T> {
        self.buffer.pop_front()
    }

    /// Replay from beginning (returns cloned records).
    pub fn replay(&self) -> Vec<T> {
        let mut all = self.flushed.clone();
        all.extend(self.buffer.iter().cloned());
        all
    }
}

/// Sensor data writer for binary output.
#[derive(Debug)]
pub struct SensorDataWriter {
    /// Output buffer.
    pub buffer: Vec<u8>,
    /// Number of records written.
    pub n_records: usize,
    /// Sensor type tag.
    pub sensor_type: SensorType,
}

/// Sensor type tag.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum SensorType {
    /// IMU sensor.
    Imu,
    /// GPS sensor.
    Gps,
    /// LiDAR sensor.
    Lidar,
    /// Force gauge sensor.
    ForceGauge,
}

impl SensorDataWriter {
    /// Create a new sensor data writer.
    pub fn new(sensor_type: SensorType) -> Self {
        Self {
            buffer: Vec::new(),
            n_records: 0,
            sensor_type,
        }
    }

    /// Write IMU record.
    pub fn write_imu(&mut self, record: &ImuDataRecord) {
        self.buffer.extend_from_slice(&record.to_bytes());
        self.n_records += 1;
    }

    /// Write GPS record.
    pub fn write_gps(&mut self, record: &GpsDataRecord) {
        self.buffer.extend_from_slice(&record.to_bytes());
        self.n_records += 1;
    }

    /// Write force gauge record.
    pub fn write_force(&mut self, record: &ForceGaugeRecord) {
        self.buffer.extend_from_slice(&record.to_bytes());
        self.n_records += 1;
    }

    /// Write LiDAR frame.
    pub fn write_lidar(&mut self, frame: &LidarFrame) {
        self.buffer.extend_from_slice(&frame.to_bytes());
        self.n_records += 1;
    }

    /// Clear buffer.
    pub fn clear(&mut self) {
        self.buffer.clear();
        self.n_records = 0;
    }

    /// Total bytes written.
    pub fn bytes_written(&self) -> usize {
        self.buffer.len()
    }
}

/// Sensor data reader with timestamp synchronization and interpolation.
#[derive(Debug)]
pub struct SensorDataReader {
    /// Raw data buffer.
    pub data: Vec<u8>,
    /// Sensor type.
    pub sensor_type: SensorType,
    /// Read position.
    pub pos: usize,
}

impl SensorDataReader {
    /// Create a new sensor data reader.
    pub fn new(data: Vec<u8>, sensor_type: SensorType) -> Self {
        Self {
            data,
            sensor_type,
            pos: 0,
        }
    }

    /// Read next IMU record.
    pub fn read_imu(&mut self) -> Option<ImuDataRecord> {
        let record = ImuDataRecord::from_bytes(&self.data[self.pos..])?;
        self.pos += 80;
        Some(record)
    }

    /// Check if more data available.
    pub fn has_more(&self) -> bool {
        self.pos < self.data.len()
    }

    /// Seek to position.
    pub fn seek(&mut self, pos: usize) {
        self.pos = pos.min(self.data.len());
    }

    /// Remaining bytes.
    pub fn remaining(&self) -> usize {
        self.data.len().saturating_sub(self.pos)
    }
}

/// Multi-column CSV telemetry writer.
///
/// Efficient append with header row.
#[derive(Debug, Clone)]
pub struct CsvTelemetry {
    /// CSV header columns.
    pub headers: Vec<String>,
    /// Data rows.
    pub rows: Vec<Vec<String>>,
}

impl CsvTelemetry {
    /// Create a new CSV telemetry writer.
    pub fn new(headers: Vec<String>) -> Self {
        Self {
            headers,
            rows: Vec::new(),
        }
    }

    /// Append a row.
    pub fn append_row(&mut self, values: Vec<String>) {
        self.rows.push(values);
    }

    /// Append IMU record as row.
    pub fn append_imu(&mut self, r: &ImuDataRecord) {
        self.rows.push(vec![
            r.timestamp.to_string(),
            r.accel[0].to_string(),
            r.accel[1].to_string(),
            r.accel[2].to_string(),
            r.gyro[0].to_string(),
            r.gyro[1].to_string(),
            r.gyro[2].to_string(),
        ]);
    }

    /// Serialize to CSV string.
    pub fn to_csv_string(&self) -> String {
        let mut out = self.headers.join(",") + "\n";
        for row in &self.rows {
            out += &(row.join(",") + "\n");
        }
        out
    }

    /// Number of data rows.
    pub fn n_rows(&self) -> usize {
        self.rows.len()
    }

    /// Write to a BufWriter.
    pub fn write_to<W: Write>(&self, writer: &mut BufWriter<W>) -> std::io::Result<()> {
        writer.write_all(self.to_csv_string().as_bytes())
    }
}

// --- Helper functions ---

/// Convert Euler angles (rad) to quaternion \[w, x, y, z\].
fn euler_to_quat(roll: f64, pitch: f64, yaw: f64) -> [f64; 4] {
    let (cr, sr) = ((roll * 0.5).cos(), (roll * 0.5).sin());
    let (cp, sp) = ((pitch * 0.5).cos(), (pitch * 0.5).sin());
    let (cy, sy) = ((yaw * 0.5).cos(), (yaw * 0.5).sin());
    [
        cr * cp * cy + sr * sp * sy,
        sr * cp * cy - cr * sp * sy,
        cr * sp * cy + sr * cp * sy,
        cr * cp * sy - sr * sp * cy,
    ]
}

/// Convert quaternion to Euler angles \[roll, pitch, yaw\] in degrees.
fn quat_to_euler(q: [f64; 4]) -> [f64; 3] {
    let [w, x, y, z] = q;
    let roll = (2.0 * (w * x + y * z))
        .atan2(1.0 - 2.0 * (x * x + y * y))
        .to_degrees();
    let sin_pitch = (2.0 * (w * y - z * x)).clamp(-1.0, 1.0);
    let pitch = sin_pitch.asin().to_degrees();
    let yaw = (2.0 * (w * z + x * y))
        .atan2(1.0 - 2.0 * (y * y + z * z))
        .to_degrees();
    [roll, pitch, yaw]
}

/// Spherical linear interpolation between two quaternions.
fn quat_slerp(q1: [f64; 4], q2: [f64; 4], t: f64) -> [f64; 4] {
    let dot = q1[0] * q2[0] + q1[1] * q2[1] + q1[2] * q2[2] + q1[3] * q2[3];
    let dot = dot.clamp(-1.0, 1.0);
    if dot.abs() > 0.9995 {
        // Linear interpolation for nearly identical quaternions
        let s = 1.0 - t;
        let q = [
            s * q1[0] + t * q2[0],
            s * q1[1] + t * q2[1],
            s * q1[2] + t * q2[2],
            s * q1[3] + t * q2[3],
        ];
        let norm = (q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3])
            .sqrt()
            .max(1e-10);
        return [q[0] / norm, q[1] / norm, q[2] / norm, q[3] / norm];
    }
    let theta = dot.acos();
    let sin_theta = theta.sin().max(1e-10);
    let s1 = ((1.0 - t) * theta).sin() / sin_theta;
    let s2 = (t * theta).sin() / sin_theta;
    [
        s1 * q1[0] + s2 * q2[0],
        s1 * q1[1] + s2 * q2[1],
        s1 * q1[2] + s2 * q2[2],
        s1 * q1[3] + s2 * q2[3],
    ]
}

// ============================================================
// Tests
// ============================================================
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_imu_serialize_deserialize() {
        let r = ImuDataRecord::new(
            1.0,
            [0.1, 0.2, 9.81],
            [0.01, 0.02, 0.03],
            [20.0, 30.0, 40.0],
        );
        let bytes = r.to_bytes();
        assert_eq!(bytes.len(), 80);
        let r2 = ImuDataRecord::from_bytes(&bytes).unwrap();
        assert!((r2.timestamp - 1.0).abs() < 1e-15);
        assert!((r2.accel[2] - 9.81).abs() < 1e-10);
    }

    #[test]
    fn test_imu_from_bytes_too_short() {
        let result = ImuDataRecord::from_bytes(&[0u8; 10]);
        assert!(result.is_none());
    }

    #[test]
    fn test_gps_serialize() {
        let r = GpsDataRecord::new(0.5, 35.0, 139.0, 100.0, [1.0, 2.0, 0.0], 1.2);
        let bytes = r.to_bytes();
        assert_eq!(bytes.len(), 64);
    }

    #[test]
    fn test_gps_ground_speed() {
        let r = GpsDataRecord::new(0.0, 0.0, 0.0, 0.0, [3.0, 4.0, 0.0], 1.0);
        assert!((r.ground_speed() - 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_lidar_frame_add_point() {
        let mut f = LidarFrame::new(1.0);
        f.add_point([1.0, 2.0, 3.0], 0.5);
        assert_eq!(f.n_points, 1);
    }

    #[test]
    fn test_lidar_frame_serialize() {
        let mut f = LidarFrame::new(1.0);
        f.add_point([1.0, 0.0, 0.0], 0.8);
        let bytes = f.to_bytes();
        assert!(!bytes.is_empty());
    }

    #[test]
    fn test_lidar_bounding_box() {
        let mut f = LidarFrame::new(0.0);
        f.add_point([1.0, 2.0, 3.0], 1.0);
        f.add_point([-1.0, -2.0, -3.0], 0.5);
        let bb = f.bounding_box();
        assert!((bb[0][0] + 1.0).abs() < 1e-6);
        assert!((bb[0][1] - 1.0).abs() < 1e-6);
    }

    #[test]
    fn test_force_gauge_magnitude() {
        let r = ForceGaugeRecord::new(0.0, [3.0, 4.0, 0.0], [0.0; 3], 25.0);
        assert!((r.force_magnitude() - 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_force_gauge_serialize() {
        let r = ForceGaugeRecord::new(1.0, [1.0, 2.0, 3.0], [0.1, 0.2, 0.3], 22.5);
        let bytes = r.to_bytes();
        assert_eq!(bytes.len(), 64);
    }

    #[test]
    fn test_telemetry_encode_decode() {
        let proto = TelemetryProtocol::new(1, 1);
        let msg = TelemetryMessage {
            message_id: 100,
            seq: 0,
            payload: vec![1u8, 2, 3, 4],
        };
        let encoded = proto.encode(&msg);
        let decoded = proto.decode(&encoded).unwrap();
        assert_eq!(decoded.message_id, 100);
        assert_eq!(decoded.payload, vec![1u8, 2, 3, 4]);
    }

    #[test]
    fn test_telemetry_wrong_magic() {
        let proto = TelemetryProtocol::new(1, 1);
        let mut data = vec![0u8; 20];
        data[0] = 0x00; // wrong magic
        assert!(proto.decode(&data).is_none());
    }

    #[test]
    fn test_sensor_fusion_complementary() {
        let mut sf = SensorFusion::new(0.98, 0.1);
        let accel = [0.0, 0.0, 9.81];
        let gyro = [0.0, 0.0, 0.0];
        sf.update_complementary(accel, gyro, 0.01);
        assert!(sf.roll().is_finite());
        assert!(sf.pitch().is_finite());
    }

    #[test]
    fn test_sensor_fusion_madgwick() {
        let mut sf = SensorFusion::new(0.98, 0.1);
        let accel = [0.0, 0.0, 9.81];
        let gyro = [0.01, 0.02, 0.0];
        sf.update_madgwick(accel, gyro, 0.01);
        let norm = sf.quaternion.iter().map(|v| v * v).sum::<f64>().sqrt();
        assert!((norm - 1.0).abs() < 0.01);
    }

    #[test]
    fn test_data_logger_log_and_flush() {
        let mut logger: DataLogger<i32> = DataLogger::new(3);
        logger.log(1);
        logger.log(2);
        logger.log(3);
        logger.log(4); // triggers flush
        assert_eq!(logger.total_count(), 4);
    }

    #[test]
    fn test_data_logger_replay() {
        let mut logger: DataLogger<i32> = DataLogger::new(10);
        for i in 0..5 {
            logger.log(i);
        }
        let all = logger.replay();
        assert_eq!(all.len(), 5);
    }

    #[test]
    fn test_sensor_data_writer_imu() {
        let mut w = SensorDataWriter::new(SensorType::Imu);
        let r = ImuDataRecord::new(0.0, [0.0; 3], [0.0; 3], [0.0; 3]);
        w.write_imu(&r);
        assert_eq!(w.n_records, 1);
        assert_eq!(w.bytes_written(), 80);
    }

    #[test]
    fn test_sensor_data_writer_clear() {
        let mut w = SensorDataWriter::new(SensorType::Gps);
        let r = GpsDataRecord::new(0.0, 0.0, 0.0, 0.0, [0.0; 3], 1.0);
        w.write_gps(&r);
        w.clear();
        assert_eq!(w.n_records, 0);
        assert_eq!(w.bytes_written(), 0);
    }

    #[test]
    fn test_sensor_data_reader_imu() {
        let r = ImuDataRecord::new(1.5, [1.0, 2.0, 3.0], [0.1, 0.2, 0.3], [10.0, 20.0, 30.0]);
        let bytes = r.to_bytes();
        let mut reader = SensorDataReader::new(bytes, SensorType::Imu);
        let r2 = reader.read_imu().unwrap();
        assert!((r2.timestamp - 1.5).abs() < 1e-10);
    }

    #[test]
    fn test_csv_telemetry_to_string() {
        let mut csv = CsvTelemetry::new(vec!["t".into(), "ax".into(), "ay".into()]);
        csv.append_row(vec!["0.0".into(), "1.0".into(), "2.0".into()]);
        let s = csv.to_csv_string();
        assert!(s.contains("t,ax,ay"));
        assert!(s.contains("0.0,1.0,2.0"));
    }

    #[test]
    fn test_csv_telemetry_n_rows() {
        let mut csv = CsvTelemetry::new(vec!["x".into()]);
        csv.append_row(vec!["1".into()]);
        csv.append_row(vec!["2".into()]);
        assert_eq!(csv.n_rows(), 2);
    }

    #[test]
    fn test_quat_to_euler_identity() {
        let euler = quat_to_euler([1.0, 0.0, 0.0, 0.0]);
        for &v in &euler {
            assert!(v.abs() < 1e-8, "Expected near zero, got {v}");
        }
    }
}

// ============================================================
// New sensor time-series I/O  (added block)
// ============================================================

/// Single sensor reading with uncertainty.
#[derive(Debug, Clone, PartialEq)]
pub struct SensorReading {
    /// Timestamp in seconds.
    pub timestamp: f64,
    /// Sensor ID.
    pub sensor_id: usize,
    /// Measured value (arbitrary units).
    pub value: f64,
    /// One-sigma measurement uncertainty.
    pub uncertainty: f64,
}

impl SensorReading {
    /// Create a new sensor reading.
    pub fn new(timestamp: f64, sensor_id: usize, value: f64, uncertainty: f64) -> Self {
        Self {
            timestamp,
            sensor_id,
            value,
            uncertainty,
        }
    }
}

/// Time series of readings for one sensor.
#[derive(Debug, Clone)]
pub struct SensorTimeSeries {
    /// Sensor identifier.
    pub id: usize,
    /// Human-readable sensor name.
    pub name: String,
    /// Ordered readings (ascending timestamp assumed).
    pub readings: Vec<SensorReading>,
    /// Nominal sample rate in Hz.
    pub sample_rate_hz: f64,
}

impl SensorTimeSeries {
    /// Create a new empty time series.
    pub fn new(id: usize, name: impl Into<String>, sample_rate_hz: f64) -> Self {
        Self {
            id,
            name: name.into(),
            readings: Vec::new(),
            sample_rate_hz,
        }
    }

    /// Append a reading.
    pub fn push(&mut self, r: SensorReading) {
        self.readings.push(r);
    }

    /// Number of readings.
    pub fn len(&self) -> usize {
        self.readings.len()
    }

    /// True when no readings.
    pub fn is_empty(&self) -> bool {
        self.readings.is_empty()
    }
}

/// Write a `SensorTimeSeries` to a two-column CSV file `timestamp,value`.
///
/// Returns an error string on I/O failure.
pub fn write_sensor_csv(series: &SensorTimeSeries, path: &str) -> Result<(), String> {
    use std::io::Write as IoWrite;
    let file = std::fs::File::create(path).map_err(|e| e.to_string())?;
    let mut writer = std::io::BufWriter::new(file);
    writeln!(writer, "timestamp,value,uncertainty").map_err(|e| e.to_string())?;
    for r in &series.readings {
        writeln!(writer, "{},{},{}", r.timestamp, r.value, r.uncertainty)
            .map_err(|e| e.to_string())?;
    }
    Ok(())
}

/// Read a CSV file written by [`write_sensor_csv`] into a `SensorTimeSeries`.
pub fn read_sensor_csv(path: &str, sensor_id: usize) -> Result<SensorTimeSeries, String> {
    let content = std::fs::read_to_string(path).map_err(|e| e.to_string())?;
    let mut lines = content.lines();
    let _header = lines.next().ok_or("missing header")?;
    let mut series = SensorTimeSeries::new(sensor_id, path, 0.0);
    for line in lines {
        let parts: Vec<&str> = line.split(',').collect();
        if parts.len() < 2 {
            continue;
        }
        let ts: f64 = parts[0]
            .trim()
            .parse()
            .map_err(|_| format!("bad timestamp: {}", parts[0]))?;
        let val: f64 = parts[1]
            .trim()
            .parse()
            .map_err(|_| format!("bad value: {}", parts[1]))?;
        let unc: f64 = if parts.len() >= 3 {
            parts[2].trim().parse().unwrap_or(0.0)
        } else {
            0.0
        };
        series.push(SensorReading::new(ts, sensor_id, val, unc));
    }
    Ok(series)
}

/// Compute a causal moving average of sensor values with the given window size.
///
/// The first `window-1` values use a partial window.
pub fn sensor_moving_average(series: &SensorTimeSeries, window: usize) -> SensorTimeSeries {
    let w = window.max(1);
    let mut out = SensorTimeSeries::new(series.id, series.name.clone(), series.sample_rate_hz);
    for i in 0..series.readings.len() {
        let start = (i + 1).saturating_sub(w);
        let slice = &series.readings[start..=i];
        let avg = slice.iter().map(|r| r.value).sum::<f64>() / slice.len() as f64;
        let mut r = series.readings[i].clone();
        r.value = avg;
        out.push(r);
    }
    out
}

/// Downsample a series by keeping every `factor`-th reading.
///
/// A factor of 1 returns a clone of the original series.
pub fn sensor_downsample(series: &SensorTimeSeries, factor: usize) -> SensorTimeSeries {
    let f = factor.max(1);
    let mut out = SensorTimeSeries::new(
        series.id,
        series.name.clone(),
        series.sample_rate_hz / f as f64,
    );
    for (i, r) in series.readings.iter().enumerate() {
        if i % f == 0 {
            out.push(r.clone());
        }
    }
    out
}

/// Find indices of local maxima above `threshold` with at least `min_gap` samples between peaks.
pub fn sensor_peak_detection(
    series: &SensorTimeSeries,
    threshold: f64,
    min_gap: usize,
) -> Vec<usize> {
    let vals: Vec<f64> = series.readings.iter().map(|r| r.value).collect();
    let n = vals.len();
    let mut peaks = Vec::new();
    let mut last_peak: Option<usize> = None;
    for i in 1..n.saturating_sub(1) {
        if vals[i] > threshold && vals[i] >= vals[i - 1] && vals[i] >= vals[i + 1] {
            if let Some(lp) = last_peak
                && i - lp < min_gap
            {
                continue;
            }
            peaks.push(i);
            last_peak = Some(i);
        }
    }
    peaks
}

/// Compute the DFT magnitude spectrum of a sensor time series.
///
/// Returns a vector of `(frequency_hz, magnitude)` pairs for bins 0..N/2.
pub fn sensor_fft_spectrum(series: &SensorTimeSeries) -> Vec<(f64, f64)> {
    let n = series.readings.len();
    if n == 0 {
        return Vec::new();
    }
    let vals: Vec<f64> = series.readings.iter().map(|r| r.value).collect();
    let dt = if series.sample_rate_hz > 0.0 {
        1.0 / series.sample_rate_hz
    } else {
        1.0
    };
    let half = n / 2 + 1;
    let mut spectrum = Vec::with_capacity(half);
    for k in 0..half {
        let mut re = 0.0f64;
        let mut im = 0.0f64;
        for (t, &v) in vals.iter().enumerate() {
            let angle = -2.0 * std::f64::consts::PI * (k as f64) * (t as f64) / (n as f64);
            re += v * angle.cos();
            im += v * angle.sin();
        }
        let mag = (re * re + im * im).sqrt();
        let freq = (k as f64) / (n as f64 * dt);
        spectrum.push((freq, mag));
    }
    spectrum
}

/// Apply a linear calibration to a raw sensor value: `gain * raw + offset`.
pub fn sensor_calibration_linear(raw: f64, gain: f64, offset: f64) -> f64 {
    gain * raw + offset
}

/// Resample all series in `series` to the same `target_rate` using nearest-neighbour interpolation.
///
/// Returns one `SensorTimeSeries` per input, all with the same length.
pub fn sensor_synchronize(series: &[SensorTimeSeries], target_rate: f64) -> Vec<SensorTimeSeries> {
    if series.is_empty() || target_rate <= 0.0 {
        return Vec::new();
    }
    // Determine common time span.
    let t_start = series
        .iter()
        .filter_map(|s| s.readings.first().map(|r| r.timestamp))
        .fold(f64::INFINITY, f64::min);
    let t_end = series
        .iter()
        .filter_map(|s| s.readings.last().map(|r| r.timestamp))
        .fold(f64::NEG_INFINITY, f64::max);
    if t_start >= t_end {
        return series
            .iter()
            .map(|s| SensorTimeSeries::new(s.id, s.name.clone(), target_rate))
            .collect();
    }
    let n_samples = ((t_end - t_start) * target_rate).ceil() as usize + 1;
    let timestamps: Vec<f64> = (0..n_samples)
        .map(|i| t_start + i as f64 / target_rate)
        .collect();

    series
        .iter()
        .map(|s| {
            let mut out = SensorTimeSeries::new(s.id, s.name.clone(), target_rate);
            for &ts in &timestamps {
                // Nearest-neighbour lookup.
                let val = if s.readings.is_empty() {
                    0.0
                } else {
                    s.readings
                        .iter()
                        .min_by(|a, b| {
                            (a.timestamp - ts)
                                .abs()
                                .partial_cmp(&(b.timestamp - ts).abs())
                                .expect("operation should succeed")
                        })
                        .map(|r| r.value)
                        .unwrap_or(0.0)
                };
                out.push(SensorReading::new(ts, s.id, val, 0.0));
            }
            out
        })
        .collect()
}

/// Build a matrix (row = time step, col = channel) from synchronised channels.
///
/// All channels must have the same length; shorter channels are zero-padded.
pub fn merge_sensor_channels(channels: &[SensorTimeSeries]) -> Vec<Vec<f64>> {
    if channels.is_empty() {
        return Vec::new();
    }
    let n_rows = channels.iter().map(|c| c.len()).max().unwrap_or(0);
    (0..n_rows)
        .map(|i| {
            channels
                .iter()
                .map(|c| c.readings.get(i).map(|r| r.value).unwrap_or(0.0))
                .collect()
        })
        .collect()
}

// ── sensor_io tests ──────────────────────────────────────────────────────────
#[cfg(test)]
mod sensor_ts_tests {
    use super::*;

    fn make_series(n: usize, sample_rate_hz: f64) -> SensorTimeSeries {
        let mut s = SensorTimeSeries::new(1, "test", sample_rate_hz);
        for i in 0..n {
            s.push(SensorReading::new(
                i as f64 / sample_rate_hz,
                1,
                i as f64,
                0.1,
            ));
        }
        s
    }

    #[test]
    fn test_sensor_reading_new() {
        let r = SensorReading::new(1.0, 2, 3.5, 0.1);
        assert!((r.timestamp - 1.0).abs() < 1e-15);
        assert_eq!(r.sensor_id, 2);
        assert!((r.value - 3.5).abs() < 1e-15);
    }

    #[test]
    fn test_series_len() {
        let s = make_series(10, 100.0);
        assert_eq!(s.len(), 10);
    }

    #[test]
    fn test_series_is_empty() {
        let s = SensorTimeSeries::new(0, "empty", 100.0);
        assert!(s.is_empty());
    }

    #[test]
    fn test_write_read_csv_roundtrip() {
        let s = make_series(5, 10.0);
        let path = "/tmp/oxiphysics_sensor_test.csv";
        write_sensor_csv(&s, path).unwrap();
        let s2 = read_sensor_csv(path, 1).unwrap();
        assert_eq!(s2.len(), 5);
        for (a, b) in s.readings.iter().zip(s2.readings.iter()) {
            assert!((a.timestamp - b.timestamp).abs() < 1e-10);
            assert!((a.value - b.value).abs() < 1e-10);
        }
    }

    #[test]
    fn test_write_csv_creates_file() {
        let s = make_series(3, 1.0);
        let path = "/tmp/oxiphysics_sensor_create_test.csv";
        assert!(write_sensor_csv(&s, path).is_ok());
    }

    #[test]
    fn test_read_csv_nonexistent() {
        let res = read_sensor_csv("/tmp/nonexistent_oxiphysics_xyz.csv", 0);
        assert!(res.is_err());
    }

    #[test]
    fn test_moving_average_same_length() {
        let s = make_series(10, 10.0);
        let avg = sensor_moving_average(&s, 3);
        assert_eq!(avg.len(), s.len());
    }

    #[test]
    fn test_moving_average_smooths() {
        // Noisy signal: alternating 0, 10
        let mut s = SensorTimeSeries::new(0, "noisy", 1.0);
        for i in 0..10 {
            s.push(SensorReading::new(
                i as f64,
                0,
                if i % 2 == 0 { 0.0 } else { 10.0 },
                0.0,
            ));
        }
        let avg = sensor_moving_average(&s, 3);
        // Peak value should be less than original max.
        let max_avg = avg
            .readings
            .iter()
            .map(|r| r.value)
            .fold(f64::NEG_INFINITY, f64::max);
        assert!(max_avg < 10.0, "max_avg={max_avg}");
    }

    #[test]
    fn test_moving_average_window_1_identity() {
        let s = make_series(5, 10.0);
        let avg = sensor_moving_average(&s, 1);
        for (a, b) in s.readings.iter().zip(avg.readings.iter()) {
            assert!((a.value - b.value).abs() < 1e-12);
        }
    }

    #[test]
    fn test_downsample_reduces_count() {
        let s = make_series(10, 10.0);
        let ds = sensor_downsample(&s, 2);
        assert_eq!(ds.len(), 5);
    }

    #[test]
    fn test_downsample_factor_1_same_count() {
        let s = make_series(8, 10.0);
        let ds = sensor_downsample(&s, 1);
        assert_eq!(ds.len(), s.len());
    }

    #[test]
    fn test_downsample_sample_rate_reduced() {
        let s = make_series(10, 100.0);
        let ds = sensor_downsample(&s, 4);
        assert!((ds.sample_rate_hz - 25.0).abs() < 1e-10);
    }

    #[test]
    fn test_peak_detection_finds_peaks() {
        let mut s = SensorTimeSeries::new(0, "peaks", 1.0);
        // Create: 0 5 0 0 7 0 0 6 0
        let vals = [0.0f64, 5.0, 0.0, 0.0, 7.0, 0.0, 0.0, 6.0, 0.0];
        for (i, &v) in vals.iter().enumerate() {
            s.push(SensorReading::new(i as f64, 0, v, 0.0));
        }
        let peaks = sensor_peak_detection(&s, 3.0, 1);
        assert!(!peaks.is_empty(), "should find at least one peak");
        // Check that all found peaks are above threshold.
        for &p in &peaks {
            assert!(s.readings[p].value > 3.0);
        }
    }

    #[test]
    fn test_peak_detection_no_peaks_below_threshold() {
        let mut s = SensorTimeSeries::new(0, "flat", 1.0);
        for i in 0..5 {
            s.push(SensorReading::new(i as f64, 0, 1.0, 0.0));
        }
        let peaks = sensor_peak_detection(&s, 5.0, 1);
        assert!(peaks.is_empty());
    }

    #[test]
    fn test_peak_detection_min_gap_enforced() {
        let mut s = SensorTimeSeries::new(0, "close_peaks", 1.0);
        let vals = [0.0f64, 5.0, 0.0, 5.0, 0.0, 5.0, 0.0];
        for (i, &v) in vals.iter().enumerate() {
            s.push(SensorReading::new(i as f64, 0, v, 0.0));
        }
        let peaks_close = sensor_peak_detection(&s, 3.0, 1);
        let peaks_far = sensor_peak_detection(&s, 3.0, 3);
        assert!(peaks_far.len() <= peaks_close.len());
    }

    #[test]
    fn test_fft_spectrum_positive_magnitude() {
        let s = make_series(16, 16.0);
        let spec = sensor_fft_spectrum(&s);
        assert!(!spec.is_empty());
        for (_, mag) in &spec {
            assert!(*mag >= 0.0);
        }
    }

    #[test]
    fn test_fft_spectrum_empty_series() {
        let s = SensorTimeSeries::new(0, "empty", 100.0);
        let spec = sensor_fft_spectrum(&s);
        assert!(spec.is_empty());
    }

    #[test]
    fn test_fft_spectrum_length() {
        let s = make_series(10, 10.0);
        let spec = sensor_fft_spectrum(&s);
        assert_eq!(spec.len(), 6); // N/2 + 1 = 5 + 1
    }

    #[test]
    fn test_calibration_linear_gain_offset() {
        let cal = sensor_calibration_linear(5.0, 2.0, 1.0);
        assert!((cal - 11.0).abs() < 1e-12);
    }

    #[test]
    fn test_calibration_linear_identity() {
        let val = sensor_calibration_linear(3.0, 1.0, 0.0);
        assert!((val - 3.0).abs() < 1e-12);
    }

    #[test]
    fn test_calibration_linear_zero_gain() {
        let val = sensor_calibration_linear(100.0, 0.0, 5.0);
        assert!((val - 5.0).abs() < 1e-12);
    }

    #[test]
    fn test_synchronize_equal_length() {
        let s1 = make_series(10, 10.0);
        let s2 = make_series(8, 8.0);
        let synced = sensor_synchronize(&[s1, s2], 10.0);
        assert_eq!(synced.len(), 2);
        assert_eq!(synced[0].len(), synced[1].len());
    }

    #[test]
    fn test_synchronize_empty_input() {
        let synced = sensor_synchronize(&[], 10.0);
        assert!(synced.is_empty());
    }

    #[test]
    fn test_synchronize_zero_rate() {
        let s = make_series(5, 10.0);
        let synced = sensor_synchronize(&[s], 0.0);
        assert!(synced.is_empty());
    }

    #[test]
    fn test_merge_sensor_channels_shape() {
        let s1 = make_series(5, 10.0);
        let s2 = make_series(5, 10.0);
        let mat = merge_sensor_channels(&[s1, s2]);
        assert_eq!(mat.len(), 5);
        for row in &mat {
            assert_eq!(row.len(), 2);
        }
    }

    #[test]
    fn test_merge_sensor_channels_empty() {
        let mat = merge_sensor_channels(&[]);
        assert!(mat.is_empty());
    }

    #[test]
    fn test_merge_sensor_channels_values() {
        let s1 = make_series(3, 1.0);
        let mat = merge_sensor_channels(&[s1]);
        assert!((mat[0][0] - 0.0).abs() < 1e-10);
        assert!((mat[1][0] - 1.0).abs() < 1e-10);
        assert!((mat[2][0] - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_csv_roundtrip_larger_series() {
        let s = make_series(100, 50.0);
        let path = "/tmp/oxiphysics_sensor_large.csv";
        write_sensor_csv(&s, path).unwrap();
        let s2 = read_sensor_csv(path, 1).unwrap();
        assert_eq!(s2.len(), 100);
    }

    #[test]
    fn test_sensor_time_series_push() {
        let mut s = SensorTimeSeries::new(0, "x", 1.0);
        s.push(SensorReading::new(0.0, 0, 42.0, 0.0));
        assert_eq!(s.len(), 1);
        assert!(!s.is_empty());
    }

    #[test]
    fn test_fft_dc_component_is_sum() {
        // DC bin (k=0) magnitude should equal |sum of values|.
        let mut s = SensorTimeSeries::new(0, "dc", 4.0);
        let vals = [2.0f64, 2.0, 2.0, 2.0];
        for (i, &v) in vals.iter().enumerate() {
            s.push(SensorReading::new(i as f64 / 4.0, 0, v, 0.0));
        }
        let spec = sensor_fft_spectrum(&s);
        // DC magnitude = |sum| = 8.0
        assert!((spec[0].1 - 8.0).abs() < 1e-8);
    }
}
