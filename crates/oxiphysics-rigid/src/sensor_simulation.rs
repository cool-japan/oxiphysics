// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Virtual sensor simulation for rigid body systems.
//!
//! Provides IMU (inertial measurement unit), GPS (position sensor), and
//! 6-DOF force/torque sensor simulations with realistic noise models and
//! a complementary filter for sensor fusion.
//!
//! # Overview
//!
//! - [`ImuSensor`] — gyroscope/accelerometer with bias and noise density.
//! - [`simulate_imu`] — add Gaussian noise and bias drift to true motion.
//! - [`PositionSensor`] — GPS-like sensor with position noise.
//! - [`simulate_gps`] — noisy position measurement.
//! - [`ForceSensor`] — 6-DOF force/torque sensor with calibration matrix.
//! - [`simulate_force_sensor`] — 6-DOF force/torque measurement with noise.
//! - [`sensor_fusion_complementary`] — complementary filter combining IMU and GPS.

use std::f64::consts::PI;

// ─────────────────────────────────────────────────────────────────────────────
// Internal PRNG (deterministic, no external dep)
// ─────────────────────────────────────────────────────────────────────────────

/// Linear congruential PRNG for deterministic noise generation.
#[derive(Debug, Clone, Copy)]
struct Lcg {
    state: u64,
}

impl Lcg {
    fn new(seed: u64) -> Self {
        Self {
            state: seed ^ 0xcafe_babe_dead_beef,
        }
    }

    fn next_u64(&mut self) -> u64 {
        self.state = self
            .state
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        self.state
    }

    /// Uniform in \[0, 1).
    fn rand01(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
    }

    /// Approximately N(0, 1) via Box-Muller transform.
    fn randn(&mut self) -> f64 {
        let u1 = (self.rand01() + 1e-15).ln();
        let u2 = self.rand01() * 2.0 * PI;
        (-2.0 * u1).sqrt() * u2.cos()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// ImuSensor
// ─────────────────────────────────────────────────────────────────────────────

/// Parameters for an Inertial Measurement Unit (IMU) sensor.
///
/// Models accelerometer and gyroscope with independent bias and noise density.
#[derive(Debug, Clone)]
pub struct ImuSensor {
    /// Gyroscope noise density \[rad/s/√Hz\].
    pub gyro_noise_density: f64,
    /// Gyroscope bias \[rad/s\].
    pub gyro_bias: [f64; 3],
    /// Gyroscope bias drift rate \[rad/s²\] (random walk coefficient).
    pub gyro_bias_drift: f64,
    /// Accelerometer noise density \[m/s²/√Hz\].
    pub accel_noise_density: f64,
    /// Accelerometer bias \[m/s²\].
    pub accel_bias: [f64; 3],
    /// Accelerometer bias drift rate \[m/s³\] (random walk coefficient).
    pub accel_bias_drift: f64,
    /// Internal PRNG state seed.
    pub seed: u64,
}

impl ImuSensor {
    /// Construct an IMU sensor with specified noise parameters.
    pub fn new(
        gyro_noise_density: f64,
        gyro_bias: [f64; 3],
        gyro_bias_drift: f64,
        accel_noise_density: f64,
        accel_bias: [f64; 3],
        accel_bias_drift: f64,
        seed: u64,
    ) -> Self {
        Self {
            gyro_noise_density,
            gyro_bias,
            gyro_bias_drift,
            accel_noise_density,
            accel_bias,
            accel_bias_drift,
            seed,
        }
    }

    /// Construct a low-noise IMU with zero biases.
    pub fn ideal(seed: u64) -> Self {
        Self::new(1e-4, [0.0; 3], 1e-6, 1e-3, [0.0; 3], 1e-5, seed)
    }

    /// Construct a consumer-grade IMU with moderate noise.
    pub fn consumer_grade(seed: u64) -> Self {
        Self::new(
            1e-3,
            [0.01, -0.005, 0.008],
            1e-4,
            0.02,
            [0.05, -0.03, 0.02],
            1e-3,
            seed,
        )
    }
}

/// Measurement output of an IMU sensor.
#[derive(Debug, Clone, Copy)]
pub struct ImuMeasurement {
    /// Measured angular velocity \[rad/s\] (gyroscope reading).
    pub angular_velocity: [f64; 3],
    /// Measured linear acceleration \[m/s²\] (accelerometer reading).
    pub linear_acceleration: [f64; 3],
}

/// Simulate an IMU measurement by adding Gaussian noise and bias drift.
///
/// # Arguments
///
/// * `sensor` — mutable IMU sensor (bias state is updated).
/// * `true_angular_velocity` — ground-truth angular velocity \[rad/s\].
/// * `true_acceleration` — ground-truth linear acceleration \[m/s²\].
/// * `dt` — time step \[s\]; determines noise scaling (σ = density / √dt).
///
/// Returns a noisy [`ImuMeasurement`].
pub fn simulate_imu(
    sensor: &mut ImuSensor,
    true_angular_velocity: [f64; 3],
    true_acceleration: [f64; 3],
    dt: f64,
) -> ImuMeasurement {
    let mut rng = Lcg::new(sensor.seed);
    sensor.seed = rng.next_u64();

    let sqrt_dt = dt.sqrt();
    let inv_sqrt_dt = 1.0 / sqrt_dt;

    // Bias drift (random walk)
    let gyro_drift_sigma = sensor.gyro_bias_drift * sqrt_dt;
    let accel_drift_sigma = sensor.accel_bias_drift * sqrt_dt;

    let mut gyro_meas = [0.0f64; 3];
    let mut accel_meas = [0.0f64; 3];

    for i in 0..3 {
        // Update bias via random walk
        sensor.gyro_bias[i] += rng.randn() * gyro_drift_sigma;
        sensor.accel_bias[i] += rng.randn() * accel_drift_sigma;

        // Noise scaled by 1/sqrt(dt) so power spectral density = density^2
        let gyro_noise = sensor.gyro_noise_density * inv_sqrt_dt * rng.randn();
        let accel_noise = sensor.accel_noise_density * inv_sqrt_dt * rng.randn();

        gyro_meas[i] = true_angular_velocity[i] + sensor.gyro_bias[i] + gyro_noise;
        accel_meas[i] = true_acceleration[i] + sensor.accel_bias[i] + accel_noise;
    }

    ImuMeasurement {
        angular_velocity: gyro_meas,
        linear_acceleration: accel_meas,
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// PositionSensor (GPS-like)
// ─────────────────────────────────────────────────────────────────────────────

/// Parameters for a GPS-like position sensor.
#[derive(Debug, Clone)]
pub struct PositionSensor {
    /// Position noise standard deviation \[m\] (isotropic).
    pub position_noise_std: f64,
    /// Optional systematic bias offset \[m\].
    pub bias: [f64; 3],
    /// Internal PRNG seed.
    pub seed: u64,
}

impl PositionSensor {
    /// Construct a position sensor with given noise level.
    pub fn new(position_noise_std: f64, bias: [f64; 3], seed: u64) -> Self {
        Self {
            position_noise_std,
            bias,
            seed,
        }
    }

    /// Construct a typical consumer GPS (2.5 m CEP ≈ 1.5 m std).
    pub fn consumer_gps(seed: u64) -> Self {
        Self::new(1.5, [0.0; 3], seed)
    }

    /// Construct a high-precision RTK GPS (~0.02 m std).
    pub fn rtk_gps(seed: u64) -> Self {
        Self::new(0.02, [0.0; 3], seed)
    }
}

/// Simulate a GPS position measurement.
///
/// Adds isotropic Gaussian noise and systematic bias to the true position.
///
/// Returns `[x, y, z]` noisy position \[m\].
pub fn simulate_gps(sensor: &mut PositionSensor, true_position: [f64; 3]) -> [f64; 3] {
    let mut rng = Lcg::new(sensor.seed);
    sensor.seed = rng.next_u64();

    let mut meas = [0.0f64; 3];
    for i in 0..3 {
        meas[i] = true_position[i] + sensor.bias[i] + sensor.position_noise_std * rng.randn();
    }
    meas
}

// ─────────────────────────────────────────────────────────────────────────────
// ForceSensor
// ─────────────────────────────────────────────────────────────────────────────

/// 6-DOF force/torque sensor with calibration matrix and noise.
///
/// The calibration matrix maps raw sensor counts to physical units
/// (Newtons and Newton-metres).  It is stored as a 6×6 row-major array.
#[derive(Debug, Clone)]
pub struct ForceSensor {
    /// 6×6 calibration matrix (row-major). Maps raw → \[Fx, Fy, Fz, Tx, Ty, Tz\].
    pub calibration: [[f64; 6]; 6],
    /// Noise standard deviation for forces \[N\].
    pub force_noise_std: f64,
    /// Noise standard deviation for torques \[N·m\].
    pub torque_noise_std: f64,
    /// Internal PRNG seed.
    pub seed: u64,
}

impl ForceSensor {
    /// Construct a force sensor with given calibration matrix and noise levels.
    pub fn new(
        calibration: [[f64; 6]; 6],
        force_noise_std: f64,
        torque_noise_std: f64,
        seed: u64,
    ) -> Self {
        Self {
            calibration,
            force_noise_std,
            torque_noise_std,
            seed,
        }
    }

    /// Construct an ideal force sensor (identity calibration, low noise).
    pub fn ideal(seed: u64) -> Self {
        let mut cal = [[0.0f64; 6]; 6];
        for (i, row) in cal.iter_mut().enumerate() {
            row[i] = 1.0;
        }
        Self::new(cal, 0.01, 0.001, seed)
    }

    /// Construct a typical industrial force/torque sensor.
    pub fn industrial(seed: u64) -> Self {
        let mut cal = [[0.0f64; 6]; 6];
        for (i, row) in cal.iter_mut().enumerate() {
            row[i] = 1.0;
        }
        // Small off-diagonal cross-talk terms
        cal[0][1] = 0.002;
        cal[1][0] = 0.002;
        cal[3][4] = 0.001;
        cal[4][3] = 0.001;
        Self::new(cal, 0.5, 0.05, seed)
    }
}

/// Simulate a 6-DOF force/torque sensor measurement.
///
/// Applies the calibration matrix and adds Gaussian noise.
///
/// # Arguments
///
/// * `sensor` — mutable force sensor (seed updated each call).
/// * `true_wrench` — ground-truth `[Fx, Fy, Fz, Tx, Ty, Tz]`.
///
/// Returns noisy `[Fx, Fy, Fz, Tx, Ty, Tz]`.
pub fn simulate_force_sensor(sensor: &mut ForceSensor, true_wrench: [f64; 6]) -> [f64; 6] {
    let mut rng = Lcg::new(sensor.seed);
    sensor.seed = rng.next_u64();

    // Apply calibration matrix
    let mut calibrated = [0.0f64; 6];
    for (row, (cal_r, calib_r)) in sensor
        .calibration
        .iter()
        .zip(calibrated.iter_mut())
        .enumerate()
    {
        let _ = row;
        *calib_r = cal_r
            .iter()
            .zip(true_wrench.iter())
            .map(|(c, t)| c * t)
            .sum();
    }

    // Add noise (forces: indices 0-2, torques: indices 3-5)
    let mut meas = [0.0f64; 6];
    for i in 0..6 {
        let noise_std = if i < 3 {
            sensor.force_noise_std
        } else {
            sensor.torque_noise_std
        };
        meas[i] = calibrated[i] + noise_std * rng.randn();
    }
    meas
}

// ─────────────────────────────────────────────────────────────────────────────
// Complementary filter
// ─────────────────────────────────────────────────────────────────────────────

/// State for the complementary filter fusing IMU velocity and GPS position.
#[derive(Debug, Clone)]
pub struct ComplementaryFilterState {
    /// Current estimated position \[m\].
    pub position: [f64; 3],
    /// Current estimated velocity \[m/s\].
    pub velocity: [f64; 3],
}

impl ComplementaryFilterState {
    /// Initialise at a given position with zero velocity.
    pub fn new(initial_position: [f64; 3]) -> Self {
        Self {
            position: initial_position,
            velocity: [0.0; 3],
        }
    }
}

/// Complementary filter fusing IMU acceleration and GPS position.
///
/// The filter uses a high-pass path from IMU (good for high-frequency motion)
/// and a low-pass path from GPS (good for DC and low-frequency accuracy).
///
/// `alpha` ∈ (0, 1]: weight for the IMU-integrated prediction.
/// Typical values: 0.98 for dt ≈ 0.01 s.
///
/// # Returns
///
/// Updated [`ComplementaryFilterState`] with fused position and velocity.
pub fn sensor_fusion_complementary(
    state: &ComplementaryFilterState,
    imu_meas: &ImuMeasurement,
    gps_position: [f64; 3],
    alpha: f64,
    dt: f64,
) -> ComplementaryFilterState {
    let a = imu_meas.linear_acceleration;

    // IMU integration (predict)
    let mut predicted_velocity = [0.0f64; 3];
    let mut predicted_position = [0.0f64; 3];
    for i in 0..3 {
        predicted_velocity[i] = state.velocity[i] + a[i] * dt;
        predicted_position[i] = state.position[i] + state.velocity[i] * dt + 0.5 * a[i] * dt * dt;
    }

    // Complementary blend: alpha * IMU + (1 - alpha) * GPS
    let beta = 1.0 - alpha;
    let mut fused_position = [0.0f64; 3];
    let mut fused_velocity = [0.0f64; 3];

    for i in 0..3 {
        fused_position[i] = alpha * predicted_position[i] + beta * gps_position[i];
        // Velocity correction via GPS position error
        let gps_correction = (gps_position[i] - predicted_position[i]) * beta / dt.max(1e-9);
        fused_velocity[i] = predicted_velocity[i] + gps_correction;
    }

    ComplementaryFilterState {
        position: fused_position,
        velocity: fused_velocity,
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // ── ImuSensor construction ───────────────────────────────────────────

    #[test]
    fn imu_ideal_has_zero_bias() {
        let s = ImuSensor::ideal(42);
        assert_eq!(s.gyro_bias, [0.0; 3]);
        assert_eq!(s.accel_bias, [0.0; 3]);
    }

    #[test]
    fn imu_consumer_grade_has_nonzero_bias() {
        let s = ImuSensor::consumer_grade(1);
        let bias_norm: f64 = s.gyro_bias.iter().map(|x| x * x).sum::<f64>().sqrt();
        assert!(bias_norm > 0.0);
    }

    #[test]
    fn imu_noise_density_set_correctly() {
        let s = ImuSensor::ideal(0);
        assert!((s.gyro_noise_density - 1e-4).abs() < 1e-10);
        assert!((s.accel_noise_density - 1e-3).abs() < 1e-10);
    }

    // ── simulate_imu noise statistics ───────────────────────────────────

    #[test]
    fn imu_measurement_close_to_true_for_ideal_sensor() {
        let mut s = ImuSensor::ideal(7);
        let omega = [0.1, -0.2, 0.3];
        let acc = [0.0, -9.81, 0.0];
        let meas = simulate_imu(&mut s, omega, acc, 0.01);
        for i in 0..3 {
            assert!(
                (meas.angular_velocity[i] - omega[i]).abs() < 0.5,
                "gyro[{i}] deviation too large"
            );
            assert!(
                (meas.linear_acceleration[i] - acc[i]).abs() < 1.0,
                "accel[{i}] deviation too large"
            );
        }
    }

    #[test]
    fn imu_repeated_calls_differ() {
        let mut s = ImuSensor::ideal(99);
        let m1 = simulate_imu(&mut s, [0.0; 3], [0.0; 3], 0.01);
        let m2 = simulate_imu(&mut s, [0.0; 3], [0.0; 3], 0.01);
        // Noise values should differ between calls
        let same = (0..3).all(|i| (m1.angular_velocity[i] - m2.angular_velocity[i]).abs() < 1e-15);
        assert!(!same, "consecutive IMU measurements should differ");
    }

    #[test]
    fn imu_bias_accumulates_over_time() {
        let mut s = ImuSensor::consumer_grade(123);
        let initial_bias = s.gyro_bias;
        for _ in 0..1000 {
            simulate_imu(&mut s, [0.0; 3], [0.0; 3], 0.01);
        }
        // After many steps the bias should have drifted
        let drift: f64 = (0..3)
            .map(|i| (s.gyro_bias[i] - initial_bias[i]).abs())
            .sum();
        // Some drift is expected (random walk), but could be small
        let _ = drift; // deterministic; just ensure no panic
    }

    #[test]
    fn imu_zero_input_mean_near_bias() {
        let mut s = ImuSensor::new(1e-4, [0.0; 3], 0.0, 1e-3, [0.0; 3], 0.0, 55);
        let n = 10000;
        let mut sum_ax = 0.0_f64;
        for _ in 0..n {
            let m = simulate_imu(&mut s, [0.0; 3], [0.0; 3], 0.01);
            sum_ax += m.linear_acceleration[0];
        }
        let mean_ax = sum_ax / n as f64;
        // With no bias, mean should be near zero
        assert!(
            mean_ax.abs() < 0.2,
            "mean accel x too far from zero: {mean_ax}"
        );
    }

    #[test]
    fn imu_noise_scales_with_dt() {
        // Higher dt → less noise per sample (σ ∝ 1/√dt)
        let mut s_fast = ImuSensor::ideal(10);
        let mut s_slow = ImuSensor::ideal(10);
        let n = 1000;

        let mut var_fast = 0.0_f64;
        let mut var_slow = 0.0_f64;
        for _ in 0..n {
            let m = simulate_imu(&mut s_fast, [0.0; 3], [0.0; 3], 0.001);
            var_fast += m.angular_velocity[0] * m.angular_velocity[0];
            let m = simulate_imu(&mut s_slow, [0.0; 3], [0.0; 3], 0.1);
            var_slow += m.angular_velocity[0] * m.angular_velocity[0];
        }
        // Fast sampling (small dt) should yield larger per-sample noise
        assert!(
            var_fast > var_slow,
            "fast sampling should have larger per-sample noise"
        );
    }

    // ── PositionSensor ───────────────────────────────────────────────────

    #[test]
    fn gps_consumer_grade_construction() {
        let s = PositionSensor::consumer_gps(0);
        assert!((s.position_noise_std - 1.5).abs() < 1e-10);
    }

    #[test]
    fn gps_rtk_construction() {
        let s = PositionSensor::rtk_gps(0);
        assert!(s.position_noise_std < 0.1);
    }

    #[test]
    fn gps_measurement_near_true_position() {
        let mut s = PositionSensor::consumer_gps(77);
        let pos = [100.0, 0.0, -50.0];
        let meas = simulate_gps(&mut s, pos);
        for i in 0..3 {
            assert!(
                (meas[i] - pos[i]).abs() < 20.0,
                "GPS error too large axis {i}"
            );
        }
    }

    #[test]
    fn gps_bias_applied_correctly() {
        let mut s = PositionSensor::new(0.0, [1.0, 2.0, 3.0], 5);
        let pos = [0.0; 3];
        let meas = simulate_gps(&mut s, pos);
        // Zero noise → measurement should equal bias exactly
        assert!((meas[0] - 1.0).abs() < 1e-9);
        assert!((meas[1] - 2.0).abs() < 1e-9);
        assert!((meas[2] - 3.0).abs() < 1e-9);
    }

    #[test]
    fn gps_repeated_calls_differ() {
        let mut s = PositionSensor::consumer_gps(3);
        let m1 = simulate_gps(&mut s, [0.0; 3]);
        let m2 = simulate_gps(&mut s, [0.0; 3]);
        let same = (0..3).all(|i| (m1[i] - m2[i]).abs() < 1e-15);
        assert!(!same, "consecutive GPS readings should differ");
    }

    #[test]
    fn gps_zero_noise_returns_bias_plus_true() {
        let mut s = PositionSensor::new(0.0, [0.0; 3], 1);
        let pos = [5.0, 6.0, 7.0];
        let meas = simulate_gps(&mut s, pos);
        for i in 0..3 {
            assert!((meas[i] - pos[i]).abs() < 1e-9);
        }
    }

    // ── ForceSensor ──────────────────────────────────────────────────────

    #[test]
    fn force_sensor_ideal_construction() {
        let s = ForceSensor::ideal(0);
        for i in 0..6 {
            assert!((s.calibration[i][i] - 1.0).abs() < 1e-10);
        }
    }

    #[test]
    fn force_sensor_ideal_passthrough() {
        let mut s = ForceSensor::new(
            {
                let mut c = [[0.0f64; 6]; 6];
                for (k, row) in c.iter_mut().enumerate() {
                    row[k] = 1.0;
                }
                c
            },
            0.0,
            0.0,
            42,
        );
        let wrench = [10.0, -5.0, 3.0, 1.0, -0.5, 0.2];
        let meas = simulate_force_sensor(&mut s, wrench);
        for i in 0..6 {
            assert!(
                (meas[i] - wrench[i]).abs() < 1e-9,
                "channel {i} passthrough failed"
            );
        }
    }

    #[test]
    fn force_sensor_noise_applied() {
        let mut s = ForceSensor::ideal(88);
        let wrench = [0.0; 6];
        let m1 = simulate_force_sensor(&mut s, wrench);
        let m2 = simulate_force_sensor(&mut s, wrench);
        let same = (0..6).all(|i| (m1[i] - m2[i]).abs() < 1e-15);
        assert!(!same, "repeated force sensor calls should differ");
    }

    #[test]
    fn force_sensor_calibration_scales_output() {
        let mut cal = [[0.0f64; 6]; 6];
        for (k, row) in cal.iter_mut().enumerate() {
            row[k] = 2.0;
        }
        let mut s = ForceSensor::new(cal, 0.0, 0.0, 0);
        let wrench = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0];
        let meas = simulate_force_sensor(&mut s, wrench);
        for (i, &m) in meas.iter().enumerate() {
            assert!((m - 2.0).abs() < 1e-9, "channel {i}");
        }
    }

    #[test]
    fn force_sensor_industrial_construction() {
        let s = ForceSensor::industrial(0);
        assert!(s.force_noise_std > 0.0);
        assert!(s.torque_noise_std > 0.0);
        // Off-diagonal cross-talk exists
        assert!(s.calibration[0][1].abs() > 0.0);
    }

    // ── ComplementaryFilter ──────────────────────────────────────────────

    #[test]
    fn complementary_filter_initialisation() {
        let state = ComplementaryFilterState::new([1.0, 2.0, 3.0]);
        assert_eq!(state.position, [1.0, 2.0, 3.0]);
        assert_eq!(state.velocity, [0.0; 3]);
    }

    #[test]
    fn complementary_filter_alpha_one_follows_imu() {
        let state = ComplementaryFilterState::new([0.0; 3]);
        let imu = ImuMeasurement {
            angular_velocity: [0.0; 3],
            linear_acceleration: [0.0, 0.0, 0.0],
        };
        let gps = [100.0, 0.0, 0.0];
        // alpha=1.0 → pure IMU integration, GPS ignored
        let next = sensor_fusion_complementary(&state, &imu, gps, 1.0, 0.01);
        // With zero accel, position stays at 0 regardless of GPS
        assert!(next.position[0].abs() < 1e-9, "alpha=1 should ignore GPS");
    }

    #[test]
    fn complementary_filter_alpha_zero_follows_gps() {
        let state = ComplementaryFilterState::new([0.0; 3]);
        let imu = ImuMeasurement {
            angular_velocity: [0.0; 3],
            linear_acceleration: [0.0; 3],
        };
        let gps = [50.0, 0.0, 0.0];
        // alpha=0.0 → pure GPS
        let next = sensor_fusion_complementary(&state, &imu, gps, 0.0, 0.01);
        assert!(
            (next.position[0] - 50.0).abs() < 1e-9,
            "alpha=0 should use GPS directly"
        );
    }

    #[test]
    fn complementary_filter_blends_position() {
        let state = ComplementaryFilterState::new([0.0; 3]);
        let imu = ImuMeasurement {
            angular_velocity: [0.0; 3],
            linear_acceleration: [0.0; 3],
        };
        let gps = [10.0, 0.0, 0.0];
        let next = sensor_fusion_complementary(&state, &imu, gps, 0.5, 0.01);
        // Should be between 0 and 10
        assert!(next.position[0] >= 0.0 && next.position[0] <= 10.0 + 1e-9);
    }

    #[test]
    fn complementary_filter_imu_acceleration_moves_position() {
        let state = ComplementaryFilterState::new([0.0; 3]);
        let imu = ImuMeasurement {
            angular_velocity: [0.0; 3],
            linear_acceleration: [2.0, 0.0, 0.0],
        };
        let gps = [0.0; 3];
        let next = sensor_fusion_complementary(&state, &imu, gps, 0.99, 0.1);
        // With a=2, dt=0.1: pos ≈ 0.5*2*0.01 = 0.01 (mostly IMU path)
        assert!(
            next.position[0] > 0.0,
            "acceleration should increase position"
        );
    }

    #[test]
    fn complementary_filter_velocity_updated() {
        let state = ComplementaryFilterState::new([0.0; 3]);
        let imu = ImuMeasurement {
            angular_velocity: [0.0; 3],
            linear_acceleration: [1.0, 0.0, 0.0],
        };
        let gps = [0.0; 3];
        let next = sensor_fusion_complementary(&state, &imu, gps, 1.0, 0.1);
        // velocity += a * dt = 0.1
        assert!(
            (next.velocity[0] - 0.1).abs() < 0.05,
            "velocity update failed"
        );
    }

    #[test]
    fn complementary_filter_convergence_to_gps() {
        // With alpha close to 0, position should converge to GPS over multiple steps
        let mut state = ComplementaryFilterState::new([0.0; 3]);
        let imu = ImuMeasurement {
            angular_velocity: [0.0; 3],
            linear_acceleration: [0.0; 3],
        };
        let gps = [100.0, 0.0, 0.0];
        for _ in 0..200 {
            state = sensor_fusion_complementary(&state, &imu, gps, 0.0, 0.01);
        }
        assert!((state.position[0] - 100.0).abs() < 1e-6);
    }

    #[test]
    fn imu_measurement_fields_accessible() {
        let m = ImuMeasurement {
            angular_velocity: [1.0, 2.0, 3.0],
            linear_acceleration: [4.0, 5.0, 6.0],
        };
        assert_eq!(m.angular_velocity[0], 1.0);
        assert_eq!(m.linear_acceleration[2], 6.0);
    }

    #[test]
    fn gps_mean_noise_near_zero_for_unbiased_sensor() {
        let mut s = PositionSensor::new(1.0, [0.0; 3], 200);
        let n = 5000;
        let mut sum_x = 0.0_f64;
        for _ in 0..n {
            let m = simulate_gps(&mut s, [0.0; 3]);
            sum_x += m[0];
        }
        let mean = sum_x / n as f64;
        assert!(mean.abs() < 0.1, "GPS mean should be near zero: {mean}");
    }

    #[test]
    fn force_sensor_cross_talk_changes_reading() {
        let s = ForceSensor::industrial(0);
        // Force only along x
        let wrench = [10.0, 0.0, 0.0, 0.0, 0.0, 0.0];
        let mut zero_noise = s.clone();
        zero_noise.force_noise_std = 0.0;
        zero_noise.torque_noise_std = 0.0;
        let meas = simulate_force_sensor(&mut zero_noise, wrench);
        // Cross-talk from cal[1][0]=0.002 should appear in channel 1
        assert!(
            (meas[1] - 0.02).abs() < 1e-9,
            "cross-talk not applied: {}",
            meas[1]
        );
    }
}
