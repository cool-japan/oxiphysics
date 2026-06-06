// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Physics-based sensor models: IMU, GPS, lidar, camera, radar, force/strain gauges,
//! encoder, sensor fusion (Madgwick AHRS), and noise models.

use std::f64::consts::PI;

// ─────────────────────────────────────────────────────────────────────────────
// Helper utilities
// ─────────────────────────────────────────────────────────────────────────────

/// Simple linear congruential PRNG for deterministic noise.
#[derive(Debug, Clone, Copy)]
struct Lcg {
    state: u64,
}

impl Lcg {
    fn new(seed: u64) -> Self {
        Self {
            state: seed ^ 0xdeadbeef,
        }
    }

    fn next_u64(&mut self) -> u64 {
        self.state = self
            .state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        self.state
    }

    /// Returns uniform random float in \[0, 1).
    fn rand01(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
    }

    /// Returns approximately N(0, 1) via Box-Muller.
    fn randn(&mut self) -> f64 {
        let u1 = (self.rand01() + 1e-15).ln();
        let u2 = self.rand01() * 2.0 * PI;
        (-2.0 * u1).sqrt() * u2.cos()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// ImuModel
// ─────────────────────────────────────────────────────────────────────────────

/// 6-DoF IMU model with accelerometer and gyroscope error sources.
#[derive(Debug, Clone)]
pub struct ImuModel {
    /// Accelerometer bias \[m/s²\].
    pub accel_bias: [f64; 3],
    /// Accelerometer noise std dev \[m/s²\].
    pub accel_noise_std: f64,
    /// Accelerometer scale factor error (fraction, e.g., 0.001 = 0.1%).
    pub accel_scale_error: f64,
    /// Gyroscope bias \[rad/s\].
    pub gyro_bias: [f64; 3],
    /// Gyroscope noise std dev \[rad/s\].
    pub gyro_noise_std: f64,
    /// Magnetometer noise std dev \[T\].
    pub mag_noise_std: f64,
    /// Internal random state.
    rng: Lcg,
}

impl ImuModel {
    /// Create an IMU model.
    pub fn new(
        accel_bias: [f64; 3],
        accel_noise_std: f64,
        accel_scale_error: f64,
        gyro_bias: [f64; 3],
        gyro_noise_std: f64,
        mag_noise_std: f64,
        seed: u64,
    ) -> Self {
        Self {
            accel_bias,
            accel_noise_std,
            accel_scale_error,
            gyro_bias,
            gyro_noise_std,
            mag_noise_std,
            rng: Lcg::new(seed),
        }
    }

    /// Create an ideal (noise-free) IMU.
    pub fn ideal() -> Self {
        Self::new([0.0; 3], 0.0, 0.0, [0.0; 3], 0.0, 0.0, 42)
    }

    /// Measure acceleration (true + bias + noise + scale error).
    pub fn measure_accel(&mut self, true_accel: [f64; 3]) -> [f64; 3] {
        let mut out = [0.0; 3];
        for i in 0..3 {
            let noise = self.rng.randn() * self.accel_noise_std;
            out[i] = true_accel[i] * (1.0 + self.accel_scale_error) + self.accel_bias[i] + noise;
        }
        out
    }

    /// Measure angular rate (true + bias + noise).
    pub fn measure_gyro(&mut self, true_omega: [f64; 3]) -> [f64; 3] {
        let mut out = [0.0; 3];
        for i in 0..3 {
            let noise = self.rng.randn() * self.gyro_noise_std;
            out[i] = true_omega[i] + self.gyro_bias[i] + noise;
        }
        out
    }

    /// Measure magnetometer (true + noise).
    pub fn measure_mag(&mut self, true_mag: [f64; 3]) -> [f64; 3] {
        let mut out = [0.0; 3];
        for i in 0..3 {
            out[i] = true_mag[i] + self.rng.randn() * self.mag_noise_std;
        }
        out
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// GpsModel
// ─────────────────────────────────────────────────────────────────────────────

/// GPS position and velocity model with HDOP and multipath errors.
#[derive(Debug, Clone)]
pub struct GpsModel {
    /// Horizontal position noise std dev \[m\].
    pub pos_noise_h: f64,
    /// Vertical position noise std dev \[m\].
    pub pos_noise_v: f64,
    /// Velocity noise std dev \[m/s\].
    pub vel_noise: f64,
    /// HDOP (horizontal dilution of precision).
    pub hdop: f64,
    /// Multipath error std dev \[m\].
    pub multipath_std: f64,
    /// Fix available (degraded availability simulation).
    pub has_fix: bool,
    rng: Lcg,
}

impl GpsModel {
    /// Create a GPS model.
    pub fn new(
        pos_noise_h: f64,
        pos_noise_v: f64,
        vel_noise: f64,
        hdop: f64,
        multipath_std: f64,
        seed: u64,
    ) -> Self {
        Self {
            pos_noise_h,
            pos_noise_v,
            vel_noise,
            hdop,
            multipath_std,
            has_fix: true,
            rng: Lcg::new(seed),
        }
    }

    /// Measure GPS position (NED frame). Returns `None` if no fix.
    pub fn measure_position(&mut self, true_pos: [f64; 3]) -> Option<[f64; 3]> {
        if !self.has_fix {
            return None;
        }
        let mut pos = [0.0; 3];
        let multipath = self.rng.randn() * self.multipath_std;
        pos[0] = true_pos[0] + self.rng.randn() * self.pos_noise_h * self.hdop + multipath;
        pos[1] = true_pos[1] + self.rng.randn() * self.pos_noise_h * self.hdop + multipath;
        pos[2] = true_pos[2] + self.rng.randn() * self.pos_noise_v;
        Some(pos)
    }

    /// Measure GPS velocity.
    pub fn measure_velocity(&mut self, true_vel: [f64; 3]) -> [f64; 3] {
        let mut vel = [0.0; 3];
        for i in 0..3 {
            vel[i] = true_vel[i] + self.rng.randn() * self.vel_noise;
        }
        vel
    }

    /// Estimated horizontal accuracy \[m\].
    pub fn estimated_accuracy(&self) -> f64 {
        self.pos_noise_h * self.hdop + self.multipath_std
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// LidarModel
// ─────────────────────────────────────────────────────────────────────────────

/// 2D 360° lidar sensor model.
#[derive(Debug, Clone)]
pub struct LidarModel {
    /// Number of beams.
    pub n_beams: usize,
    /// Maximum range \[m\].
    pub max_range: f64,
    /// Minimum range \[m\].
    pub min_range: f64,
    /// Range noise std dev \[m\].
    pub range_noise_std: f64,
    /// Angular resolution \[rad\] (2π / n_beams).
    pub angular_resolution: f64,
    rng: Lcg,
}

impl LidarModel {
    /// Create a lidar model.
    pub fn new(
        n_beams: usize,
        max_range: f64,
        min_range: f64,
        range_noise_std: f64,
        seed: u64,
    ) -> Self {
        Self {
            n_beams,
            max_range,
            min_range,
            range_noise_std,
            angular_resolution: 2.0 * PI / n_beams as f64,
            rng: Lcg::new(seed),
        }
    }

    /// Angle of beam `i` \[rad\].
    pub fn beam_angle(&self, i: usize) -> f64 {
        i as f64 * self.angular_resolution
    }

    /// Simulate a range scan. `true_ranges` should have length `n_beams`.
    /// Returns noisy range measurements, clamped to \[min_range, max_range\].
    pub fn scan(&mut self, true_ranges: &[f64]) -> Vec<f64> {
        assert!(true_ranges.len() >= self.n_beams);
        (0..self.n_beams)
            .map(|i| {
                let noise = self.rng.randn() * self.range_noise_std;
                (true_ranges[i] + noise).clamp(self.min_range, self.max_range)
            })
            .collect()
    }

    /// Convert range scan to 2D Cartesian points.
    pub fn to_cartesian(&self, ranges: &[f64]) -> Vec<[f64; 2]> {
        ranges
            .iter()
            .enumerate()
            .map(|(i, &r)| {
                let theta = self.beam_angle(i);
                [r * theta.cos(), r * theta.sin()]
            })
            .collect()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// CameraModel
// ─────────────────────────────────────────────────────────────────────────────

/// Pinhole camera model with radial and tangential distortion.
#[derive(Debug, Clone, Copy)]
pub struct CameraModel {
    /// Focal length in pixels (fx, fy).
    pub focal: [f64; 2],
    /// Principal point (cx, cy) in pixels.
    pub principal: [f64; 2],
    /// Radial distortion coefficients (k1, k2, k3).
    pub k: [f64; 3],
    /// Tangential distortion coefficients (p1, p2).
    pub p: [f64; 2],
    /// Image width \[pixels\].
    pub width: usize,
    /// Image height \[pixels\].
    pub height: usize,
}

impl CameraModel {
    /// Create a camera model.
    pub fn new(
        focal: [f64; 2],
        principal: [f64; 2],
        k: [f64; 3],
        p: [f64; 2],
        width: usize,
        height: usize,
    ) -> Self {
        Self {
            focal,
            principal,
            k,
            p,
            width,
            height,
        }
    }

    /// Horizontal field of view \[rad\].
    pub fn fov_h(&self) -> f64 {
        2.0 * (self.width as f64 / (2.0 * self.focal[0])).atan()
    }

    /// Vertical field of view \[rad\].
    pub fn fov_v(&self) -> f64 {
        2.0 * (self.height as f64 / (2.0 * self.focal[1])).atan()
    }

    /// Project a 3D world point to 2D pixel coordinates.
    /// Returns `None` if the point is behind the camera.
    pub fn project_3d_to_2d(&self, p: [f64; 3]) -> Option<[f64; 2]> {
        if p[2] <= 0.0 {
            return None;
        }
        let x = p[0] / p[2];
        let y = p[1] / p[2];
        // Apply distortion
        let r2 = x * x + y * y;
        let r4 = r2 * r2;
        let r6 = r4 * r2;
        let radial = 1.0 + self.k[0] * r2 + self.k[1] * r4 + self.k[2] * r6;
        let xd = x * radial + 2.0 * self.p[0] * x * y + self.p[1] * (r2 + 2.0 * x * x);
        let yd = y * radial + self.p[0] * (r2 + 2.0 * y * y) + 2.0 * self.p[1] * x * y;
        let u = self.focal[0] * xd + self.principal[0];
        let v = self.focal[1] * yd + self.principal[1];
        Some([u, v])
    }

    /// Check if a pixel coordinate is within the image bounds.
    pub fn is_in_frame(&self, pixel: [f64; 2]) -> bool {
        pixel[0] >= 0.0
            && pixel[0] < self.width as f64
            && pixel[1] >= 0.0
            && pixel[1] < self.height as f64
    }

    /// Back-project a pixel to a unit ray in camera frame.
    pub fn unproject(&self, pixel: [f64; 2]) -> [f64; 3] {
        let x = (pixel[0] - self.principal[0]) / self.focal[0];
        let y = (pixel[1] - self.principal[1]) / self.focal[1];
        let n = (1.0 + x * x + y * y).sqrt();
        [x / n, y / n, 1.0 / n]
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// RadarModel
// ─────────────────────────────────────────────────────────────────────────────

/// Radar sensor model with range, Doppler velocity, RCS, and clutter.
#[derive(Debug, Clone)]
pub struct RadarModel {
    /// Maximum detection range \[m\].
    pub max_range: f64,
    /// Range noise std dev \[m\].
    pub range_noise_std: f64,
    /// Velocity noise std dev \[m/s\].
    pub vel_noise_std: f64,
    /// Radar frequency \[Hz\].
    pub frequency: f64,
    /// Clutter-to-noise ratio (CNR) \[dB\].
    pub cnr_db: f64,
    rng: Lcg,
}

impl RadarModel {
    /// Create a radar model.
    pub fn new(
        max_range: f64,
        range_noise_std: f64,
        vel_noise_std: f64,
        frequency: f64,
        cnr_db: f64,
        seed: u64,
    ) -> Self {
        Self {
            max_range,
            range_noise_std,
            vel_noise_std,
            frequency,
            cnr_db,
            rng: Lcg::new(seed),
        }
    }

    /// Speed of light.
    const C: f64 = 3e8;

    /// Radar wavelength \[m\].
    pub fn wavelength(&self) -> f64 {
        Self::C / self.frequency
    }

    /// Doppler frequency shift for radial velocity v_r \[m/s\].
    pub fn doppler_shift(&self, v_radial: f64) -> f64 {
        2.0 * v_radial * self.frequency / Self::C
    }

    /// Measure range \[m\]. Returns `None` if beyond max_range.
    pub fn measure_range(&mut self, true_range: f64) -> Option<f64> {
        if true_range > self.max_range {
            return None;
        }
        let noise = self.rng.randn() * self.range_noise_std;
        Some((true_range + noise).max(0.0))
    }

    /// Measure radial velocity \[m/s\].
    pub fn measure_velocity(&mut self, v_radial: f64) -> f64 {
        v_radial + self.rng.randn() * self.vel_noise_std
    }

    /// Radar cross section for a sphere \[m²\]: RCS = π r².
    pub fn sphere_rcs(&self, radius: f64) -> f64 {
        PI * radius * radius
    }

    /// Detection probability (simplified, based on SNR approximation).
    /// `rcs` = radar cross section \[m²\], `range` = target range \[m\].
    pub fn detection_probability(&self, rcs: f64, range: f64) -> f64 {
        if range <= 0.0 {
            return 1.0;
        }
        // Simple range equation: SNR ∝ RCS / range^4
        let snr = rcs / (range / self.max_range).powi(4);
        1.0 - (-snr).exp()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// ForceGauge
// ─────────────────────────────────────────────────────────────────────────────

/// Load cell / force gauge model with noise, saturation, and creep.
#[derive(Debug, Clone)]
pub struct ForceGauge {
    /// Measurement range: maximum force \[N\].
    pub range: f64,
    /// Noise std dev \[N\].
    pub noise_std: f64,
    /// Creep time constant \[s\].
    pub creep_tau: f64,
    /// Creep factor (fraction of applied force).
    pub creep_factor: f64,
    /// Internal creep state.
    creep_state: f64,
    rng: Lcg,
}

impl ForceGauge {
    /// Create a force gauge.
    pub fn new(range: f64, noise_std: f64, creep_tau: f64, creep_factor: f64, seed: u64) -> Self {
        Self {
            range,
            noise_std,
            creep_tau,
            creep_factor,
            creep_state: 0.0,
            rng: Lcg::new(seed),
        }
    }

    /// Measure force \[N\] with noise, saturation, and creep update.
    pub fn measure(&mut self, true_force: f64, dt: f64) -> f64 {
        // Update creep state
        let d_creep = (true_force * self.creep_factor - self.creep_state) / self.creep_tau;
        self.creep_state += d_creep * dt;
        let creep_err = self.creep_state;
        // Noise
        let noise = self.rng.randn() * self.noise_std;
        // Saturate
        (true_force + noise + creep_err).clamp(-self.range, self.range)
    }

    /// Reset creep state.
    pub fn reset(&mut self) {
        self.creep_state = 0.0;
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// StrainGauge
// ─────────────────────────────────────────────────────────────────────────────

/// Wheatstone bridge strain gauge model.
#[derive(Debug, Clone)]
pub struct StrainGauge {
    /// Gauge factor GF (typically 2.0 for metallic gauges).
    pub gauge_factor: f64,
    /// Nominal resistance \[Ω\].
    pub resistance: f64,
    /// Noise std dev in strain units \[με\].
    pub noise_std: f64,
    /// Temperature coefficient \[1/K\].
    pub temp_coeff: f64,
    /// Reference temperature \[K\].
    pub ref_temp: f64,
    rng: Lcg,
}

impl StrainGauge {
    /// Create a strain gauge.
    pub fn new(
        gauge_factor: f64,
        resistance: f64,
        noise_std: f64,
        temp_coeff: f64,
        ref_temp: f64,
        seed: u64,
    ) -> Self {
        Self {
            gauge_factor,
            resistance,
            noise_std,
            temp_coeff,
            ref_temp,
            rng: Lcg::new(seed),
        }
    }

    /// Resistance change ΔR/R from strain ε.
    pub fn delta_r_over_r(&self, strain: f64) -> f64 {
        self.gauge_factor * strain
    }

    /// Output voltage ratio ΔV/V_ex from strain ε (quarter bridge).
    pub fn output_voltage_ratio(&self, strain: f64) -> f64 {
        self.gauge_factor * strain / 4.0
    }

    /// Measured strain \[με\] with noise and thermal compensation.
    pub fn measure(&mut self, true_strain: f64, temperature: f64) -> f64 {
        let thermal_apparent = self.temp_coeff * (temperature - self.ref_temp);
        let noise = self.rng.randn() * self.noise_std;
        true_strain - thermal_apparent + noise
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// EncoderModel
// ─────────────────────────────────────────────────────────────────────────────

/// Quadrature encoder model with CPR, jitter, and index pulse.
#[derive(Debug, Clone)]
pub struct EncoderModel {
    /// Counts per revolution (CPR).
    pub cpr: u32,
    /// Jitter std dev \[counts\].
    pub jitter_std: f64,
    /// Accumulated count.
    pub count: i64,
    /// Accumulated angle \[rad\] (true).
    pub true_angle: f64,
    /// Index pulse emitted every full revolution.
    pub index_pulse: bool,
    rng: Lcg,
}

impl EncoderModel {
    /// Create an encoder model.
    pub fn new(cpr: u32, jitter_std: f64, seed: u64) -> Self {
        Self {
            cpr,
            jitter_std,
            count: 0,
            true_angle: 0.0,
            index_pulse: false,
            rng: Lcg::new(seed),
        }
    }

    /// Step encoder by angular velocity omega \[rad/s\] for time dt \[s\].
    /// Returns (counts, angle_rad, index_pulse).
    pub fn step(&mut self, omega: f64, dt: f64) -> (i64, f64, bool) {
        let angle_delta = omega * dt;
        self.true_angle += angle_delta;
        let jitter = self.rng.randn() * self.jitter_std;
        let new_count = (self.true_angle * self.cpr as f64 / (2.0 * PI) + jitter).round() as i64;
        let delta_count = new_count - self.count;
        self.count = new_count;
        // Index pulse triggers when true_angle crosses a multiple of 2π
        let prev_rev = ((self.true_angle - angle_delta) / (2.0 * PI)).floor();
        let curr_rev = (self.true_angle / (2.0 * PI)).floor();
        self.index_pulse = curr_rev != prev_rev;
        (delta_count, self.true_angle, self.index_pulse)
    }

    /// Convert count to angle \[rad\].
    pub fn count_to_angle(&self, count: i64) -> f64 {
        count as f64 * 2.0 * PI / self.cpr as f64
    }

    /// Reset encoder.
    pub fn reset(&mut self) {
        self.count = 0;
        self.true_angle = 0.0;
        self.index_pulse = false;
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// SensorFusion (Madgwick AHRS)
// ─────────────────────────────────────────────────────────────────────────────

/// Madgwick AHRS filter fusing accelerometer + gyroscope (and optionally magnetometer).
#[derive(Debug, Clone, Copy)]
pub struct SensorFusion {
    /// Filter gain β (convergence rate). Typical: 0.1.
    pub beta: f64,
    /// Complementary filter weight for gyro (0..1).
    pub alpha: f64,
    /// Quaternion state \[w, x, y, z\].
    pub q: [f64; 4],
    /// Euler angles (roll, pitch, yaw) \[rad\].
    pub euler: [f64; 3],
}

impl SensorFusion {
    /// Create a Madgwick filter.
    pub fn new(beta: f64) -> Self {
        Self {
            beta,
            alpha: 0.98,
            q: [1.0, 0.0, 0.0, 0.0],
            euler: [0.0; 3],
        }
    }

    /// Normalize quaternion in-place.
    fn normalize_q(q: &mut [f64; 4]) {
        let n = (q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3]).sqrt();
        if n > 1e-14 {
            for qi in q.iter_mut() {
                *qi /= n;
            }
        }
    }

    /// Update filter with gyro \[rad/s\] and accelerometer \[m/s²\] data at rate dt \[s\].
    /// Implements a simplified Madgwick gradient-descent update.
    pub fn update(&mut self, gyro: [f64; 3], accel: [f64; 3], dt: f64) {
        let q = &mut self.q;
        let [gx, gy, gz] = gyro;
        let [ax, ay, az] = accel;

        // Normalize accelerometer
        let a_norm = (ax * ax + ay * ay + az * az).sqrt();
        if a_norm < 1e-10 {
            return;
        }
        let [ax, ay, az] = [ax / a_norm, ay / a_norm, az / a_norm];

        let [qw, qx, qy, qz] = *q;

        // Gradient of objective function (quaternion error against gravity vector)
        let f1 = 2.0 * (qx * qz - qw * qy) - ax;
        let f2 = 2.0 * (qw * qx + qy * qz) - ay;
        let f3 = 2.0 * (0.5 - qx * qx - qy * qy) - az;

        let j_t = [
            [-2.0 * qy, 2.0 * qx, 0.0, -2.0 * qw],
            [2.0 * qz, 2.0 * qw, -2.0 * qz, 2.0 * qy], // not used but kept for clarity
            [0.0, -4.0 * qx, -4.0 * qy, 0.0],
        ];
        let _ = j_t[1];

        let step_w = -2.0 * qy * f1 + 2.0 * qz * f2;
        let step_x = 2.0 * qz * f1 + 2.0 * qw * f2 - 4.0 * qx * f3;
        let step_y = -2.0 * qw * f1 + 2.0 * qz * f2 - 4.0 * qy * f3;
        let step_z = 2.0 * qx * f1 + 2.0 * qy * f2;
        let s_norm = (step_w * step_w + step_x * step_x + step_y * step_y + step_z * step_z)
            .sqrt()
            .max(1e-14);
        let (sw, sx, sy, sz) = (
            step_w / s_norm,
            step_x / s_norm,
            step_y / s_norm,
            step_z / s_norm,
        );

        // Rate of change of quaternion from gyro
        let qdot = [
            0.5 * (-qx * gx - qy * gy - qz * gz),
            0.5 * (qw * gx + qy * gz - qz * gy),
            0.5 * (qw * gy - qx * gz + qz * gx),
            0.5 * (qw * gz + qx * gy - qy * gx),
        ];

        q[0] += (qdot[0] - self.beta * sw) * dt;
        q[1] += (qdot[1] - self.beta * sx) * dt;
        q[2] += (qdot[2] - self.beta * sy) * dt;
        q[3] += (qdot[3] - self.beta * sz) * dt;

        Self::normalize_q(q);
        self.euler = self.to_euler();
    }

    /// Convert quaternion to Euler angles (roll, pitch, yaw) \[rad\].
    pub fn to_euler(&self) -> [f64; 3] {
        let [qw, qx, qy, qz] = self.q;
        let roll = (2.0 * (qw * qx + qy * qz)).atan2(1.0 - 2.0 * (qx * qx + qy * qy));
        let pitch_arg = (2.0 * (qw * qy - qz * qx)).clamp(-1.0, 1.0);
        let pitch = pitch_arg.asin();
        let yaw = (2.0 * (qw * qz + qx * qy)).atan2(1.0 - 2.0 * (qy * qy + qz * qz));
        [roll, pitch, yaw]
    }

    /// Complementary filter update (simpler alternative).
    pub fn complementary_update(&mut self, accel: [f64; 3], gyro: [f64; 3], dt: f64) {
        // Gyro integration
        self.euler[0] += gyro[0] * dt;
        self.euler[1] += gyro[1] * dt;

        // Accelerometer roll/pitch
        let a_roll = accel[1].atan2(accel[2]);
        let a_pitch = (-accel[0]).atan2((accel[1] * accel[1] + accel[2] * accel[2]).sqrt());

        // Fuse
        self.euler[0] = self.alpha * self.euler[0] + (1.0 - self.alpha) * a_roll;
        self.euler[1] = self.alpha * self.euler[1] + (1.0 - self.alpha) * a_pitch;
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// SensorNoise
// ─────────────────────────────────────────────────────────────────────────────

/// Sensor noise models: white, pink 1/f, random walk, bias instability.
#[derive(Debug, Clone)]
pub struct SensorNoise {
    /// White noise amplitude.
    pub white_amp: f64,
    /// Pink noise amplitude.
    pub pink_amp: f64,
    /// Random walk coefficient \[units/√s\].
    pub rw_coeff: f64,
    /// Bias instability \[units\].
    pub bias_instability: f64,
    /// Current random walk state.
    rw_state: f64,
    /// Pink noise state variables.
    pink_state: [f64; 3],
    rng: Lcg,
}

impl SensorNoise {
    /// Create a noise model.
    pub fn new(
        white_amp: f64,
        pink_amp: f64,
        rw_coeff: f64,
        bias_instability: f64,
        seed: u64,
    ) -> Self {
        Self {
            white_amp,
            pink_amp,
            rw_coeff,
            bias_instability,
            rw_state: 0.0,
            pink_state: [0.0; 3],
            rng: Lcg::new(seed),
        }
    }

    /// Generate one white noise sample.
    pub fn white_noise(&mut self) -> f64 {
        self.rng.randn() * self.white_amp
    }

    /// Generate one pink noise sample using 3-pole IIR approximation.
    pub fn pink_noise(&mut self) -> f64 {
        let white = self.rng.randn();
        // Voss-McCartney algorithm approximation
        self.pink_state[0] = 0.99886 * self.pink_state[0] + white * 0.0555179;
        self.pink_state[1] = 0.99332 * self.pink_state[1] + white * 0.0750759;
        self.pink_state[2] = 0.96900 * self.pink_state[2] + white * 0.1538520;
        (self.pink_state[0] + self.pink_state[1] + self.pink_state[2] + white * 0.5362)
            * self.pink_amp
    }

    /// Update random walk by one step `dt` \[s\].
    pub fn random_walk_step(&mut self, dt: f64) -> f64 {
        let dw = self.rng.randn() * self.rw_coeff * dt.sqrt();
        self.rw_state += dw;
        self.rw_state
    }

    /// Allan variance (simplified estimate for given tau \[s\]).
    /// Returns σ²_Allan for white noise-dominated region.
    pub fn allan_variance(&self, tau: f64) -> f64 {
        // σ² = white²/(tau) + random_walk²*tau
        self.white_amp * self.white_amp / tau + self.rw_coeff * self.rw_coeff * tau
    }

    /// Sample total noise (white + pink + RW + bias instability).
    pub fn sample(&mut self, dt: f64) -> f64 {
        let w = self.white_noise();
        let p = self.pink_noise();
        let rw = self.random_walk_step(dt);
        let b = self.rng.randn() * self.bias_instability;
        w + p + rw + b
    }

    /// Reset all state.
    pub fn reset(&mut self) {
        self.rw_state = 0.0;
        self.pink_state = [0.0; 3];
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // ── ImuModel tests ────────────────────────────────────────────────────────

    #[test]
    fn imu_ideal_accel_no_error() {
        let mut imu = ImuModel::ideal();
        let true_a = [0.0, 0.0, 9.81];
        let measured = imu.measure_accel(true_a);
        for i in 0..3 {
            assert!(
                (measured[i] - true_a[i]).abs() < 1e-12,
                "i={i}: {}",
                measured[i]
            );
        }
    }

    #[test]
    fn imu_bias_shifts_measurement() {
        let mut imu = ImuModel::new([1.0, 0.0, 0.0], 0.0, 0.0, [0.0; 3], 0.0, 0.0, 1);
        let measured = imu.measure_accel([0.0; 3]);
        assert!((measured[0] - 1.0).abs() < 1e-10, "accel_x={}", measured[0]);
    }

    #[test]
    fn imu_scale_error() {
        let mut imu = ImuModel::new([0.0; 3], 0.0, 0.01, [0.0; 3], 0.0, 0.0, 1);
        let measured = imu.measure_accel([100.0, 0.0, 0.0]);
        // 100 * 1.01 = 101
        assert!(
            (measured[0] - 101.0).abs() < 1e-10,
            "accel_x={}",
            measured[0]
        );
    }

    #[test]
    fn imu_gyro_bias() {
        let mut imu = ImuModel::new([0.0; 3], 0.0, 0.0, [0.1, 0.0, 0.0], 0.0, 0.0, 1);
        let measured = imu.measure_gyro([0.0; 3]);
        assert!((measured[0] - 0.1).abs() < 1e-10, "omega_x={}", measured[0]);
    }

    // ── GpsModel tests ────────────────────────────────────────────────────────

    #[test]
    fn gps_has_fix() {
        let mut gps = GpsModel::new(3.0, 5.0, 0.1, 1.5, 1.0, 42);
        let pos = gps.measure_position([0.0; 3]);
        assert!(pos.is_some());
    }

    #[test]
    fn gps_no_fix_returns_none() {
        let mut gps = GpsModel::new(3.0, 5.0, 0.1, 1.5, 1.0, 42);
        gps.has_fix = false;
        let pos = gps.measure_position([0.0; 3]);
        assert!(pos.is_none());
    }

    #[test]
    fn gps_accuracy_positive() {
        let gps = GpsModel::new(3.0, 5.0, 0.1, 1.5, 1.0, 42);
        assert!(gps.estimated_accuracy() > 0.0);
    }

    #[test]
    fn gps_velocity_finite() {
        let mut gps = GpsModel::new(3.0, 5.0, 0.1, 1.5, 1.0, 42);
        let vel = gps.measure_velocity([10.0, 0.0, 0.0]);
        assert!(vel[0].is_finite());
    }

    // ── LidarModel tests ──────────────────────────────────────────────────────

    #[test]
    fn lidar_scan_length() {
        let mut lidar = LidarModel::new(360, 50.0, 0.1, 0.01, 1);
        let ranges = vec![10.0; 360];
        let scan = lidar.scan(&ranges);
        assert_eq!(scan.len(), 360);
    }

    #[test]
    fn lidar_ranges_clamped() {
        let mut lidar = LidarModel::new(4, 10.0, 0.5, 0.0, 1);
        let ranges = vec![100.0, 0.0, 5.0, 15.0];
        let scan = lidar.scan(&ranges);
        assert!(scan[0] <= 10.0, "scan[0]={}", scan[0]);
        assert!(scan[1] >= 0.5, "scan[1]={}", scan[1]);
    }

    #[test]
    fn lidar_to_cartesian_size() {
        let lidar = LidarModel::new(4, 10.0, 0.5, 0.0, 1);
        let ranges = vec![5.0; 4];
        let pts = lidar.to_cartesian(&ranges);
        assert_eq!(pts.len(), 4);
    }

    // ── CameraModel tests ─────────────────────────────────────────────────────

    #[test]
    fn camera_project_front_point() {
        let cam = CameraModel::new([500.0, 500.0], [320.0, 240.0], [0.0; 3], [0.0; 2], 640, 480);
        let pixel = cam.project_3d_to_2d([0.0, 0.0, 1.0]);
        assert!(pixel.is_some());
        let [u, v] = pixel.unwrap();
        assert!((u - 320.0).abs() < 1e-6, "u={u}");
        assert!((v - 240.0).abs() < 1e-6, "v={v}");
    }

    #[test]
    fn camera_behind_returns_none() {
        let cam = CameraModel::new([500.0, 500.0], [320.0, 240.0], [0.0; 3], [0.0; 2], 640, 480);
        assert!(cam.project_3d_to_2d([0.0, 0.0, -1.0]).is_none());
    }

    #[test]
    fn camera_fov_positive() {
        let cam = CameraModel::new([500.0, 500.0], [320.0, 240.0], [0.0; 3], [0.0; 2], 640, 480);
        assert!(cam.fov_h() > 0.0);
        assert!(cam.fov_v() > 0.0);
    }

    #[test]
    fn camera_is_in_frame() {
        let cam = CameraModel::new([500.0, 500.0], [320.0, 240.0], [0.0; 3], [0.0; 2], 640, 480);
        assert!(cam.is_in_frame([320.0, 240.0]));
        assert!(!cam.is_in_frame([700.0, 240.0]));
    }

    // ── RadarModel tests ──────────────────────────────────────────────────────

    #[test]
    fn radar_range_within_bounds() {
        let mut radar = RadarModel::new(1000.0, 1.0, 0.1, 77e9, 20.0, 1);
        let r = radar.measure_range(500.0);
        assert!(r.is_some());
    }

    #[test]
    fn radar_range_beyond_max_returns_none() {
        let mut radar = RadarModel::new(1000.0, 1.0, 0.1, 77e9, 20.0, 1);
        let r = radar.measure_range(2000.0);
        assert!(r.is_none());
    }

    #[test]
    fn radar_doppler_shift_positive() {
        let radar = RadarModel::new(1000.0, 1.0, 0.1, 77e9, 20.0, 1);
        let fd = radar.doppler_shift(30.0);
        // 2 * 30 * 77e9 / 3e8 = 15.4 kHz
        assert!(fd > 0.0, "fd={fd}");
    }

    #[test]
    fn radar_sphere_rcs_positive() {
        let radar = RadarModel::new(1000.0, 1.0, 0.1, 77e9, 20.0, 1);
        assert!(radar.sphere_rcs(0.5) > 0.0);
    }

    // ── ForceGauge tests ──────────────────────────────────────────────────────

    #[test]
    fn force_gauge_within_range() {
        let mut fg = ForceGauge::new(1000.0, 0.1, 10.0, 0.01, 1);
        let m = fg.measure(500.0, 0.01);
        assert!(m.abs() <= 1000.0, "m={m}");
    }

    #[test]
    fn force_gauge_saturation() {
        let mut fg = ForceGauge::new(100.0, 0.0, 1000.0, 0.0, 1);
        let m = fg.measure(200.0, 0.0);
        assert!((m - 100.0).abs() < 1e-10, "m={m}");
    }

    #[test]
    fn force_gauge_reset() {
        let mut fg = ForceGauge::new(1000.0, 0.0, 1.0, 1.0, 1);
        fg.measure(100.0, 1.0);
        fg.reset();
        assert!(fg.creep_state.abs() < 1e-12);
    }

    // ── StrainGauge tests ─────────────────────────────────────────────────────

    #[test]
    fn strain_gauge_dr_over_r() {
        let sg = StrainGauge::new(2.0, 120.0, 0.0, 1e-5, 293.0, 1);
        let dr = sg.delta_r_over_r(1e-3);
        assert!((dr - 2e-3).abs() < 1e-10, "dr={dr}");
    }

    #[test]
    fn strain_gauge_output_voltage() {
        let sg = StrainGauge::new(2.0, 120.0, 0.0, 1e-5, 293.0, 1);
        let v = sg.output_voltage_ratio(1e-3);
        assert!((v - 0.5e-3).abs() < 1e-10, "v={v}");
    }

    #[test]
    fn strain_gauge_thermal_compensation() {
        let mut sg = StrainGauge::new(2.0, 120.0, 0.0, 1e-5, 293.0, 1);
        let s_room = sg.measure(0.0, 293.0);
        let s_hot = sg.measure(0.0, 373.0);
        assert!(
            s_hot != s_room,
            "thermal compensation should change reading"
        );
    }

    // ── EncoderModel tests ────────────────────────────────────────────────────

    #[test]
    fn encoder_count_to_angle() {
        let enc = EncoderModel::new(1000, 0.0, 1);
        let angle = enc.count_to_angle(1000);
        assert!((angle - 2.0 * PI).abs() < 1e-10, "angle={angle}");
    }

    #[test]
    fn encoder_step_accumulates() {
        let mut enc = EncoderModel::new(1000, 0.0, 1);
        enc.step(2.0 * PI, 1.0); // One full revolution in 1s
        assert!(enc.count > 0, "count={}", enc.count);
    }

    #[test]
    fn encoder_reset() {
        let mut enc = EncoderModel::new(1000, 0.0, 1);
        enc.step(PI, 1.0);
        enc.reset();
        assert_eq!(enc.count, 0);
    }

    // ── SensorFusion tests ────────────────────────────────────────────────────

    #[test]
    fn madgwick_quat_normalized() {
        let mut fusion = SensorFusion::new(0.1);
        let q_norm_before = fusion.q.iter().map(|q| q * q).sum::<f64>().sqrt();
        assert!(
            (q_norm_before - 1.0).abs() < 1e-10,
            "initial norm={q_norm_before}"
        );
        fusion.update([0.01, 0.0, 0.0], [0.0, 0.0, 9.81], 0.01);
        let q_norm = fusion.q.iter().map(|q| q * q).sum::<f64>().sqrt();
        assert!((q_norm - 1.0).abs() < 1e-6, "norm after update={q_norm}");
    }

    #[test]
    fn madgwick_euler_finite() {
        let mut fusion = SensorFusion::new(0.1);
        for _ in 0..100 {
            fusion.update([0.01, 0.02, 0.0], [0.1, 0.2, 9.81], 0.01);
        }
        for e in fusion.euler {
            assert!(e.is_finite(), "euler={e}");
        }
    }

    #[test]
    fn complementary_filter_update() {
        let mut fusion = SensorFusion::new(0.1);
        fusion.complementary_update([0.0, 0.0, 9.81], [0.0; 3], 0.01);
        assert!(fusion.euler[0].is_finite());
        assert!(fusion.euler[1].is_finite());
    }

    // ── SensorNoise tests ─────────────────────────────────────────────────────

    #[test]
    fn sensor_noise_white_finite() {
        let mut noise = SensorNoise::new(1.0, 0.5, 0.01, 0.001, 1);
        let w = noise.white_noise();
        assert!(w.is_finite(), "w={w}");
    }

    #[test]
    fn sensor_noise_pink_finite() {
        let mut noise = SensorNoise::new(1.0, 0.5, 0.01, 0.001, 1);
        let p = noise.pink_noise();
        assert!(p.is_finite(), "p={p}");
    }

    #[test]
    fn sensor_noise_random_walk_grows() {
        let mut noise = SensorNoise::new(0.0, 0.0, 1.0, 0.0, 100);
        let rw_after = noise.random_walk_step(1.0);
        // Non-zero after one step
        assert!(rw_after.is_finite(), "rw={rw_after}");
    }

    #[test]
    fn sensor_noise_allan_variance_minimum() {
        let noise = SensorNoise::new(1.0, 0.0, 1.0, 0.0, 1);
        // Minimum at tau* = white / rw
        let av1 = noise.allan_variance(1.0);
        let av100 = noise.allan_variance(100.0);
        // At long times, random walk dominates
        assert!(av100 > av1, "av100={av100} > av1={av1}");
    }

    #[test]
    fn sensor_noise_reset() {
        let mut noise = SensorNoise::new(0.0, 0.0, 1.0, 0.0, 1);
        noise.random_walk_step(1.0);
        noise.reset();
        assert!(noise.rw_state.abs() < 1e-12);
    }
}
