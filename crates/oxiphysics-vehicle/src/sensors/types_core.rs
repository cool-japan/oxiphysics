//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
use rand::RngExt;
use std::f64::consts::PI;

/// Pinhole camera model.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct PinholeCamera {
    /// Focal length in pixels along the x axis.
    pub fx: f64,
    /// Focal length in pixels along the y axis.
    pub fy: f64,
    /// Principal point x coordinate in pixels.
    pub cx: f64,
    /// Principal point y coordinate in pixels.
    pub cy: f64,
    /// Image width in pixels.
    pub width: u32,
    /// Image height in pixels.
    pub height: u32,
}
impl PinholeCamera {
    /// Construct from a horizontal field-of-view angle and image dimensions.
    ///
    /// Assumes square pixels and centred principal point.
    #[allow(dead_code)]
    pub fn new(fov_deg: f64, width: u32, height: u32) -> Self {
        let fx = (width as f64) / (2.0 * (fov_deg.to_radians() / 2.0).tan());
        let fy = fx;
        let cx = width as f64 / 2.0;
        let cy = height as f64 / 2.0;
        Self {
            fx,
            fy,
            cx,
            cy,
            width,
            height,
        }
    }
    /// Project a 3-D point (in camera frame, Z forward) to pixel coordinates.
    ///
    /// Returns `None` if the point is behind the camera or outside the image.
    #[allow(dead_code)]
    pub fn project(&self, point_3d: [f64; 3]) -> Option<[f64; 2]> {
        let (x, y, z) = (point_3d[0], point_3d[1], point_3d[2]);
        if z <= 0.0 {
            return None;
        }
        let u = self.fx * x / z + self.cx;
        let v = self.fy * y / z + self.cy;
        if u < 0.0 || u >= self.width as f64 || v < 0.0 || v >= self.height as f64 {
            return None;
        }
        Some([u, v])
    }
    /// Back-project a pixel plus depth to a 3-D point in camera frame.
    #[allow(dead_code)]
    pub fn unproject(&self, pixel: [f64; 2], depth: f64) -> [f64; 3] {
        let x = (pixel[0] - self.cx) * depth / self.fx;
        let y = (pixel[1] - self.cy) * depth / self.fy;
        [x, y, depth]
    }
}
/// A steering angle sensor measurement.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct SteeringAngleSensor {
    /// 1-σ noise in radians.
    pub noise_rad: f64,
    /// Constant bias in radians.
    pub bias_rad: f64,
}
impl SteeringAngleSensor {
    /// Create a high-resolution optical encoder model.
    #[allow(dead_code)]
    pub fn optical_encoder() -> Self {
        Self {
            noise_rad: 0.001,
            bias_rad: 0.0,
        }
    }
    /// Measure the steering angle.
    ///
    /// Returns the noisy steering angle in radians.
    #[allow(dead_code)]
    pub fn measure(&self, true_angle_rad: f64, noise_sample: f64) -> f64 {
        true_angle_rad + self.bias_rad + self.noise_rad * noise_sample
    }
}
/// Estimates front and rear tyre slip angles from vehicle kinematic state.
///
/// Uses the bicycle model:
/// ```text
/// alpha_f = delta - atan((vy + a * yaw_rate) / vx)
/// alpha_r =       - atan((vy - b * yaw_rate) / vx)
/// ```
/// where `a` is the distance from CoM to front axle and `b` to the rear axle.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct SlipAngleEstimator {
    /// Distance from CoM to front axle (m).
    pub a: f64,
    /// Distance from CoM to rear axle (m).
    pub b: f64,
}
impl SlipAngleEstimator {
    /// Create a new estimator from axle distances.
    #[allow(dead_code)]
    pub fn new(a: f64, b: f64) -> Self {
        Self { a, b }
    }
    /// Estimate slip angles.
    ///
    /// # Arguments
    /// * `vx`        – longitudinal velocity (m/s, must be non-zero)
    /// * `vy`        – lateral velocity (m/s)
    /// * `yaw_rate`  – yaw rate (rad/s)
    /// * `delta`     – front steering angle (rad)
    ///
    /// Returns `(alpha_front_rad, alpha_rear_rad)`.
    #[allow(dead_code)]
    pub fn estimate(&self, vx: f64, vy: f64, yaw_rate: f64, delta: f64) -> (f64, f64) {
        if vx.abs() < 0.5 {
            return (0.0, 0.0);
        }
        let alpha_f = delta - (vy + self.a * yaw_rate).atan2(vx);
        let alpha_r = -(vy - self.b * yaw_rate).atan2(vx);
        (alpha_f, alpha_r)
    }
    /// Compute front tyre lateral force coefficient from slip angle using the
    /// linearised cornering stiffness model.
    ///
    /// `Fy = C_alpha * alpha`  (N, positive = left)
    #[allow(dead_code)]
    pub fn lateral_force(&self, c_alpha: f64, slip_angle_rad: f64) -> f64 {
        c_alpha * slip_angle_rad
    }
}
/// A single LiDAR return point.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct LidarPoint {
    /// X coordinate in sensor frame (metres).
    pub x: f64,
    /// Y coordinate in sensor frame (metres).
    pub y: f64,
    /// Z coordinate in sensor frame (metres).
    pub z: f64,
    /// Intensity / reflectance (0–1 normalised).
    pub intensity: f64,
    /// Ring (beam) index.
    pub ring: usize,
}
/// IMU sensor model with noise parameters and constant bias.
///
/// (Note: the existing `ImuSensor` is a stateless helper; this struct holds
/// configuration state.)
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct ImuUnit {
    /// 1-σ accelerometer noise (m/s²).
    pub noise_std_accel: f64,
    /// 1-σ gyroscope noise (rad/s).
    pub noise_std_gyro: f64,
    /// Constant accelerometer bias `[bx, by, bz]` (m/s²).
    pub bias: [f64; 3],
}
#[allow(dead_code)]
impl ImuUnit {
    /// Create a new IMU unit.
    pub fn new(noise_std_accel: f64, noise_std_gyro: f64, bias: [f64; 3]) -> Self {
        Self {
            noise_std_accel,
            noise_std_gyro,
            bias,
        }
    }
    fn box_muller(rng: &mut impl rand::Rng) -> f64 {
        let u1: f64 = rng.random::<f64>().max(1e-15);
        let u2: f64 = rng.random::<f64>();
        (-2.0 * u1.ln()).sqrt() * (2.0 * PI * u2).cos()
    }
    /// Measure acceleration with noise and bias.
    pub fn measure_accel(&self, true_accel: [f64; 3], rng: &mut impl rand::Rng) -> [f64; 3] {
        [
            true_accel[0] + self.bias[0] + self.noise_std_accel * Self::box_muller(rng),
            true_accel[1] + self.bias[1] + self.noise_std_accel * Self::box_muller(rng),
            true_accel[2] + self.bias[2] + self.noise_std_accel * Self::box_muller(rng),
        ]
    }
    /// Measure angular velocity with noise.
    pub fn measure_gyro(&self, true_gyro: [f64; 3], rng: &mut impl rand::Rng) -> [f64; 3] {
        [
            true_gyro[0] + self.noise_std_gyro * Self::box_muller(rng),
            true_gyro[1] + self.noise_std_gyro * Self::box_muller(rng),
            true_gyro[2] + self.noise_std_gyro * Self::box_muller(rng),
        ]
    }
}
/// A scalar Kalman filter tracking position with a constant-velocity model.
///
/// State vector: `[position, velocity]`
/// Process model: position += velocity * dt
/// Observation: position only
///
/// This provides a principled fusion of a noisy position sensor (e.g. GPS)
/// with dead-reckoning (integrated velocity).
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct KalmanFilter1D {
    /// Estimated position (m).
    pub position: f64,
    /// Estimated velocity (m/s).
    pub velocity: f64,
    /// Position variance.
    pub p_pos: f64,
    /// Velocity variance.
    pub p_vel: f64,
    /// Covariance between position and velocity.
    pub p_pv: f64,
    /// Process noise (position variance per second).
    pub q_pos: f64,
    /// Process noise (velocity variance per second).
    pub q_vel: f64,
    /// Observation noise variance (m²).
    pub r_obs: f64,
}
impl KalmanFilter1D {
    /// Create a new filter at position 0, velocity 0.
    ///
    /// * `q_pos` – process noise for position (m² s⁻¹)
    /// * `q_vel` – process noise for velocity (m² s⁻¹)
    /// * `r_obs` – observation noise variance (m²)
    #[allow(dead_code)]
    pub fn new(q_pos: f64, q_vel: f64, r_obs: f64) -> Self {
        Self {
            position: 0.0,
            velocity: 0.0,
            p_pos: 1.0,
            p_vel: 1.0,
            p_pv: 0.0,
            q_pos,
            q_vel,
            r_obs,
        }
    }
    /// Predict step: propagate state forward by `dt` seconds.
    #[allow(dead_code)]
    pub fn predict(&mut self, dt: f64) {
        self.position += self.velocity * dt;
        let p_pos_new =
            self.p_pos + dt * (self.p_pv + self.p_pv) + dt * dt * self.p_vel + self.q_pos * dt;
        let p_vel_new = self.p_vel + self.q_vel * dt;
        let p_pv_new = self.p_pv + dt * self.p_vel;
        self.p_pos = p_pos_new;
        self.p_vel = p_vel_new;
        self.p_pv = p_pv_new;
    }
    /// Update step: incorporate a position measurement `z`.
    ///
    /// Uses an H = \[1, 0\] observation matrix (position-only sensor).
    #[allow(dead_code)]
    pub fn update(&mut self, z: f64) {
        let s = self.p_pos + self.r_obs;
        if s.abs() < 1e-30 {
            return;
        }
        let k_pos = self.p_pos / s;
        let k_vel = self.p_pv / s;
        let innovation = z - self.position;
        self.position += k_pos * innovation;
        self.velocity += k_vel * innovation;
        let p_pos_new = (1.0 - k_pos) * self.p_pos;
        let p_vel_new = self.p_vel - k_vel * self.p_pv;
        let p_pv_new = (1.0 - k_pos) * self.p_pv;
        self.p_pos = p_pos_new.max(0.0);
        self.p_vel = p_vel_new.max(0.0);
        self.p_pv = p_pv_new;
    }
    /// Position standard deviation (m).
    #[allow(dead_code)]
    pub fn position_std(&self) -> f64 {
        self.p_pos.max(0.0).sqrt()
    }
    /// Velocity standard deviation (m/s).
    #[allow(dead_code)]
    pub fn velocity_std(&self) -> f64 {
        self.p_vel.max(0.0).sqrt()
    }
}
/// Simple GPS position sensor with isotropic Gaussian noise.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct GpsUnit {
    /// 1-σ position noise in each axis (m).
    pub position_noise: f64,
}
#[allow(dead_code)]
impl GpsUnit {
    /// Create a new GPS unit.
    pub fn new(position_noise: f64) -> Self {
        Self { position_noise }
    }
    fn box_muller(rng: &mut impl rand::Rng) -> f64 {
        let u1: f64 = rng.random::<f64>().max(1e-15);
        let u2: f64 = rng.random::<f64>();
        (-2.0 * u1.ln()).sqrt() * (2.0 * PI * u2).cos()
    }
    /// Return a noisy 3-D position measurement.
    pub fn measure_position(&self, true_pos: [f64; 3], rng: &mut impl rand::Rng) -> [f64; 3] {
        [
            true_pos[0] + self.position_noise * Self::box_muller(rng),
            true_pos[1] + self.position_noise * Self::box_muller(rng),
            true_pos[2] + self.position_noise * Self::box_muller(rng),
        ]
    }
}
/// Simplified optical-flow sensor estimating 2-D ground-plane velocity.
///
/// Uses the pinhole projection formula: flow = v * f / h
/// where `v` is the lateral/longitudinal velocity, `f` is the focal length,
/// and `h` is the camera height above the ground plane.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct OpticalFlowSensor {
    /// Camera focal length (pixels).
    pub focal_length_px: f64,
    /// Camera height above ground plane (m).
    pub height_m: f64,
    /// 1-σ measurement noise (pixels/frame).
    pub noise_px: f64,
}
impl OpticalFlowSensor {
    /// Create a typical downward-facing camera for a ground vehicle.
    #[allow(dead_code)]
    pub fn ground_vehicle_camera() -> Self {
        Self {
            focal_length_px: 600.0,
            height_m: 0.3,
            noise_px: 1.5,
        }
    }
    /// Estimate ground-plane velocity `[vx, vy]` (m/s) from optical flow.
    ///
    /// * `flow_x_px`, `flow_y_px` – measured optical flow (pixels/frame)
    /// * `dt`                      – frame period (s)
    #[allow(dead_code)]
    pub fn velocity_from_flow(&self, flow_x_px: f64, flow_y_px: f64, dt: f64) -> [f64; 2] {
        if dt < 1e-12 || self.focal_length_px < 1.0 {
            return [0.0; 2];
        }
        let scale = self.height_m / (self.focal_length_px * dt);
        [flow_x_px * scale, flow_y_px * scale]
    }
    /// Simulate optical-flow measurement from true velocity.
    ///
    /// Returns `(flow_x_px, flow_y_px)` per frame.
    #[allow(dead_code)]
    pub fn simulate_flow(
        &self,
        vx: f64,
        vy: f64,
        dt: f64,
        noise_sample_x: f64,
        noise_sample_y: f64,
    ) -> (f64, f64) {
        let scale = self.focal_length_px * dt / self.height_m;
        let fx = vx * scale + self.noise_px * noise_sample_x;
        let fy = vy * scale + self.noise_px * noise_sample_y;
        (fx, fy)
    }
}
/// Three independent scalar Kalman filters, one per axis, providing a
/// 3-D position and velocity estimate.
#[allow(dead_code)]
pub struct Kalman3D {
    /// Filter for the X axis.
    pub kf_x: KalmanFilter1D,
    /// Filter for the Y axis.
    pub kf_y: KalmanFilter1D,
    /// Filter for the Z axis.
    pub kf_z: KalmanFilter1D,
}
impl Kalman3D {
    /// Create three independent Kalman filters with the same noise parameters.
    #[allow(dead_code)]
    pub fn new(q_pos: f64, q_vel: f64, r_obs: f64) -> Self {
        Self {
            kf_x: KalmanFilter1D::new(q_pos, q_vel, r_obs),
            kf_y: KalmanFilter1D::new(q_pos, q_vel, r_obs),
            kf_z: KalmanFilter1D::new(q_pos, q_vel, r_obs),
        }
    }
    /// Predict all three filters by `dt` seconds.
    #[allow(dead_code)]
    pub fn predict(&mut self, dt: f64) {
        self.kf_x.predict(dt);
        self.kf_y.predict(dt);
        self.kf_z.predict(dt);
    }
    /// Update all three filters with a 3-D position observation.
    #[allow(dead_code)]
    pub fn update(&mut self, z: [f64; 3]) {
        self.kf_x.update(z[0]);
        self.kf_y.update(z[1]);
        self.kf_z.update(z[2]);
    }
    /// Fused position estimate `[x, y, z]`.
    #[allow(dead_code)]
    pub fn position(&self) -> [f64; 3] {
        [self.kf_x.position, self.kf_y.position, self.kf_z.position]
    }
    /// Fused velocity estimate `[vx, vy, vz]`.
    #[allow(dead_code)]
    pub fn velocity(&self) -> [f64; 3] {
        [self.kf_x.velocity, self.kf_y.velocity, self.kf_z.velocity]
    }
    /// RMS position standard deviation across all three axes.
    #[allow(dead_code)]
    pub fn position_std_rms(&self) -> f64 {
        let v = (self.kf_x.p_pos + self.kf_y.p_pos + self.kf_z.p_pos) / 3.0;
        v.max(0.0).sqrt()
    }
}
/// Dead-reckoning navigator that integrates IMU acceleration to propagate
/// position and velocity.
///
/// Operates in 2-D (X, Y) in the vehicle's local horizontal plane.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct ImuDeadReckoning {
    /// Current position `[x, y]` (m).
    pub position: [f64; 2],
    /// Current velocity `[vx, vy]` (m/s).
    pub velocity: [f64; 2],
    /// Current heading (yaw angle, rad from +X).
    pub heading: f64,
    /// Accumulated yaw angle from gyro integration (rad).
    pub integrated_yaw: f64,
}
impl ImuDeadReckoning {
    /// Create a new dead-reckoning navigator starting at the origin.
    #[allow(dead_code)]
    pub fn new() -> Self {
        Self {
            position: [0.0; 2],
            velocity: [0.0; 2],
            heading: 0.0,
            integrated_yaw: 0.0,
        }
    }
    /// Propagate the navigator by one IMU sample.
    ///
    /// * `ax`, `ay` – body-frame accelerations (m/s²), X = forward, Y = left
    /// * `yaw_rate` – yaw rate from gyroscope (rad/s)
    /// * `dt`       – time step (s)
    #[allow(dead_code)]
    pub fn step(&mut self, ax: f64, ay: f64, yaw_rate: f64, dt: f64) {
        self.heading += yaw_rate * dt;
        self.integrated_yaw += yaw_rate * dt;
        let cos_h = self.heading.cos();
        let sin_h = self.heading.sin();
        let ax_w = ax * cos_h - ay * sin_h;
        let ay_w = ax * sin_h + ay * cos_h;
        self.velocity[0] += ax_w * dt;
        self.velocity[1] += ay_w * dt;
        self.position[0] += self.velocity[0] * dt;
        self.position[1] += self.velocity[1] * dt;
    }
    /// Speed magnitude (m/s).
    #[allow(dead_code)]
    pub fn speed(&self) -> f64 {
        (self.velocity[0].powi(2) + self.velocity[1].powi(2)).sqrt()
    }
    /// Reset position to the given value (e.g. after a GPS correction).
    #[allow(dead_code)]
    pub fn reset_position(&mut self, pos: [f64; 2]) {
        self.position = pos;
    }
    /// Reset velocity to zero (e.g. at startup).
    #[allow(dead_code)]
    pub fn reset_velocity(&mut self) {
        self.velocity = [0.0; 2];
    }
}
/// A single IMU measurement sample.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct ImuMeasurement {
    /// Measured specific force (acceleration minus gravity) `[ax, ay, az]` in m/s².
    pub accel: [f64; 3],
    /// Measured angular velocity `[ωx, ωy, ωz]` in rad/s.
    pub gyro: [f64; 3],
}
/// IMU bias instability model using a Gauss-Markov (first-order Markov) process.
///
/// Models slowly varying accelerometer and gyroscope biases that are
/// not captured by the constant-bias model in `ImuUnit`.
///
/// The bias evolves as:
/// `b(t+dt) = b(t) * exp(-dt/T) + σ_d * w`
///
/// where T is the correlation time and σ_d is the diffusion coefficient.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct ImuBiasModel {
    /// Accelerometer bias `[bx, by, bz]` (m/s²).
    pub accel_bias: [f64; 3],
    /// Gyroscope bias `[bx, by, bz]` (rad/s).
    pub gyro_bias: [f64; 3],
    /// Accelerometer bias correlation time (s).
    pub accel_time_const: f64,
    /// Gyroscope bias correlation time (s).
    pub gyro_time_const: f64,
    /// Accelerometer bias instability σ (m/s²).
    pub accel_instability: f64,
    /// Gyroscope bias instability σ (rad/s).
    pub gyro_instability: f64,
}
impl ImuBiasModel {
    /// Create a new IMU bias model.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        accel_time_const: f64,
        gyro_time_const: f64,
        accel_instability: f64,
        gyro_instability: f64,
    ) -> Self {
        Self {
            accel_bias: [0.0; 3],
            gyro_bias: [0.0; 3],
            accel_time_const: accel_time_const.max(1.0),
            gyro_time_const: gyro_time_const.max(1.0),
            accel_instability,
            gyro_instability,
        }
    }
    /// Automotive MEMS bias model (short correlation time, moderate instability).
    pub fn automotive_mems() -> Self {
        Self::new(100.0, 500.0, 0.003, 0.0002)
    }
    /// Tactical-grade IMU (long correlation, low instability).
    pub fn tactical_grade() -> Self {
        Self::new(3600.0, 7200.0, 0.0001, 0.000005)
    }
    /// Advance bias states by `dt` seconds with white noise inputs.
    ///
    /// * `accel_noise` – 3-D unit-normal noise sample for accelerometer
    /// * `gyro_noise`  – 3-D unit-normal noise sample for gyroscope
    pub fn step(&mut self, dt: f64, accel_noise: [f64; 3], gyro_noise: [f64; 3]) {
        let decay_a = (-dt / self.accel_time_const).exp();
        let sigma_a = self.accel_instability * (1.0 - decay_a * decay_a).max(0.0).sqrt();
        let decay_g = (-dt / self.gyro_time_const).exp();
        let sigma_g = self.gyro_instability * (1.0 - decay_g * decay_g).max(0.0).sqrt();
        for i in 0..3 {
            self.accel_bias[i] = self.accel_bias[i] * decay_a + sigma_a * accel_noise[i];
            self.gyro_bias[i] = self.gyro_bias[i] * decay_g + sigma_g * gyro_noise[i];
        }
    }
    /// Apply current bias to a true accelerometer measurement.
    pub fn corrupt_accel(
        &self,
        true_accel: [f64; 3],
        white_noise: [f64; 3],
        std_dev: f64,
    ) -> [f64; 3] {
        [
            true_accel[0] + self.accel_bias[0] + std_dev * white_noise[0],
            true_accel[1] + self.accel_bias[1] + std_dev * white_noise[1],
            true_accel[2] + self.accel_bias[2] + std_dev * white_noise[2],
        ]
    }
    /// Apply current bias to a true gyroscope measurement.
    pub fn corrupt_gyro(
        &self,
        true_gyro: [f64; 3],
        white_noise: [f64; 3],
        std_dev: f64,
    ) -> [f64; 3] {
        [
            true_gyro[0] + self.gyro_bias[0] + std_dev * white_noise[0],
            true_gyro[1] + self.gyro_bias[1] + std_dev * white_noise[1],
            true_gyro[2] + self.gyro_bias[2] + std_dev * white_noise[2],
        ]
    }
    /// Reset both biases to zero.
    pub fn reset(&mut self) {
        self.accel_bias = [0.0; 3];
        self.gyro_bias = [0.0; 3];
    }
    /// RMS bias magnitude for accelerometer (m/s²).
    pub fn accel_bias_rms(&self) -> f64 {
        let s: f64 = self.accel_bias.iter().map(|b| b * b).sum();
        (s / 3.0).sqrt()
    }
    /// RMS bias magnitude for gyroscope (rad/s).
    pub fn gyro_bias_rms(&self) -> f64 {
        let s: f64 = self.gyro_bias.iter().map(|b| b * b).sum();
        (s / 3.0).sqrt()
    }
}
/// Simple resistive temperature sensor model.
///
/// Measures temperature with a constant bias and Gaussian noise.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct TemperatureSensor {
    /// 1-σ noise in °C.
    pub noise_c: f64,
    /// Constant bias in °C.
    pub bias_c: f64,
    /// Update rate in Hz.
    pub update_rate_hz: f64,
}
impl TemperatureSensor {
    /// Typical thermocouple sensor for engine bay use.
    #[allow(dead_code)]
    pub fn thermocouple() -> Self {
        Self {
            noise_c: 1.0,
            bias_c: 0.5,
            update_rate_hz: 10.0,
        }
    }
    /// High-accuracy RTD (Pt100) for tyre temperature.
    #[allow(dead_code)]
    pub fn rtd_tyre() -> Self {
        Self {
            noise_c: 0.2,
            bias_c: 0.0,
            update_rate_hz: 50.0,
        }
    }
    /// Measure temperature (°C).
    #[allow(dead_code)]
    pub fn measure(&self, true_temp_c: f64, noise_sample: f64) -> f64 {
        true_temp_c + self.bias_c + self.noise_c * noise_sample
    }
}
/// IMU noise model parameterised by Allan deviation coefficients.
///
/// Allan deviation characterises sensor noise over different averaging times τ:
/// - Quantisation noise:  σ_Q = N / √(3·τ)   (not modelled here explicitly)
/// - Angle random walk:   σ_ARW = N / √τ
/// - Bias instability:    σ_BI ≈ 0.664 * B   (constant plateau)
/// - Rate random walk:    σ_RRW = K · √τ
///
/// All coefficients are for the accelerometer (in m/s²); gyro uses the same
/// structure but in rad/s.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct AllenDeviationImu {
    /// Accelerometer white-noise coefficient N (m/s² / √Hz).
    pub accel_n: f64,
    /// Accelerometer bias instability B (m/s²).
    pub accel_b: f64,
    /// Accelerometer random-walk coefficient K (m/s³ / √Hz, i.e. rate random walk).
    pub accel_k: f64,
    /// Gyroscope white-noise coefficient (rad/s / √Hz).
    pub gyro_n: f64,
}
impl AllenDeviationImu {
    /// Create a new Allan-deviation IMU noise model.
    #[allow(dead_code)]
    pub fn new(accel_n: f64, accel_b: f64, accel_k: f64, gyro_n: f64) -> Self {
        Self {
            accel_n,
            accel_b,
            accel_k,
            gyro_n,
        }
    }
    /// Automotive MEMS preset.
    ///
    /// Typical values for a consumer-grade automotive accelerometer.
    #[allow(dead_code)]
    pub fn automotive_mems() -> Self {
        Self {
            accel_n: 0.05,
            accel_b: 0.002,
            accel_k: 0.001,
            gyro_n: 0.001,
        }
    }
    /// Compute the accelerometer noise standard deviation at averaging time τ (s).
    ///
    /// Combined Allan deviation (simplified):
    /// `σ(τ) = sqrt((N/√τ)² + (0.664·B)² + (K·√τ)²)`
    #[allow(dead_code)]
    pub fn accel_noise_at_tau(&self, tau: f64) -> f64 {
        if tau <= 0.0 {
            return self.accel_n;
        }
        let arw = self.accel_n / tau.sqrt();
        let bi = 0.664 * self.accel_b;
        let rrw = self.accel_k * tau.sqrt();
        (arw * arw + bi * bi + rrw * rrw).sqrt()
    }
    /// Compute the gyroscope noise standard deviation at averaging time τ (s).
    #[allow(dead_code)]
    pub fn gyro_noise_at_tau(&self, tau: f64) -> f64 {
        if tau <= 0.0 {
            return self.gyro_n;
        }
        self.gyro_n / tau.sqrt()
    }
}
/// Constant False Alarm Rate (CFAR) detection parameters.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct CfarParams {
    /// Number of guard cells on each side of the cell under test.
    pub n_guard: usize,
    /// Number of training cells on each side.
    pub n_training: usize,
    /// False alarm rate target (e.g. 1e-6).
    pub pfa: f64,
}
impl CfarParams {
    /// Default CFAR for automotive radar.
    #[allow(dead_code)]
    pub fn automotive() -> Self {
        Self {
            n_guard: 2,
            n_training: 8,
            pfa: 1e-5,
        }
    }
    /// Detection threshold factor from PFA and number of training cells.
    ///
    /// For CA-CFAR: `T = N * (PFA^(-1/N) - 1)`
    pub fn threshold_factor(&self) -> f64 {
        let n = self.n_training as f64;
        if n < 1.0 || self.pfa <= 0.0 {
            return 10.0;
        }
        n * (self.pfa.powf(-1.0 / n) - 1.0)
    }
}
/// A GPS measurement sample.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct GpsMeasurement {
    /// Estimated position `[x, y, z]` in metres (local ENU frame).
    pub position: [f64; 3],
    /// Estimated velocity `[vx, vy, vz]` in m/s.
    pub velocity: [f64; 3],
    /// Number of satellites in view.
    pub satellites: u8,
    /// Horizontal Dilution of Precision (lower = better).
    pub hdop: f64,
}
/// Extended complementary filter for roll and pitch attitude estimation.
///
/// Blends gyroscope integration (high-pass) with accelerometer tilt estimate
/// (low-pass):
///
/// `angle_new = alpha * (angle + omega * dt) + (1 - alpha) * angle_accel`
///
/// Unlike the basic `SensorFusion`, this struct is self-contained and does
/// not delegate to `ImuSensor::accel_attitude`.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct ComplementaryFilter {
    /// Roll estimate (radians).
    pub roll: f64,
    /// Pitch estimate (radians).
    pub pitch: f64,
    /// High-pass coefficient for gyro (typical 0.95–0.99).
    pub alpha: f64,
}
impl ComplementaryFilter {
    /// Create a new filter at zero attitude.
    #[allow(dead_code)]
    pub fn new(alpha: f64) -> Self {
        Self {
            roll: 0.0,
            pitch: 0.0,
            alpha: alpha.clamp(0.0, 1.0),
        }
    }
    /// Update attitude with a new IMU sample.
    ///
    /// # Arguments
    /// * `accel` – specific force `[ax, ay, az]` (m/s²)
    /// * `gyro`  – angular velocity `[ωx, ωy, ωz]` (rad/s), ωx=roll rate, ωy=pitch rate
    /// * `dt`    – time step (s)
    #[allow(dead_code)]
    pub fn update(&mut self, accel: [f64; 3], gyro: [f64; 3], dt: f64) {
        let roll_gyro = self.roll + gyro[0] * dt;
        let pitch_gyro = self.pitch + gyro[1] * dt;
        let roll_accel = accel[1].atan2(accel[2]);
        let pitch_accel = (-accel[0]).atan2((accel[1] * accel[1] + accel[2] * accel[2]).sqrt());
        self.roll = self.alpha * roll_gyro + (1.0 - self.alpha) * roll_accel;
        self.pitch = self.alpha * pitch_gyro + (1.0 - self.alpha) * pitch_accel;
    }
    /// Reset to zero attitude.
    #[allow(dead_code)]
    pub fn reset(&mut self) {
        self.roll = 0.0;
        self.pitch = 0.0;
    }
}
/// Fuses four individual wheel speed measurements into a single vehicle speed
/// estimate, handling outliers from spinning/locked wheels.
///
/// Uses a median-of-four approach: excludes the highest and lowest readings,
/// then averages the two middle values.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct WheelSpeedFusion {
    /// Effective wheel rolling radius (m).
    pub wheel_radius: f64,
    /// Maximum allowed individual-wheel deviation from mean before exclusion (m/s).
    pub outlier_threshold: f64,
}
impl WheelSpeedFusion {
    /// Create a new wheel speed fusion object.
    pub fn new(wheel_radius: f64, outlier_threshold: f64) -> Self {
        Self {
            wheel_radius: wheel_radius.max(0.01),
            outlier_threshold: outlier_threshold.abs(),
        }
    }
    /// Typical passenger-car ABS wheel speed fusion.
    pub fn typical() -> Self {
        Self::new(0.31, 3.0)
    }
    /// Fuse four wheel angular velocities (rad/s) into a vehicle speed estimate (m/s).
    ///
    /// Applies outlier rejection: wheels whose speed deviates more than
    /// `outlier_threshold` from the mean are excluded from the final average.
    pub fn fuse(&self, omega: [f64; 4]) -> f64 {
        let speeds: Vec<f64> = omega.iter().map(|&w| w * self.wheel_radius).collect();
        let mean = speeds.iter().sum::<f64>() / 4.0;
        let valid: Vec<f64> = speeds
            .iter()
            .copied()
            .filter(|&s| (s - mean).abs() <= self.outlier_threshold)
            .collect();
        if valid.is_empty() {
            return mean;
        }
        valid.iter().sum::<f64>() / valid.len() as f64
    }
    /// Detect wheel slip: returns `true` for each wheel that is likely spinning
    /// (speed more than `outlier_threshold` above the fused speed).
    pub fn detect_slip(&self, omega: [f64; 4]) -> [bool; 4] {
        let fused = self.fuse(omega);
        let mut slip = [false; 4];
        for (i, &w) in omega.iter().enumerate() {
            let linear = w * self.wheel_radius;
            slip[i] = linear - fused > self.outlier_threshold;
        }
        slip
    }
}
/// Stateless GPS sensor functions.
#[allow(dead_code)]
pub struct GpsSensor;
impl GpsSensor {
    /// Simulate a GPS position measurement.
    ///
    /// # Arguments
    /// * `config`       – GPS configuration
    /// * `true_pos`     – true position `[x, y, z]`
    /// * `true_vel`     – true velocity `[vx, vy, vz]`
    /// * `noise_h`      – unit-normal horizontal noise sample
    /// * `noise_v`      – unit-normal vertical noise sample
    /// * `satellites`   – number of satellites in view
    #[allow(dead_code)]
    pub fn measure(
        config: &GpsConfig,
        true_pos: [f64; 3],
        true_vel: [f64; 3],
        noise_h: f64,
        noise_v: f64,
        satellites: u8,
    ) -> GpsMeasurement {
        let hdop = 1.0 + (12_u8.saturating_sub(satellites)) as f64 * 0.2;
        GpsMeasurement {
            position: [
                true_pos[0] + config.horizontal_accuracy_m * noise_h,
                true_pos[1] + config.horizontal_accuracy_m * noise_h,
                true_pos[2] + config.vertical_accuracy_m * noise_v,
            ],
            velocity: true_vel,
            satellites,
            hdop,
        }
    }
}
/// Configuration for a GPS receiver.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct GpsConfig {
    /// 1-σ horizontal position error in metres.
    pub horizontal_accuracy_m: f64,
    /// 1-σ vertical position error in metres.
    pub vertical_accuracy_m: f64,
    /// Update rate in Hz.
    pub update_rate_hz: f64,
}
impl GpsConfig {
    /// Typical automotive-grade GPS (e.g. u-blox M8).
    #[allow(dead_code)]
    pub fn automotive() -> Self {
        Self {
            horizontal_accuracy_m: 2.5,
            vertical_accuracy_m: 3.0,
            update_rate_hz: 10.0,
        }
    }
    /// High-precision RTK GPS.
    #[allow(dead_code)]
    pub fn rtk() -> Self {
        Self {
            horizontal_accuracy_m: 0.02,
            vertical_accuracy_m: 0.04,
            update_rate_hz: 20.0,
        }
    }
}
/// Wheel speed sensor based on a toothed ring (ABS/encoder type).
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct WheelPulse {
    /// Number of teeth on the tone wheel.
    pub n_teeth: u32,
    /// 1-σ angular velocity noise (rad/s), applied to speed estimate.
    pub noise: f64,
}
#[allow(dead_code)]
impl WheelPulse {
    /// Create a new wheel pulse sensor.
    pub fn new(n_teeth: u32, noise: f64) -> Self {
        Self { n_teeth, noise }
    }
    /// Count the number of pulses generated in time `dt` (s) at angular
    /// velocity `omega` (rad/s).
    ///
    /// pulses = floor(omega * dt * n_teeth / (2π))
    pub fn pulse_count(&self, omega: f64, dt: f64) -> u32 {
        let raw = omega.abs() * dt * self.n_teeth as f64 / (2.0 * PI);
        raw.floor() as u32
    }
    /// Estimate vehicle speed from `pulses` counted over `dt` seconds and
    /// tyre circumference `circumference` (m).
    ///
    /// v = (pulses / n_teeth) * circumference / dt
    pub fn velocity_from_pulses(&self, pulses: u32, dt: f64, circumference: f64) -> f64 {
        if dt < 1e-15 || self.n_teeth == 0 {
            return 0.0;
        }
        (pulses as f64 / self.n_teeth as f64) * circumference / dt
    }
}
/// Barometric pressure sensor for altitude estimation.
///
/// Uses the international barometric formula:
/// `p = p0 * (1 - L*h / T0)^(g*M/(R*L))`
/// Simplified to the linear approximation for small altitude differences:
/// `h ≈ (p0 - p) / (rho * g)`
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct BarometricSensor {
    /// Sea-level pressure in Pascals (e.g. 101 325 Pa).
    pub p0: f64,
    /// Air density at sea level (kg/m³, ~1.225).
    pub rho: f64,
    /// Gravitational acceleration (m/s²).
    pub g: f64,
    /// 1-σ pressure noise in Pascals.
    pub noise_pa: f64,
}
impl BarometricSensor {
    /// Standard atmosphere sensor preset.
    #[allow(dead_code)]
    pub fn standard() -> Self {
        Self {
            p0: 101_325.0,
            rho: 1.225,
            g: 9.81,
            noise_pa: 5.0,
        }
    }
    /// Simulate a pressure measurement at altitude `h_true` (m).
    ///
    /// Returns the measured pressure (Pa).
    #[allow(dead_code)]
    pub fn measure_pressure(&self, h_true: f64, noise_sample: f64) -> f64 {
        let true_p = self.p0 - self.rho * self.g * h_true;
        true_p + self.noise_pa * noise_sample
    }
    /// Estimate altitude from a pressure measurement.
    ///
    /// Returns altitude in metres above the reference level.
    #[allow(dead_code)]
    pub fn altitude_from_pressure(&self, pressure: f64) -> f64 {
        (self.p0 - pressure) / (self.rho * self.g)
    }
}
/// A wheel-speed-based vehicle speed sensor with Gaussian noise.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct SpeedSensor {
    /// Effective wheel rolling radius (m).
    pub wheel_radius: f64,
    /// 1-σ noise standard deviation (m/s).
    pub noise_std: f64,
}
#[allow(dead_code)]
impl SpeedSensor {
    /// Create a new speed sensor.
    pub fn new(wheel_radius: f64, noise_std: f64) -> Self {
        Self {
            wheel_radius,
            noise_std,
        }
    }
    /// Measure vehicle speed from wheel angular velocity `omega` (rad/s).
    ///
    /// Returns the noisy speed measurement (m/s).
    pub fn measure(&self, omega: f64, rng: &mut impl rand::Rng) -> f64 {
        let true_speed = omega * self.wheel_radius;
        let noise: f64 = rng.random::<f64>() * 2.0 - 1.0;
        let u1: f64 = rng.random::<f64>().max(1e-15);
        let u2: f64 = rng.random::<f64>();
        let gaussian = (-2.0 * u1.ln()).sqrt() * (2.0 * PI * u2).cos();
        let _ = noise;
        true_speed + self.noise_std * gaussian
    }
}
/// Sensor fusion combining GPS position with IMU acceleration via a simple
/// complementary filter.
///
/// State:  fused_velocity (scalar, m/s),  fused_position (\[x, y, z\], m).
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct PositionFusion {
    /// Current fused speed estimate (m/s).
    pub fused_velocity: f64,
    /// Current fused position estimate `[x, y, z]` (m).
    pub fused_position: [f64; 3],
    /// Complementary filter weight for GPS (0 = trust only IMU, 1 = trust only GPS).
    pub gps_weight: f64,
}
#[allow(dead_code)]
impl PositionFusion {
    /// Create a new fusion object, starting at the origin.
    pub fn new(gps_weight: f64) -> Self {
        Self {
            fused_velocity: 0.0,
            fused_position: [0.0; 3],
            gps_weight: gps_weight.clamp(0.0, 1.0),
        }
    }
    /// Update the fused state with a new GPS position and IMU acceleration.
    ///
    /// Returns the updated fused position.
    ///
    /// The IMU is integrated to predict position; the GPS provides a correction.
    pub fn update(&mut self, gps: [f64; 3], imu_accel: [f64; 3], dt: f64) -> [f64; 3] {
        let accel_mag = (imu_accel[0] * imu_accel[0]
            + imu_accel[1] * imu_accel[1]
            + imu_accel[2] * imu_accel[2])
            .sqrt();
        self.fused_velocity += accel_mag * dt;
        let predicted = [
            self.fused_position[0] + imu_accel[0] * dt * dt * 0.5,
            self.fused_position[1] + imu_accel[1] * dt * dt * 0.5,
            self.fused_position[2] + imu_accel[2] * dt * dt * 0.5,
        ];
        let w = self.gps_weight;
        self.fused_position = [
            (1.0 - w) * predicted[0] + w * gps[0],
            (1.0 - w) * predicted[1] + w * gps[1],
            (1.0 - w) * predicted[2] + w * gps[2],
        ];
        self.fused_position
    }
}
